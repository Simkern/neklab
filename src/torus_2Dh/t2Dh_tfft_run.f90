      module t2Dh_tfft_run
      !! The userchk-facing wrappers around the accumulator in t2Dh_tfft.
      !!
      !!   tfft_f2d          : temporal Fourier series of the converged
      !!                       baseflow orbit, truncated at a user tolerance,
      !!                       written as f2r/f2i coefficient files.
      !!
      !!   tfft_p2d          : temporal Fourier series of one or more Floquet
      !!                       modes (periodic part), truncated at a user
      !!                       tolerance, written as p2r/p2i coefficient files.
      !!
      !!   tfft_baseflow     : legacy: full-M transform of the orbit, spectrum
      !!   tfft_perturbation   file and Nek outposts, no coefficient files.
      !!
      !! The baseflow transforms sweep the recorded orbit (no solves, no time
      !! stepping: the buffer is a complete record of the trajectory). The mode
      !! transforms re-run the linear solver over one period on the replayed
      !! baseflow and transform the mode on the fly, because the perturbation
      !! fields are not stored anywhere.
      !!
      !! All wrappers work in TWO modes:
      !!
      !!   in-solver (default): the buffer is already live from the Newton or
      !!      stability run that just finished; pass no extra arguments.
      !!
      !!   post-processing: call from a fresh job by passing base= (and prefix=
      !!      if it is not 'n') to point at an existing set of chunk files. The
      !!      wrapper calls bf_init (a no-op if the buffer is already up) and
      !!      then bf_read_t2Dh to restore nsteps, nchunk, T and omega from the
      !!      sidecar <prefix><base>.t2Dh, after which the chunk files are
      !!      readable. The default prefix is 'n', which is where
      !!      unsteady_flowrate_newton leaves the converged orbit.
      !!
      !!--------------------------------------------------------------------
      !! COEFFICIENT FILES
      !!--------------------------------------------------------------------
      !!
      !! Same binary format as the n2d chunk files (t2Dh_bf_buffer_io), one
      !! file for the real part and one for the imaginary part of the complex
      !! coefficients. The leading '2d' of the buffer base is replaced by
      !! '2r' / '2i', so for the case torus2Dh (default base '2dtorus2Dh'):
      !!
      !!   f2rtorus2Dh001.fld , f2itorus2Dh001.fld     baseflow
      !!   p2rtorus2Dh00j.fld , p2itorus2Dh00j.fld     Floquet mode j = 1, 2, ...
      !!
      !! (a base not starting with '2d' gets 'r' / 'i' inserted after the
      !! prefix instead). Inside each file:
      !!
      !!   header time   period T, so omega = 2 pi / T
      !!   nsave         number of records
      !!   record k      (u_Z, u_R, u_phi) of Re c_m or Im c_m
      !!   dt(k)         relative energy of that record, as a fraction of the
      !!                 energy of the computed spectrum (see t2Dh_tfft)
      !!
      !! with c_m = (1/T) int_0^T f exp(-i m w t) dt, i.e. numpy's
      !! norm="forward" convention, so that
      !!
      !!   f(x,t) = sum_m c_m(x) exp(i m w t)
      !!
      !! with NO rescaling -- the cheapest form to reconstruct from.
      !!
      !!   f2d : records m = 0, 1, ..., M (one-sided, the field is real and
      !!         c_{-m} = conj(c_m)); reconstruct with
      !!         f = c_0 + 2 Re sum_{m>=1} c_m exp(i m w t)
      !!         or numpy irfft(c, n, norm="forward"). dt(k) is the energy of
      !!         the harmonic PAIR +-m, so the column sums to 1 - dropped.
      !!
      !!   p2d : records in numpy fft order m = 0, 1, ..., M, -M, ..., -1
      !!         (2M+1 records, c_{-m} independent of c_m for alpha /= 0);
      !!         reconstruct with the plain sum or numpy
      !!         ifft(c, norm="forward") after zero-padding the middle.
      !!         dt(k) is the energy of the single coefficient m.
      !!
      !! M is the smallest order whose discarded energy is below tol (see
      !! tfft_truncate). Only the retained harmonics are written, both to the
      !! coefficient files and to the optional Nek outposts.
      !!
      !!--------------------------------------------------------------------
      !! WHY THE PERTURBATION NEEDS THE MULTIPLIER
      !!--------------------------------------------------------------------
      !!
      !! A Floquet solution is u(t) = exp(sigma t) utilde(t) with utilde
      !! T-periodic and mu = exp(sigma T). The raw perturbation grows (or
      !! decays, or rotates) by mu over the period, so it is NOT periodic and
      !! its Fourier coefficients on [0,T] are leakage. What has harmonics is
      !! utilde = u * mu**(-t/T), and removing that modulation is the whole
      !! difference between the mode and the baseflow transforms.
      !!
      !! Pass the multiplier the eigensolver returned. The wrapper reports
      !! |utilde(T) - utilde(0)|/|utilde(0)| for every mode: for a converged
      !! Floquet pair it sits at the eigensolver's tolerance, and if it does
      !! not, the mode, the multiplier or the baseflow record do not belong
      !! together and no spectrum computed from them means anything.
      !!
      !!--------------------------------------------------------------------
      !! ASSUMPTION TO VERIFY ONCE
      !!--------------------------------------------------------------------
      !!
      !! The replay loop below mirrors what the Floquet propagator does:
      !!
      !!    bf_replay_start ; do k = 1, n ; bf_set(k) ; advance ; end do
      !!
      !! with istep incremented by the caller (settime uses nab = min(istep,3)
      !! for the order ramp) and dt taken from param(12), which bf_set forces
      !! to the recorded value. Check this against floquet_operator_torus_2Dh
      !! and align the setup_linear_solver arguments with it: if the propagator
      !! configures the solver differently, the mode this wrapper marches is
      !! not the mode the eigensolver found.
         use stdlib_optval, only: optval
         use LightKrylov, only: dp, atol_dp
         use neklab_vectors, only: nek_zvector
         use neklab_utils, only: vec2nek
         use neklab_t2Dh, only: t2Dh, lfc
         use neklab_2Dh, only: nek_advance_2Dh_planar, nek_advance_2Dh_axisym
         use neklab_nek_setup, only: setup_linear_solver, nek_log_message,
     &                               nek_log_information, nek_log_warning,
     &                               nek_log_debug, nek_stop_error, set_fldindex
         use t2Dh_bf_buffer, only: bf_init, bf_get_nsteps, bf_get_time, bf_get_dt,
     &                             bf_set, bf_set_prefix, bf_get_prefix, bf_get_base,
     &                             bf_replay_start, bf_replay_end, bf_is_recording,
     &                             bf_read_t2Dh, bf_write_fields
         use t2Dh_tfft, only: tfft_start, tfft_add, tfft_close, tfft_free,
     &                        tfft_demodulate, tfft_outpost, tfft_spectrum,
     &                        tfft_coef, tfft_truncate, tfft_get_period,
     &                        mfft_max
         implicit none
         include "SIZE"
         include "TOTAL"
         private
         character(len=*), parameter, private :: this_module = 't2Dh_tfft_run'

         integer, parameter :: lv = lx1*ly1*lz1*lelv
         integer, parameter :: lp2 = lx2*ly2*lz2*lelv

         public :: tfft_f2d, tfft_p2d
         public :: tfft_baseflow, tfft_perturbation

      contains

      !====================================================================
      !     f2d: BASEFLOW COEFFICIENTS
      !====================================================================

         subroutine tfft_f2d(tol, mmax, prefix, base, if_outpost, spec_file)
      !! Temporal Fourier series of the recorded orbit, truncated at tol, and
      !! written to f2r/f2i (see COEFFICIENT FILES above). Call after
      !! unsteady_flowrate_newton (in-solver) or from a fresh job with base=.
      !!
      !!   call tfft_f2d(1.0e-8_dp)              ! after the Newton solve
      !!   call tfft_f2d(1.0e-8_dp, base='2dtorus2Dh')  ! post-processing
            real(dp), intent(in) :: tol
      !! Relative energy below which the remaining harmonics are dropped.
            integer, optional, intent(in) :: mmax
      !! Highest harmonic computed before truncating. Default mfft_max.
            character(len=1), optional, intent(in) :: prefix
      !! Prefix of the orbit's chunk files. Default 'n'.
            character(len=*), optional, intent(in) :: base
      !! Filename stem. Only needed in post-processing mode.
            logical, optional, intent(in) :: if_outpost
      !! Also write the retained harmonics as Nek field files (r__ / i__
      !! series, file index = m+1). Default .true.
            character(len=*), optional, intent(in) :: spec_file
      !! Amplitude table of the full computed spectrum. Default
      !! 'spec_baseflow.dat'.
      ! internal
            character(len=*), parameter :: this_procedure = 'tfft_f2d'
            real(dp), dimension(-mfft_max:mfft_max) :: efrac
            real(dp), dimension(:), allocatable :: tags
            integer, dimension(:), allocatable :: mlist
            integer :: M, k

            call transform_baseflow(optval(mmax, mfft_max), optval(prefix, 'n'), base,
     &                              this_procedure)
            call tfft_truncate(tol, M, efrac)

      ! one-sided: m = 0..M, energy of the +-m pair
            allocate (mlist(M + 1), tags(M + 1))
            do k = 0, M
               mlist(k + 1) = k
               if (k == 0) then
                  tags(k + 1) = efrac(0)
               else
                  tags(k + 1) = efrac(k) + efrac(-k)
               end if
            end do
            call write_coef_files('f', 1, mlist, tags, this_procedure)

            call tfft_spectrum(filename=optval(spec_file, 'spec_baseflow.dat'))
            if (optval(if_outpost, .true.)) call tfft_outpost(mlist=mlist)
            call tfft_free()
            deallocate (mlist, tags)
         end subroutine tfft_f2d

      !====================================================================
      !     p2d: FLOQUET MODE COEFFICIENTS
      !====================================================================

         subroutine tfft_p2d(tol, zmodes, alpha, mu, mmax, prefix, base, vtol, ptol,
     &                       if_normalise, if_axisym, if_outpost, spec_file)
      !! For each mode j: one more linear period on the replayed orbit, the
      !! Floquet modulation removed, the periodic part transformed, truncated
      !! at tol, and written to p2r00j / p2i00j in numpy fft order (see
      !! COEFFICIENT FILES above).
      !!
      !!   call tfft_p2d(1.0e-8_dp, X(1:nev), alpha, mu(1:nev))
      !!
      !! zmodes and mu are the eigenvectors and multipliers as the stability
      !! analysis returned them, in the same order. alpha is the azimuthal
      !! wavenumber: alpha = 0 needs lpert = 1, alpha /= 0 needs lpert = 2.
            real(dp), intent(in) :: tol
      !! Relative energy below which the remaining harmonics are dropped.
            type(nek_zvector), dimension(:), intent(in) :: zmodes
            real(dp), intent(in) :: alpha
            complex(dp), dimension(:), intent(in) :: mu
      !! Floquet multiplier of each mode. Required: only the periodic part
      !! has a Fourier series.
            integer, optional, intent(in) :: mmax
      !! Highest harmonic computed before truncating. Default mfft_max.
            character(len=1), optional, intent(in) :: prefix
      !! Prefix of the orbit's chunk files. Default 'n'.
            character(len=*), optional, intent(in) :: base
      !! Filename stem. Only needed in post-processing mode.
            real(dp), optional, intent(in) :: vtol
            real(dp), optional, intent(in) :: ptol
            logical, optional, intent(in) :: if_normalise
      !! Scale each mode to unit norm at t = 0 first. Default .true.
            logical, optional, intent(in) :: if_axisym
      !! Torus (axisymmetric) linear step; .false. for the planar one.
      !! Default .true.
            logical, optional, intent(in) :: if_outpost
      !! Also write the retained coefficients as Nek field files: series
      !! r<jj> / i<jj> (Re / Im of c_m), file index = record position in the
      !! p2d file, time stamp = m. Modes 1..99 only. Default .true.
            character(len=*), optional, intent(in) :: spec_file
      !! Stem of the per-mode amplitude tables, <stem>_00j.dat. Default
      !! 'spec_mode'.
      ! internal
            character(len=*), parameter :: this_procedure = 'tfft_p2d'
            character(len=256) :: msg, sfile
            character(len=3) :: tr, ti
            real(dp), dimension(-mfft_max:mfft_max) :: efrac
            real(dp), dimension(:), allocatable :: tags
            integer, dimension(:), allocatable :: mlist
            integer :: nmode, j, M, k, nrec

            nmode = size(zmodes)
            if (nmode < 1) then
               call nek_stop_error('No modes given.', this_module, this_procedure)
            end if
            if (size(mu) < nmode) then
               write (msg, '(A,I0,A,I0,A)') 'mu has ', size(mu), ' entries for ', nmode, ' modes.'
               call nek_stop_error(msg, this_module, this_procedure)
            end if
            if (nmode > 999) then
               call nek_stop_error('At most 999 modes: the file index has three digits.',
     &            this_module, this_procedure)
            end if

            do j = 1, nmode
               write (msg, '(A,I0,A,I0,A,2(1X,E16.8))') 'mode ', j, ' of ', nmode,
     &            ': mu=', real(mu(j), dp), aimag(mu(j))
               call nek_log_message(msg, this_module, this_procedure)
               write (tr, '(A,I2.2)') 'm', min(j, 99)
               call transform_mode(optval(mmax, mfft_max), zmodes(j), alpha, mu(j),
     &                             optval(prefix, 'n'), base, vtol, ptol, if_normalise,
     &                             if_axisym, tr, this_procedure)
               call tfft_truncate(tol, M, efrac)

      ! numpy fft order: 0, 1, ..., M, -M, ..., -1
               nrec = 2*M + 1
               allocate (mlist(nrec), tags(nrec))
               do k = 0, M
                  mlist(k + 1) = k
               end do
               do k = 1, M
                  mlist(M + 1 + k) = k - M - 1
               end do
               do k = 1, nrec
                  tags(k) = efrac(mlist(k))
               end do
               call write_coef_files('p', j, mlist, tags, this_procedure)

               write (sfile, '(A,A,I3.3,A)') trim(optval(spec_file, 'spec_mode')), '_', j, '.dat'
               call tfft_spectrum(filename=trim(sfile))

               if (optval(if_outpost, .true.)) then
                  if (j <= 99) then
                     write (tr, '(A,I2.2)') 'r', j
                     write (ti, '(A,I2.2)') 'i', j
                     call outpost_coefs(tr, ti, mlist)
                  else
                     call nek_log_warning('Nek outposts are limited to modes 1..99; skipped.',
     &                  this_module, this_procedure)
                  end if
               end if

               call tfft_free()
               deallocate (mlist, tags)
            end do
         end subroutine tfft_p2d

      !====================================================================
      !     LEGACY WRAPPERS (no truncation, no coefficient files)
      !====================================================================

         subroutine tfft_baseflow(mmax, prefix, base, if_outpost, spec_file)
      !! Full-M transform of the recorded orbit: spectrum file and (default)
      !! r__/i__ Nek outposts of every harmonic. Prefer tfft_f2d.
            integer, intent(in) :: mmax
            character(len=1), optional, intent(in) :: prefix
      !! One-character filename prefix. Default 'n'.
            character(len=*), optional, intent(in) :: base
            logical, optional, intent(in) :: if_outpost
            character(len=*), optional, intent(in) :: spec_file
      ! internal
            character(len=*), parameter :: this_procedure = 'tfft_baseflow'
            call transform_baseflow(mmax, optval(prefix, 'n'), base, this_procedure)
            call tfft_spectrum(filename=optval(spec_file, 'spec_baseflow.dat'))
            if (optval(if_outpost, .true.)) call tfft_outpost()
         end subroutine tfft_baseflow

         subroutine tfft_perturbation(mmax, zmode, alpha, mu, prefix, base, vtol, ptol,
     &                                if_outpost, if_normalise, if_axisym, spec_file)
      !! Full-M transform of one Floquet mode: spectrum file and (default)
      !! rRe/iRe/rIm/iIm Nek outposts. Prefer tfft_p2d.
            integer, intent(in) :: mmax
            type(nek_zvector), intent(in) :: zmode
            real(dp), intent(in) :: alpha
            complex(dp), optional, intent(in) :: mu
      !! Floquet multiplier. Default (1,0), i.e. no demodulation, which is only
      !! right for a forced linear response with no homogeneous growth.
            character(len=1), optional, intent(in) :: prefix
      !! One-character filename prefix. Default 'n'.
            character(len=*), optional, intent(in) :: base
            real(dp), optional, intent(in) :: vtol
            real(dp), optional, intent(in) :: ptol
            logical, optional, intent(in) :: if_outpost
            logical, optional, intent(in) :: if_normalise
            logical, optional, intent(in) :: if_axisym
            character(len=*), optional, intent(in) :: spec_file
      ! internal
            character(len=*), parameter :: this_procedure = 'tfft_perturbation'
            character(len=3), dimension(4) :: tg
            integer :: nblk
            call transform_mode(mmax, zmode, alpha, optval(mu, cmplx(1.0_dp, 0.0_dp, kind=dp)),
     &                          optval(prefix, 'n'), base, vtol, ptol, if_normalise,
     &                          if_axisym, 'mode', this_procedure)
            call tfft_spectrum(filename=optval(spec_file, 'spec_mode.dat'))
            if (optval(if_outpost, .true.)) then
               nblk = merge(1, 2, alpha == 0.0_dp)
               tg = [ 'rRe', 'iRe', 'rIm', 'iIm' ]
               call tfft_outpost(tg(1:2*nblk))
            end if
         end subroutine tfft_perturbation

      !====================================================================
      !     TRANSFORM CORES (leave a closed accumulator behind)
      !====================================================================

         subroutine transform_baseflow(mmax, pfx, base, caller)
      !! Sweeps the recorded orbit into the accumulator and closes it.
            integer, intent(in) :: mmax
            character(len=1), intent(in) :: pfx
            character(len=*), optional, intent(in) :: base
            character(len=*), intent(in) :: caller
      ! internal
            character(len=256) :: msg
            real(dp), dimension(:, :), allocatable :: w, w0, vsave
            character(len=1) :: pfx_save
            real(dp) :: tk, dtk, T
            integer :: k, n, nbf

            call bf_attach(pfx, base, caller)
            n = bf_get_nsteps()
            T = bf_get_time()
            nbf = lx1*ly1*lz1*nelv
            write (msg, '(A,I0,A,E16.8,A,A,A)') 'transforming ', n, ' snapshots, T= ', T,
     &         ', prefix ''', pfx, ''''
            call nek_log_message(msg, this_module, caller)

            allocate (w(lv, 3), w0(lv, 3), vsave(lv, 3))
      ! The sweep overwrites vx/vy/t with recorded snapshots, so put back what
      ! the caller had. userchk is not a place to leave the state disturbed.
            call pack_base(vsave, nbf)

            pfx_save = bf_get_prefix()
            call bf_set_prefix(pfx)
            call bf_replay_start()

      ! bf_begin_step stores the state ENTERING step k, so snapshot k is
      ! u(t_{k-1}): snapshot 1 IS u(0), and u(T) is not stored at all --
      ! snapshot n is u(T - dt_n). Interval k therefore runs from snapshot k
      ! to snapshot k+1, and the last one closes on u(T) = u(0), whose error
      ! is |F(X)| carried by a single dt/T weight.
            call bf_set(1, if_lag=.false.)
            call pack_base(w0, nbf)
            call tfft_start(mmax, 3, w0, 'baseflow')

            tk = 0.0_dp
            do k = 1, n
               dtk = bf_get_dt(k)
               tk = tk + dtk
               if (k < n) then
                  call bf_set(k + 1, if_lag=.false.)
                  call pack_base(w, nbf)
               else
                  call copy(w, w0, lv*3)
               end if
               call tfft_add(w, tk, dtk)
            end do

            call bf_replay_end()
            call bf_set_prefix(pfx_save)
            call tfft_close(T)

            call unpack_base(vsave, nbf)
            deallocate (w, w0, vsave)
         end subroutine transform_baseflow

         subroutine transform_mode(mmax, zmode, alpha, mu, pfx, base, vtol, ptol,
     &                             if_normalise, if_axisym, tag, caller)
      !! One linear period on the replayed orbit, demodulated by mu, into the
      !! accumulator, which is left closed.
            integer, intent(in) :: mmax
            type(nek_zvector), intent(in) :: zmode
            real(dp), intent(in) :: alpha
            complex(dp), intent(in) :: mu
            character(len=1), intent(in) :: pfx
            character(len=*), optional, intent(in) :: base
            real(dp), optional, intent(in) :: vtol
            real(dp), optional, intent(in) :: ptol
            logical, optional, intent(in) :: if_normalise
            logical, optional, intent(in) :: if_axisym
            character(len=*), intent(in) :: tag
            character(len=*), intent(in) :: caller
      ! internal
            character(len=256) :: msg
            real(dp), dimension(:, :), allocatable :: w, w0
            character(len=1) :: pfx_save
            real(dp) :: T, dtk, mu_abs, mu_arg, nrm0, nrmT, dev, scl
            integer :: k, n, nbf, nblk, nf, i
            logical :: axisym_

            call bf_attach(pfx, base, caller)
            n = bf_get_nsteps()
            T = bf_get_time()
            nbf = lx1*ly1*lz1*nelv
            axisym_ = optval(if_axisym, .true.)
            mu_abs = abs(mu)
            mu_arg = atan2(aimag(mu), real(mu, kind=dp))
            if (mu_abs <= 0.0_dp) then
               call nek_stop_error('|mu| must be positive.', this_module, caller)
            end if

      ! --- geometry of the mode
            if (alpha == 0.0_dp) then
               nblk = 1
               if (lpert /= 1) call nek_stop_error('alpha = 0 requires lpert = 1.',
     &            this_module, caller)
            else
               nblk = 2
               if (lpert /= 2) call nek_stop_error('alpha /= 0 requires lpert = 2.',
     &            this_module, caller)
            end if
            nf = 3*nblk

            write (msg, '(A,I0,A,E16.8,A,E16.8,A,E16.8)') 'linear replay over ', n,
     &         ' steps, alpha= ', alpha, ', |mu|= ', mu_abs, ', arg(mu)= ', mu_arg
            call nek_log_message(msg, this_module, caller)
            if (abs(mu_abs - 1.0_dp) < atol_dp .and. abs(mu_arg) < atol_dp) then
               call nek_log_warning('mu = 1: no demodulation will be applied. That is correct '//
     &            'only for a forced response with no homogeneous growth, not for a Floquet mode.',
     &            this_module, caller)
            end if

      ! --- configure the linear solver exactly as the propagator does. dt is
      !     irrelevant here (bf_set forces the recorded step every step) but it
      !     must be sane before setup_linear_solver derives nsteps from it.
            dt = T/real(n, dp)
            call setup_linear_solver(transpose=.false., solve_baseflow=.false.,
     &                               recompute_dt=.false., variable_dt=.false.,
     &                               endtime=T, solve_temperature=.true.,
     &                               vtol=optval(vtol, 1.0e-10_dp),
     &                               ptol=optval(ptol, 1.0e-10_dp), silent=.true.)

      ! --- install the mode
            call vec2nek(vxp, vyp, vzp, prp, tp, zmode)
            allocate (w(lv, nf), w0(lv, nf))
            call pack_pert(w, nbf, nblk)
            nrm0 = sqrt(fnorm2(w, nbf, nf))
            if (nrm0 <= atol_dp) then
               call nek_stop_error('The mode has zero norm.', this_module, caller)
            end if
      ! A Floquet mode is defined up to a complex constant. Normalising here
      ! makes the harmonic amplitudes comparable between modes and between runs
      ! instead of inheriting whatever scaling the eigensolver left behind.
            if (optval(if_normalise, .true.)) then
               scl = 1.0_dp/nrm0
               do i = 1, lpert
                  call cmult(vxp(1, i), scl, nbf)
                  call cmult(vyp(1, i), scl, nbf)
                  call cmult(tp(1, 1, i), scl, nbf)
               end do
               call pack_pert(w, nbf, nblk)
               nrm0 = sqrt(fnorm2(w, nbf, nf))
            end if
            call copy(w0, w, lv*nf)

      ! --- replay
            pfx_save = bf_get_prefix()
            call bf_set_prefix(pfx)
            call bf_replay_start()

            time = 0.0_dp
            istep = 0
            lastep = 0
            call tfft_start(mmax, nf, w, tag)

            do k = 1, n
               istep = istep + 1
      ! if_lag defaults to .true.: the linear solve needs the baseflow lag
      ! levels, and bf_set builds them from the traversal order.
               call bf_set(k)
               if (axisym_) then
                  call nek_advance_2Dh_axisym(alpha)
               else
                  call nek_advance_2Dh_planar(alpha)
               end if
               dtk = dt
               call pack_pert(w, nbf, nblk)
      ! Strip exp(sigma t) BEFORE accumulating. Everything downstream assumes
      ! a T-periodic integrand.
               call tfft_demodulate(w, time, T, mu_abs, mu_arg)
               call tfft_add(w, time, dtk)
            end do

            call bf_replay_end()
            call bf_set_prefix(pfx_save)

      ! --- how periodic was it really
            nrmT = sqrt(fnorm2(w, nbf, nf))
            call w_minus_w0(w, w0, nbf, nf)
            dev = sqrt(fnorm2(w, nbf, nf))/max(nrm0, atol_dp)
            write (msg, '(A,E16.8,A,E16.8)') 'Floquet closure |utilde(T)-utilde(0)|/|utilde(0)|= ',
     &         dev, ', |utilde(T)|/|utilde(0)|= ', nrmT/max(nrm0, atol_dp)
            call nek_log_message(msg, this_module, caller)
            if (dev > 1.0e-04_dp) then
               call nek_log_warning('The demodulated mode does not return to itself over the '//
     &            'period. The multiplier, the mode and the recorded baseflow are inconsistent, '//
     &            'and the harmonics below are leakage rather than a spectrum.',
     &            this_module, caller)
            end if

            call tfft_close(T)
            deallocate (w, w0)
         end subroutine transform_mode

      !====================================================================
      !     COEFFICIENT OUTPUT
      !====================================================================

         function coef_fname(p, ri, idx) result(fname)
      !! '<p>2<ri><rest of base><idx>.fld', e.g. f2rtorus2Dh001.fld for base
      !! '2dtorus2Dh'; '<p><ri><base><idx>.fld' for a base not starting with '2d'.
            character(len=1), intent(in) :: p
            character(len=1), intent(in) :: ri
            integer, intent(in) :: idx
            character(len=132) :: fname
      ! internal
            character(len=132) :: b
            b = bf_get_base()
            fname = ' '
            if (b(1:2) == '2d') then
               write (fname, '(A,A,A,A,I3.3,A)') p, '2', ri, trim(b(3:)), idx, '.fld'
            else
               write (fname, '(A,A,A,I3.3,A)') p, ri, trim(b), idx, '.fld'
            end if
         end function coef_fname

         subroutine write_coef_files(p, idx, mlist, tags, caller)
      !! Writes Re c_m and Im c_m for m in mlist, in that order, to the two
      !! coefficient files of prefix p and index idx. tags go to the dt slot.
            character(len=1), intent(in) :: p
            integer, intent(in) :: idx
            integer, dimension(:), intent(in) :: mlist
            real(dp), dimension(:), intent(in) :: tags
            character(len=*), intent(in) :: caller
      ! internal
            character(len=256) :: msg
            character(len=132) :: fr, fi
            real(dp), dimension(:, :), allocatable :: rx, ry, rs, ix, iy, iss
            real(dp), dimension(:, :), allocatable :: cre, cim
            real(dp) :: T
            integer :: nrec, nbf, k

            nrec = size(mlist)
            nbf = lx1*ly1*lz1*nelv
            T = tfft_get_period()
            allocate (rx(nbf, nrec), ry(nbf, nrec), rs(nbf, nrec))
            allocate (ix(nbf, nrec), iy(nbf, nrec), iss(nbf, nrec))
            allocate (cre(nbf, 3), cim(nbf, 3))
            do k = 1, nrec
               call tfft_coef(mlist(k), cre, cim)
               call copy(rx(1, k), cre(1, 1), nbf)
               call copy(ry(1, k), cre(1, 2), nbf)
               call copy(rs(1, k), cre(1, 3), nbf)
               call copy(ix(1, k), cim(1, 1), nbf)
               call copy(iy(1, k), cim(1, 2), nbf)
               call copy(iss(1, k), cim(1, 3), nbf)
            end do

            fr = coef_fname(p, 'r', idx)
            fi = coef_fname(p, 'i', idx)
            call bf_write_fields(fr, nrec, T, tags, rx, ry, rs)
            call bf_write_fields(fi, nrec, T, tags, ix, iy, iss)
            write (msg, '(A,A,A,A,A,I0,A,E12.5)') 'wrote ', trim(fr), ' / ', trim(fi),
     &         ': ', nrec, ' records, retained energy= ', sum(tags(1:nrec))
            call nek_log_message(msg, this_module, caller)

            deallocate (rx, ry, rs, ix, iy, iss, cre, cim)
         end subroutine write_coef_files

         subroutine outpost_coefs(tr, ti, mlist)
      !! Nek field files of Re c_m (series tr) and Im c_m (series ti) for m in
      !! mlist. File index = position in mlist, time stamp = m. Field slots
      !! (u_z, u_R, --, --, u_phi) as in tfft_outpost.
            character(len=3), intent(in) :: tr, ti
            integer, dimension(:), intent(in) :: mlist
      ! internal
            real(dp), dimension(:, :), allocatable :: cre, cim
            real(dp), dimension(:), allocatable :: vzero, pzero
            real(dp) :: time_save
            logical :: ifto_save, ifpo_save, ifxyo_save
            integer :: k

            allocate (cre(lv, 3), cim(lv, 3), vzero(lv), pzero(lp2))
            call rzero(vzero, lv)
            call rzero(pzero, lp2)
            time_save = time
            ifto_save = ifto
            ifpo_save = ifpo
            ifxyo_save = ifxyo
            ifto = .true.
            ifpo = .false.
            do k = 1, size(mlist)
               call tfft_coef(mlist(k), cre, cim)
               time = real(mlist(k), dp)
               ifxyo = (k == 1)
               call set_fldindex(tr, k)
               call outpost(cre(1, 1), cre(1, 2), vzero, pzero, cre(1, 3), tr)
               call set_fldindex(ti, k)
               call outpost(cim(1, 1), cim(1, 2), vzero, pzero, cim(1, 3), ti)
            end do
            time = time_save
            ifto = ifto_save
            ifpo = ifpo_save
            ifxyo = ifxyo_save
            deallocate (cre, cim, vzero, pzero)
         end subroutine outpost_coefs

      !====================================================================
      !     PRIVATE HELPERS
      !====================================================================

         subroutine bf_attach(prefix, base, caller)
      !! Ensures the buffer module has its metadata. In-solver (buffer already
      !! live, nsteps_rec > 0): a no-op beyond the recording guard.
      !! Post-processing (fresh job or buffer empty): calls bf_init (a silent
      !! no-op if already initialised) then bf_read_t2Dh with the given prefix
      !! and base, which restores nsteps_rec, nchunk and trec from the sidecar.
      !! The sidecar also carries omega and the forcing dpds; we read them only
      !! to check that t2Dh has been initialised (it must happen separately so
      !! the caller controls which dpds and Wo end up in t2Dh).
            character(len=1), intent(in) :: prefix
            character(len=*), optional, intent(in) :: base
            character(len=*), intent(in) :: caller
      ! internal
            real(dp), dimension(lfc) :: dpds_dummy
            real(dp) :: omega_scar, period_scar
            integer :: nf_scar, ierr
            if (bf_is_recording()) then
               call nek_stop_error('The buffer is still recording. Close the nonlinear pass first.',
     &            this_module, caller)
            end if
            if (bf_get_nsteps() == 0) then
      ! Post-processing path: init buffer (no-op if already done) and read
      ! the sidecar for the metadata (nsteps, nchunk, T). write_chunks must be
      ! .true. or load_chunk refuses to read; nothing is recorded here, so
      ! nothing is written.
               if (.not. present(base)) then
                  call nek_stop_error('Buffer is empty and no base= was given. In post-processing '//
     &               'mode supply base= (e.g. base=''2dtorus2Dh'') and prefix= so the sidecar '//
     &               '<prefix><base>.t2Dh can be located.', this_module, caller)
               end if
               call bf_set_prefix(prefix)
               call bf_init(base=base, write_chunks=.true.)
               call bf_read_t2Dh(dpds_dummy, nf_scar, omega_scar, period_scar, ierr)
               if (ierr /= 0) then
                  call nek_stop_error('Could not read the t2Dh sidecar. Run the Newton solve '//
     &               'first (it writes <prefix><base>.t2Dh), or check the prefix and base '//
     &               'arguments.', this_module, caller)
               end if
               if (bf_get_nsteps() == 0) then
                  call nek_stop_error('Sidecar was read but nsteps is still zero. The sidecar '//
     &               'and the chunk files are inconsistent.', this_module, caller)
               end if
            end if
            if (.not. t2Dh%is_initialised()) then
               call nek_stop_error('t2Dh is not initialised: no omega to project onto. '//
     &            'In post-processing mode call t2Dh%init_flow with the dpds from the sidecar '//
     &            '(use bf_read_t2Dh) before calling the transform wrapper.',
     &            this_module, caller)
            end if
         end subroutine bf_attach

         subroutine pack_base(w, nbf)
      !! (u_z, u_R, u_phi) from the Nek baseflow fields. u_phi lives in the
      !! temperature slot on the 2Dh mesh, which is why the buffer stores t.
            real(dp), dimension(:, :), intent(out) :: w
            integer, intent(in) :: nbf
            call copy(w(1, 1), vx, nbf)
            call copy(w(1, 2), vy, nbf)
            call copy(w(1, 3), t(1, 1, 1, 1, 1), nbf)
         end subroutine pack_base

         subroutine unpack_base(w, nbf)
            real(dp), dimension(:, :), intent(in) :: w
            integer, intent(in) :: nbf
            call copy(vx, w(1, 1), nbf)
            call copy(vy, w(1, 2), nbf)
            call copy(t(1, 1, 1, 1, 1), w(1, 3), nbf)
         end subroutine unpack_base

         subroutine pack_pert(w, nbf, nblk)
      !! Block 1 is the spatial real part (jp = 1), block 2 the imaginary one
      !! (jp = 2). Same ordering the accumulator and tfft_outpost assume.
            real(dp), dimension(:, :), intent(out) :: w
            integer, intent(in) :: nbf, nblk
            integer :: jp, i0
            do jp = 1, nblk
               i0 = 3*(jp - 1)
               call copy(w(1, i0 + 1), vxp(1, jp), nbf)
               call copy(w(1, i0 + 2), vyp(1, jp), nbf)
               call copy(w(1, i0 + 3), tp(1, 1, jp), nbf)
            end do
         end subroutine pack_pert

         real(dp) function fnorm2(w, nbf, nf) result(s)
      !! Mass-weighted square norm summed over all field columns.
            real(dp), dimension(:, :), intent(in) :: w
            integer, intent(in) :: nbf, nf
            real(dp), external :: glsc3
            integer :: i
            s = 0.0_dp
            do i = 1, nf
               s = s + glsc3(w(1, i), bm1, w(1, i), nbf)
            end do
         end function fnorm2

         subroutine w_minus_w0(w, w0, nbf, nf)
      !! w <- w - w0, column by column.
            real(dp), dimension(:, :), intent(inout) :: w
            real(dp), dimension(:, :), intent(in) :: w0
            integer, intent(in) :: nbf, nf
            integer :: i
            do i = 1, nf
               call sub2(w(1, i), w0(1, i), nbf)
            end do
         end subroutine w_minus_w0

      end module t2Dh_tfft_run