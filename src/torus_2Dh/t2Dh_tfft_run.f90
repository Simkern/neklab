      module t2Dh_tfft_run
      !! The two userchk-facing wrappers around the accumulator in t2Dh_tfft.
      !!
      !!   tfft_baseflow     : sweeps the recorded orbit and transforms it.
      !!                       No solves, no time stepping -- the buffer is a
      !!                       complete record of the trajectory, so this is
      !!                       one pass of file reads and axpys.
      !!
      !!   tfft_perturbation : re-runs the linear solver over one period on the
      !!                       replayed baseflow and transforms the mode on the
      !!                       fly, because the perturbation fields are not
      !!                       stored anywhere.
      !!
      !! Both are meant to be called from userchk AFTER the Newton solve (and,
      !! for the second, after the stability analysis has produced a mode and
      !! its multiplier). Both leave the Nek state as they found it.
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
      !! difference between this wrapper and the baseflow one.
      !!
      !! Pass the multiplier the eigensolver returned. If you pass mu = 1 the
      !! wrapper transforms the raw response, which is the right thing only
      !! for a periodically forced linear problem with no homogeneous growth.
      !!
      !! The wrapper reports |utilde(T) - utilde(0)|/|utilde(0)| at the end.
      !! That number is a direct check on mu: for a converged Floquet pair it
      !! sits at the eigensolver's tolerance, and if it does not, the mode, the
      !! multiplier or the baseflow record do not belong together and no
      !! spectrum computed from them means anything.
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
      !! to the recorded value. Check this against floquet_operator_torus.f90
      !! and align the setup_linear_solver arguments with it: if the propagator
      !! configures the solver differently, the mode this wrapper marches is
      !! not the mode the eigensolver found.
         use stdlib_optval, only: optval
         use LightKrylov, only: dp, atol_dp
         use neklab_vectors, only: nek_zvector
         use neklab_utils, only: vec2nek
         use neklab_t2Dh, only: t2Dh
         use neklab_2Dh, only: nek_advance_2Dh
         use neklab_2Dh_axisym, only: nek_advance_2Dh_axisym
         use neklab_nek_setup, only: setup_linear_solver, nek_log_message,
     &                               nek_log_information, nek_log_warning,
     &                               nek_log_debug, nek_stop_error
         use t2Dh_bf_buffer, only: bf_get_nsteps, bf_get_time, bf_get_dt,
     &                             bf_set, bf_set_prefix, bf_get_prefix,
     &                             bf_replay_start, bf_replay_end, bf_is_recording
         use t2Dh_tfft, only: tfft_start, tfft_add, tfft_close, tfft_free,
     &                        tfft_demodulate, tfft_outpost, tfft_spectrum,
     &                        tfft_check_ubar, tfft_mharm, mfft_max
         implicit none
         include "SIZE"
         include "TOTAL"
         private
         character(len=*), parameter, private :: this_module = 't2Dh_tfft_run'

         integer, parameter :: lv = lx1*ly1*lz1*lelv

         public :: tfft_baseflow, tfft_perturbation

      contains

      !====================================================================
      !     BASEFLOW
      !====================================================================

         subroutine tfft_baseflow(mmax, prefix, if_outpost, if_check)
      !! Temporal Fourier transform of the recorded orbit.
      !!
      !! Reads the chunk files under `prefix` (default 'b', the converged orbit
      !! the driver re-records for the Floquet run) and projects every snapshot
      !! onto the temporal basis. Nothing is integrated: the recording already
      !! holds u(x,t_k) and dt_k for every step, which is exactly and only what
      !! the quadrature needs.
      !!
      !! Because the transform is a post-processing sweep rather than something
      !! bolted into the time loop, the number of harmonics is a choice you
      !! make here and can change without re-running anything.
            integer, intent(in) :: mmax
            character(len=1), optional, intent(in) :: prefix
            logical, optional, intent(in) :: if_outpost
            logical, optional, intent(in) :: if_check
      ! internal
            character(len=*), parameter :: this_procedure = 'tfft_baseflow'
            character(len=256) :: msg
            real(dp), dimension(:, :), allocatable :: w, vsave
            real(dp), dimension(0:mfft_max) :: spec
            character(len=1) :: pfx_save
            real(dp) :: tk, dtk, T
            integer :: k, n, nbf

            call guard_buffer(this_procedure)
            n = bf_get_nsteps()
            T = bf_get_time()
            nbf = lx1*ly1*lz1*nelv
            write (msg, '(A,I0,A,E16.8,A,A,A)') 'transforming ', n, ' snapshots, T= ', T,
     &         ', prefix ''', optval(prefix, 'b'), ''''
            call nek_log_message(msg, this_module, this_procedure)

            allocate (w(lv, 3), vsave(lv, 3))
      ! The sweep overwrites vx/vy/t with recorded snapshots, so put back what
      ! the caller had. userchk is not a place to leave the state disturbed.
            call pack_base(vsave, nbf)

            pfx_save = bf_get_prefix()
            call bf_set_prefix(optval(prefix, 'b'))
            call bf_replay_start()

      ! Left endpoint of the first interval. For a converged orbit u(0) = u(T),
      ! so the last snapshot is it; the error of that identity is |F(X)|.
            call bf_set(n, if_lag=.false.)
            call pack_base(w, nbf)
            call tfft_start(mmax, 3, w, 'baseflow')

            tk = 0.0_dp
            do k = 1, n
               dtk = bf_get_dt(k)
               call bf_set(k, if_lag=.false.)
               tk = tk + dtk
               call pack_base(w, nbf)
               call tfft_add(w, tk, dtk)
            end do

            call bf_replay_end()
            call bf_set_prefix(pfx_save)
            call tfft_close(T)

            call unpack_base(vsave, nbf)
            deallocate (w, vsave)

            call tfft_spectrum(spec)
      ! The bulk average of the field harmonics IS the flow-rate harmonic.
      ! Cheap, and it catches every sign and factor-of-two in the chain.
            if (optval(if_check, .true.)) call tfft_check_ubar(3)
            if (optval(if_outpost, .true.)) call tfft_outpost('ri')
         end subroutine tfft_baseflow

      !====================================================================
      !     PERTURBATION
      !====================================================================

         subroutine tfft_perturbation(mmax, zmode, alpha, mu, prefix, vtol, ptol,
     &                                if_outpost, if_normalise, if_axisym)
      !! Re-runs the linear solver over one period on the replayed baseflow and
      !! transforms the Floquet mode on the fly.
      !!
      !! zmode is the eigenvector as the stability analysis returned it, mu its
      !! multiplier. alpha is the azimuthal wavenumber: alpha = 0 needs npert=1
      !! and gives a real mode (3 field columns), alpha /= 0 needs npert=2 and
      !! gives a spatially complex one (6 columns, Re block then Im block).
            integer, intent(in) :: mmax
            type(nek_zvector), intent(in) :: zmode
            real(dp), intent(in) :: alpha
            complex(dp), optional, intent(in) :: mu
      !! Floquet multiplier. Default (1,0), i.e. no demodulation, which is only
      !! right for a forced linear response with no homogeneous growth.
            character(len=1), optional, intent(in) :: prefix
            real(dp), optional, intent(in) :: vtol
            real(dp), optional, intent(in) :: ptol
            logical, optional, intent(in) :: if_outpost
            logical, optional, intent(in) :: if_normalise
            logical, optional, intent(in) :: if_axisym
      ! internal
            character(len=*), parameter :: this_procedure = 'tfft_perturbation'
            character(len=256) :: msg
            real(dp), dimension(:, :), allocatable :: w, w0
            real(dp), dimension(0:mfft_max) :: spec
            complex(dp) :: mu_
            character(len=1) :: pfx_save
            real(dp) :: T, dtk, mu_abs, mu_arg, nrm0, nrmT, dev, scl
            integer :: k, n, nbf, nblk, nf, i
            logical :: axisym_

            call guard_buffer(this_procedure)
            n = bf_get_nsteps()
            T = bf_get_time()
            nbf = lx1*ly1*lz1*nelv
            axisym_ = optval(if_axisym, .true.)
            mu_ = optval(mu, cmplx(1.0_dp, 0.0_dp, kind=dp))
            mu_abs = abs(mu_)
            mu_arg = atan2(aimag(mu_), real(mu_, kind=dp))
            if (mu_abs <= 0.0_dp) then
               call nek_stop_error('|mu| must be positive.', this_module, this_procedure)
            end if

      ! --- geometry of the mode
            if (alpha == 0.0_dp) then
               nblk = 1
               if (lpert /= 1) call nek_stop_error('alpha = 0 requires lpert = 1.',
     &            this_module, this_procedure)
            else
               nblk = 2
               if (lpert /= 2) call nek_stop_error('alpha /= 0 requires lpert = 2.',
     &            this_module, this_procedure)
            end if
            nf = 3*nblk

            write (msg, '(A,I0,A,E16.8,A,E16.8,A,E16.8)') 'linear replay over ', n,
     &         ' steps, alpha= ', alpha, ', |mu|= ', mu_abs, ', arg(mu)= ', mu_arg
            call nek_log_message(msg, this_module, this_procedure)
            if (abs(mu_abs - 1.0_dp) < atol_dp .and. abs(mu_arg) < atol_dp) then
               call nek_log_warning('mu = 1: no demodulation will be applied. That is correct '//
     &            'only for a forced response with no homogeneous growth, not for a Floquet mode.',
     &            this_module, this_procedure)
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
               call nek_stop_error('The mode has zero norm.', this_module, this_procedure)
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
            call bf_set_prefix(optval(prefix, 'b'))
            call bf_replay_start()

            time = 0.0_dp
            istep = 0
            lastep = 0
            call tfft_start(mmax, nf, w, 'mode')

            do k = 1, n
               istep = istep + 1
      ! if_lag defaults to .true.: the linear solve needs the baseflow lag
      ! levels, and bf_set builds them from the traversal order.
               call bf_set(k)
               if (axisym_) then
                  call nek_advance_2Dh_axisym(alpha)
               else
                  call nek_advance_2Dh(alpha)
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
            call nek_log_message(msg, this_module, this_procedure)
            if (dev > 1.0e-04_dp) then
               call nek_log_warning('The demodulated mode does not return to itself over the '//
     &            'period. The multiplier, the mode and the recorded baseflow are inconsistent, '//
     &            'and the harmonics below are leakage rather than a spectrum.',
     &            this_module, this_procedure)
            end if

            call tfft_close(T)
            deallocate (w, w0)

            call tfft_spectrum(spec)
            if (optval(if_outpost, .true.)) then
               if (nblk == 1) then
                  call tfft_outpost('ri')
               else
                  call tfft_outpost('abcd')
               end if
            end if
         end subroutine tfft_perturbation

      !====================================================================
      !     PRIVATE HELPERS
      !====================================================================

         subroutine guard_buffer(caller)
            character(len=*), intent(in) :: caller
            if (bf_is_recording()) then
               call nek_stop_error('The buffer is still recording. Close the nonlinear pass first.',
     &            this_module, caller)
            end if
            if (bf_get_nsteps() == 0) then
               call nek_stop_error('Nothing recorded: run the Newton solve (or attach an '//
     &            'existing recording) before asking for its transform.', this_module, caller)
            end if
            if (.not. t2Dh%is_initialised()) then
               call nek_stop_error('t2Dh is not initialised: no omega to project onto.',
     &            this_module, caller)
            end if
         end subroutine guard_buffer

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
