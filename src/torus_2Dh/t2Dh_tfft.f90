      module t2Dh_tfft
      !! Temporal Fourier transform of a T-periodic field set, accumulated on
      !! the SAME quadrature and the SAME time grid as the flow-rate
      !! accumulator in t2Dh_flowrate.
      !!
      !!--------------------------------------------------------------------
      !! CONVENTION
      !!--------------------------------------------------------------------
      !!
      !! Native, i.e. identical to neklab_t2Dh:
      !!
      !!   f(x,t) = a_0(x) + sum_m [ a_m(x) cos(m w t) + b_m(x) sin(m w t) ]
      !!
      !!   a_0 = (1/T) int_0^T f dt
      !!   a_m = (2/T) int_0^T f cos(m w t) dt
      !!   b_m = (2/T) int_0^T f sin(m w t) dt
      !!
      !! and the packing of the third index of `fh` is exactly the packing of
      !! dpds / qfour:  j = 1 -> a_0,  j = 2m -> a_m,  j = 2m+1 -> b_m.
      !!
      !! The COMPLEX coefficient handed out by tfft_get is
      !!
      !!   fhat_m = a_m - i b_m      so that   f = Re[ fhat_m exp(i m w t) ]
      !!
      !! which makes |fhat_m| = hypot(a_m,b_m) agree with t2Dh%mf, but
      !!
      !!   arg(fhat_m) = - t2Dh%mf_phase
      !!
      !! because mf_phase is stored in the LAG convention A cos(m w t - phi).
      !! Amplitudes match, phases are negatives of each other. This is the one
      !! sign in the module worth reading twice.
      !!
      !! With two_sided = .true. the result is halved, giving the two-sided
      !! c_m = (1/T) int f exp(-i m w t) dt that an FFT library would return.
      !! The mean is never halved in either convention.
      !!
      !!--------------------------------------------------------------------
      !! QUADRATURE
      !!--------------------------------------------------------------------
      !!
      !! Composite trapezoid on the recorded (CFL-adaptive) grid, O(dt**2).
      !! Relative error on harmonic m is about (m w dt)**2/12, which for a
      !! torus run with O(10**4) steps per period sits many orders below the
      !! Newton residual. It is NOT spectrally accurate: that would need a
      !! uniform grid, and buying it is not worth changing the discrete map
      !! the orbit was converged on.
      !!
      !! A Parseval-style tail is accumulated alongside the harmonics (see
      !! tsq below): sum_i <<f_i^2>_V>_T on the same trapezoidal rule, so
      !! that the truncation loss beyond M can be reported as a fraction of
      !! the time-averaged spatial energy. This is a bilinear consistency
      !! that catches sign errors, packing errors and undersampling in one
      !! number; look at it before believing any per-harmonic amplitude.
      !!
      !!--------------------------------------------------------------------
      !! GENERALITY
      !!--------------------------------------------------------------------
      !!
      !! The accumulator knows nothing about what the nfld field columns mean.
      !! The drivers in t2Dh_tfft_run use
      !!
      !!   nfld = 3 : baseflow          (u_z, u_R, u_phi)
      !!   nfld = 6 : complex mode      (Re u_z, Re u_R, Re u_phi,
      !!                                 Im u_z, Im u_R, Im u_phi)
      !!
      !! with u_phi living at column 3 in the first block and (if nfld = 6)
      !! at column 6 in the second: tfft_spectrum uses this convention to
      !! print the cross-sectional average <u_phi>_A per harmonic, which is
      !! the physically meaningful flow-rate amplitude.
      !!
      !! In the second case the caller must have removed the Floquet
      !! modulation BEFORE calling tfft_add: exp(sigma t) u is not periodic
      !! and its transform is meaningless. See tfft_demodulate.
         use stdlib_optval, only: optval
         use LightKrylov, only: dp, atol_dp
         use neklab_t2Dh, only: t2Dh, t2Dh_basis
         use neklab_nek_setup, only: nek_log_message, nek_log_information,
     &                               nek_log_warning, nek_stop_error,
     &                               set_fldindex
         implicit none
         include "SIZE"
         include "TOTAL"
         private
         character(len=*), parameter, private :: this_module = 't2Dh_tfft'

         integer, parameter :: lv = lx1*ly1*lz1*lelv
         integer, parameter :: lp2 = lx2*ly2*lz2*lelv

      !--------------------------------------------------------------------
      !-----     COMPILE-TIME SIZES                                   -----
      !--------------------------------------------------------------------
         integer, parameter, public :: mfft_max = 32
      !! Largest temporal harmonic the accumulator will hold. Memory is
      !! nbf*nfld*(2M+1)*8 bytes, so M = 32 with nfld = 6 and 64k points is
      !! about 200 MB per rank-share; keep M at what you will actually look at.
         integer, parameter, public :: nfld_max = 6
      !! 3 for a real field set, 6 for a spatially complex one.

      !--------------------------------------------------------------------
      !-----     STATE                                                -----
      !--------------------------------------------------------------------
         real(dp), allocatable :: fh(:, :, :)
      !! (nbf, nfld, 2M+1). Native packing, see the header.
         real(dp), allocatable :: flg(:, :)
      !! (nbf, nfld). Left endpoint of the current trapezoidal interval.
         real(dp) :: omega = 0.0_dp
         real(dp) :: tlagv = 0.0_dp
         real(dp) :: tacc = 0.0_dp
         real(dp) :: period = 0.0_dp
         real(dp) :: tsq = 0.0_dp
      !! Trapezoidal accumulator of sum_i int_0^T (int f_i^2 dV) dt. After
      !! close, tsq/period is the time-averaged spatial energy, and the
      !! Parseval identity says it should equal
      !!   sum_i [ int a_0^i^2 dV + 0.5 sum_{m>=1} (int a_m^i^2 + b_m^i^2) dV ]
      !! for M -> inf. The relative gap is the truncation tail.
         real(dp) :: flg_sq = 0.0_dp
      !! Cached sum_i int flg_i^2 dV, so that tfft_add computes one norm per
      !! step rather than two -- the right endpoint of interval k is the left
      !! endpoint of interval k+1.
         integer :: nbf = 0
         integer :: mharm = 0
         integer :: nfld = 0
         integer :: nacc = 0
         logical :: accumulating = .false.
         logical :: closed = .false.
         character(len=16) :: label = 'field'

         public :: tfft_start, tfft_add, tfft_close, tfft_free
         public :: tfft_get, tfft_demodulate
         public :: tfft_outpost, tfft_spectrum

      contains

      !====================================================================
      !     ACCUMULATION
      !====================================================================

         subroutine tfft_start(mmax, nfields, f0, tag)
      !! Opens an accumulation. f0 is the field set at t = 0: the trapezoidal
      !! rule needs a left endpoint for the first interval, and supplying it
      !! is the caller's job because only the caller knows where t = 0 is.
      !!
      !! For a converged orbit read back from the buffer, u(0) = u(T) and the
      !! last snapshot is the natural left endpoint; the error of that identity
      !! is |F(X)|, i.e. the Newton tolerance.
            integer, intent(in) :: mmax
            integer, intent(in) :: nfields
            real(dp), dimension(:, :), intent(in) :: f0
            character(len=*), optional, intent(in) :: tag
      ! internal
            character(len=*), parameter :: this_procedure = 'tfft_start'
            character(len=256) :: msg
            real(dp), external :: glsc3
            integer :: i
            if (mmax < 0 .or. mmax > mfft_max) then
               write (msg, '(A,I0,A,I0,A)') 'mmax= ', mmax, ' outside [0, ', mfft_max,
     &            ']. Change mfft_max and recompile.'
               call nek_stop_error(msg, this_module, this_procedure)
            end if
            if (nfields < 1 .or. nfields > nfld_max) then
               write (msg, '(A,I0,A,I0)') 'nfields= ', nfields, ' outside [1, ', nfld_max
               call nek_stop_error(msg, this_module, this_procedure)
            end if
            if (.not. t2Dh%is_initialised()) then
               call nek_stop_error('t2Dh is not initialised: no omega to project onto.',
     &            this_module, this_procedure)
            end if
            mharm = mmax
            nfld = nfields
            nbf = lx1*ly1*lz1*nelv
            omega = t2Dh%get_omega()
            if (size(f0, 1) < nbf .or. size(f0, 2) < nfld) then
               call nek_stop_error('f0 is too small for (nbf, nfld).', this_module, this_procedure)
            end if
            call tfft_free()
            allocate (fh(nbf, nfld, 2*mharm + 1), flg(nbf, nfld))
            fh = 0.0_dp
            do i = 1, nfld
               call copy(flg(1, i), f0(1, i), nbf)
            end do
            tlagv = 0.0_dp
            tacc = 0.0_dp
            period = 0.0_dp
            nacc = 0
            tsq = 0.0_dp
            flg_sq = 0.0_dp
            do i = 1, nfld
               flg_sq = flg_sq + glsc3(f0(1, i), bm1, f0(1, i), nbf)
            end do
            accumulating = .true.
            closed = .false.
            label = optval(tag, 'field')
            write (msg, '(A,A,A,I0,A,I0,A,E16.8)') 'open [', trim(label), '] M= ', mharm,
     &         ', nfld= ', nfld, ', omega= ', omega
            call nek_log_information(msg, this_module, this_procedure)
         end subroutine tfft_start

         subroutine tfft_add(f, tval, dtn)
      !! One trapezoidal interval [tlagv, tval]. Call once per step, AFTER the
      !! advance, with tval the time at the END of the step and dtn the step
      !! just taken -- the same arguments accumulate_mflow takes, and
      !! deliberately so: the two transforms then share a grid and a rule.
      !!
      !! Cost is 2*(2M+1)*nfld axpy of length nbf per step, plus one glsc3
      !! per field column for the Parseval sum (the second endpoint is cached
      !! and reused as the next interval's first).
            real(dp), dimension(:, :), intent(in) :: f
            real(dp), intent(in) :: tval
            real(dp), intent(in) :: dtn
      ! internal
            character(len=*), parameter :: this_procedure = 'tfft_add'
            real(dp), external :: glsc3
            real(dp) :: wo, wn, rsq
            integer :: i, j, nb
            if (.not. accumulating) then
               call nek_stop_error('Accumulator is not open. Call tfft_start first.',
     &            this_module, this_procedure)
            end if
            nb = 2*mharm + 1
            do j = 1, nb
               wo = 0.5_dp*dtn*t2Dh_basis(j, omega, tlagv)
               wn = 0.5_dp*dtn*t2Dh_basis(j, omega, tval)
               do i = 1, nfld
                  call add2s2(fh(1, i, j), flg(1, i), wo, nbf)
                  call add2s2(fh(1, i, j), f(1, i), wn, nbf)
               end do
            end do
      ! Parseval tail: sum_i int f_i^2 dV, trapezoidal in time.
            rsq = 0.0_dp
            do i = 1, nfld
               rsq = rsq + glsc3(f(1, i), bm1, f(1, i), nbf)
            end do
            tsq = tsq + 0.5_dp*dtn*(flg_sq + rsq)
            flg_sq = rsq
            do i = 1, nfld
               call copy(flg(1, i), f(1, i), nbf)
            end do
            tlagv = tval
            tacc = tacc + dtn
            nacc = nacc + 1
         end subroutine tfft_add

         subroutine tfft_close(period_in)
      !! Normalises: 1/T on the mean, 2/T on every harmonic. After this the
      !! coefficients are in native convention and tfft_get works.
            real(dp), optional, intent(in) :: period_in
      ! internal
            character(len=*), parameter :: this_procedure = 'tfft_close'
            character(len=256) :: msg
            real(dp) :: T, s
            integer :: i, j
            if (.not. accumulating) then
               call nek_stop_error('Accumulator is not open.', this_module, this_procedure)
            end if
            if (tacc <= 0.0_dp) then
               call nek_stop_error('Empty accumulator.', this_module, this_procedure)
            end if
            T = optval(period_in, tacc)
      ! Same guard as close_mflow: a mismatch here means the sampled interval
      ! is not the period, and every harmonic is projected onto the wrong
      ! basis rather than being merely inaccurate.
            if (present(period_in)) then
               if (abs(tacc - T) > 1.0e-08_dp*max(abs(T), 1.0_dp)) then
                  write (msg, '(A,E16.8,A,E16.8)') 'Accumulated time ', tacc,
     &               ' does not match the period ', T
                  call nek_stop_error(msg, this_module, this_procedure)
               end if
            end if
            do j = 1, 2*mharm + 1
               s = merge(1.0_dp, 2.0_dp, j == 1)/T
               do i = 1, nfld
                  call cmult(fh(1, i, j), s, nbf)
               end do
            end do
            period = T
            accumulating = .false.
            closed = .true.
            write (msg, '(A,A,A,I0,A,E16.8)') 'closed [', trim(label), '] over ', nacc,
     &         ' steps, T= ', period
            call nek_log_information(msg, this_module, this_procedure)
            if (nacc < 8*max(mharm, 1)) then
               call nek_log_warning('Fewer than 8 samples per period of the highest harmonic: '//
     &            'the top of the spectrum is not resolved.', this_module, this_procedure)
            end if
         end subroutine tfft_close

         subroutine tfft_free()
            if (allocated(fh)) deallocate (fh)
            if (allocated(flg)) deallocate (flg)
            tsq = 0.0_dp
            flg_sq = 0.0_dp
            accumulating = .false.
            closed = .false.
         end subroutine tfft_free

      !====================================================================
      !     FLOQUET DEMODULATION
      !====================================================================

         subroutine tfft_demodulate(f, tval, T, mu_abs, mu_arg)
      !! Removes the Floquet modulation IN PLACE, so that what is handed to
      !! tfft_add is the T-periodic part of the mode.
      !!
      !! A Floquet solution is u(x,t) = exp(sigma t) utilde(x,t) with utilde
      !! T-periodic and mu = exp(sigma T) the multiplier. The raw perturbation
      !! is therefore NOT periodic -- it is multiplied by mu over one period --
      !! and projecting it onto exp(i m w t) returns leakage, not harmonics.
      !! What is periodic is
      !!
      !!    utilde(t) = u(t) * mu**(-t/T)
      !!              = u(t) * |mu|**(-t/T) * exp(-i arg(mu) t/T)
      !!
      !! and since mu is complex in general the phase rotation MIXES the two
      !! spatial blocks: this cannot be done one block at a time.
      !!
      !! f(:,1:3) is the spatial real part, f(:,4:6) the spatial imaginary
      !! part. With nfld = 3 the mode is real (alpha = 0) and only the
      !! amplitude factor applies; a negative real mu then describes a
      !! period-doubled mode which a single perturbation slot cannot carry,
      !! and that case is rejected rather than silently mangled.
            real(dp), dimension(:, :), intent(inout) :: f
            real(dp), intent(in) :: tval
            real(dp), intent(in) :: T
            real(dp), intent(in) :: mu_abs
            real(dp), intent(in) :: mu_arg
      ! internal
            character(len=*), parameter :: this_procedure = 'tfft_demodulate'
            real(dp) :: s, ang, c, sn, re, im
            integer :: i, k, ntot
            ntot = min(size(f, 1), lx1*ly1*lz1*nelv)
            if (mu_abs <= 0.0_dp) then
               call nek_stop_error('|mu| must be positive.', this_module, this_procedure)
            end if
            if (T <= 0.0_dp) then
               call nek_stop_error('Period must be positive.', this_module, this_procedure)
            end if
            s = mu_abs**(-tval/T)
            ang = mu_arg*tval/T
            c = cos(ang)
            sn = sin(ang)
            if (size(f, 2) < 6) then
               if (abs(mu_arg) > atol_dp) then
                  call nek_stop_error('A complex mu needs both perturbation slots (alpha /= 0). '//
     &               'A real mode cannot represent a rotating or period-doubled Floquet mode.',
     &               this_module, this_procedure)
               end if
               do i = 1, size(f, 2)
                  call cmult(f(1, i), s, ntot)
               end do
            else
      ! (re + i im) * s * (cos ang - i sin ang)
               do i = 1, 3
                  do k = 1, ntot
                     re = f(k, i)
                     im = f(k, i + 3)
                     f(k, i) = s*(re*c + im*sn)
                     f(k, i + 3) = s*(im*c - re*sn)
                  end do
               end do
            end if
         end subroutine tfft_demodulate

      !====================================================================
      !     EXTRACTION
      !====================================================================

         subroutine tfft_get(m, fre, fim, two_sided)
      !! Complex coefficient of harmonic m as a pair of real field sets:
      !!
      !!    f(x,t) = Re[ (fre + i fim) exp(i m w t) ]
      !!
      !! so |fre + i fim| matches t2Dh%mf and arg(fre + i fim) is MINUS
      !! t2Dh%mf_phase. See the module header.
      !!
      !! For a spatially complex field set (nfld = 6) this is applied to both
      !! blocks independently, which gives the four real objects
      !!
      !!    A_m = fre(:,1:3) + i fre(:,4:6)     coefficient of cos(m w t)
      !!   -B_m = fim(:,1:3) + i fim(:,4:6)     minus that of sin(m w t)
      !!
      !! and the two-sided complex spectrum of the complex signal is
      !! c_{+-m} = (A_m -+ i B_m)/2. Note c_{-m} is NOT the conjugate of c_m
      !! for a complex signal, so both signs carry information.
            integer, intent(in) :: m
            real(dp), dimension(:, :), intent(out) :: fre
            real(dp), dimension(:, :), intent(out) :: fim
            logical, optional, intent(in) :: two_sided
      ! internal
            character(len=*), parameter :: this_procedure = 'tfft_get'
            character(len=256) :: msg
            real(dp) :: s
            integer :: i
            call check_closed(this_procedure)
            if (m < 0 .or. m > mharm) then
               write (msg, '(A,I0,A,I0,A)') 'Harmonic ', m, ' outside [0, ', mharm, '].'
               call nek_stop_error(msg, this_module, this_procedure)
            end if
            if (size(fre, 1) < nbf .or. size(fre, 2) < nfld .or.
     &          size(fim, 1) < nbf .or. size(fim, 2) < nfld) then
               call nek_stop_error('Output arrays are too small.', this_module, this_procedure)
            end if
            s = merge(0.5_dp, 1.0_dp, optval(two_sided, .false.))
            do i = 1, nfld
               if (m == 0) then
      ! The mean is real, and it is not halved in either convention.
                  call copy(fre(1, i), fh(1, i, 1), nbf)
                  call rzero(fim(1, i), nbf)
               else
                  call copy(fre(1, i), fh(1, i, 2*m), nbf)
                  call copy(fim(1, i), fh(1, i, 2*m + 1), nbf)
                  call chsign(fim(1, i), nbf)
                  call cmult(fre(1, i), s, nbf)
                  call cmult(fim(1, i), s, nbf)
               end if
            end do
         end subroutine tfft_get

      !====================================================================
      !     DIAGNOSTICS AND OUTPUT
      !====================================================================

         subroutine tfft_spectrum(amp, filename)
      !! For every harmonic m and every field column i, computes the spatial
      !! L2 amplitude
      !!
      !!   amp(m,i)^2 = int_V a_m^i(x)^2 dV + int_V b_m^i(x)^2 dV      (m>=1)
      !!   amp(0,i)^2 = int_V a_0^i(x)^2 dV                            (m=0)
      !!
      !! using Nek's mass matrix bm1. The physically interesting quantity for
      !! the torus / pulsatile run is the CROSS-SECTIONAL AVERAGE of the
      !! azimuthal velocity, which is the flow-rate amplitude at m; it is
      !! computed with t2Dh%ubar_arr (the same routine that measures the
      !! mean flow rate elsewhere) and logged per harmonic. For nfld = 6 the
      !! two spatial blocks of the mode are printed side by side.
      !!
      !! A Parseval consistency line is also logged: with the same measure,
      !!
      !!   E_time = (1/T) int_0^T sum_i int f_i^2 dV dt   (accumulated in tsq)
      !!   E_spec = sum_i [ amp(0,i)^2 + 0.5 sum_{m>=1} amp(m,i)^2 ]
      !!
      !! agree in the M -> infinity limit. The tail (E_time - E_spec)/E_time
      !! is the fraction of energy the accumulator threw away, which is the
      !! honest measure of whether M is big enough.
      !!
      !! If filename is present, a two-block text file is written by rank 0:
      !! one block of absolute amplitudes, one block normalised by the bulk
      !! mean velocity U0 = <u_phi>_A(m=0) (first phi block).
            real(dp), dimension(0:, :), intent(out), optional :: amp
            character(len=*), intent(in), optional :: filename
      ! internal
            character(len=*), parameter :: this_procedure = 'tfft_spectrum'
            character(len=256) :: msg
            real(dp), dimension(0:mharm, nfld) :: amp_
            real(dp), dimension(0:mharm) :: E_frac
            real(dp), dimension(0:mharm, max(nfld/3, 1)) :: mean_c, mean_s
            real(dp), external :: glsc3
            real(dp) :: E_spec, E_time, tail, U0
            integer :: m, i, ib, nblk, iphi

            call check_closed(this_procedure)
            nblk = nfld/3

      ! --- per-column L2 amplitudes
            do i = 1, nfld
               amp_(0, i) = sqrt(max(glsc3(fh(1, i, 1), bm1, fh(1, i, 1), nbf), 0.0_dp))
               do m = 1, mharm
                  amp_(m, i) = sqrt(max(
     &               glsc3(fh(1, i, 2*m), bm1, fh(1, i, 2*m), nbf)
     &             + glsc3(fh(1, i, 2*m + 1), bm1, fh(1, i, 2*m + 1), nbf), 0.0_dp))
               end do
            end do

      ! --- bulk cross-sectional averages of the u_phi block(s). Native pack
      ! places u_phi at column 3*ib for ib = 1..nblk. c_m = <a_m>_A -i<b_m>_A.
            do ib = 1, nblk
               iphi = 3*ib
               mean_c(0, ib) = t2Dh%ubar_arr(fh(1, iphi, 1))
               mean_s(0, ib) = 0.0_dp
               do m = 1, mharm
                  mean_c(m, ib) = t2Dh%ubar_arr(fh(1, iphi, 2*m))
                  mean_s(m, ib) = t2Dh%ubar_arr(fh(1, iphi, 2*m + 1))
               end do
            end do

      ! --- cumulative energy fraction: sum_{k<=m} energy(k) / E_spec
            E_spec = 0.0_dp
            do i = 1, nfld
               E_spec = E_spec + amp_(0, i)**2
            end do
            E_frac(0) = E_spec
            do m = 1, mharm
               E_spec = 0.0_dp
               do i = 1, nfld
                  E_spec = E_spec + 0.5_dp*amp_(m, i)**2
               end do
               E_frac(m) = E_spec
            end do
            E_spec = 0.0_dp
            do m = 0, mharm
               E_spec = E_spec + E_frac(m)
            end do
            if (E_spec > atol_dp) then
               E_frac = E_frac/E_spec
            else
               E_frac = 0.0_dp
            end if

      ! --- Parseval tail
            E_time = tsq/max(period, atol_dp)
            tail = 1.0_dp - E_spec/max(E_time, atol_dp)

      ! --- reference velocity: bulk mean of u_phi at m=0 (first phi block).
      ! For the perturbation this may be near zero, in which case the
      ! normalised block below is meaningless -- the value of U0 is logged
      ! so the user can see it.
            U0 = mean_c(0, 1)

      ! --- log
            write (msg, '(A,A,A)') 'temporal spectrum [', trim(label), ']:'
            call nek_log_message(msg, this_module, this_procedure)
            write (msg, '(5X,A,E16.8,A,E16.8,A,E10.3)')
     &         'Parseval  E_spec= ', E_spec, ' E_time= ', E_time, ' tail= ', tail
            call nek_log_message(msg, this_module, this_procedure)
            write (msg, '(5X,A,E16.8)') 'U0 = <u_phi>_A(m=0) = ', U0
            call nek_log_message(msg, this_module, this_procedure)
            if (nblk == 1) then
               call nek_log_message(
     &            '     m       |<u_phi>|_A     phase[rad]     E rel',
     &            this_module, this_procedure)
               do m = 0, mharm
                  write (msg, '(3X,I3,3(1X,E14.6))') m,
     &               hypot(mean_c(m, 1), mean_s(m, 1)),
     &               atan2(-mean_s(m, 1), mean_c(m, 1)),
     &               E_frac(m)
                  call nek_log_message(msg, this_module, this_procedure)
               end do
            else
               call nek_log_message(
     &            '     m    |<uphi^Re>|_A   phase^Re      |<uphi^Im>|_A   phase^Im      E rel',
     &            this_module, this_procedure)
               do m = 0, mharm
                  write (msg, '(3X,I3,5(1X,E14.6))') m,
     &               hypot(mean_c(m, 1), mean_s(m, 1)),
     &               atan2(-mean_s(m, 1), mean_c(m, 1)),
     &               hypot(mean_c(m, 2), mean_s(m, 2)),
     &               atan2(-mean_s(m, 2), mean_c(m, 2)),
     &               E_frac(m)
                  call nek_log_message(msg, this_module, this_procedure)
               end do
            end if
            if (tail > 1.0e-03_dp) then
               call nek_log_warning('Spectral tail > 0.1%: energy beyond M is not negligible. '//
     &            'Increase mmax or expect truncation error at that level.',
     &            this_module, this_procedure)
            end if

      ! --- optional outputs
            if (present(amp)) then
               if (size(amp, 1) < mharm + 1 .or. size(amp, 2) < nfld) then
                  call nek_stop_error('amp array too small for (0:mharm, 1:nfld).',
     &               this_module, this_procedure)
               end if
               amp(0:mharm, 1:nfld) = amp_
            end if

            if (present(filename)) then
               call write_spectrum_file(filename, amp_, U0, E_frac)
            end if
         end subroutine tfft_spectrum

         subroutine tfft_outpost(tags, mlist, if_mesh)
      !! ONE file series per (spatial block, real/imaginary part). The harmonic
      !! is the FILE NUMBER inside the series, not part of the name:
      !!
      !!   nfld = 3   r__  Re fhat_m        i__  Im fhat_m
      !!   nfld = 6   rRe  Re_t of Re_x     iRe  Im_t of Re_x
      !!              rIm  Re_t of Im_x     iIm  Im_t of Im_x
      !!
      !! so r__2dtorus0.f00001 is the mean and r__2dtorus0.f00017 is m = 16,
      !! and a post-processor opens the series once and steps through m.
      !! The harmonic also rides on the time stamp, so m is the time axis.
      !!
      !! Field slots inside each file are (u_z, u_R, --, --, u_phi): the third
      !! velocity component and the pressure are written as zeros because the
      !! 2Dh mesh has no third direction and the buffer stores no pressure.
            character(len=3), dimension(:), optional, intent(in) :: tags
      !! 2*nblk series names, overriding the defaults above.
            integer, dimension(:), optional, intent(in) :: mlist
            logical, optional, intent(in) :: if_mesh
      !! Write the geometry into the first file of each series. Default .true.
      ! internal
            character(len=*), parameter :: this_procedure = 'tfft_outpost'
            character(len=256) :: msg
            character(len=3), dimension(4) :: name
            real(dp), dimension(:, :), allocatable :: fre, fim
            real(dp), dimension(:), allocatable :: vzero, pzero
            real(dp) :: time_save
            logical :: ifto_save, ifpo_save, ifxyo_save, mesh_
            integer :: m, mm, nblk, nser, ib, i0, nout

            call check_closed(this_procedure)
            nblk = nfld/3
            if (nblk*3 /= nfld) then
               call nek_stop_error('tfft_outpost expects nfld to be a multiple of 3.',
     &            this_module, this_procedure)
            end if
            nser = 2*nblk
            mesh_ = optval(if_mesh, .true.)

            if (present(tags)) then
               if (size(tags) < nser) then
                  write (msg, '(A,I0,A)') 'tags must supply ', nser, ' names.'
                  call nek_stop_error(msg, this_module, this_procedure)
               end if
               do ib = 1, nser
                  name(ib) = tags(ib)
               end do
            else if (nblk == 1) then
               name(1) = 'r__'
               name(2) = 'i__'
            else
               name(1) = 'rRe'
               name(2) = 'iRe'
               name(3) = 'rIm'
               name(4) = 'iIm'
            end if

            allocate (fre(lv, nfld), fim(lv, nfld), vzero(lv), pzero(lp2))
            call rzero(vzero, lv)
            call rzero(pzero, lp2)
            time_save = time
            ifto_save = ifto
            ifpo_save = ifpo
            ifxyo_save = ifxyo
            ifto = .true.
            ifpo = .false.

            if (present(mlist)) then
               nout = size(mlist)
            else
               nout = mharm + 1
            end if

            do mm = 1, nout
               if (present(mlist)) then
                  m = mlist(mm)
               else
                  m = mm - 1
               end if
               if (m < 0 .or. m > mharm) cycle
               call tfft_get(m, fre, fim)
               time = real(m, dp)
      ! Geometry in the first file of each series only: it is the same mesh
      ! every time, and 2*nblk copies of it is already one too many.
               ifxyo = mesh_ .and. (mm == 1)
               do ib = 1, nblk
                  i0 = 3*(ib - 1)
      ! Name fixed, harmonic carried by the file index: m -> f<m+1>.
                  call set_fldindex(name(2*ib - 1), m + 1)
                  call outpost(fre(1, i0 + 1), fre(1, i0 + 2), vzero, pzero,
     &                         fre(1, i0 + 3), name(2*ib - 1))
                  call set_fldindex(name(2*ib), m + 1)
                  call outpost(fim(1, i0 + 1), fim(1, i0 + 2), vzero, pzero,
     &                         fim(1, i0 + 3), name(2*ib))
               end do
            end do

            time = time_save
            ifto = ifto_save
            ifpo = ifpo_save
            ifxyo = ifxyo_save
            deallocate (fre, fim, vzero, pzero)
            write (msg, '(A,I0,A,I0,A,A,A)') 'wrote ', nser, ' series x ', nout,
     &         ' harmonics for [', trim(label), ']'
            call nek_log_message(msg, this_module, this_procedure)
         end subroutine tfft_outpost

      !====================================================================
      !     PRIVATE HELPERS
      !====================================================================

         subroutine check_closed(caller)
            character(len=*), intent(in) :: caller
            if (.not. closed) then
               call nek_stop_error('No closed transform available. Run a driver first.',
     &            this_module, caller)
            end if
         end subroutine check_closed

         subroutine write_spectrum_file(filename, amp_, U0, E_frac)
      !! Two-block text file: one block absolute, one block normalised by U0.
      !! Header carries T, omega, M, nfld, label and U0 so plots need no
      !! external context. Rank 0 only.
            character(len=*), intent(in) :: filename
            real(dp), intent(in) :: amp_(0:, :)
            real(dp), intent(in) :: U0
            real(dp), intent(in) :: E_frac(0:)
      ! internal
            character(len=*), parameter :: this_procedure = 'write_spectrum_file'
            character(len=256) :: msg
            real(dp) :: denom
            integer :: iunit, ierr, m, i

            if (nid /= 0) return
            open (newunit=iunit, file=filename, status='replace', action='write',
     &            iostat=ierr)
            if (ierr /= 0) then
               write (msg, '(A,A,A)') 'Could not open ', trim(filename),
     &            ' for writing; skipping.'
               call nek_log_warning(msg, this_module, this_procedure)
               return
            end if

            write (iunit, '(A,A,A,I0,A,I0,A,E16.8,A,E16.8)')
     &         '# t2Dh spectrum  label=', trim(label), '   nfld=', nfld,
     &         '   M=', mharm, '   T=', period, '   omega=', omega
            write (iunit, '(A,E16.8)') '# U0 = <u_phi>_A(m=0) = ', U0
            write (iunit, '(A)') '# Amplitude convention:'
            write (iunit, '(A)') '#   A_m^(i) = sqrt( int_V a_m^i(x)^2 dV + int_V b_m^i(x)^2 dV )   for m>=1'
            write (iunit, '(A)') '#   A_0^(i) = sqrt( int_V a_0^i(x)^2 dV )'
            write (iunit, '(A)') '#'

            write (iunit, '(A)') '# BLOCK 1 -- ABSOLUTE AMPLITUDES'
            call write_header(iunit, .false.)
            do m = 0, mharm
               write (iunit, '(I5,1X,10(E16.8,1X))') m,
     &            (amp_(m, i), i=1, nfld), E_frac(m)
            end do

            write (iunit, '(A)') ''
            write (iunit, '(A)') '# BLOCK 2 -- NORMALISED BY U0'
            call write_header(iunit, .true.)
            denom = max(abs(U0), atol_dp)
            do m = 0, mharm
               write (iunit, '(I5,1X,10(E16.8,1X))') m,
     &            (amp_(m, i)/denom, i=1, nfld), E_frac(m)
            end do

            close (iunit)
            write (msg, '(A,A,A,A,A)') 'wrote spectrum file ', trim(filename),
     &         ' for [', trim(label), ']'
            call nek_log_message(msg, this_module, this_procedure)
         end subroutine write_spectrum_file

         subroutine write_header(iunit, normalised)
      !! Column header line for the spectrum file, matching the block below.
            integer, intent(in) :: iunit
            logical, intent(in) :: normalised
            character(len=8) :: suf
            suf = merge('/U0     ', '        ', normalised)
            if (nfld == 3) then
               write (iunit, '(A,A,A,A,A,A,A)')
     &            '#   m      A(u_z)', trim(suf),
     &            '   A(u_R)', trim(suf),
     &            '   A(u_phi)', trim(suf),
     &            '   E rel'
            else
               write (iunit, '(A,A,A,A,A,A,A,A,A,A,A,A,A)')
     &            '#   m      A(u_z^Re)', trim(suf),
     &            '   A(u_R^Re)', trim(suf),
     &            '   A(u_phi^Re)', trim(suf),
     &            '   A(u_z^Im)', trim(suf),
     &            '   A(u_R^Im)', trim(suf),
     &            '   A(u_phi^Im)', trim(suf),
     &            '   E rel'
            end if
         end subroutine write_header

      end module t2Dh_tfft