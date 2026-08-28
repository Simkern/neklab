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
      !! The same rule applied to a linear functional of the same integrand
      !! must reproduce the flow-rate harmonics exactly, which is what
      !! tfft_check_ubar tests.
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
      !! and in the second case the caller must have removed the Floquet
      !! modulation BEFORE calling tfft_add: exp(sigma t) u is not periodic
      !! and its transform is meaningless. See tfft_demodulate.
         use stdlib_optval, only: optval
         use LightKrylov, only: dp, atol_dp
         use neklab_t2Dh, only: t2Dh, t2Dh_basis, lfc
         use neklab_nek_setup, only: nek_log_message, nek_log_information,
     &                               nek_log_warning, nek_log_debug, nek_stop_error,
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
         integer :: nbf = 0
         integer :: nb = 0
         integer :: mharm = 0
         integer :: nfld = 0
         integer :: nacc = 0
         logical :: accumulating = .false.
         logical :: closed = .false.
         character(len=16) :: label = 'field'

         public :: tfft_start, tfft_add, tfft_close, tfft_free
         public :: tfft_get, tfft_get_raw, tfft_demodulate
         public :: tfft_outpost, tfft_spectrum, tfft_check_ubar
         public :: tfft_mharm, tfft_nfld, tfft_nbf, tfft_is_closed, tfft_omega

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
            nb = 2*mmax + 1
            nfld = nfields
            nbf = lx1*ly1*lz1*nelv
            omega = t2Dh%get_omega()
            if (size(f0, 1) < nbf .or. size(f0, 2) < nfld) then
               call nek_stop_error('f0 is too small for (nbf, nfld).', this_module, this_procedure)
            end if
            call tfft_free()
            allocate (fh(nbf, nfld, nb), flg(nbf, nfld))
            fh = 0.0_dp
            do i = 1, nfld
               call copy(flg(1, i), f0(1, i), nbf)
            end do
            tlagv = 0.0_dp
            tacc = 0.0_dp
            period = 0.0_dp
            nacc = 0
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
      !! Cost is 2*(2M+1)*nfld axpy of length nbf per step. No communication:
      !! every operation is pointwise, nothing is reduced.
            real(dp), dimension(:, :), intent(in) :: f
            real(dp), intent(in) :: tval
            real(dp), intent(in) :: dtn
      ! internal
            character(len=*), parameter :: this_procedure = 'tfft_add'
            real(dp) :: wo, wn
            integer :: i, j
            if (.not. accumulating) then
               call nek_stop_error('Accumulator is not open. Call tfft_start first.',
     &            this_module, this_procedure)
            end if
            do j = 1, nb
               wo = 0.5_dp*dtn*t2Dh_basis(j, omega, tlagv)
               wn = 0.5_dp*dtn*t2Dh_basis(j, omega, tval)
               do i = 1, nfld
                  call add2s2(fh(1, i, j), flg(1, i), wo, nbf)
                  call add2s2(fh(1, i, j), f(1, i), wn, nbf)
               end do
            end do
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
            do j = 1, nb
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

         subroutine tfft_get_raw(j, f)
      !! Raw native slot j: 1 -> a_0, 2m -> a_m, 2m+1 -> b_m. For anyone who
      !! would rather do their own bookkeeping than trust the sign in tfft_get.
            integer, intent(in) :: j
            real(dp), dimension(:, :), intent(out) :: f
      ! internal
            character(len=*), parameter :: this_procedure = 'tfft_get_raw'
            integer :: i
            call check_closed(this_procedure)
            if (j < 1 .or. j > nb) then
               call nek_stop_error('Slot out of range.', this_module, this_procedure)
            end if
            do i = 1, nfld
               call copy(f(1, i), fh(1, i, j), nbf)
            end do
         end subroutine tfft_get_raw

      !====================================================================
      !     DIAGNOSTICS AND OUTPUT
      !====================================================================

         subroutine tfft_spectrum(nrm)
      !! L2 norm over the cross-section of |fhat_m|, per harmonic, summed over
      !! all field columns. This is the table to look at before believing any
      !! picture: if it has not decayed by several decades at m = M, the
      !! transform is truncated, not converged.
            real(dp), dimension(0:), intent(out) :: nrm
      ! internal
            character(len=*), parameter :: this_procedure = 'tfft_spectrum'
            character(len=256) :: msg
            real(dp), external :: glsc3
            real(dp) :: acc
            integer :: m, i, j
            call check_closed(this_procedure)
            nrm = 0.0_dp
            do m = 0, min(mharm, ubound(nrm, 1))
               acc = 0.0_dp
               do i = 1, nfld
                  if (m == 0) then
                     acc = acc + glsc3(fh(1, i, 1), bm1, fh(1, i, 1), nbf)
                  else
                     do j = 2*m, 2*m + 1
                        acc = acc + glsc3(fh(1, i, j), bm1, fh(1, i, j), nbf)
                     end do
                  end if
               end do
               nrm(m) = sqrt(max(acc, 0.0_dp))
            end do
            call nek_log_message('temporal spectrum |f_m|:', this_module, this_procedure)
            do m = 0, min(mharm, ubound(nrm, 1))
               write (msg, '(5X,A,I2,A,E16.8)') 'm= ', m, ' : ', nrm(m)
               call nek_log_message(msg, this_module, this_procedure)
            end do
         end subroutine tfft_spectrum

         subroutine tfft_check_ubar(ifld)
      !! Consistency check against the flow-rate accumulator. The bulk average
      !! of the field harmonics IS the flow-rate harmonic, because it is the
      !! same quadrature applied to a linear functional of the same integrand.
      !! Agreement to roundoff means the packing, the normalisation and the
      !! time grid are all right; a factor of two here is the failure mode this
      !! kind of code actually has.
      !!
      !! Only meaningful for the baseflow transform, with ifld the column
      !! holding u_phi (3 in the drivers).
            integer, intent(in) :: ifld
      ! internal
            character(len=*), parameter :: this_procedure = 'tfft_check_ubar'
            character(len=256) :: msg
            real(dp), dimension(lfc) :: qf
            real(dp) :: qh, err, scl
            integer :: j
            call check_closed(this_procedure)
            call t2Dh%get_mflow(qf)
            scl = max(abs(qf(1)), atol_dp)
            call nek_log_message('field harmonics vs flow-rate harmonics:', this_module, this_procedure)
            do j = 1, min(nb, t2Dh%get_nf())
               qh = t2Dh%ubar_arr(fh(1, ifld, j))
               err = abs(qh - qf(j))/scl
               write (msg, '(5X,A,I2,A,E16.8,A,E16.8,A,E10.3)') 'j= ', j, ' : ubar= ', qh,
     &            ' qfour= ', qf(j), ' rel.diff= ', err
               call nek_log_message(msg, this_module, this_procedure)
               if (err > 1.0e-06_dp) then
                  call nek_log_warning('Field and flow-rate harmonics disagree well above '//
     &               'roundoff: check the packing, the normalisation and the time grid.',
     &               this_module, this_procedure)
               end if
            end do
         end subroutine tfft_check_ubar

         subroutine tfft_outpost(letters, mlist)
      !! One field file per (harmonic, spatial block, real/imaginary part).
      !!
      !! The 3-character name is <letter><mm>, with mm the harmonic index and
      !! the letter chosen by the caller:
      !!
      !!   nfld = 3   letters = 'ri'    r = Re fhat_m,      i = Im fhat_m
      !!   nfld = 6   letters = 'abcd'  a = Re_t of Re_x,   b = Im_t of Re_x
      !!                                c = Re_t of Im_x,   d = Im_t of Im_x
      !!
      !! Field slots inside each file are (u_z, u_R, --, --, u_phi): the third
      !! velocity component and the pressure are written as zeros because the
      !! 2Dh mesh has no third direction and the buffer stores no pressure.
      !!
      !! The harmonic index is also written into the file's time stamp, so a
      !! post-processor that reads the series sees m on the time axis.
            character(len=*), intent(in) :: letters
            integer, dimension(:), optional, intent(in) :: mlist
      ! internal
            character(len=*), parameter :: this_procedure = 'tfft_outpost'
            character(len=256) :: msg
            character(len=3) :: nam3
            real(dp), dimension(:, :), allocatable :: fre, fim
            real(dp), dimension(:), allocatable :: vzero, pzero
            real(dp) :: time_save
            logical :: ifto_save, ifpo_save
            integer :: m, mm, nblk, ib, i0, nout
            call check_closed(this_procedure)
            nblk = nfld/3
            if (nblk*3 /= nfld) then
               call nek_stop_error('tfft_outpost expects nfld to be a multiple of 3.',
     &            this_module, this_procedure)
            end if
            if (len_trim(letters) < 2*nblk) then
               write (msg, '(A,I0,A)') 'letters must supply ', 2*nblk, ' characters.'
               call nek_stop_error(msg, this_module, this_procedure)
            end if
            allocate (fre(lv, nfld), fim(lv, nfld), vzero(lv), pzero(lp2))
            call rzero(vzero, lv)
            call rzero(pzero, lp2)
            time_save = time
            ifto_save = ifto
            ifpo_save = ifpo
            ifto = .true.
            ifpo = .false.
            if (present(mlist)) then
               nout = size(mlist)
            else
               nout = mharm + 1
            end if
            do mm = 1, nout
               m = merge_index(mm, mlist)
               if (m < 0 .or. m > mharm) cycle
               call tfft_get(m, fre, fim)
      ! The harmonic rides on the time stamp rather than on the name, which
      ! only has room for one integer.
               time = real(m, dp)
               do ib = 1, nblk
                  i0 = 3*(ib - 1)
                  write (nam3, '(A1,I2.2)') letters(2*ib - 1:2*ib - 1), m
                  call set_fldindex(nam3, m + 1)
                  call outpost(fre(1, i0 + 1), fre(1, i0 + 2), vzero, pzero, fre(1, i0 + 3), nam3)
                  write (nam3, '(A1,I2.2)') letters(2*ib:2*ib), m
                  call set_fldindex(nam3, m + 1)
                  call outpost(fim(1, i0 + 1), fim(1, i0 + 2), vzero, pzero, fim(1, i0 + 3), nam3)
               end do
            end do
            time = time_save
            ifto = ifto_save
            ifpo = ifpo_save
            deallocate (fre, fim, vzero, pzero)
            write (msg, '(A,I0,A,A,A)') 'wrote ', nout*2*nblk, ' harmonic files for [',
     &         trim(label), ']'
            call nek_log_message(msg, this_module, this_procedure)
         end subroutine tfft_outpost

      !====================================================================
      !     ACCESSORS
      !====================================================================

         pure integer function tfft_mharm() result(m)
            m = mharm
         end function tfft_mharm

         pure integer function tfft_nfld() result(n)
            n = nfld
         end function tfft_nfld

         pure integer function tfft_nbf() result(n)
            n = nbf
         end function tfft_nbf

         pure logical function tfft_is_closed() result(l)
            l = closed
         end function tfft_is_closed

         pure real(dp) function tfft_omega() result(w)
            w = omega
         end function tfft_omega

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

         pure integer function merge_index(mm, mlist) result(m)
            integer, intent(in) :: mm
            integer, dimension(:), optional, intent(in) :: mlist
            if (present(mlist)) then
               m = mlist(mm)
            else
               m = mm - 1
            end if
         end function merge_index

      end module t2Dh_tfft