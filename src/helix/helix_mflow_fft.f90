      submodule (neklab_helix) helix_mflow_fft
         implicit none
      
      contains

         module procedure reset_mflow_fft
            self%ubar_lag = self%compute_ubar(vx, vy, vz) ! compute ubar at t = 0
            ! zero out data arrays
            self%fftv = 0.0_dp
            self%fft_time = 0.0_dp               ! reset integration time
            self%fft_is_extracted = .false.
            call self%set_save_fft(.true.)
         end procedure reset_mflow_fft

         module procedure compute_mflow_fft
            character(len=*), parameter :: this_procedure = 'compute_mflow_fft'
            integer :: i, j
            real(dp) :: ubar, tau, dtau, twopi, pd
            logical :: var_dt_
            real(dp) :: ubar_old, tau_old, dt0, dfftv1, dfftv2
            var_dt_ = optval(var_dt, .false.)
            pd = optval(period, self%pulse_T)
            if (pd == 0.0_dp) call nek_stop_error('Period not set or zero.', this_module, this_procedure)
            if (self%is_save_fft()) then
               call lk_timer%start('neklab_helix_'//this_procedure)
               twopi = 8.0_dp*atan(1.0_dp)
               ! compute period, current ubar and time constants
               ubar = self%compute_ubar(vx,vy,vz)
               tau  = time/pd
               dtau = dt/pd
               if (var_dt_) then ! variable timestep integration
                  ! get ubar and time of previous timestep
                  ubar_old = self%ubar_lag
                  tau_old = (time - dt)/pd
                  if (ubar*ubar_old < 0.0) then ! zero crossing
                     ! find zero crossing
                     dt0 = -(ubar - ubar_old)/ubar_old
                     ! fill up fft array
                     self%fftv(1) = self%fftv(1) + (ubar_old*dt0 + ubar*(1 - dt0))*0.5*dtau
                     j = 1
                     do i = 2, 2*nfft, 2
                        dfftv1 = ubar_old*cos(j*twopi*tau_old)*   dt0
                        dfftv2 = ubar    *cos(j*twopi*tau    )*(1-dt0)
                        self%fftv(i)   = self%fftv(i)   + (dfftv1 + dfftv2)*dtau
                        dfftv1 = ubar_old*sin(j*twopi*tau_old)*   dt0
                        dfftv2 = ubar    *sin(j*twopi*tau    )*(1-dt0)
                        self%fftv(i+1) = self%fftv(i+1) + (dfftv1 + dfftv2)*dtau
                        j = j + 1
                     end do
                  else
                     ! fill up fft array
                     self%fftv(1) = self%fftv(1) + (ubar_old + ubar)*0.5*dtau
                     j = 1
                     do i = 2, 2*nfft, 2
                        dfftv1 = ubar_old*cos(j*twopi*tau_old)
                        dfftv2 = ubar    *cos(j*twopi*tau    )
                        self%fftv(i)   = self%fftv(i)   + (dfftv1 + dfftv2)*dtau
                        dfftv1 = ubar_old*sin(j*twopi*tau_old)
                        dfftv2 = ubar    *sin(j*twopi*tau    )
                        self%fftv(i+1) = self%fftv(i+1) + (dfftv1 + dfftv2)*dtau
                        j = j + 1
                     end do
                  end if
                  ! update lagged ubar
                  self%ubar_lag = ubar
               else ! constant timestep
                  ! fill up fft array
                  self%fftv(1) = self%fftv(1) + ubar*dtau
                  j = 1
                  do i = 2, 2*nfft, 2
                     self%fftv(i)   = self%fftv(i)   + ubar*cos(j*twopi*tau)*dtau
                     self%fftv(i+1) = self%fftv(i+1) + ubar*sin(j*twopi*tau)*dtau
                     j = j + 1
                  end do
               end if
               ! increment integration time
               self%fft_time = self%fft_time + dt
               if (io_rank()) print '(A,2(F18.12),F12.6)', 'neklab_helix: Compute mflow fft ', self%fft_time, pd, self%fft_time/pd
               call lk_timer%stop('neklab_helix_'//this_procedure)
            end if
         end procedure compute_mflow_fft

         module procedure extract_mflow_fft
            character(len=*), parameter :: this_procedure = 'extract_mflow_fft'
            real(dp) :: pd, pd_chk
            integer :: i, j, nprint
            logical :: if_amplitude_
            character(len=1024) :: msg
            character(len=128), parameter :: fmt = '(A,1X,F16.8,1X,A,*(1X,F16.8))'
            pd = optval(period, self%pulse_T)
            if_amplitude_ = optval(if_amplitude, .true.)
            if (pd == 0.0_dp) call nek_stop_error('Period not set or zero.', this_module, this_procedure)
            ! extract the computed FFT data, compute amplitudes and phases
            self%fft_rtime = self%fft_time ! total integration time since last call
            call copy(self%mflow, self%fftv, 2*nfft+1)
            self%mflow_amplitude(1) = self%mflow(1)
            self%mflow_phase(1) = 0.0_dp
            j = 1
            do i = 2, 2*nfft, 2
               j = j + 1
               self%mflow_amplitude(j) = sqrt(self%mflow(i)**2 + self%mflow(i+1)**2)
               self%mflow_phase(j)     = atan2(self%mflow(i+1),self%mflow(i))
            end do
            self%fftv = 0.0_dp
            self%fft_time = 0.0_dp               ! reset integration time
            ! sanity period check
            pd_chk = self%fft_rtime/pd
            ! print result
            if (if_amplitude_) then
               nprint = (self%nf+1)/2
               write(msg,fmt) 'Period',self%fft_rtime,'massflow FT amplitude  ',self%mflow_amplitude(:nprint)
               call nek_log_message(msg, this_module, this_procedure)
               write(msg,fmt) 'Period',self%fft_rtime,'massflow FT phase angle',self%mflow_phase(:nprint)
               call nek_log_message(msg, this_module, this_procedure)
               if (self%omega /= 0.0_dp) then
                  write(msg,fmt) 'Period',self%fft_rtime,'massflow FT t-shift    ',self%mflow_phase(:nprint)/self%omega
                  call nek_log_debug(msg, this_module, this_procedure)
               end if
            else
               write(msg,fmt) 'Period',self%fft_rtime,'massflow FT cmplx',self%mflow(:self%nf)
               call nek_log_message(msg, this_module, this_procedure)
            end if
            if (abs(pd_chk - 1.0_dp) > 1.0e-06_dp) then
               write(msg, '(A,E15.8,A,2(F16.8,1X))') 'Period check: ', pd_chk - 1.0_dp, ': ', pd, self%fft_rtime
               call nek_log_message(msg, this_module, this_procedure)
               call nek_stop_error('Period check failed. Maybe the integration time does not equal the period precisely.')
            end if
            self%fft_is_extracted = .true.
         end procedure extract_mflow_fft

         module procedure get_mflow_fft
            character(len=*), parameter :: this_procedure = 'get_mflow_fft'
            logical :: if_amplitude_
            character(len=128) :: msg
            if_amplitude_ = optval(if_amplitude, .true.)
            if (self%is_extracted_fft()) then
               if (if_amplitude_) then
                  allocate(mflow(nfft+1))
                  mflow = self%mflow_amplitude
                  if (present(phase)) then
                     allocate(phase(nfft+1))
                     phase = self%mflow_phase
                  end if
               else
                  allocate(mflow(2*nfft+1))
                  mflow = self%mflow
                  if (present(phase)) then
                     msg = 'To obtain phase information, use if_amplitude = .true. No phase information returned.'
                     call nek_log_message(msg, this_module, this_procedure)
                  end if
               end if
            else
               call nek_stop_error('mflow FT data has not been extracted.', this_module, this_procedure)
            end if
         end procedure get_mflow_fft
      
      end submodule helix_mflow_fft