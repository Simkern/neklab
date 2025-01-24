      submodule (neklab_helix) helix_getters_setters
         implicit none

      contains

         !! MODE SWITCHES

         module procedure set_save_base
            self%save_2d_base = log_flag_change(if_save, self%save_2d_base, 'save_2d_base')
         end procedure set_save_base

         module procedure set_save_fft
            self%if_fft = log_flag_change(if_save_fft, self%if_fft, 'if_fft')
         end procedure set_save_fft
         
         module procedure set_newton
            self%if_newton = log_flag_change(if_newton, self%if_newton, 'if_newton')
         end procedure

         module procedure set_floquet
            self%if_floquet = log_flag_change(if_floquet, self%if_floquet, 'if_floquet')
         end procedure

         module procedure set_symmetry
            self%if_sym = log_flag_change(if_sym, self%if_sym, 'if_sym')
         end procedure

         !! LOGICAL FUNCTIONS

         module procedure is_steady
            if_steady = self%if_steady
         end procedure is_steady

         module procedure is_newton
            if_newton = self%if_newton
         end procedure is_newton

         module procedure is_floquet
            if_floquet = self%if_floquet
         end procedure is_floquet

         module procedure is_save_2d
            if_save_2d_base = self%save_2d_base
         end procedure is_save_2d

         module procedure is_save_fft
            if_save_fft = self%if_fft
         end procedure is_save_fft

         module procedure is_extracted_fft
            is_extracted = self%fft_is_extracted
         end procedure is_extracted_fft

         module procedure is_sym
            mesh_is_sym = self%if_sym
         end procedure is_sym

         module procedure is_torus
            mesh_is_torus = self%if_torus
         end procedure is_torus

         module procedure is_helix
            mesh_is_helix = self%if_helix
         end procedure is_helix

         module procedure is_lowner
            is_owner = .false.
            if (ie <= nelv) is_owner = self%lowner(ie)
         end procedure is_lowner

         module procedure is_gowner
            is_owner = .false.
            if (ie <= nelv) is_owner = self%gowner(ie)
         end procedure is_gowner

         !! GETTERS

         module procedure get_dpds
            integer :: i, j
            dpds = self%dpds
            if (present(phase)) then
               allocate(phase((lf+1)/2))
               phase = 0.0_dp
               j = 1
               do i = 2, lf, 2
                  j = j + 1
                  phase(j) = atan2(dpds(i+1),dpds(i))
               end do
            end if
         end procedure get_dpds

         module procedure get_fshape
            call copy(fshape, self%fshape, lv)
         end procedure get_fshape

         module procedure get_angle_s
            call copy(angle_s, self%as, lv)
         end procedure get_angle_s

         module procedure get_alpha
            call copy(alpha, self%alpha, lv)
         end procedure get_alpha

         module procedure get_length
            length = self%length
         end procedure get_length

         module procedure get_delta
            delta = self%delta
         end procedure get_delta

         module procedure get_diameter
            diameter = self%diameter
         end procedure get_diameter

         module procedure get_pitch_s
            pitch_s = self%pitch_s
         end procedure get_pitch_s

         module procedure get_radius
            radius = self%radius
         end procedure get_radius

         module procedure get_curv_radius
            curv_radius = self%curv_radius
         end procedure get_curv_radius

         module procedure get_phi
            phi = self%phi
         end procedure get_phi

         module procedure get_sweep
            sweep = self%sweep
         end procedure get_sweep

         module procedure get_period
            T = self%pulse_T
         end procedure get_period

         module procedure get_omega
            omega = self%omega
         end procedure get_omega

         module procedure get_nf
            n = self%nf
         end procedure get_nf

         module procedure get_Wo
            Wo = self%womersley
         end procedure get_Wo

         module procedure get_nsteps
            ns = 0
            if (self%nsteps /= 0) then
               ns = self%nsteps
            else
               call nek_stop_error('nsteps not computed.', this_module, 'get_nsteps')
            end if
         end procedure get_nsteps

         module procedure get_dt_minmax
            if (self%min_dt == 100.0_dp) then
               call nek_log_message('min_dt not computed.', this_module, 'get_dt_minmax')
            end if
            if (self%max_dt == 0.0_dp) then
               call nek_log_message('max_dt not computed.', this_module, 'get_dt_minmax')
            end if
            if (self%min_dt /= 100.0_dp .and. self%max_dt /= 0.0_dp) then
               dt_minmax(1) = self%min_dt
               dt_minmax(2) = self%max_dt
            end if
         end procedure get_dt_minmax

         module procedure get_ubar_lag
            ubar_lag = 0.0_dp
            if (self%ubar_lag /= 0.0_dp) then
               ubar_lag = self%ubar_lag
            else
               call nek_stop_error('ubar_lag not computed.', this_module, 'get_ubar_lag')
            end if
         end procedure get_ubar_lag

         module procedure get_lsegment
            local_segment = 0
            if (ie <= nelv) local_segment = self%lsegment(ie)
         end procedure get_lsegment
      
         module procedure get_gsegment
            global_segment = 0
            if (ie <= nelv) global_segment = self%gsegment(ie)
         end procedure get_gsegment

         module procedure get_v2d
            v2d = 0.0_dp
            if (ix <= lx1) then
               if (iy <= ly1) then
                  if (iseg <= self%n2d_lown) then
                     if (ifld <= self%nload) then
                        if (icomp == 1) then
                           v2d = self%vx2d(ix,iy,iseg,ifld)
                        else if (icomp == 2) then
                           v2d = self%vy2d(ix,iy,iseg,ifld)
                        else if (icomp == 2) then
                           v2d = self%vz2d(ix,iy,iseg,ifld)
                        end if
                     end if
                  end if
               end if
            end if
         end procedure get_v2d

         !! SETTERS

         module procedure set_dpds
            logical :: reset_dpds
            reset_dpds = optval(reset, .true.)
            if (reset_dpds) self%dpds = 0.0_dp
            self%dpds = self%dpds + dpds
         end procedure set_dpds

         module procedure set_nsteps
            if (ns /= 0) then
               self%nsteps = ns
            else
               call nek_log_message('input is zero. nsteps not set.', this_module, 'set_nsteps')
            end if
         end procedure set_nsteps

         !! HELPER ROUTINE

         logical function log_flag_change(new_flag, old_flag, flag_name) result(out_flag)
            logical, intent(in) :: new_flag
            logical, intent(in) :: old_flag
            character(len=*), intent(in) :: flag_name
            if (new_flag .neqv. old_flag) then
               if (new_flag) then
                  call nek_log_message(trim(flag_name)//' switched ON.', this_module, 'set_logical')
               else
                  call nek_log_message(trim(flag_name)//' switched OFF.', this_module, 'set_logical')
               end if
            else
               if (new_flag) then
                  call nek_log_information(trim(flag_name)//' switched ON. (unchanged)', this_module, 'set_logical')
               else
                  call nek_log_information(trim(flag_name)//' switched OFF. (unchanged)', this_module, 'set_logical')
               end if
            end if
            out_flag = new_flag
         end function log_flag_change
      
      end submodule helix_getters_setters