      submodule (neklab_helix) helix_utils
         implicit none
      
      contains

         module procedure init_geom
            character(len=*), parameter :: this_procedure = 'init_geom'
            real(dp) :: minv, maxv, torus_r, s_angle
            real(dp) :: x_torus, y_torus, z_torus, sweep, r
            real(dp), dimension(lx1,ly1,lz1,lelv) :: tmp, pipe_r
            integer :: ix, iy, iz, ie, itmp, nxy
            character(len=128) :: msg
            ! functions
            real(dp), external :: glmax, glmin

            if (self%is_initialized) call stop_error('Attempting to reinitialize the mesh', this_module, this_procedure)
            call lk_timer%start('neklab_helix_'//this_procedure)

            !  Geometry modification for helical pipe

            pi = 4.0_dp*atan(1.0_dp)
            nxy = lx1*ly1

            if (self%is_sym()) then
               ! half mesh
               call rescale_x(xm1,0.0_dp,self%radius)       ! x in [  0, r ]
            else
               call rescale_x(xm1,-self%radius,self%radius) ! x in [ -r, r ]
            end if
            call rescale_x(ym1,-self%radius,self%radius)    ! y in [ -r, r ]
            call rescale_x(zm1,0.0_dp,1.0_dp)               ! z in [  0, 1 ]
            
            !  rotate mesh to set the center of the pipe along x-axis
            call copy(tmp,  xm1, lv)
            call copy(xm1,  zm1, lv)   ! x <--  z
            call copy(zm1, -tmp, lv)   ! z <-- -x
            call copy(self%zax,zm1,lv) ! zax set before curvature in z is added!

            call copy(pipe_r, ym1,lv) ! local distance from pipe center
            minv = glmin(xm1,lv); maxv = glmax(xm1,lv)
            write(msg,'(2(A,F16.12),A)') 'x: min ', minv, ' max ', maxv, ' (streamwise)'
            call nek_log_message(msg, this_module, this_procedure)
            minv = glmin(ym1,lv); maxv = glmax(ym1,lv)
            write(msg,'(2(A,F16.12))') 'y: min ', minv, ' max ', maxv
            call nek_log_message(msg, this_module, this_procedure)
            minv = glmin(zm1,lv); maxv = glmax(zm1,lv)
            write(msg,'(2(A,F16.12))') 'z: min ', minv, ' max ', maxv
            call nek_log_message(msg, this_module, this_procedure)
            call nek_log_message('Mesh rescaled and rotated.', this_module, this_procedure)

            ! Set up and extract 2D geometry
            call self%init_2d_geom(if_debug)
            call nek_log_message('2D geometry extracted.', this_module, this_procedure)

            ! Morph the mesh into a torus
            if (self%is_torus()) then
               ! rescale the x axis
               minv = glmin(xm1,lv)
               maxv = glmax(xm1,lv)
               xm1 = self%sweep/(maxv-minv) * xm1 ! x in [ 0, max_sweep_angle ]
               call copy(self%sweep_angle,xm1,lv) ! save sweep angle
               torus_r = self%curv_radius
               do ie = 1, nelv
               do iz = 1, lz1
               do iy = 1, ly1
               do ix = 1, lx1
                  r     = pipe_r(ix,iy,iz,ie)
                  sweep = self%sweep_angle(ix,iy,iz,ie)
                  self%ox(ix,iy,iz,ie) = torus_r * sin(sweep)
                  self%oy(ix,iy,iz,ie) = torus_r * cos(sweep)
                  xm1(ix,iy,iz,ie)     = r * sin(sweep) + self%ox(ix,iy,iz,ie)
                  ym1(ix,iy,iz,ie)     = r * cos(sweep) + self%oy(ix,iy,iz,ie)
               end do
               end do
               end do
               end do
               write(msg,'(2(A,F16.12))') 'radius: ', torus_r
               call nek_log_message(msg, this_module, this_procedure)
               minv = glmin(self%sweep_angle,lv); maxv = glmax(self%sweep_angle,lv)
               write(msg,'(2(A,F16.12))') 'sweep: min ', minv, ' max ', maxv
               call nek_log_message(msg, this_module, this_procedure)
               call nek_log_message('Mesh morphed into torus.', this_module, this_procedure)
            end if
            call copy(self%xax, xm1, lv) ! xax set before curvature in z is added!
            call copy(self%yax, ym1, lv) ! yax set before curvature in z is added!

            ! Morph the torus into a helix
            if (self%is_helix()) then
               do ie = 1, nelv
               do iz = 1, lz1
               do iy = 1, ly1
               do ix = 1, lx1
                  x_torus = xm1(ix,iy,iz,ie)
                  y_torus = ym1(ix,iy,iz,ie)
                  z_torus = self%zax(ix,iy,iz,ie)
                  sweep   = self%sweep_angle(ix,iy,iz,ie)
                  xm1(ix,iy,iz,ie) = x_torus - z_torus*sin(self%phi)*cos(sweep)
                  ym1(ix,iy,iz,ie) = y_torus + z_torus*sin(self%phi)*sin(sweep)
                  zm1(ix,iy,iz,ie) = sweep*self%pitch_s + z_torus*cos(self%phi)
               enddo
               enddo
               enddo
               enddo
               minv = glmin(zm1,lv); maxv = glmax(zm1,lv)
               write(msg,'(2(A,F16.12))') 'z: min ', minv, ' max ', maxv
               call nek_log_message(msg, this_module, this_procedure)
               call nek_log_message('Mesh morphed into helix.', this_module, this_procedure)
            end if
            param(59) = 1.   !  All elements deformed

            ! Streamwise angle in the equatorial plane
            if (self%is_helix() .or. self%is_torus()) then
               self%as = atan2(self%xax, self%yax) ! clockwise from y axis
            else
               self%as = 0.0_dp ! straight pipe
            end if
            ! Angle within cross-sectional plane
            self%alpha = atan2(self%zax, pipe_r)

            !do ie = 1, nelv
            !do iz = 1, lz1
            !do iy = 1, ly1
            !do ix = 1, lx1
            !   x_torus = xm1(ix,iy,iz,ie)
            !   y_torus = ym1(ix,iy,iz,ie)
            !   z_torus = zm1(ix,iy,iz,ie)
            !   self%r(ix,iy,iz,ie) = sqrt((sqrt(x_torus**2+y_torus**2)-self%curv_radius)**2 + z_torus**2)
            !enddo
            !enddo
            !enddo
            !enddo

            itmp = istep
            istep = 0
            call comment() ! set internal variable ifcour for standard timestep logging (--> needs to be called at istep == 0)
            istep = itmp
            self%is_initialized = .true.
            call lk_timer%stop('neklab_helix_'//this_procedure)
            
         end procedure init_geom

         module procedure init_flow
            character(len=*), parameter :: this_procedure = 'init_flow'
            integer :: i, n
            logical :: reset_nf_
            character(len=128) :: msg
            pi = 4.0_dp*atan(1.0_dp)
            self%womersley = optval(womersley, 0.0_dp)
            reset_nf_ = optval(reset_nf, .false.)
            n = size(dpds)
            write(msg,'(A,I0,A)') 'nf = ', n, ' forcing components provided.'
            call nek_log_information(msg, this_module, this_procedure)
            if (self%nf == 0) then
               self%nf = n
               write(msg,'(A,I0)') 'Number of considered forcing components set to nf = ', self%nf
               call nek_log_message(msg, this_module, this_procedure)
            else if (reset_nf_) then
               if (n /= self%nf) then
                  self%nf = n
                  write(msg,'(A,I0)') 'Number of considered forcing components reset to nf = ', self%nf
                  call nek_log_warning(msg, this_module, this_procedure)
               else
                  write(msg,'(A,I0)') 'Number of considered forcing components unchanged. nf = ', self%nf
                  call nek_log_information(msg, this_module, this_procedure)
               end if
            end if
            if (self%womersley /= 0.0_dp) then
               self%if_steady = .false.
               call nek_log_information('Running unsteady case.', this_module, this_procedure)
               self%omega     = (self%womersley**2)*cpfld(1,1)    ! pulsation frequency
               self%pulse_T   = 2.0_dp*pi/self%omega               ! pulsation period
               if (self%nf == 1) then
                  write(msg,'(A,I0,A)') 'nf > ', 1, ' forcing components required for unsteady case.'
                  call nek_stop_error(msg, this_module, this_procedure)
               else if (mod(self%nf,2)==0) then
                  msg = 'Unsteady case requires an uneven number of forcing components'
                  call nek_stop_error(msg, this_module, this_procedure)
               end if
               call nek_log_message('Unsteady flow parameters set.', this_module, this_procedure)
            else
               self%if_steady = .true.
               call nek_log_message('Running steady case.', this_module, this_procedure)
               self%omega     = 0.0_dp   
               self%pulse_T   = 0.0_dp  
               if (self%nf > 1) then
                  write(msg,'(A,I0,A)') 'Components nf > ', 1, ' will be ignored.'
                  call nek_log_warning(msg, this_module, this_procedure)
               end if
               call nek_log_message('Steady flow parameters set.', this_module, this_procedure)
            end if
            self%dpds = 0.0_dp
            self%dpds(:self%nf) = dpds
            call pipe%compute_bf_forcing(0.0_dp) ! ensure that the forcing is set (in particular for steady flows)
            call nek_log_message('Baseflow forcing set.', this_module, this_procedure)
            call self%parameter_summary()
         end procedure init_flow

         module procedure compute_fshape
            integer :: ix, iy, iz, ie
            real(dp) :: helix_r2, r, rr, alpha
            self%fshape = 0.0_dp
            if (self%is_helix() .or. self%is_torus()) then
               do ie = 1, lelv
               do iz = 1, lz1
               do iy = 1, ly1
               do ix = 1, lx1
                  ! Distance from the origin in the equatorial plane
                  helix_r2 = self%xax(ix,iy,iz,ie)**2 + self%yax(ix,iy,iz,ie)**2
                  ! Distance from the pipe center in the equatorial plane
                  r = sqrt(helix_r2) - self%curv_radius
                  ! Azimuthal angle in the cross-sectional plane
                  alpha = atan2(r, self%zax(ix,iy,iz,ie))
                  ! Radial position in the cross-sectional plane
                  rr = sqrt(r**2 + self%zax(ix,iy,iz,ie)**2)
                  ! Compute fshape
                  self%fshape(ix,iy,iz,ie) = 1.0_dp / abs(1.0_dp + self%delta * rr * sin(alpha))
               end do
               end do
               end do
               end do
            else
               self%fshape = 1.0_dp
            end if
         end procedure compute_fshape

         module procedure forcing_amplitude
            integer :: i
            complex(dp) :: eiwt, dpds
            f = self%dpds(1)
            if (.not.self%if_steady) then
               eiwt = cexp(imag * self%omega * t)
               do i = 2, self%nf, 2
                  dpds = self%dpds(i) + imag*self%dpds(i+1)
                  f = f + 2.0_dp * real(dpds * eiwt)
               end do
            end if
         end procedure forcing_amplitude

         module procedure compute_bf_forcing
            integer :: ix, iy, iz, ie
            real(dp) :: fs, phi
            real(dp), dimension(lx1,ly1,lz1,lelv) :: ffx, ffy, ffz

            if (self%is_torus()) then
               fs = self%forcing_amplitude(t) / self%curv_radius
            else
               fs = self%forcing_amplitude(t)
            end if

            if (self%is_helix()) then
               phi = self%phi
               do ie = 1, nelv
               do iz = 1, lz1
               do iy = 1, ly1
               do ix = 1, lx1
                  ffx(ix,iy,iz,ie) =  fs * self%fshape(ix,iy,iz,ie) * cos(phi) * cos(self%as(ix,iy,iz,ie))
                  ffy(ix,iy,iz,ie) = -fs * self%fshape(ix,iy,iz,ie) * cos(phi) * sin(self%as(ix,iy,iz,ie))
                  ffz(ix,iy,iz,ie) =  fs * self%fshape(ix,iy,iz,ie) * sin(phi)
               end do
               end do
               end do
               end do
            else if (self%is_torus()) then
               do ie = 1, nelv
               do iz = 1, lz1
               do iy = 1, ly1
               do ix = 1, lx1
                  ffx(ix,iy,iz,ie) =  fs * self%fshape(ix,iy,iz,ie) * cos(self%as(ix,iy,iz,ie))
                  ffy(ix,iy,iz,ie) = -fs * self%fshape(ix,iy,iz,ie) * sin(self%as(ix,iy,iz,ie))
               end do
               end do
               end do
               end do
               ffz = 0.0_dp
            else
               ffx = fs
               ffy = 0.0_dp
               ffz = 0.0_dp
            end if

            ! set baseflow forcing
            call set_neklab_forcing(ffx, ffy, ffz, ipert=0)

         end procedure compute_bf_forcing

         module procedure compute_usrt
            integer :: ix, iy, iz, ie
            real(dp) :: phi, a, s, ux, uy, uz, utmp, vtmp
            if (self%is_helix()) then
               phi = self%phi
               do ie = 1, nelv
               do iz = 1, lz1
               do iy = 1, ly1
               do ix = 1, lx1
                  s  = self%as(ix,iy,iz,ie)
                  a  = self%alpha(ix,iy,iz,ie)
                  ux = u(ix,iy,iz,ie)
                  uy = v(ix,iy,iz,ie)
                  uz = w(ix,iy,iz,ie)
                  utmp            = sin(s)*ux + cos(s)*uy
                  vtmp            = sin(phi) * (-cos(s)*ux - sin(s)*uy) + cos(phi)*uz
                  us(ix,iy,iz,ie) = cos(phi) * ( cos(s)*ux - sin(s)*uy) + sin(phi)*uz
                  ur(ix,iy,iz,ie) = cos(a) * utmp + sin(a) * vtmp
                  ut(ix,iy,iz,ie) = sin(a) * utmp - cos(a) * vtmp
               end do
               end do
               end do
               end do
            else if (self%is_torus()) then
               do ie = 1, nelv
               do iz = 1, lz1
               do iy = 1, ly1
               do ix = 1, lx1
                  s  = self%as(ix,iy,iz,ie)
                  a  = self%alpha(ix,iy,iz,ie)
                  ux = u(ix,iy,iz,ie)
                  uy = v(ix,iy,iz,ie)
                  utmp = sin(s)*ux + cos(s)*uy
                  vtmp = w(ix,iy,iz,ie)
                  us(ix,iy,iz,ie) = cos(s)*ux - sin(s)*uy
                  ur(ix,iy,iz,ie) = cos(a) * utmp + sin(a) * vtmp
                  ut(ix,iy,iz,ie) = sin(a) * utmp - cos(a) * vtmp
               end do
               end do
               end do
               end do
            else
               do ie = 1, nelv
               do iz = 1, lz1
               do iy = 1, ly1
               do ix = 1, lx1
                  a = self%alpha(ix,iy,iz,ie)
                  utmp = v(ix,iy,iz,ie)
                  vtmp = w(ix,iy,iz,ie)
                  us(ix,iy,iz,ie) = u(ix,iy,iz,ie)
                  ur(ix,iy,iz,ie) = cos(a) * utmp + sin(a) * vtmp
                  ut(ix,iy,iz,ie) = sin(a) * utmp - cos(a) * vtmp
               end do
               end do
               end do
               end do
            end if
         end procedure compute_usrt

         module procedure compute_uxyz
            integer :: ix, iy, iz, ie
            real(dp) :: phi, a, s, us_loc, ur_loc, ut_loc, utmp, vtmp
            if (self%is_helix()) then
               phi = self%phi
               do ie = 1, nelv
               do iz = 1, lz1
               do iy = 1, ly1
               do ix = 1, lx1
                  s  = self%as(ix,iy,iz,ie)
                  a  = self%alpha(ix,iy,iz,ie)
                  us_loc = us(ix,iy,iz,ie)
                  ur_loc = ur(ix,iy,iz,ie)
                  ut_loc = ut(ix,iy,iz,ie)
                  ! Invert (ur, ut) => (utmp, vtmp)
                  utmp =  cos(a) * ur_loc + sin(a) * ut_loc
                  vtmp = -sin(a) * ur_loc + cos(a) * ut_loc
                  ! Now solve for ux, uy, uz from utmp, vtmp, us
                  u(ix,iy,iz,ie) =  cos(s) * (cos(phi) * us_loc - sin(phi) * vtmp) + sin(s) * utmp
                  v(ix,iy,iz,ie) = -sin(s) * (cos(phi) * us_loc - sin(phi) * vtmp) + cos(s) * utmp
                  w(ix,iy,iz,ie) = sin(phi) * us_loc + cos(phi) * vtmp
               end do
               end do
               end do
               end do
            else if (self%is_torus()) then
               do ie = 1, nelv
               do iz = 1, lz1
               do iy = 1, ly1
               do ix = 1, lx1
                  s  = self%as(ix,iy,iz,ie)
                  a  = self%alpha(ix,iy,iz,ie)
                  us_loc = us(ix,iy,iz,ie)
                  ur_loc = ur(ix,iy,iz,ie)
                  ut_loc = ut(ix,iy,iz,ie)
                  ! Invert (ur, ut) => (utmp, vtmp)
                  utmp =  cos(a) * ur_loc + sin(a) * ut_loc
                  vtmp = -sin(a) * ur_loc + cos(a) * ut_loc
                  ! Invert us and utmp to get ux, uy
                  u(ix,iy,iz,ie) =  cos(s) * us_loc + sin(s) * utmp
                  v(ix,iy,iz,ie) = -sin(s) * us_loc + cos(s) * utmp
                  w(ix,iy,iz,ie) = vtmp
               end do
               end do
               end do
               end do
            else
               do ie = 1, nelv
               do iz = 1, lz1
               do iy = 1, ly1
               do ix = 1, lx1
                  a = self%alpha(ix,iy,iz,ie)
                  utmp =  cos(a) * ur(ix,iy,iz,ie) + sin(a) * ut(ix,iy,iz,ie)
                  vtmp = -sin(a) * ur(ix,iy,iz,ie) + cos(a) * ut(ix,iy,iz,ie)
                  u(ix,iy,iz,ie) = us(ix,iy,iz,ie)
                  v(ix,iy,iz,ie) = utmp
                  w(ix,iy,iz,ie) = vtmp
               end do
               end do
               end do
               end do
            end if
         end procedure compute_uxyz
   
         module procedure compute_ubar
            integer :: ix, iy, iz, ie
            real(dp) :: num, den, us, us_r, ux, uy, uz, phi, s, a, fs
            real(dp), external :: glsum
            call lk_timer%start('neklab_helix_compute_ubar')
            num = 0.0_dp
            den = 0.0_dp
            phi = self%phi
            do ie = 1, nelv
            do iz = 1, lz1
            do iy = 1, ly1
            do ix = 1, lx1
               ux = u(ix,iy,iz,ie)
               uy = v(ix,iy,iz,ie)
               uz = w(ix,iy,iz,ie)
               s  = self%as(ix,iy,iz,ie)
               a  = self%alpha(ix,iy,iz,ie)
               fs = self%fshape(ix,iy,iz,ie)
               us = cos(phi)*(cos(s)*ux - sin(s)*uy) + sin(phi)*uz
               us_r = us * fs ! u/r
               num = num + us_r*bm1(ix,iy,iz,ie)
               den = den + fs  *bm1(ix,iy,iz,ie)
            end do
            end do
            end do
            end do
            num = glsum(num,1)
            den = glsum(den,1)
            ubar = num/den  ! "1/r"-weighted volumetric average of streamwise velocity
            call lk_timer%stop('neklab_helix_compute_ubar')
         end procedure compute_ubar
            
         module procedure shift_mflow_phase
            character(len=*), parameter :: this_procedure = 'shift_mflow_phase'
            integer :: i
            real(dp) :: dpdsr, dpdsi, alpha, dalpha, dt_phase, prop
            real(dp) :: dpds(lf)
            character(len=128) :: msg
            ! extract current forcing components
            if (.not.self%is_steady()) then
               call self%get_dpds(dpds)
               i = 2*(icomp-1)
               dpdsr = self%dpds(i)
               dpdsi = self%dpds(i+1)
               ! get current mflow phase angle
               dalpha = self%mflow_phase(icomp) - target_mflow_phase
               dt_phase = dalpha/pipe%get_omega()
               prop = dt_phase/pipe%get_period()*100
               write(msg,'(A,I0)') 'adjusting forcing component: ', icomp
               call nek_log_message(msg, this_module, this_procedure)
               write(msg,'(3X,A,F16.8)') 'dalpha  = ', dalpha
               call nek_log_message(msg, this_module, this_procedure)
               write(msg,'(3X,A,F16.8,A,F10.5,A)') 'dt_phase= ', dt_phase , '  (', prop, ' % T)'
               call nek_log_message(msg, this_module, this_procedure)
               ! update forcing (rotation) to remove shift
               self%dpds(i  ) = dpdsr*cos(dalpha) - dpdsi*sin(dalpha)
               self%dpds(i+1) = dpdsr*sin(dalpha) + dpdsi*cos(dalpha)
            end if
         end procedure shift_mflow_phase

         module procedure gfldr_torus
            ! internal
            character(len=132) :: fname
            fname = 'r2dtorus001.fld'
            call gfldr(rstfname)
            call self%set_save_base(.true.)
	         call self%save_2d_fields(vx, vy, vz)
            call self%write_2d(fname, only_mesh=.false.)
            call self%load_baseflow(vx, vy, vz, fname, 1)
         end procedure

         module procedure setup_summary
            character(len=128) :: msg
            if (self%is_initialized) then
               call nek_log_message('##  HELIX SETUP ##', module=this_module)
               call nek_log_message('Geometry:', module=this_module)
               write (msg, '(A,F15.8)') padl('length:', 20), pipe%length
               call nek_log_message(msg, module=this_module, fmt='(5X,A)')
               write (msg, '(A,F15.8)') padl('diameter:', 20), pipe%diameter
               call nek_log_message(msg, module=this_module, fmt='(5X,A)')
               write (msg, '(A,F15.8)') padl('radius:', 20), pipe%radius
               call nek_log_message(msg, module=this_module, fmt='(5X,A)')
               write (msg, '(A,F15.8)') padl('pitch_s:', 20), pipe%pitch_s
               call nek_log_message(msg, module=this_module, fmt='(5X,A)')
               write (msg, '(A,F15.8)') padl('delta:', 20), pipe%delta
               call nek_log_message(msg, module=this_module, fmt='(5X,A)')
               write (msg, '(A,F15.8)') padl('curv_radius:', 20), pipe%curv_radius
               call nek_log_message(msg, module=this_module, fmt='(5X,A)')
               call nek_log_message('Angles:', module=this_module)
               write (msg, '(A,F15.8)') padl('rise angle:', 20), pipe%phi
               call nek_log_message(msg, module=this_module, fmt='(5X,A)')
               write (msg, '(A,F15.8)') padl('sweep angle:', 20), pipe%sweep
               call nek_log_message(msg, module=this_module, fmt='(5X,A)')
               call nek_log_message('Mesh:', module=this_module)
               write (msg, '(A,I8)') padl('slices:', 20), pipe%nslices
               call nek_log_message(msg, module=this_module, fmt='(5X,A)')
               write (msg, '(A,I8)') padl('nel/slice:', 20), pipe%nelf
               call nek_log_message(msg, module=this_module, fmt='(5X,A)')
               write (msg, '(A,L8)') padl('symmetry:', 20), pipe%is_sym()
               call nek_log_message(msg, module=this_module, fmt='(5X,A)')
               write (msg, '(A,L8)') padl('toroidal mesh', 20), pipe%is_torus()
               call nek_log_message(msg, module=this_module, fmt='(5X,A)')
               write (msg, '(A,L8)') padl('helical mesh', 20), pipe%is_helix()
               call nek_log_message(msg, module=this_module, fmt='(5X,A)')
            else
               call nek_log_warning('helix instance not initialized', module=this_module, fmt='(A)')
            end if
         end procedure setup_summary
         
         module procedure parameter_summary
            character(len=128) :: msg
            if (self%is_initialized) then
               call nek_log_message('##  HELIX PARAMETERS ##', module=this_module)
               call nek_log_message('Flow:', module=this_module)
               write (msg, '(A,L8)') padl('steady:', 20), self%if_steady
               call nek_log_message(msg, module=this_module, fmt='(5X,A)')
               if (.not.self%is_steady()) then
                  write (msg, '(A,F15.8)') padl('Wo:', 20), self%womersley
                  call nek_log_message(msg, module=this_module, fmt='(5X,A)')
                  write (msg, '(A,F15.8)') padl('omega:', 20), self%omega
                  call nek_log_message(msg, module=this_module, fmt='(5X,A)')
                  write (msg, '(A,F15.8)') padl('T:', 20), self%pulse_T
                  call nek_log_message(msg, module=this_module, fmt='(5X,A)')
               end if
               call nek_log_message('Forcing:', module=this_module)
               call self%forcing_summary()
            else
               call nek_log_warning('helix instance not initialized', module=this_module, fmt='(A)')
            end if
         end procedure parameter_summary

         module procedure forcing_summary
            integer :: i
            real(dp) :: dpds_norm, dpds_angle_rad
            character(len=128) :: msg, fmt
            if (self%is_initialized) then
               write (msg, '(4(A,F16.12))') padl('dpds_00:', 20), self%dpds(1), ' ', 0.0_dp,
     &               ' | ', self%dpds(1), ' | ', 0.0_dp 
               call nek_log_message(msg, module=this_module, fmt='(5X,A)')
               do i = 2, self%nf, 2
                  write(fmt,'("dpds_",I2.2,":")') i/2
                  dpds_norm      = sqrt(self%dpds(i)**2 + self%dpds(i+1)**2)
                  dpds_angle_rad = atan2(self%dpds(i+1),self%dpds(i))
                  write (msg, '(4(A,F16.12))') padl(trim(fmt), 20), self%dpds(i), ' ', self%dpds(i+1), 
     &               ' | ', dpds_norm, ' | ', dpds_angle_rad
                  call nek_log_message(msg, module=this_module, fmt='(5X,A)')
               end do
            else
               call nek_log_warning('helix instance not initialized', module=this_module, fmt='(A)')
            end if
         end procedure forcing_summary
      
      end submodule helix_utils