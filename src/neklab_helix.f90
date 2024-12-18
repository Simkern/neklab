      module neklab_helix
      !---------------------------------------
      !-----     LightKrylov Imports     -----
      !---------------------------------------
         use stdlib_strings, only: padl
         use stdlib_optval, only: optval
      ! Default real kind.
         use LightKrylov, only: dp
         use LightKrylov_Constants, only: imag => one_im_cdp
      ! Logging
         use LightKrylov_Logger
      ! Extensions of the abstract vector types to nek data format.
         use neklab_vectors
         use neklab_nek_forcing, only: neklab_forcing, set_neklab_forcing
         use neklab_nek_setup, only: nek_log_message, nek_log_warning
         implicit none
         include "SIZE"
         include "TOTAL"
         include "ADJOINT"
         include 'RESTART'
         private
         character(len=*), parameter, private :: this_module = 'neklab_helix'
      
         integer, parameter :: lv = lx1*ly1*lz1*lelv
      !! Local number of grid points for the velocity mesh.
         integer, parameter :: lp = lx2*ly2*lz2*lelv
      !! Local number of grid points for the pressure mesh.
         integer, parameter :: nf = 3
      !! Number of forcing components
         integer, parameter :: lbuf = 1
      !! Maximum number of 2d fields to save before outposting

         public :: pipe
         public :: helix_pipe ! constructor for the pipe instance of the helix type

         type, public :: helix
            !! Type containing the basic geometrical and dynamical properties of a helix
            private
            ! geometry inputs
            real(dp) :: length
            real(dp) :: delta
            real(dp) :: diameter
            real(dp) :: pitch_s
            ! derived geometry quantities
            real(dp) :: radius
            real(dp) :: curv_radius ! distance from helix center to pipe center
            real(dp) :: phi         ! Rise angle of helix w.r.t. equatorial plane
            real(dp) :: sweep       ! Total angle swept by helix in streamwise direction
            ! flow
            logical :: if_steady
            real(dp) :: pulse_T
            real(dp) :: omega
            real(dp) :: womersley
            ! forcing
            real(dp), dimension(lv) :: fshape
            real(dp), dimension(nf) :: dpds
            ! mesh inputs
            integer :: nslices
            integer :: nelf
            ! sanity check
            logical :: is_initialized = .false.
            ! data
            ! Sweep angle (radians)
            real(dp), dimension(lv) :: sweep_angle
            ! Angle of the cross-sectional plane around the helix, clockwise from the (positive) y axis
            real(dp), dimension(lv) :: as
            ! Angle within the cross-sectional plane, from the inside to the outside of the helix starting from the negative z axis
            real(dp), dimension(lv) :: alpha 
            ! cylindrical coordinates w.r.t the equatorial plane of the helix (with zax)
            real(dp), dimension(lv) :: ox, oy
            ! cartesian coordinates in torus (without torsion!)
            real(dp), dimension(lv) :: xax, yax, zax
            ! 2D data
            integer :: n2d     ! global number of 2d elements
            integer :: n2d_own ! number of 2d elements owned by current processor
            integer :: nsave   ! buffer fill counter
            integer :: noutc   ! number of data files in cartesian coordinates written to disk
            integer :: noutt   ! number of data files in toroidal coordinates written to disk
            logical :: save_2d_usrt = .true.! save us,ur,ut in addition to vx,vy,vz?
            ! save 2D fields
            logical, dimension(lelv)               :: gowner
            integer, dimension(lelv,2)             :: id2d
            real(dp), dimension(lbuf)              :: t2d
            real(dp), dimension(lx1,ly1,lelv)      :: x2d, y2d
            real(dp), dimension(lx1,ly1,lelv,lbuf) :: vx2d, vy2d, vz2d
         contains
            ! initialization
            procedure, pass(self), public :: init_geom
            procedure, pass(self), public :: init_flow
            ! public computation routines
            procedure, pass(self), public :: compute_fshape
            procedure, pass(self), public :: compute_bf_forcing
            procedure, pass(self), public :: compute_usrt
            procedure, pass(self), public :: compute_ubar
            ! saving 2D data
            procedure, pass(self), public :: save_2d_fields
            procedure, pass(self), public :: compute_2dusrt
            procedure, pass(self), public :: outpost_2d_fields
            procedure, pass(self), public :: load_2d_fields
            ! helper routines
            procedure, pass(self), public :: get_forcing
            procedure, pass(self), public :: get_period
            procedure, pass(self), public :: add_dpds
            procedure, pass(self), public :: is_steady
            procedure, pass(self), public :: setup_summary
            procedure, pass(self), public :: get_fshape
            procedure, pass(self), public :: get_angle_s
            procedure, pass(self), public :: get_alpha
            ! getter/setter for dpds
            procedure, pass(self), public :: get_dpds_all
            procedure, pass(self), public :: get_dpds
            procedure, pass(self), public :: set_dpds
         end type helix

         type(helix) :: pipe
      
      contains

         ! Constructor for the module level instance of helix
         subroutine helix_pipe(delta, diameter, pitch_s, length, nslices, nelf)
            real(dp), intent(in) :: delta
            real(dp), intent(in) :: diameter
            real(dp), intent(in) :: pitch_s
            real(dp), intent(in) :: length
            integer, intent(in) :: nslices
            integer, intent(in) :: nelf
            ! internal
            character(len=128) :: msg

            ! Geometry
            pipe%delta    = delta
            pipe%diameter = diameter
            pipe%pitch_s  = pitch_s
            pipe%length   = length

            ! Mesh specifics
            pipe%nslices  = nslices
            pipe%nelf     = nelf
            
            !  Derived quantities
            pipe%radius      = pipe%diameter/2.0_dp
            pipe%curv_radius = 1.0_dp/pipe%delta
            pipe%phi         = atan2(pipe%pitch_s,pipe%curv_radius)
            pipe%sweep       = pipe%length*cos(pipe%phi)/pipe%curv_radius ! sweep angle in radians

            ! intialize geometry
            call pipe%init_geom()
            ! compute forcing distribution
            call pipe%compute_fshape()

         end subroutine helix_pipe

         subroutine setup_summary(self)
            class(helix), intent(in) :: self
            ! internal
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
               write (msg, '(A,I8)') padl('ntot:', 20), nx1*ny1*nz1*nelv
               call nek_log_message(msg, module=this_module, fmt='(5X,A)')
               write (msg, '(A,I8)') padl('nelv:', 20), nelv
               call nek_log_message(msg, module=this_module, fmt='(5X,A)')
               write (msg, '(A,I8)') padl('slices:', 20), pipe%nslices
               call nek_log_message(msg, module=this_module, fmt='(5X,A)')
               write (msg, '(A,I8)') padl('nel/slice:', 20), pipe%nelf
               call nek_log_message(msg, module=this_module, fmt='(5X,A)')
            else
               call nek_log_warning('helix instance not initialized', module=this_module, fmt='(A)')
            end if
         end subroutine setup_summary

         subroutine init_geom(self)
            class(helix), intent(inout) :: self
            ! internal
            real(dp) :: xmin, xmax, helix_r, s_angle, invnv
            real(dp) :: x_torus, y_torus, z_torus
            real(dp), dimension(lv) :: tmp, pipe_r
            real(dp), dimension(2,lelv) :: xycenter
            !real(dp), dimension(nelv) :: zcenter
            integer :: ix, iy, iz, ie, i, iel, ieg, isl, iseg, nseg
            integer, dimension(lelv) :: islice, isegment, unique_segments, segment_owner
            integer, dimension(lelv) :: segment_count
            integer, dimension(lelv,2) :: cell
            !logical, dimension(lelv) :: gowner
            integer, dimension(lelv) :: lowner
            integer, allocatable :: idx(:)
            logical, allocatable :: segment_found(:)

            real(dp), external :: glmax, glmin

            integer, external :: iglsum

            if (self%is_initialized) call stop_error('Attempting to reinitialize the mesh', this_module, 'init_geom')

         !  Geometry modification for helical pipe

            pi = 4.0_dp*atan(1.0_dp)

            call rescale_x(xm1,-self%radius,self%radius) ! x in [ -r, r ]
            call rescale_x(ym1,-self%radius,self%radius) ! y in [ -r, r ]
            call rescale_x(zm1,0.0_dp,1.0_dp)            ! z in [  0, 1 ]

         !  rotate mesh to set the center of the pipe along x-axis
            call copy(tmp,  xm1, lv)
            call copy(xm1,  zm1, lv)   ! x <--  z
            call copy(zm1, -tmp, lv)   ! z <-- -x
            call copy(self%zax,zm1,lv) ! zax set before curvature in z is added!

            ! rescale the new x axis
            xmin = glmin(xm1,lv)
            xmax = glmax(xm1,lv)
            xm1 = self%sweep/(xmax-xmin) * xm1 ! x in [ 0, max_sweep_angle ]
            call copy(self%sweep_angle,xm1,lv) ! save sweep angle
            call copy(pipe_r,          ym1,lv) ! local distance from pipe center

            ! Sort local elements according to 2D mesh. 
            ! Here we use a trick that relies on the particular structure of meshes extruded
            ! from a 2D mesh using n2to3.
            call izero(islice, nelv)
            call izero(isegment, nelv)
            call izero(isegment, nelv)
            call izero(cell, nelv)
            isl = 0
            call izero(segment_count, nelv)
            allocate(segment_found(self%nelf), source=.false.)
            nseg = 0
            self%n2d_own = 0
            do ie = 1, nelv
               self%gowner(ie) = .false.
               ieg = lglel(ie)
               isl = ieg/self%nelf + 1
               if (mod(ieg,self%nelf)==0) isl = isl - 1
               islice(ie) = isl ! global numbering fills up slices in order
               isegment(ie) = ieg - (isl-1)*self%nelf
               if (isl == 1) then
                  ! the element in the first slice is the global segment owner
                  self%gowner(ie) = .true.
                  self%n2d_own = self%n2d_own + 1
               end if
               !if (nid==i) print *, nid, 'ie, ieg, islice, iseg', ie, ieg, islice(ie), isegment(ie), self%gowner(ie)
               if (.not. segment_found(isegment(ie))) then
                  nseg = nseg + 1
                  unique_segments(nseg) = isegment(ie)
                  segment_owner(nseg) = ie
                  segment_found(isegment(ie)) = .true.
               end if
               idx = findloc(unique_segments, isegment(ie))
               iseg = idx(1)
               segment_count(iseg) = segment_count(iseg) + 1
               ! establish local ownership: the smallest local slice is the local segement owner
               if (islice(ie) < islice(segment_owner(iseg))) segment_owner(iseg) = ie
            end do
            ! update local ownership map for convenience
            do ie = 1, nelv
               idx = findloc(unique_segments, isegment(ie))
               iseg = idx(1)
               lowner(ie) = segment_owner(iseg)
            end do
      
            print '(I2,1X,I3,1X,A7,*(I3,1X))', nid, nseg, 'segs:  ', unique_segments(:nseg)
            print '(I2,1X,I3,1X,A7,*(I3,1X))', nid, nseg, 'counts:', segment_count(:nseg)
            print '(I2,1X,I3,1X,A7,*(I3,1X))', nid, nseg, 'owner: ', segment_owner(:nseg)
            print *, nid, nseg, 'owned: ', self%n2d_own

            do ie = 1, nelv
               call cfill(vz(1,1,1,ie), 1.0_dp*nid, lx1*ly1*lz1)
               call cfill(vx(1,1,1,ie), 1.0_dp*islice(ie),lx1*ly1*lz1)
               call cfill(vy(1,1,1,ie), 1.0_dp*isegment(ie),lx1*ly1*lz1)
               if (lowner(ie) == ie) call cfill(vz(1,1,1,ie), -1.0_dp,lx1*ly1*lz1) 
               if (self%gowner(ie)) call cfill(vz(1,1,1,ie), np*1.0_dp,lx1*ly1*lz1)
            end do
            call outpost(vx,vy,vz,pr,t,'tst')

            ! Extract 2D mesh
            call izero(self%id2d, 2*nelv)
            iel = 0 
		      do ie = 1, nelv
		      	! extract boundary points from the global segment owners
               if (self%gowner(ie)) then
                  iel = iel + 1
		      	   do i = 1, 2*ndim
		      	   	if (cbc(i,ie,1) == 'P') then
		      	   		call ftovec(self%x2d(1,1,iel), zm1, ie, i, nx1, ny1, nz1) ! z --> x
		      	   		call ftovec(self%y2d(1,1,iel), ym1, ie, i, nx1, ny1, nz1)
                        self%id2d(iel,1) = ie
                        self%id2d(iel,2) = i
                     end if
		      	   end do
               end if
		      end do
            ! this is the number of elements that the current processor owns
            call rzero(self%vx2d, nx1*ny1*nelv*lbuf)
            call rzero(self%vy2d, nx1*ny1*nelv*lbuf)
            call rzero(self%vz2d, nx1*ny1*nelv*lbuf)
            call rzero(self%t2d,  lbuf)
            self%nsave = 0
            self%noutc = 0
            self%noutt = 0
            self%n2d   = iglsum(self%n2d_own,1)
            if (self%n2d /= self%nelf) call stop_error('Inconsistent elements in 2D mesh!', module=this_module, procedure='init_geom')

            ! Morph the mesh into a torus
            i = 0
            helix_r = self%curv_radius
            do ie = 1, nelv
            do iz = 1, lz1
            do iy = 1, ly1
            do ix = 1, lx1
               i = i + 1
               self%ox(i)       = helix_r   * sin(self%sweep_angle(i))
               self%oy(i)       = helix_r   * cos(self%sweep_angle(i))
               xm1(ix,iy,iz,ie) = pipe_r(i) * sin(self%sweep_angle(i)) + self%ox(i)
               ym1(ix,iy,iz,ie) = pipe_r(i) * cos(self%sweep_angle(i)) + self%oy(i)
            end do
            end do
            end do
            end do
            call copy(self%xax, xm1, lv) ! xax set before curvature in z is added!
            call copy(self%yax, ym1, lv) ! yax set before curvature in z is added!

            ! Morph the torus into a helix
            if (self%phi /= 0.0_dp) then
               i = 0
               do ie = 1, nelv
               do iz = 1, lz1
               do iy = 1, ly1
               do ix = 1, lx1
                  i = i + 1
                  x_torus = xm1(ix,iy,iz,ie)
                  y_torus = ym1(ix,iy,iz,ie)
                  z_torus = self%zax(i)
                  s_angle = self%sweep_angle(i)
                  xm1(ix,iy,iz,ie) = x_torus - z_torus*sin(self%phi)*cos(s_angle)
                  ym1(ix,iy,iz,ie) = y_torus + z_torus*sin(self%phi)*sin(s_angle)
                  zm1(ix,iy,iz,ie) = s_angle*self%pitch_s + z_torus*cos(self%phi)
               enddo
               enddo
               enddo
               enddo
            end if
            param(59) = 1.   !  All elements deformed

            ! Streamwise angle in the equatorial plane & angle within cross-sectional plane
            self%as    = atan2(self%xax, self%yax) ! clockwise from y axis
            self%alpha = atan2(self%zax, pipe_r)

            self%is_initialized = .true.
            
         end subroutine init_geom

         subroutine init_flow(self, dpds, womersley)
            class(helix), intent(inout) :: self
            real(dp), intent(in) :: dpds(:)
            real(dp), optional, intent(in) :: womersley
            ! internal
            integer :: i, n
            character(len=128) :: msg, fmt
            pi = 4.0_dp*atan(1.0_dp)
            n = size(dpds)
            if (present(womersley)) then
               self%if_steady = .false.
               self%womersley = womersley
               self%omega     = (self%womersley**2)*cpfld(1,1)    ! pulsation frequency
               self%pulse_T   = 2.0d0*pi/self%omega               ! pulsation period
               if (n == 1) then
                  call stop_error('Unsteady case requires more than one forcing component',module=this_module,procedure='init_flow')
               else if (mod(n,2)==0) then
                  call stop_error('Unsteady case requires an uneven number of forcing components',module=this_module,procedure='init_flow')
               end if
            else
               self%if_steady = .true.
               self%womersley = 0.0_dp   
               self%omega     = 0.0_dp   
               self%pulse_T   = 0.0_dp  
               if (n > 1) then
                  call stop_error('Steady case requires only one forcing component',module=this_module,procedure='init_flow')
               end if
            end if
            self%dpds = dpds
            call nek_log_message('Flow parameter initialization:', module=this_module)
            write (msg, '(A,L8)') padl('steady:', 20), self%if_steady
            call nek_log_message(msg, module=this_module, fmt='(5X,A)')
            write (msg, '(A,F15.8)') padl('Wo:', 20), self%womersley
            call nek_log_message(msg, module=this_module, fmt='(5X,A)')
            write (msg, '(A,F15.8)') padl('omega:', 20), self%omega
            call nek_log_message(msg, module=this_module, fmt='(5X,A)')
            write (msg, '(A,F15.8)') padl('T:', 20), self%pulse_T
            call nek_log_message(msg, module=this_module, fmt='(5X,A)')
            write (msg, '(A,F15.8,1X,F15.8)') padl('dpds_00:', 20), self%dpds(1), 0.0_dp
            call nek_log_message(msg, module=this_module, fmt='(5X,A)')
            do i = 2, n, 2
               write(fmt,'("dpds_",I2.2)') i/2
               print *, fmt
               write (msg, '(A,F15.8,1X,F15.8)') padl(trim(fmt), 20), self%dpds(i), self%dpds(i+1)
               call nek_log_message(msg, module=this_module, fmt='(5X,A)')
            end do
                  
         end subroutine init_flow

         subroutine compute_fshape(self)
            class(helix), intent(inout) :: self
            ! internals
            integer :: i, ix, iy, iz, ie
            real(dp) :: helix_r2, r, rr, alpha

            self%fshape = 0.0_dp

            i = 0
            do ie = 1, lelv
            do iz = 1, lz1
            do iy = 1, ly1
            do ix = 1, lx1
               i = i+1

               ! Distance from the origin in the equatorial plane
               helix_r2 = self%xax(i)**2 + self%yax(i)**2
               ! Distance from the pipe center in the equatorial plane
               r = sqrt(helix_r2) - self%curv_radius
               ! Azimuthal angle in the cross-sectional plane
               alpha = atan2(r, self%zax(i))
               ! Radial position in the cross-sectional plane
               rr = sqrt(r**2 + self%zax(i)**2)

               ! Compute fshape
               self%fshape(i) = 1.0_dp / abs(1.0_dp + self%delta * rr * sin(alpha))
            end do
            end do
            end do
            end do
         end subroutine compute_fshape

         real(dp) pure function get_forcing(self, t) result(f)
            class(helix), intent(in) :: self
            real(dp), intent(in) :: t
            !! time
            ! internal
            integer :: i
            complex(dp) :: eiwt, dpds

            f = self%dpds(1)
            if (.not.self%if_steady) then
               eiwt = cexp(imag * self%omega * t)
               do i = 2, nf, 2
                  dpds = self%dpds(i) + imag*self%dpds(i+1)
                  f = f + 2.0_dp * real(dpds * eiwt)
               end do
            end if
         end function get_forcing

         real(dp) pure function get_period(self) result(T)
            class(helix), intent(in) :: self
            T = self%pulse_T
         end function get_period

         subroutine compute_bf_forcing(self, t)
            class(helix), intent(in) :: self
            real(dp) :: t
            !! time
            ! internal
            integer :: i
            real(dp) :: fs
            real(dp), dimension(lv) :: ffx, ffy, ffz
            fs = self%get_forcing(t) / self%curv_radius

            do i = 1, lv
               ffx(i) =  fs * self%fshape(i) * cos(self%phi) * cos(self%as(i))
               ffy(i) = -fs * self%fshape(i) * cos(self%phi) * sin(self%as(i))
               ffz(i) =  fs * self%fshape(i) * sin(self%phi)
            end do

            ! set baseflow forcing
            call set_neklab_forcing(ffx, ffy, ffz, ipert=0)

         end subroutine compute_bf_forcing

         subroutine compute_usrt(self, u, v, w, us, ur, ut)
            class(helix), intent(in) :: self
            real(dp), dimension(lx1,ly1,lz1,lelv), intent(in) :: u
            real(dp), dimension(lx1,ly1,lz1,lelv), intent(in) :: v
            real(dp), dimension(lx1,ly1,lz1,lelv), intent(in) :: w
            real(dp), dimension(lx1,ly1,lz1,lelv), intent(out) :: us
            real(dp), dimension(lx1,ly1,lz1,lelv), intent(out) :: ur
            real(dp), dimension(lx1,ly1,lz1,lelv), intent(out) :: ut
            ! internal
            real(dp) :: x, y, z, phi, a, s
            integer :: ix, iy, iz, ie, i
            real(dp) :: utmp, vtmp
            phi = self%phi
            do ie = 1, nelv
            do iz = 1, lz1
            do iy = 1, ly1
            do ix = 1, lx1
               i = i+1
               x = self%xax(i)
               y = self%yax(i)
               z = self%zax(i)
               s = self%as(i)
               a = self%alpha(i)
               us(ix,iy,iz,ie) = cos(phi) * (
     $                          + cos(s)   * u(ix,iy,iz,ie)
     $                          - sin(s)   * v(ix,iy,iz,ie) )
     $                          + sin(phi) * w(ix,iy,iz,ie)
               utmp            = sin(s)    * u(ix,iy,iz,ie)
     $                          + cos(s)   * v(ix,iy,iz,ie)
               vtmp            = sin(phi) * (
     $                          - cos(s)   * u(ix,iy,iz,ie)
     $                          - sin(s)   * v(ix,iy,iz,ie) )
     $                          + cos(phi) * w(ix,iy,iz,ie)
               ur(ix,iy,iz,ie) = cos(a) * utmp + sin(a) * vtmp
               ut(ix,iy,iz,ie) = sin(a) * utmp - cos(a) * vtmp
            end do
            end do
            end do
            end do
         end subroutine compute_usrt
            
         real(dp) function compute_ubar(self,u,v,w) result(ubar)
            class(helix) :: self
            real(dp), dimension(lx1,ly1,lz1,lelv) :: u
            real(dp), dimension(lx1,ly1,lz1,lelv) :: v
            real(dp), dimension(lx1,ly1,lz1,lelv) :: w
            ! internal
            integer :: ix, iy, iz, ie, i
            real(dp) :: num, den, us, us_r
            real(dp), external :: glsum

            pi = atan(1.0_dp)*4.0_dp

            i = 0
            num = 0.0_dp
            den = 0.0_dp
            do ie = 1, nelv
            do iz = 1, lz1
            do iy = 1, ly1
            do ix = 1, lx1
               i = i+1
               us = cos(self%phi)*(
     $               cos( self%as(i) ) * u(ix,iy,iz,ie)
     $             - sin( self%as(i) ) * v(ix,iy,iz,ie))
     $             + sin( self%phi ) * w(ix,iy,iz,ie)
               us_r = us * self%fshape(i) ! u/r
               num = num + us_r*bm1(ix,iy,iz,ie)
               den = den + self%fshape(i)*bm1(ix,iy,iz,ie)
            end do
            end do
            end do
            end do
            num = glsum(num,1)
            den = glsum(den,1)
            ubar = num/den  ! "1/r"-weighted volumetric average of streamwise velocity
            
         end function compute_ubar

         subroutine save_2d_fields(self, u, v, w)
            class(helix), intent(inout) :: self
            real(dp), dimension(lx1,ly1,lz1,lelv), intent(in) :: u
            real(dp), dimension(lx1,ly1,lz1,lelv), intent(in) :: v
            real(dp), dimension(lx1,ly1,lz1,lelv), intent(in) :: w
            ! internal
            integer :: i, iel, ifc
            self%nsave = self%nsave + 1
            ! save data to buffer
            do i = 1, self%n2d_own
               iel = self%id2d(i, 1)
               ifc = self%id2d(i, 2)
               call ftovec(self%vx2d(1,1,i,self%nsave), u, iel, ifc, nx1, ny1, nz1)
		      	call ftovec(self%vy2d(1,1,i,self%nsave), v, iel, ifc, nx1, ny1, nz1)
		      	call ftovec(self%vz2d(1,1,i,self%nsave), w, iel, ifc, nx1, ny1, nz1)
            end do
            self%t2d(self%nsave) = time
            ! save data to file when buffer is full
            if (self%nsave == lbuf .or. lastep == 1) then
               call self%outpost_2d_fields('c')
               if (self%save_2d_usrt) then
                  call self%compute_2dusrt() ! self%v[xyz]2d are overwritten
                  call self%outpost_2d_fields('t')
               end if
               self%nsave = 0
            end if
         end subroutine save_2d_fields

         subroutine compute_2dusrt(self)
            ! this routine will overwrite self%v[xyz]2d
            class(helix), intent(inout) :: self
            ! internal
            integer, parameter :: iz = 1
            real(dp) :: x, y, z, phi, sweep, a, s
            integer :: ix, iy, ie, iel, ibuf, i
            real(dp) :: utmp, vtmp, ux, uy, uz
            phi = self%phi
            do iel = 1, self%n2d_own
               ie = self%id2d(iel, 1)
               ! add offset
               i = (ie-1)*lx1*ly1*lz1
               ! the xy plane corresponds to lz1 = 1 so nothing needs to be done for z
               do iy = 1, ly1
               do ix = 1, lx1
                  i = i+1
                  x = self%xax(i)
                  y = self%yax(i)
                  z = self%zax(i)
                  s = self%as(i)
                  a = self%alpha(i)
                  ! iterate over buffer
                  do ibuf = 1, lbuf
                     ux = self%vx2d(ix,iy,iel,ibuf)
                     uy = self%vy2d(ix,iy,iel,ibuf)
                     uz = self%vz2d(ix,iy,iel,ibuf)
                     ! overwrite v[xyz]2d with u[srt]2d
                     self%vx2d(ix,iy,iel,ibuf) = cos(phi)*(
     $                                          + cos(s) * ux 
     $                                          - sin(s) * uy )
     $                                          + sin(phi) * uz
                     utmp                      = sin(s) * ux 
     $                                          + cos(s) * uy
                     vtmp                      = sin(phi) * (
     $                                          - cos(s) * ux 
     $                                          - sin(s) * uy )
     $                                          + cos(phi) * uz
                     self%vy2d(ix,iy,iel,ibuf) = cos(a)*utmp
     $                                          + sin(a)*vtmp
                     self%vz2d(ix,iy,iel,ibuf) = sin(a)*utmp
     $                                          - cos(a)*vtmp
                  end do ! lbuf
               end do    ! lx1
               end do    ! ly1
            end do       ! self%n2d_own
         end subroutine compute_2dusrt
         
         subroutine outpost_2d_fields(self, iname)
            class(helix), intent(inout) :: self
            character(len=1), intent(in) :: iname
            ! internals
            character(len=128)  :: fname
            character(len=1024) :: head, ftm
            real rtmpv1(lx1*ly1*lelv), rtmpv(lx1*ly1*lelv)
            real*4 rtmpv2(2*lx1*ly1*lelv)
            equivalence (rtmpv1,rtmpv2)
            integer :: ierr, itmp, i, nxy, ip, ibuf, iel, len, i_own
            integer, allocatable :: n2d_own(:)
            integer, allocatable :: n2d_elmap(:)
            integer :: isend(lelv)
            integer :: wdsl, isl
            real*4 test
            parameter (test=6.54321)
            nxy = lx1*ly1
            wdsl = wdsize/4
            isl  = isize/4
            write(fname,'(A,I3.3,A)') iname, self%noutc, '.fld'

            if (nid == 0) then
               call byte_open(fname, ierr)

               ! write file's header
               ftm="('#tor',1x,i1,1x,'(lx1, ly1 =',2i9,"//          ! 18 + 18
     $         "') (nelf =',i9,"//             ! 9 + 9
     $         "') (time =',e17.9,"//          ! 9 + 17
     $         "') (nsave, lbuf =', 2i9,')')"  ! 16 + 18 + 1(1) = 116
               write(head,ftm) wdsize,lx1,ly1,self%nelf,time,self%nsave,lbuf
               call byte_write(head,116/4,ierr)

               ! write big/little endian test
               call byte_write(test,1,ierr)

               ! write metadata
               call byte_write(lx1,isl,ierr)
               call byte_write(ly1,isl,ierr)
               call byte_write(self%nelf,isl,ierr)
               call byte_write(time,wdsl,ierr)
               call byte_write(self%nsave,isl,ierr)
               call byte_write(lbuf,isl,ierr)
            end if

            ! gather information about elements on other procs  
            allocate(n2d_own(np))           ! number of elements owned by each proc
            call izero(n2d_own,np)
            allocate(n2d_elmap(self%nelf))  ! global element number of owned elements
            call izero(n2d_elmap,self%nelf)
            ! determine how many elements to dump
            if (nid == 0) then
               ! first for the master node 
               n2d_own(1) = self%n2d_own
               do i = 1, self%n2d_own
                  n2d_elmap(i) = lglel(self%id2d(i,1))     ! get global element number
               end do 
               iel = self%n2d_own
               ! then gather info from other procs
               do ip = 1, np-1
                  call csend(ip,itmp,isize,ip,0)           ! hand shake
                  call crecv(ip,i_own,isize)               ! recv number of elements
                  n2d_own(ip+1) = i_own
                  call crecv(ip,isend(:i_own),i_own*isize) ! recv global element map
                  n2d_elmap(iel+1:iel+i_own) = isend(:i_own)
                  iel = iel + i_own
               enddo
               if (iel /= self%nelf) then
                  call stop_error('Not all elements in slice found!', this_module, 'outpost_2d_fields')
               end if
               !print '(A,*(1X,I0))', 'elmap:', n2d_elmap
               ! write it to file
               call byte_write(n2d_elmap,self%nelf*isl,ierr)
            else
               call crecv(nid,itmp,isize)                  ! hand shake
               call csend(nid,self%n2d_own,isize,0,0)      ! send number of elements
               len = self%n2d_own
               do i = 1, self%n2d_own
                  isend(i) = lglel(self%id2d(i,1))       
               end do
               call csend(nid,isend(:len),len*isize,0,0)   ! send global element map
            endif
            call bcast(n2d_own, np*isize) ! broadcast to all procs
            call bcast(n2d_elmap, self%nelf*isize) ! broadcast to all procs

            !if (nid == 0) print *, 'nid', nid, ':', np
            !print *, 'nid', nid, ':', self%n2d_own, n2d_own

            ! coordinates
            call gather_and_write_slice(n2d_own, self%x2d)
            call gather_and_write_slice(n2d_own, self%y2d)
            ! velocity data
            do ibuf = 1, self%nsave
               call gather_and_write_slice(n2d_own, self%vx2d(:,:,:,ibuf))
               call gather_and_write_slice(n2d_own, self%vy2d(:,:,:,ibuf))
               call gather_and_write_slice(n2d_own, self%vz2d(:,:,:,ibuf))
            end do
            ! master closes the file
            if (nid == 0) then
               call byte_close(ierr)
            endif
            !call nek_end()
         end subroutine outpost_2d_fields

         subroutine gather_and_write_slice(n2d_own, slicedata)
            integer, intent(in) :: n2d_own(:)
            real(dp), intent(in) :: slicedata(:,:,:)
            ! internal
            integer :: nxy, itmp, wdsl, isl, len, ierr, ip
            real rtmpv1(lx1*ly1*lelv), rtmpv(lx1*ly1*lelv)
            real*4 rtmpv2(2*lx1*ly1*lelv)
            equivalence (rtmpv1,rtmpv2)
            nxy  = lx1*ly1
            wdsl = wdsize/4
            isl  = isize/4
            if (nid == 0) then
               ! master writes if there are data
               len = nxy*n2d_own(nid+1)
               if (wdsl.eq.2) then
                  call copy(rtmpv1,slicedata,len)
                  call byte_write(rtmpv2,len*wdsl,ierr)
               else
                  call copyX4(rtmpv2,slicedata,len)
                  call byte_write(rtmpv2,len,ierr)
               end if
               ! get data from other procs and write to file
               do ip = 1, np-1
                  len = nxy*n2d_own(ip+1)
                  call csend(ip,itmp,isize,ip,0) ! hand shake
                  call crecv2(ip,rtmpv,len*wdsize,ip)
                  ! write data
                  if (wdsl.eq.2) then
                     call copy(rtmpv1,rtmpv,len)
                     call byte_write(rtmpv2,len*wdsl,ierr)
                  else
                     call copyX4(rtmpv2,rtmpv,len)
                     call byte_write(rtmpv2,len,ierr)
                  endif
               end do
            else 
               ! send data to master
               call crecv2(nid,itmp,isize,0) ! hand shake
               len = nxy*n2d_own(nid+1)
               call csend(nid,slicedata,len*wdsize,0,itmp)
            end if
         end subroutine gather_and_write_slice

         subroutine load_2d_fields(self, idx)
            ! only nid 0 will read
            class(helix), intent(inout) :: self
            integer, intent(in) :: idx
            ! internal
            integer ierr, hdrsize
            real*4 test_pattern
            integer :: nxr, nyr, nelfr, nsaver, lbufr, wdsizer
            integer :: wdsl, isl, itmp, ip, nxy, len, i, lfield, ie, iel
            real rtmpv1(lx1*ly1*lelv), rtmpv(lx1*ly1)
            integer, allocatable :: global_map(:)
            real(dp) :: timer
            real(dp), allocatable :: fldata(:,:,:,:)
            real(dp), allocatable :: eldata(:,:)
            character(len=132) :: hdr
            character(len=132) :: fname
            character(len=4)   :: sdummy
            character(len=10)  :: rdc
            integer, allocatable :: n2d_own(:)
            integer, allocatable :: n2d_elmap(:)
            logical :: logend
            ! functions
            logical, external :: if_byte_swap_test
            integer, external :: igl_running_sum
            write(fname,'("c",I3.3,".fld")') idx
            print *, nid, fname
            hdrsize = 116
            nxy = lx1*ly1
            allocate(global_map(self%nelf))
            if (nid == 0) then
               call byte_open(fname,ierr)
               ! read header
               if (ierr == 0) then
                  call blank     (hdr,hdrsize)
                  call byte_read (hdr,hdrsize/4,ierr)
               endif
               if (ierr == 0) then
                  call byte_read (test_pattern,1,ierr)
                  if_byte_sw = if_byte_swap_test(test_pattern,ierr) ! determine endianess
               endif
               print *, 'header:', hdr
               read(hdr,*) sdummy, wdsizer
               wdsl = wdsizer/4
               isl  = isize/4
               call byte_read(nxr, isl, ierr)
               call byte_read(nyr, isl, ierr)
               call byte_read(nelfr, isl, ierr)
               call byte_read(timer, wdsl, ierr)
               call byte_read(nsaver, isl, ierr)
               call byte_read(lbufr, isl, ierr)
               print *, 'data:', nxr, nyr, nelfr, timer, nsaver, lbufr
               ! read global element mapping
               call byte_read(global_map, self%nelf*isl, ierr)
            end if
            !call bcast(n2d_own, self%nelf*isize) ! broadcast to all procs
            !print *, 'n2d_own:', n2d_own


            !!!!! CHECK IF THIS IS THE INFO I NEED

            ! gather information about elements on other procs  
            !allocate(n2d_own(np))           ! number of elements owned by each proc
            !call izero(n2d_own,np)
            !allocate(n2d_elmap(self%nelf))  ! global element number of owned elements
            !call izero(n2d_elmap,self%nelf)
            !! determine how many elements to dump
            !if (nid == 0) then
            !   ! first for the master node 
            !   n2d_own(1) = self%n2d_own
            !   do i = 1, self%n2d_own
            !      n2d_elmap(i) = lglel(self%id2d(i,1))     ! get global element number
            !   end do 
            !   iel = self%n2d_own
            !   ! then gather info from other procs
            !   do ip = 1, np-1
            !      call csend(ip,itmp,isize,ip,0)           ! hand shake
            !      call crecv(ip,i_own,isize)               ! recv number of elements
            !      n2d_own(ip+1) = i_own
            !      call crecv(ip,isend(:i_own),i_own*isize) ! recv global element map
            !      n2d_elmap(iel+1:iel+i_own) = isend(:i_own)
            !      iel = iel + i_own
            !   enddo
            !   if (iel /= self%nelf) then
            !      call stop_error('Not all elements in slice found!', this_module, 'outpost_2d_fields')
            !   end if
            !   !print '(A,*(1X,I0))', 'elmap:', n2d_elmap
            !   ! write it to file
            !   call byte_write(n2d_elmap,self%nelf*isl,ierr)
            !else
            !   call crecv(nid,itmp,isize)                  ! hand shake
            !   call csend(nid,self%n2d_own,isize,0,0)      ! send number of elements
            !   len = self%n2d_own
            !   do i = 1, self%n2d_own
            !      isend(i) = lglel(self%id2d(i,1))       
            !   end do
            !   call csend(nid,isend(:len),len*isize,0,0)   ! send global element map
            !endif
            !call bcast(n2d_own, np*isize) ! broadcast to all procs
            !call bcast(n2d_elmap, self%nelf*isize) ! broadcast to all procs

            allocate(fldata(lx1,ly1,self%nelf,lbufr))
            allocate(eldata(lx1,ly1))
            lfield = nxy*self%nelf
            if (nid == 0) then
               ! read x
               call byte_read(fldata, lfield, ierr)
               ! read y
               call byte_read(fldata, lfield, ierr)
            end if
            do i = 1, nsaver
               do iel = 1, self%nelf
                  ie = gllel(global_map(iel)) ! get local element number
                  if (nid == 0) then
                     ! read vx
                     call byte_read(eldata, nxy, ierr)
                  end if
                  ! distribute
                  if (nid == 0) then
                     if (self%gowner(ie)) then
                        print *, 'nid',nid, ': element ', iel, ', global:', global_map(iel), 'local:', ie
                        call copy4r(self%vx2d(1,1,ie,i),eldata,nxy)
                     else
                        !
                        ! send to owner
                        !
                        ! nid not correct: find better
                        !
                        call crecv(nid,itmp,isize) ! hand shake
                        call csend(nid,eldata,nxy*wdsize,0,0)
                     end if
                  else 
                     if (self%gowner(ie)) then
                        ! get from master
                        print *, 'nid',nid, ': element ', iel, ', global:', global_map(iel), 'local:', ie
                        call csend(0,itmp,isize,ip,0) ! hand shake
                        call crecv2(0,rtmpv,nxy*wdsize,0)
                        call copy4r(self%vx2d(1,1,ie,i),rtmpv,nxy)
                     end if
                  end if
               end do
               ! vy
               ! vz
            end do ! 1, nsaver
            ! close the file
            if (nid == 0) call byte_close(ierr)
            if (ierr /= 0) call stop_error('Error reading file '//trim(fname), procedure='load')
         end subroutine load_2d_fields

         logical pure function is_steady(self) result(steady)
            class(helix), intent(in) :: self
            steady = self%if_steady
         end function is_steady

         subroutine set_dpds(self, dpds, reset)
            class(helix), intent(inout) :: self
            real(dp), dimension(nf), intent(in) :: dpds
            logical, optional, intent(in) :: reset
            ! internal
            logical :: reset_dpds
            reset_dpds = optval(reset, .true.)
            if (reset_dpds) self%dpds = 0.0_dp
            self%dpds = self%dpds + dpds
         end subroutine set_dpds

         subroutine get_dpds_all(self, dpds)
            class(helix), intent(in) :: self
            real(dp), dimension(nf), intent(out) :: dpds
            dpds = self%dpds
         end subroutine get_dpds_all

         real(dp) pure function get_dpds(self, i) result(dpds_i)
            class(helix), intent(in) :: self
            integer, intent(in) :: i
            dpds_i = self%dpds(i)
         end function get_dpds

         subroutine add_dpds(self, df, i)
            class(helix), intent(inout) :: self
            real(dp), intent(in) :: df
            integer, intent(in) :: i
            self%dpds(i) = self%dpds(i) + df
         end subroutine add_dpds

         subroutine get_fshape(self, fshape)
            class(helix), intent(in) :: self
            real, dimension(lx1,ly1,lz1,lelv), intent(out) :: fshape
            call copy(fshape, self%fshape, lv)
         end subroutine get_fshape

         subroutine get_angle_s(self, angle_s)
            class(helix), intent(in) :: self
            real, dimension(lx1,ly1,lz1,lelv), intent(out) :: angle_s
            call copy(angle_s, self%as, lv)
         end subroutine get_angle_s

         subroutine get_alpha(self, alpha)
            class(helix), intent(in) :: self
            real, dimension(lx1,ly1,lz1,lelv), intent(out) :: alpha
            call copy(alpha, self%alpha, lv)
         end subroutine get_alpha
      
      end module neklab_helix
