      module neklab_helix
      !---------------------------------------
      !-----     LightKrylov Imports     -----
      !---------------------------------------
         use stdlib_strings, only: padl
         use stdlib_optval, only: optval
         use stdlib_logger, only: debug_level
      ! Default real kind.
         use LightKrylov, only: dp
         use LightKrylov_Constants, only: imag => one_im_cdp
      ! Logging
         use LightKrylov_Logger
      ! Extensions of the abstract vector types to nek data format.
         use neklab_vectors
         use neklab_nek_forcing, only: neklab_forcing, set_neklab_forcing
         use neklab_nek_setup, only: nek_log_message, nek_log_information, nek_log_warning, nek_log_debug
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
            integer, dimension(lelv)               :: lowner
            logical, dimension(lelv), public               :: gowner
            integer, dimension(lelv,2)             :: id2d
            integer, dimension(lelv)               :: global2local
            real(dp), dimension(lbuf)              :: t2d
            real(dp), dimension(lx1,ly1,lelv)      :: x2d, y2d
            real(dp), dimension(lx1,ly1,lelv,lbuf), public :: vx2d, vy2d, vz2d
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
            integer :: ix, iy, iz, ie, i, iel, ieg, iseg, nseg, iface, isl, level
            integer, dimension(lelv) :: islice, isegment
            integer, dimension(:), allocatable :: unique_segments, segment_owner, segment_count
            integer, allocatable :: idx(:) ! for findloc
            logical, allocatable :: segment_found(:)
            ! functions
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
            allocate(unique_segments(self%nelf)); call izero(unique_segments, self%nelf)
            allocate(segment_count  (self%nelf)); call izero(segment_count,   self%nelf)
            allocate(segment_owner  (self%nelf)); call izero(segment_owner,   self%nelf)
            allocate(segment_found  (self%nelf), source=.false.)
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
               self%lowner(ie) = segment_owner(iseg)
            end do
      
            call logger%configuration(level=level)
            if (level <= 20) then
               print '(A,2(I0,1X),A,*(1X,I0))', 'DEBUG 2dmap: ', nid, nseg, 'unique streamwise segments: ', unique_segments(:nseg)
               print '(A,2(I0,1X),A,*(1X,I0))', 'DEBUG 2dmap: ', nid, nseg, '# of element in segment:    ', segment_count(:nseg)
               print '(A,2(I0,1X),A,*(1X,I0))', 'DEBUG 2dmap: ', nid, nseg, 'local element segment owner:', segment_owner(:nseg)
               print '(A,2(I0,1X),A,I0)'      , 'DEBUG 2dmap: ', nid, nseg, 'owned: ', self%n2d_own
               call nekgsync()
               do ie = 1, nelv
                  call cfill(vz(1,1,1,ie), 1.0_dp*nid, lx1*ly1*lz1)
                  call cfill(vx(1,1,1,ie), 1.0_dp*islice(ie),lx1*ly1*lz1)
                  call cfill(vy(1,1,1,ie), 1.0_dp*isegment(ie),lx1*ly1*lz1)
                  if (self%lowner(ie) == ie) call cfill(vz(1,1,1,ie), -1.0_dp,lx1*ly1*lz1) 
                  if (self%gowner(ie)) call cfill(vz(1,1,1,ie), np*1.0_dp,lx1*ly1*lz1)
               end do
               call outpost(vx,vy,vz,pr,t,'m2d')
               print '(A,A3,3X,3(A5),2X,2(A20))', 'DEBUG: 2dmap ', 'nid', 'glob', 'locl', 'l2d', 'xavg', 'yavg'
            end if

            ! Extract 2D mesh
            call izero(self%id2d, 2*nelv)
            call ifill(self%global2local, -1, nelv)
            iel = 0 
		      do ie = 1, nelv
		      	! extract boundary points from the global segment owners
               if (self%gowner(ie)) then
                  iel = iel + 1
		      	   do iface = 1, 2*ndim
		      	   	if (cbc(iface,ie,1) == 'P') then
		      	   		call ftovec(self%x2d(1,1,iel), zm1, ie, iface, nx1, ny1, nz1) ! z --> x
		      	   		call ftovec(self%y2d(1,1,iel), ym1, ie, iface, nx1, ny1, nz1)
                        self%id2d(iel,1) = ie
                        self%id2d(iel,2) = iface
                        self%global2local(lglel(ie)) = ie
                     end if
		      	   end do
                  if (level <= 20) then
                     print '(A,I3,A,3(1X,I4),A,3X,F17.8,3x,F17.8)', 'DEBUG 2dmap: ', nid, ' el', lglel(ie), ie, iel, ': ', 
     &                           sum(self%x2d(:,:,iel))/(lx1*ly1), sum(self%y2d(:,:,iel))/(lx1*ly1)
                  end if
               end if
		      end do
            if (level <= 20) then
               call nekgsync()
               print '(A,I3,A,*(1x,I0))', 'DEBUG 2dmap: ', nid, ', global2local: ', self%global2local(:nelv)
            end if
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
            integer :: i, iel, ifc, level
            character(len=128) :: msg
            real(dp) :: xavg, yavg, vxavg, vyavg, vzavg
            call logger%configuration(level=level)
            self%nsave = self%nsave + 1
            ! save data to buffer
            write(msg,'(A,I3,A,I3)') 'Save 2D data: ', self%nsave, '/', lbuf
            call nek_log_message(msg, this_module, 'save_2d_fields')
            do i = 1, self%n2d_own
               iel = self%id2d(i, 1)
               ifc = self%id2d(i, 2)
               call ftovec(self%vx2d(1,1,i,self%nsave), u, iel, ifc, nx1, ny1, nz1)
		      	call ftovec(self%vy2d(1,1,i,self%nsave), v, iel, ifc, nx1, ny1, nz1)
		      	call ftovec(self%vz2d(1,1,i,self%nsave), w, iel, ifc, nx1, ny1, nz1)
               if (level <= 20) then
                  xavg = sum(self%x2d(:,:,iel))/(lx1*ly1)
                  yavg = sum(self%y2d(:,:,iel))/(lx1*ly1)
                  vxavg = sum(self%vx2d(:,:,iel,self%nsave))/(lx1*ly1)
                  vyavg = sum(self%vy2d(:,:,iel,self%nsave))/(lx1*ly1)
                  vzavg = sum(self%vz2d(:,:,iel,self%nsave))/(lx1*ly1)
                  print '(A,I8,I8,A,5(3X,F16.8))', 'DEBUG: save el', i, iel, ': ', xavg, yavg, vxavg, vyavg, vzavg
               end if
            end do
            self%t2d(self%nsave) = time
            ! save data to file when buffer is full
            if (self%nsave == lbuf .or. lastep == 1) then
               self%noutc = self%noutc + 1
               call self%outpost_2d_fields('c')
               if (self%save_2d_usrt) then
                  call self%compute_2dusrt() ! self%v[xyz]2d are overwritten
                  self%noutt = self%noutt + 1
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
            character(len=128)  :: fname, msg
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
            write(fname,'(A,A,I3.3,A)') iname, '2dtorus', self%noutc, '.fld'
            write(msg,'(A,I3,4X,A,A)') 'Outpost 2D data: ', self%nsave, 'fname: ', trim(fname)
            call nek_log_message(msg, this_module, 'outpost_2d_fields')
            if (nid == 0) then
               call byte_open(fname, ierr)
               if (ierr /= 0) call stop_error('Error opening file '//trim(fname), procedure='outpost_2d_fields')

               ! write file's header
               ftm="('#tor',1x,i1,1x,'(lx1, ly1 =',2i9,') (nelf =',i9,') (time =',e17.9,') (nsave, lbuf =', 2i9,')')"
               write(head,ftm) wdsize,lx1,ly1,self%nelf,time,self%nsave,lbuf
               call byte_write(head,116/4,ierr)
               if (ierr /= 0) call stop_error('Error writing header in file '//trim(fname), procedure='outpost_2d_fields')  

               ! write big/little endian test
               call byte_write(test,1,ierr)

               ! write metadata
               call byte_write(lx1,isl,ierr)
               call byte_write(ly1,isl,ierr)
               call byte_write(self%nelf,isl,ierr)
               call byte_write(time,wdsl,ierr)
               call byte_write(self%nsave,isl,ierr)
               call byte_write(lbuf,isl,ierr)
               if (ierr /= 0) call stop_error('Error writing metadata in file '//trim(fname), procedure='outpost_2d_fields')
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
               if (iel /= self%nelf) call stop_error('Not all elements in slice found!', this_module, 'outpost_2d_fields')
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
            call bcast(n2d_own, np*isize)          ! broadcast to all procs
            call bcast(n2d_elmap, self%nelf*isize) ! broadcast to all procs

            ! coordinates
            call nek_log_message('   '//trim(fname)//': write x2d ...', this_module, 'outpost_2d_fields')
            call gather_and_write_slice(self%x2d, n2d_own)
            call nek_log_message('   '//trim(fname)//': write y2d ...', this_module, 'outpost_2d_fields')
            call gather_and_write_slice(self%y2d, n2d_own)
            ! velocity data
            do i = 1, self%nsave
               write(msg,'(3X,A,A,1X,I3)') trim(fname),': write v[xyz]2d', i
               call nek_log_message(msg, this_module, 'outpost_2d_fields')
               call gather_and_write_slice(self%vx2d(:,:,:,i), n2d_own)
               call gather_and_write_slice(self%vy2d(:,:,:,i), n2d_own)
               call gather_and_write_slice(self%vz2d(:,:,:,i), n2d_own)
            end do
            ! master closes the file
            if (nid == 0) then 
               call byte_close(ierr)
               if (ierr /= 0) call stop_error('Error closing file '//trim(fname), procedure='outpost_2d_fields')
            end if
         end subroutine outpost_2d_fields

         subroutine load_2d_fields(self, idx)
            implicit none
            ! only nid 0 will read
            class(helix), intent(inout) :: self
            integer, intent(in) :: idx
            ! internal
            integer ierr, hdrsize
            real*4 test_pattern
            integer :: nxr, nyr, nelfr, nsaver, lbufr, wdsizr
            integer :: wdsl, isl, itmp, ip, nxy, i, ie, ieg, iel, nelf
            real rtmpv(lx1*ly1)
            integer, allocatable :: global_map(:)
            real(dp) :: timer
            real(dp) :: eldata(lx1,ly1)
            character(len=132) :: hdr
            character(len=132) :: fname
            character(len=4)   :: sdummy
            character(len=128) :: msg
            common /CTMP1/ fldum(lx1*ly1*lelv)
            real fldum
            ! functions
            logical, external :: if_byte_swap_test
            write(fname,'("c2dtorus",I3.3,".fld")') idx
            call nek_log_information('Load 2D data from file '//trim(fname), this_module, 'load_2d_fields')
            hdrsize = 116
            nxy = lx1*ly1
            nelf = self%nelf
            allocate(global_map(nelf))
            if (nid == 0) then
               call byte_open(fname,ierr)
               if (ierr /= 0) call stop_error('Error opening file '//trim(fname), procedure='load_2d_fields')
               ! read header
               if (ierr == 0) then
                  call blank     (hdr,hdrsize)
                  call byte_read (hdr,hdrsize/4,ierr)
               endif
               if (ierr == 0) then
                  call byte_read (test_pattern,1,ierr)
                  if_byte_sw = if_byte_swap_test(test_pattern,ierr) ! determine endianess
               endif
               call nek_log_debug('header: '//trim(hdr), this_module, 'load_2d_fields')
               ! read wdsize from header
               read(hdr,*) sdummy, wdsizr
               wdsl = wdsizr/4
               isl  = isize/4
               ! read metadata
               call byte_read(nxr,    isl, ierr)
               call byte_read(nyr,    isl, ierr)
               call byte_read(nelfr,  isl, ierr)
               call byte_read(timer, wdsl, ierr)
               call byte_read(nsaver, isl, ierr)
               call byte_read(lbufr,  isl, ierr)
               write(msg,'(A,3(1X,I0),1X,E15.7,2(1X,I0))') 'metadata: ', nxr, nyr, nelfr, timer, nsaver, lbufr
               call nek_log_debug(msg, this_module, 'load_2d_fields')
               ! read global element mapping
               call byte_read(global_map, nelf*isl, ierr)
               if (ierr /= 0) call stop_error('Error reading gloabl element map from file '//trim(fname), procedure='load_2d_fields')
               ! read coords but skip them
               call byte_read(fldum, nxy*nelf*wdsl, ierr)
               call byte_read(fldum, nxy*nelf*wdsl, ierr)
               if (ierr /= 0) call stop_error('Error reading coordinates from file '//trim(fname), procedure='load_2d_fields')
            end if
            call bcast(nsaver, isize)               ! broadcast number of saved snapshots
            call bcast(global_map, nelf*isize) ! broadcast global element map
            call bcast(wdsizr, isize)               ! broadcast word size
            ! load data
            do i = 1, nsaver
               call nek_log_information('    Read vx2d ...', this_module, 'load_2d_fields')
               call load_and_distribute_slice(self%vx2d(:,:,:,i), global_map, self%global2local, nelf, if_byte_sw)
               call nek_log_information('    Read vy2d ...', this_module, 'load_2d_fields')
               call load_and_distribute_slice(self%vy2d(:,:,:,i), global_map, self%global2local, nelf, if_byte_sw)
               call nek_log_information('    Read vz2d ...', this_module, 'load_2d_fields')
               call load_and_distribute_slice(self%vz2d(:,:,:,i), global_map, self%global2local, nelf, if_byte_sw)
            end do ! 1, nsaver
            ! master closes the file
            if (nid == 0) then 
               call byte_close(ierr)
               if (ierr /= 0) call stop_error('Error closing file '//trim(fname), procedure='load_2d_fields')
            end if
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

      ! Helper functions

         subroutine gather_and_write_slice(slicedata, n2d_own)
            real(dp), intent(in) :: slicedata(:,:,:)
            integer, intent(in) :: n2d_own(:)
            ! internal
            integer :: nxy, idum, wdsl, isl, len, ierr, ip
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
                  call csend(ip,idum,isize,ip,0) ! hand shake
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
               if (ierr /= 0) call stop_error('Error writing slice data', procedure='gather_and_write_slice')
            else 
               ! send data to master
               call crecv2(nid,idum,isize,0) ! hand shake
               len = nxy*n2d_own(nid+1)
               call csend(nid,slicedata,len*wdsize,0,0)
            end if
         end subroutine gather_and_write_slice

         subroutine load_and_distribute_slice(slicedata, global_map, global2local, nelf, if_byte_sw)
            real(dp), intent(out) :: slicedata(:,:,:)
            integer, intent(in) :: global_map(:)
            integer, intent(in) :: global2local(:)
            integer, intent(in) :: nelf
            logical, intent(in) :: if_byte_sw
            ! internal
            real(dp) :: eldata(lx1,ly1)
            integer :: nxy, ierr, iel, iel_local, ieg, ie, ip, mtype, wdsl
            real dummy, rtmpv(lx1,ly1)
            character(len=128) :: msg
            real(dp) :: elavg
            nxy = lx1*ly1
            wdsl = wdsize/4
            do iel = 1, nelf
               if (nid == 0) then
                  call byte_read(eldata, nxy*wdsl, ierr)
                  if (if_byte_sw) call byte_reverse(eldata, nxy, ierr)
                  if (ierr /= 0) call stop_error('Error reading element data', procedure='load_and_distribute_slice')
               end if
               ieg = global_map(iel) ! get global element number
               ie  = gllel(ieg)      ! get local element number
               ip  = gllnid(ieg)     ! get processor owning element
               iel_local = global2local(ieg) ! local iel for global element
               mtype = 5000+ie
               ! distribute
               if (nid == 0 .and. ip /= 0) then ! send data from reader
                  call csend (mtype,eldata,isize,ip,nullpid)
                  call crecv2(mtype,dummy,isize,ip) ! hand shake, sync procs
                  call csend (mtype,eldata,wdsize*nxy,ip,nullpid)
               else if (nid /= 0 .and. ip == nid) then ! get data from reader
                  call crecv2(mtype,dummy,isize,0)
                  call csend (mtype,rtmpv,isize,0,nullpid) ! hand shake, sync procs
                  call crecv2(mtype,rtmpv,nxy*wdsize,0)
                  if (iel_local /= -1) then
                     call copy(slicedata(1,1,iel_local),rtmpv,nxy)
                  else
                     write(msg, '(A,I0,A,I0)') 'Inconsistent global2local mapping for el ', ieg, ' on nid ', ip
                     call stop_error(msg, this_module, 'load_and_distribute_slice')
                  end if
               end if
               if (nid == 0) then
                  elavg = sum(eldata)/nxy
                  write(msg,'(5(A,1x,I3),1X,A,1X,E15.7)') 'element ', iel, ', global:', ieg, ', local:', ie, ', iel2d:', 
     &                     iel_local,', own:', ip, ', elavg:', elavg
                  call nek_log_debug(msg, this_module, 'load_and_distribute_slice')
                  if (ip == 0) then ! reader owns element
                     if (iel_local /= -1) then
                        call copy(slicedata(1,1,iel_local),eldata,nxy)
                     else
                        write(msg, '(A,I0,A,I0)') 'Inconsistent global2local mapping for el ', ieg, ' on nid ', ip
                        call stop_error(msg, this_module, 'load_and_distribute_slice')
                     end if
                  end if
               end if
            end do
         end subroutine load_and_distribute_slice
      
      end module neklab_helix
