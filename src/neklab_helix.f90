      module neklab_helix
      !---------------------------------------
      !-----     LightKrylov Imports     -----
      !---------------------------------------
         use stdlib_strings, only: padl
         use stdlib_optval, only: optval
         use stdlib_sorting, only: sort_index
         use stdlib_logger, only: debug_level, information_level
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
            real(dp), dimension(nf) :: dpds
            real(dp), dimension(lx1,ly1,lz1,lelv) :: fshape
            ! mesh inputs
            integer :: nslices
            integer :: nelf
            ! sanity check
            logical :: is_initialized = .false.
            ! data
            ! Sweep angle (radians)
            real(dp), dimension(lx1,ly1,lz1,lelv) :: sweep_angle
            ! Angle of the cross-sectional plane around the helix, clockwise from the (positive) y axis
            real(dp), dimension(lx1,ly1,lz1,lelv) :: as
            ! Angle within the cross-sectional plane, from the inside to the outside of the helix starting from the negative z axis
            real(dp), dimension(lx1,ly1,lz1,lelv) :: alpha 
            ! cylindrical coordinates w.r.t the equatorial plane of the helix (with zax)
            real(dp), dimension(lx1,ly1,lz1,lelv) :: ox, oy
            ! cartesian coordinates in torus (without torsion!)
            real(dp), dimension(lx1,ly1,lz1,lelv) :: xax, yax, zax
            ! 2D data
            integer :: n2d      ! global number of 2d elements
            integer :: n2d_gown ! number of 2d elements globally owned by current processor
            integer :: n2d_lown ! number of 2d elements locally  owned by current processor
            integer :: nsave    ! buffer fill counter
            integer :: noutc    ! number of data files in cartesian coordinates written to disk
            integer :: noutt    ! number of data files in toroidal coordinates written to disk
            integer :: nload    ! number of loaded 2D baseflow fields
            logical :: save_2d_usrt = .true.! save us,ur,ut in addition to vx,vy,vz?
            ! save 2D fields
            logical, dimension(lelv)   :: lowner   ! is the local element the local segment owner?
            logical, dimension(lelv)   :: gowner   ! is the local element the global segmet owner? (first slice)
            integer, dimension(lelv)   :: lsegment ! pointer to the local  segment the element belongs to
            integer, dimension(lelv)   :: gsegment ! pointer to the global segment the element belongs to
            integer, dimension(lelv,3) :: id2d     ! characterisation f the locally owned segments
            real(dp), dimension(lbuf)              :: dt2d ! timestep information for the saved 2d snapshots
            real(dp), dimension(lx1,ly1,lelv)      :: x2d, y2d ! coordinates of the reference 2d slice
            real(dp), dimension(lx1,ly1,lelv,lbuf) :: vx2d, vy2d, vz2d ! 2D velocity fields
         contains
            ! initialization
            procedure, pass(self), public :: init_geom
            procedure, pass(self), public :: init_flow
            ! public computation routines
            procedure, pass(self), public :: compute_fshape
            procedure, pass(self), public :: compute_bf_forcing
            procedure, pass(self), public :: compute_usrt
            procedure, pass(self), public :: compute_ubar
            ! 2D data manipulation
            procedure, pass(self), public :: save_2d_fields
            procedure, pass(self), public :: compute_2dusrt
            procedure, pass(self), public :: outpost_2d_fields
            procedure, pass(self), public :: load_2d_fields
            procedure, pass(self), public :: set_baseflow
            ! helper routines
            procedure, pass(self), public :: get_forcing
            procedure, pass(self), public :: get_period
            procedure, pass(self), public :: add_dpds
            procedure, pass(self), public :: is_steady
            procedure, pass(self), public :: setup_summary
            procedure, pass(self), public :: get_fshape
            procedure, pass(self), public :: get_angle_s
            procedure, pass(self), public :: get_alpha
            procedure, pass(self), public :: is_lowner
            procedure, pass(self), public :: is_gowner
            procedure, pass(self), public :: get_lsegment
            procedure, pass(self), public :: get_gsegment
            procedure, pass(self), public :: get_v2d
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
            real(dp) :: x_torus, y_torus, z_torus, sweep, r
            real(dp), dimension(lx1,ly1,lz1,lelv) :: tmp, pipe_r
            integer :: ix, iy, iz, ie, iel, ieg, iseg, iface, isl, level, nxy
            integer :: fileid
            integer, dimension(lelv) :: islice
            integer, dimension(:), allocatable :: unique_segments, segment_owner, segment_count
            integer, dimension(:), allocatable :: idx ! for findloc
            logical, dimension(:), allocatable :: segment_found
            character(len=3) :: fid
            ! functions
            real(dp), external :: glmax, glmin
            integer, external :: iglsum

            if (self%is_initialized) call stop_error('Attempting to reinitialize the mesh', this_module, 'init_geom')

         !  Geometry modification for helical pipe

            pi = 4.0_dp*atan(1.0_dp)
            nxy = lx1*ly1

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
            call izero(islice, nelv)
            call izero(self%id2d, 3*nelv)
            self%n2d_lown = 0
            self%n2d_gown = 0
            do ie = 1, nelv
               ieg = lglel(ie)
               isl = ieg/self%nelf + 1
               if (mod(ieg,self%nelf)==0) isl = isl - 1
               islice(ie) = isl
               self%gsegment(ie) = ieg - (isl-1)*self%nelf
               if (isl == 1) then
                  ! the element in the first slice is the global segment owner
                  self%gowner(ie) = .true.
                  self%n2d_gown = self%n2d_gown + 1
               end if
               ! gather all unique segments on proc
               if (.not. segment_found(self%gsegment(ie))) then
                  self%n2d_lown = self%n2d_lown + 1
                  segment_owner(self%n2d_lown) = ie
                  unique_segments(self%n2d_lown) = self%gsegment(ie)
                  segment_found(self%gsegment(ie)) = .true.
                  self%id2d(self%n2d_lown,1) = ie
                  self%id2d(self%n2d_lown,3) = self%gsegment(ie)
                  self%lowner(ie) = .true.
               end if
               idx = findloc(unique_segments, self%gsegment(ie))
               segment_count(idx(1)) = segment_count(idx(1)) + 1
               self%lsegment(ie) = idx(1)
            end do
      
            call logger%configuration(level=level)
            if (level <= debug_level) then
               print '(A,2(I0,1X),A,*(1X,I3))', 'DEBUG 2dmap: ', nid, self%n2d_lown, 'unique streamwise segments:   ', unique_segments(:self%n2d_lown)
               print '(A,2(I0,1X),A,*(1X,I3))', 'DEBUG 2dmap: ', nid, self%n2d_lown, '# of elements in segment:     ', segment_count(:self%n2d_lown)
               print '(A,2(I0,1X),A,*(1X,I3))', 'DEBUG 2dmap: ', nid, self%n2d_lown, 'local element local  s. owner:', segment_owner(:self%n2d_lown)
               print '(A,2(I0,1X),A,*(1X,L3))', 'DEBUG 2dmap: ', nid, self%n2d_lown, 'local element global s. owner:', self%gowner(:self%n2d_lown)
               print '(A,2(I0,1X),A,I0)'      , 'DEBUG 2dmap: ', nid, self%n2d_lown, 'globally owned: ', self%n2d_gown
               call nekgsync()
               do ie = 1, nelv
                  call cfill(vz(1,1,1,ie), 1.0_dp*nid, lx1*ly1*lz1)
                  call cfill(vx(1,1,1,ie), 1.0_dp*islice(ie),lx1*ly1*lz1)
                  call cfill(vy(1,1,1,ie), 1.0_dp*self%gsegment(ie),lx1*ly1*lz1)
                  if (self%lowner(ie)) call cfill(vz(1,1,1,ie), -1.0_dp,lx1*ly1*lz1) 
                  if (self%gowner(ie)) call cfill(vz(1,1,1,ie), np*1.0_dp,lx1*ly1*lz1)
               end do
               call outpost(vx,vy,vz,pr,t,'m2d')
               print '(A,A3,3X,3(A5),2X,2(A20))', 'DEBUG: 2dmap ', 'nid', 'glob', 'locl', 'l2d', 'xavg', 'yavg'
            end if

            ! Extract 2D mesh
            iseg = 0 
		      do ie = 1, nelv
               if (self%lowner(ie)) then
		      	   ! extract boundary points from the global segment owners
                  iseg = iseg + 1
                  if (self%gowner(ie)) then
		      	      do iface = 1, 2*ndim
		      	      	if (cbc(iface,ie,1) == 'P') then
		      	   	   	call ftovec(self%x2d(1,1,iseg), zm1, ie, iface, nx1, ny1, nz1) ! z --> x
		      	   		   call ftovec(self%y2d(1,1,iseg), ym1, ie, iface, nx1, ny1, nz1)
                           self%id2d(iseg,2) = iface
                       end if
		      	      end do
                  end if
                  if (level <= debug_level) then
                     print '(A,I3,A,4(1X,I4),A,3X,F17.8,3x,F17.8)', 'DEBUG 2dmap: ', nid, ' el', lglel(ie), ie, iseg, self%gsegment(ie),  
     &                           ': ', sum(self%x2d(:,:,iseg))/nxy, sum(self%y2d(:,:,iseg))/nxy
                  end if
               end if
		      end do
            
            ! initialize data
            call rzero(self%vx2d, nx1*ny1*nelv*lbuf)
            call rzero(self%vy2d, nx1*ny1*nelv*lbuf)
            call rzero(self%vz2d, nx1*ny1*nelv*lbuf)
            call rzero(self%dt2d, lbuf)
            self%nsave = 0
            self%noutc = 0
            self%noutt = 0
            self%n2d   = iglsum(self%n2d_gown,1)
            self%nload = 0
            if (self%n2d /= self%nelf) call stop_error('Inconsistent elements in 2D mesh!', module=this_module, procedure='init_geom')
            if (level <= debug_level) then
               call nekgsync()
               print '(A,I3,A,*(1x,I0))', 'DEBUG 2dmap: ', nid, ', nelv2iseg: ', self%lsegment(:nelv)
               write(fid,'(I3.3)') nid
               fileid = 2000+nid
               open (fileid, file='torus_map'//fid//'.txt', status='replace', action='write')
               write(fileid, *) 'nelv = ', nelv
               write(fileid, '(6(1X,A11),2(1X,A7),A12)') 'ie','ieg','slice','self%gsegment','gllel','gllnid','s%lowner','s%gowner','s%n2iseg'
               do ie = 1, nelv
                  ieg = lglel(ie)
                  write(fileid, '(6(I12),2(1X,L7),I12)') ie, ieg, islice(ie), self%gsegment(ie), gllel(ieg), gllnid(ieg), self%lowner(ie), self%gowner(ie), self%lsegment(ie)
               end do
               write(fileid, *) 'n2d_lown = ', self%n2d_lown
               write(fileid, '(*(1X,A11))') 'iel','ieg','s%id2d:ie', 's%id2d:ifc', 's%id2d%iseg'
               do ie = 1, self%n2d_lown
                  write(fileid, *) ie, lglel(self%id2d(ie,1)), self%id2d(ie,:)
               end do
               close (fileid)
               call nekgsync()
            end if

            ! Morph the mesh into a torus
            helix_r = self%curv_radius
            do ie = 1, nelv
            do iz = 1, lz1
            do iy = 1, ly1
            do ix = 1, lx1
               r     = pipe_r(ix,iy,iz,ie)
               sweep = self%sweep_angle(ix,iy,iz,ie)
               self%ox(ix,iy,iz,ie) = helix_r * sin(sweep)
               self%oy(ix,iy,iz,ie) = helix_r * cos(sweep)
               xm1(ix,iy,iz,ie)     = r * sin(sweep) + self%ox(ix,iy,iz,ie)
               ym1(ix,iy,iz,ie)     = r * cos(sweep) + self%oy(ix,iy,iz,ie)
            end do
            end do
            end do
            end do
            call copy(self%xax, xm1, lv) ! xax set before curvature in z is added!
            call copy(self%yax, ym1, lv) ! yax set before curvature in z is added!

            ! Morph the torus into a helix
            if (self%phi /= 0.0_dp) then
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
            integer :: ix, iy, iz, ie
            real(dp) :: helix_r2, r, rr, alpha
            self%fshape = 0.0_dp
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
            integer :: ix, iy, iz, ie
            real(dp) :: fs, phi
            real(dp), dimension(lx1,ly1,lz1,lelv) :: ffx, ffy, ffz
            fs = self%get_forcing(t) / self%curv_radius

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
            integer :: ix, iy, iz, ie
            real(dp) :: phi, a, s, ux, uy, uz, utmp, vtmp
            phi = self%phi
            do ie = 1, nelv
            do iz = 1, lz1
            do iy = 1, ly1
            do ix = 1, lx1
               s  = self%as(ix,iy,iy,ie)
               a  = self%alpha(ix,iy,iy,ie)
               ux = u(ix,iy,iy,ie)
               uy = v(ix,iy,iy,ie)
               uz = w(ix,iy,iy,ie)
               utmp            = sin(s)*ux + cos(s)*uy
               vtmp            = sin(phi) * (-cos(s)*ux - sin(s)*uy) + cos(phi)*uz
               us(ix,iy,iz,ie) = cos(phi) * ( cos(s)*ux - sin(s)*uy) + sin(phi)*uz
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
            integer :: ix, iy, iz, ie
            real(dp) :: num, den, us, us_r, ux, uy, uz, phi, s, a, fs
            real(dp), external :: glsum
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
         end function compute_ubar

         subroutine save_2d_fields(self, u, v, w)
            class(helix), intent(inout) :: self
            real(dp), dimension(lx1,ly1,lz1,lelv), intent(in) :: u
            real(dp), dimension(lx1,ly1,lz1,lelv), intent(in) :: v
            real(dp), dimension(lx1,ly1,lz1,lelv), intent(in) :: w
            ! internal
            integer :: iseg, ie, ifc, level, nxy
            real(dp) :: xavg, yavg, vxavg, vyavg, vzavg
            character(len=128) :: msg
            nxy = lx1*ly1
            call logger%configuration(level=level)
            self%nsave = self%nsave + 1
            ! save data to buffer
            write(msg,'(A,I3,A,I3)') 'Save 2D data: ', self%nsave, '/', lbuf
            call nek_log_message(msg, this_module, 'save_2d_fields')
            do iseg = 1, self%n2d_gown
               ie  = self%id2d(iseg, 1)
               ifc = self%id2d(iseg, 2)
               call ftovec(self%vx2d(1,1,iseg,self%nsave), u, ie, ifc, nx1, ny1, nz1)
		      	call ftovec(self%vy2d(1,1,iseg,self%nsave), v, ie, ifc, nx1, ny1, nz1)
		      	call ftovec(self%vz2d(1,1,iseg,self%nsave), w, ie, ifc, nx1, ny1, nz1)
               if (level <= debug_level) then
                  xavg  = sum(self%x2d (:,:,iseg))/nxy
                  yavg  = sum(self%y2d (:,:,iseg))/nxy
                  vxavg = sum(self%vx2d(:,:,iseg,self%nsave))/nxy
                  vyavg = sum(self%vy2d(:,:,iseg,self%nsave))/nxy
                  vzavg = sum(self%vz2d(:,:,iseg,self%nsave))/nxy
                  print '(A,I8,I8,A,5(3X,F16.8))', 'DEBUG: save el', ie, iseg, ': ', xavg, yavg, vxavg, vyavg, vzavg
               end if
            end do
            self%dt2d(self%nsave) = dt
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
            integer :: ix, iy, ie, is, ib
            real(dp) :: phi, sweep, a, s
            real(dp) :: utmp, vtmp, ux, uy, uz
            phi = self%phi
            do is = 1, self%n2d_gown ! only for the first slice
               ie = self%id2d(is, 1)
               do iy = 1, ly1
               do ix = 1, lx1
                  s = self%as(ix,iy,iz,ie)
                  a = self%alpha(ix,iy,iz,ie)
                  ! iterate over buffer
                  do ib = 1, lbuf
                     ux = self%vx2d(ix,iy,is,ib)
                     uy = self%vy2d(ix,iy,is,ib)
                     uz = self%vz2d(ix,iy,is,ib)
                     ! overwrite v[xyz]2d with u[srt]2d
                     self%vx2d(ix,iy,is,ib) = cos(phi)*( cos(s)*ux -sin(s)*uy) + sin(phi)*uz
                     utmp                   = sin(s)*ux + cos(s)*uy
                     vtmp                   = sin(phi)*(-cos(s)*ux -sin(s)*uy) + cos(phi)*uz
                     self%vy2d(ix,iy,is,ib) = cos(a)*utmp + sin(a)*vtmp
                     self%vz2d(ix,iy,is,ib) = sin(a)*utmp - cos(a)*vtmp
                  end do ! lbuf
               end do    ! lx1
               end do    ! ly1
            end do       ! self%n2d_gown
         end subroutine compute_2dusrt
         
         subroutine outpost_2d_fields(self, iname)
            class(helix), intent(inout) :: self
            character(len=1), intent(in) :: iname
            ! internals
            integer, allocatable :: n2d_gown(:)
            integer, allocatable :: n2d_elmap(:)
            integer :: ierr, itmp, i, nxy, ip, ibuf, iseg, len, i_own
            integer :: wdsl, isl, isend(lelv)
            character(len=128)  :: fname, msg
            character(len=1024) :: head, ftm
            real rtmpv1(lx1*ly1*lelv), rtmpv(lx1*ly1*lelv)
            real*4 rtmpv2(2*lx1*ly1*lelv)
            equivalence (rtmpv1,rtmpv2)
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
            allocate(n2d_gown(np))           ! number of elements owned by each proc
            call izero(n2d_gown,np)
            allocate(n2d_elmap(self%nelf))  ! global element number of owned elements
            call izero(n2d_elmap,self%nelf)
            ! determine how many elements to dump
            if (nid == 0) then
               ! first for the master node 
               n2d_gown(1) = self%n2d_gown
               do i = 1, self%n2d_gown
                  n2d_elmap(i) = lglel(self%id2d(i,1))     ! get global element number
               end do 
               iseg = self%n2d_gown
               ! then gather info from other procs
               do ip = 1, np-1
                  call csend(ip,itmp,isize,ip,0)           ! hand shake
                  call crecv(ip,i_own,isize)               ! recv number of elements
                  n2d_gown(ip+1) = i_own
                  call crecv(ip,isend(:i_own),i_own*isize) ! recv global element map
                  n2d_elmap(iseg+1:iseg+i_own) = isend(:i_own)
                  iseg = iseg + i_own
               enddo
               if (iseg /= self%nelf) call stop_error('Not all elements in slice found!', this_module, 'outpost_2d_fields')
               ! write it to file
               call byte_write(n2d_elmap,self%nelf*isl,ierr)
               ! write timestep information to file
               call byte_write(self%dt2d(:self%nsave),self%nsave*wdsl,ierr)
            else
               call crecv(nid,itmp,isize)                  ! hand shake
               call csend(nid,self%n2d_gown,isize,0,0)     ! send number of elements
               len = self%n2d_gown
               do i = 1, self%n2d_gown
                  isend(i) = lglel(self%id2d(i,1))       
               end do
               call csend(nid,isend(:len),len*isize,0,0)   ! send global element map
            endif
            call bcast(n2d_gown, np*isize)         ! broadcast to all procs
            call bcast(n2d_elmap, self%nelf*isize) ! broadcast to all procs

            ! coordinates
            call nek_log_message('   '//trim(fname)//': write x2d ...', this_module, 'outpost_2d_fields')
            call gather_and_write_slice(self%x2d, n2d_gown)
            call nek_log_message('   '//trim(fname)//': write y2d ...', this_module, 'outpost_2d_fields')
            call gather_and_write_slice(self%y2d, n2d_gown)
            ! velocity data
            do i = 1, self%nsave
               write(msg,'(3X,A,A,1X,I3)') trim(fname),': write v[xyz]2d', i
               call nek_log_message(msg, this_module, 'outpost_2d_fields')
               call gather_and_write_slice(self%vx2d(:,:,:,i), n2d_gown)
               call gather_and_write_slice(self%vy2d(:,:,:,i), n2d_gown)
               call gather_and_write_slice(self%vz2d(:,:,:,i), n2d_gown)
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
            integer :: nxr, nyr, nelfr, nsaver, lbufr, wdsizr, len
            integer :: wdsl, isl, itmp, ip, nxy, i, ie, ieg, iel, iseg, gseg, nelf
            real rtmpv(lx1*ly1)
            integer, allocatable :: global_map(:)
            integer, allocatable :: gmap_index(:)
            real(dp), allocatable :: slicedata(:,:,:,:)
            real(dp) :: dt2dr(lbuf)
            real(dp) :: timer
            character(len=132) :: hdr, fname, msg
            character(len=4)   :: sdummy
            character(len=3)   :: fid
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
            allocate(gmap_index(nelf))
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
               ! read timestep information
               call byte_read(dt2dr(:nsaver), nsaver*wdsl, ierr)
               if (ierr /= 0) call stop_error('Error reading timestep information from file '//trim(fname), procedure='load_2d_fields')
               ! read coords but skip them
               call byte_read(fldum, nxy*nelf*wdsl, ierr)
               call byte_read(fldum, nxy*nelf*wdsl, ierr)
               if (ierr /= 0) call stop_error('Error reading coordinates from file '//trim(fname), procedure='load_2d_fields')
            end if
            call bcast(nsaver, isize)          ! broadcast number of saved snapshots
            call bcast(global_map, nelf*isize) ! broadcast global element map
            call sort_index(global_map, gmap_index)
            ! load data one timestep at a time
            allocate(slicedata(lx1,ly1,nelf,3))
            call rzero(self%vx2d, nxy*nelv)
            len = 3*nxy*nelf
            do i = 1, nsaver
               if (nid == 0) then ! read v[xyz]2d for all elements at the current timestep
                  call byte_read(slicedata, len*wdsl, ierr)
                  if (if_byte_sw) call byte_reverse(slicedata, len, ierr)
                  if (ierr /= 0) call stop_error('Error reading element data', procedure='load_and_distribute_slice')
               end if
               call bcast(slicedata, len*wdsize) ! broadcast 2D data to all procs
               ! distribute to local segment owners
               do iseg = 1, self%n2d_lown
                  gseg = self%id2d(iseg,3)  ! global segment
                  iel  = gmap_index(gseg)   ! get element that is read
                  call copy(self%vx2d(1,1,iseg,i), slicedata(1,1,iel,1), nxy)
                  call copy(self%vy2d(1,1,iseg,i), slicedata(1,1,iel,2), nxy)
                  call copy(self%vz2d(1,1,iseg,i), slicedata(1,1,iel,3), nxy)
               end do
            end do ! 1, nsaver
            ! master closes the file
            if (nid == 0) then 
               call byte_close(ierr)
               if (ierr /= 0) call stop_error('Error closing file '//trim(fname), procedure='load_2d_fields')
            end if
            self%nload = nsaver
         end subroutine load_2d_fields

         subroutine set_baseflow(self, basex, basey, basez, ifld)
            class(helix), intent(in) :: self
            real(dp), dimension(lx1,ly1,lz1,lelv), intent(out) :: basex
            real(dp), dimension(lx1,ly1,lz1,lelv), intent(out) :: basey
            real(dp), dimension(lx1,ly1,lz1,lelv), intent(out) :: basez
            integer, intent(in) :: ifld
            ! internal
            integer  :: ie, ix, iy, iz, iseg
            real(dp) :: a, s, phi, u, v, w
            phi = self%phi
            if (ifld > self%nload) call stop_error('Inconsistent ifld!', this_module, 'set_baseflow')
            do ie = 1, nelv
            iseg = self%lsegment(ie) ! local segment
            do iz = 1, lz1
            do ix = 1, ly1
            do iy = 1, lx1
               u = self%vx2d(ix,iy,iseg,ifld)
               v = self%vy2d(ix,iy,iseg,ifld)
               w = self%vz2d(ix,iy,iseg,ifld)
               s = self%as(ix,iy,iz,ie)
               a = self%alpha(ix,iy,iz,ie)
               basex(ix,iy,iz,ie) = cos(phi)*( cos(s)*u - sin(s)*v) + sin(phi)*w
               basey(ix,iy,iz,ie) =            sin(s)*u + cos(s)*v
               basez(ix,iy,iz,ie) = sin(phi)*(-cos(s)*u - sin(s)*v) + cos(phi)*w
            end do
            end do
            end do
            end do
         end subroutine set_baseflow

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

         logical pure function is_lowner(self, ie) result(is_owner)
            class(helix), intent(in) :: self
            integer, intent(in) :: ie
            is_owner = .false.
            if (ie <= nelv) is_owner = self%lowner(ie)
         end function is_lowner

         logical pure function is_gowner(self, ie) result(is_owner)
            class(helix), intent(in) :: self
            integer, intent(in) :: ie
            is_owner = .false.
            if (ie <= nelv) is_owner = self%gowner(ie)
         end function is_gowner
         
         integer pure function get_lsegment(self, ie) result(local_segment)
            class(helix), intent(in) :: self
            integer, intent(in) :: ie
            local_segment = 0
            if (ie <= nelv) local_segment = self%lsegment(ie)
         end function get_lsegment
      
         integer pure function get_gsegment(self, ie) result(global_segment)
            class(helix), intent(in) :: self
            integer, intent(in) :: ie
            global_segment = 0
            if (ie <= nelv) global_segment = self%gsegment(ie)
         end function get_gsegment

         real(dp) pure function get_v2d(self,ix,iy,iseg,ifld,icomp) result(v2d)
            class(helix), intent(in) :: self
            integer, intent(in) :: ix
            integer, intent(in) :: iy
            integer, intent(in) :: iseg
            integer, intent(in) :: ifld
            integer, intent(in) :: icomp
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
         end function get_v2d

      ! Helper functions

         subroutine gather_and_write_slice(slicedata, n2d_gown)
            real(dp), intent(in) :: slicedata(:,:,:)
            integer, intent(in) :: n2d_gown(:)
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
               len = nxy*n2d_gown(nid+1)
               if (wdsl.eq.2) then
                  call copy(rtmpv1,slicedata,len)
                  call byte_write(rtmpv2,len*wdsl,ierr)
               else
                  call copyX4(rtmpv2,slicedata,len)
                  call byte_write(rtmpv2,len,ierr)
               end if
               ! get data from other procs and write to file
               do ip = 1, np-1
                  len = nxy*n2d_gown(ip+1)
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
               len = nxy*n2d_gown(nid+1)
               call csend(nid,slicedata,len*wdsize,0,0)
            end if
         end subroutine gather_and_write_slice
      
      end module neklab_helix