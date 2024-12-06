      module neklab_helix
      !---------------------------------------
      !-----     LightKrylov Imports     -----
      !---------------------------------------
         use stdlib_sorting, only: sort, sort_index
         use stdlib_linalg, only: eig
         use stdlib_optval, only: optval
      ! Default real kind.
         use LightKrylov, only: dp
         use LightKrylov_Constants, only: imag => one_im_cdp
      ! Abstract types for real-valued linear operators and vectors.
         use LightKrylov, only: abstract_linop_rdp, abstract_vector_rdp
         use LightKrylov, only: orthonormalize_basis, zero_basis, rand_basis
         use LightKrylov_Utils, only: abstract_opts
      ! Logging
         use LightKrylov_Logger
      ! Extensions of the abstract vector types to nek data format.
         use neklab_vectors
         use neklab_vectors, only: nf => n_forcing
         use neklab_nek_forcing, only: neklab_forcing, set_neklab_forcing
         use neklab_utils, only: nek2vec, vec2nek
         use neklab_nek_setup, only: setup_linear_solver
         use neklab_linops, only: apply_L
         implicit none
         include "SIZE"
         include "TOTAL"
         include "ADJOINT"
         private
         character(len=*), parameter, private :: this_module = 'neklab_helix'
      
         integer, parameter :: lv = lx1*ly1*lz1*lelv
      !! Local number of grid points for the velocity mesh.
         integer, parameter :: lp = lx2*ly2*lz2*lelv
      !! Local number of grid points for the pressure mesh.

         public :: pipe
         public :: helix_pipe ! constructor for the pipe instance of the helix type
         public :: dpds_from_vector
         public :: dpds_to_vector

         type, public :: helix
            !! Type containing the basic geometrical and dynamical properties of a helix
            private
            ! geometry inputs
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
            real(dp), dimension(lv):: fshape
            real(dp), dimension(nf), public :: dpds
            ! mesh inputs
            real(dp) :: length
            real(dp) :: nslices
            real(dp) :: nelf
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
            ! dF/df
            type(nek_ext_dvector_forcing), dimension(nf), public :: dFdf
            ! dG/df
            real(dp), dimension(nf), public :: dGdf
            ! flag to trigger computation
            logical, public :: to_compute_df = .false.
         contains
            procedure, pass(self), public :: init_geom
            procedure, pass(self), public :: init_flow
            procedure, pass(self), public :: compute_fshape
            procedure, pass(self), public :: compute_bf_forcing
            procedure, pass(self), public :: compute_usrt
            procedure, pass(self), public :: compute_ubar
            procedure, pass(self), public :: get_forcing
            procedure, pass(self), public :: get_period
            procedure, pass(self), public :: get_dpds
            procedure, pass(self), public :: set_dpds
            procedure, pass(self), public :: add_dpds_perturbation
            procedure, pass(self), public :: is_steady
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

            if (nid.eq.0) write(6,*) 'USERDAT2: phi         = ', pipe%phi*180./pi
            if (nid.eq.0) write(6,*) 'USERDAT2: curv_radius = ', pipe%curv_radius

         end subroutine helix_pipe

         subroutine init_geom(self)
            class(helix), intent(inout) :: self
            ! internal
            real(dp) :: xmin, xmax, helix_r, s_angle
            real(dp) :: x_torus, y_torus, z_torus
            real(dp), dimension(lv) :: tmp, pipe_r
            integer :: ix, iy, iz, ie, i
            real(dp), external :: glmax, glmin

            if (self%is_initialized) call stop_error('Attempting to reinitialize the mesh', this_module, 'helix%init_mesh')

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

            ! Morph the mesh into a torus
            i = 0
            helix_r = self%curv_radius
            do ie = 1, nelv
            do iz = 1, lz1
            do iy = 1, ly1
            do ix = 1, lx1
               i = i + 1
               self%ox(i)          = helix_r   * sin(self%sweep_angle(i))
               self%oy(i)          = helix_r   * cos(self%sweep_angle(i))
               xm1(ix,iy,iz,ie)    = pipe_r(i) * sin(self%sweep_angle(i)) + self%ox(i)
               ym1(ix,iy,iz,ie)    = pipe_r(i) * cos(self%sweep_angle(i)) + self%oy(i)
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
            integer :: n           
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
         end subroutine init_flow

         subroutine compute_fshape(self)
            class(helix), intent(inout) :: self
            ! internals
            integer :: i, ix, iy, iz, ie
            real(dp) :: helix_r2, r, rr, alpha

            self%fshape = 0.0_dp

            i = 0
            print *, self%delta, self%curv_radius
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
            real(dp) :: x, y, z, alpha
            integer :: ix, iy, iz, ie, i
            real(dp) :: utmp, vtmp
            do ie = 1, nelv
            do iz = 1, lz1
            do iy = 1, ly1
            do ix = 1, lx1
               i = i+1
               x = self%xax(i)
               y = self%yax(i)
               z = self%zax(i)
               ! Compute fshape
               us(ix,iy,iz,ie) = cos( self%phi ) * (
     $                          + cos( self%as(i)   * u(ix,iy,iz,ie) )
     $                          - sin( self%as(i)   * v(ix,iy,iz,ie) ))
     $                          + sin( self%phi )   * w(ix,iy,iz,ie)
               utmp            = sin( self%as(i) )  * u(ix,iy,iz,ie)
     $                          + cos( self%as(i) ) * v(ix,iy,iz,ie)
               vtmp            = sin( self%phi ) * (
     $                          - cos( self%as(i) ) * u(ix,iy,iz,ie)
     $                          - sin( self%as(i) ) * v(ix,iy,iz,ie) )
     $                          + cos( self%phi ) * w(ix,iy,iz,ie)
               ! angle within cross-sectional plane, zero at the outside edge of the helix
               alpha = atan2(z,sqrt(x**2 + y**2) - self%curv_radius)
               ur(ix,iy,iz,ie) = cos( self%alpha(i) ) * utmp + sin( self%alpha(i) ) * vtmp
               ut(ix,iy,iz,ie) = sin( self%alpha(i) ) * utmp - cos( self%alpha(i) ) * vtmp
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
               ! Compute fshape
               us = cos(self%phi)*(
     $               cos( self%as(i) * u(ix,iy,iz,ie) )
     $             - sin( self%as(i) * v(ix,iy,iz,ie) ))
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
            !s_L = self%sweep*sqrt(self%curv_radius**2 + self%pitch_s**2)
            ubar = num/den  ! "1/r"-weighted volumetric average of streamwise velocity
            
         end function compute_ubar

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

         subroutine get_dpds(self, dpds)
            class(helix), intent(in) :: self
            real(dp), dimension(nf), intent(out) :: dpds
            dpds = self%dpds
         end subroutine get_dpds

         logical pure function is_steady(self) result(steady)
            class(helix), intent(in) :: self
            steady = self%if_steady
         end function is_steady

         subroutine add_dpds_perturbation(self, df, i)
            class(helix), intent(inout) :: self
            real(dp), intent(in) :: df
            integer, intent(in) :: i
            ! internal
            self%dpds(i) = self%dpds(i) + df
         end subroutine add_dpds_perturbation

         subroutine dpds_from_vector(dpds, vec)
            real(dp), dimension(nf), intent(out) :: dpds
            type(nek_ext_dvector_forcing), intent(in) :: vec
            dpds = vec%f
         end subroutine dpds_from_vector

         subroutine dpds_to_vector(dpds, vec)
            real(dp), dimension(nf), intent(in) :: dpds
            type(nek_ext_dvector_forcing), intent(out) :: vec
            vec%f = dpds
         end subroutine dpds_to_vector
      
      end module neklab_helix
