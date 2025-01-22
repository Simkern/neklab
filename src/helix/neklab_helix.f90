      module neklab_helix
      !---------------------------------------
      !-----     LightKrylov Imports     -----
      !---------------------------------------
         use stdlib_strings, only: padl
         use stdlib_optval, only: optval
         use stdlib_sorting, only: sort_index
         use stdlib_logger, only: all_level, debug_level, information_level
      ! Default real kind.
         use LightKrylov, only: dp
         use LightKrylov_Constants, only: imag => one_im_cdp
      ! Logging & timing
         use LightKrylov_Logger
         use LightKrylov_Timing, only: lk_timer => global_lightkrylov_timer
      ! Extensions of the abstract vector types to nek data format.
         use neklab_vectors
         use neklab_nek_forcing, only: neklab_forcing, set_neklab_forcing
         use neklab_nek_setup, only: nek_log_message, nek_log_information, nek_log_warning, nek_log_debug, nek_stop_error

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
      !! Maximum number of forcing components 1 x steady + 2 x (# unsteady)
         integer, parameter :: lbuf = 1000
      !! Maximum number of 2d fields to save before outposting
         integer, parameter, public :: nfft = 16
      !! Number of FT components to compute

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
            logical :: if_sym    ! is the mesh symmetric (only half the pipe)
            logical :: if_torus  ! is the mesh curved (toroidal)
            logical :: if_helix  ! is the mesh helical?
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
            integer :: n2d       ! global number of 2d elements
            integer :: n2d_gown  ! number of 2d elements globally owned by current processor
            integer :: n2d_lown  ! number of 2d elements locally  owned by current processor
            integer :: nsave = 0 ! buffer fill counter
            integer :: noutc = 0 ! number of data files in cartesian coordinates written to disk
            integer :: noutt = 0 ! number of data files in toroidal coordinates written to disk
            integer :: noutn = 0 ! number of data files for newton
            integer :: nload = 0 ! number of loaded 2D baseflow fields
            integer :: nsteps = 0 ! number of steps in period for nonlinear simulation with variable dt
            logical :: save_2d_usrt = .true.! save us,ur,ut in addition to vx,vy,vz?
            logical :: save_2d_base = .false.
            logical :: if_newton    = .false. ! are we in newton mode?
            logical :: if_floquet   = .false. ! are we in floquet mode?
            logical :: if_fft = .false. ! compute the FT of the streamwise mass flow us on the fly
            logical :: fft_is_extracted = .false. ! has the FT data been collected?
            real(dp), dimension(2*nfft + 1) :: fftv ! temporary array for mass flow FT computation
            real(dp), dimension(2*nfft + 1) :: mflow ! FT of the streamwise mass flow 
            real(dp), dimension(nfft + 1) :: mflow_amplitude ! FT amplitude of the streamwise mass flow 
            real(dp), dimension(nfft + 1) :: mflow_phase     ! FT phase of the streamwise mass flow 
            real(dp) :: fft_time  = 0.0_dp ! Current integration time
            real(dp) :: fft_rtime = 0.0_dp ! Integration time of the FT record
            real(dp) :: min_dt = 100.0_dp ! minimum dt over a period using variable dt to determine smallest dt posssible for constant dt nonlinear run
            real(dp) :: max_dt = 0.0_dp   ! maximum dt over a period using variable dt
            real(dp) :: ubar_lag = 0.0_dp ! lagged ubar for variable dt mass flow integration
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
            ! helix_utils
            procedure, pass(self), public :: init_geom
            procedure, pass(self), public :: init_flow
            procedure, pass(self), public :: compute_fshape
            procedure, pass(self), public :: compute_bf_forcing
            procedure, pass(self), public :: compute_usrt
            procedure, pass(self), public :: compute_ubar
            procedure, pass(self), public :: forcing_amplitude
            procedure, pass(self), public :: shift_mflow_phase
            procedure, pass(self), public :: setup_summary
            procedure, pass(self), public :: parameter_summary
            procedure, pass(self), public :: forcing_summary
            ! helix_2d
            procedure, pass(self) :: init_2d_geom
            procedure, pass(self), public :: save_2d_fields
            procedure, pass(self), public :: outpost_2d
            procedure, pass(self), public :: outpost_2d_fields
            procedure, pass(self), public :: load_2d_fields
            procedure, pass(self) :: get_nsteps_from_header
            procedure, pass(self), public :: set_baseflow
            procedure, pass(self), public :: compute_2d_usrt
            procedure, pass(self), public :: set_2d_mode
            ! helix_mflow_fft
            procedure, pass(self), public :: reset_mflow_fft
            procedure, pass(self), public :: compute_mflow_fft
            procedure, pass(self), public :: extract_mflow_fft
            procedure, pass(self), public :: get_mflow_fft
            ! helix_getters_setters
            procedure, pass(self), public :: set_save_base
            procedure, pass(self), public :: set_save_fft
            procedure, pass(self), public :: set_newton
            procedure, pass(self), public :: set_floquet
            procedure, pass(self), public :: set_symmetry
            procedure, pass(self), public :: is_steady
            procedure, pass(self), public :: is_newton
            procedure, pass(self), public :: is_floquet
            procedure, pass(self), public :: is_save_2d
            procedure, pass(self), public :: is_save_fft
            procedure, pass(self), public :: is_extracted_fft
            procedure, pass(self), public :: is_sym
            procedure, pass(self), public :: is_torus
            procedure, pass(self), public :: is_helix
            procedure, pass(self), public :: is_lowner
            procedure, pass(self), public :: is_gowner
            procedure, pass(self), public :: get_period
            procedure, pass(self), public :: get_omega
            procedure, pass(self), public :: get_Wo
            procedure, pass(self), public :: get_dpds
            procedure, pass(self), public :: get_nf
            procedure, pass(self), public :: get_fshape
            procedure, pass(self), public :: get_angle_s
            procedure, pass(self), public :: get_length
            procedure, pass(self), public :: get_delta
            procedure, pass(self), public :: get_diameter
            procedure, pass(self), public :: get_pitch_s
            procedure, pass(self), public :: get_radius
            procedure, pass(self), public :: get_curv_radius
            procedure, pass(self), public :: get_phi
            procedure, pass(self), public :: get_sweep
            procedure, pass(self), public :: get_nsteps
            procedure, pass(self), public :: get_dt_minmax
            procedure, pass(self), public :: get_ubar_lag
            procedure, pass(self), public :: get_lsegment
            procedure, pass(self), public :: get_gsegment
            procedure, pass(self), public :: get_v2d
            procedure, pass(self), public :: set_dpds
            procedure, pass(self), public :: set_nsteps
         end type helix
         
         interface
         
            !-----------------------------------------------------
            ! neklab_helix % helix_utils
            !
            ! Type-bound procedures

            module subroutine init_geom(self)
               class(helix), intent(inout) :: self
            end subroutine init_geom

            module subroutine init_flow(self, dpds, womersley)
               class(helix), intent(inout) :: self
               real(dp), intent(in) :: dpds(:)
               real(dp), optional, intent(in) :: womersley
            end subroutine init_flow

            module subroutine reset_newton(self)
               class(helix), intent(inout) :: self
            end subroutine

            module subroutine compute_fshape(self)
               class(helix), intent(inout) :: self
            end subroutine compute_fshape

            module subroutine compute_bf_forcing(self, t)
               class(helix), intent(in) :: self
               real(dp) :: t
               !! time
            end subroutine compute_bf_forcing

            module subroutine compute_usrt(self, u, v, w, us, ur, ut)
               class(helix), intent(in) :: self
               real(dp), dimension(lx1,ly1,lz1,lelv), intent(in) :: u
               real(dp), dimension(lx1,ly1,lz1,lelv), intent(in) :: v
               real(dp), dimension(lx1,ly1,lz1,lelv), intent(in) :: w
               real(dp), dimension(lx1,ly1,lz1,lelv), intent(out) :: us
               real(dp), dimension(lx1,ly1,lz1,lelv), intent(out) :: ur
               real(dp), dimension(lx1,ly1,lz1,lelv), intent(out) :: ut
            end subroutine compute_usrt
            
            module function compute_ubar(self,u,v,w) result(ubar)
               class(helix), intent(in) :: self
               real(dp), dimension(lx1,ly1,lz1,lelv), intent(in) :: u
               real(dp), dimension(lx1,ly1,lz1,lelv), intent(in) :: v
               real(dp), dimension(lx1,ly1,lz1,lelv), intent(in) :: w
               real(dp) :: ubar
            end function compute_ubar

            module pure function forcing_amplitude(self, t) result(f)
               class(helix), intent(in) :: self
               real(dp), intent(in) :: t
               real(dp) :: f
               !! time
            end function forcing_amplitude

            module subroutine shift_mflow_phase(self, icomp, target_mflow_phase)
               class(helix), intent(inout) :: self
               integer, intent(in) :: icomp
               real(dp), intent(in) :: target_mflow_phase
            end subroutine shift_mflow_phase

            module subroutine setup_summary(self)
               class(helix), intent(in) :: self
            end subroutine setup_summary
            
            module subroutine parameter_summary(self)
               class(helix), intent(in) :: self
            end subroutine parameter_summary

            module subroutine forcing_summary(self)
               class(helix), intent(in) :: self
            end subroutine forcing_summary

            !-----------------------------------------------------
            ! neklab_helix % helix_2d
            !
            ! Type-bound procedures           

            module subroutine init_2d_geom(self)
               class(helix), intent(inout) :: self
            end subroutine init_2d_geom

            module subroutine save_2d_fields(self, u, v, w)
               class(helix), intent(inout) :: self
               real(dp), dimension(lx1,ly1,lz1,lelv), intent(in) :: u
               real(dp), dimension(lx1,ly1,lz1,lelv), intent(in) :: v
               real(dp), dimension(lx1,ly1,lz1,lelv), intent(in) :: w
            end subroutine save_2d_fields
         
            module subroutine outpost_2d(self)
               class(helix), intent(inout) :: self
            end subroutine outpost_2d

            module subroutine outpost_2d_fields(self, iname, iout)
               class(helix), intent(inout) :: self
               character(len=1), intent(in) :: iname
               integer, intent(in) :: iout
            end subroutine outpost_2d_fields            

            module subroutine load_2d_fields(self, idx)
               ! only nid 0 will read
               class(helix), intent(inout) :: self
               integer, intent(in) :: idx
            end subroutine load_2d_fields

            module subroutine get_nsteps_from_header(self, fname, nsaver)
               ! only nid 0 will read
               class(helix), intent(in) :: self
               character(len=132), intent(in) :: fname
               integer, intent(out) :: nsaver
            end subroutine get_nsteps_from_header

            module subroutine set_baseflow(self, basex, basey, basez, ifld)
               class(helix), intent(inout) :: self
               real(dp), dimension(lx1,ly1,lz1,lelv), intent(out) :: basex
               real(dp), dimension(lx1,ly1,lz1,lelv), intent(out) :: basey
               real(dp), dimension(lx1,ly1,lz1,lelv), intent(out) :: basez
               integer, intent(in) :: ifld
            end subroutine set_baseflow

            module subroutine compute_2d_usrt(self)
               ! this routine will overwrite self%v[xyz]2d
               class(helix), intent(inout) :: self
            end subroutine compute_2d_usrt

            module subroutine set_2d_mode(self, mode)
               class(helix), intent(inout) :: self
               character(len=*), intent(in) :: mode
            end subroutine set_2d_mode

            !-----------------------------------------------------
            ! neklab_helix % helix_mflow_fft
            !
            ! Type-bound procedures 

            module subroutine reset_mflow_fft(self)
               class(helix), intent(inout) :: self
            end subroutine reset_mflow_fft

            module subroutine compute_mflow_fft(self, period, var_dt)
               ! only for constant dt
               class(helix), intent(inout) :: self
               real(dp), optional, intent(in) :: period
               logical, optional, intent(in) :: var_dt
            end subroutine compute_mflow_fft
            
            module subroutine extract_mflow_fft(self, period, if_amplitude)
               ! only for constant dt
               class(helix), intent(inout) :: self
               real(dp), optional, intent(in) :: period
               logical, optional, intent(in) :: if_amplitude
            end subroutine extract_mflow_fft

            module subroutine get_mflow_fft(self, mflow, phase, if_amplitude)
               class(helix), intent(in) :: self
               real(dp), allocatable, intent(out) :: mflow(:)
               real(dp), optional, allocatable, intent(out) :: phase(:)
               logical, optional, intent(in) :: if_amplitude
            end subroutine get_mflow_fft

            !-----------------------------------------------------
            ! neklab_helix % helix_gs
            !
            ! Type-bound procedures

            ! mode switches

            module subroutine set_save_base(self, if_save)
               class(helix), intent(inout) :: self
               logical, intent(in) :: if_save
            end subroutine set_save_base

            module subroutine set_save_fft(self, if_save_fft)
               class(helix), intent(inout) :: self
               logical, intent(in) :: if_save_fft
            end subroutine set_save_fft

            module subroutine set_newton(self, if_newton)
               class(helix), intent(inout) :: self
               logical, intent(in) :: if_newton
            end subroutine set_newton

            module subroutine set_floquet(self, if_floquet)
               class(helix), intent(inout) :: self
               logical, intent(in) :: if_floquet
            end subroutine set_floquet

            module subroutine set_symmetry(self, if_sym)
               class(helix), intent(inout) :: self
               logical, intent(in) :: if_sym
            end subroutine set_symmetry
            
            ! logicals

            module pure function is_steady(self) result(if_steady)
               class(helix), intent(in) :: self
               logical :: if_steady
            end function is_steady

            module pure function is_newton(self) result(if_newton)
               class(helix), intent(in) :: self
               logical :: if_newton
            end function is_newton

            module pure function is_floquet(self) result(if_floquet)
               class(helix), intent(in) :: self
               logical :: if_floquet
            end function is_floquet

            module pure function is_save_2d(self) result(if_save_2d_base)
               class(helix), intent(in) :: self
               logical :: if_save_2D_base
            end function is_save_2d

            module pure function is_save_fft(self) result(if_save_fft)
               class(helix), intent(in) :: self
               logical :: if_save_fft
            end function is_save_fft

            module pure function is_extracted_fft(self) result(is_extracted)
               class(helix), intent(in) :: self
               logical :: is_extracted
            end function is_extracted_fft

            module pure function is_sym(self) result(mesh_is_sym)
               class(helix), intent(in) :: self
               logical :: mesh_is_sym
            end function is_sym

            module pure function is_torus(self) result(mesh_is_torus)
               class(helix), intent(in) :: self
               logical :: mesh_is_torus
            end function is_torus

            module pure function is_helix(self) result(mesh_is_helix)
               class(helix), intent(in) :: self
               logical :: mesh_is_helix
            end function is_helix

            module pure function is_lowner(self, ie) result(is_owner)
               class(helix), intent(in) :: self
               integer, intent(in) :: ie
               logical :: is_owner
            end function is_lowner
            
            module pure function is_gowner(self, ie) result(is_owner)
               class(helix), intent(in) :: self
               integer, intent(in) :: ie
               logical :: is_owner
            end function is_gowner

            ! getters

            module pure function get_period(self) result(T)
               class(helix), intent(in) :: self
               real(dp) :: T
            end function get_period

            module pure function get_omega(self) result(omega)
               class(helix), intent(in) :: self
               real(dp) :: omega
            end function get_omega

            module pure function get_Wo(self) result(Wo)
               class(helix), intent(in) :: self
               real(dp) :: Wo
            end function get_Wo
            
            module subroutine get_dpds(self, dpds, phase)
               class(helix), intent(in) :: self
               real(dp), dimension(nf), intent(out) :: dpds
               real(dp), optional, allocatable, intent(out) :: phase(:)
            end subroutine get_dpds

            module pure function get_nf(self) result(n)
               class(helix), intent(in) :: self
               integer :: n
            end function get_nf

            module subroutine get_fshape(self, fshape)
               class(helix), intent(in) :: self
               real(dp), dimension(lx1,ly1,lz1,lelv), intent(out) :: fshape
            end subroutine get_fshape

            module subroutine get_angle_s(self, angle_s)
               class(helix), intent(in) :: self
               real(dp), dimension(lx1,ly1,lz1,lelv), intent(out) :: angle_s
            end subroutine get_angle_s

            module subroutine get_alpha(self, alpha)
               class(helix), intent(in) :: self
               real(dp), dimension(lx1,ly1,lz1,lelv), intent(out) :: alpha
            end subroutine get_alpha
            
            module pure function get_length(self) result(length)
               class(helix), intent(in) :: self
               real(dp) :: length
            end function get_length

            module pure function get_delta(self) result(delta)
               class(helix), intent(in) :: self
               real(dp) :: delta
            end function get_delta

            module pure function get_diameter(self) result(diameter)
               class(helix), intent(in) :: self
               real(dp) :: diameter
            end function get_diameter

            module pure function get_pitch_s(self) result(pitch_s)
               class(helix), intent(in) :: self
               real(dp) :: pitch_s
            end function get_pitch_s

            module pure function get_radius(self) result(radius)
               class(helix), intent(in) :: self
               real(dp) :: radius
            end function get_radius

            module pure function get_curv_radius(self) result(curv_radius)
               class(helix), intent(in) :: self
               real(dp) :: curv_radius
            end function get_curv_radius

            module pure function get_phi(self) result(phi)
               class(helix), intent(in) :: self
               real(dp) :: phi
            end function get_phi

            module pure function get_sweep(self) result(sweep)
               class(helix), intent(in) :: self
               real(dp) :: sweep
            end function get_sweep

            module function get_nsteps(self) result(ns)
               class(helix), intent(in) :: self
               integer :: ns
            end function get_nsteps

            module subroutine get_dt_minmax(self, dt_minmax)
               class(helix), intent(in) :: self
               real(dp), dimension(2), intent(out) :: dt_minmax
            end subroutine get_dt_minmax

            module function get_ubar_lag(self) result(ubar_lag)
               class(helix), intent(in) :: self
               real(dp) :: ubar_lag
            end function get_ubar_lag

            module pure function get_lsegment(self, ie) result(local_segment)
               class(helix), intent(in) :: self
               integer, intent(in) :: ie
               integer :: local_segment
            end function get_lsegment
         
            module pure function get_gsegment(self, ie) result(global_segment)
               class(helix), intent(in) :: self
               integer, intent(in) :: ie
               integer :: global_segment
            end function get_gsegment

            module pure function get_v2d(self,ix,iy,iseg,ifld,icomp) result(v2d)
               class(helix), intent(in) :: self
               integer, intent(in) :: ix
               integer, intent(in) :: iy
               integer, intent(in) :: iseg
               integer, intent(in) :: ifld
               integer, intent(in) :: icomp
               real(dp) :: v2d
            end function get_v2d

            ! setters

            module subroutine set_dpds(self, dpds, reset)
               class(helix), intent(inout) :: self
               real(dp), dimension(nf), intent(in) :: dpds
               logical, optional, intent(in) :: reset
            end subroutine set_dpds

            module subroutine set_nsteps(self, ns)
               class(helix), intent(inout) :: self
               integer, intent(in) :: ns
            end subroutine set_nsteps
            
         end interface

         type(helix) :: pipe

      contains

         ! Constructor for the module level instance of helix
         subroutine helix_pipe(delta, diameter, pitch_s, length, nslices, nelf, if_sym, if_debug)
            real(dp), intent(in) :: delta
            real(dp), intent(in) :: diameter
            real(dp), intent(in) :: pitch_s
            real(dp), intent(in) :: length
            integer, intent(in) :: nslices
            integer, intent(in) :: nelf
            logical, optional, intent(in) :: if_sym
            logical, optional, intent(in) :: if_debug
            ! internal
            logical :: debug
            character(len=128) :: msg
            debug = optval(if_debug, .false.)

            ! Geometry
            pipe%delta    = delta
            pipe%diameter = diameter
            pipe%pitch_s  = pitch_s
            pipe%length   = length

            ! Mesh specifics
            pipe%nslices  = nslices
            pipe%nelf     = nelf
            call pipe%set_symmetry(optval(if_sym, .false.))
            
            !  Derived quantities
            pipe%radius      = pipe%diameter*0.5_dp
            if (pipe%delta /= 0.0_dp) then
               pipe%curv_radius = 1.0_dp/pipe%delta
               pipe%if_torus = .true.
            else
               pipe%curv_radius = 0.0_dp
               pipe%if_torus = .false.
            end if
            pipe%phi         = atan2(pipe%pitch_s,pipe%curv_radius)
            if (pipe%delta /= 0.0_dp) then
               pipe%sweep    = pipe%length*cos(pipe%phi)/pipe%curv_radius ! sweep angle in radians
               pipe%if_helix = .true.
            else
               pipe%sweep    = 0.0_dp
               pipe%if_helix = .false.
            end if

            ! add timers
            call lk_timer%initialize() ! in case it has not been done
            call lk_timer%add_timer('neklab_helix_init_geom', start=.false.)
            call lk_timer%add_timer('neklab_helix_save_2d', start=.false.)
            call lk_timer%add_timer('neklab_helix_load_2d', start=.false.)
            call lk_timer%add_timer('neklab_helix_outpost_2d', start=.false.)
            call lk_timer%add_timer('neklab_helix_set_baseflow', start=.false.)
            call lk_timer%add_timer('neklab_helix_compute_ubar', start=.false.)
            call lk_timer%add_timer('neklab_helix_compute_mflow_fft', start=.false.)

            ! intialize geometry
            call pipe%init_geom()
            ! compute forcing distribution
            call pipe%compute_fshape()

            ! switch on FT in the unsteady case
            if (nf > 1) then
               call pipe%set_save_fft(.true.)
            end if

            if (debug) then
               call outpost(pipe%xax, pipe%yax, pipe%zax, pr, t, 'cax')
               call outpost(pipe%alpha, pipe%as, pipe%fshape, pr, t, 'geo')
               call nek_end()
            end if

         end subroutine helix_pipe
      
      end module neklab_helix
