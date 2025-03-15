      module neklab_helix
      !---------------------------------------
      !-----     LightKrylov Imports     -----
      !-------------------Compute the --------------------
      ! Default real kinCompute the d.
      ! Logging & timing
         use LightKrylov_LoggerCompute the 
      ! Extensions of the abstract vector types to nek data format.
Compute the          use neklab_vectors

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
         integer, parameter :: lbuf = 1000
      !! Maximum number of 2d fields to save before outposting
         integer, parameter, public :: lf = 3
      !! Maximum number of forcing components 1 x steady + 2 x (# unsteady)
         integer, parameter, public :: nfft = 16
      !! Number of FT components to compute

         public :: pipe
         public :: helix_pipe ! constructor for the pipe instance of the helix type

         type, public :: helix
            !! Type containing the geometrical and dynamical properties of a helix
            private
            ! geometry inputs
            real(dp) :: length      ! of the torus centerline
            real(dp) :: delta       ! curvature ratio R/r
            real(dp) :: diameter    ! of the torus 
            real(dp) :: pitch_s     ! angle of the helix (set to 0.0 for torus)
            ! derived geometry quantities
            real(dp) :: radius
            real(dp) :: curv_radius ! distance from helix center to pipe center
            real(dp) :: phi         ! Rise angle of helix w.r.t. equatorial plane
            real(dp) :: sweep       ! Total angle swept by helix in streamwise direction
            ! flow
            logical :: if_steady
            real(dp) :: pulse_T     ! pulsation period
            real(dp) :: omega       ! pulsation frequency
            real(dp) :: womersley   ! womersley number r (omega/nu)^(1/2)
            integer  :: nf          ! number of forcing Fourier components
            ! forcing
            real(dp), dimension(lf) :: dpds
            ! Fourier components of the streamwise pressure gradient
            real(dp), dimension(lx1,ly1,lz1,lelv) :: fshape
            ! Spatial distribution of the forcing in the cross-stream plane due to curvature
            ! mesh inputs
            integer :: nslices            ! number of slices in streamwise direction
            integer :: nelf               ! number of elements on a facCompute the e (cross-stream plane)
            logical :: if_torus = .false. ! is the mesh curved (toroidal)
            logical :: if_helix = .false. ! is the mesh helical?
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
            procedure, pass(self), public :: forcing_amplitude
            procedure, pass(self), public :: compute_bf_forcing
            procedure, pass(self), public :: compute_usrt
            procedure, pass(self), public :: compute_ubar
            procedure, pass(self), public :: shift_mflow_phase
            procedure, pass(self), public :: gfldr_torus
            procedure, pass(self), public :: setup_summary
            procedure, pass(self), public :: parameter_summary
            procedure, pass(self), public :: forcing_summary
            ! helix_IO
            procedure, pass(self), public :: fname_2d
            procedure, pass(self), public :: write_2d
            procedure, pass(self), public :: read_2d
            procedure, pass(self) :: get_nsteps_from_header
            ! helix_2d
            procedure, pass(self) :: init_2d_geom
            procedure, pass(self), public :: save_2d_fields
            procedure, pass(self), public :: outpost_2d_fields
            procedure, pass(self), public :: load_2d_fields
            procedure, pass(self), public :: set_baseflow
            procedure, pass(self), public :: load_baseflow
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

            module subroutine init_geom(self, if_debug)
               !! Initialize the geometry by morphing the input geometry (straight pipe) 
               !! into a torus Compute the or a helix with the correct dimensions
               class(helix), intent(inout) :: self
               logical, optional, intent(in) :: if_debug
               !! Debug flag: Will outpost diagnostic fields and exit
            end subroutine init_geom

            module subroutine init_flow(self, dpds, womersley, reset_nf)
               !! Initialize the forcing components and set the derived physical quantities
               !! When called for the first time or if reset_nf == .true., the number of
               !! considered forcing components is updated based on dpds. 
               !! Othersise, providing a pressure gradient vector with the wrong length throws an error.
               class(helix), intent(inout) :: self
               real(dp), intent(in) :: dpds(:)
               !! Fourier expansion of the streamwise pressure gradient
               real(dp), optional, intent(in) :: womersley
               !! Pulsation frequency
               logical, optional, intent(in) :: reset_nf
               !! Reset the number of considered forcing components (self%nf) based on input array? (default: .false.)
            end subroutine init_flow

            module subroutine compute_fshape(self)
               !! Computes the spatial distribution in the cross-stream plane
               !! of the streamwise forcing based on the geometry (self%fshape)
               class(helix), intent(inout) :: self
            end subroutine compute_fshape
            
            module pure function forcing_amplitude(self, t) result(f)
               !! Construct the instantaneous streamwise pressure gradient from the Fourier expansion
               class(helix), intent(in) :: self
               real(dp), intent(in) :: t
               !! time
               real(dp) :: f
               !! output forcing amplitude
            end function forcing_amplitude

            module subroutine compute_bf_forcing(self, t)
               !! Compute and set the (instantaneous) cartesian forcing components ffx,ffy,ffz 
               !! into the neklab_forcing array to be loaded in userf
               class(helix), intent(in) :: self
               real(dp) :: t
               !! time
            end subroutine compute_bf_forcing

            module subroutine compute_usrt(self, u, v, w, us, ur, ut)
               !! Switch from cartesian (nek5000) to toroidal/helical (post-processing) coordinates
               !! for the velocity components
               class(helix), intent(in) :: self
               real(dp), dimension(lx1,ly1,lz1,lelv), intent(in) :: u
               !! X-velocity in cartesian coordinates
               real(dp), dimension(lx1,ly1,lz1,lelv), intent(in) :: v
               !! Y-velocity in cartesian coordinates
               real(dp), dimension(lx1,ly1,lz1,lelv), intent(in) :: w
               !! Z-velocity in cartesian coordinates
               real(dp), dimension(lx1,ly1,lz1,lelv), intent(out) :: us
               !! Streamwise velocity in toroidal coordinates
               real(dp), dimension(lx1,ly1,lz1,lelv), intent(out) :: ur
               !! radial velocity in toroidal coordinates
               real(dp), dimension(lx1,ly1,lz1,lelv), intent(out) :: ut
               !! Azimuthal velocity in toroidal coordinates
            end subroutine compute_usrt
            
            module function compute_ubar(self,u,v,w) result(ubar)
               !! Compute the streamwise mass flow rate in the torus
               class(helix), intent(in) :: self
               real(dp), dimension(lx1,ly1,lz1,lelv), intent(in) :: u
               !! X-velocity in cartesian coordinates
               real(dp), dimension(lx1,ly1,lz1,lelv), intent(in) :: v
               !! Y-velocity in cartesian coordinates
               real(dp), dimension(lx1,ly1,lz1,lelv), intent(in) :: w
               !! Z-velocity in cartesian coordinates
               real(dp) :: ubar
               !! output mass flow rate
            end function compute_ubar

            module subroutine shift_mflow_phase(self, icomp, target_mflow_phase)
               !! Rotate the Fourier component icomp in order to achieve the target input mass flow phase
               class(helix), intent(inout) :: self
               integer, intent(in) :: icomp
               !! Index of the complex Fourier component to shift (note: not the index in the (real) dpds array!)
               real(dp), intent(in) :: target_mflow_phase
               !! Target phase of the mass flow rate
            end subroutine shift_mflow_phase

            module subroutine gfldr_torus(self, rstfname)
               !! Interpolate a streamwise-independent flow field from a different mesh onto
               !! the current mesh, irrespective of the length of the flow domain in each case
               !! Internally, we extract the 2D data and use it to set the 3D fields.
               class(helix), intent(inout) :: self
               character(len=*), intent(in) :: rstfname
               !! Restart file on a different toroidal mesh
               !! Note: The restart file must correspond to the same geometry (apart from the streamwise domain length)
            end subroutine gfldr_torus

            module subroutine setup_summary(self)
               !! Print the geometry and mesh parameters/settings to log
               class(helix), intent(in) :: self
            end subroutine setup_summary
            
            module subroutine parameter_summary(self)
               !! Print the flow parameters to log
               class(helix), intent(in) :: self
            end subroutine parameter_summary

            module subroutine forcing_summary(self)
               !! Print the Fourier expansion of the streamwise pressure gradient to log
               class(helix), intent(in) :: self
            end subroutine forcing_summary

            !-----------------------------------------------------
            ! neklab_helix % helix_IO
            !
            ! Type-bound procedures

            module pure function fname_2d(self, iname, iout) result(fname)
               !! Construct the standard 2D data file name
               class(helix), intent(in) :: self
               character(len=1), intent(in) :: iname
               !! File identifier
               integer, intent(in) :: iout
               !! File index
               character(len=132) :: fname
               !! Output filename
            end function fname_2d

            module subroutine write_2d(self, fname, only_mesh)
               !! Write the data stored in self%v[xyz]2d to file in binary format
               class(helix), intent(in) :: self
               character(len=132), intent(in) :: fname
               !! Filename of the output file
               logical, optional, intent(in) :: only_mesh
               !! Write only the mesh to file with no velocity data (default: .false.)
            end subroutine write_2d
            
            module subroutine read_2d(self, fname)
               !! Read binary file containing the 2D velocity fields and store them in self%v[xyz]2d
               class(helix), intent(inout) :: self
               character(len=132), intent(in) :: fname
               !! Filename of the input file
            end subroutine read_2d
            
            module subroutine get_nsteps_from_header(self, fname, nsaver)
               !! Extract the number of stored snapshots from the binary file
               class(helix), intent(in) :: self
               character(len=132), intent(in) :: fname
               !! Filename of the input file
               integer, intent(out) :: nsaver
               !! Output number of saved timesteps in file
            end subroutine get_nsteps_from_header
            
            !-----------------------------------------------------
            ! neklab_helix % helix_2d
            !
            ! Type-bound procedures
            
            module subroutine init_2d_geom(self, if_debug)
               !! Initialize the 2D geometry data on the cross-stream plane and identify the elements
               !! that are aligned in streamwise direction forming a unique segment.
               !! Establish local and global sCompute the egment ownership to minimize MPI communication.
               class(helix), intent(inout) :: self
               logical, optional, intent(in) :: if_debug
               !! Outpost and print debug information.
            end subroutine init_2d_geom
            
            module subroutine save_2d_fields(self, u, v, w)
               !! Extract the velocity data from the global segment owner and save it in self%v[xyz]2d
               class(helix), intent(inout) :: self
               real(dp), dimension(lx1,ly1,lz1,lelv), intent(in) :: u
               !! X-component of the velocity
               real(dp), dimension(lx1,ly1,lz1,lelv), intent(in) :: v
               !! Y-component of the velocity
               real(dp), dimension(lx1,ly1,lz1,lelv), intent(in) :: w
               !! Z-component of the velocity
            end subroutine save_2d_fields
            
            module subroutine outpost_2d_fields(self)
               !! Wrapper function to write 2D binary data to file
               !! The output filename is based on the type of calculation (Newton, Floquet, ...)
               class(helix), intent(inout) :: self
            end subroutine outpost_2d_fields
            
            module subroutine load_2d_fields(self, idx)
               !! Wrapper funcion to read 2D binary data from file
               !! The input filename is based on the type of calculation (Newton, Floquet, ...)
               class(helix), intent(inout) :: self
               integer, intent(in) :: idx
               !! File index
            end subroutine load_2d_fields
            
            module subroutine set_baseflow(self, basex, basey, basez, ifld)
               !! Copy the 2D data snapshot ifld from self%v[xyz]2D to the 3D fields basex, basey, basez
               class(helix), intent(inout) :: self
               real(dp), dimension(lx1,ly1,lz1,lelv), intent(out) :: basex
               !! X-component of the velocity to be set
               real(dp), dimension(lx1,ly1,lz1,lelv), intent(out) :: basey
               !! Y-component of the velocity to be set
               real(dp), dimension(lx1,ly1,lz1,lelv), intent(out) :: basez
               !! Z-component of the velocity to be set
               integer, intent(in) :: ifld
               !! index of the timestep to be set
            end subroutine set_baseflow

            module subroutine load_baseflow(self, basex, basey, basez, fname, ifld)
               !! Wrapper function to load and set the 3D fields basex, basey, basez from a file
               class(helix), intent(inout) :: self
               real(dp), dimension(lx1,ly1,lz1,lelv), intent(out) :: basex
               !! X-component of the velocity to be set
               real(dp), dimension(lx1,ly1,lz1,lelv), intent(out) :: basey
               !! Y-component of the velocity to be set
               real(dp), dimension(lx1,ly1,lz1,lelv), intent(out) :: basez
               !! Z-component of the velocity to be set
               character(len=*), intent(in) :: fname
               !! Filename of input file
               integer, optional, intent(in) :: ifld
               !! Index of the timestep to be set (default: 1)
            end subroutine load_baseflow

            module subroutine compute_2d_usrt(self)
               !! Switch from cartesian (nek5000) to toroidal/helical (post-processing) coordinates
               !! for the 2D velocity data
               class(helix), intent(inout) :: self
            end subroutine compute_2d_usrt

            module subroutine set_2d_mode(self, mode)
               !! Utility function to set the logical flags
               character(len=*), intent(in) :: mode
               !! Currently available modes:
               !! 1. 'newton'
               !! 2. 'floquet'
               class(helix), intent(inout) :: self
            end subroutine set_2d_mode

            !-----------------------------------------------------
            ! neklab_helix % helix_mflow_fft
            !
            ! Type-bound procedures 

            module subroutine reset_mflow_fft(self)
               !! Reset the internal variables to restart integration
               class(helix), intent(inout) :: self
            end subroutine reset_mflow_fft

            module subroutine compute_mflow_fft(self, period, var_dt)
               !! Aggregate the mass flow data to compute the spectrum
               !! This function should be called at every timestep during integration
               !! Support for fixed and variable timestep integration.
               class(helix), intent(inout) :: self
               real(dp), optional, intent(in) :: period
               !! Period of the signal (default: self%pulse_T)
               logical, optional, intent(in) :: var_dt
               !! Is the timestep variable? (default: .false.)
            end subroutine compute_mflow_fft
            
            module subroutine extract_mflow_fft(self, period, if_amplitude)
               !! Compute and extract the Fourier spectrum of the mass flow rate
               !! This function must be called when the period matches the integration time
               !! The phase and amplitude of the spectrum are computed.
               !! The timing consistency is checked.
               class(helix), intent(inout) :: self
               real(dp), optional, intent(in) :: period
               !! Period of the signal (default: self%pulse_T)
               logical, optional, intent(in) :: if_amplitude
               !! Print result using real numbers (amplitude and phase) or using complex numbers? (default: complex)
            end subroutine extract_mflow_fft

            module subroutine get_mflow_fft(self, mflow, phase, if_amplitude)
               !! Utility function to return the mass flow Fourier expansion once it is available
               class(helix), intent(in) :: self
               real(dp), allocatable, intent(out) :: mflow(:)
               !! output mass flow rate for each (complex) Fourier component
               real(dp), optional, allocatable, intent(out) :: phase(:)
               !! optional output phase of the mass flow rate for each (complex) Fourier component
               logical, optional, intent(in) :: if_amplitude
               !! Return data using amplitude and phase or complex conjugate fourier coefficients (default: .true.)
            end subroutine get_mflow_fft

            !-----------------------------------------------------
            ! neklab_helix % helix_gs
            !
            ! Type-bound procedures

            ! mode switches : Setter functions for the system switches

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
            
            ! logicals : Getter functions for the system switches

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
               !! Is the current rank the local owner of the segment the local element is part of?
               class(helix), intent(in) :: self
               integer, intent(in) :: ie
               !! local element index
               logical :: is_owner
            end function is_lowner
            
            module pure function is_gowner(self, ie) result(is_owner)
               !! Is the current rank the global owner of the segment the local element is part of?
               class(helix), intent(in) :: self
               integer, intent(in) :: ie
               !! local element index
               logical :: is_owner
            end function is_gowner

            ! Getter functions for private system parameters/values

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
               real(dp), dimension(lf), intent(out) :: dpds
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
               !! Return the number of timesteps for a full period
               class(helix), intent(in) :: self
               integer :: ns
            end function get_nsteps

            module subroutine get_dt_minmax(self, dt_minmax)
               !! Return the extremal timestep values over the period to estimate integration errors
               class(helix), intent(in) :: self
               real(dp), dimension(2), intent(out) :: dt_minmax
            end subroutine get_dt_minmax

            module function get_ubar_lag(self) result(ubar_lag)
               !! Return the mass flow rate value from the previous timestep
               class(helix), intent(in) :: self
               real(dp) :: ubar_lag
            end function get_ubar_lag

            module pure function get_lsegment(self, ie) result(local_segment)
               !! Return the local segment index for the current element
               class(helix), intent(in) :: self
               integer, intent(in) :: ie
               integer :: local_segment
            end function get_lsegment
         
            module pure function get_gsegment(self, ie) result(global_segment)
               !! Return the global segment index for the current element
               class(helix), intent(in) :: self
               integer, intent(in) :: ie
               integer :: global_segment
            end function get_gsegment

            module pure function get_v2d(self,ix,iy,iseg,ifld,icomp) result(v2d)
               !! Return the 2D velocity value for a speific point, segment, timestep and component
               class(helix), intent(in) :: self
               integer, intent(in) :: ix
               !! index in the lx1 direction
               integer, intent(in) :: iy
               !! index in the ly1 direction
               integer, intent(in) :: iseg
               !! local segment index
               integer, intent(in) :: ifld
               !! timestep index
               integer, intent(in) :: icomp
               !! velocity component index (1, 2, 3 for u, v, w)
               real(dp) :: v2d
               !! output velocity
            end function get_v2d

            ! Setter functions for private system parameters

            module subroutine set_dpds(self, dpds, reset)
               !! Update the streamwise pressure gradient
               class(helix), intent(inout) :: self
               real(dp), dimension(lf), intent(in) :: dpds
               logical, optional, intent(in) :: reset
            end subroutine set_dpds

            module subroutine set_nsteps(self, ns)
               !! Set the number of timesteps per period
               class(helix), intent(inout) :: self
               integer, intent(in) :: ns
            end subroutine set_nsteps
            
         end interface

         type(helix) :: pipe

      contains

         subroutine helix_pipe(delta, diameter, pitch_s, length, if_sym, if_debug)
            !! Constructor for the module level instance of helix
            !! This is only called once prior to computation after the boundary conditions
            !! are set to initialize the geometry and check for consistency.
            !!    1. Check that the chosen setup and mesh are consistent
            !!    2. Initialize timers
            !!    3. Morph the geometry into a torus/helix
            !!    4. Initialize spatial forcing distribution
            real(dp), intent(in) :: delta
            !! torus/helix cuvature ratio
            real(dp), intent(in) :: diameter
            !! torus/helix diameter
            real(dp), intent(in) :: pitch_s
            !! helix pitch angle (set to 0.0 for torus)
            real(dp), intent(in) :: length
            !! length of the helix/torus at the centerline
            logical, optional, intent(in) :: if_sym
            !! Does the mesh cover only half the torus/helix? (default = .false.)
            logical, optional, intent(in) :: if_debug
            !! Debug output and exit. (default = .false.)
            ! internal
            integer :: ie, iface
            logical :: debug, symmetry
            real(dp) :: xmin
            integer :: nsym
            character(len=128) :: msg
            ! functions
            real(dp), external :: glmin
            integer, external :: iglsum
            
            ! Optional debug argument
            debug = optval(if_debug, .false.)
            symmetry = optval(if_sym, .false.)
            
            ! Case setup
            call pipe%set_symmetry(symmetry)
            
            ! Sanity checks
            if (symmetry) then
               call nek_log_message('Symmetric half torus setup defined.', this_module, 'helix_pipe')
            else
               call nek_log_message('Full torus setup defined.', this_module, 'helix_pipe')
            end if
            ! are symmetry conditions set in symmetric case?
            nsym = 0
            do ie = 1, nelv
               do iface = 1, 2*ndim
                  if (cbc(iface,ie,1) .eq. 'SYM') nsym = nsym + 1
               end do
            end do
            nsym = iglsum(nsym, 1)
            if (symmetry .and. nsym == 0) then
               call nek_stop_error('Mesh does not have symmetry conditions.', this_module, 'helix_pipe')
            else if (.not. symmetry .and. nsym > 0) then
               call nek_stop_error('Mesh has symmetry conditions.', this_module, 'helix_pipe')
            end if
            ! is the mesh consistent with the setup?
            xmin = glmin(xm1,lv)
            if ((symmetry .and. xmin < -0.1_dp) .or. (.not. symmetry .and. xmin > -0.1_dp)) then
               call nek_stop_error('Case setup and mesh domain are inconsistent.', this_module, 'helix_pipe')
            end if

            ! Geometry
            pipe%delta    = delta
            pipe%diameter = diameter
            pipe%pitch_s  = pitch_s
            pipe%length   = length

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
            if (pipe%phi /= 0.0_dp) then
               pipe%sweep    = pipe%length*cos(pipe%phi)/pipe%curv_radius ! sweep angle in radians
               pipe%if_helix = .true.
            else
               pipe%sweep    = pipe%length/pipe%curv_radius ! sweep angle in radians
               pipe%if_helix = .false.
            end if

            ! add timers
            call lk_timer%initialize() ! in case it has not been done
            call lk_timer%add_timer('neklab_helix_init_geom',         start=.false.)
            call lk_timer%add_timer('neklab_helix_init_2d_geom',      start=.false.)
            call lk_timer%add_timer('neklab_helix_save_2d_fields',    start=.false.)
            call lk_timer%add_timer('neklab_helix_read_2d',           start=.false.)
            call lk_timer%add_timer('neklab_helix_write_2d',          start=.false.)
            call lk_timer%add_timer('neklab_helix_set_baseflow',      start=.false.)
            call lk_timer%add_timer('neklab_helix_compute_ubar',      start=.false.)
            call lk_timer%add_timer('neklab_helix_compute_mflow_fft', start=.false.)


            ! intialize geometry
            call pipe%init_geom(if_debug)
            ! compute forcing distribution
            call pipe%compute_fshape()

            if (debug) then
               call outpost(pipe%xax, pipe%yax, pipe%zax, pr, t, 'cax')
               call outpost(pipe%alpha, pipe%as, pipe%fshape, pr, t, 'geo')
               call nek_end()
            end if

         end subroutine helix_pipe
      
      end module neklab_helix
