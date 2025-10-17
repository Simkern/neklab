      module neklab_otd
      !---------------------------------------
      !-----     LightKrylov Imports     -----
      !---------------------------------------
         use stdlib_sorting, only: sort, sort_index
         use stdlib_linalg, only: eig, eye
         use stdlib_optval, only: optval
         use stdlib_logger, only: debug_level
      ! Default real kind.
         use LightKrylov, only: dp
      ! Abstract types for real-valued linear operators and vectors.
         use LightKrylov, only: abstract_linop_rdp, abstract_vector_rdp
         use LightKrylov, only: orthonormalize_basis, zero_basis, rand_basis, innerprod, Gram
         use LightKrylov_Utils, only: abstract_opts
      ! Logging
         use LightKrylov_Logger
      ! Extensions of the abstract vector types to nek data format.
         use neklab_vectors
         use neklab_linops, only: apply_L
         use neklab_utils, only: nek2vec, vec2nek, outpost_nek
         use neklab_nek_forcing, only: neklab_forcing, set_neklab_forcing
         use neklab_nek_setup, only: setup_linear_solver, setup_nonlinear_solver
         use neklab_nek_setup, only: nek_log_debug, nek_log_message, nek_stop_error
         implicit none
         include "SIZE"
         include "TOTAL"
         include "ADJOINT"
         private
         character(len=*), parameter, private :: this_module = 'neklab_linops'
         character(len=*), parameter, private :: logfile_Ls = 'Ls.dat'
         character(len=*), parameter, private :: logfile_Lr = 'Lr.dat'
      
         integer, parameter :: lv = lx1*ly1*lz1*lelv
      !! Local number of grid points for the velocity mesh.
         integer, parameter :: lp = lx2*ly2*lz2*lelv
      !! Local number of grid points for the pressure mesh.

         public :: new_OTD

         type, extends(abstract_linop_rdp), public :: nek_otd
            type(nek_dvector) :: baseflow
            type(nek_dvector), allocatable :: basis(:)
            real(dp),    dimension(lpert, lpert) :: Lr  = 0.0_dp
            !! Reduced operator
            real(dp),    dimension(lpert, lpert) :: Phi = 0.0_dp
            !! Internal rotation matrix
            ! Spectral decompositions
            ! Ls (symmetrized reduced operator)
            real(dp),    dimension(lpert, lpert) :: svec   = 0.0_dp
            !! eigenvectors
            real(dp),    dimension(lpert)        :: sigma  = 0.0_dp
            !! eigenvalues
            ! Lr
            complex(dp), dimension(lpert,lpert)  :: evec   = 0.0_dp
            !! eigenvectors
            complex(dp), dimension(lpert)        :: lambda = 0.0_dp
            !! eigenvalues
            ! Miscellaneous
            integer :: r = lpert
            !! Size of the reduced operator
            logical, private :: initialized = .false.
            !! Initialization flag
         contains
            private
            procedure, pass(self), public :: matvec => apply_LNS
            procedure, pass(self), public :: rmatvec => apply_adjLNS
            procedure, pass(self), public :: init
            !! Initialize basis and set nek status
            procedure, pass(self), public :: is_initialized
            !! Getter function for the initalisation status
            procedure, pass(self), public :: spectral_analysis
            !! Perform the spectral analysis of the reduced operator
            procedure, pass(self), public :: outpost_OTDmodes
            !! Project basis vectors onto eigendirectionf of reduced operator and outpost
            procedure, pass(self), public :: outpost_basis
            !! Outpost basis vectors directly (for restart)
            procedure, pass(self), public :: generate_forcing
            !! Generate the additional forcing for each basis vector for OTD formalism
            procedure, pass(self), public :: step
            !! Single step of the OTD algorithm following the time integration step
         end type nek_otd
      
         type, extends(abstract_opts), public :: otd_opts
            integer :: printstep = 5
            !! Output frequency for spectral analysis
            integer :: orthostep = 10
            !! Reorthogonalization frequency
            integer :: iostep = 100
            !! Output frequency for the projected basis vectors
            integer :: iorststep = 100
            !! Output frequency for the basis (for restarts)
            integer :: n_usrIC = 0
            !! Number of user-defined initial conditions for the perturbations
            real(dp) :: t_init = 0.0_dp
            !! Initial time
            logical :: trans = .false.
            !! Direct of adjoint?
            logical :: solve_baseflow = .true.
            !! Solve nonlinear equations alongside the linear equations?
            logical :: output_initial_conditions = .false.
            !! Outpost initial state (OTDinit)
            character(len=128) :: OTDIC_basename = 'OTDIC_'
            !! Base filename for initial conditions
         end type
      
      contains

         function new_OTD(bf) result(OTD)
            type(nek_dvector), intent(in) :: bf
            type(nek_otd) :: OTD
            OTD%baseflow = bf
            allocate(OTD%basis(lpert), source=bf)
            call zero_basis(OTD%basis)
            OTD%Lr     = 0.0_dp
            OTD%Phi    = 0.0_dp
            OTD%svec   = 0.0_dp
            OTD%sigma  = 0.0_dp
            OTD%evec   = (0.0_dp, 0.0_dp)
            OTD%lambda = (0.0_dp, 0.0_dp)
            OTD%r = lpert
            OTD%initialized = .false.
         end function new_OTD
      
         subroutine apply_LNS(self, vec_in, vec_out)
      ! Linear Operator.
            class(nek_otd), intent(inout) :: self
      ! Input vector.
            class(abstract_vector_rdp), intent(in) :: vec_in
      ! Output vector.
            class(abstract_vector_rdp), intent(out) :: vec_out
            select type (vec_in)
            type is (nek_dvector)
               select type (vec_out)
               type is (nek_dvector)
                  ! force baseflow
                  call vec2nek(vx, vy, vz, pr, t, self%baseflow)
                  ! apply linear operator
                  call apply_L(vec_out%vx, vec_out%vy, vec_out%vz,
     &                         vec_in%vx,  vec_in%vy,  vec_in%vz, 
     &                         vec_in%pr, trans = .false.)
               class default
                  call type_error('vec_out','nek_dvector','OUT', this_module, 'apply_LNS')
               end select
            class default
               call type_error('vec_in','nek_dvector','IN', this_module, 'apply_LNS')
            end select
            return
         end subroutine apply_LNS
      
         subroutine apply_adjLNS(self, vec_in, vec_out)
      ! Linear Operator.
            class(nek_otd), intent(inout) :: self
      ! Input vector.
            class(abstract_vector_rdp), intent(in) :: vec_in
      ! Output vector.
            class(abstract_vector_rdp), intent(out) :: vec_out
            select type (vec_in)
            type is (nek_dvector)
               select type (vec_out)
               type is (nek_dvector)
                  ! force baseflow
                  call vec2nek(vx, vy, vz, pr, t, self%baseflow)
                  ! apply linear operator
                  call apply_L(vec_out%vx, vec_out%vy, vec_out%vz,
     &                         vec_in%vx,  vec_in%vy,  vec_in%vz, 
     &                         vec_in%pr, trans = .true.)
               class default
                  call type_error('vec_out','nek_dvector','OUT', this_module, 'apply_adjLNS')
               end select
            class default
               call type_error('vec_in','nek_dvector','IN', this_module, 'apply_adjLNS')
            end select
            return
         end subroutine apply_adjLNS
      
         subroutine init(self, opts)
      ! Linear Operator.
            class(nek_otd), intent(inout) :: self
            type(otd_opts), intent(in) :: opts
      ! internal
            integer :: i, r
            real(dp) :: err
            real(dp), allocatable :: G(:,:)
            logical :: exist_IC
            character(len=128) :: ifile
            character(len=128) :: msg
            
      ! number of OTD modes
            r = self%r
            
            if (.not. self%is_initialized()) then
      ! switch to perturbation mode (needed for the BCs in the random initialization)
               jp = 1
      
      ! initialize random basis
               call zero_basis(self%basis); call rand_basis(self%basis)

      ! switch back to baseflow mode
               jp = 0
      
      ! check for user-defined ICs
               if (.not. opts%n_usrIC == 0) then
                  ! sanity check
                  if (opts%n_usrIC < 0) then
                     write (msg, *) 'Incorrect number of IC fields to load. nIC=', opts%n_usrIC
                     call nek_stop_error(msg, this_module, 'init_OTD')
                  else if (opts%n_usrIC > r) then
                     write (msg, *) 'Inconsistent number of IC fields to load. nIC=', opts%n_usrIC, ' r=', self%r
                     call nek_stop_error(msg, this_module, 'init_OTD')
                  end if
                  ! load ICs
                  do i = 1, opts%n_usrIC
                     write (ifile, '(A,I2.2,".fld")') trim(opts%OTDIC_basename), i
                     inquire (file=ifile, exist=exist_IC)
                     if (exist_IC) then
                        write (msg, *) 'Loading IC file: ', trim(ifile)
                        call nek_log_message(msg, this_module, 'init_OTD')
                        call load_fld(ifile)
                        call nek2vec(self%basis(i), vx, vy, vz, pr, t)
                     else
                        write (msg, *) 'Cannot find IC file: ', trim(ifile)
                        call nek_stop_error(msg, this_module, 'init_OTD')
                     end if
                  end do
               end if
      
      ! orthonormalize
               G = Gram(self%basis); err = maxval(abs(G - eye(r, mold=1.0_dp)))
               write (msg, '(A,I5,A,1X,E10.3)') 'Step ', istep, ': max. orthonormality error: ', err
               call nek_log_debug(msg, this_module, 'OTD init')
         
               call orthonormalize_basis(self%basis)
         
               G = Gram(self%basis); err = maxval(abs(G - eye(r, mold=1.0_dp)))
               write (msg, '(A,I5,A,1X,E10.3)') 'Step ', istep, ': max. orthonormality error: ', err
               call nek_log_debug(msg, this_module, 'OTD init')
      
      ! force baseflow
               call vec2nek(vx, vy, vz, pr, t, self%baseflow)

      ! set perturbations
               do i = 1, r
                  call vec2nek(vxp(:, i:i), vyp(:, i:i), vzp(:, i:i), prp(:, i:i), tp(:, :, i:i), self%basis(i))
               end do
      
      ! outpost initial conditions (for diagnostics)
               if (opts%output_initial_conditions) then
                  call outpost(vx, vy, vz, pr, t, 'bfi')
                  do i = 1, r
                     call outpost(vxp(1, i), vyp(1, i), vzp(1, i), prp(1, i), tp(1, :, i), 'pri')
                  end do
               end if
      
      ! setup nek5000
               call setup_linear_solver(recompute_dt=.true., cfl_limit=0.4_dp, solve_baseflow=opts%solve_baseflow)

      ! Set initial time
               time = opts%t_init

      ! Prepare logfiles
               if (nid == 0) then
                  open (1234, file=logfile_Ls, status='replace', action='write'); close (1234)
                  open (1234, file=logfile_Lr, status='replace', action='write'); close (1234)
               end if

      ! Mark initialization
               self%initialized = .true.
            end if
            return
         end subroutine init

         logical pure function is_initialized(self) result(OTD_is_initialized)
            class(nek_otd), intent(in) :: self
            OTD_is_initialized = self%initialized
         end function is_initialized
      
         subroutine spectral_analysis(self, ifprint)
      ! Linear Operator.
            class(nek_otd), intent(inout) :: self
            logical, intent(in) :: ifprint
      ! internals
            real(dp), dimension(self%r) :: s
            real(dp), dimension(self%r, self%r) :: Lsym
            complex(dp), dimension(self%r) :: l
            complex(dp), dimension(self%r, self%r) :: v
            integer :: i, r
            integer :: idx(self%r)
            character(len=128) :: msg, fmt_Lr
      
            r = self%r
            write (fmt_Lr, '("(I8,1X,F15.8,A,",I0,"(1X,E15.8),A,",I0,"(1X,E15.8))")') r, r
      
      ! compute eigenvalues of the symmetrized operator
            Lsym = 0.5*(self%Lr + transpose(self%Lr))
            call eig(Lsym, l, right=v)
            s = real(l)
            call sort_index(s, idx, reverse=.true.)
            do i = 1, r
               self%sigma(i) = l(idx(i))
               self%svec(:, i) = v(:, idx(i))
            end do
            call sort(self%sigma, reverse=.true.)
      ! stamp logfile
            if (ifprint) then
               write (msg, '(I10,1X,F15.8,*(1X,E15.8))') istep, time, self%sigma
               call nek_log_message(msg, this_module, 'OTD Ls')
               if (nid == 0) then
                  open (1234, file=logfile_Ls, status='old', action='write', position='append')
                  write (1234, '(I8,1X,F15.8,A,*(1X,E15.8))') istep, time, ' Ls ', self%sigma
                  close (1234)
               end if
            end if
      ! outpost projected modes?
      
      ! compute eigenvalues of Lr
            call eig(self%Lr, l, right=v)
            s = real(l)
            call sort_index(s, idx, reverse=.true.)
            do i = 1, r
               self%lambda(i) = l(idx(i))
               self%evec(:, i) = v(:, idx(i))
            end do
      ! stamp logfile
            if (ifprint) then
               write (msg, '(I7,1X,F15.8,*(1X,E15.8))') istep, time, real(self%lambda)
               call nek_log_message(msg, this_module, 'OTD Lr%Re')
               write (msg, '(I7,1X,F15.8,*(1X,E15.8))') istep, time, aimag(self%lambda)
               call nek_log_message(msg, this_module, 'OTD Lr%Im')
               if (nid == 0) then
                  open (1234, file=logfile_Lr, status='old', action='write', position='append')
                  write (1234, fmt_Lr) istep, time, ' Lr%Re ', real(self%lambda), ' Lr%Im ', aimag(self%lambda)
                  close (1234)
               end if
            end if
            return
         end subroutine spectral_analysis
      
         subroutine outpost_OTDmodes(self, outpost_imaginary)
      !! We assume the OTD basis has already been copied to v[xyz]p
      ! Linear Operator.
            class(nek_otd), intent(inout) :: self
            logical, optional, intent(in) :: outpost_imaginary
      ! internal
            real(dp), dimension(lv, self%r) :: vxr
            real(dp), dimension(lv, self%r) :: vyr
            real(dp), dimension(lv, self%r) :: vzr
      
            integer :: i, r
            logical :: outpost_imaginary_
            character(len=3) :: file_prefix
      
            r = self%r
            outpost_imaginary_ = optval(outpost_imaginary, .false.)

      ! project modes (real part)
            call           mxm(vxp, lv, real(self%evec), r, vxr, r)
            call           mxm(vyp, lv, real(self%evec), r, vyr, r)
            if (if3d) call mxm(vzp, lv, real(self%evec), r, vzr, r)
            do i = 1, self%r
               write (file_prefix, '(A,I2.2)') 'm', i
               call outpost(vxr(1, i), vyr(1, i), vzr(1, i), prp(1, i), tp(1, :, i), file_prefix)
            end do

      ! project modes (imaginary part)
            if (outpost_imaginary_) then
               call           mxm(vxp, lv, aimag(self%evec), r, vxr, r)
               call           mxm(vyp, lv, aimag(self%evec), r, vyr, r)
               if (if3d) call mxm(vzp, lv, aimag(self%evec), r, vzr, r)
               do i = 1, self%r
                  write (file_prefix, '(A,I2.2)') 'i', i
                  call outpost(vxr(1, i), vyr(1, i), vzr(1, i), prp(1, i), tp(1, :, i), file_prefix)
               end do
            end if
            return
         end subroutine outpost_OTDmodes

         subroutine outpost_basis(self)
      ! Linear Operator.
            class(nek_otd), intent(inout) :: self
            integer :: i
            character(len=3) :: file_prefix
      
            do i = 1, self%r
               write (file_prefix, '(A,I2.2)') 'r', i
               associate (bvec => self%basis(i))
                  call outpost_nek(bvec, file_prefix)
               end associate
            end do
            return
         end subroutine outpost_basis
      
         subroutine generate_forcing(self)
      !! We assume the OTD basis has already been copied to v[xyz]p
      ! Linear Operator.
            class(nek_otd), intent(inout) :: self
      ! internal
            integer :: r, i
            real(dp), dimension(lv, self%r) :: OTDfx
            real(dp), dimension(lv, self%r) :: OTDfy
            real(dp), dimension(lv, self%r) :: OTDfz
            r = self%r
            call mxm(vxp, lv, self%Lr - self%Phi, r, OTDfx, r)
            call mxm(vyp, lv, self%Lr - self%Phi, r, OTDfy, r)
            if (if3d) call mxm(vzp, lv, self%Lr - self%Phi, r, OTDfz, r)
            do i = 1, r
               call set_neklab_forcing(OTDfx(:, i), OTDfy(:, i), OTDfy(:, i), ipert=i)
            end do
            return
         end subroutine generate_forcing

         subroutine step(self, opts, istp, Lu)
            class(nek_otd), intent(inout) :: self
            type(otd_opts), intent(in) :: opts
            integer, intent(in) :: istp
            type(nek_dvector), dimension(:), intent(inout) :: Lu
            ! internal
            character(len=*), parameter :: this_procedure = 'OTD_step'
            real(dp),    dimension(lpert, lpert) :: Phi, G
            ! Misc
            integer :: i, j, r, log_level
            real(dp) :: err
            character(len=128) :: msg

            r = self%r

            ! Initialize if not done
            call self%init(opts)

            ! Update baseflow if solved
            if (opts%solve_baseflow) call nek2vec(self%baseflow, vx, vy, vz, pr, t)
      
      ! load perturbations into neklab
            do i = 1, r
               call nek2vec(self%basis(i), vxp(:, i:i), vyp(:, i:i), vzp(:, i:i), prp(:, i:i), tp(:, :, i:i))
            end do
      
      ! orthonormalize
            if ((istp <= 10) .or.
     &              mod(istp, opts%orthostep) == 0 .or.
     &              mod(istp, opts%printstep) == 0 .or.
     &              mod(istp, opts%iostep) == 0) then
   
               if (log_level <= debug_level) then
                  G = Gram(self%basis); err = maxval(abs(G - eye(r, mold=1.0_dp)))
                  write (msg, '(A,I5,A,1X,E10.3)') 'Step ', istp, ': max. orthonormality error: ', err
                  call nek_log_debug(msg, this_module, this_procedure)
               end if
                  
               call orthonormalize_basis(self%basis)
                  
               if (log_level <= debug_level) then
                  G = Gram(self%basis); err = maxval(abs(G - eye(r, mold=1.0_dp)))
                  write (msg, '(A,I5,A,1X,E10.3)') 'Step ', istp, ': max. orthonormality error: ', err
                  call nek_log_debug(msg, this_module, this_procedure)
               end if
            end if

      ! compute Lu
            do i = 1, r
               if (opts%trans) then
                  call self%apply_rmatvec(self%basis(i), Lu(i))
               else
                  call self%apply_matvec(self%basis(i), Lu(i))
               end if
            end do

      ! compute reduced operator
            self%Lr = innerprod(self%basis, Lu)

      ! compute internal rotation
            self%Phi = 0.0_dp
            do i = 1, r
               do j = i + 1, r
                  self%Phi(i, j) =  self%Lr(i, j)
                  self%Phi(j, 1) = -self%Lr(i, j)
               end do
            end do

      ! output projected modes
            if (mod(istp, opts%printstep) == 0) then
               call self%spectral_analysis(ifprint=.true.)
            end if

      ! project basis vectors and output modes
            if (mod(istp, opts%iostep) == 0) then
               if (.not. mod(istp, opts%printstep) == 0) then ! ensure spectral data are up to date
                  call self%spectral_analysis(ifprint=.false.)
               end if
               call self%outpost_OTDmodes()
               if (opts%solve_baseflow) call outpost_nek(self%baseflow, 'bf_')
            end if
            
      ! output basis vectors for restart
            if (mod(istp, opts%iorststep) == 0) then
               call self%outpost_basis()
               if (opts%solve_baseflow) call outpost_nek(self%baseflow, 'rbf')
            end if
               
      ! set the forcing
            call self%generate_forcing()

      ! copy data back to nek
            do i = 1, r
               call vec2nek(vxp(:, i:i), vyp(:, i:i), vzp(:, i:i), prp(:, i:i), tp(:, :, i:i), self%basis(i))
            end do

         end subroutine
      
      end module neklab_otd
