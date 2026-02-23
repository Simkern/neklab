      module neklab_analysis_torus
         use stdlib_stats_distribution_normal, only: normal => rvs_normal
         use stdlib_strings, only: padl, padr
         use stdlib_optval, only: optval
         use stdlib_linalg, only: diag, eye, det, inv
         use stdlib_logger, only: information_level, warning_level, debug_level, error_level, all_level, success
         use LightKrylov, only: atol_dp, dp, eigs, svds, save_eigenspectrum
         use LightKrylov, only: kexpm, gmres_rdp
         use LightKrylov, only: initialize_krylov_subspace, orthonormalize_basis, zero_basis, rand_basis
         use LightKrylov, only: linear_combination, innerprod
         use LightKrylov, only: newton, newton_dp_opts
         use LightKrylov_Logger
         use LightKrylov_Constants, only: io_rank
         use LightKrylov_Timing, only: timer => global_lightkrylov_timer
         use LightKrylov_AbstractSystems, only: abstract_system_rdp
         use neklab_vectors
         use neklab_linops
         use neklab_utils
         use neklab_nek_setup
         use neklab_otd
         use neklab_systems
         use neklab_helix
         use neklab_analysis
      
         implicit none
         include "SIZE"
         include "TOTAL"
         include "ADJOINT"
      
         private
         character(len=*), parameter, private :: this_module = 'neklab_analysis'
      
         public :: linear_stability_analysis_periodic_orbit
         public :: mflow_newton
         public :: shift_mflow_phase_torus
         public :: compute_nonlinear_period_torus
         public :: compute_monodromy_period_torus
         public :: compute_energy_budgets_period_torus
      
      contains

         subroutine linear_stability_analysis_periodic_orbit(floquet_operator, kdim, nev, adjoint, X0, tol)
            type(floquet_linop), intent(inout) :: floquet_operator
      !! Floquet operator whose stability properties are to be investigated.
            integer, intent(in) :: kdim
      !! Maximum dimension of the Krylov subspace.
            integer, intent(in) :: nev
      !! Desired number of eigenpairs to converge.
            logical, intent(in), optional :: adjoint
      !! Whether direct or adjoint analysis should be conducted.
		   type(nek_dvector), optional, intent(in) :: X0
	   !! Initial guess for the eigenvectors
            real(dp), optional, intent(in) :: tol
      !! Tolerance for the eigenvalue convergence
      
      ! Eigenvalue computation related variables.
            type(nek_dvector), allocatable :: eigvecs(:)
            complex(kind=dp), allocatable :: eigvals(:)
            real(kind=dp), allocatable :: residuals(:)
            integer :: info
      
      ! Miscellaneous.
            character(len=*), parameter :: this_procedure = 'stability_PO_main'
            real(kind=dp) :: alpha
            integer :: i
            logical :: adjoint_
            character(len=3) :: file_prefix
      
      ! Set up logging
            call logger_setup(nio=0, log_level=information_level, log_stdout=.false., log_timestamp=.true.)
      
      ! Optional parameters.
            if (present(adjoint)) then
               adjoint_ = adjoint
            else
               adjoint_ = .false.
            end if
      
      ! Allocate eigenvectors and initialize Krylov basis.
            allocate (eigvecs(nev)); call zero_basis(eigvecs)
      
      ! Run the eigenvalue analysis.
		      call eigs(floquet_operator, eigvecs, eigvals, residuals, info, 
     &                  x0=X0, kdim=kdim, transpose=adjoint_, write_intermediate=.true.)
      
      ! Transform eigenspectrum to continuous-time representation.
            eigvals = log(eigvals)/floquet_operator%tau
      
      ! Determine the file prefix.
            file_prefix = merge("adj", "dir", adjoint_)
      
      ! Save eigenspectrum to disk.
            call save_eigenspectrum(eigvals, residuals, trim(file_prefix)//"_eigenspectrum.npy")
      
      ! Export eigenfunctions to disk.
            call outpost_dnek(eigvecs(:nev), file_prefix)

		      call log_message('Exiting eigenvalue computation.', this_module, this_procedure)

      ! Finalize exptA timings
            call floquet_operator%finalize_timer()
      ! Finalize timing
            call logger_setup(logfile='lightkrylov_tmr.log', nio=0, log_level=warning_level, log_stdout=.false., log_timestamp=.true.)
            call timer%finalize()
      
         end subroutine linear_stability_analysis_periodic_orbit

         subroutine mflow_newton(sys, bf, mflow_target, tol, tol_mf, tol_mode, maxiter_newton)
            class(abstract_system_rdp), intent(inout) :: sys
      !! System for which a fixed point is sought
            type(nek_dvector), intent(inout) :: bf
      !! Initial guess for the fixed point
            real(dp), dimension(:), intent(in) :: mflow_target
      !! Target values for the mass flow rate for each Fourier component
            real(dp), intent(inout) :: tol
      !! Absolute tolerance for the Newton solver
            real(dp), intent(inout) :: tol_mf
      !! Absolute tolerance for mass flow rate error
            integer, optional, intent(in) :: tol_mode
      !! Use constant or dynamic tolerances in the Newton-Krylov solver
            integer, optional, intent(in) :: maxiter_newton
      !! Maximum number of newton steps to converge the mass flow rate
      ! internal
            character(len=*), parameter :: this_procedure = 'mflow_newton_main'
            type(nek_dvector), allocatable :: ref
            logical :: is_fp
            integer :: tol_mode_, maxiter_newton_
            integer :: nmf, nf, inwt, i, j, icnt
            real(dp) :: Wo, df0
            real(dp) :: dpds(lf)
            real(dp), allocatable :: phase(:), dpds_tmp(:)
            real(dp), allocatable :: mflow_old(:), mflow_new(:)
            real(dp), allocatable :: dmf(:), mf_err(:)
            real(dp), allocatable :: deltaf(:), fpert(:)
            real(dp), allocatable :: jac(:,:)
				real(dp) :: dt_minmax(2)
            real(dp) :: tol_mf_inexact, tol_df
            character(len=128) :: msg, fmt
            character(len=2048) :: long_msg
            character(len=10) :: step_id
            character(len=18) :: coef_id
            integer, parameter :: pad = 10

      ! optional arguments
            tol_mode_       = optval(tol_mode, 1)
            maxiter_newton_ = optval(maxiter_newton, 10)
      ! preparation & checks
            Wo = pipe%get_Wo()
            nmf = size(mflow_target) ! number of mass flow Fourier components to converge
	      write(fmt,'("(A,",I0,"(1X,F16.10),A,E16.8)")') nmf
            nf = pipe%get_nf() ! number of real forcing components (real and imaginary parts counted individually)
            if (nmf > nfft) then
               write(msg,'(A,I0,A,I0,A)') 'nmf= ', nmf, ' > nfft= ', nfft,'. Increase nfft in neklab_helix.'
               call nek_stop_error(msg, this_module, this_procedure)
            else if (2*nmf - 1 /= nf) then
               write(msg,'(A,I0,A,I0,A)') 'nmf= ', nmf, ' and nf= ', nf, ' incompatible.'
               call nek_stop_error(msg, this_module, this_procedure)
            else
               write(msg,'(A,I0,A)') 'Starting Newton iteration for nmf = ', nmf, ' mass flow components.'
               call nek_log_information(msg, this_module, this_procedure)
            end if
            allocate(dpds_tmp(nf))
            allocate(dmf(nmf), mf_err(nmf), deltaf(nmf), fpert(nmf))
            allocate(jac(nmf,nmf))
            df0 = min(1.0e-05_dp,100*tol) ! amplitude of forcing perturbation for finite difference approximation of gradient
            tol_df = tol
            fpert(1)  = df0
            fpert(2:) = 10*df0
      ! stamp logs
            call nek_log_message('Newton configuration:', this_module, this_procedure)
            write(msg,'(3X,A,1X,E16.8)')     padr('target tol:',     18), tol
            call nek_log_message(msg, this_module, this_procedure)
            write(msg,'(3X,A,A)') padr('tol. scheduling:',18), padl(merge('constant', 'dynamic ', tol_mode_==1), 16)
            call nek_log_message(msg, this_module, this_procedure)
            write(msg,'(3X,A,*(1X,F16.12))') padr('forcing |df|:',   18), fpert
            call nek_log_message(msg, this_module, this_procedure)
            write(msg,'(3X,A,1X,E16.8)')     padr('df newton tol:',  18), tol_df
            call nek_log_message(msg, this_module, this_procedure)
            write(msg,'(3X,A,*(1X,F16.12))') padr('target mflow:',   18), mflow_target
            call nek_log_message(msg, this_module, this_procedure)
            write(msg,'(3X,A,1X,E16.8)')     padr('mflow tol:',      18), tol_mf
            call nek_log_message(msg, this_module, this_procedure)
      ! get reference forcing and phase
            call pipe%get_dpds(dpds, phase)
      !
      ! Baseline Newton iteration on the initial guess of the forcing parameters & baseflow
      !
            call newton_fixed_point_iteration(sys, bf, tol, tol_mode)
      ! extract reference values for mass flow rate, initial mass flow error and flow solution
            call pipe%get_mflow_fft(mflow_old, if_amplitude=.true.)
            allocate(ref); call nek2vec(ref, vx, vy, vz, pr, t)
            call nek_log_message('Reference solution set.', this_module, this_procedure)
            mf_err = mflow_old(:nmf) - mflow_target
		! determine an approximation of the integration error for the mass flux
            call pipe%get_dt_minmax(dt_minmax)
            tol_mf_inexact = (sum(dt_minmax)*0.5)**2/100.0
		      write(msg,'(A,1X,E16.8)') 'Approximate mass flow computation error: ', tol_mf_inexact
            call nek_log_message(msg, this_module, this_procedure)
            if (tol_mf_inexact > tol_mf) call nek_log_warning('Requested mass flow tol is below estimated error.', 
     &            this_module, this_procedure)

      ! stamp logs
            call pipe%parameter_summary()
            call nek_log_information('Initial state:', this_module, this_procedure)
            write(msg,'(A,*(1X,F16.10))') 'mf_state  = ', mflow_old(:nmf)
            call nek_log_information(msg, this_module, this_procedure)
            write(msg,fmt) 'mf_target = ', mflow_target, ' | tol= ', tol_mf
            call nek_log_information(msg, this_module, this_procedure)
            write(msg,fmt) 'mf_error  = ', mf_err, ' | sum= ', sum(abs(mf_err))
            call nek_log_information(msg, this_module, this_procedure)
      !
      ! Main Newton iteration to converge the mass flow rate for all considered Fourier components
      !
            call nek_log_message('Begin mass flow Newton iteration', this_module, this_procedure)
            df_loop: do inwt = 1, maxiter_newton_
               if (sum(abs(mf_err)) < tol_mf) then ! lucky convergence
                  call nek_log_message('The intial condition is a fixed poinf of the system!', this_module, this_procedure)
                  exit df_loop ! converged
               end if
               write(step_id,'("Step ",I3,": ")') inwt
               write(msg,'(A,I0,A)') 'Begin mflow Newton step ', inwt, ' ...'
               call nek_log_information(msg, this_module, this_procedure)
               do i = 1, nmf ! we need to compute dmf/df for all components
                  write(coef_id,'("Fourier coef. ",I2,": ")') i
                  write(msg,'(A,A,I0,A)') step_id, 'compute mflow gradient for Fourier coefficient ', i, ' ...'
                  call nek_log_information(msg, this_module, this_procedure)
                  is_fp = .true.
                  icnt = 0
                  do while (is_fp)
      ! iterate in case the forcing perturbation is too low
                     if (icnt > 0) call nek_log_message('Repeat finite difference step.',this_module, 'mflow_newton')
      ! Set flow parameters
                     dpds_tmp = 0.0_dp
                     fpert(i) = -sign(fpert(i), mf_err(i)) ! take the step in the direction of the root for better accuracy
                     if (i == 1) then
                        dpds_tmp(1) = fpert(1)
                     else
                        j = 2*(i-1)
                        dpds_tmp(j  ) = cos(phase(i))*fpert(i)
                        dpds_tmp(j+1) = sin(phase(i))*fpert(i)
                     end if
                     write(msg,'(A,A,A,*(1X,F16.10))') step_id, coef_id, padl('old frc:',pad), dpds(:nf)
                     call nek_log_message(msg, this_module, this_procedure)
                     write(msg,'(A,A,A,*(1X,F16.10))') step_id, coef_id, padl('prt frc:',pad), dpds_tmp
                     call nek_log_message(msg, this_module, this_procedure)
                     dpds_tmp = dpds_tmp + dpds(:nf)
                     write(msg,'(A,A,A,*(1X,F16.10))') step_id, coef_id, padl('new frc:',pad), dpds_tmp
                     call nek_log_message(msg, this_module, this_procedure)
                     call pipe%init_flow(dpds_tmp, Wo)
                  
      ! reset baseflow
                     call bf%zero(); call bf%add(ref)
      ! Run Newton-Krylov solver to find baseflow of perturbed system
                     call newton_fixed_point_iteration(sys, bf, tol_df, tol_mode, input_is_fixed_point=is_fp)

                     icnt = icnt + 1
      ! increase perturbation amplitude if newton exited without iteration
                     if (is_fp) then
                        fpert(i) = 10*fpert(i)
                        write(msg,'(A,I0,A,F16.10)') 'Perturbation is too small for component ', i, 
     &                        ': Reset |df| = ', fpert(i)
                        call nek_log_message(msg, this_module, this_procedure)
                     end if
                  end do
                  call pipe%get_mflow_fft(mflow_new, if_amplitude=.true.)

      ! get difference, compute gradient for current component and update mass flow Jacobian
                  dmf = mflow_new(:nmf) - mflow_old(:nmf)
                  write(msg,'(A,A,A,*(1X,F16.10))') step_id, coef_id, padl('dmflow:',pad), dmf
                  call nek_log_message(msg, this_module, this_procedure)
                  do j = 1, nmf
                     jac(i,j) = dmf(j)/fpert(i)
                  end do
               end do ! i = 1, nmf
               write(msg,'(A,A)') step_id, 'mflow jacobian'
               call nek_log_message(msg, this_module, this_procedure)
               do i = 1, nmf
                  write(msg,'(A,4X,*(1X,F16.10))') step_id, jac(i,:)
                  call nek_log_message(msg, this_module, this_procedure)
               end do
               if (det(jac) == 0.0_dp) call nek_stop_error('Jacobian is singular!', this_module, this_procedure)
      !
      ! update
      !
               deltaf = -matmul(inv(jac), mf_err)
               write(msg,'(A,A,*(1X,F16.10))') step_id, padl('|f_step|:', pad), deltaf
               call nek_log_message(msg, this_module, this_procedure)
               
      ! Set flow parameters
               dpds(1) = dpds(1) + deltaf(1) ! update reference forcing
               do i = 2, nmf
                  j = 2*(i-1)
                  dpds(j  ) = dpds(j  ) + cos(phase(i))*deltaf(i)
                  dpds(j+1) = dpds(j+1) + sin(phase(i))*deltaf(i)
               end do
      ! reset baseflow
               call bf%zero(); call bf%add(ref)
      ! Update forcing
               call nek_log_message('Forcing prior to Newton step:', module='neklab_helix')
               call pipe%forcing_summary()
               call pipe%init_flow(dpds(:nf), Wo)
      !
      ! Take Newton step for the forcing and compute new periodic orbit
      !
               call newton_fixed_point_iteration(sys, bf, tol, tol_mode)
      ! extract reference values for mass flow rate, initial mass flow error, 
      ! flow solution and forcing amplitudes and phase
               call pipe%get_mflow_fft(mflow_old, if_amplitude=.true.)
               call nek2vec(ref, vx, vy, vz, pr, t)
               call pipe%get_dpds(dpds, phase)
               mf_err = mflow_old(:nmf) - mflow_target
      ! stamp logfile
               call nek_log_information(step_id//'Final state:', this_module, this_procedure)
               write(msg,'(A,*(1X,F16.10))') step_id//'mf_state = ', mflow_old(:nmf)
               call nek_log_information(msg, this_module, this_procedure)
               write(msg,fmt) step_id//'mf_target= ', mflow_target, ' | tol= ', tol_mf
               call nek_log_information(msg, this_module, this_procedure)
               write(msg,fmt) step_id//'mf_error = ', mf_err, ' | sum= ', sum(abs(mf_err))
               call nek_log_information(msg, this_module, this_procedure)
      ! convergence check
               if (sum(abs(mf_err)) < tol_mf) then
		         write(msg,'(A,I0,A)') 'Newton iteration converged after ', inwt, ' iterations.'
                  call nek_log_message(msg, this_module, this_procedure)
                  exit df_loop ! converged
	         else
      ! save intermediate solution
                  call set_fldindex('nwf', 1)
                  call outpost_dnek(bf, 'nwf')
               end if
            end do df_loop
            call nek_log_message('Exiting mass flow Newton iteration.', this_module, this_procedure)
            call set_fldindex('BFN', 1)
            call outpost_dnek(bf, 'BFN')
            if (sum(abs(mf_err)) > tol_mf) then
               write(msg,'(A,I0,A)') 'Mass flux not converged after ', maxiter_newton_, 'steps.'
               call nek_stop_error(msg, this_module, this_procedure)
            end if
      ! stamp logs
            call nek_log_message('OUTPUT:', this_module, this_procedure)
            call pipe%setup_summary()
            call pipe%parameter_summary()
            call pipe%forcing_summary()
            write(msg,'(3X,A,*(1X,F16.12))') 'mf_state    = ', mflow_old(:nmf)
            call nek_log_message(msg, this_module, this_procedure)
            write(msg,'(3X,A,*(1X,F16.12))') 'mf_target   = ', mflow_target
            call nek_log_message(msg, this_module, this_procedure)
            write(msg,'(3X,A,1X,F16.12)') 'mflow error = ', sum(abs(mf_err))
            call nek_log_message(msg, this_module, this_procedure)
            if (size(mflow_old) > nmf) then
               write(long_msg,'(3X,A,*(1X,F16.12))') 'ext. mflow  = ', mflow_old(nmf+1:)
               call nek_log_message(long_msg, this_module, this_procedure)
            end if
         end subroutine mflow_newton

         subroutine shift_mflow_phase_torus(bf, mflow_target)
            type(nek_dvector), intent(inout) :: bf
      !! Current baseflow to be shifted
            real(dp), dimension(:), intent(in) :: mflow_target
      !! Target values for the mass flow rate for each Fourier component
            ! internal
            character(len=*), parameter :: this_procedure = 'shift_mflow_phase'
            integer :: i, nmf
            real(dp) :: phase_dt, Tend
            real(dp), allocatable :: mflow(:), phase(:)
            logical :: save_base_old, save_fft_old
            character(len=128) :: msg
            real(dp), parameter :: tol_dt = 1.0e-04_dp
            nmf = size(mflow_target)
            call nek_log_message('Current forcing:', this_module, this_procedure)
            call pipe%forcing_summary()
            call pipe%get_mflow_fft(mflow, phase, if_amplitude=.true.)
            write(msg,'(A,*(1X,F16.8))') 'mflow phase: ', phase(:nmf)
            call nek_log_message(msg, this_module, this_procedure)
            write(msg,'(A,*(1X,F16.8))') 'phase target:', 0.0_dp*phase(:nmf) 
            call nek_log_message(msg, this_module, this_procedure)
      ! save old logical flags
            save_base_old = pipe%is_save_2d(); call pipe%set_save_base(.false.)
            save_fft_old = pipe%is_save_fft(); call pipe%set_save_fft(.false.)
            do i = 2, nmf            ! the first component is purely real, phase is zero by construction
               phase_dt = phase(i)/pipe%get_omega()
               write(msg,'(A,I0,A,F16.8)') 'Component ', i, ': phase_delta (dt) = ', phase_dt 
               call nek_log_message(msg, this_module, this_procedure)
               if (abs(phase_dt) < tol_dt) then
      ! We adjust dpds without running the solver
                  msg = 'phase delta is very small. Baseflow not shifted but dpds rotated.'
                  call nek_log_message(msg, this_module, this_procedure)
               else
      ! Set the initial condition for nonlinear solve
                  call vec2nek(vx, vy, vz, pr, t, bf)
                  if (phase_dt > 0.0_dp) then ! we need to reduce the phase -> find new initial time and adjust dpds
                     Tend = phase_dt ! Set final time
                  else                             ! we need to increase the phase
                     Tend = pipe%get_period() + phase_dt ! Set final time
                  end if
      ! Set appropriate tolerances and Nek status
                  call setup_nonlinear_solver(variable_dt=.true., endtime=Tend, cfl_limit=0.4_dp)
                  time = 0.0_dp
                  istep = 0
                  do while (lastep == 0)
                     istep = istep + 1
                     call pipe%compute_bf_forcing(time) ! --> set neklab_forcing data
                     call nek_advance()
                  end do
      ! Retrieve baseflow 
                  call nek2vec(bf, vx, vy, vz, pr, t)
               end if
      ! adjust dpds
               call pipe%shift_mflow_phase(i, 0.0_dp)
            end do
            call nek_log_message('Updated forcing:', this_module, this_procedure)
            call pipe%forcing_summary()
      ! reset logical flags
            call pipe%set_save_base(save_base_old)
            call pipe%set_save_fft(save_fft_old)
         end subroutine shift_mflow_phase_torus

         subroutine compute_nonlinear_period_torus(bf_out, bf_in, save_2d, variable_dt, recompute_dt, if_fft, if_res, cfl_limit, tstart)
            type(nek_dvector), intent(out) :: bf_out
      !! Output of the nonlinear solver after a period
            type(nek_dvector), intent(in) :: bf_in
      !! Initial condition for the nonlinear solver
            logical, optional, intent(in) :: save_2d
      !! Save 2D fields?
            logical, optional, intent(in) :: variable_dt
      !! Compute with fixed or variable timestep?
            logical, optional, intent(in) :: recompute_dt
      !! Recompute timestep
            logical, optional, intent(in) :: if_fft
      !! Compute mflow fft?
            logical, optional, intent(in) :: if_res
      !! Compute periodic residual?
            real(dp), optional, intent(in) :: cfl_limit
      !! CFL limit for calculation
            real(dp), optional, intent(in) :: tstart
      !! time setting at start
            ! internal
            character(len=*), parameter :: this_procedure = 'nonlinear_period'
            logical :: get_2d, get_fft, var_dt, get_res, recompute_dt_
            logical :: get_2d_old, get_fft_old, newton_old, floquet_old
            real(dp) :: pd, ubar, rnorm, cfl
            character(len=128) :: msg
            character(len=*), parameter :: fmt = '(3(F16.8,1X),A,F16.8)'
      ! set optional arguments
            get_2d  = optval(save_2d, .false.)
            var_dt  = optval(variable_dt, .true.)     
            get_fft = optval(if_fft, .true.)
            get_res = optval(if_res, .false.)
            cfl     = optval(cfl_limit, 0.5_dp)
            time    = optval(tstart, 0.0_dp)
            recompute_dt_ = optval(recompute_dt, .true.)
      ! save toolbox status
            get_2d_old  = pipe%is_save_2d();  call pipe%set_save_base(get_2d)
            get_fft_old = pipe%is_save_fft(); call pipe%set_save_fft(get_fft)
            newton_old  = pipe%is_newton();   call pipe%set_newton(.false.)  ! in case we save 2d fields
            floquet_old = pipe%is_floquet();  call pipe%set_floquet(.false.) ! in case we save 2d fields
      ! set period
            pd = pipe%get_period()
            if (pd == 0.0_dp) pd = param(10) ! for the steady case
      ! set baseflow intial condition
            call vec2nek(vx, vy, vz, pr, t, bf_in)
      ! set nek status
            call setup_nonlinear_solver(recompute_dt = recompute_dt_,
     $                                  endtime      = pd, 
     $                                  variable_dt  = var_dt,
     $                                  cfl_limit    = cfl)
            call nek_status(full_summary=.true.)
            call pipe%reset_mflow_fft() ! in case we compute the FFT
      ! compute period
            if (var_dt) then ! variable timestep
               istep = 0
               do while (lastep == 0)
                  istep = istep + 1
                  call pipe%compute_bf_forcing(time) ! --> set neklab_forcing data
                  if (get_2d)  call pipe%save_2d_fields(vx,vy,vz) ! outposts automatically at lastep == 1
                  call nek_advance()
                  if (get_fft) call pipe%compute_mflow_fft(period = pd, var_dt = .true.) ! integrate Fourier coefficients
                  ubar = pipe%compute_ubar(vx,vy,vz)
                  write(msg,fmt) time, time/pd, mod(time,pd), 'massflow UBAR: ', ubar
                  call nek_log_information(msg, this_module, this_procedure)
               end do
            else             ! fixed dt
               do istep = 1, nsteps
                  call pipe%compute_bf_forcing(time) ! --> set neklab_forcing data
                  if (get_2d)  call pipe%save_2d_fields(vx,vy,vz) ! outposts automatically at lastep == 1
                  call nek_advance()
                  if (get_fft) call pipe%compute_mflow_fft(period = pd)      ! integrate Fourier coefficients
                  ubar = pipe%compute_ubar(vx,vy,vz)
                  write(msg,fmt) time, time/pd, mod(time,pd), 'massflow UBAR: ', ubar
                  call nek_log_information(msg, this_module, this_procedure)
               end do
            end if
            call pipe%outpost_2d_fields()
            if (get_fft) call pipe%extract_mflow_fft(period = pd)
      ! extract output
            call nek2vec(bf_out, vx, vy, vz, pr, t)
      ! compute periodic residual if requested
            if (get_res) then
               call bf_out%sub(bf_in)
               rnorm = bf_out%norm()
               write(msg,'(A,E15.8)') 'Periodic residual norm= ', rnorm
               call nek_log_message(msg, this_module, this_procedure)
      ! replace output
               call nek2vec(bf_out, vx, vy, vz, pr, t)
            end if
      ! reset logical flags
            call pipe%set_save_base(get_2d_old)
            call pipe%set_save_fft(get_fft_old)
            call pipe%set_newton(newton_old)
            call pipe%set_floquet(floquet_old)
         end subroutine compute_nonlinear_period_torus

         subroutine compute_monodromy_period_torus(pert_out, pert_in, nout, nperiod)
            type(nek_dvector), intent(out) :: pert_out
      !! Output of the linear solver after a period of the monodromy operator
            type(nek_dvector), intent(inout) :: pert_in
      !! Initial condition for the linear solver
            integer, optional, intent(in) :: nout
      !! Number of output fields per period (default = 1)
            integer, optional, intent(in) :: nperiod
      !! Number of periods to compute the linear solution across (default = 1)
		! internal
            character(len=*), parameter :: this_procedure = 'monodromy_period'
            integer :: nout_, nperiod_
            integer :: i, idx, ns, outstep, nsaver
            real(dp) :: norm, gr, Tend, FTLE, tper, tpern, pd
            logical :: existfile
            character(len=128) :: msg
            character(len=132) :: fname
            character(len=*), parameter :: fmt = '(A,1X,I4,1X,A,1X,2(F11.6),1X,A,1X,I2,3(1X,A,1X,E15.8))'
            nout_ = optval(nout, 1)
            nperiod_ = optval(nperiod, 1)
            ns = 0
            idx = 1
            write(fname,'("f2dtorus",I3.3,".fld")') idx
            inquire(file=fname, exist=existfile)
            if (existfile) then
               do while (existfile)
                  ! read first file and get nsteps
                  call pipe%get_nsteps_from_header(fname, nsaver)
                  ns = ns + nsaver              
                  idx = idx + 1
                  write(fname,'("f2dtorus",I3.3,".fld")') idx
                  inquire(file=fname, exist=existfile)
               end do
               call bcast(ns, isize)          ! broadcast number of saved snapshots
               call pipe%set_nsteps(ns)
               write(msg,'(A,I0,A,I0,A)') 'Found ', idx-1, ' baseflow files: ', ns, ' timesteps per period.'
               call nek_log_message(msg, this_module, this_procedure)
            else
               msg = "No 2d baseflow files in the format f2dtorus???.fld found. Abort."
               call nek_stop_error(msg, this_module, this_procedure)
            end if

            ! run a period to get GR and FTLE data
            ns = pipe%get_nsteps()
            pd = pipe%get_period()
            outstep = floor(1.0*ns/nout_)

            ! Initial perturbation
            norm = pert_in%norm()
            call pert_in%scal(1.0/norm)
            norm = pert_in%norm()
            call vec2nek(vxp, vyp, vzp, prp, tp, pert_in)

            do i = 1, nperiod_
               Tend = i*pd
               call setup_linear_solver(variable_dt = .true.,
     &                                  endtime     = Tend,
     &                                  cfl_limit   = 0.5_dp)

               time = (i-1)*pd

               tper = 0.0_dp
               FTLE = 0.0_dp
               call pipe%set_2d_mode('floquet') ! reset output counter to load baseflow files in order
               do istep = 1, ns
               ! update baseflow
               call pipe%set_baseflow(vx, vy, vz, istep)
               
                        ! compute linear step
               call nek_advance()
               
                        ! compute growth rate
               call nek2vec(pert_in, vxp, vyp, vzp, prp, tp)
               gr = (pert_in%norm() - norm)/(norm*dt)
               
                        ! FTLE
               FTLE = FTLE + gr*dt
               tper = tper + dt
                        tpern = tper/pd
               if (io_rank() .and. .not. istep == ns) then
                  write(msg,fmt) 'istep', istep, 't', time, tper, tpern, 'P', i,
     &                          'norm', norm, 'gr', gr, 'FTLE', FTLE/tper
                  call nek_log_message(msg, this_module, this_procedure)
               end if

               ! update norm
               norm = pert_in%norm()

               ! outpost
               if (istep == 1 .or. mod(istep,outstep) == 0) then
                  call outpost_dnek(pert_in, 'prt')
               end if
               end do
               if (io_rank()) then
                  write(msg,'(A,I3,A,E15.8)') 'Period ', i, ' FTLE ', FTLE/tper
                  call nek_log_message(msg, this_module, this_procedure)
               end if
            end do

         end subroutine compute_monodromy_period_torus
            
         subroutine compute_energy_budgets_period_torus(pert_in, nout)
            type(nek_dvector), intent(inout) :: pert_in
      !! Initial condition for the linear solver
            integer, optional, intent(in) :: nout
      !! Number of output fields per period (default = 1)
		! internal
            character(len=*), parameter :: this_procedure = 'energy_budgets_period'
            integer :: nout_, nperiod_
            integer :: i, j, k, e, ijke, lv
            integer :: idx, ns, outstep, nsaver
            real(dp) :: norm, gr, FTLE, tper, pd, t_
            real(dp), dimension(lx1,ly1,lz1,lelv) :: dUdx_g, dUdy_g, dUdz_g
            real(dp), dimension(lx1,ly1,lz1,lelv) :: dVdx_g, dVdy_g, dVdz_g
            real(dp), dimension(lx1,ly1,lz1,lelv) :: dWdx_g, dWdy_g, dWdz_g
            real(dp), dimension(lx1,ly1,lz1,lelv) :: prod_g, diss_g, tmpfld
            real(dp) :: dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz
            real(dp) :: ux, uy, uz
            real(dp) :: prod_i, diss_i
            real(dp) :: prod, diss
            type(helix), allocatable :: prtpipe
            logical :: existfile
            character(len=256) :: msg
            character(len=132) :: fname
            character(len=*), parameter :: fmt = '(A,1X,I4,1X,A,1X,2(F11.6),6(1X,A,1X,E15.8))'
            nout_ = optval(nout, 1)
            ns = 0
            idx = 1
            write(fname,'("f2dtorus",I3.3,".fld")') idx
            inquire(file=fname, exist=existfile)
            if (existfile) then
               do while (existfile)
                  ! read first file and get nsteps
                  call pipe%get_nsteps_from_header(fname, nsaver)
                  ns = ns + nsaver              
                  idx = idx + 1
                  write(fname,'("f2dtorus",I3.3,".fld")') idx
                  inquire(file=fname, exist=existfile)
               end do
               call bcast(ns, isize)          ! broadcast number of saved snapshots
               call pipe%set_nsteps(ns)
               write(msg,'(A,I0,A,I0,A)') 'Found ', idx-1, ' baseflow files: ', ns, ' timesteps per period.'
               call nek_log_message(msg, this_module, this_procedure)
            else
               msg = "No 2d baseflow files in the format f2dtorus???.fld found. Abort."
               call nek_stop_error(msg, this_module, this_procedure)
            end if

            ! run a period to get GR and FTLE data
            ns = pipe%get_nsteps()
            pd = pipe%get_period()
            outstep = floor(1.0*ns/nout_)

            ! Initial perturbation
            norm = pert_in%norm()
            call pert_in%scal(1.0/norm)
            norm = pert_in%norm()
            call vec2nek(vxp, vyp, vzp, prp, tp, pert_in)

            time = 0.0_dp
            t_ = 0.0_dp ! lagged time
            call pipe%set_2d_mode('floquet') ! reset output counter to load baseflow files in order
            call setup_linear_solver(variable_dt = .true.,
     &                               endtime     = pd,
     &                               cfl_limit   = 0.5_dp)

            allocate(prtpipe, source=pipe)
            call prtpipe%set_newton(.false.)  ! in case we save 2d fields
            call prtpipe%set_floquet(.false.) ! in case we save 2d fields
            lv = lx1*ly1*lz1*nelv
            call rzero(prod_g, lv); prod = 0.0_dp
            call rzero(diss_g, lv); diss = 0.0_dp
            do istep = 1, ns
               ! update baseflow
               call pipe%set_baseflow(vx, vy, vz, istep)
               
               ! compute linear step
               call nek_advance()
               
               ! compute growth rate
               call nek2vec(pert_in, vxp, vyp, vzp, prp, tp)

               ! compute baseflow gradients
               call gradm1(dUdx_g, dUdy_g, dUdz_g, vx)
               call gradm1(dVdx_g, dVdy_g, dVdz_g, vy)
               call gradm1(dWdx_g, dWdy_g, dWdz_g, vz)

               ! compute turbulent kinetic energy production
               call rzero(tmpfld, lv); prod_i = 0.0_dp
               do e = 1, nelv
                  do i = 1, lz1
                     do j = 1, ly1
                        do k = 1, lx1
                           ijke = i+lx1*((j-1)+ly1*((k-1) + lz1*(e-1)))
                           ux = vxp(ijke,1)
                           uy = vyp(ijke,1)
                           uz = vzp(ijke,1)
                           dUdx = dUdx_g(i,j,k,e)
                           dUdy = dUdy_g(i,j,k,e)
                           dUdz = dUdz_g(i,j,k,e)
                           dVdx = dVdx_g(i,j,k,e)
                           dVdy = dVdy_g(i,j,k,e)
                           dVdz = dVdz_g(i,j,k,e)
                           dWdx = dWdx_g(i,j,k,e)
                           dWdy = dWdy_g(i,j,k,e)
                           dWdz = dWdz_g(i,j,k,e)
                           tmpfld(i, j, k, e) = - ux*(ux*dUdx + uy*dUdy + uz*dUdz)
     &                                          - uy*(ux*dVdx + uy*dVdy + uz*dVdz)
     &                                          - uz*(ux*dWdx + uy*dWdy + uz*dWdz)
                           prod_i = tmpfld(i, j, k, e)*bm1(i, j, k, e)
                           prod_g(i, j, k, e) = (prod_g(i, j, k, e)*t_ + prod_i*dt)/time
                        end do
                     end do
                  end do
               end do
               
               ! compute perturbation gradients
               call gradm1(dUdx_g, dUdy_g, dUdz_g, vxp)
               call gradm1(dVdx_g, dVdy_g, dVdz_g, vyp)
               call gradm1(dWdx_g, dWdy_g, dWdz_g, vzp)

               ! compute turbulent kinetic energy dissipation
               call rzero(tmpfld, lx1*ly1*lz1*lelv); diss_i = 0.0_dp
               do e = 1, nelv
                  do i = 1, lz1
                     do j = 1, ly1
                        do k = 1, lx1
                           ijke = i+lx1*((j-1)+ly1*((k-1) + lz1*(e-1)))
                           dUdx = dUdx_g(i,j,k,e)
                           dUdy = dUdy_g(i,j,k,e)
                           dUdz = dUdz_g(i,j,k,e)
                           dVdx = dVdx_g(i,j,k,e)
                           dVdy = dVdy_g(i,j,k,e)
                           dVdz = dVdz_g(i,j,k,e)
                           dWdx = dWdx_g(i,j,k,e)
                           dWdy = dWdy_g(i,j,k,e)
                           dWdz = dWdz_g(i,j,k,e)
                           tmpfld(i, j, k, e) = vdiff(i, j, k, e,1)*(
     &                                        + dUdx**2 + dUdy**2 + dUdz**2
     &                                        + dVdx**2 + dVdy**2 + dVdz**2
     &                                        + dWdx**2 + dWdy**2 + dWdz**2)
                           diss_i = tmpfld(i, j, k, e)*bm1(i, j, k, e)
                           diss_g(i, j, k, e) = (diss_g(i, j, k, e)*t_ + diss_i*dt)/time
                        end do
                     end do
                  end do
               end do

               call prtpipe%save_2d_fields(prod_g, diss_g, vxp)

               gr = (pert_in%norm() - norm)/(norm*dt)
               
               ! FTLE
               FTLE = FTLE + gr*dt
               tper = time/pd
               if (io_rank() .and. .not. istep == ns) then
                  write(msg,fmt) 'istep', istep, 't', time, tper,
     &                           'norm', norm, 'E', norm**2,
     &                           'gr', gr, 'FTLE', FTLE/tper,
     &                           'P', prod_i, 'D', diss_i
                  call nek_log_message(msg, this_module, this_procedure)
               end if

               ! update norm and lagged time
               norm = pert_in%norm()
               t_ = time
            end do
            if (io_rank()) then
               write(msg,'(A,E15.8)') 'FTLE ', FTLE/tper
               call nek_log_message(msg, this_module, this_procedure)
               call outpost(prod_g, diss_g, vx, pr, t, 'bdg')
            end if
         end subroutine compute_energy_budgets_period_torus
      
         end module neklab_analysis_torus
