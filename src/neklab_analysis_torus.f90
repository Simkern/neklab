      module neklab_analysis_torus
         use stdlib_stats_distribution_normal, only: normal => rvs_normal
         use stdlib_strings, only: padl
         use stdlib_optval, only: optval
         use stdlib_linalg, only: diag, eye, det, inv
         use stdlib_logger, only: information_level, warning_level, debug_level, error_level, all_level, success
         use LightKrylov, only: atol_dp, dp, eigs, svds, save_eigenspectrum
         use LightKrylov, only: kexpm, gmres_rdp
         use LightKrylov, only: initialize_krylov_subspace, orthonormalize_basis, zero_basis, rand_basis
         use LightKrylov, only: linear_combination, innerprod
         use LightKrylov, only: newton, newton_dp_opts
         use LightKrylov_Logger
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
      
      contains

         subroutine linear_stability_analysis_periodic_orbit(floquet_operator, kdim, nev, adjoint, X0)
            type(floquet_linop), intent(inout) :: floquet_operator
      !! Floquet perator whose stability properties are to be investigated.
            integer, intent(in) :: kdim
      !! Maximum dimension of the Krylov subspace.
            integer, intent(in) :: nev
      !! Desired number of eigenpairs to converge.
            logical, intent(in), optional :: adjoint
      !! Whether direct or adjoint analysis should be conducted.
		type(nek_dvector), optional, intent(in) :: X0
	!! Initial guess for the eigenvectors
      
      ! Eigenvalue computation related variables.
            type(nek_dvector), allocatable :: eigvecs(:)
            complex(kind=dp), allocatable :: eigvals(:)
            real(kind=dp), allocatable :: residuals(:)
            integer :: info
      
      ! Miscellaneous.
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
		call eigs(floquet_operator, eigvecs, eigvals, residuals, info, x0=X0, kdim=kdim, transpose=adjoint_)
      
      ! Transform eigenspectrum to continuous-time representation.
            eigvals = log(eigvals)/floquet_operator%tau
      
      ! Determine the file prefix.
            file_prefix = merge("adj", "dir", adjoint_)
      
      ! Save eigenspectrum to disk.
            call save_eigenspectrum(eigvals, residuals, trim(file_prefix)//"_eigenspectrum.npy")
      
      ! Export eigenfunctions to disk.
            call outpost_dnek(eigvecs(:nev), file_prefix)

		call logger%log_message('Exiting eigenvalue computation.', module=this_module)

      ! Finalize exptA timings
            call floquet_operator%finalize_timer()
      ! Finalize timing
            call logger_setup(logfile='lightkrylov_tmr.log', nio=0, log_level=warning_level, log_stdout=.false., log_timestamp=.true.)
            call timer%finalize()
      
            return
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
            type(nek_dvector) :: ref
            logical :: is_new_solution
            integer :: tol_mode_, maxiter_newton_
            integer :: nmf, nf, inwt, i, j
            real(dp) :: Wo, df0
            real(dp) :: dpds(lf)
            real(dp), allocatable :: phase(:), dpds_tmp(:)
            real(dp), allocatable :: mflow_old(:), mflow_new(:)
            real(dp), allocatable :: dmf(:), mf_err(:), deltaf(:), fpert(:)
            real(dp), allocatable :: jac(:,:)
				real(dp) :: dt_minmax(2)
            real(dp) :: tol_mf_inexact, tol_df
            character(len=128) :: msg, fmt
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
            nf = pipe%get_nf()      ! number of real forcing components (real and imaginary parts counted individually)
            if (nmf > nfft) then
               write(msg,'(A,I0,A,I0,A)') 'nmf= ', nmf, ' > nfft= ', nfft,'. Increase nfft in neklab_helix.'
               call nek_stop_error(msg, module=this_module, procedure='mflow_newton')
            else if (2*nmf - 1 /= nf) then
               write(msg,'(A,I0,A,I0,A)') 'nmf= ', nmf, ' and nf= ', nf, ' incompatible.'
               call nek_stop_error(msg, module=this_module, procedure='mflow_newton')
            else
               write(msg,'(A,I0,A)') 'Starting Newton iteration for nmf = ', nmf, ' mass flow components.'
               call nek_log_information(msg, module=this_module, procedure='mflow_newton')
            end if
            allocate(dpds_tmp(nf))
            allocate(dmf(nmf), mf_err(nmf), deltaf(nmf), fpert(nmf))
            allocate(jac(nmf,nmf))
            df0 = min(1.0e-06_dp,100*tol) ! amplitude of forcing perturbation for finite difference approximation of gradient
            tol_df = tol
            fpert(1)  = df0
            fpert(2:) = 20*df0
            ! get reference forcing and phase
            call pipe%get_dpds(dpds, phase)
      !
      ! Baseline Newton iteration on the initial guess of the forcing parameters & baseflow
      !
            call newton_fixed_point_iteration(sys, bf, tol, tol_mode)
      ! extract reference values for mass flow rate, initial mass flow error and flow solution
            call pipe%get_mflow_fft(mflow_old, if_amplitude=.true.)
            call nek2vec(ref, vx, vy, vz, pr, t)
            call nek_log_message('Reference solution set.', module=this_module, procedure='mflow_newton')
            mf_err = mflow_old(:nmf) - mflow_target
		! determine an approximation of the integration error for the mass flux
				call pipe%get_dt_minmax(dt_minmax)
				tol_mf_inexact = (sum(dt_minmax)*0.5)**2/100.0
				write(msg,'(A,1X,E16.8)') 'Approximate mass flow computation error: ', tol_mf_inexact
            call nek_log_message(msg, module=this_module, procedure='mflow_newton')
            if (tol_mf_inexact > tol_mf) then
               write(msg,'(A,1X,E16.8)') 'Reset tolerance for mass flow to tol= ', tol_mf_inexact
               call nek_log_message(msg, module=this_module, procedure='mflow_newton')
               tol_mf = tol_mf_inexact
            end if

      ! stamp logs
            call pipe%parameter_summary()
            write(msg,'(A,*(1X,F16.10))') padl('|df|:',   pad), fpert
            call nek_log_information('Initial state:', module=this_module, procedure='mflow_newton')
            write(msg,'(A,*(1X,F16.10))') 'mf_state = ', mflow_old(:nmf)
            call nek_log_information(msg, module=this_module, procedure='mflow_newton')
            write(msg,fmt) 'mf_target= ', mflow_target, ' | tol= ', tol_mf
            call nek_log_information(msg, module=this_module, procedure='mflow_newton')
            write(msg,fmt) 'mf_error = ', mf_err, ' | sum= ', sum(abs(mf_err))
            call nek_log_information(msg, module=this_module, procedure='mflow_newton')
            !
            ! Main Newton iteration to converge the mass flow rate for each Fourier component
            !
            call nek_log_message('Begin mass flow Newton iteration', module=this_module, procedure='mflow_newton')
            df_loop: do inwt = 1, maxiter_newton_
               if (sum(abs(mf_err)) < tol_mf) then ! lucky convergence
                  call nek_log_message('The intial condition is a fixed poinf of the system!', module=this_module, procedure='mflow_newton')
                  exit df_loop ! converged
               end if
               write(step_id,'("Step ",I3,": ")') inwt
               write(msg,'(A,I0,A)') 'Begin mflow Newton step ', inwt, ' ...'
               call nek_log_information(msg, module=this_module, procedure='mflow_newton')
               do i = 1, nmf ! we need to compute dmf/df for all components
                  write(coef_id,'("Fourier coef. ",I2,": ")') i
                  write(msg,'(A,A,I0,A)') step_id, 'compute mflow gradient for Fourier coefficient ', i, ' ...'
                  call nek_log_information(msg, module=this_module, procedure='mflow_newton')
                  ! iterate in case the forcing is too low the first time around
                  is_new_solution = .false.
                  do while (.not. is_new_solution)
                     ! Set flow parameters
                     dpds_tmp = 0.0_dp
                     fpert(i) = -sign(fpert(i), mf_err(i)) ! take the step in the direction of the root
                     if (i == 1) then
                        dpds_tmp(1) = fpert(1)
                     else
                        j = 2*(i-1)
                        dpds_tmp(j  ) = cos(phase(i))*fpert(i)
                        dpds_tmp(j+1) = sin(phase(i))*fpert(i)
                     end if
                     write(msg,'(A,A,A,*(1X,F16.10))') step_id, coef_id, padl('old frc:',pad), dpds(:nf)
                     call nek_log_message(msg, module=this_module, procedure='mflow_newton')
                     write(msg,'(A,A,A,*(1X,F16.10))') step_id, coef_id, padl('prt frc:',pad), dpds_tmp
                     call nek_log_message(msg, module=this_module, procedure='mflow_newton')
                     dpds_tmp = dpds_tmp + dpds(:nf)
                     write(msg,'(A,A,A,*(1X,F16.10))') step_id, coef_id, padl('new frc:',pad), dpds_tmp
                     call nek_log_message(msg, module=this_module, procedure='mflow_newton')
                     call pipe%init_flow(dpds_tmp, Wo)
                  
                     ! reset baseflow
                     call bf%zero(); call bf%add(ref)
                     ! Run Newton-Krylov solver to find baseflow of perturbed system
                     call newton_fixed_point_iteration(sys, bf, tol_df, tol_mode, is_new_solution=is_new_solution)

                     ! increase perturbation amplitude if newton exited without iteration
                     if (.not. is_new_solution) then
                        fpert(i) = 10*fpert(i)
                        write(msg,'(A,I0,A,F16.10)') 'Perturbation is too small for component ', i, ': Reset |df| = ', fpert(i)
                        call nek_log_message(msg, module=this_module, procedure='mflow_newton')
                     end if
                  end do
                  call pipe%get_mflow_fft(mflow_new, if_amplitude=.true.)

                  ! get difference and compute gradient
                  dmf = mflow_new(:nmf) - mflow_old(:nmf)
			         write(msg,'(A,A,A,*(1X,F16.10))') step_id, coef_id, padl('dmflow:',pad), dmf
			         call nek_log_message(msg, module=this_module, procedure='mflow_newton')
                  do j = 1, nmf
                     jac(i,j) = dmf(j)/fpert(i)
                  end do
               end do ! i = 1, nmf
               write(msg,'(A,A)') step_id, 'mflow jacobian'
               call nek_log_message(msg, module=this_module, procedure='mflow_newton')
               do i = 1, nmf
                  write(msg,'(A,4X,*(1X,F16.10))') step_id, jac(i,:)
                  call nek_log_message(msg, module=this_module, procedure='mflow_newton')
               end do
               if (det(jac) == 0.0_dp) call nek_stop_error('Jacobian is singular!', module=this_module, procedure='mflow_newton')
               !
               ! update
               !
               deltaf = -matmul(inv(jac), mf_err)
               write(msg,'(A,A,*(1X,F16.10))') step_id, padl('|f_step|:', pad), deltaf
               call nek_log_message(msg, module=this_module, procedure='mflow_newton')
               
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
               ! Take Newton step for the forcing and compute periodic orbit
               !
               call newton_fixed_point_iteration(sys, bf, tol, tol_mode)
               ! extract reference values for mass flow rate, initial mass flow error, flow solution and forcing amplitudes and phase
               call pipe%get_mflow_fft(mflow_old, if_amplitude=.true.)
               call nek2vec(ref, vx, vy, vz, pr, t)
               call pipe%get_dpds(dpds, phase)
               mf_err = mflow_old(:nmf) - mflow_target
               ! stamp logfile
               call nek_log_information(step_id//'Final state:', module=this_module, procedure='mflow_newton')
               write(msg,'(A,*(1X,F16.10))') step_id//'mf_state = ', mflow_old(:nmf)
               call nek_log_information(msg, module=this_module, procedure='mflow_newton')
               write(msg,fmt) step_id//'mf_target= ', mflow_target, ' | tol= ', tol_mf
               call nek_log_information(msg, module=this_module, procedure='mflow_newton')
               write(msg,fmt) step_id//'mf_error = ', mf_err, ' | sum= ', sum(abs(mf_err))
               call nek_log_information(msg, module=this_module, procedure='mflow_newton')
               ! convergence check
               if (sum(abs(mf_err)) < tol_mf) then
		         write(msg,'(A,I0,A)') 'Newton iteration converged after ', inwt, ' iterations.'
                  call nek_log_message(msg, module=this_module, procedure='mflow_newton')
                  exit df_loop ! converged
					else
                  ! save intermediate solution
                  call set_fldindex('nwf', 1)
                  call outpost_dnek(bf, 'nwf')
               end if
            end do df_loop
            call nek_log_message('Exiting mass flow Newton iteration.', module=this_module, procedure='mflow_newton')
            call set_fldindex('BFN', 1)
            call outpost_dnek(bf, 'BFN')
            if (sum(abs(mf_err)) > tol_mf) then
               write(msg,'(A,I0,A)') 'Mass flux not converged after ', maxiter_newton_, 'steps.'
               call nek_stop_error(msg, module=this_module, procedure='mflow_newton')
            end if
         end subroutine mflow_newton

         subroutine shift_mflow_phase_torus(bf, nmf)
            type(nek_dvector), intent(inout) :: bf
      !! Current baseflow to be shifted
            integer, intent(in) :: nmf
      !! number of mass flow components to consider
            ! internal
            integer :: i
            real(dp) :: phase_dt, Tend
            real(dp), allocatable :: mflow(:), phase(:)
            logical :: save_base_old, save_fft_old
            character(len=128) :: msg
            real(dp), parameter :: tol_dt = 1.0e-04_dp
            call nek_log_message('Current forcing:', this_module, 'shift_mflow_phase_torus')
            call pipe%forcing_summary()
            call pipe%get_mflow_fft(mflow, phase, if_amplitude=.true.)
            write(msg,'(A,*(1X,F16.8))') 'mflow phase: ', phase(:nmf)
            call nek_log_message(msg, this_module, 'shift_mflow_phase_torus')
            write(msg,'(A,*(1X,F16.8))') 'phase target:', 0.0_dp*phase(:nmf) 
            call nek_log_message(msg, this_module, 'shift_mflow_phase_torus')
            ! save old logical flags
            save_base_old = pipe%is_save_2d(); call pipe%set_save_base(.false.)
            save_fft_old = pipe%is_save_fft(); call pipe%set_save_fft(.false.)
            do i = 2, nmf            ! the first component is purely real, phase is zero by construction
               phase_dt = phase(i)/pipe%get_omega()
               write(msg,'(A,I0,A,F16.8)') 'Component ', i, ': phase_delta (dt) = ', phase_dt 
               call nek_log_message(msg, this_module, 'shift_mflow_phase_torus')
               if (abs(phase_dt) < tol_dt) then
      ! We adjust dpds without running the solver
                  msg = 'phase delta is very small. Baseflow not shifted but dpds rotated.'
                  call nek_log_message(msg, this_module, 'shift_mflow_phase_torus')
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
            call nek_log_message('Updated forcing:', this_module, 'shift_mflow_phase_torus')
            call pipe%forcing_summary()
            ! reset logical flags
            call pipe%set_save_base(save_base_old)
            call pipe%set_save_fft(save_fft_old)
         end subroutine shift_mflow_phase_torus

         subroutine compute_nonlinear_period_torus(bf_out, bf_in, save_2d, variable_dt, if_fft, if_res, cfl_limit, tstart)
            type(nek_dvector), intent(out) :: bf_out
      !! Output of the nonlinear solver after a period
            type(nek_dvector), intent(in) :: bf_in
      !! Initial condition for the nonlinear solver
            logical, optional, intent(in) :: save_2d
      !! Save 2D fields?
            logical, optional, intent(in) :: variable_dt
      !! Compute with fixed or variable timestep?
            logical, optional, intent(in) :: if_fft
      !! Compute mflow fft?
            logical, optional, intent(in) :: if_res
      !! Compute periodic residual?
            real(dp), optional, intent(in) :: cfl_limit
      !! CFL limit for calculation
            real(dp), optional, intent(in) :: tstart
      !! time setting at start
            ! internal
            logical :: get_2d, get_fft, var_dt, get_res
            logical :: get_2d_old, get_fft_old
            real(dp) :: pd, ubar, rnorm, cfl
            character(len=128) :: msg
      ! set optional arguments
            get_2d  = optval(save_2d, .false.)
            var_dt  = optval(variable_dt, .true.)     
            get_fft = optval(if_fft, .true.)
            get_res = optval(if_res, .false.)
            cfl     = optval(cfl_limit, 0.5_dp)
            time    = optval(tstart, 0.0_dp)
      ! set period
            pd = pipe%get_period()
            if (pd == 0.0_dp) pd = param(10) ! for the steady case
      ! set baseflow intial condition
            call vec2nek(vx, vy, vz, pr, t, bf_in)
      ! set nek status
            call setup_nonlinear_solver(recompute_dt = .true.,
     $                                  endtime      = pd, 
     $                                  variable_dt  = var_dt,
     $                                  cfl_limit    = cfl)
            call nek_status(full_summary=.true.)
            call pipe%reset_mflow_fft() ! in case we compute the FFT
            get_2d_old  = pipe%is_save_2d();  call pipe%set_save_base(get_2d)
            get_fft_old = pipe%is_save_fft(); call pipe%set_save_fft(get_fft)
      ! compute period
            if (var_dt) then ! variable timestep
               istep = 0
               do while (lastep == 0)
                  istep = istep + 1
                  call pipe%compute_bf_forcing(time) ! --> set neklab_forcing data
                  call nek_advance()
                  if (get_2d)  call pipe%save_2d_fields(vx,vy,vz) ! outposts automatically at lastep == 1
                  if (get_fft) call pipe%compute_mflow_fft(period = pd, var_dt = .true.) ! integrate Fourier coefficients
                  ubar = pipe%compute_ubar(vx,vy,vz)
                  write(msg,'(3(F16.8,1X),A,F16.8)') time, time/pd, mod(time,pd), 'massflow UBAR: ', ubar
                  call nek_log_information(msg, this_module, 'compute_nonlinear_period')
               end do
            else             ! fixed dt
               do istep = 1, nsteps
                  call pipe%compute_bf_forcing(time) ! --> set neklab_forcing data
                  call nek_advance()
                  if (get_2d)  call pipe%save_2d_fields(vx,vy,vz) ! outposts automatically at lastep == 1
                  if (get_fft) call pipe%compute_mflow_fft(period = pd)      ! integrate Fourier coefficients
                  ubar = pipe%compute_ubar(vx,vy,vz)
                  write(msg,'(3(F16.8,1X),A,F16.8)') time, time/pd, mod(time,pd), 'massflow UBAR: ', ubar
                  call nek_log_information(msg, this_module, 'compute_nonlinear_period')
               end do
            end if
            if (get_fft) call pipe%extract_mflow_fft(period = pd)
      ! extract output
            call nek2vec(bf_out, vx, vy, vz, pr, t)
      ! compute periodic residual if requested
            if (get_res) then
               call bf_out%sub(bf_in)
               rnorm = bf_out%norm()
               write(msg,'(A,E15.8)') 'Periodic residual norm= ', rnorm
               call nek_log_message(msg, this_module, 'compute_nonlinear_period')
      ! replace output
               call nek2vec(bf_out, vx, vy, vz, pr, t)
            end if
      ! reset logical flags
            call pipe%set_save_base(get_2d_old)
            call pipe%set_save_fft(get_fft_old)
         end subroutine compute_nonlinear_period_torus
      
         end module neklab_analysis_torus
