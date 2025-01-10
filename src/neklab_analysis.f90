      module neklab_analysis
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
         use neklab_vectors
         use neklab_linops
         use neklab_utils
         use neklab_nek_setup
         use neklab_otd
         use neklab_systems
         use neklab_helix
      
         implicit none
         include "SIZE"
         include "TOTAL"
         include "ADJOINT"
      
         private
         character(len=*), parameter, private :: this_module = 'neklab_analysis'
      
         public :: linear_stability_analysis_fixed_point
         public :: transient_growth_analysis_fixed_point
         public :: newton_fixed_point_iteration
         public :: newton_periodic_orbit
         public :: newton_forced_periodic_orbit_torus
         public :: mflow_newton_periodic_orbit_torus
         public :: otd_analysis
      
      contains
      
         subroutine linear_stability_analysis_fixed_point(exptA, kdim, nev, adjoint, X0)
            type(exptA_linop), intent(inout) :: exptA
      !! Operator whose stability properties are to be investigated.
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
		      call eigs(exptA, eigvecs, eigvals, residuals, info, x0=X0, kdim=kdim, transpose=adjoint_)
      
      ! Transform eigenspectrum to continuous-time representation.
            eigvals = log(eigvals)/exptA%tau
      
      ! Determine the file prefix.
            file_prefix = merge("adj", "dir", adjoint_)
      
      ! Save eigenspectrum to disk.
            call save_eigenspectrum(eigvals, residuals, trim(file_prefix)//"_eigenspectrum.npy")
      
      ! Export eigenfunctions to disk.
            call outpost_dnek(eigvecs(:nev), file_prefix)

		      call logger%log_message('Exiting eigenvalue computation.', module=this_module)

      ! Finalize exptA timings
            call exptA%finalize_timer()
      ! Finalize timing
            call logger_setup(logfile='lightkrylov_tmr.log', nio=0, log_level=warning_level, log_stdout=.false., log_timestamp=.true.)
            call timer%finalize()
      
            return
         end subroutine linear_stability_analysis_fixed_point
      
         subroutine transient_growth_analysis_fixed_point(exptA, nsv, kdim)
            type(exptA_linop), intent(inout) :: exptA
      !! Operator whose singular value decomposition needs to be computed.
            integer, intent(in) :: nsv
      !! Desired number of singular triplets.
            integer, intent(in) :: kdim
      !! Maximum dimension of the Krylov subspace in LightKrylov.
      
      ! Singular value decomposition.
            type(nek_dvector), allocatable :: U(:), V(:)
            real(kind=dp), allocatable :: S(:), residuals(:)
            integer :: info
      
      ! Miscellaneous.
            integer :: i, j
            character(len=3) :: file_prefix
      
      ! Set up logging
            call logger_setup(nio=0, log_level=information_level, log_stdout=.false., log_timestamp=.true.)
      
      ! Allocate singular vectors.
            allocate (U(nsv)); call initialize_krylov_subspace(U)
            allocate (V(nsv)); call initialize_krylov_subspace(V)
      
      ! Call to LightKrylov.
            call svds(exptA, U, S, V, residuals, info, kdim=kdim)
      
      ! Save singular spectrum to disk.
            if (nid == 0) then
               open (unit=1234, file="singular_spectrum.dat")
               write (1234, *) S
               close (1234)
            end if
      
      ! Export optimal perturbations and optimal responses.
            file_prefix = "prt"; call outpost_dnek(V(:nsv), file_prefix)
            file_prefix = "rsp"; call outpost_dnek(U(:nsv), file_prefix)
      
            return
         end subroutine transient_growth_analysis_fixed_point
      
         subroutine newton_fixed_point_iteration(sys, bf, tol, tol_mode)
            type(nek_system), intent(inout) :: sys
      !! System for which a fixed point is sought
            type(nek_dvector), intent(inout) :: bf
      !! Initial guess for the fixed point
            real(dp), intent(inout) :: tol
      !! Absolute tolerance for the Newton solver
            integer, optional, intent(in) :: tol_mode
      
      ! Misc
            integer :: info, tol_mode_
            type(newton_dp_opts) :: opts
            character(len=3) :: file_prefix
      
		      tol_mode_ = optval(tol_mode, 1)

      	   call logger%log_message('Starting newton iteration.', module=this_module)
      
      ! Define options for the Newton solver
            opts = newton_dp_opts(maxiter=40, ifbisect=.true.)
      
      ! Call to LightKrylov.
            if (tol_mode_ == 1) then
               call newton(sys, bf, gmres_rdp, info, atol=tol, options=opts, scheduler=nek_constant_tol)
            else
		         call newton(sys, bf, gmres_rdp, info, atol=tol, options=opts, scheduler=nek_dynamic_tol)
		      end if
      
      ! Outpost initial condition.
            file_prefix = 'nwt'
            call outpost_dnek(bf, file_prefix)

		      call logger%log_message('Exiting newton iteration.', module=this_module)
      
            return
         end subroutine newton_fixed_point_iteration
      
         subroutine newton_periodic_orbit(sys, bf, tol, tol_mode)
            type(nek_system_upo), intent(inout) :: sys
      !! System for which a fixed point is sought
            type(nek_ext_dvector), intent(inout) :: bf
      !! Initial guess for the fixed point
            real(dp), intent(inout) :: tol
      !! Absolute tolerance for the Newton solver
            integer, optional, intent(in) :: tol_mode
      
      ! Misc
            integer :: info, tol_mode_
            type(newton_dp_opts) :: opts
      !type(gmres_dp_opts)  :: gmres_opts
            character(len=3) :: file_prefix

            ! Set up logging
            call logger%log_message('Starting newton iteration.', module=this_module)
      
      ! Define options for the Newton solver
            opts = newton_dp_opts(maxiter=40, ifbisect=.true.)
      
      ! Call to LightKrylov.
            if (tol_mode_ == 1) then
               call newton(sys, bf, gmres_rdp, info, atol=tol, options=opts, scheduler=nek_constant_tol)
            else
		         call newton(sys, bf, gmres_rdp, info, atol=tol, options=opts, scheduler=nek_dynamic_tol)
		      end if
      
      ! Outpost initial condition.
            file_prefix = 'nwt'
            call outpost_ext_dnek(bf, file_prefix)

            call logger%log_message('Exiting newton iteration.', module=this_module)
      
            return
         end subroutine newton_periodic_orbit

         subroutine newton_forced_periodic_orbit_torus(sys, bf, tol, tol_mode)
            type(nek_system_torus), intent(inout) :: sys
      !! System for which a fixed point is sought
            type(nek_dvector), intent(inout) :: bf
      !! Initial guess for the fixed point
            real(dp), intent(inout) :: tol
      !! Absolute tolerance for the Newton solver
            integer, optional, intent(in) :: tol_mode
      
      ! Misc
            integer :: info, tol_mode_
            type(newton_dp_opts) :: opts
            character(len=3) :: file_prefix
      
		      tol_mode_ = optval(tol_mode, 1)

      ! Set up logging
            call logger%log_message('Starting newton iteration.', module=this_module)
      
      ! Define options for the Newton solver
            opts = newton_dp_opts(maxiter=40, ifbisect=.true.)
      
      ! Call to LightKrylov.
            if (tol_mode_ == 1) then
               call newton(sys, bf, gmres_rdp, info, atol=tol, options=opts, scheduler=nek_constant_tol)
            else
		         call newton(sys, bf, gmres_rdp, info, atol=tol, options=opts, scheduler=nek_dynamic_tol)
		      end if
      
      ! Outpost initial condition.
            file_prefix = 'nwt'
            call outpost_dnek(bf, file_prefix)

		      call logger%log_message('Exiting newton iteration.', module=this_module)
      
            return
         end subroutine newton_forced_periodic_orbit_torus

         subroutine mflow_newton_periodic_orbit_torus(sys, bf, mflow_target, tol, tol_mf, tol_mode, maxiter_newton)
            type(nek_system_torus), intent(inout) :: sys
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
            integer :: tol_mode_, maxiter_newton_
            integer :: nmf, nf, inwt, i, j
            real(dp) :: Wo, df0
            real(dp), allocatable :: dpds(:), phase(:), dpds_tmp(:)
            real(dp), allocatable :: mflow_old(:), mflow_new(:)
            real(dp), allocatable :: dmf(:), mf_err(:), deltaf(:), fpert(:)
            real(dp), allocatable :: jac(:,:)
				real(dp) :: dt_minmax(2)
            real(dp) :: tol_mf_inexact
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
				write(fmt,'("(A,",I0,"(1X,F16.10),A,*(F16.10,1X))")') nmf
            nf  = pipe%get_nf()      ! bumber of (complex!) forcing components
            if (nmf > nfft) then
               write(msg,'(A,I0,A,I0,A)') 'nmf= ', nmf, ' > nfft= ', nfft,'. Increase nfft in neklab_helix.'
               call nek_stop_error(msg, module=this_module, procedure='mflow_newton')
            else if (2*nmf - 1 /= nf) then
               write(msg,'(A,I0,A,I0,A)') 'nmf= ', nmf, ' and nf= ', nf, ' incompatible.'
               call nek_stop_error(msg, module=this_module, procedure='mflow_newton')
            end if
            allocate(dpds(nf), dpds_tmp(nf))
            allocate(dmf(nmf), mf_err(nmf), deltaf(nmf), fpert(nmf))
            allocate(jac(nmf,nmf))
            df0 = min(1.0e-6,100*tol) ! amplitude of forcing perturbation for finite difference approximation of gradient
            fpert(1)  = df0
            fpert(2:) = 20*df0
            ! get reference forcing and phase
            call pipe%get_dpds(dpds, phase)
      !
      ! Baseline Newton iteration on the initial guess of the forcing parameters & baseflow
      !
            call newton_forced_periodic_orbit_torus(sys, bf, tol, tol_mode)
      ! extract reference values for mass flow rate, initial mass flow error and flow solution
            call pipe%get_mflow_fft(mflow_old, if_amplitude=.true.)
            call nek2vec(ref, vx, vy, vz, pr, t)
            call nek_log_message('Reference solution set.', module=this_module, procedure='mflow_newton')
            mf_err = mflow_old(:nmf) - mflow_target
		! determine an approximation of the integration error for the mass flux
				call pipe%get_dt_minmax(dt_minmax)
				tol_mf_inexact = (sum(dt_minmax)*0.5)**2/100.0
				write(msg,'(A,1X,E16.8)') 'approximate mass flow computation error: ', tol_mf_inexact
            call nek_log_message(msg, module=this_module, procedure='mflow_newton')

      ! stamp logs
            call pipe%parameter_summary()
            write(msg,'(A,*(1X,F16.10))') padl('|df|:',   pad), fpert
            call nek_log_message(msg, module=this_module, procedure='mflow_newton')
            write(msg,'(A,*(1X,F16.10))') padl('mf_old:', pad), mflow_old(:nmf)
            call nek_log_message(msg, module=this_module, procedure='mflow_newton')
            write(msg,'(A,*(1X,F16.10))') padl('mf_err:', pad), mf_err
            call nek_log_message(msg, module=this_module, procedure='mflow_newton')
            !
            ! Main Newton iteration to converge the mass flow rate for each Fourier component
            !
            call nek_log_message('Begin mass flow Newton iteration', module=this_module, procedure='mflow_newton')
            write(msg,'(A,1X,E16.8)')    padl('mf_tol:', pad), tol_mf
            call nek_log_message(msg, module=this_module, procedure='mflow_newton')
            df_loop: do inwt = 1, maxiter_newton_
               write(step_id,'("Step ",I3,": ")') inwt
               write(msg,'(A,I0,A)') 'Begin mflow Newton step ', inwt, ' ...'
               call nek_log_information(msg, module=this_module, procedure='mflow_newton')
               do i = 1, nmf ! we need to compute dmf/df for all components
                  write(coef_id,'("Fourier coef. ",I2,": ")') i
                  write(msg,'(A,A,I0,A)') step_id, 'compute mflow gradient for Fourier coefficient ', i, ' ...'
                  call nek_log_information(msg, module=this_module, procedure='mflow_newton')
                  ! Set flow parameters
                  dpds_tmp = 0.0_dp
                  if (i == 1) then
                     dpds_tmp(1) = fpert(1)
                  else
                     j = 2*(i-1)
                     dpds_tmp(j  ) = cos(phase(i))*fpert(i)
                     dpds_tmp(j+1) = sin(phase(i))*fpert(i)
                  end if
						write(msg,'(A,A,A,*(1X,F16.10))') step_id, coef_id, padl('old frc:',pad), dpds
						call nek_log_message(msg, module=this_module, procedure='mflow_newton')
						write(msg,'(A,A,A,*(1X,F16.10))') step_id, coef_id, padl('prt frc:',pad), dpds_tmp
						call nek_log_message(msg, module=this_module, procedure='mflow_newton')
						dpds_tmp = dpds_tmp + dpds
						write(msg,'(A,A,A,*(1X,F16.10))') step_id, coef_id, padl('new frc:',pad), dpds_tmp
                  call nek_log_message(msg, module=this_module, procedure='mflow_newton')
                  call pipe%init_flow(dpds_tmp, Wo)
                  
                  ! reset baseflow
                  call bf%zero(); call bf%add(ref)
                  ! Run Newton-Krylov solver to find baseflow of perturbed system
                  call newton_forced_periodic_orbit_torus(sys, bf, tol, tol_mode)
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
               call pipe%init_flow(dpds, Wo)
               !
               ! Take Newton step for the forcing and compute periodic orbit
               !
               call newton_forced_periodic_orbit_torus(sys, bf, tol, tol_mode)
               ! extract reference values for mass flow rate, initial mass flow error, flow solution and forcing amplitudes and phase
               call pipe%get_mflow_fft(mflow_old, if_amplitude=.true.)
               call nek2vec(ref, vx, vy, vz, pr, t)
               call pipe%get_dpds(dpds, phase)
               mf_err = mflow_old(:nmf) - mflow_target
               call outpost_dnek(bf, 'nwf')
               call nek_log_information(step_id//'Final state:', module=this_module, procedure='mflow_newton')
               write(msg,fmt) step_id//'mf_state = ', mflow_old(:nmf) , ' | ext= ', mflow_old(:nmf+1:5) 
               call nek_log_information(msg, module=this_module, procedure='mflow_newton')
               write(msg,fmt) step_id//'mf_target= ', mflow_target, ' | tol= ', tol_mf
               call nek_log_information(msg, module=this_module, procedure='mflow_newton')
               write(msg,fmt) step_id//'mf_error = ', mf_err, ' | sum= ', sum(abs(mf_err))
               call nek_log_information(msg, module=this_module, procedure='mflow_newton')
               if (sum(abs(mf_err)) < tol_mf) then
						write(msg,'(A,I0,A)') 'Newton iteration converged after ', inwt, ' iterations.'
               	call nek_log_message(msg, module=this_module, procedure='mflow_newton')
                  exit df_loop ! converged
					end if
            end do df_loop
            call nek_log_message('Exiting mass flow Newton iteration.', module=this_module, procedure='mflow_newton')
            if (sum(abs(mf_err)) > tol_mf) then
               write(msg,'(A,I0,A)') 'Mass flux not converged after ', maxiter_newton_, 'steps.'
               call nek_stop_error(msg, module=this_module, procedure='mflow_newton')
            end if
         end subroutine mflow_newton_periodic_orbit_torus
      
         subroutine otd_analysis(OTD, opts_)
            type(nek_otd), intent(inout) :: OTD
            type(otd_opts), optional, intent(in) :: opts_
            type(otd_opts) :: opts
      ! internal
            real(dp), dimension(:), allocatable :: sigma
            real(dp), dimension(:, :), allocatable :: Lr, Phi, svec, G
            complex(dp), dimension(:), allocatable :: lambda
            complex(dp), dimension(:, :), allocatable :: eigvec
            type(nek_dvector), allocatable :: Lu(:)
      ! Misc
            integer :: i, j, r, log_level
            character(len=3) :: file_prefix
            character(len=128) :: msg
      
            if (present(opts_)) then
               opts = opts_
            else
               opts = otd_opts()
            end if
      
      ! Set up logging
            call logger_setup(nio=0, log_level=information_level, log_stdout=.false., log_timestamp=.true.)

            call logger%configuration(level=log_level)
      
      ! initialize OTD structure
            call OTD%init(opts)
      
      ! Allocate memory
            r = OTD%r
            allocate (sigma(r), svec(r, r)); sigma = 0.0_dp; svec = 0.0_dp
            allocate (lambda(r), eigvec(r, r)); lambda = 0.0_dp; eigvec = 0.0_dp
            allocate (Lr(r, r), Phi(r, r)); Lr = 0.0_dp; Phi = 0.0_dp
            allocate (Lu(r), source=OTD%baseflow); call zero_basis(Lu)
      
      ! Intgrate the nonlinear equations forward
            time = 0.0_dp
            do istep = 1, nsteps
               call nek_advance()
               if (istep >= opts%startstep) then
      ! load perturbations
                  do i = 1, r
                     call nek2vec(OTD%basis(i), vxp(:, i:i), vyp(:, i:i), vzp(:, i:i), prp(:, i:i), tp(:, :, i:i))
                  end do
      ! orthonormalize
                  if ((istep <= opts%startstep + 10) .or.
     $               mod(istep, opts%orthostep) == 0 .or.
     $               mod(istep, opts%printstep) == 0 .or.
     $               mod(istep, opts%iostep) == 0) then
      
                     if (log_level <= debug_level) then
                        allocate (G(r, r)); G = 0.0_dp
                        call innerprod(G, OTD%basis, OTD%basis)
                        write (msg, '(A,I5,A,*(1X,E10.3))') 'Step ', istep, ': norm.  err pre: ',  (G(i,i) - 1.0_dp, i=1, r)
                        call logger%log_information(msg, module=this_module, procedure='OTD main')
                        write (msg, '(A,I5,A,*(1X,E10.3))') 'Step ', istep, ': ortho. err pre: ', ((G(i,j), j=i+1, r), i=1, r)
                        call logger%log_information(msg, module=this_module, procedure='OTD main')
                     end if
         
                     call orthonormalize_basis(OTD%basis)

                     if (log_level <= debug_level) then
                        write (msg, '(A,I5,A,*(1X,E10.3))') 'Step ', istep, ': norm.  err post:',  (G(i,i) - 1.0_dp, i=1, r)
                        call logger%log_debug(msg, module=this_module, procedure='OTD main')
                        write (msg, '(A,I5,A,*(1X,E10.3))') 'Step ', istep, ': ortho. err post:', ((G(i,j), j=i+1, r), i=1, r)
                        call logger%log_debug(msg, module=this_module, procedure='OTD main')
                     end if
                  end if
      ! compute Lu
                  do i = 1, r
                     if (opts%trans) then
                        call OTD%apply_rmatvec(OTD%basis(i), Lu(i))
                     else
                        call OTD%apply_matvec(OTD%basis(i), Lu(i))
                     end if
                  end do
      ! compute reduced operator
                  call innerprod(Lr, OTD%basis, Lu)
      
                  Phi = 0.0_dp
                  do i = 1, r
                     do j = i + 1, r
                        Phi(i, j) = Lr(i, j)
                        Phi(j, 1) = -Lr(i, j)
                     end do
                  end do
      ! output projected modes
                  if (mod(istep, opts%printstep) == 0) then
                     call OTD%spectral_analysis(Lr, sigma, svec, lambda, eigvec, ifprint=.true.)
                  end if
      ! at the end of the step we copy data back to nek2vec
                  do i = 1, r
                     call vec2nek(vxp(:, i:i), vyp(:, i:i), vzp(:, i:i), prp(:, i:i), tp(:, :, i:i), OTD%basis(i))
                  end do
      ! project basis vectors and output modes
                  if (mod(istep, opts%iostep) == 0) then
                     if (mod(istep, opts%printstep) /= 0) then
                        call OTD%spectral_analysis(Lr, sigma, svec, lambda, eigvec, ifprint=.false.)
                     end if
                     call OTD%outpost_OTDmodes(eigvec)
                  end if
      ! output basis vectors
                  if (mod(istep, opts%iorststep) == 0) then
                     write (file_prefix, '(A)') 'rst'
                     call outpost_dnek(OTD%basis, file_prefix)
                  end if
      ! set the forcing
                  call OTD%generate_forcing(Lr, Phi)
               end if ! istep >= otd_startstep
            end do ! istep ... nsteps
            return
         end subroutine otd_analysis
      
         end module neklab_analysis
