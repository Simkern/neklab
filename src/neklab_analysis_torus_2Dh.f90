      module neklab_analysis_torus_2Dh
         use stdlib_stats_distribution_normal, only: normal => rvs_normal
         use stdlib_strings, only: padl, padr
         use stdlib_optval, only: optval
         use stdlib_linalg, only: diag, eye, det, inv
         use stdlib_logger, only: information_level, warning_level, debug_level, error_level, all_level, success
         use LightKrylov, only: atol_dp, dp, eigs, svds, save_eigenspectrum
         use LightKrylov, only: initialize_krylov_subspace, orthonormalize_basis, zero_basis, rand_basis
         use LightKrylov, only: linear_combination, innerprod
         use LightKrylov, only: newton, newton_dp_opts, gmres_rdp, gmres_dp_opts
         use LightKrylov_Logger
         use LightKrylov_Timing, only: timer => global_lightkrylov_timer
         use LightKrylov_AbstractSystems, only: abstract_system_rdp
         use LightKrylov_AbstractVectors, only: abstract_vector_rdp
         use neklab_2Dh_axisym, only: nek_advance_2Dh_axisym
         use neklab_vectors
         use neklab_linops
         use neklab_utils
         use neklab_nek_setup
         use neklab_otd
         use neklab_systems
         use neklab_analysis
         use neklab_newton_control
         use neklab_bf_buffer, only: bf_init, bf_finalize_module,
     &                               bf_set_prefix, bf_get_dt_minmax,
     &                               bf_get_nsteps, bf_get_time, bf_summary
      
         implicit none
         include "SIZE"
         include "TOTAL"
         include "ADJOINT"

         private
         character(len=*), parameter, private :: this_module = 'neklab_analysis_torus_2Dh'
         
         integer, parameter, private :: lv = lx1*ly1*lz1*lelv
         real(dp), private :: dpds_ = 0.0_dp       ! streamwise forcing
      
         public :: solve_fixed_point
         public :: steady_flowrate_newton
         public :: unsteady_flowrate_newton
         public :: shift_mflow_phase_upo, upo_taylor_test
      
      contains
   
         real(dp) function get_dpds() result(f)
      !! Current value of the (constant) streamwise forcing.
            f = dpds_
         end function get_dpds

         subroutine set_dpds(f)
      !! Sets the (constant) streamwise forcing and broadcasts it.
            real(dp), intent(in) :: f
            dpds_ = f
            call bcast(dpds_, wdsize)
         end subroutine set_dpds

      !==========================================================================
      !     INNER SOLVER WRAPPER
      !==========================================================================

         subroutine solve_fixed_point(sys, X, tol, tol_mode, info, maxiter)
      !! Thin wrapper around the (working) baseline Newton-Krylov fixed-point
      !! solver. Replace the body with the exact call used by your baseline
      !! driver if the signature differs -- nothing else in this file changes.
            class(abstract_system_rdp), intent(inout) :: sys
            class(abstract_vector_rdp), intent(inout) :: X
            real(dp), intent(in) :: tol
            integer, intent(in) :: tol_mode
            integer, intent(out) :: info
            integer, optional, intent(in) :: maxiter
      ! internal
            character(len=*), parameter :: this_procedure = 'solve_fixed_point'
            type(newton_dp_opts) :: opts
            opts = newton_dp_opts(maxiter=optval(maxiter, 40), ifbisect=.false.)
            if (tol_mode == 1) then
               call newton(sys, X, gmres_rdp, info, atol=tol, options=opts, scheduler=nek_constant_tol)
            else
               call newton(sys, X, gmres_rdp, info, atol=tol, options=opts, scheduler=nek_dynamic_tol)
            end if
         end subroutine solve_fixed_point

      !==========================================================================
      !     MAIN DRIVER
      !==========================================================================

         subroutine steady_flowrate_newton(sys, bf, Q_target, tol, tol_Q, tol_mode, maxiter, if_inexact, dQdf_guess)
            class(abstract_system_rdp), intent(inout) :: sys
      !! Steady 2Dh torus system (nek_system_torus_2Dh) whose fixed point is sought
            type(nek_dvector), intent(inout) :: bf
      !! In: initial guess. Out: fixed point at the prescribed flow rate
            real(dp), intent(in) :: Q_target
      !! Target streamwise volume flux, Q = int_A u_phi dA
            real(dp), intent(in) :: tol
      !! Absolute tolerance of the inner Newton-Krylov solver
            real(dp), intent(in) :: tol_Q
      !! Absolute tolerance on the flow-rate residual |Q - Q_target|
            integer, optional, intent(in) :: tol_mode
      !! Constant (1, default) or dynamic (2) tolerances in the inner solver
            integer, optional, intent(in) :: maxiter
      !! Maximum number of outer (flow-rate) Newton steps. Default: 10
            logical, optional, intent(in) :: if_inexact
      !! Loosen the inner tolerance while far from the target. Default: .true.
            real(dp), optional, intent(in) :: dQdf_guess
      !! Optional initial slope dQ/d(dpds). Default: Stokes estimate Q/dpds 
            ! internal
            character(len=*), parameter :: this_procedure = 'flowrate_newton'
            type(nek_dvector) :: Xold, dX
            integer :: tol_mode_, maxiter_, inwt, info
            logical :: inexact_
            real(dp) :: f, df, df_old, dfmax, ratio
            real(dp) :: Q, Q_old, res, dQdf, dQdf_new
            real(dp) :: tol_inner, noise_floor
            character(len=256) :: msg
            character(len=10) :: step_id
            integer, parameter :: pad = 18
            real(dp), parameter :: step_limit = 0.5_dp    ! max |df|/|f| per step
            real(dp), parameter :: pred_limit = 2.0_dp    ! max state extrapolation
            real(dp), parameter :: maxtol = 1.0e-04_dp    ! ceiling for inexact tol

      ! ---- optional arguments
            tol_mode_ = optval(tol_mode, 1)
            maxiter_ = optval(maxiter, 10)
            inexact_ = optval(if_inexact, .true.)

      ! ---- geometry
            call build_area_weights()
            if (area_cs <= 0.0_dp) then
               call nek_stop_error('Cross-section area is not positive.', this_module, this_procedure)
            end if

      ! ---- initial forcing
            f = get_dpds()
            if (abs(f) <= atol_dp) then
      ! Q(0) = 0 for a steady flow, so no slope can be inferred: take a probe value
               f = 1.0e-03_dp
               call set_dpds(f)
               write (msg, '(A,E16.8)') 'Initial forcing is zero. Using probe value dpds= ', f
               call nek_log_warning(msg, this_module, this_procedure)
            end if

      ! ---- stamp logs
            call nek_log_message('Flow-rate Newton configuration:', this_module, this_procedure)
            write (msg, '(3X,A,1X,E16.8)') padr('inner tol:', pad), tol
            call nek_log_message(msg, this_module, this_procedure)
            write (msg, '(3X,A,A)') padr('tol. scheduling:', pad), padl(merge('constant', 'dynamic ', tol_mode_ == 1), 16) 
            call nek_log_message(msg, this_module, this_procedure)
            write (msg, '(3X,A,A)') padr('inexact outer:', pad), padl(merge('yes', 'no ', inexact_), 16)
            call nek_log_message(msg, this_module, this_procedure)
            write (msg, '(3X,A,1X,E16.8)') padr('target Q:', pad), Q_target
            call nek_log_message(msg, this_module, this_procedure)
            write (msg, '(3X,A,1X,E16.8)') padr('Q tol:', pad), tol_Q
            call nek_log_message(msg, this_module, this_procedure)
            write (msg, '(3X,A,1X,E16.8)') padr('initial dpds:', pad), f
            call nek_log_message(msg, this_module, this_procedure)

      !
      ! ---- Baseline solve: fixed point for the initial forcing
      !
            call nek_log_message('Baseline fixed-point solve ...', this_module, this_procedure)
            tol_inner = tol
            call solve_fixed_point(sys, bf, tol_inner, tol_mode_, info)
            Q = get_flowrate_arr(bf%theta)
            res = Q - Q_target
            call log_state(0, f, Q, res)

      ! ---- initial slope.
      !      For creeping flow Q is exactly proportional to dpds, so Q/f is the
      !      exact derivative; with inertia the response is sub-linear, hence
      !      Q/f OVERestimates dQ/df and the first Newton step is a safe
      !      under-shoot. No extra nonlinear solve is needed to get started.
            dQdf = optval(dQdf_guess, Q/f)
            if (dQdf <= 0.0_dp) then
               dQdf = abs(Q/f)
               call nek_log_warning('Non-positive slope guess. Using |Q/f|.', this_module, this_procedure)
            end if
            write (msg, '(3X,A,1X,E16.8)') padr('initial dQ/df:', pad), dQdf
            call nek_log_message(msg, this_module, this_procedure)

      ! ---- reference state for the secant predictor
            call Xold%zero(); call Xold%add(bf)
            df_old = 0.0_dp
            df = 0.0_dp

      !
      ! ---- Outer (secant-Newton) iteration on the scalar forcing
      !
            call nek_log_message('Begin flow-rate Newton iteration', this_module, this_procedure)
            newton_loop: do inwt = 1, maxiter_ + 1
               if (abs(res) < tol_Q) then
                  if (tol_inner > tol) then
      ! the last state solve was inexact: confirm at the target tolerance
                     call nek_log_message('Polishing state at target tol ...', this_module, this_procedure)
                     tol_inner = tol
                     call solve_fixed_point(sys, bf, tol_inner, tol_mode_, info)
                     Q = get_flowrate_arr(bf%theta)
                     res = Q - Q_target
                  end if
                  if (abs(res) < tol_Q) then
                     write (msg, '(A,I0,A)') 'Flow-rate Newton converged after ', inwt - 1, ' step(s).' 
                     call nek_log_message(msg, this_module, this_procedure)
                     exit newton_loop
                  end if
               end if
               if (inwt > maxiter_) exit newton_loop
               write (step_id, '("Step ",I3,": ")') inwt

      ! ---- Newton step on the forcing, with a relative trust region
               df = -res/dQdf
               dfmax = step_limit*abs(f)
               if (abs(df) > dfmax) then
                  write (msg, '(A,A,E16.8,A,E16.8)') step_id, 'step limited: df= ', df, ' -> ', sign(dfmax, df)
                  call nek_log_message(msg, this_module, this_procedure)
                  df = sign(dfmax, df)
               end if

      ! ---- state predictor: extrapolate along the solution branch.
      !      dX is the state increment produced by the previous forcing step;
      !      scaling it by df/df_old gives a first-order guess of the new fixed
      !      point and typically saves one or two inner Newton iterations.
      !      On the first pass dX = 0, so this is a no-op.
               call dX%zero(); call dX%add(bf); call dX%sub(Xold)
               call Xold%zero(); call Xold%add(bf)
               if (abs(df_old) > atol_dp) then
                  ratio = max(-1.0_dp, min(pred_limit, df/df_old))
                  call dX%scal(ratio)
                  call bf%add(dX)
               end if

      ! ---- update forcing
               Q_old = Q
               f = f + df
               call set_dpds(f)
               write (msg, '(A,A,1X,E16.8,A,E16.8)') step_id, padr('dpds:', pad), f, '   df= ', df
               call nek_log_message(msg, this_module, this_procedure)

      ! ---- inner tolerance: no point resolving the state far below the
      !      accuracy needed to see the current flow-rate error.
               tol_inner = tol
               if (inexact_ .and. abs(res) > 10.0_dp*tol_Q) then
                  tol_inner = max(tol, min(maxtol, 0.05_dp*abs(res)/area_cs))
               end if
               write (msg, '(A,A,1X,E16.8)') step_id, padr('inner tol:', pad), tol_inner
               call nek_log_information(msg, this_module, this_procedure)

      ! ---- converge the fixed point at the new forcing
               call solve_fixed_point(sys, bf, tol_inner, tol_mode_, info)
               Q = get_flowrate_arr(bf%theta)
               res = Q - Q_target

      ! ---- secant update of dQ/df.
      !      Only accept it if the measured change in Q is well above the noise
      !      floor set by the inner tolerance, and if it keeps the physically
      !      required sign (Q is monotone increasing in dpds).
               noise_floor = max(10.0_dp*tol_inner*area_cs, atol_dp)
               if (abs(Q - Q_old) > noise_floor) then
                  dQdf_new = (Q - Q_old)/df
                  if (dQdf_new > 0.0_dp) then
                     dQdf = dQdf_new
                  else
                     call nek_log_warning(step_id//'Secant slope has the wrong sign. Slope kept.', this_module, this_procedure)
                  end if
               else
                  call nek_log_warning(step_id//'dQ below noise floor. Slope kept.',this_module, this_procedure)
               end if
               write (msg, '(A,A,1X,E16.8)') step_id, padr('dQ/df:', pad), dQdf
               call nek_log_information(msg, this_module, this_procedure)

               df_old = df
               call log_state(inwt, f, Q, res)

      ! ---- save intermediate solution
               call outpost_dnek(bf, 'nwq')
            end do newton_loop

      !
      ! ---- Output
      !
            call nek_log_message('Exiting flow-rate Newton iteration.', this_module, this_procedure)
            if (abs(res) > tol_Q) then
               write (msg, '(A,I0,A)') 'Flow rate not converged after ', maxiter_, ' steps.'
               call nek_stop_error(msg, this_module, this_procedure)
            end if
            call set_fldindex('BFQ', 1)
            call outpost_dnek(bf, 'BFQ')

            call nek_log_message('OUTPUT:', this_module, this_procedure)
            write (msg, '(3X,A,1X,E16.8)') padr('dpds:', pad), f
            call nek_log_message(msg, this_module, this_procedure)
            write (msg, '(3X,A,1X,E16.8)') padr('dQ/dpds:', pad), dQdf
            call nek_log_message(msg, this_module, this_procedure)
            write (msg, '(3X,A,1X,E16.8)') padr('Q:', pad), Q
            call nek_log_message(msg, this_module, this_procedure)
            write (msg, '(3X,A,1X,E16.8)') padr('Q target:', pad), Q_target
            call nek_log_message(msg, this_module, this_procedure)
            write (msg, '(3X,A,1X,E16.8)') padr('Q error:', pad), res
            call nek_log_message(msg, this_module, this_procedure)

         contains

            subroutine log_state(istp, f_, Q_, res_)
               integer, intent(in) :: istp
               real(dp), intent(in) :: f_, Q_, res_
               character(len=256) :: lmsg
               write (lmsg, '(A,I3,A,3(1X,E16.8))') 'FLOWRATE-NEWTON  it= ', istp, ' | dpds, Q, q_err = ', f_, Q_, res_
               call nek_log_message(lmsg, this_module, this_procedure)
            end subroutine log_state

         end subroutine flowrate_newton

         subroutine unsteady_flowrate_newton(sys, bf, Wo, kharm, dpds, mflow_target,
     &                                  tol, tol_mf, tol_mode, maxiter, maxiter_inner,
     &                                  buffer_base, df0, if_save_orbit)
            class(abstract_system_rdp), intent(inout) :: sys
      !! Pulsatile 2Dh system (nek_system_torus_upo_2Dh). Its jacobian is set here.
            type(nek_dvector), intent(inout) :: bf
      !! In: initial guess for the orbit. Out: converged orbit at the target
      !! flow rate. Runs are always restarts, so this is expected to arrive
      !! from a field file rather than from rest.
            real(dp), intent(in) :: Wo
      !! Womersley number of the fundamental. Fixed for the whole run:
      !! omega = Wo**2 * nu and T = 2*pi/omega are set once.
            integer, intent(in) :: kharm
      !! Number of active harmonics K. nf = 2K+1 forcing components,
      !! nmf = K+1 amplitude unknowns.
            real(dp), dimension(:), intent(inout) :: dpds
      !! In: initial forcing, nf values, HELIX convention (see the conversion
      !! note in neklab_newton_control). Out: converged forcing, same units.
            real(dp), dimension(:), intent(in) :: mflow_target
      !! Target flow-rate amplitudes, nmf values, helix normalisation: index 1
      !! is the mean, index k+1 is |mflow_k|, i.e. HALF the peak excursion of
      !! harmonic k. These are the same numbers a helix deck uses.
            real(dp), intent(in) :: tol
      !! Absolute tolerance of the inner Newton-Krylov solver.
            real(dp), intent(in) :: tol_mf
      !! Absolute tolerance on sum|mflow - target|.
            integer, optional, intent(in) :: tol_mode
      !! Constant (1, default) or dynamic (2) inner tolerance scheduling.
            integer, optional, intent(in) :: maxiter
      !! Maximum number of OUTER (flow-rate) Newton steps. Default 10.
            integer, optional, intent(in) :: maxiter_inner
      !! Maximum number of inner Newton iterations. Default 40.
            character(len=*), optional, intent(in) :: buffer_base
      !! Filename stem for the baseflow buffer. Default '2dtorus'.
            real(dp), optional, intent(in) :: df0
      !! Initial finite-difference step on the mean forcing. Default
      !! min(1e-5, 100*tol). Harmonics start ten times larger.
            logical, optional, intent(in) :: if_save_orbit
      !! Re-record the converged orbit under prefix 'b', so the Floquet run has
      !! a stable input that the next Newton solve will not overwrite.
      !! Default .true.
      ! internal
            character(len=*), parameter :: this_procedure = 'upo_flowrate_newton'
            type(nek_dvector) :: ref, res
            integer :: tol_mode_, maxiter_, maxiter_inner_
            integer :: nmf, nf, inwt, i, j, k, icnt, info
            logical :: save_orbit_, too_small
            real(dp) :: df0_, tol_df, tol_mf_inexact, mf_noise, detj, jscale
            real(dp) :: dt_minmax(2), period
            real(dp), dimension(lg) :: dpds_w, dpds_try, mflow_v
            real(dp), dimension(:), allocatable :: phase, mflow_old, mflow_new
            real(dp), dimension(:), allocatable :: dmf, mf_err, deltaf, fpert
            real(dp), dimension(:, :), allocatable :: jac
            character(len=256) :: msg, fmt
            character(len=10) :: step_id
            character(len=18) :: coef_id
            integer, parameter :: pad = 18
            integer, parameter :: max_fd_retry = 6

      ! ---- optional arguments
            tol_mode_ = optval(tol_mode, 1)
            maxiter_ = optval(maxiter, 10)
            maxiter_inner_ = optval(maxiter_inner, 40)
            save_orbit_ = optval(if_save_orbit, .true.)

      ! ---- sizes and checks
            nmf = kharm + 1
            nf = 2*kharm + 1
            if (nf > lg) then
               write (msg, '(A,I0,A,I0,A)') 'nf= ', nf, ' > lg= ', lg, '. Increase kmax_ctrl.'
               call nek_stop_error(msg, this_module, this_procedure)
            end if
            if (size(dpds) < nf) then
               write (msg, '(A,I0,A,I0)') 'dpds has ', size(dpds), ' entries but nf= ', nf
               call nek_stop_error(msg, this_module, this_procedure)
            end if
            if (size(mflow_target) < nmf) then
               write (msg, '(A,I0,A,I0)') 'mflow_target has ', size(mflow_target), ' entries but nmf= ', nmf
               call nek_stop_error(msg, this_module, this_procedure)
            end if
            write (fmt, '("(A,",I0,"(1X,F16.10),A,E16.8)")') nmf

            allocate (dmf(nmf), mf_err(nmf), deltaf(nmf), fpert(nmf))
            allocate (jac(nmf, nmf))

      ! ---- configure the pulsatile control and the baseflow buffer.
      !      The DRIVER owns the buffer lifetime: the system's nonlinear map
      !      only consumes it, and fails through check_init if this is skipped.
            call init_pulsatile(Wo, kharm, dpds)
            period = get_period()
            call bf_init(base=optval(buffer_base, '2dtorus'), write_chunks=.true., min_steps=50)
            call bf_set_prefix('n')
            sys%jacobian = nek_jacobian_torus_upo_2Dh()

      ! ---- finite-difference step sizes. The mean responds most strongly, so
      !      it gets the smallest probe; harmonics start an order up.
            df0_ = optval(df0, min(1.0e-05_dp, 100.0_dp*tol))
            tol_df = tol
            fpert(1) = df0_
            if (nmf > 1) fpert(2:) = 10.0_dp*df0_

      ! ---- stamp logs
            call nek_log_message('Pulsatile flow-rate Newton configuration:', this_module, this_procedure)
            write (msg, '(3X,A,1X,E16.8)') padr('Womersley:', pad), Wo
            call nek_log_message(msg, this_module, this_procedure)
            write (msg, '(3X,A,1X,E16.8)') padr('period T:', pad), period
            call nek_log_message(msg, this_module, this_procedure)
            write (msg, '(3X,A,1X,I0,A,I0,A)') padr('harmonics K:', pad), kharm,
     &         '  (nmf= ', nmf, ' amplitude unknowns)'
            call nek_log_message(msg, this_module, this_procedure)
            write (msg, '(3X,A,1X,E16.8)') padr('inner tol:', pad), tol
            call nek_log_message(msg, this_module, this_procedure)
            write (msg, '(3X,A,A)') padr('tol. scheduling:', pad), padl(merge('constant', 'dynamic ', tol_mode_ == 1), 16)
            call nek_log_message(msg, this_module, this_procedure)
            write (msg, '(3X,A,*(1X,F16.12))') padr('forcing |df|:', pad), (fpert(i), i=1, nmf)
            call nek_log_message(msg, this_module, this_procedure)
            write (msg, '(3X,A,*(1X,F16.12))') padr('target mflow:', pad), (mflow_target(i), i=1, nmf)
            call nek_log_message(msg, this_module, this_procedure)
            write (msg, '(3X,A,1X,E16.8)') padr('mflow tol:', pad), tol_mf
            call nek_log_message(msg, this_module, this_procedure)
            if (kharm >= 2) then
               call nek_log_warning('K >= 2: the amplitude spectrum is controlled but the '//
     &            'relative phases between harmonics are not, so the waveform shape depends '//
     &            'on the incoming forcing phases.', this_module, this_procedure)
            end if

      ! ---- reference forcing and phases
            call get_dpds_fourier(dpds_w, phase)

      !
      ! ---- Baseline: converge the orbit at the initial forcing
      !
            call nek_log_message('Baseline periodic-orbit solve ...', this_module, this_procedure)
            call solve_fixed_point(sys, bf, tol, tol_mode_, info, maxiter=maxiter_inner_)
            call extract_mflow(mflow_v, mflow_old)
            call ref%zero(); call ref%add(bf)
            mf_err = mflow_old(:nmf) - mflow_target(:nmf)

      ! ---- quadrature error floor. The Fourier coefficients come from a
      !      trapezoidal accumulation along the trajectory, so they cannot be
      !      more accurate than O(dt**2). Chasing a flow-rate tolerance below
      !      that is chasing the quadrature, not the physics.
            call bf_get_dt_minmax(dt_minmax)
            tol_mf_inexact = (0.5_dp*sum(dt_minmax))**2/100.0_dp
            mf_noise = max(tol_mf_inexact, 10.0_dp*tol_df)
            write (msg, '(A,1X,E16.8,A,I0,A)') 'Approximate mflow quadrature error: ', tol_mf_inexact,
     &         '  (', bf_get_nsteps(), ' steps/period)'
            call nek_log_message(msg, this_module, this_procedure)
            if (tol_mf_inexact > tol_mf) then
               call nek_log_warning('Requested mflow tol is below the estimated quadrature error. '//
     &            'Lower the CFL target to tighten it.', this_module, this_procedure)
            end if

            call control_summary()
            call forcing_summary()
            call nek_log_information('Initial state:', this_module, this_procedure)
            write (msg, '(A,*(1X,F16.10))') 'mf_state  = ', (mflow_old(i), i=1, nmf)
            call nek_log_information(msg, this_module, this_procedure)
            write (msg, fmt) 'mf_target = ', (mflow_target(i), i=1, nmf), ' | tol= ', tol_mf
            call nek_log_information(msg, this_module, this_procedure)
            write (msg, fmt) 'mf_error  = ', (mf_err(i), i=1, nmf), ' | sum= ', sum(abs(mf_err))
            call nek_log_information(msg, this_module, this_procedure)

      !
      ! ---- Outer Newton on the forcing amplitudes
      !
            call nek_log_message('Begin mass-flow Newton iteration', this_module, this_procedure)
            df_loop: do inwt = 1, maxiter_
               if (sum(abs(mf_err)) < tol_mf) then
                  call nek_log_message('Initial forcing already meets the flow-rate target.',
     &               this_module, this_procedure)
                  exit df_loop
               end if
               write (step_id, '("Step ",I3,": ")') inwt
               write (msg, '(A,I0,A)') 'Begin mflow Newton step ', inwt, ' ...'
               call nek_log_information(msg, this_module, this_procedure)

      ! ---- finite-difference Jacobian, one inner solve per amplitude
               do i = 1, nmf
                  write (coef_id, '("Fourier coef. ",I2,": ")') i
                  write (msg, '(A,A,I0,A)') step_id, 'compute mflow gradient for component ', i, ' ...'
                  call nek_log_information(msg, this_module, this_procedure)
                  icnt = 0
                  fd_retry: do
      ! Perturb amplitude i along its FROZEN phase direction, in the direction
      ! of the root so that the secant is taken on the side we are heading for.
                     dpds_try = 0.0_dp
                     fpert(i) = -sign(fpert(i), mf_err(i))
                     if (i == 1) then
                        dpds_try(1) = fpert(1)
                     else
                        j = 2*(i - 1)
                        dpds_try(j) = cos(phase(i))*fpert(i)
                        dpds_try(j + 1) = sin(phase(i))*fpert(i)
                     end if
                     write (msg, '(A,A,A,*(1X,F16.10))') step_id, coef_id, padl('prt frc:', 10),
     &                  (dpds_try(k), k=1, nf)
                     call nek_log_message(msg, this_module, this_procedure)
                     dpds_try(:nf) = dpds_try(:nf) + dpds_w(:nf)
                     call set_dpds_fourier(dpds_try(:nf))

      ! Restart the inner solve from the reference orbit every time, so the
      ! measured difference is a property of the forcing and not of where the
      ! previous solve happened to stop.
                     call bf%zero(); call bf%add(ref)
                     call solve_fixed_point(sys, bf, tol_df, tol_mode_, info, maxiter=maxiter_inner_)
                     call extract_mflow(mflow_v, mflow_new)
                     dmf = mflow_new(:nmf) - mflow_old(:nmf)
                     icnt = icnt + 1

      ! A response buried in the quadrature/solver noise carries no gradient
      ! information. Detect that from the MEASUREMENT rather than from whether
      ! the inner Newton happened to iterate: it is the quantity that actually
      ! matters and it does not depend on the inner solver's exit convention.
                     too_small = maxval(abs(dmf)) < mf_noise
                     if (.not. too_small) exit fd_retry
                     if (icnt >= max_fd_retry) then
                        write (msg, '(A,I0,A,I0,A)') 'Component ', i, ': no measurable response after ',
     &                     icnt, ' attempts. The forcing may be decoupled from this harmonic.'
                        call nek_stop_error(msg, this_module, this_procedure)
                     end if
                     fpert(i) = 10.0_dp*fpert(i)
                     write (msg, '(A,I0,A,E16.8)') 'Response below noise floor for component ', i,
     &                  ': reset |df| = ', abs(fpert(i))
                     call nek_log_message(msg, this_module, this_procedure)
                  end do fd_retry

                  write (msg, '(A,A,A,*(1X,F16.10))') step_id, coef_id, padl('dmflow:', 10), (dmf(k), k=1, nmf)
                  call nek_log_message(msg, this_module, this_procedure)

      ! jac(j,i) = d(mflow_j)/d(f_i). Note the index order: the Newton step
      ! solves jac * deltaf = -mf_err, so the ROW index must be the constraint
      ! and the COLUMN index the unknown. (mflow_newton in
      ! neklab_analysis_torus fills jac(i,j) here and then applies inv(jac)
      ! directly, i.e. it uses the transpose. The response is close to diagonal
      ! -- linear Womersley theory makes it exactly diagonal -- so the error is
      ! small and the iteration still converges, just linearly rather than
      ! quadratically. Worth fixing there too.)
                  do j = 1, nmf
                     jac(j, i) = dmf(j)/fpert(i)
                  end do
               end do

               write (msg, '(A,A)') step_id, 'mflow jacobian  d(mflow_row)/d(f_col)'
               call nek_log_message(msg, this_module, this_procedure)
               do i = 1, nmf
                  write (msg, '(A,4X,*(1X,F16.10))') step_id, (jac(i, j), j=1, nmf)
                  call nek_log_message(msg, this_module, this_procedure)
               end do

      ! Scaled singularity test. An exact-zero determinant test is useless in
      ! floating point: compare against the product of the row magnitudes.
               jscale = 1.0_dp
               do i = 1, nmf
                  jscale = jscale*max(maxval(abs(jac(i, :))), atol_dp)
               end do
               detj = det(jac)
               if (abs(detj) <= 1.0e-12_dp*jscale) then
                  write (msg, '(A,E16.8,A,E16.8)') 'mflow Jacobian is numerically singular: det= ',
     &               detj, ', scale= ', jscale
                  call nek_stop_error(msg, this_module, this_procedure)
               end if

      ! ---- Newton step on the amplitudes
               deltaf = -matmul(inv(jac), mf_err)
               write (msg, '(A,A,*(1X,F16.10))') step_id, padl('|f_step|:', pad), (deltaf(i), i=1, nmf)
               call nek_log_message(msg, this_module, this_procedure)

               dpds_w(1) = dpds_w(1) + deltaf(1)
               do i = 2, nmf
                  j = 2*(i - 1)
                  dpds_w(j) = dpds_w(j) + cos(phase(i))*deltaf(i)
                  dpds_w(j + 1) = dpds_w(j + 1) + sin(phase(i))*deltaf(i)
               end do

               call nek_log_message('Forcing prior to Newton step:', this_module, this_procedure)
               call forcing_summary()
               call set_dpds_fourier(dpds_w(:nf))

      ! ---- converge the orbit at the updated forcing
               call bf%zero(); call bf%add(ref)
               call solve_fixed_point(sys, bf, tol, tol_mode_, info, maxiter=maxiter_inner_)
               call extract_mflow(mflow_v, mflow_old)
               call ref%zero(); call ref%add(bf)
               call get_dpds_fourier(dpds_w, phase)
               mf_err = mflow_old(:nmf) - mflow_target(:nmf)

               call nek_log_information(step_id//'Final state:', this_module, this_procedure)
               write (msg, '(A,*(1X,F16.10))') step_id//'mf_state = ', (mflow_old(i), i=1, nmf)
               call nek_log_information(msg, this_module, this_procedure)
               write (msg, fmt) step_id//'mf_target= ', (mflow_target(i), i=1, nmf), ' | tol= ', tol_mf
               call nek_log_information(msg, this_module, this_procedure)
               write (msg, fmt) step_id//'mf_error = ', (mf_err(i), i=1, nmf), ' | sum= ', sum(abs(mf_err))
               call nek_log_information(msg, this_module, this_procedure)

               if (sum(abs(mf_err)) < tol_mf) then
                  write (msg, '(A,I0,A)') 'Mass-flow Newton converged after ', inwt, ' step(s).'
                  call nek_log_message(msg, this_module, this_procedure)
                  exit df_loop
               else
                  call set_fldindex('nwf', 1)
                  call outpost_dnek(bf, 'nwf')
               end if
            end do df_loop

      !
      ! ---- Output
      !
            call nek_log_message('Exiting mass-flow Newton iteration.', this_module, this_procedure)
            if (sum(abs(mf_err)) > tol_mf) then
               write (msg, '(A,I0,A)') 'Mass flow not converged after ', maxiter_, ' steps.'
               call nek_stop_error(msg, this_module, this_procedure)
            end if

            dpds = 0.0_dp
            dpds(1:nf) = dpds_w(1:nf)

      ! ---- record the converged orbit under its own prefix, so the Floquet
      !      run reads a set of chunk files that the next Newton solve will not
      !      overwrite. This costs one extra nonlinear pass.
            if (save_orbit_) then
               call nek_log_message('Recording the converged orbit under prefix b ...',
     &            this_module, this_procedure)
               call bf_set_prefix('b')
               call sys%response(bf, res, tol)
               write (msg, '(3X,A,1X,E16.8)') padr('|F(X)| :', pad), res%norm()
               call nek_log_message(msg, this_module, this_procedure)
               call bf_summary()
               call bf_set_prefix('n')
            end if

            call set_fldindex('BFP', 1)
            call outpost_dnek(bf, 'BFP')

            call nek_log_message('OUTPUT:', this_module, this_procedure)
            call forcing_summary()
            write (msg, '(3X,A,*(1X,F16.12))') padr('mf_state:', pad), (mflow_old(i), i=1, nmf)
            call nek_log_message(msg, this_module, this_procedure)
            write (msg, '(3X,A,*(1X,F16.12))') padr('mf_target:', pad), (mflow_target(i), i=1, nmf)
            call nek_log_message(msg, this_module, this_procedure)
            write (msg, '(3X,A,1X,E16.8)') padr('mflow error:', pad), sum(abs(mf_err))
            call nek_log_message(msg, this_module, this_procedure)
            write (msg, '(3X,A,1X,I0)') padr('steps/period:', pad), bf_get_nsteps()
            call nek_log_message(msg, this_module, this_procedure)

         end subroutine unsteady_flowrate_newton

      !====================================================================
      !     PHASE GAUGE
      !====================================================================

         subroutine shift_mflow_phase_upo(bf, tol, cfl_limit)
      !! Uses the time-origin freedom to put the FUNDAMENTAL flow-rate harmonic
      !! on a pure cosine, phase_1 = 0, matching the helix convention.
      !!
      !! Advancing the state by s and re-referencing the forcing gives the same
      !! physical orbit seen from a shifted origin. Under t -> t + s:
      !!
      !!    forcing   phase_k -> phase_k + k*omega*s   (dpds rotates by +k*w*s)
      !!    flow rate phase_k -> phase_k - k*omega*s
      !!
      !! -- the opposite signs are a consequence of helix projecting the forcing
      !! onto exp(-i k w t) and the flow rate onto exp(+i k w t). Choosing
      !! s = phase_1/omega therefore zeroes the fundamental's flow-rate phase.
      !!
      !! ONE shift buys ONE phase. For K >= 2 the higher harmonics are rotated
      !! consistently but their phases are NOT zeroed, and that is the honest
      !! outcome: their relative phases are physical, not gauge.
      !!
      !! (neklab_analysis_torus's shift_mflow_phase_torus loops over harmonics
      !! applying phase(i)/omega each time, cumulatively, from phases sampled
      !! once before the loop. That is exact for nmf = 2 and does not do what
      !! its name says for nmf >= 3.)
            type(nek_dvector), intent(inout) :: bf
      !! Orbit to be shifted, in place.
            real(dp), intent(in) :: tol
            real(dp), optional, intent(in) :: cfl_limit
      ! internal
            character(len=*), parameter :: this_procedure = 'shift_mflow_phase_upo'
            character(len=256) :: msg
            real(dp), dimension(lg) :: dpds_w, dpds_new, mflow_v
            real(dp), dimension(:), allocatable :: phase
            real(dp) :: omega, period, shift, tend, ang, dc, ds
            integer :: k, i, kc, nf
            real(dp), parameter :: tol_shift = 1.0e-08_dp

            omega = get_omega()
            period = get_period()
            kc = get_kctrl()
            nf = get_nf()
            if (kc < 1) then
               call nek_log_warning('No harmonics: nothing to gauge.', this_module, this_procedure)
               return
            end if

            call extract_mflow(mflow_v, phase=phase)
            write (msg, '(A,*(1X,F16.10))') 'mflow phase (before) =', (phase(k), k=1, kc + 1)
            call nek_log_message(msg, this_module, this_procedure)

            shift = phase(2)/omega
            if (abs(shift) < tol_shift) then
               call nek_log_message('Fundamental is already in phase. Nothing to do.',
     &            this_module, this_procedure)
               return
            end if
      ! Time only runs forward, so a negative shift is taken modulo the period.
            tend = shift
            if (tend < 0.0_dp) tend = period + tend
            write (msg, '(A,E16.8,A,E16.8)') 'Shifting the orbit by s= ', shift, ', integrating for ', tend
            call nek_log_message(msg, this_module, this_procedure)

      ! ---- advance the state by tend at the CURRENT forcing. No recording:
      !      this is a re-gauge, not a residual evaluation.
            call vec2nek(vx, vy, vz, pr, t, bf)
            call setup_nonlinear_solver(variable_dt = .true.,
     &                                  endtime     = tend,
     &                                  solve_temperature = .true.,
     &                                  cfl_limit   = optval(cfl_limit, 0.4_dp),
     &                                  vtol        = tol*0.1,
     &                                  ptol        = tol*0.1)
            time = 0.0_dp
            istep = 0
            do while (lastep == 0)
               istep = istep + 1
               call nek_advance()
            end do
            call nek2vec(bf, vx, vy, vz, pr, t)

      ! ---- rotate the forcing by the same shift
            call get_dpds_fourier(dpds_w)
            dpds_new = dpds_w
            do k = 1, kc
               i = 2*k
               ang = k*omega*shift
               dc = dpds_w(i); ds = dpds_w(i + 1)
               dpds_new(i) = dc*cos(ang) - ds*sin(ang)
               dpds_new(i + 1) = dc*sin(ang) + ds*cos(ang)
            end do
            call set_dpds_fourier(dpds_new(:nf))

            call nek_log_message('Updated forcing:', this_module, this_procedure)
            call forcing_summary()
            if (kc >= 2) then
               call nek_log_warning('Only the fundamental has been gauged. Higher harmonics keep '//
     &            'their relative phases, which are physical.', this_module, this_procedure)
            end if
         end subroutine shift_mflow_phase_upo

      !====================================================================
      !     CONSISTENCY CHECK
      !====================================================================

         subroutine upo_taylor_test(sys, X0, tol, nlevels)
      !! Verifies that the replayed Jacobian really is the derivative of the
      !! recorded nonlinear map:
      !!
      !!    e(eps) = || F(X0 + eps*dX) - F(X0) - eps*J*dX ||  ~  O(eps**2)
      !!
      !! This is THE check to run before trusting anything else here. Clean
      !! second order means the buffer replay reproduces the nonlinear
      !! trajectory step for step. A stall at first order, or a plateau, means
      !! the two passes are running on different time grids -- which is exactly
      !! the failure mode a CFL-adaptive dt can produce if the step sequence
      !! drifts between the residual evaluation and the matvec.
      !!
      !! ORDERING MATTERS: J*dX must be evaluated immediately after F(X0),
      !! because the perturbed residual evaluations overwrite the buffer.
            class(abstract_system_rdp), intent(inout) :: sys
            type(nek_dvector), intent(in) :: X0
            real(dp), intent(in) :: tol
            integer, optional, intent(in) :: nlevels
      ! internal
            character(len=*), parameter :: this_procedure = 'upo_taylor_test'
            character(len=256) :: msg
            type(nek_dvector) :: F0, Fp, dX, JdX, Xp, err
            integer :: nlev, l
            real(dp) :: eps, e, e_prev, order, nrm

            nlev = optval(nlevels, 5)
            sys%jacobian = nek_jacobian_torus_upo_2Dh()

            call nek_log_message('Taylor test on the 2Dh periodic-orbit system:', this_module, this_procedure)

      ! ---- unit-norm random direction
            call dX%rand()
            nrm = dX%norm()
            if (nrm <= atol_dp) call nek_stop_error('Degenerate random direction.', this_module, this_procedure)
            call dX%scal(1.0_dp/nrm)

      ! ---- base residual, then the matvec while the buffer still holds this
      !      trajectory
            call sys%response(X0, F0, tol)
            sys%jacobian%X = X0
            call sys%jacobian%matvec(dX, JdX)
            write (msg, '(3X,A,1X,E16.8,A,I0)') '|F(X0)| = ', F0%norm(), ',  steps/period = ', bf_get_nsteps()
            call nek_log_message(msg, this_module, this_procedure)

            e_prev = 0.0_dp
            do l = 1, nlev
               eps = 10.0_dp**(-l)
               call Xp%zero(); call Xp%add(X0)
               call Xp%axpby(1.0_dp, dX, eps)
               call sys%response(Xp, Fp, tol)
      ! err = Fp - F0 - eps*J*dX
               call err%zero(); call err%add(Fp); call err%sub(F0)
               call err%axpby(1.0_dp, JdX, -eps)
               e = err%norm()
               if (l == 1) then
                  write (msg, '(3X,A,E12.4,A,E16.8)') 'eps= ', eps, '   err= ', e
               else
                  order = 0.0_dp
                  if (e > 0.0_dp .and. e_prev > 0.0_dp) order = log10(e_prev/e)
                  write (msg, '(3X,A,E12.4,A,E16.8,A,F8.3)') 'eps= ', eps, '   err= ', e,
     &               '   observed order= ', order
               end if
               call nek_log_message(msg, this_module, this_procedure)
               e_prev = e
            end do
            call nek_log_message('Expected observed order: 2.0 until the inner tolerance floors it.',
     &         this_module, this_procedure)
         end subroutine upo_taylor_test
      
         end module neklab_analysis_torus_2Dh
