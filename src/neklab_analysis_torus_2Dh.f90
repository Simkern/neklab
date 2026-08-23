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
      
         implicit none
         include "SIZE"
         include "TOTAL"
         include "ADJOINT"

         private
         character(len=*), parameter, private :: this_module = 'neklab_analysis_torus_2Dh'
         
         integer, parameter, private :: lv = lx1*ly1*lz1*lelv
         real(dp), dimension(lv), private :: bm_area      ! area mass matrix
         real(dp), private :: area_cs      = 0.0_dp       ! cross-section area
         logical,  private :: area_defined = .false.
         real(dp), private :: dpds_        = 0.0_dp       ! streamwise forcing
      
         !public :: build_area_weights
         !public :: get_area, get_flowrate, get_flowrate_nek
         !public :: get_dpds, set_dpds
         public :: solve_fixed_point
         public :: flowrate_newton
         public :: combined_flowrate_newton
      
      contains
      !--------------------------------------------------------------------------
      !-----                                                                -----
      !-----     FLOW-RATE CONTROLLED NEWTON SOLVER  --  STEADY, 2Dh        -----
      !-----                                                                -----
      !--------------------------------------------------------------------------
      !
      !  Drop-in replacement for mflow_newton, restricted (for now) to STEADY
      !  states of the 2Dh torus formulation:
      !
      !     * the state is a nek_dvector (vx,vy,pr,t) where t carries u_phi,
      !     * the control is the single scalar streamwise forcing dpds,
      !     * the constraint is  Q(dpds) = int_A t dA = Q_target.
      !
      !  Everything the old routine got from 'pipe' (surface integral, forcing
      !  get/set, area) is defined below, so this file has no dependency on
      !  neklab_helix.
      !
      !  ---------------------------------------------------------------------
      !  ADD TO THE MODULE HEADER (neklab_analysis_torus, or its 2Dh successor)
      !  ---------------------------------------------------------------------
      !
      !     use LightKrylov, only: newton, newton_dp_opts, gmres_rdp
      !     use LightKrylov, only: atol_dp, dp
      !     use LightKrylov_Logger, only: check_info
      !     use neklab_2Dh_axisym, only: nek_advance_2Dh_axisym
      !
      !     integer, parameter, private :: lv = lx1*ly1*lz1*lelv
      !
      !     real(dp), dimension(lv), private :: bm_area      ! area mass matrix
      !     real(dp), private :: area_cs      = 0.0_dp       ! cross-section area
      !     logical,  private :: area_defined = .false.
      !     real(dp), private :: dpds_        = 0.0_dp       ! streamwise forcing
      !
      !     public :: flowrate_newton
      !     public :: get_flowrate, get_flowrate_nek, get_area
      !     public :: get_dpds, set_dpds
      !
      !  userq (or userf, depending on how you inject the forcing) then simply
      !  does:   qvol = get_dpds()
      !
      !--------------------------------------------------------------------------

      !==========================================================================
      !     GEOMETRIC / DIAGNOSTIC HELPERS  (former pipe%... routines)
      !==========================================================================

!         subroutine build_area_weights(force)
!      !! Builds the mass matrix for surface integrals over the torus
!      !! cross-section (phi = const plane) and the cross-section area.
!      !!
!      !! The cross-section of the torus is a plane surface with area element
!      !! dA = dR dz, i.e. WITHOUT the radial (Jacobian) weight that Nek's bm1
!      !! carries when the mesh is run in axisymmetric mode. If your setup uses
!      !! ifaxis = .false. and handles the torus metric through the coefficient
!      !! arrays of neklab_2Dh_axisym (alphaR_coef, diag_shift, ...), bm1 is
!      !! already the plain area weight and the branch below is a no-op.
!            logical, optional, intent(in) :: force
!      ! internal
!            character(len=*), parameter :: this_procedure = 'build_area_weights'
!            real(dp), external :: glsum
!            integer :: n
!            character(len=128) :: msg
!            if (area_defined .and. .not. optval(force, .false.)) return
!            n = lx1*ly1*lz1*nelv
!            if (ifaxis) then
!               call invcol3(bm_area, bm1, ym1, n)   ! strip the radius weight
!            else
!               call copy(bm_area, bm1, n)
!            end if
!            area_cs = glsum(bm_area, n)
!            area_defined = .true.
!            write (msg, '(A,E16.8)') 'Cross-section area A = ', area_cs
!            call nek_log_information(msg, this_module, this_procedure)
!         end subroutine build_area_weights
!
!         real(dp) function get_area() result(area)
!      !! Cross-sectional area of the torus, A = int_A dA.
!            call build_area_weights()
!            area = area_cs
!         end function get_area
!
!         real(dp) function get_flowrate_nek() result(Q)
!      !! Streamwise volume flux from the current Nek fields, Q = int_A t dA / pi.
!            real(dp), external :: glsc2
!            integer :: n
!            call build_area_weights()
!            n = lx1*ly1*lz1*nelv
!            Q = glsc2(t(1, 1, 1, 1, 1), bm_area, n) / area_cs
!         end function get_flowrate_nek
!         
!         real(dp) function get_flowrate(vec) result(Q)
!      !! Streamwise volume flux carried by a nek_dvector, Q = int_A theta dA.
!            type(nek_dvector), intent(in) :: vec
!            ! internal
!            real(dp), external :: glsc2
!            integer :: n
!            call build_area_weights()
!            n = lx1*ly1*lz1*nelv
!            Q = glsc2(vec%theta(1, 1), bm_area, n) / area_cs
!         end function get_flowrate

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

         subroutine flowrate_newton(sys, bf, Q_target, tol, tol_Q, tol_mode, maxiter, if_inexact, dQdf_guess)
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

         subroutine combined_flowrate_newton(sys, bf, Q_target, tol, kharm, omega, tol_mode, maxiter, n_relax, g_scale, if_verify)
                   class(abstract_system_rdp), intent(inout) :: sys
             !! Bordered system (nek_system_bordered). Its jacobian is set here.
                   type(nek_bordered_dvector), intent(inout) :: bf
             !! In:  initial guess. bf%g(1:nctrl) MUST hold the initial forcing.
             !! Out: fixed point and forcing satisfying the flow-rate constraint.
                   real(dp), dimension(:), intent(in) :: Q_target
             !! Target flow-rate Fourier coefficients (nctrl = 2*kharm+1 values).
                   real(dp), intent(in) :: tol
             !! Target absolute tolerance for the Newton solver.
                   integer, optional, intent(in) :: kharm
             !! Number of harmonics K. Default 0 (steady).
                   real(dp), optional, intent(in) :: omega
             !! Fundamental angular frequency. Required for K > 0.
                   integer, optional, intent(in) :: tol_mode
             !! Constant (1, default) or dynamic (2) inner tolerance scheduling.
                   integer, optional, intent(in) :: maxiter
             !! Maximum number of Newton iterations. Default 40.
                   integer, optional, intent(in) :: n_relax
             !! Number of relaxation horizons run before the Newton solve, to remove
             !! transients and to calibrate the control weights. Default 1.
                   real(dp), dimension(:), optional, intent(in) :: g_scale
             !! Override for the control scaling. Default: |g_i|, falling back to
             !! |g_1|/(1+k) for harmonics that start at zero.
                   logical, optional, intent(in) :: if_verify
             !! Re-evaluate the residual after convergence and report it. Default .true.
             ! internal
                   character(len=*), parameter :: this_procedure = 'combined_flowrate_newton'
                   type(nek_bordered_dvector) :: res
                   type(newton_dp_opts) :: opts
                   integer :: kharm_, tol_mode_, maxiter_, n_relax_, info, i, k, n
                   logical :: verify_
                   real(dp) :: xnorm, area, dtn
                   real(dp), dimension(lg) :: g0, gsc, qbar, qerr
                   character(len=256) :: msg
                   integer, parameter :: pad = 18
                   real(dp), parameter :: g_probe = 1.0e-03_dp
       
             ! ---- optional arguments
                   kharm_ = optval(kharm, 0)
                   tol_mode_ = optval(tol_mode, 1)
                   maxiter_ = optval(maxiter, 40)
                   n_relax_ = optval(n_relax, 1)
                   verify_ = optval(if_verify, .true.)
                   n = 2*kharm_ + 1
                   if (size(Q_target) < n) then
                      write (msg, '(A,I0,A,I0)') 'Q_target has ', size(Q_target), ' entries but nctrl= ', n
                      call nek_stop_error(msg, this_module, this_procedure)
                   end if
       
             ! ---- initial forcing, taken from the vector itself
                   g0 = bf%g
                   if (abs(g0(1)) <= atol_dp) then
             ! Q(0) = 0, so a zero mean forcing gives no scale to work from.
                      g0(1) = g_probe
                      write (msg, '(A,E16.8)') 'Zero initial mean forcing. Using probe value g= ', g0(1) 
                      call nek_log_warning(msg, this_module, this_procedure)
                   end if
       
             ! ---- configure the control block
                   call init_control(kharm_, omega=optval(omega, 0.0_dp), target=Q_target(1:n), g0=g0(1:n))
                   bf%g = 0.0_dp
                   bf%g(1:n) = g0(1:n)
                   area = get_area()
       
             ! ---- stamp logs
                   call nek_log_message('Bordered flow-rate Newton:', this_module, this_procedure)
                   write (msg, '(3X,A,1X,I0)') padr('harmonics K:', pad), kharm_
                   call nek_log_message(msg, this_module, this_procedure)
                   write (msg, '(3X,A,1X,I0)') padr('control dofs:', pad), n
                   call nek_log_message(msg, this_module, this_procedure)
                   write (msg, '(3X,A,1X,E16.8)') padr('target tol:', pad), tol
                   call nek_log_message(msg, this_module, this_procedure)
                   write (msg, '(3X,A,A)') padr('tol. scheduling:', pad), padl(merge('constant', 'dynamic ', tol_mode_ == 1), 16)
                   call nek_log_message(msg, this_module, this_procedure)
                   write (msg, '(3X,A,1X,E16.8)') padr('area:', pad), area
                   call nek_log_message(msg, this_module, this_procedure)
       
             !
             ! ---- Relaxation: a few horizons at fixed forcing.
             !      This is not a Newton step. It kills the worst transients (which
             !      would otherwise show up as a huge first residual) and gives the
             !      state norm used to calibrate the control weights.
             !
                   if (n_relax_ > 0) then
                      write (msg, '(A,I0,A)') 'Relaxing ', n_relax_, ' horizon(s) at fixed forcing ...'
                      call nek_log_message(msg, this_module, this_procedure)
                      call set_control_base(bf%g)
                      call clear_control_pert()
                      call brd_vec2nek(vx, vy, vz, pr, t, bf)
                      call setup_nonlinear_solver(recompute_dt = .true.,
     &                                            solve_temperature = .true.,
     &                                            cfl_limit = 0.4_dp,
     &                                            vtol = tol*0.1, ptol = tol*0.1)
                      do k = 1, n_relax_
                         call reset_qfft()
                         time = 0.0_dp
                         do istep = 1, nsteps
                            call nek_advance()
                            dtn = dt
                            call accumulate_qfft(get_flowrate_nek(), time, dtn)
                         end do
                         call extract_qfft(qbar)
                         write (msg, '(A,I3,A,*(1X,E14.6))') '   relax ', k, ': Qbar =', (qbar(i), i=1, n)
                         call nek_log_message(msg, this_module, this_procedure)
                      end do
                      call nek2brd_vec(bf, vx, vy, vz, pr, t)
                      bf%g = 0.0_dp
                      bf%g(1:n) = g0(1:n)
                   end if
       
             !
             ! ---- Control scaling.
             !      wg(i) = (||X|| / gscale(i))**2 makes a unit control perturbation
             !      contribute to the norm like a unit state perturbation. Without
             !      this the two blocks differ by many orders of magnitude and GMRES
             !      spends its whole budget on one of them.
             !
                   xnorm = bf%norm()
                   gsc = 0.0_dp
                   if (present(g_scale)) then
                      gsc(1:min(size(g_scale), n)) = g_scale(1:min(size(g_scale), n))
                   else
                      gsc(1) = abs(g0(1))
                      do i = 2, n
                         k = i/2
                         if (abs(g0(i)) > atol_dp) then
                            gsc(i) = abs(g0(i))
                         else
             ! Harmonic response rolls off with k (Womersley): a common scale for all
             ! harmonics is wrong by an order of magnitude by k = 3. This is a crude
             ! stand-in -- override via g_scale once you know the actual response.
                            gsc(i) = abs(g0(1))/real(1 + k, dp)
                         end if
                      end do
                   end if
                   call set_control_weights_auto(xnorm, gsc(1:n))
                   write (msg, '(3X,A,1X,E16.8)') padr('state norm:', pad), xnorm
                   call nek_log_message(msg, this_module, this_procedure)
                   call control_summary()
       
             !
             ! ---- Single Newton solve on the bordered system.
             !
                   call nek_log_message('Begin bordered Newton iteration', this_module, this_procedure)
                   sys%jacobian = nek_jacobian_bordered()
                   opts = newton_dp_opts(maxiter=maxiter_, ifbisect=.false.)
                   if (tol_mode_ == 1) then
                      call newton(sys, bf, gmres_rdp, info, atol=tol, options=opts, scheduler=nek_constant_tol)
                   else
                      call newton(sys, bf, gmres_rdp, info, atol=tol, options=opts, scheduler=nek_dynamic_tol)
                   end if
                   !call check_info(info, 'newton', this_module, this_procedure)
       
             !
             ! ---- Output
             !
                   call nek_log_message('OUTPUT:', this_module, this_procedure)
                   write (msg, '(3X,A,*(1X,E16.8))') padr('forcing g:', pad), (bf%g(i), i=1, n)
                   call nek_log_message(msg, this_module, this_procedure)
       
                   if (verify_) then
             ! Independent check: recompute the residual from scratch. The control
             ! rows are the flow-rate errors scaled by 1/A, i.e. bulk-velocity errors.
                      call sys%response(bf, res, tol)
                      qerr = res%g
                      write (msg, '(3X,A,1X,E16.8)') padr('|R| :', pad), res%norm()
                      call nek_log_message(msg, this_module, this_procedure)
                      write (msg, '(3X,A,*(1X,E16.8))') padr('dUbar err:', pad), (qerr(i), i=1, n)
                      call nek_log_message(msg, this_module, this_procedure)
                      write (msg, '(3X,A,*(1X,E16.8))') padr('Q err:', pad), (qerr(i)/get_res_scale(), i=1, n)
                      call nek_log_message(msg, this_module, this_procedure)
                      write (msg, '(3X,A,*(1X,E16.8))') padr('Q :', pad), (qerr(i)/get_res_scale() + Q_target(i), i=1, n)
                      call nek_log_message(msg, this_module, this_procedure)
                   end if
       
                   call set_fldindex('BFQ', 1)
                   call outpost_brd_dnek(bf, 'BFQ')
       
                end subroutine combined_flowrate_newton
      
         end module neklab_analysis_torus_2Dh
