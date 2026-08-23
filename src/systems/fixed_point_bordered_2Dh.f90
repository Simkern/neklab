      !--------------------------------------------------------------------------
      !
      !  The residual is
      !
      !     R(X, g) = [ Phi_tau(X; g) - X              ]
      !               [ wres * ( Qbar_k(X; g) - Qt_k ) ]   k = 1 .. nctrl
      !
      !  where Qbar_k is the k-th temporal Fourier coefficient of the flow rate
      !  accumulated ALONG the trajectory. For K = 0 that is the mean of Q over
      !  the horizon, which equals Q once X is a fixed point.
      !
      !  The Jacobian matvec on (dX, dg) is ONE linearised integration: dX is the
      !  initial condition, dg drives a source in the tp equation, and the border
      !  row is accumulated from the same trajectory. Note that the (2,2) block is
      !  therefore NOT zero -- it is the direct sensitivity of Qbar to dg through
      !  the trajectory, which is exactly what makes the border well scaled.

      submodule(neklab_systems) fixed_point_bordered_2Dh
         implicit none
      contains

      !====================================================================
      !     NONLINEAR RESIDUAL
      !====================================================================

         module procedure nonlinear_map_bordered
      ! internal
         character(len=*), parameter :: this_procedure = 'nonlinear_map_bordered'
         real(dp), dimension(lg) :: qbar, qt
         real(dp) :: dtn
         integer :: i
         select type (vec_in)
         type is (nek_bordered_dvector)
            select type (vec_out)
            type is (nek_bordered_dvector)
      ! Control: the forcing is part of the unknown vector.
               call set_control_base(vec_in%g)
               call clear_control_pert()
      ! Set the initial condition (state block only).
               call brd_vec2nek(vx, vy, vz, pr, t, vec_in)
      ! Set appropriate tolerances and Nek status
               call setup_nonlinear_solver(recompute_dt = .true.,
     &                                     solve_temperature = .true.,
     &                                     cfl_limit    = 0.4_dp,
     &                                     vtol         = atol*0.1,
     &                                     ptol         = atol*0.1)
      ! Integrate the nonlinear equations forward, accumulating the flow rate.
               call reset_qfft()
               time = 0.0_dp
               do istep = 1, nsteps
                  call nek_advance()
                  dtn = dt
                  call accumulate_qfft(get_flowrate_nek(), time, dtn)
               end do
      ! State block: F(X) - X.
               call nek2brd_vec(vec_out, vx, vy, vz, pr, t)
               call vec_out%sub(vec_in)
      ! Control block: OVERWRITE (sub above left -g there).
               call extract_qfft(qbar)
               qt = get_control_target()
               vec_out%g = 0.0_dp
               do i = 1, nctrl
                  vec_out%g(i) = get_res_scale()*(qbar(i) - qt(i))
               end do
            class default
               call nek_stop_error("The intent [OUT] argument 'vec_out' must be of type 'nek_bordered_dvector'",
     & this_module, this_procedure)
            end select
         class default
            call nek_stop_error("The intent [IN] argument 'vec_in' must be of type 'nek_bordered_dvector'",
     & this_module, this_procedure)
         end select
         end procedure nonlinear_map_bordered

      !====================================================================
      !     JACOBIAN -- DIRECT
      !====================================================================

         module procedure jac_direct_map_bordered
      ! internal
         character(len=*), parameter :: this_procedure = 'jac_direct_map_bordered'
         real(dp) :: atol
         real(dp), dimension(lg) :: qbar
         real(dp) :: dtn
         integer :: i
         select type (vec_in)
         type is (nek_bordered_dvector)
            select type (vec_out)
            type is (nek_bordered_dvector)
               atol = param(22)
      ! Base state and base forcing come from the current Newton iterate.
               call abs_brd_vec2nek(vx, vy, vz, pr, t, self%X)
               call set_control_base(get_control_abs(self%X))
      ! Forcing perturbation: drives a source in the tp equation via userq(jp>0).
      ! This is what makes the whole border cost ONE integration instead of
      ! one per control direction -- the dg directions superpose inside it.
               call set_control_pert(vec_in%g)
      ! Ensure correct nek status
               call setup_linear_solver(solve_baseflow = .false.,
     &                                  solve_temperature = .true.,
     &                                  recompute_dt   = .true.,
     &                                  cfl_limit      = 0.4_dp,
     &                                  vtol           = atol*0.5,
     &                                  ptol           = atol*0.5)
      ! Initial condition for Nek5000's linearised solver.
               call brd_vec2nek(vxp, vyp, vzp, prp, tp, vec_in)
      ! Integrate, accumulating the linearised flow rate along the trajectory.
               call reset_qfft()
               time = 0.0_dp
               do istep = 1, nsteps
                  call nek_advance_2Dh_axisym(0.0_dp)
                  dtn = dt
                  call accumulate_qfft(get_flowrate_pert_nek(1), time, dtn)
               end do
      ! State block: [ M_tau - I ] dX + S_tau dg
               call nek2brd_vec(vec_out, vxp, vyp, vzp, prp, tp)
               call vec_out%sub(vec_in)
      ! Control block: the linearised constraint rows.
               call extract_qfft(qbar)
               vec_out%g = 0.0_dp
               do i = 1, nctrl
                  vec_out%g(i) = get_res_scale()*qbar(i)
               end do
      ! MANDATORY: leave no perturbation forcing behind, or the next nonlinear
      ! solve inherits a spurious source.
               call clear_control_pert()
               param(22) = atol
            class default
               call nek_stop_error("The intent [OUT] argument 'vec_out' must be of type 'nek_bordered_dvector'",
     & this_module, this_procedure)
            end select
         class default
            call nek_stop_error("The intent [IN] argument 'vec_in' must be of type 'nek_bordered_dvector'",
     & this_module, this_procedure)
         end select
         end procedure jac_direct_map_bordered

      !====================================================================
      !     JACOBIAN -- ADJOINT
      !====================================================================

         module procedure jac_adjoint_map_bordered
      ! internal
         character(len=*), parameter :: this_procedure = 'jac_adjoint_map_bordered'
      ! GMRES never calls this. Transposing the border correctly requires the
      ! adjoint of the forcing injection (i.e. the time-integral of the adjoint
      ! field against each basis function), which is a separate piece of work.
      ! Failing loudly is better than returning a silently wrong operator.
         call nek_stop_error('rmatvec is not implemented for the bordered system.'//
     &                       ' Use a non-transpose solver (gmres).', this_module, this_procedure)
         end procedure jac_adjoint_map_bordered

      end submodule fixed_point_bordered_2Dh