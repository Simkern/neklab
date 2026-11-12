      submodule(neklab_systems) fixed_point_2Dh_axisym
         implicit none
         real(dp), parameter :: cfl_limit = 0.4_dp
         real(dp), parameter :: tol_fact  = 0.1_dp
         real(dp), parameter :: alpha     = 0.0_dp
      contains
         module procedure nonlinear_map_2Dh_axisym
         use neklab_2Dh, only: nek_advance_2Dh_axisym_nl
      ! internal
         character(len=*), parameter :: this_procedure = "nonlinear_map_2Dh_axisym"
         select type (vec_in)
         type is (nek_dvector)
            select type (vec_out)
            type is (nek_dvector)
      ! Set the initial condition
               call vec2nek(vx, vy, vz, pr, t, vec_in)
      ! Set appropriate tolerances and Nek status
               call setup_nonlinear_solver(recompute_dt      = .false., 
     &                                     solve_temperature = .true.,
     &                                     cfl_limit         = cfl_limit,
     &                                     vtol              = atol*tol_fact, 
     &                                     ptol              = atol*tol_fact)
      ! Intgrate the nonlinear equations forward
               time = 0.0_dp
               do istep = 1, nsteps
                  call nek_advance_2Dh_axisym_nl()
               end do
      ! Copy the final solution to vector.
               call nek2vec(vec_out, vx, vy, vz, pr, t)
      ! Evaluate residual F(X) - X.
               call vec_out%sub(vec_in)
            class default
               call nek_stop_error("The intent [OUT] argument 'vec_out' must be of type 'nek_dvector'",
     & this_module, this_procedure)
            end select
         class default
            call nek_stop_error("The intent [IN] argument 'vec_in' must be of type 'nek_dvector'",
     & this_module, this_procedure)
         end select
         end procedure nonlinear_map_2Dh_axisym
      
         module procedure jac_direct_map_2Dh_axisym
      ! internal
         character(len=*), parameter :: this_procedure = "jac_direct_map_2Dh_axisym"
         real(dp) :: atol
         select type (vec_in)
         type is (nek_dvector)
            select type (vec_out)
            type is (nek_dvector)
               atol = param(22)
      ! Set the baseflow initial condition
               call abs_vec2nek(vx, vy, vz, pr, t, self%X)
      ! Ensure correct nek status
               call setup_linear_solver(solve_baseflow    = .false.,
     &                                  solve_temperature = .true.,
     &                                  recompute_dt      = .false.,
     &                                  cfl_limit         = cfl_limit, 
     &                                  vtol              = atol*tol_fact,
     &                                  ptol              = atol*tol_fact)
      ! Set the initial condition for Nek5000's linearized solver.
               call vec2nek(vxp, vyp, vzp, prp, tp, vec_in)
      ! Integrate the equations forward in time.
               time = 0.0_dp
               do istep = 1, nsteps
                  call nek_advance_2Dh_axisym(alpha)
               end do
      ! Extract the final solution to vector.
               call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)
      ! Evaluate [ exp(tau*J) - I ] @ dx.
               call vec_out%sub(vec_in)
               param(22) = atol
            class default
               call nek_stop_error("The intent [OUT] argument 'vec_out' must be of type 'nek_dvector'",
     & this_module, this_procedure)
            end select
         class default
            call nek_stop_error("The intent [IN] argument 'vec_in' must be of type 'nek_dvector'",
     & this_module, this_procedure)
         end select
         end procedure jac_direct_map_2Dh_axisym
      
         module procedure jac_adjoint_map_2Dh_axisym
      ! internal
         character(len=*), parameter :: this_procedure = "jac_adjoint_map_2Dh_axisym"
         real(dp) :: atol
         select type (vec_in)
         type is (nek_dvector)
            select type (vec_out)
            type is (nek_dvector)
               atol = param(22)
      ! Set the baseflow initial condition
               call abs_vec2nek(vx, vy, vz, pr, t, self%X)
      ! Ensure correct nek status
               call setup_linear_solver(transpose         = .true., 
     &                                  solve_baseflow    = .false.,
     &                                  solve_temperature = .true.,
     &                                  recompute_dt      = .false.,
     &                                  cfl_limit         = cfl_limit, 
     &                                  vtol              = atol*tol_fact, 
     &                                  ptol              = atol*tol_fact)
      ! Set the initial condition for Nek5000's linearized solver.
               call vec2nek(vxp, vyp, vzp, prp, tp, vec_in)
      ! Integrate the equations forward in time.
               time = 0.0_dp
               do istep = 1, nsteps
                  call nek_advance_2Dh_axisym(alpha)
               end do
      ! Extract the final solution to vector.
               call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)
      ! Evaluate [ exp(tau*J) - I ] @ dx.
               call vec_out%sub(vec_in)
               param(22) = atol
            class default
               call nek_stop_error("The intent [OUT] argument 'vec_out' must be of type 'nek_dvector'",
     & this_module, this_procedure)
            end select
         class default
            call nek_stop_error("The intent [IN] argument 'vec_in' must be of type 'nek_dvector'",
     & this_module, this_procedure)
         end select
         end procedure jac_adjoint_map_2Dh_axisym
      end submodule