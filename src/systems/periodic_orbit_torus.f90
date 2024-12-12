      submodule(neklab_systems) periodic_orbit_torus
         implicit none
      contains
         module procedure nonlinear_map_torus
      ! internal
         character(len=128) :: msg
         select type (vec_in)
         type is (nek_dvector)
            select type (vec_out)
            type is (nek_dvector)
      ! Set the initial condition
               call vec2nek(vx, vy, vz, pr, t, vec_in)
      ! Set appropriate tolerances and Nek status
               call setup_nonlinear_solver(recompute_dt=.true.,
     $   cfl_limit = 0.4_dp, vtol = atol/10.0, ptol = atol/10.0)
      ! Intgrate the nonlinear equations forward
               time = 0.0_dp
               do istep = 1, nsteps
                  call pipe%compute_bf_forcing(time) ! --> set neklab_forcing data
                  call nek_advance()
               end do
      ! Copy the final solution to vector.
               call nek2vec(vec_out, vx, vy, vz, pr, t)
      ! Evaluate residual F(X) - X.
               call vec_out%axpby(1.0_dp, vec_in, -1.0_dp)
            end select
         end select
         end procedure nonlinear_map_torus
      
         module procedure jac_direct_map_torus
      ! internal
         real(dp) :: atol
         select type (vec_in)
         type is (nek_dvector)
            select type (vec_out)
            type is (nek_dvector)
               atol = param(22)
      ! Set the baseflow initial condition
               call abs_vec2nek(vx, vy, vz, pr, t, self%X)
      ! Ensure correct nek status
               call setup_linear_solver(solve_baseflow=.true.,
     $   recompute_dt = .true., cfl_limit = 0.5_dp, 
     $   vtol = atol/2.0, ptol = atol/2.0)
      ! Set the initial condition for Nek5000's linearized solver.
               call vec2nek(vxp, vyp, vzp, prp, tp, vec_in)
      ! Integrate the equations forward in time.
               time = 0.0_dp
               do istep = 1, nsteps
                  call pipe%compute_bf_forcing(time) ! --> set neklab_forcing data
                  call nek_advance()
               end do
      ! Extract the final solution to vector.
               call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)
      ! Evaluate [ exp(tau*J) - I ] @ dx.
               call vec_out%axpby(1.0_dp, vec_in, -1.0_dp)
               param(22) = atol
            end select
         end select
         end procedure jac_direct_map_torus
      
         module procedure jac_adjoint_map_torus
      ! internal
         real(dp) :: atol
         select type (vec_in)
         type is (nek_dvector)
            select type (vec_out)
            type is (nek_dvector)
               atol = param(22)
      ! Set the baseflow initial condition
               call abs_vec2nek(vx, vy, vz, pr, t, self%X)
      ! Ensure correct nek status
               call setup_linear_solver(transpose=.true., solve_baseflow=.true.,
     $   recompute_dt = .true., cfl_limit = 0.5_dp, 
     $   vtol = atol/2.0, ptol = atol/2.0)
      ! Set the initial condition for Nek5000's linearized solver.
               call vec2nek(vxp, vyp, vzp, prp, tp, vec_in)
      ! Integrate the equations forward in time.
               time = 0.0_dp
               do istep = 1, nsteps
                  call pipe%compute_bf_forcing(time) ! --> set neklab_forcing data
                  call nek_advance()
               end do
      ! Extract the final solution to vector.
               call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)
      ! Evaluate [ exp(tau*J) - I ] @ dx.
               call vec_out%axpby(1.0_dp, vec_in, -1.0_dp)
               param(22) = atol
            end select
         end select
         end procedure jac_adjoint_map_torus
      end submodule
