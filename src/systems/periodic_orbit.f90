      submodule(neklab_systems) periodic_orbit
         implicit none
      contains
         module procedure nonlinear_map_UPO
      ! internal
         character(len=128) :: msg
         select type (vec_in)
         type is (nek_ext_dvector)
            select type (vec_out)
            type is (nek_ext_dvector)
      ! Set the initial condition
               call ext_vec2nek(vx, vy, vz, pr, t, vec_in)
      ! Set appropriate tolerances and Nek status
               call setup_nonlinear_solver(recompute_dt=.true., endtime=vec_in%T,
     $   cfl_limit = 0.4_dp, vtol = atol/10.0, ptol = atol/10.0)
               write (msg, '(A,F9.6)') 'Current period estimate, T = ', vec_in%T
               if (nid == 0) print *, msg
               call logger%log_message(msg, module=this_module, procedure='nonlinear_map_UPO')
      ! Intgrate the nonlinear equations forward
               time = 0.0_dp
               do istep = 1, nsteps
                  call nek_advance()
               end do
      ! Copy the final solution to vector.
               call nek2ext_vec(vec_out, vx, vy, vz, pr, t)
               vec_out%T = vec_in%T
      ! Evaluate residual F(X) - X.
               call vec_out%sub(vec_in)
            end select
         end select
         end procedure nonlinear_map_UPO
      
         module procedure jac_direct_map
      ! internal
         type(nek_ext_dvector) :: vec
         select type (vec_in)
         type is (nek_ext_dvector)
            select type (vec_out)
            type is (nek_ext_dvector)
      ! Set the baseflow initial condition
               call abs_ext_vec2nek(vx, vy, vz, pr, t, self%X)
      ! Ensure correct nek status -> set end time
               call setup_linear_solver(solve_baseflow=.true., transpose=.false.,
     $   recompute_dt = .true., endtime = get_period_abs(self%X), cfl_limit = 0.4_dp)
      ! Set the perturbation initial condition
               call ext_vec2nek(vxp, vyp, vzp, prp, tp, vec_in)
      ! Intgrate the coupled equations forward
               time = 0.0_dp
               do istep = 1, nsteps
                  call nek_advance()
               end do
      ! Copy the final solution to vector.
               call nek2ext_vec(vec_out, vxp, vyp, vzp, prp, tp)
      ! Evaluate [ exp(tau*J) - I ] @ dx.
               call vec_out%sub(vec_in)
      ! Evaluate f'(X(T), T) * dT and add it to the position residual
      ! Here we assume that vx,vy,vz contains the endpoint of the nonlinear trajectory
               call compute_fdot(vec)
               call vec_out%axpby(1.0_dp, vec, vec_in%T)
      ! Evaluate f'(X(0), 0).T @ dx and add phase condition
      ! Set the initial point of the nonlinear trajectory
               call abs_ext_vec2nek(vx, vy, vz, pr, t, self%X)
               call compute_fdot(vec)
               vec_out%T = vec_in%dot(vec)
            end select
         end select
         end procedure jac_direct_map
      
         module procedure jac_adjoint_map
      ! internal
         type(nek_ext_dvector) :: vec
         select type (vec_in)
         type is (nek_ext_dvector)
            select type (vec_out)
            type is (nek_ext_dvector)
      ! Set the baseflow initial condition
               call abs_ext_vec2nek(vx, vy, vz, pr, t, self%X)
      ! Ensure correct nek status -> set end time
               call setup_linear_solver(solve_baseflow=.true., transpose=.true.,
     $   recompute_dt = .true., endtime = get_period_abs(self%X), cfl_limit = 0.4_dp)
      ! Set the perturbation initial condition
               call ext_vec2nek(vxp, vyp, vzp, prp, tp, vec_in)
      ! Integrate the equations forward in time.
               time = 0.0_dp
               do istep = 1, nsteps
                  call nek_advance()
               end do
      ! Copy the final solution to vector.
               call nek2ext_vec(vec_out, vxp, vyp, vzp, prp, tp)
      ! Evaluate [ exp(tau*J) - I ] @ dx.
               call vec_out%sub(vec_in)
      ! Evaluate f'(X(T), T) * dT and add it to the position residual
               call compute_fdot(vec)
               call vec_out%axpby(1.0_dp, vec, vec_in%T)
      ! Evaluate f'(X(0), 0).T @ dx and add phase condition
      ! Set the initial point of the orbit
               call abs_ext_vec2nek(vx, vy, vz, pr, t, self%X)
               call compute_fdot(vec)
               vec_out%T = vec_in%dot(vec)
            end select
         end select
         end procedure jac_adjoint_map
      end submodule
