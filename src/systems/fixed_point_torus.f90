      submodule(neklab_systems) fixed_point_torus
         implicit none
      contains
         module procedure nonlinear_map_torus
      ! internal
         real(dp) :: pd
         select type (vec_in)
         type is (nek_dvector)
            select type (vec_out)
            type is (nek_dvector)
               pd = param(10)
      ! Set the initial condition
               call vec2nek(vx, vy, vz, pr, t, vec_in)
      ! Set appropriate tolerances and Nek status
               call setup_nonlinear_solver(recompute_dt = .true., 
     &                                     endtime      = pd,
     &                                     cfl_limit    = 0.4_dp,
     &                                     vtol         = atol*0.1, 
     &                                     ptol         = atol*0.1)
      ! Intgrate the nonlinear equations forward
               time = 0.0_dp
               call pipe%reset_mflow_fft()
               do istep = 1, nsteps
                  call pipe%compute_bf_forcing(time) ! --> set neklab_forcing data
                  call nek_advance()
                  call pipe%compute_mflow_fft(period=pd) ! integrate Fourier coefficients
               end do
      ! Record the mass flow rate and number of timesteps per period for subsequent linear runs
               call pipe%extract_mflow_fft(period=pd)
      ! Copy the final solution to vector.
               call nek2vec(vec_out, vx, vy, vz, pr, t)
      ! Evaluate residual F(X) - X.
               call vec_out%sub(vec_in)
            class default
               call nek_stop_error("The intent [OUT] argument 'vec_out' must be of type 'nek_dvector'",
     & this_module, 'nonlinear_map_torus')
            end select
         class default
            call nek_stop_error("The intent [IN] argument 'vec_in' must be of type 'nek_dvector'",
     & this_module, 'nonlinear_map_torus')
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
               call setup_linear_solver(solve_baseflow = .false.,
     &                                  recompute_dt   = .true.,
     &                                  cfl_limit      = 0.5_dp, 
     &                                  vtol           = atol*0.5,
     &                                  ptol           = atol*0.5)
      ! Set the initial condition for Nek5000's linearized solver.
               call vec2nek(vxp, vyp, vzp, prp, tp, vec_in)
      ! Integrate the equations forward in time.
               time = 0.0_dp
               do istep = 1, nsteps
                  call nek_advance()
               end do
      ! Extract the final solution to vector.
               call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)
      ! Evaluate [ exp(tau*J) - I ] @ dx.
               call vec_out%sub(vec_in)
               param(22) = atol
            class default
               call nek_stop_error("The intent [OUT] argument 'vec_out' must be of type 'nek_dvector'",
     & this_module, 'jac_direct_map_torus')
            end select
         class default
            call nek_stop_error("The intent [IN] argument 'vec_in' must be of type 'nek_dvector'",
     & this_module, 'jac_direct_map_torus')
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
               call setup_linear_solver(transpose      = .true., 
     &                                  solve_baseflow = .false.,
     &                                  recompute_dt   = .true.,
     &                                  cfl_limit      = 0.5_dp, 
     &                                  vtol           = atol*0.5, 
     &                                  ptol           = atol*0.5)
      ! Set the initial condition for Nek5000's linearized solver.
               call vec2nek(vxp, vyp, vzp, prp, tp, vec_in)
      ! Integrate the equations forward in time.
               time = 0.0_dp
               do istep = 1, nsteps
                  call nek_advance()
               end do
      ! Extract the final solution to vector.
               call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)
      ! Evaluate [ exp(tau*J) - I ] @ dx.
               call vec_out%sub(vec_in)
               param(22) = atol
            class default
               call nek_stop_error("The intent [OUT] argument 'vec_out' must be of type 'nek_dvector'",
     & this_module, 'jac_adjoint_map_torus')
            end select
         class default
            call nek_stop_error("The intent [IN] argument 'vec_in' must be of type 'nek_dvector'",
     & this_module, 'jac_adjoint_map_torus')
         end select
         end procedure jac_adjoint_map_torus
      end submodule