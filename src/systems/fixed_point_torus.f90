      submodule(neklab_systems) fixed_point_torus
         implicit none
      contains
         module procedure nonlinear_map_torus
         select type (vec_in)
         type is (nek_dvector)
            select type (vec_out)
            type is (nek_dvector)
         ! Set the initial condition
               call vec2nek(vx, vy, vz, pr, t, vec_in)
         ! Set appropriate tolerances
               call setup_nonlinear_solver(recompute_dt=.true., cfl_limit=0.5_dp,
     $   vtol = atol/10.0, ptol = atol/10.0)
         ! Intgrate the nonlinear equations forward
               time = 0.0_dp
               do istep = 1, nsteps
                  call pipe%compute_bf_forcing(time) ! --> set neklab_forcing data
                  call nek_advance()
               end do
         ! Extract the final solution to vector.
               call nek2vec(vec_out, vx, vy, vz, pr, t)
         ! Evaluate residual F(X) - X.
               call vec_out%sub(vec_in)
            end select
         end select
         end procedure nonlinear_map_torus

      end submodule