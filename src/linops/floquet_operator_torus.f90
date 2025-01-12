      submodule(neklab_linops) floquet_operator_torus
         implicit none
      contains
      
         module procedure floquet_matvec
         select type (vec_in)
         type is (nek_dvector)
            select type (vec_out)
            type is (nek_dvector)
      ! Ensure correct nek status
               call setup_linear_solver(solve_baseflow = .false., 
     $                                  variable_dt    = .true.)
      ! Set the initial condition for Nek5000's linearized solver.
               call vec2nek(vxp, vyp, vzp, prp, tp, vec_in)
      ! Integrate the equations forward in time.
               time = 0.0_dp
               call pipe%reset_newton()         ! reset output counter to load baseflow files in order
               do istep = 1, pipe%get_nsteps()
                  call pipe%set_baseflow(vx, vy, vz, istep) ! sets the baseflow field and the appropriate timestep
                  call nek_advance()
               end do
      ! Extract the final solution to vector.
               call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)
            end select
         end select
         end procedure
      
         module procedure floquet_rmatvec
         select type (vec_in)
         type is (nek_dvector)
            select type (vec_out)
            type is (nek_dvector)
      ! Ensure correct nek status
               call setup_linear_solver(transpose      = .true.,
     $                                  solve_baseflow = .false., 
     $                                  variable_dt    = .true.)
      ! Set the initial condition for Nek5000's linearized solver.
               call vec2nek(vxp, vyp, vzp, prp, tp, vec_in)
      ! Integrate the equations forward in time.
               time = 0.0_dp
               call pipe%reset_newton()         ! reset output counter to load baseflow files in order
               do istep = 1, pipe%get_nsteps()
                  call pipe%set_baseflow(vx, vy, vz, istep) ! sets the baseflow field and the appropriate timestep
                  call nek_advance()
               end do
      ! Extract the final solution to vector.
               call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)
            end select
         end select
         end procedure
      end submodule
