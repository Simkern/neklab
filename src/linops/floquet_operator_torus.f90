      submodule(neklab_linops) floquet_operator_torus
         implicit none
      contains
         module procedure floquet_init
      ! internal
         logical :: existfile
         character(len=128) :: msg
      ! Determine whether to solve for the baseflow or not
         inquire(file='f2dtorus001.fld', exist=existfile)
         if (existfile) then
            self%baseflow_computed = .true.
            msg = "Files f2dtorus001.fld exists. No not solve for baseflow"
            call logger%log_message(msg, module=this_module, procedure="init_floquet")
         else
            self%baseflow_computed = .false.
            msg = "Solve for baseflow at first iteration."
            call logger%log_message(msg, module=this_module, procedure="init_floquet")
      ! Set the baseflow for computation
            call vec2nek(vx, vy, vz, pr, t, self%baseflow)
         end if
         self%tau = pipe%get_period()
         self%is_initialized = .true.
         end procedure floquet_init
      
         module procedure floquet_matvec
         select type (vec_in)
         type is (nek_dvector)
            select type (vec_out)
            type is (nek_dvector)
               if (.not. self%is_initialized) call self%init()
               call pipe%set_2d_mode('floquet') ! reset output counter to load baseflow files in order
      ! Ensure correct nek status
               if (.not. self%baseflow_computed) then
                  call pipe%set_save_base(.true.)
                  call setup_linear_solver(solve_baseflow = .true.,
     $                                     endtime        = self%tau, 
     $                                     cfl_limit      = 0.4_dp,
     $                                     variable_dt    = .true.) ! -> solve for baseflow and save to f2dtorus***.fld
                  self%baseflow_computed = .true. ! we only need to do this once
               else
                  call pipe%set_save_base(.false.)
                  call setup_linear_solver(solve_baseflow = .false., 
     $                                     variable_dt    = .true.) ! -> load baseflow from 2d files
               end if
      ! Set the initial condition for Nek5000's linearized solver.
               call vec2nek(vxp, vyp, vzp, prp, tp, vec_in)
      ! Integrate the equations forward in time.
               time = 0.0_dp
               do istep = 1, pipe%get_nsteps()
                  call pipe%set_baseflow(vx, vy, vz, istep) ! sets the baseflow field and the appropriate timestep
                  call nek_advance()
               end do
      ! Extract the final solution to vector.
               call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)
            end select
         end select
         end procedure floquet_matvec
      
         module procedure floquet_rmatvec
         select type (vec_in)
         type is (nek_dvector)
            select type (vec_out)
            type is (nek_dvector)
               if (.not. self%is_initialized) call self%init()
               call pipe%set_2d_mode('floquet') ! reset output counter to load baseflow files in order
      ! Ensure correct nek status
               if (.not. self%baseflow_computed) then
                  call pipe%set_save_base(.true.)
                  call setup_linear_solver(transpose      = .true.,
     $                                     solve_baseflow = .true.,
     $                                     endtime        = self%tau, 
     $                                     cfl_limit      = 0.4_dp,
     $                                     variable_dt    = .true.) ! -> solve for baseflow and save to f2dtorus***.fld
                  self%baseflow_computed = .true. ! we only need to do this once
               else
                  call pipe%set_save_base(.false.)
                  call setup_linear_solver(transpose      = .true.,
     $                                     solve_baseflow = .false., 
     $                                     variable_dt    = .true.) ! -> load baseflow from 2d files
               end if
      ! Set the initial condition for Nek5000's linearized solver.
               call vec2nek(vxp, vyp, vzp, prp, tp, vec_in)
      ! Integrate the equations forward in time.
               time = 0.0_dp
               do istep = 1, pipe%get_nsteps()
                  call pipe%set_baseflow(vx, vy, vz, istep) ! sets the baseflow field and the appropriate timestep
                  call nek_advance()
               end do
      ! Extract the final solution to vector.
               call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)
            end select
         end select
         end procedure floquet_rmatvec
      end submodule
