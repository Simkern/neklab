      submodule(neklab_linops) exponential_propagator
         implicit none
      contains
         module procedure init_exptA
      ! For the baseflow field for dt/nsteps/cfl computation.
         call vec2nek(vx, vy, vz, pr, t, self%baseflow)
         call nek_log_information("Set self%baseflow -> vx, vy, vz, pr, t", this_module, "init_exptA")
      ! Setup Nek5000 for perturbation solver.
         call setup_linear_solver(solve_baseflow = .false.,
     &                            endtime        = self%tau,
     &                            recompute_dt   = .true.,
     &                            cfl_limit      = 0.5_dp)
         end procedure

         module procedure exptA_matvec
         character(len=*), parameter :: this_procedure = 'exptA_matvec'
         integer :: nrst
         select type (vec_in)
         type is (nek_dvector)
            select type (vec_out)
            type is (nek_dvector)

               nrst = abs(param(27)) - 1
      ! Set baseflow.
               call vec2nek(vx, vy, vz, pr, t, self%baseflow)
      
      ! Set nek configuration (after the v[xzy] and v[xyz]p fields are updated)   
               call setup_linear_solver(transpose    = .false.,
     &                                  silent       = .true.,
     &                                  endtime      = self%tau,
     &                                  recompute_dt = .true.,
     &                                  cfl_limit    = 0.5_dp)

      ! Set initial condition for the linearized solver.
               call vec2nek(vxp, vyp, vzp, prp, tp, vec_in)
      
      ! Integrate the equations forward in time.
               time = 0.0_dp
               do istep = 1, nsteps

                  call nek_advance()
                  
                  ! Set restart fields if present.
                  if (istep <= nrst) call get_rst_dnek(vec_in, istep)

               end do
      
      ! Copy the final solution to vector.
               call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)
      
      ! Compute restart fields.
               call compute_rst_dnek(vec_out, nrst)

            class default
               call type_error('vec_out','nek_dvector','OUT',this_module, this_procedure)
            end select
         class default
            call type_error('vec_in','nek_dvector','IN',this_module, this_procedure)
         end select
         end procedure

         module procedure exptA_rmatvec
         character(len=*), parameter :: this_procedure = 'exptA_rmatvec'
         integer :: nrst
         select type (vec_in)
         type is (nek_dvector)
            select type (vec_out)
            type is (nek_dvector)
            
               nrst = abs(param(27)) - 1
      ! Set baseflow.
               call vec2nek(vx, vy, vz, pr, t, self%baseflow)
      
      ! Set nek configuration (after the v[xzy] and v[xyz]p fields are updated)
               call setup_linear_solver(transpose    = .true., 
     &                                  silent       = .true.,  
     &                                  endtime      = self%tau,  
     &                                  recompute_dt = .true.,  
     &                                  cfl_limit    = 0.5_dp)

      ! Set initial condition for the linearized solver.
               call vec2nek(vxp, vyp, vzp, prp, tp, vec_in)
      
      ! Integrate the equations forward in time.
               time = 0.0_dp
               do istep = 1, nsteps

                  call nek_advance()
                  
                  ! Set restart fields if present.
                  if (istep <= nrst) call get_rst_dnek(vec_in, istep)

               end do
      
      ! Copy the final solution to vector.
               call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)
      
      ! Compute restart fields.
               call compute_rst_dnek(vec_out, nrst)
       
            class default
               call type_error('vec_out','nek_dvector','OUT',this_module, this_procedure)
            end select
         class default
            call type_error('vec_in','nek_dvector','IN',this_module, this_procedure)
         end select
         end procedure
      end submodule
