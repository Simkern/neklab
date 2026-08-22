      submodule(neklab_linops) exponential_propagator_2Dh
         implicit none
      contains
         module procedure init_exptA_2Dh
         character(len=*), parameter :: this_procedure = "init_exptA_2Dh"
         ! For the baseflow field for dt/nsteps/cfl computation.
         call vec2nek(vx, vy, vz, pr, t, self%baseflow)
         call nek_log_message("Set self%baseflow -> vx, vy, vz, pr, t", this_module, this_procedure)
         ! Setup Nek5000 for perturbation solver.
         call setup_linear_solver(solve_baseflow = .false.,
     &                            endtime        = self%tau,
     &                            recompute_dt   = .true.,
     &                            cfl_limit      = 0.5_dp)
         end procedure init_exptA_2Dh

         module procedure exptA_2Dh_matvec
         character(len=*), parameter :: this_procedure = "exptA_2Dh_matvec"
         select type (vec_in)
         type is (nek_zvector)
            select type (vec_out)
            type is (nek_zvector)
         ! Force baseflow.
               call vec2nek(vx, vy, vz, pr, t, self%baseflow)
         ! Ensure correct Nek status
               call setup_linear_solver(transpose     = .false.,
     &                                  silent        = .false.,
     &                                  endtime       = self%tau,
     &                                  recompute_dt  = .true.,
     &                                  cfl_limit     = 0.5_dp)
         ! Set initial condition for the linearized solver.
               call vec2nek(vxp, vyp, vzp, prp, tp, vec_in)
         ! Integrate the equations forward in time.
               time = 0.0_dp
               do istep = 1, nsteps
                  call nek_advance_2Dh(self%betaz)
               end do
         ! Copy the final solution to vector.
               call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)
            class default
               call type_error('vec_out','nek_zvector','OUT',this_module,this_procedure)
            end select
         class default
            call type_error('vec_in','nek_zvector','IN',this_module,this_procedure)
         end select
         end procedure exptA_2Dh_matvec

         module procedure exptA_2Dh_rmatvec
         character(len=*), parameter :: this_procedure = "exptA_2Dh_rmatvec"
         select type (vec_in)
         type is (nek_zvector)
            select type (vec_out)
            type is (nek_zvector)
         ! Force baseflow.
               call vec2nek(vx, vy, vz, pr, t, self%baseflow)
         ! Ensure correct Nek status
               call setup_linear_solver(transpose     = .true.,
     &                                  silent        = .false.,
     &                                  endtime       = self%tau,
     &                                  recompute_dt  = .true.,
     &                                  cfl_limit     = 0.5_dp)
         ! Set initial condition for the linearized solver.
               call vec2nek(vxp, vyp, vzp, prp, tp, vec_in)
         ! Integrate the equations forward in time.
               time = 0.0_dp
               do istep = 1, nsteps
                  call nek_advance_2Dh(-self%betaz) ! sign flip for adjoint integration
               end do
         ! Copy the final solution to vector.
               call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)
            class default
               call type_error('vec_out','nek_zvector','OUT',this_module,this_procedure)
            end select
         class default
            call type_error('vec_in','nek_zvector','IN',this_module,this_procedure)
         end select
         end procedure exptA_2Dh_rmatvec

         module procedure init_exptA_2Dh_axisym
         character(len=*), parameter :: this_procedure = "init_exptA_2Dh_axisym"
         ! For the baseflow field for dt/nsteps/cfl computation.
         call vec2nek(vx, vy, vz, pr, t, self%baseflow)
         call nek_log_message("Set self%baseflow -> vx, vy, vz, pr, t", this_module, this_procedure)
         ! Setup Nek5000 for perturbation solver.
         call setup_linear_solver(solve_baseflow = .false.,
     &                            endtime        = self%tau,
     &                            recompute_dt   = .true.,
     &                            cfl_limit      = 0.5_dp)
         end procedure init_exptA_2Dh_axisym

         module procedure exptA_2Dh_axisym_matvec
         character(len=*), parameter :: this_procedure = "exptA_2Dh_axisym_matvec"
         select type (vec_in)
         type is (nek_zvector)
            select type (vec_out)
            type is (nek_zvector)
         ! Force baseflow.
               call vec2nek(vx, vy, vz, pr, t, self%baseflow)
         ! Ensure correct Nek status
               call setup_linear_solver(transpose     = .false.,
     &                                  silent        = .false.,
     &                                  endtime       = self%tau,
     &                                  recompute_dt  = .true.,
     &                                  cfl_limit     = 0.5_dp)
         ! Set initial condition for the linearized solver.
               call vec2nek(vxp, vyp, vzp, prp, tp, vec_in)
         ! Integrate the equations forward in time.
               time = 0.0_dp
               do istep = 1, nsteps
                  call nek_advance_2Dh_axisym(self%alpha_s)
               end do
         ! Copy the final solution to vector.
               call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)
            class default
               call type_error('vec_out','nek_zvector','OUT',this_module,this_procedure)
            end select
         class default
            call type_error('vec_in','nek_zvector','IN',this_module,this_procedure)
         end select
         end procedure exptA_2Dh_axisym_matvec

         module procedure exptA_2Dh_axisym_rmatvec
         character(len=*), parameter :: this_procedure = "exptA_2Dh_axisym_rmatvec"
         select type (vec_in)
         type is (nek_zvector)
            select type (vec_out)
            type is (nek_zvector)
         ! Force baseflow.
               call vec2nek(vx, vy, vz, pr, t, self%baseflow)
         ! Ensure correct Nek status
               call setup_linear_solver(transpose     = .true.,
     &                                  silent        = .false.,
     &                                  endtime       = self%tau,
     &                                  recompute_dt  = .true.,
     &                                  cfl_limit     = 0.5_dp)
         ! Set initial condition for the linearized solver.
               call vec2nek(vxp, vyp, vzp, prp, tp, vec_in)
         ! Integrate the equations forward in time.
               time = 0.0_dp
               do istep = 1, nsteps
                  call nek_advance_2Dh_axisym(-self%alpha_s) ! sign flip for adjoint integration
               end do
         ! Copy the final solution to vector.
               call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)
            class default
               call type_error('vec_out','nek_zvector','OUT',this_module,this_procedure)
            end select
         class default
            call type_error('vec_in','nek_zvector','IN',this_module,this_procedure)
         end select
         end procedure exptA_2Dh_axisym_rmatvec

      end submodule exponential_propagator_2Dh
