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
         real(dp), parameter :: diff_tol = 1.0e-08_dp
         character(len=256) :: msg
         ! Guard: the 2.5D (axisymmetric + swirl + heat) formulation is only
         ! well-posed if the run is actually set up as:
         !    - axisymmetric (IFAXIS)
         !    - azimuthal velocity enabled (IFAZIV)
         !    - heat transport is active with matching thermal and momentum diffusivities (rhoCp == rho and conductivity == viscosity)
         if (.not. ifaxis) call nek_stop_error("2Dh axisym exptA requires ifaxis = .true.", this_module, this_procedure)
         if (.not. ifaziv) call nek_stop_error("2Dh axisym exptA requires ifaziv = .true.", this_module, this_procedure)
         if (.not. ifheat) call nek_stop_error("2Dh axisym exptA requires ifheat = .true.", this_module, this_procedure)
         if (abs(cpfld(1,1) - cpfld(2,1)) > diff_tol) then
            write(msg,'(A,E18.7)') '           viscosity = ', cpfld(1,1)
            call nek_log_message(msg, this_module, this_procedure)
            write(msg,'(A,E18.7)') 'thermal conductivity = ', cpfld(2,1)
            call nek_log_message(msg, this_module, this_procedure)
            write(msg,'(A,E18.7)') '          difference = ', abs(cpfld(2,1) - cpfld(1,1))
            call nek_log_message(msg, this_module, this_procedure)
            call nek_stop_error("2Dh axisym exptA requires viscosity = conductivity.", this_module, this_procedure)
         end if
         if (abs(cpfld(1,2) - cpfld(2,2)) > diff_tol) then
            write(msg,'(A,E18.7)') '       rho = ', cpfld(1,1)
            call nek_log_message(msg, this_module, this_procedure)
            write(msg,'(A,E18.7)') '     rhoCp = ', cpfld(1,2)
            call nek_log_message(msg, this_module, this_procedure)
            write(msg,'(A,E18.7)') 'difference = ', abs(cpfld(1,1) - cpfld(1,2))
            call nek_log_message(msg, this_module, this_procedure)
            call nek_stop_error("2Dh axisym exptA requires rho = rhoCp.", this_module, this_procedure)
         end if
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
                  call nek_advance_2Dh_axisym(self%alphas)
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
                  call nek_advance_2Dh_axisym(-self%alphas) ! sign flip for adjoint integration
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
