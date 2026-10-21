      submodule(neklab_linops) exponential_propagator_2Dh
         implicit none
      contains
         module procedure init_exptA_2Dh_planar
         character(len=*), parameter :: this_procedure = "init_exptA_2Dh_planar"
         character(len=256) :: msg
         logical :: if_advW
         real(dp) :: wmax
         real, external :: glamax
         msg = 'call init_2Dh_planar() in userchk at istep = 0 before init.'
         if (.not. is_initialized_2Dh_planar()) call nek_stop_error(msg, this_module, this_procedure)
         ! For the baseflow field for dt/nsteps/cfl computation.
         call vec2nek(vx, vy, vz, pr, t, self%baseflow)
         call nek_log_message("Set self%baseflow -> vx, vy, vz, pr, t", this_module, this_procedure)
         ! Setup Nek5000 for perturbation solver.
         call setup_linear_solver(solve_baseflow = .false.,
     &                            endtime        = self%tau,
     &                            recompute_dt   = .true.,
     &                            cfl_limit      = 0.5_dp)
         ! Check whether we need to consider spanwise advection (W)
         if_advW = check_W(self%if_advW, self%betaz, wtol)
         msg = 'W advection is '//merge('active.  ','inactive.',if_advW)
         call nek_log_message(msg, this_module, this_procedure)
         if (if_advW) then
            wmax = glamax(t(1,1,1,1,1), lx1*ly1*lz1*nelv)
            write(msg,'(A,1pE12.4)') '    max|W| = ', wmax
            call nek_log_message(msg, this_module, this_procedure)
            write(msg,'(A,1pE12.4)') '     betaz = ', self%betaz
            call nek_log_message(msg, this_module, this_procedure)
         end if
         end procedure init_exptA_2Dh_planar

         module procedure exptA_2Dh_planar_matvec
         character(len=*), parameter :: this_procedure = "exptA_2Dh_planar_matvec"
         logical :: if_advW
         select type (vec_in)
         type is (nek_zvector)
            select type (vec_out)
            type is (nek_zvector)
         ! Force baseflow.
               call vec2nek(vx, vy, vz, pr, t, self%baseflow)
         ! Check whether we need to consider spanwise advection (W)
               if_advW = check_W(self%if_advW, self%betaz, wtol)
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
                  call nek_advance_2Dh_planar(self%betaz, if_advW)
               end do
         ! Copy the final solution to vector.
               call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)
            class default
               call type_error('vec_out','nek_zvector','OUT',this_module,this_procedure)
            end select
         class default
            call type_error('vec_in','nek_zvector','IN',this_module,this_procedure)
         end select
         end procedure exptA_2Dh_planar_matvec

         module procedure exptA_2Dh_planar_rmatvec
         character(len=*), parameter :: this_procedure = "exptA_2Dh_planar_rmatvecr"
         logical :: if_advW
         select type (vec_in)
         type is (nek_zvector)
            select type (vec_out)
            type is (nek_zvector)
         ! Force baseflow.
               call vec2nek(vx, vy, vz, pr, t, self%baseflow)
         ! Check whether we need to consider spanwise advection (W)
               if_advW = check_W(self%if_advW, self%betaz, wtol)
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
                  call nek_advance_2Dh_planar(-self%betaz, if_advW) ! sign flip for adjoint integration
               end do
         ! Copy the final solution to vector.
               call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)
            class default
               call type_error('vec_out','nek_zvector','OUT',this_module,this_procedure)
            end select
         class default
            call type_error('vec_in','nek_zvector','IN',this_module,this_procedure)
         end select
         end procedure exptA_2Dh_planar_rmatvec

         module procedure init_exptA_2Dh_axisym
         character(len=*), parameter :: this_procedure = "init_exptA_2Dh_axisym"
         real(dp), parameter :: diff_tol = 1.0e-08_dp
         character(len=256) :: msg
         ! Guard: the 2.5D (axisymmetric + swirl + heat) formulation is only
         ! well-posed if the run is actually set up as:
         !    - axisymmetric (IFAXIS)
         !    - azimuthal velocity enabled (IFAZIV)
         !    - heat transport is active with matching thermal and momentum 
         !      diffusivities (rhoCp == rho and conductivity == viscosity)
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

         logical function check_W(if_advW_in, betaz, wtol) result(if_advW)
            ! Cross-check whether W advection is enabled and whether the W field is nonzero.
            ! Rerun after each reset of the baseflow (W)
            implicit none
            logical, intent(in) :: if_advW_in
            real(dp), intent(in) :: betaz
            real(dp), intent(in) :: wtol
            ! internal
            character(len=*), parameter :: this_procedure = "check_W"
            character(len=128) :: msg
            real(dp) :: wmax
            real, external :: glamax
            wmax  = glamax(t(1,1,1,1,1), lx1*ly1*lz1*nelv)
            if_advW = if_advW_in
            if (wmax < wtol .and. if_advW_in) then
               msg = "W advection is enabled, but the W field is zero."
               call nek_stop_error(msg, this_module, this_procedure)
            else if (wmax >= wtol) then
               if (.not. if_advW_in) then
                  msg = "The W field is nonzero but W advection is disabled."
                  call nek_log_warning(msg, this_module, this_procedure)
               else if (betaz == 0.0_dp) then
                  if_advW = .false.
                  msg = "betaz is zero. W advection is disabled."
                  call nek_log_warning(msg, this_module, this_procedure)
               end if
            end if
         end function check_W
      end submodule exponential_propagator_2Dh
