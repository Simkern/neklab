      submodule(neklab_linops) floquet_operator_torus
         implicit none
      contains
         module procedure floquet_init
         character(len=*), parameter :: this_procedure = 'floquet_init'
         integer :: idx, ns
         integer :: nsaver, wdsizr, ierr
         logical :: existfile
         character(len=128) :: msg
         character(len=132) :: fname
      ! Determine whether to solve for the baseflow or not
         ns = 0
         idx = 1
         write(fname,'("f2dtorus",I3.3,".fld")') idx
         inquire(file=fname, exist=existfile)
         if (existfile) then
            self%baseflow_computed = .true.
            msg = "2d baseflow files exist. No not solve for baseflow"
            call nek_log_message(msg, this_module, this_procedure)
            do while (existfile)
               ! read first file and get nsteps
               call pipe%get_nsteps_from_header(fname, nsaver)
               ns = ns + nsaver              
               idx = idx + 1
               write(fname,'("f2dtorus",I3.3,".fld")') idx
               inquire(file=fname, exist=existfile)
            end do
            call bcast(ns, isize)          ! broadcast number of saved snapshots
            call pipe%set_nsteps(ns)
            write(msg,'(A,I0,A,I0,A)') 'Found ', idx-1, ' baseflow files: ', ns, ' timesteps per period.'
            call nek_log_message(msg, this_module, this_procedure)
         else
            self%baseflow_computed = .false.
            msg = "No 2d bseflow files found. Solve for baseflow at first iteration."
            call nek_log_message(msg, this_module, this_procedure)
      ! Set the baseflow for computation
            call vec2nek(vx, vy, vz, pr, t, self%baseflow)
         end if
         self%tau = pipe%get_period()
         self%is_initialized = .true.
         end procedure floquet_init
      
         module procedure floquet_matvec
         character(len=*), parameter :: this_procedure = 'floquet_matvec'
         integer :: nrst
         select type (vec_in)
         type is (nek_dvector)
            select type (vec_out)
            type is (nek_dvector)

               if (.not. self%is_initialized) call self%init()
               ! reset output counter to load baseflow files in order
               call pipe%set_2d_mode('floquet')

               nrst = abs(param(27)) - 1

      ! Ensure correct nek status
               if (.not. self%baseflow_computed) then
                  call pipe%set_save_base(.true.)
                  call setup_linear_solver(solve_baseflow = .true.,
     &                                     endtime        = self%tau, 
     &                                     cfl_limit      = 0.4_dp,
     &                                     variable_dt    = .true.) ! -> solve for baseflow and save to f2dtorus***.fld
               else
                  call pipe%set_save_base(.false.)
                  call setup_linear_solver(solve_baseflow = .false., 
     &                                     endtime        = self%tau, 
     &                                     variable_dt    = .true.) ! -> load baseflow from 2d files
               end if

      ! Set the initial condition for Nek5000's linearized solver.
               call vec2nek(vxp, vyp, vzp, prp, tp, vec_in)

      ! Integrate the equations forward in time.
               time = 0.0_dp
               if (.not. self%baseflow_computed) then
                  istep = 0
                  do while (lastep == 0)
                     istep = istep + 1
                     call pipe%compute_bf_forcing(time) ! --> set neklab_forcing data
                     call pipe%save_2d_fields(vx,vy,vz)
                     call nek_advance()
                     ! Set restart fields if present.
                     if (istep <= nrst) call self%get_rst(vec_in, istep)
                  end do
                  call pipe%outpost_2d_fields()
                  ! Record the mass flow rate and number of timesteps per period for subsequent linear runs
                  call pipe%set_nsteps(istep)
                  self%baseflow_computed = .true. ! we only need to do this once
               else
                  do istep = 1, pipe%get_nsteps()
                     ! sets the baseflow field and the appropriate timestep
                     call pipe%set_baseflow(vx, vy, vz, istep)
                     call nek_advance()
                     ! Set restart fields if present.
                     if (istep <= nrst) call self%get_rst(vec_in, istep)
                  end do
               end if

      ! Extract the final solution to vector.
               call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)

      ! Compute restart fields.
               call self%compute_rst(vec_out, nrst)

               self%baseflow_computed = .true. ! we only need to do this once
            
            class default
               call type_error('vec_out','nek_dvector','OUT',this_module, this_procedure)
            end select
         class default
            call type_error('vec_in','nek_dvector','IN',this_module, this_procedure)
         end select
         end procedure floquet_matvec
      
         module procedure floquet_rmatvec
         character(len=*), parameter :: this_procedure = 'floquet_rmatvec'
         integer :: nrst
         select type (vec_in)
         type is (nek_dvector)
            select type (vec_out)
            type is (nek_dvector)

               if (.not. self%is_initialized) call self%init()
               ! reset output counter to load baseflow files in order
               call pipe%set_2d_mode('floquet')

               nrst = abs(param(27)) - 1

      ! Ensure correct nek status
               if (.not. self%baseflow_computed) then
                  call pipe%set_save_base(.true.)
                  call setup_linear_solver(transpose      = .true.,
     &                                     solve_baseflow = .true.,
     &                                     endtime        = self%tau, 
     &                                     cfl_limit      = 0.4_dp,
     &                                     variable_dt    = .true.) ! -> solve for baseflow and save to f2dtorus***.fld
                  self%baseflow_computed = .true. ! we only need to do this once
               else
                  call pipe%set_save_base(.false.)
                  call setup_linear_solver(transpose      = .true.,
     &                                     solve_baseflow = .false.,
     &                                     endtime        = self%tau,  
     &                                     variable_dt    = .true.) ! -> load baseflow from 2d files
               end if

      ! Set the initial condition for Nek5000's linearized solver.
               call vec2nek(vxp, vyp, vzp, prp, tp, vec_in)

      ! Integrate the equations forward in time.
               time = 0.0_dp
               if (.not. self%baseflow_computed) then
                  istep = 0
                  do while (lastep == 0)
                     istep = istep + 1
                     call pipe%compute_bf_forcing(time) ! --> set neklab_forcing data
                     call pipe%save_2d_fields(vx,vy,vz)
                     call nek_advance()
                     ! Set restart fields if present.
                     if (istep <= nrst) call self%get_rst(vec_in, istep)
                  end do
                  call pipe%outpost_2d_fields()
                  ! Record the mass flow rate and number of timesteps per period for subsequent linear runs
                  call pipe%set_nsteps(istep)
                  self%baseflow_computed = .true. ! we only need to do this once
               else
                  do istep = 1, pipe%get_nsteps()
                     ! sets the baseflow field and the appropriate timestep
                     call pipe%set_baseflow(vx, vy, vz, istep)
                     call nek_advance()
                     ! Set restart fields if present.
                     if (istep <= nrst) call self%get_rst(vec_in, istep)
                  end do
               end if

      ! Extract the final solution to vector.
               call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)
      
      ! Compute restart fields.
               call self%compute_rst(vec_out, nrst)

               self%baseflow_computed = .true. ! we only need to do this once

            class default
               call type_error('vec_out','nek_dvector','OUT',this_module, this_procedure)
            end select
         class default
            call type_error('vec_in','nek_dvector','IN',this_module, this_procedure)
         end select
         end procedure floquet_rmatvec

         module procedure floquet_compute_rst
            character(len=*), parameter :: this_procedure = 'floquet_compute_rst'
            type(nek_dvector) :: vec_rst
            character(len=128) :: msg
            integer :: nstp
            select type(vec_out)
            type is (nek_dvector)
               nstp = pipe%get_nsteps()
               write(msg,'(A,I0,A)') 'Run ', nrst, ' extra step(s) to fill up restart arrays.'
               call nek_log_information(msg, this_module, this_procedure)
               ! reset output counter to load baseflow files in order
               call pipe%set_2d_mode('floquet')
               fintim = fintim + nrst*dt
               do istep = nstp + 1, nstp + nrst
                  ! sets the baseflow field and the appropriate timestep
                  call pipe%set_baseflow(vx, vy, vz, istep - nstp)
                  call nek_advance()
                  call nek2vec(vec_rst, vxp, vyp, vzp, prp, tp)
                  call vec_out%save_rst(vec_rst, istep - nstp)
               end do
            class default
               call type_error('vec_out','nek_dvector','OUT',this_module, this_procedure)
            end select
         end procedure floquet_compute_rst

         module procedure floquet_get_rst
            character(len=*), parameter :: this_procedure = 'floquet_get_rst'
            type(nek_dvector) :: vec_rst
            character(len=128) :: msg
            select type(vec_in)
            type is (nek_dvector)
               if (vec_in%has_rst_fields()) then
                  call vec_in%get_rst(vec_rst, istep)
                  call vec2nek(vxp, vyp, vzp, prp, tp, vec_rst)
               end if
            class default
               call type_error('vec_in','nek_dvector','IN',this_module, this_procedure)
            end select
         end procedure floquet_get_rst

      end submodule
