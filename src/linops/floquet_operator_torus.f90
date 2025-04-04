      submodule(neklab_linops) floquet_operator_torus
         implicit none
      contains
         module procedure floquet_init
      ! internal
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
            call nek_log_message(msg, module=this_module, procedure="floquet_init")
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
            call nek_log_message(msg, module=this_module, procedure="floquet_init")
         else
            self%baseflow_computed = .false.
            msg = "No 2d bseflow files found. Solve for baseflow at first iteration."
            call nek_log_message(msg, module=this_module, procedure="floquet_init")
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
     &                                     endtime        = self%tau, 
     &                                     cfl_limit      = 0.4_dp,
     &                                     variable_dt    = .true.) ! -> solve for baseflow and save to f2dtorus***.fld
               else
                  call pipe%set_save_base(.false.)
                  call setup_linear_solver(solve_baseflow = .false., 
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
                     call nek_advance()
                     call pipe%save_2d_fields(vx,vy,vz) ! outposts automatically at lastep == 1
                  end do
         ! Record the mass flow rate and number of timesteps per period for subsequent linear runs
                  call pipe%set_nsteps(istep)
                  self%baseflow_computed = .true. ! we only need to do this once
               else
                  do istep = 1, pipe%get_nsteps()
                     call nek_advance()
                     call pipe%set_baseflow(vx, vy, vz, istep) ! sets the baseflow field and the appropriate timestep
                  end do
               end if
      ! Extract the final solution to vector.
               call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)

               self%baseflow_computed = .true. ! we only need to do this once
               class default
               call stop_error("The intent [OUT] argument 'vec_out' must be of type 'nek_dvector'",
     & this_module, 'floquet_matvec')
            end select
         class default
            call stop_error("The intent [IN] argument 'vec_in' must be of type 'nek_dvector'",
     & this_module, 'floquet_matvec')
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
     &                                     solve_baseflow = .true.,
     &                                     endtime        = self%tau, 
     &                                     cfl_limit      = 0.4_dp,
     &                                     variable_dt    = .true.) ! -> solve for baseflow and save to f2dtorus***.fld
                  self%baseflow_computed = .true. ! we only need to do this once
               else
                  call pipe%set_save_base(.false.)
                  call setup_linear_solver(transpose      = .true.,
     &                                     solve_baseflow = .false., 
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
                     call nek_advance()
                     call pipe%save_2d_fields(vx,vy,vz) ! outposts automatically at lastep == 1
                  end do
         ! Record the mass flow rate and number of timesteps per period for subsequent linear runs
                  call pipe%set_nsteps(istep)
                  self%baseflow_computed = .true. ! we only need to do this once
               else
                  do istep = 1, pipe%get_nsteps()
                     call pipe%set_baseflow(vx, vy, vz, istep) ! sets the baseflow field and the appropriate timestep
                     call nek_advance()
                  end do
               end if
      ! Extract the final solution to vector.
               call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)
            class default
               call stop_error("The intent [OUT] argument 'vec_out' must be of type 'nek_dvector'",
     & this_module, 'floquet_rmatvec')
            end select
         class default
            call stop_error("The intent [IN] argument 'vec_in' must be of type 'nek_dvector'",
     & this_module, 'floquet_rmatvec')
         end select
         end procedure floquet_rmatvec

      end submodule
