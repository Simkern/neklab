      submodule(neklab_systems) periodic_orbit_torus_upo
         implicit none
      contains
         module procedure nonlinear_map_torus_upo
         character(len=*), parameter :: this_procedure = 'nonlinear_map_torus_upo'
         real(dp), parameter :: eps = 1e-3_dp
         real(dp) :: pd
         select type (vec_in)
         type is (nek_dvector)
            select type (vec_out)
            type is (nek_dvector)

               pd = pipe%get_period()

      ! Set the initial condition for the nonlinear solver.
               call vec2nek(vx, vy, vz, pr, t, vec_in)

      ! Set appropriate tolerances and Nek status
               call setup_nonlinear_solver(variable_dt = .true., 
     &                                     endtime     = pd, 
     &                                     cfl_limit   = 0.4_dp,
     &                                     vtol        = atol*0.1, 
     &                                     ptol        = atol*0.1)

      ! Intgrate the nonlinear equations forward
               time = 0.0_dp

               ! reset output counter to overwrite output files, compute ubar_lag
               call pipe%set_2d_mode('newton')

               istep = 0
               do while (lastep == 0)
                  istep = istep + 1
                  call pipe%compute_bf_forcing(time) ! --> set neklab_forcing data
                  if (abs(time + 10*dt - pd) <= dt) dt = (pd-time)/(10.0_dp - eps)
                  call pipe%save_2d_fields(vx,vy,vz)
                  call nek_advance()
                  call pipe%compute_mflow_fft(var_dt = .true.) ! integrate Fourier coefficients
               end do
               
               call pipe%outpost_2d_fields()
      ! Record the mass flow rate and number of timesteps per period for subsequent linear runs
               call pipe%extract_mflow_fft()
               call pipe%set_nsteps(istep)
      
      ! Copy the final solution to vector.
               call nek2vec(vec_out, vx, vy, vz, pr, t)
      
      ! Evaluate residual F(X) - X.
               call vec_out%sub(vec_in)
      
            class default
               call type_error('vec_out','nek_dvector','OUT',this_module, this_procedure)
            end select
         class default
            call type_error('vec_in','nek_dvector','IN',this_module, this_procedure)
         end select
         end procedure nonlinear_map_torus_upo
      
         module procedure jac_direct_map_torus_upo
         character(len=*), parameter :: this_procedure = 'jac_direct_map_torus_upo'
         integer :: nrst
         real(dp) :: atol
         select type (vec_in)
         type is (nek_dvector)
            select type (vec_out)
            type is (nek_dvector)

               nrst = abs(param(27)) - 1
               atol = param(22)

      ! Set the baseflow initial condition
               call abs_vec2nek(vx, vy, vz, pr, t, self%X)

      ! Ensure correct nek status
               call setup_linear_solver(solve_baseflow = .false., 
     &                                  variable_dt    = .false.,
     &                                  vtol           = atol*0.5, 
     &                                  ptol           = atol*0.5)

      ! Set the initial condition for Nek5000's linearized solver.
               call vec2nek(vxp, vyp, vzp, prp, tp, vec_in)

      ! Integrate the equations forward in time.
               time = 0.0_dp
               call pipe%set_2d_mode('newton')         ! reset output counter to load baseflow files in order
               do istep = 1, pipe%get_nsteps()

                  call pipe%set_baseflow(vx, vy, vz, istep) ! sets the baseflow field and the appropriate timestep
                  call nek_advance()

                  ! Set restart fields if present.
                  if (istep <= nrst) call get_rst_dnek(vec_in, istep)

               end do

      ! Extract the final solution to vector.
               call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)

      ! Compute restart fields.
               call self%compute_rst(vec_out, nrst)

      ! Evaluate [ exp(tau*J) - I ] @ dx.
               call vec_out%sub(vec_in)

               param(22) = atol
            
            class default
               call type_error('vec_out','nek_dvector','OUT',this_module, this_procedure)
            end select
         class default
            call type_error('vec_in','nek_dvector','IN',this_module, this_procedure)
         end select
         end procedure jac_direct_map_torus_upo
      
         module procedure jac_adjoint_map_torus_upo
         character(len=*), parameter :: this_procedure = 'jac_adjoint_map_torus_upo'
         integer :: nrst
         real(dp) :: atol
         select type (vec_in)
         type is (nek_dvector)
            select type (vec_out)
            type is (nek_dvector)

               nrst = abs(param(27)) - 1
               atol = param(22)

      ! Set the baseflow initial condition
               call abs_vec2nek(vx, vy, vz, pr, t, self%X)

      ! Ensure correct nek status
               call setup_linear_solver(transpose      = .true., 
     &                                  solve_baseflow = .false.,
     &                                  variable_dt    = .false.,
     &                                  vtol           = atol*0.5, 
     &                                  ptol           = atol*0.5)

      ! Set the initial condition for Nek5000's linearized solver.
               call vec2nek(vxp, vyp, vzp, prp, tp, vec_in)

      ! Integrate the equations forward in time.
               time = 0.0_dp
               ! reset output counter to load baseflow files in order
               call pipe%set_2d_mode('newton')
               do istep = 1, pipe%get_nsteps()

                  ! sets the baseflow field and the appropriate timestep
                  call pipe%set_baseflow(vx, vy, vz, istep)
                  call nek_advance()

                  ! Set restart fields if present.
                  if (istep <= nrst) call get_rst_dnek(vec_in, istep)

               end do

      ! Extract the final solution to vector.
               call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)

      ! Compute restart fields.
               call self%compute_rst(vec_out, nrst)

      ! Evaluate [ exp(tau*J) - I ] @ dx.
               call vec_out%sub(vec_in)

               param(22) = atol

            class default
               call type_error('vec_out','nek_dvector','OUT',this_module, this_procedure)
            end select
         class default
            call type_error('vec_in','nek_dvector','IN',this_module, this_procedure)
         end select
         end procedure jac_adjoint_map_torus_upo

         module procedure compute_rst_torus_upo
            character(len=*), parameter :: this_procedure = 'compute_rst_torus_upo'
            type(nek_dvector), allocatable :: vec_rst
            character(len=128) :: msg
            real(dp) :: rtmp
            select type(vec_out)
            type is (nek_dvector)
               write(msg,'(A,I0,A)') 'Run ', nrst, ' extra step(s) to fill up restart arrays.'
               call nek_log_debug(msg, this_module, this_procedure)
               fintim = fintim + nrst*dt
               allocate(vec_rst)
               ! reset output counter to overwrite output files, compute ubar_lag
               call pipe%set_2d_mode('newton')
               do istep = nsteps + 1, nsteps + nrst
                  call pipe%set_baseflow(vx, vy, vz, istep) ! sets the baseflow field and the appropriate timestep
                  call nek_advance()
                  call nek2vec(vec_rst, vxp, vyp, vzp, prp, tp)
                  call vec_out%save_rst(vec_rst, istep - nsteps)
               end do
            class default
               call type_error('vec_out','nek_dvector','INOUT',this_module, this_procedure)
            end select
         end procedure compute_rst_torus_upo
      end submodule
