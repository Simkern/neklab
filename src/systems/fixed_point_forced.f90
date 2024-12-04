      submodule(neklab_systems) fixed_point_forced
         implicit none
      contains
         module procedure nonlinear_map_forced
      ! internal
            character(len=128) :: msg
            real(dp) :: Tend
            complex(dp), dimension(:), allocatable :: dpds
            select type (vec_in)
            type is (nek_ext_dvector_forcing)
               select type (vec_out)
               type is (nek_ext_dvector_forcing)
                  Tend = param(10)
      ! Set the initial condition
                  call ext_vec_f2nek(vx, vy, vz, pr, t, vec_in)
                  call dpds_from_vector(dpds, vec_in)
                  call pipe%set_dpds(dpds)
                  if (.not. pipe%is_steady()) Tend = pipe%get_period()
                  ! Set appropriate tolerances and Nek status
                  call setup_nonlinear_solver(recompute_dt=.true., endtime=Tend,
     $   cfl_limit = 0.4_dp, vtol = atol/10.0, ptol = atol/10.0)
                  write (msg, '(A,*(F9.6,1X))') 'Current forcing estimate, f = ', vec_in%f
                  if (nid == 0) print *, msg
                  call logger%log_message(msg, module=this_module, procedure='nonlinear_map_forced')
      ! Intgrate the nonlinear equations forward
                  time = 0.0_dp
                  do istep = 1, nsteps
                     call pipe%compute_bf_forcing(time) ! --> set neklab_forcing data
                     call nek_advance()
                  end do
                  ! Copy the final solution to vector.
                  call nek2ext_vec_f(vec_out, vx, vy, vz, pr, t)
                  allocate(vec_out%f(vec_in%nf))
                  vec_out%nf = vec_in%nf
                  vec_out%f = vec_in%f
      ! Evaluate residual F(X) - X.
                  call vec_out%sub(vec_in)
               end select
            end select
         end procedure nonlinear_map_forced
   
         module procedure jac_direct_map_forced
      ! internal
            integer :: i
            real(dp) :: atol, df, g_in, g_out
            complex(dp), dimension(:), allocatable :: dpds
            type(nek_ext_dvector_forcing) :: vec
            select type (vec_in)
            type is (nek_ext_dvector_forcing)
               select type (vec_out)
               type is (nek_ext_dvector_forcing)
                  atol = param(22)
      ! Set the baseflow initial condition
                  call abs_ext_vec_f2nek(vx, vy, vz, pr, t, self%X)
                  call dpds_from_vector(dpds, vec_in)
                  call pipe%set_dpds(dpds)
      ! Ensure correct nek status -> set end time
                  call setup_linear_solver(solve_baseflow=.true., transpose=.false.,
     $   recompute_dt = .true., cfl_limit = 0.4_dp, vtol = atol/2.0, ptol = atol/2.0)
      ! Set the perturbation initial condition
                  call ext_vec_f2nek(vxp, vyp, vzp, prp, tp, vec_in)
      ! Intgrate the coupled equations forward
                  time = 0.0_dp
                  do istep = 1, nsteps
                     call pipe%compute_bf_forcing(time) ! --> set neklab_forcing data
                     call nek_advance()
                  end do
      ! Copy the final solution to vector.
                  call nek2ext_vec_f(vec_out, vxp, vyp, vzp, prp, tp)
      ! Evaluate [ exp(tau*J) - I ] @ dx.
                  call vec_out%sub(vec_in)
                  df = 1e-4_dp
      ! Set up nonlinear solver
                  call setup_nonlinear_solver(recompute_dt=.true.,
     $   cfl_limit = 0.4_dp, vtol = atol/10.0, ptol = atol/10.0)
                  do i = 1, vec_in%nf
      ! Evaluate G(du)
                     vec_out%f(i) = pipe%compute_ubar(vxp,vyp,vzp)
      ! Reset baseflow
                     call abs_ext_vec_f2nek(vx, vy, vz, pr, t, self%X)
                     g_in = pipe%compute_ubar(vx,vy,vz)
                     call compute_fdot_forcing(vec, df, i)
                     g_out = pipe%compute_ubar(vx,vy,vz)
      ! Add dF/df*deltaf to output vector
                     call vec_out%axpby(1.0_dp, vec, vec_in%f(i))
      ! Add dG/df*deltaf to the output forcing
                     vec_out%f(i) = vec_out%f(i) + (g_out-g_in)/df*vec_in%f(i)
                  end do
               end select
            end select
         end procedure jac_direct_map_forced
      
         module procedure jac_adjoint_map_forced
      ! internal
            integer :: i
            real(dp) :: atol, df, g_in, g_out
            complex(dp), dimension(:), allocatable :: dpds
            type(nek_ext_dvector_forcing) :: vec
            select type (vec_in)
            type is (nek_ext_dvector_forcing)
               select type (vec_out)
               type is (nek_ext_dvector_forcing)
                  atol = param(22)
      ! Set the baseflow initial condition
                  call abs_ext_vec_f2nek(vx, vy, vz, pr, t, self%X)
                  call dpds_from_vector(dpds, vec_in)
                  call pipe%set_dpds(dpds)
      ! Ensure correct nek status -> set end time
                  call setup_linear_solver(solve_baseflow=.true., transpose=.true.,
     $   recompute_dt = .true., cfl_limit = 0.4_dp, vtol = atol/2.0, ptol = atol/2.0)
      ! Set the perturbation initial condition
                  call ext_vec_f2nek(vxp, vyp, vzp, prp, tp, vec_in)
      ! Intgrate the coupled equations forward
                  time = 0.0_dp
                  do istep = 1, nsteps
                     call pipe%compute_bf_forcing(time) ! --> set neklab_forcing data
                     call nek_advance()
                  end do
      ! Copy the final solution to vector.
                  call nek2ext_vec_f(vec_out, vxp, vyp, vzp, prp, tp)
      ! Evaluate [ exp(tau*J) - I ] @ dx.
                  call vec_out%sub(vec_in)
                  df = 1e-4_dp
      ! Set up nonlinear solver
                  call setup_nonlinear_solver(recompute_dt=.true.,
     $   cfl_limit = 0.4_dp, vtol = atol/10.0, ptol = atol/10.0)
                  do i = 1, vec_in%nf
      ! Evaluate G(du)
                     vec_out%f(i) = pipe%compute_ubar(vxp,vyp,vzp)
      ! Reset baseflow
                     call abs_ext_vec_f2nek(vx, vy, vz, pr, t, self%X)
                     g_in = pipe%compute_ubar(vx,vy,vz)
                     call compute_fdot_forcing(vec, df, i)
                     g_out = pipe%compute_ubar(vx,vy,vz)
      ! Add dF/df*deltaf to output vector
                     call vec_out%axpby(1.0_dp, vec, vec_in%f(i))
      ! Add dG/df*deltaf to the output forcing
                     vec_out%f(i) = vec_out%f(i) + (g_out-g_in)/df*vec_in%f(i)
                  end do
               end select
            end select
         end procedure jac_adjoint_map_forced
      end submodule
