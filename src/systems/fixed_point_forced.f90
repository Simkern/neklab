      submodule(neklab_systems) fixed_point_forced
         implicit none
      contains
         module procedure nonlinear_map_forced
      ! internal
            character(len=128) :: msg
            real(dp) :: Tend
            select type (vec_in)
            type is (nek_ext_dvector_forcing)
               select type (vec_out)
               type is (nek_ext_dvector_forcing)
                  Tend = param(10)
      ! Set the initial condition
                  call ext_vec_f2nek(vx, vy, vz, pr, t, vec_in)
                  call pipe%set_dpds(vec_in%f)
                  if (.not. pipe%is_steady()) Tend = pipe%get_period()
                  ! Set appropriate tolerances and Nek status
                  call setup_nonlinear_solver(recompute_dt=.true., endtime=Tend,
     $   cfl_limit = 0.5_dp, vtol = atol/10.0, ptol = atol/10.0)
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
                  vec_out%f = vec_in%f
      ! Evaluate residual F(X) - X.
                  call vec_out%sub(vec_in)
      ! Prepare computation of dFdf
                  pipe%to_compute_df = .true.
               end select
            end select
         end procedure nonlinear_map_forced
   
         module procedure jac_direct_map_forced
      ! internal
            integer :: i
            real(dp) :: atol, df, g_in, g_out
            type(nek_ext_dvector_forcing) :: F_in, F_out
            character(len=128) :: msg
            select type (vec_in)
            type is (nek_ext_dvector_forcing)
               select type (vec_out)
               type is (nek_ext_dvector_forcing)
                  atol = param(22)
      ! Set the baseflow initial condition
                  call abs_ext_vec_f2nek(vx, vy, vz, pr, t, self%X)
      ! Ensure correct nek status
                  call setup_linear_solver(solve_baseflow=.true., transpose=.false.,
     $   recompute_dt = .true., cfl_limit = 0.5_dp, vtol = atol/2.0, ptol = atol/2.0)
      ! Set the perturbation initial condition
                  call ext_vec_f2nek(vxp, vyp, vzp, prp, tp, vec_in)
                  call outpost(vxp,vyp,vzp,prp,tp,'lin')
      ! Intgrate the coupled equations forward
                  time = 0.0_dp
                  write(msg,'(A,E15.8)') 'dpds= ', real(pipe%dpds)
                  call logger%log_message(msg, module=this_module, procedure='jac_eval')
                  do istep = 1, nsteps
                     call pipe%compute_bf_forcing(time) ! --> set neklab_forcing data
                     call nek_advance()
                  end do
      ! Copy the final solution to vector.
                  call nek2ext_vec_f(vec_out, vxp, vyp, vzp, prp, tp)
                  call outpost(vxp,vyp,vzp,prp,tp,'lin')
      ! Evaluate [ exp(tau*J) - I ] @ dx.
                  call vec_out%sub(vec_in)
      ! Evaluate approximate derivative w.r.t forcing amplitude
                  if (pipe%to_compute_df) then
                     write(msg,'(A)') 'Computing dFdf and dGdf'
                     call logger%log_message(msg, module=this_module, procedure='jac_eval')
      ! Define perturbation amplitude
                     df = 1e-4_dp
                     write(msg,'(2X,A,*(E12.5))') 'df   = ', df
                     call logger%log_message(msg, module=this_module, procedure='jac_eval') 
      ! Save F(X) and G(X)
                     call nek2ext_vec_f(F_in, vx, vy, vz, pr, t)
                     g_in = pipe%compute_ubar(vx, vy, vz)
                     write(msg,'(2X,A,*(E12.5))') 'g_in = ', g_in
                     call logger%log_message(msg, module=this_module, procedure='jac_eval')
      ! Ensure correct nek status   
                     call setup_nonlinear_solver(recompute_dt=.true.,
     $   cfl_limit = 0.5_dp, vtol = atol/10.0, ptol = atol/10.0)
                     do i = 1, nf
                        write(msg,'(2X,A,I0,A)') 'Component ', i, ':'
                        call logger%log_message(msg, module=this_module, procedure='jac_eval')
      ! Set the baseflow initial condition
                        call abs_ext_vec_f2nek(vx, vy, vz, pr, t, self%X)
      ! Set perturbation
                        call pipe%add_dpds(df, i)
                        write(msg,'(4X,I2,A,*(E12.5))') i, ': dpds  = ', pipe%dpds
                        call logger%log_message(msg, module=this_module, procedure='jac_eval')
      ! Evaluate F(X+df) and G(X+df)
                        do istep = 1, nsteps
                           call pipe%compute_bf_forcing(time) ! --> set neklab_forcing data
                           call nek_advance()
                        end do
      ! Extract result and evaluate (F(X+df) - F(X))/df and (G(X+df) - G(X))/df
                        call nek2ext_vec_f(F_out, vx, vy, vz, pr, t)
                        call F_out%sub(F_in)
                        call F_out%scal(1.0_dp/df)
                        F_out%f = 0.0_dp
                        call pipe%set_dFdf(F_out, i)
                        write(msg,'(4X,I2,A,*(E12.5))') i, ':|dFdf| = ', F_out%norm()
                        call logger%log_message(msg, module=this_module, procedure='jac_eval') 
                        g_out = pipe%compute_ubar(vx,vy,vz)
                        call pipe%set_dGdf((g_out - g_in)/df, i)
                        write(msg,'(4X,I2,A,*(E12.5))') i, ': g_out = ', g_out
                        call logger%log_message(msg, module=this_module, procedure='jac_eval') 
                        write(msg,'(4X,I2,A,*(E12.5))') i, ': dGdf  = ', pipe%get_dGdf(i)
                        call logger%log_message(msg, module=this_module, procedure='jac_eval') 
      ! Reset perturbation
                        call pipe%add_dpds(-df, i)
                     end do
                     pipe%to_compute_df = .false.
                  end if
      ! Set the forcing contribution for F
                  do i = 1, nf
                     call pipe%get_dFdf(F_out,i)
                     call vec_out%axpby(1.0_dp, F_out, vec_in%f(i))
                  end do
      ! Evaluate G(du_in)
                  do i = 1, nf
                     vec_out%f(i) = pipe%compute_ubar(vec_in%vx, vec_in%vy, vec_in%vz)
                     write(msg,'(A,E15.8)') 'G(du) = ', vec_out%f(i)
                     call logger%log_message(msg, module=this_module, procedure='jac_eval')
                  end do
      ! Set the forcing contribution for G
                  vec_out%f = vec_out%f + pipe%dGdf
                  write(msg,'(A,E15.8)') 'dGdf  = ', pipe%dGdf
                  call logger%log_message(msg, module=this_module, procedure='jac_eval')
                  param(22) = atol
                  param(21) = atol
               end select
            end select
         end procedure jac_direct_map_forced
      
         module procedure jac_adjoint_map_forced
      ! internal
            integer :: i
            real(dp) :: atol, df, g_in, g_out
            type(nek_ext_dvector_forcing) :: F_in, F_out
            character(len=128) :: msg
            select type (vec_in)
            type is (nek_ext_dvector_forcing)
               select type (vec_out)
               type is (nek_ext_dvector_forcing)
                  atol = param(22)
      ! Set the baseflow initial condition
                  call abs_ext_vec_f2nek(vx, vy, vz, pr, t, self%X)
      ! Ensure correct nek status
                  call setup_linear_solver(solve_baseflow=.true., transpose=.true.,
     $   recompute_dt = .true., cfl_limit = 0.5_dp, vtol = atol/2.0, ptol = atol/2.0)
      ! Set the perturbation initial condition
                  call ext_vec_f2nek(vxp, vyp, vzp, prp, tp, vec_in)
                  call outpost(vxp,vyp,vzp,prp,tp,'lin')
      ! Intgrate the coupled equations forward
                  time = 0.0_dp
                  write(msg,'(A,E15.8)') 'dpds= ', real(pipe%dpds)
                  call logger%log_message(msg, module=this_module, procedure='jac_eval')
                  do istep = 1, nsteps
                     call pipe%compute_bf_forcing(time) ! --> set neklab_forcing data
                     call nek_advance()
                  end do
      ! Copy the final solution to vector.
                  call nek2ext_vec_f(vec_out, vxp, vyp, vzp, prp, tp)
                  call outpost(vxp,vyp,vzp,prp,tp,'lin')
      ! Evaluate [ exp(tau*J) - I ] @ dx.
                  call vec_out%sub(vec_in)
      ! Evaluate approximate derivative w.r.t forcing amplitude
                  if (pipe%to_compute_df) then
                     write(msg,'(A)') 'Computing dFdf and dGdf'
                     call logger%log_message(msg, module=this_module, procedure='jac_eval')
      ! Define perturbation amplitude
                     df = 1e-4_dp
                     write(msg,'(2X,A,*(E12.5))') 'df   = ', df
                     call logger%log_message(msg, module=this_module, procedure='jac_eval') 
      ! Save F(X) and G(X)
                     call nek2ext_vec_f(F_in, vx, vy, vz, pr, t)
                     g_in = pipe%compute_ubar(vx, vy, vz)
                     write(msg,'(2X,A,*(E12.5))') 'g_in = ', g_in
                     call logger%log_message(msg, module=this_module, procedure='jac_eval')
      ! Ensure correct nek status   
                     call setup_nonlinear_solver(recompute_dt=.true.,
     $   cfl_limit = 0.5_dp, vtol = atol/10.0, ptol = atol/10.0)
                     do i = 1, nf
                        write(msg,'(2X,A,I0,A)') 'Component ', i, ':'
                        call logger%log_message(msg, module=this_module, procedure='jac_eval')
      ! Set the baseflow initial condition
                        call abs_ext_vec_f2nek(vx, vy, vz, pr, t, self%X)
      ! Set perturbation
                        call pipe%add_dpds(df, i)
                        write(msg,'(4X,I2,A,*(E12.5))') i, ': dpds  = ', pipe%dpds
                        call logger%log_message(msg, module=this_module, procedure='jac_eval')
      ! Evaluate F(X+df) and G(X+df)
                        do istep = 1, nsteps
                           call pipe%compute_bf_forcing(time) ! --> set neklab_forcing data
                           call nek_advance()
                        end do
      ! Extract result and evaluate (F(X+df) - F(X))/df and (G(X+df) - G(X))/df
                        call nek2ext_vec_f(F_out, vx, vy, vz, pr, t)
                        call F_out%sub(F_in)
                        call F_out%scal(1.0_dp/df)
                        F_out%f = 0.0_dp
                        call pipe%set_dFdf(F_out, i)
                        write(msg,'(4X,I2,A,*(E12.5))') i, ':|dFdf| = ', F_out%norm()
                        call logger%log_message(msg, module=this_module, procedure='jac_eval') 
                        g_out = pipe%compute_ubar(vx,vy,vz)
                        call pipe%set_dGdf((g_out - g_in)/df, i)
                        write(msg,'(4X,I2,A,*(E12.5))') i, ': g_out = ', g_out
                        call logger%log_message(msg, module=this_module, procedure='jac_eval') 
                        write(msg,'(4X,I2,A,*(E12.5))') i, ': dGdf  = ', pipe%get_dGdf(i)
                        call logger%log_message(msg, module=this_module, procedure='jac_eval') 
      ! Reset perturbation
                        call pipe%add_dpds(-df, i)
                     end do
                     pipe%to_compute_df = .false.
                  end if
      ! Set the forcing contribution for F
                  do i = 1, nf
                     call pipe%get_dFdf(F_out,i)
                     call vec_out%axpby(1.0_dp, F_out, vec_in%f(i))
                  end do
      ! Evaluate G(du_in)
                  do i = 1, nf
                     vec_out%f(i) = pipe%compute_ubar(vec_in%vx, vec_in%vy, vec_in%vz)
                     write(msg,'(A,E15.8)') 'G(du) = ', vec_out%f(i)
                     call logger%log_message(msg, module=this_module, procedure='jac_eval')
                  end do
      ! Set the forcing contribution for G
                  vec_out%f = vec_out%f + pipe%dGdf
                  write(msg,'(A,E15.8)') 'dGdf  = ', pipe%dGdf
                  call logger%log_message(msg, module=this_module, procedure='jac_eval')
                  param(22) = atol
                  param(21) = atol
               end select
            end select
         end procedure jac_adjoint_map_forced
      end submodule
