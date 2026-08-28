      !--------------------------------------------------------------------------
      !
      !  PERIODIC ORBITS ON THE 2Dh (AXISYMMETRIC TORUS) MESH
      !
      !  The residual is
      !
      !     F(X) = Phi_T(X) - X
      !
      !  with T = 2*pi/omega imposed by the pulsatile forcing. The problem is
      !  NON-AUTONOMOUS, so there is no phase freedom and no unknown period:
      !  the Poincare section is t = 0 and the state vector is a plain
      !  nek_dvector, exactly as for nek_system_torus_upo. The forcing is NOT
      !  part of the unknown vector -- it is held in t2Dh and driven by the
      !  segregated outer Newton in neklab_analysis_torus_2Dh. There is no
      !  bordered vector and no perturbation forcing: the linearised equations
      !  carry no source, so userq must return zero for jp > 0.
      !
      !  The Jacobian is one linearised integration about the recorded
      !  trajectory. Because the baseflow is time dependent it cannot be held in
      !  a single field: the nonlinear pass records it snapshot by snapshot into
      !  neklab_bf_buffer, and each matvec replays it. This mirrors the 3D helix
      !  path (periodic_orbit_torus.f90), with the buffer standing in for the
      !  pipe object's 2D slice store and nek_advance_2Dh_axisym standing in for
      !  nek_advance.
      !
      !  The nonlinear pass also accumulates the flow-rate Fourier coefficients
      !  and closes them at the horizon (t2Dh%close_mflow). That measurement
      !  belongs to the trajectory, not to its endpoint, which is why it is made
      !  here rather than in the driver: the outer Newton simply reads it back
      !  through t2Dh%measure_mflow.
      !
      !  TIME GRID
      !
      !  dt is CFL-adaptive: over a pulsatile period the flow changes enough
      !  that a fixed step would either violate the CFL condition at peak flow
      !  or waste an order of magnitude of steps at the trough. The horizon is
      !  therefore endtime-driven, and the step count is whatever the CFL
      !  target produces.
      !
      !  Two consequences are handled explicitly:
      !
      !   1. The last step would otherwise be truncated to land on T, giving a
      !      sliver whose length jumps discontinuously with X and whose BDF3
      !      local error dominates. The landing logic below splits the remaining
      !      interval over the last N_LAND steps instead. See land_horizon.
      !
      !   2. A horizon resolved by a handful of steps is not a discretisation of
      !      the dynamics. bf_close_record fails below bf_get_nsteps_min.
      !
      !--------------------------------------------------------------------------

      submodule(neklab_systems) periodic_orbit_torus_2Dh
         implicit none

         integer, parameter :: n_land = 2
      !! Number of equal steps used to land exactly on T.
         real(dp), parameter :: land_factor = 2.5_dp
      !! Engage the landing when the remaining interval falls below
      !! land_factor*dt. setdt caps the growth of dt at 1.2 per step
      !! (subs1.f:293-297), so n_land further steps can consume at most
      !! 1.2 + 1.44 = 2.64 times the current dt; 2.5 leaves no room for a
      !! surprise overshoot while engaging as late as possible.
         real(dp), parameter :: cfl_upo = 0.4_dp

         logical, parameter :: if_upo_adjoint = .false.
      !! Master switch for the adjoint replay. The body is written out below so
      !! that turning it on is a one-line change here, but two other things must
      !! be done first:
      !!   * neklab_2Dh_axisym.f90 refuses ifadj outright (line ~435), and the
      !!     adjoint of the torus curvature coupling has to be written;
      !!   * bf_set_window cannot reach across a chunk boundary, so a backward
      !!     replay of a trajectory spanning more than one chunk needs the
      !!     two-snapshot side cache described in neklab_bf_buffer.

      contains

      !====================================================================
      !     NONLINEAR RESIDUAL
      !====================================================================

         module procedure nonlinear_map_torus_upo_2Dh
      ! internal
         character(len=*), parameter :: this_procedure = 'nonlinear_map_torus_upo_2Dh'
         character(len=256) :: msg
         real(dp) :: period, trem, dtland
         logical :: landing
         integer :: iland
         select type (vec_in)
         type is (nek_dvector)
            select type (vec_out)
            type is (nek_dvector)
               period = t2Dh%get_period()
               if (period <= 0.0_dp) then
                  call nek_stop_error('Period is not set. The driver must call t2Dh%init_flow '//
     &               'with a non-zero Womersley number first.', this_module, this_procedure)
               end if
      ! Set the initial condition
               call vec2nek(vx, vy, vz, pr, t, vec_in)
      ! Set appropriate tolerances and Nek status
               call setup_nonlinear_solver(variable_dt = .true.,
     &                                     endtime     = period,
     &                                     solve_temperature = .true.,
     &                                     cfl_limit   = cfl_upo,
     &                                     vtol        = atol*0.1,
     &                                     ptol        = atol*0.1)
      ! Start the recording and the flow-rate accumulator. The prefix belongs to
      ! the driver ('n' during the Newton solve, 'b' for the converged orbit),
      ! so bf_reset is called without one. The accumulator is seeded with the
      ! flow rate of the initial condition, so the trapezoidal rule has a left
      ! endpoint for the first step.
               call bf_reset()
               call t2Dh%reset_mflow(t2Dh%ubar())
      ! Integrate the nonlinear equations forward over exactly one period
               time = 0.0_dp
               istep = 0
               landing = .false.
               iland = 0
               do 
                  istep = istep + 1
                  if (.not. landing .and. istep > 2) then
                     trem = period - time
                     if (trem <= land_factor*dt) then
                        landing = .true.
                        dtland = trem/real(n_land, dp)
                        param(12) = -dtland
                        write (msg, '(A,I0,A,E16.8,A,E16.8,A)') 'Landing over ', n_land,
     &                     ' steps: dt= ', dtland, ' (was ', dt, ')'
                        call nek_log_debug(msg, this_module, this_procedure)
                     end if
                  end if
                  call bf_begin_step()
                  call neklab_timer_start(t_nl_step)
                  call nek_advance()
                  call neklab_timer_stop(t_nl_step)
      ! The landing steps are sized by the period constraint, not by the CFL
      ! condition, so they must not pollute the dt statistics the driver uses
      ! to report the resolution of the orbit.
                  call bf_end_step(count_stats = .not. landing)
                  call t2Dh%accumulate_mflow(t2Dh%ubar(), time, dt)
      ! The landing consumes exactly n_land steps by construction. Do not let
      ! setdt's TIME+DT >= FINTIM test decide: n_land equal steps sum to T only
      ! in exact arithmetic, and a shortfall of one ulp costs a whole extra
      ! step because subs1.f:306 undoes the clipping of dt to the remainder.
                  if (landing) then
                     iland = iland + 1
                     if (iland == n_land) exit
                  else if (lastep /= 0) then
                     exit   ! safety net: landing never engaged
                  end if
               end do
               lastep = 1
               write (msg, '(A,I0,A,E16.8)') 'Landed in ', iland,
     &            ' steps, period - time= ', period - time
               call nek_log_debug(msg, this_module, this_procedure)
      ! Close the recording. This is where sum(dt) == T and the minimum step
      ! count are enforced: both failures would otherwise show up much later as
      ! an inconsistency between F and dF.
               call bf_close_record(period=period)
      ! Close the flow-rate accumulation. This normalises both quadrature rules,
      ! forms the amplitudes and phases and measures the quadrature error, which
      ! is what floors the outer Newton's tolerance.
               call neklab_timer_start(t_close_mflow)
               call t2Dh%close_mflow(period=period)
               call neklab_timer_stop(t_close_mflow)
      ! Copy the final solution to vector.
               call nek2vec(vec_out, vx, vy, vz, pr, t)
      ! Evaluate residual F(X) - X.
               call vec_out%sub(vec_in)
      ! One CSV row for this F evaluation. Every timer read there is idle at
      ! this point; get_data would STOP a running one.
               call neklab_timer_dump('F', eval=self%get_eval_counter(),
     &                                nsteps=bf_get_nsteps())
            class default
               call nek_stop_error("The intent [OUT] argument 'vec_out' must be of type 'nek_dvector'",
     & this_module, this_procedure)
            end select
         class default
            call nek_stop_error("The intent [IN] argument 'vec_in' must be of type 'nek_dvector'",
     & this_module, this_procedure)
         end select
         end procedure nonlinear_map_torus_upo_2Dh

      !====================================================================
      !     JACOBIAN -- DIRECT
      !====================================================================

         module procedure jac_direct_map_torus_upo_2Dh
      ! internal
         character(len=*), parameter :: this_procedure = 'jac_direct_map_torus_upo_2Dh'
         real(dp) :: atol
         integer :: nsteps_bf
         select type (vec_in)
         type is (nek_dvector)
            select type (vec_out)
            type is (nek_dvector)
               atol = param(22)
      ! Install the Newton iterate as the baseflow. The trajectory itself comes
      ! from the buffer step by step; this is here so that setup_linear_solver
      ! sees a representative field when it sizes things from the CFL.
               call abs_vec2nek(vx, vy, vz, pr, t, self%X)
      ! Ensure correct nek status. variable_dt because the replay overrides dt
      ! per step from the buffer; endtime only sets fintim, the loop below is
      ! bounded by the recorded step count.
               call setup_linear_solver(solve_baseflow = .false.,
     &                                  solve_temperature = .true.,
     &                                  variable_dt    = .true.,
     &                                  endtime        = t2Dh%get_period(),
     &                                  cfl_limit      = cfl_upo,
     &                                  vtol           = atol*0.5,
     &                                  ptol           = atol*0.5)
      ! Set the initial condition for Nek5000's linearized solver.
               call vec2nek(vxp, vyp, vzp, prp, tp, vec_in)
      ! Replay the recorded trajectory. bf_set installs snapshot k and forces
      ! the timestep the nonlinear pass used at step k, so the linearised
      ! operator is the exact Jacobian of the discrete map that produced F --
      ! not of some nearby one.
               nsteps_bf = bf_get_nsteps()
               call bf_replay_start()
               time = 0.0_dp
               do istep = 1, nsteps_bf
                  call bf_set(istep)
                  call nek_advance_2Dh_axisym(0.0_dp)
               end do
               call bf_replay_end()
      ! Extract the final solution to vector.
               call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)
      ! Evaluate [ M_T - I ] @ dx.
               call vec_out%sub(vec_in)
               param(22) = atol
               call neklab_timer_dump('J', matvec=self%get_counter(.false.),
     &                                nsteps=nsteps_bf)
            class default
               call nek_stop_error("The intent [OUT] argument 'vec_out' must be of type 'nek_dvector'",
     & this_module, this_procedure)
            end select
         class default
            call nek_stop_error("The intent [IN] argument 'vec_in' must be of type 'nek_dvector'",
     & this_module, this_procedure)
         end select
         end procedure jac_direct_map_torus_upo_2Dh

      !====================================================================
      !     JACOBIAN -- ADJOINT
      !====================================================================
      !
      !  The adjoint of M_T = A_n ... A_1 is A_1^T ... A_n^T, so the baseflow
      !  must be presented in reverse order. bf_set cannot do that: it builds
      !  the lag levels by shifting the current field down one level, which
      !  going backwards would put buf(k+1) where buf(k-1) belongs. Hence
      !  bf_set_window, which loads slots k, k-1, k-2 explicitly.
      !
      !  The body is written out so the path is short when this goes live, but
      !  it is unreachable until if_upo_adjoint is flipped -- see the note at
      !  the top of this submodule for what else has to happen first.

         module procedure jac_adjoint_map_torus_upo_2Dh
      ! internal
         character(len=*), parameter :: this_procedure = 'jac_adjoint_map_torus_upo_2Dh'
         real(dp) :: atol
         integer :: k, nsteps_bf
         if (.not. if_upo_adjoint) then
            call nek_stop_error('rmatvec is not available for the 2Dh periodic-orbit system: '//
     &         'nek_advance_2Dh_axisym refuses adjoint mode. Use a non-transpose solver (gmres).',
     &         this_module, this_procedure)
         end if
         select type (vec_in)
         type is (nek_dvector)
            select type (vec_out)
            type is (nek_dvector)
               atol = param(22)
               call abs_vec2nek(vx, vy, vz, pr, t, self%X)
               call setup_linear_solver(transpose      = .true.,
     &                                  solve_baseflow = .false.,
     &                                  solve_temperature = .true.,
     &                                  variable_dt    = .true.,
     &                                  endtime        = t2Dh%get_period(),
     &                                  cfl_limit      = cfl_upo,
     &                                  vtol           = atol*0.5,
     &                                  ptol           = atol*0.5)
               call vec2nek(vxp, vyp, vzp, prp, tp, vec_in)
      ! Backward replay: step istep of the adjoint integration linearises about
      ! step k = nsteps - istep + 1 of the forward trajectory.
               nsteps_bf = bf_get_nsteps()
               call bf_replay_start(reverse=.true.)
               time = 0.0_dp
               do istep = 1, nsteps_bf
                  k = nsteps_bf - istep + 1
                  call bf_set_window(k)
                  call nek_advance_2Dh_axisym(0.0_dp)
               end do
               call bf_replay_end()
               call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)
      ! Evaluate [ M_T^T - I ] @ dx.
               call vec_out%sub(vec_in)
               param(22) = atol
               call neklab_timer_dump('JT', rmatvec=self%get_counter(.true.),
     &                                nsteps=nsteps_bf)
            class default
               call nek_stop_error("The intent [OUT] argument 'vec_out' must be of type 'nek_dvector'",
     & this_module, this_procedure)
            end select
         class default
            call nek_stop_error("The intent [IN] argument 'vec_in' must be of type 'nek_dvector'",
     & this_module, this_procedure)
         end select
         end procedure jac_adjoint_map_torus_upo_2Dh

      end submodule periodic_orbit_torus_2Dh