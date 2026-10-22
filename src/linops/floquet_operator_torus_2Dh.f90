      submodule(neklab_linops) floquet_operator_2Dh
         use t2Dh_bf_buffer, only: bf_init, bf_is_initialised, bf_is_recording,
     &                             bf_get_nsteps, bf_get_nchunk, bf_get_time,
     &                             bf_get_dt, bf_set, bf_replay_start, bf_replay_end,
     &                             bf_set_prefix, bf_read_t2Dh, bf_fname
 
         use neklab_t2Dh, only: t2Dh, lfc_t2Dh => lfc
         implicit none

         real(dp), parameter :: cfl_floquet = 0.4_dp
      !! Only sizes the solver setup: the replay forces the recorded dt on
      !! every step. Same value as the nonlinear pass (cfl_upo).

      contains

      !====================================================================
      !     INIT
      !====================================================================

         module procedure init_floquet_2Dh_axisym
            character(len=*), parameter :: this_procedure = 'init_floquet_2Dh_axisym'
            character(len=256) :: msg
            real(dp), dimension(lfc_t2Dh) :: dpds_s
            real(dp) :: omega_s, period_s
            integer :: nf_s, ierr
            character(len=132) :: base_

            call init_2Dh_axisym(.true.) ! no-op if already done

               ! ---- attach the orbit written by the Newton run. The buffer's base is
      !      fixed at bf_init: in a fresh job it is set here, in the Newton job
      !      it is the driver's, and the two must agree. Blank %base: case name.
            base_ = self%base
            if (len_trim(base_) == 0) base_ = '2d'//session
            if (.not. bf_is_initialised()) then
      ! write_chunks = .true.: load_chunk refuses to read when spilling is off.
      ! Nothing is recorded here, so nothing is written.
               call bf_init(base=trim(base_), write_chunks=.true.)
            else if (trim(bf_fname(1, self%prefix)) /= self%prefix//trim(base_)//'001.fld') then
               call nek_stop_error('The buffer was initialised with another file base than '//
     &         self%prefix//trim(base_)//'. Set %base to the driver''s buffer_base.',
     &         this_module, this_procedure)
            end if
            call bf_set_prefix(self%prefix)
         ! Reads nsteps, nchunk and T from <prefix><base>.t2Dh and drops the
         ! resident chunk, so the replay reads the chunk files.
            call bf_read_t2Dh(dpds_s, nf_s, omega_s, period_s, ierr)
            if (ierr /= 0 .or. bf_get_nsteps() == 0) then
               call nek_stop_error('Could not read the orbit '//self%prefix//trim(base_)//
     &         ': the .t2Dh sidecar or the chunk files are missing. Run the Newton solve '//
     &         '(unsteady_flowrate_newton) first.', this_module, this_procedure)
            end if
    
         ! ---- the period is the recorded one
            self%tau = bf_get_time()
            if (self%tau <= 0.0_dp) then
               call nek_stop_error('The recorded orbit has zero length.', this_module, this_procedure)
            end if
            if (t2Dh%is_initialised()) then
               if (t2Dh%is_unsteady()) then
                  if (abs(self%tau - t2Dh%get_period()) > 1.0e-08_dp*self%tau) then
                     write (msg, '(A,E16.8,A,E16.8)') 'Recorded horizon ', self%tau,
     &               ' differs from the t2Dh period ', t2Dh%get_period()
                     call nek_log_warning(msg, this_module, this_procedure)
                  end if
               end if
            end if
    
            call nek_log_message('', this_module, this_procedure)
            call nek_log_message('Floquet operator (2Dh axisym):', this_module, this_procedure)
            write (msg, '(3X,A,E16.8)') 'alpha  = ', self%alpha
            call nek_log_message(msg, this_module, this_procedure)
            write (msg, '(3X,A,E16.8)') 'T      = ', self%tau
            call nek_log_message(msg, this_module, this_procedure)
            write (msg, '(3X,A,I0,A,I0,A)') 'nsteps = ', bf_get_nsteps(), ' (', bf_get_nchunk(), ' chunk(s))'
            call nek_log_message(msg, this_module, this_procedure)
            write (msg, '(3X,A,A)') 'files  = ', trim(bf_fname(1))
            call nek_log_message(msg, this_module, this_procedure)
            call nek_log_message('', this_module, this_procedure)
    
            self%is_initialized = .true.
         end procedure init_floquet_2Dh_axisym

      !====================================================================
      !     DIRECT
      !====================================================================

         module procedure floquet_2Dh_axisym_matvec
            character(len=*), parameter :: this_procedure = 'floquet_2Dh_axisym_matvec'
            integer :: n
            select type (vec_in)
            type is (nek_zvector)
               select type (vec_out)
               type is (nek_zvector)
                  if (.not. self%is_initialized) call self%init()
                  n = bf_get_nsteps()
         ! A representative baseflow for the solver setup.
                  call install_snapshot(1)
         ! Ensure correct Nek status
                  call setup_linear_solver(transpose         = .false.,
     &                                     solve_baseflow    = .false.,
     &                                     solve_temperature = .true.,
     &                                     variable_dt       = .true.,
     &                                     endtime           = self%tau,
     &                                     cfl_limit         = cfl_floquet)
         ! Set initial condition for the linearized solver.
                  call vec2nek(vxp, vyp, vzp, prp, tp, vec_in)
         ! Replay the orbit: snapshot k and the recorded dt_k at step k.
                  call bf_replay_start()
                  time = 0.0_dp
                  do istep = 1, n
                     call bf_set(istep)
                     call nek_advance_2Dh_axisym(self%alpha)
                  end do
                  call bf_replay_end()
         ! Copy the final solution to vector.
                  call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)
               class default
                  call type_error('vec_out','nek_zvector','OUT',this_module,this_procedure)
               end select
            class default
               call type_error('vec_in','nek_zvector','IN',this_module,this_procedure)
            end select
         end procedure floquet_2Dh_axisym_matvec

      !====================================================================
      !     ADJOINT
      !====================================================================

         module procedure floquet_2Dh_axisym_rmatvec
            character(len=*), parameter :: this_procedure = 'floquet_2Dh_axisym_rmatvec'
            real(dp) :: dtk
            integer :: n, k
            select type (vec_in)
            type is (nek_zvector)
               select type (vec_out)
               type is (nek_zvector)
                  if (.not. self%is_initialized) call self%init()
                  n = bf_get_nsteps()
         ! A representative baseflow for the solver setup.
                  call install_snapshot(1)
         ! Ensure correct Nek status (transpose -> ifadj)
                  call setup_linear_solver(transpose         = .true.,
     &                                     solve_baseflow    = .false.,
     &                                     solve_temperature = .true.,
     &                                     variable_dt       = .true.,
     &                                     endtime           = self%tau,
     &                                     cfl_limit         = cfl_floquet)
         ! Set initial condition for the linearized solver.
                  call vec2nek(vxp, vyp, vzp, prp, tp, vec_in)
         ! Backward replay: adjoint step istep undoes forward step k, starting
         ! from U(t_k) = snapshot mod(k,n)+1, with the step length dt_k.
                  call bf_replay_start(reverse=.true.)
                  time = 0.0_dp
                  do istep = 1, n
                     k = n - istep + 1
                     dtk = bf_get_dt(k)
                     call bf_set(mod(k, n) + 1, if_lag=.false.)
                     param(12) = -abs(dtk)
                     call nek_advance_2Dh_axisym(self%alpha)   ! same alpha; adjoint terms switch on ifadj
                  end do
                  call bf_replay_end()
         ! Copy the final solution to vector.
                  call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)
               class default
                  call type_error('vec_out','nek_zvector','OUT',this_module,this_procedure)
               end select
            class default
               call type_error('vec_in','nek_zvector','IN',this_module,this_procedure)
            end select
         end procedure floquet_2Dh_axisym_rmatvec

      !====================================================================
      !     HELPER
      !====================================================================

         subroutine install_snapshot(k)
      !! Puts snapshot k into vx, vy, t without touching the lag levels, and
      !! leaves param(12) as it found it.
            integer, intent(in) :: k
            call bf_replay_start()
            call bf_set(k, if_lag=.false.)
            call bf_replay_end()
         end subroutine install_snapshot

      end submodule floquet_operator_2Dh