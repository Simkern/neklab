      !--------------------------------------------------------------------------
      !
      !  NEKLAB TIMING
      !
      !  A neklab-local watch derived from LightKrylov's `abstract_watch`, in the
      !  same way LightROM derives its own. It owns two things LightKrylov cannot
      !  see:
      !
      !    * the INSIDE of a 2Dh timestep (nek_advance_2Dh_axisym and the solvers
      !      and kernels it calls), and
      !    * the driver-level passes that are not system evaluations at all
      !      (the finite-difference Jacobian, the gauge shift, the save-orbit
      !      pass).
      !
      !  It deliberately does NOT time `response`/`matvec`/`rmatvec`: those carry
      !  their own timers inside `abstract_system_rdp` / `abstract_linop_rdp` and
      !  duplicating them here would only produce two numbers to reconcile. What
      !  this module does instead is read their COUNTERS (`get_eval_counter`,
      !  `get_counter`, both pure) so every CSV row can be joined to LightKrylov's
      !  own accounting.
      !
      !  NAMING CONVENTION
      !
      !    1. A timer that wraps exactly one procedure takes that procedure's
      !       name, lowercased. `add_timer` lowercases anyway, so `this_procedure`
      !       can be handed straight through, exactly as LightKrylov does in
      !       NewtonKrylov.f90.
      !
      !    2. A timer that wraps a REGION inside a procedure takes
      !       `<procedure>@<region>`. The `@` cannot collide with a Fortran name
      !       and sorts each region next to its parent in the summary. The set
      !       `nek_advance_2dh_axisym@*` is by construction a partition of
      !       `nek_advance_2dh_axisym`, so "do the parts add up" is a one-line
      !       check.
      !
      !  All timer names are declared as parameters below and used through those
      !  parameters at the call sites. A typo then fails at compile time rather
      !  than inside `stop_error` at run time, and the CSV column list stays in
      !  sync with the registry by construction.
      !
      !  THE CLOCK
      !
      !  LightKrylov's timers use `cpu_time`. With one thread per rank and a
      !  busy-waiting MPI that tracks wall time closely, but it is not the same
      !  quantity. `neklab_clock_tic`/`toc` therefore accumulate an INDEPENDENT
      !  wall (dnekclock) and cpu total over the body of nek_advance_2Dh_axisym,
      !  and finalize prints both. If cpu/wall drifts away from 1 on rank 0,
      !  every number in the table is suspect and that is worth seeing.
      !
      !  TWO OUTPUT CHANNELS
      !
      !    * neklab_timer.log  -- the watch's own summary, written by `finalize`.
      !      A soft `reset_all` at each outer Newton iteration turns this into a
      !      readable per-iteration table (one row per reset per timer).
      !
      !    * neklab_timings_2Dh.csv -- optional (off by default), one row per map
      !      evaluation, INCREMENTAL values (seconds and calls since the previous
      !      row). Incremental rather than cumulative because a soft reset zeroes
      !      the local elapsed time; differencing cumulative columns across a
      !      reset would give negative deltas.
      !
      !  CRITICAL CONSTRAINT ON THE CSV
      !
      !  `lightkrylov_timer%get_data` STOPS a running timer (Timer_Utils.f90:352).
      !  Polling a timer that is running at a dump point would silently truncate
      !  its interval, and `restart=.true.` would additionally inflate its call
      !  count. `csv_names` below therefore contains ONLY timers that are
      !  guaranteed idle at a dump point, i.e. step-level ones. Never add a
      !  driver-level timer to that list: `solve_fixed_point`, `flowrate_newton@iter`,
      !  `fd_jacobian`, `shift_mflow_phase_upo` and the save-orbit timer are all
      !  running while a map evaluation is in flight.
      !
      !--------------------------------------------------------------------------

      module neklab_timing
         use stdlib_optval, only: optval
         use LightKrylov, only: dp, abstract_system_rdp, abstract_linop_rdp
         use LightKrylov_Logger, only: log_message, log_information, log_debug
         use LightKrylov_Timer_Utils, only: abstract_watch
         use LightKrylov_Timing, only: global_lightkrylov_timer

         implicit none
         include "SIZE"
         include "TOTAL"

         private
         character(len=*), parameter, private :: this_module = 'neklab_timing'

      !--------------------------------------------------------------
      !     TIMER NAMES
      !--------------------------------------------------------------

      ! --- group: 2Dh timestepper. The `@` timers partition the routine.
         character(len=*), parameter, public :: t_advance = 'nek_advance_2dh_axisym'
         character(len=*), parameter, public :: t_adv_setup = 'nek_advance_2dh_axisym@setup'
         character(len=*), parameter, public :: t_adv_advection = 'nek_advance_2dh_axisym@advection'
         character(len=*), parameter, public :: t_adv_rhs = 'nek_advance_2dh_axisym@rhs'
         character(len=*), parameter, public :: t_adv_uz = 'nek_advance_2dh_axisym@uz'
         character(len=*), parameter, public :: t_adv_urphi = 'nek_advance_2dh_axisym@urphi'
         character(len=*), parameter, public :: t_adv_pressure = 'nek_advance_2dh_axisym@pressure'
         character(len=*), parameter, public :: t_adv_correction = 'nek_advance_2dh_axisym@correction'
         character(len=*), parameter, public :: t_makefp = 'makefp_2dh_axisym'
         character(len=*), parameter, public :: t_nl_step = 'nek_advance@nonlinear_2dh'
      !! Nek's own nonlinear step, timed at the call site because it is F77.

      ! --- group: 2Dh solvers
         character(len=*), parameter, public :: t_solve_hmh = 'solve_helmholtz_2dh_axisym'
         character(len=*), parameter, public :: t_solve_cpl = 'solve_coupled_helmholtz_2dh_axisym'
         character(len=*), parameter, public :: t_solve_pres = 'solve_pressure_2dh_axisym'
         character(len=*), parameter, public :: t_pres_precond = 'solve_pressure_2dh_axisym@precond'
         character(len=*), parameter, public :: t_pres_proj = 'pressure_projection_2dh_axisym'
      !! setrhs_ + gensoln_ + econj_, i.e. the whole mxprev projection.

      ! --- group: 2Dh kernels. Inclusive, and by design overlapping the two
      !     groups above: these are the drill-down, not part of any partition.
         character(len=*), parameter, public :: t_hmh_matvec = 'helmholtz_matvec_2dh_axisym'
         character(len=*), parameter, public :: t_cpl_matvec = 'coupled_helmholtz_matvec_2dh_axisym'
         character(len=*), parameter, public :: t_pres_matvec = 'pressure_matvec_2dh_axisym'
         character(len=*), parameter, public :: t_gs_comm = 'neklab_2dh_gs_comm'

      ! --- group: baseflow buffer
         character(len=*), parameter, public :: t_bf_push = 'neklab_bf_push'
         character(len=*), parameter, public :: t_bf_set = 'neklab_bf_set'
         character(len=*), parameter, public :: t_bf_window = 'neklab_bf_set_window'
         character(len=*), parameter, public :: t_bf_write = 'neklab_bf_write_chunk'
         character(len=*), parameter, public :: t_bf_read = 'neklab_bf_read_chunk'
         character(len=*), parameter, public :: t_close_mflow = 'close_mflow'

      ! --- group: drivers. NONE of these may appear in csv_names: they are
      !     running while a map evaluation is in flight.
         character(len=*), parameter, public :: t_fp_solve = 'solve_fixed_point'
         character(len=*), parameter, public :: t_fr_iter = 'flowrate_newton@iter'
         character(len=*), parameter, public :: t_fd_jac = 'fd_jacobian'
         character(len=*), parameter, public :: t_shift_phase = 'shift_mflow_phase_upo'
         character(len=*), parameter, public :: t_save_orbit = 'unsteady_flowrate_newton@save_orbit'

      !--------------------------------------------------------------
      !     THE WATCH
      !--------------------------------------------------------------

         type, extends(abstract_watch), public :: neklab_watch
      !! Global timing structure for all neklab-internal timers.
         contains
            private
            procedure, pass(self), public :: set_private_timers_and_name => set_neklab_timers
         end type neklab_watch

         type(neklab_watch), public :: global_neklab_timer

         logical, private :: if_time = .false.

      !--------------------------------------------------------------
      !     CSV CHANNEL STATE
      !--------------------------------------------------------------

         integer, parameter, private :: lname = 48
         integer, parameter, private :: ncsv = 25

      ! Step-level timers only. See the note at the top of this file: every name
      ! here MUST be idle whenever neklab_timer_dump is called.
         character(len=lname), parameter, private :: csv_names(ncsv) =
     &      [character(len=lname) ::
     &       t_advance, t_adv_setup, t_adv_advection, t_adv_rhs, t_adv_uz,
     &       t_adv_urphi, t_adv_pressure, t_adv_correction, t_makefp, t_nl_step,
     &       t_solve_hmh, t_solve_cpl, t_solve_pres, t_pres_precond, t_pres_proj,
     &       t_hmh_matvec, t_cpl_matvec, t_pres_matvec, t_gs_comm,
     &       t_bf_push, t_bf_set, t_bf_window, t_bf_write, t_bf_read,
     &       t_close_mflow]

         logical, private :: if_csv = .false.
         integer, private :: csv_unit = -1
         integer, private :: csv_row = 0
         real(dp), private :: csv_prev_s(ncsv) = 0.0_dp
         integer, private :: csv_prev_n(ncsv) = 0
         character(len=64), private :: csv_file = 'neklab_timings_2Dh.csv'

      ! --- row tag, set by the driver
         character(len=16), private :: tag_phase = 'init'
         integer, private :: tag_outer = 0

      !--------------------------------------------------------------
      !     CLOCK CROSS-CHECK STATE
      !--------------------------------------------------------------

         real(dp), private :: wall_acc = 0.0_dp
         real(dp), private :: cpu_acc = 0.0_dp
         real(dp), private :: wall_t0 = 0.0_dp
         real(dp), private :: cpu_t0 = 0.0_dp
         real(dp), private :: wall_ref = 0.0_dp
         logical, private :: clk_running = .false.

      !--------------------------------------------------------------
      !     PUBLIC INTERFACE
      !--------------------------------------------------------------

         public :: time_neklab, time_neklab_csv
         public :: neklab_initialize_timers, neklab_finalize_timers
         public :: neklab_reset_timers, neklab_enumerate_timers
         public :: neklab_timer_start, neklab_timer_stop
         public :: neklab_timer_tag, neklab_timer_dump
         public :: neklab_clock_tic, neklab_clock_toc

      contains

      !====================================================================
      !     SWITCHES
      !====================================================================

         pure logical function time_neklab() result(l)
            l = if_time
         end function time_neklab

         pure logical function time_neklab_csv() result(l)
            l = if_csv
         end function time_neklab_csv

         subroutine set_neklab_timer_switch(value)
            logical, intent(in) :: value
            if (if_time .neqv. value) then
               if_time = value
               if (if_time) then
                  call log_message('neklab timing enabled.', module=this_module)
               else
                  call log_message('neklab timing disabled.', module=this_module)
               end if
            else
               call log_debug('neklab timing switch unchanged.', module=this_module)
            end if
         end subroutine set_neklab_timer_switch

      !====================================================================
      !     GATED WRAPPERS
      !====================================================================
      !
      !  Used in preference to the inline `if (time_neklab()) call ...%start()`
      !  idiom of LightKrylov: one call site instead of two, and the branch cost
      !  is a few nanoseconds against regions that are microseconds at the very
      !  fastest. The inline idiom remains available through `time_neklab` and
      !  `global_neklab_timer`, both public.

         subroutine neklab_timer_start(name)
            character(len=*), intent(in) :: name
            if (if_time) call global_neklab_timer%start(name)
         end subroutine neklab_timer_start

         subroutine neklab_timer_stop(name)
            character(len=*), intent(in) :: name
            if (if_time) call global_neklab_timer%stop(name)
         end subroutine neklab_timer_stop

      !====================================================================
      !     TIMER REGISTRY
      !====================================================================

         subroutine set_neklab_timers(self)
      !! Define the private timers and their groups. Called once, from
      !! `initialize`. Registering here rather than at the point of first use
      !! is what makes `bf_init` re-entrant: `add_timer` aborts on a duplicate
      !! name (Timer_Utils.f90:442), so a second `bf_init` in the same run used
      !! to kill the job.
            class(neklab_watch), intent(inout) :: self
      ! internal
            integer :: istart, iend

            call self%set_watch_name('neklab_timer')

      ! --- 2Dh timestepper
            call self%add_timer(t_advance, count=istart)
            call self%add_timer(t_adv_setup)
            call self%add_timer(t_adv_advection)
            call self%add_timer(t_adv_rhs)
            call self%add_timer(t_adv_uz)
            call self%add_timer(t_adv_urphi)
            call self%add_timer(t_adv_pressure)
            call self%add_timer(t_adv_correction)
            call self%add_timer(t_makefp)
            call self%add_timer(t_nl_step, count=iend)
            call self%add_group('2Dh_timestepper', istart=istart, iend=iend)

      ! --- 2Dh solvers
            call self%add_timer(t_solve_hmh, count=istart)
            call self%add_timer(t_solve_cpl)
            call self%add_timer(t_solve_pres)
            call self%add_timer(t_pres_precond)
            call self%add_timer(t_pres_proj, count=iend)
            call self%add_group('2Dh_solvers', istart=istart, iend=iend)

      ! --- 2Dh kernels
            call self%add_timer(t_hmh_matvec, count=istart)
            call self%add_timer(t_cpl_matvec)
            call self%add_timer(t_pres_matvec)
            call self%add_timer(t_gs_comm, count=iend)
            call self%add_group('2Dh_kernels', istart=istart, iend=iend)

      ! --- baseflow buffer
            call self%add_timer(t_bf_push, count=istart)
            call self%add_timer(t_bf_set)
            call self%add_timer(t_bf_window)
            call self%add_timer(t_bf_write)
            call self%add_timer(t_bf_read)
            call self%add_timer(t_close_mflow, count=iend)
            call self%add_group('bf_buffer', istart=istart, iend=iend)

      ! --- drivers
            call self%add_timer(t_fp_solve, count=istart)
            call self%add_timer(t_fr_iter)
            call self%add_timer(t_fd_jac)
            call self%add_timer(t_shift_phase)
            call self%add_timer(t_save_orbit, count=iend)
            call self%add_group('drivers_2Dh', istart=istart, iend=iend)

            call set_neklab_timer_switch(.true.)
         end subroutine set_neklab_timers

      !====================================================================
      !     LIFECYCLE
      !====================================================================

         subroutine neklab_initialize_timers(if_csv_out, csv_filename)
      !! Switch on timing for both LightKrylov and neklab, and optionally open
      !! the CSV channel. Call from the .usr file (userchk, istep == 0) BEFORE
      !! any solve and before the first `bf_init`.
      !!
      !! Note the asymmetry inherited from LightKrylov: `set_lightkrylov_timer_switch`
      !! is private, so timing switches on the moment `initialize` runs and
      !! cannot be switched off again. "Timing off" means "do not call this".
            logical, optional, intent(in) :: if_csv_out
      !! Write the per-map-evaluation CSV? Default .false.
            character(len=*), optional, intent(in) :: csv_filename
      ! internal
            real(dp), external :: dnekclock

            call global_lightkrylov_timer%initialize()
            call global_neklab_timer%initialize()

            wall_acc = 0.0_dp
            cpu_acc = 0.0_dp
            clk_running = .false.
            wall_ref = dnekclock()

            if (present(csv_filename)) csv_file = csv_filename
            if (optval(if_csv_out, .false.)) call open_csv()
         end subroutine neklab_initialize_timers

         subroutine neklab_finalize_timers(system, linop)
      !! Close the CSV, report the clock cross-check, and print all summaries.
      !! Pass the system and/or the Jacobian to have LightKrylov's own
      !! per-evaluation and per-matvec timers printed alongside, so the cost of
      !! a map and the cost of the steps inside it appear in the same output.
            class(abstract_system_rdp), optional, intent(inout) :: system
            class(abstract_linop_rdp), optional, intent(inout) :: linop
      ! internal
            character(len=128) :: msg

            if (if_csv) call close_csv()

            if (if_time) then
               call log_message('#########   neklab clock cross-check   ##########', this_module)
               write (msg, '(2X,A30," : ",F14.6," s")') 'nek_advance_2Dh_axisym (wall)', wall_acc
               call log_message(msg, this_module)
               write (msg, '(2X,A30," : ",F14.6," s")') 'nek_advance_2Dh_axisym (cpu)', cpu_acc
               call log_message(msg, this_module)
               if (wall_acc > 0.0_dp) then
                  write (msg, '(2X,A30," : ",F14.6)') 'cpu/wall (rank 0)', cpu_acc/wall_acc
                  call log_message(msg, this_module)
                  call log_message('  A ratio far from 1 means cpu_time is not a proxy for wall '//
     &               'time here and every number below is suspect.', this_module)
               end if
            end if

            if (present(system)) call system%finalize_timer()
            if (present(linop)) call linop%finalize_timer()

            call global_lightkrylov_timer%finalize()
            call global_neklab_timer%finalize()
         end subroutine neklab_finalize_timers

         subroutine neklab_reset_timers(soft, clean, also_lightkrylov)
      !! Soft reset of the neklab watch: snapshots the current data into the
      !! per-reset history and zeroes the live counters. Called at the end of
      !! each outer Newton iteration, which is what turns the finalize summary
      !! into a per-iteration table.
      !!
      !! A residual CSV row is emitted first so that whatever happened between
      !! the last map evaluation and the reset is not lost. `reset_all` refuses
      !! to reset a RUNNING timer and logs a message when it meets one, which is
      !! why the driver-level timers are stopped before this is called.
            logical, optional, intent(in) :: soft
            logical, optional, intent(in) :: clean
            logical, optional, intent(in) :: also_lightkrylov
      !! Reset LightKrylov's watch too? Default .false. -- only safe where no
      !! LightKrylov timer is running, i.e. outside `newton`.
            if (.not. if_time) return
            if (if_csv) call neklab_timer_dump('-')
            call global_neklab_timer%reset_all(soft, clean)
            if (optval(also_lightkrylov, .false.)) then
               call global_lightkrylov_timer%reset_all(soft, clean)
            end if
      ! The live elapsed times have just been zeroed, so the CSV baseline must
      ! be zeroed with them or the next row would show a negative increment.
            csv_prev_s = 0.0_dp
            csv_prev_n = 0
         end subroutine neklab_reset_timers

         subroutine neklab_enumerate_timers(only_user)
            logical, optional, intent(in) :: only_user
            call global_lightkrylov_timer%enumerate(only_user)
            call global_neklab_timer%enumerate(only_user)
         end subroutine neklab_enumerate_timers

      !====================================================================
      !     CSV CHANNEL
      !====================================================================

         subroutine neklab_timer_tag(phase, outer)
      !! Label the rows that follow. Without this the finite-difference,
      !! gauge and save-orbit passes are indistinguishable from the Newton
      !! passes they happen to sit next to.
            character(len=*), optional, intent(in) :: phase
            integer, optional, intent(in) :: outer
            if (present(phase)) tag_phase = phase
            if (present(outer)) tag_outer = outer
         end subroutine neklab_timer_tag

         subroutine neklab_timer_dump(evaltype, eval, matvec, rmatvec, nsteps)
      !! Write one CSV row: INCREMENTAL seconds and calls for every timer in
      !! `csv_names` since the previous row.
      !!
      !! Call at the END of a map evaluation, where every timer in `csv_names`
      !! is idle. `get_data` stops a running timer, so a dump from anywhere else
      !! would corrupt the very data it reads.
            character(len=*), intent(in) :: evaltype
      !! 'F', 'J', 'JT' or '-' for the residual row emitted before a reset.
            integer, optional, intent(in) :: eval
      !! system%get_eval_counter() -- the join key to LightKrylov's accounting.
            integer, optional, intent(in) :: matvec
      !! linop%get_counter(.false.)
            integer, optional, intent(in) :: rmatvec
      !! linop%get_counter(.true.)
            integer, optional, intent(in) :: nsteps
      !! Steps in this pass. Cost PER STEP is the invariant; cost per map is
      !! not, because dt is CFL-adaptive.
      ! internal
            real(dp) :: etm, ds(ncsv), wnow
            real(dp), external :: dnekclock
            integer :: iw, lcnt, dn(ncsv)

            if (.not. if_csv) return
            if (nid /= 0) return

            do iw = 1, ncsv
               call global_neklab_timer%get_data(csv_names(iw), etime=etm, lcount=lcnt)
               ds(iw) = etm - csv_prev_s(iw)
               dn(iw) = lcnt - csv_prev_n(iw)
      ! Defensive: a reset between two rows zeroes the live counters. The
      ! baseline is zeroed with it in neklab_reset_timers, so this should never
      ! fire -- but a silently negative increment would be worse than a clamp.
               if (ds(iw) < 0.0_dp) ds(iw) = etm
               if (dn(iw) < 0) dn(iw) = lcnt
               csv_prev_s(iw) = etm
               csv_prev_n(iw) = lcnt
            end do

            wnow = dnekclock() - wall_ref
            csv_row = csv_row + 1

      ! The colon edit descriptor is what stops the trailing separator: without
      ! it the unlimited group re-enters after the last I0, writes the literal
      ! comma, and only then finds the item list exhausted -- leaving a dangling
      ! comma that makes every reader see one extra column.
            write (csv_unit, '(I0,2(",",A),4(",",I0),",",I0,",",F0.6,",",*(ES14.7,",",I0,:,","))')
     &         csv_row, trim(tag_phase), trim(evaltype), tag_outer,
     &         optval(eval, -1), optval(matvec, -1), optval(rmatvec, -1),
     &         optval(nsteps, -1), wnow, (ds(iw), dn(iw), iw=1, ncsv)
            flush (csv_unit)
         end subroutine neklab_timer_dump

         subroutine open_csv()
            character(len=4096) :: hdr
            character(len=lname) :: cstem
            integer :: iw, io_stat

            if (nid == 0) then
               open (newunit=csv_unit, file=trim(csv_file), status='replace',
     &               action='write', iostat=io_stat)
               if (io_stat /= 0) then
                  call log_message('Could not open '//trim(csv_file)//'. CSV channel disabled.',
     &               this_module, 'open_csv')
                  if_csv = .false.
                  return
               end if

      ! metadata comment: pandas.read_csv(..., comment='#') skips it
      ! lpert as well as npert: this runs from userchk before any solve, so the
      ! run-time npert is usually still 0 and only the compiled slot count is
      ! meaningful at this point.
               write (hdr, '(A,I0,3(A,I0),3(A,I0))')
     &            '# neklab 2Dh timings | clock=cpu_time | np=', np,
     &            ' | lx1=', lx1, ' ly1=', ly1, ' lx2=', lx2,
     &            ' | nelv=', nelv, ' npert=', npert, ' lpert=', lpert
               write (csv_unit, '(A)') trim(hdr)
               write (csv_unit, '(A)') '# values are INCREMENTAL: seconds and calls since the previous row'

      ! header
               hdr = 'row,phase,evaltype,outer,eval,matvec,rmatvec,nsteps,wall'
               do iw = 1, ncsv
                  cstem = csv_col(csv_names(iw))
                  hdr = trim(hdr)//','//trim(cstem)//'_s,'//trim(cstem)//'_n'
               end do
               write (csv_unit, '(A)') trim(hdr)
               flush (csv_unit)
            end if

            if_csv = .true.
            csv_row = 0
            csv_prev_s = 0.0_dp
            csv_prev_n = 0
            call log_message('neklab CSV timing channel open: '//trim(csv_file),
     &         this_module, 'open_csv')
         end subroutine open_csv

         subroutine close_csv()
            if (nid == 0 .and. csv_unit > 0) close (csv_unit)
            csv_unit = -1
            if_csv = .false.
         end subroutine close_csv

         pure function csv_col(name) result(col)
      !! Timer name -> CSV column stem. The `@` of a region timer becomes `__`,
      !! which keeps the column a legal identifier for whatever reads it.
            character(len=*), intent(in) :: name
            character(len=lname) :: col
            integer :: iw, jw
            col = ' '
            jw = 0
            do iw = 1, len_trim(name)
               if (name(iw:iw) == '@') then
                  col(jw + 1:jw + 2) = '__'
                  jw = jw + 2
               else
                  jw = jw + 1
                  col(jw:jw) = name(iw:iw)
               end if
            end do
         end function csv_col

      !====================================================================
      !     WALL/CPU CROSS-CHECK
      !====================================================================

         subroutine neklab_clock_tic()
      !! Independent wall and cpu accumulation, used ONLY around the body of
      !! nek_advance_2Dh_axisym. Not affected by resets, so the totals reported
      !! at finalize cover the whole run.
            real(dp), external :: dnekclock
            if (.not. if_time) return
            if (clk_running) return
            wall_t0 = dnekclock()
            call cpu_time(cpu_t0)
            clk_running = .true.
         end subroutine neklab_clock_tic

         subroutine neklab_clock_toc()
            real(dp) :: cnow
            real(dp), external :: dnekclock
            if (.not. if_time) return
            if (.not. clk_running) return
            call cpu_time(cnow)
            wall_acc = wall_acc + (dnekclock() - wall_t0)
            cpu_acc = cpu_acc + (cnow - cpu_t0)
            clk_running = .false.
         end subroutine neklab_clock_toc

      end module neklab_timing