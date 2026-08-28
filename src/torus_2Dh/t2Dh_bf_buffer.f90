      module t2Dh_bf_buffer
      !! Time-resolved baseflow buffer for periodic-orbit Newton solves on the
      !! 2Dh (axisymmetric torus) mesh.
      !!
      !! During the NONLINEAR pass the trajectory is recorded snapshot by
      !! snapshot; during the subsequent LINEAR passes (one per GMRES matvec)
      !! it is replayed. Because the baseflow and the perturbation live on the
      !! SAME 2D mesh, a snapshot is simply (vx, vy, u_phi, dt) and replaying it
      !! is a plain copy -- no streamwise reconstruction is needed, unlike the
      !! 3D helix case this module descends from.
      !!
      !! TIMING CONVENTION (this is the part that is easy to get wrong)
      !!
      !!    slot k  holds the state at t^{k-1}, i.e. the state ENTERING step k
      !!    dt(k)   holds the timestep that step k actually consumed
      !!
      !! Nek's settime() shifts dtlag, calls setdt (which sets dt for the step
      !! about to be taken) and only then advances time. Reading dt before
      !! nek_advance therefore returns the PREVIOUS step's value. The recording
      !! is consequently split in two:
      !!
      !!    call bf_begin_step()   ! state, before nek_advance
      !!    call nek_advance()
      !!    call bf_end_step()     ! dt, after nek_advance
      !!
      !! so that sum(dt(1:nsteps)) == T to round-off. (helix_2d.f90 stores dt
      !! inside save_2d_fields, which runs BEFORE nek_advance and therefore
      !! records step k-1's timestep in slot k. That off-by-one is fixed here.)
      !!
      !! MEMORY
      !!
      !! lbuf is a compile-time constant, as for lelv / lpert. The per-rank
      !! arrays are allocated on nelv (not lelv) at bf_init, so the footprint
      !! follows the actual partition rather than the serial worst case. At most
      !! ONE chunk is ever resident: long trajectories simply produce more
      !! files.
      !!
      !! FILE FORMAT
      !!
      !! Byte-identical to the helix 2D format (fileversion '1'), so the
      !! existing Python suite (read_2d_data.py / write_2d_data.py) reads these
      !! files unchanged. The three field slots written per snapshot are
      !!
      !!    slot 1 = vx  = u_Z    (axial,     along the torus axis)
      !!    slot 2 = vy  = u_R    (radial,    from the torus axis)
      !!    slot 3 = vz  = u_phi  = t(:,1)  (azimuthal / streamwise)
      !!
      !! i.e. the streamwise component is stored in the 'vz' slot. Pressure is
      !! NOT stored: the linearised operator does not consume it, and pressure
      !! for the Newton iterates is written separately with standard outposting.
      !!
      !! Writing is serial (gather to nid 0). The element map written to file is
      !! the GLOBAL element numbering, so files are portable across a change in
      !! the number of MPI ranks.
      !!
      !! NOTE ON NAMING: the parameter lbuf below shadows the lbuf of
      !! neklab_helix. Any scoping unit that uses both modules must import at
      !! least one of them with an explicit only-list.
         use stdlib_optval, only: optval
         use stdlib_strings, only: padl, padr
         use stdlib_sorting, only: sort_index
         use LightKrylov, only: dp
         use LightKrylov_Logger
         use neklab_timing, only: neklab_timer_start, neklab_timer_stop,
     &                            t_bf_push, t_bf_set, t_bf_window,
     &                            t_bf_write, t_bf_read
         use neklab_nek_setup, only: nek_log_message, nek_log_information,
     &                               nek_log_warning, nek_log_debug, nek_stop_error
         implicit none
         include "SIZE"
         include "TOTAL"
         include "RESTART"
      !! RESTART supplies if_byte_sw, used by the reader in the submodule.
         private
         character(len=*), parameter, private :: this_module = 'neklab_bf_buffer'

         integer, parameter, private :: lv = lx1*ly1*lz1*lelv
      !! Worst-case per-rank field length. Only used to size the static gather
      !! scratch in the IO submodule; the buffer itself is allocated on nelv.

         integer, parameter, public :: lbuf = 1000
      !! Maximum number of snapshots held in memory before the chunk is spilled
      !! to disk. Compile-time, same convention as lelv / lpert.

      !--------------------------------------------------------------------
      !-----     BUFFER STATE                                         -----
      !--------------------------------------------------------------------

         real(dp), allocatable, private :: bufx(:, :), bufy(:, :), bufs(:, :)
      !! (nbf, lbuf) snapshot storage for u_Z, u_R, u_phi.
         real(dp), allocatable, private :: bufdt(:)
      !! (lbuf) timestep consumed by each buffered step.

         real(dp), allocatable, private :: x2d(:), y2d(:)
      !! (nbf) reference coordinates, written to every file for verification.

         integer, private :: nbf = 0
      !! Per-rank field length, lx1*ly1*lz1*nelv. Set at bf_init.

         integer, private :: nsave = 0
      !! Number of snapshots currently held in the memory chunk.
         integer, private :: nsteps_rec = 0
      !! Total number of steps recorded over the whole trajectory.
         integer, private :: nchunk = 0
      !! Number of chunks produced by the current recording.
         integer, private :: ichunk = 0
      !! Index of the chunk currently resident in memory (0 = none).
         integer, private :: nelf = 0
      !! Global number of 2D elements (= nelgv).
         real(dp), private :: trec = 0.0_dp
      !! Running sum of the recorded timesteps.

         real(dp), private :: dt_min = huge(1.0_dp)
         real(dp), private :: dt_max = 0.0_dp
      !! Extremes of the recorded timestep, EXCLUDING the steps the caller
      !! flags as non-physical (the landing steps, whose length is set by the
      !! period constraint rather than by the CFL condition).

         integer, private :: nsteps_min = 50
      !! A period resolved by fewer steps than this does not represent the
      !! dynamics. bf_close_record fails rather than continuing quietly.

         real(dp), private :: p12_save = 0.0_dp
         logical, private :: replaying = .false.
         logical, private :: replay_rev = .false.

         logical, private :: is_initialized = .false.
         logical, private :: is_recording = .false.
         logical, private :: if_write = .true.
      !! Spill every chunk to disk even when the whole trajectory fits in
      !! memory. Costs one write per nonlinear solve and buys restartability.
         logical, private :: if_sym = .false.
      !! Half (symmetric) cross-section mesh? Auto-detected from the SYM bcs.

         character(len=32), private :: fbase = '2dtorus'
         character(len=1), private :: fprefix = 'n'

         public :: bf_init, bf_finalize_module
         public :: bf_reset, bf_begin_step, bf_end_step, bf_close_record
         public :: bf_replay_start, bf_replay_end, bf_set, bf_set_window
         public :: bf_get_nsteps, bf_get_dt, bf_get_dt_minmax, bf_get_nchunk
         public :: bf_get_time, bf_is_recording, bf_is_initialised, bf_get_nsteps_min
         public :: bf_set_prefix, bf_get_prefix, bf_fname
         public :: bf_write_mesh, bf_summary
         public :: bf_write_t2Dh, bf_read_t2Dh

      !--------------------------------------------------------------------
      !-----     IO INTERFACES (implemented in submodule bf_buffer_io) -----
      !--------------------------------------------------------------------
      ! The submodule has host access to every private entity above, so no
      ! accessors are needed.

         interface

            module subroutine bf_write_chunk(fname, only_mesh)
            !! Write the resident memory chunk to a binary file. Serial IO:
            !! all ranks ship their elements to nid 0, which writes.
               character(len=132), intent(in) :: fname
               logical, optional, intent(in) :: only_mesh
            !! Write coordinates only, with no field data (default .false.).
            end subroutine bf_write_chunk

            module subroutine bf_read_chunk(fname, nread)
            !! Read a binary chunk file into the memory buffer. nid 0 reads and
            !! broadcasts; each rank extracts the elements it owns using the
            !! global element map, so the file is rank-count portable.
               character(len=132), intent(in) :: fname
               integer, intent(out) :: nread
            !! Number of snapshots found in the file.
            end subroutine bf_read_chunk

            module subroutine bf_peek_chunk(fname, nread, ierr)
            !! Read only the header of a chunk file (no field data).
               character(len=132), intent(in) :: fname
               integer, intent(out) :: nread
               integer, intent(out) :: ierr
            !! Non-zero if the file cannot be opened or is malformed.
            end subroutine bf_peek_chunk

         end interface

      contains

      !====================================================================
      !     SETUP
      !====================================================================

         subroutine bf_init(base, write_chunks, min_steps)
      !! Allocates the buffer and caches the mesh metadata. Call once, after
      !! the mesh and the boundary conditions are available. The DRIVER owns
      !! this call (and bf_finalize_module): the system's nonlinear map only
      !! consumes the buffer, and fails loudly through check_init if the driver
      !! forgot.
            character(len=*), optional, intent(in) :: base
      !! Filename stem, appended to the one-character prefix. Default '2dtorus'.
            logical, optional, intent(in) :: write_chunks
      !! Spill chunks to disk even when the trajectory fits in memory.
      !! Default .true. (needed for restarts and for the stability runs).
            integer, optional, intent(in) :: min_steps
      !! Minimum acceptable number of steps per recorded horizon. Default 50.
      ! internal
            character(len=*), parameter :: this_procedure = 'bf_init'
            character(len=256) :: msg
            integer :: nsym, ie, iface
            real(dp) :: mbytes
            integer, parameter :: pad = 20
            integer, external :: iglsum

            if (is_initialized) then
               call nek_log_warning('Buffer already initialized. Ignoring.', this_module, this_procedure)
               return
            end if

            if (present(base)) fbase = base
            if_write = optval(write_chunks, .true.)
            nsteps_min = optval(min_steps, 50)

      ! This buffer is for the 2Dh reduction: a genuinely 2D mesh.
            if (if3d) then
               call nek_stop_error('neklab_bf_buffer is for 2D (2Dh) meshes only.', this_module, this_procedure)
            end if

            nelf = nelgv
            nbf = lx1*ly1*lz1*nelv

      ! Half or full cross-section? Detected exactly as helix_pipe does.
            nsym = 0
            do ie = 1, nelv
               do iface = 1, 2*ndim
                  if (cbc(iface, ie, 1) == 'SYM') nsym = nsym + 1
               end do
            end do
            nsym = iglsum(nsym, 1)
            if_sym = (nsym > 0)

      ! Allocated on nelv, not lelv: with lelv = lelg/lpmin + 3 the static size
      ! is the SERIAL worst case and overshoots by an order of magnitude on a
      ! realistic partition.
            allocate (bufx(nbf, lbuf)); call rzero(bufx, nbf*lbuf)
            allocate (bufy(nbf, lbuf)); call rzero(bufy, nbf*lbuf)
            allocate (bufs(nbf, lbuf)); call rzero(bufs, nbf*lbuf)
            allocate (bufdt(lbuf)); call rzero(bufdt, lbuf)
            allocate (x2d(nbf)); call copy(x2d, xm1, nbf)
            allocate (y2d(nbf)); call copy(y2d, ym1, nbf)

            is_initialized = .true.

            mbytes = 3.0_dp*real(nbf, dp)*real(lbuf, dp)*8.0_dp/1024.0_dp**2
            call nek_log_message('', this_module, this_procedure)
            call nek_log_message('Baseflow buffer initialized:', this_module, this_procedure)
            write (msg, '(3X,A,2X,3(1X,I4))') padl('lx1, ly1, nelv:', pad), lx1, ly1, nelv
            call nek_log_message(msg, this_module, this_procedure)
            write (msg, '(3X,A,1X,I16)') padl('snapshots per chunk:', pad), lbuf
            call nek_log_message(msg, this_module, this_procedure)
            write (msg, '(3X,A,1X,F13.1,A)') padl('footprint per rank:', pad), mbytes, ' MB'
            call nek_log_message(msg, this_module, this_procedure)
            write (msg, '(3X,A,1X,I16)') padl('global 2D elements:', pad), nelf
            call nek_log_message(msg, this_module, this_procedure)
            write (msg, '(3X,A,1X,I16)') padl('min steps / horizon:', pad), nsteps_min
            call nek_log_message(msg, this_module, this_procedure)
            write (msg, '(3X,A,1X,A16)') padl('cross-section:', pad), merge('half', 'full', if_sym)
            call nek_log_message(msg, this_module, this_procedure)
            write (msg, '(3X,A,1X,A16)') padl('chunk spilling:', pad), merge(' on', 'off', if_write)
            call nek_log_message(msg, this_module, this_procedure)
            call nek_log_message('', this_module, this_procedure)
         end subroutine bf_init

         subroutine bf_finalize_module()
            if (allocated(bufx)) deallocate (bufx)
            if (allocated(bufy)) deallocate (bufy)
            if (allocated(bufs)) deallocate (bufs)
            if (allocated(bufdt)) deallocate (bufdt)
            if (allocated(x2d)) deallocate (x2d)
            if (allocated(y2d)) deallocate (y2d)
            is_initialized = .false.
            is_recording = .false.
            replaying = .false.
         end subroutine bf_finalize_module

         subroutine bf_set_prefix(prefix)
      !! One-character file prefix. Convention inherited from the helix runs:
      !!    'n' Newton (working buffer, overwritten every nonlinear solve)
      !!    'b' converged orbit (written once, stable input for Floquet)
      !!    'f' Floquet
            character(len=1), intent(in) :: prefix
            fprefix = prefix
         end subroutine bf_set_prefix

         pure function bf_get_prefix() result(prefix)
            character(len=1) :: prefix
            prefix = fprefix
         end function bf_get_prefix

         pure function bf_fname(iout, prefix) result(fname)
      !! Standard chunk filename, e.g. 'n2dtorus001.fld'.
            integer, intent(in) :: iout
            character(len=1), optional, intent(in) :: prefix
            character(len=132) :: fname
      ! internal
            character(len=1) :: p
            p = fprefix
            if (present(prefix)) p = prefix
            fname = ' '
            write (fname, '(A,A,I3.3,A)') p, trim(fbase), iout, '.fld'
         end function bf_fname

      !====================================================================
      !     RECORDING
      !====================================================================

         subroutine bf_reset(prefix)
      !! Starts a new recording. Chunk files from a previous recording with the
      !! same prefix are overwritten as the new chunks are produced.
            character(len=1), optional, intent(in) :: prefix
      ! internal
            character(len=*), parameter :: this_procedure = 'bf_reset'
            if (.not. is_initialized) return
            if (present(prefix)) fprefix = prefix
            nsave = 0
            nsteps_rec = 0
            nchunk = 0
            ichunk = 0
            trec = 0.0_dp
            dt_min = huge(1.0_dp)
            dt_max = 0.0_dp
            call rzero(bufdt, lbuf)
            is_recording = .true.
            call nek_log_debug('Recording started.', this_module, this_procedure)
         end subroutine bf_reset

         subroutine bf_begin_step()
      !! Stores the state ENTERING the current step. Call immediately before
      !! nek_advance(). Spills the chunk first if the buffer is full.
            character(len=*), parameter :: this_procedure = 'bf_begin_step'
            if (.not. is_initialized) return
            if (.not. is_recording) then
               call nek_stop_error('bf_begin_step called outside a recording.', this_module, this_procedure)
            end if
            call neklab_timer_start(t_bf_push)

      ! Buffer full: flush before overwriting. bufdt is complete at this point
      ! because bf_end_step filled slot lbuf on the previous step.
            if (nsave == lbuf) then
               call flush_chunk(keep_resident=.false.)
            end if

            nsave = nsave + 1
            call copy(bufx(1, nsave), vx, nbf)
            call copy(bufy(1, nsave), vy, nbf)
            call copy(bufs(1, nsave), t(1, 1, 1, 1, 1), nbf)

            call neklab_timer_stop(t_bf_push)
         end subroutine bf_begin_step

         subroutine bf_end_step(count_stats)
      !! Stores the timestep the step just took. Call immediately after
      !! nek_advance(). Reading dt any earlier gives the PREVIOUS step's value.
            logical, optional, intent(in) :: count_stats
      !! Include this step in the min/max dt statistics. Default .true. Pass
      !! .false. for the landing steps, whose length is dictated by the period
      !! constraint and not by the CFL condition, so that the statistics keep
      !! describing the physical time resolution. The final step is excluded
      !! unconditionally for the same reason.
      ! internal
            character(len=*), parameter :: this_procedure = 'bf_end_step'
            if (.not. is_initialized) return
            if (nsave == 0) then
               call nek_stop_error('bf_end_step without a matching bf_begin_step.', this_module, this_procedure)
            end if
            bufdt(nsave) = dt
            nsteps_rec = nsteps_rec + 1
            trec = trec + dt
            if (optval(count_stats, .true.) .and. lastep == 0) then
               dt_min = min(dt_min, dt)
               dt_max = max(dt_max, dt)
            end if
         end subroutine bf_end_step

         subroutine bf_close_record(period)
      !! Ends the recording and flushes the partial chunk. If the horizon is
      !! supplied, the recorded timesteps are checked against it: a mismatch
      !! means the replay would run on a different time grid than the
      !! trajectory it is meant to linearise.
            real(dp), optional, intent(in) :: period
      ! internal
            character(len=*), parameter :: this_procedure = 'bf_close_record'
            character(len=256) :: msg
            real(dp) :: terr
            if (.not. is_initialized) return
            if (.not. is_recording) then
               call nek_log_warning('bf_close_record without an active recording.', this_module, this_procedure)
               return
            end if

            if (nsave > 0) call flush_chunk(keep_resident=.true.)
            is_recording = .false.

            write (msg, '(A,I0,A,I0,A,E16.8)') 'Recorded ', nsteps_rec, ' steps in ', nchunk,
     &         ' chunk(s), sum(dt)= ', trec
            call nek_log_message(msg, this_module, this_procedure)

      ! A horizon resolved by a handful of steps is not a discretisation of the
      ! dynamics, it is noise with a period attached. Fail rather than let a
      ! Newton solve chase it.
            if (nsteps_rec < nsteps_min) then
               write (msg, '(A,I0,A,I0,A)') 'Only ', nsteps_rec, ' steps in the horizon, minimum is ',
     &            nsteps_min, '. Lower the CFL target (param(26) / cfl_limit) or the period.'
               call nek_stop_error(msg, this_module, this_procedure)
            end if

            if (present(period)) then
               terr = abs(trec - period)
               if (terr > 1.0e-10_dp*max(abs(period), 1.0_dp)) then
                  write (msg, '(A,E16.8,A,E16.8)') 'sum(dt)= ', trec, ' vs T= ', period
                  call nek_log_message(msg, this_module, this_procedure)
                  call nek_stop_error('sum(dt) does not match the period: the replay would run '//
     &               'on a different time grid than the trajectory.', this_module, this_procedure)
               end if
            end if
         end subroutine bf_close_record

      !====================================================================
      !     REPLAY
      !====================================================================

         subroutine bf_replay_start(reverse)
      !! Prepares the buffer for a linear pass and saves the caller's param(12),
      !! which bf_set / bf_set_window overwrite to force the recorded timestep.
            logical, optional, intent(in) :: reverse
      !! Traverse the trajectory backwards (adjoint). Default .false.
      ! internal
            character(len=*), parameter :: this_procedure = 'bf_replay_start'
            character(len=256) :: msg
            call check_init(this_procedure)
            if (is_recording) then
               call nek_stop_error('Cannot replay while recording. Call bf_close_record first.',
     &            this_module, this_procedure)
            end if
            if (nsteps_rec == 0) then
               call nek_stop_error('Nothing recorded. Run the nonlinear map first.', this_module, this_procedure)
            end if
            replay_rev = optval(reverse, .false.)
            replaying = .true.
            p12_save = param(12)
            if (replay_rev) then
               call load_chunk(nchunk)
            else
               call load_chunk(1)
            end if
            write (msg, '(A,I0,A,I0,A,A)') 'Replay ', nsteps_rec, ' steps, ', nchunk, ' chunk(s), ',
     &         merge('backward', 'forward ', replay_rev)
            call nek_log_debug(msg, this_module, this_procedure)
         end subroutine bf_replay_start

         subroutine bf_replay_end()
      !! Restores param(12). Without this the fixed-timestep value left behind
      !! by the last bf_set leaks into whatever runs next.
            if (.not. replaying) return
            param(12) = p12_save
            replaying = .false.
         end subroutine bf_replay_end

         subroutine bf_set(k, if_lag)
      !! Installs snapshot k as the baseflow for step k and forces the step to
      !! use the timestep the nonlinear pass used. FORWARD traversal only.
      !!
      !! The lag levels are handled by nek's own shift routines rather than by
      !! storing them: on entry vx still holds snapshot k-1, so lagvel/lagscal
      !! push it down one level before the new snapshot is loaded, giving
      !!
      !!    vx = buf(k),  vxlag(1) = buf(k-1),  vxlag(2) = buf(k-2)
      !!
      !! exactly. At k = 1 the shift is skipped: settime sets nab = min(istep,3)
      !! so the (stale) lag levels are unreachable on the first two steps, and
      !! the nonlinear pass ramped its order in the same way.
      !!
      !! This is correct only because nek_advance_2Dh_axisym touches the
      !! PERTURBATION lag levels (lagfieldp / lagscalp) and never lagvel /
      !! lagscal, so there is no double shift. For backward traversal use
      !! bf_set_window, which builds the lag levels explicitly.
            integer, intent(in) :: k
      !! Step index, 1 .. nsteps.
            logical, optional, intent(in) :: if_lag
      !! Perform the lag shift (default .true.). Set .false. only if the caller
      !! manages the lag levels itself.
      ! internal
            character(len=*), parameter :: this_procedure = 'bf_set'
            character(len=256) :: msg
            integer :: c, j, ifield_save
            logical :: if_lag_
            call check_init(this_procedure)
            if (k < 1 .or. k > nsteps_rec) then
               write (msg, '(A,I0,A,I0,A)') 'Requested step ', k, ' outside [1, ', nsteps_rec, '].'
               call nek_stop_error(msg, this_module, this_procedure)
            end if
            if_lag_ = optval(if_lag, .true.)
            call neklab_timer_start(t_bf_set)

            c = (k - 1)/lbuf + 1
            j = k - (c - 1)*lbuf
            if (c /= ichunk) call load_chunk(c)

      ! Shift the CURRENT baseflow into the lag levels before overwriting it.
            if (if_lag_ .and. k > 1) then
               call lagvel
               ifield_save = ifield
               ifield = 2
               call lagscal
               ifield = ifield_save
            end if

            call copy(vx, bufx(1, j), nbf)
            call copy(vy, bufy(1, j), nbf)
            call copy(t(1, 1, 1, 1, 1), bufs(1, j), nbf)

      ! Force the step size. setdt latches its internal iffxdt as soon as it
      ! sees param(12) < 0, after which dt = abs(param(12)) on every step.
            param(12) = -abs(bufdt(j))

            call neklab_timer_stop(t_bf_set)
         end subroutine bf_set

         subroutine bf_set_window(k)
      !! Direction-agnostic alternative to bf_set: installs snapshot k AND
      !! builds the lag levels explicitly from slots k-1 and k-2, instead of
      !! relying on the traversal order to have left the right thing in vx.
      !!
      !! This is what the adjoint (backward) replay needs, since shifting the
      !! current field down a level would put buf(k+1) where buf(k-1) belongs.
      !!
      !! LIMITATION: slots k-1 and k-2 must live in the same resident chunk,
      !! because at most one chunk is ever in memory. At a chunk boundary this
      !! fails rather than silently using stale lag levels. The cheap fix, when
      !! the backward path goes live, is a two-snapshot side cache carried
      !! across the boundary (~2/lbuf overhead) rather than a second chunk.
            integer, intent(in) :: k
      ! internal
            character(len=*), parameter :: this_procedure = 'bf_set_window'
            character(len=256) :: msg
            integer :: c, j
            call check_init(this_procedure)
            if (k < 1 .or. k > nsteps_rec) then
               write (msg, '(A,I0,A,I0,A)') 'Requested step ', k, ' outside [1, ', nsteps_rec, '].'
               call nek_stop_error(msg, this_module, this_procedure)
            end if
            call neklab_timer_start(t_bf_window)

            c = (k - 1)/lbuf + 1
            j = k - (c - 1)*lbuf
            if (c /= ichunk) call load_chunk(c)

            if ((k >= 2 .and. j < 2) .or. (k >= 3 .and. j < 3)) then
               write (msg, '(A,I0,A,I0,A)') 'Step ', k, ' needs lag slots from chunk ', c - 1,
     &            '. Cross-chunk lag access is not implemented (see bf_set_window).'
               call nek_stop_error(msg, this_module, this_procedure)
            end if

            call copy(vx, bufx(1, j), nbf)
            call copy(vy, bufy(1, j), nbf)
            call copy(t(1, 1, 1, 1, 1), bufs(1, j), nbf)

            if (k >= 2) then
               call copy(vxlag(1, 1, 1, 1, 1), bufx(1, j - 1), nbf)
               call copy(vylag(1, 1, 1, 1, 1), bufy(1, j - 1), nbf)
               call copy(tlag(1, 1, 1, 1, 1, 1), bufs(1, j - 1), nbf)
            end if
            if (k >= 3) then
               call copy(vxlag(1, 1, 1, 1, 2), bufx(1, j - 2), nbf)
               call copy(vylag(1, 1, 1, 1, 2), bufy(1, j - 2), nbf)
               call copy(tlag(1, 1, 1, 1, 1, 2), bufs(1, j - 2), nbf)
            end if

            param(12) = -abs(bufdt(j))

            call neklab_timer_stop(t_bf_window)
         end subroutine bf_set_window

      !====================================================================
      !     ACCESSORS
      !====================================================================

         pure function bf_get_nsteps() result(n)
            integer :: n
            n = nsteps_rec
         end function bf_get_nsteps

         pure function bf_get_nsteps_min() result(n)
            integer :: n
            n = nsteps_min
         end function bf_get_nsteps_min

         pure function bf_get_nchunk() result(n)
            integer :: n
            n = nchunk
         end function bf_get_nchunk

         pure function bf_get_time() result(tt)
      !! Total time spanned by the recorded trajectory.
            real(dp) :: tt
            tt = trec
         end function bf_get_time

         pure function bf_is_recording() result(rec)
            logical :: rec
            rec = is_recording
         end function bf_is_recording

         pure function bf_is_initialised() result(ini)
            logical :: ini
            ini = is_initialized
         end function bf_is_initialised

         subroutine bf_get_dt_minmax(dtmm)
      !! Extremes of the physical (non-landing) timesteps of the last recording.
      !! Mirrors helix's get_dt_minmax; the flow-rate driver uses the mean to
      !! estimate the quadrature error of the Fourier accumulator.
            real(dp), dimension(2), intent(out) :: dtmm
            dtmm(1) = dt_min
            dtmm(2) = dt_max
            if (dt_max == 0.0_dp) dtmm = 0.0_dp
         end subroutine bf_get_dt_minmax

         function bf_get_dt(k) result(dtk)
      !! Timestep consumed by step k. Loads the containing chunk if needed, so
      !! it must not be called while a recording is in progress.
            integer, intent(in) :: k
            real(dp) :: dtk
      ! internal
            integer :: c, j
            c = (k - 1)/lbuf + 1
            j = k - (c - 1)*lbuf
            if (c /= ichunk) call load_chunk(c)
            dtk = bufdt(j)
         end function bf_get_dt

         subroutine bf_summary()
            character(len=*), parameter :: this_procedure = 'bf_summary'
            character(len=256) :: msg
            write (msg, '(A,I0,A,I0,A,I0,A,E16.8)') 'buffer: nsteps= ', nsteps_rec,
     &         ', chunks= ', nchunk, ', resident= ', ichunk, ', T= ', trec
            call nek_log_message(msg, this_module, this_procedure)
         end subroutine bf_summary

      !====================================================================
      !     MESH DUMP AND RESTART METADATA
      !====================================================================

         subroutine bf_write_mesh(prefix)
      !! Writes a coordinates-only file. Needed as the reference full-mesh file
      !! for the Python half -> full reconstruction used by the stability runs.
            character(len=1), optional, intent(in) :: prefix
      ! internal
            character(len=*), parameter :: this_procedure = 'bf_write_mesh'
            character(len=132) :: fname
            call check_init(this_procedure)
            fname = bf_fname(0, prefix)
            call bf_write_chunk(fname, only_mesh=.true.)
            call nek_log_message('Mesh written to '//trim(fname), this_module, this_procedure)
         end subroutine bf_write_mesh

         subroutine bf_write_t2Dh(dpds, nf, omega, period)
      !! Sidecar text file holding everything a restart needs that is not in the
      !! field files. Kept out of the binary so the format stays byte-identical
      !! to the helix files and the existing Python readers keep working.
      !! Not currently consumed by the driver: partial-run restart is not wired.
            real(dp), dimension(:), intent(in) :: dpds
            integer, intent(in) :: nf
            real(dp), intent(in) :: omega
            real(dp), intent(in) :: period
      ! internal
            character(len=*), parameter :: this_procedure = 'bf_write_t2Dh'
            character(len=132) :: fname
            integer :: iunit, i
            if (nid /= 0) return
            fname = ' '
            write (fname, '(A,A,A)') fprefix, trim(fbase), '.t2Dh'
            open (newunit=iunit, file=trim(fname), status='replace', action='write')
            write (iunit, '(A)') '# neklab baseflow buffer t2Dh state'
            write (iunit, '(A,I0)') 'nf     = ', nf
            write (iunit, '(A,I0)') 'nsteps = ', nsteps_rec
            write (iunit, '(A,I0)') 'nchunk = ', nchunk
            write (iunit, '(A,I0)') 'lbuf   = ', lbuf
            write (iunit, '(A,E24.16)') 'omega  = ', omega
            write (iunit, '(A,E24.16)') 'period = ', period
            write (iunit, '(A,*(1X,E24.16))') 'dpds   =', (dpds(i), i=1, nf)
            close (iunit)
            call nek_log_information('Control state written to '//trim(fname), this_module, this_procedure)
         end subroutine bf_write_t2Dh

         subroutine bf_read_t2Dh(dpds, nf, omega, period, ierr)
      !! Reads the sidecar file back. ierr /= 0 if it is missing.
            real(dp), dimension(:), intent(out) :: dpds
            integer, intent(out) :: nf
            real(dp), intent(out) :: omega
            real(dp), intent(out) :: period
            integer, intent(out) :: ierr
      ! internal
            character(len=*), parameter :: this_procedure = 'bf_read_t2Dh'
            character(len=132) :: fname
            character(len=1024) :: line
            integer :: iunit, i, ipos
            logical :: exists
            dpds = 0.0_dp; nf = 0; omega = 0.0_dp; period = 0.0_dp; ierr = 0
            fname = ' '
            write (fname, '(A,A,A)') fprefix, trim(fbase), '.t2Dh'
            if (nid == 0) then
               inquire (file=trim(fname), exist=exists)
               if (.not. exists) then
                  ierr = 1
               else
      ! Two passes: nf must be known before the dpds line can be parsed, and
      ! the file is small enough that rewinding is cheaper than buffering.
                  open (newunit=iunit, file=trim(fname), status='old', action='read')
                  do
                     read (iunit, '(A)', end=100) line
                     ipos = index(line, '=')
                     if (ipos == 0 .or. line(1:1) == '#') cycle
                     if (index(line, 'nf') > 0 .and. index(line, 'nchunk') == 0) read (line(ipos + 1:), *) nf
                     if (index(line, 'nsteps') > 0) read (line(ipos + 1:), *) nsteps_rec
                     if (index(line, 'nchunk') > 0) read (line(ipos + 1:), *) nchunk
                     if (index(line, 'omega') > 0) read (line(ipos + 1:), *) omega
                     if (index(line, 'period') > 0) read (line(ipos + 1:), *) period
                  end do
100               rewind (iunit)
                  if (nf > size(dpds)) ierr = 2
                  if (ierr == 0) then
                     do
                        read (iunit, '(A)', end=200) line
                        ipos = index(line, '=')
                        if (ipos == 0 .or. index(line, 'dpds') == 0) cycle
                        read (line(ipos + 1:), *) (dpds(i), i=1, nf)
                     end do
                  end if
200               close (iunit)
               end if
            end if
            call bcast(ierr, isize)
            if (ierr /= 0) then
               call nek_log_warning('Could not read t2Dh state from '//trim(fname), this_module, this_procedure)
               return
            end if
            call bcast(nf, isize)
            call bcast(nsteps_rec, isize)
            call bcast(nchunk, isize)
            call bcast(omega, wdsize)
            call bcast(period, wdsize)
            call bcast(dpds, size(dpds)*wdsize)
            ichunk = 0
            trec = period
            call nek_log_message('t2Dh state read from '//trim(fname), this_module, this_procedure)
         end subroutine bf_read_t2Dh

      !====================================================================
      !     INTERNAL HELPERS
      !====================================================================

         subroutine check_init(caller)
            character(len=*), intent(in) :: caller
            if (.not. is_initialized) then
               call nek_stop_error('Buffer not initialized. The driver must call bf_init '//
     &            'before the first residual evaluation.', this_module, caller)
            end if
         end subroutine check_init

         subroutine flush_chunk(keep_resident)
      !! Writes the resident chunk. keep_resident = .true. leaves the snapshots
      !! in memory (used at the end of a recording so that a trajectory fitting
      !! in one chunk replays with no IO at all).
            logical, intent(in) :: keep_resident
      ! internal
            character(len=*), parameter :: this_procedure = 'flush_chunk'
            character(len=132) :: fname
            character(len=256) :: msg
            integer :: k0, k1
            nchunk = nchunk + 1
            k0 = nsteps_rec - nsave + 1
            k1 = nsteps_rec
            if (if_write) then
               fname = bf_fname(nchunk)
               call bf_write_chunk(fname)
      ! One terse line per file actually written, with the step range, so a
      ! recording can be cross-checked at a glance.
               write (msg, '(A,A,I5,A)') trim(fname), ': ', nsave, ' snapshots'
               call nek_log_information(msg, this_module, this_procedure)
            else
               write (msg, '(A,I0,A,I0,A,I0,A)') 'Chunk ', nchunk, ' held in memory (', nsave,
     &            ' snapshots, steps ', k0, '-', k1, '); spilling is off.'
               call nek_log_debug(msg, this_module, this_procedure)
               if (nchunk > 1) then
                  call nek_stop_error('Trajectory exceeds lbuf but chunk spilling is off. '//
     &               'Increase lbuf or enable write_chunks in bf_init.', this_module, this_procedure)
               end if
            end if
            if (keep_resident) then
               ichunk = nchunk
            else
      ! The buffer now starts filling chunk nchunk+1, so nothing is resident.
      ! Claiming otherwise would let load_chunk short-circuit and hand back a
      ! half-filled buffer.
               ichunk = 0
               nsave = 0
            end if
         end subroutine flush_chunk

         subroutine load_chunk(c)
      !! Makes chunk c resident. At most one chunk is ever in memory.
            integer, intent(in) :: c
      ! internal
            character(len=*), parameter :: this_procedure = 'load_chunk'
            character(len=132) :: fname
            character(len=256) :: msg
            integer :: nread, k0, k1
            if (c == ichunk) return
            if (is_recording) then
               call nek_stop_error('Cannot load a chunk while recording: it would '//
     &            'overwrite the snapshots being collected.', this_module, this_procedure)
            end if
            if (.not. if_write) then
               write (msg, '(A,I0,A)') 'Chunk ', c, ' is not resident and spilling is off.'
               call nek_stop_error(msg, this_module, this_procedure)
            end if
            fname = bf_fname(c)
            call bf_read_chunk(fname, nread)
            nsave = nread
            ichunk = c
            k0 = (c - 1)*lbuf + 1
            k1 = k0 + nread - 1
            write (msg, '(A,A,I0,A,I0,A,I0)') trim(fname), ': ', nread, ' snapshots, steps ', k0, '-', k1
      ! Reads happen once per chunk per matvec, i.e. O(nchunk * nmatvec) times
      ! per Newton step. Writes are the ones worth seeing by default.
            call nek_log_debug(msg, this_module, this_procedure)
         end subroutine load_chunk

      end module t2Dh_bf_buffer