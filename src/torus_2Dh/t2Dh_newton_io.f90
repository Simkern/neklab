      module t2Dh_newton_io
      !! Restart record for the flowrate Newton iteration.
      !!
      !! ONE file per run, plain text, written by nid 0:
      !!
      !!   <tag>_final.nwt   the state the iteration exited on. Read it back in
      !!                     userchk with nwt_restore.
      !!
      !! It is written unconditionally, converged or not: a run that produced a
      !! converged orbit and no record of how to reproduce it is a wasted run,
      !! and a run that diverged wants its record more than a converged one
      !! does.
      !!
      !! WHAT A RESTART ACTUALLY NEEDS is the forcing vector and, optionally,
      !! the steady resistance that seeds the outer Jacobian. Everything else in
      !! the file is there to be read by a person, or to be checked against the
      !! run doing the reading. The record deliberately does NOT carry the inner
      !! Newton-Krylov residual history, the step sizes or the timings: those
      !! belong in the log, they differ every run, and nothing consumes them.
      !!
      !! NUMBER FORMAT. Anything between 1e-4 and 1e7 is written in fixed point
      !! with trailing zeros trimmed, so a Womersley number of 40 reads as 40.0
      !! rather than 0.4000000000000000E+02. Outside that window it falls back
      !! to exponential, because a fixed-point rendering of 7e-12 is a screen of
      !! zeros. By default enough digits are kept to round-trip a double
      !! exactly; the quantities that are only ever read by eye or compared to a
      !! tolerance -- the geometry, Re, Wo, the targets -- are printed to
      !! nsig_short digits instead, so a target of 0.02 reads as 0.02 and not as
      !! 0.020000000000000004.
      !!
      !! Note that this is a PRINTING choice, not a rounding of the stored
      !! values: nothing in the record is mutated on its way out. dpds, which is
      !! what a restart actually consumes, always goes out at full precision.
      !!
      !! WHY NOT EXTEND bf_write_t2Dh: that sidecar is parsed by substring match
      !! (index(line, 'nf') > 0), so any key merely CONTAINING 'nf' would
      !! silently corrupt it. It is left untouched. The format below splits on
      !! the first '=', normalises the key at the first blank or '(' -- which is
      !! what lets a key be labelled 'dpds  (native)' for the reader's benefit
      !! without the parser caring -- and skips keys it does not know, so fields
      !! can be added later without breaking a binary built before them.
         use stdlib_optval, only: optval
         use stdlib_strings, only: padr
         use LightKrylov, only: dp
         use neklab_nek_setup, only: nek_log_message, nek_log_information,
     &                               nek_log_warning, nek_stop_error
         use neklab_t2Dh, only: t2Dh, lfc, lmfc
         implicit none
         include "SIZE"
         include "TOTAL"
         private
         character(len=*), parameter, private :: this_module = 't2Dh_newton_io'

         real(dp), parameter, private :: check_tol = 1.0e-06_dp
      !! Relative agreement demanded of the geometry and the Reynolds number
      !! when a record is restored. Loose enough to survive a recompile, tight
      !! enough to catch a different mesh or a different .par file.

         integer, parameter, private :: nsig_short = 8
      !! Significant digits for the quantities printed for the eye. Must stay
      !! well inside check_tol: printing the geometry to 8 digits perturbs it by
      !! 1e-8 relative, which the 1e-6 check absorbs without noticing.

      !--------------------------------------------------------------------
      !-----     THE RECORD                                           -----
      !--------------------------------------------------------------------

         type, public :: t2Dh_nwt_record
      !! Fixed-size arrays throughout (lfc, lmfc are compile-time), so the type
      !! is trivially broadcastable and a record can be copied by assignment.
      ! --- sizes
            integer :: nf = 0
            integer :: nmf = 0
            integer :: kharm = 0
      ! --- geometry. Checked on restore, never restored.
            real(dp) :: radius = 0.0_dp
            real(dp) :: curv_radius = 0.0_dp
            real(dp) :: delta = 0.0_dp
      ! --- physics
            real(dp) :: Re = 0.0_dp        ! 1/nu, from cpfld(1,1). Checked on restore.
            real(dp) :: Wo = 0.0_dp
            real(dp) :: omega = 0.0_dp
            real(dp) :: period = 0.0_dp
      ! --- forcing. dpds IS the restart.
            real(dp), dimension(lfc) :: dpds = 0.0_dp
      !! Native convention. What init_flow consumes.
            real(dp), dimension(lfc) :: dpds_helix = 0.0_dp
      !! Legacy helix convention, for comparison with a helix deck. Never read.
            real(dp), dimension(lmfc) :: amp = 0.0_dp
            real(dp), dimension(lmfc) :: phase = 0.0_dp
      ! --- flow rate
            real(dp), dimension(lmfc) :: mflow = 0.0_dp
            real(dp), dimension(lmfc) :: mf_target = 0.0_dp
            real(dp), dimension(lmfc) :: mf_err = 0.0_dp
            real(dp), dimension(lmfc) :: qerr = 0.0_dp
            real(dp) :: enorm = 0.0_dp
            real(dp) :: eabs = 0.0_dp
            real(dp) :: noise = 0.0_dp
      ! --- outer Jacobian. Only the (1,1) entry is restored, as the seed.
            real(dp), dimension(lmfc, lmfc) :: jac = 0.0_dp
      ! --- discretisation, so the record identifies the run it came from
            integer :: lx1_run = 0
            integer :: nelg_run = 0
            integer :: np_run = 0
      ! --- outcome
            logical :: converged = .false.
         contains
            procedure, pass(self), public :: from_t2Dh => nwt_from_t2Dh
         end type t2Dh_nwt_record

      !--------------------------------------------------------------------
      !-----     RUN STATE                                            -----
      !--------------------------------------------------------------------
      ! Captured once at nwt_open, stamped into the record as comments.

         character(len=8), private :: tag = 'nwt'
         logical, private :: opened = .false.

         real(dp), private :: s_tol = 0.0_dp
         real(dp), private :: s_rtol_mf = 0.0_dp
         integer, private :: s_tol_mode = 1
         integer, private :: s_maxiter = 0
         integer, private :: s_maxiter_inner = 0
         logical, private :: s_inexact = .true.
         character(len=8), private :: s_jac0 = 'fd'
         real(dp), dimension(lmfc), private :: s_scal = 1.0_dp

         public :: nwt_open, nwt_close
         public :: nwt_read, nwt_restore, nwt_fname, nwt_is_open

      contains

      !====================================================================
      !     SETUP
      !====================================================================

         subroutine nwt_open(run_tag, tol, rtol_mf, tol_mode, maxiter,
     &                       maxiter_inner, if_inexact, jac0, scal)
      !! Captures the solver settings and the file stem. Touches no file: the
      !! record is written in one shot at nwt_close.
            character(len=*), intent(in) :: run_tag
      !! File stem. The outpost prefix of the run is the natural choice, so the
      !! steady solve writes nwq_final.nwt and the pulsatile one nwf_final.nwt
      !! and two solves in a session do not overwrite each other.
            real(dp), intent(in) :: tol
            real(dp), intent(in) :: rtol_mf
            integer, intent(in) :: tol_mode
            integer, intent(in) :: maxiter
            integer, intent(in) :: maxiter_inner
            logical, intent(in) :: if_inexact
            character(len=*), intent(in) :: jac0
            real(dp), dimension(:), intent(in) :: scal
      !! Norm scaling. Without it enorm cannot be reconstructed from mf_err.
      ! internal
            integer :: n

            tag = run_tag
            s_tol = tol
            s_rtol_mf = rtol_mf
            s_tol_mode = tol_mode
            s_maxiter = maxiter
            s_maxiter_inner = maxiter_inner
            s_inexact = if_inexact
            s_jac0 = jac0
            s_scal = 1.0_dp
            n = min(size(scal), lmfc)
            s_scal(1:n) = scal(1:n)
            opened = .true.
         end subroutine nwt_open

         pure function nwt_is_open() result(l)
            logical :: l
            l = opened
         end function nwt_is_open

         pure function nwt_fname() result(fname)
      !! The one filename this module produces.
            character(len=132) :: fname
            fname = ' '
            write (fname, '(A,A)') trim(tag), '_final.nwt'
         end function nwt_fname

      !====================================================================
      !     FILLING A RECORD
      !====================================================================

         subroutine nwt_from_t2Dh(self)
      !! Pulls everything the control object already knows. The caller fills the
      !! flow rate, the errors and the Jacobian, which live in the driver.
            class(t2Dh_nwt_record), intent(inout) :: self
      ! internal
            real(dp), dimension(lfc) :: d
            real(dp), dimension(lmfc) :: a, p

            self%nf = t2Dh%get_nf()
            self%nmf = t2Dh%get_nmf()
            self%kharm = t2Dh%get_kharm()

            self%radius = t2Dh%get_radius()
            self%curv_radius = t2Dh%get_curv_radius()
            self%delta = t2Dh%get_delta()

            self%Re = 1.0_dp/cpfld(1,1)
            self%Wo = t2Dh%get_womersley()
            self%omega = t2Dh%get_omega()
            self%period = t2Dh%get_period()

            call t2Dh%get_dpds(d, amp=a, phase=p)
            self%dpds = d
            self%amp = a
            self%phase = p
            self%dpds_helix = 0.0_dp
            call t2Dh%get_dpds_helix(self%dpds_helix)

            self%mf_target = t2Dh%get_target()

            self%lx1_run = lx1
            self%nelg_run = nelgv
            self%np_run = np
         end subroutine nwt_from_t2Dh

      !====================================================================
      !     WRITE
      !====================================================================

         subroutine nwt_close(rec, converged)
      !! Writes the record. UNCONDITIONAL: it does not check `opened` and does
      !! not care whether the iteration converged. If nwt_open was never called
      !! the settings block says so rather than printing its defaults as fact.
            type(t2Dh_nwt_record), intent(inout) :: rec
            logical, intent(in) :: converged
      ! internal
            character(len=*), parameter :: this_procedure = 'nwt_close'
            character(len=132) :: fname
            integer :: iunit

            rec%converged = converged
            fname = nwt_fname()

            if (nid == 0) then
               open (newunit=iunit, file=trim(fname), status='replace', action='write')
               write (iunit, '(A)') '# neklab flowrate Newton final state'
               write (iunit, '(A)') '# Read back with nwt_restore in userchk.'
               call write_settings(iunit, '# ')
               call write_record(iunit, rec)
               close (iunit)
            end if

            call nek_log_message('Final Newton record -> '//trim(fname),
     &         this_module, this_procedure)
            opened = .false.
         end subroutine nwt_close

      !--------------------------------------------------------------------
      !-----     WRITE HELPERS (nid 0 only, unit already open)        -----
      !--------------------------------------------------------------------

         subroutine write_settings(iunit, pfx)
      !! Provenance, not state. Nothing reads it back; anything a restart NEEDS
      !! is in write_record.
            integer, intent(in) :: iunit
            character(len=*), intent(in) :: pfx
      ! internal
            integer :: i
            if (.not. opened) then
               write (iunit, '(A,A)') pfx, 'solver settings NOT CAPTURED (nwt_open was never called)'
               return
            end if
            write (iunit, '(A,A)') pfx, 'solver settings'
            write (iunit, '(A,A,A)') pfx, '  inner target tol   = ', trim(adjustl(rstr(s_tol)))
            write (iunit, '(A,A,A)') pfx, '  outer rel. tol     = ', trim(adjustl(rstr(s_rtol_mf)))
            write (iunit, '(A,A,A)') pfx, '  tol. scheduling    = ',
     &         merge('constant', ' dynamic', s_tol_mode == 1)
            write (iunit, '(A,A,I0)') pfx, '  max outer steps    = ', s_maxiter
            write (iunit, '(A,A,I0)') pfx, '  max inner iters    = ', s_maxiter_inner
            write (iunit, '(A,A,L1)') pfx, '  inexact outer      = ', s_inexact
            write (iunit, '(A,A,A)') pfx, '  initial jacobian   = ', trim(s_jac0)
            write (iunit, '(A,A,*(2X,A))') pfx, '  norm scale         =',
     &         (trim(adjustl(rstr(s_scal(i), nsig_short))), i=1, lmfc)
         end subroutine write_settings

         subroutine write_record(iunit, rec)
      !! The state block: key = value, one key per line, first '=' separates,
      !! key normalised at the first blank or '('.
      !!
      !! nsig_short is passed for the quantities that are read by eye or only
      !! ever compared to a tolerance. dpds and the residual quantities are not
      !! rounded: the first is the restart, the second is what tells you whether
      !! to trust it.
            integer, intent(in) :: iunit
            type(t2Dh_nwt_record), intent(in) :: rec
      ! internal
            character(len=24) :: key
            integer :: i, n, m

            n = max(rec%nf, 1)
            m = max(rec%nmf, 1)

            write (iunit, '(A)') '#'
            write (iunit, '(A)') '# --- identification'
            call put_i(iunit, 'nf', rec%nf)
            call put_i(iunit, 'nmf', rec%nmf)
            call put_i(iunit, 'kharm', rec%kharm)
            write (iunit, '(A)') '# --- geometry'
            call put_r(iunit, 'radius', rec%radius, nsig_short)
            call put_r(iunit, 'curv_radius', rec%curv_radius, nsig_short)
            call put_r(iunit, 'delta', rec%delta, nsig_short)
            write (iunit, '(A)') '# --- physics'
            call put_r(iunit, 'Re', rec%Re, nsig_short)
            call put_r(iunit, 'Wo', rec%Wo, nsig_short)
            call put_r(iunit, 'omega', rec%omega)
            call put_r(iunit, 'period', rec%period)
            write (iunit, '(A)') '# --- forcing'
            call put_v(iunit, 'dpds  (native)', rec%dpds, n)
            call put_v(iunit, 'dpds_helix  (legacy)', rec%dpds_helix, n)
            call put_v(iunit, 'amplitude', rec%amp, m)
            call put_v(iunit, 'phase', rec%phase, m)
            write (iunit, '(A)') '# --- flow rate'
            call put_v(iunit, 'mflow', rec%mflow, m)
            call put_v(iunit, 'mf_target', rec%mf_target, m, nsig_short)
            call put_v(iunit, 'mf_err', rec%mf_err, m)
            call put_v(iunit, 'qerr', rec%qerr, m)
            call put_r(iunit, 'enorm', rec%enorm)
            call put_r(iunit, 'eabs', rec%eabs)
            call put_r(iunit, 'noise', rec%noise)
            write (iunit, '(A)') '# --- outer jacobian d(mf_i)/d(a_j)'
            do i = 1, m
               key = ' '
               write (key, '(A,I0)') 'jac_row', i
               call put_v(iunit, trim(key), rec%jac(i, :), m)
            end do
            write (iunit, '(A)') '# --- discretisation'
            call put_i(iunit, 'lx1', rec%lx1_run)
            call put_i(iunit, 'nelg', rec%nelg_run)
            call put_i(iunit, 'np', rec%np_run)
            write (iunit, '(A)') '# --- outcome'
            write (iunit, '(A,A,L1)') padr('converged', 22), '= ', rec%converged
         end subroutine write_record

         subroutine put_i(iunit, key, ival)
            integer, intent(in) :: iunit
            character(len=*), intent(in) :: key
            integer, intent(in) :: ival
            write (iunit, '(A,A,I0)') padr(key, 22), '= ', ival
         end subroutine put_i

         subroutine put_r(iunit, key, v, nsig)
            integer, intent(in) :: iunit
            character(len=*), intent(in) :: key
            real(dp), intent(in) :: v
            integer, optional, intent(in) :: nsig
            write (iunit, '(A,A,A)') padr(key, 22), '= ', trim(adjustl(rstr(v, nsig)))
         end subroutine put_r

         subroutine put_v(iunit, key, v, n, nsig)
      !! Vector values are right-justified in a fixed field so the components of
      !! successive lines sit under each other and can be compared by eye.
            integer, intent(in) :: iunit
            character(len=*), intent(in) :: key
            real(dp), dimension(:), intent(in) :: v
            integer, intent(in) :: n
            integer, optional, intent(in) :: nsig
            integer :: i
            write (iunit, '(A,A,*(2X,A))') padr(key, 22), '=', (rstr(v(i), nsig), i=1, n)
         end subroutine put_v

         function rstr(v, nsig) result(s)
      !! Readable rendering of a double, right-justified in 24 characters.
      !!
      !! Fixed point with trailing zeros trimmed inside [1e-4, 1e7), which is
      !! where the quantities a person actually reads live -- radii, Womersley
      !! numbers, flow rates, forcing amplitudes. Exponential outside it, where
      !! fixed point would be a screen of zeros.
      !!
      !! nsig is the number of SIGNIFICANT digits, defaulting to the 17 a double
      !! needs to round-trip exactly. The decimal count is derived from the
      !! exponent, so the digits requested are the digits delivered whichever
      !! branch is taken. Rounding here rather than in the record keeps the
      !! stored values untouched -- and avoids the trap of rounding a value to
      !! 8 digits and then printing it to 17, which renders 3.3333323 as
      !! 3.3333322999999999.
            real(dp), intent(in) :: v
            integer, optional, intent(in) :: nsig
            character(len=24) :: s
      ! internal
            character(len=40) :: w
            character(len=16) :: fmtstr
            real(dp) :: a
            integer :: nd, e, ip, ns

            ns = 17
            if (present(nsig)) ns = max(min(nsig, 17), 1)

            a = abs(v)
            if (a == 0.0_dp) then
               w = '0.0'
            else if (a >= 1.0e-4_dp .and. a < 1.0e7_dp) then
               e = floor(log10(a))
               nd = min(max(ns - 1 - e, 1), 24)
               write (fmtstr, '(A,I0,A)') '(F0.', nd, ')'
               write (w, fmtstr) v
               ip = len_trim(w)
               do while (ip > 1 .and. w(ip:ip) == '0')
                  ip = ip - 1
               end do
      ! Never leave a bare trailing '.': '40.' is legal Fortran input but reads
      ! as a typo to everyone else.
               if (w(ip:ip) == '.') then
                  w = w(1:ip)//'0'
               else
                  w = w(1:ip)
               end if
      ! F0.d drops the leading zero of a number below one: '.30000009'. Legal
      ! input, but it reads as a typo and misaligns against its neighbours.
               if (w(1:1) == '.') then
                  w = '0'//trim(w)
               else if (w(1:2) == '-.') then
                  w = '-0'//trim(w(2:))
               end if
            else
               write (fmtstr, '(A,I0,A,I0,A)') '(E', ns + 8, '.', ns - 1, ')'
               write (w, fmtstr) v
               w = adjustl(w)
            end if
      ! Belt and braces: if the fixed-point form somehow overran the field, fall
      ! back rather than truncate a number that has to be read back.
            if (len_trim(w) > len(s)) then
               write (w, '(E24.16)') v
               w = adjustl(w)
            end if
            s = ' '
            ip = len_trim(w)
            s(len(s) - ip + 1:) = w(1:ip)
         end function rstr

      !====================================================================
      !     READ
      !====================================================================

         subroutine nwt_read(fname, rec, ierr)
      !! Reads a record file written by nwt_close.
      !!
      !! nid 0 parses and broadcasts. Unknown keys are skipped silently, which
      !! is what lets a file written by a newer build be read by an older one.
      !! ierr: 0 fine, 1 file missing, 2 no recognisable content.
            character(len=*), intent(in) :: fname
            type(t2Dh_nwt_record), intent(out) :: rec
            integer, intent(out) :: ierr
      ! internal
            character(len=*), parameter :: this_procedure = 'nwt_read'
            character(len=1024) :: line, val
            character(len=32) :: key
            character(len=256) :: msg
            integer :: iunit, i, n, m, irow, nkey
            logical :: exists, ok

            ierr = 0
            nkey = 0

            if (nid == 0) then
               inquire (file=trim(fname), exist=exists)
               if (.not. exists) then
                  ierr = 1
               else
                  open (newunit=iunit, file=trim(fname), status='old', action='read')
      ! Single pass: nf and nmf are written before every array they size.
                  do
                     read (iunit, '(A)', end=100) line
                     call split_kv(line, key, val, ok)
                     if (.not. ok) cycle
                     nkey = nkey + 1
                     n = max(rec%nf, 1)
                     m = max(rec%nmf, 1)
                     select case (trim(key))
                     case ('nf'); read (val, *) rec%nf
                     case ('nmf'); read (val, *) rec%nmf
                     case ('kharm'); read (val, *) rec%kharm
                     case ('radius'); read (val, *) rec%radius
                     case ('curv_radius'); read (val, *) rec%curv_radius
                     case ('delta'); read (val, *) rec%delta
                     case ('Re'); read (val, *) rec%Re
                     case ('Wo'); read (val, *) rec%Wo
                     case ('omega'); read (val, *) rec%omega
                     case ('period'); read (val, *) rec%period
                     case ('dpds'); read (val, *) (rec%dpds(i), i=1, n)
                     case ('dpds_helix'); read (val, *) (rec%dpds_helix(i), i=1, n)
                     case ('amplitude'); read (val, *) (rec%amp(i), i=1, m)
                     case ('phase'); read (val, *) (rec%phase(i), i=1, m)
                     case ('mflow'); read (val, *) (rec%mflow(i), i=1, m)
                     case ('mf_target'); read (val, *) (rec%mf_target(i), i=1, m)
                     case ('mf_err'); read (val, *) (rec%mf_err(i), i=1, m)
                     case ('qerr'); read (val, *) (rec%qerr(i), i=1, m)
                     case ('enorm'); read (val, *) rec%enorm
                     case ('eabs'); read (val, *) rec%eabs
                     case ('noise'); read (val, *) rec%noise
                     case ('lx1'); read (val, *) rec%lx1_run
                     case ('nelg'); read (val, *) rec%nelg_run
                     case ('np'); read (val, *) rec%np_run
                     case ('converged'); rec%converged = (val(1:1) == 'T')
                     case default
                        if (key(1:7) == 'jac_row') then
                           read (key(8:), *) irow
                           if (irow >= 1 .and. irow <= lmfc) then
                              read (val, *) (rec%jac(irow, i), i=1, m)
                           end if
                        else
                           nkey = nkey - 1
                        end if
                     end select
                  end do
100               close (iunit)
                  if (nkey == 0 .or. rec%nf == 0) ierr = 2
               end if
            end if

            call bcast(ierr, isize)
            if (ierr /= 0) then
               write (msg, '(A,A,A,I0,A)') 'Could not read ', trim(fname), ' (ierr= ', ierr, ')'
               call nek_log_warning(msg, this_module, this_procedure)
               return
            end if
            call bcast_record(rec)
         end subroutine nwt_read

         subroutine nwt_restore(fname, ierr, if_womersley, womersley,
     &                          if_seed_jacobian, if_check, rec_out)
      !! Reads a record and configures t2Dh from it: the one call a userchk
      !! needs to pick a run up where the last one stopped.
      !!
      !! It restores the FORCING, and optionally the steady resistance that
      !! seeds the outer Jacobian. It does NOT restore the flowrate targets --
      !! the deck says what it is aiming for -- and it does not restore the
      !! state: the velocity field still comes from the field file the way it
      !! always does.
      !!
      !! Geometry and Reynolds number are CHECKED, not restored. A forcing
      !! converged on one mesh at one Re is not a good guess on another, and
      !! finding that out from a diverging Newton three hours in is worse than
      !! finding it out here.
            character(len=*), intent(in) :: fname
            integer, intent(out) :: ierr
            logical, optional, intent(in) :: if_womersley
      !! Take the Womersley number from the file. Default .true. Set .false. to
      !! continue this forcing at a DIFFERENT frequency, in which case the
      !! womersley argument must be present: a restart that quietly picked up a
      !! Wo from somewhere else would change the period without saying so.
            real(dp), optional, intent(in) :: womersley
      !! Womersley number to use when if_womersley is .false.
            logical, optional, intent(in) :: if_seed_jacobian
      !! Hand the stored (1,1) Jacobian entry over as the steady resistance,
      !! which saves the next run its finite-difference probe of the same
      !! quantity. Default .true.
            logical, optional, intent(in) :: if_check
      !! Enforce the geometry and Reynolds checks. Default .true. Turn it off
      !! only to deliberately carry a forcing onto a different mesh.
            type(t2Dh_nwt_record), optional, intent(out) :: rec_out
      ! internal
            character(len=*), parameter :: this_procedure = 'nwt_restore'
            character(len=256) :: msg
            type(t2Dh_nwt_record) :: rec
            real(dp) :: Wo, Re_run
            logical :: use_file_Wo, seed_, check_

            use_file_Wo = optval(if_womersley, .true.)
            seed_ = optval(if_seed_jacobian, .true.)
            check_ = optval(if_check, .true.)

            call nwt_read(fname, rec, ierr)
            if (ierr /= 0) return

      ! ---- Womersley number
            if (use_file_Wo) then
               Wo = rec%Wo
               if (present(womersley)) then
                  call nek_log_warning('womersley supplied but if_womersley is .true.: '//
     &               'the file value is used. Pass if_womersley = .false. to override.',
     &               this_module, this_procedure)
               end if
            else
               if (.not. present(womersley)) then
                  call nek_stop_error('if_womersley = .false. requires the womersley argument.',
     &               this_module, this_procedure)
               end if
               Wo = womersley
               write (msg, '(A,A,A,A)') 'Womersley overridden: file ',
     &            trim(adjustl(rstr(rec%Wo, nsig_short))), ' -> ',
     &            trim(adjustl(rstr(Wo, nsig_short)))
               call nek_log_warning(msg, this_module, this_procedure)
            end if

      ! ---- geometry and Reynolds number, measured from THIS run's mesh and
      !      then compared with the record. Neither is restored. All four checks
      !      fire before any forcing is installed, so a mismatched restart
      !      leaves t2Dh untouched rather than half-configured.
            Re_run = 1.0_dp/cpfld(1,1)
            if (check_) then
               call check_scalar('Re', rec%Re, Re_run)
               call check_scalar('radius', rec%radius, t2Dh%get_radius())
               call check_scalar('curv_radius', rec%curv_radius, t2Dh%get_curv_radius())
               call check_scalar('delta', rec%delta, t2Dh%get_delta())
            else if (.not. agree(rec%Re, Re_run)) then
               write (msg, '(A,A,A,A)') 'Re differs from the record (file ',
     &            trim(adjustl(rstr(rec%Re, nsig_short))), ', run ',
     &            trim(adjustl(rstr(Re_run, nsig_short)))//'). Check disabled.'
               call nek_log_warning(msg, this_module, this_procedure)
            end if

      ! ---- forcing
            if (Wo > 0.0_dp) then
               call t2Dh%init_flow(rec%dpds(1:rec%nf), womersley=Wo)
            else
               call t2Dh%init_flow(rec%dpds(1:rec%nf))
            end if

            if (seed_ .and. rec%jac(1, 1) > 0.0_dp) then
               call t2Dh%set_slope(rec%jac(1, 1))
               write (msg, '(A,A)') 'Outer Jacobian seeded from the record: dQ0/da0= ',
     &            trim(adjustl(rstr(rec%jac(1, 1))))
               call nek_log_message(msg, this_module, this_procedure)
            end if

            call nek_log_message('Forcing restored from '//trim(fname),
     &         this_module, this_procedure)
            if (.not. rec%converged) then
               call nek_log_warning('That record is NOT a converged state.',
     &            this_module, this_procedure)
            end if
            call t2Dh%summary()
            if (present(rec_out)) rec_out = rec
         end subroutine nwt_restore

      !--------------------------------------------------------------------
      !-----     READ HELPERS                                         -----
      !--------------------------------------------------------------------

         pure function agree(a, b) result(l)
      !! Relative comparison, with an absolute floor so two zeros agree.
            real(dp), intent(in) :: a, b
            logical :: l
            l = (abs(a - b) <= check_tol*max(1.0_dp, abs(a), abs(b)))
         end function agree

         subroutine check_scalar(name, from_file, from_run)
      !! A mismatch is fatal. There is no useful way to continue: the forcing
      !! being restored was converged for the other one.
            character(len=*), intent(in) :: name
            real(dp), intent(in) :: from_file, from_run
      ! internal
            character(len=256) :: msg
            if (agree(from_file, from_run)) return
            write (msg, '(A,A,A,A,A)') 'Restart mismatch in ', name, ': record has ',
     &         trim(adjustl(rstr(from_file))), ', this run has '//trim(adjustl(rstr(from_run)))
            call nek_stop_error(msg, this_module, 'check_scalar')
         end subroutine check_scalar

         subroutine split_kv(line, key, val, ok)
      !! Splits 'key = value' at the FIRST '=', then normalises the key at the
      !! first blank or '(' so that a label written for a human to read --
      !! 'dpds  (native)' -- parses as 'dpds'.
            character(len=*), intent(in) :: line
            character(len=*), intent(out) :: key
            character(len=*), intent(out) :: val
            logical, intent(out) :: ok
      ! internal
            character(len=len(line)) :: work
            integer :: ip
            key = ' '; val = ' '; ok = .false.
            if (len_trim(line) == 0) return
            work = adjustl(line)
            if (work(1:1) == '#' .or. work(1:1) == '!') return
            ip = index(work, '=')
            if (ip <= 1) return
            key = trim(work(1:ip - 1))
            val = adjustl(work(ip + 1:))
            if (len_trim(val) == 0) return
            ip = scan(trim(key), ' (')
            if (ip > 1) key = key(1:ip - 1)
            ok = .true.
         end subroutine split_kv

         subroutine bcast_record(rec)
      !! Component by component: a derived type with default initialisation is
      !! not guaranteed contiguous and nek's bcast is byte-based.
            type(t2Dh_nwt_record), intent(inout) :: rec
            call bcast(rec%nf, isize)
            call bcast(rec%nmf, isize)
            call bcast(rec%kharm, isize)
            call bcast(rec%lx1_run, isize)
            call bcast(rec%nelg_run, isize)
            call bcast(rec%np_run, isize)
            call bcast(rec%converged, lsize)
            call bcast(rec%radius, wdsize)
            call bcast(rec%curv_radius, wdsize)
            call bcast(rec%delta, wdsize)
            call bcast(rec%Re, wdsize)
            call bcast(rec%Wo, wdsize)
            call bcast(rec%omega, wdsize)
            call bcast(rec%period, wdsize)
            call bcast(rec%enorm, wdsize)
            call bcast(rec%eabs, wdsize)
            call bcast(rec%noise, wdsize)
            call bcast(rec%dpds, lfc*wdsize)
            call bcast(rec%dpds_helix, lfc*wdsize)
            call bcast(rec%amp, lmfc*wdsize)
            call bcast(rec%phase, lmfc*wdsize)
            call bcast(rec%mflow, lmfc*wdsize)
            call bcast(rec%mf_target, lmfc*wdsize)
            call bcast(rec%mf_err, lmfc*wdsize)
            call bcast(rec%qerr, lmfc*wdsize)
            call bcast(rec%jac, lmfc*lmfc*wdsize)
         end subroutine bcast_record

      end module t2Dh_newton_io