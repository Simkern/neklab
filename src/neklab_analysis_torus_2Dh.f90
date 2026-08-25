      module neklab_analysis_torus_2Dh
      !! Flow-rate constrained solves on the 2Dh (axisymmetric torus) mesh.
      !!
      !! ONE outer Newton, two wrappers. The steady and pulsatile problems are
      !! not analogues of each other: the steady problem is the nmf = 1 instance
      !! of the pulsatile one. Once the unknown is the amplitude vector
      !! a(1:nmf), nmf = K+1, and the constraint is the flow-rate amplitude
      !! vector mf(1:nmf), the two share everything -- and a Broyden rank-1
      !! update of a 1x1 matrix IS the secant update the steady driver used to
      !! do by hand. flowrate_newton below is that shared core; it contains no
      !! bf_* call and no regime branch, which is what makes it reusable and is
      !! the precondition for the 3D helix driver to adopt it later.
      !!
      !! The regime-specific setup -- buffer lifetime, jacobian type, phase
      !! gauge, converged-orbit recording, unit conversion -- stays visible in
      !! the two wrappers rather than being buried in if (unsteady) branches.
      !!
      !!--------------------------------------------------------------------
      !! THE NORM
      !!--------------------------------------------------------------------
      !!
      !!    ||e|| = max_i |mf_i - target_i| / max(|target_i|, eps_scale*|target_1|)
      !!
      !! Max-norm, because a sum lets component errors cancel -- "converged"
      !! could mean the mean is high by d and the fundamental low by d -- and
      !! because a sum makes the effective per-component requirement depend on
      !! K, so that adding a harmonic silently tightens the tolerance. Scaled by
      !! the target, because the mean and the harmonic amplitudes routinely
      !! differ by an order of magnitude and an unscaled norm would let the mean
      !! dominate. Floored at eps_scale*|target_1|, because a near-zero harmonic
      !! target otherwise divides by nothing; the mean is the natural reference.
      !!
      !! The tolerance argument is therefore RELATIVE and is named rtol_*. The
      !! absolute max-norm is kept alongside it for the quantities that live in
      !! physical units -- the finite-difference noise floor, the quadrature
      !! floor and the Broyden update gate -- and both are logged every step.
         use stdlib_stats_distribution_normal, only: normal => rvs_normal
         use stdlib_strings, only: padl, padr
         use stdlib_optval, only: optval
         use stdlib_linalg, only: diag, eye
         use stdlib_logger, only: information_level, warning_level, debug_level, error_level, all_level, success
         use LightKrylov, only: atol_dp, dp, eigs, svds, save_eigenspectrum
         use LightKrylov, only: initialize_krylov_subspace, orthonormalize_basis, zero_basis, rand_basis
         use LightKrylov, only: linear_combination, innerprod
         use LightKrylov, only: newton, newton_dp_opts, gmres_rdp, gmres_dp_opts
         use LightKrylov_Logger
         use LightKrylov_Timing, only: timer => global_lightkrylov_timer
         use LightKrylov_AbstractSystems, only: abstract_system_rdp
         use LightKrylov_AbstractVectors, only: abstract_vector_rdp
         use neklab_2Dh_axisym, only: nek_advance_2Dh_axisym
         use neklab_vectors
         use neklab_linops
         use neklab_utils
         use neklab_nek_setup
         use neklab_otd
         use neklab_systems
         use neklab_analysis
         use neklab_newton_control
         use neklab_bf_buffer, only: bf_init, bf_finalize_module,
     &                               bf_set_prefix, bf_get_dt_minmax,
     &                               bf_get_nsteps, bf_get_time, bf_summary,
     &                               bf_write_control

         implicit none
         include "SIZE"
         include "TOTAL"
         include "ADJOINT"

         private
         character(len=*), parameter, private :: this_module = 'neklab_analysis_torus_2Dh'

         integer, parameter, private :: lv = lx1*ly1*lz1*lelv

      ! --- outer Newton parameters
         real(dp), parameter, private :: eps_scale = 1.0e-02_dp
      !! Floor of the norm scaling, as a fraction of the mean target.
         real(dp), parameter, private :: step_rel = 0.5_dp
      !! Trust region: maximum |da_i| relative to the current amplitude.
         real(dp), parameter, private :: pred_limit = 2.0_dp
      !! Maximum extrapolation factor of the state predictor.
         real(dp), parameter, private :: inexact_frac = 0.05_dp
      !! Inner tolerance as a fraction of the absolute flow-rate error.
         real(dp), parameter, private :: maxtol = 1.0e-04_dp
      !! Ceiling on the inexact inner tolerance.
         real(dp), parameter, private :: noise_tol_factor = 10.0_dp
      !! Noise floor as a multiple of the inner tolerance.
         integer, parameter, private :: max_fd_retry = 6
      !! Attempts at enlarging a finite-difference step before giving up.
         real(dp), parameter, private :: sing_tol = 1.0e-12_dp
      !! Scaled-determinant threshold for declaring the outer Jacobian singular.

         public :: solve_fixed_point
         public :: flowrate_newton
         public :: steady_flowrate_newton
         public :: unsteady_flowrate_newton
         public :: shift_mflow_phase_upo, upo_taylor_test

      contains

      !==========================================================================
      !     INNER SOLVER WRAPPER
      !==========================================================================

         subroutine solve_fixed_point(sys, X, tol, tol_mode, info, maxiter)
      !! Thin wrapper around the baseline Newton-Krylov fixed-point solver.
            class(abstract_system_rdp), intent(inout) :: sys
            class(abstract_vector_rdp), intent(inout) :: X
            real(dp), intent(in) :: tol
            integer, intent(in) :: tol_mode
            integer, intent(out) :: info
            integer, optional, intent(in) :: maxiter
      ! internal
            type(newton_dp_opts) :: opts
            opts = newton_dp_opts(maxiter=optval(maxiter, 40), ifbisect=.false.)
            if (tol_mode == 1) then
               call newton(sys, X, gmres_rdp, info, atol=tol, options=opts, scheduler=nek_constant_tol)
            else
               call newton(sys, X, gmres_rdp, info, atol=tol, options=opts, scheduler=nek_dynamic_tol)
            end if
         end subroutine solve_fixed_point

      !==========================================================================
      !     THE SHARED OUTER NEWTON
      !==========================================================================

         subroutine flowrate_newton(sys, bf, tol, rtol_mf, tol_mode, maxiter,
     &                              maxiter_inner, if_inexact, jac0, prefix, info)
      !! Segregated outer Newton on the forcing amplitudes, holding the flow
      !! rate at its target. The state is converged by an inner Newton-Krylov
      !! solve at every outer step; the forcing is never part of the inner
      !! unknown.
      !!
      !! ctrl must be configured (init_flow) and its targets set before this is
      !! called, and any harmonic that is exactly zero must have been seeded --
      !! a zero harmonic has no phase, and the frozen-phase directions this
      !! iteration steps along would then be arbitrary.
            class(abstract_system_rdp), intent(inout) :: sys
      !! Inner system. Its jacobian must already be assigned by the caller.
            type(nek_dvector), intent(inout) :: bf
      !! In: initial guess. Out: converged state at the target flow rate.
            real(dp), intent(in) :: tol
      !! Absolute tolerance of the inner Newton-Krylov solver.
            real(dp), intent(in) :: rtol_mf
      !! RELATIVE tolerance on the flow-rate error, in the scaled max-norm.
            integer, optional, intent(in) :: tol_mode
      !! Constant (1, default) or dynamic (2) inner tolerance scheduling.
            integer, optional, intent(in) :: maxiter
      !! Maximum number of outer steps. Default 10.
            integer, optional, intent(in) :: maxiter_inner
      !! Maximum number of inner Newton iterations. Default 40.
            logical, optional, intent(in) :: if_inexact
      !! Loosen the inner tolerance while far from the target. Default .true.
            character(len=*), optional, intent(in) :: jac0
      !! Initial Jacobian: 'seed' (analytic diagonal, no extra solves) or 'fd'
      !! (one inner solve per amplitude). Default 'fd'.
            character(len=3), optional, intent(in) :: prefix
      !! Outpost prefix for the intermediate states. Default 'nwq'.
            integer, optional, intent(out) :: info
      ! internal
            character(len=*), parameter :: this_procedure = 'flowrate_newton'
            character(len=256) :: msg
            character(len=10) :: step_id
            character(len=8) :: jac0_
            character(len=3) :: prefix_
            type(nek_dvector) :: Xold, dX
            real(dp), dimension(lmfc) :: tgt, mf, mf_old, mf_err, dmf, qerr, scal, amp, phase
            real(dp), dimension(:, :), allocatable :: Jac, Jnew
            real(dp), dimension(:), allocatable :: da, da_old, Jda, resid
            integer :: nmf, tol_mode_, maxiter_, maxiter_inner_, inwt, i, j, info_, ierr
            logical :: inexact_, have_pred, first_update, ok
            real(dp) :: enorm, eabs, noise, tol_inner, tol_prev
            real(dp) :: dscale, fac, dmax, nda, nda_old, ratio, denom
            integer, parameter :: pad = 20

      ! ---- optional arguments
            tol_mode_ = optval(tol_mode, 1)
            maxiter_ = optval(maxiter, 10)
            maxiter_inner_ = optval(maxiter_inner, 40)
            inexact_ = optval(if_inexact, .false.)
            jac0_ = optval(jac0, 'fd')
            prefix_ = optval(prefix, 'nwq')
            if (present(info)) info = 0

      ! ---- sizes, targets and the norm scaling
            if (.not. ctrl%is_initialised()) then
               call nek_stop_error('ctrl is not configured. Call ctrl%init_flow first.',
     &            this_module, this_procedure)
            end if
            nmf = ctrl%get_nmf()
            tgt = ctrl%get_target()
            if (abs(tgt(1)) <= atol_dp) then
               call nek_stop_error('The mean flow-rate target is zero: it sets the reference scale '//
     &            'of the norm and cannot vanish.', this_module, this_procedure)
            end if
            scal = 1.0_dp
            do i = 1, nmf
               scal(i) = max(abs(tgt(i)), eps_scale*abs(tgt(1)))
            end do
            allocate (Jac(nmf, nmf), Jnew(nmf, nmf))
            allocate (da(nmf), da_old(nmf), Jda(nmf), resid(nmf))
            Jac = 0.0_dp; Jnew = 0.0_dp
            da = 0.0_dp; da_old = 0.0_dp

      ! ---- stamp logs
            call nek_log_message('Flow-rate Newton configuration:', this_module, this_procedure)
            write (msg, '(3X,A,1X,I0)') padr('unknowns nmf:', pad), nmf
            call nek_log_message(msg, this_module, this_procedure)
            write (msg, '(3X,A,1X,E16.8)') padr('inner tol:', pad), tol
            call nek_log_message(msg, this_module, this_procedure)
            write (msg, '(3X,A,1X,E16.8)') padr('rel. mflow tol:', pad), rtol_mf
            call nek_log_message(msg, this_module, this_procedure)
            write (msg, '(3X,A,A)') padr('tol. scheduling:', pad), padl(merge('constant', 'dynamic ', tol_mode_ == 1), 16)
            call nek_log_message(msg, this_module, this_procedure)
            write (msg, '(3X,A,A)') padr('inexact outer:', pad), padl(merge('yes', 'no ', inexact_), 16)
            call nek_log_message(msg, this_module, this_procedure)
            write (msg, '(3X,A,A)') padr('initial jacobian:', pad), padl(trim(jac0_), 16)
            call nek_log_message(msg, this_module, this_procedure)
            write (msg, '(3X,A,*(1X,F16.10))') padr('target:', pad), (tgt(i), i=1, nmf)
            call nek_log_message(msg, this_module, this_procedure)
            write (msg, '(3X,A,*(1X,F16.10))') padr('norm scale:', pad), (scal(i), i=1, nmf)
            call nek_log_message(msg, this_module, this_procedure)

      !
      ! ---- Baseline solve at the incoming forcing
      !
            call nek_log_message('Baseline solve ...', this_module, this_procedure)
            tol_inner = tol
            tol_prev = tol
            call solve_fixed_point(sys, bf, tol_inner, tol_mode_, info_, maxiter=maxiter_inner_)
            call ctrl%measure_mflow(bf%theta(:, 1), mf, qerr=qerr)
            call errors()
            noise = max(maxval(qerr(1:nmf)), noise_tol_factor*tol)
            write (msg, '(A,1X,E16.8)') padr('flow-rate noise floor:', pad + 4), noise
            call nek_log_message(msg, this_module, this_procedure)
            if (ctrl%is_unsteady()) then
               if (maxval(qerr(1:nmf)) > rtol_mf*minval(scal(1:nmf))) then
                  call nek_log_warning('The requested tolerance is below the measured quadrature '//
     &               'error of the flow-rate accumulator. Lower the CFL target to tighten it.',
     &               this_module, this_procedure)
               end if
            end if
            call ctrl%mflow_summary()
            call log_state(0)

      !
      ! ---- Initial Jacobian
      !
      ! Nothing to build if the incoming forcing already hits the target: the
      ! loop below will exit on its first test, and an FD Jacobian would cost
      ! nmf inner solves for nothing.
            if (enorm < rtol_mf) then
               call nek_log_message('Initial forcing already meets the flow-rate target.',
     &            this_module, this_procedure)
            else if (trim(jac0_) == 'seed') then
               call ctrl%seed_jacobian(Jac, mf)
            else
               call fd_jacobian()
            end if

            call Xold%zero(); call Xold%add(bf)
            have_pred = .false.
            first_update = .true.

      !
      ! ---- Outer iteration
      !
            call nek_log_message('Begin flow-rate Newton iteration', this_module, this_procedure)
            newton_loop: do inwt = 1, maxiter_ + 1
               if (enorm < rtol_mf) then
                  if (tol_inner > tol) then
      ! the last state solve was inexact: confirm at the target tolerance
                     call nek_log_message('Polishing state at target tol ...', this_module, this_procedure)
                     tol_inner = tol
                     call solve_fixed_point(sys, bf, tol_inner, tol_mode_, info_, maxiter=maxiter_inner_)
                     call ctrl%measure_mflow(bf%theta(:, 1), mf, qerr=qerr)
                     call errors()
                  end if
                  if (enorm < rtol_mf) then
                     write (msg, '(A,I0,A)') 'Flow-rate Newton converged after ', inwt - 1, ' step(s).'
                     call nek_log_message(msg, this_module, this_procedure)
                     exit newton_loop
                  end if
               end if
               if (inwt > maxiter_) exit newton_loop
               write (step_id, '("Step ",I3,": ")') inwt

      ! ---- Newton step, Jac da = -mf_err
               call lusolve(nmf, Jac, -mf_err(1:nmf), da, dscale, ierr)
               if (ierr /= 0 .or. abs(dscale) <= sing_tol) then
                  write (msg, '(A,A,E16.8)') step_id, 'outer Jacobian is numerically singular: '//
     &               'scaled det= ', dscale
                  if (trim(jac0_) == 'seed') then
                     call nek_log_warning(step_id//'Rebuilding the Jacobian by finite differences.',
     &                  this_module, this_procedure)
                     call fd_jacobian()
                     call lusolve(nmf, Jac, -mf_err(1:nmf), da, dscale, ierr)
                  end if
                  if (ierr /= 0 .or. abs(dscale) <= sing_tol) then
                     call nek_stop_error(msg, this_module, this_procedure)
                  end if
               end if

      ! ---- trust region. Scale the WHOLE step rather than clipping component
      !      by component, which would rotate the Newton direction.
               call ctrl%get_amp_phase(amp, phase)
               fac = 1.0_dp
               do i = 1, nmf
                  dmax = step_rel*max(abs(amp(i)), eps_scale*abs(amp(1)))
                  if (abs(da(i)) > dmax) fac = min(fac, dmax/abs(da(i)))
               end do
               if (fac < 1.0_dp) then
                  write (msg, '(A,A,E16.8)') step_id, 'step limited by factor ', fac
                  call nek_log_message(msg, this_module, this_procedure)
                  da = fac*da
               end if

      ! ---- state predictor: extrapolate along the solution branch. dX is the
      !      state increment produced by the previous forcing step; scaling it
      !      by the ratio of successive parameter steps, signed by their
      !      alignment, gives a first-order guess of the new solution and
      !      typically saves one or two inner Newton iterations.
               call dX%zero(); call dX%add(bf); call dX%sub(Xold)
               call Xold%zero(); call Xold%add(bf)
               if (have_pred) then
                  nda = norm2(da)
                  nda_old = norm2(da_old)
                  if (nda_old > atol_dp) then
                     ratio = nda/nda_old
                     if (dot_product(da, da_old) < 0.0_dp) ratio = -ratio
                     ratio = max(-1.0_dp, min(pred_limit, ratio))
                     call dX%scal(ratio)
                     call bf%add(dX)
                  end if
               end if

      ! ---- apply the step to the forcing amplitudes
               call ctrl%add_amplitude_step(da)
               write (msg, '(A,A,*(1X,F16.10))') step_id, padl('da:', 10), (da(i), i=1, nmf)
               call nek_log_message(msg, this_module, this_procedure)
               call ctrl%forcing_summary()

      ! ---- inexact inner tolerance: no point resolving the state far below
      !      the accuracy needed to see the current flow-rate error.
               tol_inner = tol
               if (inexact_ .and. enorm > 10.0_dp*rtol_mf) then
                  tol_inner = min(maxtol, max(tol, inexact_frac*eabs))
               end if
               write (msg, '(A,A,1X,E16.8)') step_id, padr('inner tol:', pad), tol_inner
               call nek_log_information(msg, this_module, this_procedure)

      ! ---- converge the state at the new forcing
               mf_old = mf
               call solve_fixed_point(sys, bf, tol_inner, tol_mode_, info_, maxiter=maxiter_inner_)
               call ctrl%measure_mflow(bf%theta(:, 1), mf, qerr=qerr)
               call errors()
               dmf = mf - mf_old
               noise = max(maxval(qerr(1:nmf)), noise_tol_factor*tol_inner)

      ! ---- Broyden ("good") rank-1 update of the outer Jacobian. Same three
      !      safeguards the scalar secant used: refuse a response buried in the
      !      noise, refuse across a tolerance change, and refuse anything that
      !      would break the physically required positive diagonal (each flow
      !      rate is monotone increasing in its own forcing amplitude).
               Jda = matmul(Jac, da)
               if (first_update) then
                  call seed_diagnostic()
                  first_update = .false.
               end if
               if (abs(tol_inner - tol_prev) > 0.1_dp*tol_prev) then
                  call nek_log_warning(step_id//'Tolerance changed. Jacobian kept.',
     &               this_module, this_procedure)
               else if (maxval(abs(dmf(1:nmf))) <= noise) then
                  call nek_log_warning(step_id//'Response below the noise floor. Jacobian kept.',
     &               this_module, this_procedure)
               else
                  denom = dot_product(da, da)
                  if (denom > atol_dp) then
                     resid = dmf(1:nmf) - Jda
                     do i = 1, nmf
                        do j = 1, nmf
                           Jnew(i, j) = Jac(i, j) + resid(i)*da(j)/denom
                        end do
                     end do
                     ok = .true.
                     do i = 1, nmf
                        if (Jnew(i, i) <= 0.0_dp) ok = .false.
                     end do
                     if (ok) then
                        Jac = Jnew
                     else
                        call nek_log_warning(step_id//'Broyden update breaks the diagonal sign. '//
     &                     'Jacobian kept.', this_module, this_procedure)
                     end if
                  end if
               end if

               tol_prev = tol_inner
               da_old = da
               have_pred = .true.
               call log_state(inwt)
               call outpost_dnek(bf, prefix_)
            end do newton_loop

      !
      ! ---- Output
      !
            call nek_log_message('Exiting flow-rate Newton iteration.', this_module, this_procedure)
            if (enorm > rtol_mf) then
               write (msg, '(A,I0,A)') 'Flow rate not converged after ', maxiter_, ' steps.'
               if (present(info)) then
                  info = -1
                  call nek_log_warning(msg, this_module, this_procedure)
               else
                  call nek_stop_error(msg, this_module, this_procedure)
               end if
            end if
            call nek_log_message('OUTPUT:', this_module, this_procedure)
            call ctrl%forcing_summary()
            call ctrl%mflow_summary()
            write (msg, '(3X,A,*(1X,F16.10))') padr('target:', pad), (tgt(i), i=1, nmf)
            call nek_log_message(msg, this_module, this_procedure)
            write (msg, '(3X,A,1X,E16.8,A,E16.8,A)') padr('error:', pad), enorm, '  (abs ', eabs, ')'
            call nek_log_message(msg, this_module, this_procedure)
            write (msg, '(3X,A)') padr('outer jacobian:', pad)
            call nek_log_message(msg, this_module, this_procedure)
            do i = 1, nmf
               write (msg, '(6X,*(1X,E16.8))') (Jac(i, j), j=1, nmf)
               call nek_log_message(msg, this_module, this_procedure)
            end do
      ! the (1,1) entry is the steady resistance an unsteady run will inherit
            if (Jac(1, 1) > 0.0_dp) call ctrl%set_slope(Jac(1, 1))

         contains

            subroutine errors()
      !! Refreshes the scaled and absolute flow-rate errors from mf.
               integer :: ii
               mf_err = 0.0_dp
               do ii = 1, nmf
                  mf_err(ii) = mf(ii) - tgt(ii)
               end do
               enorm = maxval(abs(mf_err(1:nmf))/scal(1:nmf))
               eabs = maxval(abs(mf_err(1:nmf)))
            end subroutine errors

            subroutine log_state(istp)
               integer, intent(in) :: istp
               character(len=512) :: lmsg
               real(dp), dimension(lfc) :: d
               integer :: ii
               call ctrl%get_dpds(d)
               write (lmsg, '(A,I3,A,*(1X,E16.8))') 'FLOWRATE-NEWTON it= ', istp, ' | dpds =',
     &            (d(ii), ii=1, ctrl%get_nf())
               call nek_log_message(lmsg, this_module, this_procedure)
               write (lmsg, '(A,I3,A,*(1X,E16.8))') 'FLOWRATE-NEWTON it= ', istp, ' | mflow=',
     &            (mf(ii), ii=1, nmf)
               call nek_log_message(lmsg, this_module, this_procedure)
               write (lmsg, '(A,I3,A,2(1X,E16.8))') 'FLOWRATE-NEWTON it= ', istp,
     &            ' | err (rel, abs) =', enorm, eabs
               call nek_log_message(lmsg, this_module, this_procedure)
            end subroutine log_state

            subroutine seed_diagnostic()
      !! Compares the measured directional response with the one the initial
      !! Jacobian predicted. This is the self-check on c_inertial: if the
      !! analytic seed is systematically off by a constant factor -- the open
      !! question about whether 0.3 is physics or a definition offset -- it
      !! shows up here, in every run, on the first outer step.
               character(len=512) :: lmsg
               real(dp) :: rn, rd
               integer :: ii
               rn = norm2(dmf(1:nmf))
               rd = norm2(Jda)
               if (rd > atol_dp) then
                  write (lmsg, '(A,F12.5)') 'Jacobian check: |measured|/|predicted| = ', rn/rd
                  call nek_log_message(lmsg, this_module, this_procedure)
               end if
               do ii = 1, nmf
                  if (abs(Jda(ii)) > atol_dp) then
                     write (lmsg, '(A,I0,A,F12.5,A,2(1X,E16.8))') '   component ', ii, ': ratio= ',
     &                  dmf(ii)/Jda(ii), ', (measured, predicted)=', dmf(ii), Jda(ii)
                     call nek_log_message(lmsg, this_module, this_procedure)
                  end if
               end do
            end subroutine seed_diagnostic

            subroutine fd_jacobian()
      !! Full finite-difference Jacobian, one inner solve per amplitude.
      !!
      !!    Jac(j,i) = d(mf_j)/d(a_i)
      !!
      !! ROW = constraint, COLUMN = unknown, because the step solves
      !! Jac da = -mf_err. (mflow_newton in neklab_analysis_torus fills the
      !! transpose and then applies inv(jac) directly. The response is close to
      !! diagonal, so it still converges -- linearly rather than quadratically.
      !! Worth fixing there too.)
      !!
      !! Each column restarts from the same reference state, so the measured
      !! difference is a property of the forcing and not of where the previous
      !! inner solve happened to stop. On exit the forcing and the state are
      !! restored; the measurement cached in ctrl is NOT, but nothing reads it
      !! before the next inner solve refreshes it.
               character(len=256) :: lmsg
               character(len=18) :: coef_id
               type(nek_dvector) :: refv
               real(dp), dimension(lfc) :: d0, dtry
               real(dp), dimension(lmfc) :: mfp, eps
               integer :: ii, jj, icnt, linfo
               logical :: too_small

               call nek_log_message('Building the outer Jacobian by finite differences ...',
     &            this_module, this_procedure)
               call refv%zero(); call refv%add(bf)
               call ctrl%get_dpds(d0)
      ! The mean responds most strongly, so it gets the smallest probe; the
      ! harmonics start an order up.
               eps = 0.0_dp
               eps(1) = min(1.0e-05_dp, 100.0_dp*tol)
               if (nmf > 1) eps(2:nmf) = 10.0_dp*eps(1)

               do ii = 1, nmf
                  write (coef_id, '("Amplitude ",I2,": ")') ii
                  icnt = 0
                  fd_retry: do
      ! Probe in the direction of the root, so the secant is taken on the side
      ! we are heading for.
                     eps(ii) = -sign(eps(ii), mf_err(ii))
                     dtry = ctrl%probe_dpds(ii, eps(ii))
                     call ctrl%set_dpds(dtry)
                     write (lmsg, '(A,A,E16.8)') coef_id, 'probe = ', eps(ii)
                     call nek_log_information(lmsg, this_module, this_procedure)
                     call bf%zero(); call bf%add(refv)
                     call solve_fixed_point(sys, bf, tol, tol_mode_, linfo, maxiter=maxiter_inner_)
                     call ctrl%measure_mflow(bf%theta(:, 1), mfp)
                     dmf = mfp - mf
                     icnt = icnt + 1
      ! A response buried in the quadrature/solver noise carries no gradient
      ! information. Detect that from the MEASUREMENT rather than from whether
      ! the inner Newton happened to iterate: it is the quantity that actually
      ! matters and it does not depend on the inner solver's exit convention.
                     too_small = maxval(abs(dmf(1:nmf))) <= noise
                     if (.not. too_small) exit fd_retry
                     if (icnt >= max_fd_retry) then
                        write (lmsg, '(A,I0,A,I0,A)') 'Amplitude ', ii, ': no measurable response '//
     &                     'after ', icnt, ' attempts. The forcing may be decoupled from it.'
                        call nek_stop_error(lmsg, this_module, this_procedure)
                     end if
                     eps(ii) = 10.0_dp*eps(ii)
                     write (lmsg, '(A,A,E16.8)') coef_id, 'response below the noise floor. '//
     &                  'Enlarged probe = ', abs(eps(ii))
                     call nek_log_message(lmsg, this_module, this_procedure)
                  end do fd_retry
                  write (lmsg, '(A,A,*(1X,E16.8))') coef_id, ' dmflow =', (dmf(jj), jj=1, nmf)
                  call nek_log_message(lmsg, this_module, this_procedure)
                  do jj = 1, nmf
                     Jac(jj, ii) = dmf(jj)/eps(ii)
                  end do
               end do

      ! restore the reference forcing and state
               call ctrl%set_dpds(d0)
               call bf%zero(); call bf%add(refv)
               call nek_log_message('Outer jacobian  d(mflow_row)/d(a_col):', this_module, this_procedure)
               do ii = 1, nmf
                  write (lmsg, '(6X,*(1X,E16.8))') (Jac(ii, jj), jj=1, nmf)
                  call nek_log_message(lmsg, this_module, this_procedure)
               end do
            end subroutine fd_jacobian

         end subroutine flowrate_newton

      !==========================================================================
      !     WRAPPER: STEADY
      !==========================================================================

         subroutine steady_flowrate_newton(sys, bf, dpds, Q_target, tol, rtol_Q,
     &                                     tol_mode, maxiter, if_inexact, dQdf_guess, radius)
      !! Steady fixed point at a prescribed bulk velocity.
      !!
      !! This is the nmf = 1 instance of flowrate_newton: one unknown (the mean
      !! forcing) and one constraint (the mean flow rate). Nothing regime
      !! specific happens here beyond configuring ctrl.
            class(abstract_system_rdp), intent(inout) :: sys
      !! Steady 2Dh torus system (nek_system_torus_2Dh).
            type(nek_dvector), intent(inout) :: bf
      !! In: initial guess. Out: fixed point at the prescribed flow rate.
            real(dp), intent(inout) :: dpds
      !! In: initial forcing. Out: converged forcing. Same in both conventions.
            real(dp), intent(in) :: Q_target
      !! Target bulk velocity, Q = (1/A) int_A u_phi dA.
            real(dp), intent(in) :: tol
      !! Absolute tolerance of the inner Newton-Krylov solver.
            real(dp), intent(in) :: rtol_Q
      !! Relative tolerance on the flow-rate error, |Q - Q_target|/|Q_target|.
            integer, optional, intent(in) :: tol_mode
            integer, optional, intent(in) :: maxiter
            logical, optional, intent(in) :: if_inexact
            real(dp), optional, intent(in) :: dQdf_guess
      !! Initial slope dQ/d(dpds). Default: the Stokes estimate Q/dpds, taken
      !! from the baseline solve.
            real(dp), optional, intent(in) :: radius
      !! Cross-section radius. Measured from the mesh if absent.
      ! internal
            character(len=*), parameter :: this_procedure = 'steady_flowrate_newton'
            character(len=256) :: msg
            real(dp), dimension(lfc) :: d
            real(dp), dimension(1) :: d0, tgt

            d0(1) = dpds
            tgt(1) = Q_target
            call ctrl%init_flow(d0, radius=radius)
            call ctrl%set_target(tgt)
      ! Q(0) = 0 for a steady flow, so a zero forcing yields no slope: probe.
            call ctrl%ensure_nonzero_mean(1.0e-03_dp)
            if (present(dQdf_guess)) then
               call ctrl%set_slope(dQdf_guess)
               write (msg, '(A,E16.8)') 'Using the supplied slope dQ/da= ', dQdf_guess
               call nek_log_message(msg, this_module, this_procedure)
            end if

            call flowrate_newton(sys, bf, tol, rtol_Q, tol_mode=tol_mode, maxiter=maxiter,
     &                           if_inexact=if_inexact, jac0='seed', prefix='nwq')

            call ctrl%get_dpds(d)
            dpds = d(1)
            call set_fldindex('BFQ', 1)
            call outpost_dnek(bf, 'BFQ')
         end subroutine steady_flowrate_newton

      !==========================================================================
      !     WRAPPER: UNSTEADY (PULSATILE)
      !==========================================================================

         subroutine unsteady_flowrate_newton(sys, bf, Wo, kharm, dpds, mflow_target,
     &                                       tol, rtol_mf, tol_mode, maxiter, maxiter_inner,
     &                                       buffer_base, if_save_orbit, if_gauge, jac0, radius)
      !! Pulsatile periodic orbit at a prescribed flow-rate spectrum.
      !!
      !! dpds and mflow_target are in HELIX units on the way in and on the way
      !! out; everything inside is native. Concretely: the forcing amplitudes
      !! agree between the two conventions, while a helix flow-rate amplitude is
      !! HALF the native one, so mflow_target doubles on entry. See the header
      !! of neklab_newton_control.
            class(abstract_system_rdp), intent(inout) :: sys
      !! Pulsatile 2Dh system (nek_system_torus_upo_2Dh). Its jacobian is set here.
            type(nek_dvector), intent(inout) :: bf
      !! In: initial guess for the orbit. Out: converged orbit.
            real(dp), intent(in) :: Wo
      !! Womersley number of the fundamental. Fixed for the whole run.
            integer, intent(in) :: kharm
      !! Number of active harmonics K. nf = 2K+1, nmf = K+1.
            real(dp), dimension(:), intent(inout) :: dpds
      !! In: initial forcing, nf values, helix convention. Out: converged forcing.
            real(dp), dimension(:), intent(in) :: mflow_target
      !! Target flow-rate amplitudes, nmf values, helix normalisation: index 1
      !! is the mean, index k+1 is |mflow_k|, HALF the peak excursion of
      !! harmonic k. These are the numbers a helix deck uses.
            real(dp), intent(in) :: tol
            real(dp), intent(in) :: rtol_mf
      !! RELATIVE tolerance in the scaled max-norm. NOTE: this argument changed
      !! meaning from the absolute sum-norm tolerance of the previous driver.
            integer, optional, intent(in) :: tol_mode
            integer, optional, intent(in) :: maxiter
            integer, optional, intent(in) :: maxiter_inner
            character(len=*), optional, intent(in) :: buffer_base
      !! Filename stem for the baseflow buffer. Default '2dtorus'.
            logical, optional, intent(in) :: if_save_orbit
      !! Re-record the converged orbit under prefix 'b', so the Floquet run has
      !! a stable input the next Newton solve will not overwrite. Default .true.
            logical, optional, intent(in) :: if_gauge
      !! Put the fundamental flow-rate harmonic on a pure cosine before
      !! returning. Default .true.
            character(len=*), optional, intent(in) :: jac0
      !! 'fd' (default) or 'seed'.
            real(dp), optional, intent(in) :: radius
      ! internal
            character(len=*), parameter :: this_procedure = 'unsteady_flowrate_newton'
            character(len=256) :: msg
            type(nek_dvector) :: res
            real(dp), dimension(lfc) :: d_native
            real(dp), dimension(lmfc) :: amp, phase
            logical :: save_orbit_, gauge_
            integer :: nmf, nf, i, k

            save_orbit_ = optval(if_save_orbit, .true.)
            gauge_ = optval(if_gauge, .true.)

      ! ---- sizes and checks
            nf = 2*kharm + 1
            nmf = kharm + 1
            if (kharm < 1) then
               call nek_stop_error('unsteady_flowrate_newton requires K >= 1. Use '//
     &            'steady_flowrate_newton for the steady problem.', this_module, this_procedure)
            end if
            if (nf > lfc) then
               write (msg, '(A,I0,A,I0,A)') 'nf= ', nf, ' > lfc= ', lfc, '. Increase kmax_ctrl.'
               call nek_stop_error(msg, this_module, this_procedure)
            end if
            if (size(dpds) < nf) then
               write (msg, '(A,I0,A,I0)') 'dpds has ', size(dpds), ' entries but nf= ', nf
               call nek_stop_error(msg, this_module, this_procedure)
            end if
            if (size(mflow_target) < nmf) then
               write (msg, '(A,I0,A,I0)') 'mflow_target has ', size(mflow_target), ' entries but nmf= ', nmf
               call nek_stop_error(msg, this_module, this_procedure)
            end if

      ! ---- configure. helix -> native on the forcing and on the targets.
            d_native = helix2native(dpds, nf)
            call ctrl%init_flow(d_native(1:nf), womersley=Wo, radius=radius)
            call ctrl%set_target_helix(mflow_target)

      ! ---- a cold start from a steady solve arrives with every harmonic at
      !      zero, which has no phase to freeze. Seed it.
            call ctrl%get_amp_phase(amp, phase)
            do k = 2, nmf
               if (amp(k) <= atol_dp) then
                  call nek_log_message('Zero harmonic detected: seeding from the inertial law.',
     &               this_module, this_procedure)
                  call ctrl%seed_harmonics()
                  exit
               end if
            end do
            call ctrl%summary()
            if (kharm >= 2) then
               call nek_log_warning('K >= 2: the amplitude spectrum is controlled but the relative '//
     &            'phases between harmonics are not. There are 2K+1 forcing components and one '//
     &            'gauge freedom, hence 2K genuine unknowns against only K+1 amplitude '//
     &            'constraints; the K-1 relative phases are frozen at whatever came in.',
     &            this_module, this_procedure)
            end if

      ! ---- the DRIVER owns the buffer lifetime: the nonlinear map only
      !      consumes it and fails through check_init if this is skipped.
            call bf_init(base=optval(buffer_base, '2dtorus'), write_chunks=.true., min_steps=50)
            call bf_set_prefix('n')
            sys%jacobian = nek_jacobian_torus_upo_2Dh()

            call flowrate_newton(sys, bf, tol, rtol_mf, tol_mode=tol_mode, maxiter=maxiter,
     &                           maxiter_inner=maxiter_inner, jac0=optval(jac0, 'fd'), prefix='nwf')

      ! ---- gauge the time origin
            if (gauge_) call shift_mflow_phase_upo(bf, tol)

      ! ---- record the converged orbit under its own prefix, so the Floquet run
      !      reads a set of chunk files the next Newton solve will not
      !      overwrite. This costs one extra nonlinear pass.
            if (save_orbit_) then
               call nek_log_message('Recording the converged orbit under prefix b ...',
     &            this_module, this_procedure)
               call bf_set_prefix('b')
               call sys%response(bf, res, tol)
               write (msg, '(3X,A,1X,E16.8)') 'converged |F(X)| :', res%norm()
               call nek_log_message(msg, this_module, this_procedure)
               call bf_summary()
               call bf_set_prefix('n')
            end if

      ! ---- hand the forcing back in helix units and leave a restart sidecar
            call ctrl%get_dpds_helix(dpds)
            call ctrl%get_dpds(d_native)
            call bf_write_control(d_native, nf, ctrl%get_omega(), ctrl%get_period())

            call set_fldindex('BFP', 1)
            call outpost_dnek(bf, 'BFP')

            call nek_log_message('OUTPUT:', this_module, this_procedure)
            call ctrl%forcing_summary()
            call ctrl%mflow_summary()
            write (msg, '(3X,A,*(1X,F16.10))') 'helix dpds  :', (dpds(i), i=1, nf)
            call nek_log_message(msg, this_module, this_procedure)
            write (msg, '(3X,A,1X,I0)') 'steps/period:', bf_get_nsteps()
            call nek_log_message(msg, this_module, this_procedure)
         end subroutine unsteady_flowrate_newton

      !====================================================================
      !     PHASE GAUGE
      !====================================================================

         subroutine shift_mflow_phase_upo(bf, tol, cfl_limit)
      !! Uses the time-origin freedom to put the FUNDAMENTAL flow-rate harmonic
      !! on a pure cosine, phase_1 = 0.
      !!
      !! Advancing the state by s and re-referencing the forcing gives the same
      !! physical orbit seen from a shifted origin. In the native convention
      !! BOTH phases rotate the same way,
      !!
      !!    t -> t + s :  phase_k -> phase_k - k*omega*s ,
      !!
      !! so one rotation applied to the forcing array is all that is needed;
      !! s = phase_1/omega zeroes the fundamental's flow-rate phase. (In the
      !! helix convention the two rotate in opposite directions, which is a
      !! convention artefact and the reason this used to need a paragraph of
      !! sign bookkeeping.)
      !!
      !! ONE shift buys ONE phase. For K >= 2 the higher harmonics are rotated
      !! consistently but their phases are NOT zeroed, and that is the honest
      !! outcome: their relative phases are physical, not gauge.
      !! (shift_mflow_phase_torus in neklab_analysis_torus loops over harmonics
      !! applying phase(i)/omega cumulatively, from phases sampled once before
      !! the loop. That is exact for nmf = 2 and does not do what its name says
      !! for nmf >= 3.)
            type(nek_dvector), intent(inout) :: bf
            real(dp), intent(in) :: tol
            real(dp), optional, intent(in) :: cfl_limit
      ! internal
            character(len=*), parameter :: this_procedure = 'shift_mflow_phase_upo'
            character(len=256) :: msg
            real(dp), dimension(lfc) :: mflow
            real(dp), dimension(lmfc) :: amp, phase
            real(dp) :: omega, period, shift, tend
            integer :: k, kc
            real(dp), parameter :: tol_shift = 1.0e-08_dp

            omega = ctrl%get_omega()
            period = ctrl%get_period()
            kc = ctrl%get_kharm()
            if (kc < 1) then
               call nek_log_warning('No harmonics: nothing to gauge.', this_module, this_procedure)
               return
            end if

            call ctrl%get_mflow(mflow, amp, phase)
            write (msg, '(A,*(1X,F16.10))') 'mflow phase (before) =', (phase(k), k=1, kc + 1)
            call nek_log_message(msg, this_module, this_procedure)

            shift = phase(2)/omega
            if (abs(shift) < tol_shift) then
               call nek_log_message('Fundamental is already in phase. Nothing to do.',
     &            this_module, this_procedure)
               return
            end if
      ! Time only runs forward, so a negative shift is taken modulo the period.
            tend = shift
            if (tend < 0.0_dp) tend = period + tend
            write (msg, '(A,E16.8,A,E16.8)') 'Shifting the orbit by s= ', shift, ', integrating for ', tend
            call nek_log_message(msg, this_module, this_procedure)

      ! ---- advance the state by tend at the CURRENT forcing. No recording:
      !      this is a re-gauge, not a residual evaluation.
            call vec2nek(vx, vy, vz, pr, t, bf)
            call setup_nonlinear_solver(variable_dt = .true.,
     &                                  endtime     = tend,
     &                                  solve_temperature = .true.,
     &                                  cfl_limit   = optval(cfl_limit, 0.4_dp),
     &                                  vtol        = tol*0.1,
     &                                  ptol        = tol*0.1)
            time = 0.0_dp
            istep = 0
            do while (lastep == 0)
               istep = istep + 1
               call nek_advance()
            end do
            call nek2vec(bf, vx, vy, vz, pr, t)

      ! ---- re-reference the forcing to the same origin
            call ctrl%rotate_in_time(shift)
            call nek_log_message('Updated forcing:', this_module, this_procedure)
            call ctrl%forcing_summary()
            if (kc >= 2) then
               call nek_log_warning('Only the fundamental has been gauged. Higher harmonics keep '//
     &            'their relative phases, which are physical.', this_module, this_procedure)
            end if
         end subroutine shift_mflow_phase_upo

      !====================================================================
      !     CONSISTENCY CHECK
      !====================================================================

         subroutine upo_taylor_test(sys, X0, tol, nlevels)
      !! Verifies that the replayed Jacobian really is the derivative of the
      !! recorded nonlinear map:
      !!
      !!    e(eps) = || F(X0 + eps*dX) - F(X0) - eps*J*dX ||  ~  O(eps**2)
      !!
      !! This is THE check to run before trusting anything else here. Clean
      !! second order means the buffer replay reproduces the nonlinear
      !! trajectory step for step. A stall at first order, or a plateau, means
      !! the two passes are running on different time grids -- exactly the
      !! failure mode a CFL-adaptive dt can produce if the step sequence drifts
      !! between the residual evaluation and the matvec.
      !!
      !! ORDERING MATTERS: J*dX must be evaluated immediately after F(X0),
      !! because the perturbed residual evaluations overwrite the buffer.
            class(abstract_system_rdp), intent(inout) :: sys
            type(nek_dvector), intent(in) :: X0
            real(dp), intent(in) :: tol
            integer, optional, intent(in) :: nlevels
      ! internal
            character(len=*), parameter :: this_procedure = 'upo_taylor_test'
            character(len=256) :: msg
            type(nek_dvector) :: F0, Fp, dX, JdX, Xp, err
            integer :: nlev, l
            real(dp) :: eps, e, e_prev, order, nrm

            nlev = optval(nlevels, 5)
            sys%jacobian = nek_jacobian_torus_upo_2Dh()

            call nek_log_message('Taylor test on the 2Dh periodic-orbit system:', this_module, this_procedure)

      ! ---- unit-norm random direction
            call dX%rand()
            nrm = dX%norm()
            if (nrm <= atol_dp) call nek_stop_error('Degenerate random direction.', this_module, this_procedure)
            call dX%scal(1.0_dp/nrm)

      ! ---- base residual, then the matvec while the buffer still holds this
      !      trajectory
            call sys%response(X0, F0, tol)
            sys%jacobian%X = X0
            call sys%jacobian%matvec(dX, JdX)
            write (msg, '(3X,A,1X,E16.8,A,I0)') '|F(X0)| = ', F0%norm(), ',  steps/period = ', bf_get_nsteps()
            call nek_log_message(msg, this_module, this_procedure)

            e_prev = 0.0_dp
            do l = 1, nlev
               eps = 10.0_dp**(-l)
               call Xp%zero(); call Xp%add(X0)
               call Xp%axpby(1.0_dp, dX, eps)
               call sys%response(Xp, Fp, tol)
      ! err = Fp - F0 - eps*J*dX
               call err%zero(); call err%add(Fp); call err%sub(F0)
               call err%axpby(1.0_dp, JdX, -eps)
               e = err%norm()
               if (l == 1) then
                  write (msg, '(3X,A,E12.4,A,E16.8)') 'eps= ', eps, '   err= ', e
               else
                  order = 0.0_dp
                  if (e > 0.0_dp .and. e_prev > 0.0_dp) order = log10(e_prev/e)
                  write (msg, '(3X,A,E12.4,A,E16.8,A,F8.3)') 'eps= ', eps, '   err= ', e,
     &               '   observed order= ', order
               end if
               call nek_log_message(msg, this_module, this_procedure)
               e_prev = e
            end do
            call nek_log_message('Expected observed order: 2.0 until the inner tolerance floors it.',
     &         this_module, this_procedure)
         end subroutine upo_taylor_test

      !====================================================================
      !     PRIVATE HELPERS
      !====================================================================

         function helix2native(d, nf) result(a)
      !! Forcing conversion, helix -> native:  a_0 = d_0, a_ck = 2 d_ck,
      !! a_sk = -2 d_sk. The sign flip on the sine is what makes the forcing and
      !! the flow-rate phases rotate together natively.
            real(dp), dimension(:), intent(in) :: d
            integer, intent(in) :: nf
            real(dp), dimension(lfc) :: a
      ! internal
            integer :: k, i
            a = 0.0_dp
            a(1) = d(1)
            do k = 1, (nf - 1)/2
               i = 2*k
               a(i) = 2.0_dp*d(i)
               a(i + 1) = -2.0_dp*d(i + 1)
            end do
         end function helix2native

         subroutine lusolve(n, A, b, x, dscale, ierr)
      !! Dense solve by Gaussian elimination with partial pivoting, plus the
      !! determinant scaled by the product of the row magnitudes.
      !!
      !! Local rather than stdlib's inv/det: n is at most kmax_ctrl+1, the LU is
      !! free once the elimination has run, and an exact-zero determinant test
      !! is useless in floating point -- what matters is the determinant
      !! relative to the size of the rows, which is what dscale is.
            integer, intent(in) :: n
            real(dp), dimension(:, :), intent(in) :: A
            real(dp), dimension(:), intent(in) :: b
            real(dp), dimension(:), intent(out) :: x
            real(dp), intent(out) :: dscale
            integer, intent(out) :: ierr
      ! internal
            real(dp), dimension(n, n) :: M
            real(dp), dimension(n) :: y, rowmax
            real(dp) :: piv, fct, det, scl
            integer :: i, j, k, ip

            ierr = 0
            dscale = 0.0_dp
            x = 0.0_dp
            M = A(1:n, 1:n)
            y = b(1:n)
            do i = 1, n
               rowmax(i) = max(maxval(abs(M(i, :))), atol_dp)
            end do
            det = 1.0_dp
            do k = 1, n - 1
               ip = k
               piv = abs(M(k, k))
               do i = k + 1, n
                  if (abs(M(i, k)) > piv) then
                     piv = abs(M(i, k))
                     ip = i
                  end if
               end do
               if (piv <= 0.0_dp) then
                  ierr = 1
                  return
               end if
               if (ip /= k) then
                  do j = 1, n
                     fct = M(k, j); M(k, j) = M(ip, j); M(ip, j) = fct
                  end do
                  fct = y(k); y(k) = y(ip); y(ip) = fct
                  fct = rowmax(k); rowmax(k) = rowmax(ip); rowmax(ip) = fct
                  det = -det
               end if
               do i = k + 1, n
                  fct = M(i, k)/M(k, k)
                  M(i, k) = 0.0_dp
                  do j = k + 1, n
                     M(i, j) = M(i, j) - fct*M(k, j)
                  end do
                  y(i) = y(i) - fct*y(k)
               end do
            end do
            if (abs(M(n, n)) <= 0.0_dp) then
               ierr = 1
               return
            end if
      ! determinant, scaled by the product of the row magnitudes
            scl = 1.0_dp
            do i = 1, n
               det = det*M(i, i)
               scl = scl*rowmax(i)
            end do
            dscale = det/scl
      ! back substitution
            do i = n, 1, -1
               fct = y(i)
               do j = i + 1, n
                  fct = fct - M(i, j)*x(j)
               end do
               x(i) = fct/M(i, i)
            end do
         end subroutine lusolve

      end module neklab_analysis_torus_2Dh