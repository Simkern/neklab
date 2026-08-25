      module neklab_newton_control
      !! Scalar control unknowns for bordered (flow-rate constrained) systems.
      !!
      !! Holds everything the old 'pipe' object provided for the flow-rate
      !! problem: the control basis, the streamwise forcing, the cross-section
      !! surface integral and the running Fourier accumulator.
      !!
      !! DEPENDENCIES: LightKrylov only. This module sits BELOW neklab_vectors
      !! in the dependency graph (neklab_vectors uses it for the inner-product
      !! weights), so it must NOT use neklab_nek_setup / neklab_utils.
         use stdlib_optval, only: optval
         use LightKrylov, only: dp, atol_dp
         use LightKrylov_Logger
         implicit none
         include "SIZE"
         include "TOTAL"
         private
         character(len=*), parameter, private :: this_module = 'neklab_newton_control'

         integer, parameter, private :: lv = lx1*ly1*lz1*lelv

      !--------------------------------------------------------------------
      !-----     COMPILE-TIME SIZE OF THE CONTROL VECTOR              -----
      !--------------------------------------------------------------------
      ! Maximum number of harmonics. Change and recompile (same convention as
      ! lpert / lelv / nfft). lg is the number of real control unknowns:
      ! one mean + two (cos, sin) per harmonic.
         integer, parameter, public :: kmax_ctrl = 4
         integer, parameter, public :: lg = 2*kmax_ctrl + 1

      !--------------------------------------------------------------------
      !-----     ACTIVE CONFIGURATION (runtime)                       -----
      !--------------------------------------------------------------------
         integer, public :: nctrl = 1
      !! Number of ACTIVE control unknowns, nctrl = 2*kctrl + 1. Public
      !! because neklab_vectors needs it in dot / get_size / rand.
         real(dp), dimension(lg), public :: wg = 1.0_dp
      !! Inner-product weights for the control block. Public for the same reason.

         integer, private :: kctrl = 0
         real(dp), private :: omega_ctrl = 0.0_dp
         real(dp), private :: wres = 1.0_dp
         logical, private :: ctrl_defined = .false.

      ! --- control values
         real(dp), dimension(lg), private :: g_base = 0.0_dp
         real(dp), dimension(lg), private :: g_pert = 0.0_dp
         real(dp), dimension(lg), private :: q_target = 0.0_dp

      ! --- forcing profile
         logical, public :: if_pgrad_profile = .true.
      !! .true.  : f_phi = g * R0/R   (a true mean pressure gradient, irrotational)
      !! .false. : f_phi = g          (a uniform body force)
         real(dp), public :: R0_ctrl = 1.0_dp
      !! Reference radius used to normalise the 1/R profile.

      ! --- surface-integral machinery
         real(dp), dimension(lv), private :: bm_area = 0.0_dp
         real(dp), private :: area_cs = 0.0_dp
         logical, private :: area_defined = .false.

      ! --- running Fourier accumulator
         real(dp), dimension(lg), private :: q_acc = 0.0_dp
         real(dp), private :: t_acc = 0.0_dp

         integer, private :: nf_ = 1
      !! Number of active real forcing components, nf = 2*kctrl + 1.
         real(dp), dimension(lg), private :: dpds_f = 0.0_dp
      !! Forcing coefficients in the HELIX convention. Authoritative; g_base is
      !! derived from these.
         real(dp), private :: womersley_ = 0.0_dp
         real(dp), private :: pulse_T = 0.0_dp
      !! Fundamental period, 2*pi/omega. Fixed for the whole run.
         real(dp), private :: q_lag = 0.0_dp
      !! Previous flow rate, for the trapezoidal accumulator.
   

         public :: init_control, control_summary
         public :: set_control_base, get_control_base
         public :: set_control_pert, clear_control_pert
         public :: set_control_weights, set_control_weights_auto
         public :: get_control_target, set_control_target
         public :: get_nctrl, get_kctrl, get_omega, get_res_scale
         public :: control_basis, control_forcing
         public :: build_area_weights, get_area
         public :: get_flowrate_nek, get_flowrate_pert_nek, get_flowrate_arr
         public :: reset_qfft, accumulate_qfft, extract_qfft
         public :: init_pulsatile, forcing_summary
         public :: set_dpds_fourier, get_dpds_fourier, get_nf
         public :: get_pulsation_period, get_womersley, forcing_amplitude
         public :: reset_qfft_trap, accumulate_qfft_trap
         public :: extract_mflow

      contains

      !====================================================================
      !     CONFIGURATION
      !====================================================================

         subroutine init_control(kharm, omega, target, g0)
      !! Configures the control vector. kharm = 0 gives the steady problem
      !! (one unknown: the mean forcing; one constraint: the mean flow rate).
            integer, intent(in) :: kharm
      !! Number of active harmonics K. nctrl = 2*K + 1.
            real(dp), optional, intent(in) :: omega
      !! Fundamental angular frequency. Required for kharm > 0.
            real(dp), dimension(:), optional, intent(in) :: target
      !! Target flow-rate coefficients (nctrl values).
            real(dp), dimension(:), optional, intent(in) :: g0
      !! Initial forcing coefficients (nctrl values).
      ! internal
            character(len=*), parameter :: this_procedure = 'init_control'
            character(len=128) :: msg
            if (kharm < 0 .or. kharm > kmax_ctrl) then
               write (msg, '(A,I0,A,I0)') 'Invalid kharm= ', kharm, '. Must be in [0, kmax_ctrl] with kmax_ctrl= ', kmax_ctrl
               call stop_error(msg, this_module, this_procedure)
            end if
            kctrl = kharm
            nctrl = 2*kctrl + 1
            omega_ctrl = optval(omega, 0.0_dp)
            if (kctrl > 0 .and. abs(omega_ctrl) <= atol_dp) then
               call stop_error('kharm > 0 requires a non-zero omega.', this_module, this_procedure)
            end if
            g_base = 0.0_dp
            g_pert = 0.0_dp
            q_target = 0.0_dp
            wg = 1.0_dp
            if (present(target)) call set_control_target(target)
            if (present(g0)) call set_control_base(g0)
      ! constraint rows are scaled to read a bulk-velocity error
            call build_area_weights()
            wres = 1.0_dp/area_cs
            ctrl_defined = .true.
            call control_summary()
         end subroutine init_control

         subroutine control_summary()
            character(len=*), parameter :: this_procedure = 'control_summary'
            character(len=256) :: msg
            integer :: i
            write (msg, '(A,I0,A,I0,A,I0)') 'control: K= ', kctrl, ', nctrl= ', nctrl, ', lg= ', lg
            call logger%log_message(msg, this_module, this_procedure)
            write (msg, '(A,E16.8,A,E16.8)') '   omega= ', omega_ctrl, ', area= ', area_cs
            call logger%log_message(msg, this_module, this_procedure)
            write (msg, '(A,*(1X,E14.6))') '   g      =', (g_base(i), i=1, nctrl)
            call logger%log_message(msg, this_module, this_procedure)
            write (msg, '(A,*(1X,E14.6))') '   Qtarget=', (q_target(i), i=1, nctrl)
            call logger%log_message(msg, this_module, this_procedure)
            write (msg, '(A,*(1X,E14.6))') '   weights=', (wg(i), i=1, nctrl)
            call logger%log_message(msg, this_module, this_procedure)
         end subroutine control_summary

         integer function get_nctrl() result(n)
            n = nctrl
         end function get_nctrl

         integer function get_kctrl() result(k)
            k = kctrl
         end function get_kctrl

         real(dp) function get_omega() result(w)
            w = omega_ctrl
         end function get_omega

         real(dp) function get_res_scale() result(w)
      !! Scaling applied to the constraint rows of the residual.
            w = wres
         end function get_res_scale

      !====================================================================
      !     CONTROL VALUES
      !====================================================================

         subroutine set_control_base(g)
            real(dp), dimension(:), intent(in) :: g
            g_base = 0.0_dp
            g_base(1:min(size(g), nctrl)) = g(1:min(size(g), nctrl))
            call bcast(g_base, lg*wdsize)
         end subroutine set_control_base

         function get_control_base() result(g)
            real(dp), dimension(lg) :: g
            g = g_base
         end function get_control_base

         subroutine set_control_pert(g)
      !! Sets the FORCING PERTURBATION used by the linearised solver. Must be
      !! cleared after every Jacobian matvec, or the next nonlinear solve
      !! inherits a spurious source.
            real(dp), dimension(:), intent(in) :: g
            g_pert = 0.0_dp
            g_pert(1:min(size(g), nctrl)) = g(1:min(size(g), nctrl))
            call bcast(g_pert, lg*wdsize)
         end subroutine set_control_pert

         subroutine clear_control_pert()
            g_pert = 0.0_dp
         end subroutine clear_control_pert

         subroutine set_control_target(target)
            real(dp), dimension(:), intent(in) :: target
            q_target = 0.0_dp
            q_target(1:min(size(target), nctrl)) = target(1:min(size(target), nctrl))
         end subroutine set_control_target

         function get_control_target() result(target)
            real(dp), dimension(lg) :: target
            target = q_target
         end function get_control_target

         subroutine set_control_weights(w)
            real(dp), dimension(:), intent(in) :: w
            wg = 1.0_dp
            wg(1:min(size(w), nctrl)) = w(1:min(size(w), nctrl))
         end subroutine set_control_weights

         subroutine set_control_weights_auto(xnorm, gscale)
      !! Sets wg so that a unit perturbation of each control produces a control
      !! block contribution to the norm comparable to the state block.
      !!
      !! The norm contribution is wg(i)*g(i)**2, so wg = (xnorm/gscale)**2.
      !! gscale(i) should be the magnitude of control i that produces an O(1)
      !! relative change of the state -- for the mean forcing that is simply
      !! |g_1|; for harmonic k the response rolls off with k (Womersley), so a
      !! single common value will be badly wrong at k = 3.
            real(dp), intent(in) :: xnorm
            real(dp), dimension(:), intent(in) :: gscale
      ! internal
            character(len=*), parameter :: this_procedure = 'set_control_weights_auto'
            character(len=256) :: msg
            integer :: i
            wg = 1.0_dp
            do i = 1, min(size(gscale), nctrl)
               if (abs(gscale(i)) > atol_dp) then
                  wg(i) = (xnorm/gscale(i))**2
               else
                  wg(i) = 1.0_dp
                  write (msg, '(A,I0,A)') 'Zero scale for control ', i, '. Weight set to unity.'
                  call logger%log_warning(msg, this_module, this_procedure)
               end if
            end do
            write (msg, '(A,*(1X,E14.6))') 'control weights =', (wg(i), i=1, nctrl)
            call logger%log_message(msg, this_module, this_procedure)
         end subroutine set_control_weights_auto

      !====================================================================
      !     CONTROL BASIS AND FORCING
      !====================================================================

         real(dp) function control_basis(i, tval) result(phi)
      !! i-th temporal basis function: 1, cos(wt), sin(wt), cos(2wt), ...
            integer, intent(in) :: i
            real(dp), intent(in) :: tval
      ! internal
            integer :: k
            if (i == 1) then
               phi = 1.0_dp
            else
               k = i/2
               if (mod(i, 2) == 0) then
                  phi = cos(k*omega_ctrl*tval)
               else
                  phi = sin(k*omega_ctrl*tval)
               end if
            end if
         end function control_basis

         real(dp) function control_forcing(tval, jp_, y) result(f)
      !! Streamwise source term for userq. Call as
      !!
      !!    qvol = control_forcing(time, jp, y)
      !!
      !! jp = 0 returns the base forcing, jp > 0 the perturbation forcing.
      !! Note that Nek's makeqp evaluates userq at t_n (it shifts 'time' by
      !! -dt internally), which is what the linearisation requires.
            real(dp), intent(in) :: tval
            integer, intent(in) :: jp_
            real(dp), intent(in) :: y
      ! internal
            integer :: i
            real(dp) :: amp
            amp = 0.0_dp
            if (jp_ == 0) then
               do i = 1, nctrl
                  amp = amp + g_base(i)*control_basis(i, tval)
               end do
            else
               do i = 1, nctrl
                  amp = amp + g_pert(i)*control_basis(i, tval)
               end do
            end if
            if (if_pgrad_profile) then
      ! f_phi = -(1/R) dp/dphi with p = -g*R0*phi. Irrotational, i.e. a genuine
      ! mean pressure gradient. A uniform f_phi is NOT (curl f = f/R).
               f = amp*R0_ctrl/y
            else
               f = amp
            end if
         end function control_forcing

      !====================================================================
      !     SURFACE INTEGRALS OVER THE CROSS-SECTION
      !====================================================================

         subroutine build_area_weights(force)
      !! Mass matrix for integrals over the phi = const cross-section, whose
      !! area element is dA = dR dz. If the mesh runs with ifaxis = .true.,
      !! bm1 carries an extra radial weight which is stripped here; if the
      !! torus metric lives entirely in the coefficient arrays of
      !! neklab_2Dh_axisym on a plain 2D mesh, bm1 is already correct.
            logical, optional, intent(in) :: force
      ! internal
            character(len=*), parameter :: this_procedure = 'build_area_weights'
            real(dp), external :: glsum
            integer :: n
            character(len=128) :: msg
            if (area_defined .and. .not. optval(force, .false.)) return
            n = lx1*ly1*lz1*nelv
            if (ifaxis) then
               call invcol3(bm_area, bm1, ym1, n)
            else
               call copy(bm_area, bm1, n)
            end if
            area_cs = glsum(bm_area, n)
            area_defined = .true.
            write (msg, '(A,E16.8)') 'Cross-section area A = ', area_cs
            call logger%log_information(msg, this_module, this_procedure)
         end subroutine build_area_weights

         real(dp) function get_area() result(area)
            call build_area_weights()
            area = area_cs
         end function get_area

         real(dp) function get_flowrate_arr(theta) result(Q)
      !! Q = int_A theta dA for an arbitrary field array.
            real(dp), dimension(lx1*ly1*lz1*lelt), intent(in) :: theta
      ! internal
            real(dp), external :: glsc2
            integer :: n
            call build_area_weights()
            n = lx1*ly1*lz1*nelv
            Q = glsc2(theta, bm_area, n) / area_cs
         end function get_flowrate_arr

         real(dp) function get_flowrate_nek() result(Q)
      !! Flow rate of the current nonlinear field t (which carries u_phi).
            Q = get_flowrate_arr(t(1, 1, 1, 1, 1))
         end function get_flowrate_nek

         real(dp) function get_flowrate_pert_nek(jp_) result(Q)
      !! Flow rate of the current perturbation field tp.
            integer, optional, intent(in) :: jp_
            Q = get_flowrate_arr(tp(1, 1, optval(jp_, 1)))
         end function get_flowrate_pert_nek

      !====================================================================
      !     RUNNING FOURIER ACCUMULATOR
      !====================================================================
      !
      !  The constraint is the time average of the flow rate against each basis
      !  function over the integration horizon:
      !
      !     c_1     = (1/T) int_0^T Q(t) dt
      !     c_{2k}  = (2/T) int_0^T Q(t) cos(k w t) dt
      !     c_{2k+1}= (2/T) int_0^T Q(t) sin(k w t) dt
      !
      !  For K = 0 this is just the mean of Q over the horizon, which equals
      !  Q itself once the state is a genuine fixed point. Away from the
      !  solution the two differ, which is fine: the residual and its Jacobian
      !  are consistent with each other, and the root is unchanged.

         subroutine reset_qfft()
            q_acc = 0.0_dp
            t_acc = 0.0_dp
         end subroutine reset_qfft

         subroutine accumulate_qfft(Q, tval, dtn)
      !! Call once per timestep, after nek_advance, with the flow rate of the
      !! current field, the current time and the timestep just taken.
            real(dp), intent(in) :: Q
            real(dp), intent(in) :: tval
            real(dp), intent(in) :: dtn
      ! internal
            integer :: i
            t_acc = t_acc + dtn
            do i = 1, nctrl
               q_acc(i) = q_acc(i) + dtn*Q*control_basis(i, tval)
            end do
         end subroutine accumulate_qfft

         subroutine extract_qfft(coeffs)
      !! Normalises the accumulator. Inactive slots are returned as zero.
            real(dp), dimension(lg), intent(out) :: coeffs
      ! internal
            character(len=*), parameter :: this_procedure = 'extract_qfft'
            integer :: i
            coeffs = 0.0_dp
            if (t_acc <= 0.0_dp) then
               call stop_error('Empty accumulator: reset_qfft/accumulate_qfft were not called.', this_module, this_procedure)
            end if
            coeffs(1) = q_acc(1)/t_acc
            do i = 2, nctrl
               coeffs(i) = 2.0_dp*q_acc(i)/t_acc
            end do
         end subroutine extract_qfft

               !====================================================================
      !     PULSATILE CONFIGURATION
      !====================================================================
 
         subroutine init_pulsatile(womersley, kharm, dpds)
      !! Configures the pulsatile control from a Womersley number, exactly as
      !! neklab_helix does (helix_utils.f90:170-171):
      !!
      !!    omega = Wo**2 * nu ,   T = 2*pi/omega
      !!
      !! with nu = cpfld(1,1). Wo is fixed for the whole run, so this is called
      !! once by the driver and the period never changes afterwards.
            real(dp), intent(in) :: womersley
            integer, intent(in) :: kharm
      !! Number of active harmonics K. nf = 2*K + 1.
            real(dp), dimension(:), intent(in) :: dpds
      !! Initial forcing coefficients in the helix convention (nf values).
      ! internal
            character(len=*), parameter :: this_procedure = 'init_pulsatile'
            character(len=256) :: msg
            real(dp) :: omega, pi
            pi = 4.0_dp*atan(1.0_dp)
            womersley_ = womersley
            if (abs(womersley_) <= atol_dp) then
               call stop_error('init_pulsatile requires a non-zero Womersley number. '//
     &            'Use init_control directly for the steady problem.', this_module, this_procedure)
            end if
            omega = (womersley_**2)*cpfld(1, 1)
            pulse_T = 2.0_dp*pi/omega
            nf_ = 2*kharm + 1
            if (size(dpds) < nf_) then
               write (msg, '(A,I0,A,I0)') 'dpds has ', size(dpds), ' entries but nf = ', nf_
               call stop_error(msg, this_module, this_procedure)
            end if
            call init_control(kharm, omega=omega)
            call set_dpds_fourier(dpds)
            write (msg, '(A,E16.8,A,E16.8,A,E16.8)') 'pulsatile: Wo= ', womersley_,
     &         ', nu= ', cpfld(1, 1), ', omega= ', omega
            call logger%log_message(msg, this_module, this_procedure)
            write (msg, '(A,E24.16)') '   period T = ', pulse_T
            call logger%log_message(msg, this_module, this_procedure)
         end subroutine init_pulsatile
 
         real(dp) function get_pulsation_period() result(T)
            T = pulse_T
         end function get_pulsation_period
 
         real(dp) function get_womersley() result(Wo)
            Wo = womersley_
         end function get_womersley
 
         integer function get_nf() result(n)
            n = nf_
         end function get_nf
 
      !====================================================================
      !     FORCING IN THE HELIX CONVENTION
      !====================================================================
 
         subroutine set_dpds_fourier(dpds)
      !! Stores the forcing in the helix convention and pushes the equivalent
      !! internal coefficients to the base control, which is what userq sees
      !! through control_forcing.
            real(dp), dimension(:), intent(in) :: dpds
      ! internal
            real(dp), dimension(lg) :: g
            integer :: i, k, n
            n = min(size(dpds), nf_)
            dpds_f = 0.0_dp
            dpds_f(1:n) = dpds(1:n)
            call bcast(dpds_f, lg*wdsize)
            g = 0.0_dp
            g(1) = dpds_f(1)
            do k = 1, kctrl
               i = 2*k
               g(i) = 2.0_dp*dpds_f(i)
               g(i + 1) = -2.0_dp*dpds_f(i + 1)
            end do
            call set_control_base(g)
         end subroutine set_dpds_fourier
 
         subroutine get_dpds_fourier(dpds, phase)
      !! Returns the forcing coefficients and, optionally, the phase of each
      !! harmonic in the helix convention, phase_k = atan2(dpds(2k+1), dpds(2k))
      !! (helix_gs.f90:78-90). Index 1 of phase is the mean, whose phase is zero
      !! by construction.
            real(dp), dimension(:), intent(out) :: dpds
            real(dp), dimension(:), allocatable, optional, intent(out) :: phase
      ! internal
            integer :: i, k, n
            n = min(size(dpds), lg)
            dpds = 0.0_dp
            dpds(1:n) = dpds_f(1:n)
            if (present(phase)) then
               allocate (phase(kctrl + 1))
               phase = 0.0_dp
               do k = 1, kctrl
                  i = 2*k
                  phase(k + 1) = atan2(dpds_f(i + 1), dpds_f(i))
               end do
            end if
         end subroutine get_dpds_fourier
 
         real(dp) function forcing_amplitude(tval) result(f)
      !! Instantaneous streamwise forcing amplitude in the helix convention.
      !! Diagnostic only: the forcing actually applied comes from
      !! control_forcing, which evaluates the equivalent internal coefficients.
            real(dp), intent(in) :: tval
      ! internal
            integer :: k, i
            f = dpds_f(1)
            do k = 1, kctrl
               i = 2*k
               f = f + 2.0_dp*(dpds_f(i)*cos(k*omega_ctrl*tval) - dpds_f(i + 1)*sin(k*omega_ctrl*tval))
            end do
         end function forcing_amplitude
 
         subroutine forcing_summary()
            character(len=*), parameter :: this_procedure = 'forcing_summary'
            character(len=512) :: msg
            real(dp), dimension(lg) :: d
            real(dp), dimension(:), allocatable :: ph
            integer :: i, k
            call get_dpds_fourier(d, ph)
            write (msg, '(A,*(1X,F16.10))') 'dpds      =', (d(i), i=1, nf_)
            call logger%log_message(msg, this_module, this_procedure)
            write (msg, '(A,*(1X,F16.10))') 'amplitude =', d(1), (2.0_dp*sqrt(d(2*k)**2 + d(2*k + 1)**2), k=1, kctrl)
            call logger%log_message(msg, this_module, this_procedure)
            write (msg, '(A,*(1X,F16.10))') 'phase     =', (ph(k), k=1, kctrl + 1)
            call logger%log_message(msg, this_module, this_procedure)
         end subroutine forcing_summary
 
      !====================================================================
      !     TRAPEZOIDAL FLOW-RATE ACCUMULATOR
      !====================================================================
      !
      !  accumulate_qfft uses the right-endpoint rule, which is O(dt) on a
      !  non-uniform grid. With a variable timestep and as few as 50 steps per
      !  period that is a percent-level error on the harmonics -- far above the
      !  Newton tolerance, and it would put a floor under the flow-rate solve.
      !  The trapezoidal variant below is O(dt**2) and mirrors what helix does
      !  (helix_mflow_fft.f90:30-65).
      !
      !  REVERSE FLOW. When Q changes sign within a step, a node is placed at
      !  the zero crossing of the linear interpolant and the rule is applied on
      !  each half. For the mean this makes no difference (the trapezoid of a
      !  linear function is exact either way, and the split form is
      !  algebraically identical); for the harmonics it is a genuine refinement,
      !  since Q*phi is not linear and the extra node helps.
      !
      !  Two corrections versus helix_mflow_fft.f90:34-49, which this replaces:
      !
      !    * the crossing fraction. Helix uses dt0 = -(Q - Q_old)/Q_old, which
      !      is not a fraction of the interval at all: for Q_old = 1, Q = -1 it
      !      returns 2. The linear interpolant vanishes at
      !      s = Q_old/(Q_old - Q), which is 0.5 for that case.
      !
      !    * the factor 1/2. Helix applies it to the mean but not to the
      !      harmonic terms, so every step containing a crossing contributes
      !      twice what it should to every harmonic.
      !
      !  Both only fire on a sign change, which is why they have survived: with
      !  no reverse flow the branch is never taken.
 
         subroutine reset_qfft_trap(Q0)
      !! Starts a trapezoidal accumulation. Q0 is the flow rate at t = 0, i.e.
      !! of the initial condition, before the first step is taken.
            real(dp), intent(in) :: Q0
            q_acc = 0.0_dp
            t_acc = 0.0_dp
            q_lag = Q0
         end subroutine reset_qfft_trap
 
         subroutine accumulate_qfft_trap(Q, tval, dtn)
      !! Call once per timestep, AFTER nek_advance, with the flow rate of the
      !! current field, the current time and the timestep just taken. tval is
      !! the time at the END of the step, so the interval is [tval-dtn, tval].
            real(dp), intent(in) :: Q
            real(dp), intent(in) :: tval
            real(dp), intent(in) :: dtn
      ! internal
            integer :: i
            real(dp) :: told, s
            told = tval - dtn
            t_acc = t_acc + dtn
            if (Q*q_lag < 0.0_dp) then
      ! Sign change: split at the zero of the linear interpolant. The integrand
      ! vanishes there, so only the two outer endpoints contribute and the
      ! weights are the sub-interval lengths.
               s = q_lag/(q_lag - Q)
 
               do i = 1, nctrl
                  q_acc(i) = q_acc(i) + 0.5_dp*dtn*(s*q_lag*control_basis(i, told)
     &                                              + (1.0_dp - s)*Q*control_basis(i, tval))
               end do
            else
               do i = 1, nctrl
                  q_acc(i) = q_acc(i) + 0.5_dp*dtn*(q_lag*control_basis(i, told) + Q*control_basis(i, tval))
               end do
            end if
            q_lag = Q
         end subroutine accumulate_qfft_trap
 
      !====================================================================
      !     FLOW RATE IN THE HELIX CONVENTION
      !====================================================================
 
         subroutine extract_mflow(mflow, amplitude, phase)
      !! Normalises the accumulator and converts to the helix convention.
      !! Note that amplitude(k+1) is |mflow_k|, i.e. HALF the peak-to-mean
      !! excursion of harmonic k -- the same quantity helix reports and the
      !! same one its targets are given in.
            real(dp), dimension(lg), intent(out) :: mflow
            real(dp), dimension(:), allocatable, optional, intent(out) :: amplitude
            real(dp), dimension(:), allocatable, optional, intent(out) :: phase
      ! internal
            real(dp), dimension(lg) :: coeffs
            integer :: i, k
            call extract_qfft(coeffs)
            mflow = 0.0_dp
            mflow(1) = coeffs(1)
            do i = 2, nctrl
               mflow(i) = 0.5_dp*coeffs(i)
            end do
            if (present(amplitude)) then
               allocate (amplitude(kctrl + 1))
               amplitude = 0.0_dp
               amplitude(1) = mflow(1)
               do k = 1, kctrl
                  i = 2*k
                  amplitude(k + 1) = sqrt(mflow(i)**2 + mflow(i + 1)**2)
               end do
            end if
            if (present(phase)) then
               allocate (phase(kctrl + 1))
               phase = 0.0_dp
               do k = 1, kctrl
                  i = 2*k
                  phase(k + 1) = atan2(mflow(i + 1), mflow(i))
               end do
            end if
         end subroutine extract_mflow


      end module neklab_newton_control