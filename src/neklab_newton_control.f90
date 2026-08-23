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

      end module neklab_newton_control