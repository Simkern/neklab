      module neklab_newton_control
      !! Owner of the flow-rate control state for the 2Dh torus.
      !!
      !! Everything the flow-rate problem needs that is not a Nek field lives
      !! in the single instance `ctrl` of the type below: the streamwise
      !! forcing, the pulsation time scales, the cross-section geometry, the
      !! flow-rate targets, the running Fourier accumulator and the steady
      !! resistance inherited by an unsteady run. This mirrors the way
      !! neklab_helix keeps all of the helix state in `pipe`.
      !!
      !! There is no bordered vector: the forcing is NOT part of the Newton
      !! unknown. It is driven by the segregated outer iteration in
      !! neklab_analysis_torus_2Dh, which is why nothing here is exposed to
      !! neklab_vectors and why this module is free to sit above
      !! neklab_nek_setup in the dependency graph.
      !!
      !!--------------------------------------------------------------------
      !! FOURIER CONVENTION  (native)
      !!--------------------------------------------------------------------
      !!
      !!   f(t) = a_0 + sum_k [ a_ck cos(k w t) + a_sk sin(k w t) ]
      !!   Q(t) = q_0 + sum_k [ q_ck cos(k w t) + q_sk sin(k w t) ]
      !!
      !!   amp_k   = hypot(._c, ._s)      ! PEAK EXCURSION, both quantities
      !!   phase_k = atan2(._s, ._c)      ! same sign, both quantities
      !!   t -> t+s :  phase_k -> phase_k - k w s   ! same direction, both
      !!
      !! ONE convention for the forcing and the flow rate, unlike the helix
      !! code this was ported from, which projects the forcing onto exp(-ikwt)
      !! and the flow rate onto exp(+ikwt) and reports the forcing amplitude
      !! as a peak excursion but the flow-rate amplitude as half of one. The
      !! packing is dpds = (a_0, a_c1, a_s1, a_c2, a_s2, ...), stored exactly
      !! as it is evaluated: there is no separate internal 'g' representation.
      !!
      !! Helix decks are still readable through set_dpds_helix / get_dpds_helix
      !! and mflow_from_helix / mflow_to_helix, which are the ONLY places the
      !! old convention appears:
      !!
      !!   a_0 = d_0 ,  a_ck =  2 d_ck ,  a_sk = -2 d_sk      (forcing)
      !!   q_0 = m_0 ,  q_ck =  2 m_ck ,  q_sk =  2 m_sk      (flow rate)
      !!
      !! so a helix flow-rate amplitude is HALF the native one, while the two
      !! forcing amplitudes agree.
      !!
      !!--------------------------------------------------------------------
      !! FORCING PROFILE
      !!--------------------------------------------------------------------
      !!
      !!   f_phi = amplitude(t) / R
      !!
      !! which is irrotational, i.e. a genuine mean pressure gradient (a
      !! uniform f_phi is not: curl f = f/R). This matches helix exactly --
      !! helix_utils divides by curv_radius and then multiplies by
      !! fshape = 1/(1 + delta*rr*sin(alpha)), whose product is 1/R -- so the
      !! same dpds means the same physical forcing in both codes. The earlier
      !! 2Dh port carried a spurious reference radius R0 in front; it is gone.
         use stdlib_optval, only: optval
         use LightKrylov, only: dp, atol_dp
         use LightKrylov_Logger
         use neklab_nek_setup, only: nek_log_message, nek_log_information,
     &                               nek_log_warning, nek_log_debug, nek_stop_error
         implicit none
         include "SIZE"
         include "TOTAL"
         private
         character(len=*), parameter, private :: this_module = 'neklab_newton_control'

         integer, parameter :: lv = lx1*ly1*lz1*lelv

      !--------------------------------------------------------------------
      !-----     COMPILE-TIME SIZES                                   -----
      !--------------------------------------------------------------------
      ! Change and recompile, same convention as lpert / lelv / nfft.
         integer, parameter, public :: kmax_ctrl = 4
      !! Maximum number of harmonics K.
         integer, parameter, public :: lfc = 2*kmax_ctrl + 1
      !! Number of real forcing components, one mean + two per harmonic.
         integer, parameter, public :: lmfc = kmax_ctrl + 1
      !! Number of amplitude unknowns / constraints, one mean + one per harmonic.

         real(dp), parameter, public :: c_inertial = 0.3_dp
      !! Empirical prefactor of the inertial forcing law used to seed the
      !! outer Jacobian,
      !!
      !!    |a_k| = c_inertial * k*omega*|q_k| / delta
      !!
      !! calibrated at Q_mean = 1 over a range of omega and delta, for which
      !! all amplitude curves collapse onto a line of this slope once scaled
      !! by omega/delta. NOTE: the angular-momentum balance
      !! R_c d(u_phi)/dt = f_phi with f_phi = a/R gives a prefactor of ONE, and
      !! viscosity can only push the required forcing UP (impedances add in
      !! quadrature), so a measured value below one is almost certainly a
      !! definition offset between the plotted axes and the symbols above --
      !! peak versus half excursion, or r versus D as the length scale. The
      !! value is therefore treated as a calibration, not as physics, and the
      !! first Broyden update logs the measured-to-seeded ratio so that a
      !! systematic factor shows up in every run instead of in a plot. See
      !! seed_jacobian.

      !--------------------------------------------------------------------
      !-----     THE CONTROL TYPE                                     -----
      !--------------------------------------------------------------------

         type, public :: nek_control
            private
      ! --- regime. NOTE the polarity: this module has if_unsteady, neklab_helix
      !     has is_steady(). Do not mix them up.
            logical :: if_unsteady = .false.
            integer :: kharm = 0
      !! Number of active harmonics K.
            integer :: nf = 1
      !! Number of active forcing components, 2K+1.
            integer :: nmf = 1
      !! Number of amplitude unknowns / constraints, K+1.
      ! --- time scales
            real(dp) :: womersley = 0.0_dp
            real(dp) :: omega = 0.0_dp
            real(dp) :: period = 0.0_dp
      ! --- forcing, native convention
            real(dp), dimension(lfc) :: dpds = 0.0_dp
      ! --- flow rate: coefficients, amplitudes, phases and targets
            real(dp), dimension(lfc) :: qfour = 0.0_dp
            real(dp), dimension(lmfc) :: mf = 0.0_dp
            real(dp), dimension(lmfc) :: mf_phase = 0.0_dp
            real(dp), dimension(lmfc) :: mf_qerr = 0.0_dp
            real(dp), dimension(lmfc) :: mf_target = 0.0_dp
            logical :: mf_extracted = .false.
            logical :: target_defined = .false.
      ! --- Fourier accumulator. q_acc is the trapezoidal rule on every step,
      !     q_crs the same rule on pairs of steps; their difference is a
      !     Richardson estimate of the quadrature error, which is what floors
      !     the outer tolerance (see close_mflow).
            real(dp), dimension(lfc) :: q_acc = 0.0_dp
            real(dp), dimension(lfc) :: q_crs = 0.0_dp
            real(dp) :: t_acc = 0.0_dp
            real(dp) :: q_lag = 0.0_dp
            real(dp) :: t_lag = 0.0_dp
            real(dp) :: q_lag2 = 0.0_dp
            real(dp) :: t_lag2 = 0.0_dp
            integer :: nacc = 0
            logical :: accumulating = .false.
      ! --- cross-section geometry
            real(dp), dimension(lv) :: bm_area = 0.0_dp
            real(dp) :: area = 0.0_dp
            real(dp) :: curv_radius = 0.0_dp
            real(dp) :: radius = 0.0_dp
            real(dp) :: delta = 0.0_dp
            logical :: area_defined = .false.
      ! --- steady resistance dQ_0/da_0, inherited by an unsteady run
            real(dp) :: gslope = 0.0_dp
            logical :: gslope_defined = .false.
      ! --- sanity
            logical :: is_initialized = .false.
         contains
            private
      ! neklab_newton_control (this file)
            procedure, pass(self), public :: init_flow
            procedure, pass(self), public :: amplitude
            procedure, pass(self), public :: forcing
            procedure, pass(self), public :: is_unsteady
            procedure, pass(self), public :: is_initialised
      ! control_gs
            procedure, pass(self), public :: get_dpds
            procedure, pass(self), public :: set_dpds
            procedure, pass(self), public :: get_dpds_helix
            procedure, pass(self), public :: set_dpds_helix
            procedure, pass(self), public :: get_amp_phase
            procedure, pass(self), public :: add_amplitude_step
            procedure, pass(self), public :: probe_dpds
            procedure, pass(self), public :: ensure_nonzero_mean
            procedure, pass(self), public :: seed_harmonics
            procedure, pass(self), public :: rotate_in_time
            procedure, pass(self), public :: get_target
            procedure, pass(self), public :: set_target
            procedure, pass(self), public :: set_target_helix
            procedure, pass(self), public :: get_slope
            procedure, pass(self), public :: set_slope
            procedure, pass(self), public :: has_slope
            procedure, pass(self), public :: get_mflow
            procedure, pass(self), public :: get_mflow_helix
            procedure, pass(self), public :: get_omega
            procedure, pass(self), public :: get_period
            procedure, pass(self), public :: get_womersley
            procedure, pass(self), public :: get_nf
            procedure, pass(self), public :: get_nmf
            procedure, pass(self), public :: get_kharm
            procedure, pass(self), public :: summary
            procedure, pass(self), public :: forcing_summary
            procedure, pass(self), public :: mflow_summary
      ! control_flowrate
            procedure, pass(self), public :: build_area_weights
            procedure, pass(self), public :: get_area
            procedure, pass(self), public :: get_delta
            procedure, pass(self), public :: get_curv_radius
            procedure, pass(self), public :: get_radius
            procedure, pass(self), public :: ubar_arr
            procedure, pass(self), public :: ubar
            procedure, pass(self), public :: reset_mflow
            procedure, pass(self), public :: accumulate_mflow
            procedure, pass(self), public :: close_mflow
            procedure, pass(self), public :: measure_mflow
            procedure, pass(self), public :: seed_jacobian
         end type nek_control

         type(nek_control), public :: ctrl
      !! The one instance. Everything talks to this.

         public :: ctrl_basis
      !! Temporal basis function, exposed for the analysis driver's diagnostics.

      !--------------------------------------------------------------------
      !-----     SUBMODULE INTERFACES                                 -----
      !--------------------------------------------------------------------

      ! --- control_gs
         interface
            module subroutine get_dpds(self, dpds, amp, phase)
               class(nek_control), intent(in) :: self
               real(dp), dimension(lfc), intent(out) :: dpds
               real(dp), dimension(lmfc), optional, intent(out) :: amp
               real(dp), dimension(lmfc), optional, intent(out) :: phase
            end subroutine get_dpds

            module subroutine set_dpds(self, dpds)
               class(nek_control), intent(inout) :: self
               real(dp), dimension(:), intent(in) :: dpds
            end subroutine set_dpds

            module subroutine get_dpds_helix(self, dpds)
               class(nek_control), intent(in) :: self
               real(dp), dimension(:), intent(out) :: dpds
            end subroutine get_dpds_helix

            module subroutine set_dpds_helix(self, dpds)
               class(nek_control), intent(inout) :: self
               real(dp), dimension(:), intent(in) :: dpds
            end subroutine set_dpds_helix

            module subroutine get_amp_phase(self, amp, phase)
               class(nek_control), intent(in) :: self
               real(dp), dimension(lmfc), intent(out) :: amp
               real(dp), dimension(lmfc), intent(out) :: phase
            end subroutine get_amp_phase

            module subroutine add_amplitude_step(self, da)
               class(nek_control), intent(inout) :: self
               real(dp), dimension(:), intent(in) :: da
            end subroutine add_amplitude_step

            module function probe_dpds(self, i, eps) result(d)
               class(nek_control), intent(in) :: self
               integer, intent(in) :: i
               real(dp), intent(in) :: eps
               real(dp), dimension(lfc) :: d
            end function probe_dpds

            module subroutine ensure_nonzero_mean(self, probe)
               class(nek_control), intent(inout) :: self
               real(dp), intent(in) :: probe
            end subroutine ensure_nonzero_mean

            module subroutine seed_harmonics(self, force)
               class(nek_control), intent(inout) :: self
               logical, optional, intent(in) :: force
            end subroutine seed_harmonics

            module subroutine rotate_in_time(self, shift)
               class(nek_control), intent(inout) :: self
               real(dp), intent(in) :: shift
            end subroutine rotate_in_time

            module function get_target(self) result(tgt)
               class(nek_control), intent(in) :: self
               real(dp), dimension(lmfc) :: tgt
            end function get_target

            module subroutine set_target(self, tgt)
               class(nek_control), intent(inout) :: self
               real(dp), dimension(:), intent(in) :: tgt
            end subroutine set_target

            module subroutine set_target_helix(self, tgt)
               class(nek_control), intent(inout) :: self
               real(dp), dimension(:), intent(in) :: tgt
            end subroutine set_target_helix

            module function get_slope(self) result(g)
               class(nek_control), intent(in) :: self
               real(dp) :: g
            end function get_slope

            module subroutine set_slope(self, g)
               class(nek_control), intent(inout) :: self
               real(dp), intent(in) :: g
            end subroutine set_slope

            module function has_slope(self) result(l)
               class(nek_control), intent(in) :: self
               logical :: l
            end function has_slope

            module subroutine get_mflow(self, mflow, amp, phase, qerr)
               class(nek_control), intent(in) :: self
               real(dp), dimension(lfc), intent(out) :: mflow
               real(dp), dimension(lmfc), optional, intent(out) :: amp
               real(dp), dimension(lmfc), optional, intent(out) :: phase
               real(dp), dimension(lmfc), optional, intent(out) :: qerr
            end subroutine get_mflow

            module subroutine get_mflow_helix(self, amp, phase)
               class(nek_control), intent(in) :: self
               real(dp), dimension(:), intent(out) :: amp
               real(dp), dimension(:), optional, intent(out) :: phase
            end subroutine get_mflow_helix

            module function get_omega(self) result(w)
               class(nek_control), intent(in) :: self
               real(dp) :: w
            end function get_omega

            module function get_period(self) result(T)
               class(nek_control), intent(in) :: self
               real(dp) :: T
            end function get_period

            module function get_womersley(self) result(Wo)
               class(nek_control), intent(in) :: self
               real(dp) :: Wo
            end function get_womersley

            module function get_nf(self) result(n)
               class(nek_control), intent(in) :: self
               integer :: n
            end function get_nf

            module function get_nmf(self) result(n)
               class(nek_control), intent(in) :: self
               integer :: n
            end function get_nmf

            module function get_kharm(self) result(k)
               class(nek_control), intent(in) :: self
               integer :: k
            end function get_kharm

            module subroutine summary(self)
               class(nek_control), intent(in) :: self
            end subroutine summary

            module subroutine forcing_summary(self)
               class(nek_control), intent(in) :: self
            end subroutine forcing_summary

            module subroutine mflow_summary(self)
               class(nek_control), intent(in) :: self
            end subroutine mflow_summary
         end interface

      ! --- control_flowrate
         interface
            module subroutine build_area_weights(self, force, radius)
               class(nek_control), intent(inout) :: self
               logical, optional, intent(in) :: force
               real(dp), optional, intent(in) :: radius
            end subroutine build_area_weights

            module function get_area(self) result(a)
               class(nek_control), intent(in) :: self
               real(dp) :: a
            end function get_area

            module function get_delta(self) result(d)
               class(nek_control), intent(in) :: self
               real(dp) :: d
            end function get_delta

            module function get_curv_radius(self) result(r)
               class(nek_control), intent(in) :: self
               real(dp) :: r
            end function get_curv_radius

            module function get_radius(self) result(r)
               class(nek_control), intent(in) :: self
               real(dp) :: r
            end function get_radius

            module function ubar_arr(self, theta) result(Q)
               class(nek_control), intent(in) :: self
               real(dp), dimension(lv), intent(in) :: theta
               real(dp) :: Q
            end function ubar_arr

            module function ubar(self) result(Q)
               class(nek_control), intent(in) :: self
               real(dp) :: Q
            end function ubar

            module subroutine reset_mflow(self, Q0)
               class(nek_control), intent(inout) :: self
               real(dp), intent(in) :: Q0
            end subroutine reset_mflow

            module subroutine accumulate_mflow(self, Q, tval, dtn)
               class(nek_control), intent(inout) :: self
               real(dp), intent(in) :: Q
               real(dp), intent(in) :: tval
               real(dp), intent(in) :: dtn
            end subroutine accumulate_mflow

            module subroutine close_mflow(self, period)
               class(nek_control), intent(inout) :: self
               real(dp), optional, intent(in) :: period
            end subroutine close_mflow

            module subroutine measure_mflow(self, theta, mf, phase, qerr)
               class(nek_control), intent(inout) :: self
               real(dp), dimension(lv), intent(in) :: theta
               real(dp), dimension(:), intent(out) :: mf
               real(dp), dimension(:), optional, intent(out) :: phase
               real(dp), dimension(:), optional, intent(out) :: qerr
            end subroutine measure_mflow

            module subroutine seed_jacobian(self, J, mf)
               class(nek_control), intent(inout) :: self
               real(dp), dimension(:, :), intent(out) :: J
               real(dp), dimension(:), intent(in) :: mf
            end subroutine seed_jacobian
         end interface

      contains

      !====================================================================
      !     CONFIGURATION
      !====================================================================

         subroutine init_flow(self, dpds, womersley, radius)
      !! Configures the control from a forcing vector in the NATIVE convention
      !! and a Womersley number, mirroring neklab_helix % init_flow.
      !!
      !! The number of harmonics follows from the length of dpds: nf = 2K+1, so
      !! a steady deck passes one component and a K = 1 deck passes three. A
      !! zero (or absent) Womersley number selects the steady problem, in which
      !! case only the mean is kept.
      !!
      !!    omega = Wo**2 * nu ,   T = 2*pi/omega ,   nu = cpfld(1,1)
      !!
      !! Wo is fixed for the whole run: the driver calls this once and the
      !! period never changes afterwards.
            class(nek_control), intent(inout) :: self
            real(dp), dimension(:), intent(in) :: dpds
      !! Forcing components, native convention, nf = 2K+1 of them.
            real(dp), optional, intent(in) :: womersley
      !! Womersley number of the fundamental. Zero or absent: steady problem.
            real(dp), optional, intent(in) :: radius
      !! Cross-section radius. Measured from the mesh if absent.
      ! internal
            character(len=*), parameter :: this_procedure = 'init_flow'
            character(len=256) :: msg
            real(dp) :: Wo, pi
            integer :: n

            pi = 4.0_dp*atan(1.0_dp)
            n = size(dpds)
            Wo = optval(womersley, 0.0_dp)

            if (n < 1) then
               call nek_stop_error('init_flow requires at least one forcing component.',
     &            this_module, this_procedure)
            end if
            if (n > lfc) then
               write (msg, '(A,I0,A,I0,A)') 'nf= ', n, ' > lfc= ', lfc, '. Increase kmax_ctrl and recompile.'
               call nek_stop_error(msg, this_module, this_procedure)
            end if
            if (mod(n, 2) == 0) then
               write (msg, '(A,I0,A)') 'nf= ', n, ' is even. nf must be 2K+1.'
               call nek_stop_error(msg, this_module, this_procedure)
            end if

            self%if_unsteady = (abs(Wo) > atol_dp)
            self%womersley = Wo

            if (self%if_unsteady) then
               if (n == 1) then
                  call nek_stop_error('A non-zero Womersley number requires nf > 1 forcing components.',
     &               this_module, this_procedure)
               end if
               self%nf = n
               self%kharm = (n - 1)/2
               self%omega = (Wo**2)*cpfld(1, 1)
               self%period = 2.0_dp*pi/self%omega
               call nek_log_message('Running unsteady (pulsatile) case.', this_module, this_procedure)
            else
               if (n > 1) then
                  call nek_log_warning('Zero Womersley number: harmonic forcing components are ignored.',
     &               this_module, this_procedure)
               end if
               self%nf = 1
               self%kharm = 0
               self%omega = 0.0_dp
               self%period = 0.0_dp
               call nek_log_message('Running steady case.', this_module, this_procedure)
            end if
            self%nmf = self%kharm + 1

            self%dpds = 0.0_dp
            self%dpds(1:self%nf) = dpds(1:self%nf)
            call bcast(self%dpds, lfc*wdsize)

      ! measurement state belongs to a trajectory, not to a configuration
            self%qfour = 0.0_dp
            self%mf = 0.0_dp
            self%mf_phase = 0.0_dp
            self%mf_qerr = 0.0_dp
            self%mf_extracted = .false.
            self%accumulating = .false.

      ! geometry. Forced, because a restart may have re-read the mesh.
            call self%build_area_weights(force=.true., radius=radius)

            self%is_initialized = .true.
            call self%summary()
         end subroutine init_flow

      !====================================================================
      !     FORCING
      !====================================================================

         real(dp) function amplitude(self, tval) result(f)
      !! Instantaneous streamwise forcing amplitude in the native convention.
      !! This is the analogue of neklab_helix % forcing_amplitude and it is the
      !! ONLY place the temporal shape of the forcing is defined.
            class(nek_control), intent(in) :: self
            real(dp), intent(in) :: tval
      ! internal
            integer :: k, i
            real(dp) :: wt
            f = self%dpds(1)
            if (self%if_unsteady) then
               do k = 1, self%kharm
                  i = 2*k
                  wt = k*self%omega*tval
                  f = f + self%dpds(i)*cos(wt) + self%dpds(i + 1)*sin(wt)
               end do
            end if
         end function amplitude

         real(dp) function forcing(self, tval, y) result(f)
      !! Streamwise source term for a single grid point of the 2Dh mesh, where
      !! u_phi is carried by the temperature field. Call from userq as
      !!
      !!    qvol = ctrl%forcing(time, ym1(ix, iy, iz, gllel(ieg)))
      !!
      !! Note that Nek's makeq evaluates userq at t_n (it shifts 'time' by -dt
      !! internally), which is what the linearisation requires. There is no
      !! perturbation branch: the forcing is not a Newton unknown, so the
      !! linearised equations carry no source and userq must return zero for
      !! jp > 0.
      !!
      !! f_phi = amplitude/R is irrotational; a uniform f_phi is not.
            class(nek_control), intent(in) :: self
            real(dp), intent(in) :: tval
            real(dp), intent(in) :: y
            f = self%amplitude(tval)/y
         end function forcing

      !====================================================================
      !     SMALL PREDICATES AND THE TEMPORAL BASIS
      !====================================================================

         logical function is_unsteady(self) result(l)
            class(nek_control), intent(in) :: self
            l = self%if_unsteady
         end function is_unsteady

         logical function is_initialised(self) result(l)
            class(nek_control), intent(in) :: self
            l = self%is_initialized
         end function is_initialised

         real(dp) function ctrl_basis(i, omega, tval) result(phi)
      !! i-th temporal basis function: 1, cos(wt), sin(wt), cos(2wt), ...
      !! Module-level rather than type-bound so that the accumulator can call
      !! it without a self reference in an inner loop.
            integer, intent(in) :: i
            real(dp), intent(in) :: omega
            real(dp), intent(in) :: tval
      ! internal
            integer :: k
            if (i == 1) then
               phi = 1.0_dp
            else
               k = i/2
               if (mod(i, 2) == 0) then
                  phi = cos(k*omega*tval)
               else
                  phi = sin(k*omega*tval)
               end if
            end if
         end function ctrl_basis

      end module neklab_newton_control