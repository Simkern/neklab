      submodule(neklab_newton_control) control_gs
      !! Getters, setters, convention converters and summaries for nek_control.
      !!
      !! The helix converters here are the ONLY place the old (helix) Fourier
      !! convention appears. Everything else in neklab is native:
      !!
      !!   forcing    a_0 = d_0 ,  a_ck =  2 d_ck ,  a_sk = -2 d_sk
      !!   flow rate  q_0 = m_0 ,  q_ck =  2 m_ck ,  q_sk =  2 m_sk
      !!
      !! so the two forcing amplitudes agree while a helix flow-rate amplitude
      !! is HALF the native one. The sign flip on the forcing sine is what makes
      !! helix's forcing and flow-rate phases rotate in opposite directions
      !! under a time shift; natively they rotate together.
         implicit none
      contains

      !====================================================================
      !     FORCING
      !====================================================================

         module procedure get_dpds
         integer :: k, i
         dpds = self%dpds
         if (present(amp) .or. present(phase)) then
            if (present(amp)) then
               amp = 0.0_dp
               amp(1) = self%dpds(1)
            end if
            if (present(phase)) phase = 0.0_dp
            do k = 1, self%kharm
               i = 2*k
               if (present(amp)) amp(k + 1) = hypot(self%dpds(i), self%dpds(i + 1))
               if (present(phase)) phase(k + 1) = atan2(self%dpds(i + 1), self%dpds(i))
            end do
         end if
         end procedure get_dpds

         module procedure set_dpds
         integer :: n
         n = min(size(dpds), self%nf)
         self%dpds = 0.0_dp
         self%dpds(1:n) = dpds(1:n)
         call bcast(self%dpds, lfc*wdsize)
         end procedure set_dpds

         module procedure get_dpds_helix
         integer :: k, i, n
         n = min(size(dpds), self%nf)
         dpds = 0.0_dp
         dpds(1) = self%dpds(1)
         do k = 1, self%kharm
            i = 2*k
            if (i + 1 > n) exit
            dpds(i) = 0.5_dp*self%dpds(i)
            dpds(i + 1) = -0.5_dp*self%dpds(i + 1)
         end do
         end procedure get_dpds_helix

         module procedure set_dpds_helix
         character(len=*), parameter :: this_procedure = 'set_dpds_helix'
         character(len=128) :: msg
         real(dp), dimension(lfc) :: d
         integer :: k, i, n
         n = min(size(dpds), self%nf)
         if (n < self%nf) then
            write (msg, '(A,I0,A,I0,A)') 'Only ', n, ' of nf= ', self%nf, ' components supplied.'
            call nek_log_warning(msg, this_module, this_procedure)
         end if
         d = 0.0_dp
         d(1) = dpds(1)
         do k = 1, self%kharm
            i = 2*k
            if (i + 1 > n) exit
            d(i) = 2.0_dp*dpds(i)
            d(i + 1) = -2.0_dp*dpds(i + 1)
         end do
         call self%set_dpds(d)
         end procedure set_dpds_helix

         module procedure get_amp_phase
         real(dp), dimension(lfc) :: d
         call self%get_dpds(d, amp, phase)
         end procedure get_amp_phase

         module procedure add_amplitude_step
      !! Adds da(i) to the amplitude of component i along its CURRENT (frozen)
      !! phase direction. The phase of a harmonic is therefore never changed by
      !! the outer Newton -- only its amplitude is an unknown. A harmonic that
      !! is exactly zero has no phase to freeze, which is a configuration error
      !! rather than something to paper over: call seed_harmonics first.
         character(len=*), parameter :: this_procedure = 'add_amplitude_step'
         character(len=256) :: msg
         real(dp), dimension(lmfc) :: amp, phase
         integer :: k, i, n
         n = min(size(da), self%nmf)
         call self%get_amp_phase(amp, phase)
         do k = 1, self%kharm
            if (k + 1 > n) exit
            if (amp(k + 1) <= atol_dp .and. abs(da(k + 1)) > atol_dp) then
               write (msg, '(A,I0,A)') 'Harmonic ', k, ' is zero: its phase direction is undefined. '//
     &            'Call ctrl%seed_harmonics before stepping the amplitudes.'
               call nek_stop_error(msg, this_module, this_procedure)
            end if
         end do
         self%dpds(1) = self%dpds(1) + da(1)
         do k = 1, self%kharm
            if (k + 1 > n) exit
            i = 2*k
            self%dpds(i) = self%dpds(i) + cos(phase(k + 1))*da(k + 1)
            self%dpds(i + 1) = self%dpds(i + 1) + sin(phase(k + 1))*da(k + 1)
         end do
         call bcast(self%dpds, lfc*wdsize)
         end procedure add_amplitude_step

         module procedure probe_dpds
      !! Forcing vector with amplitude i displaced by eps along its frozen
      !! phase direction. Used to build the finite-difference columns of the
      !! outer Jacobian without disturbing the stored state.
         real(dp), dimension(lmfc) :: amp, phase
         integer :: j
         d = self%dpds
         if (i < 1 .or. i > self%nmf) return
         if (i == 1) then
            d(1) = d(1) + eps
         else
            call self%get_amp_phase(amp, phase)
            j = 2*(i - 1)
            d(j) = d(j) + cos(phase(i))*eps
            d(j + 1) = d(j + 1) + sin(phase(i))*eps
         end if
         end procedure probe_dpds

         module procedure ensure_nonzero_mean
      !! A zero mean forcing gives a zero mean flow rate, hence no slope to
      !! infer: displace it to a probe value before the first solve.
         character(len=*), parameter :: this_procedure = 'ensure_nonzero_mean'
         character(len=256) :: msg
         if (abs(self%dpds(1)) > atol_dp) return
         self%dpds(1) = probe
         call bcast(self%dpds, lfc*wdsize)
         write (msg, '(A,E16.8)') 'Mean forcing is zero. Using probe value dpds(1)= ', probe
         call nek_log_warning(msg, this_module, this_procedure)
         end procedure ensure_nonzero_mean

         module procedure seed_harmonics
      !! Fills harmonics that are exactly zero -- the state a cold start from a
      !! steady solve arrives in -- using the lumped RL model of seed_jacobian:
      !!
      !!    |a_k|  = target_k / J_kk ,   J_kk = 1/sqrt( (1/G)**2 + (c k w/d)**2 )
      !!    phase  = atan( c k w G / d )
      !!
      !! The phase matters more than the amplitude here. A zero harmonic has no
      !! phase, atan2(0,0) returns zero, and the frozen-phase machinery would
      !! then probe along an arbitrary direction and build a meaningless
      !! Jacobian column. In the inertial limit the forcing leads the flow by
      !! exactly pi/2, and the gauge puts the fundamental flow-rate phase at
      !! zero, so the expression above is the forcing phase measured from that
      !! gauge. It interpolates to zero lead in the resistive limit.
         character(len=*), parameter :: this_procedure = 'seed_harmonics'
         character(len=256) :: msg
         real(dp), dimension(lmfc) :: amp, phase, tgt
         real(dp) :: G, wk, xk, ak, pk
         integer :: k, i
         logical :: force_, seeded
         force_ = optval(force, .false.)
         if (.not. self%if_unsteady) return
         if (self%kharm < 1) return
         if (.not. self%target_defined) then
            call nek_stop_error('seed_harmonics needs the flow-rate targets. Call set_target first.',
     &         this_module, this_procedure)
         end if
         if (.not. self%gslope_defined) then
      ! No steady run on record. The mean forcing and the mean target are both
      ! known and the steady problem is nearly linear, so Q_target/a_0 is a
      ! serviceable estimate of the resistance -- it is exactly the Stokes
      ! estimate evaluated at the target rather than at the current state. Only
      ! the seed depends on it, and Broyden removes the error within a step or
      ! two, but running the steady solve first is still the better path.
            if (abs(self%dpds(1)) <= atol_dp) then
               call nek_stop_error('seed_harmonics needs either the steady resistance dQ/da '//
     &            '(run the steady solve first, or set it with ctrl%set_slope) or a non-zero '//
     &            'mean forcing to estimate it from.', this_module, this_procedure)
            end if
            call self%set_slope(abs(self%mf_target(1)/self%dpds(1)))
            write (msg, '(A,E16.8)') 'No steady resistance on record. Estimating dQ/da= ', self%gslope
            call nek_log_warning(msg, this_module, this_procedure)
         end if
         if (self%delta <= atol_dp) then
            call nek_stop_error('Curvature ratio delta is not set.', this_module, this_procedure)
         end if
         G = self%gslope
         tgt = self%mf_target
         call self%get_amp_phase(amp, phase)
         seeded = .false.
         do k = 1, self%kharm
            i = 2*k
            if (amp(k + 1) > atol_dp .and. .not. force_) cycle
            wk = c_inertial*k*self%omega/self%delta
            xk = sqrt((1.0_dp/G)**2 + wk**2)
            ak = tgt(k + 1)*xk
            pk = atan(wk*G)
            self%dpds(i) = ak*cos(pk)
            self%dpds(i + 1) = ak*sin(pk)
            seeded = .true.
            write (msg, '(A,I0,A,E16.8,A,E16.8)') 'Seeded harmonic ', k, ': |a|= ', ak, ', phase= ', pk
            call nek_log_message(msg, this_module, this_procedure)
         end do
         if (seeded) then
            call bcast(self%dpds, lfc*wdsize)
            if (self%kharm >= 2) then
               call nek_log_warning('The inertial seed is calibrated at k = 1 only; its linear '//
     &            'extrapolation to k >= 2 is untested.', this_module, this_procedure)
            end if
            call self%forcing_summary()
         end if
         end procedure seed_harmonics

         module procedure rotate_in_time
      !! Re-references the forcing to a time origin displaced by `shift`, i.e.
      !! replaces f(t) by f(t + shift). Under that map phase_k -> phase_k - k w s
      !! for BOTH the forcing and the flow rate, which is the whole point of the
      !! native convention: the same rotation applies to both arrays.
         integer :: k, i
         real(dp) :: ang, c, s, ca, sa
         if (.not. self%if_unsteady) return
         do k = 1, self%kharm
            i = 2*k
            ang = k*self%omega*shift
            c = self%dpds(i); s = self%dpds(i + 1)
            ca = cos(ang); sa = sin(ang)
            self%dpds(i) = c*ca + s*sa
            self%dpds(i + 1) = -c*sa + s*ca
         end do
         call bcast(self%dpds, lfc*wdsize)
         end procedure rotate_in_time

      !====================================================================
      !     TARGETS AND THE STEADY RESISTANCE
      !====================================================================

         module procedure get_target
         tgt = self%mf_target
         end procedure get_target

         module procedure set_target
         integer :: n
         n = min(size(tgt), self%nmf)
         self%mf_target = 0.0_dp
         self%mf_target(1:n) = tgt(1:n)
         self%target_defined = .true.
         end procedure set_target

         module procedure set_target_helix
      !! Helix reports the flow-rate amplitude as HALF the peak excursion, so a
      !! target taken from a helix deck doubles on the way in. The mean does not.
         real(dp), dimension(lmfc) :: t_
         integer :: n, i
         n = min(size(tgt), self%nmf)
         t_ = 0.0_dp
         t_(1) = tgt(1)
         do i = 2, n
            t_(i) = 2.0_dp*tgt(i)
         end do
         call self%set_target(t_)
         end procedure set_target_helix

         module procedure get_slope
         g = self%gslope
         end procedure get_slope

         module procedure set_slope
         character(len=*), parameter :: this_procedure = 'set_slope'
         character(len=128) :: msg
         if (g <= 0.0_dp) then
            write (msg, '(A,E16.8,A)') 'Non-positive resistance dQ/da= ', g, '. Ignored.'
            call nek_log_warning(msg, this_module, this_procedure)
            return
         end if
         self%gslope = g
         self%gslope_defined = .true.
         end procedure set_slope

         module procedure has_slope
         l = self%gslope_defined
         end procedure has_slope

      !====================================================================
      !     MEASURED FLOW RATE
      !====================================================================

         module procedure get_mflow
         mflow = self%qfour
         if (present(amp)) amp = self%mf
         if (present(phase)) phase = self%mf_phase
         if (present(qerr)) qerr = self%mf_qerr
         end procedure get_mflow

         module procedure get_mflow_helix
         integer :: n, i
         n = min(size(amp), self%nmf)
         amp = 0.0_dp
         amp(1) = self%mf(1)
         do i = 2, n
            amp(i) = 0.5_dp*self%mf(i)
         end do
         if (present(phase)) then
            phase = 0.0_dp
            phase(1:min(size(phase), self%nmf)) = self%mf_phase(1:min(size(phase), self%nmf))
         end if
         end procedure get_mflow_helix

      !====================================================================
      !     SCALARS
      !====================================================================

         module procedure get_omega
         w = self%omega
         end procedure get_omega

         module procedure get_period
         T = self%period
         end procedure get_period

         module procedure get_womersley
         Wo = self%womersley
         end procedure get_womersley

         module procedure get_nf
         n = self%nf
         end procedure get_nf

         module procedure get_nmf
         n = self%nmf
         end procedure get_nmf

         module procedure get_kharm
         k = self%kharm
         end procedure get_kharm

      !====================================================================
      !     SUMMARIES
      !====================================================================

         module procedure summary
         character(len=*), parameter :: this_procedure = 'summary'
         character(len=256) :: msg
         integer :: i
         call nek_log_message('Newton control configuration:', this_module, this_procedure)
         write (msg, '(3X,A,L1,A,I0,A,I0,A,I0)') 'unsteady= ', self%if_unsteady,
     &      ', K= ', self%kharm, ', nf= ', self%nf, ', nmf= ', self%nmf
         call nek_log_message(msg, this_module, this_procedure)
         if (self%if_unsteady) then
            write (msg, '(3X,A,E16.8,A,E16.8,A,E16.8)') 'Wo= ', self%womersley,
     &         ', omega= ', self%omega, ', T= ', self%period
            call nek_log_message(msg, this_module, this_procedure)
         end if
         write (msg, '(3X,A,E16.8,A,E16.8,A,E16.8,A,E16.8)') 'area= ', self%area,
     &      ', R_c= ', self%curv_radius, ', r= ', self%radius, ', delta= ', self%delta
         call nek_log_message(msg, this_module, this_procedure)
         if (self%gslope_defined) then
            write (msg, '(3X,A,E16.8)') 'dQ/da (steady resistance)= ', self%gslope
            call nek_log_message(msg, this_module, this_procedure)
         end if
         if (self%target_defined) then
            write (msg, '(3X,A,*(1X,F16.10))') 'target  =', (self%mf_target(i), i=1, self%nmf)
            call nek_log_message(msg, this_module, this_procedure)
         end if
         call self%forcing_summary()
         end procedure summary

         module procedure forcing_summary
         character(len=*), parameter :: this_procedure = 'forcing_summary'
         character(len=512) :: msg
         real(dp), dimension(lfc) :: d
         real(dp), dimension(lmfc) :: amp, phase
         integer :: i
         call self%get_dpds(d, amp, phase)
         write (msg, '(3X,A,*(1X,F16.10))') 'dpds     =', (d(i), i=1, self%nf)
         call nek_log_message(msg, this_module, this_procedure)
         if (self%if_unsteady) then
            write (msg, '(3X,A,*(1X,F16.10))') 'amplitude=', (amp(i), i=1, self%nmf)
            call nek_log_message(msg, this_module, this_procedure)
            write (msg, '(3X,A,*(1X,F16.10))') 'phase    =', (phase(i), i=1, self%nmf)
            call nek_log_message(msg, this_module, this_procedure)
         end if
         end procedure forcing_summary

         module procedure mflow_summary
         character(len=*), parameter :: this_procedure = 'mflow_summary'
         character(len=512) :: msg
         integer :: i
         write (msg, '(3X,A,*(1X,F16.10))') 'mflow    =', (self%mf(i), i=1, self%nmf)
         call nek_log_message(msg, this_module, this_procedure)
         if (self%if_unsteady) then
            write (msg, '(3X,A,*(1X,F16.10))') 'phase    =', (self%mf_phase(i), i=1, self%nmf)
            call nek_log_message(msg, this_module, this_procedure)
            write (msg, '(3X,A,*(1X,E16.8))') 'quad err =', (self%mf_qerr(i), i=1, self%nmf)
            call nek_log_message(msg, this_module, this_procedure)
         end if
         end procedure mflow_summary

      end submodule control_gs
