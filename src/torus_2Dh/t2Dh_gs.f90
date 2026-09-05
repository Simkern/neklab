      submodule(neklab_t2Dh) t2Dh_gs
      !! Getters, setters, convention converters and summaries for nek_t2Dh.
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
     &            'Call t2Dh%seed_harmonics before stepping the amplitudes.'
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
     &            '(run the steady solve first, or set it with t2Dh%set_slope) or a non-zero '//
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

         module procedure set_target_ratio
      !! Targets in the input convention:
      !!
      !!    tgt(1)   = Q(0)          mean bulk velocity, ABSOLUTE
      !!    tgt(k+1) = 2*Q(k)/Q(0)   amplitude RATIO of harmonic k
      !!
      !! The native amplitude is the peak excursion, which is exactly 2*Q(k),
      !! so the harmonic conversion is a multiplication by the mean. At the
      !! usual Q(0) = 1 this is the identity and the input passes through
      !! unchanged; it only bites if the mean target is ever moved off one,
      !! which is precisely why the ratio is stored rather than assumed away.
      !!
      !! NOTE the factor 2: a helix deck reports Q(k) itself (half the peak
      !! excursion), so a helix number must be doubled before it is handed to
      !! this routine.
         character(len=*), parameter :: this_procedure = 'set_target_ratio'
         real(dp), dimension(lmfc) :: t_
         integer :: n, i
         n = min(size(tgt), self%nmf)
         if (abs(tgt(1)) <= atol_dp) then
            call nek_stop_error('The mean flow-rate target is zero: the harmonic targets are '//
     &         'ratios to it and cannot be resolved.', this_module, this_procedure)
         end if
         t_ = 0.0_dp
         t_(1) = tgt(1)
         do i = 2, n
            t_(i) = tgt(i)*abs(tgt(1))
         end do
         call self%set_target(t_)
         end procedure set_target_ratio

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

         module procedure get_mflow_ratio
      !! The measured flow rate in the same convention as set_target_ratio:
      !! index 1 is the absolute mean, index k+1 is 2*Q(k)/Q(0). Directly
      !! comparable with what the deck asked for.
         integer :: n, i
         n = min(size(amp), self%nmf)
         amp = 0.0_dp
         amp(1) = self%mf(1)
         if (abs(self%mf(1)) > atol_dp) then
            do i = 2, n
               amp(i) = self%mf(i)/abs(self%mf(1))
            end do
         end if
         if (present(phase)) then
            phase = 0.0_dp
            phase(1:min(size(phase), self%nmf)) = self%mf_phase(1:min(size(phase), self%nmf))
         end if
         end procedure get_mflow_ratio

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
         integer, parameter :: pad = 16
         call nek_log_message('', this_module, this_procedure)
         call nek_log_message('Newton configuration:', this_module, this_procedure)
         call nek_log_message('', this_module, this_procedure)
         write (msg, '(3X,A,1X,L16)') padl('unsteady:',pad), self%if_unsteady
         call nek_log_message(msg, this_module, this_procedure)
         call nek_log_message('Geometry:', this_module, this_procedure)
         write (msg, '(3X,A,1X,E16.8)') padl('area:',pad), self%area
         call nek_log_message(msg, this_module, this_procedure)
         write (msg, '(3X,A,1X,E16.8)') padl('R_c:',pad), self%curv_radius
         call nek_log_message(msg, this_module, this_procedure)
         write (msg, '(3X,A,1X,E16.8)') padl('r:',pad), self%radius
         call nek_log_message(msg, this_module, this_procedure)
         write (msg, '(3X,A,1X,E16.8)') padl('delta:',pad), self%delta
         call nek_log_message(msg, this_module, this_procedure)
         if (self%torsion_defined) then
            write (msg, '(3X,A,1X,E16.8)') padl('lambda:',pad), self%lambda
            call nek_log_message(msg, this_module, this_procedure)
            write (msg, '(3X,A,1X,E16.8)') padl('axial centre:',pad), self%axial_centre
            call nek_log_message(msg, this_module, this_procedure)
         end if
         call nek_log_message('Dynamics:', this_module, this_procedure)
         write (msg, '(3X,A,2X,3(1X,I4))') padl('K, nf, nmf:',pad), self%kharm, self%nf, self%nmf
         call nek_log_message(msg, this_module, this_procedure)
         if (self%if_unsteady) then
            write (msg, '(3X,A,1X,E16.8)') padl('Wo:',pad), self%womersley
            call nek_log_message(msg, this_module, this_procedure)
            write (msg, '(3X,A,1X,E16.8)') padl('omega:',pad), self%omega
            call nek_log_message(msg, this_module, this_procedure)
            write (msg, '(3X,A,1X,E16.8)') padl('T:',pad), self%period
            call nek_log_message(msg, this_module, this_procedure)
         end if
         if (self%gslope_defined) then
            write (msg, '(3X,A,1X,E16.8)') padl('dQ/da:',pad), self%gslope
            call nek_log_message(msg, this_module, this_procedure)
         end if
         if (self%target_defined) then
            write (msg, '(3X,A,*(1X,F16.10))') padl('Q_target:',pad), (self%mf_target(i), i=1, self%nmf)
            call nek_log_message(msg, this_module, this_procedure)
         end if
         call nek_log_message('', this_module, this_procedure)
         call self%forcing_summary()
         end procedure summary

         module procedure forcing_summary
         character(len=*), parameter :: this_procedure = 'forcing_summary'
         character(len=512) :: msg
         real(dp), dimension(lfc) :: d
         real(dp), dimension(lmfc) :: amp, phase
         integer :: i
         call self%get_dpds(d, amp, phase)
         write (msg, '(3X,A,*(1X,F16.10))')    'dpds      =', (d(i), i=1, self%nf)
         call nek_log_message(msg, this_module, this_procedure)
         if (self%if_unsteady) then
            write (msg, '(3X,A,*(1X,F16.10))') 'amplitude =', (amp(i), i=1, self%nmf)
            call nek_log_message(msg, this_module, this_procedure)
            write (msg, '(3X,A,*(1X,F16.10))') 'phase     =', (phase(i), i=1, self%nmf)
            call nek_log_message(msg, this_module, this_procedure)
         end if
         end procedure forcing_summary

         module procedure mflow_summary
         character(len=*), parameter :: this_procedure = 'mflow_summary'
         character(len=512) :: msg
         integer :: i
         write (msg, '(5X,A,*(1X,F16.10))')    'mflow     =', (self%mf(i), i=1, self%nmf)
         call nek_log_message(msg, this_module, this_procedure)
         if (self%if_unsteady) then
            write (msg, '(5X,A,*(1X,F16.10))') 'phase     =', (self%mf_phase(i), i=1, self%nmf)
            call nek_log_message(msg, this_module, this_procedure)
            write (msg, '(5X,A,*(1X,E16.8))')  'quad err  =', (self%mf_qerr(i), i=1, self%nmf)
            call nek_log_message(msg, this_module, this_procedure)
         end if
         end procedure mflow_summary

      end submodule t2Dh_gs