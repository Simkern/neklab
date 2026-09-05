      submodule(neklab_t2Dh) t2Dh_flowrate
      !! Cross-section geometry, the flow rate, its running Fourier transform
      !! and the seed for the outer Jacobian.
      !!
      !!--------------------------------------------------------------------
      !! THE CONSTRAINT
      !!--------------------------------------------------------------------
      !!
      !! The measured quantity is the BULK VELOCITY, Q = (1/A) int_A u_phi dA,
      !! not the volume flux. It is O(1), it is what the helix targets are given
      !! in, and using it consistently removes every stray factor of A that the
      !! earlier driver carried around its noise floors and tolerances.
      !!
      !! For the unsteady problem the constraint is the Fourier content of Q
      !! over one period,
      !!
      !!    q_1     = (1/T) int_0^T Q dt
      !!    q_{2k}  = (2/T) int_0^T Q cos(k w t) dt
      !!    q_{2k+1}= (2/T) int_0^T Q sin(k w t) dt
      !!
      !! accumulated along the trajectory. For the steady problem it is simply
      !! the instantaneous bulk velocity of the converged state, read from the
      !! state vector rather than from the Nek field -- the last residual
      !! evaluation of the inner Newton does not necessarily leave the iterate
      !! in vx/vy/t.
      !!
      !!--------------------------------------------------------------------
      !! QUADRATURE ERROR
      !!--------------------------------------------------------------------
      !!
      !! The accumulator is a plain trapezoidal rule on the (CFL-adaptive) time
      !! grid, O(dt**2), and there is no point chasing a flow-rate tolerance
      !! below its own error. That error is MEASURED rather than guessed: a
      !! second accumulator applies the same rule on PAIRS of steps, and
      !! Richardson gives
      !!
      !!    err ~ |fine - coarse| / 3 .
      !!
      !! This replaces the earlier dt**2/100, which was dimensionally a time
      !! squared with an invented constant in front.
      !!
      !! The zero-crossing split that helix applies (helix_mflow_fft.f90:35-49)
      !! is deliberately NOT reproduced. It changes nothing for the mean -- the
      !! trapezoid of a linear function is exact either way -- and for the
      !! harmonics it is an arbitrary partial refinement on the one step where
      !! Q happens to change sign, which is not where the quadrature error
      !! lives. Both of its bugs (a crossing "fraction" that can exceed one, and
      !! a factor 1/2 applied to the mean but not to the harmonics) become moot.
         implicit none
      contains

      !====================================================================
      !     CROSS-SECTION GEOMETRY
      !====================================================================

         module procedure build_area_weights
      !! Mass matrix for integrals over the phi = const cross-section, whose
      !! area element is dA = dR dz. If the mesh runs with ifaxis = .true., bm1
      !! carries an extra radial weight which is stripped here; if the torus
      !! metric lives entirely in the coefficient arrays of neklab_2Dh_axisym on
      !! a plain 2D mesh, bm1 is already correct.
      !!
         character(len=*), parameter :: this_procedure = 'build_area_weights'
         character(len=256) :: msg
         real(dp), external :: glsum, glsc2, glmax, glmin
         real(dp) :: r_area, r_extent, pi
         integer :: n
         integer, parameter :: pad = 10
         if (self%area_defined .and. .not. optval(force, .false.)) return
         pi = 4.0_dp*atan(1.0_dp)
         n = lx1*ly1*lz1*nelv
         if (ifaxis) then
            call invcol3(self%bm_area, bm1, ym1, n)
         else
            call copy(self%bm_area, bm1, n)
         end if
         self%area = glsum(self%bm_area, n)
         if (self%area <= 0.0_dp) then
            call nek_stop_error('Cross-section area is not positive.', this_module, this_procedure)
         end if
         self%area_defined = .true.
         end procedure build_area_weights

         module procedure get_area
         a = self%area
         end procedure get_area

         module procedure get_delta
         d = self%delta
         end procedure get_delta

         module procedure get_curv_radius
         r = self%curv_radius
         end procedure get_curv_radius

         module procedure get_radius
         r = self%radius
         end procedure get_radius

         module procedure get_lambda
         lam = self%lambda
         end procedure get_lambda

         module procedure get_axial_centre
         zc = self%axial_centre
         end procedure get_axial_centre

         module procedure is_torsion_zero
         l = (.not. self%torsion_defined) .or. (self%lambda == 0.0_dp)
         end procedure is_torsion_zero

      !====================================================================
      !     FLOW RATE
      !====================================================================

         module procedure ubar_arr
      !! Bulk velocity of an arbitrary field array, Q = (1/A) int_A theta dA.
         character(len=*), parameter :: this_procedure = 'ubar_arr'
         real(dp), external :: glsc2
         integer :: n
         if (.not. self%area_defined) then
            call nek_stop_error('Area weights are not built. Call t2Dh%init_geom first.',
     &         this_module, this_procedure)
         end if
         n = lx1*ly1*lz1*nelv
         Q = glsc2(theta, self%bm_area, n)/self%area
         end procedure ubar_arr

         module procedure ubar
      !! Bulk velocity of the current Nek field t, which carries u_phi on the
      !! 2Dh mesh.
         Q = self%ubar_arr(t(1, 1, 1, 1, 1))
         end procedure ubar

      !====================================================================
      !     RUNNING FOURIER ACCUMULATOR
      !====================================================================

         module procedure reset_mflow
      !! Starts an accumulation. Q0 is the flow rate at t = 0, i.e. that of the
      !! initial condition, before the first step is taken: the trapezoidal rule
      !! needs a left endpoint for the first interval.
         self%q_acc = 0.0_dp
         self%q_crs = 0.0_dp
         self%t_acc = 0.0_dp
         self%q_lag = Q0
         self%t_lag = 0.0_dp
         self%q_lag2 = Q0
         self%t_lag2 = 0.0_dp
         self%nacc = 0
         self%accumulating = .true.
         self%mf_extracted = .false.
         end procedure reset_mflow

         module procedure accumulate_mflow
      !! Call once per timestep, AFTER nek_advance, with the flow rate of the
      !! current field, the current time and the timestep just taken. tval is
      !! the time at the END of the step.
         character(len=*), parameter :: this_procedure = 'accumulate_mflow'
         integer :: i
         real(dp) :: dtc
         if (.not. self%accumulating) then
            call nek_stop_error('Accumulator is not open. Call t2Dh%reset_mflow first.',
     &         this_module, this_procedure)
         end if
         self%t_acc = self%t_acc + dtn
         self%nacc = self%nacc + 1
      ! fine rule: trapezoid on [t_lag, tval]
         do i = 1, self%nf
            self%q_acc(i) = self%q_acc(i) + 0.5_dp*dtn*
     &         (self%q_lag*t2Dh_basis(i, self%omega, self%t_lag) + Q*t2Dh_basis(i, self%omega, tval))
         end do
      ! coarse rule: trapezoid on [t_lag2, tval], closed every second step
         if (mod(self%nacc, 2) == 0) then
            dtc = tval - self%t_lag2
            do i = 1, self%nf
               self%q_crs(i) = self%q_crs(i) + 0.5_dp*dtc*
     &            (self%q_lag2*t2Dh_basis(i, self%omega, self%t_lag2) + Q*t2Dh_basis(i, self%omega, tval))
            end do
            self%q_lag2 = Q
            self%t_lag2 = tval
         end if
         self%q_lag = Q
         self%t_lag = tval
         end procedure accumulate_mflow

         module procedure close_mflow
      !! Normalises both accumulators, converts to amplitudes and phases and
      !! forms the Richardson error estimate. Call once at the end of the
      !! nonlinear pass, after the horizon has been reached.
         character(len=*), parameter :: this_procedure = 'close_mflow'
         character(len=256) :: msg
         real(dp), dimension(lfc) :: cfine, ccrs, cerr
         real(dp) :: T, dtc
         integer :: i, k
         if (.not. self%accumulating) then
            call nek_stop_error('Accumulator is not open.', this_module, this_procedure)
         end if
         if (self%t_acc <= 0.0_dp) then
            call nek_stop_error('Empty accumulator.', this_module, this_procedure)
         end if
         T = optval(period, self%t_acc)
         if (present(period)) then
            if (abs(self%t_acc - period) > 1.0e-08_dp*period) then
               write (msg, '(A,E16.8,A,E16.8)') 'Integration time ', self%t_acc,
     &            ' does not match the period ', period
               call nek_stop_error(msg, this_module, this_procedure)
            end if
         end if
      ! close an incomplete pair so that both rules span the same interval
         if (mod(self%nacc, 2) == 1) then
            dtc = self%t_lag - self%t_lag2
            do i = 1, self%nf
               self%q_crs(i) = self%q_crs(i) + 0.5_dp*dtc*
     &            (self%q_lag2*t2Dh_basis(i, self%omega, self%t_lag2)
     &             + self%q_lag*t2Dh_basis(i, self%omega, self%t_lag))
            end do
         end if
      ! normalise: mean, then twice the projection for each harmonic
         cfine = 0.0_dp; ccrs = 0.0_dp
         cfine(1) = self%q_acc(1)/self%t_acc
         ccrs(1) = self%q_crs(1)/self%t_acc
         do i = 2, self%nf
            cfine(i) = 2.0_dp*self%q_acc(i)/self%t_acc
            ccrs(i) = 2.0_dp*self%q_crs(i)/self%t_acc
         end do
         self%qfour = cfine
      ! Richardson estimate of the trapezoidal error
         cerr = abs(cfine - ccrs)/3.0_dp
         self%mf_qerr = 0.0_dp
         self%mf = 0.0_dp
         self%mf_phase = 0.0_dp
         self%mf(1) = cfine(1)
         self%mf_qerr(1) = cerr(1)
         do k = 1, self%kharm
            i = 2*k
            self%mf(k + 1) = hypot(cfine(i), cfine(i + 1))
            self%mf_phase(k + 1) = atan2(cfine(i + 1), cfine(i))
            self%mf_qerr(k + 1) = hypot(cerr(i), cerr(i + 1))
         end do
         if (self%nacc < 4) then
            self%mf_qerr = 0.0_dp
            call nek_log_warning('Fewer than four steps accumulated: no quadrature error estimate.',
     &         this_module, this_procedure)
         end if
         self%accumulating = .false.
         self%mf_extracted = .true.
         end procedure close_mflow

         module procedure measure_mflow
      !! The single measurement entry point. Dispatches on the regime so that
      !! the outer Newton never has to know which problem it is solving:
      !!
      !!   unsteady : the Fourier amplitudes accumulated along the last orbit,
      !!              which the nonlinear map closed with close_mflow;
      !!   steady   : the instantaneous bulk velocity of the state vector.
      !!
      !! theta is the temperature block of the Newton iterate. It is used in the
      !! steady branch only; in the unsteady branch the measurement belongs to
      !! the trajectory, not to its endpoint, and theta is ignored.
         character(len=*), parameter :: this_procedure = 'measure_mflow'
         integer :: n
         n = min(size(mf), self%nmf)
         mf = 0.0_dp
         if (present(phase)) phase = 0.0_dp
         if (present(qerr)) qerr = 0.0_dp
         if (self%if_unsteady) then
            if (.not. self%mf_extracted) then
               call nek_stop_error('No flow-rate measurement available: the nonlinear map must '//
     &            'call t2Dh%close_mflow before the driver measures.', this_module, this_procedure)
            end if
            mf(1:n) = self%mf(1:n)
            if (present(phase)) phase(1:min(size(phase), self%nmf)) = self%mf_phase(1:min(size(phase), self%nmf))
            if (present(qerr)) qerr(1:min(size(qerr), self%nmf)) = self%mf_qerr(1:min(size(qerr), self%nmf))
         else
            mf(1) = self%ubar_arr(theta)
      ! cache it, so that t2Dh holds the last measurement in both regimes and
      ! mflow_summary has something to print
            self%qfour = 0.0_dp
            self%mf = 0.0_dp
            self%mf_phase = 0.0_dp
            self%mf_qerr = 0.0_dp
            self%qfour(1) = mf(1)
            self%mf(1) = mf(1)
            self%mf_extracted = .true.
         end if
         end procedure measure_mflow

      !====================================================================
      !     JACOBIAN SEED
      !====================================================================

         module procedure seed_jacobian
      !! Diagonal seed for the outer Jacobian J_ij = d(mf_i)/d(a_j), from a
      !! lumped RL model of each harmonic:
      !!
      !!    resistance  1/G      G = dQ_0/da_0, measured by the steady solve
      !!    inductance  c k w/d  the inertial law, f_phi = a/R and R_c du/dt = f
      !!
      !! Impedances add in quadrature, so
      !!
      !!    J_kk = 1 / sqrt( (1/G)**2 + (c k w/delta)**2 )
      !!
      !! which is exact in both limits and monotone between them. k = 0 gives
      !! J_11 = G, so the steady case is the same formula. The true Womersley
      !! admittance would need Bessel functions for a correction that Broyden
      !! removes within one or two outer steps.
      !!
      !! If no steady resistance is on record, the Stokes estimate G = Q_0/a_0
      !! is taken from the current measurement. For creeping flow that is exact;
      !! with inertia it OVERestimates dQ/da, so the first Newton step is a safe
      !! undershoot and no extra nonlinear solve is needed to get started.
         character(len=*), parameter :: this_procedure = 'seed_jacobian'
         character(len=256) :: msg
         real(dp) :: G, wk
         integer :: k, n
         n = self%nmf
         J = 0.0_dp
         if (.not. self%gslope_defined) then
            if (abs(self%dpds(1)) <= atol_dp) then
               call nek_stop_error('Cannot infer the resistance from a zero mean forcing.',
     &            this_module, this_procedure)
            end if
            G = mf(1)/self%dpds(1)
            if (G <= 0.0_dp) then
               G = abs(mf(1)/self%dpds(1))
               call nek_log_warning('Non-positive Stokes slope. Using its magnitude.',
     &            this_module, this_procedure)
            end if
            if (G <= atol_dp) then
               call nek_stop_error('Degenerate Stokes slope: the flow rate does not respond to the '//
     &            'mean forcing.', this_module, this_procedure)
            end if
            call self%set_slope(G)
            write (msg, '(A,E16.8)') 'No steady resistance on record. Stokes estimate dQ/da= ', G
            call nek_log_message(msg, this_module, this_procedure)
         end if
         G = self%gslope
         J(1, 1) = G
         do k = 1, min(self%kharm, n - 1)
            wk = c_inertial*k*self%omega/self%delta
            J(k + 1, k + 1) = 1.0_dp/sqrt((1.0_dp/G)**2 + wk**2)
         end do
         write (msg, '(A,*(1X,E16.8))') 'Seeded outer Jacobian diagonal:', (J(k, k), k=1, n)
         call nek_log_message(msg, this_module, this_procedure)
         end procedure seed_jacobian

      end submodule t2Dh_flowrate