      module neklab_2Dh_axisym
         use LightKrylov, only: dp
         use neklab_2Dh_common, only: wmask, build_wmask,
     &                                check_2Dh_properties, check_2Dh_masks,
     &                                bm2_map12_2Dh, bm2_map12T_2Dh,
     &                                solve_helmholtz_2Dh
         use neklab_2Dh_torsion, only: if_torsion_zero, check_torsion_cache,
     &                                 build_torsion_coeffs,
     &                                 dphi_apply, dphi_apply_t,
     &                                 dphi_alpha_apply, dphi_alpha_apply_t,
     &                                 add_torsion_forcing_RZ, add_torsion_forcing_phi,
     &                                 add_torsion_coupling
         use neklab_nek_setup, only: nek_log_debug, nek_log_message, nek_stop_error
         use neklab_timing, only: neklab_timer_start, neklab_timer_stop,
     &                            t_axisym_makefp, t_axisym_solve_cpl, t_axisym_solve_pres_cpl,
     &                            t_axisym_pres_prec_cpl, t_axisym_helm_matvec, t_axisym_cpl_matvec,
     &                            t_axisym_pres_matvec
         use neklab_nek_forcing, only: add_neklab_forcing_term, add_neklab_source_term
         implicit none
         include "SIZE"
         include "TOTAL"
         include "ADJOINT" !ifadj
         
         private
         character(len=*), parameter, private :: this_module = 'neklab_2Dh_axisym'
         integer, parameter :: lv = lx1*ly1*lz1*lelv
         integer, parameter :: lp = lx2*ly2*lz2*lelv
         integer, parameter, public :: maxit_p_axisym = 1000 ! max iterations for pressure solve
         
         real(dp), dimension(:,:,:,:), pointer, private :: shift_op => null() ! shift seen by helm_matvec_axisym
         
         real(dp), dimension(lx1,ly1,lz1,lelv), public, target :: h2z_shift
         real(dp), dimension(lx1,ly1,lz1,lelv), public, target :: diag_shift
         real(dp), dimension(lx1,ly1,lz1,lelv), public :: couple_coef
      ! diag_shift is S_PHI = (alpha^2+1)/R^2 and diag_shift_r the u_R shift.
      ! They are equal: the torsion shifts lambda^2/R^2 (Z, R) are explicit
      ! (route A, neklab_2Dh_torsion). diag_shift_r is kept so that the
      ! call sites need not change.
         real(dp), dimension(lx1,ly1,lz1,lelv), public, target :: diag_shift_r
         real(dp), dimension(lx1,ly1,lz1,lelv), public  :: alphaR_coef   ! alpha/R, continuity/pressure coupling
         ! Adjoint: base-flow gradients for (grad U)^T v (prepare_step, ifadj only)
         real(dp), dimension(lv),               private :: gUZ_Z, gUZ_R, gUR_Z, gUR_R
         logical,                               public  :: axisym_coeffs_defined = .false.
         logical,                               public  :: alpha_eq_zero = .false.
         real(dp),                              private :: alpha_cached = huge(1.0_dp)
         logical,                               private :: init_2Dh_axisym_done = .false.

      ! Fields produced and consumed inside the model
         real(dp), dimension(lv,lpert), private :: advZ, advR, advPhi ! slot-swapping explicit terms
         real(dp), dimension(lv),       private :: tmp_m1a, tmp_m1b   ! ophx/ophinv, make_coupling, makefp, incomprp_div/dw, pressure matvecs
         real(dp), dimension(lv),       private :: tmp_m1c            ! helmholtz_matvec_2Dh_axisym and coupled_helmholtz_matvec_2Dh_axisym 
         real(dp), dimension(lp),       private :: tmp_m2             ! extrapolated pressure, divergence, pressure matvecs
         
         public ::             init_2Dh_axisym
         public ::   is_initialized_2Dh_axisym
         public ::     prepare_step_2Dh_axisym
         public ::    make_coupling_2Dh_axisym
         public ::           makefp_2Dh_axisym
         public ::             ophx_2Dh_axisym
         public ::        ophinv_uz_2Dh_axisym
         public ::     ophinv_urphi_2Dh_axisym
         public ::     incomprp_div_2Dh_axisym
         public ::      incomprp_dw_2Dh_axisym
         public :: incomprp_coupled_2Dh_axisym
         public ::          pres_matvec_axisym

         public :: build_axisym_coeffs
         public :: pressure_matvec_coupled_2Dh_axisym
         public :: solve_pressure_coupled_2Dh_axisym
         public :: enable_axisym_curvature_terms
         
      contains

      !-----------------------------------------------------------------!
      !                                                                 !
      !     PUBLIC STAGE WRAPPERS FOR THE AXISYMMETRIC 2Dh SOLVER       !
      !                                                                 !
      !-----------------------------------------------------------------!

         subroutine init_2Dh_axisym(verbose)
      !! Set up the axisymmetric 2Dh solver. Idempotent. Called by
      !! init_exptA_2Dh_axisym and, as a fallback, by prepare_step_2Dh_axisym;
      !! may also be called by hand from userchk (istep = 0, after usrdat2).
      !! For a torus/helix, t2Dh%init_geom must come first: lambda is read and
      !! cached here (check_torsion_cache stops on the wrong order).
            logical, optional, intent(in) :: verbose
            ! internal
            character(len=*), parameter :: this_procedure = 'init_2Dh_axisym'
            character(len=256) :: msg
            real(dp) :: rmin
            real(dp), external :: glmin

            if (init_2Dh_axisym_done) return

            if (.not. ifaxis) call nek_stop_error('init_2Dh_axisym requires ifaxis = .true.',
     &                                            this_module, this_procedure)
            if (.not. ifaziv) call nek_stop_error('init_2Dh_axisym requires ifaziv = .true.',
     &                                            this_module, this_procedure)
            call check_2Dh_properties(this_procedure)

      ! alpha/R, 1/R^2 and (alpha^2+1)/R^2 are built pointwise and are singular
      ! at R = 0. Meshes touching the axis are not supported yet.
            rmin = glmin(ym1, lx1*ly1*lz1*nelv)
            if (rmin <= 0.0_dp) then
               write(msg,'(A,1pE12.4,A)') 'min(R) = ', rmin,
     &            ': meshes touching the axis are not supported.'
               call nek_stop_error(msg, this_module, this_procedure)
            end if

            call build_wmask(verbose)
            call check_2Dh_masks(this_procedure)

      ! lambda is fixed for the run: build the torsion coefficients once, here.
      ! Without t2Dh%init_geom (pure axisymmetric case) torsion is off.
            call build_torsion_coeffs()
            if (.not. if_torsion_zero) call check_full_section_axisym(this_procedure)
            msg = 'torsion terms '//merge('inactive.', 'active.  ', if_torsion_zero)
            call nek_log_message(msg, this_module, this_procedure)

            init_2Dh_axisym_done = .true.
            call nek_log_message('axisymmetric 2Dh solver initialised.', this_module, this_procedure)
         end subroutine init_2Dh_axisym

         logical function is_initialized_2Dh_axisym()
            is_initialized_2Dh_axisym = init_2Dh_axisym_done
         end function is_initialized_2Dh_axisym

         subroutine prepare_step_2Dh_axisym(alpha)
      !! Coefficients and checks for the current wavenumber. Called first by
      !! nek_advance_2Dh_axisym.
            real(dp), intent(in) :: alpha
            ! internal
            character(len=*), parameter :: this_procedure = 'prepare_step_2Dh_axisym'
            call init_2Dh_axisym() ! for safety; normally done already
            call check_torsion_cache(this_procedure)
            if (.not. axisym_coeffs_defined .or. alpha /= alpha_cached) call build_axisym_coeffs(alpha)
            if (ifadj) then
               if (.not. if_torsion_zero) call nek_stop_error(
     &            'adjoint mode is not implemented with torsion (lambda /= 0).',
     &            this_module, this_procedure)
               call gradm1(dTdx, dTdy, dTdz, t)          ! grad(U_phi)
               call gradm1(gUZ_Z, gUZ_R, tmp_m1c, vx)    ! grad(U_Z)
               call gradm1(gUR_Z, gUR_R, tmp_m1c, vy)    ! grad(U_R)
            end if
            ! Re/Im pairing is only meaningful for alpha /= 0; npert must match
            ! or the caller has set up the wrong number of perturbation slots.
            alpha_eq_zero = (alpha == 0.0_dp)
            if (.not. alpha_eq_zero .and. .not. if_torsion_zero) call nek_stop_error(
     &         'alpha /= 0 with lambda /= 0 is not supported yet (Step 7 of the 2Dh plan).',
     &         this_module, this_procedure)
            if (alpha_eq_zero) then
               if (lpert /= 1) call nek_stop_error('alpha = 0 requires npert = 1.', this_module, this_procedure)
            else
               if (lpert /= 2) call nek_stop_error('alpha /= 0 requires npert = 2.', this_module, this_procedure)
            end if
         end subroutine prepare_step_2Dh_axisym

         subroutine make_coupling_2Dh_axisym(gradw)
      !! Slot-coupled explicit terms of slot jp (read jp, write ipert): the
      !! alpha/lambda pressure gradient of the u_phi equation, the swirl advection
      !! -i alpha (U_phi/R) u, and the torsion alpha*lambda coupling. Runs for all
      !! slots before any slot is advanced.
            real(dp), dimension(lv,lpert), intent(inout) :: gradw
            call extrapprp(tmp_m2)
            call compute_gradp_axisym(gradw, tmp_m2, alphaR_coef)
            call compute_swirl_advection_axisym(advZ, advR, advPhi, jp)
            if (.not. alpha_eq_zero .and. .not. if_torsion_zero) then
               call add_torsion_coupling(advZ, advR, advPhi, alpha_cached, jp)
            end if
         end subroutine make_coupling_2Dh_axisym

         subroutine ophx_2Dh_axisym(resv_jp, h1, h2)
      !! resv -= H u for all slots (cf. ophx in cresvipp). u_Z per slot; u_R of
      !! slot jp couples to u_phi of the other slot at alpha /= 0.
            real(dp), dimension(lv,3,lpert), intent(inout) :: resv_jp
            real(dp), dimension(lv),         intent(in)    :: h1, h2
            ! internal
            real(dp) :: csign_jp
            integer :: ipert, ntot1
            ntot1 = lx1*ly1*lz1*nelv
            do jp = 1, npert
               call helmholtz_matvec_2Dh_axisym(tmp_m1a, vxp(1,jp), h1, h2, h2z_shift, 1)
               call sub2(resv_jp(1,1,jp), tmp_m1a, ntot1)
            end do
            if (alpha_eq_zero) then
               do jp = 1, npert
                  call helmholtz_matvec_2Dh_axisym(tmp_m1a, vyp(1,jp),  h1, h2, diag_shift_r, 1)
                  call sub2(resv_jp(1,2,jp), tmp_m1a, ntot1)
                  call helmholtz_matvec_2Dh_axisym(tmp_m1b, tp(1,1,jp), h1, h2, diag_shift,   1)
                  call sub2(resv_jp(1,3,jp), tmp_m1b, ntot1)
               end do
            else
               do jp = 1, npert
                  ipert = npert + 1 - jp
                  csign_jp = merge(1.0_dp, -1.0_dp, jp == 1)
                  call coupled_helmholtz_matvec_2Dh_axisym(tmp_m1a, tmp_m1b, vyp(1,jp), tp(1,1,ipert),
     &                                                     h1, h2, diag_shift_r, diag_shift, couple_coef, csign_jp)
                  call sub2(resv_jp(1,2,jp),    tmp_m1a, ntot1)
                  call sub2(resv_jp(1,3,ipert), tmp_m1b, ntot1)
               end do
            end if
         end subroutine ophx_2Dh_axisym

         subroutine ophinv_uz_2Dh_axisym(resv_jp, h1, h2)
      !! u_Z increments of all slots: assemble, mask, solve, update vxp.
            real(dp), dimension(lv,3,lpert), intent(inout) :: resv_jp
            real(dp), dimension(lv),         intent(in)    :: h1, h2
            ! internal
            character(len=2) :: info_str
            character(len=7) :: helm_info
            integer :: ntot1
            ntot1 = lx1*ly1*lz1*nelv
            do jp = 1, npert
               info_str = merge('Re', 'Im', jp == 1)
               call dssum(resv_jp(1,1,jp), lx1, ly1, lz1)
               call col2 (resv_jp(1,1,jp), v1mask, ntot1)
               helm_info = info_str//' V_Z  '
               if (istep < 10) call chktcg1(tolhv, resv_jp(1,1,jp), h1, h2, v1mask, vmult, imesh, 1)
               shift_op => h2z_shift
               call solve_helmholtz_2Dh(tmp_m1a, resv_jp(1,1,jp), h1, h2, v1mask, vmult, imesh,
     &                                  tolhv, nmxv, 1, binvm1, helm_info, helm_matvec_axisym)
               call add2(vxp(1,jp), tmp_m1a, ntot1)
            end do
         end subroutine ophinv_uz_2Dh_axisym

         subroutine ophinv_urphi_2Dh_axisym(resv_jp, h1, h2)
      !! u_R/u_phi increments of all slots, decoupled at alpha = 0, coupled across
      !! the slots otherwise; updates vyp and tp.
            real(dp), dimension(lv,3,lpert), intent(inout) :: resv_jp
            real(dp), dimension(lv),         intent(in)    :: h1, h2
            ! internal
            character(len=2) :: info_str
            character(len=7) :: helm_info
            real(dp) :: csign_jp
            integer :: ipert, ntot1
            ntot1 = lx1*ly1*lz1*nelv
            if (alpha_eq_zero) then
               do jp = 1, npert
                  info_str = merge('Re', 'Im', jp == 1)
                  call dssum(resv_jp(1,2,jp), lx1, ly1, lz1)
                  call col2 (resv_jp(1,2,jp), v2mask, ntot1)
                  helm_info = info_str//' V_R  '
                  if (istep < 10) call chktcg1(tolhv, resv_jp(1,2,jp), h1, h2, v2mask, vmult, imesh, 1)
                  shift_op => diag_shift_r
                  call solve_helmholtz_2Dh(tmp_m1a, resv_jp(1,2,jp), h1, h2, v2mask, vmult, imesh,
     &                                     tolhv, nmxv, 1, binvm1, helm_info, helm_matvec_axisym)
                  
                  call dssum(resv_jp(1,3,jp), lx1, ly1, lz1)
                  call col2 (resv_jp(1,3,jp), wmask, ntot1)
                  helm_info = info_str//' V_PHI'
                  if (istep < 10) call chktcg1(tolhv, resv_jp(1,3,jp), h1, h2, wmask, vmult, imesh, 1)
                  shift_op => diag_shift
                  call solve_helmholtz_2Dh(tmp_m1b, resv_jp(1,3,jp), h1, h2, wmask, vmult, imesh,
     &                                     tolhv, nmxv, 1, binvm1, helm_info, helm_matvec_axisym)
                  
                  call add2(vyp(1,jp),  tmp_m1a, ntot1)
                  call add2(tp(1,1,jp), tmp_m1b, ntot1)
               end do
            else
               do jp = 1, npert
                  ipert = npert + 1 - jp
                  csign_jp = merge(1.0_dp, -1.0_dp, jp == 1)
                  helm_info = merge('BLK A  ', 'BLK B  ', jp == 1)
                  call solve_coupled_helmholtz_2Dh_axisym(tmp_m1a, tmp_m1b, resv_jp(1,2,jp), resv_jp(1,3,ipert),
     &                                                    h1, h2, v2mask, wmask, vmult, imesh, tolhv, nmxv, binvm1,
     &                                                    helm_info, diag_shift_r, diag_shift, couple_coef, csign_jp)
                  call add2(vyp(1,jp),     tmp_m1a, ntot1)
                  call add2(tp(1,1,ipert), tmp_m1b, ntot1)
               end do
            end if
         end subroutine ophinv_urphi_2Dh_axisym

         subroutine incomprp_div_2Dh_axisym(dpr_jp)
      !! Right-hand side of the pressure step of slot jp: -(D u + D_alpha u_phi).
      !! The alpha part reads the other slot's u_phi.
            real(dp), dimension(lp), intent(out) :: dpr_jp
            ! internal
            integer :: ntot2
            ntot2 = lx2*ly2*lz2*nelv
            call opdiv (dpr_jp, vxp(1,jp), vyp(1,jp), vzp(1,jp))
            call compute_frc_div_axisym(tmp_m2, tp, alphaR_coef)
            call add2  (dpr_jp, tmp_m2, ntot2)
            call chsign(dpr_jp, ntot2)
         end subroutine incomprp_div_2Dh_axisym

         subroutine incomprp_dw_2Dh_axisym(dw, dpr, h2inv)
      !! u_phi correction of slot jp from the pressure corrections.
            real(dp), dimension(lv), intent(out) :: dw
            real(dp), dimension(lp,lpert), intent(in) :: dpr
            real(dp), dimension(lv), intent(in) :: h2inv
            call compute_dw_axisym(dw, dpr, h2inv, alphaR_coef)
         end subroutine incomprp_dw_2Dh_axisym

         subroutine incomprp_coupled_2Dh_axisym(dpr, h1, h2, h2inv, alpha, ebar)
      !! alpha /= 0 and lambda /= 0: both slots' pressure corrections in one
      !! stacked solve. No projection. Deleted in Step 7.
            real(dp), dimension(lp,lpert), intent(inout) :: dpr
            real(dp), dimension(lv), intent(in) :: h1, h2, h2inv
            real(dp), intent(in) :: alpha, ebar
            ! internal
            integer :: iter
            do jp = 1, npert
               call incomprp_div_2Dh_axisym(dpr(1,jp))
            end do
            call solve_pressure_coupled_2Dh_axisym(dpr(1,1), dpr(1,2), h1, h2, h2inv,
     &                                             alphaR_coef, alpha, ebar, 1, iter)
         end subroutine incomprp_coupled_2Dh_axisym

      !-----------------------------------!
      !                                   !
      !     PRIVATE HELPER ROUTINES       !
      !                                   !
      !-----------------------------------!
         
         subroutine build_axisym_coeffs(alpha)
            implicit none
            real(dp), intent(in) :: alpha
            integer :: ix, iy, iz, ie
            real(dp) :: r
            
            do ie = 1, nelv
               do iz = 1, lz1
               do iy = 1, ly1
               do ix = 1, lx1
                  r = ym1(ix,iy,iz,ie)
                  alphaR_coef (ix,iy,iz,ie) = alpha / r
                  h2z_shift   (ix,iy,iz,ie) = alpha**2 / r**2
                  diag_shift  (ix,iy,iz,ie) = h2z_shift(ix,iy,iz,ie) + 1.0_dp / r**2
                  diag_shift_r(ix,iy,iz,ie) = diag_shift(ix,iy,iz,ie)
                  couple_coef (ix,iy,iz,ie) = 2.0_dp * alpha / r**2
               end do
               end do
               end do
            end do
            
            alpha_cached          = alpha
            axisym_coeffs_defined  = .true.
         end subroutine build_axisym_coeffs
         
         subroutine helmholtz_matvec_2Dh_axisym(Au, u, h1, h2, shift, isd)
            implicit none
            ! Au = (h1*A + h2*B) u + vdiff*shift*bm1*u: Nek's Helmholtz operator
            ! plus the diagonal wavenumber/curvature shift (field).
            real, dimension(lx1,ly1,lz1,1), intent(out) :: Au
            real, dimension(lx1,ly1,lz1,1), intent(in)  :: u, h1, h2, shift
            integer, intent(in) :: isd
            integer :: imesh, ntot1
            call neklab_timer_start(t_axisym_helm_matvec)
            imesh = 1
            ntot1 = lx1*ly1*lz1*nelv
            call axhelm(Au, u, h1, h2, imesh, isd)
            call col3 (tmp_m1c, u, vdiff(1,1,1,1,1), ntot1)
            call col2 (tmp_m1c, shift, ntot1)
            call col2 (tmp_m1c, bm1, ntot1)
            call add2 (Au, tmp_m1c, ntot1)
            call neklab_timer_stop(t_axisym_helm_matvec)
         end subroutine helmholtz_matvec_2Dh_axisym
         
         subroutine compute_gradp_axisym(gradp, prextr, alphaR)
      !! Third-component pressure gradient, accumulated over the make_coupling
      !! loop jp = 1, npert (zeroed at jp = 1): the alpha part, from slot jp
      !! into slot ipert (skew, slot swap), exact transpose of the alpha block
      !! of compute_frc_div_axisym. The lambda part is explicit, in
      !! add_torsion_forcing_phi (route A).
            implicit none
            include 'SIZE'
            include 'SOLN'
            include 'MASS'
            real(dp), dimension(lv,lpert), intent(inout) :: gradp
            real(dp), dimension(lp), intent(in) :: prextr
            real(dp), dimension(lx1,ly1,lz1,lelv), intent(in) :: alphaR
            ! internal
            integer :: ipert, spert, kp, ntot1

            ntot1 = lx1*ly1*lz1*nelv
            ipert = npert + 1 - jp

            if (jp == 1) then
               do kp = 1, npert
                  call rzero(gradp(1,kp), ntot1)
               end do
            end if

            if (alpha_eq_zero) return

      ! --- alpha part: skew in slot space (swaps jp <-> ipert). Absent at
      !     alpha = 0.
            if (.not. alpha_eq_zero) then
               spert = merge(1, -1, jp == 1)
               call bm2_map12T_2Dh(tmp_m1a, prextr)          ! I12^T (bm2 p)
               call col2          (tmp_m1a, alphaR, ntot1)
               if (spert == -1) call chsign(tmp_m1a, ntot1)
               call add2          (gradp(1,ipert), tmp_m1a, ntot1)
            end if

      ! dssum is applied to the velocity residual later
         end subroutine compute_gradp_axisym
         
         subroutine compute_dw_axisym(dw, dpr, h2inv, alphaR)
            implicit none
            include 'SIZE'
            include 'SOLN'
            include 'MASS'
            real(dp), dimension(lv), intent(out) :: dw
            real(dp), dimension(lp,lpert), intent(in) :: dpr
            real(dp), dimension(lx1,ly1,lz1,lelv), intent(in) :: h2inv, alphaR
            ! internal
            integer :: ipert, spert, ntot1
            
            ntot1 = lx1*ly1*lz1*nelv
            call rzero(dw, ntot1)
            
            if (alpha_eq_zero) return
            
            if (.not. alpha_eq_zero) then
               ipert = npert + 1 - jp
               spert = merge(1, -1, jp == 1)
               call bm2_map12T_2Dh(dw, dpr(1,ipert))      ! I12^T (bm2 dp)
               call col2          (dw, alphaR, ntot1)
               if (spert == 1) call chsign(dw, ntot1)     ! matches -spert*beta_z sign convention
            end if
            
            call col2 (dw, wmask, ntot1)
            call dssum(dw, lx1, ly1, lz1)
            call col2 (dw, binvm1, ntot1)
            call col2 (dw, h2inv, ntot1)
         end subroutine compute_dw_axisym
         
         subroutine compute_frc_div_axisym(frc_div, w, alphaR)
            implicit none
            include 'SIZE'
            include 'SOLN'
            include 'MASS'
            real(dp), dimension(lx2,ly2,lz2,lelv), intent(out) :: frc_div
            real(dp), dimension(lv,ldimt,lpert), intent(in) :: w
            real(dp), dimension(lx1,ly1,lz1,lelv), intent(in) :: alphaR
            integer :: ipert, spert, ntot1, ntot2
            
            ntot1 = lx1*ly1*lz1*nelv
            ntot2 = lx2*ly2*lz2*nelv

            call rzero(frc_div, ntot2)
            
            if (alpha_eq_zero .and. if_torsion_zero) return

            call rzero(tmp_m1a, ntot1)
            
            ! --- alpha part, from the opposite slot
            if (.not. alpha_eq_zero) then
               ipert = npert + 1 - jp
               spert = merge(1, -1, jp == 1)
               call copy(tmp_m1a, w(1,1,ipert), ntot1)
               call col2(tmp_m1a, alphaR, ntot1)       ! scale on M1 BEFORE mapping to M2
               if (spert == -1) call chsign(tmp_m1a, ntot1)
            end if
            
            ! --- lambda part, from THIS slot: -(lambda/R)*dtheta(u_phi). Accumulated
            !     on M1 alongside the alpha part so a single map12 below carries both.
            if (.not. if_torsion_zero) then
               call dphi_apply(tmp_m1b, w(1,1,jp))
               call add2      (tmp_m1a, tmp_m1b, ntot1)
            end if
            
            call bm2_map12_2Dh(frc_div, tmp_m1a)           
         end subroutine compute_frc_div_axisym
         
         subroutine makefp_2Dh_axisym()
            ! Perturbation RHS assembly of the axisymmetric 2Dh solver.
            !
            ! Replaces the core `makefp` + `makeqp` pair. Structure is identical
            ! to those routines (perturb.f), with the swirl, curvature and torsion
            ! explicit terms
            ! injected AFTER the convective step and BEFORE the extrapolation
            ! step, so that makextp/makeabqp apply EXT3 to them and roll them
            ! into the lag arrays (exx1p/exx2p, vgradt1p/vgradt2p) automatically.
            !
            ! Adjoint mode (ifadj): advabp_adjoint/convabp give Nek's adjoint
            ! transport; the swirl advection flips sign (compute_swirl_advection_axisym)
            ! and the curvature block is transposed (add_axisym_curv_R/_phi).
            !
            ! Operates on the current global `jp` (the Nek routines called here
            ! all read it), and restores `ifield` on exit.
            implicit none
            include 'SIZE'
            include 'SOLN'
            include 'INPUT'
            include 'TSTEP'
            include 'ADJOINT'
            ! internal
            character(len=*), parameter :: this_procedure = 'makefp_2Dh_axisym'
            integer :: ifield_bak, ntot1
            
            call neklab_timer_start(t_axisym_makefp)
            ntot1 = lx1*ly1*lz1*nelv
            ifield_bak = ifield
            
            ! ---------------- velocity: u_Z, u_R ----------------
            ifield = 1
            call makeufp
            if (ifnav .and. .not. ifchar) then
               if (ifadj) then
                  call advabp_adjoint_axisym           ! +(U.grad)u - (grad U)^T u - u_phi grad(U_phi)
                  ! we don't use the nek core advabp_adjoint because it ignores the R weight of bm1 under ifaxis
               else
                  call advabp
               end if
            end if
            call add2(bfxp(1,jp), advZ(1,jp), ntot1)   ! -+i*alpha*(U_phi/R)*u_Z  (direct/adjoint)
            call add2(bfyp(1,jp), advR(1,jp), ntot1)   ! -+i*alpha*(U_phi/R)*u_R
            call add_axisym_curv_R(bfyp(1,jp), jp)     ! direct: +2 U_phi u_phi/R, adjoint: -U_phi u_phi/R
            if (.not. if_torsion_zero) call add_torsion_forcing_RZ(bfxp(1,jp), bfyp(1,jp), jp) ! viscous advection
            if (iftran) call makextp                   ! <-- EXT3
            call makebdfp
            call lagfieldp
            
            ! ---------------- scalar: u_phi ----------------
            ifield = 2
            if (ifadvc(ifield) .and. ifchar) call nek_stop_error(
     &         'characteristics (ifchar) not supported for the swirl field',
     &         this_module, this_procedure)
            
            call makeuqp
            if (ifadvc(ifield) .and. .not. ifchar) call convabp
            call add2(bqp(1,ifield-1,jp), advPhi(1,jp), ntot1)   ! -+i*alpha*(U_phi/R)*u_phi
            call add_axisym_curv_phi(bqp(1,ifield-1,jp), jp)     ! direct: -(U_phi u_R + U_R u_phi)/R
                                                                 ! adjoint: (2 U_phi u_R - U_R u_phi)/R
            if (.not. if_torsion_zero) call add_torsion_forcing_phi(bqp(1,ifield-1,jp), jp) ! viscous advection D_lam^T p
            if (iftran) call makeabqp                            ! <-- EXT3
            if ((iftran .and. .not. ifchar) .or.
     &          (iftran .and. .not. ifadvc(ifield) .and. ifchar)) call makebdqp
            call lagscalp
            
            ifield = ifield_bak
            call neklab_timer_stop(t_axisym_makefp)
         end subroutine makefp_2Dh_axisym
         
         subroutine add_axisym_curv_R(rhs, jp_)
            implicit none
            include 'SIZE'
            include 'SOLN'   ! vx,vy,t (baseflow); vxp,vyp,tp (perturbation); vtrans
            include 'MASS'   ! bm1
            real(dp), dimension(lv), intent(inout) :: rhs
            integer, intent(in) :: jp_
            ! internal
            integer :: ntot1
            
            ntot1 = lx1*ly1*lz1*nelv
            
            ! --- R-momentum: + 2*U_phi*u'_phi/R * rho * bm1
            !     linearization of the missing -U_phi^2/R term in (U.grad U)_R
            call col3   (tmp_m1a, t(1,1,1,1,1), tp(1,1,jp_), ntot1)     ! U_phi * u'_phi
            call invcol2(tmp_m1a, ym1, ntot1)                           ! / R
            call cmult  (tmp_m1a, merge(-1.0_dp, 2.0_dp, ifadj), ntot1) ! adjoint: transposed block
            call col2   (tmp_m1a, vtrans(1,1,1,1,1), ntot1)             ! rho
            call col2   (tmp_m1a, bm1, ntot1)                           ! bm1
            call add2   (rhs, tmp_m1a, ntot1)
            
         end subroutine add_axisym_curv_R

         subroutine advabp_adjoint_axisym()
            !! Adjoint transport of the u_Z, u_R equations, slot jp:
            !!   + (U.grad) v  -  (grad U)^T v  -  v_phi grad(U_phi)
            !! Pointwise replacement for Nek's advabp_adjoint, whose convop_adj ignores
            !! the R weight of bm1 under ifaxis and whose grad(T) term reads slot 1 only.
            !! Curvature parts are handled separately (add_axisym_curv_R/_phi).
            integer :: ntot1
            ntot1 = lx1*ly1*lz1*nelv

      ! --- Z equation
            call convop  (tmp_m1a, vxp(1,jp))                  ! U.grad v_Z
            call subcol3 (tmp_m1a, vxp(1,jp),  gUZ_Z, ntot1)   ! - v_Z dU_Z/dZ
            call subcol3 (tmp_m1a, vyp(1,jp),  gUR_Z, ntot1)   ! - v_R dU_R/dZ
            call subcol3 (tmp_m1a, tp(1,1,jp), dTdx,  ntot1)   ! - v_phi dU_phi/dZ
            call col2    (tmp_m1a, vtrans(1,1,1,1,1), ntot1)   ! rho
            call col2    (tmp_m1a, bm1, ntot1)                 ! bm1
            call add2    (bfxp(1,jp), tmp_m1a, ntot1)

      ! --- R equation
            call convop  (tmp_m1b, vyp(1,jp))                  ! U.grad v_R
            call subcol3 (tmp_m1b, vxp(1,jp),  gUZ_R, ntot1)   ! - v_Z dU_Z/dR
            call subcol3 (tmp_m1b, vyp(1,jp),  gUR_R, ntot1)   ! - v_R dU_R/dR
            call subcol3 (tmp_m1b, tp(1,1,jp), dTdy,  ntot1)   ! - v_phi dU_phi/dR
            call col2    (tmp_m1b, vtrans(1,1,1,1,1), ntot1)   ! rho
            call col2    (tmp_m1b, bm1, ntot1)                 ! bm1
            call add2    (bfyp(1,jp), tmp_m1b, ntot1)
         end subroutine advabp_adjoint_axisym
         
         subroutine add_axisym_curv_phi(rhs, jp_)
            implicit none
            include 'SIZE'
            include 'SOLN'   ! vx,vy,t (baseflow); vxp,vyp,tp (perturbation); vtrans
            include 'MASS'   ! bm1
            real(dp), dimension(lv), intent(inout) :: rhs
            integer, intent(in) :: jp_
            ! internal
            integer :: ntot1
            
            ntot1 = lx1*ly1*lz1*nelv
            
            ! --- swirl ("phi") equation: - (U_phi*u'_R + U_R*u'_phi)/R * rhocp * bm1
            !     linearization of the missing +U_R*U_phi/R term in (U.grad U)_phi
            call col3   (tmp_m1a, t(1,1,1,1,1), vyp(1,jp_),  ntot1)    ! U_phi * u'_R
            if (ifadj) call cmult(tmp_m1a, -2.0_dp, ntot1)             ! adjoint: -2 U_phi u_R
            call col3   (tmp_m1b, vy,           tp(1,1,jp_), ntot1)    ! U_R   * u'_phi  (drop if U_R deemed negligible)
            call add2   (tmp_m1a, tmp_m1b, ntot1)
            call invcol2(tmp_m1a, ym1, ntot1)                          ! / R
            call chsign (tmp_m1a, ntot1)
            call col2   (tmp_m1a, vtrans(1,1,1,1,2), ntot1)            ! rhocp
            call col2   (tmp_m1a, bm1, ntot1)                          ! bm1
            call add2   (rhs, tmp_m1a, ntot1)
            
         end subroutine add_axisym_curv_phi
         
         subroutine compute_swirl_advection_axisym(advZ, advR, advPhi, jp_)
            ! Swirl advection -i*alpha*(U_phi/R)*u, present in ALL THREE momentum
            ! equations because the base flow has U_phi != 0 in the ignorable
            ! direction. Planar analogue: compute_W_advection_planar.
            implicit none
            include 'SIZE'
            include 'SOLN'   ! vxp,vyp,tp; t (baseflow); vtrans
            include 'MASS'   ! bm1
            real(dp), dimension(lv,lpert), intent(inout) :: advZ, advR, advPhi
            integer, intent(in) :: jp_
            ! internal
            integer :: ipert, ntot1
            
            ntot1 = lx1*ly1*lz1*nelv
            ipert = npert + 1 - jp_
            call rzero(advZ(1,ipert),   ntot1)
            call rzero(advR(1,ipert),   ntot1)
            call rzero(advPhi(1,ipert), ntot1)
            
            if (alpha_eq_zero) return
            
            call col3(tmp_m1a, alphaR_coef, t(1,1,1,1,1), ntot1)   ! (alpha/R) * U_phi
            if (ifadj) call chsign(tmp_m1a, ntot1)                 ! antisymmetric in slot space: adjoint = -direct
            if (ipert == 1) call chsign(tmp_m1a, ntot1)            ! sign keyed to the READING jp, not the caller
            
            call col3 (advZ(1,ipert), tmp_m1a, vxp(1,jp_), ntot1)
            call col2 (advZ(1,ipert), vtrans(1,1,1,1,1), ntot1)      ! rho
            call col2 (advZ(1,ipert), bm1, ntot1)
            
            call col3 (advR(1,ipert), tmp_m1a, vyp(1,jp_), ntot1)
            call col2 (advR(1,ipert), vtrans(1,1,1,1,1), ntot1)      ! rho
            call col2 (advR(1,ipert), bm1, ntot1)
            
            call col3 (advPhi(1,ipert), tmp_m1a, tp(1,1,jp_), ntot1)
            call col2 (advPhi(1,ipert), vtrans(1,1,1,1,2), ntot1)    ! rho*cp
            call col2 (advPhi(1,ipert), bm1, ntot1)
         end subroutine compute_swirl_advection_axisym

         subroutine enable_axisym_curvature_terms()
      !! Adds the curvature terms of (U.grad)U that Nek's axisymmetric advection
      !! leaves out, to the nonlinear solver (ipert = 0):
      !!   R  : + U_phi^2 / R
      !!   phi: - U_R U_phi / R
      !! Call once from userchk (istep = 0). Idempotent.
            call add_neklab_forcing_term(axisym_curv_forcing, 'axisym_curvature')
            call add_neklab_source_term (axisym_curv_source,  'axisym_curvature')
         end subroutine enable_axisym_curvature_terms

         subroutine axisym_curv_forcing(fx, fy, fz, ix, iy, iz, e)
            real(dp), intent(inout) :: fx, fy, fz
            integer, intent(in) :: ix, iy, iz, e
            real(dp) :: r
            r = ym1(ix, iy, iz, e)
            if (r > 0.0_dp) fy = fy + t(ix, iy, iz, e, 1)**2 / r
         end subroutine axisym_curv_forcing

         subroutine axisym_curv_source(q, ix, iy, iz, e)
            real(dp), intent(inout) :: q
            integer, intent(in) :: ix, iy, iz, e
            real(dp) :: r
            r = ym1(ix, iy, iz, e)
            if (r > 0.0_dp) q = q - vy(ix, iy, iz, e)*t(ix, iy, iz, e, 1) / r
         end subroutine axisym_curv_source
         
         subroutine pressure_matvec_2Dh_axisym(ap, wp, h1, h2, h2inv, alphaR, intype)
            implicit none
            real(dp), dimension(lx2,ly2,lz2,lelv), intent(out) :: Ap
            real(dp), dimension(lx2,ly2,lz2,lelv), intent(in)  :: wp
            real(dp), dimension(lx1,ly1,lz1,lelv), intent(in)  :: h1, h2, h2inv, alphaR
            integer, intent(in) :: intype
            ! internal
            integer :: ntot1, ntot2
            
            call neklab_timer_start(t_axisym_pres_matvec)
            call cdabdtp(Ap, wp, h1, h2, h2inv, intype)   ! core Nek, already ifaxis-correct globally
            
            if (alpha_eq_zero) then
               call neklab_timer_stop(t_axisym_pres_matvec)
               return
            end if
            
            ntot1 = lx1*ly1*lz1*nelv
            ntot2 = lx2*ly2*lz2*nelv
            
            ! --- alpha block: (alpha/R) H^-1 (alpha/R). The two slot swaps of
            !     D_phi and D_phi^T cancel, so this block is diagonal in slot space
            !     and each perturbation's pressure can be solved on its own.
            if (.not. alpha_eq_zero) then
               call bm2_map12T_2Dh(tmp_m1a, wp)                    ! I12^T (bm2 p)
               call col2          (tmp_m1a, alphaR, ntot1)         ! alpha/R
               call col2          (tmp_m1a, wmask,  ntot1)
               call dssum         (tmp_m1a, lx1, ly1, lz1)
               call col2          (tmp_m1a, binvm1, ntot1)
               call col2          (tmp_m1a, h2inv,  ntot1)
               call col2          (tmp_m1a, alphaR, ntot1)         ! alpha/R
               call chsign        (tmp_m1a, ntot1)
               call bm2_map12_2Dh (tmp_m2, tmp_m1a)                ! bm2 * I12 (.)
               call sub2(ap, tmp_m2, ntot2)
            end if
            
            call neklab_timer_stop(t_axisym_pres_matvec)
         end subroutine pressure_matvec_2Dh_axisym
         
         subroutine pmatvec_cross_block(out2, in2, h2inv, alphaR, left_is_alpha)
            !! Applies ONE cross block of the coupled pressure Schur complement to a
            !! single M2 pressure field. The Schur complement is
            !!   E = D_phi B^-1 D_phi^T,   D_phi = A S + L,  A = alpha/R, L = -(lambda/R) dtheta,
            !! whose slot-off-diagonal ("cross") part is
            !!   term2 = A S B^-1 L^T      (left op = A, right op = L^T)
            !!   term3 = -L B^-1 A S = term2^T   (left op = L, right op = A)
            !! left_is_alpha selects which: .true. builds A(.)B^-1 L^T(.), .false.
            !! builds L(.)B^-1 A(.). The slot swap S is NOT applied here; the caller
            !! routes inputs/outputs between slots and supplies the swap sign.
            !!
            !! The B^-1 sandwich (bm2_map12T, wmask, dssum, binvm1, h2inv, bm2_map12)
            !! is identical to the diagonal blocks in pressure_matvec_2Dh_axisym;
            !! only the left/right operator halves differ, which keeps E symmetric
            !! (term3 = term2^T).
            implicit none
            real(dp), dimension(lx2,ly2,lz2,lelv), intent(out) :: out2
            real(dp), dimension(lx2,ly2,lz2,lelv), intent(in)  :: in2
            real(dp), dimension(lx1,ly1,lz1,lelv), intent(in)  :: h2inv, alphaR
            logical, intent(in) :: left_is_alpha
            integer :: ntot1

            ntot1 = lx1*ly1*lz1*nelv

            call bm2_map12T_2Dh       (tmp_m1a, in2)              ! I12^T (bm2 p)

            ! right operator half (applied first): L^T if left is alpha, else A
            if (left_is_alpha) then
               call dphi_apply_t      (tmp_m1b, tmp_m1a)          ! L^T
            else
               call dphi_alpha_apply_t(tmp_m1b, tmp_m1a, alphaR)  ! A^T = A
            end if

            call col2 (tmp_m1b, wmask, ntot1)
            call dssum(tmp_m1b, lx1, ly1, lz1)
            call col2 (tmp_m1b, binvm1, ntot1)
            call col2 (tmp_m1b, h2inv, ntot1)

            ! left operator half (applied last): A if left is alpha, else L
            if (left_is_alpha) then
               call dphi_alpha_apply  (tmp_m1a, tmp_m1b, alphaR)    ! A
            else
               call dphi_apply        (tmp_m1a, tmp_m1b)            ! L
            end if
            call chsign               (tmp_m1a, ntot1)
            call bm2_map12_2Dh        (out2, tmp_m1a)               ! bm2 * I12 (.)
         end subroutine pmatvec_cross_block

         subroutine pressure_matvec_coupled_2Dh_axisym(apR, apI, wpR, wpI, h1, h2, h2inv, alphaR, intype)
      !! Coupled pressure Schur-complement matvec for alpha /= 0 AND lambda /= 0.
      !! Produces both slots at once:
      !!   [apR]   [ E_diag   E_cross ] [wpR]
      !!   [apI] = [ E_cross^T E_diag ] [wpI]
      !! The diagonal blocks are exactly pressure_matvec_2Dh_axisym (per slot).
      !! The cross blocks come from D_phi's alpha part carrying the slot swap S
      !! while its lambda part does not:
      !!   E_cross acting between slots = A S B^-1 L^T  +  (-L B^-1 A S)
      !! With the swap S sending slot R<->I and contributing the sign pattern of
      !! multiplication by i (Re[i w]=-w^I, Im[i w]=+w^R):
      !!   apR += (term2 reads wpI) - (term3 reads wpI)   [both from the I slot]
      !!   apI += (term2 reads wpR) - (term3 reads wpR)
      !! with the appropriate swap sign. The two cross contributions are exact
      !! transposes, so the full 2x2 operator is symmetric.
            implicit none
            real(dp), dimension(lx2,ly2,lz2,lelv), intent(out) :: apR, apI
            real(dp), dimension(lx2,ly2,lz2,lelv), intent(in)  :: wpR, wpI
            real(dp), dimension(lx1,ly1,lz1,lelv), intent(in)  :: h1, h2, h2inv, alphaR
            integer, intent(in) :: intype
            real(dp), dimension(lx2,ly2,lz2,lelv) :: c2a, c2b
            integer :: ntot2
            
            ntot2 = lx2*ly2*lz2*nelv
            
            ! diagonal blocks (alpha^2 + lambda^2 parts), reuse existing routine
            call pressure_matvec_2Dh_axisym(apR, wpR, h1, h2, h2inv, alphaR, intype)
            call pressure_matvec_2Dh_axisym(apI, wpI, h1, h2, h2inv, alphaR, intype)
            
      ! cross blocks. S swaps slots and is ANTISYMMETRIC (S = [[0,1],[-1,0]]),
      ! so the two off-diagonal blocks carry OPPOSITE signs. Writing
      ! K = term2 - term3 = A B^-1 L^T - L B^-1 A (with term3 = term2^T):
      !   apR += K(pI)        (block R<-I, S-sign +1)
      !   apI -= K(pR)        (block I<-R, S-sign -1)
      ! Then block(R<-I)^T = K^T = -K = block(I<-R), so the full 2x2 E is
      ! SYMMETRIC, as a Schur complement must be. (The equal-sign version is
      ! antisymmetric and wrong; caught by the symmetry check.)
            call pmatvec_cross_block(c2a, wpI, h2inv, alphaR, .true.)
            call pmatvec_cross_block(c2b, wpI, h2inv, alphaR, .false.)
            call add2   (apR, c2a, ntot2)
            call sub2   (apR, c2b, ntot2)
            
            call pmatvec_cross_block(c2a, wpR, h2inv, alphaR, .true.)
            call pmatvec_cross_block(c2b, wpR, h2inv, alphaR, .false.)
            call sub2   (apI, c2a, ntot2)
            call add2   (apI, c2b, ntot2)
         end subroutine pressure_matvec_coupled_2Dh_axisym
         
         subroutine coupled_helmholtz_matvec_2Dh_axisym(y1, y2, x1, x2, h1, h2, shift1, shift2, couple_coef, csign)
            implicit none
            ! Matvec for the coupled 2x2 block (x1 = u_R, x2 = u_phi):
            !   [ D1  csign*c ] [x1]   [y1]
            !   [csign*c   D2 ] [x2] = [y2]
            ! D1 = axhelm + shift1 (u_R  shift, = diag_shift_r: gains lambda^2/R^2)
            ! D2 = axhelm + shift2 (u_phi shift, = diag_shift: no lambda term)
            ! c  = couple_coef (== 2*alpha/R^2, field-weighted)
            ! shift1 == shift2 at lambda = 0 (no torsion).
            real(dp), dimension(lx1,ly1,lz1,lelv), intent(out) :: y1, y2
            real(dp), dimension(lx1,ly1,lz1,lelv), intent(in)  :: x1, x2, h1, h2, shift1, shift2, couple_coef
            real(dp), intent(in) :: csign
            integer :: ntot1
            ntot1 = lx1*ly1*lz1*nelv
            
            call neklab_timer_start(t_axisym_cpl_matvec)
            call helmholtz_matvec_2Dh_axisym(y1, x1, h1, h2, shift1, 1)
            call helmholtz_matvec_2Dh_axisym(y2, x2, h1, h2, shift2, 1)
            
            call col3   (tmp_m1c, x2, couple_coef, ntot1)
            call col2   (tmp_m1c, vdiff(1,1,1,1,1), ntot1)
            call col2   (tmp_m1c, bm1, ntot1)
            call add2s2 (y1, tmp_m1c, csign, ntot1)
            
            call col3   (tmp_m1c, x1, couple_coef, ntot1)
            call col2   (tmp_m1c, vdiff(1,1,1,1,1), ntot1)
            call col2   (tmp_m1c, bm1, ntot1)
            call add2s2 (y2, tmp_m1c, csign, ntot1)
            call neklab_timer_stop(t_axisym_cpl_matvec)
         end subroutine coupled_helmholtz_matvec_2Dh_axisym

         subroutine helm_matvec_axisym(Au, u, h1, h2)
            ! Wrapper around helmholtz_matvec_2Dh_axisym matching the abstract interface.
            ! The shift field is selected through shift_op.
            real(dp), dimension(lv), intent(out) :: Au
            real(dp), dimension(lv), intent(in)  :: u, h1, h2
            call helmholtz_matvec_2Dh_axisym(Au, u, h1, h2, shift_op, 1)
         end subroutine helm_matvec_axisym

         subroutine pres_matvec_axisym(Ap, p, h1, h2, h2inv, intype)
            ! Wrapper around pressure_matvec_2Dh_axisym matching the abstract interface.
            real(dp), dimension(lp), intent(out) :: Ap
            real(dp), dimension(lp), intent(in)  :: p
            real(dp), dimension(lv), intent(in)  :: h1, h2, h2inv
            integer, intent(in) :: intype
            call pressure_matvec_2Dh_axisym(Ap, p, h1, h2, h2inv, alphaR_coef, intype)
         end subroutine pres_matvec_axisym
         
         subroutine solve_coupled_helmholtz_2Dh_axisym(x1,x2,f1,f2,h1,h2,mask1,mask2,mult,imsh,tin,maxit,
     &                                                 binv,name,shift1,shift2,couple_coef,csign)
            implicit none
             ! Same CG loop as solve_helmholtz_2Dh (common), but operating on the
            ! stacked 2-field vector (x1;x2). dssum/mask applied to x1 and x2
            ! SEPARATELY (mirrors the planar-case pattern), never combined
            ! before dssum -- this is the ordering question from before, tested
            ! directly by construction.
            ! There are two separate BC masks:
            ! mask1: field 1 (u_R,  v2mask)
            ! mask2: field 2 (u_phi, wmask)
            include 'FDMH1'
            real, intent(out) :: x1(1), x2(1)
            real, intent(inout)  :: f1(1), f2(1)
            real, intent(in)  :: h1(1), h2(1), mask1(1), mask2(1), mult(1), binv(1)
            integer, intent(in) :: imsh, maxit
            real, intent(in) :: tin
            character(len=7), intent(in) :: name
            real, intent(in) :: shift1(1), shift2(1), couple_coef(1)
            real(dp), intent(in) :: csign
            
            integer :: iter, krylov, n, nel, niter, nxyz
            real :: alpha, beta, alphm, fmax
            real :: rbn0, rbn2, rho, rho0, rtz1, rtz2
            real :: tol, vol
            LOGICAL IFPRINT, IFHZPC
            COMMON /CPRINT/ IFPRINT, IFHZPC
            logical ifdfrm, iffast, ifh2, ifsolv
            common /fastmd/ ifdfrm(lelt), iffast(lelt), ifh2, ifsolv
            logical ifprint_helm
            integer, parameter :: lg = lx1*ly1*lz1*lelt
            ! stacked scratch: [1:n] = field 1, [n+1:2n] = field 2
            real, dimension(2*lg) :: d, r, w, p, z
            real :: scalar(2)
            integer :: niterhm
            common /iterhm/ niterhm
            real, external :: glamax, glsum, glsc2, glsc3, vlsc3, vlsc32
            
            call neklab_timer_start(t_axisym_solve_cpl)
            rho = 0.0
            NXYZ = lx1*ly1*lz1
            NEL  = NELV
            VOL  = VOLVM1
            n = NEL*NXYZ            ! per-field length
            tol = abs(tin)
            if (restol(ifield).ne.0) tol = abs(restol(ifield))
            if (tin.lt.0) tol = abs(tin)
            niter = min(maxit,900)
            
            if (.not.ifsolv) then
            call setfast(h1,h2,imsh)
            ifsolv = .true.
            endif
            ! diagonal preconditioner: reuse setprec on the diag_shift-weighted
            ! operator (ignores the off-diagonal coupling, same spirit as the
            ! existing solver ignoring `shift` in its preconditioner)
            call setprec(d(1), h1,h2,imsh,1)
            call copy   (d(n+1), d(1), n)
            
            call dssum(f1, lx1, ly1, lz1)
            call dssum(f2, lx1, ly1, lz1)
            call copy(r(1),   f1, n);  call col2(r(1),   mask1, n)
            call copy(r(n+1), f2, n);  call col2(r(n+1), mask2, n)
            
            call rzero(x1,n); call rzero(x2,n)
            call rzero(p,2*n)
            
            fmax = max(glamax(r(1),n), glamax(r(n+1),n))
            if (fmax == 0.0) then
      ! Zero RHS: nothing to solve. The stop MUST happen here -- a timer left
      ! running is not restarted by the next `start` (Timer_Utils.f90:185), so
      ! a missed stop silently loses every subsequent call as well.
               call neklab_timer_stop(t_axisym_solve_cpl)
               return
            end if
            
            krylov = 0
            rtz1 = 1.0
            niterhm = 0
            
            do iter=1,niter
               call col3(z(1),   r(1),   d(1),   n)
               call col3(z(n+1), r(n+1), d(n+1), n)
               rtz2 = rtz1
               scalar(1) = vlsc3(z(1),r(1),mult,n) + vlsc3(z(n+1),r(n+1),mult,n)
               if (param(18).eq.1) then
                  scalar(2) = vlsc3(r(1),r(1),mult,n) + vlsc3(r(n+1),r(n+1),mult,n)
               else
                  scalar(2) = vlsc32(r(1),mult,binv,n) + vlsc32(r(n+1),mult,binv,n)
               endif
               call gop(scalar,w,'+  ',2)
               rtz1 = scalar(1)
               rbn2 = sqrt(scalar(2)/(2.0*vol))
               if (iter.eq.1) rbn0 = rbn2
               if (param(22).lt.0) tol = abs(param(22))*rbn0
               if (tin.lt.0)       tol = abs(tin)*rbn0
               
               ifprint_helm = .false.
               if (nio.eq.0.and.ifprint.and.param(74).ne.0) ifprint_helm=.true.
               if (nio.eq.0.and.istep.eq.1)                 ifprint_helm=.true.
               if (ifprint_helm) write(6,3002) istep,'  CplHmh ' // name, iter, rbn2, h1(1), tol, h2(1)
               
               IF (rbn2.LE.TOL.and.(iter.gt.1 .or. istep.le.5)) THEN
                  NITER = ITER-1
                  if (nio.eq.0) write(6,3000) istep,'  CplHmh ' // name, niter, rbn2, rbn0, tol
                  goto 9999
               endif
               
               beta = rtz1/rtz2
               if (iter.eq.1) beta = 0.0
               call add2s1(p, z, beta, 2*n)
               
               call coupled_helmholtz_matvec_2Dh_axisym(w(1), w(n+1), p(1), p(n+1),
     &                                                  h1, h2, shift1, shift2, couple_coef, csign)
               call dssum(w(1),   lx1, ly1, lz1); call col2(w(1),   mask1, n)
               call dssum(w(n+1), lx1, ly1, lz1); call col2(w(n+1), mask2, n)
               
               rho0 = rho
               rho  = glsc3(w(1),p(1),mult,n) + glsc3(w(n+1),p(n+1),mult,n)
               alpha = rtz1/rho
               alphm = -alpha
               call add2s2(x1, p(1),   alpha, n)
               call add2s2(x2, p(n+1), alpha, n)
               call add2s2(r(1),   w(1),   alphm, n)
               call add2s2(r(n+1), w(n+1), alphm, n)
            enddo
            niter = iter-1
            if (nio.eq.0) write(6,3001) istep,'  CplHmh Error ' // name, niter, rbn2, rbn0, tol
 3000       format(i11,a,1x,I7,1p4E13.4)
 3001       format(i11,a,1x,I7,1p4E13.4)
 3002       format(i11,a,1x,I7,1p4E13.4)
 9999       continue
            niterhm = niter
            ifsolv = .false.
            call neklab_timer_stop(t_axisym_solve_cpl)
         end subroutine solve_coupled_helmholtz_2Dh_axisym
         
         subroutine solve_pressure_coupled_2Dh_axisym(resR, resI, h1, h2, h2inv, alphaR, alpha, ebar, intype, iter)
            !! Stacked GMRES for the coupled pressure system at alpha /= 0 AND
            !! lambda /= 0. Solves E [pR;pI] = [resR;resI] where E is the symmetric
            !! 2x2 Schur complement assembled by pressure_matvec_coupled_2Dh_axisym.
            !!
            !! Self-contained: it does NOT use the shared /cgmres1/, /cgmres2/... Nek
            !! GMRES commons (those are sized for a single ntot2 field and would be
            !! overrun by the 2*ntot2 stacked vector). Krylov storage is local.
            !!
            !! Preconditioner. hsmg_solve is a scalar single-field multigrid; applied
            !! BLOCK-DIAGONALLY, once per slot. The cross blocks are left to the
            !! Krylov iteration. This converges cheaply when lambda/alpha <~ 1, which
            !! covers the physical range (alpha = k R_c ~ O(1-10), lambda ~ O(0.1-5)).
            !!
            !! The alpha /= 0 branch never calls ortho (the constant-pressure mode is
            !! not in the null space when alpha /= 0), matching solve_pressure_2Dh.
            implicit none
            include 'GMRES'
            logical          ifprint
            common  /cprint/ ifprint
            integer, parameter :: n2 = lx2*ly2*lz2*lelv
            integer, parameter :: kmax = 200   ! restart length for the stacked solve
            real, dimension(n2), intent(inout) :: resR, resI
            real, dimension(lx1,ly1,lz1,lelv), intent(in) :: h1, h2, h2inv, alphaR
            real(dp), intent(in) :: alpha
            real(dp), intent(in) :: ebar
            integer, intent(in) :: intype
            integer, intent(inout) :: iter
            ! local Krylov storage (stacked: [R half ; I half])
            real, dimension(n2) :: xR, xI, rR, rI, wR, wI, zR, zI, tR, tI
            real, dimension(n2, kmax+1) :: vR, vI    ! store Arnoldi basis
            real, dimension(n2, kmax) :: zRk, zIk    ! store preconditioned basis for reconstruction
            real, dimension(kmax+1, kmax) :: hh
            real, dimension(kmax+1) :: gam
            real, dimension(kmax)   :: cc, ss
            real :: nrm, div0, tol_local, beta, alph, l, temp, norm_fac
            integer :: i, j, k, m, iconv, ntot2
            real, external :: glsc2, glsum
            real*8 :: etime1, dnekclock
            logical, save :: iflag_cpl = .false.
            
            call neklab_timer_start(t_axisym_solve_pres_cpl)
            
            if (.not. iflag_cpl) then
            iflag_cpl = .true.
            call uzawa_gmres_split0(ml_gmres, mu_gmres, bm2, bm2inv, n2)
            end if
            
            ntot2 = n2
            norm_fac = 1./sqrt(volvm2)
            etime1 = dnekclock()
            m = kmax
            iter = 0
            iconv = 0
            
            ! tolerance: mirror solve_pressure_2Dh's chktcg2-based choice
            call chktcg2(tolps, resR, i)
            if (param(21).gt.0.and.tolps.gt.abs(param(21))) tolps = abs(param(21))
            if (istep.eq.0) tolps = 1.e-4
            tol_local = tolps
            
            call rzero(xR, ntot2)
            call rzero(xI, ntot2)
            
            do while (iconv.eq.0 .and. iter.lt.2000)
               ! r = res - E x  (x = 0 on first pass)
               if (iter.eq.0) then
                  call col3(rR, ml_gmres, resR, ntot2)     ! r = L^-1 res
                  call col3(rI, ml_gmres, resI, ntot2)
               else
                  call pressure_matvec_coupled_2Dh_axisym(wR, wI, xR, xI, h1, h2, h2inv, alphaR, intype)
                  call copy(rR, resR, ntot2); call sub2(rR, wR, ntot2)
                  call copy(rI, resI, ntot2); call sub2(rI, wI, ntot2)
                  call col2(rR, ml_gmres, ntot2)           ! r = L^-1 r
                  call col2(rI, ml_gmres, ntot2)
               end if
               beta = sqrt(glsc2(rR, rR, ntot2) + glsc2(rI, rI, ntot2))
               if (iter.eq.0) div0 = beta*norm_fac
               if (beta.eq.0.) goto 900
               gam(1) = beta
               temp = 1./beta
               call cmult2(vR(1,1), rR, temp, ntot2)
               call cmult2(vI(1,1), rI, temp, ntot2)
               
               do j = 1, m
                  iter = iter + 1
                  ! block-diagonal preconditioner: hsmg on each slot separately
                  call neklab_timer_start(t_axisym_pres_prec_cpl)
                  call col3(wR, mu_gmres, vR(1,j), ntot2)   ! w = U^-1 v
                  call col3(wI, mu_gmres, vI(1,j), ntot2)
                  call hsmg_solve(zR, wR)     ! z = M^-1 w
                  call hsmg_solve(zI, wI)
                  if (ebar /= 0.0_dp) then
                     call cadd(zR, glsum(wR,ntot2)/ebar, ntot2)
                     call cadd(zI, glsum(wI,ntot2)/ebar, ntot2)
                  end if
                  call copy(zRk(1,j), zR, ntot2)   ! save z_j
                  call copy(zIk(1,j), zI, ntot2)
                  call neklab_timer_stop(t_axisym_pres_prec_cpl)
                  ! w = E z = E M^-1 v
                  call pressure_matvec_coupled_2Dh_axisym(wR, wI, zR, zI, h1, h2, h2inv, alphaR, intype)
                  call col2(wR, ml_gmres, ntot2)            ! w = L^-1 w
                  call col2(wI, ml_gmres, ntot2)
                  ! Arnoldi (modified Gram-Schmidt) over the stacked inner product
                  do i = 1, j
                     hh(i,j) = glsc2(wR, vR(1,i), ntot2) + glsc2(wI, vI(1,i), ntot2)
                     call add2s2(wR, vR(1,i), -hh(i,j), ntot2)
                     call add2s2(wI, vI(1,i), -hh(i,j), ntot2)
                  end do
                  alph = sqrt(glsc2(wR, wR, ntot2) + glsc2(wI, wI, ntot2))
                  hh(j+1,j) = alph
                  ! apply previous Givens rotations to column j
                  do i = 1, j-1
                     temp      =  cc(i)*hh(i,j) + ss(i)*hh(i+1,j)
                     hh(i+1,j) = -ss(i)*hh(i,j) + cc(i)*hh(i+1,j)
                     hh(i,j)   =  temp
                  end do
                  l = sqrt(hh(j,j)*hh(j,j) + alph*alph)
                  if (l.eq.0.) goto 900
                  temp = 1./l
                  cc(j) = hh(j,j)*temp
                  ss(j) = alph*temp
                  hh(j,j) = l
                  gam(j+1) = -ss(j)*gam(j)
                  gam(j)   =  cc(j)*gam(j)
                  nrm = abs(gam(j+1))*norm_fac
                  if (ifprint.and.nio.eq.0) write(6,66) iter, tol_local, nrm, div0, istep
   66             format(i5,1p3e12.5,i8,' P-cpl')
                  if (nrm.lt.tol_local) goto 900
                  if (j.eq.m) goto 800
                  if (alph.eq.0.) goto 900
                  temp = 1./alph
                  call cmult2(vR(1,j+1), wR, temp, ntot2)
                  call cmult2(vI(1,j+1), wI, temp, ntot2)
               end do
  800          continue
               ! back-substitute for the Krylov coefficients (store in cc)
               do k = j, 1, -1
                  temp = gam(k)
                  do i = j, k+1, -1
                     temp = temp - hh(k,i)*cc(i)
                  end do
                  cc(k) = temp/hh(k,k)
               end do
               ! form the solution update from the preconditioned vectors
               do i = 1, j
                  call add2s2(xR, zRk(1,i), cc(i), ntot2)   ! x += c_i z_i
                  call add2s2(xI, zIk(1,i), cc(i), ntot2)
               end do
               if (nrm.lt.tol_local) iconv = 1
            end do
  900       continue
            call copy(resR, xR, ntot2)
            call copy(resI, xI, ntot2)
            etime1 = dnekclock()-etime1
            if (nio.eq.0) write(6,9998) istep, '  U-PRES-cpl gmres  alpha= ', alpha, iter, nrm, div0, tol_local, etime1
 9998       format(i11,a,F5.1,1X,I6,1p4e13.4)
            call neklab_timer_stop(t_axisym_solve_pres_cpl)
         end subroutine solve_pressure_coupled_2Dh_axisym

         subroutine check_full_section_axisym(caller)
      !! Torsion breaks the Z -> -Z symmetry (it maps a helix to its mirror
      !! image), and dtheta is wall-normal on a plane Z = z_c. A helix therefore
      !! needs the full cross-section: no SYM or ON face.
            character(len=*), intent(in) :: caller
            character(len=3) :: cb
            integer :: ie, ifc, nbad
            integer, external :: iglsum
            nbad = 0
            do ie = 1, nelv
               do ifc = 1, 2*ndim
                  cb = cbc(ifc, ie, 1)
                  if (cb == 'SYM' .or. cb == 'ON ' .or. cb == 'on ') nbad = nbad + 1
               end do
            end do
            nbad = iglsum(nbad, 1)
            if (nbad > 0) call nek_stop_error('lambda /= 0 (helix) requires the full '//
     &         'cross-section: SYM/ON faces found.', this_module, caller)
         end subroutine check_full_section_axisym
         
      end module neklab_2Dh_axisym