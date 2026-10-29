      module neklab_2Dh_axisym
      !! Axisymmetric 2Dh linear stepper stages (torus, helix).
      !!
      !! Torsion (helix, lambda /= 0) enters through neklab_2Dh_torsion only, at
      !! five places, each the exact linearisation of the nonlinear step
      !! (neklab_2Dh%nek_advance_2Dh_axisym_nl):
      !!   makefp_2Dh_axisym         G2 + G3 (advection/frame, viscous), before EXT3
      !!   compute_gradp_axisym      +D_lam^T p'_ext in the u_s right-hand side
      !!   compute_frc_div_axisym    D_lam u'_s in the divergence
      !!   pressure_matvec_2Dh_axisym  + D_lam (h2inv B^-1) D_lam^T
      !!   compute_dw_axisym         u_s correction h2inv B^-1 D_lam^T dp
      !! At lambda = 0 every torsion call returns at once (torus unchanged).
      !! alpha /= 0 with lambda /= 0 is stage 6 and stops here.
         use LightKrylov, only: dp
         use neklab_2Dh_common, only: wmask, build_wmask,
     &                                check_2Dh_properties, check_2Dh_masks,
     &                                bm2_map12_2Dh, bm2_map12T_2Dh,
     &                                solve_helmholtz_2Dh
         use neklab_2Dh_torsion, only: init_torsion_2Dh, torsion_active,
     &                                 dlam_div_add, dlam_grad_t, lambda_block_E_add,
     &                                 torsion_add_linear_uv, torsion_add_linear_s
         use neklab_nek_setup, only: nek_log_debug, nek_log_message, nek_stop_error
         use neklab_timing, only: neklab_timer_start, neklab_timer_stop,
     &                            t_axisym_makefp, t_axisym_solve_cpl,
     &                            t_axisym_helm_matvec, t_axisym_cpl_matvec,
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
      ! diag_shift is S_PHI = (alpha^2+1)/R^2, diag_shift_r the u_R shift. They are
      ! equal: the torsion shifts lambda^2/R^2 (Z, R) are explicit (G3).
         real(dp), dimension(lx1,ly1,lz1,lelv), public, target :: diag_shift_r
         real(dp), dimension(lx1,ly1,lz1,lelv), public  :: alphaR_coef   ! alpha/R
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
         public ::          pres_matvec_axisym
 
         public :: build_axisym_coeffs
         public :: enable_axisym_curvature_terms
 
      contains

      !-----------------------------------------------------------------!
      !                                                                 !
      !     PUBLIC STAGE WRAPPERS FOR THE AXISYMMETRIC 2Dh SOLVER       !
      !                                                                 !
      !-----------------------------------------------------------------!
 
         subroutine init_2Dh_axisym(verbose)
      !! Set up the axisymmetric 2Dh solver. Idempotent. Call from userchk
      !! (istep = 0) after usrdat2; the nonlinear and linear steppers call it
      !! too. For a torus/helix, t2Dh%init_geom must come first: lambda is read
      !! here, once.
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
 
            call init_torsion_2Dh()
 
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
            call init_2Dh_axisym()
            if (.not. axisym_coeffs_defined .or. alpha /= alpha_cached) call build_axisym_coeffs(alpha)
            if (ifadj) then
               if (torsion_active()) call nek_stop_error(
     &            'adjoint mode is not implemented with torsion (lambda /= 0).',
     &            this_module, this_procedure)
               call gradm1(dTdx, dTdy, dTdz, t)          ! grad(U_phi)
               call gradm1(gUZ_Z, gUZ_R, tmp_m1c, vx)    ! grad(U_Z)
               call gradm1(gUR_Z, gUR_R, tmp_m1c, vy)    ! grad(U_R)
            end if
            ! Re/Im pairing is only meaningful for alpha /= 0; npert must match
            ! or the caller has set up the wrong number of perturbation slots.
            alpha_eq_zero = (alpha == 0.0_dp)
            if (.not. alpha_eq_zero .and. torsion_active()) call nek_stop_error(
     &         'alpha /= 0 with lambda /= 0 is not available yet (stage 6).',
     &         this_module, this_procedure)
            if (alpha_eq_zero) then
               if (lpert /= 1) call nek_stop_error('alpha = 0 requires npert = 1.', this_module, this_procedure)
            else
               if (lpert /= 2) call nek_stop_error('alpha /= 0 requires npert = 2.', this_module, this_procedure)
            end if
         end subroutine prepare_step_2Dh_axisym
 
         subroutine make_coupling_2Dh_axisym(gradw)
      !! Slot-coupled explicit terms of slot jp (read jp, write ipert): the
      !! third-component pressure gradient and the swirl advection
      !! -i alpha (U_phi/R) u. Runs for all slots before any slot is advanced.
            real(dp), dimension(lv,lpert), intent(inout) :: gradw
            call extrapprp(tmp_m2)
            call compute_gradp_axisym(gradw, tmp_m2, alphaR_coef)
            call compute_swirl_advection_axisym(advZ, advR, advPhi, jp)
         end subroutine make_coupling_2Dh_axisym
 
         subroutine ophx_2Dh_axisym(resv_jp, h1, h2)
      !! resv -= H u for all slots. u_Z per slot; u_R of slot jp couples to
      !! u_phi of the other slot at alpha /= 0.
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
      !! Right-hand side of the pressure step of slot jp:
      !! -(D u + D_alpha u_phi + D_lam u_phi). The alpha part reads the other slot.
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
      !! loop jp = 1, npert (zeroed at jp = 1):
      !!   alpha part: from slot jp into slot ipert (skew, slot swap), exact
      !!               transpose of the alpha block of compute_frc_div_axisym;
      !!   lambda part: +D_lam^T p_ext into slot jp, with the same extrapolated
      !!               pressure as the in-plane gradient (cresvipp_2Dh).
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
 
            if (torsion_active()) then
               call dlam_grad_t(tmp_m1a, prextr)
               call add2       (gradp(1,jp), tmp_m1a, ntot1)
            end if
 
            if (alpha_eq_zero) return
 
            spert = merge(1, -1, jp == 1)
            call bm2_map12T_2Dh(tmp_m1a, prextr)          ! I12^T (bm2 p)
            call col2          (tmp_m1a, alphaR, ntot1)
            if (spert == -1) call chsign(tmp_m1a, ntot1)
            call add2          (gradp(1,ipert), tmp_m1a, ntot1)
         end subroutine compute_gradp_axisym
 
         subroutine compute_dw_axisym(dw, dpr, h2inv, alphaR)
      !! u_phi correction of slot jp: mask h2inv B^-1 (alpha part + D_lam^T dp_jp).
            real(dp), dimension(lv), intent(out) :: dw
            real(dp), dimension(lp,lpert), intent(in) :: dpr
            real(dp), dimension(lx1,ly1,lz1,lelv), intent(in) :: h2inv, alphaR
            ! internal
            integer :: ipert, spert, ntot1
 
            ntot1 = lx1*ly1*lz1*nelv
            call rzero(dw, ntot1)
            if (alpha_eq_zero .and. .not. torsion_active()) return
 
            if (.not. alpha_eq_zero) then
               ipert = npert + 1 - jp
               spert = merge(1, -1, jp == 1)
               call bm2_map12T_2Dh(dw, dpr(1,ipert))      ! I12^T (bm2 dp)
               call col2          (dw, alphaR, ntot1)
               if (spert == 1) call chsign(dw, ntot1)     ! matches -spert*beta_z sign convention
            end if
            if (torsion_active()) then
               call dlam_grad_t(tmp_m1a, dpr(1,jp))
               call add2       (dw, tmp_m1a, ntot1)
            end if
 
            call col2 (dw, wmask, ntot1)
            call dssum(dw, lx1, ly1, lz1)
            call col2 (dw, binvm1, ntot1)
            call col2 (dw, h2inv, ntot1)
         end subroutine compute_dw_axisym
 
         subroutine compute_frc_div_axisym(frc_div, w, alphaR)
      !! frc_div = D_alpha u_phi (other slot) + D_lam u_phi (this slot).
            real(dp), dimension(lp), intent(out) :: frc_div
            real(dp), dimension(lv,ldimt,lpert), intent(in) :: w
            real(dp), dimension(lx1,ly1,lz1,lelv), intent(in) :: alphaR
            integer :: ipert, spert, ntot1, ntot2
 
            ntot1 = lx1*ly1*lz1*nelv
            ntot2 = lx2*ly2*lz2*nelv
            call rzero(frc_div, ntot2)
            if (alpha_eq_zero .and. .not. torsion_active()) return
 
            if (.not. alpha_eq_zero) then
               ipert = npert + 1 - jp
               spert = merge(1, -1, jp == 1)
               call copy(tmp_m1a, w(1,1,ipert), ntot1)
               call col2(tmp_m1a, alphaR, ntot1)       ! scale on M1 BEFORE mapping to M2
               if (spert == -1) call chsign(tmp_m1a, ntot1)
               call bm2_map12_2Dh(frc_div, tmp_m1a)
            end if
            call dlam_div_add(frc_div, w(1,1,jp))
         end subroutine compute_frc_div_axisym
 
         subroutine makefp_2Dh_axisym()
      !! Perturbation RHS of slot jp, replacing core makefp + makeqp: the swirl,
      !! curvature and torsion explicit terms are injected after the convective
      !! step and before EXT3 (makextp/makeabqp), as in the nonlinear step.
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
            call torsion_add_linear_uv(bfxp(1,jp), bfyp(1,jp), jp)   ! G2 + G3 (Z, R)
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
            call torsion_add_linear_s(bqp(1,ifield-1,jp))             ! G2 + G3 (s)
            if (iftran) call makeabqp                            ! <-- EXT3
            if ((iftran .and. .not. ifchar) .or.
     &          (iftran .and. .not. ifadvc(ifield) .and. ifchar)) call makebdqp
            call lagscalp
 
            ifield = ifield_bak
            call neklab_timer_stop(t_axisym_makefp)
         end subroutine makefp_2Dh_axisym
 
         subroutine add_axisym_curv_R(rhs, jp_)
      !! R-momentum: + 2*U_phi*u'_phi/R * rho * bm1 (linearisation of +U_phi^2/R)
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
      !! swirl equation: - (U_phi*u'_R + U_R*u'_phi)/R * rhocp * bm1
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
            call invcol2(tmp_m1a, ym1, ntot1)
            call chsign (tmp_m1a, ntot1)
            call col2   (tmp_m1a, vtrans(1,1,1,1,2), ntot1)
            call col2   (tmp_m1a, bm1, ntot1)
            call add2   (rhs, tmp_m1a, ntot1)
         end subroutine add_axisym_curv_phi
 
         subroutine compute_swirl_advection_axisym(advZ, advR, advPhi, jp_)
      !! Swirl advection -i*alpha*(U_phi/R)*u in all three momentum equations.
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
      !! Curvature terms of (U.grad)U missing from Nek's axisymmetric advection,
      !! for the nonlinear solver: R: +U_phi^2/R, phi: -U_R U_phi/R. Idempotent.
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
      !! E p = D (h2inv B^-1) D^T p with D = [D_in, D_alpha + D_lam].
            real(dp), dimension(lx2,ly2,lz2,lelv), intent(out) :: Ap
            real(dp), dimension(lx2,ly2,lz2,lelv), intent(in)  :: wp
            real(dp), dimension(lx1,ly1,lz1,lelv), intent(in)  :: h1, h2, h2inv, alphaR
            integer, intent(in) :: intype
            ! internal
            integer :: ntot1, ntot2
 
            call neklab_timer_start(t_axisym_pres_matvec)
            call cdabdtp(Ap, wp, h1, h2, h2inv, intype)   ! core Nek, already ifaxis-correct globally
            call lambda_block_E_add(ap, wp, h2inv, wmask)       ! no-op without torsion
            if (alpha_eq_zero) then
               call neklab_timer_stop(t_axisym_pres_matvec)
               return
            end if
            
            ntot1 = lx1*ly1*lz1*nelv
            ntot2 = lx2*ly2*lz2*nelv
      ! alpha block: (alpha/R) H^-1 (alpha/R), diagonal in slot space
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
            call neklab_timer_stop(t_axisym_pres_matvec)
         end subroutine pressure_matvec_2Dh_axisym
 
         subroutine coupled_helmholtz_matvec_2Dh_axisym(y1, y2, x1, x2, h1, h2, shift1, shift2, couple_coef, csign)
      !! [ D1  csign*c ; csign*c  D2 ] [x1; x2],  x1 = u_R, x2 = u_phi, c = 2 alpha/R^2
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
      !! CG on the stacked (u_R; u_phi) vector at alpha /= 0; separate masks.
            include 'FDMH1'
            real, intent(out) :: x1(1), x2(1)
            real, intent(inout)  :: f1(1), f2(1)
            real, intent(in)  :: h1(1), h2(1), mask1(1), mask2(1), mult(1), binv(1)
            integer, intent(in) :: imsh, maxit
            real, intent(in) :: tin
            character(len=7), intent(in) :: name
            real, intent(in) :: shift1(1), shift2(1), couple_coef(1)
            real(dp), intent(in) :: csign
            integer :: iter, n, nel, niter, nxyz
            real :: alpha, beta, alphm, fmax
            real :: rbn0, rbn2, rho, rho0, rtz1, rtz2
            real :: tol, vol
            LOGICAL IFPRINT, IFHZPC
            COMMON /CPRINT/ IFPRINT, IFHZPC
            logical ifdfrm, iffast, ifh2, ifsolv
            common /fastmd/ ifdfrm(lelt), iffast(lelt), ifh2, ifsolv
            logical ifprint_helm
            integer, parameter :: lg = lx1*ly1*lz1*lelt
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
            n = NEL*NXYZ
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
 
      end module neklab_2Dh_axisym