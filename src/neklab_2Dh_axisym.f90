      module neklab_2Dh_axisym
      !! Axisymmetric 2Dh linear stepper stages (torus, helix).
      !!
      !! Torsion (helix, lambda /= 0) enters through neklab_2Dh_torsion only, at
      !! five places, each the exact linearisation of the nonlinear step
      !! (neklab_2Dh%nek_advance_2Dh_axisym_nl):
      !!   makefp_2Dh_axisym           G2 + G3 (advection/frame, viscous), before EXT3
      !!   compute_gradp_axisym        +D_lam^T p'_ext in the u_s right-hand side
      !!   compute_frc_div_axisym      D_lam u'_s in the divergence
      !!   pressure_matvec_2Dh_axisym  + D_lam (h2inv B^-1) D_lam^T
      !!   compute_dw_axisym           u_s correction h2inv B^-1 D_lam^T dp
      !! At lambda = 0 every torsion call returns at once (torus unchanged).
      !!
      !! alpha /= 0 with lambda /= 0: the alpha*lambda viscous coupling G6 is
      !! added by torsion_add_linear_uv (slot swap), and the pressure system
      !! couples the slots: E = blockdiag(E_in) + F M F^T with
      !! F = [[D_lam, C], [-C, D_lam]] on (u_s Re, u_s Im) (C: the alpha block).
      !! It is assembled from the same divergence (compute_frc_div_axisym) and
      !! correction (compute_dw_axisym) routines that the step uses, so the
      !! corrected field is divergence-free by construction, and solved by the
      !! block GMRES of neklab_2Dh_common (incomprp_block_2Dh_axisym).
      !!
      !! Adjoint (ifadj, same alpha): Nek's adjoint transport (advabp_adjoint_axisym,
      !! convabp), the swirl advection with its sign flipped, the transposed
      !! curvature block, and from neklab_2Dh_torsion the exact transposes of
      !! G2/G3/G6 plus the (div_in U) v term. The D_lam pieces are part of the
      !! self-adjoint projection and unchanged.
         use LightKrylov, only: dp
         use neklab_2Dh_common, only: wmask, build_wmask,
     &                                check_2Dh_properties, check_2Dh_masks,
     &                                bm2_map12_2Dh, bm2_map12T_2Dh,
     &                                solve_helmholtz_2Dh, incomprp_block_solve_2Dh
         use neklab_2Dh_torsion, only: init_torsion_2Dh, lambda_eq_zero,
     &                                 add_dlam_div, dlam_grad, add_lambda_block_E,
     &                                 torsion_add_linear_uv, torsion_add_linear_s,
     &                                 set_torsion_alpha, torsion_prepare_adjoint
         use neklab_nek_setup, only: nek_log_message, nek_stop_error
         use neklab_timing, only: neklab_timer_start, neklab_timer_stop,
     &                            t_axisym_makefp, t_axisym_solve_cpl,
     &                            t_axisym_helm_matvec, t_axisym_cpl_matvec,
     &                            t_axisym_pres_matvec, t_axisym_pres_block_matvec
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
         real(dp), dimension(lv),       private :: tmp_m1a, tmp_m1b
         real(dp), dimension(lv),       private :: tmp_m1c
         real(dp), dimension(lp),       private :: tmp_m2
         real(dp), dimension(lv,ldimt,lpert), private :: blk_us    ! block E: u_s = M G p
         real(dp), dimension(lp,2),     private :: blk_one, blk_e1  ! block ebar

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
         public ::   incomprp_block_2Dh_axisym
         public ::          pres_matvec_axisym
         public ::    pres_block_matvec_axisym

         public :: build_axisym_coeffs
         public :: enable_axisym_curvature_terms

      contains

         subroutine init_2Dh_axisym(verbose)
      !! Set up the axisymmetric 2Dh solver. Idempotent. Call from userchk
      !! (istep = 0) after usrdat2; the nonlinear and linear steppers call it
      !! too. For a torus/helix, t2Dh%init_geom must come first: lambda is read
      !! here, once.
            logical, optional, intent(in) :: verbose
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
            character(len=*), parameter :: this_procedure = 'prepare_step_2Dh_axisym'
            call init_2Dh_axisym()
            if (.not. axisym_coeffs_defined .or. alpha /= alpha_cached) call build_axisym_coeffs(alpha)
            if (ifadj) then
               call gradm1(dTdx, dTdy, dTdz, t)          ! grad(U_phi)
               call gradm1(gUZ_Z, gUZ_R, tmp_m1c, vx)    ! grad(U_Z)
               call gradm1(gUR_Z, gUR_R, tmp_m1c, vy)    ! grad(U_R)
               call torsion_prepare_adjoint()            ! div_in U; no-op at lambda = 0
            end if
            alpha_eq_zero = (alpha == 0.0_dp)
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
            character(len=2) :: info_str
            character(len=7) :: helm_info
            integer :: ntot1
            ntot1 = lx1*ly1*lz1*nelv
            do jp = 1, npert
               info_str = merge('Re', 'Im', jp == 1)
               call dssum(resv_jp(1,1,jp), lx1, ly1, lz1)
               call col2 (resv_jp(1,1,jp), v1mask, ntot1)
               helm_info = info_str//' V_Z '
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
                  helm_info = info_str//' V_R '
                  if (istep < 10) call chktcg1(tolhv, resv_jp(1,2,jp), h1, h2, v2mask, vmult, imesh, 1)
                  shift_op => diag_shift_r
                  call solve_helmholtz_2Dh(tmp_m1a, resv_jp(1,2,jp), h1, h2, v2mask, vmult, imesh,
     &                                     tolhv, nmxv, 1, binvm1, helm_info, helm_matvec_axisym)

                  call dssum(resv_jp(1,3,jp), lx1, ly1, lz1)
                  call col2 (resv_jp(1,3,jp), wmask, ntot1)
                  helm_info = info_str//' VPHI'
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
            real(dp), intent(in) :: alpha
            integer :: ntot1
            ntot1 = lx1*ly1*lz1*nelv
            call invers2(alphaR_coef, ym1, ntot1)                      ! 1/R
            call col3   (diag_shift, alphaR_coef, alphaR_coef, ntot1)  ! 1/R^2
            call cmult2 (h2z_shift, diag_shift, alpha**2, ntot1)       ! alpha^2/R^2
            call cmult2 (couple_coef, diag_shift, 2.0_dp*alpha, ntot1) ! 2 alpha/R^2
            call cmult  (diag_shift, 1.0_dp + alpha**2, ntot1)         ! (alpha^2+1)/R^2
            call copy   (diag_shift_r, diag_shift, ntot1)
            call cmult  (alphaR_coef, alpha, ntot1)                    ! alpha/R
            call set_torsion_alpha(alpha)       ! G6; no-op without torsion
            alpha_cached          = alpha
            axisym_coeffs_defined  = .true.
         end subroutine build_axisym_coeffs

         subroutine helmholtz_matvec_2Dh_axisym(Au, u, h1, h2, shift, isd)
      !! Au = (h1*A + h2*B) u + vdiff*shift*bm1*u
            real(dp), dimension(lv), intent(out) :: Au
            real(dp), dimension(lv), intent(in)  :: u, h1, h2, shift
            integer, intent(in) :: isd
            integer :: ntot1
            call neklab_timer_start(t_axisym_helm_matvec)
            ntot1 = lx1*ly1*lz1*nelv
            call axhelm (Au, u, h1, h2, 1, isd)
            call col3   (tmp_m1c, u, vdiff(1,1,1,1,1), ntot1)
            call col2   (tmp_m1c, shift, ntot1)
            call addcol3(Au, tmp_m1c, bm1, ntot1)
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
            real(dp), dimension(lv), intent(in) :: alphaR
            integer :: ipert, kp, ntot1

            ntot1 = lx1*ly1*lz1*nelv
            ipert = npert + 1 - jp

            if (jp == 1) then
               do kp = 1, npert
                  call rzero(gradp(1,kp), ntot1)
               end do
            end if

            if (.not. lambda_eq_zero) then
               call dlam_grad(tmp_m1a, prextr)
               call add2     (gradp(1,jp), tmp_m1a, ntot1)
            end if

            if (alpha_eq_zero) return

            call bm2_map12T_2Dh(tmp_m1a, prextr)          ! I12^T (bm2 p)
            call col2          (tmp_m1a, alphaR, ntot1)
            if (jp == 2) call chsign(tmp_m1a, ntot1)
            call add2          (gradp(1,ipert), tmp_m1a, ntot1)
         end subroutine compute_gradp_axisym

         subroutine compute_dw_axisym(dw, dpr, h2inv, alphaR)
      !! u_phi correction of slot jp: mask h2inv B^-1 (alpha part + D_lam^T dp_jp).
            real(dp), dimension(lv), intent(out) :: dw
            real(dp), dimension(lp,lpert), intent(in) :: dpr
            real(dp), dimension(lv), intent(in) :: h2inv, alphaR
            integer :: ipert, ntot1

            ntot1 = lx1*ly1*lz1*nelv
            call rzero(dw, ntot1)
            if (alpha_eq_zero .and. lambda_eq_zero) return

            if (.not. alpha_eq_zero) then
               ipert = npert + 1 - jp
               call bm2_map12T_2Dh(dw, dpr(1,ipert))      ! I12^T (bm2 dp)
               call col2          (dw, alphaR, ntot1)
               if (jp == 1) call chsign(dw, ntot1)
            end if
            if (.not. lambda_eq_zero) then
               call dlam_grad(tmp_m1a, dpr(1,jp))
               call add2     (dw, tmp_m1a, ntot1)
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
            real(dp), dimension(lv), intent(in) :: alphaR
            integer :: ipert, ntot1

            ntot1 = lx1*ly1*lz1*nelv
            call rzero(frc_div, lx2*ly2*lz2*nelv)
            if (alpha_eq_zero .and. lambda_eq_zero) return

            if (.not. alpha_eq_zero) then
               ipert = npert + 1 - jp
               call col3(tmp_m1a, w(1,1,ipert), alphaR, ntot1)   ! scale on M1 BEFORE mapping to M2
               if (jp == 2) call chsign(tmp_m1a, ntot1)
               call bm2_map12_2Dh(frc_div, tmp_m1a)
            end if
            call add_dlam_div(frc_div, w(1,1,jp))
         end subroutine compute_frc_div_axisym

         subroutine makefp_2Dh_axisym()
      !! Perturbation RHS of slot jp, replacing core makefp + makeqp: the swirl,
      !! curvature and torsion explicit terms are injected after the convective
      !! step and before EXT3 (makextp/makeabqp), as in the nonlinear step.
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
                  call advabp_adjoint_axisym
               else
                  call advabp
               end if
            end if
            call add2(bfxp(1,jp), advZ(1,jp), ntot1)
            call add2(bfyp(1,jp), advR(1,jp), ntot1)
            call add_axisym_curv_R(bfyp(1,jp), jp)
            call torsion_add_linear_uv(bfxp(1,jp), bfyp(1,jp), jp)   ! torsion (Z, R)
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
            call add2(bqp(1,ifield-1,jp), advPhi(1,jp), ntot1)
            call add_axisym_curv_phi(bqp(1,ifield-1,jp), jp)
            call torsion_add_linear_s(bqp(1,ifield-1,jp))             ! torsion (s)
            if (iftran) call makeabqp                            ! <-- EXT3
            if ((iftran .and. .not. ifchar) .or.
     &          (iftran .and. .not. ifadvc(ifield) .and. ifchar)) call makebdqp
            call lagscalp

            ifield = ifield_bak
            call neklab_timer_stop(t_axisym_makefp)
         end subroutine makefp_2Dh_axisym

         subroutine add_axisym_curv_R(rhs, jp_)
      !! R-momentum, rho*bm1 times:  direct  +2 U_phi u_phi/R  (linearisation of +U_phi^2/R)
      !!                             adjoint  -U_phi u_phi/R   (transposed curvature block)
            real(dp), dimension(lv), intent(inout) :: rhs
            integer, intent(in) :: jp_
            integer :: ntot1
            ntot1 = lx1*ly1*lz1*nelv
            call col3   (tmp_m1a, t(1,1,1,1,1), tp(1,1,jp_), ntot1)
            call invcol2(tmp_m1a, ym1, ntot1)
            call cmult  (tmp_m1a, merge(-1.0_dp, 2.0_dp, ifadj), ntot1)
            call col2   (tmp_m1a, vtrans(1,1,1,1,1), ntot1)
            call addcol3(rhs, tmp_m1a, bm1, ntot1)
         end subroutine add_axisym_curv_R

         subroutine advabp_adjoint_axisym()
      !! Adjoint transport of the u_Z, u_R equations, slot jp:
      !!   + (U.grad) v  -  (grad U)^T v  -  v_phi grad(U_phi)
            integer :: ntot1
            ntot1 = lx1*ly1*lz1*nelv
            call convop  (tmp_m1a, vxp(1,jp))
            call subcol3 (tmp_m1a, vxp(1,jp),  gUZ_Z, ntot1)
            call subcol3 (tmp_m1a, vyp(1,jp),  gUR_Z, ntot1)
            call subcol3 (tmp_m1a, tp(1,1,jp), dTdx,  ntot1)
            call col2    (tmp_m1a, vtrans(1,1,1,1,1), ntot1)
            call addcol3 (bfxp(1,jp), tmp_m1a, bm1, ntot1)
            call convop  (tmp_m1b, vyp(1,jp))
            call subcol3 (tmp_m1b, vxp(1,jp),  gUZ_R, ntot1)
            call subcol3 (tmp_m1b, vyp(1,jp),  gUR_R, ntot1)
            call subcol3 (tmp_m1b, tp(1,1,jp), dTdy,  ntot1)
            call col2    (tmp_m1b, vtrans(1,1,1,1,1), ntot1)
            call addcol3 (bfyp(1,jp), tmp_m1b, bm1, ntot1)
         end subroutine advabp_adjoint_axisym

         subroutine add_axisym_curv_phi(rhs, jp_)
      !! swirl equation, rho cp*bm1 times:  direct  -(U_phi u_R + U_R u_phi)/R
      !!                                     adjoint  (2 U_phi u_R - U_R u_phi)/R
            real(dp), dimension(lv), intent(inout) :: rhs
            integer, intent(in) :: jp_
            integer :: ntot1
            ntot1 = lx1*ly1*lz1*nelv
            call col3   (tmp_m1a, t(1,1,1,1,1), vyp(1,jp_),  ntot1)
            if (ifadj) call cmult(tmp_m1a, -2.0_dp, ntot1)
            call addcol3(tmp_m1a, vy, tp(1,1,jp_), ntot1)
            call invcol2(tmp_m1a, ym1, ntot1)
            call chsign (tmp_m1a, ntot1)
            call col2   (tmp_m1a, vtrans(1,1,1,1,2), ntot1)
            call addcol3(rhs, tmp_m1a, bm1, ntot1)
         end subroutine add_axisym_curv_phi

         subroutine compute_swirl_advection_axisym(advZ, advR, advPhi, jp_)
      !! Swirl advection -i*alpha*(U_phi/R)*u in all three momentum equations
      !! (adjoint: sign flipped, antisymmetric in slot space).
            real(dp), dimension(lv,lpert), intent(inout) :: advZ, advR, advPhi
            integer, intent(in) :: jp_
            integer :: ipert, ntot1
            ntot1 = lx1*ly1*lz1*nelv
            ipert = npert + 1 - jp_
            call rzero(advZ(1,ipert),   ntot1)
            call rzero(advR(1,ipert),   ntot1)
            call rzero(advPhi(1,ipert), ntot1)
            if (alpha_eq_zero) return
            call col3(tmp_m1a, alphaR_coef, t(1,1,1,1,1), ntot1)
            call col2(tmp_m1a, bm1, ntot1)
            if (ifadj) call chsign(tmp_m1a, ntot1)
            if (ipert == 1) call chsign(tmp_m1a, ntot1)
            call col3(tmp_m1b, tmp_m1a, vtrans(1,1,1,1,1), ntot1)
            call col3(advZ(1,ipert),   tmp_m1b, vxp(1,jp_),  ntot1)
            call col3(advR(1,ipert),   tmp_m1b, vyp(1,jp_),  ntot1)
            call col2(tmp_m1a, vtrans(1,1,1,1,2), ntot1)
            call col3(advPhi(1,ipert), tmp_m1a, tp(1,1,jp_), ntot1)
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
            real(dp), dimension(lp), intent(out) :: ap
            real(dp), dimension(lp), intent(in)  :: wp
            real(dp), dimension(lv), intent(in)  :: h1, h2, h2inv, alphaR
            integer, intent(in) :: intype
            integer :: ntot1

            call neklab_timer_start(t_axisym_pres_matvec)
            call cdabdtp(ap, wp, h1, h2, h2inv, intype)
            call add_lambda_block_E(ap, wp, h2inv, wmask)       ! no-op without torsion
            if (alpha_eq_zero) then
               call neklab_timer_stop(t_axisym_pres_matvec)
               return
            end if
            ntot1 = lx1*ly1*lz1*nelv
      ! alpha block: (alpha/R) H^-1 (alpha/R), diagonal in slot space
            call bm2_map12T_2Dh(tmp_m1a, wp)
            call col2          (tmp_m1a, alphaR, ntot1)
            call col2          (tmp_m1a, wmask,  ntot1)
            call dssum         (tmp_m1a, lx1, ly1, lz1)
            call col2          (tmp_m1a, binvm1, ntot1)
            call col2          (tmp_m1a, h2inv,  ntot1)
            call col2          (tmp_m1a, alphaR, ntot1)
            call bm2_map12_2Dh (tmp_m2, tmp_m1a)
            call add2          (ap, tmp_m2, lx2*ly2*lz2*nelv)
            call neklab_timer_stop(t_axisym_pres_matvec)
         end subroutine pressure_matvec_2Dh_axisym

         subroutine coupled_helmholtz_matvec_2Dh_axisym(y1, y2, x1, x2, h1, h2, shift1, shift2, couple_coef, csign)
      !! [ D1  csign*c ; csign*c  D2 ] [x1; x2],  x1 = u_R, x2 = u_phi, c = 2 alpha/R^2
            real(dp), dimension(lv), intent(out) :: y1, y2
            real(dp), dimension(lv), intent(in)  :: x1, x2, h1, h2, shift1, shift2, couple_coef
            real(dp), intent(in) :: csign
            integer :: ntot1
            ntot1 = lx1*ly1*lz1*nelv
            call neklab_timer_start(t_axisym_cpl_matvec)
            call helmholtz_matvec_2Dh_axisym(y1, x1, h1, h2, shift1, 1)
            call helmholtz_matvec_2Dh_axisym(y2, x2, h1, h2, shift2, 1)
            call col3   (tmp_m1c, couple_coef, vdiff(1,1,1,1,1), ntot1)
            call col2   (tmp_m1c, bm1, ntot1)
            call cmult  (tmp_m1c, csign, ntot1)                ! csign * c * mu * bm1
            call addcol3(y1, tmp_m1c, x2, ntot1)
            call addcol3(y2, tmp_m1c, x1, ntot1)
            call neklab_timer_stop(t_axisym_cpl_matvec)
         end subroutine coupled_helmholtz_matvec_2Dh_axisym

         subroutine helm_matvec_axisym(Au, u, h1, h2)
            real(dp), dimension(lv), intent(out) :: Au
            real(dp), dimension(lv), intent(in)  :: u, h1, h2
            call helmholtz_matvec_2Dh_axisym(Au, u, h1, h2, shift_op, 1)
         end subroutine helm_matvec_axisym

         subroutine pres_block_matvec_axisym(Ap, p, h1, h2, h2inv)
      !! Two-slot pressure operator E p = D (h2inv B^-1) G p, G = D^T:
      !!   in-plane, per slot:  cdabdtp (Nek's D_in h2inv B^-1 D_in^T);
      !!   u_s part:            u_s = compute_dw (mask, dssum, B^-1, h2inv of G_s p),
      !!                        Ap += compute_frc_div (u_s),
      !! with exactly the operators of the step's correction and divergence.
            real(dp), dimension(lp,2), intent(out) :: Ap
            real(dp), dimension(lp,2), intent(in)  :: p
            real(dp), dimension(lv),   intent(in)  :: h1, h2, h2inv
            integer :: jp_save, k, ntot2
            call neklab_timer_start(t_axisym_pres_block_matvec)
            ntot2 = lx2*ly2*lz2*nelv
            jp_save = jp
            do k = 1, 2
               call cdabdtp(Ap(1,k), p(1,k), h1, h2, h2inv, 1)
            end do
            do jp = 1, 2
               call compute_dw_axisym(blk_us(1,1,jp), p, h2inv, alphaR_coef)
            end do
            do jp = 1, 2
               call compute_frc_div_axisym(tmp_m2, blk_us, alphaR_coef)
               call add2(Ap(1,jp), tmp_m2, ntot2)
            end do
            jp = jp_save
            call neklab_timer_stop(t_axisym_pres_block_matvec)
         end subroutine pres_block_matvec_axisym

         subroutine incomprp_block_2Dh_axisym(dpr, h1, h2, h2inv, alpha, maxit)
      !! Coupled pressure step of both slots: dpr(:,k) holds the right-hand
      !! sides -(D u + F u_s) on entry (incomprp_div_2Dh_axisym), the pressure
      !! corrections on exit. ebar (coarse correction of the hsmg preconditioner)
      !! from the slot-diagonal part: sum(E (1, 0))_Re = sum(E (0, 1))_Im.
            real(dp), dimension(lp,lpert), intent(inout) :: dpr
            real(dp), dimension(lv), intent(in) :: h1, h2, h2inv
            real(dp), intent(in) :: alpha
            integer,  intent(in) :: maxit
            real(dp) :: ebar
            real(dp), external :: glsum
            integer :: ntot2
            ntot2 = lx2*ly2*lz2*nelv
            call rzero(blk_one, 2*lp)
            call rone (blk_one(1,1), ntot2)
            call pres_block_matvec_axisym(blk_e1, blk_one, h1, h2, h2inv)
            ebar = glsum(blk_e1(1,1), ntot2)
            call incomprp_block_solve_2Dh(dpr, h1, h2, h2inv, pres_block_matvec_axisym,
     &                                    ebar, alpha, maxit, 'alpha')
         end subroutine incomprp_block_2Dh_axisym

         subroutine pres_matvec_axisym(Ap, p, h1, h2, h2inv, intype)
            real(dp), dimension(lp), intent(out) :: Ap
            real(dp), dimension(lp), intent(in)  :: p
            real(dp), dimension(lv), intent(in)  :: h1, h2, h2inv
            integer, intent(in) :: intype
            call pressure_matvec_2Dh_axisym(Ap, p, h1, h2, h2inv, alphaR_coef, intype)
         end subroutine pres_matvec_axisym

         subroutine solve_coupled_helmholtz_2Dh_axisym(x1,x2,f1,f2,h1,h2,mask1,mask2,mult,imsh,tin,maxit,
     &                                                 binv,name,shift1,shift2,couple_coef,csign)
      !! CG on the stacked (u_R; u_phi) vector at alpha /= 0; separate masks.
      !! Nek's cggo (hmholtz.f) with the coupled matvec: cggo hard-wires axhelm.
            real(dp), intent(out) :: x1(lv), x2(lv)
            real(dp), intent(inout)  :: f1(lv), f2(lv)
            real(dp), intent(in)  :: h1(lv), h2(lv), mask1(lv), mask2(lv), mult(lv), binv(lv)
            integer, intent(in) :: imsh, maxit
            real(dp), intent(in) :: tin
            character(len=7), intent(in) :: name
            real(dp), intent(in) :: shift1(lv), shift2(lv), couple_coef(lv)
            real(dp), intent(in) :: csign
            integer :: iter, n, niter
            real(dp) :: alpha, beta, alphm, fmax
            real(dp) :: rbn0, rbn2, rho, rho0, rtz1, rtz2
            real(dp) :: tol, vol
            logical          ifprint, ifhzpc
            common  /cprint/ ifprint, ifhzpc
            logical ifdfrm, iffast, ifh2, ifsolv
            common /fastmd/ ifdfrm(lelt), iffast(lelt), ifh2, ifsolv
            logical ifprint_helm
            integer, parameter :: lg = lx1*ly1*lz1*lelt
            real(dp), dimension(2*lg) :: d, r, w, p, z
            real(dp) :: scalar(2)
            integer :: niterhm
            common /iterhm/ niterhm
            real(dp), external :: glamax, glsc3, vlsc3, vlsc32

            call neklab_timer_start(t_axisym_solve_cpl)
            rho = 0.0_dp
            vol = volvm1
            n   = lx1*ly1*lz1*nelv
            tol = abs(tin)
            if (restol(ifield) /= 0) tol = abs(restol(ifield))
            if (tin < 0) tol = abs(tin)
            niter = min(maxit, 900)
            if (.not. ifsolv) then
               call setfast(h1, h2, imsh)
               ifsolv = .true.
            end if
            call setprec(d(1), h1, h2, imsh, 1)
            call copy   (d(n+1), d(1), n)
            call dssum(f1, lx1, ly1, lz1)
            call dssum(f2, lx1, ly1, lz1)
            call col3(r(1),   f1, mask1, n)
            call col3(r(n+1), f2, mask2, n)
            call rzero(x1, n); call rzero(x2, n)
            call rzero(p, 2*n)
            fmax = max(glamax(r(1), n), glamax(r(n+1), n))
            if (fmax == 0.0_dp) then
               call neklab_timer_stop(t_axisym_solve_cpl)
               return
            end if
            rtz1 = 1.0_dp
            niterhm = 0
            do iter = 1, niter
               call col3(z(1),   r(1),   d(1),   n)
               call col3(z(n+1), r(n+1), d(n+1), n)
               rtz2 = rtz1
               scalar(1) = vlsc3(z(1), r(1), mult, n) + vlsc3(z(n+1), r(n+1), mult, n)
               if (param(18) == 1) then
                  scalar(2) = vlsc3(r(1), r(1), mult, n) + vlsc3(r(n+1), r(n+1), mult, n)
               else
                  scalar(2) = vlsc32(r(1), mult, binv, n) + vlsc32(r(n+1), mult, binv, n)
               end if
               call gop(scalar, w, '+  ', 2)
               rtz1 = scalar(1)
               rbn2 = sqrt(scalar(2)/(2.0_dp*vol))
               if (iter == 1) rbn0 = rbn2
               if (param(22) < 0) tol = abs(param(22))*rbn0
               if (tin < 0)       tol = abs(tin)*rbn0
               ifprint_helm = .false.
               if (nio == 0 .and. ifprint .and. param(74) /= 0) ifprint_helm = .true.
               if (nio == 0 .and. istep == 1)                   ifprint_helm = .true.
               if (ifprint_helm) write(6,3002) istep, '  CplHmh '//name, iter, rbn2, h1(1), tol, h2(1)
               if (rbn2 <= tol .and. (iter > 1 .or. istep <= 5)) then
                  niter = iter - 1
                  if (nio == 0) write(6,3000) istep, '  CplHmh '//name, niter, rbn2, rbn0, tol
                  goto 9999
               end if
               beta = rtz1/rtz2
               if (iter == 1) beta = 0.0_dp
               call add2s1(p, z, beta, 2*n)
               call coupled_helmholtz_matvec_2Dh_axisym(w(1), w(n+1), p(1), p(n+1),
     &                                                  h1, h2, shift1, shift2, couple_coef, csign)
               call dssum(w(1),   lx1, ly1, lz1); call col2(w(1),   mask1, n)
               call dssum(w(n+1), lx1, ly1, lz1); call col2(w(n+1), mask2, n)
               rho0  = rho
               rho   = glsc3(w(1), p(1), mult, n) + glsc3(w(n+1), p(n+1), mult, n)
               alpha = rtz1/rho
               alphm = -alpha
               call add2s2(x1, p(1),   alpha, n)
               call add2s2(x2, p(n+1), alpha, n)
               call add2s2(r(1),   w(1),   alphm, n)
               call add2s2(r(n+1), w(n+1), alphm, n)
            end do
            niter = iter - 1
            if (nio == 0) write(6,3001) istep, '  CplHmh Error '//name, niter, rbn2, rbn0, tol
 3000       format(i11,a,1x,I7,1p4E13.4)
 3001       format(i11,a,1x,I7,1p4E13.4)
 3002       format(i11,a,1x,I7,1p4E13.4)
 9999       continue
            niterhm = niter
            ifsolv = .false.
            call neklab_timer_stop(t_axisym_solve_cpl)
         end subroutine solve_coupled_helmholtz_2Dh_axisym

      end module neklab_2Dh_axisym