      module neklab_2Dh_torsion
      !! Torsion terms for helical pipes in the Frenet frame of the centreline
      !! (derivation and symbolic verification: claude/helix_2Dh_derivation.tex).
      !!
      !!--------------------------------------------------------------------
      !! OPERATORS
      !!--------------------------------------------------------------------
      !!
      !!   streamwise derivative   T.grad = (1/R)(d_phi - lambda*dtheta)
      !!   rotation generator      dtheta = zhat*d/dR - rhat*d/dZ ,
      !!                           zhat = xm1 - z_c, rhat = ym1 - R_c
      !!   rotation of a vector    S(v) = (dtheta v_Z + v_R, dtheta v_R - v_Z, dtheta v_phi)
      !!
      !! lambda = tau/kappa and the axial centre z_c are set by t2Dh%init_geom;
      !! lambda = 0 (default) is the torus and every entry point returns at once.
      !! dtheta is tangent to the wall only for a circular wall centred on
      !! (z_c, R_c): the mesh must be registered accordingly, and must be the
      !! FULL cross-section (torsion breaks the Z -> -Z symmetry).
      !!
      !!--------------------------------------------------------------------
      !! ROUTE A  (claude/neklab_2Dh_plan.md, section 3.1)
      !!--------------------------------------------------------------------
      !!
      !! Every torsion term is explicit, in the nonlinear AND in the linear
      !! solver, built from the same kernels below, at the same time level:
      !!
      !!   advection + frame rotation   RHS +(lambda U_phi/R) S(U)
      !!                                (linear: + c(U) S(u) + c(u) S(U))
      !!   viscous                      nu * [torsion part of the vector
      !!                                Laplacian], INCLUDING -lambda^2/R^2 on
      !!                                the Z and R components
      !!   pressure (u_phi equation)    +D_lam^T p   (weak, exact transpose of D_lam)
      !!   continuity                   D u + D_lam u_phi = 0,
      !!                                D_lam u = bm2 * I12( -(lambda/R) dtheta u )
      !!
      !! Nonlinear: refresh_torsion_forcing_nl (before every step, from U^n,
      !! P^n) feeds registered pointwise terms; set_torsion_usrdiv (after heat(2),
      !! with U_phi^{n+1}) sets Nek's usrdiv. Both are called by
      !! nek_advance_2Dh_axisym_nl.
      !!
      !! Linear (alpha = 0): add_torsion_forcing_RZ/_phi inject the same terms
      !! before EXT3; the continuity term stays in compute_frc_div_axisym; the
      !! pressure operator E and the u_phi correction carry NO lambda block.
      !!
      !! alpha /= 0 with lambda /= 0 is refused by prepare_step_2Dh_axisym until
      !! Step 7 (add_torsion_coupling and dphi_alpha_apply* are kept for it).
         use LightKrylov, only: dp
         use neklab_t2Dh, only: t2Dh
         use neklab_nek_setup, only: nek_log_debug, nek_log_message,
     &                               nek_log_warning, nek_stop_error
         use neklab_nek_forcing, only: add_neklab_forcing_term, add_neklab_source_term
         use neklab_2Dh_common, only: map12_2Dh, bm2_map12T_2Dh
         implicit none
         include "SIZE"
         include "TOTAL"

         private
         character(len=*), parameter, private :: this_module = 'neklab_2Dh_torsion'

         logical, public :: if_torsion_zero = .true.
      !! Mirrors t2Dh%is_torsion_zero(), cached at build_torsion_coeffs time.

         integer, parameter, private :: lv = lx1*ly1*lz1*lelv
         integer, parameter, private :: lxyz = lx1*ly1*lz1

      ! --- coefficient fields, built once (lambda does not change at runtime)
         real(dp), dimension(lv), private :: zhat_f, rhat_f
         real(dp), dimension(lv), private :: lamR_f     ! lambda/R
         real(dp), dimension(lv), private :: lamR2_f    ! lambda/R^2
         real(dp), dimension(lv), private :: lam2R2_f   ! lambda^2/R^2
         real(dp), dimension(lv), private :: lamzR3_f   ! lambda*zhat/R^3
         real(dp), dimension(lv), private :: invR_f     ! 1/R
         logical,  private :: torsion_coeffs_defined = .false.
         real(dp), private :: lambda_built = 0.0_dp     ! lambda the coefficients were built with

      ! --- nonlinear forcing (accelerations), refreshed before every step and
      !     read pointwise by the registered terms below
         real(dp), dimension(lv), private :: fZ_nl = 0.0_dp, fR_nl = 0.0_dp, fP_nl = 0.0_dp

      ! --- scratch. Each group belongs to one level of the call tree, so no
      !     routine overwrites an array its caller still reads.
         real(dp), dimension(lv), private :: g_x, g_y, g_z       ! apply_dtheta
         real(dp), dimension(lv), private :: l_g                 ! apply_Ltheta
         real(dp), dimension(lv), private :: c_t                 ! add_cxy
         real(dp), dimension(lv), private :: p_t                 ! dphi_apply_t
         real(dp), dimension(lv), private :: m_t                 ! dlam_grad_t
         real(dp), dimension(lv), private :: k_dZ, k_dR, k_dP, k_L ! kernels
         real(dp), dimension(lv), private :: acc_a, acc_b, acc_c, cf_1 ! drivers

         public :: build_torsion_coeffs
         public :: check_torsion_cache
         public :: apply_dtheta
         public :: apply_dtheta_t
         public :: dphi_apply
         public :: dphi_apply_t
         public :: dphi_alpha_apply
         public :: dphi_alpha_apply_t
         public :: add_torsion_forcing_RZ
         public :: add_torsion_forcing_phi
         public :: add_torsion_coupling
         public :: refresh_torsion_forcing_nl
         public :: set_torsion_usrdiv

      contains

      !====================================================================
      !     SETUP
      !====================================================================

         subroutine build_torsion_coeffs()
      !! Builds the coefficient fields from t2Dh%init_geom and, for lambda /= 0,
      !! registers the nonlinear torsion terms. Idempotent. Called by
      !! init_2Dh_axisym.
            real(dp) :: r, zc, Rc
            integer :: i, ntot1

            if (torsion_coeffs_defined) return

            if_torsion_zero = t2Dh%is_torsion_zero()
            if (if_torsion_zero) then
               torsion_coeffs_defined = .true.
               return
            end if

            ntot1 = lx1*ly1*lz1*nelv
            lambda_built = t2Dh%get_lambda()
            zc = t2Dh%get_axial_centre()
            Rc = t2Dh%get_curv_radius()

            call copy(zhat_f, xm1, ntot1)
            call copy(rhat_f, ym1, ntot1)
            do i = 1, ntot1
               r = rhat_f(i)
               zhat_f  (i) = zhat_f(i) - zc
               rhat_f  (i) = r - Rc
               invR_f  (i) = 1.0_dp/r
               lamR_f  (i) = lambda_built/r
               lamR2_f (i) = lambda_built/r**2
               lam2R2_f(i) = lambda_built**2/r**2
               lamzR3_f(i) = lambda_built*zhat_f(i)/r**3
            end do

      ! nonlinear solver: pointwise terms reading fZ_nl, fR_nl, fP_nl (ipert = 0)
            call add_neklab_forcing_term(torsion_forcing_pw, 'torsion')
            call add_neklab_source_term (torsion_source_pw,  'torsion')

            torsion_coeffs_defined = .true.
         end subroutine build_torsion_coeffs

         subroutine check_torsion_cache(caller)
      !! Stop if lambda has changed since build_torsion_coeffs ran, i.e. if
      !! t2Dh%init_geom was called after init_2Dh_axisym.
            character(len=*), intent(in) :: caller
            real(dp) :: lam_now
            if (.not. torsion_coeffs_defined) return
            lam_now = merge(0.0_dp, t2Dh%get_lambda(), t2Dh%is_torsion_zero())
            if (lam_now /= lambda_built) call nek_stop_error(
     &         'lambda changed after the torsion coefficients were built: '//
     &         'call t2Dh%init_geom before init_2Dh_axisym.', this_module, caller)
         end subroutine check_torsion_cache

      !====================================================================
      !     THE dtheta KERNEL AND ITS EXACT TRANSPOSE
      !====================================================================

         subroutine apply_dtheta(dout, din)
      !! dout = dtheta(din) = zhat*d(din)/dR - rhat*d(din)/dZ (Nek x = Z, y = R).
      !! Element-local collocation derivative, exactly as gradm1 returns it, so
      !! that apply_dtheta_t is its exact discrete transpose.
            real(dp), dimension(lv), intent(out) :: dout
            real(dp), dimension(lv), intent(in)  :: din
            integer :: i, ntot1
            ntot1 = lx1*ly1*lz1*nelv
            call gradm1(g_x, g_y, g_z, din)
            do i = 1, ntot1
               dout(i) = zhat_f(i)*g_y(i) - rhat_f(i)*g_x(i)
            end do
         end subroutine apply_dtheta

         subroutine apply_dtheta_t(dout, din)
      !! dout = dtheta^T(din), the exact discrete transpose of apply_dtheta under
      !! the Euclidean inner product: gradm1 run backwards. gradm1 pairs
      !! local_grad2(ur,us,u,N,e,dxm1,dytm1); its transpose is
      !! local_grad2_t(u,ur,us,N,e,dym1,dxtm1).
            real(dp), dimension(lv), intent(out) :: dout
            real(dp), dimension(lv), intent(in)  :: din
            real(dp), dimension(lxyz) :: ur, us, wtmp, ax, ay
            integer :: e, i, i0, ix, iy, N

            if (if3d) call nek_stop_error('apply_dtheta_t is 2D only.',
     &         this_module, 'apply_dtheta_t')

            N = lx1 - 1
            do e = 1, nelv
               i0 = (e - 1)*lxyz
               do i = 1, lxyz
                  ax(i) = -rhat_f(i0 + i)*din(i0 + i)
                  ay(i) =  zhat_f(i0 + i)*din(i0 + i)
               end do
               do iy = 1, ly1
               do ix = 1, lx1
                  i = ix + (iy - 1)*lx1
                  ur(i) = jacmi(i, e)*(rxm1(ix,iy,1,e)*ax(i) + rym1(ix,iy,1,e)*ay(i))
                  us(i) = jacmi(i, e)*(sxm1(ix,iy,1,e)*ax(i) + sym1(ix,iy,1,e)*ay(i))
               end do
               end do
               if (ifaxis) call setaxdy(ifrzer(e))
               call local_grad2_t(dout(i0 + 1), ur, us, N, 1, dym1, dxtm1, wtmp)
            end do
         end subroutine apply_dtheta_t

      !====================================================================
      !     THE D_lambda PAIR  (continuity / u_phi pressure)
      !====================================================================

         subroutine dphi_apply(dout, uin)
      !! dout = -(lambda/R)*dtheta(uin) on M1, unweighted. D_lam u = bm2*I12(dout).
            real(dp), dimension(lv), intent(out) :: dout
            real(dp), dimension(lv), intent(in)  :: uin
            integer :: ntot1
            ntot1 = lx1*ly1*lz1*nelv
            if (if_torsion_zero) then
               call rzero(dout, ntot1)
               return
            end if
            call apply_dtheta(dout, uin)
            call col2  (dout, lamR_f, ntot1)
            call chsign(dout, ntot1)
         end subroutine dphi_apply

         subroutine dphi_apply_t(dout, win)
      !! dout = [-(lambda/R)*dtheta]^T(win) = -dtheta^T[(lambda/R)*win], the exact
      !! discrete transpose of dphi_apply.
            real(dp), dimension(lv), intent(out) :: dout
            real(dp), dimension(lv), intent(in)  :: win
            integer :: ntot1
            ntot1 = lx1*ly1*lz1*nelv
            if (if_torsion_zero) then
               call rzero(dout, ntot1)
               return
            end if
            call col3          (p_t, win, lamR_f, ntot1)
            call apply_dtheta_t(dout, p_t)
            call chsign        (dout, ntot1)
         end subroutine dphi_apply_t

         subroutine dlam_grad_t(dout, p)
      !! dout = D_lam^T p = dphi_apply_t( I12^T (bm2 p) ): the weak, element-local
      !! u_phi pressure term (Nek's sign convention: the momentum RHS gets
      !! +D^T p). p lives on M2.
            real(dp), dimension(lv), intent(out) :: dout
            real(dp), dimension(lx2*ly2*lz2*lelv), intent(in) :: p
            call bm2_map12T_2Dh(m_t, p)
            call dphi_apply_t  (dout, m_t)
         end subroutine dlam_grad_t

         subroutine dphi_alpha_apply(dout, uin, alphaR)
      !! The alpha half of the streamwise divergence, without the slot swap:
      !! dout = (alpha/R) * uin. Used by the coupled pressure path (Step 7).
            real(dp), dimension(lv), intent(out) :: dout
            real(dp), dimension(lv), intent(in)  :: uin
            real(dp), dimension(lx1,ly1,lz1,lelv), intent(in) :: alphaR
            call col3(dout, uin, alphaR, lx1*ly1*lz1*nelv)
         end subroutine dphi_alpha_apply

         subroutine dphi_alpha_apply_t(dout, win, alphaR)
      !! Transpose of dphi_alpha_apply (a diagonal multiply).
            real(dp), dimension(lv), intent(out) :: dout
            real(dp), dimension(lv), intent(in)  :: win
            real(dp), dimension(lx1,ly1,lz1,lelv), intent(in) :: alphaR
            call col3(dout, win, alphaR, lx1*ly1*lz1*nelv)
         end subroutine dphi_alpha_apply_t

      !====================================================================
      !     KERNELS  (strong form, unweighted; shared by both solvers)
      !====================================================================

         subroutine adv_RZ(aZ, aR, cf, uZ, uR)
      !! a += cf * S(u), Z and R components (RHS sign):
      !!    Z: cf*(dtheta u_Z + u_R) ,  R: cf*(dtheta u_R - u_Z)
            real(dp), dimension(lv), intent(inout) :: aZ, aR
            real(dp), dimension(lv), intent(in)    :: cf, uZ, uR
            integer :: ntot1
            ntot1 = lx1*ly1*lz1*nelv
            call apply_dtheta(k_dZ, uZ)
            call apply_dtheta(k_dR, uR)
            call add2   (k_dZ, uR, ntot1)
            call sub2   (k_dR, uZ, ntot1)
            call addcol3(aZ, cf, k_dZ, ntot1)
            call addcol3(aR, cf, k_dR, ntot1)
         end subroutine adv_RZ

         subroutine adv_P(aP, cf, uP)
      !! a += cf * S(u), phi component: cf*dtheta u_phi.
            real(dp), dimension(lv), intent(inout) :: aP
            real(dp), dimension(lv), intent(in)    :: cf, uP
            integer :: ntot1
            ntot1 = lx1*ly1*lz1*nelv
            call apply_dtheta(k_dP, uP)
            call addcol3(aP, cf, k_dP, ntot1)
         end subroutine adv_P

         subroutine visc_RZ(aZ, aR, uZ, uR, uP)
      !! a += torsion part of the vector Laplacian (without nu), Z and R:
      !!    Z: lambda^2 Lth u_Z + (2 lambda^2/R^2) dth u_R + (lambda/R^2) u_phi
      !!       - (lambda^2 zhat/R^3) u_R - (lambda^2/R^2) u_Z
      !!    R: lambda^2 Lth u_R - (2 lambda^2/R^2) dth u_Z + (2 lambda/R^2) dth u_phi
      !!       + (lambda^2 zhat/R^3) u_Z - (lambda zhat/R^3) u_phi - (lambda^2/R^2) u_R
            real(dp), dimension(lv), intent(inout) :: aZ, aR
            real(dp), dimension(lv), intent(in)    :: uZ, uR, uP
            real(dp) :: lam, lam2
            integer :: ntot1
            ntot1 = lx1*ly1*lz1*nelv
            lam  = lambda_built
            lam2 = lam*lam
            call apply_dtheta(k_dZ, uZ)
            call apply_dtheta(k_dR, uR)
            call apply_dtheta(k_dP, uP)
      ! Z
            call apply_Ltheta(k_L, uZ)
            call add2s2 (aZ, k_L, lam2, ntot1)
            call add_cxy(aZ, lam2R2_f, k_dR,  2.0_dp, ntot1)
            call add_cxy(aZ, lamR2_f,  uP,    1.0_dp, ntot1)
            call add_cxy(aZ, lamzR3_f, uR,   -lam,    ntot1)
            call add_cxy(aZ, lam2R2_f, uZ,   -1.0_dp, ntot1)
      ! R
            call apply_Ltheta(k_L, uR)
            call add2s2 (aR, k_L, lam2, ntot1)
            call add_cxy(aR, lam2R2_f, k_dZ, -2.0_dp, ntot1)
            call add_cxy(aR, lamR2_f,  k_dP,  2.0_dp, ntot1)
            call add_cxy(aR, lamzR3_f, uZ,    lam,    ntot1)
            call add_cxy(aR, lamzR3_f, uP,   -1.0_dp, ntot1)
            call add_cxy(aR, lam2R2_f, uR,   -1.0_dp, ntot1)
         end subroutine visc_RZ

         subroutine visc_P(aP, uZ, uR, uP)
      !! a += torsion part of the vector Laplacian (without nu), phi:
      !!    lambda^2 Lth u_phi - (2 lambda/R^2) dth u_R + (lambda/R^2) u_Z
      !!    + (lambda zhat/R^3) u_R
            real(dp), dimension(lv), intent(inout) :: aP
            real(dp), dimension(lv), intent(in)    :: uZ, uR, uP
            integer :: ntot1
            ntot1 = lx1*ly1*lz1*nelv
            call apply_dtheta(k_dR, uR)
            call apply_Ltheta(k_L, uP)
            call add2s2 (aP, k_L, lambda_built**2, ntot1)
            call add_cxy(aP, lamR2_f,  k_dR, -2.0_dp, ntot1)
            call add_cxy(aP, lamR2_f,  uZ,    1.0_dp, ntot1)
            call add_cxy(aP, lamzR3_f, uR,    1.0_dp, ntot1)
         end subroutine visc_P

      !====================================================================
      !     LINEAR SOLVER  (alpha = 0 terms; called from makefp_2Dh_axisym)
      !====================================================================

         subroutine add_torsion_forcing_RZ(rhs_z, rhs_r, jp_)
      !! Z and R momentum RHS of slot jp_, injected before EXT3:
      !!    rho*bm1*[ c(U) S(u) + c(u) S(U) ]  +  vdiff*bm1*[viscous(u)],
      !! c(v) = lambda v_phi/R. Exact linearisation of refresh_torsion_forcing_nl.
            real(dp), dimension(lv), intent(inout) :: rhs_z, rhs_r
            integer, intent(in) :: jp_
            integer :: ntot1
            if (if_torsion_zero) return
            ntot1 = lx1*ly1*lz1*nelv
      ! ---- advection + frame rotation (rho)
            call rzero(acc_a, ntot1)
            call rzero(acc_b, ntot1)
            call col3  (cf_1, lamR_f, t(1,1,1,1,1), ntot1)       ! c(U)
            call adv_RZ(acc_a, acc_b, cf_1, vxp(1,jp_), vyp(1,jp_))
            call col3  (cf_1, lamR_f, tp(1,1,jp_), ntot1)        ! c(u)
            call adv_RZ(acc_a, acc_b, cf_1, vx, vy)
            call weight_forcing(acc_a, vtrans(1,1,1,1,1))
            call weight_forcing(acc_b, vtrans(1,1,1,1,1))
            call add2(rhs_z, acc_a, ntot1)
            call add2(rhs_r, acc_b, ntot1)
      ! ---- viscous (nu), including -lambda^2/R^2 on Z and R
            call rzero  (acc_a, ntot1)
            call rzero  (acc_b, ntot1)
            call visc_RZ(acc_a, acc_b, vxp(1,jp_), vyp(1,jp_), tp(1,1,jp_))
            call weight_forcing(acc_a, vdiff(1,1,1,1,1))
            call weight_forcing(acc_b, vdiff(1,1,1,1,1))
            call add2(rhs_z, acc_a, ntot1)
            call add2(rhs_r, acc_b, ntot1)
         end subroutine add_torsion_forcing_RZ

         subroutine add_torsion_forcing_phi(rhs_phi, jp_)
      !! u_phi RHS of slot jp_, injected before EXT3:
      !!    rhoCp*bm1*[ c(U) S(u) + c(u) S(U) ]_phi + vdiff*bm1*[viscous(u)]_phi
      !!    + D_lam^T p'   (p' = prp(:,jp_), the current perturbation pressure)
            real(dp), dimension(lv), intent(inout) :: rhs_phi
            integer, intent(in) :: jp_
            integer :: ntot1
            if (if_torsion_zero) return
            ntot1 = lx1*ly1*lz1*nelv
      ! ---- advection (rhoCp)
            call rzero(acc_c, ntot1)
            call col3 (cf_1, lamR_f, t(1,1,1,1,1), ntot1)        ! c(U)
            call adv_P(acc_c, cf_1, tp(1,1,jp_))
            call col3 (cf_1, lamR_f, tp(1,1,jp_), ntot1)         ! c(u)
            call adv_P(acc_c, cf_1, t(1,1,1,1,1))
            call weight_forcing(acc_c, vtrans(1,1,1,1,2))
            call add2(rhs_phi, acc_c, ntot1)
      ! ---- viscous (nu)
            call rzero (acc_c, ntot1)
            call visc_P(acc_c, vxp(1,jp_), vyp(1,jp_), tp(1,1,jp_))
            call weight_forcing(acc_c, vdiff(1,1,1,1,1))
            call add2(rhs_phi, acc_c, ntot1)
      ! ---- pressure (weak, already weighted)
            call dlam_grad_t(acc_c, prp(1,jp_))
            call add2(rhs_phi, acc_c, ntot1)
         end subroutine add_torsion_forcing_phi

         subroutine add_torsion_coupling(rhs_z, rhs_r, rhs_phi, alpha, jp_)
      !! alpha*lambda viscous terms that couple the Re/Im slots (alpha /= 0 only,
      !! currently refused by prepare_step_2Dh_axisym; kept for Step 7). With
      !! D = -i alpha - lambda dtheta, the single-i terms of the vector Laplacian:
      !!   T1: +(2 alpha lambda/R^2) dtheta(u_i)
      !!   T2: -(alpha lambda zhat/R^3) u_i
      !!   T3: R-eq -(2 alpha lambda/R^2) u_Z, Z-eq +(2 alpha lambda/R^2) u_R
      !! all times +i: read slot jp_, write ipert with Re[i w] = -w^I, Im[i w] = w^R.
            real(dp), dimension(lv,lpert), intent(inout) :: rhs_z, rhs_r, rhs_phi
            real(dp), intent(in) :: alpha
            integer, intent(in) :: jp_
            real(dp) :: lam, sgn
            integer :: ipert, ntot1

            if (if_torsion_zero) return
            if (alpha == 0.0_dp) return
            ntot1 = lx1*ly1*lz1*nelv
            lam = lambda_built
            ipert = npert + 1 - jp_
            sgn = merge(-1.0_dp, 1.0_dp, ipert == 1)

            call rzero(acc_a, ntot1)
            call rzero(acc_b, ntot1)
            call rzero(acc_c, ntot1)
            call apply_dtheta(k_dZ, vxp(1,jp_))
            call apply_dtheta(k_dR, vyp(1,jp_))
            call apply_dtheta(k_dP, tp(1,1,jp_))
      ! T1
            call add_cxy(acc_a, lamR2_f, k_dZ, 2.0_dp*alpha*lam, ntot1)
            call add_cxy(acc_b, lamR2_f, k_dR, 2.0_dp*alpha*lam, ntot1)
            call add_cxy(acc_c, lamR2_f, k_dP, 2.0_dp*alpha*lam, ntot1)
      ! T2
            call add_cxy(acc_a, lamzR3_f, vxp(1,jp_),  -alpha, ntot1)
            call add_cxy(acc_b, lamzR3_f, vyp(1,jp_),  -alpha, ntot1)
            call add_cxy(acc_c, lamzR3_f, tp(1,1,jp_), -alpha, ntot1)
      ! T3
            call add_cxy(acc_b, lamR2_f, vxp(1,jp_), -2.0_dp*alpha*lam, ntot1)
            call add_cxy(acc_a, lamR2_f, vyp(1,jp_),  2.0_dp*alpha*lam, ntot1)

            call weight_forcing(acc_a, vdiff(1,1,1,1,1))
            call weight_forcing(acc_b, vdiff(1,1,1,1,1))
            call weight_forcing(acc_c, vdiff(1,1,1,1,1))
            call add2s2(rhs_z(1,ipert),   acc_a, sgn, ntot1)
            call add2s2(rhs_r(1,ipert),   acc_b, sgn, ntot1)
            call add2s2(rhs_phi(1,ipert), acc_c, sgn, ntot1)
         end subroutine add_torsion_coupling

      !====================================================================
      !     NONLINEAR SOLVER  (called from nek_advance_2Dh_axisym_nl)
      !====================================================================

         subroutine refresh_torsion_forcing_nl()
      !! Torsion forcing of the nonlinear step, as ACCELERATIONS (Nek multiplies
      !! userf by rho*bm1; neklab_source multiplies by rhoCp, Nek by bm1), from
      !! the current U^n and P^n. Same kernels as the linear solver:
      !!    Z, R : c(U) S(U) + nu*viscous(U)
      !!    phi  : c(U) S(U) + nu*viscous(U) + D_lam^T P/(bm1*rhoCp)
      !! with nu = vdiff/vtrans, so that Nek rebuilds exactly the weights of the
      !! linear solver.
            integer :: ntot1
            if (if_torsion_zero) return
            ntot1 = lx1*ly1*lz1*nelv
            call rzero(fZ_nl, ntot1)
            call rzero(fR_nl, ntot1)
            call rzero(fP_nl, ntot1)
      ! ---- advection + frame rotation
            call col3  (cf_1, lamR_f, t(1,1,1,1,1), ntot1)       ! c(U)
            call adv_RZ(fZ_nl, fR_nl, cf_1, vx, vy)
            call adv_P (fP_nl, cf_1, t(1,1,1,1,1))
      ! ---- viscous
            call rzero  (acc_a, ntot1)
            call rzero  (acc_b, ntot1)
            call rzero  (acc_c, ntot1)
            call visc_RZ(acc_a, acc_b, vx, vy, t(1,1,1,1,1))
            call visc_P (acc_c, vx, vy, t(1,1,1,1,1))
            call col2   (acc_a, vdiff(1,1,1,1,1),  ntot1)        ! nu = vdiff/rho
            call invcol2(acc_a, vtrans(1,1,1,1,1), ntot1)
            call add2   (fZ_nl, acc_a, ntot1)
            call col2   (acc_b, vdiff(1,1,1,1,1),  ntot1)
            call invcol2(acc_b, vtrans(1,1,1,1,1), ntot1)
            call add2   (fR_nl, acc_b, ntot1)
            call col2   (acc_c, vdiff(1,1,1,1,1),  ntot1)        ! vdiff/rhoCp
            call invcol2(acc_c, vtrans(1,1,1,1,2), ntot1)
            call add2   (fP_nl, acc_c, ntot1)
      ! ---- pressure (weak -> divide by bm1 node by node, Nek multiplies back)
            call dlam_grad_t(cf_1, pr)
            call invcol2(cf_1, bm1, ntot1)
            call invcol2(cf_1, vtrans(1,1,1,1,2), ntot1)
            call add2   (fP_nl, cf_1, ntot1)
         end subroutine refresh_torsion_forcing_nl

         subroutine set_torsion_usrdiv()
      !! Continuity of the nonlinear step. Nek's incomprn enforces
      !! D u^{n+1} = (dt/bd(1)) bm2 usrdiv; requiring D u + D_lam U_phi^{n+1} = 0
      !! with D_lam U = bm2 I12(dphi_apply U) gives
      !!    usrdiv = -(bd(1)/dt) I12( dphi_apply(U_phi^{n+1}) ).
      !! Call after heat(2) (U_phi^{n+1} in t) and before fluid(2).
            character(len=*), parameter :: this_procedure = 'set_torsion_usrdiv'
            character(len=128) :: msg
            real(dp), external :: glsc2
            real(dp) :: dmean
            integer :: e, ntot2
            if (if_torsion_zero) return
            ntot2 = lx2*ly2*lz2*nelv
            call dphi_apply(acc_a, t(1,1,1,1,1))
            do e = 1, nelv
               call map12_2Dh(usrdiv(1,1,1,e), acc_a((e - 1)*lxyz + 1), e)
            end do
            call cmult(usrdiv, -bd(1)/dt, ntot2)
      ! ortho removes the mean of the pressure RHS; with U_phi = 0 on the wall the
      ! continuous mean is zero, so this must sit at discretisation level.
            dmean = glsc2(usrdiv, bm2, ntot2)/volvm2
            write (msg, '(A,1pE12.4)') 'mean(usrdiv) = ', dmean
            call nek_log_debug(msg, this_module, this_procedure)
         end subroutine set_torsion_usrdiv

         subroutine torsion_forcing_pw(fx, fy, fz, ix, iy, iz, e)
      !! Registered pointwise forcing (ipert = 0): Z and R torsion accelerations.
            real(dp), intent(inout) :: fx, fy, fz
            integer, intent(in) :: ix, iy, iz, e
            integer :: i
            i = ix + lx1*((iy - 1) + ly1*((iz - 1) + lz1*(e - 1)))
            fx = fx + fZ_nl(i)
            fy = fy + fR_nl(i)
         end subroutine torsion_forcing_pw

         subroutine torsion_source_pw(q, ix, iy, iz, e)
      !! Registered pointwise source (ipert = 0): u_phi torsion acceleration.
            real(dp), intent(inout) :: q
            integer, intent(in) :: ix, iy, iz, e
            integer :: i
            i = ix + lx1*((iy - 1) + ly1*((iz - 1) + lz1*(e - 1)))
            q = q + fP_nl(i)
         end subroutine torsion_source_pw

      !====================================================================
      !     HELPERS
      !====================================================================

         subroutine apply_Ltheta(dout, din)
      !! Ltheta(u) = (1/R) dtheta[ (1/R) dtheta u ], with the intermediate field
      !! mass-averaged across elements before the second derivative.
            real(dp), dimension(lv), intent(out) :: dout
            real(dp), dimension(lv), intent(in)  :: din
            integer :: ntot1
            ntot1 = lx1*ly1*lz1*nelv
            call apply_dtheta(l_g, din)
            call col2 (l_g, invR_f, ntot1)
            call col2 (l_g, bm1, ntot1)
            call dssum(l_g, lx1, ly1, lz1)
            call col2 (l_g, binvm1, ntot1)
            call apply_dtheta(dout, l_g)
            call col2 (dout, invR_f, ntot1)
         end subroutine apply_Ltheta

         subroutine add_cxy(dst, a, b, s, ntot1)
      !! dst += s*a*b, elementwise.
            real(dp), dimension(lv), intent(inout) :: dst
            real(dp), dimension(lv), intent(in)    :: a, b
            real(dp), intent(in) :: s
            integer, intent(in)  :: ntot1
            call col3  (c_t, a, b, ntot1)
            call add2s2(dst, c_t, s, ntot1)
         end subroutine add_cxy

         subroutine weight_forcing(rhs, prop)
      !! Strong form -> the weighting bfxp/bfyp/bqp expect: times the property
      !! (rho for advection, vdiff for viscous terms) and bm1.
            real(dp), dimension(lv), intent(inout) :: rhs
            real(dp), dimension(lv), intent(in)    :: prop
            integer :: ntot1
            ntot1 = lx1*ly1*lz1*nelv
            call col2(rhs, prop, ntot1)
            call col2(rhs, bm1, ntot1)
         end subroutine weight_forcing

      end module neklab_2Dh_torsion
