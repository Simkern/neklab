      module neklab_2Dh_torsion
      !! Torsion terms of the 2Dh axisymmetric solver: helical pipe in the
      !! Germano (Frenet) frame, from helix_2Dh_derivation.tex. Nothing here
      !! depends on either stepper. The steppers (neklab_2Dh, neklab_2Dh_axisym)
      !! call:
      !!
      !!   init_torsion_2Dh           once, from init_2Dh_axisym
      !!   torsion_refresh_nl         nonlinear step, from U^n, before heat(1)
      !!   torsion_add_pressure_nl    nonlinear step, after heat(1): bq += D_lam^T p_ext
      !!   torsion_add_linear_uv/_s   linear step, makefp: bfxp/bfyp/bqp += L(U; u)
      !!                              (adjoint mode: += L^T(U; v), see ADJOINT below)
      !!   torsion_prepare_adjoint    adjoint linear step, from prepare_step_2Dh_axisym
      !!   torsion_transport_velocity CFL hook (rotation advection c(U) dth)
      !!   add_dlam_div, dlam_grad    the continuity/pressure pair, both steppers
      !!   add_lambda_block_E         the lambda block of the pressure operator
      !!   lambda_correction          the u_s pressure correction
      !!
      !!--------------------------------------------------------------------
      !! NOTATION  (helix_2Dh_derivation.tex, Sec. 2-3)
      !!--------------------------------------------------------------------
      !!
      !!   Z = xm1, R = ym1, zeta = Z - z_c, rho = R - R_c, lambda = tau/kappa
      !!   dth f = zeta df/dR - rho df/dZ        (rotation about (z_c, R_c))
      !!   D_s   = (1/R)(P - lambda dth),  P = -i alpha (P = 0 for the base flow)
      !!   u = (u_Z, u_R, u_s) = (vx, vy, t(:,1))  base;  (vxp, vyp, tp) perturbation
      !!
      !! Terms Nek does not provide, as right-hand-side force densities (the G
      !! numbers are those of the derivation):
      !!
      !!   G2  add_torsion_advection     rho c(a) S(b),   c(a) = lambda a_s / R
      !!                                 S(b) = (dth b_Z + b_R, dth b_R - b_Z, dth b_s)
      !!       nonlinear: G2(U,U);   linear: G2(U,u) + G2(u,U)   (exact by construction)
      !!   G3  add_torsion_viscous       mu * [torsion part of (D_s + M)^2 v] at P = 0
      !!                                 (lambda^2/R^2 shifts on Z, R included)
      !!   G6  add_alpha_torsion_viscous the P*lambda part of mu (D_s + M)^2 u, linear
      !!       step at alpha /= 0:  mu i alpha lambda [(2/R^2) dth u - (zeta/R^3) u
      !!       + (2/R^2)(u_R, -u_Z, 0)]; odd in P, so it swaps slots (+i rule)
      !!   G4  pressure, u_s eq.         +(lambda/R) dth p   =  weak  D_lam^T p
      !!   G5  continuity                -(lambda/R) dth u_s =  weak  D_lam u_s
      !!
      !! G2 and G3 are computed ONCE, as weak (bm1-weighted) force vectors, by the
      !! same routines for both solvers. The linear solver adds them to
      !! bfxp/bfyp/bqp; the nonlinear one passes F/(rho bm1) through userf/userq,
      !! which Nek multiplies back by bm1 (makeuf/makeuq) and rho (makeabf /
      !! neklab_source). Both go through EXT3 at the same time level.
      !!
      !! G4/G5 are IMPLICIT in both solvers: D_lam is part of the projection
      !! (E = D B^-1 D^T with D = [D_in, D_lam]) and the correction updates u_s.
      !! An explicit D_lam^T P has a dt-independent pressure-error instability
      !! for (lambda delta/(1-delta))^2 > 1/7 with EXT3.
      !!
      !! The D_lam pair is made exactly mean-free (1^T D_lam u = 0, D_lam^T 1 = 0),
      !! so that E 1 = 0 holds exactly: ortho and the hsmg preconditioner assume it.
      !!
      !!--------------------------------------------------------------------
      !! ADJOINT (ifadj, linear step only)
      !!--------------------------------------------------------------------
      !!
      !! The adjoint step runs at the SAME alpha and reuses the step unchanged
      !! except for the explicit force: with F(u) the weak torsion force of the
      !! direct step (local, linear in u), the adjoint force is its exact
      !! algebraic transpose,  v^T F(u) = F^T(v)^T u  to round-off. For M = rho B
      !! the adjoint of L = M^-1 F in the M-inner product is M^-1 F^T, so F^T is
      !! the right weak force while rho = rho cp is one constant (2Dh requires it).
      !!
      !!   G2 lin.  add_torsion_advection_lin_adj
      !!   G3       add_torsion_viscous_adj (the lambda^2 L_theta block is symmetric)
      !!   G6       slot block [[0, -a K], [a K, 0]] -> [[0, a K^T], [-a K^T, 0]]:
      !!            slot 1 += +alpha K^T(v_Im), slot 2 += -alpha K^T(v_Re)
      !!   D_lam    part of the projection P = I - B^-1 D^T E^-1 D: self-adjoint
      !!
      !! add_adjoint_transport_div adds the part of the adjoint of Nek's transport
      !! -(U.grad) u that Nek's adjoint (+(U.grad) v) leaves out because it
      !! assumes div U = 0. With torsion the in-plane divergence of the base is
      !!   div_in U = dU_Z/dZ + (1/R) d(R U_R)/dR = (lambda/R) dth U_s   (G5),
      !! so the adjoint transport is +(U.grad) v + (div_in U) v. div_in U is taken
      !! from the base flow by gradm1 at every adjoint step (torsion_prepare_adjoint).
      !!
      !!--------------------------------------------------------------------
      !! DISCRETE dth
      !!--------------------------------------------------------------------
      !!
      !! Element-local collocation derivative, the r,s derivatives of gradm1:
      !!   dth f = A_r f_r + A_s f_s ,
      !!   A_r = (zeta ry - rho rx)/J ,  A_s = (zeta sy - rho sx)/J
      !! with Nek's local_grad2 (dxm1, dytm1); dtheta_t is its exact transpose
      !! with local_grad2_t, never "-dth". The meshes must not touch the axis
      !! (R > 0 is checked by init_2Dh_axisym), so no element uses setaxdy.
         use LightKrylov, only: dp
         use stdlib_optval, only: optval
         use neklab_nek_setup, only: nek_log_message, nek_log_warning,
     &                               nek_stop_error, set_cfl_velocity_hook
         use neklab_nek_forcing, only: add_neklab_forcing_term, add_neklab_source_term
         use neklab_2Dh_common, only: wmask, bm2_map12_2Dh, bm2_map12T_2Dh
         use neklab_t2Dh, only: t2Dh
         implicit none
         include "SIZE"
         include "TOTAL"
         include "ADJOINT"

         private
         character(len=*), parameter, private :: this_module = 'neklab_2Dh_torsion'
         integer, parameter :: lv = lx1*ly1*lz1*lelv
         integer, parameter :: lp = lx2*ly2*lz2*lelv
         integer, parameter :: lxy = lx1*ly1

      ! --- state
         logical,  public, protected :: lambda_eq_zero = .true.   ! no torsion (torus), or init not done
         logical,  private :: torsion_built = .false.
         real(dp), private :: lambda_t  = 0.0_dp        ! lambda = tau/kappa
         real(dp), private :: alpha_t   = 0.0_dp        ! alpha of the linear step (G6)
         real(dp), private :: z_centre  = 0.0_dp
         real(dp), private :: R_centre  = 0.0_dp
         real(dp), private :: vol2      = 0.0_dp        ! sum(bm2), mean removal

      ! --- coefficient fields (M1), built once by init_torsion_2Dh
         real(dp), dimension(lv), private :: zeta_f, rho_f, invR_f
         real(dp), dimension(lv), private :: iR2_f      ! 1/R^2
         real(dp), dimension(lv), private :: zR3_f      ! zeta/R^3
         real(dp), dimension(lv), private :: mlamR_f    ! -lambda/R   (D_lam)
         real(dp), dimension(lv), private :: Ar_f, As_f ! dth = Ar d/dr + As d/ds
         real(dp), dimension(lv), private :: divU_f = 0.0_dp  ! div_in U of the base (adjoint)

      ! --- nonlinear forcing as accelerations, read pointwise by the registered terms
         real(dp), dimension(lv), private :: fZ_nl = 0.0_dp, fR_nl = 0.0_dp, fs_nl = 0.0_dp

      ! --- linear forcing of the current slot, computed at the velocity stage of
      !     makefp and consumed (u_s part) at its scalar stage
         real(dp), dimension(lv), private :: lin_FZ, lin_FR, lin_Fs

      ! --- scratch, one group per call level
         real(dp), dimension(lv), private :: k_dZ, k_dR, k_ds, k_t1, k_t2   ! kernels
         real(dp), dimension(lv), private :: k_aZ, k_aR, k_as               ! kernels (2nd arg)
         real(dp), dimension(lv), private :: k_gZ, k_gR, k_gs, k_h          ! adjoint kernels
         real(dp), dimension(lv), private :: p_m1a, p_m1b                   ! pressure pair
         real(dp), dimension(lp), private :: p_m2a                          ! pressure pair
         real(dp), dimension(lv), private :: e_m1                           ! E block, correction
         real(dp), dimension(lp), private :: n_pext                         ! nonlinear p_ext
         real(dp), dimension(lxy), private :: e_r, e_s, e_w                 ! dth element work

         public :: init_torsion_2Dh
         public :: torsion_lambda
         public :: set_torsion_alpha
         public :: check_wall_centering
         public :: dtheta, dtheta_t
         public :: add_dlam_div, dlam_grad
         public :: add_lambda_block_E, lambda_correction
         public :: torsion_force_nl, torsion_force_lin, torsion_force_lin_adj
         public :: torsion_refresh_nl, torsion_add_pressure_nl
         public :: torsion_add_linear_uv, torsion_add_linear_s
         public :: torsion_prepare_adjoint
         public :: torsion_transport_velocity
      ! kernels (also called by the tests)
         public :: add_torsion_advection, add_torsion_viscous, add_alpha_torsion_viscous
         public :: add_torsion_advection_lin_adj, add_torsion_viscous_adj
         public :: add_alpha_torsion_viscous_adj, add_adjoint_transport_div

      contains

      !====================================================================
      !     SETUP
      !====================================================================

         subroutine set_torsion_alpha(alpha)
      !! Wavenumber of the linear step, for the alpha*lambda viscous coupling G6.
      !! Set by build_axisym_coeffs (neklab_2Dh_axisym).
            real(dp), intent(in) :: alpha
            alpha_t = alpha
         end subroutine set_torsion_alpha

         real(dp) function torsion_lambda()
            torsion_lambda = merge(0.0_dp, lambda_t, lambda_eq_zero)
         end function torsion_lambda

         subroutine init_torsion_2Dh()
      !! Build the coefficients from t2Dh%init_geom and register the nonlinear
      !! terms and the CFL velocity. Idempotent. Called by init_2Dh_axisym, after
      !! the masks are built. A torus (lambda = 0, or no init_geom) leaves
      !! everything off: every entry point below then returns at once.
            character(len=*), parameter :: this_procedure = 'init_torsion_2Dh'
            character(len=256) :: msg
            real(dp), external :: glsum
            integer :: n
            real(dp) :: dev

            if (torsion_built) return
            torsion_built = .true.
            lambda_eq_zero = .true.
            if (t2Dh%is_torsion_zero()) then
               call nek_log_message('lambda = 0: torsion terms inactive.', this_module, this_procedure)
               return
            end if
            if (ldim /= 2) call nek_stop_error('torsion requires a 2D mesh.', this_module, this_procedure)

            lambda_t = t2Dh%get_lambda()
            z_centre = t2Dh%get_axial_centre()
            R_centre = t2Dh%get_curv_radius()
            n = lx1*ly1*lz1*nelv

            call copy   (zeta_f, xm1, n)
            call cadd   (zeta_f, -z_centre, n)
            call copy   (rho_f, ym1, n)
            call cadd   (rho_f, -R_centre, n)
            call invers2(invR_f, ym1, n)
            call col3   (iR2_f, invR_f, invR_f, n)
            call col3   (zR3_f, zeta_f, iR2_f, n)
            call col2   (zR3_f, invR_f, n)
            call cmult2 (mlamR_f, invR_f, -lambda_t, n)
      ! A_r = (zeta ry - rho rx)/J,  A_s = (zeta sy - rho sx)/J
            call col3   (Ar_f, zeta_f, rym1, n)
            call subcol3(Ar_f, rho_f, rxm1, n)
            call col2   (Ar_f, jacmi, n)
            call col3   (As_f, zeta_f, sym1, n)
            call subcol3(As_f, rho_f, sxm1, n)
            call col2   (As_f, jacmi, n)
            vol2 = glsum(bm2, lx2*ly2*lz2*nelv)

            call check_full_section(this_procedure)
            call check_masks_equal(this_procedure)
            call check_wall_centering(dev, verbose=.true.)

            call add_neklab_forcing_term(torsion_forcing_pw, 'torsion')
            call add_neklab_source_term (torsion_source_pw,  'torsion')
            call set_cfl_velocity_hook(torsion_transport_velocity)

            lambda_eq_zero = .false.
            write (msg, '(A,ES14.6,A,ES14.6,A,ES14.6)') 'torsion active: lambda =', lambda_t,
     &         ', z_c =', z_centre, ', R_c =', R_centre
            call nek_log_message(msg, this_module, this_procedure)
         end subroutine init_torsion_2Dh

         subroutine check_full_section(caller)
      !! Torsion breaks Z -> -Z symmetry and dth is wall-normal on a plane
      !! Z = z_c: a helix needs the full cross-section (no SYM/ON faces).
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
            if (nbad > 0) call nek_stop_error('lambda /= 0 requires the full cross-section: '//
     &         'SYM/ON faces found.', this_module, caller)
         end subroutine check_full_section

         subroutine check_masks_equal(caller)
      !! The nonlinear step masks u_s with tmask (Nek's scalar), the linear one
      !! with wmask. With torsion the two steps must be the same operator.
            character(len=*), intent(in) :: caller
            integer :: i, nbad
            integer, external :: iglsum
            nbad = 0
            do i = 1, lx1*ly1*lz1*nelv
               if (wmask(i,1,1,1) /= tmask(i,1,1,1,1)) nbad = nbad + 1
            end do
            nbad = iglsum(nbad, 1)
            if (nbad > 0) call nek_stop_error('lambda /= 0 requires tmask = wmask '//
     &         '(temperature BCs must match the velocity BCs).', this_module, caller)
         end subroutine check_masks_equal

         subroutine check_wall_centering(maxdev, verbose)
      !! dth is tangent to the wall only if the wall is the circle of radius
      !! t2Dh%radius about (z_c, R_c). Reports the largest deviation of a
      !! velocity-Dirichlet boundary node from that circle. Midpoint ('m') curved
      !! sides give O(h^4) deviations; 'C' arcs give round-off.
            real(dp), intent(out) :: maxdev
            logical, optional, intent(in) :: verbose
            character(len=*), parameter :: this_procedure = 'check_wall_centering'
            character(len=256) :: msg
            real(dp) :: a, dist, zcl, Rcl
            integer :: ie, ifc, ix, iy, kx1, kx2, ky1, ky2, kz1, kz2, iz
            real(dp), external :: glmax
            a = t2Dh%get_radius()
            zcl = t2Dh%get_axial_centre()
            Rcl = t2Dh%get_curv_radius()
            maxdev = 0.0_dp
            do ie = 1, nelv
               do ifc = 1, 2*ndim
                  if (cbc(ifc, ie, 1) /= 'W  ' .and. cbc(ifc, ie, 1) /= 'v  ') cycle
                  call facind(kx1, kx2, ky1, ky2, kz1, kz2, lx1, ly1, lz1, ifc)
                  do iz = kz1, kz2
                  do iy = ky1, ky2
                  do ix = kx1, kx2
                     dist = hypot(xm1(ix,iy,iz,ie) - zcl, ym1(ix,iy,iz,ie) - Rcl)
                     maxdev = max(maxdev, abs(dist - a))
                  end do
                  end do
                  end do
               end do
            end do
            maxdev = glmax(maxdev, 1)
            write (msg, '(A,ES12.4)') 'max wall deviation from the circle about (z_c, R_c): ', maxdev
            if (maxdev > 1.0e-2_dp*a) then
               call nek_stop_error(trim(msg)//' > 1e-2*radius: the wall is not centred on '//
     &            '(z_c, R_c). Mesh/geometry mismatch (e.g. coordinates read by load_fld: '//
     &            'use load_fld(''file U P T'')).', this_module, this_procedure)
            end if
            if (optval(verbose, .false.)) then
               if (maxdev > 1.0e-3_dp*a) then
                  call nek_log_warning(trim(msg)//' -- (z_c, R_c) is not the wall centre?',
     &               this_module, this_procedure)
               else
                  call nek_log_message(msg, this_module, this_procedure)
               end if
            end if
         end subroutine check_wall_centering

      !====================================================================
      !     dth AND ITS EXACT TRANSPOSE
      !====================================================================

         subroutine dtheta(out, f)
      !! out = dth f = A_r f_r + A_s f_s, element-local (no dssum).
            real(dp), dimension(lv), intent(out) :: out
            real(dp), dimension(lv), intent(in)  :: f
            integer :: ie, i, i0
            do ie = 1, nelv
               call local_grad2(e_r, e_s, f, lx1-1, ie, dxm1, dytm1)
               i0 = lxy*(ie - 1)
               do i = 1, lxy
                  out(i0+i) = Ar_f(i0+i)*e_r(i) + As_f(i0+i)*e_s(i)
               end do
            end do
         end subroutine dtheta

         subroutine dtheta_t(out, g)
      !! out = dth^T g = D_r^T (A_r g) + D_s^T (A_s g), element-local.
            real(dp), dimension(lv), intent(out) :: out
            real(dp), dimension(lv), intent(in)  :: g
            integer :: ie, i, i0
            do ie = 1, nelv
               i0 = lxy*(ie - 1)
               do i = 1, lxy
                  e_r(i) = Ar_f(i0+i)*g(i0+i)
                  e_s(i) = As_f(i0+i)*g(i0+i)
               end do
               call local_grad2_t(out, e_r, e_s, lx1-1, ie, dym1, dxtm1, e_w)
            end do
         end subroutine dtheta_t

      !====================================================================
      !     THE CONTINUITY / PRESSURE PAIR (mean-free)
      !====================================================================
      !
      !   A u     = bm2 * I12( -(lambda/R) dth u )                 (M1 -> M2)
      !   D_lam   = P0 A,    P0 d = d - bm2 (1^T d)/sum(bm2)
      !   D_lam^T = A^T P0^T,  P0^T p = p - 1 (bm2^T p)/sum(bm2)
      !   A^T q   = dth^T( -(lambda/R) I12^T(bm2 q) )
      !
      ! so 1^T D_lam u = 0 and D_lam^T 1 = 0 exactly, and the pair is an exact
      ! transpose. The removed mean is O(discretisation): the continuous
      ! integral of dth u vanishes on a circular wall centred at (z_c, R_c).

         subroutine add_dlam_div(d, u)
      !! d += D_lam u   (d on M2, u on M1). No-op without torsion.
            real(dp), dimension(lp), intent(inout) :: d
            real(dp), dimension(lv), intent(in)    :: u
            real(dp), external :: glsum
            real(dp) :: s
            integer :: ntot2
            if (lambda_eq_zero) return
            ntot2 = lx2*ly2*lz2*nelv
            call dtheta(p_m1a, u)
            call col2  (p_m1a, mlamR_f, lx1*ly1*lz1*nelv)
            call bm2_map12_2Dh(p_m2a, p_m1a)
            s = glsum(p_m2a, ntot2)/vol2
            call add2s2(p_m2a, bm2, -s, ntot2)
            call add2  (d, p_m2a, ntot2)
         end subroutine add_dlam_div

         subroutine dlam_grad(g, p)
      !! g = D_lam^T p   (g on M1, element-local, unassembled; p on M2).
      !! Zero without torsion.
            real(dp), dimension(lv), intent(out) :: g
            real(dp), dimension(lp), intent(in)  :: p
            real(dp), external :: glsc2
            real(dp) :: pm
            integer :: ntot2
            ntot2 = lx2*ly2*lz2*nelv
            if (lambda_eq_zero) then
               call rzero(g, lx1*ly1*lz1*nelv)
               return
            end if
            pm = glsc2(bm2, p, ntot2)/vol2
            call copy(p_m2a, p, ntot2)
            call cadd(p_m2a, -pm, ntot2)
            call bm2_map12T_2Dh(p_m1b, p_m2a)                 ! I12^T (bm2 (p - pm))
            call col2(p_m1b, mlamR_f, lx1*ly1*lz1*nelv)
            call dtheta_t(g, p_m1b)
         end subroutine dlam_grad

         subroutine add_lambda_block_E(ap, p, h2inv, mask)
      !! ap += D_lam (mask h2inv B^-1) D_lam^T p : the lambda block of the
      !! pressure operator E = D (h2inv B^-1) D^T. Symmetric, positive
      !! semi-definite, annihilates constants.
            real(dp), dimension(lp), intent(inout) :: ap
            real(dp), dimension(lp), intent(in)    :: p
            real(dp), dimension(lv), intent(in)    :: h2inv, mask
            if (lambda_eq_zero) return
            call lambda_correction(e_m1, p, h2inv, mask)
            call add_dlam_div(ap, e_m1)
         end subroutine add_lambda_block_E

         subroutine lambda_correction(dw, dp_, h2inv, mask)
      !! dw = mask h2inv B^-1 D_lam^T dp : the u_s pressure correction, same
      !! sign and same B^-1 as the in-plane one (opgradt -> opbinv).
            real(dp), dimension(lv), intent(out) :: dw
            real(dp), dimension(lp), intent(in)  :: dp_
            real(dp), dimension(lv), intent(in)  :: h2inv, mask
            integer :: ntot1
            call dlam_grad(dw, dp_)
            if (lambda_eq_zero) return
            ntot1 = lx1*ly1*lz1*nelv
            call col2 (dw, mask, ntot1)
            call dssum(dw, lx1, ly1, lz1)
            call col2 (dw, binvm1, ntot1)
            call col2 (dw, h2inv, ntot1)
         end subroutine lambda_correction

      !====================================================================
      !     WEAK FORCE KERNELS (G2, G3, G6)
      !====================================================================
      ! One fused loop per kernel, written as the formula of the derivation;
      ! the R-dependent factors are the coefficient fields of init.

         subroutine add_torsion_advection(FZ, FR, Fs, aZ, aR, as, bZ, bR, bs)
      !! G2(a, b):  F += rho bm1 c(a) S(b),  c(a) = lambda a_s/R,
      !! S(b) = (dth b_Z + b_R, dth b_R - b_Z, dth b_s).
            real(dp), dimension(lv), intent(inout) :: FZ, FR, Fs
            real(dp), dimension(lv), intent(in)    :: aZ, aR, as, bZ, bR, bs
            integer :: i
            real(dp) :: c, w1, w2
            call dtheta(k_aZ, bZ)
            call dtheta(k_aR, bR)
            call dtheta(k_as, bs)
            do i = 1, lx1*ly1*lz1*nelv
               c  = lambda_t*as(i)*invR_f(i)*bm1(i,1,1,1)
               w1 = c*vtrans(i,1,1,1,1)
               w2 = c*vtrans(i,1,1,1,2)
               FZ(i) = FZ(i) + w1*(k_aZ(i) + bR(i))
               FR(i) = FR(i) + w1*(k_aR(i) - bZ(i))
               Fs(i) = Fs(i) + w2*k_as(i)
            end do
         end subroutine add_torsion_advection

         subroutine add_torsion_viscous(FZ, FR, Fs, vZ, vR, vs)
      !! G3(v):  F += mu * [torsion part of (D_s + M)^2 v] at P = 0, weak
      !! (l = lambda):
      !!   bm1 * ( 2 l^2/R^2 (dth v_R, -dth v_Z, 0) + 2 l/R^2 (0, dth v_s, -dth v_R)
      !!         + l/R^2 (v_s, 0, v_Z) + l zeta/R^3 (-l v_R, l v_Z - v_s, v_R)
      !!         - l^2/R^2 (v_Z, v_R, 0) )
      !!   + l^2 L_theta v, weak symmetric form  - l^2 dth^T[(bm1/R^2) dth v]
            real(dp), dimension(lv), intent(inout) :: FZ, FR, Fs
            real(dp), dimension(lv), intent(in)    :: vZ, vR, vs
            integer :: i
            real(dp) :: l, l2, mu1, mu2
            l  = lambda_t
            l2 = l*l
            call dtheta(k_dZ, vZ)
            call dtheta(k_dR, vR)
            call dtheta(k_ds, vs)
            do i = 1, lx1*ly1*lz1*nelv
               mu1 = vdiff(i,1,1,1,1)*bm1(i,1,1,1)
               mu2 = vdiff(i,1,1,1,2)*bm1(i,1,1,1)
               FZ(i) = FZ(i) + mu1*( 2.0_dp*l2*iR2_f(i)*k_dR(i) + l*iR2_f(i)*vs(i)
     &                             - l2*zR3_f(i)*vR(i) - l2*iR2_f(i)*vZ(i) )
               FR(i) = FR(i) + mu1*(-2.0_dp*l2*iR2_f(i)*k_dZ(i) + 2.0_dp*l*iR2_f(i)*k_ds(i)
     &                             + l2*zR3_f(i)*vZ(i) - l*zR3_f(i)*vs(i) - l2*iR2_f(i)*vR(i) )
               Fs(i) = Fs(i) + mu2*(-2.0_dp*l*iR2_f(i)*k_dR(i) + l*iR2_f(i)*vZ(i)
     &                             + l*zR3_f(i)*vR(i) )
            end do
            call add_ltheta(FZ, k_dZ, 1)
            call add_ltheta(FR, k_dR, 1)
            call add_ltheta(Fs, k_ds, 2)
         end subroutine add_torsion_viscous

         subroutine add_ltheta(F, dv, ifld)
      !! F += -lambda^2 dth^T[(mu bm1/R^2) dv],  dv = dth v: the lambda^2 L_theta
      !! block of G3, symmetric (so also its own transpose). ifld: 1 for Z, R,
      !! 2 for s (vdiff of the field).
            real(dp), dimension(lv), intent(inout) :: F
            real(dp), dimension(lv), intent(in)    :: dv
            integer, intent(in) :: ifld
            integer :: i
            do i = 1, lx1*ly1*lz1*nelv
               k_t1(i) = -lambda_t*lambda_t*vdiff(i,1,1,1,ifld)*bm1(i,1,1,1)*iR2_f(i)*dv(i)
            end do
            call dtheta_t(k_t2, k_t1)
            call add2(F, k_t2, lx1*ly1*lz1*nelv)
         end subroutine add_ltheta

         subroutine add_alpha_torsion_viscous(FZ, FR, Fs, vZ, vR, vs, coef)
      !! G6: the alpha*lambda viscous coupling, the P*lambda part of
      !! mu (D_s + M)^2 v with P = -i alpha (helix_2Dh_derivation.tex, eq. lapexp):
      !!   mu * i * alpha lambda [ (2/R^2) dth v - (zeta/R^3) v + (2/R^2)(v_R, -v_Z, 0) ].
      !! Adds the REAL operator, F += coef * mu * bm1 * lambda [...], with
      !! coef = +-alpha chosen by the caller from the +i slot rule
      !! (Re[i w] = -Im w, Im[i w] = Re w).
            real(dp), dimension(lv), intent(inout) :: FZ, FR, Fs
            real(dp), dimension(lv), intent(in)    :: vZ, vR, vs
            real(dp), intent(in) :: coef
            integer :: i
            real(dp) :: c1, c2
            if (coef == 0.0_dp) return
            call dtheta(k_dZ, vZ)
            call dtheta(k_dR, vR)
            call dtheta(k_ds, vs)
            do i = 1, lx1*ly1*lz1*nelv
               c1  = coef*lambda_t*vdiff(i,1,1,1,1)*bm1(i,1,1,1)
               c2  = coef*lambda_t*vdiff(i,1,1,1,2)*bm1(i,1,1,1)
               FZ(i) = FZ(i) + c1*( 2.0_dp*iR2_f(i)*(k_dZ(i) + vR(i)) - zR3_f(i)*vZ(i) )
               FR(i) = FR(i) + c1*( 2.0_dp*iR2_f(i)*(k_dR(i) - vZ(i)) - zR3_f(i)*vR(i) )
               Fs(i) = Fs(i) + c2*( 2.0_dp*iR2_f(i)*k_ds(i)           - zR3_f(i)*vs(i) )
            end do
         end subroutine add_alpha_torsion_viscous

         subroutine torsion_force_nl(FZ, FR, Fs, UZ, UR, Us)
      !! Nonlinear torsion force (weak): G2(U,U) + G3(U).
            real(dp), dimension(lv), intent(out) :: FZ, FR, Fs
            real(dp), dimension(lv), intent(in)  :: UZ, UR, Us
            integer :: ntot1
            ntot1 = lx1*ly1*lz1*nelv
            call rzero(FZ, ntot1); call rzero(FR, ntot1); call rzero(Fs, ntot1)
            if (lambda_eq_zero) return
            call add_torsion_advection(FZ, FR, Fs, UZ, UR, Us, UZ, UR, Us)
            call add_torsion_viscous  (FZ, FR, Fs, UZ, UR, Us)
         end subroutine torsion_force_nl

         subroutine torsion_force_lin(FZ, FR, Fs, UZ, UR, Us, uZp, uRp, usp)
      !! Linearised torsion force (weak): G2(U,u) + G2(u,U) + G3(u). Exact
      !! derivative of torsion_force_nl by construction (G2 bilinear, G3 linear).
            real(dp), dimension(lv), intent(out) :: FZ, FR, Fs
            real(dp), dimension(lv), intent(in)  :: UZ, UR, Us, uZp, uRp, usp
            integer :: ntot1
            ntot1 = lx1*ly1*lz1*nelv
            call rzero(FZ, ntot1); call rzero(FR, ntot1); call rzero(Fs, ntot1)
            if (lambda_eq_zero) return
            call add_torsion_advection(FZ, FR, Fs, UZ, UR, Us, uZp, uRp, usp)
            call add_torsion_advection(FZ, FR, Fs, uZp, uRp, usp, UZ, UR, Us)
            call add_torsion_viscous  (FZ, FR, Fs, uZp, uRp, usp)
         end subroutine torsion_force_lin

      !====================================================================
      !     ADJOINT: EXACT TRANSPOSES OF THE WEAK FORCES
      !====================================================================
      !
      ! Each routine adds F^T(w) for the corresponding direct kernel:
      !   sum_i w_i F_i(u) = sum_i F^T(w)_i u_i      (all local nodes, all components)
      ! The row weight of the direct force (rho bm1, mu bm1) multiplies w first;
      ! dth in a row becomes dth^T of that product in the column.

         subroutine add_torsion_advection_lin_adj(FZ, FR, Fs, UZ, UR, Us, wZ, wR, ws)
      !! F += [G2(U,.) + G2(.,U)]^T w. Transpose of
      !!   u -> rho b c(U) S(u):   g = rho b c(U) w,
      !!        F_Z += dth^T g_Z - g_R,  F_R += dth^T g_R + g_Z,  F_s += dth^T g_s
      !!   u -> rho b (lambda u_s/R) S(U):
      !!        F_s += (lambda/R) b sum_k rho_k w_k S_k(U)
            real(dp), dimension(lv), intent(inout) :: FZ, FR, Fs
            real(dp), dimension(lv), intent(in)    :: UZ, UR, Us, wZ, wR, ws
            integer :: i, ntot1
            real(dp) :: c
            ntot1 = lx1*ly1*lz1*nelv
      ! transpose of c(U) S(u)
            do i = 1, ntot1
               c = lambda_t*Us(i)*invR_f(i)*bm1(i,1,1,1)
               k_gZ(i) = c*vtrans(i,1,1,1,1)*wZ(i)
               k_gR(i) = c*vtrans(i,1,1,1,1)*wR(i)
               k_gs(i) = c*vtrans(i,1,1,1,2)*ws(i)
               FZ(i) = FZ(i) - k_gR(i)
               FR(i) = FR(i) + k_gZ(i)
            end do
            call add_dtheta_t(FZ, k_gZ)
            call add_dtheta_t(FR, k_gR)
            call add_dtheta_t(Fs, k_gs)
      ! transpose of c(u) S(U)
            call dtheta(k_aZ, UZ)
            call dtheta(k_aR, UR)
            call dtheta(k_as, Us)
            do i = 1, ntot1
               c = lambda_t*invR_f(i)*bm1(i,1,1,1)
               Fs(i) = Fs(i) + c*( vtrans(i,1,1,1,1)*( wZ(i)*(k_aZ(i) + UR(i))
     &                                               + wR(i)*(k_aR(i) - UZ(i)) )
     &                           + vtrans(i,1,1,1,2)*ws(i)*k_as(i) )
            end do
         end subroutine add_torsion_advection_lin_adj

         subroutine add_torsion_viscous_adj(FZ, FR, Fs, wZ, wR, ws)
      !! F += G3^T w. With a = (mu1 b w_Z, mu1 b w_R, mu2 b w_s), l = lambda:
      !!   F_Z += -l^2/R^2 a_Z + l^2 zeta/R^3 a_R + l/R^2 a_s   + dth^T(-2 l^2/R^2 a_R)
      !!   F_R += -l^2 zeta/R^3 a_Z - l^2/R^2 a_R + l zeta/R^3 a_s
      !!                                                     + dth^T( 2 l^2/R^2 a_Z - 2 l/R^2 a_s)
      !!   F_s +=  l/R^2 a_Z - l zeta/R^3 a_R               + dth^T( 2 l/R^2 a_R)
      !! plus the lambda^2 L_theta block (symmetric).
            real(dp), dimension(lv), intent(inout) :: FZ, FR, Fs
            real(dp), dimension(lv), intent(in)    :: wZ, wR, ws
            integer :: i
            real(dp) :: l, l2, aZ, aR, as
            l  = lambda_t
            l2 = l*l
            do i = 1, lx1*ly1*lz1*nelv
               aZ = vdiff(i,1,1,1,1)*bm1(i,1,1,1)*wZ(i)
               aR = vdiff(i,1,1,1,1)*bm1(i,1,1,1)*wR(i)
               as = vdiff(i,1,1,1,2)*bm1(i,1,1,1)*ws(i)
               FZ(i) = FZ(i) - l2*iR2_f(i)*aZ + l2*zR3_f(i)*aR + l*iR2_f(i)*as
               FR(i) = FR(i) - l2*zR3_f(i)*aZ - l2*iR2_f(i)*aR + l*zR3_f(i)*as
               Fs(i) = Fs(i) + l*iR2_f(i)*aZ  - l*zR3_f(i)*aR
               k_gZ(i) = -2.0_dp*l2*iR2_f(i)*aR
               k_gR(i) =  2.0_dp*l2*iR2_f(i)*aZ - 2.0_dp*l*iR2_f(i)*as
               k_gs(i) =  2.0_dp*l*iR2_f(i)*aR
            end do
            call add_dtheta_t(FZ, k_gZ)
            call add_dtheta_t(FR, k_gR)
            call add_dtheta_t(Fs, k_gs)
            call dtheta(k_dZ, wZ)
            call dtheta(k_dR, wR)
            call dtheta(k_ds, ws)
            call add_ltheta(FZ, k_dZ, 1)
            call add_ltheta(FR, k_dR, 1)
            call add_ltheta(Fs, k_ds, 2)
         end subroutine add_torsion_viscous_adj

         subroutine add_alpha_torsion_viscous_adj(FZ, FR, Fs, wZ, wR, ws, coef)
      !! F += coef * K^T w, K the real G6 operator of add_alpha_torsion_viscous.
      !! With a = coef lambda (mu1 b w_Z, mu1 b w_R, mu2 b w_s):
      !!   F_Z += dth^T(2/R^2 a_Z) - zeta/R^3 a_Z - 2/R^2 a_R
      !!   F_R += dth^T(2/R^2 a_R) - zeta/R^3 a_R + 2/R^2 a_Z
      !!   F_s += dth^T(2/R^2 a_s) - zeta/R^3 a_s
      !! The caller picks coef from the TRANSPOSED slot block (see header).
            real(dp), dimension(lv), intent(inout) :: FZ, FR, Fs
            real(dp), dimension(lv), intent(in)    :: wZ, wR, ws
            real(dp), intent(in) :: coef
            integer :: i
            real(dp) :: aZ, aR, as
            if (coef == 0.0_dp) return
            do i = 1, lx1*ly1*lz1*nelv
               aZ = coef*lambda_t*vdiff(i,1,1,1,1)*bm1(i,1,1,1)*wZ(i)
               aR = coef*lambda_t*vdiff(i,1,1,1,1)*bm1(i,1,1,1)*wR(i)
               as = coef*lambda_t*vdiff(i,1,1,1,2)*bm1(i,1,1,1)*ws(i)
               FZ(i) = FZ(i) - zR3_f(i)*aZ - 2.0_dp*iR2_f(i)*aR
               FR(i) = FR(i) - zR3_f(i)*aR + 2.0_dp*iR2_f(i)*aZ
               Fs(i) = Fs(i) - zR3_f(i)*as
               k_gZ(i) = 2.0_dp*iR2_f(i)*aZ
               k_gR(i) = 2.0_dp*iR2_f(i)*aR
               k_gs(i) = 2.0_dp*iR2_f(i)*as
            end do
            call add_dtheta_t(FZ, k_gZ)
            call add_dtheta_t(FR, k_gR)
            call add_dtheta_t(Fs, k_gs)
         end subroutine add_alpha_torsion_viscous_adj

         subroutine add_dtheta_t(F, g)
      !! F += dth^T g
            real(dp), dimension(lv), intent(inout) :: F
            real(dp), dimension(lv), intent(in)    :: g
            call dtheta_t(k_h, g)
            call add2(F, k_h, lx1*ly1*lz1*nelv)
         end subroutine add_dtheta_t

         subroutine torsion_force_lin_adj(FZ, FR, Fs, UZ, UR, Us, wZ, wR, ws)
      !! Exact transpose of torsion_force_lin (same base U).
            real(dp), dimension(lv), intent(out) :: FZ, FR, Fs
            real(dp), dimension(lv), intent(in)  :: UZ, UR, Us, wZ, wR, ws
            integer :: ntot1
            ntot1 = lx1*ly1*lz1*nelv
            call rzero(FZ, ntot1); call rzero(FR, ntot1); call rzero(Fs, ntot1)
            if (lambda_eq_zero) return
            call add_torsion_advection_lin_adj(FZ, FR, Fs, UZ, UR, Us, wZ, wR, ws)
            call add_torsion_viscous_adj      (FZ, FR, Fs, wZ, wR, ws)
         end subroutine torsion_force_lin_adj

         subroutine torsion_prepare_adjoint()
      !! Adjoint step setup, from the base flow in (vx, vy): the in-plane
      !! divergence div_in U = dU_Z/dZ + dU_R/dR + U_R/R, collocated (gradm1).
      !! Called by prepare_step_2Dh_axisym at every adjoint step. No-op without
      !! torsion.
            integer :: ntot1
            if (lambda_eq_zero) return
            ntot1 = lx1*ly1*lz1*nelv
            call gradm1 (divU_f, k_h, k_gs, vx)    ! dU_Z/dZ
            call gradm1 (k_gZ, k_t2, k_gs, vy)     ! dU_R/dR
            call add2   (divU_f, k_t2, ntot1)
            call addcol3(divU_f, vy, invR_f, ntot1)
         end subroutine torsion_prepare_adjoint

         subroutine add_adjoint_transport_div(FZ, FR, Fs, wZ, wR, ws)
      !! F_k += rho_k bm1 (div_in U) w_k: the part of the adjoint of Nek's
      !! transport that Nek's adjoint leaves out (see the ADJOINT section).
            real(dp), dimension(lv), intent(inout) :: FZ, FR, Fs
            real(dp), dimension(lv), intent(in)    :: wZ, wR, ws
            integer :: i
            real(dp) :: q
            if (lambda_eq_zero) return
            do i = 1, lx1*ly1*lz1*nelv
               q = divU_f(i)*bm1(i,1,1,1)
               FZ(i) = FZ(i) + q*vtrans(i,1,1,1,1)*wZ(i)
               FR(i) = FR(i) + q*vtrans(i,1,1,1,1)*wR(i)
               Fs(i) = Fs(i) + q*vtrans(i,1,1,1,2)*ws(i)
            end do
         end subroutine add_adjoint_transport_div

      !====================================================================
      !     NONLINEAR STEP HOOKS
      !====================================================================

         subroutine torsion_refresh_nl()
      !! From U^n (vx, vy, t) before heat(1): the accelerations read by the
      !! registered userf/userq terms. Nek rebuilds F exactly: makeuf multiplies
      !! by bm1 and makeabf by rho; makeuq by bm1 and neklab_source by rhoCp.
            integer :: ntot1
            if (lambda_eq_zero) return
            ntot1 = lx1*ly1*lz1*nelv
            call torsion_force_nl(fZ_nl, fR_nl, fs_nl, vx, vy, t(1,1,1,1,1))
            call col3   (k_t1, vtrans(1,1,1,1,1), bm1, ntot1)
            call invcol2(fZ_nl, k_t1, ntot1)
            call invcol2(fR_nl, k_t1, ntot1)
            call col3   (k_t1, vtrans(1,1,1,1,2), bm1, ntot1)
            call invcol2(fs_nl, k_t1, ntot1)
         end subroutine torsion_refresh_nl

         subroutine torsion_forcing_pw(fx, fy, fz, ix, iy, iz, e)
            real(dp), intent(inout) :: fx, fy, fz
            integer, intent(in) :: ix, iy, iz, e
            integer :: i
            if (lambda_eq_zero) return
            i = ix + lx1*((iy - 1) + ly1*((iz - 1) + lz1*(e - 1)))
            fx = fx + fZ_nl(i)
            fy = fy + fR_nl(i)
         end subroutine torsion_forcing_pw

         subroutine torsion_source_pw(q, ix, iy, iz, e)
            real(dp), intent(inout) :: q
            integer, intent(in) :: ix, iy, iz, e
            integer :: i
            if (lambda_eq_zero) return
            i = ix + lx1*((iy - 1) + ly1*((iz - 1) + lz1*(e - 1)))
            q = q + fs_nl(i)
         end subroutine torsion_source_pw

         subroutine torsion_add_pressure_nl()
      !! After heat(1), before heat(2): bq(u_s) += D_lam^T p_ext, with p_ext the
      !! extrapolation that cresvif/extrapp will apply to the in-plane gradient
      !! in fluid(2) of this step (pr/prlag untouched). Added AFTER makeabq, so
      !! it is not extrapolated again.
            real(dp) :: cst
            integer :: ntot2
            if (lambda_eq_zero) return
            ntot2 = lx2*ly2*lz2*nelv
            if (nbd == 3) then
               cst = dtlag(1)/dtlag(2)
               call add3s2(n_pext, pr, prlag, 1.0_dp + cst, -cst, ntot2)   ! p + cst (p - p_lag)
            else if (nbd > 3) then
               call nek_stop_error('pressure extrapolation order > 3.', this_module,
     &            'torsion_add_pressure_nl')
            else
               call copy(n_pext, pr, ntot2)
            end if
            call dlam_grad(e_m1, n_pext)
            call add2(bq(1,1,1,1,1), e_m1, lx1*ly1*lz1*nelv)
         end subroutine torsion_add_pressure_nl

      !====================================================================
      !     LINEAR STEP HOOKS (slot jp)
      !====================================================================

         subroutine torsion_add_linear_uv(bfx_jp, bfy_jp, jp_)
      !! Velocity stage of makefp: bfxp/bfyp += L_Z, L_R of slot jp_; the u_s
      !! component is kept for torsion_add_linear_s.
      !!   direct:  L = G2/G3 linearised + G6 from the other slot
      !!            (slot 1 gets -alpha G6(u_Im), slot 2 gets +alpha G6(u_Re))
      !!   adjoint: L^T + (div_in U) v + the transposed G6 block
      !!            (slot 1 gets +alpha G6^T(v_Im), slot 2 gets -alpha G6^T(v_Re))
            real(dp), dimension(lv), intent(inout) :: bfx_jp, bfy_jp
            integer, intent(in) :: jp_
            integer :: ntot1, ipert
            if (lambda_eq_zero) return
            ntot1 = lx1*ly1*lz1*nelv
            ipert = npert + 1 - jp_
            if (.not. ifadj) then
               call torsion_force_lin(lin_FZ, lin_FR, lin_Fs, vx, vy, t(1,1,1,1,1),
     &                                vxp(1,jp_), vyp(1,jp_), tp(1,1,jp_))
               if (alpha_t /= 0.0_dp) call add_alpha_torsion_viscous(lin_FZ, lin_FR, lin_Fs,
     &            vxp(1,ipert), vyp(1,ipert), tp(1,1,ipert), merge(-alpha_t, alpha_t, jp_ == 1))
            else
               call torsion_force_lin_adj(lin_FZ, lin_FR, lin_Fs, vx, vy, t(1,1,1,1,1),
     &                                    vxp(1,jp_), vyp(1,jp_), tp(1,1,jp_))
               call add_adjoint_transport_div(lin_FZ, lin_FR, lin_Fs,
     &                                        vxp(1,jp_), vyp(1,jp_), tp(1,1,jp_))
               if (alpha_t /= 0.0_dp) call add_alpha_torsion_viscous_adj(lin_FZ, lin_FR, lin_Fs,
     &            vxp(1,ipert), vyp(1,ipert), tp(1,1,ipert), merge(alpha_t, -alpha_t, jp_ == 1))
            end if
            call add2(bfx_jp, lin_FZ, ntot1)
            call add2(bfy_jp, lin_FR, ntot1)
         end subroutine torsion_add_linear_uv

         subroutine torsion_add_linear_s(bq_jp)
      !! Scalar stage of makefp: bqp += L_s (computed by torsion_add_linear_uv).
            real(dp), dimension(lv), intent(inout) :: bq_jp
            if (lambda_eq_zero) return
            call add2(bq_jp, lin_Fs, lx1*ly1*lz1*nelv)
         end subroutine torsion_add_linear_s

      !====================================================================
      !     CFL
      !====================================================================

         subroutine torsion_transport_velocity(cx, cy, cz)
      !! In-plane velocity that sets the explicit advective limit:
      !!   W = (U_Z + c rho, U_R - c zeta, 0),  c = lambda U_s/R.
      !! The rotation part c dth is explicit and invisible to Nek's CFL. The
      !! adjoint is transported by the same W.
            real(dp), dimension(lv), intent(out) :: cx, cy, cz
            integer :: i
            real(dp) :: c
            do i = 1, lx1*ly1*lz1*nelv
               c = merge(0.0_dp, lambda_t*t(i,1,1,1,1)*invR_f(i), lambda_eq_zero)
               cx(i) = vx(i,1,1,1) + c*rho_f(i)
               cy(i) = vy(i,1,1,1) - c*zeta_f(i)
               cz(i) = 0.0_dp
            end do
         end subroutine torsion_transport_velocity

      end module neklab_2Dh_torsion