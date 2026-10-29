      module neklab_2Dh_torsion
      !! Torsion terms of the 2Dh axisymmetric solver: helical pipe in the
      !! Germano (Frenet) frame. Written from scratch (restart, 28 Oct 2026)
      !! from helix_2Dh_derivation.tex; nothing here depends on either
      !! stepper. The steppers (neklab_2Dh, neklab_2Dh_axisym) call:
      !!
      !!   init_torsion_2Dh           once, from init_2Dh_axisym
      !!   torsion_refresh_nl         nonlinear step, from U^n, before heat(1)
      !!   torsion_add_pressure_nl    nonlinear step, after heat(1): bq += D_lam^T p_ext
      !!   torsion_add_linear_uv/_s   linear step, makefp: bfxp/bfyp/bqp += L(U; u)
      !!   torsion_transport_velocity CFL hook (rotation advection c(U) dth)
      !!   dlam_div_add, dlam_grad_t  the continuity/pressure pair, both steppers
      !!   lambda_block_E_add         the lambda block of the pressure operator
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
      !! Terms Nek does not provide, as right-hand-side force densities:
      !!
      !!   G2  advection + frame   rho_d c(a) S(b),   c(a) = lambda a_s / R
      !!                           S(b) = (dth b_Z + b_R, dth b_R - b_Z, dth b_s)
      !!       nonlinear: G2(U,U);   linear: G2(U,u) + G2(u,U)   (exact by construction)
      !!   G3  viscous             mu * [torsion part of (D_s + M)^2 v] at P = 0
      !!                           (linear in v; lambda^2/R^2 shifts on Z, R included)
      !!   G6  alpha*lambda (linear step, alpha /= 0): the P*lambda part of
      !!       mu (D_s + M)^2 u,  mu i alpha lambda [(2/R^2) dth u - (zeta/R^3) u
      !!       + (2/R^2)(u_R, -u_Z, 0)]; odd in P, so it swaps slots (+i rule)
      !!   G4  pressure, u_s eq.   +(lambda/R) dth p   =  weak  D_lam^T p
      !!   G5  continuity          -(lambda/R) dth u_s =  weak  D_lam u_s
      !!
      !! G2 and G3 are computed ONCE, as weak (bm1-weighted) force vectors, by the
      !! same routines for both solvers. The linear solver adds them to
      !! bfxp/bfyp/bqp; the nonlinear one passes F/(rho bm1) through userf/userq,
      !! which Nek multiplies back by bm1 (makeuf/makeuq) and rho (makeabf /
      !! neklab_source). Both go through EXT3 at the same time level.
      !!
      !! G4/G5 are IMPLICIT in both solvers: D_lam is part of the projection
      !! (E = D B^-1 D^T with D = [D_in, D_lam]) and the correction updates u_s.
      !! An explicit D_lam^T P (route A) has a dt-independent pressure-error
      !! instability for (lambda delta/(1-delta))^2 > 1/7 with EXT3.
      !!
      !! The D_lam pair is made exactly mean-free (1^T D_lam u = 0, D_lam^T 1 = 0),
      !! so that E 1 = 0 holds exactly: ortho and the hsmg preconditioner assume it.
      !!
      !!--------------------------------------------------------------------
      !! DISCRETE dth
      !!--------------------------------------------------------------------
      !!
      !! Element-local collocation derivative (as gradm1):
      !!   dth f = A_r f_r + A_s f_s ,
      !!   A_r = (zeta ry - rho rx)/J ,  A_s = (zeta sy - rho sx)/J
      !! dth_t is its exact algebraic transpose, never "-dth". The meshes must not
      !! touch the axis (R > 0 is checked by init_2Dh_axisym).
      !!
      !!--------------------------------------------------------------------
      !! TERM-GROUP SWITCHES
      !!--------------------------------------------------------------------
      !!
      !! set_torsion_groups(adv, visc, pres) acts on BOTH solvers, so that a
      !! partial configuration is still a consistent nonlinear/linear pair (the
      !! Jacobian tests remain valid). 'pres' removes D_lam everywhere at once:
      !! continuity, E, gradient and correction. Call before init_2Dh_axisym.
         use LightKrylov, only: dp
         use stdlib_optval, only: optval
         use neklab_nek_setup, only: nek_log_message, nek_log_warning, nek_log_debug,
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

      ! --- state
         logical,  private :: torsion_on    = .false.   ! lambda /= 0 and init done
         logical,  private :: torsion_built = .false.
         real(dp), private :: lamt  = 0.0_dp
         real(dp), private :: alph  = 0.0_dp            ! alpha of the linear step (G6)
         real(dp), private :: zcen   = 0.0_dp
         real(dp), private :: Rcen   = 0.0_dp
         real(dp), private :: vol2 = 0.0_dp            ! sum(bm2), mean removal
         logical,  private :: use_adv  = .true.
         logical,  private :: use_visc = .true.
         logical,  private :: use_pres = .true.

      ! --- coefficient fields (M1)
         real(dp), dimension(lv), private :: zeta_f, rho_f, invR_f
         real(dp), dimension(lv), private :: Ar_f, As_f   ! dth = Ar d/dr + As d/ds

      ! --- nonlinear forcing as accelerations, read pointwise by the registered terms
         real(dp), dimension(lv), private :: fZ_nl = 0.0_dp, fR_nl = 0.0_dp, fs_nl = 0.0_dp

      ! --- linear forcing of the current slot, computed at the velocity stage of
      !     makefp and consumed (u_s part) at its scalar stage
         real(dp), dimension(lv), private :: lin_FZ, lin_FR, lin_Fs

      ! --- scratch, one group per call level
         real(dp), dimension(lv), private :: k_dZ, k_dR, k_ds, k_t1, k_t2   ! kernels
         real(dp), dimension(lv), private :: k_aZ, k_aR, k_as               ! kernels (2nd arg)
         real(dp), dimension(lv), private :: p_m1a, p_m1b                   ! pressure pair
         real(dp), dimension(lp), private :: p_m2a                          ! pressure pair
         real(dp), dimension(lv), private :: e_m1                           ! E block, correction
         real(dp), dimension(lp), private :: e_m2                           ! E block
         real(dp), dimension(lp), private :: n_pext                         ! nonlinear p_ext
         real(dp), dimension(lx1,ly1), private :: t_r, t_s                  ! dth element work

         public :: init_torsion_2Dh
         public :: set_torsion_groups
         public :: torsion_active
         public :: torsion_lambda
         public :: dtheta, dtheta_t
         public :: dlam_div_add, dlam_grad_t
         public :: lambda_block_E_add, lambda_correction
         public :: torsion_force_nl, torsion_force_lin
         public :: torsion_refresh_nl, torsion_add_pressure_nl
         public :: torsion_add_linear_uv, torsion_add_linear_s
         public :: torsion_transport_velocity
         public :: check_wall_centering
         public :: set_torsion_alpha
         public :: add_g2, add_g3, add_g6  ! kernels, exposed for the unit tests

      contains

      !====================================================================
      !     SETUP
      !====================================================================

         subroutine set_torsion_groups(adv, visc, pres)
      !! Switch term groups on/off in BOTH solvers. Default: all on.
            logical, optional, intent(in) :: adv, visc, pres
            character(len=128) :: msg
            if (torsion_built) call nek_stop_error('set_torsion_groups must be called before '//
     &         'init_2Dh_axisym.', this_module, 'set_torsion_groups')
            use_adv  = optval(adv,  .true.)
            use_visc = optval(visc, .true.)
            use_pres = optval(pres, .true.)
            write (msg, '(A,3L2)') 'torsion groups (adv, visc, pres):', use_adv, use_visc, use_pres
            call nek_log_message(msg, this_module, 'set_torsion_groups')
         end subroutine set_torsion_groups

         subroutine set_torsion_alpha(alpha)
      !! Wavenumber of the linear step, for the alpha*lambda viscous coupling G6.
      !! Set by build_axisym_coeffs (neklab_2Dh_axisym).
            real(dp), intent(in) :: alpha
            alph = alpha
         end subroutine set_torsion_alpha

         logical function torsion_active()
            torsion_active = torsion_on
         end function torsion_active

         real(dp) function torsion_lambda()
            torsion_lambda = merge(lamt, 0.0_dp, torsion_on)
         end function torsion_lambda

         subroutine init_torsion_2Dh()
      !! Build the coefficients from t2Dh%init_geom and register the nonlinear
      !! terms and the CFL velocity. Idempotent. Called by init_2Dh_axisym, after
      !! the masks are built. A torus (lambda = 0, or no init_geom) leaves
      !! everything off: every entry point below then returns at once.
            character(len=*), parameter :: this_procedure = 'init_torsion_2Dh'
            character(len=256) :: msg
            real(dp), external :: glsum
            integer :: i, ie, ix, iy, ntot1
            real(dp) :: dev

            if (torsion_built) return
            torsion_built = .true.
            torsion_on = .false.
            if (t2Dh%is_torsion_zero()) then
               call nek_log_message('lambda = 0: torsion terms inactive.', this_module, this_procedure)
               return
            end if
            if (ldim /= 2) call nek_stop_error('torsion requires a 2D mesh.', this_module, this_procedure)
            if (ifadj) call nek_stop_error('torsion: adjoint not implemented.', this_module, this_procedure)

            lamt = t2Dh%get_lambda()
            zcen  = t2Dh%get_axial_centre()
            Rcen  = t2Dh%get_curv_radius()
            ntot1 = lx1*ly1*lz1*nelv

            do i = 1, ntot1
               zeta_f(i) = xm1(i,1,1,1) - zcen
               rho_f(i)  = ym1(i,1,1,1) - Rcen
               invR_f(i) = 1.0_dp/ym1(i,1,1,1)
            end do
            do ie = 1, nelv
               do iy = 1, ly1
               do ix = 1, lx1
                  i = ix + lx1*(iy - 1) + lx1*ly1*(ie - 1)
                  Ar_f(i) = jacmi(ix + lx1*(iy - 1),ie)*(zeta_f(i)*rym1(ix,iy,1,ie) - rho_f(i)*rxm1(ix,iy,1,ie))
                  As_f(i) = jacmi(ix + lx1*(iy - 1),ie)*(zeta_f(i)*sym1(ix,iy,1,ie) - rho_f(i)*sxm1(ix,iy,1,ie))
               end do
               end do
            end do
            vol2 = glsum(bm2, lx2*ly2*lz2*nelv)

            call check_full_section(this_procedure)
            call check_masks_equal(this_procedure)
            call check_wall_centering(dev, verbose=.true.)

            call add_neklab_forcing_term(torsion_forcing_pw, 'torsion')
            call add_neklab_source_term (torsion_source_pw,  'torsion')
            call set_cfl_velocity_hook(torsion_transport_velocity)

            torsion_on = .true.
            write (msg, '(A,ES14.6,A,ES14.6,A,ES14.6)') 'torsion active: lambda =', lamt,
     &         ', z_c =', zcen, ', R_c =', Rcen
            call nek_log_message(msg, this_module, this_procedure)
            write (msg, '(A,3L2)') 'torsion groups (adv, visc, pres):', use_adv, use_visc, use_pres
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
            real(dp), dimension(lx1,ly1,lelv), intent(out) :: out
            real(dp), dimension(lx1,ly1,lelv), intent(in)  :: f
            integer :: ie, ix, iy, i0
            do ie = 1, nelv
               call mxm(dxm1, lx1, f(1,1,ie), lx1, t_r, ly1)      ! f_r
               call mxm(f(1,1,ie), lx1, dytm1, ly1, t_s, ly1)     ! f_s
               i0 = lx1*ly1*(ie - 1)
               do iy = 1, ly1
               do ix = 1, lx1
                  out(ix,iy,ie) = Ar_f(i0 + ix + lx1*(iy - 1))*t_r(ix,iy)
     &                          + As_f(i0 + ix + lx1*(iy - 1))*t_s(ix,iy)
               end do
               end do
            end do
         end subroutine dtheta

         subroutine dtheta_t(out, g)
      !! out = dth^T g = D_r^T (A_r g) + D_s^T (A_s g), element-local.
            real(dp), dimension(lx1,ly1,lelv), intent(out) :: out
            real(dp), dimension(lx1,ly1,lelv), intent(in)  :: g
            real(dp), dimension(lx1,ly1) :: w2
            integer :: ie, ix, iy, i0
            do ie = 1, nelv
               i0 = lx1*ly1*(ie - 1)
               do iy = 1, ly1
               do ix = 1, lx1
                  t_r(ix,iy) = Ar_f(i0 + ix + lx1*(iy - 1))*g(ix,iy,ie)
                  t_s(ix,iy) = As_f(i0 + ix + lx1*(iy - 1))*g(ix,iy,ie)
               end do
               end do
               call mxm(dxtm1, lx1, t_r, lx1, out(1,1,ie), ly1)   ! D_r^T
               call mxm(t_s, lx1, dym1, ly1, w2, ly1)             ! D_s^T
               out(:,:,ie) = out(:,:,ie) + w2
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

         subroutine dlam_div_add(d, u)
      !! d += D_lam u   (d on M2, u on M1). No-op without torsion or with pres off.
            real(dp), dimension(lp), intent(inout) :: d
            real(dp), dimension(lv), intent(in)    :: u
            real(dp), external :: glsum
            real(dp) :: s
            integer :: ntot1, ntot2
            if (.not. (torsion_on .and. use_pres)) return
            ntot1 = lx1*ly1*lz1*nelv
            ntot2 = lx2*ly2*lz2*nelv
            call dtheta(p_m1a, u)
            p_m1a(1:ntot1) = -lamt*invR_f(1:ntot1)*p_m1a(1:ntot1)
            call bm2_map12_2Dh(p_m2a, p_m1a)
            s = glsum(p_m2a, ntot2)/vol2
            call add2s2(p_m2a, bm2, -s, ntot2)
            call add2(d, p_m2a, ntot2)
         end subroutine dlam_div_add

         subroutine dlam_grad_t(g, p)
      !! g = D_lam^T p   (g on M1, element-local, unassembled; p on M2).
      !! Zero without torsion or with pres off.
            real(dp), dimension(lv), intent(out) :: g
            real(dp), dimension(lp), intent(in)  :: p
            real(dp), external :: glsc2
            real(dp) :: pm
            integer :: ntot1, ntot2
            ntot1 = lx1*ly1*lz1*nelv
            ntot2 = lx2*ly2*lz2*nelv
            if (.not. (torsion_on .and. use_pres)) then
               call rzero(g, ntot1)
               return
            end if
            pm = glsc2(bm2, p, ntot2)/vol2
            call copy(p_m2a, p, ntot2)
            call cadd(p_m2a, -pm, ntot2)
            call bm2_map12T_2Dh(p_m1b, p_m2a)                 ! I12^T (bm2 (p - pm))
            p_m1b(1:ntot1) = -lamt*invR_f(1:ntot1)*p_m1b(1:ntot1)
            call dtheta_t(g, p_m1b)
         end subroutine dlam_grad_t

         subroutine lambda_block_E_add(ap, p, h2inv, mask)
      !! ap += D_lam (mask h2inv B^-1) D_lam^T p : the lambda block of the
      !! pressure operator E = D (h2inv B^-1) D^T. Symmetric, positive
      !! semi-definite, annihilates constants.
            real(dp), dimension(lp), intent(inout) :: ap
            real(dp), dimension(lp), intent(in)    :: p
            real(dp), dimension(lv), intent(in)    :: h2inv, mask
            integer :: ntot1
            if (.not. (torsion_on .and. use_pres)) return
            call lambda_correction(e_m1, p, h2inv, mask)
            call dlam_div_add(ap, e_m1)
         end subroutine lambda_block_E_add

         subroutine lambda_correction(dw, dp_, h2inv, mask)
      !! dw = mask h2inv B^-1 D_lam^T dp : the u_s pressure correction, same
      !! sign and same B^-1 as the in-plane one (opgradt -> opbinv).
            real(dp), dimension(lv), intent(out) :: dw
            real(dp), dimension(lp), intent(in)  :: dp_
            real(dp), dimension(lv), intent(in)  :: h2inv, mask
            integer :: ntot1
            ntot1 = lx1*ly1*lz1*nelv
            call dlam_grad_t(dw, dp_)
            if (.not. (torsion_on .and. use_pres)) return
            call col2 (dw, mask, ntot1)
            call dssum(dw, lx1, ly1, lz1)
            call col2 (dw, binvm1, ntot1)
            call col2 (dw, h2inv, ntot1)
         end subroutine lambda_correction

      !====================================================================
      !     G2 AND G3: WEAK FORCE VECTORS
      !====================================================================

         subroutine add_g2(FZ, FR, Fs, aZ, aR, as, bZ, bR, bs)
      !! F += rho bm1 c(a) S(b),  c(a) = lambda a_s/R,
      !! S(b) = (dth b_Z + b_R, dth b_R - b_Z, dth b_s).
      !! dth b is element-local: fine, the result is a local weak contribution.
            real(dp), dimension(lv), intent(inout) :: FZ, FR, Fs
            real(dp), dimension(lv), intent(in)    :: aZ, aR, as, bZ, bR, bs
            integer :: i, ntot1
            real(dp) :: c, w1, w2
            if (.not. use_adv) return
            ntot1 = lx1*ly1*lz1*nelv
            call dtheta(k_aZ, bZ)
            call dtheta(k_aR, bR)
            call dtheta(k_as, bs)
            do i = 1, ntot1
               c  = lamt*as(i)*invR_f(i)*bm1(i,1,1,1)
               w1 = c*vtrans(i,1,1,1,1)
               w2 = c*vtrans(i,1,1,1,2)
               FZ(i) = FZ(i) + w1*(k_aZ(i) + bR(i))
               FR(i) = FR(i) + w1*(k_aR(i) - bZ(i))
               Fs(i) = Fs(i) + w2*k_as(i)
            end do
         end subroutine add_g2

         subroutine add_g3(FZ, FR, Fs, vZ, vR, vs)
      !! F += mu * [torsion part of (D_s + M)^2 v] at P = 0, weak:
      !!   pointwise: bm1 * ( 2 l^2/R^2 (dth v_R, -dth v_Z, 0) + 2 l/R^2 (0, dth v_s, -dth v_R)
      !!                    + l/R^2 (v_s, 0, v_Z) + l z/R^3 (-l v_R, l v_Z - v_s, v_R)
      !!                    - l^2/R^2 (v_Z, v_R, 0) )
      !!   l^2 L_theta, weak symmetric form:  - l^2 dth^T[ (bm1/R^2) dth v_c ]
            real(dp), dimension(lv), intent(inout) :: FZ, FR, Fs
            real(dp), dimension(lv), intent(in)    :: vZ, vR, vs
            integer :: i, ntot1
            real(dp) :: iR, iR2, l2R2, lR2, lzR3, mu1, mu2, b
            if (.not. use_visc) return
            ntot1 = lx1*ly1*lz1*nelv
            call dtheta(k_dZ, vZ)
            call dtheta(k_dR, vR)
            call dtheta(k_ds, vs)
            do i = 1, ntot1
               iR   = invR_f(i)
               iR2  = iR*iR
               l2R2 = lamt*lamt*iR2
               lR2  = lamt*iR2
               lzR3 = lamt*zeta_f(i)*iR2*iR
               b    = bm1(i,1,1,1)
               mu1  = vdiff(i,1,1,1,1)*b
               mu2  = vdiff(i,1,1,1,2)*b
               FZ(i) = FZ(i) + mu1*( 2.0_dp*l2R2*k_dR(i) + lR2*vs(i)
     &                             - lamt*lzR3*vR(i) - l2R2*vZ(i) )
               FR(i) = FR(i) + mu1*(-2.0_dp*l2R2*k_dZ(i) + 2.0_dp*lR2*k_ds(i)
     &                             + lamt*lzR3*vZ(i) - lzR3*vs(i) - l2R2*vR(i) )
               Fs(i) = Fs(i) + mu2*(-2.0_dp*lR2*k_dR(i) + lR2*vZ(i) + lzR3*vR(i) )
            end do
      ! lambda^2 L_theta: - lambda^2 mu dth^T[(bm1/R^2) dth v]
            do i = 1, ntot1
               b = -lamt*lamt*bm1(i,1,1,1)*invR_f(i)*invR_f(i)
               k_t1(i) = b*vdiff(i,1,1,1,1)*k_dZ(i)
            end do
            call dtheta_t(k_t2, k_t1)
            call add2(FZ, k_t2, ntot1)
            do i = 1, ntot1
               b = -lamt*lamt*bm1(i,1,1,1)*invR_f(i)*invR_f(i)
               k_t1(i) = b*vdiff(i,1,1,1,1)*k_dR(i)
            end do
            call dtheta_t(k_t2, k_t1)
            call add2(FR, k_t2, ntot1)
            do i = 1, ntot1
               b = -lamt*lamt*bm1(i,1,1,1)*invR_f(i)*invR_f(i)
               k_t1(i) = b*vdiff(i,1,1,1,2)*k_ds(i)
            end do
            call dtheta_t(k_t2, k_t1)
            call add2(Fs, k_t2, ntot1)
         end subroutine add_g3

         subroutine add_g6(FZ, FR, Fs, vZ, vR, vs, coef)
      !! alpha*lambda viscous coupling, the P*lambda part of mu (D_s + M)^2 v with
      !! P = -i alpha (helix_2Dh_derivation.tex, eq. lapexp):
      !!   mu * i * alpha lambda [ (2/R^2) dth v - (zeta/R^3) v + (2/R^2)(v_R, -v_Z, 0) ].
      !! This routine adds the REAL coefficient, F += coef * mu * bm1 * lambda [...],
      !! with coef = +-alpha chosen by the caller from the +i slot rule
      !! (Re[i w] = -Im w, Im[i w] = Re w).
            real(dp), dimension(lv), intent(inout) :: FZ, FR, Fs
            real(dp), dimension(lv), intent(in)    :: vZ, vR, vs
            real(dp), intent(in) :: coef
            integer :: i, ntot1
            real(dp) :: iR2, zR3, c1, c2
            if (.not. use_visc) return
            if (coef == 0.0_dp) return
            ntot1 = lx1*ly1*lz1*nelv
            call dtheta(k_dZ, vZ)
            call dtheta(k_dR, vR)
            call dtheta(k_ds, vs)
            do i = 1, ntot1
               iR2 = invR_f(i)*invR_f(i)
               zR3 = zeta_f(i)*iR2*invR_f(i)
               c1  = coef*lamt*vdiff(i,1,1,1,1)*bm1(i,1,1,1)
               c2  = coef*lamt*vdiff(i,1,1,1,2)*bm1(i,1,1,1)
               FZ(i) = FZ(i) + c1*( 2.0_dp*iR2*k_dZ(i) - zR3*vZ(i) + 2.0_dp*iR2*vR(i) )
               FR(i) = FR(i) + c1*( 2.0_dp*iR2*k_dR(i) - zR3*vR(i) - 2.0_dp*iR2*vZ(i) )
               Fs(i) = Fs(i) + c2*( 2.0_dp*iR2*k_ds(i) - zR3*vs(i) )
            end do
         end subroutine add_g6

         subroutine torsion_force_nl(FZ, FR, Fs, UZ, UR, Us)
      !! Nonlinear torsion force (weak): G2(U,U) + G3(U).
            real(dp), dimension(lv), intent(out) :: FZ, FR, Fs
            real(dp), dimension(lv), intent(in)  :: UZ, UR, Us
            integer :: ntot1
            ntot1 = lx1*ly1*lz1*nelv
            call rzero(FZ, ntot1); call rzero(FR, ntot1); call rzero(Fs, ntot1)
            if (.not. torsion_on) return
            call add_g2(FZ, FR, Fs, UZ, UR, Us, UZ, UR, Us)
            call add_g3(FZ, FR, Fs, UZ, UR, Us)
         end subroutine torsion_force_nl

         subroutine torsion_force_lin(FZ, FR, Fs, UZ, UR, Us, uZp, uRp, usp)
      !! Linearised torsion force (weak): G2(U,u) + G2(u,U) + G3(u). Exact
      !! derivative of torsion_force_nl by construction (G2 bilinear, G3 linear).
            real(dp), dimension(lv), intent(out) :: FZ, FR, Fs
            real(dp), dimension(lv), intent(in)  :: UZ, UR, Us, uZp, uRp, usp
            integer :: ntot1
            ntot1 = lx1*ly1*lz1*nelv
            call rzero(FZ, ntot1); call rzero(FR, ntot1); call rzero(Fs, ntot1)
            if (.not. torsion_on) return
            call add_g2(FZ, FR, Fs, UZ, UR, Us, uZp, uRp, usp)
            call add_g2(FZ, FR, Fs, uZp, uRp, usp, UZ, UR, Us)
            call add_g3(FZ, FR, Fs, uZp, uRp, usp)
         end subroutine torsion_force_lin

      !====================================================================
      !     NONLINEAR STEP HOOKS
      !====================================================================

         subroutine torsion_refresh_nl()
      !! From U^n (vx, vy, t) before heat(1): the accelerations read by the
      !! registered userf/userq terms. Nek rebuilds F exactly: makeuf multiplies
      !! by bm1 and makeabf by rho; makeuq by bm1 and neklab_source by rhoCp.
            integer :: i, ntot1
            if (.not. torsion_on) return
            ntot1 = lx1*ly1*lz1*nelv
            call torsion_force_nl(fZ_nl, fR_nl, fs_nl, vx, vy, t(1,1,1,1,1))
            do i = 1, ntot1
               fZ_nl(i) = fZ_nl(i)/(vtrans(i,1,1,1,1)*bm1(i,1,1,1))
               fR_nl(i) = fR_nl(i)/(vtrans(i,1,1,1,1)*bm1(i,1,1,1))
               fs_nl(i) = fs_nl(i)/(vtrans(i,1,1,1,2)*bm1(i,1,1,1))
            end do
         end subroutine torsion_refresh_nl

         subroutine torsion_forcing_pw(fx, fy, fz, ix, iy, iz, e)
            real(dp), intent(inout) :: fx, fy, fz
            integer, intent(in) :: ix, iy, iz, e
            integer :: i
            if (.not. torsion_on) return
            i = ix + lx1*((iy - 1) + ly1*((iz - 1) + lz1*(e - 1)))
            fx = fx + fZ_nl(i)
            fy = fy + fR_nl(i)
         end subroutine torsion_forcing_pw

         subroutine torsion_source_pw(q, ix, iy, iz, e)
            real(dp), intent(inout) :: q
            integer, intent(in) :: ix, iy, iz, e
            integer :: i
            if (.not. torsion_on) return
            i = ix + lx1*((iy - 1) + ly1*((iz - 1) + lz1*(e - 1)))
            q = q + fs_nl(i)
         end subroutine torsion_source_pw

         subroutine torsion_add_pressure_nl()
      !! After heat(1), before heat(2): bq(u_s) += D_lam^T p_ext, with p_ext the
      !! extrapolation that cresvif/extrapp will apply to the in-plane gradient
      !! in fluid(2) of this step (copied, pr/prlag untouched). Added AFTER
      !! makeabq, so it is not extrapolated again.
            integer :: i, ntot1, ntot2
            real(dp) :: cst
            integer :: igeom
            common /cgeom/ igeom
            if (.not. (torsion_on .and. use_pres)) return
            ntot1 = lx1*ly1*lz1*nelv
            ntot2 = lx2*ly2*lz2*nelv
            if (nbd == 3) then
               cst = dtlag(1)/dtlag(2)
               call extrap_copy(n_pext, pr, prlag, cst, ntot2)
            else if (nbd > 3) then
               call nek_stop_error('pressure extrapolation order > 3.', this_module,
     &            'torsion_add_pressure_nl')
            else
               call copy(n_pext, pr, ntot2)
            end if
            call dlam_grad_t(e_m1, n_pext)
            call add2(bq(1,1,1,1,1), e_m1, ntot1)
         end subroutine torsion_add_pressure_nl

         subroutine extrap_copy(pe, p, pl, cst, n)
            integer, intent(in) :: n
            real(dp), intent(out) :: pe(n)
            real(dp), intent(in)  :: p(n), pl(n), cst
            pe = p + cst*(p - pl)
         end subroutine extrap_copy

      !====================================================================
      !     LINEAR STEP HOOKS (slot jp)
      !====================================================================

         subroutine torsion_add_linear_uv(bfx_jp, bfy_jp, jp_)
      !! Velocity stage of makefp: bfxp/bfyp += L_Z, L_R of slot jp_; the u_s
      !! component is kept for torsion_add_linear_s.
            real(dp), dimension(lv), intent(inout) :: bfx_jp, bfy_jp
            integer, intent(in) :: jp_
            integer :: ntot1, ipert
            if (.not. torsion_on) return
            ntot1 = lx1*ly1*lz1*nelv
            call torsion_force_lin(lin_FZ, lin_FR, lin_Fs, vx, vy, t(1,1,1,1,1),
     &                             vxp(1,jp_), vyp(1,jp_), tp(1,1,jp_))
      ! alpha*lambda coupling (G6) from the other slot: slot 1 (Re) receives
      ! -alpha*G6(u_Im), slot 2 (Im) receives +alpha*G6(u_Re)
            if (alph /= 0.0_dp) then
               ipert = npert + 1 - jp_
               call add_g6(lin_FZ, lin_FR, lin_Fs, vxp(1,ipert), vyp(1,ipert), tp(1,1,ipert),
     &                     merge(-alph, alph, jp_ == 1))
            end if
            call add2(bfx_jp, lin_FZ, ntot1)
            call add2(bfy_jp, lin_FR, ntot1)
         end subroutine torsion_add_linear_uv

         subroutine torsion_add_linear_s(bq_jp)
      !! Scalar stage of makefp: bqp += L_s (computed by torsion_add_linear_uv).
            real(dp), dimension(lv), intent(inout) :: bq_jp
            if (.not. torsion_on) return
            call add2(bq_jp, lin_Fs, lx1*ly1*lz1*nelv)
         end subroutine torsion_add_linear_s

      !====================================================================
      !     CFL
      !====================================================================

         subroutine torsion_transport_velocity(cx, cy, cz)
      !! In-plane velocity that sets the explicit advective limit:
      !!   a = (U_Z + c rho, U_R - c zeta, 0),  c = lambda U_s/R.
      !! The rotation part c dth is explicit and invisible to Nek's CFL.
            real(dp), dimension(lv), intent(out) :: cx, cy, cz
            integer :: i, ntot1
            real(dp) :: c
            ntot1 = lx1*ly1*lz1*nelv
            do i = 1, ntot1
               c = merge(lamt*t(i,1,1,1,1)*invR_f(i), 0.0_dp, torsion_on)
               cx(i) = vx(i,1,1,1) + c*rho_f(i)
               cy(i) = vy(i,1,1,1) - c*zeta_f(i)
               cz(i) = 0.0_dp
            end do
         end subroutine torsion_transport_velocity

      end module neklab_2Dh_torsion