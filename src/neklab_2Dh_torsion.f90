      module neklab_2Dh_torsion
      !! Torsion terms for helical-pipe base flows in the Germano frame,
      !! scoped to alpha = 0 (the nonlinear/base-flow path). See the
      !! derivation notes for the full alpha /= 0 formulation; this module
      !! deliberately implements only the subset of it that survives at
      !! alpha = 0, and is structured so that extending it later is additive.
      !!
      !!--------------------------------------------------------------------
      !! GEOMETRY
      !!--------------------------------------------------------------------
      !!
      !! lambda = tau/kappa (torsion-to-curvature ratio) and the axial centre
      !! of the cross-section are properties of the MESH, set once via
      !! t2Dh%init_geom, never measured: a 2D cross-sectional mesh carries no
      !! information about pitch, so there is no mesh-based fallback the way
      !! there is for delta. lambda = 0 (the default) recovers the torus
      !! exactly; every entry point below returns immediately in that case.
      !!
      !!--------------------------------------------------------------------
      !! THE dtheta OPERATOR
      !!--------------------------------------------------------------------
      !!
      !!    dtheta = zhat*d/dR - rhat*d/dZ ,   zhat = xm1-z_c, rhat = ym1-R_c
      !!
      !! generates rotation of the cross-section about the registered centre
      !! (z_c, R_c). It is exactly tangent to a circular wall centred there,
      !! so it needs no boundary data -- PROVIDED (z_c, R_c) match the mesh.
      !! If they don't, dtheta picks up a spurious wall-normal component
      !! proportional to the mis-registration and to the near-wall gradient of
      !! whatever it is applied to; check_wall_centering below is the
      !! diagnostic for that failure mode, and is a DIFFERENT check from
      !! check_dtheta_adjoint (which only verifies the transpose pair is
      !! self-consistent, and is blind to both being wrong by the same
      !! amount).
      !!
      !!--------------------------------------------------------------------
      !! WHAT SURVIVES AT alpha = 0
      !!--------------------------------------------------------------------
      !!
      !! Every term proportional to alpha (the ones that swap the real and
      !! imaginary perturbation slots in the general formulation) is simply
      !! ABSENT here, not guarded -- there is no slot bookkeeping in this
      !! module at all. What remains:
      !!
      !!   - two new diagonal Helmholtz shifts: +lambda^2/R^2 on the R and Z
      !!     equations (the phi shift is unchanged: diag(Omega^2)_phi = -1
      !!     regardless of lambda);
      !!   - explicit advection terms (linearised (U.grad)u + (u.grad)U along
      !!     the torsion-corrected streamwise direction), added exactly where
      !!     add_torus_curv_R/phi are added, so they ride EXT3 for free;
      !!   - explicit viscous terms (the alpha-free part of the vector
      !!     Laplacian's torsion correction);
      !!   - a genuinely new pressure/continuity coupling,
      !!     -(lambda/R)*dtheta(u_phi) in the divergence and its exact adjoint
      !!     in the streamwise momentum equation. This is the one piece with
      !!     no precedent in the torus code: compute_gradp_axisym,
      !!     compute_dw_axisym, compute_frc_div_axisym and
      !!     pressure_matvec_2Dh_axisym all currently return unconditionally
      !!     at alpha = 0, so this is a code path with no prior coverage.
      !!
      !! Deferred to the alpha /= 0 phase, and NOT implemented here: the
      !! alpha-proportional advection/viscous terms, the second shift argument
      !! of the coupled Helmholtz solve, and the slot-off-diagonal block of
      !! the pressure Schur complement. neklab_2Dh_axisym refuses alpha /= 0
      !! with lambda /= 0 rather than silently dropping those terms.
      !!
      !! The kernel layer (apply_dtheta, apply_dtheta_t, dphi_apply,
      !! dphi_apply_t) takes and returns plain field arrays with no jp/slot
      !! argument anywhere. That is deliberate: it is unchanged by the
      !! alpha /= 0 extension, which only adds a slot-aware wrapper around it.
         use LightKrylov, only: dp
         use neklab_2Dh, only: wmask, wmask_defined, build_wmask
         use neklab_t2Dh, only: t2Dh
         use neklab_nek_setup, only: nek_log_debug, nek_log_message,
     &                               nek_log_warning, nek_stop_error
         use neklab_nek_forcing, only: set_neklab_forcing
         implicit none
         include "SIZE"
         include "TOTAL"

         private
         character(len=*), parameter, private :: this_module = 'neklab_2Dh_torsion'

         logical, public :: if_torsion_zero = .true.
      !! Mirrors t2Dh%is_torsion_zero(), cached at build_torsion_coeffs time so
      !! every routine below can branch on a plain logical instead of going
      !! through the type-bound function on every call.

      ! --- coefficient fields, built once (lambda does not change at runtime).
      !     Rank-1 (flat) throughout, matching every point of use: they are
      !     always indexed by a single flat node number or passed to
      !     col2/col3/add2-style routines, never used as rank-4 arrays.
         integer, parameter, private :: lv = lx1*ly1*lz1*lelv
         real(dp), dimension(lv), private :: zhat_f, rhat_f
         real(dp), dimension(lv), private :: lamR_f     ! lambda/R
         real(dp), dimension(lv), private :: lamR2_f    ! lambda/R^2
         real(dp), dimension(lv), private :: lam2R2_f   ! lambda^2/R^2
         real(dp), dimension(lv), private :: lamzR3_f   ! lambda*zhat/R^3
         real(dp), dimension(lv), private :: invR_f     ! 1/R
         logical, private :: torsion_coeffs_defined = .false.

         public :: build_torsion_coeffs
         public :: torsion_shift_add
         public :: apply_dtheta
         public :: apply_dtheta_t
         public :: dphi_apply
         public :: dphi_apply_t
         public :: add_torsion_forcing_RZ
         public :: add_torsion_forcing_phi
         public :: build_torsion_forcing_nl

      contains

      !====================================================================
      !     SETUP
      !====================================================================

         subroutine build_torsion_coeffs()
      !! Builds the coefficient fields from whatever t2Dh%init_geom registered.
      !! Idempotent: does nothing after the first successful call, since
      !! lambda/axial centre are fixed for the run. Safe to call every step.
            real(dp), dimension(lv) :: xflat, yflat
            real(dp) :: r, lam, zc, Rc
            integer :: i, ntot1

            if (torsion_coeffs_defined) return

            if_torsion_zero = t2Dh%is_torsion_zero()
            if (if_torsion_zero) then
               torsion_coeffs_defined = .true.
               return
            end if

            ntot1 = lx1*ly1*lz1*nelv
            lam = t2Dh%get_lambda()
            zc  = t2Dh%get_axial_centre()
            Rc  = t2Dh%get_curv_radius()

      ! xm1/ym1 are rank-4 (lx1,ly1,lz1,lelt); flatten once via copy (the
      ! same sequence-association idiom used throughout this codebase, e.g.
      ! compute_frc_div_axisym's copy(wtmp, w(1,1,ipert), ...)) so everything
      ! below is a plain single-index loop.
            call copy(xflat, xm1, ntot1)
            call copy(yflat, ym1, ntot1)

            do i = 1, ntot1
               r = yflat(i)
               zhat_f  (i) = xflat(i) - zc
               rhat_f  (i) = r - Rc
               invR_f  (i) = 1.0_dp / r
               lamR_f  (i) = lam / r
               lamR2_f (i) = lam / r**2
               lam2R2_f(i) = lam*lam / r**2
               lamzR3_f(i) = lam*zhat_f(i) / r**3
            end do

            torsion_coeffs_defined = .true.
         end subroutine build_torsion_coeffs

         subroutine torsion_shift_add(shift_r, shift_z)
      !! Adds the torsion part of diag(Omega^2)/R^2 to the R and Z Helmholtz
      !! shifts:
      !!    S_R   = (1 + lambda^2)/R^2      (torus value 1/R^2, already in shift_r)
      !!    S_Z   =      lambda^2 /R^2      (torus value 0,     already in shift_z)
      !!    S_PHI =  unchanged -- not touched here, diag(Omega^2)_phi = -1
      !!             regardless of lambda.
            real(dp), dimension(lx1,ly1,lz1,lelv), intent(inout) :: shift_r, shift_z
            integer :: ntot1
            if (if_torsion_zero) return
            ntot1 = lx1*ly1*lz1*nelv
            call add2(shift_r, lam2R2_f, ntot1)
            call add2(shift_z, lam2R2_f, ntot1)
         end subroutine torsion_shift_add

      !====================================================================
      !     THE dtheta KERNEL AND ITS EXACT TRANSPOSE
      !====================================================================

         subroutine apply_dtheta(dout, din)
      !! dout = dtheta(din) = zhat*d(din)/dR - rhat*d(din)/dZ. Nek's x is the
      !! axial direction Z, Nek's y the radial direction R (r = ym1), matching
      !! neklab_2Dh_axisym's convention throughout.
      !!
      !! Element-local collocation derivative, exactly as gradm1 returns it --
      !! deliberately not averaged across element boundaries, so that
      !! apply_dtheta_t below is its exact discrete transpose.
            real(dp), dimension(lx1*ly1*lz1*lelv), intent(out) :: dout
            real(dp), dimension(lx1*ly1*lz1*lelv), intent(in) :: din
            real(dp), dimension(lx1*ly1*lz1*lelv) :: gx, gy, gz
            integer :: i, ntot1
            ntot1 = lx1*ly1*lz1*nelv
            call gradm1(gx, gy, gz, din)
            do i = 1, ntot1
               dout(i) = zhat_f(i)*gy(i) - rhat_f(i)*gx(i)
            end do
         end subroutine apply_dtheta

         subroutine apply_dtheta_t(dout, din)
      !! dout = dtheta^T(din), the EXACT discrete transpose of apply_dtheta
      !! under the plain Euclidean inner product (no mass weight). Needed
      !! because the divergence gains -(lambda/R)*dtheta(u_phi), and the
      !! streamwise momentum equation must carry its exact adjoint or the
      !! projection loses the discrete divergence-free property. dtheta is
      !! skew only up to the weight term dtheta(R) = zhat, so "apply dtheta
      !! and flip the sign" is NOT its transpose.
      !!
      !! Runs gradm1 backwards. gradm1 calls
      !! local_grad2(ur,us,u,N,e,dxm1,dytm1); the transpose of THAT specific
      !! pairing needs (D,Dt) = (dym1,dxtm1) in local_grad2_t, not the
      !! (dxm1,dxtm1) pair used elsewhere in Nek. The two coincide only before
      !! setaxdy has been called for a deformed (ifrzer) element.
            real(dp), dimension(lx1*ly1*lz1*lelv), intent(out) :: dout
            real(dp), dimension(lx1*ly1*lz1*lelv), intent(in) :: din
            integer, parameter :: lxyz = lx1*ly1*lz1
            real(dp), dimension(lxyz) :: ur, us, wtmp, ax, ay
            integer :: e, i, i0, N

            if (if3d) call nek_stop_error('apply_dtheta_t is 2D only.',
     &         this_module, 'apply_dtheta_t')

            N = lx1 - 1
            do e = 1, nelv
               i0 = (e-1)*lxyz
               do i = 1, lxyz
                  ax(i) = -rhat_f(i0+i) * din(i0+i)
                  ay(i) =  zhat_f(i0+i) * din(i0+i)
               end do
               do i = 1, lxyz
                  ur(i) = jacmi(i,e)*(rxm1(i,1,1,e)*ax(i) + rym1(i,1,1,e)*ay(i))
                  us(i) = jacmi(i,e)*(sxm1(i,1,1,e)*ax(i) + sym1(i,1,1,e)*ay(i))
               end do
               if (ifaxis) call setaxdy(ifrzer(e))
               call local_grad2_t(dout(i0+1), ur, us, N, 1, dym1, dxtm1, wtmp)
            end do
         end subroutine apply_dtheta_t

      !====================================================================
      !     THE D_phi PAIR  (continuity/pressure coupling)
      !====================================================================

         subroutine dphi_apply(dout, uin)
      !! dout = -(lambda/R)*dtheta(uin), on M1, unweighted -- the lambda part
      !! of the streamwise divergence operator D_phi at alpha = 0.
            real(dp), dimension(lx1*ly1*lz1*lelv), intent(out) :: dout
            real(dp), dimension(lx1*ly1*lz1*lelv), intent(in) :: uin
            integer :: ntot1
            ntot1 = lx1*ly1*lz1*nelv
            if (if_torsion_zero) then
               call rzero(dout, ntot1)
               return
            end if
            call apply_dtheta(dout, uin)
            call col2(dout, lamR_f, ntot1)
            call chsign(dout, ntot1)
         end subroutine dphi_apply

         subroutine dphi_apply_t(dout, win)
      !! dout = [-(lambda/R)*dtheta]^T(win) = -dtheta^T[(lambda/R)*win], the
      !! EXACT discrete transpose of dphi_apply. The diagonal scaling is
      !! applied BEFORE dtheta^T -- the reverse order of dphi_apply -- which
      !! is what makes this the transpose of that operator rather than of a
      !! different one that merely looks similar.
            real(dp), dimension(lx1*ly1*lz1*lelv), intent(out) :: dout
            real(dp), dimension(lx1*ly1*lz1*lelv), intent(in) :: win
            real(dp), dimension(lx1*ly1*lz1*lelv) :: tmp
            integer :: ntot1
            ntot1 = lx1*ly1*lz1*nelv
            if (if_torsion_zero) then
               call rzero(dout, ntot1)
               return
            end if
            call col3(tmp, win, lamR_f, ntot1)
            call apply_dtheta_t(dout, tmp)
            call chsign(dout, ntot1)
         end subroutine dphi_apply_t

      !====================================================================
      !     EXPLICIT FORCING  (advection + viscous, alpha-free part only)
      !====================================================================
      !
      !  Both routines below accumulate into LOCAL buffers first and only
      !  ADD the (weighted) result onto the caller's rhs at the very end --
      !  exactly mirroring add_torus_curv_R/phi. rhs arrives already carrying
      !  other, already-weighted contributions (advZ/advR and the torus
      !  curvature terms); weighting the torsion contribution in isolation
      !  before the final add2 is what keeps those earlier contributions from
      !  being weighted a second time.

         subroutine add_torsion_forcing_RZ(rhs_z, rhs_r, jp_)
      !! Adds the torsion RHS to the Z- and R-momentum equations. Injected at
      !! the same point in makefp_2Dh_axisym as add_torus_curv_R, i.e. after
      !! the convective step and before EXT3, so it rides makextp automatically.
      !!
      !! ADVECTION -- from (U.grad)u + (u.grad)U with the torsion-corrected
      !! streamwise derivative and frame rotation, real (alpha-free) part only:
      !!
      !!   Z += +(lambda U_phi/R) dtheta(u_Z) + (lambda U_phi/R) u_R
      !!      +  (lambda u_phi/R) dtheta(U_Z) + (lambda u_phi/R) U_R
      !!   R += +(lambda U_phi/R) dtheta(u_R) - (lambda U_phi/R) u_Z
      !!      +  (lambda u_phi/R) dtheta(U_R) - (lambda u_phi/R) U_Z
      !!
      !! VISCOUS -- the alpha-free groups of the vector Laplacian's torsion
      !! correction (see module header):
      !!
      !!   Z += lambda^2*Ltheta(u_Z) + (2*lambda^2/R^2)*dtheta(u_R)
      !!      + (lambda/R^2)*u_phi - (lambda^2*zhat/R^3)*u_R
      !!   R += lambda^2*Ltheta(u_R) - (2*lambda^2/R^2)*dtheta(u_Z)
      !!      + (2*lambda/R^2)*dtheta(u_phi) + (lambda^2*zhat/R^3)*u_Z
      !!      - (lambda*zhat/R^3)*u_phi
      !!
      !! WEIGHTING. The advection group multiplies rho (vtrans) and the viscous
      !! group multiplies nu (vdiff); the two are NOT interchangeable
      !! (nu = 1/Re, so a viscous term wrongly scaled by rho = 1 would be a
      !! factor Re too large). They are therefore accumulated into separate
      !! buffers and weighted independently before being summed onto the rhs.
      !! This mirrors how the code already treats the two contributions
      !! elsewhere: the advection RHS is weighted by vtrans in
      !! add_torus_curv_R, while the viscous shift terms are weighted by vdiff
      !! inside helmholtz_matvec_2Dh_axisym.
      !!
      !! STABILITY NOTE: the advective group above is a solid-body rotation
      !! of the cross-section at angular rate lambda*U_phi/R. Treated
      !! explicitly, it can become the limiting CFL constraint at large
      !! lambda; if it does, fold it into the advecting velocity passed to
      !! advabp/makeufp instead of adding it here. Not done in this version.
            real(dp), dimension(lx1*ly1*lz1*lelv), intent(inout) :: rhs_z, rhs_r
            integer, intent(in) :: jp_
            real(dp), dimension(lx1*ly1*lz1*lelv) :: aZadv, aRadv, aZvis, aRvis
            real(dp), dimension(lx1*ly1*lz1*lelv) :: coef_base, coef_pert, dthR, dthZ, dthP
            real(dp), dimension(lx1*ly1*lz1*lelv) :: LthR, LthZ, tmp
            real(dp) :: lam
            integer :: ntot1

            if (if_torsion_zero) return
            ntot1 = lx1*ly1*lz1*nelv
            lam = t2Dh%get_lambda()
            call rzero(aZadv, ntot1)
            call rzero(aRadv, ntot1)
            call rzero(aZvis, ntot1)
            call rzero(aRvis, ntot1)

      ! ---- advection (weighted by rho) ----
            call col3(coef_base, lamR_f, t(1,1,1,1,1), ntot1)   ! lambda*U_phi/R
            call col3(coef_pert, lamR_f, tp(1,1,jp_),  ntot1)   ! lambda*u_phi/R

            call apply_dtheta(dthZ, vxp(1,jp_))
            call apply_dtheta(dthR, vyp(1,jp_))
            call apply_dtheta(dthP, tp(1,1,jp_))

            call add_cxy(aZadv, coef_base, dthZ,        1.0_dp, ntot1)
            call add_cxy(aZadv, coef_base, vyp(1,jp_),  1.0_dp, ntot1)
            call add_cxy(aRadv, coef_base, dthR,        1.0_dp, ntot1)
            call add_cxy(aRadv, coef_base, vxp(1,jp_), -1.0_dp, ntot1)

            call apply_dtheta(tmp, vx)
            call add_cxy(aZadv, coef_pert, tmp, 1.0_dp, ntot1)
            call add_cxy(aZadv, coef_pert, vy,  1.0_dp, ntot1)
            call apply_dtheta(tmp, vy)
            call add_cxy(aRadv, coef_pert, tmp,  1.0_dp, ntot1)
            call add_cxy(aRadv, coef_pert, vx,  -1.0_dp, ntot1)

      ! ---- viscous (weighted by nu) ----
            call apply_Ltheta(LthZ, vxp(1,jp_))
            call apply_Ltheta(LthR, vyp(1,jp_))
            call add2s2(aZvis, LthZ, lam*lam, ntot1)
            call add2s2(aRvis, LthR, lam*lam, ntot1)

            call add_cxy(aZvis, lam2R2_f, dthR,         2.0_dp, ntot1)
            call add_cxy(aRvis, lam2R2_f, dthZ,        -2.0_dp, ntot1)
            call add_cxy(aRvis, lamR2_f,  dthP,         2.0_dp, ntot1)

            call add_cxy(aZvis, lamR2_f,  tp(1,1,jp_),  1.0_dp, ntot1)
            call add_cxy(aZvis, lamzR3_f, vyp(1,jp_),  -lam,    ntot1)
            call add_cxy(aRvis, lamzR3_f, vxp(1,jp_),   lam,    ntot1)
            call add_cxy(aRvis, lamzR3_f, tp(1,1,jp_), -1.0_dp, ntot1)

            call weight_forcing(aZadv, vtrans(1,1,1,1,1))
            call weight_forcing(aRadv, vtrans(1,1,1,1,1))
            call weight_forcing(aZvis, vdiff(1,1,1,1,1))
            call weight_forcing(aRvis, vdiff(1,1,1,1,1))
            call add2(rhs_z, aZadv, ntot1)
            call add2(rhs_z, aZvis, ntot1)
            call add2(rhs_r, aRadv, ntot1)
            call add2(rhs_r, aRvis, ntot1)
         end subroutine add_torsion_forcing_RZ

         subroutine add_torsion_forcing_phi(rhs_phi, jp_)
      !! Adds the torsion RHS to the phi ("swirl") equation. Injected in
      !! makefp_2Dh_axisym at the same point as add_torus_curv_phi.
      !!
      !!   phi += (lambda U_phi/R) dtheta(u_phi) + (lambda u_phi/R) dtheta(U_phi)
      !!        + lambda^2*Ltheta(u_phi) - (2*lambda/R^2)*dtheta(u_R)
      !!        + (lambda/R^2)*u_Z + (lambda*zhat/R^3)*u_R
      !!
      !! As in add_torsion_forcing_RZ, the advection line (first two terms)
      !! multiplies rho and the viscous line (last four) multiplies nu, so the
      !! two are weighted separately. The phi ("swirl") field lives in the
      !! temperature slot, so its transport property is vtrans(...,2); its
      !! diffusivity is still vdiff(...,1), the momentum viscosity, because
      !! u_phi is a velocity component carried in the scalar slot, not an
      !! independent scalar with its own conductivity.
            real(dp), dimension(lx1*ly1*lz1*lelv), intent(inout) :: rhs_phi
            integer, intent(in) :: jp_
            real(dp), dimension(lx1*ly1*lz1*lelv) :: aPadv, aPvis
            real(dp), dimension(lx1*ly1*lz1*lelv) :: coef_base, coef_pert, dthR, dthP, LthP, tmp
            real(dp) :: lam
            integer :: ntot1

            if (if_torsion_zero) return
            ntot1 = lx1*ly1*lz1*nelv
            lam = t2Dh%get_lambda()
            call rzero(aPadv, ntot1)
            call rzero(aPvis, ntot1)

            call col3(coef_base, lamR_f, t(1,1,1,1,1), ntot1)
            call col3(coef_pert, lamR_f, tp(1,1,jp_),  ntot1)

            call apply_dtheta(dthP, tp(1,1,jp_))
            call apply_dtheta(dthR, vyp(1,jp_))

      ! ---- advection (weighted by rho) ----
            call add_cxy(aPadv, coef_base, dthP, 1.0_dp, ntot1)
            call apply_dtheta(tmp, t(1,1,1,1,1))
            call add_cxy(aPadv, coef_pert, tmp, 1.0_dp, ntot1)

      ! ---- viscous (weighted by nu) ----
            call apply_Ltheta(LthP, tp(1,1,jp_))
            call add2s2(aPvis, LthP, lam*lam, ntot1)
            call add_cxy(aPvis, lamR2_f,  dthR,       -2.0_dp, ntot1)
            call add_cxy(aPvis, lamR2_f,  vxp(1,jp_),  1.0_dp, ntot1)
            call add_cxy(aPvis, lamzR3_f, vyp(1,jp_),  1.0_dp, ntot1)

            call weight_forcing(aPadv, vtrans(1,1,1,1,2))
            call weight_forcing(aPvis, vdiff(1,1,1,1,1))
            call add2(rhs_phi, aPadv, ntot1)
            call add2(rhs_phi, aPvis, ntot1)
         end subroutine add_torsion_forcing_phi

         subroutine build_torsion_forcing_nl()
      !! Fills the jp = 0 (nonlinear/DNS) neklab forcing vector with the helix
      !! torsion terms, evaluated on the CURRENT base flow (vx, vy, t). Call
      !! once per timestep from userchk; userf/userq then read it back
      !! pointwise via neklab_forcing.
      !!
      !! This is the nonlinear counterpart of add_torsion_forcing_RZ/_phi. Two
      !! differences from those linear routines, both essential:
      !!
      !!   (1) The advection is the genuine quadratic (U.grad)U, not the
      !!       linearised (U.grad)u + (u.grad)U. Concretely, every u_i and its
      !!       dtheta is replaced by the base field U_i: the two contributions
      !!       that were distinct in the linear routine coincide here, so each
      !!       advective group is written ONCE, not twice.
      !!
      !!   (2) The weighting is that of a userf/userq BODY FORCE, not of a
      !!       pre-integrated bfxp/bqp RHS. Nek multiplies ffx/ffy/qvol by the
      !!       mass matrix and (for viscous-looking pieces, none here) by no
      !!       diffusivity, in makef/makeq. So we emit a BARE force DENSITY:
      !!       no bm1, and viscous groups carry an EXPLICIT vdiff factor
      !!       (= nu), since a body force is not diffusivity-scaled by Nek the
      !!       way the Helmholtz operator is. Advective groups carry vtrans
      !!       (= rho) explicitly, matching the existing torus userf force
      !!       ffy = temp**2/y (bare density, rho = 1). Signs follow the same
      !!       LHS -> RHS flip as add_torus_curv_R.
      !!
      !! The phi ("swirl") force is carried in u_phi's slot, i.e. temperature;
      !! it is written into the fz component of the jp = 0 vector, to be read
      !! in userq as qvol = qvol + <fz at this point>.
            real(dp), dimension(lx1*ly1*lz1*lelv) :: fZ, fR, fPhi
            real(dp), dimension(lx1*ly1*lz1*lelv) :: coefU, dthR, dthZ, dthP
            real(dp), dimension(lx1*ly1*lz1*lelv) :: LthR, LthZ, LthP, tmp
            real(dp) :: lam
            integer :: ntot1

            if (if_torsion_zero) then
               if (nid == 0) write(*,*) 'build_torsion_forcing_nl: torsion forcing skipped for jp=0'
               return
            end if
            call build_torsion_coeffs()
            ntot1 = lx1*ly1*lz1*nelv
            lam = t2Dh%get_lambda()
            call rzero(fZ, ntot1)
            call rzero(fR, ntot1)
            call rzero(fPhi, ntot1)

      ! ============ ADVECTION  (nonlinear (U.grad)U, weighted by rho) ============
      ! coef = lambda*U_phi/R
            call col3(coefU, lamR_f, t(1,1,1,1,1), ntot1)
            call apply_dtheta(dthZ, vx)
            call apply_dtheta(dthR, vy)
            call apply_dtheta(dthP, t(1,1,1,1,1))
      ! LHS: (lambda U_phi/R)(-dtheta + Omega')U ; RHS = -LHS
            call add_cxy(fZ,  coefU, dthZ,       -1.0_dp, ntot1)
            call add_cxy(fZ,  coefU, vy,         -1.0_dp, ntot1)  ! -(lambda U_phi/R) U_R
            call add_cxy(fR,  coefU, dthR,       -1.0_dp, ntot1)
            call add_cxy(fR,  coefU, vx,          1.0_dp, ntot1)  ! +(lambda U_phi/R) U_Z
            call add_cxy(fPhi,coefU, dthP,       -1.0_dp, ntot1)
            call col2(fZ,  vtrans(1,1,1,1,1), ntot1)
            call col2(fR,  vtrans(1,1,1,1,1), ntot1)
            call col2(fPhi,vtrans(1,1,1,1,2), ntot1)

      ! ============ VISCOUS  (alpha-free vector-Laplacian torsion, weighted by nu) ============
      ! accumulate the viscous part in tmp-style buffers, then add with vdiff.
      ! reuse dth* on the BASE flow (already computed above).
            call rzero(tmp, ntot1)
      ! (a) lambda^2 * Ltheta(U_i)
            call apply_Ltheta(LthZ, vx)
            call apply_Ltheta(LthR, vy)
            call apply_Ltheta(LthP, t(1,1,1,1,1))
      ! --- build viscous R/Z/phi into separate accumulators via a local helper block
      ! Z viscous:
            call add2s2_nu(fZ, LthZ, lam*lam, ntot1)
            call add_cxy_nu(fZ, lam2R2_f, dthR,        2.0_dp, ntot1)
            call add_cxy_nu(fZ, lamR2_f,  t(1,1,1,1,1),1.0_dp, ntot1)   ! (lambda/R^2) U_phi
            call add_cxy_nu(fZ, lamzR3_f, vy,         -lam,    ntot1)   ! -(lambda^2 zhat/R^3) U_R
      ! R viscous:
            call add2s2_nu(fR, LthR, lam*lam, ntot1)
            call add_cxy_nu(fR, lam2R2_f, dthZ,       -2.0_dp, ntot1)
            call add_cxy_nu(fR, lamR2_f,  dthP,        2.0_dp, ntot1)
            call add_cxy_nu(fR, lamzR3_f, vx,          lam,    ntot1)   ! +(lambda^2 zhat/R^3) U_Z
            call add_cxy_nu(fR, lamzR3_f, t(1,1,1,1,1),-1.0_dp, ntot1)  ! -(lambda zhat/R^3) U_phi
      ! phi viscous:
            call add2s2_nu(fPhi, LthP, lam*lam, ntot1)
            call add_cxy_nu(fPhi, lamR2_f,  dthR,     -2.0_dp, ntot1)
            call add_cxy_nu(fPhi, lamR2_f,  vx,        1.0_dp, ntot1)   ! (lambda/R^2) U_Z
            call add_cxy_nu(fPhi, lamzR3_f, vy,        1.0_dp, ntot1)   ! (lambda zhat/R^3) U_R

            call set_neklab_forcing(fZ, fR, fPhi, 0)
         end subroutine build_torsion_forcing_nl

      !====================================================================
      !     HELPERS
      !====================================================================

         subroutine apply_Ltheta(dout, din)
      !! Ltheta(u) = (1/R) dtheta[ (1/R) dtheta u ]. The intermediate field is
      !! direct-stiffness-averaged before the second derivative: gradm1
      !! returns an element-local quantity, and differentiating a
      !! discontinuous field a second time is meaningless.
            real(dp), dimension(lx1*ly1*lz1*lelv), intent(out) :: dout
            real(dp), dimension(lx1*ly1*lz1*lelv), intent(in) :: din
            real(dp), dimension(lx1*ly1*lz1*lelv) :: g
            integer :: ntot1
            ntot1 = lx1*ly1*lz1*nelv
            call apply_dtheta(g, din)
            call col2(g, invR_f, ntot1)
            call col2(g, bm1, ntot1)
            call dssum(g, lx1, ly1, lz1)
            call col2(g, binvm1, ntot1)
            call apply_dtheta(dout, g)
            call col2(dout, invR_f, ntot1)
         end subroutine apply_Ltheta

         subroutine add_cxy(dst, a, b, s, ntot1)
      !! dst += s*a*b, elementwise.
            real(dp), dimension(lx1*ly1*lz1*lelv), intent(inout) :: dst
            real(dp), dimension(lx1*ly1*lz1*lelv), intent(in) :: a, b
            real(dp), intent(in) :: s
            integer, intent(in) :: ntot1
            real(dp), dimension(lx1*ly1*lz1*lelv) :: tmp
            call col3(tmp, a, b, ntot1)
            call add2s2(dst, tmp, s, ntot1)
         end subroutine add_cxy

         subroutine add_cxy_nu(dst, a, b, s, ntot1)
      !! dst += s * nu * a * b, elementwise. The viscous body-force weighting:
      !! a body force is not diffusivity-scaled by Nek, so nu is applied here.
            real(dp), dimension(lx1*ly1*lz1*lelv), intent(inout) :: dst
            real(dp), dimension(lx1*ly1*lz1*lelv), intent(in) :: a, b
            real(dp), intent(in) :: s
            integer, intent(in) :: ntot1
            real(dp), dimension(lx1*ly1*lz1*lelv) :: tmp
            call col3(tmp, a, b, ntot1)
            call col2(tmp, vdiff(1,1,1,1,1), ntot1)
            call add2s2(dst, tmp, s, ntot1)
         end subroutine add_cxy_nu

         subroutine add2s2_nu(dst, src, s, ntot1)
      !! dst += s * nu * src, elementwise (scalar-coefficient viscous term).
            real(dp), dimension(lx1*ly1*lz1*lelv), intent(inout) :: dst
            real(dp), dimension(lx1*ly1*lz1*lelv), intent(in) :: src
            real(dp), intent(in) :: s
            integer, intent(in) :: ntot1
            real(dp), dimension(lx1*ly1*lz1*lelv) :: tmp
            call col3(tmp, src, vdiff(1,1,1,1,1), ntot1)
            call add2s2(dst, tmp, s, ntot1)
         end subroutine add2s2_nu

         subroutine weight_forcing(rhs, prop)
      !! Strong-form residual -> the weighting bfxp/bfyp/bqp expect: multiply
      !! by the supplied transport property (rho = vtrans for advection,
      !! nu = vdiff for viscous) and by bm1. The caller chooses which property
      !! so that advection and viscous groups, which differ by a factor of Re,
      !! are never scaled by the same field.
            real(dp), dimension(lx1*ly1*lz1*lelv), intent(inout) :: rhs
            real(dp), dimension(lx1*ly1*lz1*lelv), intent(in) :: prop
            integer :: ntot1
            ntot1 = lx1*ly1*lz1*nelv
            call col2(rhs, prop, ntot1)
            call col2(rhs, bm1, ntot1)
         end subroutine weight_forcing

      end module neklab_2Dh_torsion