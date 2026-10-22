      module neklab_2Dh_axisym
         use LightKrylov, only: dp
         use neklab_2Dh_common, only: wmask, wmask_defined, build_wmask,
     &                                bm2_map12_2Dh, bm2_map12T_2Dh,
     &                                solve_helmholtz_2Dh, solve_pressure_2Dh,
     &                                setrhs_pressure_2Dh, gensoln_pressure_2Dh
         use neklab_2Dh_torsion, only: if_torsion_zero,
     &                                 build_torsion_coeffs, torsion_shift_add,
     &                                 dphi_apply, dphi_apply_t,
     &                                 dphi_alpha_apply, dphi_alpha_apply_t,
     &                                 add_torsion_forcing_RZ, add_torsion_forcing_phi,
     &                                 add_torsion_coupling
         use neklab_nek_setup, only: nek_log_debug, nek_log_message, nek_stop_error
         use neklab_timing, only: neklab_timer_start, neklab_timer_stop, neklab_clock_tic, neklab_clock_toc,
     &                            t_axisym_advance, t_axisym_adv_setup, t_axisym_adv_advection, t_axisym_adv_rhs,
     &                            t_axisym_adv_uz, t_axisym_adv_urphi, t_axisym_adv_pressure, t_axisym_adv_correction,
     &                            t_axisym_makefp, t_axisym_solve_helm, t_axisym_solve_cpl, t_axisym_solve_pres,
     &                            t_axisym_pres_precond, t_axisym_pres_proj, t_axisym_helm_matvec, t_axisym_cpl_matvec,
     &                            t_axisym_pres_matvec
         implicit none
         include "SIZE"
         include "TOTAL"
         
         private
         character(len=*), parameter, private :: this_module = 'neklab_2Dh_axisym'
         integer, parameter :: maxit_p = 1000 ! max iterations for pressure solve
         
         real(dp), private :: pbasis(lx2*ly2*lz2*lelv,mxprev,lpert)
         real(dp), private :: dtbdlast = 0.0_dp
         real(dp), private :: betaz_op = 0.0_dp ! module variable to pass to the abstract operators
         real(dp), dimension(:,:,:,:), pointer, private :: shift_op => null() ! shift seen by helm_matvec_axisym
         integer,  private :: nprev(lpert) = 0
         
         real(dp), dimension(lx1,ly1,lz1,lelv), public, target :: h2z_shift
         real(dp), dimension(lx1,ly1,lz1,lelv), public, target :: diag_shift
         real(dp), dimension(lx1,ly1,lz1,lelv), public :: couple_coef
      ! diag_shift is S_PHI = (alpha^2+1)/R^2. Torsion splits the R equation's
      ! shift away from it: S_R = (alpha^2+1+lambda^2)/R^2. In the torus the
      ! two coincide and one field sufficed; diag_shift_r carries the R value
      ! whenever lambda /= 0, and is byte-identical to diag_shift otherwise.
         real(dp), dimension(lx1,ly1,lz1,lelv), public, target :: diag_shift_r
         real(dp), dimension(lx1,ly1,lz1,lelv), public :: alphaR_coef   ! alpha/R, continuity/pressure coupling
         logical, public :: torus_coeffs_defined = .false.
         logical, public :: alpha_eq_zero = .false.
         real(dp), private :: alpha_cached = huge(1.0_dp)
         
         public :: build_torus_coeffs
         public :: nek_advance_2Dh_axisym
         public :: pressure_matvec_coupled_2Dh_axisym
         public :: solve_pressure_coupled_2Dh_axisym
         
      contains
         
         subroutine nek_advance_2Dh_axisym(alpha)
            implicit none
            real(dp), intent(in) :: alpha
            
            real(dp), dimension(lx1,ly1,lz1,lelv) :: w1, w2, w3
            real(dp), dimension(lx1,ly1,lz1,lelv) :: resv1, resv2, resv3
            real(dp), dimension(lx1,ly1,lz1,lelv) :: dv1, dv2, dv3
            real(dp), dimension(lx1,ly1,lz1,lelv) :: h1, h2, h2inv
            common /SCRNS/  w1, w2, w3, resv1, resv2, resv3, dv1, dv2, dv3
            common /SCRVH/  h1
            common /SCRCH/  h2
            common /scrhi/  h2inv
            
            real(dp), dimension(lx1,ly1,lz1,lelv,lpert) :: resv2_jp, resv3_jp
            real(dp), dimension(lx1*ly1*lz1*lelv,lpert) :: gradp, advZ, advR, advPhi
            
            real(dp), dimension(lx2*ly2*lz2*lelv,lpert) :: dpr
            real(dp), dimension(lx2*ly2*lz2*lelv) :: prextr, frc_div, onep, ep
            real(dp), dimension(lx1*ly1*lz1*lelv) :: dw
            real(dp) :: ebar, csign_jp
            
            character(len=*), parameter :: this_procedure = 'nek_advance_2Dh_axisym'
            character(len=256) :: msg
            character(len=2) :: info_str
            character(len=7) :: helm_info
            real(dp) :: dtbd
            logical :: ifprjp
            integer :: igeom, iter, intype, ipert
            integer :: ntot1, ntot2, istart
            real, external :: glsum
            
            ntot1 = lx1*ly1*lz1*nelv
            ntot2 = lx2*ly2*lz2*nelv
            
            if (.not. wmask_defined) call build_wmask(.true.)
            if (.not. torus_coeffs_defined .or. alpha /= alpha_cached) call build_torus_coeffs(alpha)
            
            ! Re/Im pairing is only meaningful for alpha /= 0; npert must match
            ! or the caller has set up the wrong number of perturbation slots.
            alpha_eq_zero = (alpha == 0.0_dp)
            if (alpha_eq_zero) then
               if (lpert /= 1) call nek_stop_error('alpha = 0 requires npert = 1.', this_module, this_procedure)
            else
               if (lpert /= 2) call nek_stop_error('alpha /= 0 requires npert = 2.', this_module, this_procedure)
            end if
            if (istep == 1) nprev(:) = 0
            
            ! Timing. 
            call neklab_timer_start(t_axisym_advance)
            call neklab_clock_tic()
            
            call neklab_timer_start(t_axisym_adv_setup)
            call nekgsync
            call settime
            call setup_convect(2)
            call setsolv
            call comment
            call setprop
            
            ifield = 1
            imesh  = 1
            call unorm
            call settolv
            call neklab_timer_stop(t_axisym_adv_setup)
            
            ! --- pressure-gradient / alpha-R coupling forcing (analog of gradz_p)
            
            msg = 'Compute explicit pressure gradient term for w'
            call nek_log_debug(msg,this_module,this_procedure)
            
            call neklab_timer_start(t_axisym_adv_advection)
            do jp = 1, npert
               call extrapprp(prextr)
               call compute_gradp_axisym(gradp, prextr, alphaR_coef)
               call compute_torus_s_advection(advZ, advR, advPhi, jp)
            end do
            ! Torsion alpha*lambda coupling. Like compute_torus_s_advection it reads
            ! slot jp and writes the OPPOSITE slot, so it runs in this same pre-loop,
            ! augmenting (not overwriting) the advZ/advR/advPhi swap-slot buffers,
            ! which makefp_2Dh_axisym then adds into bfxp/bfyp/bqp per jp. Guarded
            ! internally on alpha /= 0 and lambda /= 0.
            if (.not. alpha_eq_zero .and. .not. if_torsion_zero) then
               do jp = 1, npert
                  call add_torsion_coupling(advZ, advR, advPhi, alpha, jp)
               end do
            end if
            call neklab_timer_stop(t_axisym_adv_advection)
            
            ! --- assemble RHS for both jp; solve u_Z immediately (uncoupled across jp);
            !     stash u_R/u_phi RHS for the deferred, jointly-solved transform step.
            
            msg = 'Solve momentum equations for u/v/w'
            call nek_log_debug(msg, this_module, this_procedure)
            
            do jp = 1, npert
               info_str = merge('Re', 'Im', jp == 1)
               call neklab_timer_start(t_axisym_adv_rhs)
               igeom = 1
               !
               !  Old geometry, old velocity
               !
               call makefp_2Dh_axisym(advZ, advR, advPhi)
               !
               igeom = 2
               !
               !  New geometry, new velocity
               !
               intype = -1
               ifield = 1
               call sethlm(h1, h2, intype)
               ! cresvipp
               call bcdirvc(vxp(1,jp), vyp(1,jp), vzp(1,jp), v1mask, v2mask, v3mask)
               call col2   (tp(1,1,jp), wmask, ntot1)
               call extrapprp(prextr)
               call opgradt(resv1, resv2, resv3, prextr)
               call opadd2 (resv1, resv2, resv3, bfxp(1,jp), bfyp(1,jp), bfzp(1,jp))
               
               ! compute w term
               call copy(resv3, gradp(1,jp),  ntot1)
               call add2(resv3, bqp(1,1,jp),  ntot1)
               
               ! hv1
               call helmholtz_matvec_2Dh_axisym(w1, vxp(1,jp), h1, h2, h2z_shift, 1)
               call sub2 (resv1, w1, ntot1)
               
               call dssum(resv1, lx1, ly1, lz1)
               call col2 (resv1, v1mask, ntot1)
               helm_info = info_str//' V_Z  '
               if (istep < 10) call chktcg1(tolhv, resv1, h1, h2, v1mask, vmult, imesh, 1)
               call neklab_timer_stop(t_axisym_adv_rhs)
               
               call neklab_timer_start(t_axisym_adv_uz)
               shift_op => h2z_shift
               call solve_helmholtz_2Dh(dv1, resv1, h1, h2, v1mask, vmult, imesh,
     &                                  tolhv, nmxv, 1, binvm1, helm_info, helm_matvec_axisym)
               call add2(vxp(1,jp), dv1, ntot1)
               call neklab_timer_stop(t_axisym_adv_uz)
               
               ! stash u_R/u_phi RHS for the joint transform+solve below
               call neklab_timer_start(t_axisym_adv_rhs)
               call copy(resv2_jp(1,1,1,1,jp), resv2, ntot1)
               call copy(resv3_jp(1,1,1,1,jp), resv3, ntot1)
               call neklab_timer_stop(t_axisym_adv_rhs)
            end do ! jp
            
            call neklab_timer_start(t_axisym_adv_urphi)
            if (alpha_eq_zero) then
               do jp = 1, npert
                  call helmholtz_matvec_2Dh_axisym(w2, vyp(1,jp), h1, h2, diag_shift_r, 1)
                  call sub2(resv2_jp(1,1,1,1,jp), w2, ntot1)
                  call helmholtz_matvec_2Dh_axisym(w3, tp(1,1,jp), h1, h2, diag_shift, 1)
                  call sub2(resv3_jp(1,1,1,1,jp), w3, ntot1)
               end do
            else
               do jp = 1, npert
                  ipert = npert + 1 - jp
                  csign_jp = merge(1.0_dp, -1.0_dp, jp == 1)
                  call coupled_helmholtz_matvec_2Dh_axisym(w2, w3, vyp(1,jp), tp(1,1,ipert), h1, h2, diag_shift_r, diag_shift, couple_coef, csign_jp)
                  call sub2(resv2_jp(1,1,1,1,jp),    w2, ntot1)
                  call sub2(resv3_jp(1,1,1,1,ipert), w3, ntot1)
               end do
            end if
            
            msg = 'Solve u_R/u_phi momentum equations'
            call nek_log_debug(msg, this_module, this_procedure)
            
            if (alpha_eq_zero) then
               do jp = 1, npert
                  info_str = merge('Re', 'Im', jp == 1)
                  call dssum(resv2_jp(1,1,1,1,jp), lx1, ly1, lz1)
                  call col2 (resv2_jp(1,1,1,1,jp), v2mask, ntot1)
                  helm_info = info_str//' V_R  '
                  if (istep < 10) call chktcg1(tolhv, resv2_jp(1,1,1,1,jp), h1, h2, v2mask, vmult, imesh, 1)
                  shift_op => diag_shift_r
                  call solve_helmholtz_2Dh(dv2, resv2_jp(1,1,1,1,jp), h1, h2, v2mask, vmult, imesh,
     &                                     tolhv, nmxv, 1, binvm1, helm_info, helm_matvec_axisym)
                  
                  call dssum(resv3_jp(1,1,1,1,jp), lx1, ly1, lz1)
                  call col2 (resv3_jp(1,1,1,1,jp), wmask, ntot1)
                  helm_info = info_str//' V_PHI'
                  if (istep < 10) call chktcg1(tolhv, resv3_jp(1,1,1,1,jp), h1, h2, wmask, vmult, imesh, 1)
                  shift_op => diag_shift
                  call solve_helmholtz_2Dh(dv3, resv3_jp(1,1,1,1,jp), h1, h2, wmask, vmult, imesh,
     &                                     tolhv, nmxv, 1, binvm1, helm_info, helm_matvec_axisym)
                  
                  call add2(vyp(1,jp),  dv2, ntot1)
                  call add2(tp(1,1,jp), dv3, ntot1)
               end do
            else
               do jp = 1, npert
                  ipert = npert + 1 - jp
                  csign_jp = merge(1.0_dp, -1.0_dp, jp == 1)
                  helm_info = merge('BLK A  ', 'BLK B  ', jp == 1)
                  call solve_coupled_helmholtz_2Dh_axisym(dv2, dv3, resv2_jp(1,1,1,1,jp), resv3_jp(1,1,1,1,ipert),
     &                                                    h1, h2, v2mask, wmask, vmult, imesh, tolhv, nmxv, binvm1, 
     &                                                    helm_info, diag_shift_r, diag_shift, couple_coef, csign_jp)
                  
                  call add2(vyp(1,jp),     dv2, ntot1)
                  call add2(tp(1,1,ipert), dv3, ntot1)
               end do
            end if
            call neklab_timer_stop(t_axisym_adv_urphi)
            
            ! --- pressure correction stage: identical structure to the beta case,
            !     field routines substituted for the beta_z ones
            msg = 'Compute pressure correction to enforce mass balance'
            call nek_log_debug(msg, this_module, this_procedure)
            
            call neklab_timer_start(t_axisym_adv_pressure)
            ifield = 1
            intype = 1
            dtbd   = bd(1)/dt
            if (dtbd /= dtbdlast) then
               nprev(:) = 0
               dtbdlast = dtbd
            end if
            
            call rzero  (h1, ntot1)
            call cmult2 (h2, vtrans(1,1,1,1,ifield), dtbd, ntot1)
            call invers2(h2inv, h2, ntot1)
            
            ebar = 0.0_dp
            if (.not. alpha_eq_zero .and. ifvcor) then
               call rone(onep, ntot2)
               call pressure_matvec_2Dh_axisym(ep, onep, h1, h2, h2inv, alphaR_coef, intype)
               ebar = glsum(ep, ntot2)
            end if
            
            ifprjp = .false.
            istart = param(95)
            if (istep >= istart .and. istart /= 0) ifprjp = .true.
      ! Pressure projection is a per-slot optimisation and is not applied in the
      ! coupled (alpha /= 0 & lambda /= 0) path, whose two slots are solved
      ! jointly. Correctness is unaffected; only the initial guess quality.
            if (.not. alpha_eq_zero .and. .not. if_torsion_zero) ifprjp = .false.
            
            if (.not. alpha_eq_zero .and. .not. if_torsion_zero) then
      ! ---- coupled path: build both RHS, then one stacked solve ----
               do jp = 1, npert
                  call opdiv(dpr(1,jp), vxp(1,jp), vyp(1,jp), vzp(1,jp))
                  call compute_frc_div_axisym(frc_div, tp, alphaR_coef)
                  call add2 (dpr(1,jp), frc_div, ntot2)
                  call chsign(dpr(1,jp), ntot2)
               end do
               call solve_pressure_coupled_2Dh_axisym(dpr(1,1), dpr(1,2), h1, h2, h2inv,
     &                                              alphaR_coef, alpha, ebar, intype, iter)
            else
      ! ---- decoupled path: per-slot solve (torus or base-flow) ----
               do jp = 1, npert
                  info_str = merge('Re', 'Im', jp == 1)
                  call opdiv(dpr(1,jp), vxp(1,jp), vyp(1,jp), vzp(1,jp))
                  call compute_frc_div_axisym(frc_div, tp, alphaR_coef)
                  call add2 (dpr(1,jp), frc_div, ntot2)
                  call chsign(dpr(1,jp), ntot2)
                  
                  if (alpha_eq_zero) call ortho(dpr(1,jp))
                  
                  if (ifprjp) then
                     call  setrhs_pressure_2Dh(dpr(1,jp), h1, h2, h2inv, pbasis(1,1,jp), nprev(jp),
     &                                         pres_matvec_axisym, info_str)
                  end if
                  call      solve_pressure_2Dh(dpr(1,jp), h1, h2, h2inv, pres_matvec_axisym, intype, iter,
     &                                         ebar, alpha_eq_zero, maxit_p, 'alpha', alpha, info_str)
                  if (ifprjp) then
                     call gensoln_pressure_2Dh(dpr(1,jp), h1, h2, h2inv, pbasis(1,1,jp), nprev(jp),
     &                                         pres_matvec_axisym)
                  end if
               end do
            end if
            call neklab_timer_stop(t_axisym_adv_pressure)
            
            msg = 'Compute velocity correction based on pressure correction'
            call nek_log_debug(msg, this_module, this_procedure)
            
            call neklab_timer_start(t_axisym_adv_correction)
            do jp = 1, npert
               call opgradt(w1, w2, w3, dpr(1,jp))
               call opbinv (dv1, dv2, dv3, w1, w2, w3, h2inv)
               call compute_dw_axisym(dw, dpr, h2inv, alphaR_coef)
               
               call opadd2(vxp(1,jp), vyp(1,jp), vzp(1,jp), dv1, dv2, dv3)
               call add2  (tp(1,1,jp), dw, ntot1)
               
               call extrapprp(prextr)
               call lagpresp
               call add3(prp(1,jp), prextr, dpr(1,jp), ntot2)
            end do
            call neklab_timer_stop(t_axisym_adv_correction)
            
            call neklab_clock_toc()
            call neklab_timer_stop(t_axisym_advance)
            
            ! reset jp = 0 in case we switch to the nonlinear solver next
            jp = 0
            
         end subroutine nek_advance_2Dh_axisym
         
         subroutine build_torus_coeffs(alpha)
            implicit none
            real(dp), intent(in) :: alpha
            integer :: ix, iy, iz, ie
            real(dp) :: r
            
            do ie = 1, nelv
               do iz = 1, lz1
               do iy = 1, ly1
               do ix = 1, lx1
                  r = ym1(ix,iy,iz,ie)
                  alphaR_coef(ix,iy,iz,ie) = alpha / r
                  h2z_shift  (ix,iy,iz,ie) = alpha**2 / r**2
                  diag_shift (ix,iy,iz,ie) = h2z_shift(ix,iy,iz,ie) + 1.0_dp / r**2
                  diag_shift_r(ix,iy,iz,ie) = diag_shift(ix,iy,iz,ie)
                  couple_coef(ix,iy,iz,ie) = 2.0_dp * alpha / r**2
               end do
               end do
               end do
            end do
            
            ! Torsion: adds lambda^2/R^2 to the R and Z shifts, leaves the phi shift
            ! (diag_shift) alone. Idempotent and cheap after the first call, so it is
            ! safe to call unconditionally here rather than threading a second cache
            ! flag through this routine.
            call build_torsion_coeffs()
            call torsion_shift_add(diag_shift_r, h2z_shift)
            
            alpha_cached          = alpha
            torus_coeffs_defined  = .true.
         end subroutine build_torus_coeffs
         
         subroutine helmholtz_matvec_2Dh_axisym(Au, u, h1, h2, shift, isd)
            implicit none
            ! Same CG loop as solve_helmholtz_2Dh_axisym, but operating on the
            ! stacked 2-field vector (x1;x2). dssum/mask applied to x1 and x2
            real, dimension(lx1,ly1,lz1,1), intent(out) :: Au
            real, dimension(lx1,ly1,lz1,1), intent(in)  :: u, h1, h2, shift
            integer, intent(in) :: isd
            real, dimension(lx1,ly1,lz1,lelv) :: tmp
            integer :: imesh, ntot1
            call neklab_timer_start(t_axisym_helm_matvec)
            imesh = 1
            ntot1 = lx1*ly1*lz1*nelv
            call axhelm(Au, u, h1, h2, imesh, isd)
            call col3 (tmp, u, vdiff(1,1,1,1,1), ntot1)
            call col2 (tmp, shift, ntot1)              ! field multiply, replaces beta_z**2 scalar
            call col2 (tmp, bm1, ntot1)
            call add2 (Au, tmp, ntot1)
            call neklab_timer_stop(t_axisym_helm_matvec)
         end subroutine helmholtz_matvec_2Dh_axisym
         
         subroutine compute_gradp_axisym(gradp, prextr, alphaR)
            implicit none
            include 'SIZE'
            include 'SOLN'
            include 'MASS'
            real(dp), dimension(lx1*ly1*lz1*lelv,lpert), intent(inout) :: gradp
            real(dp), dimension(lx2*ly2*lz2*lelv), intent(in) :: prextr
            real(dp), dimension(lx1,ly1,lz1,lelv), intent(in) :: alphaR
            integer :: ipert, spert, ntot1
            real(dp), dimension(lx1*ly1*lz1) :: wrk1, wrk2
            real(dp), dimension(lx1*ly1*lz1*lelv) :: wl1, wl2
            
            ntot1 = lx1*ly1*lz1*nelv
            ipert = npert + 1 - jp
            
            if (alpha_eq_zero .and. if_torsion_zero) return
            
      ! --- alpha part: skew in slot space (swaps jp <-> ipert). Torus term,
      !     unchanged, absent entirely at alpha = 0.
            if (.not. alpha_eq_zero) then
               spert = merge(1, -1, jp == 1)
               call bm2_map12T_2Dh(gradp(1,ipert), prextr)           ! I12^T (bm2 p)
               call col2 (gradp(1,ipert), alphaR, ntot1)
               if (spert == -1) call chsign(gradp(1,ipert), ntot1)
               ! dssum will be called on the vel residual later
            end if
            
      ! --- lambda part: real, so it stays in slot jp (never ipert). This is
      !     the exact transpose of the lambda block compute_frc_div_axisym
      !     adds to the divergence: dphi_apply_t is deliberately "scale THEN
      !     transpose", the reverse order of dphi_apply's "transpose THEN
      !     scale", which is what makes the pair exact adjoints of each other
      !     rather than merely similar-looking operators.
            if (.not. if_torsion_zero) then
               call mappr(wl1, prextr, wrk1, wrk2)
               call col2 (wl1, bm1, ntot1)
               call dphi_apply_t(wl2, wl1)
               call add2 (gradp(1,jp), wl2, ntot1)
               ! dssum will be called on the vel residual later
            end if
         end subroutine compute_gradp_axisym
         
         subroutine compute_dw_axisym(dw, dpr, h2inv, alphaR)
            implicit none
            include 'SIZE'
            include 'SOLN'
            include 'MASS'
            real(dp), dimension(lx1*ly1*lz1*lelv), intent(out) :: dw
            real(dp), dimension(lx2*ly2*lz2*lelv,lpert), intent(in) :: dpr
            real(dp), dimension(lx1,ly1,lz1,lelv), intent(in) :: h2inv, alphaR
            integer :: ipert, spert, ntot1
            real(dp), dimension(lx1*ly1*lz1) :: wrk1, wrk2
            real(dp), dimension(lx1*ly1*lz1*lelv) :: wl1, wl2
            
            ntot1 = lx1*ly1*lz1*nelv
            call rzero(dw, ntot1)
            
            if (alpha_eq_zero .and. if_torsion_zero) return
            
            if (.not. alpha_eq_zero) then
               ipert = npert + 1 - jp
               spert = merge(1, -1, jp == 1)
               call bm2_map12T_2Dh(dw, dpr(1,ipert))      ! I12^T (bm2 dp)
               call col2 (dw, alphaR, ntot1)
               if (spert == 1) call chsign(dw, ntot1)     ! matches -spert*beta_z sign convention
            end if
            
      ! Lambda part, from THIS slot's own pressure correction. Same operator
      ! as compute_gradp_axisym's lambda block; the overall correction here
      ! must UNDO the gradient contribution, hence sub2 rather than add2.
            if (.not. if_torsion_zero) then
               call mappr(wl1, dpr(1,jp), wrk1, wrk2)
               call col2 (wl1, bm1, ntot1)
               call dphi_apply_t(wl2, wl1)
               call sub2 (dw, wl2, ntot1)
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
            real(dp), dimension(lx1*ly1*lz1*lelv,ldimt,lpert), intent(in) :: w
            real(dp), dimension(lx1,ly1,lz1,lelv), intent(in) :: alphaR
            integer :: ipert, spert, ie, ie1, nxyz1, ntot1, ntot2
            real(dp), dimension(lx1,ly1,lz1,lelv) :: wtmp
            real(dp), dimension(lx1*ly1*lz1*lelv) :: wl1
            
            ntot1 = lx1*ly1*lz1*nelv
            ntot2 = lx2*ly2*lz2*nelv
            
            if (alpha_eq_zero .and. if_torsion_zero) return
            
            ! --- alpha part, from the opposite slot
            if (.not. alpha_eq_zero) then
               ipert = npert + 1 - jp
               spert = merge(1, -1, jp == 1)
               call copy(wtmp, w(1,1,ipert), ntot1)
               call col2(wtmp, alphaR, ntot1)       ! scale on M1 BEFORE mapping to M2
               if (spert == -1) call chsign(wtmp, ntot1)
            end if
            
            ! --- lambda part, from THIS slot: -(lambda/R)*dtheta(u_phi). Accumulated
            !     on M1 alongside the alpha part so a single map12 below carries both.
            if (.not. if_torsion_zero) then
               call dphi_apply(wl1, w(1,1,jp))
               call add2(wtmp, wl1, ntot1)
            end if
            
            call bm2_map12_2Dh(frc_div, wtmp)           
         end subroutine compute_frc_div_axisym
         
         subroutine makefp_2Dh_axisym(advZ, advR, advPhi)
            ! Perturbation RHS assembly for the torus formulation.
            !
            ! Replaces the core `makefp` + `makeqp` pair. Structure is identical
            ! to those routines (perturb.f), with the toroidal explicit terms
            ! injected AFTER the convective step and BEFORE the extrapolation
            ! step, so that makextp/makeabqp apply EXT3 to them and roll them
            ! into the lag arrays (exx1p/exx2p, vgradt1p/vgradt2p) automatically.
            !
            ! Operates on the current global `jp` (the Nek routines called here
            ! all read it), and restores `ifield` on exit.
            implicit none
            include 'SIZE'
            include 'SOLN'
            include 'INPUT'
            include 'TSTEP'
            include 'ADJOINT'
            real(dp), dimension(lx1*ly1*lz1*lelv,lpert), intent(in) :: advZ, advR, advPhi
            character(len=*), parameter :: this_procedure = 'makefp_2Dh_axisym'
            integer :: ifield_bak, ntot1
            
            call neklab_timer_start(t_axisym_makefp)
            ntot1 = lx1*ly1*lz1*nelv
            ifield_bak = ifield
            
            ! The toroidal terms are linearised about the forward base flow and
            ! are not adjoint-aware; refuse rather than silently solve the wrong
            ! system.
            if (ifadj) call nek_stop_error(
     &         'adjoint mode not supported by the torus formulation',
     &         this_module, this_procedure)
            
            ! ---------------- velocity: u_Z, u_R ----------------
            ifield = 1
            call makeufp
            if (ifnav .and. .not. ifchar) call advabp
            call add2(bfxp(1,jp), advZ(1,jp), ntot1)   ! -i*alpha*(U_phi/R)*u_Z
            call add2(bfyp(1,jp), advR(1,jp), ntot1)   ! -i*alpha*(U_phi/R)*u_R
            call add_torus_curv_R(bfyp(1,jp), jp)      ! +2*U_phi*u'_phi/R
            if (.not. if_torsion_zero) call add_torsion_forcing_RZ(bfxp(1,jp), bfyp(1,jp), jp)
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
            call add2(bqp(1,ifield-1,jp), advPhi(1,jp), ntot1)   ! -i*alpha*(U_phi/R)*u_phi
            call add_torus_curv_phi(bqp(1,ifield-1,jp), jp)      ! -(U_phi*u'_R + U_R*u'_phi)/R
            if (.not. if_torsion_zero) call add_torsion_forcing_phi(bqp(1,ifield-1,jp), jp)
            if (iftran) call makeabqp                            ! <-- EXT3
            if ((iftran .and. .not. ifchar) .or.
     &          (iftran .and. .not. ifadvc(ifield) .and. ifchar)) call makebdqp
            call lagscalp
            
            ifield = ifield_bak
            call neklab_timer_stop(t_axisym_makefp)
         end subroutine makefp_2Dh_axisym
         
         subroutine add_torus_curv_R(rhs, jp_)
            implicit none
            include 'SIZE'
            include 'SOLN'   ! vx,vy,t (baseflow); vxp,vyp,tp (perturbation); vtrans
            include 'MASS'   ! bm1
            real(dp), dimension(lx1*ly1*lz1*lelv), intent(inout) :: rhs
            integer, intent(in) :: jp_
            real(dp), dimension(lx1,ly1,lz1,lelv) :: tmp
            integer :: ntot1
            
            ntot1 = lx1*ly1*lz1*nelv
            
            ! --- R-momentum: + 2*U_phi*u'_phi/R * rho * bm1
            !     linearization of the missing -U_phi^2/R term in (U.grad U)_R
            call col3   (tmp, t(1,1,1,1,1), tp(1,1,jp_), ntot1)   ! U_phi * u'_phi
            call invcol2(tmp, ym1, ntot1)                          ! / R
            call cmult  (tmp, 2.0_dp, ntot1)
            call col2   (tmp, vtrans(1,1,1,1,1), ntot1)            ! rho  (density = 1.0 in your .par)
            call col2   (tmp, bm1, ntot1)
            call add2   (rhs, tmp, ntot1)
            
         end subroutine add_torus_curv_R
         
         subroutine add_torus_curv_phi(rhs, jp_)
            implicit none
            include 'SIZE'
            include 'SOLN'   ! vx,vy,t (baseflow); vxp,vyp,tp (perturbation); vtrans
            include 'MASS'   ! bm1
            real(dp), dimension(lx1*ly1*lz1*lelv), intent(inout) :: rhs
            integer, intent(in) :: jp_
            real(dp), dimension(lx1,ly1,lz1,lelv) :: tmp1, tmp2
            integer :: ntot1
            
            ntot1 = lx1*ly1*lz1*nelv
            
            ! --- swirl ("phi") equation: - (U_phi*u'_R + U_R*u'_phi)/R * rhocp * bm1
            !     linearization of the missing +U_R*U_phi/R term in (U.grad U)_phi
            call col3   (tmp1, t(1,1,1,1,1), vyp(1,jp_),  ntot1)    ! U_phi * u'_R
            call col3   (tmp2, vy,           tp(1,1,jp_), ntot1)    ! U_R   * u'_phi  (drop if U_R deemed negligible)
            call add2   (tmp1, tmp2, ntot1)
            call invcol2(tmp1, ym1, ntot1)                           ! / R
            call chsign (tmp1, ntot1)
            call col2   (tmp1, vtrans(1,1,1,1,2), ntot1)             ! rhocp (= 1.0 in your .par)
            call col2   (tmp1, bm1, ntot1)
            call add2   (rhs, tmp1, ntot1)
            
         end subroutine add_torus_curv_phi
         
         subroutine compute_torus_s_advection(advZ, advR, advPhi, jp_)
            ! -i*alpha*(U_phi/R)*u term, present in ALL THREE momentum equations
            ! because the base flow has U_phi != 0 in the ignorable direction
            ! (no analogue in the planar beta_z case, where W=0 by assumption).
            implicit none
            include 'SIZE'
            include 'SOLN'   ! vxp,vyp,tp; t (baseflow); vtrans
            include 'MASS'   ! bm1
            real(dp), dimension(lx1*ly1*lz1*lelv,lpert), intent(inout) :: advZ, advR, advPhi
            integer, intent(in) :: jp_
            real(dp), dimension(lx1,ly1,lz1,lelv) :: coef
            integer :: ipert, ntot1
            
            ntot1 = lx1*ly1*lz1*nelv
            ipert = npert + 1 - jp_
            call rzero(advZ(1,ipert),   ntot1)
            call rzero(advR(1,ipert),   ntot1)
            call rzero(advPhi(1,ipert), ntot1)
            
            if (alpha_eq_zero) return
            
            call col3(coef, alphaR_coef, t(1,1,1,1,1), ntot1)   ! (alpha/R) * U_phi
            if (ipert == 1) call chsign(coef, ntot1)             ! sign keyed to the READING jp, not the caller
            
            call col3 (advZ(1,ipert), coef, vxp(1,jp_), ntot1)
            call col2 (advZ(1,ipert), vtrans(1,1,1,1,1), ntot1)
            call col2 (advZ(1,ipert), bm1, ntot1)
            
            call col3 (advR(1,ipert), coef, vyp(1,jp_), ntot1)
            call col2 (advR(1,ipert), vtrans(1,1,1,1,1), ntot1)
            call col2 (advR(1,ipert), bm1, ntot1)
            
            call col3 (advPhi(1,ipert), coef, tp(1,1,jp_), ntot1)
            call col2 (advPhi(1,ipert), vtrans(1,1,1,1,2), ntot1)
            call col2 (advPhi(1,ipert), bm1, ntot1)
         end subroutine compute_torus_s_advection
         
         subroutine pressure_matvec_2Dh_axisym(ap, wp, h1, h2, h2inv, alphaR, intype)
            implicit none
            real(dp), dimension(lx2,ly2,lz2,lelv), intent(out) :: Ap
            real(dp), dimension(lx2,ly2,lz2,lelv), intent(in)  :: wp
            real(dp), dimension(lx1,ly1,lz1,lelv), intent(in)  :: h1, h2, h2inv, alphaR
            integer, intent(in) :: intype
            real(dp), dimension(lx1*ly1*lz1) :: wrk1, wrk2
            real(dp), dimension(lx1,ly1,lz1,lelv) :: wdivm1
            real(dp), dimension(lx1*ly1*lz1*lelv) :: wl1
            real(dp), dimension(lx2,ly2,lz2,lelv) :: wdivm2
            integer :: ie, ntot1, ntot2
            
            call neklab_timer_start(t_axisym_pres_matvec)
            call cdabdtp(Ap, wp, h1, h2, h2inv, intype)   ! core Nek, already ifaxis-correct globally
            
            if (alpha_eq_zero .and. if_torsion_zero) then
               call neklab_timer_stop(t_axisym_pres_matvec)
               return
            end if
            
            ntot1 = lx1*ly1*lz1*nelv
            ntot2 = lx2*ly2*lz2*nelv
            
            ! --- alpha block: (alpha/R) H^-1 (alpha/R). The two slot swaps of
            !     D_phi and D_phi^T cancel, so this block is diagonal in slot space
            !     and each perturbation's pressure can be solved on its own.
            if (.not. alpha_eq_zero) then
               call bm2_map12T_2Dh(wdivm1, wp)            ! I12^T (bm2 p), local
               call col2 (wdivm1, alphaR, ntot1)          ! alpha/R
               call col2 (wdivm1, wmask, ntot1)
               call dssum(wdivm1, lx1, ly1, lz1)
               call col2 (wdivm1, binvm1, ntot1)
               call col2 (wdivm1, h2inv, ntot1)
               call col2 (wdivm1, alphaR, ntot1)          ! alpha/R
               call chsign(wdivm1, ntot1)
               call bm2_map12_2Dh(wdivm2, wdivm1)         ! bm2 * I12 (.)
               call sub2(ap, wdivm2, ntot2)
            end if
            
      ! --- lambda block: [-(lambda/R)*dtheta] H^-1 [-(lambda/R)*dtheta]^T.
      !     Carries no slot swap either, so at alpha = 0 (the only case
      !     implemented so far) this is simply a second, additive, symmetric
      !     positive semi-definite contribution to the same pressure operator.
            if (.not. if_torsion_zero) then
               call mappr(wl1, wp, wrk1, wrk2)
               call col2 (wl1, bm1, ntot1)
               call dphi_apply_t(wdivm1, wl1)
               call col2 (wdivm1, wmask, ntot1)
               call dssum(wdivm1, lx1, ly1, lz1)
               call col2 (wdivm1, binvm1, ntot1)
               call col2 (wdivm1, h2inv, ntot1)
               call dphi_apply(wl1, wdivm1)
               call chsign(wl1, ntot1)
               do ie = 1, nelv
                  call map12(wdivm2(1,1,1,ie), wl1((ie-1)*lx1*ly1*lz1+1), ie)
               end do
               call col2(wdivm2, bm2, ntot2)
               call sub2(ap, wdivm2, ntot2)
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
            !! routes inputs/outputs between slots and supplies the swap sign, because
            !! the swap couples the two M2 fields and must be tracked at that level.
            !!
            !! The B^-1 sandwich (mappr, bm1, wmask, dssum, binvm1, h2inv, map12, bm2)
            !! is byte-identical to the diagonal blocks in pressure_matvec_2Dh_axisym;
            !! only the left/right operator halves differ. That identity is what makes
            !! the cross block consistent with the diagonal ones and the whole E
            !! symmetric (term3 = term2^T).
            implicit none
            real(dp), dimension(lx2,ly2,lz2,lelv), intent(out) :: out2
            real(dp), dimension(lx2,ly2,lz2,lelv), intent(in)  :: in2
            real(dp), dimension(lx1,ly1,lz1,lelv), intent(in)  :: h2inv, alphaR
            logical, intent(in) :: left_is_alpha
            real(dp), dimension(lx1*ly1*lz1) :: wrk1, wrk2
            real(dp), dimension(lx1,ly1,lz1,lelv) :: wm1
            real(dp), dimension(lx1*ly1*lz1*lelv) :: wf
            real(dp), dimension(lx2,ly2,lz2,lelv) :: wm2
            integer :: ie, ntot1
            
            ntot1 = lx1*ly1*lz1*nelv
            
            call mappr(wm1, in2, wrk1, wrk2)
            call col2 (wm1, bm1, ntot1)
            
            ! right operator half (applied first): L^T if left is alpha, else A
            if (left_is_alpha) then
               call dphi_apply_t(wf, wm1)                 ! L^T = -dtheta^T[(lambda/R) .]
               call copy(wm1, wf, ntot1)
            else
               call dphi_alpha_apply_t(wf, wm1, alphaR)   ! A^T = A = (alpha/R) .
               call copy(wm1, wf, ntot1)
            end if
            
            call col2 (wm1, wmask, ntot1)
            call dssum(wm1, lx1, ly1, lz1)
            call col2 (wm1, binvm1, ntot1)
            call col2 (wm1, h2inv, ntot1)
            
            ! left operator half (applied last): A if left is alpha, else L
            if (left_is_alpha) then
               call dphi_alpha_apply(wf, wm1, alphaR)      ! A = (alpha/R) .
            else
               call dphi_apply(wf, wm1)                    ! L = -(lambda/R) dtheta .
            end if
            call chsign(wf, ntot1)
            
            do ie = 1, nelv
               call map12(wm2(1,1,1,ie), wf((ie-1)*lx1*ly1*lz1+1), ie)
            end do
            call col2(wm2, bm2, lx2*ly2*lz2*nelv)
            call copy(out2, wm2, lx2*ly2*lz2*nelv)
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
            ! shift1 == shift2 at lambda = 0, recovering the torus block exactly.
            real(dp), dimension(lx1,ly1,lz1,lelv), intent(out) :: y1, y2
            real(dp), dimension(lx1,ly1,lz1,lelv), intent(in)  :: x1, x2, h1, h2, shift1, shift2, couple_coef
            real(dp), intent(in) :: csign
            real(dp), dimension(lx1,ly1,lz1,lelv) :: tmp
            integer :: ntot1
            ntot1 = lx1*ly1*lz1*nelv
            
            call neklab_timer_start(t_axisym_cpl_matvec)
            call helmholtz_matvec_2Dh_axisym(y1, x1, h1, h2, shift1, 1)
            call helmholtz_matvec_2Dh_axisym(y2, x2, h1, h2, shift2, 1)
            
            call col3   (tmp, x2, couple_coef, ntot1)
            call col2   (tmp, vdiff(1,1,1,1,1), ntot1)
            call col2   (tmp, bm1, ntot1)
            call add2s2 (y1, tmp, csign, ntot1)
            
            call col3   (tmp, x1, couple_coef, ntot1)
            call col2   (tmp, vdiff(1,1,1,1,1), ntot1)
            call col2   (tmp, bm1, ntot1)
            call add2s2 (y2, tmp, csign, ntot1)
            call neklab_timer_stop(t_axisym_cpl_matvec)
         end subroutine coupled_helmholtz_matvec_2Dh_axisym

         subroutine helm_matvec_axisym(Au, u, h1, h2)
            ! Wrapper around helmholtz_matvec_2Dh_axisym matching the abstract interface.
            ! The shift field is selected through shift_op.
            real(dp), dimension(lx1*ly1*lz1*lelv), intent(out) :: Au
            real(dp), dimension(lx1*ly1*lz1*lelv), intent(in)  :: u, h1, h2
            call helmholtz_matvec_2Dh_axisym(Au, u, h1, h2, shift_op, 1)
         end subroutine helm_matvec_axisym

         subroutine pres_matvec_axisym(Ap, p, h1, h2, h2inv, intype)
            ! Wrapper around pressure_matvec_2Dh_axisym matching the abstract interface.
            real(dp), dimension(lx2*ly2*lz2*lelv), intent(out) :: Ap
            real(dp), dimension(lx2*ly2*lz2*lelv), intent(in)  :: p
            real(dp), dimension(lx1*ly1*lz1*lelv), intent(in)  :: h1, h2, h2inv
            integer, intent(in) :: intype
            call pressure_matvec_2Dh_axisym(Ap, p, h1, h2, h2inv, alphaR_coef, intype)
         end subroutine pres_matvec_axisym
         
         subroutine solve_coupled_helmholtz_2Dh_axisym(x1,x2,f1,f2,h1,h2,mask1,mask2,mult,imsh,tin,maxit,
     &                                                 binv,name,shift1,shift2,couple_coef,csign)
            implicit none
            ! Same CG loop as solve_helmholtz_2Dh_axisym, but operating on the
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
            
            integer :: i, j, iter, krylov, n, nel, niter, nxyz
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
            
            real, dimension(lg) :: w1a, w2a, p1a, p2a, z1a, z2a
            
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
            call setprec(d(1),   h1,h2,imsh,1)
            call setprec(d(n+1), h1,h2,imsh,1)
            
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
               
               call copy(p1a, p(1),   n)
               call copy(p2a, p(n+1), n)
               call coupled_helmholtz_matvec_2Dh_axisym(w1a, w2a, p1a, p2a, h1, h2, shift1, shift2, couple_coef, csign)
               call dssum(w1a, lx1, ly1, lz1); call col2(w1a, mask1, n)
               call dssum(w2a, lx1, ly1, lz1); call col2(w2a, mask2, n)
               call copy(w(1),   w1a, n)
               call copy(w(n+1), w2a, n)
               
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
            !! not in the null space when alpha /= 0), matching solve_pressure_2Dh_axisym.
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
            
            call neklab_timer_start(t_axisym_solve_pres)
            
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
            
            ! tolerance: mirror solve_pressure_2Dh_axisym's chktcg2-based choice
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
                  call neklab_timer_start(t_axisym_pres_precond)
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
                  call neklab_timer_stop(t_axisym_pres_precond)
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
            call neklab_timer_stop(t_axisym_solve_pres)
         end subroutine solve_pressure_coupled_2Dh_axisym
         
      end module neklab_2Dh_axisym