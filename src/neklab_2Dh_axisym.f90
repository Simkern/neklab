      module neklab_2Dh_axisym
         use LightKrylov, only: dp
         use neklab_2Dh, only: wmask, wmask_defined, build_wmask
         use neklab_2Dh_torsion, only: if_torsion_zero,
     &                                 build_torsion_coeffs, torsion_shift_add,
     &                                 dphi_apply, dphi_apply_t,
     &                                 add_torsion_forcing_RZ, add_torsion_forcing_phi
         use neklab_nek_setup, only: nek_log_debug, nek_log_message, nek_stop_error
         use neklab_timing, only: neklab_timer_start, neklab_timer_stop, neklab_clock_tic, neklab_clock_toc,
     &                            t_advance, t_adv_setup, t_adv_advection, t_adv_rhs,
     &                            t_adv_uz, t_adv_urphi, t_adv_pressure, t_adv_correction,
     &                            t_makefp, t_solve_hmh, t_solve_cpl, t_solve_pres,
     &                            t_pres_precond, t_pres_proj, t_hmh_matvec, t_cpl_matvec,
     &                            t_pres_matvec, t_gs_comm
         implicit none
         include "SIZE"
         include "TOTAL"
      
         private
         character(len=*), parameter, private :: this_module = 'neklab_2Dh_axisym'
         
         real(dp), dimension(lx1,ly1,lz1,lelv), public :: h2z_shift
         real(dp), dimension(lx1,ly1,lz1,lelv), public :: diag_shift, couple_coef
      ! diag_shift is S_PHI = (alpha^2+1)/R^2. Torsion splits the R equation's
      ! shift away from it: S_R = (alpha^2+1+lambda^2)/R^2. In the torus the
      ! two coincide and one field sufficed; diag_shift_r carries the R value
      ! whenever lambda /= 0, and is byte-identical to diag_shift otherwise.
         real(dp), dimension(lx1,ly1,lz1,lelv), public :: diag_shift_r
         real(dp), dimension(lx1,ly1,lz1,lelv), public :: alphaR_coef   ! alpha/R, continuity/pressure coupling
         logical, public :: torus_coeffs_defined = .false.
         logical, public :: if_alpha_zero = .false.
         real(dp), private :: alpha_cached = huge(1.0_dp)

         real(dp), dimension(lx2*ly2*lz2*lelv,mxprev,lpert), private :: pbasis
         integer, dimension(lpert), private :: nprev = 0
         real(dp), private :: dtbdlast = 0.0_dp
      
         public :: build_torus_coeffs
         public :: nek_advance_2Dh_axisym
      
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
            character(len=7) :: hmh_info
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
            if_alpha_zero = (alpha == 0.0_dp)
            if (if_alpha_zero) then
               if (lpert /= 1) call nek_stop_error('alpha = 0 requires npert = 1.', this_module, this_procedure)
            else
               if (lpert /= 2) call nek_stop_error('alpha /= 0 requires npert = 2.', this_module, this_procedure)
            end if
            if (istep == 1) nprev(:) = 0

      ! Timing. The `@` regions below partition this routine exactly; the
      ! solver and kernel timers they contain are the drill-down and overlap
      ! them on purpose. neklab_clock_tic/toc accumulate an independent
      ! wall/cpu pair over the same span as the cross-check on cpu_time.
            call neklab_timer_start(t_advance)
            call neklab_clock_tic()

            call neklab_timer_start(t_adv_setup)
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
            call neklab_timer_stop(t_adv_setup)
         
            ! --- pressure-gradient / alpha-R coupling forcing (analog of gradz_p)

            msg = 'Compute explicit pressure gradient term for w'
            call nek_log_debug(msg,this_module,this_procedure)

            call neklab_timer_start(t_adv_advection)
            do jp = 1, npert
               call extrapprp(prextr)
               call compute_gradp_axisym(gradp, prextr, alphaR_coef)
               call compute_torus_s_advection(advZ, advR, advPhi, jp)
            end do
            call neklab_timer_stop(t_adv_advection)
         
            ! --- assemble RHS for both jp; solve u_Z immediately (uncoupled across jp);
            !     stash u_R/u_phi RHS for the deferred, jointly-solved transform step.

            msg = 'Solve momentum equations for u/v/w'
            call nek_log_debug(msg, this_module, this_procedure)
         
            do jp = 1, npert
               info_str = merge('Re', 'Im', jp == 1)
               call neklab_timer_start(t_adv_rhs)
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
               hmh_info = info_str//' V_Z  '
               if (istep < 10) call chktcg1(tolhv, resv1, h1, h2, v1mask, vmult, imesh, 1)
               call neklab_timer_stop(t_adv_rhs)

               call neklab_timer_start(t_adv_uz)
               call solve_helmholtz_2Dh_axisym(dv1, resv1, h1, h2, v1mask, vmult, imesh,
     &                                              tolhv, nmxv, 1, binvm1, hmh_info, h2z_shift)
               call add2(vxp(1,jp), dv1, ntot1)
               call neklab_timer_stop(t_adv_uz)

               ! stash u_R/u_phi RHS for the joint transform+solve below
               call neklab_timer_start(t_adv_rhs)
               call copy(resv2_jp(1,1,1,1,jp), resv2, ntot1)
               call copy(resv3_jp(1,1,1,1,jp), resv3, ntot1)
               call neklab_timer_stop(t_adv_rhs)
            end do ! jp

            call neklab_timer_start(t_adv_urphi)
            if (if_alpha_zero) then
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
                  call coupled_helmholtz_matvec_2Dh_axisym(w2, w3, vyp(1,jp), tp(1,1,ipert), h1, h2, diag_shift, couple_coef, csign_jp)
                  call sub2(resv2_jp(1,1,1,1,jp),    w2, ntot1)
                  call sub2(resv3_jp(1,1,1,1,ipert), w3, ntot1)
               end do
            end if
         
            msg = 'Solve u_R/u_phi momentum equations'
            call nek_log_debug(msg, this_module, this_procedure)
 
            if (if_alpha_zero) then
               do jp = 1, npert
                  info_str = merge('Re', 'Im', jp == 1)
                  call dssum(resv2_jp(1,1,1,1,jp), lx1, ly1, lz1)
                  call col2 (resv2_jp(1,1,1,1,jp), v2mask, ntot1)
                  hmh_info = info_str//' V_R  '
                  if (istep < 10) call chktcg1(tolhv, resv2_jp(1,1,1,1,jp), h1, h2, v2mask, vmult, imesh, 1)
                  call solve_helmholtz_2Dh_axisym(dv2, resv2_jp(1,1,1,1,jp), h1, h2, v2mask, vmult, imesh,
     &                                             tolhv, nmxv, 1, binvm1, hmh_info, diag_shift_r)

                  call dssum(resv3_jp(1,1,1,1,jp), lx1, ly1, lz1)
                  call col2 (resv3_jp(1,1,1,1,jp), wmask, ntot1)
                  hmh_info = info_str//' V_PHI'
                  if (istep < 10) call chktcg1(tolhv, resv3_jp(1,1,1,1,jp), h1, h2, wmask, vmult, imesh, 1)
                  call solve_helmholtz_2Dh_axisym(dv3, resv3_jp(1,1,1,1,jp), h1, h2, wmask, vmult, imesh,
     &                                             tolhv, nmxv, 1, binvm1, hmh_info, diag_shift)

                  call add2(vyp(1,jp),  dv2, ntot1)
                  call add2(tp(1,1,jp), dv3, ntot1)
               end do
            else
               do jp = 1, npert
                  ipert = npert + 1 - jp
                  csign_jp = merge(1.0_dp, -1.0_dp, jp == 1)
                  hmh_info = merge('BLK A  ', 'BLK B  ', jp == 1)
                  call solve_coupled_helmholtz_2Dh_axisym(dv2, dv3, resv2_jp(1,1,1,1,jp), resv3_jp(1,1,1,1,ipert), 
     &                                                     h1, h2, wmask, vmult, imesh, tolhv, nmxv, binvm1, hmh_info, 
     &                                                     diag_shift, couple_coef, csign_jp)
         
                  call add2(vyp(1,jp),     dv2, ntot1)
                  call add2(tp(1,1,ipert), dv3, ntot1)
               end do
            end if
            call neklab_timer_stop(t_adv_urphi)
         
            ! --- pressure correction stage: identical structure to the beta case,
            !     field routines substituted for the beta_z ones
            msg = 'Compute pressure correction to enforce mass balance'
            call nek_log_debug(msg, this_module, this_procedure)
         
            call neklab_timer_start(t_adv_pressure)
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
            if (.not. if_alpha_zero .and. ifvcor) then
               call rone(onep, ntot2)
               call pressure_matvec_2Dh_axisym(ep, onep, h1, h2, h2inv, alphaR_coef, intype)
               ebar = glsum(ep, ntot2)
            end if
         
            ifprjp = .false.
            istart = param(95)
            if (istep >= istart .and. istart /= 0) ifprjp = .true.
         
            do jp = 1, npert
               info_str = merge('Re', 'Im', jp == 1)
               call opdiv(dpr(1,jp), vxp(1,jp), vyp(1,jp), vzp(1,jp))
               call compute_frc_div_axisym(frc_div, tp, alphaR_coef)
               call add2 (dpr(1,jp), frc_div, ntot2)
               call chsign(dpr(1,jp), ntot2)

               if (if_alpha_zero) call ortho(dpr(1,jp))

               if (ifprjp) then
                  call neklab_timer_start(t_pres_proj)
                  call setrhs_pressure_2Dh_axisym(dpr(1,jp), h1, h2, h2inv, pbasis(1,1,jp), nprev(jp), alphaR_coef, info_str)
                  call neklab_timer_stop(t_pres_proj)
               end if
               call solve_pressure_2Dh_axisym(dpr(1,jp), h1, h2, h2inv, alphaR_coef, alpha, intype, iter, ebar, info_str)
               if (ifprjp) then
                  call neklab_timer_start(t_pres_proj)
                  call gensoln_pressure_2Dh_axisym(dpr(1,jp), h1, h2, h2inv, pbasis(1,1,jp), nprev(jp), alphaR_coef)
                  call neklab_timer_stop(t_pres_proj)
               end if
            end do
            call neklab_timer_stop(t_adv_pressure)
         
            msg = 'Compute velocity correction based on pressure correction'
            call nek_log_debug(msg, this_module, this_procedure)
         
            call neklab_timer_start(t_adv_correction)
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
            call neklab_timer_stop(t_adv_correction)

            call neklab_clock_toc()
            call neklab_timer_stop(t_advance)

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
            call neklab_timer_start(t_hmh_matvec)
            imesh = 1
            ntot1 = lx1*ly1*lz1*nelv
            call axhelm(Au, u, h1, h2, imesh, isd)
            call col3 (tmp, u, vdiff(1,1,1,1,1), ntot1)
            call col2 (tmp, shift, ntot1)              ! field multiply, replaces beta_z**2 scalar
            call col2 (tmp, bm1, ntot1)
            call add2 (Au, tmp, ntot1)
            call neklab_timer_stop(t_hmh_matvec)
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
            call rzero(gradp(1,ipert), ntot1)
            
            if (if_alpha_zero .and. if_torsion_zero) return

      ! --- alpha part: skew in slot space (swaps jp <-> ipert). Torus term,
      !     unchanged, absent entirely at alpha = 0.
            if (.not. if_alpha_zero) then
               spert = merge(1, -1, jp == 1)
               call mappr(gradp(1,ipert), prextr, wrk1, wrk2)
               call col2 (gradp(1,ipert), bm1, ntot1)
               call col2 (gradp(1,ipert), alphaR, ntot1)
               if (spert == -1) call chsign(gradp(1,ipert), ntot1)
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
            
            if (if_alpha_zero .and. if_torsion_zero) return
            
            if (.not. if_alpha_zero) then
               ipert = npert + 1 - jp
               spert = merge(1, -1, jp == 1)
               call mappr(dw, dpr(1,ipert), wrk1, wrk2)
               call col2 (dw, bm1, ntot1)
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
            call rzero(frc_div, ntot2)
            
            if (if_alpha_zero .and. if_torsion_zero) return
            
            nxyz1 = lx1*ly1*lz1
            call rzero(wtmp, ntot1)

      ! --- alpha part, from the opposite slot
            if (.not. if_alpha_zero) then
               ipert = npert + 1 - jp
               spert = merge(1, -1, jp == 1)
               call copy(wtmp, w(1,1,ipert), nxyz1*nelv)
               call col2(wtmp, alphaR, nxyz1*nelv)       ! scale on M1 BEFORE mapping to M2 -- see caveat below
               if (spert == -1) call chsign(wtmp, nxyz1*nelv)
            end if

      ! --- lambda part, from THIS slot: -(lambda/R)*dtheta(u_phi). Accumulated
      !     on M1 alongside the alpha part so a single map12 below carries both.
            if (.not. if_torsion_zero) then
               call dphi_apply(wl1, w(1,1,jp))
               call add2(wtmp, wl1, ntot1)
            end if

            do ie = 1, nelv
               ie1 = (ie-1)*nxyz1 + 1
               call map12(frc_div(1,1,1,ie), wtmp(ie1,1,1,1), ie)
            end do
            call col2(frc_div, bm2, ntot2)
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

            call neklab_timer_start(t_makefp)
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
            call neklab_timer_stop(t_makefp)
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
            
            if (if_alpha_zero) return
            
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
         
            call neklab_timer_start(t_pres_matvec)
            call cdabdtp(Ap, wp, h1, h2, h2inv, intype)   ! core Nek, already ifaxis-correct globally

            if (if_alpha_zero .and. if_torsion_zero) then
               call neklab_timer_stop(t_pres_matvec)
               return
            end if
         
            ntot1 = lx1*ly1*lz1*nelv
            ntot2 = lx2*ly2*lz2*nelv

      ! --- alpha block: (alpha/R) H^-1 (alpha/R). The two slot swaps of
      !     D_phi and D_phi^T cancel, so this block is diagonal in slot space
      !     and each perturbation's pressure can be solved on its own.
            if (.not. if_alpha_zero) then
               call mappr(wdivm1, wp, wrk1, wrk2)
               call col2 (wdivm1, bm1, ntot1)
               call col2 (wdivm1, alphaR, ntot1)          ! alpha/R
               call col2 (wdivm1, wmask, ntot1)
               call dssum(wdivm1, lx1, ly1, lz1)
               call col2 (wdivm1, binvm1, ntot1)
               call col2 (wdivm1, h2inv, ntot1)
               call col2 (wdivm1, alphaR, ntot1)          ! alpha/R
               call chsign(wdivm1, ntot1)
               do ie = 1, nelv
                  call map12(wdivm2(1,1,1,ie), wdivm1(1,1,1,ie), ie)
               end do
               call col2(wdivm2, bm2, ntot2)
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
            call neklab_timer_stop(t_pres_matvec)
         end subroutine pressure_matvec_2Dh_axisym

         subroutine coupled_helmholtz_matvec_2Dh_axisym(y1, y2, x1, x2, h1, h2, diag_shift, couple_coef, csign)
            implicit none
            ! Matvec for one of the two coupled 2x2 blocks:
            !   [ D   csign*c ] [x1]   [y1]
            !   [csign*c   D  ] [x2] = [y2]
            ! D = axhelm + diag_shift (== h2z_shift + 1/R^2, field-weighted)
            ! c = couple_coef (== 2*alpha/R^2, field-weighted)
            real(dp), dimension(lx1,ly1,lz1,lelv), intent(out) :: y1, y2
            real(dp), dimension(lx1,ly1,lz1,lelv), intent(in)  :: x1, x2, h1, h2, diag_shift, couple_coef
            real(dp), intent(in) :: csign
            real(dp), dimension(lx1,ly1,lz1,lelv) :: tmp
            integer :: ntot1
            ntot1 = lx1*ly1*lz1*nelv
   
            call neklab_timer_start(t_cpl_matvec)
            call helmholtz_matvec_2Dh_axisym(y1, x1, h1, h2, diag_shift, 1)
            call helmholtz_matvec_2Dh_axisym(y2, x2, h1, h2, diag_shift, 1)
   
            call col3   (tmp, x2, couple_coef, ntot1)
            call col2   (tmp, vdiff(1,1,1,1,1), ntot1)
            call col2   (tmp, bm1, ntot1)
            call add2s2 (y1, tmp, csign, ntot1)
   
            call col3   (tmp, x1, couple_coef, ntot1)
            call col2   (tmp, vdiff(1,1,1,1,1), ntot1)
            call col2   (tmp, bm1, ntot1)
            call add2s2 (y2, tmp, csign, ntot1)
            call neklab_timer_stop(t_cpl_matvec)
         end subroutine coupled_helmholtz_matvec_2Dh_axisym
   
         subroutine solve_coupled_helmholtz_2Dh_axisym(x1,x2,f1,f2,h1,h2,mask,mult,imsh,tin,maxit,binv,name,diag_shift,couple_coef,csign)
            implicit none
            ! Same CG loop as solve_helmholtz_2Dh_axisym, but operating on the
            ! stacked 2-field vector (x1;x2). dssum/mask applied to x1 and x2
            ! SEPARATELY (mirrors the planar-case pattern), never combined
            ! before dssum -- this is the ordering question from before, tested
            ! directly by construction.
            include 'FDMH1'
            real, intent(out) :: x1(1), x2(1)
            real, intent(inout)  :: f1(1), f2(1)
            real, intent(in)  :: h1(1), h2(1), mask(1), mult(1), binv(1)
            integer, intent(in) :: imsh, maxit
            real, intent(in) :: tin
            character(len=7), intent(in) :: name
            real, intent(in) :: diag_shift(1), couple_coef(1)
            real(dp), intent(in) :: csign
   
            integer :: i, j, iter, krylov, n, nel, niter, nxyz
            real :: alpha, beta, alphm, fmax
            real :: rbn0, rbn2, rho, rho0, rtz1, rtz2
            real :: tol, vol
            LOGICAL IFPRINT, IFHZPC
            COMMON /CPRINT/ IFPRINT, IFHZPC
            logical ifdfrm, iffast, ifh2, ifsolv
            common /fastmd/ ifdfrm(lelt), iffast(lelt), ifh2, ifsolv
            logical ifprint_hmh
            integer, parameter :: lg = lx1*ly1*lz1*lelt
            ! stacked scratch: [1:n] = field 1, [n+1:2n] = field 2
            real, dimension(2*lg) :: d, r, w, p, z
            real :: scalar(2)
            integer :: niterhm
            common /iterhm/ niterhm
            real, external :: glamax, glsum, glsc2, glsc3, vlsc3, vlsc32
   
            real, dimension(lg) :: w1a, w2a, p1a, p2a, z1a, z2a
   
            call neklab_timer_start(t_solve_cpl)
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
   
            call neklab_timer_start(t_gs_comm)
            call dssum(f1, lx1, ly1, lz1)
            call dssum(f2, lx1, ly1, lz1)
            call neklab_timer_stop(t_gs_comm)
            call copy(r(1),   f1, n);  call col2(r(1),   mask, n)
            call copy(r(n+1), f2, n);  call col2(r(n+1), mask, n)
   
            call rzero(x1,n); call rzero(x2,n)
            call rzero(p,2*n)
   
            fmax = max(glamax(r(1),n), glamax(r(n+1),n))
            if (fmax == 0.0) then
      ! Zero RHS: nothing to solve. The stop MUST happen here -- a timer left
      ! running is not restarted by the next `start` (Timer_Utils.f90:185), so
      ! a missed stop silently loses every subsequent call as well.
               call neklab_timer_stop(t_solve_cpl)
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
               call neklab_timer_start(t_gs_comm)
               call gop(scalar,w,'+  ',2)
               call neklab_timer_stop(t_gs_comm)
               rtz1 = scalar(1)
               rbn2 = sqrt(scalar(2)/(2.0*vol))
               if (iter.eq.1) rbn0 = rbn2
               if (param(22).lt.0) tol = abs(param(22))*rbn0
               if (tin.lt.0)       tol = abs(tin)*rbn0
   
               ifprint_hmh = .false.
               if (nio.eq.0.and.ifprint.and.param(74).ne.0) ifprint_hmh=.true.
               if (nio.eq.0.and.istep.eq.1)                 ifprint_hmh=.true.
               if (ifprint_hmh) write(6,3002) istep,'  CplHmh ' // name, iter, rbn2, h1(1), tol, h2(1)
   
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
               call coupled_helmholtz_matvec_2Dh_axisym(w1a, w2a, p1a, p2a, h1, h2, diag_shift, couple_coef, csign)
               call neklab_timer_start(t_gs_comm)
               call dssum(w1a, lx1, ly1, lz1); call col2(w1a, mask, n)
               call dssum(w2a, lx1, ly1, lz1); call col2(w2a, mask, n)
               call neklab_timer_stop(t_gs_comm)
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
            call neklab_timer_stop(t_solve_cpl)
         end subroutine solve_coupled_helmholtz_2Dh_axisym

         subroutine solve_helmholtz_2Dh_axisym(x,f,h1,h2,mask,mult,imsh,tin,maxit,isd,binv,name,shift)
            implicit none
            include 'FDMH1'
            real, intent(out) :: x(1)
            real, intent(in) :: f(1), h1(1), h2(1), mask(1), mult(1), binv(1)
            integer, intent(in) :: imsh, maxit, isd
            real, intent(in) :: tin
            character(len=7), intent(in) :: name
            real, intent(in) :: shift(1)
            integer :: i, j, iter, krylov, n, nel, niter, nxyz
            real :: alpha, beta, alphm, fmax, h2max
            real :: rbn0, rbn2, rho, rho0, rmean, rtz1, rtz2
            real :: skmin, smean, tol, vol, div0, ratio, divex
            real :: etime2, etime_p
            LOGICAL IFPRINT, IFHZPC
            COMMON /CPRINT/ IFPRINT, IFHZPC
            logical ifdfrm, iffast, ifh2, ifsolv
            common /fastmd/ ifdfrm(lelt), iffast(lelt), ifh2, ifsolv
            logical ifmcor,ifprint_hmh
            integer, parameter :: lg = lx1*ly1*lz1*lelt
            real, dimension(lg) :: d, r, w, p, z
            real :: scalar(2)
            COMMON /SCRCG/ d, scalar
            common /SCRMG/ r, w, p, z
            integer, parameter :: maxcg = 900
            real, dimension(maxcg) :: diagt, upper
            common /tdarray/ diagt, upper
            integer :: niterhm
            common /iterhm/ niterhm
            real, external :: glamax, glmax, glmin, glsum, glsc2, glsc3, vlsc3, vlsc32
         
            call neklab_timer_start(t_solve_hmh)
            call rzero(diagt,maxcg)
            call rzero(upper,maxcg)
            rho = 0.00
            NXYZ = lx1*ly1*lz1
            NEL  = NELV
            VOL  = VOLVM1
            IF (IMSH.EQ.2) NEL=NELT
            IF (IMSH.EQ.2) VOL=VOLTM1
            n = NEL*NXYZ
            tol = abs(tin)
            if (restol(ifield).ne.0) tol = abs(restol(ifield))
            if (tin.lt.0) tol=abs(tin)
            niter = min(maxit,maxcg)
         
            if (.not.ifsolv) then
               call setfast(h1,h2,imsh)
               ifsolv = .true.
            endif
            call setprec(D,h1,h2,imsh,isd)   ! NOTE: preconditioner still ignores `shift` -- see caveat below
            call copy (r,f,n)
            call rzero(x,n)
            call rzero(p,n)
            fmax = glamax(f,n)
            if (fmax == 0.0) then
               call neklab_timer_stop(t_solve_hmh)
               return
            end if
         
            krylov = 0
            rtz1=1.0
            niterhm = 0
         
            do iter=1,niter
               call col3(z,r,d,n)
               rtz2=rtz1
               scalar(1)=vlsc3(z,r,mult,n)
               if(param(18).eq.1) then
                 scalar(2)=vlsc3(r,r,mult,n)
               else
                scalar(2)=vlsc32(r,mult,binv,n)
               endif
               call neklab_timer_start(t_gs_comm)
               call gop(scalar,w,'+  ',2)
               call neklab_timer_stop(t_gs_comm)
               rtz1=scalar(1)
               rbn2=sqrt(scalar(2)/vol)
               if (iter.eq.1) rbn0 = rbn2
               if (param(22).lt.0) tol=abs(param(22))*rbn0
               if (tin.lt.0)       tol=abs(tin)*rbn0
         
               ifprint_hmh = .false.
               if (nio.eq.0.and.ifprint.and.param(74).ne.0) ifprint_hmh=.true.
               if (nio.eq.0.and.istep.eq.1)                 ifprint_hmh=.true.
               if (ifprint_hmh) write(6,3002) istep,'  Hmholtz ' // name,  iter,rbn2,h1(1),tol,h2(1),ifmcor
         
               IF (rbn2.LE.TOL.and.(iter.gt.1 .or. istep.le.5)) THEN
               NITER = ITER-1
               if (nio.eq.0) write(6,3000) istep,'  Hmholtz ' // name, niter,rbn2,rbn0,tol
                  goto 9999
               endif
         
               beta = rtz1/rtz2
               if (iter.eq.1) beta=0.0
               call add2s1 (p,z,beta,n)
               call helmholtz_matvec_2Dh_axisym(w,p,h1,h2,shift,isd)
               call neklab_timer_start(t_gs_comm)
               call dssum  (w,lx1,ly1,lz1)
               call neklab_timer_stop(t_gs_comm)
               call col2   (w,mask,n)
         
               rho0 = rho
               rho  = glsc3(w,p,mult,n)
               alpha=rtz1/rho
               alphm=-alpha
               call add2s2(x,p ,alpha,n)
               call add2s2(r,w ,alphm,n)
         
               if (iter.eq.1) then
                  krylov = krylov+1
                  diagt(iter) = rho/rtz1
               elseif (iter.le.maxcg) then
                  krylov = krylov+1
                  diagt(iter)    = (beta**2 * rho0 + rho ) / rtz1
                  upper(iter-1)  = -beta * rho0 / sqrt(rtz2 * rtz1)
               endif
            enddo
            niter = iter-1
            if (nio.eq.0) write (6,3001) istep, '  Error Hmholtz ' // name, niter,rbn2,rbn0,tol
 3000       format(i11,a,1x,I7,1p4E13.4)
 3001       format(i11,a,1x,I7,1p4E13.4)
 3002       format(i11,a,1x,I7,1p4E13.4,l4)
 9999       continue
            niterhm = niter
            ifsolv = .false.
            call neklab_timer_stop(t_solve_hmh)
         end subroutine solve_helmholtz_2Dh_axisym

         subroutine solve_pressure_2Dh_axisym(res,h1,h2,h2inv,alphaR,alpha,intype,iter,ebar,info_str)
            implicit none
            include 'GMRES'
            real :: divex
            common  /ctolpr/ divex
            logical          ifprint
            common  /cprint/ ifprint
            real, dimension(lx2*ly2*lz2*lelv), intent(inout) :: res
            real, dimension(lx1,ly1,lz1,lelv), intent(in)    :: h1, h2, h2inv, alphaR
            real(dp), intent(in) :: alpha        ! the toroidal wavenumber itself, for the ortho guard
            integer, intent(in) :: intype
            integer, intent(inout) :: iter
            real(dp), intent(in) :: ebar
            character(len=2), optional, intent(in) :: info_str
            real, dimension(lx2,ly2,lz2,lelv) :: wp
            common /scrmg/   wp
            real, dimension(lgmres) :: y, wk1, wk2
            common /ctmp0/   wk1, wk2
            common /cgmres1/ y
            real alph, l, temp, div0, ratio, rnorm, tolpss
            real :: etime2, etime_p
            integer i, j, k, m, iconv, ntot2
            logical iflag
            save    iflag
            data    iflag /.false./
            real    norm_fac
            save    norm_fac
            real*8 etime1,dnekclock
            real, external :: vlsc2, glsc2, glsum
            integer, parameter :: gmres_imax = 1000
         
            call neklab_timer_start(t_solve_pres)
            if(.not.iflag) then
               iflag=.true.
               call uzawa_gmres_split0(ml_gmres,mu_gmres,bm2,bm2inv,lx2*ly2*lz2*nelv)
               norm_fac = 1./sqrt(volvm2)
            endif
         
            etime1 = dnekclock()
            etime_p = 0.
            divex = 0.
            iter  = 0
            m = lgmres
            call chktcg2(tolps,res,iconv)
            if (param(21).gt.0.and.tolps.gt.abs(param(21))) tolps = abs(param(21))
            if (istep.eq.0) tolps = 1.e-4
            tolpss = tolps
            ntot2  = lx2*ly2*lz2*nelv
            iconv = 0
            call rzero(x_gmres,ntot2)
         
            do while(iconv.eq.0.and.iter.lt.gmres_imax)
               if(iter.eq.0) then
                  call col3(r_gmres,ml_gmres,res,ntot2)
               else
                  call copy(r_gmres,res,ntot2)
                  call pressure_matvec_2Dh_axisym(w_gmres,x_gmres,h1,h2,h2inv,alphaR,intype)
                  call add2s2(r_gmres,w_gmres,-1.,ntot2)
                  call col2(r_gmres,ml_gmres,ntot2)
               endif
               gamma_gmres(1) = sqrt(glsc2(r_gmres,r_gmres,ntot2))
               if(iter.eq.0) then
                  div0 = gamma_gmres(1)*norm_fac
                  if (param(21).lt.0) tolpss=abs(param(21))*div0
               endif
               rnorm = 0.
               if(gamma_gmres(1) .eq. 0.) goto 9000
               temp = 1./gamma_gmres(1)
               call cmult2(v_gmres(1,1),r_gmres,temp,ntot2)
               do j=1,m
                  iter = iter+1
                  call col3(w_gmres,mu_gmres,v_gmres(1,j),ntot2)
                  etime2 = dnekclock()
                  call neklab_timer_start(t_pres_precond)
                  call hsmg_solve(z_gmres(1,j),w_gmres)
                  if (ebar /= 0.0_dp) then
                     call cadd(z_gmres(1,j), glsum(w_gmres,ntot2)/ebar, ntot2)
                  end if
                  call neklab_timer_stop(t_pres_precond)
                  etime_p = etime_p + dnekclock()-etime2
                  call pressure_matvec_2Dh_axisym(w_gmres,z_gmres(1,j),h1,h2,h2inv,alphaR,intype)
                  call col2(w_gmres,ml_gmres,ntot2)
                  do i=1,j
                     h_gmres(i,j)=vlsc2(w_gmres,v_gmres(1,i),ntot2)
                  enddo
                  call gop(h_gmres(1,j),wk1,'+  ',j)
                  do i=1,j
                     call add2s2(w_gmres,v_gmres(1,i),-h_gmres(i,j),ntot2)
                  enddo
                  do i=1,j-1
                     temp = h_gmres(i,j)
                     h_gmres(i  ,j)=  c_gmres(i)*temp + s_gmres(i)*h_gmres(i+1,j)
                     h_gmres(i+1,j)= -s_gmres(i)*temp + c_gmres(i)*h_gmres(i+1,j)
                  enddo
                  alph = sqrt(glsc2(w_gmres,w_gmres,ntot2))
                  rnorm = 0.
                  if(alph.eq.0.) goto 900
                  l = sqrt(h_gmres(j,j)*h_gmres(j,j)+alph*alph)
                  temp = 1./l
                  c_gmres(j) = h_gmres(j,j) * temp
                  s_gmres(j) = alph  * temp
                  h_gmres(j,j) = l
                  gamma_gmres(j+1) = -s_gmres(j) * gamma_gmres(j)
                  gamma_gmres(j)   =  c_gmres(j) * gamma_gmres(j)
                  rnorm = abs(gamma_gmres(j+1))*norm_fac
                  ratio = rnorm/div0
                  if (ifprint.and.nio.eq.0) write (6,66) iter,tolpss,rnorm,div0,ratio,istep
   66             format(i5,1p4e12.5,i8,' Divergence')
                  if (rnorm .lt. tolpss) goto 900
                  if (j.eq.m) goto 1000
                  temp = 1./alph
                  call cmult2(v_gmres(1,j+1),w_gmres,temp,ntot2)
               enddo
  900          iconv = 1
 1000          continue
               do k=j,1,-1
                  temp = gamma_gmres(k)
                  do i=j,k+1,-1
                     temp = temp - h_gmres(k,i)*c_gmres(i)
                  enddo
                  c_gmres(k) = temp/h_gmres(k,k)
               enddo
               do i=1,j
                  call add2s2(x_gmres,z_gmres(1,i),c_gmres(i),ntot2)
               enddo
            enddo
 9000       continue
            divex = rnorm
         
            if (alpha == 0.0_dp) call ortho  (w_gmres)
            call copy(res,x_gmres,ntot2)
            if (alpha == 0.0_dp) call ortho (res)
         
            etime1 = dnekclock()-etime1
            if (present(info_str)) then
               if (nio.eq.0) write(6,9998) istep,'  U-PRES  '//info_str//' gmres  alpha= ', alpha, iter,divex,div0,tolpss,etime_p,etime1
            else
               if (nio.eq.0) write(6,9999) istep,'  U-PRES  gmres  ', iter,divex,div0,tolpss,etime_p,etime1
            end if
 9998       format(i11,a,F5.1,1X,I6,1p5e13.4)
 9999       format(i11,a,I6,1p5e13.4)
            call neklab_timer_stop(t_solve_pres)
         end subroutine solve_pressure_2Dh_axisym

         subroutine setrhs_pressure_2Dh_axisym(dp,h1,h2,h2inv,proj_set,niprev,alphaR,info_str)
            implicit none
            !
            !     Project soln onto best fit in the "E" norm.
            !
            real, dimension(lx2,ly2,lz2,lelv), intent(inout) :: dp
            real, dimension(lx1,ly1,lz1,lelv), intent(in) :: h1
            real, dimension(lx1,ly1,lz1,lelv), intent(in) :: h2
            real, dimension(lx1,ly1,lz1,lelv), intent(in) :: h2inv
            ! projection basis
            real, dimension(lx2*ly2*lz2*lelv,mxprev), intent(inout) :: proj_set
            integer, intent(inout) :: niprev
            ! alpha/R coupling field (replaces scalar beta_z of the planar case)
            real, dimension(lx1,ly1,lz1,lelv), intent(in) :: alphaR
            character(len=2), intent(in) :: info_str
            ! common block
            integer, parameter :: ltot2 = lx2*ly2*lz2*lelv
            real :: pbar(ltot2), pnew(ltot2)
            real :: alpha(mxprev), work(mxprev)
            common /orthox/ pbar, pnew
            common /orthos/ alpha, work
            ! internal
            integer :: ntot2, n10, intype, i
            real, external :: glsc3, vlsc2
            real :: alpha1, alpha2, ratio

            ntot2 = lx2*ly2*lz2*nelv
            call rzero(pbar,ntot2)
            if (niprev == 0) return

            ! Diag to see how much reduction in the residual is attained.
            alpha1 = glsc3(dp,dp,bm2inv,ntot2)
            if (alpha1.gt.0) then
               alpha1 = sqrt(alpha1/volvm2)
            else
               return
            endif

            do i = 1, niprev  ! Perform Gram-Schmidt for previous soln's.
               alpha(i) = vlsc2(dp,proj_set(1,i),ntot2)
            enddo
            call gop(alpha,work,'+  ',niprev)

            do i = 1, niprev
               call add2s2(pbar,proj_set(1,i),alpha(i),ntot2)
            enddo

            intype = 1
            call pressure_matvec_2Dh_axisym(pnew,pbar,h1,h2,h2inv,alphaR,intype)
            call sub2(dp,pnew,ntot2)

            alpha2 = glsc3(dp,dp,bm2inv,ntot2) ! Diagnostics
            if (alpha2.gt.0) then
               alpha2 = sqrt(alpha2/volvm2)
               ratio  = alpha1/alpha2
               n10=min(10,niprev)
               if (nio.eq.0) write(6,13) istep,'  Project '//info_str//' PRES',
     &                         alpha2,alpha1,ratio,niprev,mxprev
            endif
 13         format(i11,a,6x,1p3e13.4,i4,i4)
         end subroutine setrhs_pressure_2Dh_axisym

         subroutine gensoln_pressure_2Dh_axisym(dp,h1,h2,h2inv,proj_set,niprev,alphaR)
            implicit none
            !
            !     Reconstruct the solution to the original problem by adding back
            !     the previous solutions
            !
            real, dimension(lx2,ly2,lz2,lelv), intent(inout) :: dp
            real, dimension(lx1,ly1,lz1,lelv), intent(in) :: h1
            real, dimension(lx1,ly1,lz1,lelv), intent(in) :: h2
            real, dimension(lx1,ly1,lz1,lelv), intent(in) :: h2inv
            ! projection basis
            real, dimension(lx2*ly2*lz2*lelv,mxprev), intent(inout) :: proj_set
            integer, intent(inout) :: niprev
            ! alpha/R coupling field (replaces scalar beta_z of the planar case)
            real, dimension(lx1,ly1,lz1,lelv), intent(in) :: alphaR
            ! common block
            integer, parameter :: ltot2 = lx2*ly2*lz2*lelv
            real :: pbar(ltot2), pnew(ltot2)
            real :: alpha(mxprev), work(mxprev)
            common /orthox/ pbar, pnew
            common /orthos/ alpha, work
            ! internal
            integer :: ntot2, mprev, ierr

            mprev=param(93)
            mprev=min(mprev,mxprev)

            ntot2=lx2*ly2*lz2*nelv

            if (niprev.lt.mprev) then
               niprev = niprev+1
               call copy (proj_set(1,niprev),dp,ntot2)        ! Save current solution
               call add2 (dp,pbar,ntot2)                      ! Reconstruct solution.
               call econj_pressure_2Dh_axisym(proj_set,niprev,h1,h2,h2inv,alphaR,ierr) ! Orthonormalize set
               if (ierr.eq.1) then
                  niprev = 1
                  call copy (proj_set(1,niprev),dp,ntot2)     ! Save current solution
                  call econj_pressure_2Dh_axisym(proj_set,niprev,h1,h2,h2inv,alphaR,ierr) ! and orthonormalize.
               endif
            else                                              !          (uses pnew).
               niprev = 1
               call add2 (dp,pbar,ntot2)                      ! Reconstruct solution.
               call copy (proj_set(1,niprev),dp,ntot2)        ! Save current solution
               call econj_pressure_2Dh_axisym(proj_set,niprev,h1,h2,h2inv,alphaR,ierr) ! and orthonormalize.
            endif
         end subroutine gensoln_pressure_2Dh_axisym

         subroutine econj_pressure_2Dh_axisym(proj_set,nprev,h1,h2,h2inv,alphaR,ierr)
            implicit none
            !     Orthogonalize the soln wrt previous soln's for which we already
            !     know the soln.
            real, dimension(lx2*ly2*lz2*lelv,mxprev), intent(inout) :: proj_set
            integer, intent(inout) :: nprev
            real, dimension(lx1,ly1,lz1,lelv), intent(in) :: h1
            real, dimension(lx1,ly1,lz1,lelv), intent(in) :: h2
            real, dimension(lx1,ly1,lz1,lelv), intent(in) :: h2inv
            ! alpha/R coupling field (replaces scalar beta_z of the planar case)
            real, dimension(lx1,ly1,lz1,lelv), intent(in) :: alphaR
            integer, intent(out) :: ierr
            ! common block
            integer, parameter :: ltot2 = lx2*ly2*lz2*lelv
            real :: pbar(ltot2), pnew(ltot2)
            real :: alpha(mxprev), work(mxprev)
            common /orthox/ pbar,pnew
            common /orthos/ alpha,work
            ! internal
            integer :: ntot2, i, ipass, npass, nprev1, intype
            real :: alphad, alpham
            real, external :: vlsc2, glsc2

            ierr  = 0
            ntot2 = lx2*ly2*lz2*nelv

            !
            !     Gram Schmidt, w re-orthogonalization
            !
            npass = 1
            do ipass = 1, npass

               intype = 1
               call pressure_matvec_2Dh_axisym(pnew,proj_set(1,nprev),h1,h2,h2inv,alphaR,intype)
               alphad = glsc2(pnew,proj_set(1,nprev),ntot2) ! compute part of the norm

               nprev1 = nprev - 1
               do i = 1, nprev1   !   Gram-Schmidt
                  alpha(i) = vlsc2(pnew,proj_set(1,i),ntot2)
               enddo
               if (nprev1.gt.0) call gop(alpha,work,'+  ',nprev1)

               do i = 1, nprev1
                  alpham = -alpha(i)
                  call add2s2(proj_set(1,nprev),proj_set(1,i),alpham,ntot2)
                  alphad = alphad - alpha(i)**2
               enddo

            enddo
            !
            !    Normalize new element in P~
            !
            if (alphad.le.0) then
               write(6,*) 'ERROR:  alphad .le. 0 in econj_pressure_2Dh_axisym',alphad,nprev
               ierr = 1
               return
            endif
            alphad = 1./sqrt(alphad)
            call cmult(proj_set(1,nprev),alphad,ntot2)

         end subroutine econj_pressure_2Dh_axisym

      end module neklab_2Dh_axisym