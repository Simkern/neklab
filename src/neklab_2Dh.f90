      module neklab_2Dh
      !! The two 2Dh linear steppers.
      !!
      !! Each stepper is a short sequence of stage calls named after Nek's
      !! perturbation solver (perturb.f: makefp, cresvipp, ophx, ophinv,
      !! incomprp). Model-specific stages live in neklab_2Dh_planar and
      !! neklab_2Dh_axisym, shared stages in neklab_2Dh_common.
      !!
      !! A loop over the perturbation slots jp is split only where a stage reads
      !! the other slot ipert = npert + 1 - jp:
      !!   makecpl_*       reads jp, writes ipert  (slot-coupled explicit terms)
      !!   ophx/ophinv     u_R/u_phi coupling at alpha /= 0 (axisym only)
      !!   incomprp_div_*  reads the other slot's third velocity component
      !!   incomprp_dw_*   reads the other slot's pressure correction
         use LightKrylov, only: dp
         use neklab_2Dh_common, only: setup_step_2Dh, cresvipp_2Dh, setup_incomprp_2Dh,
     &                                ebar_2Dh, incomprp_solve_2Dh, incomprp_update_2Dh,
     &                                solve_pressure_2Dh, setrhs_pressure_2Dh,
     &                                gensoln_pressure_2Dh, bm2_map12_2Dh
         use neklab_2Dh_planar, only: prepare_step_2Dh_planar, make_coupling_2Dh_planar,
     &                                makefp_2Dh_planar, ophx_2Dh_planar, ophinv_2Dh_planar,
     &                                incomprp_div_2Dh_planar, incomprp_dw_2Dh_planar,
     &                                pres_matvec_planar, maxit_p_planar
         use neklab_2Dh_axisym, only: init_2Dh_axisym,
     &                                prepare_step_2Dh_axisym, make_coupling_2Dh_axisym,
     &                                makefp_2Dh_axisym, ophx_2Dh_axisym,
     &                                ophinv_uz_2Dh_axisym, ophinv_urphi_2Dh_axisym,
     &                                incomprp_div_2Dh_axisym, incomprp_dw_2Dh_axisym,
     &                                incomprp_coupled_2Dh_axisym,
     &                                pres_matvec_axisym, maxit_p_axisym, alpha_eq_zero
         use neklab_2Dh_torsion, only: if_torsion_zero, refresh_torsion_forcing_nl,
     &                                 add_torsion_pressure_phi_nl,
     &                                 dphi_apply, lambda_correction, add_lambda_block_E,
     &                                 torsion_lambda
         use neklab_nek_setup, only: nek_stop_error
         use neklab_timing, only: neklab_timer_start, neklab_timer_stop,
     &                            neklab_clock_tic, neklab_clock_toc,
     &                            t_planar_advance, t_planar_adv_setup, t_planar_adv_advection,
     &                            t_planar_adv_rhs, t_planar_adv_velocity, t_planar_adv_pressure,
     &                            t_planar_adv_correction,
     &                            t_axisym_advance, t_axisym_adv_setup, t_axisym_adv_advection,
     &                            t_axisym_adv_rhs, t_axisym_adv_uz, t_axisym_adv_urphi,
     &                            t_axisym_adv_pressure, t_axisym_adv_correction
         implicit none
         include "SIZE"
         include "TOTAL"
         include "ADJOINT" ! ifadj

         private
         character(len=*), parameter, private :: this_module = 'neklab_2Dh'
         integer, parameter :: lv = lx1*ly1*lz1*lelv
         integer, parameter :: lp = lx2*ly2*lz2*lelv

      ! Fields handed from one stage to the next.
      ! m1
         real(dp), dimension(lv,lpert),   private :: gradw    ! third-component pressure gradient (make_coupling -> cresvipp)
         real(dp), dimension(lv,3,lpert), private :: resv_jp  ! momentum residuals per slot (axisym: ophx after the slot loop)
         real(dp), dimension(lv),         private :: dw       ! third-component pressure correction
      ! m2
         real(dp), dimension(lp,lpert),   private :: dpr      ! pressure correction
      ! Nonlinear helix step (fluid_2Dh_helix, incomprn_2Dh_helix)
         real(dp), dimension(lv,3), private :: nl_resv, nl_dv
         real(dp), dimension(lv),   private :: nl_h1, nl_h2, nl_h2inv, nl_w
         real(dp), dimension(lp),   private :: nl_dp, nl_w2
         real(dp), private :: pbasis_nl(lp,mxprev)
         integer,  private :: nprev_nl = 0
         real(dp), private :: dtbdlast_nl = 0.0_dp

         public :: nek_advance_2Dh_planar
         public :: nek_advance_2Dh_axisym
         public :: nek_advance_2Dh_axisym_nl

      contains

         subroutine nek_advance_2Dh_planar(betaz, if_advW_in)
      !! One linear step of the cartesian 2Dh problem, spanwise wavenumber betaz.
            real(dp), intent(in) :: betaz
            logical, optional, intent(in) :: if_advW_in
            real(dp), dimension(lx1,ly1,lz1,lelv) :: w1, w2, w3, resv1, resv2, resv3, dv1, dv2, dv3
            real(dp), dimension(lx1,ly1,lz1,lelv) :: h1, h2, h2inv
            common /SCRNS/  w1, w2, w3, resv1, resv2, resv3, dv1, dv2, dv3
            common /SCRVH/  h1
            common /SCRCH/  h2
            common /scrhi/  h2inv
            real(dp) :: ebar
            integer :: igeom, intype, ntot1

            ntot1 = lx1*ly1*lz1*nelv
            call prepare_step_2Dh_planar(betaz, if_advW_in)    ! init check, betaz, W-advection switch

            call neklab_timer_start(t_planar_advance)
            call neklab_clock_tic()

            call neklab_timer_start(t_planar_adv_setup)
            call setup_step_2Dh()
            call neklab_timer_stop(t_planar_adv_setup)

      ! --- explicit slot-coupled terms: read jp, write ipert
            call neklab_timer_start(t_planar_adv_advection)
            do jp = 1, npert
               call make_coupling_2Dh_planar(gradw)              ! gradz_p, W-advection
            end do
            call neklab_timer_stop(t_planar_adv_advection)

      ! --- momentum, per slot
            do jp = 1, npert
               call neklab_timer_start(t_planar_adv_rhs)
               igeom = 1                                         ! old geometry, old velocity
               call makefp_2Dh_planar()                          ! makefp + lagfieldp + makeqp + lagscalp
               igeom = 2                                         ! new geometry, new velocity
               intype = -1
               ifield = 1
               call sethlm(h1, h2, intype)
               if (ifadj) call seth2bc(h2)   ! adjoint O/ON boundary term
               call cresvipp_2Dh(resv1, resv2, resv3, gradw(1,jp))
               call ophx_2Dh_planar(resv1, resv2, resv3, h1, h2)
               call neklab_timer_stop(t_planar_adv_rhs)

               call neklab_timer_start(t_planar_adv_velocity)
               call ophinv_2Dh_planar(dv1, dv2, dv3, resv1, resv2, resv3, h1, h2)
               call opadd2(vxp(1,jp), vyp(1,jp), vzp(1,jp), dv1, dv2, dv3)
               call add2  (tp(1,1,jp), dv3, ntot1)
               call neklab_timer_stop(t_planar_adv_velocity)
            end do

      ! --- incompressibility
            call neklab_timer_start(t_planar_adv_pressure)
            call setup_incomprp_2Dh(h1, h2, h2inv)
            call ebar_2Dh(ebar, pres_matvec_planar, h1, h2, h2inv, betaz)
            do jp = 1, npert
               call incomprp_div_2Dh_planar(dpr(1,jp))           ! -(D u + D_beta w), reads ipert
               call incomprp_solve_2Dh(dpr(1,jp), h1, h2, h2inv, pres_matvec_planar, ebar,
     &                                 betaz, maxit_p_planar, 'beta')
            end do
            call neklab_timer_stop(t_planar_adv_pressure)

            call neklab_timer_start(t_planar_adv_correction)
            do jp = 1, npert
               call incomprp_dw_2Dh_planar(dw, dpr, h2inv)       ! w correction, reads ipert
               call incomprp_update_2Dh(dpr(1,jp), dw, h2inv)    ! u, v, w and p updates
            end do
            call neklab_timer_stop(t_planar_adv_correction)

            call neklab_clock_toc()
            call neklab_timer_stop(t_planar_advance)
            jp = 0   ! back to the nonlinear solver's slot
         end subroutine nek_advance_2Dh_planar

         subroutine nek_advance_2Dh_axisym(alpha)
      !! One linear step of the axisymmetric 2Dh problem, azimuthal wavenumber alpha.
            real(dp), intent(in) :: alpha
            real(dp), dimension(lx1,ly1,lz1,lelv) :: h1, h2, h2inv
            common /SCRVH/  h1
            common /SCRCH/  h2
            common /scrhi/  h2inv
            real(dp) :: ebar
            integer :: igeom, intype

            call prepare_step_2Dh_axisym(alpha)                 ! coefficients, npert check

            call neklab_timer_start(t_axisym_advance)
            call neklab_clock_tic()

            call neklab_timer_start(t_axisym_adv_setup)
            call setup_step_2Dh()
            call neklab_timer_stop(t_axisym_adv_setup)

      ! --- explicit slot-coupled terms: read jp, write ipert
            call neklab_timer_start(t_axisym_adv_advection)
            do jp = 1, npert
               call make_coupling_2Dh_axisym(gradw)               ! gradp, swirl advection, torsion alpha*lambda
            end do
            call neklab_timer_stop(t_axisym_adv_advection)

      ! --- momentum right-hand sides, per slot
            call neklab_timer_start(t_axisym_adv_rhs)
            do jp = 1, npert
               igeom = 1                                         ! old geometry, old velocity
               call makefp_2Dh_axisym()                          ! incl. curvature terms
               igeom = 2                                         ! new geometry, new velocity
               intype = -1
               ifield = 1
               call sethlm(h1, h2, intype)
               if (ifadj) call seth2bc(h2)   ! adjoint O/ON boundary term
               call cresvipp_2Dh(resv_jp(1,1,jp), resv_jp(1,2,jp), resv_jp(1,3,jp), gradw(1,jp))
            end do
      ! --- H u: u_Z per slot; u_R/u_phi coupled across the slots at alpha /= 0
            call ophx_2Dh_axisym(resv_jp, h1, h2)
            call neklab_timer_stop(t_axisym_adv_rhs)

            call neklab_timer_start(t_axisym_adv_uz)
            call ophinv_uz_2Dh_axisym(resv_jp, h1, h2)           ! solves and updates vxp
            call neklab_timer_stop(t_axisym_adv_uz)

            call neklab_timer_start(t_axisym_adv_urphi)
            call ophinv_urphi_2Dh_axisym(resv_jp, h1, h2)        ! solves and updates vyp, tp
            call neklab_timer_stop(t_axisym_adv_urphi)

      ! --- incompressibility
            call neklab_timer_start(t_axisym_adv_pressure)
            call setup_incomprp_2Dh(h1, h2, h2inv)
            call ebar_2Dh(ebar, pres_matvec_axisym, h1, h2, h2inv, alpha)
            if (.not. alpha_eq_zero .and. .not. if_torsion_zero) then
               call incomprp_coupled_2Dh_axisym(dpr, h1, h2, h2inv, alpha, ebar)   ! deleted in Step 7
            else
               do jp = 1, npert
                  call incomprp_div_2Dh_axisym(dpr(1,jp))        ! -(D u + D_alpha u_phi), reads ipert
                  call incomprp_solve_2Dh(dpr(1,jp), h1, h2, h2inv, pres_matvec_axisym, ebar,
     &                                    alpha, maxit_p_axisym, 'alpha')
               end do
            end if
            call neklab_timer_stop(t_axisym_adv_pressure)

            call neklab_timer_start(t_axisym_adv_correction)
            do jp = 1, npert
               call incomprp_dw_2Dh_axisym(dw, dpr, h2inv)       ! u_phi correction, reads ipert
               call incomprp_update_2Dh(dpr(1,jp), dw, h2inv)    ! u_Z, u_R, u_phi and p updates
            end do
            call neklab_timer_stop(t_axisym_adv_correction)

            call neklab_clock_toc()
            call neklab_timer_stop(t_axisym_advance)
            jp = 0   ! back to the nonlinear solver's slot
         end subroutine nek_advance_2Dh_axisym

         subroutine nek_advance_2Dh_axisym_nl()
      !! One NONLINEAR step of the axisymmetric 2Dh problem (torus, helix).
      !! A copy of the P_N-P_{N-2} branch of Nek's nek_advance (drive1.f,
      !! v19) with the torsion hooks:
      !!   - the explicit torsion forcing, refreshed from U^n before the step;
      !!   - +D_lam^T p_ext into the u_phi right-hand side, after heat(1);
      !!   - fluid(2) replaced by fluid_2Dh_helix: Nek's in-plane predictor,
      !!     then ONE projection with the full E, correcting u_Z, u_R AND u_phi.
      !! At lambda = 0 every hook returns at once and fluid(2) is Nek's: this
      !! is nek_advance.
            character(len=*), parameter :: this_procedure = 'nek_advance_2Dh_axisym_nl'
            integer :: igeom
            common /cgeom/ igeom

            call init_2Dh_axisym()          ! idempotent; builds and registers the torsion terms
            if (ifsplit .or. ifmhd .or. ifpert .or. ifmvbd) call nek_stop_error(
     &         'only the standard P_N-P_{N-2} nonlinear step is supported.',
     &         this_module, this_procedure)
            if (param(103) > 0.0_dp) call nek_stop_error('filtering is not supported.',
     &         this_module, this_procedure)
            if (.not. if_torsion_zero .and. (param(54) /= 0.0_dp .or. param(55) /= 0.0_dp))
     &         call nek_stop_error('vol_flow (param 54/55) projects in-plane only: '//
     &         'not supported with torsion.', this_module, this_procedure)

            call nekgsync
            if (iftran) call settime
            call setsolv
            call comment
            call setprop
            call refresh_torsion_forcing_nl()            ! from U^n
            do igeom = 1, ngeom
               if (ifgeom) then
                  call gengeom(igeom)
                  call geneig (igeom)
               end if
               if (ifheat) call heat(igeom)
               if (igeom == 1) call add_torsion_pressure_phi_nl()   ! +D_lam^T p_ext
               if (ifflow) then
                  if (igeom == 2 .and. .not. if_torsion_zero) then
                     call fluid_2Dh_helix()
                  else
                     call fluid(igeom)
                  end if
               end if
            end do
         end subroutine nek_advance_2Dh_axisym_nl

         subroutine fluid_2Dh_helix()
      !! fluid(2) for the helix: Nek's fluid -> plan3(igeom=2)
      !! up to and including the in-plane predictor (sethlm, cresvif, ophinv,
      !! opadd2), then incomprn_2Dh_helix instead of incomprn.
      !! CHECK against fluid (drive2.f) and plan3 (plan3.f) of your Nek: this
      !! copies the v19 sequence.
            ifield = 1
            imesh  = 1
            call unorm
            call settolv
            call sethlm (nl_h1, nl_h2, -1)
            call cresvif(nl_resv(1,1), nl_resv(1,2), nl_resv(1,3), nl_h1, nl_h2)   ! lagvel, BCs, extrapp(pr), opgradt, bf, -H u
            call ophinv  (nl_dv(1,1), nl_dv(1,2), nl_dv(1,3),
     &                       nl_resv(1,1), nl_resv(1,2), nl_resv(1,3), nl_h1, nl_h2, tolhv, nmxv)
            call opadd2(vx, vy, vz, nl_dv(1,1), nl_dv(1,2), nl_dv(1,3))
            call incomprn_2Dh_helix()
         end subroutine fluid_2Dh_helix

         subroutine incomprn_2Dh_helix()
      !! Full projection of the nonlinear helix step (cf. incomprn, navier1/plan3):
      !!    E dp = -(D u* + D_lam U_phi*),   E = D_in (h2inv B^-1) D_in^T + D_lam (h2inv B^-1) D_lam^T,
      !!    u += h2inv B^-1 D_in^T dp,  U_phi += h2inv B^-1 D_lam^T dp,  p = p_ext + dp,
      !! with h2inv = dt/(bd(1) rho), as in the linear solver (setup_incomprp_2Dh).
      !! Nek's incomprn solves the same system scaled by bd(1)/dt with the same
      !! residual, so tolerances and iteration counts are comparable. U_phi* is
      !! from heat(2) of this step; pr holds p_ext (cresvif ran extrapp).
            real(dp) :: dtbd
            logical :: ifprjp
            integer :: istart, iter, ntot1, ntot2
            ntot1 = lx1*ly1*lz1*nelv
            ntot2 = lx2*ly2*lz2*nelv
            ifield = 1

            dtbd = bd(1)/dt
            if (dtbd /= dtbdlast_nl .or. istep <= 1) then
               nprev_nl = 0                                  ! operator changed: basis is stale
               dtbdlast_nl = dtbd
            end if
            call rzero  (nl_h1, ntot1)
            call cmult2 (nl_h2, vtrans(1,1,1,1,1), dtbd, ntot1)
            call invers2(nl_h2inv, nl_h2, ntot1)

      ! right-hand side -(D u* + D_lam U_phi*)
            call opdiv        (nl_dp, vx, vy, vz)
            call dphi_apply   (nl_w, t(1,1,1,1,1))
            call bm2_map12_2Dh(nl_w2, nl_w)
            call add2  (nl_dp, nl_w2, ntot2)
            call chsign(nl_dp, ntot2)
            call ortho (nl_dp)

            istart = param(95)
            ifprjp = (istep >= istart .and. istart /= 0)
            if (ifprjp) call setrhs_pressure_2Dh(nl_dp, nl_h1, nl_h2, nl_h2inv, pbasis_nl, nprev_nl,
     &                                           pres_matvec_helix_nl, 'NL')
            call solve_pressure_2Dh(nl_dp, nl_h1, nl_h2, nl_h2inv, pres_matvec_helix_nl, 1, iter,
     &                              0.0_dp, .true., maxit_p_axisym, 'lambda', torsion_lambda(), 'NL')
            if (ifprjp) call gensoln_pressure_2Dh(nl_dp, nl_h1, nl_h2, nl_h2inv, pbasis_nl, nprev_nl,
     &                                            pres_matvec_helix_nl)

      ! corrections: in-plane exactly as incomprn, u_phi with the same B^-1 (tmask)
            call opgradt(nl_resv(1,1), nl_resv(1,2), nl_resv(1,3), nl_dp)
            call opbinv (nl_dv(1,1), nl_dv(1,2), nl_dv(1,3), nl_resv(1,1), nl_resv(1,2), nl_resv(1,3), nl_h2inv)
            call opadd2 (vx, vy, vz, nl_dv(1,1), nl_dv(1,2), nl_dv(1,3))
            call lambda_correction(nl_w, nl_dp, nl_h2inv, tmask(1,1,1,1,1))
            call add2   (t(1,1,1,1,1), nl_w, ntot1)
            call add2   (pr, nl_dp, ntot2)
         end subroutine incomprn_2Dh_helix

         subroutine pres_matvec_helix_nl(Ap, p, h1, h2, h2inv, intype)
      !! Full pressure operator of the nonlinear helix step. The u_phi block uses
      !! the fluid h2inv: rhoCp = rho is enforced by check_2Dh_properties.
            real(dp), dimension(lp), intent(out) :: Ap
            real(dp), dimension(lp), intent(in)  :: p
            real(dp), dimension(lv), intent(in)  :: h1, h2, h2inv
            integer, intent(in) :: intype
            call cdabdtp(Ap, p, h1, h2, h2inv, intype)
            call add_lambda_block_E(Ap, p, h2inv, tmask(1,1,1,1,1))
         end subroutine pres_matvec_helix_nl

         subroutine seth2bc(h2)
      ! THIS SHOULD BE IN THE NEK CORE LATER BUT IS CURRENTLY NOT
      !
      !     sets diagonal terms for boundary conditions if required
      !
         include 'SIZE'
         include 'MASS'
         include 'PARALLEL'
         include 'SOLN'
         include 'ADJOINT'
         include 'TSTEP'
         include 'INPUT'
         include 'NEKUSE'
         include 'GEOM'
         
         real      h2(lx1,ly1,lz1,1)
         character cb*3
         integer   nfaces,nxyz,nel,ntot,iel,iface,ieg
         integer   kx1,kx2,ky1,ky2,kz1,kz2,ix,iy,iz,ia
         logical   ifalgn,ifnorx,ifnory,ifnorz
         real      unorm
         
         nfaces=2*ldim
         nxyz  =nx1*ny1*nz1
         nel   =nelfld(ifield)
         ntot  =nxyz*nel
         
         do iel=1,nelv
            ieg=lglel(iel)
            do iface=1,nfaces
               cb = cbc(iface,iel,ifield)
      !
      !     Adjoint O/o/ON/on conditions
      !
               if ((cb.eq.'O  '.or.cb.eq.'o  '.or.cb.eq.'ON '.or.cb.eq.'on ').and.ifadj) then
                  ia = 0
      !     Loop on the points of the face
                  call facind(kx1,kx2,ky1,ky2,kz1,kz2,nx1,ny1,nz1,iface)
                  do iz=kz1,kz2
                     do iy=ky1,ky2
                        do ix=kx1,kx2
                           ia = ia + 1
      !     add an U_n (baseflow velocity) coefficient to H2
                           h2(ix,iy,iz,iel)=h2(ix,iy,iz,iel)+
     &                     area(ia,1,iface,iel)/bm1(ix,iy,iz,iel)*
     &                     (unx(ia,1,iface,iel)*vx(ix,iy,iz,iel)
     &                     +uny(ia,1,iface,iel)*vy(ix,iy,iz,iel)
     &                     +unz(ia,1,iface,iel)*vz(ix,iy,iz,iel))
                        end do
                     end do
                  end do
               end if
      !
      !     bn condition
      !
               if (cb.eq.'bn ') then
                  ia = 0
                  call chknord (ifalgn,ifnorx,ifnory,ifnorz,iface,iel)
                  call facind(kx1,kx2,ky1,ky2,kz1,kz2,nx1,ny1,nz1,iface)
                  do iz=kz1,kz2
                     do iy=ky1,ky2
                        do ix=kx1,kx2
                           ia = ia + 1
      !     add an un (from userbc) coefficient to H2
                           if (optlevel.le.2) call nekasgn (ix,iy,iz,iel)
                           call userbc(ix,iy,iz,iface,ieg)
                           if(ifnorx) unorm=ux*unx(ia,1,iface,iel)
                           if(ifnory) unorm=uy*uny(ia,1,iface,iel)
                           if(ifnorz) unorm=uz*unz(ia,1,iface,iel)
                           h2(ix,iy,iz,iel)=h2(ix,iy,iz,iel)-param(1)*unorm
     &                       *area(ia,1,iface,iel)/bm1(ix,iy,iz,iel)
                        end do
                     end do
                  end do
               end if
            enddo
         enddo

         return
         end

      end module neklab_2Dh
