      module neklab_2Dh
      !! The 2Dh steppers.
      !!
      !!   nek_advance_2Dh_planar(betaz)   linear, cartesian, spanwise wavenumber
      !!   nek_advance_2Dh_axisym(alpha)   linear, axisymmetric (torus, helix)
      !!   nek_advance_2Dh_axisym_nl()     NONLINEAR, axisymmetric (torus, helix)
      !!
      !! The linear steppers are short sequences of stage calls named after Nek's
      !! perturbation solver (perturb.f). A loop over the perturbation slots jp
      !! is split only where a stage reads the other slot ipert = npert + 1 - jp.
      !!
      !! THE NONLINEAR STEP. Every nonlinear loop of the 2Dh axisymmetric drivers
      !! (steady map, orbit map, phase gauge, warmup) must call
      !! nek_advance_2Dh_axisym_nl, never nek_advance: it is the P_N-P_{N-2}
      !! branch of Nek's nek_advance (drive1.f) with three torsion hooks, and
      !! the linear stepper is its exact linearisation.
      !!   1. torsion_refresh_nl         G2 + G3 from U^n, read by userf/userq
      !!   2. torsion_add_pressure_nl    after heat(1): bq(u_s) += D_lam^T p_ext
      !!   3. fluid_2Dh_full             fluid(2) with ONE projection on the full
      !!                                 E = D B^-1 D^T, D = [D_in, D_lam],
      !!                                 correcting u_Z, u_R and u_s.
      !! Without torsion the hooks return at once and fluid(2) is Nek's: the
      !! step is nek_advance.
         use LightKrylov, only: dp
         use neklab_2Dh_common, only: setup_step_2Dh, cresvipp_2Dh, setup_incomprp_2Dh,
     &                                ebar_2Dh, incomprp_solve_2Dh, incomprp_update_2Dh,
     &                                solve_pressure_2Dh, setrhs_pressure_2Dh,
     &                                gensoln_pressure_2Dh
         use neklab_2Dh_planar, only: prepare_step_2Dh_planar, make_coupling_2Dh_planar,
     &                                makefp_2Dh_planar, ophx_2Dh_planar, ophinv_2Dh_planar,
     &                                incomprp_div_2Dh_planar, incomprp_dw_2Dh_planar,
     &                                pres_matvec_planar, maxit_p_planar
         use neklab_2Dh_axisym, only: init_2Dh_axisym,
     &                                prepare_step_2Dh_axisym, make_coupling_2Dh_axisym,
     &                                makefp_2Dh_axisym, ophx_2Dh_axisym,
     &                                ophinv_uz_2Dh_axisym, ophinv_urphi_2Dh_axisym,
     &                                incomprp_div_2Dh_axisym, incomprp_dw_2Dh_axisym,
     &                                incomprp_block_2Dh_axisym, alpha_eq_zero,
     &                                pres_matvec_axisym, maxit_p_axisym
         use neklab_2Dh_torsion, only: lambda_eq_zero, torsion_lambda,
     &                                 torsion_refresh_nl, torsion_add_pressure_nl,
     &                                 add_dlam_div, lambda_correction, add_lambda_block_E
         use neklab_nek_setup, only: nek_stop_error, neklab_cfl_scale
         use neklab_timing, only: neklab_timer_start, neklab_timer_stop, t_incomprn_full
         implicit none
         include "SIZE"
         include "TOTAL"
         include "ADJOINT" ! ifadj

         private
         character(len=*), parameter, private :: this_module = 'neklab_2Dh'
         integer, parameter :: lv = lx1*ly1*lz1*lelv
         integer, parameter :: lp = lx2*ly2*lz2*lelv

      ! Fields handed from one stage to the next (linear steppers)
         real(dp), dimension(lv,lpert),   private :: gradw    ! third-component pressure gradient
         real(dp), dimension(lv,3,lpert), private :: resv_jp  ! momentum residuals per slot
         real(dp), dimension(lv),         private :: dw       ! third-component pressure correction
         real(dp), dimension(lp,lpert),   private :: dpr      ! pressure correction

      ! Nonlinear step (fluid_2Dh_full, incomprn_2Dh_full)
         real(dp), dimension(lv,3), private :: nl_resv, nl_dv
         real(dp), dimension(lv),   private :: nl_h1, nl_h2, nl_h2inv, nl_ws
         real(dp), dimension(lp),   private :: nl_dp
         real(dp), private :: pbasis_nl(lp,mxprev)
         integer,  private :: nprev_nl = 0
         real(dp), private :: dtbdlast_nl = 0.0_dp
         real(dp), dimension(lp),   private :: dv_m2

         public :: nek_advance_2Dh_planar
         public :: nek_advance_2Dh_axisym
         public :: nek_advance_2Dh_axisym_nl
         public :: divergence_2Dh_axisym

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
            call prepare_step_2Dh_planar(betaz, if_advW_in)

            call setup_step_2Dh()

      ! --- explicit slot-coupled terms: read jp, write ipert
            do jp = 1, npert
               call make_coupling_2Dh_planar(gradw)
            end do

      ! --- per slot: right-hand side and velocity increments
            do jp = 1, npert
               igeom = 1
               call makefp_2Dh_planar()
               igeom = 2
               intype = -1
               ifield = 1
               call sethlm(h1, h2, intype)
               if (ifadj) call seth2bc(h2)
               call cresvipp_2Dh(resv1, resv2, resv3, gradw(1,jp))
               call ophx_2Dh_planar(resv1, resv2, resv3, h1, h2)

               call ophinv_2Dh_planar(dv1, dv2, dv3, resv1, resv2, resv3, h1, h2)
               call opadd2(vxp(1,jp), vyp(1,jp), vzp(1,jp), dv1, dv2, dv3)
               call add2  (tp(1,1,jp), dv3, ntot1)
            end do

      ! --- incompressibility
            call setup_incomprp_2Dh(h1, h2, h2inv)
            call ebar_2Dh(ebar, pres_matvec_planar, h1, h2, h2inv, betaz)
            do jp = 1, npert
               call incomprp_div_2Dh_planar(dpr(1,jp))
               call incomprp_solve_2Dh(dpr(1,jp), h1, h2, h2inv, pres_matvec_planar, ebar,
     &                                 betaz, maxit_p_planar, 'beta')
            end do

      ! --- correction
            do jp = 1, npert
               call incomprp_dw_2Dh_planar(dw, dpr, h2inv)
               call incomprp_update_2Dh(dpr(1,jp), dw, h2inv)
            end do

            jp = 0
         end subroutine nek_advance_2Dh_planar

         subroutine nek_advance_2Dh_axisym(alpha)
      !! One linear step of the axisymmetric 2Dh problem, azimuthal wavenumber
      !! alpha: the exact linearisation of nek_advance_2Dh_axisym_nl about the
      !! base flow in (vx, vy, t).
            real(dp), intent(in) :: alpha
            real(dp), dimension(lx1,ly1,lz1,lelv) :: h1, h2, h2inv
            common /SCRVH/  h1
            common /SCRCH/  h2
            common /scrhi/  h2inv
            real(dp) :: ebar
            integer :: igeom, intype

            call prepare_step_2Dh_axisym(alpha)

            call setup_step_2Dh()

      ! --- explicit slot-coupled terms: read jp, write ipert
            do jp = 1, npert
               call make_coupling_2Dh_axisym(gradw)       ! gradp (incl. D_lam^T p_ext), swirl advection
            end do

      ! --- momentum right-hand sides, per slot
            do jp = 1, npert
               igeom = 1
               call makefp_2Dh_axisym()                   ! incl. curvature and torsion
               igeom = 2
               intype = -1
               ifield = 1
               call sethlm(h1, h2, intype)
               if (ifadj) call seth2bc(h2)
               call cresvipp_2Dh(resv_jp(1,1,jp), resv_jp(1,2,jp), resv_jp(1,3,jp), gradw(1,jp))
            end do
            call ophx_2Dh_axisym(resv_jp, h1, h2)

      ! --- velocity increments: u_Z per slot, (u_R, u_phi) coupled at alpha /= 0
            call ophinv_uz_2Dh_axisym(resv_jp, h1, h2)
            call ophinv_urphi_2Dh_axisym(resv_jp, h1, h2)

      ! --- incompressibility (E includes the lambda block when torsion is on;
      !     at alpha /= 0 with torsion the two slots are solved together)
            call setup_incomprp_2Dh(h1, h2, h2inv)
            if (.not. (alpha_eq_zero .or. lambda_eq_zero)) then
               do jp = 1, npert
                  call incomprp_div_2Dh_axisym(dpr(1,jp))  ! -(D u + D_alpha u_s + D_lam u_s)
               end do
               call incomprp_block_2Dh_axisym(dpr, h1, h2, h2inv, alpha, maxit_p_axisym)
            else
               call ebar_2Dh(ebar, pres_matvec_axisym, h1, h2, h2inv, alpha)
               do jp = 1, npert
                  call incomprp_div_2Dh_axisym(dpr(1,jp))  ! -(D u + D_alpha u_phi + D_lam u_phi)
                  call incomprp_solve_2Dh(dpr(1,jp), h1, h2, h2inv, pres_matvec_axisym, ebar,
     &                                    alpha, maxit_p_axisym, 'alpha')
               end do
            end if

      ! --- correction
            do jp = 1, npert
               call incomprp_dw_2Dh_axisym(dw, dpr, h2inv)
               call incomprp_update_2Dh(dpr(1,jp), dw, h2inv)
            end do

            jp = 0
         end subroutine nek_advance_2Dh_axisym

      !=====================================================================
      !     NONLINEAR STEP
      !=====================================================================

         subroutine nek_advance_2Dh_axisym_nl()
      !! One NONLINEAR step of the axisymmetric 2Dh problem (torus, helix): the
      !! P_N-P_{N-2} branch of Nek's nek_advance (drive1.f), with the torsion
      !! hooks listed in the module header.
            character(len=*), parameter :: this_procedure = 'nek_advance_2Dh_axisym_nl'
            integer :: igeom
            common /cgeom/ igeom
            real(dp) :: p26, ct0

            call init_2Dh_axisym()
            if (ifsplit) call nek_stop_error('P_N-P_N is not supported by the 2Dh steppers.',
     &         this_module, this_procedure)
            if (ifmhd .or. ifpert .or. ifmvbd .or. ifcvode .or. ifrich .or. ifneknekc)
     &         call nek_stop_error('unsupported Nek option (mhd/pert/moving mesh/cvode/rich/neknek).',
     &         this_module, this_procedure)
            if (.not. lambda_eq_zero) then
               if (param(54) /= 0.0_dp .or. param(55) /= 0.0_dp) call nek_stop_error(
     &            'vol_flow (param 54/55) is in-plane only: not supported with torsion.',
     &            this_module, this_procedure)
               if (filterType /= 0) call nek_stop_error('filtering is not supported with torsion '//
     &            '(the linear step has no filter).', this_module, this_procedure)
            end if

            call nekgsync
            call setup_convect(2)                     ! convecting field from U^n (dealiasing)
            if (iftran) then
      ! variable dt: Nek's setdt sees only (vx, vy, vz); scale its target so that
      ! the full transport velocity (incl. the torsion rotation) meets param(26)
               p26 = param(26)
               ct0 = ctarg
               if (.not. lambda_eq_zero .and. param(12) > 0.0_dp) then
                  param(26) = p26*neklab_cfl_scale()
                  ctarg = param(26)
               end if
               call settime
               param(26) = p26
               ctarg = ct0
            end if
            call setsolv
            call comment
            call torsion_refresh_nl()                 ! G2 + G3 from U^n
            call setprop
            do igeom = 1, ngeom
               if (ifgeom) then
                  call gengeom(igeom)
                  call geneig (igeom)
               end if
               if (ifheat) call heat(igeom)
               if (igeom == 1) call torsion_add_pressure_nl()   ! bq(u_s) += D_lam^T p_ext
               if (ifflow) then
                  if (igeom == 2 .and. .not. lambda_eq_zero) then
                     call fluid_2Dh_full()
                  else
                     call fluid(igeom)
                  end if
               end if
               if (igeom == ngeom .and. filterType == 1) call q_filter(param(103))
            end do
         end subroutine nek_advance_2Dh_axisym_nl

         subroutine fluid_2Dh_full()
      !! fluid(2) of Nek (drive2.f) -> plan3(igeom = 2) (planx.f) up to and
      !! including the in-plane predictor, then incomprn_2Dh_full.
            integer :: intype
            ifield = 1
            imesh  = 1
            call unorm
            call settolv
            intype = -1
            call sethlm (nl_h1, nl_h2, intype)
            call cresvif(nl_resv(1,1), nl_resv(1,2), nl_resv(1,3), nl_h1, nl_h2)
            call ophinv (nl_dv(1,1), nl_dv(1,2), nl_dv(1,3),
     &                   nl_resv(1,1), nl_resv(1,2), nl_resv(1,3), nl_h1, nl_h2, tolhv, nmxv)
            call opadd2 (vx, vy, vz, nl_dv(1,1), nl_dv(1,2), nl_dv(1,3))
            call incomprn_2Dh_full()
            call chkptol
         end subroutine fluid_2Dh_full

         subroutine incomprn_2Dh_full()
      !! Projection with the full operator (cf. incomprn, induct.f):
      !!   E dp = -(D u* + D_lam U_s*),  E = D (h2inv B^-1) D^T + D_lam (h2inv B^-1) D_lam^T
      !!   u += h2inv B^-1 D^T dp,  U_s += h2inv B^-1 D_lam^T dp,  p = p_ext + dp
      !! with h2inv = dt/(bd(1) rho), the scaling of the linear solver. Nek's
      !! incomprn solves the same system scaled by bd(1)/dt with the same
      !! residual, so dp, tolerances and iteration counts are comparable. U_s* is
      !! from heat(2) of this step; pr holds p_ext (cresvif ran extrapp).
            real(dp) :: dtbd
            logical :: ifprjp
            integer :: istart, iter, ntot1, ntot2
            call neklab_timer_start(t_incomprn_full)
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

            call opdiv       (nl_dp, vx, vy, vz)
            call add_dlam_div(nl_dp, t(1,1,1,1,1))
            call chsign      (nl_dp, ntot2)
            call ortho       (nl_dp)

            istart = param(95)
            ifprjp = (istep >= istart .and. istart /= 0)
            if (ifprjp) call setrhs_pressure_2Dh(nl_dp, nl_h1, nl_h2, nl_h2inv, pbasis_nl, nprev_nl,
     &                                           pres_matvec_full_nl, 'NL')
            call solve_pressure_2Dh(nl_dp, nl_h1, nl_h2, nl_h2inv, pres_matvec_full_nl, 1, iter,
     &                              0.0_dp, .true., maxit_p_axisym, 'lambda', torsion_lambda(), 'NL')
            if (ifprjp) call gensoln_pressure_2Dh(nl_dp, nl_h1, nl_h2, nl_h2inv, pbasis_nl, nprev_nl,
     &                                            pres_matvec_full_nl)

            call opgradt(nl_resv(1,1), nl_resv(1,2), nl_resv(1,3), nl_dp)
            call opbinv (nl_dv(1,1), nl_dv(1,2), nl_dv(1,3),
     &                   nl_resv(1,1), nl_resv(1,2), nl_resv(1,3), nl_h2inv)
            call opadd2 (vx, vy, vz, nl_dv(1,1), nl_dv(1,2), nl_dv(1,3))
            call lambda_correction(nl_ws, nl_dp, nl_h2inv, tmask(1,1,1,1,1))
            call add2(t(1,1,1,1,1), nl_ws, ntot1)
            call add2(pr, nl_dp, ntot2)
            call neklab_timer_stop(t_incomprn_full)
         end subroutine incomprn_2Dh_full

         subroutine pres_matvec_full_nl(Ap, p, h1, h2, h2inv, intype)
      !! Full pressure operator of the nonlinear step (tmask = wmask is enforced
      !! by the torsion init).
            real(dp), dimension(lp), intent(out) :: Ap
            real(dp), dimension(lp), intent(in)  :: p
            real(dp), dimension(lv), intent(in)  :: h1, h2, h2inv
            integer, intent(in) :: intype
            call cdabdtp(Ap, p, h1, h2, h2inv, intype)
            call add_lambda_block_E(Ap, p, h2inv, tmask(1,1,1,1,1))
         end subroutine pres_matvec_full_nl

         real(dp) function divergence_2Dh_axisym(nonlinear, slot) result(rdiv)
      !! ||D u + D_lam u_s|| in the bm2^-1 norm of the pressure GMRES residual,
      !! for the nonlinear fields (nonlinear = .true.) or perturbation slot
      !! `slot`, alpha = 0. Should sit at the pressure tolerance after a step.
            logical, intent(in) :: nonlinear
            integer, optional, intent(in) :: slot
            real(dp), external :: glsc3
            integer :: ntot2, k
            ntot2 = lx2*ly2*lz2*nelv
            if (nonlinear) then
               call opdiv       (dv_m2, vx, vy, vz)
               call add_dlam_div(dv_m2, t(1,1,1,1,1))
            else
               k = 1
               if (present(slot)) k = slot
               call opdiv       (dv_m2, vxp(1,k), vyp(1,k), vzp(1,k))
               call add_dlam_div(dv_m2, tp(1,1,k))
            end if
            rdiv = sqrt(glsc3(dv_m2, dv_m2, bm2inv, ntot2)/volvm2)
         end function divergence_2Dh_axisym

         subroutine seth2bc(h2)
      ! Sets diagonal terms for adjoint O/ON and bn boundary conditions.
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
               if ((cb.eq.'O  '.or.cb.eq.'o  '.or.cb.eq.'ON '.or.cb.eq.'on ').and.ifadj) then
                  ia = 0
                  call facind(kx1,kx2,ky1,ky2,kz1,kz2,nx1,ny1,nz1,iface)
                  do iz=kz1,kz2
                     do iy=ky1,ky2
                        do ix=kx1,kx2
                           ia = ia + 1
                           h2(ix,iy,iz,iel)=h2(ix,iy,iz,iel)+
     &                     area(ia,1,iface,iel)/bm1(ix,iy,iz,iel)*
     &                     (unx(ia,1,iface,iel)*vx(ix,iy,iz,iel)
     &                     +uny(ia,1,iface,iel)*vy(ix,iy,iz,iel)
     &                     +unz(ia,1,iface,iel)*vz(ix,iy,iz,iel))
                        end do
                     end do
                  end do
               end if
               if (cb.eq.'bn ') then
                  ia = 0
                  call chknord (ifalgn,ifnorx,ifnory,ifnorz,iface,iel)
                  call facind(kx1,kx2,ky1,ky2,kz1,kz2,nx1,ny1,nz1,iface)
                  do iz=kz1,kz2
                     do iy=ky1,ky2
                        do ix=kx1,kx2
                           ia = ia + 1
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