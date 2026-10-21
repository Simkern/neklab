      module neklab_2Dh_planar
         use stdlib_optval, only: optval
         use LightKrylov, only: dp
         use neklab_nek_setup, only: nek_log_debug, nek_log_message, nek_log_information,
     &                               nek_log_warning, nek_stop_error
         use neklab_2Dh_common, only: wmask, is_initialized_2Dh_planar,
     &                                bm2_map12_2Dh, bm2_map12T_2Dh,
     &                                solve_helmholtz_2Dh, solve_pressure_2Dh,
     &                                setrhs_pressure_2Dh, gensoln_pressure_2Dh
         use neklab_timing, only: neklab_timer_start, neklab_timer_stop,
     &                            neklab_clock_tic, neklab_clock_toc,
     &                            t_planar_advance, t_planar_adv_setup, t_planar_adv_advection,
     &                            t_planar_adv_rhs, t_planar_adv_velocity, t_planar_adv_pressure,
     &                            t_planar_adv_correction, t_planar_makefp, t_planar_makeqp,
     &                            t_planar_solve_hmh, t_planar_solve_pres, t_planar_pres_precond,
     &                            t_planar_pres_proj, t_planar_hmh_matvec, t_planar_pres_matvec
         implicit none
         include "SIZE"
         include "TOTAL"
         include "ADJOINT" ! ifadj
         
         private
         character(len=*), parameter, private :: this_module = 'neklab_2Dh_planar'
         integer, parameter :: maxit_p = 100 ! max iterations for pressure solve
         
         real(dp), private :: pbasis(lx2*ly2*lz2*lelv,mxprev,lpert)
         real(dp), private :: dtbdlast = 0.0_dp
         real(dp), private :: betaz_op = 0.0_dp ! module variable to pass to the abstract operators
         integer,  private :: nprev(lpert) = 0

         ! Scratch arrays, at module level to keep them off the stack
         ! m1
         real(dp), dimension(lx1*ly1*lz1*lelv,lpert), private :: gradz_p               ! spanwise pressure gradient
         real(dp), dimension(lx1*ly1*lz1*lelv,lpert), private :: avdW_x,avdW_y,avdW_z  ! spanwise advection terms
         real(dp), dimension(lx1*ly1*lz1*lelv),       private :: dw                    ! spanwise velocity correction
         real(dp), dimension(lx1*ly1*lz1*lelv),       private :: advW_coef
         real(dp), dimension(lx1,ly1,lz1,lelv),       private :: hmv_tmp, pmv_tmp1     ! scratch
         ! m2
         real(dp), dimension(lx2*ly2*lz2*lelv,lpert), private :: dpr                   ! prefactor
         real(dp), dimension(lx2,ly2,lz2,lelv),       private :: frc_div               ! spanwise divergence
         real(dp), dimension(lx2*ly2*lz2*lelv),       private :: prextr, onep, ep      ! temp flds for p. recovery / extrap.
         real(dp), dimension(lx2,ly2,lz2,lelv),       private :: pmv_tmp2              ! scratch
         
         public :: nek_advance_2Dh_planar
      
      contains
      
         subroutine nek_advance_2Dh_planar(betaz, if_advW_in)
            implicit none
            real(dp), intent(in) :: betaz
            logical, optional, intent(in) :: if_advW_in
            ! spanwise wavenumber
            ! common
            real(dp), dimension(lx1,ly1,lz1,lelv) :: w1, w2, w3
            ! temporary velocity fields
            real(dp), dimension(lx1,ly1,lz1,lelv) :: resv1, resv2, resv3
            ! velocity residuals for Helmholtz solve
            real(dp), dimension(lx1,ly1,lz1,lelv) :: dv1, dv2, dv3
            ! velocity correction
            real(dp), dimension(lx1,ly1,lz1,lelv) :: h1, h2
            ! weights for the Helmholtz operator
            real(dp), dimension(lx1,ly1,lz1,lelv) :: h2inv
            ! utility array
            common /SCRNS/  w1, w2, w3, resv1, resv2, resv3, dv1, dv2, dv3
            common /SCRVH/  h1
            common /SCRCH/  h2
            common /scrhi/ h2inv
            real(dp) :: ebar
            ! mean pressure
            logical :: if_advW
            ! W advection active?
            ! Miscellaneous
            character(len=*), parameter :: this_procedure = 'nek_advance_2Dh_planar'
            character(len=*), parameter :: fmt = '(A,I6)'
            character(len=256) :: msg
            character(len=2) :: info_str
            character(len=7) :: hmh_info
            real(dp) :: dtbd
            logical :: ifprjp, if_ortho
            integer :: igeom, iter, intype, kfldfdm
            integer :: ntot1, ntot2, istart
            real, external :: glsum
            
            ntot1 = lx1*ly1*lz1*nelv
            ntot2 = lx2*ly2*lz2*nelv
            msg = 'call init_2Dh_planar() in userchk at istep = 0 before using the 2Dh solver.'
            if (.not. is_initialized_2Dh_planar()) call nek_stop_error(msg, this_module, this_procedure)
            betaz_op = betaz ! set the abstract operator's betaz
            if_ortho = (betaz == 0.0_dp)
            
            if_advW = optval(if_advW_in, .false.)
            msg = 'Spanwise advection '//merge('active.  ', 'inactive.', if_advW)
            call nek_log_debug(msg, this_module, this_procedure)
            
            write(msg,fmt) 'Step', istep
            call nek_log_debug(msg,this_module,this_procedure)
            
            if (istep == 1) nprev(:) = 0 ! reset pressure projection basis
            
            call neklab_timer_start(t_planar_advance)
            call neklab_clock_tic()
            
            call neklab_timer_start(t_planar_adv_setup)
            call nekgsync
            call settime
            call setup_convect(2)
            call setsolv
            call comment
            call setprop
            
            ifield = 1
            imesh = 1
            call unorm
            call settolv
            call neklab_timer_stop(t_planar_adv_setup)
            
            ! compute best estimate for d/dz (p)
            
            msg = 'Compute explicit pressure gradient term for w'
            call nek_log_debug(msg,this_module,this_procedure)
            
            call neklab_timer_start(t_planar_adv_advection)
            do jp = 1, npert
               ! we need to do this outside of the other jp loop for extrapolation of the correct pressure term
               call extrapprp(prextr)
               call compute_gradz_p(gradz_p, prextr, betaz)
               ! spanwise advection by the base flow: active iff betaz /= 0 and W /= 0
               if (if_advW) call compute_W_advection(avdW_x, avdW_y, avdW_z, betaz)
            end do
            call neklab_timer_stop(t_planar_adv_advection)
            
            ! Compute intermediate velocities
            
            msg = 'Solve momentum equations for u/v/w'
            call nek_log_debug(msg,this_module,this_procedure)
            
            do jp = 1, npert
               info_str = merge('Re', 'Im', jp == 1)
               call neklab_timer_start(t_planar_adv_rhs)
               igeom = 1
               !
               !  Old geometry, old velocity
               !
               ifield = 1     !  velocity
               call makefp_2Dh_planar(avdW_x, avdW_y, if_advW) ! W-advection added conditionally
               call lagfieldp
               ifield = 2     !  temperature (w)
               call makeqp_2Dh_planar(avdW_z, if_advW)        ! W-advection added conditionally
               call lagscalp
               !
               igeom = 2
               !
               !  New geometry, new velocity
               !
               intype = -1
               ifield = 1     !  velocity
               call sethlm   (h1,h2,intype)
               ! cresvipp
               call bcdirvc(vxp(1,jp), vyp(1,jp), vzp(1,jp),v1mask, v2mask, v3mask)
               call col2   (tp(1,1,jp), wmask, ntot1)
               call extrapprp(prextr)
               call opgradt(resv1,resv2,resv3,prextr)      ! d/dx pr, d/dy pr
               call opadd2 (resv1,resv2,resv3,bfxp(1,jp), bfyp(1,jp),bfzp(1,jp))
               call copy   (resv3,gradz_p(1,jp),ntot1)     ! d/dz pr
               call   add2 (resv3, bqp(1,1,jp), ntot1)
               ! ophx
               call helmholtz_matvec_2Dh(w1, vxp(1,jp),h1,h2,betaz)
               call helmholtz_matvec_2Dh(w2, vyp(1,jp),h1,h2,betaz)
               call helmholtz_matvec_2Dh(w3,tp(1,1,jp),h1,h2,betaz)
               call opsub2(resv1,resv2,resv3,w1,w2,w3)
               call   sub2(            resv3,      w3,ntot1)
               call neklab_timer_stop(t_planar_adv_rhs)
               
               call neklab_timer_start(t_planar_adv_velocity)
               ! ophinv
               tolhv = abs(tolhv)
               ifield = 1 !  velocity
               
               call dssum  (resv1,lx1,ly1,lz1)
               call col2   (resv1,v1mask,ntot1)
               hmh_info = info_str//' VELX'
               if (istep < 10) call chktcg1 (tolhv,resv1,h1,h2,v1mask,vmult,imesh,1)
               call solve_helmholtz_2Dh(dv1,resv1,h1,h2,v1mask,vmult,imesh,
     &                                  tolhv,nmxv,1,binvm1,hmh_info,
     &                                  hmh_matvec_planar,t_planar_solve_hmh)
               
               call dssum  (resv2,lx1,ly1,lz1)
               call col2   (resv2,v2mask,ntot1)
               hmh_info = info_str//' VELY'
               if (istep < 10) call chktcg1 (tolhv,resv2,h1,h2,v2mask,vmult,imesh,2)
               call solve_helmholtz_2Dh(dv2,resv2,h1,h2,v2mask,vmult,imesh,
     &                                  tolhv,nmxv,2,binvm1,hmh_info,
     &                                  hmh_matvec_planar,t_planar_solve_hmh)
               
               call dssum  (resv3,lx1,ly1,lz1)
               call col2   (resv3,wmask,ntot1)
               hmh_info = info_str//' VELZ'
               if (istep < 10) call chktcg1 (tolhv,resv3,h1,h2, wmask,vmult,imesh,3)
               call solve_helmholtz_2Dh(dv3,resv3,h1,h2,wmask,vmult,imesh,
     &                                  tolhv,nmxv,3,binvm1,hmh_info,
     &                                  hmh_matvec_planar,t_planar_solve_hmh)
               
               call opadd2 (vxp(1,jp),vyp(1,jp),vzp(1,jp),dv1,dv2,dv3)
               call add2   (tp(1,1,jp),dv3,ntot1)
               call neklab_timer_stop(t_planar_adv_velocity)
            end do ! jp
            
            msg = 'Compute pressure correction to enforce mass balance'
            call nek_log_debug(msg,this_module,this_procedure)
            
            call neklab_timer_start(t_planar_adv_pressure)
            ! prepare incomprp
            ifield = 1 ! velocity
            intype = 1
            dtbd   = bd(1)/dt
            if (dtbd /= dtbdlast) then
               nprev(:) = 0          ! operator changed: basis is stale
               dtbdlast = dtbd
            end if
            
            call rzero   (h1,ntot1)
            call cmult2  (h2,vtrans(1,1,1,1,ifield),dtbd,ntot1)
            call invers2 (h2inv,h2,ntot1)

            ! if beta /= 0, the E-operator is non-singular and the mean is NOT in the nullspace.
            ! But since the preconditioner automatically projects it out, we compute mean pressure
            ! to add it back after the pressure solve
            ebar = 0.0_dp
            if (.not. if_ortho .and. ifvcor) then
               call rone(onep, ntot2)
               call pressure_matvec_2Dh(ep,onep,h1,h2,h2inv,betaz,intype)
               ebar = glsum(ep, ntot2)
            end if

            ! project out previous pressure solutions?
            ifprjp=.false.
            istart=param(95)
            if (istep.ge.istart.and.istart.ne.0) ifprjp=.true.
            
            do jp = 1, npert
               info_str = merge('Re', 'Im', jp == 1)
               
               call opdiv   (dpr(1,jp),vxp(1,jp),vyp(1,jp),vzp(1,jp))
               
               call compute_frc_div(frc_div, tp, betaz)
               call add2    (dpr(1,jp),frc_div,ntot2)
               
               call chsign  (dpr(1,jp),ntot2)
               if (if_ortho) call ortho (dpr(1,jp))
               
               if (ifprjp) then
                  call  setrhs_pressure_2Dh(dpr(1,jp), h1, h2, h2inv, pbasis(1,1,jp), nprev(jp),
     &                                      pres_matvec_planar, t_planar_pres_proj, info_str)
               end if
               call      solve_pressure_2Dh(dpr(1,jp), h1, h2, h2inv, pres_matvec_planar, intype, iter,
     &                                      ebar, if_ortho, maxit_p, 'beta', betaz, info_str,
     &                                      t_planar_solve_pres, t_planar_pres_precond)
               if (ifprjp) then
                  call gensoln_pressure_2Dh(dpr(1,jp), h1, h2, h2inv, pbasis(1,1,jp), nprev(jp),
     &                                      pres_matvec_planar, t_planar_pres_proj)
               end if
            end do ! jp
            
            call neklab_timer_stop(t_planar_adv_pressure)
            
            ! Reconstruct pressure and add pressure correction
            
            msg = 'Compute velocity correction based on pressure correction'
            call nek_log_debug(msg,this_module,this_procedure)
            
            call neklab_timer_start(t_planar_adv_correction)
            do jp = 1, npert
               call opgradt (w1 ,w2 ,w3 ,dpr(1,jp))
               call opbinv  (dv1,dv2,dv3,w1 ,w2 ,w3 ,h2inv)
               call compute_dw(dw, dpr, h2inv, betaz) ! dw
               
               call opadd2 (vxp(1,jp),vyp(1,jp),vzp(1,jp), dv1,dv2,dv3)
               call   add2 (tp(1,1,jp), dw, ntot1)
               
               call extrapprp(prextr)
               call lagpresp
               call add3(prp(1,jp),prextr,dpr(1,jp), ntot2)
            end do ! jp
            
            call neklab_timer_stop(t_planar_adv_correction)
            
            call neklab_clock_toc()
            call neklab_timer_stop(t_planar_advance)
            
            ! reset jp = 0 in case we run the nonlinear step next
            jp = 0
         
         end subroutine nek_advance_2Dh_planar
      
         subroutine compute_dw(dw, dpr, h2inv, betaz)
            implicit none
            include 'SIZE'
            include 'SOLN' ! jp
            include 'MASS' ! bm1
            real, dimension(lx1*ly1*lz1*lelv), intent(out) :: dw
            real, dimension(lx2*ly2*lz2*lelv,lpert), intent(in) :: dpr
            real, dimension(lx1,ly1,lz1,lelv), intent(in) :: h2inv
            real, intent(in) :: betaz
            ! internal
            integer :: ipert, spert, ntot1
            real, dimension(lx1*ly1*lz1) :: wrk1, wrk2
            
            ipert = npert + 1 - jp
            spert = merge(1, -1, jp == 1)  ! +1 for jp = 1, -1 for jp = 2
            ntot1 = lx1*ly1*lz1*nelv
            call rzero  (dw, ntot1)
            call mappr  (dw, dpr(1,ipert), wrk1, wrk2)   ! map to vmesh
            call cmult  (dw, -spert*betaz, ntot1)       ! collate coupling
            !
            call col2   (dw, bm1, ntot1)                 ! collate mass matrix
            call col2   (dw, wmask, ntot1)               ! mask Dirichlet conditions
            call dssum  (dw, lx1, ly1, lz1)              ! make continuous
            call col2   (dw, binvm1, ntot1)              ! collate inverse mass matrix
            !
            call col2   (dw, h2inv, ntot1)               ! collate inverse h2
         end
      
         subroutine compute_frc_div(frc_div, w, betaz)
            implicit none
            include 'SIZE'
            include 'SOLN' ! jp
            include 'MASS' ! bm2
            real, dimension(lx2,ly2,lz2,lelv), intent(out) :: frc_div
            real, dimension(lx1*ly1*lz1*lelv,ldimt,lpert), intent(in) :: w
            real, intent(in) :: betaz
            ! internal
            integer :: ipert, spert, ie, ie1, nxyz1, ntot2
            
            ipert = npert + 1 - jp
            spert = merge(1, -1, jp == 1)  ! +1 for jp = 1, -1 for jp = 2
            nxyz1 = lx1*ly1*lz1
            ntot2 = lx2*ly2*lz2*nelv
            call rzero(frc_div, ntot2)
            do ie = 1, nelv
               ie1 = (ie-1)*nxyz1+1
               call map12 (frc_div(1,1,1,ie), w(ie1,1,ipert), ie) ! map to pmesh
            end do
            call col2c(frc_div, bm2, spert*betaz, ntot2) ! this is to be consistent with opdiv output
         end
      
         subroutine compute_gradz_p(gradz_p, prextr, betaz)
            implicit none
            include 'SIZE'
            include 'SOLN' ! jp
            include 'MASS' ! bm1
            real, dimension(lx1*ly1*lz1*lelv,lpert), intent(inout) :: gradz_p
            real, dimension(lx2*ly2*lz2*lelv), intent(in) :: prextr
            real, intent(in) :: betaz
            ! internal
            integer :: ipert, spert, ntot1
            real, dimension(lx1*ly1*lz1) :: wrk1, wrk2
            
            ipert = npert + 1 - jp
            spert = merge(1, -1, jp == 1)  ! +1 for jp = 1, -1 for jp = 2
            ntot1 = lx1*ly1*lz1*nelv
            call rzero(gradz_p(1,ipert), ntot1)
            call mappr(gradz_p(1,ipert), prextr, wrk1, wrk2)  ! map to vmesh
            call col2c(gradz_p(1,ipert), bm1, spert*betaz, ntot1)
            ! dssum will be called on the vel residual later
         end
      
         subroutine compute_W_advection(avdW_x, avdW_y, avdW_z, betaz)
            ! + i * betaz * W*v[xy]p: read/write switches jp.
            implicit none
            include 'SIZE'
            include 'SOLN'   ! jp, vxp, vyp, tp
            include 'MASS'   ! bm1
            real(dp), dimension(lx1*ly1*lz1*lelv,lpert), intent(inout) :: avdW_x, avdW_y, avdW_z
            real(dp), intent(in) :: betaz
            ! internal
            integer :: ipert, spert, ntot1
            
            ntot1 = lx1*ly1*lz1*nelv
            ipert = npert + 1 - jp
            spert = merge(1, -1, jp == 1)  ! +1 for jp = 1, -1 for jp = 2
            
            call col3 (advW_coef, t(1,1,1,1,1), bm1, ntot1)         ! W*bm1
            call cmult(advW_coef, spert*betaz, ntot1)              ! +-beta*W*bm1
            
            call col3 (avdW_x(1,ipert), advW_coef, vxp(1,jp),  ntot1)
            call col2 (avdW_x(1,ipert), vtrans(1,1,1,1,1), ntot1)   ! rho
            call col3 (avdW_y(1,ipert), advW_coef, vyp(1,jp),  ntot1)
            call col2 (avdW_y(1,ipert), vtrans(1,1,1,1,1), ntot1)   ! rho
            call col3 (avdW_z(1,ipert), advW_coef, tp(1,1,jp), ntot1)
            call col2 (avdW_z(1,ipert), vtrans(1,1,1,1,2), ntot1)   ! rho*cp
         end subroutine compute_W_advection
      
      
         subroutine makefp_2Dh_planar(avdW_x, avdW_y, if_advW)
            ! Copy of the core makefp (perturb.f)
            ! Spanwise advection of vxp and vyp injected before makextp if W /= 0.
            implicit none
            include 'SIZE'
            include 'SOLN'
            include 'INPUT'
            include 'TSTEP'
            include 'ADJOINT'
            real(dp), dimension(lx1*ly1*lz1*lelv,lpert), intent(in) :: avdW_x, avdW_y
            logical, intent(in) :: if_advW
            integer :: ntot1
            
            call neklab_timer_start(t_planar_makefp)
            ntot1 = lx1*ly1*lz1*nelv
            ! user-specified forcing
                                                               call makeufp
            ! baseflow advection
            if (ifnav .and. (.not.ifchar) .and. (.not.ifadj))  call advabp
            if (ifnav .and. (.not.ifchar) .and. (     ifadj))  call advabp_adjoint
            !! W addition
            if (if_advW) then
               ! advection via the W component
               call add2(bfxp(1,jp), avdW_x(1,jp), ntot1)
               call add2(bfyp(1,jp), avdW_y(1,jp), ntot1)
            end if
            ! add extrapolation terms
            if (iftran)                                        call makextp
            ! add bd contribution
                                                               call makebdfp
            call neklab_timer_stop(t_planar_makefp)
         end subroutine makefp_2Dh_planar
         
         subroutine makeqp_2Dh_planar(avdW_z, if_advW)
            ! Copy of the core makeqp (perturb.f)
            ! Spanwise advection of vzp injected before makeabqp if W /= 0.
            ! Note: Expects ifield = 2.
            implicit none
            include 'SIZE'
            include 'SOLN'
            include 'INPUT'
            include 'TSTEP'
            real(dp), dimension(lx1*ly1*lz1*lelv,lpert), intent(in) :: avdW_z
            logical, intent(in) :: if_advW
            integer :: ntot1
            
            call neklab_timer_start(t_planar_makeqp)
            ntot1 = lx1*ly1*lz1*nelv
            ! user-specified forcing
                                                                  call makeuqp
            ! baseflow advection
            if (ifadvc(ifield) .and. (.not.ifchar))               call convabp
            !! W additions
            if (if_advW) then
               ! advection via the W component
               call add2(bqp(1,ifield-1,jp), avdW_z(1,jp), ntot1)
            end if
            ! add extrapolation terms
            if (iftran)                                           call makeabqp
            if ((iftran .and. .not.ifchar) .or.
     &          (iftran .and. .not.ifadvc(ifield) .and. ifchar))  call makebdqp
            call neklab_timer_stop(t_planar_makeqp)
         end subroutine makeqp_2Dh_planar
         
         subroutine pressure_matvec_2Dh(ap,wp,h1,h2,h2inv,betaz,intype)
            implicit none
            
         ! INTYPE= 0  Compute the matrix-vector product    DA(-1)DT*p
         ! INTYPE= 1  Compute the matrix-vector product    D(B/DT)(-1)DT*p
         ! INTYPE=-1  Compute the matrix-vector product    D(A+B/DT)(-1)DT*p
            
            real, dimension(lx2,ly2,lz2,lelv), intent(out) :: Ap
            real, dimension(lx2,ly2,lz2,lelv), intent(in) :: wp
            real, dimension(lx1,ly1,lz1,lelv), intent(in) :: h1
            real, dimension(lx1,ly1,lz1,lelv), intent(in) :: h2
            real, dimension(lx1,ly1,lz1,lelv), intent(in) :: h2inv
            !!! add the contribution from the 3rd perturbation component
            real, intent(in) :: betaz
            integer, intent(in) :: intype
            ! internal
            real, dimension(lx1*ly1*lz1) :: wrk1, wrk2
            
            integer :: ie, ntot1, ntot2
            
            call neklab_timer_start(t_planar_pres_matvec)
            ! Regular Ax = D (h2*B)^(-1) D^T
            call cdabdtp(Ap, wp, h1, h2, h2inv, intype)
            ! add 3rd perturbation component.
            ! we want to subtract the contribution
            !
            ! B * d/dz (h2*B)^(-1) d/dz wp = (-i betaz) (h2)^(-1) (-i betaz) wp = - betaz^2 (h2)^(-1) wp
            !
            ! Note: we are only considering the ifanls = .false. case
            ntot1 = lx1*ly1*lz1*nelv
            ntot2 = lx2*ly2*lz2*nelv
            call mappr  (pmv_tmp1, wp, wrk1, wrk2)                    ! map to vmesh
            !
            call col2   (pmv_tmp1, bm1,   ntot1)                      ! collate mass matrix
            call col2   (pmv_tmp1, wmask, ntot1)                      ! mask Dirichlet conditions
            call dssum  (pmv_tmp1, lx1, ly1, lz1)                     ! make continuous
            call col2   (pmv_tmp1, binvm1, ntot1)                     ! collate inverse mass matrix
            !
            call col2   (pmv_tmp1, h2inv, ntot1)                      ! collate (h2)^-1
            call cmult  (pmv_tmp1, -betaz**2, ntot1)                 ! collate -betaz^2
            do ie = 1, nelv
               call map12 (pmv_tmp2(1,1,1,ie), pmv_tmp1(1,1,1,ie), ie)! map wdiv to pmesh
            end do
            call col2(pmv_tmp2, bm2, ntot2)                           ! collate mass matrix (pressure mesh)
            call sub2(ap, pmv_tmp2, ntot2)
            call neklab_timer_stop(t_planar_pres_matvec)
         end subroutine pressure_matvec_2Dh
         
         subroutine helmholtz_matvec_2Dh(Au, u, h1, h2, betaz)
            implicit none
            real, dimension(lx1,ly1,lz1,1), intent(out) :: Au
            real, dimension(lx1,ly1,lz1,1), intent(in) :: u
            real, dimension(lx1,ly1,lz1,1), intent(in) :: h1
            real, dimension(lx1,ly1,lz1,1), intent(in) :: h2
            real betaz
            ! internal
            integer :: imesh, ntot1
            call neklab_timer_start(t_planar_hmh_matvec)
            imesh = 1
            ntot1 = lx1*ly1*lz1*nelv
            ! regular Ax
            call axhelm (Au, u, h1, h2, imesh, 1)
            ! added spanwise viscous diffusion term
            call col3 (hmv_tmp, u, vdiff(1,1,1,1,1), ntot1)
            call col2c(hmv_tmp, bm1, betaz**2, ntot1)
            call add2 (Au, hmv_tmp, ntot1)
            call neklab_timer_stop(t_planar_hmh_matvec)
         end subroutine helmholtz_matvec_2Dh

         subroutine hmh_matvec_planar(Au, u, h1, h2)
            ! Wrapper around helmholtz_matvec_2Dh to match the abstract interface while passing betaz
            real(dp), dimension(lx1*ly1*lz1*lelv), intent(out) :: Au
            real(dp), dimension(lx1*ly1*lz1*lelv), intent(in)  :: u, h1, h2
            call helmholtz_matvec_2Dh(Au, u, h1, h2, betaz_op)
         end subroutine hmh_matvec_planar

         subroutine pres_matvec_planar(Ap, p, h1, h2, h2inv, intype)
            ! Wrapper around the pressure_matvec_2Dh to match the abstract interface while passing betaz
            real(dp), dimension(lx2*ly2*lz2*lelv), intent(out) :: Ap
            real(dp), dimension(lx2*ly2*lz2*lelv), intent(in)  :: p
            real(dp), dimension(lx1*ly1*lz1*lelv), intent(in)  :: h1, h2, h2inv
            integer, intent(in) :: intype
            call pressure_matvec_2Dh(Ap, p, h1, h2, h2inv, betaz_op, intype)
         end subroutine pres_matvec_planar
      
      end module neklab_2Dh_planar
