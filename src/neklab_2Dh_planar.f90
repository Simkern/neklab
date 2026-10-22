      module neklab_2Dh_planar
         use stdlib_optval, only: optval
         use LightKrylov, only: dp
         use neklab_nek_setup, only: nek_log_debug, nek_log_message, nek_log_information,
     &                               nek_log_warning, nek_stop_error
         use neklab_2Dh_common, only: wmask, build_wmask, 
     &                                check_2Dh_properties, check_2Dh_masks,
     &                                bm2_map12_2Dh, bm2_map12T_2Dh,
     &                                solve_helmholtz_2Dh
         use neklab_timing, only: neklab_timer_start, neklab_timer_stop,
     &                            t_planar_makefp, t_planar_makeqp,
     &                            t_planar_helm_matvec, t_planar_pres_matvec
         implicit none
         include "SIZE"
         include "TOTAL"
         include "ADJOINT" ! ifadj
         
         private
         character(len=*), parameter, private :: this_module = 'neklab_2Dh_planar'
         integer, parameter :: lv = lx1*ly1*lz1*lelv
         integer, parameter :: lp = lx2*ly2*lz2*lelv
         integer, parameter, public :: maxit_p_planar = 100 ! max iterations for pressure solve
         
         real(dp), private :: betaz_op = 0.0_dp                ! betaz of the current step (prepare_step_2Dh_planar)
         logical,  private :: if_advW_op = .false.             ! W-advection switch of the current step
         logical,  private :: init_2Dh_planar_done = .false.   ! initialization flag
         logical,  public  :: betaz_eq_zero = .false.          ! flag for beta = 0 (no Re/Im pairing)

      ! Fields produced and consumed inside the model
         real(dp), dimension(lv,lpert), private :: advW_x, advW_y, advW_z  ! W advection (make_coupling -> makefp)
      ! Scratch, module level to keep it off the stack. Holds nothing between stage calls.
         real(dp), dimension(lv), private :: tmp_m1a ! ophx (H u), W-advection coefficient, pressure matvec
         real(dp), dimension(lv), private :: tmp_m1b ! helmholtz_matvec_2Dh_planar (ophx passes tmp_m1a as its output)
         real(dp), dimension(lp), private :: tmp_m2  ! extrapolated pressure (make_coupling), divergence (incomprp_div),
                                                     ! pressure matvec (runs only after both are consumed)
         
         public ::           init_2Dh_planar
         public :: is_initialized_2Dh_planar
         public ::   prepare_step_2Dh_planar
         public ::  make_coupling_2Dh_planar
         public ::         makefp_2Dh_planar
         public ::           ophx_2Dh_planar
         public ::         ophinv_2Dh_planar
         public ::   incomprp_div_2Dh_planar
         public ::    incomprp_dw_2Dh_planar
         public ::        pres_matvec_planar
      
      contains

      !-----------------------------------------------------------!
      !                                                           !
      !     PUBLIC STAGE WRAPPERS FOR THE PLANAR 2Dh SOLVER       !
      !                                                           !
      !-----------------------------------------------------------!

         subroutine init_2Dh_planar(verbose)
      !! Set up the cartesian 2Dh solver. Idempotent. Called by
      !! init_exptA_2Dh_planar and, as a fallback, by prepare_step_2Dh_planar;
      !! may also be called by hand from userchk (istep = 0, after usrdat2).
            logical, optional, intent(in) :: verbose
            ! internal
            character(len=*), parameter :: this_procedure = 'init_2Dh_planar'

            if (init_2Dh_planar_done) return

            if (ifaxis) call nek_stop_error(
     &         'init_2Dh_planar is for the cartesian 2Dh solver, incompatible with ifaxis.',
     &         this_module, this_procedure)
            call check_2Dh_properties(this_procedure)

            call build_wmask(verbose)
            call check_2Dh_masks(this_procedure)

            init_2Dh_planar_done = .true.
            call nek_log_message('cartesian 2Dh solver initialised.', this_module, this_procedure)
         end subroutine init_2Dh_planar

         logical function is_initialized_2Dh_planar()
            is_initialized_2Dh_planar = init_2Dh_planar_done
         end function is_initialized_2Dh_planar 

         subroutine prepare_step_2Dh_planar(betaz, if_advW_in)
      !! Checks and per-step switches. Called first by nek_advance_2Dh_planar.
            real(dp), intent(in) :: betaz
            logical, optional, intent(in) :: if_advW_in
            ! internal
            character(len=*), parameter :: this_procedure = 'prepare_step_2Dh_planar'
            character(len=256) :: msg
            call init_2Dh_planar() ! for safety; normally done already
            betaz_op   = betaz
            if_advW_op = optval(if_advW_in, .false.)
            msg = 'Spanwise advection '//merge('active.  ', 'inactive.', if_advW_op)
            call nek_log_debug(msg, this_module, this_procedure)
            write(msg,'(A,I6)') 'Step', istep
            call nek_log_debug(msg, this_module, this_procedure)
            ! Re/Im pairing is only meaningful for betaz /= 0; npert must match
            ! or the caller has set up the wrong number of perturbation slots.
            betaz_eq_zero = (betaz == 0.0_dp)
            if (betaz_eq_zero) then
               if (lpert /= 1) call nek_stop_error('betaz = 0 requires npert = 1.', this_module, this_procedure)
            else
               if (lpert /= 2) call nek_stop_error('betaz /= 0 requires npert = 2.', this_module, this_procedure)
            end if
         end subroutine prepare_step_2Dh_planar

         subroutine make_coupling_2Dh_planar(gradw)
      !! Slot-coupled explicit terms of slot jp (read jp, write ipert): the
      !! spanwise pressure gradient and, if active, the W-advection. Runs for all
      !! slots before any slot is advanced.
            real(dp), dimension(lv,lpert), intent(inout) :: gradw
            call extrapprp(tmp_m2)
            call compute_gradp_planar(gradw, tmp_m2, betaz_op)
            if (if_advW_op) call compute_W_advection_planar(advW_x, advW_y, advW_z, betaz_op)
         end subroutine make_coupling_2Dh_planar

         subroutine makefp_2Dh_planar()
      !! Old-geometry part of slot jp (cf. perturbv, igeom = 1): makefp + lagfieldp
      !! for u, v and makeqp + lagscalp for w, with the W-advection injected.
            ifield = 1
            call makefp_uv_2Dh_planar(advW_x, advW_y, if_advW_op)
            call lagfieldp
            ifield = 2
            call makeqp_w_2Dh_planar(advW_z, if_advW_op)
            call lagscalp
         end subroutine makefp_2Dh_planar

         subroutine ophx_2Dh_planar(resv1, resv2, resv3, h1, h2)
      !! resv -= H u for the three components of slot jp (cf. ophx in cresvipp).
            real(dp), dimension(lv), intent(inout) :: resv1, resv2, resv3
            real(dp), dimension(lv), intent(in)    :: h1, h2
            ! internal
            integer :: ntot1
            ntot1 = lx1*ly1*lz1*nelv
            call helmholtz_matvec_2Dh_planar(tmp_m1a, vxp(1,jp),  h1, h2, betaz_op)
            call sub2(resv1, tmp_m1a, ntot1)
            call helmholtz_matvec_2Dh_planar(tmp_m1a, vyp(1,jp),  h1, h2, betaz_op)
            call sub2(resv2, tmp_m1a, ntot1)
            call helmholtz_matvec_2Dh_planar(tmp_m1a, tp(1,1,jp), h1, h2, betaz_op)
            call sub2(resv3, tmp_m1a, ntot1)
         end subroutine ophx_2Dh_planar

         subroutine ophinv_2Dh_planar(dv1, dv2, dv3, resv1, resv2, resv3, h1, h2)
      !! Velocity increments of slot jp (cf. ophinv): assemble, mask, solve.
            real(dp), dimension(lv), intent(out)   :: dv1, dv2, dv3
            real(dp), dimension(lv), intent(inout) :: resv1, resv2, resv3
            real(dp), dimension(lv), intent(in)    :: h1, h2
            ! internal
            character(len=2) :: info_str
            character(len=7) :: helm_info
            integer :: ntot1
            ntot1    = lx1*ly1*lz1*nelv
            info_str = merge('Re', 'Im', jp == 1)
            tolhv    = abs(tolhv)
            ifield   = 1
            
            call dssum  (resv1,lx1,ly1,lz1)
            call col2   (resv1,v1mask,ntot1)
            helm_info = info_str//' VELX'
            if (istep < 10) call chktcg1 (tolhv,resv1,h1,h2,v1mask,vmult,imesh,1)
            call solve_helmholtz_2Dh(dv1,resv1,h1,h2,v1mask,vmult,imesh,
     &                               tolhv,nmxv,1,binvm1,helm_info,helm_matvec_planar)
            
            call dssum  (resv2,lx1,ly1,lz1)
            call col2   (resv2,v2mask,ntot1)
            helm_info = info_str//' VELY'
            if (istep < 10) call chktcg1 (tolhv,resv2,h1,h2,v2mask,vmult,imesh,2)
            call solve_helmholtz_2Dh(dv2,resv2,h1,h2,v2mask,vmult,imesh,
     &                               tolhv,nmxv,2,binvm1,helm_info,helm_matvec_planar)
            
            call dssum  (resv3,lx1,ly1,lz1)
            call col2   (resv3,wmask,ntot1)
            helm_info = info_str//' VELZ'
            if (istep < 10) call chktcg1 (tolhv,resv3,h1,h2,wmask,vmult,imesh,3)
            call solve_helmholtz_2Dh(dv3,resv3,h1,h2,wmask,vmult,imesh,
     &                               tolhv,nmxv,3,binvm1,helm_info,helm_matvec_planar)
         end subroutine ophinv_2Dh_planar

         subroutine incomprp_div_2Dh_planar(dpr_jp)
      !! Right-hand side of the pressure step of slot jp: -(D u + D_beta w). The
      !! beta part reads the other slot's w.
            real(dp), dimension(lp), intent(out) :: dpr_jp
            ! internal
            integer :: ntot2
            ntot2 = lx2*ly2*lz2*nelv
            call opdiv  (dpr_jp, vxp(1,jp), vyp(1,jp), vzp(1,jp))
            call compute_frc_div_planar(tmp_m2, tp, betaz_op)
            call add2   (dpr_jp, tmp_m2, ntot2)
            call chsign (dpr_jp, ntot2)
         end subroutine incomprp_div_2Dh_planar

         subroutine incomprp_dw_2Dh_planar(dw, dpr, h2inv)
      !! w correction of slot jp from the other slot's pressure correction.
            real(dp), dimension(lv), intent(out) :: dw
            real(dp), dimension(lp,lpert), intent(in) :: dpr
            real(dp), dimension(lv), intent(in) :: h2inv
            call compute_dw_planar(dw, dpr, h2inv, betaz_op)
         end subroutine incomprp_dw_2Dh_planar

      !-----------------------------------!
      !                                   !
      !     PRIVATE HELPER ROUTINES       !
      !                                   !
      !-----------------------------------!
      
         subroutine compute_dw_planar(dw, dpr, h2inv, betaz)
            implicit none
            include 'SIZE'
            include 'SOLN' ! jp
            include 'MASS' ! bm1
            real, dimension(lv), intent(out) :: dw
            real, dimension(lp,lpert), intent(in) :: dpr
            real, dimension(lx1,ly1,lz1,lelv), intent(in) :: h2inv
            real, intent(in) :: betaz
            ! internal
            integer :: ipert, spert, ntot1
            
            ipert = npert + 1 - jp
            spert = merge(1, -1, jp == 1)  ! +1 for jp = 1, -1 for jp = 2
            ntot1 = lx1*ly1*lz1*nelv
            call bm2_map12t_2Dh(dw, dpr(1,ipert))            ! I12^T (bm2 dp)
            call cmult         (dw, -spert*betaz, ntot1)     ! collate coupling coefficient
            call col2          (dw, wmask, ntot1)            ! mask Dirichlet conditions
            call dssum         (dw, lx1, ly1, lz1)           ! make continuous
            call col2          (dw, binvm1, ntot1)           ! collate inverse mass matrix
            call col2          (dw, h2inv, ntot1)            ! collate inverse h2
         end
      
         subroutine compute_frc_div_planar(frc_div, w, betaz)
            implicit none
            include 'SIZE'
            include 'SOLN' ! jp
            include 'MASS' ! bm2
            real, dimension(lx2,ly2,lz2,lelv), intent(out) :: frc_div
            real, dimension(lv,ldimt,lpert), intent(in) :: w
            real, intent(in) :: betaz
            ! internal
            integer :: ipert, spert, ntot2
            
            ipert = npert + 1 - jp
            spert = merge(1, -1, jp == 1)  ! +1 for jp = 1, -1 for jp = 2
            ntot2 = lx2*ly2*lz2*nelv
            call bm2_map12_2Dh(frc_div, w(1,1,ipert))       ! bm2 * I12 w
            call cmult        (frc_div, spert*betaz, ntot2) ! this is to be consistent with opdiv output
         end
      
         subroutine compute_gradp_planar(gradz_p, prextr, betaz)
            implicit none
            include 'SIZE'
            include 'SOLN' ! jp
            include 'MASS' ! bm1
            real, dimension(lv,lpert), intent(inout) :: gradz_p
            real, dimension(lp), intent(in) :: prextr
            real, intent(in) :: betaz
            ! internal
            integer :: ipert, spert, ntot1
            
            ipert = npert + 1 - jp
            spert = merge(1, -1, jp == 1)  ! +1 for jp = 1, -1 for jp = 2
            ntot1 = lx1*ly1*lz1*nelv
            call bm2_map12t_2Dh(gradz_p(1,ipert), prextr)             ! I12^T (bm2 p)
            call cmult         (gradz_p(1,ipert), spert*betaz, ntot1)
            ! dssum will be called on the vel residual later
         end
      
         subroutine compute_W_advection_planar(advW_x, advW_y, advW_z, betaz)
            ! + i * betaz * W*v[xy]p: read/write switches jp.
            implicit none
            include 'SIZE'
            include 'SOLN'   ! jp, vxp, vyp, tp
            include 'MASS'   ! bm1
            real(dp), dimension(lv,lpert), intent(inout) :: advW_x, advW_y, advW_z
            real(dp), intent(in) :: betaz
            ! internal
            integer :: ipert, spert, ntot1
            real(dp) :: sadj
            
            ntot1 = lx1*ly1*lz1*nelv
            ipert = npert + 1 - jp
            spert = merge(1, -1, jp == 1)  ! +1 for jp = 1, -1 for jp = 2
            sadj  = merge(-1.0_dp, 1.0_dp, ifadj)   ! antisymmetric in slot space: adjoint = -direct
            
            call col3 (tmp_m1a, t(1,1,1,1,1), bm1, ntot1)            ! W*bm1
            call cmult(tmp_m1a, sadj*spert*betaz, ntot1)             ! +-beta*W*bm1
            
            call col3 (advW_x(1,ipert), tmp_m1a, vxp(1,jp),  ntot1)
            call col2 (advW_x(1,ipert), vtrans(1,1,1,1,1),   ntot1)  ! rho
            call col3 (advW_y(1,ipert), tmp_m1a, vyp(1,jp),  ntot1)
            call col2 (advW_y(1,ipert), vtrans(1,1,1,1,1),   ntot1)  ! rho
            call col3 (advW_z(1,ipert), tmp_m1a, tp(1,1,jp), ntot1)
            call col2 (advW_z(1,ipert), vtrans(1,1,1,1,2),   ntot1)  ! rho*cp
         end subroutine compute_W_advection_planar
      
      
         subroutine makefp_uv_2Dh_planar(advW_x, advW_y, if_advW)
            ! Copy of the core makefp (perturb.f)
            ! Spanwise advection of vxp and vyp injected before makextp if W /= 0.
            implicit none
            include 'SIZE'
            include 'SOLN'
            include 'INPUT'
            include 'TSTEP'
            include 'ADJOINT'
            real(dp), dimension(lv,lpert), intent(in) :: advW_x, advW_y
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
               call add2(bfxp(1,jp), advW_x(1,jp), ntot1)
               call add2(bfyp(1,jp), advW_y(1,jp), ntot1)
            end if
            ! add extrapolation terms
            if (iftran)                                        call makextp
            ! add bd contribution
                                                               call makebdfp
            call neklab_timer_stop(t_planar_makefp)
         end subroutine makefp_uv_2Dh_planar
         
         subroutine makeqp_w_2Dh_planar(advW_z, if_advW)
            ! Copy of the core makeqp (perturb.f)
            ! Spanwise advection of vzp injected before makeabqp if W /= 0.
            ! Note: Expects ifield = 2.
            implicit none
            include 'SIZE'
            include 'SOLN'
            include 'INPUT'
            include 'TSTEP'
            real(dp), dimension(lv,lpert), intent(in) :: advW_z
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
               call add2(bqp(1,ifield-1,jp), advW_z(1,jp), ntot1)
            end if
            ! add extrapolation terms
            if (iftran)                                           call makeabqp
            if ((iftran .and. .not.ifchar) .or.
     &          (iftran .and. .not.ifadvc(ifield) .and. ifchar))  call makebdqp
            call neklab_timer_stop(t_planar_makeqp)
         end subroutine makeqp_w_2Dh_planar
         
         subroutine pressure_matvec_2Dh_planar(ap,wp,h1,h2,h2inv,betaz,intype)
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
            integer :: ntot1, ntot2
            
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
            call bm2_map12t_2Dh(tmp_m1a, wp)                     ! I12^T (bm2 p)
            call col2          (tmp_m1a, wmask, ntot1)           ! mask Dirichlet conditions
            call dssum         (tmp_m1a, lx1, ly1, lz1)          ! make continuous
            call col2          (tmp_m1a, binvm1, ntot1)          ! collate inverse mass matrix
            call col2          (tmp_m1a, h2inv, ntot1)           ! collate (h2)^-1
            call cmult         (tmp_m1a, -betaz**2, ntot1)       ! collate -betaz^2
            call bm2_map12_2Dh (tmp_m2, tmp_m1a)               ! bm2 * I12 (.)
            call sub2          (ap, tmp_m2, ntot2)
            call neklab_timer_stop(t_planar_pres_matvec)
         end subroutine pressure_matvec_2Dh_planar
         
         subroutine helmholtz_matvec_2Dh_planar(Au, u, h1, h2, betaz)
            implicit none
            real, dimension(lx1,ly1,lz1,1), intent(out) :: Au
            real, dimension(lx1,ly1,lz1,1), intent(in) :: u
            real, dimension(lx1,ly1,lz1,1), intent(in) :: h1
            real, dimension(lx1,ly1,lz1,1), intent(in) :: h2
            real betaz
            ! internal
            integer :: imesh, ntot1
            call neklab_timer_start(t_planar_helm_matvec)
            imesh = 1
            ntot1 = lx1*ly1*lz1*nelv
            ! regular Ax
            call axhelm (Au, u, h1, h2, imesh, 1)
            ! added spanwise viscous diffusion term
            call col3 (tmp_m1b, u, vdiff(1,1,1,1,1), ntot1)
            call col2c(tmp_m1b, bm1, betaz**2, ntot1)
            call add2 (Au, tmp_m1b, ntot1)
            call neklab_timer_stop(t_planar_helm_matvec)
         end subroutine helmholtz_matvec_2Dh_planar

         subroutine helm_matvec_planar(Au, u, h1, h2)
            ! Wrapper around helmholtz_matvec_2Dh_planar to match the abstract interface while passing betaz
            real(dp), dimension(lv), intent(out) :: Au
            real(dp), dimension(lv), intent(in)  :: u, h1, h2
            call helmholtz_matvec_2Dh_planar(Au, u, h1, h2, betaz_op)
         end subroutine helm_matvec_planar

         subroutine pres_matvec_planar(Ap, p, h1, h2, h2inv, intype)
            ! Wrapper around the pressure_matvec_2Dh_planar to match the abstract interface while passing betaz
            real(dp), dimension(lp), intent(out) :: Ap
            real(dp), dimension(lp), intent(in)  :: p
            real(dp), dimension(lv), intent(in)  :: h1, h2, h2inv
            integer, intent(in) :: intype
            call pressure_matvec_2Dh_planar(Ap, p, h1, h2, h2inv, betaz_op, intype)
         end subroutine pres_matvec_planar
      
      end module neklab_2Dh_planar
