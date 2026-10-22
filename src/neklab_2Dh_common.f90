      module neklab_2Dh_common
      !! Infrastructure shared by the 2Dh solvers (cartesian: neklab_2Dh_planar,
      !! axisymmetric: neklab_2Dh_axisym):
      !!   - wmask: Dirichlet mask of the third velocity component, which is
      !!     carried in the temperature slot (w, or u_phi);
      !!   - the run checks common to both formulations, called by the model
      !!     inits (init_2Dh_planar, init_2Dh_axisym);
      !!   - the shared stage wrappers and generic solvers of the steppers.
      !! The solvers stop if init has not been called, instead of building
      !! their data lazily on the first time step.
         use stdlib_optval, only: optval
         use LightKrylov, only: dp
         use neklab_nek_setup, only: nek_log_message, nek_log_warning, nek_stop_error
         use neklab_timing, only: neklab_timer_start, neklab_timer_stop,
     &                            t_helm_solve, t_pres_solve, t_pres_prec, t_pres_proj
         implicit none
         include "SIZE"
         include "TOTAL"

         private
         character(len=*), parameter, private :: this_module = 'neklab_2Dh_common'
         integer, parameter :: lv = lx1*ly1*lz1*lelv
         integer, parameter :: lp = lx2*ly2*lz2*lelv

         real(dp), dimension(lx1, ly1, lz1, lelv), public :: wmask
         logical, public :: wmask_defined = .false.
         ! shared stage wrappers of the 2Dh steppers
         public :: setup_step_2Dh, cresvipp_2Dh
         public :: setup_incomprp_2Dh, ebar_2Dh
         public :: incomprp_solve_2Dh, incomprp_update_2Dh

         real(dp), private :: pbasis(lp,mxprev,lpert)
         integer,  private :: nprev(lpert) = 0
         real(dp), private :: dtbdlast = 0.0_dp
         real(dp), dimension(lp), private :: prextr_c, onep, ep
         
         ! abstract interfactes for the helmholtz and pressure matvecs passed to the solvers
         public :: abstract_helmholtz_matvec, abstract_pressure_matvec

         ! BC masks for the w component and checks used by the model inits
         public :: build_wmask
         public :: check_2Dh_properties
         public :: check_2Dh_masks

         ! weighted/unweighted mappings from m1 -> m2 and their discrete transposes
         public ::     map12_2Dh
         public ::     map12T_2Dh
         public :: bm2_map12_2Dh
         public :: bm2_map12T_2Dh
         ! generic solvers for the helmholtz (pCG) and pressure system (pGMRES); concrete matvecs are inputs
         public :: solve_helmholtz_2Dh
         public :: solve_pressure_2Dh
         ! genergic projection / reconstruction routines for the pressure solver; contrete matvecs are inputs
         public ::  setrhs_pressure_2Dh
         public :: gensoln_pressure_2Dh

      !-----------------------------------------------------------------------------------
      !
      !  ABASTRACT INTERFACES for the helmholtz and pressure matvecs passed to the solvers
      !
      !-----------------------------------------------------------------------------------
         
         abstract interface
            subroutine abstract_helmholtz_matvec(Au, u, h1, h2)
               ! based on the core nek helmholtz matvec
               import :: dp, lx1, ly1, lz1, lelv
               real(dp), dimension(lx1*ly1*lz1*lelv), intent(out) :: Au
               real(dp), dimension(lx1*ly1*lz1*lelv), intent(in)  :: u, h1, h2
            end subroutine abstract_helmholtz_matvec
            subroutine abstract_pressure_matvec(Ap, p, h1, h2, h2inv, intype)
               ! based on the core nek pressure matvec
               import :: dp, lx1, ly1, lz1, lx2, ly2, lz2, lelv
               real(dp), dimension(lx2*ly2*lz2*lelv), intent(out) :: Ap
               real(dp), dimension(lx2*ly2*lz2*lelv), intent(in)  :: p
               real(dp), dimension(lx1*ly1*lz1*lelv), intent(in)  :: h1, h2, h2inv
               integer, intent(in) :: intype
            end subroutine abstract_pressure_matvec
         end interface

      contains

         subroutine build_wmask(verbose)
      !! Dirichlet mask of the third velocity component: zero wherever the full
      !! velocity vector is prescribed.
            logical, optional, intent(in) :: verbose
            ! internal
            character(len=*), parameter :: this_procedure = 'build_wmask'
            character(len=128) :: msg
            character(len=3) :: cb
            integer :: ie, ifc, nfaces, ntot, ifield_bak, nzero
            real(dp), external :: vlsum
            logical :: verb

            nfaces = 2*ndim
            ntot = lx1*ly1*lz1*nelv
            verb = optval(verbose, .false.)

            call rone(wmask, ntot)

            do ie = 1, nelv
               do ifc = 1, nfaces
                  cb = cbc(ifc, ie, 1)   ! velocity boundary conditions
                  select case (cb)
      ! Full velocity vector prescribed -> w is prescribed too.
                  case ('W  ', 'v  ', 'V  ', 'vl ', 'VL ', 'mv ', 'MV ', 'ON ', 'on ')
                     call facev(wmask, ie, ifc, 0.0_dp, lx1, ly1, lz1)
      ! Everything else leaves w free:
      !   'SYM'          - normal lies in the x-y plane, w is tangential
      !   'O  ','o  '    - outflow, natural BC
      !   'SL ','sl '    - slip, tangential components free
      !   'shl'          - stress-free tangential
      !   'P  ','p  '    - periodic, handled by dssum
      !   'E  '          - interior
                  end select
               end do
            end do

      ! Propagate zeros across element interfaces. Called directly rather than
      ! through opdsop so that the 2D dimension guard does not apply. 'MUL' is
      ! the correct operator for a 0/1 mask: any node that is Dirichlet in one
      ! element becomes Dirichlet in all elements sharing it.
            ifield_bak = ifield
            ifield = 1
            call dsop(wmask, 'MUL', lx1, ly1, lz1)
            ifield = ifield_bak

            wmask_defined = .true.

            if (verb) then
               nzero = ntot - nint(vlsum(wmask, ntot))
               write (msg, '(A,I10,A,I10,A)')
     &            ' neklab_masks: wmask built, ', nzero, ' of ', ntot,
     &            ' local nodes constrained.'
               call nek_log_message(msg, this_module, this_procedure)
            end if

         end subroutine build_wmask

         subroutine check_2Dh_properties(caller)
      !! The third velocity component lives in the temperature slot, so the
      !! scalar must exist and diffuse/convect exactly like momentum.
            character(len=*), intent(in) :: caller
            real(dp), parameter :: diff_tol = 1.0e-08_dp
            character(len=128) :: msg
            if (.not. ifheat)  then
               msg = '2Dh requires a temperature field (ifheat): it carries the vzp component.'
               call nek_stop_error(msg, this_module, caller)
            end if
            if (abs(cpfld(1,1) - cpfld(2,1)) > diff_tol) then
               write(msg,'(A,2E18.7)') 'viscosity, conductivity = ', cpfld(1,1), cpfld(2,1)
               call nek_log_message(msg, this_module, caller)
               call nek_stop_error('2Dh requires conductivity = viscosity.', this_module, caller)
            end if
            if (abs(cpfld(1,2) - cpfld(2,2)) > diff_tol) then
               write(msg,'(A,2E18.7)') 'rho, rhoCp = ', cpfld(1,2), cpfld(2,2)
               call nek_log_message(msg, this_module, caller)
               call nek_stop_error('2Dh requires rhoCp = rho.', this_module, caller)
            end if
         end subroutine check_2Dh_properties

         subroutine check_2Dh_masks(caller)
      !! The solver masks the third component with wmask (from the velocity
      !! BCs); the deck sets the temperature BCs separately. Warn if the two
      !! disagree. Counts are node instances summed over all ranks.
            character(len=*), intent(in) :: caller
            integer :: ix, iy, iz, ie, nbad(2), nsum(2)
            character(len=256) :: msg

            nbad = 0
            do ie = 1, nelv
               do iz = 1, lz1
               do iy = 1, ly1
               do ix = 1, lx1
                  if (wmask(ix,iy,iz,ie) == 0.0_dp .and. tmask(ix,iy,iz,ie,1) /= 0.0_dp)
     &               nbad(1) = nbad(1) + 1
                  if (wmask(ix,iy,iz,ie) /= 0.0_dp .and. tmask(ix,iy,iz,ie,1) == 0.0_dp)
     &               nbad(2) = nbad(2) + 1
               end do
               end do
               end do
            end do
            call igop(nbad, nsum, '+  ', 2)
            if (nbad(1) > 0) then
               write(msg,'(A,I10,A)') 'u3 Dirichlet (wmask) but temperature free at ', nbad(1), ' nodes.'
               call nek_log_warning(msg, this_module, caller)
            end if
            if (nbad(2) > 0) then
               write(msg,'(A,I10,A)') 'temperature Dirichlet but u3 free (wmask) at ', nbad(2), ' nodes.'
               call nek_log_warning(msg, this_module, caller)
            end if
            if (nbad(1) > 0 .or. nbad(2) > 0) then
               call nek_log_message('Ignore this warning when running the linear solver for antisymmetric modes.', this_module, caller)
            end if
         end subroutine check_2Dh_masks

      !-------------------------------------------------------------------------------
      !
      !     COMMON ROUTINES for planar (cartesian) and axisymmetric 2Dh solvers
      !
      !-------------------------------------------------------------------------------

      !-------------------------------------------------------------------------------
      !  Shared stage wrappers of the 2Dh steppers (names follow perturb.f)
      !-------------------------------------------------------------------------------

         subroutine setup_step_2Dh()
      !! Per-step solver setup, as in Nek's nek_advance.
            if (istep == 1) nprev(:) = 0     ! new linear pass: discard the pressure projection basis
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
         end subroutine setup_step_2Dh

         subroutine cresvipp_2Dh(resv1, resv2, resv3, gradw_jp)
      !! Momentum right-hand side of slot jp without the H u term (cf. cresvipp):
      !! Dirichlet values, extrapolated pressure gradient (in-plane: opgradt,
      !! third component: gradw_jp from makecpl), forcing.
            real(dp), dimension(lv), intent(out) :: resv1, resv2, resv3
            real(dp), dimension(lv), intent(in)  :: gradw_jp
            ! internal
            integer :: ntot1
            ntot1 = lx1*ly1*lz1*nelv
            call bcdirvc  (vxp(1,jp), vyp(1,jp), vzp(1,jp), v1mask, v2mask, v3mask)
            call col2     (tp(1,1,jp), wmask, ntot1)
            call extrapprp(prextr_c)
            call opgradt  (resv1, resv2, resv3, prextr_c)
            call opadd2   (resv1, resv2, resv3, bfxp(1,jp), bfyp(1,jp), bfzp(1,jp))
            call copy     (resv3, gradw_jp, ntot1)
            call add2     (resv3, bqp(1,1,jp), ntot1)
         end subroutine cresvipp_2Dh

         subroutine setup_incomprp_2Dh(h1, h2, h2inv)
      !! Helmholtz coefficients of the pressure step (intype = 1). The projection
      !! basis is discarded when dt/bd(1) changes.
            real(dp), dimension(lv), intent(out) :: h1, h2, h2inv
            ! internal
            real(dp) :: dtbd
            integer :: ntot1
            ntot1  = lx1*ly1*lz1*nelv
            ifield = 1
            dtbd   = bd(1)/dt
            if (dtbd /= dtbdlast) then
               nprev(:) = 0                  ! operator changed: basis is stale
               dtbdlast = dtbd
            end if
            call rzero  (h1, ntot1)
            call cmult2 (h2, vtrans(1,1,1,1,ifield), dtbd, ntot1)
            call invers2(h2inv, h2, ntot1)
         end subroutine setup_incomprp_2Dh

         subroutine ebar_2Dh(ebar, matvec, h1, h2, h2inv, wavenumber)
      !! Mean of E*1. For a nonzero wavenumber E is not singular, but the hsmg
      !! preconditioner removes the mean; the solver adds it back with ebar.
            real(dp), intent(out) :: ebar
            procedure(abstract_pressure_matvec) :: matvec
            real(dp), dimension(lv), intent(in) :: h1, h2, h2inv
            real(dp), intent(in) :: wavenumber
            ! internal
            logical :: if_singular
            real(dp), external :: glsum
            integer :: ntot2
            ntot2 = lx2*ly2*lz2*nelv
            ebar = 0.0_dp
            if_singular = (wavenumber == 0.0_dp)
            if (.not. if_singular .and. ifvcor) then
               call rone(onep, ntot2)
               call matvec(ep, onep, h1, h2, h2inv, 1)
               ebar = glsum(ep, ntot2)
            end if
         end subroutine ebar_2Dh

         subroutine incomprp_solve_2Dh(dpr_jp, h1, h2, h2inv, matvec, ebar, wavenumber, maxit, label)
      !! Pressure correction of slot jp: E dp = rhs, with projection (param(95)).
            real(dp), dimension(lp), intent(inout) :: dpr_jp
            real(dp), dimension(lv), intent(in)    :: h1, h2, h2inv
            procedure(abstract_pressure_matvec) :: matvec
            real(dp), intent(in) :: ebar, wavenumber
            integer, intent(in)  :: maxit
            character(len=*), intent(in) :: label
            ! internal
            logical :: if_singular
            character(len=2) :: info_str
            logical :: ifprjp
            integer :: istart, iter
            info_str = merge('Re', 'Im', jp == 1)
            if_singular = (wavenumber == 0.0_dp)
            if (if_singular) call ortho(dpr_jp)
            istart = param(95)
            ifprjp = (istep >= istart .and. istart /= 0)
            if (ifprjp) then
               call  setrhs_pressure_2Dh(dpr_jp, h1, h2, h2inv, pbasis(1,1,jp), nprev(jp),
     &                                   matvec, info_str)
            end if
            call      solve_pressure_2Dh(dpr_jp, h1, h2, h2inv, matvec, 1, iter, ebar,
     &                                   if_singular, maxit, label, wavenumber, info_str)
            if (ifprjp) then
               call gensoln_pressure_2Dh(dpr_jp, h1, h2, h2inv, pbasis(1,1,jp), nprev(jp),
     &                                   matvec)
            end if
         end subroutine incomprp_solve_2Dh

         subroutine incomprp_update_2Dh(dpr_jp, dw, h2inv)
      !! Velocity and pressure update of slot jp from the pressure correction
      !! (cf. the end of incomprp). dw is the third-component correction, computed
      !! by the model (incomprp_dw_*).
            real(dp), dimension(lp), intent(in) :: dpr_jp
            real(dp), dimension(lv), intent(in) :: dw, h2inv
            real(dp), dimension(lx1,ly1,lz1,lelv) :: w1, w2, w3, resv1, resv2, resv3, dv1, dv2, dv3
            common /SCRNS/ w1, w2, w3, resv1, resv2, resv3, dv1, dv2, dv3
            integer :: ntot1, ntot2
            ntot1 = lx1*ly1*lz1*nelv
            ntot2 = lx2*ly2*lz2*nelv
            call opgradt  (w1, w2, w3, dpr_jp)
            call opbinv   (dv1, dv2, dv3, w1, w2, w3, h2inv)
            call opadd2   (vxp(1,jp), vyp(1,jp), vzp(1,jp), dv1, dv2, dv3)
            call add2     (tp(1,1,jp), dw, ntot1)
            call extrapprp(prextr_c)
            call lagpresp
            call add3     (prp(1,jp), prextr_c, dpr_jp, ntot2)
         end subroutine incomprp_update_2Dh
      
      !-------------------------------------------------------------------------------
      !  Weighted/unweighted mappings from m1 -> m2 and their discrete transposes

      !--------------------------------------------------------------------
      !  Nek core:
      !--------------------------------------------------------------------
      !  D : M1 -> I12 -> M2
      !
      ! navier1.f, multd:
      !
      ! CARTESIAN:
      ! 624:      call mxm (dxm12,lx2,x(1,e),lx1,ta1,nyz1)      ! r: derivative GLL -> GL
      ! 625:      call mxm (ta1,lx2,iytm12,ly1,dx(1,e),ly2)     ! s: interpolation GLL -> GL
      ! 626:      call col2 (dx(1,e),rm2(1,e),nxyz2)            ! * dr/dx (incl. Jacobian)
      ! 627:      call mxm (ixm12,lx2,x(1,e),lx1,ta1,nyz1)      ! r: interpolation GLL -> GL
      ! 628:      call mxm (ta1,lx2,dytm12,ly1,ta3,ly2)         ! s: derivative GLL -> GL
      ! 629:      call addcol3 (dx(1,e),ta3,sm2(1,e),nxyz2)     ! * ds/dx
      !           ...
      ! 693:      call col2 (dx(1,e),w3m2,nxyz2)                ! GL weights
      !
      ! CYLINDRICAL (ifaxis, added for isd=2; elements on the axis are treated with ifrzer):
      ! 725:      call mxm     (ixm12,lx2,x(1,e),lx1,ta1,ly1)   ! I12, r direction
      ! 726:      call mxm     (ta1,lx2,iytm12,ly1,ta2,ly2)     ! I12, s direction (iytm12 <- ictm12)
      ! 727:      call col3    (ta3,bm2(1,1,1,e),ta2,nxyz2)     ! * bm2
      ! 728:      call invcol2 (ta3,ym2(1,1,1,e),nxyz2)         ! / R  (on M2)
      ! 729:      call add2    (dx(1,e),ta3,nxyz2)              ! 
      !
      !-------------------------------------------------------------------
      !  D^T : M2 -> I12^T -> M1
      !
      ! navier1.f, cdtp:
      !
      ! CARTESIAN:
      ! 402:      call col3 (wx,w3m2,x(1,e),nxyz2)              ! weights first
      !           ...
      ! 430:      call col3 (ta1,wx,rm2(1,e),nxyz2)
      ! 431:      call mxm  (dxtm12,lx1,ta1,lx2,ta2,nyz2)       ! dxm12^T
      ! 432:      call mxm  (ta2,lx1,iym12,ly2,dtx(1,e),ly1)    ! (iytm12)^T = iym12
      !           
      ! 434:      call col3 (ta1,wx,sm2(1,e),nxyz2)
      ! 435:      call mxm  (ixtm12,lx1,ta1,lx2,ta2,nyz2)       ! ixm12^T
      ! 436:      call mxm  (ta2,lx1,dym12,ly2,ta1,ly1)         ! (dytm12)^T
      ! 
      ! CYLINDRICAL (ifaxis, added for isd=2; elements on the axis are treated with ifrzer):
      ! 521:      call col3    (ta1,x(1,e),bm2(1,1,1,e),nxyz2)  ! * bm2
      ! 522:      call invcol2 (ta1,ym2(1,1,1,e),nxyz2)         ! / R  (on M2)
      ! 523:      call mxm     (ixtm12,lx1,ta1,lx2,ta2,ly2)     ! I12^T, r direction
      ! 524:      call mxm     (ta2,lx1,iym12,ly2,ta1,ly1)      ! I12^T, s direction
      ! 525:      call add2    (dtx(1,e),ta1,nxyz1)
      !      
      ! -------------------------------------------------------------------

         subroutine map12_2Dh(y_m2, x_m1, ie)
            ! y_m2 = I12 @ x_m1
            ! element local operation
            real(dp), dimension(lx2,ly2), intent(out) :: y_m2
            real(dp), dimension(lx1,ly1), intent(in)  :: x_m1
            integer, intent(in) :: ie
            real(dp), dimension(lx2,ly1) :: ta
            call mxm(ixm12, lx2, x_m1, lx1, ta, ly1)
            if (ifaxis) then
               if (ifrzer(ie)) then ! element on the axis
                  call mxm(ta, lx2, iatm12, ly1, y_m2, ly2)
               else
                  call mxm(ta, lx2, ictm12, ly1, y_m2, ly2)
               end if
            else
               call mxm(ta, lx2, iytm12, ly1, y_m2, ly2)
            end if
         end subroutine map12_2Dh

         subroutine map12T_2Dh(y_m1, x_m2, ie)
            ! y_m1 = I12^T @ x_m2
            ! element local operation
            real(dp), dimension(lx1,ly1), intent(out) :: y_m1
            real(dp), dimension(lx2,ly2), intent(in)  :: x_m2
            integer, intent(in) :: ie
            real(dp), dimension(lx1,ly2) :: ta
            call mxm(ixtm12, lx1, x_m2, lx2, ta, ly2)
            if (ifaxis) then
               if (ifrzer(ie)) then ! element on the axis
                  call mxm(ta, lx1, iam12, ly2, y_m1, ly1)
               else
                  call mxm(ta, lx1, icm12, ly2, y_m1, ly1)
               end if
            else
               call mxm(ta, lx1, iym12, ly2, y_m1, ly1)
            end if
         end subroutine map12T_2Dh
         
         subroutine bm2_map12_2Dh(y_m2, x_m1, const_fld)
            ! y_m2 = [const_fld *] bm2 * I12 @ x_m1
            ! apply to all elements
            real(dp), dimension(lx2*ly2*lz2,lelv), intent(out) :: y_m2
            real(dp), dimension(lx1*ly1*lz1,lelv), intent(in)  :: x_m1
            real(dp), dimension(lx2*ly2*lz2,lelv), intent(in), optional :: const_fld
            integer :: ie, n2
            n2 = lx2*ly2*lz2
            do ie = 1, nelv
               call map12_2Dh(y_m2(1,ie), x_m1(1,ie), ie)
               call col2 (y_m2(1,ie), bm2(1,1,1,ie), n2)
               if (present(const_fld)) call col2(y_m2(1,ie), const_fld(1,ie), n2)
            end do
         end subroutine bm2_map12_2Dh
         
         subroutine bm2_map12T_2Dh(y_m1, x_m2, const_fld)
            ! y_m1 = I12^T @ ( [const_fld *] bm2 * x_m2 )
            ! apply to all elements
            real(dp), dimension(lx1*ly1*lz1,lelv), intent(out) :: y_m1
            real(dp), dimension(lx2*ly2*lz2,lelv), intent(in)  :: x_m2
            real(dp), dimension(lx2*ly2*lz2,lelv), intent(in), optional :: const_fld
            real(dp), dimension(lx2*ly2*lz2) :: wk
            integer :: ie, n2
            n2 = lx2*ly2*lz2
            do ie = 1, nelv
               call col3(wk, x_m2(1,ie), bm2(1,1,1,ie), n2)
               if (present(const_fld)) call col2(wk, const_fld(1,ie), n2)
               call map12T_2Dh(y_m1(1,ie), wk, ie)
            end do
         end subroutine bm2_map12T_2Dh

      !-------------------------------------------------------------------------------
      !  Generic solvers for the helmholtz (pCG) and pressure system (pGMRES)

      ! -------------------------------------------------------------------
      !  Essentially a copy of the helmholtz solver in navier1.f, but with the abstract matvec as input.
      ! -------------------------------------------------------------------
         subroutine solve_helmholtz_2Dh(x, f, h1, h2, mask, mult, imsh, tin, maxit, isd,
     &                                  binv, name, matvec)
            include 'FDMH1'
            real(dp), dimension(lv), intent(out) :: x
            real(dp), dimension(lv), intent(in)  :: f, h1, h2, mask, mult, binv
            integer, intent(in) :: imsh, maxit, isd
            real(dp), intent(in) :: tin
            character(len=7), intent(in) :: name
            procedure(abstract_helmholtz_matvec) :: matvec
            ! internal
            character(len=*), parameter :: this_procedure = 'solve_helmholtz_2Dh'
            integer :: iter, krylov, n, nel, niter, nxyz
            real(dp) :: alpha, beta, alphm, fmax
            real(dp) :: rbn0, rbn2, rho, rho0, rtz1, rtz2, tol, vol
            logical          ifprint, ifhzpc
            common  /cprint/ ifprint, ifhzpc
            logical ifdfrm, iffast, ifh2, ifsolv
            common /fastmd/ ifdfrm(lelt), iffast(lelt), ifh2, ifsolv
            logical :: ifmcor, ifprint_helm
            integer, parameter :: lg = lx1*ly1*lz1*lelt
            real(dp), dimension(lg) :: d, r, w, p, z
            real(dp) :: scalar(2)
            common /scrcg/ d, scalar
            common /scrmg/ r, w, p, z
            integer, parameter :: maxcg = 900
            real(dp), dimension(maxcg) :: diagt, upper
            common /tdarray/ diagt, upper
            integer :: niterhm
            common /iterhm/ niterhm
            real(dp), external :: glamax, glsc3, vlsc3, vlsc32

            call neklab_timer_start(t_helm_solve)
            ! Initialization
            call rzero(diagt, maxcg)
            call rzero(upper, maxcg)
            rho  = 0.0_dp
            nxyz = lx1*ly1*lz1
            nel  = nelv
            vol  = volvm1
            if (imsh == 2) nel = nelt
            if (imsh == 2) vol = voltm1
            n = nel*nxyz

            ! Set tolerances
            tol = abs(tin)
            if (restol(ifield) /= 0) tol = abs(restol(ifield))
            if (tin < 0) tol = abs(tin)
            niter = min(maxit, maxcg)

            if (.not. ifsolv) then
               call setfast(h1, h2, imsh)
               ifsolv = .true.
            end if

            ! Setup the Jacobi preconditioner: d = 1/diag(H_nek), assembled
            call setprec(d, h1, h2, imsh, isd)

            ! Check for trivial solution
            call copy (r, f, n)
            call rzero(x, n)
            call rzero(p, n)
            fmax = glamax(f, n)
            if (fmax == 0.0_dp) then
               call neklab_timer_stop(t_helm_solve)
               return
            end if

            krylov  = 0
            rtz1    = 1.0_dp
            niterhm = 0
            do iter = 1, niter

               ! apply the Jacobi preconditioner
               call col3(z, r, d, n)                              ! z = D^-1 r

               rtz2 = rtz1
               scalar(1) = vlsc3(z, r, mult, n)
               if (param(18) == 1) then
                  scalar(2) = vlsc3(r, r, mult, n)
               else
                  scalar(2) = vlsc32(r, mult, binv, n)
               end if
               call gop(scalar, w, '+  ', 2)
               rtz1 = scalar(1)
               rbn2 = sqrt(scalar(2)/vol)                         ! residual norm
               if (iter == 1) rbn0 = rbn2                         ! record initial residual norm
               if (param(22) < 0) tol = abs(param(22))*rbn0
               if (tin < 0)       tol = abs(tin)*rbn0

               ifprint_helm = .false.
               if (nio == 0 .and. ifprint .and. param(74) /= 0) ifprint_helm = .true.
               if (nio == 0 .and. istep == 1)                   ifprint_helm = .true.
               if (ifprint_helm) write(6,3002) istep, '  Hmholtz '//name, iter, rbn2, h1(1), tol, h2(1), ifmcor

               if (rbn2 <= tol .and. (iter > 1 .or. istep <= 5)) then
                  niter = iter - 1
                  if (nio == 0) write(6,3000) istep, '  Hmholtz '//name, niter, rbn2, rbn0, tol
                  goto 9999
               end if

               ! Lanczos iteration
               beta = rtz1/rtz2
               if (iter == 1) beta = 0.0_dp

               ! Update search direction
               call add2s1(p, z, beta, n)                         ! p = z + beta*p

               call matvec(w, p, h1, h2)                          ! matvec
               call dssum (w, lx1, ly1, lz1)
               call col2  (w, mask, n)

               ! Compute step size
               rho0  = rho
               rho   = glsc3(w, p, mult, n)
               alpha = rtz1/rho
               alphm = -alpha

               ! Update solution and residual
               call add2s2(x, p, alpha, n)                        ! x = x + alpha*p
               call add2s2(r, w, alphm, n)                        ! r = r - alpha*w

               ! Update the tridiagonal coefficient matrix
               if (iter == 1) then
                  krylov = krylov + 1
                  diagt(iter) = rho/rtz1
               else if (iter <= maxcg) then
                  krylov = krylov + 1
                  diagt(iter)   = (beta**2*rho0 + rho)/rtz1
                  upper(iter-1) = -beta*rho0/sqrt(rtz2*rtz1)
               end if
            end do
            niter = iter - 1
            if (nio == 0) write(6,3001) istep, '  Error Hmholtz '//name, niter, rbn2, rbn0, tol
            call nek_log_warning('Helmholtz '//name//' not converged.', this_module, this_procedure)

 3000       format(i11,a,1x,I7,1p4E13.4)
 3001       format(i11,a,1x,I7,1p4E13.4)
 3002       format(i11,a,1x,I7,1p4E13.4,l4)
 9999       continue
            niterhm = niter
            ifsolv  = .false.

            call neklab_timer_stop(t_helm_solve)

         end subroutine solve_helmholtz_2Dh
        
      ! -------------------------------------------------------------------
      !  Essentially a copy of the uzawa_gmres in gmres.f, but with the abstract matvec as input.
      !      if_ortho : project out the nullspace if the operator is singular (wavenumber = 0 only)
      ! -------------------------------------------------------------------
         subroutine solve_pressure_2Dh(res, h1, h2, h2inv, matvec, intype, iter, ebar,
     &                                 if_ortho, maxit, label, wavenumber, info_str)
            include 'GMRES'
            real(dp), dimension(lp), intent(inout) :: res
            real(dp), dimension(lv), intent(in)    :: h1, h2, h2inv
            procedure(abstract_pressure_matvec) :: matvec
            integer, intent(in)  :: intype, maxit
            integer, intent(out) :: iter
            real(dp), intent(in) :: ebar, wavenumber
            logical, intent(in)  :: if_ortho
            character(len=*), intent(in) :: label, info_str
            ! internal
            character(len=*), parameter :: this_procedure = 'solve_pressure_2Dh'
            real(dp) :: divex
            common  /ctolpr/ divex
            logical          ifprint
            common  /cprint/ ifprint
            real(dp), dimension(lgmres) :: wk1
            real(dp) :: alpha, l, temp, div0, ratio, rnorm, tolpss
            real(dp) :: etime_p
            real(dp) :: etime1, etime2, dnekclock
            integer  :: i, j, k, m, iconv, ntot2
            logical, save  :: iflag = .false.
            real(dp), save :: norm_fac
            real(dp), external :: vlsc2, glsc2, glsum
            character(len=128) :: msg

            call neklab_timer_start(t_pres_solve)
            ntot2 = lx2*ly2*lz2*nelv

            ! compute the lower and upper triangular factors of the preconditioner
            if (.not. iflag) then
               iflag = .true.
               call uzawa_gmres_split0(ml_gmres, mu_gmres, bm2, bm2inv, ntot2)
               norm_fac = 1.0_dp/sqrt(volvm2)
            end if

            etime1  = dnekclock()
            etime_p = 0.0_dp
            divex   = 0.0_dp
            iter    = 0
            m       = lgmres

            call chktcg2(tolps, res, iconv)
            if (param(21) > 0 .and. tolps > abs(param(21))) tolps = abs(param(21))
            if (istep == 0) tolps = 1.e-4
            tolpss = tolps

            iconv = 0
            call rzero(x_gmres, ntot2)

            do while (iconv == 0 .and. iter < maxit)

               if (iter == 0) then
                  ! initial residual
                  call col3(r_gmres, ml_gmres, res, ntot2)                 ! r = L^-1 res
               else
                  ! update residual
                  call copy(r_gmres, res, ntot2)
                  call matvec(w_gmres, x_gmres, h1, h2, h2inv, intype)     ! matvec
                  call add2s2(r_gmres, w_gmres, -1.0_dp, ntot2)
                  call col2(r_gmres, ml_gmres, ntot2)                      ! r = L^-1 (res - A x)
               end if

               gamma_gmres(1) = sqrt(glsc2(r_gmres, r_gmres, ntot2))       ! prec. residual norm

               ! record initial residual norm for divergence check
               if (iter == 0) then
                  div0 = gamma_gmres(1)*norm_fac
                  if (param(21) < 0) tolpss = abs(param(21))*div0
               end if

               ! check for lucky convergence
               rnorm = 0.0_dp
               if (gamma_gmres(1) == 0.0_dp) goto 9000
               temp = 1.0_dp/gamma_gmres(1)
               call cmult2(v_gmres(1,1), r_gmres, temp, ntot2)             ! v(:,1) = r / gamma(1)

               do j = 1, m
                  iter = iter + 1

                  ! apply inverse of the lower triangular factor U
                  call col3(w_gmres, mu_gmres, v_gmres(1,j), ntot2)        ! w = U^-1 v

                  ! apply hsmg preconditioner (internally assumes a singular operator)
                  etime2 = dnekclock()
                  call neklab_timer_start(t_pres_prec)
                  call hsmg_solve(z_gmres(1,j), w_gmres)                   ! z = M^-1 w
                  call neklab_timer_stop(t_pres_prec)
                  etime_p = etime_p + dnekclock() - etime2

                  ! add back the mean if the operator is not singular (betaz /= 0)
                  if (ebar /= 0.0_dp) call cadd(z_gmres(1,j), glsum(w_gmres,ntot2)/ebar, ntot2)

                  call matvec(w_gmres, z_gmres(1,j), h1, h2, h2inv, intype)! matvec

                  ! apply inverse of the upper triangular factor L
                  call col2(w_gmres, ml_gmres, ntot2)

                  ! Gram-Schmidt orthogonalization
                  do i = 1, j
                     h_gmres(i,j) = vlsc2(w_gmres, v_gmres(1,i), ntot2)
                  end do
                  call gop(h_gmres(1,j), wk1, '+  ', j)
                  do i = 1, j
                     call add2s2(w_gmres, v_gmres(1,i), -h_gmres(i,j), ntot2)
                  end do

                  ! apply Givens rotations
                  do i = 1, j-1
                     temp = h_gmres(i,j)
                     h_gmres(i  ,j) =  c_gmres(i)*temp + s_gmres(i)*h_gmres(i+1,j)
                     h_gmres(i+1,j) = -s_gmres(i)*temp + c_gmres(i)*h_gmres(i+1,j)
                  end do

                  alpha = sqrt(glsc2(w_gmres, w_gmres, ntot2))             ! new residual norm
                  rnorm = 0.0_dp
                  if (alpha == 0.0_dp) goto 900                            ! converged?

                  ! compute new Givens rotation
                  l = sqrt(h_gmres(j,j)*h_gmres(j,j) + alpha**2)
                  temp = 1.0_dp/l
                  c_gmres(j) = h_gmres(j,j)*temp
                  s_gmres(j) = alpha*temp
                  h_gmres(j,j) = l
                  gamma_gmres(j+1) = -s_gmres(j)*gamma_gmres(j)
                  gamma_gmres(j)   =  c_gmres(j)*gamma_gmres(j)

                  rnorm = abs(gamma_gmres(j+1))*norm_fac
                  ratio = rnorm/div0
                  if (ifprint .and. nio == 0) write(6,66) iter, tolpss, rnorm, div0, ratio, istep
   66             format(i5,1p4e12.5,i8,' Divergence')

                  if (rnorm < tolpss) goto 900                             ! converged?
                  if (j == m) goto 1000                                    ! not converged, restart

                  temp = 1.0_dp/alpha
                  call cmult2(v_gmres(1,j+1), w_gmres, temp, ntot2)        ! normalize
               end do
  900          iconv = 1
 1000          continue

               ! solution update through back substitution
               ! compute coefficients ...
               do k = j, 1, -1
                  temp = gamma_gmres(k)
                  do i = j, k+1, -1
                     temp = temp - h_gmres(k,i)*c_gmres(i)
                  end do
                  c_gmres(k) = temp/h_gmres(k,k)
               end do
               ! ... for the linear combination of the basis vectors
               do i = 1, j
                  call add2s2(x_gmres, z_gmres(1,i), c_gmres(i), ntot2)
               end do
            end do
 9000       continue
            divex = rnorm

            if (iconv == 0 .and. rnorm >= tolpss) then
               write(msg,'(A,I6,A,1pE12.4)') 'pressure GMRES not converged after ', iter,
     &              ' iterations, residual ', rnorm
               call nek_log_warning(msg, this_module, this_procedure)
            end if

            call copy(res, x_gmres, ntot2)
            if (if_ortho) call ortho(res)

            etime1 = dnekclock() - etime1
            if (nio == 0) write(6,9998) istep, '  U-PRES gmres  '//label//' = ', wavenumber,
     &         info_str, iter, divex, div0, tolpss, etime_p, etime1
 9998       format(i11,a,F5.1,1X,A,1X,I6,1p5e13.4)

            call neklab_timer_stop(t_pres_solve)

         end subroutine solve_pressure_2Dh

      !-------------------------------------------------------------------------------
      !  Generic pressure projection routines; matvecs are inputs

      !-------------------------------------------------------------------------------
      !  Essentially a copy of setrhs_/gensoln_/econj
      !-------------------------------------------------------------------------------
         subroutine setrhs_pressure_2Dh(dpr, h1, h2, h2inv, proj_set, niprev, matvec, info_str)
            real(dp), dimension(lp), intent(inout) :: dpr
            real(dp), dimension(lv), intent(in)    :: h1, h2, h2inv
            real(dp), dimension(lp,mxprev), intent(inout) :: proj_set
            integer, intent(inout) :: niprev
            procedure(abstract_pressure_matvec) :: matvec
            character(len=2), intent(in) :: info_str
            ! internal
            integer, parameter :: ltot2 = lx2*ly2*lz2*lelv
            real(dp) :: pbar(ltot2), pnew(ltot2)
            real(dp) :: alpha(mxprev), work(mxprev)
            common /orthox/ pbar, pnew
            common /orthos/ alpha, work
            integer :: ntot2, i
            real(dp), external :: glsc3, vlsc2
            real(dp) :: alpha1, alpha2, ratio
 
            call neklab_timer_start(t_pres_proj)

            ntot2 = lx2*ly2*lz2*nelv
            call rzero(pbar, ntot2)

            ! Initial residual.
            if (niprev == 0) then
               call neklab_timer_stop(t_pres_proj)
               return
            end if
            alpha1 = glsc3(dpr, dpr, bm2inv, ntot2)
            if (alpha1 > 0) then
               alpha1 = sqrt(alpha1/volvm2)
            else
               call neklab_timer_stop(t_pres_proj)
               return
            end if

            ! Build pbar = proj_set @ (proj_set^T @ dpr)
            do i = 1, niprev
               alpha(i) = vlsc2(dpr, proj_set(1,i), ntot2)
            end do
            call gop(alpha, work, '+  ', niprev)
            do i = 1, niprev
               call add2s2(pbar, proj_set(1,i), alpha(i), ntot2)
            end do

            ! Compute the new residual dpr = dpr - E @ pbar
            call matvec(pnew, pbar, h1, h2, h2inv, 1)
            call sub2(dpr, pnew, ntot2)

            ! Check the reduction in the residual.
            alpha2 = glsc3(dpr, dpr, bm2inv, ntot2)
            if (alpha2 > 0) then
               alpha2 = sqrt(alpha2/volvm2)
               ratio  = alpha1/alpha2
               if (nio == 0) write(6,13) istep, '  Project PRES '//info_str,
     &                                   alpha2, alpha1, ratio, niprev, mxprev
            end if
 13         format(i11,a,6x,1p3e13.4,i4,i4)

            call neklab_timer_stop(t_pres_proj)

         end subroutine setrhs_pressure_2Dh
 
         subroutine gensoln_pressure_2Dh(dpr, h1, h2, h2inv, proj_set, niprev, matvec)
            real(dp), dimension(lp), intent(inout) :: dpr
            real(dp), dimension(lv), intent(in)    :: h1, h2, h2inv
            real(dp), dimension(lp,mxprev), intent(inout) :: proj_set
            integer, intent(inout) :: niprev
            procedure(abstract_pressure_matvec) :: matvec
            ! internal
            integer, parameter :: ltot2 = lx2*ly2*lz2*lelv
            real(dp) :: pbar(ltot2), pnew(ltot2)
            real(dp) :: alpha(mxprev), work(mxprev)
            common /orthox/ pbar, pnew
            common /orthos/ alpha, work
            integer :: ntot2, mprev, ierr

            call neklab_timer_start(t_pres_proj)
 
            mprev = min(int(param(93)), mxprev)
            ntot2 = lx2*ly2*lz2*nelv                                                       

            ! add new vector to proj_set and E-conjugate it wrt. the previous vectors
            if (niprev < mprev) then
               niprev = niprev + 1
               call copy(proj_set(1,niprev), dpr, ntot2)                                  ! add to proj_set
               call add2(dpr, pbar, ntot2)                                                ! reconstruct solution
               call econj_pressure_2Dh(proj_set, niprev, h1, h2, h2inv, matvec, ierr)     ! E-conjugate new vector
               if (ierr == 1) then
                  niprev = 1
                  call copy(proj_set(1,niprev), dpr, ntot2)                               ! restart proj_set
                  call econj_pressure_2Dh(proj_set, niprev, h1, h2, h2inv, matvec, ierr)  ! E-conjugate new vector
               end if
            else
               niprev = 1
               call add2(dpr, pbar, ntot2)                                                ! reconstruct solution
               call copy(proj_set(1,niprev), dpr, ntot2)                                  ! restart proj_set
               call econj_pressure_2Dh(proj_set, niprev, h1, h2, h2inv, matvec, ierr)     ! E-conjugate new vector
            end if

            call neklab_timer_stop(t_pres_proj)

         end subroutine gensoln_pressure_2Dh
 
         subroutine econj_pressure_2Dh(proj_set, nprev, h1, h2, h2inv, matvec, ierr)
            ! E-conjugate the last vector in proj_set wrt. the previous vectors
            real(dp), dimension(lp,mxprev), intent(inout) :: proj_set
            integer, intent(inout) :: nprev
            real(dp), dimension(lv), intent(in) :: h1, h2, h2inv
            procedure(abstract_pressure_matvec) :: matvec
            integer, intent(out) :: ierr
            ! internal
            integer, parameter :: ltot2 = lx2*ly2*lz2*lelv
            real(dp) :: pbar(ltot2), pnew(ltot2)
            real(dp) :: alpha(mxprev), work(mxprev)
            common /orthox/ pbar, pnew
            common /orthos/ alpha, work
            integer :: ntot2, i, nprev1
            real(dp) :: alphad
            real(dp), external :: vlsc2, glsc2
 
            ierr  = 0
            ntot2 = lx2*ly2*lz2*nelv

            ! E-matvec
            call matvec(pnew, proj_set(1,nprev), h1, h2, h2inv, 1)

            ! proj_set(:,nprev) -= proj_set(1:nprev-1) @ (proj_set(1:nprev-1)^T @ pnew)
            alphad = glsc2(pnew, proj_set(1,nprev), ntot2)
            nprev1 = nprev - 1
            do i = 1, nprev1
               alpha(i) = vlsc2(pnew, proj_set(1,i), ntot2)
            end do
            if (nprev1 > 0) call gop(alpha, work, '+  ', nprev1)
            do i = 1, nprev1
               call add2s2(proj_set(1,nprev), proj_set(1,i), -alpha(i), ntot2)
               alphad = alphad - alpha(i)**2
            end do

            ! Identify linear dependence (= breakdown)
            if (alphad <= 0) then
               if (alphad < -1e-12) write(6,*) 'ERROR: alphad <= 0 in econj_pressure_2Dh', alphad, nprev
               ierr = 1
               return
            end if

            ! Normalize
            call cmult(proj_set(1,nprev), 1.0_dp/sqrt(alphad), ntot2)
         end subroutine econj_pressure_2Dh
       
      end module neklab_2Dh_common
