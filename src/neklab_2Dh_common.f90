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

         public :: abstract_helmholtz_matvec, abstract_pressure_matvec

         public :: build_wmask
         public :: check_2Dh_properties
         public :: check_2Dh_masks

         public ::     map12_2Dh
         public ::     map12T_2Dh
         public :: bm2_map12_2Dh
         public :: bm2_map12T_2Dh
         public :: solve_helmholtz_2Dh
         public :: solve_pressure_2Dh
         public :: solve_pressure_block_2Dh
         public :: abstract_pressure_block_matvec
         public ::  setrhs_pressure_2Dh
         public :: gensoln_pressure_2Dh

         abstract interface
            subroutine abstract_helmholtz_matvec(Au, u, h1, h2)
               import :: dp, lx1, ly1, lz1, lelv
               real(dp), dimension(lx1*ly1*lz1*lelv), intent(out) :: Au
               real(dp), dimension(lx1*ly1*lz1*lelv), intent(in)  :: u, h1, h2
            end subroutine abstract_helmholtz_matvec
            subroutine abstract_pressure_matvec(Ap, p, h1, h2, h2inv, intype)
               import :: dp, lx1, ly1, lz1, lx2, ly2, lz2, lelv
               real(dp), dimension(lx2*ly2*lz2*lelv), intent(out) :: Ap
               real(dp), dimension(lx2*ly2*lz2*lelv), intent(in)  :: p
               real(dp), dimension(lx1*ly1*lz1*lelv), intent(in)  :: h1, h2, h2inv
               integer, intent(in) :: intype
            end subroutine abstract_pressure_matvec
            subroutine abstract_pressure_block_matvec(Ap, p, h1, h2, h2inv)
      !! Two-slot (Re, Im) pressure operator: Ap(:,k) = sum_l E_kl p(:,l).
               import :: dp, lx1, ly1, lz1, lx2, ly2, lz2, lelv
               real(dp), dimension(lx2*ly2*lz2*lelv,2), intent(out) :: Ap
               real(dp), dimension(lx2*ly2*lz2*lelv,2), intent(in)  :: p
               real(dp), dimension(lx1*ly1*lz1*lelv),   intent(in)  :: h1, h2, h2inv
            end subroutine abstract_pressure_block_matvec
         end interface

      ! Block GMRES work space (solve_pressure_block_2Dh)
         real(dp), dimension(lp,2,lgmres), private :: bv_gm, bz_gm
         real(dp), dimension(lp,2),        private :: bw_gm, br_gm, bx_gm
      ! Block residual projection (incomprp_block_solve_2Dh): E-orthonormal
      ! basis of previous two-slot solutions. Reset when bd(1)/dt or the
      ! wavenumber changes (E depends on both).
         real(dp), dimension(lp,2,mxprev), private :: bpbasis
         integer,  private :: bnprev = 0
         real(dp), private :: bwav_last = -1.0_dp
         real(dp), dimension(lp,2), private :: bpbar, bpnew
         public :: incomprp_block_solve_2Dh

      contains

         subroutine build_wmask(verbose)
            logical, optional, intent(in) :: verbose
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
                  cb = cbc(ifc, ie, 1)
                  select case (cb)
                  case ('W  ', 'v  ', 'V  ', 'vl ', 'VL ', 'mv ', 'MV ', 'ON ', 'on ')
                     call facev(wmask, ie, ifc, 0.0_dp, lx1, ly1, lz1)
                  end select
               end do
            end do

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

         subroutine setup_step_2Dh()
            if (istep == 1) nprev(:) = 0
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
            real(dp), dimension(lv), intent(out) :: resv1, resv2, resv3
            real(dp), dimension(lv), intent(in)  :: gradw_jp
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
            real(dp), dimension(lv), intent(out) :: h1, h2, h2inv
            real(dp) :: dtbd
            integer :: ntot1
            ntot1  = lx1*ly1*lz1*nelv
            ifield = 1
            dtbd   = bd(1)/dt
            if (dtbd /= dtbdlast) then
               nprev(:) = 0
               bnprev   = 0
               dtbdlast = dtbd
            end if
            call rzero  (h1, ntot1)
            call cmult2 (h2, vtrans(1,1,1,1,ifield), dtbd, ntot1)
            call invers2(h2inv, h2, ntot1)
         end subroutine setup_incomprp_2Dh

         subroutine ebar_2Dh(ebar, matvec, h1, h2, h2inv, wavenumber)
            real(dp), intent(out) :: ebar
            procedure(abstract_pressure_matvec) :: matvec
            real(dp), dimension(lv), intent(in) :: h1, h2, h2inv
            real(dp), intent(in) :: wavenumber
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
            real(dp), dimension(lp), intent(inout) :: dpr_jp
            real(dp), dimension(lv), intent(in)    :: h1, h2, h2inv
            procedure(abstract_pressure_matvec) :: matvec
            real(dp), intent(in) :: ebar, wavenumber
            integer, intent(in)  :: maxit
            character(len=*), intent(in) :: label
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

         subroutine map12_2Dh(y_m2, x_m1, ie)
            real(dp), dimension(lx2,ly2), intent(out) :: y_m2
            real(dp), dimension(lx1,ly1), intent(in)  :: x_m1
            integer, intent(in) :: ie
            real(dp), dimension(lx2,ly1) :: ta
            call mxm(ixm12, lx2, x_m1, lx1, ta, ly1)
            if (ifaxis) then
               if (ifrzer(ie)) then
                  call mxm(ta, lx2, iatm12, ly1, y_m2, ly2)
               else
                  call mxm(ta, lx2, ictm12, ly1, y_m2, ly2)
               end if
            else
               call mxm(ta, lx2, iytm12, ly1, y_m2, ly2)
            end if
         end subroutine map12_2Dh

         subroutine map12T_2Dh(y_m1, x_m2, ie)
            real(dp), dimension(lx1,ly1), intent(out) :: y_m1
            real(dp), dimension(lx2,ly2), intent(in)  :: x_m2
            integer, intent(in) :: ie
            real(dp), dimension(lx1,ly2) :: ta
            call mxm(ixtm12, lx1, x_m2, lx2, ta, ly2)
            if (ifaxis) then
               if (ifrzer(ie)) then
                  call mxm(ta, lx1, iam12, ly2, y_m1, ly1)
               else
                  call mxm(ta, lx1, icm12, ly2, y_m1, ly1)
               end if
            else
               call mxm(ta, lx1, iym12, ly2, y_m1, ly1)
            end if
         end subroutine map12T_2Dh

         subroutine bm2_map12_2Dh(y_m2, x_m1, const_fld)
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

         subroutine solve_helmholtz_2Dh(x, f, h1, h2, mask, mult, imsh, tin, maxit, isd,
     &                                  binv, name, matvec)
            include 'FDMH1'
            real(dp), dimension(lv), intent(out) :: x
            real(dp), dimension(lv), intent(in)  :: f, h1, h2, mask, mult, binv
            integer, intent(in) :: imsh, maxit, isd
            real(dp), intent(in) :: tin
            character(len=7), intent(in) :: name
            procedure(abstract_helmholtz_matvec) :: matvec
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
            call rzero(diagt, maxcg)
            call rzero(upper, maxcg)
            rho  = 0.0_dp
            nxyz = lx1*ly1*lz1
            nel  = nelv
            vol  = volvm1
            if (imsh == 2) nel = nelt
            if (imsh == 2) vol = voltm1
            n = nel*nxyz

            tol = abs(tin)
            if (restol(ifield) /= 0) tol = abs(restol(ifield))
            if (tin < 0) tol = abs(tin)
            niter = min(maxit, maxcg)

            if (.not. ifsolv) then
               call setfast(h1, h2, imsh)
               ifsolv = .true.
            end if

            call setprec(d, h1, h2, imsh, isd)

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
               call col3(z, r, d, n)
               rtz2 = rtz1
               scalar(1) = vlsc3(z, r, mult, n)
               if (param(18) == 1) then
                  scalar(2) = vlsc3(r, r, mult, n)
               else
                  scalar(2) = vlsc32(r, mult, binv, n)
               end if
               call gop(scalar, w, '+  ', 2)
               rtz1 = scalar(1)
               rbn2 = sqrt(scalar(2)/vol)
               if (iter == 1) rbn0 = rbn2
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

               beta = rtz1/rtz2
               if (iter == 1) beta = 0.0_dp
               call add2s1(p, z, beta, n)
               call matvec(w, p, h1, h2)
               call dssum (w, lx1, ly1, lz1)
               call col2  (w, mask, n)
               rho0  = rho
               rho   = glsc3(w, p, mult, n)
               alpha = rtz1/rho
               alphm = -alpha
               call add2s2(x, p, alpha, n)
               call add2s2(r, w, alphm, n)
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
                  call col3(r_gmres, ml_gmres, res, ntot2)
               else
                  call copy(r_gmres, res, ntot2)
                  call matvec(w_gmres, x_gmres, h1, h2, h2inv, intype)
                  call add2s2(r_gmres, w_gmres, -1.0_dp, ntot2)
                  call col2(r_gmres, ml_gmres, ntot2)
               end if
               gamma_gmres(1) = sqrt(glsc2(r_gmres, r_gmres, ntot2))
               if (iter == 0) then
                  div0 = gamma_gmres(1)*norm_fac
                  if (param(21) < 0) tolpss = abs(param(21))*div0
               end if
               rnorm = 0.0_dp
               if (gamma_gmres(1) == 0.0_dp) goto 9000
               temp = 1.0_dp/gamma_gmres(1)
               call cmult2(v_gmres(1,1), r_gmres, temp, ntot2)

               do j = 1, m
                  iter = iter + 1
                  call col3(w_gmres, mu_gmres, v_gmres(1,j), ntot2)
                  etime2 = dnekclock()
                  call neklab_timer_start(t_pres_prec)
                  call hsmg_solve(z_gmres(1,j), w_gmres)
                  call neklab_timer_stop(t_pres_prec)
                  etime_p = etime_p + dnekclock() - etime2
                  if (ebar /= 0.0_dp) call cadd(z_gmres(1,j), glsum(w_gmres,ntot2)/ebar, ntot2)
                  call matvec(w_gmres, z_gmres(1,j), h1, h2, h2inv, intype)
                  call col2(w_gmres, ml_gmres, ntot2)
                  do i = 1, j
                     h_gmres(i,j) = vlsc2(w_gmres, v_gmres(1,i), ntot2)
                  end do
                  call gop(h_gmres(1,j), wk1, '+  ', j)
                  do i = 1, j
                     call add2s2(w_gmres, v_gmres(1,i), -h_gmres(i,j), ntot2)
                  end do
                  do i = 1, j-1
                     temp = h_gmres(i,j)
                     h_gmres(i  ,j) =  c_gmres(i)*temp + s_gmres(i)*h_gmres(i+1,j)
                     h_gmres(i+1,j) = -s_gmres(i)*temp + c_gmres(i)*h_gmres(i+1,j)
                  end do
                  alpha = sqrt(glsc2(w_gmres, w_gmres, ntot2))
                  rnorm = 0.0_dp
                  if (alpha == 0.0_dp) goto 900
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
                  if (rnorm < tolpss) goto 900
                  if (j == m) goto 1000
                  temp = 1.0_dp/alpha
                  call cmult2(v_gmres(1,j+1), w_gmres, temp, ntot2)
               end do
  900          iconv = 1
 1000          continue
               do k = j, 1, -1
                  temp = gamma_gmres(k)
                  do i = j, k+1, -1
                     temp = temp - h_gmres(k,i)*c_gmres(i)
                  end do
                  c_gmres(k) = temp/h_gmres(k,k)
               end do
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

         subroutine solve_pressure_block_2Dh(res, h1, h2, h2inv, matvec, iter, ebar, maxit, label, wavenumber)
      !! GMRES(lgmres) for the two-slot pressure system E dp = res, E symmetric
      !! positive definite in the real (Re, Im) form [[A, -K], [K, A]]. Same
      !! algorithm, scaling (ml/mu, norm_fac) and tolerance logic as
      !! solve_pressure_2Dh; the vectors are (lp, 2), inner products sum over both
      !! slots, and the preconditioner is block-diagonal: hsmg + ebar correction on
      !! each slot (it is exact for the slot-diagonal part A). No ortho
      !! (alpha /= 0: E is non-singular). Residual projection: see
      !! incomprp_block_solve_2Dh, which wraps this solver.
            include 'GMRES'
            real(dp), dimension(lp,2), intent(inout) :: res
            real(dp), dimension(lv),   intent(in)    :: h1, h2, h2inv
            procedure(abstract_pressure_block_matvec) :: matvec
            integer, intent(in)  :: maxit
            integer, intent(out) :: iter
            real(dp), intent(in) :: ebar, wavenumber
            character(len=*), intent(in) :: label
            character(len=*), parameter :: this_procedure = 'solve_pressure_block_2Dh'
            real(dp) :: divex
            common  /ctolpr/ divex
            logical          ifprint
            common  /cprint/ ifprint
            real(dp), dimension(lgmres) :: wk1
            real(dp) :: alpha, l, temp, div0, ratio, rnorm, tolpss, tol1, tol2
            real(dp) :: etime_p, etime1, etime2, dnekclock
            integer  :: i, j, k, m, iconv, ntot2, kk
            real(dp), save :: norm_fac
            logical,  save :: iflag = .false.
            real(dp), external :: vlsc2, glsum
            character(len=128) :: msg

            call neklab_timer_start(t_pres_solve)
            ntot2 = lx2*ly2*lz2*nelv
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

            call chktcg2(tol1, res(1,1), iconv)
            call chktcg2(tol2, res(1,2), iconv)
            tolps = min(tol1, tol2)
            if (param(21) > 0 .and. tolps > abs(param(21))) tolps = abs(param(21))
            if (istep == 0) tolps = 1.e-4
            tolpss = tolps

            iconv = 0
            call rzero(bx_gm, 2*lp)
            div0  = 0.0_dp
            rnorm = 0.0_dp

            do while (iconv == 0 .and. iter < maxit)
               if (iter == 0) then
                  do kk = 1, 2
                     call col3(br_gm(1,kk), ml_gmres, res(1,kk), ntot2)
                  end do
               else
                  call matvec(bw_gm, bx_gm, h1, h2, h2inv)
                  do kk = 1, 2
                     call copy  (br_gm(1,kk), res(1,kk), ntot2)
                     call add2s2(br_gm(1,kk), bw_gm(1,kk), -1.0_dp, ntot2)
                     call col2  (br_gm(1,kk), ml_gmres, ntot2)
                  end do
               end if
               gamma_gmres(1) = sqrt(bdot(br_gm, br_gm))
               if (iter == 0) then
                  div0 = gamma_gmres(1)*norm_fac
                  if (param(21) < 0) tolpss = abs(param(21))*div0
               end if
               rnorm = 0.0_dp
               if (gamma_gmres(1) == 0.0_dp) goto 9000
               temp = 1.0_dp/gamma_gmres(1)
               do kk = 1, 2
                  call cmult2(bv_gm(1,kk,1), br_gm(1,kk), temp, ntot2)
               end do

               do j = 1, m
                  iter = iter + 1
                  etime2 = dnekclock()
                  call neklab_timer_start(t_pres_prec)
                  do kk = 1, 2
                     call col3(bw_gm(1,kk), mu_gmres, bv_gm(1,kk,j), ntot2)
                     call hsmg_solve(bz_gm(1,kk,j), bw_gm(1,kk))
                     if (ebar /= 0.0_dp) call cadd(bz_gm(1,kk,j), glsum(bw_gm(1,kk),ntot2)/ebar, ntot2)
                  end do
                  call neklab_timer_stop(t_pres_prec)
                  etime_p = etime_p + dnekclock() - etime2
                  call matvec(bw_gm, bz_gm(1,1,j), h1, h2, h2inv)
                  do kk = 1, 2
                     call col2(bw_gm(1,kk), ml_gmres, ntot2)
                  end do
                  do i = 1, j
                     h_gmres(i,j) = vlsc2(bw_gm(1,1), bv_gm(1,1,i), ntot2)
     &                            + vlsc2(bw_gm(1,2), bv_gm(1,2,i), ntot2)
                  end do
                  call gop(h_gmres(1,j), wk1, '+  ', j)
                  do i = 1, j
                     do kk = 1, 2
                        call add2s2(bw_gm(1,kk), bv_gm(1,kk,i), -h_gmres(i,j), ntot2)
                     end do
                  end do
                  do i = 1, j-1
                     temp = h_gmres(i,j)
                     h_gmres(i  ,j) =  c_gmres(i)*temp + s_gmres(i)*h_gmres(i+1,j)
                     h_gmres(i+1,j) = -s_gmres(i)*temp + c_gmres(i)*h_gmres(i+1,j)
                  end do
                  alpha = sqrt(bdot(bw_gm, bw_gm))
                  rnorm = 0.0_dp
                  if (alpha == 0.0_dp) goto 900
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
   66             format(i5,1p4e12.5,i8,' Divergence (block)')
                  if (rnorm < tolpss) goto 900
                  if (j == m) goto 1000
                  temp = 1.0_dp/alpha
                  do kk = 1, 2
                     call cmult2(bv_gm(1,kk,j+1), bw_gm(1,kk), temp, ntot2)
                  end do
               end do
  900          iconv = 1
 1000          continue
               do k = j, 1, -1
                  temp = gamma_gmres(k)
                  do i = j, k+1, -1
                     temp = temp - h_gmres(k,i)*c_gmres(i)
                  end do
                  c_gmres(k) = temp/h_gmres(k,k)
               end do
               do i = 1, j
                  do kk = 1, 2
                     call add2s2(bx_gm(1,kk), bz_gm(1,kk,i), c_gmres(i), ntot2)
                  end do
               end do
            end do
 9000       continue
            divex = rnorm

            if (iconv == 0 .and. rnorm >= tolpss) then
               write(msg,'(A,I6,A,1pE12.4)') 'block pressure GMRES not converged after ', iter,
     &              ' iterations, residual ', rnorm
               call nek_log_warning(msg, this_module, this_procedure)
            end if

            do kk = 1, 2
               call copy(res(1,kk), bx_gm(1,kk), ntot2)
            end do

            etime1 = dnekclock() - etime1
            if (nio == 0) write(6,9998) istep, '  U-PRES gmres  '//label//' = ', wavenumber,
     &         'BL', iter, divex, div0, tolpss, etime_p, etime1
 9998       format(i11,a,F5.1,1X,A,1X,I6,1p5e13.4)
            call neklab_timer_stop(t_pres_solve)
         contains
            real(dp) function bdot(a, b)
               real(dp), dimension(lp,2), intent(in) :: a, b
               real(dp), external :: glsc2
               bdot = glsc2(a(1,1), b(1,1), ntot2) + glsc2(a(1,2), b(1,2), ntot2)
            end function bdot
         end subroutine solve_pressure_block_2Dh

      !--------------------------------------------------------------------
      !     Block (two-slot) residual projection
      !--------------------------------------------------------------------
      ! Same algorithm as setrhs/gensoln/econj_pressure_2Dh (Fischer's
      ! projection, Nek's param(93) basis size, param(95) start step) on (lp,2)
      ! vectors, with the E-inner product of the coupled operator:
      !   <a, b>_E = sum_k a(:,k).(E b)(:,k).
      ! E is symmetric positive definite at alpha /= 0 (U10a/b), so the basis is
      ! well defined. The projected part is exact; GMRES solves for the rest.

         subroutine incomprp_block_solve_2Dh(dpr, h1, h2, h2inv, matvec, ebar, wavenumber, maxit, label)
      !! Coupled pressure solve of both slots with residual projection:
      !! dpr holds the right-hand sides on entry, the corrections on exit.
            real(dp), dimension(lp,2), intent(inout) :: dpr
            real(dp), dimension(lv),   intent(in)    :: h1, h2, h2inv
            procedure(abstract_pressure_block_matvec) :: matvec
            real(dp), intent(in) :: ebar, wavenumber
            integer,  intent(in) :: maxit
            character(len=*), intent(in) :: label
            logical :: ifprjp
            integer :: istart, iter
            if (wavenumber /= bwav_last) then
               bnprev    = 0
               bwav_last = wavenumber
            end if
            istart = param(95)
            ifprjp = (istep >= istart .and. istart /= 0)
            if (ifprjp) call setrhs_pressure_block_2Dh(dpr, h1, h2, h2inv, matvec)
            call solve_pressure_block_2Dh(dpr, h1, h2, h2inv, matvec, iter, ebar, maxit,
     &                                    label, wavenumber)
            if (ifprjp) call gensoln_pressure_block_2Dh(dpr, h1, h2, h2inv, matvec)
         end subroutine incomprp_block_solve_2Dh

         subroutine setrhs_pressure_block_2Dh(dpr, h1, h2, h2inv, matvec)
      !! pbar = sum_i <dpr, q_i> q_i (q_i E-orthonormal), dpr <- dpr - E pbar.
            real(dp), dimension(lp,2), intent(inout) :: dpr
            real(dp), dimension(lv),   intent(in)    :: h1, h2, h2inv
            procedure(abstract_pressure_block_matvec) :: matvec
            real(dp) :: alpha(mxprev), work(mxprev)
            real(dp) :: alpha1, alpha2, ratio
            real(dp), external :: glsc3, vlsc2
            integer :: ntot2, i, k
            call neklab_timer_start(t_pres_proj)
            ntot2 = lx2*ly2*lz2*nelv
            call rzero(bpbar, 2*lp)
            if (bnprev == 0) then
               call neklab_timer_stop(t_pres_proj)
               return
            end if
            alpha1 = glsc3(dpr(1,1), dpr(1,1), bm2inv, ntot2) + glsc3(dpr(1,2), dpr(1,2), bm2inv, ntot2)
            if (alpha1 <= 0.0_dp) then
               call neklab_timer_stop(t_pres_proj)
               return
            end if
            alpha1 = sqrt(alpha1/volvm2)
            do i = 1, bnprev
               alpha(i) = vlsc2(dpr(1,1), bpbasis(1,1,i), ntot2) + vlsc2(dpr(1,2), bpbasis(1,2,i), ntot2)
            end do
            call gop(alpha, work, '+  ', bnprev)
            do i = 1, bnprev
               do k = 1, 2
                  call add2s2(bpbar(1,k), bpbasis(1,k,i), alpha(i), ntot2)
               end do
            end do
            call matvec(bpnew, bpbar, h1, h2, h2inv)
            do k = 1, 2
               call sub2(dpr(1,k), bpnew(1,k), ntot2)
            end do
            alpha2 = glsc3(dpr(1,1), dpr(1,1), bm2inv, ntot2) + glsc3(dpr(1,2), dpr(1,2), bm2inv, ntot2)
            if (alpha2 > 0.0_dp) then
               alpha2 = sqrt(alpha2/volvm2)
               ratio  = alpha1/alpha2
               if (nio == 0) write(6,13) istep, '  Project PRES BL', alpha2, alpha1, ratio, bnprev, mxprev
            end if
 13         format(i11,a,6x,1p3e13.4,i4,i4)
            call neklab_timer_stop(t_pres_proj)
         end subroutine setrhs_pressure_block_2Dh

         subroutine gensoln_pressure_block_2Dh(dpr, h1, h2, h2inv, matvec)
      !! dpr <- dpr + pbar; add the new solution to the basis (restart when full).
            real(dp), dimension(lp,2), intent(inout) :: dpr
            real(dp), dimension(lv),   intent(in)    :: h1, h2, h2inv
            procedure(abstract_pressure_block_matvec) :: matvec
            integer :: ntot2, mprev, ierr, k
            call neklab_timer_start(t_pres_proj)
            mprev = min(int(param(93)), mxprev)
            ntot2 = lx2*ly2*lz2*nelv
            if (bnprev < mprev) then
               bnprev = bnprev + 1
               do k = 1, 2
                  call copy(bpbasis(1,k,bnprev), dpr(1,k), ntot2)
                  call add2(dpr(1,k), bpbar(1,k), ntot2)
               end do
               call econj_pressure_block_2Dh(h1, h2, h2inv, matvec, ierr)
               if (ierr == 1) then
                  bnprev = 1
                  do k = 1, 2
                     call copy(bpbasis(1,k,1), dpr(1,k), ntot2)
                  end do
                  call econj_pressure_block_2Dh(h1, h2, h2inv, matvec, ierr)
                  if (ierr == 1) bnprev = 0
               end if
            else
               bnprev = 1
               do k = 1, 2
                  call add2(dpr(1,k), bpbar(1,k), ntot2)
                  call copy(bpbasis(1,k,1), dpr(1,k), ntot2)
               end do
               call econj_pressure_block_2Dh(h1, h2, h2inv, matvec, ierr)
               if (ierr == 1) bnprev = 0
            end if
            call neklab_timer_stop(t_pres_proj)
         end subroutine gensoln_pressure_block_2Dh

         subroutine econj_pressure_block_2Dh(h1, h2, h2inv, matvec, ierr)
      !! E-orthonormalise basis vector bnprev against 1..bnprev-1 (classical
      !! Gram-Schmidt, as econj_pressure_2Dh).
            real(dp), dimension(lv), intent(in) :: h1, h2, h2inv
            procedure(abstract_pressure_block_matvec) :: matvec
            integer, intent(out) :: ierr
            real(dp) :: alpha(mxprev), work(mxprev), alphad
            real(dp), external :: vlsc2, glsc2
            integer :: ntot2, i, k, n1
            ierr  = 0
            ntot2 = lx2*ly2*lz2*nelv
            call matvec(bpnew, bpbasis(1,1,bnprev), h1, h2, h2inv)
            alphad = glsc2(bpnew(1,1), bpbasis(1,1,bnprev), ntot2)
     &             + glsc2(bpnew(1,2), bpbasis(1,2,bnprev), ntot2)
            n1 = bnprev - 1
            do i = 1, n1
               alpha(i) = vlsc2(bpnew(1,1), bpbasis(1,1,i), ntot2) + vlsc2(bpnew(1,2), bpbasis(1,2,i), ntot2)
            end do
            if (n1 > 0) call gop(alpha, work, '+  ', n1)
            do i = 1, n1
               do k = 1, 2
                  call add2s2(bpbasis(1,k,bnprev), bpbasis(1,k,i), -alpha(i), ntot2)
               end do
               alphad = alphad - alpha(i)**2
            end do
            if (alphad <= 0.0_dp) then
               if (alphad < -1e-12 .and. nio == 0)
     &            write(6,*) 'WARNING: alphad <= 0 in econj_pressure_block_2Dh', alphad, bnprev
               ierr = 1
               return
            end if
            do k = 1, 2
               call cmult(bpbasis(1,k,bnprev), 1.0_dp/sqrt(alphad), ntot2)
            end do
         end subroutine econj_pressure_block_2Dh

         subroutine setrhs_pressure_2Dh(dpr, h1, h2, h2inv, proj_set, niprev, matvec, info_str)
            real(dp), dimension(lp), intent(inout) :: dpr
            real(dp), dimension(lv), intent(in)    :: h1, h2, h2inv
            real(dp), dimension(lp,mxprev), intent(inout) :: proj_set
            integer, intent(inout) :: niprev
            procedure(abstract_pressure_matvec) :: matvec
            character(len=2), intent(in) :: info_str
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
            do i = 1, niprev
               alpha(i) = vlsc2(dpr, proj_set(1,i), ntot2)
            end do
            call gop(alpha, work, '+  ', niprev)
            do i = 1, niprev
               call add2s2(pbar, proj_set(1,i), alpha(i), ntot2)
            end do
            call matvec(pnew, pbar, h1, h2, h2inv, 1)
            call sub2(dpr, pnew, ntot2)
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
            integer, parameter :: ltot2 = lx2*ly2*lz2*lelv
            real(dp) :: pbar(ltot2), pnew(ltot2)
            real(dp) :: alpha(mxprev), work(mxprev)
            common /orthox/ pbar, pnew
            common /orthos/ alpha, work
            integer :: ntot2, mprev, ierr

            call neklab_timer_start(t_pres_proj)
            mprev = min(int(param(93)), mxprev)
            ntot2 = lx2*ly2*lz2*nelv
            if (niprev < mprev) then
               niprev = niprev + 1
               call copy(proj_set(1,niprev), dpr, ntot2)
               call add2(dpr, pbar, ntot2)
               call econj_pressure_2Dh(proj_set, niprev, h1, h2, h2inv, matvec, ierr)
               if (ierr == 1) then
                  niprev = 1
                  call copy(proj_set(1,niprev), dpr, ntot2)
                  call econj_pressure_2Dh(proj_set, niprev, h1, h2, h2inv, matvec, ierr)
               end if
            else
               niprev = 1
               call add2(dpr, pbar, ntot2)
               call copy(proj_set(1,niprev), dpr, ntot2)
               call econj_pressure_2Dh(proj_set, niprev, h1, h2, h2inv, matvec, ierr)
            end if
            call neklab_timer_stop(t_pres_proj)
         end subroutine gensoln_pressure_2Dh

         subroutine econj_pressure_2Dh(proj_set, nprev, h1, h2, h2inv, matvec, ierr)
            real(dp), dimension(lp,mxprev), intent(inout) :: proj_set
            integer, intent(inout) :: nprev
            real(dp), dimension(lv), intent(in) :: h1, h2, h2inv
            procedure(abstract_pressure_matvec) :: matvec
            integer, intent(out) :: ierr
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
            call matvec(pnew, proj_set(1,nprev), h1, h2, h2inv, 1)
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
            if (alphad <= 0) then
               if (alphad < -1e-12) write(6,*) 'ERROR: alphad <= 0 in econj_pressure_2Dh', alphad, nprev
               ierr = 1
               return
            end if
            call cmult(proj_set(1,nprev), 1.0_dp/sqrt(alphad), ntot2)
         end subroutine econj_pressure_2Dh

      end module neklab_2Dh_common