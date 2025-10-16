      module neklab_2Dh
         use LightKrylov, only: dp
         use neklab_nek_setup, only: nek_log_debug, nek_log_message, nek_stop_error
         implicit none
         include "SIZE"
         include "TOTAL"
      
         private
         character(len=*), parameter, private :: this_module = 'neklab_2Dh'
      
         public :: nek_advance_2Dh
         public :: pressure_matvec_2Dh, helmholtz_matvec_2Dh
         public :: solve_pressure_2Dh, solve_helmholtz_2Dh
         public :: compute_dw, compute_frc_div, compute_gradz_p
      
      contains

         subroutine nek_advance_2dh(beta_z)
            !implicit none
            !include 'SIZE'
            !include 'SOLN' ! jp
            !include 'MASS' ! bm1
            real(dp), intent(in) :: beta_z
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
            real(dp), dimension(lx1*ly1*lz1*lelv,lpert) :: gradz_p
            ! spanwise pressure gradient
            real(dp), dimension(lx2,ly2,lz2,lelv,lpert) :: frc_div
            ! spanwise divergence
            real(dp), dimension(lx2*ly2*lz2*lelv) :: prextr
            ! extrapolated pressure
            real(dp), dimension(lx1*ly1*lz1*lelv,lpert) :: dw
            ! spanwise velocity correction
            real(dp), dimension(lx2*ly2*lz2*lelv,lpert) :: dpr
            ! pressure correction
            real(dp), dimension(lx1*ly1*lz1*lelv,lpert) :: bm1h2inv
            ! prefactor
            ! Miscellaneous
            character(len=*), parameter :: this_procedure = 'nek_advance'
            character(len=*), parameter :: fmt = '(A,I6)'
            character(len=256) :: msg
            real(dp) :: dtbd
            integer :: igeom, iter, intype, kfldfdm
            integer :: ntot1, ntot2

            ntot1 = lx1*ly1*lz1*nelv
            ntot2 = lx2*ly2*lz2*nelv
            write(msg,fmt) 'Step', istep
            call nek_log_debug(msg,this_module,this_procedure)

            call nekgsync
            call setup_convect(2)
            call settime
            call setsolv
            call comment
            call setprop

            ! compute best estimate for d/dz (p)
            
            msg = 'Compute explicit pressure gradient term for w'
            call nek_log_debug(msg,this_module,this_procedure)
            ifield = 1
            do jp = 1, npert
               ! we need to do this outside of the other jp loop for extrapolation of the correct pressure term
               call extrapprp(prextr)
               call compute_gradz_p(gradz_p, prextr, beta_z)
            end do

            ! Compute intermediate velocities

            msg = 'Solve momentum equations for u/v/w'
            call nek_log_debug(msg,this_module,this_procedure)
            do jp = 1, npert
               do igeom = 1,ngeom
                  if (igeom == 1) then
                  !
                  !  Old geometry, old velocity
                  !
                     ifield = 1     !  velocity
                     call makefp
                     call lagfieldp
                     ifield = 2     !  temperature (w)
                     call makeqp
                     call lagscalp
                     !
                  else
                  !
                  !  New geometry, new velocity
                  !
                     intype = -1
                     ifield = 1     !  velocity
                     call sethlm   (h1,h2,intype)
                     ! cresvipp 
                     call bcdirvc (vxp(1,jp), vyp(1,jp), vzp(1,jp),v1mask,v2mask,v3mask)
                     call bcdirvc (tp(1,1,jp),vyp(1,jp), vzp(1,jp),v3mask,v2mask,v3mask)

                     call extrapprp(prextr)
                     call opgradt(resv1,resv2,resv3,prextr)      ! d/dx pr, d/dy pr
                     call opadd2 (resv1,resv2,resv3,bfxp(1,jp), bfyp(1,jp),bfzp(1,jp))
                     call copy   (resv3,gradz_p(1,jp),ntot1)     ! d/dz pr
                     call   add2 (resv3, bqp(1,1,jp), ntot1)
                     ! ophx
                     call helmholtz_matvec_2Dh(w1, vxp(1,jp),h1,h2,beta_z)
                     call helmholtz_matvec_2Dh(w2, vyp(1,jp),h1,h2,beta_z)
                     call helmholtz_matvec_2Dh(w3,tp(1,1,jp),h1,h2,beta_z)
                     call opsub2(resv1,resv2,resv3,w1,w2,w3)
                     call   sub2(            resv3,      w3,ntot1)
                
                     ! ophinv
                     kfldfdm = -1
                     tolhv = abs(tolhv)
                     ifield = 1 !  velocity

                     call dssum  (resv1,lx1,ly1,lz1)
                     call col2   (resv1,v1mask,ntot1)
                     if (istep < 10) call chktcg1 (tolhv,resv1,h1,h2,v1mask,vmult,imesh,1)
                     call solve_helmholtz_2Dh(dv1,resv1,h1,h2,v1mask,vmult,imesh,tolhv,nmxv,1,binvm1,'VELX',beta_z)

                     call dssum  (resv2,lx1,ly1,lz1)
                     call col2   (resv2,v2mask,ntot1)
                     if (istep < 10) call chktcg1 (tolhv,resv2,h1,h2,v2mask,vmult,imesh,2)
                     call solve_helmholtz_2Dh(dv2,resv2,h1,h2,v2mask,vmult,imesh,tolhv,nmxv,2,binvm1,'VELY',beta_z)
                     
                     call dssum  (resv3,lx1,ly1,lz1)
                     call col2   (resv3,v3mask,ntot1)
                     if (istep < 10) call chktcg1 (tolhv,resv3,h1,h2,v3mask,vmult,imesh,3)
                     call solve_helmholtz_2Dh(dv3,resv3,h1,h2,v3mask,vmult,imesh,tolhv,nmxv,3,binvm1,'VELZ',beta_z)
                     
                     call opadd2 (vxp(1,jp),vyp(1,jp),vzp(1,jp),dv1,dv2,dv3)
                     call add2   (tp(1,1,jp),dv3,ntot1)
                  end if
               end do ! igeom
            end do ! jp
                  
            msg = 'Compute pressure correction to enforce mass balance'
            call nek_log_debug(msg,this_module,this_procedure)
            do jp = 1, npert
               ! incomprp
               ifield = 1 ! velocity
               intype = 1
               dtbd   = bd(1)/dt
               
               call rzero   (h1,ntot1)
               call cmult2  (h2,vtrans(1,1,1,1,ifield),dtbd,ntot1)
               call invers2 (h2inv,h2,ntot1)
                  
               call opdiv   (dpr(1,jp),vxp(1,jp),vyp(1,jp),vzp(1,jp))

               call compute_frc_div(frc_div, tp, beta_z)
               call add2(dpr(1,jp),frc_div(1,1,1,1,jp),ntot2)
                  
               call chsign  (dpr(1,jp),ntot2)
               call ortho   (dpr(1,jp))
             
               call solve_pressure_2Dh(dpr(1,jp),h1,h2,h2inv,beta_z,intype,iter)
            end do ! jp                

            ! Reconstruct pressure and add pressure correction

            msg = 'Compute velocity correction based on pressure correction'
            call nek_log_debug(msg,this_module,this_procedure)
            ! prepare
            call invcol3 (bm1h2inv,bm1,h2inv,ntot1)    ! = bm1/h2inv = B*h2
            call dssum   (bm1h2inv,lx1,ly1,lz1)

            do jp = 1, npert
               call opgradt (w1 ,w2 ,w3 ,dpr(1,jp))
               call opbinv  (dv1,dv2,dv3,w1 ,w2 ,w3 ,h2inv)
               call compute_dw(dw, dpr, bm1h2inv, beta_z) ! dw

               call opadd2 (vxp(1,jp),vyp(1,jp),vzp(1,jp), dv1,dv2,dv3)
               call   add2 (tp(1,1,jp), dw(1,jp), ntot1)
             
               call extrapprp(prextr)
               call lagpresp
               call add3(prp(1,jp),prextr,dpr(1,jp), ntot2)
            end do ! jp
         end subroutine nek_advance_2Dh

         subroutine compute_dw(dw, dpr, bm1h2inv, beta_z)
            implicit none
            include 'SIZE'
            include 'SOLN' ! v3mask, jp
            include 'MASS' ! bm1
            real, dimension(lx1*ly1*lz1*lelv,lpert), intent(out) :: dw
            real, dimension(lx2*ly2*lz2*lelv,lpert), intent(in) :: dpr
            real, dimension(lx1*ly1*lz1*lelv,lpert), intent(in) :: bm1h2inv
            real, intent(in) :: beta_z
            ! internal
            integer :: ipert, spert, ntot1
            real, dimension(lx1*ly1*lz1) :: wrk1, wrk2
   
            ipert = npert + 1 - jp
            spert = merge(1, -1, jp == 1)  ! +1 for jp = 1, -1 for jp = 2
            ntot1 = lx1*ly1*lz1*nelv
            call rzero  (dw(1,jp), ntot1)
            call mappr  (dw(1,jp), dpr(1,ipert), wrk1, wrk2) ! map to vmesh
            call col2c  (dw(1,jp), bm1, -spert*beta_z, ntot1)
            ! call dssum on dw
            call col2   (dw(1,jp), v3mask, ntot1)
            call dssum  (dw(1,jp), lx1, ly1, lz1)
            call invcol2(dw(1,jp), bm1h2inv,ntot1)
         end
   
         subroutine compute_frc_div(frc_div, w, beta_z)
            implicit none
            include 'SIZE'
            include 'SOLN' ! jp
            include 'MASS' ! bm2
            real, dimension(lx2,ly2,lz2,lelv,lpert), intent(out) :: frc_div
            real, dimension(lx1*ly1*lz1*lelv,ldimt,lpert), intent(in) :: w
            real, intent(in) :: beta_z
            ! internal
            integer :: ipert, spert, ie, ie1, nxyz1, ntot2
   
            ipert = npert + 1 - jp
            spert = merge(1, -1, jp == 1)  ! +1 for jp = 1, -1 for jp = 2
            nxyz1 = lx1*ly1*lz1
            ntot2 = lx2*ly2*lz2*nelv
            call rzero(frc_div(1,1,1,1,jp), ntot2)
            do ie = 1, nelv
               ie1 = (ie-1)*nxyz1+1
               call map12 (frc_div(1,1,1,ie,jp), w(ie1,1,ipert), ie) ! map to pmesh
            end do
            call col2c(frc_div(1,1,1,1,jp), bm2, spert*beta_z, ntot2) ! this is to be consistent with opdiv output
         end
   
         subroutine compute_gradz_p(gradz_p, prextr, beta_z)
            implicit none
            include 'SIZE'
            include 'SOLN' ! jp
            include 'MASS' ! bm1
            real, dimension(lx1*ly1*lz1*lelv,lpert), intent(out) :: gradz_p
            real, dimension(lx2*ly2*lz2*lelv), intent(in) :: prextr
            real, intent(in) :: beta_z
            ! internal
            integer :: ipert, spert, ntot1
            real, dimension(lx1,ly1,lz1,lelv) :: wrk
            real, dimension(lx1*ly1*lz1) :: wrk1, wrk2
   
            ipert = npert + 1 - jp
            spert = merge(1, -1, jp == 1)  ! +1 for jp = 1, -1 for jp = 2
            ntot1 = lx1*ly1*lz1*nelv
            call rzero(gradz_p(1,ipert), ntot1)
            call mappr(gradz_p(1,ipert), prextr, wrk1, wrk2)  ! map to vmesh
            call col2c(gradz_p(1,ipert), bm1, spert*beta_z, ntot1)
            ! dssum will be called on the vel residual later
         end

         subroutine solve_pressure_2Dh(res,h1,h2,h2inv,beta_z,intype,iter)
            implicit none
            include 'GMRES'

            ! Solve the pressure equation by right-preconditioned GMRES iteration.
            ! intype =  0  (steady)
            ! intype =  1  (explicit)
            ! intype = -1  (implicit)

            ! This is essentially a copy of uzawa_gmres with updated matvecs to reflect the
            ! updated operator.

            real :: divex
            common  /ctolpr/ divex
            logical          ifprint
            common  /cprint/ ifprint
            real, dimension(lx2*ly2*lz2*lelv), intent(inout) :: res  
            real, dimension(lx1,ly1,lz1,lelv), intent(in)    :: h1   
            real, dimension(lx1,ly1,lz1,lelv), intent(in)    :: h2   
            real, dimension(lx1,ly1,lz1,lelv), intent(in)    :: h2inv
            integer, intent(in) :: intype
            integer, intent(inout) :: iter
            !!! add the contribution from the 3rd perturbation component
            real             beta_z
            !!!
            ! internal
            real, dimension(lx2,ly2,lz2,lelv) :: wp
            common /scrmg/   wp
            real, dimension(lgmres) :: y, wk1, wk2
            common /ctmp0/   wk1, wk2
            common /cgmres1/ y

            real alpha, l, temp, div0, ratio, rnorm, tolpss
            real :: etime2, etime_p
            integer i, j, k, m, iconv, ntot2

            logical iflag
            save    iflag
            data    iflag /.false./
            real    norm_fac
            save    norm_fac

            real*8 etime1,dnekclock

            real, external :: vlsc2, glsc2

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

            do while(iconv.eq.0.and.iter.lt.100)
            
               if(iter.eq.0) then
                                                              !      -1
                  call col3(r_gmres,ml_gmres,res,ntot2)       ! r = L  res
               
               else
                  !update residual
                  call copy(r_gmres,res,ntot2)                                ! r = res
                  call pressure_matvec_2Dh(w_gmres,x_gmres,h1,h2,h2inv,beta_z,intype) ! w = A x
                  call add2s2(r_gmres,w_gmres,-1.,ntot2)                      ! r = r - w
                                                                              !      -1
                  call col2(r_gmres,ml_gmres,ntot2)                           ! r = L   r
               endif
                                                                  !            ______
               gamma_gmres(1) = sqrt(glsc2(r_gmres,r_gmres,ntot2))! gamma  = \/ (r,r) 
                                                                  !      1
               if(iter.eq.0) then
                  div0 = gamma_gmres(1)*norm_fac
                  if (param(21).lt.0) tolpss=abs(param(21))*div0
               endif
            
               !check for lucky convergence
               rnorm = 0.
               if(gamma_gmres(1) .eq. 0.) goto 9000
               temp = 1./gamma_gmres(1)
               call cmult2(v_gmres(1,1),r_gmres,temp,ntot2)! v  = r / gamma
                                                           !  1            1
               do j=1,m
                  iter = iter+1
                                                                 !       -1
                  call col3(w_gmres,mu_gmres,v_gmres(1,j),ntot2) ! w  = U   v
                                                                 !           j
                  etime2 = dnekclock()
                                                           !       -1                                  
                  call hsmg_solve(z_gmres(1,j),w_gmres)    ! z  = M   w
                  etime_p = etime_p + dnekclock()-etime2
               
                  call pressure_matvec_2Dh(w_gmres,z_gmres(1,j),      ! w = A z
     $                                h1,h2,h2inv,beta_z,intype)      !        j
                                                               !      -1
                  call col2(w_gmres,ml_gmres,ntot2)            ! w = L   w
                                                               !      i,j  i
                  ! 2-PASS GS, 1st pass:
              
                  do i=1,j
                     h_gmres(i,j)=vlsc2(w_gmres,v_gmres(1,i),ntot2) ! h    = (w,v )
                  enddo                                             !  i,j       i
               
                  call gop(h_gmres(1,j),wk1,'+  ',j)          ! sum over P procs
               
                  do i=1,j
                     call add2s2(w_gmres,v_gmres(1,i),-h_gmres(i,j),ntot2) ! w = w - h    v
                  enddo                                                    !          i,j  i
               
                  ! apply Givens rotations to new column
                  do i=1,j-1
                     temp = h_gmres(i,j)                   
                     h_gmres(i  ,j)=  c_gmres(i)*temp + s_gmres(i)*h_gmres(i+1,j)  
                     h_gmres(i+1,j)= -s_gmres(i)*temp + c_gmres(i)*h_gmres(i+1,j)
                  enddo
                                                                    !            ______
                  alpha = sqrt(glsc2(w_gmres,w_gmres,ntot2))        ! alpha =  \/ (w,w)
                  rnorm = 0.
                  if(alpha.eq.0.) goto 900  !converged
                  l = sqrt(h_gmres(j,j)*h_gmres(j,j)+alpha*alpha)
                  temp = 1./l
                  c_gmres(j) = h_gmres(j,j) * temp
                  s_gmres(j) = alpha  * temp
                  h_gmres(j,j) = l
                  gamma_gmres(j+1) = -s_gmres(j) * gamma_gmres(j)
                  gamma_gmres(j)   =  c_gmres(j) * gamma_gmres(j)

                  rnorm = abs(gamma_gmres(j+1))*norm_fac
                  ratio = rnorm/div0
                  if (ifprint.and.nio.eq.0) write (6,66) iter,tolpss,rnorm,div0,ratio,istep
   66             format(i5,1p4e12.5,i8,' Divergence')
               
                  if (rnorm .lt. tolpss) goto 900  !converged
                  if (j.eq.m) goto 1000 !not converged, restart
                  temp = 1./alpha
                  call cmult2(v_gmres(1,j+1),w_gmres,temp,ntot2) ! v    = w / alpha
                                                           !  j+1            
               enddo
  900          iconv = 1
 1000          continue
               !back substitution
               !     -1
               !c = H   gamma
               do k=j,1,-1
                  temp = gamma_gmres(k)
                  do i=j,k+1,-1
                     temp = temp - h_gmres(k,i)*c_gmres(i)
                  enddo
                  c_gmres(k) = temp/h_gmres(k,k)
               enddo
               !sum up Arnoldi vectors
               do i=1,j
                  call add2s2(x_gmres,z_gmres(1,i),c_gmres(i),ntot2) 
                             ! x = x + c  z
                             !          i  i
               enddo
            enddo
 9000       continue

            divex = rnorm

            ! DIAGNOSTICS
             call ortho  (w_gmres) ! Orthogonalize wrt null space, if present
            ! DIAGNOSTICS
            call copy(res,x_gmres,ntot2)

            call ortho (res)  ! Orthogonalize wrt null space, if present

            etime1 = dnekclock()-etime1
            if (nio.eq.0) write(6,9999) istep,'  U-PRES gmres  ', 
     $                                  iter,divex,div0,tolpss,etime_p,etime1
 9999       format(i11,a,I6,1p5e13.4)
         
         end subroutine solve_pressure_2Dh

         subroutine pressure_matvec_2Dh(ap,wp,h1,h2,h2inv,beta_z,intype)
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
            real, intent(in) :: beta_z
            integer, intent(in) :: intype
            ! internal
            real, dimension(lx1*ly1*lz1) :: wrk1, wrk2
            real, dimension(lx1,ly1,lz1,lelv) :: wdivm1
            real, dimension(lx1,ly1,lz1,lelv) :: h2B
            real, dimension(lx2,ly2,lz2,lelv) :: wdivm2

            integer :: ie, ntot1, ntot2

            ! Regular Ax = D (h2*B)^(-1) D^T
            call cdabdtp(Ap, wp, h1, h2, h2inv, intype)
            ! add 3rd perturbation component. 
            ! we want to subtract the contribution
            !
            ! B * d/dz (h2*B)^(-1) d/dz wp = (-i beta_z) (h2)^(-1) (-i beta_z) wp = - beta_z^2 (h2)^(-1) wp
            !
            ! Note: we are only considering the ifanls = .false. case
            ntot1 = lx1*ly1*lz1*nelv
            ntot2 = lx2*ly2*lz2*nelv
            call mappr  (wdivm1, wp, wrk1, wrk2)                    ! map to vmesh
            call col2   (wdivm1, v3mask, ntot1)
            call col2c  (wdivm1, h2inv, -beta_z**2, ntot1)          ! collate -beta_z^2 * (h2)^-1
            call dssum  (wdivm1, lx1, ly1, lz1)
            do ie = 1, nelv
               call map12 (wdivm2(1,1,1,ie), wdivm1(1,1,1,ie), ie)  ! map wdiv to pmesh
            end do
            call col2(wdivm2, bm2, ntot2)
            call sub2(ap, wdivm2, ntot2)
         end subroutine pressure_matvec_2Dh

         subroutine helmholtz_matvec_2Dh(Au, u, h1, h2, beta_z)
            implicit none
            real, dimension(lx1,ly1,lz1,1), intent(in) :: Au   
            real, dimension(lx1,ly1,lz1,1), intent(in) :: u   
            real, dimension(lx1,ly1,lz1,1), intent(in) :: h1
            real, dimension(lx1,ly1,lz1,1), intent(in) :: h2
            real beta_z
            ! internal
            real, dimension(lx1,ly1,lz1,lelv) :: tmp
            integer :: imesh, ntot1
            imesh = 1
            ntot1 = lx1*ly1*lz1*nelv
            ! regular Ax
            call axhelm (Au, u, h1, h2, imesh, 1)
            ! added spanwise viscous diffusion term
            call col3 (tmp, u, vdiff(1,1,1,1,1), ntot1)
            call col2c(tmp, bm1, -beta_z**2, ntot1)
            call add2 (Au, tmp, ntot1)
         end subroutine helmholtz_matvec_2Dh

         subroutine solve_helmholtz_2Dh(x,f,h1,h2,mask,mult,imsh,tin,maxit,isd,binv,name,beta_z)
            implicit none
            include 'FDMH1'
            !-------------------------------------------------------------------------
            !
            !     Solve the Helmholtz equation, H*U = RHS,
            !     using preconditioned conjugate gradient iteration.
            !     Preconditioner: diag(H).
            !
            !------------------------------------------------------------------------
            real, intent(out) :: x(1)
            real, intent(in) :: f(1)
            real, intent(in) :: h1(1)
            real, intent(in) :: h2(1)
            real, intent(in) :: mask(1)
            real, intent(in) :: mult(1)
            real, intent(in) :: binv(1)
            integer, intent(in) :: imsh
            real, intent(in) :: tin
            integer, intent(in) :: maxit
            integer, intent(in) :: isd
            character(len=4), intent(in) :: name
            real, intent(in) :: beta_z
            ! internal 
            integer :: i, j, iter, krylov, n, nel, niter, nxyz
            real :: alpha, beta, alphm, fmax, h2max
            real :: rbn0, rbn2
            real :: rho, rho0
            real :: rmean, rtz1, rtz2
            real :: skmin, smean
            real :: tol, vol, div0, ratio, divex
            real :: etime2, etime_p
            LOGICAL          IFPRINT, IFHZPC
            COMMON  /CPRINT/ IFPRINT, IFHZPC
         
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

            ! ** zero out stuff for Lanczos eigenvalue estimator
            call rzero(diagt,maxcg)
            call rzero(upper,maxcg)
            rho = 0.00
            !  
            ! Initialization
            !  
            NXYZ   = lx1*ly1*lz1
            NEL    = NELV
            VOL    = VOLVM1
            IF (IMSH.EQ.2) NEL=NELT
            IF (IMSH.EQ.2) VOL=VOLTM1
            n      = NEL*NXYZ

            tol=abs(tin)

            ! overrule input tolerance
            if (restol(ifield).ne.0) tol=restol(ifield)

            if (tin.lt.0) tol=abs(tin)
            niter = min(maxit,maxcg)

            if (.not.ifsolv) then
               call setfast(h1,h2,imsh)
               ifsolv = .true.
            endif
            !  
            ! Set up diag preconditioner.
            !  
            call setprec(D,h1,h2,imsh,isd)

            call copy (r,f,n)
            call rzero(x,n)
            call rzero(p,n)

            fmax = glamax(f,n)
            if (fmax == 0.0) call nek_stop_error(name//': residual is zero.', this_module, 'solve_helmholtz_2Dh')

            ! Check for non-trivial null-space

            krylov = 0
            rtz1=1.0
            niterhm = 0

            do iter=1,niter
            
               call col3(z,r,d,n) ! Joacobi
            
               rtz2=rtz1
               scalar(1)=vlsc3 (z,r,mult,n)
               if(param(18).eq.1) then
                 scalar(2)=vlsc3(r,r,mult,n)
               else 
                scalar(2)=vlsc32(r,mult,binv,n)
               endif
               call gop(scalar,w,'+  ',2)
               rtz1=scalar(1)
               rbn2=sqrt(scalar(2)/vol)
               if (iter.eq.1) rbn0 = rbn2
               if (param(22).lt.0) tol=abs(param(22))*rbn0
               if (tin.lt.0)       tol=abs(tin)*rbn0

               ifprint_hmh = .false.
               if (nio.eq.0.and.ifprint.and.param(74).ne.0) ifprint_hmh=.true.
               if (nio.eq.0.and.istep.eq.1)                 ifprint_hmh=.true.

               if (ifprint_hmh) write(6,3002) istep,'  Hmholtz ' // name,  iter,rbn2,h1(1),tol,h2(1),ifmcor


               ! Always take at least one iteration   (for projection) pff 11/23/98
               IF (rbn2.LE.TOL.and.(iter.gt.1 .or. istep.le.5)) THEN
               NITER = ITER-1
               if (nio.eq.0) write(6,3000) istep,'  Hmholtz ' // name, niter,rbn2,rbn0,tol
                  goto 9999
               endif
            
               beta = rtz1/rtz2
               if (iter.eq.1) beta=0.0
               call add2s1 (p,z,beta,n)
               call helmholtz_matvec_2Dh (w,p,h1,h2,beta)
               call dssum  (w,lx1,ly1,lz1)
               call col2   (w,mask,n)
            
               rho0 = rho
               rho  = glsc3(w,p,mult,n)
               alpha=rtz1/rho
               alphm=-alpha
               call add2s2(x,p ,alpha,n)
               call add2s2(r,w ,alphm,n)
            
               ! Generate tridiagonal matrix for Lanczos scheme
               if (iter.eq.1) then
                  krylov = krylov+1
                  diagt(iter) = rho/rtz1
               elseif (iter.le.maxcg) then
                  krylov = krylov+1
                  diagt(iter)    = (beta**2 * rho0 + rho ) / rtz1
                  upper(iter-1)  = -beta * rho0 / sqrt(rtz2 * rtz1)
               endif
 1000       enddo
            niter = iter-1
!     
            if (nio.eq.0) write (6,3001) istep, '  Error Hmholtz ' // name, niter,rbn2,rbn0,tol


 3000       format(i11,a,1x,I7,1p4E13.4)
 3001       format(i11,a,1x,I7,1p4E13.4)
 3002       format(i11,a,1x,I7,1p4E13.4,l4)
 9999       continue
            niterhm = niter
            ifsolv = .false.
         end subroutine solve_helmholtz_2Dh

      end module neklab_2Dh
