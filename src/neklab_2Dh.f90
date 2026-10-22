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
     &                                ebar_2Dh, incomprp_solve_2Dh, incomprp_update_2Dh
         use neklab_2Dh_planar, only: prepare_step_2Dh_planar, make_coupling_2Dh_planar,
     &                                makefp_2Dh_planar, ophx_2Dh_planar, ophinv_2Dh_planar,
     &                                incomprp_div_2Dh_planar, incomprp_dw_2Dh_planar,
     &                                pres_matvec_planar, maxit_p_planar
         use neklab_2Dh_axisym, only: prepare_step_2Dh_axisym, make_coupling_2Dh_axisym,
     &                                makefp_2Dh_axisym, ophx_2Dh_axisym,
     &                                ophinv_uz_2Dh_axisym, ophinv_urphi_2Dh_axisym,
     &                                incomprp_div_2Dh_axisym, incomprp_dw_2Dh_axisym,
     &                                incomprp_coupled_2Dh_axisym,
     &                                pres_matvec_axisym, maxit_p_axisym, alpha_eq_zero
         use neklab_2Dh_torsion, only: if_torsion_zero
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

         private
         character(len=*), parameter, private :: this_module = 'neklab_2Dh'

      ! Fields handed from one stage to the next.
      ! m1
         real(dp), dimension(lx1*ly1*lz1*lelv,lpert),   private :: gradw    ! third-component pressure gradient (make_coupling -> cresvipp)
         real(dp), dimension(lx1*ly1*lz1*lelv,3,lpert), private :: resv_jp  ! momentum residuals per slot (axisym: ophx after the slot loop)
         real(dp), dimension(lx1*ly1*lz1*lelv),         private :: dw       ! third-component pressure correction
      ! m2
         real(dp), dimension(lx2*ly2*lz2*lelv,lpert),   private :: dpr      ! pressure correction

         public :: nek_advance_2Dh_planar
         public :: nek_advance_2Dh_axisym

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
            real(dp), dimension(lx1,ly1,lz1,lelv) :: w1, w2, w3, resv1, resv2, resv3, dv1, dv2, dv3
            real(dp), dimension(lx1,ly1,lz1,lelv) :: h1, h2, h2inv
            common /SCRNS/  w1, w2, w3, resv1, resv2, resv3, dv1, dv2, dv3
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

      end module neklab_2Dh
