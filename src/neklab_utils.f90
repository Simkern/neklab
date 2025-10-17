      module neklab_utils
         use stdlib_strings, only: padl
         use stdlib_optval, only: optval
      !---------------------------------------
      !-----     LightKrylov Imports     -----
      !---------------------------------------
      ! Default real kind.
         use LightKrylov, only: dp
         use LightKrylov_Logger
      ! Abstract types for real-valued vectors.
         use LightKrylov, only: abstract_vector, abstract_vector_rdp
      ! Neklab vectors
         use neklab_vectors
      
         implicit none
         include "SIZE"
         include "TOTAL"
         include "ADJOINT"
         private
         character(len=*), parameter, private :: this_module = 'neklab_utils'
      
         integer, parameter :: lv = lx1*ly1*lz1*lelv
      !! Local number of grid points for the velocity mesh.
         integer, parameter :: lp = lx2*ly2*lz2*lelv
      !! Local number of grid points for the pressure mesh.
         integer, parameter :: lt = lx1*ly1*lz1*lelt
      !! Local number of grid points for the temperature/passive scalar mesh.
      
      ! utilities for regular nek vectors
         public :: nek2vec, vec2nek, abs_vec2nek
      ! utilities for extended nek vectors
         public :: get_period, get_period_abs
      ! utility for outposting
         public :: outpost_nek
      ! miscellaneous
         public :: nopcopy
      
      ! Nek vector utilities
         interface nek2vec
            ! nek_dvector
            module procedure nek2vec_std      ! vx, vy, vz, pr, t
            module procedure nek2vec_prt_i    ! vxp, vyp, vzp, prp, tp         pert i
            module procedure nek2vec_prt      ! wrapper for single pert mode
            ! nek_zvector
            module procedure nek2vec_2Dh      ! vxp, vyp, vzp, prp, tp         pert 1:2
            ! nek_ext_dvector
            module procedure nek2ext_vec_std
            module procedure nek2ext_vec_prt_i
            module procedure nek2ext_vec_prt
            ! nekv
            module procedure nek2v_std
            module procedure nek2v_prt_i
            module procedure nek2v_prt
            ! nekp
            module procedure nek2p_std
            module procedure nek2p_prt_i
            module procedure nek2p_prt
         end interface
      
         interface vec2nek
            ! nek_dvector
            module procedure vec2nek_std
            module procedure vec2nek_prt_i
            module procedure vec2nek_prt
            ! nek_zvector
            module procedure vec2nek_2Dh
            ! nek_ext_dvector
            module procedure ext_vec2nek_std
            module procedure ext_vec2nek_prt_i
            module procedure ext_vec2nek_prt
            ! nekv
            module procedure v2nek_std
            module procedure v2nek_prt_i
            module procedure v2nek_prt
            ! nekp
            module procedure p2nek_std
            module procedure p2nek_prt_i
            module procedure p2nek_prt
         end interface
      
         interface abs_vec2nek
            module procedure abstract_vec2nek_std
            module procedure abstract_vec2nek_prt_i
            module procedure abstract_vec2nek_prt
         end interface
      
      ! Outposting
         interface outpost_nek
            module procedure outpost_vector
            module procedure outpost_basis
         end interface
      
      contains

      !--------------------------
      ! Nek5000 ---> LightKrylov
      !--------------------------

      ! nek baseflow -> nek_dvector

         subroutine nek2vec_std(vec, vx_, vy_, vz_, pr_, t_)
            type(nek_dvector), intent(out) :: vec
            real(kind=dp), dimension(lx1, ly1, lz1, lelv), intent(in) :: vx_
            real(kind=dp), dimension(lx1, ly1, lz1, lelv), intent(in) :: vy_
            real(kind=dp), dimension(lx1, ly1, lz1, lelv), intent(in) :: vz_
            real(kind=dp), dimension(lx2, ly2, lz2, lelv), intent(in) :: pr_
            real(kind=dp), dimension(lx1, ly1, lz1, lelt, ldimt), intent(in) :: t_
         
            call nopcopy(vec%vx, vec%vy, vec%vz, vec%pr, vec%theta, vx_, vy_, vz_, pr_, t_)
         
         end subroutine nek2vec_std

      ! nek perturbation -> nek_dvector

         subroutine nek2vec_prt_i(vec, vx_, vy_, vz_, pr_, t_, ipert)
            type(nek_dvector), intent(out) :: vec
            real(kind=dp), dimension(lv, lpert), intent(in) :: vx_
            real(kind=dp), dimension(lv, lpert), intent(in) :: vy_
            real(kind=dp), dimension(lv, lpert), intent(in) :: vz_
            real(kind=dp), dimension(lp, lpert), intent(in) :: pr_
            real(kind=dp), dimension(lx1*ly1*lz1*lelt, ldimt, lpert), intent(in) :: t_
            integer, intent(in) :: ipert
   
            if (ipert > npert) call stop_error('The chosen perturbation index is not defined.',this_module,'nek2vec_prt_i')
            call nopcopy(vec%vx, vec%vy, vec%vz, vec%pr, vec%theta, vx_(:, ipert), vy_(:, ipert), vz_(:, ipert), pr_(:, ipert), t_(:, :, ipert))
            
         end subroutine nek2vec_prt_i
      
         subroutine nek2vec_prt(vec, vx_, vy_, vz_, pr_, t_)
            type(nek_dvector), intent(out) :: vec
            real(kind=dp), dimension(lv, lpert), intent(in) :: vx_
            real(kind=dp), dimension(lv, lpert), intent(in) :: vy_
            real(kind=dp), dimension(lv, lpert), intent(in) :: vz_
            real(kind=dp), dimension(lp, lpert), intent(in) :: pr_
            real(kind=dp), dimension(lx1*ly1*lz1*lelt, ldimt, lpert), intent(in) :: t_
      
            call nek2vec_prt_i(vec, vx_, vy_, vz_, pr_, t_, 1)
      
         end subroutine nek2vec_prt

      ! nek perturbation (1:2) -> nek_zvector

         subroutine nek2vec_2Dh(vec, vx_, vy_, vz_, pr_, t_)
            type(nek_zvector), intent(out) :: vec
            real(kind=dp), dimension(lv, lpert), intent(in) :: vx_
            real(kind=dp), dimension(lv, lpert), intent(in) :: vy_
            real(kind=dp), dimension(lv, lpert), intent(in) :: vz_
            real(kind=dp), dimension(lp, lpert), intent(in) :: pr_
            real(kind=dp), dimension(lx1*ly1*lz1*lelt, ldimt, lpert), intent(in) :: t_
      
            call nek2vec_prt_i(vec%re, vx_, vy_, vz_, pr_, t_, 1)
            call nek2vec_prt_i(vec%im, vx_, vy_, vz_, pr_, t_, 2)
      
         end subroutine nek2vec_2Dh
      
      ! nek baseflow --> nek_ext_dvector

         subroutine nek2ext_vec_std(vec, vx_, vy_, vz_, pr_, t_)
            type(nek_ext_dvector), intent(out) :: vec
            real(kind=dp), dimension(lx1, ly1, lz1, lelv), intent(in) :: vx_
            real(kind=dp), dimension(lx1, ly1, lz1, lelv), intent(in) :: vy_
            real(kind=dp), dimension(lx1, ly1, lz1, lelv), intent(in) :: vz_
            real(kind=dp), dimension(lx2, ly2, lz2, lelv), intent(in) :: pr_
            real(kind=dp), dimension(lx1, ly1, lz1, lelt, ldimt), intent(in) :: t_
      
            call nopcopy(vec%vx, vec%vy, vec%vz, vec%pr, vec%theta, vx_, vy_, vz_, pr_, t_)
      
         end subroutine nek2ext_vec_std

      ! nek perturbation --> nek_ext_dvector

         subroutine nek2ext_vec_prt_i(vec, vx_, vy_, vz_, pr_, t_, ipert)
            type(nek_ext_dvector), intent(out) :: vec
            real(kind=dp), dimension(lv, lpert), intent(in) :: vx_
            real(kind=dp), dimension(lv, lpert), intent(in) :: vy_
            real(kind=dp), dimension(lv, lpert), intent(in) :: vz_
            real(kind=dp), dimension(lp, lpert), intent(in) :: pr_
            real(kind=dp), dimension(lt, ldimt, lpert), intent(in) :: t_
            integer, intent(in) :: ipert
      
            call nopcopy(vec%vx, vec%vy, vec%vz, vec%pr, vec%theta, vx_(:, ipert), vy_(:, ipert), vz_(:, ipert), pr_(:, ipert), t_(:, :, ipert))
      
         end subroutine nek2ext_vec_prt_i

         subroutine nek2ext_vec_prt(vec, vx_, vy_, vz_, pr_, t_)
            type(nek_ext_dvector), intent(out) :: vec
            real(kind=dp), dimension(lv, lpert), intent(in) :: vx_
            real(kind=dp), dimension(lv, lpert), intent(in) :: vy_
            real(kind=dp), dimension(lv, lpert), intent(in) :: vz_
            real(kind=dp), dimension(lp, lpert), intent(in) :: pr_
            real(kind=dp), dimension(lt, ldimt, lpert), intent(in) :: t_
      
            call nek2ext_vec_prt_i(vec, vx_, vy_, vz_, pr_, t_, 1)
      
         end subroutine nek2ext_vec_prt

      ! nek baseflow --> nekv_dvector

         subroutine nek2v_std(vec, v_, mask_, vmult_, isd_)
            type(nekv_dvector), intent(out) :: vec
            real(kind=dp), dimension(lx1, ly1, lz1, lelv), intent(in) :: v_
            real(kind=dp), dimension(lx1, ly1, lz1, lelv), intent(in) :: mask_
            real(kind=dp), dimension(lx1, ly1, lz1, lelv), intent(in) :: vmult_
            integer, intent(in) :: isd_
            ! internal
            integer :: n
            n = lx1*ly1*lz1*nelv
            call copy(vec%v,     v_,       n)
            call copy(vec%mask,  mask_,    n)
            call copy(vec%vmult, vmult_,   n)
            vec%isd = isd_
      
         end subroutine nek2v_std

      ! nek perturbation --> nekv_dvector

         subroutine nek2v_prt_i(vec, v_, mask_, vmult_, isd_, ipert)
            type(nekv_dvector), intent(out) :: vec
            real(kind=dp), dimension(lv, lpert), intent(in) :: v_
            real(kind=dp), dimension(lx1, ly1, lz1, lelv), intent(in) :: mask_
            real(kind=dp), dimension(lx1, ly1, lz1, lelv), intent(in) :: vmult_
            integer, intent(in) :: isd_
            integer, intent(in) :: ipert
            ! internal
            integer :: n
            n = lx1*ly1*lz1*nelv
            call copy(vec%v,     v_(:, ipert), n)
            call copy(vec%mask,  mask_,        n)
            call copy(vec%vmult, vmult_,       n)
            vec%isd = isd_
      
         end subroutine nek2v_prt_i

         subroutine nek2v_prt(vec, v_, mask_, vmult_, isd_)
            type(nekv_dvector), intent(out) :: vec
            real(kind=dp), dimension(lv, lpert), intent(in) :: v_
            real(kind=dp), dimension(lx1, ly1, lz1, lelv), intent(in) :: mask_
            real(kind=dp), dimension(lx1, ly1, lz1, lelv), intent(in) :: vmult_
            integer, intent(in) :: isd_
            call nek2v_prt_i(vec, v_, mask_, vmult_, isd_, 1)
         end subroutine nek2v_prt

      ! nek baseflow --> nekv_dvector

         subroutine nek2p_std(vec, pr_)
            type(nekp_dvector), intent(out) :: vec
            real(kind=dp), dimension(lx2, ly2, lz2, lelv), intent(in) :: pr_
            ! internal
            integer :: n
            n = lx2*ly2*lz2*nelv
            call copy(vec%pr, pr_, n)
      
         end subroutine nek2p_std

      ! nek perturbation --> nekp_dvector

         subroutine nek2p_prt_i(vec, pr_, ipert)
            type(nekp_dvector), intent(out) :: vec
            real(kind=dp), dimension(lp, lpert), intent(in) :: pr_
            integer, intent(in) :: ipert
            ! internal
            integer :: n
            n = lx2*ly2*lz2*nelv
            call copy(vec%pr, pr_(:, ipert), n)
      
         end subroutine nek2p_prt_i

         subroutine nek2p_prt(vec, pr_)
            type(nekp_dvector), intent(out) :: vec
            real(kind=dp), dimension(lp, lpert), intent(in) :: pr_
            call nek2p_prt_i(vec, pr_, 1)
      
         end subroutine nek2p_prt

      !--------------------------
      ! LightKrylov ---> Nek5000
      !--------------------------

      ! nek_dvector -> nek baseflow
      
         subroutine vec2nek_std(vx_, vy_, vz_, pr_, t_, vec)
            type(nek_dvector), intent(in) :: vec
            real(kind=dp), dimension(lx1, ly1, lz1, lelv), intent(out) :: vx_
            real(kind=dp), dimension(lx1, ly1, lz1, lelv), intent(out) :: vy_
            real(kind=dp), dimension(lx1, ly1, lz1, lelv), intent(out) :: vz_
            real(kind=dp), dimension(lx2, ly2, lz2, lelv), intent(out) :: pr_
            real(kind=dp), dimension(lx1, ly1, lz1, lelt, ldimt), intent(out) :: t_
      
            call nopcopy(vx_, vy_, vz_, pr_, t_, vec%vx, vec%vy, vec%vz, vec%pr, vec%theta)
      
         end subroutine vec2nek_std

      ! nek_dvector -> nek perturbation

         subroutine vec2nek_prt_i(vx_, vy_, vz_, pr_, t_, vec, ipert)
            type(nek_dvector), intent(in) :: vec
            integer, intent(in) :: ipert
            real(kind=dp), dimension(lv, lpert), intent(out) :: vx_
            real(kind=dp), dimension(lv, lpert), intent(out) :: vy_
            real(kind=dp), dimension(lv, lpert), intent(out) :: vz_
            real(kind=dp), dimension(lp, lpert), intent(out) :: pr_
            real(kind=dp), dimension(lt, ldimt, lpert), intent(out) :: t_
      
            if (ipert > npert) call stop_error('The chosen perturbation index is not defined.',this_module,'vec2nek_prt_i')
            call nopcopy(vx_(:, ipert), vy_(:, ipert), vz_(:, ipert), pr_(:, ipert), t_(:, :, ipert), vec%vx, vec%vy, vec%vz, vec%pr, vec%theta)
      
         end subroutine vec2nek_prt_i
      
         subroutine vec2nek_prt(vx_, vy_, vz_, pr_, t_, vec)
            type(nek_dvector), intent(in) :: vec
            real(kind=dp), dimension(lv, lpert), intent(out) :: vx_
            real(kind=dp), dimension(lv, lpert), intent(out) :: vy_
            real(kind=dp), dimension(lv, lpert), intent(out) :: vz_
            real(kind=dp), dimension(lp, lpert), intent(out) :: pr_
            real(kind=dp), dimension(lt, ldimt, lpert), intent(out) :: t_
      
            call vec2nek_prt_i(vx_, vy_, vz_, pr_, t_, vec, 1)
      
         end subroutine vec2nek_prt

      ! nek_zvector -> nek perturbation (1:2) 

         subroutine vec2nek_2Dh(vx_, vy_, vz_, pr_, t_, vec)
            type(nek_zvector), intent(in) :: vec
            real(kind=dp), dimension(lv, lpert), intent(out) :: vx_
            real(kind=dp), dimension(lv, lpert), intent(out) :: vy_
            real(kind=dp), dimension(lv, lpert), intent(out) :: vz_
            real(kind=dp), dimension(lp, lpert), intent(out) :: pr_
            real(kind=dp), dimension(lx1*ly1*lz1*lelt, ldimt, lpert), intent(out) :: t_
      
            call vec2nek_prt_i(vx_, vy_, vz_, pr_, t_, vec%re, 1)
            call vec2nek_prt_i(vx_, vy_, vz_, pr_, t_, vec%im, 2)
      
         end subroutine vec2nek_2Dh

      ! nek_ext_dvector --> nek baseflow

         subroutine ext_vec2nek_std(vx_, vy_, vz_, pr_, t_, vec)
            type(nek_ext_dvector), intent(in) :: vec
            real(kind=dp), dimension(lx1, ly1, lz1, lelv), intent(out) :: vx_
            real(kind=dp), dimension(lx1, ly1, lz1, lelv), intent(out) :: vy_
            real(kind=dp), dimension(lx1, ly1, lz1, lelv), intent(out) :: vz_
            real(kind=dp), dimension(lx2, ly2, lz2, lelv), intent(out) :: pr_
            real(kind=dp), dimension(lx1, ly1, lz1, lelt, ldimt), intent(out) :: t_
      
            call nopcopy(vx_, vy_, vz_, pr_, t_, vec%vx, vec%vy, vec%vz, vec%pr, vec%theta)
      
         end subroutine ext_vec2nek_std

      ! nek_ext_dvector --> nek perturbation

         subroutine ext_vec2nek_prt_i(vx_, vy_, vz_, pr_, t_, vec, ipert)
            type(nek_ext_dvector), intent(in) :: vec
            real(kind=dp), dimension(lv, lpert), intent(out) :: vx_
            real(kind=dp), dimension(lv, lpert), intent(out) :: vy_
            real(kind=dp), dimension(lv, lpert), intent(out) :: vz_
            real(kind=dp), dimension(lp, lpert), intent(out) :: pr_
            real(kind=dp), dimension(lt, ldimt, lpert), intent(out) :: t_
            integer, intent(in) :: ipert
      
            call nopcopy(vx_(:, ipert), vy_(:, ipert), vz_(:, ipert), pr_(:, ipert), t_(:, :, ipert), vec%vx, vec%vy, vec%vz, vec%pr, vec%theta)
      
         end subroutine ext_vec2nek_prt_i

         subroutine ext_vec2nek_prt(vx_, vy_, vz_, pr_, t_, vec)
            type(nek_ext_dvector), intent(in) :: vec
            real(kind=dp), dimension(lv, lpert), intent(out) :: vx_
            real(kind=dp), dimension(lv, lpert), intent(out) :: vy_
            real(kind=dp), dimension(lv, lpert), intent(out) :: vz_
            real(kind=dp), dimension(lp, lpert), intent(out) :: pr_
            real(kind=dp), dimension(lt, ldimt, lpert), intent(out) :: t_
      
            call ext_vec2nek_prt_i(vx_, vy_, vz_, pr_, t_, vec, 1)
      
         end subroutine ext_vec2nek_prt

      ! nekv_dvector --> nek baseflow

         subroutine v2nek_std(v_, vec)
            type(nekv_dvector), intent(in) :: vec
            real(kind=dp), dimension(lx1, ly1, lz1, lelv), intent(out) :: v_
            ! internal
            integer :: n
            n = lx1*ly1*lz1*nelv
            call copy(v_, vec%v, n)
      
         end subroutine v2nek_std

      ! nekv_dvector --> nek perturbation
      
         subroutine v2nek_prt_i(v_, vec, ipert)
            type(nekv_dvector), intent(in) :: vec
            real(kind=dp), dimension(lv, lpert), intent(out) :: v_
            integer, intent(in) :: ipert
            ! internal
            integer :: n
            n = lx1*ly1*lz1*nelv
            call copy(v_(:, ipert), vec%v, n)
      
         end subroutine v2nek_prt_i

         subroutine v2nek_prt(v_, vec)
            type(nekv_dvector), intent(in) :: vec
            real(kind=dp), dimension(lv, lpert), intent(out) :: v_
            call v2nek_prt_i(v_, vec, 1)
      
         end subroutine v2nek_prt

      ! nekp_dvector --> nek baseflow

         subroutine p2nek_std(pr_, vec)
            type(nekp_dvector), intent(in) :: vec
            real(kind=dp), dimension(lx2, ly2, lz2, lelv), intent(out) :: pr_
            ! internal
            integer :: n
            n = lx2*ly2*lz2*nelv
            call copy(pr_, vec%pr, n)
      
         end subroutine p2nek_std

      ! nekp_dvector --> nek perturbation
      
         subroutine p2nek_prt_i(pr_, vec, ipert)
            type(nekp_dvector), intent(in) :: vec
            real(kind=dp), dimension(lp, lpert), intent(out) :: pr_
            integer, intent(in) :: ipert
            ! internal
            integer :: n
            n = lx2*ly2*lz2*nelv
            call copy(pr_(:,ipert), vec%pr, n)
      
         end subroutine p2nek_prt_i

         subroutine p2nek_prt(pr_, vec)
            type(nekp_dvector), intent(in) :: vec
            real(kind=dp), dimension(lp, lpert), intent(out) :: pr_
            call p2nek_prt_i(pr_, vec, 1)
      
         end subroutine p2nek_prt

      ! abstract_vector_rdp -> nek baseflow
      
         subroutine abstract_vec2nek_std(vx_, vy_, vz_, pr_, t_, vec)
            class(abstract_vector_rdp), intent(in) :: vec
            real(kind=dp), dimension(lx1, ly1, lz1, lelv), intent(out) :: vx_
            real(kind=dp), dimension(lx1, ly1, lz1, lelv), intent(out) :: vy_
            real(kind=dp), dimension(lx1, ly1, lz1, lelv), intent(out) :: vz_
            real(kind=dp), dimension(lx2, ly2, lz2, lelv), intent(out) :: pr_
            real(kind=dp), dimension(lx1, ly1, lz1, lelt, ldimt), intent(out) :: t_
            select type (vec)
            type is (nek_dvector)
               call vec2nek(vx_, vy_, vz_, pr_, t_, vec)
            type is (nek_ext_dvector)
               call vec2nek(vx_, vy_, vz_, pr_, t_, vec)
            class default
               call type_error('vec','nek_dvector/nek_ext_dvector','IN',this_module,'abstract_vec2nek_std')
            end select
         end subroutine abstract_vec2nek_std

      ! abstract_vector_rdp -> nek perturbation
      
         subroutine abstract_vec2nek_prt_i(vx_, vy_, vz_, pr_, t_, vec, ipert)
            class(abstract_vector_rdp), intent(in) :: vec
            real(kind=dp), dimension(lv, lpert), intent(out) :: vx_
            real(kind=dp), dimension(lv, lpert), intent(out) :: vy_
            real(kind=dp), dimension(lv, lpert), intent(out) :: vz_
            real(kind=dp), dimension(lp, lpert), intent(out) :: pr_
            real(kind=dp), dimension(lt, ldimt, lpert), intent(out) :: t_
            integer, intent(in) :: ipert
            select type (vec)
            type is (nek_dvector)
               call vec2nek(vx_, vy_, vz_, pr_, t_, vec, ipert)
            type is (nek_ext_dvector)
               call vec2nek(vx_, vy_, vz_, pr_, t_, vec, ipert)
            class default
               call type_error('vec','nek_dvector/nek_ext_dvector','IN',this_module,'abstract_vec2nek_prt')
            end select
         end subroutine abstract_vec2nek_prt_i

         subroutine abstract_vec2nek_prt(vx_, vy_, vz_, pr_, t_, vec)
            class(abstract_vector_rdp), intent(in) :: vec
            real(kind=dp), dimension(lv, lpert), intent(out) :: vx_
            real(kind=dp), dimension(lv, lpert), intent(out) :: vy_
            real(kind=dp), dimension(lv, lpert), intent(out) :: vz_
            real(kind=dp), dimension(lp, lpert), intent(out) :: pr_
            real(kind=dp), dimension(lt, ldimt, lpert), intent(out) :: t_
            call abstract_vec2nek_prt_i(vx_, vy_, vz_, pr_, t_, vec, 1)
         end subroutine abstract_vec2nek_prt
      
      ! EXTENDED vector utils

         real(dp) function get_period_abs(vec) result(period)
            class(abstract_vector_rdp), intent(in) :: vec
            select type (vec)
            type is (nek_ext_dvector)
               
               period = vec%T

            class default
               call type_error('vec','nek_ext_dvector','IN',this_module,'get_period_abs')
            end select
         end function get_period_abs
      
         pure real(dp) function get_period(vec) result(period)
            class(nek_ext_dvector), intent(in) :: vec
            period = vec%T
         end function get_period

         subroutine nopcopy(a1, a2, a3, a4, a5, b1, b2, b3, b4, b5)
            implicit none
            include 'SIZE'
            include 'TOTAL'
            integer :: n, k
            real(kind=dp), intent(inout) :: a1(1), a2(1), a3(1), a4(1), a5(lx1*ly1*lz1*lelt, 1)
            real(kind=dp), intent(in) :: b1(1), b2(1), b3(1), b4(1), b5(lx1*ly1*lz1*lelt, 1)
            n = nx1*ny1*nz1*nelv
            call copy(a1, b1, n)
            call copy(a2, b2, n)
            if (if3D) call copy(a3, b3, n)
            if (ifpo) call copy(a4, b4, nx2*ny2*nz2*nelv)
            if (ifto) call copy(a5(1, 1), b5(1, 1), lx1*ly1*lz1*nelfld(2))
            if (ldimt > 1) then
            do k = 1, npscal
               if (ifpsco(k)) call copy(a5(1, k + 1), b5(1, k + 1), lx1*ly1*lz1*nelfld(k + 2))
            end do
            end if
         end subroutine nopcopy
      
      ! OUTPOSTING

      ! abstract_vector

         subroutine outpost_vector(vec, prefix)
            class(abstract_vector), intent(in) :: vec
            character(len=3), intent(in) :: prefix
            ! internal
            real(dp) :: t_tmp
            character(len=3) :: prefix2
            select type (vec)
            type is (nek_dvector)
               call outpost(vec%vx, vec%vy, vec%vz, vec%pr, vec%theta, prefix)
            type is (nek_zvector)
               write(prefix2,'(A2,A1)') prefix(1:2), 'r' 
               associate (v => vec%re)
                  call outpost(v%vx, v%vy, v%vz, v%pr, v%theta, prefix2)
               end associate
               write(prefix2,'(A2,A1)') prefix(1:2), 'i' 
               associate (v => vec%im)
                  call outpost(v%vx, v%vy, v%vz, v%pr, v%theta, prefix2)
               end associate
            type is (nek_ext_dvector)
               t_tmp = time 
               time = vec%T
               call outpost(vec%vx, vec%vy, vec%vz, vec%pr, vec%theta, prefix)
               time = t_tmp
            class default
               call type_error('vec','nek_dvector/zvector/ext_dvector','IN',this_module,'outpost_vector')
            end select
         end subroutine outpost_vector

         subroutine outpost_basis(vec, prefix)
            class(abstract_vector), intent(in) :: vec(:)
            character(len=3), intent(in) :: prefix
            ! internal
            integer :: i
            do i = 1, size(vec)
               call outpost_vector(vec(i), prefix)
            end do
         end subroutine outpost_basis      

      end module neklab_utils
