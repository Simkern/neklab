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
         use LightKrylov, only: abstract_vector_rdp
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
         public :: nek2vec, vec2nek, abs_vec2nek, outpost_dnek, outpost_dnek_abs_vector
      ! utilities for extended nek vectors
         public :: nek2ext_vec, ext_vec2nek, abs_ext_vec2nek, outpost_ext_dnek
         public :: get_period, get_period_abs
      ! Utilities for velocity components
         public :: nek2vcomp, vcomp2nek
      ! miscellaneous
         public :: nopcopy
      
      ! Nek vector utilities
         interface nek2vec
            module procedure nek2vec_std
            module procedure nek2vec_prt
         end interface
      
         interface vec2nek
            module procedure vec2nek_std
            module procedure vec2nek_prt
         end interface
      
         interface abs_vec2nek
            module procedure abstract_vec2nek_std
            module procedure abstract_vec2nek_prt
         end interface
      
         interface outpost_dnek
            module procedure outpost_dnek_vector
            module procedure outpost_dnek_basis
         end interface
      
      ! Extended nek vector utilities
         interface nek2ext_vec
            module procedure nek2ext_vec_std
            module procedure nek2ext_vec_prt
         end interface
      
         interface ext_vec2nek
            module procedure ext_vec2nek_std
            module procedure ext_vec2nek_prt
         end interface
      
         interface abs_ext_vec2nek
            module procedure abstract_ext_vec2nek_std
            module procedure abstract_ext_vec2nek_prt
         end interface
      
         interface outpost_ext_dnek
            module procedure outpost_ext_dnek_vector
            module procedure outpost_ext_dnek_basis
         end interface

      ! V component utilities
         interface nek2vcomp
            module procedure nek2vcomp_std
            module procedure nek2vcomp_prt
         end interface
      
         interface vcomp2nek
            module procedure vcomp2nek_std
            module procedure vcomp2nek_prt
         end interface
      
      contains
      
         subroutine nek2vec_prt(vec, vx_, vy_, vz_, pr_, t_)
            include "SIZE"
            type(nek_dvector), intent(out) :: vec
            real(kind=dp), dimension(lv, lpert), intent(in) :: vx_
            real(kind=dp), dimension(lv, lpert), intent(in) :: vy_
            real(kind=dp), dimension(lv, lpert), intent(in) :: vz_
            real(kind=dp), dimension(lp, lpert), intent(in) :: pr_
            real(kind=dp), dimension(lx1*ly1*lz1*lelt, ldimt, lpert), intent(in) :: t_
      
            call nopcopy(vec%vx, vec%vy, vec%vz, vec%pr, vec%theta, vx_(:, 1), vy_(:, 1), vz_(:, 1), pr_(:, 1), t_(:, :, 1))
      
         end subroutine nek2vec_prt
      
         subroutine nek2vec_std(vec, vx_, vy_, vz_, pr_, t_)
            include "SIZE"
            type(nek_dvector), intent(out) :: vec
            real(kind=dp), dimension(lx1, ly1, lz1, lelv), intent(in) :: vx_
            real(kind=dp), dimension(lx1, ly1, lz1, lelv), intent(in) :: vy_
            real(kind=dp), dimension(lx1, ly1, lz1, lelv), intent(in) :: vz_
            real(kind=dp), dimension(lx2, ly2, lz2, lelv), intent(in) :: pr_
            real(kind=dp), dimension(lx1, ly1, lz1, lelt, ldimt), intent(in) :: t_
      
            call nopcopy(vec%vx, vec%vy, vec%vz, vec%pr, vec%theta, vx_, vy_, vz_, pr_, t_)
      
         end subroutine nek2vec_std
      
         subroutine vec2nek_std(vx_, vy_, vz_, pr_, t_, vec)
            include "SIZE"
            type(nek_dvector), intent(in) :: vec
            real(kind=dp), dimension(lx1, ly1, lz1, lelv), intent(out) :: vx_
            real(kind=dp), dimension(lx1, ly1, lz1, lelv), intent(out) :: vy_
            real(kind=dp), dimension(lx1, ly1, lz1, lelv), intent(out) :: vz_
            real(kind=dp), dimension(lx2, ly2, lz2, lelv), intent(out) :: pr_
            real(kind=dp), dimension(lx1, ly1, lz1, lelt, ldimt), intent(out) :: t_
      
            call nopcopy(vx_, vy_, vz_, pr_, t_, vec%vx, vec%vy, vec%vz, vec%pr, vec%theta)
      
         end subroutine vec2nek_std
      
         subroutine vec2nek_prt(vx_, vy_, vz_, pr_, t_, vec)
            include "SIZE"
            type(nek_dvector), intent(in) :: vec
            real(kind=dp), dimension(lv, 1), intent(out) :: vx_
            real(kind=dp), dimension(lv, 1), intent(out) :: vy_
            real(kind=dp), dimension(lv, 1), intent(out) :: vz_
            real(kind=dp), dimension(lp, 1), intent(out) :: pr_
            real(kind=dp), dimension(lt, ldimt, 1), intent(out) :: t_
      
            call nopcopy(vx_(:, 1), vy_(:, 1), vz_(:, 1), pr_(:, 1), t_(:, :, 1), vec%vx, vec%vy, vec%vz, vec%pr, vec%theta)
      
         end subroutine vec2nek_prt
      
         subroutine abstract_vec2nek_std(vx_, vy_, vz_, pr_, t_, vec)
            include "SIZE"
            class(abstract_vector_rdp), intent(in) :: vec
            real(kind=dp), dimension(lx1, ly1, lz1, lelv), intent(out) :: vx_
            real(kind=dp), dimension(lx1, ly1, lz1, lelv), intent(out) :: vy_
            real(kind=dp), dimension(lx1, ly1, lz1, lelv), intent(out) :: vz_
            real(kind=dp), dimension(lx2, ly2, lz2, lelv), intent(out) :: pr_
            real(kind=dp), dimension(lx1, ly1, lz1, lelt, ldimt), intent(out) :: t_
            select type (vec)
            type is (nek_dvector)
               
               call nopcopy(vx_, vy_, vz_, pr_, t_, vec%vx, vec%vy, vec%vz, vec%pr, vec%theta)

            class default
               call type_error('vec','nek_dvector','IN',this_module,'abstract_vec2nek_std')
            end select
         end subroutine abstract_vec2nek_std
      
         subroutine abstract_vec2nek_prt(vx_, vy_, vz_, pr_, t_, vec)
            include "SIZE"
            class(abstract_vector_rdp), intent(in) :: vec
            real(kind=dp), dimension(lv, 1), intent(out) :: vx_
            real(kind=dp), dimension(lv, 1), intent(out) :: vy_
            real(kind=dp), dimension(lv, 1), intent(out) :: vz_
            real(kind=dp), dimension(lp, 1), intent(out) :: pr_
            real(kind=dp), dimension(lt, ldimt, 1), intent(out) :: t_
            select type (vec)
            type is (nek_dvector)
               
               call nopcopy(vx_(:, 1), vy_(:, 1), vz_(:, 1), pr_(:, 1), t_(:, :, 1), vec%vx, vec%vy, vec%vz, vec%pr, vec%theta)

            class default
               call type_error('vec','nek_dvector','IN',this_module,'abstract_vec2nek_prt')
            end select
         end subroutine abstract_vec2nek_prt
      
      ! EXTENDED
      
         subroutine nek2ext_vec_prt(vec, vx_, vy_, vz_, pr_, t_)
            include "SIZE"
            type(nek_ext_dvector), intent(out) :: vec
            real(kind=dp), dimension(lv, 1), intent(in) :: vx_
            real(kind=dp), dimension(lv, 1), intent(in) :: vy_
            real(kind=dp), dimension(lv, 1), intent(in) :: vz_
            real(kind=dp), dimension(lp, 1), intent(in) :: pr_
            real(kind=dp), dimension(lt, ldimt, 1), intent(in) :: t_
      
            call nopcopy(vec%vx, vec%vy, vec%vz, vec%pr, vec%theta, vx_(:, 1), vy_(:, 1), vz_(:, 1), pr_(:, 1), t_(:, :, 1))
      
         end subroutine nek2ext_vec_prt
      
         subroutine nek2ext_vec_std(vec, vx_, vy_, vz_, pr_, t_)
            include "SIZE"
            type(nek_ext_dvector), intent(out) :: vec
            real(kind=dp), dimension(lx1, ly1, lz1, lelv), intent(in) :: vx_
            real(kind=dp), dimension(lx1, ly1, lz1, lelv), intent(in) :: vy_
            real(kind=dp), dimension(lx1, ly1, lz1, lelv), intent(in) :: vz_
            real(kind=dp), dimension(lx2, ly2, lz2, lelv), intent(in) :: pr_
            real(kind=dp), dimension(lx1, ly1, lz1, lelt, ldimt), intent(in) :: t_
      
            call nopcopy(vec%vx, vec%vy, vec%vz, vec%pr, vec%theta, vx_, vy_, vz_, pr_, t_)
      
         end subroutine nek2ext_vec_std
      
         subroutine ext_vec2nek_std(vx_, vy_, vz_, pr_, t_, vec)
            include "SIZE"
            type(nek_ext_dvector), intent(in) :: vec
            real(kind=dp), dimension(lx1, ly1, lz1, lelv), intent(out) :: vx_
            real(kind=dp), dimension(lx1, ly1, lz1, lelv), intent(out) :: vy_
            real(kind=dp), dimension(lx1, ly1, lz1, lelv), intent(out) :: vz_
            real(kind=dp), dimension(lx2, ly2, lz2, lelv), intent(out) :: pr_
            real(kind=dp), dimension(lx1, ly1, lz1, lelt, ldimt), intent(out) :: t_
      
            call nopcopy(vx_, vy_, vz_, pr_, t_, vec%vx, vec%vy, vec%vz, vec%pr, vec%theta)
      
         end subroutine ext_vec2nek_std
      
         subroutine ext_vec2nek_prt(vx_, vy_, vz_, pr_, t_, vec)
            include "SIZE"
            type(nek_ext_dvector), intent(in) :: vec
            real(kind=dp), dimension(lv, 1), intent(out) :: vx_
            real(kind=dp), dimension(lv, 1), intent(out) :: vy_
            real(kind=dp), dimension(lv, 1), intent(out) :: vz_
            real(kind=dp), dimension(lp, 1), intent(out) :: pr_
            real(kind=dp), dimension(lt, ldimt, 1), intent(out) :: t_
      
            call nopcopy(vx_(:, 1), vy_(:, 1), vz_(:, 1), pr_(:, 1), t_(:, :, 1), vec%vx, vec%vy, vec%vz, vec%pr, vec%theta)
      
         end subroutine ext_vec2nek_prt
      
         subroutine abstract_ext_vec2nek_std(vx_, vy_, vz_, pr_, t_, vec)
            include "SIZE"
            class(abstract_vector_rdp), intent(in) :: vec
            real(kind=dp), dimension(lx1, ly1, lz1, lelv), intent(out) :: vx_
            real(kind=dp), dimension(lx1, ly1, lz1, lelv), intent(out) :: vy_
            real(kind=dp), dimension(lx1, ly1, lz1, lelv), intent(out) :: vz_
            real(kind=dp), dimension(lx2, ly2, lz2, lelv), intent(out) :: pr_
            real(kind=dp), dimension(lx1, ly1, lz1, lelt, ldimt), intent(out) :: t_
            select type (vec)
            type is (nek_ext_dvector)
               
               call nopcopy(vx_, vy_, vz_, pr_, t_, vec%vx, vec%vy, vec%vz, vec%pr, vec%theta)

            class default
               call type_error('vec','nek_ext_dvector','IN',this_module,'abstract_ext_vec2nek_std')
            end select
         end subroutine abstract_ext_vec2nek_std
      
         subroutine abstract_ext_vec2nek_prt(vx_, vy_, vz_, pr_, t_, vec)
            include "SIZE"
            class(abstract_vector_rdp), intent(in) :: vec
            real(kind=dp), dimension(lv, 1), intent(out) :: vx_
            real(kind=dp), dimension(lv, 1), intent(out) :: vy_
            real(kind=dp), dimension(lv, 1), intent(out) :: vz_
            real(kind=dp), dimension(lp, 1), intent(out) :: pr_
            real(kind=dp), dimension(lt, ldimt, 1), intent(out) :: t_
            select type (vec)
            type is (nek_ext_dvector)
               
               call nopcopy(vx_(:, 1), vy_(:, 1), vz_(:, 1), pr_(:, 1), t_(:, :, 1), vec%vx, vec%vy, vec%vz, vec%pr, vec%theta)

            class default
               call type_error('vec','nek_ext_dvector','IN',this_module,'abstract_ext_vec2nek_prt')
            end select
         end subroutine abstract_ext_vec2nek_prt

      ! V component

         subroutine nek2vcomp_prt(vec, v_, mask_, vmult_, isd_)
            include "SIZE"
            type(nekv_dvector), intent(out) :: vec
            real(kind=dp), dimension(lv, lpert), intent(in) :: v_
            real(kind=dp), dimension(lx1, ly1, lz1, lelv), intent(in) :: mask_
            real(kind=dp), dimension(lx1, ly1, lz1, lelv), intent(in) :: vmult_
            integer, intent(in) :: isd_
            ! internal
            integer :: n
            n = lx1*ly1*lz1*nelv
            call copy(vec%v,     v_(:, 1), n)
            call copy(vec%mask,  mask_,    n)
            call copy(vec%vmult, vmult_,   n)
            vec%isd = isd_
      
         end subroutine nek2vcomp_prt
      
         subroutine nek2vcomp_std(vec, v_, mask_, vmult_, isd_)
            include "SIZE"
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
      
         end subroutine nek2vcomp_std
      
         subroutine vcomp2nek_std(v_, vec)
            include "SIZE"
            type(nekv_dvector), intent(in) :: vec
            real(kind=dp), dimension(lx1, ly1, lz1, lelv), intent(out) :: v_
            ! internal
            integer :: n
            n = lx1*ly1*lz1*nelv
            call copy(v_,        vec%v, n)
      
         end subroutine vcomp2nek_std
      
         subroutine vcomp2nek_prt(v_, vec)
            include "SIZE"
            type(nekv_dvector), intent(in) :: vec
            real(kind=dp), dimension(lv, 1), intent(out) :: v_
            ! internal
            integer :: n
            n = lx1*ly1*lz1*nelv
            call copy(v_,        vec%v,    n)
      
         end subroutine vcomp2nek_prt
      
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
         
         subroutine outpost_dnek_vector(vec, prefix)
            type(nek_dvector), intent(in) :: vec
            character(len=3), intent(in) :: prefix
            call outpost(vec%vx, vec%vy, vec%vz, vec%pr, vec%theta, prefix)
         end subroutine outpost_dnek_vector
      
         subroutine outpost_dnek_abs_vector(vec, prefix)
            class(abstract_vector_rdp), intent(in) :: vec
            character(len=3), intent(in) :: prefix
            select type (vec)
            type is (nek_dvector)
               
               call outpost(vec%vx, vec%vy, vec%vz, vec%pr, vec%theta, prefix)

            class default
               call type_error('vec','nek_dvector','IN',this_module,'outpost_dnek_abs_vector')
            end select
         end subroutine outpost_dnek_abs_vector
      
         subroutine outpost_dnek_basis(vec, prefix)
            type(nek_dvector), intent(in) :: vec(:)
            character(len=3), intent(in) :: prefix
            integer :: i
            do i = 1, size(vec)
               call outpost_dnek_vector(vec(i), prefix)
            end do
         end subroutine outpost_dnek_basis
      
         subroutine outpost_ext_dnek_vector(vec, prefix)
            type(nek_ext_dvector), intent(in) :: vec
            character(len=3), intent(in) :: prefix
            call outpost(vec%vx, vec%vy, vec%vz, vec%pr, vec%theta, prefix)
         end subroutine outpost_ext_dnek_vector
      
         subroutine outpost_ext_dnek_basis(vec, prefix)
            type(nek_ext_dvector), intent(in) :: vec(:)
            character(len=3), intent(in) :: prefix
            integer :: i
            do i = 1, size(vec)
               call outpost_ext_dnek_vector(vec(i), prefix)
            end do
         end subroutine outpost_ext_dnek_basis

      end module neklab_utils
