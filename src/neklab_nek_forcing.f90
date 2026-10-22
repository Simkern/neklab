      module neklab_nek_forcing
         use LightKrylov, only: dp
         use LightKrylov_Logger
         use neklab_nek_setup, only: nek_log_message, nek_log_debug, nek_stop_error

         implicit none
         include "SIZE"
         include "TOTAL"
         include "ADJOINT"
         private
         character(len=*), parameter, private :: this_module = 'neklab_forcing'
      
         integer, parameter :: lv = lx1*ly1*lz1*lelv
         integer, parameter :: lp = lx2*ly2*lz2*lelv
      
         public :: get_neklab_forcing, set_neklab_forcing
         public :: zero_neklab_forcing, zero_neklab_forcing_ipert
         public :: neklab_forcing
         public :: get_neklab_source, set_neklab_source
         public :: neklab_source
         public :: add_neklab_forcing_term, add_neklab_source_term, remove_neklab_term
         public :: pointwise_forcing_term, pointwise_source_term
      
         real(dp), private, dimension(lv, lpert + 1) :: neklab_ffx = 0.0_dp
         real(dp), private, dimension(lv, lpert + 1) :: neklab_ffy = 0.0_dp
         real(dp), private, dimension(lv, lpert + 1) :: neklab_ffz = 0.0_dp
         real(dp), private, dimension(lv, lpert + 1) :: neklab_fq  = 0.0_dp

      ! Pointwise terms for the nonlinear solver (ipert = 0 only), added on read.
      ! They return ACCELERATIONS and ADD to their arguments. e is the LOCAL element.
         abstract interface
            subroutine pointwise_forcing_term(fx, fy, fz, ix, iy, iz, e)
               import :: dp
               real(dp), intent(inout) :: fx, fy, fz
               integer, intent(in) :: ix, iy, iz, e
            end subroutine pointwise_forcing_term
            subroutine pointwise_source_term(q, ix, iy, iz, e)
               import :: dp
               real(dp), intent(inout) :: q
               integer, intent(in) :: ix, iy, iz, e
            end subroutine pointwise_source_term
         end interface

         integer, parameter, private :: max_terms = 8
         type, private :: forcing_term_t
            character(len=32) :: name = ''
            procedure(pointwise_forcing_term), pointer, nopass :: fn => null()
         end type forcing_term_t
         type, private :: source_term_t
            character(len=32) :: name = ''
            procedure(pointwise_source_term), pointer, nopass :: fn => null()
         end type source_term_t
         type(forcing_term_t), dimension(max_terms), private :: f_terms
         type(source_term_t),  dimension(max_terms), private :: q_terms
         integer, private :: nf_terms = 0
         integer, private :: nq_terms = 0
      
      contains
      
         subroutine zero_neklab_forcing()
      ! internal
            integer :: ipert
            character(len=128) :: msg
            write (msg, *) 'Setting all neklab_forcing vectors to zero.'
            call nek_log_debug(msg, this_module, procedure='zero_neklab_forcing')
            neklab_ffx = 0.0_dp
            neklab_ffy = 0.0_dp
            neklab_ffz = 0.0_dp
            neklab_fq  = 0.0_dp
         end subroutine zero_neklab_forcing
      
         subroutine get_neklab_forcing(fx, fy, fz, ipert)
            real(dp), dimension(lv), intent(out) :: fx
            real(dp), dimension(lv), intent(out) :: fy
            real(dp), dimension(lv), intent(out) :: fz
            integer, intent(in) :: ipert
      ! internal
            character(len=128) :: msg
            if (ipert < 0 .or. ipert > lpert) then
               write (msg, '(A,I0)') 'Invalid value for ipert specified! ipert = ', ipert
               call nek_log_message(msg, this_module, procedure='get_neklab_forcing')
               if (nid == 0) print *, trim(msg)
               call nek_end()
            else
               write (msg, '(A,I0)') 'Retrieving value of the neklab forcing. ipert = ', ipert
               call nek_log_debug(msg, this_module, procedure='get_neklab_forcing')
            end if
            fx = neklab_ffx(:, ipert + 1)
            fy = neklab_ffy(:, ipert + 1)
            fz = neklab_ffz(:, ipert + 1)
         end subroutine get_neklab_forcing
      
         subroutine set_neklab_forcing(fx, fy, fz, ipert)
            real(dp), dimension(lv), intent(in) :: fx
            real(dp), dimension(lv), intent(in) :: fy
            real(dp), dimension(lv), intent(in) :: fz
            integer, intent(in) :: ipert
      ! internal
            character(len=128) :: msg
            if (ipert < 0 .or. ipert > lpert) then
               write (msg, '(A,I0)') 'Invalid value for ipert specified! ipert = ', ipert
               call nek_log_message(msg, this_module, procedure='set_neklab_forcing')
               if (nid == 0) print *, trim(msg)
               call nek_end()
            else
               write (msg, '(A,I0)') 'Setting value of the neklab forcing. ipert = ', ipert
               call nek_log_debug(msg, this_module, procedure='get_neklab_forcing')
            end if
            neklab_ffx(:, ipert + 1) = fx
            neklab_ffy(:, ipert + 1) = fy
            neklab_ffz(:, ipert + 1) = fz
         end subroutine set_neklab_forcing
      
         subroutine zero_neklab_forcing_ipert(ipert)
            integer, intent(in) :: ipert
      ! internal
            character(len=128) :: msg
            if (ipert < 0 .or. ipert > lpert) then
               write (msg, '(A,I0)') 'Invalid value for ipert specified! ipert = ', ipert
               call nek_log_message(msg, this_module, procedure='zero_neklab_forcing_ipert')
               call nek_end()
            else
               write (msg, '(A,I0)') 'Setting value of the neklab forcing to zero. ipert = ', ipert
               call nek_log_debug(msg, this_module, procedure='zero_neklab_forcing_ipert')
            end if
            neklab_ffx(:, ipert + 1) = 0.0_dp
            neklab_ffy(:, ipert + 1) = 0.0_dp
            neklab_ffz(:, ipert + 1) = 0.0_dp
            neklab_fq (:, ipert + 1) = 0.0_dp
         end subroutine zero_neklab_forcing_ipert
      
         subroutine neklab_forcing(ffx, ffy, ffz, ix, iy, iz, ieg, ipert)
            real(dp), intent(inout) :: ffx
            real(dp), intent(inout) :: ffy
            real(dp), intent(inout) :: ffz
            integer, intent(in) :: ix
            integer, intent(in) :: iy
            integer, intent(in) :: iz
            integer, intent(in) :: ieg
            integer, intent(in) :: ipert
      ! internal
            integer :: e, ijke, i
            e = gllel(ieg)
            ijke = ix + lx1*((iy - 1) + ly1*((iz - 1) + lz1*(e - 1)))
      ! note: ipert/jp == 0 : forcing for nonlinear solver
      ! note: ipert/jp >  0 : forcing for linear solver
            ffx = ffx + neklab_ffx(ijke, ipert + 1)
            ffy = ffy + neklab_ffy(ijke, ipert + 1)
            ffz = ffz + neklab_ffz(ijke, ipert + 1)
      ! registered pointwise terms: nonlinear solver only
            if (ipert == 0) then
               do i = 1, nf_terms
                  call f_terms(i)%fn(ffx, ffy, ffz, ix, iy, iz, e)
               end do
            end if
         end subroutine neklab_forcing

         subroutine get_neklab_source(fq, ipert)
            real(dp), dimension(lv), intent(out) :: fq
            integer, intent(in) :: ipert
            call check_ipert(ipert, 'get_neklab_source')
            fq = neklab_fq(:, ipert + 1)
         end subroutine get_neklab_source

         subroutine set_neklab_source(fq, ipert)
            real(dp), dimension(lv), intent(in) :: fq
            integer, intent(in) :: ipert
            call check_ipert(ipert, 'set_neklab_source')
            neklab_fq(:, ipert + 1) = fq
         end subroutine set_neklab_source

         subroutine neklab_source(qvol, ix, iy, iz, ieg, ipert)
      !! Third velocity component (temperature slot). Call from userq. Nek weights
      !! qvol by bm1 only (makeuq), while ffx/ffy get rho*bm1 (makeuf); the sum is
      !! therefore multiplied by rhoCp here, so that every term is an acceleration.
            real(dp), intent(inout) :: qvol
            integer, intent(in) :: ix, iy, iz, ieg, ipert
      ! internal
            integer :: e, ijke, i
            real(dp) :: q
            e = gllel(ieg)
            ijke = ix + lx1*((iy - 1) + ly1*((iz - 1) + lz1*(e - 1)))
            q = neklab_fq(ijke, ipert + 1)
            if (ipert == 0) then
               do i = 1, nq_terms
                  call q_terms(i)%fn(q, ix, iy, iz, e)
               end do
            end if
            qvol = qvol + q*vtrans(ix, iy, iz, e, 2)
         end subroutine neklab_source

         subroutine add_neklab_forcing_term(fn, name)
      !! Register a pointwise forcing term (ipert = 0 only). A term with the same
      !! name is replaced, so registration is idempotent.
            procedure(pointwise_forcing_term) :: fn
            character(len=*), intent(in) :: name
            integer :: i
            do i = 1, nf_terms
               if (f_terms(i)%name == name) then
                  f_terms(i)%fn => fn
                  return
               end if
            end do
            if (nf_terms == max_terms) call nek_stop_error('too many forcing terms; raise max_terms.',
     &                                                     this_module, 'add_neklab_forcing_term')
            nf_terms = nf_terms + 1
            f_terms(nf_terms)%name = name
            f_terms(nf_terms)%fn => fn
            call nek_log_message('forcing term registered: '//trim(name), this_module,
     &                           procedure='add_neklab_forcing_term')
         end subroutine add_neklab_forcing_term

         subroutine add_neklab_source_term(fn, name)
      !! Register a pointwise source term for the temperature slot (ipert = 0 only).
            procedure(pointwise_source_term) :: fn
            character(len=*), intent(in) :: name
            integer :: i
            do i = 1, nq_terms
               if (q_terms(i)%name == name) then
                  q_terms(i)%fn => fn
                  return
               end if
            end do
            if (nq_terms == max_terms) call nek_stop_error('too many source terms; raise max_terms.',
     &                                                     this_module, 'add_neklab_source_term')
            nq_terms = nq_terms + 1
            q_terms(nq_terms)%name = name
            q_terms(nq_terms)%fn => fn
            call nek_log_message('source term registered: '//trim(name), this_module,
     &                           procedure='add_neklab_source_term')
         end subroutine add_neklab_source_term

         subroutine remove_neklab_term(name)
      !! Remove the forcing and/or source term registered under this name.
            character(len=*), intent(in) :: name
            integer :: i, j
            j = 0
            do i = 1, nf_terms
               if (f_terms(i)%name /= name) then
                  j = j + 1
                  f_terms(j) = f_terms(i)
               end if
            end do
            nf_terms = j
            j = 0
            do i = 1, nq_terms
               if (q_terms(i)%name /= name) then
                  j = j + 1
                  q_terms(j) = q_terms(i)
               end if
            end do
            nq_terms = j
         end subroutine remove_neklab_term

         subroutine check_ipert(ipert, caller)
            integer, intent(in) :: ipert
            character(len=*), intent(in) :: caller
            character(len=128) :: msg
            if (ipert < 0 .or. ipert > lpert) then
               write (msg, '(A,I0)') 'Invalid value for ipert specified! ipert = ', ipert
               call nek_stop_error(msg, this_module, caller)
            end if
         end subroutine check_ipert
      end module neklab_nek_forcing
