      module neklab_2Dh_common
      !! Infrastructure shared by the 2Dh solvers (cartesian: neklab_2Dh_planar,
      !! axisymmetric: neklab_2Dh_axisym):
      !!   - wmask: Dirichlet mask of the third velocity component, which is
      !!     carried in the temperature slot (w, or u_phi);
      !!   - the run checks common to both formulations;
      !!   - init_2Dh_planar(), called once from userchk at istep = 0, after the case's
      !!     own geometry setup.
      !! The solvers stop if init has not been called, instead of building
      !! their data lazily on the first time step.
         use stdlib_optval, only: optval
         use LightKrylov, only: dp
         use neklab_nek_setup, only: nek_log_message, nek_log_warning, nek_stop_error
         use neklab_timing, only: neklab_timer_start, neklab_timer_stop, t_gs_comm
         implicit none
         include "SIZE"
         include "TOTAL"

         private
         character(len=*), parameter, private :: this_module = 'neklab_2Dh_common'

         real(dp), dimension(lx1, ly1, lz1, lelv), public :: wmask
         logical, public :: wmask_defined = .false.
         logical, private :: init_2Dh_planar_done = .false.

         ! BC masks for the w component
         public :: build_wmask
         ! init and check routines
         public :: init_2Dh_planar
         public :: is_initialized_2Dh_planar
         public :: check_2Dh_properties
         public :: check_2Dh_masks

      contains

         subroutine init_2Dh_planar(verbose)
      !! Set up the cartesian 2Dh solver. Call from userchk at istep = 0, after
      !! usrdat2 (final mesh and mass matrices) and before any 2Dh linop is
      !! initialised. Idempotent.
            logical, optional, intent(in) :: verbose
            character(len=*), parameter :: this_procedure = 'init_2Dh_planar'

            if (init_2Dh_planar_done) return

            if (ifaxis) call nek_stop_error(
     &         'init_2Dh_planar is for the cartesian 2Dh solver, incompatible with ifaxis.',
     &         this_module, this_procedure)
            
            call build_wmask(verbose)
            
            call check_2Dh_masks(this_procedure)
            call check_2Dh_properties(this_procedure)

            init_2Dh_planar_done = .true.
            call nek_log_message('cartesian 2Dh solver initialised.', this_module, this_procedure)
         end subroutine init_2Dh_planar

         logical function is_initialized_2Dh_planar()
            is_initialized_2Dh_planar = init_2Dh_planar_done
         end function is_initialized_2Dh_planar

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
                  case ('W  ', 'v  ', 'V  ', 'vl ', 'VL ', 'mv ', 'MV ')
                     call facev(wmask, ie, ifc, 0.0_dp, lx1, ly1, lz1)
      ! Everything else leaves w free:
      !   'SYM'          - normal lies in the x-y plane, w is tangential
      !   'O  ','o  '    - outflow, natural BC
      !   'ON ','on '    - only the normal component is constrained
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
      !! BCs), while the deck sets the temperature BCs separately. Warn if the
      !! two disagree: the scalar BCs in the deck are then not the ones solved.
            character(len=*), intent(in) :: caller
            integer :: ix, iy, iz, ie, nbad(1)
            integer, external :: iglsum
            character(len=256) :: msg

            nbad(1) = 0
            do ie = 1, nelv
               do iz = 1, lz1
               do iy = 1, ly1
               do ix = 1, lx1
                  if (tmask(ix,iy,iz,ie,1) /= wmask(ix,iy,iz,ie)) nbad(1) = nbad(1) + 1
               end do
               end do
               end do
            end do
            nbad(1) = iglsum(nbad, 1)
            if (nbad(1) > 0) then
               write(msg,'(A,I10,A)') 'temperature mask differs from wmask at ', nbad(1),
     &            ' nodes: check the temperature BCs of the third velocity component.'
               call nek_log_warning(msg, this_module, caller)
            end if
         end subroutine check_2Dh_masks

      end module neklab_2Dh_common
