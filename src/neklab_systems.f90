      module neklab_systems
      ! Fortran standard Library
         use stdlib_stats_distribution_normal, only: normal => rvs_normal
         use stdlib_optval, only: optval
         use stdlib_linalg, only: svd, diag
      !---------------------------------------
      !-----     LightKrylov Imports     -----
      !---------------------------------------
      ! Default real kind.
         use LightKrylov, only: dp, qr, apply_inverse_permutation_matrix
      ! Abstract types for real-valued vectors and utilities
         use LightKrylov, only: abstract_linop_rdp, abstract_vector_rdp
         use LightKrylov, only: orthonormalize_basis, orthogonalize_against_basis, zero_basis, copy, rand_basis, linear_combination
         use LightKrylov, only: abstract_system_rdp, abstract_jacobian_linop_rdp
         use LightKrylov_Logger
         use LightKrylov_utils, only: assert_shape
      ! Extensions of the abstract vector types to nek data format.
         use neklab_vectors
         use neklab_linops
         use neklab_utils
         use neklab_nek_setup
         implicit none
         include "SIZE"
         include "TOTAL"
         include "ADJOINT"
         private
         character(len=*), parameter, private :: this_module = 'neklab_systems'
      
         integer, parameter :: lv = lx1*ly1*lz1*lelv
      !! Local number of grid points for the velocity mesh.
         integer, parameter :: lp = lx2*ly2*lz2*lelv
      !! Local number of grid points for the pressure mesh.

         type(nek_ext_dvector_forcing) :: F_residual
         real(dp)                      :: G_residual
         real(dp), parameter           :: delta = 0.3_dp

         public :: get_mass_flux
      
      !--------------------------------------------------------------
      !-----     NEK SYSTEM FOR REAL VECTORS (fixed-points)   -------
      !--------------------------------------------------------------
      
         type, extends(abstract_system_rdp), public :: nek_system
         contains
            private
            procedure, pass(self), public :: response => nonlinear_map
         end type nek_system

         type, extends(abstract_jacobian_linop_rdp), public :: nek_jacobian
         contains
            private
            procedure, pass(self), public :: matvec => jac_exptA_matvec
            procedure, pass(self), public :: rmatvec => jac_exptA_rmatvec
         end type nek_jacobian

      !--------------------------------------------------------------------------
      !-----     NEK SYSTEM FOR REAL EXTENDED VECTORS (periodic orbits)   -------
      !--------------------------------------------------------------------------

         type, extends(abstract_system_rdp), public :: nek_system_upo
         contains
            private
            procedure, pass(self), public :: response => nonlinear_map_upo
         end type nek_system_upo
      
         type, extends(abstract_jacobian_linop_rdp), public :: nek_jacobian_upo
         contains
            private
            procedure, pass(self), public :: matvec => jac_direct_map
            procedure, pass(self), public :: rmatvec => jac_adjoint_map
         end type nek_jacobian_upo

      !-----------------------------------------------------------------------------
      !-----     NEK SYSTEM FOR REAL EXTENDED VECTORS (forced fixed-point)   -------
      !-----------------------------------------------------------------------------

         type, extends(abstract_system_rdp), public :: nek_system_forcing
         contains
            private
            procedure, pass(self), public :: response => nonlinear_map_forcing
         end type nek_system_forcing

         type, extends(abstract_jacobian_linop_rdp), public :: nek_jacobian_forcing
         contains
            private
            procedure, pass(self), public :: matvec => jac_direct_map_forcing
            procedure, pass(self), public :: rmatvec => jac_adjoint_map_forcing
         end type nek_jacobian_forcing
      
      contains
      
      !-------------------------------------------------------
      !-----     TYPE BOUND PROCEDURES FOR NEK_SYSTEM    -----
      !-------------------------------------------------------
      
         subroutine nonlinear_map(self, vec_in, vec_out, atol)
      ! Dynamical system.
            class(nek_system), intent(inout) :: self
      ! Input vector.
            class(abstract_vector_rdp), intent(in) :: vec_in
      ! Output vector.
            class(abstract_vector_rdp), intent(out) :: vec_out
      ! Solver tolerances if needed
            real(dp), intent(in) :: atol
      
            select type (vec_in)
            type is (nek_dvector)
               select type (vec_out)
               type is (nek_dvector)
      ! Set the initial condition
                  call vec2nek(vx, vy, vz, pr, t, vec_in)
      
      ! Set appropriate tolerances
                  call setup_nonlinear_solver(recompute_dt=.true., cfl_limit=0.4_dp,
     $   vtol = atol/10.0, ptol = atol/10.0)
      
      ! Intgrate the nonlinear equations forward
                  time = 0.0_dp
                  do istep = 1, nsteps
                     call nek_advance()
                  end do
      
      ! Extract the final solution to vector.
                  call nek2vec(vec_out, vx, vy, vz, pr, t)
      
      ! Evaluate residual F(X) - X.
                  call vec_out%sub(vec_in)
               end select
            end select
      
            return
         end subroutine nonlinear_map
      
      !---------------------------------------------------------
      !-----     TYPE BOUND PROCEDURES FOR NEK_JACOBIAN    -----
      !---------------------------------------------------------
      
         subroutine jac_exptA_matvec(self, vec_in, vec_out)
            class(nek_jacobian), intent(inout) :: self
            class(abstract_vector_rdp), intent(in) :: vec_in
            class(abstract_vector_rdp), intent(out) :: vec_out
      
            select type (vec_in)
            type is (nek_dvector)
               select type (vec_out)
               type is (nek_dvector)
      ! Set the baseflow initial condition
                  call abs_vec2nek(vx, vy, vz, pr, t, self%X)
      
      ! Ensure correct nek status
                  call setup_linear_solver(solve_baseflow=.false.,
     $   recompute_dt = .true., cfl_limit = 0.5_dp)
      
      ! Set the initial condition for Nek5000's linearized solver.
                  call vec2nek(vxp, vyp, vzp, prp, tp, vec_in)
      
      ! Integrate the equations forward in time.
                  time = 0.0_dp
                  do istep = 1, nsteps
                     call nek_advance()
                  end do
      
      ! Extract the final solution to vector.
                  call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)
      
      ! Evaluate [ exp(tau*J) - I ] @ dx.
                  call vec_out%sub(vec_in)
               end select
            end select
      
            return
         end subroutine jac_exptA_matvec
      
         subroutine jac_exptA_rmatvec(self, vec_in, vec_out)
            class(nek_jacobian), intent(inout) :: self
            class(abstract_vector_rdp), intent(in) :: vec_in
            class(abstract_vector_rdp), intent(out) :: vec_out
      
            select type (vec_in)
            type is (nek_dvector)
               select type (vec_out)
               type is (nek_dvector)
      ! Set the baseflow initial condition
                  call abs_vec2nek(vx, vy, vz, pr, t, self%X)
      
      ! Ensure correct nek status
                  call setup_linear_solver(transpose=.true., solve_baseflow=.false.,
     $   recompute_dt = .true., cfl_limit = 0.5_dp)
      
      ! Set the initial condition for Nek5000's linearized solver.
                  call vec2nek(vxp, vyp, vzp, prp, tp, vec_in)
      
      ! Integrate the equations forward in time.
                  time = 0.0_dp
                  do istep = 1, nsteps
                     call nek_advance()
                  end do
      
      ! Extract the final solution to vector.
                  call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)
      
      ! Evaluate [ exp(tau*J) - I ] @ dx.
                  call vec_out%sub(vec_in)
               end select
            end select
      
            return
         end subroutine jac_exptA_rmatvec
      
      !-----------------------------------------------------------
      !-----     TYPE BOUND PROCEDURES FOR NEK_SYSTEM_UPO    -----
      !-----------------------------------------------------------
      
         subroutine nonlinear_map_UPO(self, vec_in, vec_out, atol)
      ! Dynamical system.
            class(nek_system_upo), intent(inout) :: self
      ! Input vector.
            class(abstract_vector_rdp), intent(in) :: vec_in
      ! Output vector.
            class(abstract_vector_rdp), intent(out) :: vec_out
      ! Solver tolerances if needed
            real(dp), intent(in) :: atol
      ! internal
            character(len=128) :: msg
      
            select type (vec_in)
            type is (nek_ext_dvector)
               select type (vec_out)
               type is (nek_ext_dvector)
      ! Set the initial condition
                  call ext_vec2nek(vx, vy, vz, pr, t, vec_in)
      
      ! Set appropriate tolerances and Nek status
                  call setup_nonlinear_solver(recompute_dt=.true., endtime=vec_in%T,
     $   cfl_limit = 0.4_dp, vtol = atol/10.0, ptol = atol/10.0)
      
                  write (msg, '(A,F9.6)') 'Current period estimate, T = ', vec_in%T
                  if (nid == 0) print *, msg
                  call logger%log_message(msg, module=this_module, procedure='nonlinear_map_UPO')
      
      ! Intgrate the nonlinear equations forward
                  time = 0.0_dp
                  do istep = 1, nsteps
                     call nek_advance()
                  end do
      
      ! Copy the final solution to vector.
                  call nek2ext_vec(vec_out, vx, vy, vz, pr, t)
                  vec_out%T = vec_in%T
      
      ! Evaluate residual F(X) - X.
                  call vec_out%sub(vec_in)
               end select
            end select
      
            return
         end subroutine nonlinear_map_UPO
      
      !-------------------------------------------------------------
      !-----     TYPE BOUND PROCEDURES FOR NEK_JACOBIAN_UPO    -----
      !-------------------------------------------------------------
      
         subroutine jac_direct_map(self, vec_in, vec_out)
      ! Dynamical system.
            class(nek_jacobian_upo), intent(inout) :: self
      ! Input vector.
            class(abstract_vector_rdp), intent(in) :: vec_in
      ! Output vector.
            class(abstract_vector_rdp), intent(out) :: vec_out
      ! internal
            type(nek_ext_dvector) :: vec
      
            select type (vec_in)
            type is (nek_ext_dvector)
               select type (vec_out)
               type is (nek_ext_dvector)
      ! Set the baseflow initial condition
                  call abs_ext_vec2nek(vx, vy, vz, pr, t, self%X)
      
      ! Ensure correct nek status -> set end time
                  call setup_linear_solver(solve_baseflow=.true., transpose=.false.,
     $   recompute_dt = .true., endtime = get_period_abs(self%X), cfl_limit = 0.4_dp)
      
      ! Set the perturbation initial condition
                  call ext_vec2nek(vxp, vyp, vzp, prp, tp, vec_in)
      
      ! Intgrate the coupled equations forward
                  time = 0.0_dp
                  do istep = 1, nsteps
                     call nek_advance()
                  end do
      
      ! Copy the final solution to vector.
                  call nek2ext_vec(vec_out, vxp, vyp, vzp, prp, tp)
      
      ! Evaluate [ exp(tau*J) - I ] @ dx.
                  call vec_out%sub(vec_in)
      
      ! Evaluate f'(X(T), T) * dT and add it to the position residual
      ! Here we assume that vx,vy,vz contains the endpoint of the nonlinear trajectory
                  call compute_fdot(vec)
                  call vec_out%axpby(1.0_dp, vec, vec_in%T)
      
      ! Evaluate f'(X(0), 0).T @ dx and add phase condition
      ! Set the initial point of the nonlinear trajectory
                  call abs_ext_vec2nek(vx, vy, vz, pr, t, self%X)
                  call compute_fdot(vec)
                  vec_out%T = vec_in%dot(vec)
               end select
            end select
      
            return
         end subroutine jac_direct_map
      
         subroutine jac_adjoint_map(self, vec_in, vec_out)
            class(nek_jacobian_upo), intent(inout) :: self
            class(abstract_vector_rdp), intent(in) :: vec_in
            class(abstract_vector_rdp), intent(out) :: vec_out
      ! internal
            type(nek_ext_dvector) :: vec
      
            select type (vec_in)
            type is (nek_ext_dvector)
               select type (vec_out)
               type is (nek_ext_dvector)
      ! Set the baseflow initial condition
                  call abs_ext_vec2nek(vx, vy, vz, pr, t, self%X)
      
      ! Ensure correct nek status -> set end time
                  call setup_linear_solver(solve_baseflow=.true., transpose=.true.,
     $   recompute_dt = .true., endtime = get_period_abs(self%X), cfl_limit = 0.4_dp)
      
      ! Set the perturbation initial condition
                  call ext_vec2nek(vxp, vyp, vzp, prp, tp, vec_in)
      
      ! Integrate the equations forward in time.
                  time = 0.0_dp
                  do istep = 1, nsteps
                     call nek_advance()
                  end do
      
      ! Copy the final solution to vector.
                  call nek2ext_vec(vec_out, vxp, vyp, vzp, prp, tp)
      
      ! Evaluate [ exp(tau*J) - I ] @ dx.
                  call vec_out%sub(vec_in)
      
      ! Evaluate f'(X(T), T) * dT and add it to the position residual
                  call compute_fdot(vec)
                  call vec_out%axpby(1.0_dp, vec, vec_in%T)
      
      ! Evaluate f'(X(0), 0).T @ dx and add phase condition
      ! Set the initial point of the orbit
                  call abs_ext_vec2nek(vx, vy, vz, pr, t, self%X)
                  call compute_fdot(vec)
                  vec_out%T = vec_in%dot(vec)
               end select
            end select
      
            return
         end subroutine jac_adjoint_map
      
         subroutine compute_fdot(vec)
            class(nek_ext_dvector), intent(out) :: vec
      ! internal
            type(nek_ext_dvector) :: vec_in
      ! copy initial condition
            call nek2ext_vec(vec_in, vx, vy, vz, pr, t)
      ! Integrate the nonlinear equations forward in time.
            time = 0.0_dp
            do istep = 1, 1
               call nek_advance()
            end do
      ! Extract f(X(t0+dt))
            call nek2ext_vec(vec, vx, vy, vz, pr, t)
      ! Approximate derivative at t = t0:
      !
      !    f'(X(t0)) ~ ( f(X(t0+dt)) - f(X(t0)) ) / dt
      !
            call vec%sub(vec_in)
            call vec%scal(1.0/dt)
            vec%T = 0.0_dp ! ensure that the period shift vec_out%T is zero
            return
         end subroutine compute_fdot

      !---------------------------------------------------------------
      !-----     TYPE BOUND PROCEDURES FOR NEK_SYSTEM_FORCING    -----
      !---------------------------------------------------------------
      
         subroutine nonlinear_map_forcing(self, vec_in, vec_out, atol)
      ! Dynamical system.
            class(nek_system_forcing), intent(inout) :: self
      ! Input vector.
            class(abstract_vector_rdp), intent(in) :: vec_in
      ! Output vector.
            class(abstract_vector_rdp), intent(out) :: vec_out
      ! Solver tolerances if needed
            real(dp), intent(in) :: atol
      ! internal
            character(len=128) :: msg
      
            select type (vec_in)
            type is (nek_ext_dvector_forcing)
               select type (vec_out)
               type is (nek_ext_dvector_forcing)
      ! Set the initial condition
                  call ext_vec_f2nek(vx, vy, vz, pr, t, vec_in)
      
      ! Set appropriate tolerances and Nek status
                  call setup_nonlinear_solver(recompute_dt=.true.,
     $   cfl_limit = 0.4_dp, vtol = atol/10.0, ptol = atol/10.0)
      
      ! Intgrate the nonlinear equations forward
                  time = 0.0_dp
                  do istep = 1, nsteps
                     call nek_advance()
                  end do
      
      ! Copy the final solution to vector.
                  call nek2ext_vec_f(vec_out, vx, vy, vz, pr, t)
                  vec_out%f = vec_in%f

      ! Copy to scratch
                  call copy(F_residual, vec_out)
                  F_residual%f = vec_in%f
      
      ! Evaluate residual F(X) - X.
                  call vec_out%sub(vec_in)

      ! Copy to scratch
                  call F_residual%sub(vec_in)
      ! Compute G_residual
                  G_residual = get_mass_flux(F_residual)
                  
               end select
            end select
      
            return
         end subroutine nonlinear_map_forcing
      
      !-------------------------------------------------------------
      !-----     TYPE BOUND PROCEDURES FOR NEK_JACOBIAN_UPO    -----
      !-------------------------------------------------------------
      
         subroutine jac_direct_map_forcing(self, vec_in, vec_out)
      ! Dynamical system.
            class(nek_jacobian_forcing), intent(inout) :: self
      ! Input vector.
            class(abstract_vector_rdp), intent(in) :: vec_in
      ! Output vector.
            class(abstract_vector_rdp), intent(out) :: vec_out
      ! internal
            type(nek_ext_dvector_forcing) :: dFdf
            real(dp) :: dGdf, df
            real(dp), parameter :: mean = 0.0_dp
            real(dp), parameter :: var  = 1.0_dp
      
            select type (vec_in)
            type is (nek_ext_dvector_forcing)
               select type (vec_out)
               type is (nek_ext_dvector_forcing)
      ! Set the baseflow initial condition
                  call abs_ext_vec_f2nek(vx, vy, vz, pr, t, self%X)
      
      ! Ensure correct nek status -> set end time
                  call setup_linear_solver(solve_baseflow=.true., transpose=.false.,
     $   recompute_dt = .true., cfl_limit = 0.4_dp)
      
      ! Set the perturbation initial condition
                  call ext_vec_f2nek(vxp, vyp, vzp, prp, tp, vec_in)
      
      ! Intgrate the equations forward
                  time = 0.0_dp
                  do istep = 1, nsteps
                     call nek_advance()
                  end do
      
      ! Copy the final solution to vector.
                  call nek2ext_vec_f(vec_out, vxp, vyp, vzp, prp, tp)
      
      ! Evaluate [ exp(tau*J) - I ] @ dx.
                  call vec_out%sub(vec_in)
           
      ! Evaluate G'(X) @ dx and add to forcing term
                  vec_out%f = get_mass_flux(vec_in)

      ! Evaluate dF/df and dG/df
      ! Set the initial guess for fixed point
                  call abs_ext_vec_f2nek(vx, vy, vz, pr, t, self%X)
      ! Add perturbation to f
                  df = normal(mean, var)*1e-6
                  call compute_f_derivatives(dFdf, dGdf, df)
                  call vec_out%axpby(1.0_dp, dFdf, vec_in%f)

                  vec_out%f = vec_out%f + dGdf*vec_in%f

               end select
            end select
      
            return
         end subroutine jac_direct_map_forcing
      
         subroutine jac_adjoint_map_forcing(self, vec_in, vec_out)
            class(nek_jacobian_forcing), intent(inout) :: self
            class(abstract_vector_rdp), intent(in) :: vec_in
            class(abstract_vector_rdp), intent(out) :: vec_out
      ! internal
            type(nek_ext_dvector_forcing) :: dFdf
            real(dp) :: dGdf, df
            real(dp), parameter :: mean = 0.0_dp
            real(dp), parameter :: var  = 1.0_dp
      
            select type (vec_in)
            type is (nek_ext_dvector_forcing)
               select type (vec_out)
               type is (nek_ext_dvector_forcing)
      ! Set the baseflow initial condition
                  call abs_ext_vec_f2nek(vx, vy, vz, pr, t, self%X)
      
      ! Ensure correct nek status -> set end time
                  call setup_linear_solver(solve_baseflow=.true., transpose=.true.,
     $   recompute_dt = .true., cfl_limit = 0.4_dp)
      
      ! Set the perturbation initial condition
                  call ext_vec_f2nek(vxp, vyp, vzp, prp, tp, vec_in)
      
      ! Integrate the equations forward in time.
                  time = 0.0_dp
                  do istep = 1, nsteps
                     call nek_advance()
                  end do
      
      ! Copy the final solution to vector.
                  call nek2ext_vec_f(vec_out, vxp, vyp, vzp, prp, tp)
      
      ! Evaluate [ exp(tau*J) - I ] @ dx.
                  call vec_out%sub(vec_in)
           
      ! Evaluate G'(X) @ dx and add to forcing term
                  vec_out%f = get_mass_flux(vec_in)

      ! Evaluate dF/df and dG/df
      ! Set the initial guess for fixed point
                  call abs_ext_vec_f2nek(vx, vy, vz, pr, t, self%X)
      ! Add perturbation to f
                  df = normal(mean, var)*1e-6
                  ! add forcing to baseflow
                  call compute_f_derivatives(dFdf, dGdf, df)
                  call vec_out%axpby(1.0_dp, dFdf, vec_in%f)

                  vec_out%f = vec_out%f + dGdf*vec_in%f
               end select
            end select
      
            return
         end subroutine jac_adjoint_map_forcing

         subroutine compute_f_derivatives(dFdf, dGdf, df)
            class(nek_ext_dvector_forcing), intent(out) :: dFdf
            real(dp),                       intent(out) :: dGdf
            real(dp),                       intent(in)  :: df
      ! internal
            type(nek_ext_dvector_forcing) :: vec_in
      ! copy initial condition
            call nek2ext_vec_f(vec_in, vx, vy, vz, pr, t)
      ! Set nonlinear solver
            call setup_nonlinear_solver(recompute_dt=.true., cfl_limit = 0.4_dp)
      ! Integrate the nonlinear equations forward in time.
            time = 0.0_dp
            do istep = 1, nsteps
               call nek_advance()
            end do
      ! Extract F(X,f+df)
            call nek2ext_vec_f(dFdf, vx, vy, vz, pr, t)
            call dFdf%sub(vec_in)
      !
      ! Approximate mass flux derivative
      !
      !    dG/df ~ ( G(X(U,f + df)) - G(X(U,f)) ) / df
      !
            dGdf = ( get_mass_flux(dFdf) - G_residual ) / df
      !
      ! Approximate derivative
      !
      !    dF/df ~ ( F(X(U,f + df)) - F(X(U,f)) ) / df
      !
            call dFdf%sub(F_residual)
            call dFdf%scal(1.0/df)
            dFdf%f = 0.0_dp ! ensure that the forcing shift vec_out%f is zero
            return
         end subroutine compute_f_derivatives

         real(dp) function get_mass_flux(vec) result(mass_flux)
            implicit none
            type(nek_ext_dvector_forcing), intent(in) :: vec
            ! internal
            real(dp), dimension(lx1, ly1, lz1, lelt) :: ux, uy, uz
            real(dp), dimension(lx2, ly2, lz2, lelt) :: pres
            real(dp), dimension(lx1, ly1, lz1, lelt, ldimt) :: tmp
            real(dp) :: num, den, angle_t, an_phi
            real(dp) :: x, y, z, p, pp, alpha, rr, rinv, us, us_rinv
            integer  :: i, j, k, e, eg
            ! functions
            real(dp), external :: glsum

            pi = atan(1.0_dp)*4.0_dp

            num = 0.0_dp
            den = 0.0_dp
            an_phi = 0.0_dp

            call ext_vec_f2nek(ux, uy, uz, pres, tmp, vec)

            do e = 1, nelv
               eg = lglel(e)
               do k = 1,lz1
                  do j = 1,ly1
                     do i = 1,lx1
                        x  = xm1(i,j,k,e)
                        y  = ym1(i,j,k,e)
                        z  = zm1(i,j,k,e)
                        pp = x*x + y*y
                        p  = sqrt(pp) - 1.0_dp/delta
                        alpha = atan2(p,z)
                        rr = p*p + z*z
                        angle_t = atan2(x,y) ! clockwise from y axis
                        rinv = 1.0_dp/(1+delta*sqrt(rr)*sin(alpha))

                        us = cos(an_phi)*(cos(angle_t)*ux(i,j,k,e)
     $                       -sin(angle_t)*uy(i,j,k,e))+sin(an_phi)*uz(i,j,k,e)

                        us_rinv = us*rinv      ! Streamwise u/r
                        num = num + us_rinv*bm1(i,j,k,e)
                        den = den +    rinv*bm1(i,j,k,e)
                     enddo
                  enddo
               enddo
            enddo
            num = glsum(num,1)
            den = glsum(den,1)
            mass_flux = num/den  ! "1/r"-weighted volumetric average of streamwise velocity
         end function get_mass_flux
      
      end module neklab_systems
