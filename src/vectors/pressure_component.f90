      submodule(neklab_vectors) pressure_component
         implicit none
      contains
      
      !-------------------------------
      !-----     CONSTRUCTOR     -----
      !-------------------------------
      
         module procedure construct_nekp_dvector
      ! Pressure.
         out%pr = pr
         end procedure
      
      !-----------------------------------------
      !-----     TYPE-BOUND PROCEDURES     -----
      !-----------------------------------------
      
         module procedure nekp_dzero
         call self%scal(0.0_dp)
         end procedure
      
         module procedure nekp_drand
         logical :: normalize
         integer :: i
         real(kind=dp) :: xl(ldim), fcoeff(3), alpha
         normalize = optval(ifnorm, .false.)
         call self%scal(0.0_dp)
         do i = 1, lp
            xl(1) = xm2(i, 1, 1, 1)
            xl(2) = ym2(i, 1, 1, 1)
            if (if3D) xl(ldim) = zm2(i, 1, 1, 1)
            call random_number(fcoeff); fcoeff = fcoeff*1.0e4_dp
            self%pr(i) = self%pr(i) + mth_rand(i, 1, 1, 1, xl, fcoeff)
         end do
         if (normalize) then
            alpha = self%norm()
            call self%scal(1.0_dp/alpha)
         end if
         end procedure
      
         module procedure nekp_dscal
         integer :: n2
         n2 = nx2*ny2*nz2*nelv
         call dscal(n2, alpha, self%pr, 1)
         end procedure
      
         module procedure nekp_daxpby
         integer :: n2
         n2 = nx2*ny2*nz2*nelv
         call self%scal(beta)
         select type (vec)
         type is (nekp_dvector)
            call add2s2(self%pr, vec%pr, alpha, n2)
         class default
            call type_error('vec','nekp_dvector','IN',this_module,'nekp_daxpby')
         end select
         end procedure
      
         module procedure nekp_ddot
         real(kind=dp), external :: glsc2
         integer :: n
         n = nx2*ny2*nz2*nelv
         select type (vec)
         type is (nekp_dvector)
            !alpha = glsc3(self%pr, vec%pr, bm2, lp)
            !alpha = glsc3(self%pr, bm2inv, vec%pr, lp)/volvm2    ! rnorm from convprn in navier1.f
            alpha = glsc2(self%pr, vec%pr, n)                    ! convprn in navier1.f
         class default
            call type_error('vec','nekp_dvector','IN',this_module,'nekp_ddot')
         end select
         end procedure
      
         module procedure nekp_dsize
         N = nx2*ny2*nz2*nelv
         end procedure
      
      end submodule
