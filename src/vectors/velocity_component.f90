      submodule(neklab_vectors) velocity_component
         implicit none
      contains
      
      !-------------------------------
      !-----     CONSTRUCTOR     -----
      !-------------------------------
      
         module procedure construct_nekv_dvector
      ! Velocity array.
         out%v     = v
         out%mask  = mask
         out%vmult = vmult
         out%isd   = isd
         end procedure
      
      !-----------------------------------------
      !-----     TYPE-BOUND PROCEDURES     -----
      !-----------------------------------------
      
         module procedure nekv_dzero
         call self%scal(0.0_dp)
         end procedure
      
         module procedure nekv_drand
         logical :: normalize
         integer :: ix, iy, iz, iel, ieg, ijke
         integer :: iface, kx1, kx2, ky1, ky2, kz1, kz2, ntot
         real(kind=dp) :: xl(ldim), fcoeff(3), alpha
         ! internal
         real(kind=dp), dimension(lx1,ly1,lz1,lelv) :: dum1, dum2
         normalize = optval(ifnorm, .false.)
      
         ifield = 1 ! for bcdirvc
         ntot = lx1*ly1*lz1*nelv

         call self%scal(0.0_dp)
         do iel = 1, nelv
         do iz = 1, lz1
         do iy = 1, ly1
         do ix = 1, lx1
            ieg = lglel(iel)
            xl(1) = xm1(ix, iy, iz, iel)
            xl(2) = ym1(ix, iy, iz, iel)
            if (if3d) xl(3) = zm1(ix, iy, iz, iel)
            ijke = ix + lx1*((iy-1) + ly1*((iz-1) + lz1*(iel-1)))
            
            call random_number(fcoeff); fcoeff = fcoeff*1.0e4_dp
            self%v(ijke) = self%v(ijke) + mth_rand(ix, iy, iz, ieg, xl, fcoeff)
   
         end do
         end do
         end do
         end do

         ! zero out noise on the cyclic boundary to avoid spurious modes
         if (ifcyclic) then
            do iel = 1, nelv
               do iface = 1, 2*ndim
                  if (cbc(iface,iel,1) .eq. 'P  ') then
                     call facind(kx1, kx2, ky1, ky2, kz1, kz2, lx1, ly1, lz1, iface)
                     do iz = kz1, kz2
                     do iy = ky1, ky2
                     do ix = kx1, kx2
                        ijke = ix + lx1*((iy-1) + ly1*((iz-1) + lz1*(iel-1)))
                        self%v(ijke) = 0.0_dp
                     end do
                     end do
                     end do
                  end if
               end do
            end do
         end if
      
      ! Face averaging.
         call dssum(self%v, lx1, ly1, lz1)
         call col2(self%v, self%vmult, ntot)
         call dsavg(self%v)
         call bcdirvc(self%v, dum1, dum2, self%mask, self%mask, self%mask)
      
         if (normalize) then
            alpha = self%norm()
            call self%scal(1.0_dp/alpha)
         end if
         end procedure
      
         module procedure nekv_dscal
         integer :: n1, n2
         n1 = nx1*ny1*nz1*nelv
         call cmult(self%v, alpha, n1)
         end procedure
      
         module procedure nekv_daxpby
         integer :: n1
         n1 = nx1*ny1*nz1*nelv
         call self%scal(beta)
         select type (vec)
         type is (nekv_dvector)
            call add2s2(self%v, vec%v, alpha, n1)
         class default
            call type_error('vec','nekv_dvector','IN',this_module,'nekv_daxpby')
         end select
         end procedure
      
         module procedure nekv_ddot
         real(kind=dp), external :: glsc3
         integer :: i, n
         n = nx1*ny1*nz1*nelv
         select type (vec)
         type is (nekv_dvector)
            alpha = glsc3(self%v, vec%v, self%vmult, n)
         class default
            call type_error('vec','nekv_dvector','IN',this_module,'nekv_ddot')
         end select
         end procedure
      
         module procedure nekv_dsize
         n = nx1*ny1*nz1*nelv
         end procedure
      
      end submodule
