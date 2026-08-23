      submodule(neklab_vectors) bordered_vectors
         implicit none
      contains

      !-------------------------------
      !-----     CONSTRUCTOR     -----
      !-------------------------------

         module procedure construct_nek_brd_dvector
         out%vx = vx; out%vy = vy
         if (present(vz)) then
            out%vz = vz
         else
            out%vz = 0.0_dp
         end if
         if (present(pr)) then
            out%pr = pr
         else
            out%pr = 0.0_dp
         end if
         if (present(theta)) then
            out%theta = theta
         else
            out%theta = 0.0_dp
         end if
         if (present(g)) then
            out%g = g
         else
            out%g = 0.0_dp
         end if
         end procedure

      !-----------------------------------------
      !-----     TYPE-BOUND PROCEDURES     -----
      !-----------------------------------------

         module procedure nek_brd_dzero
         call self%scal(0.0_dp)
         self%g = 0.0_dp
         end procedure

         module procedure nek_brd_drand
      !! Randomises the state block exactly as nek_drand does, and the ACTIVE
      !! control slots only. Padding slots must stay identically zero.
         logical :: normalize
         integer :: ix, iy, iz, iel, ieg, ijke, i
         integer :: iface, kx1, kx2, ky1, ky2, kz1, kz2
         real(kind=dp) :: xl(ldim), fcoeff(3), alpha, r
         normalize = optval(ifnorm, .false.)

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
            self%vx(ijke) = self%vx(ijke) + mth_rand(ix, iy, iz, ieg, xl, fcoeff)

            call random_number(fcoeff); fcoeff = fcoeff*1.0e4_dp
            self%vy(ijke) = self%vy(ijke) + mth_rand(ix, iy, iz, ieg, xl, fcoeff)

            if (if3d) then
               call random_number(fcoeff); fcoeff = fcoeff*1.0e4_dp
               self%vz(ijke) = self%vz(ijke) + mth_rand(ix, iy, iz, ieg, xl, fcoeff)
            end if

            if (ifto) then
               call random_number(fcoeff); fcoeff = fcoeff*1.0e4_dp
               self%theta(ijke,1) = self%theta(ijke,1) + mth_rand(ix, iy, iz, ieg, xl, fcoeff)
            end if
         end do
         end do
         end do
         end do

         if (ifcyclic) then
            do iel = 1, nelv
               do iface = 1, 2*ndim
                  if (cbc(iface,iel,1) .eq. 'P  ') then
                     call facind(kx1, kx2, ky1, ky2, kz1, kz2, lx1, ly1, lz1, iface)
                     do iz = kz1, kz2
                     do iy = ky1, ky2
                     do ix = kx1, kx2
                        ijke = ix + lx1*((iy-1) + ly1*((iz-1) + lz1*(iel-1)))
                        self%vx(ijke) = 0.0_dp
                        self%vy(ijke) = 0.0_dp
                        if (if3d) self%vz(ijke) = 0.0_dp
                        if (ifto) self%theta(ijke,1) = 0.0_dp
                     end do
                     end do
                     end do
                  end if
               end do
            end do
         end if

         call opdssum(self%vx, self%vy, self%vz)
         call opcolv (self%vx, self%vy, self%vz, vmult)
         call dsavg  (self%vx)
         call dsavg  (self%vy)
         if (if3d) call dsavg(self%vz)
         if (ifto) then
            call dssum(self%theta, lx1, ly1, lz1)
            call col2 (self%theta, vmult, lx1*ly1*lz1*lelv)
            call dsavg(self%theta)
         end if

         ifield = 1
         call bcdirvc(self%vx, self%vy, self%vz, v1mask, v2mask, v3mask)
         if (ifto .and. ifaxis .and. ifaziv) then
            do iel = 1, nelv
               do iface = 1, 2*ndim
                  select case (cbc(iface,iel,1))
                  case ('W  ','v  ','V  ','vl ','VL ','mv ','MV ')
                     call facind(kx1, kx2, ky1, ky2, kz1, kz2, lx1, ly1, lz1, iface)
                     do iz = kz1, kz2
                     do iy = ky1, ky2
                     do ix = kx1, kx2
                        ijke = ix + lx1*((iy-1) + ly1*((iz-1) + lz1*(iel-1)))
                        self%theta(ijke,1) = 0.0_dp
                     end do
                     end do
                     end do
                  end select
               end do
            end do
         end if

      ! control block: active slots only, and identical on all ranks.
         self%g = 0.0_dp
         if (nid == 0) then
            do i = 1, nctrl
               call random_number(r)
               self%g(i) = 2.0_dp*r - 1.0_dp
            end do
         end if
         call bcast(self%g, lg*wdsize)

         if (normalize) then
            alpha = self%norm()
            call self%scal(1.0_dp/alpha)
         end if
         end procedure

         module procedure nek_brd_dscal
         integer :: n1, n2
         n1 = nx1*ny1*nz1*nelv
         n2 = nx2*ny2*nz2*nelv
         call           cmult(self%vx,          alpha, n1)
         call           cmult(self%vy,          alpha, n1)
         if (if3d) call cmult(self%vz,          alpha, n1)
         call           cmult(self%pr,          alpha, n2)
         if (ifto) call cmult(self%theta(:, 1), alpha, n1)
         self%g = alpha*self%g
         end procedure

         module procedure nek_brd_daxpby
         integer :: n1, n2
         n1 = nx1*ny1*nz1*nelv
         n2 = nx2*ny2*nz2*nelv
         call self%scal(beta)
         select type (vec)
         type is (nek_bordered_dvector)
            call           add2s2(self%vx,          vec%vx,          alpha, n1)
            call           add2s2(self%vy,          vec%vy,          alpha, n1)
            if (if3d) call add2s2(self%vz,          vec%vz,          alpha, n1)
            call           add2s2(self%pr,          vec%pr,          alpha, n2)
            if (ifto) call add2s2(self%theta(:, 1), vec%theta(:, 1), alpha, n1)
            call opdssum(self%vx, self%vy, self%vz)
            call opcolv (self%vx, self%vy, self%vz, vmult)
            call dsavg  (self%vx)
            call dsavg  (self%vy)
            if (if3d) call dsavg(self%vz)
            if (ifto) then
               call dssum(self%theta, lx1, ly1, lz1)
               call col2 (self%theta, vmult, lx1*ly1*lz1*lelv)
               call dsavg(self%theta)
            end if
            self%g = self%g + alpha*vec%g
         class default
            call type_error('vec','nek_bordered_dvector','IN',this_module,'nek_brd_daxpby')
         end select
         end procedure

         module procedure nek_brd_ddot
      !! Weighted inner product. The control block carries wg so that the two
      !! blocks are commensurate -- see set_control_weights_auto. Getting this
      !! wrong does not change the root, but it does change every Krylov
      !! convergence test and every Newton residual norm.
         real(kind=dp), external :: glsc3
         integer :: i, n
         n = nx1*ny1*nz1*nelv
         select type (vec)
         type is (nek_bordered_dvector)
            alpha =         glsc3(self%vx, vec%vx, bm1, n)
            alpha = alpha + glsc3(self%vy, vec%vy, bm1, n)
            if (if3d) alpha = alpha + glsc3(self%vz, vec%vz, bm1, n)
            if (ifto) then
               alpha = alpha + glsc3(self%theta(:, 1), vec%theta(:, 1), bm1, n)
            end if
            if (ldimt > 1) then
            do i = 2, ldimt
               if (ifpsco(i - 1)) alpha = alpha + glsc3(self%theta(:, i), vec%theta(:, i), bm1, n)
            end do
            end if
            do i = 1, nctrl
               alpha = alpha + wg(i)*self%g(i)*vec%g(i)
            end do
         class default
            call type_error('vec','nek_bordered_dvector','IN',this_module,'nek_brd_ddot')
         end select
         end procedure

         module procedure nek_brd_dsize
      !! ACTIVE controls only -- LightKrylov uses this for the
      !! "Krylov dimension exceeds problem size" guard.
         integer :: i, n1
         n1 = nx1*ny1*nz1*nelv
         n = 2*n1 + nx2*ny2*nz2*nelv
         if (if3d) n = n + n1
         if (ifto) n = n + n1
         if (ldimt > 1) then
         do i = 2, ldimt
            if (ifpsco(i - 1)) n = n + n1
         end do
         end if
         n = n + nctrl
         end procedure

      end submodule bordered_vectors