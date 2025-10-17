      submodule(neklab_linops) helmholtz_operator
         implicit none
      contains
         module procedure helmholtz_init
         ! internal
         character(len=*), parameter :: this_procedure = 'helmholtz_init'
         character(len=128) :: msg
         self%baseflow         = bf
         self%betaz            = betaz
         self%integration_type = integration_type
         self%imesh            = optval(imesh, 1)
         self%isd              = optval(isd, 1)
         call sethlm(self%h1,self%h2,self%integration_type)
         self%is_initialized = .true.
         write(msg,'(A,I6)')'Setup complete for intype =', integration_type
         call nek_log_message(msg, this_module, this_procedure)
         end procedure helmholtz_init
      
         module procedure helmholtz_matvec
         ! internal
         character(len=*), parameter :: this_procedure = 'helmholtz_matvec'
         real(dp), dimension(lx1,ly1,lz1,lelv) :: tmp
         integer :: ntot
         ntot = lx1*ly1*lz1*nelv
         select type (vec_in)
         type is (nekv_dvector)
            select type (vec_out)
            type is (nekv_dvector)
               if (.not.self%is_initialized) then
                  call nek_stop_error('Operator not set up.', this_module, this_procedure)
               end if
               ! set baseflow
               call vec2nek(vx, vy, vz, pr, t, self%baseflow)
               ! regular Ax
               call axhelm (vec_out%v,vec_in%v,self%h1,self%h2,self%imesh,self%isd)
               ! added spanwise diffusion
               call col3   (tmp, vec_in%v, vdiff(1,1,1,1,1), ntot)
               call col2c  (tmp, bm1, -self%betaz**2, ntot)
               call add2   (vec_out%v, tmp, ntot)
               ! make continuous
               call dssum  (vec_out%v, lx1, ly1, lz1) ! make continuous
               call col2   (vec_out%v, vec_in%mask, ntot)
            class default
               call type_error('vec_out','nekv_dvector','OUT',this_module,this_procedure)
            end select
         class default
            call type_error('vec_in','nekv_dvector','IN',this_module,this_procedure)
         end select
         end procedure helmholtz_matvec

         module procedure construct_jacobi_preconditioner
         call setprec(self%D,h1,h2,optval(imsh,1),optval(isd,1))
         self%is_initialized = .true.
         end procedure construct_jacobi_preconditioner

         module procedure apply_jacobi_preconditioner
         ! internal
         character(len=*), parameter :: this_procedure = 'apply_jacobi_preconditioner'  
         real(dp), dimension(lx1,ly1,lz1,lelv) :: tmp
         integer :: ntot
         ntot = lx1*ly1*lz1*nelv
         select type (vec)
         type is (nekv_dvector)
            if (.not.self%is_initialized) then
               call nek_stop_error('Preconditioner not set up.', this_module, this_procedure)
            end if
            call col2(vec%v,self%D,ntot)
         class default
            call type_error('vec','nekv_dvector','INOUT',this_module,this_procedure)
         end select
         end procedure apply_jacobi_preconditioner
      end submodule
