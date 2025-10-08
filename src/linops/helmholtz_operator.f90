      submodule(neklab_linops) helmholtz_operator
         implicit none
      contains
         module procedure helmholtz_init
            self%integration_type = integration_type
            call sethlm(self%h1,self%h2,self%integration_type)
            self%is_initialized = .true.
         end procedure helmholtz_init
      
         module procedure helmholtz_matvec
         ! internal
         integer :: ntot
         ntot = lx1*ly1*lz1*nelv
         select type (vec_in)
         type is (nekv_dvector)
            select type (vec_out)
            type is (nekv_dvector)
               if (.not.self%is_initialized) call self%init(-1)
               call vec2nek(vx, vy, vz, pr, t, self%baseflow)
               call axhelm (vec_out%v,vec_in%v,self%h1,self%h2,self%imesh,self%isd)
               call dssum  (vec_out%v,lx1,ly1,lz1)
               call col2   (vec_out%v,vec_in%mask,ntot)
            class default
               call type_error('vec_out','nekv_dvector','OUT',this_module,'helmholtz_matvec')
            end select
         class default
            call type_error('vec_in','nekv_dvector','IN',this_module,'helmholtz_matvec')
         end select
         end procedure helmholtz_matvec

         module procedure construct_jacobi_preconditioner
            call setprec(self%D,h1,h2,imsh,isd)
            self%is_initialized = .true.
         end procedure construct_jacobi_preconditioner

         module procedure apply_jacobi_preconditioner
         ! internal
            real(dp), dimension(lx1,ly1,lz1,lelv) :: tmp
            integer :: ntot
            ntot = lx1*ly1*lz1*nelv
            select type (vec)
            type is (nekv_dvector)
               call col3(tmp,vec%v,self%D,ntot)
               call copy(vec%v,tmp,ntot)
            class default
               call type_error('vec','nekv_dvector','INOUT',this_module,'apply_jacobi_preconditioner')
            end select
         end procedure apply_jacobi_preconditioner
      end submodule
