      submodule(neklab_linops) pressure_projection_operator
         implicit none
      contains
         module procedure pressure_projection_init
         ! internal
         character(len=*), parameter :: this_procedure = 'pressure_projection_init'
         character(len=128) :: msg
         integer :: ntot
         ntot = lx1*ly1*lz1*nelv
         self%baseflow         = bf
         self%betaz            = betaz
         self%integration_type = integration_type
         call sethlm(self%h1,self%h2,self%integration_type)
         call copy(self%h2inv, h2inv, ntot)
         self%is_initialized = .true.
         write(msg,'(A,I6)')'Setup complete for intype =', integration_type
         call nek_log_message(msg, this_module, this_procedure)
         end procedure pressure_projection_init
      
         module procedure pressure_projection_matvec
         ! internal
         character(len=*), parameter :: this_procedure = 'pressure_projection_matvec'
         real(dp), dimension(lx1*ly1*lz1) :: wrk1, wrk2
         real(dp), dimension(lx1,ly1,lz1,lelv) :: wdivm1
         real(dp), dimension(lx2,ly2,lz2,lelv) :: wdivm2
         integer :: ie, ntot1, ntot2
         ntot1 = lx1*ly1*lz1*nelv
         ntot2 = lx2*ly2*lz2*nelv
         select type (vec_in)
         type is (nekp_dvector)
            select type (vec_out)
            type is (nekp_dvector)
               if (.not.self%is_initialized) then
                  call nek_stop_error('pressure_projection operator not set up.', this_module, this_procedure)
               end if
               ! Regular Ax = D (h2*B)^(-1) D^T
               call cdabdtp(vec_out%pr, vec_in%pr, self%h1, self%h2, self%h2inv, self%integration_type)
               ! add 3rd perturbation component. 
               ! we want to subtract the contribution
               !
               ! B * d/dz (h2*B)^(-1) d/dz wp = (-i betaz) (h2)^(-1) (-i betaz) wp = - betaz^2 (h2)^(-1) wp
               !
               ! Note: we are only considering the ifanls = .false. case
               call mappr  (wdivm1, vec_in%pr, wrk1, wrk2)
               call col2   (wdivm1, bm1,   ntot1)
               call col2   (wdivm1, wmask, ntot1)     ! not v3mask
               call dssum  (wdivm1, lx1, ly1, lz1)
               call col2   (wdivm1, binvm1, ntot1)
               call col2   (wdivm1, self%h2inv, ntot1)
               call cmult  (wdivm1, -self%betaz**2, ntot1)                    ! make continuous
               do ie = 1, nelv
                  call map12 (wdivm2(1,1,1,ie), wdivm1(1,1,1,ie), ie)  ! map wdiv to pmesh
               end do
               call col2(wdivm2, bm2, ntot2)
               call sub2(vec_out%pr, wdivm2, ntot2)
            class default
               call type_error('vec_out','nekp_dvector','OUT',this_module,this_procedure)
            end select
         class default
            call type_error('vec_in','nekp_dvector','IN',this_module,this_procedure)
         end select
         end procedure pressure_projection_matvec
   
         !module procedure apply_hsmg_preconditioner
         !! internal
         !character(len=*), parameter :: this_procedure = 'apply_hsmg_preconditioner'  
         !integer :: ntot
         !ntot = lx1*ly1*lz1*nelv
         !select type (vec)
         !type is (nekp_dvector)
         !   if (.not.self%is_initialized) then
         !      call nek_stop_error('Preconditioner not set up.', this_module, this_procedure)
         !   end if
         !   call col2(vec%v,self%D,ntot)
         !   !call copy(vec%v,tmp,ntot)
         !class default
         !   call type_error('vec','nekv_dvector','INOUT',this_module,this_procedure)
         !end select
         !end procedure apply_hsmg_preconditioner
      end submodule
