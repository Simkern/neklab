      submodule (neklab_helix) helix_io
         implicit none
            
      contains

         module procedure fname_2d
            ! internal
            character(len=*), parameter :: fmt = '(A,A,I3.3,A)'
            write(fname,fmt) iname, '2dtorus', iout, '.fld'
         end procedure fname_2d

         module procedure write_2d
            integer, allocatable :: n2d_gown(:)
            integer, allocatable :: n2d_elmap(:)
            integer :: ierr, itmp, i, nxy, ip, ibuf, iseg, length, i_own, nsave
            integer :: wdsl, isl, isend(lelv)
            character(len=128)  :: msg
            character(len=1024) :: head, ftm
            character(len=2) :: id
            character(len=1), parameter :: fileversion = '1'
            logical :: only_mesh_
            real rtmpv1(lx1*ly1*lelv), rtmpv(lx1*ly1*lelv)
            real*4 rtmpv2(2*lx1*ly1*lelv)
            equivalence (rtmpv1,rtmpv2)
            real*4 test
            parameter (test=6.54321)
            nxy = lx1*ly1
            wdsl = wdsize/4
            isl  = isize/4
            only_mesh_ = optval(only_mesh, .false.)
            if (only_mesh_) then
               write(msg,'(A,A)') 'Outpost 2D mesh: ', trim(fname)
               nsave = 0
            else
               write(msg,'(A,I5,4X,A,A)') 'Outpost 2D data: ', self%nsave, 'fname: ', trim(fname)
               nsave = self%nsave
            end if
            call nek_log_information(msg, this_module, 'outpost_2d_fields')
            if (nid == 0) then
               call byte_open(fname, ierr)
               if (ierr /= 0) call nek_stop_error('Error opening file '//trim(fname), procedure='outpost_2d_fields')

               ! write file's header
               ftm="('#',A1,A2,1x,i1,1x,'(lx1, ly1 =',2i9,') (nelf =',i9,') (time =',e17.9,') (nsave, lbuf =', 2i9,')')"
               if (self%is_sym()) then
                  id = 'th'
               else
                  id = 'tf'
               end if
               write(head,ftm) fileversion, id, wdsize,lx1,ly1,self%nelf,time,nsave,lbuf
               call byte_write(head,116/4,ierr)
               if (ierr /= 0) call nek_stop_error('Error writing header in file '//trim(fname), procedure='outpost_2d_fields')  

               ! write big/little endian test
               call byte_write(test,1,ierr)

               ! write metadata
               call byte_write(lx1,isl,ierr)
               call byte_write(ly1,isl,ierr)
               call byte_write(self%nelf,isl,ierr)
               call byte_write(time,wdsl,ierr)
               call byte_write(nsave,isl,ierr)
               call byte_write(lbuf,isl,ierr)
               if (ierr /= 0) call nek_stop_error('Error writing metadata in file '//trim(fname), procedure='outpost_2d_fields')
            end if

            ! gather information about elements on other procs  
            allocate(n2d_gown(np))           ! number of elements owned by each proc
            call izero(n2d_gown,np)
            allocate(n2d_elmap(self%nelf))  ! global element number of owned elements
            call izero(n2d_elmap,self%nelf)
            ! determine how many elements to dump
            if (nid == 0) then
               ! first for the master node 
               n2d_gown(1) = self%n2d_gown
               do i = 1, self%n2d_gown
                  n2d_elmap(i) = lglel(self%id2d(i,1))     ! get global element number
               end do 
               iseg = self%n2d_gown
               ! then gather info from other procs
               do ip = 1, np-1
                  call csend(ip,itmp,isize,ip,0)           ! hand shake
                  call crecv(ip,i_own,isize)               ! recv number of elements
                  n2d_gown(ip+1) = i_own
                  call crecv(ip,isend(:i_own),i_own*isize) ! recv global element map
                  n2d_elmap(iseg+1:iseg+i_own) = isend(:i_own)
                  iseg = iseg + i_own
               enddo
               if (iseg /= self%nelf) call nek_stop_error('Not all elements in slice found!', this_module, 'outpost_2d_fields')
               ! write it to file
               call byte_write(n2d_elmap,self%nelf*isl,ierr)
               ! write timestep information to file
               call byte_write(self%dt2d(:nsave),nsave*wdsl,ierr)
            else
               call crecv(nid,itmp,isize)                  ! hand shake
               call csend(nid,self%n2d_gown,isize,0,0)     ! send number of elements
               length = self%n2d_gown
               do i = 1, self%n2d_gown
                  isend(i) = lglel(self%id2d(i,1))       
               end do
               call csend(nid,isend(:length),length*isize,0,0)   ! send global element map
            endif
            call bcast(n2d_gown, np*isize)         ! broadcast to all procs
            call bcast(n2d_elmap, self%nelf*isize) ! broadcast to all procs

            ! coordinates
            call nek_log_debug('   '//trim(fname)//': write x2d ...', this_module, 'outpost_2d_fields')
            call gather_and_write_slice(self%x2d, n2d_gown)
            call nek_log_debug('   '//trim(fname)//': write y2d ...', this_module, 'outpost_2d_fields')
            call gather_and_write_slice(self%y2d, n2d_gown)
            ! velocity data
            if (.not. only_mesh_) then
               write(msg,'(3X,A,A,1X,I5)') trim(fname),': write v[xyz]2d', nsave
               call nek_log_debug(msg, this_module, 'outpost_2d_fields')
               do i = 1, nsave
                  call gather_and_write_slice(self%vx2d(:,:,:,i), n2d_gown)
                  call gather_and_write_slice(self%vy2d(:,:,:,i), n2d_gown)
                  call gather_and_write_slice(self%vz2d(:,:,:,i), n2d_gown)
               end do
            end if
            ! master closes the file
            if (nid == 0) then 
               call byte_close(ierr)
               if (ierr /= 0) call nek_stop_error('Error closing file '//trim(fname), procedure='outpost_2d_fields')
            end if
         end procedure write_2d

         module procedure read_2d
            ! only nid 0 will read
            integer ierr, hdrsize
            real*4 test_pattern
            integer :: nxr, nyr, nelfr, nsaver, lbufr, wdsizr, length
            integer :: wdsl, isl, nxy, i, ie, ieg, iel, iseg, gseg, nelf
            integer, allocatable :: global_map(:)
            integer, allocatable :: gmap_index(:)
            real(dp), allocatable :: slicedata(:,:,:,:)
            real(dp) :: dt2dr(lbuf)
            real(dp) :: timer
            character(len=132) :: hdr
            character(len=256) :: msg
            character(len=4)   :: sdummy
            common /CTMP1/ fldum(lx1*ly1*lelv)
            real fldum
            ! functions
            logical, external :: if_byte_swap_test
            hdrsize = 116
            nxy = lx1*ly1
            nelf = self%nelf
            allocate(global_map(nelf))
            allocate(gmap_index(nelf))
            if (nid == 0) then
               call byte_open(fname,ierr)
               if (ierr /= 0) call nek_stop_error('Error opening file '//trim(fname), procedure='load_2d_fields')
               ! read header
               if (ierr == 0) then
                  call blank     (hdr,hdrsize)
                  call byte_read (hdr,hdrsize/4,ierr)
               endif
               if (ierr == 0) then
                  call byte_read (test_pattern,1,ierr)
                  if_byte_sw = if_byte_swap_test(test_pattern,ierr) ! determine endianess
               endif
               call nek_log_debug('header: '//trim(hdr), this_module, 'load_2d_fields')
               ! read wdsize from header
               read(hdr,*) sdummy, wdsizr
               wdsl = wdsizr/4
               isl  = isize/4
               ! read metadata
               call byte_read(nxr,    isl, ierr)
               call byte_read(nyr,    isl, ierr)
               call byte_read(nelfr,  isl, ierr)
               call byte_read(timer, wdsl, ierr)
               call byte_read(nsaver, isl, ierr)
               call byte_read(lbufr,  isl, ierr)
               write(msg,'(A,3(1X,I0),1X,E15.7,2(1X,I0))') 'metadata: ', nxr, nyr, nelfr, timer, nsaver, lbufr
               call nek_log_debug(msg, this_module, 'load_2d_fields')
               if (lx1 /= nxr .or. ly1 /= nyr) call nek_stop_error('Reading '//trim(fname)//': Inconsistent lx1/ly1', procedure='read_2d')
               if (nelfr /= nelf) call nek_stop_error('Reading '//trim(fname)//': Inconsistent nelf', procedure='read_2d')
               if (lbufr /= lbuf) call nek_stop_error('Reading '//trim(fname)//': Inconsistent lbuf', procedure='read_2d')
               ! read global element mapping
               call byte_read(global_map, nelf*isl, ierr)
               if (ierr /= 0) call nek_stop_error('Error reading gloabl element map from file '//trim(fname), procedure='load_2d_fields')
               ! read timestep information
               call byte_read(dt2dr(:nsaver), nsaver*wdsl, ierr)
               self%dt2d = dt2dr
               if (ierr /= 0) call nek_stop_error('Error reading timestep information from file '//trim(fname), procedure='load_2d_fields')
               ! read coords but skip them
               call byte_read(fldum, nxy*nelf*wdsl, ierr)
               call byte_read(fldum, nxy*nelf*wdsl, ierr)
               if (ierr /= 0) call nek_stop_error('Error reading coordinates from file '//trim(fname), procedure='load_2d_fields')
            end if
            call bcast(nsaver, isize)          ! broadcast number of saved snapshots
            call bcast(self%dt2d, lbuf*wdsize) ! broadcast timestep data
            call bcast(global_map, nelf*isize) ! broadcast global element map
            call sort_index(global_map, gmap_index)
            ! initialize data and prepare arrays
            length = 3*nxy*nelf
            allocate(slicedata(lx1,ly1,nelf,3))
            call rzero(slicedata, length)
            ! load data one timestep at a time
            do i = 1, nsaver
               if (nid == 0) then ! read v[xyz]2d for all elements at the current timestep
                  call byte_read(slicedata, length*wdsl, ierr)
                  if (if_byte_sw) call byte_reverse(slicedata, length, ierr)
                  if (ierr /= 0) call stop_error('Error reading element data', procedure='load_and_distribute_slice')
               end if
               call bcast(slicedata, length*wdsize) ! broadcast 2D data to all procs
               ! distribute to local segment owners
               do iseg = 1, self%n2d_lown
                  gseg = self%id2d(iseg,3)  ! global segment
                  iel  = gmap_index(gseg)   ! get element that is read
                  call copy(self%vx2d(1,1,iseg,i), slicedata(1,1,iel,1), nxy)
                  call copy(self%vy2d(1,1,iseg,i), slicedata(1,1,iel,2), nxy)
                  call copy(self%vz2d(1,1,iseg,i), slicedata(1,1,iel,3), nxy)
               end do
            end do ! 1, nsaver
            ! master closes the file
            if (nid == 0) then 
               call byte_close(ierr)
               if (ierr /= 0) call nek_stop_error('Error closing file '//trim(fname), procedure='load_2d_fields')
            end if
            self%nload = nsaver
            write(msg,'(A,A,A,I0)') 'Loaded 2D data from file ', trim(fname), ': ', self%nload
            call nek_log_information(msg, this_module, 'load_2d_fields')
            call lk_timer%stop('neklab_helix_load_2d')
         end procedure read_2d

         module procedure get_nsteps_from_header
            ! only nid 0 will read
            integer ierr, hdrsize
            real*4 test_pattern
            integer :: nxr, nyr, nelfr, lbufr
            integer :: wdsl, isl
            real(dp) :: timer
            character(len=132) :: hdr, msg
            character(len=4)   :: sdummy
            ! functions
            logical, external :: if_byte_swap_test
            hdrsize = 116
            if (nid == 0) then
               call byte_open(fname,ierr)
               if (ierr /= 0) call nek_stop_error('Error opening file '//trim(fname), procedure='get_nsteps_from_header')
               ! read header
               if (ierr == 0) then
                  call blank     (hdr,hdrsize)
                  call byte_read (hdr,hdrsize/4,ierr)
               endif
               if (ierr == 0) then
                  call byte_read (test_pattern,1,ierr)
                  if_byte_sw = if_byte_swap_test(test_pattern,ierr) ! determine endianess
               endif
               call nek_log_debug('header: '//trim(hdr), this_module, 'get_nsteps_from_header')
               ! read wdsize from header
               read(hdr,*) sdummy, wdsizr
               wdsl = wdsizr/4
               isl  = isize/4
               ! read metadata
               call byte_read(nxr,    isl, ierr)
               call byte_read(nyr,    isl, ierr)
               call byte_read(nelfr,  isl, ierr)
               call byte_read(timer, wdsl, ierr)
               call byte_read(nsaver, isl, ierr)
               call byte_read(lbufr,  isl, ierr)
               write(msg,'(A,3(1X,I0),1X,E15.7,2(1X,I0))') 'metadata: ', nxr, nyr, nelfr, timer, nsaver, lbufr
               call nek_log_debug(msg, this_module, 'get_nsteps_from_header')
               call byte_close(ierr)
               if (ierr /= 0) call nek_stop_error('Error closing file '//trim(fname), procedure='get_nsteps_from_header')
            end if
            write(msg,'(3X,A,A,I0,A)') trim(fname), ': ', nsaver, ' timesteps.'
            call nek_log_message(msg, this_module, 'get_nsteps_from_header')
         end procedure get_nsteps_from_header

         ! Helper functions

         subroutine gather_and_write_slice(slicedata, n2d_gown)
            real(dp), intent(in) :: slicedata(:,:,:)
            integer, intent(in) :: n2d_gown(:)
            ! internal
            integer :: nxy, idum, wdsl, isl, length, ierr, ip
            real rtmpv1(lx1*ly1*lelv), rtmpv(lx1*ly1*lelv)
            real*4 rtmpv2(2*lx1*ly1*lelv)
            equivalence (rtmpv1,rtmpv2)
            nxy  = lx1*ly1
            wdsl = wdsize/4
            isl  = isize/4
            if (nid == 0) then
               ! master writes if there are data
               length = nxy*n2d_gown(nid+1)
               if (wdsl.eq.2) then
                  call copy(rtmpv1,slicedata,length)
                  call byte_write(rtmpv2,length*wdsl,ierr)
               else
                  call copyX4(rtmpv2,slicedata,length)
                  call byte_write(rtmpv2,length,ierr)
               end if
               ! get data from other procs and write to file
               do ip = 1, np-1
                  length = nxy*n2d_gown(ip+1)
                  call csend(ip,idum,isize,ip,0) ! hand shake
                  call crecv2(ip,rtmpv,length*wdsize,ip)
                  ! write data
                  if (wdsl.eq.2) then
                     call copy(rtmpv1,rtmpv,length)
                     call byte_write(rtmpv2,length*wdsl,ierr)
                  else
                     call copyX4(rtmpv2,rtmpv,length)
                     call byte_write(rtmpv2,length,ierr)
                  endif
               end do
               if (ierr /= 0) call nek_stop_error('Error writing slice data', procedure='gather_and_write_slice')
            else 
               ! send data to master
               call crecv2(nid,idum,isize,0) ! hand shake
               length = nxy*n2d_gown(nid+1)
               call csend(nid,slicedata,length*wdsize,0,0)
            end if
         end subroutine gather_and_write_slice

      end submodule helix_io