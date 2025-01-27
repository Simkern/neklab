      submodule (neklab_helix) helix_2d
         implicit none
      
      contains

         module procedure init_2d_geom
            integer :: ie, iel, ieg, iseg, iface, isl, nxy
            integer, dimension(lelv) :: islice
            integer, dimension(:), allocatable :: unique_segments, segment_owner, segment_count
            integer, dimension(:), allocatable :: idx ! for findloc
            logical, dimension(:), allocatable :: segment_found
            character(len=128) :: msg
            ! for debug
            logical :: debug
            character(len=3) :: fid
            integer :: fileid
            ! functions
            integer, external :: iglsum

            debug = optval(if_debug, .false.)
            
            nxy = lx1*ly1
            call nek_log_message('start extraction', this_module, 'init_2d_geom')

            ! Sort local elements according to 2D mesh. 
            ! Here we use a trick that relies on the particular structure of meshes extruded
            ! from a 2D mesh using n2to3.
            allocate(unique_segments(self%nelf)); call izero(unique_segments, self%nelf)
            allocate(segment_count  (self%nelf)); call izero(segment_count,   self%nelf)
            allocate(segment_owner  (self%nelf)); call izero(segment_owner,   self%nelf)
            allocate(segment_found  (self%nelf), source=.false.)
            call izero(islice, nelv)
            call izero(self%id2d, 3*nelv)
            self%n2d_lown = 0
            self%n2d_gown = 0
            do ie = 1, nelv
               ieg = lglel(ie)
               isl = ieg/self%nelf + 1
               if (mod(ieg,self%nelf)==0) isl = isl - 1
               islice(ie) = isl
               self%gsegment(ie) = ieg - (isl-1)*self%nelf
               if (isl == 1) then
                  ! the element in the first slice is the global segment owner
                  self%gowner(ie) = .true.
                  self%n2d_gown = self%n2d_gown + 1
               end if
               ! gather all unique segments on proc
               if (.not. segment_found(self%gsegment(ie))) then
                  self%n2d_lown = self%n2d_lown + 1
                  segment_owner(self%n2d_lown) = ie
                  unique_segments(self%n2d_lown) = self%gsegment(ie)
                  segment_found(self%gsegment(ie)) = .true.
                  self%id2d(self%n2d_lown,1) = ie
                  self%id2d(self%n2d_lown,3) = self%gsegment(ie)
                  self%lowner(ie) = .true.
               end if
               idx = findloc(unique_segments, self%gsegment(ie))
               segment_count(idx(1)) = segment_count(idx(1)) + 1
               self%lsegment(ie) = idx(1)
            end do
      
            if (debug) then
               print '(A,2(I0,1X),A,*(1X,I3))', 'DEBUG 2dmap: ', nid, self%n2d_lown, 'unique streamwise segments:   ', unique_segments(:self%n2d_lown)
               print '(A,2(I0,1X),A,*(1X,I3))', 'DEBUG 2dmap: ', nid, self%n2d_lown, '# of elements in segment:     ', segment_count(:self%n2d_lown)
               print '(A,2(I0,1X),A,*(1X,I3))', 'DEBUG 2dmap: ', nid, self%n2d_lown, 'local element local  s. owner:', segment_owner(:self%n2d_lown)
               print '(A,2(I0,1X),A,*(1X,L3))', 'DEBUG 2dmap: ', nid, self%n2d_lown, 'local element global s. owner:', self%gowner(:self%n2d_lown)
               print '(A,2(I0,1X),A,I0)'      , 'DEBUG 2dmap: ', nid, self%n2d_lown, 'globally owned: ', self%n2d_gown
               call nekgsync()
               do ie = 1, nelv
                  call cfill(vz(1,1,1,ie), 1.0_dp*nid, lx1*ly1*lz1)
                  call cfill(vx(1,1,1,ie), 1.0_dp*islice(ie),lx1*ly1*lz1)
                  call cfill(vy(1,1,1,ie), 1.0_dp*self%gsegment(ie),lx1*ly1*lz1)
                  if (self%lowner(ie)) call cfill(vz(1,1,1,ie), -1.0_dp,lx1*ly1*lz1) 
                  if (self%gowner(ie)) call cfill(vz(1,1,1,ie), np*1.0_dp,lx1*ly1*lz1)
               end do
               call outpost(vx,vy,vz,pr,t,'m2d')
               print '(A,A3,3X,3(A5),2X,2(A20))', 'DEBUG: 2dmap ', 'nid', 'glob', 'locl', 'l2d', 'xavg', 'yavg'
            end if

            ! Extract 2D mesh
            iseg = 0 
		      do ie = 1, nelv
               if (self%lowner(ie)) then
		      	   ! extract boundary points from the global segment owners
                  iseg = iseg + 1
                  if (self%gowner(ie)) then
		      	      do iface = 1, 2*ndim
		      	      	if (cbc(iface,ie,1) == 'P') then
		      	   	   	call ftovec(self%x2d(1,1,iseg), zm1, ie, iface, nx1, ny1, nz1) ! z --> x
		      	   		   call ftovec(self%y2d(1,1,iseg), ym1, ie, iface, nx1, ny1, nz1)
                           self%id2d(iseg,2) = iface
                       end if
		      	      end do
                  end if
                  if (debug) then
                     print '(A,I3,A,4(1X,I4),A,3X,F17.8,3x,F17.8)', 'DEBUG 2dmap: ', nid, ' el', lglel(ie), ie, iseg, self%gsegment(ie),  
     &                           ': ', sum(self%x2d(:,:,iseg))/nxy, sum(self%y2d(:,:,iseg))/nxy
                  end if
               end if
		      end do
            
            ! initialize data
            call rzero(self%vx2d, nx1*ny1*nelv*lbuf)
            call rzero(self%vy2d, nx1*ny1*nelv*lbuf)
            call rzero(self%vz2d, nx1*ny1*nelv*lbuf)
            call rzero(self%dt2d, lbuf)
            self%nsave = 0
            self%noutc = 0
            self%noutt = 0
            self%n2d   = iglsum(self%n2d_gown,1)
            self%nload = 0
            if (self%n2d /= self%nelf) then
               call nek_stop_error('Inconsistent elements in 2D mesh!', module=this_module, procedure='init_2d_geom')
            else
               call nek_log_message('global 2D element ownership establised', this_module, 'init_2d_geom')
               msg = 'neklab_helix % init_2d_geom :'
               do ie = 0, np-1
                  if (nid == ie) print '(A,4X,A,I3,A,I3,A)', trim(msg), 'proc ', ie, ': ', self%n2d_gown, ' 2D elements'
                  call nekgsync()
               end do
               write(msg,'(A,I3,A)') 'total: ', self%n2d, ' 2D elements'
               call nek_log_message(msg, this_module, 'init_2d_geom')
            end if
            if (debug) then
               call nekgsync()
               print '(A,I3,A,*(1x,I0))', 'DEBUG 2dmap: ', nid, ', nelv2iseg: ', self%lsegment(:nelv)
               write(fid,'(I3.3)') nid
               fileid = 2000+nid
               open (fileid, file='torus_map'//fid//'.txt', status='replace', action='write')
               write(fileid, *) 'nelv = ', nelv
               write(fileid, '(6(1X,A11),2(1X,A7),A12)') 'ie','ieg','slice','self%gsegment','gllel','gllnid','s%lowner','s%gowner','s%n2iseg'
               do ie = 1, nelv
                  ieg = lglel(ie)
                  write(fileid, '(6(I12),2(1X,L7),I12)') ie, ieg, islice(ie), self%gsegment(ie), gllel(ieg), gllnid(ieg), self%lowner(ie), self%gowner(ie), self%lsegment(ie)
               end do
               write(fileid, *) 'n2d_lown = ', self%n2d_lown
               write(fileid, '(*(1X,A11))') 'iel','ieg','s%id2d:ie', 's%id2d:ifc', 's%id2d%iseg'
               do ie = 1, self%n2d_lown
                  write(fileid, *) ie, lglel(self%id2d(ie,1)), self%id2d(ie,:)
               end do
               close (fileid)
               call nekgsync()
            end if
            call nek_log_message('extraction complete', this_module, 'init_2d_geom')
         end procedure init_2d_geom

         module procedure save_2d_fields
            integer :: iseg, ie, ifc, level, nxy
            real(dp) :: xavg, yavg, vxavg, vyavg, vzavg
            character(len=128) :: msg
            nxy = lx1*ly1
            call logger%configuration(level=level)
            if (self%save_2d_base) then
               call lk_timer%start('neklab_helix_save_2d')
               self%nsave = self%nsave + 1
               ! save data to buffer
               write(msg,'(A,I5,A,I5,A,E12.5,A,F12.8)') 'Save 2D field ', self%nsave, '/', lbuf, ', time=', time, ', dt=', dt
               call logger%log_debug(msg, this_module, 'save_2d_fields')
               if (nid == 0) print '(A,A)', 'neklab_helix: ', trim(msg)
               do iseg = 1, self%n2d_gown
                  ie  = self%id2d(iseg, 1)
                  ifc = self%id2d(iseg, 2)
                  call ftovec(self%vx2d(1,1,iseg,self%nsave), u, ie, ifc, nx1, ny1, nz1)
		         	call ftovec(self%vy2d(1,1,iseg,self%nsave), v, ie, ifc, nx1, ny1, nz1)
		         	call ftovec(self%vz2d(1,1,iseg,self%nsave), w, ie, ifc, nx1, ny1, nz1)
                  if (level == all_level) then
                     xavg  = sum(self%x2d (:,:,iseg))/nxy
                     yavg  = sum(self%y2d (:,:,iseg))/nxy
                     vxavg = sum(self%vx2d(:,:,iseg,self%nsave))/nxy
                     vyavg = sum(self%vy2d(:,:,iseg,self%nsave))/nxy
                     vzavg = sum(self%vz2d(:,:,iseg,self%nsave))/nxy
                     print '(A,I8,I8,A,5(3X,F16.8))', 'DEBUG: save el', ie, iseg, ': ', xavg, yavg, vxavg, vyavg, vzavg
                  end if
                  call lk_timer%stop('neklab_helix_save_2d')
               end do
               ! save timestep information and record minimum dt
               self%dt2d(self%nsave) = dt
               if (lastep == 0) then ! exclude the potentially very short last step
                  self%min_dt = min(dt, self%min_dt)
                  self%max_dt = max(dt, self%max_dt)
               end if
               ! save data to file when buffer is full
               if (self%nsave == lbuf .or. lastep == 1) call self%outpost_2d()
            else
               call nek_log_debug('Baseflow saving turned off', this_module, 'save_2d_fields')
            end if
         end procedure save_2d_fields
         
         module procedure outpost_2d
            if (self%nsave > 0) then
               call lk_timer%start('neklab_helix_outpost_2d')
               if (self%is_newton()) then
                  self%noutn = self%noutn + 1
                  call self%outpost_2d_fields(iname='n', iout=self%noutn)
               else if (self%is_floquet()) then
                  self%noutn = self%noutn + 1
                  call self%outpost_2d_fields(iname='f', iout=self%noutn)
               else
                  self%noutc = self%noutc + 1
                  call self%outpost_2d_fields(iname='c', iout=self%noutc)
                  if (self%save_2d_usrt) then
                     call self%compute_2d_usrt() ! self%v[xyz]2d are overwritten
                     self%noutt = self%noutt + 1
                     call self%outpost_2d_fields(iname='t', iout=self%noutt)
                  end if
               end if
               self%nsave = 0
               call lk_timer%stop('neklab_helix_outpost_2d')
            else
               call nek_log_message('No 2D data to outpost.', this_module, 'outpost')
            end if
         end procedure outpost_2d
         
         module procedure outpost_2d_fields
            integer, allocatable :: n2d_gown(:)
            integer, allocatable :: n2d_elmap(:)
            integer :: ierr, itmp, i, nxy, ip, ibuf, iseg, length, i_own, nsave
            integer :: wdsl, isl, isend(lelv)
            character(len=128)  :: fname, msg
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
            write(fname,'(A,A,I3.3,A)') iname, '2dtorus', iout, '.fld'
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
         end procedure outpost_2d_fields

         module procedure load_2d_fields
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
            character(len=132) :: hdr, fname, msg
            character(len=4)   :: sdummy
            common /CTMP1/ fldum(lx1*ly1*lelv)
            real fldum
            ! functions
            logical, external :: if_byte_swap_test
            call lk_timer%start('neklab_helix_load_2d')
            if (self%is_newton()) then
               write(fname,'("n2dtorus",I3.3,".fld")') idx
            else if (self%is_floquet()) then
               write(fname,'("f2dtorus",I3.3,".fld")') idx
            else
               write(fname,'("c2dtorus",I3.3,".fld")') idx
            end if
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
         end procedure load_2d_fields

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

         module procedure set_baseflow
            integer  :: ie, ix, iy, iz, iseg, ifld_
            real(dp) :: s, phi, u, v, w
            character(len=128) :: msg
            phi = self%phi
            if (self%is_newton() .or. (self%is_floquet() .and. .not. self%is_save_2d()) ) then
               ifld_ = ifld - (self%noutn-1)*lbuf
               if (ifld_ > self%nload) then
                  ! load next file
                  self%noutn = self%noutn + 1
                  write(msg,'(A,I5)') 'Load file: ', self%noutn
                  call nek_log_debug(msg, this_module, 'set_baseflow')
                  call self%load_2d_fields(self%noutn)
               end if
               ifld_ = ifld - (self%noutn-1)*lbuf
            else
               ifld_ = ifld
               if (ifld_ > self%nload) call nek_stop_error('Inconsistent ifld!', this_module, 'set_baseflow')
            end if
            call lk_timer%start('neklab_helix_set_baseflow')
            ! set dt
            param(12) = -abs(self%dt2d(ifld_)) ! negative to force the stepsize in settime
            write(msg,'(A,I5,"/",I5,A,I5,A,F10.6)') 'Set field ', ifld_, lbuf, ' (', ifld, '), dt= ', -param(12)
            call logger%log_debug(msg, this_module, 'set_baseflow')
            if (nid == 0) print '(A,A)', 'neklab_helix: ', trim(msg)
            do ie = 1, nelv
            iseg = self%lsegment(ie) ! local segment
            do iz = 1, lz1
            do iy = 1, ly1
            do ix = 1, lx1
               s = self%as(ix,iy,iz,ie)
               u = self%vx2d(ix,iy,iseg,ifld_)
               v = self%vy2d(ix,iy,iseg,ifld_)
               w = self%vz2d(ix,iy,iseg,ifld_)
               basex(ix,iy,iz,ie) = cos(phi)*( cos(s)*u + sin(s)*v) + sin(phi)*w
               basey(ix,iy,iz,ie) =           -sin(s)*u + cos(s)*v
               basez(ix,iy,iz,ie) = sin(phi)*(-cos(s)*u - sin(s)*v) + cos(phi)*w
            end do
            end do
            end do
            end do
            call lk_timer%stop('neklab_helix_set_baseflow')
         end procedure set_baseflow

         module procedure compute_2d_usrt
            ! this routine will overwrite self%v[xyz]2d
            integer, parameter :: iz = 1
            integer :: ix, iy, ie, is, ib
            real(dp) :: phi, a, s
            real(dp) :: utmp, vtmp, ux, uy, uz
            phi = self%phi
            do is = 1, self%n2d_gown ! only for the first slice
               ie = self%id2d(is, 1)
               do iy = 1, ly1
               do ix = 1, lx1
                  s = self%as(ix,iy,iz,ie)
                  a = self%alpha(ix,iy,iz,ie)
                  ! iterate over buffer
                  do ib = 1, lbuf
                     ux = self%vx2d(ix,iy,is,ib)
                     uy = self%vy2d(ix,iy,is,ib)
                     uz = self%vz2d(ix,iy,is,ib)
                     ! overwrite v[xyz]2d with u[srt]2d
                     self%vx2d(ix,iy,is,ib) = cos(phi)*( cos(s)*ux -sin(s)*uy) + sin(phi)*uz
                     utmp                   = sin(s)*ux + cos(s)*uy
                     vtmp                   = sin(phi)*(-cos(s)*ux -sin(s)*uy) + cos(phi)*uz
                     self%vy2d(ix,iy,is,ib) = cos(a)*utmp + sin(a)*vtmp
                     self%vz2d(ix,iy,is,ib) = sin(a)*utmp - cos(a)*vtmp
                  end do ! lbuf
               end do    ! lx1
               end do    ! ly1
            end do       ! self%n2d_gown
         end procedure compute_2d_usrt
         
         module procedure set_2d_mode
            if (trim(mode)=='newton') then
               self%noutn = 0
               self%nload = 0
               self%min_dt = 100.0_dp
               self%max_dt = 0.0_dp
               call self%reset_mflow_fft()
               call self%set_newton(.true.)
               call self%set_floquet(.false.)
               call self%set_save_base(.true.)
            else if (trim(mode)=='floquet') then
               self%noutn = 0
               self%nload = 0
               call self%set_floquet(.true.)
               call self%set_newton(.false.)
            else
               call nek_stop_error('Selected mode '//trim(mode)//' is invalid.', this_module, 'set_2d_mode')
            end if
         end procedure set_2d_mode

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
      
      end submodule helix_2d