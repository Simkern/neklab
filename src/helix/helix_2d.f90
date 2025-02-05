      submodule (neklab_helix) helix_2d
         implicit none
      
      contains

         module procedure init_2d_geom
            integer :: ie, iel, ieg, iseg, iface, isl, nxy, nelf, nslices, nown
            integer, dimension(lelv) :: islice
            integer, dimension(:), allocatable :: unique_segments, segment_owner, segment_count
            integer, dimension(:), allocatable :: idx ! for findloc
            logical, dimension(:), allocatable :: segment_found
            real(dp) :: xmin
            logical :: has_P
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

            ! Find number of elements in the first slice of the mesh
            nelf = 0
            do ie = 1, nelv
               ! the element should have a face on the yz plane
               xmin = abs(minval(xm1(:,:,:,ie)))
               has_P = .false.
               do iface = 1, 2*ndim
               ! and have a periodic bc
                  if (cbc(iface,ie,1) == 'P  ') has_P = .true.
               end do
               if (xmin < 1e-6_dp .and. has_P) nelf = nelf + 1
            end do
            ! gather info from all procs
            nelf = iglsum(nelf, 1)
            ! deduce number of slices
            nslices = nelgv/nelf
            ! stamp logs
            write(msg,'(A,I8)') 'Elements in cross-stream plane:   ', nelf
            call nek_log_message(msg, this_module, 'init_2d_geom')
            write(msg,'(A,I8)') 'Elements in streamwise direction: ', nslices
            call nek_log_message(msg, this_module, 'init_2d_geom')
            write(msg,'(A,I8)') 'Total number of elements mesh:    ', nelgv
            call nek_log_message(msg, this_module, 'init_2d_geom')

            ! sanity check
            if (nelf*nslices /= nelgv) then
               call nek_stop_error('Inconsistent mesh partitioning!', module=this_module, procedure='init_2d_geom')
            else
               self%nelf = nelf
               self%nslices = nslices
            end if

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
               write(msg,'(A,I0)') 'Number of globally owned 2D elements: ', self%n2d
               call nek_log_message(msg, this_module, 'init_2d_geom')
               write(msg,'(A,I0)') 'Total number of 2D elements:          ', self%nelf
               call nek_log_message(msg, this_module, 'init_2d_geom')
               call nek_stop_error('Inconsistent element ownership!', module=this_module, procedure='init_2d_geom')
            else
               call nek_log_message('global 2D element ownership established', this_module, 'init_2d_geom')
               ! this is quite ugly, but it's just once for information purposes ...
               do ie = 0, np-1
                  nown = 0
                  if (nid == ie) nown = self%n2d_gown
                  nown = iglsum(nown, 1)
                  write(msg,'(A,I3,A,I3,A)') 'proc ', ie, ': ', nown, ' 2D elements'
                  call nek_log_message(msg, this_module, 'init_2d_geom')
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
               if (self%nsave == lbuf .or. lastep == 1) call self%outpost_2d_fields()
            else
               call nek_log_debug('Baseflow saving turned off', this_module, 'save_2d_fields')
            end if
         end procedure save_2d_fields
         
         module procedure outpost_2d_fields
            if (self%nsave > 0) then
               call lk_timer%start('neklab_helix_outpost_2d')
               if (self%is_newton()) then
                  self%noutn = self%noutn + 1
                  call self%write_2d(self%fname_2d('n', self%noutn))
               else if (self%is_floquet()) then
                  self%noutn = self%noutn + 1
                  call self%write_2d(self%fname_2d('f', self%noutn))
               else
                  self%noutc = self%noutc + 1
                  call self%write_2d(self%fname_2d('c', self%noutc))
                  if (self%save_2d_usrt) then
                     call self%compute_2d_usrt() ! self%v[xyz]2d are overwritten
                     self%noutt = self%noutt + 1
                     call self%write_2d(self%fname_2d('t', self%noutt))
                  end if
               end if
               self%nsave = 0
               call lk_timer%stop('neklab_helix_outpost_2d')
            else
               call nek_log_message('No 2D data to outpost.', this_module, 'outpost')
            end if
         end procedure outpost_2d_fields

         module procedure load_2d_fields
            ! only nid 0 will read
            character(len=132) :: fname
            call lk_timer%start('neklab_helix_load_2d')
            if (self%is_newton()) then
               fname = self%fname_2d('n', idx)
            else if (self%is_floquet()) then
               fname = self%fname_2d('f', idx)
            else
               fname = self%fname_2d('c', idx)
            end if
            call self%read_2d(fname)
            call lk_timer%stop('neklab_helix_load_2d')
         end procedure load_2d_fields

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

         module procedure load_baseflow
            ! internal
            integer :: ifld_, nchar
            character(len=132) :: filename
            ifld_ = optval(ifld, 1)
            call blank(filename, 132)
            nchar = min(len(fname), 132)
            filename(1:nchar) = fname(1:nchar)
            call self%read_2d(filename)
            if (ifld > self%nload) call nek_stop_error('Inconsistent ifld', this_module, 'load_baseflow')
            call self%set_baseflow(vx, vy, vz, ifld)
         end procedure load_baseflow

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
      
      end submodule helix_2d