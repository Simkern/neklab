      submodule(t2Dh_bf_buffer) t2Dh_bf_buffer_io
      !! Binary IO for the baseflow snapshot buffer.
      !!
      !! Wire format is byte-identical to the helix 2D format, so
      !! read_2d_data.py / write_2d_data.py handle these files unchanged:
      !!
      !!   [116 B]  ascii header  '#1tf <wdsize> (lx1, ly1 = ..) (nelf = ..)
      !!                           (time = ..) (nsave, lbuf = ..)'
      !!   [  4 B]  real*4 endian test pattern (6.54321)
      !!   [ int ]  lx1, ly1, nelf
      !!   [ flt ]  time
      !!   [ int ]  nsave, lbuf
      !!   [ int ]  elmap(nelf)        global element number of each written el.
      !!   [ flt ]  dt(nsave)
      !!   [ flt ]  x(lx1,ly1,nelf)
      !!   [ flt ]  y(lx1,ly1,nelf)
      !!   repeated nsave times:
      !!   [ flt ]  vx(lx1,ly1,nelf)   u_Z
      !!   [ flt ]  vy(lx1,ly1,nelf)   u_R
      !!   [ flt ]  vz(lx1,ly1,nelf)   u_phi   (the streamwise component)
      !!
      !! Elements are written in rank order, and elmap records the global number
      !! of each one, so a reader on a different number of ranks recovers the
      !! correct ordering by sorting elmap. Data are written at full wdsize
      !! precision: the Jacobian consumes this baseflow, and truncating to
      !! real*4 would put an O(1e-7) inconsistency between F and dF.
      !!
      !! Simplification versus helix_io: there are no streamwise segments here,
      !! so every local element is a 2D element and every rank contributes
      !! nelv of them. No ownership bookkeeping is needed.
      !!
      !! LOGGING: everything in here is debug. The one line per file worth
      !! seeing by default is emitted by flush_chunk in the parent module.
         implicit none

      contains

      !====================================================================
      !     WRITE
      !====================================================================

         module procedure bf_write_chunk
            character(len=*), parameter :: this_procedure = 'bf_write_chunk'
            integer, allocatable :: ngown(:), elmap(:)
            integer :: ierr, itmp, i, ip, iseg, i_own, nsave_
            integer :: wdsl, isl
            integer :: isend(lelv)
            character(len=256) :: msg
            character(len=1024) :: head, ftm
            character(len=2) :: id
            character(len=1), parameter :: fileversion = '1'
            logical :: only_mesh_
            real*4 :: test
            parameter(test=6.54321)

            call neklab_timer_start(t_bf_write)
            wdsl = wdsize/4
            isl = isize/4
            only_mesh_ = optval(only_mesh, .false.)

            if (only_mesh_) then
               nsave_ = 0
               write (msg, '(A,A)') 'Write 2D mesh: ', trim(fname)
            else
               nsave_ = nsave
               write (msg, '(A,I0,A,A)') 'Write 2D data: ', nsave_, ' snapshots -> ', trim(fname)
            end if
            call nek_log_debug(msg, this_module, this_procedure)

            ierr = 0
            if (nid == 0) call byte_open(fname, ierr)
            call bcast(ierr, isize)
            if (ierr /= 0) call nek_stop_error('Error opening '//trim(fname), this_module, this_procedure)

      ! --- header
            if (nid == 0) then
               ftm = "('#',A1,A2,1x,i1,1x,'(lx1, ly1 =',2i9,') (nelf =',i9,"//
     &               "') (time =',e17.9,') (nsave, lbuf =', 2i9,')')"
               id = merge('th', 'tf', if_sym)
               write (head, ftm) fileversion, id, wdsize, lx1, ly1, nelf, time, nsave_, lbuf
               call byte_write(head, 116/4, ierr)
            end if
            call bcast(ierr, isize)
            if (ierr /= 0) call nek_stop_error('Error writing header in '//trim(fname), this_module, this_procedure)

      ! --- endian test pattern and metadata
            if (nid == 0) then
               call byte_write(test, 1, ierr)
               call byte_write(lx1, isl, ierr)
               call byte_write(ly1, isl, ierr)
               call byte_write(nelf, isl, ierr)
               call byte_write(time, wdsl, ierr)
               call byte_write(nsave_, isl, ierr)
               call byte_write(lbuf, isl, ierr)
            end if
            call bcast(ierr, isize)
            if (ierr /= 0) call nek_stop_error('Error writing metadata in '//trim(fname), this_module, this_procedure)

      ! --- element counts and the global element map
            allocate (ngown(np)); call izero(ngown, np)
            allocate (elmap(nelf)); call izero(elmap, nelf)

            iseg = 0
            if (nid == 0) then
               ngown(1) = nelv
               do i = 1, nelv
                  elmap(i) = lglel(i)
               end do
               iseg = nelv
               do ip = 1, np - 1
                  call csend(ip, itmp, isize, ip, 0)            ! hand shake
                  call crecv(ip, i_own, isize)                  ! element count
                  ngown(ip + 1) = i_own
                  call crecv(ip, isend(:i_own), i_own*isize)    ! global numbers
                  elmap(iseg + 1:iseg + i_own) = isend(:i_own)
                  iseg = iseg + i_own
               end do
               call byte_write(elmap, nelf*isl, ierr)
               if (nsave_ > 0) call byte_write(bufdt(:nsave_), nsave_*wdsl, ierr)
            else
               call crecv(nid, itmp, isize)                     ! hand shake
               call csend(nid, nelv, isize, 0, 0)               ! element count
               do i = 1, nelv
                  isend(i) = lglel(i)
               end do
               call csend(nid, isend(:nelv), nelv*isize, 0, 0)  ! global numbers
            end if
            call bcast(iseg, isize)
            if (iseg /= nelf) then
               write (msg, '(A,I0,A,I0)') 'Element count mismatch: gathered ', iseg, ' of ', nelf
               call nek_stop_error(msg, this_module, this_procedure)
            end if
            call bcast(ngown, np*isize)

      ! --- coordinates
            call nek_log_debug('   '//trim(fname)//': write x, y ...', this_module, this_procedure)
            call gather_and_write(x2d, ngown)
            call gather_and_write(y2d, ngown)

      ! --- fields, one snapshot at a time
            if (.not. only_mesh_) then
               write (msg, '(3X,A,A,I0)') trim(fname), ': write u_Z, u_R, u_phi for snapshots: ', nsave_
               call nek_log_debug(msg, this_module, this_procedure)
               do i = 1, nsave_
                  call gather_and_write(bufx(:, i), ngown)
                  call gather_and_write(bufy(:, i), ngown)
                  call gather_and_write(bufs(:, i), ngown)
               end do
            end if

            if (nid == 0) call byte_close(ierr)
            call bcast(ierr, isize)
            if (ierr /= 0) call nek_stop_error('Error closing '//trim(fname), this_module, this_procedure)

            deallocate (ngown, elmap)
            call neklab_timer_stop(t_bf_write)
         end procedure bf_write_chunk

      !====================================================================
      !     READ
      !====================================================================

         module procedure bf_read_chunk
            character(len=*), parameter :: this_procedure = 'bf_read_chunk'
            integer :: ierr, hdrsize
            real*4 :: test_pattern
            integer :: nxr, nyr, nelfr, nsaver, lbufr, wdsizr
            integer :: wdsl, isl, nxy, i, ie, ieg, iel, length
            integer, allocatable :: elmap(:), gmap_index(:)
            real(dp), allocatable :: fdata(:, :, :, :)
            real(dp), allocatable :: xr(:, :, :), yr(:, :, :)
            real(dp) :: timer(1), xyavg, xyavgr
            real(dp), parameter :: mesh_tol = 1.0e-10_dp
            character(len=132) :: hdr
            character(len=256) :: msg
            character(len=4) :: sdummy
            logical, external :: if_byte_swap_test
            integer, external :: iglsum

            call neklab_timer_start(t_bf_read)
            hdrsize = 116
            nxy = lx1*ly1
            isl = isize/4
            nread = 0
            call nekgsync()   ! avoid false positives in the error checks below

            allocate (elmap(nelf)); call izero(elmap, nelf)
            allocate (gmap_index(nelf)); call izero(gmap_index, nelf)

            ierr = 0
            wdsizr = wdsize
            if (nid == 0) call byte_open(fname, ierr)
            call bcast(ierr, isize)
            if (ierr /= 0) call nek_stop_error('Error opening '//trim(fname), this_module, this_procedure)

      ! --- header and metadata
            if (nid == 0) then
               call blank(hdr, hdrsize)
               call byte_read(hdr, hdrsize/4, ierr)
               if (ierr == 0) then
                  call byte_read(test_pattern, 1, ierr)
                  if_byte_sw = if_byte_swap_test(test_pattern, ierr)
               end if
            end if
            call bcast(ierr, isize)
            if (ierr /= 0) call nek_stop_error('Error reading header from '//trim(fname), this_module, this_procedure)

            if (nid == 0) then
               call nek_log_debug('header: '//trim(hdr), this_module, this_procedure)
      ! wdsl comes from the FILE, not from the run: reading the following
      ! real with the wrong width would misalign everything after it and the
      ! resulting error would point at the wrong thing.
               read (hdr, *) sdummy, wdsizr
               wdsl = wdsizr/4
               call byte_read(nxr, isl, ierr)
               call byte_read(nyr, isl, ierr)
               call byte_read(nelfr, isl, ierr)
               call byte_read(timer, wdsl, ierr)
               call byte_read(nsaver, isl, ierr)
               call byte_read(lbufr, isl, ierr)
               if (if_byte_sw) then
                  call byte_reverse(nxr, 1, ierr)
                  call byte_reverse(nyr, 1, ierr)
                  call byte_reverse(nelfr, 1, ierr)
                  call byte_reverse(nsaver, 1, ierr)
                  call byte_reverse(lbufr, 1, ierr)
                  call reverse_real(timer(1), 1, wdsizr, ierr)
               end if
            end if
            call bcast(ierr, isize)
            if (ierr /= 0) call nek_stop_error('Error reading metadata from '//trim(fname), this_module, this_procedure)
            call bcast(nxr, isize); call bcast(nyr, isize); call bcast(nelfr, isize)
            call bcast(nsaver, isize); call bcast(lbufr, isize); call bcast(wdsizr, isize)

      ! --- consistency of the discretisation
            if (wdsizr /= wdsize) then
               write (msg, '(A,I0,A,I0)') 'Word size mismatch: file ', wdsizr, ', run ', wdsize
               call nek_stop_error(msg, this_module, this_procedure)
            end if
            wdsl = wdsize/4
            if (nxr /= lx1 .or. nyr /= ly1) then
               write (msg, '(A,2(1X,I0),A,2(1X,I0))') 'Polynomial order mismatch: file', nxr, nyr, ', run', lx1, ly1
               call nek_stop_error(msg, this_module, this_procedure)
            end if
            if (nelfr /= nelf) then
               write (msg, '(A,I0,A,I0)') 'Element count mismatch: file ', nelfr, ', run ', nelf
               call nek_stop_error(msg, this_module, this_procedure)
            end if
            if (nsaver > lbuf) then
               write (msg, '(A,I0,A,I0)') 'File holds ', nsaver, ' snapshots but lbuf = ', lbuf
               call nek_stop_error(msg, this_module, this_procedure)
            end if

      ! --- element map and timesteps
            call rzero(bufdt, lbuf)   ! do not let a previous chunk's tail survive
            if (nid == 0) then
               call byte_read(elmap, nelf*isl, ierr)
               if (if_byte_sw) call byte_reverse(elmap, nelf, ierr)
               if (nsaver > 0) then
                  call byte_read(bufdt, nsaver*wdsl, ierr)
                  if (if_byte_sw) call reverse_real(bufdt, nsaver, wdsize, ierr)
               end if
            end if
            call bcast(ierr, isize)
            if (ierr /= 0) call nek_stop_error('Error reading element map from '//trim(fname), this_module, this_procedure)
            call bcast(elmap, nelf*isize)
            call bcast(bufdt, lbuf*wdsize)

      ! Position of global element ieg within the file is gmap_index(ieg).
            call sort_index(elmap, gmap_index)

      ! --- coordinates: read and verify against the running mesh
            length = nxy*nelf
            allocate (xr(lx1, ly1, nelf)); allocate (yr(lx1, ly1, nelf))
            if (nid == 0) then
               call byte_read(xr, length*wdsl, ierr)
               if (if_byte_sw) call reverse_real(xr, length, wdsize, ierr)
               call byte_read(yr, length*wdsl, ierr)
               if (if_byte_sw) call reverse_real(yr, length, wdsize, ierr)
            end if
            call bcast(ierr, isize)
            if (ierr /= 0) call nek_stop_error('Error reading coordinates from '//trim(fname), this_module, this_procedure)
            call bcast(xr, length*wdsize)
            call bcast(yr, length*wdsize)

            ierr = 0
            do ie = 1, nelv
               ieg = lglel(ie)
               iel = gmap_index(ieg)
               xyavg = sum(xm1(:, :, 1, ie)) + sum(ym1(:, :, 1, ie))
               xyavgr = sum(xr(:, :, iel)) + sum(yr(:, :, iel))
      ! Relative: xyavg is a SUM over lx1*ly1 coordinates, so on a torus at
      ! R0 = O(10) an absolute 1e-10 sits at the double-precision floor.
               if (abs(xyavg - xyavgr) > mesh_tol*max(1.0_dp, abs(xyavg))) ierr = ierr + 1
            end do
            ierr = iglsum(ierr, 1)
            if (ierr > 0) then
               write (msg, '(A,I0,A,A)') 'Mesh mismatch on ', ierr, ' element(s) in ', trim(fname)
               call nek_stop_error(msg, this_module, this_procedure)
            end if

      ! --- fields, one snapshot at a time
            length = 3*nxy*nelf
            allocate (fdata(lx1, ly1, nelf, 3)); call rzero(fdata, length)
            do i = 1, nsaver
               if (nid == 0) then
                  call byte_read(fdata, length*wdsl, ierr)
                  if (if_byte_sw) call reverse_real(fdata, length, wdsize, ierr)
               end if
               call bcast(ierr, isize)
               if (ierr /= 0) call nek_stop_error('Error reading snapshot data', this_module, this_procedure)
               call bcast(fdata, length*wdsize)
               do ie = 1, nelv
                  ieg = lglel(ie)
                  iel = gmap_index(ieg)
                  call copy(bufx(1 + (ie - 1)*nxy, i), fdata(1, 1, iel, 1), nxy)
                  call copy(bufy(1 + (ie - 1)*nxy, i), fdata(1, 1, iel, 2), nxy)
                  call copy(bufs(1 + (ie - 1)*nxy, i), fdata(1, 1, iel, 3), nxy)
               end do
            end do

            if (nid == 0) call byte_close(ierr)
            call bcast(ierr, isize)
            if (ierr /= 0) call nek_stop_error('Error closing '//trim(fname), this_module, this_procedure)

            nread = nsaver
            write (msg, '(A,I0,A,A)') 'Read ', nread, ' snapshots from ', trim(fname)
            call nek_log_debug(msg, this_module, this_procedure)

            deallocate (elmap, gmap_index, xr, yr, fdata)
            call neklab_timer_stop(t_bf_read)
         end procedure bf_read_chunk

      !====================================================================
      !     HEADER PEEK
      !====================================================================

         module procedure bf_peek_chunk
            character(len=*), parameter :: this_procedure = 'bf_peek_chunk'
            integer :: hdrsize, wdsl, isl, wdsizr
            integer :: nxr, nyr, nelfr, lbufr, nsaver
            real*4 :: test_pattern
            real(dp) :: timer(1)
            character(len=132) :: hdr
            character(len=256) :: msg
            character(len=4) :: sdummy
            logical, external :: if_byte_swap_test
            logical :: exists

            hdrsize = 116
            isl = isize/4
            nread = 0
            nsaver = 0
            ierr = 0

            if (nid == 0) then
               inquire (file=trim(fname), exist=exists)
               if (.not. exists) then
                  ierr = 1
               else
                  call byte_open(fname, ierr)
               end if
            end if
            call bcast(ierr, isize)
            if (ierr /= 0) return

            if (nid == 0) then
               call blank(hdr, hdrsize)
               call byte_read(hdr, hdrsize/4, ierr)
               if (ierr == 0) then
                  call byte_read(test_pattern, 1, ierr)
                  if_byte_sw = if_byte_swap_test(test_pattern, ierr)
               end if
            end if
            call bcast(ierr, isize)
      ! A failed header read leaves hdr blank; a list-directed read of that
      ! aborts, so bail out before touching it.
            if (ierr /= 0) then
               if (nid == 0) call byte_close(ierr)
               ierr = 1
               return
            end if

            if (nid == 0) then
               read (hdr, *) sdummy, wdsizr
               wdsl = wdsizr/4
               call byte_read(nxr, isl, ierr)
               call byte_read(nyr, isl, ierr)
               call byte_read(nelfr, isl, ierr)
               call byte_read(timer(1), wdsl, ierr)
               call byte_read(nsaver, isl, ierr)
               call byte_read(lbufr, isl, ierr)
               if (if_byte_sw) then
                  call byte_reverse(nxr, 1, ierr)
                  call byte_reverse(nyr, 1, ierr)
                  call byte_reverse(nelfr, 1, ierr)
                  call byte_reverse(nsaver, 1, ierr)
                  call byte_reverse(lbufr, 1, ierr)
                  call reverse_real(timer(1), 1, wdsizr, ierr)
               end if
               write (msg, '(A,A,3(1X,I0),1X,E15.7,2(1X,I0))') trim(fname), ':',
     &            nxr, nyr, nelfr, timer(1), nsaver, lbufr
               call nek_log_debug(msg, this_module, this_procedure)
               call byte_close(ierr)
            end if
            call bcast(ierr, isize)
            call bcast(nsaver, isize)
            nread = nsaver
         end procedure bf_peek_chunk

      !====================================================================
      !     HELPERS
      !====================================================================

         subroutine reverse_real(buf, nwords, wds, ierr)
      !! Endian-swap nwords reals of width wds bytes.
      !!
      !! Nek's byte_reverse reverses 4-BYTE words, so applying it to real*8
      !! data scrambles the doubles no matter what count is passed:
      !! byte_reverse(x, n) touches only half the array and byte_reverse(x, 2n)
      !! swaps each half of every double independently. byte_reverse8 is the
      !! routine for 8-byte quantities. (helix_io has the first form of this
      !! bug; it only shows up on cross-endian files.)
            real(dp), intent(inout) :: buf(*)
            integer, intent(in) :: nwords, wds
            integer, intent(inout) :: ierr
            if (wds == 8) then
               call byte_reverse8(buf, nwords, ierr)
            else
               call byte_reverse(buf, nwords, ierr)
            end if
         end subroutine reverse_real

         subroutine gather_and_write(fld, ngown)
      !! Ships every rank's elements to nid 0, which appends them to the open
      !! file in rank order. The element map written earlier records which
      !! global element each one is, so the ordering here is arbitrary and the
      !! result is independent of the rank count.
            real(dp), intent(in) :: fld(:)
            integer, intent(in) :: ngown(:)
      ! internal
            integer :: nxy, idum, wdsl, length, ierr, ip
      ! Sized on lelv, not nelv: nid 0 receives from every rank and any of them
      ! may own more elements than it does.
            real(dp) :: rtmp(lx1*ly1*lelv)

            nxy = lx1*ly1
            wdsl = wdsize/4
            ierr = 0

            if (nid == 0) then
               length = nxy*ngown(1)
               call byte_write(fld, length*wdsl, ierr)
               do ip = 1, np - 1
                  length = nxy*ngown(ip + 1)
                  call csend(ip, idum, isize, ip, 0)         ! hand shake
                  call crecv2(ip, rtmp, length*wdsize, ip)
                  call byte_write(rtmp, length*wdsl, ierr)
               end do
            else
               call crecv2(nid, idum, isize, 0)              ! hand shake
               length = nxy*ngown(nid + 1)
               call csend(nid, fld, length*wdsize, 0, 0)
            end if
            call bcast(ierr, isize)
            if (ierr /= 0) call nek_stop_error('Error writing field data', this_module, 'gather_and_write')
         end subroutine gather_and_write

      end submodule t2Dh_bf_buffer_io