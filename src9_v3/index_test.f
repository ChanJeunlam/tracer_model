      program index_test
      implicit none
      integer i,j,k,nrecl, npad, ia, ja
      integer, parameter :: nbdy=2,ii=1475,jj=1950,ms=99
      integer, dimension(1-nbdy:ii+nbdy,1-nbdy:jj+nbdy) :: ip,iu,iv
      integer, dimension(1-nbdy:jj+nbdy,ms) :: ifp,ilp,ifu,ilu,ifv,ilv
      integer, dimension(1-nbdy:jj+nbdy) :: isp, isu, isv
      real, dimension(1-nbdy:ii+nbdy,1-nbdy:jj+nbdy) :: depths
      real, allocatable :: pad(:)
c
      npad=4096-mod((ii+nbdy)*(jj+nbdy),4096)
      allocate (pad(npad))
      inquire(iolength=nrecl)depths,pad
      open(13,file='../regional.depth.a',access='direct',
     $     recl=nrecl,form='unformatted',status='unknown',
     $     convert='big_endian')
      read(13,rec=1) depths
c
      do j=1-nbdy,jj+nbdy
        do i=1-nbdy,ii+nbdy
          ip(i,j)=0
          iu(i,j)=0
          iv(i,j)=0
        enddo
      enddo
c
c --- mass points are defined where water depth is greater than zero
c --- but for GLBa0.08 there are large numbers over land to be masked later.
c
      do j=1-nbdy,jj+nbdy
         do i=1-nbdy,ii+nbdy
            if (depths(i,j).gt.0.0) then
               ip(i,j)=1
            endif
         enddo
      enddo
c
c --- u,v points are located halfway between any 2 adjoining mass points
c
      do j=1-nbdy,jj+nbdy
        do i=1-nbdy,ii+nbdy
           ia=max(1-nbdy,i-1)
           if (ip(ia,j).gt.0.and.ip(i,j).gt.0) then
              iu(i,j)=1
           endif
           ja=max(1-nbdy,j-1)
           if (ip(i,ja).gt.0.and.ip(i,j).gt.0) then
              iv(i,j)=1
           endif
        enddo
      enddo
c
      do j= 1-nbdy,jj+nbdy
        do i= 1-nbdy,ii+nbdy
           if     (depths(i,j).gt.1.e10) depths(i,j) = 0.0
        enddo
      enddo
c recompute ip, iu,iv
      do j=1-nbdy,jj+nbdy
         do i=1-nbdy,ii+nbdy
            ip(i,j)=0
            iu(i,j)=0
            iv(i,j)=0
            if (depths(i,j).gt.0.)then
               ip(i,j)=1
            endif
         enddo
      enddo
c --- u,v points are located halfway between any 2 adjoining mass points
      do j=1-nbdy,jj+nbdy
        do i=1-nbdy,ii+nbdy
           ia=max(1-nbdy,i-1)
           if (ip(ia,j).gt.0.and.ip(i,j).gt.0) then
              iu(i,j)=1
           endif
           ja=max(1-nbdy,j-1)
           if (j.gt.1.and.ip(i,ja).gt.0.and.ip(i,j).gt.0) then
              iv(i,j)=1
           endif
        enddo
      enddo
c --- determine loop indices for mass and velocity points
      call indxi(ip,ifp,ilp,isp)
c      call indxj(ip,jfp,jlp,jsp)
      call indxi(iu,ifu,ilu,isu)
c      call indxj(iu,jfu,jlu,jsu)
      call indxi(iv,ifv,ilv,isv)
c      call indxj(iv,jfv,jlv,jsv)
c
      j=1000
      print *,'ifp,ilp,ifu,ilu,ifv,ilv j=1000'
      print *, ifp(1000,1:isp(j))
      print *, ilp(1000,1:isp(j))
      print *, ifu(1000,1:isu(j))
      print *, ilu(1000,1:isu(j))
      print *, ifv(1000,1:isv(j))
      print *, ilv(1000,1:isv(j))
c
      stop
      end


      subroutine indxi(ipt,if,il,is)
      implicit none
      integer i,j,k
      integer, parameter :: nbdy=2,ii=1475,jj=1950,ms=99
c
      integer, dimension (1-nbdy:ii+nbdy,1-nbdy:jj+nbdy) :: ipt
      integer, dimension (1-nbdy:jj+nbdy,ms) :: if,il
      integer, dimension (1-nbdy:jj+nbdy) :: is
c
c --- input array ipt contains 1 at grid point locations, 0 elsewhere
c --- output is arrays if, il, is  where
c --- if(j,k) gives row index of first point in column j for k-th section
c --- il(j,k) gives row index of last point
c --- is(j) gives number of sections in column j (maximum: ms)
c
      integer last,i0,j0
c
      do j=1-nbdy,jj+nbdy
         is(j) = 0
         do k=1,ms
            if(j,k) = 0
            il(j,k) = 0
         end do
c
         k=1
         last = ipt(1-nbdy,j)
         if     (last .eq. 1) then
            if(j,k) = 1-nbdy
         endif
         do i=2-nbdy,ii+nbdy
            if      (last .eq. 1 .and. ipt(i,j) .eq. 0) then
               il(j,k) = i-1
               k = k+1
            elseif (last .eq. 0 .and. ipt(i,j) .eq. 1) then
               if     (k .gt. ms) then
c                  write(lp,'(a,i5)')  'indxi problem on proc ',mnproc
c                  write(lp,'(a,2i5)')
c     &   ' error in indxi -- ms too small at i,j =',i0+i,j0+j
c                  call xchalt('(indxi)')
                  stop '(indxi)'
               endif
               if(j,k) = i
            endif
            last = ipt(i,j)
         enddo
         if     (last .eq. 1) then
            il(j,k) = ii+nbdy
            is(j) = k
         else
            is(j) = k-1
         endif
      enddo
c      call xcsync(no_flush)
      return
      end
