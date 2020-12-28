!     ifort -g -O3 -convert big_endian -openmp -o calc_mean_uvdp_HRA.f
      program calc_mean_uvdp_HRA
!     
!     The program reads the variables needed for the offline model
!     and saves them in the forcing files.

!     The forcing fields are calc by time-averaging and spatial filtering
!     'make_uvdp_HRA.f' only performs time-averaging.
!
      implicit none
      character*4 var_name
c
      var_name = 'dpii'
      call process_var(var_name)
      var_name='uflx' 
      call process_var(var_name)
      var_name='vflx'
      call process_var(var_name)
      var_name='dpmm'
      call process_var(var_name)
      stop
      end
c
c=====================================================================
c          SUBTOUTINE: process_var
c=====================================================================
      subroutine process_var(var_name)
      implicit none
      integer nday,i,j,k,days,daye
      integer nrec,npad,nrecl
      character*3 dayn
      character*2 hourn
      integer, parameter :: isz=1573,jsz=1073,ksz=30,nvar=3
      real, parameter :: lflag=1.0e10
      integer, parameter :: nfilx=25, nfily=25
      real, dimension(isz,jsz) :: var, varc, sc2, scux, scuy, scvx,
     $     scvy, scpx, scpy
      real, allocatable :: pad(:)
      character*4 var_name
      character*100 foc_path, grd_path, sav_path
      integer lp
c
      npad=4096-mod(isz*jsz,4096)
      allocate (pad(npad))
      inquire(iolength=nrecl)var,pad
c
      grd_path = 'MODEL_GRID_PATH'
      atl_path = 'MODEL_FORCING_PATH'
      sav_path = 'OUTPUT_PATH'
c
      lp = 6
      days = daysSh
      daye = dayeSh
c
c     ================ Read grid (fid=31)
c
      open(31,file=trim(grd_path)//'regional.grid.a', access='direct',
     $     recl=nrecl, form='unformatted', status='old',
     $     convert='big_endian')
      write(lp,*)  'READING '//trim(grd_path)//'regional.grid.a'
c
      read(31,rec=10) scpx
      read(31,rec=11) scpy
      read(31,rec=14) scux
      read(31,rec=15) scuy
      read(31,rec=16) scvx
      read(31,rec=17) scvy
      close(31)
c
c     ================== Read UVDP (fid=21)
c
      print*,'PROCESSING ',var_name 
c
c     ---- loop over time
c
      DO nday = 2*days,2*daye-1

        write(dayn,'(i3.3)') (nday+1)/2 
        write(hourn,'(i2.2)') 12*(nday+1-2*((nday+1)/2))
c
c       ---- open file to be READ
        write(lp,*) 'READING '//trim(foc_path)//'uvdp_offline.'//dayn//
     $  '_'//hourn//'.a'
        open(21, file=trim(foc_path)//'uvdp_offline.'//dayn//'_'//
     $  hourn//'.a', access='direct', recl=nrecl, convert='big_endian',
     $  form='unformatted',status='unknown')
c
c       ---- open file to be SAVED
        write(lp,*) 'TO SAVE '//trim(sav_path)//'uvdp_offline.'//dayn//
     $    '_'//hourn//'.a' 
        open(12,file=trim(sav_path)//'uvdp_offline.'//dayn//'_'//
     $  hourn//'.a', access='direct', recl=nrecl, convert='big_endian',
     $  form='unformatted', status='unknown')
c
c       ---- layer by layer
c
!OMP PARALLEL DO PRIVATE(k)
        do k = 1,ksz
          print*,'K level ', k
          if (var_name.eq.'dpii') nrec = k
          if (var_name.eq.'uflx') nrec = ksz + nvar*(k-1) + 1
          if (var_name.eq.'vflx') nrec = ksz + nvar*(k-1) + 2
          if (var_name.eq.'dpmm') nrec = ksz + nvar*(k-1) + 3
          read(21,rec=nrec) var
c
c         ---- Spatially running average and save (fid=12)
c
          if (var_name.eq.'uflx') then
            sc2=scux*scuy
          elseif (var_name.eq.'vflx') then
            sc2=scvx*scvy
          else
            sc2=scpx*scpy
          endif
          call coarsen(var,nfilx,nfily,sc2,varc)  
          write(12,rec=nrec) ((varc(i,j),i=1,isz),j=1,jsz)
        enddo   ! end of k-loop
!OMP END PARALLEL DO
         close(12)
         close(21)
      ENDDO       ! end of nday-loop
      return
      end
c=====================================================================
c          SUBTOUTINE: spatially coarsen (similar to 'smooth_subgrid')
c===================================================================== 
      subroutine coarsen(varin,nfilx,nfily,sc2,varout)
      implicit none
      integer, parameter :: isz=1573,jsz=1073 ! include buffer points
      integer nfilx, nfily, i,j,ic,jc,is,ie,js,je
      real, dimension(isz,jsz) :: varin, varout, sc2
      real varsum, arsum

!OMP PARALLEL DO PRIVATE(i,j,ic,jc,varsum,arsum,is,js,ie,je)
      do i=1,isz
         is=max(1,  i-nfilx)
         ie=min(isz,i+nfilx)
         do j=1,jsz
c       make sure you do not create extra points
            varout(i,j)=0.0
            if  (varin(i,j).lt.1.0e10) then
              js=max(1,  j-nfily)
              je=min(jsz,j+nfily)
              varsum=0.0
              arsum =0.0
              do ic=is,ie
                do jc=js,je
                  if (varin(ic,jc).lt.1.0e10) then
                     varsum = varsum + varin(ic,jc)*sc2(ic,jc)
                     arsum = arsum + sc2(ic,jc)
                  endif
                enddo
              enddo
              if (arsum.gt.0.0) varout(i,j) = varsum/arsum
            endif
         enddo
      enddo
!OMP END PARALLEL DO
      return
      end
