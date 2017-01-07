c-----------------------------------------------------------------------
      subroutine PCGUARD32L(f,scs,kstrt,nvpy,nvpz,nxv,nypmx,nzpmx,mblok,
     1nblok,kyp,kzp,ngds,tag1,tag2,rid,sid,ierr)
c this subroutine copies data from field to particle partitions, copying
c data to guard cells, where the field and particle partitions are 
c assumed to be the same.
c f(j,k,l,m) = real data for grid j,k,l in particle partition m.  the
c number of grids per partition is uniform and includes one extra guard
c cell.
c scs(j,k,m) = scratch array for particle partition m
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nxv = first dimension of f, must be >= nx
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c mblok/nblok = number of particle partitions in y/z
c kyp/kzp = number of complex grids in y/z for each field partition.
c ngds = number of guard cells
c linear interpolation, for distributed data,
c with 2D spatial decomposition
      implicit none
      integer kstrt, nvpy, nvpz, nxv, nypmx, nzpmx, mblok, nblok, ngds
      integer kyp, kzp, tag1, tag2, rid, sid, ierr
      real f, scs
      dimension f(nxv,nypmx,nzpmx,mblok*nblok)
      dimension scs(nxv,nzpmx,2*ngds,mblok*nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus
      integer ky, kz, js, ks, moff, noff, kr, kl, mnblok
      integer nxvz, nxvy, m, my, mz, j, k
      dimension istatus(lstat)
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      noff = nypmx*nzpmx
      mnblok = mblok*nblok
      nxvz = nxv*nzpmx
      nxvy = nxv*nypmx
c copy to guard cells in z
      do 130 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 120 my = 1, mblok
      m = my + moff
      ky = my + js + 1
      kz = mz + ks
      kr = kz + 1
c      if (kr.ge.nvpz) kr = kr - nvpz
      kl = kz - 1
c      if (kl.lt.0) kl = kl + nvpz
c     kr = ky + nvpy*kr
c     kl = ky + nvpy*kl
c this segment is used for mpi computers
      if (kr.lt.nvpz) then
         call MPI_IRECV(f(1,1,kzp+1,m),nxvy,mreal,ky+nvpy*kr-1,tag1,lwor
     1ld,rid,ierr)
      endif
      if (kl.ge.0) then
         call MPI_ISEND(f(1,1,1,m),nxvy,mreal,ky+nvpy*kl-1,tag2,lworld,s
     1id,ierr)
      endif
  120 continue
  130 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PACGUARD32L(f,scs,kstrt,nvpy,nvpz,nx,nxv,nypmx,nzpmx,mb
     1lok,nblok,kyp,kzp,ngds)
c this subroutine copies data from particle to field partitions, adding
c data from guard cells, where the field and particle partitions are 
c assumed to be the same.
c f(3,j,k,l,m) = real data for grid j,k,l in particle partition m.  the
c number of grids per partition is uniform and includes one extra guard
c cell.
c scs = scratch array for particle partition m
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx+1
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c mblok/nblok = number of particle partitions in y/z
c kyp/kzp = number of complex grids in y/z for each field partition.
c ngds = number of guard cells
c linear interpolation, for distributed data
c with 2D spatial decomposition
      implicit none
      integer kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx, mblok, nblok
      integer ngds, kyp, kzp
      real f, scs
      dimension f(3,nxv,nypmx,nzpmx,mblok*nblok)
      dimension scs(3,nxv,nzpmx,2*ngds,mblok*nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer ierr, msid, istatus
      integer ky, kz, js, ks, moff, noff, kr, kl, mnblok
      integer nx1, kyp1, kzp1, nxvz, nxvy, m, my, mz, j, k, n
      dimension istatus(lstat)
      nx1 = nx + 1
      kyp1 = kyp + 1
      kzp1 = kzp + 1
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      noff = nypmx*nzpmx
      mnblok = mblok*nblok
      nxvz = nxv*nzpmx
      nxvy = nxv*nypmx
c special case for one processor in y
      if (nvpy.eq.1) then
         do 40 m = 1, mnblok
         do 30 k = 1, kzp1
         do 20 j = 1, nx1
         do 10 n = 1, 3
         f(n,j,1,k,m) = f(n,j,1,k,m) + f(n,j,kyp+1,k,m)
   10    continue
   20    continue
   30    continue
   40    continue
         go to 170
      endif
c buffer data in y
      do 80 m = 1, mnblok
      do 70 k = 1, nzpmx
      do 60 j = 1, nxv
      do 50 n = 1, 3
      scs(n,j,k,1,m) = f(n,j,kyp+1,k,m)
   50 continue
   60 continue
   70 continue
   80 continue
c add guard cells in y
      do 160 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 150 my = 1, mblok
      m = my + moff
      ky = my + js
      kz = nvpy*(mz + ks) + 1
      kr = ky + 1
c      if (kr.ge.nvpy) kr = kr - nvpy
      kl = ky - 1
c      if (kl.lt.0) kl = kl + nvpy
c      kr = kr + kz
c      kl = kl + kz
c this segment is used for shared memory computers
c     do 110 k = 1, nzpmx
c     do 100 j = 1, nxv
c     do 90 n = 1, 3
c     scs(n,j,k,2,m) = scs(n,j,k,1,kl)
c  90 continue
c 100 continue
c 110 continue
c this segment is used for mpi computers
      if (kl.ge.0) then
         call MPI_IRECV(scs(1,1,1,2,m),3*nxvz,mreal,kl+kz-1,noff+1,lgrp,
     1msid,ierr)
      endif
      if (kr.lt.nvpy) then
         call MPI_SEND(scs(1,1,1,1,m),3*nxvz,mreal,kr+kz-1,noff+1,lgrp,i
     1err)
      endif
      if (kl.ge.0) then
         call MPI_WAIT(msid,istatus,ierr)
         do 140 k = 1, kzp1
         do 130 j = 1, nx1
         do 120 n = 1, 3
         f(n,j,1,k,m) = f(n,j,1,k,m) + scs(n,j,k,2,m)
  120    continue
  130    continue
  140    continue
      endif
  150 continue
  160 continue
  170 return
      end
c-----------------------------------------------------------------------
      subroutine PAGUARD32L(f,scs,scr,kstrt,nvpy,nvpz,nx,nxv,nypmx,nzpmx
     1,mblok,nblok,kyp,kzp,ngds,tag1,tag2,id,ierr)
c this subroutine copies data from particle to field partitions, adding
c data from guard cells, where the field and particle partitions are 
c assumed to be the same.
c f(j,k,l,m) = real data for grid j,k,l in particle partition m.  the
c number of grids per partition is uniform and includes one extra guard
c cell.
c scs/scr = scratch array for particle partition m
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx+1
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c mblok/nblok = number of particle partitions in y/z
c kyp/kzp = number of complex grids in y/z for each field partition.
c ngds = number of guard cells
c linear interpolation, for distributed data
c with 2D spatial decomposition
      implicit none
      integer kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx, mblok, nblok
      integer ngds, kyp, kzp, tag1, tag2, id, ierr
      real f, scs, scr
      dimension f(nxv,nypmx,nzpmx,mblok*nblok)
      dimension scs(nxv,nzpmx,2*ngds,mblok*nblok)
      dimension scr(nxv,nypmx,mblok*nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer msid, istatus
      integer ky, kz, js, ks, moff, noff, kr, kl, mnblok
      integer nx1, kyp1, kzp1, nxvz, nxvy, m, my, mz, j, k
      dimension istatus(lstat)
      nx1 = nx + 1
      kyp1 = kyp + 1
      kzp1 = kzp + 1
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      noff = nypmx*nzpmx
      mnblok = mblok*nblok
      nxvz = nxv*nzpmx
      nxvy = nxv*nypmx
c special case for one processor in y
      if (nvpy.eq.1) then
         do 30 m = 1, mnblok
         do 20 k = 1, kzp1
         do 10 j = 1, nx1
         f(j,1,k,m) = f(j,1,k,m) + f(j,kyp+1,k,m)
   10    continue
   20    continue
   30    continue
         go to 130
      endif
c buffer data in y
      do 60 m = 1, mnblok
      do 50 k = 1, nzpmx
      do 40 j = 1, nxv
      scs(j,k,1,m) = f(j,kyp+1,k,m)
   40 continue
   50 continue
   60 continue
c add guard cells in y
      do 120 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 110 my = 1, mblok
      m = my + moff
      ky = my + js
      kz = nvpy*(mz + ks) + 1
      kr = ky + 1
c      if (kr.ge.nvpy) kr = kr - nvpy
      kl = ky - 1
c      if (kl.lt.0) kl = kl + nvpy
c      kr = kr + kz
c      kl = kl + kz
c this segment is used for shared memory computers
c     do 80 k = 1, nzpmx
c     do 70 j = 1, nxv
c     scs(j,k,2,m) = scs(j,k,1,kl)
c  70 continue
c  80 continue
c this segment is used for mpi computers
      if (kl.ge.0) then
         call MPI_IRECV(scs(1,1,2,m),nxvz,mreal,kl+kz-1,noff+1,lworld,ms
     1id,ierr)
      endif
      if (kr.lt.nvpy) then
         call MPI_SEND(scs(1,1,1,m),nxvz,mreal,kr+kz-1,noff+1,lworld,ier
     1r)
      endif
      if (kl.ge.0) then
         call MPI_WAIT(msid,istatus,ierr)
         do 100 k = 1, kzp1
         do 90 j = 1, nx1
         f(j,1,k,m) = f(j,1,k,m) + scs(j,k,2,m)
   90    continue
  100    continue
      endif
  110 continue
  120 continue
c special case for one processor in z
  130 if (nvpz.eq.1) then
         do 160 m = 1, mnblok
         do 150 k = 1, kyp1
         do 140 j = 1, nx1
         f(j,k,1,m) = f(j,k,1,m) + f(j,k,kzp+1,m)
  140    continue
  150    continue
  160    continue
         return
      endif
c add guard cells in z
      do 220 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 210 my = 1, mblok
      m = my + moff
      ky = my + js + 1
      kz = mz + ks
      kr = kz + 1
c      if (kr.ge.nvpz) kr = kr - nvpz
      kl = kz - 1
c      if (kl.lt.0) kl = kl + nvpz
c      kr = ky + nvpy*kr
c      kl = ky + nvpy*kl
c this segment is used for shared memory computers
c     do 180 k = 1, nypmx
c     do 170 j = 1, nxv
c     scr(j,k,m) = f(j,k,kzp+1,kl)
c 170 continue
c 180 continue
c this segment is used for mpi computers
      if (kl.ge.0) then
         call MPI_IRECV(scr,nxvy,mreal,ky+nvpy*kl-1,tag1,lworld,msid,ier
     1r)
      endif
      if (kr.lt.nvpz) then
         call MPI_ISEND(f(1,1,kzp+1,m),nxvy,mreal,ky+nvpy*kr-1,tag2,lwor
     1ld,id,ierr)
      endif
      if (kl.ge.0) then
         call MPI_WAIT(msid,istatus,ierr)
         do 200 k = 1, kyp1
         do 190 j = 1, nx1
         f(j,k,1,m) = f(j,k,1,m) + scr(j,k,m)
  190    continue
  200    continue
      endif
  210 continue
  220 continue
      return
      end
c-----------------------------------------------------------------------
      