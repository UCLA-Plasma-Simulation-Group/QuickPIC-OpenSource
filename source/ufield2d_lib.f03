! ufield2d_lib module for QuickPIC Open Source 1.0
! update: 04/18/2016

      module ufield2d_lib

      use mpi
         
      implicit none

      public

      contains
!
      subroutine PPNCGUARD2L(f,nyp,kstrt,nvp,nxv,nypmx,lgrp,mreal)
! this subroutine copies data to guard cells in non-uniform partitions
! f(j,k) = real data for grid j,k in particle partition.
! the grid is non-uniform and includes one extra guard cell.
! output: f
! nyp = number of primary gridpoints in field partition
! it is assumed the nyp > 0.
! kstrt = starting data block number
! nvp = number of real or virtual processors
! nxv = first dimension of f, must be >= nx
! nypmx = maximum size of field partition, including guard cell.
! linear interpolation, for distributed data
! lgrp = current communicator
! mreal = default datatype for reals
      implicit none
      integer, intent(in) :: lgrp, mreal
      integer, intent(in) :: nyp, kstrt, nvp, nxv, nypmx
      real, dimension(nxv,nypmx), intent(inout) :: f
! local data
      integer :: j, ks, moff, kl, kr
      integer :: msid, ierr
      integer, dimension(MPI_STATUS_SIZE) :: istatus
! special case for one processor
      if (nvp==1) then
         do j = 1, nxv
            f(j,nyp+1) = f(j,1)
         enddo
         return
      endif
      ks = kstrt - 1
      moff = nypmx*nvp + 2
! copy to guard cells
      kr = ks + 1
      kl = ks - 1
      ks = nyp + 1
! this segment is used for mpi computers
      if (kr < nvp) then
         call MPI_IRECV(f(1,ks),nxv,mreal,kr,moff,lgrp,msid,ierr)
      end if
      if (kl >= 0) then
         call MPI_SEND(f,nxv,mreal,kl,moff,lgrp,ierr)
      end if
      if (kr < nvp) then
         call MPI_WAIT(msid,istatus,ierr)
      end if
      end subroutine
!      
      subroutine PPNACGUARD2L(f,scr,nyp,nx,ndim,kstrt,nvp,nxv,nypmx,lgrp,mreal)
! this subroutine adds data from guard cells in non-uniform partitions
! f(ndim,j,k) = real data for grid j,k in particle partition.
! the grid is non-uniform and includes one extra guard cell.
! output: f, scr
! scr(ndim,j) = scratch array for particle partition
! nyp = number of primary gridpoints in particle partition
! it is assumed the nyp > 0.
! kstrt = starting data block number
! nvp = number of real or virtual processors
! nx = system length in x direction
! ndim = leading dimension of array f
! nxv = first dimension of f, must be >= nx
! nypmx = maximum size of field partition, including guard cells.
! linear interpolation, for distributed data
! lgrp = current communicator
! mreal = default datatype for reals
      implicit none
      integer, intent(in) :: lgrp, mreal
      integer, intent(in) ::  nyp, kstrt, nvp, nx, ndim, nxv, nypmx
      real, dimension(ndim,nxv,nypmx), intent(inout) :: f
      real, dimension(ndim,nxv), intent(inout) :: scr
! local data
      integer :: j, n, nx1, ks, moff, kl, kr
      integer :: nnxv
      integer :: msid, ierr
      integer, dimension(MPI_STATUS_SIZE) :: istatus
      nx1 = nx + 1
! special case for one processor
      if (nvp==1) then
         do j = 1, nx1
            do n = 1, ndim
               f(n,j,1) = f(n,j,1) + f(n,j,nyp+1)
               f(n,j,nyp+1) = 0.0
            enddo
         enddo
         return
      endif
      ks = kstrt - 1
      moff = nypmx*nvp + 1
      nnxv = ndim*nxv
! add guard cells
      kr = ks + 1
      kl = ks - 1
      ks = nyp + 1
! this segment is used for mpi computers
      if (kl >= 0) then
         call MPI_IRECV(scr,nnxv,mreal,kl,moff,lgrp,msid,ierr)
      end if
      if (kr < nvp) then
         call MPI_SEND(f(1,1,ks),nnxv,mreal,kr,moff,lgrp,ierr)
      end if
      if (kl >= 0) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         scr(:,:) = 0.0
      end if
! add up the guard cells
      do j = 1, nx1
         do n = 1, ndim
            f(n,j,1) = f(n,j,1) + scr(n,j)
         enddo
      enddo

      if (kr < nvp) then
         do j = 1, nx1
            do n = 1, ndim
               f(n,j,nyp+1) = 0.0
            enddo
         enddo
      end if
      end subroutine 
!      
      end module ufield2d_lib
