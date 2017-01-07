c-----------------------------------------------------------------------
      subroutine PWRITE2(f,nx,kyp,nxv,kypmx,nblok,iunit,nrec,lrec,name)
c this subroutine collects distributed real 2d data f and writes to a
c direct access binary file
c f = input data to be written, modified on node 0
c nx/kyp = length of data f in x/y on each processor to write
c nxv = first dimension of data array f, must be >= nx
c kypmx = second dimension of data array f, must be >= kyp
c nblok = number of parallel partitions
c iunit = fortran unit number
c nrec = current record number for write, if nrec > 0
c if nrec < 0, open new file and write first record
c if nrec = 0, open old file, do not write
c lrec = record length (used only if nrec <= 0)
c name = file name (used only if nrec <= 0)
c input: f, nx, kyp, nxv, kypmx, nblok, iunit, nrec, lrec, fname
c output: nrec
      implicit none
      integer nx, kyp, nxv, kypmx, nblok, iunit, nrec, lrec
      real f
      character*(*) name
      dimension f(nxv,kypmx,nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c nproc = number of real or virtual processors obtained
c lgrp = current communicator
c lworld = MPI_COMM_WORLD communicator
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, nvp, idproc, np, ioff, id, nrec0, i, j, k, l
      integer ierr
      dimension istatus(lstat)
c determine the rank of the calling process in the communicator
      call MPI_COMM_RANK(lgrp,idproc,ierr)
c determine the size of the group associated with a communicator
      call MPI_COMM_SIZE(lgrp,nvp,ierr)
c node 0 receives messages from other nodes
      if (idproc.eq.0) then
         if (nrec.lt.0) then
            open(unit=iunit,file=name,form='unformatted',access='direct'
     1,recl=lrec,status='replace')
            nrec = 1
c open old file
         else if (nrec.eq.0) then
            open(unit=iunit,file=name,form='unformatted',access='direct'
     1,recl=lrec,status='old')
         endif
c no special diagnostic node
         if (nvp.eq.nproc) then
            np = nvp
            ioff = 1
c special diagnostic node present
         else
            np = nvp - nproc
            ioff = 0
            id = 1
            call MPI_RECV(f,nxv*kyp,mreal,id,99,lgrp,istatus,ierr)
         endif
c first write data for node 0
         nrec0 = nrec
         write (unit=iunit,rec=nrec) (((f(j,k,l),j=1,nx),k=1,kyp),l=1,nb
     1lok)
         nrec = nrec + 1
c then write data from remaining nodes
         do 10 i = 2, np
            id = i - ioff
            call MPI_RECV(f,nxv*kyp,mreal,id,99,lgrp,istatus,ierr)
            write (unit=iunit,rec=nrec) (((f(j,k,l),j=1,nx),k=1,kyp),l=1
     1,nblok)
            nrec = nrec + 1
   10    continue
c read data back for node 0
         read (unit=iunit,rec=nrec0) (((f(j,k,l),j=1,nx),k=1,kyp),l=1,nb
     1lok)
c other nodes send data to node 0
      elseif (idproc.le.(nproc+1)) then
         call MPI_SEND(f,nxv*kyp,mreal,0,99,lgrp,ierr)
      endif
      return
      end
