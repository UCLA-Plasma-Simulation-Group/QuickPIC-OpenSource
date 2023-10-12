! part2d_lib_h module for QuickPIC Open Source 1.0
! update: 04/18/2016

      module part2d_lib

      use mpi
         
      implicit none

      public
!
      integer :: nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
!
      interface
         subroutine PPDBLKP2L(part,kpic,npp,noff,nppmx,idimp,npmax,mx,my&
     &,mx1,mxyp1,irc)
         implicit none
         integer, intent(in) :: idimp, npmax, mx, my, mx1, mxyp1, npp
         integer, intent(in) :: noff
         integer, intent(inout) :: nppmx, irc
         real, dimension(idimp,npmax), intent(in) :: part
         integer, dimension(mxyp1), intent(inout) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPPMOVIN2L(part,ppart,kpic,npp,noff,nppmx,idimp,    &
     &npmax,mx,my,mx1,mxyp1,irc)
         implicit none
         integer, intent(in) :: nppmx, idimp, npmax, mx, my, mx1, mxyp1
         integer, intent(in) :: npp, noff
         integer, intent(inout) :: irc
         real, dimension(idimp,npmax), intent(in) :: part
         real, dimension(idimp,nppmx,mxyp1), intent(inout) :: ppart
         integer, dimension(mxyp1), intent(inout) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPPCHECK2L(ppart,kpic,noff,nyp,idimp,nppmx,nx,mx,my,&
     &mx1,myp1,irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, mx, my, mx1, myp1
         integer, intent(in) :: noff, nyp
         integer, intent(inout) :: irc
         real, dimension(idimp,nppmx,mx1*myp1), intent(in) :: ppart
         integer, dimension(mx1*myp1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPGPPOST2L(ppart,q,kpic,noff,idimp,nppmx,mx,my,  &
     &nxv,nypmx,mx1,mxyp1)
         implicit none
         integer, intent(in) :: idimp, nppmx, mx, my, nxv, nypmx
         integer, intent(in) :: mx1, mxyp1
         integer, intent(in) :: noff
         real, dimension(idimp,nppmx,mxyp1), intent(in) :: ppart
         real, dimension(nxv,nypmx), intent(inout) :: q
         integer, dimension(mxyp1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
         subroutine PPGRDCJPPOST2L_QP(ppart,fxy,bxy,psit,cu,dcu,amu,kpic,noff,&
     &nyp,qbm,dt,ci,idimp,nppmx,nx,mx,my,nxv,nypmx,mx1,mxyp1,dex)
         implicit none
         integer, intent(in) :: noff, nyp, idimp, nppmx, nx, mx, my
         integer, intent(in) :: nxv, nypmx, mx1, mxyp1
         real, intent(in) ::  qbm, dt, ci
         real, dimension(idimp,nppmx,mxyp1), intent(inout) :: ppart
         real, dimension(2,nxv,nypmx), intent(in) :: fxy
         real, dimension(3,nxv,nypmx), intent(in) :: bxy
         real, dimension(nxv,nypmx), intent(in) :: psit
         real, dimension(3,nxv,nypmx), intent(inout) :: cu
         real, dimension(2,nxv,nypmx), intent(inout) :: dcu
         real, dimension(3,nxv,nypmx), intent(inout) :: amu
         integer, dimension(mxyp1), intent(in) :: kpic
         real, intent(in) :: dex
         end subroutine
      end interface
!
      interface
         subroutine PPGRBPPUSHF23L_QP(ppart,fxy,bxy,psit,kpic,ncl,ihole,&
     &noff,nyp,qbm,dt,dtc,ci,ek,idimp,nppmx,nx,ny,mx,my,nxv,nypmx,mx1,mx&
     &yp1,ntmax,irc,dex)
         implicit none
         integer noff, nyp, idimp, nppmx, nx, ny, mx, my, nxv, nypmx
         integer mx1, mxyp1, ntmax, irc
         real qbm, dt, dtc, ci, ek, dex
         real ppart, fxy, bxy, psit
         integer kpic, ncl, ihole
         dimension ppart(idimp,nppmx,mxyp1)
         dimension fxy(2,nxv,nypmx), bxy(3,nxv,nypmx)
         dimension psit(nxv,nypmx)
         dimension kpic(mxyp1), ncl(8,mxyp1)
         dimension ihole(2,ntmax+1,mxyp1)     
         end subroutine     
      end interface
!
      interface
         subroutine PPPORDER2LA(ppart,ppbuff,sbufl,sbufr,kpic,ncl,ihole,&  
     &ncll,nclr,noff,nyp,idimp,nppmx,nx,ny,mx,my,mx1,myp1,npbmx,ntmax,  &
     &nbmax,irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, ny, mx, my, mx1, myp1
         integer, intent(in) :: npbmx, ntmax, nbmax, noff, nyp
         integer, intent(inout) :: irc
         real, dimension(idimp,nppmx,mx1*myp1), intent(inout) :: ppart
         real, dimension(idimp,npbmx,mx1*myp1), intent(inout) :: ppbuff
         real, dimension(idimp,nbmax), intent(inout) :: sbufl, sbufr
         integer, dimension(mx1*myp1), intent(in) :: kpic
         integer, dimension(8,mx1*myp1), intent(inout) :: ncl
         integer, dimension(2,ntmax+1,mx1*myp1), intent(inout) :: ihole
         integer, dimension(3,mx1), intent(inout) :: ncll, nclr
         end subroutine
      end interface
!
      interface
         subroutine PPPORDERF2LA(ppart,ppbuff,sbufl,sbufr,ncl,ihole,ncll&
     &,nclr,idimp,nppmx,mx1,myp1,npbmx,ntmax,nbmax,irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, mx1, myp1, npbmx, ntmax
         integer, intent(in) :: nbmax
         integer, intent(inout) :: irc
         real, dimension(idimp,nppmx,mx1*myp1), intent(inout) :: ppart
         real, dimension(idimp,npbmx,mx1*myp1), intent(inout) :: ppbuff
         real, dimension(idimp,nbmax), intent(inout) :: sbufl, sbufr
         integer, dimension(8,mx1*myp1), intent(inout) :: ncl
         integer, dimension(2,ntmax+1,mx1*myp1), intent(in) :: ihole
         integer, dimension(3,mx1), intent(inout) :: ncll, nclr
         end subroutine
      end interface
!
      interface
         subroutine PPPORDER2LB(ppart,ppbuff,rbufl,rbufr,kpic,ncl,ihole,&
     &mcll,mclr,idimp,nppmx,mx1,myp1,npbmx,ntmax,nbmax,irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, mx1, myp1, npbmx, ntmax
         integer, intent(in) :: nbmax
         integer, intent(inout) :: irc
         real, dimension(idimp,nppmx,mx1*myp1), intent(inout) :: ppart
         real, dimension(idimp,npbmx,mx1*myp1), intent(in) :: ppbuff
         real, dimension(idimp,nbmax), intent(in) :: rbufl, rbufr
         integer, dimension(mx1*myp1), intent(inout) :: kpic
         integer, dimension(8,mx1*myp1), intent(in) :: ncl
         integer, dimension(2,ntmax+1,mx1*myp1), intent(in) :: ihole
         integer, dimension(3,mx1), intent(in) :: mcll, mclr
         end subroutine
      end interface
!
      interface
         subroutine WPGPSIPOST2L_QP(ppart,psi,kpic,qbm,noff,nyp,idimp,np&
     &pmx,nx,mx,my,nxv,nypmx,mx1,mxyp1,dex)
         implicit none
         integer, intent(in) :: noff, nyp, idimp, nppmx, nx, mx, my
         integer, intent(in) :: mx1, mxyp1, nxv, nypmx
         real, intent(in) :: dex,qbm
         real, dimension(idimp,nppmx,mxyp1), intent(inout) :: ppart
         real, dimension(nxv,nypmx), intent(in) :: psi
         integer, dimension(mxyp1), intent(in) :: kpic
         end subroutine
      end interface      
!
      interface
         subroutine PPPCOPYOUT2(part,ppart,kpic,npp,npmax,nppmx,idimp,  &
     &mxyp1,irc)
         implicit none
         integer, intent(in) :: npmax, nppmx, idimp, mxyp1
         integer, intent(inout) :: npp, irc
         real, dimension(idimp,npmax), intent(inout) :: part
         real, dimension(idimp,nppmx,mxyp1), intent(in) :: ppart
         integer, dimension(mxyp1), intent(in) :: kpic
         end subroutine
      end interface
!
      contains
!      
      subroutine PPPMOVE2(sbufr,sbufl,rbufr,rbufl,ncll,nclr,mcll,mclr,  &
     &kstrt,nvp,idimp,nbmax,mx1)
! this subroutine moves particles into appropriate spatial regions
! for distributed data, with 1d domain decomposition in y.
! tiles are assumed to be arranged in 2D linear memory
! output: rbufr, rbufl, mcll, mclr
! sbufl = buffer for particles being sent to lower processor
! sbufr = buffer for particles being sent to upper processor
! rbufl = buffer for particles being received from lower processor
! rbufr = buffer for particles being received from upper processor
! ncll = particle number being sent to lower processor
! nclr = particle number being sent to upper processor
! mcll = particle number being received from lower processor
! mclr = particle number being received from upper processor
! kstrt = starting data block number
! nvp = number of real or virtual processors
! idimp = size of phase space = 4 or 5
! nbmax =  size of buffers for passing particles between processors
! mx1 = (system length in x direction - 1)/mx + 1
      implicit none
      integer, intent(in) :: kstrt, nvp, idimp, nbmax, mx1
      real, dimension(idimp,nbmax), intent(in) :: sbufl, sbufr
      real, dimension(idimp,nbmax), intent(inout) :: rbufl, rbufr
      integer, dimension(3,mx1), intent(in) :: ncll, nclr
      integer, dimension(3,mx1), intent(inout) :: mcll, mclr
! lgrp = current communicator
! mint = default datatype for integers
! mreal = default datatype for reals
! local data
      integer :: ierr, ks, kl, kr, i, j, jsl, jsr
      integer :: nbsize, ncsize
      integer, dimension(8) :: msid
      integer, dimension(4) :: itg
      integer, dimension(10) :: istatus
      data itg /3,4,5,6/
      ks = kstrt - 1
      nbsize = idimp*nbmax
      ncsize = 3*mx1
! copy particle buffers: update rbufl, rbufr, mcll, mclr
! special case for one processor
      if (nvp==1) then
         do j = 1, mx1
            do i = 1, 3
               mcll(i,j) = nclr(i,j)
            enddo
         continue
         enddo
         do j = 1, mx1
            do i = 1, 3
               mclr(i,j) = ncll(i,j)
            enddo
         enddo
         do j = 1, nclr(3,mx1)
            do i = 1, idimp
               rbufl(i,j) = sbufr(i,j)
            enddo
         enddo
         do j = 1, ncll(3,mx1)
            do i = 1, idimp
               rbufr(i,j) = sbufl(i,j)
            enddo
         enddo
! this segment is used for mpi computers
      else
! get particles from below and above
         kr = ks + 1
         if (kr >= nvp) kr = kr - nvp
         kl = ks - 1
         if (kl < 0) kl = kl + nvp
! post receives
         call MPI_IRECV(mcll,ncsize,mint,kl,itg(1),lgrp,msid(1),ierr)
         call MPI_IRECV(mclr,ncsize,mint,kr,itg(2),lgrp,msid(2),ierr)
         call MPI_IRECV(rbufl,nbsize,mreal,kl,itg(3),lgrp,msid(3),ierr)
         call MPI_IRECV(rbufr,nbsize,mreal,kr,itg(4),lgrp,msid(4),ierr)
! send particle number offsets
         call MPI_ISEND(nclr,ncsize,mint,kr,itg(1),lgrp,msid(5),ierr)
         call MPI_ISEND(ncll,ncsize,mint,kl,itg(2),lgrp,msid(6),ierr)
         call MPI_WAIT(msid(1),istatus,ierr)
         call MPI_WAIT(msid(2),istatus,ierr)
! send particles
         jsr = idimp*nclr(3,mx1)
         call MPI_ISEND(sbufr,jsr,mreal,kr,itg(3),lgrp,msid(7),ierr)
         jsl = idimp*ncll(3,mx1)
         call MPI_ISEND(sbufl,jsl,mreal,kl,itg(4),lgrp,msid(8),ierr)
         call MPI_WAIT(msid(3),istatus,ierr)
         call MPI_WAIT(msid(4),istatus,ierr)
      endif
! make sure sbufr, sbufl, ncll, and nclr have been sent
      if (nvp /= 1) then
         do i = 1, 4
            call MPI_WAIT(msid(i+4),istatus,ierr)
         enddo
      endif
      end subroutine PPPMOVE2

subroutine PPGRBPPUSHF23L_RB( ppart, fxy, bxy, psit, kpic, ncl, ihole, noff, &
  nyp, qbm, dt, dtc, ci, ek, idimp, nppmx, nx, ny, mx, my, nxv, nypmx, mx1, &
  mxyp1, ntmax, irc, dex )

  ! see PPGRBPPUSHF23L_QP for the meaning of the input and output parameters

  implicit none

  real, intent(inout), dimension(idimp, nppmx, mxyp1) :: ppart
  real, intent(in), dimension(2, nxv, nypmx) :: fxy
  real, intent(in), dimension(3, nxv, nypmx) :: bxy
  real, intent(in), dimension(nxv, nypmx) :: psit
  integer, intent(in), dimension(mxyp1) :: kpic
  integer, intent(inout), dimension(8, mxyp1) :: ncl
  integer, intent(inout), dimension(2, ntmax+1, mxyp1) :: ihole
  integer, intent(in) :: noff, nyp, idimp, nppmx, nx, ny, mx, my, nxv, nypmx
  integer, intent(in) :: mx1, mxyp1, ntmax
  integer, intent(inout) :: irc
  real, intent(in) :: qbm, dt, dtc, ci, dex
  real, intent(inout) :: ek

  ! local variables
  real :: idex, qtmh, qtmh1, qtmh2, ex, ey, ez, bx, by, bz, px, py, pz, x, y
  real :: dxp, dyp, amx, amy, tmp_x, tmp_y, tmp_z, edgelx, edgely, edgerx, edgery, dtc1
  real :: sum1, sum2, anx, any, gam, ostq, p6, p7
  integer :: i, j, k, noffp, moffp, nppp, nn, mm, ih, nh, mnoff
  real, dimension(2, mx+1, my+1) :: sfxy
  real, dimension(3, mx+1, my+1) :: sbxy
  real, dimension(mx+1, my+1) :: spsit

  idex = 1.0 / dex
  qtmh = 0.5 * qbm * dt
  anx = real(nx)
  any = real(ny)
  sum2 = 0.0

! loop over tiles
!$omp parallel do &
!$omp& private(i, j, k, noffp, moffp, nppp, nn, mm, ih, nh, mnoff, x, y, dxp) &
!$omp& private(dyp, amx, amy, ex, ey, ez, bx, by, bz, tmp_x, tmp_y, tmp_z, edgelx, edgely, edgerx) &
!$omp& private(edgery, sum1, sfxy, sbxy, spsit, qtmh1, qtmh2, p6, p7, dtc1, gam, ostq) &
! $omp& shared(idex, qtmh, anx, any) &
! $omp& shared(ppart, fxy, bxy, psit, kpic, ncl, ihole, noff, nyp, dt, dtc, ek) &
! $omp& shared(nx, ny, mx, my, mx1, mxyp1, ntmax, irc, dex) &
!$omp& reduction(+:sum2)
  do k = 1, mxyp1

    noffp = (k - 1) / mx1
    moffp = my * noffp
    noffp = mx * (k - mx1 * noffp - 1)
    nppp = kpic(k)
    nn = min(mx, nx - noffp)
    mm = min(my, nyp - moffp)
    edgelx = noffp
    edgerx = noffp + nn
    edgely = noff + moffp
    edgery = noff + moffp + mm
    ih = 0
    nh = 0
    mnoff = moffp + noff - 1

    ! load local fields from global arrays
    do j = 1, mm + 1
      do i = 1, nn + 1
        sfxy(1,i,j) = fxy(1, i + noffp, j + moffp)
        sfxy(2,i,j) = fxy(2, i + noffp, j + moffp)
      enddo
    enddo

    do j = 1, mm + 1
      do i = 1, nn + 1
        sbxy(1,i,j) = bxy(1, i + noffp, j + moffp)
        sbxy(2,i,j) = bxy(2, i + noffp, j + moffp)
        sbxy(3,i,j) = bxy(3, i + noffp, j + moffp)
      enddo
    enddo

    do j = 1, mm + 1
      do i = 1, nn + 1
        spsit(i,j) = psit(i + noffp, j + moffp)
      enddo
    enddo

    ! clear counters
    do j = 1, 8
      ncl(j,k) = 0
    enddo
    sum1 = 0.0d0

    ! loop over particles in tile
    do j = 1, nppp

      ! find interpolation weights
      x = ppart(1,j,k)
      y = ppart(2,j,k)
      nn = x
      mm = y

      p6 = ppart(6,j,k)
      p7 = ppart(7,j,k)

      qtmh1 = qtmh / p7 
      qtmh2 = qtmh1 * p6

      dxp = x - real(nn)
      dyp = y - real(mm)
      nn = nn - noffp + 1
      mm = mm - mnoff
      amx = 1.0 - dxp
      amy = 1.0 - dyp

      ! find electric field
      ex = amx * sfxy(1,nn, mm)
      ey = amx * sfxy(2,nn, mm)
      ez = amx *  spsit(nn, mm)
      ex = amy * (dxp * sfxy(1,nn+1, mm) + ex)
      ey = amy * (dxp * sfxy(2,nn+1, mm) + ey)
      ez = amy * (dxp *  spsit(nn+1, mm) + ez)
      tmp_x = amx * sfxy(1,nn, mm+1)
      tmp_y = amx * sfxy(2,nn, mm+1)
      tmp_z = amx *  spsit(nn, mm+1)
      ex = ex + dyp * (dxp * sfxy(1,nn+1, mm+1) + tmp_x)
      ey = ey + dyp * (dxp * sfxy(2,nn+1, mm+1) + tmp_y)
      ez = ez + dyp * (dxp *  spsit(nn+1, mm+1) + tmp_z)

      ! E-field is in the unit of dex, i.e., (ex, ey, ez) / dex is normalized.
      ez = ez * idex

      ! find magnetic field
      bx = amx * sbxy(1,nn,mm)
      by = amx * sbxy(2,nn,mm)
      bz = amx * sbxy(3,nn,mm)
      bx = amy * (dxp * sbxy(1, nn+1, mm) + bx)
      by = amy * (dxp * sbxy(2, nn+1, mm) + by)
      bz = amy * (dxp * sbxy(3, nn+1, mm) + bz)
      tmp_x = amx * sbxy(1, nn, mm+1)
      tmp_y = amx * sbxy(2, nn, mm+1)
      tmp_z = amx * sbxy(3, nn, mm+1)
      bx = bx + dyp * (dxp * sbxy(1, nn+1, mm+1) + tmp_x)
      by = by + dyp * (dxp * sbxy(2, nn+1, mm+1) + tmp_y)
      bz = bz + dyp * (dxp * sbxy(3, nn+1, mm+1) + tmp_z)

      ! B-field is in normalized unit
      bx = bx * dex
      by = by * dex

      ! calculate half impulse
      ex = qtmh2 * ex
      ey = qtmh2 * ey
      ez = qtmh2 * ez

      ! half acceleration
      px = ppart(3,j,k) + ex
      py = ppart(4,j,k) + ey
      pz = ppart(5,j,k) + ez

      ! rotate about magnetic field
      bx = qtmh1 * bx
      by = qtmh1 * by
      bz = qtmh1 * bz

      tmp_x = px + py * bz - pz * by
      tmp_y = py + pz * bx - px * bz
      tmp_z = pz + px * by - py * bx

      ostq = 2.0 / ( 1.0 + bx**2 + by**2 + bz**2 )
      bx = bx * ostq
      by = by * ostq
      bz = bz * ostq

      ppart(3,j,k) = px + tmp_y * bz - tmp_z * by
      ppart(4,j,k) = py + tmp_z * bx - tmp_x * bz
      ppart(5,j,k) = pz + tmp_x * by - tmp_y * bx

      ! half acceleration
      ppart(3,j,k) = ppart(3,j,k) + ex
      ppart(4,j,k) = ppart(4,j,k) + ey
      ppart(5,j,k) = ppart(5,j,k) + ez

      ! advance particle position
      gam = sqrt( 1.0 + (ppart(3,j,k)**2 + ppart(4,j,k)**2 + ppart(5,j,k)**2) * dex**2 )
      dtc1 = dt / (gam - ppart(5,j,k) * dex)
      x = x + ppart(3,j,k) * dtc1
      y = y + ppart(4,j,k) * dtc1

      ! find particles going out of bounds
      mm = 0
      ! count how many particles are going in each direction in ncl
      ! save their address and destination in ihole
      ! use periodic boundary conditions and check for roundoff error
      ! mm = direction particle is going
      if (x >= edgerx) then
         if (x >= anx) x = x - anx
         mm = 2
      else if (x < edgelx) then
         if (x < 0.0) then
            x = x + anx
            if (x < anx) then
               mm = 1
            else
               x = 0.0
            endif
         else
            mm = 1
         endif
      endif
      if (y >= edgery) then
         if (y >= any) y = y - any
         mm = mm + 6
      else if (y < edgely) then
         if (y < 0.0) then
            y = y + any
            if (y < any) then
               mm = mm + 3
            else
               y = 0.0
            endif
         else
            mm = mm + 3
         endif
      endif
      ! set new position
      ppart(1,j,k) = x
      ppart(2,j,k) = y
      ! increment counters
      if (mm.gt.0) then
         ncl(mm,k) = ncl(mm,k) + 1
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(1,ih+1,k) = j
            ihole(2,ih+1,k) = mm
         else
            nh = 1
         endif
      endif

    enddo

    sum2 = sum2 + sum1

    ! set error and end of file flag
    ! ihole overflow
    if (nh.gt.0) then
       irc = ih
       ih = -ih
    endif
    ihole(1,1,k) = ih

  enddo
!$omp end parallel do

! normalize kinetic energy
  ek = ek + sum2

end subroutine PPGRBPPUSHF23L_RB

subroutine PPGRDCJPPOST2L_RB( ppart, fxy, bxy, psit, cu, dcu, amu, kpic, noff, &
  nyp, qbm, dt, ci, idimp, nppmx, nx, mx, my, nxv, nypmx, mx1, mxyp1, dex )

  implicit none

  integer,  intent(in) :: noff, nyp, idimp, nppmx, nx, mx, my
  integer, intent(in) :: nxv, nypmx, mx1, mxyp1
  real, intent(in) ::  qbm, dt, ci
  real, dimension(idimp, nppmx, mxyp1), intent(inout) :: ppart
  real, dimension(2, nxv, nypmx), intent(in) :: fxy
  real, dimension(3, nxv, nypmx), intent(in) :: bxy
  real, dimension(nxv, nypmx), intent(in) :: psit
  real, dimension(3, nxv, nypmx), intent(inout) :: cu
  real, dimension(2, nxv, nypmx), intent(inout) :: dcu
  real, dimension(3, nxv, nypmx), intent(inout) :: amu
  integer, dimension(mxyp1), intent(in) :: kpic
  real, intent(in) :: dex
  ! local variables
  real :: idex, dex2, qtmh, qtmh1, qtmh2, ex, ey, ez, wx, wy, wz, bx, by, bz, px, py, pz, x, y
  real :: dxp, dyp, amx, amy, tmp_x, tmp_y, tmp_z, dpx, dpy, dpz, pxx, pyy, pxy
  real :: ostq, p6, p7, dti, dx, dy, qm, inv_p7
  integer :: i, j, k, noffp, moffp, nppp, nn, mm, mnoff
  real, dimension(2, mx+1, my+1) :: sfxy, sdcu
  real, dimension(3, mx+1, my+1) :: sbxy, scu, samu
  real, dimension(mx+1, my+1) :: spsit

  qtmh = 0.5 * qbm * dt
  dti = 1.0 / dt
  idex = 1.0 / dex
  dex2 = dex * dex
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV)) return
! loop over tiles
!$omp parallel do &
!$omp& private(i, j, k, noffp, moffp, nppp, nn, mm, mnoff) &
!$omp& private(qtmh1, qtmh2, ex, ey, ez, wx, wy, wz, bx, by, bz, px, py, pz, x, y) &
!$omp& private(dxp, dyp, amx, amy, tmp_x, tmp_y, tmp_z, dpx, dpy, dpz, pxx, pyy, pxy) &
!$omp& private(ostq, p6, p7, dx, dy, qm, inv_p7) &
!$omp& private(sfxy, sdcu, sbxy, scu, samu, spsit)
  do k = 1, mxyp1
    noffp = (k - 1) / mx1
    moffp = my * noffp
    noffp = mx * (k - mx1 * noffp - 1)
    nppp = kpic(k)
    mnoff = moffp + noff - 1

    ! load local fields from global arrays
    nn = min(mx,nx-noffp) + 1
    mm = min(my,nyp-moffp) + 1
    do j = 1, mm
      do i = 1, nn
        sfxy(1,i,j) = fxy(1, i + noffp, j + moffp)
        sfxy(2,i,j) = fxy(2, i + noffp, j + moffp)
      enddo
    enddo
    do j = 1, mm
      do i = 1, nn
        sbxy(1,i,j) = bxy(1, i + noffp, j + moffp)
        sbxy(2,i,j) = bxy(2, i + noffp, j + moffp)
        sbxy(3,i,j) = bxy(3, i + noffp, j + moffp)
      enddo
    enddo
    do j = 1, mm
      do i = 1, nn
        spsit(i,j) = psit(i + noffp, j + moffp)
      enddo
    enddo

    ! zero out local accumulators
    do j = 1, my + 1
      do i = 1, mx + 1
        samu(1,i,j) = 0.0
        samu(2,i,j) = 0.0
        samu(3,i,j) = 0.0
        sdcu(1,i,j) = 0.0
        sdcu(2,i,j) = 0.0
        scu(1,i,j) = 0.0
        scu(2,i,j) = 0.0
        scu(3,i,j) = 0.0
      enddo
    enddo

    ! loop over particles in tile
    do j = 1, nppp

      ! find interpolation weights
      x = ppart(1,j,k)
      y = ppart(2,j,k)
      ! p6 = ppart(6,j,k)
      ! p7 = ppart(7,j,k)
      p6 = sqrt(1.0 + (ppart(3,j,k)**2 + ppart(4,j,k)**2 + ppart(5,j,k)**2) * dex2 )
      p7 = p6 - ppart(5,j,k) * dex
      qm = ppart(8,j,k)
      nn = x
      mm = y
      dxp = x - real(nn)
      dyp = y - real(mm)
      nn = nn - noffp + 1
      mm = mm - mnoff
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      
      inv_p7 = 1.0 / p7
      qtmh1 = qtmh * inv_p7
      qtmh2 = qtmh1 * p6
      
      ! find wakefield
      wx = amx * sfxy(1,nn,mm)
      wy = amx * sfxy(2,nn,mm)
      wz = amx * spsit(nn,mm)

      wx = amy * (dxp * sfxy(1,nn+1,mm) + wx)
      wy = amy * (dxp * sfxy(2,nn+1,mm) + wy)
      wz = amy * ( dxp * spsit(nn+1,mm) + wz)

      tmp_x = amx * sfxy(1,nn,mm+1)
      tmp_y = amx * sfxy(2,nn,mm+1)
      tmp_z = amx * spsit(nn,mm+1)

      wx = wx + dyp * (dxp * sfxy(1,nn+1,mm+1) + tmp_x) 
      wy = wy + dyp * (dxp * sfxy(2,nn+1,mm+1) + tmp_y)
      wz = wz + dyp * ( dxp * spsit(nn+1,mm+1) + tmp_z)

      ! find magnetic field
      bx = amx * sbxy(1,nn,mm)
      by = amx * sbxy(2,nn,mm)
      bz = amx * sbxy(3,nn,mm)

      bx = amy * (dxp * sbxy(1,nn+1,mm) + bx)
      by = amy * (dxp * sbxy(2,nn+1,mm) + by)
      bz = amy * (dxp * sbxy(3,nn+1,mm) + bz)

      tmp_x = amx * sbxy(1,nn,mm+1)
      tmp_y = amx * sbxy(2,nn,mm+1)
      tmp_z = amx * sbxy(3,nn,mm+1)

      bx = bx + dyp * (dxp * sbxy(1,nn+1,mm+1) + tmp_x)
      by = by + dyp * (dxp * sbxy(2,nn+1,mm+1) + tmp_y)
      bz = bz + dyp * (dxp * sbxy(3,nn+1,mm+1) + tmp_z)

      ! NEW
      ! obtain E-field (in the unit of dex)
      ex = by - wx
      ey =-bx - wy
      ez = wz * idex

      ! convert B-field to normalized unit
      bx = bx * dex
      by = by * dex

      ! half acceleration
      ex = qtmh2 * ex
      ey = qtmh2 * ey
      ez = qtmh2 * ez
      px = ppart(3,j,k) + ex
      py = ppart(4,j,k) + ey
      pz = ppart(5,j,k) + ez

      ! rotate about magnetic field
      bx = qtmh1 * bx
      by = qtmh1 * by
      bz = qtmh1 * bz

      tmp_x = px + py * bz - pz * by
      tmp_y = py + pz * bx - px * bz
      tmp_z = pz + px * by - py * bx

      ostq = 2.0 / ( 1.0 + bx**2 + by**2 + bz**2 )
      bx = bx * ostq
      by = by * ostq
      bz = bz * ostq

      px = px + tmp_y * bz - tmp_z * by
      py = py + tmp_z * bx - tmp_x * bz
      pz = pz + tmp_x * by - tmp_y * bx

      ! half acceleration
      px = px + ex
      py = py + ey
      pz = pz + ez

      dpx = (px - ppart(3,j,k)) * dti
      dpy = (py - ppart(4,j,k)) * dti

      ! time-centered momentum
      px = 0.5 * (px + ppart(3,j,k))
      py = 0.5 * (py + ppart(4,j,k))
      pz = 0.5 * (pz + ppart(5,j,k))

      ! time-centered gamma and gamma - pz
      ppart(6,j,k) = sqrt(1.0 + (px**2 + py**2 + pz**2) * dex2 )
      ppart(7,j,k) = ppart(6,j,k) - pz * dex
      inv_p7 = 1.0 / ppart(7,j,k)

      ! time-centered Ez (normalized unit)
      wz = qbm * ( wz + ( wx * px + wy * py ) * dex2 * inv_p7 )

      ! time-centered values of the time derivative of p
      dpx = dpx + px * wz * inv_p7
      dpy = dpy + py * wz * inv_p7

      pxx = px * px * inv_p7
      pyy = py * py * inv_p7
      pxy = px * py * inv_p7

      ! deposit momentum flux, acceleration density, and current density
      qm = qm * inv_p7
      amx = qm * amx
      dxp = qm * dxp

      dx = amx * amy
      dy = dxp * amy

      samu(1,nn,mm) = samu(1,nn,mm) + pxx * dx
      samu(2,nn,mm) = samu(2,nn,mm) + pyy * dx
      samu(3,nn,mm) = samu(3,nn,mm) + pxy * dx

      sdcu(1,nn,mm) = sdcu(1,nn,mm) + dpx * dx
      sdcu(2,nn,mm) = sdcu(2,nn,mm) + dpy * dx

      scu(1,nn,mm) = scu(1,nn,mm) + px * dx
      scu(2,nn,mm) = scu(2,nn,mm) + py * dx
      scu(3,nn,mm) = scu(3,nn,mm) + pz * dx

      dx = amx * dyp
      samu(1,nn+1,mm) = samu(1,nn+1,mm) + pxx * dy
      samu(2,nn+1,mm) = samu(2,nn+1,mm) + pyy * dy
      samu(3,nn+1,mm) = samu(3,nn+1,mm) + pxy * dy

      sdcu(1,nn+1,mm) = sdcu(1,nn+1,mm) + dpx * dy
      sdcu(2,nn+1,mm) = sdcu(2,nn+1,mm) + dpy * dy

      scu(1,nn+1,mm) = scu(1,nn+1,mm) + px * dy
      scu(2,nn+1,mm) = scu(2,nn+1,mm) + py * dy
      scu(3,nn+1,mm) = scu(3,nn+1,mm) + pz * dy

      dy = dxp * dyp
      samu(1,nn,mm+1) = samu(1,nn,mm+1) + pxx * dx
      samu(2,nn,mm+1) = samu(2,nn,mm+1) + pyy * dx
      samu(3,nn,mm+1) = samu(3,nn,mm+1) + pxy * dx

      sdcu(1,nn,mm+1) = sdcu(1,nn,mm+1) + dpx * dx
      sdcu(2,nn,mm+1) = sdcu(2,nn,mm+1) + dpy * dx

      scu(1,nn,mm+1) = scu(1,nn,mm+1) + px * dx
      scu(2,nn,mm+1) = scu(2,nn,mm+1) + py * dx
      scu(3,nn,mm+1) = scu(3,nn,mm+1) + pz * dx

      samu(1,nn+1,mm+1) = samu(1,nn+1,mm+1) + pxx * dy
      samu(2,nn+1,mm+1) = samu(2,nn+1,mm+1) + pyy * dy
      samu(3,nn+1,mm+1) = samu(3,nn+1,mm+1) + pxy * dy

      sdcu(1,nn+1,mm+1) = sdcu(1,nn+1,mm+1) + dpx * dy
      sdcu(2,nn+1,mm+1) = sdcu(2,nn+1,mm+1) + dpy * dy

      scu(1,nn+1,mm+1) = scu(1,nn+1,mm+1) + px * dy
      scu(2,nn+1,mm+1) = scu(2,nn+1,mm+1) + py * dy
      scu(3,nn+1,mm+1) = scu(3,nn+1,mm+1) + pz * dy
    enddo

    ! deposit currents to interior points in global array
    nn = min(mx, nxv-noffp)
    mm = min(my, nypmx-moffp)

    do j = 2, mm
      do i = 2, nn
        amu(1, i + noffp, j + moffp) = amu(1, i + noffp, j + moffp) + samu(1,i,j)
        amu(2, i + noffp, j + moffp) = amu(2, i + noffp, j + moffp) + samu(2,i,j)
        amu(3, i + noffp, j + moffp) = amu(3, i + noffp, j + moffp) + samu(3,i,j)

        dcu(1, i + noffp, j + moffp) = dcu(1, i + noffp, j + moffp) + sdcu(1,i,j)
        dcu(2, i + noffp, j + moffp) = dcu(2, i + noffp, j + moffp) + sdcu(2,i,j)

        cu(1, i + noffp, j + moffp) = cu(1, i + noffp, j + moffp) + scu(1,i,j)
        cu(2, i + noffp, j + moffp) = cu(2, i + noffp, j + moffp) + scu(2,i,j)
        cu(3, i + noffp, j + moffp) = cu(3, i + noffp, j + moffp) + scu(3,i,j)
      enddo
    enddo

    ! deposit currents to edge points in global array
    mm = min(my+1,nypmx-moffp)
    do i = 2, nn
!$omp atomic
      amu(1,i+noffp,1+moffp) = amu(1,i+noffp,1+moffp) + samu(1,i,1)
!$omp atomic
      amu(2,i+noffp,1+moffp) = amu(2,i+noffp,1+moffp) + samu(2,i,1)
!$omp atomic
      amu(3,i+noffp,1+moffp) = amu(3,i+noffp,1+moffp) + samu(3,i,1)
!$omp atomic
      dcu(1,i+noffp,1+moffp) = dcu(1,i+noffp,1+moffp) + sdcu(1,i,1)
!$omp atomic
      dcu(2,i+noffp,1+moffp) = dcu(2,i+noffp,1+moffp) + sdcu(2,i,1)
!$omp atomic
      cu(1,i+noffp,1+moffp) = cu(1,i+noffp,1+moffp) + scu(1,i,1)
!$omp atomic
      cu(2,i+noffp,1+moffp) = cu(2,i+noffp,1+moffp) + scu(2,i,1)
!$omp atomic
      cu(3,i+noffp,1+moffp) = cu(3,i+noffp,1+moffp) + scu(3,i,1)
      if (mm > my) then
!$omp atomic
         amu(1,i+noffp,mm+moffp) = amu(1,i+noffp,mm+moffp) + samu(1,i,mm)
!$omp atomic
         amu(2,i+noffp,mm+moffp) = amu(2,i+noffp,mm+moffp) + samu(2,i,mm)
!$omp atomic
         amu(3,i+noffp,mm+moffp) = amu(3,i+noffp,mm+moffp) + samu(3,i,mm)
!$omp atomic
         dcu(1,i+noffp,mm+moffp) = dcu(1,i+noffp,mm+moffp) + sdcu(1,i,mm)
!$omp atomic
         dcu(2,i+noffp,mm+moffp) = dcu(2,i+noffp,mm+moffp) + sdcu(2,i,mm)
!$omp atomic
         cu(1,i+noffp,mm+moffp) = cu(1,i+noffp,mm+moffp) + scu(1,i,mm)
!$omp atomic
         cu(2,i+noffp,mm+moffp) = cu(2,i+noffp,mm+moffp) + scu(2,i,mm)
!$omp atomic
         cu(3,i+noffp,mm+moffp) = cu(3,i+noffp,mm+moffp) + scu(3,i,mm)
      endif
    enddo

    nn = min(mx+1,nxv-noffp)
    do j = 1, mm
!$omp atomic
      amu(1,1+noffp,j+moffp) = amu(1,1+noffp,j+moffp) + samu(1,1,j)
!$omp atomic
      amu(2,1+noffp,j+moffp) = amu(2,1+noffp,j+moffp) + samu(2,1,j)
!$omp atomic
      amu(3,1+noffp,j+moffp) = amu(3,1+noffp,j+moffp) + samu(3,1,j)
!$omp atomic
      dcu(1,1+noffp,j+moffp) = dcu(1,1+noffp,j+moffp) + sdcu(1,1,j)
!$omp atomic
      dcu(2,1+noffp,j+moffp) = dcu(2,1+noffp,j+moffp) + sdcu(2,1,j)
!$omp atomic
      cu(1,1+noffp,j+moffp) = cu(1,1+noffp,j+moffp) + scu(1,1,j)
!$omp atomic
      cu(2,1+noffp,j+moffp) = cu(2,1+noffp,j+moffp) + scu(2,1,j)
!$omp atomic
      cu(3,1+noffp,j+moffp) = cu(3,1+noffp,j+moffp) + scu(3,1,j)
      if (nn > mx) then
!$omp atomic
         amu(1,nn+noffp,j+moffp) = amu(1,nn+noffp,j+moffp) + samu(1,nn,j)
!$omp atomic
         amu(2,nn+noffp,j+moffp) = amu(2,nn+noffp,j+moffp) + samu(2,nn,j)
!$omp atomic
         amu(3,nn+noffp,j+moffp) = amu(3,nn+noffp,j+moffp) + samu(3,nn,j)
!$omp atomic
         dcu(1,nn+noffp,j+moffp) = dcu(1,nn+noffp,j+moffp) + sdcu(1,nn,j)
!$omp atomic
         dcu(2,nn+noffp,j+moffp) = dcu(2,nn+noffp,j+moffp) + sdcu(2,nn,j)
!$omp atomic
         cu(1,nn+noffp,j+moffp) = cu(1,nn+noffp,j+moffp) + scu(1,nn,j)
!$omp atomic
         cu(2,nn+noffp,j+moffp) = cu(2,nn+noffp,j+moffp) + scu(2,nn,j)
!$omp atomic
         cu(3,nn+noffp,j+moffp) = cu(3,nn+noffp,j+moffp) + scu(3,nn,j)
      endif
    enddo
  enddo
!$omp end parallel do

end subroutine PPGRDCJPPOST2L_RB

      end module part2d_lib
