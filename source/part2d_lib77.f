!-----------------------------------------------------------------------
      subroutine PPDBLKP2L(part,kpic,npp,noff,nppmx,idimp,npmax,mx,my,  &
     &mx1,mxyp1,irc)
! this subroutine finds the maximum number of particles in each tile of
! mx, my to calculate size of segmented particle array ppart
! linear interpolation, spatial decomposition in y direction
! input: all except kpic, nppmx, output: kpic, nppmx
! part = input particle array
! part(1,n) = position x of particle n in partition
! part(2,n) = position y of particle n in partition
! kpic = output number of particles per tile
! nppmx = return maximum number of particles in tile
! npp = number of particles in partition
! noff = backmost global gridpoint in particle partition
! idimp = size of phase space = 4
! npmax = maximum number of particles in each partition
! mx/my = number of grids in sorting cell in x and y
! mx1 = (system length in x direction - 1)/mx + 1
! mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
! irc = maximum overflow, returned only if error occurs, when irc > 0
      implicit none
      integer nppmx, idimp, npmax, mx, my, mx1, mxyp1, irc
      integer kpic, npp, noff
      real part
      dimension part(idimp,npmax)
      dimension kpic(mxyp1)
! local data
      integer j, k, n, m, mnoff, isum, ist, npx, ierr
      mnoff = noff
      ierr = 0
! clear counter array
      do 10 k = 1, mxyp1
      kpic(k) = 0
   10 continue
! find how many particles in each tile
      do 20 j = 1, npp
      n = part(1,j)
      n = n/mx + 1
      m = part(2,j)
      m = (m - mnoff)/my
      m = n + mx1*m
      if (m.le.mxyp1) then
         kpic(m) = kpic(m) + 1
      else
         ierr = max(ierr,m-mxyp1)
      endif
   20 continue
! find maximum
      isum = 0
      npx = 0
      do 30 k = 1, mxyp1
      ist = kpic(k)
      npx = max(npx,ist)
      isum = isum + ist
   30 continue
      nppmx = npx
! check for errors
      if (ierr.gt.0) then
         irc = ierr
      else if (isum.ne.npp) then
         irc = isum
      endif
      return
      end
!-----------------------------------------------------------------------
      subroutine PPPMOVIN2L(part,ppart,kpic,npp,noff,nppmx,idimp,npmax, &
     &mx,my,mx1,mxyp1,irc)
! this subroutine sorts particles by x,y grid in tiles of
! mx, my and copies to segmented array ppart
! linear interpolation, spatial decomposition in y direction
! input: all except ppart, kpic, output: ppart, kpic
! part/ppart = input/output particle arrays
! part(1,n) = position x of particle n in partition
! part(2,n) = position y of particle n in partition
! kpic = output number of particles per tile
! nppmx = maximum number of particles in tile
! npp = number of particles in partition
! noff = backmost global gridpoint in particle partition
! idimp = size of phase space = 4
! npmax = maximum number of particles in each partition
! mx/my = number of grids in sorting cell in x and y
! mx1 = (system length in x direction - 1)/mx + 1
! mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
! irc = maximum overflow, returned only if error occurs, when irc > 0
      implicit none
      integer nppmx, idimp, npmax, mx, my, mx1, mxyp1, irc
      integer kpic, npp, noff
      real part, ppart
      dimension part(idimp,npmax), ppart(nppmx,idimp,mxyp1)
      dimension kpic(mxyp1)
! local data
      integer i, j, k, n, m, mnoff, ip, ierr
      mnoff = noff
      ierr = 0
! clear counter array
      do 10 k = 1, mxyp1
      kpic(k) = 0
   10 continue
! find addresses of particles at each tile and reorder particles
      do 30 j = 1, npp
      n = part(1,j)
      n = n/mx + 1
      m = part(2,j)
      m = (m - mnoff)/my
      m = n + mx1*m
      ip = kpic(m) + 1
      if (ip.le.nppmx) then
         do 20 i = 1, idimp
         ppart(ip,i,m) = part(i,j)
   20    continue
      else
         ierr = max(ierr,ip-nppmx)
      endif
      kpic(m) = ip
   30 continue
      if (ierr.gt.0) irc = ierr
      return
      end
!-----------------------------------------------------------------------
      subroutine PPPCHECK2L(ppart,kpic,noff,nyp,idimp,nppmx,nx,mx,my,mx1&
     &,myp1,irc)
! this subroutine performs a sanity check to make sure particles sorted
! by x,y grid in tiles of mx, my, are all within bounds.
! tiles are assumed to be arranged in 2D linear memory
! input: all except irc
! output: irc
! ppart(1,n,k) = position x of particle n in tile k
! ppart(2,n,k) = position y of particle n in tile k
! kpic(k) = number of reordered output particles in tile k
! noff = lowermost global gridpoint in particle partition.
! nyp = number of primary (complete) gridpoints in particle partition
! idimp = size of phase space = 4
! nppmx = maximum number of particles in tile
! nx = system length in x direction
! mx/my = number of grids in sorting cell in x/y
! mx1 = (system length in x direction - 1)/mx + 1
! myp1 = (partition length in y direction - 1)/my + 1
! irc = particle error, returned only if error occurs, when irc > 0
      implicit none
      integer noff, nyp, idimp, nppmx, nx, mx, my, mx1, myp1, irc
      real ppart
      integer kpic
      dimension ppart(nppmx,idimp,mx1*myp1)
      dimension kpic(mx1*myp1)
! local data
      integer mxyp1, noffp, moffp, nppp, j, k, ist, nn, mm
      real edgelx, edgely, edgerx, edgery, dx, dy
      mxyp1 = mx1*myp1
! loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(j,k,noffp,moffp,nppp,nn,mm,ist,edgelx,edgely,edgerx,
!$OMP& edgery,dx,dy)
      do 20 k = 1, mxyp1
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(k)
      nn = min(mx,nx-noffp)
      mm = min(my,nyp-moffp)
      edgelx = noffp
      edgerx = noffp + nn
      edgely = noff + moffp
      edgery = noff + moffp + mm
! loop over particles in tile
      do 10 j = 1, nppp
      dx = ppart(j,1,k)
      dy = ppart(j,2,k)
! find particles going out of bounds
      ist = 0
      if (dx.lt.edgelx) ist = 1
      if (dx.ge.edgerx) ist = 2
      if (dy.lt.edgely) ist = ist + 3
      if (dy.ge.edgery) ist = ist + 6
      if (ist.gt.0) irc = k
   10 continue
   20 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPGPPOST2L(ppart,q,kpic,noff,qm,idimp,nppmx,mx,my,nxv, &
     &nypmx,mx1,mxyp1)
! vecterization
! transposed
! for 2d code, this subroutine calculates particle charge density
! using first-order linear interpolation, periodic boundaries
! OpenMP version using guard cells, for distributed data
! data deposited in tiles
! particles stored segmented array
! 17 flops/particle, 6 loads, 4 stores
! input: all, output: q
! charge density is approximated by values at the nearest grid points
! q(n,m)=qm*(1.-dx)*(1.-dy)
! q(n+1,m)=qm*dx*(1.-dy)
! q(n,m+1)=qm*(1.-dx)*dy
! q(n+1,m+1)=qm*dx*dy
! where n,m = leftmost grid points and dx = x-n, dy = y-m
! ppart(1,n,m) = position x of particle n in partition in tile m
! ppart(2,n,m) = position y of particle n in partition in tile m
! q(j,k) = charge density at grid point (j,kk),
! where kk = k + noff - 1
! kpic = number of particles per tile
! noff = lowermost global gridpoint in particle partition.
! qm = charge on particle, in units of e
! idimp = size of phase space = 4
! nppmx = maximum number of particles in tile
! mx/my = number of grids in sorting cell in x/y
! nxv = first dimension of charge array, must be >= nx+1
! nypmx = maximum size of particle partition, including guard cells.
! mx1 = (system length in x direction - 1)/mx + 1
! mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
      implicit none
      integer noff, idimp, nppmx, mx, my, nxv, nypmx, mx1, mxyp1
      real qm
      real ppart, q
      integer kpic
      dimension ppart(nppmx,idimp,mxyp1), q(nxv*nypmx), kpic(mxyp1)
! local data
!      integer MXV, MYV
!      parameter(MXV=33,MYV=33)
      integer npblk, lvect
      parameter(npblk=32,lvect=4)
      integer ipp, joff, nps, m, lxv
      integer noffp, moffp, nppp
      integer mnoff, i, j, k, nn, mm
      real x, y, dxp, dyp, amx, amy
      real sq
!     dimension sq(MXV,MYV)
      dimension sq((mx+1)*(my+1))
c scratch arrays
      integer n
      real s
      dimension n(npblk), s(npblk,lvect)
      lxv = mx + 1      
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV)) return
! loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,noffp,moffp,nppp,mnoff,nn,mm,x,y,dxp,dyp,amx,amy,
!$OMP& ipp,joff,nps,n,m,s,
!$OMP& sq)
      do 110 k = 1, mxyp1
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(k)
      mnoff = moffp + noff
! zero out local accumulator
      do 10 j = 1, (mx+1)*(my+1)
      sq(j) = 0.0
   10 continue
c loop over particles in tile
      ipp = nppp/npblk
c outer loop over number of full blocks
      do 50 m = 1, ipp
      joff = npblk*(m - 1)
c inner loop over particles in block
      do 20 j = 1, npblk
! find interpolation weights
      x = ppart(j+joff,1,k)
      y = ppart(j+joff,2,k)
      nn = x
      mm = y
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
      n(j) = nn - noffp + lxv*(mm - mnoff)
      amx = qm - dxp
      amy = 1.0 - dyp
      s(j,1) = amx*amy
      s(j,2) = dxp*amy
      s(j,3) = amx*dyp
      s(j,4) = dxp*dyp
   20 continue
c deposit charge within tile to local accumulator
      do 40 j = 1, npblk
      nn = n(j)
      mm = nn + lxv - 2
!dir$ ivdep
      do 30 i = 1, lvect
      if (i.gt.2) nn = mm
      sq(i+nn) = sq(i+nn) + s(j,i)
   30 continue
   40 continue
   50 continue
      nps = npblk*ipp + 1
c loop over remaining particles
      do 60 j = nps, nppp
c find interpolation weights
! find interpolation weights
      x = ppart(j,1,k)
      y = ppart(j,2,k)
      nn = x
      mm = y
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
      nn = nn - noffp + 1 + lxv*(mm - mnoff)
      amx = qm - dxp
      amy = 1.0 - dyp
! deposit charge within tile to local accumulator
      x = sq(nn) + amx*amy
      y = sq(nn+1) + dxp*amy
      sq(nn) = x
      sq(nn+1) = y
      x = sq(nn+lxv) + amx*dyp
      y = sq(nn+1+lxv) + dxp*dyp
      sq(nn+lxv) = x
      sq(nn+1+lxv) = y
   60 continue
! deposit charge to interior points in global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      do 80 j = 2, mm
      do 70 i = 2, nn
      q(i+noffp+nxv*(j+moffp-1)) = q(i+noffp+nxv*(j+moffp-1)) + 
     1sq(i+lxv*(j-1))
   70 continue
   80 continue
! deposit charge to edge points in global array
      mm = min(my+1,nypmx-moffp)
      do 90 i = 2, nn
!$OMP ATOMIC
      q(i+noffp+nxv*moffp) = q(i+noffp+nxv*moffp) + sq(i)
      if (mm > my) then
!$OMP ATOMIC
         q(i+noffp+nxv*(mm+moffp-1)) = q(i+noffp+nxv*(mm+moffp-1)) + 
     1sq(i+lxv*(mm-1))
      endif
   90 continue
      nn = min(mx+1,nxv-noffp)
      do 100 j = 1, mm
!$OMP ATOMIC
      q(1+noffp+nxv*(j+moffp-1)) = q(1+noffp+nxv*(j+moffp-1)) + 
     1sq(1+lxv*(j-1))
      if (nn > mx) then
!$OMP ATOMIC
         q(nn+noffp+nxv*(j+moffp-1)) = q(nn+noffp+nxv*(j+moffp-1)) + 
     1sq(nn+lxv*(j-1))
      endif
  100 continue
  110 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPGRDCJPPOST2L_QP(ppart,fxy,bxy,psit,cu,dcu,amu,kpic,no
     1ff,nyp,qm,qbm,dt,ci,idimp,nppmx,nx,mx,my,nxv,nypmx,mx1,mxyp1,dex)
! vecterization
! transposed
! for 2-1/2d code, this subroutine calculates particle momentum flux,
! acceleration density and current density using first-order spline
! interpolation for relativistic particles.
! OpenMP version using guard cells, for distributed data
! data read/written in tiles
! particles stored segmented array
! 241 flops/particle, 2 divide, 1 sqrt, 69 loads, 40 stores
! input: all, output: cu, dcu, amu
! current density is approximated by values at the nearest grid points
! cu(i,n,m)=qci*(1.-dx)*(1.-dy)
! cu(i,n+1,m)=qci*dx*(1.-dy)
! cu(i,n,m+1)=qci*(1.-dx)*dy
! cu(i,n+1,m+1)=qci*dx*dy
! and qci = qm*pj*gami, where j = x,y,z, for i = 1, 3
! where pj = .5*(pj(t+dt/2)+pj(t-dt/2))
! where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
! acceleration density is approximated by values at the nearest grid
! points
! dcu(i,n,m)=qci*(1.-dx)*(1.-dy)
! dcu(i,n+1,m)=qci*dx*(1.-dy)
! dcu(i,n,m+1)=qci*(1.-dx)*dy
! dcu(i,n+1,m+1)=qci*dx*dy
! and qci = qm*dvj*gami/dt, where j = x,y,z, for i = 1, 3
! where dvj = dpj - pj*gami*dgamma, dpj = (pj(t+dt/2)-pj(t-dt/2)), 
! pj = .5*(pj(t+dt/2)+pj(t-dt/2)),
! dgamma = (q/m)*ci*ci*gami*(sum(pj*Ej))*dt,
! and Ej = jth component of electric field
! momentum flux is approximated by values at the nearest grid points
! amu(i,n,m)=qci*(1.-dx)*(1.-dy)
! amu(i,n+1,m)=qci*dx*(1.-dy)
! amu(i,n,m+1)=qci*(1.-dx)*dy
! amu(i,n+1,m+1)=qci*dx*dy
! and qci = qm*pj*pk*gami**2, where jk = xx-yy,xy,zx,zy, for i = 1, 4
! where pj = 0.5*(pj(t+dt/2)+pj(t-dt/2),
! and pk = 0.5*(pk(t+dt/2)+pk(t-dt/2))
! where n,m = nearest grid points and dx = x-n, dy = y-m
! momentum equations at t=t+dt/2 are calculated from:
! px(t+dt/2) = rot(1)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
!    rot(2)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
!    rot(3)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
!    .5*(q/m)*fx(x(t),y(t))*dt)
! py(t+dt/2) = rot(4)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
!    rot(5)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
!    rot(6)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
!    .5*(q/m)*fy(x(t),y(t))*dt)
! pz(t+dt/2) = rot(7)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
!    rot(8)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
!    rot(9)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
!    .5*(q/m)*fz(x(t),y(t))*dt)
! where q/m is charge/mass, and the rotation matrix is given by:
!    rot(1) = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
!    rot(2) = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
!    rot(3) = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
!    rot(4) = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
!    rot(5) = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
!    rot(6) = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
!    rot(7) = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
!    rot(8) = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
!    rot(9) = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
! and om**2 = omx**2 + omy**2 + omz**2
! the rotation matrix is determined by:
! omx = (q/m)*bx(x(t),y(t))*gami, omy = (q/m)*by(x(t),y(t))*gami, and
! omz = (q/m)*bz(x(t),y(t))*gami.
! fx(x(t),y(t)), fy(x(t),y(t)), and fz(x(t),y(t))
! bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
! are approximated by interpolation from the nearest grid points:
! fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
!    + dx*fx(n+1,m+1))
! where n,m = leftmost grid points and dx = x-n, dy = y-m
! similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
! ppart(1,n,m) = position x of particle n in partition in tile m at t
! ppart(2,n,m) = position y of particle n in partition in tile m at t
! ppart(3,n,m) = momentum px of particle n in partition in tile m
! at t - dt/2
! ppart(4,n,m) = momentum py of particle n in partition in tile m
! at t - dt/2
! ppart(5,n,m) = momentum pz of particle n in partition in tile m
! at t - dt/2
! fxy(1,j,k) = x component of force/charge at grid (j,kk)
! fxy(2,j,k) = y component of force/charge at grid (j,kk)
! fxy(3,j,k) = z component of force/charge at grid (j,kk)
! that is, convolution of electric field over particle shape,
! where kk = k + noff - 1
! bxy(1,j,k) = x component of magnetic field at grid (j,kk)
! bxy(2,j,k) = y component of magnetic field at grid (j,kk)
! bxy(3,j,k) = z component of magnetic field at grid (j,kk)
! that is, the convolution of magnetic field over particle shape,
! where kk = k + noff - 1
! cu(i,j,k) = ith component of current density
! at grid point j,kk for i = 1, 3
! dcu(i,j,k) = ith component of acceleration density
! at grid point j,kk for i = 1, 3
! amu(i,j,k) = ith component of momentum flux
! at grid point j,kk for i = 1, 4
! where kk = k + noff - 1
! kpic = number of particles per tile
! noff = lowermost global gridpoint in particle partition.
! nyp = number of primary (complete) gridpoints in particle partition
! qm = charge on particle, in units of e
! qbm = particle charge/mass ratio
! dt = time interval between successive calculations
! ci = reciprical of velocity of light
! idimp = size of phase space = 5
! nppmx = maximum number of particles in tile
! nx = system length in x direction
! mx/my = number of grids in sorting cell in x/y
! nxv = second dimension of field arrays, must be >= nx+1
! nypmx = maximum size of particle partition, including guard cells.
! mx1 = (system length in x direction - 1)/mx + 1
! mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
      implicit none
      integer noff, nyp, idimp, nppmx, nx, mx, my, nxv, nypmx
      integer mx1, mxyp1
      real qm, qbm, dt, ci,dex
      real ppart, fxy, bxy, cu, dcu, amu, psit
      integer kpic
      dimension ppart(nppmx,idimp,mxyp1)
      dimension fxy(2,nxv*nypmx), bxy(3,nxv*nypmx), psit(nxv*nypmx)
      dimension cu(3,nxv*nypmx), dcu(2,nxv*nypmx), amu(3,nxv*nypmx)
      dimension kpic(mxyp1)
! local data
!      integer MXV, MYV
!      parameter(MXV=33,MYV=33)
      integer npblk, lvect
      parameter(npblk=32,lvect=4)
      integer ipp, joff, nps, m, lxv
      integer noffp, moffp, nppp
      integer mnoff, i, j, k, nn, mm
      real qtmh, dti, ci2, gami, qtmg, gh, dxp, dyp, amx, amy
      real dx, dy, dz, ox, oy, oz
      real acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real x, y, vx, vy, vz, p2, v1, v2, v3, v4
      real sfxy, sbxy, scu, sdcu, samu
      dimension sfxy(3,(mx+1)*(my+1)), sbxy(3,(mx+1)*(my+1))
      dimension scu(3,(mx+1)*(my+1)), sdcu(2,(mx+1)*(my+1))
      dimension samu(3,(mx+1)*(my+1))
!      dimension sfxy(2,MXV,MYV), sbxy(3,MXV,MYV)
!      dimension spsit(MXV,MYV)
!      dimension scu(3,MXV,MYV), sdcu(2,MXV,MYV), samu(3,MXV,MYV)
      real qm1, qtmh1,qtmh2,idex,inv_part_7,ddx,ddy,p6,p7
!     dimension sfxy(3,mx+1,my+1), sbxy(3,mx+1,my+1)
!     dimension scu(3,mx+1,my+1,), sdcu(3,mx+1,my+1), samu(4,mx+1,my+1)
c scratch arrays
      integer n
      real s1, s2, s3, t
      dimension n(npblk), s1(npblk,lvect), s2(npblk,lvect), t(npblk,4)
      dimension s3(npblk,lvect)
      lxv = mx + 1
      qtmh = 0.5*qbm*dt
      dti = 1.0/dt
      ci2 = ci*ci
      idex = 1.0/dex
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV)) return
! loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,noffp,moffp,nppp,nn,mm,mnoff,x,y,vx,vy,v1,v2,v3,
!$OMP& dxp,dyp,amx,amy,dx,dy,dz,ox,oy,oz,acx,acy,acz,omzt, 
!$OMP& omt,anorm,rot1,rot2,  
!$OMP& sfxy,sbxy,scu,sdcu,samu,qm1,qtmh1,qtmh2,inv_part_7,ddx,
!$OMP& m,n,ipp,joff,nps,s1,s2,s3,t,
!$OMP& ddy,p6,p7)
      do 120 k = 1, mxyp1
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(k)
      mnoff = moffp + noff
! load local fields from global arrays
      nn = min(mx,nx-noffp) + 1
      mm = min(my,nyp-moffp) + 1
      do 20 j = 1, mm
      do 10 i = 1, nn
      sfxy(1,i+lxv*(j-1)) = fxy(1,i+noffp+nxv*(j+moffp-1))
      sfxy(2,i+lxv*(j-1)) = fxy(2,i+noffp+nxv*(j+moffp-1))
      sfxy(3,i+lxv*(j-1)) = psit(i+noffp+nxv*(j+moffp-1))
   10 continue
   20 continue
      do 40 j = 1, mm
      do 30 i = 1, nn
      sbxy(1,i+lxv*(j-1)) = bxy(1,i+noffp+nxv*(j+moffp-1))
      sbxy(2,i+lxv*(j-1)) = bxy(2,i+noffp+nxv*(j+moffp-1))
      sbxy(3,i+lxv*(j-1)) = bxy(3,i+noffp+nxv*(j+moffp-1))
   30 continue
   40 continue
! zero out local accumulators
      nn = lxv*(my + 1)
      do 50 i = 1, nn
      samu(1,i) = 0.0
      samu(2,i) = 0.0
      samu(3,i) = 0.0
      sdcu(1,i) = 0.0
      sdcu(2,i) = 0.0
      scu(1,i) = 0.0
      scu(2,i) = 0.0
      scu(3,i) = 0.0
   50 continue
      ipp = nppp/npblk
c outer loop over number of full blocks
      do 110 m = 1, ipp
      joff = npblk*(m - 1)
c inner loop over particles in block
      do 60 j = 1, npblk
! find interpolation weights
      x = ppart(j+joff,1,k)
      y = ppart(j+joff,2,k)
      p6 = ppart(j+joff,6,k)
      p7 = ppart(j+joff,7,k)
      nn = x
      mm = y
      dxp = x - real(nn)
      dyp = y - real(mm)
      n(j) = nn - noffp + lxv*(mm - mnoff)
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      
      inv_part_7 = 1.0 / p7
      qtmh1 = qtmh * inv_part_7
      qtmh2 = qtmh1 * p6
      s2(j,1) = amx*amy
      s2(j,2) = dxp*amy
      s2(j,3) = amx*dyp
      s2(j,4) = dxp*dyp
      t(j,1) = inv_part_7
      t(j,2) = qtmh1
      t(j,3) = qtmh2
      t(j,4) = p7
      qm1 = qm*inv_part_7
      amx = qm1*amx
      dxp = qm1*dxp
      s1(j,1) = amx*amy
      s1(j,2) = dxp*amy
      s1(j,3) = amx*dyp
      s1(j,4) = dxp*dyp      
   60 continue
      do 70 j = 1, npblk
      nn = n(j)
      mm = nn + lxv - 2
      dx = 0.0
      dy = 0.0
      dz = 0.0
      ox = 0.0
      oy = 0.0
      oz = 0.0
!dir$ ivdep
      do 65 i = 1, lvect
      if (i.gt.2) nn = mm
      dx = dx + sfxy(1,i+nn)*s2(j,i)
      dy = dy + sfxy(2,i+nn)*s2(j,i)
      dz = dz + sfxy(3,i+nn)*s2(j,i)
      ox = ox + sbxy(1,i+nn)*s2(j,i)
      oy = oy + sbxy(2,i+nn)*s2(j,i)
      oz = oz + sbxy(3,i+nn)*s2(j,i)
   65 continue
      s2(j,1) = dx
      s2(j,2) = dy
      s2(j,3) = dz
      s3(j,1) = ox
      s3(j,2) = oy
      s3(j,3) = oz
   70 continue
      do 80 j = 1, npblk
      dx = s2(j,1)
      dy = s2(j,2)
      dz = s2(j,3)
      ox = s3(j,1)
      oy = s3(j,2)
      oz = s3(j,3)
      inv_part_7 = t(j,1)
      qtmh1 = t(j,2)
      qtmh2 = t(j,3)
      p7 = t(j,4)
! calculate half impulse
      ddx = (-1.0)*qtmh2*dx + qtmh*oy
      ddy = (-1.0)*qtmh2*dy - qtmh*ox
! half acceleration
      vx = ppart(j+joff,3,k)
      vy = ppart(j+joff,4,k)
      acx = vx + ddx
      acy = vy + ddy
! find inverse gamma
! renormalize magnetic field
! calculate cyclotron frequency
      omzt = qtmh1*oz
! calculate rotation matrix
      omt = omzt*omzt
      anorm = 2.0/(1.0 + omt)
      rot1 = 0.5*(1.0 - omt)
      rot2 = omzt
! new momentum
      v1 = (rot1*acx + rot2*acy)*anorm + ddx
      v2 = (rot1*acy - rot2*acx)*anorm + ddy
      
! deposit momentum flux, acceleration density, and current density
      ox = 0.5*(v1 + vx)
      oy = 0.5*(v2 + vy)
      oz = 0.5*(1.0+dex*dex*(ox*ox+oy*oy))*inv_part_7 - 0.5*p7

      ppart(j+joff,6,k) = p7 + oz

      oz = oz*idex

      vx = (v1 - vx)*dti
      vy = (v2 - vy)*dti

      dz = qbm*(dz + (dx*ox+dy*oy)*dex*dex*inv_part_7)

      v1 = ox*ox*inv_part_7
      v2 = oy*oy*inv_part_7
      v3 = ox*oy*inv_part_7

      vx = vx+ox*dz*inv_part_7
      vy = vy+oy*dz*inv_part_7

      s2(j,1) = v1
      s2(j,2) = v2
      s2(j,3) = v3
      s2(j,4) = vx
      s3(j,1) = vy
      s3(j,2) = ox
      s3(j,3) = oy
      s3(j,4) = oz
   80 continue
      do 90 j = 1, npblk
      nn = n(j)
      mm = nn + lxv - 2
      v1 = s2(j,1)
      v2 = s2(j,2)
      v3 = s2(j,3)
      vx = s2(j,4)
      vy = s3(j,1)
      ox = s3(j,2)
      oy = s3(j,3)
      oz = s3(j,4)
!dir$ ivdep
      do 85 i = 1, lvect
      if (i.gt.2) nn = mm
      samu(1,i+nn) = samu(1,i+nn) + v1*s1(j,i)
      samu(2,i+nn) = samu(2,i+nn) + v2*s1(j,i)
      samu(3,i+nn) = samu(3,i+nn) + v3*s1(j,i)
      sdcu(1,i+nn) = sdcu(1,i+nn) + vx*s1(j,i)
      sdcu(2,i+nn) = sdcu(2,i+nn) + vy*s1(j,i)
      scu(1,i+nn) = scu(1,i+nn) + ox*s1(j,i)
      scu(2,i+nn) = scu(2,i+nn) + oy*s1(j,i)
      scu(3,i+nn) = scu(3,i+nn) + oz*s1(j,i)
   85 continue
   90 continue
  110 continue
      nps = npblk*ipp + 1
c loop over remaining particles
      do 111 j = nps, nppp
! find interpolation weights
      x = ppart(j,1,k)
      y = ppart(j,2,k)
      p6 = ppart(j,6,k)
      p7 = ppart(j,7,k)
      nn = x
      mm = y
      dxp = x - real(nn)
      dyp = y - real(mm)
      nn = nn - noffp + 1 + lxv*(mm - mnoff)
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      
      inv_part_7 = 1.0 / p7
      qtmh1 = qtmh * inv_part_7
      qtmh2 = qtmh1 * p6
      
! find electric field
      dx = amx*sfxy(1,nn)
      dy = amx*sfxy(2,nn)
      dz = amx*sfxy(3,nn)

      dx = amy*(dxp*sfxy(1,nn+1) + dx)
      dy = amy*(dxp*sfxy(2,nn+1) + dy)
      dz = amy*(dxp*sfxy(3,nn+1) + dz)

      acx = amx*sfxy(1,nn+lxv)
      acy = amx*sfxy(2,nn+lxv)
      acz = amx*sfxy(3,nn+lxv)

      dx = dx + dyp*(dxp*sfxy(1,nn+1+lxv) + acx) 
      dy = dy + dyp*(dxp*sfxy(2,nn+1+lxv) + acy)
      dz = dz + dyp*(dxp*sfxy(3,nn+1+lxv) + acz)

! find magnetic field
      ox = amx*sbxy(1,nn)
      oy = amx*sbxy(2,nn)
      oz = amx*sbxy(3,nn)

      ox = amy*(dxp*sbxy(1,nn+1) + ox)
      oy = amy*(dxp*sbxy(2,nn+1) + oy)
      oz = amy*(dxp*sbxy(3,nn+1) + oz)

      acx = amx*sbxy(1,nn+lxv)
      acy = amx*sbxy(2,nn+lxv)
      acz = amx*sbxy(3,nn+lxv)

      ox = ox + dyp*(dxp*sbxy(1,nn+1+lxv) + acx) 
      oy = oy + dyp*(dxp*sbxy(2,nn+1+lxv) + acy)
      oz = oz + dyp*(dxp*sbxy(3,nn+1+lxv) + acz)

! calculate half impulse
      ddx = (-1.0)*qtmh2*dx + qtmh*oy
      ddy = (-1.0)*qtmh2*dy - qtmh*ox

! half acceleration
      vx = ppart(j,3,k)
      vy = ppart(j,4,k)

      acx = vx + ddx
      acy = vy + ddy

! find inverse gamma
! renormalize magnetic field
! calculate cyclotron frequency
      omzt = qtmh1*oz
! calculate rotation matrix
      omt = omzt*omzt
      anorm = 2.0/(1.0 + omt)
      rot1 = 0.5*(1.0 - omt)
      rot2 = omzt
! new momentum
      v1 = (rot1*acx + rot2*acy)*anorm + ddx
      v2 = (rot1*acy - rot2*acx)*anorm + ddy

! deposit momentum flux, acceleration density, and current density
      qm1 = qm*inv_part_7

      amx = qm1*amx
      dxp = qm1*dxp

      ox = 0.5*(v1 + vx)
      oy = 0.5*(v2 + vy)
      oz = 0.5*(1.0+dex*dex*(ox*ox+oy*oy))*inv_part_7 - 0.5*p7

      ppart(j,6,k) = p7 + oz

      oz = oz*idex

      vx = (v1 - vx)*dti
      vy = (v2 - vy)*dti

      dz = qbm*(dz + (dx*ox+dy*oy)*dex*dex*inv_part_7)


      dx = amx*amy
      dy = dxp*amy

      v1 = ox*ox*inv_part_7
      v2 = oy*oy*inv_part_7
      v3 = ox*oy*inv_part_7

      vx = vx+ox*dz*inv_part_7
      vy = vy+oy*dz*inv_part_7

      samu(1,nn) = samu(1,nn) + v1*dx
      samu(2,nn) = samu(2,nn) + v2*dx
      samu(3,nn) = samu(3,nn) + v3*dx

      sdcu(1,nn) = sdcu(1,nn) + vx*dx
      sdcu(2,nn) = sdcu(2,nn) + vy*dx

      scu(1,nn) = scu(1,nn) + ox*dx
      scu(2,nn) = scu(2,nn) + oy*dx
      scu(3,nn) = scu(3,nn) + oz*dx

      dx = amx*dyp
      samu(1,nn+1) = samu(1,nn+1) + v1*dy
      samu(2,nn+1) = samu(2,nn+1) + v2*dy
      samu(3,nn+1) = samu(3,nn+1) + v3*dy

      sdcu(1,nn+1) = sdcu(1,nn+1) + vx*dy
      sdcu(2,nn+1) = sdcu(2,nn+1) + vy*dy

      scu(1,nn+1) = scu(1,nn+1) + ox*dy
      scu(2,nn+1) = scu(2,nn+1) + oy*dy
      scu(3,nn+1) = scu(3,nn+1) + oz*dy

      dy = dxp*dyp
      samu(1,nn+lxv) = samu(1,nn+lxv) + v1*dx
      samu(2,nn+lxv) = samu(2,nn+lxv) + v2*dx
      samu(3,nn+lxv) = samu(3,nn+lxv) + v3*dx

      sdcu(1,nn+lxv) = sdcu(1,nn+lxv) + vx*dx
      sdcu(2,nn+lxv) = sdcu(2,nn+lxv) + vy*dx

      scu(1,nn+lxv) = scu(1,nn+lxv) + ox*dx
      scu(2,nn+lxv) = scu(2,nn+lxv) + oy*dx
      scu(3,nn+lxv) = scu(3,nn+lxv) + oz*dx

      samu(1,nn+1+lxv) = samu(1,nn+1+lxv) + v1*dy
      samu(2,nn+1+lxv) = samu(2,nn+1+lxv) + v2*dy
      samu(3,nn+1+lxv) = samu(3,nn+1+lxv) + v3*dy

      sdcu(1,nn+1+lxv) = sdcu(1,nn+1+lxv) + vx*dy
      sdcu(2,nn+1+lxv) = sdcu(2,nn+1+lxv) + vy*dy

      scu(1,nn+1+lxv) = scu(1,nn+1+lxv) + ox*dy
      scu(2,nn+1+lxv) = scu(2,nn+1+lxv) + oy*dy
      scu(3,nn+1+lxv) = scu(3,nn+1+lxv) + oz*dy
  111 continue
! deposit currents to interior points in global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      do 113 j = 2, mm
      do 112 i = 2, nn
      amu(1,i+noffp+nxv*(j+moffp-1)) = amu(1,i+noffp+nxv*(j+moffp-1)) + 
     1samu(1,i+lxv*(j-1))
      amu(2,i+noffp+nxv*(j+moffp-1)) = amu(2,i+noffp+nxv*(j+moffp-1)) + 
     1samu(2,i+lxv*(j-1))
      amu(3,i+noffp+nxv*(j+moffp-1)) = amu(3,i+noffp+nxv*(j+moffp-1)) + 
     1samu(3,i+lxv*(j-1))

      dcu(1,i+noffp+nxv*(j+moffp-1)) = dcu(1,i+noffp+nxv*(j+moffp-1)) + 
     1sdcu(1,i+lxv*(j-1))
      dcu(2,i+noffp+nxv*(j+moffp-1)) = dcu(2,i+noffp+nxv*(j+moffp-1)) + 
     1sdcu(2,i+lxv*(j-1))

      cu(1,i+noffp+nxv*(j+moffp-1)) = cu(1,i+noffp+nxv*(j+moffp-1)) + 
     1scu(1,i+lxv*(j-1))
      cu(2,i+noffp+nxv*(j+moffp-1)) = cu(2,i+noffp+nxv*(j+moffp-1)) + 
     1scu(2,i+lxv*(j-1))
      cu(3,i+noffp+nxv*(j+moffp-1)) = cu(3,i+noffp+nxv*(j+moffp-1)) + 
     1scu(3,i+lxv*(j-1))
  112 continue
  113 continue
! deposit currents to edge points in global array
      mm = min(my+1,nypmx-moffp)
      do 114 i = 2, nn
!$OMP ATOMIC
      amu(1,i+noffp+nxv*moffp) = amu(1,i+noffp+nxv*moffp) + samu(1,i)
!$OMP ATOMIC
      amu(2,i+noffp+nxv*moffp) = amu(2,i+noffp+nxv*moffp) + samu(2,i)
!$OMP ATOMIC
      amu(3,i+noffp+nxv*moffp) = amu(3,i+noffp+nxv*moffp) + samu(3,i)
!$OMP ATOMIC
      dcu(1,i+noffp+nxv*moffp) = dcu(1,i+noffp+nxv*moffp) + sdcu(1,i)
!$OMP ATOMIC
      dcu(2,i+noffp+nxv*moffp) = dcu(2,i+noffp+nxv*moffp) + sdcu(2,i)
!$OMP ATOMIC
      cu(1,i+noffp+nxv*moffp) = cu(1,i+noffp+nxv*moffp) + scu(1,i)
!$OMP ATOMIC
      cu(2,i+noffp+nxv*moffp) = cu(2,i+noffp+nxv*moffp) + scu(2,i)
!$OMP ATOMIC
      cu(3,i+noffp+nxv*moffp) = cu(3,i+noffp+nxv*moffp) + scu(3,i)
      if (mm > my) then
!$OMP ATOMIC
      amu(1,i+noffp+nxv*(mm+moffp-1)) = amu(1,i+noffp+nxv*(mm+moffp-1))
     1 + samu(1,i+lxv*(mm-1))
!$OMP ATOMIC
      amu(2,i+noffp+nxv*(mm+moffp-1)) = amu(2,i+noffp+nxv*(mm+moffp-1))
     1 + samu(2,i+lxv*(mm-1))
!$OMP ATOMIC
      amu(3,i+noffp+nxv*(mm+moffp-1)) = amu(3,i+noffp+nxv*(mm+moffp-1))
     1 + samu(3,i+lxv*(mm-1))
!$OMP ATOMIC
      dcu(1,i+noffp+nxv*(mm+moffp-1)) = dcu(1,i+noffp+nxv*(mm+moffp-1))
     1 + sdcu(1,i+lxv*(mm-1))
!$OMP ATOMIC
      dcu(2,i+noffp+nxv*(mm+moffp-1)) = dcu(2,i+noffp+nxv*(mm+moffp-1))
     1 + sdcu(2,i+lxv*(mm-1))
!$OMP ATOMIC
      cu(1,i+noffp+nxv*(mm+moffp-1)) = cu(1,i+noffp+nxv*(mm+moffp-1))
     1 + scu(1,i+lxv*(mm-1))
!$OMP ATOMIC
      cu(2,i+noffp+nxv*(mm+moffp-1)) = cu(2,i+noffp+nxv*(mm+moffp-1))
     1 + scu(2,i+lxv*(mm-1))
!$OMP ATOMIC
      cu(3,i+noffp+nxv*(mm+moffp-1)) = cu(3,i+noffp+nxv*(mm+moffp-1))
     1 + scu(3,i+lxv*(mm-1))
      endif
  114 continue
      nn = min(mx+1,nxv-noffp)
      do 115 j = 1, mm
!$OMP ATOMIC
      amu(1,1+noffp+nxv*(j+moffp-1)) = amu(1,1+noffp+nxv*(j+moffp-1)) + 
     1samu(1,1+lxv*(j-1))
!$OMP ATOMIC
      amu(2,1+noffp+nxv*(j+moffp-1)) = amu(2,1+noffp+nxv*(j+moffp-1)) + 
     1samu(2,1+lxv*(j-1))
!$OMP ATOMIC
      amu(3,1+noffp+nxv*(j+moffp-1)) = amu(3,1+noffp+nxv*(j+moffp-1)) + 
     1samu(3,1+lxv*(j-1))
!$OMP ATOMIC
      dcu(1,1+noffp+nxv*(j+moffp-1)) = dcu(1,1+noffp+nxv*(j+moffp-1)) + 
     1sdcu(1,1+lxv*(j-1))
!$OMP ATOMIC
      dcu(2,1+noffp+nxv*(j+moffp-1)) = dcu(2,1+noffp+nxv*(j+moffp-1)) + 
     1sdcu(2,1+lxv*(j-1))
!$OMP ATOMIC
      cu(1,1+noffp+nxv*(j+moffp-1)) = cu(1,1+noffp+nxv*(j+moffp-1)) + 
     1scu(1,1+lxv*(j-1))
!$OMP ATOMIC
      cu(2,1+noffp+nxv*(j+moffp-1)) = cu(2,1+noffp+nxv*(j+moffp-1)) + 
     1scu(2,1+lxv*(j-1))
!$OMP ATOMIC
      cu(3,1+noffp+nxv*(j+moffp-1)) = cu(3,1+noffp+nxv*(j+moffp-1)) + 
     1scu(3,1+lxv*(j-1))
      if (nn > mx) then
!$OMP ATOMIC
      amu(1,nn+noffp+nxv*(j+moffp-1)) = amu(1,nn+noffp+nxv*(j+moffp-1))
     1 + samu(1,nn+lxv*(j-1))
!$OMP ATOMIC
      amu(2,nn+noffp+nxv*(j+moffp-1)) = amu(2,nn+noffp+nxv*(j+moffp-1))
     1 + samu(2,nn+lxv*(j-1))
!$OMP ATOMIC
      amu(3,nn+noffp+nxv*(j+moffp-1)) = amu(3,nn+noffp+nxv*(j+moffp-1))
     1+ samu(3,nn+lxv*(j-1))
!$OMP ATOMIC
      dcu(1,nn+noffp+nxv*(j+moffp-1)) = dcu(1,nn+noffp+nxv*(j+moffp-1))
     1 + sdcu(1,nn+lxv*(j-1))
!$OMP ATOMIC
      dcu(2,nn+noffp+nxv*(j+moffp-1)) = dcu(2,nn+noffp+nxv*(j+moffp-1))
     1 + sdcu(2,nn+lxv*(j-1))
!$OMP ATOMIC
      cu(1,nn+noffp+nxv*(j+moffp-1)) = cu(1,nn+noffp+nxv*(j+moffp-1)) +
     1scu(1,nn+lxv*(j-1))
!$OMP ATOMIC
      cu(2,nn+noffp+nxv*(j+moffp-1)) = cu(2,nn+noffp+nxv*(j+moffp-1)) + 
     1scu(2,nn+lxv*(j-1))
!$OMP ATOMIC
      cu(3,nn+noffp+nxv*(j+moffp-1)) = cu(3,nn+noffp+nxv*(j+moffp-1)) + 
     1scu(3,nn+lxv*(j-1))
      endif
  115 continue
  120 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPGRBPPUSHF23L_QP(ppart,fxy,bxy,psit,kpic,ncl,ihole,nof&
     &f,nyp,qbm,dt,dtc,ci,ek,idimp,nppmx,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1&
     &,ntmax,irc,dex)
! vectorization 
! transposed
! for 2-1/2d code, this subroutine updates particle co-ordinates and
! velocities using leap-frog scheme in time and first-order linear
! interpolation in space, for relativistic particles with magnetic field
! Using the Boris Mover.
! with periodic boundary conditions.
! also determines list of particles which are leaving this tile
! OpenMP version using guard cells, for distributed data
! data read in tiles
! particles stored segmented array
! 131 flops/particle, 4 divides, 2 sqrts, 25 loads, 5 stores
! input: all except ncl, ihole, irc, output: ppart, ncl, ihole, ek, irc
! momentum equations used are:
! px(t+dt/2) = rot(1)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
!    rot(2)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
!    rot(3)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
!    .5*(q/m)*fx(x(t),y(t))*dt)
! py(t+dt/2) = rot(4)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
!    rot(5)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
!    rot(6)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
!    .5*(q/m)*fy(x(t),y(t))*dt)
! pz(t+dt/2) = rot(7)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
!    rot(8)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
!    rot(9)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
!    .5*(q/m)*fz(x(t),y(t))*dt)
! where q/m is charge/mass, and the rotation matrix is given by:
!    rot(1) = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
!    rot(2) = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
!    rot(3) = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
!    rot(4) = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
!    rot(5) = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
!    rot(6) = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
!    rot(7) = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
!    rot(8) = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
!    rot(9) = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
! and om**2 = omx**2 + omy**2 + omz**2
! the rotation matrix is determined by:
! omx = (q/m)*bx(x(t),y(t))*gami, omy = (q/m)*by(x(t),y(t))*gami, and
! omz = (q/m)*bz(x(t),y(t))*gami,
! where gami = 1./sqrt(1.+(px(t)*px(t)+py(t)*py(t)+pz(t)*pz(t))*ci*ci)
! position equations used are:
! x(t+dt) = x(t) + px(t+dt/2)*dtg
! y(t+dt) = y(t) + py(t+dt/2)*dtg
! where dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2)+
! pz(t+dt/2)*pz(t+dt/2))*ci*ci)
! fx(x(t),y(t)), fy(x(t),y(t)), and fz(x(t),y(t))
! bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
! are approximated by interpolation from the nearest grid points:
! fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
!    + dx*fx(n+1,m+1))
! where n,m = leftmost grid points and dx = x-n, dy = y-m
! similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
! ppart(1,n,m) = position x of particle n in partition in tile m
! ppart(2,n,m) = position y of particle n in partition in tile m
! ppart(3,n,m) = momentum vx of particle n in partition in tile m
! ppart(4,n,m) = momentum vy of particle n in partition in tile m
! ppart(5,n,m) = momentum vz of particle n in partition in tile m
! fxy(1,j,k) = x component of force/charge at grid (j,kk)
! fxy(2,j,k) = y component of force/charge at grid (j,kk)
! fxy(3,j,k) = z component of force/charge at grid (j,kk)
! that is, convolution of electric field over particle shape,
! where kk = k + noff - 1
! bxy(1,j,k) = x component of magnetic field at grid (j,kk)
! bxy(2,j,k) = y component of magnetic field at grid (j,kk)
! bxy(3,j,k) = z component of magnetic field at grid (j,kk)
! that is, the convolution of magnetic field over particle shape,
! where kk = k + noff - 1
! kpic(k) = number of particles in tile k
! ncl(i,k) = number of particles going to destination i, tile k
! ihole(1,:,k) = location of hole in array left by departing particle
! ihole(2,:,k) = destination of particle leaving hole
! ihole(1,1,k) = ih, number of holes left (error, if negative)
! noff = lowermost global gridpoint in particle partition.
! nyp = number of primary (complete) gridpoints in particle partition
! qbm = particle charge/mass ratio
! dt = time interval between successive calculations
! dtc = time interval between successive co-ordinate calculations
! ci = reciprical of velocity of light
! kinetic energy/mass at time t is also calculated, using
! ek = gami*sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
!      (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 +
!      (pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)/(1. + gami)
! idimp = size of phase space = 5
! nppmx = maximum number of particles in tile
! nx/ny = system length in x/y direction
! mx/my = number of grids in sorting cell in x/y
! nxv = first dimension of field arrays, must be >= nx+1
! nypmx = maximum size of particle partition, including guard cells.
! mx1 = (system length in x direction - 1)/mx + 1
! mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
! ntmax = size of hole array for particles leaving tiles
! irc = maximum overflow, returned only if error occurs, when irc > 0
! optimized version
      implicit none
      integer noff, nyp, idimp, nppmx, nx, ny, mx, my, nxv, nypmx
      integer mx1, mxyp1, ntmax, irc
      real qbm, dt, dtc, ci, ek, dex
      real ppart, fxy, bxy, psit
      integer kpic, ncl, ihole
      dimension ppart(nppmx,idimp,mxyp1)
!      dimension fxy(2,nxv,nypmx), bxy(3,nxv,nypmx)
      dimension fxy(2,nxv*nypmx), bxy(3,nxv*nypmx)
!      dimension psit(nxv,nypmx)
      dimension psit(nxv*nypmx)
      dimension kpic(mxyp1), ncl(8,mxyp1)
      dimension ihole(2,ntmax+1,mxyp1)
! local data
!      integer MXV, MYV
!      parameter(MXV=33,MYV=33)
      integer npblk, lvect
      parameter(npblk=32,lvect=4)
      integer ipp, joff, nps, m, lxv
      integer noffp, moffp, nppp
      integer mnoff, i, j, k, ih, nh, nn, mm
      real qtmh, ci2, dxp, dyp, amx, amy
      real dx, dy, dz, ox, oy, oz, acx, acy, acz, p2, gami, qtmg, dtg
      real omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real anx, any, edgelx, edgely, edgerx, edgery
      real x, y, vx, vy, vz
!      real sfxy, sbxy, spsit
      real sfxy, sbxy
      dimension sfxy(3,(mx+1)*(my+1)), sbxy(3,(mx+1)*(my+1))
!      dimension sfxy(2,mx+1,my+1), sbxy(3,mx+1,my+1), spsit(mx+1,my+1)
!     dimension sfxy(3,mx+1,my+1), sbxy(3,mx+1,my+1)
! scratch arrays
      integer n
      real s1, s2, t
      dimension n(npblk), s1(npblk,lvect), s2(npblk,lvect), t(npblk,4)
      double precision sum1, sum2
      real idex,dtc1,qtmh1,qtmh2,p6,p7
      
      lxv = mx + 1
      idex = 1.0/dex
      qtmh = 0.5*qbm*dt
      ci2 = ci*ci
      anx = real(nx)
      any = real(ny)
      sum2 = 0.0d0
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV)) return
! loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,noffp,moffp,nppp,nn,mm,ih,nh,mnoff,x,y,dxp,dyp,amx,
!$OMP& amy,dx,dy,dz,ox,oy,oz,acx,acy,acz,omxt,omyt,omzt,omt,anorm,rot1,
!$OMP& rot2,rot3,rot4,rot5,rot6,rot7,rot8,rot9,edgelx,edgely,edgerx,
!$OMP& edgery,sum1,sfxy,sbxy,vx,vy,vz,
!$OMP& ipp,joff,nps,n,s1,s2,t,m,
!$OMP& qtmh1,qtmh2,p6,p7,dtc1)
!$OMP& REDUCTION(+:sum2)
      do 170 k = 1, mxyp1
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(k)
      nn = min(mx,nx-noffp)
      mm = min(my,nyp-moffp)
      edgelx = noffp
      edgerx = noffp + nn
      edgely = noff + moffp
      edgery = noff + moffp + mm
      ih = 0
      nh = 0
      mnoff = moffp + noff
! load local fields from global arrays
      do 20 j = 1, mm+1
      do 10 i = 1, nn+1
      sfxy(1,i+lxv*(j-1)) = fxy(1,i+noffp+nxv*(j+moffp-1))
      sfxy(2,i+lxv*(j-1)) = fxy(2,i+noffp+nxv*(j+moffp-1))
      sfxy(3,i+lxv*(j-1)) = psit(i+noffp+nxv*(j+moffp-1))
   10 continue
   20 continue
      do 40 j = 1, mm+1
      do 30 i = 1, nn+1
      sbxy(1,i+lxv*(j-1)) = bxy(1,i+noffp+nxv*(j+moffp-1))
      sbxy(2,i+lxv*(j-1)) = bxy(2,i+noffp+nxv*(j+moffp-1))
      sbxy(3,i+lxv*(j-1)) = bxy(3,i+noffp+nxv*(j+moffp-1))
   30 continue
   40 continue
! clear counters
      do 50 j = 1, 8
      ncl(j,k) = 0
   50 continue
      sum1 = 0.0d0
      ipp = nppp/npblk
c outer loop over number of full blocks
      do 110 m = 1, ipp
      joff = npblk*(m - 1)
c inner loop over particles in block
      do 60 j = 1, npblk
! find interpolation weights
      x = ppart(j+joff,1,k)
      y = ppart(j+joff,2,k)
      nn = x
      mm = y

      p6 = ppart(j+joff,6,k)
      p7 = ppart(j+joff,7,k)

      qtmh1 = qtmh / p7 
      qtmh2 = qtmh1 * p6
      t(j,3) = qtmh1
      t(j,4) = qtmh2

      dxp = x - real(nn)
      dyp = y - real(mm)
      n(j) = nn - noffp + lxv*(mm - mnoff)
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      s1(j,1) = amx*amy
      s1(j,2) = dxp*amy
      s1(j,3) = amx*dyp
      s1(j,4) = dxp*dyp
      t(j,1) = x
      t(j,2) = y
   60 continue
c find acceleration
      do 80 j = 1, npblk
      nn = n(j)
      mm = nn + lxv - 2
      dx = 0.0
      dy = 0.0
      dz = 0.0
      ox = 0.0
      oy = 0.0
      oz = 0.0
      do 70 i = 1, lvect
      if (i.gt.2) nn = mm
      dx = dx + sfxy(1,i+nn)*s1(j,i)
      dy = dy + sfxy(2,i+nn)*s1(j,i)
      dz = dz + sfxy(3,i+nn)*s1(j,i)
      ox = ox + sbxy(1,i+nn)*s1(j,i)
      oy = oy + sbxy(2,i+nn)*s1(j,i)
      oz = oz + sbxy(3,i+nn)*s1(j,i)
   70 continue
      s1(j,1) = dx
      s1(j,2) = dy
      s1(j,3) = dz*idex
      s2(j,1) = ox*dex
      s2(j,2) = oy*dex
      s2(j,3) = oz
   80 continue
c new velocity
      do 90 j = 1, npblk
      x = t(j,1)
      y = t(j,2)
      qtmh1 = t(j,3)
      qtmh2 = t(j,4)
      
! calculate half impulse
      dx = qtmh2*s1(j,1)
      dy = qtmh2*s1(j,2)
      dz = qtmh2*s1(j,3)
! half acceleration
      acx = ppart(j+joff,3,k) + dx
      acy = ppart(j+joff,4,k) + dy
      acz = ppart(j+joff,5,k) + dz
c calculate cyclotron frequency
      omxt = qtmh1*s2(j,1)
      omyt = qtmh1*s2(j,2)
      omzt = qtmh1*s2(j,3)
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2.0/(1.0 + omt)
      omt = 0.5*(1.0 - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c new velocity
      vx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      vy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      vz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
! update inverse gamma
      dtc1 = dt/(sqrt(1+(vx*vx+vy*vy+vz*vz)*dex*dex)-vz*dex)
c new position
      s1(j,1) = x + vx*dtc1
      s1(j,2) = y + vy*dtc1
      s2(j,1) = vx
      s2(j,2) = vy
      s2(j,3) = vz
   90 continue
! check boundary conditions
!dir$ novector
      do 100 j = 1, npblk
      dx = s1(j,1)
      dy = s1(j,2)
! find particles going out of bounds
      mm = 0
! count how many particles are going in each direction in ncl
! save their address and destination in ihole
! use periodic boundary conditions and check for roundoff error
! mm = direction particle is going
      if (dx.ge.edgerx) then
         if (dx.ge.anx) dx = dx - anx
         mm = 2
      else if (dx.lt.edgelx) then
         if (dx.lt.0.0) then
            dx = dx + anx
            if (dx.lt.anx) then
               mm = 1
            else
               dx = 0.0
            endif
         else
            mm = 1
         endif
      endif
      if (dy.ge.edgery) then
         if (dy.ge.any) dy = dy - any
         mm = mm + 6
      else if (dy.lt.edgely) then
         if (dy.lt.0.0) then
            dy = dy + any
            if (dy.lt.any) then
               mm = mm + 3
            else
               dy = 0.0
            endif
         else
            mm = mm + 3
         endif
      endif
! set new position
      ppart(j+joff,1,k) = dx
      ppart(j+joff,2,k) = dy
c set new velocity
      ppart(j+joff,3,k) = s2(j,1)
      ppart(j+joff,4,k) = s2(j,2)
      ppart(j+joff,5,k) = s2(j,3)
! increment counters
      if (mm.gt.0) then
         ncl(mm,k) = ncl(mm,k) + 1
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(1,ih+1,k) = j + joff
            ihole(2,ih+1,k) = mm
         else
            nh = 1
         endif
      endif
  100 continue
  110 continue
      nps = npblk*ipp + 1
! loop over remaining particles
      do 120 j = nps, nppp
! find interpolation weights
      x = ppart(j,1,k)
      y = ppart(j,2,k)
      nn = x
      mm = y

      p6 = ppart(j,6,k)
      p7 = ppart(j,7,k)

      qtmh1 = qtmh / p7 
      qtmh2 = qtmh1 * p6

      dxp = x - real(nn)
      dyp = y - real(mm)
      nn = nn - noffp + 1 + lxv*(mm - mnoff)
      amx = 1.0 - dxp
      amy = 1.0 - dyp
c find electric field
      dx = amx*sfxy(1,nn)
      dy = amx*sfxy(2,nn)
      dz = amx*sfxy(3,nn)
      dx = amy*(dxp*sfxy(1,nn+1) + dx)
      dy = amy*(dxp*sfxy(2,nn+1) + dy)
      dz = amy*(dxp*sfxy(3,nn+1) + dz)
      acx = amx*sfxy(1,nn+lxv)
      acy = amx*sfxy(2,nn+lxv)
      acz = amx*sfxy(3,nn+lxv)
      dx = dx + dyp*(dxp*sfxy(1,nn+1+lxv) + acx) 
      dy = dy + dyp*(dxp*sfxy(2,nn+1+lxv) + acy)
      dz = dz + dyp*(dxp*sfxy(3,nn+1+lxv) + acz)
      dz = dz * idex
c find magnetic field
      ox = amx*sbxy(1,nn)
      oy = amx*sbxy(2,nn)
      oz = amx*sbxy(3,nn)
      ox = amy*(dxp*sbxy(1,nn+1) + ox)
      oy = amy*(dxp*sbxy(2,nn+1) + oy)
      oz = amy*(dxp*sbxy(3,nn+1) + oz)
      acx = amx*sbxy(1,nn+lxv)
      acy = amx*sbxy(2,nn+lxv)
      acz = amx*sbxy(3,nn+lxv)
      ox = ox + dyp*(dxp*sbxy(1,nn+1+lxv) + acx) 
      oy = oy + dyp*(dxp*sbxy(2,nn+1+lxv) + acy)
      oz = oz + dyp*(dxp*sbxy(3,nn+1+lxv) + acz)
      ox = ox * dex
      oy = oy * dex
! calculate half impulse
      dx = qtmh2*dx
      dy = qtmh2*dy
      dz = qtmh2*dz
! half acceleration
      acx = ppart(j,3,k) + dx
      acy = ppart(j,4,k) + dy
      acz = ppart(j,5,k) + dz
! find inverse gamma
!      p2 = acx*acx + acy*acy + acz*acz
!      gami = 1.0/sqrt(1.0 + p2*ci2)
! renormalize magnetic field
!      qtmg = qtmh*gami
! time-centered kinetic energy
!      sum1 = sum1 + gami*p2/(1.0 + gami)
! calculate cyclotron frequency
      omxt = qtmh1*ox
      omyt = qtmh1*oy
      omzt = qtmh1*oz
! calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2.0/(1.0 + omt)
      omt = 0.5*(1.0 - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
! new momentum
      vx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      vy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      vz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
! update inverse gamma
      dtc1 = dt/(sqrt(1+(vx*vx+vy*vy+vz*vz)*dex*dex)-vz*dex)
c new position
      dx = x + vx*dtc1
      dy = y + vy*dtc1

!      p2 = dx*dx + dy*dy + dz*dz
!      dtg = dtc/sqrt(1.0 + p2*ci2)
! find particles going out of bounds
      mm = 0
! count how many particles are going in each direction in ncl
! save their address and destination in ihole
! use periodic boundary conditions and check for roundoff error
! mm = direction particle is going
      if (dx.ge.edgerx) then
         if (dx.ge.anx) dx = dx - anx
         mm = 2
      else if (dx.lt.edgelx) then
         if (dx.lt.0.0) then
            dx = dx + anx
            if (dx.lt.anx) then
               mm = 1
            else
               dx = 0.0
            endif
         else
            mm = 1
         endif
      endif
      if (dy.ge.edgery) then
         if (dy.ge.any) dy = dy - any
         mm = mm + 6
      else if (dy.lt.edgely) then
         if (dy.lt.0.0) then
            dy = dy + any
            if (dy.lt.any) then
               mm = mm + 3
            else
               dy = 0.0
            endif
         else
            mm = mm + 3
         endif
      endif
! set new position
      ppart(j,1,k) = dx
      ppart(j,2,k) = dy
c set new velocity
      ppart(j,3,k) = vx
      ppart(j,4,k) = vy
      ppart(j,5,k) = vz      
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
  120 continue
      sum2 = sum2 + sum1
! set error and end of file flag
! ihole overflow
      if (nh.gt.0) then
         irc = ih
         ih = -ih
      endif
      ihole(1,1,k) = ih
  170 continue
!$OMP END PARALLEL DO
! normalize kinetic energy
      ek = ek + sum2
      return
      end
!-----------------------------------------------------------------------
      subroutine PPPORDER2LA(ppart,ppbuff,sbufl,sbufr,kpic,ncl,ihole,   &
     &ncll,nclr,noff,nyp,idimp,nppmx,nx,ny,mx,my,mx1,myp1,npbmx,ntmax,  &
     &nbmax,irc)
! this subroutine performs first part of a particle sort by x,y grid
! in tiles of mx, my
! linear interpolation, with periodic boundary conditions
! for distributed data, with 1d domain decomposition in y.
! tiles are assumed to be arranged in 2D linear memory
! this part of the algorithm has 3 steps.  first, one finds particles
! leaving tile and stores their number in each directon, location, and
! destination in ncl and ihole.  then, a prefix scan of ncl is performed
! and departing particles are buffered in ppbuff in direction order.
! finally, we buffer particles leaving the processor in sbufl and sbufr,
! and store particle number offsets in ncll and nclr.
! input: all except ppbuff, sbufl, sbufr, ncl, ihole, ncll, nclr, irc
! output: ppart, ppbuff, sbufl, sbufr, ncl, ihole, ncll, nclr, irc
! ppart(1,n,k) = position x of particle n in tile k
! ppart(2,n,k) = position y of particle n in tile k
! ppbuff(i,n,k) = i co-ordinate of particle n in tile k
! sbufl = buffer for particles being sent to lower processor
! sbufr = buffer for particles being sent to upper processor
! kpic(k) = number of particles in tile k
! ncl(i,k) = number of particles going to destination i, tile k
! ihole(1,:,k) = location of hole in array left by departing particle
! ihole(2,:,k) = direction destination of particle leaving hole
! all for tile k
! ihole(1,1,k) = ih, number of holes left (error, if negative)
! ncll = number offset being sent to lower processor
! nclr = number offset being sent to upper processor
! noff = lowermost global gridpoint in particle partition.
! nyp = number of primary (complete) gridpoints in particle partition
! idimp = size of phase space = 4
! nppmx = maximum number of particles in tile
! nx/ny = system length in x/y direction
! mx/my = number of grids in sorting cell in x/y
! mx1 = (system length in x direction - 1)/mx + 1
! myp1 = (partition length in y direction - 1)/my + 1
! npbmx = size of buffer array ppbuff
! ntmax = size of hole array for particles leaving tiles
! nbmax =  size of buffers for passing particles between processors
! irc = maximum overflow, returned only if error occurs, when irc > 0
      implicit none
      integer noff, nyp, idimp, nppmx, nx, ny, mx, my, mx1, myp1, npbmx
      integer ntmax, nbmax, irc
      real ppart, ppbuff, sbufl, sbufr
      integer kpic, ncl, ihole, ncll, nclr
      dimension ppart(nppmx,idimp,mx1*myp1)
      dimension ppbuff(idimp,npbmx,mx1*myp1)
      dimension sbufl(idimp,nbmax), sbufr(idimp,nbmax)
      dimension kpic(mx1*myp1), ncl(8,mx1*myp1)
      dimension ihole(2,ntmax+1,mx1*myp1)
      dimension ncll(3,mx1), nclr(3,mx1)
! local data
      integer mxyp1, noffp, moffp, nppp
      integer i, j, k, ii, ih, nh, ist, nn, mm, isum, ip, j1, kk
      real anx, any, edgelx, edgely, edgerx, edgery, dx, dy
      mxyp1 = mx1*myp1
      anx = real(nx)
      any = real(ny)
! find and count particles leaving tiles and determine destination
! update ppart, ihole, ncl
! loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(j,k,noffp,moffp,nppp,nn,mm,ih,nh,ist,dx,dy,edgelx,edgely,
!$OMP& edgerx,edgery)
      do 30 k = 1, mxyp1
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(k)
      nn = min(mx,nx-noffp)
      mm = min(my,nyp-moffp)
      ih = 0
      nh = 0
      edgelx = noffp
      edgerx = noffp + nn
      edgely = noff + moffp
      edgery = noff + moffp + mm
! clear counters
      do 10 j = 1, 8
      ncl(j,k) = 0
   10 continue
! loop over particles in tile
      do 20 j = 1, nppp
      dx = ppart(j,1,k)
      dy = ppart(j,2,k)
! find particles going out of bounds
      ist = 0
! count how many particles are going in each direction in ncl
! save their address and destination in ihole
! use periodic boundary conditions and check for roundoff error
! ist = direction particle is going
      if (dx.ge.edgerx) then
         if (dx.ge.anx) ppart(j,1,k) = dx - anx
         ist = 2
      else if (dx.lt.edgelx) then
         if (dx.lt.0.0) then
            dx = dx + anx
            if (dx.lt.anx) then
               ist = 1
            else
               dx = 0.0
            endif
            ppart(j,1,k) = dx
         else
            ist = 1
         endif
      endif
      if (dy.ge.edgery) then
         if (dy.ge.any) ppart(j,2,k) = dy - any
         ist = ist + 6
      else if (dy.lt.edgely) then
         if (dy.lt.0.0) then
            dy = dy + any
            if (dy.lt.any) then
               ist = ist + 3
            else
               dy = 0.0
            endif
            ppart(j,2,k) = dy
         else
            ist = ist + 3
         endif
      endif
      if (ist.gt.0) then
         ncl(ist,k) = ncl(ist,k) + 1
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(1,ih+1,k) = j
            ihole(2,ih+1,k) = ist
         else
            nh = 1
         endif
      endif
   20 continue
! set error and end of file flag
      if (nh.gt.0) then
         irc = ih
         ih = -ih
      endif
      ihole(1,1,k) = ih
   30 continue
!$OMP END PARALLEL DO
! ihole overflow
      if (irc.gt.0) return
!
! buffer particles that are leaving tile: update ppbuff, ncl
! loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,isum,ist,nh,ip,j1,ii)
      do 70 k = 1, mxyp1
! find address offset for ordered ppbuff array
      isum = 0
      do 40 j = 1, 8
      ist = ncl(j,k)
      ncl(j,k) = isum
      isum = isum + ist
   40 continue
      nh = ihole(1,1,k)
      ip = 0
! loop over particles leaving tile
      do 60 j = 1, nh
! buffer particles that are leaving tile, in direction order
      j1 = ihole(1,j+1,k)
      ist = ihole(2,j+1,k)
      ii = ncl(ist,k) + 1
      if (ii.le.npbmx) then
         do 50 i = 1, idimp
         ppbuff(i,ii,k) = ppart(j1,i,k)
   50    continue
      else
         ip = 1
      endif
      ncl(ist,k) = ii
   60 continue
! set error
      if (ip.gt.0) irc = ncl(8,k)
   70 continue
!$OMP END PARALLEL DO
! ppbuff overflow
      if (irc.gt.0) return
!
! buffer particles and their number leaving the node:
! update sbufl, sbufr, ncll, nclr
      kk = mx1*(myp1 - 1)
!$OMP PARALLEL DO PRIVATE(k)
      do 80 k = 1, mx1
      ncll(1,k) = ncl(5,k) - ncl(2,k)
      nclr(1,k) = ncl(8,k+kk) - ncl(5,k+kk)
   80 continue
!$OMP END PARALLEL DO
! perform prefix scan
      kk = 1
   90 if (kk.ge.mx1) go to 110
!$OMP PARALLEL DO PRIVATE(k,ii,nn,mm)
      do 100 k = 1, mx1
      ii = (k - 1)/kk
      nn = kk*ii
      mm = 2*nn + kk - 1
      nn = nn + k + kk
      if (nn.le.mx1) then
         ncll(1,nn) = ncll(1,nn) + ncll(1,mm+1)
         nclr(1,nn) = nclr(1,nn) + nclr(1,mm+1)
      endif
  100 continue
!$OMP END PARALLEL DO
      kk = kk + kk
      go to 90
  110 kk = mx1*(myp1 - 1)
!$OMP PARALLEL DO PRIVATE(i,j,k,ii,nn,mm)
      do 180 k = 1, mx1
      ii = ncl(5,k) - ncl(2,k)
      nn = ncll(1,k) - ii
      do 130 j = 1, min(ii,nbmax-nn)
      do 120 i = 1, idimp
      sbufl(i,j+nn) = ppbuff(i,j+ncl(2,k),k)
  120 continue
  130 continue
      do 140 i = 1, 3
      ncll(i,k) = ncl(i+2,k) - ncl(2,k) + nn
  140 continue
      ii = ncl(8,k+kk) - ncl(5,k+kk)
      mm = nclr(1,k) - ii
      do 160 j = 1, min(ii,nbmax-mm)
      do 150 i = 1, idimp
      sbufr(i,j+mm) = ppbuff(i,j+ncl(5,k+kk),k+kk)
  150 continue
  160 continue
      do 170 i = 1, 3
      nclr(i,k) = ncl(i+5,k+kk) - ncl(5,k+kk) + mm
  170 continue
  180 continue
!$OMP END PARALLEL DO
! sbufl or sbufr overflow
      ii = max(ncll(3,mx1),nclr(3,mx1))
      if (ii.gt.nbmax) then
         irc = ii
      endif
      return
      end      
!-----------------------------------------------------------------------
      subroutine PPPORDERF2LA(ppart,ppbuff,sbufl,sbufr,ncl,ihole,ncll,  &
     &nclr,idimp,nppmx,mx1,myp1,npbmx,ntmax,nbmax,irc)
! this subroutine performs first part of a particle sort by x,y grid
! in tiles of mx, my
! linear interpolation, with periodic boundary conditions
! for distributed data, with 1d domain decomposition in y.
! tiles are assumed to be arranged in 2D linear memory
! this part of the algorithm has 2 steps.  first, a prefix scan of ncl
! is performed and departing particles are buffered in ppbuff in
! direction order. then, we buffer particles leaving the processor in
! sbufl and sbufr, and store particle number offsets in ncll and nclr.
! it assumes that the number, location, and destination of particles 
! leaving a tile have been previously stored in ncl and ihole by the
! PPGPPUSHF2L subroutine.
! input: all except ppbuff, sbufl, sbufr, ncll, nclr, irc
! output: ppart, ppbuff, sbufl, sbufr, ncl, ncll, nclr, irc
! ppart(1,n,k) = position x of particle n in tile k
! ppart(2,n,k) = position y of particle n in tile k
! ppbuff(i,n,k) = i co-ordinate of particle n in tile k
! sbufl = buffer for particles being sent to lower processor
! sbufr = buffer for particles being sent to upper processor
! ncl(i,k) = number of particles going to destination i, tile k
! ihole(1,:,k) = location of hole in array left by departing particle
! ihole(2,:,k) = direction destination of particle leaving hole
! all for tile k
! ihole(1,1,k) = ih, number of holes left (error, if negative)
! ncll = number offset being sent to lower processor
! nclr = number offset being sent to upper processor
! idimp = size of phase space = 4
! nppmx = maximum number of particles in tile
! mx1 = (system length in x direction - 1)/mx + 1
! myp1 = (partition length in y direction - 1)/my + 1
! npbmx = size of buffer array ppbuff
! ntmax = size of hole array for particles leaving tiles
! nbmax =  size of buffers for passing particles between processors
! irc = maximum overflow, returned only if error occurs, when irc > 0
      implicit none
      integer idimp, nppmx, mx1, myp1, npbmx, ntmax, nbmax, irc
      real ppart, ppbuff, sbufl, sbufr
      integer ncl, ihole, ncll, nclr
      dimension ppart(nppmx,idimp,mx1*myp1)
      dimension ppbuff(idimp,npbmx,mx1*myp1)
      dimension sbufl(idimp,nbmax), sbufr(idimp,nbmax)
      dimension ncl(8,mx1*myp1)
      dimension ihole(2,ntmax+1,mx1*myp1)
      dimension ncll(3,mx1), nclr(3,mx1)
! local data
      integer mxyp1
      integer i, j, k, ii, nh, ist, nn, mm, isum, ip, j1, kk
      mxyp1 = mx1*myp1
! buffer particles that are leaving tile: update ppbuff, ncl
! loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,isum,ist,nh,ip,j1,ii)
      do 40 k = 1, mxyp1
! find address offset for ordered ppbuff array
      isum = 0
      do 10 j = 1, 8
      ist = ncl(j,k)
      ncl(j,k) = isum
      isum = isum + ist
   10 continue
      nh = ihole(1,1,k)
      ip = 0
! loop over particles leaving tile
      do 30 j = 1, nh
! buffer particles that are leaving tile, in direction order
      j1 = ihole(1,j+1,k)
      ist = ihole(2,j+1,k)
      ii = ncl(ist,k) + 1
      if (ii.le.npbmx) then
         do 20 i = 1, idimp
         ppbuff(i,ii,k) = ppart(j1,i,k)
   20    continue
      else
         ip = 1
      endif
      ncl(ist,k) = ii
   30 continue
! set error
      if (ip.gt.0) irc = ncl(8,k)
   40 continue
!$OMP END PARALLEL DO
! ppbuff overflow
      if (irc.gt.0) return
!
! buffer particles and their number leaving the node:
! update sbufl, sbufr, ncll, nclr
      kk = mx1*(myp1 - 1)
!$OMP PARALLEL DO PRIVATE(k)
      do 50 k = 1, mx1
      ncll(1,k) = ncl(5,k) - ncl(2,k)
      nclr(1,k) = ncl(8,k+kk) - ncl(5,k+kk)
   50 continue
!$OMP END PARALLEL DO
! perform prefix scan
      kk = 1
   60 if (kk.ge.mx1) go to 80
!$OMP PARALLEL DO PRIVATE(k,ii,nn,mm)
      do 70 k = 1, mx1
      ii = (k - 1)/kk
      nn = kk*ii
      mm = 2*nn + kk - 1
      nn = nn + k + kk
      if (nn.le.mx1) then
         ncll(1,nn) = ncll(1,nn) + ncll(1,mm+1)
         nclr(1,nn) = nclr(1,nn) + nclr(1,mm+1)
      endif
   70 continue
!$OMP END PARALLEL DO
      kk = kk + kk
      go to 60
   80 kk = mx1*(myp1 - 1)
!$OMP PARALLEL DO PRIVATE(i,j,k,ii,nn,mm)
      do 150 k = 1, mx1
      ii = ncl(5,k) - ncl(2,k)
      nn = ncll(1,k) - ii
      do 100 j = 1, min(ii,nbmax-nn)
      do 90 i = 1, idimp
      sbufl(i,j+nn) = ppbuff(i,j+ncl(2,k),k)
   90 continue
  100 continue
      do 110 i = 1, 3
      ncll(i,k) = ncl(i+2,k) - ncl(2,k) + nn
  110 continue
      ii = ncl(8,k+kk) - ncl(5,k+kk)
      mm = nclr(1,k) - ii
      do 130 j = 1, min(ii,nbmax-mm)
      do 120 i = 1, idimp
      sbufr(i,j+mm) = ppbuff(i,j+ncl(5,k+kk),k+kk)
  120 continue
  130 continue
      do 140 i = 1, 3
      nclr(i,k) = ncl(i+5,k+kk) - ncl(5,k+kk) + mm
  140 continue
  150 continue
!$OMP END PARALLEL DO
! sbufl or sbufr overflow
      ii = max(ncll(3,mx1),nclr(3,mx1))
      if (ii.gt.nbmax) then
         irc = ii
      endif
      return
      end
!-----------------------------------------------------------------------
      subroutine PPPORDER2LB(ppart,ppbuff,rbufl,rbufr,kpic,ncl,ihole,   &
     &mcll,mclr,idimp,nppmx,mx1,myp1,npbmx,ntmax,nbmax,irc)
! this subroutine performs second part of a particle sort by x,y grid
! in tiles of mx, my
! linear interpolation, with periodic boundary conditions
! for distributed data, with 1d domain decomposition in y.
! tiles are assumed to be arranged in 2D linear memory
! incoming particles from other tiles are copied from ppbuff, rbufl, and
! rbufr into ppart
! input: all except ppart, kpic, irc
! output: ppart, kpic, irc
! ppart(1,n,k) = position x of particle n in tile k
! ppart(2,n,k) = position y of particle n in tile k
! ppbuff(i,n,k) = i co-ordinate of particle n in tile k
! rbufl = buffer for particles being received from lower processor
! rbufr = buffer for particles being received from upper processor
! kpic(k) = number of particles in tile k
! ncl(i,k) = number of particles going to destination i, tile k
! ihole(1,:,k) = location of hole in array left by departing particle
! ihole(2,:,k) = direction destination of particle leaving hole
! all for tile k
! ihole(1,1,k) = ih, number of holes left (error, if negative)
! mcll = number offset being received from lower processor
! mclr = number offset being received from upper processor
! idimp = size of phase space = 4
! nppmx = maximum number of particles in tile
! mx1 = (system length in x direction - 1)/mx + 1
! myp1 = (partition length in y direction - 1)/my + 1
! npbmx = size of buffer array ppbuff
! ntmax = size of hole array for particles leaving tiles
! nbmax =  size of buffers for passing particles between processors
! irc = maximum overflow, returned only if error occurs, when irc > 0
      implicit none
      integer idimp, nppmx, mx1, myp1, npbmx
      integer ntmax, nbmax, irc
      real ppart, ppbuff, rbufl, rbufr
      integer kpic, ncl, ihole, mcll, mclr
      dimension ppart(nppmx,idimp,mx1*myp1)
      dimension ppbuff(idimp,npbmx,mx1*myp1)
      dimension rbufl(idimp,nbmax), rbufr(idimp,nbmax)
      dimension kpic(mx1*myp1), ncl(8,mx1*myp1)
      dimension ihole(2,ntmax+1,mx1*myp1)
      dimension mcll(3,mx1), mclr(3,mx1)
! local data
      integer mxyp1, nppp, ncoff, noff, moff
      integer i, j, k, ii, kx, ky, ih, nh, ist
      integer ip, j1, j2, kxl, kxr, kk, kl, kr
      integer ks
      dimension ks(8)
      mxyp1 = mx1*myp1
! copy incoming particles from buffer into ppart: update ppart, kpic
! loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,ii,kk,nppp,kx,ky,kl,kr,kxl,kxr,ih,nh,ncoff,noff,
!$OMP& moff,ist,j1,j2,ip,ks)
      do 200 k = 1, mxyp1
      nppp = kpic(k)
      ky = (k - 1)/mx1 + 1
! loop over tiles in y
      kk = (ky - 1)*mx1
! find tile above
      kl = (ky - 2)*mx1
! find tile below
      kr = ky*mx1
! loop over tiles in x, assume periodic boundary conditions
      kx = k - (ky - 1)*mx1
      kxl = kx - 1 
      if (kxl.lt.1) kxl = kxl + mx1
      kxr = kx + 1
      if (kxr.gt.mx1) kxr = kxr - mx1
! find tile number for different directions
      ks(1) = kxr + kk
      ks(2) = kxl + kk
      ks(3) = kx + kr
      ks(4) = kxr + kr
      ks(5) = kxl + kr
      ks(6) = kx + kl
      ks(7) = kxr + kl
      ks(8) = kxl + kl
! loop over directions
      nh = ihole(1,1,k)
      noff = 0
      moff = 0
      if (ky.eq.1) then
         if (kx.gt.1) noff = mcll(3,kx-1)
      endif
      if (ky.eq.myp1) then
         if (kx.gt.1) moff = mclr(3,kx-1)
      endif
      ncoff = 0
      ih = 0
      ist = 0
      j1 = 0
      do 170 ii = 1, 8
! ip = number of particles coming from direction ii
      if (ks(ii).le.0) then
         if (ii.gt.6) noff = mcll(ii-6,ks(ii)+mx1)
         ip = mcll(ii-5,ks(ii)+mx1) - noff
      else if (ks(ii).gt.mxyp1) then
         if (ii.gt.3) moff = mclr(ii-3,ks(ii)-mxyp1)
         ip = mclr(ii-2,ks(ii)-mxyp1) - moff
      else
         if (ii.gt.1) ncoff = ncl(ii-1,ks(ii))
         ip = ncl(ii,ks(ii)) - ncoff
      endif
      do 160 j = 1, ip
      ih = ih + 1
! insert incoming particles into holes
      if (ih.le.nh) then
         j1 = ihole(1,ih+1,k)
! place overflow at end of array
      else
         j1 = nppp + 1
         nppp = j1
      endif
      if (j1.le.nppmx) then
         if (ks(ii).le.0) then
            do 130 i = 1, idimp
            ppart(j1,i,k) = rbufl(i,j+noff)
  130       continue
         else if (ks(ii).gt.mxyp1) then
            do 140 i = 1, idimp
            ppart(j1,i,k) = rbufr(i,j+moff)
  140       continue
         else
            do 150 i = 1, idimp
            ppart(j1,i,k) = ppbuff(i,j+ncoff,ks(ii))
  150       continue
         endif
      else
         ist = 1
      endif
  160 continue
  170 continue
! set error
      if (ist.gt.0) irc = j1
! fill up remaining holes in particle array with particles from bottom
      if (ih.lt.nh) then
         ip = nh - ih
         do 190 j = 1, ip
         j1 = nppp - j + 1
         j2 = ihole(1,nh-j+2,k)
         if (j1.gt.j2) then
! move particle only if it is below current hole
            do 180 i = 1, idimp
            ppart(j2,i,k) = ppart(j1,i,k)
  180       continue
         endif
  190    continue
         nppp = nppp - ip
      endif
      kpic(k) = nppp
  200 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
! deposit psi on particles
      subroutine WPGPSIPOST2L_QP(ppart,psi,kpic,qbm,noff,nyp,idimp,nppmx
     1,nx,mx,my,nxv,nypmx,mx1,mxyp1,dex)
! vecterization
! transposed
      implicit none
      integer noff, nyp, idimp, nppmx, nx, mx, my, nxv, nypmx
      integer mx1, mxyp1
      real dex,qbm
      real ppart, psi
      integer kpic
      dimension ppart(nppmx,idimp,mxyp1)
      dimension psi(nxv*nypmx)
      dimension kpic(mxyp1)
! local data
!      integer MXV, MYV
!      parameter(MXV=33,MYV=33)
      integer npblk, lvect
      parameter(npblk=32,lvect=4)
      integer ipp, joff, nps, m, lxv
      integer noffp, moffp, nppp
      integer mnoff, i, j, k, nn, mm
      real dxp, dyp, amx, amy, acx
      real dx, dy, dz
      real x, y, vx, vy, vz
      real spsi
      dimension spsi((mx+1)*(my+1))
      real dx2
c scratch arrays
      integer n
      real s1, s2, t
      dimension n(npblk), s1(npblk,lvect), s2(npblk,lvect), t(npblk,2)

      lxv = mx + 1
      dx2 = dex * dex
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV)) return
! loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,noffp,moffp,nppp,nn,mm,mnoff,x,y,vx,vy,
!$OMP& ipp,joff,nps,n,s1,s2,t,m,
!$OMP& dxp,dyp,amx,amy,dx,acx,spsi)
      do 120 k = 1, mxyp1
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(k)
      mnoff = moffp + noff
! load local fields from global arrays
      nn = min(mx,nx-noffp) + 1
      mm = min(my,nyp-moffp) + 1
      do 20 j = 1, mm
      do 10 i = 1, nn
      spsi(i+lxv*(j-1)) = psi(i+noffp+nxv*(j+moffp-1))
   10 continue
   20 continue
      ipp = nppp/npblk
c outer loop over number of full blocks
      do 110 m = 1, ipp
      joff = npblk*(m - 1)
c inner loop over particles in block
      do 60 j = 1, npblk
! find interpolation weights
      x = ppart(j+joff,1,k)
      y = ppart(j+joff,2,k)
      vx = ppart(j+joff,3,k)
      vy = ppart(j+joff,4,k)
      nn = x
      mm = y
      dxp = x - real(nn)
      dyp = y - real(mm)
      n(j) = nn - noffp + lxv*(mm - mnoff)
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      s1(j,1) = amx*amy
      s1(j,2) = dxp*amy
      s1(j,3) = amx*dyp
      s1(j,4) = dxp*dyp
      t(j,1) = vx
      t(j,2) = vy
   60 continue
c find acceleration
      do 80 j = 1, npblk
      nn = n(j)
      mm = nn + lxv - 2
      dx = 0.0
      do 70 i = 1, lvect
      if (i.gt.2) nn = mm
      dx = dx + spsi(i+nn)*s1(j,i)
   70 continue
      s1(j,1) = -dx*dx2*qbm
   80 continue
      do 90 j = 1, npblk
      vx = t(j,1)
      vy = t(j,2)
      dx = s1(j,1)
      ppart(j+joff,7,k) = 1.0 + dx
      ppart(j+joff,6,k) = dx+(dx2*(vx**2+vy**2)-dx**2-2.*dx)/(2.*(1.+dx)&
     &)+1.0
   90 continue
  110 continue
! loop over particles in tile
      nps = npblk*ipp + 1
      do 115 j = nps, nppp
! find interpolation weights
      x = ppart(j,1,k)
      y = ppart(j,2,k)
      vx = ppart(j,3,k)
      vy = ppart(j,4,k)
      nn = x
      mm = y
      dxp = x - real(nn)
      dyp = y - real(mm)
      nn = nn - noffp + 1 + lxv*(mm - mnoff)
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      
      
! find electric field
      dx = amx*spsi(nn)
      dx = amy*(dxp*spsi(nn+1) + dx)
      acx = amx*spsi(nn+lxv)
      dx = dx + dyp*(dxp*spsi(nn+1+lxv) + acx) 

      dx = - dx*dx2*qbm
      
      ppart(j,7,k) = 1.0 + dx
      ppart(j,6,k) = dx+(dx2*(vx**2+vy**2)-dx**2-2.*dx)/(2.*(1.+dx))+1.0
  115 continue
  120 continue
!$OMP END PARALLEL DO
      return
      end      
!-----------------------------------------------------------------------
      subroutine PPPCOPYOUT2(part,ppart,kpic,npp,npmax,nppmx,idimp,mxyp1&
     &,irc)
! for 2d code, this subroutine copies segmented particle data ppart to
! the array part with original tiled layout
! spatial decomposition in y direction
! input: all except part, npp, irc, output: part, npp, irc
! part(i,j) = i-th coordinate for particle j in partition
! ppart(i,j,k) = i-th coordinate for particle j in partition in tile k
! kpic = number of particles per tile
! npp = number of particles in partition
! npmax = maximum number of particles in each partition
! nppmx = maximum number of particles in tile
! idimp = size of phase space = 5
! mxyp1 = total number of tiles in partition
! irc = maximum overflow, returned only if error occurs, when irc > 0
      implicit none
      integer npp, npmax, nppmx, idimp, mxyp1, irc
      real part, ppart
      integer kpic
      dimension part(idimp,npmax), ppart(nppmx,idimp,mxyp1)
      dimension kpic(mxyp1)
! local data
      integer i, j, k, npoff, nppp, ne, ierr
      npoff = 0
      ierr = 0
! loop over tiles
      do 30 k = 1, mxyp1
      nppp = kpic(k)
      ne = nppp + npoff
      if (ne.gt.npmax) ierr = max(ierr,ne-npmax)
      if (ierr.gt.0) nppp = 0
! loop over particles in tile
      do 20 j = 1, nppp
      do 10 i = 1, idimp
      part(i,j+npoff) = ppart(j,i,k)
   10 continue
   20 continue
      npoff = npoff + nppp
   30 continue
      npp = npoff
      if (ierr.gt.0) irc = ierr
      return
      end
!-----------------------------------------------------------------------      