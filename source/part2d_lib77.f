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
      dimension part(idimp,npmax), ppart(idimp,nppmx,mxyp1)
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
         ppart(i,ip,m) = part(i,j)
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
      dimension ppart(idimp,nppmx,mx1*myp1)
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
      dx = ppart(1,j,k)
      dy = ppart(2,j,k)
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
      dimension ppart(idimp,nppmx,mxyp1), q(nxv,nypmx), kpic(mxyp1)
! local data
!      integer MXV, MYV
!      parameter(MXV=33,MYV=33)
      integer noffp, moffp, nppp
      integer mnoff, i, j, k, nn, mm
      real x, y, dxp, dyp, amx, amy
      real sq
!     dimension sq(MXV,MYV)
      dimension sq(mx+1,my+1)
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV)) return
! loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,noffp,moffp,nppp,mnoff,nn,mm,x,y,dxp,dyp,amx,amy,
!$OMP& sq)
      do 80 k = 1, mxyp1
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(k)
      mnoff = moffp + noff - 1
! zero out local accumulator
      do 20 j = 1, my+1
      do 10 i = 1, mx+1
      sq(i,j) = 0.0
   10 continue
   20 continue
! loop over particles in tile
      do 30 j = 1, nppp
! find interpolation weights
      x = ppart(1,j,k)
      y = ppart(2,j,k)
      nn = x
      mm = y
      dxp = qm*(x - real(nn))
      dyp = y - real(mm)
      nn = nn - noffp + 1
      mm = mm - mnoff
      amx = qm - dxp
      amy = 1.0 - dyp
! deposit charge within tile to local accumulator
      x = sq(nn,mm) + amx*amy
      y = sq(nn+1,mm) + dxp*amy
      sq(nn,mm) = x
      sq(nn+1,mm) = y
      x = sq(nn,mm+1) + amx*dyp
      y = sq(nn+1,mm+1) + dxp*dyp
      sq(nn,mm+1) = x
      sq(nn+1,mm+1) = y
   30 continue
! deposit charge to interior points in global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      do 50 j = 2, mm
      do 40 i = 2, nn
      q(i+noffp,j+moffp) = q(i+noffp,j+moffp) + sq(i,j)
   40 continue
   50 continue
! deposit charge to edge points in global array
      mm = min(my+1,nypmx-moffp)
      do 60 i = 2, nn
!$OMP ATOMIC
      q(i+noffp,1+moffp) = q(i+noffp,1+moffp) + sq(i,1)
      if (mm > my) then
!$OMP ATOMIC
         q(i+noffp,mm+moffp) = q(i+noffp,mm+moffp) + sq(i,mm)
      endif
   60 continue
      nn = min(mx+1,nxv-noffp)
      do 70 j = 1, mm
!$OMP ATOMIC
      q(1+noffp,j+moffp) = q(1+noffp,j+moffp) + sq(1,j)
      if (nn > mx) then
!$OMP ATOMIC
         q(nn+noffp,j+moffp) = q(nn+noffp,j+moffp) + sq(nn,j)
      endif
   70 continue
   80 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPGRDCJPPOST2L_QP(ppart,fxy,bxy,psit,cu,dcu,amu,kpic,no&
     &ff,nyp,qm,qbm,dt,ci,idimp,nppmx,nx,mx,my,nxv,nypmx,mx1,mxyp1,dex)
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
      dimension ppart(idimp,nppmx,mxyp1)
      dimension fxy(2,nxv,nypmx), bxy(3,nxv,nypmx), psit(nxv,nypmx)
      dimension cu(3,nxv,nypmx), dcu(2,nxv,nypmx), amu(3,nxv,nypmx)
      dimension kpic(mxyp1)
! local data
!      integer MXV, MYV
!      parameter(MXV=33,MYV=33)
      integer noffp, moffp, nppp
      integer mnoff, i, j, k, nn, mm
      real qtmh, dti, ci2, gami, qtmg, gh, dxp, dyp, amx, amy
      real dx, dy, dz, ox, oy, oz
      real acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real x, y, vx, vy, vz, p2, v1, v2, v3, v4
      real sfxy, sbxy, scu, sdcu, samu, spsit
      dimension sfxy(2,mx+1,my+1), sbxy(3,mx+1,my+1)
      dimension spsit(mx+1,my+1)
      dimension scu(3,mx+1,my+1), sdcu(2,mx+1,my+1), samu(3,mx+1,my+1)
!      dimension sfxy(2,MXV,MYV), sbxy(3,MXV,MYV)
!      dimension spsit(MXV,MYV)
!      dimension scu(3,MXV,MYV), sdcu(2,MXV,MYV), samu(3,MXV,MYV)
      real qm1, qtmh1,qtmh2,idex,inv_part_7,ddx,ddy,p6,p7
!     dimension sfxy(3,mx+1,my+1), sbxy(3,mx+1,my+1)
!     dimension scu(3,mx+1,my+1,), sdcu(3,mx+1,my+1), samu(4,mx+1,my+1)
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
!$OMP& omt,anorm,rot1,rot2,spsit,  
!$OMP& sfxy,sbxy,scu,sdcu,samu,qm1,qtmh1,qtmh2,inv_part_7,ddx,
!$OMP& ddy,p6,p7)
      do 120 k = 1, mxyp1
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(k)
      mnoff = moffp + noff - 1
! load local fields from global arrays
      nn = min(mx,nx-noffp) + 1
      mm = min(my,nyp-moffp) + 1
      do 20 j = 1, mm
      do 10 i = 1, nn
      sfxy(1,i,j) = fxy(1,i+noffp,j+moffp)
      sfxy(2,i,j) = fxy(2,i+noffp,j+moffp)
   10 continue
   20 continue
      do 40 j = 1, mm
      do 30 i = 1, nn
      sbxy(1,i,j) = bxy(1,i+noffp,j+moffp)
      sbxy(2,i,j) = bxy(2,i+noffp,j+moffp)
      sbxy(3,i,j) = bxy(3,i+noffp,j+moffp)
   30 continue
   40 continue
      do 46 j = 1, mm
      do 42 i = 1, nn
      spsit(i,j) = psit(i+noffp,j+moffp)
   42 continue
   46 continue
! zero out local accumulators
      do 60 j = 1, my+1
      do 50 i = 1, mx+1
      samu(1,i,j) = 0.0
      samu(2,i,j) = 0.0
      samu(3,i,j) = 0.0
      sdcu(1,i,j) = 0.0
      sdcu(2,i,j) = 0.0
      scu(1,i,j) = 0.0
      scu(2,i,j) = 0.0
      scu(3,i,j) = 0.0
   50 continue
   60 continue
! loop over particles in tile
      do 70 j = 1, nppp
! find interpolation weights
      x = ppart(1,j,k)
      y = ppart(2,j,k)
      p6 = ppart(6,j,k)
      p7 = ppart(7,j,k)
      nn = x
      mm = y
      dxp = x - real(nn)
      dyp = y - real(mm)
      nn = nn - noffp + 1
      mm = mm - mnoff
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      
      inv_part_7 = 1.0 / p7
      qtmh1 = qtmh * inv_part_7
      qtmh2 = qtmh1 * p6
      
! find electric field
      dx = amx*sfxy(1,nn,mm)
      dy = amx*sfxy(2,nn,mm)
      dz = amx*spsit(nn,mm)

      dx = amy*(dxp*sfxy(1,nn+1,mm) + dx)
      dy = amy*(dxp*sfxy(2,nn+1,mm) + dy)
      dz = amy*(dxp*spsit(nn+1,mm) + dz)

      acx = amx*sfxy(1,nn,mm+1)
      acy = amx*sfxy(2,nn,mm+1)
      acz = amx*spsit(nn,mm+1)

      dx = dx + dyp*(dxp*sfxy(1,nn+1,mm+1) + acx) 
      dy = dy + dyp*(dxp*sfxy(2,nn+1,mm+1) + acy)
      dz = dz + dyp*(dxp*spsit(nn+1,mm+1) + acz)

! find magnetic field
      ox = amx*sbxy(1,nn,mm)
      oy = amx*sbxy(2,nn,mm)
      oz = amx*sbxy(3,nn,mm)

      ox = amy*(dxp*sbxy(1,nn+1,mm) + ox)
      oy = amy*(dxp*sbxy(2,nn+1,mm) + oy)
      oz = amy*(dxp*sbxy(3,nn+1,mm) + oz)

      acx = amx*sbxy(1,nn,mm+1)
      acy = amx*sbxy(2,nn,mm+1)
      acz = amx*sbxy(3,nn,mm+1)

      ox = ox + dyp*(dxp*sbxy(1,nn+1,mm+1) + acx) 
      oy = oy + dyp*(dxp*sbxy(2,nn+1,mm+1) + acy)
      oz = oz + dyp*(dxp*sbxy(3,nn+1,mm+1) + acz)

! calculate half impulse
      ddx = (-1.0)*qtmh2*dx + qtmh*oy
      ddy = (-1.0)*qtmh2*dy - qtmh*ox

! half acceleration
      vx = ppart(3,j,k)
      vy = ppart(4,j,k)

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

      ppart(6,j,k) = p7 + oz

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

      samu(1,nn,mm) = samu(1,nn,mm) + v1*dx
      samu(2,nn,mm) = samu(2,nn,mm) + v2*dx
      samu(3,nn,mm) = samu(3,nn,mm) + v3*dx

      sdcu(1,nn,mm) = sdcu(1,nn,mm) + vx*dx
      sdcu(2,nn,mm) = sdcu(2,nn,mm) + vy*dx

      scu(1,nn,mm) = scu(1,nn,mm) + ox*dx
      scu(2,nn,mm) = scu(2,nn,mm) + oy*dx
      scu(3,nn,mm) = scu(3,nn,mm) + oz*dx

      dx = amx*dyp
      samu(1,nn+1,mm) = samu(1,nn+1,mm) + v1*dy
      samu(2,nn+1,mm) = samu(2,nn+1,mm) + v2*dy
      samu(3,nn+1,mm) = samu(3,nn+1,mm) + v3*dy

      sdcu(1,nn+1,mm) = sdcu(1,nn+1,mm) + vx*dy
      sdcu(2,nn+1,mm) = sdcu(2,nn+1,mm) + vy*dy

      scu(1,nn+1,mm) = scu(1,nn+1,mm) + ox*dy
      scu(2,nn+1,mm) = scu(2,nn+1,mm) + oy*dy
      scu(3,nn+1,mm) = scu(3,nn+1,mm) + oz*dy

      dy = dxp*dyp
      samu(1,nn,mm+1) = samu(1,nn,mm+1) + v1*dx
      samu(2,nn,mm+1) = samu(2,nn,mm+1) + v2*dx
      samu(3,nn,mm+1) = samu(3,nn,mm+1) + v3*dx

      sdcu(1,nn,mm+1) = sdcu(1,nn,mm+1) + vx*dx
      sdcu(2,nn,mm+1) = sdcu(2,nn,mm+1) + vy*dx

      scu(1,nn,mm+1) = scu(1,nn,mm+1) + ox*dx
      scu(2,nn,mm+1) = scu(2,nn,mm+1) + oy*dx
      scu(3,nn,mm+1) = scu(3,nn,mm+1) + oz*dx

      samu(1,nn+1,mm+1) = samu(1,nn+1,mm+1) + v1*dy
      samu(2,nn+1,mm+1) = samu(2,nn+1,mm+1) + v2*dy
      samu(3,nn+1,mm+1) = samu(3,nn+1,mm+1) + v3*dy

      sdcu(1,nn+1,mm+1) = sdcu(1,nn+1,mm+1) + vx*dy
      sdcu(2,nn+1,mm+1) = sdcu(2,nn+1,mm+1) + vy*dy

      scu(1,nn+1,mm+1) = scu(1,nn+1,mm+1) + ox*dy
      scu(2,nn+1,mm+1) = scu(2,nn+1,mm+1) + oy*dy
      scu(3,nn+1,mm+1) = scu(3,nn+1,mm+1) + oz*dy
   70 continue
! deposit currents to interior points in global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      do 90 j = 2, mm
      do 80 i = 2, nn
      amu(1,i+noffp,j+moffp) = amu(1,i+noffp,j+moffp) + samu(1,i,j)
      amu(2,i+noffp,j+moffp) = amu(2,i+noffp,j+moffp) + samu(2,i,j)
      amu(3,i+noffp,j+moffp) = amu(3,i+noffp,j+moffp) + samu(3,i,j)

      dcu(1,i+noffp,j+moffp) = dcu(1,i+noffp,j+moffp) + sdcu(1,i,j)
      dcu(2,i+noffp,j+moffp) = dcu(2,i+noffp,j+moffp) + sdcu(2,i,j)

      cu(1,i+noffp,j+moffp) = cu(1,i+noffp,j+moffp) + scu(1,i,j)
      cu(2,i+noffp,j+moffp) = cu(2,i+noffp,j+moffp) + scu(2,i,j)
      cu(3,i+noffp,j+moffp) = cu(3,i+noffp,j+moffp) + scu(3,i,j)
   80 continue
   90 continue
! deposit currents to edge points in global array
      mm = min(my+1,nypmx-moffp)
      do 100 i = 2, nn
!$OMP ATOMIC
      amu(1,i+noffp,1+moffp) = amu(1,i+noffp,1+moffp) + samu(1,i,1)
!$OMP ATOMIC
      amu(2,i+noffp,1+moffp) = amu(2,i+noffp,1+moffp) + samu(2,i,1)
!$OMP ATOMIC
      amu(3,i+noffp,1+moffp) = amu(3,i+noffp,1+moffp) + samu(3,i,1)
!$OMP ATOMIC
      dcu(1,i+noffp,1+moffp) = dcu(1,i+noffp,1+moffp) + sdcu(1,i,1)
!$OMP ATOMIC
      dcu(2,i+noffp,1+moffp) = dcu(2,i+noffp,1+moffp) + sdcu(2,i,1)
!$OMP ATOMIC
      cu(1,i+noffp,1+moffp) = cu(1,i+noffp,1+moffp) + scu(1,i,1)
!$OMP ATOMIC
      cu(2,i+noffp,1+moffp) = cu(2,i+noffp,1+moffp) + scu(2,i,1)
!$OMP ATOMIC
      cu(3,i+noffp,1+moffp) = cu(3,i+noffp,1+moffp) + scu(3,i,1)
      if (mm > my) then
!$OMP ATOMIC
         amu(1,i+noffp,mm+moffp) = amu(1,i+noffp,mm+moffp)              &
     & + samu(1,i,mm)
!$OMP ATOMIC
         amu(2,i+noffp,mm+moffp) = amu(2,i+noffp,mm+moffp)              &
     & + samu(2,i,mm)
!$OMP ATOMIC
         amu(3,i+noffp,mm+moffp) = amu(3,i+noffp,mm+moffp)              &
     & + samu(3,i,mm)
!$OMP ATOMIC
         dcu(1,i+noffp,mm+moffp) = dcu(1,i+noffp,mm+moffp)              &
     & + sdcu(1,i,mm)
!$OMP ATOMIC
         dcu(2,i+noffp,mm+moffp) = dcu(2,i+noffp,mm+moffp)              &
     & + sdcu(2,i,mm)
!$OMP ATOMIC
         cu(1,i+noffp,mm+moffp) = cu(1,i+noffp,mm+moffp) + scu(1,i,mm)
!$OMP ATOMIC
         cu(2,i+noffp,mm+moffp) = cu(2,i+noffp,mm+moffp) + scu(2,i,mm)
!$OMP ATOMIC
         cu(3,i+noffp,mm+moffp) = cu(3,i+noffp,mm+moffp) + scu(3,i,mm)
      endif
  100 continue
      nn = min(mx+1,nxv-noffp)
      do 110 j = 1, mm
!$OMP ATOMIC
      amu(1,1+noffp,j+moffp) = amu(1,1+noffp,j+moffp) + samu(1,1,j)
!$OMP ATOMIC
      amu(2,1+noffp,j+moffp) = amu(2,1+noffp,j+moffp) + samu(2,1,j)
!$OMP ATOMIC
      amu(3,1+noffp,j+moffp) = amu(3,1+noffp,j+moffp) + samu(3,1,j)
!$OMP ATOMIC
      dcu(1,1+noffp,j+moffp) = dcu(1,1+noffp,j+moffp) + sdcu(1,1,j)
!$OMP ATOMIC
      dcu(2,1+noffp,j+moffp) = dcu(2,1+noffp,j+moffp) + sdcu(2,1,j)
!$OMP ATOMIC
      cu(1,1+noffp,j+moffp) = cu(1,1+noffp,j+moffp) + scu(1,1,j)
!$OMP ATOMIC
      cu(2,1+noffp,j+moffp) = cu(2,1+noffp,j+moffp) + scu(2,1,j)
!$OMP ATOMIC
      cu(3,1+noffp,j+moffp) = cu(3,1+noffp,j+moffp) + scu(3,1,j)
      if (nn > mx) then
!$OMP ATOMIC
         amu(1,nn+noffp,j+moffp) = amu(1,nn+noffp,j+moffp)              &
     & + samu(1,nn,j)
!$OMP ATOMIC
         amu(2,nn+noffp,j+moffp) = amu(2,nn+noffp,j+moffp)              &
     & + samu(2,nn,j)
!$OMP ATOMIC
         amu(3,nn+noffp,j+moffp) = amu(3,nn+noffp,j+moffp)              &
     & + samu(3,nn,j)
!$OMP ATOMIC
         dcu(1,nn+noffp,j+moffp) = dcu(1,nn+noffp,j+moffp)              &
     & + sdcu(1,nn,j)
!$OMP ATOMIC
         dcu(2,nn+noffp,j+moffp) = dcu(2,nn+noffp,j+moffp)              &
     & + sdcu(2,nn,j)
!$OMP ATOMIC
         cu(1,nn+noffp,j+moffp) = cu(1,nn+noffp,j+moffp) + scu(1,nn,j)
!$OMP ATOMIC
         cu(2,nn+noffp,j+moffp) = cu(2,nn+noffp,j+moffp) + scu(2,nn,j)
!$OMP ATOMIC
         cu(3,nn+noffp,j+moffp) = cu(3,nn+noffp,j+moffp) + scu(3,nn,j)
      endif
  110 continue
  120 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPGRBPPUSHF23L_QP(ppart,fxy,bxy,psit,kpic,ncl,ihole,nof&
     &f,nyp,qbm,dt,dtc,ci,ek,idimp,nppmx,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1&
     &,ntmax,irc,dex)
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
      dimension ppart(idimp,nppmx,mxyp1)
      dimension fxy(2,nxv,nypmx), bxy(3,nxv,nypmx)
      dimension psit(nxv,nypmx)
      dimension kpic(mxyp1), ncl(8,mxyp1)
      dimension ihole(2,ntmax+1,mxyp1)
! local data
!      integer MXV, MYV
!      parameter(MXV=33,MYV=33)
      integer noffp, moffp, nppp
      integer mnoff, i, j, k, ih, nh, nn, mm
      real qtmh, ci2, dxp, dyp, amx, amy
      real dx, dy, dz, ox, oy, oz, acx, acy, acz, p2, gami, qtmg, dtg
      real omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real anx, any, edgelx, edgely, edgerx, edgery
      real x, y
      real sfxy, sbxy, spsit
      dimension sfxy(2,mx+1,my+1), sbxy(3,mx+1,my+1), spsit(mx+1,my+1)
!     dimension sfxy(3,mx+1,my+1), sbxy(3,mx+1,my+1)
      double precision sum1, sum2
      real idex,dtc1,ddx,ddy,ddz,qtmh1,qtmh2,p6,p7
      
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
!$OMP& edgery,sum1,sfxy,sbxy,spsit,
!$OMP& ddx,ddy,ddz,qtmh1,qtmh2,p6,p7,dtc1)
!$OMP& REDUCTION(+:sum2)
      do 70 k = 1, mxyp1
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
      mnoff = moffp + noff - 1
! load local fields from global arrays
      do 20 j = 1, mm+1
      do 10 i = 1, nn+1
      sfxy(1,i,j) = fxy(1,i+noffp,j+moffp)
      sfxy(2,i,j) = fxy(2,i+noffp,j+moffp)
   10 continue
   20 continue
      do 40 j = 1, mm+1
      do 30 i = 1, nn+1
      sbxy(1,i,j) = bxy(1,i+noffp,j+moffp)
      sbxy(2,i,j) = bxy(2,i+noffp,j+moffp)
      sbxy(3,i,j) = bxy(3,i+noffp,j+moffp)
   30 continue
   40 continue
      do 46 j = 1, mm+1
      do 42 i = 1, nn+1
      spsit(i,j) = psit(i+noffp,j+moffp)
   42 continue
   46 continue
! clear counters
      do 50 j = 1, 8
      ncl(j,k) = 0
   50 continue
      sum1 = 0.0d0
! loop over particles in tile
      do 60 j = 1, nppp
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
      dx = amx*sfxy(1,nn,mm)
      dy = amx*sfxy(2,nn,mm)
      dz = amx*spsit(nn,mm)
      dx = amy*(dxp*sfxy(1,nn+1,mm) + dx)
      dy = amy*(dxp*sfxy(2,nn+1,mm) + dy)
      dz = amy*(dxp*spsit(nn+1,mm) + dz)
      acx = amx*sfxy(1,nn,mm+1)
      acy = amx*sfxy(2,nn,mm+1)
      acz = amx*spsit(nn,mm+1)
      dx = dx + dyp*(dxp*sfxy(1,nn+1,mm+1) + acx) 
      dy = dy + dyp*(dxp*sfxy(2,nn+1,mm+1) + acy)
      dz = dz + dyp*(dxp*spsit(nn+1,mm+1) + acz)
      
      dz = dz * idex

! find magnetic field
      ox = amx*sbxy(1,nn,mm)
      oy = amx*sbxy(2,nn,mm)
      oz = amx*sbxy(3,nn,mm)
      ox = amy*(dxp*sbxy(1,nn+1,mm) + ox)
      oy = amy*(dxp*sbxy(2,nn+1,mm) + oy)
      oz = amy*(dxp*sbxy(3,nn+1,mm) + oz)
      acx = amx*sbxy(1,nn,mm+1)
      acy = amx*sbxy(2,nn,mm+1)
      acz = amx*sbxy(3,nn,mm+1)
      ox = ox + dyp*(dxp*sbxy(1,nn+1,mm+1) + acx) 
      oy = oy + dyp*(dxp*sbxy(2,nn+1,mm+1) + acy)
      oz = oz + dyp*(dxp*sbxy(3,nn+1,mm+1) + acz)

      ox = ox * dex
      oy = oy * dex

! calculate half impulse
      dx = qtmh2*dx
      dy = qtmh2*dy
      dz = qtmh2*dz
! half acceleration
      acx = ppart(3,j,k) + dx
      acy = ppart(4,j,k) + dy
      acz = ppart(5,j,k) + dz
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
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
      ppart(3,j,k) = dx
      ppart(4,j,k) = dy
      ppart(5,j,k) = dz
! update inverse gamma
      dtc1 = dt/(sqrt(1+(dx*dx+dy*dy+dz*dz)*dex*dex)-dz*dex)

!      p2 = dx*dx + dy*dy + dz*dz
!      dtg = dtc/sqrt(1.0 + p2*ci2)
! new position
      dx = x + dx*dtc1
      dy = y + dy*dtc1
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
      ppart(1,j,k) = dx
      ppart(2,j,k) = dy
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
   60 continue
      sum2 = sum2 + sum1
! set error and end of file flag
! ihole overflow
      if (nh.gt.0) then
         irc = ih
         ih = -ih
      endif
      ihole(1,1,k) = ih
   70 continue
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
      dimension ppart(idimp,nppmx,mx1*myp1)
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
      dx = ppart(1,j,k)
      dy = ppart(2,j,k)
! find particles going out of bounds
      ist = 0
! count how many particles are going in each direction in ncl
! save their address and destination in ihole
! use periodic boundary conditions and check for roundoff error
! ist = direction particle is going
      if (dx.ge.edgerx) then
         if (dx.ge.anx) ppart(1,j,k) = dx - anx
         ist = 2
      else if (dx.lt.edgelx) then
         if (dx.lt.0.0) then
            dx = dx + anx
            if (dx.lt.anx) then
               ist = 1
            else
               dx = 0.0
            endif
            ppart(1,j,k) = dx
         else
            ist = 1
         endif
      endif
      if (dy.ge.edgery) then
         if (dy.ge.any) ppart(2,j,k) = dy - any
         ist = ist + 6
      else if (dy.lt.edgely) then
         if (dy.lt.0.0) then
            dy = dy + any
            if (dy.lt.any) then
               ist = ist + 3
            else
               dy = 0.0
            endif
            ppart(2,j,k) = dy
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
         ppbuff(i,ii,k) = ppart(i,j1,k)
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
      dimension ppart(idimp,nppmx,mx1*myp1)
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
         ppbuff(i,ii,k) = ppart(i,j1,k)
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
      dimension ppart(idimp,nppmx,mx1*myp1)
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
            ppart(i,j1,k) = rbufl(i,j+noff)
  130       continue
         else if (ks(ii).gt.mxyp1) then
            do 140 i = 1, idimp
            ppart(i,j1,k) = rbufr(i,j+moff)
  140       continue
         else
            do 150 i = 1, idimp
            ppart(i,j1,k) = ppbuff(i,j+ncoff,ks(ii))
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
            ppart(i,j2,k) = ppart(i,j1,k)
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
!-deposit psi on particles
      subroutine WPGPSIPOST2L_QP(ppart,psi,kpic,qbm,noff,nyp,idimp,nppmx
     1,nx,mx,my,nxv,nypmx,mx1,mxyp1,dex)

      implicit none
      integer noff, nyp, idimp, nppmx, nx, mx, my, nxv, nypmx
      integer mx1, mxyp1
      real dex,qbm
      real ppart, psi
      integer kpic
      dimension ppart(idimp,nppmx,mxyp1)
      dimension psi(nxv,nypmx)
      dimension kpic(mxyp1)
! local data
!      integer MXV, MYV
!      parameter(MXV=33,MYV=33)
      integer noffp, moffp, nppp
      integer mnoff, i, j, k, nn, mm
      real dxp, dyp, amx, amy, acx
      real dx, dy, dz
      real x, y, vx, vy, vz
      real spsi
      dimension spsi(mx+1,my+1)
      real dx2
!     dimension sfxy(3,mx+1,my+1), sbxy(3,mx+1,my+1)
!     dimension scu(3,mx+1,my+1,), sdcu(3,mx+1,my+1), samu(4,mx+1,my+1)
      dx2 = dex * dex
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV)) return
! loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,noffp,moffp,nppp,nn,mm,mnoff,x,y,vx,vy,
!$OMP& dxp,dyp,amx,amy,dx,acx,spsi)
      do 120 k = 1, mxyp1
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(k)
      mnoff = moffp + noff - 1
! load local fields from global arrays
      nn = min(mx,nx-noffp) + 1
      mm = min(my,nyp-moffp) + 1
      do 20 j = 1, mm
      do 10 i = 1, nn
      spsi(i,j) = psi(i+noffp,j+moffp)
   10 continue
   20 continue
! loop over particles in tile
      do 70 j = 1, nppp
! find interpolation weights
      x = ppart(1,j,k)
      y = ppart(2,j,k)
      vx = ppart(3,j,k)
      vy = ppart(4,j,k)
      nn = x
      mm = y
      dxp = x - real(nn)
      dyp = y - real(mm)
      nn = nn - noffp + 1
      mm = mm - mnoff
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      
      
! find electric field
      dx = amx*spsi(nn,mm)
      dx = amy*(dxp*spsi(nn+1,mm) + dx)
      acx = amx*spsi(nn,mm+1)
      dx = dx + dyp*(dxp*spsi(nn+1,mm+1) + acx) 

      dx = - dx*dx2*qbm
      
      ppart(7,j,k) = 1.0 + dx
      ppart(6,j,k) = dx+(dx2*(vx**2+vy**2)-dx**2-2.*dx)/(2.*(1.+dx))+1.0
   70 continue
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
      dimension part(idimp,npmax), ppart(idimp,nppmx,mxyp1)
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
      part(i,j+npoff) = ppart(i,j,k)
   10 continue
   20 continue
      npoff = npoff + nppp
   30 continue
      npp = npoff
      if (ierr.gt.0) irc = ierr
      return
      end
!-----------------------------------------------------------------------      