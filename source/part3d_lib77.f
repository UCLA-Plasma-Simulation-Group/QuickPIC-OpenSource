c-----------------------------------------------------------------------
      subroutine PRVDIST32_TWISS(part,qm,edges,npp,nps,alpha_x,alpha_y,
     1beta_x,beta_y,emt_x,emt_y,sigz,vdx,vdy,vdz,vtz,npx,npy,npz,idimp ,
     1npmax, nx, ny,nz,x0,y0,z0,mblok,nblok,idps,ierr,gamma,lquiet)
c Cut-off at 3 sigma
c cut-off by a circle
c keep 1 + p^2 = gamma     
c for 3d code, this subroutine calculates initial particle co-ordinates
c and velocities using twiss parameters with 2D spatial decomposition.
c part(1,n,m) = position x of particle n
c part(2,n,m) = position y of particle n
c part(3,n,m) = position z of particle n
c part(4,n,m) = velocity vx of particle n
c part(5,n,m) = velocity vy of particle n
c part(6,n,m) = velocity vz of particle n
c edges(1,m) = lower boundary in y of particle partition m
c edges(2,m) = upper boundary in y of particle partition m
c edges(3,m) = lower boundary in z of particle partition m
c edges(4,m) = upper boundary in z of particle partition m
c alpha, beta = twiss parameters
c emt_x emt_y = normalized emittances
c npx/npy/npz = initial number of particles distributed in x/y/z
c direction
c idimp = size of phase space = 6
c npmax = max number of particles
c nx/ny/nz = system length in x/y/z direction
      implicit none
      include "mpif.h"

c      common /f77_common/ f77_log_unit, f77_output_unit
c      integer f77_log_unit, f77_output_unit
      character(len=60) strMessage

c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mint = default datatype for integers
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld

      real qm,x0,y0,z0,sigz,vdx,vdy,vdz,vtz, gamma
      real part, edges, alpha_x, alpha_y, beta_x, beta_y, emt_x, emt_y
      integer npx,npy,npz,idimp,nx,ny,nz
      integer nps, npp, npmax, mblok, nblok, idps, ierr
      dimension part(idimp,npmax,nblok)
      dimension nps(nblok), npp(nblok)
      dimension edges(idps,nblok)
             
      double precision randum, ranorm
      real sigx,sigy,tvtx,tvty,tvtz,vtx,vty
      integer isum2,iwork2
      double precision sum0, sum1, sum2
      real sum3,work3
      dimension sum3(3), work3(3), isum2(2), iwork2(2)
      integer one
      integer pt,j,k,l, m, mnblok, npxy, npxyz
      double precision r, tempz, tempx, tempy, dnpxyz
      real tx,ty
      real borderlx,borderly,borderlz, borderx, bordery, borderz
      logical lquiet
      integer my,mz,moff,k1,npt
      double precision at1
      real x2,y2

      ierr = 0
      npt = 1
      dnpxyz = npx
      dnpxyz = dnpxyz*npy
      dnpxyz = dnpxyz*npz
      x2 = 2.0 * x0
      y2 = 2.0 * y0
      
      sigx = sqrt(beta_x*emt_x/gamma)
      sigy = sqrt(beta_y*emt_y/gamma)
      vtx = emt_x/sigx
      vty = emt_y/sigy

      borderlx = max((x0-3.0*sigx),1.0)
      borderly = max((y0-3.0*sigy),1.0)
      borderlz = max((z0-3.0*sigz),1.0)
      borderx = min((x0+3.0*sigx),float(nx-1)) 
      bordery = min((y0+3.0*sigy),float(ny-1))
      borderz = min((z0+3.0*sigz),float(nz-1))
      
      j = 0
      l = npz
      do while (j<npx)
      k = 0
      do while (k<npy)
      l = l - npz
      do while (l<npz)

  10    tempz = z0+sigz*ranorm()  
        if (tempz>=(borderz) .or. tempz<=borderlz) goto 10
        
  20    tx = ranorm()
        tempx = x0+sigx*tx
        if (tempx>=(borderx) .or. tempx<=borderlx) goto 20

  30    ty = ranorm()
        tempy = y0+sigy*ty
        if (tempy>=(bordery) .or. tempy<=borderly) goto 30

        if ((tx*tx+ty*ty) > 9.0) goto 20
c  generate velocity 
        tvtx = vtx*ranorm()
        tvty = vty*ranorm()
        tvtz = vtz*ranorm() + vdz

        tvtx = tvtx - gamma*alpha_x/beta_x*(tempx-x0)
        tvty = tvty - gamma*alpha_y/beta_y*(tempy-y0)
!        tvtz = sqrt(tvtz*tvtz-1-tvtx*tvtx-tvty*tvty)
        do 110 mz = 1, nblok
        moff = mblok*(mz - 1)
        do 100 my = 1, mblok
        m = my + moff   
c  check if particle belongs to this partition
        if ((tempy.ge.edges(1,m)) .and. (tempy.lt.edges(2,m)) .and.     &
     &      (tempz.ge.edges(3,m)) .and. (tempz.lt.edges(4,m)) ) then 
          if (npt.le.npmax) then        
            npt = npp(m) + 1
            part(3,npt,m) = tempz 
            part(1,npt,m) = tempx
            part(2,npt,m) = tempy
            part(4,npt,m) = tvtx
            part(5,npt,m) = tvty
            part(6,npt,m) = tvtz  
            part(7,npt,m) = qm
            npp(m) = npt
          else
            ierr = ierr + 1
          endif
c quiet start          
          if (lquiet) then
             if (npt.le.npmax) then
                npt = npp(m) + 1
                part(3,npt,m) = tempz
                part(1,npt,m) = x2 - tempx
                part(2,npt,m) = y2 - tempy
                part(4,npt,m) = -tvtx
                part(5,npt,m) = -tvty
                part(6,npt,m) = tvtz  
                part(7,npt,m) = qm
                npp(m) = npt
             else
                ierr = ierr + 1
             endif
          endif
        endif   
  100   continue        
  110   continue        
        l = l + 1
        if (lquiet) l = l + 1
      enddo     
      k = k + 1
      enddo
      j = j + 1
      enddo
c add correct drift
      sum3(1) = 0.
      sum3(2) = 0.
      sum3(3) = 0.
      mnblok = mblok*nblok
      do 210 m = 1, mnblok
      sum0 = 0.0d0
      sum1 = 0.0d0
      sum2 = 0.0d0
      do 200 j = nps(m), npp(m)
c      npxyz = npxyz + 1
      sum0 = sum0 + part(4,j,m)
      sum1 = sum1 + part(5,j,m)
      sum2 = sum2 + part(6,j,m)
  200 continue
      sum3(1) = sum3(1) + sum0
      sum3(2) = sum3(2) + sum1
      sum3(3) = sum3(3) + sum2
  210 continue
      isum2(1) = ierr
      isum2(2) = 0
c use MPI call instead of PISUM      
      iwork2 = isum2
      call MPI_ALLREDUCE(iwork2,isum2,2,mint,MPI_SUM,lgrp,ierr) 
      ierr = isum2(1)
c process errors
      if (ierr.gt.0) then
         write (2,*) 'particle overflow error, ierr = ', ierr
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PRVDIST32_RAN_PFL(part,qm,edges,npp,nps,x0,y0,z0,sigx,s
     1igy,vtx,vty,vtz,vdx,vdy,vdz,cx,cy,npx,npy,npz,nx,ny,nz,ipbc,idimp,
     2npmax,mblok,nblok,idps,dp,lquiet,ierr)
c new quiet start
c keep 1 + p^2 = gamma
c for 3d code, this subroutine calculates initial particle co-ordinates
c and velocities with 2D spatial decomposition. The longitudinal density 
c profile is described by 1D array dp and transverse profile is gaussian.
c velocity is maxwellian with drift. The method used is acceptance-
c rejection method.
c part(1,n,m) = position x of particle n in partition m
c part(2,n,m) = position y of particle n in partition m
c part(3,n,m) = position z of particle n in partition m
c part(4,n,m) = velocity vx of particle n in partition m
c part(5,n,m) = velocity vy of particle n in partition m
c part(6,n,m) = velocity vz of particle n in partition m
c edges(1,m) = lower boundary in y of particle partition m
c edges(2,m) = upper boundary in y of particle partition m
c edges(3,m) = lower boundary in z of particle partition m
c edges(4,m) = upper boundary in z of particle partition m
c npp(m) = number of particles in partition m
c nps(m) = starting address of particles in partition m
c vtx/vty/vtz = thermal velocity of electrons in x/y/z direction
c vdx/vdy/vdz = drift velocity of beam electrons in x/y/z direction
c npx/npy/npz = initial number of particles distributed in x/y/z
c direction
c nx/ny/nz = system length in x/y/z direction
c idimp = size of phase space = 6
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,xy periodic,xy reflecting, x reflecting/y periodic)
c npmax = maximum number of particles in each partition
c mblok/nblok = number of particle partitions in y/z
c idps = number of partition boundaries
c ierr = (0,1) = (no,yes) error condition exists
c ranorm = gaussian random number with zero mean and unit variance
c with spatial decomposition
      implicit none

c      common /f77_common/ f77_log_unit, f77_output_unit
c      integer f77_log_unit, f77_output_unit
      character(len=60) strMessage

      integer nps,npp,npmax,nblok,npx,npy,npz,idimp,nx,ny,nz,idps,ierr
      integer mblok, ipbc
      real part,edges,x0,y0,z0,sigx,sigy,vtx,vty,vtz,vdx,vdy,vdz,dp
      double precision random,ranorm
      real cx, cy, qm
      real cdth,sdth
      dimension part(idimp,npmax,nblok),dp(nz)
      dimension edges(idps,nblok), npp(nblok), nps(nblok)
      dimension cx(0:2),cy(0:2)
      logical lquiet

c local variables
      real tmpu(1024),tmpz(1024)
      integer xi, yi, zi 
      integer npxyz,isum2,iwork2
      integer j,k,l,m,my,mz,moff,mnblok,k1,npt
      real xf, xs, yf, ys, zf, zs, u, amax, dpi
      double precision sum0, sum1, sum2, dnpxyz
      real at1,sum3,work3,borderlx,borderly, borderx, bordery
      real tempx,tempy,tempxx,tempyy, x2, y2,tempz,tvtx,tvty,tvtz
      real tempx0,tempy0,tvtx0,tvty0
      dimension sum3(3), work3(3), isum2(2), iwork2(2)
c borderlx(yz), lower bound of border, borderx(yz), upper bound.      
      integer nz1 
      integer cnt

      ierr = 0

      cdth = sqrt(2.0)/2.0
      sdth = sqrt(2.0)/2.0
      
      npt = 1
!      npxyz = npx*npy*npz
      dnpxyz = npx
      dnpxyz = dnpxyz*npy
      dnpxyz = dnpxyz*npz
     
      x2 = 2.0 * x0
      y2 = 2.0 * y0
      borderlx = max((x0-5.0*sigx),1.0)
      borderly = max((y0-5.0*sigy),1.0)
      borderx = min((x0+5.0*sigx),float(nx-1)) 
      bordery = min((y0+5.0*sigy),float(ny-1))

      nz1 = nz -1
      dp = abs(dp)
      amax = maxval(dp)*1.2
      dp = dp/amax
      j=0
      cnt = 0
      l = npz
      do while (j<npx)
      k = 0
      do while (k<npy)
      l = l - npz
      do while (l<npz)
  495    if (cnt<=0) then
            call RANDOM_NUMBER(tmpu)
            call RANDOM_NUMBER(tmpz)
            cnt=1024
         endif
         u=tmpu(cnt)
         tempz=tmpz(cnt)*nz1
         cnt=cnt-1
         zi = int(tempz)
         zf = tempz-zi
         zs = 1.-zf
         zi = zi+1
        dpi = dp(zi)*zs + dp(zi+1)*zf    
        if (u<=dpi) then  
c particle is accepted        

  20    tempx = x0+sigx*ranorm()
        if (tempx>=(borderx) .or. tempx<=borderlx) goto 20

  30    tempy = y0+sigy*ranorm()
        if (tempy>=(bordery) .or. tempy<=borderly) goto 30

c  check if particle belongs to this partition
             do m = 1, nblok
              if ((tempy<edges(2,m)) .and. (tempy>edges(1,m)) .and.     &
     & (tempz<edges(4,m)) .and. (tempz>edges(3,m))) then  
                tvtx = vtx*ranorm()
                tvty = vty*ranorm()
                tvtz = vtz*ranorm() + vdz
                tvtz = sqrt(tvtz*tvtz-1-tvtx*tvtx-tvty*tvty)
                if (npt.le.npmax) then        
                     npt = npp(m) + 1
                     part(3,npt,m) = tempz
                     tempxx = -cx(0)*(part(3,npt,m)-z0)**2-cx(1)*(part(3&
     &,npt,m)-z0)-cx(2)                
                     part(1,npt,m) = tempx + tempxx
                     tempyy = -cy(0)*(part(3,npt,m)-z0)**2-cy(1)*(part(3&
     &,npt,m)-z0)-cy(2)        
                     part(2,npt,m) = tempy + tempyy 
                     part(4,npt,m) = tvtx
                     part(5,npt,m) = tvty
                     part(6,npt,m) = tvtz  
                     part(7,npt,m) = qm
                     npp(m) = npt
                else
                     ierr = ierr + 1
                endif     
c quiet start
                if (lquiet) then
                   if (npt.le.npmax) then
                      npt = npp(m) + 1
                      part(3,npt,m) = tempz
                      part(1,npt,m) = x2 - tempx+tempxx
                      part(2,npt,m) = y2 - tempy+tempyy
                      part(4,npt,m) = -tvtx
                      part(5,npt,m) = -tvty
                      part(6,npt,m) = tvtz  
                      part(7,npt,m) = qm
                      npp(m) = npt
                   else
                      ierr = ierr + 1
                   endif
                endif
              endif
             enddo 
        l = l + 1
        if (lquiet) l = l + 1
        endif
      enddo 
      k = k + 1
      enddo
      j = j + 1
      enddo

      write (2,*) "part. gen. done, now add drift" 
      
c      npxyz = 0
c add correct drift
      sum3(1) = 0.
      sum3(2) = 0.
      sum3(3) = 0.
      mnblok = mblok*nblok
      do 260 m = 1, mnblok
      sum0 = 0.0d0
      sum1 = 0.0d0
      sum2 = 0.0d0
      do 250 j = nps(m), npp(m)
c      npxyz = npxyz + 1
      sum0 = sum0 + part(4,j,m)
      sum1 = sum1 + part(5,j,m)
      sum2 = sum2 + part(6,j,m)
  250 continue
      sum3(1) = sum3(1) + sum0
      sum3(2) = sum3(2) + sum1
      sum3(3) = sum3(3) + sum2
  260 continue
      isum2(1) = ierr
c      isum2(2) = npxyz
c      call PISUM(isum2,iwork2,2,1)
      ierr = isum2(1)
c      npxyz = isum2(2)
c      call PSUM(sum3,work3,3,1)
c      at1 = 1./float(npxyz)
c      sum3(1) = at1*sum3(1) - vdx
c      sum3(2) = at1*sum3(2) - vdy
c      sum3(3) = at1*sum3(3) - vdz
c      do 280 m = 1, nblok
c      do 270 j = nps(m), npp(m)
c      part(4,j,m) = part(4,j,m) - sum3(1)
c      part(5,j,m) = part(5,j,m) - sum3(2)
c      part(6,j,m) = part(6,j,m) - sum3(3)
c  270 continue
c  280 continue
c process errors
      if (ierr.gt.0) then
         write (2,*) 'particle overflow error, ierr = ', ierr
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PRVDIST32_RANDOM(part,qm,edges,npp,nps,vtx,vty,vtz,vdx,
     1vdy,vdz,npx,npy,npz,nx,ny,nz,ipbc,idimp,npmax,mblok,nblok,idps,sig
     1x,sigy,sigz,x0,y0,z0,cx,cy,lquiet,ierr)
c old quiet start     
c keep 1 + p^2 = gamma
c for 3d code, this subroutine calculates initial particle co-ordinates
c and velocities with tri-gaussian density and maxwellian velocity with 
c drift for distributed data with 2D spatial decomposition
c part(1,n,m) = position x of particle n in partition m
c part(2,n,m) = position y of particle n in partition m
c part(3,n,m) = position z of particle n in partition m
c part(4,n,m) = velocity vx of particle n in partition m
c part(5,n,m) = velocity vy of particle n in partition m
c part(6,n,m) = velocity vz of particle n in partition m
c edges(1,m) = lower boundary in y of particle partition m
c edges(2,m) = upper boundary in y of particle partition m
c edges(3,m) = lower boundary in z of particle partition m
c edges(4,m) = upper boundary in z of particle partition m
c npp(m) = number of particles in partition m
c nps(m) = starting address of particles in partition m
c vtx/vty/vtz = thermal velocity of electrons in x/y/z direction
c vdx/vdy/vdz = drift velocity of beam electrons in x/y/z direction
c npx/npy/npz = initial number of particles distributed in x/y/z
c direction
c nx/ny/nz = system length in x/y/z direction
c idimp = size of phase space = 6
c npmax = maximum number of particles in each partition
c mblok/nblok = number of particle partitions in y/z
c idps = number of partition boundaries
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,xy periodic,xy reflecting, x reflecting/y periodic)
c ierr = (0,1) = (no,yes) error condition exists
c ranorm = gaussian random number with zero mean and unit variance
      implicit none

      integer nps,npp,npmax,nblok,npx,npy,npz,idimp,nx,ny,nz,idps,ierr
      integer mblok,ipbc
      real qm,sigx,sigy,sigz,x0,y0,z0
      real cx,cy
      double precision random,ranorm
      real part,edges,vtx,vty,vtz,vdx,vdy,vdz
      dimension part(idimp,npmax,nblok)
      dimension edges(idps,nblok), npp(nblok), nps(nblok)
      dimension cx(0:2),cy(0:2)
      logical lquiet

c local variables
      integer j,k,l,m,my,mz,moff,mnblok,k1,npt
      real tempx,tempy,tempxx,tempyy,x2,y2,tempz,tvtx,tvty,tvtz
c borderlx(yz), lower bound; borderx(yz), upper bound.      
      real borderlx,borderly,borderlz, borderx, bordery, borderz 

      ierr = 0
      
      npt = 1

      x2 = 2.0 * x0
      y2 = 2.0 * y0
      borderlx = max((x0-5.0*sigx),1.0)
      borderly = max((y0-5.0*sigy),1.0)
      borderlz = max((z0-5.0*sigz),1.0)
      borderx = min((x0+5.0*sigx),float(nx-1)) 
      bordery = min((y0+5.0*sigy),float(ny-1))
      borderz = min((z0+5.0*sigz),float(nz-1))
      

      j = 0
      l = npz
      do while (j<npx)
      k = 0
      do while (k<npy)
      l = l - npz
      do while (l<npz)
  10    tempz = z0+sigz*ranorm()  
        if (tempz>=(borderz) .or. tempz<=borderlz) goto 10
        
  20    tempx = x0+sigx*ranorm()
        if (tempx>=(borderx) .or. tempx<=borderlx) goto 20

  30    tempy = y0+sigy*ranorm()
        if (tempy>=(bordery) .or. tempy<=borderly) goto 30
           
c  generate velocity 
        tvtx = vtx*ranorm() + vdx
        tvty = vty*ranorm() + vdy
        tvtz = vtz*ranorm() + vdz
        tvtz = sqrt(tvtz*tvtz-1-tvtx*tvtx-tvty*tvty)
        do 110 mz = 1, nblok
        moff = mblok*(mz - 1)
        do 100 my = 1, mblok
        m = my + moff   
c  check if particle belongs to this partition
        if ((tempy.ge.edges(1,m)) .and. (tempy.lt.edges(2,m)) .and.     &
     &      (tempz.ge.edges(3,m)) .and. (tempz.lt.edges(4,m)) ) then 
          if (npt.lt.npmax) then        
            npt = npp(m) + 1
            part(3,npt,m) = tempz 
c calculate offset in x            
            tempxx = -cx(0)*(part(3,npt,m)-z0)**2-cx(1)*(part(3,npt,m)  &
     & - z0)-cx(2)                
            part(1,npt,m) = tempx + tempxx
c calculate offset in y            
            tempyy = -cy(0)*(part(3,npt,m)-z0)**2-cy(1)*(part(3,npt,m)  &
     & - z0)-cy(2)        
            part(2,npt,m) = tempy + tempyy 
            part(4,npt,m) = tvtx
            part(5,npt,m) = tvty
            part(6,npt,m) = tvtz 
            part(7,npt,m) = qm
            npp(m) = npt
          else
            ierr = ierr + 1
          endif
c quiet start          
          if (lquiet) then
             if (npt.lt.npmax) then
                npt = npp(m) + 1
                part(3,npt,m) = tempz
                part(1,npt,m) = x2 - tempx + tempxx
                part(2,npt,m) = y2 - tempy + tempyy
                part(4,npt,m) = -tvtx
                part(5,npt,m) = -tvty
                part(6,npt,m) = tvtz 
                part(7,npt,m) = qm
                npp(m) = npt
             else
                ierr = ierr + 1
             endif
          endif
        endif   
  100   continue        
  110   continue        
        l = l + 1
        if (lquiet) l = l + 1
      enddo     
      k = k + 1
      enddo
      j = j + 1
      enddo
      
      return
      end
c-----------------------------------------------------------------------
      subroutine PGPOST32L(part,q,npp,noff,idimp,npmax,mnblok,nxv,nypmx,
     1nzpmx,idds)
c for 3d code, this subroutine calculates particle charge density
c using first-order linear interpolation, and distributed data
c with 2D spatial decomposition
c scalar version using guard cells, for distributed data
c 33 flops/particle, 11 loads, 8 stores
c input: all, output: q
c charge density is approximated by values at the nearest grid points
c q(n,m,l)=qm*(1.-dx)*(1.-dy)*(1.-dz)
c q(n+1,m,l)=qm*dx*(1.-dy)*(1.-dz)
c q(n,m+1,l)=qm*(1.-dx)*dy*(1.-dz)
c q(n+1,m+1,l)=qm*dx*dy*(1.-dz)
c q(n,m,l+1)=qm*(1.-dx)*(1.-dy)*dz
c q(n+1,m,l+1)=qm*dx*(1.-dy)*dz
c q(n,m+1,l+1)=qm*(1.-dx)*dy*dz
c q(n+1,m+1,l+1)=qm*dx*dy*dz
c where n,m,l = nearest grid points and dx = x-n, dy = y-m, dz = z-l
c part(1,n,m) = position x of particle n in partition m
c part(2,n,m) = position y of particle n in partition m
c part(3,n,m) = position z of particle n in partition m
c q(j,k,l,m) = charge density at grid point (j,kk,ll),
c where kk = k + noff(1,m) - 1, and ll = l + noff(2,m) - 1
c npp(m) = number of particles in partition m
c noff(1,m) = lowermost global gridpoint in y in particle partition m
c noff(2,m) = backmost global gridpoint in z in particle partition m
c qm = charge on particle, in units of e
c idimp = size of phase space = 6
c npmax = maximum number of particles in each partition
c mnblok = number of particle partitions.
c nxv = first dimension of charge array, must be >= nx+1
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c idds = dimensionality of domain decomposition
      dimension part(idimp,npmax,mnblok), q(nxv,nypmx,nzpmx,mnblok)
      dimension npp(mnblok), noff(idds,mnblok)
      real qm
      do 20 m = 1, mnblok
      mnoff = noff(1,m) - 1
      lnoff = noff(2,m) - 1
c find interpolation weights
      do 10 j = 1, npp(m)
      nn = part(1,j,m)
      mm = part(2,j,m)
      ll = part(3,j,m)
      qm = part(7,j,m)
      dxp = qm*(part(1,j,m) - float(nn))
      dyp = part(2,j,m) - float(mm)
      dzp = part(3,j,m) - float(ll)
      nn = nn + 1
      amx = qm - dxp
      amy = 1. - dyp
      np = nn + 1
      mm = mm - mnoff
      dx1 = dxp*dyp
      dyp = amx*dyp
      mp = mm + 1
      amx = amx*amy
      amz = 1. - dzp
      ll = ll - lnoff
      amy = dxp*amy
      lp = ll + 1
c deposit charge
      q(nn,mm,ll,m) = q(nn,mm,ll,m) + amx*amz
      q(np,mm,ll,m) = q(np,mm,ll,m) + amy*amz
      q(nn,mp,ll,m) = q(nn,mp,ll,m) + dyp*amz
      q(np,mp,ll,m) = q(np,mp,ll,m) + dx1*amz
      q(nn,mm,lp,m) = q(nn,mm,lp,m) + amx*dzp
      q(np,mm,lp,m) = q(np,mm,lp,m) + amy*dzp
      q(nn,mp,lp,m) = q(nn,mp,lp,m) + dyp*dzp
      q(np,mp,lp,m) = q(np,mp,lp,m) + dx1*dzp
   10 continue
   20 continue
      return
      end      
c-----------------------------------------------------------------------
      subroutine PGBPUSH32L_QP(part,fxyz,bxyz,npp,noff,qbm,dt,dtc,ek,nx,
     1ny,nz,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc,deltax,deltaz ,
     2cofd)

      double precision sum1
      dimension part(idimp,npmax,mnblok)
      dimension fxyz(3,nxv,nypmx,nzpmx,mnblok)
      dimension bxyz(3,nxv,nypmx,nzpmx,mnblok)
      dimension npp(mnblok), noff(idds,mnblok)

      real qtmg,dtgx,dtgy,cofd
      real delx, dely, delz, ngamma
      real p2t,p2l, cofd1 
      real dtc_over_deltax, dtc_over_deltaz
      real one_minus_vz0, one_minus_vz

      qtmh = qbm*dt
      sum1 = 0.0d0
      dtc_over_deltax = dtc/deltax
      dtc_over_deltaz = dtc/deltaz

      one_minus_vz0 = 0. 
      cofd1 = cofd/dt
c set boundary values
      if (ipbc.eq.1) then
         edgelx = 1.
         edgely = 1.
         edgelz = 1.
         edgerx = float(nx-1)
         edgery = float(ny-1)
         edgerz = float(nz-1)
      endif

      do 20 m = 1, mnblok
      mnoff = noff(1,m) - 1
      lnoff = noff(2,m) - 1
c find interpolation weights
      inpp = npp(m)
      do 10 j = 1, npp(m)
      if (j.gt.inpp) exit
   11 nn = part(1,j,m)
      mm = part(2,j,m)
      ll = part(3,j,m)
      dxp = part(1,j,m) - float(nn)
      dyp = part(2,j,m) - float(mm)
      dzp = part(3,j,m) - float(ll)
      nn = nn + 1
      amx = 1. - dxp
      amy = 1. - dyp
      np = nn + 1
      mm = mm - mnoff
      dx1 = dxp*dyp
      dyp = amx*dyp
      mp = mm + 1
      amx = amx*amy
      amz = 1. - dzp
      ll = ll - lnoff
      amy = dxp*amy
      lp = ll + 1
c find electric field
      dx = amz*(amx*fxyz(1,nn,mm,ll,m) + amy*fxyz(1,np,mm,ll,m) + dyp*fx
     1yz(1,nn,mp,ll,m) + dx1*fxyz(1,np,mp,ll,m)) + dzp*(amx*fxyz(1,nn,mm
     2,lp,m) + amy*fxyz(1,np,mm,lp,m) + dyp*fxyz(1,nn,mp,lp,m) + dx1*fxy
     3z(1,np,mp,lp,m))
      dy = amz*(amx*fxyz(2,nn,mm,ll,m) + amy*fxyz(2,np,mm,ll,m) + dyp*fx
     1yz(2,nn,mp,ll,m) + dx1*fxyz(2,np,mp,ll,m)) + dzp*(amx*fxyz(2,nn,mm
     2,lp,m) + amy*fxyz(2,np,mm,lp,m) + dyp*fxyz(2,nn,mp,lp,m) + dx1*fxy
     3z(2,np,mp,lp,m))
      dz = amz*(amx*fxyz(3,nn,mm,ll,m) + amy*fxyz(3,np,mm,ll,m) + dyp*fx
     1yz(3,nn,mp,ll,m) + dx1*fxyz(3,np,mp,ll,m)) + dzp*(amx*fxyz(3,nn,mm
     2,lp,m) + amy*fxyz(3,np,mm,lp,m) + dyp*fxyz(3,nn,mp,lp,m) + dx1*fxy
     3z(3,np,mp,lp,m))
c find magnetic field
      ox = amz*(amx*bxyz(1,nn,mm,ll,m) + amy*bxyz(1,np,mm,ll,m) + dyp*bx
     1yz(1,nn,mp,ll,m) + dx1*bxyz(1,np,mp,ll,m)) + dzp*(amx*bxyz(1,nn,mm
     2,lp,m) + amy*bxyz(1,np,mm,lp,m) + dyp*bxyz(1,nn,mp,lp,m) + dx1*bxy
     3z(1,np,mp,lp,m))
      oy = amz*(amx*bxyz(2,nn,mm,ll,m) + amy*bxyz(2,np,mm,ll,m) + dyp*bx
     1yz(2,nn,mp,ll,m) + dx1*bxyz(2,np,mp,ll,m)) + dzp*(amx*bxyz(2,nn,mm
     2,lp,m) + amy*bxyz(2,np,mm,lp,m) + dyp*bxyz(2,nn,mp,lp,m) + dx1*bxy
     3z(2,np,mp,lp,m))
     
      dx = dx - oy
      dy = dy + ox

      delx = qtmh*dx
      dely = qtmh*dy
      delz = qtmh*dz
c half acceleration
c momentums are normalized to c temporarily
      dx = part(4,j,m) + delx
      dy = part(5,j,m) + dely
      dz = part(6,j,m) + delz
      part(4,j,m) = dx
      part(5,j,m) = dy
c update gamma and inverse gamma
      p2t = dx*dx + dy*dy
      p2l = dz*dz
      p2 = p2t + p2l 
      ngamma = sqrt(1.0 + p2)
      dtgx = dtc_over_deltax/ngamma
      dtgy = dtc_over_deltax/ngamma
      one_minus_vz = (1.0+p2t)/(dz*(dz+ngamma))

c include radiation damping in the z direction      
c note: this is not time-centered
c      print *,"cofd1,ngamma,dz,delx,dely=",cofd1,ngamma, dz, delx, dely,&
c     &cofd1*ngamma*((delx*delx+dely*dely)*dz-delz*(delx*dx+dely*dy))  
      dz = dz - cofd1*ngamma*((delx*delx+dely*dely)*dz-delz*(delx*dx+del&
     &y*dy))  
c comment this line to shut off longitudinal push     
      part(6,j,m) = dz
c new position in grid unit
      dx = part(1,j,m) + dx*dtgx
      dy = part(2,j,m) + dy*dtgy
c      dz = part(3,j,m) + (one_minus_vz - one_minus_vz0)*dtc_over_deltaz 
      dz = part(3,j,m) + one_minus_vz*dtc_over_deltaz 
c dropping boundary conditions in x and y
      if (ipbc.eq.1) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            if (j.eq.inpp) then
               inpp = inpp - 1
               exit
            end if
            part(:,j,m) = part(:,inpp,m)
            inpp = inpp -1
            goto 11
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            if (j.eq.inpp) then
               inpp = inpp - 1
               exit
            end if
            part(:,j,m) = part(:,inpp,m)
            inpp = inpp -1
            goto 11
         endif
         if ((dz.lt.edgelz).or.(dz.ge.edgerz)) then
            if (j.eq.inpp) then
               inpp = inpp - 1
               exit
            end if
            part(:,j,m) = part(:,inpp,m)
            inpp = inpp -1
            goto 11
         endif
      endif
c set new position
      part(1,j,m) = dx
      part(2,j,m) = dy
c comment this line to shut off longitudinal push     
      part(3,j,m) = dz
   10 continue
      npp(m) = inpp
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PMOVE32(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,pb
     1uff,jsr,jsl,jss,ny,nz,kstrt,nvpy,nvpz,idimp,npmax,mblok,nblok,idps
     2,nbmax,idds,ntmax,tag1,tag2,id,info)
c modified version with particle moving in y only    
c this subroutine moves particles into appropriate spatial regions
c periodic boundary conditions with 2D spatial decomposition
c part(1,n,m) = position x of particle n in partition m
c part(2,n,m) = position y of particle n in partition m
c part(3,n,m) = position z of particle n in partition m
c part(4,n,m) = velocity vx of particle n in partition m
c part(5,n,m) = velocity vy of particle n in partition m
c part(6,n,m) = velocity vz of particle n in partition m
c edges(1,m) = lower boundary in y of particle partition m
c edges(2,m) = upper boundary in y of particle partition m
c edges(3,m) = back boundary in z of particle partition m
c edges(4,m) = front boundary in z of particle partition m
c npp(m) = number of particles in partition m
c sbufl = buffer for particles being sent to back processor
c sbufr = buffer for particles being sent to front processor
c rbufl = buffer for particles being received from back processor
c rbufr = buffer for particles being received from front processor
c pbuff = buffer for particles being sent to next pipeline stage
c ihole = location of holes left in particle arrays
c jsl(idds,m) = number of particles going back in particle partition m
c jsr(idds,m) = number of particles going front in particle partition m
c jss(idds,m) = scratch array for particle partition m
c ny/nz = system length in y/z direction
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c idimp = size of phase space = 6
c npmax = maximum number of particles in each partition
c mblok/nblok = number of particle partitions in y/z
c idps = number of particle partition boundaries
c nbmax =  size of buffers for passing particles between processors
c idds = dimensionality of domain decomposition
c ntmax =  size of hole array for particles leaving processors
c info = status information
c info(1) = ierr = (0,N) = (no,yes) error condition exists
c info(2) = maximum number of particles per processor
c info(3) = minimum number of particles per processor
c info(4) = maximum number of buffer overflows in y
c info(5) = maximum number of buffer overflows in z
c info(6) = maximum number of particle passes required in y
c info(7) = maximum number of particle passes required in z
c info(8) = total number of particles on entry
c info(9) = difference of total number of particles on exit
      implicit none
      include "mpif.h"
      
      real part, edges, sbufr, sbufl, rbufr, rbufl, pbuff
      integer npp, ihole, jsr, jsl, jss, info
      integer ny, nz, kstrt, nvpy, nvpz, idimp, npmax, mblok, nblok
      integer idps, nbmax, idds, ntmax, tag1, tag2, id
      dimension part(idimp,npmax,mblok*nblok)
      dimension edges(idps,mblok*nblok), npp(mblok*nblok)
      dimension sbufl(idimp,nbmax,mblok*nblok)
      dimension sbufr(idimp,nbmax,mblok*nblok)
      dimension rbufl(idimp,nbmax,mblok*nblok)
      dimension rbufr(idimp,nbmax,mblok*nblok)
      dimension pbuff(idimp,nbmax)
      dimension jsl(idds,mblok*nblok), jsr(idds,mblok*nblok)
      dimension jss(idds,mblok*nblok)
      dimension ihole(ntmax,mblok*nblok)
      dimension info(9)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mint = default datatype for integers
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer iy, iz
      parameter(iy=2,iz=3)
      integer ierr, ic, js, ks, mnblok, i, n, m, my, mz, moff, nvp, iter
      integer npr, nps, npt, kb, kl, kr, j, j1, j2, nbsize, nter, mter
      integer itermax
      integer msid, istatus
      integer ibflg, iwork
      double precision bflg, work
      real an, xt
      dimension msid(4), istatus(lstat)
      dimension ibflg(4), iwork(4)
      dimension bflg(2), work(2)
      dimension kb(2)
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      mnblok = mblok*nblok
      nbsize = idimp*nbmax
      do 5 j = 1, 9
      info(j) = 0
    5 continue
c buffer outgoing particles, first in y then in z direction
      ic = iz
      nvp = nvpy*nvpz
      an = float(nz)
      n = 2

      kl = kstrt - nvpy
      if (kl.lt.1) go to 21
      call MPI_IRECV(rbufl,nbsize,mreal,kl-1,tag1,lworld,msid(1),ierr)
      call MPI_WAIT(msid(1),istatus,ierr)
      call MPI_GET_COUNT(istatus,mreal,nps,ierr)
      jsl(2,1) = nps/idimp
      jss(2,1) = npp(1) + jsl(2,1)
      if (jss(2,1).le.npmax) then
         do 11 j = 1,jsl(2,1)
         do 8  i = 1, idimp
         part(i,j+npp(1),1) = rbufl(i,j,1)
    8    continue
   11    continue
         npp(1) = jss(2,1)         
      else
         write (2,*) 'particle overflow', jss(2,1)
         info(1) = jss(2,1)
         return
      endif
      
   21 iter = 2
      nter = 0
      mter = 0
      do 61 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 51 my = 1, mblok
      m = my + moff
      jsl(1,m) = 0
      jsr(1,m) = 0
      jss(2,m) = 0
      do 31 j = 1, npp(m)
      xt = part(ic,j,m)
c particles going backward, not going to happen
      if (xt.lt.edges(2*n-1,m)) then
         jss(2,m) = 1
         write (2,*) 'Error: particles move to the previous stage'
         go to 41
c particles going forward
      else if (xt.ge.edges(2*n,m)) then
         if (jsr(1,m).lt.nbmax) then
            jsr(1,m) = jsr(1,m) + 1
            do 28 i = 1, idimp
            pbuff(i,jsr(1,m)) = part(i,j,m)
   28       continue
            ihole(jsr(1,m),m) = j
         else
            jss(2,m) = 1
            go to 41
         endif
      endif
   31 continue
   41 jss(1,m) = jsl(1,m) + jsr(1,m)
   51 continue
   61 continue
c check for full buffer condition
      nps = 0
      do 101 m = 1, mnblok
      nps = max0(nps,jss(2,m))
  101 continue
      if (nps.gt.0) then
         info(1) = nps
         write (2,*) 'particle buffer full'
         return
      endif
c copy particle buffers
  111 iter = iter + 2
      mter = mter + 1
      do 151 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 141 my = 1, mblok
      m = my + moff
      kr = kstrt + nvpy
      if (kr.gt.nvp) go to 155
c send particles
      call MPI_ISEND(pbuff,idimp*jsr(1,m),mreal,kr-1,tag2,lworld,id,ierr
     1)
  141 continue
  151 continue
c fill the holes
  155 npt = npp(1)
      do 291 j = 1, jsr(1,1)
      j1 = ihole((jsr(1,1)-j+1),1)
      if (j1.lt.npt) then
         do 274 i = 1, idimp
         part(i,j1,1) = part(i,npt,m)
  274    continue
         npt = npt - 1
      else
         npt = npt - 1
      endif
  291 continue
      npp(1) = npt
      itermax = 20000
c buffer outgoing particles, first in y then in z direction
      do 300 n = 1, 1
      if (n.eq.1) then
         ic = iy
         nvp = nvpy
         an = float(ny)
      elseif (n.eq.2) then
         ic = iz
         nvp = nvpz
         an = float(nz)
      endif
      iter = 2
      nter = 0
   20 mter = 0
      do 60 mz = 1, nblok
      moff = mblok*(mz - 1)
      kb(2) = mz + ks
      do 50 my = 1, mblok
      m = my + moff
      kb(1) = my + js
      jsl(1,m) = 0
      jsr(1,m) = 0
      jss(2,m) = 0
      do 30 j = 1, npp(m)
      xt = part(ic,j,m)
c particles going down or backward
      if (xt.lt.edges(2*n-1,m)) then
         if (jsl(1,m).lt.nbmax) then
            jsl(1,m) = jsl(1,m) + 1
            if (kb(n).eq.0) xt = xt + an
            do 23 i = 1, idimp
            sbufl(i,jsl(1,m),m) = part(i,j,m)
   23       continue
            sbufl(ic,jsl(1,m),m) = xt
            ihole(jsl(1,m)+jsr(1,m),m) = j
         else
            jss(2,m) = 1
            go to 40
         endif
c particles going up or forward
      else if (xt.ge.edges(2*n,m)) then
         if (jsr(1,m).lt.nbmax) then
            jsr(1,m) = jsr(1,m) + 1
            if ((kb(n)+1).eq.nvp) xt = xt - an
            do 27 i = 1, idimp
            sbufr(i,jsr(1,m),m) = part(i,j,m)
   27       continue
            sbufr(ic,jsr(1,m),m) = xt
            ihole(jsl(1,m)+jsr(1,m),m) = j
         else
            jss(2,m) = 1
            go to 40
         endif
      endif
   30 continue
   40 jss(1,m) = jsl(1,m) + jsr(1,m)
   50 continue
   60 continue
c check for full buffer condition
      nps = 0
      do 100 m = 1, mnblok
      nps = max0(nps,jss(2,m))
  100 continue
      ibflg(3) = nps
c copy particle buffers
  110 iter = iter + 2
      mter = mter + 1
      do 150 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 140 my = 1, mblok
      m = my + moff
      kb(1) = my + js
      kb(2) = mz + ks
c get particles from below and above or back and front
      kl = kb(n)
      kb(n) = kl + 1
      if (kb(n).ge.nvp) kb(n) = kb(n) - nvp
      kr = kb(1) + nvpy*kb(2) + 1
      kb(n) = kl - 1
      if (kb(n).lt.0) kb(n) = kb(n) + nvp
      kl = kb(1) + nvpy*kb(2) + 1
c this segment is used for shared memory computers
c     jsl(2,m) = jsr(1,kl)
c     do 120 j = 1, jsl(2,m)
c     do 115 i = 1, idimp
c     rbufl(i,j,m) = sbufr(i,j,kl)
c 115 continue
c 120 continue
c     jsr(2,m) = jsl(1,kr)
c     do 130 j = 1, jsr(2,m)
c     do 125 i = 1, idimp
c     rbufr(i,j,m) = sbufl(i,j,kr)
c 125 continue
c 130 continue
c this segment is used for mpi computers
c post receive
      call MPI_IRECV(rbufl,nbsize,mreal,kl-1,iter-1,lworld,msid(1),ierr)
      call MPI_IRECV(rbufr,nbsize,mreal,kr-1,iter,lworld,msid(2),ierr)
c send particles
      call MPI_ISEND(sbufr,idimp*jsr(1,m),mreal,kr-1,iter-1,lworld,msid(
     13),ierr)
      call MPI_ISEND(sbufl,idimp*jsl(1,m),mreal,kl-1,iter,lworld,msid(4)
     1,ierr)
c wait for particles to arrive
      call MPI_WAIT(msid(1),istatus,ierr)
      call MPI_GET_COUNT(istatus,mreal,nps,ierr)
      jsl(2,m) = nps/idimp
      call MPI_WAIT(msid(2),istatus,ierr)
      call MPI_GET_COUNT(istatus,mreal,nps,ierr)
      jsr(2,m) = nps/idimp
  140 continue
  150 continue
c check if particles must be passed further
      nps = 0
      do 180 m = 1, mnblok
c check if any particles coming from above or front belong here
      jsl(1,m) = 0
      jsr(1,m) = 0
      jss(2,m) = 0
      do 160 j = 1, jsr(2,m)
      if (rbufr(ic,j,m).lt.edges(2*n-1,m)) jsl(1,m) = jsl(1,m) + 1
      if (rbufr(ic,j,m).ge.edges(2*n,m)) jsr(1,m) = jsr(1,m) + 1
  160 continue
      if (jsr(1,m).ne.0) then
         if (n.eq.1) then
            write (2,*) 'Info:',jsr(1,m),' particles returning above'
         elseif (n.eq.2) then
            write (2,*) 'Info:',jsr(1,m),' particles returning front'
         endif
      endif
c check if any particles coming from below or back belong here
      do 170 j = 1, jsl(2,m)
      if (rbufl(ic,j,m).ge.edges(2*n,m)) jsr(1,m) = jsr(1,m) + 1
      if (rbufl(ic,j,m).lt.edges(2*n-1,m)) jss(2,m) = jss(2,m) + 1
  170 continue
      if (jss(2,m).ne.0) then
         if (n.eq.1) then
            write (2,*) 'Info:',jss(2,m),' particles returning below'
         elseif (n.eq.2) then
            write (2,*) 'Info:',jss(2,m),' particles returning back'
         endif
      endif
      jsl(1,m) = jsl(1,m) + jss(2,m)
      nps = max0(nps,jsl(1,m)+jsr(1,m))
  180 continue
      ibflg(2) = nps
c make sure sbufr and sbufl have been sent
      call MPI_WAIT(msid(3),istatus,ierr)
      call MPI_WAIT(msid(4),istatus,ierr)
      if (nps.eq.0) go to 240
c remove particles which do not belong here
      do 230 mz = 1, nblok
      moff = mblok*(mz - 1)
      kb(2) = mz + ks
      do 220 my = 1, mblok
      m = my + moff
      kb(1) = my + js
c first check particles coming from above or front
      jsl(1,m) = 0
      jsr(1,m) = 0
      jss(2,m) = 0
      do 190 j = 1, jsr(2,m)
      xt = rbufr(ic,j,m)
c particles going down or back
      if (xt.lt.edges(2*n-1,m)) then
         jsl(1,m) = jsl(1,m) + 1
         if (kb(n).eq.0) xt = xt + an
         rbufr(ic,j,m) = xt
         do 183 i = 1, idimp
         sbufl(i,jsl(1,m),m) = rbufr(i,j,m)
  183    continue
c particles going up or front, should not happen
      elseif (xt.ge.edges(2*n,m)) then
         jsr(1,m) = jsr(1,m) + 1
         if ((kb(n)+1).eq.nvp) xt = xt - an
         rbufr(ic,j,m) = xt
         do 185 i = 1, idimp
         sbufr(i,jsr(1,m),m) = rbufr(i,j,m)
  185    continue
c particles staying here
      else
         jss(2,m) = jss(2,m) + 1
         do 187 i = 1, idimp
         rbufr(i,jss(2,m),m) = rbufr(i,j,m)
  187    continue
      endif
  190 continue
      jsr(2,m) = jss(2,m)
c next check particles coming from below or back
      jss(2,m) = 0
      do 200 j = 1, jsl(2,m)
      xt = rbufl(ic,j,m)
c particles going up or front
      if (xt.ge.edges(2*n,m)) then
         if (jsr(1,m).lt.nbmax) then
            jsr(1,m) = jsr(1,m) + 1
            if ((kb(n)+1).eq.nvp) xt = xt - an
            rbufl(ic,j,m) = xt
            do 193 i = 1, idimp
            sbufr(i,jsr(1,m),m) = rbufl(i,j,m)
  193       continue 
         else
            jss(2,m) = 2*npmax
            go to 210
         endif
c particles going down back, should not happen
      elseif (xt.lt.edges(2*n-1,m)) then
         if (jsl(1,m).lt.nbmax) then
            jsl(1,m) = jsl(1,m) + 1
            if (kb(n).eq.0) xt = xt + an
            rbufl(ic,j,m) = xt
            do 195 i = 1, idimp
            sbufl(i,jsl(1,m),m) = rbufl(i,j,m)
  195       continue
         else
            jss(2,m) = 2*npmax
            go to 210
         endif
c particles staying here
      else
         jss(2,m) = jss(2,m) + 1
         do 197 i = 1, idimp
         rbufl(i,jss(2,m),m) = rbufl(i,j,m)
  197    continue
      endif
  200 continue
  210 jsl(2,m) = jss(2,m)
  220 continue
  230 continue
c check if move would overflow particle array
  240 nps = 0
      npt = npmax
      do 250 m = 1, mnblok
      jss(2,m) = npp(m) + jsl(2,m) + jsr(2,m) - jss(1,m)
      nps = max0(nps,jss(2,m))
      npt = min0(npt,jss(2,m))
  250 continue
      ibflg(1) = nps
      ibflg(4) = -npt
      iwork = ibflg
      call MPI_ALLREDUCE(iwork,ibflg,4,mint,MPI_MAX,lgrp,ierr) 
      info(2) = ibflg(1)
      info(3) = -ibflg(4)
      ierr = ibflg(1) - npmax
      if (ierr.gt.0) then
         write (2,*) 'particle overflow error, ierr = ', ierr
         info(1) = ierr
         return
      endif
c distribute incoming particles from buffers
      do 290 m = 1, mnblok
c distribute particles coming from below or back into holes
      jss(2,m) = min0(jss(1,m),jsl(2,m))
      do 260 j = 1, jss(2,m)
      do 255 i = 1, idimp
      part(i,ihole(j,m),m) = rbufl(i,j,m)
  255 continue
  260 continue
      if (jss(1,m).gt.jsl(2,m)) then
         jss(2,m) = min0(jss(1,m)-jsl(2,m),jsr(2,m))
      else
         jss(2,m) = jsl(2,m) - jss(1,m)
      endif
      do 270 j = 1, jss(2,m)
c no more particles coming from below or back
c distribute particles coming from above or front into holes
      if (jss(1,m).gt.jsl(2,m)) then
         do 263 i = 1, idimp
         part(i,ihole(j+jsl(2,m),m),m) = rbufr(i,j,m)
  263    continue
      else
c no more holes
c distribute remaining particles from below or back into bottom
         do 267 i = 1, idimp
         part(i,j+npp(m),m) = rbufl(i,j+jss(1,m),m)
  267    continue
      endif
  270 continue
      if (jss(1,m).le.jsl(2,m)) then
         npp(m) = npp(m) + (jsl(2,m) - jss(1,m))
         jss(1,m) = jsl(2,m)
      endif
      jss(2,m) = jss(1,m) - (jsl(2,m) + jsr(2,m))
      if (jss(2,m).gt.0) then
         jss(1,m) = (jsl(2,m) + jsr(2,m))
         jsr(2,m) = jss(2,m)
      else
         jss(1,m) = jss(1,m) - jsl(2,m)
         jsr(2,m) = -jss(2,m)
      endif
      do 280 j = 1, jsr(2,m)
c holes left over
c fill up remaining holes in particle array with particles from bottom
      if (jss(2,m).gt.0) then
         j1 = npp(m) - j + 1
         j2 = jss(1,m) + jss(2,m) - j + 1
         if (j1.gt.ihole(j2,m)) then
c move particle only if it is below current hole
            do 273 i = 1, idimp
            part(i,ihole(j2,m),m) = part(i,j1,m)
  273       continue
         endif
      else
c no more holes
c distribute remaining particles from above or front into bottom
         do 277 i = 1, idimp
         part(i,j+npp(m),m) = rbufr(i,j+jss(1,m),m)
  277    continue
      endif
  280 continue
      if (jss(2,m).gt.0) then
         npp(m) = npp(m) - jsr(2,m)
      else
         npp(m) = npp(m) + jsr(2,m)
      endif
      jss(1,m) = 0
  290 continue
c check if any particles have to be passed further
      info(5+n) = max0(info(5+n),mter)
      if (ibflg(2).gt.0) then
         write (2,*) 'Info: particles being passed further = ', ibflg(2)
         if (ibflg(3).gt.0) ibflg(3) = 1
         if (iter.lt.itermax) go to 110
         ierr = -((iter-2)/2)
         write (2,*) 'Iteration overflow, iter = ', ierr
         info(1) = ierr
         go to 320
      endif
c check if buffer overflowed and more particles remain to be checked
      if (ibflg(3).gt.0) then
         nter = nter + 1
         info(3+n) = nter
         write (2,*) "new loop, nter=", nter 
         go to 20
      endif
  300 continue
c information
  320 if (nter.gt.0) then
         write (2,*) 'Info: ', nter, ' buffer overflows, nbmax=', nbmax
      endif
      return
      end
c-----------------------------------------------------------------------
      function ranorm()
c this program calculates a random number y from a gaussian distribution
c with zero mean and unit variance, according to the method of
c mueller and box:
c    y(k) = (-2*ln(x(k)))**1/2*sin(2*pi*x(k+1))
c    y(k+1) = (-2*ln(x(k)))**1/2*cos(2*pi*x(k+1)),
c where x is a random number uniformly distributed on (0,1).
c written for the ibm by viktor k. decyk, ucla
      integer r1,r2,r4,r5
      double precision ranorm,h1l,h1u,h2l,r0,r3,asc,bsc,temp
      save iflg,r1,r2,r4,r5,h1l,h1u,h2l,r0
      data r1,r2,r4,r5 /885098780,1824280461,1396483093,55318673/
      data h1l,h1u,h2l /65531.0d0,32767.0d0,65525.0d0/
      data iflg,r0 /0,0.0d0/
      if (iflg.eq.0) go to 10
      ranorm = r0
      r0 = 0.0d0
      iflg = 0
      return
   10 isc = 65536
      asc = dble(isc)
      bsc = asc*asc
      i1 = r1 - (r1/isc)*isc
      r3 = h1l*dble(r1) + asc*h1u*dble(i1)
      i1 = r3/bsc
      r3 = r3 - dble(i1)*bsc
      bsc = 0.5d0*bsc
      i1 = r2/isc
      isc = r2 - i1*isc
      r0 = h1l*dble(r2) + asc*h1u*dble(isc)
      asc = 1.0d0/bsc
      isc = r0*asc
      r2 = r0 - dble(isc)*bsc
      r3 = r3 + (dble(isc) + 2.0d0*h1u*dble(i1))
      isc = r3*asc
      r1 = r3 - dble(isc)*bsc
      temp = dsqrt(-2.0d0*dlog((dble(r1) + dble(r2)*asc)*asc))
      isc = 65536
      asc = dble(isc)
      bsc = asc*asc
      i1 = r4 - (r4/isc)*isc
      r3 = h2l*dble(r4) + asc*h1u*dble(i1)
      i1 = r3/bsc
      r3 = r3 - dble(i1)*bsc
      bsc = 0.5d0*bsc
      i1 = r5/isc
      isc = r5 - i1*isc
      r0 = h2l*dble(r5) + asc*h1u*dble(isc)
      asc = 1.0d0/bsc
      isc = r0*asc
      r5 = r0 - dble(isc)*bsc
      r3 = r3 + (dble(isc) + 2.0d0*h1u*dble(i1))
      isc = r3*asc
      r4 = r3 - dble(isc)*bsc
      r0 = 6.28318530717959d0*((dble(r4) + dble(r5)*asc)*asc)
      ranorm = temp*dsin(r0)
      r0 = temp*dcos(r0)
      iflg = 1
      return
      end
c-----------------------------------------------------------------------
