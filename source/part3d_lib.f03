! part3d_lib module for QuickPIC Open Source 1.0
! update: 04/18/2016

      module part3d_lib
         
      implicit none

!
      interface
         subroutine PRVDIST32_TWISS(part,qm,edges,npp,nps,alpha_x,alpha_y,&
         &beta_x,beta_y,emt_x,emt_y,sigz,vdx,vdy,vdz,vtz,npx,npy,npz,idimp,&
         &npmax, nx, ny,nz,x0,y0,z0,mblok,nblok,idps,ierr,gamma,lquiet)
         implicit none
         integer, intent(in) :: npmax,nblok,npx,npy,npz,idimp,nx,ny,nz,i&
         &dps,mblok
         integer, intent(inout) :: nps,npp,ierr
         real, intent(in) :: qm,sigz,x0,y0,z0,edges,vtz,vdx,vdy,vdz,gamma
         real, intent(in) :: alpha_x, alpha_y, beta_x, beta_y, emt_x, emt_y
         real, intent(inout) :: part
         logical, intent(in) :: lquiet
         dimension part(idimp,npmax,nblok)
         dimension edges(idps,nblok)
         end subroutine
      end interface
!
      interface
         subroutine PRVDIST32_RAN_PFL(part,qm,edges,npp,nps,x0,y0,z0,sigx,s&
         &igy,vtx,vty,vtz,vdx,vdy,vdz,cx,cy,npx,npy,npz,nx,ny,nz,ipbc,idimp,&
         &npmax,mblok,nblok,idps,dp,lquiet,ierr)
         implicit none
         integer, intent(in) :: npmax,nblok,npx,npy,npz,idimp,nx,ny,nz,i&
         &dps,mblok,ipbc
         integer, intent(inout) :: nps,npp,ierr
         real, intent(in) :: qm,sigx,sigy,x0,y0,z0,cx,cy,edges,vtx,vty&
         &,vtz,vdx,vdy,vdz,dp
         real, intent(inout) :: part
         logical, intent(in) :: lquiet
         dimension part(idimp,npmax,nblok)
         dimension edges(idps,nblok)
         dimension dp(nz)
         dimension cx(0:2),cy(0:2)
         end subroutine
      end interface
!
      interface
         subroutine PRVDIST32_RANDOM(part,qm,edges,npp,nps,vtx,vty,vtz,vdx,&
         &vdy,vdz,npx,npy,npz,nx,ny,nz,ipbc,idimp,npmax,mblok,nblok,idps&
         &,sigx,sigy,sigz,x0,y0,z0,cx,cy,lquiet,ierr)
         implicit none
         integer, intent(in) :: npmax,nblok,npx,npy,npz,idimp,nx,ny,nz,i&
         &dps,mblok,ipbc
         integer, intent(inout) :: nps,npp,ierr
         real, intent(in) :: qm,sigx,sigy,sigz,x0,y0,z0,cx,cy,edges,vtx,vty&
         &,vtz,vdx,vdy,vdz
         real, intent(inout) :: part
         logical, intent(in) :: lquiet
         dimension part(idimp,npmax,nblok)
         dimension edges(idps,nblok)
         dimension cx(0:2),cy(0:2)
         end subroutine
      end interface
!
      interface
         subroutine PGPOST32L(part,q,npp,noff,idimp,npmax,mnblok,nxv,&
         &nypmx,nzpmx,idds)
         implicit none
         integer, intent(in) :: npp,noff,idimp,npmax,mnblok,nxv,nypmx,nz&
         &pmx,idds
         real, intent(inout) :: part
         real, intent(inout) :: q
         dimension part(idimp,npmax,mnblok), q(nxv,nypmx,nzpmx,mnblok)
         dimension noff(idds,mnblok)
         end subroutine
      end interface
!
      interface
         subroutine PGBPUSH32L_QP(part,fxyz,bxyz,npp,noff,qbm,dt,dtc,ek,&
         &nx,ny,nz,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc,deltax,d&
         &eltaz,cofd)
         implicit none
         real, intent(inout) :: part
         integer, intent(inout) :: npp
         real, intent(inout) :: fxyz,bxyz
         real, intent(in) :: qbm,dt,dtc,ek,deltax,deltaz,cofd
         integer, intent(in) :: noff,nx,ny,nz,idimp,npmax,mnblok,nxv,nyp&
         &mx,nzpmx,idds,ipbc
         dimension part(idimp,npmax,mnblok)
         dimension fxyz(3,nxv,nypmx,nzpmx,mnblok)
         dimension bxyz(3,nxv,nypmx,nzpmx,mnblok)
         dimension noff(idds,mnblok)
         end subroutine
      end interface
!      
      interface
         subroutine PMOVE32(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole&
         &,pbuff,jsr,jsl,jss,ny,nz,kstrt,nvpy,nvpz,idimp,npmax,mblok,nbl&
         &ok,idps,nbmax,idds,ntmax,tag1,tag2,id,info)
         implicit none
         real, intent(inout) :: part, pbuff
         real, intent(in) :: edges, sbufr, sbufl, rbufr, rbufl
         integer, intent(inout) :: npp, id, info
         integer, intent(in) :: ihole, jsr, jsl, jss
         integer, intent(in) :: ny, nz, kstrt, nvpy, nvpz, idimp, npmax
         integer, intent(in) :: idps, nbmax, idds, ntmax, mblok, nblok
         integer, intent(in) :: tag1, tag2
         dimension part(idimp,npmax,mblok*nblok)
         dimension pbuff(idimp,nbmax)
         dimension edges(idps,mblok*nblok)
         dimension sbufl(idimp,nbmax,mblok*nblok)
         dimension sbufr(idimp,nbmax,mblok*nblok)
         dimension rbufl(idimp,nbmax,mblok*nblok)
         dimension rbufr(idimp,nbmax,mblok*nblok)
         dimension jsl(idds,mblok*nblok), jsr(idds,mblok*nblok)
         dimension jss(idds,mblok*nblok)
         dimension ihole(ntmax,mblok*nblok)
         dimension info(9)      
         end subroutine
      end interface
!      
      end module part3d_lib