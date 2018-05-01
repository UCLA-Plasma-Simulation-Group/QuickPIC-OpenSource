! ufield3d_lib module for QuickPIC Open Source 1.0
! update: 04/18/2016

      module ufield3d_lib
         
      implicit none

!
      interface
         subroutine PCGUARD32L(f,scs,kstrt,nvpy,nvpz,nxv,nypmx,nzpmx,mbl&
         &ok,nblok,kyp,kzp,ngds,tag1,tag2,rid,sid,ierr)
         implicit none
         integer, intent(in) :: kstrt, nvpy, nvpz, nxv, nypmx, nzpmx, mb&
         &lok, nblok, ngds, kyp, kzp, tag1, tag2
         real, intent(inout) :: f, scs
         integer, intent(inout) :: rid, sid, ierr
         dimension f(nxv,nypmx,nzpmx)
         dimension scs(nxv,nzpmx,2*ngds)
         end subroutine
      end interface
!
      interface
         subroutine PACGUARD32L(f,scs,kstrt,nvpy,nvpz,nx,nxv,nypmx,nzpmx&
         &,mblok,nblok,kyp,kzp,ngds)
         implicit none
         integer, intent(in) :: kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx&
         &, mblok, nblok, ngds, kyp, kzp
         real, intent(inout) :: f, scs
         dimension f(3,nxv,nypmx,nzpmx,mblok*nblok)
         dimension scs(3,nxv,nzpmx,2*ngds,mblok*nblok)
         end subroutine
      end interface     
!
      interface
         subroutine PAGUARD32L(f,scs,scr,kstrt,nvpy,nvpz,nx,nxv,nypmx,nz&
         &pmx,mblok,nblok,kyp,kzp,ngds,tag1,tag2,id,ierr)
         implicit none
         integer, intent(in) :: kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx&
         &, mblok, nblok, ngds, kyp, kzp, tag1, tag2
         integer, intent(inout) :: id, ierr
         real, intent(inout) :: f, scs, scr
         dimension f(nxv,nypmx,nzpmx,mblok*nblok)
         dimension scs(nxv,nzpmx,2*ngds,mblok*nblok)
         dimension scr(nxv,nypmx,ngds,mblok*nblok)
         end subroutine
      end interface     
!
      end module ufield3d_lib