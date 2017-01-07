      module fpois2d_lib
!
      interface
         subroutine PPOISDX2(q,fx,fy,isign,ffd,ax,ay,affp,we,nx,ny,kstrt&
     &,ny2d,kxp2,j2blok,nyd)
         implicit none
         real, intent(in) :: ax, ay, affp, we
         integer, intent(in) :: isign, nx, ny, kstrt, ny2d, kxp2, j2blok, nyd
         complex, dimension(ny2d,kxp2,j2blok), intent(inout) :: q, fx, fy
         complex, dimension(nyd,kxp2,j2blok), intent(inout) :: ffd
         end subroutine
      end interface
!
      interface
         subroutine PPOISD2(q,fx,fy,isign,ffd,ax,ay,affp,we,nx,ny,kstrt,&
     &nyv,kxp2,j2blok,nyd)
         implicit none
         real, intent(in) :: ax, ay, affp, we
         integer, intent(in) :: isign, nx, ny, kstrt, nyv, kxp2, j2blok, nyd
         real, dimension(nyv,kxp2+1,j2blok), intent(inout) :: q, fx, fy
         complex, dimension(nyd,kxp2,j2blok), intent(in) :: ffd
         end subroutine
      end interface
!
      interface
         subroutine PPOISD22(q,fxy,isign,ffd,ax,ay,affp,we,nx,ny,kstrt,n&
     &yv,kxp2,j2blok,nyd)
         implicit none
         real, intent(in) :: ax, ay, affp, we
         integer, intent(in) :: isign, nx, ny, kstrt, nyv, kxp2, j2blok, nyd
         real, dimension(nyv,kxp2+1,j2blok), intent(inout) :: q
         real, dimension(2,nyv,kxp2+1,j2blok), intent(inout) :: fxy
         complex, dimension(nyd,kxp2,j2blok), intent(in) :: ffd
         end subroutine
      end interface
!
      interface
         subroutine PPOISD23(q,fxy,isign,ffd,ax,ay,affp,we,nx,ny,kstrt,n&
     &yv,kxp2,j2blok,nyd)
         implicit none
         real, intent(in) :: ax, ay, affp, we
         integer, intent(in) :: isign, nx, ny, kstrt, nyv, kxp2, j2blok, nyd
         real, dimension(nyv,kxp2+1,j2blok), intent(inout) :: q
         real, dimension(3,nyv,kxp2+1,j2blok), intent(inout) :: fxy
         complex, dimension(nyd,kxp2,j2blok), intent(in) :: ffd
         end subroutine
      end interface
!
      interface
         subroutine PBPOISD22(cu,bxy,bz,isign,ffd,ax,ay,affp,ci,wm,nx,ny&
     &,kstrt,nyv,kxp2,j2blok,nyd)
         implicit none
         real, intent(in) :: ax, ay, affp, ci, wm
         integer, intent(in) :: isign, nx, ny, kstrt, nyv, kxp2, j2blok, nyd
         real, dimension(2,nyv,kxp2+1,j2blok), intent(inout) :: cu, bxy
         real, dimension(nyv,kxp2+1,j2blok), intent(inout) :: bz
         complex, dimension(nyd,kxp2,j2blok), intent(in) :: ffd
         end subroutine
      end interface
!      
      interface
         subroutine PBPOISD23(cu,bxy,isign,ffd,ax,ay,affp,ci,wm,nx,ny,ks&
     &trt,nyv,kxp2,j2blok,nyd)
         implicit none
         real, intent(in) :: ax, ay, affp, ci, wm
         integer, intent(in) :: isign, nx, ny, kstrt, nyv, kxp2, j2blok, nyd
         real, dimension(3,nyv,kxp2+1,j2blok), intent(inout) :: cu, bxy
         complex, dimension(nyd,kxp2,j2blok), intent(in) :: ffd
         end subroutine
      end interface
!
      interface
         subroutine PBPOISD22N_QP(cu,dcu,amu,bxy,bz,isign,ffd,ax,ay,affp&
     &,ci,wm,nx,ny,kstrt,nyv,kxp2,j2blok,nyd,aa,dex)
		 implicit none
         real, intent(in) :: ax, ay, affp, ci, wm, aa,dex
         integer, intent(in) :: isign, nx, ny, kstrt, nyv, kxp2, j2blok, nyd
         real, dimension(2,nyv,kxp2+1,j2blok), intent(inout) :: dcu, bxy
         real, dimension(3,nyv,kxp2+1,j2blok), intent(inout) :: cu, amu
         real, dimension(nyv,kxp2+1,j2blok), intent(inout) :: bz
         complex, dimension(nyd,kxp2,j2blok), intent(in) :: ffd
		 end subroutine
      end interface
!
      end module fpois2d_lib
