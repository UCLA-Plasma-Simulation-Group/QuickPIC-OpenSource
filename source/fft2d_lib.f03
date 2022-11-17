      module fft2d_lib

      use mpi
         
      implicit none

      public
!
      integer :: nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld     
!      
      interface
         subroutine WPFST2RINIT(mixup,sctd,indx,indy,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: indx, indy
         integer, intent(in) :: nxhyd, nxyd
         integer, dimension(nxhyd), intent(inout) :: mixup
         complex, dimension(nxyd), intent(inout) :: sctd
         end subroutine
      end interface
!
      interface
         subroutine WPFSST2R(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,&
     &indy,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: isign, ntpose, indx, indy, kstrt, nxvh, nyv
         integer, intent(in) :: kxp2, kyp, kypd, kxp2d, jblok, kblok, nxhyd, nxyd
         real, intent(inout) :: ttp
         real, intent(inout) :: f
         real, dimension(nyv,kxp2d,jblok), intent(inout) :: g
         real, dimension(kxp2+1,kyp+1,kblok), intent(inout) :: bs
         real, dimension(kxp2+1,kyp+1,jblok), intent(inout) :: br
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyd), intent(in) :: sctd
         end subroutine
      end interface
!
      interface
         subroutine WPFSCT2R(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,&
     &indy,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: isign, ntpose, indx, indy, kstrt, nxvh, nyv
         integer, intent(in) :: kxp2, kyp, kypd, kxp2d, jblok, kblok, nxhyd, nxyd
         real, intent(inout) :: ttp
         real, intent(inout) :: f
         real, dimension(nyv,kxp2d,jblok), intent(inout) :: g
         real, dimension(kxp2+1,kyp+1,kblok), intent(inout) :: bs
         real, dimension(kxp2+1,kyp+1,jblok), intent(inout) :: br
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyd), intent(in) :: sctd
         end subroutine
      end interface
!
      interface
         subroutine WPFCST2R(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,&
     &indy,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: isign, ntpose, indx, indy, kstrt, nxvh, nyv
         integer, intent(in) :: kxp2, kyp, kypd, kxp2d, jblok, kblok, nxhyd, nxyd
         real, intent(inout) :: ttp
         real, intent(inout) :: f
         real, dimension(nyv,kxp2d,jblok), intent(inout) :: g
         real, dimension(kxp2+1,kyp+1,kblok), intent(inout) :: bs
         real, dimension(kxp2+1,kyp+1,jblok), intent(inout) :: br
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyd), intent(in) :: sctd
         end subroutine
      end interface
!
      interface
         subroutine WPFCCT2R(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,&
     &indy,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: isign, ntpose, indx, indy, kstrt, nxvh, nyv
         integer, intent(in) :: kxp2, kyp, kypd, kxp2d, jblok, kblok, nxhyd, nxyd
         real, intent(inout) :: ttp
         real, intent(inout) :: f
         real, dimension(nyv,kxp2d,jblok), intent(inout) :: g
         real, dimension(kxp2+1,kyp+1,kblok), intent(inout) :: bs
         real, dimension(kxp2+1,kyp+1,jblok), intent(inout) :: br
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyd), intent(in) :: sctd
         end subroutine
      end interface
!
      interface
         subroutine WPFCST2R2(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx&
     &,indy,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: isign, ntpose, indx, indy, kstrt, nxvh, nyv
         integer, intent(in) :: kxp2, kyp, kypd, kxp2d, jblok, kblok, nxhyd, nxyd
         real, intent(inout) :: ttp
         real, intent(inout) :: f
         real, dimension(2,nyv,kxp2d,jblok), intent(inout) :: g
         real, dimension(2,kxp2+1,kyp+1,kblok), intent(inout) :: bs
         real, dimension(2,kxp2+1,kyp+1,jblok), intent(inout) :: br
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyd), intent(in) :: sctd
         end subroutine
      end interface
!
      interface
         subroutine WPFSCT2R2(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx&
     &,indy,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: isign, ntpose, indx, indy, kstrt, nxvh, nyv
         integer, intent(in) :: kxp2, kyp, kypd, kxp2d, jblok, kblok, nxhyd, nxyd
         real, intent(inout) :: ttp
         real, intent(inout) :: f
         real, dimension(2,nyv,kxp2d,jblok), intent(inout) :: g
         real, dimension(2,kxp2+1,kyp+1,kblok), intent(inout) :: bs
         real, dimension(2,kxp2+1,kyp+1,jblok), intent(inout) :: br
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyd), intent(in) :: sctd
         end subroutine
      end interface
!
      interface
         subroutine WPFCST2R3(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx&
     &,indy,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: isign, ntpose, indx, indy, kstrt, nxvh, nyv
         integer, intent(in) :: kxp2, kyp, kypd, kxp2d, jblok, kblok, nxhyd, nxyd
         real, intent(inout) :: ttp
         real, intent(inout) :: f
         real, dimension(3,nyv,kxp2d,jblok), intent(inout) :: g
         real, dimension(3,kxp2+1,kyp+1,kblok), intent(inout) :: bs
         real, dimension(3,kxp2+1,kyp+1,jblok), intent(inout) :: br
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyd), intent(in) :: sctd
         end subroutine
      end interface
!
      interface
         subroutine WPFSCT2R3(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx&
     &,indy,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: isign, ntpose, indx, indy, kstrt, nxvh, nyv
         integer, intent(in) :: kxp2, kyp, kypd, kxp2d, jblok, kblok, nxhyd, nxyd
         real, intent(inout) :: ttp
         real, intent(inout) :: f
         real, dimension(3,nyv,kxp2d,jblok), intent(inout) :: g
         real, dimension(3,kxp2+1,kyp+1,kblok), intent(inout) :: bs
         real, dimension(3,kxp2+1,kyp+1,jblok), intent(inout) :: br
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyd), intent(in) :: sctd
         end subroutine
      end interface
!
      interface
         subroutine WPFS3T2R3(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx&
     &,indy,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: isign, ntpose, indx, indy, kstrt, nxvh, nyv
         integer, intent(in) :: kxp2, kyp, kypd, kxp2d, jblok, kblok, nxhyd, nxyd
         real, intent(inout) :: ttp
         real, intent(inout) :: f
         real, dimension(3,nyv,kxp2d,jblok), intent(inout) :: g
         real, dimension(3,kxp2+1,kyp+1,kblok), intent(inout) :: bs
         real, dimension(3,kxp2+1,kyp+1,jblok), intent(inout) :: br
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyd), intent(in) :: sctd
         end subroutine
      end interface
!
      interface
         subroutine PDIVFD2(f,df,nx,ny,kstrt,ndim,nyv,kxp2,j2blok)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, ndim, nyv, kxp2, j2blok
         real, dimension(ndim,nyv,kxp2+1,j2blok), intent(inout) :: f
         real, dimension(nyv,kxp2+1,j2blok), intent(inout) :: df
         end subroutine
      end interface
!
      interface
         subroutine PGRADFD2(df,f,nx,ny,kstrt,ndim,nyv,kxp2,j2blok)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, ndim, nyv, kxp2, j2blok
         real, dimension(nyv,kxp2+1,j2blok), intent(inout) :: df
         real, dimension(ndim,nyv,kxp2+1,j2blok), intent(inout) :: f
         end subroutine
      end interface
!
      interface
         subroutine PCURLFD2(f,g,nx,ny,kstrt,nyv,kxp2,j2blok)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, nyv, kxp2, j2blok
         real, dimension(3,nyv,kxp2+1,j2blok), intent(inout) :: f, g
         end subroutine
      end interface
!
      interface
         subroutine PCURLFD22(f,g,nx,ny,kstrt,nyv,kxp2,j2blok)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, nyv, kxp2, j2blok
         real, dimension(2,nyv,kxp2+1,j2blok), intent(inout) :: f
         real, dimension(nyv,kxp2+1,j2blok), intent(inout) :: g
         end subroutine
      end interface
!
      interface
         subroutine WPPFSST2RM(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,   &
     &indx,indy,kstrt,nvp,nxvh,nyv,kxp2,kyp,kypd,kxp2d,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: isign, ntpose, indx, indy, kstrt, nvp
         integer, intent(in) :: nxvh, nyv, kxp2, kyp, kypd, kxp2d
         integer, intent(in) :: nxhyd, nxyd
         real, intent(inout) :: ttp
         real, dimension(2*nxvh,kypd), intent(inout) :: f
         real, dimension(nyv,kxp2d), intent(inout) :: g
         real, dimension(kxp2+1,kyp+1), intent(inout) :: bs, br
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyd), intent(in) :: sctd
         end subroutine
      end interface
!
      interface
         subroutine WPPFSCT2RM(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,   &
     &indx,indy,kstrt,nvp,nxvh,nyv,kxp2,kyp,kypd,kxp2d,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: isign, ntpose, indx, indy, kstrt, nvp
         integer, intent(in) :: nxvh, nyv, kxp2, kyp, kypd, kxp2d
         integer, intent(in) :: nxhyd, nxyd
         real, intent(inout) :: ttp
         real, dimension(2*nxvh,kypd), intent(inout) :: f
         real, dimension(nyv,kxp2d), intent(inout) :: g
         real, dimension(kxp2+1,kyp+1), intent(inout) :: bs, br
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyd), intent(in) :: sctd
         end subroutine
      end interface
!
      interface
         subroutine WPPFCST2RM(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,   &
     &indx,indy,kstrt,nvp,nxvh,nyv,kxp2,kyp,kypd,kxp2d,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: isign, ntpose, indx, indy, kstrt, nvp
         integer, intent(in) :: nxvh, nyv, kxp2, kyp, kypd, kxp2d
         integer, intent(in) :: nxhyd, nxyd
         real, intent(inout) :: ttp
         real, dimension(2*nxvh,kypd), intent(inout) :: f
         real, dimension(nyv,kxp2d), intent(inout) :: g
         real, dimension(kxp2+1,kyp+1), intent(inout) :: bs, br
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyd), intent(in) :: sctd
         end subroutine
      end interface
!
      interface
         subroutine WPPFCCT2RM(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,   &
     &indx,indy,kstrt,nvp,nxvh,nyv,kxp2,kyp,kypd,kxp2d,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: isign, ntpose, indx, indy, kstrt, nvp
         integer, intent(in) :: nxvh, nyv, kxp2, kyp, kypd, kxp2d
         integer, intent(in) :: nxhyd, nxyd
         real, intent(inout) :: ttp
         real, dimension(2*nxvh,kypd), intent(inout) :: f
         real, dimension(nyv,kxp2d), intent(inout) :: g
         real, dimension(kxp2+1,kyp+1), intent(inout) :: bs, br
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyd), intent(in) :: sctd
         end subroutine
      end interface
!
      interface
         subroutine PPFST2RMXX(f,isign,mixup,sctd,indx,indy,kstrt,kypi, &
     &kypp,nxvh,kypd,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: isign, indx, indy, kstrt
         integer, intent(in) :: kypi, kypp, nxvh, kypd, nxhyd, nxyd
         real, dimension(2*nxvh,kypd), intent(inout) :: f
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyd), intent(in) :: sctd
         end subroutine
      end interface
!
      interface
         subroutine PPFCT2RMXX(f,isign,mixup,sctd,indx,indy,kstrt,kypi, &
     &kypp,nxvh,kypd,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: isign, indx, indy, kstrt
         integer, intent(in) :: kypi, kypp, nxvh, kypd, nxhyd, nxyd
         real, dimension(2*nxvh,kypd), intent(inout) :: f
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyd), intent(in) :: sctd
         end subroutine
      end interface
!
      interface
         subroutine PPFST2RMXY(g,isign,mixup,sctd,indx,indy,kstrt,kxpi, &
     &kxpp,nyv,kxpd,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: isign, indx, indy, kstrt
         integer, intent(in) :: kxpi, kxpp, nyv, kxpd, nxhyd, nxyd
         real, dimension(nyv,kxpd), intent(inout) :: g
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyd), intent(in) :: sctd
         end subroutine
      end interface
!
      interface
         subroutine PPFCT2RMXY(g,isign,mixup,sctd,indx,indy,kstrt,kxpi, &
     &kxpp,nyv,kxpd,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: isign, indx, indy, kstrt
         integer, intent(in) :: kxpi, kxpp, nyv, kxpd, nxhyd, nxyd
         real, dimension(nyv,kxpd), intent(inout) :: g
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyd), intent(in) :: sctd
         end subroutine
      end interface
!
      interface
         subroutine WPPFCST2RM2(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,  &
     &indx,indy,kstrt,nvp,nxvh,nyv,kxp2,kyp,kypd,kxp2d,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: isign, ntpose, indx, indy, kstrt, nvp
         integer, intent(in) :: nxvh, nyv, kxp2, kyp, kypd, kxp2d
         integer, intent(in) :: nxhyd, nxyd
         real, intent(inout) :: ttp
         real, dimension(2,2*nxvh,kypd), intent(inout) :: f
         real, dimension(2,nyv,kxp2d), intent(inout) :: g
         real, dimension(2,kxp2+1,kyp+1), intent(inout) :: bs, br
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyd), intent(in) :: sctd
         end subroutine
      end interface
!
      interface
         subroutine WPPFSCT2RM2(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,  &
     &indx,indy,kstrt,nvp,nxvh,nyv,kxp2,kyp,kypd,kxp2d,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: isign, ntpose, indx, indy, kstrt, nvp
         integer, intent(in) :: nxvh, nyv, kxp2, kyp, kypd, kxp2d
         integer, intent(in) :: nxhyd, nxyd
         real, intent(inout) :: ttp
         real, dimension(2,2*nxvh,kypd), intent(inout) :: f
         real, dimension(2,nyv,kxp2d), intent(inout) :: g
         real, dimension(2,kxp2+1,kyp+1), intent(inout) :: bs, br
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyd), intent(in) :: sctd
         end subroutine
      end interface
!
      interface
         subroutine PPFCST2RM2X(f,isign,mixup,sctd,indx,indy,kstrt,kypi,&
     &kypp,nxvh,kypd,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: isign, indx, indy, kstrt
         integer, intent(in) :: kypi, kypp, nxvh, kypd, nxhyd, nxyd
         real, dimension(2,2*nxvh,kypd), intent(inout) :: f
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyd), intent(in) :: sctd
         end subroutine
      end interface
!
      interface
         subroutine PPFSCT2RM2X(f,isign,mixup,sctd,indx,indy,kstrt,kypi,&
     &kypp,nxvh,kypd,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: isign, indx, indy, kstrt
         integer, intent(in) :: kypi, kypp, nxvh, kypd, nxhyd, nxyd
         real, dimension(2,2*nxvh,kypd), intent(inout) :: f
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyd), intent(in) :: sctd
         end subroutine
      end interface
!
      interface
         subroutine PPFSCT2RM2Y(g,isign,mixup,sctd,indx,indy,kstrt,kxpi,&
     &kxpp,nyv,kxpd,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: isign, indx, indy, kstrt
         integer, intent(in) :: kxpi, kxpp, nyv, kxpd, nxhyd, nxyd
         real, dimension(2,nyv,kxpd), intent(inout) :: g
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyd), intent(in) :: sctd
         end subroutine
      end interface
!
      interface
         subroutine PPFCST2RM2Y(g,isign,mixup,sctd,indx,indy,kstrt,kxpi,&
     &kxpp,nyv,kxpd,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: isign, indx, indy, kstrt
         integer, intent(in) :: kxpi, kxpp, nyv, kxpd, nxhyd, nxyd
         real, dimension(2,nyv,kxpd), intent(inout) :: g
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyd), intent(in) :: sctd
         end subroutine
      end interface
!
      interface
         subroutine WPPFCST2RM3(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,  &
     &indx,indy,kstrt,nvp,nxvh,nyv,kxp2,kyp,kypd,kxp2d,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: isign, ntpose, indx, indy, kstrt, nvp
         integer, intent(in) :: nxvh, nyv, kxp2, kyp, kypd, kxp2d
         integer, intent(in) :: nxhyd, nxyd
         real, intent(inout) :: ttp
         real, dimension(3,2*nxvh,kypd), intent(inout) :: f
         real, dimension(3,nyv,kxp2d), intent(inout) :: g
         real, dimension(3,kxp2+1,kyp+1), intent(inout) :: bs, br
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyd), intent(in) :: sctd
         end subroutine
      end interface
!
      interface
         subroutine WPPFSCT2RM3(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,  &
     &indx,indy,kstrt,nvp,nxvh,nyv,kxp2,kyp,kypd,kxp2d,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: isign, ntpose, indx, indy, kstrt, nvp
         integer, intent(in) :: nxvh, nyv, kxp2, kyp, kypd, kxp2d
         integer, intent(in) :: nxhyd, nxyd
         real, intent(inout) :: ttp
         real, dimension(3,2*nxvh,kypd), intent(inout) :: f
         real, dimension(3,nyv,kxp2d), intent(inout) :: g
         real, dimension(3,kxp2+1,kyp+1), intent(inout) :: bs, br
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyd), intent(in) :: sctd
         end subroutine
      end interface
!
      interface
         subroutine PPFCSST2RM3X(f,isign,mixup,sctd,indx,indy,kstrt,kypi&
     &,kypp,nxvh,kypd,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: isign, indx, indy, kstrt
         integer, intent(in) :: kypi, kypp, nxvh, kypd, nxhyd, nxyd
         real, dimension(3,2*nxvh,kypd), intent(inout) :: f
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyd), intent(in) :: sctd
         end subroutine
      end interface
!
      interface
         subroutine PPFSCCT2RM3X(f,isign,mixup,sctd,indx,indy,kstrt,kypi&
     &,kypp,nxvh,kypd,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: isign, indx, indy, kstrt
         integer, intent(in) :: kypi, kypp, nxvh, kypd, nxhyd, nxyd
         real, dimension(3,2*nxvh,kypd), intent(inout) :: f
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyd), intent(in) :: sctd
         end subroutine
      end interface
!
      interface
         subroutine PPFSCST2RM3Y(g,isign,mixup,sctd,indx,indy,kstrt,kxpi&
     &,kxpp,nyv,kxpd,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: isign, indx, indy, kstrt
         integer, intent(in) :: kxpi, kxpp, nyv, kxpd, nxhyd, nxyd
         real, dimension(3,nyv,kxpd), intent(inout) :: g
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyd), intent(in) :: sctd
         end subroutine
      end interface
!
      interface
         subroutine PPFCSCT2RM3Y(g,isign,mixup,sctd,indx,indy,kstrt,kxpi&
     &,kxpp,nyv,kxpd,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: isign, indx, indy, kstrt
         integer, intent(in) :: kxpi, kxpp, nyv, kxpd, nxhyd, nxyd
         real, dimension(3,nyv,kxpd), intent(inout) :: g
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyd), intent(in) :: sctd
         end subroutine
      end interface
!
      interface
         subroutine WPPFSCT2RM4(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,  &
     &indx,indy,kstrt,nvp,nxvh,nyv,kxp2,kyp,kypd,kxp2d,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: isign, ntpose, indx, indy, kstrt, nvp
         integer, intent(in) :: nxvh, nyv, kxp2, kyp, kypd, kxp2d
         integer, intent(in) :: nxhyd, nxyd
         real, intent(inout) :: ttp
         real, dimension(4,2*nxvh,kypd), intent(inout) :: f
         real, dimension(4,nyv,kxp2d), intent(inout) :: g
         real, dimension(4,kxp2+1,kyp+1), intent(inout) :: bs, br
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyd), intent(in) :: sctd
         end subroutine
      end interface
!
      interface
         subroutine PPFSCCST2RM4X(f,isign,mixup,sctd,indx,indy,kstrt,   &
     &kypi,kypp,nxvh,kypd,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: isign, indx, indy, kstrt
         integer, intent(in) :: kypi, kypp, nxvh, kypd, nxhyd, nxyd
         real, dimension(4,2*nxvh,kypd), intent(inout) :: f
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyd), intent(in) :: sctd
         end subroutine
      end interface
!
      interface
         subroutine PPFSCSCT2RM4Y(g,isign,mixup,sctd,indx,indy,kstrt,   &
     &kxpi,kxpp,nyv,kxpd,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: isign, indx, indy, kstrt
         integer, intent(in) :: kxpi, kxpp, nyv, kxpd, nxhyd, nxyd
         real, dimension(4,nyv,kxpd), intent(inout) :: g
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyd), intent(in) :: sctd
         end subroutine
      end interface
!
      interface
         subroutine WPPFSCT2RM22(f,g,bs,br,isign,ntpose,mixup,sctd,ttp, &
     &indx,indy,kstrt,nvp,nxvh,nyv,kxp2,kyp,kypd,kxp2d,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: isign, ntpose, indx, indy, kstrt, nvp
         integer, intent(in) :: nxvh, nyv, kxp2, kyp, kypd, kxp2d
         integer, intent(in) :: nxhyd, nxyd
         real, intent(inout) :: ttp
         real, dimension(2,2*nxvh,kypd), intent(inout) :: f
         real, dimension(2,nyv,kxp2d), intent(inout) :: g
         real, dimension(2,kxp2+1,kyp+1), intent(inout) :: bs, br
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyd), intent(in) :: sctd
         end subroutine
      end interface
!
      interface
         subroutine PPFSCCST2RM22X(f,isign,mixup,sctd,indx,indy,kstrt,  &
     &kypi,kypp,nxvh,kypd,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: isign, indx, indy, kstrt
         integer, intent(in) :: kypi, kypp, nxvh, kypd, nxhyd, nxyd
         real, dimension(2,2*nxvh,kypd), intent(inout) :: f
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyd), intent(in) :: sctd
         end subroutine
      end interface
!
      interface
         subroutine PPFSCSCT2RM22Y(g,isign,mixup,sctd,indx,indy,kstrt,  &
     &kxpi,kxpp,nyv,kxpd,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: isign, indx, indy, kstrt
         integer, intent(in) :: kxpi, kxpp, nyv, kxpd, nxhyd, nxyd
         real, dimension(2,nyv,kxpd), intent(inout) :: g
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyd), intent(in) :: sctd
         end subroutine
      end interface
!
      interface
         subroutine WPPFSST2RM23(f,g,bs,br,isign,ntpose,mixup,sctd,ttp, &
     &indx,indy,kstrt,nvp,nxvh,nyv,kxp2,kyp,kypd,kxp2d,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: isign, ntpose, indx, indy, kstrt, nvp
         integer, intent(in) :: nxvh, nyv, kxp2, kyp, kypd, kxp2d
         integer, intent(in) :: nxhyd, nxyd
         real, intent(inout) :: ttp
         real, dimension(3,2*nxvh,kypd), intent(inout) :: f
         real, dimension(3,nyv,kxp2d), intent(inout) :: g
         real, dimension(3,kxp2+1,kyp+1), intent(inout) :: bs, br
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyd), intent(in) :: sctd
         end subroutine
      end interface
!
      interface
         subroutine PPFSSCT2RM32X(f,isign,mixup,sctd,indx,indy,kstrt,   &
     &kypi,kypp,nxvh,kypd,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: isign, indx, indy, kstrt
         integer, intent(in) :: kypi, kypp, nxvh, kypd, nxhyd, nxyd
         real, dimension(3,2*nxvh,kypd), intent(inout) :: f
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyd), intent(in) :: sctd
         end subroutine
      end interface
!
      interface
         subroutine PPFSSCT2RM23Y(g,isign,mixup,sctd,indx,indy,kstrt,   &
     &kxpi,kxpp,nyv,kxpd,nxhyd,nxyd)
         implicit none
         integer, intent(in) :: isign, indx, indy, kstrt
         integer, intent(in) :: kxpi, kxpp, nyv, kxpd, nxhyd, nxyd
         real, dimension(3,nyv,kxpd), intent(inout) :: g
         integer, dimension(nxhyd), intent(in) :: mixup
         complex, dimension(nxyd), intent(in) :: sctd
         end subroutine
      end interface
!
      interface
         subroutine MPPDIVFD2(f,df,nx,ny,kstrt,ndim,nyv,kxp2)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, ndim, nyv, kxp2
         real, dimension(3,nyv,kxp2+1), intent(in) :: f
         real, dimension(nyv,kxp2+1), intent(inout) :: df
         end subroutine
      end interface
!
      interface
         subroutine MPPGRADFD2(df,f,nx,ny,kstrt,ndim,nyv,kxp2)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, ndim, nyv, kxp2
         real, dimension(nyv,kxp2+1), intent(in) :: df
         real, dimension(3,nyv,kxp2+1), intent(inout) :: f
         end subroutine
      end interface
!
      interface
         subroutine MPPCURLFD2(f,g,nx,ny,kstrt,nyv,kxp2)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, nyv, kxp2
         real, dimension(3,nyv,kxp2+1), intent(in) :: f
         real, dimension(3,nyv,kxp2+1), intent(inout) :: g
         end subroutine
      end interface
!
      contains
!-----------------------------------------------------------------------
      subroutine PPRTPOSE(f,g,s,t,nx,ny,kxp,kyp,kstrt,nvp,nxv,nyv,kxpd, &
     &kypd)
! this subroutine performs a transpose of a real matrix f, distributed
! in y, to a real matrix g, distributed in x, that is,
! g(k+kyp*(m-1),j,l) = f(j+kxp*(l-1),k,m), where
! 1 <= j <= kxp, 1 <= k <= kyp, 1 <= l <= nx/kxp, 1 <= m <= ny/kyp
! and where indices l and m can be distributed across processors.
! includes an extra guard cell for last row and column
! this subroutine sends and receives one message at a time, either
! synchronously or asynchronously. it uses a minimum of system resources
! f = real input array
! g = real output array
! s, t = real scratch arrays
! nx/ny = number of points in x/y
! kxp/kyp = number of data values per block in x/y
! kstrt = starting data block number
! nvp = number of real or virtual processors
! nxv = first dimension of f, nxv >= nx+1
! nyv = first dimension of g, nyv >= ny+1
! kypd = second dimension of f, kypd >= kyp+1
! kxpd = second dimension of g, kxpd >= kxp+1
      implicit none
      integer, intent(in) :: nx, ny, kxp, kyp, kstrt, nvp, nxv, nyv
      integer, intent(in) :: kxpd, kypd
      real, dimension(nxv,kypd), intent(in) :: f
      real, dimension(nyv,kxpd), intent(inout) :: g
      real, dimension((kxp+1)*(kyp+1)), intent(inout) :: s, t
! lgrp = current communicator
! mreal = default datatype for reals
! local data
      integer :: n, j, k, ks, kxps, kyps, kxyp, id, joff, koff, ld
      integer :: nx1, ny1, kxb, kyb
      integer :: ierr, msid
      integer, dimension(10) :: istatus
      nx1 = nx + 1
      ny1 = ny + 1
! ks = processor id
      ks = kstrt - 1
! kxps = actual size used in x direction
      kxps = min(kxp,max(0,nx-kxp*ks))
! kyps = actual size used in y direction
      kyps = min(kyp,max(0,ny-kyp*ks))
! kxb = minimum number of processors needed in x direction
      kxb = (nx - 1)/kxp + 1
! kyb = minimum number of processors needed in y direction
      kyb = (ny - 1)/kyp + 1
! add extra word for last processor in x
      if (ks==(kxb-1)) kxps = kxps + 1
! add extra word for last processor in y
      if (ks==(kyb-1)) kyps = kyps + 1
! kxyp = maximum amount of data to be received
      kxyp = (kxp + 1)*(kyp + 1)
! special case for one processor
      if (nvp==1) then
!$OMP PARALLEL DO PRIVATE(j,k)
         do k = 1, ny1
            do j = 1, nx1
               g(k,j) = f(j,k)
            enddo
         enddo
!$OMP END PARALLEL DO
         return
      endif
! this segment is used for shared memory computers
!     do m = 1, min(ny,nvp)
!        koff = kyp*(m - 1)
!        kyps = min(kyp,max(0,ny-koff))
!        if (m==kyb) kyps = kyps + 1
!        do k = 1, kyps
!           do l = 1, min(nx,nvp)
!              joff = kxp*(l - 1)
!              kxps = min(kxp,max(0,nx-joff))
!              if (l==kxb) kxps = kxps + 1
!              do j = 1, kxps
!                 g(k+koff,j+joff) = f(j+joff,k+koff)
!              enddo
!           enddo
!        enddo
!     enddo
! this segment is used for mpi computers
      do n = 1, nvp
         id = n - ks - 1
         if (id < 0) id = id + nvp
! extract data to send
         joff = kxp*id
         ld = min(kxp,max(0,nx-joff))
! add extra word for last processor in x
         if (id==(kxb-1)) ld = ld + 1
!$OMP PARALLEL DO PRIVATE(j,k)
         do k = 1, kyps
            do j = 1, ld
               s(j+ld*(k-1)) = f(j+joff,k)
            enddo
         enddo
!$OMP END PARALLEL DO
         ld = ld*kyps
! post receive
         call MPI_IRECV(t,kxyp,mreal,id,n,lgrp,msid,ierr)
! send data
         call MPI_SEND(s,ld,mreal,id,n,lgrp,ierr)
! receive data
         call MPI_WAIT(msid,istatus,ierr)
! insert data received
         koff = kyp*id
         ld = min(kyp,max(0,ny-koff))
! add extra word for last processor in y
         if (id==(kyb-1)) ld = ld + 1
!$OMP PARALLEL DO PRIVATE(j,k)
         do k = 1, ld
            do j = 1, kxps
               g(k+koff,j) = t(j+kxps*(k-1))
            enddo
         enddo
!$OMP END PARALLEL DO
      enddo
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPRNTPOSE(f,g,s,t,nx,ny,kxp,kyp,kstrt,nvp,ndim,nxv,nyv,&
     &kxpd,kypd)
! this subroutine performs a transpose of a real matrix f, distributed
! in y, to a real matrix g, distributed in x, that is,
! g(1:ndim,k+kyp*(m-1),j,l) = f(1:ndim,j+kxp*(l-1),k,m), where
! 1 <= j <= kxp, 1 <= k <= kyp, 1 <= l <= nx/kxp, 1 <= m <= ny/kyp
! and where indices l and m can be distributed across processors.
! includes an extra guard cell for last row and column
! this subroutine sends and receives one message at a time, either
! synchronously or asynchronously. it uses a minimum of system resources
! f = real input array
! g = real output array
! s, t = real scratch arrays
! nx/ny = number of points in x/y
! kxp/kyp = number of data values per block in x/y
! kstrt = starting data block number
! nvp = number of real or virtual processors
! ndim = leading dimension of arrays f and g
! nxv = second dimension of f, nxv >= nx+1
! nyv = second dimension of g, nyv >= ny+1
! kypd = third dimension of f, kypd >= kyp+1
! kxpd = third dimension of g, kxpd >= kxp+1
      implicit none
      integer, intent(in) :: nx, ny, kxp, kyp, kstrt, nvp, ndim
      integer, intent(in) :: nxv, nyv, kxpd, kypd
      real, dimension(ndim,nxv,kypd), intent(in) :: f
      real, dimension(ndim,nyv,kxpd), intent(inout) :: g
      real, dimension(ndim,(kxp+1)*(kyp+1)), intent(inout) :: s, t
! lgrp = current communicator
! mreal = default datatype for reals
! local data
      integer :: i, n, j, k, ks, kxps, kyps, kxyp, id, joff, koff, ld
      integer :: nx1, ny1, kxb, kyb
      integer :: ierr, msid
      integer, dimension(10) :: istatus
      nx1 = nx + 1
      ny1 = ny + 1
! ks = processor id
      ks = kstrt - 1
! kxps = actual size used in x direction
      kxps = min(kxp,max(0,nx-kxp*ks))
! kyps = actual size used in y direction
      kyps = min(kyp,max(0,ny-kyp*ks))
! kxb = minimum number of processors needed in x direction
      kxb = (nx - 1)/kxp + 1
! kyb = minimum number of processors needed in y direction
      kyb = (ny - 1)/kyp + 1
! add extra word for last processor in x
      if (ks==(kxb-1)) kxps = kxps + 1
! add extra word for last processor in y
      if (ks==(kyb-1)) kyps = kyps + 1
! kxyp = maximum amount of data to be received
      kxyp = ndim*(kxp + 1)*(kyp + 1)
! special case for one processor
      if (nvp==1) then
!$OMP PARALLEL DO PRIVATE(i,j,k)
         do k = 1, ny1
            do j = 1, nx1
               do i = 1, ndim
                  g(i,k,j) = f(i,j,k)
               enddo
            enddo
         enddo
!$OMP END PARALLEL DO
         return
      endif
! this segment is used for shared memory computers
!     do m = 1, min(ny,nvp)
!        koff = kyp*(m - 1)
!        kyps = min(kyp,max(0,ny-koff))
!        if (m==kyb) kyps = kyps + 1
!        do k = 1, kyps
!           do l = 1, min(nx,nvp)
!              joff = kxp*(l - 1)
!              kxps = min(kxp,max(0,nx-joff))
!              if (l==kxb) kxps = kxps + 1
!              do j = 1, kxps
!                 do i = 1, ndim
!                    g(i,k+koff,j+joff) = f(i,j+joff,k+koff)
!                 enddo
!              enddo
!           enddo
!        enddo
!     enddo
! this segment is used for mpi computers
      do n = 1, nvp
         id = n - ks - 1
         if (id.lt.0) id = id + nvp
! extract data to send
         joff = kxp*id
         ld = min(kxp,max(0,nx-joff))
! add extra word for last processor in x
         if (id==(kxb-1)) ld = ld + 1
!$OMP PARALLEL DO PRIVATE(i,j,k)
         do k = 1, kyps
            do j = 1, ld
               do i = 1, ndim
                  s(i,j+ld*(k-1)) = f(i,j+joff,k)
               enddo
            enddo
         enddo
!$OMP END PARALLEL DO
         ld = ndim*ld*kyps
! post receive
         call MPI_IRECV(t,kxyp,mreal,id,n,lgrp,msid,ierr)
! send data
         call MPI_SEND(s,ld,mreal,id,n,lgrp,ierr)
! receive data
         call MPI_WAIT(msid,istatus,ierr)
! insert data received
         koff = kyp*id
         ld = min(kyp,max(0,ny-koff))
! add extra word for last processor in y
         if (id==(kyb-1)) ld = ld + 1
!$OMP PARALLEL DO PRIVATE(i,j,k)
         do k = 1, ld
            do j = 1, kxps
               do i = 1, ndim
                  g(i,k+koff,j) = t(i,j+kxps*(k-1))
               enddo
            enddo
         enddo
!$OMP END PARALLEL DO
      enddo
      end subroutine
!
      end module fft2d_lib
!-----------------------------------------------------------------------
      subroutine PPRTPOSE(f,g,s,t,nx,ny,kxp,kyp,kstrt,nvp,nxv,nyv,kxpd, &
     &kypd)
      use fft2d_lib, only: SUB => PPRTPOSE
      implicit none
      integer, intent(in) :: nx, ny, kxp, kyp, kstrt, nvp, nxv, nyv
      integer, intent(in) :: kxpd, kypd
      real, dimension(nxv,kypd), intent(in) :: f
      real, dimension(nyv,kxpd), intent(inout) :: g
      real, dimension((kxp+1)*(kyp+1)), intent(inout) :: s, t
      call SUB(f,g,s,t,nx,ny,kxp,kyp,kstrt,nvp,nxv,nyv,kxpd,kypd)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine PPRNTPOSE(f,g,s,t,nx,ny,kxp,kyp,kstrt,nvp,ndim,nxv,nyv,&
     &kxpd,kypd)
      use fft2d_lib, only: SUB => PPRNTPOSE
      implicit none
      integer, intent(in) :: nx, ny, kxp, kyp, kstrt, nvp, ndim
      integer, intent(in) :: nxv, nyv, kxpd, kypd
      real, dimension(ndim,nxv,kypd), intent(in) :: f
      real, dimension(ndim,nyv,kxpd), intent(inout) :: g
      real, dimension(ndim,kxp*kyp), intent(inout) :: s, t
      call SUB(f,g,s,t,nx,ny,kxp,kyp,kstrt,nvp,ndim,nxv,nyv,kxpd,kypd)
      end subroutine

      