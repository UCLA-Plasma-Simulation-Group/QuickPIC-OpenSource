      module fft2d_lib
      
      implicit none
      
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
      end module fft2d_lib