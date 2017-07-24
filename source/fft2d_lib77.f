c 2d parallel PIC library for fast fourier transforms
c-----------------------------------------------------------------------
      subroutine WPFST2RINIT(mixup,sctd,indx,indy,nxhyd,nxyd)
c this subroutine calculates tables needed by a two dimensional
c fast real sine and cosine transforms and their inverses.
c input: indx, indy, nxhyd, nxyd
c output: mixup, sctd
c mixup = array of bit reversed addresses
c sctd = sine/cosine table
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c nxhyd = maximum of (nx/2,ny)
c nxyd = maximum of (nx,ny)
c written by viktor k. decyk, ucla
      implicit none
      integer indx, indy, nxhyd, nxyd
      integer mixup
      complex sctd
      dimension mixup(nxhyd), sctd(nxyd)
c local data
      integer indx1, indx1y, nx, ny, nxy, nxhy
      integer j, k, lb, ll, jb, it
      real dnxy, arg
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
c bit-reverse index table: mixup(j) = 1 + reversed bits of (j - 1)
      do 20 j = 1, nxhy
      lb = j - 1
      ll = 0
      do 10 k = 1, indx1y
      jb = lb/2
      it = lb - 2*jb
      lb = jb
      ll = 2*ll + it
   10 continue
      mixup(j) = ll + 1
   20 continue
c sine/cosine table for the angles n*pi/nxy
      dnxy = 0.5*6.28318530717959/float(nxy)
      do 30 j = 1, nxy
      arg = dnxy*float(j - 1)
      sctd(j) = cmplx(cos(arg),-sin(arg))
   30 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine WPFSST2R(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,ind
     1y,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhyd,nxyd)
c wrapper function for parallel real sine/sine transform
      implicit none
      integer isign, ntpose, mixup, indx, indy, kstrt, nxvh, nyv
      integer kxp2, kyp, kypd, kxp2d, jblok, kblok, nxhyd, nxyd
      real f, g, bs, br, ttp
      complex sctd
      dimension f(2*nxvh,kypd,kblok), g(nyv,kxp2d,jblok)
      dimension bs(kxp2+1,kyp+1,kblok), br(kxp2+1,kyp+1,jblok)
      dimension mixup(nxhyd), sctd(nxyd)
c local data
      integer nx, ny, kxpi, kypi
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
c inverse fourier transform
      if (isign.lt.0) then
c perform x sine transform
         call PFST2RXX(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kyp,n
     1xvh,kypd,kblok,nxhyd,nxyd)
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PRTPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,nyv,kxp2,kyp,kxp2d,ky
     1pd,jblok,kblok)
         call PWTIMERA(1,ttp,dtime)
c perform y sine transform
         call PFST2RXY(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxp2
     1,nyv,kxp2d,jblok,nxhyd,nxyd)
c transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PRTPOSE(g,f,br,bs,ny,nx,kstrt,nyv,2*nxvh,kyp,kxp2,kypd,
     1kxp2d,kblok,jblok)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PRTPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,nyv,kxp2,kyp,kxp2d
     1,kypd,jblok,kblok)
            call PWTIMERA(1,tf,dtime)
         endif
c perform y sine transform
         call PFST2RXY(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxp2
     1,nyv,kxp2d,jblok,nxhyd,nxyd)
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PRTPOSE(g,f,br,bs,ny,nx,kstrt,nyv,2*nxvh,kyp,kxp2,kypd,kxp
     12d,kblok,jblok)
         call PWTIMERA(1,ttp,dtime)
c perform x sine transform
         call PFST2RXX(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kyp,n
     1xvh,kypd,kblok,nxhyd,nxyd)
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine WPFCCT2R(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,ind
     1y,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhyd,nxyd)
c wrapper function for parallel real cosine/cosine transform
      implicit none
      integer isign, ntpose, mixup, indx, indy, kstrt, nxvh, nyv
      integer kxp2, kyp, kypd, kxp2d, jblok, kblok, nxhyd, nxyd
      real f, g, bs, br, ttp
      complex sctd
      dimension f(2*nxvh,kypd,kblok), g(nyv,kxp2d,jblok)
      dimension bs(kxp2+1,kyp+1,kblok), br(kxp2+1,kyp+1,jblok)
      dimension mixup(nxhyd), sctd(nxyd)
c local data
      integer nx, ny, kxpi, kypi
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
c inverse fourier transform
      if (isign.lt.0) then
c perform x cosine transform
         call PFCT2RXX(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kyp,n
     1xvh,kypd,kblok,nxhyd,nxyd)
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PRTPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,nyv,kxp2,kyp,kxp2d,ky
     1pd,jblok,kblok)
         call PWTIMERA(1,ttp,dtime)
c perform y cosine transform
         call PFCT2RXY(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxp2
     1,nyv,kxp2d,jblok,nxhyd,nxyd)
c transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PRTPOSE(g,f,br,bs,ny,nx,kstrt,nyv,2*nxvh,kyp,kxp2,kypd,
     1kxp2d,kblok,jblok)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PRTPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,nyv,kxp2,kyp,kxp2d
     1,kypd,jblok,kblok)
            call PWTIMERA(1,tf,dtime)
         endif
c perform y cosine transform
         call PFCT2RXY(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxp2
     1,nyv,kxp2d,jblok,nxhyd,nxyd)
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PRTPOSE(g,f,br,bs,ny,nx,kstrt,nyv,2*nxvh,kyp,kxp2,kypd,kxp
     12d,kblok,jblok)
         call PWTIMERA(1,ttp,dtime)
c perform x cosine transform
         call PFCT2RXX(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kyp,n
     1xvh,kypd,kblok,nxhyd,nxyd)
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine PFST2RXX(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,ky
     1pp,nxvh,kypd,kblok,nxhyd,nxyd)
c this subroutine performs the x part of a two dimensional fast real
c sine transform and its inverse, for a subset of y,
c using real arithmetic, for data which is distributed in blocks
c algorithm is described in Numerical Recipies in Fortran, Second Ed.,
c by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
c [Cambridge Univ. Press, 1992], p. 508.
c for isign = (-1,1), input: all, output: f
c approximate flop count: N*(5*log2(N) + 18)/nvp
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, an inverse sine transform is performed
c f(n,k,i) = (1/nx*ny)*sum(f(j,k,i)*sin(pi*n*j/nx))
c if isign = 1, a forward sine transform is performed
c f(j,k,i) = sum(f(n,k,i)*sin(pi*n*j/nx))
c mixup = array of bit reversed addresses
c sctd = sine/cosine table
c kstrt = starting data block number
c kyp = number of data values per block in y
c kypi = initial y index used
c kypp = number of y indices used
c nxvh = first dimension of f >= nx/2 + 1
c kypd = second dimension of f >= kyp + 1
c kblok = number of data blocks in y
c nxhyd = maximum of (nx/2,ny)
c nxyd = maximum of (nx,ny)
c written by viktor k. decyk, ucla
      implicit none
      integer isign, mixup, indx, indy, kstrt, kyp, kypi, kypp
      integer nxvh, kypd, kblok, nxhyd, nxyd
      real f
      complex sctd
      dimension f(2*nxvh,kypd,kblok)
      dimension mixup(nxhyd), sctd(nxyd)
c local data
      integer indx1, indx1y, nx, nxh, nxhh, nx3, ny, nxy, nxhy, ks, kypt
      integer i, j, k, l, m, km, kmr, nrx, j1, j2, ns, ns2, k1, k2, kyps
      real at1, at2, t2, t3, t4, t5, t6, ani, sum1
      complex t1
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nx3 = nx + 3
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 2
      kyps = kypi + kypp - 1
      if (kstrt.gt.ny) return
      if (isign.eq.0) return
c create auxiliary array in x
      kmr = nxy/nx
      kypt = kyps
      do 30 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 20 k = kypi, kypt
      do 10 j = 2, nxh
      j1 = 1 + kmr*(j - 1)
      at2 = f(nx+2-j,k,l)
      at1 = f(j,k,l) + at2
      at2 = f(j,k,l) - at2
      at1 = -aimag(sctd(j1))*at1
      at2 = .5*at2
      f(j,k,l) = at1 + at2
      f(nx+2-j,k,l) = at1 - at2
   10 continue
      f(1,k,l) = 0.0
      f(nxh+1,k,l) = 2.0*f(nxh+1,k,l)
   20 continue
   30 continue
c bit-reverse array elements in x
      nrx = nxhy/nxh
      kypt = kyps
      do 60 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 50 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 50
      do 40 k = kypi, kypt
      t2 = f(2*j1-1,k,l)
      t3 = f(2*j1,k,l)
      f(2*j1-1,k,l) = f(2*j-1,k,l)
      f(2*j1,k,l) = f(2*j,k,l)
      f(2*j-1,k,l) = t2
      f(2*j,k,l) = t3
   40 continue
   50 continue
   60 continue
c first transform in x
      nrx = nxy/nxh
      do 110 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = 2*km*nrx
      kypt = kyps
      do 100 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 90 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 80 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 70 i = kypi, kypt
      t2 = real(t1)*f(2*j2-1,i,l) - aimag(t1)*f(2*j2,i,l)
      t3 = aimag(t1)*f(2*j2-1,i,l) + real(t1)*f(2*j2,i,l)
      f(2*j2-1,i,l) = f(2*j1-1,i,l) - t2
      f(2*j2,i,l) = f(2*j1,i,l) - t3
      f(2*j1-1,i,l) = f(2*j1-1,i,l) + t2
      f(2*j1,i,l) = f(2*j1,i,l) + t3
   70 continue
   80 continue
   90 continue
  100 continue
  110 continue
c unscramble coefficients and normalize
c inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nxh
         ani = 1./float(2*nx*ny)
         kypt = kyps
         do 150 l = 1, kblok
         if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
         do 130 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 120 k = kypi, kypt
         t4 = f(nx3-2*j,k,l)
         t5 = -f(nx3-2*j+1,k,l)
         t2 = f(2*j-1,k,l) + t4
         t3 = f(2*j,k,l) + t5
         t6 = f(2*j-1,k,l) - t4
         t5 = f(2*j,k,l) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(2*j-1,k,l) = ani*(t2 + t4)
         f(2*j,k,l) = ani*(t3 + t5)
         f(nx3-2*j,k,l) = ani*(t2 - t4)
         f(nx3-2*j+1,k,l) = ani*(t5 - t3)
  120    continue
  130    continue
         ani = 2.*ani
         do 140 k = kypi, kypt
         f(nxh+1,k,l) = ani*f(nxh+1,k,l)
         f(nxh+2,k,l) = -ani*f(nxh+2,k,l)
         t2 = ani*(f(1,k,l) + f(2,k,l))
         f(2,k,l) = ani*(f(1,k,l) - f(2,k,l))
         f(1,k,l) = t2
         f(nx+1,k,l) = ani*f(nx+1,k,l)
  140    continue
  150    continue
c forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nxh
         kypt = kyps
         do 190 l = 1, kblok
         if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
         do 170 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 160 k = kypi, kypt
         t4 = f(nx3-2*j,k,l)
         t5 = -f(nx3-2*j+1,k,l)
         t2 = f(2*j-1,k,l) + t4
         t3 = f(2*j,k,l) + t5
         t6 = f(2*j-1,k,l) - t4
         t5 = f(2*j,k,l) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(2*j-1,k,l) = t2 + t4
         f(2*j,k,l) = t3 + t5
         f(nx3-2*j,k,l) = t2 - t4
         f(nx3-2*j+1,k,l) = t5 - t3
  160    continue
  170    continue
         do 180 k = kypi, kypt
         f(nxh+1,k,l) = 2.0*f(nxh+1,k,l)
         f(nxh+2,k,l) = -2.0*f(nxh+2,k,l)
         t2 = 2.0*(f(1,k,l) + f(2,k,l))
         f(2,k,l) = 2.0*(f(1,k,l) - f(2,k,l))
         f(1,k,l) = t2
         f(nx+1,k,l) = 2.0*f(nx+1,k,l)
  180    continue
  190    continue
      endif
c perform recursion for sine transform
      kypt = kyps
      do 220 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 210 k = kypi, kypt
      sum1 = .5*f(1,k,l)
      f(1,k,l) = 0.0
      f(2,k,l) = sum1
      do 200 j = 2, nxh
      sum1 = sum1 + f(2*j-1,k,l)
      f(2*j-1,k,l) = -f(2*j,k,l)
      f(2*j,k,l) = sum1
  200 continue
      f(nx+1,k,l) = 0.0
  210 continue
  220 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PFCT2RXX(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,ky
     1pp,nxvh,kypd,kblok,nxhyd,nxyd)
c this subroutine performs the x part of a two dimensional fast real
c cosine transform and its inverse, for a subset of y,
c using real arithmetic, for data which is distributed in blocks
c algorithm is described in Numerical Recipies in Fortran, Second Ed.,
c by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
c [Cambridge Univ. Press, 1992], p. 508.
c for isign = (-1,1), input: all, output: f
c approximate flop count: N*(5*log2(N) + 18)/nvp
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, an inverse cosine transform is performed
c f(n,k,i) = (1/nx*ny)*(.5*f(1,k,i) + ((-1)**n)*f(nx+1,k,i)
c            + sum(f(j,k,i)*cos(pi*n*j/nx)))
c if isign = 1, a forward cosine transform is performed
c f(j,k,i) = 2*(.5*f(1,k,i) + ((-1)**j)*f(n+1,k,i) + sum(f(n,k,i)*
c            cos(pi*n*j/nx))
c mixup = array of bit reversed addresses
c sctd = sine/cosine table
c kstrt = starting data block number
c kyp = number of data values per block in y
c kypi = initial y index used
c kypp = number of y indices used
c nxvh = first dimension of f >= nx/2 + 1
c kypd = second dimension of f >= kyp+1
c kblok = number of data blocks in y
c nxhyd = maximum of (nx/2,ny)
c nxyd = maximum of (nx,ny)
c written by viktor k. decyk, ucla
      implicit none
      integer isign, mixup, indx, indy, kstrt, kyp, kypi, kypp
      integer nxvh, kypd, kblok, nxhyd, nxyd
      real f
      complex sctd
      dimension f(2*nxvh,kypd,kblok)
      dimension mixup(nxhyd), sctd(nxyd)
c local data
      integer indx1, indx1y, nx, nxh, nxhh, nx3, ny, nxy, nxhy, ks, kypt
      integer i, j, k, l, m, km, kmr, nrx, j1, j2, ns, ns2, k1, k2, kyps
      real at1, at2, t2, t3, t4, t5, t6, ani, sum1
      complex t1
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nx3 = nx + 3
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 2
      kyps = kypi + kypp - 1
      if (kstrt.gt.ny) return
      if (isign.eq.0) return
c create auxiliary array in x
      kmr = nxy/nx
      kypt = kyps
      do 30 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 20 k = kypi, kypt
      sum1 = .5*(f(1,k,l) - f(nx+1,k,l))
      do 10 j = 2, nxh
      j1 = 1 + kmr*(j - 1)
      at2 = f(nx+2-j,k,l)
      at1 = f(j,k,l) + at2
      at2 = f(j,k,l) - at2
      sum1 = sum1 + real(sctd(j1))*at2
      at2 = -aimag(sctd(j1))*at2
      at1 = .5*at1
      f(j,k,l) = at1 - at2
      f(nx+2-j,k,l) = at1 + at2
   10 continue
      f(1,k,l) = .5*(f(1,k,l) + f(nx+1,k,l))
      f(nx+1,k,l) = sum1
   20 continue
   30 continue
c bit-reverse array elements in x
      nrx = nxhy/nxh
      kypt = kyps
      do 60 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 50 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 50
      do 40 k = kypi, kypt
      t2 = f(2*j1-1,k,l)
      t3 = f(2*j1,k,l)
      f(2*j1-1,k,l) = f(2*j-1,k,l)
      f(2*j1,k,l) = f(2*j,k,l)
      f(2*j-1,k,l) = t2
      f(2*j,k,l) = t3
   40 continue
   50 continue
   60 continue
c first transform in x
      nrx = nxy/nxh
      do 110 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = 2*km*nrx
      kypt = kyps
      do 100 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 90 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 80 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 70 i = kypi, kypt
      t2 = real(t1)*f(2*j2-1,i,l) - aimag(t1)*f(2*j2,i,l)
      t3 = aimag(t1)*f(2*j2-1,i,l) + real(t1)*f(2*j2,i,l)
      f(2*j2-1,i,l) = f(2*j1-1,i,l) - t2
      f(2*j2,i,l) = f(2*j1,i,l) - t3
      f(2*j1-1,i,l) = f(2*j1-1,i,l) + t2
      f(2*j1,i,l) = f(2*j1,i,l) + t3
   70 continue
   80 continue
   90 continue
  100 continue
  110 continue
c unscramble coefficients and normalize
c inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nxh
         ani = 1./float(2*nx*ny)
         kypt = kyps
         do 150 l = 1, kblok
         if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
         do 130 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 120 k = kypi, kypt
         t4 = f(nx3-2*j,k,l)
         t5 = -f(nx3-2*j+1,k,l)
         t2 = f(2*j-1,k,l) + t4
         t3 = f(2*j,k,l) + t5
         t6 = f(2*j-1,k,l) - t4
         t5 = f(2*j,k,l) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(2*j-1,k,l) = ani*(t2 + t4)
         f(2*j,k,l) = ani*(t3 + t5)
         f(nx3-2*j,k,l) = ani*(t2 - t4)
         f(nx3-2*j+1,k,l) = ani*(t5 - t3)
  120    continue
  130    continue
         ani = 2.*ani
         do 140 k = kypi, kypt
         f(nxh+1,k,l) = ani*f(nxh+1,k,l)
         f(nxh+2,k,l) = -ani*f(nxh+2,k,l)
         t2 = ani*(f(1,k,l) + f(2,k,l))
         f(2,k,l) = ani*(f(1,k,l) - f(2,k,l))
         f(1,k,l) = t2
         f(nx+1,k,l) = ani*f(nx+1,k,l)
  140    continue
  150    continue
c forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nxh
         kypt = kyps
         do 190 l = 1, kblok
         if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
         do 170 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 160 k = kypi, kypt
         t4 = f(nx3-2*j,k,l)
         t5 = -f(nx3-2*j+1,k,l)
         t2 = f(2*j-1,k,l) + t4
         t3 = f(2*j,k,l) + t5
         t6 = f(2*j-1,k,l) - t4
         t5 = f(2*j,k,l) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(2*j-1,k,l) = t2 + t4
         f(2*j,k,l) = t3 + t5
         f(nx3-2*j,k,l) = t2 - t4
         f(nx3-2*j+1,k,l) = t5 - t3
  160    continue
  170    continue
         do 180 k = kypi, kypt
         f(nxh+1,k,l) = 2.0*f(nxh+1,k,l)
         f(nxh+2,k,l) = -2.0*f(nxh+2,k,l)
         t2 = 2.0*(f(1,k,l) + f(2,k,l))
         f(2,k,l) = 2.0*(f(1,k,l) - f(2,k,l))
         f(1,k,l) = t2
         f(nx+1,k,l) = 2.0*f(nx+1,k,l)
  180    continue
  190    continue
      endif
c perform recursion for cosine transform
      kypt = kyps
      do 220 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 210 k = kypi, kypt
      sum1 = f(nx+1,k,l)
      f(nx+1,k,l) = f(2,k,l)
      f(2,k,l) = sum1
      do 200 j = 2, nxh
      sum1 = sum1 - f(2*j,k,l)
      f(2*j,k,l) = sum1
  200 continue
  210 continue
  220 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PFST2RXY(g,isign,mixup,sctd,indx,indy,kstrt,kxp,kxpi,kx
     1pp,nyv,kxpd,jblok,nxhyd,nxyd)
c this subroutine performs the y part of a two dimensional fast real
c sine transform and its inverse, for a subset of x,
c using real arithmetic, for data which is distributed in blocks
c algorithm is described in Numerical Recipies in Fortran, Second Ed.,
c by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
c [Cambridge Univ. Press, 1992], p. 508.
c for isign = (-1,1), input: all, output: f
c approximate flop count: N*(5*log2(N) + 18)/nvp
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, an inverse sine transform is performed
c g(m,n,i) = sum(g(k,n,i)*sin(pi*m*k/ny))
c if isign = 1, a forward sine transform is performed
c g(k,n,i) = sum(g(m,n,i)*sin(pi*m*k/ny))
c mixup = array of bit reversed addresses
c sctd = sine/cosine table
c kstrt = starting data block number
c kxp = number of data values per block in x
c kxpi = initial x index used
c kxpp = number of x indices used
c nyv = first dimension of g >= ny + 1
c kxpd = second dimension of g >= kxp + 1
c jblok = number of data blocks in x
c nxhyd = maximum of (nx/2,ny)
c nxyd = maximum of (nx,ny)
c written by viktor k. decyk, ucla
      implicit none
      integer isign, mixup, indx, indy, kstrt, kxp, kxpi, kxpp
      integer nyv, kxpd, jblok, nxhyd, nxyd
      real g
      complex sctd
      dimension g(nyv,kxpd,jblok)
      dimension mixup(nxhyd), sctd(nxyd)
c local data
      integer indx1, indy1, indx1y, nx, ny, nyh, nyhh, ny3, nxy, nxhy
      integer i, j, k, l, m, ks, km, kmr, nry, j1, j2, ns, ns2, k1, k2
      integer kxps, kxpt
      real at1, at2, t2, t3, t4, t5, t6, ani, sum1
      complex t1
      indx1 = indx - 1
      indy1 = indy - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      nyhh = ny/4
      ny3 = ny + 3
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 2
      kxps = kxpi + kxpp - 1
      if (kstrt.gt.nx) return
      if (isign.eq.0) return
c create auxiliary array in y
      kmr = nxy/ny
      kxpt = kxps
      do 30 l = 1, jblok
      if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
      do 20 j = kxpi, kxpt
      do 10 k = 2, nyh
      k1 = 1 + kmr*(k - 1)
      at2 = g(ny+2-k,j,l)
      at1 = g(k,j,l) + at2
      at2 = g(k,j,l) - at2
      at1 = -aimag(sctd(k1))*at1
      at2 = .5*at2
      g(k,j,l) = at1 + at2
      g(ny+2-k,j,l) = at1 - at2
   10 continue
      g(1,j,l) = 0.0
      g(nyh+1,j,l) = 2.0*g(nyh+1,j,l)
   20 continue
   30 continue
c bit-reverse array elements in y
      nry = nxhy/nyh
      kxpt = kxps
      do 60 l = 1, jblok
      if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
      do 50 k = 1, nyh
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 50
      do 40 j = kxpi, kxpt
      t2 = g(2*k1-1,j,l)
      t3 = g(2*k1,j,l)
      g(2*k1-1,j,l) = g(2*k-1,j,l)
      g(2*k1,j,l) = g(2*k,j,l)
      g(2*k-1,j,l) = t2
      g(2*k,j,l) = t3
   40 continue
   50 continue
   60 continue
c first transform in y
      nry = nxy/nyh
      do 110 m = 1, indy1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyhh/ns
      kmr = 2*km*nry
      kxpt = kxps
      do 100 l = 1, jblok
      if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
      do 90 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 80 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 70 i = kxpi, kxpt
      t2 = real(t1)*g(2*j2-1,i,l) - aimag(t1)*g(2*j2,i,l)
      t3 = aimag(t1)*g(2*j2-1,i,l) + real(t1)*g(2*j2,i,l)
      g(2*j2-1,i,l) = g(2*j1-1,i,l) - t2
      g(2*j2,i,l) = g(2*j1,i,l) - t3
      g(2*j1-1,i,l) = g(2*j1-1,i,l) + t2
      g(2*j1,i,l) = g(2*j1,i,l) + t3
   70 continue
   80 continue
   90 continue
  100 continue
  110 continue
c unscramble coefficients and normalize
c inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nyh
         ani = 0.5
         kxpt = kxps
         do 150 l = 1, jblok
         if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
         do 130 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         do 120 j = kxpi, kxpt
         t4 = g(ny3-2*k,j,l)
         t5 = -g(ny3-2*k+1,j,l)
         t2 = g(2*k-1,j,l) + t4
         t3 = g(2*k,j,l) + t5
         t6 = g(2*k-1,j,l) - t4
         t5 = g(2*k,j,l) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         g(2*k-1,j,l) = ani*(t2 + t4)
         g(2*k,j,l) = ani*(t3 + t5)
         g(ny3-2*k,j,l) = ani*(t2 - t4)
         g(ny3-2*k+1,j,l) = ani*(t5 - t3)
  120    continue
  130    continue
         do 140 j = kxpi, kxpt
         g(nyh+1,j,l) = g(nyh+1,j,l)
         g(nyh+2,j,l) = -g(nyh+2,j,l)
         t2 = g(1,j,l) + g(2,j,l)
         g(2,j,l) = g(1,j,l) - g(2,j,l)
         g(1,j,l) = t2
         g(ny+1,j,l) = g(ny+1,j,l)
  140    continue
  150    continue
c forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nyh
         kxpt = kxps
         do 190 l = 1, jblok
         if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
         do 170 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         do 160 j = kxpi, kxpt
         t4 = g(ny3-2*k,j,l)
         t5 = -g(ny3-2*k+1,j,l)
         t2 = g(2*k-1,j,l) + t4
         t3 = g(2*k,j,l) + t5
         t6 = g(2*k-1,j,l) - t4
         t5 = g(2*k,j,l) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         g(2*k-1,j,l) = t2 + t4
         g(2*k,j,l) = t3 + t5
         g(ny3-2*k,j,l) = t2 - t4
         g(ny3-2*k+1,j,l) = t5 - t3
  160    continue
  170    continue
         do 180 j = kxpi, kxpt
         g(nyh+1,j,l) = 2.0*g(nyh+1,j,l)
         g(nyh+2,j,l) = -2.0*g(nyh+2,j,l)
         t2 = 2.0*(g(1,j,l) + g(2,j,l))
         g(2,j,l) = 2.0*(g(1,j,l) - g(2,j,l))
         g(1,j,l) = t2
         g(ny+1,j,l) = 2.0*g(ny+1,j,l)
  180    continue
  190    continue
      endif
c perform recursion for sine transform
      kxpt = kxps
      do 220 l = 1, jblok
      if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
      do 210 j = kxpi, kxpt
      sum1 = .5*g(1,j,l)
      g(1,j,l) = 0.0
      g(2,j,l) = sum1
      do 200 k = 2, nyh
      sum1 = sum1 + g(2*k-1,j,l)
      g(2*k-1,j,l) = -g(2*k,j,l)
      g(2*k,j,l) = sum1
  200 continue
      g(ny+1,j,l) = 0.0
  210 continue
  220 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PFCT2RXY(g,isign,mixup,sctd,indx,indy,kstrt,kxp,kxpi,kx
     1pp,nyv,kxpd,jblok,nxhyd,nxyd)
c this subroutine performs the y part of a two dimensional fast real
c cosine transform and its inverse, for a subset of x,
c using real arithmetic, for data which is distributed in blocks
c algorithm is described in Numerical Recipies in Fortran, Second Ed.,
c by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
c [Cambridge Univ. Press, 1992], p. 508.
c for isign = (-1,1), input: all, output: f
c approximate flop count: N*(5*log2(N) + 18)/nvp
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, an inverse cosine transform is performed
c g(m,n,i) = (.5*g(1,n,i) + ((-1)**m)*g(ny+1,n,i)
c            + sum(g(k,n,i)*cos(pi*m*k/ny))
c if isign = 1, a forward cosine transform is performed
c g(k,n,i) = 2*(.5*g(1,n,i) + ((-1)**m)*g(ny+1,n,i) + sum(g(m,n,i)*
c            cos(pi*m*k/ny))
c mixup = array of bit reversed addresses
c sctd = sine/cosine table
c kstrt = starting data block number
c kxp = number of data values per block in x
c kxpi = initial x index used
c kxpp = number of x indices used
c nyv = first dimension of g >= ny + 1
c kxpd = second dimension of g >= kxp + 1
c jblok = number of data blocks in x
c nxhyd = maximum of (nx/2,ny)
c nxyd = maximum of (nx,ny)
c written by viktor k. decyk, ucla
      implicit none
      integer isign, mixup, indx, indy, kstrt, kxp, kxpi, kxpp
      integer nyv, kxpd, jblok, nxhyd, nxyd
      real g
      complex sctd
      dimension g(nyv,kxpd,jblok)
      dimension mixup(nxhyd), sctd(nxyd)
c local data
      integer indx1, indy1, indx1y, nx, ny, nyh, nyhh, ny3, nxy, nxhy
      integer i, j, k, l, m, ks, km, kmr, nry, j1, j2, ns, ns2, k1, k2
      integer kxps, kxpt
      real at1, at2, t2, t3, t4, t5, t6, ani, sum1
      complex t1
      indx1 = indx - 1
      indy1 = indy - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      nyhh = ny/4
      ny3 = ny + 3
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 2
      kxps = kxpi + kxpp - 1
      if (kstrt.gt.nx) return
      if (isign.eq.0) return
c create auxiliary array in y
      kmr = nxy/ny
      kxpt = kxps
      do 30 l = 1, jblok
      if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
      do 20 j = kxpi, kxpt
      sum1 = .5*(g(1,j,l) - g(ny+1,j,l))
      do 10 k = 2, nyh
      k1 = 1 + kmr*(k - 1)
      at2 = g(ny+2-k,j,l)
      at1 = g(k,j,l) + at2
      at2 = g(k,j,l) - at2
      sum1 = sum1 + real(sctd(k1))*at2
      at2 = -aimag(sctd(k1))*at2
      at1 = .5*at1
      g(k,j,l) = at1 - at2
      g(ny+2-k,j,l) = at1 + at2
   10 continue
      g(1,j,l) = .5*(g(1,j,l) + g(ny+1,j,l))
      g(ny+1,j,l) = sum1
   20 continue
   30 continue
c bit-reverse array elements in y
      nry = nxhy/nyh
      kxpt = kxps
      do 60 l = 1, jblok
      if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
      do 50 k = 1, nyh
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 50
      do 40 j = kxpi, kxpt
      t2 = g(2*k1-1,j,l)
      t3 = g(2*k1,j,l)
      g(2*k1-1,j,l) = g(2*k-1,j,l)
      g(2*k1,j,l) = g(2*k,j,l)
      g(2*k-1,j,l) = t2
      g(2*k,j,l) = t3
   40 continue
   50 continue
   60 continue
c first transform in y
      nry = nxy/nyh
      do 110 m = 1, indy1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyhh/ns
      kmr = 2*km*nry
      kxpt = kxps
      do 100 l = 1, jblok
      if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
      do 90 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 80 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 70 i = kxpi, kxpt
      t2 = real(t1)*g(2*j2-1,i,l) - aimag(t1)*g(2*j2,i,l)
      t3 = aimag(t1)*g(2*j2-1,i,l) + real(t1)*g(2*j2,i,l)
      g(2*j2-1,i,l) = g(2*j1-1,i,l) - t2
      g(2*j2,i,l) = g(2*j1,i,l) - t3
      g(2*j1-1,i,l) = g(2*j1-1,i,l) + t2
      g(2*j1,i,l) = g(2*j1,i,l) + t3
   70 continue
   80 continue
   90 continue
  100 continue
  110 continue
c unscramble coefficients and normalize
c inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nyh
         ani = 0.5
         kxpt = kxps
         do 150 l = 1, jblok
         if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
         do 130 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         do 120 j = kxpi, kxpt
         t4 = g(ny3-2*k,j,l)
         t5 = -g(ny3-2*k+1,j,l)
         t2 = g(2*k-1,j,l) + t4
         t3 = g(2*k,j,l) + t5
         t6 = g(2*k-1,j,l) - t4
         t5 = g(2*k,j,l) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         g(2*k-1,j,l) = ani*(t2 + t4)
         g(2*k,j,l) = ani*(t3 + t5)
         g(ny3-2*k,j,l) = ani*(t2 - t4)
         g(ny3-2*k+1,j,l) = ani*(t5 - t3)
  120    continue
  130    continue
         do 140 j = kxpi, kxpt
         g(nyh+1,j,l) = g(nyh+1,j,l)
         g(nyh+2,j,l) = -g(nyh+2,j,l)
         t2 = g(1,j,l) + g(2,j,l)
         g(2,j,l) = g(1,j,l) - g(2,j,l)
         g(1,j,l) = t2
         g(ny+1,j,l) = g(ny+1,j,l)
  140    continue
  150    continue
c forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nyh
         kxpt = kxps
         do 190 l = 1, jblok
         if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
         do 170 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         do 160 j = kxpi, kxpt
         t4 = g(ny3-2*k,j,l)
         t5 = -g(ny3-2*k+1,j,l)
         t2 = g(2*k-1,j,l) + t4
         t3 = g(2*k,j,l) + t5
         t6 = g(2*k-1,j,l) - t4
         t5 = g(2*k,j,l) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         g(2*k-1,j,l) = t2 + t4
         g(2*k,j,l) = t3 + t5
         g(ny3-2*k,j,l) = t2 - t4
         g(ny3-2*k+1,j,l) = t5 - t3
  160    continue
  170    continue
         do 180 j = kxpi, kxpt
         g(nyh+1,j,l) = 2.0*g(nyh+1,j,l)
         g(nyh+2,j,l) = -2.0*g(nyh+2,j,l)
         t2 = 2.0*(g(1,j,l) + g(2,j,l))
         g(2,j,l) = 2.0*(g(1,j,l) - g(2,j,l))
         g(1,j,l) = t2
         g(ny+1,j,l) = 2.0*g(ny+1,j,l)
  180    continue
  190    continue
      endif
c perform recursion for cosine transform
      kxpt = kxps
      do 220 l = 1, jblok
      if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
      do 210 j = kxpi, kxpt
      sum1 = g(ny+1,j,l)
      g(ny+1,j,l) = g(2,j,l)
      g(2,j,l) = sum1
      do 200 k = 2, nyh
      sum1 = sum1 - g(2*k,j,l)
      g(2*k,j,l) = sum1
  200 continue
  210 continue
  220 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine WPFCST2R2(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,in
     1dy,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhyd,nxyd)
c wrapper function for 2 parallel real sine/sine transforms
c for the electric field with dirichlet or magnetic field with neumann
c boundary conditions
      implicit none
      integer isign, ntpose, mixup, indx, indy, kstrt, nxvh, nyv
      integer kxp2, kyp, kypd, kxp2d, jblok, kblok, nxhyd, nxyd
      real f, g, bs, br, ttp
      complex sctd
      dimension f(2,2*nxvh,kypd,kblok), g(2,nyv,kxp2d,jblok)
      dimension bs(2,kxp2+1,kyp+1,kblok), br(2,kxp2+1,kyp+1,jblok)
      dimension mixup(nxhyd), sctd(nxyd)
c local data
      integer nx, ny, kxpi, kypi
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
c inverse fourier transform
      if (isign.lt.0) then
c perform x cosine-sine transform
         call PFCST2R2X(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kyp,
     1nxvh,kypd,kblok,nxhyd,nxyd)
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PR2TPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,nyv,kxp2,kyp,kxp2d,k
     1ypd,jblok,kblok)
         call PWTIMERA(1,ttp,dtime)
c perform y sine-cosine transform
         call PFSCT2R2Y(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxp
     12,nyv,kxp2d,jblok,nxhyd,nxyd)
c transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PR2TPOSE(g,f,br,bs,ny,nx,kstrt,nyv,2*nxvh,kyp,kxp2,kypd
     1,kxp2d,kblok,jblok)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PR2TPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,nyv,kxp2,kyp,kxp2
     1d,kypd,jblok,kblok)
            call PWTIMERA(1,tf,dtime)
         endif
c perform y sine-cosine transform
         call PFSCT2R2Y(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxp
     12,nyv,kxp2d,jblok,nxhyd,nxyd)
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PR2TPOSE(g,f,br,bs,ny,nx,kstrt,nyv,2*nxvh,kyp,kxp2,kypd,kx
     1p2d,kblok,jblok)
         call PWTIMERA(1,ttp,dtime)
c perform x cosine-sine transform
         call PFCST2R2X(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kyp,
     1nxvh,kypd,kblok,nxhyd,nxyd)
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine WPFSCT2R2(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,in
     1dy,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhyd,nxyd)
c wrapper function for 2 parallel real sine/sine transforms
c for the magnetic field with dirichlet or electric field with neumann
c boundary conditions
      implicit none
      integer isign, ntpose, mixup, indx, indy, kstrt, nxvh, nyv
      integer kxp2, kyp, kypd, kxp2d, jblok, kblok, nxhyd, nxyd
      real f, g, bs, br, ttp
      complex sctd
      dimension f(2,2*nxvh,kypd,kblok), g(2,nyv,kxp2d,jblok)
      dimension bs(2,kxp2+1,kyp+1,kblok), br(2,kxp2+1,kyp+1,jblok)
      dimension mixup(nxhyd), sctd(nxyd)
c local data
      integer nx, ny, kxpi, kypi
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
c inverse fourier transform
      if (isign.lt.0) then
c perform x sine-cosine transform
         call PFSCT2R2X(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kyp,
     1nxvh,kypd,kblok,nxhyd,nxyd)
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PR2TPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,nyv,kxp2,kyp,kxp2d,k
     1ypd,jblok,kblok)
         call PWTIMERA(1,ttp,dtime)
c perform y cosine-sine transform
         call PFCST2R2Y(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxp
     12,nyv,kxp2d,jblok,nxhyd,nxyd)
c transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PR2TPOSE(g,f,br,bs,ny,nx,kstrt,nyv,2*nxvh,kyp,kxp2,kypd
     1,kxp2d,kblok,jblok)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PR2TPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,nyv,kxp2,kyp,kxp2
     1d,kypd,jblok,kblok)
            call PWTIMERA(1,tf,dtime)
         endif
c perform y cosine-sine transform
         call PFCST2R2Y(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxp
     12,nyv,kxp2d,jblok,nxhyd,nxyd)
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PR2TPOSE(g,f,br,bs,ny,nx,kstrt,nyv,2*nxvh,kyp,kxp2,kypd,kx
     1p2d,kblok,jblok)
         call PWTIMERA(1,ttp,dtime)
c perform x sine-cosine transform
         call PFSCT2R2X(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kyp,
     1nxvh,kypd,kblok,nxhyd,nxyd)
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine PFCST2R2X(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,k
     1ypp,nxvh,kypd,kblok,nxhyd,nxyd)
c this subroutine performs the x part of 2 two dimensional fast real
c sine and cosine transforms and their inverses, for a subset of y,
c using real arithmetic, for data which is distributed in blocks
c algorithm is described in Numerical Recipies in Fortran, Second Ed.,
c by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
c [Cambridge Univ. Press, 1992], p. 508.
c for isign = (-1,1), input: all, output: f
c approximate flop count: N*(5*log2(N) + 18)/nvp
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, inverse sine-cosine transforms are performed
c f(1,n,k,i) = (1/nx*ny)*(.5*f(1,1,k,i) + ((-1)**n)*f(1,nx+1,k,i)
c              + sum(f(1,j,k,i)*cos(pi*n*j/nx)))
c f(2,n,k,i) = (1/nx*ny)*sum(f(2,j,k,i)*sin(pi*n*j/nx))
c if isign = 1, forward sine transforms are performed
c f(1,j,k,i) = 2*(.5*f(1,1,k,i) + ((-1)**j)*f(1,n+1,k,i)
c              + sum(f(1,n,k,i)*cos(pi*n*j/nx))
c f(2,j,k,i) = sum(f(2,n,k,i)*sin(pi*n*j/nx))
c mixup = array of bit reversed addresses
c sctd = sine/cosine table
c kstrt = starting data block number
c kyp = number of data values per block in y
c kypi = initial y index used
c kypp = number of y indices used
c nxvh = first dimension of f >= nx/2 + 1
c kypd = second dimension of f >= kyp + 1
c kblok = number of data blocks in y
c nxhyd = maximum of (nx/2,ny)
c nxyd = maximum of (nx,ny)
c written by viktor k. decyk, ucla
      implicit none
      integer isign, mixup, indx, indy, kstrt, kyp, kypi, kypp
      integer nxvh, kypd, kblok, nxhyd, nxyd
      real f
      complex sctd
      dimension f(2,2*nxvh,kypd,kblok)
      dimension mixup(nxhyd), sctd(nxyd)
c local data
      integer indx1, indx1y, nx, nxh, nxhh, nx3, ny, nxy, nxhy, ks, kypt
      integer i, j, k, l, m, km, kmr, nrx, j1, j2, ns, ns2, k1, k2, kyps
      integer jj
      real at1, at2, at3, t2, t3, t4, t5, t6, ani, sum1, sum2
      complex t1
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nx3 = nx + 3
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 2
      kyps = kypi + kypp - 1
      if (kstrt.gt.ny) return
      if (isign.eq.0) return
c create auxiliary array in x
      kmr = nxy/nx
      kypt = kyps
      do 30 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 20 k = kypi, kypt
      sum1 = .5*(f(1,1,k,l) - f(1,nx+1,k,l))
      do 10 j = 2, nxh
      j1 = 1 + kmr*(j - 1)
      at3 = -aimag(sctd(j1))
      at2 = f(1,nx+2-j,k,l)
      at1 = f(1,j,k,l) + at2
      at2 = f(1,j,k,l) - at2
      sum1 = sum1 + real(sctd(j1))*at2
      at2 = at3*at2
      at1 = .5*at1
      f(1,j,k,l) = at1 - at2
      f(1,nx+2-j,k,l) = at1 + at2
      at2 = f(2,nx+2-j,k,l)
      at1 = f(2,j,k,l) + at2
      at2 = f(2,j,k,l) - at2
      at1 = at3*at1
      at2 = .5*at2
      f(2,j,k,l) = at1 + at2
      f(2,nx+2-j,k,l) = at1 - at2
   10 continue
      f(1,1,k,l) = .5*(f(1,1,k,l) + f(1,nx+1,k,l))
      f(1,nx+1,k,l) = sum1
      f(2,1,k,l) = 0.0
      f(2,nxh+1,k,l) = 2.0*f(2,nxh+1,k,l)
   20 continue
   30 continue
c bit-reverse array elements in x
      nrx = nxhy/nxh
      kypt = kyps
      do 70 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 60 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 60
      do 50 k = kypi, kypt
      do 40 jj = 1, 2
      t2 = f(jj,2*j1-1,k,l)
      t3 = f(jj,2*j1,k,l)
      f(jj,2*j1-1,k,l) = f(jj,2*j-1,k,l)
      f(jj,2*j1,k,l) = f(jj,2*j,k,l)
      f(jj,2*j-1,k,l) = t2
      f(jj,2*j,k,l) = t3
   40 continue
   50 continue
   60 continue
   70 continue
c first transform in x
      nrx = nxy/nxh
      do 130 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = 2*km*nrx
      kypt = kyps
      do 120 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 110 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 100 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 90 i = kypi, kypt
      do 80 jj = 1, 2
      t2 = real(t1)*f(jj,2*j2-1,i,l) - aimag(t1)*f(jj,2*j2,i,l)
      t3 = aimag(t1)*f(jj,2*j2-1,i,l) + real(t1)*f(jj,2*j2,i,l)
      f(jj,2*j2-1,i,l) = f(jj,2*j1-1,i,l) - t2
      f(jj,2*j2,i,l) = f(jj,2*j1,i,l) - t3
      f(jj,2*j1-1,i,l) = f(jj,2*j1-1,i,l) + t2
      f(jj,2*j1,i,l) = f(jj,2*j1,i,l) + t3
   80 continue
   90 continue
  100 continue
  110 continue
  120 continue
  130 continue
c unscramble coefficients and normalize
c inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nxh
         ani = 1./float(2*nx*ny)
         kypt = kyps
         do 190 l = 1, kblok
         if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
         do 160 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 150 k = kypi, kypt
         do 140 jj = 1, 2
         t4 = f(jj,nx3-2*j,k,l)
         t5 = -f(jj,nx3-2*j+1,k,l)
         t2 = f(jj,2*j-1,k,l) + t4
         t3 = f(jj,2*j,k,l) + t5
         t6 = f(jj,2*j-1,k,l) - t4
         t5 = f(jj,2*j,k,l) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,2*j-1,k,l) = ani*(t2 + t4)
         f(jj,2*j,k,l) = ani*(t3 + t5)
         f(jj,nx3-2*j,k,l) = ani*(t2 - t4)
         f(jj,nx3-2*j+1,k,l) = ani*(t5 - t3)
  140    continue
  150    continue
  160    continue
         ani = 2.*ani
         do 180 k = kypi, kypt
         do 170 jj = 1, 2
         f(jj,nxh+1,k,l) = ani*f(jj,nxh+1,k,l)
         f(jj,nxh+2,k,l) = -ani*f(jj,nxh+2,k,l)
         t2 = ani*(f(jj,1,k,l) + f(jj,2,k,l))
         f(jj,2,k,l) = ani*(f(jj,1,k,l) - f(jj,2,k,l))
         f(jj,1,k,l) = t2
         f(jj,nx+1,k,l) = ani*f(jj,nx+1,k,l)
  170    continue
  180    continue
  190    continue
c forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nxh
         kypt = kyps
         do 250 l = 1, kblok
         if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
         do 220 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 210 k = kypi, kypt
         do 200 jj = 1, 2
         t4 = f(jj,nx3-2*j,k,l)
         t5 = -f(jj,nx3-2*j+1,k,l)
         t2 = f(jj,2*j-1,k,l) + t4
         t3 = f(jj,2*j,k,l) + t5
         t6 = f(jj,2*j-1,k,l) - t4
         t5 = f(jj,2*j,k,l) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,2*j-1,k,l) = t2 + t4
         f(jj,2*j,k,l) = t3 + t5
         f(jj,nx3-2*j,k,l) = t2 - t4
         f(jj,nx3-2*j+1,k,l) = t5 - t3
  200    continue
  210    continue
  220    continue
         do 240 k = kypi, kypt
         do 230 jj = 1, 2
         f(jj,nxh+1,k,l) = 2.0*f(jj,nxh+1,k,l)
         f(jj,nxh+2,k,l) = -2.0*f(jj,nxh+2,k,l)
         t2 = 2.0*(f(jj,1,k,l) + f(jj,2,k,l))
         f(jj,2,k,l) = 2.0*(f(jj,1,k,l) - f(jj,2,k,l))
         f(jj,1,k,l) = t2
         f(jj,nx+1,k,l) = 2.0*f(jj,nx+1,k,l)
  230    continue
  240    continue
  250    continue
      endif
c perform recursion for cosine-sine transform
      kypt = kyps
      do 280 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 270 k = kypi, kypt
      sum1 = f(1,nx+1,k,l)
      f(1,nx+1,k,l) = f(1,2,k,l)
      f(1,2,k,l) = sum1
      sum2 = .5*f(2,1,k,l)
      f(2,1,k,l) = 0.0
      f(2,2,k,l) = sum2
      do 260 j = 2, nxh
      sum1 = sum1 - f(1,2*j,k,l)
      f(1,2*j,k,l) = sum1
      sum2 = sum2 + f(2,2*j-1,k,l)
      f(2,2*j-1,k,l) = -f(2,2*j,k,l)
      f(2,2*j,k,l) = sum2
  260 continue
      f(2,nx+1,k,l) = 0.0
  270 continue
  280 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PFSCT2R2X(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,k
     1ypp,nxvh,kypd,kblok,nxhyd,nxyd)
c this subroutine performs the x part of 2 two dimensional fast real
c sine and cosine transforms and their inverses, for a subset of y,
c using real arithmetic, for data which is distributed in blocks
c algorithm is described in Numerical Recipies in Fortran, Second Ed.,
c by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
c [Cambridge Univ. Press, 1992], p. 508.
c for isign = (-1,1), input: all, output: f
c approximate flop count: N*(5*log2(N) + 18)/nvp
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, inverse sine-cosine transforms are performed
c f(1,n,k,i) = (1/nx*ny)*sum(f(1,j,k,i)*sin(pi*n*j/nx))
c f(2,n,k,i) = (1/nx*ny)*(.5*f(2,1,k,i) + ((-1)**n)*f(2,nx+1,k,i)
c              + sum(f(2,j,k,i)*cos(pi*n*j/nx)))
c if isign = 1, forward sine transforms are performed
c f(1,j,k,i) = sum(f(1,n,k,i)*sin(pi*n*j/nx))
c f(2,j,k,i) = 2*(.5*f(2,1,k,i) + ((-1)**j)*f(2,n+1,k,i)
c              + sum(f(2,n,k,i)*cos(pi*n*j/nx))
c mixup = array of bit reversed addresses
c sctd = sine/cosine table
c kstrt = starting data block number
c kyp = number of data values per block in y
c kypi = initial y index used
c kypp = number of y indices used
c nxvh = first dimension of f >= nx/2 + 1
c kypd = second dimension of f >= kyp + 1
c kblok = number of data blocks in y
c nxhyd = maximum of (nx/2,ny)
c nxyd = maximum of (nx,ny)
c written by viktor k. decyk, ucla
      implicit none
      integer isign, mixup, indx, indy, kstrt, kyp, kypi, kypp
      integer nxvh, kypd, kblok, nxhyd, nxyd
      real f
      complex sctd
      dimension f(2,2*nxvh,kypd,kblok)
      dimension mixup(nxhyd), sctd(nxyd)
c local data
      integer indx1, indx1y, nx, nxh, nxhh, nx3, ny, nxy, nxhy, ks, kypt
      integer i, j, k, l, m, km, kmr, nrx, j1, j2, ns, ns2, k1, k2, kyps
      integer jj
      real at1, at2, at3, t2, t3, t4, t5, t6, ani, sum1, sum2
      complex t1
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nx3 = nx + 3
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 2
      kyps = kypi + kypp - 1
      if (kstrt.gt.ny) return
      if (isign.eq.0) return
c create auxiliary array in x
      kmr = nxy/nx
      kypt = kyps
      do 30 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 20 k = kypi, kypt
      sum1 = .5*(f(2,1,k,l) - f(2,nx+1,k,l))
      do 10 j = 2, nxh
      j1 = 1 + kmr*(j - 1)
      at3 = -aimag(sctd(j1))
      at2 = f(1,nx+2-j,k,l)
      at1 = f(1,j,k,l) + at2
      at2 = f(1,j,k,l) - at2
      at1 = at3*at1
      at2 = .5*at2
      f(1,j,k,l) = at1 + at2
      f(1,nx+2-j,k,l) = at1 - at2
      at2 = f(2,nx+2-j,k,l)
      at1 = f(2,j,k,l) + at2
      at2 = f(2,j,k,l) - at2
      sum1 = sum1 + real(sctd(j1))*at2
      at2 = at3*at2
      at1 = .5*at1
      f(2,j,k,l) = at1 - at2
      f(2,nx+2-j,k,l) = at1 + at2
   10 continue
      f(1,1,k,l) = 0.0
      f(1,nxh+1,k,l) = 2.0*f(1,nxh+1,k,l)
      f(2,1,k,l) = .5*(f(2,1,k,l) + f(2,nx+1,k,l))
      f(2,nx+1,k,l) = sum1
   20 continue
   30 continue
c bit-reverse array elements in x
      nrx = nxhy/nxh
      kypt = kyps
      do 70 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 60 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 60
      do 50 k = kypi, kypt
      do 40 jj = 1, 2
      t2 = f(jj,2*j1-1,k,l)
      t3 = f(jj,2*j1,k,l)
      f(jj,2*j1-1,k,l) = f(jj,2*j-1,k,l)
      f(jj,2*j1,k,l) = f(jj,2*j,k,l)
      f(jj,2*j-1,k,l) = t2
      f(jj,2*j,k,l) = t3
   40 continue
   50 continue
   60 continue
   70 continue
c first transform in x
      nrx = nxy/nxh
      do 130 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = 2*km*nrx
      kypt = kyps
      do 120 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 110 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 100 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 90 i = kypi, kypt
      do 80 jj = 1, 2
      t2 = real(t1)*f(jj,2*j2-1,i,l) - aimag(t1)*f(jj,2*j2,i,l)
      t3 = aimag(t1)*f(jj,2*j2-1,i,l) + real(t1)*f(jj,2*j2,i,l)
      f(jj,2*j2-1,i,l) = f(jj,2*j1-1,i,l) - t2
      f(jj,2*j2,i,l) = f(jj,2*j1,i,l) - t3
      f(jj,2*j1-1,i,l) = f(jj,2*j1-1,i,l) + t2
      f(jj,2*j1,i,l) = f(jj,2*j1,i,l) + t3
   80 continue
   90 continue
  100 continue
  110 continue
  120 continue
  130 continue
c unscramble coefficients and normalize
c inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nxh
         ani = 1./float(2*nx*ny)
         kypt = kyps
         do 190 l = 1, kblok
         if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
         do 160 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 150 k = kypi, kypt
         do 140 jj = 1, 2
         t4 = f(jj,nx3-2*j,k,l)
         t5 = -f(jj,nx3-2*j+1,k,l)
         t2 = f(jj,2*j-1,k,l) + t4
         t3 = f(jj,2*j,k,l) + t5
         t6 = f(jj,2*j-1,k,l) - t4
         t5 = f(jj,2*j,k,l) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,2*j-1,k,l) = ani*(t2 + t4)
         f(jj,2*j,k,l) = ani*(t3 + t5)
         f(jj,nx3-2*j,k,l) = ani*(t2 - t4)
         f(jj,nx3-2*j+1,k,l) = ani*(t5 - t3)
  140    continue
  150    continue
  160    continue
         ani = 2.*ani
         do 180 k = kypi, kypt
         do 170 jj = 1, 2
         f(jj,nxh+1,k,l) = ani*f(jj,nxh+1,k,l)
         f(jj,nxh+2,k,l) = -ani*f(jj,nxh+2,k,l)
         t2 = ani*(f(jj,1,k,l) + f(jj,2,k,l))
         f(jj,2,k,l) = ani*(f(jj,1,k,l) - f(jj,2,k,l))
         f(jj,1,k,l) = t2
         f(jj,nx+1,k,l) = ani*f(jj,nx+1,k,l)
  170    continue
  180    continue
  190    continue
c forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nxh
         kypt = kyps
         do 250 l = 1, kblok
         if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
         do 220 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 210 k = kypi, kypt
         do 200 jj = 1, 2
         t4 = f(jj,nx3-2*j,k,l)
         t5 = -f(jj,nx3-2*j+1,k,l)
         t2 = f(jj,2*j-1,k,l) + t4
         t3 = f(jj,2*j,k,l) + t5
         t6 = f(jj,2*j-1,k,l) - t4
         t5 = f(jj,2*j,k,l) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,2*j-1,k,l) = t2 + t4
         f(jj,2*j,k,l) = t3 + t5
         f(jj,nx3-2*j,k,l) = t2 - t4
         f(jj,nx3-2*j+1,k,l) = t5 - t3
  200    continue
  210    continue
  220    continue
         do 240 k = kypi, kypt
         do 230 jj = 1, 2
         f(jj,nxh+1,k,l) = 2.0*f(jj,nxh+1,k,l)
         f(jj,nxh+2,k,l) = -2.0*f(jj,nxh+2,k,l)
         t2 = 2.0*(f(jj,1,k,l) + f(jj,2,k,l))
         f(jj,2,k,l) = 2.0*(f(jj,1,k,l) - f(jj,2,k,l))
         f(jj,1,k,l) = t2
         f(jj,nx+1,k,l) = 2.0*f(jj,nx+1,k,l)
  230    continue
  240    continue
  250    continue
      endif
c perform recursion for cosine-sine transform
      kypt = kyps
      do 280 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 270 k = kypi, kypt
      sum1 = .5*f(1,1,k,l)
      f(1,1,k,l) = 0.0
      f(1,2,k,l) = sum1
      sum2 = f(2,nx+1,k,l)
      f(2,nx+1,k,l) = f(2,2,k,l)
      f(2,2,k,l) = sum2
      do 260 j = 2, nxh
      sum1 = sum1 + f(1,2*j-1,k,l)
      f(1,2*j-1,k,l) = -f(1,2*j,k,l)
      f(1,2*j,k,l) = sum1
      sum2 = sum2 - f(2,2*j,k,l)
      f(2,2*j,k,l) = sum2
  260 continue
      f(1,nx+1,k,l) = 0.0
  270 continue
  280 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PFSCT2R2Y(g,isign,mixup,sctd,indx,indy,kstrt,kxp,kxpi,k
     1xpp,nyv,kxpd,jblok,nxhyd,nxyd)
c this subroutine performs the y part of 2 two dimensional fast real
c sine and cosine transforms and their inverses, for a subset of x,
c using real arithmetic, for data which is distributed in blocks
c algorithm is described in Numerical Recipies in Fortran, Second Ed.,
c by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
c [Cambridge Univ. Press, 1992], p. 508.
c for isign = (-1,1), input: all, output: f
c approximate flop count: N*(5*log2(N) + 18)/nvp
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, inverse sine-cosine transform are performed
c g(1,m,n,i) = sum(g(1,k,n,i)*sin(pi*m*k/ny))
c g(2,m,n,i) = (.5*g(2,1,n,i) + ((-1)**m)*g(2,ny+1,n,i)
c              + sum(g(2,k,n,i)*cos(pi*m*k/ny))
c if isign = 1, a forward sine-cosine transforms are performed
c g(1,k,n,i) = sum(g(1,m,n,i)*sin(pi*m*k/ny))
c g(2,k,n,i) = 2*(.5*g(2,1,n,i) + ((-1)**m)*g(2,ny+1,n,i)
c              + sum(g(2,m,n,i)*cos(pi*m*k/ny))
c mixup = array of bit reversed addresses
c sctd = sine/cosine table
c kstrt = starting data block number
c kxp = number of data values per block in x
c kxpi = initial x index used
c kxpp = number of x indices used
c nyv = first dimension of g >= ny + 1
c kxpd = second dimension of g >= kxp + 1
c jblok = number of data blocks in x
c nxhyd = maximum of (nx/2,ny)
c nxyd = maximum of (nx,ny)
c written by viktor k. decyk, ucla
      implicit none
      integer isign, mixup, indx, indy, kstrt, kxp, kxpi, kxpp
      integer nyv, kxpd, jblok, nxhyd, nxyd
      real g
      complex sctd
      dimension g(2,nyv,kxpd,jblok)
      dimension mixup(nxhyd), sctd(nxyd)
c local data
      integer indx1, indy1, indx1y, nx, ny, nyh, nyhh, ny3, nxy, nxhy
      integer i, j, k, l, m, ks, km, kmr, nry, j1, j2, ns, ns2, k1, k2
      integer kxps, kxpt, jj
      real at1, at2, at3, t2, t3, t4, t5, t6, ani, sum1, sum2
      complex t1
      indx1 = indx - 1
      indy1 = indy - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      nyhh = ny/4
      ny3 = ny + 3
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 2
      kxps = kxpi + kxpp - 1
      if (kstrt.gt.nx) return
      if (isign.eq.0) return
c create auxiliary array in y
      kmr = nxy/ny
      kxpt = kxps
      do 30 l = 1, jblok
      if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
      do 20 j = kxpi, kxpt
      sum1 = .5*(g(2,1,j,l) - g(2,ny+1,j,l))
      do 10 k = 2, nyh
      k1 = 1 + kmr*(k - 1)
      at3 = -aimag(sctd(k1))
      at2 = g(1,ny+2-k,j,l)
      at1 = g(1,k,j,l) + at2
      at2 = g(1,k,j,l) - at2
      at1 = at3*at1
      at2 = .5*at2
      g(1,k,j,l) = at1 + at2
      g(1,ny+2-k,j,l) = at1 - at2
      at2 = g(2,ny+2-k,j,l)
      at1 = g(2,k,j,l) + at2
      at2 = g(2,k,j,l) - at2
      sum1 = sum1 + real(sctd(k1))*at2
      at2 = at3*at2
      at1 = .5*at1
      g(2,k,j,l) = at1 - at2
      g(2,ny+2-k,j,l) = at1 + at2
   10 continue
      g(1,1,j,l) = 0.0
      g(1,nyh+1,j,l) = 2.0*g(1,nyh+1,j,l)
      g(2,1,j,l) = .5*(g(2,1,j,l) + g(2,ny+1,j,l))
      g(2,ny+1,j,l) = sum1
   20 continue
   30 continue
c bit-reverse array elements in y
      nry = nxhy/nyh
      kxpt = kxps
      do 70 l = 1, jblok
      if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
      do 60 k = 1, nyh
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 60
      do 50 j = kxpi, kxpt
      do 40 jj = 1, 2
      t2 = g(jj,2*k1-1,j,l)
      t3 = g(jj,2*k1,j,l)
      g(jj,2*k1-1,j,l) = g(jj,2*k-1,j,l)
      g(jj,2*k1,j,l) = g(jj,2*k,j,l)
      g(jj,2*k-1,j,l) = t2
      g(jj,2*k,j,l) = t3
   40 continue
   50 continue
   60 continue
   70 continue
c first transform in y
      nry = nxy/nyh
      do 130 m = 1, indy1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyhh/ns
      kmr = 2*km*nry
      kxpt = kxps
      do 120 l = 1, jblok
      if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
      do 110 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 100 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 90 i = kxpi, kxpt
      do 80 jj = 1, 2
      t2 = real(t1)*g(jj,2*j2-1,i,l) - aimag(t1)*g(jj,2*j2,i,l)
      t3 = aimag(t1)*g(jj,2*j2-1,i,l) + real(t1)*g(jj,2*j2,i,l)
      g(jj,2*j2-1,i,l) = g(jj,2*j1-1,i,l) - t2
      g(jj,2*j2,i,l) = g(jj,2*j1,i,l) - t3
      g(jj,2*j1-1,i,l) = g(jj,2*j1-1,i,l) + t2
      g(jj,2*j1,i,l) = g(jj,2*j1,i,l) + t3
   80 continue
   90 continue
  100 continue
  110 continue
  120 continue
  130 continue
c unscramble coefficients and normalize
c inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nyh
         ani = 0.5
         kxpt = kxps
         do 190 l = 1, jblok
         if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
         do 160 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         do 150 j = kxpi, kxpt
         do 140 jj = 1, 2
         t4 = g(jj,ny3-2*k,j,l)
         t5 = -g(jj,ny3-2*k+1,j,l)
         t2 = g(jj,2*k-1,j,l) + t4
         t3 = g(jj,2*k,j,l) + t5
         t6 = g(jj,2*k-1,j,l) - t4
         t5 = g(jj,2*k,j,l) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         g(jj,2*k-1,j,l) = ani*(t2 + t4)
         g(jj,2*k,j,l) = ani*(t3 + t5)
         g(jj,ny3-2*k,j,l) = ani*(t2 - t4)
         g(jj,ny3-2*k+1,j,l) = ani*(t5 - t3)
  140    continue
  150    continue
  160    continue
         do 180 j = kxpi, kxpt
         do 170 jj = 1, 2
         g(jj,nyh+1,j,l) = g(jj,nyh+1,j,l)
         g(jj,nyh+2,j,l) = -g(jj,nyh+2,j,l)
         t2 = g(jj,1,j,l) + g(jj,2,j,l)
         g(jj,2,j,l) = g(jj,1,j,l) - g(jj,2,j,l)
         g(jj,1,j,l) = t2
         g(jj,ny+1,j,l) = g(jj,ny+1,j,l)
  170    continue
  180    continue
  190    continue
c forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nyh
         kxpt = kxps
         do 250 l = 1, jblok
         if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
         do 220 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         do 210 j = kxpi, kxpt
         do 200 jj = 1, 2
         t4 = g(jj,ny3-2*k,j,l)
         t5 = -g(jj,ny3-2*k+1,j,l)
         t2 = g(jj,2*k-1,j,l) + t4
         t3 = g(jj,2*k,j,l) + t5
         t6 = g(jj,2*k-1,j,l) - t4
         t5 = g(jj,2*k,j,l) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         g(jj,2*k-1,j,l) = t2 + t4
         g(jj,2*k,j,l) = t3 + t5
         g(jj,ny3-2*k,j,l) = t2 - t4
         g(jj,ny3-2*k+1,j,l) = t5 - t3
  200    continue
  210    continue
  220    continue
         do 240 j = kxpi, kxpt
         do 230 jj = 1, 2
         g(jj,nyh+1,j,l) = 2.0*g(jj,nyh+1,j,l)
         g(jj,nyh+2,j,l) = -2.0*g(jj,nyh+2,j,l)
         t2 = 2.0*(g(jj,1,j,l) + g(jj,2,j,l))
         g(jj,2,j,l) = 2.0*(g(jj,1,j,l) - g(jj,2,j,l))
         g(jj,1,j,l) = t2
         g(jj,ny+1,j,l) = 2.0*g(jj,ny+1,j,l)
  230    continue
  240    continue
  250    continue
      endif
c perform recursion for sine-cosine transform
      kxpt = kxps
      do 280 l = 1, jblok
      if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
      do 270 j = kxpi, kxpt
      sum1 = .5*g(1,1,j,l)
      g(1,1,j,l) = 0.0
      g(1,2,j,l) = sum1
      sum2 = g(2,ny+1,j,l)
      g(2,ny+1,j,l) = g(2,2,j,l)
      g(2,2,j,l) = sum2
      do 260 k = 2, nyh
      sum1 = sum1 + g(1,2*k-1,j,l)
      g(1,2*k-1,j,l) = -g(1,2*k,j,l)
      g(1,2*k,j,l) = sum1
      sum2 = sum2 - g(2,2*k,j,l)
      g(2,2*k,j,l) = sum2
  260 continue
      g(1,ny+1,j,l) = 0.0
  270 continue
  280 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PFCST2R2Y(g,isign,mixup,sctd,indx,indy,kstrt,kxp,kxpi,k
     1xpp,nyv,kxpd,jblok,nxhyd,nxyd)
c this subroutine performs the y part of 2 two dimensional fast real
c sine and cosine transforms and their inverses, for a subset of x,
c using real arithmetic, for data which is distributed in blocks
c algorithm is described in Numerical Recipies in Fortran, Second Ed.,
c by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
c [Cambridge Univ. Press, 1992], p. 508.
c for isign = (-1,1), input: all, output: f
c approximate flop count: N*(5*log2(N) + 18)/nvp
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, inverse sine-cosine transform are performed
c g(1,m,n,i) = (.5*g(1,1,n,i) + ((-1)**m)*g(1,ny+1,n,i)
c              + sum(g(1,k,n,i)*cos(pi*m*k/ny))
c g(2,m,n,i) = sum(g(2,k,n,i)*sin(pi*m*k/ny))
c if isign = 1, a forward sine-cosine transforms are performed
c g(1,k,n,i) = 2*(.5*g(1,1,n,i) + ((-1)**m)*g(1,ny+1,n,i)
c              + sum(g(1,m,n,i)*cos(pi*m*k/ny))
c g(2,k,n,i) = sum(g(2,m,n,i)*sin(pi*m*k/ny))
c mixup = array of bit reversed addresses
c sctd = sine/cosine table
c kstrt = starting data block number
c kxp = number of data values per block in x
c kxpi = initial x index used
c kxpp = number of x indices used
c nyv = first dimension of g >= ny + 1
c kxpd = second dimension of g >= kxp + 1
c jblok = number of data blocks in x
c nxhyd = maximum of (nx/2,ny)
c nxyd = maximum of (nx,ny)
c written by viktor k. decyk, ucla
      implicit none
      integer isign, mixup, indx, indy, kstrt, kxp, kxpi, kxpp
      integer nyv, kxpd, jblok, nxhyd, nxyd
      real g
      complex sctd
      dimension g(2,nyv,kxpd,jblok)
      dimension mixup(nxhyd), sctd(nxyd)
c local data
      integer indx1, indy1, indx1y, nx, ny, nyh, nyhh, ny3, nxy, nxhy
      integer i, j, k, l, m, ks, km, kmr, nry, j1, j2, ns, ns2, k1, k2
      integer kxps, kxpt, jj
      real at1, at2, at3, t2, t3, t4, t5, t6, ani, sum1, sum2
      complex t1
      indx1 = indx - 1
      indy1 = indy - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      nyhh = ny/4
      ny3 = ny + 3
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 2
      kxps = kxpi + kxpp - 1
      if (kstrt.gt.nx) return
      if (isign.eq.0) return
c create auxiliary array in y
      kmr = nxy/ny
      kxpt = kxps
      do 30 l = 1, jblok
      if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
      do 20 j = kxpi, kxpt
      sum1 = .5*(g(1,1,j,l) - g(1,ny+1,j,l))
      do 10 k = 2, nyh
      k1 = 1 + kmr*(k - 1)
      at3 = -aimag(sctd(k1))
      at2 = g(1,ny+2-k,j,l)
      at1 = g(1,k,j,l) + at2
      at2 = g(1,k,j,l) - at2
      sum1 = sum1 + real(sctd(k1))*at2
      at2 = at3*at2
      at1 = .5*at1
      g(1,k,j,l) = at1 - at2
      g(1,ny+2-k,j,l) = at1 + at2
      at2 = g(2,ny+2-k,j,l)
      at1 = g(2,k,j,l) + at2
      at2 = g(2,k,j,l) - at2
      at1 = at3*at1
      at2 = .5*at2
      g(2,k,j,l) = at1 + at2
      g(2,ny+2-k,j,l) = at1 - at2
   10 continue
      g(1,1,j,l) = .5*(g(1,1,j,l) + g(1,ny+1,j,l))
      g(1,ny+1,j,l) = sum1
      g(2,1,j,l) = 0.0
      g(2,nyh+1,j,l) = 2.0*g(2,nyh+1,j,l)
   20 continue
   30 continue
c bit-reverse array elements in y
      nry = nxhy/nyh
      kxpt = kxps
      do 70 l = 1, jblok
      if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
      do 60 k = 1, nyh
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 60
      do 50 j = kxpi, kxpt
      do 40 jj = 1, 2
      t2 = g(jj,2*k1-1,j,l)
      t3 = g(jj,2*k1,j,l)
      g(jj,2*k1-1,j,l) = g(jj,2*k-1,j,l)
      g(jj,2*k1,j,l) = g(jj,2*k,j,l)
      g(jj,2*k-1,j,l) = t2
      g(jj,2*k,j,l) = t3
   40 continue
   50 continue
   60 continue
   70 continue
c first transform in y
      nry = nxy/nyh
      do 130 m = 1, indy1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyhh/ns
      kmr = 2*km*nry
      kxpt = kxps
      do 120 l = 1, jblok
      if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
      do 110 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 100 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 90 i = kxpi, kxpt
      do 80 jj = 1, 2
      t2 = real(t1)*g(jj,2*j2-1,i,l) - aimag(t1)*g(jj,2*j2,i,l)
      t3 = aimag(t1)*g(jj,2*j2-1,i,l) + real(t1)*g(jj,2*j2,i,l)
      g(jj,2*j2-1,i,l) = g(jj,2*j1-1,i,l) - t2
      g(jj,2*j2,i,l) = g(jj,2*j1,i,l) - t3
      g(jj,2*j1-1,i,l) = g(jj,2*j1-1,i,l) + t2
      g(jj,2*j1,i,l) = g(jj,2*j1,i,l) + t3
   80 continue
   90 continue
  100 continue
  110 continue
  120 continue
  130 continue
c unscramble coefficients and normalize
c inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nyh
         ani = 0.5
         kxpt = kxps
         do 190 l = 1, jblok
         if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
         do 160 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         do 150 j = kxpi, kxpt
         do 140 jj = 1, 2
         t4 = g(jj,ny3-2*k,j,l)
         t5 = -g(jj,ny3-2*k+1,j,l)
         t2 = g(jj,2*k-1,j,l) + t4
         t3 = g(jj,2*k,j,l) + t5
         t6 = g(jj,2*k-1,j,l) - t4
         t5 = g(jj,2*k,j,l) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         g(jj,2*k-1,j,l) = ani*(t2 + t4)
         g(jj,2*k,j,l) = ani*(t3 + t5)
         g(jj,ny3-2*k,j,l) = ani*(t2 - t4)
         g(jj,ny3-2*k+1,j,l) = ani*(t5 - t3)
  140    continue
  150    continue
  160    continue
         do 180 j = kxpi, kxpt
         do 170 jj = 1, 2
         g(jj,nyh+1,j,l) = g(jj,nyh+1,j,l)
         g(jj,nyh+2,j,l) = -g(jj,nyh+2,j,l)
         t2 = g(jj,1,j,l) + g(jj,2,j,l)
         g(jj,2,j,l) = g(jj,1,j,l) - g(jj,2,j,l)
         g(jj,1,j,l) = t2
         g(jj,ny+1,j,l) = g(jj,ny+1,j,l)
  170    continue
  180    continue
  190    continue
c forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nyh
         kxpt = kxps
         do 250 l = 1, jblok
         if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
         do 220 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         do 210 j = kxpi, kxpt
         do 200 jj = 1, 2
         t4 = g(jj,ny3-2*k,j,l)
         t5 = -g(jj,ny3-2*k+1,j,l)
         t2 = g(jj,2*k-1,j,l) + t4
         t3 = g(jj,2*k,j,l) + t5
         t6 = g(jj,2*k-1,j,l) - t4
         t5 = g(jj,2*k,j,l) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         g(jj,2*k-1,j,l) = t2 + t4
         g(jj,2*k,j,l) = t3 + t5
         g(jj,ny3-2*k,j,l) = t2 - t4
         g(jj,ny3-2*k+1,j,l) = t5 - t3
  200    continue
  210    continue
  220    continue
         do 240 j = kxpi, kxpt
         do 230 jj = 1, 2
         g(jj,nyh+1,j,l) = 2.0*g(jj,nyh+1,j,l)
         g(jj,nyh+2,j,l) = -2.0*g(jj,nyh+2,j,l)
         t2 = 2.0*(g(jj,1,j,l) + g(jj,2,j,l))
         g(jj,2,j,l) = 2.0*(g(jj,1,j,l) - g(jj,2,j,l))
         g(jj,1,j,l) = t2
         g(jj,ny+1,j,l) = 2.0*g(jj,ny+1,j,l)
  230    continue
  240    continue
  250    continue
      endif
c perform recursion for sine-cosine transform
      kxpt = kxps
      do 280 l = 1, jblok
      if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
      do 270 j = kxpi, kxpt
      sum1 = g(1,ny+1,j,l)
      g(1,ny+1,j,l) = g(1,2,j,l)
      g(1,2,j,l) = sum1
      sum2 = .5*g(2,1,j,l)
      g(2,1,j,l) = 0.0
      g(2,2,j,l) = sum2
      do 260 k = 2, nyh
      sum1 = sum1 - g(1,2*k,j,l)
      g(1,2*k,j,l) = sum1
      sum2 = sum2 + g(2,2*k-1,j,l)
      g(2,2*k-1,j,l) = -g(2,2*k,j,l)
      g(2,2*k,j,l) = sum2
  260 continue
      g(2,ny+1,j,l) = 0.0
  270 continue
  280 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine WPFCST2R3(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,in
     1dy,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhyd,nxyd)
c wrapper function for 3 parallel real sine/sine transforms
c for the electric field with dirichlet or magnetic field with neumann
c boundary conditions
      implicit none
      integer isign, ntpose, mixup, indx, indy, kstrt, nxvh, nyv
      integer kxp2, kyp, kypd, kxp2d, jblok, kblok, nxhyd, nxyd
      real f, g, bs, br, ttp
      complex sctd
      dimension f(3,2*nxvh,kypd,kblok), g(3,nyv,kxp2d,jblok)
      dimension bs(3,kxp2+1,kyp+1,kblok), br(3,kxp2+1,kyp+1,jblok)
      dimension mixup(nxhyd), sctd(nxyd)
c local data
      integer nx, ny, kxpi, kypi
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
c inverse fourier transform
      if (isign.lt.0) then
c perform x cosine-sine transform
         call PFCSST2R3X(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kyp
     1,nxvh,kypd,kblok,nxhyd,nxyd)
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PR3TPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,nyv,kxp2,kyp,kxp2d,k
     1ypd,jblok,kblok)
         call PWTIMERA(1,ttp,dtime)
c perform y sine-cosine transform
         call PFSCST2R3Y(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kx
     1p2,nyv,kxp2d,jblok,nxhyd,nxyd)
c transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PR3TPOSE(g,f,br,bs,ny,nx,kstrt,nyv,2*nxvh,kyp,kxp2,kypd
     1,kxp2d,kblok,jblok)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PR3TPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,nyv,kxp2,kyp,kxp2
     1d,kypd,jblok,kblok)
            call PWTIMERA(1,tf,dtime)
         endif
c perform y sine-cosine transform
         call PFSCST2R3Y(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kx
     1p2,nyv,kxp2d,jblok,nxhyd,nxyd)
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PR3TPOSE(g,f,br,bs,ny,nx,kstrt,nyv,2*nxvh,kyp,kxp2,kypd,kx
     1p2d,kblok,jblok)
         call PWTIMERA(1,ttp,dtime)
c perform x cosine-sine transform
         call PFCSST2R3X(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kyp
     1,nxvh,kypd,kblok,nxhyd,nxyd)
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine WPFSCT2R3(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,in
     1dy,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhyd,nxyd)
c wrapper function for 3 parallel real sine/sine transforms
c for the magnetic field with dirichlet or electric field with neumann
c boundary conditions
      implicit none
      integer isign, ntpose, mixup, indx, indy, kstrt, nxvh, nyv
      integer kxp2, kyp, kypd, kxp2d, jblok, kblok, nxhyd, nxyd
      real f, g, bs, br, ttp
      complex sctd
      dimension f(3,2*nxvh,kypd,kblok), g(3,nyv,kxp2d,jblok)
      dimension bs(3,kxp2+1,kyp+1,kblok), br(3,kxp2+1,kyp+1,jblok)
      dimension mixup(nxhyd), sctd(nxyd)
c local data
      integer nx, ny, kxpi, kypi
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
c inverse fourier transform
      if (isign.lt.0) then
c perform x sine-cosine transform
         call PFSCCT2R3X(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kyp
     1,nxvh,kypd,kblok,nxhyd,nxyd)
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PR3TPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,nyv,kxp2,kyp,kxp2d,k
     1ypd,jblok,kblok)
         call PWTIMERA(1,ttp,dtime)
c perform y cosine-sine transform
         call PFCSCT2R3Y(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kx
     1p2,nyv,kxp2d,jblok,nxhyd,nxyd)
c transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PR3TPOSE(g,f,br,bs,ny,nx,kstrt,nyv,2*nxvh,kyp,kxp2,kypd
     1,kxp2d,kblok,jblok)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PR3TPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,nyv,kxp2,kyp,kxp2
     1d,kypd,jblok,kblok)
            call PWTIMERA(1,tf,dtime)
         endif
c perform y cosine-sine transform
         call PFCSCT2R3Y(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kx
     1p2,nyv,kxp2d,jblok,nxhyd,nxyd)
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PR3TPOSE(g,f,br,bs,ny,nx,kstrt,nyv,2*nxvh,kyp,kxp2,kypd,kx
     1p2d,kblok,jblok)
         call PWTIMERA(1,ttp,dtime)
c perform x sine-cosine transform
         call PFSCCT2R3X(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kyp
     1,nxvh,kypd,kblok,nxhyd,nxyd)
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine PFCSST2R3X(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,
     1kypp,nxvh,kypd,kblok,nxhyd,nxyd)
c this subroutine performs the x part of 2 two dimensional fast real
c sine and cosine transforms and their inverses, for a subset of y,
c using real arithmetic, for data which is distributed in blocks
c algorithm is described in Numerical Recipies in Fortran, Second Ed.,
c by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
c [Cambridge Univ. Press, 1992], p. 508.
c for isign = (-1,1), input: all, output: f
c approximate flop count: N*(5*log2(N) + 18)/nvp
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, inverse sine-cosine transforms are performed
c f(1,n,k,i) = (1/nx*ny)*(.5*f(1,1,k,i) + ((-1)**n)*f(1,nx+1,k,i)
c              + sum(f(1,j,k,i)*cos(pi*n*j/nx)))
c f(2:3,n,k,i) = (1/nx*ny)*sum(f(2:3,j,k,i)*sin(pi*n*j/nx))
c if isign = 1, forward sine transforms are performed
c f(1,j,k,i) = 2*(.5*f(1,1,k,i) + ((-1)**j)*f(1,n+1,k,i)
c              + sum(f(1,n,k,i)*cos(pi*n*j/nx))
c f(2:3,j,k,i) = sum(f(2:3,n,k,i)*sin(pi*n*j/nx))
c mixup = array of bit reversed addresses
c sctd = sine/cosine table
c kstrt = starting data block number
c kyp = number of data values per block in y
c kypi = initial y index used
c kypp = number of y indices used
c nxvh = first dimension of f >= nx/2 + 1
c kypd = second dimension of f >= kyp + 1
c kblok = number of data blocks in y
c nxhyd = maximum of (nx/2,ny)
c nxyd = maximum of (nx,ny)
c written by viktor k. decyk, ucla
      implicit none
      integer isign, mixup, indx, indy, kstrt, kyp, kypi, kypp
      integer nxvh, kypd, kblok, nxhyd, nxyd
      real f
      complex sctd
      dimension f(3,2*nxvh,kypd,kblok)
      dimension mixup(nxhyd), sctd(nxyd)
c local data
      integer indx1, indx1y, nx, nxh, nxhh, nx3, ny, nxy, nxhy, ks, kypt
      integer i, j, k, l, m, km, kmr, nrx, j1, j2, ns, ns2, k1, k2, kyps
      integer jj
      real at1, at2, at3, t2, t3, t4, t5, t6, ani, sum1, sum2, sum3
      complex t1
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nx3 = nx + 3
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 2
      kyps = kypi + kypp - 1
      if (kstrt.gt.ny) return
      if (isign.eq.0) return
c create auxiliary array in x
      kmr = nxy/nx
      kypt = kyps
      do 30 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 20 k = kypi, kypt
      sum1 = .5*(f(1,1,k,l) - f(1,nx+1,k,l))
      do 10 j = 2, nxh
      j1 = 1 + kmr*(j - 1)
      at3 = -aimag(sctd(j1))
      at2 = f(1,nx+2-j,k,l)
      at1 = f(1,j,k,l) + at2
      at2 = f(1,j,k,l) - at2
      sum1 = sum1 + real(sctd(j1))*at2
      at2 = at3*at2
      at1 = .5*at1
      f(1,j,k,l) = at1 - at2
      f(1,nx+2-j,k,l) = at1 + at2
      at2 = f(2,nx+2-j,k,l)
      at1 = f(2,j,k,l) + at2
      at2 = f(2,j,k,l) - at2
      at1 = at3*at1
      at2 = .5*at2
      f(2,j,k,l) = at1 + at2
      f(2,nx+2-j,k,l) = at1 - at2
      at2 = f(3,nx+2-j,k,l)
      at1 = f(3,j,k,l) + at2
      at2 = f(3,j,k,l) - at2
      at1 = at3*at1
      at2 = .5*at2
      f(3,j,k,l) = at1 + at2
      f(3,nx+2-j,k,l) = at1 - at2
   10 continue
      f(1,1,k,l) = .5*(f(1,1,k,l) + f(1,nx+1,k,l))
      f(1,nx+1,k,l) = sum1
      f(2,1,k,l) = 0.0
      f(2,nxh+1,k,l) = 2.0*f(2,nxh+1,k,l)
      f(3,1,k,l) = 0.0
      f(3,nxh+1,k,l) = 2.0*f(3,nxh+1,k,l)
   20 continue
   30 continue
c bit-reverse array elements in x
      nrx = nxhy/nxh
      kypt = kyps
      do 70 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 60 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 60
      do 50 k = kypi, kypt
      do 40 jj = 1, 3
      t2 = f(jj,2*j1-1,k,l)
      t3 = f(jj,2*j1,k,l)
      f(jj,2*j1-1,k,l) = f(jj,2*j-1,k,l)
      f(jj,2*j1,k,l) = f(jj,2*j,k,l)
      f(jj,2*j-1,k,l) = t2
      f(jj,2*j,k,l) = t3
   40 continue
   50 continue
   60 continue
   70 continue
c first transform in x
      nrx = nxy/nxh
      do 130 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = 2*km*nrx
      kypt = kyps
      do 120 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 110 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 100 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 90 i = kypi, kypt
      do 80 jj = 1, 3
      t2 = real(t1)*f(jj,2*j2-1,i,l) - aimag(t1)*f(jj,2*j2,i,l)
      t3 = aimag(t1)*f(jj,2*j2-1,i,l) + real(t1)*f(jj,2*j2,i,l)
      f(jj,2*j2-1,i,l) = f(jj,2*j1-1,i,l) - t2
      f(jj,2*j2,i,l) = f(jj,2*j1,i,l) - t3
      f(jj,2*j1-1,i,l) = f(jj,2*j1-1,i,l) + t2
      f(jj,2*j1,i,l) = f(jj,2*j1,i,l) + t3
   80 continue
   90 continue
  100 continue
  110 continue
  120 continue
  130 continue
c unscramble coefficients and normalize
c inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nxh
         ani = 1./float(2*nx*ny)
         kypt = kyps
         do 190 l = 1, kblok
         if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
         do 160 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 150 k = kypi, kypt
         do 140 jj = 1, 3
         t4 = f(jj,nx3-2*j,k,l)
         t5 = -f(jj,nx3-2*j+1,k,l)
         t2 = f(jj,2*j-1,k,l) + t4
         t3 = f(jj,2*j,k,l) + t5
         t6 = f(jj,2*j-1,k,l) - t4
         t5 = f(jj,2*j,k,l) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,2*j-1,k,l) = ani*(t2 + t4)
         f(jj,2*j,k,l) = ani*(t3 + t5)
         f(jj,nx3-2*j,k,l) = ani*(t2 - t4)
         f(jj,nx3-2*j+1,k,l) = ani*(t5 - t3)
  140    continue
  150    continue
  160    continue
         ani = 2.*ani
         do 180 k = kypi, kypt
         do 170 jj = 1, 3
         f(jj,nxh+1,k,l) = ani*f(jj,nxh+1,k,l)
         f(jj,nxh+2,k,l) = -ani*f(jj,nxh+2,k,l)
         t2 = ani*(f(jj,1,k,l) + f(jj,2,k,l))
         f(jj,2,k,l) = ani*(f(jj,1,k,l) - f(jj,2,k,l))
         f(jj,1,k,l) = t2
         f(jj,nx+1,k,l) = ani*f(jj,nx+1,k,l)
  170    continue
  180    continue
  190    continue
c forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nxh
         kypt = kyps
         do 250 l = 1, kblok
         if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
         do 220 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 210 k = kypi, kypt
         do 200 jj = 1, 3
         t4 = f(jj,nx3-2*j,k,l)
         t5 = -f(jj,nx3-2*j+1,k,l)
         t2 = f(jj,2*j-1,k,l) + t4
         t3 = f(jj,2*j,k,l) + t5
         t6 = f(jj,2*j-1,k,l) - t4
         t5 = f(jj,2*j,k,l) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,2*j-1,k,l) = t2 + t4
         f(jj,2*j,k,l) = t3 + t5
         f(jj,nx3-2*j,k,l) = t2 - t4
         f(jj,nx3-2*j+1,k,l) = t5 - t3
  200    continue
  210    continue
  220    continue
         do 240 k = kypi, kypt
         do 230 jj = 1, 3
         f(jj,nxh+1,k,l) = 2.0*f(jj,nxh+1,k,l)
         f(jj,nxh+2,k,l) = -2.0*f(jj,nxh+2,k,l)
         t2 = 2.0*(f(jj,1,k,l) + f(jj,2,k,l))
         f(jj,2,k,l) = 2.0*(f(jj,1,k,l) - f(jj,2,k,l))
         f(jj,1,k,l) = t2
         f(jj,nx+1,k,l) = 2.0*f(jj,nx+1,k,l)
  230    continue
  240    continue
  250    continue
      endif
c perform recursion for cosine-sine transform
      kypt = kyps
      do 280 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 270 k = kypi, kypt
      sum1 = f(1,nx+1,k,l)
      f(1,nx+1,k,l) = f(1,2,k,l)
      f(1,2,k,l) = sum1
      sum2 = .5*f(2,1,k,l)
      f(2,1,k,l) = 0.0
      f(2,2,k,l) = sum2
      sum3 = .5*f(3,1,k,l)
      f(3,1,k,l) = 0.0
      f(3,2,k,l) = sum3
      do 260 j = 2, nxh
      sum1 = sum1 - f(1,2*j,k,l)
      f(1,2*j,k,l) = sum1
      sum2 = sum2 + f(2,2*j-1,k,l)
      f(2,2*j-1,k,l) = -f(2,2*j,k,l)
      f(2,2*j,k,l) = sum2
      sum3 = sum3 + f(3,2*j-1,k,l)
      f(3,2*j-1,k,l) = -f(3,2*j,k,l)
      f(3,2*j,k,l) = sum3
  260 continue
      f(2,nx+1,k,l) = 0.0
      f(3,nx+1,k,l) = 0.0
  270 continue
  280 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PFSCCT2R3X(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,
     1kypp,nxvh,kypd,kblok,nxhyd,nxyd)
c this subroutine performs the x part of 3 two dimensional fast real
c sine and cosine transforms and their inverses, for a subset of y,
c using real arithmetic, for data which is distributed in blocks
c algorithm is described in Numerical Recipies in Fortran, Second Ed.,
c by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
c [Cambridge Univ. Press, 1992], p. 508.
c for isign = (-1,1), input: all, output: f
c approximate flop count: N*(5*log2(N) + 18)/nvp
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, inverse sine-cosine transforms are performed
c f(1,n,k,i) = (1/nx*ny)*sum(f(1,j,k,i)*sin(pi*n*j/nx))
c f(2:3,n,k,i) = (1/nx*ny)*(.5*f(2:3,1,k,i) + ((-1)**n)*f(2:3,nx+1,k,i)
c              + sum(f(2:3,j,k,i)*cos(pi*n*j/nx)))
c if isign = 1, forward sine transforms are performed
c f(1,j,k,i) = sum(f(1,n,k,i)*sin(pi*n*j/nx))
c f(2:3,j,k,i) = 2*(.5*f(2:3,1,k,i) + ((-1)**j)*f(2:3,n+1,k,i)
c              + sum(f(2:3,n,k,i)*cos(pi*n*j/nx))
c mixup = array of bit reversed addresses
c sctd = sine/cosine table
c kstrt = starting data block number
c kyp = number of data values per block in y
c kypi = initial y index used
c kypp = number of y indices used
c nxvh = first dimension of f >= nx/2 + 1
c kypd = second dimension of f >= kyp + 1
c kblok = number of data blocks in y
c nxhyd = maximum of (nx/2,ny)
c nxyd = maximum of (nx,ny)
c written by viktor k. decyk, ucla
      implicit none
      integer isign, mixup, indx, indy, kstrt, kyp, kypi, kypp
      integer nxvh, kypd, kblok, nxhyd, nxyd
      real f
      complex sctd
      dimension f(3,2*nxvh,kypd,kblok)
      dimension mixup(nxhyd), sctd(nxyd)
c local data
      integer indx1, indx1y, nx, nxh, nxhh, nx3, ny, nxy, nxhy, ks, kypt
      integer i, j, k, l, m, km, kmr, nrx, j1, j2, ns, ns2, k1, k2, kyps
      integer jj
      real at1, at2, at3, t2, t3, t4, t5, t6, ani, sum1, sum2, sum3
      complex t1
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nx3 = nx + 3
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 2
      kyps = kypi + kypp - 1
      if (kstrt.gt.ny) return
      if (isign.eq.0) return
c create auxiliary array in x
      kmr = nxy/nx
      kypt = kyps
      do 30 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 20 k = kypi, kypt
      sum1 = .5*(f(2,1,k,l) - f(2,nx+1,k,l))
      sum2 = .5*(f(3,1,k,l) - f(3,nx+1,k,l))
      do 10 j = 2, nxh
      j1 = 1 + kmr*(j - 1)
      at3 = -aimag(sctd(j1))
      at2 = f(1,nx+2-j,k,l)
      at1 = f(1,j,k,l) + at2
      at2 = f(1,j,k,l) - at2
      at1 = at3*at1
      at2 = .5*at2
      f(1,j,k,l) = at1 + at2
      f(1,nx+2-j,k,l) = at1 - at2
      at2 = f(2,nx+2-j,k,l)
      at1 = f(2,j,k,l) + at2
      at2 = f(2,j,k,l) - at2
      sum1 = sum1 + real(sctd(j1))*at2
      at2 = at3*at2
      at1 = .5*at1
      f(2,j,k,l) = at1 - at2
      f(2,nx+2-j,k,l) = at1 + at2
      at2 = f(3,nx+2-j,k,l)
      at1 = f(3,j,k,l) + at2
      at2 = f(3,j,k,l) - at2
      sum2 = sum2 + real(sctd(j1))*at2
      at2 = at3*at2
      at1 = .5*at1
      f(3,j,k,l) = at1 - at2
      f(3,nx+2-j,k,l) = at1 + at2
   10 continue
      f(1,1,k,l) = 0.0
      f(1,nxh+1,k,l) = 2.0*f(1,nxh+1,k,l)
      f(2,1,k,l) = .5*(f(2,1,k,l) + f(2,nx+1,k,l))
      f(2,nx+1,k,l) = sum1
      f(3,1,k,l) = .5*(f(3,1,k,l) + f(3,nx+1,k,l))
      f(3,nx+1,k,l) = sum2
   20 continue
   30 continue
c bit-reverse array elements in x
      nrx = nxhy/nxh
      kypt = kyps
      do 70 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 60 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 60
      do 50 k = kypi, kypt
      do 40 jj = 1, 3
      t2 = f(jj,2*j1-1,k,l)
      t3 = f(jj,2*j1,k,l)
      f(jj,2*j1-1,k,l) = f(jj,2*j-1,k,l)
      f(jj,2*j1,k,l) = f(jj,2*j,k,l)
      f(jj,2*j-1,k,l) = t2
      f(jj,2*j,k,l) = t3
   40 continue
   50 continue
   60 continue
   70 continue
c first transform in x
      nrx = nxy/nxh
      do 130 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = 2*km*nrx
      kypt = kyps
      do 120 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 110 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 100 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 90 i = kypi, kypt
      do 80 jj = 1, 3
      t2 = real(t1)*f(jj,2*j2-1,i,l) - aimag(t1)*f(jj,2*j2,i,l)
      t3 = aimag(t1)*f(jj,2*j2-1,i,l) + real(t1)*f(jj,2*j2,i,l)
      f(jj,2*j2-1,i,l) = f(jj,2*j1-1,i,l) - t2
      f(jj,2*j2,i,l) = f(jj,2*j1,i,l) - t3
      f(jj,2*j1-1,i,l) = f(jj,2*j1-1,i,l) + t2
      f(jj,2*j1,i,l) = f(jj,2*j1,i,l) + t3
   80 continue
   90 continue
  100 continue
  110 continue
  120 continue
  130 continue
c unscramble coefficients and normalize
c inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nxh
         ani = 1./float(2*nx*ny)
         kypt = kyps
         do 190 l = 1, kblok
         if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
         do 160 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 150 k = kypi, kypt
         do 140 jj = 1, 3
         t4 = f(jj,nx3-2*j,k,l)
         t5 = -f(jj,nx3-2*j+1,k,l)
         t2 = f(jj,2*j-1,k,l) + t4
         t3 = f(jj,2*j,k,l) + t5
         t6 = f(jj,2*j-1,k,l) - t4
         t5 = f(jj,2*j,k,l) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,2*j-1,k,l) = ani*(t2 + t4)
         f(jj,2*j,k,l) = ani*(t3 + t5)
         f(jj,nx3-2*j,k,l) = ani*(t2 - t4)
         f(jj,nx3-2*j+1,k,l) = ani*(t5 - t3)
  140    continue
  150    continue
  160    continue
         ani = 2.*ani
         do 180 k = kypi, kypt
         do 170 jj = 1, 3
         f(jj,nxh+1,k,l) = ani*f(jj,nxh+1,k,l)
         f(jj,nxh+2,k,l) = -ani*f(jj,nxh+2,k,l)
         t2 = ani*(f(jj,1,k,l) + f(jj,2,k,l))
         f(jj,2,k,l) = ani*(f(jj,1,k,l) - f(jj,2,k,l))
         f(jj,1,k,l) = t2
         f(jj,nx+1,k,l) = ani*f(jj,nx+1,k,l)
  170    continue
  180    continue
  190    continue
c forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nxh
         kypt = kyps
         do 250 l = 1, kblok
         if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
         do 220 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 210 k = kypi, kypt
         do 200 jj = 1, 3
         t4 = f(jj,nx3-2*j,k,l)
         t5 = -f(jj,nx3-2*j+1,k,l)
         t2 = f(jj,2*j-1,k,l) + t4
         t3 = f(jj,2*j,k,l) + t5
         t6 = f(jj,2*j-1,k,l) - t4
         t5 = f(jj,2*j,k,l) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,2*j-1,k,l) = t2 + t4
         f(jj,2*j,k,l) = t3 + t5
         f(jj,nx3-2*j,k,l) = t2 - t4
         f(jj,nx3-2*j+1,k,l) = t5 - t3
  200    continue
  210    continue
  220    continue
         do 240 k = kypi, kypt
         do 230 jj = 1, 3
         f(jj,nxh+1,k,l) = 2.0*f(jj,nxh+1,k,l)
         f(jj,nxh+2,k,l) = -2.0*f(jj,nxh+2,k,l)
         t2 = 2.0*(f(jj,1,k,l) + f(jj,2,k,l))
         f(jj,2,k,l) = 2.0*(f(jj,1,k,l) - f(jj,2,k,l))
         f(jj,1,k,l) = t2
         f(jj,nx+1,k,l) = 2.0*f(jj,nx+1,k,l)
  230    continue
  240    continue
  250    continue
      endif
c perform recursion for cosine-sine transform
      kypt = kyps
      do 280 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 270 k = kypi, kypt
      sum1 = .5*f(1,1,k,l)
      f(1,1,k,l) = 0.0
      f(1,2,k,l) = sum1
      sum2 = f(2,nx+1,k,l)
      f(2,nx+1,k,l) = f(2,2,k,l)
      f(2,2,k,l) = sum2
      sum3 = f(3,nx+1,k,l)
      f(3,nx+1,k,l) = f(3,2,k,l)
      f(3,2,k,l) = sum3
      do 260 j = 2, nxh
      sum1 = sum1 + f(1,2*j-1,k,l)
      f(1,2*j-1,k,l) = -f(1,2*j,k,l)
      f(1,2*j,k,l) = sum1
      sum2 = sum2 - f(2,2*j,k,l)
      f(2,2*j,k,l) = sum2
      sum3 = sum3 - f(3,2*j,k,l)
      f(3,2*j,k,l) = sum3
  260 continue
      f(1,nx+1,k,l) = 0.0
  270 continue
  280 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PFSCST2R3Y(g,isign,mixup,sctd,indx,indy,kstrt,kxp,kxpi,
     1kxpp,nyv,kxpd,jblok,nxhyd,nxyd)
c this subroutine performs the y part of 3 two dimensional fast real
c sine and cosine transforms and their inverses, for a subset of x,
c using real arithmetic, for data which is distributed in blocks
c algorithm is described in Numerical Recipies in Fortran, Second Ed.,
c by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
c [Cambridge Univ. Press, 1992], p. 508.
c for isign = (-1,1), input: all, output: f
c approximate flop count: N*(5*log2(N) + 18)/nvp
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, inverse sine-cosine transform are performed
c g(1,m,n,i) = sum(g(1,k,n,i)*sin(pi*m*k/ny))
c g(2,m,n,i) = (.5*g(2,1,n,i) + ((-1)**m)*g(2,ny+1,n,i)
c              + sum(g(2,k,n,i)*cos(pi*m*k/ny))
c g(3,m,n,i) = sum(g(3,k,n,i)*sin(pi*m*k/ny))
c if isign = 1, a forward sine-cosine transforms are performed
c g(1,k,n,i) = sum(g(1,m,n,i)*sin(pi*m*k/ny))
c g(2,k,n,i) = 2*(.5*g(2,1,n,i) + ((-1)**m)*g(2,ny+1,n,i)
c              + sum(g(2,m,n,i)*cos(pi*m*k/ny))
c g(3,k,n,i) = sum(g(3,m,n,i)*sin(pi*m*k/ny))
c mixup = array of bit reversed addresses
c sctd = sine/cosine table
c kstrt = starting data block number
c kxp = number of data values per block in x
c kxpi = initial x index used
c kxpp = number of x indices used
c nyv = first dimension of g >= ny + 1
c kxpd = second dimension of g >= kxp + 1
c jblok = number of data blocks in x
c nxhyd = maximum of (nx/2,ny)
c nxyd = maximum of (nx,ny)
c written by viktor k. decyk, ucla
      implicit none
      integer isign, mixup, indx, indy, kstrt, kxp, kxpi, kxpp
      integer nyv, kxpd, jblok, nxhyd, nxyd
      real g
      complex sctd
      dimension g(3,nyv,kxpd,jblok)
      dimension mixup(nxhyd), sctd(nxyd)
c local data
      integer indx1, indy1, indx1y, nx, ny, nyh, nyhh, ny3, nxy, nxhy
      integer i, j, k, l, m, ks, km, kmr, nry, j1, j2, ns, ns2, k1, k2
      integer kxps, kxpt, jj
      real at1, at2, at3, t2, t3, t4, t5, t6, ani, sum1, sum2, sum3
      complex t1
      indx1 = indx - 1
      indy1 = indy - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      nyhh = ny/4
      ny3 = ny + 3
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 2
      kxps = kxpi + kxpp - 1
      if (kstrt.gt.nx) return
      if (isign.eq.0) return
c create auxiliary array in y
      kmr = nxy/ny
      kxpt = kxps
      do 30 l = 1, jblok
      if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
      do 20 j = kxpi, kxpt
      sum1 = .5*(g(2,1,j,l) - g(2,ny+1,j,l))
      do 10 k = 2, nyh
      k1 = 1 + kmr*(k - 1)
      at3 = -aimag(sctd(k1))
      at2 = g(1,ny+2-k,j,l)
      at1 = g(1,k,j,l) + at2
      at2 = g(1,k,j,l) - at2
      at1 = at3*at1
      at2 = .5*at2
      g(1,k,j,l) = at1 + at2
      g(1,ny+2-k,j,l) = at1 - at2
      at2 = g(2,ny+2-k,j,l)
      at1 = g(2,k,j,l) + at2
      at2 = g(2,k,j,l) - at2
      sum1 = sum1 + real(sctd(k1))*at2
      at2 = at3*at2
      at1 = .5*at1
      g(2,k,j,l) = at1 - at2
      g(2,ny+2-k,j,l) = at1 + at2
      at2 = g(3,ny+2-k,j,l)
      at1 = g(3,k,j,l) + at2
      at2 = g(3,k,j,l) - at2
      at1 = at3*at1
      at2 = .5*at2
      g(3,k,j,l) = at1 + at2
      g(3,ny+2-k,j,l) = at1 - at2
   10 continue
      g(1,1,j,l) = 0.0
      g(1,nyh+1,j,l) = 2.0*g(1,nyh+1,j,l)
      g(2,1,j,l) = .5*(g(2,1,j,l) + g(2,ny+1,j,l))
      g(2,ny+1,j,l) = sum1
      g(3,1,j,l) = 0.0
      g(3,nyh+1,j,l) = 2.0*g(3,nyh+1,j,l)
   20 continue
   30 continue
c bit-reverse array elements in y
      nry = nxhy/nyh
      kxpt = kxps
      do 70 l = 1, jblok
      if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
      do 60 k = 1, nyh
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 60
      do 50 j = kxpi, kxpt
      do 40 jj = 1, 3
      t2 = g(jj,2*k1-1,j,l)
      t3 = g(jj,2*k1,j,l)
      g(jj,2*k1-1,j,l) = g(jj,2*k-1,j,l)
      g(jj,2*k1,j,l) = g(jj,2*k,j,l)
      g(jj,2*k-1,j,l) = t2
      g(jj,2*k,j,l) = t3
   40 continue
   50 continue
   60 continue
   70 continue
c first transform in y
      nry = nxy/nyh
      do 130 m = 1, indy1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyhh/ns
      kmr = 2*km*nry
      kxpt = kxps
      do 120 l = 1, jblok
      if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
      do 110 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 100 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 90 i = kxpi, kxpt
      do 80 jj = 1, 3
      t2 = real(t1)*g(jj,2*j2-1,i,l) - aimag(t1)*g(jj,2*j2,i,l)
      t3 = aimag(t1)*g(jj,2*j2-1,i,l) + real(t1)*g(jj,2*j2,i,l)
      g(jj,2*j2-1,i,l) = g(jj,2*j1-1,i,l) - t2
      g(jj,2*j2,i,l) = g(jj,2*j1,i,l) - t3
      g(jj,2*j1-1,i,l) = g(jj,2*j1-1,i,l) + t2
      g(jj,2*j1,i,l) = g(jj,2*j1,i,l) + t3
   80 continue
   90 continue
  100 continue
  110 continue
  120 continue
  130 continue
c unscramble coefficients and normalize
c inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nyh
         ani = 0.5
         kxpt = kxps
         do 190 l = 1, jblok
         if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
         do 160 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         do 150 j = kxpi, kxpt
         do 140 jj = 1, 3
         t4 = g(jj,ny3-2*k,j,l)
         t5 = -g(jj,ny3-2*k+1,j,l)
         t2 = g(jj,2*k-1,j,l) + t4
         t3 = g(jj,2*k,j,l) + t5
         t6 = g(jj,2*k-1,j,l) - t4
         t5 = g(jj,2*k,j,l) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         g(jj,2*k-1,j,l) = ani*(t2 + t4)
         g(jj,2*k,j,l) = ani*(t3 + t5)
         g(jj,ny3-2*k,j,l) = ani*(t2 - t4)
         g(jj,ny3-2*k+1,j,l) = ani*(t5 - t3)
  140    continue
  150    continue
  160    continue
         do 180 j = kxpi, kxpt
         do 170 jj = 1, 3
         g(jj,nyh+1,j,l) = g(jj,nyh+1,j,l)
         g(jj,nyh+2,j,l) = -g(jj,nyh+2,j,l)
         t2 = g(jj,1,j,l) + g(jj,2,j,l)
         g(jj,2,j,l) = g(jj,1,j,l) - g(jj,2,j,l)
         g(jj,1,j,l) = t2
         g(jj,ny+1,j,l) = g(jj,ny+1,j,l)
  170    continue
  180    continue
  190    continue
c forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nyh
         kxpt = kxps
         do 250 l = 1, jblok
         if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
         do 220 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         do 210 j = kxpi, kxpt
         do 200 jj = 1, 3
         t4 = g(jj,ny3-2*k,j,l)
         t5 = -g(jj,ny3-2*k+1,j,l)
         t2 = g(jj,2*k-1,j,l) + t4
         t3 = g(jj,2*k,j,l) + t5
         t6 = g(jj,2*k-1,j,l) - t4
         t5 = g(jj,2*k,j,l) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         g(jj,2*k-1,j,l) = t2 + t4
         g(jj,2*k,j,l) = t3 + t5
         g(jj,ny3-2*k,j,l) = t2 - t4
         g(jj,ny3-2*k+1,j,l) = t5 - t3
  200    continue
  210    continue
  220    continue
         do 240 j = kxpi, kxpt
         do 230 jj = 1, 3
         g(jj,nyh+1,j,l) = 2.0*g(jj,nyh+1,j,l)
         g(jj,nyh+2,j,l) = -2.0*g(jj,nyh+2,j,l)
         t2 = 2.0*(g(jj,1,j,l) + g(jj,2,j,l))
         g(jj,2,j,l) = 2.0*(g(jj,1,j,l) - g(jj,2,j,l))
         g(jj,1,j,l) = t2
         g(jj,ny+1,j,l) = 2.0*g(jj,ny+1,j,l)
  230    continue
  240    continue
  250    continue
      endif
c perform recursion for sine-cosine transform
      kxpt = kxps
      do 280 l = 1, jblok
      if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
      do 270 j = kxpi, kxpt
      sum1 = .5*g(1,1,j,l)
      g(1,1,j,l) = 0.0
      g(1,2,j,l) = sum1
      sum2 = g(2,ny+1,j,l)
      g(2,ny+1,j,l) = g(2,2,j,l)
      g(2,2,j,l) = sum2
      sum3 = .5*g(3,1,j,l)
      g(3,1,j,l) = 0.0
      g(3,2,j,l) = sum3
      do 260 k = 2, nyh
      sum1 = sum1 + g(1,2*k-1,j,l)
      g(1,2*k-1,j,l) = -g(1,2*k,j,l)
      g(1,2*k,j,l) = sum1
      sum2 = sum2 - g(2,2*k,j,l)
      g(2,2*k,j,l) = sum2
      sum3 = sum3 + g(3,2*k-1,j,l)
      g(3,2*k-1,j,l) = -g(3,2*k,j,l)
      g(3,2*k,j,l) = sum3
  260 continue
      g(1,ny+1,j,l) = 0.0
      g(3,ny+1,j,l) = 0.0
  270 continue
  280 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PFCSCT2R3Y(g,isign,mixup,sctd,indx,indy,kstrt,kxp,kxpi,
     1kxpp,nyv,kxpd,jblok,nxhyd,nxyd)
c this subroutine performs the y part of 3 two dimensional fast real
c sine and cosine transforms and their inverses, for a subset of x,
c using real arithmetic, for data which is distributed in blocks
c algorithm is described in Numerical Recipies in Fortran, Second Ed.,
c by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
c [Cambridge Univ. Press, 1992], p. 508.
c for isign = (-1,1), input: all, output: f
c approximate flop count: N*(5*log2(N) + 18)/nvp
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, inverse sine-cosine transform are performed
c g(1,m,n,i) = (.5*g(1,1,n,i) + ((-1)**m)*g(1,ny+1,n,i)
c              + sum(g(1,k,n,i)*cos(pi*m*k/ny))
c g(2,m,n,i) = sum(g(2,k,n,i)*sin(pi*m*k/ny))
c g(3,m,n,i) = (.5*g(3,1,n,i) + ((-1)**m)*g(3,ny+1,n,i)
c              + sum(g(3,k,n,i)*cos(pi*m*k/ny))
c if isign = 1, a forward sine-cosine transforms are performed
c g(1,k,n,i) = 2*(.5*g(1,1,n,i) + ((-1)**m)*g(1,ny+1,n,i)
c              + sum(g(1,m,n,i)*cos(pi*m*k/ny))
c g(2,k,n,i) = sum(g(2,m,n,i)*sin(pi*m*k/ny))
c g(3,k,n,i) = 2*(.5*g(3,1,n,i) + ((-1)**m)*g(3,ny+1,n,i)
c              + sum(g(3,m,n,i)*cos(pi*m*k/ny))
c mixup = array of bit reversed addresses
c sctd = sine/cosine table
c kstrt = starting data block number
c kxp = number of data values per block in x
c kxpi = initial x index used
c kxpp = number of x indices used
c nyv = first dimension of g >= ny + 1
c kxpd = second dimension of g >= kxp + 1
c jblok = number of data blocks in x
c nxhyd = maximum of (nx/2,ny)
c nxyd = maximum of (nx,ny)
c written by viktor k. decyk, ucla
      implicit none
      integer isign, mixup, indx, indy, kstrt, kxp, kxpi, kxpp
      integer nyv, kxpd, jblok, nxhyd, nxyd
      real g
      complex sctd
      dimension g(3,nyv,kxpd,jblok)
      dimension mixup(nxhyd), sctd(nxyd)
c local data
      integer indx1, indy1, indx1y, nx, ny, nyh, nyhh, ny3, nxy, nxhy
      integer i, j, k, l, m, ks, km, kmr, nry, j1, j2, ns, ns2, k1, k2
      integer kxps, kxpt, jj
      real at1, at2, at3, t2, t3, t4, t5, t6, ani, sum1, sum2, sum3
      complex t1
      indx1 = indx - 1
      indy1 = indy - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      nyhh = ny/4
      ny3 = ny + 3
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 2
      kxps = kxpi + kxpp - 1
      if (kstrt.gt.nx) return
      if (isign.eq.0) return
c create auxiliary array in y
      kmr = nxy/ny
      kxpt = kxps
      do 30 l = 1, jblok
      if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
      do 20 j = kxpi, kxpt
      sum1 = .5*(g(1,1,j,l) - g(1,ny+1,j,l))
      sum2 = .5*(g(3,1,j,l) - g(3,ny+1,j,l))
      do 10 k = 2, nyh
      k1 = 1 + kmr*(k - 1)
      at3 = -aimag(sctd(k1))
      at2 = g(1,ny+2-k,j,l)
      at1 = g(1,k,j,l) + at2
      at2 = g(1,k,j,l) - at2
      sum1 = sum1 + real(sctd(k1))*at2
      at2 = at3*at2
      at1 = .5*at1
      g(1,k,j,l) = at1 - at2
      g(1,ny+2-k,j,l) = at1 + at2
      at2 = g(2,ny+2-k,j,l)
      at1 = g(2,k,j,l) + at2
      at2 = g(2,k,j,l) - at2
      at1 = at3*at1
      at2 = .5*at2
      g(2,k,j,l) = at1 + at2
      g(2,ny+2-k,j,l) = at1 - at2
      at2 = g(3,ny+2-k,j,l)
      at1 = g(3,k,j,l) + at2
      at2 = g(3,k,j,l) - at2
      sum2 = sum2 + real(sctd(k1))*at2
      at2 = at3*at2
      at1 = .5*at1
      g(3,k,j,l) = at1 - at2
      g(3,ny+2-k,j,l) = at1 + at2
   10 continue
      g(1,1,j,l) = .5*(g(1,1,j,l) + g(1,ny+1,j,l))
      g(1,ny+1,j,l) = sum1
      g(2,1,j,l) = 0.0
      g(2,nyh+1,j,l) = 2.0*g(2,nyh+1,j,l)
      g(3,1,j,l) = .5*(g(3,1,j,l) + g(3,ny+1,j,l))
      g(3,ny+1,j,l) = sum2
   20 continue
   30 continue
c bit-reverse array elements in y
      nry = nxhy/nyh
      kxpt = kxps
      do 70 l = 1, jblok
      if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
      do 60 k = 1, nyh
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 60
      do 50 j = kxpi, kxpt
      do 40 jj = 1, 3
      t2 = g(jj,2*k1-1,j,l)
      t3 = g(jj,2*k1,j,l)
      g(jj,2*k1-1,j,l) = g(jj,2*k-1,j,l)
      g(jj,2*k1,j,l) = g(jj,2*k,j,l)
      g(jj,2*k-1,j,l) = t2
      g(jj,2*k,j,l) = t3
   40 continue
   50 continue
   60 continue
   70 continue
c first transform in y
      nry = nxy/nyh
      do 130 m = 1, indy1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyhh/ns
      kmr = 2*km*nry
      kxpt = kxps
      do 120 l = 1, jblok
      if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
      do 110 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 100 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 90 i = kxpi, kxpt
      do 80 jj = 1, 3
      t2 = real(t1)*g(jj,2*j2-1,i,l) - aimag(t1)*g(jj,2*j2,i,l)
      t3 = aimag(t1)*g(jj,2*j2-1,i,l) + real(t1)*g(jj,2*j2,i,l)
      g(jj,2*j2-1,i,l) = g(jj,2*j1-1,i,l) - t2
      g(jj,2*j2,i,l) = g(jj,2*j1,i,l) - t3
      g(jj,2*j1-1,i,l) = g(jj,2*j1-1,i,l) + t2
      g(jj,2*j1,i,l) = g(jj,2*j1,i,l) + t3
   80 continue
   90 continue
  100 continue
  110 continue
  120 continue
  130 continue
c unscramble coefficients and normalize
c inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nyh
         ani = 0.5
         kxpt = kxps
         do 190 l = 1, jblok
         if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
         do 160 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         do 150 j = kxpi, kxpt
         do 140 jj = 1, 3
         t4 = g(jj,ny3-2*k,j,l)
         t5 = -g(jj,ny3-2*k+1,j,l)
         t2 = g(jj,2*k-1,j,l) + t4
         t3 = g(jj,2*k,j,l) + t5
         t6 = g(jj,2*k-1,j,l) - t4
         t5 = g(jj,2*k,j,l) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         g(jj,2*k-1,j,l) = ani*(t2 + t4)
         g(jj,2*k,j,l) = ani*(t3 + t5)
         g(jj,ny3-2*k,j,l) = ani*(t2 - t4)
         g(jj,ny3-2*k+1,j,l) = ani*(t5 - t3)
  140    continue
  150    continue
  160    continue
         do 180 j = kxpi, kxpt
         do 170 jj = 1, 3
         g(jj,nyh+1,j,l) = g(jj,nyh+1,j,l)
         g(jj,nyh+2,j,l) = -g(jj,nyh+2,j,l)
         t2 = g(jj,1,j,l) + g(jj,2,j,l)
         g(jj,2,j,l) = g(jj,1,j,l) - g(jj,2,j,l)
         g(jj,1,j,l) = t2
         g(jj,ny+1,j,l) = g(jj,ny+1,j,l)
  170    continue
  180    continue
  190    continue
c forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nyh
         kxpt = kxps
         do 250 l = 1, jblok
         if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
         do 220 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         do 210 j = kxpi, kxpt
         do 200 jj = 1, 3
         t4 = g(jj,ny3-2*k,j,l)
         t5 = -g(jj,ny3-2*k+1,j,l)
         t2 = g(jj,2*k-1,j,l) + t4
         t3 = g(jj,2*k,j,l) + t5
         t6 = g(jj,2*k-1,j,l) - t4
         t5 = g(jj,2*k,j,l) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         g(jj,2*k-1,j,l) = t2 + t4
         g(jj,2*k,j,l) = t3 + t5
         g(jj,ny3-2*k,j,l) = t2 - t4
         g(jj,ny3-2*k+1,j,l) = t5 - t3
  200    continue
  210    continue
  220    continue
         do 240 j = kxpi, kxpt
         do 230 jj = 1, 3
         g(jj,nyh+1,j,l) = 2.0*g(jj,nyh+1,j,l)
         g(jj,nyh+2,j,l) = -2.0*g(jj,nyh+2,j,l)
         t2 = 2.0*(g(jj,1,j,l) + g(jj,2,j,l))
         g(jj,2,j,l) = 2.0*(g(jj,1,j,l) - g(jj,2,j,l))
         g(jj,1,j,l) = t2
         g(jj,ny+1,j,l) = 2.0*g(jj,ny+1,j,l)
  230    continue
  240    continue
  250    continue
      endif
c perform recursion for sine-cosine transform
      kxpt = kxps
      do 280 l = 1, jblok
      if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
      do 270 j = kxpi, kxpt
      sum1 = g(1,ny+1,j,l)
      g(1,ny+1,j,l) = g(1,2,j,l)
      g(1,2,j,l) = sum1
      sum2 = .5*g(2,1,j,l)
      g(2,1,j,l) = 0.0
      g(2,2,j,l) = sum2
      sum3 = g(3,ny+1,j,l)
      g(3,ny+1,j,l) = g(3,2,j,l)
      g(3,2,j,l) = sum3
      do 260 k = 2, nyh
      sum1 = sum1 - g(1,2*k,j,l)
      g(1,2*k,j,l) = sum1
      sum2 = sum2 + g(2,2*k-1,j,l)
      g(2,2*k-1,j,l) = -g(2,2*k,j,l)
      g(2,2*k,j,l) = sum2
      sum3 = sum3 - g(3,2*k,j,l)
      g(3,2*k,j,l) = sum3
  260 continue
      g(2,ny+1,j,l) = 0.0
  270 continue
  280 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PFSSCT2R3Y(g,isign,mixup,sctd,indx,indy,kstrt,kxp,kxpi,
     1kxpp,nyv,kxpd,jblok,nxhyd,nxyd)
c this subroutine performs the y part of 3 two dimensional fast real
c sine and cosine transforms and their inverses, for a subset of x,
c using real arithmetic, for data which is distributed in blocks
c algorithm is described in Numerical Recipies in Fortran, Second Ed.,
c by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
c [Cambridge Univ. Press, 1992], p. 508.
c for isign = (-1,1), input: all, output: f
c approximate flop count: N*(5*log2(N) + 18)/nvp
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, inverse sine-cosine transform are performed
c g(1,m,n,i) = (.5*g(1,1,n,i) + ((-1)**m)*g(1,ny+1,n,i)
c              + sum(g(1,k,n,i)*cos(pi*m*k/ny))
c g(2,m,n,i) = sum(g(2,k,n,i)*sin(pi*m*k/ny))
c g(3,m,n,i) = (.5*g(3,1,n,i) + ((-1)**m)*g(3,ny+1,n,i)
c              + sum(g(3,k,n,i)*cos(pi*m*k/ny))
c if isign = 1, a forward sine-cosine transforms are performed
c g(1,k,n,i) = 2*(.5*g(1,1,n,i) + ((-1)**m)*g(1,ny+1,n,i)
c              + sum(g(1,m,n,i)*cos(pi*m*k/ny))
c g(2,k,n,i) = sum(g(2,m,n,i)*sin(pi*m*k/ny))
c g(3,k,n,i) = 2*(.5*g(3,1,n,i) + ((-1)**m)*g(3,ny+1,n,i)
c              + sum(g(3,m,n,i)*cos(pi*m*k/ny))
c mixup = array of bit reversed addresses
c sctd = sine/cosine table
c kstrt = starting data block number
c kxp = number of data values per block in x
c kxpi = initial x index used
c kxpp = number of x indices used
c nyv = first dimension of g >= ny + 1
c kxpd = second dimension of g >= kxp + 1
c jblok = number of data blocks in x
c nxhyd = maximum of (nx/2,ny)
c nxyd = maximum of (nx,ny)
c written by viktor k. decyk, ucla
      implicit none
      integer isign, mixup, indx, indy, kstrt, kxp, kxpi, kxpp
      integer nyv, kxpd, jblok, nxhyd, nxyd
      real g
      complex sctd
      dimension g(3,nyv,kxpd,jblok)
      dimension mixup(nxhyd), sctd(nxyd)
c local data
      integer indx1, indy1, indx1y, nx, ny, nyh, nyhh, ny3, nxy, nxhy
      integer i, j, k, l, m, ks, km, kmr, nry, j1, j2, ns, ns2, k1, k2
      integer kxps, kxpt, jj
      real at1, at2, at3, t2, t3, t4, t5, t6, ani, sum1, sum2, sum3
      complex t1
      indx1 = indx - 1
      indy1 = indy - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      nyhh = ny/4
      ny3 = ny + 3
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 2
      kxps = kxpi + kxpp - 1
      if (kstrt.gt.nx) return
      if (isign.eq.0) return
c create auxiliary array in y
      kmr = nxy/ny
      kxpt = kxps
      do 30 l = 1, jblok
      if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
      do 20 j = kxpi, kxpt
      sum2 = .5*(g(3,1,j,l) - g(3,ny+1,j,l))
      do 10 k = 2, nyh
      k1 = 1 + kmr*(k - 1)
      at3 = -aimag(sctd(k1))
      at2 = g(1,ny+2-k,j,l)
      at1 = g(1,k,j,l) + at2
      at2 = g(1,k,j,l) - at2
      at1 = at3*at1
      at2 = .5*at2
      g(1,k,j,l) = at1 + at2
      g(1,ny+2-k,j,l) = at1 - at2
      at2 = g(2,ny+2-k,j,l)
      at1 = g(2,k,j,l) + at2
      at2 = g(2,k,j,l) - at2
      at1 = at3*at1
      at2 = .5*at2
      g(2,k,j,l) = at1 + at2
      g(2,ny+2-k,j,l) = at1 - at2
      at2 = g(3,ny+2-k,j,l)
      at1 = g(3,k,j,l) + at2
      at2 = g(3,k,j,l) - at2
      sum2 = sum2 + real(sctd(k1))*at2
      at2 = at3*at2
      at1 = .5*at1
      g(3,k,j,l) = at1 - at2
      g(3,ny+2-k,j,l) = at1 + at2
   10 continue
      g(1,1,j,l) = 0.0
      g(1,nyh+1,j,l) = 2.0*g(1,nyh+1,j,l)
      g(2,1,j,l) = 0.0
      g(2,nyh+1,j,l) = 2.0*g(2,nyh+1,j,l)
      g(3,1,j,l) = .5*(g(3,1,j,l) + g(3,ny+1,j,l))
      g(3,ny+1,j,l) = sum2
   20 continue
   30 continue
c bit-reverse array elements in y
      nry = nxhy/nyh
      kxpt = kxps
      do 70 l = 1, jblok
      if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
      do 60 k = 1, nyh
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 60
      do 50 j = kxpi, kxpt
      do 40 jj = 1, 3
      t2 = g(jj,2*k1-1,j,l)
      t3 = g(jj,2*k1,j,l)
      g(jj,2*k1-1,j,l) = g(jj,2*k-1,j,l)
      g(jj,2*k1,j,l) = g(jj,2*k,j,l)
      g(jj,2*k-1,j,l) = t2
      g(jj,2*k,j,l) = t3
   40 continue
   50 continue
   60 continue
   70 continue
c first transform in y
      nry = nxy/nyh
      do 130 m = 1, indy1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyhh/ns
      kmr = 2*km*nry
      kxpt = kxps
      do 120 l = 1, jblok
      if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
      do 110 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 100 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 90 i = kxpi, kxpt
      do 80 jj = 1, 3
      t2 = real(t1)*g(jj,2*j2-1,i,l) - aimag(t1)*g(jj,2*j2,i,l)
      t3 = aimag(t1)*g(jj,2*j2-1,i,l) + real(t1)*g(jj,2*j2,i,l)
      g(jj,2*j2-1,i,l) = g(jj,2*j1-1,i,l) - t2
      g(jj,2*j2,i,l) = g(jj,2*j1,i,l) - t3
      g(jj,2*j1-1,i,l) = g(jj,2*j1-1,i,l) + t2
      g(jj,2*j1,i,l) = g(jj,2*j1,i,l) + t3
   80 continue
   90 continue
  100 continue
  110 continue
  120 continue
  130 continue
c unscramble coefficients and normalize
c inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nyh
         ani = 0.5
         kxpt = kxps
         do 190 l = 1, jblok
         if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
         do 160 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         do 150 j = kxpi, kxpt
         do 140 jj = 1, 3
         t4 = g(jj,ny3-2*k,j,l)
         t5 = -g(jj,ny3-2*k+1,j,l)
         t2 = g(jj,2*k-1,j,l) + t4
         t3 = g(jj,2*k,j,l) + t5
         t6 = g(jj,2*k-1,j,l) - t4
         t5 = g(jj,2*k,j,l) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         g(jj,2*k-1,j,l) = ani*(t2 + t4)
         g(jj,2*k,j,l) = ani*(t3 + t5)
         g(jj,ny3-2*k,j,l) = ani*(t2 - t4)
         g(jj,ny3-2*k+1,j,l) = ani*(t5 - t3)
  140    continue
  150    continue
  160    continue
         do 180 j = kxpi, kxpt
         do 170 jj = 1, 3
         g(jj,nyh+1,j,l) = g(jj,nyh+1,j,l)
         g(jj,nyh+2,j,l) = -g(jj,nyh+2,j,l)
         t2 = g(jj,1,j,l) + g(jj,2,j,l)
         g(jj,2,j,l) = g(jj,1,j,l) - g(jj,2,j,l)
         g(jj,1,j,l) = t2
         g(jj,ny+1,j,l) = g(jj,ny+1,j,l)
  170    continue
  180    continue
  190    continue
c forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nyh
         kxpt = kxps
         do 250 l = 1, jblok
         if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
         do 220 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         do 210 j = kxpi, kxpt
         do 200 jj = 1, 3
         t4 = g(jj,ny3-2*k,j,l)
         t5 = -g(jj,ny3-2*k+1,j,l)
         t2 = g(jj,2*k-1,j,l) + t4
         t3 = g(jj,2*k,j,l) + t5
         t6 = g(jj,2*k-1,j,l) - t4
         t5 = g(jj,2*k,j,l) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         g(jj,2*k-1,j,l) = t2 + t4
         g(jj,2*k,j,l) = t3 + t5
         g(jj,ny3-2*k,j,l) = t2 - t4
         g(jj,ny3-2*k+1,j,l) = t5 - t3
  200    continue
  210    continue
  220    continue
         do 240 j = kxpi, kxpt
         do 230 jj = 1, 3
         g(jj,nyh+1,j,l) = 2.0*g(jj,nyh+1,j,l)
         g(jj,nyh+2,j,l) = -2.0*g(jj,nyh+2,j,l)
         t2 = 2.0*(g(jj,1,j,l) + g(jj,2,j,l))
         g(jj,2,j,l) = 2.0*(g(jj,1,j,l) - g(jj,2,j,l))
         g(jj,1,j,l) = t2
         g(jj,ny+1,j,l) = 2.0*g(jj,ny+1,j,l)
  230    continue
  240    continue
  250    continue
      endif
c perform recursion for sine-cosine transform
      kxpt = kxps
      do 280 l = 1, jblok
      if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
      do 270 j = kxpi, kxpt
      sum1 = .5*g(1,1,j,l)
      g(1,1,j,l) = 0.0
      g(1,2,j,l) = sum1
      sum2 = .5*g(2,1,j,l)
      g(2,1,j,l) = 0.0
      g(2,2,j,l) = sum2
      sum3 = g(3,ny+1,j,l)
      g(3,ny+1,j,l) = g(3,2,j,l)
      g(3,2,j,l) = sum3
      do 260 k = 2, nyh
      sum1 = sum1 + g(1,2*k-1,j,l)
      g(1,2*k-1,j,l) = -g(1,2*k,j,l)
      g(1,2*k,j,l) = sum1
      sum2 = sum2 + g(2,2*k-1,j,l)
      g(2,2*k-1,j,l) = -g(2,2*k,j,l)
      g(2,2*k,j,l) = sum2
      sum3 = sum3 - g(3,2*k,j,l)
      g(3,2*k,j,l) = sum3
  260 continue
      g(1,ny+1,j,l) = 0.0
      g(2,ny+1,j,l) = 0.0
  270 continue
  280 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PFSSCT2R3X(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,
     1kypp,nxvh,kypd,kblok,nxhyd,nxyd)
c this subroutine performs the x part of 3 two dimensional fast real
c sine and cosine transforms and their inverses, for a subset of y,
c using real arithmetic, for data which is distributed in blocks
c algorithm is described in Numerical Recipies in Fortran, Second Ed.,
c by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
c [Cambridge Univ. Press, 1992], p. 508.
c for isign = (-1,1), input: all, output: f
c approximate flop count: N*(5*log2(N) + 18)/nvp
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, inverse sine-cosine transforms are performed
c f(1,n,k,i) = (1/nx*ny)*sum(f(1,j,k,i)*sin(pi*n*j/nx))
c f(2:3,n,k,i) = (1/nx*ny)*(.5*f(2:3,1,k,i) + ((-1)**n)*f(2:3,nx+1,k,i)
c              + sum(f(2:3,j,k,i)*cos(pi*n*j/nx)))
c if isign = 1, forward sine transforms are performed
c f(1,j,k,i) = sum(f(1,n,k,i)*sin(pi*n*j/nx))
c f(2:3,j,k,i) = 2*(.5*f(2:3,1,k,i) + ((-1)**j)*f(2:3,n+1,k,i)
c              + sum(f(2:3,n,k,i)*cos(pi*n*j/nx))
c mixup = array of bit reversed addresses
c sctd = sine/cosine table
c kstrt = starting data block number
c kyp = number of data values per block in y
c kypi = initial y index used
c kypp = number of y indices used
c nxvh = first dimension of f >= nx/2 + 1
c kypd = second dimension of f >= kyp + 1
c kblok = number of data blocks in y
c nxhyd = maximum of (nx/2,ny)
c nxyd = maximum of (nx,ny)
c written by viktor k. decyk, ucla
      implicit none
      integer isign, mixup, indx, indy, kstrt, kyp, kypi, kypp
      integer nxvh, kypd, kblok, nxhyd, nxyd
      real f
      complex sctd
      dimension f(3,2*nxvh,kypd,kblok)
      dimension mixup(nxhyd), sctd(nxyd)
c local data
      integer indx1, indx1y, nx, nxh, nxhh, nx3, ny, nxy, nxhy, ks, kypt
      integer i, j, k, l, m, km, kmr, nrx, j1, j2, ns, ns2, k1, k2, kyps
      integer jj
      real at1, at2, at3, t2, t3, t4, t5, t6, ani, sum1, sum2, sum3
      complex t1
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nx3 = nx + 3
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 2
      kyps = kypi + kypp - 1
      if (kstrt.gt.ny) return
      if (isign.eq.0) return
c create auxiliary array in x
      kmr = nxy/nx
      kypt = kyps
      do 30 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 20 k = kypi, kypt
      sum2 = .5*(f(3,1,k,l) - f(3,nx+1,k,l))
      do 10 j = 2, nxh
      j1 = 1 + kmr*(j - 1)
      at3 = -aimag(sctd(j1))
      at2 = f(1,nx+2-j,k,l)
      at1 = f(1,j,k,l) + at2
      at2 = f(1,j,k,l) - at2
      at1 = at3*at1
      at2 = .5*at2
      f(1,j,k,l) = at1 + at2
      f(1,nx+2-j,k,l) = at1 - at2
      at2 = f(2,nx+2-j,k,l)
      at1 = f(2,j,k,l) + at2
      at2 = f(2,j,k,l) - at2
      at1 = at3*at1
      at2 = .5*at2
      f(2,j,k,l) = at1 + at2
      f(2,nx+2-j,k,l) = at1 - at2
      at2 = f(3,nx+2-j,k,l)
      at1 = f(3,j,k,l) + at2
      at2 = f(3,j,k,l) - at2
      sum2 = sum2 + real(sctd(j1))*at2
      at2 = at3*at2
      at1 = .5*at1
      f(3,j,k,l) = at1 - at2
      f(3,nx+2-j,k,l) = at1 + at2
   10 continue
      f(1,1,k,l) = 0.0
      f(1,nxh+1,k,l) = 2.0*f(1,nxh+1,k,l)
      f(2,1,k,l) = 0.0
      f(2,nxh+1,k,l) = 2.0*f(2,nxh+1,k,l)
      f(3,1,k,l) = .5*(f(3,1,k,l) + f(3,nx+1,k,l))
      f(3,nx+1,k,l) = sum2
   20 continue
   30 continue
c bit-reverse array elements in x
      nrx = nxhy/nxh
      kypt = kyps
      do 70 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 60 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 60
      do 50 k = kypi, kypt
      do 40 jj = 1, 3
      t2 = f(jj,2*j1-1,k,l)
      t3 = f(jj,2*j1,k,l)
      f(jj,2*j1-1,k,l) = f(jj,2*j-1,k,l)
      f(jj,2*j1,k,l) = f(jj,2*j,k,l)
      f(jj,2*j-1,k,l) = t2
      f(jj,2*j,k,l) = t3
   40 continue
   50 continue
   60 continue
   70 continue
c first transform in x
      nrx = nxy/nxh
      do 130 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = 2*km*nrx
      kypt = kyps
      do 120 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 110 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 100 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 90 i = kypi, kypt
      do 80 jj = 1, 3
      t2 = real(t1)*f(jj,2*j2-1,i,l) - aimag(t1)*f(jj,2*j2,i,l)
      t3 = aimag(t1)*f(jj,2*j2-1,i,l) + real(t1)*f(jj,2*j2,i,l)
      f(jj,2*j2-1,i,l) = f(jj,2*j1-1,i,l) - t2
      f(jj,2*j2,i,l) = f(jj,2*j1,i,l) - t3
      f(jj,2*j1-1,i,l) = f(jj,2*j1-1,i,l) + t2
      f(jj,2*j1,i,l) = f(jj,2*j1,i,l) + t3
   80 continue
   90 continue
  100 continue
  110 continue
  120 continue
  130 continue
c unscramble coefficients and normalize
c inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nxh
         ani = 1./float(2*nx*ny)
         kypt = kyps
         do 190 l = 1, kblok
         if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
         do 160 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 150 k = kypi, kypt
         do 140 jj = 1, 3
         t4 = f(jj,nx3-2*j,k,l)
         t5 = -f(jj,nx3-2*j+1,k,l)
         t2 = f(jj,2*j-1,k,l) + t4
         t3 = f(jj,2*j,k,l) + t5
         t6 = f(jj,2*j-1,k,l) - t4
         t5 = f(jj,2*j,k,l) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,2*j-1,k,l) = ani*(t2 + t4)
         f(jj,2*j,k,l) = ani*(t3 + t5)
         f(jj,nx3-2*j,k,l) = ani*(t2 - t4)
         f(jj,nx3-2*j+1,k,l) = ani*(t5 - t3)
  140    continue
  150    continue
  160    continue
         ani = 2.*ani
         do 180 k = kypi, kypt
         do 170 jj = 1, 3
         f(jj,nxh+1,k,l) = ani*f(jj,nxh+1,k,l)
         f(jj,nxh+2,k,l) = -ani*f(jj,nxh+2,k,l)
         t2 = ani*(f(jj,1,k,l) + f(jj,2,k,l))
         f(jj,2,k,l) = ani*(f(jj,1,k,l) - f(jj,2,k,l))
         f(jj,1,k,l) = t2
         f(jj,nx+1,k,l) = ani*f(jj,nx+1,k,l)
  170    continue
  180    continue
  190    continue
c forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nxh
         kypt = kyps
         do 250 l = 1, kblok
         if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
         do 220 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 210 k = kypi, kypt
         do 200 jj = 1, 3
         t4 = f(jj,nx3-2*j,k,l)
         t5 = -f(jj,nx3-2*j+1,k,l)
         t2 = f(jj,2*j-1,k,l) + t4
         t3 = f(jj,2*j,k,l) + t5
         t6 = f(jj,2*j-1,k,l) - t4
         t5 = f(jj,2*j,k,l) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,2*j-1,k,l) = t2 + t4
         f(jj,2*j,k,l) = t3 + t5
         f(jj,nx3-2*j,k,l) = t2 - t4
         f(jj,nx3-2*j+1,k,l) = t5 - t3
  200    continue
  210    continue
  220    continue
         do 240 k = kypi, kypt
         do 230 jj = 1, 3
         f(jj,nxh+1,k,l) = 2.0*f(jj,nxh+1,k,l)
         f(jj,nxh+2,k,l) = -2.0*f(jj,nxh+2,k,l)
         t2 = 2.0*(f(jj,1,k,l) + f(jj,2,k,l))
         f(jj,2,k,l) = 2.0*(f(jj,1,k,l) - f(jj,2,k,l))
         f(jj,1,k,l) = t2
         f(jj,nx+1,k,l) = 2.0*f(jj,nx+1,k,l)
  230    continue
  240    continue
  250    continue
      endif
c perform recursion for cosine-sine transform
      kypt = kyps
      do 280 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 270 k = kypi, kypt
      sum1 = .5*f(1,1,k,l)
      f(1,1,k,l) = 0.0
      f(1,2,k,l) = sum1
      sum2 = .5*f(2,1,k,l)
      f(2,1,k,l) = 0.0
      f(2,2,k,l) = sum2
      sum3 = f(3,nx+1,k,l)
      f(3,nx+1,k,l) = f(3,2,k,l)
      f(3,2,k,l) = sum3
      do 260 j = 2, nxh
      sum1 = sum1 + f(1,2*j-1,k,l)
      f(1,2*j-1,k,l) = -f(1,2*j,k,l)
      f(1,2*j,k,l) = sum1
      sum2 = sum2 + f(2,2*j-1,k,l)
      f(2,2*j-1,k,l) = -f(2,2*j,k,l)
      f(2,2*j,k,l) = sum2
      sum3 = sum3 - f(3,2*j,k,l)
      f(3,2*j,k,l) = sum3
  260 continue
      f(1,nx+1,k,l) = 0.0
      f(2,nx+1,k,l) = 0.0
  270 continue
  280 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine WPFS3T2R3(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,in
     1dy,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhyd,nxyd)
c wrapper function for 3 parallel real sine/sine transforms
c for the electric field with dirichlet or magnetic field with neumann
c boundary conditions
      implicit none
      integer isign, ntpose, mixup, indx, indy, kstrt, nxvh, nyv
      integer kxp2, kyp, kypd, kxp2d, jblok, kblok, nxhyd, nxyd
      real f, g, bs, br, ttp
      complex sctd
      dimension f(3,2*nxvh,kypd,kblok), g(3,nyv,kxp2d,jblok)
      dimension bs(3,kxp2+1,kyp+1,kblok), br(3,kxp2+1,kyp+1,jblok)
      dimension mixup(nxhyd), sctd(nxyd)
c local data
      integer nx, ny, kxpi, kypi
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
c inverse fourier transform
      if (isign.lt.0) then
c perform x cosine-sine transform
         call PFSSCT2R3X(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kyp
     1,nxvh,kypd,kblok,nxhyd,nxyd)
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PR3TPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,nyv,kxp2,kyp,kxp2d,k
     1ypd,jblok,kblok)
         call PWTIMERA(1,ttp,dtime)
c perform y sine-cosine transform
         call PFSSCT2R3Y(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kx
     1p2,nyv,kxp2d,jblok,nxhyd,nxyd)
c transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PR3TPOSE(g,f,br,bs,ny,nx,kstrt,nyv,2*nxvh,kyp,kxp2,kypd
     1,kxp2d,kblok,jblok)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PR3TPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,nyv,kxp2,kyp,kxp2
     1d,kypd,jblok,kblok)
            call PWTIMERA(1,tf,dtime)
         endif
c perform y sine-cosine transform
         call PFSSCT2R3Y(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kx
     1p2,nyv,kxp2d,jblok,nxhyd,nxyd)
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PR3TPOSE(g,f,br,bs,ny,nx,kstrt,nyv,2*nxvh,kyp,kxp2,kypd,kx
     1p2d,kblok,jblok)
         call PWTIMERA(1,ttp,dtime)
c perform x cosine-sine transform
         call PFSSCT2R3X(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kyp
     1,nxvh,kypd,kblok,nxhyd,nxyd)
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine PRTPOSE(f,g,s,t,nx,ny,kstrt,nxv,nyv,kxp,kyp,kxpd,kypd,j
     1blok,kblok)
c this subroutine performs a transpose of a matrix f, distributed in y,
c to a matrix g, distributed in x, that is,
c g(k+kyp*(m-1),j,l) = f(j+kxp*(l-1),k,m), where
c 1 <= j <= kxp, 1 <= k <= kyp, 1 <= l <= nx/kxp, 1 <= m <= ny/kyp
c and where indices l and m can be distributed across processors.
c includes an extra guard cell for last row and column
c this subroutine sends and receives one message at a time, either
c synchronously or asynchronously. it uses a minimum of system resources
c f = real input array
c g = real output array
c s, t = real scratch arrays
c nx/ny = number of points in x/y
c kstrt = starting data block number
c nxv = first dimension of f >= nx+1
c nyv = first dimension of g >= ny+1
c kypd = second dimension of f >= kyp+1
c kxpd = second dimension of g >= kxp+1
c kxp/kyp = number of data values per block in x/y
c jblok/kblok = number of data blocks in x/y
      implicit none
      integer nx, ny, kstrt, nxv, nyv, kxp, kyp, kxpd, kypd
      integer jblok, kblok
      real f, g, s, t
      dimension f(nxv,kypd,kblok), g(nyv,kxpd,jblok)
      dimension s(kxp+1,kyp+1,kblok), t(kxp+1,kyp+1,jblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer ks, kxb, kyb, kxp1, kyp1, kxpt, kypt
      integer jkblok, kxym, mtr, ntr, mntr
      integer l, i, joff, koff, k, j
      integer ir0, is0, ii, ir, is, ierr, msid, istatus
      dimension istatus(lstat)
      ks = kstrt - 2
      kxb = nx/kxp
      kyb = ny/kyp
c set constants to receive extra guard cells
      kxp1 = kxp + 1
      kyp1 = kyp + 1
      kxpt = kxp
      if (kstrt.eq.kxb) kxpt = kxp1
c this segment is used for shared memory computers
c     if (kstrt.gt.nx) return
c     kypt = kyp
c     do 40 l = 1, jblok
c     joff = kxp*(l + ks)
c     if ((l+ks).eq.(kxb-1)) kxpt = kxp1
c     do 30 i = 1, kyb
c     koff = kyp*(i - 1)
c     if (i.eq.kyb) kypt = kyp1
c     do 20 k = 1, kypt
c     do 10 j = 1, kxpt
c     g(k+koff,j,l) = f(j+joff,k,i)
c  10 continue
c  20 continue
c  30 continue
c  40 continue
c this segment is used for mpi computers
      jkblok = max0(jblok,kblok)
      kxym = min0(kxb,kyb)
      mtr = kyb/kxym
      ntr = kxb/kxym
      mntr = max0(mtr,ntr)
      do 70 l = 1, jkblok
      do 60 i = 1, kxym
      ir0 = iand(kxym-1,ieor(l+ks,i-1)) + 1
      is0 = ir0
      do 50 ii = 1, mntr
c post receive
      if ((kstrt.le.nx).and.(ii.le.mtr)) then
         ir = ir0 + kxym*(ii - 1)
         kypt = kyp
         if (ir.eq.kyb) kypt = kyp1
         call MPI_IRECV(t(1,1,l),kxp1*kyp1,mreal,ir-1,ir+kxym+1,lgrp,msi
     1d,ierr)
      endif
c send data
      if ((kstrt.le.ny).and.(ii.le.ntr)) then
         is = is0 + kxym*(ii - 1)
         joff = kxp*(is - 1)
         do 20 k = 1, kyp1
         do 10 j = 1, kxp1
         s(j,k,l) = f(j+joff,k,l)
   10    continue
   20    continue
         call MPI_SEND(s(1,1,l),kxp1*kyp1,mreal,is-1,l+ks+kxym+2,lgrp,ie
     1rr)
      endif
c receive data
      if ((kstrt.le.nx).and.(ii.le.mtr)) then
         koff = kyp*(ir - 1)
         call MPI_WAIT(msid,istatus,ierr)
         do 40 k = 1, kypt
         do 30 j = 1, kxpt
         g(k+koff,j,l) = t(j,k,l)
   30    continue
   40    continue
      endif
   50 continue
   60 continue
   70 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PR2TPOSE(f,g,s,t,nx,ny,kstrt,nxv,nyv,kxp,kyp,kxpd,kypd,
     1jblok,kblok)
c this subroutine performs a transpose of a matrix f, distributed in y,
c to a matrix g, distributed in x, that is,
c g(1:2,k+kyp*(m-1),j,l) = f(1:2,j+kxp*(l-1),k,m), where
c 1 <= j <= kxp, 1 <= k <= kyp, 1 <= l <= nx/kxp, 1 <= m <= ny/kyp
c and where indices l and m can be distributed across processors.
c includes an extra guard cell for last row and column
c this subroutine sends and receives one message at a time, either
c synchronously or asynchronously. it uses a minimum of system resources
c f = real input array
c g = real output array
c s, t = real scratch arrays
c nx/ny = number of points in x/y
c kstrt = starting data block number
c nxv = first dimension of f >= nx+1
c nyv = first dimension of g >= ny+1
c kypd = second dimension of f >= kyp+1
c kxpd = second dimension of g >= kxp+1
c kxp/kyp = number of data values per block in x/y
c jblok/kblok = number of data blocks in x/y
      implicit none
      integer nx, ny, kstrt, nxv, nyv, kxp, kyp, kxpd, kypd
      integer jblok, kblok
      real f, g, s, t
      dimension f(2,nxv,kypd,kblok), g(2,nyv,kxpd,jblok)
      dimension s(2,kxp+1,kyp+1,kblok), t(2,kxp+1,kyp+1,jblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer ks, kxb, kyb, kxp1, kyp1, kxpt, kypt
      integer jkblok, kxym, mtr, ntr, mntr
      integer l, i, joff, koff, k, j
      integer ir0, is0, ii, ir, is, ierr, msid, istatus
      dimension istatus(lstat)
      ks = kstrt - 2
      kxb = nx/kxp
      kyb = ny/kyp
c set constants to receive extra guard cells
      kxp1 = kxp + 1
      kyp1 = kyp + 1
      kxpt = kxp
      if (kstrt.eq.kxb) kxpt = kxp1
c this segment is used for shared memory computers
c     if (kstrt.gt.nx) return
c     kypt = kyp
c     do 40 l = 1, jblok
c     joff = kxp*(l + ks)
c     if ((l+ks).eq.(kxb-1)) kxpt = kxp1
c     do 30 i = 1, kyb
c     koff = kyp*(i - 1)
c     if (i.eq.kyb) kypt = kyp1
c     do 20 k = 1, kypt
c     do 10 j = 1, kxpt
c     g(1,k+koff,j,l) = f(1,j+joff,k,i)
c     g(2,k+koff,j,l) = f(2,j+joff,k,i)
c  10 continue
c  20 continue
c  30 continue
c  40 continue
c this segment is used for mpi computers
      jkblok = max0(jblok,kblok)
      kxym = min0(kxb,kyb)
      mtr = kyb/kxym
      ntr = kxb/kxym
      mntr = max0(mtr,ntr)
      do 70 l = 1, jkblok
      do 60 i = 1, kxym
      ir0 = iand(kxym-1,ieor(l+ks,i-1)) + 1
      is0 = ir0
      do 50 ii = 1, mntr
c post receive
      if ((kstrt.le.nx).and.(ii.le.mtr)) then
         ir = ir0 + kxym*(ii - 1)
         kypt = kyp
         if (ir.eq.kyb) kypt = kyp1
         call MPI_IRECV(t(1,1,1,l),2*kxp1*kyp1,mreal,ir-1,ir+kxym+1,lgrp
     1,msid,ierr)
      endif
c send data
      if ((kstrt.le.ny).and.(ii.le.ntr)) then
         is = is0 + kxym*(ii - 1)
         joff = kxp*(is - 1)
         do 20 k = 1, kyp1
         do 10 j = 1, kxp1
         s(1,j,k,l) = f(1,j+joff,k,l)
         s(2,j,k,l) = f(2,j+joff,k,l)
   10    continue
   20    continue
         call MPI_SEND(s(1,1,1,l),2*kxp1*kyp1,mreal,is-1,l+ks+kxym+2,lgr
     1p,ierr)
      endif
c receive data
      if ((kstrt.le.nx).and.(ii.le.mtr)) then
         koff = kyp*(ir - 1)
         call MPI_WAIT(msid,istatus,ierr)
         do 40 k = 1, kypt
         do 30 j = 1, kxpt
         g(1,k+koff,j,l) = t(1,j,k,l)
         g(2,k+koff,j,l) = t(2,j,k,l)
   30    continue
   40    continue
      endif
   50 continue
   60 continue
   70 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PR3TPOSE(f,g,s,t,nx,ny,kstrt,nxv,nyv,kxp,kyp,kxpd,kypd,
     1jblok,kblok)
c this subroutine performs a transpose of a matrix f, distributed in y,
c to a matrix g, distributed in x, that is,
c g(1:3,k+kyp*(m-1),j,l) = f(1:3,j+kxp*(l-1),k,m), where
c 1 <= j <= kxp, 1 <= k <= kyp, 1 <= l <= nx/kxp, 1 <= m <= ny/kyp
c and where indices l and m can be distributed across processors.
c includes an extra guard cell for last row and column
c this subroutine sends and receives one message at a time, either
c synchronously or asynchronously. it uses a minimum of system resources
c f = real input array
c g = real output array
c s, t = real scratch arrays
c nx/ny = number of points in x/y
c kstrt = starting data block number
c nxv = first dimension of f >= nx+1
c nyv = first dimension of g >= ny+1
c kypd = second dimension of f >= kyp+1
c kxpd = second dimension of g >= kxp+1
c kxp/kyp = number of data values per block in x/y
c jblok/kblok = number of data blocks in x/y
      implicit none
      integer nx, ny, kstrt, nxv, nyv, kxp, kyp, kxpd, kypd
      integer jblok, kblok
      real f, g, s, t
      dimension f(3,nxv,kypd,kblok), g(3,nyv,kxpd,jblok)
      dimension s(3,kxp+1,kyp+1,kblok), t(3,kxp+1,kyp+1,jblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer ks, kxb, kyb, kxp1, kyp1, kxpt, kypt
      integer jkblok, kxym, mtr, ntr, mntr
      integer l, i, joff, koff, k, j
      integer ir0, is0, ii, ir, is, ierr, msid, istatus
      dimension istatus(lstat)
      ks = kstrt - 2
      kxb = nx/kxp
      kyb = ny/kyp
c set constants to receive extra guard cells
      kxp1 = kxp + 1
      kyp1 = kyp + 1
      kxpt = kxp
      if (kstrt.eq.kxb) kxpt = kxp1
c this segment is used for shared memory computers
c     if (kstrt.gt.nx) return
c     kypt = kyp
c     do 40 l = 1, jblok
c     joff = kxp*(l + ks)
c     if ((l+ks).eq.(kxb-1)) kxpt = kxp1
c     do 30 i = 1, kyb
c     koff = kyp*(i - 1)
c     if (i.eq.kyb) kypt = kyp1
c     do 20 k = 1, kypt
c     do 10 j = 1, kxpt
c     g(1,k+koff,j,l) = f(1,j+joff,k,i)
c     g(2,k+koff,j,l) = f(2,j+joff,k,i)
c     g(3,k+koff,j,l) = f(3,j+joff,k,i)
c  10 continue
c  20 continue
c  30 continue
c  40 continue
c this segment is used for mpi computers
      jkblok = max0(jblok,kblok)
      kxym = min0(kxb,kyb)
      mtr = kyb/kxym
      ntr = kxb/kxym
      mntr = max0(mtr,ntr)
      do 70 l = 1, jkblok
      do 60 i = 1, kxym
      ir0 = iand(kxym-1,ieor(l+ks,i-1)) + 1
      is0 = ir0
      do 50 ii = 1, mntr
c post receive
      if ((kstrt.le.nx).and.(ii.le.mtr)) then
         ir = ir0 + kxym*(ii - 1)
         kypt = kyp
         if (ir.eq.kyb) kypt = kyp1
         call MPI_IRECV(t(1,1,1,l),3*kxp1*kyp1,mreal,ir-1,ir+kxym+1,lgrp
     1,msid,ierr)
      endif
c send data
      if ((kstrt.le.ny).and.(ii.le.ntr)) then
         is = is0 + kxym*(ii - 1)
         joff = kxp*(is - 1)
         do 20 k = 1, kyp1
         do 10 j = 1, kxp1
         s(1,j,k,l) = f(1,j+joff,k,l)
         s(2,j,k,l) = f(2,j+joff,k,l)
         s(3,j,k,l) = f(3,j+joff,k,l)
   10    continue
   20    continue
         call MPI_SEND(s(1,1,1,l),3*kxp1*kyp1,mreal,is-1,l+ks+kxym+2,lgr
     1p,ierr)
      endif
c receive data
      if ((kstrt.le.nx).and.(ii.le.mtr)) then
         koff = kyp*(ir - 1)
         call MPI_WAIT(msid,istatus,ierr)
         do 40 k = 1, kypt
         do 30 j = 1, kxpt
         g(1,k+koff,j,l) = t(1,j,k,l)
         g(2,k+koff,j,l) = t(2,j,k,l)
         g(3,k+koff,j,l) = t(3,j,k,l)
   30    continue
   40    continue
      endif
   50 continue
   60 continue
   70 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PWTIMERA(icntrl,time,dtime)
c this subroutine performs local wall clock timing
c input: icntrl, dtime
c icntrl = (-1,0,1) = (initialize,ignore,read) clock
c clock should be initialized before it is read!
c time = elapsed time in seconds
c dtime = current time
c written for mpi
      implicit none
      integer icntrl
      real time
      double precision dtime
c local data
      double precision jclock
      double precision MPI_WTIME
      external MPI_WTIME
c initialize clock
      if (icntrl.eq.(-1)) then
         dtime = MPI_WTIME()
c read clock and write time difference from last clock initialization
      else if (icntrl.eq.1) then
         jclock = dtime
         dtime = MPI_WTIME()
         time = real(dtime - jclock)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PDIVFD2(f,df,nx,ny,kstrt,ndim,nyv,kxp2,j2blok)
c this subroutine calculates the divergence in fourier space
c with dirichlet boundary conditions (zero potential)
c intended for calculating the charge density from the electric field
c input: all except df, output: df
c approximate flop count is: 6*nxc*nyc
c where nxc = (nx/2-1)/nvp, nyc = ny - 1, and nvp = number of procs
c the divergence is calculated using the equation:
c df(kx,ky) = sqrt(-1)*(kx*fx(kx,ky)+ky*fy(kx,ky))
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c except for df(kx=pi) = df(ky=pi) = df(kx=0,ky=0) = 0.
c nx/ny = system length in x/y direction
c ndim = number of field arrays, must be >= 2
c kstrt = starting data block number
c nyv = first dimension of field arrays, must be >= ny+1
c kxp2 = number of data values per block
c j2blok = number of data blocks
      real f, df
      dimension f(ndim,nyv,kxp2+1,j2blok), df(nyv,kxp2+1,j2blok)
      if (ndim.lt.2) return
      ks = kstrt - 2
      ny1 = ny + 1
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
c calculate the divergence
      if (kstrt.gt.nx) return
      do 50 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 20 j = 1, kxp2
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 10 k = 2, ny
         dky = dny*float(k - 1)
         df(k,j,l) = -(dkx*f(1,k,j,l) + dky*f(2,k,j,l))
   10    continue
      endif
c mode numbers ky = 0, ny
      df(1,j,l) = 0.
      df(ny+1,j,l) = 0.
   20 continue
c mode numbers kx = 0, nx
      if ((l+ks).eq.0) then
         do 30 k = 2, ny
         df(k,1,l) = 0.
   30    continue
      endif
      do 40 k = 1, ny1
      df(k,kxp2+1,l) = 0.
   40 continue
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PGRADFD2(df,f,nx,ny,kstrt,ndim,nyv,kxp2,j2blok)
c this subroutine calculates the gradient in fourier space
c with dirichlet boundary conditions (zero potential)
c intended for calculating the electric field from the potential
c input: all except f, output: f
c approximate flop count is: 4*nxc*nyc
c where nxc = (nx/2-1)/nvp, nyc = ny - 1, and nvp = number of procs
c the gradient is calculated using the equations:
c fx(kx,ky) = sqrt(-1)*kx*df(kx,ky)
c fy(kx,ky) = sqrt(-1)*ky*df(kx,ky)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c except for fx(kx=pi) = fy(kx=pi) = 0, fx(ky=pi) = fy(ky=pi) = 0,
c and fx(kx=0,ky=0) = fy(kx=0,ky=0) = 0.
c nx/ny = system length in x/y direction
c ndim = number of field arrays, must be >= 2
c kstrt = starting data block number
c nyv = first dimension of field arrays, must be >= ny+1
c kxp2 = number of data values per block
c j2blok = number of data blocks
      real df, f
      dimension df(nyv,kxp2+1,j2blok), f(ndim,nyv,kxp2+1,j2blok)
      ks = kstrt - 2
      ny1 = ny + 1
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
c calculate the gradient
      if (kstrt.gt.nx) return
      do 50 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 20 j = 1, kxp2
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 10 k = 2, ny
         dky = dny*float(k - 1)
         f(1,k,j,l) = dkx*df(k,j,l)
         f(2,k,j,l) = dky*df(k,j,l)
   10    continue
      endif
c mode numbers ky = 0, ny
      f(1,1,j,l) = 0.
      f(2,1,j,l) = 0.
      f(1,ny+1,j,l) = 0.
      f(2,ny+1,j,l) = 0.
   20 continue
c mode numbers kx = 0, nx
      if ((l+ks).eq.0) then
         do 30 k = 2, ny
         f(1,k,1,l) = 0.
         f(2,k,1,l) = 0.
   30    continue
      endif
      do 40 k = 1, ny1
      f(1,k,kxp2+1,l) = 0.
      f(2,k,kxp2+1,l) = 0.
   40 continue
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PCURLFD22(f,g,nx,ny,kstrt,nyv,kxp2,j2blok)
c this subroutine calculates the curl in fourier space
c with dirichlet boundary conditions (zero potential)
c intended for calculating the magnetic field from the vector potential
c input: all except g, output: g
c approximate flop count is: 32*nxc*nyc
c where nxc = (nx/2-1)/nvp, nyc = ny - 1, and nvp = number of procs
c the curl is calculated using the equations:
c g(kx,ky) = sqrt(-1)*(kx*fy(kx,ky)-ky*fx(kx,ky))
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c nx/ny = system length in x/y direction
c kstrt = starting data block number
c nyv = first dimension of field arrays, must be >= ny+1
c kxp2 = number of data values per block
c j2blok = number of data blocks
      real f, g
      dimension f(2,nyv,kxp2+1,j2blok), g(nyv,kxp2+1,j2blok)
      ks = kstrt - 2
      ny1 = ny + 1
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
c calculate the curl
      if (kstrt.gt.nx) return
      do 50 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 20 j = 1, kxp2
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 10 k = 2, ny
         dky = dny*float(k - 1)
         g(k,j,l) = dkx*f(2,k,j,l) - dky*f(1,k,j,l)
   10    continue
c mode numbers ky = 0, ny
         g(1,j,l) = dkx*f(2,1,j,l)
      endif
      g(ny+1,j,l) = 0.
   20 continue
c mode numbers kx = 0, nx
      if ((l+ks).eq.0) then
         do 30 k = 2, ny
         dky = dny*float(k - 1)
         g(k,1,l) = -dky*f(1,k,1,l)
   30    continue
         g(1,1,l) = 0.
      endif
      do 40 k = 1, ny1
      g(k,kxp2+1,l) = 0.
   40 continue
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PCURLFD2(f,g,nx,ny,kstrt,nyv,kxp2,j2blok)
c this subroutine calculates the curl in fourier space
c with dirichlet boundary conditions (zero potential)
c intended for calculating the magnetic field from the vector potential
c input: all except g, output: g
c approximate flop count is: 8*nxc*nyc
c where nxc = (nx/2-1)/nvp, nyc = ny - 1, and nvp = number of procs
c the curl is calculated using the equations:
c gx(kx,ky) = sqrt(-1)*ky*fz(kx,ky)
c gy(kx,ky) = -sqrt(-1)*kx*fz(kx,ky)
c gz(kx,ky) = sqrt(-1)*(kx*fy(kx,ky)-ky*fx(kx,ky))
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c except for gx(kx=pi) = gy(kx=pi) = 0, gx(ky=pi) = gy(ky=pi) = 0,
c and gx(kx=0,ky=0) = gy(kx=0,ky=0) = 0.
c nx/ny = system length in x/y direction
c kstrt = starting data block number
c nyv = first dimension of field arrays, must be >= ny+1
c kxp2 = number of data values per block
c j2blok = number of data blocks
      real f, g
      dimension f(3,nyv,kxp2+1,j2blok), g(3,nyv,kxp2+1,j2blok)
      ks = kstrt - 2
      ny1 = ny + 1
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
c calculate the curl
      if (kstrt.gt.nx) return
      do 50 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 20 j = 1, kxp2
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 10 k = 2, ny
         dky = dny*float(k - 1)
         g(1,k,j,l) = dky*f(3,k,j,l)
         g(2,k,j,l) = -dkx*f(3,k,j,l)
         g(3,k,j,l) = dkx*f(2,k,j,l) - dky*f(1,k,j,l)
   10    continue
c mode numbers ky = 0, ny
         g(1,1,j,l) = 0.
         g(2,1,j,l) = 0.
         g(3,1,j,l) = dkx*f(2,1,j,l)
      endif
      g(1,ny+1,j,l) = 0.
      g(2,ny+1,j,l) = 0.
      g(3,ny+1,j,l) = 0.
   20 continue
c mode numbers kx = 0, nx
      if ((l+ks).eq.0) then
         do 30 k = 2, ny
         dky = dny*float(k - 1)
         g(1,k,1,l) = 0.
         g(2,k,1,l) = 0.
         g(3,k,1,l) = -dky*f(1,k,1,l)
   30    continue
         g(1,1,1,l) = 0.
         g(2,1,1,l) = 0.
         g(3,1,1,l) = 0.
      endif
      do 40 k = 1, ny1
      g(1,k,kxp2+1,l) = 0.
      g(2,k,kxp2+1,l) = 0.
      g(3,k,kxp2+1,l) = 0.
   40 continue
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine WPPFSST2RM(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx, &
     &indy,kstrt,nvp,nxvh,nyv,kxp2,kyp,kypd,kxp2d,nxhyd,nxyd)
! wrapper function for parallel real sine/sine transform
      implicit none
      integer isign, ntpose, indx, indy, kstrt, nvp, nxvh, nyv
      integer kxp2, kyp, kypd, kxp2d, nxhyd, nxyd, mixup
      real f, g, bs, br, ttp
      complex sctd
      dimension f(2*nxvh,kypd), g(nyv,kxp2d)
      dimension bs(kxp2+1,kyp+1), br(kxp2+1,kyp+1)
      dimension mixup(nxhyd), sctd(nxyd)
! local data
      integer nx, ny, kxpi, kypi, ks, kxp2p, kypp, kxb2, kyb
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
! calculate range of indices
      nx = 2**indx
      ny = 2**indy
! ks = processor id
      ks = kstrt - 1
! kxp2p = actual size used in x direction
      kxp2p = min(kxp2,max(0,nx-kxp2*ks))
! kypp = actual size used in y direction
      kypp = min(kyp,max(0,ny-kyp*ks))
! kxb2 = minimum number of processors needed in x direction
      kxb2 = (nx - 1)/kxp2 + 1
! kyb = minimum number of processors needed in y direction
      kyb = (ny - 1)/kyp + 1
! add extra word for last processor in x
      if (ks==(kxb2-1)) kxp2p = kxp2p + 1
! add extra word for last processor in y
      if (ks==(kyb-1)) kypp = kypp + 1
! inverse fourier transform
      if (isign.lt.0) then
! perform x sine transform
         call PPFST2RMXX(f,isign,mixup,sctd,indx,indy,kstrt,kypi,kypp,  &
     &nxvh,kypd,nxhyd,nxyd)
! transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PPRTPOSE(f,g,bs,br,nx,ny,kxp2,kyp,kstrt,nvp,2*nxvh,nyv,   &
     &kxp2d,kypd)
         call PWTIMERA(1,ttp,dtime)
! perform y sine transform
         call PPFST2RMXY(g,isign,mixup,sctd,indx,indy,kstrt,kxpi,kxp2p, &
     &nyv,kxp2d,nxhyd,nxyd)
! transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPRTPOSE(g,f,br,bs,ny,nx,kyp,kxp2,kstrt,nvp,nyv,2*nxvh,&
     &kypd,kxp2d)
            call PWTIMERA(1,tf,dtime)
         endif
! forward fourier transform
      else if (isign.gt.0) then
! transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPRTPOSE(f,g,bs,br,nx,ny,kxp2,kyp,kstrt,nvp,2*nxvh,nyv,&
     &kxp2d,kypd)
            call PWTIMERA(1,tf,dtime)
         endif
! perform y sine transform
         call PPFST2RMXY(g,isign,mixup,sctd,indx,indy,kstrt,kxpi,kxp2p, &
     &nyv,kxp2d,nxhyd,nxyd)
! transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PPRTPOSE(g,f,br,bs,ny,nx,kyp,kxp2,kstrt,nvp,nyv,2*nxvh,   &
     &kypd,kxp2d)
         call PWTIMERA(1,ttp,dtime)
! perform x sine transform
         call PPFST2RMXX(f,isign,mixup,sctd,indx,indy,kstrt,kypi,kypp,  &
     &nxvh,kypd,nxhyd,nxyd)
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
!-----------------------------------------------------------------------
      subroutine WPPFSCT2RM(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx, &
     &indy,kstrt,nvp,nxvh,nyv,kxp2,kyp,kypd,kxp2d,nxhyd,nxyd)
! wrapper function for parallel real sine/cosine transform
      implicit none
      integer isign, ntpose, indx, indy, kstrt, nvp, nxvh, nyv
      integer kxp2, kyp, kypd, kxp2d, nxhyd, nxyd, mixup
      real f, g, bs, br, ttp
      complex sctd
      dimension f(2*nxvh,kypd), g(nyv,kxp2d)
      dimension bs(kxp2+1,kyp+1), br(kxp2+1,kyp+1)
      dimension mixup(nxhyd), sctd(nxyd)
! local data
      integer nx, ny, kxpi, kypi, ks, kxp2p, kypp, kxb2, kyb
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
! calculate range of indices
      nx = 2**indx
      ny = 2**indy
! ks = processor id
      ks = kstrt - 1
! kxp2p = actual size used in x direction
      kxp2p = min(kxp2,max(0,nx-kxp2*ks))
! kypp = actual size used in y direction
      kypp = min(kyp,max(0,ny-kyp*ks))
! kxb2 = minimum number of processors needed in x direction
      kxb2 = (nx - 1)/kxp2 + 1
! kyb = minimum number of processors needed in y direction
      kyb = (ny - 1)/kyp + 1
! add extra word for last processor in x
      if (ks==(kxb2-1)) kxp2p = kxp2p + 1
! add extra word for last processor in y
      if (ks==(kyb-1)) kypp = kypp + 1
! inverse fourier transform
      if (isign.lt.0) then
! perform x sine transform
         call PPFST2RMXX(f,isign,mixup,sctd,indx,indy,kstrt,kypi,kypp,  &
     &nxvh,kypd,nxhyd,nxyd)
! transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PPRTPOSE(f,g,bs,br,nx,ny,kxp2,kyp,kstrt,nvp,2*nxvh,nyv,   &
     &kxp2d,kypd)
         call PWTIMERA(1,ttp,dtime)
! perform y cosine transform
         call PPFCT2RMXY(g,isign,mixup,sctd,indx,indy,kstrt,kxpi,kxp2p, &
     &nyv,kxp2d,nxhyd,nxyd)
! transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPRTPOSE(g,f,br,bs,ny,nx,kyp,kxp2,kstrt,nvp,nyv,2*nxvh,&
     &kypd,kxp2d)
            call PWTIMERA(1,tf,dtime)
         endif
! forward fourier transform
      else if (isign.gt.0) then
! transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPRTPOSE(f,g,bs,br,nx,ny,kxp2,kyp,kstrt,nvp,2*nxvh,nyv,&
     &kxp2d,kypd)
            call PWTIMERA(1,tf,dtime)
         endif
! perform y cosine transform
         call PPFCT2RMXY(g,isign,mixup,sctd,indx,indy,kstrt,kxpi,kxp2p, &
     &nyv,kxp2d,nxhyd,nxyd)
! transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PPRTPOSE(g,f,br,bs,ny,nx,kyp,kxp2,kstrt,nvp,nyv,2*nxvh,   &
     &kypd,kxp2d)
         call PWTIMERA(1,ttp,dtime)
! perform x sine transform
         call PPFST2RMXX(f,isign,mixup,sctd,indx,indy,kstrt,kypi,kypp,  &
     &nxvh,kypd,nxhyd,nxyd)
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
!-----------------------------------------------------------------------
      subroutine WPPFCST2RM(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx, &
     &indy,kstrt,nvp,nxvh,nyv,kxp2,kyp,kypd,kxp2d,nxhyd,nxyd)
! wrapper function for parallel real cosine/sine transform
      implicit none
      integer isign, ntpose, indx, indy, kstrt, nvp, nxvh, nyv
      integer kxp2, kyp, kypd, kxp2d, nxhyd, nxyd, mixup
      real f, g, bs, br, ttp
      complex sctd
      dimension f(2*nxvh,kypd), g(nyv,kxp2d)
      dimension bs(kxp2+1,kyp+1), br(kxp2+1,kyp+1)
      dimension mixup(nxhyd), sctd(nxyd)
! local data
      integer nx, ny, kxpi, kypi, ks, kxp2p, kypp, kxb2, kyb
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
! calculate range of indices
      nx = 2**indx
      ny = 2**indy
! ks = processor id
      ks = kstrt - 1
! kxp2p = actual size used in x direction
      kxp2p = min(kxp2,max(0,nx-kxp2*ks))
! kypp = actual size used in y direction
      kypp = min(kyp,max(0,ny-kyp*ks))
! kxb2 = minimum number of processors needed in x direction
      kxb2 = (nx - 1)/kxp2 + 1
! kyb = minimum number of processors needed in y direction
      kyb = (ny - 1)/kyp + 1
! add extra word for last processor in x
      if (ks==(kxb2-1)) kxp2p = kxp2p + 1
! add extra word for last processor in y
      if (ks==(kyb-1)) kypp = kypp + 1
! inverse fourier transform
      if (isign.lt.0) then
! perform x cosine transform
         call PPFCT2RMXX(f,isign,mixup,sctd,indx,indy,kstrt,kypi,kypp,  &
     &nxvh,kypd,nxhyd,nxyd)
! transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PPRTPOSE(f,g,bs,br,nx,ny,kxp2,kyp,kstrt,nvp,2*nxvh,nyv,   &
     &kxp2d,kypd)
         call PWTIMERA(1,ttp,dtime)
! perform y sine transform
         call PPFST2RMXY(g,isign,mixup,sctd,indx,indy,kstrt,kxpi,kxp2p, &
     &nyv,kxp2d,nxhyd,nxyd)
! transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPRTPOSE(g,f,br,bs,ny,nx,kyp,kxp2,kstrt,nvp,nyv,2*nxvh,&
     &kypd,kxp2d)
            call PWTIMERA(1,tf,dtime)
         endif
! forward fourier transform
      else if (isign.gt.0) then
! transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPRTPOSE(f,g,bs,br,nx,ny,kxp2,kyp,kstrt,nvp,2*nxvh,nyv,&
     &kxp2d,kypd)
            call PWTIMERA(1,tf,dtime)
         endif
! perform y sine transform
         call PPFST2RMXY(g,isign,mixup,sctd,indx,indy,kstrt,kxpi,kxp2p, &
     &nyv,kxp2d,nxhyd,nxyd)
! transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PPRTPOSE(g,f,br,bs,ny,nx,kyp,kxp2,kstrt,nvp,nyv,2*nxvh,   &
     &kypd,kxp2d)
         call PWTIMERA(1,ttp,dtime)
! perform x cosine transform
         call PPFCT2RMXX(f,isign,mixup,sctd,indx,indy,kstrt,kypi,kypp,  &
     &nxvh,kypd,nxhyd,nxyd)
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
!-----------------------------------------------------------------------
      subroutine WPPFCCT2RM(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx, &
     &indy,kstrt,nvp,nxvh,nyv,kxp2,kyp,kypd,kxp2d,nxhyd,nxyd)
! wrapper function for parallel real cosine/cosine transform
      implicit none
      integer isign, ntpose, indx, indy, kstrt, nvp, nxvh, nyv
      integer kxp2, kyp, kypd, kxp2d, nxhyd, nxyd, mixup
      real f, g, bs, br, ttp
      complex sctd
      dimension f(2*nxvh,kypd), g(nyv,kxp2d)
      dimension bs(kxp2+1,kyp+1), br(kxp2+1,kyp+1)
      dimension mixup(nxhyd), sctd(nxyd)
! local data
      integer nx, ny, kxpi, kypi, ks, kxp2p, kypp, kxb2, kyb
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
! calculate range of indices
      nx = 2**indx
      ny = 2**indy
! ks = processor id
      ks = kstrt - 1
! kxp2p = actual size used in x direction
      kxp2p = min(kxp2,max(0,nx-kxp2*ks))
! kypp = actual size used in y direction
      kypp = min(kyp,max(0,ny-kyp*ks))
! kxb2 = minimum number of processors needed in x direction
      kxb2 = (nx - 1)/kxp2 + 1
! kyb = minimum number of processors needed in y direction
      kyb = (ny - 1)/kyp + 1
! add extra word for last processor in x
      if (ks==(kxb2-1)) kxp2p = kxp2p + 1
! add extra word for last processor in y
      if (ks==(kyb-1)) kypp = kypp + 1
! inverse fourier transform
      if (isign.lt.0) then
! perform x cosine transform
         call PPFCT2RMXX(f,isign,mixup,sctd,indx,indy,kstrt,kypi,kypp,  &
     &nxvh,kypd,nxhyd,nxyd)
! transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PPRTPOSE(f,g,bs,br,nx,ny,kxp2,kyp,kstrt,nvp,2*nxvh,nyv,   &
     &kxp2d,kypd)
         call PWTIMERA(1,ttp,dtime)
! perform y cosine transform
         call PPFCT2RMXY(g,isign,mixup,sctd,indx,indy,kstrt,kxpi,kxp2p, &
     &nyv,kxp2d,nxhyd,nxyd)
! transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPRTPOSE(g,f,br,bs,ny,nx,kyp,kxp2,kstrt,nvp,nyv,2*nxvh,&
     &kypd,kxp2d)
            call PWTIMERA(1,tf,dtime)
         endif
! forward fourier transform
      else if (isign.gt.0) then
! transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPRTPOSE(f,g,bs,br,nx,ny,kxp2,kyp,kstrt,nvp,2*nxvh,nyv,&
     &kxp2d,kypd)
            call PWTIMERA(1,tf,dtime)
         endif
! perform y cosine transform
         call PPFCT2RMXY(g,isign,mixup,sctd,indx,indy,kstrt,kxpi,kxp2p, &
     &nyv,kxp2d,nxhyd,nxyd)
! transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PPRTPOSE(g,f,br,bs,ny,nx,kyp,kxp2,kstrt,nvp,nyv,2*nxvh,   &
     &kypd,kxp2d)
         call PWTIMERA(1,ttp,dtime)
! perform x cosine transform
         call PPFCT2RMXX(f,isign,mixup,sctd,indx,indy,kstrt,kypi,kypp,  &
     &nxvh,kypd,nxhyd,nxyd)
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
!-----------------------------------------------------------------------
      subroutine PPFST2RMXX(f,isign,mixup,sctd,indx,indy,kstrt,kypi,kypp&
     &,nxvh,kypd,nxhyd,nxyd)
! this subroutine performs the x part of a two dimensional fast real
! sine transform and its inverse, for a subset of y,
! using real arithmetic, with OpenMP,
! for data which is distributed in blocks
! algorithm is described in Numerical Recipies in Fortran, Second Ed.,
! by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
! [Cambridge Univ. Press, 1992], p. 508.
! for isign = (-1,1), input: all, output: f
! approximate flop count: N*(5*log2(N) + 18)/nvp
! where N = (nx/2)*ny
! indx/indy = exponent which determines length in x/y direction,
! where nx=2**indx, ny=2**indy
! if isign = -1, an inverse sine transform is performed
! f(n,k) = (1/nx*ny)*sum(f(j,k)*sin(pi*n*j/nx))
! if isign = 1, a forward sine transform is performed
! f(j,k) = sum(f(n,k)*sin(pi*n*j/nx))
! mixup = array of bit reversed addresses
! sctd = sine/cosine table
! kstrt = starting data block number
! kypi = initial y index used
! kypp = number of y indices used
! nxvh = first dimension of f >= nx/2 + 1
! kypd = second dimension of f >= kyp + 1
! nxhyd = maximum of (nx/2,ny)
! nxyd = maximum of (nx,ny)
! written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, kstrt, kypi, kypp, nxvh, kypd
      integer nxhyd, nxyd, mixup
      real f
      complex sctd
      dimension f(2*nxvh,kypd)
      dimension mixup(nxhyd), sctd(nxyd)
! local data
      integer indx1, indx1y, nx, nxh, nxhh, nx3, ny, nxy, nxhy, ks
      integer i, j, k, m, km, kmr, nrx, j1, j2, ns, ns2, k1, k2, kyps
      integer nrxb
      real at1, at2, t2, t3, t4, t5, t6, ani
      complex t1
      double precision sum1
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nx3 = nx + 3
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 1
      kyps = kypi + kypp - 1
      if (kstrt.gt.ny) return
      if (isign.eq.0) return
      ani = 0.5/(real(nx)*real(ny))
      nrxb = nxhy/nxh
      nrx = nxy/nxh
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,m,ns,ns2,km,kmr,k1,k2,j1,j2,at1,at2,t2,t3,t4,t5,t6,
!$OMP& t1,sum1)
      do 90 i = kypi, kyps
! create auxiliary array in x
      kmr = nxy/nx
      do 10 j = 2, nxh
      j1 = 1 + kmr*(j - 1)
      at2 = f(nx+2-j,i)
      at1 = f(j,i) + at2
      at2 = f(j,i) - at2
      at1 = -aimag(sctd(j1))*at1
      at2 = 0.5*at2
      f(j,i) = at1 + at2
      f(nx+2-j,i) = at1 - at2
   10 continue
      f(1,i) = 0.0
      f(nxh+1,i) = 2.0*f(nxh+1,i)
! bit-reverse array elements in x
      do 20 j = 1, nxh
      j1 = (mixup(j) - 1)/nrxb + 1
      if (j.lt.j1) then
         t2 = f(2*j1-1,i)
         t3 = f(2*j1,i)
         f(2*j1-1,i) = f(2*j-1,i)
         f(2*j1,i) = f(2*j,i)
         f(2*j-1,i) = t2
         f(2*j,i) = t3
      endif
   20 continue
! then transform in x
      do 50 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = 2*km*nrx
      do 40 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 30 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      t2 = real(t1)*f(2*j2-1,i) - aimag(t1)*f(2*j2,i)
      t3 = aimag(t1)*f(2*j2-1,i) + real(t1)*f(2*j2,i)
      f(2*j2-1,i) = f(2*j1-1,i) - t2
      f(2*j2,i) = f(2*j1,i) - t3
      f(2*j1-1,i) = f(2*j1-1,i) + t2
      f(2*j1,i) = f(2*j1,i) + t3
   30 continue
   40 continue
   50 continue
! unscramble coefficients and normalize
! inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nxh
         do 60 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         t4 = f(nx3-2*j,i)
         t5 = -f(nx3-2*j+1,i)
         t2 = f(2*j-1,i) + t4
         t3 = f(2*j,i) + t5
         t6 = f(2*j-1,i) - t4
         t5 = f(2*j,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(2*j-1,i) = ani*(t2 + t4)
         f(2*j,i) = ani*(t3 + t5)
         f(nx3-2*j,i) = ani*(t2 - t4)
         f(nx3-2*j+1,i) = ani*(t5 - t3)
   60    continue
         f(nxh+1,i) = 2.0*ani*f(nxh+1,i)
         f(nxh+2,i) = -2.0*ani*f(nxh+2,i)
         t2 = 2.0*ani*(f(1,i) + f(2,i))
         f(2,i) = 2.0*ani*(f(1,i) - f(2,i))
         f(1,i) = t2
         f(nx+1,i) = 2.0*ani*f(nx+1,i)
! forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nxh
         do 70 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         t4 = f(nx3-2*j,i)
         t5 = -f(nx3-2*j+1,i)
         t2 = f(2*j-1,i) + t4
         t3 = f(2*j,i) + t5
         t6 = f(2*j-1,i) - t4
         t5 = f(2*j,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(2*j-1,i) = t2 + t4
         f(2*j,i) = t3 + t5
         f(nx3-2*j,i) = t2 - t4
         f(nx3-2*j+1,i) = t5 - t3
   70    continue
         f(nxh+1,i) = 2.0*f(nxh+1,i)
         f(nxh+2,i) = -2.0*f(nxh+2,i)
         t2 = 2.0*(f(1,i) + f(2,i))
         f(2,i) = 2.0*(f(1,i) - f(2,i))
         f(1,i) = t2
         f(nx+1,i) = 2.0*f(nx+1,i)
      endif
! perform recursion for sine transform
      sum1 = 0.5*f(1,i)
      f(1,i) = 0.0
      f(2,i) = sum1
      do 80 j = 2, nxh
      sum1 = sum1 + f(2*j-1,i)
      f(2*j-1,i) = -f(2*j,i)
      f(2*j,i) = sum1
   80 continue
      f(nx+1,i) = 0.0
   90 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPFCT2RMXX(f,isign,mixup,sctd,indx,indy,kstrt,kypi,kypp&
     &,nxvh,kypd,nxhyd,nxyd)
! this subroutine performs the x part of a two dimensional fast real
! cosine transform and its inverse, for a subset of y,
! using real arithmetic, with OpenMP,
! for data which is distributed in blocks
! algorithm is described in Numerical Recipies in Fortran, Second Ed.,
! by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
! [Cambridge Univ. Press, 1992], p. 508.
! for isign = (-1,1), input: all, output: f
! approximate flop count: N*(5*log2(N) + 18)/nvp
! where N = (nx/2)*ny
! indx/indy = exponent which determines length in x/y direction,
! where nx=2**indx, ny=2**indy
! if isign = -1, an inverse cosine transform is performed
! f(n,k) = (1/nx*ny)*(.5*f(1,k) + ((-1)**n)*f(nx+1,k)
!            + sum(f(j,k)*cos(pi*n*j/nx)))
! if isign = 1, a forward cosine transform is performed
! f(j,k) = 2*(.5*f(1,k) + ((-1)**j)*f(n+1,k) + sum(f(n,k)*
!            cos(pi*n*j/nx))
! mixup = array of bit reversed addresses
! sctd = sine/cosine table
! kstrt = starting data block number
! kypi = initial y index used
! kypp = number of y indices used
! nxvh = first dimension of f >= nx/2 + 1
! kypd = second dimension of f >= kyp+1
! nxhyd = maximum of (nx/2,ny)
! nxyd = maximum of (nx,ny)
! written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, kstrt, kypi, kypp, nxvh, kypd
      integer nxhyd, nxyd, mixup
      real f
      complex sctd
      dimension f(2*nxvh,kypd)
      dimension mixup(nxhyd), sctd(nxyd)
! local data
      integer indx1, indx1y, nx, nxh, nxhh, nx3, ny, nxy, nxhy, ks
      integer i, j, k, m, km, kmr, nrx, j1, j2, ns, ns2, k1, k2, kyps
      integer nrxb
      real at1, at2, t2, t3, t4, t5, t6, ani
      double precision sum1
      complex t1
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nx3 = nx + 3
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 1
      kyps = kypi + kypp - 1
      if (kstrt.gt.ny) return
      if (isign.eq.0) return
      ani = 0.5/(real(nx)*real(ny))
      nrxb = nxhy/nxh
      nrx = nxy/nxh
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,m,ns,ns2,km,kmr,k1,k2,j1,j2,at1,at2,t2,t3,t4,t5,t6,
!$OMP& t1,sum1)
      do 90 i = kypi, kyps
! create auxiliary array in x
      kmr = nxy/nx
      sum1 = 0.5*(f(1,i) - f(nx+1,i))
      do 10 j = 2, nxh
      j1 = 1 + kmr*(j - 1)
      at2 = f(nx+2-j,i)
      at1 = f(j,i) + at2
      at2 = f(j,i) - at2
      sum1 = sum1 + real(sctd(j1))*at2
      at2 = -aimag(sctd(j1))*at2
      at1 = 0.5*at1
      f(j,i) = at1 - at2
      f(nx+2-j,i) = at1 + at2
   10 continue
      f(1,i) = 0.5*(f(1,i) + f(nx+1,i))
      f(nx+1,i) = sum1
! bit-reverse array elements in x
      do 20 j = 1, nxh
      j1 = (mixup(j) - 1)/nrxb + 1
      if (j.lt.j1) then
         t2 = f(2*j1-1,i)
         t3 = f(2*j1,i)
         f(2*j1-1,i) = f(2*j-1,i)
         f(2*j1,i) = f(2*j,i)
         f(2*j-1,i) = t2
         f(2*j,i) = t3
      endif
   20 continue
! then transform in x
      do 50 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = 2*km*nrx
      do 40 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 30 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      t2 = real(t1)*f(2*j2-1,i) - aimag(t1)*f(2*j2,i)
      t3 = aimag(t1)*f(2*j2-1,i) + real(t1)*f(2*j2,i)
      f(2*j2-1,i) = f(2*j1-1,i) - t2
      f(2*j2,i) = f(2*j1,i) - t3
      f(2*j1-1,i) = f(2*j1-1,i) + t2
      f(2*j1,i) = f(2*j1,i) + t3
   30 continue
   40 continue
   50 continue
! unscramble coefficients and normalize
! inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nxh
         do 60 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         t4 = f(nx3-2*j,i)
         t5 = -f(nx3-2*j+1,i)
         t2 = f(2*j-1,i) + t4
         t3 = f(2*j,i) + t5
         t6 = f(2*j-1,i) - t4
         t5 = f(2*j,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(2*j-1,i) = ani*(t2 + t4)
         f(2*j,i) = ani*(t3 + t5)
         f(nx3-2*j,i) = ani*(t2 - t4)
         f(nx3-2*j+1,i) = ani*(t5 - t3)
   60    continue
         f(nxh+1,i) = 2.0*ani*f(nxh+1,i)
         f(nxh+2,i) = -2.0*ani*f(nxh+2,i)
         t2 = 2.0*ani*(f(1,i) + f(2,i))
         f(2,i) = 2.0*ani*(f(1,i) - f(2,i))
         f(1,i) = t2
         f(nx+1,i) = 2.0*ani*f(nx+1,i)
! forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nxh
         do 70 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         t4 = f(nx3-2*j,i)
         t5 = -f(nx3-2*j+1,i)
         t2 = f(2*j-1,i) + t4
         t3 = f(2*j,i) + t5
         t6 = f(2*j-1,i) - t4
         t5 = f(2*j,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(2*j-1,i) = t2 + t4
         f(2*j,i) = t3 + t5
         f(nx3-2*j,i) = t2 - t4
         f(nx3-2*j+1,i) = t5 - t3
   70    continue
         f(nxh+1,i) = 2.0*f(nxh+1,i)
         f(nxh+2,i) = -2.0*f(nxh+2,i)
         t2 = 2.0*(f(1,i) + f(2,i))
         f(2,i) = 2.0*(f(1,i) - f(2,i))
         f(1,i) = t2
         f(nx+1,i) = 2.0*f(nx+1,i)
      endif
! perform recursion for cosine transform
      sum1 = f(nx+1,i)
      f(nx+1,i) = f(2,i)
      f(2,i) = sum1
      do 80 j = 2, nxh
      sum1 = sum1 - f(2*j,i)
      f(2*j,i) = sum1
   80 continue
   90 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPFST2RMXY(g,isign,mixup,sctd,indx,indy,kstrt,kxpi,kxpp&
     &,nyv,kxpd,nxhyd,nxyd)
! this subroutine performs the y part of a two dimensional fast real
! sine transform and its inverse, for a subset of x,
! using real arithmetic, with OpenMP,
! for data which is distributed in blocks
! algorithm is described in Numerical Recipies in Fortran, Second Ed.,
! by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
! [Cambridge Univ. Press, 1992], p. 508.
! for isign = (-1,1), input: all, output: f
! approximate flop count: N*(5*log2(N) + 18)/nvp
! where N = (nx/2)*ny
! indx/indy = exponent which determines length in x/y direction,
! where nx=2**indx, ny=2**indy
! if isign = -1, an inverse sine transform is performed
! g(m,n) = sum(g(k,n)*sin(pi*m*k/ny))
! if isign = 1, a forward sine transform is performed
! g(k,n) = sum(g(m,n)*sin(pi*m*k/ny))
! mixup = array of bit reversed addresses
! sctd = sine/cosine table
! kstrt = starting data block number
! kxpi = initial x index used
! kxpp = number of x indices used
! nyv = first dimension of g >= ny + 1
! kxpd = second dimension of g >= kxp + 1
! nxhyd = maximum of (nx/2,ny)
! nxyd = maximum of (nx,ny)
! written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, kstrt, kxpi, kxpp, nyv, kxpd
      integer nxhyd, nxyd, mixup
      real g
      complex sctd
      dimension g(nyv,kxpd)
      dimension mixup(nxhyd), sctd(nxyd)
! local data
      integer indx1, indy1, indx1y, nx, ny, nyh, nyhh, ny3, nxy, nxhy
      integer i, j, k, m, ks, km, kmr, nry, j1, j2, ns, ns2, k1, k2
      integer kxps, nryb
      real at1, at2, t2, t3, t4, t5, t6
      complex t1
      double precision sum1
      indx1 = indx - 1
      indy1 = indy - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      nyhh = ny/4
      ny3 = ny + 3
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 1
      kxps = kxpi + kxpp - 1
      if (kstrt.gt.nx) return
      if (isign.eq.0) return
      nryb = nxhy/nyh
      nry = nxy/nyh
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,m,ns,ns2,km,kmr,k1,k2,j1,j2,at1,at2,t2,t3,t4,t5,t6,
!$OMP& t1,sum1)
      do 90 i = kxpi, kxps
! create auxiliary array in y
      kmr = nxy/ny
      do 10 k = 2, nyh
      k1 = 1 + kmr*(k - 1)
      at2 = g(ny+2-k,i)
      at1 = g(k,i) + at2
      at2 = g(k,i) - at2
      at1 = -aimag(sctd(k1))*at1
      at2 = 0.5*at2
      g(k,i) = at1 + at2
      g(ny+2-k,i) = at1 - at2
   10 continue
      g(1,i) = 0.0
      g(nyh+1,i) = 2.0*g(nyh+1,i)
! bit-reverse array elements in y
      do 20 k = 1, nyh
      k1 = (mixup(k) - 1)/nryb + 1
      if (k.lt.k1) then
         t2 = g(2*k1-1,i)
         t3 = g(2*k1,i)
         g(2*k1-1,i) = g(2*k-1,i)
         g(2*k1,i) = g(2*k,i)
         g(2*k-1,i) = t2
         g(2*k,i) = t3
      endif
   20 continue
! then transform in y
      do 50 m = 1, indy1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyhh/ns
      kmr = 2*km*nry
      do 40 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 30 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      t2 = real(t1)*g(2*j2-1,i) - aimag(t1)*g(2*j2,i)
      t3 = aimag(t1)*g(2*j2-1,i) + real(t1)*g(2*j2,i)
      g(2*j2-1,i) = g(2*j1-1,i) - t2
      g(2*j2,i) = g(2*j1,i) - t3
      g(2*j1-1,i) = g(2*j1-1,i) + t2
      g(2*j1,i) = g(2*j1,i) + t3
   30 continue
   40 continue
   50 continue
! unscramble coefficients and normalize
! inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nyh
         do 60 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         t4 = g(ny3-2*k,i)
         t5 = -g(ny3-2*k+1,i)
         t2 = g(2*k-1,i) + t4
         t3 = g(2*k,i) + t5
         t6 = g(2*k-1,i) - t4
         t5 = g(2*k,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         g(2*k-1,i) = 0.5*(t2 + t4)
         g(2*k,i) = 0.5*(t3 + t5)
         g(ny3-2*k,i) = 0.5*(t2 - t4)
         g(ny3-2*k+1,i) = 0.5*(t5 - t3)
   60    continue
         g(nyh+1,i) = g(nyh+1,i)
         g(nyh+2,i) = -g(nyh+2,i)
         t2 = g(1,i) + g(2,i)
         g(2,i) = g(1,i) - g(2,i)
         g(1,i) = t2
         g(ny+1,i) = g(ny+1,i)
! forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nyh
         do 70 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         t4 = g(ny3-2*k,i)
         t5 = -g(ny3-2*k+1,i)
         t2 = g(2*k-1,i) + t4
         t3 = g(2*k,i) + t5
         t6 = g(2*k-1,i) - t4
         t5 = g(2*k,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         g(2*k-1,i) = t2 + t4
         g(2*k,i) = t3 + t5
         g(ny3-2*k,i) = t2 - t4
         g(ny3-2*k+1,i) = t5 - t3
   70    continue
         g(nyh+1,i) = 2.0*g(nyh+1,i)
         g(nyh+2,i) = -2.0*g(nyh+2,i)
         t2 = 2.0*(g(1,i) + g(2,i))
         g(2,i) = 2.0*(g(1,i) - g(2,i))
         g(1,i) = t2
         g(ny+1,i) = 2.0*g(ny+1,i)
      endif
! perform recursion for sine transform
      sum1 = 0.5*g(1,i)
      g(1,i) = 0.0
      g(2,i) = sum1
      do 80 k = 2, nyh
      sum1 = sum1 + g(2*k-1,i)
      g(2*k-1,i) = -g(2*k,i)
      g(2*k,i) = sum1
   80 continue
      g(ny+1,i) = 0.0
   90 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPFCT2RMXY(g,isign,mixup,sctd,indx,indy,kstrt,kxpi,kxpp&
     &,nyv,kxpd,nxhyd,nxyd)
! this subroutine performs the y part of a two dimensional fast real
! cosine transform and its inverse, for a subset of x,
! using real arithmetic, with OpenMP,
! for data which is distributed in blocks
! algorithm is described in Numerical Recipies in Fortran, Second Ed.,
! by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
! [Cambridge Univ. Press, 1992], p. 508.
! for isign = (-1,1), input: all, output: f
! approximate flop count: N*(5*log2(N) + 18)/nvp
! where N = (nx/2)*ny
! indx/indy = exponent which determines length in x/y direction,
! where nx=2**indx, ny=2**indy
! if isign = -1, an inverse cosine transform is performed
! g(m,n) = (.5*g(1,n) + ((-1)**m)*g(ny+1,n)
!            + sum(g(k,n)*cos(pi*m*k/ny))
! if isign = 1, a forward cosine transform is performed
! g(k,n) = 2*(.5*g(1,n) + ((-1)**m)*g(ny+1,n) + sum(g(m,n)*
!            cos(pi*m*k/ny))
! mixup = array of bit reversed addresses
! sctd = sine/cosine table
! kstrt = starting data block number
! kxpi = initial x index used
! kxpp = number of x indices used
! nyv = first dimension of g >= ny + 1
! kxpd = second dimension of g >= kxp + 1
! nxhyd = maximum of (nx/2,ny)
! nxyd = maximum of (nx,ny)
! written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, kstrt, kxpi, kxpp, nyv, kxpd
      integer nxhyd, nxyd, mixup
      real g
      complex sctd
      dimension g(nyv,kxpd)
      dimension mixup(nxhyd), sctd(nxyd)
! local data
      integer indx1, indy1, indx1y, nx, ny, nyh, nyhh, ny3, nxy, nxhy
      integer i, j, k, m, ks, km, kmr, nry, j1, j2, ns, ns2, k1, k2
      integer kxps, nryb
      real at1, at2, t2, t3, t4, t5, t6
      complex t1
      double precision sum1
      indx1 = indx - 1
      indy1 = indy - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      nyhh = ny/4
      ny3 = ny + 3
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 1
      kxps = kxpi + kxpp - 1
      if (kstrt.gt.nx) return
      if (isign.eq.0) return
      nryb = nxhy/nyh
      nry = nxy/nyh
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,m,ns,ns2,km,kmr,k1,k2,j1,j2,at1,at2,t2,t3,t4,t5,t6,
!$OMP& t1,sum1)
      do 90 i = kxpi, kxps
! create auxiliary array in y
      kmr = nxy/ny
      sum1 = 0.5*(g(1,i) - g(ny+1,i))
      do 10 k = 2, nyh
      k1 = 1 + kmr*(k - 1)
      at2 = g(ny+2-k,i)
      at1 = g(k,i) + at2
      at2 = g(k,i) - at2
      sum1 = sum1 + real(sctd(k1))*at2
      at2 = -aimag(sctd(k1))*at2
      at1 = 0.5*at1
      g(k,i) = at1 - at2
      g(ny+2-k,i) = at1 + at2
   10 continue
      g(1,i) = 0.5*(g(1,i) + g(ny+1,i))
      g(ny+1,i) = sum1
! bit-reverse array elements in y
      do 20 k = 1, nyh
      k1 = (mixup(k) - 1)/nryb + 1
      if (k.lt.k1) then
         t2 = g(2*k1-1,i)
         t3 = g(2*k1,i)
         g(2*k1-1,i) = g(2*k-1,i)
         g(2*k1,i) = g(2*k,i)
         g(2*k-1,i) = t2
         g(2*k,i) = t3
      endif
   20 continue
! then transform in y
      do 50 m = 1, indy1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyhh/ns
      kmr = 2*km*nry
      do 40 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 30 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      t2 = real(t1)*g(2*j2-1,i) - aimag(t1)*g(2*j2,i)
      t3 = aimag(t1)*g(2*j2-1,i) + real(t1)*g(2*j2,i)
      g(2*j2-1,i) = g(2*j1-1,i) - t2
      g(2*j2,i) = g(2*j1,i) - t3
      g(2*j1-1,i) = g(2*j1-1,i) + t2
      g(2*j1,i) = g(2*j1,i) + t3
   30 continue
   40 continue
   50 continue
! unscramble coefficients and normalize
! inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nyh
         do 60 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         t4 = g(ny3-2*k,i)
         t5 = -g(ny3-2*k+1,i)
         t2 = g(2*k-1,i) + t4
         t3 = g(2*k,i) + t5
         t6 = g(2*k-1,i) - t4
         t5 = g(2*k,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         g(2*k-1,i) = 0.5*(t2 + t4)
         g(2*k,i) = 0.5*(t3 + t5)
         g(ny3-2*k,i) = 0.5*(t2 - t4)
         g(ny3-2*k+1,i) = 0.5*(t5 - t3)
   60    continue
         g(nyh+1,i) = g(nyh+1,i)
         g(nyh+2,i) = -g(nyh+2,i)
         t2 = g(1,i) + g(2,i)
         g(2,i) = g(1,i) - g(2,i)
         g(1,i) = t2
         g(ny+1,i) = g(ny+1,i)
! forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nyh
         do 70 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         t4 = g(ny3-2*k,i)
         t5 = -g(ny3-2*k+1,i)
         t2 = g(2*k-1,i) + t4
         t3 = g(2*k,i) + t5
         t6 = g(2*k-1,i) - t4
         t5 = g(2*k,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         g(2*k-1,i) = t2 + t4
         g(2*k,i) = t3 + t5
         g(ny3-2*k,i) = t2 - t4
         g(ny3-2*k+1,i) = t5 - t3
   70    continue
         g(nyh+1,i) = 2.0*g(nyh+1,i)
         g(nyh+2,i) = -2.0*g(nyh+2,i)
         t2 = 2.0*(g(1,i) + g(2,i))
         g(2,i) = 2.0*(g(1,i) - g(2,i))
         g(1,i) = t2
         g(ny+1,i) = 2.0*g(ny+1,i)
      endif
! perform recursion for cosine transform
      sum1 = g(ny+1,i)
      g(ny+1,i) = g(2,i)
      g(2,i) = sum1
      do 80 k = 2, nyh
      sum1 = sum1 - g(2*k,i)
      g(2*k,i) = sum1
   80 continue
   90 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine WPPFCST2RM2(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,&
     &indy,kstrt,nvp,nxvh,nyv,kxp2,kyp,kypd,kxp2d,nxhyd,nxyd)
! wrapper function for 2 parallel real cosine/sine transforms
! for the electric field with dirichlet or magnetic field with neumann
! boundary conditions
! x component has a cosine/sine transform in x and y, respectively
! y component has a sine/cosine transform in x and y, respectively
      implicit none
      integer isign, ntpose, indx, indy, kstrt, nvp, nxvh, nyv
      integer kxp2, kyp, kypd, kxp2d, nxhyd, nxyd, mixup
      real f, g, bs, br, ttp
      complex sctd
      dimension f(2,2*nxvh,kypd), g(2,nyv,kxp2d)
      dimension bs(2,kxp2+1,kyp+1), br(2,kxp2+1,kyp+1)
      dimension mixup(nxhyd), sctd(nxyd)
! local data
      integer nx, ny, kxpi, kypi, ks, kxp2p, kypp, kxb2, kyb
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
! calculate range of indices
      nx = 2**indx
      ny = 2**indy
! ks = processor id
      ks = kstrt - 1
! kxp2p = actual size used in x direction
      kxp2p = min(kxp2,max(0,nx-kxp2*ks))
! kypp = actual size used in y direction
      kypp = min(kyp,max(0,ny-kyp*ks))
! kxb2 = minimum number of processors needed in x direction
      kxb2 = (nx - 1)/kxp2 + 1
! kyb = minimum number of processors needed in y direction
      kyb = (ny - 1)/kyp + 1
! add extra word for last processor in x
      if (ks==(kxb2-1)) kxp2p = kxp2p + 1
! add extra word for last processor in y
      if (ks==(kyb-1)) kypp = kypp + 1
! inverse fourier transform
      if (isign.lt.0) then
! perform x cosine-sine transform
         call PPFCST2RM2X(f,isign,mixup,sctd,indx,indy,kstrt,kypi,kypp, &
     &nxvh,kypd,nxhyd,nxyd)
! transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PPRNTPOSE(f,g,bs,br,nx,ny,kxp2,kyp,kstrt,nvp,2,2*nxvh,nyv,&
     &kxp2d,kypd)
         call PWTIMERA(1,ttp,dtime)
! perform y sine-cosine transform
         call PPFSCT2RM2Y(g,isign,mixup,sctd,indx,indy,kstrt,kxpi,kxp2p,&
     &nyv,kxp2d,nxhyd,nxyd)
! transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPRNTPOSE(g,f,br,bs,ny,nx,kyp,kxp2,kstrt,nvp,2,nyv,    &
     &2*nxvh,kypd,kxp2d)
            call PWTIMERA(1,tf,dtime)
         endif
! forward fourier transform
      else if (isign.gt.0) then
! transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPRNTPOSE(f,g,bs,br,nx,ny,kxp2,kyp,kstrt,nvp,2,2*nxvh, &
     &nyv,kxp2d,kypd)
            call PWTIMERA(1,tf,dtime)
         endif
! perform y sine-cosine transform
         call PPFSCT2RM2Y(g,isign,mixup,sctd,indx,indy,kstrt,kxpi,kxp2p,&
     &nyv,kxp2d,nxhyd,nxyd)
! transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PPRNTPOSE(g,f,br,bs,ny,nx,kyp,kxp2,kstrt,nvp,2,nyv,2*nxvh,&
     &kypd,kxp2d)
         call PWTIMERA(1,ttp,dtime)
! perform x cosine-sine transform
         call PPFCST2RM2X(f,isign,mixup,sctd,indx,indy,kstrt,kypi,kypp, &
     &nxvh,kypd,nxhyd,nxyd)
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
!-----------------------------------------------------------------------
      subroutine WPPFSCT2RM2(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,&
     &indy,kstrt,nvp,nxvh,nyv,kxp2,kyp,kypd,kxp2d,nxhyd,nxyd)
! wrapper function for 2 parallel real sine/cosine transforms
! for the magnetic field with dirichlet or electric field with neumann
! boundary conditions
! x component has a sine/cosine transform in x and y, respectively
! y component has a cosine/sine transform in x and y, respectively
      implicit none
      integer isign, ntpose, indx, indy, kstrt, nvp, nxvh, nyv
      integer kxp2, kyp, kypd, kxp2d, nxhyd, nxyd, mixup
      real f, g, bs, br, ttp
      complex sctd
      dimension f(2,2*nxvh,kypd), g(2,nyv,kxp2d)
      dimension bs(2,kxp2+1,kyp+1), br(2,kxp2+1,kyp+1)
      dimension mixup(nxhyd), sctd(nxyd)
! local data
      integer nx, ny, kxpi, kypi, ks, kxp2p, kypp, kxb2, kyb
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
! calculate range of indices
      nx = 2**indx
      ny = 2**indy
! ks = processor id
      ks = kstrt - 1
! kxp2p = actual size used in x direction
      kxp2p = min(kxp2,max(0,nx-kxp2*ks))
! kypp = actual size used in y direction
      kypp = min(kyp,max(0,ny-kyp*ks))
! kxb2 = minimum number of processors needed in x direction
      kxb2 = (nx - 1)/kxp2 + 1
! kyb = minimum number of processors needed in y direction
      kyb = (ny - 1)/kyp + 1
! add extra word for last processor in x
      if (ks==(kxb2-1)) kxp2p = kxp2p + 1
! add extra word for last processor in y
      if (ks==(kyb-1)) kypp = kypp + 1
! inverse fourier transform
      if (isign.lt.0) then
! perform x sine-cosine transform
         call PPFSCT2RM2X(f,isign,mixup,sctd,indx,indy,kstrt,kypi,kypp, &
     &nxvh,kypd,nxhyd,nxyd)
! transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PPRNTPOSE(f,g,bs,br,nx,ny,kxp2,kyp,kstrt,nvp,2,2*nxvh,nyv,&
     &kxp2d,kypd)
         call PWTIMERA(1,ttp,dtime)
! perform y cosine-sine transform
         call PPFCST2RM2Y(g,isign,mixup,sctd,indx,indy,kstrt,kxpi,kxp2p,&
     &nyv,kxp2d,nxhyd,nxyd)
! transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPRNTPOSE(g,f,br,bs,ny,nx,kyp,kxp2,kstrt,nvp,2,nyv,    &
     &2*nxvh,kypd,kxp2d)
            call PWTIMERA(1,tf,dtime)
         endif
! forward fourier transform
      else if (isign.gt.0) then
! transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPRNTPOSE(f,g,bs,br,nx,ny,kxp2,kyp,kstrt,nvp,2,2*nxvh, &
     &nyv,kxp2d,kypd)
            call PWTIMERA(1,tf,dtime)
         endif
! perform y cosine-sine transform
         call PPFCST2RM2Y(g,isign,mixup,sctd,indx,indy,kstrt,kxpi,kxp2p,&
     &nyv,kxp2d,nxhyd,nxyd)
! transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PPRNTPOSE(g,f,br,bs,ny,nx,kyp,kxp2,kstrt,nvp,2,nyv,2*nxvh,&
     &kypd,kxp2d)
         call PWTIMERA(1,ttp,dtime)
! perform x sine-cosine transform
         call PPFSCT2RM2X(f,isign,mixup,sctd,indx,indy,kstrt,kypi,kypp, &
     &nxvh,kypd,nxhyd,nxyd)
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
!-----------------------------------------------------------------------
      subroutine PPFCST2RM2X(f,isign,mixup,sctd,indx,indy,kstrt,kypi,   &
     &kypp,nxvh,kypd,nxhyd,nxyd)
! this subroutine performs the x part of 2 two dimensional fast real
! sine and cosine transforms and their inverses, for a subset of y,
! using real arithmetic, with OpenMP,
! for data which is distributed in blocks
! x component has a cosine transform, y component a sine transform
! algorithm is described in Numerical Recipies in Fortran, Second Ed.,
! by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
! [Cambridge Univ. Press, 1992], p. 508.
! for isign = (-1,1), input: all, output: f
! approximate flop count: N*(5*log2(N) + 18)/nvp
! where N = (nx/2)*ny
! indx/indy = exponent which determines length in x/y direction,
! where nx=2**indx, ny=2**indy
! if isign = -1, inverse cosine-sine transforms are performed
! f(1,n,k) = (1/nx*ny)*(.5*f(1,1,k) + ((-1)**n)*f(1,nx+1,k)
!              + sum(f(1,j,k)*cos(pi*n*j/nx)))
! f(2,n,k) = (1/nx*ny)*sum(f(2,j,k)*sin(pi*n*j/nx))
! if isign = 1, forward cosine-sine transforms are performed
! f(1,j,k) = 2*(.5*f(1,1,k) + ((-1)**j)*f(1,n+1,k)
!              + sum(f(1,n,k)*cos(pi*n*j/nx))
! f(2,j,k) = sum(f(2,n,k)*sin(pi*n*j/nx))
! mixup = array of bit reversed addresses
! sctd = sine/cosine table
! kstrt = starting data block number
! kypi = initial y index used
! kypp = number of y indices used
! nxvh = second dimension of f >= nx/2 + 1
! kypd = third dimension of f >= kyp + 1
! nxhyd = maximum of (nx/2,ny)
! nxyd = maximum of (nx,ny)
! written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, kstrt, kypi, kypp
      integer nxvh, kypd, nxhyd, nxyd, mixup
      real f
      complex sctd
      dimension f(2,2*nxvh,kypd)
      dimension mixup(nxhyd), sctd(nxyd)
! local data
      integer indx1, indx1y, nx, nxh, nxhh, nx3, ny, nxy, nxhy, ks
      integer i, j, k, m, km, kmr, nrx, j1, j2, ns, ns2, k1, k2, kyps
      integer nrxb, jj
      real at1, at2, at3, t2, t3, t4, t5, t6, ani
      complex t1
      double precision sum1, sum2
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nx3 = nx + 3
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 1
      kyps = kypi + kypp - 1
      if (kstrt.gt.ny) return
      if (isign.eq.0) return
      ani = 0.5/(real(nx)*real(ny))
      nrxb = nxhy/nxh
      nrx = nxy/nxh
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,m,jj,ns,ns2,km,kmr,k1,k2,j1,j2,at1,at2,at3,t2,t3,t4
!$OMP& ,t5,t6,t1,sum1,sum2)
      do 150 i = kypi, kyps
! create auxiliary array in x
      kmr = nxy/nx
      sum1 = 0.5*(f(1,1,i) - f(1,nx+1,i))
      do 10 j = 2, nxh
      j1 = 1 + kmr*(j - 1)
      at3 = -aimag(sctd(j1))
      at2 = f(1,nx+2-j,i)
      at1 = f(1,j,i) + at2
      at2 = f(1,j,i) - at2
      sum1 = sum1 + real(sctd(j1))*at2
      at2 = at3*at2
      at1 = 0.5*at1
      f(1,j,i) = at1 - at2
      f(1,nx+2-j,i) = at1 + at2
      at2 = f(2,nx+2-j,i)
      at1 = f(2,j,i) + at2
      at2 = f(2,j,i) - at2
      at1 = at3*at1
      at2 = 0.5*at2
      f(2,j,i) = at1 + at2
      f(2,nx+2-j,i) = at1 - at2
   10 continue
      f(1,1,i) = 0.5*(f(1,1,i) + f(1,nx+1,i))
      f(1,nx+1,i) = sum1
      f(2,1,i) = 0.0
      f(2,nxh+1,i) = 2.0*f(2,nxh+1,i)
! bit-reverse array elements in x
      do 30 j = 1, nxh
      j1 = (mixup(j) - 1)/nrxb + 1
      if (j.lt.j1) then
         do 20 jj = 1, 2
         t2 = f(jj,2*j1-1,i)
         t3 = f(jj,2*j1,i)
         f(jj,2*j1-1,i) = f(jj,2*j-1,i)
         f(jj,2*j1,i) = f(jj,2*j,i)
         f(jj,2*j-1,i) = t2
         f(jj,2*j,i) = t3
   20    continue
      endif
   30 continue
! then transform in x
      do 70 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = 2*km*nrx
      do 60 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 50 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 40 jj = 1, 2
      t2 = real(t1)*f(jj,2*j2-1,i) - aimag(t1)*f(jj,2*j2,i)
      t3 = aimag(t1)*f(jj,2*j2-1,i) + real(t1)*f(jj,2*j2,i)
      f(jj,2*j2-1,i) = f(jj,2*j1-1,i) - t2
      f(jj,2*j2,i) = f(jj,2*j1,i) - t3
      f(jj,2*j1-1,i) = f(jj,2*j1-1,i) + t2
      f(jj,2*j1,i) = f(jj,2*j1,i) + t3
   40 continue
   50 continue
   60 continue
   70 continue
! unscramble coefficients and normalize
! inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nxh
         do 90 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 80 jj = 1, 2
         t4 = f(jj,nx3-2*j,i)
         t5 = -f(jj,nx3-2*j+1,i)
         t2 = f(jj,2*j-1,i) + t4
         t3 = f(jj,2*j,i) + t5
         t6 = f(jj,2*j-1,i) - t4
         t5 = f(jj,2*j,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,2*j-1,i) = ani*(t2 + t4)
         f(jj,2*j,i) = ani*(t3 + t5)
         f(jj,nx3-2*j,i) = ani*(t2 - t4)
         f(jj,nx3-2*j+1,i) = ani*(t5 - t3)
   80    continue
   90    continue
         do 100 jj = 1, 2
         f(jj,nxh+1,i) = 2.0*ani*f(jj,nxh+1,i)
         f(jj,nxh+2,i) = -2.0*ani*f(jj,nxh+2,i)
         t2 = 2.0*ani*(f(jj,1,i) + f(jj,2,i))
         f(jj,2,i) = 2.0*ani*(f(jj,1,i) - f(jj,2,i))
         f(jj,1,i) = t2
         f(jj,nx+1,i) = 2.0*ani*f(jj,nx+1,i)
  100    continue
! forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nxh
         do 120 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 110 jj = 1, 2
         t4 = f(jj,nx3-2*j,i)
         t5 = -f(jj,nx3-2*j+1,i)
         t2 = f(jj,2*j-1,i) + t4
         t3 = f(jj,2*j,i) + t5
         t6 = f(jj,2*j-1,i) - t4
         t5 = f(jj,2*j,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,2*j-1,i) = t2 + t4
         f(jj,2*j,i) = t3 + t5
         f(jj,nx3-2*j,i) = t2 - t4
         f(jj,nx3-2*j+1,i) = t5 - t3
  110    continue
  120    continue
         do 130 jj = 1, 2
         f(jj,nxh+1,i) = 2.0*f(jj,nxh+1,i)
         f(jj,nxh+2,i) = -2.0*f(jj,nxh+2,i)
         t2 = 2.0*(f(jj,1,i) + f(jj,2,i))
         f(jj,2,i) = 2.0*(f(jj,1,i) - f(jj,2,i))
         f(jj,1,i) = t2
         f(jj,nx+1,i) = 2.0*f(jj,nx+1,i)
  130    continue
      endif
! perform recursion for cosine-sine transform
      sum1 = f(1,nx+1,i)
      f(1,nx+1,i) = f(1,2,i)
      f(1,2,i) = sum1
      sum2 = 0.5*f(2,1,i)
      f(2,1,i) = 0.0
      f(2,2,i) = sum2
      do 140 j = 2, nxh
      sum1 = sum1 - f(1,2*j,i)
      f(1,2*j,i) = sum1
      sum2 = sum2 + f(2,2*j-1,i)
      f(2,2*j-1,i) = -f(2,2*j,i)
      f(2,2*j,i) = sum2
  140 continue
      f(2,nx+1,i) = 0.0
  150 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPFSCT2RM2X(f,isign,mixup,sctd,indx,indy,kstrt,kypi,   &
     &kypp,nxvh,kypd,nxhyd,nxyd)
! this subroutine performs the x part of 2 two dimensional fast real
! sine and cosine transforms and their inverses, for a subset of y,
! using real arithmetic, with OpenMP,
! for data which is distributed in blocks
! x component has a sine transform, y component a cosine transform
! algorithm is described in Numerical Recipies in Fortran, Second Ed.,
! by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
! [Cambridge Univ. Press, 1992], p. 508.
! for isign = (-1,1), input: all, output: f
! approximate flop count: N*(5*log2(N) + 18)/nvp
! where N = (nx/2)*ny
! indx/indy = exponent which determines length in x/y direction,
! where nx=2**indx, ny=2**indy
! if isign = -1, inverse sine-cosine transforms are performed
! f(1,n,k) = (1/nx*ny)*sum(f(1,j,k)*sin(pi*n*j/nx))
! f(2,n,k) = (1/nx*ny)*(.5*f(2,1,k) + ((-1)**n)*f(2,nx+1,k)
!              + sum(f(2,j,k)*cos(pi*n*j/nx)))
! if isign = 1, forward sine-cosine transforms are performed
! f(1,j,k) = sum(f(1,n,k)*sin(pi*n*j/nx))
! f(2,j,k) = 2*(.5*f(2,1,k) + ((-1)**j)*f(2,n+1,k)
!              + sum(f(2,n,k)*cos(pi*n*j/nx))
! mixup = array of bit reversed addresses
! sctd = sine/cosine table
! kstrt = starting data block number
! kypi = initial y index used
! kypp = number of y indices used
! nxvh = second dimension of f >= nx/2 + 1
! kypd = third dimension of f >= kyp + 1
! nxhyd = maximum of (nx/2,ny)
! nxyd = maximum of (nx,ny)
! written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, kstrt, kypi, kypp
      integer nxvh, kypd, nxhyd, nxyd, mixup
      real f
      complex sctd
      dimension f(2,2*nxvh,kypd)
      dimension mixup(nxhyd), sctd(nxyd)
! local data
      integer indx1, indx1y, nx, nxh, nxhh, nx3, ny, nxy, nxhy, ks
      integer i, j, k, m, km, kmr, nrx, j1, j2, ns, ns2, k1, k2, kyps
      integer nrxb, jj
      real at1, at2, at3, t2, t3, t4, t5, t6, ani
      complex t1
      double precision sum1, sum2
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nx3 = nx + 3
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 1
      kyps = kypi + kypp - 1
      if (kstrt.gt.ny) return
      if (isign.eq.0) return
      ani = 0.5/(real(nx)*real(ny))
      nrxb = nxhy/nxh
      nrx = nxy/nxh
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,m,jj,ns,ns2,km,kmr,k1,k2,j1,j2,at1,at2,at3,t2,t3,t4
!$OMP& ,t5,t6,t1,sum1,sum2)
      do 150 i = kypi, kyps
! create auxiliary array in x
      kmr = nxy/nx
      sum1 = 0.5*(f(2,1,i) - f(2,nx+1,i))
      do 10 j = 2, nxh
      j1 = 1 + kmr*(j - 1)
      at3 = -aimag(sctd(j1))
      at2 = f(1,nx+2-j,i)
      at1 = f(1,j,i) + at2
      at2 = f(1,j,i) - at2
      at1 = at3*at1
      at2 = 0.5*at2
      f(1,j,i) = at1 + at2
      f(1,nx+2-j,i) = at1 - at2
      at2 = f(2,nx+2-j,i)
      at1 = f(2,j,i) + at2
      at2 = f(2,j,i) - at2
      sum1 = sum1 + real(sctd(j1))*at2
      at2 = at3*at2
      at1 = 0.5*at1
      f(2,j,i) = at1 - at2
      f(2,nx+2-j,i) = at1 + at2
   10 continue
      f(1,1,i) = 0.0
      f(1,nxh+1,i) = 2.0*f(1,nxh+1,i)
      f(2,1,i) = 0.5*(f(2,1,i) + f(2,nx+1,i))
      f(2,nx+1,i) = sum1
! bit-reverse array elements in x
      do 30 j = 1, nxh
      j1 = (mixup(j) - 1)/nrxb + 1
      if (j.lt.j1) then
         do 20 jj = 1, 2
         t2 = f(jj,2*j1-1,i)
         t3 = f(jj,2*j1,i)
         f(jj,2*j1-1,i) = f(jj,2*j-1,i)
         f(jj,2*j1,i) = f(jj,2*j,i)
         f(jj,2*j-1,i) = t2
         f(jj,2*j,i) = t3
   20    continue
      endif
   30 continue
! then transform in x
      do 70 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = 2*km*nrx
      do 60 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 50 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 40 jj = 1, 2
      t2 = real(t1)*f(jj,2*j2-1,i) - aimag(t1)*f(jj,2*j2,i)
      t3 = aimag(t1)*f(jj,2*j2-1,i) + real(t1)*f(jj,2*j2,i)
      f(jj,2*j2-1,i) = f(jj,2*j1-1,i) - t2
      f(jj,2*j2,i) = f(jj,2*j1,i) - t3
      f(jj,2*j1-1,i) = f(jj,2*j1-1,i) + t2
      f(jj,2*j1,i) = f(jj,2*j1,i) + t3
   40 continue
   50 continue
   60 continue
   70 continue
! unscramble coefficients and normalize
! inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nxh
         do 90 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 80 jj = 1, 2
         t4 = f(jj,nx3-2*j,i)
         t5 = -f(jj,nx3-2*j+1,i)
         t2 = f(jj,2*j-1,i) + t4
         t3 = f(jj,2*j,i) + t5
         t6 = f(jj,2*j-1,i) - t4
         t5 = f(jj,2*j,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,2*j-1,i) = ani*(t2 + t4)
         f(jj,2*j,i) = ani*(t3 + t5)
         f(jj,nx3-2*j,i) = ani*(t2 - t4)
         f(jj,nx3-2*j+1,i) = ani*(t5 - t3)
   80    continue
   90    continue
         do 100 jj = 1, 2
         f(jj,nxh+1,i) = 2.0*ani*f(jj,nxh+1,i)
         f(jj,nxh+2,i) = -2.0*ani*f(jj,nxh+2,i)
         t2 = 2.0*ani*(f(jj,1,i) + f(jj,2,i))
         f(jj,2,i) = 2.0*ani*(f(jj,1,i) - f(jj,2,i))
         f(jj,1,i) = t2
         f(jj,nx+1,i) = 2.0*ani*f(jj,nx+1,i)
  100    continue
! forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nxh
         do 120 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 110 jj = 1, 2
         t4 = f(jj,nx3-2*j,i)
         t5 = -f(jj,nx3-2*j+1,i)
         t2 = f(jj,2*j-1,i) + t4
         t3 = f(jj,2*j,i) + t5
         t6 = f(jj,2*j-1,i) - t4
         t5 = f(jj,2*j,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,2*j-1,i) = t2 + t4
         f(jj,2*j,i) = t3 + t5
         f(jj,nx3-2*j,i) = t2 - t4
         f(jj,nx3-2*j+1,i) = t5 - t3
  110    continue
  120    continue
         do 130 jj = 1, 2
         f(jj,nxh+1,i) = 2.0*f(jj,nxh+1,i)
         f(jj,nxh+2,i) = -2.0*f(jj,nxh+2,i)
         t2 = 2.0*(f(jj,1,i) + f(jj,2,i))
         f(jj,2,i) = 2.0*(f(jj,1,i) - f(jj,2,i))
         f(jj,1,i) = t2
         f(jj,nx+1,i) = 2.0*f(jj,nx+1,i)
  130    continue
      endif
! perform recursion for cosine-sine transform
      sum1 = 0.5*f(1,1,i)
      f(1,1,i) = 0.0
      f(1,2,i) = sum1
      sum2 = f(2,nx+1,i)
      f(2,nx+1,i) = f(2,2,i)
      f(2,2,i) = sum2
      do 140 j = 2, nxh
      sum1 = sum1 + f(1,2*j-1,i)
      f(1,2*j-1,i) = -f(1,2*j,i)
      f(1,2*j,i) = sum1
      sum2 = sum2 - f(2,2*j,i)
      f(2,2*j,i) = sum2
  140 continue
      f(1,nx+1,i) = 0.0
  150 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPFSCT2RM2Y(g,isign,mixup,sctd,indx,indy,kstrt,kxpi,   &
     &kxpp,nyv,kxpd,nxhyd,nxyd)
! this subroutine performs the y part of 2 two dimensional fast real
! sine and cosine transforms and their inverses, for a subset of x,
! using real arithmetic, with OpenMP,
! for data which is distributed in blocks
! x component has a sine transform, y component a cosine transform
! algorithm is described in Numerical Recipies in Fortran, Second Ed.,
! by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
! [Cambridge Univ. Press, 1992], p. 508.
! for isign = (-1,1), input: all, output: f
! approximate flop count: N*(5*log2(N) + 18)/nvp
! where N = (nx/2)*ny
! indx/indy = exponent which determines length in x/y direction,
! where nx=2**indx, ny=2**indy
! if isign = -1, inverse sine-cosine transform are performed
! g(1,m,n) = sum(g(1,k,n)*sin(pi*m*k/ny))
! g(2,m,n) = (.5*g(2,1,n) + ((-1)**m)*g(2,ny+1,n)
!              + sum(g(2,k,n)*cos(pi*m*k/ny))
! if isign = 1, a forward sine-cosine transforms are performed
! g(1,k,n) = sum(g(1,m,n)*sin(pi*m*k/ny))
! g(2,k,n) = 2*(.5*g(2,1,n) + ((-1)**m)*g(2,ny+1,n)
!              + sum(g(2,m,n)*cos(pi*m*k/ny))
! mixup = array of bit reversed addresses
! sctd = sine/cosine table
! kstrt = starting data block number
! kxpi = initial x index used
! kxpp = number of x indices used
! nyv = second dimension of g >= ny + 1
! kxpd = third dimension of g >= kxp + 1
! nxhyd = maximum of (nx/2,ny)
! nxyd = maximum of (nx,ny)
! written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, kstrt, kxpi, kxpp
      integer nyv, kxpd, nxhyd, nxyd, mixup
      real g
      complex sctd
      dimension g(2,nyv,kxpd)
      dimension mixup(nxhyd), sctd(nxyd)
! local data
      integer indx1, indy1, indx1y, nx, ny, nyh, nyhh, ny3, nxy, nxhy
      integer i, j, k, m, ks, km, kmr, nry, j1, j2, ns, ns2, k1, k2
      integer kxps, nryb, jj
      real at1, at2, at3, t2, t3, t4, t5, t6
      complex t1
      double precision sum1, sum2
      indx1 = indx - 1
      indy1 = indy - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      nyhh = ny/4
      ny3 = ny + 3
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 1
      kxps = kxpi + kxpp - 1
      if (kstrt.gt.nx) return
      if (isign.eq.0) return
      nryb = nxhy/nyh
      nry = nxy/nyh
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,m,jj,ns,ns2,km,kmr,k1,k2,j1,j2,at1,at2,at3,t2,t3,t4
!$OMP& ,t5,t6,t1,sum1,sum2)
      do 150 i = kxpi, kxps
! create auxiliary array in y
      kmr = nxy/ny
      sum1 = 0.5*(g(2,1,i) - g(2,ny+1,i))
      do 10 k = 2, nyh
      k1 = 1 + kmr*(k - 1)
      at3 = -aimag(sctd(k1))
      at2 = g(1,ny+2-k,i)
      at1 = g(1,k,i) + at2
      at2 = g(1,k,i) - at2
      at1 = at3*at1
      at2 = 0.5*at2
      g(1,k,i) = at1 + at2
      g(1,ny+2-k,i) = at1 - at2
      at2 = g(2,ny+2-k,i)
      at1 = g(2,k,i) + at2
      at2 = g(2,k,i) - at2
      sum1 = sum1 + real(sctd(k1))*at2
      at2 = at3*at2
      at1 = 0.5*at1
      g(2,k,i) = at1 - at2
      g(2,ny+2-k,i) = at1 + at2
   10 continue
      g(1,1,i) = 0.0
      g(1,nyh+1,i) = 2.0*g(1,nyh+1,i)
      g(2,1,i) = 0.5*(g(2,1,i) + g(2,ny+1,i))
      g(2,ny+1,i) = sum1
! bit-reverse array elements in y
      do 30 k = 1, nyh
      k1 = (mixup(k) - 1)/nryb + 1
      if (k.lt.k1) then
         do 20 jj = 1, 2
         t2 = g(jj,2*k1-1,i)
         t3 = g(jj,2*k1,i)
         g(jj,2*k1-1,i) = g(jj,2*k-1,i)
         g(jj,2*k1,i) = g(jj,2*k,i)
         g(jj,2*k-1,i) = t2
         g(jj,2*k,i) = t3
   20    continue
      endif
   30 continue
! then transform in y
      do 70 m = 1, indy1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyhh/ns
      kmr = 2*km*nry
      do 60 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 50 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 40 jj = 1, 2
      t2 = real(t1)*g(jj,2*j2-1,i) - aimag(t1)*g(jj,2*j2,i)
      t3 = aimag(t1)*g(jj,2*j2-1,i) + real(t1)*g(jj,2*j2,i)
      g(jj,2*j2-1,i) = g(jj,2*j1-1,i) - t2
      g(jj,2*j2,i) = g(jj,2*j1,i) - t3
      g(jj,2*j1-1,i) = g(jj,2*j1-1,i) + t2
      g(jj,2*j1,i) = g(jj,2*j1,i) + t3
   40 continue
   50 continue
   60 continue
   70 continue
! unscramble coefficients and normalize
! inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nyh
         do 90 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         do 80 jj = 1, 2
         t4 = g(jj,ny3-2*k,i)
         t5 = -g(jj,ny3-2*k+1,i)
         t2 = g(jj,2*k-1,i) + t4
         t3 = g(jj,2*k,i) + t5
         t6 = g(jj,2*k-1,i) - t4
         t5 = g(jj,2*k,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         g(jj,2*k-1,i) = 0.5*(t2 + t4)
         g(jj,2*k,i) = 0.5*(t3 + t5)
         g(jj,ny3-2*k,i) = 0.5*(t2 - t4)
         g(jj,ny3-2*k+1,i) = 0.5*(t5 - t3)
   80    continue
   90    continue
         do 100 jj = 1, 2
         g(jj,nyh+1,i) = g(jj,nyh+1,i)
         g(jj,nyh+2,i) = -g(jj,nyh+2,i)
         t2 = g(jj,1,i) + g(jj,2,i)
         g(jj,2,i) = g(jj,1,i) - g(jj,2,i)
         g(jj,1,i) = t2
         g(jj,ny+1,i) = g(jj,ny+1,i)
  100    continue
! forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nyh
         do 120 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         do 110 jj = 1, 2
         t4 = g(jj,ny3-2*k,i)
         t5 = -g(jj,ny3-2*k+1,i)
         t2 = g(jj,2*k-1,i) + t4
         t3 = g(jj,2*k,i) + t5
         t6 = g(jj,2*k-1,i) - t4
         t5 = g(jj,2*k,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         g(jj,2*k-1,i) = t2 + t4
         g(jj,2*k,i) = t3 + t5
         g(jj,ny3-2*k,i) = t2 - t4
         g(jj,ny3-2*k+1,i) = t5 - t3
  110    continue
  120    continue
         do 130 jj = 1, 2
         g(jj,nyh+1,i) = 2.0*g(jj,nyh+1,i)
         g(jj,nyh+2,i) = -2.0*g(jj,nyh+2,i)
         t2 = 2.0*(g(jj,1,i) + g(jj,2,i))
         g(jj,2,i) = 2.0*(g(jj,1,i) - g(jj,2,i))
         g(jj,1,i) = t2
         g(jj,ny+1,i) = 2.0*g(jj,ny+1,i)
  130    continue
      endif
! perform recursion for sine-cosine transform
      sum1 = 0.5*g(1,1,i)
      g(1,1,i) = 0.0
      g(1,2,i) = sum1
      sum2 = g(2,ny+1,i)
      g(2,ny+1,i) = g(2,2,i)
      g(2,2,i) = sum2
      do 140 k = 2, nyh
      sum1 = sum1 + g(1,2*k-1,i)
      g(1,2*k-1,i) = -g(1,2*k,i)
      g(1,2*k,i) = sum1
      sum2 = sum2 - g(2,2*k,i)
      g(2,2*k,i) = sum2
  140 continue
      g(1,ny+1,i) = 0.0
  150 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPFCST2RM2Y(g,isign,mixup,sctd,indx,indy,kstrt,kxpi,   &
     &kxpp,nyv,kxpd,nxhyd,nxyd)
! this subroutine performs the y part of 2 two dimensional fast real
! sine and cosine transforms and their inverses, for a subset of x,
! using real arithmetic, with OpenMP,
! for data which is distributed in blocks
! x component has a cosine transform, y component a sine transform
! algorithm is described in Numerical Recipies in Fortran, Second Ed.,
! by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
! [Cambridge Univ. Press, 1992], p. 508.
! for isign = (-1,1), input: all, output: f
! approximate flop count: N*(5*log2(N) + 18)/nvp
! where N = (nx/2)*ny
! indx/indy = exponent which determines length in x/y direction,
! where nx=2**indx, ny=2**indy
! if isign = -1, inverse cosine-sine transform are performed
! g(1,m,n) = (.5*g(1,1,n) + ((-1)**m)*g(1,ny+1,n)
!              + sum(g(1,k,n)*cos(pi*m*k/ny))
! g(2,m,n) = sum(g(2,k,n)*sin(pi*m*k/ny))
! if isign = 1, a forward cosine-sine transforms are performed
! g(1,k,n) = 2*(.5*g(1,1,n) + ((-1)**m)*g(1,ny+1,n)
!              + sum(g(1,m,n)*cos(pi*m*k/ny))
! g(2,k,n) = sum(g(2,m,n)*sin(pi*m*k/ny))
! mixup = array of bit reversed addresses
! sctd = sine/cosine table
! kstrt = starting data block number
! kxpi = initial x index used
! kxpp = number of x indices used
! nyv = second dimension of g >= ny + 1
! kxpd = third dimension of g >= kxp + 1
! nxhyd = maximum of (nx/2,ny)
! nxyd = maximum of (nx,ny)
! written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, kstrt, kxpi, kxpp
      integer nyv, kxpd, nxhyd, nxyd, mixup
      real g
      complex sctd
      dimension g(2,nyv,kxpd)
      dimension mixup(nxhyd), sctd(nxyd)
! local data
      integer indx1, indy1, indx1y, nx, ny, nyh, nyhh, ny3, nxy, nxhy
      integer i, j, k, m, ks, km, kmr, nry, j1, j2, ns, ns2, k1, k2
      integer kxps, nryb, jj
      real at1, at2, at3, t2, t3, t4, t5, t6
      complex t1
      double precision sum1, sum2
      indx1 = indx - 1
      indy1 = indy - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      nyhh = ny/4
      ny3 = ny + 3
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 1
      kxps = kxpi + kxpp - 1
      if (kstrt.gt.nx) return
      if (isign.eq.0) return
      nryb = nxhy/nyh
      nry = nxy/nyh
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,m,jj,ns,ns2,km,kmr,k1,k2,j1,j2,at1,at2,at3,t2,t3,t4
!$OMP& ,t5,t6,t1,sum1,sum2)
      do 150 i = kxpi, kxps
! create auxiliary array in y
      kmr = nxy/ny
      sum1 = 0.5*(g(1,1,i) - g(1,ny+1,i))
      do 10 k = 2, nyh
      k1 = 1 + kmr*(k - 1)
      at3 = -aimag(sctd(k1))
      at2 = g(1,ny+2-k,i)
      at1 = g(1,k,i) + at2
      at2 = g(1,k,i) - at2
      sum1 = sum1 + real(sctd(k1))*at2
      at2 = at3*at2
      at1 = 0.5*at1
      g(1,k,i) = at1 - at2
      g(1,ny+2-k,i) = at1 + at2
      at2 = g(2,ny+2-k,i)
      at1 = g(2,k,i) + at2
      at2 = g(2,k,i) - at2
      at1 = at3*at1
      at2 = 0.5*at2
      g(2,k,i) = at1 + at2
      g(2,ny+2-k,i) = at1 - at2
   10 continue
      g(1,1,i) = 0.5*(g(1,1,i) + g(1,ny+1,i))
      g(1,ny+1,i) = sum1
      g(2,1,i) = 0.0
      g(2,nyh+1,i) = 2.0*g(2,nyh+1,i)
! bit-reverse array elements in y
      do 30 k = 1, nyh
      k1 = (mixup(k) - 1)/nryb + 1
      if (k.lt.k1) then
         do 20 jj = 1, 2
         t2 = g(jj,2*k1-1,i)
         t3 = g(jj,2*k1,i)
         g(jj,2*k1-1,i) = g(jj,2*k-1,i)
         g(jj,2*k1,i) = g(jj,2*k,i)
         g(jj,2*k-1,i) = t2
         g(jj,2*k,i) = t3
   20    continue
      endif
   30 continue
! then transform in y
      do 70 m = 1, indy1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyhh/ns
      kmr = 2*km*nry
      do 60 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 50 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 40 jj = 1, 2
      t2 = real(t1)*g(jj,2*j2-1,i) - aimag(t1)*g(jj,2*j2,i)
      t3 = aimag(t1)*g(jj,2*j2-1,i) + real(t1)*g(jj,2*j2,i)
      g(jj,2*j2-1,i) = g(jj,2*j1-1,i) - t2
      g(jj,2*j2,i) = g(jj,2*j1,i) - t3
      g(jj,2*j1-1,i) = g(jj,2*j1-1,i) + t2
      g(jj,2*j1,i) = g(jj,2*j1,i) + t3
   40 continue
   50 continue
   60 continue
   70 continue
! unscramble coefficients and normalize
! inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nyh
         do 90 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         do 80 jj = 1, 2
         t4 = g(jj,ny3-2*k,i)
         t5 = -g(jj,ny3-2*k+1,i)
         t2 = g(jj,2*k-1,i) + t4
         t3 = g(jj,2*k,i) + t5
         t6 = g(jj,2*k-1,i) - t4
         t5 = g(jj,2*k,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         g(jj,2*k-1,i) = 0.5*(t2 + t4)
         g(jj,2*k,i) = 0.5*(t3 + t5)
         g(jj,ny3-2*k,i) = 0.5*(t2 - t4)
         g(jj,ny3-2*k+1,i) = 0.5*(t5 - t3)
   80    continue
   90    continue
         do 100 jj = 1, 2
         g(jj,nyh+1,i) = g(jj,nyh+1,i)
         g(jj,nyh+2,i) = -g(jj,nyh+2,i)
         t2 = g(jj,1,i) + g(jj,2,i)
         g(jj,2,i) = g(jj,1,i) - g(jj,2,i)
         g(jj,1,i) = t2
         g(jj,ny+1,i) = g(jj,ny+1,i)
  100    continue
! forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nyh
         do 120 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         do 110 jj = 1, 2
         t4 = g(jj,ny3-2*k,i)
         t5 = -g(jj,ny3-2*k+1,i)
         t2 = g(jj,2*k-1,i) + t4
         t3 = g(jj,2*k,i) + t5
         t6 = g(jj,2*k-1,i) - t4
         t5 = g(jj,2*k,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         g(jj,2*k-1,i) = t2 + t4
         g(jj,2*k,i) = t3 + t5
         g(jj,ny3-2*k,i) = t2 - t4
         g(jj,ny3-2*k+1,i) = t5 - t3
  110    continue
  120    continue
         do 130 jj = 1, 2
         g(jj,nyh+1,i) = 2.0*g(jj,nyh+1,i)
         g(jj,nyh+2,i) = -2.0*g(jj,nyh+2,i)
         t2 = 2.0*(g(jj,1,i) + g(jj,2,i))
         g(jj,2,i) = 2.0*(g(jj,1,i) - g(jj,2,i))
         g(jj,1,i) = t2
         g(jj,ny+1,i) = 2.0*g(jj,ny+1,i)
  130    continue
      endif
! perform recursion for sine-cosine transform
      sum1 = g(1,ny+1,i)
      g(1,ny+1,i) = g(1,2,i)
      g(1,2,i) = sum1
      sum2 = 0.5*g(2,1,i)
      g(2,1,i) = 0.0
      g(2,2,i) = sum2
      do 140 k = 2, nyh
      sum1 = sum1 - g(1,2*k,i)
      g(1,2*k,i) = sum1
      sum2 = sum2 + g(2,2*k-1,i)
      g(2,2*k-1,i) = -g(2,2*k,i)
      g(2,2*k,i) = sum2
  140 continue
      g(2,ny+1,i) = 0.0
  150 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine WPPFCST2RM3(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,&
     &indy,kstrt,nvp,nxvh,nyv,kxp2,kyp,kypd,kxp2d,nxhyd,nxyd)
! wrapper function for 3 parallel real cosine/sine transforms
! for the electric field with dirichlet or magnetic field with neumann
! boundary conditions
! x component has a cosine/sine transform in x and y, respectively
! y/z component has a sine/cosine transform in x and y, respectively
      implicit none
      integer isign, ntpose, mixup, indx, indy, kstrt, nvp, nxvh, nyv
      integer kxp2, kyp, kypd, kxp2d, nxhyd, nxyd
      real f, g, bs, br, ttp
      complex sctd
      dimension f(3,2*nxvh,kypd), g(3,nyv,kxp2d)
      dimension bs(3,kxp2+1,kyp+1), br(3,kxp2+1,kyp+1)
      dimension mixup(nxhyd), sctd(nxyd)
! local data
      integer nx, ny, kxpi, kypi, ks, kxp2p, kypp, kxb2, kyb
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
! calculate range of indices
      nx = 2**indx
      ny = 2**indy
! ks = processor id
      ks = kstrt - 1
! kxp2p = actual size used in x direction
      kxp2p = min(kxp2,max(0,nx-kxp2*ks))
! kypp = actual size used in y direction
      kypp = min(kyp,max(0,ny-kyp*ks))
! kxb2 = minimum number of processors needed in x direction
      kxb2 = (nx - 1)/kxp2 + 1
! kyb = minimum number of processors needed in y direction
      kyb = (ny - 1)/kyp + 1
! add extra word for last processor in x
      if (ks==(kxb2-1)) kxp2p = kxp2p + 1
! add extra word for last processor in y
      if (ks==(kyb-1)) kypp = kypp + 1
! inverse fourier transform
      if (isign.lt.0) then
! perform x cosine-sine-sine transform
         call PPFCSST2RM3X(f,isign,mixup,sctd,indx,indy,kstrt,kypi,kypp,&
     &nxvh,kypd,nxhyd,nxyd)
! transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PPRNTPOSE(f,g,bs,br,nx,ny,kxp2,kyp,kstrt,nvp,3,2*nxvh,nyv,&
     &kxp2d,kypd)
         call PWTIMERA(1,ttp,dtime)
! perform y sine-cosine-sine transform
         call PPFSCST2RM3Y(g,isign,mixup,sctd,indx,indy,kstrt,kxpi,kxp2p&
     &,nyv,kxp2d,nxhyd,nxyd)
! transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPRNTPOSE(g,f,br,bs,ny,nx,kyp,kxp2,kstrt,nvp,3,nyv,    &
     &2*nxvh,kypd,kxp2d)
            call PWTIMERA(1,tf,dtime)
         endif
! forward fourier transform
      else if (isign.gt.0) then
! transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPRNTPOSE(f,g,bs,br,nx,ny,kxp2,kyp,kstrt,nvp,3,2*nxvh, &
     &nyv,kxp2d,kypd)
            call PWTIMERA(1,tf,dtime)
         endif
! perform y sine-cosine-sine transform
         call PPFSCST2RM3Y(g,isign,mixup,sctd,indx,indy,kstrt,kxpi,kxp2p&
     &,nyv,kxp2d,nxhyd,nxyd)
! transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PPRNTPOSE(g,f,br,bs,ny,nx,kyp,kxp2,kstrt,nvp,3,nyv,2*nxvh,&
     &kypd,kxp2d)
         call PWTIMERA(1,ttp,dtime)
! perform x cosine-sine-sine transform
         call PPFCSST2RM3X(f,isign,mixup,sctd,indx,indy,kstrt,kypi,kypp,&
     &nxvh,kypd,nxhyd,nxyd)
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
!-----------------------------------------------------------------------
      subroutine WPPFSCT2RM3(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,&
     &indy,kstrt,nvp,nxvh,nyv,kxp2,kyp,kypd,kxp2d,nxhyd,nxyd)
! wrapper function for 3 parallel real sine/cosine transforms
! for the magnetic field with dirichlet or electric field with neumann
! boundary conditions
! x component has a sine/cosine transform in x and y, respectively
! y component has a cosine/sine transform in x and y, respectively
! z component has a cosine/cosine transform in x and y, respectively
      implicit none
      integer isign, ntpose, mixup, indx, indy, kstrt, nvp, nxvh, nyv
      integer kxp2, kyp, kypd, kxp2d, nxhyd, nxyd
      real f, g, bs, br, ttp
      complex sctd
      dimension f(3,2*nxvh,kypd), g(3,nyv,kxp2d)
      dimension bs(3,kxp2+1,kyp+1), br(3,kxp2+1,kyp+1)
      dimension mixup(nxhyd), sctd(nxyd)
! local data
      integer nx, ny, kxpi, kypi, ks, kxp2p, kypp, kxb2, kyb
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
! calculate range of indices
      nx = 2**indx
      ny = 2**indy
! ks = processor id
      ks = kstrt - 1
! kxp2p = actual size used in x direction
      kxp2p = min(kxp2,max(0,nx-kxp2*ks))
! kypp = actual size used in y direction
      kypp = min(kyp,max(0,ny-kyp*ks))
! kxb2 = minimum number of processors needed in x direction
      kxb2 = (nx - 1)/kxp2 + 1
! kyb = minimum number of processors needed in y direction
      kyb = (ny - 1)/kyp + 1
! add extra word for last processor in x
      if (ks==(kxb2-1)) kxp2p = kxp2p + 1
! add extra word for last processor in y
      if (ks==(kyb-1)) kypp = kypp + 1
! inverse fourier transform
      if (isign.lt.0) then
! perform x sine-cosine-cosine transform
         call PPFSCCT2RM3X(f,isign,mixup,sctd,indx,indy,kstrt,kypi,kypp,&
     &nxvh,kypd,nxhyd,nxyd)
! transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PPRNTPOSE(f,g,bs,br,nx,ny,kxp2,kyp,kstrt,nvp,3,2*nxvh,nyv,&
     &kxp2d,kypd)
         call PWTIMERA(1,ttp,dtime)
! perform y cosine-sine-cosine transform
         call PPFCSCT2RM3Y(g,isign,mixup,sctd,indx,indy,kstrt,kxpi,kxp2p&
     &,nyv,kxp2d,nxhyd,nxyd)
! transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPRNTPOSE(g,f,br,bs,ny,nx,kyp,kxp2,kstrt,nvp,3,nyv,    &
     &2*nxvh,kypd,kxp2d)
            call PWTIMERA(1,tf,dtime)
         endif
! forward fourier transform
      else if (isign.gt.0) then
! transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPRNTPOSE(f,g,bs,br,nx,ny,kxp2,kyp,kstrt,nvp,3,2*nxvh, &
     &nyv,kxp2d,kypd)
            call PWTIMERA(1,tf,dtime)
         endif
! perform y cosine-sine-cosine transform
         call PPFCSCT2RM3Y(g,isign,mixup,sctd,indx,indy,kstrt,kxpi,kxp2p&
     &,nyv,kxp2d,nxhyd,nxyd)
! transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PPRNTPOSE(g,f,br,bs,ny,nx,kyp,kxp2,kstrt,nvp,3,nyv,2*nxvh,&
     &kypd,kxp2d)
         call PWTIMERA(1,ttp,dtime)
! perform x sine-cosine-cosine transform
         call PPFSCCT2RM3X(f,isign,mixup,sctd,indx,indy,kstrt,kypi,kypp,&
     &nxvh,kypd,nxhyd,nxyd)
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
!-----------------------------------------------------------------------
      subroutine PPFCSST2RM3X(f,isign,mixup,sctd,indx,indy,kstrt,kypi,  &
     &kypp,nxvh,kypd,nxhyd,nxyd)
! this subroutine performs the x part of 3 two dimensional fast real
! sine and cosine transforms and their inverses, for a subset of y,
! using real arithmetic, with OpenMP,
! for data which is distributed in blocks
! x component has a cosine transform, y/z component a sine transform
! algorithm is described in Numerical Recipies in Fortran, Second Ed.,
! by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
! [Cambridge Univ. Press, 1992], p. 508.
! for isign = (-1,1), input: all, output: f
! approximate flop count: N*(5*log2(N) + 18)/nvp
! where N = (nx/2)*ny
! indx/indy = exponent which determines length in x/y direction,
! where nx=2**indx, ny=2**indy
! if isign = -1, inverse cosine-sine-sine transforms are performed
! f(1,n,k) = (1/nx*ny)*(.5*f(1,1,k) + ((-1)**n)*f(1,nx+1,k)
!              + sum(f(1,j,k)*cos(pi*n*j/nx)))
! f(2:3,n,k) = (1/nx*ny)*sum(f(2:3,j,k)*sin(pi*n*j/nx))
! if isign = 1, forward cosine-sine-sine transforms are performed
! f(1,j,k) = 2*(.5*f(1,1,k) + ((-1)**j)*f(1,n+1,k)
!              + sum(f(1,n,k)*cos(pi*n*j/nx))
! f(2:3,j,k) = sum(f(2:3,n,k)*sin(pi*n*j/nx))
! mixup = array of bit reversed addresses
! sctd = sine/cosine table
! kstrt = starting data block number
! kypi = initial y index used
! kypp = number of y indices used
! nxvh = first dimension of f >= nx/2 + 1
! kypd = second dimension of f >= kyp + 1
! nxhyd = maximum of (nx/2,ny)
! nxyd = maximum of (nx,ny)
! written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, kstrt, kypi, kypp
      integer nxvh, kypd, nxhyd, nxyd, mixup
      real f
      complex sctd
      dimension f(3,2*nxvh,kypd)
      dimension mixup(nxhyd), sctd(nxyd)
! local data
      integer indx1, indx1y, nx, nxh, nxhh, nx3, ny, nxy, nxhy, ks
      integer i, j, k, m, km, kmr, nrx, j1, j2, ns, ns2, k1, k2, kyps
      integer nrxb, jj
      real at1, at2, at3, t2, t3, t4, t5, t6, ani
      complex t1
      double precision sum1, sum2, sum3
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nx3 = nx + 3
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 1
      kyps = kypi + kypp - 1
      if (kstrt.gt.ny) return
      if (isign.eq.0) return
      ani = 0.5/(real(nx)*real(ny))
      nrxb = nxhy/nxh
      nrx = nxy/nxh
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,m,jj,ns,ns2,km,kmr,k1,k2,j1,j2,at1,at2,at3,t2,t3,t4
!$OMP& ,t5,t6,t1,sum1,sum2,sum3)
      do 150 i = kypi, kyps
! create auxiliary array in x
      kmr = nxy/nx
      sum1 = 0.5*(f(1,1,i) - f(1,nx+1,i))
      do 10 j = 2, nxh
      j1 = 1 + kmr*(j - 1)
      at3 = -aimag(sctd(j1))
      at2 = f(1,nx+2-j,i)
      at1 = f(1,j,i) + at2
      at2 = f(1,j,i) - at2
      sum1 = sum1 + real(sctd(j1))*at2
      at2 = at3*at2
      at1 = 0.5*at1
      f(1,j,i) = at1 - at2
      f(1,nx+2-j,i) = at1 + at2
      at2 = f(2,nx+2-j,i)
      at1 = f(2,j,i) + at2
      at2 = f(2,j,i) - at2
      at1 = at3*at1
      at2 = 0.5*at2
      f(2,j,i) = at1 + at2
      f(2,nx+2-j,i) = at1 - at2
      at2 = f(3,nx+2-j,i)
      at1 = f(3,j,i) + at2
      at2 = f(3,j,i) - at2
      at1 = at3*at1
      at2 = 0.5*at2
      f(3,j,i) = at1 + at2
      f(3,nx+2-j,i) = at1 - at2
   10 continue
      f(1,1,i) = 0.5*(f(1,1,i) + f(1,nx+1,i))
      f(1,nx+1,i) = sum1
      f(2,1,i) = 0.0
      f(2,nxh+1,i) = 2.0*f(2,nxh+1,i)
      f(3,1,i) = 0.0
      f(3,nxh+1,i) = 2.0*f(3,nxh+1,i)
! bit-reverse array elements in x
      do 30 j = 1, nxh
      j1 = (mixup(j) - 1)/nrxb + 1
      if (j.lt.j1) then
         do 20 jj = 1, 3
         t2 = f(jj,2*j1-1,i)
         t3 = f(jj,2*j1,i)
         f(jj,2*j1-1,i) = f(jj,2*j-1,i)
         f(jj,2*j1,i) = f(jj,2*j,i)
         f(jj,2*j-1,i) = t2
         f(jj,2*j,i) = t3
   20    continue
      endif
   30 continue
! then transform in x
      do 70 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = 2*km*nrx
      do 60 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 50 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 40 jj = 1, 3
      t2 = real(t1)*f(jj,2*j2-1,i) - aimag(t1)*f(jj,2*j2,i)
      t3 = aimag(t1)*f(jj,2*j2-1,i) + real(t1)*f(jj,2*j2,i)
      f(jj,2*j2-1,i) = f(jj,2*j1-1,i) - t2
      f(jj,2*j2,i) = f(jj,2*j1,i) - t3
      f(jj,2*j1-1,i) = f(jj,2*j1-1,i) + t2
      f(jj,2*j1,i) = f(jj,2*j1,i) + t3
   40 continue
   50 continue
   60 continue
   70 continue
! unscramble coefficients and normalize
! inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nxh
         do 90 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 80 jj = 1, 3
         t4 = f(jj,nx3-2*j,i)
         t5 = -f(jj,nx3-2*j+1,i)
         t2 = f(jj,2*j-1,i) + t4
         t3 = f(jj,2*j,i) + t5
         t6 = f(jj,2*j-1,i) - t4
         t5 = f(jj,2*j,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,2*j-1,i) = ani*(t2 + t4)
         f(jj,2*j,i) = ani*(t3 + t5)
         f(jj,nx3-2*j,i) = ani*(t2 - t4)
         f(jj,nx3-2*j+1,i) = ani*(t5 - t3)
   80    continue
   90    continue
         do 100 jj = 1, 3
         f(jj,nxh+1,i) = 2.0*ani*f(jj,nxh+1,i)
         f(jj,nxh+2,i) = -2.0*ani*f(jj,nxh+2,i)
         t2 = 2.0*ani*(f(jj,1,i) + f(jj,2,i))
         f(jj,2,i) = 2.0*ani*(f(jj,1,i) - f(jj,2,i))
         f(jj,1,i) = t2
         f(jj,nx+1,i) = 2.0*ani*f(jj,nx+1,i)
  100    continue
! forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nxh
         do 120 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 110 jj = 1, 3
         t4 = f(jj,nx3-2*j,i)
         t5 = -f(jj,nx3-2*j+1,i)
         t2 = f(jj,2*j-1,i) + t4
         t3 = f(jj,2*j,i) + t5
         t6 = f(jj,2*j-1,i) - t4
         t5 = f(jj,2*j,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,2*j-1,i) = t2 + t4
         f(jj,2*j,i) = t3 + t5
         f(jj,nx3-2*j,i) = t2 - t4
         f(jj,nx3-2*j+1,i) = t5 - t3
  110    continue
  120    continue
         do 130 jj = 1, 3
         f(jj,nxh+1,i) = 2.0*f(jj,nxh+1,i)
         f(jj,nxh+2,i) = -2.0*f(jj,nxh+2,i)
         t2 = 2.0*(f(jj,1,i) + f(jj,2,i))
         f(jj,2,i) = 2.0*(f(jj,1,i) - f(jj,2,i))
         f(jj,1,i) = t2
         f(jj,nx+1,i) = 2.0*f(jj,nx+1,i)
  130    continue
      endif
! perform recursion for cosine-sine transform
      sum1 = f(1,nx+1,i)
      f(1,nx+1,i) = f(1,2,i)
      f(1,2,i) = sum1
      sum2 = 0.5*f(2,1,i)
      f(2,1,i) = 0.0
      f(2,2,i) = sum2
      sum3 = 0.5*f(3,1,i)
      f(3,1,i) = 0.0
      f(3,2,i) = sum3
      do 140 j = 2, nxh
      sum1 = sum1 - f(1,2*j,i)
      f(1,2*j,i) = sum1
      sum2 = sum2 + f(2,2*j-1,i)
      f(2,2*j-1,i) = -f(2,2*j,i)
      f(2,2*j,i) = sum2
      sum3 = sum3 + f(3,2*j-1,i)
      f(3,2*j-1,i) = -f(3,2*j,i)
      f(3,2*j,i) = sum3
  140 continue
      f(2,nx+1,i) = 0.0
      f(3,nx+1,i) = 0.0
  150 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPFSCCT2RM3X(f,isign,mixup,sctd,indx,indy,kstrt,kypi,  &
     &kypp,nxvh,kypd,nxhyd,nxyd)
! this subroutine performs the x part of 3 two dimensional fast real
! sine and cosine transforms and their inverses, for a subset of y,
! using real arithmetic, with OpenMP,
! for data which is distributed in blocks
! x component has a sine transform, y/z component a cosine transform
! algorithm is described in Numerical Recipies in Fortran, Second Ed.,
! by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
! [Cambridge Univ. Press, 1992], p. 508.
! for isign = (-1,1), input: all, output: f
! approximate flop count: N*(5*log2(N) + 18)/nvp
! where N = (nx/2)*ny
! indx/indy = exponent which determines length in x/y direction,
! where nx=2**indx, ny=2**indy
! if isign = -1, inverse sine-cosine-cosine transforms are performed
! f(1,n,k) = (1/nx*ny)*sum(f(1,j,k)*sin(pi*n*j/nx))
! f(2:3,n,k) = (1/nx*ny)*(.5*f(2:3,1,k) + ((-1)**n)*f(2:3,nx+1,k)
!              + sum(f(2:3,j,k)*cos(pi*n*j/nx)))
! if isign = 1, forward sine-cosine-cosine transforms are performed
! f(1,j,k) = sum(f(1,n,k)*sin(pi*n*j/nx))
! f(2:3,j,k) = 2*(.5*f(2:3,1,k) + ((-1)**j)*f(2:3,n+1,k)
!              + sum(f(2:3,n,k)*cos(pi*n*j/nx))
! mixup = array of bit reversed addresses
! sctd = sine/cosine table
! kstrt = starting data block number
! kypi = initial y index used
! kypp = number of y indices used
! nxvh = first dimension of f >= nx/2 + 1
! kypd = second dimension of f >= kyp + 1
! nxhyd = maximum of (nx/2,ny)
! nxyd = maximum of (nx,ny)
! written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, kstrt, kypi, kypp
      integer nxvh, kypd, nxhyd, nxyd, mixup
      real f
      complex sctd
      dimension f(3,2*nxvh,kypd)
      dimension mixup(nxhyd), sctd(nxyd)
! local data
      integer indx1, indx1y, nx, nxh, nxhh, nx3, ny, nxy, nxhy, ks
      integer i, j, k, m, km, kmr, nrx, j1, j2, ns, ns2, k1, k2, kyps
      integer nrxb, jj
      real at1, at2, at3, t2, t3, t4, t5, t6, ani
      complex t1
      double precision sum1, sum2, sum3
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nx3 = nx + 3
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 1
      kyps = kypi + kypp - 1
      if (kstrt.gt.ny) return
      if (isign.eq.0) return
      ani = 0.5/(real(nx)*real(ny))
      nrxb = nxhy/nxh
      nrx = nxy/nxh
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,m,jj,ns,ns2,km,kmr,k1,k2,j1,j2,at1,at2,at3,t2,t3,t4
!$OMP& ,t5,t6,t1,sum1,sum2,sum3)
      do 150 i = kypi, kyps
! create auxiliary array in x
      kmr = nxy/nx
      sum1 = 0.5*(f(2,1,i) - f(2,nx+1,i))
      sum2 = 0.5*(f(3,1,i) - f(3,nx+1,i))
      do 10 j = 2, nxh
      j1 = 1 + kmr*(j - 1)
      at3 = -aimag(sctd(j1))
      at2 = f(1,nx+2-j,i)
      at1 = f(1,j,i) + at2
      at2 = f(1,j,i) - at2
      at1 = at3*at1
      at2 = 0.5*at2
      f(1,j,i) = at1 + at2
      f(1,nx+2-j,i) = at1 - at2
      at2 = f(2,nx+2-j,i)
      at1 = f(2,j,i) + at2
      at2 = f(2,j,i) - at2
      sum1 = sum1 + real(sctd(j1))*at2
      at2 = at3*at2
      at1 = 0.5*at1
      f(2,j,i) = at1 - at2
      f(2,nx+2-j,i) = at1 + at2
      at2 = f(3,nx+2-j,i)
      at1 = f(3,j,i) + at2
      at2 = f(3,j,i) - at2
      sum2 = sum2 + real(sctd(j1))*at2
      at2 = at3*at2
      at1 = 0.5*at1
      f(3,j,i) = at1 - at2
      f(3,nx+2-j,i) = at1 + at2
   10 continue
      f(1,1,i) = 0.0
      f(1,nxh+1,i) = 2.0*f(1,nxh+1,i)
      f(2,1,i) = 0.5*(f(2,1,i) + f(2,nx+1,i))
      f(2,nx+1,i) = sum1
      f(3,1,i) = 0.5*(f(3,1,i) + f(3,nx+1,i))
      f(3,nx+1,i) = sum2
! bit-reverse array elements in x
      do 30 j = 1, nxh
      j1 = (mixup(j) - 1)/nrxb + 1
      if (j.lt.j1) then
         do 20 jj = 1, 3
         t2 = f(jj,2*j1-1,i)
         t3 = f(jj,2*j1,i)
         f(jj,2*j1-1,i) = f(jj,2*j-1,i)
         f(jj,2*j1,i) = f(jj,2*j,i)
         f(jj,2*j-1,i) = t2
         f(jj,2*j,i) = t3
   20    continue
      endif
   30 continue
! then transform in x
      do 70 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = 2*km*nrx
      do 60 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 50 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 40 jj = 1, 3
      t2 = real(t1)*f(jj,2*j2-1,i) - aimag(t1)*f(jj,2*j2,i)
      t3 = aimag(t1)*f(jj,2*j2-1,i) + real(t1)*f(jj,2*j2,i)
      f(jj,2*j2-1,i) = f(jj,2*j1-1,i) - t2
      f(jj,2*j2,i) = f(jj,2*j1,i) - t3
      f(jj,2*j1-1,i) = f(jj,2*j1-1,i) + t2
      f(jj,2*j1,i) = f(jj,2*j1,i) + t3
   40 continue
   50 continue
   60 continue
   70 continue
! unscramble coefficients and normalize
! inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nxh
         do 90 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 80 jj = 1, 3
         t4 = f(jj,nx3-2*j,i)
         t5 = -f(jj,nx3-2*j+1,i)
         t2 = f(jj,2*j-1,i) + t4
         t3 = f(jj,2*j,i) + t5
         t6 = f(jj,2*j-1,i) - t4
         t5 = f(jj,2*j,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,2*j-1,i) = ani*(t2 + t4)
         f(jj,2*j,i) = ani*(t3 + t5)
         f(jj,nx3-2*j,i) = ani*(t2 - t4)
         f(jj,nx3-2*j+1,i) = ani*(t5 - t3)
   80    continue
   90    continue
         do 100 jj = 1, 3
         f(jj,nxh+1,i) = 2.0*ani*f(jj,nxh+1,i)
         f(jj,nxh+2,i) = -2.0*ani*f(jj,nxh+2,i)
         t2 = 2.0*ani*(f(jj,1,i) + f(jj,2,i))
         f(jj,2,i) = 2.0*ani*(f(jj,1,i) - f(jj,2,i))
         f(jj,1,i) = t2
         f(jj,nx+1,i) = 2.0*ani*f(jj,nx+1,i)
  100    continue
! forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nxh
         do 120 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 110 jj = 1, 3
         t4 = f(jj,nx3-2*j,i)
         t5 = -f(jj,nx3-2*j+1,i)
         t2 = f(jj,2*j-1,i) + t4
         t3 = f(jj,2*j,i) + t5
         t6 = f(jj,2*j-1,i) - t4
         t5 = f(jj,2*j,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,2*j-1,i) = t2 + t4
         f(jj,2*j,i) = t3 + t5
         f(jj,nx3-2*j,i) = t2 - t4
         f(jj,nx3-2*j+1,i) = t5 - t3
  110    continue
  120    continue
         do 130 jj = 1, 3
         f(jj,nxh+1,i) = 2.0*f(jj,nxh+1,i)
         f(jj,nxh+2,i) = -2.0*f(jj,nxh+2,i)
         t2 = 2.0*(f(jj,1,i) + f(jj,2,i))
         f(jj,2,i) = 2.0*(f(jj,1,i) - f(jj,2,i))
         f(jj,1,i) = t2
         f(jj,nx+1,i) = 2.0*f(jj,nx+1,i)
  130    continue
      endif
! perform recursion for cosine-sine transform
      sum1 = 0.5*f(1,1,i)
      f(1,1,i) = 0.0
      f(1,2,i) = sum1
      sum2 = f(2,nx+1,i)
      f(2,nx+1,i) = f(2,2,i)
      f(2,2,i) = sum2
      sum3 = f(3,nx+1,i)
      f(3,nx+1,i) = f(3,2,i)
      f(3,2,i) = sum3
      do 140 j = 2, nxh
      sum1 = sum1 + f(1,2*j-1,i)
      f(1,2*j-1,i) = -f(1,2*j,i)
      f(1,2*j,i) = sum1
      sum2 = sum2 - f(2,2*j,i)
      f(2,2*j,i) = sum2
      sum3 = sum3 - f(3,2*j,i)
      f(3,2*j,i) = sum3
  140 continue
      f(1,nx+1,i) = 0.0
  150 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPFSCST2RM3Y(g,isign,mixup,sctd,indx,indy,kstrt,kxpi,  &
     &kxpp,nyv,kxpd,nxhyd,nxyd)
! this subroutine performs the y part of 3 two dimensional fast real
! sine and cosine transforms and their inverses, for a subset of x,
! using real arithmetic, with OpenMP,
! for data which is distributed in blocks
! x/z component has a sine transform, y component a cosine transform
! algorithm is described in Numerical Recipies in Fortran, Second Ed.,
! by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
! [Cambridge Univ. Press, 1992], p. 508.
! for isign = (-1,1), input: all, output: f
! approximate flop count: N*(5*log2(N) + 18)/nvp
! where N = (nx/2)*ny
! indx/indy = exponent which determines length in x/y direction,
! where nx=2**indx, ny=2**indy
! if isign = -1, inverse sine-cosine-sine transform are performed
! g(1,m,n) = sum(g(1,k,n)*sin(pi*m*k/ny))
! g(2,m,n) = (.5*g(2,1,n) + ((-1)**m)*g(2,ny+1,n)
!              + sum(g(2,k,n)*cos(pi*m*k/ny))
! g(3,m,n) = sum(g(3,k,n)*sin(pi*m*k/ny))
! if isign = 1, a forward sine-cosine-sine transforms are performed
! g(1,k,n) = sum(g(1,m,n)*sin(pi*m*k/ny))
! g(2,k,n) = 2*(.5*g(2,1,n) + ((-1)**m)*g(2,ny+1,n)
!              + sum(g(2,m,n)*cos(pi*m*k/ny))
! g(3,k,n) = sum(g(3,m,n)*sin(pi*m*k/ny))
! mixup = array of bit reversed addresses
! sctd = sine/cosine table
! kstrt = starting data block number
! kxpi = initial x index used
! kxpp = number of x indices used
! nyv = first dimension of g >= ny + 1
! kxpd = second dimension of g >= kxp + 1
! nxhyd = maximum of (nx/2,ny)
! nxyd = maximum of (nx,ny)
! written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, kstrt, kxpi, kxpp
      integer nyv, kxpd, nxhyd, nxyd, mixup
      real g
      complex sctd
      dimension g(3,nyv,kxpd)
      dimension mixup(nxhyd), sctd(nxyd)
! local data
      integer indx1, indy1, indx1y, nx, ny, nyh, nyhh, ny3, nxy, nxhy
      integer i, j, k, m, ks, km, kmr, nry, j1, j2, ns, ns2, k1, k2
      integer kxps, nryb, jj
      real at1, at2, at3, t2, t3, t4, t5, t6
      complex t1
      double precision sum1, sum2, sum3
      indx1 = indx - 1
      indy1 = indy - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      nyhh = ny/4
      ny3 = ny + 3
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 1
      kxps = kxpi + kxpp - 1
      if (kstrt.gt.nx) return
      if (isign.eq.0) return
      nryb = nxhy/nyh
      nry = nxy/nyh
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,m,jj,ns,ns2,km,kmr,k1,k2,j1,j2,at1,at2,at3,t2,t3,t4
!$OMP& ,t5,t6,t1,sum1,sum2,sum3)
      do 150 i = kxpi, kxps
! create auxiliary array in y
      kmr = nxy/ny
      sum1 = 0.5*(g(2,1,i) - g(2,ny+1,i))
      do 10 k = 2, nyh
      k1 = 1 + kmr*(k - 1)
      at3 = -aimag(sctd(k1))
      at2 = g(1,ny+2-k,i)
      at1 = g(1,k,i) + at2
      at2 = g(1,k,i) - at2
      at1 = at3*at1
      at2 = 0.5*at2
      g(1,k,i) = at1 + at2
      g(1,ny+2-k,i) = at1 - at2
      at2 = g(2,ny+2-k,i)
      at1 = g(2,k,i) + at2
      at2 = g(2,k,i) - at2
      sum1 = sum1 + real(sctd(k1))*at2
      at2 = at3*at2
      at1 = 0.5*at1
      g(2,k,i) = at1 - at2
      g(2,ny+2-k,i) = at1 + at2
      at2 = g(3,ny+2-k,i)
      at1 = g(3,k,i) + at2
      at2 = g(3,k,i) - at2
      at1 = at3*at1
      at2 = 0.5*at2
      g(3,k,i) = at1 + at2
      g(3,ny+2-k,i) = at1 - at2
   10 continue
      g(1,1,i) = 0.0
      g(1,nyh+1,i) = 2.0*g(1,nyh+1,i)
      g(2,1,i) = 0.5*(g(2,1,i) + g(2,ny+1,i))
      g(2,ny+1,i) = sum1
      g(3,1,i) = 0.0
      g(3,nyh+1,i) = 2.0*g(3,nyh+1,i)
! bit-reverse array elements in y
      do 30 k = 1, nyh
      k1 = (mixup(k) - 1)/nryb + 1
      if (k.lt.k1) then
         do 20 jj = 1, 3
         t2 = g(jj,2*k1-1,i)
         t3 = g(jj,2*k1,i)
         g(jj,2*k1-1,i) = g(jj,2*k-1,i)
         g(jj,2*k1,i) = g(jj,2*k,i)
         g(jj,2*k-1,i) = t2
         g(jj,2*k,i) = t3
   20    continue
      endif
   30 continue
! then transform in y
      do 70 m = 1, indy1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyhh/ns
      kmr = 2*km*nry
      do 60 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 50 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 40 jj = 1, 3
      t2 = real(t1)*g(jj,2*j2-1,i) - aimag(t1)*g(jj,2*j2,i)
      t3 = aimag(t1)*g(jj,2*j2-1,i) + real(t1)*g(jj,2*j2,i)
      g(jj,2*j2-1,i) = g(jj,2*j1-1,i) - t2
      g(jj,2*j2,i) = g(jj,2*j1,i) - t3
      g(jj,2*j1-1,i) = g(jj,2*j1-1,i) + t2
      g(jj,2*j1,i) = g(jj,2*j1,i) + t3
   40 continue
   50 continue
   60 continue
   70 continue
! unscramble coefficients and normalize
! inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nyh
         do 90 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         do 80 jj = 1, 3
         t4 = g(jj,ny3-2*k,i)
         t5 = -g(jj,ny3-2*k+1,i)
         t2 = g(jj,2*k-1,i) + t4
         t3 = g(jj,2*k,i) + t5
         t6 = g(jj,2*k-1,i) - t4
         t5 = g(jj,2*k,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         g(jj,2*k-1,i) = 0.5*(t2 + t4)
         g(jj,2*k,i) = 0.5*(t3 + t5)
         g(jj,ny3-2*k,i) = 0.5*(t2 - t4)
         g(jj,ny3-2*k+1,i) = 0.5*(t5 - t3)
   80    continue
   90    continue
         do 100 jj = 1, 3
         g(jj,nyh+1,i) = g(jj,nyh+1,i)
         g(jj,nyh+2,i) = -g(jj,nyh+2,i)
         t2 = g(jj,1,i) + g(jj,2,i)
         g(jj,2,i) = g(jj,1,i) - g(jj,2,i)
         g(jj,1,i) = t2
         g(jj,ny+1,i) = g(jj,ny+1,i)
  100    continue
! forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nyh
         do 120 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         do 110 jj = 1, 3
         t4 = g(jj,ny3-2*k,i)
         t5 = -g(jj,ny3-2*k+1,i)
         t2 = g(jj,2*k-1,i) + t4
         t3 = g(jj,2*k,i) + t5
         t6 = g(jj,2*k-1,i) - t4
         t5 = g(jj,2*k,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         g(jj,2*k-1,i) = t2 + t4
         g(jj,2*k,i) = t3 + t5
         g(jj,ny3-2*k,i) = t2 - t4
         g(jj,ny3-2*k+1,i) = t5 - t3
  110    continue
  120    continue
         do 130 jj = 1, 3
         g(jj,nyh+1,i) = 2.0*g(jj,nyh+1,i)
         g(jj,nyh+2,i) = -2.0*g(jj,nyh+2,i)
         t2 = 2.0*(g(jj,1,i) + g(jj,2,i))
         g(jj,2,i) = 2.0*(g(jj,1,i) - g(jj,2,i))
         g(jj,1,i) = t2
         g(jj,ny+1,i) = 2.0*g(jj,ny+1,i)
  130    continue
      endif
! perform recursion for sine-cosine transform
      sum1 = 0.5*g(1,1,i)
      g(1,1,i) = 0.0
      g(1,2,i) = sum1
      sum2 = g(2,ny+1,i)
      g(2,ny+1,i) = g(2,2,i)
      g(2,2,i) = sum2
      sum3 = 0.5*g(3,1,i)
      g(3,1,i) = 0.0
      g(3,2,i) = sum3
      do 140 k = 2, nyh
      sum1 = sum1 + g(1,2*k-1,i)
      g(1,2*k-1,i) = -g(1,2*k,i)
      g(1,2*k,i) = sum1
      sum2 = sum2 - g(2,2*k,i)
      g(2,2*k,i) = sum2
      sum3 = sum3 + g(3,2*k-1,i)
      g(3,2*k-1,i) = -g(3,2*k,i)
      g(3,2*k,i) = sum3
  140 continue
      g(1,ny+1,i) = 0.0
      g(3,ny+1,i) = 0.0
  150 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPFCSCT2RM3Y(g,isign,mixup,sctd,indx,indy,kstrt,kxpi,  &
     &kxpp,nyv,kxpd,nxhyd,nxyd)
! this subroutine performs the y part of 3 two dimensional fast real
! sine and cosine transforms and their inverses, for a subset of x,
! using real arithmetic, with OpenMP,
! for data which is distributed in blocks
! x/z component has a cosine transform, y component a sine transform
! algorithm is described in Numerical Recipies in Fortran, Second Ed.,
! by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
! [Cambridge Univ. Press, 1992], p. 508.
! for isign = (-1,1), input: all, output: f
! approximate flop count: N*(5*log2(N) + 18)/nvp
! where N = (nx/2)*ny
! indx/indy = exponent which determines length in x/y direction,
! where nx=2**indx, ny=2**indy
! if isign = -1, inverse cosine-sine-cosine transform are performed
! g(1,m,n) = (.5*g(1,1,n) + ((-1)**m)*g(1,ny+1,n)
!              + sum(g(1,k,n)*cos(pi*m*k/ny))
! g(2,m,n) = sum(g(2,k,n)*sin(pi*m*k/ny))
! g(3,m,n) = (.5*g(3,1,n) + ((-1)**m)*g(3,ny+1,n)
!              + sum(g(3,k,n)*cos(pi*m*k/ny))
! if isign = 1, a forward cosine-sine-cosine transforms are performed
! g(1,k,n) = 2*(.5*g(1,1,n) + ((-1)**m)*g(1,ny+1,n)
!              + sum(g(1,m,n)*cos(pi*m*k/ny))
! g(2,k,n) = sum(g(2,m,n)*sin(pi*m*k/ny))
! g(3,k,n) = 2*(.5*g(3,1,n) + ((-1)**m)*g(3,ny+1,n)
!              + sum(g(3,m,n)*cos(pi*m*k/ny))
! mixup = array of bit reversed addresses
! sctd = sine/cosine table
! kstrt = starting data block number
! kxpi = initial x index used
! kxpp = number of x indices used
! nyv = first dimension of g >= ny + 1
! kxpd = second dimension of g >= kxp + 1
! nxhyd = maximum of (nx/2,ny)
! nxyd = maximum of (nx,ny)
! written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, kstrt, kxpi, kxpp
      integer nyv, kxpd, nxhyd, nxyd, mixup
      real g
      complex sctd
      dimension g(3,nyv,kxpd)
      dimension mixup(nxhyd), sctd(nxyd)
! local data
      integer indx1, indy1, indx1y, nx, ny, nyh, nyhh, ny3, nxy, nxhy
      integer i, j, k, m, ks, km, kmr, nry, j1, j2, ns, ns2, k1, k2
      integer kxps, nryb, jj
      real at1, at2, at3, t2, t3, t4, t5, t6
      complex t1
      double precision sum1, sum2, sum3
      indx1 = indx - 1
      indy1 = indy - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      nyhh = ny/4
      ny3 = ny + 3
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 1
      kxps = kxpi + kxpp - 1
      if (kstrt.gt.nx) return
      if (isign.eq.0) return
      nryb = nxhy/nyh
      nry = nxy/nyh
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,m,jj,ns,ns2,km,kmr,k1,k2,j1,j2,at1,at2,at3,t2,t3,t4
!$OMP& ,t5,t6,t1,sum1,sum2,sum3)
      do 150 i = kxpi, kxps
! create auxiliary array in y
      kmr = nxy/ny
      sum1 = 0.5*(g(1,1,i) - g(1,ny+1,i))
      sum2 = 0.5*(g(3,1,i) - g(3,ny+1,i))
      do 10 k = 2, nyh
      k1 = 1 + kmr*(k - 1)
      at3 = -aimag(sctd(k1))
      at2 = g(1,ny+2-k,i)
      at1 = g(1,k,i) + at2
      at2 = g(1,k,i) - at2
      sum1 = sum1 + real(sctd(k1))*at2
      at2 = at3*at2
      at1 = 0.5*at1
      g(1,k,i) = at1 - at2
      g(1,ny+2-k,i) = at1 + at2
      at2 = g(2,ny+2-k,i)
      at1 = g(2,k,i) + at2
      at2 = g(2,k,i) - at2
      at1 = at3*at1
      at2 = 0.5*at2
      g(2,k,i) = at1 + at2
      g(2,ny+2-k,i) = at1 - at2
      at2 = g(3,ny+2-k,i)
      at1 = g(3,k,i) + at2
      at2 = g(3,k,i) - at2
      sum2 = sum2 + real(sctd(k1))*at2
      at2 = at3*at2
      at1 = 0.5*at1
      g(3,k,i) = at1 - at2
      g(3,ny+2-k,i) = at1 + at2
   10 continue
      g(1,1,i) = 0.5*(g(1,1,i) + g(1,ny+1,i))
      g(1,ny+1,i) = sum1
      g(2,1,i) = 0.0
      g(2,nyh+1,i) = 2.0*g(2,nyh+1,i)
      g(3,1,i) = 0.5*(g(3,1,i) + g(3,ny+1,i))
      g(3,ny+1,i) = sum2
! bit-reverse array elements in y
      do 30 k = 1, nyh
      k1 = (mixup(k) - 1)/nryb + 1
      if (k.lt.k1) then
         do 20 jj = 1, 3
         t2 = g(jj,2*k1-1,i)
         t3 = g(jj,2*k1,i)
         g(jj,2*k1-1,i) = g(jj,2*k-1,i)
         g(jj,2*k1,i) = g(jj,2*k,i)
         g(jj,2*k-1,i) = t2
         g(jj,2*k,i) = t3
   20    continue
      endif
   30 continue
! then transform in y
      do 70 m = 1, indy1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyhh/ns
      kmr = 2*km*nry
      do 60 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 50 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 40 jj = 1, 3
      t2 = real(t1)*g(jj,2*j2-1,i) - aimag(t1)*g(jj,2*j2,i)
      t3 = aimag(t1)*g(jj,2*j2-1,i) + real(t1)*g(jj,2*j2,i)
      g(jj,2*j2-1,i) = g(jj,2*j1-1,i) - t2
      g(jj,2*j2,i) = g(jj,2*j1,i) - t3
      g(jj,2*j1-1,i) = g(jj,2*j1-1,i) + t2
      g(jj,2*j1,i) = g(jj,2*j1,i) + t3
   40 continue
   50 continue
   60 continue
   70 continue
! unscramble coefficients and normalize
! inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nyh
         do 90 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         do 80 jj = 1, 3
         t4 = g(jj,ny3-2*k,i)
         t5 = -g(jj,ny3-2*k+1,i)
         t2 = g(jj,2*k-1,i) + t4
         t3 = g(jj,2*k,i) + t5
         t6 = g(jj,2*k-1,i) - t4
         t5 = g(jj,2*k,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         g(jj,2*k-1,i) = 0.5*(t2 + t4)
         g(jj,2*k,i) = 0.5*(t3 + t5)
         g(jj,ny3-2*k,i) = 0.5*(t2 - t4)
         g(jj,ny3-2*k+1,i) = 0.5*(t5 - t3)
   80    continue
   90    continue
         do 100 jj = 1, 3
         g(jj,nyh+1,i) = g(jj,nyh+1,i)
         g(jj,nyh+2,i) = -g(jj,nyh+2,i)
         t2 = g(jj,1,i) + g(jj,2,i)
         g(jj,2,i) = g(jj,1,i) - g(jj,2,i)
         g(jj,1,i) = t2
         g(jj,ny+1,i) = g(jj,ny+1,i)
  100    continue
! forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nyh
         do 120 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         do 110 jj = 1, 3
         t4 = g(jj,ny3-2*k,i)
         t5 = -g(jj,ny3-2*k+1,i)
         t2 = g(jj,2*k-1,i) + t4
         t3 = g(jj,2*k,i) + t5
         t6 = g(jj,2*k-1,i) - t4
         t5 = g(jj,2*k,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         g(jj,2*k-1,i) = t2 + t4
         g(jj,2*k,i) = t3 + t5
         g(jj,ny3-2*k,i) = t2 - t4
         g(jj,ny3-2*k+1,i) = t5 - t3
  110    continue
  120    continue
         do 130 jj = 1, 3
         g(jj,nyh+1,i) = 2.0*g(jj,nyh+1,i)
         g(jj,nyh+2,i) = -2.0*g(jj,nyh+2,i)
         t2 = 2.0*(g(jj,1,i) + g(jj,2,i))
         g(jj,2,i) = 2.0*(g(jj,1,i) - g(jj,2,i))
         g(jj,1,i) = t2
         g(jj,ny+1,i) = 2.0*g(jj,ny+1,i)
  130    continue
      endif
! perform recursion for sine-cosine transform
      sum1 = g(1,ny+1,i)
      g(1,ny+1,i) = g(1,2,i)
      g(1,2,i) = sum1
      sum2 = 0.5*g(2,1,i)
      g(2,1,i) = 0.0
      g(2,2,i) = sum2
      sum3 = g(3,ny+1,i)
      g(3,ny+1,i) = g(3,2,i)
      g(3,2,i) = sum3
      do 140 k = 2, nyh
      sum1 = sum1 - g(1,2*k,i)
      g(1,2*k,i) = sum1
      sum2 = sum2 + g(2,2*k-1,i)
      g(2,2*k-1,i) = -g(2,2*k,i)
      g(2,2*k,i) = sum2
      sum3 = sum3 - g(3,2*k,i)
      g(3,2*k,i) = sum3
  140 continue
      g(2,ny+1,i) = 0.0
  150 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine WPPFSCT2RM4(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,&
     &indy,kstrt,nvp,nxvh,nyv,kxp2,kyp,kypd,kxp2d,nxhyd,nxyd)
! wrapper function for 4 parallel real sine/cosine transforms
! for the momentum flux with dirichlet boundary conditions
! x component has a sine/sine transform in x and y, respectively
! y component has a cosine/cosine transform in x and y, respectively
! z component has a cosine/sine transform in x and y, respectively
! w component has a sine/cosine transform in x and y, respectively
      implicit none
      integer isign, ntpose, mixup, indx, indy, kstrt, nvp, nxvh, nyv
      integer kxp2, kyp, kypd, kxp2d, nxhyd, nxyd
      real f, g, bs, br, ttp
      complex sctd
      dimension f(4,2*nxvh,kypd), g(4,nyv,kxp2d)
      dimension bs(4,kxp2+1,kyp+1), br(4,kxp2+1,kyp+1)
      dimension mixup(nxhyd), sctd(nxyd)
! local data
      integer nx, ny, kxpi, kypi, ks, kxp2p, kypp, kxb2, kyb
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
! calculate range of indices
      nx = 2**indx
      ny = 2**indy
! ks = processor id
      ks = kstrt - 1
! kxp2p = actual size used in x direction
      kxp2p = min(kxp2,max(0,nx-kxp2*ks))
! kypp = actual size used in y direction
      kypp = min(kyp,max(0,ny-kyp*ks))
! kxb2 = minimum number of processors needed in x direction
      kxb2 = (nx - 1)/kxp2 + 1
! kyb = minimum number of processors needed in y direction
      kyb = (ny - 1)/kyp + 1
! add extra word for last processor in x
      if (ks==(kxb2-1)) kxp2p = kxp2p + 1
! add extra word for last processor in y
      if (ks==(kyb-1)) kypp = kypp + 1
! inverse fourier transform
      if (isign.lt.0) then
! perform x sine-cosine-cosine-sine transform
         call PPFSCCST2RM4X(f,isign,mixup,sctd,indx,indy,kstrt,kypi,kypp&
     &,nxvh,kypd,nxhyd,nxyd)
! transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PPRNTPOSE(f,g,bs,br,nx,ny,kxp2,kyp,kstrt,nvp,4,2*nxvh,nyv,&
     &kxp2d,kypd)
         call PWTIMERA(1,ttp,dtime)
! perform y sine-cosine-sine-cosine transform
         call PPFSCSCT2RM4Y(g,isign,mixup,sctd,indx,indy,kstrt,kxpi,    &
     &kxp2p,nyv,kxp2d,nxhyd,nxyd)
! transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPRNTPOSE(g,f,br,bs,ny,nx,kyp,kxp2,kstrt,nvp,4,nyv,    &
     &2*nxvh,kypd,kxp2d)
            call PWTIMERA(1,tf,dtime)
         endif
! forward fourier transform
      else if (isign.gt.0) then
! transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPRNTPOSE(f,g,bs,br,nx,ny,kxp2,kyp,kstrt,nvp,4,2*nxvh, &
     &nyv,kxp2d,kypd)
            call PWTIMERA(1,tf,dtime)
         endif
! perform y sine-cosine-sine-cosine transform
         call PPFSCSCT2RM4Y(g,isign,mixup,sctd,indx,indy,kstrt,kxpi,    &
     &kxp2p,nyv,kxp2d,nxhyd,nxyd)
! transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PPRNTPOSE(g,f,br,bs,ny,nx,kyp,kxp2,kstrt,nvp,4,nyv,2*nxvh,&
     &kypd,kxp2d)
         call PWTIMERA(1,ttp,dtime)
! perform y sine-cosine-sine-cosine transform
         call PPFSCCST2RM4X(f,isign,mixup,sctd,indx,indy,kstrt,kypi,kypp&
     &,nxvh,kypd,nxhyd,nxyd)
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
!-----------------------------------------------------------------------
      subroutine PPFSCCST2RM4X(f,isign,mixup,sctd,indx,indy,kstrt,kypi, &
     &kypp,nxvh,kypd,nxhyd,nxyd)
! this subroutine performs the x part of 4 two dimensional fast real
! sine and cosine transforms and their inverses, for a subset of y,
! using real arithmetic, with OpenMP,
! for data which is distributed in blocks
! x/w component has a sine transform, y/z component a cosine transform
! algorithm is described in Numerical Recipies in Fortran, Second Ed.,
! by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
! [Cambridge Univ. Press, 1992], p. 508.
! for isign = (-1,1), input: all, output: f
! approximate flop count: N*(5*log2(N) + 18)/nvp
! where N = (nx/2)*ny
! indx/indy = exponent which determines length in x/y direction,
! where nx=2**indx, ny=2**indy
! if isign = -1, inverse sine-cosine-cosine-sine transforms are
! performed
! f(1,n,k) = (1/nx*ny)*sum(f(1,j,k)*sin(pi*n*j/nx))
! f(2:3,n,k) = (1/nx*ny)*(.5*f(2:3,1,k) + ((-1)**n)*f(2:3,nx+1,k)
!              + sum(f(2:3,j,k)*cos(pi*n*j/nx)))
! f(4,n,k) = (1/nx*ny)*sum(f(4,j,k)*sin(pi*n*j/nx))
! if isign = 1, forward sine-cosine-cosine-sine transforms are performed
! f(1,j,k) = sum(f(1,n,k)*sin(pi*n*j/nx))
! f(2:3,j,k) = 2*(.5*f(2:3,1,k) + ((-1)**j)*f(2:3,n+1,k)
!              + sum(f(2:3,n,k)*cos(pi*n*j/nx))
! f(4,j,k) = sum(f(4,n,k)*sin(pi*n*j/nx))
! mixup = array of bit reversed addresses
! sctd = sine/cosine table
! kstrt = starting data block number
! kypi = initial y index used
! kypp = number of y indices used
! nxvh = first dimension of f >= nx/2 + 1
! kypd = second dimension of f >= kyp + 1
! nxhyd = maximum of (nx/2,ny)
! nxyd = maximum of (nx,ny)
! written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, kstrt, kypi, kypp
      integer nxvh, kypd, nxhyd, nxyd, mixup
      real f
      complex sctd
      dimension f(4,2*nxvh,kypd)
      dimension mixup(nxhyd), sctd(nxyd)
! local data
      integer indx1, indx1y, nx, nxh, nxhh, nx3, ny, nxy, nxhy, ks
      integer i, j, k, m, km, kmr, nrx, j1, j2, ns, ns2, k1, k2, kyps
      integer nrxb, jj
      real at1, at2, at3, t2, t3, t4, t5, t6, ani
      complex t1
      double precision sum1, sum2, sum3, sum4
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nx3 = nx + 3
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 1
      kyps = kypi + kypp - 1
      if (kstrt.gt.ny) return
      if (isign.eq.0) return
      ani = 0.5/(real(nx)*real(ny))
      nrxb = nxhy/nxh
      nrx = nxy/nxh
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,m,jj,ns,ns2,km,kmr,k1,k2,j1,j2,at1,at2,at3,t2,t3,t4
!$OMP& ,t5,t6,t1,sum1,sum2,sum3,sum4)
      do 150 i = kypi, kyps
! create auxiliary array in x
      kmr = nxy/nx
      sum2 = 0.5*(f(2,1,i) - f(2,nx+1,i))
      sum3 = 0.5*(f(3,1,i) - f(3,nx+1,i))
      do 10 j = 2, nxh
      j1 = 1 + kmr*(j - 1)
      at3 = -aimag(sctd(j1))
      at2 = f(1,nx+2-j,i)
      at1 = f(1,j,i) + at2
      at2 = f(1,j,i) - at2
      at1 = at3*at1
      at2 = 0.5*at2
      f(1,j,i) = at1 + at2
      f(1,nx+2-j,i) = at1 - at2
      at2 = f(2,nx+2-j,i)
      at1 = f(2,j,i) + at2
      at2 = f(2,j,i) - at2
      sum2 = sum2 + real(sctd(j1))*at2
      at2 = at3*at2
      at1 = 0.5*at1
      f(2,j,i) = at1 - at2
      f(2,nx+2-j,i) = at1 + at2
      at2 = f(3,nx+2-j,i)
      at1 = f(3,j,i) + at2
      at2 = f(3,j,i) - at2
      sum3 = sum3 + real(sctd(j1))*at2
      at2 = at3*at2
      at1 = 0.5*at1
      f(3,j,i) = at1 - at2
      f(3,nx+2-j,i) = at1 + at2
      at2 = f(4,nx+2-j,i)
      at1 = f(4,j,i) + at2
      at2 = f(4,j,i) - at2
      at1 = at3*at1
      at2 = 0.5*at2
      f(4,j,i) = at1 + at2
      f(4,nx+2-j,i) = at1 - at2
   10 continue
      f(1,1,i) = 0.0
      f(1,nxh+1,i) = 2.0*f(1,nxh+1,i)
      f(2,1,i) = 0.5*(f(2,1,i) + f(2,nx+1,i))
      f(2,nx+1,i) = sum2
      f(3,1,i) = 0.5*(f(3,1,i) + f(3,nx+1,i))
      f(3,nx+1,i) = sum3
      f(4,1,i) = 0.0
      f(4,nxh+1,i) = 2.0*f(4,nxh+1,i)
! bit-reverse array elements in x
      do 30 j = 1, nxh
      j1 = (mixup(j) - 1)/nrxb + 1
      if (j.lt.j1) then
         do 20 jj = 1, 4
         t2 = f(jj,2*j1-1,i)
         t3 = f(jj,2*j1,i)
         f(jj,2*j1-1,i) = f(jj,2*j-1,i)
         f(jj,2*j1,i) = f(jj,2*j,i)
         f(jj,2*j-1,i) = t2
         f(jj,2*j,i) = t3
   20    continue
      endif
   30 continue
! then transform in x
      do 70 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = 2*km*nrx
      do 60 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 50 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 40 jj = 1, 4
      t2 = real(t1)*f(jj,2*j2-1,i) - aimag(t1)*f(jj,2*j2,i)
      t3 = aimag(t1)*f(jj,2*j2-1,i) + real(t1)*f(jj,2*j2,i)
      f(jj,2*j2-1,i) = f(jj,2*j1-1,i) - t2
      f(jj,2*j2,i) = f(jj,2*j1,i) - t3
      f(jj,2*j1-1,i) = f(jj,2*j1-1,i) + t2
      f(jj,2*j1,i) = f(jj,2*j1,i) + t3
   40 continue
   50 continue
   60 continue
   70 continue
! unscramble coefficients and normalize
! inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nxh
         do 90 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 80 jj = 1, 4
         t4 = f(jj,nx3-2*j,i)
         t5 = -f(jj,nx3-2*j+1,i)
         t2 = f(jj,2*j-1,i) + t4
         t3 = f(jj,2*j,i) + t5
         t6 = f(jj,2*j-1,i) - t4
         t5 = f(jj,2*j,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,2*j-1,i) = ani*(t2 + t4)
         f(jj,2*j,i) = ani*(t3 + t5)
         f(jj,nx3-2*j,i) = ani*(t2 - t4)
         f(jj,nx3-2*j+1,i) = ani*(t5 - t3)
   80    continue
   90    continue
         do 100 jj = 1, 4
         f(jj,nxh+1,i) = 2.0*ani*f(jj,nxh+1,i)
         f(jj,nxh+2,i) = -2.0*ani*f(jj,nxh+2,i)
         t2 = 2.0*ani*(f(jj,1,i) + f(jj,2,i))
         f(jj,2,i) = 2.0*ani*(f(jj,1,i) - f(jj,2,i))
         f(jj,1,i) = t2
         f(jj,nx+1,i) = 2.0*ani*f(jj,nx+1,i)
  100    continue
! forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nxh
         do 120 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 110 jj = 1, 4
         t4 = f(jj,nx3-2*j,i)
         t5 = -f(jj,nx3-2*j+1,i)
         t2 = f(jj,2*j-1,i) + t4
         t3 = f(jj,2*j,i) + t5
         t6 = f(jj,2*j-1,i) - t4
         t5 = f(jj,2*j,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,2*j-1,i) = t2 + t4
         f(jj,2*j,i) = t3 + t5
         f(jj,nx3-2*j,i) = t2 - t4
         f(jj,nx3-2*j+1,i) = t5 - t3
  110    continue
  120    continue
         do 130 jj = 1, 4
         f(jj,nxh+1,i) = 2.0*f(jj,nxh+1,i)
         f(jj,nxh+2,i) = -2.0*f(jj,nxh+2,i)
         t2 = 2.0*(f(jj,1,i) + f(jj,2,i))
         f(jj,2,i) = 2.0*(f(jj,1,i) - f(jj,2,i))
         f(jj,1,i) = t2
         f(jj,nx+1,i) = 2.0*f(jj,nx+1,i)
  130    continue
      endif
! perform recursion for cosine-sine transform
      sum1 = 0.5*f(1,1,i)
      f(1,1,i) = 0.0
      f(1,2,i) = sum1
      sum2 = f(2,nx+1,i)
      f(2,nx+1,i) = f(2,2,i)
      f(2,2,i) = sum2
      sum3 = f(3,nx+1,i)
      f(3,nx+1,i) = f(3,2,i)
      f(3,2,i) = sum3
      sum4 = 0.5*f(4,1,i)
      f(4,1,i) = 0.0
      f(4,2,i) = sum4
      do 140 j = 2, nxh
      sum1 = sum1 + f(1,2*j-1,i)
      f(1,2*j-1,i) = -f(1,2*j,i)
      f(1,2*j,i) = sum1
      sum2 = sum2 - f(2,2*j,i)
      f(2,2*j,i) = sum2
      sum3 = sum3 - f(3,2*j,i)
      f(3,2*j,i) = sum3
      sum4 = sum4 + f(4,2*j-1,i)
      f(4,2*j-1,i) = -f(4,2*j,i)
      f(4,2*j,i) = sum4
  140 continue
      f(1,nx+1,i) = 0.0
      f(4,nx+1,i) = 0.0
  150 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPFSCSCT2RM4Y(g,isign,mixup,sctd,indx,indy,kstrt,kxpi, &
     &kxpp,nyv,kxpd,nxhyd,nxyd)
! this subroutine performs the y part of 4 two dimensional fast real
! sine and cosine transforms and their inverses, for a subset of x,
! using real arithmetic, with OpenMP,
! for data which is distributed in blocks
! x/z component has a sine transform, y/w component a cosine transform
! algorithm is described in Numerical Recipies in Fortran, Second Ed.,
! by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
! [Cambridge Univ. Press, 1992], p. 508.
! for isign = (-1,1), input: all, output: f
! approximate flop count: N*(5*log2(N) + 18)/nvp
! where N = (nx/2)*ny
! indx/indy = exponent which determines length in x/y direction,
! where nx=2**indx, ny=2**indy
! if isign = -1, inverse sine-cosine-sine-cosine transform are performed
! g(1,m,n) = sum(g(1,k,n)*sin(pi*m*k/ny))
! g(2,m,n) = (.5*g(2,1,n) + ((-1)**m)*g(2,ny+1,n)
!              + sum(g(2,k,n)*cos(pi*m*k/ny))
! g(3,m,n) = sum(g(3,k,n)*sin(pi*m*k/ny))
! g(4,m,n) = (.5*g(4,1,n) + ((-1)**m)*g(4,ny+1,n)
!              + sum(g(4,k,n)*cos(pi*m*k/ny))
! if isign = 1, forward sine-cosine-sine-cosine transform are performed
! g(1,k,n) = sum(g(1,m,n)*sin(pi*m*k/ny))
! g(2,k,n) = 2*(.5*g(2,1,n) + ((-1)**m)*g(2,ny+1,n)
!              + sum(g(2,m,n)*cos(pi*m*k/ny))
! g(3,k,n) = sum(g(3,m,n)*sin(pi*m*k/ny))
! g(4,k,n) = 2*(.5*g(4,1,n) + ((-1)**m)*g(4,ny+1,n)
!              + sum(g(4,m,n)*cos(pi*m*k/ny))
! mixup = array of bit reversed addresses
! sctd = sine/cosine table
! kstrt = starting data block number
! kxpi = initial x index used
! kxpp = number of x indices used
! nyv = first dimension of g >= ny + 1
! kxpd = second dimension of g >= kxp + 1
! nxhyd = maximum of (nx/2,ny)
! nxyd = maximum of (nx,ny)
! written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, kstrt, kxpi, kxpp
      integer nyv, kxpd, nxhyd, nxyd, mixup
      real g
      complex sctd
      dimension g(4,nyv,kxpd)
      dimension mixup(nxhyd), sctd(nxyd)
! local data
      integer indx1, indy1, indx1y, nx, ny, nyh, nyhh, ny3, nxy, nxhy
      integer i, j, k, m, ks, km, kmr, nry, j1, j2, ns, ns2, k1, k2
      integer kxps, nryb, jj
      real at1, at2, at3, t2, t3, t4, t5, t6
      complex t1
      double precision sum1, sum2, sum3, sum4
      indx1 = indx - 1
      indy1 = indy - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      nyhh = ny/4
      ny3 = ny + 3
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 1
      kxps = kxpi + kxpp - 1
      if (kstrt.gt.nx) return
      if (isign.eq.0) return
      nryb = nxhy/nyh
      nry = nxy/nyh
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,m,jj,ns,ns2,km,kmr,k1,k2,j1,j2,at1,at2,at3,t2,t3,t4
!$OMP& ,t5,t6,t1,sum1,sum2,sum3,sum4)
      do 150 i = kxpi, kxps
! create auxiliary array in y
      kmr = nxy/ny
      sum2 = 0.5*(g(2,1,i) - g(2,ny+1,i))
      sum4 = 0.5*(g(4,1,i) - g(4,ny+1,i))
      do 10 k = 2, nyh
      k1 = 1 + kmr*(k - 1)
      at3 = -aimag(sctd(k1))
      at2 = g(1,ny+2-k,i)
      at1 = g(1,k,i) + at2
      at2 = g(1,k,i) - at2
      at1 = at3*at1
      at2 = 0.5*at2
      g(1,k,i) = at1 + at2
      g(1,ny+2-k,i) = at1 - at2
      at2 = g(2,ny+2-k,i)
      at1 = g(2,k,i) + at2
      at2 = g(2,k,i) - at2
      sum2 = sum2 + real(sctd(k1))*at2
      at2 = at3*at2
      at1 = 0.5*at1
      g(2,k,i) = at1 - at2
      g(2,ny+2-k,i) = at1 + at2
      at2 = g(3,ny+2-k,i)
      at1 = g(3,k,i) + at2
      at2 = g(3,k,i) - at2
      at1 = at3*at1
      at2 = 0.5*at2
      g(3,k,i) = at1 + at2
      g(3,ny+2-k,i) = at1 - at2
      at2 = g(4,ny+2-k,i)
      at1 = g(4,k,i) + at2
      at2 = g(4,k,i) - at2
      sum4 = sum4 + real(sctd(k1))*at2
      at2 = at3*at2
      at1 = 0.5*at1
      g(4,k,i) = at1 - at2
      g(4,ny+2-k,i) = at1 + at2
   10 continue
      g(1,1,i) = 0.0
      g(1,nyh+1,i) = 2.0*g(1,nyh+1,i)
      g(2,1,i) = 0.5*(g(2,1,i) + g(2,ny+1,i))
      g(2,ny+1,i) = sum2
      g(3,1,i) = 0.0
      g(3,nyh+1,i) = 2.0*g(3,nyh+1,i)
      g(4,1,i) = 0.5*(g(4,1,i) + g(4,ny+1,i))
      g(4,ny+1,i) = sum4
! bit-reverse array elements in y
      do 30 k = 1, nyh
      k1 = (mixup(k) - 1)/nryb + 1
      if (k.lt.k1) then
         do 20 jj = 1, 4
         t2 = g(jj,2*k1-1,i)
         t3 = g(jj,2*k1,i)
         g(jj,2*k1-1,i) = g(jj,2*k-1,i)
         g(jj,2*k1,i) = g(jj,2*k,i)
         g(jj,2*k-1,i) = t2
         g(jj,2*k,i) = t3
   20    continue
      endif
   30 continue
! then transform in y
      do 70 m = 1, indy1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyhh/ns
      kmr = 2*km*nry
      do 60 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 50 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 40 jj = 1, 4
      t2 = real(t1)*g(jj,2*j2-1,i) - aimag(t1)*g(jj,2*j2,i)
      t3 = aimag(t1)*g(jj,2*j2-1,i) + real(t1)*g(jj,2*j2,i)
      g(jj,2*j2-1,i) = g(jj,2*j1-1,i) - t2
      g(jj,2*j2,i) = g(jj,2*j1,i) - t3
      g(jj,2*j1-1,i) = g(jj,2*j1-1,i) + t2
      g(jj,2*j1,i) = g(jj,2*j1,i) + t3
   40 continue
   50 continue
   60 continue
   70 continue
! unscramble coefficients and normalize
! inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nyh
         do 90 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         do 80 jj = 1, 4
         t4 = g(jj,ny3-2*k,i)
         t5 = -g(jj,ny3-2*k+1,i)
         t2 = g(jj,2*k-1,i) + t4
         t3 = g(jj,2*k,i) + t5
         t6 = g(jj,2*k-1,i) - t4
         t5 = g(jj,2*k,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         g(jj,2*k-1,i) = 0.5*(t2 + t4)
         g(jj,2*k,i) = 0.5*(t3 + t5)
         g(jj,ny3-2*k,i) = 0.5*(t2 - t4)
         g(jj,ny3-2*k+1,i) = 0.5*(t5 - t3)
   80    continue
   90    continue
         do 100 jj = 1, 4
         g(jj,nyh+1,i) = g(jj,nyh+1,i)
         g(jj,nyh+2,i) = -g(jj,nyh+2,i)
         t2 = g(jj,1,i) + g(jj,2,i)
         g(jj,2,i) = g(jj,1,i) - g(jj,2,i)
         g(jj,1,i) = t2
         g(jj,ny+1,i) = g(jj,ny+1,i)
  100    continue
! forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nyh
         do 120 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         do 110 jj = 1, 4
         t4 = g(jj,ny3-2*k,i)
         t5 = -g(jj,ny3-2*k+1,i)
         t2 = g(jj,2*k-1,i) + t4
         t3 = g(jj,2*k,i) + t5
         t6 = g(jj,2*k-1,i) - t4
         t5 = g(jj,2*k,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         g(jj,2*k-1,i) = t2 + t4
         g(jj,2*k,i) = t3 + t5
         g(jj,ny3-2*k,i) = t2 - t4
         g(jj,ny3-2*k+1,i) = t5 - t3
  110    continue
  120    continue
         do 130 jj = 1, 4
         g(jj,nyh+1,i) = 2.0*g(jj,nyh+1,i)
         g(jj,nyh+2,i) = -2.0*g(jj,nyh+2,i)
         t2 = 2.0*(g(jj,1,i) + g(jj,2,i))
         g(jj,2,i) = 2.0*(g(jj,1,i) - g(jj,2,i))
         g(jj,1,i) = t2
         g(jj,ny+1,i) = 2.0*g(jj,ny+1,i)
  130    continue
      endif
! perform recursion for sine-cosine transform
      sum1 = 0.5*g(1,1,i)
      g(1,1,i) = 0.0
      g(1,2,i) = sum1
      sum2 = g(2,ny+1,i)
      g(2,ny+1,i) = g(2,2,i)
      g(2,2,i) = sum2
      sum3 = 0.5*g(3,1,i)
      g(3,1,i) = 0.0
      g(3,2,i) = sum3
      sum4 = g(4,ny+1,i)
      g(4,ny+1,i) = g(4,2,i)
      g(4,2,i) = sum4
      do 140 k = 2, nyh
      sum1 = sum1 + g(1,2*k-1,i)
      g(1,2*k-1,i) = -g(1,2*k,i)
      g(1,2*k,i) = sum1
      sum2 = sum2 - g(2,2*k,i)
      g(2,2*k,i) = sum2
      sum3 = sum3 + g(3,2*k-1,i)
      g(3,2*k-1,i) = -g(3,2*k,i)
      g(3,2*k,i) = sum3
      sum4 = sum4 - g(4,2*k,i)
      g(4,2*k,i) = sum4
  140 continue
      g(1,ny+1,i) = 0.0
      g(3,ny+1,i) = 0.0
  150 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine WPPFSCT2RM22(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx&
     &,indy,kstrt,nvp,nxvh,nyv,kxp2,kyp,kypd,kxp2d,nxhyd,nxyd)
! wrapper function for 2 parallel real sine/cosine transforms
! for the momentum flux with dirichlet boundary conditions
! x component has a sine/sine transform in x and y, respectively
! y component has a cosine/cosine transform in x and y, respectively
      implicit none
      integer isign, ntpose, mixup, indx, indy, kstrt, nvp, nxvh, nyv
      integer kxp2, kyp, kypd, kxp2d, nxhyd, nxyd
      real f, g, bs, br, ttp
      complex sctd
      dimension f(2,2*nxvh,kypd), g(2,nyv,kxp2d)
      dimension bs(2,kxp2+1,kyp+1), br(2,kxp2+1,kyp+1)
      dimension mixup(nxhyd), sctd(nxyd)
! local data
      integer nx, ny, kxpi, kypi, ks, kxp2p, kypp, kxb2, kyb
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
! calculate range of indices
      nx = 2**indx
      ny = 2**indy
! ks = processor id
      ks = kstrt - 1
! kxp2p = actual size used in x direction
      kxp2p = min(kxp2,max(0,nx-kxp2*ks))
! kypp = actual size used in y direction
      kypp = min(kyp,max(0,ny-kyp*ks))
! kxb2 = minimum number of processors needed in x direction
      kxb2 = (nx - 1)/kxp2 + 1
! kyb = minimum number of processors needed in y direction
      kyb = (ny - 1)/kyp + 1
! add extra word for last processor in x
      if (ks==(kxb2-1)) kxp2p = kxp2p + 1
! add extra word for last processor in y
      if (ks==(kyb-1)) kypp = kypp + 1
! inverse fourier transform
      if (isign.lt.0) then
! perform x sine-cosine transform
         call PPFSCCST2RM22X(f,isign,mixup,sctd,indx,indy,kstrt,kypi,   &
     &kypp,nxvh,kypd,nxhyd,nxyd)
! transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PPRNTPOSE(f,g,bs,br,nx,ny,kxp2,kyp,kstrt,nvp,2,2*nxvh,nyv,&
     &kxp2d,kypd)
         call PWTIMERA(1,ttp,dtime)
! perform y sine-cosine transform
         call PPFSCSCT2RM22Y(g,isign,mixup,sctd,indx,indy,kstrt,kxpi,   &
     &kxp2p,nyv,kxp2d,nxhyd,nxyd)
! transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPRNTPOSE(g,f,br,bs,ny,nx,kyp,kxp2,kstrt,nvp,2,nyv,    &
     &2*nxvh,kypd,kxp2d)
            call PWTIMERA(1,tf,dtime)
         endif
! forward fourier transform
      else if (isign.gt.0) then
! transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPRNTPOSE(f,g,bs,br,nx,ny,kxp2,kyp,kstrt,nvp,2,2*nxvh, &
     &nyv,kxp2d,kypd)
            call PWTIMERA(1,tf,dtime)
         endif
! perform y cosine-sine transform
         call PPFSCSCT2RM22Y(g,isign,mixup,sctd,indx,indy,kstrt,kxpi,   &
     &kxp2p,nyv,kxp2d,nxhyd,nxyd)
! transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PPRNTPOSE(g,f,br,bs,ny,nx,kyp,kxp2,kstrt,nvp,2,nyv,2*nxvh,&
     &kypd,kxp2d)
         call PWTIMERA(1,ttp,dtime)
! perform y sine-cosine transform
         call PPFSCCST2RM22X(f,isign,mixup,sctd,indx,indy,kstrt,kypi,   &
     &kypp,nxvh,kypd,nxhyd,nxyd)
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
!-----------------------------------------------------------------------
      subroutine PPFSCCST2RM22X(f,isign,mixup,sctd,indx,indy,kstrt,kypi,&
     &kypp,nxvh,kypd,nxhyd,nxyd)
! this subroutine performs the x part of 2 two dimensional fast real
! sine and cosine transforms and their inverses, for a subset of y,
! using real arithmetic, with OpenMP,
! for data which is distributed in blocks
! x component has a sine transform, y component a cosine transform
! algorithm is described in Numerical Recipies in Fortran, Second Ed.,
! by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
! [Cambridge Univ. Press, 1992], p. 508.
! for isign = (-1,1), input: all, output: f
! approximate flop count: N*(5*log2(N) + 18)/nvp
! where N = (nx/2)*ny
! indx/indy = exponent which determines length in x/y direction,
! where nx=2**indx, ny=2**indy
! if isign = -1, inverse sine-cosine transforms are
! performed
! f(1,n,k) = (1/nx*ny)*sum(f(1,j,k)*sin(pi*n*j/nx))
! f(2,n,k) = (1/nx*ny)*(.5*f(2,1,k) + ((-1)**n)*f(2,nx+1,k)
!              + sum(f(2,j,k)*cos(pi*n*j/nx)))
! if isign = 1, forward sine-cosine transforms are performed
! f(1,j,k) = sum(f(1,n,k)*sin(pi*n*j/nx))
! f(2,j,k) = 2*(.5*f(2,1,k) + ((-1)**j)*f(2,n+1,k)
!              + sum(f(2:3,n,k)*cos(pi*n*j/nx))
! mixup = array of bit reversed addresses
! sctd = sine/cosine table
! kstrt = starting data block number
! kypi = initial y index used
! kypp = number of y indices used
! nxvh = first dimension of f >= nx/2 + 1
! kypd = second dimension of f >= kyp + 1
! nxhyd = maximum of (nx/2,ny)
! nxyd = maximum of (nx,ny)
! written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, kstrt, kypi, kypp
      integer nxvh, kypd, nxhyd, nxyd, mixup
      real f
      complex sctd
      dimension f(2,2*nxvh,kypd)
      dimension mixup(nxhyd), sctd(nxyd)
! local data
      integer indx1, indx1y, nx, nxh, nxhh, nx3, ny, nxy, nxhy, ks
      integer i, j, k, m, km, kmr, nrx, j1, j2, ns, ns2, k1, k2, kyps
      integer nrxb, jj
      real at1, at2, at3, t2, t3, t4, t5, t6, ani
      complex t1
      double precision sum1, sum2
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nx3 = nx + 3
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 1
      kyps = kypi + kypp - 1
      if (kstrt.gt.ny) return
      if (isign.eq.0) return
      ani = 0.5/(real(nx)*real(ny))
      nrxb = nxhy/nxh
      nrx = nxy/nxh
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,m,jj,ns,ns2,km,kmr,k1,k2,j1,j2,at1,at2,at3,t2,t3,t4
!$OMP& ,t5,t6,t1,sum1,sum2)
      do 150 i = kypi, kyps
! create auxiliary array in x
      kmr = nxy/nx
      sum2 = 0.5*(f(2,1,i) - f(2,nx+1,i))
      do 10 j = 2, nxh
      j1 = 1 + kmr*(j - 1)
      at3 = -aimag(sctd(j1))
      at2 = f(1,nx+2-j,i)
      at1 = f(1,j,i) + at2
      at2 = f(1,j,i) - at2
      at1 = at3*at1
      at2 = 0.5*at2
      f(1,j,i) = at1 + at2
      f(1,nx+2-j,i) = at1 - at2
      at2 = f(2,nx+2-j,i)
      at1 = f(2,j,i) + at2
      at2 = f(2,j,i) - at2
      sum2 = sum2 + real(sctd(j1))*at2
      at2 = at3*at2
      at1 = 0.5*at1
      f(2,j,i) = at1 - at2
      f(2,nx+2-j,i) = at1 + at2
   10 continue
      f(1,1,i) = 0.0
      f(1,nxh+1,i) = 2.0*f(1,nxh+1,i)
      f(2,1,i) = 0.5*(f(2,1,i) + f(2,nx+1,i))
      f(2,nx+1,i) = sum2
! bit-reverse array elements in x
      do 30 j = 1, nxh
      j1 = (mixup(j) - 1)/nrxb + 1
      if (j.lt.j1) then
         do 20 jj = 1, 2
         t2 = f(jj,2*j1-1,i)
         t3 = f(jj,2*j1,i)
         f(jj,2*j1-1,i) = f(jj,2*j-1,i)
         f(jj,2*j1,i) = f(jj,2*j,i)
         f(jj,2*j-1,i) = t2
         f(jj,2*j,i) = t3
   20    continue
      endif
   30 continue
! then transform in x
      do 70 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = 2*km*nrx
      do 60 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 50 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 40 jj = 1, 2
      t2 = real(t1)*f(jj,2*j2-1,i) - aimag(t1)*f(jj,2*j2,i)
      t3 = aimag(t1)*f(jj,2*j2-1,i) + real(t1)*f(jj,2*j2,i)
      f(jj,2*j2-1,i) = f(jj,2*j1-1,i) - t2
      f(jj,2*j2,i) = f(jj,2*j1,i) - t3
      f(jj,2*j1-1,i) = f(jj,2*j1-1,i) + t2
      f(jj,2*j1,i) = f(jj,2*j1,i) + t3
   40 continue
   50 continue
   60 continue
   70 continue
! unscramble coefficients and normalize
! inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nxh
         do 90 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 80 jj = 1, 2
         t4 = f(jj,nx3-2*j,i)
         t5 = -f(jj,nx3-2*j+1,i)
         t2 = f(jj,2*j-1,i) + t4
         t3 = f(jj,2*j,i) + t5
         t6 = f(jj,2*j-1,i) - t4
         t5 = f(jj,2*j,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,2*j-1,i) = ani*(t2 + t4)
         f(jj,2*j,i) = ani*(t3 + t5)
         f(jj,nx3-2*j,i) = ani*(t2 - t4)
         f(jj,nx3-2*j+1,i) = ani*(t5 - t3)
   80    continue
   90    continue
         do 100 jj = 1, 2
         f(jj,nxh+1,i) = 2.0*ani*f(jj,nxh+1,i)
         f(jj,nxh+2,i) = -2.0*ani*f(jj,nxh+2,i)
         t2 = 2.0*ani*(f(jj,1,i) + f(jj,2,i))
         f(jj,2,i) = 2.0*ani*(f(jj,1,i) - f(jj,2,i))
         f(jj,1,i) = t2
         f(jj,nx+1,i) = 2.0*ani*f(jj,nx+1,i)
  100    continue
! forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nxh
         do 120 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 110 jj = 1, 2
         t4 = f(jj,nx3-2*j,i)
         t5 = -f(jj,nx3-2*j+1,i)
         t2 = f(jj,2*j-1,i) + t4
         t3 = f(jj,2*j,i) + t5
         t6 = f(jj,2*j-1,i) - t4
         t5 = f(jj,2*j,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,2*j-1,i) = t2 + t4
         f(jj,2*j,i) = t3 + t5
         f(jj,nx3-2*j,i) = t2 - t4
         f(jj,nx3-2*j+1,i) = t5 - t3
  110    continue
  120    continue
         do 130 jj = 1, 2
         f(jj,nxh+1,i) = 2.0*f(jj,nxh+1,i)
         f(jj,nxh+2,i) = -2.0*f(jj,nxh+2,i)
         t2 = 2.0*(f(jj,1,i) + f(jj,2,i))
         f(jj,2,i) = 2.0*(f(jj,1,i) - f(jj,2,i))
         f(jj,1,i) = t2
         f(jj,nx+1,i) = 2.0*f(jj,nx+1,i)
  130    continue
      endif
! perform recursion for cosine-sine transform
      sum1 = 0.5*f(1,1,i)
      f(1,1,i) = 0.0
      f(1,2,i) = sum1
      sum2 = f(2,nx+1,i)
      f(2,nx+1,i) = f(2,2,i)
      f(2,2,i) = sum2
      do 140 j = 2, nxh
      sum1 = sum1 + f(1,2*j-1,i)
      f(1,2*j-1,i) = -f(1,2*j,i)
      f(1,2*j,i) = sum1
      sum2 = sum2 - f(2,2*j,i)
      f(2,2*j,i) = sum2
  140 continue
      f(1,nx+1,i) = 0.0
  150 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPFSCSCT2RM22Y(g,isign,mixup,sctd,indx,indy,kstrt,kxpi,&
     &kxpp,nyv,kxpd,nxhyd,nxyd)
! this subroutine performs the y part of 2 two dimensional fast real
! sine and cosine transforms and their inverses, for a subset of x,
! using real arithmetic, with OpenMP,
! for data which is distributed in blocks
! x component has a sine transform, y component a cosine transform
! algorithm is described in Numerical Recipies in Fortran, Second Ed.,
! by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
! [Cambridge Univ. Press, 1992], p. 508.
! for isign = (-1,1), input: all, output: f
! approximate flop count: N*(5*log2(N) + 18)/nvp
! where N = (nx/2)*ny
! indx/indy = exponent which determines length in x/y direction,
! where nx=2**indx, ny=2**indy
! if isign = -1, inverse sine-cosine transform are performed
! g(1,m,n) = sum(g(1,k,n)*sin(pi*m*k/ny))
! g(2,m,n) = (.5*g(2,1,n) + ((-1)**m)*g(2,ny+1,n)
!              + sum(g(2,k,n)*cos(pi*m*k/ny))
! if isign = 1, forward sine-cosine transform are performed
! g(1,k,n) = sum(g(1,m,n)*sin(pi*m*k/ny))
! g(2,k,n) = 2*(.5*g(2,1,n) + ((-1)**m)*g(2,ny+1,n)
!              + sum(g(2,m,n)*cos(pi*m*k/ny))
! mixup = array of bit reversed addresses
! sctd = sine/cosine table
! kstrt = starting data block number
! kxpi = initial x index used
! kxpp = number of x indices used
! nyv = first dimension of g >= ny + 1
! kxpd = second dimension of g >= kxp + 1
! nxhyd = maximum of (nx/2,ny)
! nxyd = maximum of (nx,ny)
! written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, kstrt, kxpi, kxpp
      integer nyv, kxpd, nxhyd, nxyd, mixup
      real g
      complex sctd
      dimension g(2,nyv,kxpd)
      dimension mixup(nxhyd), sctd(nxyd)
! local data
      integer indx1, indy1, indx1y, nx, ny, nyh, nyhh, ny3, nxy, nxhy
      integer i, j, k, m, ks, km, kmr, nry, j1, j2, ns, ns2, k1, k2
      integer kxps, nryb, jj
      real at1, at2, at3, t2, t3, t4, t5, t6
      complex t1
      double precision sum1, sum2
      indx1 = indx - 1
      indy1 = indy - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      nyhh = ny/4
      ny3 = ny + 3
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 1
      kxps = kxpi + kxpp - 1
      if (kstrt.gt.nx) return
      if (isign.eq.0) return
      nryb = nxhy/nyh
      nry = nxy/nyh
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,m,jj,ns,ns2,km,kmr,k1,k2,j1,j2,at1,at2,at3,t2,t3,t4
!$OMP& ,t5,t6,t1,sum1,sum2)
      do 150 i = kxpi, kxps
! create auxiliary array in y
      kmr = nxy/ny
      sum2 = 0.5*(g(2,1,i) - g(2,ny+1,i))
      do 10 k = 2, nyh
      k1 = 1 + kmr*(k - 1)
      at3 = -aimag(sctd(k1))
      at2 = g(1,ny+2-k,i)
      at1 = g(1,k,i) + at2
      at2 = g(1,k,i) - at2
      at1 = at3*at1
      at2 = 0.5*at2
      g(1,k,i) = at1 + at2
      g(1,ny+2-k,i) = at1 - at2
      at2 = g(2,ny+2-k,i)
      at1 = g(2,k,i) + at2
      at2 = g(2,k,i) - at2
      sum2 = sum2 + real(sctd(k1))*at2
      at2 = at3*at2
      at1 = 0.5*at1
      g(2,k,i) = at1 - at2
      g(2,ny+2-k,i) = at1 + at2
   10 continue
      g(1,1,i) = 0.0
      g(1,nyh+1,i) = 2.0*g(1,nyh+1,i)
      g(2,1,i) = 0.5*(g(2,1,i) + g(2,ny+1,i))
      g(2,ny+1,i) = sum2
! bit-reverse array elements in y
      do 30 k = 1, nyh
      k1 = (mixup(k) - 1)/nryb + 1
      if (k.lt.k1) then
         do 20 jj = 1, 2
         t2 = g(jj,2*k1-1,i)
         t3 = g(jj,2*k1,i)
         g(jj,2*k1-1,i) = g(jj,2*k-1,i)
         g(jj,2*k1,i) = g(jj,2*k,i)
         g(jj,2*k-1,i) = t2
         g(jj,2*k,i) = t3
   20    continue
      endif
   30 continue
! then transform in y
      do 70 m = 1, indy1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyhh/ns
      kmr = 2*km*nry
      do 60 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 50 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 40 jj = 1, 2
      t2 = real(t1)*g(jj,2*j2-1,i) - aimag(t1)*g(jj,2*j2,i)
      t3 = aimag(t1)*g(jj,2*j2-1,i) + real(t1)*g(jj,2*j2,i)
      g(jj,2*j2-1,i) = g(jj,2*j1-1,i) - t2
      g(jj,2*j2,i) = g(jj,2*j1,i) - t3
      g(jj,2*j1-1,i) = g(jj,2*j1-1,i) + t2
      g(jj,2*j1,i) = g(jj,2*j1,i) + t3
   40 continue
   50 continue
   60 continue
   70 continue
! unscramble coefficients and normalize
! inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nyh
         do 90 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         do 80 jj = 1, 2
         t4 = g(jj,ny3-2*k,i)
         t5 = -g(jj,ny3-2*k+1,i)
         t2 = g(jj,2*k-1,i) + t4
         t3 = g(jj,2*k,i) + t5
         t6 = g(jj,2*k-1,i) - t4
         t5 = g(jj,2*k,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         g(jj,2*k-1,i) = 0.5*(t2 + t4)
         g(jj,2*k,i) = 0.5*(t3 + t5)
         g(jj,ny3-2*k,i) = 0.5*(t2 - t4)
         g(jj,ny3-2*k+1,i) = 0.5*(t5 - t3)
   80    continue
   90    continue
         do 100 jj = 1, 2
         g(jj,nyh+1,i) = g(jj,nyh+1,i)
         g(jj,nyh+2,i) = -g(jj,nyh+2,i)
         t2 = g(jj,1,i) + g(jj,2,i)
         g(jj,2,i) = g(jj,1,i) - g(jj,2,i)
         g(jj,1,i) = t2
         g(jj,ny+1,i) = g(jj,ny+1,i)
  100    continue
! forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nyh
         do 120 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         do 110 jj = 1, 2
         t4 = g(jj,ny3-2*k,i)
         t5 = -g(jj,ny3-2*k+1,i)
         t2 = g(jj,2*k-1,i) + t4
         t3 = g(jj,2*k,i) + t5
         t6 = g(jj,2*k-1,i) - t4
         t5 = g(jj,2*k,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         g(jj,2*k-1,i) = t2 + t4
         g(jj,2*k,i) = t3 + t5
         g(jj,ny3-2*k,i) = t2 - t4
         g(jj,ny3-2*k+1,i) = t5 - t3
  110    continue
  120    continue
         do 130 jj = 1, 2
         g(jj,nyh+1,i) = 2.0*g(jj,nyh+1,i)
         g(jj,nyh+2,i) = -2.0*g(jj,nyh+2,i)
         t2 = 2.0*(g(jj,1,i) + g(jj,2,i))
         g(jj,2,i) = 2.0*(g(jj,1,i) - g(jj,2,i))
         g(jj,1,i) = t2
         g(jj,ny+1,i) = 2.0*g(jj,ny+1,i)
  130    continue
      endif
! perform recursion for sine-cosine transform
      sum1 = 0.5*g(1,1,i)
      g(1,1,i) = 0.0
      g(1,2,i) = sum1
      sum2 = g(2,ny+1,i)
      g(2,ny+1,i) = g(2,2,i)
      g(2,2,i) = sum2
      do 140 k = 2, nyh
      sum1 = sum1 + g(1,2*k-1,i)
      g(1,2*k-1,i) = -g(1,2*k,i)
      g(1,2*k,i) = sum1
      sum2 = sum2 - g(2,2*k,i)
      g(2,2*k,i) = sum2
  140 continue
      g(1,ny+1,i) = 0.0
  150 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine WPPFSST2RM23(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx&
     &,indy,kstrt,nvp,nxvh,nyv,kxp2,kyp,kypd,kxp2d,nxhyd,nxyd)
! wrapper function for 3 parallel real sine transforms
! x/y component has a sine/sine transform in x and y, respectively
! z component has a cosine/cosine transform in x and y, respectively
      implicit none
      integer isign, ntpose, mixup, indx, indy, kstrt, nvp, nxvh, nyv
      integer kxp2, kyp, kypd, kxp2d, nxhyd, nxyd
      real f, g, bs, br, ttp
      complex sctd
      dimension f(3,2*nxvh,kypd), g(3,nyv,kxp2d)
      dimension bs(3,kxp2+1,kyp+1), br(3,kxp2+1,kyp+1)
      dimension mixup(nxhyd), sctd(nxyd)
! local data
      integer nx, ny, kxpi, kypi, ks, kxp2p, kypp, kxb2, kyb
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
! calculate range of indices
      nx = 2**indx
      ny = 2**indy
! ks = processor id
      ks = kstrt - 1
! kxp2p = actual size used in x direction
      kxp2p = min(kxp2,max(0,nx-kxp2*ks))
! kypp = actual size used in y direction
      kypp = min(kyp,max(0,ny-kyp*ks))
! kxb2 = minimum number of processors needed in x direction
      kxb2 = (nx - 1)/kxp2 + 1
! kyb = minimum number of processors needed in y direction
      kyb = (ny - 1)/kyp + 1
! add extra word for last processor in x
      if (ks==(kxb2-1)) kxp2p = kxp2p + 1
! add extra word for last processor in y
      if (ks==(kyb-1)) kypp = kypp + 1
! inverse fourier transform
      if (isign.lt.0) then
! perform x sine transforms
         call PPFSSCT2RM23X(f,isign,mixup,sctd,indx,indy,kstrt,kypi,kypp&
     &,nxvh,kypd,nxhyd,nxyd)
! transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PPRNTPOSE(f,g,bs,br,nx,ny,kxp2,kyp,kstrt,nvp,3,2*nxvh,nyv,&
     &kxp2d,kypd)
         call PWTIMERA(1,ttp,dtime)
! perform y sine transforms
         call PPFSSCT2RM23Y(g,isign,mixup,sctd,indx,indy,kstrt,kxpi,    &
     &kxp2p,nyv,kxp2d,nxhyd,nxyd)
! transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPRNTPOSE(g,f,br,bs,ny,nx,kyp,kxp2,kstrt,nvp,3,nyv,    &
     &2*nxvh,kypd,kxp2d)
            call PWTIMERA(1,tf,dtime)
         endif
! forward fourier transform
      else if (isign.gt.0) then
! transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PPRNTPOSE(f,g,bs,br,nx,ny,kxp2,kyp,kstrt,nvp,3,2*nxvh, &
     &nyv,kxp2d,kypd)
            call PWTIMERA(1,tf,dtime)
         endif
! perform y sine transforms
         call PPFSSCT2RM23Y(g,isign,mixup,sctd,indx,indy,kstrt,kxpi,    &
     &kxp2p,nyv,kxp2d,nxhyd,nxyd)
! transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PPRNTPOSE(g,f,br,bs,ny,nx,kyp,kxp2,kstrt,nvp,3,nyv,2*nxvh,&
     &kypd,kxp2d)
         call PWTIMERA(1,ttp,dtime)
! perform x sine transforms
         call PPFSSCT2RM23X(f,isign,mixup,sctd,indx,indy,kstrt,kypi,kypp&
     &,nxvh,kypd,nxhyd,nxyd)
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
!-----------------------------------------------------------------------
      subroutine PPFSSCT2RM23X(f,isign,mixup,sctd,indx,indy,kstrt,kypi, &
     &kypp,nxvh,kypd,nxhyd,nxyd)
! this subroutine performs the x part of 3 two dimensional fast real
! sine transforms and their inverses, for a subset of y,
! using real arithmetic, with OpenMP,
! for data which is distributed in blocks
! x/y component has a sine transform, z component a cosine transform
! algorithm is described in Numerical Recipies in Fortran, Second Ed.,
! by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
! [Cambridge Univ. Press, 1992], p. 508.
! for isign = (-1,1), input: all, output: f
! approximate flop count: N*(5*log2(N) + 18)/nvp
! where N = (nx/2)*ny
! indx/indy = exponent which determines length in x/y direction,
! where nx=2**indx, ny=2**indy
! if isign = -1, inverse sine-sine-cosine transforms are performed
! f(1:2,n,k) = (1/nx*ny)*sum(f(2:3,j,k)*sin(pi*n*j/nx))
! f(3,n,k) = (1/nx*ny)*(.5*f(1,1,k) + ((-1)**n)*f(1,nx+1,k)
!              + sum(f(1,j,k)*cos(pi*n*j/nx)))
! if isign = 1, forward sine-sine-cosine transforms are performed
! f(1:2,j,k) = sum(f(2:3,n,k)*sin(pi*n*j/nx))
! f(3,j,k) = 2*(.5*f(1,1,k) + ((-1)**j)*f(1,n+1,k)
!              + sum(f(1,n,k)*cos(pi*n*j/nx))
! mixup = array of bit reversed addresses
! sctd = sine/cosine table
! kstrt = starting data block number
! kypi = initial y index used
! kypp = number of y indices used
! nxvh = first dimension of f >= nx/2 + 1
! kypd = second dimension of f >= kyp + 1
! nxhyd = maximum of (nx/2,ny)
! nxyd = maximum of (nx,ny)
! written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, kstrt, kypi, kypp
      integer nxvh, kypd, nxhyd, nxyd, mixup
      real f
      complex sctd
      dimension f(3,2*nxvh,kypd)
      dimension mixup(nxhyd), sctd(nxyd)
! local data
      integer indx1, indx1y, nx, nxh, nxhh, nx3, ny, nxy, nxhy, ks
      integer i, j, k, m, km, kmr, nrx, j1, j2, ns, ns2, k1, k2, kyps
      integer nrxb, jj
      real at1, at2, at3, t2, t3, t4, t5, t6, ani
      complex t1
      double precision sum1, sum2, sum3
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nx3 = nx + 3
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 1
      kyps = kypi + kypp - 1
      if (kstrt.gt.ny) return
      if (isign.eq.0) return
      ani = 0.5/(real(nx)*real(ny))
      nrxb = nxhy/nxh
      nrx = nxy/nxh
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,m,jj,ns,ns2,km,kmr,k1,k2,j1,j2,at1,at2,at3,t2,t3,t4
!$OMP& ,t5,t6,t1,sum1,sum2,sum3)
      do 150 i = kypi, kyps
! create auxiliary array in x
      kmr = nxy/nx
      sum3 = 0.5*(f(3,1,i) - f(3,nx+1,i))
      do 10 j = 2, nxh
      j1 = 1 + kmr*(j - 1)
      at3 = -aimag(sctd(j1))
      at2 = f(1,nx+2-j,i)
      at1 = f(1,j,i) + at2
      at2 = f(1,j,i) - at2
      at1 = at3*at1
      at2 = 0.5*at2
      f(1,j,i) = at1 + at2
      f(1,nx+2-j,i) = at1 - at2
      at2 = f(2,nx+2-j,i)
      at1 = f(2,j,i) + at2
      at2 = f(2,j,i) - at2
      at1 = at3*at1
      at2 = 0.5*at2
      f(2,j,i) = at1 + at2
      f(2,nx+2-j,i) = at1 - at2
      at2 = f(3,nx+2-j,i)
      at1 = f(3,j,i) + at2
      at2 = f(3,j,i) - at2
      sum3 = sum3 + real(sctd(j1))*at2
      at2 = at3*at2
      at1 = 0.5*at1
      f(3,j,i) = at1 - at2
      f(3,nx+2-j,i) = at1 + at2
   10 continue
      f(1,1,i) = 0.0
      f(1,nxh+1,i) = 2.0*f(1,nxh+1,i)
      f(2,1,i) = 0.0
      f(2,nxh+1,i) = 2.0*f(2,nxh+1,i)
      f(3,1,i) = 0.5*(f(3,1,i) + f(3,nx+1,i))
      f(3,nx+1,i) = sum3
! bit-reverse array elements in x
      do 30 j = 1, nxh
      j1 = (mixup(j) - 1)/nrxb + 1
      if (j.lt.j1) then
         do 20 jj = 1, 3
         t2 = f(jj,2*j1-1,i)
         t3 = f(jj,2*j1,i)
         f(jj,2*j1-1,i) = f(jj,2*j-1,i)
         f(jj,2*j1,i) = f(jj,2*j,i)
         f(jj,2*j-1,i) = t2
         f(jj,2*j,i) = t3
   20    continue
      endif
   30 continue
! then transform in x
      do 70 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = 2*km*nrx
      do 60 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 50 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 40 jj = 1, 3
      t2 = real(t1)*f(jj,2*j2-1,i) - aimag(t1)*f(jj,2*j2,i)
      t3 = aimag(t1)*f(jj,2*j2-1,i) + real(t1)*f(jj,2*j2,i)
      f(jj,2*j2-1,i) = f(jj,2*j1-1,i) - t2
      f(jj,2*j2,i) = f(jj,2*j1,i) - t3
      f(jj,2*j1-1,i) = f(jj,2*j1-1,i) + t2
      f(jj,2*j1,i) = f(jj,2*j1,i) + t3
   40 continue
   50 continue
   60 continue
   70 continue
! unscramble coefficients and normalize
! inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nxh
         do 90 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 80 jj = 1, 3
         t4 = f(jj,nx3-2*j,i)
         t5 = -f(jj,nx3-2*j+1,i)
         t2 = f(jj,2*j-1,i) + t4
         t3 = f(jj,2*j,i) + t5
         t6 = f(jj,2*j-1,i) - t4
         t5 = f(jj,2*j,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,2*j-1,i) = ani*(t2 + t4)
         f(jj,2*j,i) = ani*(t3 + t5)
         f(jj,nx3-2*j,i) = ani*(t2 - t4)
         f(jj,nx3-2*j+1,i) = ani*(t5 - t3)
   80    continue
   90    continue
         do 100 jj = 1, 3
         f(jj,nxh+1,i) = 2.0*ani*f(jj,nxh+1,i)
         f(jj,nxh+2,i) = -2.0*ani*f(jj,nxh+2,i)
         t2 = 2.0*ani*(f(jj,1,i) + f(jj,2,i))
         f(jj,2,i) = 2.0*ani*(f(jj,1,i) - f(jj,2,i))
         f(jj,1,i) = t2
         f(jj,nx+1,i) = 2.0*ani*f(jj,nx+1,i)
  100    continue
! forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nxh
         do 120 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 110 jj = 1, 3
         t4 = f(jj,nx3-2*j,i)
         t5 = -f(jj,nx3-2*j+1,i)
         t2 = f(jj,2*j-1,i) + t4
         t3 = f(jj,2*j,i) + t5
         t6 = f(jj,2*j-1,i) - t4
         t5 = f(jj,2*j,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,2*j-1,i) = t2 + t4
         f(jj,2*j,i) = t3 + t5
         f(jj,nx3-2*j,i) = t2 - t4
         f(jj,nx3-2*j+1,i) = t5 - t3
  110    continue
  120    continue
         do 130 jj = 1, 3
         f(jj,nxh+1,i) = 2.0*f(jj,nxh+1,i)
         f(jj,nxh+2,i) = -2.0*f(jj,nxh+2,i)
         t2 = 2.0*(f(jj,1,i) + f(jj,2,i))
         f(jj,2,i) = 2.0*(f(jj,1,i) - f(jj,2,i))
         f(jj,1,i) = t2
         f(jj,nx+1,i) = 2.0*f(jj,nx+1,i)
  130    continue
      endif
! perform recursion for sine-sine transform
      sum1 = 0.5*f(1,1,i)
      f(1,1,i) = 0.0
      f(1,2,i) = sum1
      sum2 = 0.5*f(2,1,i)
      f(2,1,i) = 0.0
      f(2,2,i) = sum2
      sum3 = f(3,nx+1,i)
      f(3,nx+1,i) = f(3,2,i)
      f(3,2,i) = sum3
      do 140 j = 2, nxh
      sum1 = sum1 + f(1,2*j-1,i)
      f(1,2*j-1,i) = -f(1,2*j,i)
      f(1,2*j,i) = sum1
      sum2 = sum2 + f(2,2*j-1,i)
      f(2,2*j-1,i) = -f(2,2*j,i)
      f(2,2*j,i) = sum2
      sum3 = sum3 - f(3,2*j,i)
      f(3,2*j,i) = sum3
  140 continue
      f(1,nx+1,i) = 0.0
      f(2,nx+1,i) = 0.0
  150 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine PPFSSCT2RM23Y(g,isign,mixup,sctd,indx,indy,kstrt,kxpi, &
     &kxpp,nyv,kxpd,nxhyd,nxyd)
! this subroutine performs the y part of 3 two dimensional fast real
! sine transforms and their inverses, for a subset of x,
! using real arithmetic, with OpenMP,
! for data which is distributed in blocks
! x/y component has a sine transform, z component a cosine transform
! algorithm is described in Numerical Recipies in Fortran, Second Ed.,
! by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
! [Cambridge Univ. Press, 1992], p. 508.
! for isign = (-1,1), input: all, output: f
! approximate flop count: N*(5*log2(N) + 18)/nvp
! where N = (nx/2)*ny
! indx/indy = exponent which determines length in x/y direction,
! where nx=2**indx, ny=2**indy
! if isign = -1, inverse sine-sine-cosine transform are performed
! g(1:2,m,n) = sum(g(1:2,k,n)*sin(pi*m*k/ny))
! g(3,m,n) = (.5*g(3,1,n) + ((-1)**m)*g(3,ny+1,n)
!              + sum(g(3,k,n)*cos(pi*m*k/ny))
! if isign = 1, a forward sine-sine-cosine transforms are performed
! g(1:2,k,n) = sum(g(1:2,m,n)*sin(pi*m*k/ny))
! g(3,k,n) = 2*(.5*g(3,1,n) + ((-1)**m)*g(3,ny+1,n)
!              + sum(g(2,m,n)*cos(pi*m*k/ny))
! mixup = array of bit reversed addresses
! sctd = sine/cosine table
! kstrt = starting data block number
! kxpi = initial x index used
! kxpp = number of x indices used
! nyv = first dimension of g >= ny + 1
! kxpd = second dimension of g >= kxp + 1
! nxhyd = maximum of (nx/2,ny)
! nxyd = maximum of (nx,ny)
! written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, kstrt, kxpi, kxpp
      integer nyv, kxpd, nxhyd, nxyd, mixup
      real g
      complex sctd
      dimension g(3,nyv,kxpd)
      dimension mixup(nxhyd), sctd(nxyd)
! local data
      integer indx1, indy1, indx1y, nx, ny, nyh, nyhh, ny3, nxy, nxhy
      integer i, j, k, m, ks, km, kmr, nry, j1, j2, ns, ns2, k1, k2
      integer kxps, nryb, jj
      real at1, at2, at3, t2, t3, t4, t5, t6
      complex t1
      double precision sum1, sum2, sum3
      indx1 = indx - 1
      indy1 = indy - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      nyhh = ny/4
      ny3 = ny + 3
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 1
      kxps = kxpi + kxpp - 1
      if (kstrt.gt.nx) return
      if (isign.eq.0) return
      nryb = nxhy/nyh
      nry = nxy/nyh
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,m,jj,ns,ns2,km,kmr,k1,k2,j1,j2,at1,at2,at3,t2,t3,t4
!$OMP& ,t5,t6,t1,sum1,sum2,sum3)
      do 150 i = kxpi, kxps
! create auxiliary array in y
      kmr = nxy/ny
      sum3 = 0.5*(g(3,1,i) - g(3,ny+1,i))
      do 10 k = 2, nyh
      k1 = 1 + kmr*(k - 1)
      at3 = -aimag(sctd(k1))
      at2 = g(1,ny+2-k,i)
      at1 = g(1,k,i) + at2
      at2 = g(1,k,i) - at2
      at1 = at3*at1
      at2 = 0.5*at2
      g(1,k,i) = at1 + at2
      g(1,ny+2-k,i) = at1 - at2
      at2 = g(2,ny+2-k,i)
      at1 = g(2,k,i) + at2
      at2 = g(2,k,i) - at2
      at1 = at3*at1
      at2 = 0.5*at2
      g(2,k,i) = at1 + at2
      g(2,ny+2-k,i) = at1 - at2
      at2 = g(3,ny+2-k,i)
      at1 = g(3,k,i) + at2
      at2 = g(3,k,i) - at2
      sum3 = sum3 + real(sctd(k1))*at2
      at2 = at3*at2
      at1 = 0.5*at1
      g(3,k,i) = at1 - at2
      g(3,ny+2-k,i) = at1 + at2
   10 continue
      g(1,1,i) = 0.0
      g(1,nyh+1,i) = 2.0*g(1,nyh+1,i)
      g(2,1,i) = 0.0
      g(2,nyh+1,i) = 2.0*g(2,nyh+1,i)
      g(3,1,i) = 0.5*(g(3,1,i) + g(3,ny+1,i))
      g(3,ny+1,i) = sum3
! bit-reverse array elements in y
      do 30 k = 1, nyh
      k1 = (mixup(k) - 1)/nryb + 1
      if (k.lt.k1) then
         do 20 jj = 1, 3
         t2 = g(jj,2*k1-1,i)
         t3 = g(jj,2*k1,i)
         g(jj,2*k1-1,i) = g(jj,2*k-1,i)
         g(jj,2*k1,i) = g(jj,2*k,i)
         g(jj,2*k-1,i) = t2
         g(jj,2*k,i) = t3
   20    continue
      endif
   30 continue
! then transform in y
      do 70 m = 1, indy1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyhh/ns
      kmr = 2*km*nry
      do 60 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 50 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 40 jj = 1, 3
      t2 = real(t1)*g(jj,2*j2-1,i) - aimag(t1)*g(jj,2*j2,i)
      t3 = aimag(t1)*g(jj,2*j2-1,i) + real(t1)*g(jj,2*j2,i)
      g(jj,2*j2-1,i) = g(jj,2*j1-1,i) - t2
      g(jj,2*j2,i) = g(jj,2*j1,i) - t3
      g(jj,2*j1-1,i) = g(jj,2*j1-1,i) + t2
      g(jj,2*j1,i) = g(jj,2*j1,i) + t3
   40 continue
   50 continue
   60 continue
   70 continue
! unscramble coefficients and normalize
! inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nyh
         do 90 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         do 80 jj = 1, 3
         t4 = g(jj,ny3-2*k,i)
         t5 = -g(jj,ny3-2*k+1,i)
         t2 = g(jj,2*k-1,i) + t4
         t3 = g(jj,2*k,i) + t5
         t6 = g(jj,2*k-1,i) - t4
         t5 = g(jj,2*k,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         g(jj,2*k-1,i) = 0.5*(t2 + t4)
         g(jj,2*k,i) = 0.5*(t3 + t5)
         g(jj,ny3-2*k,i) = 0.5*(t2 - t4)
         g(jj,ny3-2*k+1,i) = 0.5*(t5 - t3)
   80    continue
   90    continue
         do 100 jj = 1, 3
         g(jj,nyh+1,i) = g(jj,nyh+1,i)
         g(jj,nyh+2,i) = -g(jj,nyh+2,i)
         t2 = g(jj,1,i) + g(jj,2,i)
         g(jj,2,i) = g(jj,1,i) - g(jj,2,i)
         g(jj,1,i) = t2
         g(jj,ny+1,i) = g(jj,ny+1,i)
  100    continue
! forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nyh
         do 120 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         do 110 jj = 1, 3
         t4 = g(jj,ny3-2*k,i)
         t5 = -g(jj,ny3-2*k+1,i)
         t2 = g(jj,2*k-1,i) + t4
         t3 = g(jj,2*k,i) + t5
         t6 = g(jj,2*k-1,i) - t4
         t5 = g(jj,2*k,i) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         g(jj,2*k-1,i) = t2 + t4
         g(jj,2*k,i) = t3 + t5
         g(jj,ny3-2*k,i) = t2 - t4
         g(jj,ny3-2*k+1,i) = t5 - t3
  110    continue
  120    continue
         do 130 jj = 1, 3
         g(jj,nyh+1,i) = 2.0*g(jj,nyh+1,i)
         g(jj,nyh+2,i) = -2.0*g(jj,nyh+2,i)
         t2 = 2.0*(g(jj,1,i) + g(jj,2,i))
         g(jj,2,i) = 2.0*(g(jj,1,i) - g(jj,2,i))
         g(jj,1,i) = t2
         g(jj,ny+1,i) = 2.0*g(jj,ny+1,i)
  130    continue
      endif
! perform recursion for sine-cosine transform
      sum1 = 0.5*g(1,1,i)
      g(1,1,i) = 0.0
      g(1,2,i) = sum1
      sum2 = 0.5*g(2,1,i)
      g(2,1,i) = 0.0
      g(2,2,i) = sum2
      sum3 = g(3,ny+1,i)
      g(3,ny+1,i) = g(3,2,i)
      g(3,2,i) = sum3
      do 140 k = 2, nyh
      sum1 = sum1 + g(1,2*k-1,i)
      g(1,2*k-1,i) = -g(1,2*k,i)
      g(1,2*k,i) = sum1
      sum2 = sum2 + g(2,2*k-1,i)
      g(2,2*k-1,i) = -g(2,2*k,i)
      g(2,2*k,i) = sum2
      sum3 = sum3 - g(3,2*k,i)
      g(3,2*k,i) = sum3
  140 continue
      g(1,ny+1,i) = 0.0
      g(2,ny+1,i) = 0.0
  150 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPDIVFD2(f,df,nx,ny,kstrt,ndim,nyv,kxp2)
! this subroutine calculates the divergence in fourier space
! with dirichlet boundary conditions (zero potential)
! using fast sine/cosine transforms for distributed data.
! intended for calculating the charge density from the electric field
! input: all except df, output: df
! approximate flop count is: 6*nx*ny
! the divergence is calculated using the equation:
! df(kx,ky) = -(kx*fx(kx,ky)+ky*fy(kx,ky))
! where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
! all for fourier mode (jj-1,k-1), where jj = j + kxp2*(kstrt - 1)
! modes nx and ny are zeroed out
! nx/ny = system length in x/y direction
! ndim = number of field arrays, must be >= 2
! kstrt = starting data block number
! nyv = first dimension of field arrays, must be >= ny+1
! kxp2 = number of data values per block
      implicit none
      integer nx, ny, kstrt, ndim, nyv ,kxp2
      real f, df
      dimension f(ndim,nyv,kxp2+1), df(nyv,kxp2+1)
! local data
      integer j, k, ks, ny1, joff, kxp2s
      real dnx, dny, dkx, dky
      if (ndim.lt.2) return
      ks = kstrt - 1
      ny1 = ny + 1
      joff = kxp2*ks
      kxp2s = min(kxp2,max(0,nx-joff))
      joff = joff - 1
      dnx = 6.28318530717959/real(nx + nx)
      dny = 6.28318530717959/real(ny + ny)
! calculate the divergence
      if (kstrt.gt.nx) return
! mode numbers 0 < kx < nx and 0 < ky < ny
!$OMP PARALLEL DO PRIVATE(j,k,dkx,dky)
      do 20 j = 1, kxp2s
      dkx = dnx*real(j + joff)
      if ((j+joff).gt.0) then
         do 10 k = 2, ny
         dky = dny*real(k - 1)
         df(k,j) = -(dkx*f(1,k,j) + dky*f(2,k,j))
   10    continue
      endif
! mode numbers ky = 0, ny
      df(1,j) = 0.0
      df(ny+1,j) = 0.0
   20 continue
!$OMP END PARALLEL DO
! mode numbers kx = 0, nx
      if (ks.eq.0) then
         do 30 k = 2, ny
         df(k,1) = 0.0
   30    continue
      endif
      do 40 k = 1, ny1
      df(k,kxp2s+1) = 0.0
   40 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPGRADFD2(df,f,nx,ny,kstrt,ndim,nyv,kxp2)
! this subroutine calculates the gradient in fourier space
! with dirichlet boundary conditions (zero potential)
! using fast sine/cosine transforms for distributed data.
! intended for calculating the electric field from the potential
! input: all except f, output: f
! approximate flop count is: 4*nx*ny
! the gradient is calculated using the equations:
! fx(kx,ky) = kx*df(kx,ky)
! fy(kx,ky) = ky*df(kx,ky)
! where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
! all for fourier mode (jj-1,k-1), where jj = j + kxp2*(kstrt - 1)
! modes nx and ny are zeroed out
! nx/ny = system length in x/y direction
! ndim = number of field arrays, must be >= 2
! kstrt = starting data block number
! nyv = first dimension of field arrays, must be >= ny+1
! kxp2 = number of data values per block
      implicit none
      integer nx, ny, kstrt, ndim, nyv, kxp2
      real df, f
      dimension df(nyv,kxp2+1), f(ndim,nyv,kxp2+1)
! local data
      integer j, k, ks, ny1, joff, kxp2s
      real dnx, dny, dkx, dky
      ks = kstrt - 1
      ny1 = ny + 1
      joff = kxp2*ks
      kxp2s = min(kxp2,max(0,nx-joff))
      joff = joff - 1
      dnx = 6.28318530717959/real(nx + nx)
      dny = 6.28318530717959/real(ny + ny)
! calculate the gradient
      if (kstrt.gt.nx) return
! mode numbers 0 < kx < nx and 0 < ky < ny
!$OMP PARALLEL DO PRIVATE(j,k,dkx,dky)
      do 20 j = 1, kxp2s
      dkx = dnx*real(j + joff)
      if ((j+joff).gt.0) then
         do 10 k = 2, ny
         dky = dny*real(k - 1)
         f(1,k,j) = dkx*df(k,j)
         f(2,k,j) = dky*df(k,j)
   10    continue
      endif
! mode numbers ky = 0, ny
      f(1,1,j) = 0.0
      f(2,1,j) = 0.0
      f(1,ny+1,j) = 0.0
      f(2,ny+1,j) = 0.0
   20 continue
!$OMP END PARALLEL DO
! mode numbers kx = 0, nx
      if (ks.eq.0) then
         do 30 k = 2, ny
         f(1,k,1) = 0.0
         f(2,k,1) = 0.0
   30    continue
      endif
      do 40 k = 1, ny1
      f(1,k,kxp2s+1) = 0.0
      f(2,k,kxp2s+1) = 0.0
   40 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPCURLFD2(f,g,nx,ny,kstrt,nyv,kxp2)
! this subroutine calculates the curl in fourier space
! with dirichlet boundary conditions (zero potential)
! using fast sine/cosine transforms for distributed data.
! intended for calculating the magnetic field from the vector potential
! input: all except g, output: g
! approximate flop count is: 8*nx*ny
! the curl is calculated using the equations:
! gx(kx,ky) = ky*fz(kx,ky)
! gy(kx,ky) = -kx*fz(kx,ky)
! gz(kx,ky) = (kx*fy(kx,ky)-ky*fx(kx,ky))
! where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
! all for fourier mode (jj-1,k-1), where jj = j + kxp2*(kstrt - 1)
! nx/ny = system length in x/y direction
! kstrt = starting data block number
! nyv = first dimension of field arrays, must be >= ny+1
! kxp2 = number of data values per block
      implicit none
      integer nx, ny, kstrt, nyv, kxp2
      real f, g
      dimension f(3,nyv,kxp2+1), g(3,nyv,kxp2+1)
! local data
      integer j, k, ks, ny1, joff, kxp2s
      real dnx, dny, dkx, dky
      ks = kstrt - 1
      ny1 = ny + 1
      joff = kxp2*ks
      kxp2s = min(kxp2,max(0,nx-joff))
      joff = joff - 1
      dnx = 6.28318530717959/real(nx + nx)
      dny = 6.28318530717959/real(ny + ny)
! calculate the curl
      if (kstrt.gt.nx) return
! mode numbers 0 < kx < nx and 0 < ky < ny
!$OMP PARALLEL DO PRIVATE(j,k,dkx,dky)
      do 20 j = 1, kxp2s
      dkx = dnx*real(j + joff)
      if ((j+joff).gt.0) then
         do 10 k = 2, ny
         dky = dny*real(k - 1)
         g(1,k,j) = dky*f(3,k,j)
         g(2,k,j) = -dkx*f(3,k,j)
         g(3,k,j) = dkx*f(2,k,j) - dky*f(1,k,j)
   10    continue
! mode numbers ky = 0, ny
         g(1,1,j) = 0.0
         g(2,1,j) = 0.0
         g(3,1,j) = dkx*f(2,1,j)
      endif
      g(1,ny+1,j) = 0.0
      g(2,ny+1,j) = 0.0
      g(3,ny+1,j) = 0.0
   20 continue
!$OMP END PARALLEL DO
! mode numbers kx = 0, nx
      if (ks.eq.0) then
         do 30 k = 2, ny
         dky = dny*real(k - 1)
         g(1,k,1) = 0.0
         g(2,k,1) = 0.0
         g(3,k,1) = -dky*f(1,k,1)
   30    continue
         g(1,1,1) = 0.0
         g(2,1,1) = 0.0
         g(3,1,1) = 0.0
      endif
      do 40 k = 1, ny1
      g(1,k,kxp2s+1) = 0.0
      g(2,k,kxp2s+1) = 0.0
      g(3,k,kxp2s+1) = 0.0
   40 continue
      return
      end            