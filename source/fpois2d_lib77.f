c-----------------------------------------------------------------------
      subroutine PPOISDX2(q,fx,fy,isign,ffd,ax,ay,affp,we,nx,ny,kstrt,ny
     12d,kxp2,j2blok,nyd)
c this subroutine solves 2d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c or for potential, or provides a smoothing function, with dirichlet
c boundary conditions (zero potential), for distributed data.
c fourier coefficients are constructed so that a real to complex fft
c will perform the appropriate sin-sin, sin-cos, or cos-sin transform
c for isign = 0,input: isign,ax,ay,affp,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd
c               output: ffd
c for isign = -1, input: q,ffd,isign,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd
c                output: fx,fy,we
c approximate flop count is: 11*nx*ny
c for isign = 1, input: q,ffd,isign,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd
c                output: fx,we
c approximate flop count is: 6*nx*ny
c for isign = 2, input: q,ffd,isign,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd
c                output: fy
c approximate flop count is: 2*nx*ny
c if isign < 0, force/charge is calculated using the equations:
c fx(kx,ky) = -sqrt(-1)*kx*g(kx,ky)*s(kx,ky)*q(kx,ky),
c fy(kx,ky) = -sqrt(-1)*ky*g(kx,ky)*s(kx,ky)*q(kx,ky),
c where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2)
c if isign = 1, potential is calculated using the equation:
c fx(kx,ky) = g(kx,ky)*q(kx,ky)*s(kx,ky)
c if isign = 2, smoothing is calculated using the equation:
c fy(kx,ky) = q(kx,ky)*s(kx,ky)
c q(k,j,l) = complex charge density for fourier mode (jj-1,k-1)
c fx(k,j,l) = x component of complex force/charge,
c fy(k,j,l) = y component of complex force/charge,
c all for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c if isign = 0, form factor array is prepared
c ffd(k,2*j,l) = finite-size particle shape factor s
c ffd(k,2*j-1,l) = potential green's function g
c all for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c j2blok = number of data blocks
c kxp2 = number of data values per block
c kstrt = starting data block number
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c electric field energy is also calculated, using
c we = 2*nx*ny*sum((affp/(kx**2+ky**2))*|q(kx,ky)*s(kx,ky)|**2)
c nx/ny = system length in x/y direction
c ny2d = first dimension of field arrays, must be >= 2*ny
c nyd = first dimension of form factor array, must be >= ny
      double precision wp
      complex q, fx, fy, ffd, zero
      dimension q(ny2d,kxp2,j2blok)
      dimension fx(ny2d,kxp2,j2blok), fy(ny2d,kxp2,j2blok)
      dimension ffd(nyd,kxp2,j2blok)
      ny2 = 2*ny + 2
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
      zero = cmplx(0.,0.)
      if (isign.ne.0) go to 40
      if (kstrt.gt.nx) return
c prepare form factor array
      do 30 l = 1, j2blok
      joff = kxp2*(l + ks) - 1
      do 20 j = 1, kxp2
      dkx = dnx*float(j + joff)
      at1 = dkx*dkx
      at2 = (dkx*ax)**2
      do 10 k = 1, ny
      dky = dny*float(k - 1)
      at3 = dky*dky + at1
      at4 = exp(-.5*((dky*ay)**2 + at2))
      if (at3.eq.0.) then
         ffd(k,j,l) = cmplx(affp,1.)
      else
         ffd(k,j,l) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
   30 continue
      return
   40 if (isign.gt.0) go to 100
c calculate force/charge and sum field energy
      wp = 0.0d0
      if (kstrt.gt.nx) go to 90
      do 80 l = 1, j2blok
c mode numbers kx > 0 and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 60 j = 1, kxp2
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 50 k = 2, ny
         k1 = ny2 - k
         at1 = real(ffd(k,j,l))*aimag(ffd(k,j,l))
         at3 = -at1*real(q(k,j,l))
         at2 = dkx*at3
         at3 = dny*float(k - 1)*at3
         fx(k,j,l) = cmplx(0.,at2)
         fx(k1,j,l) = cmplx(0.,-at2)
         fy(k,j,l) = cmplx(0.,at3)
         fy(k1,j,l) = cmplx(0.,at3)
         wp = wp + at1*real(q(k,j,l))**2
   50    continue
      endif
c mode numbers ky = 0, ny
      fx(1,j,l) = zero
      fx(ny+1,j,l) = zero
      fy(1,j,l) = zero
      fy(ny+1,j,l) = zero
   60 continue
c mode number kx = 0
      if ((l+ks).eq.0) then
         do 70 k = 2, ny
         k1 = ny2 - k
         fx(k,1,l) = zero
         fx(k1,1,l) = zero
         fy(k,1,l) = zero
         fy(k1,1,l) = zero
   70    continue
      endif
   80 continue
   90 continue
      we = 2.0*float(nx*ny)*wp
      return
c calculate potential and sum field energy
  100 if (isign.gt.1) go to 160
      wp = 0.0d0
      if (kstrt.gt.nx) go to 150
      do 140 l = 1, j2blok
c mode numbers kx > 0 and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 120 j = 1, kxp2
      if ((j+joff).gt.0) then
         do 110 k = 2, ny
         k1 = ny2 - k
         at2 = real(ffd(k,j,l))
         at1 = at2*aimag(ffd(k,j,l))
         at3 = at2*real(q(k,j,l))
         fx(k,j,l) = cmplx(at3,0.)
         fx(k1,j,l) = cmplx(-at3,0.)
         wp = wp + at1*real(q(k,j,l))**2
  110    continue
      endif
c mode numbers ky = 0, ny
      fx(1,j,l) = zero
      fx(ny+1,j,l) = zero
  120 continue
c mode number kx = 0
      if ((l+ks).eq.0) then
         do 130 k = 2, ny
         k1 = ny2 - k
         fx(k,1,l) = zero
         fx(k1,1,l) = zero
  130    continue
      endif
  140 continue
  150 continue
      we = 2.0*float(nx*ny)*wp
      return
c calculate smoothing
  160 if (kstrt.gt.nx) go to 210
      do 200 l = 1, j2blok
c mode numbers kx > 0 and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 180 j = 1, kxp2
      if ((j+joff).gt.0) then
         do 170 k = 2, ny
         k1 = ny2 - k
         at1 = aimag(ffd(k,j,l))
         at2 = at1*real(q(k,j,l))
         fy(k,j,l) = cmplx(at2,0.)
         fy(k1,j,l) = cmplx(-at2,0.)
  170    continue
      endif
c mode numbers ky = 0, ny
      fy(1,j,l) = zero
      fy(ny+1,j,l) = zero
  180 continue
c mode number kx = 0
      if ((l+ks).eq.0) then
         do 190 k = 2, ny
         k1 = ny2 - k
         fy(k,1,l) = zero
         fy(k1,1,l) = zero
  190    continue
      endif
  200 continue
  210 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PPOISD2(q,fx,fy,isign,ffd,ax,ay,affp,we,nx,ny,kstrt,ny
     1v,kxp2,j2blok,nyd)
c this subroutine solves 2d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c or for potential, or provides a smoothing function, with dirichlet
c boundary conditions (zero potential), for distributed data.
c fourier coefficients are constructed so that a real to complex fft
c will perform the appropriate sin-sin, sin-cos, or cos-sin transform
c for isign = 0,input: isign,ax,ay,affp,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd
c               output: ffd
c for isign = -1, input: q,ffd,isign,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd
c                output: fx,fy,we
c approximate flop count is: 10*nx*ny
c for isign = 1, input: q,ffd,isign,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd
c                output: fx,we
c approximate flop count is: 5*nx*ny
c for isign = 2, input: q,ffd,isign,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd
c                output: fy
c approximate flop count is: 1*nx*ny
c if isign < 0, force/charge is calculated using the equations:
c fx(kx,ky) = -sqrt(-1)*kx*g(kx,ky)*s(kx,ky)*q(kx,ky),
c fy(kx,ky) = -sqrt(-1)*ky*g(kx,ky)*s(kx,ky)*q(kx,ky),
c where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2)
c if isign = 1, potential is calculated using the equation:
c fx(kx,ky) = g(kx,ky)*q(kx,ky)*s(kx,ky)
c if isign = 2, smoothing is calculated using the equation:
c fy(kx,ky) = q(kx,ky)*s(kx,ky)
c q(k,j,l) = complex charge density for fourier mode (jj-1,k-1)
c fx(k,j,l) = x component of complex force/charge,
c fy(k,j,l) = y component of complex force/charge,
c all for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c if isign = 0, form factor array is prepared
c ffd(k,2*j,l) = finite-size particle shape factor s
c ffd(k,2*j-1,l) = potential green's function g
c all for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c j2blok = number of data blocks
c kxp2 = number of data values per block
c kstrt = starting data block number
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c electric field energy is also calculated, using
c we = 2*nx*ny*sum((affp/(kx**2+ky**2))*|q(kx,ky)*s(kx,ky)|**2)
c nx/ny = system length in x/y direction
c nyv = first dimension of field arrays, must be >= ny+1
c nyd = first dimension of form factor array, must be >= ny
      double precision wp
      complex ffd
      dimension q(nyv,kxp2+1,j2blok)
      dimension fx(nyv,kxp2+1,j2blok), fy(nyv,kxp2+1,j2blok)
      dimension ffd(nyd,kxp2,j2blok)
      ks = kstrt - 2
      ny1 = ny + 1
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
      if (isign.ne.0) go to 40
      if (kstrt.gt.nx) return
c prepare form factor array
      do 30 l = 1, j2blok
      joff = kxp2*(l + ks) - 1
      do 20 j = 1, kxp2
      dkx = dnx*float(j + joff)
      at1 = dkx*dkx
      at2 = (dkx*ax)**2
      do 10 k = 1, ny
      dky = dny*float(k - 1)
      at3 = dky*dky + at1
      at4 = exp(-.5*((dky*ay)**2 + at2))
      if (at3.eq.0.) then
         ffd(k,j,l) = cmplx(affp,1.)
      else
         ffd(k,j,l) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
   30 continue
      return
   40 if (isign.gt.0) go to 110
c calculate force/charge and sum field energy
      wp = 0.0d0
      if (kstrt.gt.nx) go to 100
      do 90 l = 1, j2blok
c mode numbers kx > 0 and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 60 j = 1, kxp2
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 50 k = 2, ny
         at1 = real(ffd(k,j,l))*aimag(ffd(k,j,l))
         at3 = -at1*q(k,j,l)
         at2 = dkx*at3
         at3 = dny*float(k - 1)*at3
         fx(k,j,l) = at2
         fy(k,j,l) = at3
         wp = wp + at1*q(k,j,l)**2
   50    continue
      endif
c mode numbers ky = 0, ny
      fx(1,j,l) = 0.
      fx(ny+1,j,l) = 0.
      fy(1,j,l) = 0.
      fy(ny+1,j,l) = 0.
   60 continue
c mode number kx = 0
      if ((l+ks).eq.0) then
         do 70 k = 2, ny
         fx(k,1,l) = 0.
         fy(k,1,l) = 0.
   70    continue
      endif
      do 80 k = 1, ny1
      fx(k,kxp2+1,l) = 0.
      fy(k,kxp2+1,l) = 0.
   80 continue
   90 continue
  100 continue
      we = 2.0*float(nx*ny)*wp
      return
c calculate potential and sum field energy
  110 if (isign.gt.1) go to 180
      wp = 0.0d0
      if (kstrt.gt.nx) go to 170
      do 160 l = 1, j2blok
c mode numbers kx > 0 and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 130 j = 1, kxp2
      if ((j+joff).gt.0) then
         do 120 k = 2, ny
         at2 = real(ffd(k,j,l))
         at1 = at2*aimag(ffd(k,j,l))
c         at3 = at2*q(k,j,l)
         at3 = at1*q(k,j,l)
         fx(k,j,l) = at3
         wp = wp + at1*q(k,j,l)**2
  120    continue
      endif
c mode numbers ky = 0, ny
      fx(1,j,l) = 0.
      fx(ny+1,j,l) = 0.
  130 continue
c mode number kx = 0
      if ((l+ks).eq.0) then
         do 140 k = 2, ny
         fx(k,1,l) = 0.
  140    continue
      endif
      do 150 k = 1, ny1
      fx(k,kxp2+1,l) = 0.
  150 continue
  160 continue
  170 continue
      we = 2.0*float(nx*ny)*wp
      return
c calculate smoothing
  180 if (kstrt.gt.nx) go to 240
      do 230 l = 1, j2blok
c mode numbers kx > 0 and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 200 j = 1, kxp2
      if ((j+joff).gt.0) then
         do 190 k = 2, ny
         at1 = aimag(ffd(k,j,l))
         at2 = at1*q(k,j,l)
         fy(k,j,l) = at2
  190    continue
      endif
c mode numbers ky = 0, ny
      fy(1,j,l) = 0.
      fy(ny+1,j,l) = 0.
  200 continue
c mode number kx = 0
      if ((l+ks).eq.0) then
         do 210 k = 2, ny
         fy(k,1,l) = 0.
  210    continue
      endif
      do 220 k = 1, ny1
      fy(k,kxp2+1,l) = 0.
  220 continue
  230 continue
  240 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PPOISD22(q,fxy,isign,ffd,ax,ay,affp,we,nx,ny,kstrt,nyv,
     1kxp2,j2blok,nyd)
c this subroutine solves 2d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c with dirichlet boundary conditions (zero potential),
c for distributed data.
c fourier coefficients are constructed so that a real to complex fft
c will perform the appropriate sin-cos, or cos-sin transform
c for isign = 0,input: isign,ax,ay,affp,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd
c               output: ffd
c for isign /= 0, input: q,ffd,isign,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd
c                output: fxy,we
c approximate flop count is: 10*nx*ny
c equation used is:
c fx(kx,ky) = -sqrt(-1)*kx*g(kx,ky)*s(kx,ky)*q(kx,ky),
c fy(kx,ky) = -sqrt(-1)*ky*g(kx,ky)*s(kx,ky)*q(kx,ky),
c where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2)
c q(k,j,l) = complex charge density for fourier mode (jj-1,k-1)
c fxy(1,k,j,l) = x component of complex force/charge,
c fxy(2,k,j,l) = y component of complex force/charge,
c all for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c if isign = 0, form factor array is prepared
c ffd(k,2*j,l) = finite-size particle shape factor s
c ffd(k,2*j-1,l) = potential green's function g
c all for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c j2blok = number of data blocks
c kxp2 = number of data values per block
c kstrt = starting data block number
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c electric field energy is also calculated, using
c we = 2*nx*ny*sum((affp/(kx**2+ky**2))*|q(kx,ky)*s(kx,ky)|**2)
c nx/ny = system length in x/y direction
c nyv = first dimension of field arrays, must be >= ny+1
c nyd = first dimension of form factor array, must be >= ny
      double precision wp
      complex ffd
      dimension q(nyv,kxp2+1,j2blok), fxy(2,nyv,kxp2+1,j2blok)
      dimension ffd(nyd,kxp2,j2blok)
      ks = kstrt - 2
      ny1 = ny + 1
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
      if (isign.ne.0) go to 40
      if (kstrt.gt.nx) return
c prepare form factor array
      do 30 l = 1, j2blok
      joff = kxp2*(l + ks) - 1
      do 20 j = 1, kxp2
      dkx = dnx*float(j + joff)
      at1 = dkx*dkx
      at2 = (dkx*ax)**2
      do 10 k = 1, ny
      dky = dny*float(k - 1)
      at3 = dky*dky + at1
      at4 = exp(-.5*((dky*ay)**2 + at2))
      if (at3.eq.0.) then
         ffd(k,j,l) = cmplx(affp,1.)
      else
         ffd(k,j,l) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
   30 continue
      return
c calculate force/charge and sum field energy
   40 wp = 0.0d0
      if (kstrt.gt.nx) go to 100
      do 90 l = 1, j2blok
c mode numbers kx > 0 and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 60 j = 1, kxp2
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 50 k = 2, ny
         at1 = real(ffd(k,j,l))*aimag(ffd(k,j,l))
         at3 = -at1*q(k,j,l)
         at2 = dkx*at3
         at3 = dny*float(k - 1)*at3
         fxy(1,k,j,l) = at2
         fxy(2,k,j,l) = at3
         wp = wp + at1*q(k,j,l)**2
   50    continue
      endif
c mode numbers ky = 0, ny
      fxy(1,1,j,l) = 0.
      fxy(2,1,j,l) = 0.
      fxy(1,ny+1,j,l) = 0.
      fxy(2,ny+1,j,l) = 0.
   60 continue
c mode number kx = 0
      if ((l+ks).eq.0) then
         do 70 k = 2, ny
         fxy(1,k,1,l) = 0.
         fxy(2,k,1,l) = 0.
   70    continue
      endif
      do 80 k = 1, ny1
      fxy(1,k,kxp2+1,l) = 0.
      fxy(2,k,kxp2+1,l) = 0.
   80 continue
   90 continue
  100 continue
      we = 2.0*float(nx*ny)*wp
      return
      end
c-----------------------------------------------------------------------
      subroutine PPOISDX23(q,fxy,isign,ffd,ax,ay,affp,we,nx,ny,kstrt,ny2
     1d,kxp2,j2blok,nyd)
c this subroutine solves 2d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c with dirichlet boundary conditions (zero potential),
c for distributed data.  Zeros out z component
c fourier coefficients are constructed so that a real to complex fft
c will perform the appropriate sin-cos, or cos-sin transform
c for isign = 0,input: isign,ax,ay,affp,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd
c               output: ffd
c for isign /= 0, input: q,ffd,isign,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd
c                output: fxy,we
c approximate flop count is: 11*nx*ny
c equation used is:
c fx(kx,ky) = -sqrt(-1)*kx*g(kx,ky)*s(kx,ky)*q(kx,ky),
c fy(kx,ky) = -sqrt(-1)*ky*g(kx,ky)*s(kx,ky)*q(kx,ky),
c fz(kx,ky) = zero,
c where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2)
c q(k,j,l) = complex charge density for fourier mode (jj-1,k-1)
c fxy(1,k,j,l) = x component of complex force/charge,
c fxy(2,k,j,l) = y component of complex force/charge,
c fxy(3,k,j,l) = zero,
c all for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c if isign = 0, form factor array is prepared
c ffd(k,2*j,l) = finite-size particle shape factor s
c ffd(k,2*j-1,l) = potential green's function g
c all for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c j2blok = number of data blocks
c kxp2 = number of data values per block
c kstrt = starting data block number
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c electric field energy is also calculated, using
c we = 2*nx*ny*sum((affp/(kx**2+ky**2))*|q(kx,ky)*s(kx,ky)|**2)
c nx/ny = system length in x/y direction
c ny2d = first dimension of field arrays, must be >= 2*ny
c nyd = first dimension of form factor array, must be >= ny
      double precision wp
      complex q, fxy, ffd, zero
      dimension q(ny2d,kxp2,j2blok), fxy(3,ny2d,kxp2,j2blok)
      dimension ffd(nyd,kxp2,j2blok)
      ny2 = 2*ny + 2
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
      zero = cmplx(0.,0.)
      if (isign.ne.0) go to 40
      if (kstrt.gt.nx) return
c prepare form factor array
      do 30 l = 1, j2blok
      joff = kxp2*(l + ks) - 1
      do 20 j = 1, kxp2
      dkx = dnx*float(j + joff)
      at1 = dkx*dkx
      at2 = (dkx*ax)**2
      do 10 k = 1, ny
      dky = dny*float(k - 1)
      at3 = dky*dky + at1
      at4 = exp(-.5*((dky*ay)**2 + at2))
      if (at3.eq.0.) then
         ffd(k,j,l) = cmplx(affp,1.)
      else
         ffd(k,j,l) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
   30 continue
      return
c calculate force/charge and sum field energy
   40 wp = 0.0d0
      if (kstrt.gt.nx) go to 90
      do 80 l = 1, j2blok
c mode numbers kx > 0 and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 60 j = 1, kxp2
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 50 k = 2, ny
         k1 = ny2 - k
         at1 = real(ffd(k,j,l))*aimag(ffd(k,j,l))
         at3 = -at1*real(q(k,j,l))
         at2 = dkx*at3
         at3 = dny*float(k - 1)*at3
         fxy(1,k,j,l) = cmplx(0.,at2)
         fxy(2,k,j,l) = cmplx(0.,at3)
         fxy(3,k,j,l) = zero
         fxy(1,k1,j,l) = cmplx(0.,-at2)
         fxy(2,k1,j,l) = cmplx(0.,at3)
         fxy(3,k1,j,l) = zero
         wp = wp + at1*real(q(k,j,l))**2
   50    continue
      endif
c mode numbers ky = 0, ny
      fxy(1,1,j,l) = zero
      fxy(2,1,j,l) = zero
      fxy(3,1,j,l) = zero
      fxy(1,ny+1,j,l) = zero
      fxy(2,ny+1,j,l) = zero
      fxy(3,ny+1,j,l) = zero
   60 continue
c mode number kx = 0
      if ((l+ks).eq.0) then
         do 70 k = 2, ny
         k1 = ny2 - k
         fxy(1,k,1,l) = zero
         fxy(2,k,1,l) = zero
         fxy(3,k,1,l) = zero
         fxy(1,k1,1,l) = zero
         fxy(2,k1,1,l) = zero
         fxy(3,k1,1,l) = zero
   70    continue
      endif
   80 continue
   90 continue
      we = 2.0*float(nx*ny)*wp
      return
      end
c-----------------------------------------------------------------------
      subroutine PPOISD23(q,fxy,isign,ffd,ax,ay,affp,we,nx,ny,kstrt,nyv,
     1kxp2,j2blok,nyd)
c this subroutine solves 2d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c with dirichlet boundary conditions (zero potential),
c for distributed data.  Zeros out z component
c fourier coefficients are constructed so that a real to complex fft
c will perform the appropriate sin-cos, or cos-sin transform
c for isign = 0,input: isign,ax,ay,affp,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd
c               output: ffd
c for isign /= 0, input: q,ffd,isign,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd
c                output: fxy,we
c approximate flop count is: 10*nx*ny
c equation used is:
c fx(kx,ky) = -sqrt(-1)*kx*g(kx,ky)*s(kx,ky)*q(kx,ky),
c fy(kx,ky) = -sqrt(-1)*ky*g(kx,ky)*s(kx,ky)*q(kx,ky),
c fz(kx,ky) = zero,
c where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2)
c q(k,j,l) = complex charge density for fourier mode (jj-1,k-1)
c fxy(1,k,j,l) = x component of complex force/charge,
c fxy(2,k,j,l) = y component of complex force/charge,
c fxy(3,k,j,l) = zero,
c all for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c if isign = 0, form factor array is prepared
c ffd(k,2*j,l) = finite-size particle shape factor s
c ffd(k,2*j-1,l) = potential green's function g
c all for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c j2blok = number of data blocks
c kxp2 = number of data values per block
c kstrt = starting data block number
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c electric field energy is also calculated, using
c we = 2*nx*ny*sum((affp/(kx**2+ky**2))*|q(kx,ky)*s(kx,ky)|**2)
c nx/ny = system length in x/y direction
c nyv = first dimension of field arrays, must be >= ny+1
c nyd = first dimension of form factor array, must be >= ny
      double precision wp
      complex ffd
      dimension q(nyv,kxp2+1,j2blok), fxy(3,nyv,kxp2+1,j2blok)
      dimension ffd(nyd,kxp2,j2blok)
      ks = kstrt - 2
      ny1 = ny + 1
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
      if (isign.ne.0) go to 40
      if (kstrt.gt.nx) return
c prepare form factor array
      do 30 l = 1, j2blok
      joff = kxp2*(l + ks) - 1
      do 20 j = 1, kxp2
      dkx = dnx*float(j + joff)
      at1 = dkx*dkx
      at2 = (dkx*ax)**2
      do 10 k = 1, ny
      dky = dny*float(k - 1)
      at3 = dky*dky + at1
      at4 = exp(-.5*((dky*ay)**2 + at2))
      if (at3.eq.0.) then
         ffd(k,j,l) = cmplx(affp,1.)
      else
         ffd(k,j,l) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
   30 continue
      return
c calculate force/charge and sum field energy
   40 wp = 0.0d0
      if (kstrt.gt.nx) go to 100
      do 90 l = 1, j2blok
c mode numbers kx > 0 and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 60 j = 1, kxp2
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 50 k = 2, ny
         at1 = real(ffd(k,j,l))*aimag(ffd(k,j,l))
         at3 = -at1*q(k,j,l)
         at2 = dkx*at3
         at3 = dny*float(k - 1)*at3
         fxy(1,k,j,l) = at2
         fxy(2,k,j,l) = at3
         fxy(3,k,j,l) = 0.
         wp = wp + at1*q(k,j,l)**2
   50    continue
      endif
c mode numbers ky = 0, ny
      fxy(1,1,j,l) = 0.
      fxy(2,1,j,l) = 0.
      fxy(3,1,j,l) = 0.
      fxy(1,ny+1,j,l) = 0.
      fxy(2,ny+1,j,l) = 0.
      fxy(3,ny+1,j,l) = 0.
   60 continue
c mode number kx = 0
      if ((l+ks).eq.0) then
         do 70 k = 2, ny
         fxy(1,k,1,l) = 0.
         fxy(2,k,1,l) = 0.
         fxy(3,k,1,l) = 0.
   70    continue
      endif
      do 80 k = 1, ny1
      fxy(1,k,kxp2+1,l) = 0.
      fxy(2,k,kxp2+1,l) = 0.
      fxy(3,k,kxp2+1,l) = 0.
   80 continue
   90 continue
  100 continue
      we = 2.0*float(nx*ny)*wp
      return
      end
c-----------------------------------------------------------------------
      subroutine PBPOISD22(cu,bxy,bz,isign,ffd,ax,ay,affp,ci,wm,nx,ny,ks
     1trt,nyv,kxp2,j2blok,nyd)
c this subroutine solves 2d poisson's equation in fourier space for
c magnetic field (or convolution of magnetic field over particle shape)
c or for vector potential, or provides a smoothing function,
c with dirichlet boundary conditions (zero potential),
c for distributed data.
c fourier coefficients are constructed so that a real to complex fft
c will perform the appropriate sin-cos, or cos-sin transform
c for isign = 0, input: isign,ax,ay,affp,nx,ny,kstrt,ny2d,kxp,j2blok,nyd
c output: ffd
c for isign = -1, input:cu,ffd,isign,ci,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd
c output: bz,wm
c approximate flop count is: 15*nxc*nyc + 7*(nxc + nyc)
c for isign = 1, input: cu,ffd,isign,ci,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd
c output: bxy,wm
c approximate flop count is: 10*nxc*nyc + 6*(nxc + nyc)
c for isign = 2, input: cu,ffd,isign,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd
c output: bxy
c approximate flop count is: 2*nxc*nyc + 1*(nxc + nyc)
c where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c if isign = 0, form factor array is prepared
c if isign < 0, magnetic field is calculated using the equations:
c bz(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*(kx*cuy(kx,ky)-ky*cux(kx,ky))*
c             s(kx,ky),
c where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2)
c if isign = 1, vector potential is calculated using the equation:
c bx(kx,ky) = ci*ci*g(kx,ky)*cux(kx,ky)
c by(kx,ky) = ci*ci*g(kx,ky)*cuy(kx,ky)
c if isign = 2, smoothing is calculated using the equation:
c bx(kx,ky) = cux(kx,ky)*s(kx,ky)
c by(kx,ky) = cuy(kx,ky)*s(kx,ky)
c cu(i,k,j,l) = i-th component of complex current density and
c bxy(i,k,j,l) = i-th component of complex vector potential,
c bz(k,j,l) = i-th component of complex magnetic field,
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c j2blok = number of data blocks
c kxp2 = number of data values per block
c kstrt = starting data block number
c aimag(ffd(k,j,l)) = finite-size particle shape factor s
c real(ffd(k,j,l)) = potential green's function g
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c ci = reciprical of velocity of light
c magnetic field energy is also calculated, using
c wm = nx*ny*nz*sum((affp/(kx**2+ky**2+kz**2))*ci*ci
c    |cu(kx,ky,kz)*s(kx,ky,kz)|**2)
c this expression is valid only if the current is divergence-free
c nx/ny = system length in x/y direction
c nyv = first dimension of field arrays, must be >= ny+1
c nyd = first dimension of form factor array, must be >= ny
      double precision wp
      complex ffd
      dimension cu(2,nyv,kxp2+1,j2blok), bxy(2,nyv,kxp2+1,j2blok)
      dimension bz(nyv,kxp2+1,j2blok)
      dimension ffd(nyd,kxp2,j2blok)
      ks = kstrt - 2
      ny1 = ny + 1
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
      ci2 = ci*ci
      if (isign.ne.0) go to 40
      if (kstrt.gt.nx) return
c prepare form factor array
      do 30 l = 1, j2blok
      joff = kxp2*(l + ks) - 1
      do 20 j = 1, kxp2
      dkx = dnx*float(j + joff)
      at1 = dkx*dkx
      at2 = (dkx*ax)**2
      do 10 k = 1, ny
      dky = dny*float(k - 1)
      at3 = dky*dky + at1
      at4 = exp(-.5*((dky*ay)**2 + at2))
      if (at3.eq.0.) then
         ffd(k,j,l) = cmplx(affp,1.)
      else
         ffd(k,j,l) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
   30 continue
      return
   40 if (isign.gt.0) go to 110
c calculate magnetic field and sum field energy
      wp = 0.0d0
      if (kstrt.gt.nx) go to 100
      do 90 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 60 j = 1, kxp2
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 50 k = 2, ny
         dky = dny*float(k - 1)
         at1 = ci2*real(ffd(k,j,l))*aimag(ffd(k,j,l))
         at2 = dky*at1
         at3 = dkx*at1
         bz(k,j,l) = at3*cu(2,k,j,l) - at2*cu(1,k,j,l)
         wp = wp + 2.0*at1*(cu(1,k,j,l)**2 + cu(2,k,j,l)**2)
   50    continue
c mode numbers ky = 0, ny
         at1 = ci2*real(ffd(1,j,l))*aimag(ffd(1,j,l))
         at2 = dkx*at1
         bz(1,j,l) = at2*cu(2,1,j,l)
         wp = wp + at1*cu(2,1,j,l)**2
      endif
      bz(ny+1,j,l) = 0.
   60 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 70 k = 2, ny
         dky = dny*float(k - 1)
         at1 = ci2*real(ffd(k,1,l))*aimag(ffd(k,1,l))
         at2 = dky*at1
         bz(k,1,l) = -at2*cu(1,k,1,l)
         wp = wp + at1*cu(1,k,1,l)**2
   70    continue
         bz(1,1,l) = 0.
      endif
      do 80 k = 1, ny1
      bz(k,kxp2+1,l) = 0.
   80 continue
   90 continue
  100 continue
      wm = float(nx*ny)*wp
      return
c calculate vector potential and sum field energy
  110 if (isign.gt.1) go to 180
      wp = 0.0d0
      if (kstrt.gt.nx) go to 170
      do 160 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 130 j = 1, kxp2
      if ((j+joff).gt.0) then
         do 120 k = 2, ny
         at2 = ci2*real(ffd(k,j,l))
         at1 = at2*aimag(ffd(k,j,l))
c         bxy(1,k,j,l) = at2*cu(1,k,j,l)
c         bxy(2,k,j,l) = at2*cu(2,k,j,l)
         bxy(1,k,j,l) = at1*cu(1,k,j,l)
         bxy(2,k,j,l) = at1*cu(2,k,j,l)
         wp = wp + 2.0*at1*(cu(1,k,j,l)**2 + cu(2,k,j,l)**2)
  120    continue
c mode numbers ky = 0, ny
         at2 = ci2*real(ffd(1,j,l))
         at1 = at2*aimag(ffd(1,j,l))
         bxy(1,1,j,l) = 0.
c         bxy(2,1,j,l) = at2*cu(2,1,j,l)
         bxy(2,1,j,l) = at1*cu(2,1,j,l)
         wp = wp + at1*cu(2,1,j,l)**2
      endif
      bxy(1,ny+1,j,l) = 0.
      bxy(2,ny+1,j,l) = 0.
  130 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 140 k = 2, ny
         at2 = ci2*real(ffd(k,1,l))
         at1 = at2*aimag(ffd(k,1,l))
c         bxy(1,k,1,l) = at2*cu(1,k,1,l)
         bxy(1,k,1,l) = at1*cu(1,k,1,l)
         bxy(2,k,1,l) = 0.
         wp = wp + at1*cu(1,k,1,l)**2
  140    continue
         bxy(1,1,1,l) = 0.
         bxy(2,1,1,l) = 0.
      endif
      do 150 k = 1, ny1
      bxy(1,k,kxp2+1,l) = 0.
      bxy(2,k,kxp2+1,l) = 0.
  150 continue
  160 continue
  170 continue
      wm = float(nx*ny)*wp
      return
c calculate smoothing
  180 if (kstrt.gt.nx) go to 240
      do 230 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 200 j = 1, kxp2
      if ((j+joff).gt.0) then
         do 190 k = 2, ny
         at1 = aimag(ffd(k,j,l))
         bxy(1,k,j,l) = at1*cu(1,k,j,l)
         bxy(2,k,j,l) = at1*cu(2,k,j,l)
  190    continue
c mode numbers ky = 0, ny
         at1 = aimag(ffd(1,j,l))
         bxy(1,1,j,l) = 0.
         bxy(2,1,j,l) = at1*cu(2,1,j,l)
      endif
      bxy(1,ny+1,j,l) = 0.
      bxy(2,ny+1,j,l) = 0.
  200 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 210 k = 2, ny
         at1 = aimag(ffd(k,1,l))
         bxy(1,k,1,l) = at1*cu(1,k,1,l)
         bxy(2,k,1,l) = 0.
  210    continue
         bxy(1,1,1,l) = 0.
         bxy(2,1,1,l) = 0.
      endif
      do 220 k = 1, ny1
      bxy(1,k,kxp2+1,l) = 0.
      bxy(2,k,kxp2+1,l) = 0.
  220 continue
  230 continue
  240 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PBPOISD23(cu,bxy,isign,ffd,ax,ay,affp,ci,wm,nx,ny,kstrt
     1,nyv,kxp2,j2blok,nyd)
c this subroutine solves 2-1/2d poisson's equation in fourier space for
c magnetic field (or convolution of magnetic field over particle shape)
c or for vector potential, or provides a smoothing function,
c with dirichlet boundary conditions (zero potential),
c for distributed data.
c fourier coefficients are constructed so that a real to complex fft
c will perform the appropriate sin-cos, or cos-sin transform
c for isign = 0, input: isign,ax,ay,affp,nx,ny,kstrt,ny2d,kxp,j2blok,nyd
c output: ffd
c for isign = -1, input:cu,ffd,isign,ci,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd
c output: bxy,wm
c approximate flop count is: 20*nxc*nyc + 8*(nxc + nyc)
c for isign = 1, input: cu,ffd,isign,ci,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd
c output: bxy,wm
c approximate flop count is: 13*nxc*nyc + 8*(nxc + nyc)
c for isign = 2, input: cu,ffd,isign,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd
c output: bxy
c approximate flop count is: 3*nxc*nyc + 1*(nxc + nyc)
c where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c if isign = 0, form factor array is prepared
c if isign < 0, magnetic field is calculated using the equations:
c bx(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*ky*cuz(kx,ky)*s(kx,ky),
c by(kx,ky) = -ci*ci*sqrt(-1)*g(kx,ky)*kx*cuz(kx,ky)*s(kx,ky),
c bz(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*(kx*cuy(kx,ky)-ky*cux(kx,ky))*
c             s(kx,ky),
c where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2)
c if isign = 1, vector potential is calculated using the equation:
c bx(kx,ky) = ci*ci*g(kx,ky)*cux(kx,ky)
c by(kx,ky) = ci*ci*g(kx,ky)*cuy(kx,ky)
c bz(kx,ky) = ci*ci*g(kx,ky)*cuz(kx,ky)
c if isign = 2, smoothing is calculated using the equation:
c bx(kx,ky) = cux(kx,ky)*s(kx,ky)
c by(kx,ky) = cuy(kx,ky)*s(kx,ky)
c bz(kx,ky) = cuz(kx,ky)*s(kx,ky)
c cu(i,k,j,l) = i-th component of complex current density and
c bxy(i,k,j,l) = i-th component of complex magnetic field,
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c j2blok = number of data blocks
c kxp2 = number of data values per block
c kstrt = starting data block number
c aimag(ffd(k,j,l)) = finite-size particle shape factor s
c real(ffd(k,j,l)) = potential green's function g
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c ci = reciprical of velocity of light
c magnetic field energy is also calculated, using
c wm = nx*ny*nz*sum((affp/(kx**2+ky**2+kz**2))*ci*ci
c    |cu(kx,ky,kz)*s(kx,ky,kz)|**2)
c this expression is valid only if the current is divergence-free
c nx/ny = system length in x/y direction
c nyv = first dimension of field arrays, must be >= ny+1
c nyd = first dimension of form factor array, must be >= ny
      double precision wp
      complex ffd
      dimension cu(3,nyv,kxp2+1,j2blok), bxy(3,nyv,kxp2+1,j2blok)
      dimension ffd(nyd,kxp2,j2blok)
      ks = kstrt - 2
      ny1 = ny + 1
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
      ci2 = ci*ci
      if (isign.ne.0) go to 40
      if (kstrt.gt.nx) return
c prepare form factor array
      do 30 l = 1, j2blok
      joff = kxp2*(l + ks) - 1
      do 20 j = 1, kxp2
      dkx = dnx*float(j + joff)
      at1 = dkx*dkx
      at2 = (dkx*ax)**2
      do 10 k = 1, ny
      dky = dny*float(k - 1)
      at3 = dky*dky + at1
      at4 = exp(-.5*((dky*ay)**2 + at2))
      if (at3.eq.0.) then
         ffd(k,j,l) = cmplx(affp,1.)
      else
         ffd(k,j,l) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
   30 continue
      return
   40 if (isign.gt.0) go to 110
c calculate magnetic field and sum field energy
      wp = 0.0d0
      if (kstrt.gt.nx) go to 100
      do 90 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 60 j = 1, kxp2
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 50 k = 2, ny
         dky = dny*float(k - 1)
         at1 = ci2*real(ffd(k,j,l))*aimag(ffd(k,j,l))
         at2 = dky*at1
         at3 = dkx*at1
         bxy(1,k,j,l) = at2*cu(3,k,j,l)
         bxy(2,k,j,l) = -at3*cu(3,k,j,l)
         bxy(3,k,j,l) = at3*cu(2,k,j,l) - at2*cu(1,k,j,l)
         wp = wp + 2.0*at1*(cu(1,k,j,l)**2 + cu(2,k,j,l)**2 + cu(3,k,j,l
     1)**2)
   50    continue
c mode numbers ky = 0, ny
         at1 = ci2*real(ffd(1,j,l))*aimag(ffd(1,j,l))
         at2 = dkx*at1
         bxy(1,1,j,l) = 0.
         bxy(2,1,j,l) = 0.
         bxy(3,1,j,l) = at2*cu(2,1,j,l)
         wp = wp + at1*(cu(2,1,j,l)**2 + cu(3,1,j,l)**2)
      endif
      bxy(1,ny+1,j,l) = 0.
      bxy(2,ny+1,j,l) = 0.
      bxy(3,ny+1,j,l) = 0.
   60 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 70 k = 2, ny
         dky = dny*float(k - 1)
         at1 = ci2*real(ffd(k,1,l))*aimag(ffd(k,1,l))
         at2 = dky*at1
         bxy(1,k,1,l) = 0.
         bxy(2,k,1,l) = 0.
         bxy(3,k,1,l) = -at2*cu(1,k,1,l)
         wp = wp + at1*(cu(1,k,1,l)**2 + cu(3,k,1,l)**2)
   70    continue
         bxy(1,1,1,l) = 0.
         bxy(2,1,1,l) = 0.
         bxy(3,1,1,l) = 0.
      endif
      do 80 k = 1, ny1
      bxy(1,k,kxp2+1,l) = 0.
      bxy(2,k,kxp2+1,l) = 0.
      bxy(3,k,kxp2+1,l) = 0.
   80 continue
   90 continue
  100 continue
      wm = float(nx*ny)*wp
      return
c calculate vector potential and sum field energy
  110 if (isign.gt.1) go to 180
      wp = 0.0d0
      if (kstrt.gt.nx) go to 170
      do 160 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 130 j = 1, kxp2
      if ((j+joff).gt.0) then
         do 120 k = 2, ny
         at2 = ci2*real(ffd(k,j,l))
         at1 = at2*aimag(ffd(k,j,l))
c         bxy(1,k,j,l) = at2*cu(1,k,j,l)
c         bxy(2,k,j,l) = at2*cu(2,k,j,l)
c         bxy(3,k,j,l) = at2*cu(3,k,j,l)
         bxy(1,k,j,l) = at1*cu(1,k,j,l)
         bxy(2,k,j,l) = at1*cu(2,k,j,l)
         bxy(3,k,j,l) = at1*cu(3,k,j,l)
         wp = wp + 2.0*at1*(cu(1,k,j,l)**2 + cu(2,k,j,l)**2 + cu(3,k,j,l
     1)**2)
  120    continue
c mode numbers ky = 0, ny
         at2 = ci2*real(ffd(1,j,l))
         at1 = at2*aimag(ffd(1,j,l))
         bxy(1,1,j,l) = 0.
c         bxy(2,1,j,l) = at2*cu(2,1,j,l)
         bxy(2,1,j,l) = at1*cu(2,1,j,l)
         bxy(3,1,j,l) = 0.
         wp = wp + at1*(cu(2,1,j,l)**2 + cu(3,1,j,l)**2)
      endif
      bxy(1,ny+1,j,l) = 0.
      bxy(2,ny+1,j,l) = 0.
      bxy(3,ny+1,j,l) = 0.
  130 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 140 k = 2, ny
         at2 = ci2*real(ffd(k,1,l))
         at1 = at2*aimag(ffd(k,1,l))
c         bxy(1,k,1,l) = at2*cu(1,k,1,l)
         bxy(1,k,1,l) = at1*cu(1,k,1,l)
         bxy(2,k,1,l) = 0.
         bxy(3,k,1,l) = 0.
         wp = wp + at1*(cu(1,k,1,l)**2 + cu(3,k,1,l)**2)
  140    continue
         bxy(1,1,1,l) = 0.
         bxy(2,1,1,l) = 0.
         bxy(3,1,1,l) = 0.
      endif
      do 150 k = 1, ny1
      bxy(1,k,kxp2+1,l) = 0.
      bxy(2,k,kxp2+1,l) = 0.
      bxy(3,k,kxp2+1,l) = 0.
  150 continue
  160 continue
  170 continue
      wm = float(nx*ny)*wp
      return
c calculate smoothing
  180 if (kstrt.gt.nx) go to 240
      do 230 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 200 j = 1, kxp2
      if ((j+joff).gt.0) then
         do 190 k = 2, ny
         at1 = aimag(ffd(k,j,l))
         bxy(1,k,j,l) = at1*cu(1,k,j,l)
         bxy(2,k,j,l) = at1*cu(2,k,j,l)
         bxy(3,k,j,l) = at1*cu(3,k,j,l)
  190    continue
c mode numbers ky = 0, ny
         at1 = aimag(ffd(1,j,l))
         bxy(1,1,j,l) = 0.
         bxy(2,1,j,l) = at1*cu(2,1,j,l)
         bxy(3,1,j,l) = 0.
      endif
      bxy(1,ny+1,j,l) = 0.
      bxy(2,ny+1,j,l) = 0.
      bxy(3,ny+1,j,l) = 0.
  200 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 210 k = 2, ny
         at1 = aimag(ffd(k,1,l))
         bxy(1,k,1,l) = at1*cu(1,k,1,l)
         bxy(2,k,1,l) = 0.
         bxy(3,k,1,l) = 0.
  210    continue
         bxy(1,1,1,l) = 0.
         bxy(2,1,1,l) = 0.
         bxy(3,1,1,l) = 0.
      endif
      do 220 k = 1, ny1
      bxy(1,k,kxp2+1,l) = 0.
      bxy(2,k,kxp2+1,l) = 0.
      bxy(3,k,kxp2+1,l) = 0.
  220 continue
  230 continue
  240 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PBPOISD22N_QP(cu,dcu,amu,bxy,bz,isign,ffd,ax,ay,affp,ci
     1,wm,nx,ny,kstrt,nyv,kxp2,j2blok,nyd,aa,dex)
c this subroutine solves 2d poisson's equation in fourier space for
c magnetic field (or convolution of magnetic field over particle shape)
c or for vector potential, or provides a smoothing function,
c with dirichlet boundary conditions (zero potential),
c for distributed data.
c fourier coefficients are constructed so that a real to complex fft
c will perform the appropriate sin-cos, or cos-sin transform
c for isign = 0, input: isign,ax,ay,affp,nx,ny,kstrt,ny2d,kxp,j2blok,nyd
c output: ffd
c for isign = -1, input:cu,ffd,isign,ci,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd
c output: bz,wm
c approximate flop count is: 15*nxc*nyc + 7*(nxc + nyc)
c for isign = 1, input: cu,ffd,isign,ci,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd
c output: bxy,wm
c approximate flop count is: 10*nxc*nyc + 6*(nxc + nyc)
c for isign = 2, input: cu,ffd,isign,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd
c output: bxy
c approximate flop count is: 2*nxc*nyc + 1*(nxc + nyc)
c where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c if isign = 0, form factor array is prepared
c if isign < 0, magnetic field is calculated using the equations:
c bz(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*(kx*cuy(kx,ky)-ky*cux(kx,ky))*
c             s(kx,ky),
c where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2)
c if isign = 1, vector potential is calculated using the equation:
c bx(kx,ky) = ci*ci*g(kx,ky)*cux(kx,ky)
c by(kx,ky) = ci*ci*g(kx,ky)*cuy(kx,ky)
c if isign = 2, smoothing is calculated using the equation:
c bx(kx,ky) = cux(kx,ky)*s(kx,ky)
c by(kx,ky) = cuy(kx,ky)*s(kx,ky)
c cu(i,k,j,l) = i-th component of complex current density and
c bxy(i,k,j,l) = i-th component of complex vector potential,
c bz(k,j,l) = i-th component of complex magnetic field,
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c j2blok = number of data blocks
c kxp2 = number of data values per block
c kstrt = starting data block number
c aimag(ffd(k,j,l)) = finite-size particle shape factor s
c real(ffd(k,j,l)) = potential green's function g
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c ci = reciprical of velocity of light
c magnetic field energy is also calculated, using
c wm = nx*ny*nz*sum((affp/(kx**2+ky**2+kz**2))*ci*ci
c    |cu(kx,ky,kz)*s(kx,ky,kz)|**2)
c this expression is valid only if the current is divergence-free
c nx/ny = system length in x/y direction
c nyv = first dimension of field arrays, must be >= ny+1
c nyd = first dimension of form factor array, must be >= ny
      double precision wp
      complex ffd
      real cu, dcu, amu, dex,idex
      dimension cu(3,nyv,kxp2+1,j2blok), bxy(2,nyv,kxp2+1,j2blok)
      dimension dcu(2,nyv,kxp2+1,j2blok), amu(3,nyv,kxp2+1,j2blok)
      dimension bz(nyv,kxp2+1,j2blok)
      dimension ffd(nyd,kxp2,j2blok)
      idex = 1.0/dex
      ks = kstrt - 2
      ny1 = ny + 1
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
      ci2 = ci*ci
      if (isign.ne.0) go to 40
      if (kstrt.gt.nx) return
c prepare form factor array
      do 30 l = 1, j2blok
      joff = kxp2*(l + ks) - 1
      do 20 j = 1, kxp2
      dkx = dnx*float(j + joff)
      at1 = dkx*dkx
      at2 = (dkx*ax)**2
      do 10 k = 1, ny
      dky = dny*float(k - 1)
      at3 = dky*dky + at1
      at4 = exp(-.5*((dky*ay)**2 + at2))
      if (at3.eq.0.) then
         ffd(k,j,l) = cmplx(affp,1.)
      else
         ffd(k,j,l) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
   30 continue
      return
   40 if (isign.gt.0) go to 110
c calculate magnetic field and sum field energy
      wp = 0.0d0
      if (kstrt.gt.nx) go to 100
      do 90 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 60 j = 1, kxp2
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 50 k = 2, ny
         dky = dny*float(k - 1)
         at1 = ci2*real(ffd(k,j,l))*aimag(ffd(k,j,l))
         at2 = dky*at1
         at3 = dkx*at1
         bz(k,j,l) = at3*cu(2,k,j,l) - at2*cu(1,k,j,l)
         wp = wp + 2.0*at1*(cu(1,k,j,l)**2 + cu(2,k,j,l)**2)
   50    continue
c mode numbers ky = 0, ny
         at1 = ci2*real(ffd(1,j,l))*aimag(ffd(1,j,l))
         at2 = dkx*at1
         bz(1,j,l) = at2*cu(2,1,j,l)
         wp = wp + at1*cu(2,1,j,l)**2
      endif
      bz(ny+1,j,l) = 0.
   60 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 70 k = 2, ny
         dky = dny*float(k - 1)
         at1 = ci2*real(ffd(k,1,l))*aimag(ffd(k,1,l))
         at2 = dky*at1
         bz(k,1,l) = -at2*cu(1,k,1,l)
         wp = wp + at1*cu(1,k,1,l)**2
   70    continue
         bz(1,1,l) = 0.
      endif
      do 80 k = 1, ny1
      bz(k,kxp2+1,l) = 0.
   80 continue
   90 continue
  100 continue
      wm = float(nx*ny)*wp
      return
c calculate vector potential and sum field energy
  110 if (isign.gt.1) go to 180
      wp = 0.0d0
      if (kstrt.gt.nx) go to 170
      do 160 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 130 j = 1, kxp2
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 120 k = 2, ny
         dky = dny*float(k - 1)
         dcu(1,k,j,l) = -dcu(1,k,j,l)+dkx*amu(1,k,j,l)-dky*amu(3,k,j,l)
         aat1 = dcu(1,k,j,l) - dkx*cu(3,k,j,l)*idex         
         dcu(2,k,j,l) = -dcu(2,k,j,l)-dkx*amu(3,k,j,l)+dky*amu(2,k,j,l)
         aat2 = dcu(2,k,j,l) - dky*cu(3,k,j,l)*idex         
c-------------------------------------------------------
c         at2 = real(ffd(k,j,l))
c         at2 = ci2*at2/(1+aa*at2)
c         at1 = at2*aimag(ffd(k,j,l))
c         bxy(1,k,j,l) = at2*cu(1,k,j,l)
c         bxy(2,k,j,l) = at2*cu(2,k,j,l)
c-------------------------------------------------------
c         at2 = real(ffd(k,j,l))
c         at1 = at2*aimag(ffd(k,j,l))
c         at2 = at2/aimag(ffd(k,j,l))*aa
c         bxy(1,k,j,l) = (at1*cu(1,k,j,l)+at2*bxy(1,k,j,l))/(1+at2)
c         bxy(2,k,j,l) = (at1*cu(2,k,j,l)+at2*bxy(2,k,j,l))/(1+at2)
c-------------------------------------------------------
         at2 = real(ffd(k,j,l))
         at2 = ci2*at2/(1+aa*at2)
         at1 = at2*aimag(ffd(k,j,l))
         bxy(1,k,j,l) = at1*aat1+aa*at2*bxy(1,k,j,l)
         bxy(2,k,j,l) = at1*aat2+aa*at2*bxy(2,k,j,l)
         wp = wp + 2.0*at1*(dcu(1,k,j,l)**2 + dcu(2,k,j,l)**2)
  120    continue
c mode numbers ky = 0, ny
c-----------------------------------------------------------------
c         at2 = real(ffd(k,j,l))
c         at2 = ci2*at2/(1+aa*at2)
c         at1 = at2*aimag(ffd(1,j,l))
c         bxy(1,1,j,l) = 0.
c         bxy(2,1,j,l) = at2*cu(2,1,j,l)
c-----------------------------------------------------------------
c         at2 = real(ffd(k,j,l))
c         at1 = at2*aimag(ffd(k,j,l))
c         at2 = at2/aimag(ffd(k,j,l))*aa
c         bxy(1,1,j,l) = 0.
c         bxy(2,1,j,l) = (at1*cu(2,1,j,l)+at2*bxy(2,1,j,l))/(1+at2)
c-----------------------------------------------------------------
         dcu(1,1,j,l) = -dcu(1,1,j,l)+dkx*amu(1,1,j,l)
         dcu(2,1,j,l) = -dcu(2,1,j,l)-dkx*amu(3,1,j,l)         

         at2 = real(ffd(1,j,l))
         at2 = ci2*at2/(1+aa*at2)
         at1 = at2*aimag(ffd(1,j,l))
         bxy(1,1,j,l) = 0.
         bxy(2,1,j,l) = at1*dcu(2,1,j,l)+aa*at2*bxy(2,1,j,l)
         wp = wp + at1*dcu(2,1,j,l)**2
      endif
      bxy(1,ny+1,j,l) = 0.
      bxy(2,ny+1,j,l) = 0.
  130 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 140 k = 2, ny
c-----------------------------------------------------------------
c         at2 = real(ffd(k,j,l))
c         at2 = ci2*at2/(1+aa*at2)
c         at1 = at2*aimag(ffd(k,1,l))
c         bxy(1,k,1,l) = at2*cu(1,k,1,l)
c         bxy(2,k,1,l) = 0.
c-----------------------------------------------------------------
c         at2 = real(ffd(k,j,l))
c         at1 = at2*aimag(ffd(k,j,l))
c         at2 = at2/aimag(ffd(k,j,l))*aa
c         bxy(1,k,1,l) = (at1*cu(1,k,1,l)+at2*bxy(1,k,1,l))/(1+at2)
c         bxy(2,k,1,l) = 0.
c-----------------------------------------------------------------
         dky = dny*float(k - 1)
         dcu(1,k,1,l) = -dcu(1,k,1,l)-dky*amu(3,k,1,l)
         dcu(2,k,1,l) = -dcu(2,k,1,l)+dky*amu(2,k,1,l)

         at2 = real(ffd(k,1,l))
         at2 = ci2*at2/(1+aa*at2)
         at1 = at2*aimag(ffd(k,1,l))
         bxy(1,k,1,l) = at1*dcu(1,k,1,l)+aa*at2*bxy(1,k,1,l)
         bxy(2,k,1,l) = 0.
         wp = wp + at1*dcu(1,k,1,l)**2
  140    continue
         bxy(1,1,1,l) = 0.
         bxy(2,1,1,l) = 0.
      endif
      do 150 k = 1, ny1
      bxy(1,k,kxp2+1,l) = 0.
      bxy(2,k,kxp2+1,l) = 0.
  150 continue
  160 continue
  170 continue
      wm = float(nx*ny)*wp
      return
c calculate smoothing
  180 if (kstrt.gt.nx) go to 240
      do 230 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 200 j = 1, kxp2
      if ((j+joff).gt.0) then
         do 190 k = 2, ny
         at1 = aimag(ffd(k,j,l))
         bxy(1,k,j,l) = at1*cu(1,k,j,l)
         bxy(2,k,j,l) = at1*cu(2,k,j,l)
  190    continue
c mode numbers ky = 0, ny
         at1 = aimag(ffd(1,j,l))
         bxy(1,1,j,l) = 0.
         bxy(2,1,j,l) = at1*cu(2,1,j,l)
      endif
      bxy(1,ny+1,j,l) = 0.
      bxy(2,ny+1,j,l) = 0.
  200 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 210 k = 2, ny
         at1 = aimag(ffd(k,1,l))
         bxy(1,k,1,l) = at1*cu(1,k,1,l)
         bxy(2,k,1,l) = 0.
  210    continue
         bxy(1,1,1,l) = 0.
         bxy(2,1,1,l) = 0.
      endif
      do 220 k = 1, ny1
      bxy(1,k,kxp2+1,l) = 0.
      bxy(2,k,kxp2+1,l) = 0.
  220 continue
  230 continue
  240 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPOISD22(q,fxy,isign,ffd,ax,ay,affp,we,nx,ny,kstrt,nyv&
     &,kxp2,nyd)
! this subroutine solves 2d poisson's equation in fourier space for
! force/charge (or convolution of electric field over particle shape)
! with dirichlet boundary conditions (zero potential),
! using fast sine/cosine transforms for distributed data.
! for isign = 0,input: isign,ax,ay,affp,nx,ny,kstrt,ny2d,kxp2,nyd
!               output: ffd
! for isign /= 0, input: q,ffd,isign,nx,ny,kstrt,ny2d,kxp2,nyd
!                 output: fxy,we
! approximate flop count is: 10*nx*ny
! equation used is:
! fx(kx,ky) = -kx*g(kx,ky)*s(kx,ky)*q(kx,ky),
! fy(kx,ky) = -ky*g(kx,ky)*s(kx,ky)*q(kx,ky),
! where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
! modes nx and ny are zeroed out
! g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
! s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2)
! q(k,j) = transformed charge density for fourier mode (jj-1,k-1)
! fxy(1,k,j) = x component of transformed force/charge,
! fxy(2,k,j) = y component of transformed force/charge,
! all for fourier mode (jj-1,k-1), where jj = j + kxp2*(kstrt - 1)
! if isign = 0, form factor array is prepared
! aimag(ffd(k,j)= finite-size particle shape factor s
! real(ffd(k,j)) = potential green's function g
! all for fourier mode (jj-1,k-1), where jj = j + kxp2*(kstrt - 1)
! kxp2 = number of data values per block
! kstrt = starting data block number
! ax/ay = half-width of particle in x/y direction
! affp = normalization constant = nx*ny/np, where np=number of particles
! electric field energy is also calculated, using
! we = 2*nx*ny*sum((affp/(kx**2+ky**2))*|q(kx,ky)*s(kx,ky)|**2)
! nx/ny = system length in x/y direction
! nyv = first dimension of field arrays, must be >= ny+1
! nyd = first dimension of form factor array, must be >= ny
      implicit none
      integer isign, nx, ny, kstrt, nyv, kxp2, nyd
      real ax, ay, affp, we
      real q, fxy
      dimension q(nyv,kxp2+1), fxy(2,nyv,kxp2+1)
      complex ffd
      dimension ffd(nyd,kxp2)
! local data
      integer j, k, ks, ny1, joff, kxp2s
      real dnx, dny, dkx, dky, at1, at2, at3, at4
      double precision wp, sum1
      ks = kstrt - 1
      ny1 = ny + 1
      joff = kxp2*ks
      kxp2s = min(kxp2,max(0,nx-joff))
      joff = joff - 1
      dnx = 6.28318530717959/real(nx + nx)
      dny = 6.28318530717959/real(ny + ny)
      if (isign.ne.0) go to 30
! prepare form factor array
      do 20 j = 1, kxp2s
      dkx = dnx*real(j + joff)
      at1 = dkx*dkx
      at2 = (dkx*ax)**2
      do 10 k = 1, ny
      dky = dny*real(k - 1)
      at3 = dky*dky + at1
      at4 = exp(-0.5*((dky*ay)**2 + at2))
      if (at3.eq.0.0) then
         ffd(k,j) = cmplx(affp,1.0)
      else
         ffd(k,j) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
      return
! calculate force/charge and sum field energy
   30 sum1 = 0.0d0
      if (kstrt.gt.nx) go to 80
! mode numbers kx > 0 and 0 < ky < ny
!$OMP PARALLEL DO PRIVATE(j,k,dkx,at1,at2,at3,wp)
!$OMP& REDUCTION(+:sum1)
      do 50 j = 1, kxp2s
      dkx = dnx*real(j + joff)
      wp = 0.0d0
      if ((j+joff).gt.0) then
         do 40 k = 2, ny
         at1 = real(ffd(k,j))*aimag(ffd(k,j))
         at3 = -at1*q(k,j)
         at2 = dkx*at3
         at3 = dny*real(k - 1)*at3
         fxy(1,k,j) = at2
         fxy(2,k,j) = at3
         wp = wp + at1*q(k,j)**2
   40    continue
      endif
! mode numbers ky = 0, ny
      fxy(1,1,j) = 0.0
      fxy(2,1,j) = 0.0
      fxy(1,ny+1,j) = 0.0
      fxy(2,ny+1,j) = 0.0
      sum1 = sum1 + wp
   50 continue
!$OMP END PARALLEL DO
! mode number kx = 0
      if (ks.eq.0) then
         do 60 k = 2, ny
         fxy(1,k,1) = 0.0
         fxy(2,k,1) = 0.0
   60    continue
      endif
! zero out kx = nx mode and unused extra cells
      do 70 k = 1, ny1
      fxy(1,k,kxp2s+1) = 0.0
      fxy(2,k,kxp2s+1) = 0.0
   70 continue
   80 continue
      we = 2.0*real(nx)*real(ny)*sum1
      return
      end   
!-----------------------------------------------------------------------
      subroutine MPPOISD23(q,fxy,isign,ffd,ax,ay,affp,we,nx,ny,kstrt,nyv&
     &,kxp2,nyd)
! this subroutine solves 2-1/2d poisson's equation in fourier space for
! force/charge (or convolution of electric field over particle shape)
! with dirichlet boundary conditions (zero potential),
! using fast sine/cosine transforms for distributed data.
! Zeros out z component
! for isign = 0,input: isign,ax,ay,affp,nx,ny,kstrt,ny2d,kxp2,nyd
!               output: ffd
! for isign /= 0, input: q,ffd,isign,nx,ny,kstrt,ny2d,kxp2,nyd
!                 output: fxy,we
! approximate flop count is: 10*nx*ny
! equation used is:
! fx(kx,ky) = -kx*g(kx,ky)*s(kx,ky)*q(kx,ky),
! fy(kx,ky) = -ky*g(kx,ky)*s(kx,ky)*q(kx,ky),
! fz(kx,ky) = zero,
! where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
! modes nx and ny are zeroed out
! g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
! s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2)
! q(k,j) = transformed charge density for fourier mode (jj-1,k-1)
! fxy(1,k,j) = x component of transformed force/charge,
! fxy(2,k,j) = y component of transformed force/charge,
! fxy(3,k,j) = zero,
! all for fourier mode (jj-1,k-1), where jj = j + kxp2*(kstrt - 1)
! if isign = 0, form factor array is prepared
! aimag(ffd(k,j)= finite-size particle shape factor s
! real(ffd(k,j)) = potential green's function g
! all for fourier mode (jj-1,k-1), where jj = j + kxp2*(kstrt - 1)
! kxp2 = number of data values per block
! kstrt = starting data block number
! ax/ay = half-width of particle in x/y direction
! affp = normalization constant = nx*ny/np, where np=number of particles
! electric field energy is also calculated, using
! we = 2*nx*ny*sum((affp/(kx**2+ky**2))*|q(kx,ky)*s(kx,ky)|**2)
! nx/ny = system length in x/y direction
! nyv = first dimension of field arrays, must be >= ny+1
! nyd = first dimension of form factor array, must be >= ny
      implicit none
      integer isign, nx, ny, kstrt, nyv, kxp2, nyd
      real ax, ay, affp, we
      real q, fxy
      dimension q(nyv,kxp2+1), fxy(3,nyv,kxp2+1)
      complex ffd
      dimension ffd(nyd,kxp2)
! local data
      integer j, k, ks, ny1, joff, kxp2s
      real dnx, dny, dkx, dky, at1, at2, at3, at4
      double precision wp, sum1
      ks = kstrt - 1
      ny1 = ny + 1
      joff = kxp2*ks
      kxp2s = min(kxp2,max(0,nx-joff))
      joff = joff - 1
      dnx = 6.28318530717959/real(nx + nx)
      dny = 6.28318530717959/real(ny + ny)
      if (isign.ne.0) go to 30
! prepare form factor array
      do 20 j = 1, kxp2s
      dkx = dnx*real(j + joff)
      at1 = dkx*dkx
      at2 = (dkx*ax)**2
      do 10 k = 1, ny
      dky = dny*real(k - 1)
      at3 = dky*dky + at1
      at4 = exp(-0.5*((dky*ay)**2 + at2))
      if (at3.eq.0.0) then
         ffd(k,j) = cmplx(affp,1.0)
      else
         ffd(k,j) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
      return
! calculate force/charge and sum field energy
   30 sum1 = 0.0d0
      if (kstrt.gt.nx) go to 80
! mode numbers kx > 0 and 0 < ky < ny
!$OMP PARALLEL DO PRIVATE(j,k,dkx,at1,at2,at3,wp)
!$OMP& REDUCTION(+:sum1)
      do 50 j = 1, kxp2s
      dkx = dnx*real(j + joff)
      wp = 0.0d0
      if ((j+joff).gt.0) then
         do 40 k = 2, ny
         at1 = real(ffd(k,j))*aimag(ffd(k,j))
         at3 = -at1*q(k,j)
         at2 = dkx*at3
         at3 = dny*real(k - 1)*at3
         fxy(1,k,j) = at2
         fxy(2,k,j) = at3
         fxy(3,k,j) = 0.0
         wp = wp + at1*q(k,j)**2
   40    continue
      endif
! mode numbers ky = 0, ny
      fxy(1,1,j) = 0.0
      fxy(2,1,j) = 0.0
      fxy(3,1,j) = 0.0
      fxy(1,ny+1,j) = 0.0
      fxy(2,ny+1,j) = 0.0
      fxy(3,ny+1,j) = 0.0
      sum1 = sum1 + wp
   50 continue
!$OMP END PARALLEL DO
! mode number kx = 0
      if (ks.eq.0) then
         do 60 k = 2, ny
         fxy(1,k,1) = 0.0
         fxy(2,k,1) = 0.0
         fxy(3,k,1) = 0.0
   60    continue
      endif
! zero out kx = nx mode and unused extra cells
      do 70 k = 1, ny1
      fxy(1,k,kxp2s+1) = 0.0
      fxy(2,k,kxp2s+1) = 0.0
      fxy(3,k,kxp2s+1) = 0.0
   70 continue
   80 continue
      we = 2.0*real(nx)*real(ny)*sum1
      return
      end         
!-----------------------------------------------------------------------
      subroutine MPPOTPD2(q,pot,ffd,we,nx,ny,kstrt,nyv,kxp2,nyd)
! this subroutine solves 2d poisson's equation in fourier space for
! potential, with dirichlet boundary conditions (zero potential),
! using fast sine transforms for distributed data.
! input: q,ffd,nx,ny,kstrt,nyv,kxp2,nyd, output: pot,we
! approximate flop count is: 5*nx*ny
! potential is calculated using the equation:
! fx(kx,ky) = g(kx,ky)*q(kx,ky)*s(kx,ky)
! where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
! modes nx and ny are zeroed out
! g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
! s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2)
! q(k,j) = transformed charge density for fourier mode (jj-1,k-1)
! pot(k,j) = transformed potential,
! all for fourier mode (jj-1,k-1), where jj = j + kxp2*(kstrt - 1)
! kxp2 = number of data values per block
! kstrt = starting data block number
! aimag(ffd(k,j)= finite-size particle shape factor s
! real(ffd(k,j)) = potential green's function g
! electric field energy is also calculated, using
! we = 2*nx*ny*sum((affp/(kx**2+ky**2))*|q(kx,ky)*s(kx,ky)|**2)
! where affp = normalization constant = nx*ny/np,
! where np=number of particles
! nx/ny = system length in x/y direction
! nyv = second dimension of field arrays, must be >= ny+1
! nyd = first dimension of form factor array, must be >= ny
      implicit none
      integer nx, ny, kstrt, nyv, kxp2, nyd
      real we
      real q, pot
      dimension q(nyv,kxp2+1), pot(nyv,kxp2+1)
      complex ffd
      dimension ffd(nyd,kxp2)
! local data
      integer j, k, ks, ny1, joff, kxp2s
      real at1, at2, at3
      double precision wp, sum1
      ks = kstrt - 1
      ny1 = ny + 1
      joff = kxp2*ks
      kxp2s = min(kxp2,max(0,nx-joff))
      joff = joff - 1
! calculate potential and sum field energy
      sum1 = 0.0d0
      if (kstrt.gt.nx) go to 50
! mode numbers kx > 0 and 0 < ky < ny
!$OMP PARALLEL DO PRIVATE(j,k,at1,at2,at3,wp)
!$OMP& REDUCTION(+:sum1)
      do 20 j = 1, kxp2s
      wp = 0.0d0
      if ((j+joff).gt.0) then
         do 10 k = 2, ny
         at2 = real(ffd(k,j))
         at1 = at2*aimag(ffd(k,j))
         at3 = at2*q(k,j)
         pot(k,j) = at3
         wp = wp + at1*q(k,j)**2
  10     continue
      endif
! mode numbers ky = 0, ny
      pot(1,j) = 0.0
      pot(ny+1,j) = 0.0
      sum1 = sum1 + wp
   20 continue
!$OMP END PARALLEL DO
! mode number kx = 0
      if (ks.eq.0) then
         do 30 k = 2, ny
         pot(k,1) = 0.0
   30    continue
      endif
! zero out kx = nx mode and unused extra cells
      do 40 k = 1, ny1
      pot(k,kxp2s+1) = 0.0
   40 continue
   50 continue
      we = 2.0*real(nx)*real(ny)*sum1
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPSMOOTHD2(q,qs,ffd,nx,ny,kstrt,nyv,kxp2,nyd)
! this subroutine provides a 2d scalar smoothing function
! in fourier space, with dirichlet boundary conditions (zero potential),
! using fast sine transforms for distributed data.
! input: q,ffd,nx,ny,kstrt,nyv,kxp2,nyd, output: qs
! approximate flop count is: 1*nx*ny
! smoothing is calculated using the equation:
! qs(kx,ky) = q(kx,ky)*s(kx,ky)
! where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
! modes nx and ny are zeroed out
! s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2)
! q(k,j) = transformed charge density for fourier mode (jj-1,k-1)
! qs(k,j) = transformed smoothed charge density,
! all for fourier mode (jj-1,k-1), where jj = j + kxp2*(kstrt - 1)
! kxp2 = number of data values per block
! kstrt = starting data block number
! aimag(ffd(k,j)) = finite-size particle shape factor s
! nx/ny = system length in x/y direction
! nyv = first dimension of field arrays, must be >= ny+1
! nyd = first dimension of form factor array, must be >= ny
      implicit none
      integer nx, ny, kstrt, nyv, kxp2, nyd
      real q, qs
      dimension q(nyv,kxp2+1), qs(nyv,kxp2+1)
      complex ffd
      dimension ffd(nyd,kxp2)
! local data
      integer j, k, ks, ny1, joff, kxp2s
      real at1, at2
      ks = kstrt - 1
      ny1 = ny + 1
      joff = kxp2*ks
      kxp2s = min(kxp2,max(0,nx-joff))
      joff = joff - 1
! calculate smoothing
      if (kstrt.gt.nx) return
! mode numbers kx > 0 and 0 < ky < ny
!$OMP PARALLEL DO PRIVATE(j,k,at1,at2)
      do 20 j = 1, kxp2s
      if ((j+joff).gt.0) then
         do 10 k = 2, ny
         at1 = aimag(ffd(k,j))
         at2 = at1*q(k,j)
         qs(k,j) = at2
  10     continue
      endif
! mode numbers ky = 0, ny
      qs(1,j) = 0.0
      qs(ny+1,j) = 0.0
   20 continue
!$OMP END PARALLEL DO
! mode number kx = 0
      if (ks.eq.0) then
         do 30 k = 2, ny
         qs(k,1) = 0.0
   30    continue
      endif
! zero out kx = nx mode and unused extra cells
      do 40 k = 1, ny1
      qs(k,kxp2s+1) = 0.0
   40 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPBBPOISD23(cu,bxy,ffd,ci,wm,nx,ny,kstrt,nyv,kxp2,nyd)
! this subroutine solves 2-1/2d poisson's equation in fourier space for
! smoothed magnetic field (or convolution of magnetic field over 
! particle shape) with dirichlet boundary conditions (zero potential),
! using fast sine/cosine transforms for distributed data.
! input: cu,ffd,ci,nx,ny,kstrt,ny2d,kxp2,nyd, output: bxy,wm
! approximate flop count is: 20*nx*ny
! magnetic field is calculated using the equations:
! bx(kx,ky) = ci*ci**g(kx,ky)*ky*cuz(kx,ky)*s(kx,ky),
! by(kx,ky) = -ci*ci*sg(kx,ky)*kx*cuz(kx,ky)*s(kx,ky),
! bz(kx,ky) = ci*ci*g(kx,ky)*(kx*cuy(kx,ky)-ky*cux(kx,ky))*s(kx,ky),
! where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
! modes nx and ny are zeroed out
! g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
! s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2)
! cu(i,k,j) = i-th component of transformed current density and
! bxy(i,k,j) = i-th component of transformed magnetic field,
! all for fourier mode (jj-1,k-1), where jj = j + kxp2*(kstrt - 1)
! aimag(ffd(k,j)) = finite-size particle shape factor s
! real(ffd(k,j)) = potential green's function g
! all for fourier mode (jj-1,k-1), where jj = j + kxp2*(kstrt - 1)
! kxp2 = number of data values per block
! kstrt = starting data block number=
! ci = reciprocal of velocity of light
! magnetic field energy is also calculated, using
! wm = nx*ny*nz*sum((affp/(kx**2+ky**2+kz**2))*ci*ci
!    |cu(kx,ky,kz)*s(kx,ky,kz)|**2)
! where affp = normalization constant = nx*ny/np,
! where np=number of particles
! this expression is valid only if the current is divergence-free
! nx/ny = system length in x/y direction
! nyv = second dimension of field arrays, must be >= ny+1
! nyd = first dimension of form factor array, must be >= ny
      implicit none
      integer nx, ny, kstrt, nyv, kxp2, nyd
      real ci, wm
      real cu, bxy
      dimension cu(3,nyv,kxp2+1), bxy(3,nyv,kxp2+1)
      complex ffd
      dimension ffd(nyd,kxp2)
! local data
      integer j, k, ks, ny1, joff, kxp2s
      real dnx, dny, dkx, dky, ci2, at1, at2, at3
      double precision wp, sum1
      ks = kstrt - 1
      ny1 = ny + 1
      joff = kxp2*ks
      kxp2s = min(kxp2,max(0,nx-joff))
      joff = joff - 1
      dnx = 6.28318530717959/real(nx + nx)
      dny = 6.28318530717959/real(ny + ny)
      ci2 = ci*ci
! calculate smoothed magnetic field and sum field energy
      sum1 = 0.0d0
      if (kstrt.gt.nx) go to 50
! mode numbers kx > 0 and 0 < ky < ny
!$OMP PARALLEL DO PRIVATE(j,k,dkx,dky,at1,at2,at3,wp)
!$OMP& REDUCTION(+:sum1)
      do 20 j = 1, kxp2s
      dkx = dnx*real(j + joff)
      wp = 0.0d0
      if ((j+joff).gt.0) then
         do 10 k = 2, ny
         dky = dny*real(k - 1)
         at1 = ci2*real(ffd(k,j))*aimag(ffd(k,j))
         at2 = dky*at1
         at3 = dkx*at1
         bxy(1,k,j) = at2*cu(3,k,j)
         bxy(2,k,j) = -at3*cu(3,k,j)
         bxy(3,k,j) = at3*cu(2,k,j) - at2*cu(1,k,j)
         wp = wp + 2.0*at1*(cu(1,k,j)**2 + cu(2,k,j)**2 + cu(3,k,j)**2)
   10    continue
! mode numbers ky = 0, ny
         at1 = ci2*real(ffd(1,j))*aimag(ffd(1,j))
         at2 = dkx*at1
         bxy(1,1,j) = 0.0
         bxy(2,1,j) = 0.0
         bxy(3,1,j) = at2*cu(2,1,j)
         wp = wp + at1*(cu(2,1,j)**2 + cu(3,1,j)**2)
      endif
      bxy(1,ny+1,j) = 0.0
      bxy(2,ny+1,j) = 0.0
      bxy(3,ny+1,j) = 0.0
      sum1 = sum1 + wp
   20 continue
!$OMP END PARALLEL DO
      wp = 0.0d0
! mode numbers kx = 0, nx
      if (ks.eq.0) then
         do 30 k = 2, ny
         dky = dny*real(k - 1)
         at1 = ci2*real(ffd(k,1))*aimag(ffd(k,1))
         at2 = dky*at1
         bxy(1,k,1) = 0.0
         bxy(2,k,1) = 0.0
         bxy(3,k,1) = -at2*cu(1,k,1)
         wp = wp + at1*(cu(1,k,1)**2 + cu(3,k,1)**2)
   30    continue
         bxy(1,1,1) = 0.0
         bxy(2,1,1) = 0.0
         bxy(3,1,1) = 0.0
      endif
      sum1 = sum1 + wp
! zero out kx = nx mode and unused extra cells
      do 40 k = 1, ny1
      bxy(1,k,kxp2s+1) = 0.0
      bxy(2,k,kxp2s+1) = 0.0
      bxy(3,k,kxp2s+1) = 0.0
   40 continue
   50 continue
      wm = real(nx)*real(ny)*sum1
      return
      end      