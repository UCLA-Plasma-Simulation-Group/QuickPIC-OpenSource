module integral
    use perrors_class
    use parallel_pipe_class

    class(perrors), pointer, public :: err => null()
!         class(spect3d), pointer, public :: sp => null()
    class(parallel_pipe), pointer, public :: p => null()
abstract interface 
    subroutine event(k1,k2,dt,re,im)
    implicit none
    REAL(8),  INTENT(inout) :: re,im
    real(8),  INTENT(IN) :: k1, k2
    real(8),  INTENT(IN):: dt
    ! real(8) :: dtt
   end subroutine event
end interface

contains

subroutine calculate_inte0(k1,k2,dt,re,im)
    implicit none
    REAL(8),  INTENT(inout) :: re,im
    real(8),  INTENT(IN) :: k1, k2
    real(8),  INTENT(IN):: dt
    real(8) :: dtt
    re = dt
    im = 0
    end subroutine

subroutine calculate_inte1(k1,k2,dt,re,im)
    implicit none
    REAL(8),  INTENT(inout) :: re,im
    real(8),  INTENT(IN) :: k1, k2
    real(8),  INTENT(IN):: dt
    real(8) :: dtt
    call calculate_inte_temp1(k1,k2,dt,re,im)
    end subroutine

subroutine calculate_inte2(k1,k2,dt,re,im)
    implicit none
    REAL(8),  INTENT(inout) :: re,im
    real(8),  INTENT(IN) :: k1, k2
    real(8),  INTENT(IN):: dt
    real(8) :: dtt
    dtt =0.5*dt
    re = 2*sin(k1*dtt)/k1
    im = k2*(2*dtt**2*sin(k1*dtt)/k1 + 4*dtt*cos(k1*dtt)/k1**2 - 4*sin(k1*dtt)/k1**3)
    end subroutine

subroutine calculate_inte_temp1(k1,k2,dt,re,im)
    implicit none
    real(8), INTENT(IN) :: k1, k2
    real(8) ,INTENT(IN) :: dt
    real(8), INTENT(INOUT)  :: re,im

    real(8) :: dtt,pi,xi1,xi2,u1,u2,v1,v2,d1,d2,temp_c1,temp_c2
    real(8) :: co,si,temp1_re,temp1_im,temp2_re,temp2_im,xi11,xi22
    Logical :: flag
    pi=3.141592653589793d0
    dtt =0.5*dt

    if (k2<0) then 
        d1 = k1
        d2 = -k2

        temp_c1 = (2*pi/d2)**0.5/2
        temp_c2 = d1**2/4/d2
        co = cos(temp_c2)
        si = sin(temp_c2)
        temp1_re = temp_c1*(co-si)
        temp1_im = temp_c1*(co+si)
        ! print*,"temp1:",temp1_re,temp1_im

        xi1 =d1/d2**0.5
        xi2 = 2*d2**0.5*dtt
        xi11 = (xi1-xi2)/2.8284271247461903d0
        xi22 = (xi1+xi2)/2.8284271247461903d0
        call erfc_c(xi11, xi11, u1, v1, flag)
        call erfc_c(xi22, xi22, u2, v2, flag)

    else 
        d1 = k1
        d2 = k2

        temp_c1 = (2*pi/d2)**0.5/2
        temp_c2 = -d1**2/(4*d2)
        co = cos(temp_c2)
        si = sin(temp_c2)
        temp1_re = temp_c1*(co+si)
        temp1_im = temp_c1*(si-co)

        
        
        xi1 =d1/d2**0.5
        xi2 = 2*d2**0.5*dtt
        xi11 = (xi1+xi2)/2.8284271247461903d0
        xi22 = (xi1-xi2)/2.8284271247461903d0            
        call erfc_c(xi11, -xi11, u1, v1, flag)
        call erfc_c(xi22, -xi22, u2, v2, flag)
    end if
        temp2_re = (v1-v2)/2
        temp2_im = (u2-u1)/2
        ! print*,"temp2:",temp2_re,temp2_im
        ! print*,u1,v1
        ! print*,temp2_re,temp2_im
        re = temp1_re*temp2_re - temp1_im*temp2_im
        im = temp1_re*temp2_im + temp1_im*temp2_re 
    
    end subroutine 
subroutine calculate_inte_temp2(k1,k2,dt,re,im)
    implicit none
    REAL(8),  INTENT(inout) :: re,im
    real(8),  INTENT(IN) :: k1, k2
    real(8),  INTENT(IN):: dt
    real(8) :: dtt
    dtt =0.5*dt
    re = 2*sin(k1*dtt)/k1
    im = k2*(2*dtt**2*sin(k1*dtt)/k1 + 4*dtt*cos(k1*dtt)/k1**2 - 4*sin(k1*dtt)/k1**3)
    end SUBROUTINE 

subroutine erfc_c(xi,yi,u,v,flag)
  real(8),INTENT(in) ::xi,yi
  real(8),INTENT(inout) :: u,v
  real(8) :: wx,wy,temp1,temp2_co,temp2_si
  Logical,INTENT(inout):: flag
  call wofz(-yi,xi,wx,wy,flag)
  temp1 = exp(yi**2-xi**2)
  temp2_co = cos(2*xi*yi)
  temp2_si = sin(2*xi*yi)
  u = (wx*temp2_co+wy*temp2_si)*temp1
  v = (wy*temp2_co-wx*temp2_si)*temp1
end subroutine


Subroutine wofz(xi, yi, u, v, flag)
  Implicit Double Precision (A-H, O-Z)

  Logical a, b, flag
  integer :: kapn,n,np1,nu
  Parameter (factor=1.12837916709551257388D0, rmaxreal=0.5D+154, rmaxexp=708.503061461606D0, rmaxgoni=3.53711887601422D+15)
!
  flag = .False.
!
  xabs = dabs(xi)
  yabs = dabs(yi)
  x = xabs/6.3
  y = yabs/4.4
!
!
!     THE FOLLOWING IF-STATEMENT PROTECTS
!     QRHO = (X**2 + Y**2) AGAINST OVERFLOW
!
  If ((xabs>rmaxreal) .Or. (yabs>rmaxreal)) Goto 100
!
  qrho = x**2 + y**2
!
  xabsq = xabs**2
  xquad = xabsq - yabs**2
  yquad = 2*xabs*yabs
!
  a = qrho < 0.085264D0
!
  If (a) Then
!
!  IF (QRHO.LT.0.085264D0) THEN THE FADDEEVA-FUNCTION IS EVALUATED
!  USING A POWER-SERIES (ABRAMOWITZ/STEGUN, EQUATION (7.1.5), P.297)
!  N IS THE MINIMUM NUMBER OF TERMS NEEDED TO OBTAIN THE REQUIRED
!  ACCURACY
!
    qrho = (1-0.85*y)*dsqrt(qrho)
    n = idnint(6+72*qrho)
    j = 2*n + 1
    xsum = 1.0/j
    ysum = 0.0D0
    Do i = n, 1, -1
      j = j - 2
      xaux = (xsum*xquad-ysum*yquad)/i
      ysum = (xsum*yquad+ysum*xquad)/i
      xsum = xaux + 1.0/j
    End Do
    u1 = -factor*(xsum*yabs+ysum*xabs) + 1.0
    v1 = factor*(xsum*xabs-ysum*yabs)
    daux = dexp(-xquad)
    u2 = daux*dcos(yquad)
    v2 = -daux*dsin(yquad)
!
    u = u1*u2 - v1*v2
    v = u1*v2 + v1*u2
!
  Else
!
!  IF (QRHO.GT.1.O) THEN W(Z) IS EVALUATED USING THE LAPLACE
!  CONTINUED FRACTION
!  NU IS THE MINIMUM NUMBER OF TERMS NEEDED TO OBTAIN THE REQUIRED
!  ACCURACY
!
!  IF ((QRHO.GT.0.085264D0).AND.(QRHO.LT.1.0)) THEN W(Z) IS EVALUATED
!  BY A TRUNCATED TAYLOR EXPANSION, WHERE THE LAPLACE CONTINUED FRACTION
!  IS USED TO CALCULATE THE DERIVATIVES OF W(Z)
!  KAPN IS THE MINIMUM NUMBER OF TERMS IN THE TAYLOR EXPANSION NEEDED
!  TO OBTAIN THE REQUIRED ACCURACY
!  NU IS THE MINIMUM NUMBER OF TERMS OF THE CONTINUED FRACTION NEEDED
!  TO CALCULATE THE DERIVATIVES WITH THE REQUIRED ACCURACY
!
!
    If (qrho>1.0) Then
      h = 0.0D0
      kapn = 0
      qrho = dsqrt(qrho)
      nu = idint(3+(1442/(26*qrho+77)))
    Else
      qrho = (1-y)*dsqrt(1-qrho)
      h = 1.88*qrho
      h2 = 2*h
      kapn = idnint(7+34*qrho)
      nu = idnint(16+26*qrho)
    End If
!
    b = (h>0.0)
!
    If (b) qlambda = h2**kapn
!
    rx = 0.0
    ry = 0.0
    sx = 0.0
    sy = 0.0
!
    Do n = nu, 0, -1
      np1 = n + 1
      tx = yabs + h + np1*rx
      ty = xabs - np1*ry
      c = 0.5/(tx**2+ty**2)
      rx = c*tx
      ry = c*ty
      If ((b) .And. (n<=kapn)) Then
        tx = qlambda + sx
        sx = rx*tx - ry*sy
        sy = ry*tx + rx*sy
        qlambda = qlambda/h2
      End If
    End Do
!
    If (h==0.0) Then
      u = factor*rx
      v = factor*ry
    Else
      u = factor*sx
      v = factor*sy
    End If
!
    If (yabs==0.0) u = dexp(-xabs**2)
!
  End If
!
!
!
!  EVALUATION OF W(Z) IN THE OTHER QUADRANTS
!
!
  If (yi<0.0) Then
!
    If (a) Then
      u2 = 2*u2
      v2 = 2*v2
    Else
      xquad = -xquad
!
!
!         THE FOLLOWING IF-STATEMENT PROTECTS 2*EXP(-Z**2)
!         AGAINST OVERFLOW
!
      If ((yquad>rmaxgoni) .Or. (xquad>rmaxexp)) Goto 100
!
      w1 = 2*dexp(xquad)
      u2 = w1*dcos(yquad)
      v2 = -w1*dsin(yquad)
    End If
!
    u = u2 - u
    v = v2 - v
    If (xi>0.0) v = -v
  Else
    If (xi<0.0) v = -v
  End If
!
  Return
!
  100 flag = .True.
  Return
!
End Subroutine wofz
end module integral