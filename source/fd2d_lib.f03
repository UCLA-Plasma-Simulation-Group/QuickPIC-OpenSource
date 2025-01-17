module fd2d_lib
   
implicit none

public

contains

subroutine divf_2order( f, df, m, n, dx, kstrt, nvp )

  implicit none

  real, intent(inout), dimension(:,:,:), pointer :: f, df
  integer, intent(in) :: m, n, kstrt, nvp
  real, intent(in) :: dx

  ! local variables
  integer :: i, j, k, kl, ku
  real :: inv_dxh, a1, a2, a3

  inv_dxh = 0.5 / dx
  a1 = -3.0 * inv_dxh
  a2 =  4.0 * inv_dxh
  a3 = -1.0 * inv_dxh

  kl = 1; ku = n
  if ( kstrt == 1 ) kl = 2
  if ( kstrt == nvp ) ku = n-1

  ! interior
  do k = kl, ku
    do j = 2, m-1
      df(1,j,k) = inv_dxh * ( f(1,j+1,k) - f(1,j-1,k) + f(2,j,k+1) - f(2,j,k-1) )
    enddo
  enddo

  ! top and bottom boundaries
  do k = kl, ku
    df(1,1,k) =  a1 * f(1,1,k) + a2 * f(1,2,k)   + a3 * f(1,3,k)   + inv_dxh * ( f(2,1,k+1) - f(2,1,k-1) )
    df(1,m,k) = -a1 * f(1,m,k) - a2 * f(1,m-1,k) - a3 * f(1,m-2,k) + inv_dxh * ( f(2,m,k+1) - f(2,m,k-1) )
  enddo

  ! left boundary
  if ( kstrt == 1 ) then
    do j = 2, m-1
      df(1,j,1) = inv_dxh * ( f(1,j+1,1) - f(1,j-1,1) ) + a1 * f(2,j,1) + a2 * f(2,j,2) + a3 * f(2,j,3)
    enddo
    df(1,1,1) =  a1 * ( f(1,1,1) + f(2,1,1) ) + a2 * ( f(1,2,1)   + f(2,1,2) ) + a3 * ( f(1,3,1)   + f(2,1,3) )
    df(1,m,1) = -a1 * ( f(1,m,1) - f(2,m,1) ) - a2 * ( f(1,m-1,1) - f(2,m,2) ) - a3 * ( f(1,m-2,1) - f(2,m,3) )
  endif

  ! right boundary
  if ( kstrt == nvp ) then
    do j = 2, m-1
      df(1,j,n) = inv_dxh * ( f(1,j+1,n) - f(1,j-1,n) ) - a1 * f(2,j,n) - a2 * f(2,j,n-1) - a3 * f(2,j,n-2)
    enddo
    df(1,1,n) =  a1 * ( f(1,1,1) - f(2,1,n) ) + a2 * ( f(1,2,1)   - f(2,1,n-1) ) + a3 * ( f(1,3,1)   - f(2,1,n-2) )
    df(1,m,n) = -a1 * ( f(1,m,1) + f(2,m,n) ) - a2 * ( f(1,m-1,1) + f(2,m,n-1) ) - a3 * ( f(1,m-2,1) + f(2,m,n-2) )
  endif

end subroutine divf_2order

subroutine divf_4order( f, df, m, n, dx, kstrt, nvp )

  implicit none

  real, intent(inout), dimension(:,:,:), pointer :: f, df
  integer, intent(in) :: m, n, kstrt, nvp
  real, intent(in) :: dx

  ! local variables
  integer :: i, j, k, kl, ku
  real :: inv_dx, a1, a2, a3, a4, a5, b1, b2, b3, b4, b5, c1, c2

  inv_dx = 1.0 / dx
  c1 =  0.666666666666667 * inv_dx
  c2 = -0.083333333333333 * inv_dx
  a1 = -0.25 * inv_dx
  a2 = -0.833333333333333 * inv_dx
  a3 =  1.5 * inv_dx
  a4 = -0.5 * inv_dx
  a5 =  0.083333333333333 * inv_dx
  b1 = -2.083333333333333 * inv_dx
  b2 =  4.0 * inv_dx
  b3 = -3.0 * inv_dx
  b4 =  1.333333333333333 * inv_dx
  b5 = -0.25 * inv_dx

  kl = 1; ku = n
  if ( kstrt == 1 ) kl = 3
  if ( kstrt == nvp ) ku = n-2

  ! interior
  do k = kl, ku
    do j = 3, m-2
      df(1,j,k) = c1 * ( f(1,j+1,k) - f(1,j-1,k) + f(2,j,k+1) - f(2,j,k-1) ) + &
                  c2 * ( f(1,j+2,k) - f(1,j-2,k) + f(2,j,k+2) - f(2,j,k-2) )
    enddo
  enddo

  ! top and bottom boundaries
  do k = kl, ku
    df(1,1,k)   =  b1 * f(1,1,k) + b2 * f(1,2,k) + b3 * f(1,3,k) + b4 * f(1,4,k) + b5 * f(1,5,k) + &
                   c1 * ( f(2,1,k+1) - f(2,1,k-1) ) + c2 * ( f(2,1,k+2) - f(2,1,k-2) )
    df(1,2,k)   =  a1 * f(1,1,k) + a2 * f(1,2,k) + a3 * f(1,3,k) + a4 * f(1,4,k) + a5 * f(1,5,k) + &
                   c1 * ( f(2,2,k+1) - f(2,2,k-1) ) + c2 * ( f(2,2,k+2) - f(2,2,k-2) )
    df(1,m-1,k) = -a1 * f(1,m,k) - a2 * f(1,m-1,k) - a3 * f(1,m-2,k) - a4 * f(1,m-3,k) - a5 * f(1,m-4,k) + &
                   c1 * ( f(2,m-1,k+1) - f(2,m-1,k-1) ) + c2 * ( f(2,m-1,k+2) - f(2,m-1,k-2) )
    df(1,m,k)   = -b1 * f(1,m,k) - b2 * f(1,m-1,k) - b3 * f(1,m-2,k) - b4 * f(1,m-3,k) - b5 * f(1,m-4,k) + &
                   c1 * ( f(2,m,k+1) - f(2,m,k-1) ) + c2 * ( f(2,m,k+2) - f(2,m,k-2) )
  enddo

  ! left boundary
  if ( kstrt == 1 ) then
    do j = 3, m-2
      df(1,j,1) = c1 * ( f(1,j+1,1) - f(1,j-1,1) ) + c2 * ( f(1,j+2,1) - f(1,j-2,1) ) + &
                  b1 * f(2,j,1) + b2 * f(2,j,2) + b3 * f(2,j,3) + b4 * f(2,j,4) + b5 * f(2,j,5)
      df(1,j,2) = c1 * ( f(1,j+1,2) - f(1,j-1,2) ) + c2 * ( f(1,j+2,2) - f(1,j-2,2) ) + &
                  a1 * f(2,j,1) + a2 * f(2,j,2) + a3 * f(2,j,3) + a4 * f(2,j,4) + a5 * f(2,j,5)
    enddo
    df(1,1,1)   =  b1 * f(1,1,1)   + b2 * f(1,2,1)   + b3 * f(1,3,1)   + b4 * f(1,4,1)   + b5 * f(1,5,1) + &
                   b1 * f(2,1,1)   + b2 * f(2,1,2)   + b3 * f(2,1,3)   + b4 * f(2,1,4)   + b5 * f(2,1,5)
    df(1,2,1)   =  a1 * f(1,1,1)   + a2 * f(1,2,1)   + a3 * f(1,3,1)   + a4 * f(1,4,1)   + a5 * f(1,5,1) + &
                   b1 * f(2,2,1)   + b2 * f(2,2,2)   + b3 * f(2,2,3)   + b4 * f(2,2,4)   + b5 * f(2,2,5)
    df(1,m-1,1) = -a1 * f(1,m,1)   - a2 * f(1,m-1,1) - a3 * f(1,m-2,1) - a4 * f(1,m-3,1) - a5 * f(1,m-4,1) + &
                   b1 * f(2,m-1,1) + b2 * f(2,m-1,2) + b3 * f(2,m-1,3) + b4 * f(2,m-1,4) + b5 * f(2,m-1,5)
    df(1,m,1)   = -b1 * f(1,m,1)   - b2 * f(1,m-1,1) - b3 * f(1,m-2,1) - b4 * f(1,m-3,1) - b5 * f(1,m-4,1) + &
                   b1 * f(2,m,1)   + b2 * f(2,m,2)   + b3 * f(2,m,3)   + b4 * f(2,m,4)   + b5 * f(2,m,5)

    df(1,1,2)   =  b1 * f(1,1,2)   + b2 * f(1,2,2)   + b3 * f(1,3,2)   + b4 * f(1,4,2)   + b5 * f(1,5,2) + &
                   a1 * f(2,1,1)   + a2 * f(2,1,2)   + a3 * f(2,1,3)   + a4 * f(2,1,4)   + a5 * f(2,1,5)
    df(1,2,2)   =  a1 * f(1,1,2)   + a2 * f(1,2,2)   + a3 * f(1,3,2)   + a4 * f(1,4,2)   + a5 * f(1,5,2) + &
                   a1 * f(2,2,1)   + a2 * f(2,2,2)   + a3 * f(2,2,3)   + a4 * f(2,2,4)   + a5 * f(2,2,5)
    df(1,m-1,2) = -a1 * f(1,m,2)   - a2 * f(1,m-1,2) - a3 * f(1,m-2,2) - a4 * f(1,m-3,2) - a5 * f(1,m-4,2) + &
                   a1 * f(2,m-1,1) + a2 * f(2,m-1,2) + a3 * f(2,m-1,3) + a4 * f(2,m-1,4) + a5 * f(2,m-1,5)
    df(1,m,2)   = -b1 * f(1,m,2)   - b2 * f(1,m-1,2) - b3 * f(1,m-2,2) - b4 * f(1,m-3,2) - b5 * f(1,m-4,2) + &
                   a1 * f(2,m,1)   + a2 * f(2,m,2)   + a3 * f(2,m,3)   + a4 * f(2,m,4)   + a5 * f(2,m,5)
  endif

  ! right boundary
  if ( kstrt == nvp ) then
    do j = 3, m-2
      df(1,j,n-1) = c1 * ( f(1,j+1,n-1) - f(1,j-1,n-1) ) + c2 * ( f(1,j+2,n-1) - f(1,j-2,n-1) ) - &
                    a1 * f(2,j,n) - a2 * f(2,j,n-1) - a3 * f(2,j,n-2) - a4 * f(2,j,n-3) - a5 * f(2,j,n-4)
      df(1,j,n)   = c1 * ( f(1,j+1,n) - f(1,j-1,n) ) + c2 * ( f(1,j+2,n) - f(1,j-2,n) ) - &
                    b1 * f(2,j,n) - b2 * f(2,j,n-1) - b3 * f(2,j,n-2) - b4 * f(2,j,n-3) - b5 * f(2,j,n-4)
    enddo
    df(1,1,n-1)   =  b1 * f(1,1,n-1) + b2 * f(1,2,n-1)   + b3 * f(1,3,n-1)   + b4 * f(1,4,n-1)   + b5 * f(1,5,n-1) - &
                     a1 * f(2,1,n)   - a2 * f(2,1,n-1)   - a3 * f(2,1,n-2)   - a4 * f(2,1,n-3)   - a5 * f(2,1,n-4)
    df(1,2,n-1)   =  a1 * f(1,1,n-1) + a2 * f(1,2,n-1)   + a3 * f(1,3,n-1)   + a4 * f(1,4,n-1)   + a5 * f(1,5,n-1) - &
                     a1 * f(2,2,n)   - a2 * f(2,2,n-1)   - a3 * f(2,2,n-2)   - a4 * f(2,2,n-3)   - a5 * f(2,2,n-4)
    df(1,m-1,n-1) = -a1 * f(1,m,n-1) - a2 * f(1,m-1,n-1) - a3 * f(1,m-2,n-1) - a4 * f(1,m-3,n-1) - a5 * f(1,m-4,n-1) - &
                     a1 * f(2,m-1,n) - a2 * f(2,m-1,n-1) - a3 * f(2,m-1,n-2) - a4 * f(2,m-1,n-3) - a5 * f(2,m-1,n-4)
    df(1,m,n-1)   = -b1 * f(1,m,n-1) - b2 * f(1,m-1,n-1) - b3 * f(1,m-2,n-1) - b4 * f(1,m-3,n-1) - b5 * f(1,m-4,n-1) - &
                     a1 * f(2,m,n)   - a2 * f(2,m,n-1)   - a3 * f(2,m,n-2)   - a4 * f(2,m,n-3)   - a5 * f(2,m,n-4)

    df(1,1,n)   =  b1 * f(1,1,n)   + b2 * f(1,2,n)     + b3 * f(1,3,n)     + b4 * f(1,4,n)     + b5 * f(1,5,n) - &
                   b1 * f(2,1,n)   - b2 * f(2,1,n-1)   - b3 * f(2,1,n-2)   - b4 * f(2,1,n-3)   - b5 * f(2,1,n-4)
    df(1,2,n)   =  a1 * f(1,1,n)   + a2 * f(1,2,n)     + a3 * f(1,3,n)     + a4 * f(1,4,n)     + a5 * f(1,5,n) - &
                   b1 * f(2,2,n)   - b2 * f(2,2,n-1)   - b3 * f(2,2,n-2)   - b4 * f(2,2,n-3)   - b5 * f(2,2,n-4)
    df(1,m-1,n) = -a1 * f(1,m,n)   - a2 * f(1,m-1,n)   - a3 * f(1,m-2,n)   - a4 * f(1,m-3,n)   - a5 * f(1,m-4,n) - &
                   b1 * f(2,m-1,n) - b2 * f(2,m-1,n-1) - b3 * f(2,m-1,n-2) - b4 * f(2,m-1,n-3) - b5 * f(2,m-1,n-4)
    df(1,m,n)   = -b1 * f(1,m,n)   - b2 * f(1,m-1,n)   - b3 * f(1,m-2,n)   - b4 * f(1,m-3,n)   - b5 * f(1,m-4,n) - &
                   b1 * f(2,m,n)   - b2 * f(2,m,n-1)   - b3 * f(2,m,n-2)   - b4 * f(2,m,n-3)   - b5 * f(2,m,n-4) 
  endif

end subroutine divf_4order

subroutine divf_dyad_2order( f, df, m, n, dx, kstrt, nvp )

  implicit none

  real, intent(inout), dimension(:,:,:), pointer :: f, df
  integer, intent(in) :: m, n, kstrt, nvp
  real, intent(in) :: dx

  ! local variables
  integer :: i, j, k, kl, ku
  real :: inv_dxh, a1, a2, a3

  inv_dxh = 0.5 / dx
  a1 = -3.0 * inv_dxh
  a2 =  4.0 * inv_dxh
  a3 = -1.0 * inv_dxh

  kl = 1; ku = n
  if ( kstrt == 1 ) kl = 2
  if ( kstrt == nvp ) ku = n-1

  ! interior
  do k = kl, ku
    do j = 2, m-1
      df(1,j,k) = inv_dxh * ( f(1,j+1,k) - f(1,j-1,k) + f(3,j,k+1) - f(3,j,k-1) )
      df(2,j,k) = inv_dxh * ( f(3,j+1,k) - f(3,j-1,k) + f(2,j,k+1) - f(2,j,k-1) )
    enddo
  enddo

  ! top and bottom boundaries
  do k = kl, ku
    df(1,1,k) =  a1 * f(1,1,k) + a2 * f(1,2,k)   + a3 * f(1,3,k)   + inv_dxh * ( f(3,1,k+1) - f(3,1,k-1) )
    df(2,1,k) =  a1 * f(3,1,k) + a2 * f(3,2,k)   + a3 * f(3,3,k)   + inv_dxh * ( f(2,1,k+1) - f(2,1,k-1) )
    df(1,m,k) = -a1 * f(1,m,k) - a2 * f(1,m-1,k) - a3 * f(1,m-2,k) + inv_dxh * ( f(3,m,k+1) - f(3,m,k-1) )
    df(2,m,k) = -a1 * f(3,m,k) - a2 * f(3,m-1,k) - a3 * f(3,m-2,k) + inv_dxh * ( f(2,m,k+1) - f(2,m,k-1) )
  enddo

  ! left boundary
  if ( kstrt == 1 ) then
    do j = 2, m-1
      df(1,j,1) = inv_dxh * ( f(1,j+1,1) - f(1,j-1,1) ) + a1 * f(3,j,1) + a2 * f(3,j,2) + a3 * f(3,j,3)
      df(2,j,1) = inv_dxh * ( f(3,j+1,1) - f(3,j-1,1) ) + a1 * f(2,j,1) + a2 * f(2,j,2) + a3 * f(2,j,3)
    enddo
    df(1,1,1) =  a1 * f(1,1,1) + a2 * f(1,2,1)   + a3 * f(1,3,1) + &
                 a1 * f(3,1,1) + a2 * f(3,1,2)   + a3 * f(3,1,3)
    df(2,1,1) =  a1 * f(3,1,1) + a2 * f(3,2,1)   + a3 * f(3,3,1) + &
                 a1 * f(2,1,1) + a2 * f(2,1,2)   + a3 * f(2,1,3)
    df(1,m,1) = -a1 * f(1,m,1) - a2 * f(1,m-1,1) - a3 * f(1,m-2,1) + &
                 a1 * f(3,m,1) + a2 * f(3,m,2)   + a3 * f(3,m,3)
    df(2,m,1) = -a1 * f(3,m,1) - a2 * f(3,m-1,1) - a3 * f(3,m-2,1) + &
                 a1 * f(2,m,1) + a2 * f(2,m,2)   + a3 * f(2,m,3)
  endif

  ! right boundary
  if ( kstrt == nvp ) then
    do j = 2, m-1
      df(1,j,m) = inv_dxh * ( f(1,j+1,m) - f(1,j-1,m) ) - a1 * f(3,j,m) - a2 * f(3,j,m-1) - a3 * f(3,j,m-2)
      df(2,j,m) = inv_dxh * ( f(3,j+1,m) - f(3,j-1,m) ) - a1 * f(2,j,m) - a2 * f(2,j,m-1) - a3 * f(2,j,m-2)
    enddo
    df(1,1,m) =  a1 * f(1,1,m) + a2 * f(1,2,m)   + a3 * f(1,3,m) - &
                 a1 * f(3,1,m) - a2 * f(3,1,m-1) - a3 * f(3,1,m-2)
    df(2,1,m) =  a1 * f(3,1,m) + a2 * f(3,2,m)   + a3 * f(3,3,m) - &
                 a1 * f(2,1,m) - a2 * f(2,1,m-1) - a3 * f(2,1,m-2)
    df(1,m,m) = -a1 * f(1,m,m) - a2 * f(1,m-1,m) - a3 * f(1,m-2,m) - &
                 a1 * f(3,m,m) - a2 * f(3,m,m-1) - a3 * f(3,m,m-2)
    df(2,m,m) = -a1 * f(3,m,m) - a2 * f(3,m-1,m) - a3 * f(3,m-2,m) - &
                 a1 * f(2,m,m) - a2 * f(2,m,m-1) - a3 * f(2,m,m-2)
  endif

end subroutine divf_dyad_2order

subroutine gradf_2order( f, df, m, n, dx, kstrt, nvp )

  implicit none

  real, intent(inout), dimension(:,:,:), pointer :: f, df
  integer, intent(in) :: m, n, kstrt, nvp
  real, intent(in) :: dx

  ! local variables
  integer :: i, j, k, kl, ku
  real :: inv_dxh, a1, a2, a3

  inv_dxh = 0.5 / dx
  a1 = -3.0 * inv_dxh
  a2 =  4.0 * inv_dxh
  a3 = -1.0 * inv_dxh

  kl = 1; ku = n
  if ( kstrt == 1 ) kl = 2
  if ( kstrt == nvp ) ku = n-1

  ! interior
  do k = kl, ku
    do j = 2, m-1
      df(1,j,k) = inv_dxh * ( f(1,j+1,k) - f(1,j-1,k) )
      df(2,j,k) = inv_dxh * ( f(1,j,k+1) - f(1,j,k-1) )
    enddo
  enddo

  ! top and bottom boundaries
  do k = kl, ku
    df(1,1,k) =  a1 * f(1,1,k) + a2 * f(1,2,k) + a3 * f(1,3,k)
    df(2,1,k) =  inv_dxh * ( f(1,1,k+1) - f(1,1,k-1) )
    df(1,m,k) = -a1 * f(1,m,k) - a2 * f(1,m-1,k) - a3 * f(1,m-2,k)
    df(2,m,k) =  inv_dxh * ( f(1,m,k+1) - f(1,m,k-1) )
  enddo

  ! left boundary
  if ( kstrt == 1 ) then
    do j = 2, m-1
      df(1,j,1) = inv_dxh * ( f(1,j+1,1) - f(1,j-1,1) )
      df(2,j,1) = a1 * f(1,j,1) + a2 * f(1,j,2) + a3 * f(1,j,3)
    enddo
    df(1,1,1) =  a1 * f(1,1,1) + a2 * f(1,2,1)   + a3 * f(1,3,1)
    df(2,1,1) =  a1 * f(1,1,1) + a2 * f(1,1,2)   + a3 * f(1,1,3)
    df(1,m,1) = -a1 * f(1,m,1) - a2 * f(1,m-1,1) - a3 * f(1,m-2,1)
    df(2,m,1) =  a1 * f(1,m,1) + a2 * f(1,m,2)   + a3 * f(1,m,3)
  endif

  ! right boundary
  if ( kstrt == nvp ) then
    do j = 2, m-1
      df(1,j,n) =  inv_dxh * ( f(1,j+1,n) - f(1,j-1,n) )
      df(2,j,n) = -a1 * f(1,j,n) - a2 * f(1,j,n-1) - a3 * f(1,j,n-2)
    enddo
    df(1,1,n) =  a1 * f(1,1,n) + a2 * f(1,2,n)   + a3 * f(1,3,n)
    df(2,1,n) = -a1 * f(1,1,n) - a2 * f(1,1,n-1) - a3 * f(1,1,n-2)
    df(1,m,n) = -a1 * f(1,m,n) - a2 * f(1,m-1,n) - a3 * f(1,m-2,n)
    df(2,m,n) = -a1 * f(1,m,n) - a2 * f(1,m,n-1) - a3 * f(1,m,n-2)
  endif

end subroutine gradf_2order

subroutine gradf_4order( f, df, m, n, dx, kstrt, nvp )

  implicit none

  real, intent(inout), dimension(:,:,:) :: f, df
  integer, intent(in) :: m, n, kstrt, nvp
  real, intent(in) :: dx

  ! local variables
  integer :: i, j, k, kl, ku
  real :: inv_dx, c1, c2, a1, a2, a3, a4, a5, b1, b2, b3, b4, b5

  inv_dx = 1.0 / dx
  c1 =  0.666666666666667 * inv_dx
  c2 = -0.083333333333333 * inv_dx
  a1 = -0.25 * inv_dx
  a2 = -0.833333333333333 * inv_dx
  a3 =  1.5 * inv_dx
  a4 = -0.5 * inv_dx
  a5 =  0.083333333333333 * inv_dx
  b1 = -2.083333333333333 * inv_dx
  b2 =  4.0 * inv_dx
  b3 = -3.0 * inv_dx
  b4 =  1.333333333333333 * inv_dx
  b5 = -0.25 * inv_dx

  kl = 1; ku = n
  if ( kstrt == 1 ) kl = 3
  if ( kstrt == nvp ) ku = n-2

  ! interior
  do k = kl, ku
    do j = 3, m-2
      df(1,j,k) = c1 * ( f(1,j+1,k) - f(1,j-1,k) ) + c2 * ( f(1,j+2,k) - f(1,j-2,k) ) 
      df(2,j,k) = c1 * ( f(1,j,k+1) - f(1,j,k-1) ) + c2 * ( f(1,j,k+2) - f(1,j,k-2) )
    enddo
  enddo

  ! top and bottom boundaries
  do k = kl, ku
    df(1,1,k)   =  b1 * f(1,1,k) + b2 * f(1,2,k) + b3 * f(1,3,k) + b4 * f(1,4,k) + b5 * f(1,5,k)
    df(2,1,k)   =  c1 * ( f(1,1,k+1) - f(1,1,k-1) ) + c2 * ( f(1,1,k+2) - f(1,1,k-2) )
    df(1,2,k)   =  a1 * f(1,1,k) + a2 * f(1,2,k) + a3 * f(1,3,k) + a4 * f(1,4,k) + a5 * f(1,5,k)
    df(2,2,k)   =  c1 * ( f(1,2,k+1) - f(1,2,k-1) ) + c2 * ( f(1,2,k+2) - f(1,2,k-2) )
    df(1,m-1,k) = -a1 * f(1,m,k) - a2 * f(1,m-1,k) - a3 * f(1,m-2,k) - a4 * f(1,m-3,k) - a5 * f(1,m-4,k)
    df(2,m-1,k) =  c1 * ( f(1,m-1,k+1) - f(1,m-1,k-1) ) + c2 * ( f(1,m-1,k+2) - f(1,m-1,k-2) )          
    df(1,m,k)   = -b1 * f(1,m,k) - b2 * f(1,m-1,k) - b3 * f(1,m-2,k) - b4 * f(1,m-3,k) - b5 * f(1,m-4,k)
    df(2,m,k)   =  c1 * ( f(1,m,k+1) - f(1,m,k-1) ) + c2 * ( f(1,m,k+2) - f(1,m,k-2) ) 
  enddo

  ! left boundary
  if ( kstrt == 1 ) then
    do j = 3, m-2
      df(1,j,1) = c1 * ( f(1,j+1,1) - f(1,j-1,1) ) + c2 * ( f(1,j+2,1) - f(1,j-2,1) )
      df(2,j,1) = b1 * f(1,j,1) + b2 * f(1,j,2) + b3 * f(1,j,3) + b4 * f(1,j,4) + b5 * f(1,j,5)
      df(1,j,2) = c1 * ( f(1,j+1,2) - f(1,j-1,2) ) + c2 * ( f(1,j+2,2) - f(1,j-2,2) )
      df(2,j,2) = a1 * f(1,j,1) + a2 * f(1,j,2) + a3 * f(1,j,3) + a4 * f(1,j,4) + a5 * f(1,j,5)
    enddo
    df(1,1,1)   =  b1 * f(1,1,1)   + b2 * f(1,2,1)   + b3 * f(1,3,1)   + b4 * f(1,4,1)   + b5 * f(1,5,1)
    df(2,1,1)   =  b1 * f(1,1,1)   + b2 * f(1,1,2)   + b3 * f(1,1,3)   + b4 * f(1,1,4)   + b5 * f(1,1,5)
    df(1,2,1)   =  a1 * f(1,1,1)   + a2 * f(1,2,1)   + a3 * f(1,3,1)   + a4 * f(1,4,1)   + a5 * f(1,5,1)
    df(2,2,1)   =  b1 * f(1,2,1)   + b2 * f(1,2,2)   + b3 * f(1,2,3)   + b4 * f(1,2,4)   + b5 * f(1,2,5)
    df(1,m-1,1) = -a1 * f(1,m,1)   - a2 * f(1,m-1,1) - a3 * f(1,m-2,1) - a4 * f(1,m-3,1) - a5 * f(1,m-4,1)
    df(2,m-1,1) =  b1 * f(1,m-1,1) + b2 * f(1,m-1,2) + b3 * f(1,m-1,3) + b4 * f(1,m-1,4) + b5 * f(1,m-1,5)
    df(1,m,1)   = -b1 * f(1,m,1)   - b2 * f(1,m-1,1) - b3 * f(1,m-2,1) - b4 * f(1,m-3,1) - b5 * f(1,m-4,1)
    df(2,m,1)   =  b1 * f(1,m,1)   + b2 * f(1,m,2)   + b3 * f(1,m,3)   + b4 * f(1,m,4)   + b5 * f(1,m,5)
    df(1,1,2)   =  b1 * f(1,1,2)   + b2 * f(1,2,2)   + b3 * f(1,3,2)   + b4 * f(1,4,2)   + b5 * f(1,5,2)
    df(2,1,2)   =  a1 * f(1,1,1)   + a2 * f(1,1,2)   + a3 * f(1,1,3)   + a4 * f(1,1,4)   + a5 * f(1,1,5)
    df(1,2,2)   =  a1 * f(1,1,2)   + a2 * f(1,2,2)   + a3 * f(1,3,2)   + a4 * f(1,4,2)   + a5 * f(1,5,2)
    df(2,2,2)   =  a1 * f(1,2,1)   + a2 * f(1,2,2)   + a3 * f(1,2,3)   + a4 * f(1,2,4)   + a5 * f(1,2,5)
    df(1,m-1,2) = -a1 * f(1,m,2)   - a2 * f(1,m-1,2) - a3 * f(1,m-2,2) - a4 * f(1,m-3,2) - a5 * f(1,m-4,2)
    df(2,m-1,2) =  a1 * f(1,m-1,1) + a2 * f(1,m-1,2) + a3 * f(1,m-1,3) + a4 * f(1,m-1,4) + a5 * f(1,m-1,5)
    df(1,m,2)   = -b1 * f(1,m,2)   - b2 * f(1,m-1,2) - b3 * f(1,m-2,2) - b4 * f(1,m-3,2) - b5 * f(1,m-4,2)
    df(2,m,2)   =  a1 * f(1,m,1)   + a2 * f(1,m,2)   + a3 * f(1,m,3)   + a4 * f(1,m,4)   + a5 * f(1,m,5)

  endif

  ! right boundary
  if ( kstrt == nvp ) then
    do j = 3, m-2
      df(1,j,n-1) =  c1 * ( f(1,j+1,n-1) - f(1,j-1,n-1) ) + c2 * ( f(1,j+2,n-1) - f(1,j-2,n-1) )
      df(2,j,n-1) = -a1 * f(1,j,n) - a2 * f(1,j,n-1) - a3 * f(1,j,n-2) - a4 * f(1,j,n-3) - a5 * f(1,j,n-4)
      df(1,j,n)   =  c1 * ( f(1,j+1,n) - f(1,j-1,n) ) + c2 * ( f(1,j+2,n) - f(1,j-2,n) )
      df(2,j,n)   = -b1 * f(1,j,n) - b2 * f(1,j,n-1) - b3 * f(1,j,n-2) - b4 * f(1,j,n-3) - b5 * f(1,j,n-4)
    enddo
    df(1,1,n-1)   =  b1 * f(1,1,n-1) + b2 * f(1,2,n-1)   + b3 * f(1,3,n-1)   + b4 * f(1,4,n-1)   + b5 * f(1,5,n-1)
    df(2,1,n-1)   = -a1 * f(1,1,n)   - a2 * f(1,1,n-1)   - a3 * f(1,1,n-2)   - a4 * f(1,1,n-3)   - a5 * f(1,1,n-4)
    df(1,2,n-1)   =  a1 * f(1,1,n-1) + a2 * f(1,2,n-1)   + a3 * f(1,3,n-1)   + a4 * f(1,4,n-1)   + a5 * f(1,5,n-1)
    df(2,2,n-1)   = -a1 * f(1,2,n)   - a2 * f(1,2,n-1)   - a3 * f(1,2,n-2)   - a4 * f(1,2,n-3)   - a5 * f(1,2,n-4)
    df(1,m-1,n-1) = -a1 * f(1,m,n-1) - a2 * f(1,m-1,n-1) - a3 * f(1,m-2,n-1) - a4 * f(1,m-3,n-1) - a5 * f(1,m-4,n-1)
    df(2,m-1,n-1) = -a1 * f(1,m-1,n) - a2 * f(1,m-1,n-1) - a3 * f(1,m-1,n-2) - a4 * f(1,m-1,n-3) - a5 * f(1,m-1,n-4)
    df(1,m,n-1)   = -b1 * f(1,m,n-1) - b2 * f(1,m-1,n-1) - b3 * f(1,m-2,n-1) - b4 * f(1,m-3,n-1) - b5 * f(1,m-4,n-1)
    df(2,m,n-1)   = -a1 * f(1,m,n)   - a2 * f(1,m,n-1)   - a3 * f(1,m,n-2)   - a4 * f(1,m,n-3)   - a5 * f(1,m,n-4)
    df(1,1,n)     =  b1 * f(1,1,n)   + b2 * f(1,2,n)     + b3 * f(1,3,n)     + b4 * f(1,4,n)     + b5 * f(1,5,n)
    df(2,1,n)     = -b1 * f(1,1,n)   - b2 * f(1,1,n-1)   - b3 * f(1,1,n-2)   - b4 * f(1,1,n-3)   - b5 * f(1,1,n-4)
    df(1,2,n)     =  a1 * f(1,1,n)   + a2 * f(1,2,n)     + a3 * f(1,3,n)     + a4 * f(1,4,n)     + a5 * f(1,5,n)
    df(2,2,n)     = -b1 * f(1,2,n)   - b2 * f(1,2,n-1)   - b3 * f(1,2,n-2)   + b4 * f(1,2,n-3)   - b5 * f(1,2,n-4)
    df(1,m-1,n)   = -a1 * f(1,m,n)   - a2 * f(1,m-1,n)   - a3 * f(1,m-2,n)   - a4 * f(1,m-3,n)   - a5 * f(1,m-4,n)
    df(2,m-1,n)   = -b1 * f(1,m-1,n) - b2 * f(1,m-1,n-1) - b3 * f(1,m-1,n-2) - b4 * f(1,m-1,n-3) - b5 * f(1,m-1,n-4)
    df(1,m,n)     = -b1 * f(1,m,n)   - b2 * f(1,m-1,n)   - b3 * f(1,m-2,n)   - b4 * f(1,m-3,n)   - b5 * f(1,m-4,n)
    df(2,m,n)     = -b1 * f(1,m,n)   - b2 * f(1,m,n-1)   - b3 * f(1,m,n-2)   - b4 * f(1,m,n-3)   - b5 * f(1,m,n-4)
  endif

end subroutine gradf_4order

subroutine curlf_2order( f, df, m, n, dx, kstrt, nvp )

  implicit none

  real, intent(inout), dimension(:,:,:), pointer :: f, df
  integer, intent(in) :: m, n, kstrt, nvp
  real, intent(in) :: dx

  ! local variables
  integer :: i, j, k, kl, ku
  real :: inv_dxh, a1, a2, a3

  inv_dxh = 0.5 / dx
  a1 = -3.0 * inv_dxh
  a2 =  4.0 * inv_dxh
  a3 = -1.0 * inv_dxh

  kl = 1; ku = n
  if ( kstrt == 1 ) kl = 2
  if ( kstrt == nvp ) ku = n-1

  ! interior
  do k = kl, ku
    do j = 2, m-1
      df(1,j,k) =  inv_dxh * ( f(3,j,k+1) - f(3,j,k-1) )
      df(2,j,k) = -inv_dxh * ( f(3,j+1,k) - f(3,j-1,k) )
      df(3,j,k) =  inv_dxh * ( f(2,j+1,k) - f(2,j-1,k) - f(1,j,k+1) + f(1,j,k-1) )
    enddo
  enddo

  ! top and bottom boundaries
  do k = kl, ku
    df(1,1,k) =  inv_dxh * ( f(3,1,k+1) - f(3,1,k-1) )
    df(2,1,k) = -a1 * f(3,1,k) - a2 * f(3,2,k) - a3 * f(3,3,k)
    df(3,1,k) =  a1 * f(2,1,k) + a2 * f(2,2,k) + a3 * f(2,3,k) - inv_dxh * ( f(1,1,k+1) - f(1,1,k-1) )
    df(1,m,k) =  inv_dxh * ( f(3,m,k+1) - f(3,m,k-1) )
    df(2,m,k) =  a1 * f(3,m,k) + a2 * f(3,m-1,k) + a3 * f(3,m-2,k)
    df(3,m,k) = -a1 * f(2,m,k) - a2 * f(2,m-1,k) - a3 * f(2,m-2,k) - inv_dxh * ( f(1,m,k+1) - f(1,m,k-1) )
  enddo

  ! left boundary
  if ( kstrt == 1 ) then
    do j = 2, m-1
      df(1,j,1) =  a1 * f(3,j,1) + a2 * f(3,j,2) + a3 * f(3,j,3)
      df(2,j,1) = -inv_dxh * ( f(3,j+1,1) - f(3,j-1,1) )
      df(3,j,1) =  inv_dxh * ( f(2,j+1,1) - f(2,j-1,1) ) - a1 * f(1,j,1) - a2 * f(1,j,2) - a3 * f(1,j,3)
    enddo
    df(1,1,1) =  a1 * f(3,1,1) + a2 * f(3,1,2)   + a3 * f(3,1,3)
    df(2,1,1) = -a1 * f(3,1,1) - a2 * f(3,2,1)   - a3 * f(3,3,1)
    df(3,1,1) =  a1 * f(2,1,1) + a2 * f(2,2,1)   + a3 * f(2,3,1)   - a1 * f(1,1,1) - a2 * f(1,1,2) - a3 * f(1,1,3)
    df(1,m,1) =  a1 * f(3,m,1) + a2 * f(3,m,2)   + a3 * f(3,m,3)
    df(2,m,1) =  a1 * f(3,m,1) + a2 * f(3,m-1,1) + a3 * f(3,m-2,1)
    df(3,m,1) = -a1 * f(2,m,1) - a2 * f(2,m-1,1) - a3 * f(2,m-2,1) - a1 * f(1,m,1) - a2 * f(1,m,2) - a3 * f(1,m,3)
  endif

  ! right boundary
  if ( kstrt == nvp ) then
    do j = 2, m-1
      df(1,j,n) = -a1 * f(3,j,n) - a2 * f(3,j,n-1) - a3 * f(3,j,n-2)
      df(2,j,n) = -inv_dxh * ( f(3,j+1,n) - f(3,j-1,n) )
      df(3,j,n) =  inv_dxh * ( f(2,j+1,n) - f(2,j-1,n) ) + a1 * f(1,j,n) + a2 * f(1,j,n-1) + a3 * f(1,j,n-2)
    enddo
    df(1,1,n) = -a1 * f(3,1,n) - a2 * f(3,1,n-1) - a3 * f(3,1,n-2)
    df(2,1,n) = -a1 * f(3,1,n) - a2 * f(3,2,n)   - a3 * f(3,3,n)
    df(3,1,n) =  a1 * f(2,1,n) + a2 * f(2,2,n)   + a3 * f(2,3,n)   + a1 * f(1,1,n) + a2 * f(1,1,n-1) + a3 * f(1,1,n-2)
    df(1,m,n) = -a1 * f(3,m,n) - a2 * f(3,m,n-1) - a3 * f(3,m,n-2)
    df(2,m,n) =  a1 * f(3,m,n) + a2 * f(3,m-1,n) + a3 * f(3,m-2,n)
    df(3,m,n) = -a1 * f(2,m,n) - a2 * f(2,m-1,n) - a3 * f(2,m-2,n) + a1 * f(1,m,n) + a2 * f(1,m,n-1) + a3 * f(1,m,n-2)
  endif

end subroutine curlf_2order

subroutine curlfz_2order( f, df, m, n, dx, kstrt, nvp )

  implicit none

  real, intent(inout), dimension(:,:,:), pointer :: f, df
  integer, intent(in) :: m, n, kstrt, nvp
  real, intent(in) :: dx

  ! local variables
  integer :: i, j, k, kl, ku
  real :: inv_dxh, a1, a2, a3

  inv_dxh = 0.5 / dx
  a1 = -3.0 * inv_dxh
  a2 =  4.0 * inv_dxh
  a3 = -1.0 * inv_dxh

  kl = 1; ku = n
  if ( kstrt == 1 ) kl = 2
  if ( kstrt == nvp ) ku = n-1

  ! interior
  do k = kl, ku
    do j = 2, m-1
      df(1,j,k) =  inv_dxh * ( f(2,j+1,k) - f(2,j-1,k) - f(1,j,k+1) + f(1,j,k-1) )
    enddo
  enddo

  ! top and bottom boundaries
  do k = kl, ku
    df(1,1,k) =  a1 * f(2,1,k) + a2 * f(2,2,k) + a3 * f(2,3,k) - inv_dxh * ( f(1,1,k+1) - f(1,1,k-1) )
    df(1,m,k) = -a1 * f(2,m,k) - a2 * f(2,m-1,k) - a3 * f(2,m-2,k) - inv_dxh * ( f(1,m,k+1) - f(1,m,k-1) )
  enddo

  ! left boundary
  if ( kstrt == 1 ) then
    do j = 2, m-1
      df(1,j,1) =  inv_dxh * ( f(2,j+1,1) - f(2,j-1,1) ) - a1 * f(1,j,1) - a2 * f(1,j,2) - a3 * f(1,j,3)
    enddo
    df(1,1,1) =  a1 * f(2,1,1) + a2 * f(2,2,1)   + a3 * f(2,3,1)   - a1 * f(1,1,1) - a2 * f(1,1,2) - a3 * f(1,1,3)
    df(1,m,1) = -a1 * f(2,m,1) - a2 * f(2,m-1,1) - a3 * f(2,m-2,1) - a1 * f(1,m,1) - a2 * f(1,m,2) - a3 * f(1,m,3)
  endif

  ! right boundary
  if ( kstrt == nvp ) then
    do j = 2, m-1
      df(1,j,n) =  inv_dxh * ( f(2,j+1,n) - f(2,j-1,n) ) + a1 * f(1,j,n) + a2 * f(1,j,n-1) + a3 * f(1,j,n-2)
    enddo
    df(1,1,n) =  a1 * f(2,1,n) + a2 * f(2,2,n)   + a3 * f(2,3,n)   + a1 * f(1,1,n) + a2 * f(1,1,n-1) + a3 * f(1,1,n-2)
    df(1,m,n) = -a1 * f(2,m,n) - a2 * f(2,m-1,n) - a3 * f(2,m-2,n) + a1 * f(1,m,n) + a2 * f(1,m,n-1) + a3 * f(1,m,n-2)
  endif

end subroutine curlfz_2order

subroutine curlf_4order( f, df, m, n, dx, kstrt, nvp )

  implicit none

  real, intent(inout), dimension(:,:,:), pointer :: f, df
  integer, intent(in) :: m, n, kstrt, nvp
  real, intent(in) :: dx

  ! local variables
  integer :: i, j, k, kl, ku
  real :: inv_dx, c1, c2, a1, a2, a3, a4, a5, b1, b2, b3, b4, b5

  inv_dx = 1.0 / dx
  c1 =  0.666666666666667 * inv_dx
  c2 = -0.083333333333333 * inv_dx
  a1 = -0.25 * inv_dx
  a2 = -0.833333333333333 * inv_dx
  a3 =  1.5 * inv_dx
  a4 = -0.5 * inv_dx
  a5 =  0.083333333333333 * inv_dx
  b1 = -2.083333333333333 * inv_dx
  b2 =  4.0 * inv_dx
  b3 = -3.0 * inv_dx
  b4 =  1.333333333333333 * inv_dx
  b5 = -0.25 * inv_dx

  kl = 1; ku = n
  if ( kstrt == 1 ) kl = 3
  if ( kstrt == nvp ) ku = n-2

  ! interior
  do k = kl, ku
    do j = 3, m-2
      df(1,j,k) =  c1 * ( f(3,j,k+1) - f(3,j,k-1) ) + c2 * ( f(3,j,k+2) - f(3,j,k-2) )
      df(2,j,k) = -c1 * ( f(3,j+1,k) - f(3,j-1,k) ) - c2 * ( f(3,j+2,k) - f(3,j-2,k) )
      df(3,j,k) =  c1 * ( f(2,j+1,k) - f(2,j-1,k) ) + c2 * ( f(2,j+2,k) - f(2,j-2,k) ) - &
                   c1 * ( f(1,j,k+1) - f(1,j,k-1) ) - c2 * ( f(1,j,k+2) - f(1,j,k-2) )
    enddo
  enddo

  ! top and bottom boundaries
  do k = kl, ku
    df(1,1,k)   =  c1 * ( f(3,1,k+1) - f(3,1,k-1) ) + c2 * ( f(3,1,k+2) - f(3,1,k-2) )
    df(2,1,k)   = -b1 * f(3,1,k) - b2 * f(3,2,k) - b3 * f(3,3,k) - b4 * f(3,4,k) - b5 * f(3,5,k)
    df(3,1,k)   =  b1 * f(2,1,k) + b2 * f(2,2,k) + b3 * f(2,3,k) + b4 * f(2,4,k) + b5 * f(2,5,k) - &
                   c1 * ( f(1,1,k+1) - f(1,1,k-1) ) - c2 * ( f(1,1,k+2) - f(1,1,k-2) )
    df(1,2,k)   =  c1 * ( f(3,2,k+1) - f(3,2,k-1) ) + c2 * ( f(3,2,k+2) - f(3,2,k-2) )
    df(2,2,k)   = -a1 * f(3,1,k) - a2 * f(3,2,k) - a3 * f(3,3,k) - a4 * f(3,4,k) - a5 * f(3,5,k)
    df(3,2,k)   =  a1 * f(2,1,k) + a2 * f(2,2,k) + a3 * f(2,3,k) + a4 * f(2,4,k) + a5 * f(2,5,k) - &
                   c1 * ( f(1,2,k+1) - f(1,2,k-1) ) - c2 * ( f(1,2,k+2) - f(1,2,k-2) )
    df(1,m-1,k) =  c1 * ( f(3-1,m,k+1) - f(3,m-1,k-1) ) + c2 * ( f(3,m-1,k+2) - f(3,m-1,k-2) )
    df(2,m-1,k) =  a1 * f(3,m,k) + a2 * f(3,m-1,k) + a3 * f(3,m-2,k) + a4 * f(3,m-3,k) + a5 * f(3,m-4,k)
    df(3,m-1,k) = -a1 * f(2,m,k) - a2 * f(2,m-1,k) - a3 * f(2,m-2,k) - a4 * f(2,m-3,k) - a5 * f(2,m-4,k) - &
                   c1 * ( f(1,m-1,k+1) - f(1,m-1,k-1) ) - c2 * ( f(1,m-1,k+2) - f(1,m-1,k-2) )
    df(1,m,k)   =  c1 * ( f(3,m,k+1) - f(3,m,k-1) ) + c2 * ( f(3,m,k+2) - f(3,m,k-2) )
    df(2,m,k)   =  b1 * f(3,m,k) + b2 * f(3,m-1,k) + b3 * f(3,m-2,k) + b4 * f(3,m-3,k) + b5 * f(3,m-4,k)
    df(3,m,k)   = -b1 * f(2,m,k) - b2 * f(2,m-1,k) - b3 * f(2,m-2,k) - b4 * f(2,m-3,k) - b5 * f(2,m-4,k) - &
                   c1 * ( f(1,m,k+1) - f(1,m,k-1) ) - c2 * ( f(1,m,k+2) - f(1,m,k-2) )
  enddo

  ! left boundary
  if ( kstrt == 1 ) then
    do j = 3, m-2
      df(1,j,1) =  b1 * f(3,j,1) + b2 * f(3,j,2) + b3 * f(3,j,3) + b4 * f(3,j,4) + b5 * f(3,j,5)
      df(2,j,1) = -c1 * ( f(3,j+1,1) - f(3,j-1,1) ) - c2 * ( f(3,j+2,1) - f(3,j-2,1) )
      df(3,j,1) =  c1 * ( f(2,j+1,1) - f(2,j-1,1) ) + c2 * ( f(2,j+2,1) - f(2,j-2,1) ) - &
                   b1 * f(1,j,1) - b2 * f(1,j,2) - b3 * f(1,j,3) - b4 * f(1,j,4) - b5 * f(1,j,5)
      df(1,j,2) =  a1 * f(3,j,1) + a2 * f(3,j,2) + a3 * f(3,j,3) + a4 * f(3,j,4) + a5 * f(3,j,5)
      df(2,j,2) = -c1 * ( f(3,j+1,2) - f(3,j-1,2) ) - c2 * ( f(3,j+2,2) - f(3,j-2,2) )
      df(3,j,2) =  c1 * ( f(2,j+1,2) - f(2,j-1,2) ) + c2 * ( f(2,j+2,2) - f(2,j-2,2) ) - &
                   a1 * f(1,j,1) - a2 * f(1,j,2) - a3 * f(1,j,3) - a4 * f(1,j,4) - a5 * f(1,j,5)
    enddo
    df(1,1,1) =  b1 * f(3,1,1) + b2 * f(3,1,2)   + b3 * f(3,1,3) + b4 * f(3,1,4) + b5 * f(3,1,5)
    df(2,1,1) = -b1 * f(3,1,1) - b2 * f(3,2,1)   - b3 * f(3,3,1) - b4 * f(3,4,1) - b5 * f(3,5,1)
    df(3,1,1) =  b1 * f(2,1,1) + b2 * f(2,2,1)   + b3 * f(2,3,1) + b4 * f(2,4,1) + b5 * f(2,5,1) - &
                 b1 * f(1,1,1) - b2 * f(1,1,2)   - b3 * f(1,1,3) - b4 * f(1,1,4) - b5 * f(1,1,5)
    df(1,2,1) =  b1 * f(3,2,1) + b2 * f(3,2,2)   + b3 * f(3,2,3) + b4 * f(3,2,4) + b5 * f(3,2,5)
    df(2,2,1) = -a1 * f(3,1,1) - a2 * f(3,2,1)   - a3 * f(3,3,1) - a4 * f(3,4,1) - a5 * f(3,5,1)
    df(3,2,1) =  a1 * f(2,1,1) + a2 * f(2,2,1)   + a3 * f(2,3,1) + a4 * f(2,4,1) + a5 * f(2,5,1) - &
                 b1 * f(1,2,1) - b2 * f(1,2,2)   - b3 * f(1,2,3) - b4 * f(1,2,4) - b5 * f(1,2,5)

    df(1,m-1,1) =  b1 * f(3,m-1,1) + b2 * f(3,m-1,2) + b3 * f(3,m-1,3) + b4 * f(3,m-1,4) + b5 * f(3,m-1,5)
    df(2,m-1,1) =  a1 * f(3,m,1)   + a2 * f(3,m-1,1) + a3 * f(3,m-2,1) + a4 * f(3,m-3,1) + a5 * f(3,m-4,1)
    df(3,m-1,1) = -a1 * f(2,m,1)   - a2 * f(2,m-1,1) - a3 * f(2,m-2,1) - a4 * f(2,m-3,1) - a5 * f(2,m-4,1) - &
                   b1 * f(1,m-1,1) - b2 * f(1,m-1,2) - b3 * f(1,m-1,3) - b4 * f(1,m-1,4) - b5 * f(1,m-1,5)
    df(1,m,1) =    b1 * f(3,m,1)   + b2 * f(3,m,2)   + b3 * f(3,m,3)   + b4 * f(3,m,4)   + b5 * f(3,m,5)
    df(2,m,1) =    b1 * f(3,m,1)   + b2 * f(3,m-1,1) + b3 * f(3,m-2,1) + b4 * f(3,m-3,1) + b5 * f(3,m-4,1)
    df(3,m,1) =   -b1 * f(2,m,1)   - b2 * f(2,m-1,1) - b3 * f(2,m-2,1) - b4 * f(2,m-3,1) - b5 * f(2,m-4,1) - &
                   b1 * f(1,m,1)   - b2 * f(1,m,2)   - b3 * f(1,m,3)   - b4 * f(1,m,4)   - b5 * f(1,m,5)

    df(1,1,2) =  a1 * f(3,1,1) + a2 * f(3,1,2)   + a3 * f(3,1,3) + a4 * f(3,1,4) + a5 * f(3,1,5)
    df(2,1,2) = -b1 * f(3,1,2) - b2 * f(3,2,2)   - b3 * f(3,3,2) - b4 * f(3,4,2) - b5 * f(3,5,2)
    df(3,1,2) =  b1 * f(2,1,2) + b2 * f(2,2,2)   + b3 * f(2,3,2) + b4 * f(2,4,2) + b5 * f(2,5,2) - &
                 a1 * f(1,1,1) - a2 * f(1,1,2)   - a3 * f(1,1,3) - a4 * f(1,1,4) - a5 * f(1,1,5)
    df(1,2,2) =  a1 * f(3,2,1) + a2 * f(3,2,2)   + a3 * f(3,2,3) + a4 * f(3,2,4) + a5 * f(3,2,5)
    df(2,2,2) = -a1 * f(3,1,2) - a2 * f(3,2,2)   - a3 * f(3,3,2) - a4 * f(3,4,2) - a5 * f(3,5,2)
    df(3,2,2) =  a1 * f(2,1,2) + a2 * f(2,2,2)   + a3 * f(2,3,2) + a4 * f(2,4,2) + a5 * f(2,5,2) - &
                 a1 * f(1,2,1) - a2 * f(1,2,2)   - a3 * f(1,2,3) - a4 * f(1,2,4) - a5 * f(1,2,5)

    df(1,m-1,2) =  a1 * f(3,m-1,1) + a2 * f(3,m-1,2) + a3 * f(3,m-1,3) + a4 * f(3,m-1,4) + a5 * f(3,m-1,5)
    df(2,m-1,2) =  a1 * f(3,m,2)   + a2 * f(3,m-1,2) + a3 * f(3,m-2,2) + a4 * f(3,m-3,2) + a5 * f(3,m-4,2)
    df(3,m-1,2) = -a1 * f(2,m,2)   - a2 * f(2,m-1,2) - a3 * f(2,m-2,2) - a4 * f(2,m-3,2) - a5 * f(2,m-4,2) - &
                   a1 * f(1,m-1,1) - a2 * f(1,m-1,2) - a3 * f(1,m-1,3) - a4 * f(1,m-1,4) - a5 * f(1,m-1,5)
    df(1,m,2) =    a1 * f(3,m,1)   + a2 * f(3,m,2)   + a3 * f(3,m,3)   + a4 * f(3,m,4)   + a5 * f(3,m,5)
    df(2,m,2) =    b1 * f(3,m,2)   + b2 * f(3,m-1,2) + b3 * f(3,m-2,2) + b4 * f(3,m-3,2) + b5 * f(3,m-4,2)
    df(3,m,2) =   -b1 * f(2,m,2)   - b2 * f(2,m-1,2) - b3 * f(2,m-2,2) - b4 * f(2,m-3,2) - b5 * f(2,m-4,2) - &
                   a1 * f(1,m,1)   - a2 * f(1,m,2)   - a3 * f(1,m,3)   - a4 * f(1,m,4)   - a5 * f(1,m,5)
  endif

  ! right boundary
  if ( kstrt == nvp ) then
    do j = 3, m-2
      df(1,j,n-1) = -a1 * f(3,j,n) - a2 * f(3,j,n-1) - a3 * f(3,j,n-2) - a4 * f(3,j,n-3) - a5 * f(3,j,n-4)
      df(2,j,n-1) = -c1 * ( f(3,j+1,n-1) - f(3,j-1,n-1) ) - c2 * ( f(3,j+2,n-1) - f(3,j-2,n-1) )
      df(3,j,n-1) =  c1 * ( f(2,j+1,n-1) - f(2,j-1,n-1) ) + c2 * ( f(2,j+2,n-1) - f(2,j-2,n-1) ) + &
                     a1 * f(1,j,n) + a2 * f(1,j,n-1) + a3 * f(1,j,n-2) + a4 * f(1,j,n-3) + a5 * f(1,j,n-4)
      df(1,j,n)   = -b1 * f(3,j,n) - b2 * f(3,j,n-1) - b3 * f(3,j,n-2) - b4 * f(3,j,n-3) - b5 * f(3,j,n-4)
      df(2,j,n)   = -c1 * ( f(3,j+1,n) - f(3,j-1,n) ) - c2 * ( f(3,j+2,n) - f(3,j-2,n) )
      df(3,j,n)   =  c1 * ( f(2,j+1,n) - f(2,j-1,n) ) + c2 * ( f(2,j+2,n) - f(2,j-2,n) ) + &
                     b1 * f(1,j,n) + b2 * f(1,j,n-1) + b3 * f(1,j,n-2) + b4 * f(1,j,n-3) + b5 * f(1,j,n-4)
    enddo
    df(1,1,n-1) = -a1 * f(3,1,n)   - a2 * f(3,1,n-1) - a3 * f(3,1,n-2) - a4 * f(3,1,n-3) - a5 * f(3,1,n-4)
    df(2,1,n-1) = -b1 * f(3,1,n-1) - b2 * f(3,2,n-1) - b3 * f(3,3,n-1) - b4 * f(3,4,n-1) - b5 * f(3,5,n-1)
    df(3,1,n-1) =  b1 * f(2,1,n-1) + b2 * f(2,2,n-1) + b3 * f(2,3,n-1) + b4 * f(2,4,n-1) + b5 * f(2,5,n-1) + &
                   a1 * f(1,1,n)   + a2 * f(1,1,n-1) + a3 * f(1,1,n-2) + a4 * f(1,1,n-3) + a5 * f(1,1,n-4)
    df(1,2,n-1) = -a1 * f(3,2,n)   - a2 * f(3,2,n-1) - a3 * f(3,2,n-2) - a4 * f(3,2,n-3) - a5 * f(3,2,n-4)
    df(2,2,n-1) = -a1 * f(3,1,n-1) - a2 * f(3,2,n-1) - a3 * f(3,3,n-1) - a4 * f(3,4,n-1) - a5 * f(3,5,n-1)
    df(3,2,n-1) =  a1 * f(2,1,n-1) + a2 * f(2,2,n-1) + a3 * f(2,3,n-1) + a4 * f(2,4,n-1) + a5 * f(2,5,n-1) + &
                   a1 * f(1,2,n)   + a2 * f(1,2,n-1) + a3 * f(1,2,n-2) + a4 * f(1,2,n-3) + a5 * f(1,2,n-4)

    df(1,m-1,n-1) = -a1 * f(3,m-1,n) - a2 * f(3,m-1,n-1) - a3 * f(3,m-1,n-2) - a4 * f(3,m-1,n-3) - a5 * f(3,m-1,n-4)
    df(2,m-1,n-1) =  a1 * f(3,m,n-1) + a2 * f(3,m-1,n-1) + a3 * f(3,m-2,n-1) + a4 * f(3,m-3,n-1) + a5 * f(3,m-4,n-1)
    df(3,m-1,n-1) = -a1 * f(2,m,n-1) - a2 * f(2,m-1,n-1) - a3 * f(2,m-2,n-1) - a4 * f(2,m-3,n-1) - a5 * f(2,m-4,n-1) + &
                     a1 * f(1,m-1,n) + a2 * f(1,m-1,n-1) + a3 * f(1,m-1,n-2) + a4 * f(1,m-1,n-3) + a5 * f(1,m-1,n-4)
    df(1,m,n-1) =   -a1 * f(3,m,n)   - a2 * f(3,m,n-1)   - a3 * f(3,m,n-2)   - a4 * f(3,m,n-3)   - a5 * f(3,m,n-4)
    df(2,m,n-1) =    b1 * f(3,m,n-1) + b2 * f(3,m-1,n-1) + b3 * f(3,m-2,n-1) + b4 * f(3,m-3,n-1) + b5 * f(3,m-4,n-1)
    df(3,m,n-1) =   -b1 * f(2,m,n-1) - b2 * f(2,m-1,n-1) - b3 * f(2,m-2,n-1) - b4 * f(2,m-3,n-1) - b5 * f(2,m-4,n-1) + &
                     a1 * f(1,m,n)   + a2 * f(1,m,n-1)   + a3 * f(1,m,n-2)   + a4 * f(1,m,n-3)   + a5 * f(1,m,n-4)

    df(1,1,n) = -b1 * f(3,1,n) - b2 * f(3,1,n-1) - b3 * f(3,1,n-2) - b4 * f(3,1,n-3) - b5 * f(3,1,n-4)
    df(2,1,n) = -b1 * f(3,1,n) - b2 * f(3,2,n)   - b3 * f(3,3,n)   - b4 * f(3,4,n)   - b5 * f(3,5,n)  
    df(3,1,n) =  b1 * f(2,1,n) + b2 * f(2,2,n)   + b3 * f(2,3,n)   + b4 * f(2,4,n)   + b5 * f(2,5,n) + &
                 b1 * f(1,1,n) + b2 * f(1,1,n-1) + b3 * f(1,1,n-2) + b4 * f(1,1,n-3) + b5 * f(1,1,n-4)
    df(1,2,n) = -b1 * f(3,2,n) - b2 * f(3,2,n-1) - b3 * f(3,2,n-2) - b4 * f(3,2,n-3) - b5 * f(3,2,n-4)
    df(2,2,n) = -a1 * f(3,1,n) - a2 * f(3,2,n)   - a3 * f(3,3,n)   - a4 * f(3,4,n)   - a5 * f(3,5,n)  
    df(3,2,n) =  a1 * f(2,1,n) + a2 * f(2,2,n)   + a3 * f(2,3,n)   + a4 * f(2,4,n)   + a5 * f(2,5,n) + &
                 b1 * f(1,2,n) + b2 * f(1,2,n-1) + b3 * f(1,2,n-2) + b4 * f(1,2,n-3) + b5 * f(1,2,n-4)

    df(1,m-1,n) = -b1 * f(3,m-1,n) - b2 * f(3,m-1,n-1) - b3 * f(3,m-1,n-2) - b4 * f(3,m-1,n-3) - b5 * f(3,m-1,n-4)
    df(2,m-1,n) =  a1 * f(3,m,n)   + a2 * f(3,m-1,n)   + a3 * f(3,m-2,n)   + a4 * f(3,m-3,n)   + a5 * f(3,m-4,n)  
    df(3,m-1,n) = -a1 * f(2,m,n)   - a2 * f(2,m-1,n)   - a3 * f(2,m-2,n)   - a4 * f(2,m-3,n)   - a5 * f(2,m-4,n) + &
                   b1 * f(1,m-1,n) + b2 * f(1,m-1,n-1) + b3 * f(1,m-1,n-2) + b4 * f(1,m-1,n-3) + b5 * f(1,m-1,n-4)
    df(1,m,n) =   -b1 * f(3,m,n)   - b2 * f(3,m,n-1)   - b3 * f(3,m,n-2)   - b4 * f(3,m,n-3)   - b5 * f(3,m,n-4)
    df(2,m,n) =    b1 * f(3,m,n)   + b2 * f(3,m-1,n)   + b3 * f(3,m-2,n)   + b4 * f(3,m-3,n)   + b5 * f(3,m-4,n)
    df(3,m,n) =   -b1 * f(2,m,n)   - b2 * f(2,m-1,n)   - b3 * f(2,m-2,n)   - b4 * f(2,m-3,n)   - b5 * f(2,m-4,n) + &
                   b1 * f(1,m,n)   + b2 * f(1,m,n-1)   + b3 * f(1,m,n-2)   + b4 * f(1,m,n-3)   + b5 * f(1,m,n-4)
    

  endif

end subroutine curlf_4order

subroutine bfqp_src_2order( cu, dcu, amu, bxy, src, m, n, dx, kstrt, nvp )

  implicit none

  real, intent(in), dimension(:,:,:), pointer :: cu, dcu, amu, bxy
  real, intent(inout), dimension(:,:,:), pointer :: src
  integer, intent(in) :: m, n, kstrt, nvp
  real, intent(in) :: dx

  ! local variables
  integer :: i, j, k, kl, ku
  real :: inv_dxh, a1, a2, a3

  inv_dxh = 0.5 / dx
  a1 = -3.0 * inv_dxh
  a2 =  4.0 * inv_dxh
  a3 = -1.0 * inv_dxh

  kl = 1; ku = n
  if ( kstrt == 1 ) kl = 2
  if ( kstrt == nvp ) ku = n-1

  ! interior
  do k = kl, ku
    do j = 2, m-1
      src(1,j,k) = -dcu(2,j,k) - bxy(1,j,k) + inv_dxh * ( amu(3,j+1,k) - amu(3,j-1,k) + amu(2,j,k+1) - &
                    amu(2,j,k-1) ) - inv_dxh * ( cu(3,j,k+1) - cu(3,j,k-1) )
      src(2,j,k) =  dcu(1,j,k) - bxy(2,j,k) - inv_dxh * ( amu(1,j+1,k) - amu(1,j-1,k) + amu(3,j,k+1) - &
                    amu(3,j,k-1) ) + inv_dxh * ( cu(3,j+1,k) - cu(3,j-1,k) )
    enddo
  enddo

  ! top and bottom boundaries
  do k = kl, ku
    src(1,1,k) = -dcu(2,1,k) - bxy(1,1,k) + a1 * amu(3,1,k) + a2 * amu(3,2,k)   + a3 * amu(3,3,k)   + &
                  inv_dxh * ( amu(2,1,k+1) - amu(2,1,k-1) ) - inv_dxh * ( cu(3,1,k+1) - cu(3,1,k-1) )
    src(2,1,k) =  dcu(1,1,k) - bxy(2,1,k) - a1 * amu(1,1,k) - a2 * amu(1,2,k)   - a3 * amu(1,3,k)   - &
                  inv_dxh * ( amu(3,1,k+1) - amu(3,1,k-1) ) + a1 * cu(3,1,k) + a2 * cu(3,2,k) + a3 * cu(3,3,k)
    src(1,m,k) = -dcu(2,m,k) - bxy(1,m,k) - a1 * amu(3,m,k) - a2 * amu(3,m-1,k) - a3 * amu(3,m-2,k) + &
                  inv_dxh * ( amu(2,m,k+1) - amu(2,m,k-1) ) - inv_dxh * ( cu(3,m,k+1) - cu(3,m,k-1) )
    src(2,m,k) =  dcu(1,m,k) - bxy(2,m,k) + a1 * amu(1,m,k) + a2 * amu(1,m-1,k) + a3 * amu(1,m-2,k) - &
                  inv_dxh * ( amu(3,m,k+1) - amu(3,m,k-1) ) - a1 * cu(3,m,k) - a2 * cu(3,m-1,k) - a3 * cu(3,m-2,k)
  enddo

  ! left boundary
  if ( kstrt == 1 ) then
    do j = 2, m-1
      src(1,j,1) = -dcu(2,j,1) - bxy(1,j,1) + inv_dxh * ( amu(3,j+1,1) - amu(3,j-1,1) ) + &
                    a1 * amu(2,j,1) + a2 * amu(2,j,2) + a3 * amu(2,j,3) - a1 * cu(3,j,1) - a2 * cu(3,j,2) - a3 * cu(3,j,3)
      src(2,j,1) =  dcu(1,j,1) - bxy(2,j,1) - inv_dxh * ( amu(1,j+1,1) - amu(1,j-1,1) ) - &
                    a1 * amu(3,j,1) - a2 * amu(3,j,2) - a3 * amu(3,j,3) + inv_dxh * ( cu(3,j+1,1) - cu(3,j-1,1) )
    enddo
    src(1,1,1) = -dcu(2,1,1) - bxy(1,1,1) + &
                  a1 * ( amu(3,1,1) + amu(2,1,1) - cu(3,1,1) ) + &
                  a2 * ( amu(3,2,1) + amu(2,1,2) - cu(3,1,2) ) + &
                  a3 * ( amu(3,3,1) + amu(2,1,3) - cu(3,1,3) )
    src(2,1,1) =  dcu(1,1,1) - bxy(2,1,1) - &
                  a1 * ( amu(1,1,1) + amu(3,1,1) - cu(3,1,1) ) - &
                  a2 * ( amu(1,2,1) + amu(3,1,2) - cu(3,2,1) ) - &
                  a3 * ( amu(1,3,1) + amu(3,1,3) - cu(3,3,1) )
    src(1,m,1) = -dcu(2,m,1) - bxy(1,m,1) - &
                  a1 * ( amu(3,m,1)   - amu(2,m,1) + cu(3,m,1) ) - &
                  a2 * ( amu(3,m-1,1) - amu(2,m,2) + cu(3,m,2) ) - &
                  a3 * ( amu(3,m-2,1) - amu(2,m,3) + cu(3,m,3) )
    src(2,m,1) =  dcu(1,m,1) - bxy(2,m,1) + &
                  a1 * ( amu(1,m,1)   - amu(3,m,1) - cu(3,m,1) ) + &
                  a2 * ( amu(1,m-1,1) - amu(3,m,2) - cu(3,m-1,1) ) + &
                  a3 * ( amu(1,m-2,1) - amu(3,m,3) - cu(3,m-2,1) )
  endif

  ! right boundary
  if ( kstrt == nvp ) then
    do j = 2, m-1
      src(1,j,n) = -dcu(2,j,n) - bxy(1,j,n) + inv_dxh * ( amu(3,j+1,n) - amu(3,j-1,n) ) - &
                    a1 * amu(2,j,n) - a2 * amu(2,j,n-1) - a3 * amu(2,j,n-2) + a1 * cu(3,j,n) + a2 * cu(3,j,n-1) + a3 * cu(3,j,n-2)
      src(2,j,n) =  dcu(1,j,n) - bxy(2,j,n) - inv_dxh * ( amu(1,j+1,n) - amu(1,j-1,n) ) + &
                    a1 * amu(3,j,n) + a2 * amu(3,j,n-1) + a3 * amu(3,j,n-2) + inv_dxh * ( cu(3,j+1,n) - cu(3,j-1,n) )
    enddo

    src(1,1,n) = -dcu(2,1,n) - bxy(1,1,n) + &
                  a1 * ( amu(3,1,n) - amu(2,1,n)   + cu(3,1,n) ) + &
                  a2 * ( amu(3,2,n) - amu(2,1,n-1) + cu(3,1,n-1) ) + &
                  a3 * ( amu(3,3,n) - amu(2,1,n-2) + cu(3,1,n-2) )
    src(2,1,n) =  dcu(1,1,n) - bxy(2,1,n) - &
                  a1 * ( amu(1,1,n) - amu(3,1,n) - cu(3,1,n) ) - &
                  a2 * ( amu(1,2,n) - amu(3,1,n-1) - cu(3,2,n) ) - &
                  a3 * ( amu(1,3,n) - amu(3,1,n-2) - cu(3,3,n) )
    src(1,m,n) = -dcu(2,m,n) - bxy(1,m,n) - &
                  a1 * ( amu(3,m,n)   + amu(2,m,n)   - cu(3,m,n) ) - &
                  a2 * ( amu(3,m-1,n) + amu(2,m,n-1) - cu(3,m,n-1) ) - &
                  a3 * ( amu(3,m-2,n) + amu(2,m,n-2) - cu(3,m,n-2) )
    src(2,m,n) =  dcu(1,m,n) - bxy(2,m,n) + &
                  a1 * ( amu(1,m,n)   + amu(3,m,n)   - cu(3,m,n) ) + &
                  a2 * ( amu(1,m-1,n) + amu(3,m,n-1) - cu(3,m-1,n) ) + &
                  a3 * ( amu(1,m-2,n) + amu(3,m,n-2) - cu(3,m-2,n) )
  endif

end subroutine bfqp_src_2order

subroutine bfqp_src_2order_new( cu, dcu, amu, src, m, n, dx, kstrt, nvp )

  implicit none

  real, intent(in), dimension(:,:,:), pointer :: cu, dcu, amu
  real, intent(inout), dimension(:,:,:), pointer :: src
  integer, intent(in) :: m, n, kstrt, nvp
  real, intent(in) :: dx

  ! local variables
  integer :: i, j, k, kl, ku
  real :: inv_dxh, a1, a2, a3

  inv_dxh = 0.5 / dx
  a1 = -3.0 * inv_dxh
  a2 =  4.0 * inv_dxh
  a3 = -1.0 * inv_dxh

  kl = 1; ku = n
  if ( kstrt == 1 ) kl = 2
  if ( kstrt == nvp ) ku = n-1
  ! write(*,*) "m=", m
  ! write(*,*) "n=", n
  ! write(*,*) "kstrt=", kstrt
  ! write(*,*) "ku=", ku
  ! interior
  do k = kl, ku
    do j = 2, m-1
      src(1,j,k) = -dcu(2,j,k) + inv_dxh * ( amu(3,j+1,k) - amu(3,j-1,k) + amu(2,j,k+1) - &
                    amu(2,j,k-1) ) - inv_dxh * ( cu(3,j,k+1) - cu(3,j,k-1) )
      src(2,j,k) =  dcu(1,j,k) - inv_dxh * ( amu(1,j+1,k) - amu(1,j-1,k) + amu(3,j,k+1) - &
                    amu(3,j,k-1) ) + inv_dxh * ( cu(3,j+1,k) - cu(3,j-1,k) )
    enddo
  enddo

  ! top and bottom boundaries
  do k = kl, ku
    src(1,1,k) = -dcu(2,1,k) + a1 * amu(3,1,k) + a2 * amu(3,2,k)   + a3 * amu(3,3,k)   + &
                  inv_dxh * ( amu(2,1,k+1) - amu(2,1,k-1) ) - inv_dxh * ( cu(3,1,k+1) - cu(3,1,k-1) )
    src(2,1,k) =  dcu(1,1,k) - a1 * amu(1,1,k) - a2 * amu(1,2,k)   - a3 * amu(1,3,k)   - &
                  inv_dxh * ( amu(3,1,k+1) - amu(3,1,k-1) ) + a1 * cu(3,1,k) + a2 * cu(3,2,k) + a3 * cu(3,3,k)
    src(1,m,k) = -dcu(2,m,k) - a1 * amu(3,m,k) - a2 * amu(3,m-1,k) - a3 * amu(3,m-2,k) + &
                  inv_dxh * ( amu(2,m,k+1) - amu(2,m,k-1) ) - inv_dxh * ( cu(3,m,k+1) - cu(3,m,k-1) )
    src(2,m,k) =  dcu(1,m,k) + a1 * amu(1,m,k) + a2 * amu(1,m-1,k) + a3 * amu(1,m-2,k) - &
                  inv_dxh * ( amu(3,m,k+1) - amu(3,m,k-1) ) - a1 * cu(3,m,k) - a2 * cu(3,m-1,k) - a3 * cu(3,m-2,k)
  enddo

  ! left boundary
  if ( kstrt == 1 ) then
    do j = 2, m-1
      src(1,j,1) = -dcu(2,j,1) + inv_dxh * ( amu(3,j+1,1) - amu(3,j-1,1) ) + &
                    a1 * amu(2,j,1) + a2 * amu(2,j,2) + a3 * amu(2,j,3) - a1 * cu(3,j,1) - a2 * cu(3,j,2) - a3 * cu(3,j,3)
      src(2,j,1) =  dcu(1,j,1) - inv_dxh * ( amu(1,j+1,1) - amu(1,j-1,1) ) - &
                    a1 * amu(3,j,1) - a2 * amu(3,j,2) - a3 * amu(3,j,3) + inv_dxh * ( cu(3,j+1,1) - cu(3,j-1,1) )
    enddo
    src(1,1,1) = -dcu(2,1,1) + &
                  a1 * ( amu(3,1,1) + amu(2,1,1) - cu(3,1,1) ) + &
                  a2 * ( amu(3,2,1) + amu(2,1,2) - cu(3,1,2) ) + &
                  a3 * ( amu(3,3,1) + amu(2,1,3) - cu(3,1,3) )
    src(2,1,1) =  dcu(1,1,1) - &
                  a1 * ( amu(1,1,1) + amu(3,1,1) - cu(3,1,1) ) - &
                  a2 * ( amu(1,2,1) + amu(3,1,2) - cu(3,2,1) ) - &
                  a3 * ( amu(1,3,1) + amu(3,1,3) - cu(3,3,1) )
    src(1,m,1) = -dcu(2,m,1) - &
                  a1 * ( amu(3,m,1)   - amu(2,m,1) + cu(3,m,1) ) - &
                  a2 * ( amu(3,m-1,1) - amu(2,m,2) + cu(3,m,2) ) - &
                  a3 * ( amu(3,m-2,1) - amu(2,m,3) + cu(3,m,3) )
    src(2,m,1) =  dcu(1,m,1) + &
                  a1 * ( amu(1,m,1)   - amu(3,m,1) - cu(3,m,1) ) + &
                  a2 * ( amu(1,m-1,1) - amu(3,m,2) - cu(3,m-1,1) ) + &
                  a3 * ( amu(1,m-2,1) - amu(3,m,3) - cu(3,m-2,1) )
  endif

  ! right boundary
  if ( kstrt == nvp ) then
    do j = 2, m-1
      src(1,j,n) = -dcu(2,j,n) + inv_dxh * ( amu(3,j+1,n) - amu(3,j-1,n) ) - &
                    a1 * amu(2,j,n) - a2 * amu(2,j,n-1) - a3 * amu(2,j,n-2) + a1 * cu(3,j,n) + a2 * cu(3,j,n-1) + a3 * cu(3,j,n-2)
      src(2,j,n) =  dcu(1,j,n) - inv_dxh * ( amu(1,j+1,n) - amu(1,j-1,n) ) + &
                    a1 * amu(3,j,n) + a2 * amu(3,j,n-1) + a3 * amu(3,j,n-2) + inv_dxh * ( cu(3,j+1,n) - cu(3,j-1,n) )
    enddo

    src(1,1,n) = -dcu(2,1,n) + &
                  a1 * ( amu(3,1,n) - amu(2,1,n)   + cu(3,1,n) ) + &
                  a2 * ( amu(3,2,n) - amu(2,1,n-1) + cu(3,1,n-1) ) + &
                  a3 * ( amu(3,3,n) - amu(2,1,n-2) + cu(3,1,n-2) )
    src(2,1,n) =  dcu(1,1,n) - &
                  a1 * ( amu(1,1,n) - amu(3,1,n) - cu(3,1,n) ) - &
                  a2 * ( amu(1,2,n) - amu(3,1,n-1) - cu(3,2,n) ) - &
                  a3 * ( amu(1,3,n) - amu(3,1,n-2) - cu(3,3,n) )
    src(1,m,n) = -dcu(2,m,n) - &
                  a1 * ( amu(3,m,n)   + amu(2,m,n)   - cu(3,m,n) ) - &
                  a2 * ( amu(3,m-1,n) + amu(2,m,n-1) - cu(3,m,n-1) ) - &
                  a3 * ( amu(3,m-2,n) + amu(2,m,n-2) - cu(3,m,n-2) )
    src(2,m,n) =  dcu(1,m,n) + &
                  a1 * ( amu(1,m,n)   + amu(3,m,n)   - cu(3,m,n) ) + &
                  a2 * ( amu(1,m-1,n) + amu(3,m,n-1) - cu(3,m-1,n) ) + &
                  a3 * ( amu(1,m-2,n) + amu(3,m,n-2) - cu(3,m-2,n) )
  endif

end subroutine bfqp_src_2order_new

end module fd2d_lib
