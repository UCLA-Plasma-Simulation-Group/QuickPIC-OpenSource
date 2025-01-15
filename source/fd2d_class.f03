!> author: Fei Li
!> date: 05/17/2018
!> summary: Finite difference class of two dimensional fields
module fd2d_class

use perrors_class
use parallel_pipe_class
use ufield2d_class
use fd2d_lib

implicit none

private

integer, parameter :: p_fd_2order = 1, p_fd_4order = 2

public :: fd2d
public :: p_fd_2order, p_fd_4order

type fd2d

  private

  class(perrors), pointer, public :: err => null()
  class(parallel_pipe), pointer, public :: p => null()
  integer :: kind

  contains

  generic :: new => init_fd2d
  generic :: del => end_fd2d
  generic :: divf => select_divf
  generic :: divf_dyad => select_divf_dyad
  generic :: gradf => select_gradf
  generic :: curlf => select_curlf
  generic :: curlfz => select_curlfz
  generic :: bfqp_src => select_bfqp_src
  generic :: bfqp_src_new => select_bfqp_src_new
  generic :: smoothf => select_smoothf
  procedure, private :: init_fd2d
  procedure, private :: end_fd2d
  procedure, private :: select_gradf, select_divf, select_curlf, select_divf_dyad, select_bfqp_src, &
                        select_curlfz, select_smoothf, select_bfqp_src_new

end type fd2d

character(len=10), save :: class = 'fd2d:'
character(len=128), save :: erstr

contains

subroutine init_fd2d( this, pp, perr, kind )

  implicit none

  class(fd2d), intent(inout) :: this
  class(parallel_pipe), intent(in), pointer :: pp
  class(perrors), intent(in), pointer :: perr
  integer, intent(in) :: kind

  ! local variables
  character(len=18), save :: sname = 'init_fd2d:'
  
  call perr%werrfl2( class//sname//' started' )

  this%err => perr
  this%p => pp
  this%kind = kind

  call perr%werrfl2( class//sname//' ended' )

end subroutine init_fd2d

subroutine end_fd2d( this )

  implicit none

  class(fd2d), intent(inout) :: this

  ! local variables
  character(len=18), save :: sname = 'end_fd2d:'

  call this%err%werrfl2( class//sname//' started' )
  ! do nothing
  call this%err%werrfl2( class//sname//' ended' )

end subroutine end_fd2d

subroutine select_divf( this, rspace, dspace, dx, order )

  implicit none

  class(fd2d), intent(in) :: this
  class(ufield2d), intent(inout) :: rspace, dspace
  real, intent(in) :: dx
  integer, intent(in) :: order

  ! local variables
  integer :: kstrt, nvp, nd1p, nd2p
  real, dimension(:,:,:), pointer :: f => null(), df => null()
  character(len=18), save :: sname = 'select_divf:'

  call this%err%werrfl2(class//sname//' started')

  kstrt = this%p%getlkstrt()
  nvp = this%p%getlnvp()
  f => rspace%getrf_gc()
  df => dspace%getrf_gc()
  nd1p = rspace%getnd1p()
  nd2p = rspace%getnd2p()

  select case ( order )
  case ( p_fd_2order )
    call divf_2order( f, df, nd1p, nd2p, dx, kstrt, nvp )
  case ( p_fd_4order )
    call divf_4order( f, df, nd1p, nd2p, dx, kstrt, nvp )
  case default
    call this%err%equit(class//sname//' invalid order')
  end select

  call this%err%werrfl2(class//sname//' ended')

end subroutine select_divf

subroutine select_divf_dyad( this, rspace, dspace, dx, order )

  implicit none

  class(fd2d), intent(in) :: this
  class(ufield2d), intent(inout) :: rspace, dspace
  real, intent(in) :: dx
  integer, intent(in) :: order

  ! local variables
  integer :: kstrt, nvp, nd1p, nd2p
  real, dimension(:,:,:), pointer :: f => null(), df => null()
  character(len=18), save :: sname = 'select_divf_dyad:'

  call this%err%werrfl2(class//sname//' started')

  kstrt = this%p%getlkstrt()
  nvp = this%p%getlnvp()
  f => rspace%getrf_gc()
  df => dspace%getrf_gc()
  nd1p = rspace%getnd1p()
  nd2p = rspace%getnd2p()

  select case ( order )
  case ( p_fd_2order )
    call divf_dyad_2order( f, df, nd1p, nd2p, dx, kstrt, nvp )
  case ( p_fd_4order )
    ! call divf_dyad_4order( f, df, nd1p, nd2p, dx, kstrt, nvp )
    call this%err%equit(class//sname//' divf_dyad_4order not implemented')
  case default
    call this%err%equit(class//sname//' invalid order')
  end select

  call this%err%werrfl2(class//sname//' ended')

end subroutine select_divf_dyad

subroutine select_gradf( this, rspace, dspace, dx, order )

  implicit none

  class(fd2d), intent(in) :: this
  class(ufield2d), intent(inout) :: rspace, dspace
  real, intent(in) :: dx
  integer, intent(in) :: order

  ! local variables
  integer :: kstrt, nvp, nd1p, nd2p
  real, dimension(:,:,:), pointer :: f => null(), df => null()
  character(len=18), save :: sname = 'select_gradf:'

  call this%err%werrfl2(class//sname//' started')

  kstrt = this%p%getlkstrt()
  nvp = this%p%getlnvp()
  f => rspace%getrf_gc()
  df => dspace%getrf_gc()
  nd1p = rspace%getnd1p()
  nd2p = rspace%getnd2p()

  select case ( order )
  case ( p_fd_2order )
    call gradf_2order( f, df, nd1p, nd2p, dx, kstrt, nvp )
  case ( p_fd_4order )
    call gradf_4order( f, df, nd1p, nd2p, dx, kstrt, nvp )
  case default
    call this%err%equit(class//sname//' invalid order')
  end select

  call this%err%werrfl2(class//sname//' ended')

end subroutine select_gradf

subroutine select_curlf( this, rspace, dspace, dx, order )

  implicit none

  class(fd2d), intent(in) :: this
  class(ufield2d), intent(inout) :: rspace, dspace
  real, intent(in) :: dx
  integer, intent(in) :: order

  ! local variables
  integer :: kstrt, nvp, nd1p, nd2p
  real, dimension(:,:,:), pointer :: f => null(), df => null()
  character(len=18), save :: sname = 'select_curlf:'

  call this%err%werrfl2(class//sname//' started')

  kstrt = this%p%getlkstrt()
  nvp = this%p%getlnvp()
  f => rspace%getrf_gc()
  df => dspace%getrf_gc()
  nd1p = rspace%getnd1p()
  nd2p = rspace%getnd2p()

  select case ( order )
  case ( p_fd_2order )
    call curlf_2order( f, df, nd1p, nd2p, dx, kstrt, nvp )
  case ( p_fd_4order )
    call curlf_4order( f, df, nd1p, nd2p, dx, kstrt, nvp )
  case default
    call this%err%equit(class//sname//' invalid order')
  end select

  call this%err%werrfl2(class//sname//' ended')

end subroutine select_curlf

subroutine select_curlfz( this, rspace, dspace, dx, order )

  implicit none

  class(fd2d), intent(in) :: this
  class(ufield2d), intent(inout) :: rspace, dspace
  real, intent(in) :: dx
  integer, intent(in) :: order

  ! local variables
  integer :: kstrt, nvp, nd1p, nd2p
  real, dimension(:,:,:), pointer :: f => null(), df => null()
  character(len=18), save :: sname = 'select_curlfz:'

  call this%err%werrfl2(class//sname//' started')

  kstrt = this%p%getlkstrt()
  nvp = this%p%getlnvp()
  f => rspace%getrf_gc()
  df => dspace%getrf_gc()
  nd1p = rspace%getnd1p()
  nd2p = rspace%getnd2p()

  select case ( order )
  case ( p_fd_2order )
    call curlfz_2order( f, df, nd1p, nd2p, dx, kstrt, nvp )
  case ( p_fd_4order )
    ! call curlf_4order( f, df, nd1p, nd2p, dx, kstrt, nvp )
  case default
    call this%err%equit(class//sname//' invalid order')
  end select

  call this%err%werrfl2(class//sname//' ended')

end subroutine select_curlfz

subroutine select_bfqp_src( this, cu, dcu, amu, bxy, src, dx, order )

  implicit none

  class(fd2d), intent(in) :: this
  class(ufield2d), intent(in) :: cu, dcu, amu, bxy
  class(ufield2d), intent(inout) :: src
  real, intent(in) :: dx
  integer, intent(in) :: order

  ! local variables
  integer :: kstrt, nvp, nd1p, nd2p
  real, dimension(:,:,:), pointer :: pcu => null(), pdcu => null(), &
                                     pamu => null(), pbxy => null(), psrc => null()

  character(len=18), save :: sname = 'select_bfqp_src:'

  call this%err%werrfl2(class//sname//' started')

  kstrt = this%p%getlkstrt()
  nvp   = this%p%getlnvp()
  pcu  => cu%getrf_gc()
  pdcu => dcu%getrf_gc()
  pamu => amu%getrf_gc()
  pbxy => bxy%getrf_gc()
  psrc => src%getrf_gc()
  nd1p = cu%getnd1p()
  nd2p = cu%getnd2p()

  select case ( order )
  case ( p_fd_2order )
    call bfqp_src_2order( pcu, pdcu, pamu, pbxy, psrc, nd1p, nd2p, dx, kstrt, nvp )
  case ( p_fd_4order )
    ! call divf_dyad_4order( f, df, nd1p, nd2p, dx, kstrt, nvp )
    call this%err%equit(class//sname//' bfqp_src_4order not implemented')
  case default
    call this%err%equit(class//sname//' invalid order')
  end select

  call this%err%werrfl2(class//sname//' ended')

end subroutine select_bfqp_src

subroutine select_bfqp_src_new( this, cu, dcu, amu, src, dx, order )

  implicit none

  class(fd2d), intent(in) :: this
  class(ufield2d), intent(in) :: cu, dcu, amu
  class(ufield2d), intent(inout) :: src
  real, intent(in) :: dx
  integer, intent(in) :: order

  ! local variables
  integer :: kstrt, nvp, nd1p, nd2p
  real, dimension(:,:,:), pointer :: pcu => null(), pdcu => null(), &
                                     pamu => null(), psrc => null()

  character(len=18), save :: sname = 'select_bfqp_src_new:'
  
  call this%err%werrfl2(class//sname//' started')

  kstrt = this%p%getlkstrt()
  nvp   = this%p%getlnvp()
  pcu  => cu%getrf_gc()
  pdcu => dcu%getrf_gc()
  pamu => amu%getrf_gc()
  psrc => src%getrf_gc()
  nd1p = cu%getnd1p()
  nd2p = cu%getnd2p()

  select case ( order )
  case ( p_fd_2order )
    ! call bfqp_src_2order( pcu, pdcu, pamu, pbxy, psrc, nd1p, nd2p, dx, kstrt, nvp )
    call bfqp_src_2order_new( pcu, pdcu, pamu, psrc, nd1p, nd2p, dx, kstrt, nvp )
  case ( p_fd_4order )
    ! call divf_dyad_4order( f, df, nd1p, nd2p, dx, kstrt, nvp )
    call this%err%equit(class//sname//' bfqp_src_4order not implemented')
  case default
    call this%err%equit(class//sname//' invalid order')
  end select

  call this%err%werrfl2(class//sname//' ended')

end subroutine select_bfqp_src_new

subroutine select_smoothf( this, rspace, dspace, order )

  implicit none

  class(fd2d), intent(in) :: this
  class(ufield2d), intent(inout) :: rspace, dspace
  integer, intent(in) :: order

  ! local variables
  integer :: kstrt, nvp, nd1p, nd2p
  real, dimension(:,:,:), pointer :: f => null(), df => null()
  character(len=18), save :: sname = 'select_smoothf:'

  call this%err%werrfl2(class//sname//' started')

  kstrt = this%p%getlkstrt()
  nvp = this%p%getlnvp()
  f => rspace%getrf_gc()
  df => dspace%getrf_gc()
  nd1p = rspace%getnd1p()
  nd2p = rspace%getnd2p()

  select case ( order )
  case ( p_fd_2order )
    ! call smoothf_2order( f, df, nd1p, nd2p, kstrt, nvp )
    call this%err%equit(class//sname//'smoothf_2order not implemented')
  case ( p_fd_4order )
    ! call smoothf_4order( f, df, nd1p, nd2p, dx, kstrt, nvp )
    call this%err%equit(class//sname//'smoothf_4order not implemented')
  case default
    call this%err%equit(class//sname//' invalid order')
  end select

  call this%err%werrfl2(class//sname//' ended')

end subroutine select_smoothf

end module fd2d_class