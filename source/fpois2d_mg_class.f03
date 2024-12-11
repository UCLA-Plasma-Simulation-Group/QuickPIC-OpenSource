module fpois2d_mg_class

use perrors_class
use parallel_pipe_class
use spect2d_class
use fd2d_class, only : p_fd_2order, p_fd_4order
use ufield2d_class

implicit none

private

integer, parameter :: HYPRE_TYPE = 8
integer, parameter :: p_xdim = 2
integer, parameter :: p_bnd_dd = 1, &   ! Dirichlet-Dirichlet boundary condition
                      p_bnd_dn = 2, &   ! Dirichlet-Neumann boundary condition
                      p_bnd_nd = 3, &   ! Neumann-Dirichlet boundary condition
                      p_bnd_nn = 4      ! Neumann-Neumann boundary condition
integer, parameter :: p_lower = 1, p_upper = 2
! integer, parameter :: p_fd_2order = 1, p_fd_4order = 2
integer, parameter :: p_solver_smg = 1, &
                      p_solver_pfmg = 2, &
                      p_solver_pcg = 3
integer, parameter :: p_precond_smg = 1, &
                      p_precond_pfmg = 2, &
                      p_precond_ds = 3, &
                      p_precond_none = 99
integer, parameter :: p_stencil_poisson = 1, &
                      p_stencil_pseudo = 2, &
                      p_stencil_explicit = 3

character(len=16), save :: class = 'fpois2d_mg:'
character(len=128), save :: erstr

real, dimension(:), pointer :: values => null()

public :: fpois2d_mg, get_poismg2table

type :: fpois2d_mg

  private

  class( perrors ), pointer, public :: err => null()
  class( parallel_pipe ), pointer, public :: p => null()
  logical, dimension(4) :: phys_bnd_flag
  integer :: order

  ! HYPRE parameters
  integer, dimension(:,:), pointer :: offsets => null()
  integer, dimension(:), pointer :: stencil_idx => null()
  integer :: stencil_num
  integer :: solver_type
  integer :: precond_type
  ! integer :: maxiter, n_pre, n_post
  real :: tol
  ! matrices and solvers for different types of boundary conditions and stencil kind
  ! (1,*) - Dirichlet-Dirichlet
  ! (2,*) - Dirichlet-Neumann
  ! (3,*) - Neumann-Dirichlet
  ! (4,*) - Neumann-Neumann
  ! (*,1) - Poisson stencil
  ! (*,2) - pseudo-Poisson stencil
  ! (*,3) - explicit stencil
  ! integer(HYPRE_TYPE), dimension(4,2) :: A, solver, precond
  integer(HYPRE_TYPE), dimension(4,3) :: A, solver, precond
  integer(HYPRE_TYPE) :: b, x, grid, stencil
  integer, dimension(p_xdim) :: iupper, ilower
  
  contains

  generic :: new => init_fpois2d_mg_explicit,init_fpois2d_mg
  generic :: renew => renew_fpois2d_mg_explicit
  generic :: get_sol => get_hypre_sol
  generic :: set_src => set_hypre_src
  generic :: solve => solve_hypre_poisson
  generic :: del => end_fpois2d_mg
  procedure, private :: init_fpois2d_mg
  procedure, private :: init_fpois2d_mg_explicit
  procedure, private :: renew_fpois2d_mg_explicit
  procedure, private :: end_fpois2d_mg
  ! procedure :: print_matrix
  procedure, private :: get_hypre_sol
  procedure, private :: set_hypre_src
  procedure, private :: solve_hypre_poisson
  ! procedure, private :: read_hypre_sol_ptrcopy

end type fpois2d_mg

type fpois2d_mg_link
  type( fpois2d_mg_link ), pointer :: next => null()
  type( fpois2d_mg ), pointer :: table => null()
  integer :: refcount
end type fpois2d_mg_link

integer, save :: numtables = 0
type( fpois2d_mg_link ), target, save :: table_list

contains

subroutine init_fpois2d_mg( this, pp, perr, noff, nd, ndp, order, dex, solver_type, &
  precond_type, tol )

  implicit none

  class( fpois2d_mg ), intent(inout) :: this
  class( parallel_pipe ), intent(in), pointer :: pp
  class( perrors ), intent(in), pointer :: perr
  integer, intent(in), dimension(2) :: nd, ndp, noff
  integer, intent(in) :: order, solver_type
  integer, intent(in), optional :: precond_type
  ! integer, intent(in), dimension(4) :: gc
  real, intent(in) :: dex, tol

  ! local variables
  integer :: i, j, ierr
  integer, dimension(4), parameter :: p_bnd = (/p_bnd_dd, p_bnd_dn, p_bnd_nd, p_bnd_nn/)
  integer, dimension(2), parameter :: p_stencil = (/p_stencil_poisson, p_stencil_pseudo/)
  integer :: n_pre = 1, n_post = 1, maxiter = 1000
  character(len=18), save :: sname = 'init_fpois2d_mg:'

  call perr%werrfl2( class//sname//' started' )
  this%p => pp
  this%err => perr

  this%phys_bnd_flag = is_phys_bnd( noff, nd, ndp )
  this%order = order
  this%solver_type = solver_type
  this%tol = tol
  if ( present(precond_type) ) then
    this%precond_type = precond_type
  else
    this%precond_type = p_precond_none
  endif

  ! setup HYPRE grid
  call HYPRE_StructGridCreate( this%p%getlgrp(), p_xdim, this%grid, ierr )

  this%ilower = (/ noff(1) + 1, noff(2) + 1 /)
  this%iupper = (/ noff(1) + ndp(1), noff(2) + ndp(2) /)

  call HYPRE_StructGridSetExtents( this%grid, this%ilower, this%iupper, ierr )
  call HYPRE_StructGridAssemble( this%grid, ierr )

  ! set up HYPRE stencil
  call setup_hypre_stencil( this%stencil, this%stencil_num, this%offsets, this%stencil_idx, order )

  if (.not. associated(values)) allocate( values( ndp(1)*ndp(2)*this%stencil_num ) )


  ! set up HYPRE matrix
  do j = 1, 2
    do i = 1, 4
      call HYPRE_StructMatrixCreate( this%p%getlgrp(), this%grid, this%stencil, this%A(i,j), ierr )
      call HYPRE_StructMatrixInitialize( this%A(i,j), ierr )
      ! call HYPRE_StructMatrixSetNumGhost( this%A(i), gc, ierr )
      call setvalue_hypre_matrix_inter( this%A(i,j), this%order, this%ilower, this%iupper, dex, p_stencil(j) )
      call setvalue_hypre_matrix_bnd( this%A(i,j), this%order, this%ilower, this%iupper, p_bnd(i), this%phys_bnd_flag, dex )
      call HYPRE_StructMatrixAssemble( this%A(i,j), ierr )
    enddo
  enddo

  ! set up HYPRE vector
  call HYPRE_StructVectorCreate( this%p%getlgrp(), this%grid, this%b, ierr )
  call HYPRE_StructVectorInitialize( this%b, ierr )

  call HYPRE_StructVectorCreate( this%p%getlgrp(), this%grid, this%x, ierr )
  call HYPRE_StructVectorInitialize( this%x, ierr )

  ! set up solver
  do j = 1, 2
    do i = 1, 4
      call setup_hypre_solver( this%p%getlgrp(), this%solver(i,j), this%A(i,j), this%b, this%x, &
      solver_type, maxiter, n_pre, n_post, tol, this%precond(i,j), this%precond_type )
    enddo
  enddo

  call perr%werrfl2( class//sname//' ended' )

end subroutine init_fpois2d_mg

subroutine init_fpois2d_mg_explicit(this, pp, perr, noff, nd, ndp, order, dex, solver_type, &
  precond_type, tol, pqe, ppsi )

  implicit none

  class( fpois2d_mg ), intent(inout) :: this
  class( parallel_pipe ), intent(in), pointer :: pp
  class( perrors ), intent(in), pointer :: perr
  integer, intent(in), dimension(2) :: nd, ndp, noff
  integer, intent(in) :: order, solver_type
  integer, intent(in), optional :: precond_type
  real, dimension(:,:,:), intent(in), pointer :: pqe
  real, dimension(:,:,:), intent(in), pointer :: ppsi 
  ! real, intent(in) :: dex, tol
  real, intent(in) :: tol
  real, intent(inout) :: dex


  ! local variables
  integer :: i, j, ierr
  integer, dimension(4), parameter :: p_bnd = (/p_bnd_dd, p_bnd_dn, p_bnd_nd, p_bnd_nn/)
  integer, dimension(3), parameter :: p_stencil = (/p_stencil_poisson, p_stencil_pseudo, p_stencil_explicit/)
  integer :: n_pre = 1, n_post = 1, maxiter = 1000
  character(len=18), save :: sname = 'init_fpois2d_mg_explicit:'

  call perr%werrfl2( class//sname//' started' )
  this%p => pp
  this%err => perr

  this%phys_bnd_flag = is_phys_bnd( noff, nd, ndp )
  this%order = order
  this%solver_type = solver_type
  this%tol = tol
  dex = 0.046875
  if ( present(precond_type) ) then
    this%precond_type = precond_type
  else
    this%precond_type = p_precond_none
  endif

  ! setup HYPRE grid
  call HYPRE_StructGridCreate( this%p%getlgrp(), p_xdim, this%grid, ierr )

  this%ilower = (/ noff(1) + 1, noff(2) + 1 /)
  this%iupper = (/ noff(1) + ndp(1), noff(2) + ndp(2) /)

  call HYPRE_StructGridSetExtents( this%grid, this%ilower, this%iupper, ierr )
  call HYPRE_StructGridAssemble( this%grid, ierr )

  ! set up HYPRE stencil
  call setup_hypre_stencil( this%stencil, this%stencil_num, this%offsets, this%stencil_idx, order )

  ! if (.not. associated(values)) allocate( values( ndp(1)*ndp(2)*this%stencil_num ) )
  if (.not. associated(values)) allocate( values( (ndp(1)+1)*(ndp(2)+1)*this%stencil_num ) )
  ! set up HYPRE matrix
  j = 3
    do i = 1, 4
      call HYPRE_StructMatrixCreate( this%p%getlgrp(), this%grid, this%stencil, this%A(i,j), ierr )
      call HYPRE_StructMatrixInitialize( this%A(i,j), ierr )
      ! call HYPRE_StructMatrixSetNumGhost( this%A(i), gc, ierr )
      ! call setvalue_hypre_matrix_inter( this%A(i,j), this%order, this%ilower, this%iupper, dex, p_stencil(j) )
      call setvalue_hypre_matrix_inter_explicit( this%A(i,j), this%order, this%ilower, this%iupper, dex, p_stencil(j), pqe, ppsi )
      call setvalue_hypre_matrix_bnd( this%A(i,j), this%order, this%ilower, this%iupper, p_bnd(i), this%phys_bnd_flag, dex )
      call HYPRE_StructMatrixAssemble( this%A(i,j), ierr )
    enddo
  

  ! set up HYPRE vector
  call HYPRE_StructVectorCreate( this%p%getlgrp(), this%grid, this%b, ierr )
  call HYPRE_StructVectorInitialize( this%b, ierr )

  call HYPRE_StructVectorCreate( this%p%getlgrp(), this%grid, this%x, ierr )
  call HYPRE_StructVectorInitialize( this%x, ierr )

  ! set up solver
  j = 3
    do i = 1, 4
      call setup_hypre_solver( this%p%getlgrp(), this%solver(i,j), this%A(i,j), this%b, this%x, &
      solver_type, maxiter, n_pre, n_post, tol, this%precond(i,j), this%precond_type )
    enddo

  call perr%werrfl2( class//sname//' ended' )

end subroutine init_fpois2d_mg_explicit

! ----------------------------------------------------------------------------------------
subroutine renew_fpois2d_mg_explicit(this, qe, psi, dex )

  implicit none

  class( fpois2d_mg ), intent(inout) :: this
  class( ufield2d ), intent(in), pointer :: qe, psi
  real, intent(in) :: dex
  integer :: i, j, ierr
  ! real :: dex
  real, dimension(:,:,:), pointer :: pqe  
  real, dimension(:,:,:), pointer :: ppsi  
  integer, dimension(4), parameter :: p_bnd = (/p_bnd_dd, p_bnd_dn, p_bnd_nd, p_bnd_nn/)
  integer, dimension(3), parameter :: p_stencil = (/p_stencil_poisson, p_stencil_pseudo, p_stencil_explicit/)
  character(len=18), save :: sname = 'renew_fpois2d_mg_explicit:'

  call this%err%werrfl2(class//sname//' started')

  pqe => qe%getrf()
  ppsi => psi%getrf()
  !write(*,*) "renew_fpois2d_mg_explicit_dex is", dex
! renew HYPRE matrix
  j = 3
  do i = 1, 4
    call setvalue_hypre_matrix_inter_explicit( this%A(i,j), this%order, this%ilower, this%iupper, dex, p_stencil(j), pqe, ppsi)
    call setvalue_hypre_matrix_bnd( this%A(i,j), this%order, this%ilower, this%iupper, p_bnd(i), this%phys_bnd_flag, dex )
    call HYPRE_StructMatrixAssemble( this%A(i,j), ierr )
  enddo
  call this%err%werrfl2(class//sname//' ended')

end subroutine renew_fpois2d_mg_explicit

! ----------------------------------------------------------------------------------------
subroutine end_fpois2d_mg( this )

  implicit none

  class( fpois2d_mg ), intent(inout) :: this

  integer :: i, j, ierr
  character(len=18), save :: sname = 'end_fpois2d_mg:'

  call this%err%werrfl2( class//sname//' started' )

  if ( associated(values) ) deallocate(values)

  call HYPRE_StructGridDestroy( this%grid, ierr )
  call HYPRE_StructStencilDestroy( this%stencil, ierr )
  call HYPRE_StructVectorDestroy( this%b, ierr )
  call HYPRE_StructVectorDestroy( this%x, ierr )

  ! do j = 1, 2
  !   do i = 1, 4
  do j = 1, 3
    do i = 1, 4
    call HYPRE_StructMatrixDestroy( this%A(i,j), ierr )

    select case ( this%solver_type )
    case ( p_solver_smg )
      call HYPRE_StructSMGDestroy( this%solver(i,j), ierr )
    case ( p_solver_pfmg )
      call HYPRE_StructPFMGDestroy( this%solver(i,j), ierr )
    case ( p_solver_pcg )
      call HYPRE_StructPCGDestroy( this%solver(i,j), ierr )
    end select

    select case ( this%precond_type )
    case ( p_precond_smg )
      call HYPRE_StructSMGDestroy( this%precond(i,j), ierr )
    case ( p_precond_pfmg )
      call HYPRE_StructPFMGDestroy( this%precond(i,j), ierr )
    case ( p_precond_ds )
      ! do nothing
    end select

    enddo
  enddo

  call this%err%werrfl2( class//sname//' ended' )

end subroutine end_fpois2d_mg

! ----------------------------------------------------------------------------------------
subroutine setup_hypre_stencil( stencil, stencil_num, offsets, stencil_idx, order )

  implicit none

  ! class( fpois2d_mg ), intent(inout) :: this
  integer(HYPRE_TYPE), intent(inout) :: stencil
  integer, intent(inout) :: stencil_num
  integer, intent(inout), dimension(:,:), pointer :: offsets
  integer, intent(inout), dimension(:), pointer :: stencil_idx
  integer, intent(in) :: order

  ! local
  integer :: i, ierr

  select case ( order )
  case ( p_fd_2order )

    stencil_num = 5
    if ( .not. associated( offsets ) ) then
      allocate( offsets( p_xdim, stencil_num ) )
    endif

    if ( .not. associated( stencil_idx ) ) then
      allocate( stencil_idx( stencil_num ) )
      do i = 1, stencil_num
        stencil_idx(i) = i-1
      enddo
    endif

    offsets(:,1) = (/ 0, 0/)
    offsets(:,2) = (/-1, 0/)
    offsets(:,3) = (/ 1, 0/)
    offsets(:,4) = (/ 0,-1/)
    offsets(:,5) = (/ 0, 1/)

  case ( p_fd_4order )

    print *, "HYPRE solver not implemented for 4order algorithm"
    stop

  case default

    print *, "HYPRE solver not implemented for specified algorithm"
    stop

  end select

  call HYPRE_StructStencilCreate( p_xdim, stencil_num, stencil, ierr )

  do i = 1, stencil_num
    call HYPRE_StructStencilSetElement( stencil, i-1, offsets(1,i), ierr )
  enddo

end subroutine setup_hypre_stencil

! ! ----------------------------------------------------------------------------------------
subroutine setvalue_hypre_matrix_inter( matrix, order, ilower, iupper, dex, kind )

  implicit none

  integer(HYPRE_TYPE), intent(inout) :: matrix
  integer, intent(in) :: order, kind
  integer, intent(in), dimension(2) :: ilower, iupper
  real, intent(in) :: dex

  ! local variables
  integer :: i, nd1p, nd2p, stencil_num, local_vol, ierr
  integer, dimension(:), pointer :: stencil_idx
  real :: idex2

  idex2 = 1.0 / (dex*dex)
  nd1p = iupper(1) - ilower(1) + 1
  nd2p = iupper(2) - ilower(2) + 1
  
  select case ( order )
  case ( p_fd_2order )

    stencil_num = 5
    local_vol = nd1p * nd2p
    allocate( stencil_idx( stencil_num ) )
    stencil_idx = (/0,1,2,3,4/)

    select case ( kind )

    case ( p_stencil_poisson )

      do i = 1, local_vol*stencil_num, stencil_num
        values(i)   =  4.0 * idex2
        values(i+1) = -1.0 * idex2
        values(i+2) = -1.0 * idex2
        values(i+3) = -1.0 * idex2
        values(i+4) = -1.0 * idex2
      enddo

    case ( p_stencil_pseudo )

      do i = 1, local_vol*stencil_num, stencil_num
        values(i)   =  4.0 * idex2 + 1.0
        values(i+1) = -1.0 * idex2
        values(i+2) = -1.0 * idex2
        values(i+3) = -1.0 * idex2
        values(i+4) = -1.0 * idex2
      enddo

    end select

    call HYPRE_StructMatrixSetBoxValues( matrix, ilower, iupper, stencil_num, stencil_idx, values, ierr )
  
  case ( p_fd_4order )

    print *, "HYPRE solver not implemented for 4order algorithm"
    stop

  case default

    print *, "HYPRE solver not implemented for specified algorithm"
    stop

  end select

end subroutine setvalue_hypre_matrix_inter

subroutine setvalue_hypre_matrix_inter_explicit( matrix, order, ilower, iupper, dex, kind, pqe, ppsi )

  implicit none

  integer(HYPRE_TYPE), intent(inout) :: matrix
  integer, intent(in) :: order, kind
  integer, intent(in), dimension(2) :: ilower, iupper
  real, intent(in) :: dex

  ! local variables
  ! integer :: i, nd1p, nd2p, stencil_num, local_vol, ierr
  integer :: i, y, x, nd1p, nd2p, stencil_num, local_vol, ierr
  integer, dimension(:), pointer :: stencil_idx
  real, dimension(:,:,:), pointer :: pqe  
  real, dimension(:,:,:), pointer :: ppsi 
  real :: idex2

  idex2 = 1.0 / (dex*dex)
  nd1p = iupper(1) - ilower(1) + 1
  nd2p = iupper(2) - ilower(2) + 1

  ! nd1p = iupper(1) - ilower(1) + 2
  ! nd2p = iupper(2) - ilower(2) + 2
  !write(*,*) "idex2", idex2
  ! write(*,*) "iupper", iupper, ilower
  ! write(*,*) "nd1p", nd1p
  ! write(*,*) "nd2p", nd2p
  select case ( order )
  case ( p_fd_2order )

    stencil_num = 5
    local_vol = nd1p * nd2p
    allocate( stencil_idx( stencil_num ) )
    stencil_idx = (/0,1,2,3,4/)

    select case ( kind )

    case ( p_stencil_poisson )

      do i = 1, local_vol*stencil_num, stencil_num
        values(i)   =  4.0 * idex2
        values(i+1) = -1.0 * idex2
        values(i+2) = -1.0 * idex2
        values(i+3) = -1.0 * idex2
        values(i+4) = -1.0 * idex2
      enddo

    case ( p_stencil_pseudo )

      do i = 1, local_vol*stencil_num, stencil_num
        values(i)   =  4.0 * idex2 + 1.0
        values(i+1) = -1.0 * idex2
        values(i+2) = -1.0 * idex2
        values(i+3) = -1.0 * idex2
        values(i+4) = -1.0 * idex2
      enddo

    case ( p_stencil_explicit )

      i = 1
      do y = 1, nd2p
        do x = 1, nd1p 
          values(i)   = 4.0 * idex2 - pqe(1, x, y) / (1 + ppsi(1, x, y))
          values(i+1) = -1.0 * idex2
          values(i+2) = -1.0 * idex2
          values(i+3) = -1.0 * idex2
          values(i+4) = -1.0 * idex2
          i = i + 5
        enddo
      enddo
    end select
    !write(2,*) values, "values"  
    call HYPRE_StructMatrixSetBoxValues( matrix, ilower, iupper, stencil_num, stencil_idx, values, ierr )
  
  case ( p_fd_4order )

    print *, "HYPRE solver not implemented for 4order algorithm"
    stop

  case default

    print *, "HYPRE solver not implemented for specified algorithm"
    stop

  end select 

end subroutine setvalue_hypre_matrix_inter_explicit

subroutine setvalue_hypre_matrix_bnd_ori( matrix, order, ilower, iupper, bnd_type, &
  bnd_flag, dex )

  implicit none

  integer(HYPRE_TYPE), intent(inout) :: matrix
  integer, intent(in), dimension(2) :: ilower, iupper
  integer, intent(in) :: order, bnd_type
  logical, intent(in), dimension(4) :: bnd_flag
  real, intent(in) :: dex

  ! local
  integer, dimension(2) :: ilower_bnd, iupper_bnd, icorner
  integer :: nd1p, nd2p, stencil_num, ierr, i
  integer, dimension(:), pointer :: stencil_idx
  real :: idex2

  nd1p = iupper(1) - ilower(1) + 1
  nd2p = iupper(2) - ilower(2) + 1
  
  idex2 = 1.0/(dex**2)

  select case ( order )
  case ( p_fd_2order )

    !stencil_num = 5
    stencil_num = 2
    allocate( stencil_idx( stencil_num ) )

    ! top boundary
    if ( bnd_flag(1) ) then

      ilower_bnd = ilower
      iupper_bnd = (/ ilower(1), iupper(2) /)

      select case ( bnd_type )
      case ( p_bnd_dd, p_bnd_dn )

        stencil_num = 1
        stencil_idx(1) = 1
        do i = 1, nd2p*stencil_num
          values(i) = 0.0
        enddo

      case ( p_bnd_nn, p_bnd_nd )

        stencil_num = 2
        stencil_idx(1:2) = (/1,2/)
        do i = 1, nd2p*stencil_num, stencil_num
          values(i)   = 0.0
          values(i+1) = -2.0 * idex2
        enddo

      end select

      call HYPRE_StructMatrixSetBoxValues( matrix, ilower_bnd, iupper_bnd, stencil_num, &
        stencil_idx, values, ierr )

    endif

    ! bottom boundary
    if ( bnd_flag(2) ) then

      ilower_bnd = (/ iupper(1), ilower(2) /)
      iupper_bnd = iupper

      select case ( bnd_type )
      case ( p_bnd_dd, p_bnd_dn )

        stencil_num = 1
        stencil_idx(1) = 2
        do i = 1, nd2p*stencil_num
          values(i) = 0.0
        enddo

      case ( p_bnd_nn, p_bnd_nd )

        stencil_num = 2
        stencil_idx(1:2) = (/1,2/)
        do i = 1, nd2p*stencil_num, stencil_num
          values(i)   = -2.0 * idex2
          values(i+1) = 0.0
        enddo

      end select

      call HYPRE_StructMatrixSetBoxValues( matrix, ilower_bnd, iupper_bnd, stencil_num, &
        stencil_idx, values, ierr )

    endif

    ! left boundary
    if ( bnd_flag(3) ) then

      ilower_bnd = ilower
      iupper_bnd = (/ iupper(1), ilower(2) /)

      select case ( bnd_type )
      case ( p_bnd_dd, p_bnd_nd )

        stencil_num = 1
        stencil_idx(1) = 3
        do i = 1, nd1p*stencil_num
          values(i) = 0.0
        enddo

      case ( p_bnd_nn, p_bnd_dn )

        stencil_num = 2
        stencil_idx(1:2) = (/3,4/)
        do i = 1, nd1p*stencil_num, stencil_num
          values(i)   = 0.0
          values(i+1) = -2.0 * idex2
        enddo

      end select

      call HYPRE_StructMatrixSetBoxValues( matrix, ilower_bnd, iupper_bnd, stencil_num, &
        stencil_idx, values, ierr )

    endif

    ! right boundary
    if ( bnd_flag(4) ) then

      ilower_bnd = (/ ilower(1), iupper(2) /)
      iupper_bnd = iupper

      select case ( bnd_type )
      case ( p_bnd_dd, p_bnd_nd )
        
        stencil_num = 1
        stencil_idx(1) = 4
        do i = 1, nd1p*stencil_num
          values(i) = 0.0
        enddo

      case ( p_bnd_nn, p_bnd_dn )

        stencil_num = 2
        stencil_idx(1:2) = (/3,4/)
        do i = 1, nd1p*stencil_num, stencil_num
          values(i)   = -2.0 * idex2
          values(i+1) = 0.0
        enddo

      end select

      call HYPRE_StructMatrixSetBoxValues( matrix, ilower_bnd, iupper_bnd, stencil_num, &
        stencil_idx, values, ierr )

    endif

    ! deal four corners if pure Neumann boundary is used
    ! top left
    if ( bnd_type == p_bnd_nn ) then

      if ( bnd_flag(1) .and. bnd_flag(3) ) then 
        icorner = ilower
        stencil_idx = (/2,4/)
      endif
      ! top right
      if ( bnd_flag(1) .and. bnd_flag(4) ) then
        icorner = (/ ilower(1), iupper(2) /)
        stencil_idx = (/2,3/)
      endif
      ! bottom left
      if ( bnd_flag(2) .and. bnd_flag(3) ) then
        icorner = (/ iupper(1), ilower(2) /)
        stencil_idx = (/1,4/)
      endif
      ! bottom right
      if ( bnd_flag(2) .and. bnd_flag(4) ) then
        icorner = iupper
        stencil_idx = (/1,3/)
      endif

      stencil_num = 2
      values(1) = 0.0
      values(2) = 0.0
      call HYPRE_StructMatrixSetValues( matrix, icorner, stencil_num, &
        stencil_idx, values, ierr )

    endif

    deallocate( stencil_idx )

  case ( p_fd_4order )

    print *, "HYPRE solver not implemented for 4order algorithm"
    stop

  case default

    print *, "HYPRE solver not implemented for specified algorithm"
    stop

  end select

end subroutine setvalue_hypre_matrix_bnd_ori

subroutine setvalue_hypre_matrix_bnd_1( matrix, order, ilower, iupper, bnd_type, &
  bnd_flag, dex )

  implicit none

  integer(HYPRE_TYPE), intent(inout) :: matrix
  integer, intent(in), dimension(2) :: ilower, iupper
  integer, intent(in) :: order, bnd_type
  logical, intent(in), dimension(4) :: bnd_flag
  real, intent(in) :: dex

  ! local
  integer, dimension(2) :: ilower_bnd, iupper_bnd, icorner
  integer :: nd1p, nd2p, stencil_num, ierr, i
  integer, dimension(:), pointer :: stencil_idx
  real :: idex2

  nd1p = iupper(1) - ilower(1) + 1
  nd2p = iupper(2) - ilower(2) + 1
  
  idex2 = 1.0/(dex**2)

  select case ( order )
  case ( p_fd_2order )

    stencil_num = 5
    allocate( stencil_idx( stencil_num ) )

    ! top boundary
    if ( bnd_flag(1) ) then

      ilower_bnd = ilower
      iupper_bnd = (/ ilower(1), iupper(2) /)

      select case ( bnd_type )
      case ( p_bnd_dd, p_bnd_dn )

        stencil_num = 4
        stencil_idx(1:4) = (/0,1,3,4/)
        do i = 1, nd2p*stencil_num,  stencil_num
          values(i) = 0.0
          values(i+1) = 0.0
          values(i+2) = 0.0
          values(i+3) = 0.0
        enddo

      case ( p_bnd_nn, p_bnd_nd )

        stencil_num = 2
        stencil_idx(1:2) = (/1,2/)
        do i = 1, nd2p*stencil_num, stencil_num
          values(i)   = 0.0
          values(i+1) = -2.0 * idex2
        enddo

      end select

      call HYPRE_StructMatrixSetBoxValues( matrix, ilower_bnd, iupper_bnd, stencil_num, &
        stencil_idx, values, ierr )

    endif

    ! bottom boundary
    if ( bnd_flag(2) ) then

      ilower_bnd = (/ iupper(1), ilower(2) /)
      iupper_bnd = iupper

      select case ( bnd_type )
      case ( p_bnd_dd, p_bnd_dn )

        stencil_num = 4
        stencil_idx(1:4) =  (/0,2,3,4/)
        do i = 1, nd2p*stencil_num,  stencil_num
          values(i) = 0.0
          values(i+1) = 0.0
          values(i+2) = 0.0
          values(i+3) = 0.0
        enddo

      case ( p_bnd_nn, p_bnd_nd )

        stencil_num = 2
        stencil_idx(1:2) = (/1,2/)
        do i = 1, nd2p*stencil_num, stencil_num
          values(i)   = -2.0 * idex2
          values(i+1) = 0.0
        enddo

      end select

      call HYPRE_StructMatrixSetBoxValues( matrix, ilower_bnd, iupper_bnd, stencil_num, &
        stencil_idx, values, ierr )

    endif

    ! left boundary
    if ( bnd_flag(3) ) then

      ilower_bnd = ilower
      iupper_bnd = (/ iupper(1), ilower(2) /)

      select case ( bnd_type )
      case ( p_bnd_dd, p_bnd_nd )

        stencil_num = 4
        stencil_idx(1:4) = (/0,1,2,3/)
        do i = 1, nd1p*stencil_num,  stencil_num
          values(i) = 0.0
          values(i+1) = 0.0
          values(i+2) = 0.0
          values(i+3) = 0.0
        enddo

      case ( p_bnd_nn, p_bnd_dn )

        stencil_num = 2
        stencil_idx(1:2) = (/3,4/)
        do i = 1, nd1p*stencil_num, stencil_num
          values(i)   = 0.0
          values(i+1) = -2.0 * idex2
        enddo

      end select

      call HYPRE_StructMatrixSetBoxValues( matrix, ilower_bnd, iupper_bnd, stencil_num, &
        stencil_idx, values, ierr )

    endif

    ! right boundary
    if ( bnd_flag(4) ) then

      ilower_bnd = (/ ilower(1), iupper(2) /)
      iupper_bnd = iupper

      select case ( bnd_type )
      case ( p_bnd_dd, p_bnd_nd )
        
        stencil_num = 4
        stencil_idx(1:4) = (/0,1,2,4/)
        do i = 1, nd1p*stencil_num,  stencil_num
          values(i) = 0.0
          values(i+1) = 0.0
          values(i+2) = 0.0
          values(i+3) = 0.0
        enddo

      case ( p_bnd_nn, p_bnd_dn )

        stencil_num = 2
        stencil_idx(1:2) = (/3,4/)
        do i = 1, nd1p*stencil_num, stencil_num
          values(i)   = -2.0 * idex2
          values(i+1) = 0.0
        enddo

      end select

      call HYPRE_StructMatrixSetBoxValues( matrix, ilower_bnd, iupper_bnd, stencil_num, &
        stencil_idx, values, ierr )

    endif

    ! deal four corners if pure Neumann boundary is used
    ! top left
    if ( bnd_type == p_bnd_nn ) then

      if ( bnd_flag(1) .and. bnd_flag(3) ) then 
        icorner = ilower
        stencil_idx(1:2) = (/2,4/)
      endif
      ! top right
      if ( bnd_flag(1) .and. bnd_flag(4) ) then
        icorner = (/ ilower(1), iupper(2) /)
        stencil_idx(1:2) = (/2,3/)
      endif
      ! bottom left
      if ( bnd_flag(2) .and. bnd_flag(3) ) then
        icorner = (/ iupper(1), ilower(2) /)
        stencil_idx(1:2) = (/1,4/)
      endif
      ! bottom right
      if ( bnd_flag(2) .and. bnd_flag(4) ) then
        icorner = iupper
        stencil_idx(1:2) = (/1,3/)
      endif

      stencil_num = 2
      values(1) = 0.0
      values(2) = 0.0
      call HYPRE_StructMatrixSetValues( matrix, icorner, stencil_num, &
        stencil_idx, values, ierr )

    endif

    deallocate( stencil_idx )

  case ( p_fd_4order )

    print *, "HYPRE solver not implemented for 4order algorithm"
    stop

  case default

    print *, "HYPRE solver not implemented for specified algorithm"
    stop

  end select

end subroutine setvalue_hypre_matrix_bnd_1

subroutine setvalue_hypre_matrix_bnd( matrix, order, ilower, iupper, bnd_type, &
  bnd_flag, dex )

  implicit none

  integer(HYPRE_TYPE), intent(inout) :: matrix
  integer, intent(in), dimension(2) :: ilower, iupper
  integer, intent(in) :: order, bnd_type
  logical, intent(in), dimension(4) :: bnd_flag
  real, intent(in) :: dex

  ! local
  integer, dimension(2) :: ilower_bnd, iupper_bnd, icorner
  integer :: nd1p, nd2p, stencil_num, ierr, i
  integer, dimension(:), pointer :: stencil_idx
  real :: idex2

  nd1p = iupper(1) - ilower(1) + 1
  nd2p = iupper(2) - ilower(2) + 1

  ! nd1p = iupper(1) - ilower(1) + 2
  ! nd2p = iupper(2) - ilower(2) + 2
  
  idex2 = 1.0/(dex**2)

  select case ( order )
  case ( p_fd_2order )

    stencil_num = 5
    allocate( stencil_idx( stencil_num ) )

    ! top boundary
    if ( bnd_flag(1) ) then

      ilower_bnd = ilower
      iupper_bnd = (/ ilower(1), iupper(2) /)

      select case ( bnd_type )
      case ( p_bnd_dd, p_bnd_dn )

        ! stencil_num = 5
        ! stencil_idx(1:5) = (/0,1,2,3,4/)
        ! do i = 1, nd2p*stencil_num, stencil_num
        !   values(i) = 1.0 
        !   values(i+1) = 0.0
        !   values(i+2) = 0.0
        !   values(i+3) = 0.0
        !   values(i+4) = 0.0
        ! enddo

        ! stencil_num = 4
        ! stencil_idx(1:4) = (/0,1,3,4/)
        ! do i = 1, nd2p*stencil_num, stencil_num
        !   values(i) = 0.0 
        !   values(i+1) = 0.0
        !   values(i+2) = 0.0
        !   values(i+3) = 0.0
        ! enddo

        stencil_num = 1
        stencil_idx(1) = 1
        do i = 1, nd2p*stencil_num
          values(i) = 0.0
        enddo

      case ( p_bnd_nn, p_bnd_nd )

        stencil_num = 2
        stencil_idx(1:2) = (/1,2/)
        do i = 1, nd2p*stencil_num, stencil_num
          values(i)   = 0.0
          values(i+1) = -2.0 * idex2
        enddo

      end select

      call HYPRE_StructMatrixSetBoxValues( matrix, ilower_bnd, iupper_bnd, stencil_num, &
        stencil_idx, values, ierr )

    endif

    ! bottom boundary
    if ( bnd_flag(2) ) then

      ilower_bnd = (/ iupper(1), ilower(2) /)
      iupper_bnd = iupper

      select case ( bnd_type )
      case ( p_bnd_dd, p_bnd_dn )

        ! stencil_num = 5
        ! stencil_idx(1:5) = (/0,1,2,3,4/)
        ! do i = 1, nd2p*stencil_num, stencil_num
        !   values(i) = 1.0 
        !   values(i+1) = 0.0
        !   values(i+2) = 0.0
        !   values(i+3) = 0.0
        !   values(i+4) = 0.0
        ! enddo

        ! stencil_num = 4
        ! stencil_idx(1:4) = (/0,2,3,4/)
        ! do i = 1, nd2p*stencil_num, stencil_num
        !   values(i) = 0.0 
        !   values(i+1) = 0.0
        !   values(i+2) = 0.0
        !   values(i+3) = 0.0
        ! enddo

        stencil_num = 1
        stencil_idx(1) = 2
        do i = 1, nd2p*stencil_num
          values(i) = 0.0
        enddo

      case ( p_bnd_nn, p_bnd_nd )

        stencil_num = 2
        stencil_idx(1:2) = (/1,2/)
        do i = 1, nd2p*stencil_num, stencil_num
          values(i)   = -2.0 * idex2
          values(i+1) = 0.0
        enddo

      end select

      call HYPRE_StructMatrixSetBoxValues( matrix, ilower_bnd, iupper_bnd, stencil_num, &
        stencil_idx, values, ierr )

    endif

    ! left boundary
    if ( bnd_flag(3) ) then

      ilower_bnd = ilower
      iupper_bnd = (/ iupper(1), ilower(2) /)

      select case ( bnd_type )
      case ( p_bnd_dd, p_bnd_nd )

        ! stencil_num = 5
        ! stencil_idx(1:5) = (/0,1,2,3,4/)
        ! do i = 1, nd1p*stencil_num, stencil_num
        !   values(i) = 1.0
        !   values(i+1) = 0.0
        !   values(i+2) = 0.0
        !   values(i+3) = 0.0
        !   values(i+4) = 0.0
        ! enddo

        ! stencil_num = 1
        ! stencil_idx(1) = 3
        ! do i = 1, nd1p*stencil_num
        !   values(i) = 0.0
        ! enddo

        ! stencil_num = 4
        ! stencil_idx(1:4) = (/0,1,2,3/)
        ! do i = 1, nd1p*stencil_num, stencil_num
        !   values(i) = 0.0 
        !   values(i+1) = 0.0
        !   values(i+2) = 0.0
        !   values(i+3) = 0.0
        ! enddo

        stencil_num = 1
        stencil_idx(1) = 3
        do i = 1, nd1p*stencil_num
          values(i) = 0.0
        enddo

      case ( p_bnd_nn, p_bnd_dn )

        stencil_num = 2
        stencil_idx(1:2) = (/3,4/)
        do i = 1, nd1p*stencil_num, stencil_num
          values(i)   = 0.0
          values(i+1) = -2.0 * idex2
        enddo

      end select

      call HYPRE_StructMatrixSetBoxValues( matrix, ilower_bnd, iupper_bnd, stencil_num, &
        stencil_idx, values, ierr )

    endif

    ! right boundary
    if ( bnd_flag(4) ) then

      ilower_bnd = (/ ilower(1), iupper(2) /)
      iupper_bnd = iupper

      select case ( bnd_type )
      case ( p_bnd_dd, p_bnd_nd )
        
        ! stencil_num = 5
        ! stencil_idx(1:5) = (/0,1,2,3,4/)
        ! do i = 1, nd1p*stencil_num, stencil_num
        !   values(i) = 1.0 
        !   values(i+1) = 0.0
        !   values(i+2) = 0.0
        !   values(i+3) = 0.0
        !   values(i+4) = 0.0
        ! enddo

        ! stencil_num = 1
        ! stencil_idx(1) = 4
        ! do i = 1, nd1p*stencil_num
        !   values(i) = 0.0
        ! enddo

        ! stencil_num = 4
        ! stencil_idx(1:4) = (/0,1,2,4/)
        ! do i = 1, nd1p*stencil_num, stencil_num
        !   values(i) = 0.0 
        !   values(i+1) = 0.0
        !   values(i+2) = 0.0
        !   values(i+3) = 0.0
        ! enddo

        stencil_num = 1
        stencil_idx(1) = 4
        do i = 1, nd1p*stencil_num
          values(i) = 0.0
        enddo

      case ( p_bnd_nn, p_bnd_dn )

        stencil_num = 2
        stencil_idx(1:2) = (/3,4/)
        do i = 1, nd1p*stencil_num, stencil_num
          values(i)   = -2.0 * idex2
          values(i+1) = 0.0
        enddo

      end select

      call HYPRE_StructMatrixSetBoxValues( matrix, ilower_bnd, iupper_bnd, stencil_num, &
        stencil_idx, values, ierr )

    endif

    ! deal four corners if pure Neumann boundary is used
    ! top left
    if ( bnd_type == p_bnd_nn ) then

      if ( bnd_flag(1) .and. bnd_flag(3) ) then 
        icorner = ilower
        stencil_idx(1:2) = (/2,4/)
      endif
      ! top right
      if ( bnd_flag(1) .and. bnd_flag(4) ) then
        icorner = (/ ilower(1), iupper(2) /)
        stencil_idx(1:2) = (/2,3/)
      endif
      ! bottom left
      if ( bnd_flag(2) .and. bnd_flag(3) ) then
        icorner = (/ iupper(1), ilower(2) /)
        stencil_idx(1:2) = (/1,4/)
      endif
      ! bottom right
      if ( bnd_flag(2) .and. bnd_flag(4) ) then
        icorner = iupper
        stencil_idx(1:2) = (/1,3/)
      endif

      stencil_num = 2
      values(1) = 0.0
      values(2) = 0.0
      call HYPRE_StructMatrixSetValues( matrix, icorner, stencil_num, &
        stencil_idx, values, ierr )

    endif

    deallocate( stencil_idx )

  case ( p_fd_4order )

    print *, "HYPRE solver not implemented for 4order algorithm"
    stop

  case default

    print *, "HYPRE solver not implemented for specified algorithm"
    stop

  end select

end subroutine setvalue_hypre_matrix_bnd

subroutine setup_hypre_solver( comm, solver, A, b, x, type_, maxiter, n_pre, n_post, tol, &
  precond, precond_type )

  implicit none

  integer, intent(in) :: comm
  integer(HYPRE_TYPE), intent(inout) :: solver, A, b, x
  integer, intent(in) :: type_, maxiter, n_pre, n_post
  real, intent(in) :: tol
  integer(HYPRE_TYPE), intent(inout), optional :: precond
  integer, intent(in), optional :: precond_type

  ! local
  integer :: ierr

  select case ( type_ )

  case ( p_solver_smg )

    print *, 'Using SMG solver'
    call HYPRE_StructSMGCreate( comm, solver, ierr )
    call HYPRE_StructSMGSetMemoryUse( solver, 0, ierr )
    call HYPRE_StructSMGSetMaxIter( solver, maxiter, ierr )
    call HYPRE_StructSMGSetTol( solver, tol, ierr )
    ! call HYPRE_StructSMGSetPrintLevel( solver, 2, ierr )
    call HYPRE_StructSMGSetRelChange( solver, 0, ierr )
    call HYPRE_StructSMGSetNumPreRelax( solver, n_pre, ierr )
    call HYPRE_StructSMGSetNumPostRelax( solver, n_post, ierr )
    call HYPRE_StructSMGSetLogging( solver, 1, ierr )
    call HYPRE_StructSMGSetup( solver, A, b, x, ierr )

  case ( p_solver_pfmg )

    print *, 'Using PFMG solver'
    call HYPRE_StructPFMGCreate( comm, solver, ierr )
    call HYPRE_StructPFMGSetMaxIter( solver, maxiter, ierr )
    call HYPRE_StructPFMGSetTol( solver, tol, ierr )
    ! call HYPRE_StructPFMGSetPrintLevel( solver, 2, ierr )
    call HYPRE_StructPFMGSetRelChange( solver, 0, ierr )
    ! weighted Jacobi = 1; red-black GS = 2
    call HYPRE_StructPFMGSetRelaxType( solver, 1, ierr )
    call HYPRE_StructPFMGSetNumPreRelax( solver, n_pre, ierr )
    call HYPRE_StructPFMGSetNumPostRelax( solver, n_post, ierr )
    call HYPRE_StructPFMGSetLogging( solver, 1, ierr )
    call HYPRE_StructPFMGSetup( solver, A, b, x, ierr )

  case ( p_solver_pcg )

    print *, 'Using PCG solver'
    call HYPRE_StructPCGCreate( comm, solver, ierr )
    call HYPRE_StructPCGSetMaxIter( solver, maxiter, ierr )
    call HYPRE_StructPCGSetTol( solver, tol, ierr )
    ! call HYPRE_StructPCGSetPrintLevel( solver, 2, ierr )
    call HYPRE_StructPCGSetTwoNorm( solver, 1, ierr )
    call HYPRE_StructPCGSetRelChange( solver, 0, ierr )
    call HYPRE_StructPCGSetLogging( solver, 1, ierr )

    if ( present(precond) .and. present(precond_type) ) then
      call setup_hypre_precond( comm, precond, solver, precond_type, n_pre, n_post )
    endif

    call HYPRE_StructPCGSetup( solver, A, b, x, ierr )

  end select

end subroutine setup_hypre_solver

subroutine setup_hypre_precond( comm, precond, solver, type_, n_pre, n_post )

  implicit none

  integer, intent(in) :: comm
  integer(HYPRE_TYPE), intent(inout) :: precond, solver
  integer, intent(in) :: type_, n_pre, n_post

  ! local
  integer :: ierr, precond_id, maxiter = 1
  real :: tol = 0.0

  select case ( type_ )

  case ( p_precond_smg )

    precond_id = 0
    call HYPRE_StructSMGCreate( comm, precond, ierr )
    call HYPRE_StructSMGSetMemoryUse( precond, 0, ierr )
    call HYPRE_StructSMGSetMaxIter( precond, maxiter, ierr )
    call HYPRE_StructSMGSetTol( precond, tol, ierr )
    call HYPRE_StructSMGSetNumPreRelax( precond, n_pre, ierr )
    call HYPRE_StructSMGSetNumPostRelax( precond, n_post, ierr )
    call HYPRE_StructSMGSetLogging( precond, 0, ierr )

  case ( p_precond_pfmg )

    precond_id = 1
    call HYPRE_StructPFMGCreate( comm, precond, ierr )
    call HYPRE_StructPFMGSetMaxIter( precond, maxiter, ierr )
    call HYPRE_StructPFMGSetTol( precond, tol, ierr )
    ! 1 - Jacobi, 2 - red-black GS
    call HYPRE_StructPFMGSetRelaxType( precond, 1, ierr )
    call HYPRE_StructPFMGSetNumPreRelax( precond, n_pre, ierr )
    call HYPRE_StructPFMGSetNumPostRelax( precond, n_post, ierr )
    call HYPRE_StructPFMGSetLogging( precond, 0, ierr )
    
  case ( p_precond_ds )

    precond_id = 8
    precond = 0

  case ( p_precond_none )

    return

  end select

  call HYPRE_StructPCGSetPrecond( solver, precond_id, precond, ierr )

end subroutine setup_hypre_precond

subroutine get_hypre_sol( this, dspace, dim )
! this is the hard copy version

  implicit none

  class( fpois2d_mg ), intent(inout) :: this
  class( ufield2d ), intent(inout), pointer :: dspace
  integer, intent(in) :: dim
  ! integer, intent(in), dimension(2) :: nd

  ! local variables
  integer :: i, j, k, ierr, nd1p, nd2p
  real, dimension(:,:,:), pointer :: df => null()
  character(len=18), save :: sname = 'get_hypre_sol:'

  call this%err%werrfl2( class//sname//' started' )

  df => dspace%getrf()
  nd1p = dspace%getnd1p()
  nd2p = dspace%getnd2p()

  ! nd1p = dspace%getnd1p() + 1
  ! nd2p = dspace%getnd2p() + 1

  call HYPRE_StructVectorGetBoxValues( this%x, this%ilower, this%iupper, values, ierr )
  
  k = 1
  do j = 1, nd2p
    do i = 1, nd1p
      df(dim,i,j) = values(k)
      k = k + 1
    enddo
  enddo

  call this%err%werrfl2( class//sname//' ended' )

end subroutine get_hypre_sol

! subroutine read_hypre_sol_ptrcopy( this, ptr_out, n1, n2 )
! ! this is the pointer copy version

!   implicit none

!   class( fpois2d_mg ), intent(inout) :: this
!   real, intent(inout), dimension(:,:), pointer :: ptr_out
!   integer, intent(in) :: n1, n2

!   ! local variables
!   real, dimension(:), pointer :: values_

!   allocate( values_(n1*n2) )
!   deallocate( ptr_out )

!   call HYPRE_StructVectorGetBoxValues( this%x, this%ilower, this%iupper, values_ )

!   ptr_out => reshape_array( values_, (/n1, n2/) )

! end subroutine read_hypre_sol_ptrcopy

subroutine set_hypre_src( this, rspace, dim )

  implicit none

  class( fpois2d_mg ), intent(inout) :: this
  class( ufield2d ), intent(inout), pointer :: rspace
  integer, intent(in) :: dim

  ! local variables
  integer :: i, j, k, ierr, nd1p, nd2p
  real, dimension(:,:,:), pointer :: f => null()
  character(len=18), save :: sname = 'set_hypre_src:'

  call this%err%werrfl2( class//sname//' started' )

  f => rspace%getrf()
  nd1p = rspace%getnd1p()
  nd2p = rspace%getnd2p()

  ! nd1p = rspace%getnd1p() + 1
  ! nd2p = rspace%getnd2p() + 1
  ! write(*,*) "nd1p=", nd1p
  ! write(*,*) "nd2p=", nd2p
  ! f(1,1,:) = 0
  ! f(1,257,:) = 0
  ! f(1,:,1) = 0
  ! f(1,:,257) = 0
  ! f(2,1,:) = 0
  ! f(2,257,:) = 0
  ! f(2,:,1) = 0
  ! f(2,:,257) = 0
  k = 1
  do j = 1, nd2p
    do i = 1, nd1p
      values(k) = f(dim,i,j)
      k = k + 1
    enddo
  enddo
  ! write(2,*) "f(1,1,:)=", f(1,1,:) 
  ! write(2,*) "f(1,257,:)=", f(1,257,:)
  ! write(2,*) "f(1,:,1)=", f(1,:,1)  
  ! write(2,*) "f(1,:,257)=", f(1,:,257)  
  ! write(2,*) "f(2,1,:)=", f(2,1,:) 
  ! write(2,*) "f(2,257,:)=", f(2,257,:)  
  ! write(2,*) "f(2,:,1)=", f(2,:,1)
  ! write(2,*) "f(2,:,257)=", f(2,:,257) 
  call HYPRE_StructVectorSetBoxValues( this%b, this%ilower, this%iupper, values, ierr )

  call this%err%werrfl2( class//sname//' ended' )

end subroutine set_hypre_src

subroutine solve_hypre_poisson( this, src, sol, bnd_types, stencil_kind )

  implicit none

  class( fpois2d_mg ), intent(inout) :: this
  class( ufield2d ), intent(inout), pointer :: src, sol
  integer, intent(in), dimension(:) :: bnd_types
  integer, intent(in) :: stencil_kind

  ! local
  integer :: k, l, dim, ndim, ierr
  character(len=18), save :: sname = 'solve_hypre_poisson:'

  call this%err%werrfl2( class//sname//' started' )
  !write(*,*) "bnd_types",bnd_types
  ndim = size( bnd_types )
  !write(*,*) "ndim",ndim, stencil_kind
  l = stencil_kind
  !write(*,*) "solve_hypre_poisson1"
  do dim = 1, ndim

    k = bnd_types(dim) 
    call this%set_src( src, dim )
    !write(*,*) "solve_hypre_poisson2"
    select case ( this%solver_type )
    case ( p_solver_smg )
      call HYPRE_StructSMGSolve( this%solver(k,l), this%A(k,l), this%b, this%x, ierr )
    case ( p_solver_pfmg )
      call HYPRE_StructPFMGSolve( this%solver(k,l), this%A(k,l), this%b, this%x, ierr )
    case ( p_solver_pcg )
      call HYPRE_StructPCGSolve( this%solver(k,l), this%A(k,l), this%b, this%x, ierr )
    end select
    
    call this%get_sol( sol, dim )

  enddo

  call this%err%werrfl2( class//sname//' ended' )

end subroutine solve_hypre_poisson

function is_phys_bnd( noff, nd, ndp )

  implicit none

  integer, intent(in), dimension(2) :: noff, nd, ndp
  logical, dimension(4) :: is_phys_bnd

  is_phys_bnd(1) = ( noff(1) == 0 )
  is_phys_bnd(2) = ( noff(1)+ndp(1) == nd(1) )
  is_phys_bnd(3) = ( noff(2) == 0 )
  is_phys_bnd(4) = ( noff(2)+ndp(2) == nd(2) )

end function is_phys_bnd

function reshape_array( array, shape_ ) result(aptr)

  use iso_c_binding, only : C_F_POINTER, C_LOC

  implicit none

  real, intent(in), dimension(:), pointer :: array
  integer, intent(in), dimension(2) :: shape_
  real, dimension(:,:), pointer :: aptr

  call C_F_POINTER( C_LOC(array), aptr, shape_ )

end function reshape_array

! subroutine print_matrix( this, which )

!   implicit none

!   class( fpois2d_mg ), intent(in) :: this
!   integer, intent(in) :: which

!   ! local 
!   integer :: ierr

!   select case (which)
!   case (p_bnd_dd)
!     call HYPRE_StructMatrixPrint( this%A(1), 0, ierr )
!   case (p_bnd_dn)
!     call HYPRE_StructMatrixPrint( this%A(2), 0, ierr )
!   case (p_bnd_nd)
!     call HYPRE_StructMatrixPrint( this%A(3), 0, ierr )
!   case (p_bnd_nn)
!     call HYPRE_StructMatrixPrint( this%A(4), 0, ierr )
!   end select

! end subroutine print_matrix

! subroutine print_x( this )

!   implicit none

!   class( fpois2d_mg ), intent(in) :: this

!   ! local
!   integer :: ierr

!   call HYPRE_StructVectorPrint( this%x, 0, ierr )

! end subroutine print_x

function get_poismg2table( pp, perr, psp, order, uqe, upsi ) result(table)
! this function gets an fft table entry, either creates one or points
! to one if the required table already exists for real to complex ffts
  implicit none
  class(perrors), intent(in), pointer :: perr
  class(parallel_pipe), intent(in), pointer :: pp
  class(spect2d), intent(in), pointer :: psp
  class(ufield2d), intent(in), pointer :: uqe, upsi
  integer, intent(in) :: order
  type(fpois2d_mg), pointer :: table

  ! local data
  type(fpois2d_mg_link), pointer :: link => null()
  type(fpois2d_mg), pointer :: ltable => null()
  real, dimension(:,:,:), pointer :: pqe  => null()
  real, dimension(:,:,:), pointer :: ppsi => null()

  integer :: solver_type, precond_type
  
  integer, dimension(2) :: nd, ndp, noff, nvp
  real :: dex, mg_prec
  character(len=18), save :: sname = 'get_poismg2table:'

  call perr%werrfl2(class//sname//' started')

  nvp(1) = 1
  nvp(2) = pp%getlnvp()
  nd(1) = 2**psp%getindx()
  nd(2) = 2**psp%getindy()
  ndp(1) = nd(1)/nvp(1)
  ndp(2) = nd(2)/nvp(2)
  noff(1) = 0
  noff(2) = ( pp%getlkstrt()-1 ) * ndp(2)
  ! write(*,*) "nvp2", nvp(2) 2
  !write(2,*) dex, "dexa" 0
  ! write(*,*) "noff(2)", noff(2)
  dex = psp%getdex()
  !write(2,*) dex, "dexb" 0
  ! dex = 0.0468
  mg_prec = psp%getmg_prec()
  solver_type = psp%getmg_type()
  precond_type = 1 ! SMG preconditioner
  pqe => uqe%getrf()
  ppsi => upsi%getrf()
  nullify(table)         

  if ( numtables == 0 ) then
    nullify( table_list%next, table_list%table )         
    table_list%refcount = 0
  endif
  link => table_list
  table => link%table
  ! search link list of table to see if required table already exists
  do while ( associated(table) )
  ! found it
    if ( order == table%order &
      .and. solver_type == table%solver_type &
      .and. precond_type == table%precond_type &
      .and. mg_prec == table%tol ) then

      link%refcount = link%refcount + 1
      call perr%werrfl2(class//sname//' ended')
      return
  ! check next table, create new empty table if end is reached
    else
      if (associated(link%next)) then
        link => link%next
      else
        allocate(link%next)
        link => link%next
        nullify(link%next,link%table)
        link%refcount = 0
      endif
      table => link%table
    endif
  end do
  ! allocate table entries
  allocate(ltable)
  link%table => ltable
  table => link%table
  ! write(2,*) dex, "dex"
  ! dex = 0.0468
  ! call table%new_fd(pp, perr, noff, nd, ndp, order, dex, solver_type, &
    ! precond_type, mg_prec, pqe, ppsi)
  call table%new(pp, perr, noff, nd, ndp, order, dex, solver_type, &
    precond_type, mg_prec, pqe, ppsi)
  link%refcount = 1
  numtables = numtables + 1

  call perr%werrfl2(class//sname//' ended')

end function get_poismg2table

end module fpois2d_mg_class