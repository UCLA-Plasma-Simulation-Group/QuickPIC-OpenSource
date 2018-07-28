! fdist2d class for QuickPIC Open Source 1.0
! update: 04/18/2016

      module fdist2d_class

      use perrors_class
      use parallel_pipe_class
      use spect2d_class
      use ufield2d_class
      use input_class

      implicit none

      private

      public :: fdist2d, fdist2d_000

      type, abstract :: fdist2d

         private

         class(spect2d), pointer, public :: sp => null()
         class(perrors), pointer, public :: err => null()
         class(parallel_pipe), pointer, public :: p => null()
!
! ndprof = profile type 
         integer :: npf
                          
         contains
         generic :: new => init_fdist2d         
         generic :: del => end_fdist2d
         generic :: dist => dist2d
         procedure(initialize), deferred, private :: init_fdist2d
         procedure, private :: end_fdist2d
         procedure(init_prof), deferred, private :: dist2d
         procedure :: getnpf
                  
      end type fdist2d

      abstract interface
!
      subroutine init_prof(this,part2d,npp,fd,s)
         import fdist2d
         import ufield2d
         implicit none
         class(fdist2d), intent(inout) :: this
         real, dimension(:,:), pointer, intent(inout) :: part2d
         integer, intent(inout) :: npp
         class(ufield2d), intent(in), pointer :: fd         
         real, intent(in) :: s
      end subroutine init_prof
!
      subroutine initialize(this,input,i)
         import fdist2d
         import input_json
         implicit none
         class(fdist2d), intent(inout) :: this
         type(input_json), intent(inout), pointer :: input
         integer, intent(in) :: i
      end subroutine initialize
!
      end interface

      type, extends(fdist2d) :: fdist2d_000

         private

! xppc, yppc = particle per cell in x and y directions
         integer :: xppc, yppc
         real :: qm
                          
         contains
         procedure, private :: init_fdist2d => init_fdist2d_000
         procedure, private :: dist2d => dist2d_000
                  
      end type fdist2d_000


      character(len=10), save :: class = 'fdist2d:'
      character(len=128), save :: erstr
      
      contains
!
      function getnpf(this)

         implicit none

         class(fdist2d), intent(in) :: this
         integer :: getnpf
         
         getnpf = this%npf
      
      end function getnpf  
!      
      subroutine end_fdist2d(this)
          
         implicit none
         
         class(fdist2d), intent(inout) :: this
         character(len=18), save :: sname = 'end_fdist2d:'

         call this%err%werrfl2(class//sname//' started')
         call this%err%werrfl2(class//sname//' ended')
                  
      end subroutine end_fdist2d
!      
      subroutine init_fdist2d_000(this,input,i)
      
         implicit none
         
         class(fdist2d_000), intent(inout) :: this
         type(input_json), intent(inout), pointer :: input
         integer, intent(in) :: i

! local data
         integer :: npf,xppc,yppc
         real :: qm
         character(len=20) :: sn,s1
         character(len=18), save :: sname = 'init_fdist2d_000:'
         
         this%sp => input%sp
         this%err => input%err
         this%p => input%pp

         call this%err%werrfl2(class//sname//' started')
         write (sn,'(I3.3)') i
         s1 = 'species('//trim(sn)//')'
         call input%get(trim(s1)//'.profile',npf)
         call input%get(trim(s1)//'.ppc(1)',xppc)
         call input%get(trim(s1)//'.ppc(2)',yppc)
         call input%get(trim(s1)//'.q',qm)

         this%npf = npf
         this%xppc = xppc
         this%yppc = yppc
         this%qm = qm
         call this%err%werrfl2(class//sname//' ended')

      end subroutine init_fdist2d_000
!
      subroutine dist2d_000(this,part2d,npp,fd,s)
         implicit none
         class(fdist2d_000), intent(inout) :: this
         class(ufield2d), intent(in), pointer :: fd
         real, dimension(:,:), pointer, intent(inout) :: part2d
         real, intent(in) :: s 
         integer, intent(inout) :: npp
! local data
         character(len=18), save :: sname = 'dist2d_000:'
         real, dimension(:,:), pointer :: pt => null()
         integer :: nps, nx, ny, noff, xppc, yppc, i, j
         integer :: ix, iy
         real :: qm

         call this%err%werrfl2(class//sname//' started')
         
         nx = fd%getnd1p(); ny = fd%getnd2p(); noff = fd%getnoff()
         xppc = this%xppc; yppc = this%yppc
         qm = this%qm/abs(this%qm)/xppc/yppc
         nps = 1
         pt => part2d
! initialize the particle positions
         if (noff < ny) then
         do i=2, nx-1
            do j=2, ny
               do ix = 0, xppc-1
                  do iy=0, yppc-1
                     pt(1,nps) = (ix + 0.5)/xppc + i - 1
                     pt(2,nps) = (iy + 0.5)/yppc + j - 1 + noff
                     pt(3,nps) = 0.0
                     pt(4,nps) = 0.0
                     pt(5,nps) = 0.0
                     pt(6,nps) = 1.0
                     pt(7,nps) = 1.0
                     pt(8,nps) = qm
                     nps = nps + 1
                  enddo
               enddo
            enddo
         enddo
         else if (noff > (nx-ny-1)) then       
         do i=2, nx-1
            do j=1, ny-1
               do ix = 0, xppc-1
                  do iy=0, yppc-1
                     pt(1,nps) = (ix + 0.5)/xppc + i - 1
                     pt(2,nps) = (iy + 0.5)/yppc + j - 1 + noff
                     pt(3,nps) = 0.0
                     pt(4,nps) = 0.0
                     pt(5,nps) = 0.0
                     pt(6,nps) = 1.0
                     pt(7,nps) = 1.0
                     pt(8,nps) = qm
                     nps = nps + 1
                  enddo
               enddo
            enddo
         enddo
         else
         do i=2, nx-1
            do j=1, ny
               do ix = 0, xppc-1
                  do iy=0, yppc-1
                     pt(1,nps) = (ix + 0.5)/xppc + i - 1
                     pt(2,nps) = (iy + 0.5)/yppc + j - 1 + noff
                     pt(3,nps) = 0.0
                     pt(4,nps) = 0.0
                     pt(5,nps) = 0.0
                     pt(6,nps) = 1.0
                     pt(7,nps) = 1.0
                     pt(8,nps) = qm
                     nps = nps + 1
                  enddo
               enddo
            enddo
         enddo
         endif
         
         npp = nps - 1
         
         call this%err%werrfl2(class//sname//' ended')

      end subroutine dist2d_000
!
      end module fdist2d_class