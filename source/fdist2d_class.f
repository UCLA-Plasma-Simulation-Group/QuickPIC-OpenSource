! fdist2d class for QuickPIC Open Source 1.0
! update: 04/18/2016

      module fdist2d_class

      use perrors_class
      use parallel_pipe_class
      use spect2d_class
         
      implicit none

      private

      public :: fdist2d

      type fdist2d

         private

         class(spect2d), pointer, public :: sp => null()
         class(perrors), pointer, public :: err => null()
         class(parallel_pipe), pointer, public :: p => null()
!
! ndprof = profile type 
! npx/npy = number of particles distributed in x/y direction
! arg = arguments for the profile in x/y directions
!        
         integer :: npf
         integer :: npx, npy
         real, dimension(2,100) :: arg
                          
         contains
         
         procedure, private :: init_fdist2d
         procedure, private :: end_fdist2d
         generic :: new => init_fdist2d         
         generic :: del => end_fdist2d
         procedure :: getnpf,getnpx,getnpy,getarg
                  
      end type 

      character(len=10), save :: class = 'fdist2d:'
      character(len=128), save :: erstr
      
      contains
!
      subroutine init_fdist2d(this,pp,perr,psp,npf,npx,npy,arg)
      
         implicit none
         
         class(fdist2d), intent(inout) :: this
         class(spect2d), intent(in), pointer :: psp
         class(perrors), intent(in), pointer :: perr
         class(parallel_pipe), intent(in), pointer :: pp
         integer, intent(in) :: npf
         integer, intent(in) :: npx, npy
         real, dimension(2,100), intent(in) :: arg
        
! local data
         character(len=18), save :: sname = 'init_fdist2d:'
         
         this%sp => psp
         this%err => perr
         this%p => pp

         call this%err%werrfl2(class//sname//' started')
         this%npf = npf
         this%npx = npx
         this%npy = npy
         this%arg = arg

         call this%err%werrfl2(class//sname//' ended')

      end subroutine init_fdist2d
!
      subroutine end_fdist2d(this)
          
         implicit none
         
         class(fdist2d), intent(inout) :: this
         character(len=18), save :: sname = 'end_fdist2d:'

         call this%err%werrfl2(class//sname//' started')
         call this%err%werrfl2(class//sname//' ended')
                  
      end subroutine end_fdist2d
!
      function getnpf(this)

         implicit none

         class(fdist2d), intent(in) :: this
         integer :: getnpf
         
         getnpf = this%npf

      end function getnpf  
!      
      function getnpx(this)

         implicit none

         class(fdist2d), intent(in) :: this
         integer :: getnpx
         
         getnpx = this%npx

      end function getnpx
!      
      function getnpy(this)

         implicit none

         class(fdist2d), intent(in) :: this
         integer :: getnpy
         
         getnpy = this%npy

      end function getnpy
!      
      function getarg(this)

         implicit none

         class(fdist2d), intent(in) :: this
         real, dimension(2,100) :: getarg
         
         getarg = this%arg

      end function getarg
!      
      end module fdist2d_class