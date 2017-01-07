! fdist3d class for QuickPIC Open Source 1.0
! update: 04/18/2016

      module fdist3d_class

      use perrors_class
      use parallel_pipe_class
      use spect3d_class
         
      implicit none

      private

      public :: fdist3d

      type fdist3d

         private

         class(spect3d), pointer, public :: sp => null()
         class(perrors), pointer, public :: err => null()
         class(parallel_pipe), pointer, public :: p => null()
!
! ndprof = profile type
! npx/npy/npz = number of particles distributed in x/y/z direction
! arg = arguments for the profile in x/y/z directions
!        
         integer :: npf
         integer :: npx, npy, npz
         real, dimension(3,100) :: arg
                          
         contains
         
         procedure, private :: init_fdist3d
         procedure, private :: end_fdist3d
         generic :: new => init_fdist3d         
         generic :: del => end_fdist3d
         procedure :: getnpf,getnpx,getnpy,getnpz,getarg
                  
      end type 

      character(len=10), save :: class = 'fdist3d:'
      character(len=128), save :: erstr
      
      contains
!
      subroutine init_fdist3d(this,pp,perr,psp,npf,npx,npy,npz,arg)
      
         implicit none
         
         class(fdist3d), intent(inout) :: this
         class(spect3d), intent(in), pointer :: psp
         class(perrors), intent(in), pointer :: perr
         class(parallel_pipe), intent(in), pointer :: pp
         integer, intent(in) :: npf
         integer, intent(in) :: npx, npy, npz
         real, dimension(3,100), intent(in) :: arg
        
! local data
         character(len=18), save :: sname = 'init_fdist3d:'
         
         this%sp => psp
         this%err => perr
         this%p => pp

         call this%err%werrfl2(class//sname//' started')
         this%npf = npf
         this%npx = npx
         this%npy = npy
         this%npz = npz
         this%arg = arg

         call this%err%werrfl2(class//sname//' ended')

      end subroutine init_fdist3d
!
      subroutine end_fdist3d(this)
          
         implicit none
         
         class(fdist3d), intent(inout) :: this
         character(len=18), save :: sname = 'end_fdist3d:'

         call this%err%werrfl2(class//sname//' started')
         call this%err%werrfl2(class//sname//' ended')
                  
      end subroutine end_fdist3d
!
      function getnpf(this)

         implicit none

         class(fdist3d), intent(in) :: this
         integer :: getnpf
         
         getnpf = this%npf

      end function getnpf  
!      
      function getnpx(this)

         implicit none

         class(fdist3d), intent(in) :: this
         integer :: getnpx
         
         getnpx = this%npx

      end function getnpx
!      
      function getnpy(this)

         implicit none

         class(fdist3d), intent(in) :: this
         integer :: getnpy
         
         getnpy = this%npy

      end function getnpy
!      
      function getnpz(this)

         implicit none

         class(fdist3d), intent(in) :: this
         integer :: getnpz
         
         getnpz = this%npz

      end function getnpz
!      
      function getarg(this)

         implicit none

         class(fdist3d), intent(in) :: this
         real, dimension(3,100) :: getarg
         
         getarg = this%arg

      end function getarg
!      
      end module fdist3d_class