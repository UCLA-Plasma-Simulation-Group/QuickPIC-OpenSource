! Spect2d class for QuickPIC Open Source 1.0
! update: 04/18/2016

      module spect2d_class

      use perrors_class
      use parallel_pipe_class
         
      implicit none

      private

      public :: spect2d

      type spect2d

         private

! psolver = solver type = (1) = (conductive)
! inorder = (1) = (linear)         
      !    integer :: indx, indy, psolver, inorder
         integer :: indx, indy, psolver, inorder, mg_type
         real :: dex, mg_prec
         class(perrors), pointer, public :: err => null()
         class(parallel_pipe), pointer, public :: p => null()
         
         contains
         
         procedure, private :: init_spect2d
         procedure, private :: init_spect2d_fd
         procedure, private :: end_spect2d
      !    generic :: new => init_spect2
         generic :: new => init_spect2d, init_spect2d_fd
         generic :: del => end_spect2d
         procedure :: getindx
         procedure :: getindy
         procedure :: getpsolver
         procedure :: getinorder
         procedure :: getdex
         procedure :: getmg_prec
         procedure :: getmg_type

      end type spect2d
      
      contains
!
      subroutine init_spect2d(this,pp,perr,indx,indy,psolver,inorder)
      
         implicit none
         
         class(spect2d), intent(inout) :: this
         class(perrors), intent(in), pointer :: perr
         class(parallel_pipe), intent(in), pointer :: pp
         integer, intent(in) :: indx, indy, psolver, inorder
         
         this%indx = indx
         this%indy = indy
         this%psolver = psolver
         this%inorder = inorder
         this%err => perr
         this%p => pp         
      end subroutine init_spect2d
!
      subroutine init_spect2d_fd(this,pp,perr,indx,indy,dex,psolver,inorder,mg_prec,mg_type)
      
         implicit none
         
         class(spect2d), intent(inout) :: this
         class(perrors), intent(in), pointer :: perr
         class(parallel_pipe), intent(in), pointer :: pp
         integer, intent(in) :: indx, indy, psolver, inorder, mg_type
         real, intent(in) :: dex, mg_prec
         
         this%indx = indx
         this%indy = indy
         this%dex = dex
         this%psolver = psolver
         this%inorder = inorder
         this%err => perr
         this%p => pp         
         this%mg_prec = mg_prec
         this%mg_type = mg_type
      end subroutine init_spect2d_fd
!
      subroutine end_spect2d(this)
          
         implicit none
         
         class(spect2d), intent(inout) :: this
         
         return
         
      end subroutine end_spect2d
!      
      function getindx(this)

         implicit none

         class(spect2d), intent(in) :: this
         integer :: getindx
         
         getindx = this%indx

      end function getindx         
!      
      function getindy(this)

         implicit none

         class(spect2d), intent(in) :: this
         integer :: getindy
         
         getindy = this%indy

      end function getindy         
!      
      function getpsolver(this)

         implicit none

         class(spect2d), intent(in) :: this
         integer :: getpsolver
         
         getpsolver = this%psolver

      end function getpsolver         
!      
      function getinorder(this)

         implicit none

         class(spect2d), intent(in) :: this
         integer :: getinorder
         
         getinorder = this%inorder

      end function getinorder         
!     
      function getdex(this)
         
         implicit none

         class(spect2d), intent(in) :: this
         real :: getdex

         getdex = this%dex

      end function getdex
!
      function getmg_prec(this)
         
         implicit none

         class(spect2d), intent(in) :: this
         real :: getmg_prec

         getmg_prec = this%mg_prec

      end function getmg_prec
!
      function getmg_type(this)
         
         implicit none

         class(spect2d), intent(in) :: this
         integer :: getmg_type

         getmg_type = this%mg_type

      end function getmg_type
!
      end module spect2d_class