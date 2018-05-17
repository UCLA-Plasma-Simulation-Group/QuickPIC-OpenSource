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
         integer :: indx, indy, psolver, inorder
         class(perrors), pointer, public :: err => null()
         class(parallel_pipe), pointer, public :: p => null()
         
         contains
         
         procedure, private :: init_spect2d
         procedure, private :: end_spect2d
         generic :: new => init_spect2d
         generic :: del => end_spect2d
         procedure :: getindx
         procedure :: getindy
         procedure :: getpsolver
         procedure :: getinorder
                  
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
      end module spect2d_class