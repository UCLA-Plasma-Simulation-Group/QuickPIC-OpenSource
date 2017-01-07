! Spect3d class for QuickPIC Open Source 1.0
! update: 04/18/2016

      module spect3d_class

      use perrors_class
      use parallel_pipe_class
      use spect2d_class
         
      implicit none

      private

      public :: spect3d

      type, extends(spect2d) :: spect3d

         private
         
         integer :: indz
         
         contains
         
         procedure, private :: init_spect3d
         procedure, private :: end_spect2d => end_spect3d
         generic :: new => init_spect3d
!         generic :: del => end_spect3d
         procedure :: getindz
                  
      end type spect3d
      
      contains
!
      subroutine init_spect3d(this,pp,perr,indx,indy,indz,psolver,inorder)
      
         implicit none
         
         class(spect3d), intent(inout) :: this
         class(perrors), intent(in), pointer :: perr
         class(parallel_pipe), intent(in), pointer :: pp
         integer, intent(in) :: indx, indy, indz,psolver, inorder
         
         call this%spect2d%new(pp,perr,indx,indy,psolver,inorder)
         this%indz = indz

      end subroutine init_spect3d
!
      subroutine end_spect3d(this)
          
         implicit none
         
         class(spect3d), intent(inout) :: this
         
         return
         
      end subroutine end_spect3d
!      
      function getindz(this)

         implicit none

         class(spect3d), intent(in) :: this
         integer :: getindz
         
         getindz = this%indz

      end function getindz        
!      
      end module spect3d_class