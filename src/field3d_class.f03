! field3d class for QuickPIC Open Source 1.0
! update: 04/18/2016

      module field3d_class

      use perrors_class
      use parallel_pipe_class
      use spect3d_class
      use ufield3d_class
      use hdf5io_class
         
      implicit none

      private

      public :: field3d

      type field3d
! gcells = (0,1) = (no, yes) guard cell processing is performed

         private

         class(spect3d), pointer, public :: sp => null()
         class(perrors), pointer, public :: err => null()
         class(parallel_pipe), pointer, public :: p => null()
         class(ufield3d), pointer :: rs         
         integer :: gcells
         
         contains
         
         generic :: new => init_field3d
         generic :: del => end_field3d
         generic :: pcg => pipecg_field3d
         generic :: ag => acopyguard_field3d
         generic :: wr => writehdf5_3d, writehdf5_2dslice
         generic :: as => asc,asa
         generic :: add => sum
         generic :: sub => minus
         generic :: mult => multiply
         
         procedure, private :: init_field3d
         procedure, private :: end_field3d
         procedure, private :: pipecg_field3d, acopyguard_field3d
         procedure, private :: writehdf5_3d, writehdf5_2dslice
         procedure, private :: asc, asa, sum, minus, multiply
         procedure :: getgcells, getrs
                  
      end type 

      character(len=10), save :: class = 'field3d:'
      character(len=128), save :: erstr
      
      contains
!
      subroutine init_field3d(this,pp,perr,psp,dim)
      
         implicit none
         
         class(field3d), intent(inout) :: this
         class(spect3d), intent(in), pointer :: psp
         class(perrors), intent(in), pointer :: perr
         class(parallel_pipe), intent(in), pointer :: pp
         integer, intent(in) :: dim
! local data
         character(len=18), save :: sname = 'init_field3d:'
         
         this%sp => psp
         this%err => perr
         this%p => pp

         call this%err%werrfl2(class//sname//' started')

         allocate(this%rs)
         
         this%gcells = 0
         
         call this%rs%new(pp,perr,psp,dim,0,pp%getlnvp(),pp%getnstage())         

         call this%err%werrfl2(class//sname//' ended')

      end subroutine init_field3d
!
      subroutine end_field3d(this)
          
         implicit none
         
         class(field3d), intent(inout) :: this
         character(len=18), save :: sname = 'end_field3d:'

         call this%err%werrfl2(class//sname//' started')
         call this%rs%del()
         deallocate(this%rs)         
         call this%err%werrfl2(class//sname//' ended')
         
      end subroutine end_field3d
!
      subroutine pipecg_field3d(this,rtag,stag,rid,sid)
      
         implicit none
         
         class(field3d), intent(inout) :: this
         integer, intent(in) :: rtag,stag
         integer, intent(inout) :: rid,sid
         character(len=18), save :: sname = 'pipecg_field3d:'

         call this%err%werrfl2(class//sname//' started')
         
         call this%rs%pcg(rtag,stag,rid,sid)
         
         call this%err%werrfl2(class//sname//' ended')
         
      end subroutine pipecg_field3d
!      
      subroutine acopyguard_field3d(this,rtag,stag,id)
      
         implicit none
         
         class(field3d), intent(inout) :: this
         integer, intent(in) :: rtag,stag
         integer, intent(inout) :: id
         character(len=20), save :: sname = 'acopyguard_field3d:'

         call this%err%werrfl2(class//sname//' started')
         
         call this%rs%ag(rtag,stag,id)
         
         this%gcells = 1
         
         call this%err%werrfl2(class//sname//' ended')
         
      end subroutine acopyguard_field3d
!      
      subroutine writehdf5_3d(this,file,dim,rtag,stag,id)
      
         implicit none
         
         class(field3d), intent(inout) :: this
         class(hdf5file), intent(in) :: file
         integer, intent(in) :: dim, rtag, stag
         integer, intent(inout) :: id
         character(len=20), save :: sname = 'writehdf5_3d:'

         call this%err%werrfl2(class//sname//' started')
         
         call this%rs%wr(file,dim,rtag,stag,id)
                  
         call this%err%werrfl2(class//sname//' ended')
         
      end subroutine writehdf5_3d
!      
      subroutine writehdf5_2dslice(this,file,dim,slice,spos,rtag,stag,id)
      
         implicit none
         
         class(field3d), intent(inout) :: this
         class(hdf5file), intent(in) :: file
         integer, intent(in) :: dim, rtag, stag, slice, spos
         integer, intent(inout) :: id
         character(len=20), save :: sname = 'writehdf5_2dslice:'

         call this%err%werrfl2(class//sname//' started')
         
         call this%rs%wr(file,dim,slice,spos,rtag,stag,id)
                  
         call this%err%werrfl2(class//sname//' ended')
         
      end subroutine writehdf5_2dslice
!      
      subroutine asc(this,value) 

         implicit none
         
         class(field3d), intent(inout) :: this
         real, intent(in) :: value
         character(len=18), save :: sname = 'asc:'

         call this%err%werrfl2(class//sname//' started')
         
         call this%rs%as(value)         
                  
         call this%err%werrfl2(class//sname//' ended')
                 
      end subroutine asc
!
      subroutine asa(this,that) 

         implicit none
         
         class(field3d), intent(inout) :: this,that
         character(len=18), save :: sname = 'asc:'

         call this%err%werrfl2(class//sname//' started')
         
         call this%rs%as(that%rs)         
         
         call this%err%werrfl2(class//sname//' ended')
                 
      end subroutine asa
!
      subroutine sum(this,a1,a2) 

         implicit none
         
         class(field3d), intent(inout) :: this
         class(field3d), intent(in) :: a1, a2
         character(len=18), save :: sname = 'sum:'

         call this%err%werrfl2(class//sname//' started')
         
         call this%rs%add(a1%rs,a2%rs)
         
         call this%err%werrfl2(class//sname//' ended')
                 
      end subroutine sum
!
      subroutine minus(this,a1,a2) 

         implicit none
         
         class(field3d), intent(inout) :: this
         class(field3d), intent(in) :: a1, a2
         character(len=18), save :: sname = 'minus:'

         call this%err%werrfl2(class//sname//' started')
         
         call this%rs%sub(a1%rs,a2%rs)
         
         call this%err%werrfl2(class//sname//' ended')
                 
      end subroutine minus
!
      subroutine multiply(this,a1,value) 

         implicit none
         
         class(field3d), intent(inout) :: this
         class(field3d), intent(in) :: a1
         real, intent(in) :: value
         character(len=18), save :: sname = 'multiply:'

         call this%err%werrfl2(class//sname//' started')
         
         call this%rs%mult(a1%rs,value)
         
         call this%err%werrfl2(class//sname//' ended')
                 
      end subroutine multiply
!      
      function getgcells(this)

         implicit none

         class(field3d), intent(in) :: this
         integer :: getgcells
         
         getgcells = this%gcells

      end function getgcells      
!      
      function getrs(this)

         implicit none

         class(field3d), intent(in) :: this
         class(ufield3d), pointer :: getrs
         
         getrs => this%rs

      end function getrs
!      
      end module field3d_class