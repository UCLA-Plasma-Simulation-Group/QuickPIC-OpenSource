! beam3d class for QuickPIC Open Source 1.0
! update: 04/18/2016

      module beam3d_class

      use perrors_class
      use parallel_pipe_class
      use spect3d_class
      use fdist3d_class
      use field3d_class
      use field2d_class
      use part3d_class
      use hdf5io_class
      use mpi
               
      implicit none

      private

      public :: beam3d

      type beam3d

         private

         class(spect3d), pointer, public :: sp => null()
         class(perrors), pointer, public :: err => null()
         class(parallel_pipe), pointer, public :: p => null()
         class(part3d), pointer :: pd
         class(field3d), pointer :: q => null()
         class(fdist3d), pointer :: pf => null()
         logical :: evol
         contains
         
         generic :: new => init_beam3d
         generic :: del => end_beam3d
         generic :: push => push_beam3d
         generic :: pmv => pmove_beam3d
         generic :: qdp => qdeposit_beam3d, qdpcopy_beam3d  
         generic :: wr => writehdf5_beam3d
         generic :: wrq => writeq_beam3d, writeqslice_beam3d
         generic :: wrst => writerst_beam3d       
         generic :: rrst => readrst_beam3d       
         procedure, private :: init_beam3d
         procedure, private :: end_beam3d
         procedure, private :: push_beam3d
         procedure, private :: pmove_beam3d   
         procedure, private :: qdeposit_beam3d, writehdf5_beam3d
         procedure, private :: writerst_beam3d, readrst_beam3d
         procedure, private :: writeq_beam3d, writeqslice_beam3d 
         procedure, private :: qdpcopy_beam3d                 
      end type 

      save      

      character(len=10) :: class = 'beam3d:'
      character(len=128) :: erstr
      
      contains
!
      subroutine init_beam3d(this,pp,perr,psp,pf,qbm,dt,ci,xdim)
      
         implicit none
         
         class(beam3d), intent(inout) :: this
         class(spect3d), intent(in), pointer :: psp
         class(perrors), intent(in), pointer :: perr
         class(parallel_pipe), intent(in), pointer :: pp
         class(fdist3d), intent(inout), target :: pf
         real, intent(in) :: qbm, dt, ci
         integer, intent(in) :: xdim

! local data
         character(len=18), save :: sname = 'init_beam3d:'
         integer :: id, ierr
         integer, dimension(10) :: istat
                  
         this%sp => psp
         this%err => perr
         this%p => pp
         this%pf => pf

         call this%err%werrfl2(class//sname//' started')

         allocate(this%pd,this%q)
         this%evol = pf%getevol()
         call this%q%new(this%p,this%err,this%sp,dim=1)
         call this%pd%new(pp,perr,psp,pf,this%q%getrs(),qbm,dt,ci,xdim)
         call this%pmv(this%q,1,1,id)
         call MPI_WAIT(id,istat,ierr)
         
         call this%err%werrfl2(class//sname//' ended')

      end subroutine init_beam3d
!
      subroutine end_beam3d(this)
          
         implicit none
         
         class(beam3d), intent(inout) :: this
         character(len=18), save :: sname = 'end_beam3d:'

         call this%err%werrfl2(class//sname//' started')
         call this%pd%del()
         call this%q%del()
         call this%err%werrfl2(class//sname//' ended')
                  
      end subroutine end_beam3d
!      
      subroutine qdeposit_beam3d(this,id1,id2,id3,tag1,tag2)
! deposit the charge density      
      
         implicit none
         
         class(beam3d), intent(inout) :: this
         integer, intent(inout) :: id1, id2, id3, tag1, tag2
! local data
         character(len=18), save :: sname = 'qdeposit_beam3d:'
         integer, dimension(10) :: istat
         integer :: ierr
                  
         call this%err%werrfl2(class//sname//' started')
         call this%q%as(0.0)
         call MPI_WAIT(id1,istat,ierr)
         call MPI_WAIT(id3,istat,ierr)
         call this%pd%qdp(this%q%getrs())
         call this%q%ag(tag1,tag1,id1)
         call this%q%pcg(tag2,tag2,id2,id3)    
                 
         call this%err%werrfl2(class//sname//' ended')
         
      end subroutine qdeposit_beam3d
!      
      subroutine qdpcopy_beam3d(this,q,slice)
! copy and add the charge density to a 2d slice
      
         implicit none
         
         class(beam3d), intent(inout) :: this
         class(field2d), intent(inout) :: q
         integer, intent(in) :: slice
! local data
         character(len=18), save :: sname = 'qdpcopy_beam3d:'
                  
         call this%err%werrfl2(class//sname//' started')

         call q%ca(this%q,slice,(/1/),(/1/))
                
         call this%err%werrfl2(class//sname//' ended')
         
      end subroutine qdpcopy_beam3d
!      
      subroutine push_beam3d(this,ef,bf,dex,dez,rtag,stag,sid)
      
         implicit none
         
         class(beam3d), intent(inout) :: this
         class(field3d), intent(in) :: ef, bf
         real, intent(in) :: dex, dez
         integer, intent(in) :: rtag, stag
         integer, intent(inout) :: sid         
! local data
         character(len=18), save :: sname = 'partpush'

         call this%err%werrfl2(class//sname//' started')

         if (.not. this%evol) then
            call this%err%werrfl2(class//sname//' ended')
            return
         end if

         call this%pd%push(ef%getrs(),bf%getrs(),dex,dez)
         
         call this%pmv(ef,rtag,stag,sid)
         
         call this%err%werrfl2(class//sname//' ended')
         
      end subroutine push_beam3d
!
      subroutine pmove_beam3d(this,fd,rtag,stag,sid)
      
         implicit none
         
         class(beam3d), intent(inout) :: this
         class(field3d), intent(in) :: fd
         integer, intent(in) :: rtag, stag
         integer, intent(inout) :: sid
! local data
         character(len=18), save :: sname = 'pmove:'
         
         call this%err%werrfl2(class//sname//' started')
         
         call this%pd%pmv(fd%getrs(),rtag,stag,sid)
         
         call this%err%werrfl2(class//sname//' ended')
         
      end subroutine pmove_beam3d
!
      subroutine writehdf5_beam3d(this,file,dspl,delta,rtag,stag,id)

         implicit none
         
         class(beam3d), intent(inout) :: this
         class(hdf5file), intent(in) :: file
         real, dimension(3), intent(in) :: delta
         integer, intent(in) :: dspl, rtag, stag
         integer, intent(inout) :: id
! local data
         character(len=18), save :: sname = 'writehdf5_beam3d:'

         call this%err%werrfl2(class//sname//' started')                  
         call this%pd%wr(file,dspl,delta,rtag,stag,id)
         call this%err%werrfl2(class//sname//' ended')
      
      end subroutine writehdf5_beam3d
!            
      subroutine writerst_beam3d(this,file)

         implicit none
         
         class(beam3d), intent(inout) :: this
         class(hdf5file), intent(in) :: file
! local data
         character(len=18), save :: sname = 'writerst_beam3d:'

         call this%err%werrfl2(class//sname//' started')                  
         call this%pd%wrst(file)
         call this%err%werrfl2(class//sname//' ended')
      
      end subroutine writerst_beam3d
!            
      subroutine writeq_beam3d(this,file,rtag,stag,id)

         implicit none
         
         class(beam3d), intent(inout) :: this
         class(hdf5file), intent(in) :: file
         integer, intent(in) :: rtag, stag
         integer, intent(inout) :: id
! local data
         character(len=18), save :: sname = 'writeq_beam3d:'

         call this%err%werrfl2(class//sname//' started')                  
         call this%q%wr(file,1,rtag,stag,id)
         call this%err%werrfl2(class//sname//' ended')
      
      end subroutine writeq_beam3d
!            
      subroutine writeqslice_beam3d(this,file,slice,spos,rtag,stag,id)

         implicit none
         
         class(beam3d), intent(inout) :: this
         class(hdf5file), intent(in) :: file
         integer, intent(in) :: rtag, stag, slice, spos
         integer, intent(inout) :: id         
! local data
         character(len=18), save :: sname = 'writeqslice_beam3d:'

         call this%err%werrfl2(class//sname//' started')                  
         call this%q%wr(file,1,slice,spos,rtag,stag,id)
         call this%err%werrfl2(class//sname//' ended')
      
      end subroutine writeqslice_beam3d
!            
      subroutine readrst_beam3d(this,file)

         implicit none
         
         class(beam3d), intent(inout) :: this
         class(hdf5file), intent(in) :: file
! local data
         character(len=18), save :: sname = 'readrst_beam3d:'

         call this%err%werrfl2(class//sname//' started')                  
         call this%pd%rrst(file)
         call this%err%werrfl2(class//sname//' ended')
      
      end subroutine readrst_beam3d
!            
      end module beam3d_class