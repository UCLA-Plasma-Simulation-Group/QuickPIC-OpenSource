! ufield2d class for QuickPIC Open Source 1.0
! update: 04/18/2016

      module ufield2d_class

      use perrors_class
      use parallel_pipe_class
      use spect2d_class
      use ufield2d_lib
      use hdf5io_class
      use mpi
         
      implicit none

      private

      public :: ufield2d

      type ufield2d

         private
! layout values: XLOCAL = xy = 0, YLOCAL = yx = 1
! dim = dimension of the field
! nd1, nd2 = size of global array data in each dimension
! nvpx, nvpy = number of processors in each dimension
! nd1p, nd2p = size of local array data in each dimension (without guardcell)
! rf = pointer of the local 2d field array
! noff = smallest global gridpoint in y
!
         class(spect2d), pointer, public :: sp => null()
         class(perrors), pointer, public :: err => null()
         class(parallel_pipe), pointer, public :: p => null()
         integer :: layout, dim, noff
         integer :: nd1, nvpx, nd1p
         integer :: nd2, nvpy, nd2p
         real, dimension(:,:,:), pointer :: rf => null(), buff => null()

         contains
         
         generic :: new => init_ufield2d, init_ufield2d_k
         generic :: del => end_ufield2d
         generic :: cg => copyguard
         generic :: ag => acopyguard
         generic :: wr => writehdf5
         generic :: psend => pipesend_ufield2d
         generic :: precv => piperecv_ufield2d
!         generic :: assignment(=) => asc,asa
!         generic :: operator(+) => add
!         generic :: operator(*) => mult1
         generic :: as => asc, asa
         generic :: add => sum1, sum2
         generic :: sub => minus1, minus2
         generic :: mult => multiply1, multiply2
         final :: final_ufield2d
         
         procedure, private :: init_ufield2d
         procedure, private :: init_ufield2d_k
         procedure, private :: end_ufield2d
         procedure :: getlayout, getdim
         procedure :: getnd1, getnvpx, getnd1p
         procedure :: getnd2, getnvpy, getnd2p
         procedure :: getrf, getnoff
         procedure, private :: copyguard, acopyguard, writehdf5
         procedure, private :: pipesend_ufield2d, piperecv_ufield2d
         procedure, private :: asc, asa, sum1, multiply1, multiply2, minus1
         procedure, private :: sum2, minus2
!         procedure, private :: asc,asa,add,mult1
         
      end type 
      
      character(len=10), save :: class = 'ufield2d:'
      character(len=128), save :: erstr
! scr = guard cell buffer received from nearby processors
      real, dimension(:), allocatable, save  :: scr
      integer, save :: szscr = 0
      
      contains
!
      subroutine init_ufield2d(this,pp,perr,psp,dim,layout,nvpx,nvpy)
      
         implicit none
         
         class(ufield2d), intent(inout) :: this
         integer, intent(in) :: layout,dim
         integer, intent(in) :: nvpx,nvpy
         class(spect2d), intent(in), pointer :: psp
         class(perrors), intent(in), pointer :: perr
         class(parallel_pipe), intent(in), pointer :: pp
! local data
         character(len=18), save :: sname = 'init_ufield2d:'
         integer :: nd1,nd2
         
         call perr%werrfl2(class//sname//' started')
         this%sp => psp
         this%err => perr
         this%p => pp
         this%layout = layout
         this%dim = dim
         nd1 = 2**psp%getindx()
         nd2 = 2**psp%getindy()
         select case (layout)
         case (0)
            this%nd1 = nd1
            this%nd2 = nd2
         case (1)
            this%nd1 = nd2
            this%nd2 = nd1
         case default
            write (erstr,*) 'invalid layout = ', layout
            call this%err%equit(class//sname//erstr)
            return
         end select        
! make sure data is a multiple of the number of processors
         if ((((this%nd2/nvpy)*nvpy)/=this%nd2) .and. (((nvpy/this%nd2)*&
     &this%nd2)/=nvpy)) then
            write (erstr,*) 'data, proc number not multiples:', this%nd2&
     &, nvpy
            call this%err%equit(class//sname//erstr)
            return
         endif
! save number of processors in each dimension
         this%nvpx = 0
         this%nvpy = nvpy
         select case (layout)
         case (0)
            this%nd1p = nd1
            this%nd2p = nd2/nvpy
            this%noff = (pp%getlkstrt()-1)*nd2/nvpy
            allocate(this%rf(dim,this%nd1p+2,this%nd2p*2+psp%getinorder()))
            this%rf(:,:,:) = 0.0         
         case (1)
            this%nd1p = nd2/nvpy
            this%nd2p = nd1
            this%noff = (pp%getlkstrt()-1)*nd2/nvpy
            allocate(this%rf(dim,this%nd2p,this%nd1p+psp%getinorder()))         
            this%rf(:,:,:) = 0.0         
         case default
            write (erstr,*) 'invalid layout = ', layout
            call this%err%equit(class//sname//erstr)
            return
         end select        
         call perr%werrfl2(class//sname//' ended')
                  
      end subroutine init_ufield2d
!
      subroutine init_ufield2d_k(this,pp,perr,psp,dim,layout,nd1,nd2,nvpx,nvpy)
      
         implicit none
         
         class(ufield2d), intent(inout) :: this
         integer, intent(in) :: layout,dim
         integer, intent(in) :: nvpx,nvpy
         class(spect2d), intent(in), pointer :: psp
         class(perrors), intent(in), pointer :: perr
         class(parallel_pipe), intent(in), pointer :: pp
         integer, intent(in) :: nd1,nd2
! local data
         character(len=18), save :: sname = 'init_ufield2d_k:'
         
         call perr%werrfl2(class//sname//' started')
         this%sp => psp
         this%err => perr
         this%p => pp
         this%layout = layout
         this%dim = dim
         select case (layout)
         case (0)
            this%nd1 = nd1
            this%nd2 = nd2
         case (1)
            this%nd1 = nd2
            this%nd2 = nd1
         case default
            write (erstr,*) 'invalid layout = ', layout
            call this%err%equit(class//sname//erstr)
            return
         end select        
! make sure data is a multiple of the number of processors
         if ((((this%nd2/nvpy)*nvpy)/=this%nd2) .and. (((nvpy/this%nd2)*&
     &this%nd2)/=nvpy)) then
            write (erstr,*) 'data, proc number not multiples:', this%nd2&
     &, nvpy
            call this%err%equit(class//sname//erstr)
            return
         endif
! save number of processors in each dimension
         this%nvpx = 0
         this%nvpy = nvpy
         select case (layout)
         case (0)
            this%nd1p = nd1
            this%nd2p = nd2/nvpy
            allocate(this%rf(dim,this%nd1p,this%nd2p))
            this%rf(:,:,:) = 0.0         
         case (1)
            this%nd1p = nd2
            this%nd2p = (nd1-1)/nvpy + 1
            allocate(this%rf(dim,this%nd1p+1,this%nd2p+1))         
            this%rf(:,:,:) = 0.0         
         case default
            write (erstr,*) 'invalid layout = ', layout
            call this%err%equit(class//sname//erstr)
            return
         end select        
         call perr%werrfl2(class//sname//' ended')
                  
      end subroutine init_ufield2d_k
!
      subroutine end_ufield2d(this)
          
         implicit none
         
         class(ufield2d), intent(inout) :: this
! local data
         character(len=18), save :: sname = 'end_ufield2d:'

         call this%err%werrfl2(class//sname//' started')
         if (associated(this%rf)) deallocate(this%rf)
         call this%err%werrfl2(class//sname//' ended')
         
      end subroutine end_ufield2d
!      
      subroutine final_ufield2d(this)
          
         implicit none
         
         type(ufield2d), intent(inout) :: this
! local data
         character(len=18), save :: sname = 'final_ufield2d:'

         call this%err%werrfl2(class//sname//' started')
         if (associated(this%rf)) deallocate(this%rf)
         call this%err%werrfl2(class//sname//' ended')
         
      end subroutine final_ufield2d
!      
      function getlayout(this)

         implicit none

         class(ufield2d), intent(in) :: this
         integer :: getlayout
         
         getlayout = this%layout

      end function getlayout
!      
      function getnvpy(this)

         implicit none

         class(ufield2d), intent(in) :: this
         integer :: getnvpy
         
         getnvpy = this%nvpy

      end function getnvpy
!      
      function getnvpx(this)

         implicit none

         class(ufield2d), intent(in) :: this
         integer :: getnvpx
                  
         getnvpx = this%nvpx

      end function getnvpx
!      
      function getnd2(this)

         implicit none

         class(ufield2d), intent(in) :: this
         integer :: getnd2
         
         getnd2 = this%nd2

      end function getnd2
!      
      function getnd1(this)

         implicit none

         class(ufield2d), intent(in) :: this
         integer :: getnd1
         
         getnd1 = this%nd1

      end function getnd1
!      
      function getnd2p(this)

         implicit none

         class(ufield2d), intent(in) :: this
         integer :: getnd2p
         
         getnd2p = this%nd2p

      end function getnd2p
!      
      function getnd1p(this)

         implicit none

         class(ufield2d), intent(in) :: this
         integer :: getnd1p
         
         getnd1p = this%nd1p

      end function getnd1p
!      
      function getdim(this)

         implicit none

         class(ufield2d), intent(in) :: this
         integer :: getdim
         
         getdim = this%dim

      end function getdim
!      
      function getnoff(this)

         implicit none

         class(ufield2d), intent(in) :: this
         integer :: getnoff
         
         getnoff = this%noff

      end function getnoff
!      
      function getrf(this)

         implicit none

         class(ufield2d), intent(in) :: this
         real, dimension(:,:,:), pointer :: getrf
         
         getrf => this%rf

      end function getrf
!
      subroutine copyguard(this)
      
         implicit none
         
         class(ufield2d), intent(inout) :: this
         character(len=18), save :: sname = 'copyguard'
! local data         
         integer :: nxv, nypmx, nyp
         
         call this%err%werrfl2(class//sname//' started')
         nxv = size(this%rf,1)*size(this%rf,2)
         nypmx = size(this%rf,3)
         nyp = this%nd2p

         if (this%layout /= 0) then
            write (erstr,*) 'invalid layout = ', this%layout
            call this%err%equit(class//sname//erstr); return
         endif

         select case (this%sp%getpsolver())
         case (1)
            select case (this%sp%getinorder())
            case (1)
! copy data to guard cells in distributed direction
               call PPNCGUARD2L(this%rf,nyp,this%p%getlkstrt(),this%p%getlnvp(),&
               &nxv,nypmx,this%p%getlgrp(),this%p%getmreal())            
            case default
               call PPNCGUARD2L(this%rf,nyp,this%p%getlkstrt(),this%p%getlnvp(),&
               &nxv,nypmx,this%p%getlgrp(),this%p%getmreal())            
            end select
         case default
            call PPNCGUARD2L(this%rf,nyp,this%p%getlkstrt(),this%p%getlnvp(),&
            &nxv,nypmx,this%p%getlgrp(),this%p%getmreal())            
         end select
         call this%err%werrfl2(class//sname//' ended')
         
      end subroutine copyguard
!      
      subroutine acopyguard(this)
      
         implicit none
         
         class(ufield2d), intent(inout) :: this
         character(len=18), save :: sname = 'acopyguard'
! local data         
         integer :: ndim, nxv, nypmx, nyp, nx
         
         call this%err%werrfl2(class//sname//' started')
         ndim = size(this%rf,1)
         nxv = size(this%rf,2)
         nypmx = size(this%rf,3)
         nyp = this%nd2p
         nx = this%nd1

         if (szscr < ndim*nxv) then
            if (szscr /= 0) deallocate(scr)
! allocate new buffer
            allocate(scr(ndim*nxv))
            szscr = ndim*nxv
         endif
                  
         if (this%layout /= 0) then
            write (erstr,*) 'invalid layout = ', this%layout
            call this%err%equit(class//sname//erstr); return
         endif

         select case (this%sp%getpsolver())
         case (1)
            select case (this%sp%getinorder())
            case (1)
! copy data to guard cells in distributed direction
               call PPNACGUARD2L(this%rf,scr,nyp,nx,ndim,this%p%getlkstrt(),&
               &this%p%getlnvp(),nxv,nypmx,this%p%getlgrp(),this%p%getmreal())
            case default
               call PPNACGUARD2L(this%rf,scr,nyp,nx,ndim,this%p%getlkstrt(),&
               &this%p%getlnvp(),nxv,nypmx,this%p%getlgrp(),this%p%getmreal())
            end select
         case default
            call PPNACGUARD2L(this%rf,scr,nyp,nx,ndim,this%p%getlkstrt(),&
            &this%p%getlnvp(),nxv,nypmx,this%p%getlgrp(),this%p%getmreal())
         end select
         call this%err%werrfl2(class//sname//' ended')
         
      end subroutine acopyguard
!
      subroutine writehdf5(this,file,dim)

         implicit none

         class(ufield2d), intent(inout) :: this
         class(hdf5file), intent(in) :: file
         integer, intent(in) :: dim
! local data
         integer, dimension(2) :: gsize,lsize
         integer :: noff
         
         integer :: ierr
         real, dimension(:,:), pointer :: fdata
         character(len=10), save :: sname = 'writehdf5:'
         
         call this%err%werrfl2(class//sname//' started')
            noff = this%noff
            gsize =(/this%nd1,this%nd2/)
            lsize =(/this%nd1p,this%nd2p/)
            call pwfield(this%p,this%err,file,this%rf(dim,:,:),gsize,lsize,&
            &noff,ierr)
         call this%err%werrfl2(class//sname//' ended')

      end subroutine writehdf5
!
      subroutine pipesend_ufield2d(this,stag,id)
      
         implicit none
         
         class(ufield2d), intent(inout) :: this
         integer, intent(in) :: stag
         integer, intent(inout) :: id
! local data
         character(len=18), save :: sname = 'pipesend_ufield2d:'
         integer :: bsize,i,j,k,des,ierr
                  
         call this%err%werrfl2(class//sname//' started')

         des = this%p%getidproc()+this%p%getlnvp()
         
         if (des >= this%p%getnvp()) then
            id = MPI_REQUEST_NULL         
            call this%err%werrfl2(class//sname//' ended')
            return
         endif
         
         if (.not.associated(this%buff)) then
            select case (this%layout)
            case (0)
               allocate(this%buff(this%dim,this%nd1p+2,this%nd2p+this%sp%getinorder())) 
            case (1)        
               allocate(this%buff(this%dim,this%nd1p+1,this%nd2p+1)) 
            end select
         endif
         bsize = size(this%buff)

!$OMP PARALLEL DO PRIVATE(i,j,k)
         do i = 1, size(this%buff,1)
            do j = 1, size(this%buff,2)
               do k = 1, size(this%buff,3)       
                  this%buff(i,j,k) = this%rf(i,j,k)
               enddo
            enddo
         enddo 
!$OMP END PARALLEL DO

         call MPI_ISEND(this%buff,bsize,this%p%getmreal(),des,&
         &stag,this%p%getlworld(),id,ierr)


! check for errors
         if (ierr /= 0) then
            write (erstr,*) 'MPI_ISEND failed'
            call this%err%equit(class//sname//erstr); return
         endif
         
         
         call this%err%werrfl2(class//sname//' ended')
         
      end subroutine pipesend_ufield2d
!      
      subroutine piperecv_ufield2d(this,rtag)
      
         implicit none
         
         class(ufield2d), intent(inout) :: this
         integer, intent(in) :: rtag
! local data
         character(len=18), save :: sname = 'piperecv_ufield2d:'
         integer, dimension(10) :: istat
         integer :: bsize,id,i,j,k,des,ierr
         
         call this%err%werrfl2(class//sname//' started')

         des = this%p%getidproc()-this%p%getlnvp()
         
         if (des < 0) then
            this%rf(:,:,:) = 0.0
            call this%err%werrfl2(class//sname//' ended')
            return
         endif

!         if (.not.associated(this%buff)) then
!            select case (this%layout)
!            case (0)
!               allocate(this%buff(this%dim,this%nd1p+2,this%nd2p+this%sp%getinorder())) 
!            case (1)        
!               allocate(this%buff(this%dim,this%nd2p,this%nd1p+this%sp%getinorder())) 
!            end select
!         endif

         call MPI_IRECV(this%rf,size(this%rf),this%p%getmreal(),des,&
         &rtag,this%p%getlworld(),id,ierr)
         
         call MPI_WAIT(id,istat,ierr)

! check for errors
         if (ierr /= 0) then
            write (erstr,*) 'MPI failed'
            call this%err%equit(class//sname//erstr); return
         endif

         call this%err%werrfl2(class//sname//' ended')
         
      end subroutine piperecv_ufield2d
!
      subroutine sum1(this,a1,a2)
      
         implicit none
         
         class(ufield2d),intent(inout) :: this
         class(ufield2d), target, intent(in) :: a1,a2
! local data
         character(len=18), save :: sname = 'sum1:'
         integer :: i,j,k
         real, dimension(:,:,:), pointer :: rf1 => null(), rf2 => null(),&
         &rf3 => null()
         
         call this%err%werrfl2(class//sname//' started')
         
         rf1 => this%rf         
         rf2 => a1%rf
         rf3 => a2%rf
         

!$OMP PARALLEL DO PRIVATE(i,j,k)         
         do k = 1, size(rf1,3)
            do j = 1, size(rf1,2)
               do i = 1, size(rf1,1)
                  rf1(i,j,k) = rf2(i,j,k) + rf3(i,j,k)
               enddo
            enddo
         enddo
!$OMP END PARALLEL DO    

         call this%err%werrfl2(class//sname//' ended')

      end subroutine sum1
!
      subroutine sum2(this,a1,a2,dim,dim1,dim2)
      
         implicit none
         
         class(ufield2d),intent(inout) :: this
         class(ufield2d), target, intent(in) :: a1,a2
         integer, dimension(:), intent(in) :: dim,dim1, dim2
! local data
         character(len=18), save :: sname = 'sum2:'
         integer :: i,j,k
         real, dimension(:,:,:), pointer :: rf1 => null(), rf2 => null(),&
         &rf3 => null()
         
         call this%err%werrfl2(class//sname//' started')
         
         rf1 => this%rf         
         rf2 => a1%rf
         rf3 => a2%rf
         

!$OMP PARALLEL DO PRIVATE(i,j,k)         
         do k = 1, size(rf1,3)
            do j = 1, size(rf1,2)
               do i = 1, size(dim1)
                  rf1(dim(i),j,k) = rf2(dim1(i),j,k) + rf3(dim2(i),j,k)
               enddo
            enddo
         enddo
!$OMP END PARALLEL DO    

         call this%err%werrfl2(class//sname//' ended')

      end subroutine sum2
!
      subroutine minus1(this,a1,a2)
      
         implicit none
         
         class(ufield2d),intent(inout) :: this
         class(ufield2d), target, intent(in) :: a1,a2
! local data
         character(len=18), save :: sname = 'minus1:'
         integer :: i,j,k
         real, dimension(:,:,:), pointer :: rf1 => null(), rf2 => null(),&
         &rf3 => null()
         
         call this%err%werrfl2(class//sname//' started')
         
         rf1 => this%rf         
         rf2 => a1%rf
         rf3 => a2%rf
         

!$OMP PARALLEL DO PRIVATE(i,j,k)         
         do k = 1, size(rf1,3)
            do j = 1, size(rf1,2)
               do i = 1, size(rf1,1)
                  rf1(i,j,k) = rf2(i,j,k) - rf3(i,j,k)
               enddo
            enddo
         enddo
!$OMP END PARALLEL DO    

         call this%err%werrfl2(class//sname//' ended')

      end subroutine minus1
!
      subroutine minus2(this,a1,a2,dim,dim1,dim2)
      
         implicit none
         
         class(ufield2d),intent(inout) :: this
         class(ufield2d), target, intent(in) :: a1,a2
         integer, dimension(:), intent(in) :: dim, dim1, dim2         
! local data
         character(len=18), save :: sname = 'minus1:'
         integer :: i,j,k
         real, dimension(:,:,:), pointer :: rf1 => null(), rf2 => null(),&
         &rf3 => null()
         
         call this%err%werrfl2(class//sname//' started')
         
         rf1 => this%rf         
         rf2 => a1%rf
         rf3 => a2%rf
         

!$OMP PARALLEL DO PRIVATE(i,j,k)         
         do k = 1, size(rf1,3)
            do j = 1, size(rf1,2)
               do i = 1, size(dim1)
                  rf1(dim(i),j,k) = rf2(dim1(i),j,k) - rf3(dim2(i),j,k)
               enddo
            enddo
         enddo
!$OMP END PARALLEL DO    

         call this%err%werrfl2(class//sname//' ended')

      end subroutine minus2
!      
      subroutine multiply1(this,a1,value)
      
         implicit none
         
         class(ufield2d), intent(inout) :: this
         class(ufield2d), target, intent(in) :: a1
         real, intent(in) :: value
! local data
         character(len=18), save :: sname = 'multiply1:'
         integer :: i,j,k
         real, dimension(:,:,:), pointer :: rf1 => null(), rf2 => null()
         
         call this%err%werrfl2(class//sname//' started')

         rf1 => this%rf         
         rf2 => a1%rf

!$OMP PARALLEL DO PRIVATE(i,j,k)         
         do k = 1, size(rf1,3)
            do j = 1, size(rf1,2)
               do i = 1, size(rf1,1)
                  rf1(i,j,k) = rf2(i,j,k) * value
               enddo
            enddo
         enddo
!$OMP END PARALLEL DO        
         call this%err%werrfl2(class//sname//' ended')

      end subroutine multiply1
!
      subroutine multiply2(this,a1,dim,dim1,value)
      
         implicit none
         
         class(ufield2d), intent(inout) :: this
         class(ufield2d), target, intent(in) :: a1
         integer, dimension(:), intent(in) :: dim, dim1
         real, dimension(:), intent(in) :: value
! local data
         character(len=18), save :: sname = 'multiply2:'
         integer :: i,j,k
         real, dimension(:,:,:), pointer :: rf1 => null(), rf2 => null()
         
         call this%err%werrfl2(class//sname//' started')

         rf1 => this%rf         
         rf2 => a1%rf

!$OMP PARALLEL DO PRIVATE(i,j,k)         
         do k = 1, size(rf1,3)
            do j = 1, size(rf1,2)
               do i = 1, size(dim,1)
                  rf1(dim(i),j,k) = rf2(dim1(i),j,k) * value(i)
               enddo
            enddo
         enddo
!$OMP END PARALLEL DO        
         call this%err%werrfl2(class//sname//' ended')

      end subroutine multiply2
!
      subroutine asc(this,value)
      
         implicit none
         
         class(ufield2d), intent(inout) :: this
         real, intent(in) :: value
! local data
         character(len=18), save :: sname = 'asc:'
         integer :: i,j,k
         real, dimension(:,:,:), pointer :: rf => null()
         
         call this%err%werrfl2(class//sname//' started')
         
         rf => this%rf

!$OMP PARALLEL DO PRIVATE(i,j,k)         
         do k = 1, size(rf,3)
            do j = 1, size(rf,2)
               do i = 1, size(rf,1)
                  rf(i,j,k) = value
               enddo
            enddo
         enddo
!$OMP END PARALLEL DO
                    
         call this%err%werrfl2(class//sname//' ended')

      end subroutine asc
!
      subroutine asa(this,that)
      
         implicit none
         
         class(ufield2d), intent(inout) :: this
         class(ufield2d), target, intent(in) :: that
! local data
         character(len=18), save :: sname = 'asa:'
         integer :: i,j,k
         real, dimension(:,:,:), pointer :: rf1 => null(), rf2 => null()
         
         call this%err%werrfl2(class//sname//' started')
         
         rf1 => this%rf
         rf2 => that%rf

!$OMP PARALLEL DO PRIVATE(i,j,k)         
         do k = 1, size(rf1,3)
            do j = 1, size(rf1,2)
               do i = 1, size(rf1,1)
                  rf1(i,j,k) = rf2(i,j,k)
               enddo
            enddo
         enddo
!$OMP END PARALLEL DO
!         deallocate(rf2)

         call this%err%werrfl2(class//sname//' ended')

      end subroutine asa
!
      function add(this,that)
      
         implicit none
         
         class(ufield2d),intent(in) :: this,that
         class(ufield2d), allocatable :: add
! local data
         character(len=18), save :: sname = 'add:'
         integer :: i,j,k
         real, dimension(:,:,:), pointer :: rf1 => null(), rf2 => null(),&
         &rf3 => null()
         
         call this%err%werrfl2(class//sname//' started')
         
         allocate(add,source=this)
         
         rf2 => this%rf
         rf3 => that%rf
         
         allocate(add%rf(size(rf2,1),size(rf2,2),size(rf2,3)))

         rf1 => add%rf

!$OMP PARALLEL DO PRIVATE(i,j,k)         
         do k = 1, size(rf1,3)
            do j = 1, size(rf1,2)
               do i = 1, size(rf1,1)
                  rf1(i,j,k) = rf2(i,j,k) + rf3(i,j,k)
               enddo
            enddo
         enddo
!$OMP END PARALLEL DO    

!         deallocate(add%rf)

         
         call this%err%werrfl2(class//sname//' ended')

      end function add
!
      function mult1(this,value)
      
         implicit none
         
         class(ufield2d), intent(in) :: this
         real, intent(in) :: value
         class(ufield2d), allocatable :: mult1
! local data
         character(len=18), save :: sname = 'mult1:'
         integer :: i,j,k
         real, dimension(:,:,:), pointer :: rf1 => null(), rf2 => null()
         
         call this%err%werrfl2(class//sname//' started')

         allocate(mult1,source=this)
         
         rf2 => this%rf

         allocate(mult1%rf(size(rf2,1),size(rf2,2),size(rf2,3)))
         rf1 => mult1%rf

!$OMP PARALLEL DO PRIVATE(i,j,k)         
         do k = 1, size(rf1,3)
            do j = 1, size(rf1,2)
               do i = 1, size(rf1,1)
                  rf1(i,j,k) = rf2(i,j,k) * value
               enddo
            enddo
         enddo
!$OMP END PARALLEL DO        
         call this%err%werrfl2(class//sname//' ended')

      end function mult1
!
      function mult2(value,this)
      
         implicit none
         
         class(ufield2d), intent(in) :: this
         real, intent(in) :: value
         class(ufield2d), allocatable :: mult2
! local data
         character(len=18), save :: sname = 'mult2:'
         integer :: i,j,k
         real, dimension(:,:,:), pointer :: rf1 => null(), rf2 => null()
         
         call this%err%werrfl2(class//sname//' started')

         allocate(mult2,source=this)
         
         rf2 => this%rf

         allocate(mult2%rf(size(rf2,1),size(rf2,2),size(rf2,3)))
         rf1 => mult2%rf

!$OMP PARALLEL DO PRIVATE(i,j,k)         
         do k = 1, size(rf1,3)
            do j = 1, size(rf1,2)
               do i = 1, size(rf1,1)
                  rf1(i,j,k) = rf2(i,j,k) * value
               enddo
            enddo
         enddo
!$OMP END PARALLEL DO        
         call this%err%werrfl2(class//sname//' ended')

      end function mult2
!
      end module ufield2d_class