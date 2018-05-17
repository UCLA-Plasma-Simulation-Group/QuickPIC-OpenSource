! ufield3d class for QuickPIC Open Source 1.0
! update: 04/18/2016

      module ufield3d_class

      use perrors_class
      use parallel_pipe_class
      use spect3d_class
      use ufield3d_lib
      use ufield2d_class
      use hdf5io_class
      use mpi
         
      implicit none

      private

      public :: ufield3d

      type ufield3d

         private
! dim = dimension of the field
! nd1, nd2, nd3 = size of global array data in each dimension
! nvpx, nvpy, nvpz = number of processors in each dimension
! nd1p, nd2p, nd3p = size of local array data in each dimension (without guardcell)
! rf = pointer of the local 3d field array
! noff = smallest global gridpoint in y
!
         class(spect3d), pointer, public :: sp => null()
         class(perrors), pointer, public :: err => null()
         class(parallel_pipe), pointer, public :: p => null()
         integer :: dim
         integer, dimension(2) :: noff
         integer :: nd1, nvpx, nd1p
         integer :: nd2, nvpy, nd2p
         integer :: nd3, nvpz, nd3p
         real, dimension(:,:,:,:), pointer :: rf => null()

         contains
         
         generic :: new => init_ufield3d
         generic :: del => end_ufield3d
         generic :: pcg => copyguard_pipe
         generic :: ag => acopyguard
         generic :: cp => copyin
         generic :: cb => copyout
         generic :: wr => writehdf5_3d, writehdf5_2dslice
         generic :: as => asc,asa
         generic :: add => sum
         generic :: sub => minus
         generic :: mult => multiply
         final :: final_ufield3d
         procedure, private :: init_ufield3d
         procedure, private :: end_ufield3d
         procedure :: getdim
         procedure :: getnd1, getnvpx, getnd1p
         procedure :: getnd2, getnvpy, getnd2p
         procedure :: getnd3, getnvpz, getnd3p
         procedure :: getrf, getnoff
         procedure, private :: copyguard_pipe, acopyguard
         procedure, private :: copyin, copyout
         procedure, private :: writehdf5_3d, writehdf5_2dslice
         procedure, private :: asc, asa, sum, minus, multiply

      end type 
      
      character(len=10), save :: class = 'ufield3d:'
      character(len=128), save :: erstr
      
      contains
!
      subroutine init_ufield3d(this,pp,perr,psp,dim,nvpx,nvpy,nvpz)
      
         implicit none
         
         class(ufield3d), intent(inout) :: this
         integer, intent(in) :: dim
         integer, intent(in) :: nvpx,nvpy,nvpz
         class(spect3d), intent(in), pointer :: psp
         class(perrors), intent(in), pointer :: perr
         class(parallel_pipe), intent(in), pointer :: pp
! local data
         character(len=18), save :: sname = 'init_ufield3d:'
         integer :: nd1,nd2,nd3
         
         call perr%werrfl2(class//sname//' started')
         this%sp => psp
         this%err => perr
         this%p => pp
         this%dim = dim
         nd1 = 2**psp%getindx()
         nd2 = 2**psp%getindy()
         nd3 = 2**psp%getindz()
         this%nd1 = nd1
         this%nd2 = nd2
         this%nd3 = nd3
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
         this%nvpz = nvpz
         this%nd1p = nd1
         this%nd2p = nd2/nvpy
         this%nd3p = nd3/nvpz
         this%noff(1) = (pp%getlkstrt()-1)*nd2/nvpy
         this%noff(2) = pp%getstageid()*nd3/pp%getnstage()
         allocate(this%rf(dim,this%nd1p+2,this%nd2p+1,this%nd3p+1))
         this%rf(:,:,:,:) = 0.0         
         call perr%werrfl2(class//sname//' ended')
                  
      end subroutine init_ufield3d
!
      subroutine end_ufield3d(this)
          
         implicit none
         
         class(ufield3d), intent(inout) :: this
! local data
         character(len=18), save :: sname = 'end_ufield3d:'

         call this%err%werrfl2(class//sname//' started')
         if (associated(this%rf)) deallocate(this%rf)
         call this%err%werrfl2(class//sname//' ended')
         
      end subroutine end_ufield3d
!      
      subroutine final_ufield3d(this)
          
         implicit none
         
         type(ufield3d), intent(inout) :: this
! local data
         character(len=18), save :: sname = 'final_ufield3d:'

         call this%err%werrfl2(class//sname//' started')
         if (associated(this%rf)) deallocate(this%rf)
         call this%err%werrfl2(class//sname//' ended')
         
      end subroutine final_ufield3d
!      
      function getnvpy(this)

         implicit none

         class(ufield3d), intent(in) :: this
         integer :: getnvpy
         
         getnvpy = this%nvpy

      end function getnvpy
!      
      function getnvpx(this)

         implicit none

         class(ufield3d), intent(in) :: this
         integer :: getnvpx
                  
         getnvpx = this%nvpx

      end function getnvpx
!      
      function getnvpz(this)

         implicit none

         class(ufield3d), intent(in) :: this
         integer :: getnvpz
         
         getnvpz = this%nvpz

      end function getnvpz

!      
      function getnd2(this)

         implicit none

         class(ufield3d), intent(in) :: this
         integer :: getnd2
         
         getnd2 = this%nd2

      end function getnd2
!      
      function getnd1(this)

         implicit none

         class(ufield3d), intent(in) :: this
         integer :: getnd1
         
         getnd1 = this%nd1

      end function getnd1
!      
      function getnd3(this)

         implicit none

         class(ufield3d), intent(in) :: this
         integer :: getnd3
         
         getnd3 = this%nd3

      end function getnd3
!      
      function getnd3p(this)

         implicit none

         class(ufield3d), intent(in) :: this
         integer :: getnd3p
         
         getnd3p = this%nd3p

      end function getnd3p
!      
      function getnd2p(this)

         implicit none

         class(ufield3d), intent(in) :: this
         integer :: getnd2p
         
         getnd2p = this%nd2p

      end function getnd2p
!      
      function getnd1p(this)

         implicit none

         class(ufield3d), intent(in) :: this
         integer :: getnd1p
         
         getnd1p = this%nd1p

      end function getnd1p
!      
      function getdim(this)

         implicit none

         class(ufield3d), intent(in) :: this
         integer :: getdim
         
         getdim = this%dim

      end function getdim
!      
      function getnoff(this)

         implicit none

         class(ufield3d), intent(in) :: this
         integer, dimension(2) :: getnoff
         
         getnoff = this%noff

      end function getnoff
!      
      function getrf(this)

         implicit none

         class(ufield3d), intent(in) :: this
         real, dimension(:,:,:,:), pointer :: getrf
         
         getrf => this%rf

      end function getrf
!
      subroutine copyguard_pipe(this,rtag,stag,rid,sid)
      
         implicit none
         
         class(ufield3d), intent(inout) :: this
         integer, intent(in) :: rtag, stag 
         integer, intent(inout) :: rid, sid
! local data
         integer :: nvpy, nvpz, kyp, kzp, ngds = 1
         real, dimension(:,:,:,:), pointer :: f
         character(len=18), save :: sname = 'copyguard:'
         integer :: nxv, nypmx, nzpmx, ierr
         real, dimension(:,:,:), pointer :: scs
         
         call this%err%werrfl2(class//sname//' started')

         f => this%rf
         nxv = size(f,1)*size(f,2); nypmx = size(f,3); nzpmx = size(f,4)
         nvpy=this%p%getlnvp(); kyp = this%nd2p; kzp= this%nd3p
         nvpz=this%p%getnstage()
         allocate(scs(size(f,1)*size(f,2),size(f,4),2*ngds))
         
         select case (this%sp%getpsolver())
         case (1)
            select case (this%sp%getinorder())
            case (1)
               call PCGUARD32L(f,scs,this%p%getkstrt(),nvpy,nvpz,nxv,nypmx,nzpmx,&
               &1,1,kyp,kzp,ngds,rtag,stag,rid,sid,ierr)
               if (this%p%getstageid() == 0) then
                  sid = MPI_REQUEST_NULL         
               endif
               if (this%p%getstageid() == this%p%getnstage()-1) then
                  rid = MPI_REQUEST_NULL         
               endif
            case default
               call PCGUARD32L(f,scs,this%p%getkstrt(),nvpy,nvpz,nxv,nypmx,nzpmx,&
               &1,1,kyp,kzp,ngds,rtag,stag,rid,sid,ierr)
               if (this%p%getstageid() == 0) then
                  sid = MPI_REQUEST_NULL         
               endif
               if (this%p%getstageid() == this%p%getnstage()-1) then
                  rid = MPI_REQUEST_NULL         
               endif
            end select
         case default
               call PCGUARD32L(f,scs,this%p%getkstrt(),nvpy,nvpz,nxv,nypmx,nzpmx,&
               &1,1,kyp,kzp,ngds,rtag,stag,rid,sid,ierr)
               if (this%p%getstageid() == 0) then
                  sid = MPI_REQUEST_NULL         
               endif
               if (this%p%getstageid() == this%p%getnstage()-1) then
                  rid = MPI_REQUEST_NULL         
               endif
         end select
         call this%err%werrfl2(class//sname//' ended')
         deallocate(scs)

      end subroutine copyguard_pipe
!      
      subroutine acopyguard(this,rtag,stag,id)
      
         implicit none
         
         class(ufield3d), intent(inout) :: this
         integer, intent(in) :: rtag, stag 
         integer, intent(inout) :: id
! local data
         integer :: nvpy, nvpz, nx, kyp, kzp, ngds = 1
         real, dimension(:,:,:,:), pointer :: f
         character(len=18), save :: sname = 'acopyguard:'
         integer :: nxv, nypmx, nzpmx, ierr
         real, dimension(:,:,:,:), pointer :: scs, scr

         call this%err%werrfl2(class//sname//' started')
         
         f => this%rf
         nxv = size(f,2); nypmx = size(f,3); nzpmx = size(f,4)
         nvpy=this%p%getlnvp(); kyp = this%nd2p; kzp= this%nd3p
         nvpz=this%p%getnstage(); nx = this%nd1p
         if (size(f,1) == 1) then
            allocate(scs(size(f,2),size(f,4),2*ngds,1),scr(size(f,2),size(f,3),&
            &ngds,1))
         else if (size(f,1) == 3) then
            allocate(scs(size(f,1),size(f,2),size(f,4),2*ngds))
         end if
            
         select case (this%sp%getpsolver())
         case (1)
            select case (this%sp%getinorder())
            case (1)
               if (size(f,1) == 1) then
                  call PAGUARD32L(f,scs,scr,this%p%getkstrt(),nvpy,nvpz,nx,nxv,&
                  &nypmx,nzpmx,1,1,kyp,kzp,ngds,rtag,stag,id,ierr)
                  if (this%p%getstageid() == this%p%getnstage() - 1) then
                     id = MPI_REQUEST_NULL         
                  endif
               else if (size(f,1) == 3) then
                  call PACGUARD32L(f,scs,this%p%getlkstrt(),nvpy,nvpz,nx,nxv,&
                  &nypmx,nzpmx,1,1,kyp,kzp,ngds)
               end if
            case default
                  call PAGUARD32L(f,scs,scr,this%p%getkstrt(),nvpy,nvpz,nx,nxv,&
                  &nypmx,nzpmx,1,1,kyp,kzp,ngds,rtag,stag,id,ierr)               
            end select
         case default
                  call PAGUARD32L(f,scs,scr,this%p%getkstrt(),nvpy,nvpz,nx,nxv,&
                  &nypmx,nzpmx,1,1,kyp,kzp,ngds,rtag,stag,id,ierr)               
         end select
         deallocate(scs,scr)
         call this%err%werrfl2(class//sname//' ended')
         
      end subroutine acopyguard
!     
      subroutine copyin(this,fd2d,lpos,sdim,ddim) 
! copy 2d field to 3d field at the local longitudinal position lpos
         implicit none
         
         class(ufield3d), intent(inout) :: this
         class(ufield2d), intent(in), target :: fd2d
         integer, intent(in) :: lpos
         integer, dimension(:), intent(in):: sdim, ddim
! local data         
         character(len=18), save :: sname = 'copyin:'
         real, dimension(:,:,:), pointer :: rf2d
         integer :: i,j,k,rank,nd1,nd2

         call this%err%werrfl2(class//sname//' started')
         
         rf2d => fd2d%getrf()
         rank = size(sdim)
         nd1 = size(this%rf,2)
         nd2 = size(this%rf,3)

!$OMP PARALLEL DO PRIVATE(i,j,k)
         do k = 1, nd2
            do j = 1, nd1
               do i = 1, rank
                  this%rf(ddim(i),j,k,lpos) = rf2d(sdim(i),j,k)
               enddo
            enddo
         enddo
!$OMP END PARALLEL DO
         call this%err%werrfl2(class//sname//' ended')
      end subroutine copyin
!
      subroutine copyout(this,fd2d,lpos,sdim,ddim) 
! copy 3d field from the local longitudinal position lpos to a 2d field
         implicit none
         
         class(ufield3d), intent(inout) :: this
         class(ufield2d), intent(in), target :: fd2d
         integer, intent(in) :: lpos
         integer, dimension(:), intent(in):: sdim, ddim
! local data         
         character(len=18), save :: sname = 'copyout:'
         real, dimension(:,:,:), pointer :: rf2d
         integer :: i,j,k,rank,nd1,nd2

         call this%err%werrfl2(class//sname//' started')
         
         rf2d => fd2d%getrf()
         rank = size(sdim)
         nd1 = size(this%rf,2)
         nd2 = size(this%rf,3)

!$OMP PARALLEL DO PRIVATE(i,j,k)
         do k = 1, nd2
            do j = 1, nd1
               do i = 1, rank
                  rf2d(ddim(i),j,k) = this%rf(sdim(i),j,k,lpos)
               enddo
            enddo
         enddo
!$OMP END PARALLEL DO
         call this%err%werrfl2(class//sname//' ended')
      end subroutine copyout
!
      subroutine writehdf5_3d(this,file,dim,rtag,stag,id)

         implicit none
         
         class(ufield3d), intent(inout) :: this
         class(hdf5file), intent(in) :: file
         integer, intent(in) :: dim, rtag, stag
         integer, intent(inout) :: id
! local data
         integer, dimension(3) :: gsize, lsize
         integer, dimension(2) :: noff
         integer :: ierr
         character(len=18), save :: sname = 'writehdf5_3d:'

         call this%err%werrfl2(class//sname//' started')
                
         gsize = (/this%nd1,this%nd2,this%nd3/)
         lsize = (/this%nd1p,this%nd2p,this%nd3p/)
         noff = this%noff       
         call pwfield_pipe(this%p,this%err,file,this%rf(dim,:,:,:),gsize,lsize,&
         &noff,rtag,stag,id,ierr)

         call this%err%werrfl2(class//sname//' ended')
      
      end subroutine writehdf5_3d
!      
      subroutine writehdf5_2dslice(this,file,dim,slice,spos,rtag,stag,id)
! spos = global position of the slice
         implicit none
         
         class(ufield3d), intent(inout) :: this
         class(hdf5file), intent(in) :: file
         integer, intent(in) :: dim, rtag, stag, slice, spos
         integer, intent(inout) :: id
! local data
         integer, dimension(2) :: gsize, lsize
         integer, dimension(2) :: noff
         integer :: ks, lpos, ierr
         character(len=18), save :: sname = 'writehdf5_2dslice:'

         call this%err%werrfl2(class//sname//' started')

         noff = this%noff       
            
         select case (slice)
         case (1)
            gsize = (/this%nd2,this%nd3/)
            if (this%p%getnstage() == 1) then
               lsize = (/this%nd2p,this%nd3p/)
               call pwfield_pipe(this%p,this%err,file,this%rf(dim,spos,:,:),&
               &gsize,lsize,noff,rtag,stag,id,ierr)
            else
               if (this%p%getstageid() == 0) then
                  lsize = (/this%nd2p,this%nd3p+1/)
                  call pwfield_pipe(this%p,this%err,file,this%rf(dim,spos,:,:),&
                  &gsize,lsize,noff,rtag,stag,id,ierr)
               else if (this%p%getstageid() /= (this%p%getnstage()-1)) then
                  lsize = (/this%nd2p,this%nd3p/)
                  call pwfield_pipe(this%p,this%err,file,this%rf(dim,spos,:,2:),&
                  &gsize,lsize,noff+(/0,1/),rtag,stag,id,ierr)
               else
                  lsize = (/this%nd2p,this%nd3p-1/)
                  call pwfield_pipe(this%p,this%err,file,this%rf(dim,spos,:,2:),&
                  &gsize,lsize,noff+(/0,1/),rtag,stag,id,ierr)
               end if
            end if
         case (2)
            ks = (spos-1)/this%nd2p
            if (this%p%getlkstrt() == ks + 1) then
               gsize = (/this%nd1,this%nd3/)
               if (this%p%getnstage() == 1) then
                  lsize = (/this%nd1p,this%nd3p/)
                  lpos = spos - ks*this%nd2p
                  call wfield_pipe(this%p,this%err,file,this%rf(dim,:,lpos,:),&
                  &gsize,lsize,noff,rtag,stag,id,ierr)
               else
                  if (this%p%getstageid() == 0) then
                     lsize = (/this%nd1p,this%nd3p+1/)
                     lpos = spos - ks*this%nd2p
                     call wfield_pipe(this%p,this%err,file,this%rf(dim,:,lpos,:),&
                     &gsize,lsize,noff,rtag,stag,id,ierr)
                  else if (this%p%getstageid() /= (this%p%getnstage()-1)) then
                     lsize = (/this%nd1p,this%nd3p/)
                     lpos = spos - ks*this%nd2p
                     call wfield_pipe(this%p,this%err,file,this%rf(dim,:,lpos,2:),&
                     &gsize,lsize,noff+(/0,1/),rtag,stag,id,ierr)               
                  else
                     lsize = (/this%nd1p,this%nd3p-1/)
                     lpos = spos - ks*this%nd2p
                     call wfield_pipe(this%p,this%err,file,this%rf(dim,:,lpos,2:),&
                     &gsize,lsize,noff+(/0,1/),rtag,stag,id,ierr)               
                  end if
               end if
            else         
               id = MPI_REQUEST_NULL         
            endif   
         case (3)
            ks = (spos-1)/this%nd3p
            if (this%p%getstageid() == ks) then
               gsize = (/this%nd1,this%nd2/)
               lsize = (/this%nd1p,this%nd2p/)
               lpos = spos - ks*this%nd3p
               call pwfield(this%p,this%err,file,this%rf(dim,:,:,lpos),&
               &gsize,lsize,noff(1),ierr)
            endif
            id = MPI_REQUEST_NULL         
         end select

         call this%err%werrfl2(class//sname//' ended')
      
      end subroutine writehdf5_2dslice
!
      subroutine sum(this,a1,a2)
      
         implicit none
         
         class(ufield3d),intent(inout) :: this
         class(ufield3d), target, intent(in) :: a1,a2
! local data
         character(len=18), save :: sname = 'sum:'
         integer :: i,j,k,l
         real, dimension(:,:,:,:), pointer :: rf1 => null(), rf2 => null(),&
         &rf3 => null()
         
         call this%err%werrfl2(class//sname//' started')
         
         rf1 => this%rf         
         rf2 => a1%rf
         rf3 => a2%rf
         

!$OMP PARALLEL DO PRIVATE(i,j,k,l)         
         do l = 1, size(rf1,4)  
            do k = 1, size(rf1,3)
               do j = 1, size(rf1,2)
                  do i = 1, size(rf1,1)             
                     rf1(i,j,k,l) = rf2(i,j,k,l) + rf3(i,j,k,l)
                  enddo
               enddo
            enddo
         enddo
!$OMP END PARALLEL DO    

         call this%err%werrfl2(class//sname//' ended')

      end subroutine sum
!
      subroutine minus(this,a1,a2)
      
         implicit none
         
         class(ufield3d),intent(inout) :: this
         class(ufield3d), target, intent(in) :: a1,a2
! local data
         character(len=18), save :: sname = 'minus:'
         integer :: i,j,k,l
         real, dimension(:,:,:,:), pointer :: rf1 => null(), rf2 => null(),&
         &rf3 => null()
         
         call this%err%werrfl2(class//sname//' started')
         
         rf1 => this%rf         
         rf2 => a1%rf
         rf3 => a2%rf
         

!$OMP PARALLEL DO PRIVATE(i,j,k)         
         do l = 1, size(rf1,4)  
            do k = 1, size(rf1,3)
               do j = 1, size(rf1,2)
                  do i = 1, size(rf1,1)             
                     rf1(i,j,k,l) = rf2(i,j,k,l) - rf3(i,j,k,l)
                  enddo
               enddo
            enddo
         enddo
!$OMP END PARALLEL DO    

         call this%err%werrfl2(class//sname//' ended')

      end subroutine minus
!
      subroutine multiply(this,a1,value)
      
         implicit none
         
         class(ufield3d), intent(inout) :: this
         class(ufield3d), target, intent(in) :: a1
         real, intent(in) :: value
! local data
         character(len=18), save :: sname = 'multiply:'
         integer :: i,j,k,l
         real, dimension(:,:,:,:), pointer :: rf1 => null(), rf2 => null()
         
         call this%err%werrfl2(class//sname//' started')

         rf1 => this%rf         
         rf2 => a1%rf

!$OMP PARALLEL DO PRIVATE(i,j,k)         
         do l = 1, size(rf1,4)
            do k = 1, size(rf1,3)
               do j = 1, size(rf1,2)
                  do i = 1, size(rf1,1)               
                     rf1(i,j,k,l) = rf2(i,j,k,l) * value
                  enddo
               enddo
            enddo
         enddo
!$OMP END PARALLEL DO        
         call this%err%werrfl2(class//sname//' ended')

      end subroutine multiply
!
      subroutine asc(this,value)
      
         implicit none
         
         class(ufield3d), intent(inout) :: this
         real, intent(in) :: value
! local data
         character(len=18), save :: sname = 'asc:'
         integer :: i,j,k,l
         real, dimension(:,:,:,:), pointer :: rf => null()
         
         call this%err%werrfl2(class//sname//' started')
         
         rf => this%rf

!$OMP PARALLEL DO PRIVATE(i,j,k)         
         do l = 1, size(rf,4)
            do k = 1, size(rf,3)
               do j = 1, size(rf,2)
                  do i = 1, size(rf,1)               
                     rf(i,j,k,l) = value
                  enddo
               enddo
            enddo
         enddo
!$OMP END PARALLEL DO
                    
         call this%err%werrfl2(class//sname//' ended')

      end subroutine asc
!
      subroutine asa(this,that)
      
         implicit none
         
         class(ufield3d), intent(inout) :: this
         class(ufield3d), target, intent(in) :: that
! local data
         character(len=18), save :: sname = 'asa:'
         integer :: i,j,k,l
         real, dimension(:,:,:,:), pointer :: rf1 => null(), rf2 => null()
         
         call this%err%werrfl2(class//sname//' started')
         
         rf1 => this%rf
         rf2 => that%rf

!$OMP PARALLEL DO PRIVATE(i,j,k)         
         do l = 1, size(rf1,4) 
            do k = 1, size(rf1,3)
               do j = 1, size(rf1,2)
                  do i = 1, size(rf1,1)              
                     rf1(i,j,k,l) = rf2(i,j,k,l)
                  enddo
               enddo
            enddo
         enddo
!$OMP END PARALLEL DO
!         deallocate(rf2)

         call this%err%werrfl2(class//sname//' ended')

      end subroutine asa
!      
      end module ufield3d_class