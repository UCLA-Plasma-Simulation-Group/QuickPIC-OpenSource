! Laser class for QuickPIC Open Source 1.0
! update: 08/10/2016

      module laser_class

      use perrors_class
      use parallel_pipe_class
      use spect3d_class
      use fdist3d_class
      use field3d_class
      use field2d_class
      use ufield3d_class
      use ufield2d_class
      use hdf5io_class
      use mpi
               
      implicit none

      private

      public :: laser

      type laser

         private

         class(spect3d), pointer, public :: sp => null()
         class(perrors), pointer, public :: err => null()
         class(parallel_pipe), pointer, public :: p => null()
         class(field3d), pointer :: lp => null()
         class(field3d), pointer :: a0x => null(),a0y => null(),a0z => null()
         class(field3d), pointer :: a1x => null(),a1y => null(),a1z => null()
         class(field3d), pointer :: a2x => null(),a2y => null(),a2z => null()
         class(fdist3d), pointer :: pf => null()
         class(field2d), pointer :: a0r => null(),a0i => null()
         real :: k0, ds, dz, dx
         contains
         
         generic :: new => init_laser
         generic :: del => end_laser
         generic :: push => push_laser
         generic :: wrlp => writelp_laser, writelpslice_laser
         generic :: wrst => writerst_laser       
         generic :: rrst => readrst_laser       
         procedure :: getlp
         procedure, private :: init_laser
         procedure, private :: end_laser
         procedure, private :: push_laser
         procedure, private :: writerst_laser, readrst_laser
         procedure, private :: writelp_laser, writelpslice_laser 
      end type 

      save      

      character(len=10) :: class = 'laser:'
      character(len=128) :: erstr
      
      contains
!
      subroutine init_laser(this,pp,perr,psp,pf,k0,ds,dz,dx)
      
         implicit none
         
         class(beam3d), intent(inout) :: this
         class(spect3d), intent(in), pointer :: psp
         class(perrors), intent(in), pointer :: perr
         class(parallel_pipe), intent(in), pointer :: pp
         class(fdist3d), intent(inout), target :: pf
         real, intent(in) :: k0, ds, dx

! local data
         character(len=18), save :: sname = 'init_laser:'
         integer :: id, ierr
         integer, dimension(10) :: istat
                  
         this%sp => psp
         this%err => perr
         this%p => pp
         this%pf => pf

         call this%err%werrfl2(class//sname//' started')

         this%k0 = k0
         this%ds = ds
         this%dz = dz
         this%dx = dx
         allocate(this%lp,this%a0x,this%a0y,this%a0z)
         allocate(this%a1x,this%a1y,this%a1z)
         allocate(this%a2x,this%a2y,this%a2z)
         allocate(this%a0r,this%a0i)
         call this%lp%new(this%p,this%err,this%sp,dim=1)
         call this%a0x%new(this%p,this%err,this%sp,dim=2)
         call this%a0y%new(this%p,this%err,this%sp,dim=2)
         call this%a0z%new(this%p,this%err,this%sp,dim=2)
         call this%a1x%new(this%p,this%err,this%sp,dim=2)
         call this%a1y%new(this%p,this%err,this%sp,dim=2)
         call this%a1z%new(this%p,this%err,this%sp,dim=2)
         call this%a2x%new(this%p,this%err,this%sp,dim=2)
         call this%a2y%new(this%p,this%err,this%sp,dim=2)
         call this%a2z%new(this%p,this%err,this%sp,dim=2)
         call this%a0r%new(this%p,this%err,this%sp,dim=1,fftflag=.true.)
         call this%a0i%new(this%p,this%err,this%sp,dim=1,fftflag=.true.)
         call this%pf%dist%(this%a0x,this%a0y,this%a0z)
         
         call this%err%werrfl2(class//sname//' ended')

      end subroutine init_laser
!
      subroutine end_beam3d(this)
          
         implicit none
         
         class(beam3d), intent(inout) :: this
         character(len=18), save :: sname = 'end_beam3d:'

         call this%err%werrfl2(class//sname//' started')
         call this%lp%del()
         call this%a0x%del()
         call this%a0y%del()
         call this%a0z%del()
         call this%a1x%del()
         call this%a1y%del()
         call this%a1z%del()
         call this%a2x%del()
         call this%a2y%del()
         call this%a2z%del()
         call this%err%werrfl2(class//sname//' ended')
                  
      end subroutine end_beam3d
!      
      subroutine push_laser(this,chi)
      
         implicit none
         
         class(beam3d), intent(inout) :: this
         class(field3d), intent(in) :: chi
! local data
         character(len=18), save :: sname = 'push_laser:'

         type(ufield3d), pointer, :: f3d => null()
         type(ufield2d), pointer, :: f2d => null()
         real, pointer, dimension(:,:,:,:) :: a0 => null()
         real, pointer, dimension(:,:,:,:) :: a1 => null()
         real, pointer, dimension(:,:,:,:) :: a2 => null()
         real, pointer, dimension(:,:) :: a0r => null(), a0i => null()
         real :: k0, ds, dz, dx
         real, pointer, dimension(:,:,:,:) :: vp => null(), pchi => null()
         complex, dimension(:), allocatable :: cp,dp
         real :: idx2
         integer :: i, ierr, j, k, mm, ii, jj
         integer :: nx, ny, nz
         ! integer :: ntasks, nxe, kyp
         ! integer :: nx, ny, nxh, nz
         ! integer, dimension(3) :: nd_uniform, ndg_uniform
         ! real :: zero = 0.0, affp,  we, wke
         real :: dzds, ks, k2, kk
         complex :: c1, c2, c3, d1, d2, d3, cc
        
         call this%err%werrfl2(class//sname//' started')
   
         k0 = this%k0
         ds = this%ds
         dz = this%dz
         dx = this%dx
         ci = cmplx(0.0,1.0)
         idx2 = 1.0/dx/dx         
         dzds = 1.0/dz/ds
         k2 = k0*k0*0.5
         ks = -2.0*k0/ds
         f3d => this%a0x%getrs()
         nx = f3d%getnd1p()
         ny = f3d%getnd2p()
         nz = f3d%getnd3p()
         a0 => f3d%getrf()
         f3d => chi%getrs()
         pchi => f3d%getrf()
         f3d => this%vp%getrs()
         vp => f3d%getrf()
         f2d => this%a0r%getks()
         a0r => f2d%getrf()
         f2d => this%a0i%getks()
         a0i => f2d%getrf()
         allocate(cp(nz),dp(nz))
         
         do i = 1, nx
            do j = 1, ny
               cp(1) = cmplx(dzds,0.0)/cmplx(-k2*chi(1,i,j,1),ks)
               dp(1) = (cmplx(k2*chi(1,i,j,1),ks)*cmplx(a0(1,i,j,1),a0(2,i,j,1))+&
               &cmplx(dzds,0.0)*cmplx(a0(1,i,j,2),a0(2,i,j,2)))/cmplx(-k2*chi(1,i,j,1),ks)
               do k = 2, nz-1
                  cp(k) = cmplx(dzds,0.0)/(cmplx(-k2*chi(1,i,j,k),ks)+cmplx(dzds,0.0)*cp(k-1))
                  dp(k) = (cmplx(-dzds,0.0)*cmplx(a0(1,i,j,k-1),a0(2,i,j,k-1))+&
                  &cmplx(k2*chi(1,i,j,k),ks)*cmplx(a0(1,i,j,k),a0(2,i,j,k))+&
                  &cmplx(dzds,0.0)*cmplx(a0(1,i,j,k+1),a0(2,i,j,k+1))+cmplx(dzds,0.0)*dp(k-1))/&
                  &(cmplx(-k2*chi(1,i,j,k),ks)+cmplx(dzds,0.0)*cp(k-1))
               enddo
               cp(nz) = cmplx(dzds,0.0)/(cmplx(-k2*chi(1,i,j,nz),ks)+cmplx(dzds,0.0)*cp(nz-1))
               dp(nz) = (cmplx(-dzds,0.0)*cmplx(a0(1,i,j,nz-1),a0(2,i,j,nz-1))+&
               &cmplx(k2*chi(1,i,j,nz),ks)*cmplx(a0(1,i,j,nz),a0(2,i,j,nz))+&
               &cmplx(dzds,0.0)*dp(nz-1))/&
               &(cmplx(-k2*chi(1,i,j,nz),ks)+cmplx(dzds,0.0)*cp(nz-1))
               a2(1,i,j,nz) = real(dp(nz))
               a2(2,i,j,nz) = aimag(dp(nz))
               do k = nz-1,1,-1
                  c1 = dp(k) - cp(k)*cmplx(a2(1,i,j,k+1),a2(2,i,j,k+1))
                  a2(1,i,j,k) = real(c1)
                  a2(2,i,j,k) = aimag(c1)
               enddo      
            enddo
         enddo
   
         a0(:,:,:,:) = a1(:,:,:,:)
         
         do k = 1, nz      
            do i = 1, nx
               do j = 1, ny
                  a0r(i,j) = a2(1,i,j,k)
                  a0i(i,j) = a2(2,i,j,k)      
               enddo
            enddo
            call fftrk(a0r)
            call fftrk(a0i)
   
         do i = 1, nd_uniform(1)+1
         do j = 1, nd_uniform(2)+1
         
         a2(i,j,k,1) = cmplx(a0r%sfd(i,j,1),a0i%sfd(i,j,1))
   
         enddo
         enddo
         enddo
         
         do i = 1, nd_uniform(1)+1
         do j = 1, nd_uniform(2)+1
         
         if((i==1) .or. (i==(nd_uniform(1)+1)) .or. (j==(nd_uniform(2)+1))) then
            kk = cmplx(0.0,0.0)
         else if ((j==1) .and. (id0==0)) then
            kk = cmplx(0.0,0.0)
         else
            kk = aimag(a0r%pt%ffc(i,j,1))/real(a0r%pt%ffc(i,j,1))/idx2*0.25
         endif
   
   !      ace = 2.25*a0(i,j,nz,1)-1.5*a0(i,j,nz-1,1)+0.25*a0(i,j,nz-2,1)
         
         cp(1) = cmplx(dzds,0.0)/cmplx(kk,ks)
         dp(1) = (cmplx(-kk,ks)*a2(i,j,1,1)+cmplx(dzds,0.0)*a2(i,j,2,1))&
         &/cmplx(kk,ks)
         
         do k = 2, nz-1
         
         cp(k) = cmplx(dzds,0.0)/(cmplx(kk,ks)+cmplx(dzds,0.0)*cp(k-1))
         dp(k) = (cmplx(-dzds,0.0)*a2(i,j,k-1,1)+cmplx(-kk,ks)*a2(i,j,k,1)+&
         &cmplx(dzds,0.0)*a2(i,j,k+1,1)+cmplx(dzds,0.0)*dp(k-1))/&
         &(cmplx(kk,ks)+cmplx(dzds,0.0)*cp(k-1))
   
         enddo
   
         cp(nz) = cmplx(dzds,0.0)/(cmplx(kk,ks)+cmplx(dzds,0.0)*cp(nz-1))
         dp(nz) = (cmplx(-dzds,0.0)*a2(i,j,nz-1,1)+cmplx(-kk,ks)*a2(i,j,nz,1)+&
         &cmplx(dzds,0.0)*dp(nz-1))/&
         &(cmplx(kk,ks)+cmplx(dzds,0.0)*cp(nz-1))
   
         a1(i,j,nz,1) = dp(nz)
   
         do k = nz-1,1,-1
   
         a1(i,j,k,1) = dp(k) - cp(k)*a1(i,j,k+1,1)
         
         enddo      
         enddo
         enddo
              
        do k = 1, nz      
         do i = 1, ndg_uniform(1)-1
         do j = 1, ndg_uniform(2)
         a0r%sfd(i,j,1) = real(a1(i,j,k,1))
         a0i%sfd(i,j,1) = aimag(a1(i,j,k,1))      
         enddo
         enddo
         
         call smoothf(a0r)
   
        a0r%state = 1 
        a0r%gcells = 0      
   
        call fftkr(a0r,uniform=.true.)
         
        a0i%state = 1 
        a0i%gcells = 0      
   
         call smoothf(a0i)
   
        call fftkr(a0i,uniform=.true.)
   
         do i = 1, nd_uniform(1)
         do j = 1, nd_uniform(2)
         
         a1(i,j,k,1) = cmplx(a0r%sf(i,j,1),a0i%sf(i,j,1))
   
         enddo
         enddo
         enddo
   
         vp=real(a1*conjg(a1))
   
         deallocate(cp,dp)
         
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