! part3d class for QuickPIC Open Source 1.0
! update: 04/18/2016

      module part3d_class

      use perrors_class
      use parallel_pipe_class
      use spect3d_class
      use fdist3d_class
      use ufield3d_class
      use part3d_lib
      use hdf5io_class
      use mpi
               
      implicit none

      private

      public :: part3d

      type part3d

         private

         class(spect3d), pointer, public :: sp => null()
         class(perrors), pointer, public :: err => null()
         class(parallel_pipe), pointer, public :: p => null()
!
! qm = charge on particle, in units of e
! qbm = particle charge/mass ratio
! dt = time interval between successive calculations
! ci = reciprical of velocity of light
! xdim = dimension of the particle coordinates
! nbmax = size of buffer for passing particles between processors
! npp = number of particles in current partition
! npmax = maximum number of particles in each partition
! part(:,:) = initial particle coordinates
!         
         real :: qm, qbm, dt, ci
         integer :: npmax, nbmax, xdim, npp = 0
         real, dimension(:,:), pointer :: part => null(), pbuff => null()
         
         contains
         
         generic :: new => init_part3d
         generic :: del => end_part3d
         generic :: push => partpush
         generic :: pmv => pmove
         generic :: qdp => qdeposit  
         generic :: wr => writehdf5_part3d
         generic :: wrst => writerst_part3d
         generic :: rrst => readrst_part3d
         procedure, private :: init_part3d
         procedure, private :: end_part3d
         procedure, private :: partpush
         procedure, private :: pmove         
         procedure, private :: qdeposit, writehdf5_part3d
         procedure, private :: writerst_part3d, readrst_part3d
         procedure :: getnpp
                  
      end type 

      save      

      character(len=10) :: class = 'part3d:'
      character(len=128) :: erstr
!
! buffer data for particle managers
      real, dimension(:,:), allocatable :: sbufl, sbufr, rbufl, rbufr
      integer, dimension(:), allocatable :: ihole
      integer :: szbuf = 0
      
      
      contains
!
      subroutine init_part3d(this,pp,perr,psp,pf,fd,qm,qbm,dt,ci,xdim,npmax,nbmax)
      
         implicit none
         
         class(part3d), intent(inout) :: this
         class(spect3d), intent(in), pointer :: psp
         class(perrors), intent(in), pointer :: perr
         class(parallel_pipe), intent(in), pointer :: pp
         class(fdist3d), intent(in) :: pf
         class(ufield3d), target, intent(in) :: fd
         real, intent(in) :: qm, qbm, dt, ci
         integer, intent(in) :: npmax, nbmax, xdim

! local data
         character(len=18), save :: sname = 'init_part3d:'
         integer :: noff, nxyp, nx, prof
                  
         this%sp => psp
         this%err => perr
         this%p => pp

         call this%err%werrfl2(class//sname//' started')

         this%qm = qm
         this%qbm = qbm
         this%dt = dt
         this%ci = ci
         this%xdim = xdim
         this%npmax = npmax
         this%nbmax = nbmax
         prof = pf%getnpf()

         allocate(this%part(xdim,npmax),this%pbuff(xdim,nbmax))
         
                           
         select case (prof)
         case (1)
            call init_prof1(this,pf,fd)
         end select

         call this%err%werrfl2(class//sname//' ended')

      end subroutine init_part3d
!
      subroutine end_part3d(this)
          
         implicit none
         
         class(part3d), intent(inout) :: this
         character(len=18), save :: sname = 'end_part3d:'

         call this%err%werrfl2(class//sname//' started')
         deallocate(this%part,this%pbuff)
         call this%err%werrfl2(class//sname//' ended')
         
         return
         
      end subroutine end_part3d
!
      subroutine init_prof1(this,pf,fd)
      
         implicit none
         
         class(part3d), intent(inout) :: this
         class(fdist3d), intent(in) :: pf
         class(ufield3d), intent(in) :: fd
! local data1

! edges(1) = lower boundary in y of particle partition
! edges(2) = upper boundary in y of particle partition
! edges(3) = lower boundary in z of particle partition
! edges(4) = upper boundary in z of particle partition
         real, dimension(:,:), pointer :: pt => null()
         real, dimension(3,100) :: arg
         integer :: npx, npy, npz, nx, ny, nz, ipbc
         real :: vtx, vty, vtz, vdx, vdy, vdz
         real :: sigx, sigy, sigz, x0, y0, z0
         real, dimension(3) :: cx, cy
         real, dimension(4) :: edges
         integer, dimension(2) :: noff
         integer :: nps=1
         logical :: lquiet = .false.
         integer :: idimp, npmax, ierr = 0
         character(len=18), save :: sname = 'init_prof1:'

         call this%err%werrfl2(class//sname//' started')
         
         npx = pf%getnpx(); npy = pf%getnpy(); npz = pf%getnpz()
         nx = fd%getnd1(); ny = fd%getnd2(); nz = fd%getnd3()
         ipbc = this%sp%getpsolver()
         arg = pf%getarg()
         pt => this%part
         vtx = arg(1,1); vty = arg(2,1); vtz = arg(3,1)
         vdx = arg(1,2); vdy = arg(2,2); vdz = arg(3,2)
         sigx = arg(1,3); sigy = arg(2,3); sigz = arg(3,3)
         x0 = arg(1,4); y0 = arg(2,4); z0 = arg(3,4)
         cx = arg(1,5:7); cy = arg(2,5:7)
         if (arg(1,8) >= 0) lquiet = .true.
         idimp = this%xdim; npmax = this%npmax
         noff = fd%getnoff()
         edges(1) = noff(1); edges(3) = noff(2)
         edges(2) = edges(1) + fd%getnd2p()
         edges(4) = edges(3) + fd%getnd3p()         
         
         call PRVDIST32_RANDOM(pt,this%qm,edges,this%npp,nps,vtx,vty,vtz,vdx,vdy,&
         &vdz,npx,npy,npz,nx,ny,nz,ipbc,idimp,npmax,1,1,4,sigx,sigy,sigz,&
         &x0,y0,z0,cx,cy,lquiet,ierr)

         if (ierr /= 0) then
            write (erstr,*) 'PRVDIST32_RANDOM error'
            call this%err%equit(class//sname//erstr)
         endif
         
         
         call this%err%werrfl2(class//sname//' ended')
         
      end subroutine init_prof1
!      
      subroutine qdeposit(this,q)
! deposit the charge density      
      
         implicit none
         
         class(part3d), intent(in) :: this
         class(ufield3d), target, intent(inout) :: q
         character(len=18), save :: sname = 'qdeposit:'
! local data
         real, dimension(:,:,:,:), pointer :: pq => null()
         real :: qm
         real, dimension(:,:), pointer :: part
         integer :: npp
         integer, dimension(2) :: noff
         integer :: idimp, npmax, nxv, nypmx, nzpmx, nxyzp
         integer :: order, opt
                  
         call this%err%werrfl2(class//sname//' started')
         
         pq => q%getrf()
         part => this%part
         noff = q%getnoff()
         qm = this%qm
         npp = this%npp

         idimp = this%xdim; npmax = this%npmax
         nxv = size(pq,2); nypmx = size(pq,3); nzpmx = size(pq,4)
         nxyzp = nxv*nypmx*nzpmx
         
         select case (this%sp%getinorder())
         case (1)
            call PGPOST32L(part,pq(1,:,:,:),npp,noff,idimp,npmax,1,nxv,nypmx,&
            &nzpmx,2)
         case default
            call PGPOST32L(part,pq(1,:,:,:),npp,noff,idimp,npmax,1,nxv,nypmx,&
            &nzpmx,2)         
         end select
         call this%err%werrfl2(class//sname//' ended')
         
      end subroutine qdeposit
!      
      subroutine partpush(this,ef,bf,dex,dez)
      
         implicit none
         
         class(part3d), intent(inout) :: this
         class(ufield3d), target, intent(in) :: ef, bf
         real, intent(in) :: dex, dez
         character(len=18), save :: sname = 'partpush'
! local data
         real, dimension(:,:,:,:), pointer :: pef => null(), pbf => null()
         integer :: nx, ny, nz, ipbc
         real :: qbm, dt, dtc, ek
         integer, dimension(2) :: noff
         integer :: idimp, npmax, nxv, nypmx, nzpmx, nxyzp

         call this%err%werrfl2(class//sname//' started')

         pef => ef%getrf(); pbf => bf%getrf()
         qbm = this%qbm; dt = this%dt
         nx = ef%getnd1(); ny = ef%getnd2(); nz = ef%getnd3()
         idimp = this%xdim; npmax = this%npmax;
         ipbc = this%sp%getpsolver()
         noff = ef%getnoff()
         nxv = size(pef,2); nypmx = size(pef,3); nzpmx = size(pef,4)
         nxyzp = nxv*nypmx*nzpmx
         
         select case (this%sp%getinorder())
         case (1)
            call PGBPUSH32L_QP(this%part,pef,pbf,this%npp,noff,qbm,dt,dt,ek,&
            &nx,ny,nz,idimp,npmax,1,nxv,nypmx,nzpmx,2,ipbc,dex,dez,0.0)
         case default
            call PGBPUSH32L_QP(this%part,pef,pbf,this%npp,noff,qbm,dt,dt,ek,&
            &nx,ny,nz,idimp,npmax,1,nxv,nypmx,nzpmx,2,ipbc,dex,dez,0.0)
         end select

         call this%err%werrfl2(class//sname//' ended')
         
      end subroutine partpush
!
      subroutine pmove(this,fd,rtag,stag,sid)
      
         implicit none
         
         class(part3d), intent(inout) :: this
         class(ufield3d), target, intent(in) :: fd
         integer, intent(in) :: rtag, stag
         integer, intent(inout) :: sid
         character(len=18), save :: sname = 'pmove:'
! local data
         integer :: ny, nz, nvpy, nvpz, nbmax, idds = 2
         integer :: ierr
         real, dimension(4) :: edges
         integer, dimension(2) :: noff         
         integer, dimension(2) :: jsl, jsr, jss
         integer, dimension(9) :: info
         integer :: idimp, npmax, idps, ntmax
         real, dimension(:,:), pointer :: pbuff
         
         call this%err%werrfl2(class//sname//' started')
         
         idimp = this%xdim; npmax = this%npmax; nbmax = this%nbmax
         idps = size(edges,1)
         ntmax = 2*nbmax
         ny = fd%getnd2(); nz = fd%getnd3()
         nvpy = this%p%getlnvp()
         nvpz = this%p%getnstage()
         noff = fd%getnoff()         
         edges(1) = noff(1); edges(3) = noff(2)
         edges(2) = edges(1) + fd%getnd2p()
         edges(4) = edges(3) + fd%getnd3p()
         pbuff => this%pbuff         
         
! check if size of buffers has changed
         if (szbuf.ne.nbmax) then
            if (szbuf /= 0) deallocate(sbufl,sbufr,rbufl,rbufr,ihole)
! allocate buffers
            allocate(sbufl(idimp,nbmax))
            allocate(sbufr(idimp,nbmax))
            allocate(rbufl(idimp,nbmax))
            allocate(rbufr(idimp,nbmax))
            allocate(ihole(ntmax))
            szbuf = nbmax
         endif

         call PMOVE32(this%part,edges,this%npp,sbufr,sbufl,rbufr,rbufl,ihole,pbuff, &
         &jsr,jsl,jss,ny,nz,this%p%getkstrt(),nvpy,nvpz,idimp,npmax,1,1,idps,nbmax,&
         &idds,ntmax,rtag,stag,sid,info)

         if (info(1) /= 0) then
            write (erstr,*) 'PMOVE32 error'
            call this%err%equit(class//sname//erstr)
         endif
         
         if (this%p%getstageid() == this%p%getnstage() - 1) sid = MPI_REQUEST_NULL         

         call this%err%werrfl2(class//sname//' ended')
         
      end subroutine pmove
!      
      function getnpp(this)

         implicit none

         class(part3d), intent(in) :: this
         integer :: getnpp
         
         getnpp = this%npp

      end function getnpp
!
      subroutine writehdf5_part3d(this,file,dspl,delta,rtag,stag,id)

         implicit none
         
         class(part3d), intent(inout) :: this
         class(hdf5file), intent(in) :: file
         real, dimension(3), intent(in) :: delta
         integer, intent(in) :: dspl, rtag, stag
         integer, intent(inout) :: id
! local data
         character(len=18), save :: sname = 'writehdf5_part3d:'
         integer :: ierr = 0

         call this%err%werrfl2(class//sname//' started')                  
         call pwpart_pipe(this%p,this%err,file,this%part,this%npp,dspl,delta,&
         &rtag,stag,id,ierr)
         call this%err%werrfl2(class//sname//' ended')
      
      end subroutine writehdf5_part3d
!            
      subroutine writerst_part3d(this,file)

         implicit none
         
         class(part3d), intent(inout) :: this
         class(hdf5file), intent(in) :: file
! local data
         character(len=18), save :: sname = 'writerst_part3d:'
         integer :: ierr = 0

         call this%err%werrfl2(class//sname//' started')                  
         call wpart(this%p,this%err,file,this%part,this%npp,1,ierr)
         call this%err%werrfl2(class//sname//' ended')
      
      end subroutine writerst_part3d
!            
      subroutine readrst_part3d(this,file)

         implicit none
         
         class(part3d), intent(inout) :: this
         class(hdf5file), intent(in) :: file
! local data
         character(len=18), save :: sname = 'readrst_part3d:'
         integer :: ierr = 0

         call this%err%werrfl2(class//sname//' started')                  
         call rpart(this%p,this%err,file,this%part,this%npp,ierr)
         call this%err%werrfl2(class//sname//' ended')
      
      end subroutine readrst_part3d
!            
      end module part3d_class