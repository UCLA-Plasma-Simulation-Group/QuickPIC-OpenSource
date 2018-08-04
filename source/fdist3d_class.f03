! fdist3d class for QuickPIC Open Source 1.0
! update: 04/18/2016

      module fdist3d_class

      use perrors_class
      use parallel_pipe_class
      use spect3d_class
      use ufield3d_class
      use part3d_lib
      use input_class
         
      implicit none

      private

      public :: fdist3d, fdist3d_000

      type, abstract :: fdist3d

         private

         class(spect3d), pointer, public :: sp => null()
         class(perrors), pointer, public :: err => null()
         class(parallel_pipe), pointer, public :: p => null()
!
! ndprof = profile type
         integer :: npf, npmax
                         
         contains
         
         generic :: new => init_fdist3d         
         generic :: del => end_fdist3d
         generic :: dist => dist3d
         procedure(ab_init_fdist3d), deferred, private :: init_fdist3d
         procedure, private :: end_fdist3d
         procedure(ab_dist3d), deferred, private :: dist3d
         procedure :: getnpf, getnpmax
                  
      end type 

      abstract interface
!
      subroutine ab_dist3d(this,part3d,npp,fd)
         import fdist3d
         import ufield3d
         implicit none
         class(fdist3d), intent(inout) :: this
         real, dimension(:,:), pointer, intent(inout) :: part3d
         integer, intent(inout) :: npp
         class(ufield3d), intent(in), pointer :: fd         
      end subroutine ab_dist3d
!
      subroutine ab_init_fdist3d(this,input,i)
         import fdist3d
         import input_json
         implicit none
         class(fdist3d), intent(inout) :: this
         type(input_json), intent(inout), pointer :: input
         integer, intent(in) :: i
      end subroutine ab_init_fdist3d
!
      end interface

      type, extends(fdist3d) :: fdist3d_000

         private

         integer :: npx, npy, npz
         real :: qm, sigx, sigy, sigz
         real :: bcx, bcy, bcz, sigvx, sigvy, sigvz
         real :: cx1,cx2,cx3,cy1,cy2,cy3,gamma,np
         logical :: quiet

         contains
         procedure, private :: init_fdist3d => init_fdist3d_000
         procedure, private :: dist3d => dist3d_000

      end type fdist3d_000

      character(len=10), save :: class = 'fdist3d:'
      character(len=128), save :: erstr
      
      contains
!
      function getnpf(this)

         implicit none

         class(fdist3d), intent(in) :: this
         integer :: getnpf
         
         getnpf = this%npf

      end function getnpf
!      
      function getnpmax(this)

         implicit none

         class(fdist3d), intent(in) :: this
         integer :: getnpmax
         
         getnpmax = this%npmax

      end function getnpmax
!      
      subroutine end_fdist3d(this)
          
         implicit none
         
         class(fdist3d), intent(inout) :: this
         character(len=18), save :: sname = 'end_fdist3d:'

         call this%err%werrfl2(class//sname//' started')
         call this%err%werrfl2(class//sname//' ended')
                  
      end subroutine end_fdist3d
!
      subroutine init_fdist3d_000(this,input,i)
      
         implicit none
         
         class(fdist3d_000), intent(inout) :: this
         type(input_json), intent(inout), pointer :: input
         integer, intent(in) :: i        
! local data
         integer :: npf,npx,npy,npz,npmax
         real :: qm,sigx,sigy,sigz,bcx,bcy,bcz,sigvx,sigvy,sigvz
         real :: cx1,cx2,cx3,cy1,cy2,cy3,gamma,np
         logical :: quiet
         real :: min, max, cwp, n0
         real :: alx, aly, alz, dx, dy, dz
         integer :: indx, indy, indz         
         character(len=20) :: sn,s1
         character(len=18), save :: sname = 'init_fdist3d_000:'
         
         this%sp => input%sp
         this%err => input%err
         this%p => input%pp

         call this%err%werrfl2(class//sname//' started')

         write (sn,'(I3.3)') i
         s1 = 'beam('//trim(sn)//')'

         call input%get('simulation.n0',n0)
         call input%get('simulation.indx',indx)
         call input%get('simulation.indy',indy)
         call input%get('simulation.indz',indz)

         cwp=5.32150254*1e9/sqrt(n0)
         call input%get('simulation.box.x(1)',min)
         call input%get('simulation.box.x(2)',max)
         call input%get(trim(s1)//'.center(1)',bcx)
         bcx = bcx - min
         alx = (max-min)/cwp 
         dx=alx/real(2**indx)
         call input%get('simulation.box.y(1)',min)
         call input%get('simulation.box.y(2)',max)
         call input%get(trim(s1)//'.center(2)',bcy)
         bcy = bcy -min
         aly = (max-min)/cwp 
         dy=aly/real(2**indy)
         call input%get('simulation.box.z(1)',min)
         call input%get('simulation.box.z(2)',max)
         call input%get(trim(s1)//'.center(3)',bcz)
         bcz = bcz -min
         alz = (max-min)/cwp 
         dz=alz/real(2**indz)

         call input%get(trim(s1)//'.profile',npf)
         call input%get(trim(s1)//'.np(1)',npx)
         call input%get(trim(s1)//'.np(2)',npy)
         call input%get(trim(s1)//'.np(3)',npz)
         call input%get(trim(s1)//'.q',qm)
         call input%get(trim(s1)//'.sigma(1)',sigx)
         call input%get(trim(s1)//'.sigma(2)',sigy)
         call input%get(trim(s1)//'.sigma(3)',sigz)
         call input%get(trim(s1)//'.sigma_v(1)',sigvx)
         call input%get(trim(s1)//'.sigma_v(2)',sigvy)
         call input%get(trim(s1)//'.sigma_v(3)',sigvz)
         call input%get(trim(s1)//'.centroid_x(1)',cx1)
         call input%get(trim(s1)//'.centroid_x(2)',cx2)
         call input%get(trim(s1)//'.centroid_x(3)',cx3)
         call input%get(trim(s1)//'.centroid_y(1)',cy1)
         call input%get(trim(s1)//'.centroid_y(2)',cy2)
         call input%get(trim(s1)//'.centroid_y(3)',cy3)
         call input%get(trim(s1)//'.quiet_start',quiet)
         call input%get(trim(s1)//'.gamma',gamma)
         call input%get(trim(s1)//'.num_particle',np)
         call input%get(trim(s1)//'.npmax',npmax)


         this%npf = npf
         this%npx = npx
         this%npy = npy
         this%npz = npz
         this%npmax = npmax
         this%bcx = bcx/dx/cwp
         this%bcy = bcy/dy/cwp
         this%bcz = bcz/dz/cwp
         this%sigx = sigx/dx/cwp
         this%sigy = sigy/dy/cwp
         this%sigz = sigz/dz/cwp
         this%sigvx = sigvx
         this%sigvy = sigvy
         this%sigvz = sigvz
         this%cx1 = cx1*cwp*dz*dz/dx
         this%cx2 = cx2*dz/dx
         this%cx3 = cx3/cwp/dx
         this%cy1 = cy1*cwp*dz*dz/dy
         this%cy2 = cy2*dz/dy
         this%cy3 = cy3/cwp/dy
         this%gamma = gamma
         this%np = np
         this%quiet = quiet

         qm = qm*np*1e12*(2**indz)
         qm = qm*(2**indx)
         qm = qm*(2**indy)/(npx*alx*aly*alz*n0*cwp*cwp*cwp) 
         qm = qm/npy
         qm = qm/npz
         this%qm = qm

         call this%err%werrfl2(class//sname//' ended')

      end subroutine init_fdist3d_000
!
      subroutine dist3d_000(this,part3d,npp,fd)
      
         implicit none
         
         class(fdist3d_000), intent(inout) :: this
         real, dimension(:,:), pointer, intent(inout) :: part3d
         integer, intent(inout) :: npp
         class(ufield3d), intent(in) :: fd
! local data1

! edges(1) = lower boundary in y of particle partition
! edges(2) = upper boundary in y of particle partition
! edges(3) = lower boundary in z of particle partition
! edges(4) = upper boundary in z of particle partition
         real, dimension(:,:), pointer :: pt => null()
         integer :: npx, npy, npz, nx, ny, nz, ipbc
         real :: vtx, vty, vtz, vdx, vdy, vdz
         real :: sigx, sigy, sigz, x0, y0, z0
         real, dimension(3) :: cx, cy
         real, dimension(4) :: edges
         integer, dimension(2) :: noff
         integer :: nps=1
         logical :: lquiet = .false.
         integer :: idimp, npmax, ierr = 0
         character(len=18), save :: sname = 'dist3d_000:'

         call this%err%werrfl2(class//sname//' started')
         
         npx = this%npx; npy = this%npy; npz = this%npz
         nx = fd%getnd1(); ny = fd%getnd2(); nz = fd%getnd3()
         ipbc = this%sp%getpsolver()
         pt => part3d
         vtx = this%sigvx; vty = this%sigvy; vtz = this%sigvz
         vdx = 0.0; vdy = 0.0; vdz = this%gamma
         sigx = this%sigx; sigy = this%sigy; sigz = this%sigz
         x0 = this%bcx; y0 = this%bcy; z0 = this%bcz
         cx = (/this%cx1,this%cx2,this%cx3/); cy = (/this%cy1,this%cy2,this%cy3/)
         lquiet = this%quiet
         idimp = size(part3d,1); npmax = size(part3d,2)
         noff = fd%getnoff()
         edges(1) = noff(1); edges(3) = noff(2)
         edges(2) = edges(1) + fd%getnd2p()
         edges(4) = edges(3) + fd%getnd3p()         
         
         call PRVDIST32_RANDOM(pt,this%qm,edges,npp,nps,vtx,vty,vtz,vdx,vdy,&
         &vdz,npx,npy,npz,nx,ny,nz,ipbc,idimp,npmax,1,1,4,sigx,sigy,sigz,&
         &x0,y0,z0,cx,cy,lquiet,ierr)

         if (ierr /= 0) then
            write (erstr,*) 'PRVDIST32_RANDOM error'
            call this%err%equit(class//sname//erstr)
         endif
         
         call this%err%werrfl2(class//sname//' ended')
         
      end subroutine dist3d_000
!
      end module fdist3d_class