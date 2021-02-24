! fpois2d class for QuickPIC Open Source 1.0
! update: 04/18/2016

      module fpois2d_class

      use perrors_class
      use parallel_pipe_class
      use spect2d_class
      use ufield2d_class
      use fpois2d_lib
         
      implicit none

      private

      public :: fpois2d, get_pois2table

      type fpois2d

         private
! nd = system length in each direction
! a = half-width of particle in each direction
! anorm = normalization constant for poisson solver
! ffc = complex table for poisson solver
! ffg = real table for poisson solver
         class(spect2d), pointer, public :: sp => null()
         class(perrors), pointer, public :: err => null()
         class(parallel_pipe), pointer, public :: p => null()
         integer, dimension(2) :: nd
         real, dimension(2) :: a
         real :: anorm
         complex, dimension(:,:), pointer :: ffc => null()
         real, dimension(:,:,:), pointer :: ffg => null()
         
         contains
         
         generic :: new => init_fpois2d
         generic :: del => end_fpois2d
         generic :: potential => ipotd2
         generic :: vpotential => ivpotd2
         generic :: smoothf => ismoothfd2
         generic :: elfield => ippoisd23
         generic :: bfield => ibfieldd2
         generic :: bfield_qp => jpbpoisd23n_qp
         procedure, private :: init_fpois2d
         procedure, private :: end_fpois2d
         procedure, private :: ipotd2, ismoothfd2, ippoisd23, ibfieldd2
         procedure, private :: ivpotd2
         procedure, private :: jpbpoisd23n_qp                  
      end type 
!
      type fpois2d_link
         type (fpois2d_link), pointer :: next => null()
         type (fpois2d), pointer :: table => null()
         integer :: refcount
      end type fpois2d_link
!

      character(len=10), save :: class = 'fpois2d:'
      character(len=128), save :: erstr
! link list for poisson tables
      integer, save :: numtables = 0
      type (fpois2d_link), target, save :: table_list
            
      contains
!
      subroutine init_fpois2d(this,pp,perr,psp,nx,ny,ax,ay,affp)
      
         implicit none
         
         class(fpois2d), intent(inout) :: this
         class(spect2d), intent(in), pointer :: psp
         class(perrors), intent(in), pointer :: perr
         class(parallel_pipe), intent(in), pointer :: pp
         integer, intent(in) :: nx, ny
         real, intent(in) :: ax, ay, affp
         
! local data
         character(len=18), save :: sname = 'init_fpois2d:'
         
         this%sp => psp
         this%err => perr
         this%p => pp

         if ((ax < 0.) .or. (ay < 0.)) then
            write (erstr,*) 'invalid ax or ay=', ax, ay
            call this%err%equit(class//sname//erstr)
            return
         endif

         call this%err%werrfl2(class//sname//' started')
         select case (this%sp%getpsolver())
         case (1)
            this%nd = (/nx,ny/); this%a = (/ax,ay/)
            this%anorm = affp
            allocate(this%ffc(ny,(nx-1)/this%p%getlnvp()+1))
            call ippoisd2init(this,ax,ay,affp,nx,ny,this%p%getlkstrt())
         case default
            this%nd = (/nx,ny/); this%a = (/ax,ay/)
            this%anorm = affp
            allocate(this%ffc(ny,(nx-1)/this%p%getlnvp()+1))
            call ippoisd2init(this,ax,ay,affp,nx,ny,this%p%getlkstrt())
         end select

         call this%err%werrfl2(class//sname//' ended')

      end subroutine init_fpois2d
!
      subroutine end_fpois2d(this)
          
         implicit none
         
         class(fpois2d), intent(inout) :: this
         character(len=18), save :: sname = 'end_fpois2d'

         call this%err%werrfl2(class//sname//' started')

         if(associated(this%ffc)) deallocate(this%ffc)
         if(associated(this%ffg)) deallocate(this%ffg)
         
         call this%err%werrfl2(class//sname//' ended')
         
         return
         
      end subroutine end_fpois2d
!      
         subroutine ippoisd2init(this,ax,ay,affp,nx,ny,kstrt)
! this subroutine initializes ffd table
         implicit none
         class(fpois2d), intent(in) :: this
         integer, intent(in) :: nx, ny, kstrt
         real, intent(in) :: ax, ay, affp
! local data
         integer :: isign = 0, nyv, kxp2, nyd
         real :: we
         complex, dimension(:,:), pointer:: ffd
         real, dimension(1,1)  :: q
         real, dimension(2,1,1) :: fxy
         character(len=14), save :: sname = 'ippoisd2init:'

         call this%err%werrfl2(class//sname//' started')
         ffd => this%ffc
         nyv = size(q,1)
         nyd = size(ffd,1); kxp2 = size(ffd,2);
         call MPPOISD22(q,fxy,isign,ffd,ax,ay,affp,we,nx,ny,kstrt,nyv,kxp2,&
         &nyd)
         call this%err%werrfl2(class//sname//' ended')
         end subroutine ippoisd2init
!
         subroutine ipotd2(this,q,fx,we)
! this subroutine solves poisson's equation for potential with 
! dirichlet (zero potential) boundary conditions and 1d partition
! this = fpois2d descriptor
! q = input charge density array, in fourier space
! fx = output potential array, in fourier space
! we = output potential energy
         implicit none
         class(fpois2d), intent(in) :: this
         real, intent(inout) :: we
         class(ufield2d), intent(inout) :: q, fx
! local data
         integer :: isign = 1
         real, dimension(:,:,:), pointer :: pq => null(), pfx => null()
         character(len=8), save :: sname = 'ipotd2:'
         
         call this%err%werrfl2(class//sname//' started')
         select case (this%sp%getpsolver())
         case (1)
            pq => q%getrf(); pfx => fx%getrf()
            call ippoisd2(this,pq,pfx,isign,we)
         case default
            pq => q%getrf(); pfx => fx%getrf()
            call ippoisd2(this,pq,pfx,isign,we)
         end select
         call this%err%werrfl2(class//sname//' ended')
         end subroutine ipotd2
!
         subroutine ippoisd2(this,q,fx,isign,we)
! this subroutine solves poisson's equation for potential or smoothing
! with dirichlet (zero potential) boundary conditions and 1d partition
! this = poisson solver descriptor, includes table pointers
! q = input charge density array, in fourier space
! fx = output potential array, in fourier space
! isign = (1,2) = solve for (potential,smooth)
! we = output potential energy
         implicit none
         class(fpois2d), intent(in) :: this
         integer, intent(in) :: isign
         real, intent(inout) :: we
         real, dimension(:,:,:), pointer, intent(inout) :: q, fx
! local data
         integer :: nx, ny, nyv, kxp2, nyd
         complex, dimension(:,:), pointer :: ffc
         character(len=10), save :: sname = 'ippoisd2:'

         call this%err%werrfl2(class//sname//' started')
! unpack common arguments
         ffc => this%ffc
         nx = this%nd(1); ny = this%nd(2)
         nyv = size(q,2); nyd = size(ffc,1)
         kxp2 = size(q,3)-1;

         if (isign==1) then
            call MPPOTPD2(q(1,:,:),fx(1,:,:),ffc,we,nx,ny,this%p%getlkstrt(),&
            &nyv,kxp2,nyd)
         else if (isign==2) then
            call MPPSMOOTHD2(q(1,:,:),fx(1,:,:),ffc,nx,ny,this%p%getlkstrt(),&
            &nyv,kxp2,nyd)
         endif
         call this%err%werrfl2(class//sname//' ended')

         end subroutine ippoisd2
!
         subroutine ivpotd2(this,q,fx,we)
! this subroutine solves poisson's equation for vector potential with 
! dirichlet (zero potential) boundary conditions and 1d partition
! this = fpois2d descriptor
! q = input current density array, in fourier space
! fx = output potential array, in fourier space
! we = output potential energy
         implicit none
         class(fpois2d), intent(in) :: this
         real, intent(inout) :: we
         class(ufield2d), intent(inout) :: q, fx
! local data
         integer :: isign = 1
         real, dimension(:,:,:), pointer :: pq => null(), pfx => null()
         character(len=8), save :: sname = 'ivpotd2:'
         
         call this%err%werrfl2(class//sname//' started')
         select case (this%sp%getpsolver())
         case (1)
            pq => q%getrf(); pfx => fx%getrf()
            call ivppoisd2(this,pq,pfx,isign,we)
         case default
            pq => q%getrf(); pfx => fx%getrf()
            call ivppoisd2(this,pq,pfx,isign,we)
         end select
         call this%err%werrfl2(class//sname//' ended')
         end subroutine ivpotd2
!
         subroutine ivppoisd2(this,q,fx,isign,we)
! this subroutine solves poisson's equation for vetor potential
! with dirichlet (zero potential) boundary conditions and 1d partition
! this = poisson solver descriptor, includes table pointers
! q = input charge density array, in fourier space
! fx = output potential array, in fourier space
! isign = (1,2) = solve for (potential,smooth)
! we = output potential energy
         implicit none
         class(fpois2d), intent(in) :: this
         integer, intent(in) :: isign
         real, intent(inout) :: we
         real, dimension(:,:,:), pointer, intent(inout) :: q, fx
! local data
         integer :: nx, ny, nyv, kxp2, nyd
         complex, dimension(:,:), pointer :: ffc
         character(len=10), save :: sname = 'ivppoisd2:'

         call this%err%werrfl2(class//sname//' started')
! unpack common arguments
         ffc => this%ffc
         nx = this%nd(1); ny = this%nd(2)
         nyv = size(q,2); nyd = size(ffc,1)
         kxp2 = size(q,3)-1;

         if (isign==1) then
            call PBPOISD23(q,fx,isign,ffc,this%a(1),this%a(2),this%anorm,1.0,&
            &we,nx,ny,this%p%getlkstrt(),nyv,kxp2,1,nyd)
         endif
         call this%err%werrfl2(class//sname//' ended')

         end subroutine ivppoisd2
!
         subroutine ismoothfd2(this,f,fs)
! this subroutine provides smoothing for fourier transformed data with
! 1d partition
! this = fpois2d descriptor
! f = input data, in fourier space
! fs = output data, in fourier space
         implicit none
         class(fpois2d), intent(in) :: this
         class(ufield2d), intent(inout) :: f, fs
! local data
         integer :: isign = 2
         real :: we
         real, dimension(:,:,:), pointer :: pf, pfs
         character(len=12), save :: sname = 'ismoothfd2:'
         call this%err%werrfl2(class//sname//' started')
         select case (this%sp%getpsolver())
         case (1)
            pf => f%getrf(); pfs => fs%getrf()
            call ippoisd2(this,pf,pfs,isign,we)
         case default
            pf => f%getrf(); pfs => fs%getrf()
            call ippoisd2(this,pf,pfs,isign,we)
         end select
         call this%err%werrfl2(class//sname//' ended')
         end subroutine ismoothfd2
!
         subroutine ippoisd23(this,q,fxy,we)
! this subroutine solves poisson's equation for electric field
! with dirichlet (zero potential) boundary conditions and 1d partition
! this = fpois2d descriptor, includes table pointers
! q = input charge density array, in fourier space
! fxy = output electric field array, in fourier space
! we = output potential energy
! kstrt = starting data block number, a global variable
         implicit none
         class(fpois2d), intent(in) :: this
         class(ufield2d), intent(inout) :: q, fxy
         real, intent(inout) :: we
! local data
         integer :: isign = -1, nx, ny, nyv, kxp2, nyd
         real :: ax, ay, affp
         real, dimension(:,:,:), pointer :: pq
         real, dimension(:,:,:), pointer :: pfxy
         complex, dimension(:,:), pointer :: ffc
         character(len=11), save :: sname = 'ippoisd23:'

         call this%err%werrfl2(class//sname//' started')
         pq => q%getrf(); pfxy => fxy%getrf()         
         if ((size(pfxy,1) < 2) .or. (size(pfxy,1) > 3)) then
            write (erstr,*) 'invalid vector size=', size(pfxy,1)
            call this%err%equit(class//sname//erstr)
            return
         endif
! unpack common arguments
         ffc => this%ffc
         nx = this%nd(1); ny = this%nd(2)
! choose the proper solver
         select case (this%sp%getpsolver())
         case (1)
            nyv = size(pq,2); nyd = size(ffc,1)
            kxp2 = size(pq,3)-1;
            select case (size(pfxy,1))
            case (2)
               call MPPOISD22(pq(1,:,:),pfxy,isign,ffc,ax,ay,affp,we,nx,ny,&
               &this%p%getlkstrt(),nyv,kxp2,nyd)
            case (3)
               call MPPOISD23(pq(1,:,:),pfxy,isign,ffc,ax,ay,affp,we,nx,ny,&
               &this%p%getlkstrt(),nyv,kxp2,nyd)
            end select
         case default
            nyv = size(pq,2); nyd = size(ffc,1)
            kxp2 = size(ffc,2);
            select case (size(pfxy,1))
            case (2)
               call MPPOISD22(pq(1,:,:),pfxy,isign,ffc,ax,ay,affp,we,nx,ny,&
               &this%p%getlkstrt(),nyv,kxp2,nyd)
            case (3)
               call MPPOISD23(pq(1,:,:),pfxy,isign,ffc,ax,ay,affp,we,nx,ny,&
               &this%p%getlkstrt(),nyv,kxp2,nyd)
            end select
         end select
         call this%err%werrfl2(class//sname//' ended')

         end subroutine ippoisd23
!
         subroutine ibfieldd2(this,cu,bxy,ci,wm)
! this subroutine solves poisson's equation for magnetic field
! with dirichlet (zero potential) boundary conditions and 1d partition
! this = fpois2d descriptor, includes table pointers
! cu = input current density array, in fourier space
! bxy = output magnetic field array, in fourier space
! ci = reciprical of velocity of light
! wm = output magnetic energy
         implicit none
         class(fpois2d), intent(in) :: this
         real, intent(in) :: ci
         real, intent(inout) :: wm
         class(ufield2d), intent(inout) :: cu, bxy
! local data
         integer :: isign = -1
         real, dimension(:,:,:), pointer :: pcu, pbxy
         character(len=11), save :: sname = 'ibfieldd2:'

         call this%err%werrfl2(class//sname//' started')

         select case (this%sp%getpsolver())
         case (1)
            pcu => cu%getrf(); pbxy => bxy%getrf()         
            call jpbpoisd23(this,pcu,pbxy,isign,ci,wm)
         case default
            pcu => cu%getrf(); pbxy => bxy%getrf()         
            call jpbpoisd23(this,pcu,pbxy,isign,ci,wm)
         end select
         call this%err%werrfl2(class//sname//' ended')

         end subroutine ibfieldd2
!
         subroutine jpbpoisd23(this,cu,bxy,isign,ci,wm)
! this subroutine solves poisson's equation for magnetic field or
! vector potential or smoothing with dirichlet (zero potential) boundary
! conditions and 1d partition
! this = fpois2d descriptor, includes table pointers
! cu = input current density array, in fourier space
! bxy = output magnetic field array, in fourier space
! isign = (-1,1,2) = solve for (magnetic field,vector potential,smooth)
! ci = reciprical of velocity of light
! wm = output magnetic energy
! kstrt = starting data block number, a global variable
         implicit none
         class(fpois2d), intent(in) :: this
         integer, intent(in) :: isign
         real, intent(in) :: ci
         real, intent(inout) :: wm
         real, dimension(:,:,:), pointer, intent(inout) :: cu, bxy
! local data
         integer :: nx, ny, nyv, kxp2, j2blok, nyd
         real :: ax, ay, affp
         real, dimension(1,1,1) :: dxy
         complex, dimension(:,:), pointer :: ffc
         character(len=12), save :: sname = 'jpbpoisd23:'

         call this%err%werrfl2(class//sname//' started')
         if ((size(cu,1) < 2) .or. (size(cu,1) > 3)) then
            write (erstr,*) 'invalid cu vector size=', size(cu,1)
            call this%err%equit(class//sname//erstr)
            return            
         endif
         if (size(cu,1) /= size(bxy,1)) then
            if ((size(cu,1)==2).and.((size(bxy,1)/=1).or.(isign/=(-1))))&
     &then
               write (erstr,*) 'invalid bxy vector size=', size(bxy,1)
               call this%err%equit(class//sname//erstr)
               return            
            endif
         endif
! unpack common arguments
         ffc => this%ffc
         nx = this%nd(1); ny = this%nd(2)
         nyv = size(cu,2); nyd = size(ffc,1)
         kxp2 = size(cu,3)-1;
         select case (size(cu,1))
            case (3)
               call MPPBBPOISD23(cu,bxy,ffc,ci,wm,nx,ny,this%p%getlkstrt(),nyv,&
               &kxp2,nyd)
         end select
         call this%err%werrfl2(class//sname//' ended')

         end subroutine jpbpoisd23
!
         subroutine jpbpoisd23n_qp(this,cu,dcu,amu,bxy,ci,c,dex,wm)
! this = poisson solver descriptor
! cu = input current density array, in fourier space
! bxy = output magnetic field array, in fourier space
! isign = (-1,1,2) = solve for (magnetic field,vector potential,smooth)
! ci = reciprical of velocity of light
! wm = output magnetic energy
! kstrt = starting data block number, a global variable
         implicit none
         class(fpois2d), intent(in) :: this
         real, intent(in) :: ci, c, dex
         real, intent(inout) :: wm
         class(ufield2d), intent(inout) :: cu, dcu, amu, bxy
! local data
         integer :: isign = 1
         real, dimension(:,:,:), pointer :: pcu, pdcu, pamu, pbxy
         integer :: nx, ny, nyv, kxp2, j2blok, nyd
         real :: ax, ay, affp
         real, dimension(1,1,1) :: dxy
         complex, dimension(:,:), pointer :: ffc
         character(len=20), save :: sname = 'jpbpoisd23n_qp:'

         call this%err%werrfl2(class//sname//' started')
         pcu => cu%getrf()
         pdcu => dcu%getrf()
         pamu => amu%getrf()
         pbxy => bxy%getrf()
         if ((size(pcu,1) < 2) .or. (size(pcu,1) > 3)) then
            write (erstr,*) 'invalid cu vector size=', size(pcu,1)
            call this%err%equit(class//sname//erstr)
            return
         endif

! unpack common arguments
         ffc => this%ffc
         nx = this%nd(1); ny = this%nd(2)
! choose the proper solver
         select case (this%sp%getpsolver())
         case (1)
            nyv = size(pcu,2); nyd = size(ffc,1)
            kxp2 = size(ffc,2); j2blok = 1
            call PBPOISD22N_QP(pcu,pdcu,pamu,pbxy,dxy,isign,ffc,ax,ay,&
            &affp,ci,wm,nx,ny,this%p%getlkstrt(),nyv,kxp2,j2blok,nyd,c,dex)
         case default
            nyv = size(pcu,2); nyd = size(ffc,1)
            kxp2 = size(ffc,2); j2blok = 1
            call PBPOISD22N_QP(pcu,pdcu,pamu,pbxy,dxy,isign,ffc,ax,ay,&
            &affp,ci,wm,nx,ny,this%p%getlkstrt(),nyv,kxp2,j2blok,nyd,c,dex)
         end select
         call this%err%werrfl2(class//sname//' ended')

         end subroutine jpbpoisd23n_qp
!
      function get_pois2table(pp,perr,psp,ax,ay,affp) result(table)

         implicit none
         
         class(spect2d), intent(in), pointer :: psp
         class(perrors), intent(in), pointer :: perr
         class(parallel_pipe), intent(in), pointer :: pp         
         real, intent(in) :: ax, ay, affp
         type (fpois2d), pointer :: table
! local data
         type (fpois2d_link), pointer :: link => null()
         type (fpois2d), pointer :: ltable => null()
         integer :: nx, ny, ierr = 0
         character(len=18), save :: sname = 'get_pois2table:'

         call perr%werrfl2(class//sname//' started')
         nullify(table)

         if (numtables==0) then
            nullify(table_list%next,table_list%table)
            table_list%refcount = 0
         endif

         nx = 2**psp%getindx()
         ny = 2**psp%getindy()

         link => table_list
         table => link%table
! search link list of table to see if required table already exists
         do while (associated(table))
! found it
            if ((nx==table%nd(1)).and.(ny==table%nd(2)).and.(ax==table%a(1))&
            &.and.(ay==table%a(2)).and.(affp==table%anorm).and.(psp%getpsolver()&
            &==table%sp%getpsolver())) then
               link%refcount = link%refcount + 1
               call perr%werrfl2(class//sname//' ended')
               return
! check next table, create new empty table if end is reached
            else
               if (associated(link%next)) then
                  link => link%next
               else
                  allocate(link%next)
                  link => link%next
                  nullify(link%next,link%table)
                  link%refcount = 0
               endif
               table => link%table
            endif
         end do
! allocate table entries
         allocate(ltable)
         link%table => ltable
         table => link%table
         call table%new(pp,perr,psp,nx,ny,ax,ay,affp)
         link%refcount = 1
         numtables = numtables + 1
         call perr%werrfl2(class//sname//' ended')

      end function get_pois2table
!         
      end module fpois2d_class