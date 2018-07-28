! part2d class for QuickPIC Open Source 1.0
! update: 04/18/2016

      module part2d_class

      use perrors_class
      use parallel_pipe_class
      use spect2d_class
      use fdist2d_class
      use ufield2d_class
      use part2d_lib
      use hdf5io_class
      use mpi
               
      implicit none

      private

      public :: part2d

      type part2d

         private

         class(spect2d), pointer, public :: sp => null()
         class(perrors), pointer, public :: err => null()
         class(parallel_pipe), pointer, public :: p => null()
!
! qm = charge on particle, in units of e
! qbm = particle charge/mass ratio
! dt = time interval between successive calculations
! ci = reciprical of velocity of light
! xdim = dimension of the particle coordinates
! nbmax = size of buffer for passing particles between processors
! np = total number of particles
! npp = number of particles in current partition
! npmax = maximum number of particles in each partition
! part(:,:) = initial particle coordinates
! ppart(:,:,:) = particle coordinates for OpenMP
! nppmx, nppmx0, nbmaxp, ntmaxp, npbmx, irc, ncl, ihole, kpic = parameters for OpenMP
!         
         real :: qm, qbm, dt, ci
         integer :: npmax, nbmax, np, xdim, npp = 0
         real, dimension(:,:), pointer :: part => null()
         real, dimension(:,:,:), pointer :: ppart => null()
         integer :: nppmx, nppmx0, nbmaxp, ntmaxp, npbmx, irc = 0
         integer, dimension(:,:), pointer :: ncl => null()
         integer, dimension(:,:,:), pointer :: ihole => null()
         integer, dimension(:), pointer :: kpic => null()
         
         contains
         
         generic :: new => init_part2d
         generic :: renew => renew_part2d
         generic :: del => end_part2d
         generic :: qdp => qdeposit
         generic :: amjdp => amjdeposit
         generic :: push => partpush
         generic :: pmv => pmove
         generic :: extpsi => extractpsi
         generic :: pcp => partcopy
         generic :: pcb => partcopyback
         generic :: psend => pipesend_part2d
         generic :: precv => piperecv_part2d
         generic :: wr => writehdf5_part2d
         procedure, private :: init_part2d, renew_part2d
         procedure, private :: end_part2d
         procedure, private :: qdeposit
         procedure, private :: amjdeposit
         procedure, private :: partpush
         procedure, private :: pmove
         procedure, private :: extractpsi
         procedure, private :: partcopy
         procedure, private :: partcopyback
         procedure, private :: pipesend_part2d
         procedure, private :: piperecv_part2d, writehdf5_part2d
                           
      end type 

      save      

      character(len=10) :: class = 'part2d:'
      character(len=128) :: erstr
! parameters for OpenMP
      integer :: mx = 16, my = 16
      integer :: mx1, myp1, mxyp1
! ppbuff = buffer array for reordering tiled particle array
      real, dimension(:,:,:), allocatable :: ppbuff
      integer :: szpbuf = 0
! sbufl/sbufr = particle buffers sent to nearby processors
! rbufl/rbufr = particle buffers received from nearby processors
      real, dimension(:,:), allocatable :: sbufl, sbufr, rbufl, rbufr
      integer :: szbufs = 0
! ncll/nclr/mcll/mclr = number offsets send/received from processors
      integer, dimension(:,:), allocatable :: ncll, nclr, mcll, mclr
      integer :: sznbufs = 0
      
      contains
!
      subroutine init_part2d(this,pp,perr,psp,pf,fd,qm,qbm,dt,ci,xdim,npmax,nbmax,s)
      
         implicit none
         
         class(part2d), intent(inout) :: this
         class(spect2d), intent(in), pointer :: psp
         class(perrors), intent(in), pointer :: perr
         class(parallel_pipe), intent(in), pointer :: pp
         class(fdist2d), intent(inout) :: pf
         class(ufield2d), intent(in), pointer :: fd
         real, intent(in) :: qm, qbm, dt, ci, s
         integer, intent(in) :: npmax, nbmax, xdim

! local data
         character(len=18), save :: sname = 'init_part2d:'
         integer :: xtras, noff, nxyp, nx
                  
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
         noff = fd%getnoff()
         nxyp = fd%getnd2p()
         nx = fd%getnd1p()
         
         allocate(this%part(xdim,npmax))
         mx1 = (nx - 1)/mx + 1
         myp1 = (nxyp - 1)/my + 1; mxyp1 = mx1*myp1
         allocate(this%kpic(mxyp1))
         
         call pf%dist(this%part,this%npp,fd,s)

! find number of particles in each of mx, my tiles: updates kpic, nppmx
         call PPDBLKP2L(this%part,this%kpic,this%npp,noff,this%nppmx,&
         &this%xdim,this%npmax,mx,my,mx1,mxyp1,this%irc)
! check for errors
         if (this%irc /= 0) then
            write (erstr,*) 'PPDBLKP2L error, irc=', this%irc
            call this%err%equit(class//sname//erstr); return
         endif
         
         xtras = 10
!    
! allocate vector particle data
         this%nppmx0 = (1.0 + xtras)*this%nppmx
         this%ntmaxp = xtras*this%nppmx
         this%npbmx = xtras*this%nppmx
         this%nbmaxp = 0.25*mx1*this%npbmx
         allocate(this%ppart(xdim,this%nppmx0,mxyp1))
         allocate(this%ncl(8,mxyp1))
         allocate(this%ihole(2,this%ntmaxp+1,mxyp1))
!
! copy ordered particle data for OpenMP
         call PPPMOVIN2L(this%part,this%ppart,this%kpic,this%npp,noff,&
         &this%nppmx0,this%xdim,this%npmax,mx,my,mx1,mxyp1,this%irc)
! check for errors
         if (this%irc /= 0) then
            write (erstr,*) 'PPPMOVIN2L overflow error, irc=', this%irc
            call this%err%equit(class//sname//erstr); return
         endif
! 
! sanity check
!         call PPPCHECK2L(this%ppart,this%kpic,noff,nxyp,this%xdim,&
!         &this%nppmx0,nx,mx,my,mx1,myp1,this%irc)
! check error
!         if (this%irc /= 0) then
!            write (erstr,*) 'PPPCHECK2L error: irc=', this%irc
!            call this%err%equit(class//sname//erstr); return
!         endif

         call this%err%werrfl2(class//sname//' ended')

      end subroutine init_part2d
!
      subroutine end_part2d(this)
          
         implicit none
         
         class(part2d), intent(inout) :: this
         character(len=18), save :: sname = 'end_part2d:'

         call this%err%werrfl2(class//sname//' started')
         deallocate(this%part,this%ppart,this%ncl,this%ihole,this%kpic)
         call this%err%werrfl2(class//sname//' ended')
         
         return
         
      end subroutine end_part2d
!
      subroutine renew_part2d(this,pf,fd,s)
      
         implicit none
         
         class(part2d), intent(inout) :: this
         class(fdist2d), intent(inout) :: pf
         class(ufield2d), pointer, intent(in) :: fd
         real, intent(in) :: s

! local data
         character(len=18), save :: sname = 'renew_part2d:'
         integer :: noff, prof
                  
         call this%err%werrfl2(class//sname//' started')
         
         noff = fd%getnoff()
         prof = pf%getnpf()         
         
         call pf%dist(this%part,this%npp,fd,s)
         
         call PPDBLKP2L(this%part,this%kpic,this%npp,noff,this%nppmx,&
         &this%xdim,this%npmax,mx,my,mx1,mxyp1,this%irc)
! check for errors
         if (this%irc /= 0) then
            write (erstr,*) 'PPDBLKP2L error, irc=', this%irc
            call this%err%equit(class//sname//erstr); return
         endif
         
! copy ordered particle data for OpenMP
         call PPPMOVIN2L(this%part,this%ppart,this%kpic,this%npp,noff,&
         &this%nppmx0,this%xdim,this%npmax,mx,my,mx1,mxyp1,this%irc)
! check for errors
         if (this%irc /= 0) then
            write (erstr,*) 'PPPMOVIN2L overflow error, irc=', this%irc
            call this%err%equit(class//sname//erstr); return
         endif


         call this%err%werrfl2(class//sname//' ended')

      end subroutine renew_part2d
!      
      subroutine qdeposit(this,q)
! deposit the charge density      
      
         implicit none
         
         class(part2d), intent(in) :: this
         class(ufield2d), target, intent(inout) :: q
! local data
         character(len=18), save :: sname = 'qdeposit:'
         real, dimension(:,:,:), pointer :: pq => null()
         integer :: noff, nxv, nypmx
                  
         call this%err%werrfl2(class//sname//' started')
         
         pq => q%getrf()
         noff = q%getnoff()
         nxv = size(pq,2)
         nypmx = size(pq,3)
         
         call PPGPPOST2L(this%ppart,pq(1,:,:),this%kpic,noff,&
         &this%xdim,this%nppmx0,mx,my,nxv,nypmx,mx1,mxyp1)         
         
         call this%err%werrfl2(class//sname//' ended')
         
      end subroutine qdeposit
!      
      subroutine amjdeposit(this,ef,bf,psit,cu,amu,dcu,dex)
! deposit the current, acceleration and momentum flux      
      
         implicit none
         
         class(part2d), intent(inout) :: this
         class(ufield2d), target, intent(inout) :: cu, amu, dcu
         class(ufield2d), target, intent(in) :: ef, bf, psit
         real, intent(in) :: dex
         character(len=18), save :: sname = 'amjdeposit'
! local data
         real, dimension(:,:,:), pointer :: pef => null(), pbf => null()
         real, dimension(:,:,:), pointer :: ppsit => null(), pcu => null()
         real, dimension(:,:,:), pointer :: pamu => null(), pdcu => null()
         integer :: noff, nyp, nx, nxv, nypmx

         call this%err%werrfl2(class//sname//' started')
         
         pef => ef%getrf(); pbf => bf%getrf()
         ppsit => psit%getrf(); pcu => cu%getrf()
         pamu => amu%getrf(); pdcu => dcu%getrf()
         noff = ef%getnoff()
         nxv = size(pef,2); nypmx = size(pef,3)
         nx = ef%getnd1(); nyp = ef%getnd2p()
         
         call PPGRDCJPPOST2L_QP(this%ppart,pef,pbf,ppsit(1,:,:),pcu,pdcu,&
         &pamu,this%kpic,noff,nyp,this%qbm, this%dt,this%ci,this%xdim,&
         &this%nppmx0,nx,mx,my,nxv,nypmx,mx1,mxyp1,dex)
         
         call this%err%werrfl2(class//sname//' ended')
         
      end subroutine amjdeposit
!      
      subroutine partpush(this,ef,bf,psit,dex)
      
         implicit none
         
         class(part2d), intent(inout) :: this
         class(ufield2d), target, intent(in) :: ef, bf, psit
         real, intent(in) :: dex
         character(len=18), save :: sname = 'partpush'
! local data
         real, dimension(:,:,:), pointer :: pef => null(), pbf => null()
         real, dimension(:,:,:), pointer :: ppsit => null()
         integer :: noff, nyp, nx, ny, nxv, nypmx
         real :: ek
         
         call this%err%werrfl2(class//sname//' started')

         pef => ef%getrf(); pbf => bf%getrf()
         ppsit => psit%getrf()
         noff = ef%getnoff(); ny = ef%getnd2()
         nxv = size(pef,2); nypmx = size(pef,3)
         nx = ef%getnd1(); nyp = ef%getnd2p()

         call PPGRBPPUSHF23L_QP(this%ppart,pef,pbf,ppsit(1,:,:),this%kpic,&
         &this%ncl,this%ihole,noff,nyp,this%qbm,this%dt,this%dt,this%ci,ek,&
         &this%xdim,this%nppmx0,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,this%ntmaxp,&
         &this%irc,dex)

         call this%err%werrfl2(class//sname//' ended')
         
      end subroutine partpush
!
      subroutine pmove(this,fd)
      
         implicit none
         
         class(part2d), intent(inout) :: this
         class(ufield2d), target, intent(in) :: fd
         character(len=18), save :: sname = 'pmove:'
! local data
! list = (true,false) = list of particles leaving tiles found in push
         integer :: noff, nyp, nx, ny, nxv, nypmx, kstrt, nvp
         integer :: npbmx, nbmax, idimp, nppmx, ntmax, irc
         real, dimension(:,:,:), pointer :: ppart => null()
         integer, dimension(:,:), pointer :: ncl => null()
         integer, dimension(:,:,:), pointer :: ihole => null()
         integer, dimension(:), pointer :: kpic => null()
         logical :: list = .true.
         
         call this%err%werrfl2(class//sname//' started')
         
         noff = fd%getnoff(); ny = fd%getnd2()
         nx = fd%getnd1(); nyp = fd%getnd2p()
         npbmx = this%npbmx; nbmax = this%nbmaxp
         idimp = this%xdim; nppmx = this%nppmx0
         ntmax = this%ntmaxp; ppart => this%ppart
         kstrt = this%p%getlkstrt(); nvp = this%p%getlnvp()
         ncl => this%ncl; ihole => this%ihole; kpic => this%kpic
         irc = this%irc
! check if required size of buffer has increased
         if (szpbuf < idimp*npbmx*mxyp1) then
            if (szpbuf /= 0) deallocate(ppbuff)
! allocate new buffer
            allocate(ppbuff(idimp,npbmx,mxyp1))
            szpbuf = idimp*npbmx*mxyp1
         endif
! check if required size of buffers has increased
         if (szbufs < idimp*nbmax) then
            if (szbufs /= 0) deallocate(sbufl,sbufr,rbufl,rbufr)
! allocate new buffers
            allocate(sbufl(idimp,nbmax),sbufr(idimp,nbmax))
            allocate(rbufl(idimp,nbmax),rbufr(idimp,nbmax))
            szbufs = idimp*nbmax
         endif
! check if required size of buffers has increased
         if (sznbufs < 3*mx1) then
            if (sznbufs /= 0) deallocate(ncll,nclr,mcll,mclr)
! allocate new buffers
            allocate(ncll(3,mx1),nclr(3,mx1),mcll(3,mx1),mclr(3,mx1))
            sznbufs = 3*mx1
         endif
!
! first part of particle reorder on x and y cell with mx, my tiles:
! list of particles leaving tile already calculated by push
         if (list) then
! updates: ppart, ppbuff, sbufl, sbufr, ncl, ncll, nclr, irc
            call PPPORDERF2LA(ppart,ppbuff,sbufl,sbufr,ncl,ihole,ncll,nclr,&
            &idimp,nppmx,mx1,myp1,npbmx,ntmax,nbmax,irc)
            if (irc /= 0) then
               write (erstr,*) kstrt,'mporderf2a error: ntmax, irc=',ntmax,irc            
            endif
! calculate list of particles leaving tile
         else
! updates ppart, ppbuff, sbufl, sbufr, ncl, ihole, ncll, nclr, irc
            call PPPORDER2LA(ppart,ppbuff,sbufl,sbufr,kpic,ncl,ihole,ncll,nclr&
            &,noff,nyp,idimp,nppmx,nx,ny,mx,my,mx1,myp1,npbmx,ntmax,nbmax,irc)
            if (irc /= 0) then
               write (erstr,*) kstrt,'mporder2a error: ntmax, irc=',ntmax,irc
            endif
         endif
         if (irc /= 0) then
               call this%err%equit(class//sname//erstr)
               return
         endif
!
! move particles into appropriate spatial regions with MPI:
! updates rbufr, rbufl, mcll, mclr
         call PPPMOVE2(sbufr,sbufl,rbufr,rbufl,ncll,nclr,mcll,mclr,kstrt,&
         &nvp,idimp,nbmax,mx1)
!
! second part of particle reorder on x and y cell with mx, my tiles:
! updates ppart, kpic
         call PPPORDER2LB(ppart,ppbuff,rbufl,rbufr,kpic,ncl,ihole,mcll,mclr,&
         &idimp,nppmx,mx1,myp1,npbmx,ntmax,nbmax,irc)
         
         if (irc /= 0) then
            write (erstr,*) kstrt,'mporder2b error: nppmx, irc=',nppmx,irc
            call this%err%equit(class//sname//erstr)
            stop
         endif

         this%irc = irc

         call this%err%werrfl2(class//sname//' ended')
         
      end subroutine pmove
!      
      subroutine extractpsi(this,psi,dex)
      
         implicit none
         
         class(part2d), intent(inout) :: this
         class(ufield2d), target, intent(in) :: psi
         real, intent(in) :: dex
         character(len=18), save :: sname = 'extractpsi'
! local data
         real, dimension(:,:,:), pointer :: ppsi
         integer :: noff, nyp, nx, nxv, nypmx
         
         call this%err%werrfl2(class//sname//' started')

         ppsi => psi%getrf(); noff = psi%getnoff()
         nyp = psi%getnd2p(); nx = psi%getnd1()
         nxv = size(ppsi,2); nypmx = size(ppsi,3)
         
         call WPGPSIPOST2L_QP(this%ppart,ppsi(1,:,:),this%kpic,this%qbm,noff,&
         &nyp,this%xdim,this%nppmx0,nx,mx,my,nxv,nypmx,mx1,mxyp1,dex)

         call this%err%werrfl2(class//sname//' ended')
         
      end subroutine extractpsi
!
      subroutine partcopy(this,fd)
      
         implicit none
         
         class(part2d), intent(inout) :: this
         class(ufield2d), target, intent(in) :: fd
! local data
         character(len=18), save :: sname = 'partcopy:'
         integer :: noff
         
         call this%err%werrfl2(class//sname//' started')

         noff = fd%getnoff()         
         
! copy ordered particle data for OpenMP
         call PPPMOVIN2L(this%part,this%ppart,this%kpic,this%npp,noff,&
         &this%nppmx0,this%xdim,this%npmax,mx,my,mx1,mxyp1,this%irc)
! check for errors
         if (this%irc /= 0) then
            write (erstr,*) 'PPPMOVIN2L overflow error, irc=', this%irc
            call this%err%equit(class//sname//erstr); return
         endif
         
         
         call this%err%werrfl2(class//sname//' ended')
         
      end subroutine partcopy
!
      subroutine partcopyback(this)
      
         implicit none
         
         class(part2d), intent(inout) :: this
! local data
         character(len=18), save :: sname = 'partcopyback:'
         
         call this%err%werrfl2(class//sname//' started')

         this%irc = 0

         call PPPCOPYOUT2(this%part,this%ppart,this%kpic,this%npp,&
         &this%npmax,this%nppmx0,this%xdim,mxyp1,this%irc)

         call this%err%werrfl2(class//sname//' ended')
         
      end subroutine partcopyback  
!
      subroutine pipesend_part2d(this,tag,id)
      
         implicit none
         
         class(part2d), intent(inout) :: this
         integer, intent(in) :: tag
         integer, intent(inout) :: id
! local data
         character(len=18), save :: sname = 'pipesend_part2d:'
         integer :: des, ierr
         
         
         call this%err%werrfl2(class//sname//' started')
         
         des = this%p%getidproc()+this%p%getlnvp()
         
         if (des >= this%p%getnvp()) then
            id = MPI_REQUEST_NULL         
            call this%err%werrfl2(class//sname//' ended')
            return
         endif
         
         call this%pcb()
                  
         call MPI_ISEND(this%part,this%npp*this%xdim,this%p%getmreal(),&
         &des,tag,this%p%getlworld(),id,ierr)

! check for errors
         if (ierr /= 0) then
            write (erstr,*) 'MPI_ISEND failed'
            call this%err%equit(class//sname//erstr); return
         endif

         call this%err%werrfl2(class//sname//' ended')
         
      end subroutine pipesend_part2d
!      
      subroutine piperecv_part2d(this,fd,tag)
      
         implicit none
         
         class(part2d), intent(inout) :: this
         class(ufield2d), target, intent(in) :: fd
         integer, intent(in) :: tag
! local data
         character(len=18), save :: sname = 'piperecv_part2d:'
         integer, dimension(10) :: istat
         integer :: nps, id, des, ierr
         
         
         call this%err%werrfl2(class//sname//' started')

         des = this%p%getidproc()-this%p%getlnvp()
         
         if (des < 0) then
            call this%err%werrfl2(class//sname//' ended')
            return
         endif

         call MPI_IRECV(this%part,this%npmax*this%xdim,this%p%getmreal(),&
         &des,tag,this%p%getlworld(),id,ierr)

         call MPI_WAIT(id,istat,ierr)
         
         call MPI_GET_COUNT(istat,this%p%getmreal(),nps,ierr)

         this%npp = nps/this%xdim
         
         call this%pcp(fd)
         
! check for errors
         if (ierr /= 0) then
            write (erstr,*) 'MPI failed'
            call this%err%equit(class//sname//erstr); return
         endif
         
         call this%err%werrfl2(class//sname//' ended')
         
      end subroutine piperecv_part2d
!      
      subroutine writehdf5_part2d(this,file,delta)

         implicit none
         
         class(part2d), intent(inout) :: this
         class(hdf5file), intent(in) :: file
         real, dimension(2), intent(in) :: delta
! local data
         character(len=18), save :: sname = 'writehdf5_part2d:'
         integer :: ierr

         call this%err%werrfl2(class//sname//' started')                  
         call this%pcb()
         call pwpart(this%p,this%err,file,this%part,this%npp,1,delta,ierr)
         call this%err%werrfl2(class//sname//' ended')
      
      end subroutine writehdf5_part2d
!      
      end module part2d_class