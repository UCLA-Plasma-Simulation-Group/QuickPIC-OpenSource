! simulation module for QuickPIC Open Source 1.0
! update: 01/09/2018

      module simulation_class
      
      use parallel_class
      use parallel_pipe_class
      use perrors_class
      use spect3d_class
      use spect2d_class
      use field2d_class
      use fdist2d_class
      use field3d_class
      use beam3d_class
      use species2d_class
      use fdist3d_class
      use hdf5io_class
      use input_class
      use mpi

      implicit none
      
      private

      public :: simulation
!
      type sim_fields

         private

         class(parallel_pipe),pointer :: p => null()
         class(perrors),pointer :: err => null()
         class(spect3d), pointer :: sp3 => null()
         class(spect2d), pointer :: sp2 => null()
         type(field2d), allocatable :: qb, qe, qi, psit, psi, div_vpot, reg
         type(field2d), allocatable :: fxy, bxyz, cu, dcu, amu, epw, epwb         
         type(field2d), dimension(:), allocatable :: qe0, cu0, dcu0, amu0      
         type(field3d), allocatable :: bexyz, bbxyz, qeb
         type(field3d), dimension(:), allocatable :: qep
         type(field3d), allocatable :: psi3d,cu3d

         contains
         
         generic :: new => init_sim_fields
         generic :: del => end_sim_fields

         procedure, private :: init_sim_fields, end_sim_fields
         
      end type sim_fields
!
      type sim_beams

         private

         class(parallel_pipe),pointer :: p => null()
         class(perrors),pointer :: err => null()
         class(spect3d), pointer :: sp3 => null()
         class(spect2d), pointer :: sp2 => null()
         type(beam3d), dimension(:), allocatable :: beam
         type(fdist3d), dimension(:), allocatable :: pf3
         
         contains
         
         generic :: new => init_sim_beams
         generic :: del => end_sim_beams

         procedure, private :: init_sim_beams, end_sim_beams

      end type sim_beams
!      
      type sim_species

         private

         class(parallel_pipe),pointer :: p => null()
         class(perrors),pointer :: err => null()
         class(spect3d), pointer :: sp3 => null()
         class(spect2d), pointer :: sp2 => null()
         type(fdist2d), dimension(:), allocatable :: pf
         type(species2d), dimension(:), allocatable :: spe

         contains
         
         generic :: new => init_sim_species
         generic :: del => end_sim_species

         procedure, private :: init_sim_species, end_sim_species

      end type sim_species
!      
      type simulation

         private

!         type(input) :: sim
         type(input_json), pointer :: in => null()
         class(parallel_pipe),pointer :: p => null()
         class(perrors),pointer :: err => null()
         class(spect3d), pointer :: sp3 => null()
         class(spect2d), pointer :: sp2 => null()
         
         type(sim_fields) :: fields
         type(sim_beams) :: beams
         type(sim_species) :: species
         integer :: iter, nstep3d, nstep2d, start3d,nbeams,nspecies
         integer, dimension(54) :: id
         integer, dimension(8) :: tag
         integer, dimension(:), allocatable :: tag_spe, id_spe, id_br
         integer, dimension(:), allocatable :: tag_beam, id_beam, id_qep, id_qeps
         real :: dex, dxi, dex2

         contains
         
         generic :: new => init_simulation
         generic :: del => end_simulation
         generic :: go => go_simulation
         generic :: diag => diag_simulation

         procedure, private :: init_simulation, end_simulation
         procedure, private :: go_simulation, diag_simulation

      end type simulation
!
      character(len=20) :: class = 'simulation: '
      character(len=128) :: erstr
            
      contains
!
      subroutine init_simulation(this)

         implicit none
         
         class(simulation), intent(inout) :: this
! local data
         character(len=18), save :: sname = 'init_simulation:'
         real :: min, max, cwp, n0, dx, dy, dz
         integer :: indx, indy, indz
         logical :: read_rst
         
         allocate(this%in)
         call this%in%new()
         this%err => this%in%err
         this%p => this%in%pp
         this%sp3 => this%in%sp
         this%sp2 => this%in%sp

         call this%err%werrfl2(class//sname//' started')

         call this%in%get('simulation.n0',n0)
         call this%in%get('simulation.indx',indx)
         call this%in%get('simulation.indy',indy)
         call this%in%get('simulation.indz',indz)

         cwp=5.32150254*1e9/sqrt(n0)
         call this%in%get('simulation.box.x(1)',min)
         call this%in%get('simulation.box.x(2)',max)
         dx=(max-min)/cwp/real(2**indx)
         call this%in%get('simulation.box.y(1)',min)
         call this%in%get('simulation.box.y(2)',max)
         dy=(max-min)/cwp/real(2**indy)
         call this%in%get('simulation.box.z(1)',min)
         call this%in%get('simulation.box.z(2)',max)
         dz=(max-min)/cwp/real(2**indz)
         
         this%dex = dx
         this%dxi = dz
         this%dex2 = dx * dx
         this%nstep2d = 2**this%sp3%getindz()/this%p%getnstage()
         call this%in%get('simulation.time',min)
         call this%in%get('simulation.dt',max)
         this%nstep3d = min/max
         call this%in%get('simulation.read_restart',read_rst)
         if (read_rst) then
            call this%in%get('simulation.restart_timestep',this%start3d)
         else
            this%start3d = 1
         end if
         call this%in%get('simulation.iter',this%iter)
         call this%in%get('simulation.nbeams',this%nbeams)
         call this%in%get('simulation.nspecies',this%nspecies)

         call this%fields%new(this%in)
         call this%beams%new(this%in,this%fields)
         call this%species%new(this%in,this%fields)

         allocate(this%tag_spe(this%nspecies),this%id_spe(this%nspecies))
         allocate(this%id_qep(this%nspecies),this%id_qeps(this%nspecies))
         allocate(this%tag_beam(this%nbeams),this%id_beam(this%nbeams),this%id_br(this%nbeams))
         this%id(:) = MPI_REQUEST_NULL
         this%id_spe(:) = MPI_REQUEST_NULL
         this%id_qep(:) = MPI_REQUEST_NULL
         this%id_qeps(:) = MPI_REQUEST_NULL
         this%id_beam(:) = MPI_REQUEST_NULL                 
         this%id_br(:) = MPI_REQUEST_NULL                 

         call this%err%werrfl2(class//sname//' ended')

      end subroutine init_simulation
!         
      subroutine end_simulation(this)

         implicit none

         class(simulation), intent(inout) :: this
! local data
         character(len=18), save :: sname = 'end_simulation:'
         integer :: ierr
         
         call this%err%werrfl2(class//sname//' started')
         
!         call this%fields%del()
!         call this%beams%del()
!         call this%species%del()
         call MPI_FINALIZE(ierr)

         call this%err%werrfl2(class//sname//' ended')

      end subroutine end_simulation
!
      subroutine init_sim_fields(this,input)

         implicit none

         class(sim_fields), intent(inout) :: this
         type(input_json), pointer, intent(inout) :: input
! local data
         character(len=18), save :: sname = 'init_sim_fields:'
         character(len=18), save :: class = 'sim_fields:'
         character(len=20) :: sn
         character(len=:), allocatable :: ff   
         integer :: i,n,ndump

         this%err => input%err
         this%p => input%pp
         this%sp3 => input%sp
         this%sp2 => input%sp


         call this%err%werrfl2(class//sname//' started')

         allocate(this%qb, this%qe, this%qi, this%psit, this%psi)
         allocate(this%div_vpot, this%reg, this%fxy, this%bxyz, this%cu)
         allocate(this%dcu, this%amu, this%epw, this%epwb)
         allocate(this%bexyz, this%bbxyz, this%qeb)

         call this%bexyz%new(this%p,this%err,this%sp3,dim=3)
         call this%bbxyz%new(this%p,this%err,this%sp3,dim=3)
         call this%qeb%new(this%p,this%err,this%sp3,dim=1)
         
         call this%qb%new(this%p,this%err,this%sp2,dim=1,fftflag=.true.,gcells=1)
         call this%qe%new(this%p,this%err,this%sp2,dim=1,fftflag=.true.)
         call this%qi%new(this%p,this%err,this%sp2,dim=1,fftflag=.true.,gcells=1)
         call this%psit%new(this%p,this%err,this%sp2,dim=1,fftflag=.true.)
         call this%div_vpot%new(this%p,this%err,this%sp2,dim=1,fftflag=.true.)
         call this%psi%new(this%p,this%err,this%sp2,dim=1,fftflag=.true.)
         call this%reg%new(this%p,this%err,this%sp2,dim=1,fftflag=.true.)
         call this%fxy%new(this%p,this%err,this%sp2,dim=2,fftflag=.true.)
         call this%cu%new(this%p,this%err,this%sp2,dim=3,fftflag=.true.,state=1)
         call this%dcu%new(this%p,this%err,this%sp2,dim=2,fftflag=.true.)
         call this%amu%new(this%p,this%err,this%sp2,dim=3,fftflag=.true.)
         call this%epw%new(this%p,this%err,this%sp2,dim=2,fftflag=.true.,state=1)
         call this%epwb%new(this%p,this%err,this%sp2,dim=2,fftflag=.true.)
         call this%bxyz%new(this%p,this%err,this%sp2,dim=3,fftflag=.true.)

         call input%get('simulation.nspecies',n)

         allocate(this%qe0(n),this%cu0(n),this%dcu0(n),this%amu0(n))
         allocate(this%qep(n))

         do i = 1, n
            call this%qep(i)%new(this%p,this%err,this%sp3,dim=1)
            call this%qe0(i)%new(this%p,this%err,this%sp2,dim=1,fftflag=.false.)
            call this%cu0(i)%new(this%p,this%err,this%sp2,dim=3,fftflag=.false.)
            call this%dcu0(i)%new(this%p,this%err,this%sp2,dim=2,fftflag=.false.)
            call this%amu0(i)%new(this%p,this%err,this%sp2,dim=3,fftflag=.false.)
         end do
         
         call input%info('field.diag',n_children=n)
         
         do i = 1, n
            write (sn,'(I4.4)') i
            call input%get('field.diag('//trim(sn)//').ndump',ndump)
            if (ndump > 0) then
               call input%get('field.diag('//trim(sn)//').name',ff)
               if (ff == 'psi') then
                  deallocate(ff)
                  allocate(this%psi3d)
                  call this%psi3d%new(this%p,this%err,this%sp3,dim=1)
                  exit
               end if
            end if
         end do
         
         do i = 1, n
            write (sn,'(I4.4)') i
            call input%get('field.diag('//trim(sn)//').ndump',ndump)
            if (ndump > 0) then
               call input%get('field.diag('//trim(sn)//').name',ff)
               if (ff == 'jp') then
                  deallocate(ff)
                  allocate(this%cu3d)
                  call this%cu3d%new(this%p,this%err,this%sp3,dim=1)
                  exit
               end if
            end if
         end do

         call this%err%werrfl2(class//sname//' ended')

      end subroutine init_sim_fields
!
      subroutine end_sim_fields(this)

         implicit none

         class(sim_fields), intent(inout) :: this
! local data
         character(len=18), save :: sname = 'end_sim_fields:'
         character(len=18), save :: class = 'sim_fields:'
         integer :: i, n
         
         call this%err%werrfl2(class//sname//' started')

         call this%bexyz%del()
         call this%bbxyz%del()
         call this%qeb%del()
         call this%qb%del()
         call this%qe%del()
         call this%qi%del()
         call this%psit%del()
         call this%div_vpot%del()
         call this%psi%del()
         call this%reg%del()
         call this%fxy%del()
         call this%cu%del()
         call this%dcu%del()
         call this%amu%del()
         call this%epw%del()
         call this%epwb%del()
         call this%bxyz%del()

         n = size(this%qep)

         do i = 1, n
            call this%qep(i)%del()
            call this%qe0(i)%del()
            call this%cu0(i)%del()
            call this%dcu0(i)%del()
            call this%amu0(i)%del()
         end do

         if(allocated(this%psi3d)) call this%psi3d%del()
         if(allocated(this%cu3d)) call this%cu3d%del()

         call this%err%werrfl2(class//sname//' ended')

      end subroutine end_sim_fields
!
      subroutine init_sim_beams(this,input,fields)

         implicit none

         class(sim_beams), intent(inout) :: this
         type(input_json), pointer, intent(inout) :: input
         class(sim_fields), intent(inout) :: fields
! local data
         character(len=18), save :: class = 'sim_beams:'
         character(len=18), save :: sname = 'init_sim_beams:'
         integer :: i,n
         integer :: npf,npx,npy,npz,npp
         real, dimension(3,100) :: arg
         logical :: quiet
         real :: min, max, cwp, n0, gamma
         real :: qm, qbm, dt, ntemp, alx, aly, alz
         real :: dx, dy, dz
         logical :: read_rst
         integer :: rst_timestep, ierr, indx, indy, indz, npmax
         type(hdf5file) :: file_rst
         character(len=20) :: sn, sid, stime,s1

         this%err => input%err
         this%p => input%pp
         this%sp3 => input%sp
         this%sp2 => input%sp
         
         call this%err%werrfl2(class//sname//' started')

         call input%get('simulation.n0',n0)
         call input%get('simulation.indx',indx)
         call input%get('simulation.indy',indy)
         call input%get('simulation.indz',indz)

         cwp=5.32150254*1e9/sqrt(n0)
         call input%get('simulation.box.x(1)',min)
         call input%get('simulation.box.x(2)',max)
         alx = (max-min)/cwp 
         dx=alx/real(2**indx)
         call input%get('simulation.box.y(1)',min)
         call input%get('simulation.box.y(2)',max)
         aly = (max-min)/cwp 
         dy=aly/real(2**indy)
         call input%get('simulation.box.z(1)',min)
         call input%get('simulation.box.z(2)',max)
         alz = (max-min)/cwp 
         dz=alz/real(2**indz)

         call input%get('simulation.nbeams',n)         

         allocate(this%beam(n),this%pf3(n))
         
         do i = 1, n
            arg(:,:) = 0.0
            write (sn,'(I3.3)') i
            s1 = 'beam('//trim(sn)//')'
            call input%get(trim(s1)//'.profile',npf)
            call input%get(trim(s1)//'.np(1)',npx)
            call input%get(trim(s1)//'.np(2)',npy)
            call input%get(trim(s1)//'.np(3)',npz)

            call input%get(trim(s1)//'.parameter(3)(1)',arg(1,1))
            call input%get(trim(s1)//'.parameter(3)(2)',arg(2,1))
            call input%get(trim(s1)//'.parameter(3)(3)',arg(3,1))
            arg(1,2) = 0.0; arg(2,2) = 0.0
            call input%get(trim(s1)//'.gamma',arg(3,2))
            call input%get(trim(s1)//'.parameter(2)(1)',arg(1,3))
            call input%get(trim(s1)//'.parameter(2)(2)',arg(2,3))
            call input%get(trim(s1)//'.parameter(2)(3)',arg(3,3))
            call input%get(trim(s1)//'.parameter(1)(1)',arg(1,4))
            call input%get(trim(s1)//'.parameter(1)(2)',arg(2,4))
            call input%get(trim(s1)//'.parameter(1)(3)',arg(3,4))
            call input%get(trim(s1)//'.parameter(4)(1)',arg(1,5))
            call input%get(trim(s1)//'.parameter(4)(2)',arg(1,6))
            call input%get(trim(s1)//'.parameter(4)(3)',arg(1,7))
            call input%get(trim(s1)//'.parameter(5)(1)',arg(2,5))
            call input%get(trim(s1)//'.parameter(5)(2)',arg(2,6))
            call input%get(trim(s1)//'.parameter(5)(3)',arg(2,7))
            call input%get(trim(s1)//'.quiet_start',quiet)
            call input%get(trim(s1)//'.gamma',gamma)
            if (quiet) then
               arg(1,8) = 1.0
            else
               arg(1,8) = 0.0
            endif

            arg(1,1) = arg(1,1)/arg(1,3)
            arg(2,1) = arg(2,1)/arg(2,3)
            arg(3,1) = arg(3,1)*gamma
            arg(1:3,3) = arg(1:3,3)/cwp
            arg(1:3,4) = arg(1:3,4)/cwp
            arg(1,5) = arg(1,5)*cwp 
            arg(1,7) = arg(1,7)/cwp
            arg(2,5) = arg(2,5)*cwp 
            arg(2,7) = arg(2,7)/cwp
            arg(1,3) = arg(1,3)/dx
            arg(1,4) = arg(1,4)/dx
            arg(2,3) = arg(2,3)/dy
            arg(2,4) = arg(2,4)/dy
            arg(3,3) = arg(3,3)/dz
            arg(3,4) = arg(3,4)/dz
            arg(1,5) = arg(1,5)*dz*dz/dx
            arg(1,6) = arg(1,6)*dz/dx
            arg(1,7) = arg(1,7)/dz*dz/dx
            arg(2,5) = arg(2,5)*dz*dz/dy
            arg(2,6) = arg(2,6)*dz/dy
            arg(2,7) = arg(2,7)/dz*dz/dy 

            call this%pf3(i)%new(this%p,this%err,this%sp3,npf=npf,npx=npx,&
            &npy=npy,npz=npz,arg=arg)

            call input%get(trim(s1)//'.q',qm)
            call input%get(trim(s1)//'.m',qbm)
            qbm = qm/qbm
            call input%get(trim(s1)//'.num_particle',npp)
            ntemp = npp*1e12*(2**indz)
            ntemp = ntemp*(2**indx)
            ntemp = ntemp*(2**indy)/(npx*alx*aly*alz*n0*cwp*cwp*cwp) 
            ntemp = ntemp/npy
            ntemp = ntemp/npz
            qm = ntemp*qm
               
            call input%get('simulation.dt',dt)
            call input%get(trim(s1)//'.npmax',npmax)
            
            call this%beam(i)%new(this%p,this%err,this%sp3,this%pf3(i),fields%qeb,qm=qm,qbm=qbm,&
            &dt=dt,ci=1.0,xdim=6,npmax=npmax,nbmax=int(npmax*0.01))

            call input%get('simulation.read_restart',read_rst)

            if (read_rst) then
               call input%get('simulation.restart_timestep',rst_timestep)
               write (sn,'(I2.2)') i
               write (sid,'(I8.8)') this%p%getidproc()
               write (stime,'(I8.8)') rst_timestep
               call file_rst%new(filename = './RST/RST-'//trim(sn)//'-'//trim(sid)//&
               &'_'//trim(stime)//'.h5',dataname = 'RST-'//trim(sn)//'-'//trim(sid)//&
               &'_'//trim(stime))
               call this%beam(i)%rrst(file_rst)
            end if
         end do
         call MPI_BARRIER(this%p%getlworld(),ierr)

         call this%err%werrfl2(class//sname//' ended')
      end subroutine init_sim_beams
!
      subroutine end_sim_beams(this)

         implicit none

         class(sim_beams), intent(inout) :: this
         type(input_json), intent(in) :: input
! local data
         character(len=18), save :: sname = 'end_sim_beams:'
         character(len=18), save :: class = 'sim_beams:'
         integer :: i, n
         
         call this%err%werrfl2(class//sname//' started')

         n = size(this%beam)

         do i = 1, n
            call this%beam(i)%del()
         end do

         call this%err%werrfl2(class//sname//' ended')

      end subroutine end_sim_beams
!
      subroutine init_sim_species(this,input,fields)

         implicit none

         class(sim_species), intent(inout) :: this
         type(input_json), pointer, intent(inout) :: input
         class(sim_fields), intent(inout) :: fields
! local data
         character(len=18), save :: class = 'sim_species:'
         character(len=18), save :: sname = 'init_sim_species:'
         integer :: i,n,ndump
         integer :: npf,npx,npy,npz
         real, dimension(3,100) :: arg
         character(len=20) :: sn,s1
         integer :: indx, indy, indz, npmax
         real :: min, max, cwp, n0
         real :: qm, qbm, dt, dx, dy, dz

         this%err => input%err
         this%p => input%pp
         this%sp3 => input%sp
         this%sp2 => input%sp
         
         call this%err%werrfl2(class//sname//' started')

         call input%get('simulation.n0',n0)
         call input%get('simulation.indx',indx)
         call input%get('simulation.indy',indy)
         call input%get('simulation.indz',indz)

         cwp=5.32150254*1e9/sqrt(n0)
         call input%get('simulation.box.x(1)',min)
         call input%get('simulation.box.x(2)',max)
         dx=(max-min)/cwp/real(2**indx)
         call input%get('simulation.box.y(1)',min)
         call input%get('simulation.box.y(2)',max)
         dy=(max-min)/cwp/real(2**indy)
         call input%get('simulation.box.z(1)',min)
         call input%get('simulation.box.z(2)',max)
         dz=(max-min)/cwp/real(2**indz)

         call input%get('simulation.nspecies',n)

         allocate(this%spe(n),this%pf(n))
         
         do i = 1, n

            write (sn,'(I3.3)') i
            s1 = 'species('//trim(sn)//')'
            call input%get(trim(s1)//'.profile',npf)
            call input%get(trim(s1)//'.np(1)',npx)
            call input%get(trim(s1)//'.np(2)',npy)

            call this%pf(i)%new(this%p,this%err,this%sp2,npf=npf,npx=npx,&
            &npy=npy,arg=arg(1:2,:))

            call input%get(trim(s1)//'.q',qm)
            call input%get(trim(s1)//'.m',qbm)
            qbm = qm/qbm
            npmax = npx*npy/this%p%getlnvp()
            qm = qm/abs(qm)/(real(npx)/2**indx)/(real(npy)/2**indy)
            call this%spe(i)%new(this%p,this%err,this%sp2,this%pf(i),fields%qe0(i),qm=qm,&
            &qbm=qbm,dt=dz,ci=1.0,xdim=7,npmax=npmax,nbmax=int(0.01*npmax))

         end do
         call this%err%werrfl2(class//sname//' ended')
      end subroutine init_sim_species
!
      subroutine end_sim_species(this)

         implicit none

         class(sim_species), intent(inout) :: this
! local data
         character(len=18), save :: sname = 'end_sim_species:'
         character(len=18), save :: class = 'sim_species:'
         integer :: i, n
         
         call this%err%werrfl2(class//sname//' started')

         n = size(this%spe)

         do i = 1, n
            call this%spe(i)%del()
         end do

         call this%err%werrfl2(class//sname//' ended')

      end subroutine end_sim_species
!
      subroutine go_simulation(this)

         implicit none

         class(simulation), intent(inout) :: this
! local data
         character(len=18), save :: sname = 'go_simulation:'
         integer, dimension(10) :: istat
         integer :: i,j,l,m,n,ierr

         call this%err%werrfl2(class//sname//' started')

         do i = this%start3d, this%nstep3d
   
            write (erstr,*) '3D step:', i        
            call this%err%werrfl0(erstr)
            
            call MPI_WAIT(this%id(1),istat,ierr)
            call MPI_WAIT(this%id(3),istat,ierr)
            call this%fields%qeb%as(0.0)
            do m = 1, this%nbeams
               call this%beams%beam(m)%qdp(this%fields%qeb)
            end do
            this%tag(1) = ntag()
            call this%fields%qeb%ag(this%tag(1),this%tag(1),this%id(1))
            this%tag(1) = ntag()
            call this%fields%qeb%pcg(this%tag(1),this%tag(1),this%id(2),this%id(3))    
   
            do l =  1, this%nspecies
               this%tag_spe(l) = ntag()
               call this%species%spe(l)%precv(this%fields%qe0(l),this%tag_spe(l))
            end do
            this%tag(2) = ntag()
            call this%fields%qi%precv(this%tag(2))
            this%tag(3) = ntag()
            call this%fields%cu%precv(this%tag(3))
            this%tag(4) = ntag()
            call this%fields%epw%precv(this%tag(4))
            this%tag(5) = ntag()
            call this%fields%fxy%precv(this%tag(5))
            this%tag(6) = ntag()
            call this%fields%psit%precv(this%tag(6))
            this%tag(7) = ntag()
            call this%fields%bxyz%precv(this%tag(7))
   
            call this%fields%fxy%cb(this%fields%bexyz,1,(/1,2/),(/1,2/))         
            call this%fields%psit%cb(this%fields%bexyz,1,(/1/),(/3/))
            call this%fields%bxyz%cb(this%fields%bbxyz,1,(/1,2,3/),(/1,2,3/))
   
   
            if (this%p%getstageid() == 0) then
               do m = 1, this%nspecies
                  call this%fields%qe0(m)%as(0.0)           
                  call this%species%spe(m)%qdp(this%fields%qe0(m))
                  call this%fields%qi%add(this%fields%qi,this%fields%qe0(m))
                  call this%fields%qe0(m)%cb(this%fields%qep(m),1,(/1/),(/1/))         
               end do
            end if
                     
            do j = 1, this%nstep2d
               write (erstr,*) '2D step:', j
               call this%err%werrfl0(erstr)
               if (j == this%nstep2d) then
                  call MPI_WAIT(this%id(2),istat,ierr)
               endif
               call this%fields%qb%cp(this%fields%qeb,j+1,(/1/),(/1/))
               call this%fields%qb%fftrk(1)
               call this%fields%qb%elf(this%fields%epwb)
               call this%fields%qe%mult(this%fields%qi,-1.0)
               do l = 1, this%nspecies
                  call this%fields%qe0(l)%as(0.0)
                  call this%species%spe(l)%qdp(this%fields%qe0(l))                     
                  call this%fields%qe%add(this%fields%qe,this%fields%qe0(l))
               end do
               call this%fields%qe%fftrk(1)
               call this%fields%qe%pot(this%fields%psi)
               call this%fields%psi%grad(this%fields%fxy)
               call this%fields%fxy%fftkr(1)
               call this%fields%psi%fftkr(1)
               do l = 1, this%nspecies            
                  call this%species%spe(l)%extpsi(this%fields%psi,this%dex)
               end do
               call this%fields%psi%mult(this%fields%psi,this%dex*this%dex)
               if (allocated(this%fields%psi3d)) call this%fields%psi%cb(this%fields%psi3d,j+1,(/1/),(/1/))
               call this%fields%cu%div(this%fields%div_vpot)
               call this%fields%div_vpot%pot(this%fields%psit)
               call this%fields%psit%mult(this%fields%psit,-this%dex)
               call this%fields%psit%fftkr(1)
               call this%fields%cu%bf(this%fields%bxyz)
               do l = 1, this%iter
                  call this%fields%bxyz%mult(this%fields%epw,(/1,2/),(/2,1/),(/-this%dex,this%dex/))
                  call this%fields%bxyz%mult(this%fields%bxyz,(/3/),(/3/),(/this%dex/))
                  call this%fields%bxyz%sub(this%fields%bxyz,this%fields%epwb,(/1/),(/1/),(/2/))
                  call this%fields%bxyz%add(this%fields%bxyz,this%fields%epwb,(/2/),(/2/),(/1/))
                  call this%fields%bxyz%fftkr(2)
                  call this%fields%cu%as(0.0)
                  call this%fields%dcu%as(0.0)
                  call this%fields%amu%as(0.0)
                  do m = 1, this%nspecies
                     call this%fields%cu0(m)%as(0.0)
                     call this%fields%dcu0(m)%as(0.0)
                     call this%fields%amu0(m)%as(0.0)
                     call this%species%spe(m)%amjdp(this%fields%fxy,this%fields%bxyz,this%fields%psit,&
                     &this%fields%cu0(m),this%fields%amu0(m),this%fields%dcu0(m),this%dex)
                     call this%fields%cu0(m)%mult(this%fields%cu0(m),this%dex)
                     call this%fields%amu0(m)%mult(this%fields%amu0(m),this%dex)
                     call this%fields%dcu0(m)%mult(this%fields%dcu0(m),this%dex)
                     call this%fields%cu%add(this%fields%cu,this%fields%cu0(m))
                     call this%fields%amu%add(this%fields%amu,this%fields%amu0(m))
                     call this%fields%dcu%add(this%fields%dcu,this%fields%dcu0(m))
                  end do
                  if (l == this%iter) then
                     do m = 1, this%nspecies              
                        call this%fields%reg%add(this%fields%qe0(m),this%fields%cu0(m),(/1/),(/1/),(/3/))
                        call this%fields%reg%fftrk(1)
                        call this%fields%reg%smooth(this%fields%reg)
                        call this%fields%reg%fftkr(1)
                        call this%fields%qe0(m)%as(this%fields%reg)
                        call this%fields%qe0(m)%cb(this%fields%qep(m),j+1,(/1/),(/1/))                     
                     end do
                     if (allocated(this%fields%cu3d)) then
                        call this%fields%cu%cb(this%fields%cu3d,j+1,(/1,2,3/),(/1,2,3/))
                     end if               
                  endif
                  call this%fields%cu%fftrk(1)
                  call this%fields%dcu%fftrk(1)
                  call this%fields%amu%fftrk(3)
                  call this%fields%cu%bfqp(this%fields%dcu,this%fields%amu,this%fields%epw,this%dex2,&
                  &this%dex)
                  call this%fields%cu%div(this%fields%div_vpot)
                  call this%fields%div_vpot%pot(this%fields%psit)
                  call this%fields%psit%mult(this%fields%psit,-this%dex)
                  call this%fields%psit%fftkr(1)
                  call this%fields%cu%bf(this%fields%bxyz)
               enddo
               call this%fields%bxyz%mult(this%fields%epw,(/1,2/),(/2,1/),(/-this%dex,this%dex/))
               call this%fields%bxyz%mult(this%fields%bxyz,(/3/),(/3/),(/this%dex/))
               call this%fields%bxyz%sub(this%fields%bxyz,this%fields%epwb,(/1/),(/1/),(/2/))
               call this%fields%bxyz%add(this%fields%bxyz,this%fields%epwb,(/2/),(/2/),(/1/))
               call this%fields%bxyz%fftkr(2)
               call this%fields%fxy%mult(this%fields%fxy,-1.0)
               call this%fields%fxy%add(this%fields%fxy,this%fields%bxyz,(/1/),(/1/),(/2/))
               call this%fields%fxy%sub(this%fields%fxy,this%fields%bxyz,(/2/),(/2/),(/1/))
               call this%fields%dcu%mult(this%fields%dcu,this%dxi)
               call this%fields%cu%sub(this%fields%cu,this%fields%dcu,(/1,2/),(/1,2/),(/1,2/))
               do m = 1, this%nspecies              
                  call this%species%spe(m)%push(this%fields%fxy,this%fields%bxyz,this%fields%psit,&
                  &this%dex)
               end do
               call this%fields%fxy%mult(this%fields%fxy,this%dex)
               call this%fields%fxy%cb(this%fields%bexyz,j+1,(/1,2/),(/1,2/))
               call this%fields%psit%cb(this%fields%bexyz,j+1,(/1/),(/3/))
               call this%fields%bxyz%mult(this%fields%bxyz,(/1,2/),(/1,2/),(/this%dex,this%dex/))
               call this%fields%bxyz%cb(this%fields%bbxyz,j+1,(/1,2,3/),(/1,2,3/))
            enddo
            
            do m = 1, this%nspecies                       
               call this%species%spe(m)%psend(this%tag_spe(m),this%id_spe(m))
            end do
            call MPI_WAIT(this%id(4),istat,ierr)
            call this%fields%qi%psend(this%tag(2),this%id(4))
            call MPI_WAIT(this%id(5),istat,ierr)
            call this%fields%cu%psend(this%tag(3),this%id(5))
            call MPI_WAIT(this%id(6),istat,ierr)
            call this%fields%epw%psend(this%tag(4),this%id(6))
            call MPI_WAIT(this%id(7),istat,ierr)
            call this%fields%fxy%psend(this%tag(5),this%id(7))
            call MPI_WAIT(this%id(8),istat,ierr)
            call this%fields%psit%psend(this%tag(6),this%id(8))
            call MPI_WAIT(this%id(9),istat,ierr)
            call this%fields%bxyz%psend(this%tag(7),this%id(9))
   
            call this%diag()
   
            do m = 1, this%nbeams
               this%tag_beam(m) = ntag()
               call MPI_WAIT(this%id_beam(m),istat,ierr)
               call this%beams%beam(m)%push(this%fields%bexyz,this%fields%bbxyz,this%dex,this%dxi,&
               &this%tag_beam(m),this%tag_beam(m),this%id_beam(m))
            end do
            
            do m = 1, this%nspecies                       
               call MPI_WAIT(this%id_spe(m),istat,ierr)
               call this%species%spe(m)%renew(this%species%pf(m),this%fields%qe0(m))
            end do
                             
         enddo


         call this%err%werrfl2(class//sname//' ended')

      end subroutine go_simulation
!
      subroutine diag_simulation(this)

         implicit none

         class(simulation), intent(inout) :: this
! local data
         character(len=18), save :: sname = 'diag_simulation:'

         call this%err%werrfl2(class//sname//' started')


         call this%err%werrfl2(class//sname//' ended')

      end subroutine diag_simulation
!
      function ntag()
      
         implicit none
         integer, save :: tag = 0
         integer :: ntag
                 
         ntag = tag
         tag = tag + 1
      
      end function ntag
!                  
      end module simulation_class