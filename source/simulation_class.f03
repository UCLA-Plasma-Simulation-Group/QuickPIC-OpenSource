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
      type sim_diag

         private

         class(parallel_pipe),pointer :: p => null()
         class(perrors),pointer :: err => null()
         class(spect3d), pointer :: sp3 => null()
         class(spect2d), pointer :: sp2 => null()
         type(hdf5file) :: file
         class(*), pointer :: obj => null()
         integer, allocatable :: slice, slice_pos, psample, dim
         integer :: df

      end type sim_diag
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
         type(sim_diag), dimension(:), allocatable :: diag
         integer :: iter, nstep3d, nstep2d, start3d,nbeams,nspecies, tstep
         integer, dimension(8) :: tag
         integer, dimension(:), allocatable :: tag_spe, id_spe, id
         integer, dimension(:), allocatable :: tag_beam, id_beam
         real :: dex, dxi, dex2

         contains
         
         generic :: new => init_simulation
         generic :: del => end_simulation
         generic :: go => go_simulation

         procedure, private :: init_simulation, end_simulation
         procedure, private :: init_diag, diag_simulation
         procedure, private :: go_simulation

      end type simulation
!
      character(len=20), save :: class = 'simulation: '
      character(len=128), save :: erstr
            
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

         call this%init_diag()                 

         allocate(this%tag_spe(this%nspecies),this%tag_beam(this%nbeams))
         allocate(this%id_spe(this%nspecies),this%id_beam(this%nbeams))

         allocate(this%id(9+size(this%diag)))
         this%id(:) = MPI_REQUEST_NULL
         this%id_spe(:) = MPI_REQUEST_NULL
         this%id_beam(:) = MPI_REQUEST_NULL                 

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
         character(len=20) :: s1, s2, s3
         character(len=:), allocatable :: ff   
         integer :: i,n,ndump,j,k,l,m

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
         
         loop1: do i = 1, n
            write (s1, '(I4.4)') i
            call input%info('species('//trim(s1)//').diag',n_children=m)
            do j = 1, m
               write (s2, '(I4.4)') j
               call input%get('species('//trim(s1)//').diag'//'('//trim(s2)//').ndump',ndump)
               if (ndump>0) then
                  call input%info('species('//trim(s1)//').diag'//'('//trim(s2)//').name',n_children=l)
                  do k = 1, l
                     write (s3, '(I4.4)') k
                     if(allocated(ff)) deallocate(ff)
                     call input%get('species('//trim(s1)//').diag'//'('//trim(s2)//').name'&
                     &//'('//trim(s3)//')',ff)
                     if (ff == 'jx' .or. ff == 'jy' .or. ff == 'jz') then
                        allocate(this%cu3d)
                        call this%cu3d%new(this%p,this%err,this%sp3,dim=1)
                        exit loop1
                     end if
                  end do
               end if
            end do
         end do loop1

         call input%info('field.diag',n_children=n)
         
         loop2: do i = 1, n
            write (s1,'(I4.4)') i
            call input%get('field.diag('//trim(s1)//').ndump',ndump)
            if (ndump > 0) then
               call input%info('field.diag('//trim(s1)//').name',n_children=m)
               do j = 1, m
                  write (s2,'(I4.4)') j
                  if(allocated(ff)) deallocate(ff)
                  call input%get('field.diag('//trim(s1)//').name('//trim(s2)//')',ff)
                  if (ff == 'psi') then
                     allocate(this%psi3d)
                     call this%psi3d%new(this%p,this%err,this%sp3,dim=1)
                     exit loop2
                  end if
               end do
            end if
         end do loop2
         

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
               call file_rst%new(&
               &filename = './RST/Beam-'//trim(sn)//'/',&
               &dataname = 'RST-beam'//trim(sn)//'-'//trim(sid),&
               &n = rst_timestep)
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

            this%tstep = i
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

            call this%diag_simulation()
                             
         end do


         call this%err%werrfl2(class//sname//' ended')

      end subroutine go_simulation
!
      subroutine init_diag(this)

         implicit none

         class(simulation), intent(inout), target :: this
! local data
         character(len=18), save :: sname = 'init_diag:'
         integer :: n_diag = 0, ndump, slice, slice_pos, sample
         integer :: n, m, l, i, j, k, ii
         character(len=20) :: s1, s2, s3, s4, sn1, sn2, sn3, sn4
         character(len=:), allocatable :: ss,sl
         real :: min, max, cwp, n0, dt
         real :: alx1, aly1, alz1, alx2, aly2, alz2
         logical :: rst
         integer :: ierr, indx, indy, indz, dim
         class (*), pointer :: obj => null()


         call this%err%werrfl2(class//sname//' started')

         call this%in%get('simulation.n0',n0)
         call this%in%get('simulation.indx',indx)
         call this%in%get('simulation.indy',indy)
         call this%in%get('simulation.indz',indz)

         cwp=5.32150254*1e9/sqrt(n0)
         call this%in%get('simulation.box.x(1)',min)
         call this%in%get('simulation.box.x(2)',max)
         alx1 = min/cwp
         alx2 = max/cwp
         call this%in%get('simulation.box.y(1)',min)
         call this%in%get('simulation.box.y(2)',max)
         aly1 = min/cwp 
         aly2 = max/cwp 
         call this%in%get('simulation.box.z(1)',min)
         call this%in%get('simulation.box.z(2)',max)
         alz1 = min/cwp 
         alz2 = max/cwp 
         call this%in%get('simulation.dt',dt)

         do i = 1, this%nbeams
            write (s1, '(I4.4)') i
            call this%in%info('beam('//trim(s1)//').diag',n_children=m)
            do j = 1, m
               write (s2, '(I4.4)') j
               call this%in%get('beam('//trim(s1)//').diag'//'('//trim(s2)//').ndump',ndump)
               if (ndump>0) then
                  call this%in%info('beam('//trim(s1)//').diag'//'('//trim(s2)//').name',n_children=l)
                  if(this%in%found('beam('//trim(s1)//').diag'//'('//trim(s2)//').slice')) then
                     call this%in%info('beam('//trim(s1)//').diag'//'('//trim(s2)//').slice',n_children=n)
                     n_diag = n_diag + l*n
                  else
                     n_diag = n_diag + l
                  end if
               end if
            end do
         end do

         do i = 1, this%nspecies
            write (s1, '(I4.4)') i
            call this%in%info('species('//trim(s1)//').diag',n_children=m)
            do j = 1, m
               write (s2, '(I4.4)') j
               call this%in%get('species('//trim(s1)//').diag'//'('//trim(s2)//').ndump',ndump)
               if (ndump>0) then
                  call this%in%info('species('//trim(s1)//').diag'//'('//trim(s2)//').name',n_children=l)
                  if(this%in%found('species('//trim(s1)//').diag'//'('//trim(s2)//').slice')) then
                     call this%in%info('species('//trim(s1)//').diag'//'('//trim(s2)//').slice',n_children=n)
                     n_diag = n_diag + l*n
                  else
                     n_diag = n_diag + l
                  end if
               end if
            end do
         end do

         call this%in%info('field.diag',n_children=n)
         do i = 1, n
            write (s1, '(I4.4)') i
            call this%in%get('field.diag('//trim(s1)//').ndump',ndump)
            if (ndump>0) then
               call this%in%info('field.diag('//trim(s1)//').name',n_children=l)
               if(this%in%found('field.diag('//trim(s1)//').slice')) then
                  call this%in%info('field.diag('//trim(s1)//').slice',n_children=m)
                  n_diag = n_diag + l*m
               else
                  n_diag = n_diag + l
               end if
            end if
         end do

         call this%in%get('simulation.dump_restart',rst)
         if (rst) n_diag = n_diag + this%nbeams

         allocate(this%diag(n_diag))
         n_diag = 0
         do i = 1, this%nbeams
            write (s1, '(I4.4)') i
            call this%in%info('beam('//trim(s1)//').diag',n_children=m)
            do j = 1, m
               write (s2, '(I4.4)') j
               call this%in%get('beam('//trim(s1)//').diag'//'('//trim(s2)//').ndump',ndump)
               if (ndump>0) then
                  call this%in%info('beam('//trim(s1)//').diag'//'('//trim(s2)//').name',n_children=l)
                  do k = 1, l
                     write (s3, '(I4.4)') k
                     if (allocated(ss)) deallocate(ss)
                     call this%in%get('beam('//trim(s1)//').diag'//'('//trim(s2)//').name'&
                     &//'('//trim(s3)//')',ss)
                     select case (trim(ss))
                     case ('charge')
                        if (this%in%found('beam('//trim(s1)//').diag'//'('//trim(s2)//').slice')) then
                           call this%in%info('beam('//trim(s1)//').diag'//'('//trim(s2)//').slice',n_children=n)
                           do ii = 1, n
                              n_diag = n_diag + 1 
                              this%diag(n_diag)%df = ndump
                              this%diag(n_diag)%obj => this%fields%qeb
                              allocate(this%diag(n_diag)%slice,this%diag(n_diag)%slice_pos)
                              allocate(this%diag(n_diag)%dim)
                              this%diag(n_diag)%dim = 1
                              if (allocated(sl)) deallocate(sl)
                              write (s4, '(I4.4)') ii
                              call this%in%get('beam('//trim(s1)//').diag'//'('//trim(s2)//').&
                              &slice('//trim(s4)//').(1)',sl)
                              select case (sl)
                              case ('yz')
                                 this%diag(n_diag)%slice = 1
                                 call this%diag(n_diag)%file%new(&
                                 &axismin = (/aly1,alz1,alx1/),&
                                 &axismax = (/aly2,alz2,alx2/),&
                                 &axisname  = (/'y  ','\xi','x  '/),&
                                 &axislabel = (/'y  ','\xi','x  '/))
                              case ('xz')
                                 this%diag(n_diag)%slice = 2
                                 call this%diag(n_diag)%file%new(&
                                 &axismin = (/alx1,alz1,aly1/),&
                                 &axismax = (/alx2,alz2,aly2/),&
                                 &axisname  = (/'x  ','\xi','y  '/),&
                                 &axislabel = (/'x  ','\xi','y  '/))
                              case ('xy')
                                 this%diag(n_diag)%slice = 3
                                 call this%diag(n_diag)%file%new(&
                                 &axismin = (/alx1,aly1,alz1/),&
                                 &axismax = (/alx2,aly2,alz2/),&
                                 &axisname  = (/'x','y','z'/),&
                                 &axislabel = (/'x','y','z'/))
                              end select
                              call this%in%get('beam('//trim(s1)//').diag'//'('//trim(s2)//').&
                              &slice('//trim(s4)//').(2)',this%diag(n_diag)%slice_pos)
                              if (this%p%getidproc() == 0) then
                                 call system('mkdir -p ./ChargeDensity/beam'//trim(s1)//'_slice'//trim(s4)//'/')
                              end if
                              call this%diag(n_diag)%file%new(&
                              &timeunits = '1 / \omega_p',&
                              &dt = dt,&
                              &axisunits = (/'c / \omega_p','c / \omega_p','c / \omega_p'/),&
                              &rank = 2,&
                              &filename = './ChargeDensity/beam'//trim(s1)//'_slice'//trim(s4)//'/',&
                              &dataname = 'chargeslice'//sl,&
                              &units = 'n_0',&
                              &label = 'Charge Density')
                           end do
                        else
                           n_diag = n_diag + 1 
                           this%diag(n_diag)%df = ndump
                           this%diag(n_diag)%obj => this%fields%qeb
                           allocate(this%diag(n_diag)%dim)
                           this%diag(n_diag)%dim = 1
                           if (this%p%getidproc() == 0) then
                              call system('mkdir -p ./ChargeDensity/beam'//trim(s1)//'/')
                           end if
                           call this%diag(n_diag)%file%new(&
                           &timeunits = '1 / \omega_p',&
                           &dt = dt,&
                           &axisname  = (/'x  ','y  ','\xi'/),&
                           &axislabel = (/'x  ','y  ','\xi'/),&
                           &axisunits = (/'c / \omega_p','c / \omega_p','c / \omega_p'/),&
                           &axismin = (/alx1,aly1,alz1/),&
                           &axismax = (/alx2,aly2,alz2/),&
                           &rank = 3,&
                           &filename = './ChargeDensity/beam'//trim(s1)//'/',&
                           &dataname = 'charge',&
                           &units = 'n_0',&
                           &label = 'Charge Density')
                        end if
                     case ('raw')
                        n_diag = n_diag + 1 
                        this%diag(n_diag)%df = ndump
                        this%diag(n_diag)%obj => this%beams%beam(i)
                        allocate(this%diag(n_diag)%psample)
                        call this%in%get('beam('//trim(s1)//').diag'//'('//trim(s2)//').&
                        &sample',this%diag(n_diag)%psample)
                        if (this%p%getidproc() == 0) then
                           call system('mkdir -p ./Raw/beam'//trim(s1)//'/')
                        end if
                        call this%diag(n_diag)%file%new(&
                        &timeunits = '1 / \omega_p',&
                        &dt = dt,&
                        &ty = 'particles',&
                        &filename = './Raw/beam'//trim(s1)//'/',&
                        &dataname = 'raw',&
                        &units = '',&
                        &label = 'Beam RAw')
                     end select
                  end do
               end if
            end do
         end do

         do i = 1, this%nspecies
            write (s1, '(I4.4)') i
            call this%in%info('species('//trim(s1)//').diag',n_children=m)
            do j = 1, m
               write (s2, '(I4.4)') j
               call this%in%get('species('//trim(s1)//').diag'//'('//trim(s2)//').ndump',ndump)
               if (ndump>0) then
                  call this%in%info('species('//trim(s1)//').diag'//'('//trim(s2)//').name',n_children=l)
                  do k = 1, l
                     write (s3, '(I4.4)') k
                     if (allocated(ss)) deallocate(ss)
                     call this%in%get('species('//trim(s1)//').diag'//'('//trim(s2)//').name'&
                     &//'('//trim(s3)//')',ss)
                     select case (trim(ss))
                     case ('charge')
                        sn1 = 'ChargeDensity'
                        sn2 = 'charge'
                        sn3 = 'n_0'
                        sn4 = 'Charge Density'
                        dim = 1
                        obj => this%fields%qep(i)
                     case ('jx')
                        sn1 = 'Jx'
                        sn2 = 'jx'
                        sn3 = 'n_0 c'
                        sn4 = 'J_x'
                        dim = 1
                        obj => this%fields%cu3d
                     case ('jy')
                        sn1 = 'Jy'
                        sn2 = 'jy'
                        sn3 = 'n_0 c'
                        sn4 = 'J_y'
                        dim = 2
                        obj => this%fields%cu3d
                     case ('jz')
                        sn1 = 'Jz'
                        sn2 = 'jz'
                        sn3 = 'n_0 c'
                        sn4 = 'J_z'
                        dim = 3
                        obj => this%fields%cu3d
                     end select
                     if (this%in%found('species('//trim(s1)//').diag'//'('//trim(s2)//').slice')) then
                        call this%in%info('species('//trim(s1)//').diag'//'('//trim(s2)//').slice',n_children=n)
                        do ii = 1, n
                           n_diag = n_diag + 1 
                           this%diag(n_diag)%df = ndump
                           this%diag(n_diag)%obj => obj
                           allocate(this%diag(n_diag)%slice,this%diag(n_diag)%slice_pos)
                           allocate(this%diag(n_diag)%dim)
                           this%diag(n_diag)%dim = dim
                           if (allocated(sl)) deallocate(sl)
                           write (s4, '(I4.4)') ii
                           call this%in%get('species('//trim(s1)//').diag'//'('//trim(s2)//').&
                           &slice('//trim(s4)//').(1)',sl)
                           select case (sl)
                           case ('yz')
                              this%diag(n_diag)%slice = 1
                              call this%diag(n_diag)%file%new(&
                              &axismin = (/aly1,alz1,alx1/),&
                              &axismax = (/aly2,alz2,alx2/),&
                              &axisname  = (/'y  ','\xi','x  '/),&
                              &axislabel = (/'y  ','\xi','x  '/))
                           case ('xz')
                              this%diag(n_diag)%slice = 2
                              call this%diag(n_diag)%file%new(&
                              &axismin = (/alx1,alz1,aly1/),&
                              &axismax = (/alx2,alz2,aly2/),&
                              &axisname  = (/'x  ','\xi','y  '/),&
                              &axislabel = (/'x  ','\xi','y  '/))
                           case ('xy')
                              this%diag(n_diag)%slice = 3
                              call this%diag(n_diag)%file%new(&
                              &axismin = (/alx1,aly1,alz1/),&
                              &axismax = (/alx2,aly2,alz2/),&
                              &axisname  = (/'x','y','z'/),&
                              &axislabel = (/'x','y','z'/))
                           end select
                           call this%in%get('species('//trim(s1)//').diag'//'('//trim(s2)//').&
                           &slice('//trim(s4)//').(2)',this%diag(n_diag)%slice_pos)
                           if (this%p%getidproc() == 0) then
                              call system('mkdir -p ./'//trim(sn1)//'/species'//trim(s1)//'_slice'//trim(s4)//'/')
                           end if
                           call this%diag(n_diag)%file%new(&
                           &timeunits = '1 / \omega_p',&
                           &dt = dt,&
                           &axisunits = (/'c / \omega_p','c / \omega_p','c / \omega_p'/),&
                           &rank = 2,&
                           &filename = './'//trim(sn1)//'/species'//trim(s1)//'_slice'//trim(s4)//'/',&
                           &dataname = trim(sn2)//'slice'//sl,&
                           &units = trim(sn3),&
                           &label = trim(sn4))
                        end do
                     else
                        n_diag = n_diag + 1 
                        this%diag(n_diag)%df = ndump
                        this%diag(n_diag)%obj => obj
                        allocate(this%diag(n_diag)%dim)
                        this%diag(n_diag)%dim = dim
                        if (this%p%getidproc() == 0) then
                           call system('mkdir -p ./'//trim(sn1)//'/species'//trim(s1)//'/')
                        end if
                        call this%diag(n_diag)%file%new(&
                        &timeunits = '1 / \omega_p',&
                        &dt = dt,&
                        &axisname  = (/'x  ','y  ','\xi'/),&
                        &axislabel = (/'x  ','y  ','\xi'/),&
                        &axisunits = (/'c / \omega_p','c / \omega_p','c / \omega_p'/),&
                        &axismin = (/alx1,aly1,alz1/),&
                        &axismax = (/alx2,aly2,alz2/),&
                        &rank = 3,&
                        &filename = './'//trim(sn1)//'/species'//trim(s1)//'/',&
                        &dataname = trim(sn2),&
                        &units = trim(sn3),&
                        &label = trim(sn4))
                     end if
                  end do
               end if
            end do
         end do

         call this%in%info('field.diag',n_children=n)
         do i = 1, n
            write (s1, '(I4.4)') i
            call this%in%get('field.diag('//trim(s1)//').ndump',ndump)
            if (ndump>0) then
               call this%in%info('field.diag('//trim(s1)//').name',n_children=m)
               do j = 1, m
                  write (s2, '(I4.4)') j
                  if (allocated(ss)) deallocate(ss)
                  call this%in%get('field.diag('//trim(s1)//').name('//trim(s2)//')',ss)
                  select case (trim(ss))
                  case ('ex')
                     sn1 = 'Ex'
                     sn2 = 'ex'
                     sn3 = 'mc\omega_p/e'
                     sn4 = 'Electric Field'
                     dim = 1
                     obj => this%fields%bexyz
                  case ('ey')
                     sn1 = 'Ey'
                     sn2 = 'ey'
                     sn3 = 'mc\omega_p/e'
                     sn4 = 'Electric Field'
                     dim = 2
                     obj => this%fields%bexyz
                  case ('ez')
                     sn1 = 'Ez'
                     sn2 = 'ez'
                     sn3 = 'mc\omega_p/e'
                     sn4 = 'Electric Field'
                     dim = 3
                     obj => this%fields%bexyz
                  case ('bx')
                     sn1 = 'Bx'
                     sn2 = 'bx'
                     sn3 = 'mc\omega_p/e'
                     sn4 = 'Magnetic Field'
                     dim = 1
                     obj => this%fields%bbxyz
                  case ('by')
                     sn1 = 'By'
                     sn2 = 'by'
                     sn3 = 'mc\omega_p/e'
                     sn4 = 'Magnetic Field'
                     dim = 2
                     obj => this%fields%bbxyz
                  case ('bz')
                     sn1 = 'Bz'
                     sn2 = 'bz'
                     sn3 = 'mc\omega_p/e'
                     sn4 = 'Magnetic Field'
                     dim = 3
                     obj => this%fields%bbxyz
                  case ('psi')
                     sn1 = 'Psi'
                     sn2 = 'psi'
                     sn3 = 'mc^2'
                     sn4 = '\Psi'
                     dim = 1
                     obj => this%fields%psi3d
                  end select
                  if (this%in%found('field.diag('//trim(s1)//').slice')) then
                     call this%in%info('field.diag('//trim(s1)//').slice',n_children=l)
                     do k = 1, l
                        n_diag = n_diag + 1 
                        this%diag(n_diag)%df = ndump
                        this%diag(n_diag)%obj => obj
                        allocate(this%diag(n_diag)%slice,this%diag(n_diag)%slice_pos)
                        allocate(this%diag(n_diag)%dim)
                        this%diag(n_diag)%dim = dim
                        if (allocated(sl)) deallocate(sl)
                        write (s3, '(I4.4)') k
                        call this%in%get('field.diag('//trim(s1)//').&
                        &slice('//trim(s3)//').(1)',sl)
                        select case (sl)
                        case ('yz')
                           this%diag(n_diag)%slice = 1
                           call this%diag(n_diag)%file%new(&
                           &axismin = (/aly1,alz1,alx1/),&
                           &axismax = (/aly2,alz2,alx2/),&
                           &axisname  = (/'y  ','\xi','x  '/),&
                           &axislabel = (/'y  ','\xi','x  '/))
                        case ('xz')
                           this%diag(n_diag)%slice = 2
                           call this%diag(n_diag)%file%new(&
                           &axismin = (/alx1,alz1,aly1/),&
                           &axismax = (/alx2,alz2,aly2/),&
                           &axisname  = (/'x  ','\xi','y  '/),&
                           &axislabel = (/'x  ','\xi','y  '/))
                        case ('xy')
                           this%diag(n_diag)%slice = 3
                           call this%diag(n_diag)%file%new(&
                           &axismin = (/alx1,aly1,alz1/),&
                           &axismax = (/alx2,aly2,alz2/),&
                           &axisname  = (/'x','y','z'/),&
                           &axislabel = (/'x','y','z'/))
                        end select
                        call this%in%get('field.diag('//trim(s1)//').&
                        &slice('//trim(s3)//').(2)',this%diag(n_diag)%slice_pos)
                        if (this%p%getidproc() == 0) then
                           call system('mkdir -p ./Fields/'//trim(sn1)//'_slice'//trim(s3)//'/')
                        end if
                        call this%diag(n_diag)%file%new(&
                        &timeunits = '1 / \omega_p',&
                        &dt = dt,&
                        &axisunits = (/'c / \omega_p','c / \omega_p','c / \omega_p'/),&
                        &rank = 2,&
                        &filename = './Fields/'//trim(sn1)//'_slice'//trim(s3)//'/',&
                        &dataname = trim(sn2)//'slice'//sl,&
                        &units = trim(sn3),&
                        &label = trim(sn4))
                     end do
                  else
                     n_diag = n_diag + 1 
                     this%diag(n_diag)%df = ndump
                     this%diag(n_diag)%obj => obj
                     allocate(this%diag(n_diag)%dim)
                     this%diag(n_diag)%dim = dim
                     if (this%p%getidproc() == 0) then
                        call system('mkdir -p ./Fields/'//trim(sn1)//'/')
                     end if
                     call this%diag(n_diag)%file%new(&
                     &timeunits = '1 / \omega_p',&
                     &dt = dt,&
                     &axisname  = (/'x  ','y  ','\xi'/),&
                     &axislabel = (/'x  ','y  ','\xi'/),&
                     &axisunits = (/'c / \omega_p','c / \omega_p','c / \omega_p'/),&
                     &axismin = (/alx1,aly1,alz1/),&
                     &axismax = (/alx2,aly2,alz2/),&
                     &rank = 3,&
                     &filename = './Fields/'//trim(sn1)//'/',&
                     &dataname = trim(sn2),&
                     &units = trim(sn3),&
                     &label = trim(sn4))
                  end if
               end do
            end if
         end do

         if (rst) then
            call this%in%get('simulation.ndump_restart',ndump) 
            do i = 1, this%nbeams
               write (s1,'(I4.4)') i
               write (s2,'(I10.10)') this%p%getidproc()
               n_diag = n_diag + 1
               this%diag(n_diag)%df = ndump
               this%diag(n_diag)%obj => this%beams%beam(i)
               call this%diag(n_diag)%file%new(&
               &ty='rst',&
               &filename = './RST/Beam-'//trim(s1)//'/',&
               &dataname = 'RST-beam'//trim(s1)//'-'//trim(s2))
            end do
         end if

         call MPI_BARRIER(this%p%getlworld(),ierr)

         call this%err%werrfl2(class//sname//' ended')

      end subroutine init_diag
!
      subroutine diag_simulation(this)

         implicit none

         class(simulation), intent(inout) :: this
! local data
         character(len=18), save :: sname = 'diag_simulation:'
         integer :: n, m, l, i, j, k, ierr, idn = 9
         integer, dimension(10) :: istat
         real :: dt

         call this%err%werrfl2(class//sname//' started')

         call this%in%get('simulation.dt',dt)

         n = size(this%diag)

         do i = 1, n
            if (mod(this%tstep-1,this%diag(i)%df) == 0) then
               call this%diag(i)%file%new(n = this%tstep, t = this%tstep*dt)
               select type (obj => this%diag(i)%obj)
               type is (field3d)
                  if (allocated(this%diag(i)%slice)) then
                     this%tag(1) = ntag()
                     call MPI_WAIT(this%id(idn+i),istat,ierr)
                     call obj%wr(this%diag(i)%file,this%diag(i)%dim,this%diag(i)%slice,this%diag(i)%slice_pos,&
                     &this%tag(1),this%tag(1),this%id(idn+i))
                  else
                     this%tag(1) = ntag()
                     call MPI_WAIT(this%id(idn+i),istat,ierr)
                     call obj%wr(this%diag(i)%file,this%diag(i)%dim,this%tag(1),this%tag(1),this%id(idn+i))
                  end if
               type is (beam3d)
                  this%tag(1) = ntag()
                  call MPI_WAIT(this%id(idn+i),istat,ierr)
                  call obj%wr(this%diag(i)%file,this%diag(i)%psample,(/this%dex,this%dex,this%dxi/),&
                  &this%tag(1),this%tag(1),this%id(idn+i))         
               end select
            end if
         end do

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