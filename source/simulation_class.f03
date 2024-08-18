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

      type fdist2d_wrap
         class(fdist2d), allocatable :: p
      end type fdist2d_wrap
!
      type fdist3d_wrap
         class(fdist3d), allocatable :: p
      end type fdist3d_wrap
!
      type sim_fields

         private

         class(parallel_pipe),pointer :: p => null()
         class(perrors),pointer :: err => null()
         class(spect3d), pointer :: sp3 => null()
         class(spect2d), pointer :: sp2 => null()
         type(field2d), allocatable :: qb, qe, qe1, psit, psi, div_vpot, reg
         !type(field2d), allocatable :: fxy, bxyz, cu, dcu, amu, epw, epwb,vpot 
         type(field2d), allocatable :: fxy, bxyz, bxyzfd, cu, dcu, amu, dccu, damu, epw, epwb, vpot        
         type(field3d), allocatable :: bexyz, bbxyz
         type(field3d), allocatable :: psi3d,cu3d,vpot3d

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
         type(fdist3d_wrap), dimension(:), allocatable :: pf

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
         type(fdist2d_wrap), dimension(:), allocatable :: pf
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
         character(len=128) :: path2

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
         type(sim_diag), dimension(:), allocatable :: diag, diag2
         integer :: iter, nstep3d, nstep2d, start3d,nbeams,nspecies, tstep, tstep2, noffstep2d
         integer, dimension(8) :: tag
         integer, dimension(:), allocatable :: tag_spe, id_spe, id
         integer, dimension(:), allocatable :: tag_beam, id_beam
         integer, dimension(:,:), allocatable :: id_bq, tag_bq
         !real :: dex, dxi, dex2, dt
         real :: dex, dxi, dex2, dt, Bx, By, Bz, Ex, Ey, Ez, bxp, byp

         contains

         generic :: new => init_simulation
         generic :: del => end_simulation
         generic :: go => go_simulation_2
         ! generic :: go=> go_simulation_1

         procedure, private :: init_simulation, end_simulation
         procedure, private :: init_diag, diag_simulation
         procedure, private :: go_simulation_1, go_simulation_2
         procedure, private :: diag2_simulation

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
         real :: min, max, n0, dx, dy, dz, dt
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

         call this%in%get('simulation.box.x(1)',min)
         call this%in%get('simulation.box.x(2)',max)
         dx=(max-min)/real(2**indx)
         call this%in%get('simulation.box.y(1)',min)
         call this%in%get('simulation.box.y(2)',max)
         dy=(max-min)/real(2**indy)
         call this%in%get('simulation.box.z(1)',min)
         call this%in%get('simulation.box.z(2)',max)
         dz=(max-min)/real(2**indz)

         this%dex = dx
         this%dxi = dz
         this%dex2 = dx * dx
         this%nstep2d = 2**this%sp3%getindz()/this%p%getnstage()
         this%noffstep2d = this%p%getstageid()*this%nstep2d
         call this%in%get('simulation.time',min)
         call this%in%get('simulation.dt',max)
         this%nstep3d = min/max
         call this%in%get('simulation.read_restart',read_rst)
         call this%in%get('simulation.dt',dt)
         this%dt = dt
         if (read_rst) then
            call this%in%get('simulation.restart_timestep',this%start3d)
            this%start3d = this%start3d + 1
         else
            this%start3d = 1
         end if
         call this%in%get('simulation.iter',this%iter)
         call this%in%get('simulation.nbeams',this%nbeams)
         call this%in%get('simulation.nspecies',this%nspecies)

         call this%fields%new(this%in)
         call this%beams%new(this%in)
         call this%species%new(this%in,(this%start3d-1)*dt)

         call this%init_diag()

         allocate(this%tag_spe(this%nspecies),this%tag_beam(this%nbeams))
         allocate(this%id_spe(this%nspecies),this%id_beam(this%nbeams))
         allocate(this%id_bq(this%nbeams,3),this%tag_bq(this%nbeams,2))

         allocate(this%id(9+size(this%diag)))
         this%id(:) = MPI_REQUEST_NULL
         this%id_spe(:) = MPI_REQUEST_NULL
         this%id_beam(:) = MPI_REQUEST_NULL
         this%id_bq(:,:) = MPI_REQUEST_NULL

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
         write(*,*) input%sp%getdex(), "dexe"


         call this%err%werrfl2(class//sname//' started')

         allocate(this%qb, this%qe, this%qe1, this%psit, this%psi)
         !allocate(this%div_vpot, this%reg, this%fxy, this%bxyz, this%cu)
         !allocate(this%dcu, this%amu, this%epw, this%epwb)
         allocate(this%div_vpot, this%reg, this%fxy, this%bxyz, this%bxyzfd, this%cu, this%dccu)
         allocate(this%dcu, this%amu, this%damu, this%epw, this%epwb)
         allocate(this%bexyz, this%bbxyz)
         

         call this%bexyz%new(this%p,this%err,this%sp3,dim=3)
         call this%bbxyz%new(this%p,this%err,this%sp3,dim=3)

         call this%qb%new(this%p,this%err,this%sp2,dim=1,fftflag=.true.,gcells=1)
         call this%qe%new(this%p,this%err,this%sp2,dim=1,fftflag=.true.)
         call this%qe1%new(this%p,this%err,this%sp2,dim=1,fftflag=.true.)
         call this%psit%new(this%p,this%err,this%sp2,dim=1,fftflag=.true.)
         call this%div_vpot%new(this%p,this%err,this%sp2,dim=1,fftflag=.true.)
         call this%psi%new(this%p,this%err,this%sp2,dim=1,fftflag=.true.)
         call this%reg%new(this%p,this%err,this%sp2,dim=1,fftflag=.true.)
         call this%fxy%new(this%p,this%err,this%sp2,dim=2,fftflag=.true.)
         call this%cu%new(this%p,this%err,this%sp2,dim=3,fftflag=.true.)
         call this%dcu%new(this%p,this%err,this%sp2,dim=2,fftflag=.true.)
         call this%amu%new(this%p,this%err,this%sp2,dim=3,fftflag=.true.)
         call this%epw%new(this%p,this%err,this%sp2,dim=2,fftflag=.true.,state=1) !state=1 meaning in kspace.
         call this%epwb%new(this%p,this%err,this%sp2,dim=2,fftflag=.true.)
         call this%dccu%new(this%p,this%err,this%sp2,dim=2,fftflag=.true.)
         !call this%damu%new(this%p,this%err,this%sp2,dim=3,fftflag=.true.)
         !write(*,*) this%sp2%getdex(), "dexd"
         call this%damu%new_fd(this%p,this%err,this%sp2,dim=2,fftflag=.true.,qe1=this%qe1,psi=this%psi,algorithm=1)
         ! call this%damu%new(this%p,this%err,this%sp2,dim=3,fftflag=.true.)
         call this%bxyz%new(this%p,this%err,this%sp2,dim=3,fftflag=.true.)
         call this%bxyzfd%new(this%p,this%err,this%sp2,dim=2,fftflag=.false.)
         ! call this%bxyzfd%new_fd(this%p,this%err,this%sp2,dim=2,fftflag=.false.,qe1=this%qe1,psi=this%psi,algorithm=1)
         
         call input%get('simulation.nspecies',n)

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
                        call this%cu3d%new(this%p,this%err,this%sp3,dim=3)
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
         
         loop3: do i = 1, n
            write (s1,'(I4.4)') i
            call input%get('field.diag('//trim(s1)//').ndump',ndump)
            if (ndump > 0) then
               call input%info('field.diag('//trim(s1)//').name',n_children=m)
               do j = 1, m
                  write (s2,'(I4.4)') j
                  if(allocated(ff)) deallocate(ff)
                  call input%get('field.diag('//trim(s1)//').name('//trim(s2)//')',ff)
                  if (ff == 'vpotx' .or. ff == 'vpoty' .or. ff == 'vpotz' ) then
                     allocate(this%vpot,this%vpot3d)
                     call this%vpot%new(this%p,this%err,this%sp2,dim=3,fftflag=.true.)
                     call this%vpot3d%new(this%p,this%err,this%sp3,dim=3)
                     exit loop3
                  end if
               end do
            end if
         end do loop3

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
         call this%qb%del()
         call this%qe%del()
         call this%qe1%del()
         call this%psit%del()
         call this%div_vpot%del()
         call this%psi%del()
         call this%reg%del()
         call this%fxy%del()
         call this%cu%del()
         call this%dcu%del()
         call this%amu%del()
         call this%dccu%del()
         call this%damu%del()         
         call this%epw%del()
         call this%epwb%del()
         call this%bxyz%del()
         call this%bxyzfd%del()

         if(allocated(this%psi3d)) call this%psi3d%del()
         if(allocated(this%cu3d)) call this%cu3d%del()

         call this%err%werrfl2(class//sname//' ended')

      end subroutine end_sim_fields
!
      subroutine init_sim_beams(this,input)

         implicit none

         class(sim_beams), intent(inout) :: this
         type(input_json), pointer, intent(inout) :: input
! local data
         character(len=18), save :: class = 'sim_beams:'
         character(len=18), save :: sname = 'init_sim_beams:'
         integer :: i,n
         integer :: npf,npx,npy,npz,npp
         real, dimension(3,100) :: arg
         logical :: quiet
         real :: gamma
         real :: qm, qbm, dt
         logical :: read_rst
         integer :: rst_timestep, ierr
         type(hdf5file) :: file_rst
         character(len=20) :: sn, sid, stime,s1

         this%err => input%err
         this%p => input%pp
         this%sp3 => input%sp
         this%sp2 => input%sp

         call this%err%werrfl2(class//sname//' started')

         call input%get('simulation.nbeams',n)

         allocate(this%beam(n),this%pf(n))

         do i = 1, n
            arg(:,:) = 0.0
            write (sn,'(I3.3)') i
            s1 = 'beam('//trim(sn)//')'
            call input%get(trim(s1)//'.profile',npf)
            select case (npf)
            case (0)
               allocate(fdist3d_000::this%pf(i)%p)
               call this%pf(i)%p%new(input,i)
            case (1)
               allocate(fdist3d_001::this%pf(i)%p)
               call this%pf(i)%p%new(input,i)
            case (2)
               allocate(fdist3d_002::this%pf(i)%p)
               call this%pf(i)%p%new(input,i)
            case (3)
               allocate(fdist3d_003::this%pf(i)%p)
               call this%pf(i)%p%new(input,i)
            case (4)
               allocate(fdist3d_004::this%pf(i)%p)
               call this%pf(i)%p%new(input,i)
            case (100)
               allocate(fdist3d_100::this%pf(i)%p)
               call this%pf(i)%p%new(input,i)
! Add new distributions right above this line
            case default
               write (erstr,*) 'Invalid beam profile number:', npf
               call this%err%equit(class//sname//erstr)
            end select

            call input%get(trim(s1)//'.q',qm)
            call input%get(trim(s1)//'.m',qbm)
            qbm = qm/qbm

            call input%get('simulation.dt',dt)

            call this%beam(i)%new(this%p,this%err,this%sp3,this%pf(i)%p,qbm=qbm,&
            &dt=dt,ci=1.0,xdim=7)

            call input%get('simulation.read_restart',read_rst)

            if (read_rst) then
               call input%get('simulation.restart_timestep',rst_timestep)
               write (sn,'(I4.4)') i
               write (sid,'(I10.10)') this%p%getidproc()
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
      subroutine init_sim_species(this,input,s)

         implicit none

         class(sim_species), intent(inout) :: this
         type(input_json), pointer, intent(inout) :: input
         real, intent(in) :: s
! local data
         character(len=18), save :: class = 'sim_species:'
         character(len=18), save :: sname = 'init_sim_species:'
         integer :: i,n,ndump
         integer :: npf
         real, dimension(3,100) :: arg
         character(len=20) :: sn,s1
         integer :: indz, xppc, yppc
         real :: min, max, cwp, n0
         real :: qm, qbm, dz

         this%err => input%err
         this%p => input%pp
         this%sp3 => input%sp
         this%sp2 => input%sp

         call this%err%werrfl2(class//sname//' started')

         call input%get('simulation.n0',n0)
         call input%get('simulation.indz',indz)

         cwp=5.32150254*1e9/sqrt(n0)
         call input%get('simulation.box.z(1)',min)
         call input%get('simulation.box.z(2)',max)
         dz=(max-min)/real(2**indz)

         call input%get('simulation.nspecies',n)

         allocate(this%spe(n),this%pf(n))

         do i = 1, n

            write (sn,'(I3.3)') i
            s1 = 'species('//trim(sn)//')'
            call input%get(trim(s1)//'.profile',npf)
            select case (npf)
            case (0)
               allocate(fdist2d_000::this%pf(i)%p)
               call this%pf(i)%p%new(input,i)
            case (10)
               allocate(fdist2d_010::this%pf(i)%p)
               call this%pf(i)%p%new(input,i)
            case (11)
               allocate(fdist2d_011::this%pf(i)%p)
               call this%pf(i)%p%new(input,i)
            case (12)
               allocate(fdist2d_012::this%pf(i)%p)
               call this%pf(i)%p%new(input,i)
            case (200)
               allocate(fdist2d_200::this%pf(i)%p)
               call this%pf(i)%p%new(input,i)
! Add new distributions right above this line
            case default
               write (erstr,*) 'Invalid species profile number:', npf
               call this%err%equit(class//sname//erstr)
            end select

            call input%get(trim(s1)//'.q',qm)
            call input%get(trim(s1)//'.m',qbm)
            qbm = qm/qbm
            call this%spe(i)%new(this%p,this%err,this%sp3,this%pf(i)%p,&
            &qbm=qbm,dt=dz,ci=1.0,xdim=13,s=s)
            ! if (npf == 200) then
            !    call this%spe(i)%new(this%p,this%err,this%sp3,this%pf(i)%p,&
            !    &qbm=qbm,dt=dz,ci=1.0,xdim=10,s=s)
            ! else
            !    call this%spe(i)%new(this%p,this%err,this%sp3,this%pf(i)%p,&
            !    &qbm=qbm,dt=dz,ci=1.0,xdim=8,s=s)
            ! end if               

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
      subroutine go_simulation_1(this)

         implicit none

         class(simulation), intent(inout) :: this
! local data
         character(len=18), save :: sname = 'go_simulation_1:'
         integer, dimension(10) :: istat
         integer :: i,j,l,m,n,ierr

         call this%err%werrfl2(class//sname//' started')

         ! dump initialization data
         if (this%start3d == 1) then
            this%tstep = 0
            call this%diag_simulation()
         endif

         do i = this%start3d, this%nstep3d

            this%tstep = i
            write (erstr,*) '3D step:', i
            call this%err%werrfl0(erstr)

            do m = 1, this%nbeams
               this%tag_bq(m,1) = ntag()
               this%tag_bq(m,2) = ntag()
               call this%beams%beam(m)%qdp(this%id_bq(m,1),this%id_bq(m,2),&
               &this%id_bq(m,3),this%tag_bq(m,1),this%tag_bq(m,2))
            end do

            do l =  1, this%nspecies
               this%tag_spe(l) = ntag()
               call this%species%spe(l)%precv(this%tag_spe(l))
            end do
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

            do j = 1, this%nstep2d
               write (erstr,*) '2D step:', j
               call this%err%werrfl0(erstr)
               this%tstep2 = j+this%noffstep2d
               call this%diag2_simulation()
               if (j == this%nstep2d) then
                  do m = 1, this%nbeams
                     call MPI_WAIT(this%id_bq(m,2),istat,ierr)
                  end do
               endif
               call this%fields%qb%as(0.0)
               do m = 1, this%nbeams
                  call this%beams%beam(m)%qdp(this%fields%qb,j+1)
               end do
               call this%fields%qb%fftrk(1)
               call this%fields%qb%elf(this%fields%epwb)
               call this%fields%qe%as(0.0)
               do l = 1, this%nspecies
                  call this%species%spe(l)%qdp(this%fields%qe)
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
                     call this%species%spe(m)%amjdp(this%fields%fxy,this%fields%bxyz,this%fields%psit,&
                     &this%fields%cu,this%fields%amu,this%fields%dcu,this%dex)
                  end do
                  call this%fields%cu%mult(this%fields%cu,this%dex)
                  call this%fields%amu%mult(this%fields%amu,this%dex)
                  call this%fields%dcu%mult(this%fields%dcu,this%dex)
                  if (l == this%iter) then
                     do m = 1, this%nspecies
                        call this%species%spe(m)%cbq(j+1)
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
               if (allocated(this%fields%vpot)) then
                  call this%fields%cu%vpot(this%fields%vpot)
                  call this%fields%vpot%fftkr(1)
                  call this%fields%vpot%mult(this%fields%vpot,-this%dex)
                  call this%fields%vpot%cb(this%fields%vpot3d,j+1,(/1,2,3/),(/1,2,3/))
               end if               
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

            call this%diag_simulation()

            do m = 1, this%nspecies
               call MPI_WAIT(this%id_spe(m),istat,ierr)
               call this%species%spe(m)%renew(i*this%dt)
            end do
         end do

         call this%err%werrfl2(class//sname//' ended')

      end subroutine go_simulation_1
!
      subroutine go_simulation_2(this)

         implicit none

         class(simulation), intent(inout) :: this
! local data
         character(len=18), save :: sname = 'go_simulation_2:'
         integer, dimension(10) :: istat
         integer :: i,j,l,m,n,ierr

         call this%err%werrfl2(class//sname//' started')

         ! dump initialization data
         if (this%start3d == 1) then
            this%tstep = 0
            call this%diag_simulation()
         endif

         do i = this%start3d, this%nstep3d

            this%tstep = i
            write (erstr,*) '3D step:', i
            call this%err%werrfl0(erstr)

            do m = 1, this%nbeams
               this%tag_bq(m,1) = ntag()
               this%tag_bq(m,2) = ntag()
               call this%beams%beam(m)%qdp(this%id_bq(m,1),this%id_bq(m,2),&
               &this%id_bq(m,3),this%tag_bq(m,1),this%tag_bq(m,2)) !3d charge density of beam
            end do

            do l =  1, this%nspecies
               this%tag_spe(l) = ntag()
               call this%species%spe(l)%precv(this%tag_spe(l))
            end do
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

            call this%fields%fxy%cb(this%fields%bexyz,1,(/1,2/),(/1,2/)) !copy 2d field to 3d field at the local longitudinal position lpos
            call this%fields%psit%cb(this%fields%bexyz,1,(/1/),(/3/))
            call this%fields%bxyz%cb(this%fields%bbxyz,1,(/1,2,3/),(/1,2,3/))

            do j = 1, this%nstep2d
               write (erstr,*) '2D step:', j
               call this%err%werrfl0(erstr)
               this%tstep2 = j+this%noffstep2d
               call this%diag2_simulation()               
               if (j == this%nstep2d) then
                  do m = 1, this%nbeams
                     call MPI_WAIT(this%id_bq(m,2),istat,ierr)
                  end do
               endif
               call this%fields%qb%as(0.0)
               do m = 1, this%nbeams
                  call this%beams%beam(m)%qdp(this%fields%qb,j+1) !2d charge density of beam,copy and add the charge density to a 2d slice
               end do
               write(2,*) this%fields%qb%getsum() , "qb"
!qb is the normalized unit in k space.               
               call this%fields%qb%fftrk(1)
!epwb is the normalized unit/dex in k space.                
               call this%fields%qb%elf(this%fields%epwb) !E⊥ of beam
               write(2,*) this%fields%epwb%getsum() , "epwb in k"              
               call this%fields%qe%as(0.0)
               do l = 1, this%nspecies
                  call this%species%spe(l)%qdp(this%fields%qe) !deposit rho-Jz
               end do
               call this%fields%qe1%as(0.0)
               call this%fields%qe1%as(this%fields%qe)
               write(2,*) this%fields%qe%getsum() , "qe"
!qe is the normalized unit in k space.                
               call this%fields%qe%fftrk(1)
!psi is the normalized unit/dex**2 in k space.             
               call this%fields%qe%pot(this%fields%psi)
               call this%fields%fxy%as(0.0)
!fxy is the normalized unit/dex in k space.                
               call this%fields%psi%grad(this%fields%fxy)
!fxy is the normalized unit/dex in r space.               
               call this%fields%fxy%fftkr(1)
               write(2,*) this%fields%fxy%getsum() , "fxy"
!psi is the normalized unit/dex**2 in r space.               
               call this%fields%psi%fftkr(1)
               write(2,*) this%fields%psi%getsum() , "psi"
               do l = 1, this%nspecies
                  call this%species%spe(l)%extpsi(this%fields%psi,this%dex)
               end do
!psi is the normalized unit in r space.               
               call this%fields%psi%mult(this%fields%psi,this%dex*this%dex)
               if (allocated(this%fields%psi3d)) call this%fields%psi%cb(this%fields%psi3d,j+1,(/1/),(/1/))
               call this%fields%cu%as(0.0)
!J, cu is the normalized unit/dex in r space.              
               do m = 1, this%nspecies
                  call this%species%spe(m)%jdp(this%fields%cu,this%dex) !deposit J
               end do
!cu is the normalized unit in r space.               
               call this%fields%cu%mult(this%fields%cu,this%dex)
               do m = 1, this%nspecies
                  call this%species%spe(m)%cbq(j+1)
               end do  
               if (allocated(this%fields%cu3d)) then
                  call this%fields%cu%cb(this%fields%cu3d,j+1,(/1,2,3/),(/1,2,3/))
               end if
               write(2,*) this%fields%cu%getsum() , "cu-r-N"
!cu is the normalized unit in k space.
               call this%fields%cu%fftrk(1)
               write(2,*) this%fields%cu%getsum() , "cu-k-N"
!potencil A,???               
               ! if (allocated(this%fields%vpot)) then
               !    call this%fields%cu%vpot(this%fields%vpot)
               !    call this%fields%vpot%fftkr(1)
               !    call this%fields%vpot%mult(this%fields%vpot,-this%dex)
               !    call this%fields%vpot%cb(this%fields%vpot3d,j+1,(/1,2,3/),(/1,2,3/))
               ! end if 

               call this%fields%cu%div(this%fields%div_vpot)
!Ez of plasma, psit is the normalized unit/dex in k space.              
               call this%fields%div_vpot%pot(this%fields%psit)
!Ez=Ez*(-dex)，psit is the normalized unit in k space.                
               call this%fields%psit%mult(this%fields%psit,-this%dex)
!psit is the normalized unit in r space.               
               call this%fields%psit%fftkr(1)
               !write(2,*) this%fields%psit%getsum() , "psit"
               call this%fields%bxyz%as(0.0)
!Bz of plasma, bxyz is the normalized unit/dex in k space.               
               call this%fields%cu%bf(this%fields%bxyz)
!Bz is the normalized unit in k space.             
               call this%fields%bxyz%mult(this%fields%bxyz,(/3/),(/3/),(/this%dex/))
!Bx=-Ey(beam)=Bx(beam),By=Ex(beam)=By(beam),Bx and By are the normalized unit/dex in k space.               
               call this%fields%bxyz%mult(this%fields%epwb,(/1,2/),(/2,1/),(/-1.0,1.0/))
!Bx and By are the normalized unit/dex and Bz is the normalized unit in r space.
               call this%fields%bxyz%fftkr(2)
               call this%fields%dcu%as(0.0)
               call this%fields%amu%as(0.0)              
!partial J⊥/partial n; psit and Bz(bxyz) are the normalized unit. bxy、fxy, dcu is normalized unit/dex, amu is the normalized unit/dex^2.               
               do m = 1, this%nspecies
                  call this%species%spe(m)%bdp(this%fields%fxy,this%fields%bxyz,this%fields%psit,&
                  &this%fields%dcu,this%fields%amu,this%dex) !source of equation when solving B⊥ in rspace.
               end do 
               write(2,*) this%fields%amu%getsum() , "amu-r-N/dex^2"
!amu=amu*dex，the amu is the normalized unit/dex in r space.
               call this%fields%amu%mult(this%fields%amu,this%dex) 
!dcu=dcu*dex，the dcu is the normalized unit in r space.            
               call this%fields%dcu%mult(this%fields%dcu,this%dex)
               write(2,*) this%fields%dcu%getsum() , "dcu-r-N" 
!amu is the normalized unit/dex in k space.           
               call this%fields%amu%fftrk(3) 
               call this%fields%damu%as(0.0)
!cu(k) is the normalized unit,amu(k) is the normalized unit/dex, damu is the normalized unit in k space.
               !call this%fields%cu%bfqpn(this%fields%dcu,this%fields%amu,this%fields%damu,this%dex)
               call this%fields%cu%bfqpn(this%fields%amu,this%fields%damu,this%dex)
!damu is the normalized unit in r space.               
               call this%fields%damu%fftkr(2)
               call this%fields%damu%sub(this%fields%damu,this%fields%dcu,(/1/),(/1/),(/2/))
               call this%fields%damu%add(this%fields%damu,this%fields%dcu,(/2/),(/2/),(/1/))
               write(2,*) this%fields%damu%getsum() , "damu-r-N"
               !call this%fields%damu%mult(this%fields%damu,-1.0)
               call this%fields%bxyzfd%as(0.0)   
!qe1 and psi are the normalized unit, damu is the normalized unit in r space.                            
               call this%fields%damu%renew_fd(this%fields%qe1,this%fields%psi,this%dex) !update stenci l when qe and psi change.
!damu and bxyzfd are the normalized unit. The bxyzfd is the B⊥. fd does not change units.   
               call this%fields%damu%solve_mg(this%fields%bxyzfd,(/2,3/),3) !damu and bxyzfd are both in rspace.
!bxyzfd is the normalized uni/dex in r space.
               call this%fields%bxyzfd%mult(this%fields%bxyzfd,1.0/this%dex)
               write(2,*) this%fields%bxyzfd%getsum(), "bxyzfd(B_perp)-r-N/dex"
               call this%fields%bxyzfd%sum_ele((/1/),this%bxp)
               write(2,*) this%bxp , "bxp"
               call this%fields%bxyzfd%sum_ele((/2/),this%byp)
               write(2,*) this%byp , "byp"
!Bx= Bx(beam)+Bx(plasma), Bx is normalized unit/dex in rspace.     
               call this%fields%bxyz%add(this%fields%bxyz,this%fields%bxyzfd,(/1/),(/1/),(/1/))
!By= By(beam)+By(plasma), By is the normalized unit/dex in rspace.              
               call this%fields%bxyz%add(this%fields%bxyz,this%fields%bxyzfd,(/2/),(/2/),(/2/))         
!fxy is the normalized unit/dex in r space.                 
               call this%fields%fxy%mult(this%fields%fxy,-1.0)
               call this%fields%fxy%add(this%fields%fxy,this%fields%bxyz,(/1/),(/1/),(/2/)) !Ex
               call this%fields%fxy%sub(this%fields%fxy,this%fields%bxyz,(/2/),(/2/),(/1/)) !Ey
!ex、ey(fxy)、bx and by are the normalized unit/dex,bz、ez(psit) are the normalized unit.              
               do m = 1, this%nspecies
                  call this%species%spe(m)%push(this%fields%fxy,this%fields%bxyz,this%fields%psit,&
                  &this%dex) !push particle
               end do
!fxy is the normalized unit.
               call this%fields%fxy%mult(this%fields%fxy,this%dex)               
               call this%fields%fxy%cb(this%fields%bexyz,j+1,(/1,2/),(/1,2/))
!psit is the normalized unit.               
               call this%fields%psit%cb(this%fields%bexyz,j+1,(/1/),(/3/))
!Bx=Bx*dex, By=By*dex, bxyz is the normalized unit.
               call this%fields%bxyz%mult(this%fields%bxyz,(/1,2/),(/1,2/),(/this%dex,this%dex/))
!bxyz is the normalized unit.               
               call this%fields%bxyz%cb(this%fields%bbxyz,j+1,(/1,2,3/),(/1,2,3/))
               write(2,*) this%fields%fxy%getsum() , "fxy_3d"
               call this%fields%fxy%sum_ele((/1/),this%Ex)
               write(2,*) this%Ex , "Ex"
               call this%fields%fxy%sum_ele((/2/),this%Ey)
               write(2,*) this%Ey , "Ey"
               write(2,*) this%fields%psit%getsum() , "psit_3d"
               write(2,*) this%fields%bxyz%getsum() , "B_3d"
               call this%fields%bxyz%sum_ele((/1/),this%Bx)
               write(2,*) this%Bx , "Bx"
               call this%fields%bxyz%sum_ele((/2/),this%By)
               write(2,*) this%By , "By"
               call this%fields%bxyz%sum_ele((/3/),this%Bz)
               write(2,*) this%Bz , "Bz"
            enddo

            do m = 1, this%nspecies
               call this%species%spe(m)%psend(this%tag_spe(m),this%id_spe(m))
            end do
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

            call this%diag_simulation()

            do m = 1, this%nspecies
               call MPI_WAIT(this%id_spe(m),istat,ierr)
               call this%species%spe(m)%renew(i*this%dt)
            end do
         end do

         call this%err%werrfl2(class//sname//' ended')

      end subroutine go_simulation_2
!
      subroutine init_diag(this)

         implicit none

         class(simulation), intent(inout), target :: this
! local data
         character(len=18), save :: sname = 'init_diag:'
         integer :: n_diag = 0, ndump, slice, slice_pos, sample
         integer :: n, m, l, i, j, k, ii, n_diag2 = 0
         character(len=20) :: s1, s2, s3, s4, sn1, sn2, sn3, sn4
         character(len=:), allocatable :: ss,sl
         real :: min, max, n0, dt
         real :: alx1, aly1, alz1, alx2, aly2, alz2
         logical :: rst
         integer :: ierr, indx, indy, indz, dim, dim2
         class (*), pointer :: obj => null(), obj2 => null()


         call this%err%werrfl2(class//sname//' started')

         call this%in%get('simulation.n0',n0)
         call this%in%get('simulation.indx',indx)
         call this%in%get('simulation.indy',indy)
         call this%in%get('simulation.indz',indz)

         call this%in%get('simulation.box.x(1)',min)
         call this%in%get('simulation.box.x(2)',max)
         alx1 = min
         alx2 = max
         call this%in%get('simulation.box.y(1)',min)
         call this%in%get('simulation.box.y(2)',max)
         aly1 = min
         aly2 = max
         call this%in%get('simulation.box.z(1)',min)
         call this%in%get('simulation.box.z(2)',max)
         alz1 = min
         alz2 = max
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
                     if (allocated(ss)) deallocate(ss)
                     call this%in%get('species('//trim(s1)//').diag'//'('//trim(s2)//').name'&
                     &//'(1)',ss)
                     if (ss == 'raw') then
                        n_diag2 = n_diag2 + 1
                     else
                        n_diag = n_diag + l
                     end if
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
               else if(this%in%found('field.diag('//trim(s1)//').debug')) then
                  n_diag2 = n_diag2 + l
               else
                  n_diag = n_diag + l
               end if
            end if
         end do
         call this%in%get('simulation.dump_restart',rst)
         if (rst) n_diag = n_diag + this%nbeams

         allocate(this%diag(n_diag), this%diag2(n_diag2))
         n_diag = 0
         n_diag2 = 0
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
                              this%diag(n_diag)%obj => this%beams%beam(i)
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
                                 call system('mkdir -p ./Beam'//trim(s1)//'/Charge_slice_'//trim(s4)//'/')
                              end if
                              call this%diag(n_diag)%file%new(&
                              &timeunits = '1 / \omega_p',&
                              &dt = dt,&
                              &axisunits = (/'c / \omega_p','c / \omega_p','c / \omega_p'/),&
                              &rank = 2,&
                              &filename = './Beam'//trim(s1)//'/Charge_slice_'//trim(s4)//'/',&
                              &dataname = 'charge_slice_'//sl,&
                              &units = 'n_0',&
                              &label = 'Charge Density')
                           end do
                        else
                           n_diag = n_diag + 1
                           this%diag(n_diag)%df = ndump
                           this%diag(n_diag)%obj => this%beams%beam(i)
                           allocate(this%diag(n_diag)%dim)
                           this%diag(n_diag)%dim = 1
                           if (this%p%getidproc() == 0) then
                              call system('mkdir -p ./Beam'//trim(s1)//'/Charge/')
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
                           &filename = './Beam'//trim(s1)//'/Charge/',&
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
                           call system('mkdir -p ./Beam'//trim(s1)//'/Raw/')
                        end if
                        call this%diag(n_diag)%file%new(&
                        &timeunits = '1 / \omega_p',&
                        &dt = dt,&
                        &ty = 'particles',&
                        &filename = './Beam'//trim(s1)//'/Raw/',&
                        &dataname = 'raw',&
                        &units = ' ',&
                        &label = 'Beam Raw')
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
                        sn1 = 'Charge'
                        sn2 = 'charge'
                        sn3 = 'n_0'
                        sn4 = 'Charge Density'
                        dim = 1
                        obj => this%species%spe(i)
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
                     case ('raw')
                         obj => this%species%spe(i)                       
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
                              call system('mkdir -p ./Species'//trim(s1)//'/'//trim(sn1)//'_slice_'//trim(s4)//'/')
                           end if
                           call this%diag(n_diag)%file%new(&
                           &timeunits = '1 / \omega_p',&
                           &dt = dt,&
                           &axisunits = (/'c / \omega_p','c / \omega_p','c / \omega_p'/),&
                           &rank = 2,&
                           &filename = './Species'//trim(s1)//'/'//trim(sn1)//'_slice_'//trim(s4)//'/',&
                           &dataname = trim(sn2)//'_slice_'//sl,&
                           &units = trim(sn3),&
                           &label = trim(sn4))
                        end do
                     else
                        if (ss == 'raw') then
                           n_diag2 = n_diag2 + 1
                           this%diag2(n_diag2)%df = ndump
                           this%diag2(n_diag2)%obj => obj
                           this%diag2(n_diag2)%path2 = './Species'//trim(s1)//'/Raw/'
                           ! allocate(this%diag(n_diag)%psample)
                           ! call this%in%get('beam('//trim(s1)//').diag'//'('//trim(s2)//').&
                           ! &sample',this%diag(n_diag)%psample)
                           if (this%p%getidproc() == 0) then
                              call system('mkdir -p ./Species'//trim(s1)//'/Raw/')
                           end if
                           call this%diag2(n_diag2)%file%new(&
                           &timeunits = '1 / \omega_p',&
                           &dt = dt,&
                           &ty = 'particles',&
                           &filename = './Species'//trim(s1)//'/Raw/',&
                           &dataname = 'raw',&
                           &units = ' ',&
                           &label = 'Test Particles Raw Data')
                        else
                           n_diag = n_diag + 1
                           this%diag(n_diag)%df = ndump
                           this%diag(n_diag)%obj => obj
                           allocate(this%diag(n_diag)%dim)
                           this%diag(n_diag)%dim = dim
                           if (this%p%getidproc() == 0) then
                              call system('mkdir -p ./Species'//trim(s1)//'/'//trim(sn1)//'/')
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
                           &filename = './Species'//trim(s1)//'/'//trim(sn1)//'/',&
                           &dataname = trim(sn2),&
                           &units = trim(sn3),&
                           &label = trim(sn4))
                        end if
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
                     obj2 => this%fields%fxy
                     dim2 = 1
                  case ('ey')
                     sn1 = 'Ey'
                     sn2 = 'ey'
                     sn3 = 'mc\omega_p/e'
                     sn4 = 'Electric Field'
                     dim = 2
                     obj => this%fields%bexyz
                     obj2 => this%fields%fxy
                     dim2 = 2
                  case ('ez')
                     sn1 = 'Ez'
                     sn2 = 'ez'
                     sn3 = 'mc\omega_p/e'
                     sn4 = 'Electric Field'
                     dim = 3
                     obj => this%fields%bexyz
                     obj2 => this%fields%psit
                     dim2 = 1
                  case ('bx')
                     sn1 = 'Bx'
                     sn2 = 'bx'
                     sn3 = 'mc\omega_p/e'
                     sn4 = 'Magnetic Field'
                     dim = 1
                     obj => this%fields%bbxyz
                     obj2 => this%fields%bxyz
                     dim2 = 1
                  case ('by')
                     sn1 = 'By'
                     sn2 = 'by'
                     sn3 = 'mc\omega_p/e'
                     sn4 = 'Magnetic Field'
                     dim = 2
                     obj => this%fields%bbxyz
                     obj2 => this%fields%bxyz
                     dim2 = 2
                  case ('bz')
                     sn1 = 'Bz'
                     sn2 = 'bz'
                     sn3 = 'mc\omega_p/e'
                     sn4 = 'Magnetic Field'
                     dim = 3
                     obj => this%fields%bbxyz
                     obj2 => this%fields%bxyz
                     dim2 = 3
                  case ('psi')
                     sn1 = 'Psi'
                     sn2 = 'psi'
                     sn3 = 'mc^2'
                     sn4 = '\Psi'
                     dim = 1
                     obj => this%fields%psi3d
                  case ('vpotx')
                     sn1 = 'Ax'
                     sn2 = 'ax'
                     sn3 = 'mc^2/e'
                     sn4 = 'Vector Potential'
                     dim = 1
                     obj => this%fields%vpot3d
                  case ('vpoty')
                     sn1 = 'Ay'
                     sn2 = 'ay'
                     sn3 = 'mc^2/e'
                     sn4 = 'Vector Potential'
                     dim = 2
                     obj => this%fields%vpot3d
                  case ('vpotz')
                     sn1 = 'Az'
                     sn2 = 'az'
                     sn3 = 'mc^2/e'
                     sn4 = 'Vector Potential'
                     dim = 3
                     obj => this%fields%vpot3d
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
                  else if (this%in%found('field.diag('//trim(s1)//').debug')) then
                        n_diag2 = n_diag2 + 1
                        this%diag2(n_diag2)%df = ndump
                        this%diag2(n_diag2)%obj => obj2
                        ! allocate(this%diag2(n_diag2)%slice,this%diag(n_diag)%slice_pos)
                        allocate(this%diag2(n_diag2)%dim)
                        this%diag2(n_diag2)%dim = dim2
                        ! if (allocated(sl)) deallocate(sl)
                        ! write (s3, '(I4.4)') k
                        ! call this%in%get('field.diag('//trim(s1)//').&
                        ! &slice('//trim(s3)//').(1)',sl)
                        call this%diag2(n_diag2)%file%new(&
                        &axismin = (/alx1,aly1,alz1/),&
                        &axismax = (/alx2,aly2,alz2/),&
                        &axisname  = (/'x','y','z'/),&
                        &axislabel = (/'x','y','z'/))

                        ! select case (sl)
                        ! case ('yz')
                        !    this%diag(n_diag)%slice = 1
                        !    call this%diag(n_diag)%file%new(&
                        !    &axismin = (/aly1,alz1,alx1/),&
                        !    &axismax = (/aly2,alz2,alx2/),&
                        !    &axisname  = (/'y  ','\xi','x  '/),&
                        !    &axislabel = (/'y  ','\xi','x  '/))
                        ! case ('xz')
                        !    this%diag(n_diag)%slice = 2
                        !    call this%diag(n_diag)%file%new(&
                        !    &axismin = (/alx1,alz1,aly1/),&
                        !    &axismax = (/alx2,alz2,aly2/),&
                        !    &axisname  = (/'x  ','\xi','y  '/),&
                        !    &axislabel = (/'x  ','\xi','y  '/))
                        ! case ('xy')
                        !    this%diag(n_diag)%slice = 3
                        !    call this%diag(n_diag)%file%new(&
                        !    &axismin = (/alx1,aly1,alz1/),&
                        !    &axismax = (/alx2,aly2,alz2/),&
                        !    &axisname  = (/'x','y','z'/),&
                        !    &axislabel = (/'x','y','z'/))
                        ! end select
                        ! call this%in%get('field.diag('//trim(s1)//').&
                        ! &slice('//trim(s3)//').(2)',this%diag(n_diag)%slice_pos)
                        if (this%p%getidproc() == 0) then
                           call system('mkdir -p ./Fields/'//trim(sn1)//'_debug'//'/')
                        end if
                        this%diag2(n_diag2)%path2 = './Fields/'//trim(sn1)//'_debug'//'/'

                        call this%diag2(n_diag2)%file%new(&
                        &timeunits = '1 / \omega_p',&
                        &dt = dt,&
                        &axisunits = (/'c / \omega_p','c / \omega_p','c / \omega_p'/),&
                        &rank = 2,&
                        &filename = './Fields/'//trim(sn1)//'_debug'//'/',&
                        &dataname = trim(sn2),&
                        &units = trim(sn3),&
                        &label = trim(sn4))
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
               if (this%p%getidproc() == 0) then
                  call system('mkdir -p ./RST/Beam-'//trim(s1)//'/')
               end if
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
      subroutine diag2_simulation(this)

         implicit none

         class(simulation), intent(inout) :: this
! local data
         character(len=18), save :: sname = 'diag2_simulation:'
         character(len=20) :: s1
         integer :: n, m, l, i, j, k, ierr, idn = 9
         integer, dimension(10) :: istat
         real :: dt

         call this%err%werrfl2(class//sname//' started')

         call this%in%get('simulation.dt',dt)

         n = size(this%diag2)

         do i = 1, n
            if (mod(this%tstep,this%diag2(i)%df) == 0) then
               write (s1, '(I8.8)') this%tstep
               if (this%p%getidproc() == 0) then
                  call system('mkdir -p '//trim(this%diag2(i)%path2)//trim(s1)//'/')
               end if               
               select type (obj => this%diag2(i)%obj)
               type is (field2d)
                  call this%diag2(i)%file%new(n = this%tstep2, t = real(this%tstep2),filename = trim(this%diag2(i)%path2)&
                  &//trim(s1)//'/')
                  call obj%wr(this%diag2(i)%file,this%diag2(i)%dim)
               type is (species2d)
                  call this%diag2(i)%file%new(n = this%tstep2, t = real(this%tstep2),filename = trim(this%diag2(i)%path2)&
                  &//trim(s1)//'/')
                  call obj%wr(this%diag2(i)%file,(/this%dex,this%dex/))
               end select
            end if
         end do

         call this%err%werrfl2(class//sname//' ended')

      end subroutine diag2_simulation
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
            if (mod(this%tstep,this%diag(i)%df) == 0) then
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
                  if (allocated(this%diag(i)%psample)) then
                     this%tag(1) = ntag()
                     call MPI_WAIT(this%id(idn+i),istat,ierr)
                     call obj%wr(this%diag(i)%file,this%diag(i)%psample,(/this%dex,this%dex,this%dxi/),&
                     &this%tag(1),this%tag(1),this%id(idn+i))
                  else if (allocated(this%diag(i)%slice)) then
                     this%tag(1) = ntag()
                     call MPI_WAIT(this%id(idn+i),istat,ierr)
                     call obj%wrq(this%diag(i)%file,this%diag(i)%slice,this%diag(i)%slice_pos,&
                     &this%tag(1),this%tag(1),this%id(idn+i))
                  else if (allocated(this%diag(i)%dim)) then
                     this%tag(1) = ntag()
                     call MPI_WAIT(this%id(idn+i),istat,ierr)
                     call obj%wrq(this%diag(i)%file,this%tag(1),this%tag(1),this%id(idn+i))
                  else
                     call obj%wrst(this%diag(i)%file)
                  end if
               type is (species2d)
                  if (allocated(this%diag(i)%slice)) then
                     this%tag(1) = ntag()
                     call MPI_WAIT(this%id(idn+i),istat,ierr)
                     call obj%wrq(this%diag(i)%file,this%diag(i)%slice,this%diag(i)%slice_pos,&
                     &this%tag(1),this%tag(1),this%id(idn+i))
                  else
                     this%tag(1) = ntag()
                     call MPI_WAIT(this%id(idn+i),istat,ierr)
                     call obj%wrq(this%diag(i)%file,this%tag(1),this%tag(1),this%id(idn+i))
                  end if
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