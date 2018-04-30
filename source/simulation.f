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
      
      type simulation

         private

!         type(input) :: sim
         type(input_json) :: sim_json
         class(parallel_pipe),pointer :: p => null()
         class(perrors),pointer :: err => null()
         class(spect3d), pointer :: sp3 => null()
         class(spect2d), pointer :: sp2 => null()
         
         type(sim_fields) :: fields
         type(sim_beams) :: beams
         type(sim_species) :: species

         contains
         
         generic :: new => init_simulation
         generic :: del => end_simulation
         
         procedure, private :: init_simulation, end_simulation

      end type simulation
!
      type sim_fields

         private

         class(parallel_pipe),pointer :: p => null()
         class(perrors),pointer :: err => null()
         class(spect3d), pointer :: sp3 => null()
         class(spect2d), pointer :: sp2 => null()
         type(field2d), allocatable :: qb, qe, qi, psit, psi, div_vpot, reg
         type(field2d), allocatable :: fxy, bxyz, cu, dcu, amu, epw, epwb         
         type(field3d), allocatable :: bexyz, bbxyz
         type(field3d), allocatable :: qeb
         type(field3d), allocatable :: psi3d

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
         integer :: nbeams

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
      character(len=10) :: class = 'simulation:'
      character(len=128) :: erstr
      	
      type(hdf5file) :: file2d,file3d,filep,file_rst

      integer :: ierr, iter, nstep3d, nstep2d, start3d
      integer, dimension(54) :: id
      integer, dimension(8) :: tag
      integer, dimension(:), allocatable :: tag_spe, id_spe, id_br
      integer, dimension(:), allocatable :: tag_beam, id_beam, id_qep, id_qeps
      integer, dimension(10) :: istat
      integer :: i,j,l,m,n
      character(len=128) :: erstr
      character(len=20) :: stime
      
      real :: dex, dxi, dex2
            
      contains
!
      subroutine init_simulation(this)

         implicit none
         
         class(simulation), intent(inout) :: this
! local data
         character(len=18), save :: sname = 'init_simulation:'
         
         
         call this%sim_json%new()
         this%err => this%sim_json%err
         this%p => this%sim_json%pp
         this%sp3 => this%sim_json%sp
         this%sp2 => this%sim_json%sp

         call this%err%werrfl2(class//sname//' started')
         
         call this%fields%new(this%sim_json)
         call this%beams%new(this%sim_json)
         call this%species%new(this%sim_json)

         call this%err%werrfl2(class//sname//' ended')

      end subroutine init_simulation
!         
      subroutine end_simulation(this)

         implicit none

         class(simulation), intent(inout) :: this
! local data
         character(len=18), save :: sname = 'end_simulation:'
         
         call this%err%werrfl2(class//sname//' started')
         
         call this%fields%del()
         call this%beams%del()
         call this%species%del()

         call this%err%werrfl2(class//sname//' ended')

      end subroutine end_simulation
!
      subroutine init_sim_fields(this,input)

         implicit none

         class(sim_fields), intent(inout) :: this
         type(input_json), intent(in) :: input
! local data
         character(len=18), save :: sname = 'init_sim_fields:'
         character(len=20) :: sn
         character(len=:), allocatable :: ff         
         integer:: i,n,ndump

         this%err => this%input%err
         this%p => this%input%pp
         this%sp3 => this%input%sp
         this%sp2 => this%input%sp

         call this%err%werrfl2(class//sname//' started')

         allocate(this%qb, this%qe, this%qi, this%psit, this%psi)
         allocate(this%div_vpot, this%reg, this%fxy, this%bxyz, this%cu)
         allocate(this%dcu, this%amu, this%epw, this%epwb)
         allocate(this%bexyz, this%bbxyz, this%qeb)
         
         call input%info('field.diag',n_children=n)
         
         do i = 1, n
            call str(i,sn,4)
            call input%get('field.diag('//trim(sn)//').ndump',ndump)
            if (ndump > 0) then
               call input%get('field.diag('//trim(sn)//').name',ff)
               if (ff == 'psi') then
                  allocate(this%psi3d)
                  call this%psi3d%new(this%p,this%err,this%sp3,dim=1)
                  exit
               end if
            end if
         end do
         
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

         call this%err%werrfl2(class//sname//' ended')

      end subroutine init_sim_fields
!
      subroutine end_sim_fields(this)

         implicit none

         class(sim_fields), intent(inout) :: this
         type(input_json), intent(in) :: input
! local data
         character(len=18), save :: sname = 'end_sim_fields:'
         
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

         call this%err%werrfl2(class//sname//' ended')


      end subroutine end_sim_fields
!
      subroutine init_sim_beams(this)

         implicit none

         class(sim_fields), intent(inout) :: this
         type(input_json), intent(in) :: input
! local data
         character(len=18), save :: sname = 'end_sim_fields:'
         
         call this%err%werrfl2(class//sname//' started')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         
         
         
         allocate(pf(sim%sim%nspecies),spe(sim%sim%nspecies))
         allocate(pf3(sim%sim%nbeams),beam(sim%sim%nbeams))
         allocate(qep(sim%sim%nspecies),qe0(sim%sim%nspecies))
         allocate(cu0(sim%sim%nspecies),dcu0(sim%sim%nspecies),amu0(sim%sim%nspecies))
         
         do i = 1, sim%sim%nspecies
            call qep(i)%new(pp,perr,psp3,dim=1)
            call qe0(i)%new(pp,perr,psp2,dim=1,fftflag=.false.)
            call cu0(i)%new(pp,perr,psp2,dim=3,fftflag=.false.)
            call dcu0(i)%new(pp,perr,psp2,dim=2,fftflag=.false.)
            call amu0(i)%new(pp,perr,psp2,dim=3,fftflag=.false.)
            call pf(i)%new(pp,perr,psp2,npf = sim%species(i)%prof,npx=sim%species(i)%npx,&
            &npy=sim%species(i)%npy,arg=sim%species(i)%arg(1:2,:))
            call spe(i)%new(pp,perr,psp2,pf(i),qe0(i),qm=sim%species(i)%q,&
            &qbm=sim%species(i)%q/sim%species(i)%m,dt=sim%sim%dz,ci=1.0,&
            &xdim=7,npmax=sim%species(i)%npmax,nbmax=int(0.01*sim%species(i)%npmax))
         end do
         
         do i = 1, sim%sim%nbeams
            call pf3(i)%new(pp,perr,psp3,npf = sim%beam(i)%prof,npx=sim%beam(i)%npx,&
            &npy=sim%beam(i)%npy,npz=sim%beam(i)%npz,arg=sim%beam(i)%arg)
            call beam(i)%new(pp,perr,psp3,pf3(i),qeb,qm=sim%beam(i)%q,qbm=sim%beam(i)%q/&
            &sim%beam(i)%m,dt=sim%sim%dt,ci=1.0,xdim=6,npmax=sim%beam(i)%npmax,&
            &nbmax=int(sim%beam(i)%npmax*0.01))
         end do
                     
         call file3d%new(timeunits = '1 / \omega_p',dt = sim%sim%dt,&
         &axisname  = (/'x  ','y  ','\xi'/), axislabel = (/'x  ','y  ','\xi'/),&
         &axisunits = (/'c / \omega_p','c / \omega_p','c / \omega_p'/),&
         &axismax = (/sim%sim%lx,sim%sim%ly,sim%sim%lz/),rank = 3)

         call file2d%new(timeunits = '1 / \omega_p',dt = sim%sim%dt,&
         &axisname  = (/'x  ','y  ','\xi'/), axislabel = (/'x  ','y  ','\xi'/),&
         &axisunits = (/'c / \omega_p','c / \omega_p','c / \omega_p'/),&
         &axismax = (/sim%sim%lx,sim%sim%ly,sim%sim%lz/),rank = 2)

         call filep%new(timeunits = '1 / \omega_p',dt = sim%sim%dt,ty = 'particles')

         call file_rst%new(ty='rst',rank = 2)
         
         nstep2d = 2**psp3%getindz()/pp%getnstage()
         nstep3d = sim%sim%time / sim%sim%dt
         start3d = 1
         dex = sim%sim%dx
         dxi = sim%sim%dz
         dex2 = dex * dex
         iter = sim%sim%iter

         
         if (pp%getidproc() == 0) then
            if (sim%diag%dfpsi > 0) call system('mkdir ./PSI')
            if (sim%diag%dfqep > 0) then
               do i = 1, sim%sim%nspecies
                  call str(i,stime,2)                  
                  call system('mkdir ./QEP'//trim(stime))
               end do
            end if
            if (sim%diag%dfqeb > 0) then
               call system('mkdir ./QEB')
            end if
            if (sim%diag%dfjp > 0) then
               call system('mkdir ./JPX')
               call system('mkdir ./JPY')
               call system('mkdir ./JPZ')
            end if
            if (sim%diag%dfe > 0) then
               call system('mkdir ./FEX')
               call system('mkdir ./FEY')
               call system('mkdir ./FEZ')
            end if
            if (sim%diag%dfb > 0) then
               call system('mkdir ./FBX')
               call system('mkdir ./FBY')
               call system('mkdir ./FBZ')
            end if
            if (sim%diag%dfpsislice > 0) then
               call system('mkdir ./PSI-XZ')
               call system('mkdir ./PSI-YZ')
               call system('mkdir ./PSI-XY')
            end if
            if (sim%diag%dfqebslice > 0) then
               call system('mkdir ./QEB-XZ')
               call system('mkdir ./QEB-YZ')
               call system('mkdir ./QEB-XY')
            end if
            if (sim%diag%dfqepslice > 0) then
               do i = 1, sim%sim%nspecies
                  call str(i,stime,2)                  
                  call system('mkdir ./QEP'//trim(stime)//'-XZ')
                  call system('mkdir ./QEP'//trim(stime)//'-YZ')
                  call system('mkdir ./QEP'//trim(stime)//'-XY')
               end do
            end if
            if (sim%diag%dfjpslice > 0) then
               call system('mkdir ./JPX-XZ')
               call system('mkdir ./JPX-YZ')
               call system('mkdir ./JPX-XY')
               call system('mkdir ./JPY-XZ')
               call system('mkdir ./JPY-YZ')
               call system('mkdir ./JPY-XY')
               call system('mkdir ./JPZ-XZ')
               call system('mkdir ./JPZ-YZ')
               call system('mkdir ./JPZ-XY')
            end if
            if (sim%diag%dfeslice > 0) then
               call system('mkdir ./FEX-XZ')
               call system('mkdir ./FEX-YZ')
               call system('mkdir ./FEX-XY')
               call system('mkdir ./FEY-XZ')
               call system('mkdir ./FEY-YZ')
               call system('mkdir ./FEY-XY')
               call system('mkdir ./FEZ-XZ')
               call system('mkdir ./FEZ-YZ')
               call system('mkdir ./FEZ-XY')
            end if
            if (sim%diag%dfbslice > 0) then
               call system('mkdir ./FBX-XZ')
               call system('mkdir ./FBX-YZ')
               call system('mkdir ./FBX-XY')
               call system('mkdir ./FBY-XZ')
               call system('mkdir ./FBY-YZ')
               call system('mkdir ./FBY-XY')
               call system('mkdir ./FBZ-XZ')
               call system('mkdir ./FBZ-YZ')
               call system('mkdir ./FBZ-XY')
            end if
            if (sim%diag%dump_beam_raw) then
               call system('mkdir ./RAW-BEAM')
               do i = 1, sim%sim%nbeams
                  call str(i,stime,2)                  
                  call system('mkdir ./RAW-BEAM/'//trim(stime))
               end do
            end if
            if (sim%res%read_rst_file) then
               do m = 1, sim%sim%nbeams
                  call str(m,sn,2)                  
                  call str(pp%getidproc(),sid,8)
                  call str(sim%res%rst_timestep,stime,8)                  
                  call file_rst%new(filename = './RST/RST-'//trim(sn)//'-'//trim(sid)//&
                  &'_'//trim(stime)//'.h5',dataname = 'RST-'//trim(sn)//'-'//trim(sid)//&
                  &'_'//trim(stime)//'.h5')
                  call beam(m)%rrst(file_rst)
                end do
                start3d = sim%res%rst_timestep
            else if (sim%res%dump_rst_file) then
               call system('mkdir ./RST')
            end if
            
         endif
         
         call MPI_BARRIER(pp%getlworld(),ierr)

         allocate(tag_spe(sim%sim%nspecies),id_spe(sim%sim%nspecies))
         allocate(id_qep(sim%sim%nspecies),id_qeps(sim%sim%nspecies))
         allocate(tag_beam(sim%sim%nbeams),id_beam(sim%sim%nbeams),id_br(sim%sim%nbeams))
         id(:) = MPI_REQUEST_NULL
         id_spe(:) = MPI_REQUEST_NULL
         id_qep(:) = MPI_REQUEST_NULL
         id_qeps(:) = MPI_REQUEST_NULL
         id_beam(:) = MPI_REQUEST_NULL                 
         id_br(:) = MPI_REQUEST_NULL                 

      end subroutine initialization     
!
      subroutine diagnostic()
         implicit none
         
         character(len=20) :: sn, sid
         
         call str(i,stime,8)                  
         
         if ((sim%diag%dfqeb > 0) .and. (mod((i-1),sim%diag%dfqeb) == 0)) then
            call file3d%new(filename = './QEB/QEB_'//trim(stime)//'.h5',&
            &dataname = 'QEB',units = 'n_0',label = 'QEB',n = i, t = i*sim%sim%dt)
            tag(1) = ntag()
            call MPI_WAIT(id(11),istat,ierr)
            call qeb%wr(file3d,1,tag(1),tag(1),id(11))
         end if

         if ((sim%diag%dfqep > 0) .and. (mod((i-1),sim%diag%dfqep) == 0)) then
            do m = 1, sim%sim%nspecies
               call str(m,sn,2)                  
               call file3d%new(filename = './QEP'//trim(sn)//'/QEP'//trim(sn)//&
               &'_'//trim(stime)//'.h5',&
               &dataname = 'QEP'//trim(sn),units = 'n_0',label = 'QEP'//trim(sn),&
               &n = i, t = i*sim%sim%dt)
               tag(1) = ntag()
               call MPI_WAIT(id_qep(m),istat,ierr)
               call qep(m)%wr(file3d,1,tag(1),tag(1),id_qep(m))
            end do
         end if

         if ((sim%diag%dfpsi > 0) .and. (mod((i-1),sim%diag%dfpsi) == 0)) then
            call file3d%new(filename = './PSI/PSI_'//trim(stime)//'.h5',&
            &dataname = 'PSI',units = 'mc^2',label = 'PSI',n = i, t = i*sim%sim%dt)
            tag(1) = ntag()
            call MPI_WAIT(id(12),istat,ierr)
            call psi3d%wr(file3d,1,tag(1),tag(1),id(12))
         end if

         if ((sim%diag%dfjp > 0) .and. (mod((i-1),sim%diag%dfjp) == 0)) then
            call file3d%new(filename = './JPX/JPX_'//trim(stime)//'.h5',&
            &dataname = 'JPX',units = 'n_0c',label = 'JPX',n = i, t = i*sim%sim%dt)
            tag(1) = ntag()
            call MPI_WAIT(id(13),istat,ierr)
            call cu3d%wr(file3d,1,tag(1),tag(1),id(13))

            call file3d%new(filename = './JPY/JPY_'//trim(stime)//'.h5',&
            &dataname = 'JPY',units = 'n_0c',label = 'JPY',n = i, t = i*sim%sim%dt)
            tag(1) = ntag()
            call MPI_WAIT(id(14),istat,ierr)
            call cu3d%wr(file3d,2,tag(1),tag(1),id(14))
            
            call file3d%new(filename = './JPZ/JPZ_'//trim(stime)//'.h5',&
            &dataname = 'JPZ',units = 'n_0c',label = 'JPZ',n = i, t = i*sim%sim%dt)
            tag(1) = ntag()
            call MPI_WAIT(id(15),istat,ierr)
            call cu3d%wr(file3d,3,tag(1),tag(1),id(15))            
         end if

         if ((sim%diag%dfe > 0) .and. (mod((i-1),sim%diag%dfe) == 0)) then
            call file3d%new(filename = './FEX/FEX_'//trim(stime)//'.h5',&
            &dataname = 'FEX',units = 'mc\omega_p/e',label = 'FEX',n = i,t = i*sim%sim%dt)
            tag(1) = ntag()
            call MPI_WAIT(id(16),istat,ierr)
            call bexyz%wr(file3d,1,tag(1),tag(1),id(16))

            call file3d%new(filename = './FEY/FEY_'//trim(stime)//'.h5',&
            &dataname = 'FEY',units = 'mc\omega_p/e',label = 'FEY',n = i,t = i*sim%sim%dt)
            tag(1) = ntag()
            call MPI_WAIT(id(17),istat,ierr)
            call bexyz%wr(file3d,2,tag(1),tag(1),id(17))
            
            call file3d%new(filename = './FEZ/FEZ_'//trim(stime)//'.h5',&
            &dataname = 'FEZ',units = 'mc\omega_p/e',label = 'FEZ',n = i,t = i*sim%sim%dt)
            tag(1) = ntag()
            call MPI_WAIT(id(18),istat,ierr)
            call bexyz%wr(file3d,3,tag(1),tag(1),id(18))            
         end if         

         if ((sim%diag%dfb > 0) .and. (mod((i-1),sim%diag%dfb) == 0)) then
            call file3d%new(filename = './FBX/FBX_'//trim(stime)//'.h5',&
            &dataname = 'FBX',units = 'mc\omega_p/e',label = 'FBX',n = i,t = i*sim%sim%dt)
            tag(1) = ntag()
            call MPI_WAIT(id(19),istat,ierr)
            call bbxyz%wr(file3d,1,tag(1),tag(1),id(19))

            call file3d%new(filename = './FBY/FBY_'//trim(stime)//'.h5',&
            &dataname = 'FBY',units = 'mc\omega_p/e',label = 'FBY',n = i,t = i*sim%sim%dt)
            tag(1) = ntag()
            call MPI_WAIT(id(20),istat,ierr)
            call bbxyz%wr(file3d,2,tag(1),tag(1),id(20))
            
            call file3d%new(filename = './FBZ/FBZ_'//trim(stime)//'.h5',&
            &dataname = 'FBZ',units = 'mc\omega_p/e',label = 'FBZ',n = i,t = i*sim%sim%dt)
            tag(1) = ntag()
            call MPI_WAIT(id(21),istat,ierr)
            call bbxyz%wr(file3d,3,tag(1),tag(1),id(21))            
         end if         

         if ((sim%diag%dfqebslice > 0) .and. (mod((i-1),sim%diag%dfqebslice) == 0)) then
            call file2d%new(filename = './QEB-XZ/QEB-XZ_'//trim(stime)//'.h5',&
            &dataname = 'QEB-XZ',units = 'n_0',label = 'QEB-XZ',n = i, t = i*sim%sim%dt,&
            &rank = 2, axisname = (/'x  ','\xi','z  '/),&
            &axislabel = (/'x  ','\xi','z  '/),&
            &axismax = (/sim%sim%lx,sim%sim%lz,1.0/), axismin = (/0.0,0.0,0.0/))
            tag(1) = ntag()
            call MPI_WAIT(id(22),istat,ierr)
            call qeb%wr(file2d,1,2,sim%diag%qeby0,tag(1),tag(1),id(22))

            call file2d%new(filename = './QEB-YZ/QEB-YZ_'//trim(stime)//'.h5',&
            &dataname = 'QEB-YZ',units = 'n_0',label = 'QEB-YZ',n = i, t = i*sim%sim%dt,&
            &rank = 2, axisname  = (/'y  ','\xi','x3 '/),&
            &axislabel = (/'y  ','\xi','x3 '/),&
            &axismax = (/sim%sim%ly,sim%sim%lz,1.0/), axismin = (/0.0,0.0,0.0/))
            tag(1) = ntag()
            call MPI_WAIT(id(23),istat,ierr)
            call qeb%wr(file2d,1,1,sim%diag%qebx0,tag(1),tag(1),id(23))

            call file2d%new(filename = './QEB-XY/QEB-XY_'//trim(stime)//'.h5',&
            &dataname = 'QEB-XY',units = 'n_0',label = 'QEB-XY',n = i, t = i*sim%sim%dt,&
            &rank = 2, axisname  = (/'x','y','z'/),axislabel = (/'x','y','z'/),&
            &axismax = (/sim%sim%lx,sim%sim%ly,1.0/), axismin = (/0.0,0.0,0.0/))
            tag(1) = ntag()
            call MPI_WAIT(id(24),istat,ierr)
            call qeb%wr(file2d,1,3,sim%diag%qebz0,tag(1),tag(1),id(24))

         end if

         if ((sim%diag%dfqepslice > 0) .and. (mod((i-1),sim%diag%dfqepslice) == 0)) then
            do m = 1, sim%sim%nspecies
               call str(m,sn,2)                  
               call file2d%new(filename = './QEP'//trim(sn)//'-XZ/QEP'//trim(sn)//&
               &'-XZ_'//trim(stime)//'.h5',&
               &dataname = 'QEP'//trim(sn)//'-XZ',units = 'n_0',&
               &label = 'QEP'//trim(sn)//'-XZ',&
               &n = i, t = i*sim%sim%dt,rank = 2, axisname  = (/'x  ','\xi','x3 '/),&
               &axislabel = (/'x  ','\xi','x3 '/),axismax=(/sim%sim%lx,sim%sim%lz,1.0/),&
               &axismin = (/0.0,0.0,0.0/))
               tag(1) = ntag()
               call MPI_WAIT(id_qeps(m),istat,ierr)
               call qep(m)%wr(file2d,1,2,sim%diag%qepy0,tag(1),tag(1),id_qeps(m))

               call file2d%new(filename = './QEP'//trim(sn)//'-YZ/QEP'//trim(sn)//&
               &'-YZ_'//trim(stime)//'.h5',&
               &dataname = 'QEP'//trim(sn)//'-YZ',units = 'n_0',&
               &label = 'QEP'//trim(sn)//'-YZ',&
               &n = i, t = i*sim%sim%dt,rank = 2, axisname  = (/'y  ','\xi','x3 '/),&
               &axislabel = (/'y  ','\xi','x3 '/),axismax=(/sim%sim%ly,sim%sim%lz,1.0/),&
               &axismin = (/0.0,0.0,0.0/))
               tag(1) = ntag()
               call MPI_WAIT(id_qeps(m),istat,ierr)
               call qep(m)%wr(file2d,1,1,sim%diag%qepx0,tag(1),tag(1),id_qeps(m))

               call file2d%new(filename = './QEP'//trim(sn)//'-XY/QEP'//trim(sn)//&
               &'-XY_'//trim(stime)//'.h5',&
               &dataname = 'QEP'//trim(sn)//'-XY',units = 'n_0',&
               &label = 'QEP'//trim(sn)//'-XY',&
               &n = i, t = i*sim%sim%dt,rank = 2, axisname  = (/'x','y','z'/),&
               &axislabel = (/'x','y','z'/),axismax = (/sim%sim%lx,sim%sim%ly,1.0/),&
               &axismin = (/0.0,0.0,0.0/))
               tag(1) = ntag()
               call MPI_WAIT(id_qeps(m),istat,ierr)
               call qep(m)%wr(file2d,1,3,sim%diag%qepz0,tag(1),tag(1),id_qeps(m))
            end do
         end if

         if ((sim%diag%dfpsislice > 0) .and. (mod((i-1),sim%diag%dfpsislice) == 0)) then
            call file2d%new(filename = './PSI-XZ/PSI-XZ_'//trim(stime)//'.h5',&
            &dataname = 'PSI-XZ',units = 'mc^2',label = 'PSI-XZ',n = i,t = i*sim%sim%dt,&
            &rank = 2, axisname  = (/'x  ','\xi','x3 '/),&
            &axislabel = (/'x  ','\xi','x3 '/),&
            &axismax = (/sim%sim%lx,sim%sim%lz,1.0/), axismin = (/0.0,0.0,0.0/))
            tag(1) = ntag()
            call MPI_WAIT(id(25),istat,ierr)
            call psi3d%wr(file2d,1,2,sim%diag%psiy0,tag(1),tag(1),id(25))

            call file2d%new(filename = './PSI-YZ/PSI-YZ_'//trim(stime)//'.h5',&
            &dataname = 'PSI-YZ',units = 'mc^2',label = 'PSI-YZ',n = i,t = i*sim%sim%dt,&
            &rank = 2, axisname  = (/'y  ','\xi','x3 '/),&
            &axislabel = (/'y  ','\xi','x3 '/),&
            &axismax = (/sim%sim%ly,sim%sim%lz,1.0/), axismin = (/0.0,0.0,0.0/))
            tag(1) = ntag()
            call MPI_WAIT(id(26),istat,ierr)
            call psi3d%wr(file2d,1,1,sim%diag%psix0,tag(1),tag(1),id(26))

            call file2d%new(filename = './PSI-XY/PSI-XY_'//trim(stime)//'.h5',&
            &dataname = 'PSI-XY',units = 'mc^2',label = 'PSI-XY',n = i,t = i*sim%sim%dt,&
            &rank = 2, axisname  = (/'x','y','z'/),axislabel = (/'x','y','z'/),&
            &axismax = (/sim%sim%lx,sim%sim%ly,1.0/), axismin = (/0.0,0.0,0.0/))
            tag(1) = ntag()
            call MPI_WAIT(id(27),istat,ierr)
            call psi3d%wr(file2d,1,3,sim%diag%psiz0,tag(1),tag(1),id(27))
         end if

         if ((sim%diag%dfjpslice > 0) .and. (mod((i-1),sim%diag%dfjpslice) == 0)) then
            call file2d%new(filename = './JPX-XZ/JPX-XZ_'//trim(stime)//'.h5',&
            &dataname = 'JPX-XZ',units = 'n_0c',label = 'JPX-XZ',n = i,t = i*sim%sim%dt,&
            &rank = 2, axisname  = (/'x  ','\xi','x3 '/),&
            &axislabel = (/'x  ','\xi','x3 '/),&
            &axismax = (/sim%sim%lx,sim%sim%lz,1.0/), axismin = (/0.0,0.0,0.0/))
            tag(1) = ntag()
            call MPI_WAIT(id(28),istat,ierr)
            call cu3d%wr(file2d,1,2,sim%diag%jpy0,tag(1),tag(1),id(28))

            call file2d%new(filename = './JPX-YZ/JPX-YZ_'//trim(stime)//'.h5',&
            &dataname = 'JPX-YZ',units = 'n_0c',label = 'JPX-YZ',n = i,t = i*sim%sim%dt,&
            &rank = 2, axisname  = (/'y  ','\xi','x3 '/),&
            &axislabel = (/'y  ','\xi','x3 '/),&
            &axismax = (/sim%sim%ly,sim%sim%lz,1.0/), axismin = (/0.0,0.0,0.0/))
            tag(1) = ntag()
            call MPI_WAIT(id(29),istat,ierr)
            call cu3d%wr(file2d,1,1,sim%diag%jpx0,tag(1),tag(1),id(29))

            call file2d%new(filename = './JPX-XY/JPX-XY_'//trim(stime)//'.h5',&
            &dataname = 'JPX-XY',units = 'n_0c',label = 'JPX-XY',n = i,t = i*sim%sim%dt,&
            &rank = 2, axisname  = (/'x','y','z'/),axislabel = (/'x','y','z'/),&
            &axismax = (/sim%sim%lx,sim%sim%ly,1.0/), axismin = (/0.0,0.0,0.0/))
            tag(1) = ntag()
            call MPI_WAIT(id(30),istat,ierr)
            call cu3d%wr(file2d,1,3,sim%diag%jpz0,tag(1),tag(1),id(30))

            call file2d%new(filename = './JPY-XZ/JPY-XZ_'//trim(stime)//'.h5',&
            &dataname = 'JPY-XZ',units = 'n_0c',label = 'JPY-XZ',n = i,t = i*sim%sim%dt,&
            &rank = 2, axisname  = (/'x  ','\xi','x3 '/),&
            &axislabel = (/'x  ','\xi','x3 '/),&
            &axismax = (/sim%sim%lx,sim%sim%lz,1.0/), axismin = (/0.0,0.0,0.0/))
            tag(1) = ntag()
            call MPI_WAIT(id(32),istat,ierr)
            call cu3d%wr(file2d,2,2,sim%diag%jpy0,tag(1),tag(1),id(31))

            call file2d%new(filename = './JPY-YZ/JPY-YZ_'//trim(stime)//'.h5',&
            &dataname = 'JPY-YZ',units = 'n_0c',label = 'JPY-YZ',n = i,t = i*sim%sim%dt,&
            &rank = 2, axisname  = (/'y  ','\xi','x3 '/),&
            &axislabel = (/'y  ','\xi','x3 '/),&
            &axismax = (/sim%sim%ly,sim%sim%lz,1.0/), axismin = (/0.0,0.0,0.0/))
            tag(1) = ntag()
            call MPI_WAIT(id(32),istat,ierr)
            call cu3d%wr(file2d,2,1,sim%diag%jpx0,tag(1),tag(1),id(32))

            call file2d%new(filename = './JPY-XY/JPY-XY_'//trim(stime)//'.h5',&
            &dataname = 'JPY-XY',units = 'n_0c',label = 'JPY-XY',n = i,t = i*sim%sim%dt,&
            &rank = 2, axisname  = (/'x','y','z'/),axislabel = (/'x','y','z'/),&
            &axismax = (/sim%sim%lx,sim%sim%ly,1.0/), axismin = (/0.0,0.0,0.0/))
            tag(1) = ntag()
            call MPI_WAIT(id(33),istat,ierr)
            call cu3d%wr(file2d,2,3,sim%diag%jpz0,tag(1),tag(1),id(33))

            call file2d%new(filename = './JPZ-XZ/JPZ-XZ_'//trim(stime)//'.h5',&
            &dataname = 'JPZ-XZ',units = 'n_0c',label = 'JPZ-XZ',n = i,t = i*sim%sim%dt,&
            &rank = 2, axisname  = (/'x  ','\xi','x3 '/),&
            &axislabel = (/'x  ','\xi','x3 '/),&
            &axismax = (/sim%sim%lx,sim%sim%lz,1.0/), axismin = (/0.0,0.0,0.0/))
            tag(1) = ntag()
            call MPI_WAIT(id(34),istat,ierr)
            call cu3d%wr(file2d,3,2,sim%diag%jpy0,tag(1),tag(1),id(34))

            call file2d%new(filename = './JPZ-YZ/JPZ-YZ_'//trim(stime)//'.h5',&
            &dataname = 'JPZ-YZ',units = 'n_0c',label = 'JPZ-YZ',n = i,t = i*sim%sim%dt,&
            &rank = 2, axisname  = (/'y  ','\xi','x3 '/),&
            &axislabel = (/'y  ','\xi','x3 '/),&
            &axismax = (/sim%sim%ly,sim%sim%lz,1.0/), axismin = (/0.0,0.0,0.0/))
            tag(1) = ntag()
            call MPI_WAIT(id(35),istat,ierr)
            call cu3d%wr(file2d,3,1,sim%diag%jpx0,tag(1),tag(1),id(35))

            call file2d%new(filename = './JPZ-XY/JPZ-XY_'//trim(stime)//'.h5',&
            &dataname = 'JPZ-XY',units = 'n_0c',label = 'JPZ-XY',n = i,t = i*sim%sim%dt,&
            &rank = 2, axisname  = (/'x','y','z'/),axislabel = (/'x','y','z'/),&
            &axismax = (/sim%sim%lx,sim%sim%ly,1.0/), axismin = (/0.0,0.0,0.0/))
            tag(1) = ntag()
            call MPI_WAIT(id(36),istat,ierr)
            call cu3d%wr(file2d,3,3,sim%diag%jpz0,tag(1),tag(1),id(36))

         end if

         if ((sim%diag%dfeslice > 0) .and. (mod((i-1),sim%diag%dfeslice) == 0)) then
            call file2d%new(filename = './FEX-XZ/FEX-XZ_'//trim(stime)//'.h5',&
            &dataname = 'FEX-XZ',units = 'mc\omega_p/e',label = 'FEX-XZ',n = i,&
            &t = i*sim%sim%dt,rank = 2, axisname  = (/'x  ','\xi','x3 '/),&
            &axislabel = (/'x  ','\xi','x3 '/),&
            &axismax = (/sim%sim%lx,sim%sim%lz,1.0/), axismin = (/0.0,0.0,0.0/))
            tag(1) = ntag()
            call MPI_WAIT(id(37),istat,ierr)
            call bexyz%wr(file2d,1,2,sim%diag%ey0,tag(1),tag(1),id(37))

            call file2d%new(filename = './FEX-YZ/FEX-YZ_'//trim(stime)//'.h5',&
            &dataname = 'FEX-YZ',units = 'mc\omega_p/e',label = 'FEX-YZ',n = i,&
            &t = i*sim%sim%dt,rank = 2, axisname  = (/'y  ','\xi','x3 '/),&
            &axislabel = (/'y  ','\xi','x3 '/),&
            &axismax = (/sim%sim%ly,sim%sim%lz,1.0/), axismin = (/0.0,0.0,0.0/))
            tag(1) = ntag()
            call MPI_WAIT(id(38),istat,ierr)
            call bexyz%wr(file2d,1,1,sim%diag%ex0,tag(1),tag(1),id(38))

            call file2d%new(filename = './FEX-XY/FEX-XY_'//trim(stime)//'.h5',&
            &dataname = 'FEX-XY',units = 'mc\omega_p/e',label = 'FEX-XY',n = i,&
            &t = i*sim%sim%dt,&
            &rank = 2, axisname  = (/'x','y','z'/),axislabel = (/'x','y','z'/),&
            &axismax = (/sim%sim%lx,sim%sim%ly,1.0/), axismin = (/0.0,0.0,0.0/))
            tag(1) = ntag()
            call MPI_WAIT(id(39),istat,ierr)
            call bexyz%wr(file2d,1,3,sim%diag%ez0,tag(1),tag(1),id(39))

            call file2d%new(filename = './FEY-XZ/FEY-XZ_'//trim(stime)//'.h5',&
            &dataname = 'FEY-XZ',units = 'mc\omega_p/e',label = 'FEY-XZ',n = i,&
            &t = i*sim%sim%dt,rank = 2, axisname  = (/'x  ','\xi','x3 '/),&
            &axislabel = (/'x  ','\xi','x3 '/),&
            &axismax = (/sim%sim%lx,sim%sim%lz,1.0/), axismin = (/0.0,0.0,0.0/))
            tag(1) = ntag()
            call MPI_WAIT(id(40),istat,ierr)
            call bexyz%wr(file2d,2,2,sim%diag%ey0,tag(1),tag(1),id(40))

            call file2d%new(filename = './FEY-YZ/FEY-YZ_'//trim(stime)//'.h5',&
            &dataname = 'FEY-YZ',units = 'mc\omega_p/e',label = 'FEY-YZ',n = i,&
            &t = i*sim%sim%dt,rank = 2, axisname  = (/'y  ','\xi','x3 '/),&
            &axislabel = (/'y  ','\xi','x3 '/),&
            &axismax = (/sim%sim%ly,sim%sim%lz,1.0/), axismin = (/0.0,0.0,0.0/))
            tag(1) = ntag()
            call MPI_WAIT(id(41),istat,ierr)
            call bexyz%wr(file2d,2,1,sim%diag%ex0,tag(1),tag(1),id(41))

            call file2d%new(filename = './FEY-XY/FEY-XY_'//trim(stime)//'.h5',&
            &dataname = 'FEY-XY',units = 'mc\omega_p/e',label = 'FEY-XY',n = i,&
            &t = i*sim%sim%dt,&
            &rank = 2, axisname  = (/'x','y','z'/),axislabel = (/'x','y','z'/),&
            &axismax = (/sim%sim%lx,sim%sim%ly,1.0/), axismin = (/0.0,0.0,0.0/))
            tag(1) = ntag()
            call MPI_WAIT(id(42),istat,ierr)
            call bexyz%wr(file2d,2,3,sim%diag%ez0,tag(1),tag(1),id(42))

            call file2d%new(filename = './FEZ-XZ/FEZ-XZ_'//trim(stime)//'.h5',&
            &dataname = 'FEZ-XZ',units = 'mc\omega_p/e',label = 'FEZ-XZ',n = i,&
            &t = i*sim%sim%dt,rank = 2, axisname  = (/'x  ','\xi','x3 '/),&
            &axislabel = (/'x  ','\xi','x3 '/),&
            &axismax = (/sim%sim%lx,sim%sim%lz,1.0/), axismin = (/0.0,0.0,0.0/))
            tag(1) = ntag()
            call MPI_WAIT(id(43),istat,ierr)
            call bexyz%wr(file2d,3,2,sim%diag%ey0,tag(1),tag(1),id(43))

            call file2d%new(filename = './FEZ-YZ/FEZ-YZ_'//trim(stime)//'.h5',&
            &dataname = 'FEZ-YZ',units = 'mc\omega_p/e',label = 'FEZ-YZ',n = i,&
            &t = i*sim%sim%dt,rank = 2, axisname  = (/'y  ','\xi','x3 '/),&
            &axislabel = (/'y  ','\xi','x3 '/),&
            &axismax = (/sim%sim%ly,sim%sim%lz,1.0/), axismin = (/0.0,0.0,0.0/))
            tag(1) = ntag()
            call MPI_WAIT(id(44),istat,ierr)
            call bexyz%wr(file2d,3,1,sim%diag%ex0,tag(1),tag(1),id(44))

            call file2d%new(filename = './FEZ-XY/FEZ-XY_'//trim(stime)//'.h5',&
            &dataname = 'FEZ-XY',units = 'mc\omega_p/e',label = 'FEZ-XY',n = i,&
            &t = i*sim%sim%dt,&
            &rank = 2, axisname  = (/'x','y','z'/),axislabel = (/'x','y','z'/),&
            &axismax = (/sim%sim%lx,sim%sim%ly,1.0/), axismin = (/0.0,0.0,0.0/))
            tag(1) = ntag()
            call MPI_WAIT(id(45),istat,ierr)
            call bexyz%wr(file2d,3,3,sim%diag%ez0,tag(1),tag(1),id(45))
         end if

         if ((sim%diag%dfbslice > 0) .and. (mod((i-1),sim%diag%dfbslice) == 0)) then
            call file2d%new(filename = './FBX-XZ/FBX-XZ_'//trim(stime)//'.h5',&
            &dataname = 'FBX-XZ',units = 'mc\omega_p/e',label = 'FBX-XZ',n = i,&
            &t = i*sim%sim%dt,rank = 2, axisname  = (/'x  ','\xi','x3 '/),&
            &axislabel = (/'x  ','\xi','x3 '/),&
            &axismax = (/sim%sim%lx,sim%sim%lz,1.0/), axismin = (/0.0,0.0,0.0/))
            tag(1) = ntag()
            call MPI_WAIT(id(46),istat,ierr)
            call bbxyz%wr(file2d,1,2,sim%diag%by0,tag(1),tag(1),id(46))

            call file2d%new(filename = './FBX-YZ/FBX-YZ_'//trim(stime)//'.h5',&
            &dataname = 'FBX-YZ',units = 'mc\omega_p/e',label = 'FBX-YZ',n = i,&
            &t = i*sim%sim%dt,rank = 2, axisname  = (/'y  ','\xi','x3 '/),&
            &axislabel = (/'y  ','\xi','x3 '/),&
            &axismax = (/sim%sim%ly,sim%sim%lz,1.0/), axismin = (/0.0,0.0,0.0/))
            tag(1) = ntag()
            call MPI_WAIT(id(47),istat,ierr)
            call bbxyz%wr(file2d,1,1,sim%diag%bx0,tag(1),tag(1),id(47))

            call file2d%new(filename = './FBX-XY/FBX-XY_'//trim(stime)//'.h5',&
            &dataname = 'FBX-XY',units = 'mc\omega_p/e',label = 'FBX-XY',n = i,&
            &t = i*sim%sim%dt,&
            &rank = 2, axisname  = (/'x','y','z'/),axislabel = (/'x','y','z'/),&
            &axismax = (/sim%sim%lx,sim%sim%ly,1.0/), axismin = (/0.0,0.0,0.0/))
            tag(1) = ntag()
            call MPI_WAIT(id(48),istat,ierr)
            call bbxyz%wr(file2d,1,3,sim%diag%bz0,tag(1),tag(1),id(48))

            call file2d%new(filename = './FBY-XZ/FBY-XZ_'//trim(stime)//'.h5',&
            &dataname = 'FBY-XZ',units = 'mc\omega_p/e',label = 'FBY-XZ',n = i,&
            &t = i*sim%sim%dt,rank = 2, axisname  = (/'x  ','\xi','x3 '/),&
            &axislabel = (/'x  ','\xi','x3 '/),&
            &axismax = (/sim%sim%lx,sim%sim%lz,1.0/), axismin = (/0.0,0.0,0.0/))
            tag(1) = ntag()
            call MPI_WAIT(id(49),istat,ierr)
            call bbxyz%wr(file2d,2,2,sim%diag%by0,tag(1),tag(1),id(49))

            call file2d%new(filename = './FBY-YZ/FBY-YZ_'//trim(stime)//'.h5',&
            &dataname = 'FBY-YZ',units = 'mc\omega_p/e',label = 'FBY-YZ',n = i,&
            &t = i*sim%sim%dt,rank = 2, axisname  = (/'y  ','\xi','x3 '/),&
            &axislabel = (/'y  ','\xi','x3 '/),&
            &axismax = (/sim%sim%ly,sim%sim%lz,1.0/), axismin = (/0.0,0.0,0.0/))
            tag(1) = ntag()
            call MPI_WAIT(id(50),istat,ierr)
            call bbxyz%wr(file2d,2,1,sim%diag%bx0,tag(1),tag(1),id(50))

            call file2d%new(filename = './FBY-XY/FBY-XY_'//trim(stime)//'.h5',&
            &dataname = 'FBY-XY',units = 'mc\omega_p/e',label = 'FBY-XY',n = i,&
            &t = i*sim%sim%dt,&
            &rank = 2, axisname  = (/'x','y','z'/),axislabel = (/'x','y','z'/),&
            &axismax = (/sim%sim%lx,sim%sim%ly,1.0/), axismin = (/0.0,0.0,0.0/))
            tag(1) = ntag()
            call MPI_WAIT(id(51),istat,ierr)
            call bbxyz%wr(file2d,2,3,sim%diag%bz0,tag(1),tag(1),id(51))

            call file2d%new(filename = './FBZ-XZ/FBZ-XZ_'//trim(stime)//'.h5',&
            &dataname = 'FBZ-XZ',units = 'mc\omega_p/e',label = 'FBZ-XZ',n = i,&
            &t = i*sim%sim%dt,rank = 2, axisname  = (/'x  ','\xi','x3 '/),&
            &axislabel = (/'x  ','\xi','x3 '/),&
            &axismax = (/sim%sim%lx,sim%sim%lz,1.0/), axismin = (/0.0,0.0,0.0/))
            tag(1) = ntag()
            call MPI_WAIT(id(52),istat,ierr)
            call bbxyz%wr(file2d,3,2,sim%diag%by0,tag(1),tag(1),id(52))

            call file2d%new(filename = './FBZ-YZ/FBZ-YZ_'//trim(stime)//'.h5',&
            &dataname = 'FBZ-YZ',units = 'mc\omega_p/e',label = 'FBZ-YZ',n = i,&
            &t = i*sim%sim%dt,&
            &rank = 2, axisname  = (/'y  ','\xi','x3 '/),&
            &axislabel = (/'y  ','\xi','x3 '/),&
            &axismax = (/sim%sim%ly,sim%sim%lz,1.0/), axismin = (/0.0,0.0,0.0/))
            tag(1) = ntag()
            call MPI_WAIT(id(53),istat,ierr)
            call bbxyz%wr(file2d,3,1,sim%diag%bx0,tag(1),tag(1),id(53))

            call file2d%new(filename = './FBZ-XY/FBZ-XY_'//trim(stime)//'.h5',&
            &dataname = 'FBZ-XY',units = 'mc\omega_p/e',label = 'FBZ-XY',n = i,&
            &t = i*sim%sim%dt,&
            &rank = 2, axisname  = (/'x','y','z'/),axislabel = (/'x','y','z'/),&
            &axismax = (/sim%sim%lx,sim%sim%ly,1.0/), axismin = (/0.0,0.0,0.0/))
            tag(1) = ntag()
            call MPI_WAIT(id(54),istat,ierr)
            call bbxyz%wr(file2d,3,3,sim%diag%bz0,tag(1),tag(1),id(54))
         end if

         if (sim%diag%dump_beam_raw .and. (mod((i-1),sim%diag%dfbeam_raw) == 0)) then
            do m = 1, sim%sim%nbeams
               call str(m,sn,2)
               call filep%new(filename = './RAW-BEAM/'//trim(sn)//'/RAW-BEAM-'//&
               &trim(sn)//'_'//trim(stime)//'.h5',&
               &dataname = 'RAW-BEAM-'//trim(sn),units = '',&
               &label = 'RAW-BEAM-'//trim(sn),n = i, t = i*sim%sim%dt,&
               &ty='particles')
               tag(1) = ntag()
               call MPI_WAIT(id_br(m),istat,ierr)
               call beam(m)%wr(filep,sim%diag%beam_raw_factor,(/dex,dex,dxi/),&
               &tag(1),tag(1),id_br(m))         
            end do
         end if

         if (sim%res%dump_rst_file .and. (mod((i-1),sim%res%dfrst) == 0)) then
            do m = 1, sim%sim%nbeams
               call str(m,sn,2)                  
               call str(pp%getidproc(),sid,8)
               call file_rst%new(filename = './RST/RST-'//trim(sn)//'-'//trim(sid)//&
               &'_'//trim(stime)//'.h5',dataname = 'RST-'//trim(sn)//'-'//trim(sid)//&
               &'_'//trim(stime))
               call beam(m)%wrst(file_rst)
             end do
         end if

      end subroutine diagnostic
!
      subroutine str(int_in,string,ndigits)

      implicit none
      integer, intent(in) :: int_in, ndigits
      character(len=20), intent(inout) :: string
      
! local variables      
      integer  ::  izero, i, nd, m
      character(len=20)  :: chindx 
      
      m = 1 
      izero =  ichar('0')
      if (ndigits > 20) then 
         nd = 20
      else
         nd = ndigits
      endif
      chindx = ''
      do i = nd, 1, -1 
         m = 10**(i-1)
         chindx = trim(chindx) // char(  izero + mod( int_in/m , 10 ) ) 
      enddo 
      string = trim(chindx)

      end subroutine str    
!                
      function ntag()
      
         implicit none
         integer, save :: tag = 0
         integer :: ntag
                 
         ntag = tag
         tag = tag + 1
      
      end function ntag
                  
      end module simulation_class