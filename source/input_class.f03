! input class for QuickPIC Open Source 1.0
! update: 04/18/2016

      module input_class

      use perrors_class
      use parallel_class
      use parallel_pipe_class
      use spect3d_class
      use mpi
      use json_module
         
      implicit none

      private

      public :: input, input_json
      
      type simsys_in
! The cell number of the simulation box is (nx,ny,nz) = (2**indx,2**indy,2**indz)
! time = total time length of the simulation
! dt = time interval for pushing the drive beam
! ax/ay/az = half-width of particle in x/y/z direction
! The length of the simulation box is (lx,ly,lz)

         integer :: indx = 8, indy = 8, indz = 8
         real :: time, dt
         real :: ax = .912871, ay = .912871, az = .912871
         real :: lx, ly, lz
         integer :: smonitor = 0 
         integer :: psolve = 1 , num_stages = 1
         real :: cwp
         real :: n0
         real :: dx, dy, dz
         integer :: iter
         integer :: nbeams 
         integer :: nspecies       
      end type
      
      type beam_in
! npx3/npy3/npz3 = initial number of particles distributed in x/y/z
! vtx/vty/vtz = thermal velocity of beam electrons in x/y/z direction
! vdx/vdy/vdz = drift velocity of beam electrons in x/y/z direction
! q = charge for 3D particle, in unit of e
! m = mass for 3D particle, in unit of electron mass 

         integer :: npmax
         integer :: npx, npy, npz
         real, dimension(3,100) :: arg
         real :: q, m
         integer :: prof
      end type
      
      type species_in
! npx, npy = number of particles in 2d code in x/y

         integer :: npx, npy
         integer :: npmax
         real, dimension(3,100) :: arg
         integer :: prof
         real :: q, m  
      end type
      
      type diag_in

         integer :: dfpsi,dfqep,dfqeb
         integer :: dfjp,dfe,dfb
         integer :: dfpsislice,dfqebslice,dfqepslice
         integer :: dfjpslice,dfjbslice
         integer :: dfeslice,dfbslice
         integer :: psix0,psiy0,psiz0
         integer :: qebx0,qeby0,qebz0
         integer :: qepx0,qepy0,qepz0
         integer :: jpx0,jpy0,jpz0
         integer :: ex0,ey0,ez0
         integer :: bx0,by0,bz0
         logical :: dump_beam_raw
         integer :: dfbeam_raw, beam_raw_factor
      end type
       
      type restart_in

         logical :: read_rst_file, dump_rst_file
         integer :: rst_timestep, dfrst
         
      end type

      type input

         class(spect3d), pointer :: sp => null()
         class(perrors), pointer :: err => null()
         class(parallel), pointer :: p => null()
         class(parallel_pipe), pointer :: pp => null()

         type (simsys_in) :: sim
         type (restart_in) :: res
         type (beam_in), dimension(:), allocatable :: beam
         type (species_in), dimension(:), allocatable :: species
         type (diag_in) :: diag

         contains
         
         generic :: new => read_input
         procedure, private :: read_input

      end type

      type input_json

         class(spect3d), pointer :: sp => null()
         class(perrors), pointer :: err => null()
         class(parallel), pointer :: p => null()
         class(parallel_pipe), pointer :: pp => null()
         type(json_file), private :: input

         contains
         
         generic :: new => read_input_json
         generic :: get => json_file_get_object,json_file_get_integer,&
         &json_file_get_double, json_file_get_logical,&
         &json_file_get_string, json_file_get_integer_vec,&
         &json_file_get_double_vec,json_file_get_logical_vec,&
         &json_file_get_string_vec,json_file_get_alloc_string_vec,&
         &json_file_get_root
         generic :: info => json_file_variable_info

         procedure, private :: read_input_json
         procedure, private :: initialize, set_json_core_in_file
         procedure, private :: load_file, print_to_string
         procedure, private :: load_from_string
         procedure, private :: json_file_get_object,json_file_get_integer,&
         &json_file_get_double, json_file_get_logical,&
         &json_file_get_string, json_file_get_integer_vec,&
         &json_file_get_double_vec,json_file_get_logical_vec,&
         &json_file_get_string_vec,json_file_get_alloc_string_vec,&
         &json_file_get_root
         procedure, private :: json_file_variable_info
         
      end type
      
      character(len=10), save :: class = 'input:'
      character(len=128), save :: erstr
      type(spect3d), save, target  :: sp
      type(perrors), save, target :: err
      type(parallel), save, target :: p
      type(parallel_pipe), save, target :: pp
     
      contains
!
      subroutine read_input(this)
    
         implicit none
         
         class(input), intent(inout) :: this
! local data
         character(len=18), save :: sname = 'read_input:'

         integer :: indx, indy, indz, npmax
         integer, dimension(:), allocatable :: min_bparticles
         integer, dimension(:), allocatable :: npx3, npy3, npz3
         integer, dimension(:), allocatable :: np2s
         integer :: inorder, popt, dopt
         integer :: djopt
         integer :: ntw = 1, ntp = 0
         real :: tend, dt
         real :: vtx, vty, vtz
         real, dimension(:), allocatable :: charge3s
         real, dimension(:), allocatable :: mass3s
         real, dimension(:), allocatable :: load_balance_ths,vt2xs,vt2ys
         real, dimension(:), allocatable :: vdx3, vdy3, vdz3
         real :: vtdx, vtdy, vtdz
         real :: ax = .912871, ay = .912871, az = .912871
         real, dimension(:), allocatable :: bgamma
         logical, dimension(:), allocatable :: ftwiss
         real, dimension(:), allocatable :: alphax, alphay, betax, betay
         real, dimension(:), allocatable :: emtx,emty 
         integer :: seunit = 2, smonitor = 0 
         integer :: emf = 0
         integer :: psolve = 1 
         integer :: ipbc = 1
         real :: alx, aly, alz
         integer :: dfpsi,dfphi,dfqep,dfqeb,dfbc
         integer :: dfchi,dfvp,dfjp,dfjb,dfe,dfb
         integer :: dfphislice,dfpsislice,dfqebslice,dfqepslice
         integer :: dfchislice,dfvpslice,dfjpslice,dfjbslice
         integer :: dfeslice,dfbslice
         real :: phix0,phiy0,phiz0
         real :: psix0,psiy0,psiz0
         real :: qebx0,qeby0,qebz0
         real :: qepx0,qepy0,qepz0
         real :: chix0,chiy0,chiz0
         real :: vpx0,vpy0,vpz0
         real :: jpx0,jpy0,jpz0
         real :: jbx0,jby0,jbz0
         real :: ex0,ey0,ez0
         real :: bx0,by0,bz0
         real :: cwp, gamma

         real :: plasma_density
         integer :: bc_diag_res
         integer :: version
         integer :: nbeams, num_stages
         integer :: nspecies=0, nneutrals=0, ntotal=0
         integer :: neutral_gas, neutral_z
         integer :: min_beam_particle,npx, npy, npz
         real :: vdx, vdy, vdz
         real :: dx, dy, dz
         real :: box_x, box_y, box_z
         real :: num_particle, load_balance_th
         real :: charge, mass
         integer :: np2, i,k
         real :: argx1,argx2,argx3, argx4, argx5
         integer :: nsrand
         real :: vt2x,vt2y
         integer :: profile_type
         real :: non_neutral_factor
         character(len=20) :: sboundary
         real, dimension(:),allocatable :: neutral_species,ionization_level
         integer, dimension(:), allocatable :: profile_types
         real, dimension(:), allocatable :: argx1s,argx2s,argx3s,argx4s,arg&
         &x5s
         integer, dimension(:), allocatable :: nsrands
         real, dimension(:), allocatable :: charge2s, mass2s
         logical, dimension(:), allocatable :: den_vars
         integer, dimension(:), allocatable :: den_var_nsecs
         integer, dimension(:), allocatable :: prof_nsecs
         real, dimension(:,:), allocatable :: den_var_fss, den_var_ss
         real, dimension(:,:,:), allocatable :: prof_paras
         logical, dimension(:), allocatable :: bevolutions
         logical, dimension(:), allocatable :: track
         integer, dimension(:), allocatable :: track_num
         logical, dimension(:), allocatable :: bquiets
         logical, dimension(:), allocatable :: buse_shifters, buse_destroye&
         &rs, buse_radiation_dampings
         integer, dimension(:), allocatable :: bshifter_nsecs,bdestroyer_nc&
         &riterias 
         logical :: beam_evolution, quiet_start
         logical :: use_shifter, use_destroyer, use_radiation_damping
         logical :: beam_track
         integer :: beam_track_num
         integer :: shifter_nsec, destroyer_ncriteria
         character(len=80), dimension(:), allocatable :: beam_profiles
         character(len=80) :: beam_profile
         integer :: temp
         real :: temp1
         integer, dimension(:,:), allocatable :: bp_sizes
         real, dimension(:,:,:), pointer :: ptr_temp

         real, dimension(:), allocatable :: non_neutral_factors
         logical :: dump_pha_beam
         integer :: dfpha_beam, dsample_beam      
         logical :: dump_pha_plasma
         integer :: dfpha_plasma, dsample_plasma
         logical :: laser_on 
         logical :: read_rst_file, dump_rst_file
         integer :: rst_timestep, dfrst, sortime_2d, sortime_3d
         real :: omx, omy, omz
         integer :: ntv, nts, sortime, nplot, idpal
         integer :: max_iter, verbose 
         real :: fac_exy, fac_bxy, fac_az, fac_bz, c_dif, j_dif 
         integer :: init_routine
         integer, dimension(:), allocatable :: iroutine
         real :: parameter_array(5,100), shifter_parameter(3,100)
         real :: destroyer_criteria(3,100)
         real, dimension(:,:,:), allocatable :: para_arr, bshifter_paras
         real, dimension(:,:,:), allocatable :: bdestroyer_criterias
         logical :: density_variation
         integer :: density_variation_nsec, prof_nsec
         real, dimension(100) :: density_variation_fs, density_variation_s 
         real, dimension(2,100) :: prof_parameter

         
         namelist /input_file/ version
         
         namelist /pipeline/ num_stages

         namelist /simulation_sys/ box_x, box_y, box_z, indx, indy, indz

         namelist /boundary/ sboundary
        
         namelist /num_beams/ nbeams 

         namelist /plasma/ nspecies,nneutrals,plasma_density
         
         namelist /beam/  beam_evolution, min_beam_particle, npx, npy, npz, &
         &charge, mass, gamma, num_particle, vdx, vdy, vdz, init_routine,   &
         &beam_profile,quiet_start, parameter_array, use_shifter,           &
         &shifter_nsec, shifter_parameter,use_destroyer,destroyer_ncriteria,&
         &destroyer_criteria,use_radiation_damping,beam_track,beam_track_num
        
         namelist /laser_input/ laser_on

         namelist /species/ load_balance_th, np2, charge, mass,             &
         &vt2x, vt2y, non_neutral_factor, profile_type,argx1,argx2,argx3,arg&
         &x4,argx5,nsr&
         &and,prof_nsec, prof_parameter, density_variation, density_variatio&
         &n_nsec,density_variation_fs,density_variation_s
         
         namelist /neutral/ neutral_gas, neutral_z
         
         namelist /simulation_time/ tend, dt
        
         namelist /potential_diag/  dfphi, dfphislice, phix0, phiy0, phiz0, &
         &dfpsi, dfpsislice, psix0, psiy0, psiz0

         namelist /ponderomotive_potential_diag/ dfvp,dfvpslice, vpx0, vpy0,&
         &vpz0

         namelist /chi_diag/ dfchi, dfchislice, chix0, chiy0, chiz0

         namelist /current_diag/ dfjp, dfjpslice, jpx0, jpy0, jpz0, dfjb, df&
         &jbslice, jbx0, jby0, jbz0

         namelist /field_diag/ dfe, dfeslice, ex0,ey0,ez0, dfb, dfbslice, bx&
         &0, by0, bz0
        
         namelist /beam_diag/  dfqeb, dfqebslice, qebx0, qeby0, qebz0,      &
         &dfbc, bc_diag_res     
        
         namelist /plasma_diag/ dfqep, dfqepslice, qepx0, qepy0, qepz0

         namelist /beam_phase_space_diag/ dump_pha_beam, dfpha_beam,        &
         &dsample_beam
         
         namelist /plasma_phase_space_diag/ dump_pha_plasma, dfpha_plasma,  &
         &dsample_plasma

         namelist /restart_file/ read_rst_file, rst_timestep,               &
         &dump_rst_file, dfrst
        
         namelist /optimization/ inorder, popt, dopt, djopt, sortime_2d,    &
         &sortime_3d
      
         namelist /debug/ max_iter, fac_exy, fac_bxy, fac_az, fac_bz, c_dif,&
         &j_dif, verbose      
         
         integer :: ierr
         real :: ntemp
         character(len=20) :: chari
         real, dimension(500) :: ddata
         real, dimension(3,100) :: barg
         real, dimension(2,100) :: sarg
         
         call p%new()
         
         this%p => p

         call err%new(this%p,2,monitor=0)
         
         this%err => err
         
         call this%err%werrfl0(class//sname//' started')
         
         if (p%getidproc() == 0) then
             open (unit=8,file='rpinput',form='formatted',status='old',iostat&
             &=ierr)
             if (ierr /= 0 ) then
                write (erstr,*) 'error: reading input file failed'
                call this%err%equit(class//sname//erstr)
                return
             end if
             read (8,pipeline)
             read (8,simulation_sys)
             read (8,boundary)
             read (8,num_beams)
             read (8,plasma)
             read (8,simulation_time)
             read (8,potential_diag)
             read (8,current_diag)
             read (8,field_diag)
             read (8,beam_diag)
             read (8,plasma_diag)
             read (8,beam_phase_space_diag)
             read (8,restart_file)
             read (8,debug)

             ntotal = nspecies + nneutrals
             cwp=5.32150254*1e9/sqrt(plasma_density)
             alx=box_x/cwp
             aly=box_y/cwp 
             alz=box_z/cwp 
             dx = alx/real(2**indx)
             dy = aly/real(2**indy)
             dz = alz/real(2**indz)
             
             psix0=psix0/cwp
             if ((psix0 < 0.0) .or. (psix0 >= alx)) then
                write (erstr,*) 'error: psix0 out of range'
                call this%err%equit(class//sname//erstr)
                return             
             end if
             psiy0=psiy0/cwp
             if ((psiy0 < 0.0) .or. (psiy0 >= aly)) then
                write (erstr,*) 'error: psiy0 out of range'
                call this%err%equit(class//sname//erstr)
                return             
             end if
             psiz0=psiz0/cwp
             if ((psiz0 < 0.0) .or. (psiz0 >= alz)) then
                write (erstr,*) 'error: psiz0 out of range'
                call this%err%equit(class//sname//erstr)
                return             
             end if
             qebx0=qebx0/cwp
             if ((qebx0 < 0.0) .or. (qebx0 >= alx)) then
                write (erstr,*) 'error: qebx0 out of range'
                call this%err%equit(class//sname//erstr)
                return             
             end if
             qeby0=qeby0/cwp
             if ((qeby0 < 0.0) .or. (qeby0 >= aly)) then
                write (erstr,*) 'error: qeby0 out of range'
                call this%err%equit(class//sname//erstr)
                return             
             end if
             qebz0=qebz0/cwp
             if ((qebz0 < 0.0) .or. (qebz0 >= alz)) then
                write (erstr,*) 'error: qebz0 out of range'
                call this%err%equit(class//sname//erstr)
                return             
             end if
             qepx0=qepx0/cwp
             if ((qepx0 < 0.0) .or. (qepx0 >= alx)) then
                write (erstr,*) 'error: qepx0 out of range'
                call this%err%equit(class//sname//erstr)
                return             
             end if
             qepy0=qepy0/cwp
             if ((qepy0 < 0.0) .or. (qepy0 >= aly)) then
                write (erstr,*) 'error: qepy0 out of range'
                call this%err%equit(class//sname//erstr)
                return             
             end if
             qepz0=qepz0/cwp
             if ((qepz0 < 0.0) .or. (qepz0 >= alz)) then
                write (erstr,*) 'error: qepz0 out of range'
                call this%err%equit(class//sname//erstr)
                return             
             end if
             jpx0=jpx0/cwp
             if ((jpx0 < 0.0) .or. (jpx0 >= alx)) then
                write (erstr,*) 'error: jpx0 out of range'
                call this%err%equit(class//sname//erstr)
                return             
             end if
             jpy0=jpy0/cwp
             if ((jpy0 < 0.0) .or. (jpy0 >= aly)) then
                write (erstr,*) 'error: jpy0 out of range'
                call this%err%equit(class//sname//erstr)
                return             
             end if
             jpz0=jpz0/cwp
             if ((jpz0 < 0.0) .or. (jpz0 >= alz)) then
                write (erstr,*) 'error: jpz0 out of range'
                call this%err%equit(class//sname//erstr)
                return             
             end if
             ex0=ex0/cwp
             if ((ex0 < 0.0) .or. (ex0 >= alx)) then
                write (erstr,*) 'error: ex0 out of range'
                call this%err%equit(class//sname//erstr)
                return             
             end if
             ey0=ey0/cwp
             if ((ey0 < 0.0) .or. (ey0 >= aly)) then
                write (erstr,*) 'error: ey0 out of range'
                call this%err%equit(class//sname//erstr)
                return             
             end if
             ez0=ez0/cwp
             if ((ez0 < 0.0) .or. (ez0 >= alz)) then
                write (erstr,*) 'error: ez0 out of range'
                call this%err%equit(class//sname//erstr)
                return             
             end if
             bx0=bx0/cwp
             if ((bx0 < 0.0) .or. (bx0 >= alx)) then
                write (erstr,*) 'error: bx0 out of range'
                call this%err%equit(class//sname//erstr)
                return             
             end if
             by0=by0/cwp
             if ((by0 < 0.0) .or. (by0 >= aly)) then
                write (erstr,*) 'error: by0 out of range'
                call this%err%equit(class//sname//erstr)
                return             
             end if
             bz0=bz0/cwp
             if ((bz0 < 0.0) .or. (bz0 >= alz)) then
                write (erstr,*) 'error: bz0 out of range'
                call this%err%equit(class//sname//erstr)
                return             
             end if

             phix0 =int(phix0/dx) + 1
             phiy0 =int(phiy0/dy) + 1
             phiz0 =int(phiz0/dz) + 1
             
             psix0 =int(psix0/dx) + 1
             psiy0 =int(psiy0/dy) + 1
             psiz0 =int(psiz0/dz) + 1
             
             qebx0 =int(qebx0/dx) + 1
             qeby0 =int(qeby0/dy) + 1
             qebz0 =int(qebz0/dz) + 1
             
             qepx0 =int(qepx0/dx) + 1
             qepy0 =int(qepy0/dy) + 1
             qepz0 =int(qepz0/dz) + 1

             vpx0 =int(vpx0/dx) + 1
             vpy0 =int(vpy0/dy) + 1
             vpz0 =int(vpz0/dz) + 1

             chix0 =int(chix0/dx) + 1
             chiy0 =int(chiy0/dy) + 1
             chiz0 =int(chiz0/dz) + 1

             jpx0 =int(jpx0/dx) + 1
             jpy0 =int(jpy0/dy) + 1
             jpz0 =int(jpz0/dz) + 1

             jbx0 =int(jbx0/dx) + 1
             jby0 =int(jby0/dy) + 1
             jbz0 =int(jbz0/dz) + 1
             
             ex0 =int(ex0/dx) + 1
             ey0 =int(ey0/dy) + 1
             ez0 =int(ez0/dz) + 1

             bx0 =int(bx0/dx) + 1
             by0 =int(by0/dy) + 1
             bz0 =int(bz0/dz) + 1

             sortime = sortime_3d  

             select case (sboundary)
             case ("conducting")
                ipbc = 1
                psolve = 1
             case default
                ipbc = 1
                psolve = 1
             end select

             ddata(1) = indx; ddata(2) = indy; ddata(3) = indz
             ddata(7) = inorder; ddata(8) = popt; ddata(9) = dopt
             ddata(10) = ntw; ddata(11) = ntp
             ddata(12) = tend; ddata(13) = dt
             ddata(14) = vtdx; ddata(15) = vtdy; ddata(16) = vtdz
             ddata(17) = ax; ddata(18) = ay; ddata(19) = az
             ddata(20) = alx; ddata(21) = aly; ddata(22) = alz; 
             ddata(23) = num_stages
             ddata(24) = dfpsi; ddata(25) = dfphi
             ddata(26) = dfqep; ddata(27) = dfqeb
             ddata(28) = dfphislice 
             ddata(29) = phix0; ddata(30) = phiy0; ddata(31) = phiz0
             ddata(32) = dfpsislice 
             ddata(33) = psix0; ddata(34) = psiy0; ddata(35) = psiz0
             ddata(36) = dfqebslice 
             ddata(37) = qebx0; ddata(38) = qeby0; ddata(39) = qebz0
             ddata(40) = dfqepslice 
             ddata(41) = qepx0; ddata(42) = qepy0; ddata(43) = qepz0
             ddata(44) = vt2x; ddata(45) = vt2y
             ddata(46) = bc_diag_res; ddata(47) = dfbc
             ddata(48) = dx; ddata(49) = dy; ddata(50) = dz
             ddata(51) = djopt; ddata(52) = sortime
             ddata(53) = non_neutral_factor; ddata(54) = max_iter
             ddata(55) = fac_exy; ddata(56) = fac_bxy
             ddata(57) = fac_az; ddata(58) = fac_bz
             ddata(59) = c_dif; ddata(60) = j_dif 
             ddata(61) = nbeams 
             if (read_rst_file) then
                ddata(62) = 1
             else 
                ddata(62) = 0
             endif
             ddata(63) = rst_timestep
             if (dump_rst_file) then
                ddata(64) = 1
             else
                ddata(64) = 0
             endif
             ddata(65) = dfrst   
             if (laser_on) then
                ddata(66) = 1
             else
                ddata(66) = 0
             endif
             ddata(67) = cwp
             ddata(68) = dfchi; ddata(69) = dfchislice 
             ddata(70) = chix0; ddata(71) = chiy0; ddata(72) = chiz0
             ddata(73) = dfvp; ddata(74) = dfvpslice 
             ddata(75) = vpx0; ddata(76) = vpy0; ddata(77) = vpz0
             ddata(78) = dfjp; ddata(79) = dfjpslice 
             ddata(80) = jpx0; ddata(81) = jpy0; ddata(82) = jpz0
             ddata(83) = dfjb; ddata(84) = dfjbslice 
             ddata(85) = jbx0; ddata(86) = jby0; ddata(87) = jbz0
             ddata(88) = dfe; ddata(89) = dfeslice 
             ddata(90) = ex0; ddata(91) = ey0; ddata(92) = ez0
             ddata(93) = dfb; ddata(94) = dfbslice 
             ddata(95) = bx0; ddata(96) = by0; ddata(97) = bz0
             ddata(98) = ipbc
             if (dump_pha_beam) then
                ddata(99) = 1
             else
                ddata(99) = 0
             endif
             ddata(100) = dfpha_beam
             ddata(101) = dsample_beam
             ddata(102) = verbose
             ddata(103) = ntotal
             ddata(104) = psolve
             rewind(8)
         end if
         
         call MPI_BCAST(ddata, 104, this%p%getmreal(), 0, this%p%getlworld(), ierr)

         indx = ddata(1); indy = ddata(2); indz = ddata(3)
         inorder = ddata(7); popt = ddata(8); dopt = ddata(9)
         ntw = ddata(10); ntp = ddata(11)
         tend = ddata(12); dt = ddata(13)
         vtdx = ddata(14); vtdy = ddata(15); vtdz = ddata(16)
         ax = ddata(17); ay = ddata(18); az = ddata(19)
         alx = ddata(20); aly = ddata(21); alz = ddata(22); 
         num_stages = ddata(23) 
         dfpsi = ddata(24); dfphi = ddata(25) 
         dfqep = ddata(26); dfqeb = ddata(27) 
         dfphislice = ddata(28)
         phix0 = ddata(29); phiy0 = ddata(30); phiz0 = ddata(31)  
         dfpsislice = ddata(32) 
         psix0 = ddata(33); psiy0 = ddata(34); psiz0 = ddata(35)
         dfqebslice = ddata(36) 
         qebx0 = ddata(37); qeby0 = ddata(38); qebz0 = ddata(39)
         dfqepslice = ddata(40) 
         qepx0 = ddata(41); qepy0 = ddata(42); qepz0 =ddata(43)  
         vt2x = ddata(44); vt2y = ddata(45)
         bc_diag_res = ddata(46); dfbc = ddata(47)
         dx = ddata(48); dy = ddata(49); dz = ddata(50)
         djopt = ddata(51); sortime = ddata(52)
         non_neutral_factor = ddata(53); max_iter = ddata(54)
         fac_exy = ddata(55); fac_bxy = ddata(56)
         fac_az = ddata(57); fac_bz = ddata(58)
         c_dif = ddata(59); j_dif = ddata(60) 
         nbeams = ddata(61)
         if (ddata(62)>0.5) then
            read_rst_file = .true. 
         else 
            read_rst_file = .false. 
         endif
         rst_timestep = ddata(63)
         if (ddata(64)>0.5) then
            dump_rst_file = .true. 
         else
            dump_rst_file = .false. 
         endif
         dfrst = ddata(65)  
         if (ddata(66)>0.5) then
            laser_on = .true. 
         else
            laser_on = .false. 
         endif
         cwp = ddata(67)
         dfchi = ddata(68); dfchislice = ddata(69)
         chix0 = ddata(70); chiy0 = ddata(71); chiz0 = ddata(72)
         dfvp = ddata(73); dfvpslice = ddata(74)
         vpx0 = ddata(75); vpy0 = ddata(76); vpz0 = ddata(77)
         dfjp = ddata(78); dfjpslice = ddata(79)
         jpx0 = ddata(80); jpy0 = ddata(81); jpz0 = ddata(82)
         dfjb = ddata(83); dfjbslice = ddata(84)
         jbx0 = ddata(85); jby0 = ddata(86); jbz0 = ddata(87)
         dfe = ddata(88); dfeslice = ddata(89)
         ex0 = ddata(90); ey0 = ddata(91); ez0 = ddata(92)
         dfb = ddata(93); dfbslice = ddata(94)
         bx0 = ddata(95); by0 = ddata(96); bz0 = ddata(97)
         ipbc = ddata(98)
         if (ddata(99)>0.5) then
            dump_pha_beam = .true.
         else
            dump_pha_beam = .false.
         endif
         dfpha_beam = ddata(100)
         dsample_beam = ddata(101)
         verbose = ddata(102)
         ntotal = ddata(103)
         psolve = ddata(104)
         

         if (nbeams > 0) then
            allocate(bevolutions(nbeams),min_bparticles(nbeams),npx3(nbea&
            &ms),npy3(nbeams),npz3(nbeams),bgamma(nbeams),charge3s(nbeams),mass&
            &3s(nbeams),vdx3(nbeams),vdy3(nbeams),vdz3(nbeams),iroutine(nbeams)&
            &,beam_profiles(nbeams),bquiets(nbeams),para_arr(5,100,nbeams),&
            &bp_sizes(3,nbeams),buse_shifters(nbeams),bshifter_ns&
            &ecs(nbeams),bshifter_paras(3,100,nbeams),buse_destroyers(nbeams),b&
            &destroyer_ncriterias(nbeams),bdestroyer_criterias(3,100,nbeams),bu&
            &se_radiation_dampings(nbeams),track(nbeams),track_num(nbeams))
         endif

         do i = 1, nbeams
            if (p%getidproc() == 0) then
               read(8, beam)
               track(i)=beam_track
               track_num(i)=beam_track_num
               bevolutions(i)=beam_evolution
               min_bparticles(i)=min_beam_particle
               npx3(i)=npx
               npy3(i)=npy
               npz3(i)=npz
               ntemp = num_particle*1e12*(2**indz)
               ntemp = ntemp*(2**indx)
               ntemp = ntemp*(2**indy)/(npx*alx*aly*alz*plasma_density*cwp*cwp*cwp) 
               ntemp = ntemp/npy
               ntemp = ntemp/npz
               charge3s(i) = ntemp*charge
               mass3s(i) = ntemp*mass
               bgamma(i) = gamma
               vdx3(i)=vdx*gamma
               vdy3(i)=vdy*gamma
               vdz3(i)=gamma
               vtdx=0.0
               vtdy=0.0
               vtdz=0.0 
               iroutine(i) = init_routine
               beam_profiles(i) = beam_profile 
               bquiets(i) = quiet_start
               buse_radiation_dampings(i)= use_radiation_damping
               buse_shifters(i) = use_shifter
               bshifter_nsecs(i) = shifter_nsec
               if (use_shifter) then
                  shifter_parameter(:,:)=shifter_parameter(:,:)/cwp
               endif
               buse_destroyers(i) = use_destroyer
               bdestroyer_ncriterias(i) = destroyer_ncriteria
               if (use_destroyer) then
                  do k = 1, destroyer_ncriteria
                     if (destroyer_criteria(1,k)<4) then
                        destroyer_criteria(2:3,k)=destroyer_criteria(2:3,k)/cwp
                     endif
                  end do
               end if
               select case (iroutine(i))
               case (1)
! for random tri-gaussian initialization 
! para_arr(1,1:3) = centerx,centery,centerz
! para_arr(2,1:3) = sigx, sigy, sigz
! para_arr(3,1:3) = vtx, vty, vtz
! para_arr(4,1:3) = c2x,c1x,c0x
! para_arr(5,1:3) = c2y,c1y,c0y
                  parameter_array(3,1) = parameter_array(3,1)/parameter_array(2,1)
                  parameter_array(3,2) = parameter_array(3,2)/parameter_array(2,2)
                  parameter_array(3,3) = parameter_array(3,3)*gamma
                  parameter_array(2,1:3) = parameter_array(2,1:3)/cwp
                  parameter_array(1,1:3) = parameter_array(1,1:3)/cwp
                  parameter_array(4,1) = parameter_array(4,1)*cwp 
                  parameter_array(4,3) = parameter_array(4,3)/cwp
                  parameter_array(5,1) = parameter_array(5,1)*cwp 
                  parameter_array(5,3) = parameter_array(5,3)/cwp
                  parameter_array(2,1) = parameter_array(2,1)/dx
                  parameter_array(1,1) = parameter_array(1,1)/dx
                  parameter_array(2,2) = parameter_array(2,2)/dy
                  parameter_array(1,2) = parameter_array(1,2)/dy
                  parameter_array(2,3) = parameter_array(2,3)/dz
                  parameter_array(1,3) = parameter_array(1,3)/dz
                  parameter_array(4,1) = parameter_array(4,1)*dz*dz/dx
                  parameter_array(4,2) = parameter_array(4,2)*dz/dx
                  parameter_array(4,3) = parameter_array(4,3)/dz*dz/dx
                  parameter_array(5,1) = parameter_array(5,1)*dz*dz/dy
                  parameter_array(5,2) = parameter_array(5,2)*dz/dy
                  parameter_array(5,3) = parameter_array(5,3)/dz*dz/dy 
               end select
               para_arr(:,:,i) = parameter_array(:,:)
               
               ddata(1) = npx3(i); ddata(2) = npy3(i); ddata(3) = npz3(i)
               ddata(4) = charge3s(i); ddata(5) = mass3s(i)
               ddata(6) = vdx3(i); ddata(7) = vdy3(i); ddata(8) = vdz3(i)
               ddata(9)=bgamma(i); ddata(10) = iroutine(i)
               ddata(11) = min_bparticles(i)
               if (bevolutions(i)) then
                  ddata(12) = 1.0
               else
                  ddata(12) = 0.0
               endif
               if (bquiets(i)) then
                  ddata(13) = 1.0
               else
                  ddata(13) = 0.0
               endif 
               if (buse_shifters(i)) then
                  ddata(14) = 1.0
               else
                  ddata(14) = 0.0
               endif
               ddata(15) = bshifter_nsecs(i)
               if (buse_destroyers(i)) then
                  ddata(16) = 1.0
               else
                  ddata(16) = 0.0
               endif
               ddata(17) = bdestroyer_ncriterias(i)
               if (buse_radiation_dampings(i)) then
                  ddata(18) = 1.0
               else
                  ddata(18) = 0.0
               endif
               if (track(i)) then
                  ddata(19) = 1.0
               else
                  ddata(19) = 0.0
               endif
               ddata(20) = track_num(i)
            end if

            call MPI_BCAST(ddata, 20, this%p%getmreal(), 0,&
            &this%p%getlworld(), ierr)

            npx3(i) = ddata(1); npy3(i) = ddata(2); npz3(i) = ddata(3)
            charge3s(i) = ddata(4); mass3s(i) = ddata(5)
            vdx3(i) = ddata(6); vdy3(i) = ddata(7); vdz3(i) = ddata(8)
            bgamma(i) = ddata(9); iroutine(i) = ddata(10)
            min_bparticles(i) = ddata(11)
            if (ddata(12) > 0.5) then
               bevolutions(i) = .true.
            else
               bevolutions(i) = .false.
            endif
            if (ddata(13) > 0.5) then
               bquiets(i) = .true.
            else
               bquiets(i) = .false.
            endif
            if (ddata(14) > 0.5) then
               buse_shifters(i) = .true.
            else
               buse_shifters(i) = .false.
            endif
            bshifter_nsecs(i) = ddata(15)
            if (ddata(16) > 0.5) then
               buse_destroyers(i) = .true.
            else
               buse_destroyers(i) = .false.
            endif
            bdestroyer_ncriterias(i) = ddata(17)
            if (ddata(18) > 0.5) then
               buse_radiation_dampings(i) = .true.
            else
               buse_radiation_dampings(i) = .false.
            endif
            if (ddata(19) > 0.5) then
               track(i) = .true.
            else
               track(i) = .false.
            endif
            
            track_num(i) = ddata(20)
            
            if (p%getidproc() == 0) then
               ddata = reshape(para_arr(:,:,i),(/500/))
            end if
            call MPI_BCAST(ddata, 500, p%getmreal(), 0, p%getlworld(), ierr)
            para_arr(:,:,i) = reshape(ddata,(/5,100/))
         end do
         if (p%getidproc() == 0) then
               rewind(8)
         end if

         allocate(charge2s(ntotal),mass2s(ntotal),load_balance_ths(ntot&
         &al), np2s(ntotal),vt2xs(ntotal),                            &
         &vt2ys(ntotal),non_neutral_factors(ntotal),                    &
         &profile_types(ntotal),argx1s(ntotal),argx2s(ntotal),        &
         &argx3s(ntotal), argx4s(ntotal), argx5s(ntotal), nsrands(ntotal&
         &),prof_nsecs(ntotal), prof_para&
         &s(ntotal,2,100), den_vars(ntotal)                             &
         &,den_var_nsecs(ntotal),den_var_fss(ntotal,100),den_var_ss(ntot&
         &al,100))
         allocate(neutral_species(nneutrals),ionization_level(nneutrals))
 

         if (p%getidproc() == 0) then
            do i = 1, ntotal
               read (8,species)
               charge2s(i) = charge
               mass2s(i) = mass
               np2s(i)=np2
               load_balance_ths(i) = load_balance_th
               vt2xs(i)=vt2x/dx
               vt2ys(i)=vt2y/dy
               non_neutral_factors(i)=non_neutral_factor
               profile_types(i)=profile_type
               argx1s(i)=argx1
               argx3s(i)=argx3
               argx2s(i)=argx2
               argx4s(i)=argx4
               argx5s(i)=argx5
               nsrands(i)=nsrand
               prof_nsecs(i)=prof_nsec
               prof_paras(i,1:2,1:100) = prof_parameter
               prof_paras(i,2,1:100) = prof_paras(i,2,1:100)/cwp/dx 
               den_vars(i)=density_variation
               den_var_nsecs(i)=density_variation_nsec
               den_var_fss(i,1:100)=density_variation_fs(1:100)
               den_var_ss(i,1:100)=density_variation_s(1:100)/cwp/dt
               select case (profile_types(i))
               case (1)
                  charge2s(i) = charge2s(i)/(real(np2s(i))**2)*real(2**indx&
                  &)*real(2**indy)
                  mass2s(i) = mass2s(i)/(real(np2s(i))**2)*real(2**indx&
                  &)*real(2**indy)
               end select
            enddo
            rewind(8)
            do i = 1, nneutrals
               read (8,neutral)
               neutral_species(i) = neutral_gas
               ionization_level(i) = neutral_z
            enddo
         end if

         do i = 1, ntotal
            if (p%getidproc() == 0) then
               ddata(1) = charge2s(i)
               ddata(2) = mass2s(i)
               ddata(3) = np2s(i)
               ddata(4) = load_balance_ths(i)
               ddata(5) = vt2xs(i)
               ddata(6) = vt2ys(i)
               ddata(7) = non_neutral_factors(i)
               ddata(8) = profile_types(i)
               ddata(9) = argx1s(i)
               ddata(10) = argx3s(i)
               ddata(11) = argx2s(i)
               ddata(12) = argx4s(i)
               ddata(13) = argx5s(i)
               ddata(14) = nsrands(i)
               ddata(15) = prof_nsecs(i)
               if (den_vars(i)) then
                  ddata(16) = 1
               else
                  ddata(16) = 0
               end if
               ddata(17) = den_var_nsecs(i)
               ddata(18:117) = reshape(den_var_fss(i,1:100),(/100/))
               ddata(118:217) = reshape(den_var_ss(i,1:100),(/100/))
            end if
            call MPI_BCAST(ddata, 217, p%getmreal(), 0, p%getlworld(), ierr)
            charge2s(i) = ddata(1) 
            mass2s(i) = ddata(2)
            np2s(i) = ddata(3) 
            load_balance_ths(i) = ddata(4) 
            vt2xs(i) = ddata(5) 
            vt2ys(i) = ddata(6) 
            non_neutral_factors(i) = ddata(7) 
            profile_types(i) = ddata(8) 
            argx1s(i) = ddata(9) 
            argx3s(i) = ddata(10) 
            argx2s(i) = ddata(11) 
            argx4s(i) = ddata(12) 
            argx5s(i) = ddata(13) 
            nsrands(i) = ddata(14) 
            prof_nsecs(i) = ddata(15) 
            if (ddata(16) > 0.5) then
               den_vars(i) = .true.
            else
               den_vars(i) = .false.
            end if
            den_var_nsecs(i) = ddata(17)
            den_var_fss(i,1:100) = reshape(ddata(18:117),(/100/))
            den_var_ss(i,1:100) = reshape(ddata(118:217),(/100/))
         enddo
         
         do i = 1, nneutrals
            if (p%getidproc() == 0) then
               ddata(1) = neutral_species(i)
               ddata(2) = ionization_level(i)
            end if
            call MPI_BCAST(ddata, 2, p%getmreal(), 0, p%getlworld(), ierr)
            neutral_species(i) = ddata(1)
            ionization_level(i) = ddata(2)
         end do
         
         call pp%new(nst=num_stages)
         
         this%pp => pp
         
         this%sim = simsys_in(indx,indy,indz,tend,dt,0.912871,0.912871,0.912871,&
         &alx,aly,alz,verbose,psolve,num_stages,cwp,plasma_density,dx,dy,dz,max_iter,&
         &nbeams,ntotal)
         
         this%res = restart_in(read_rst_file,dump_rst_file,rst_timestep,dfrst)
                  
         this%diag = diag_in(dfpsi,dfqep,dfqeb,dfjp,dfe,dfb,dfpsislice,dfqebslice,&
         &dfqepslice,dfjpslice,dfjbslice,dfeslice,dfbslice,int(psix0),int(psiy0),&
         &int(psiz0),int(qebx0),int(qeby0),int(qebz0),int(qepx0),int(qepy0),int(qepz0),&
         &int(jpx0),int(jpy0),int(jpz0),int(ex0),int(ey0),int(ez0),int(bx0),int(by0),&
         &int(bz0),dump_pha_beam,dfpha_beam,dsample_beam)
         
         allocate(this%beam(nbeams),this%species(ntotal))
         
         do i = 1, nbeams
            select case (iroutine(i))
            case (1)
               barg(1,1) = para_arr(3,1,i); barg(2,1) = para_arr(3,2,i)
               barg(3,1) = para_arr(3,3,i)
               barg(1,2) = 0.0; barg(2,2) = 0.0; barg(3,2) = bgamma(i)
               barg(1,3) = para_arr(2,1,i); barg(2,3) = para_arr(2,2,i)
               barg(3,3) = para_arr(2,3,i)
               barg(1,4) = para_arr(1,1,i); barg(2,4) = para_arr(1,2,i)
               barg(3,4) = para_arr(1,3,i)
               barg(1,5:7) = para_arr(4,1:3,i); barg(2,5:7) = para_arr(5,1:3,i)
               if (bquiets(i)) then
                  barg(1,8) = 1.0
               else
                  barg(1,8) = 0.0
               endif
            end select
            this%beam(i) = beam_in(min_bparticles(i),npx3(i),npy3(i),npz3(i),&
            &barg,charge3s(i),mass3s(i),iroutine(i))
         end do
         
         do i = 1, ntotal
            select case (profile_types(i))
            case (1)
               npmax = 4 * np2s(i) * np2s(i)/pp%getlnvp()
            end select
            this%species(i) = species_in(np2s(i),np2s(i),npmax,barg,profile_types(i),&
            &charge2s(i),mass2s(i))
         end do         
   
         call sp%new(this%pp,this%err,indx,indy,indz,psolve,1)

         this%sp => sp
         
         call this%err%werrfl0(class//sname//' ended')
         
         call err%setmonitor(verbose)

      end subroutine read_input
!      
      subroutine read_input_json(this)
    
         implicit none
         
         class(input_json), intent(inout) :: this
! local data
         character(len=18), save :: sname = 'read_input_json:'
         logical :: found
         character(len=:), allocatable :: ff, boundary
         integer :: length, num_stages, verbose, indx, indy, indz, psolve, ierr
         
         call p%new()
         
         this%p => p

         call err%new(this%p,2,monitor=0)
         
         this%err => err
         
         call this%err%werrfl0(class//sname//' started')
         
         call this%initialize(comment_char='/')
         
         if (p%getidproc() == 0) then
            inquire(FILE='./qpinput.json', EXIST=found)
            if(found) then
! read the file
               call this%load_file(filename = './qpinput.json')
            else
                write (erstr,*) 'error: cannot find the input file'
                call this%err%equit(class//sname//erstr)
                return
            end if
            call this%print_to_string(ff)
            length = len(ff)
         end if

         call MPI_BCAST(length, 1, this%p%getmint(), 0, this%p%getlworld(), ierr)
         
         if (.not. allocated(ff)) allocate(character(len=length) :: ff)
         
         call MPI_BCAST(ff, length, this%p%getmchar(), 0, this%p%getlworld(), ierr)
         
         call this%load_from_string(ff)
         
         call this%get('simulation.nodes(2)',num_stages,found)
         if (.not. found) then
             write (erstr,*) 'error: cannot find simulation.nodes(2) in the input file'
             call this%err%equit(class//sname//erstr)
             return
         end if
         
         call pp%new(nst=num_stages)
         
         this%pp => pp
                     
         call this%get('simulation.indx',indx,found)
         if (.not. found) then
             write (erstr,*) 'error: cannot find simulation.indx in the input file'
             call this%err%equit(class//sname//erstr)
             return
         end if

         call this%get('simulation.indy',indy,found)
         if (.not. found) then
             write (erstr,*) 'error: cannot find simulation.indy in the input file'
             call this%err%equit(class//sname//erstr)
             return
         end if

         call this%get('simulation.indz',indz,found)
         if (.not. found) then
             write (erstr,*) 'error: cannot find simulation.indz in the input file'
             call this%err%equit(class//sname//erstr)
             return
         end if

         call this%get('simulation.boundary',boundary,found)
         if (.not. found) then
             write (erstr,*) 'error: cannot find simulation.boundary in the input file'
             call this%err%equit(class//sname//erstr)
             return
         end if

         select case (trim(boundary))
         case ("conducting")
            psolve = 1
         case default
            psolve = 1
         end select
         
         call sp%new(this%pp,this%err,indx,indy,indz,psolve,1)

         this%sp => sp

         call this%err%werrfl0(ff)
         
         call this%err%werrfl0(class//sname//' ended')
         
         call this%get('simulation.verbose',verbose,found)
         if (.not. found) then
             write (erstr,*) 'error: cannot find simulation.verbose in the input file'
             call this%err%equit(class//sname//erstr)
             return
         end if

         call err%setmonitor(verbose)

      end subroutine read_input_json
!
      subroutine initialize(this,verbose,compact_reals,print_signs,&
      &real_format,spaces_per_tab,strict_type_checking,trailing_spaces_significant,&
      &case_sensitive_keys,no_whitespace,unescape_strings,comment_char,path_mode,&
      &path_separator,compress_vectors,allow_duplicate_keys)
    
         implicit none
         
         class(input_json), intent(inout) :: this
         logical,intent(in),optional :: verbose
         logical,intent(in),optional :: compact_reals
         logical,intent(in),optional :: print_signs
         character(len=*),intent(in),optional :: real_format
         integer,intent(in),optional :: spaces_per_tab
         logical,intent(in),optional :: strict_type_checking
         logical,intent(in),optional :: trailing_spaces_significant
         logical,intent(in),optional :: case_sensitive_keys
         logical,intent(in),optional :: no_whitespace
         logical,intent(in),optional :: unescape_strings
         character(len=1),intent(in),optional :: comment_char
         integer,intent(in),optional :: path_mode
         character(len=1),intent(in),optional :: path_separator
         logical,intent(in),optional :: compress_vectors
         logical,intent(in),optional :: allow_duplicate_keys
! local data
         character(len=38), save :: sname = 'initialize_json_core_in_file:'
         call this%err%werrfl0(class//sname//' started')
         
         call this%input%initialize(verbose,compact_reals,print_signs,&
         &real_format,spaces_per_tab,strict_type_checking,trailing_spaces_significant,&
         &case_sensitive_keys,no_whitespace,unescape_strings,comment_char,path_mode,&
         &path_separator,compress_vectors,allow_duplicate_keys)

         call this%err%werrfl0(class//sname//' ended')
      
      end subroutine initialize
!
      subroutine set_json_core_in_file(this,core)

         implicit none

         class(input_json),intent(inout) :: this
         type(json_core),intent(in) :: core
! local data
         character(len=38), save :: sname = 'set_json_core_in_file:'
         call this%err%werrfl0(class//sname//' started')

         call this%input%initialize(core)

         call this%err%werrfl0(class//sname//' ended')
         
      end subroutine set_json_core_in_file   
!
      subroutine load_file(this, filename, unit)

         implicit none

         class(input_json),intent(inout) :: this
         character(len=*),intent(in) :: filename
         integer,intent(in),optional :: unit
! local data
         character(len=18), save :: sname = 'json_file_load:'
         call this%err%werrfl0(class//sname//' started')

         call this%input%load_file(filename, unit)

         call this%err%werrfl0(class//sname//' ended')
         
      end subroutine load_file
!
      subroutine print_to_string(this, str)

         implicit none

         class(input_json),intent(inout) :: this
         character(len=:),allocatable,intent(out) :: str
! local data
         character(len=38), save :: sname = 'json_file_print_to_string:'
         call this%err%werrfl0(class//sname//' started')

         call this%input%print_to_string(str)

         call this%err%werrfl0(class//sname//' ended')

      end subroutine print_to_string     
!
      subroutine load_from_string(this, str)

         implicit none

         class(input_json),intent(inout) :: this
         character(len=*),intent(in) :: str
! local data
         character(len=38), save :: sname = 'json_file_load_from_string:'
         call this%err%werrfl0(class//sname//' started')

         call this%input%load_from_string(str)

         call this%err%werrfl0(class//sname//' ended')

      end subroutine load_from_string
!
      subroutine json_file_get_object(this, path, p, found)

         implicit none

         class(input_json),intent(inout) :: this
         character(len=*),intent(in) :: path
         type(json_value),pointer,intent(out) :: p
         logical,intent(out),optional :: found
! local data
         character(len=38), save :: sname = 'json_file_get_object:'
         call this%err%werrfl0(class//sname//' started')

         call this%input%get(path, p, found)

         call this%err%werrfl0(class//sname//' ended')

      end subroutine json_file_get_object
!
      subroutine json_file_get_integer(this, path, val, found)

         implicit none

         class(input_json),intent(inout) :: this
         character(len=*),intent(in) :: path
         integer,intent(out) :: val 
         logical,intent(out),optional :: found
! local data
         character(len=38), save :: sname = 'json_file_get_integer:'
         call this%err%werrfl0(class//sname//' started')

         call this%input%get(path, val, found)

         call this%err%werrfl0(class//sname//' ended')

      end subroutine json_file_get_integer
!
      subroutine json_file_get_double(this, path, val, found)

         implicit none

         class(input_json),intent(inout) :: this
         character(len=*),intent(in) :: path
         real,intent(out) :: val 
         logical,intent(out),optional :: found
! local data
         character(len=38), save :: sname = 'json_file_get_double:'
         call this%err%werrfl0(class//sname//' started')

         call this%input%get(path, val, found)

         call this%err%werrfl0(class//sname//' ended')

      end subroutine json_file_get_double
!
      subroutine json_file_get_logical(this, path, val, found)

         implicit none

         class(input_json),intent(inout) :: this
         character(len=*),intent(in) :: path
         logical,intent(out) :: val 
         logical,intent(out),optional :: found
! local data
         character(len=38), save :: sname = 'json_file_get_logical:'
         call this%err%werrfl0(class//sname//' started')

         call this%input%get(path, val, found)

         call this%err%werrfl0(class//sname//' ended')

      end subroutine json_file_get_logical
!
      subroutine json_file_get_string(this, path, val, found)

         implicit none

         class(input_json),intent(inout) :: this
         character(len=*),intent(in) :: path
         character(len=:),allocatable,intent(out) :: val 
         logical,intent(out),optional :: found
! local data
         character(len=38), save :: sname = 'json_file_get_string:'
         call this%err%werrfl0(class//sname//' started')

         call this%input%get(path, val, found)

         call this%err%werrfl0(class//sname//' ended')

      end subroutine json_file_get_string
!
      subroutine json_file_get_integer_vec(this, path, vec, found)

         implicit none

         class(input_json),intent(inout) :: this
         character(len=*),intent(in) :: path
         integer,dimension(:),allocatable,intent(out) :: vec 
         logical,intent(out),optional :: found
! local data
         character(len=38), save :: sname = 'json_file_get_integer_vec:'
         call this%err%werrfl0(class//sname//' started')

         call this%input%get(path, vec, found)

         call this%err%werrfl0(class//sname//' ended')

      end subroutine json_file_get_integer_vec
!
      subroutine json_file_get_double_vec(this, path, vec, found)

         implicit none

         class(input_json),intent(inout) :: this
         character(len=*),intent(in) :: path
         real,dimension(:),allocatable,intent(out) :: vec 
         logical,intent(out),optional :: found
! local data
         character(len=38), save :: sname = 'json_file_get_double_vec:'
         call this%err%werrfl0(class//sname//' started')

         call this%input%get(path, vec, found)

         call this%err%werrfl0(class//sname//' ended')

      end subroutine json_file_get_double_vec
!
      subroutine json_file_get_logical_vec(this, path, vec, found)

         implicit none

         class(input_json),intent(inout) :: this
         character(len=*),intent(in) :: path
         logical,dimension(:),allocatable,intent(out) :: vec 
         logical,intent(out),optional :: found
! local data
         character(len=38), save :: sname = 'json_file_get_logical_vec:'
         call this%err%werrfl0(class//sname//' started')

         call this%input%get(path, vec, found)

         call this%err%werrfl0(class//sname//' ended')

      end subroutine json_file_get_logical_vec
!
      subroutine json_file_get_string_vec(this, path, vec, found)

         implicit none

         class(input_json),intent(inout) :: this
         character(len=*),intent(in) :: path
         character(len=*),dimension(:),allocatable,intent(out) :: vec 
         logical,intent(out),optional :: found
! local data
         character(len=38), save :: sname = 'json_file_get_string_vec:'
         call this%err%werrfl0(class//sname//' started')

         call this%input%get(path, vec, found)

         call this%err%werrfl0(class//sname//' ended')

      end subroutine json_file_get_string_vec
!
      subroutine json_file_get_alloc_string_vec(this, path, vec, ilen, found)

         implicit none

         class(input_json),intent(inout) :: this
         character(len=*),intent(in) :: path
         character(len=:),dimension(:),allocatable,intent(out) :: vec
         integer,dimension(:),allocatable,intent(out) :: ilen 
         logical,intent(out),optional :: found
! local data
         character(len=38), save :: sname = 'json_file_get_alloc_string_vec:'
         call this%err%werrfl0(class//sname//' started')

         call this%input%get(path, vec, ilen, found)

         call this%err%werrfl0(class//sname//' ended')

      end subroutine json_file_get_alloc_string_vec
!
      subroutine json_file_get_root(this,p)

         implicit none

         class(input_json),intent(inout) :: this
         type(json_value),pointer,intent(out) :: p
! local data
         character(len=38), save :: sname = 'json_file_get_root:'
         call this%err%werrfl0(class//sname//' started')

         call this%input%get(p)

         call this%err%werrfl0(class//sname//' ended')

      end subroutine json_file_get_root                            
!
      subroutine  json_file_variable_info(this,path, found, var_type, n_children, name)

         implicit none

         class(input_json),intent(inout) :: this
         character(len=*),intent(in) :: path
         logical,intent(out),optional :: found
         integer,intent(out),optional :: var_type
         integer,intent(out),optional :: n_children
         character(len=:),allocatable,intent(out),optional :: name

! local data
         character(len=38), save :: sname = 'json_file_variable_info:'
         call this%err%werrfl0(class//sname//' started')

         call this%input%info(path, found, var_type, n_children, name)

         call this%err%werrfl0(class//sname//' ended')

      end subroutine json_file_variable_info                           

      end module input_class