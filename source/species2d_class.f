! species2d class for QuickPIC Open Source 1.0
! update: 04/18/2016

      module species2d_class

      use perrors_class
      use parallel_pipe_class
      use spect2d_class
      use fdist2d_class
      use field2d_class
      use part2d_class
      use hdf5io_class
               
      implicit none

      private

      public :: species2d

      type species2d

         private

         class(spect2d), pointer, public :: sp => null()
         class(perrors), pointer, public :: err => null()
         class(parallel_pipe), pointer, public :: p => null()
         class(part2d), pointer :: pd => null()
         
         
         contains
         
         generic :: new => init_species2d
         generic :: renew => renew_species2d
         generic :: del => end_species2d
         generic :: qdp => qdp_species2d
         generic :: amjdp => amjdp_species2d
         generic :: push => push_species2d
         generic :: pmv => pmove_species2d
         generic :: extpsi => extpsi_species2d
         generic :: pcp => pcp_species2d
         generic :: pcb => pcb_species2d
         generic :: psend => psend_species2d
         generic :: precv => precv_species2d
         generic :: wr => writehdf5_species2d
         procedure, private :: init_species2d, renew_species2d
         procedure, private :: end_species2d
         procedure, private :: qdp_species2d
         procedure, private :: amjdp_species2d
         procedure, private :: push_species2d
         procedure, private :: pmove_species2d
         procedure, private :: extpsi_species2d
         procedure, private :: pcp_species2d
         procedure, private :: pcb_species2d
         procedure, private :: psend_species2d
         procedure, private :: precv_species2d, writehdf5_species2d
                           
      end type 

      save      

      character(len=10) :: class = 'species2d:'
      character(len=128) :: erstr
      
      contains
!
      subroutine init_species2d(this,pp,perr,psp,pf,fd,qm,qbm,dt,ci,xdim,&
      &npmax,nbmax)
      
         implicit none
         
         class(species2d), intent(inout) :: this
         class(spect2d), intent(in), pointer :: psp
         class(perrors), intent(in), pointer :: perr
         class(parallel_pipe), intent(in), pointer :: pp
         class(fdist2d), intent(in) :: pf
         class(field2d), intent(in) :: fd
         real, intent(in) :: qm, qbm, dt, ci
         integer, intent(in) :: npmax, nbmax, xdim

! local data
         character(len=18), save :: sname = 'init_species2d:'
                  
         this%sp => psp
         this%err => perr
         this%p => pp

         call this%err%werrfl2(class//sname//' started')
         
         allocate(this%pd)
         
         call this%pd%new(pp,perr,psp,pf,fd%getrs(),qm,qbm,dt,ci,xdim,npmax,nbmax)
         
         call this%err%werrfl2(class//sname//' ended')

      end subroutine init_species2d
!
      subroutine end_species2d(this)
          
         implicit none
         
         class(species2d), intent(inout) :: this
         character(len=18), save :: sname = 'end_species2d:'

         call this%err%werrfl2(class//sname//' started')
         call this%pd%del()
         deallocate(this%pd)
         call this%err%werrfl2(class//sname//' ended')
                  
      end subroutine end_species2d
!
      subroutine renew_species2d(this,pf,fd)
      
         implicit none
         
         class(species2d), intent(inout) :: this
         class(fdist2d), intent(in) :: pf
         class(field2d), intent(in) :: fd

! local data
         character(len=18), save :: sname = 'renew_species2d:'
                  
         call this%err%werrfl2(class//sname//' started')
         
         if (this%p%getstageid() == 0) then
            call this%pd%renew(pf,fd%getrs())
         end if
         
         call this%err%werrfl2(class//sname//' ended')

      end subroutine renew_species2d
!      
      subroutine qdp_species2d(this,q)
! deposit the charge density      
      
         implicit none
         
         class(species2d), intent(in) :: this
         class(field2d), intent(inout) :: q
! local data
         character(len=18), save :: sname = 'qdp_species2d:'
         
                  
         call this%err%werrfl2(class//sname//' started')
                  
         call this%pd%qdp(q%getrs())
         
         call q%ag()
                  
         call this%err%werrfl2(class//sname//' ended')
         
      end subroutine qdp_species2d
!      
      subroutine amjdp_species2d(this,ef,bf,psit,cu,amu,dcu,dex)
! deposit the current, acceleration and momentum flux      
      
         implicit none
         
         class(species2d), intent(inout) :: this
         class(field2d), intent(inout) :: cu, amu, dcu
         class(field2d), intent(in) :: ef, bf, psit
         real, intent(in) :: dex
! local data
         character(len=18), save :: sname = 'amjdp_species2d'

         call this%err%werrfl2(class//sname//' started')
         
         call this%pd%amjdp(ef%getrs(),bf%getrs(),psit%getrs(),cu%getrs(),&
         &amu%getrs(),dcu%getrs(),dex)
         
         call cu%ag()
         call dcu%ag()
         call amu%ag()

         call this%err%werrfl2(class//sname//' ended')
         
      end subroutine amjdp_species2d
!      
      subroutine push_species2d(this,ef,bf,psit,dex)
      
         implicit none
         
         class(species2d), intent(inout) :: this
         class(field2d), intent(in) :: ef, bf, psit
         real, intent(in) :: dex
! local data
         character(len=18), save :: sname = 'push_species2d'

         call this%err%werrfl2(class//sname//' started')

         call this%pd%push(ef%getrs(),bf%getrs(),psit%getrs(),dex)         
         call this%pmv(psit)
         
         call this%err%werrfl2(class//sname//' ended')
         
      end subroutine push_species2d
!
      subroutine pmove_species2d(this,fd)
      
         implicit none
         
         class(species2d), intent(inout) :: this
         class(field2d), intent(in) :: fd
! local data
         character(len=18), save :: sname = 'pmove_species2d:'
         
         call this%err%werrfl2(class//sname//' started')
         
         call this%pd%pmv(fd%getrs())

         call this%err%werrfl2(class//sname//' ended')
         
      end subroutine pmove_species2d
!      
      subroutine extpsi_species2d(this,psi,dex)
      
         implicit none
         
         class(species2d), intent(inout) :: this
         class(field2d), intent(in) :: psi
         real, intent(in) :: dex
! local data
         character(len=18), save :: sname = 'extpsi_species2d:'
         
         call this%err%werrfl2(class//sname//' started')
         
         call this%pd%extpsi(psi%getrs(),dex)

         call this%err%werrfl2(class//sname//' ended')
         
      end subroutine extpsi_species2d
!
      subroutine pcp_species2d(this,fd)
      
         implicit none
         
         class(species2d), intent(inout) :: this
         class(field2d), intent(in) :: fd
! local data
         character(len=18), save :: sname = 'pcp_species2d:'
         
         call this%err%werrfl2(class//sname//' started')
         
         call this%pd%pcp(fd%getrs())

         call this%err%werrfl2(class//sname//' ended')
         
      end subroutine pcp_species2d
!
      subroutine pcb_species2d(this)
      
         implicit none
         
         class(species2d), intent(inout) :: this
! local data
         character(len=18), save :: sname = 'pcb_species2d:'
         
         call this%err%werrfl2(class//sname//' started')
         
         call this%pd%pcb()

         call this%err%werrfl2(class//sname//' ended')
         
      end subroutine pcb_species2d
!
      subroutine psend_species2d(this,tag,id)
      
         implicit none
         
         class(species2d), intent(inout) :: this
         integer, intent(in) :: tag
         integer, intent(inout) :: id
! local data
         character(len=18), save :: sname = 'pipesend_part2d:'
                  
         call this%err%werrfl2(class//sname//' started')
         
         call this%pd%psend(tag,id)
         
         call this%err%werrfl2(class//sname//' ended')
         
      end subroutine psend_species2d
!      
      subroutine precv_species2d(this,fd,tag)
      
         implicit none
         
         class(species2d), intent(inout) :: this
         class(field2d), intent(in) :: fd
         integer, intent(in) :: tag
! local data
         character(len=18), save :: sname = 'precv_species2d:'
         
         
         call this%err%werrfl2(class//sname//' started')
         
         call this%pd%precv(fd%getrs(),tag)
                  
         call this%err%werrfl2(class//sname//' ended')
         
      end subroutine precv_species2d
!      
      subroutine writehdf5_species2d(this,file,delta)

         implicit none
         
         class(species2d), intent(inout) :: this
         class(hdf5file), intent(in) :: file
         real, dimension(2), intent(in) :: delta
! local data
         character(len=18), save :: sname = 'writehdf5_species2d:'

         call this%err%werrfl2(class//sname//' started')                  
         
         call this%pd%wr(file,delta) 

         call this%err%werrfl2(class//sname//' ended')
      
      end subroutine writehdf5_species2d
!      
      end module species2d_class