! species2d class for QuickPIC Open Source 1.0
! update: 04/18/2016

      module species2d_class

      use perrors_class
      use parallel_pipe_class
      use spect2d_class
      use spect3d_class
      use fdist2d_class
      use field2d_class
      use field3d_class
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
         class(field2d), pointer :: q => null(), qn => null(), cu => null()
         class(field2d), pointer :: amu => null(), dcu => null()
         class(field3d), pointer :: q3 => null()
         class(fdist2d), pointer :: pf => null()
         
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
         generic :: wrq => writeq_species2d, writeqslice_species2d
         generic :: cbq => cbq_species2d
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
         procedure, private :: cbq_species2d, writeq_species2d, writeqslice_species2d
                           
      end type 

      save      

      character(len=10) :: class = 'species2d:'
      character(len=128) :: erstr
      
      contains
!
      subroutine init_species2d(this,pp,perr,psp,pf,qbm,dt,ci,xdim,s)

         implicit none
         
         class(species2d), intent(inout) :: this
         class(spect3d), intent(in), pointer :: psp
         class(perrors), intent(in), pointer :: perr
         class(parallel_pipe), intent(in), pointer :: pp
         class(fdist2d), intent(inout), target :: pf
         real, intent(in) :: qbm, dt, ci, s
         integer, intent(in) :: xdim

! local data
         character(len=18), save :: sname = 'init_species2d:'
                  
         this%sp => psp
         this%err => perr
         this%p => pp
         this%pf => pf

         call this%err%werrfl2(class//sname//' started')
         
         allocate(this%pd,this%q,this%qn,this%cu,this%amu,this%dcu,this%q3)
         call this%q%new(this%p,this%err,this%sp,dim=1,fftflag=.true.)
         call this%q3%new(this%p,this%err,psp,dim=1)
         call this%qn%new(this%p,this%err,this%sp,dim=1,fftflag=.false.)
         call this%cu%new(this%p,this%err,this%sp,dim=3,fftflag=.false.)
         call this%dcu%new(this%p,this%err,this%sp,dim=2,fftflag=.false.)
         call this%amu%new(this%p,this%err,this%sp,dim=3,fftflag=.false.)
         call this%pd%new(pp,perr,this%sp,pf,this%q%getrs(),qbm,dt,ci,xdim,s)
         call this%qn%as(0.0)
         call this%cu%as(0.0)
         call this%pd%qdp(this%qn%getrs())
         call this%qn%ag()
         call this%q%as(this%qn)
         if (this%p%getstageid() == 0) then
            call this%q%fftrk(1)
            call this%q%smooth(this%q)
            call this%q%fftkr(1)
            call this%q%cb(this%q3,1,(/1/),(/1/))         
         end if
         call this%qn%mult(this%qn,-1.0)
         call this%err%werrfl2(class//sname//' ended')

      end subroutine init_species2d
!
      subroutine end_species2d(this)
          
         implicit none
         
         class(species2d), intent(inout) :: this
         character(len=18), save :: sname = 'end_species2d:'

         call this%err%werrfl2(class//sname//' started')
         call this%pd%del()
         call this%q%del()
         call this%qn%del()
         call this%cu%del()
         call this%dcu%del()
         call this%amu%del()
         call this%err%werrfl2(class//sname//' ended')
                  
      end subroutine end_species2d
!
      subroutine renew_species2d(this,s)
      
         implicit none
         
         class(species2d), intent(inout) :: this
         real, intent(in) :: s

! local data
         character(len=18), save :: sname = 'renew_species2d:'
                  
         call this%err%werrfl2(class//sname//' started')
         
         call this%pd%renew(this%pf,this%qn%getrs(),s)
         call this%qn%as(0.0)
         call this%pd%qdp(this%qn%getrs())
         call this%qn%ag()
         call this%q%as(this%qn)
         if (this%p%getstageid() == 0) then
            call this%q%fftrk(1)
            call this%q%smooth(this%q)
            call this%q%fftkr(1)
            call this%q%cb(this%q3,1,(/1/),(/1/))         
         end if
         call this%qn%mult(this%qn,-1.0)
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
         call this%q%as(0.0)
         call this%pd%qdp(this%q%getrs())
         call this%q%ag()
         call q%add(this%q,q)
         call q%add(this%qn,q)
                  
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
         call this%cu%as(0.0)
         call this%dcu%as(0.0)
         call this%amu%as(0.0)
         
         call this%pd%amjdp(ef%getrs(),bf%getrs(),psit%getrs(),this%cu%getrs(),&
         &this%amu%getrs(),this%dcu%getrs(),dex)
         
         call this%cu%ag()
         call this%dcu%ag()
         call this%amu%ag()

         call cu%add(this%cu,cu)
         call dcu%add(this%dcu,dcu)
         call amu%add(this%amu,amu)
         call this%cu%mult(this%cu,dex)

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
         
         call this%pd%pmv(this%q%getrs())

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
      subroutine precv_species2d(this,tag)
      
         implicit none
         
         class(species2d), intent(inout) :: this
         integer, intent(in) :: tag
! local data
         character(len=18), save :: sname = 'precv_species2d:'
         
         
         call this%err%werrfl2(class//sname//' started')
         
         call this%pd%precv(this%q%getrs(),tag)
                  
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
      subroutine writeq_species2d(this,file,rtag,stag,id)

         implicit none
         
         class(species2d), intent(inout) :: this
         class(hdf5file), intent(in) :: file
         integer, intent(in) :: rtag, stag
         integer, intent(inout) :: id
! local data
         character(len=18), save :: sname = 'writeq_species2d:'

         call this%err%werrfl2(class//sname//' started')                  
         
         call this%q3%wr(file,1,rtag,stag,id) 

         call this%err%werrfl2(class//sname//' ended')
      
      end subroutine writeq_species2d
!
      subroutine writeqslice_species2d(this,file,slice,spos,rtag,stag,id)

         implicit none
         
         class(species2d), intent(inout) :: this
         class(hdf5file), intent(in) :: file
         integer, intent(in) :: rtag, stag, slice, spos
         integer, intent(inout) :: id
! local data
         character(len=18), save :: sname = 'writeqslice_species2d:'

         call this%err%werrfl2(class//sname//' started')                  
         
         call this%q3%wr(file,1,slice,spos,rtag,stag,id)

         call this%err%werrfl2(class//sname//' ended')
      
      end subroutine writeqslice_species2d
!
      subroutine cbq_species2d(this,pos)

         implicit none
         
         class(species2d), intent(inout) :: this
         integer, intent(in) :: pos
! local data
         character(len=18), save :: sname = 'cpq_species2d:'

         call this%err%werrfl2(class//sname//' started')

         call this%q%add(this%q,this%cu,(/1/),(/1/),(/3/))
         call this%q%fftrk(1)
         call this%q%smooth(this%q)
         call this%q%fftkr(1)
         call this%q%cb(this%q3,pos,(/1/),(/1/))

         call this%err%werrfl2(class//sname//' ended')
      
      end subroutine cbq_species2d
!      
      end module species2d_class