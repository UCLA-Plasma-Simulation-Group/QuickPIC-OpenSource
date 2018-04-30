! field2d class for QuickPIC Open Source 1.0
! update: 04/18/2016

      module field2d_class

      use perrors_class
      use parallel_pipe_class
      use spect2d_class
      use field3d_class
      use ufield2d_class
      use ufield3d_class
      use fft2d_class
      use fpois2d_class
      use hdf5io_class
         
      implicit none

      private

      public :: field2d
      
      type field2d

         private
! state = (0,1) = current data is in (rs,ks)
! gcells = (0,1) = (no, yes) guard cell processing is performed

         class(spect2d), pointer, public :: sp => null()
         class(perrors), pointer, public :: err => null()
         class(parallel_pipe), pointer, public :: p => null()
         class(ufield2d), pointer :: rs => null(), ks => null()
         class(fft2d), pointer :: ft => null()
         class(fpois2d), pointer :: pt => null()
         integer :: state, gcells
                  
         contains
         
         generic :: new => init_field2d
         generic :: del => end_field2d
         generic :: fftrk => fftrk_field2d
         generic :: fftkr => fftkr_field2d
         generic :: div => divf_field2d
         generic :: grad => gradf_field2d
         generic :: curl => curlf_field2d
         generic :: pot => potential_field2d
         generic :: smooth => smoothf_field2d
         generic :: elf => elfield_field2d
         generic :: bf => bfield_field2d
         generic :: bfqp => bfield_qp_field2d
         generic :: cg => copyguard_field2d
         generic :: ag => acopyguard_field2d
         generic :: psend => pipesend_field2d
         generic :: precv => piperecv_field2d
         generic :: cp => copyfrom
         generic :: cb => copyto
         generic :: as => asc, asa
         generic :: add => sum1, sum2
         generic :: sub => minus1, minus2
         generic :: mult => multiply1, multiply2
         generic :: wr => writehdf5_field2d

         procedure, private :: init_field2d, end_field2d
         procedure, private :: fftrk_field2d, fftkr_field2d, divf_field2d
         procedure, private :: gradf_field2d, curlf_field2d, potential_field2d
         procedure, private :: smoothf_field2d, elfield_field2d, bfield_field2d 
         procedure, private :: bfield_qp_field2d, copyguard_field2d
         procedure, private :: acopyguard_field2d, pipesend_field2d
         procedure, private :: piperecv_field2d, asc, asa, sum1, minus1, multiply1
         procedure, private :: writehdf5_field2d, multiply2, sum2, minus2
         procedure, private :: copyfrom, copyto
         procedure :: getstate, getgcells, getrs, getks
                  
      end type 

      character(len=10), save :: class = 'field2d:'
      character(len=128), save :: erstr
      
      contains
!
      subroutine init_field2d(this,pp,perr,psp,dim,fftflag,state,gcells)

! fftflag = (.false.,.true.) = (no,yes) initial ks, fft and poisson solver table
         implicit none
         
         class(field2d), intent(inout) :: this
         class(spect2d), intent(in), pointer :: psp
         class(perrors), intent(in), pointer :: perr
         class(parallel_pipe), intent(in), pointer :: pp
         integer, intent(in) :: dim
         logical, intent(in) :: fftflag
         integer, intent(in), optional :: state,gcells
! local data
         character(len=18), save :: sname = 'init_field2d:'
         integer :: nvpy,indx,indy,nd1,nd2
         
         this%err => perr
         call this%err%werrfl2(class//sname//' started')

         this%sp => psp
         this%p => pp
         if (present(state)) then
            this%state = state
         else
            this%state = 0
         endif
         if (present(gcells)) then
            this%gcells = gcells
         else
            this%gcells = 0
         endif
         nvpy = pp%getlnvp()
         indx = psp%getindx()
         indy = psp%getindy()         
         nd1 = 2**indx
         nd2 = 2**indy
         
         allocate(this%rs)
         
         call this%rs%new(pp,perr,psp,dim,0,0,nvpy)
         
         if (fftflag) then
            allocate(this%ks)
            call this%ks%new(pp,perr,psp,dim,1,nd1,nd2,0,nvpy)
            select case (psp%getpsolver())
            case (1)
               this%ft => get_fft2table(pp,perr,psp,indx,indy)
               this%pt => get_pois2table(pp,perr,psp,ax=0.912871,ay=0.912871,affp=1.0)
            case default
               this%ft => get_fft2table(pp,perr,psp,indx,indy)
               this%pt => get_pois2table(pp,perr,psp,ax=0.912871,ay=0.912871,affp=1.0)
            end select
         endif

         call this%err%werrfl2(class//sname//' ended')

      end subroutine init_field2d
!
      subroutine end_field2d(this)
          
         implicit none
         
         class(field2d), intent(inout) :: this
         character(len=18), save :: sname = 'end_field2d:'

         call this%err%werrfl2(class//sname//' started')
         
         call this%rs%del()
         call this%ks%del()
         deallocate(this%rs,this%ks)

         call this%err%werrfl2(class//sname//' ended')
         
      end subroutine end_field2d
!
      subroutine fftrk_field2d(this,kind)
      
         implicit none
         
         class(field2d), intent(inout) :: this
         integer, intent(in) :: kind
         character(len=18), save :: sname = 'fftrk_field2d:'

         call this%err%werrfl2(class//sname//' started')
         
         if (this%state /= 0) then 
            call this%err%equit(class//sname//'data is in ks')
            return
         endif
         
         if (this%gcells == 0) call this%err%werrfl2('guard cell not processed')
         
         select case (kind)
         case (1)
            call this%ft%fsst(this%rs,this%ks,-1)
         case (2)   
            call this%ft%fcct(this%rs,this%ks,-1)
         case (3)
            if (this%rs%getdim() /= 3) then
               call this%err%equit(class//sname//'wrong fft kind')
               return
            endif
            call this%ft%fs3t(this%rs,this%ks,-1)
         end select
         
         this%state = 1
         
         call this%err%werrfl2(class//sname//' ended')
         
      end subroutine fftrk_field2d
!
      subroutine fftkr_field2d(this,kind)
      
         implicit none
         
         class(field2d), intent(inout) :: this
         integer, intent(in) :: kind
         character(len=18), save :: sname = 'fftkr_field2d:'

         call this%err%werrfl2(class//sname//' started')
         
         if (this%state /= 1) then 
            call this%err%equit(class//sname//'data is in rs')
            return
         endif
                  
         select case (kind)
         case (1)
            call this%ft%fsst(this%rs,this%ks,1)
         case (2)   
            call this%ft%fcct(this%rs,this%ks,1)
         case (3)
            if (this%rs%getdim() /= 3) then
               call this%err%equit(class//sname//'wrong fft kind')
               return
            endif
            call this%ft%fs3t(this%rs,this%ks,1)
         end select
         
         this%state = 0
         
         call this%rs%cg()
                  
         call this%err%werrfl2(class//sname//' ended')
         
      end subroutine fftkr_field2d
!
      subroutine divf_field2d(this,that)
      
         implicit none
         
         class(field2d), intent(inout) :: this,that
         character(len=18), save :: sname = 'divf_field2d:'

         call this%err%werrfl2(class//sname//' started')
         
         if (this%state /= 1) then 
            call this%err%equit(class//sname//'data is in rs')
            return
         endif
         
         call this%ft%divf(this%ks,that%ks)

         that%state = 1
         
         call this%err%werrfl2(class//sname//' ended')
         
      end subroutine divf_field2d      
!            
      subroutine gradf_field2d(this,that)
      
         implicit none
         
         class(field2d), intent(inout) :: this,that
         character(len=18), save :: sname = 'gradf_field2d:'

         call this%err%werrfl2(class//sname//' started')
         
         if (this%state /= 1) then 
            call this%err%equit(class//sname//'data is in rs')
            return
         endif
         
         call this%ft%gradf(this%ks,that%ks)

         that%state = 1
         
         call this%err%werrfl2(class//sname//' ended')
         
      end subroutine gradf_field2d      
!            
      subroutine curlf_field2d(this,that)
      
         implicit none
         
         class(field2d), intent(inout) :: this,that
         character(len=18), save :: sname = 'curlf_field2d:'

         call this%err%werrfl2(class//sname//' started')
         
         if (this%state /= 1) then 
            call this%err%equit(class//sname//'data is in rs')
            return
         endif
         
         call this%ft%curlf(this%ks,that%ks)

         that%state = 1
         
         call this%err%werrfl2(class//sname//' ended')
         
      end subroutine curlf_field2d      
!            
      subroutine potential_field2d(this,that)
      
         implicit none
         
         class(field2d), intent(inout) :: this,that
         character(len=18), save :: sname = 'potential_field2d:'
         real :: we

         call this%err%werrfl2(class//sname//' started')
         
         if (this%state /= 1) then 
            call this%err%equit(class//sname//'data is in rs')
            return
         endif
         
         call this%pt%potential(this%ks,that%ks,we)

         that%state = 1
         
         call this%err%werrfl2(class//sname//' ended')
         
      end subroutine potential_field2d      
!            
      subroutine smoothf_field2d(this,that)
      
         implicit none
         
         class(field2d), intent(inout) :: this,that
         character(len=18), save :: sname = 'smoothf_field2d:'
         real :: we

         call this%err%werrfl2(class//sname//' started')
         
         if (this%state /= 1) then 
            call this%err%equit(class//sname//'data is in rs')
            return
         endif
         
         call this%pt%smoothf(this%ks,that%ks)

         that%state = 1
         
         call this%err%werrfl2(class//sname//' ended')
         
      end subroutine smoothf_field2d      
!            
      subroutine elfield_field2d(this,that)
      
         implicit none
         
         class(field2d), intent(inout) :: this,that
         character(len=18), save :: sname = 'elfield_field2d:'
         real :: we

         call this%err%werrfl2(class//sname//' started')
         
         if (this%state /= 1) then 
            call this%err%equit(class//sname//'data is in rs')
            return
         endif
         
         call this%pt%elfield(this%ks,that%ks,we)

         that%state = 1
         
         call this%err%werrfl2(class//sname//' ended')
         
      end subroutine elfield_field2d      
!            
      subroutine bfield_field2d(this,that)
      
         implicit none
         
         class(field2d), intent(inout) :: this,that
         character(len=18), save :: sname = 'bfield_field2d:'
         real :: we

         call this%err%werrfl2(class//sname//' started')
         
         if (this%state /= 1) then 
            call this%err%equit(class//sname//'data is in rs')
            return
         endif
         
         call this%pt%bfield(this%ks,that%ks,1.0,we)

         that%state = 1
         
         call this%err%werrfl2(class//sname//' ended')
         
      end subroutine bfield_field2d      
!            
      subroutine bfield_qp_field2d(cu,dcu,amu,bxy,c,dex)
      
         implicit none
         
         class(field2d), intent(inout) :: cu,dcu,amu,bxy
         real, intent(in) :: c, dex
         character(len=18), save :: sname = 'bfield_qp_field2d:'
         real :: we

         call cu%err%werrfl2(class//sname//' started')
         
         if ((cu%state /= 1).or.(dcu%state /= 1).or.(amu%state /= 1))  then 
            call cu%err%equit(class//sname//'data is in rs')
            return
         endif
         
         call cu%pt%bfield_qp(cu%ks,dcu%ks,amu%ks,bxy%ks,1.0,c,dex,we)

         bxy%state = 1
         
         call cu%err%werrfl2(class//sname//' ended')
         
      end subroutine bfield_qp_field2d      
!
      subroutine copyguard_field2d(this) 

         implicit none
         
         class(field2d), intent(inout) :: this
         character(len=18), save :: sname = 'copyguard_field2d:'

         call this%err%werrfl2(class//sname//' started')
         if (this%state /= 0) then 
            call this%err%equit(class//sname//'data is in ks')
            return
         endif
         
         call this%rs%cg()
         
         this%gcells = 1
         this%state = 0
         
         call this%err%werrfl2(class//sname//' ended')
                 
      end subroutine copyguard_field2d
!
      subroutine acopyguard_field2d(this) 

         implicit none
         
         class(field2d), intent(inout) :: this
         character(len=18), save :: sname = 'acopyguard_field2d:'

         call this%err%werrfl2(class//sname//' started')
         if (this%state /= 0) then 
            call this%err%equit(class//sname//'data is in ks')
            return
         endif
         
         call this%rs%ag()
         
         this%gcells = 1
         this%state = 0
         
         call this%err%werrfl2(class//sname//' ended')
                 
      end subroutine acopyguard_field2d
!
      subroutine pipesend_field2d(this,stag,id) 

         implicit none
         
         class(field2d), intent(inout) :: this
         integer, intent(in) :: stag
         integer, intent(inout) :: id
         character(len=18), save :: sname = 'pipesend_field2d:'

         call this%err%werrfl2(class//sname//' started')
         
         if (this%state == 0) then
            call this%rs%psend(stag,id)
         else
            call this%ks%psend(stag,id)
         endif
         
         call this%err%werrfl2(class//sname//' ended')
                 
      end subroutine pipesend_field2d
!
      subroutine piperecv_field2d(this,rtag) 

         implicit none
         
         class(field2d), intent(inout) :: this
         integer, intent(in) :: rtag
         character(len=18), save :: sname = 'piperecv_field2d:'

         call this%err%werrfl2(class//sname//' started')
         
         if (this%state == 0) then
            call this%rs%precv(rtag)
         else
            call this%ks%precv(rtag)
         endif
         
         call this%err%werrfl2(class//sname//' ended')
                 
      end subroutine piperecv_field2d
!
      subroutine asc(this,value) 

         implicit none
         
         class(field2d), intent(inout) :: this
         real, intent(in) :: value
         character(len=18), save :: sname = 'asc:'

         call this%err%werrfl2(class//sname//' started')
         
         call this%rs%as(value)         
         
         this%state = 0
         this%gcells = 1
         
         call this%err%werrfl2(class//sname//' ended')
                 
      end subroutine asc
!
      subroutine copyto(this,that,lpos,sdim,ddim)
      
         implicit none
         
         class(field2d), intent(inout) :: this
         class(field3d), intent(in) :: that
         integer, intent(in) :: lpos
         integer, intent(in), dimension(:) :: sdim, ddim
! local data         
         character(len=20), save :: sname = 'copyto:'
         class(ufield3d), pointer :: rs3d

         call this%err%werrfl2(class//sname//' started')
         
         rs3d => that%getrs()
         
         call rs3d%cp(this%rs,lpos,sdim,ddim)
                  
         call this%err%werrfl2(class//sname//' ended')
         
      end subroutine copyto
!
      subroutine copyfrom(this,that,lpos,sdim,ddim)
      
         implicit none
         
         class(field2d), intent(inout) :: this
         class(field3d), intent(in) :: that
         integer, intent(in) :: lpos
         integer, intent(in), dimension(:) :: sdim, ddim
! local data         
         character(len=20), save :: sname = 'copyfrom:'
         class(ufield3d), pointer :: rs3d


         call this%err%werrfl2(class//sname//' started')
         
         rs3d => that%getrs()
         
         call rs3d%cb(this%rs,lpos,sdim,ddim)
         
         this%gcells = 1
         this%state = 0
                  
         call this%err%werrfl2(class//sname//' ended')
         
      end subroutine copyfrom
!      
      subroutine asa(this,that) 

         implicit none
         
         class(field2d), intent(inout) :: this,that
         character(len=18), save :: sname = 'asc:'

         call this%err%werrfl2(class//sname//' started')
         
         if (that%state == 0) then
            call this%rs%as(that%rs)         
         else
            call this%ks%as(that%ks)
         endif

         this%state = that%state
         this%gcells = that%gcells
         
         call this%err%werrfl2(class//sname//' ended')
                 
      end subroutine asa
!
      subroutine sum1(this,a1,a2) 

         implicit none
         
         class(field2d), intent(inout) :: this
         class(field2d), intent(in) :: a1, a2
         character(len=18), save :: sname = 'sum:'

         call this%err%werrfl2(class//sname//' started')
         
         if (a1%state /= a2%state) then
            call this%err%equit(class//sname//'states are different')
            return
         else if (a1%state == 0) then
            if (a1%gcells /= a2%gcells) then
               call this%err%equit(class//sname//'gcells are different')
               return
            else            
               call this%rs%add(a1%rs,a2%rs)
               this%state = 0
               this%gcells = a1%gcells
            endif
         else
            call this%ks%add(a1%ks,a2%ks)
            this%state = 1         
         endif
         
         call this%err%werrfl2(class//sname//' ended')
                 
      end subroutine sum1
!
      subroutine sum2(this,a1,a2,dim,dim1,dim2) 

         implicit none
         
         class(field2d), intent(inout) :: this
         class(field2d), intent(in) :: a1, a2
         integer, dimension(:), intent(in) :: dim, dim1, dim2
         character(len=18), save :: sname = 'sum:'

         call this%err%werrfl2(class//sname//' started')
         
         if (a1%state /= a2%state) then
            call this%err%equit(class//sname//'states are different')
            return
         else if (a1%state == 0) then
            if (a1%gcells /= a2%gcells) then
               call this%err%equit(class//sname//'gcells are different')
               return
            else            
               call this%rs%add(a1%rs,a2%rs,dim,dim1,dim2)
               this%state = 0
               this%gcells = a1%gcells
            endif
         else
            call this%ks%add(a1%ks,a2%ks,dim,dim1,dim2)
            this%state = 1         
         endif
         
         call this%err%werrfl2(class//sname//' ended')
                 
      end subroutine sum2
!
      subroutine minus1(this,a1,a2) 

         implicit none
         
         class(field2d), intent(inout) :: this
         class(field2d), intent(in) :: a1, a2
         character(len=18), save :: sname = 'minus:'

         call this%err%werrfl2(class//sname//' started')
         
         if (a1%state /= a2%state) then
            call this%err%equit(class//sname//'states are different')
            return
         else if (a1%state == 0) then
            if (a1%gcells /= a2%gcells) then
               call this%err%equit(class//sname//'gcells are different')
               return
            else            
               call this%rs%sub(a1%rs,a2%rs)
               this%state = 0
               this%gcells = a1%gcells
            endif
         else
            call this%ks%sub(a1%ks,a2%ks)
            this%state = 1         
         endif
         
         call this%err%werrfl2(class//sname//' ended')
                 
      end subroutine minus1
!
      subroutine minus2(this,a1,a2,dim,dim1,dim2) 

         implicit none
         
         class(field2d), intent(inout) :: this
         class(field2d), intent(in) :: a1, a2
         integer, dimension(:), intent(in) :: dim, dim1, dim2         
         character(len=18), save :: sname = 'minus:'

         call this%err%werrfl2(class//sname//' started')
         
         if (a1%state /= a2%state) then
            call this%err%equit(class//sname//'states are different')
            return
         else if (a1%state == 0) then
            if (a1%gcells /= a2%gcells) then
               call this%err%equit(class//sname//'gcells are different')
               return
            else            
               call this%rs%sub(a1%rs,a2%rs,dim,dim1,dim2)
               this%state = 0
               this%gcells = a1%gcells
            endif               
         else
            call this%ks%sub(a1%ks,a2%ks,dim,dim1,dim2)
            this%state = 1         
         endif
         
         call this%err%werrfl2(class//sname//' ended')
                 
      end subroutine minus2
!
      subroutine multiply1(this,a1,value) 

         implicit none
         
         class(field2d), intent(inout) :: this
         class(field2d), intent(in) :: a1
         real, intent(in) :: value
         character(len=18), save :: sname = 'multiply1:'

         call this%err%werrfl2(class//sname//' started')
         
         if (a1%state == 0) then
            call this%rs%mult(a1%rs,value)
            this%state = 0
            this%gcells = a1%gcells            
         else
            call this%ks%mult(a1%ks,value)
            this%state = 1         
         endif
         
         call this%err%werrfl2(class//sname//' ended')
                 
      end subroutine multiply1
!
      subroutine multiply2(this,a1,dim,dim1,value) 

         implicit none
         
         class(field2d), intent(inout) :: this
         class(field2d), intent(in) :: a1
         integer, dimension(:), intent(in) :: dim,dim1
         real, dimension(:), intent(in) :: value
         character(len=18), save :: sname = 'multiply2:'

         call this%err%werrfl2(class//sname//' started')
         
         if (a1%state == 0) then
            call this%rs%mult(a1%rs,dim,dim1,value)
            this%state = 0
            this%gcells = a1%gcells                        
         else
            call this%ks%mult(a1%ks,dim,dim1,value)
            this%state = 1         
         endif
         
         call this%err%werrfl2(class//sname//' ended')
                 
      end subroutine multiply2
!
      subroutine writehdf5_field2d(this,file,dim) 

         implicit none
         
         class(field2d), intent(inout) :: this
         class(hdf5file), intent(in) :: file
         integer, intent(in) :: dim
         character(len=20), save :: sname = 'writehdf5_field2d:'

         call this%err%werrfl2(class//sname//' started')
         
         call this%rs%wr(file,dim)
         
         call this%err%werrfl2(class//sname//' ended')
                 
      end subroutine writehdf5_field2d
!
      function getstate(this)

         implicit none

         class(field2d), intent(in) :: this
         integer :: getstate

         getstate = this%state

      end function getstate      
!      
      function getgcells(this)

         implicit none

         class(field2d), intent(in) :: this
         integer :: getgcells

         getgcells = this%gcells

      end function getgcells
!      
      function getrs(this)

         implicit none

         class(field2d), intent(in) :: this
         class(ufield2d), pointer :: getrs

         getrs => this%rs

      end function getrs
!      
      function getks(this)

         implicit none

         class(field2d), intent(in) :: this
         class(ufield2d), pointer :: getks

         getks => this%ks

      end function getks
!      
      end module field2d_class