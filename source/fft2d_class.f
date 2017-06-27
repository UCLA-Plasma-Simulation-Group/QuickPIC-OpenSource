! fft2d class for QuickPIC Open Source 1.0
! update: 04/18/2016

      module fft2d_class

      use perrors_class
      use parallel_pipe_class
      use spect2d_class
      use ufield2d_class
      use fft2d_lib
         
      implicit none

      private

      public :: fft2d, get_fft2table

      type fft2d

         private
!
! ind = exponent which determines length in each direction
! nrc = (0,1) = table for real to complex (0) or complex to complex (1)
! mixup = array of bit reversed addresses
! sct = sine/cosine table         
         class(spect2d), pointer, public :: sp => null()
         class(perrors), pointer, public :: err => null()
         class(parallel_pipe), pointer, public :: p => null()
         integer, dimension(2) :: ind
         integer :: nrc
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
         
         contains
         
         generic :: new => init_fft2d
         generic :: del => end_fft2d
         generic :: fsst => iwpfsst2r         
         generic :: fcct => iwpfcct2r       
         generic :: fs3t => iwpfs3t2r
         generic :: divf => ipdivfd2
         generic :: gradf => ipgradfd2
         generic :: curlf => ipcurlfd2
         procedure, private :: init_fft2d
         procedure, private :: end_fft2d
         procedure, private :: iwpfsst2r, iwpfcct2r, iwpfs3t2r 
         procedure, private :: ipdivfd2, ipgradfd2, ipcurlfd2
      end type 
!
      type fft2d_link
         type (fft2d_link), pointer :: next => null()
         type (fft2d), pointer :: table => null()
         integer :: refcount
      end type fft2d_link
!      
      character(len=10), save :: class = 'fft2d:'
      character(len=128), save :: erstr
! link list for fft tables
      integer, save :: numtables = 0
      type (fft2d_link), target, save :: table_list
      
      contains
!
      subroutine init_fft2d(this,pp,perr,psp,indx,indy,nrc)
      
         implicit none
         
         class(fft2d), intent(inout) :: this
         class(spect2d), intent(in), pointer :: psp
         class(perrors), intent(in), pointer :: perr
         class(parallel_pipe), intent(in), pointer :: pp
         integer, intent(in) :: indx,indy,nrc
! local data
         character(len=18), save :: sname = 'init_fft2d:'
         integer :: nx, ny, n1, n2
         
         this%sp => psp
         this%err => perr
         this%p => pp
         call this%err%werrfl2(class//sname//' started')

         if ((indx < 1) .or. (indy < 1)) then
            write (erstr,*) 'invalid indx or indy=', indx, indy
            call this%err%equit(class//sname//erstr)
            return
         endif
         
         this%ind = (/indx,indy/)
         this%nrc = nrc
         nx = 2**indx; ny = 2**indy
         select case (this%sp%getpsolver())
         case (1)
            n1 = max(nx/2,ny); n2 = max(nx,ny)
            allocate(this%mixup(n1),this%sct(n2))
            call WPFST2RINIT(this%mixup,this%sct,indx,indy,n1,n2)
         case default
            n1 = max(nx/2,ny); n2 = max(nx,ny)
            allocate(this%mixup(n1),this%sct(n2))
            call WPFST2RINIT(this%mixup,this%sct,indx,indy,n1,n2)
         end select

         call this%err%werrfl2(class//sname//' ended')

      end subroutine init_fft2d
!
      subroutine end_fft2d(this)
          
         implicit none
         
         class(fft2d), intent(inout) :: this
! local data
         character(len=18), save :: sname = 'end_fft2d:'

         call this%err%werrfl2(class//sname//' started')
         if (associated(this%mixup)) deallocate(this%mixup)
         if (associated(this%sct)) deallocate(this%sct)
         call this%err%werrfl2(class//sname//' ended')
                  
      end subroutine end_fft2d
!      
      subroutine iwpfsst2r(this,rspace,krspace,isign)
! this subroutine performs 2d real sine-cosine-sine transform for 1d partition
! rspace, krspace = ufield2d of input and output data
! isign = sign for transform (-1 = real to fourier, 1 = fourier to real)
         implicit none
         class(fft2d), intent(in) :: this
         class(ufield2d), intent(inout) :: rspace, krspace
         integer, intent(in) :: isign
! local data
! f = real space source or destination for transform
! g = real fourier space source or destination for transform
! inorder = interpolation order, determines starting point for transform
! kstrt = starting data block number, a global variable
         real, dimension(:,:,:), pointer :: f => null()
         real, dimension(:,:,:), pointer :: g => null()

         integer, dimension(:), pointer :: mixup => null()
         complex, dimension(:), pointer :: sctd => null()
         integer :: indx, indy
         integer :: ntpose = 1, nxvh, nyv, kyp, kxp2, kypd, kxp2d
         integer :: jblok, kblok, nxhyd, nxyd, order
         real :: ttp
         real, dimension(:,:,:), allocatable :: bs, br
         character(len=11), save :: sname = 'iwpfsst2r:'
! check for errors
         call this%err%werrfl2(class//sname//' started')
         if ((rspace%getlayout() /= 0) .or. (krspace%getlayout() /= 1)) &
     & then
            erstr = ' invalid layout'
            call this%err%equit(class//sname//erstr)
            return
         endif
         if ((rspace%getnd1() /= krspace%getnd2()).or.&
         &(rspace%getnd2() /= krspace%getnd1())) then
            erstr = ' non-conforming array'
            call this%err%equit(class//sname//erstr)
            return
         endif
! unpack arguments
         indx = this%ind(1); indy = this%ind(2)
         kyp = rspace%getnd2p(); kxp2 = krspace%getnd2p()
         kblok = 1; jblok = 1
         allocate(bs(rspace%getdim(),kxp2+1,kyp+1))
         allocate(br(krspace%getdim(),kxp2+1,kyp+1))
         f => rspace%getrf()
         g => krspace%getrf()
         nxvh = size(f,2)/2; nyv = size(g,2)
         kypd = size(f,3); kxp2d = size(g,3)
         nxhyd = size(this%mixup); nxyd = size(this%sct)
         mixup => this%mixup; sctd => this%sct
! choose the proper function
         order = this%sp%getinorder()
         select case (order)
         case (1)
            select case (rspace%getdim())
            case (1)
               call WPPFSST2RM(f(1,:,:),g,bs,br,isign,ntpose,mixup,sctd,ttp,&
               &indx,indy,this%p%getlkstrt(),this%p%getlnvp(),nxvh,nyv,kxp2,&
               &kyp,kypd,kxp2d,nxhyd,nxyd)
!               call WPFSST2R(f(1,1,1),g,bs,br,isign,ntpose,mixup,sctd,ttp,&
!               &indx,indy,this%p%getlkstrt(),nxvh,nyv,kxp2,kyp,kypd,kxp2d,&
!               &jblok,kblok,nxhyd,nxyd)
            case (2)
               call WPPFCST2RM2(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,&
               &indx,indy,this%p%getlkstrt(),this%p%getlnvp(),nxvh,nyv,kxp2,&
               &kyp,kypd,kxp2d,nxhyd,nxyd)
!               call WPFCST2R2(f(1,1,1),g,bs,br,isign,ntpose,mixup,sctd,ttp,&
!               &indx,indy,this%p%getlkstrt(),nxvh,nyv,kxp2,kyp,kypd,kxp2d,&
!               &jblok,kblok,nxhyd,nxyd)
            case (3)
               call WPPFCST2RM3(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,&
               &indx,indy,this%p%getlkstrt(),this%p%getlnvp(),nxvh,nyv,kxp2,&
               &kyp,kypd,kxp2d,nxhyd,nxyd)
!               call WPFCST2R3(f(1,1,1),g,bs,br,isign,ntpose,mixup,sctd,ttp,&
!               &indx,indy,this%p%getlkstrt(),nxvh,nyv,kxp2,kyp,kypd,kxp2d,&
!               &jblok,kblok,nxhyd,nxyd)
            end select
         case default
            select case (rspace%getdim())
            case (1)
               call WPPFSST2RM(f(1,:,:),g,bs,br,isign,ntpose,mixup,sctd,ttp,&
               &indx,indy,this%p%getlkstrt(),this%p%getlnvp(),nxvh,nyv,kxp2,&
               &kyp,kypd,kxp2d,nxhyd,nxyd)
!               call WPFSST2R(f(1,1,1),g,bs,br,isign,ntpose,mixup,sctd,ttp,&
!               &indx,indy,this%p%getlkstrt(),nxvh,nyv,kxp2,kyp,kypd,kxp2d,&
!               &jblok,kblok,nxhyd,nxyd)
            case (2)
               call WPPFCST2RM2(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,&
               &indx,indy,this%p%getlkstrt(),this%p%getlnvp(),nxvh,nyv,kxp2,&
               &kyp,kypd,kxp2d,nxhyd,nxyd)
!               call WPFCST2R2(f(1,1,1),g,bs,br,isign,ntpose,mixup,sctd,ttp,&
!               &indx,indy,this%p%getlkstrt(),nxvh,nyv,kxp2,kyp,kypd,kxp2d,&
!               &jblok,kblok,nxhyd,nxyd)
            case (3)
               call WPPFCST2RM3(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,&
               &indx,indy,this%p%getlkstrt(),this%p%getlnvp(),nxvh,nyv,kxp2,&
               &kyp,kypd,kxp2d,nxhyd,nxyd)
!               call WPFCST2R3(f(1,1,1),g,bs,br,isign,ntpose,mixup,sctd,ttp,&
!               &indx,indy,this%p%getlkstrt(),nxvh,nyv,kxp2,kyp,kypd,kxp2d,&
!               &jblok,kblok,nxhyd,nxyd)
            end select
         end select
         deallocate(bs,br)
         call this%err%werrfl2(class//sname//' ended')

      end subroutine iwpfsst2r
!
      subroutine iwpfcct2r(this,rspace,krspace,isign)
! this subroutine performs 2d real cosine-sine-cosine transform for 1d partition
! rspace, krspace = ufield2d of input and output data
! isign = sign for transform (-1 = real to fourier, 1 = fourier to real)
         implicit none
         class(fft2d), intent(in) :: this
         class(ufield2d), intent(inout) :: rspace, krspace
         integer, intent(in) :: isign
! local data
! f = real space source or destination for transform
! g = real fourier space source or destination for transform
! inorder = interpolation order, determines starting point for transform
! kstrt = starting data block number, a global variable
         real, dimension(:,:,:), pointer :: f => null()
         real, dimension(:,:,:), pointer :: g => null()

         integer, dimension(:), pointer :: mixup => null()
         complex, dimension(:), pointer :: sctd => null()
         integer :: indx, indy
         integer :: ntpose = 1, nxvh, nyv, kyp, kxp2, kypd, kxp2d
         integer :: jblok, kblok, nxhyd, nxyd, order
         real :: ttp
         real, dimension(:,:,:), allocatable :: bs, br
         character(len=11), save :: sname = 'iwpfcct2r:'
! check for errors
         call this%err%werrfl2(class//sname//' started')
         if ((rspace%getlayout() /= 0) .or. (krspace%getlayout() /= 1)) &
     & then
            erstr = ' invalid layout'
            call this%err%equit(class//sname//erstr)
            return
         endif
         if ((rspace%getnd1() /= krspace%getnd2()).or.&
         &(rspace%getnd2() /= krspace%getnd1())) then
            erstr = ' non-conforming array'
            call this%err%equit(class//sname//erstr)
            return
         endif
! unpack arguments
         indx = this%ind(1); indy = this%ind(2)
         kyp = rspace%getnd2p(); kxp2 = krspace%getnd2p()
         kblok = 1; jblok = 1
         allocate(bs(rspace%getdim(),kxp2+1,kyp+1))
         allocate(br(krspace%getdim(),kxp2+1,kyp+1))
         f => rspace%getrf()
         g => krspace%getrf()
         nxvh = size(f,2)/2; nyv = size(g,2)
         kypd = size(f,3); kxp2d = size(g,3)
         nxhyd = size(this%mixup); nxyd = size(this%sct)
         mixup => this%mixup; sctd => this%sct
! choose the proper function
         order = this%sp%getinorder()
         select case (order)
         case (1)
            select case (rspace%getdim())
            case (1)
               call WPPFCCT2RM(f(1,:,:),g,bs,br,isign,ntpose,mixup,sctd,ttp,&
               &indx,indy,this%p%getlkstrt(),this%p%getlnvp(),nxvh,nyv,kxp2,&
               &kyp,kypd,kxp2d,nxhyd,nxyd)
!               call WPFCCT2R(f(1,1,1),g,bs,br,isign,ntpose,mixup,sctd,ttp,&
!               &indx,indy,this%p%getlkstrt(),nxvh,nyv,kxp2,kyp,kypd,kxp2d,&
!               &jblok,kblok,nxhyd,nxyd)
            case (2)
               call WPPFSCT2RM2(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,&
               &indx,indy,this%p%getlkstrt(),this%p%getlnvp(),nxvh,nyv,kxp2,&
               &kyp,kypd,kxp2d,nxhyd,nxyd)
!               call WPFSCT2R2(f(1,1,1),g,bs,br,isign,ntpose,mixup,sctd,ttp,&
!               &indx,indy,this%p%getlkstrt(),nxvh,nyv,kxp2,kyp,kypd,kxp2d,&
!               &jblok,kblok,nxhyd,nxyd)
            case (3)
               call WPPFSCT2RM3(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,&
               &indx,indy,this%p%getlkstrt(),this%p%getlnvp(),nxvh,nyv,kxp2,&
               &kyp,kypd,kxp2d,nxhyd,nxyd)
!               call WPFSCT2R3(f(1,1,1),g,bs,br,isign,ntpose,mixup,sctd,ttp,&
!               &indx,indy,this%p%getlkstrt(),nxvh,nyv,kxp2,kyp,kypd,kxp2d,&
!               &jblok,kblok,nxhyd,nxyd)
            end select
         case default
            select case (rspace%getdim())
            case (1)
               call WPPFCCT2RM(f(1,:,:),g,bs,br,isign,ntpose,mixup,sctd,ttp,&
               &indx,indy,this%p%getlkstrt(),this%p%getlnvp(),nxvh,nyv,kxp2,&
               &kyp,kypd,kxp2d,nxhyd,nxyd)
!               call WPFCCT2R(f(1,1,1),g,bs,br,isign,ntpose,mixup,sctd,ttp,&
!               &indx,indy,this%p%getlkstrt(),nxvh,nyv,kxp2,kyp,kypd,kxp2d,&
!               &jblok,kblok,nxhyd,nxyd)
            case (2)
               call WPPFSCT2RM2(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,&
               &indx,indy,this%p%getlkstrt(),this%p%getlnvp(),nxvh,nyv,kxp2,&
               &kyp,kypd,kxp2d,nxhyd,nxyd)
!               call WPFSCT2R2(f(1,1,1),g,bs,br,isign,ntpose,mixup,sctd,ttp,&
!               &indx,indy,this%p%getlkstrt(),nxvh,nyv,kxp2,kyp,kypd,kxp2d,&
!               &jblok,kblok,nxhyd,nxyd)
            case (3)
               call WPPFSCT2RM3(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,&
               &indx,indy,this%p%getlkstrt(),this%p%getlnvp(),nxvh,nyv,kxp2,&
               &kyp,kypd,kxp2d,nxhyd,nxyd)
!               call WPFSCT2R3(f(1,1,1),g,bs,br,isign,ntpose,mixup,sctd,ttp,&
!               &indx,indy,this%p%getlkstrt(),nxvh,nyv,kxp2,kyp,kypd,kxp2d,&
!               &jblok,kblok,nxhyd,nxyd)
            end select
         end select
         deallocate(bs,br)
         call this%err%werrfl2(class//sname//' ended')

      end subroutine iwpfcct2r
!         
      subroutine iwpfs3t2r(this,rspace,krspace,isign)
! this subroutine performs 2d real sine-cosine-sine transform for 1d partition
! rspace, krspace = ufield2d of input and output data
! isign = sign for transform (-1 = real to fourier, 1 = fourier to real)
         implicit none
         class(fft2d), intent(in) :: this
         class(ufield2d), intent(inout) :: rspace, krspace
         integer, intent(in) :: isign
! local data
! f = real space source or destination for transform
! g = real fourier space source or destination for transform
! inorder = interpolation order, determines starting point for transform
! kstrt = starting data block number, a global variable
         real, dimension(:,:,:), pointer :: f => null()
         real, dimension(:,:,:), pointer :: g => null()

         integer, dimension(:), pointer :: mixup => null()
         complex, dimension(:), pointer :: sctd => null()
         integer :: indx, indy
         integer :: ntpose = 1, nxvh, nyv, kyp, kxp2, kypd, kxp2d
         integer :: jblok, kblok, nxhyd, nxyd, order
         real :: ttp
         real, dimension(:,:,:), allocatable :: bs, br
         character(len=11), save :: sname = 'iwpfs3t2r:'
! check for errors
         call this%err%werrfl2(class//sname//' started')
         if ((rspace%getlayout() /= 0) .or. (krspace%getlayout() /= 1)) &
     & then
            erstr = ' invalid layout'
            call this%err%equit(class//sname//erstr)
            return
         endif
         if ((rspace%getnd1() /= krspace%getnd2()).or.&
         &(rspace%getnd2() /= krspace%getnd1())) then
            erstr = ' non-conforming array'
            call this%err%equit(class//sname//erstr)
            return
         endif
! unpack arguments
         indx = this%ind(1); indy = this%ind(2)
         kyp = rspace%getnd2p(); kxp2 = krspace%getnd2p()
         allocate(bs(rspace%getdim(),kxp2+1,kyp+1))
         allocate(br(rspace%getdim(),kxp2+1,kyp+1))
         kblok = 1; jblok = 1
         f => rspace%getrf()
         g => krspace%getrf()
         nxvh = size(f,2)/2; nyv = size(g,2)
         kypd = size(f,3); kxp2d = size(g,3)
         nxhyd = size(this%mixup); nxyd = size(this%sct)
         mixup => this%mixup; sctd => this%sct
! choose the proper function
         order = this%sp%getinorder()
         
         select case (order)
         case (1)
            select case (rspace%getdim())
            case (3)
               call WPPFSST2RM23(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,&
               &indx,indy,this%p%getlkstrt(),this%p%getlnvp(),nxvh,nyv,kxp2,&
               &kyp,kypd,kxp2d,nxhyd,nxyd)
!               call WPFS3T2R3(f(1,1,1),g,bs,br,isign,ntpose,mixup,sctd,ttp,&
!               &indx,indy,this%p%getlkstrt(),nxvh,nyv,kxp2,kyp,kypd,kxp2d,&
!               &jblok,kblok,nxhyd,nxyd)
            end select
         case default
            select case (rspace%getdim())
            case (3)
               call WPPFSST2RM23(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,&
               &indx,indy,this%p%getlkstrt(),this%p%getlnvp(),nxvh,nyv,kxp2,&
               &kyp,kypd,kxp2d,nxhyd,nxyd)
!               call WPFS3T2R3(f(1,1,1),g,bs,br,isign,ntpose,mixup,sctd,ttp,&
!               &indx,indy,this%p%getlkstrt(),nxvh,nyv,kxp2,kyp,kypd,kxp2d,&
!               &jblok,kblok,nxhyd,nxyd)
            end select
         end select
         deallocate(bs,br)
         call this%err%werrfl2(class//sname//' ended')

      end subroutine iwpfs3t2r
!
      subroutine ipdivfd2(this,krspace,kdspace)
! this subroutine calculates the divergence in fourier space
! with dirichlet (zero potential) boundary conditions
! this = fft2d descriptor
! krspace = source ufield2d
! kdspace = destiny ufield2d
         implicit none
         class(fft2d), intent(in) :: this
         class(ufield2d), intent(inout) :: krspace, kdspace
! local data
         real, dimension(:,:,:), pointer :: f => null()
         real, dimension(:,:,:), pointer :: df => null()
         integer :: nx, ny, ndim, nyv, kxp2
         character(len=10), save :: sname = 'ipdivfd2:'

         call this%err%werrfl2(class//sname//' started')
! unpack arguments
         nx = 2**this%ind(1); ny = 2**this%ind(2)
         f => krspace%getrf()
         df => kdspace%getrf()
         ndim = size(f,1)
         nyv = size(f,2); kxp2 = size(f,3) - 1;
! call the operator
         call MPPDIVFD2(f,df(1,:,:),nx,ny,this%p%getlkstrt(),ndim,nyv,kxp2)
         call this%err%werrfl2(class//sname//' ended')

      end subroutine ipdivfd2
!      
      subroutine ipgradfd2(this,krspace,kdspace)
! this subroutine calculates the gradient in fourier space
! with dirichlet (zero potential) boundary conditions
! this = fft2d descriptor
! krspace = source ufield2d
! kdspace = destiny ufield2d
         implicit none
         class(fft2d), intent(in) :: this
         class(ufield2d), intent(inout) :: krspace, kdspace
! local data
         integer :: nx, ny, ndim, nyv, kxp2, j2blok
         real, dimension(:,:,:), pointer :: df => null()
         real, dimension(:,:,:), pointer :: f => null()
         character(len=11), save :: sname = 'ipgradfd2:'

         call this%err%werrfl2(class//sname//' started')
! unpack arguments
         nx = 2**this%ind(1); ny = 2**this%ind(2)
         f => krspace%getrf()
         df => kdspace%getrf()
         ndim = size(df,1)
         nyv = size(df,2); kxp2 = size(df,3) - 1;
! call the operator
         call MPPGRADFD2(f(1,:,:),df,nx,ny,this%p%getlkstrt(),ndim,nyv,kxp2)
         call this%err%werrfl2(class//sname//' ended')

      end subroutine ipgradfd2
!
      subroutine ipcurlfd2(this,krspace,kdspace)
! this subroutine calculates the curl in fourier space
! with dirichlet (zero potential) boundary conditions
! this = fft2d descriptor
! krspace = ufield2d descriptor of data
! krspace = source ufield2d
! kdspace = destiny ufield2d
         implicit none
         class(fft2d), intent(in) :: this
         class(ufield2d), intent(inout) :: krspace, kdspace
! local data
         integer :: nx, ny, nyv, kxp2, j2blok
         real, dimension(:,:,:), pointer :: pf => null(), pg => null()
         character(len=11), save :: sname = 'ipcurlfd2:'

         call this%err%werrfl2(class//sname//' started')
! unpack arguments
         nx = 2**this%ind(1); ny = 2**this%ind(2)
         pf => krspace%getrf()
         pg => kdspace%getrf()
         nyv = size(pf,2); kxp2 = size(pf,3) - 1; j2blok = 1
! choose the proper function
         select case (size(pf,1))
         case (2)
            call PCURLFD22(pf,pg,nx,ny,this%p%getlkstrt(),nyv,kxp2,j2blok)
         case (3)
            call MPPCURLFD2(pf,pg,nx,ny,this%p%getlkstrt(),nyv,kxp2)         
         end select
         call this%err%werrfl2(class//sname//' ended')

      end subroutine ipcurlfd2
!
      function get_fft2table(pp,perr,psp,indx,indy) result(table)
! this function gets an fft table entry, either creates one or points
! to one if the required table already exists for real to complex ffts
         implicit none
         class(spect2d), intent(in), pointer :: psp
         class(perrors), intent(in), pointer :: perr
         class(parallel_pipe), intent(in), pointer :: pp         
         integer, intent(in) :: indx, indy
         type (fft2d), pointer :: table
! local data
         type (fft2d_link), pointer :: link => null()
         type (fft2d), pointer :: ltable => null()
         integer :: nrc = 0, ierr = 0
         character(len=18), save :: sname = 'get_fft2table:'

         call perr%werrfl2(class//sname//' started')

         nullify(table)         
         select case (psp%getpsolver())
         case (1)
            nrc = 0
         case default
            nrc = 0
         end select
         
         if (numtables == 0) then
            nullify(table_list%next,table_list%table)         
            table_list%refcount = 0
         endif
         link => table_list
         table => link%table
! search link list of table to see if required table already exists
         do while (associated(table))
! found it
            if ((indx==table%ind(1)) .and. (indy==table%ind(2)) .and. (t&
            &able%nrc==nrc).and.(psp%getpsolver()==table%sp%getpsolver())) then
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
         call table%new(pp,perr,psp,indx,indy,nrc)
         link%refcount = 1
         numtables = numtables + 1

      end function get_fft2table
!
         
      end module fft2d_class