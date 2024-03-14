! fdist2d class for QuickPIC Open Source 1.0
! update: 04/18/2016

      module fdist2d_class

      use perrors_class
      use parallel_pipe_class
      use spect2d_class
      use ufield2d_class
      use input_class
      use m_fparser
      implicit none

      private

      public :: fdist2d, fdist2d_000, fdist2d_010, fdist2d_011, fdist2d_012, fdist2d_013

      type, abstract :: fdist2d

         private

         class(spect2d), pointer, public :: sp => null()
         class(perrors), pointer, public :: err => null()
         class(parallel_pipe), pointer, public :: p => null()
!
! ndprof = profile type 
         integer :: npf, npmax
                          
         contains
         generic :: new => init_fdist2d         
         generic :: del => end_fdist2d
         generic :: dist => dist2d
         procedure(ab_init_fdist2d), deferred, private :: init_fdist2d
         procedure, private :: end_fdist2d
         procedure(ab_dist2d), deferred, private :: dist2d
         procedure :: getnpf, getnpmax
                  
      end type fdist2d

      abstract interface
!
      subroutine ab_dist2d(this,part2d,npp,fd,s)
         import fdist2d
         import ufield2d
         implicit none
         class(fdist2d), intent(inout) :: this
         real, dimension(:,:), pointer, intent(inout) :: part2d
         integer, intent(inout) :: npp
         class(ufield2d), intent(in), pointer :: fd         
         real, intent(in) :: s
      end subroutine ab_dist2d
!
      subroutine ab_init_fdist2d(this,input,i)
         import fdist2d
         import input_json
         implicit none
         class(fdist2d), intent(inout) :: this
         type(input_json), intent(inout), pointer :: input
         integer, intent(in) :: i
      end subroutine ab_init_fdist2d
!
      end interface
!
      type, extends(fdist2d) :: fdist2d_000
! Transeversely uniform profile with uniform or piecewise longitudinal profile
         private
! xppc, yppc = particle per cell in x and y directions
         integer :: xppc, yppc
         real :: qm, den
         character(len=:), allocatable :: long_prof
         real, dimension(:), allocatable :: s, fs
                          
         contains
         procedure, private :: init_fdist2d => init_fdist2d_000
         procedure, private :: dist2d => dist2d_000
                  
      end type fdist2d_000
!
      type, extends(fdist2d) :: fdist2d_010
! sharp edge hollow channel with different radius in x and y
         private
! xppc, yppc = particle per cell in x and y directions
         integer :: xppc, yppc
         real :: qm, den
         real :: cx, cy
         real :: irx, iry, orx, ory
         character(len=:), allocatable :: long_prof
         real, dimension(:), allocatable :: s, fs
                          
         contains
         procedure, private :: init_fdist2d => init_fdist2d_010
         procedure, private :: dist2d => dist2d_010
                  
      end type fdist2d_010
!
      type, extends(fdist2d) :: fdist2d_011
!  Hollow channel with Gaussian profile
         private
! xppc, yppc = particle per cell in x and y directions
         integer :: xppc, yppc
         real :: qm, den
         real :: cx, cy
         real :: r0,sigr
         character(len=:), allocatable :: long_prof
         real, dimension(:), allocatable :: s, fs
                          
         contains
         procedure, private :: init_fdist2d => init_fdist2d_011
         procedure, private :: dist2d => dist2d_011
                  
      end type fdist2d_011
!
      type, extends(fdist2d) :: fdist2d_012
! hollow channel with f(r) profile
         private
! xppc, yppc = particle per cell in x and y directions
         integer :: xppc, yppc
         real :: qm, den
         real :: cx, cy
         real, dimension(:), allocatable :: r, fr
         character(len=:), allocatable :: long_prof
         real, dimension(:), allocatable :: s, fs
                          
         contains
         procedure, private :: init_fdist2d => init_fdist2d_012
         procedure, private :: dist2d => dist2d_012
                  
      end type fdist2d_012

      type, extends(fdist2d) :: fdist2d_013
! Transeversely uniform profile with uniform or piecewise longitudinal profile
         private
! xppc, yppc = particle per cell in x and y directions
         integer :: xppc, yppc
         real :: qm, den
         real :: cx, cy, dx, dy
         ! analytic density math function
         type(t_fparser), pointer :: math_func => null()
                          
         contains
         procedure, private :: init_fdist2d => init_fdist2d_013
         procedure, private :: dist2d => dist2d_013
                  
      end type fdist2d_013
!
      character(len=10), save :: class = 'fdist2d:'
      character(len=128), save :: erstr
      
      contains
!
      function getnpf(this)

         implicit none

         class(fdist2d), intent(in) :: this
         integer :: getnpf
         
         getnpf = this%npf
      
      end function getnpf
!      
      function getnpmax(this)

         implicit none

         class(fdist2d), intent(in) :: this
         integer :: getnpmax
         
         getnpmax = this%npmax
      
      end function getnpmax
!
      subroutine end_fdist2d(this)
          
         implicit none
         
         class(fdist2d), intent(inout) :: this
         character(len=18), save :: sname = 'end_fdist2d:'

         call this%err%werrfl2(class//sname//' started')
         call this%err%werrfl2(class//sname//' ended')
                  
      end subroutine end_fdist2d
!      
      subroutine init_fdist2d_000(this,input,i)
      
         implicit none
         
         class(fdist2d_000), intent(inout) :: this
         type(input_json), intent(inout), pointer :: input
         integer, intent(in) :: i
! local data
         integer :: npf,xppc,yppc,npmax,indx,indy
         real :: qm, den
         character(len=20) :: sn,s1
         character(len=18), save :: sname = 'init_fdist2d_000:'
         
         this%sp => input%sp
         this%err => input%err
         this%p => input%pp

         call this%err%werrfl2(class//sname//' started')
         write (sn,'(I3.3)') i
         s1 = 'species('//trim(sn)//')'
         call input%get('simulation.indx',indx)
         call input%get('simulation.indy',indy)
         call input%get(trim(s1)//'.profile',npf)
         call input%get(trim(s1)//'.ppc(1)',xppc)
         call input%get(trim(s1)//'.ppc(2)',yppc)
         call input%get(trim(s1)//'.q',qm)
         call input%get(trim(s1)//'.density',den)
         npmax = xppc*yppc*(2**indx)*(2**indy)/this%p%getlnvp()*4
         call input%get(trim(s1)//'.longitudinal_profile',this%long_prof)
         if (trim(this%long_prof) == 'piecewise') then
            call input%get(trim(s1)//'.piecewise_density',this%fs)
            call input%get(trim(s1)//'.piecewise_s',this%s)
         end if
         this%npf = npf
         this%xppc = xppc
         this%yppc = yppc
         this%qm = qm
         this%den = den
         this%npmax = npmax
         call this%err%werrfl2(class//sname//' ended')

      end subroutine init_fdist2d_000
!
      subroutine dist2d_000(this,part2d,npp,fd,s)
         implicit none
         class(fdist2d_000), intent(inout) :: this
         real, dimension(:,:), pointer, intent(inout) :: part2d
         integer, intent(inout) :: npp
         class(ufield2d), intent(in), pointer :: fd
         real, intent(in) :: s 
! local data
         character(len=18), save :: sname = 'dist2d_000:'
         real, dimension(:,:), pointer :: pt => null()
         integer :: nps, nx, ny, noff, xppc, yppc, i, j
         integer :: ix, iy
         real :: qm, den_temp
         integer :: prof_l

         call this%err%werrfl2(class//sname//' started')
         
         nx = fd%getnd1p(); ny = fd%getnd2p(); noff = fd%getnoff()
         xppc = this%xppc; yppc = this%yppc
         den_temp = 1.0
         if (trim(this%long_prof) == 'piecewise') then
            prof_l = size(this%fs)
            if (s<this%s(1) .or. s>this%s(prof_l)) then
               write (erstr,*) 'The s is out of the bound!'
               call this%err%equit(class//sname//erstr)
               return
            end if
            do i = 2, prof_l
               if (this%s(i) < this%s(i-1)) then
                  write (erstr,*) 's is not monotonically increasing!'
                  call this%err%equit(class//sname//erstr)
                  return
               end if
               if (s<=this%s(i)) then
                  den_temp = this%fs(i-1) + (this%fs(i)-this%fs(i-1))/&
                  &(this%s(i)-this%s(i-1))*(s-this%s(i-1))
                  exit
               end if
            end do
         end if
         qm = den_temp*this%den*this%qm/abs(this%qm)/real(xppc)/real(yppc)
         nps = 1
         pt => part2d
! initialize the particle positions
         if (noff < ny) then
         do i=2, nx-1
            do j=2, ny
               do ix = 0, xppc-1
                  do iy=0, yppc-1
                     pt(1,nps) = (ix + 0.5)/xppc + i - 1
                     pt(2,nps) = (iy + 0.5)/yppc + j - 1 + noff
                     pt(3,nps) = 0.0
                     pt(4,nps) = 0.0
                     pt(5,nps) = 0.0
                     pt(6,nps) = 1.0
                     pt(7,nps) = 1.0
                     pt(8,nps) = qm
                     nps = nps + 1
                  enddo
               enddo
            enddo
         enddo
         else if (noff > (nx-ny-1)) then       
         do i=2, nx-1
            do j=1, ny-1
               do ix = 0, xppc-1
                  do iy=0, yppc-1
                     pt(1,nps) = (ix + 0.5)/xppc + i - 1
                     pt(2,nps) = (iy + 0.5)/yppc + j - 1 + noff
                     pt(3,nps) = 0.0
                     pt(4,nps) = 0.0
                     pt(5,nps) = 0.0
                     pt(6,nps) = 1.0
                     pt(7,nps) = 1.0
                     pt(8,nps) = qm
                     nps = nps + 1
                  enddo
               enddo
            enddo
         enddo
         else
         do i=2, nx-1
            do j=1, ny
               do ix = 0, xppc-1
                  do iy=0, yppc-1
                     pt(1,nps) = (ix + 0.5)/xppc + i - 1
                     pt(2,nps) = (iy + 0.5)/yppc + j - 1 + noff
                     pt(3,nps) = 0.0
                     pt(4,nps) = 0.0
                     pt(5,nps) = 0.0
                     pt(6,nps) = 1.0
                     pt(7,nps) = 1.0
                     pt(8,nps) = qm
                     nps = nps + 1
                  enddo
               enddo
            enddo
         enddo
         endif
         
         npp = nps - 1
         
         call this%err%werrfl2(class//sname//' ended')

      end subroutine dist2d_000
!
      subroutine init_fdist2d_010(this,input,i)
      
         implicit none
         
         class(fdist2d_010), intent(inout) :: this
         type(input_json), intent(inout), pointer :: input
         integer, intent(in) :: i
! local data
         integer :: npf,xppc,yppc,npmax,indx,indy
         real :: qm, den
         real :: cx, cy
         real :: irx, iry, orx, ory
         real :: min, max
         real :: alx, aly, dx, dy
         character(len=20) :: sn,s1
         character(len=18), save :: sname = 'init_fdist2d_010:'
         
         this%sp => input%sp
         this%err => input%err
         this%p => input%pp

         call this%err%werrfl2(class//sname//' started')
         write (sn,'(I3.3)') i
         s1 = 'species('//trim(sn)//')'
         call input%get('simulation.indx',indx)
         call input%get('simulation.indy',indy)
         call input%get('simulation.box.x(1)',min)
         call input%get('simulation.box.x(2)',max)
         call input%get(trim(s1)//'.center(1)',cx)
         cx = cx - min
         alx = (max-min) 
         dx=alx/real(2**indx)
         call input%get('simulation.box.y(1)',min)
         call input%get('simulation.box.y(2)',max)
         call input%get(trim(s1)//'.center(2)',cy)
         cy = cy - min
         aly = (max-min) 
         dy=aly/real(2**indy)
         call input%get(trim(s1)//'.profile',npf)
         call input%get(trim(s1)//'.ppc(1)',xppc)
         call input%get(trim(s1)//'.ppc(2)',yppc)
         call input%get(trim(s1)//'.q',qm)
         call input%get(trim(s1)//'.density',den)
         call input%get(trim(s1)//'.inner_radius(1)',irx)
         call input%get(trim(s1)//'.inner_radius(2)',iry)
         call input%get(trim(s1)//'.outer_radius(1)',orx)
         call input%get(trim(s1)//'.outer_radius(2)',ory)
         npmax = xppc*yppc*(2**indx)*(2**indy)/this%p%getlnvp()*4
         call input%get(trim(s1)//'.longitudinal_profile',this%long_prof)
         if (trim(this%long_prof) == 'piecewise') then
            call input%get(trim(s1)//'.piecewise_density',this%fs)
            call input%get(trim(s1)//'.piecewise_s',this%s)
         end if

         this%npf = npf
         this%xppc = xppc
         this%yppc = yppc
         this%qm = qm
         this%den = den
         this%npmax = npmax
         this%cx = cx/dx
         this%cy = cy/dy
         this%irx = irx/dx
         this%iry = iry/dy
         this%orx = orx/dx
         this%ory = ory/dy

         call this%err%werrfl2(class//sname//' ended')

      end subroutine init_fdist2d_010
!
      subroutine dist2d_010(this,part2d,npp,fd,s)
         implicit none
         class(fdist2d_010), intent(inout) :: this
         real, dimension(:,:), pointer, intent(inout) :: part2d
         integer, intent(inout) :: npp
         class(ufield2d), intent(in), pointer :: fd
         real, intent(in) :: s 
! local data
         character(len=18), save :: sname = 'dist2d_010:'
         real, dimension(:,:), pointer :: pt => null()
         integer :: nps, nx, ny, noff, xppc, yppc, i, j
         integer :: ix, iy
         real :: qm, x, y
         real :: cx, cy
         real :: irx, iry, orx, ory         
         real :: iirx2, iiry2, iorx2, iory2   
         real :: den_temp
         integer :: prof_l

         call this%err%werrfl2(class//sname//' started')
         
         nx = fd%getnd1p(); ny = fd%getnd2p(); noff = fd%getnoff()
         xppc = this%xppc; yppc = this%yppc
         den_temp = 1.0
         if (trim(this%long_prof) == 'piecewise') then
            prof_l = size(this%fs)
            if (s<this%s(1) .or. s>this%s(prof_l)) then
               write (erstr,*) 'The s is out of the bound!'
               call this%err%equit(class//sname//erstr)
               return
            end if
            do i = 2, prof_l
               if (this%s(i) < this%s(i-1)) then
                  write (erstr,*) 's is not monotonically increasing!'
                  call this%err%equit(class//sname//erstr)
                  return
               end if
               if (s<=this%s(i)) then
                  den_temp = this%fs(i-1) + (this%fs(i)-this%fs(i-1))/&
                  &(this%s(i)-this%s(i-1))*(s-this%s(i-1))
                  exit
               end if
            end do
         end if
         qm = den_temp*this%den*this%qm/abs(this%qm)/real(xppc)/real(yppc)
         cx = this%cx; cy = this%cy
         irx = this%irx; iry = this%iry
         orx = this%orx; ory = this%ory
         iirx2 = 1.0/irx**2; iiry2 = 1.0/iry**2
         iorx2 = 1.0/orx**2; iory2 = 1.0/ory**2
         nps = 1
         pt => part2d
! initialize the particle positions
         if (noff < ny) then
         do i=2, nx-1
            do j=2, ny
               do ix = 0, xppc-1
                  do iy=0, yppc-1
                     x = (ix + 0.5)/xppc + i - 1
                     y = (iy + 0.5)/yppc + j - 1 + noff
                     if (((x-cx)**2*iirx2+(y-cy)**2*iiry2 < 1) .or. &
                     &((x-cx)**2*iorx2+(y-cy)**2*iory2 > 1)) then
                        cycle
                     else 
                        pt(1,nps) = x
                        pt(2,nps) = y
                        pt(3,nps) = 0.0
                        pt(4,nps) = 0.0
                        pt(5,nps) = 0.0
                        pt(6,nps) = 1.0
                        pt(7,nps) = 1.0
                        pt(8,nps) = qm
                        nps = nps + 1
                     end if
                  enddo
               enddo
            enddo
         enddo
         else if (noff > (nx-ny-1)) then       
         do i=2, nx-1
            do j=1, ny-1
               do ix = 0, xppc-1
                  do iy=0, yppc-1
                     x = (ix + 0.5)/xppc + i - 1
                     y = (iy + 0.5)/yppc + j - 1 + noff
                     if (((x-cx)**2*iirx2+(y-cy)**2*iiry2 < 1) .or. &
                     &((x-cx)**2*iorx2+(y-cy)**2*iory2 > 1)) then
                        cycle
                     else 
                        pt(1,nps) = x
                        pt(2,nps) = y
                        pt(3,nps) = 0.0
                        pt(4,nps) = 0.0
                        pt(5,nps) = 0.0
                        pt(6,nps) = 1.0
                        pt(7,nps) = 1.0
                        pt(8,nps) = qm
                        nps = nps + 1
                     end if
                  enddo
               enddo
            enddo
         enddo
         else
         do i=2, nx-1
            do j=1, ny
               do ix = 0, xppc-1
                  do iy=0, yppc-1
                     x = (ix + 0.5)/xppc + i - 1
                     y = (iy + 0.5)/yppc + j - 1 + noff
                     if (((x-cx)**2*iirx2+(y-cy)**2*iiry2 < 1) .or. &
                     &((x-cx)**2*iorx2+(y-cy)**2*iory2 > 1)) then
                        cycle
                     else 
                        pt(1,nps) = x
                        pt(2,nps) = y
                        pt(3,nps) = 0.0
                        pt(4,nps) = 0.0
                        pt(5,nps) = 0.0
                        pt(6,nps) = 1.0
                        pt(7,nps) = 1.0
                        pt(8,nps) = qm
                        nps = nps + 1
                     end if
                  enddo
               enddo
            enddo
         enddo
         endif
         
         npp = nps - 1
         
         call this%err%werrfl2(class//sname//' ended')

      end subroutine dist2d_010
!
      subroutine init_fdist2d_011(this,input,i)
      
         implicit none
         
         class(fdist2d_011), intent(inout) :: this
         type(input_json), intent(inout), pointer :: input
         integer, intent(in) :: i
! local data
         integer :: npf,xppc,yppc,npmax,indx,indy
         real :: qm, den
         real :: cx, cy
         real :: r0,sigr
         real :: min, max
         real :: alx, aly, dx, dy
         character(len=20) :: sn,s1
         character(len=18), save :: sname = 'init_fdist2d_011:'
         
         this%sp => input%sp
         this%err => input%err
         this%p => input%pp

         call this%err%werrfl2(class//sname//' started')
         write (sn,'(I3.3)') i
         s1 = 'species('//trim(sn)//')'
         call input%get('simulation.indx',indx)
         call input%get('simulation.indy',indy)
         call input%get('simulation.box.x(1)',min)
         call input%get('simulation.box.x(2)',max)
         call input%get(trim(s1)//'.center(1)',cx)
         cx = cx - min
         alx = (max-min) 
         dx=alx/real(2**indx)
         call input%get('simulation.box.y(1)',min)
         call input%get('simulation.box.y(2)',max)
         call input%get(trim(s1)//'.center(2)',cy)
         cy = cy - min
         aly = (max-min) 
         dy=aly/real(2**indy)
         call input%get(trim(s1)//'.profile',npf)
         call input%get(trim(s1)//'.ppc(1)',xppc)
         call input%get(trim(s1)//'.ppc(2)',yppc)
         call input%get(trim(s1)//'.q',qm)
         call input%get(trim(s1)//'.density',den)
         call input%get(trim(s1)//'.radius',r0)
         call input%get(trim(s1)//'.width',sigr)
         npmax = xppc*yppc*(2**indx)*(2**indy)/this%p%getlnvp()*4
         call input%get(trim(s1)//'.longitudinal_profile',this%long_prof)
         if (trim(this%long_prof) == 'piecewise') then
            call input%get(trim(s1)//'.piecewise_density',this%fs)
            call input%get(trim(s1)//'.piecewise_s',this%s)
         end if

         this%npf = npf
         this%xppc = xppc
         this%yppc = yppc
         this%qm = qm
         this%den = den
         this%npmax = npmax
         this%cx = cx/dx
         this%cy = cy/dy
         this%r0 = r0/dx
         this%sigr = sigr/dx

         call this%err%werrfl2(class//sname//' ended')

      end subroutine init_fdist2d_011
!
      subroutine dist2d_011(this,part2d,npp,fd,s)
         implicit none
         class(fdist2d_011), intent(inout) :: this
         real, dimension(:,:), pointer, intent(inout) :: part2d
         integer, intent(inout) :: npp
         class(ufield2d), intent(in), pointer :: fd
         real, intent(in) :: s 
! local data
         character(len=18), save :: sname = 'dist2d_011:'
         real, dimension(:,:), pointer :: pt => null()
         integer :: nps, nx, ny, noff, xppc, yppc, i, j
         integer :: ix, iy
         real :: qm, x, y
         real :: cx, cy
         real :: r0, sigr, isigr2, rr
         real :: den_temp
         integer :: prof_l

         call this%err%werrfl2(class//sname//' started')
         
         nx = fd%getnd1p(); ny = fd%getnd2p(); noff = fd%getnoff()
         xppc = this%xppc; yppc = this%yppc
         den_temp = 1.0
         if (trim(this%long_prof) == 'piecewise') then
            prof_l = size(this%fs)
            if (s<this%s(1) .or. s>this%s(prof_l)) then
               write (erstr,*) 'The s is out of the bound!'
               call this%err%equit(class//sname//erstr)
               return
            end if
            do i = 2, prof_l
               if (this%s(i) < this%s(i-1)) then
                  write (erstr,*) 's is not monotonically increasing!'
                  call this%err%equit(class//sname//erstr)
                  return
               end if
               if (s<=this%s(i)) then
                  den_temp = this%fs(i-1) + (this%fs(i)-this%fs(i-1))/&
                  &(this%s(i)-this%s(i-1))*(s-this%s(i-1))
                  exit
               end if
            end do
         end if
         qm = den_temp*this%den*this%qm/abs(this%qm)/real(xppc)/real(yppc)
         cx = this%cx; cy = this%cy
         r0 = this%r0; sigr = this%sigr
         isigr2 = 0.5/sigr**2
         nps = 1
         pt => part2d
! initialize the particle positions
         if (noff < ny) then
         do i=2, nx-1
            do j=2, ny
               do ix = 0, xppc-1
                  do iy=0, yppc-1
                     x = (ix + 0.5)/xppc + i - 1
                     y = (iy + 0.5)/yppc + j - 1 + noff
                     rr = sqrt((x-cx)**2+(y-cy)**2)
                     pt(1,nps) = x
                     pt(2,nps) = y
                     pt(3,nps) = 0.0
                     pt(4,nps) = 0.0
                     pt(5,nps) = 0.0
                     pt(6,nps) = 1.0
                     pt(7,nps) = 1.0
                     pt(8,nps) = qm*exp(-(rr-r0)**2*isigr2)
                     nps = nps + 1
                  enddo
               enddo
            enddo
         enddo
         else if (noff > (nx-ny-1)) then       
         do i=2, nx-1
            do j=1, ny-1
               do ix = 0, xppc-1
                  do iy=0, yppc-1
                     x = (ix + 0.5)/xppc + i - 1
                     y = (iy + 0.5)/yppc + j - 1 + noff
                     rr = sqrt((x-cx)**2+(y-cy)**2)
                     pt(1,nps) = x
                     pt(2,nps) = y
                     pt(3,nps) = 0.0
                     pt(4,nps) = 0.0
                     pt(5,nps) = 0.0
                     pt(6,nps) = 1.0
                     pt(7,nps) = 1.0
                     pt(8,nps) = qm*exp(-(rr-r0)**2*isigr2)
                     nps = nps + 1
                  enddo
               enddo
            enddo
         enddo
         else
         do i=2, nx-1
            do j=1, ny
               do ix = 0, xppc-1
                  do iy=0, yppc-1
                     x = (ix + 0.5)/xppc + i - 1
                     y = (iy + 0.5)/yppc + j - 1 + noff
                     rr = sqrt((x-cx)**2+(y-cy)**2)
                     pt(1,nps) = x
                     pt(2,nps) = y
                     pt(3,nps) = 0.0
                     pt(4,nps) = 0.0
                     pt(5,nps) = 0.0
                     pt(6,nps) = 1.0
                     pt(7,nps) = 1.0
                     pt(8,nps) = qm*exp(-(rr-r0)**2*isigr2)
                     nps = nps + 1
                  enddo
               enddo
            enddo
         enddo
         endif
         
         npp = nps - 1
         
         call this%err%werrfl2(class//sname//' ended')

      end subroutine dist2d_011
!
      subroutine init_fdist2d_012(this,input,i)
      
         implicit none
         
         class(fdist2d_012), intent(inout) :: this
         type(input_json), intent(inout), pointer :: input
         integer, intent(in) :: i
! local data
         integer :: npf,xppc,yppc,npmax,indx,indy
         real :: qm, den
         real :: cx, cy
         real :: min, max
         real :: alx, aly, dx, dy
         character(len=20) :: sn,s1
         character(len=18), save :: sname = 'init_fdist2d_012:'
         
         this%sp => input%sp
         this%err => input%err
         this%p => input%pp

         call this%err%werrfl2(class//sname//' started')
         write (sn,'(I3.3)') i
         s1 = 'species('//trim(sn)//')'
         call input%get('simulation.indx',indx)
         call input%get('simulation.indy',indy)
         call input%get('simulation.box.x(1)',min)
         call input%get('simulation.box.x(2)',max)
         call input%get(trim(s1)//'.center(1)',cx)
         cx = cx - min
         alx = (max-min) 
         dx=alx/real(2**indx)
         call input%get('simulation.box.y(1)',min)
         call input%get('simulation.box.y(2)',max)
         call input%get(trim(s1)//'.center(2)',cy)
         cy = cy - min
         aly = (max-min) 
         dy=aly/real(2**indy)
         call input%get(trim(s1)//'.profile',npf)
         call input%get(trim(s1)//'.ppc(1)',xppc)
         call input%get(trim(s1)//'.ppc(2)',yppc)
         call input%get(trim(s1)//'.q',qm)
         call input%get(trim(s1)//'.density',den)
         npmax = xppc*yppc*(2**indx)*(2**indy)/this%p%getlnvp()*4
         call input%get(trim(s1)//'.longitudinal_profile',this%long_prof)
         if (trim(this%long_prof) == 'piecewise') then
            call input%get(trim(s1)//'.piecewise_density',this%fs)
            call input%get(trim(s1)//'.piecewise_s',this%s)
         end if
         call input%get(trim(s1)//'.piecewise_radial_density',this%fr)
         call input%get(trim(s1)//'.piecewise_r',this%r)

         this%r = this%r/dx
         this%npf = npf
         this%xppc = xppc
         this%yppc = yppc
         this%qm = qm
         this%den = den
         this%npmax = npmax
         this%cx = cx/dx
         this%cy = cy/dy

         call this%err%werrfl2(class//sname//' ended')

      end subroutine init_fdist2d_012
!
      subroutine dist2d_012(this,part2d,npp,fd,s)
         implicit none
         class(fdist2d_012), intent(inout) :: this
         real, dimension(:,:), pointer, intent(inout) :: part2d
         integer, intent(inout) :: npp
         class(ufield2d), intent(in), pointer :: fd
         real, intent(in) :: s 
! local data
         character(len=18), save :: sname = 'dist2d_012:'
         real, dimension(:,:), pointer :: pt => null()
         integer :: nps, nx, ny, noff, xppc, yppc, i, j
         integer :: ix, iy
         real :: qm, x, y
         real :: cx, cy
         real :: den_temp, rr
         integer :: prof_l, ii

         call this%err%werrfl2(class//sname//' started')
         
         nx = fd%getnd1p(); ny = fd%getnd2p(); noff = fd%getnoff()
         xppc = this%xppc; yppc = this%yppc
         den_temp = 1.0
         if (trim(this%long_prof) == 'piecewise') then
            prof_l = size(this%fs)
            if (prof_l /= size(this%s)) then
               write (erstr,*) 'The piecewise_density and s array have different sizes!'
               call this%err%equit(class//sname//erstr)
               return
            end if
            if (s<this%s(1) .or. s>this%s(prof_l)) then
               write (erstr,*) 'The s is out of the bound!'
               call this%err%equit(class//sname//erstr)
               return
            end if
            do i = 2, prof_l
               if (this%s(i) < this%s(i-1)) then
                  write (erstr,*) 's is not monotonically increasing!'
                  call this%err%equit(class//sname//erstr)
                  return
               end if
               if (s<=this%s(i)) then
                  den_temp = this%fs(i-1) + (this%fs(i)-this%fs(i-1))/&
                  &(this%s(i)-this%s(i-1))*(s-this%s(i-1))
                  exit
               end if
            end do
         end if
         qm = den_temp*this%den*this%qm/abs(this%qm)/real(xppc)/real(yppc)
         cx = this%cx; cy = this%cy
         nps = 1
         pt => part2d
         prof_l = size(this%fr)
         if (prof_l /= size(this%r)) then
            write (erstr,*) 'The piecewise_radial_density and r array have different sizes!'
            call this%err%equit(class//sname//erstr)
            return
         end if
! initialize the particle positions
         if (noff < ny) then
         do i=2, nx-1
            do j=2, ny
               do ix = 0, xppc-1
                  do iy=0, yppc-1
                     x = (ix + 0.5)/xppc + i - 1
                     y = (iy + 0.5)/yppc + j - 1 + noff
                     rr = sqrt((x-cx)**2+(y-cy)**2)
                     if (rr<this%r(1) .or. rr>this%r(prof_l)) then
                        cycle
                     end if
                     do ii = 2, prof_l
                        if (this%r(ii) <= this%r(ii-1)) then
                           write (erstr,*) 'r is not monotonically increasing!'
                           call this%err%equit(class//sname//erstr)
                           return
                        end if
                        if (rr<=this%r(ii)) then
                           den_temp = this%fr(ii-1) + (this%fr(ii)-this%fr(ii-1))/&
                           &(this%r(ii)-this%r(ii-1))*(rr-this%r(ii-1))
                           exit
                        end if
                     end do
                     pt(1,nps) = x
                     pt(2,nps) = y
                     pt(3,nps) = 0.0
                     pt(4,nps) = 0.0
                     pt(5,nps) = 0.0
                     pt(6,nps) = 1.0
                     pt(7,nps) = 1.0
                     pt(8,nps) = qm*den_temp
                     nps = nps + 1
                  enddo
               enddo
            enddo
         enddo
         else if (noff > (nx-ny-1)) then       
         do i=2, nx-1
            do j=1, ny-1
               do ix = 0, xppc-1
                  do iy=0, yppc-1
                     x = (ix + 0.5)/xppc + i - 1
                     y = (iy + 0.5)/yppc + j - 1 + noff
                     rr = sqrt((x-cx)**2+(y-cy)**2)
                     if (rr<this%r(1) .or. rr>this%r(prof_l)) then
                        cycle
                     end if
                     do ii = 2, prof_l
                        if (this%r(ii) <=this%r(ii-1)) then
                           write (erstr,*) 'r is not monotonically increasing!'
                           call this%err%equit(class//sname//erstr)
                           return
                        end if
                        if (rr<=this%r(ii)) then
                           den_temp = this%fr(ii-1) + (this%fr(ii)-this%fr(ii-1))/&
                           &(this%r(ii)-this%r(ii-1))*(rr-this%r(ii-1))
                           exit
                        end if
                     end do
                     pt(1,nps) = x
                     pt(2,nps) = y
                     pt(3,nps) = 0.0
                     pt(4,nps) = 0.0
                     pt(5,nps) = 0.0
                     pt(6,nps) = 1.0
                     pt(7,nps) = 1.0
                     pt(8,nps) = qm*den_temp
                     nps = nps + 1
                  enddo
               enddo
            enddo
         enddo
         else
         do i=2, nx-1
            do j=1, ny
               do ix = 0, xppc-1
                  do iy=0, yppc-1
                     x = (ix + 0.5)/xppc + i - 1
                     y = (iy + 0.5)/yppc + j - 1 + noff
                     rr = sqrt((x-cx)**2+(y-cy)**2)
                     if (rr<this%r(1) .or. rr>this%r(prof_l)) then
                        cycle
                     end if
                     do ii = 2, prof_l
                        if (this%r(ii) < this%r(ii-1)) then
                           write (erstr,*) 'r is not monotonically increasing!'
                           call this%err%equit(class//sname//erstr)
                           return
                        end if
                        if (rr<=this%r(ii)) then
                           den_temp = this%fr(ii-1) + (this%fr(ii)-this%fr(ii-1))/&
                           &(this%r(ii)-this%r(ii-1))*(rr-this%r(ii-1))
                           exit
                        end if
                     end do
                     pt(1,nps) = x
                     pt(2,nps) = y
                     pt(3,nps) = 0.0
                     pt(4,nps) = 0.0
                     pt(5,nps) = 0.0
                     pt(6,nps) = 1.0
                     pt(7,nps) = 1.0
                     pt(8,nps) = qm*den_temp
                     nps = nps + 1
                  enddo
               enddo
            enddo
         enddo
         endif
         
         npp = nps - 1
         
         call this%err%werrfl2(class//sname//' ended')

      end subroutine dist2d_012
!

      subroutine init_fdist2d_013(this,input,i)
      
         implicit none
         
         class(fdist2d_013), intent(inout) :: this
         type(input_json), intent(inout), pointer :: input
         integer, intent(in) :: i
! local data
         integer :: npf,xppc,yppc,npmax,indx,indy
         real :: qm, den
         real :: cx, cy
         real :: min, max
         real :: alx, aly, dx, dy
         character(len=20) :: sn,s1
         character(len=18), save :: sname = 'init_fdist2d_013:'
         character(len=:), allocatable :: read_str
         integer :: ierr
         this%sp => input%sp
         this%err => input%err
         this%p => input%pp

         call this%err%werrfl2(class//sname//' started')
         write (sn,'(I3.3)') i
         s1 = 'species('//trim(sn)//')'
         call input%get('simulation.indx',indx)
         call input%get('simulation.indy',indy)
         call input%get('simulation.box.x(1)',min)
         call input%get('simulation.box.x(2)',max)
         cx = - min
         alx = (max-min) 
         dx=alx/real(2**indx)
         call input%get('simulation.box.y(1)',min)
         call input%get('simulation.box.y(2)',max)
         cy = - min
         aly = (max-min) 
         dy=aly/real(2**indy)
         call input%get(trim(s1)//'.profile',npf)
         call input%get(trim(s1)//'.ppc(1)',xppc)
         call input%get(trim(s1)//'.ppc(2)',yppc)
         call input%get(trim(s1)//'.q',qm)
         call input%get(trim(s1)//'.density',den)
         npmax = xppc*yppc*(2**indx)*(2**indy)/this%p%getlnvp()*4

         this%dx = dx
         this%dy = dy
         this%npf = npf
         this%xppc = xppc
         this%yppc = yppc
         this%qm = qm
         this%den = den
         this%npmax = npmax
         this%cx = cx/dx
         this%cy = cy/dy

         if ( .not. associated( this%math_func ) ) then
          allocate( t_fparser :: this%math_func )
        endif
        call input%get( trim(s1) // '.math_func', read_str )
        call setup(this%math_func, trim(read_str), (/'x','y','s'/), ierr)

        call this%err%werrfl2(class//sname//' ended')
         

      end subroutine init_fdist2d_013
!
      subroutine dist2d_013(this,part2d,npp,fd,s)
         implicit none
         class(fdist2d_013), intent(inout) :: this
         real, dimension(:,:), pointer, intent(inout) :: part2d
         integer, intent(inout) :: npp
         class(ufield2d), intent(in), pointer :: fd
         real, intent(in) :: s 
! local data
         character(len=18), save :: sname = 'dist2d_013:'
         real, dimension(:,:), pointer :: pt => null()
         integer :: nps, nx, ny, noff, xppc, yppc, i, j
         integer :: ix, iy
         real :: qm, x, y
         real :: cx, cy
         real :: den_temp, rr
         integer :: prof_l, ii
         real(p_k_fparse), dimension(3) :: fparser_arr 

         call this%err%werrfl2(class//sname//' started')
         
         nx = fd%getnd1p(); ny = fd%getnd2p(); noff = fd%getnoff()
         xppc = this%xppc; yppc = this%yppc
         den_temp = 1.0
         qm = den_temp*this%den*this%qm/abs(this%qm)/real(xppc)/real(yppc)
         cx = this%cx; cy = this%cy
         nps = 1
         pt => part2d
         fparser_arr(3) = s
! initialize the particle positions
         if (noff < ny) then
         do i=2, nx-1
            do j=2, ny
               do ix = 0, xppc-1
                  do iy=0, yppc-1
                     x = (ix + 0.5)/xppc + i - 1
                     y = (iy + 0.5)/yppc + j - 1 + noff
                     fparser_arr(1) = (x-cx) * this%dx
                     fparser_arr(2) = (y-cy) * this%dy
                     den_temp = eval(this%math_func, fparser_arr)
                     pt(1,nps) = x
                     pt(2,nps) = y
                     pt(3,nps) = 0.0
                     pt(4,nps) = 0.0
                     pt(5,nps) = 0.0
                     pt(6,nps) = 1.0
                     pt(7,nps) = 1.0
                     pt(8,nps) = qm*den_temp
                     nps = nps + 1
                  enddo
               enddo
            enddo
         enddo
         else if (noff > (nx-ny-1)) then       
         do i=2, nx-1
            do j=1, ny-1
               do ix = 0, xppc-1
                  do iy=0, yppc-1
                     x = (ix + 0.5)/xppc + i - 1
                     y = (iy + 0.5)/yppc + j - 1 + noff
                     fparser_arr(1) = (x-cx) * this%dx
                     fparser_arr(2) = (y-cy) * this%dy
                     den_temp = eval(this%math_func, fparser_arr)
                     pt(1,nps) = x
                     pt(2,nps) = y
                     pt(3,nps) = 0.0
                     pt(4,nps) = 0.0
                     pt(5,nps) = 0.0
                     pt(6,nps) = 1.0
                     pt(7,nps) = 1.0
                     pt(8,nps) = qm*den_temp
                     nps = nps + 1
                  enddo
               enddo
            enddo
         enddo
         else
         do i=2, nx-1
            do j=1, ny
               do ix = 0, xppc-1
                  do iy=0, yppc-1
                     x = (ix + 0.5)/xppc + i - 1
                     y = (iy + 0.5)/yppc + j - 1 + noff
                     fparser_arr(1) = (x-cx) * this%dx
                     fparser_arr(2) = (y-cy) * this%dy
                     den_temp = eval(this%math_func, fparser_arr)
                     pt(1,nps) = x
                     pt(2,nps) = y
                     pt(3,nps) = 0.0
                     pt(4,nps) = 0.0
                     pt(5,nps) = 0.0
                     pt(6,nps) = 1.0
                     pt(7,nps) = 1.0
                     pt(8,nps) = qm*den_temp
                     nps = nps + 1
                  enddo
               enddo
            enddo
         enddo
         endif
         
         npp = nps - 1
         
         call this%err%werrfl2(class//sname//' ended')

      end subroutine dist2d_013

      end module fdist2d_class