module ldist3d_class

    use perrors_class
    use parallel_pipe_class
    use spect3d_class
    use ufield3d_class
    use part3d_lib
    use input_class
    use mpi
           
    implicit none
  
    private
  
    public :: ldist3d, ldist3d_000
  
    interface
        real function ranorm()
        implicit none
        end function
    end interface
  
    type, abstract :: ldist3d
  
        private
  
        class(spect3d), pointer, public :: sp => null()
        class(perrors), pointer, public :: err => null()
        class(parallel_pipe), pointer, public :: p => null()
        integer :: npf
        real, dimension(3) :: origin=(/0.0,0.0,0.0/)
        logical :: evol = .true.
                           
        contains
           
        generic :: new => init_ldist3d         
        !generic :: del => end_ldist3d
        generic :: dist => dist3d
        procedure(ab_init_ldist3d), deferred, private :: init_ldist3d
        procedure(ab_dist3d), deferred, private :: dist3d
        !procedure, private :: end_ldist3d
        !procedure :: getnpf, getevol, getorigin
                    
    end type 
  
    abstract interface

    subroutine ab_init_ldist3d(this,input,i)
        import ldist3d
        import input_json
        implicit none
        class(ldist3d), intent(inout) :: this
        type(input_json), intent(inout), pointer :: input
        integer, intent(in) :: i
    end subroutine ab_init_ldist3d

    subroutine ab_dist3d(this,lrf)
        import ldist3d
        import ufield3d
        implicit none
        class(ldist3d), intent(in) :: this
        class(ufield3d), intent(inout) :: lrf        
    end subroutine ab_dist3d

    end interface

    type, extends(ldist3d) :: ldist3d_000
        ! Tri Gaussian profile (the same particle charge)
        private

        real :: sigx, sigy, sigz, lcx, lcy, lcz, nx, ny, nz
        real :: peak_density

        contains
        
        procedure, private :: init_ldist3d => init_ldist3d_000
        procedure, private :: dist3d => dist3d_000
  
    end type ldist3d_000
  
    character(len=10), save :: class = 'ldist3d:'
    character(len=128), save :: erstr
        
    contains

    subroutine init_ldist3d_000(this,input,i)
        
        implicit none
           
        class(ldist3d_000), intent(inout) :: this
        type(input_json), intent(inout), pointer :: input
        integer, intent(in) :: i        
        integer :: npf, indx, indy, indz, nx, ny, nz
        real :: sigx,sigy,sigz, alx, aly, alz, dx, dy, dz, lcx, lcy, lcz
        real :: min, max, peak_density
        logical :: evol
  
        character(len=20) :: sn,s1
        character(len=18), save :: sname = 'init_ldist3d_000:'
           
        this%sp => input%sp
        this%err => input%err
        this%p => input%pp

        call this%err%werrfl2(class//sname//' started')
        write (sn,'(I3.3)') i
        s1 = 'beam('//trim(sn)//')'
  
        call input%get('simulation.indx',indx)
        call input%get('simulation.indy',indy)
        call input%get('simulation.indz',indz)
        nx = 2**indx
        ny = 2**indy
        nz = 2**indz
        call input%get('simulation.box.x(1)',min)
        call input%get('simulation.box.x(2)',max)
        this%origin(1) = min
        call input%get(trim(s1)//'.center(1)',lcx)
        lcx = lcx - min
        alx = (max-min) 
        dx = alx/real(nx)
        call input%get('simulation.box.y(1)',min)
        call input%get('simulation.box.y(2)',max)
        this%origin(2) = min
        call input%get(trim(s1)//'.center(2)',lcy)
        lcy = lcy -min
        aly = (max-min) 
        dy = aly/real(ny)
        call input%get('simulation.box.z(1)',min)
        call input%get('simulation.box.z(2)',max)
        this%origin(3) = min
        call input%get(trim(s1)//'.center(3)',lcz)
        lcz = lcz -min
        alz = (max-min) 
        dz = alz/real(nz)
  
        call input%get(trim(s1)//'.profile',npf)
        call input%get(trim(s1)//'.sigma(1)',sigx)
        call input%get(trim(s1)//'.sigma(2)',sigy)
        call input%get(trim(s1)//'.sigma(3)',sigz)
        call input%get(trim(s1)//'.peak_density',peak_density)
        call input%get(trim(s1)//'.evolution',evol)
        if (input%found(trim(s1)//'.static')) then
            call input%get(trim(s1)//'.static',static)
            if (static) evol = .false.
        else
            static = .false.
        end if
  
        this%npf = npf
        this%nx = nx
        this%ny = ny
        this%nz = nz
        this%sigx = sigx
        this%sigy = sigy
        this%sigz = sigz
        this%lcx = lcx
        this%lcy = lcy
        this%lcz = lcz
        this%peak_density = peak_density
        this%evol = evol
  
        call this%err%werrfl2(class//sname//' ended')
  
    end subroutine init_ldist3d_000

    subroutine dist3d_000(this,lrf)
    
        implicit none
        
        class(ldist3d_000), intent(in) :: this
        class(ufield3d), intent(inout) :: lrf
        real, dimension(:,:,:,:), pointer :: lf => null()        
        integer :: nx, ny, nz, i, j, k, x, y, z
        real :: sigx, sigy, sigz, lcx, lcy, lcz

        character(len=18), save :: sname = 'dist3d_000:'
        call this%err%werrfl2(class//sname//' started')
        
        nx = this%nx
        ny = this%ny
        nz = this%nz
        sigx = this%sigx
        sigy = this%sigy
        sigz = this%sigz
        lcx = this%lcx 
        lcy = this%lcy
        lcz = this%lcz
        peak_density = this%peak_density
        lf => lrf%getrf()

        do i = 1, nx
            x = i*dx 
            do j = 1, ny
                y = j*dy
                do k = 1, nz
                    z = k*dz
                    lf(1,i,j,k) = peak_density*exp(-(x-lcx)**2/2.0/sigx**2)&
                    *exp(-(y-lcy)**2/2.0/sigy**2)*exp(-(z-lcz)**2/2.0/sigz**2)
                enddo
            enddo    
        enddo

        call this%err%werrfl2(class//sname//' ended')
        
    end subroutine dist3d_000

    function getnpf(this)
  
        implicit none
  
        class(ldist3d), intent(in) :: this
        integer :: getnpf
           
        getnpf = this%npf
  
    end function getnpf
      
    function getnpmax(this)
  
        implicit none
  
        class(ldist3d), intent(in) :: this
        integer :: getnpmax
           
        getnpmax = this%npmax
  
    end function getnpmax
     
    function getevol(this)
  
        implicit none
  
        class(ldist3d), intent(in) :: this
        logical :: getevol
           
        getevol = this%evol
  
    end function getevol
      
    function getstatic(this)
  
        implicit none
  
        class(ldist3d), intent(in) :: this
        logical :: getstatic
           
        getstatic = this%static
  
    end function getstatic

    function getorigin(this)
  
        implicit none
  
        class(ldist3d), intent(in) :: this
        real, dimension(3) :: getorigin
           
        getorigin = this%origin
  
    end function getorigin

    subroutine end_ldist3d(this)
            
        implicit none
           
        class(ldist3d), intent(inout) :: this
        character(len=18), save :: sname = 'end_ldist3d:'
  
        call this%err%werrfl2(class//sname//' started')
        call this%err%werrfl2(class//sname//' ended')
                    
    end subroutine end_ldist3d

end module ldist3d_class