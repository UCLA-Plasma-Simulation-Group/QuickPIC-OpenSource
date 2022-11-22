module laser_class

    use perrors_class
    use parallel_pipe_class
    use spect3d_class
    use ldist3d_class
    use field3d_class
    use field2d_class
    use ufield3d_class
    use ufield2d_class
    use hdf5io_class
    use mpi
                     
    implicit none
      
    private 
    public :: laser
      
    type laser
      
        private
      
        class(spect3d), pointer, public :: sp => null()
        class(perrors), pointer, public :: err => null()
        class(parallel_pipe), pointer, public :: p => null()
        class(field3d), pointer :: lp => null()
        class(ldist3d), pointer :: lf => null() !laser
        class(field3d), pointer :: a_n => null(), a_half => null() !The envelope of integer order and half-integer order, the integer order communicates with the plasma, and the half-integer order is responsible for the circulation
        class(field2d), pointer :: a_fft => null() !for FFT
        real :: k0, ds, dz, dx, beta1, beta2 !extra parameters of envelope equation
    
        contains
               
        generic :: new => init_laser 
        generic :: del => end_laser
        generic :: push => push_laser
        generic :: wrlp => writelp_laser, writelpslice_laser
        ! generic :: wrst => writerst_laser       
        ! generic :: rrst => readrst_laser       
        ! procedure :: getlp 
        procedure, private :: init_laser
        procedure, private :: end_laser
        procedure, private :: push_laser
        ! procedure, private :: writerst_laser, readrst_laser
        procedure, private :: writelp_laser, writelpslice_laser 
    
    end type laser
      
    save      
      
    character(len=10) :: class = 'laser:'
    character(len=128) :: erstr
            
    contains
    
    subroutine push_laser(this,chi)
            
        implicit none
               
        class(laser), intent(inout) :: this
        class(field3d), intent(in) :: chi
        
        character(len=18), save :: sname = 'push_laser:'
      
        type(ufield3d), pointer :: f3d => null()
        type(ufield2d), pointer :: f2d => null()
        real, pointer, dimension(:,:,:,:) :: a0, a1, pchi => null()
        real, pointer, dimension(:,:,:) :: a_fftreal => null() 
        integer :: i, j, k, nx, ny, nz
        real :: k0, ds, dz, dx, beta1, beta2
        real :: mid1, mid2, mid3, mid4, c4r, c4i, c6r, c6i!Intermediate variables
        real, dimension(:,:,:), allocatable :: aa !matrix for dual sweep 2*nz-2*3
        real, dimension(:,:), allocatable :: c5, bb, a_single
        
        call this%err%werrfl2(class//sname//' started')
         
        k0 = this%k0
        ds = this%ds
        dz = this%dz
        dx = this%dx
        beta1 = this%beta1
        beta2 = this%beta2
    
        f3d => this%a_half%getrs() 
        nx = f3d%getnd1p()
        ny = f3d%getnd2p()
        nz = f3d%getnd3p()
        a0 => f3d%getrf()

        f3d => chi%getrs()
        pchi => f3d%getrf()

        allocate(aa(2,nz-2,3),bb(2,nz-2),c5(2,nz-2),a_single(2,nz-2))

        mid1 = 1/(dz*ds)
        mid2 = beta2/(4*dz*dz) !need beta1 and beta2 in type laser
        mid3 = (k0*beta1)/(8*dz)
        mid4 = 0.5*k0*k0
        
        aa(1,:,1) = -mid1+mid2
        aa(2,:,1) = mid3
        aa(1,:,3) = mid1+mid2
        aa(2,:,3) = -mid3

        c4r = -mid1-mid2
        c4i = -mid3
        c6r = mid1-mid2
        c6i = mid3
    
        do i = 1,nx
            do j = 1,ny
                do k = 1,nz-2 
                    aa(1,k,2)=-2*mid2-mid4*pchi(1,i,j,k)!predicted chiT 
                    aa(2,k,2)=-2*k0/ds
                    c5(1,k)=2*mid2+mid4*pchi(1,i,j,k)
                    c5(2,k)=aa(2,k,2)
                    bb(1,k)=(c4r*a0(1,i,j,k)-c4i*a0(2,i,j,k))+(c5(1,k)*a0(1,i,j,k+1)&
                    -c5(2,k)*a0(2,i,j,k+1))+(c6r*a0(1,i,j,k+2)-c6i*a0(2,i,j,k+2))
                    bb(2,k)=(c4r*a0(2,i,j,k)+c4i*a0(1,i,j,k))+(c5(1,k)+a0(2,i,j,k+1)&
                    +c5(2,k)*a0(1,i,j,k+1))+(c6r*a0(2,i,j,k+2)+c6i*a0(1,i,j,k+2))
                enddo                      
                call dual_sweep(aa(1,:,:),aa(2,:,:),bb(1,:),bb(2,:),a_single(1,2:nz-1),a_single(2,2:nz-1),nz-2) !dual sweep can be further optimized
                
                a_single(1,1)=0
                a_single(1,nz)=0
                a_single(2,1)=0
                a_single(2,nz)=0

                a0(1,i,j,:)=a_single(1,:)
                a0(2,i,j,:)=a_single(2,:)
            enddo
        enddo

        f2d => this%a_fft%getrs()
        a_fftreal => f2d%getrf()
        
        do i = 1,nz 
            a_fftreal(1,:,:) = a0(1,:,:,i)
            call this%a_fft%fftrk(1)
            f2d => this%a_fft%getks()
            a_fftreal => f2d%getrf()
            a0(1,:,:,i)=a_fftreal(1,:,:)   
            
            a_fftreal(1,:,:) = a0(2,:,:,i)
            call this%a_fft%fftrk(1)
            f2d => this%a_fft%getks()
            a_fftreal => f2d%getrf()
            a0(2,:,:,i)=a_fftreal(1,:,:)
        enddo
        
        do i=1,nx
            do j=1,ny !nduniform?
                !kk = aimag(a0r%pt%ffc(i,j,1))/real(a0r%pt%ffc(i,j,1))/idx2*0.25!kk是一个数组，此处实际有查表操作
                do k=1,nz-2
                    aa(1,k,2)=-2*mid2+0.5*(kk)!kk有问题
                    c5(1,k)=2*mid2-0.5*(kk)
                    bb(1,k)=(c4r*a0(1,i,j,k)-c4i*a0(2,i,j,k))+(c5(1,k)*a0(1,i,j,k+1)&
                    -c5(2,k)*a0(2,i,j,k+1))+(c6r*a0(1,i,j,k+2)-c6i*a0(2,i,j,k+2))
                    bb(2,k)=(c4r*a0(2,i,j,k)+c4i*a0(1,i,j,k))+(c5(1,k)+a0(2,i,j,k+1)&
                    +c5(2,k)*a0(1,i,j,k+1))+(c6r*a0(2,i,j,k+2)+c6i*a0(1,i,j,k+2))
                enddo

                call dual_sweep(aa(1,:,:),aa(2,:,:),bb(1,:),bb(2,:),a_single(1,2:nz-1),a_single(2,2:nz-1),nz-2)  
                a0(1,i,j,:)=a_single(1,:)
                a0(2,i,j,:)=a_single(2,:) 
            enddo
        enddo

        do i = 1,nz 
            a_fftreal(1,:,:) = a0(1,:,:,i)
            call this%a_fft%fftkr(1)
            f2d => this%a_fft%getrs()
            a_fftreal => f2d%getrf()
            a0(1,:,:,i)=a_fftreal(1,:,:)   

            a_fftreal(1,:,:) = a0(2,:,:,i)
            call this%a_fft%fftkr(1)
            f2d => this%a_fft%getrs()
            a_fftreal => f2d%getrf()
            a0(2,:,:,i)=a_fftreal(1,:,:)
        enddo
        
        f3d => this%a_n%getrs() 
        a1 => f3d%getrf()
        
        aa(1,:,1)=aa(1,:,1)-mid1 
        aa(1,:,3)=aa(1,:,3)+mid1
        c4r=c4r-mid1
        c6r=c6r+mid1
    
        do i = 1,nx
            do j = 1,ny
                do k = 1,nz-2 
                    aa(1,k,2)=-2*mid2+mid4*pchi(1,i,j,k)
                    aa(2,k,2)=-4*k0/ds
                    c5(1,k)=2*mid2+mid4*pchi(1,i,j,k)
                    c5(2,k)=aa(2,k,2)
                    bb(1,k)=(c4r*a0(1,i,j,k)-c4i*a0(2,i,j,k))+(c5(1,k)*a0(1,i,j,k+1)&
                    -c5(2,k)*a0(2,i,j,k+1))+(c6r*a0(1,i,j,k+2)-c6i*a0(2,i,j,k+2))
                    bb(2,k)=(c4r*a0(2,i,j,k)+c4i*a0(1,i,j,k))+(c5(1,k)+a0(2,i,j,k+1)&
                    +c5(2,k)*a0(1,i,j,k+1))+(c6r*a0(2,i,j,k+2)+c6i*a0(1,i,j,k+2))
                enddo
                call dual_sweep(aa(1,:,:),aa(2,:,:),bb(1,:),bb(2,:),a_single(1,2:nz-1),a_single(2,2:nz-1),nz-2) 
                a1(1,i,j,:)=a_single(1,:)
                a1(2,i,j,:)=a_single(2,:) 
            enddo
        enddo   
        !After return a, should be able to get a real chi, and need to make a comparison, confirm the convergence through iteration      
        !deallocate() ?
        call this%err%werrfl2(class//sname//' ended')     
    end subroutine push_laser
    
    subroutine dual_sweep(aa,ai,bb,bi,xx,xi,n)
        implicit none
        integer, intent(in) :: n
        real, intent(in), dimension (n,3) :: aa,ai
        real, intent(in), dimension (n) :: bb,bi
        real, intent(out), dimension(n) :: xx,xi
        real, dimension(:), allocatable :: cc,ci,dd,di
        
        integer :: i
        real :: m,mi,w,wi,p,pi,q,qi    
        allocate(cc(n),ci(n),dd(n),di(n))
    
        cc(1)=(aa(1,3)*aa(1,2)+ai(1,3)*ai(1,2))/(aa(1,2)*aa(1,2)+ai(1,2)*ai(1,2))
        ci(1)=(ai(1,3)*aa(1,2)-aa(1,3)*ai(1,2))/(aa(1,2)*aa(1,2)+ai(1,2)*ai(1,2))
        dd(1)=(bb(1)*aa(1,2)+bi(1)*ai(1,2))/(aa(1,2)*aa(1,2)+ai(1,2)*ai(1,2))
        di(1)=(bi(1)*aa(1,2)-bb(1)*ai(1,2))/(aa(1,2)*aa(1,2)+ai(1,2)*ai(1,2))
    
        do i = 2 , n-1
            m=aa(i,1)*cc(i-1)-ai(i,1)*ci(i-1)
            mi=aa(i,1)*ci(i-1)+ai(i,1)*cc(i-1)
            w=aa(i,2)-m
            wi=ai(i,2)-mi
            cc(i)=(aa(i,3)*w+ai(i,3)*wi)/(w*w+wi*wi)
            ci(i)=(ai(i,3)*w-aa(i,3)*wi)/(w*w+wi*wi)
            p=aa(i,1)*dd(i-1)-ai(i,1)*di(i-1)
            pi=aa(i,1)*di(i-1)+ai(i,1)*dd(i-1)
            q=bb(i)-p
            qi=bi(i)-pi
            dd(i)=(q*w+qi*wi)/(w*w+wi*wi)
            di(i)=(qi*w-q*wi)/(w*w+wi*wi)
        enddo
    
        m=aa(n,1)*cc(n-1)-ai(n,1)*ci(n-1)
        mi=aa(n,1)*ci(n-1)+ai(n,1)*cc(n-1)
        w=aa(n,2)-m
        wi=ai(n,2)-mi
        p=aa(n,1)*dd(n-1)-ai(n,1)*di(n-1)
        pi=aa(n,1)*di(n-1)+ai(n,1)*dd(n-1)
        q=bb(n)-p
        qi=bi(n)-pi
        dd(n)=(q*w+qi*wi)/(w*w+wi*wi)
        di(n)=(qi*w-q*wi)/(w*w+wi*wi)
    
        xx(n)=dd(n)
        xi(n)=di(n)
        do i = n-1 , 1 , -1
            m=cc(i)*xx(i+1)-ci(i)*xi(i+1)
            mi=cc(i)*xi(i+1)+ci(i)*xx(i+1)
            xx(i)=dd(i)-m
            xi(i)=di(i)-mi
        enddo
    
    end subroutine dual_sweep
    
    subroutine init_laser(this,pp,perr,psp,lf,k0,ds,dz,dx,beta1,beta2)
            
        implicit none
               
        class(laser), intent(inout) :: this
        class(spect3d), intent(in), pointer :: psp
        class(perrors), intent(in), pointer :: perr
        class(parallel_pipe), intent(in), pointer :: pp
        class(ldist3d), intent(inout), target :: lf
        real, intent(in) :: k0, ds, dz, dx, beta1, beta2
      
        character(len=18), save :: sname = 'init_laser:'
        
        integer :: id, ierr
        integer, dimension(10) :: istat
                        
        this%sp => psp
        this%err => perr
        this%p => pp
        this%lf => lf
      
        call this%err%werrfl2(class//sname//' started')
      
        this%k0 = k0
        this%ds = ds
        this%dz = dz
        this%dx = dx
        this%beta1 = beta1
        this%beta2 = beta2

        allocate(this%lp,this%a_n,this%a_half,this%a_fft)

        call this%lp%new(this%p,this%err,this%sp,dim=1)
        call this%a_n%new(this%p,this%err,this%sp,dim=2)
        call this%a_half%new(this%p,this%err,this%sp,dim=2)
        call this%a_fft%new(this%p,this%err,this%sp,dim=2,fftflag=.true.)
        call this%lf%dist(this%a_n%getrs())
        call this%lf%dist(this%a_half%getrs())

        call this%err%werrfl2(class//sname//' ended')
      
    end subroutine init_laser
            
    subroutine end_laser(this)
                
        implicit none
               
        class(laser), intent(inout) :: this
        character(len=18), save :: sname = 'end_laser:'
      
        call this%err%werrfl2(class//sname//' started')
        
        call this%lp%del()
        call this%a_n%del()
        call this%a_half%del()
        call this%a_fft%del()

        call this%err%werrfl2(class//sname//' ended')
                        
    end subroutine end_laser
          
    !subroutine pmove_beam3d(this,fd,rtag,stag,sid)
    !subroutine init_laser(this,pp,perr,psp,lf,k0,ds,dz,dx)
            
        !     implicit none
                   
        !     class(laser), intent(inout) :: this
        !     class(spect3d), intent(in), pointer :: psp
        !     class(perrors), intent(in), pointer :: perr
        !     class(parallel_pipe), intent(in), pointer :: pp
        !     class(ldist3d), intent(inout), target :: lf
        !     real, intent(in) :: k0, ds, dx
          
        !     character(len=18), save :: sname = 'init_laser:'
        !     integer :: id, ierr
        !     integer, dimension(10) :: istat
                            
        !     this%sp => psp
        !     this%err => perr
        !     this%p => pp
        !     this%lf => lf
          
        !     call this%err%werrfl2(class//sname//' started')
          
        !     this%k0 = k0
        !     this%ds = ds
        !     this%dz = dz
        !     this%dx = dx
        !     allocate(this%lp,this%a0x,this%a0y,this%a0z)!这是声明了一个什么数组？
        !     allocate(this%a1x,this%a1y,this%a1z)
        !     allocate(this%a2x,this%a2y,this%a2z)
        !     allocate(this%a0r,this%a0i)
        !     call this%lp%new(this%p,this%err,this%sp,dim=1)
        !     call this%a0x%new(this%p,this%err,this%sp,dim=2)
        !     call this%a0y%new(this%p,this%err,this%sp,dim=2)
        !     call this%a0z%new(this%p,this%err,this%sp,dim=2)
        !     call this%a1x%new(this%p,this%err,this%sp,dim=2)
        !     call this%a1y%new(this%p,this%err,this%sp,dim=2)
        !     call this%a1z%new(this%p,this%err,this%sp,dim=2)
        !     call this%a2x%new(this%p,this%err,this%sp,dim=2)
        !     call this%a2y%new(this%p,this%err,this%sp,dim=2)
        !     call this%a2z%new(this%p,this%err,this%sp,dim=2)
        !     call this%a0r%new(this%p,this%err,this%sp,dim=1,fftflag=.true.)
        !     call this%a0i%new(this%p,this%err,this%sp,dim=1,fftflag=.true.)
        !     call this%lf%dist%(this%a0x,this%a0y,this%a0z)
                   
        !     call this%err%werrfl2(class//sname//' ended')
          
        ! end subroutine init_laser
                
        ! subroutine end_beam3d(this)
                    
        !     implicit none
                   
        !     class(beam3d), intent(inout) :: this
        !     character(len=18), save :: sname = 'end_beam3d:'
          
        !     call this%err%werrfl2(class//sname//' started')
        !     call this%lp%del()
        !     call this%a0x%del()
        !     call this%a0y%del()
        !     call this%a0z%del()
        !     call this%a1x%del()
        !     call this%a1y%del()
        !     call this%a1z%del()
        !     call this%a2x%del()
        !     call this%a2y%del()
        !     call this%a2z%del()
        !     call this%err%werrfl2(class//sname//' ended')
                            
        ! end subroutine end_beam3d
              
        ! subroutine pmove_beam3d(this,fd,rtag,stag,sid)
                
        !     implicit none
            
        !     class(beam3d), intent(inout) :: this
        !     class(field3d), intent(in) :: fd
        !     integer, intent(in) :: rtag, stag
        !     integer, intent(inout) :: sid
        
        !     character(len=18), save :: sname = 'pmove:'
            
        !     call this%err%werrfl2(class//sname//' started')    
        !     call this%pd%pmv(fd%getrs(),rtag,stag,sid)
        !     call this%err%werrfl2(class//sname//' ended')
            
        ! end subroutine pmove_beam3d
        
        ! subroutine writehdf5_beam3d(this,file,dspl,delta,rtag,stag,id)
        
        !     implicit none
            
        !     class(beam3d), intent(inout) :: this
        !     class(hdf5file), intent(in) :: file
        !     real, dimension(3), intent(in) :: delta
        !     integer, intent(in) :: dspl, rtag, stag
        !     integer, intent(inout) :: id
        
        !     character(len=18), save :: sname = 'writehdf5_beam3d:'
        
        !     call this%err%werrfl2(class//sname//' started')                  
        !     call this%pd%wr(file,dspl,delta,rtag,stag,id)
        !     call this%err%werrfl2(class//sname//' ended')
         
        ! end subroutine writehdf5_beam3d
                   
        ! subroutine writerst_beam3d(this,file)
        
        !     implicit none
            
        !     class(beam3d), intent(inout) :: this
        !     class(hdf5file), intent(in) :: file
        
        !     character(len=18), save :: sname = 'writerst_beam3d:'
        
        !     call this%err%werrfl2(class//sname//' started')                  
        !     call this%pd%wrst(file)
        !     call this%err%werrfl2(class//sname//' ended')
         
        ! end subroutine writerst_beam3d
        
        subroutine writelp_laser(this,file,rtag,stag,id)
        
            implicit none
            
            class(laser), intent(inout) :: this
            class(hdf5file), intent(in) :: file
            integer, intent(in) :: rtag, stag
            integer, intent(inout) :: id
        
            character(len=18), save :: sname = 'writelp_laser:'
        
            call this%err%werrfl2(class//sname//' started')                  
            call this%lp%wr(file,1,rtag,stag,id)
            call this%err%werrfl2(class//sname//' ended')
         
        end subroutine writelp_laser
                   
        subroutine writelpslice_laser(this,file,slice,spos,rtag,stag,id)
        
            implicit none
            
            class(laser), intent(inout) :: this
            class(hdf5file), intent(in) :: file
            integer, intent(in) :: rtag, stag, slice, spos
            integer, intent(inout) :: id         
        
            character(len=18), save :: sname = 'writeqslice_beam3d:'
        
            call this%err%werrfl2(class//sname//' started')                  
            call this%lp%wr(file,1,slice,spos,rtag,stag,id)
            call this%err%werrfl2(class//sname//' ended')
         
        end subroutine writelpslice_laser
                  
        ! subroutine readrst_beam3d(this,file)
        
        !     implicit none
            
        !     class(beam3d), intent(inout) :: this
        !     class(hdf5file), intent(in) :: file
        
        !     character(len=18), save :: sname = 'readrst_beam3d:'
        
        !     call this%err%werrfl2(class//sname//' started')                  
        !     call this%pd%rrst(file)
        !     call this%err%werrfl2(class//sname//' ended')
         
        ! end subroutine readrst_beam3d
        
    !     class(beam3d), intent(inout) :: this
    !     class(field3d), intent(in) :: fd
    !     integer, intent(in) :: rtag, stag
    !     integer, intent(inout) :: sid
    
    !     character(len=18), save :: sname = 'pmove:'
        
    !     call this%err%werrfl2(class//sname//' started')    
    !     call this%pd%pmv(fd%getrs(),rtag,stag,sid)
    !     call this%err%werrfl2(class//sname//' ended')
        
    ! end subroutine pmove_beam3d
    
    ! subroutine writehdf5_beam3d(this,file,dspl,delta,rtag,stag,id)
    
    !     implicit none
        
    !     class(beam3d), intent(inout) :: this
    !     class(hdf5file), intent(in) :: file
    !     real, dimension(3), intent(in) :: delta
    !     integer, intent(in) :: dspl, rtag, stag
    !     integer, intent(inout) :: id
    
    !     character(len=18), save :: sname = 'writehdf5_beam3d:'
    
    !     call this%err%werrfl2(class//sname//' started')                  
    !     call this%pd%wr(file,dspl,delta,rtag,stag,id)
    !     call this%err%werrfl2(class//sname//' ended')
     
    ! end subroutine writehdf5_beam3d
               
    ! subroutine writerst_beam3d(this,file)
    
    !     implicit none
        
    !     class(beam3d), intent(inout) :: this
    !     class(hdf5file), intent(in) :: file
    
    !     character(len=18), save :: sname = 'writerst_beam3d:'
    
    !     call this%err%werrfl2(class//sname//' started')                  
    !     call this%pd%wrst(file)
    !     call this%err%werrfl2(class//sname//' ended')
     
    ! end subroutine writerst_beam3d
    
    ! subroutine writeq_beam3d(this,file,rtag,stag,id)
    
    !     implicit none
        
    !     class(beam3d), intent(inout) :: this
    !     class(hdf5file), intent(in) :: file
    !     integer, intent(in) :: rtag, stag
    !     integer, intent(inout) :: id
    
    !     character(len=18), save :: sname = 'writeq_beam3d:'
    
    !     call this%err%werrfl2(class//sname//' started')                  
    !     call this%q%wr(file,1,rtag,stag,id)
    !     call this%err%werrfl2(class//sname//' ended')
     
    ! end subroutine writeq_beam3d
               
    ! subroutine writeqslice_beam3d(this,file,slice,spos,rtag,stag,id)
    
    !     implicit none
        
    !     class(beam3d), intent(inout) :: this
    !     class(hdf5file), intent(in) :: file
    !     integer, intent(in) :: rtag, stag, slice, spos
    !     integer, intent(inout) :: id         
    
    !     character(len=18), save :: sname = 'writeqslice_beam3d:'
    
    !     call this%err%werrfl2(class//sname//' started')                  
    !     call this%q%wr(file,1,slice,spos,rtag,stag,id)
    !     call this%err%werrfl2(class//sname//' ended')
     
    ! end subroutine writeqslice_beam3d
              
    ! subroutine readrst_beam3d(this,file)
    
    !     implicit none
        
    !     class(beam3d), intent(inout) :: this
    !     class(hdf5file), intent(in) :: file
    
    !     character(len=18), save :: sname = 'readrst_beam3d:'
    
    !     call this%err%werrfl2(class//sname//' started')                  
    !     call this%pd%rrst(file)
    !     call this%err%werrfl2(class//sname//' ended')
     
    ! end subroutine readrst_beam3d

    end module laser_class