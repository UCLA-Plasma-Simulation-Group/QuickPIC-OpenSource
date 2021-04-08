! radiation module for QuickPIC Open Source 

        module radiation_class
        
        implicit none

        type radiation
         real :: x0,x1,y0,y1,z
         real :: t0,t1
         real :: w0,w1
         integer::Nx,Ny,Nt,Nw
         real,dimension(:,:,:),pointer:: intensity => null()
         real,dimension(:),pointer:: spectrum => null()
         real,dimension(3,:,:,:),pointer :: E => null()
         real,dimension(3,:,:,:),pointer :: Ew_re => null()
         real,dimension(3,:,:,:),pointer :: Ew_im => null()
         real,dimension(3,:,:,:),pointer :: B => null()
         real,dimension(3,:,:,:),pointer :: Bw_re => null()
         real,dimension(3,:,:,:),pointer :: Bw_im => null()
         
         
         contains
         procedure :: new
         procedure :: calculate_intensity
         procedure :: calculate_spectrum
         procedure :: calculate_field
         procedure :: calculate_field_w
         procedure :: write_radiation
        end type radiation
        
        
        contains
        subroutine new(this,input，i)

            implicit none

            class(radiation), intent(inout) :: this
            type(input_json), intent(inout), pointer :: input
            integer, intent(in) :: i

            real :: x0,x1,y0,y1,z
            real :: min, max,
            real :: t0,t1
            real :: w0,w1
            integer::Nx,Ny,Nt,Nw
            character(len=20) :: sn,s1

            real,allocatable :: intensity(:,:,:,:) 
            real,allocatable:: spectrum(:) 
            real,allocatable :: E(:,:,:,:),Ew_re(:,:,:,:),Ew_im(:,:,:,:) 
            real,allocatable :: B(:,:,:,:),Bw_re(:,:,:,:),Bw_im(:,:,:,:) 

            write (sn,'(I3.3)') i
            s1 = 'detector('//trim(sn)//')'

            call input%get(s1//'.x(1)',min)
            call input%get(s1//'.x(2)',max)
            call input%get(s1//'.nx',Nx)
            this%x0 = min
            this%x1 = max
            this%Nx = Nx

            call input%get(s1//'.y(1)',min)
            call input%get(s1//'.y(2)',max)
            call input%get(s1//'.ny',Ny)
            this%y0 = min
            this%y1 = max
            this%yx = Nx

            call input%get(s1//'.z',z)
            this%z = z

            call input%get('radiation.time(1)',min )
            call input%get('radiation.time(1)',max )
            call input%get('radiation.nt',Nt )
            this%t0 = min
            this%t1 = max
            this%Nt = Nt

            call input%get('spectrum.frequence(1)',min )
            call input%get('spectrum.frequence(2)',max)
            call input%get('spectrum.nw',Nw )
            this%w0 = min
            this%w1 = max
            this%Nw = Nw
            
            allocate(intensity(3,Nx,Ny,Nt))
            intensity = 0
            this%intensity = intensity

            allocate(spectrum(Nw))
            spectrum(:) = 0
            this%spectrum = spectrum

            allocate(E(3,Nx,Ny,Nt))
            E = 0
            this%E = E

            allocate(B(3,this%Nx,this%Ny,this%Nt))
            B = 0
            this%B = B
            
            allocate(Ew_re(3,this%Nx,this%Ny,this%Nw))
            Ew_re = 0
            this%Ew_re =Ew_re

            real,dimension(Ew_im(3,this%Nx,this%Ny,this%Nw))
            Ew_im = 0
            this%B = Ew_im
            real,dimension(Bw_re(3,this%Nx,this%Ny,this%Nw)) 
            Bw_re = 0
            this%Bw_re = Bw_re
            real,dimension(Bw_im(3,this%Nx,this%Ny,this%Nw))
            Bw_im = 0
            this%Bw_im = Bw_im
            
        end subroutine new

        subroutine calculate_energy(this,part,t_beam)

        end subroutine calculate_energy

        subroutine calculate_spectrum(this,part,t_beam)

        end subroutine calculate_spectrum

        subroutine calculate_field(this,part,t_beam)

        end subroutine

        subroutine calculate_field_w(this,part,t_beam)
        
        end subroutine