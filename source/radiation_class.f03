! radiation module for QuickPIC Open Source 

        module radiation_class
        
        use input_class

        implicit none

        type radiation
         real :: x0,x1,y0,y1,z,dx,dy,dt
         real :: t0,t1
         real :: w0,w1
         integer::Nx,Ny,Nt,Nw
         real,dimension(:,:,:),pointer:: intensity => null()
         real,dimension(:),pointer:: spectrum => null()
         real,dimension(:,:,:,:),pointer :: E => null()
         real,dimension(:,:,:,:),pointer :: Ew_re => null()
         real,dimension(:,:,:,:),pointer :: Ew_im => null()
         real,dimension(:,:,:,:),pointer :: B => null()
         real,dimension(:,:,:,:),pointer :: Bw_re => null()
         real,dimension(:,:,:,:),pointer :: Bw_im => null()
         
         
         contains
         procedure :: new
         procedure :: calculate_radiation
         procedure :: calculate_intensity
!        procedure :: calculate_spectrum
!        procedure :: calculate_field
!        procedure :: calculate_field_w
!        procedure :: write_radiation
        end type radiation
        
        
        contains
        FUNCTION cross(a, b)
        REAL, DIMENSION(3) :: cross
        real, DIMENSION(3), INTENT(IN) :: a, b

        cross(1) = a(2) * b(3) - a(3) * b(2)
        cross(2) = a(3) * b(1) - a(1) * b(3)
        cross(3) = a(1) * b(2) - a(2) * b(1)
        END FUNCTION cross

        subroutine new(this,input,i)

            implicit none

            class(radiation), intent(inout) :: this
            type(input_json), intent(inout), pointer :: input
            integer, intent(in) :: i

            real :: x0,x1,y0,y1,z
            real :: min, max
            real :: t0,t1
            real :: w0,w1
            integer::Nx,Ny,Nt,Nw,ii
            character(len=20) :: sn,s1,diag

            real,allocatable :: intensity(:,:,:) 
            real,allocatable:: spectrum(:) 
            real,allocatable :: E(:,:,:,:),Ew_re(:,:,:,:),Ew_im(:,:,:,:) 
            real,allocatable :: B(:,:,:,:),Bw_re(:,:,:,:),Bw_im(:,:,:,:) 

            write (sn,'(I3.3)') i
            s1 = 'beam('//trim(sn)//').radiation'

            call input%get(s1//'.detector.x(1)',min)
            call input%get(s1//'.detector.x(2)',max)
            call input%get(s1//'.detector.nx',Nx)
            this%x0 = min
            this%x1 = max
            this%Nx = Nx

            call input%get(s1//'.detector.y(1)',min)
            call input%get(s1//'.detector.y(2)',max)
            call input%get(s1//'.detector.ny',Ny)
            this%y0 = min
            this%y1 = max
            this%Ny = Ny

            call input%get(s1//'.z',z)
            this%z = z

            call input%get(s1//'.time(1)',min )
            call input%get(s1//'.time(1)',max )
            call input%get(s1//'.nt',Nt )
            this%t0 = min
            this%t1 = max
            this%Nt = Nt

            call input%get(s1//'.spectrum.frequence(1)',min )
            call input%get(s1//'.spectrum.frequence(2)',max)
            call input%get(s1//'.spectrum.nw',Nw )
            this%w0 = min
            this%w1 = max
            this%Nw = Nw

            this%dx = (this%x1 - this%x0)/this%Nx
            this%dy = (this%y1 - this%y0)/this%Ny
            this%dt = (this%t1 - this%t0)/this%Nt

            call input%in%info(s1//'diag',n_children=n)

            DO ii in 1:n
                call input%get(s1//'.diag(ii)',diag)
                select case(trim(diag))
                case('intensity')
                    allocate(this%intensity(Nx,Ny,Nt))
                    intensity = 0
                case('spectrum')
                    allocate(this%spectrum(Nw))
                    spectrum(:) = 0
                case('field')
                    allocate(this%E(3,Nx,Ny,Nt))
                    this%E = 0
                    allocate(this%B(3,this%Nx,this%Ny,this%Nt))
                    this%B = 0
                case('field_omega')
                    allocate(this%Ew_re(3,this%Nx,this%Ny,this%Nw))
                    this%Ew_re = 0
                    allocate(this%Ew_im(3,this%Nx,this%Ny,this%Nw))
                    this%Ew_im = 0
                    allocate(this%Bw_re(3,this%Nx,this%Ny,this%Nw)) 
                    this%Bw_re = 0
                    allocate(this%Bw_im(3,this%Nx,this%Ny,this%Nw))
                    this%Bw_im = 0
                end select

            
        end subroutine new

        subroutine calculate_radiation(this,part,s_beam,ds)
            implicit none
            if allocated(this%intensity)
                call calculate_intensity(this,part,Np,s_beam,ds)
            end if
            if allocated(this%spectrum)
                call calculate_spectrum(this,part,Np,s_beam,ds)
            end if
        end subroutine calculate_radiation

        subroutine calculate_intensity(this,part,Np,s_beam,ds)
            implicit none

            class(radiation), intent(inout) :: this
            integer :: n,i,j,k,Np
            integer :: t_n
            real :: x,y,z,tempz
            real :: t,t_obs,m0,m1
            real :: s_beam,ds,delta_t
            real :: dx,dy,dt
            real ::r
            real :: intensity
            real,dimension(10,Np) :: part
            real,dimension(3) :: posi,vel,acc,posi_x,vec_r,vec_n,p,temp1,temp2
            real :: c
            z = s_beam
            DO n = 1,Np
                x = part(1,n)
                y = part(2,n)
                tempz = part(3,n)
                t = s_beam + tempz
                posi = (/x,y,z/)
                p = part(4:6,n)
                vel = p/(p(1)**2+p(2)**2+p(3)**2 + 1)**0.5
                acc = part(8:10,n)
                delta_t = ds/vel(3)
                DO i=1,this%Ny,1
                    DO j = 1,this%Nx,1
                        posi_x = (/this%x0+(j-1)*this%dx,this%y1-(i-1)*this%dy,z/)
                        vec_r = posi_x - posi
                        r = NORM2(vec_r)
                        vec_n = vec_r/r
                        t_obs = (t + r/c-this%t0)/dt
                        if (t_obs > 0 .AND. t_obs < this%Nt) then                            
                            t_n = t_obs
                            m0 = t_obs - t_n
                            m1 = 1 - m0
                            temp1 = cross(vec_n-vel,acc)
                            temp2 = cross(vec_n,temp1)
                            intensity =((temp2(1)**2+temp2(2)**2+temp2(3)**2)/(1 - vec_n(1)*vel(1)-vec_n(2)*vel(2)-vec_n(2)*vel(2)))**5 )*delta_t
                            this%intensity(i,j,t_n) = this%intensity(i,j,t_n) + intensity*m1
                            this%intensity(i,j,t_n + 1) = this%intensity(i,j,t_n+1) + intensity*m0
                        end if
                    end DO
                end DO
            end DO
        end subroutine calculate_intensity
!
!        subroutine calculate_spectrum(this,part,t_beam)
!
!        end subroutine calculate_spectrum
!
!        subroutine calculate_field(this,part,t_beam)
!
!        end subroutine
!
!        subroutine calculate_field_w(this,part,t_beam)
!        
!        end subroutine


        end module radiation_class