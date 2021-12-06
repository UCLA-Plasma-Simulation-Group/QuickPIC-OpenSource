! radiation module for QuickPIC Open Source 

        module radiation_class
        
        use input_class
        use perrors_class
        use parallel_pipe_class
        use hdf5io_class

        implicit none

        private

        public :: radiation

        type radiation

         private

         class(perrors), pointer, public :: err => null()
!         class(spect3d), pointer, public :: sp => null()
         class(parallel_pipe), pointer, public :: p => null()

         real :: x0,x1,y0,y1,z,dx,dy,dt,dw,n0,dx1,dy1
         real :: t0,t1
         real :: w0,w1
         real, dimension(3) :: delta
         integer :: Nx,Ny,Nt,Nw,in_n,Nx1,Ny1,type
         
         real,allocatable :: intensity(:,:,:) 
         real,allocatable :: spectrum(:,:,:) 
         real,allocatable :: field_E(:,:,:,:),field_B(:,:,:,:)
         real,allocatable :: Ew_re(:,:,:,:),Ew_im(:,:,:,:) 
         
         contains
         generic :: new => init_radiation
         generic :: del => end_radiation         
         generic :: dp => calculate_radiation
         generic :: wrr => write_radiation
         procedure, private :: init_radiation, end_radiation
         procedure, private :: calculate_radiation
         procedure, private :: calculate_intensity
         procedure, private :: write_radiation
         procedure, private :: calculate_spectrum
         procedure, private :: calculate_field
         procedure, private :: calculate_field_w
         procedure, private :: calculate_field_w2
         procedure, private :: calculate_radiation2
         procedure, private :: calculate_intensity2
!         procedure :: calculate_spectrum2
!         procedure :: calculate_field2
!         
        end type radiation
        character(len=10) :: class = "radiation"
        character(len=128) :: erstr 
        
        contains
!
        FUNCTION cross(a, b)

        implicit none
        
        REAL, DIMENSION(3) :: cross
        real, DIMENSION(3), INTENT(IN) :: a, b

        cross(1) = a(2) * b(3) - a(3) * b(2)
        cross(2) = a(3) * b(1) - a(1) * b(3)
        cross(3) = a(1) * b(2) - a(2) * b(1)
        END FUNCTION cross
!
        subroutine init_radiation(this,pp,perr,input,i)

            implicit none

            class(radiation), intent(inout) :: this
            type(input_json), intent(inout), pointer :: input
            class(perrors), intent(in), pointer :: perr
            class(parallel_pipe), intent(in), pointer :: pp
            integer, intent(in) :: i
!
            real :: x0,x1,y0,y1,z,beam_z
            real :: min, max
            real :: t0,t1
            real :: w0,w1
            real(8) :: n0
            integer::Nx,Ny,Nx1,Ny1,Nt,Nw,ii,n_children,n,in_n         
            integer :: indx, indy, indz,type
            character(len=20) :: sn,s1
            character(len=:), allocatable :: diag
            character(len=25), save :: sname = 'init_radiation:'

            this%err => perr
            this%p => pp

            call this%err%werrfl2(class//sname//' started')

            call input%get('simulation.indx',indx)
            call input%get('simulation.indy',indy)
            call input%get('simulation.indz',indz)
            call input%get('simulation.n0',n0)
            this%n0 = n0
            write (sn,'(I3.3)') i
            s1 = 'beam('//trim(sn)//').radiation'

            call input%get(s1//'.detector.x(1)',min)
            call input%get(s1//'.detector.x(2)',max)
            call input%get(s1//'.detector.nx',Nx)
            this%x0 = min
            this%x1 = max
            this%Nx = Nx
            call input%get('simulation.box.x(1)',min)
            call input%get('simulation.box.x(2)',max)
            this%delta(1) = (max-min) / real(2**indx)
            this%x0 = this%x0 - min
            this%x1 = this%x1 - min

            call input%get(s1//'.detector.y(1)',min)
            call input%get(s1//'.detector.y(2)',max)
            call input%get(s1//'.detector.ny',Ny)
            this%y0 = min
            this%y1 = max
            this%Ny = Ny
            call input%get('simulation.box.y(1)',min)
            call input%get('simulation.box.y(2)',max)
            this%delta(2) = (max-min) / real(2**indy)
            this%y0 = this%y0 - min
            this%y1 = this%y1 - min

            call input%get(s1//'.detector.z',z)
            call input%get(s1//'.detector.type',type)
            call input%get('simulation.box.z(1)',min)
            call input%get('simulation.box.z(2)',max)
            call input%get('beam('//trim(sn)//').center(3)',beam_z)
            this%delta(3) = (max-min)/real(2**indz)
            this%type = type
            if (this%type==1) then
                this%z = z
            else if (this%type==2) then
                this%z = z - beam_z
            end if 
            call input%get(s1//'.time(1)',min )
            call input%get(s1//'.time(2)',max )
            call input%get(s1//'.nt',Nt )
            this%t0 = min
            this%t1 = max
            this%Nt = Nt

            this%dx = (this%x1 - this%x0)/this%Nx
            this%dy = (this%y1 - this%y0)/this%Ny
            this%dt = (this%t1 - this%t0)/this%Nt



            call input%info(s1//'.diag',n_children=n)
            
            DO ii = 1,n
                write (sn,'(I3.3)') ii
                call input%get(s1//'.diag('//trim(sn)//')',diag)
                select case(trim(diag))
                case('intensity')
                    allocate(this%intensity(this%Nx,this%Ny,this%Nt+1))
                    this%intensity = 0
                case('spectrum')
                    call input%get(s1//'.spectrum.frequence(1)',min )
                    call input%get(s1//'.spectrum.frequence(2)',max)
                    call input%get(s1//'.spectrum.nw',Nw )
                    call input%get(s1//'.spectrum.Nx',Nx1 )
                    call input%get(s1//'.spectrum.Ny',Ny1 )
                    call input%get(s1//'.spectrum.in_n',in_n )
                    this%in_n = in_n
                    this%w0 = min
                    this%w1 = max
                    this%Nw = Nw
                    this%Nx1 = Nx1
                    this%Ny1 = Ny1
                    this%dx1 = (this%x1 - this%x0)/this%Nx1
                    this%dy1 = (this%y1 - this%y0)/this%Ny1
                    this%dw = (this%w1 - this%w0)/this%Nw
                    allocate(this%Ew_re(3,this%Nx1,this%Ny1,this%Nw+1))
                    this%Ew_re = 0
                    allocate(this%Ew_im(3,this%Nx1,this%Ny1,this%Nw+1))
                    this%Ew_im = 0
                    allocate(this%spectrum(this%Nx1,this%Ny1,this%Nw+1))
                    this%spectrum = 0
                case('field')
                    allocate(this%field_E(3,this%Nx,this%Ny,this%Nt+1))
                    this%field_E = 0
                    allocate(this%field_B(3,this%Nx,this%Ny,this%Nt+1))
                    this%field_B = 0
                case('field_w')
                    allocate(this%Ew_re(3,this%Nx,this%Ny,this%Nw+1))
                    this%Ew_re = 0
                    allocate(this%Ew_im(3,this%Nx,this%Ny,this%Nw+1))
                    this%Ew_im = 0
                end select
            end DO

            call this%err%werrfl2(class//sname//' ended')
        end subroutine init_radiation
!
        subroutine end_radiation(this)

            implicit none

            class(radiation), intent(inout) :: this

            return
        
        end subroutine end_radiation
!
        subroutine calculate_radiation(this,part,s_beam,ds,Np)

            implicit none

            class(radiation), intent(inout) :: this
            real, intent(inout) :: s_beam,ds
            real,dimension(:,:),pointer,intent(inout) :: part
            integer,intent(inout) :: Np

            character(len=18), save :: sname = 'calculate_radiation:'
            call this%err%werrfl2(class//sname//' started')
            if (allocated(this%intensity)) then
                call calculate_intensity(this,part,Np,s_beam,ds)
            end if
            if (allocated(this%spectrum)) then
                if (this%in_n == 0 )then
                    call calculate_field_w(this,part,Np,s_beam,ds)
                else 
                    call calculate_field_w2(this,part,Np,s_beam,ds)
                end if 
                call calculate_spectrum(this)
            end if 
            if (allocated(this%field_E)) then 
                call calculate_field(this,part,Np,s_beam,ds)
            end if
            call this%err%werrfl2(class//sname//' ended')
        end subroutine calculate_radiation
!
        subroutine calculate_intensity(this,part,Np,s_beam,ds)
            implicit none

            class(radiation), intent(inout) :: this
            real, intent(inout) :: s_beam,ds
            integer, intent(inout) :: Np
            real,dimension(:,:), pointer, intent(inout) :: part

            integer :: n,i,j,k,tt3
            integer :: t_n,n1,n2
            real :: x,y,z,tempz
            real :: t,t_obs,m0,m1
            real :: delta_t
            real :: deltax, deltay, deltaz
            real :: dx,dy,dt
            real :: r
            real :: intensity
            real(8) :: A,n0
            real(8),dimension(3) :: posi,vel,acc,posi_x,vec_r,vec_n,p,temp0,temp1,temp2,pnew,pold,&
                                 &vnew,vold,force
            real :: temp3,temp4
            real :: c
            character(len=25), save :: sname = 'calculate_intensity:'
            call this%err%werrfl2(class//sname//' started') 
            deltax = this%delta(1)
            deltay = this%delta(2)
            deltaz = this%delta(3)
            A = 3.4548e-33*sqrt(this%n0)
            z = s_beam
            tt3 = size(part,1)
            if (this%type == 1) then
                n1 = 1 
                n2 = 0
            else if (this%type == 2) then
                n1 = 0
                n2 = 1
            end if 
            DO n = 1,Np
               if (part(tt3,n) > 0.5) then
                    x = part(1,n)*deltax
                    y = part(2,n)*deltay
                    tempz = part(3,n)*deltaz
                    t = s_beam + tempz
!                    posi = (/x,y,n1*z-n2*tempz/)
                    posi = (/x,y,z/)
                    p = part(4:6,n)
                    vel = p/(p(1)**2+p(2)**2+p(3)**2 + 1)**0.5
                    force = part(8:10,n)
                    pold = p - force
                    pnew = p + force
                    vold = pold/(pold(1)**2+pold(2)**2+pold(3)**2 + 1)**0.5
                    vnew = pnew/(pnew(1)**2+pnew(2)**2+pnew(3)**2 + 1)**0.5
                    delta_t = ds/vel(3)
                    acc = (vnew-vold)/delta_t
                    DO j = 1,this%Ny,1
                        DO i = 1,this%Nx,1 
                            posi_x = (/this%x0+(i-1)*this%dx,this%y1-(j-1)*this%dy,this%z/)
                            vec_r = posi_x - n1*posi
                            r = (vec_r(1)**2 +vec_r(2)**2 +vec_r(3)**2)**0.5
                            vec_n = vec_r/r
!                            t_obs = (t + r - this%t0)/this%dt
                            t_obs = (n1*t + r - this%t0)/this%dt
                            if (t_obs > 0 .AND. t_obs < this%Nt) then                            
                                t_n = t_obs
                                m0 = t_obs - t_n
                                m1 = 1.0 - m0
                                t_n = t_n + 1
                                temp0 = vec_n-vel
    !                          temp1(1) = temp0(2)*acc(3) - temp0(3)*acc(2)
    !                          temp1(2) = temp0(3)*acc(1) - temp0(1)*acc(3)
    !                          temp1(3) = temp0(1)*acc(2) - temp0(2)*acc(1)
                                temp1 = cross(temp0,acc)
                                temp2 = cross(vec_n,temp1)
                                temp3 = (temp2(1)**2+temp2(2)**2+temp2(3)**2)
                                temp4 = (1 - vec_n(1)*vel(1)-vec_n(2)*vel(2)-vec_n(3)*vel(3))
                                intensity = A*(temp3/temp4**5)*delta_t
                                this%intensity(i,j,t_n) = this%intensity(i,j,t_n) + intensity*m1
                                this%intensity(i,j,t_n+1) = this%intensity(i,j,t_n+1) + intensity*m0
                            else 
                                print*,r
                            end if
                        end DO
                    end DO
                end if
            end DO
            call this%err%werrfl2(class//sname//' ended')
        end subroutine calculate_intensity
!
!        subroutine calculate_spectrum(this,part,t_beam)
!
!        end subroutine calculate_spectrum
!
        subroutine calculate_field(this,part,Np,s_beam,ds)
            implicit none

            class(radiation), intent(inout) :: this
            real,dimension(:,:),intent(inout) :: part
            integer,intent(in) :: Np
            real,intent(in) :: s_beam,ds

            integer :: n,i,j,k,tt3,n1,n2
            integer :: t_n
            real :: x,y,z,tempz
            real :: t,t_obs,m0,m1
            real :: delta_t
            real :: dx,dy,dt
            real :: r 
            real :: deltax, deltay, deltaz
            real,dimension(3) :: field_E,field_B,S
            real(8) :: A,n0
            real,dimension(3) :: posi,vel,acc,posi_x,vec_r,vec_n,p,temp1,temp2,pnew,pold,vnew,vold,force
            real :: temp3,temp4
            character(len=25), save :: sname = 'calculate_field:'
            call this%err%werrfl2(class//sname//' started')
            z = s_beam
            deltax = this%delta(1)
            deltay = this%delta(2)
            deltaz = this%delta(3)
            tt3 = size(part,1)
            A = 1.528671463266763e-08
            if (this%type == 1) then
                n1 = 1 
                n2 = 0
            else if (this%type == 2) then
                n1 = 0
                n2 = 1
            end if 
            DO n = 1,Np
            if (part(tt3,n) > 0.5) then
                x = part(1,n)*deltax
                y = part(2,n)*deltay
                tempz = part(3,n)*deltaz
                t = s_beam + tempz
                posi = (/x,y,z/)
!                posi = (/x,y,-tempz/)
                p = part(4:6,n)
                vel = p/(p(1)**2+p(2)**2+p(3)**2 + 1)**0.5
                force = part(8:10,n)
                pold = p - force
                pnew = p + force
                vold = pold/(pold(1)**2+pold(2)**2+pold(3)**2 + 1)**0.5
                vnew = pnew/(pnew(1)**2+pnew(2)**2+pnew(3)**2 + 1)**0.5
                delta_t = ds/vel(3)
                acc = (vnew-vold)/delta_t
                DO i=1,this%Nx,1
                    DO j = 1,this%Ny,1
                        posi_x = (/this%x0+(i-1)*this%dx,this%y1-(j-1)*this%dy,this%z/)
                        vec_r = posi_x - n1*posi
                        r = (vec_r(1)**2 +vec_r(2)**2 +vec_r(3)**2)**0.5
                        vec_n = vec_r/r
                        !t_obs = (t + r-this%t0)/this%dt
                        t_obs = (n1*t + r - this%t0)/this%dt
                        if (t_obs > 0 .AND. t_obs < this%Nt) then                            
                            t_n = t_obs
                            m0 = t_obs - t_n
                            m1 = 1.0 - m0
                            t_n = t_n + 1
                            temp1 = cross(vec_n-vel,acc)
                            temp2 = cross(vec_n,temp1)
                            temp4 = (1 - vec_n(1)*vel(1)-vec_n(2)*vel(2)-vec_n(3)*vel(3))
                            field_E = (temp2/(temp4**2.5)/(n1*r+n2*this%z))
                            field_B = cross(vec_n,field_E)
                            this%field_E(:,i,j,t_n) = this%field_E(:,i,j,t_n) + field_E*m1
                            this%field_E(:,i,j,t_n + 1) = this%field_E(:,i,j,t_n+1) + field_E*m0
                            this%field_B(:,i,j,t_n) = this%field_B(:,i,j,t_n) + field_B*m1
                            this%field_B(:,i,j,t_n + 1) = this%field_B(:,i,j,t_n+1) + field_B*m0
                        else 
                            print*,r
                        end if
                    end DO
                end DO
            end if
            end DO    
            call this%err%werrfl2(class//sname//' ended')      
        end subroutine calculate_field
!

        subroutine calculate_field_w(this,part,Np,s_beam,ds)
            implicit none

            class(radiation), intent(inout) :: this
            real,dimension(:,:),intent(inout) :: part
            integer,intent(in) :: Np
            real,intent(in) :: s_beam,ds

            real :: deltax, deltay, deltaz
            integer :: n,i,j,k,tt3,n1,n2
            integer :: t_n
            real :: x,y,z,tempz
            real :: t,t_obs
            real :: delta_t
            real :: dx,dy,dt
            real :: r 
            real :: w
            
            real,dimension(3) :: Ew_im,Ew_re
            real,dimension(3) :: posi,vel,acc,posi_x,vec_r,vec_n,p,temp1,temp2,pnew,pold,vnew,vold,force
            real :: temp3,temp4,temp5_im,temp5_Re
            real :: c
            character(len=25), save :: sname = 'calculate_field_w:'
            call this%err%werrfl2(class//sname//' started')
            z = s_beam
            deltax = this%delta(1)
            deltay = this%delta(2)
            deltaz = this%delta(3)
            tt3 = size(part,1)   
            if (this%type == 1) then
                n1 = 1 
                n2 = 0
            else if (this%type == 2) then
                n1 = 0
                n2 = 1
            end if 
            DO n = 1,Np,1
             if (part(tt3,n) > 1.5) then
                x = part(1,n)*deltax
                y = part(2,n)*deltay
                tempz = part(3,n)*deltaz
                t = s_beam + tempz
                posi = (/x,y,z/)
                p = part(4:6,n)
                vel = p/(p(1)**2+p(2)**2+p(3)**2 + 1)**0.5
                force = part(8:10,n)
                pold = p - force
                pnew = p + force
                vold = pold/(pold(1)**2+pold(2)**2+pold(3)**2 + 1)**0.5
                vnew = pnew/(pnew(1)**2+pnew(2)**2+pnew(3)**2 + 1)**0.5
                delta_t = ds/vel(3)
                acc = (vnew-vold)/delta_t
                DO i=1,this%Nx1,1
                    DO j = 1,this%Ny1,1
                        posi_x = (/this%x0+(i-0.5)*this%dx1,this%y1-(j-0.5)*this%dy1,this%z/)
                        vec_r = posi_x - n1*posi
                        r = (vec_r(1)**2 +vec_r(2)**2 +vec_r(3)**2)**0.5
                        vec_n = vec_r/r
                        t_obs = t + r
                        temp1 = cross(vec_n-vel,acc)
                        temp2 = cross(vec_n,temp1)
                        temp4 = (1 - vec_n(1)*vel(1)-vec_n(2)*vel(2)-vec_n(3)*vel(3))
                        Do k = 1,this%Nw,1
                            w = this%w0 + (k-1)*this%dw
                            temp5_Re = cos(w*(t_obs*n1+n2*(t-(vec_n(1)*posi(1)+vec_n(2)*posi(2)+vec_n(3)*posi(3)))))
                            temp5_im = sin(w*(t_obs*n1+n2*(t-(vec_n(1)*posi(1)+vec_n(2)*posi(2)+vec_n(3)*posi(3)))))
                            Ew_re = (temp2*temp5_re/temp4**2)*delta_t                        
                            Ew_im =(temp2*temp5_im/temp4**2)*delta_t
                            this%Ew_re(:,i,j,k) = this%Ew_re(:,i,j,k) + Ew_re
                            this%Ew_im(:,i,j,k) = this%Ew_im(:,i,j,k) + Ew_im
                        end DO
                    end DO
                end DO
             end if
            end DO    
            call this%err%werrfl2(class//sname//' ended')      
        end subroutine calculate_field_w
!        
        subroutine calculate_field_w2(this,part,Np,s_beam,ds)
            implicit none

            class(radiation), intent(inout) :: this
            real,dimension(:,:),intent(inout) :: part
            integer,intent(in) :: Np
            real,intent(in) :: s_beam,ds

            real :: deltax, deltay, deltaz
            integer :: n,i,j,k,in_n,m,tt3,n1,n2
            integer :: t_n
            real :: x,y,z,tempz
            real :: t,t_obs,t_old
            real :: delta_t,in_t
            real :: dx,dy,dt
            real :: r 
            real :: w
            real,dimension(3) :: Ew_im,Ew_re
            real,dimension(3) :: posi,posi0,vel,acc,posi_x,vec_r,vec_n,p,pnew,pold,vnew,vold,force,posi_new,vel_new
            real,dimension(3) ::temp0,temp1,temp2
            real :: temp3,temp4,temp5_im,temp5_Re
            real :: c
            character(len=25), save :: sname = 'calculate_field_w2:'
            call this%err%werrfl2(class//sname//' started')
            z = s_beam
            deltax = this%delta(1)
            deltay = this%delta(2)
            deltaz = this%delta(3)
            tt3 = size(part,1)   
            if (this%type == 1) then
                n1 = 1 
                n2 = 0
            else if (this%type == 2) then
                n1 = 0
                n2 = 1
            end if 
            DO n = 1,Np,1
            if (part(tt3,n) > 1.5) then
                x = part(1,n)*deltax
                y = part(2,n)*deltay
                tempz = part(3,n)*deltaz
                posi = (/x,y,z/)
                p = part(4:6,n)
                vel = p/(p(1)**2+p(2)**2+p(3)**2 + 1)**0.5
                force = part(8:10,n)
                pold = p - force
                pnew = p + force
                vold = pold/(pold(1)**2+pold(2)**2+pold(3)**2 + 1)**0.5
                vnew = pnew/(pnew(1)**2+pnew(2)**2+pnew(3)**2 + 1)**0.5
                delta_t = ds/vel(3)
                t_old = s_beam + tempz -0.5*delta_t
                acc = (vnew-vold)/delta_t
                posi0 = posi - vold*0.5*delta_t-0.5*acc*(0.5*delta_t)**2
                in_t =  delta_t/this%in_n
                Do m = 0,this%in_n-1,1 
                    posi_new = posi0 + vold*m*in_t+0.5*acc*(m*in_t)**2
                    vel_new = vold + acc*m*in_t
                    t = t_old + m*in_t
                    DO i=1,this%Nx1,1
                        DO j = 1,this%Ny1,1
                            posi_x = (/this%x0+(i-0.5)*this%dx1,this%y1-(j-0.5)*this%dy1,this%z/)
                            vec_r = posi_x - n1*posi_new
                            r = (vec_r(1)**2 +vec_r(2)**2 +vec_r(3)**2)**0.5
                            vec_n = vec_r/r
                            t_obs = t + r
                            temp0 = vec_n-vel_new
                            temp1 = cross(temp0,acc)
                            temp2 = cross(vec_n,temp1)
                            temp4 = (1 - vec_n(1)*vel_new(1)-vec_n(2)*vel_new(2)-vec_n(3)*vel_new(3))
                            
                            Do k = 1,this%Nw+1,1
                                w= this%w0 + (k-1)*this%dw
                            !    w = 10**(this%w0 + (k-1)*this%dw)
                                temp5_Re = cos(w*(t_obs*n1+n2*(t-(vec_n(1)*posi_new(1)+vec_n(2)*posi_new(2)+vec_n(3)*posi_new(3)))))
                                temp5_im = sin(w*(t_obs*n1+n2*(t-(vec_n(1)*posi_new(1)+vec_n(2)*posi_new(2)+vec_n(3)*posi_new(3)))))
                                Ew_re = (temp2*temp5_re/temp4**2)*in_t
                                Ew_im =(temp2*temp5_im/temp4**2)*in_t
                                this%Ew_re(:,i,j,k) = this%Ew_re(:,i,j,k) + Ew_re
                                this%Ew_im(:,i,j,k) = this%Ew_im(:,i,j,k) + Ew_im
                            end DO
                        end DO
                    end DO
                end DO
            end if 
            end DO 
            call this%err%werrfl2(class//sname//' ended')      
        end subroutine calculate_field_w2
!



        subroutine calculate_radiation2(this,input,part,Np,s_beam,ds)

            implicit none

            class(radiation), intent(inout) :: this
            type(input_json), intent(inout), pointer :: input
            real,intent(in) :: s_beam,ds
            real,dimension(:,:),intent(inout) :: part
            integer,intent(in) :: Np


            real,dimension(3,Np,this%Nx,this%Ny) :: temp_up
            real,dimension(Np,this%Nx,this%Ny) :: temp_down,temp_r,temp_t
            character(len=25), save :: sname = 'calculate_radiation2:'
            call this%err%werrfl2(class//sname//' started')
            call calculate_temp(this,temp_up,&
            temp_down,temp_r,temp_t,part,Np,s_beam,ds) 

            if (allocated(this%intensity)) then
                call calculate_intensity2(this,temp_up,&
                temp_down,temp_r,temp_t,Np,s_beam,ds)
            end if
!           if (allocated(this%spectrum)) then
!               call calculate_field_w2(this,temp_up,&
!               temp_down,temp_r,temp_t,Np,s_beam,ds)
!           end if
!           if (allocated(this%Ew_re)) then 
!               call calculate_field2(this,temp_up,&
!               temp_down,temp_r,temp_t,Np,s_beam,ds)
!           end if
!           call this%err%werrfl2(class//sname//' ended')
        end subroutine calculate_radiation2
!
        subroutine calculate_temp(this,temp_up,temp_down,temp_r,temp_t,part,Np,s_beam,ds)
            implicit none

            class(radiation), intent(inout) :: this
            real,dimension(:,:),intent(inout) :: part
            integer,intent(in) :: Np
            real,intent(in) :: s_beam,ds
            real,dimension(3,Np,this%Nx,this%Ny),intent(inout) :: temp_up
            real,dimension(Np,this%Nx,this%Ny),intent(inout) :: temp_down,temp_r,temp_t
            real :: deltax, deltay, deltaz
            integer :: n,i,j,k
            integer :: t_n
            real :: x,y,z,tempz
            real :: t,t_obs,m0,m1
            real :: delta_t
            real :: dx,dy,dt
            real :: r
            real :: intensity

            real,dimension(3) :: posi,vel,acc,posi_x,vec_r,vec_n,p,temp0,temp1,temp2,pnew,pold,vnew,vold,force
            real :: temp3,temp4
            real :: c
            character(len=25), save :: sname = 'calculate_temp:'
            call this%err%werrfl2(class//sname//' started')
            z = s_beam
            DO n = 1,Np
                x = part(1,n)
                y = part(2,n)
                tempz = part(3,n)
                t = s_beam + tempz
                posi = (/x,y,z/)
                p = part(4:6,n)
                vel = p/(p(1)**2+p(2)**2+p(3)**2 + 1)**0.5
                force = part(8:10,n)
                pold = p - acc
                pnew = p + acc
                vold = pold/(pold(1)**2+pold(2)**2+pold(3)**2 + 1)**0.5
                vnew = pnew/(pnew(1)**2+pnew(2)**2+pnew(3)**2 + 1)**0.5
                delta_t = ds/vel(3)
                acc = (vnew-vold)/delta_t
                DO i=1,this%Nx,1
                    DO j = 1,this%Ny,1
                        posi_x = (/this%x0+(i-1)*this%dx,this%y1-(j-1)*this%dy,z/)
                        vec_r = posi_x - posi
                        r = (vec_r(1)**2 +vec_r(2)**2 +vec_r(3)**2)**0.5
                        vec_n = vec_r/r
                        t_obs = (t + r)
                        temp0 = vec_n-vel
                        temp1 = cross(temp0,acc)
                        temp2 = cross(vec_n,temp1)
                        temp4 = (1 - vec_n(1)*vel(1)-vec_n(2)*vel(2)-vec_n(3)*vel(3))
                        temp_up(:,n,i,j) = temp2
                        temp_down(n,i,j) = temp4
                        temp_r(n,i,j) = r
                        temp_t(n,i,j) = t_obs
                    end DO
                end DO
            end DO      
            call this%err%werrfl2(class//sname//' ended')      
        end subroutine calculate_temp
!
        subroutine calculate_intensity2(this,temp_up,temp_down,temp_r,temp_t,Np,s_beam,ds)
            implicit none

            class(radiation), intent(inout) :: this
            integer,intent(in) :: Np
            real,intent(in) :: s_beam,ds
            real,dimension(3,Np,this%Nx,this%Ny),intent(in) :: temp_up
            real,dimension(Np,this%Nx,this%Ny),intent(in) :: temp_down,temp_r,temp_t

            integer :: n,i,j
            integer :: t_n
            real :: temp
            real :: t_obs,m0,m1
            real :: dx,dy,dt
            real :: r
            real :: intensity
            character(len=25), save :: sname = 'calculate_intensity2:'
            call this%err%werrfl2(class//sname//' started')
            DO n = 1,Np,1
                DO i=1,this%Nx,1
                    DO j = 1,this%Ny,1
                        t_obs = (temp_t(n,i,j) - this%t0)/this%dt
                        if (t_obs > 0 .AND. t_obs < this%Nt) then                            
                            t_n = t_obs
                            m0 = t_obs - t_n
                            m1 = 1 - m0
                            t_n = t_n + 1
                            temp = (temp_up(1,n,i,j)**2+temp_up(2,n,i,j)**2+temp_up(3,n,i,j)**2)
                            intensity =(temp/temp_down(n,i,j)**5)
                            this%intensity(i,j,t_n) = this%intensity(i,j,t_n) + intensity*m1
                            this%intensity(i,j,t_n+1) = this%intensity(i,j,t_n+1) + intensity*m0
                        end if
                    end DO
                end DO
            end DO
            call this%err%werrfl2(class//sname//' ended')
        end subroutine calculate_intensity2
!
!      subroutine calculate_field2(this,temp_up,temp_down,temp_r,temp_t,Np,s_beam,ds)     
!          implicit none
!
!          class(radiation), intent(inout) :: this
!          integer,intent(in) :: Np
!          real,intent(in) :: s_beam,ds
!
!
!          integer :: n,i,j,k
!          integer :: t_n
!          real :: x,y,z,tempz
!          real :: t,t_obs,m0,m1
!          real :: delta_t
!          real :: dx,dy,dt
!          real :: r 
!          
!          real,dimension(3) :: field_E,field_B
!          real,dimension(3) :: posi,vel,acc,posi_x,vec_r,vec_n,p,temp1,temp2,pnew,pold,vnew,vold,force
!          real :: temp3,temp4
!          real :: c
!
!      end subroutine calculate_field2
!
!      subroutine calculate_field_w2(this,temp_up,temp_down,temp_r,temp_t,Np,s_beam,ds)     
!          implicit none
!
!          class(radiation), intent(inout) :: this
!          integer,intent(in) :: Np
!          real,intent(in) :: s_beam,ds
!
!
!          integer :: n,i,j,k
!          integer :: t_n
!          real :: x,y,z,tempz
!          real :: t,t_obs,m0,m1
!          real :: delta_t
!          real :: dx,dy,dt
!          real :: r 
!          
!          real,dimension(3) :: field_E,field_B
!          real,dimension(3) :: posi,vel,acc,posi_x,vec_r,vec_n,p,temp1,temp2,pnew,pold,vnew,vold,force
!          real :: temp3,temp4
!          real :: c
!
!      end subroutine calculate_field_w2
!
        subroutine calculate_spectrum(this)
            implicit none
            class(radiation), intent(inout) :: this
            real(8) :: A,delta_O
            real,dimension(this%Nw+1) :: spectrum2
            
            A = 2.961528977012875e-17
            delta_O = this%dx1*this%dy1/(this%z**2)
            this%spectrum =  A*(sum(this%Ew_im**2,1)+sum(this%Ew_re**2,1))
            spectrum2 = sum(sum(this%spectrum,1),1)*delta_O

        end subroutine calculate_spectrum


        subroutine write_radiation(this,file,tag,id,rad)

            implicit none

            class(radiation), intent(inout) :: this
            class(hdf5file), intent(in) :: file
            integer, intent(in) :: tag
            integer, intent(inout) :: id
            integer, intent(in) :: rad
            real, dimension(:,:), allocatable :: trad
            real, dimension(:,:,:), allocatable :: S
            real, dimension(:,:), allocatable :: trad2
            real, dimension(:,:), allocatable :: spectrum
            integer, dimension(2) :: gs
            integer, dimension(2) :: gs1
            integer :: ierr = 0
            integer :: i,j,k
            real(8) :: delta_O,dt
            if (rad == 1) then 
                allocate(trad(this%Nx,this%Ny))
                trad = sum(this%intensity,3)
                gs(1) = this%Nx
                gs(2) = this%Ny
                call pwrad(this%p,this%err,file,trad,gs,tag,id,ierr)  
                deallocate(trad)  
            else if (rad == 2) then
                allocate(spectrum(1,this%Nw+1))
                delta_O = this%dx1*this%dy1/(this%z**2)
                spectrum = reshape(sum(sum(this%spectrum,1),1),(/1,this%Nw+1/))*delta_O
                gs1(1) = 1
                gs1(2) = this%Nw+1
                call pwrad(this%p,this%err,file,spectrum,gs1,tag,id,ierr)
                deallocate(spectrum)
            else if (rad == 3) then
                allocate(S(this%Nx,this%Ny,this%Nt+1))
                allocate(trad2(this%Nx,this%Ny))
                DO j = 1,this%Ny,1
                    DO i = 1,this%Nx,1
                        Do k = 1,this%Nt+1,1
                            S(i,j,k) = this%field_E(1,i,j,k)*this%field_B(2,i,j,k) - this%field_E(2,i,j,k)*this%field_B(1,i,j,k)
                        end DO
                    end DO
                end DO
                dt = (this%t1-this%t0)/this%nt
                !trad2 = 3.4548e-33*sqrt(this%n0)*reshape(sum(sum(S,1),1),(/1,this%Nt+1/))*(this%z)**2*dt
                trad2 = 3.4548e-33*sqrt(this%n0)*sum(S,3)*(this%z)**2*dt
                gs(1) = this%Nx
                gs(2) = this%Ny
                call pwrad(this%p,this%err,file,trad2,gs,tag,id,ierr)
                deallocate(trad2)
                deallocate(S)
            end if
        end subroutine write_radiation
!
        end module radiation_class
