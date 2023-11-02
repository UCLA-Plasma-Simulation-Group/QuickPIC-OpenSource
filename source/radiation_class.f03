! radiation module for QuickPIC Open Source 

        module radiation_class
        
        use input_class
        use perrors_class
        use parallel_pipe_class
        use hdf5io_class
        ! use integral

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
         real,allocatable :: field_E(:,:,:,:),E_old(:,:,:),t_det(:,:)
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
        !  procedure, private :: calculate_radiation2
        !  procedure, private :: calculate_intensity2
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
            ! integer::integral_type        
            integer :: indx, indy, indz,type
            character(len=20) :: sn,s1
            character(len=:), allocatable :: diag
            character(len=25), save :: sname = 'init_radiation:'
            ! procedure(event), pointer :: calculate_inte !=> null ()

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
            this%z = z + min
            call input%get(s1//'.time(1)',min )
            call input%get(s1//'.time(2)',max )
            call input%get(s1//'.nt',Nt )
            this%t0 = min
            this%t1 = max
            this%Nt = Nt

            this%dx = (this%x1 - this%x0)/(this%Nx-1)
            this%dy = (this%y1 - this%y0)/(this%Ny-1)
            this%dt = (this%t1 - this%t0)/this%Nt



            call input%info(s1//'.diag',n_children=n)
            
            DO ii = 1,n
                write (sn,'(I3.3)') ii
                call input%get(s1//'.diag('//trim(sn)//')',diag)
                select case(trim(diag))
                case('intensity' ) 
                    if ( .not. allocated(this%intensity)) then
                        allocate(this%intensity(this%Nx,this%Ny,this%Nt+1))
                        this%intensity = 0
                    end if
                case('intens_t' ) 
                    if ( .not. allocated(this%intensity)) then
                        allocate(this%intensity(this%Nx,this%Ny,this%Nt+1))
                        this%intensity = 0
                    end if
                case('spectrum')
                    if ( .not. allocated(this%Ew_re)) then
                    call input%get(s1//'.spectrum.frequence(1)',min )
                    call input%get(s1//'.spectrum.frequence(2)',max)
                    call input%get(s1//'.spectrum.nw',Nw )
                    call input%get(s1//'.spectrum.Nx',Nx1 )
                    call input%get(s1//'.spectrum.Ny',Ny1 )
                    call input%get(s1//'.spectrum.in_n',in_n )
                    ! call input%get(s1//'.spectrum.integral_type',integral_type )
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
                    ! if (integral_type==0) then 
                    ! calculate_inte => calculate_inte0
                    ! else if (integral_type==1) then 
                    ! calculate_inte => calculate_inte1
                    ! else if (integral_type==2) then 
                    ! calculate_inte => calculate_inte2
                    ! end if
                    end if
                case('spectrum_intens')
                    if ( .not. allocated(this%Ew_re)) then
                    call input%get(s1//'.spectrum.frequence(1)',min )
                    call input%get(s1//'.spectrum.frequence(2)',max)
                    call input%get(s1//'.spectrum.nw',Nw )
                    call input%get(s1//'.spectrum.Nx',Nx1 )
                    call input%get(s1//'.spectrum.Ny',Ny1 )
                    call input%get(s1//'.spectrum.in_n',in_n )
                    ! call input%get(s1//'.spectrum.integral_type',integral_type )
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
                    ! if (integral_type==0) then 
                    ! calculate_inte => calculate_inte0
                    ! else if (integral_type==1) then 
                    ! calculate_inte => calculate_inte1
                    ! else if (integral_type==2) then 
                    ! calculate_inte => calculate_inte2
                    ! end if
                    end if
                case('field')
                    allocate(this%field_E(3,this%Nx,this%Ny,this%Nt+1))
                    this%field_E = 0
                    ! allocate(this%E_field(3,this%Nt+1))
                    ! this%E_old = 0
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
            real :: x,y,z,tempz,xx,yy,zz
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
                    posi = (/x,y,n1*z+n2*tempz/)
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
                            xx = this%x0+(i-1)*this%dx
                            yy = this%y1-(j-1)*this%dy
                            zz = (this%z**2-n2*(xx**2+yy**2))**0.5
                            posi_x = (/xx,yy,zz/)
                            vec_r = posi_x - posi
                            r = (vec_r(1)**2 +vec_r(2)**2 +vec_r(3)**2)**0.5
                            vec_n = vec_r/r
!                            t_obs = (t + r - this%t0)/this%dt
                            t_obs = (n1*t +n2*tempz + r - this%t0)/this%dt
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
                                print*,(n1*t + r - this%t0)
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
            integer :: t_n,t_n_new
            real :: x,y,z,tempz,xx,yy,zz
            real :: x_new,y_new,z_new,tempz_new
            real :: t,t_obs,m
            real :: t_new,t_obs_new,t_temp
            real :: delta_t
            real :: dx,dy,dt
            real :: r ,r_new
            real :: deltax, deltay, deltaz
            real,dimension(3) :: field_E,field_B,S
            real(8) :: A,n0
            real,dimension(3) :: posi,vel,acc,posi_x,vec_r,vec_n,p,temp1,temp2,pnew,pold,vnew,vold,force
            real,dimension(3) :: posi_new,vec_r_new
            real, dimension(:,:), allocatable :: E_field
            real :: temp3,temp4
            character(len=25), save :: sname = 'calculate_field:'
            call this%err%werrfl2(class//sname//' started')
            z = s_beam
            deltax = this%delta(1)
            deltay = this%delta(2)
            deltaz = this%delta(3)
            tt3 = size(part,1)
            A = 1.528671463266763e-08
            allocate(E_field(3,this%Nt+1))
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
                posi = (/x,y,n1*z+n2*tempz/)
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
                x_new = x + vnew(1)*dt
                y_new = y+ vnew(2)*dt
                z_new = s_beam+dt
                tempz_new = z+(1-vnew(3))*dt
                posi_new = (/x_new,y_new,n1*z_new+n2*tempz_new/)
                t_new = s_beam+dt + tempz_new
                DO i=1,this%Nx,1
                    DO j = 1,this%Ny,1
                        xx = this%x0+(i-1)*this%dx
                        yy = this%y1-(j-1)*this%dy
                        zz = (this%z**2-n2*(xx**2+yy**2))**0.5
                        posi_x = (/xx,yy,zz/)
                        vec_r = posi_x - posi
                        vec_r_new = posi_x - posi_new
                        r_new = (vec_r_new(1)**2 +vec_r_new(2)**2 +vec_r_new(3)**2)**0.5
                        vec_n = vec_r/r
                        !t_obs = (t + r-this%t0)/this%dt
                        t_obs = (n1*t +n2*tempz + r - this%t0)/this%dt
                        t_obs_new = (n1*t_new +n2*tempz_new + r_new - this%t0)/this%dt
                        if (t_obs > 0 .AND. t_obs < this%Nt) then                            
                            t_n = t_obs
                            t_n_new = t_obs_new
                            temp1 = cross(vec_n-vel,acc)
                            temp2 = cross(vec_n,temp1)
                            temp4 = (1 - vec_n(1)*vel(1)-vec_n(2)*vel(2)-vec_n(3)*vel(3))
                            field_E = (temp2/(temp4**3)/(n1*r+n2*this%z))
                            t_temp = t_obs
                            if (t_n_new == t_n) then
                                E_field(:,t_n) =E_field(:,t_n)+(t_obs_new-t_obs)*field_E
                            else 
                                do while (t_n < t_n_new) 
                                    m = (t_n + 1-t_temp)
                                    E_field(:,t_n) = E_field(:,t_n) + m*field_E
                                    t_n = t_n+1
                                    t_temp = t_n
                                end do
                                m = (t_obs_new - t_n)
                                E_field(:,t_n_new) = E_field(:,t_n_new) + m*field_E
                            end if                 
                            this%field_E(:,i,j,:) = this%field_E(:,i,j,:) + E_field
                            ! this%t_det(i,j) = t_obs
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
            real,dimension(3) :: posi,vel,acc,posi_x,vec_r,vec_n,p,pnew,pold,vnew,vold,force
            real,dimension(3) ::temp1,temp2,temp7
            real :: temp3,temp4,temp5_im,temp5_Re,temp6,temp_w
            real :: k1,k2,inte_re,inte_im
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

                            temp_w =w*(t_obs*n1+n2*(t-(vec_n(1)*posi(1)+vec_n(2)*posi(2)+vec_n(3)*posi(3))))
                            temp6 =  - vec_n(1)*acc(1)-vec_n(2)*acc(2)-vec_n(3)*acc(3)
                            temp5_Re = cos(temp_w)
                            temp5_im = sin(temp_w)
                            temp6 =  - vec_n(1)*acc(1)-vec_n(2)*acc(2)-vec_n(3)*acc(3)
                            k1 = w*temp4
                            k2 = w*temp6
                            call calculate_inte2(k1,k2,delta_t,inte_re,inte_im)   
                            temp7=  temp2/temp4**2
                            Ew_re = temp7*(temp5_re*inte_re-temp5_im*inte_im)                  
                            Ew_im =temp7*(temp5_re*inte_im+temp5_im*inte_re)

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
            real,dimension(3) ::temp0,temp1,temp2,temp7
            real :: temp3,temp4,temp5_im,temp5_Re,temp_w,temp6
            real :: k1,k2,inte_re,inte_im
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
                                temp_w =w*(t_obs*n1+n2*(t-(vec_n(1)*posi_new(1)+vec_n(2)*posi_new(2)+vec_n(3)*posi_new(3))))
                                temp5_Re = cos(temp_w)
                                temp5_im = sin(temp_w)
                                
                                temp6 =  - vec_n(1)*acc(1)-vec_n(2)*acc(2)-vec_n(3)*acc(3)
                                k1 = w*temp4
                                k2 = w*temp6
                                call calculate_inte2(k1,k2,in_t,inte_re,inte_im)    
                                temp7=  temp2/temp4**2
                                Ew_re = temp7*(temp5_re*inte_re-temp5_im*inte_im)                  
                                Ew_im =temp7*(temp5_re*inte_im+temp5_im*inte_re)

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
            real, dimension(:,:), allocatable :: trad,rad_t
            real, dimension(:,:,:), allocatable :: S
            real, dimension(:,:), allocatable :: trad2
            real, dimension(:,:), allocatable :: spectrum
            integer, dimension(2) :: gs,gs3
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
            else if (rad == 3) then 
                allocate(rad_t(1,this%Nt+1))
                delta_O = this%dx*this%dy/(this%z**2)
                rad_t = reshape(sum(sum(this%intensity,1),1),(/1,this%Nt+1/))*delta_O
                gs3(1) = 1
                gs3(2) = this%Nt+1
                call pwrad(this%p,this%err,file,rad_t,gs3,tag,id,ierr)  
                deallocate(rad_t) 
            else if (rad == 2) then
                allocate(spectrum(1,this%Nw+1))
                !ALLOCATE(spectrum(this%Nx1,this%Ny1))
                delta_O = this%dx1*this%dy1/(this%z**2)
                spectrum = reshape(sum(sum(this%spectrum,1),1),(/1,this%Nw+1/))*delta_O
                !spectrum=sum(this%spectrum,3)
                gs1(1) = 1
                gs1(2) = this%Nw+1
                !gs1(1)=this%Nx1
                !gs1(2)=this%Ny1
                call pwrad(this%p,this%err,file,spectrum,gs1,tag,id,ierr)
                deallocate(spectrum)
            else if (rad == 4) then
                allocate(S(this%Nx,this%Ny,this%Nt+1))
                allocate(trad2(this%Nx,this%Ny))
                DO j = 1,this%Ny,1
                    DO i = 1,this%Nx,1
                        Do k = 1,this%Nt+1,1
                            S(i,j,k) = this%field_E(1,i,j,k)**2+this%field_E(2,i,j,k)**2+this%field_E(3,i,j,k)**2
                            !S(i,j,k) = this%field_E(1,i,j,k)*this%field_B(2,i,j,k) - this%field_E(2,i,j,k)*this%field_B(1,i,j,k)
                        end DO
                    end DO
                end DO
                dt = (this%t1-this%t0)/this%nt
                !trad2 = 3.4548e-33*sqrt(this%n0)*reshape(sum(sum(S,1),1),(/1,this%Nt+1/))*(this%z)**2*dt
                trad2 = 3.4548e-33*sqrt(this%n0)*sum(S,3)*dt
                gs(1) = this%Nx
                gs(2) = this%Ny
                call pwrad(this%p,this%err,file,trad2,gs,tag,id,ierr)
                deallocate(trad2)
                deallocate(S)
            end if
        end subroutine write_radiation
!
        ! subroutine calculate_inte1(k1,k2,dt,re,im)
        !     implicit none
        !     real(8), INTENT(IN) :: k1, k2
        !     real(8) ,INTENT(IN) :: dt
        !     real(8), INTENT(INOUT)  :: re,im

        !     real(8) :: dtt,pi,xi1,xi2,u1,u2,v1,v2,d1,d2,temp_c1,temp_c2
        !     real(8) :: co,si,temp1_re,temp1_im,temp2_re,temp2_im,xi11,xi22
        !     Logical :: flag
        !     pi=3.141592653589793d0
        !     dtt =0.5*dt

        !     if (k2<0) then 
        !         d1 = k1
        !         d2 = -k2

        !         temp_c1 = (2*pi/d2)**0.5/2
        !         temp_c2 = d1**2/4/d2
        !         co = cos(temp_c2)
        !         si = sin(temp_c2)
        !         temp1_re = temp_c1*(co-si)
        !         temp1_im = temp_c1*(co+si)
        !         ! print*,"temp1:",temp1_re,temp1_im

        !         xi1 =d1/d2**0.5
        !         xi2 = 2*d2**0.5*dtt
        !         xi11 = (xi1-xi2)/2.8284271247461903d0
        !         xi22 = (xi1+xi2)/2.8284271247461903d0
        !         call erfc_c(xi11, xi11, u1, v1, flag)
        !         call erfc_c(xi22, xi22, u2, v2, flag)

        !     else 
        !         d1 = k1
        !         d2 = k2

        !         temp_c1 = (2*pi/d2)**0.5/2
        !         temp_c2 = -d1**2/(4*d2)
        !         co = cos(temp_c2)
        !         si = sin(temp_c2)
        !         temp1_re = temp_c1*(co+si)
        !         temp1_im = temp_c1*(si-co)

                
                
        !         xi1 =d1/d2**0.5
        !         xi2 = 2*d2**0.5*dtt
        !         xi11 = (xi1+xi2)/2.8284271247461903d0
        !         xi22 = (xi1-xi2)/2.8284271247461903d0            
        !         call erfc_c(xi11, -xi11, u1, v1, flag)
        !         call erfc_c(xi22, -xi22, u2, v2, flag)
        !     end if
        !         temp2_re = (v1-v2)/2
        !         temp2_im = (u2-u1)/2
        !         ! print*,"temp2:",temp2_re,temp2_im
        !         ! print*,u1,v1
        !         ! print*,temp2_re,temp2_im
        !         re = temp1_re*temp2_re - temp1_im*temp2_im
        !         im = temp1_re*temp2_im + temp1_im*temp2_re 
            
        !     end subroutine 
        subroutine calculate_inte2(k1,k2,dt,re,im)
            implicit none
            REAL(8),  INTENT(inout) :: re,im
            real(8),  INTENT(IN) :: k1, k2
            real(8),  INTENT(IN):: dt
            real(8) :: dtt
            dtt =0.5*dt
            re = 2*sin(k1*dtt)/k1
            im = k2*(2*dtt**2*sin(k1*dtt)/k1 + 4*dtt*cos(k1*dtt)/k1**2 - 4*sin(k1*dtt)/k1**3)
        end SUBROUTINE 


            
        ! subroutine erfc_c(xi,yi,u,v,flag)
        ! real(8),INTENT(in) ::xi,yi
        ! real(8),INTENT(inout) :: u,v
        ! real(8) :: wx,wy,temp1,temp2_co,temp2_si
        ! Logical,INTENT(inout):: flag
        ! call wofz(-yi,xi,wx,wy,flag)
        ! temp1 = exp(yi**2-xi**2)
        ! temp2_co = cos(2*xi*yi)
        ! temp2_si = sin(2*xi*yi)
        ! u = (wx*temp2_co+wy*temp2_si)*temp1
        ! v = (wy*temp2_co-wx*temp2_si)*temp1
        ! end subroutine


        ! Subroutine wofz(xi, yi, u, v, flag)

        ! !  GIVEN A COMPLEX NUMBER Z = (XI,YI), THIS SUBROUTINE COMPUTES
        ! !  THE VALUE OF THE FADDEEVA-FUNCTION W(Z) = EXP(-Z**2)*ERFC(-I*Z),
        ! !  WHERE ERFC IS THE COMPLEX COMPLEMENTARY ERROR-FUNCTION AND I
        ! !  MEANS SQRT(-1).
        ! !  THE ACCURACY OF THE ALGORITHM FOR Z IN THE 1ST AND 2ND QUADRANT
        ! !  IS 14 SIGNIFICANT DIGITS; IN THE 3RD AND 4TH IT IS 13 SIGNIFICANT
        ! !  DIGITS OUTSIDE A CIRCULAR REGION WITH RADIUS 0.126 AROUND A ZERO
        ! !  OF THE FUNCTION.
        ! !  ALL REAL VARIABLES IN THE PROGRAM ARE DOUBLE PRECISION.
        ! !
        ! !
        ! !  THE CODE CONTAINS A FEW COMPILER-DEPENDENT PARAMETERS :
        ! !     RMAXREAL = THE MAXIMUM VALUE OF RMAXREAL EQUALS THE ROOT OF
        ! !                RMAX = THE LARGEST NUMBER WHICH CAN STILL BE
        ! !                IMPLEMENTED ON THE COMPUTER IN DOUBLE PRECISION
        ! !                FLOATING-POINT ARITHMETIC
        ! !     RMAXEXP  = LN(RMAX) - LN(2)
        ! !     RMAXGONI = THE LARGEST POSSIBLE ARGUMENT OF A DOUBLE PRECISION
        ! !                GONIOMETRIC FUNCTION (DCOS, DSIN, ...)
        ! !  THE REASON WHY THESE PARAMETERS ARE NEEDED AS THEY ARE DEFINED WILL
        ! !  BE EXPLAINED IN THE CODE BY MEANS OF COMMENTS
        ! !
        ! !
        ! !  PARAMETER LIST
        ! !     XI     = REAL      PART OF Z
        ! !     YI     = IMAGINARY PART OF Z
        ! !     U      = REAL      PART OF W(Z)
        ! !     V      = IMAGINARY PART OF W(Z)
        ! !     FLAG   = AN ERROR FLAG INDICATING WHETHER OVERFLOW WILL
        ! !              OCCUR OR NOT; TYPE LOGICAL;
        ! !              THE VALUES OF THIS VARIABLE HAVE THE FOLLOWING
        ! !              MEANING :
        ! !              FLAG=.FALSE. : NO ERROR CONDITION
        ! !              FLAG=.TRUE.  : OVERFLOW WILL OCCUR, THE ROUTINE
        ! !                             BECOMES INACTIVE
        ! !  XI, YI      ARE THE INPUT-PARAMETERS
        ! !  U, V, FLAG  ARE THE OUTPUT-PARAMETERS
        ! !
        ! !  FURTHERMORE THE PARAMETER FACTOR EQUALS 2/SQRT(PI)
        ! !
        ! !  THE ROUTINE IS NOT UNDERFLOW-PROTECTED BUT ANY VARIABLE CAN BE
        ! !  PUT TO 0 UPON UNDERFLOW;
        ! !
        ! !  REFERENCE - GPM POPPE, CMJ WIJERS; MORE EFFICIENT COMPUTATION OF
        ! !  THE COMPLEX ERROR-FUNCTION, ACM TRANS. MATH. SOFTWARE.
        ! !
        ! !
        ! !
        ! !
        ! !
        ! Implicit Double Precision (A-H, O-Z)
        ! !
        ! integer :: i,j,kapn,n,np1,nu
        ! Logical a, b, flag
        ! Parameter (factor=1.12837916709551257388D0, rmaxreal=0.5D+154, rmaxexp=708.503061461606D0, rmaxgoni=3.53711887601422D+15)
        ! !
        ! flag = .False.
        ! !
        ! xabs = dabs(xi)
        ! yabs = dabs(yi)
        ! x = xabs/6.3
        ! y = yabs/4.4
        ! !
        ! !
        ! !     THE FOLLOWING IF-STATEMENT PROTECTS
        ! !     QRHO = (X**2 + Y**2) AGAINST OVERFLOW
        ! !
        ! If ((xabs>rmaxreal) .Or. (yabs>rmaxreal)) Goto 100
        ! !
        ! qrho = x**2 + y**2
        ! !
        ! xabsq = xabs**2
        ! xquad = xabsq - yabs**2
        ! yquad = 2*xabs*yabs
        ! !
        ! a = qrho < 0.085264D0
        ! !
        ! If (a) Then
        ! !
        ! !  IF (QRHO.LT.0.085264D0) THEN THE FADDEEVA-FUNCTION IS EVALUATED
        ! !  USING A POWER-SERIES (ABRAMOWITZ/STEGUN, EQUATION (7.1.5), P.297)
        ! !  N IS THE MINIMUM NUMBER OF TERMS NEEDED TO OBTAIN THE REQUIRED
        ! !  ACCURACY
        ! !
        !     qrho = (1-0.85*y)*dsqrt(qrho)
        !     n = idnint(6+72*qrho)
        !     j = 2*n + 1
        !     xsum = 1.0/j
        !     ysum = 0.0D0
        !     Do i = n, 1, -1
        !     j = j - 2
        !     xaux = (xsum*xquad-ysum*yquad)/i
        !     ysum = (xsum*yquad+ysum*xquad)/i
        !     xsum = xaux + 1.0/j
        !     End Do
        !     u1 = -factor*(xsum*yabs+ysum*xabs) + 1.0
        !     v1 = factor*(xsum*xabs-ysum*yabs)
        !     daux = dexp(-xquad)
        !     u2 = daux*dcos(yquad)
        !     v2 = -daux*dsin(yquad)
        ! !
        !     u = u1*u2 - v1*v2
        !     v = u1*v2 + v1*u2
        ! !
        ! Else
        ! !
        ! !  IF (QRHO.GT.1.O) THEN W(Z) IS EVALUATED USING THE LAPLACE
        ! !  CONTINUED FRACTION
        ! !  NU IS THE MINIMUM NUMBER OF TERMS NEEDED TO OBTAIN THE REQUIRED
        ! !  ACCURACY
        ! !
        ! !  IF ((QRHO.GT.0.085264D0).AND.(QRHO.LT.1.0)) THEN W(Z) IS EVALUATED
        ! !  BY A TRUNCATED TAYLOR EXPANSION, WHERE THE LAPLACE CONTINUED FRACTION
        ! !  IS USED TO CALCULATE THE DERIVATIVES OF W(Z)
        ! !  KAPN IS THE MINIMUM NUMBER OF TERMS IN THE TAYLOR EXPANSION NEEDED
        ! !  TO OBTAIN THE REQUIRED ACCURACY
        ! !  NU IS THE MINIMUM NUMBER OF TERMS OF THE CONTINUED FRACTION NEEDED
        ! !  TO CALCULATE THE DERIVATIVES WITH THE REQUIRED ACCURACY
        ! !
        ! !
        !     If (qrho>1.0) Then
        !     h = 0.0D0
        !     kapn = 0
        !     qrho = dsqrt(qrho)
        !     nu = idint(3+(1442/(26*qrho+77)))
        !     Else
        !     qrho = (1-y)*dsqrt(1-qrho)
        !     h = 1.88*qrho
        !     h2 = 2*h
        !     kapn = idnint(7+34*qrho)
        !     nu = idnint(16+26*qrho)
        !     End If
        ! !
        !     b = (h>0.0)
        ! !
        !     If (b) qlambda = h2**kapn
        ! !
        !     rx = 0.0
        !     ry = 0.0
        !     sx = 0.0
        !     sy = 0.0
        ! !
        !     Do n = nu, 0, -1
        !     np1 = n + 1
        !     tx = yabs + h + np1*rx
        !     ty = xabs - np1*ry
        !     c = 0.5/(tx**2+ty**2)
        !     rx = c*tx
        !     ry = c*ty
        !     If ((b) .And. (n<=kapn)) Then
        !         tx = qlambda + sx
        !         sx = rx*tx - ry*sy
        !         sy = ry*tx + rx*sy
        !         qlambda = qlambda/h2
        !     End If
        !     End Do
        ! !
        !     If (h==0.0) Then
        !     u = factor*rx
        !     v = factor*ry
        !     Else
        !     u = factor*sx
        !     v = factor*sy
        !     End If
        ! !
        !     If (yabs==0.0) u = dexp(-xabs**2)
        ! !
        ! End If
        ! !
        ! !
        ! !
        ! !  EVALUATION OF W(Z) IN THE OTHER QUADRANTS
        ! !
        ! !
        ! If (yi<0.0) Then
        ! !
        !     If (a) Then
        !     u2 = 2*u2
        !     v2 = 2*v2
        !     Else
        !     xquad = -xquad
        ! !
        ! !
        ! !         THE FOLLOWING IF-STATEMENT PROTECTS 2*EXP(-Z**2)
        ! !         AGAINST OVERFLOW
        ! !
        !     If ((yquad>rmaxgoni) .Or. (xquad>rmaxexp)) Goto 100
        ! !
        !     w1 = 2*dexp(xquad)
        !     u2 = w1*dcos(yquad)
        !     v2 = -w1*dsin(yquad)
        !     End If
        ! !
        !     u = u2 - u
        !     v = v2 - v
        !     If (xi>0.0) v = -v
        ! Else
        !     If (xi<0.0) v = -v
        ! End If
        ! !
        ! Return
        ! !
        ! 100 flag = .True.
        ! Return
        ! !
        ! End Subroutine wofz
        end module radiation_class
