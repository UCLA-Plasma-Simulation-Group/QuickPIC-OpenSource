! Perrors class for QuickPIC Open Source 1.0
! update: 04/18/2016

      module perrors_class

      use mpi
      use parallel_class
      
      implicit none

      private

      public :: perrors

      type timing
           double precision :: tfft = 0
           double precision :: tpush = 0
           double precision :: tqdeposit = 0
           double precision :: tamjdeposit = 0
           double precision :: textractpsi = 0
           double precision :: t2dpmove = 0
           double precision :: t2d = 0
           double precision :: tsolve = 0
           double precision :: tcp = 0
           double precision :: tarith = 0
           integer :: nfft = 0
      end type timing

      type perrors
      
         private
         
         class(parallel), pointer, public :: p => null()
         integer :: eunit = 2
         integer :: monitor = 0
         type(timing) :: time
         
         contains
         
         generic :: new => init_perrors
         generic :: del => end_perrors
         procedure, private :: init_perrors
         procedure, private :: end_perrors
         procedure :: equit
         procedure :: werrfl0
         procedure :: werrfl1
         procedure :: werrfl2
         procedure :: setmonitor
         procedure :: add_profile
         procedure :: wprof

      end type perrors

      integer, dimension(4), save :: itime
      double precision, save :: dtime
                  
      contains
!
      subroutine init_perrors(this,prl,eunit,monitor)
         
         implicit none
         
         class(perrors), intent(inout) :: this
         class(parallel), intent(in), pointer :: prl
         integer, intent(in) :: eunit, monitor
! local data
         character(len=20) :: str
         character(len=8) :: cn
         integer :: k, l
         integer :: ierror

         write(str,*) prl%getidproc()
         cn = '00000000'
         str = trim(adjustl(str))
         l = len_trim(str)
         k = 8
         cn(k+1-l:k) = str(1:l)
                                   
         this%p => prl
         this%eunit = eunit
         this%monitor = monitor
         
         if (this%p%getidproc() == 0) then
            
            call system('mkdir ./ELOG')
            
         endif
         
         call MPI_BARRIER(this%p%getlworld(),ierror)
         
         call set_ename(eunit,'./ELOG/elog-'//cn)
         call dtimer(dtime,itime,-1)
               
      end subroutine init_perrors      
!
      subroutine set_ename(eunit,ename)
         
         implicit none
         
! this subroutine sets the name of the error file and opens it
! ename = new error file name
         character(len=*), intent(in) :: ename
         integer, intent(in) :: eunit
         
         open(unit=eunit,file=trim(ename),form='formatted',status='repla&
         &ce')

      end subroutine set_ename
!      
      subroutine end_perrors(this)
      
         implicit none
         
         class(perrors), intent(inout) :: this
         
         close(unit=this%eunit)
      
      end subroutine end_perrors
!
      subroutine equit(this,estr)
! this subroutine handles errors and optionally logs error message
! if estr is present, prepends the node number to the output string
! how = keyword on how to handle error
! estr = optional error message string to be logged in error file
         implicit none

         class(perrors), intent(in) :: this
         character(len=*), intent(in), optional :: estr      
! loccal data
         character(len=20) :: tstr
         integer :: ierror

         call dtimer(dtime,itime,1)
         write (tstr,'(f12.3)') dtime
         if (present(estr)) write (this%eunit,*) trim(tstr)//" : "//trim(adjustl(estr))
         flush(this%eunit)
         call MPI_ABORT(this%p%getlworld(),1,ierror)
!         call MPI_FINALIZE(ierror)
         print *, 'Error!'
         stop

      end subroutine equit
!      
      subroutine werrfl0(this,estr)
! this subroutine handles error message on level 0
         implicit none

         class(perrors), intent(in) :: this
         character(len=*), intent(in) :: estr
! local data
         character(len=20) :: tstr
                  
         call dtimer(dtime,itime,1)
         write (tstr,'(f12.3)') dtime
         write (this%eunit,*) trim(tstr)//" : "//trim(adjustl(estr))
         flush(this%eunit)
         
      end subroutine werrfl0
!      
      subroutine werrfl1(this,estr)
! this subroutine handles error message on level 1
         implicit none

         class(perrors), intent(in) :: this
         character(len=*), intent(in) :: estr
! local data
         character(len=20) :: tstr

         if (this%monitor < 1) then
            
            return
            
         else
					   
            call dtimer(dtime,itime,1)
            write (tstr,'(f12.3)') dtime
            write (this%eunit,*) trim(tstr)//" : "//trim(adjustl(estr))
            flush(this%eunit)
         
         endif
         
      end subroutine werrfl1
!      
      subroutine werrfl2(this,estr)
! this subroutine handles error message on level 2
         implicit none

         class(perrors), intent(in) :: this
         character(len=*), intent(in) :: estr
! local data
         character(len=20) :: tstr

         if (this%monitor < 2) then
            
            return
            
         else
					   
   	        call dtimer(dtime,itime,1)
            write (tstr,'(f12.3)') dtime
            write (this%eunit,*) trim(tstr)//" : "//trim(adjustl(estr))
            flush(this%eunit)
         
         endif
         
      end subroutine werrfl2
!      
      subroutine setmonitor(this,moniter)
 
         implicit none

         class(perrors), intent(inout) :: this
         integer, intent(in) :: moniter
         
         this%monitor = moniter
      
      end subroutine
!
      subroutine add_profile(this,event)
! this subroutine handles error message on level 1
         implicit none

         class(perrors), intent(inout) :: this
         character(len=*), intent(in) :: event
					   
         call dtimer(dtime,itime,1)
         select case (event)
         case ('fft_start')
            this%time%tfft = this%time%tfft - dtime
            this%time%nfft = this%time%nfft + 1
            return
         case ('fft_end')
            this%time%tfft = this%time%tfft + dtime
            return
         case ('push_start')
            this%time%tpush = this%time%tpush - dtime
            return
         case ('push_end')
            this%time%tpush = this%time%tpush + dtime
            return
         case ('qdeposit_start')
            this%time%tqdeposit = this%time%tqdeposit - dtime
            return
         case ('qdeposit_end')
            this%time%tqdeposit = this%time%tqdeposit + dtime
            return
         case ('amjdeposit_start')
            this%time%tamjdeposit = this%time%tamjdeposit - dtime
            return
         case ('amjdeposit_end')
            this%time%tamjdeposit = this%time%tamjdeposit + dtime
            return
         case ('extractpsi_start')
            this%time%textractpsi = this%time%textractpsi - dtime
            return
         case ('extractpsi_end')
            this%time%textractpsi = this%time%textractpsi + dtime
            return
         case ('2dpmove_start')
            this%time%t2dpmove = this%time%t2dpmove - dtime
            return
         case ('2dpmove_end')
            this%time%t2dpmove = this%time%t2dpmove + dtime
            return
         case ('solve_start')
            this%time%tsolve = this%time%tsolve - dtime
            return
         case ('solve_end')
            this%time%tsolve = this%time%tsolve + dtime
            return
         case ('cp_start')
            this%time%tcp = this%time%tcp - dtime
            return
         case ('cp_end')
            this%time%tcp = this%time%tcp + dtime
            return
         case ('arith_start')
            this%time%tarith = this%time%tarith - dtime
            return
         case ('arith_end')
            this%time%tarith = this%time%tarith + dtime
            return
         case ('2d_start')
            this%time%t2d = this%time%t2d - dtime
            return
         case ('2d_end')
            this%time%t2d = this%time%t2d + dtime
            return
         end select
      end subroutine add_profile
!      
      subroutine wprof(this)
! this subroutine handles error message on level 1
         implicit none

         class(perrors), intent(inout) :: this
! local data
         character(len=20) :: tstr
					   
         write (tstr,'(f16.4)') this%time%tfft
         call werrfl0(this,'fft time: '//trim(tstr))
         write (tstr,'(i16)') this%time%nfft
         call werrfl0(this,'total number of fft: '//trim(tstr))
         write (tstr,'(f16.4)') this%time%tfft/this%time%nfft
         call werrfl0(this,'each fft time: '//trim(tstr))
         write (tstr,'(f16.4)') this%time%tpush
         call werrfl0(this,'push time: '//trim(tstr))
         write (tstr,'(f16.4)') this%time%tqdeposit
         call werrfl0(this,'qdeposit time: '//trim(tstr))
         write (tstr,'(f16.4)') this%time%tamjdeposit
         call werrfl0(this,'amjdeposit time: '//trim(tstr))
         write (tstr,'(f16.4)') this%time%textractpsi
         call werrfl0(this,'extractpsi time: '//trim(tstr))
         write (tstr,'(f16.4)') this%time%t2dpmove
         call werrfl0(this,'2dpmove time: '//trim(tstr))
         write (tstr,'(f16.4)') this%time%tsolve
         call werrfl0(this,'solve time: '//trim(tstr))
         write (tstr,'(f16.4)') this%time%tcp
         call werrfl0(this,'cp time: '//trim(tstr))
         write (tstr,'(f16.4)') this%time%tarith
         call werrfl0(this,'arith time: '//trim(tstr))
         write (tstr,'(f16.4)') this%time%tpush+this%time%tqdeposit+&
         &this%time%tamjdeposit+this%time%textractpsi+this%time%t2dpmove
         call werrfl0(this,'2d part total time: '//trim(tstr))
         write (tstr,'(f16.4)') this%time%tsolve+this%time%tarith
         call werrfl0(this,'2d solve total time: '//trim(tstr))
         write (tstr,'(f16.4)') this%time%t2d
         call werrfl0(this,'2d time: '//trim(tstr))
                  
      end subroutine wprof
!      
      end module perrors_class
