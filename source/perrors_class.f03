! Perrors class for QuickPIC Open Source 1.0
! update: 04/18/2016

      module perrors_class

      use mpi
      use parallel_class
      
      implicit none

      private

      public :: perrors

      type perrors
      
         private
         
         class(parallel), pointer, public :: p => null()
         integer :: eunit = 2
         integer :: monitor = 0
         
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
         integer :: ierror

         write(str,'(I8.8)') prl%getidproc()
                                   
         this%p => prl
         this%eunit = eunit
         this%monitor = monitor
         
         if (this%p%getidproc() == 0) then
            
            call system('mkdir ./ELOG')
            
         endif
         
         call MPI_BARRIER(this%p%getlworld(),ierror)
         
         call set_ename(eunit,'./ELOG/elog-'//trim(str))
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
      end module perrors_class
