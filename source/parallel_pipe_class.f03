! Parallel class with pipelines for QuickPIC Open Source 1.0
! update: 04/18/2016

      module parallel_pipe_class

      use mpi
      use parallel_class
      
      implicit none
      
      private

      public :: parallel_pipe

      type, extends (parallel) :: parallel_pipe
      
         private
! nstage: number of pipeline stages
! stageid: pipeline stage id
! lidproc: local processor id
! lkstrt: idproc+1
! lgrp = local pipeline stage communicator
! lnvp = number of MPI nodes in the local pipeline stage
         
         integer :: nstage = 1
         integer :: stageid = 0
         integer :: lidproc
         integer :: lkstrt
         integer :: lgrp
         integer :: lnvp
         
         contains
            
         procedure :: getnstage
         procedure :: getstageid
         procedure :: getlidproc
         procedure :: getlkstrt
         procedure :: getlgrp
         procedure :: getlnvp
         procedure, private :: init_parallel_pipe
         generic :: new => init_parallel_pipe
      
      end type parallel_pipe
! used by f77 subroutines
      integer :: nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
      
      contains
!
      function getlidproc(this)

         implicit none

         class(parallel_pipe), intent(in) :: this
         integer :: getlidproc
         
         getlidproc = this%lidproc

      end function getlidproc         
!
      function getlkstrt(this)

         implicit none

         class(parallel_pipe), intent(in) :: this
         integer :: getlkstrt
         
         getlkstrt = this%lkstrt

      end function getlkstrt      
!
      function getlgrp(this)

         implicit none

         class(parallel_pipe), intent(in) :: this
         integer :: getlgrp
         
         getlgrp = this%lgrp

      end function getlgrp      
!
      function getlnvp(this)

         implicit none

         class(parallel_pipe), intent(in) :: this
         integer :: getlnvp
         
         getlnvp = this%lnvp

      end function getlnvp    
!
      function getnstage(this)

         implicit none

         class(parallel_pipe), intent(in) :: this
         integer :: getnstage
         
         getnstage = this%nstage

      end function getnstage         
!
      function getstageid(this)

         implicit none

         class(parallel_pipe), intent(in) :: this
         integer :: getstageid
         
         getstageid = this%stageid

      end function getstageid      
!    
      subroutine init_parallel_pipe(this,nst)
      
         implicit none

         class(parallel_pipe), intent(inout) :: this
         integer, intent(in) :: nst
! local data
         integer :: ierror         
         integer :: idproc, llworld, nvp
         integer :: lnvp, stageid, lidproc, llgrp
         
         call this%parallel%new        
         idproc = this%getidproc()
         llworld = this%getlworld()
         nvp = this%getnvp()
         lnvp = nvp / nst
         this%nstage = nst         
         lidproc = mod(idproc,lnvp)  
         stageid = int(idproc/lnvp) 
         call MPI_COMM_SPLIT(llworld,stageid,lidproc,llgrp,ierror)
         call MPI_COMM_RANK(llgrp,this%lidproc,ierror)
         call MPI_COMM_SIZE(llgrp,this%lnvp,ierror)

         this%stageid = stageid
         this%lgrp = llgrp
         this%lkstrt = this%lidproc+1
         
         nproc = lnvp
         lgrp = llgrp
         lworld = llworld
         mreal = this%getmreal()
         mint = this%getmint()
         mcplx = this%getmcplx()
         mdouble = this%getmdouble()
                      
      end subroutine init_parallel_pipe
!      
      end module parallel_pipe_class