! Main program for QuickPIC Open Source 1.0
! update: 04/18/2016
      program quickpic
      
      use simulation_class

      implicit none
            
      type(simulation) :: sim

      call sim%new()
      call sim%go()
      call sim%del()

      stop

      end program quickpic