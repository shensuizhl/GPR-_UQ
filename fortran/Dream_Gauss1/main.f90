!#####################################################################
!                                                                    #
!  The FORTRAN programs are written by Dan Lu at the Department of   #
!  Scientific Computing at the Florida State University under        #
!  supervision of Professor Ming Ye with support of NSF-EAR grant    #
!  0911074 and DOE-SBR grant DESC0002687.                            #
!  The codes implement the DREAM algorithm developed by Jasper A.    #
!  Vrugt for Markov Chain Monte Carlo simulation.                    #
!  There is no guarantee that this program will suit the user's      #
!  needs or goals, execute efficiently and without mishap on the     #
!  user's computer, exhibit no errors or bugs, or yield a            #
!  scientifically defensible result.                                 #
!  The user is welcome to make any modifications needed to suit      #
!  his/her interpretive and modeling needs.                          #
!                                                                    #
!  Dan Lu (demi.ludan@gmail.com)                                     # 
!  Florida state university, August, 2013                            #
!                                                                    #
!#####################################################################
program main

   use rand_generator
   use pdf_density
   use dream
   use run_mod
   call read_vars( )
   call import_point( )
   call import_observ()
   call Guass() 
end program
