module mod_distance
  use mod_constants
  use mod_read_parameters
  use mod_variables
  implicit none 
  
  private  

  public :: distance

  contains 

    subroutine distance(state_in, dist, temp_c, sigma_c)
      implicit none
      type(state)                                 :: state_in
      real(kind = x_precision), dimension(n_cell) :: dist
      real(kind = x_precision), dimension(n_cell) :: temp_c, sigma_c

      dist = sqrt ((abs(state_in%T-temp_c))**2._x_precision + & 
           (abs(state_in%S-sigma_c))**2._x_precision)
      
    end subroutine distance
  end module mod_distance
