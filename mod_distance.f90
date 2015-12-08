module mod_distance
  use mod_constants
  use mod_read_parameters
  use mod_s_curve
  implicit none 
  
  private  

  public :: distance

  contains 

    subroutine distance(state_in, dist, temp_c, sigma_c)
      implicit none
      type(state)                                 :: state_in
      real(kind = x_precision), dimension(n_cell) :: dist
      real(kind = x_precision), dimension(n_cell) :: temp_c, sigma_c

      dist = 1. !FIXME
      
    end subroutine distance
  end module mod_distance
