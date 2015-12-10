module mod_distance
  use mod_constants
  use mod_read_parameters
  use mod_variables
  implicit none 
  
  private  

  public :: pre_factor

  contains 

    function pre_factor(state_in, S_c, dist_crit)
      implicit none
      type(state), intent(in)                                 :: state_in
      real(kind = x_precision)                                :: pre_factor, tmp
      real(kind = x_precision), dimension(n_cell), intent(in) :: S_c, dist_crit
      tmp = 1 - (0.999_x_precision*exp(-(minval((S_c - state_in%S) / dist_crit))))
      
      pre_factor = tmp
    end function pre_factor
  end module mod_distance
