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
      real(x_precision), dimension(n_cell), intent(in) :: S_c, dist_crit
      type(state),                          intent(in) :: state_in

      real(x_precision) :: pre_factor

      pre_factor = 1 - (0.99_x_precision*exp(-(minval((S_c - state_in%S) / dist_crit))))

    end function pre_factor

end module mod_distance
