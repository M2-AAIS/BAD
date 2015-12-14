module mod_distance
  use mod_constants
  use mod_read_parameters
  use mod_variables
  implicit none

  private

  public :: pre_factor

  contains

    function pre_factor(S, S_c, dist_crit)
      implicit none
      real(x_precision), dimension(n_cell), intent(in) :: S, S_c, dist_crit

      real(x_precision) :: pre_factor, distance

      distance = maxval((S - S_c)/dist_crit)

      pre_factor = 1 - 0.90_x_precision * exp(distance)

    end function pre_factor

end module mod_distance
