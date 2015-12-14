program black_hole_s_curve
  use mod_constants
  use mod_read_parameters
  use mod_s_curve

  implicit none

  real(x_precision), dimension(n_cell) :: temperature_c
  real(x_precision), dimension(n_cell) :: sigma_c

  call get_parameters()
  call make_temperature()
  call set_conditions()

  call curve(temperature_c, sigma_c)

end program black_hole_s_curve
