program black_hole_s_curve
  use mod_constants
  use mod_read_parameters
  use mod_s_curve

  implicit none

  ! Parameters for s_curve
  real(x_precision), parameter         :: eps_in = 1.e-7_x_precision  ! Precision required for the dichotomy
  real(x_precision), parameter         :: Tmin   = 2.5e-2_x_precision / 2.
  real(x_precision), parameter         :: Tmax   = 4.49e0_x_precision * 2.
  real(x_precision), parameter         :: Smin   = 2.36e1_x_precision / 5.
  real(x_precision), parameter         :: Smax   = 2.36e3_x_precision * 5. 

  real(x_precision), dimension(n_cell) :: temperature_c
  real(x_precision), dimension(n_cell) :: sigma_c

  call get_parameters()

  call make_temperature(Tmin, Tmax)

  call set_conditions(eps_in, Smin, Smax)

  call curve(temperature_c, sigma_c)

end program black_hole_s_curve
