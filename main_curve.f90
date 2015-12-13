program black_hole_s_curve
  use mod_constants
  use mod_read_parameters
  use mod_s_curve

  implicit none

  ! Parameters for s_curve
  real(x_precision), parameter         :: eps_in = 1.e-7_x_precision  ! Precision required for the dichotomy
  real(x_precision), parameter         :: Tmin   = 10._x_precision**5.3_x_precision
  real(x_precision), parameter         :: Tmax   = 10._x_precision**7.4_x_precision
  real(x_precision), parameter         :: Smin   = 10._x_precision**1.5_x_precision
  real(x_precision), parameter         :: Smax   = 10._x_precision**3.5_x_precision

  real(x_precision), dimension(n_cell) :: temperature_c
  real(x_precision), dimension(n_cell) :: sigma_c

  call get_parameters()

  call make_temperature(Tmin / state_0%T_0, Tmax / state_0%T_0)

  call set_conditions(eps_in, Smin / state_0%S_0, Smax / state_0%S_0)

  call curve(temperature_c, sigma_c)

end program black_hole_s_curve
