program black_hole_s_curve
  use mod_constants
  use mod_read_parameters
  use mod_s_curve

  implicit none

  ! Parameters for s_curve
  real(kind = x_precision), parameter         :: eps    = 1.e-8_x_precision ! Precision required for the dichotomy
  real(kind = x_precision), parameter         :: T_min  = 2.e-1_x_precision
  real(kind = x_precision), parameter         :: T_max  = 4.49e0_x_precision
  real(kind = x_precision), parameter         :: S_min  = 4.5e1_x_precision ! Greater values lead to FPE
  real(kind = x_precision), parameter         :: S_max  = 2.5e3_x_precision

  real(kind = x_precision), dimension(n_cell) :: temperature_c
  real(kind = x_precision), dimension(n_cell) :: sigma_c

  call get_parameters()

  call make_temperature(T_min, T_max)

  call curve(eps, S_min, S_max, temperature_c, sigma_c)

end program black_hole_s_curve
