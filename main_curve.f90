program black_hole_s_curve
  use mod_constants
  use mod_read_parameters
  use mod_s_curve

  implicit none

  ! Parameters for s_curve
  real(kind = x_precision), parameter         :: eps_in = 1.e-8_x_precision  ! Precision required for the dichotomy
  real(kind = x_precision), parameter         :: Tmin   = 2.e-1_x_precision
  real(kind = x_precision), parameter         :: Tmax   = 4.49e0_x_precision
  real(kind = x_precision), parameter         :: Smin   = 4.67e1_x_precision ! Greater values lead to error below
  real(kind = x_precision), parameter         :: Smax   = 2.36e3_x_precision

  ! Error with higher value of Smin:
  !
  ! At line 228 of file mod_s_curve.f90
  ! Fortran runtime error: Index '0' of dimension 1 of array 'sigma_t_thin' below lower bound of 1
  !
  ! This is related to the fact that any greater value would be above the first solution at the first position, though I didn’t have
  ! the time to investigate why it leads to that error then. This should be fixed to get leaner curves.

  real(kind = x_precision), dimension(n_cell) :: temperature_c
  real(kind = x_precision), dimension(n_cell) :: sigma_c

  call get_parameters()

  call make_temperature(Tmin, Tmax)

  call set_conditions(eps_in, Smin, Smax)

  call curve(temperature_c, sigma_c)

end program black_hole_s_curve
