program black_hole_s_curve
  use mod_constants
  use mod_read_parameters
  use mod_s_curve


! CAREFUL !!! The output arrays are for the NON dimensioned values of the surface density S (not Sigma) and the Temperature  whereas the produced data files used when ploting are dimensioned.

  implicit none

  integer                 ,parameter         :: nb_it  = 1000
  integer                 ,parameter         :: max_it = 10000000 ! Maximum number of dichotomy iterations, max 1e9
  real(kind = x_precision),parameter         :: eps    = 1.e-8_x_precision ! Precision required for the dichotomy
  real(kind = x_precision),parameter         :: t_min  = 2.e-1_x_precision
  real(kind = x_precision),parameter         :: t_max  = 4.49e0_x_precision
  real(kind = x_precision),parameter         :: dt     = (t_max - t_min) / (nb_it - 1)
  real(kind = x_precision),parameter         :: S_min  = 4.5e1_x_precision ! Greater values lead to FPE
  real(kind = x_precision),parameter         :: S_max  = 2.5e3_x_precision

  real(kind = x_precision),dimension(nb_it)  :: temperature
  real(kind = x_precision),dimension(n_cell) :: temperature_c
  real(kind = x_precision),dimension(n_cell) :: sigma_c

  integer                                    :: i

  do i = 1, nb_it
    temperature(i) = dt * (i-1) + t_min
  enddo

  call get_parameters()

  call curve(nb_it, max_it, eps, temperature, s_min, s_max, temperature_c, sigma_c)

end program black_hole_s_curve
