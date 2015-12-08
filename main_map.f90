program black_hole_maps
  use mod_constants
  use mod_read_parameters
  use mod_maps

  implicit none

  real(kind = x_precision), parameter :: T_min   = 2.5e-2_x_precision
  real(kind = x_precision), parameter :: T_max   = 4.49e0_x_precision
  real(kind = x_precision), parameter :: S_min   = 2.36e1_x_precision
  real(kind = x_precision), parameter :: S_max   = 2.36e3_x_precision
  integer, parameter :: nbT = 10
  integer, parameter :: nbS = 10

  call get_parameters()

  call set_conditions(nbS, nbT)

  call build_grid(T_min, T_max, S_min, S_max)

end program black_hole_maps
