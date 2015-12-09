program black_hole_maps
  use mod_constants
  use mod_read_parameters
  use mod_maps

  implicit none

  real(kind = x_precision), parameter :: Tmin = 2.5e-2_x_precision
  real(kind = x_precision), parameter :: Tmax = 4.49e0_x_precision
  real(kind = x_precision), parameter :: Smin = 2.36e1_x_precision
  real(kind = x_precision), parameter :: Smax = 2.36e3_x_precision
  integer, parameter :: nbT = 400
  integer, parameter :: nbS = 400

  call get_parameters()

  call set_conditions(nbS, nbT, Tmin, Tmax, Smin, Smax)

  call build_grid()

end program black_hole_maps
