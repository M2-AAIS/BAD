program black_hole_maps
  use mod_constants
  use mod_read_parameters
  use mod_maps

  implicit none

  integer,           parameter :: nbT = 400
  integer,           parameter :: nbS = 400
  real(x_precision), parameter :: Tmin = 10._x_precision**5.3_x_precision
  real(x_precision), parameter :: Tmax = 10._x_precision**7.4_x_precision
  real(x_precision), parameter :: Smin = 10._x_precision**1.5_x_precision
  real(x_precision), parameter :: Smax = 10._x_precision**3.5_x_precision

  real(x_precision), dimension(n_cell,nbT,nbS) :: Q_res   ! The Q+-Q- grids, one for each position in the disk
  real(x_precision), dimension(n_cell,nbT,nbS) :: tau_res ! The tau grids, one for each position in the disk

  call get_parameters()

  call set_conditions(nbS, nbT, Tmin / state_0%T_0, Tmax / state_0%T_0, Smin / state_0%S_0, Smax / state_0%S_0)

  call build_grid(Q_res, tau_res)

  call save_data(Q_res, tau_res)

end program black_hole_maps
