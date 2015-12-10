program black_hole_maps
  use mod_constants
  use mod_read_parameters
  use mod_maps

  implicit none

  integer,           parameter :: nbT = 400
  integer,           parameter :: nbS = 400
  real(x_precision), parameter :: Tmin = 2.5e-2_x_precision
  real(x_precision), parameter :: Tmax = 4.49e0_x_precision
  real(x_precision), parameter :: Smin = 2.36e1_x_precision
  real(x_precision), parameter :: Smax = 2.36e3_x_precision

  real(x_precision), dimension(n_cell,nbT,nbS) :: Q_res   ! The Q+-Q- grids, one for each position in the disk
  real(x_precision), dimension(n_cell,nbT,nbS) :: tau_res ! The tau grids, one for each position in the disk

  call get_parameters()

  call set_conditions(nbS, nbT, Tmin, Tmax, Smin, Smax)

  call build_grid(Q_res, tau_res)

  call save_data(Q_res, tau_res)

end program black_hole_maps
