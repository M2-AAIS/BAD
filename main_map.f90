program black_hole_maps
  use mod_constants
  use mod_read_parameters
  use mod_maps

  implicit none

  integer, parameter :: nbT = 400
  integer, parameter :: nbS = 400

  real(x_precision), dimension(n_cell,nbT,nbS) :: Q_res   ! The Q+-Q- grids, one for each position in the disk
  real(x_precision), dimension(n_cell,nbT,nbS) :: tau_res ! The tau grids, one for each position in the disk

  call get_parameters()

  call set_conditions(nbS, nbT)

  call build_grid(Q_res, tau_res)

  call save_data(Q_res, tau_res)

end program black_hole_maps
