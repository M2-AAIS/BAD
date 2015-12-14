! Module that gives the timestep

module mod_timestep
  use mod_constants
  use mod_read_parameters
  implicit none

  private

  public :: timestep_T, timestep_nu

contains

  subroutine timestep_T(state_in, dt_T)
    implicit none

    type(state),       intent(in)  :: state_in
    real(x_precision), intent(out) :: dt_T

    real(x_precision), dimension(n_cell) :: diffQ
    real(x_precision), dimension(n_cell) :: threshold

    threshold = 1.e10_x_precision
    diffQ     = abs(state_in%Qplus - state_in%Qminus)

    dt_T = minval(state_in%Cv * state_in%T / max(diffQ, threshold)) / cst_dt_T

  end subroutine timestep_T

  subroutine timestep_nu(state_in, dt_nu)
    implicit none

    type(state),       intent(in)  :: state_in
    real(x_precision), intent(out) :: dt_nu

   ! dt_nu  = dt_T / (state_in%H / x_state%x)**(2._x_precision) / cst_dt_nu
   ! dt_nu = x_state%x**2._x_precision / state_in%nu

    dt_nu = minval(x_state%x / abs(state_in%v)) / cst_dt_nu

  end subroutine timestep_nu

end module mod_timestep

