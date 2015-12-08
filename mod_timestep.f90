! Module that gives the timestep

module mod_timestep
  use mod_constants
  use mod_read_parameters
  implicit none

  private

  public :: timestep
  
contains
  
  subroutine timestep (state_in, dt_T, dt_nu)
    implicit none
    
    type (state), intent(in)              :: state_in
    real(kind = x_precision), dimension(n_cell),  intent(out) :: dt_nu, dt_T

    dt_T = state_in%Cv * state_in%T / (state_in%Qplus - state_in%Qminus) / cst_dt
    dt_nu  = dt_T / (state_in%H / x_state%x)**(2._x_precision)
  end subroutine timestep
  
end module mod_timestep

