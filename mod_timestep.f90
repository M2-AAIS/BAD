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

    type(state),                          intent(in)  :: state_in
    real(x_precision), dimension(n_cell), intent(out) :: dt_nu, dt_T

    real(x_precision), dimension(n_cell) :: diffQ
    real(x_precision)                    :: threshold

    threshold = 1e10_x_precision
    diffQ = abs(state_in%Qplus - state_in%Qminus)

    where (diffQ <= threshold)
       dt_T = state_in%Cv * state_in%T / threshold / cst_dt_T
    elsewhere
       dt_T = state_in%Cv * state_in%T / diffQ / cst_dt_T
    end where
   ! dt_nu  = dt_T / (state_in%H / x_state%x)**(2._x_precision) / cst_dt_nu
   ! dt_nu = x_state%x**2._x_precision / state_in%nu

    dt_nu = x_state%x / abs(state_in%v) / cst_dt_nu
    where (dt_nu < dt_T*100._x_precision)
       dt_T = dt_nu/100._x_precision
    end where
   ! dt_T  = (state_in%H / x_state%x)**(2._x_precision) * dt_nu * cst_dt_nu &
   !      / cst_dt_T

   ! dt_T  =  1._x_precision / params%alpha / x_state%Omega / cst_dt_T
   ! dt_T  = state_in%Cs / params%alpha / state_in%H / x_state%Omega**2._x_precision
  end subroutine timestep

end module mod_timestep

