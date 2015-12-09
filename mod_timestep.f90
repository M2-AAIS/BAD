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
!    real(kind = x_precision), dimension(n_cell)               :: diffQ

  !  diffQ = abs(state_in%Qplus - state_in%Qminus)

  !  where (diffQ <= 0.00001)
  !     dt_T = state_in%Cv * state_in%T / 0.00001 / cst_dt_T
  !  elsewhere
  !     dt_T = state_in%Cv * state_in%T / diffQ / cst_dt_T
  !  end where
      ! dt_nu  = dt_T / (state_in%H / x_state%x)**(2._x_precision) / cst_dt_nu
      !  dt_nu = x_state%x**2._x_precision / state_in%nu 

    dt_nu = x_state%x / abs(state_in%v) / cst_dt_nu
    dt_T  = (state_in%H / x_state%x)**(2._x_precision) * dt_nu * cst_dt_nu &
         / cst_dt_T
 
  end subroutine timestep
  
end module mod_timestep






