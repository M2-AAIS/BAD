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
    real(kind = x_precision), dimension(n_cell)               :: diffQ
    real(kind = x_precision) :: threshold
    integer :: i

    threshold = 1e10_x_precision
    diffQ = abs(state_in%Qplus - state_in%Qminus)

    do i = 1, n_cell
       if (diffQ(i) <= threshold) then
          dt_T(i) = state_in%Cv(i) * state_in%T(i) / threshold / cst_dt_T
       else
          dt_T(i) = state_in%Cv(i) * state_in%T(i) / diffQ(i) / cst_dt_T
       end if
    end do
      ! dt_nu  = dt_T / (state_in%H / x_state%x)**(2._x_precision) / cst_dt_nu
      !  dt_nu = x_state%x**2._x_precision / state_in%nu 

    dt_nu = x_state%x / abs(state_in%v) / cst_dt_nu
   ! dt_T  = (state_in%H / x_state%x)**(2._x_precision) * dt_nu * cst_dt_nu &
   !      / cst_dt_T

   ! dt_T  =  1._x_precision / params%alpha / x_state%Omega / cst_dt_T
   ! dt_T  = state_in%Cs / params%alpha / state_in%H / x_state%Omega**2._x_precision
  end subroutine timestep
  
end module mod_timestep






