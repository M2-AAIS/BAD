! Module that gives the timestep

module mod_timestep
  use mod_constants
  implicit none

  private

  public :: timestep
  
contains
  
  subroutine timestep (state_in, dt)
    implicit none
    
    type (state), intent(in)              :: state_in
    real(kind = x_precision), intent(out) :: dt    
    
  end subroutine timestep
  
end module mod_timestep

