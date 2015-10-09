! Module that gives the timestep

module mod_timestep
  
  private

  public :: timestep
  
contains
  
  subroutine timestep (state_in)
    use mod_constants
    implicit none
    type (state), intent(in), dimension(:)  :: state_in
    
  end subroutine timestep
  
end module mod_timestep

