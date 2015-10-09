! Module that exposes the integrator that transforms the state_in to thet state_out

module mod_integrator
  
  private

  public :: do_timestep
  
contains
  
  subroutine do_timestep (state_in, state_out)
    use mod_constants
    implicit none
    type (state), intent(in), dimension(:)  :: state_in
    type (state), intent(out), dimension(:) :: state_out
    
  end subroutine do_timestep
  
end module mod_integrator

