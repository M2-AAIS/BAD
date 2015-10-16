! Module that exposes the integrator that transforms the state_in to the state_out

module mod_integrator
  
  private

  public :: do_timestep
  
contains
  
  subroutine do_timestep (state_in, state_out)
    use mod_constants
    use mod_timestep
    implicit none
    
    type (state), intent(in), dimension(:)  :: state_in
    type (state), intent(out), dimension(:) :: state_out

    real(kind = x_precision)                :: dt
    integer                                 :: i

    ! Get the timestep
    call timestep(state_in, dt)
    
    ! Do the integration
    do i = 0, n_iterations
       ! Integrate equation 1
       ! Integrate equation 2
       ! Integrate equation 3
    end do
    
  end subroutine do_timestep
  
end module mod_integrator

