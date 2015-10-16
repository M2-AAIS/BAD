! Module that exposes the integrator that transforms the state_in to the state_out

module mod_integrator
  
  private

  public :: do_timestep
  
contains
  
  subroutine do_timestep (state)
    use mod_constants
    use mod_timestep
    implicit none
    
    type (state), intent(in), dimension(n_cell)   :: state

    real(kind = x_precision)                      :: dt, dx2, dx = 0
    integer                                       :: i, info

    real(kind = x_precision), dimension(n_cell)   :: S, T, nu, diag
    real(kind = x_precision), dimension(n_cell-1) :: diag_low, diag_up
    
    ! Get the timestep
    call timestep(state, dt)

    dx2 = dx**2
    
    ! Create the diagonals and copy T and S
    diag     = 1/dt + 2/dx2 * state(:)%nu / state(:)%x
    diag_up  = -1/dx2 * state(2:n_cell)%nu / state(:)%x**2
    diag_low = -1/dx2 * state(1:n_cell-1)%nu / state(:)%x**2

    diag_up(1)           = -1/dx2 * state(2)%nu / state(1)%x**2
    diag(n_cell)         = state(n_cell)%nu
    diag_low(n_cell - 1) = state(n_cell-1)%nu
    
    ! Solving for S
    call dgtsv(n_cell, 1, diag_low, diag, diag_up, state_out(:)%S, n_cell, info)
    
    if (info /= 0) then
       print *, "Ooops, something bad happened!"
    end if

    ! Solve for T
    ! TODO

    ! Copy the result
  end subroutine do_timestep
  
end module mod_integrator

