! Module that exposes the integrator that transforms the state_in to the state_out

module mod_integrator
  use mod_constants
  use mod_timestep
  use mod_read_parameters

  implicit none 
  private

  public :: do_timestep
  
contains
  
  subroutine do_timestep (states)
   implicit none
    
    type (state), intent(in)                      :: states

    real(kind = x_precision)                      :: dt, dx2, dx != 0
    integer                                       :: i, info

    real(kind = x_precision), dimension(n_cell)   :: S, T, nu, diag
    real(kind = x_precision), dimension(n_cell-1) :: diag_low, diag_up
    
    ! Get the timestep and the spacestep
    call timestep(states, dt)
    call process_dx(dx)

    dx2 = dx**2
    
    ! Create the diagonals and copy T and S
    diag      = 1._x_precision/dt + 2._x_precision/dx2 * states%nu / states%x
    diag_up  = -1._x_precision/dx2 * states%nu(2:n_cell) / states%x(1:n_cell-1)**2
    diag_low = -1._x_precision/dx2 * states%nu(1:n_cell-1) / states%x(2:n_cell)**2

    diag_up(1)           = -1/dx2 * states%nu(2) / states%x(1)**2
    diag(n_cell)         = states%nu(n_cell)
    diag_low(n_cell - 1) = states%nu(n_cell-1)
    
    ! Solving for S
    call dgtsv(n_cell, 1, diag_low, diag, diag_up, states%S(:), n_cell, info)
    
    if (info /= 0) then
       print *, "Ooops, something bad happened!"
    end if

    ! Solve for T
    ! TODO

    ! Copy the result
  end subroutine do_timestep
  
end module mod_integrator

