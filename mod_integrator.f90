! Module that exposes the integrator that transforms the state_in to the state_out

module mod_integrator
  use mod_constants
  use mod_timestep
  use mod_read_parameters

  implicit none 
  private

  public :: do_timestep
  
contains

  function f (s, dt, dx)
    use mod_constants
    use mod_variables
    implicit none
    
    type(state), intent(in)                    :: s
    real (kind=x_precision), intent(in)        :: dt, dx
    real (kind=x_precision), dimension(n_cell) :: f

    real (kind=x_precision)                    :: overdx, overdt
    real (kind=x_precision), dimension(n_cell) :: dS_over_dt, dS_over_x_over_dx, S_over_x, dT_over_dx

    overdx = 1/dx
    overdt = 1/dt

    dT_over_dx(1:n_cell - 1) = (s%T(2:n_cell) - s%T(1:n_cell - 1)) * overdx
    dT_over_dx(n_cell)       = 0

    S_over_x                        = s%S / s%x
    dS_over_x_over_dx(1:n_cell - 1) = (S_over_x(2:n_cell) - S_over_x(1:n_cell-1)) * overdx
    dS_over_x_over_dx(n_cell)       = 0 !TODO

    dS_over_dt  = 0 !TODO

    dS_over_x_over_dx(1:n_cell - 1) = (S_over_x(2:n_cell) - S_over_x(1:n_cell-1)) * overdx
    dS_over_dt = 0
    
    f = (3_x_precision*state_0%v_0**2/state_0%T_0 * s%nu * s%Omega**2 - &
         s%Fz*s%x/s%S &
         + kmp/mu * (4._x_precision-3._x_precision*s%beta)/s%beta * s%T/s%S * &
         (dS_over_dt + s%v*dS_over_x_over_dx) &
         - s%cv*s%v/s%x * dT_over_dx) / s%cv
  end function f

  
  subroutine do_timestep (states, param)
    use mod_variables
    use mod_constants
    implicit none
    
    type (state), intent(inout)                   :: states
    type (parameters), intent(in)                 :: param

    type (state)                                  :: states_deriv
    real(kind = x_precision)                      :: dt, overdx, overdx2, overdt, dx
    integer                                       :: info

    real(kind = x_precision), dimension(n_cell)   :: dtemp, T_deriv, diag
    real(kind = x_precision), dimension(n_cell-1) :: diag_low, diag_up
    real(kind = x_precision), dimension(n_cell)   :: f0, f1
    
    ! Get the timestep and the spacestep
    call timestep(states, dt)
    call process_dx(dx)

    overdx  = 1_x_precision/dx
    overdx2 = overdx**2
    overdt  = 1_x_precision/dt
    
    
    ! Create the diagonals and copy T and S
    diag     =  overdt + 2*overdx2 * states%nu / states%x**2
    diag_up  = -overdx2 * states%nu(2:n_cell) / states%x(1:n_cell-1)**2
    diag_low = -overdx2 * states%nu(1:n_cell-1) / states%x(2:n_cell)**2

    diag(1)            = 1_x_precision
    diag_up(1)         = 0
    diag(n_cell)       = states%nu(n_cell)
    diag_low(n_cell-1) = states%nu(n_cell-1)
    
    ! Solving for S
    call dgtsv(n_cell, 1, diag_low, diag, diag_up, states%S, n_cell, info)
    
    if (info /= 0) then
       print *, "Ooops, something bad happened!"
    end if
    
    call compute_variables(states%x, states%Omega, T_deriv, states%S, states_deriv)
    f0 = f(states, dt, dx)
    dtemp = 0.01_x_precision * states%T
    T_deriv = states%T + dtemp

    f1 = (f(states_deriv, dt, dx) - f(states, dt, dx)) / dtemp
    
    ! Solve for T
    states%T = states%T + f0 / (1_x_precision / dt - f1)

    ! Copy the result
  end subroutine do_timestep
  
end module mod_integrator

