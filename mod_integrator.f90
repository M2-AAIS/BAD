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

    real (kind=x_precision)                    :: overdx, overdx2, overdt
    real (kind=x_precision), dimension(n_cell) :: dS_over_dt, dS_over_x_over_dx, S_over_x, dT_over_dx, nuS

    overdx  = 1/dx
    overdx2 = overdx**2
    overdt  = 1/dt
    nuS = s%nu*s%S

    ! Compute dT/dx as the mean spatial derivative (left + right)/2
    dT_over_dx(1)            = 0 ! FIXME
    dT_over_dx(2:n_cell - 1) = (s%T(3:n_cell) - s%T(1:n_cell - 2)) * overdx / 2
    dT_over_dx(n_cell)       = 0

    ! Compute d(S/x)/dx as the mean spatial derivative (left + right)/2
    S_over_x                        = s%S / s%x
    dS_over_x_over_dx(1)            = -1._x_precision / s%v(1) / s%x(1)**2 * nuS(2) / dx**2
    dS_over_x_over_dx(2:n_cell - 1) = (S_over_x(3:n_cell) - S_over_x(1:n_cell - 2)) * overdx / 2
    dS_over_x_over_dx(n_cell)       = 1._x_precision/dx**2 * & 
    ( 2._x_precision * dx + nuS(n_cell-1) - 2._x_precision * nuS(n_cell) + nuS(n_cell-1) )

    ! Compute dS/dt using the corresponding equation, with d²(nuS)/dx² as (d(left)+d(right))/2
    nuS = s%nu*s%S
    dS_over_dt(1)          = 0 ! FIXME
    dS_over_dt(2:n_cell-1) = 1/s%x(2:n_cell-1)**2 * (nuS(3:n_cell) - 2_x_precision*nuS(2:n_cell-1) + nuS(1:n_cell-2)) * overdx2
    dS_over_dt(n_cell)     = 0 ! FIXME 

    ! Second member of the dT/dt equation
    f = (3_x_precision * state_0%v_0**2 * s%nu * s%Omega**2 - &
         s%Fz * s%x / s%S &
         + state_0%T_0 * kmp/mu * (4._x_precision - 3._x_precision * s%beta) / s%beta * s%T / s%S * &
         (dS_over_dt + s%v * dS_over_x_over_dx) &
         - s%Cv * s%v / s%x * dT_over_dx) / s%Cv
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
    
    
    ! Create the diagonals
    diag     =  overdt + 2*overdx2 * states%nu / states%x**2
    diag_up  = -overdx2 * states%nu(2:n_cell) / states%x(1:n_cell-1)**2
    diag_low = -overdx2 * states%nu(1:n_cell-1) / states%x(2:n_cell)**2

    diag(1)            = 1_x_precision
    diag_up(1)         = 0
    diag(n_cell)       = states%nu(n_cell)
    diag_low(n_cell-1) = states%nu(n_cell-1)
    
    ! Solving for S, this does modify states%S directly
    call dgtsv(n_cell, 1, diag_low, diag, diag_up, states%S, n_cell, info)
    
    if (info /= 0) then
       print *, "Ooops, something bad happened!"
    end if
    
    dtemp = 0.01_x_precision * states%T
    T_deriv = states%T + dtemp
    !call compute_variables(states%x, states%Omega, T_deriv, states%S, states_deriv)
    !modification of compute_varaible--just enter a state now -- FIXME

    f0 = f(states, dt, dx)
    f1 = (f(states_deriv, dt, dx) - f(states, dt, dx)) / dtemp
    
    ! Solve for T
    states%T = states%T + f0 / (1_x_precision / dt - f1)

    ! Copy the result
  end subroutine do_timestep
  
end module mod_integrator

