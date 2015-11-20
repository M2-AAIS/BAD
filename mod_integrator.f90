! Module that exposes the integrator that transforms the state_in to the state_out

module mod_integrator
  use mod_constants
  use mod_timestep
  use mod_read_parameters
  use mod_variables

  implicit none 
  private

  public :: do_timestep_T, do_timestep_S
  
contains
  ! (T(t+dt) - T(t))/dt = f(s, dt)
  function f (s, dt)
  !process the right term of \partial T* / \partial t* = (see Recapitulatif des adimensionenemts in report)
    implicit none
    
    type(state), intent(in)                    :: s
    real (kind=x_precision), intent(in)        :: dt
    real (kind=x_precision), dimension(n_cell) :: f

    real (kind=x_precision)                    :: overdx, overdx2, overdt
    real (kind=x_precision), dimension(n_cell) :: dS_over_dt, dS_over_x_over_dx, S_over_x, dT_over_dx, nuS

    overdx  = 1/params%dx
    overdx2 = overdx**2
    overdt  = 1/dt
    nuS = s%nu*s%S

    ! Compute dT/dx as the mean spatial derivative (left + right)/2
    dT_over_dx(1:n_cell - 1) = (s%T(2:n_cell) - s%T(1:n_cell - 1)) * overdx
    dT_over_dx(n_cell)       = 0

    ! Compute d(S/x)/dx as the mean spatial derivative (left + right)/2
    S_over_x                        = s%S / x_state%x
    dS_over_x_over_dx(1:n_cell - 1) = (S_over_x(2:n_cell) - S_over_x(1:n_cell - 1)) * overdx
    dS_over_x_over_dx(n_cell)       = ( params%dx - nuS(n_cell) + nuS(n_cell-1) ) / &
                                      ( x_state%x(n_cell)**2 * s%v(n_cell) * params%dx**2 )

    ! Compute dS/dt using the corresponding equation, with d²(nuS)/dx² as (d(left)+d(right))/2
    dS_over_dt(1:n_cell-1) = 1 / x_state%x(1:n_cell-1)**2 * (nuS(2:n_cell) - 2._x_precision*nuS(1:n_cell-1) + &
                             nuS(1:n_cell-1)) * overdx2
    dS_over_dt(n_cell)     = ( params%dx - nuS(n_cell) + nuS(n_cell-1) ) / &
                             ( x_state%x(n_cell)**2 * params%dx**2 )

    ! Second member of the dT/dt equation
    f = (3_x_precision * state_0%v_0**2 * s%nu * x_state%Omega**2 - &
         s%Fz * x_state%x / s%S &
         + params%RTM * (4._x_precision - 3._x_precision * s%beta) / s%beta * s%T / s%S * &
         (dS_over_dt + s%v * dS_over_x_over_dx) &
         - s%Cv * s%v / x_state%x * dT_over_dx) / s%Cv
  end function f

  subroutine do_timestep_S (states, dt)
  !process the temporal evolution of S
    implicit none
    
    type (state), intent(inout)                   :: states

    real(kind = x_precision)                      :: overdx, overdx2, overdt
    real(kind = x_precision), intent(in)          :: dt
    integer                                       :: info

    real(kind = x_precision), dimension(n_cell)   :: diag
    real(kind = x_precision), dimension(n_cell-1) :: diag_low, diag_up
    
    overdx  = 1_x_precision/params%dx
    overdx2 = overdx**2
    overdt  = 1_x_precision/dt
    
    
    ! Create the diagonals
    diag     =  overdt + 2*overdx2 * states%nu / x_state%x**2
    diag_up  = -overdx2 * states%nu(2:n_cell) / x_state%x(1:n_cell-1)**2
    diag_low = -overdx2 * states%nu(1:n_cell-1) / x_state%x(2:n_cell)**2

    diag(1)            = 1_x_precision
    diag_up(1)         = 0
    diag(n_cell)       = states%nu(n_cell)
    diag_low(n_cell-1) = states%nu(n_cell-1)
    
    ! Solving for S, this does modify states%S directly
    call dgtsv(n_cell, 1, diag_low, diag, diag_up, states%S, n_cell, info)
    
    if (info /= 0) then
      stop "error in DGTSV call in subroutine do_timestep_S"
    end if
  end subroutine do_timestep_S

  subroutine do_timestep_S_exp (states, dt)
    implicit none
    
    type (state), intent(inout)                   :: states
    real (kind = x_precision), intent(in)         :: dt
    
    real (kind = x_precision), dimension(n_cell)  :: S_tmp

    S_tmp(2:n_cell-1) = states%S(2:n_cell-1) + ((states%nu(3:n_cell) * states%S(3:n_cell)) &
         - 2._x_precision * (states%nu(2:n_cell-1) * states%S(2:n_cell-1)) + ( states%nu(1:n_cell-2) &
         * states%S(1:n_cell-2))) / (params%dx**2 *  x_state%x(2:n_cell-1)**2) * dt

    S_tmp(1)          = states%S(1) + (states%nu(2) * states%S(2) - (2._x_precision * states%nu(1) &
         * states%S(1)) ) /  (params%dx**2 *  x_state%x(1)**2) * dt

    S_tmp(n_cell)     = states%S(n_cell) + (params%dx - states%nu(n_cell) * states%nu(n_cell) &
         + states%nu(n_cell-1) * states%S(n_cell-1)) / params%dx**2 * dt

    states%S = S_tmp
  end subroutine do_timestep_S_exp
  
  subroutine do_timestep_T(states, dt) 
  !process the temporal evolution of T
    implicit none
    
    type (state), intent(inout)                   :: states

    type (state)                                  :: states_deriv
    real(kind = x_precision), intent(in)          :: dt

    real(kind = x_precision), dimension(n_cell)   :: dtemp
    real(kind = x_precision), dimension(n_cell)   :: f0, f1

    !call timestep(states, dt) !FIXME call of function to process the viscosity timestep
    dtemp = 0.01_x_precision * states%T
    states_deriv = states
    states_deriv%T = states%T + dtemp
    
    call compute_variables(states_deriv)

    ! Let dT/dt = f0 and d²T/dt² = f1 so T(t+1) = T(t) + dt * f0 * (1 + dt * f1)
    
    f0 = f(states, dt)
    f1 = (f(states_deriv, dt) - f(states, dt)) / dtemp
    
    states%T = states%T + dt * f0 * (1_x_precision + dt * f1)

    ! Copy the result
  end subroutine do_timestep_T
  
end module mod_integrator
