! Module that exposes the integrator that transforms the state_in to the state_out

module mod_integrator
  use mod_constants
  use mod_timestep
  use mod_read_parameters

  implicit none 
  private

  public :: do_timestep_T, do_timestep_S
  
contains
  ! (T(t+dt) - T(t))/dt = f(s, dt, dx)
  function f (s, dt, dx)
!process the right term of \partial T* / \partial t* = (see Recapitulatif des adimensionenemts in report)
    use mod_constants
    use mod_variables
    implicit none
    
    type(state), intent(in)                    :: s
    type(parameters)                           :: para
    real (kind=x_precision), intent(in)        :: dt, dx
    real (kind=x_precision), dimension(n_cell) :: f

    real (kind=x_precision)                    :: overdx, overdx2, overdt
    real (kind=x_precision), dimension(n_cell) :: dS_over_dt, dS_over_x_over_dx, S_over_x, dT_over_dx, nuS

    call get_parameters(para)

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
    dS_over_x_over_dx(1)            = -2._x_precision / s%v(1) / s%x(1)**2 * nuS(2) / dx**2
    dS_over_x_over_dx(2:n_cell - 1) = (S_over_x(3:n_cell) - S_over_x(1:n_cell - 2)) * overdx / 2
    dS_over_x_over_dx(n_cell)       = -2._x_precision / s%x(n_cell)**2 / s%v(n_cell) / dx**2 * &
    ( dx + nuS(n_cell-1) - nuS(n_cell) )

    ! Compute dS/dt using the corresponding equation, with d²(nuS)/dx² as (d(left)+d(right))/2
    dS_over_dt(1)          = 2._x_precision * nuS(2) / s%x(1)**2 / dx**2
    dS_over_dt(2:n_cell-1) = 1/s%x(2:n_cell-1)**2 * &
    (nuS(3:n_cell) - 2_x_precision*nuS(2:n_cell-1) + nuS(1:n_cell-2)) * overdx2
    dS_over_dt(n_cell)     = 2._x_precision / s%x(n_cell)**2 / dx**2 * &
    ( dx + nuS(n_cell-1) - nuS(n_cell) )

    ! Second member of the dT/dt equation
    f = (3_x_precision * state_0%v_0**2 * s%nu * s%Omega**2 - &
         s%Fz * s%x / s%S &
         + para%RTM * (4._x_precision - 3._x_precision * s%beta) / s%beta * s%T / s%S * &
         (dS_over_dt + s%v * dS_over_x_over_dx) &
         - s%Cv * s%v / s%x * dT_over_dx) / s%Cv
  end function f


  subroutine do_timestep_S (states)
  !process the temporal evolution of S
    use mod_variables
    use mod_constants
    implicit none
    
    type (state), intent(inout)                   :: states

    real(kind = x_precision)                      :: dt, overdx, overdx2, overdt, dx
    integer                                       :: info

    real(kind = x_precision), dimension(n_cell)   :: dtemp, diag
    real(kind = x_precision), dimension(n_cell-1) :: diag_low, diag_up
    
    ! Get the timestep and the spacestep
    !call timestep(states, dt) !FIXME call of function to process the viscosity timestep
    dt=1d-6
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
      stop "error in DGTSV call in subroutine do_timestep_S"
    end if
  end subroutine do_timestep_S 


  subroutine do_timestep_T(states) 
  !process the temporal evolution of T
    use mod_variables
    use mod_constants
    implicit none
    
    type (state), intent(inout)                   :: states

    type (state)                                  :: states_deriv
    real(kind = x_precision)                      :: dt, dx

    real(kind = x_precision), dimension(n_cell)   :: dtemp, T_deriv 
    real(kind = x_precision), dimension(n_cell)   :: f0, f1

    !call timestep(states, dt) !FIXME call of function to process the viscosity timestep
    call process_dx(dx)
    dt=1d-6

    dtemp = 0.01_x_precision * states%T
    states_deriv = states
    states_deriv%T = states%T + dtemp
    
    call compute_variables(states_deriv)
    !modification of compute_varaible--just enter a state now -- FIXME

    ! Let dT/dt = f0 and d²T/dt² = f1 so T(t+1) = T(t) + dt * f0 * (1 + dt * f1)
    
    f0 = f(states, dt, dx)
    f1 = (f(states_deriv, dt, dx) - f(states, dt, dx)) / dtemp
    
    states%T = states%T + dt * f0 * (1_x_precision + dt * f1)

    ! Copy the result
  end subroutine do_timestep_T
  
end module mod_integrator

