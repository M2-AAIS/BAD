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
    use mod_variables
    implicit none
    
    type (state), intent(in)                      :: states

    real(kind = x_precision)                      :: dt, overdx, overdx2, overdt, dx, B
    integer                                       :: info

    real(kind = x_precision), dimension(n_cell)   :: diag, dT_over_dx, S_over_x, dS_over_x_over_dx
    real(kind = x_precision), dimension(n_cell)   :: dS_over_dt, dH_over_dT, f0, f1, dFz_over_dT
    real(kind = x_precision), dimension(n_cell-1) :: diag_low, diag_up
    
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

    dT_over_dx(1:n_cell - 1) = (states%T(2:n_cell) - states%T(1:n_cell - 1)) * overdx
    dT_over_dx(n_cell) = 0

    S_over_x = states%S / states%x
    dS_over_x_over_dx(1:n_cell - 1) = (S_over_x(2:n_cell) - S_over_x(1:n_cell-1)) * overdx
    dS_over_x_over_dx(n_cell) = 0 !TODO

    dS_over_dt = 0 !TODO
    dH_over_dT = 0 !TODO

    dFz_over_dT = 0 !TODO

    B = state_0%P_rad / state_0%P_gaz
    
    f0 = 3*v0**2/T0 * states%nu * states%Omega**2 - &
         states%Fz*states%x/states%S &
         + R/mu * (4-3*states%beta)/states%beta * states%T/states%S * ((states%S - prevS)/dt + states%v*dS_over_x_over_dx) &
         - cv*states%v/states%x * dT_over_dx

    f1 = 3*v0**2/T0 * (2*param%alpha * states%H * dH_over_dT * states%Omega**3) + &
         dFz_over_dT * states%x / states%S + &
         R/mu * (1 - 12*B*states%T**3 / states%rho) * (dS_over_dt + states%v * dS_over_x_over_dx )
    !           \____ (4-3beta)/beta * T*/S* _____/

    ! Solve for T
    states%T = states%T + f0 / (param%cv / dt - f1)

    ! Copy the result
  end subroutine do_timestep
  
end module mod_integrator

