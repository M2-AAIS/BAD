! Module that exposes the integrator that transforms the state_in to thet state_out

module mod_constants
  integer, parameter, public :: x_precision = selected_real_kind(15)
  integer, parameter, public :: n_cell = 100
  integer, parameter, public :: n_iterations = 1000

  type parameters
     real(kind = x_precision) :: M, Mdot, X, mu, alpha, rmax
     !M     : Black hole Mass
     !Mdot  : Acretion rate at rmax
     !X     : Chemical composition
     !mu    : Atomic mass
     !alpha : Viscosity parameter
     !rmax  : Maximal radius of accretion disk
  end type parameters

  type state
     real(kind = x_precision) :: Omega, x, nu, v, T, P, beta, cs, H, rho, sigma, Fz
     !Omega : Angular velocity
     !x     : Space variable
     !nu    : Viscosity
     !v     : Local speed accretion
     !T     : Temperature
     !P     : Pression
     !beta  : Pressure indicator
     !cs    : Speed of sound
     !H     : Half height of disk
     !rho   : Volume density
     !sigma : Surface density
     !Fz    : Radiative Flux
  end type state
  
contains
    
end module mod_constants

