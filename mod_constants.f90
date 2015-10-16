! Module that exposes the integrator that transforms the state_in to thet state_out

module mod_constants
  integer, parameter, public :: x_precision = selected_real_kind(15)
  integer, parameter, public :: n_cell = 256
  integer, parameter, public :: n_iterations = 1000

  !Unit of constants : cgs
  real(kind = x_precision), parameter, public :: G = 6.67384e-8 
  real(kind = x_precision), parameter, public :: c = 2.99792458e10
  real(kind = x_precision), parameter, public :: pi = 4.*datan(1.d0)
  real(kind = x_precision), parameter, public :: stefan = 5.66956e-5
  real(kind = x_precision), parameter, public :: kb_mp = 8.31434e7   !kb/mp
  real(kind = x_precision), parameter, public :: a = 7.564e-15
  real(kind = x_precision), parameter, public :: mu = 0.62
    

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
     real(kind = x_precision), dimension(n_cell) :: Omega, x, nu, v, T, P_rad, P_gaz, beta, cs, H, rho, S, Fz, M_dot
     !Omega : Angular velocity
     !x     : Space variable
     !nu    : Viscosity
     !v     : Local speed accretion
     !T     : Temperature
     !P_rad : Radiative pressure
     !P_gaz : Gaz pressure
     !beta  : Pressure indicator
     !cs    : Speed of sound
     !H     : Half height of disk
     !rho   : Volume density
     !S     : Variable of density
     !Fz    : Radiative Flux
     !M_dot : acretion rate
     !s     : S
  end type state

  type adim_state
     real(kind= x_precision)         :: Omega, x, nu, v, T, P_rad, P_gaz, cs, H, rho, S, M_dot
  end type adim_state
  
contains
    
end module mod_constants

