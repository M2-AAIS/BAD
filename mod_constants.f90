! Module that exposes the integrator that transforms the state_in to thet state_out

module mod_constants
  integer, parameter, public :: x_precision = selected_real_kind(15)
  integer, parameter, public :: n_cell = 100
  integer, parameter, public :: n_iterations = 1000

  type parameters
     real(kind = x_precision) :: M, Mdot, X, Y, alpha, rmax
  end type parameters

  type state
     real(kind = x_precision) :: Omega, x, S, nu, v, T, L
  end type state
  
contains
    
end module mod_constants

