! Module that exposes all the necessary constants, parameters and variable definitions

module mod_constants
  implicit none

  integer, parameter, public :: x_precision  = selected_real_kind(15) ! Numerical precision
  integer, parameter, public :: n_cell       = 256                    ! Number of cells along the disk radius
  integer, parameter, public :: n_iterations = 1000                   ! Number of iterations to run

  real(kind = x_precision), parameter, public :: G       = 6.67408e-8_x_precision ! Gravitationnal cst in cgs
  real(kind = x_precision), parameter, public :: c       = 2.99792458e10_x_precision ! Speed of light in cgs
  real(kind = x_precision), parameter, public :: pi      = 4.0_x_precision*atan(1.0_x_precision)
  real(kind = x_precision), parameter, public :: M_sun   = 1.98855e33_x_precision ! Mass of the sun in cgs
  real(kind = x_precision), parameter, public :: cst_rad = 7.5657308531642009e-17_x_precision ! Radiation cst in cgs
  real(kind = x_precision), parameter, public :: stefan  = (c * cst_rad) / 4.0_x_precision ! Stefan cst in cgs
  real(kind = x_precision), parameter, public :: kb      = 1.3806488e-16_x_precision ! Boltzmann cst in cgs
  real(kind = x_precision), parameter, public :: R       = 8.3144598e7_x_precision ! Gas cst in csg = Boltzmann cst over proton mass
  real(kind = x_precision), parameter, public :: gammag  = 5._x_precision / 3._x_precision ! Adiabatic coefficient
  real(kind = x_precision), parameter, public :: thomson = 6.6524587158e-25_x_precision ! Thomson cross-section in cgs

  type parameters
    real(kind = x_precision) :: M, Mdot, X, RTM, alpha, dx
    ! M     : Black hole Mass
    ! Mdot  : Acretion rate at rmax
    ! X     : Chemical composition in H
    ! RTM   : Gas constant × T_0 ÷ μ
    ! alpha : Viscosity parameter
    ! dx    : Dimensionless spacestep
  end type parameters

  type state
    real(kind = x_precision), dimension(n_cell) :: x, Omega, nu, v, cs, S, H, M_dot, rho, T, Fz, Cv, P_gaz, P_rad, beta
    ! x     : Space variable
    ! Omega : Angular velocity
    ! nu    : Viscosity
    ! v     : Local, radial accretion speed
    ! cs    : Speed of sound
    ! S     : Variable of density
    ! H     : Disk half-height
    ! M_dot : Accretion rate
    ! rho   : Volume density
    ! T     : Temperature
    ! Fz    : Radiative flux
    ! Cv    : Heat capacity at constant volume
    ! P_gaz : Gaz pressure
    ! P_rad : Radiative pressure
    ! beta  : Pressure indicator
  end type state

  type adim_state
    real(kind= x_precision) :: temps_0, Omega_0, nu_0, v_0, cs_0, S_0, H_0, M_dot_0, rho_0, T_0, Fz_0, Cv_0, P_gaz_0, P_rad_0
    ! Usefull quantities to get the dimensionless ones
  end type adim_state

contains

end module mod_constants

