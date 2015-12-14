! Module that exposes all the necessary constants, parameters and variable definitions

module mod_constants
  implicit none

  integer, parameter, public :: x_precision  = selected_real_kind(15) ! Numerical precision
  integer, parameter, public :: n_cell       = 256                    ! Number of cells along the disk radius
  ! Specific to main program
  integer, parameter, public :: n_iterations = 10000000               ! Number of iterations to run
  integer, parameter, public :: output_freq  = 8000                   ! Frequency of the outputs
  real(x_precision), parameter, public :: cst_dt_nu = 100._x_precision
  real(x_precision), parameter, public :: cst_dt_T  = 100._x_precision
  ! Specific to s_curve
  integer,           parameter, public :: nb_it  = 1000              ! Number of points for the S curve
  integer,           parameter, public :: max_it = 10000000          ! Maximum number of s_curve dichotomy iterations
  real(x_precision), parameter, public :: eps    = 1.e-6_x_precision ! Precision required for the dichotomy
  real(x_precision), parameter, public :: Tmin   = 10._x_precision**5.3_x_precision
  real(x_precision), parameter, public :: Tmax   = 10._x_precision**7.4_x_precision
  real(x_precision), parameter, public :: Smin   = 10._x_precision**1.0_x_precision
  real(x_precision), parameter, public :: Smax   = 10._x_precision**3.5_x_precision

  ! Physical constants
  real(x_precision), parameter, public :: G       = 6.67408e-8_x_precision ! Gravitationnal cst in cgs
  real(x_precision), parameter, public :: c       = 2.99792458e10_x_precision ! Speed of light in cgs
  real(x_precision), parameter, public :: pi      = 4.0_x_precision*atan(1.0_x_precision)
  real(x_precision), parameter, public :: M_sun   = 1.98855e33_x_precision ! Mass of the sun in cgs
  real(x_precision), parameter, public :: cst_rad = 7.5657308531642009e-15_x_precision ! Radiation cst in cgs
  real(x_precision), parameter, public :: stefan  = (c * cst_rad) / 4.0_x_precision ! Stefan cst in cgs
  real(x_precision), parameter, public :: kb      = 1.3806488e-16_x_precision ! Boltzmann cst in cgs
  real(x_precision), parameter, public :: R       = 8.3144598e7_x_precision ! Gas cst in csg = Boltzmann cst over proton mass
  real(x_precision), parameter, public :: gammag  = 5._x_precision / 3._x_precision ! Adiabatic coefficient
  real(x_precision), parameter, public :: thomson = 6.6524587158e-25_x_precision ! Thomson cross-section in cgs

  type parameters
    real(x_precision) :: M, Mdot, Mdot_crit, kappa_e, RTM, alpha, dx, t_nu, t_T
    ! M         : Black hole Mass
    ! Mdot      : Acretion rate at rmax
    ! Mdot_crit : Critical accretion rate
    ! kappa_e   : Thomson scattering thickness
    ! RTM       : Gas constant × T_0 ÷ μ
    ! alpha     : Viscosity parameter
    ! dx        : Dimensionless spacestep
    ! t_nu      : Viscous time
    ! t_T       : Thermal time
  end type parameters

  type state
    real(x_precision), dimension(n_cell) :: nu, v, cs, S, H, Mdot, rho, T, Fz, Cv, Pgaz, Prad, beta, Qplus, Qminus, tau, Qadv
    ! nu     : Viscosity
    ! v      : Local, radial accretion speed
    ! cs     : Speed of sound
    ! S      : Variable of density
    ! H      : Disk half-height
    ! Mdot   : Accretion rate
    ! rho    : Volume density
    ! T      : Temperature
    ! Fz     : Radiative flux
    ! Cv     : Heat capacity at constant volume
    ! Pgaz   : Gaz pressure
    ! Prad   : Radiative pressure
    ! beta   : Pressure indicator
    ! Qplus  : Heating term
    ! Qminus : Cooling term
    ! tau    : Opacity
    ! Qadv   : Advection term
  end type state

  type dim_state
    real(x_precision), dimension(n_cell) :: r, Omega_r
    ! r       : Space variable
    ! Omega_r : Angular velocity
  end type dim_state

  type adim_state
    real(x_precision), dimension(n_cell) :: x, Omega
    ! x     : Space variable
    ! Omega : Angular velocity
  end type adim_state

  type state_zero
    real(x_precision) :: temps_0, Omega_0, nu_0, v_0, cs_0, S_0, H_0, Mdot_0, rho_0, T_0, Fz_0, Cv_0, Pgaz_0, Prad_0, beta_0
    ! Usefull quantities to get the dimensionless ones
  end type state_zero

  type state_ic
     real(x_precision), dimension(n_cell) :: T, Sigma, H
     ! T     : Initial condition for Temperature
     ! Sigma : Initial condition for Sigma
     ! H     : Initial condition for H (for comparison only)
  end type state_ic

contains

end module mod_constants
