! Module that exposes the initial (physical) parameters after parsing them from a file

module mod_read_parameters
  use mod_constants

  implicit none

  type(parameters) :: params  ! Contains the initial parameters of the black hole accretion disk (see mod_constants.f90 for details)
  type(adim_state) :: x_state ! Contains the x and Omega tables
  type(dim_state)  :: r_state ! Contains the r and Omega_r tables
  type(state_zero) :: state_0
  type(state_ci)   :: CI      ! Contains the initial parmaters for Temparature end Sigma

  real(kind = x_precision), dimension(n_cell) :: f1 

  private

  public           :: get_parameters, params, x_state, r_state, state_0, f1, CI

contains

  ! Read parameter from the input file: ./input_parameter.dat
  subroutine get_parameters()
    implicit none

    ! Internal variables
    real(kind=x_precision) :: X,Y,Z   ! Chemical composition : X+Y+Z = 1
    real(kind=x_precision) :: mu      ! Mean molecular mass
    real(kind=x_precision) :: rmax    ! Maximum considered radius (in rs)
    real(kind=x_precision) :: rmin    ! Minimal radius (in rs)
    real(kind=x_precision) :: mp      ! Proton mass in cgs
    real(kind=x_precision) :: Ledd    ! Eddington luminosity
    real(kind=x_precision) :: c2      ! Light speed to the square
    real(kind=x_precision) :: rs      ! Schwarzschild radius
    real(kind=x_precision) :: t_dyn   ! Dynamic time
    integer(kind=4)        :: ios     ! I/O test variable
    integer                :: i       ! Cells iteration variable
    character(len=50)      :: line    ! String reading variable

    ! Open file
    open(unit=11, file="./input_parameter.dat", action="read", status="old", iostat=ios)
    if (ios /= 0) stop "OPENING input_parameter.dat ERROR"

    ! Read parameters
    read(11,fmt=*) line, params%M
    read(11,fmt=*) line, rmax
    read(11,fmt=*) line, params%Mdot
    read(11,fmt=*) line, params%alpha
    read(11,fmt=*) line, X
    read(11,fmt=*) line, Y

    ! Close the file
    close(11)

    ! Compute mp
    mp = kb / R

    ! Compute Ledd
    Ledd = 4._x_precision * pi * G * mp * c / thomson

    ! Process M and M_dot in cgs, M_dot_crit = 12 * Ledd * M / c^2
    params%M         = params%M * M_sun
    params%Mdot_crit =  12._x_precision * Ledd * params%M / c**2 
    params%Mdot      = params%Mdot * params%Mdot_crit

    ! Compute c2 and rs
    c2 = c**2
    rs = 2._x_precision * G * params%M / c2

    ! Compute T_0 and Omega_0
    state_0%T_0     = (params%Mdot * c2 / (sqrt(27.0) * 48._x_precision * pi * rs**2 * stefan))**0.25_x_precision
    state_0%Omega_0 = sqrt(G * params%M / rs**3)

    ! Compute mu
    Z  = 1._x_precision - X - Y
    mu = 1._x_precision / (2._x_precision*X + 3._x_precision*Y/4._x_precision + Z/2._x_precision)

    ! Process RTM
    params%RTM = R * state_0%T_0 / mu

    ! Process dx
    rmin      = 3._x_precision
    params%dx = (sqrt(rmax) - sqrt(rmin)) / n_cell

    ! Process x_state, r_state
    do i = 1, n_cell
      x_state%x(i)       = sqrt(rmin) + i * params%dx
      r_state%r(i)       = x_state%x(i)**2 * rs
      x_state%Omega(i)   = 1._x_precision / x_state%x(i)**3
      r_state%Omega_r(i) = x_state%Omega(i) * state_0%Omega_0
    end do

    ! Process t_nu
    params%t_nu = 0.28 * params%alpha**(-0.8_x_precision) * (params%Mdot / params%Mdot_crit)**(-.3_x_precision) * &
                  (params%M / M_sun)**1.2_x_precision * (r_state%r(n_cell) / (3*rs))**1.25_x_precision * &
                  (1 - (r_state%r(n_cell)/(3 * rs))**(-0.5_x_precision))**(-1.2_x_precision)

    ! Process t_dyn
    t_dyn       = 7.2 * 10**(-5._x_precision) * (params%M / M_sun) * (r_state%r(n_cell) / (3 * rs))**1.5_x_precision
    
    ! Process t_T
    params%t_T  = t_dyn / params%alpha
    
    !-----------------------------------------------------------
    !--Process of initial parmameter for Temperature and Sigma-- 
    !-----------------------------------------------------------
    ! Process f1 to compute T_ci and Sig_ci
    f1 = 1._x_precision - (sqrt(3._x_precision) / x_state%x)
    
    ! Process T_ci
    CI%T_ci = 1.4e4_x_precision * (params%alpha)**(-1._x_precision/5._x_precision) & 
         * (params%Mdot / 1.e16_x_precision)**(3._x_precision/10) * (params%M / M_sun)**(1._x_precision/4._x_precision) & 
         * (r_state%r / 1.e10_x_precision)**(-3._x_precision/4._x_precision) * f1**(3._x_precision/10._x_precision)
    ! Process Sig_ci 
    CI%Sig_ci = 5.2_x_precision * params%alpha**(-4._x_precision/5._x_precision) * &
         (params%Mdot / 1.e16_x_precision)**(7._x_precision/10._x_precision) * (params%M / M_sun)**(1._x_precision/4._x_precision) &
         * (r_state%r / 1.e10_x_precision)**(-3._x_precision/4._x_precision) * f1**(7._x_precision/10._x_precision)
    !-----------------------------------------------------------
    !-----------------------------------------------------------

    ! Process kappa_e
    params%kappa_e = 0.2_x_precision * (1._x_precision + X)
    
    ! Display parameters
    write(*,"('#           Input Parameters             ')")
    write(*,"('#****************************************')")
    write(*,"('# BH_mass     =',1p,E12.4)") params%M
    write(*,"('# rmax        =',1p,E12.4)") rmax
    write(*,"('# Mdot        =',1p,E12.4)") params%Mdot
    write(*,"('# alpha       =',1p,E12.4)") params%alpha
    write(*,"('# X           =',1p,E12.4)") X
    write(*,"('# Y           =',1p,E12.4)") Y
    write(*,"('# t_nu        =',1p,E12.4)") params%t_nu
    write(*,"('# t_T         =',1p,E12.4)") params%t_T
    write(*,"('#****************************************')")

    ! Compute the initial variables needed for dimensionless variables
    !state_0%Omega_0 = sqrt(G * params%M / rs**3)
    state_0%temps_0 = 2._x_precision / state_0%Omega_0
    state_0%nu_0    = 2._x_precision * state_0%Omega_0 * rs**2 / 3._x_precision
    state_0%v_0     = state_0%Omega_0 * rs
    state_0%cs_0    = state_0%v_0
    state_0%S_0     = params%Mdot / (3._x_precision * pi * state_0%nu_0)
    state_0%H_0     = rs
    state_0%Mdot_0  = params%Mdot
    state_0%rho_0   = state_0%S_0 / (2._x_precision * rs)
    !state_0%T_0     = (params%Mdot * c2 / (sqrt(27.0) * 48._x_precision * pi * rs**2 * stefan))**0.25_x_precision
    state_0%Fz_0    = state_0%S_0 * state_0%Omega_0 / 4._x_precision
    state_0%Cv_0    = 1._x_precision / state_0%T_0
    state_0%Pgaz_0  = state_0%rho_0 * params%RTM
    state_0%Prad_0  = cst_rad * state_0%T_0**4 / 3._x_precision
    state_0%beta_0  = state_0%Prad_0 / state_0%Pgaz_0

    write(*,"('#           Initial Variables            ')")
    write(*,"('#****************************************')")
    write(*,"('# Omega_0     =',1p,E12.4)") state_0%Omega_0
    write(*,"('# t_0         =',1p,E12.4)") state_0%temps_0
    write(*,"('# nu_0        =',1p,E12.4)") state_0%nu_0
    write(*,"('# v_0         =',1p,E12.4)") state_0%v_0
    write(*,"('# cs_0        =',1p,E12.4)") state_0%cs_0
    write(*,"('# Sigma_0     =',1p,E12.4)") state_0%S_0
    write(*,"('# H_0         =',1p,E12.4)") state_0%H_0
    write(*,"('# Mdot_0      =',1p,E12.4)") state_0%Mdot_0
    write(*,"('# rho_0       =',1p,E12.4)") state_0%rho_0
    write(*,"('# T_0         =',1p,E12.4)") state_0%T_0
    write(*,"('# Fz_0        =',1p,E12.4)") state_0%Fz_0
    write(*,"('# Cv_0        =',1p,E12.4)") state_0%Cv_0
    write(*,"('# P_gaz_0     =',1p,E12.4)") state_0%Pgaz_0
    write(*,"('# P_rad_0     =',1p,E12.4)") state_0%Prad_0
    write(*,"('#****************************************')")

  end subroutine get_parameters

end module mod_read_parameters
