! Module that exposes the initial (physical) parameters after parsing them from a file

module mod_read_parameters
  use mod_constants

  implicit none

  type(parameters) :: params  ! Contains the initial parameters of the black hole accretion disk (see mod_constants.f90 for details)
  type(adim_state) :: x_state ! Contains the x and Omega tables
  type(dim_state)  :: r_state ! Contains the r and Omega_r tables

  private

  public           :: get_parameters, params, x_state, r_state

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
    real(kind=x_precision) :: T_0     ! Order of magnitude for temperature
    real(kind=x_precision) :: Omega_0 ! Order of magnitude for angular velocity
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
    params%M    = params%M * M_sun
    params%Mdot_crit =  12._x_precision * Ledd * params%M / c**2 
    params%Mdot = params%Mdot *  params%Mdot_crit

    ! Compute c2, rs, T_0 and Omega_0
    c2      = c**2
    rs      = 2._x_precision * G * params%M / c2
    T_0     = (params%Mdot * c2 / (sqrt(27.0) * 48._x_precision * pi * rs**2 * stefan))**0.25_x_precision
    Omega_0 = sqrt(G * params%M / rs**3)

    ! Compute mu
    Z  = 1._x_precision - X - Y
    mu = 1._x_precision / (2._x_precision*X + 3._x_precision*Y/4._x_precision + Z/2._x_precision)

    ! Process RTM
    params%RTM = R * T_0 / mu

    ! Process dx
    rmin      = 3._x_precision
    params%dx = (sqrt(rmax) - sqrt(rmin)) / n_cell

    
    ! Process x_state, r_state
    do i = 1, n_cell
      x_state%x(i)       = sqrt(rmin) + i * params%dx
      r_state%r(i)       = x_state%x(i)**2 * rs
      x_state%Omega(i)   = 1._x_precision / x_state%x(i)**3
      r_state%Omega_r(i) = x_state%Omega(i) * Omega_0
    end do

    ! Process t_v
    t_v = 0.28 * params%alpha**(-4._x_precision / 5._x_precision) * (M_dot/M_dot_c)**(-3._x_precision/10._x_precision)*(M/M_sun)**(6._x_precision/5._x_precision)*(r/(3*rs))**(5._x_precision / 4._x_precision)*( 1 - (r_state%r(i)/(3 * rs))**(-1._x_precision/2._x_precision))**(-6._x_precision / 5._x_precision)
    
    ! Process kappa_e
    params%kappa_e = 0.2_x_precision * (1._x_precision + X)

    ! Display parameters
    write(*,*)'           Input Parameters             '
    write(*,*)'****************************************'
    write(*,"(' BH_mass     =',1p,E12.4)") params%M
    write(*,"(' rmax        =',1p,E12.4)") rmax
    write(*,"(' Mdot        =',1p,E12.4)") params%Mdot
    write(*,"(' alpha       =',1p,E12.4)") params%alpha
    write(*,"(' X           =',1p,E12.4)") X
    write(*,"(' Y           =',1p,E12.4)") Y
    write(*,*)'****************************************'
    read(*,*)
  end subroutine get_parameters

end module mod_read_parameters
