! Module that reads the parameter out of a file

module mod_read_parameters
  use mod_constants
  implicit none
  private

  public :: get_parameters, process_dx

contains

  subroutine get_parameters (initial_param)
    ! Read parameter in the input file: ./input_parameter.dat

    implicit none

    ! Function parameters
    type(parameters), intent(out) :: initial_param ! Contains the initial parameters of black hole accretion disk (see mod_constants.f90 for details)

    ! Internal variables
    real (kind=x_precision)       :: Y,Z           ! Chemical composition : x+y+z = 1
    real (kind=x_precision)       :: mu            ! Mean molecular mass
    real (kind=x_precision)       :: rmax          ! Maximum considered radius (in rs)
    real (kind=x_precision)       :: Ledd          ! Eddington luminosity
    real (kind=x_precision)       :: rs            ! Schwarzschild radius
    real (kind=x_precision)       :: T_0           ! Temperature order of magnitude
    integer (kind=4)              :: ios           ! i/o variable test
    character (len=50)            :: bla           ! string reading variable

    !Compute of Ledd
    Ledd = 4._x_precision * pi * G * mp * c / thomson

    !Opening file

    open(unit=11, file="./input_parameter.dat", &
         action="read", status="old", iostat=ios)
    if (ios .ne. 0) stop "OPENING input_parameter.dat ERROR"

    !Reading parameters
    read(11,fmt=*) bla, initial_param%M
    read(11,fmt=*) bla, rmax
    read(11,fmt=*) bla, initial_param%Mdot
    read(11,fmt=*) bla, initial_param%alpha
    read(11,fmt=*) bla, initial_param%X
    read(11,fmt=*) bla, Y

    !Compute M and M_dot in cgs
    initial_param%Mdot = initial_param%Mdot * ( Ledd * 12._x_precision / c**2) * initial_param%M * M_sun
    initial_param%M = initial_param%M * M_sun

    !Processing r_s and T_0
    rs  = 2._x_precision * G * initial_param%M / (c*c)
    T_0  = (1._x_precision/sqrt(27.0) * 1._x_precision/48._x_precision * &
           initial_param%Mdot * (c*c) / ( pi * rs * rs * stefan ))**0.25_x_precision

    !Processing mu
    Z = 1._x_precision - initial_param%X - Y
    mu = 1._x_precision / (2._x_precision*initial_param%X + 3._x_precision*Y/4._x_precision + Z/2._x_precision)

    !Compute dx, rmin = 3rs
    initial_param%dx = (sqrt(rmax) - sqrt(3._x_precision)) / (n_cell - 1._x_precision)

    !Compute RTM
    initial_param%RTM = kmp * T_0 / mu

    close(11)

  end subroutine get_parameters

end module mod_read_parameters
