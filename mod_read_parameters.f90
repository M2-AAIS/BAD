! Module that expose the initial (physical) parameters after parsing them from a file

module mod_read_parameters
  use mod_constants
  implicit none
  type(parameters) :: params ! Contains the initial parameters of the black hole accretion disk (see mod_constants.f90 for details)
  
  private

  public           :: get_parameters, params

contains

  subroutine get_parameters ()
    ! Read parameter in the input file: ./input_parameter.dat

    implicit none

    ! Internal variables
    real (kind=x_precision) :: Y,Z  ! Chemical composition : x+y+z = 1
    real (kind=x_precision) :: mu   ! Mean molecular mass
    real (kind=x_precision) :: rmax ! Maximum considered radius (in rs)
    real (kind=x_precision) :: Ledd ! Eddington luminosity
    real (kind=x_precision) :: rs   ! Schwarzschild radius
    real (kind=x_precision) :: T_0  ! Temperature order of magnitude
    integer (kind=4)        :: ios  ! i/o variable test
    character (len=50)      :: line ! string reading variable

    !Compute of Ledd
    Ledd = 4._x_precision * pi * G * mp * c / thomson

    !Opening file

    open(unit=11, file="./input_parameter.dat", &
         action="read", status="old", iostat=ios)
    if (ios .ne. 0) stop "OPENING input_parameter.dat ERROR"

    !Reading parameters
    read(11,fmt=*) line, params%M
    read(11,fmt=*) line, rmax
    read(11,fmt=*) line, params%Mdot
    read(11,fmt=*) line, params%alpha
    read(11,fmt=*) line, params%X
    read(11,fmt=*) line, Y

    !Compute M and M_dot in cgs
    params%Mdot = params%Mdot * ( Ledd * 12._x_precision / c**2) * params%M * M_sun
    params%M = params%M * M_sun

    !Processing r_s and T_0
    rs  = 2._x_precision * G * params%M / (c*c)
    T_0  = (1._x_precision/sqrt(27.0) * 1._x_precision/48._x_precision * &
           params%Mdot * (c*c) / ( pi * rs * rs * stefan ))**0.25_x_precision

    !Processing mu
    Z = 1._x_precision - params%X - Y
    mu = 1._x_precision / (2._x_precision*params%X + 3._x_precision*Y/4._x_precision + Z/2._x_precision)

    !Compute dx, rmin = 3rs
    params%dx = (sqrt(rmax) - sqrt(3._x_precision)) / (n_cell - 1._x_precision)

    !Compute RTM
    params%RTM = kmp * T_0 / mu

    close(11)

  end subroutine get_parameters

end module mod_read_parameters
