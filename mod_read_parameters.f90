! Module that reads the parameter out of a file

module mod_read_parameters
  use mod_constants
  implicit none
  private

  public :: get_parameters, process_dx
  
contains
  
  subroutine get_parameters (initial_param)
    !read parameter in the input file : ./input_parameter.dat

    implicit none
    

    ! function parameters
    type(parameters), intent(out) :: initial_param ! contains the initial parameters of Black hole accretion disk (see mod_constants.f90 for details)
    
    ! internal variable
    real (kind=x_precision)       :: Y,Z           ! chemical composition : x+y+z = 1 
    integer (kind=4)              :: ios           ! i/o variable test 
    character (len=50)            :: bla           ! string reading variable
    real(kind=x_precision)        :: Ledd

    !Compute pf Ledd
    Ledd = 4._x_precision * pi * G * mp * c / thomson

    !Opening file

    open(unit=11, file="./input_parameter.dat", &
         action="read", status="old", iostat=ios)
    if (ios .ne. 0) stop "OPENING input_parameter.dat ERROR"

    
    !reading parameter
    read(11,fmt=*) bla, initial_param%M
    read(11,fmt=*) bla, initial_param%rmax
    read(11,fmt=*) bla, initial_param%Mdot
    read(11,fmt=*) bla, initial_param%alpha
    read(11,fmt=*) bla, initial_param%X
    read(11,fmt=*) bla, Y

    !Compute M and M_dot in cgs
    initial_param%Mdot = initial_param%Mdot * ( Ledd * 12._x_precision / c**2) * initial_param%M * M_sol
    initial_param%M = initial_param%M * M_sol

    !processing mu
    Z = 1._x_precision-initial_param%X-Y
    initial_param%mu = 1._x_precision / (2._x_precision*initial_param%X + 3._x_precision*Y/4._x_precision + Z/2._x_precision)
    

    close(11)
    
  end subroutine get_parameters



  subroutine process_dx(dx)
  !process the space step

  !function parameter
  real (kind=x_precision), intent(out) :: dx !space step
  
  dx = (10._x_precision - sqrt(3._x_precision)) / n_cell

  end subroutine process_dx
  
end module mod_read_parameters
