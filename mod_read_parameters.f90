! Module that reads the parameter out of a file

module mod_read_parameters
  
  private

  public :: get_parameters
  
contains
  
  subroutine get_parameters (initial_param)
    use mod_constants
    !read parameter in the input file : ./input_parameter.dat

    implicit none
    

    !function parameters
    type(parameters),intent(out)::initial_param !contains the initial parameters of Black hole accretion disk (see mod_constants.f90 for details)
    
    !internal variable
    real(kind=x_precision)::Y,Z !chemical composition : x+y+z = 1 
    integer(kind=4)::ios !i/o variable test 
    character(len=50)::bla !string reading variable

    !Opening file

    open(unit=11,file="./input_parameter.dat",&
         action="read",status="old",iostat=ios)
    if(ios.ne.0) stop "OPENING input_parameter.dat ERROR"

    
    !reading parameter
    read(11,fmt=*)bla,initial_param%M
    read(11,fmt=*)bla,initial_param%rmax
    read(11,fmt=*)bla,initial_param%Mdot
    read(11,fmt=*)bla,initial_param%alpha
    read(11,fmt=*)bla,initial_param%X
    read(11,fmt=*)bla,Y

    !processing mu
    Z=1.D0-initial_param%X-Y
    initial_param%mu=1.D0/(2.D0*initial_param%X + 3.D0*Y/4.D0 + Z/2.D0)
    
    

    

    close(11)
    
  end subroutine get_parameters
  
end module mod_read_parameters
