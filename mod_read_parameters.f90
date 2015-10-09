! Module that reads the parameter out of a file

module mod_read_parameters
  
  private

  public :: get_parameters
  
contains
  
  subroutine get_parameters (M,r_max,M_0_dot,alpha,X,mu)
    use mod_constants
    !read parameter in the input file : ./input_parameter.dat
    implicit none
    
       !function parameters
    real(kind=x_precision),intent(out)::M !Black hole mass
    real(kind=x_precision),intent(out)::r_max !maximal radius of black hole
    real(kind=x_precision),intent(out)::M_0_dot !accretion rate at r_max
    real(kind=x_precision),intent(out)::alpha !viscosity parameter
    real(kind=x_precision),intent(out)::X !chemical composition
    real(kind=x_precision),intent(out)::mu !atomic mass

    !internal variable
    real(kind=x_precision)::Y,Z !chemical composition
    integer(kind=4)::ios !i/o variable
    character(len=50)::bla !string reading variable

    !Opening file

    open(unit=11,file="./input_parameter.dat",&
         action="read",status="old",iostat=ios)
    if(ios.ne.0) stop "OPENING input_parameter.dat ERROR"

    
    !reading parameter
    read(11,fmt=*)bla,M
    read(11,fmt=*)bla,r_max
    read(11,fmt=*)bla,M_0_dot
    read(11,fmt=*)bla,alpha
    read(11,fmt=*)bla,X
    read(11,fmt=*)bla,Y

    !processing mu
    Z=1.D0-X-Y
    mu=1.D0/(2.D0*X + 3.D0*Y/4.D0 + Z/2.D0)
    
    

    

    close(11)
    
  end subroutine get_parameters
  
end module mod_read_parameters
