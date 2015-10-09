! Module that reads the parameter out of a file

module mod_read_parameters
  
  private

  public :: get_parameters
  
contains
  
  subroutine get_parameters (filename, params)
    use mod_constants
    implicit none
    
    character(len=*), intent(in) :: filename
    type (parameters), intent(out) :: params
    
  end subroutine get_parameters
  
end module mod_read_parameters
