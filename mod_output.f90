! Module that provides a way to dump the state

module mod_output
  
  private

  public :: snapshot
  
contains
  
  subroutine snapshot (state_in, filename)
    use mod_constants
    implicit none
    
    type (state), intent(in), dimension(:)  :: state_in
    character(len=*), intent(in)            :: filename
    
  end subroutine snapshot
  
end module mod_output

