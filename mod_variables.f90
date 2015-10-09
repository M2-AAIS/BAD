! Module that computes the variables

module mod_variables
  
  private

  public :: compute_variables
  
contains
  
  subroutine compute_variables (T, Sigma, state_out)
    use mod_constants
    implicit none
    real (kind = x_precision), dimension(:) :: T, Sigma
    type (state), intent(in), dimension(:)  :: state_out
    
  end subroutine compute_variables
  
end module mod_variables
