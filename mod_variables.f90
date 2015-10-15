! Module that computes the variables
module mod_variables
  
  private

  public :: compute_variables, dim_adim
  
contains
  
  subroutine compute_variables (T, Sigma, state_out)
    use mod_constants
    implicit none
    real (kind = x_precision), dimension(:) :: T, Sigma
    type (state), intent(in), dimension(:)  :: state_out
    
  end subroutine compute_variables

  subroutine dim_adim(type,state_in)
    use mod_constants
    implicit none
    !select 0 or 1 to dimension or adimension your state 
    integer         , intent(in)                 :: type
    type(state)     , intent(inout)              :: state_in

    select case(type)
    case(0)
       !write here the adimension conversion
    case(1)
       !write here the dimension conversion
    case default
       write(*,*)'Enter 0 to adim and 1 to dim'
    end select
  end subroutine dim_adim
  
end module mod_variables
