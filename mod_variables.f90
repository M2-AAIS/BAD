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
    use mod_read_parameters
    implicit none
    !select 0 or 1 to adimension or dimension your state 
    real(kind = x_precision)                     :: omega_max
    real(kind = x_precision)                     :: rs
    
    integer         , intent(in)                 :: type
    type(state)     , intent(inout)              :: state_in
    type(parameters)                             :: para
    call get_parameters ( para )

    rs = 2*G*para%M/c/c
	
    select case(type)
    case(0)
       !write here the adimension conversion

       omega_max = sqrt( G*para%M/(3*rs)/(3*rs)/(3*rs) )
       state_in%omega = state_in%omega/omega_max

       state_in%H = state_in%H/rs

       state_in%cs = state_in%H * state_in%omega
       
       state_in%v = state_in%v / (omega_max * rs)

       state_in%sigma = state_in%sigma * 2 * pi * rs * rs * omega_max / para%Mdot

       state_in%x = sqrt( state_in%x / rs )

       state_in%s = state_in%sigma * state_in%x

       state_in%m_dot = state_in%v * state_in%s * state_in%x

       state_in%T = state_in%T / ( 1.0/sqrt(27.0) * 1.0/12.0 * para%mdot * c * c / ( 4 * pi * rs * rs * stefan ) )**0.25_x_precision

       state_in%rho = state_in%rho * rs * 2 * pi * rs * rs * omega_max / para%Mdot

       state_in%nu = 3.0/4.0 * state_in%nu / omega_max / rs /rs

       
 



       state_in%M_dot = state_in%M_dot       
       print*,para%m
     
    case(1)
       !write here the dimension conversion
    case default
       write(*,*)'Enter 0 to adim and 1 to dim'
    end select
  end subroutine dim_adim
  
end module mod_variables
