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

  subroutine dim_adim(mode,state_0,state_in)
    use mod_constants
    use mod_read_parameters
    implicit none
    !select 0 or 1 to adimension or dimension your state 
    integer         , intent(in)                 :: mode
    type(adim_state), intent(inout)              :: state_0
    type(state)     , intent(inout)              :: state_in
	
    select case(mode)
    case(0)
       state_in%H = state_in%H / state_0%H
       state_in%v = state_in%v / state_0%v
       state_in%cs = state_in%cs / state_0%cs
       state_in%S = state_in%S / state_0%S
       state_in%T = state_in%T / state_0%T
       state_in%rho = state_in%rho / state_0%rho
       state_in%nu = state_in%nu / state_0%nu
       state_in%M_dot = state_in%M_dot / state_0%M_dot
       state_in%P_rad = state_in%P_rad / state_0%P_rad
       state_in%P_gaz = state_in%P_gaz / state_0%P_gaz
    case(1)
       state_in%H = state_in%H * state_0%H
       state_in%v = state_in%v * state_0%v
       state_in%cs = state_in%cs * state_0%cs
       state_in%S = state_in%S * state_0%S
       state_in%T = state_in%T * state_0%T
       state_in%rho = state_in%rho * state_0%rho
       state_in%nu = state_in%nu * state_0%nu
       state_in%M_dot = state_in%M_dot * state_0%M_dot
       state_in%P_rad = state_in%P_rad * state_0%P_rad
       state_in%P_gaz = state_in%P_gaz * state_0%P_gaz
    case default
       write(*,*)'Enter 0 to adim and 1 to dim'
    end select
  end subroutine dim_adim

  subroutine init_variable_0(state_0)
    !Compute the initial adimention parameters
    use mod_constants
    use mod_read_parameters
    real(kind = x_precision)         :: omega_max
    real(kind = x_precision)         :: rs
    type(parameters)                 :: para
    type(adim_state), intent(inout)  :: state_0

    call get_parameters(para)
    rs = 2*G*para%M/c/c
    omega_max = sqrt( G*para%M/(3*rs)/(3*rs)/(3*rs) )

    state_0%x = sqrt(rs)
    state_0%H = rs
    state_0%omega = omega_max
    state_0%v = omega_max * rs
    state_0%cs = state_0%v
    state_0%S = para%Mdot / (2 * pi * rs * rs * omega_max)
    state_0%T = (1.0/sqrt(27.0) * 1.0/12.0 * para%Mdot * c * c &
         / ( 4 * pi * rs * rs * stefan ))**0.25_x_precision
    state_0%rho = state_0%S / rs
    state_0%nu = 4.0/3.0 * omega_max * rs * rs
    state_0%M_dot = para%Mdot
    state_0%P_rad = 1.0/3.0 * a * (state_0%T **4)
    state_0%P_gaz = state_0%rho * kb_mp * state_0%T / mu
  
  end subroutine init_variable_0
  
end module mod_variables
