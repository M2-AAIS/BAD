! Module that computes the variables
module mod_variables
  use mod_constants
  use mod_read_parameters
  implicit none
  type(adim_state)               :: state_0

  private

  public :: compute_variables, dim_adim, init_variable_0, state_0
  
contains
  !--------------------------------------------------------------
  subroutine compute_variables(x, Omega, T_in, S_in, state_out)
    use mod_constants
    use mod_read_parameters
    implicit none
    integer                                                   :: i
    type(parameters)                                          :: para
    real (kind = x_precision),              dimension(n_cell) :: a1, b1, c1
    real (kind = x_precision),              dimension(n_cell) :: Delta
    real (kind = x_precision),              dimension(n_cell) :: kappa_ff, tau, epsil
    real (kind = x_precision)                                 :: kappa_e !cgs
    
    real (kind = x_precision), intent(in),  dimension(n_cell) :: x, Omega
    real (kind = x_precision), intent(in),  dimension(n_cell) :: T_in, S_in
    type(state),               intent(out)                    :: state_out

    call get_parameters(para)
    kappa_e = 0.2_x_precision * (1._x_precision + para%X)
    !-------------Compute trinome coeff for H----------------
    a1 = ((Omega**2) * (state_0%Omega_0**2) * S_in * state_0%S_0) / (2._x_precision * x)
    b1 = - (cst_rad * (T_in**4) * (state_0%T_0**4)) / (3._x_precision)
    c1 = - (kmp * T_in * S_in * state_0%S_0 * state_0%T_0) / (2._x_precision * mu * x)
    Delta = b1**2 - (4._x_precision * a1 * c1)
    !--------------------------------------------------------
    do i=1,n_cell
       state_out%H(i)     = - 0.5_x_precision * (b1(i) + sign(sqrt(Delta(i)),b1(i))) / a1(i) 
       state_out%P_rad(i) = T_in(i)**4
       state_out%cs(i)    = Omega(i) * state_out%H(i)
       state_out%rho(i)   = S_in(i) / (state_out%H(i) * x(i))
       state_out%nu(i)    = para%alpha * state_out%cs(i) * state_out%H(i)
       state_out%P_gaz(i) = state_out%rho(i) * (T_in(i)**4)
       state_out%beta(i)  = state_out%P_gaz(i) / (state_out%P_gaz(i) + state_out%P_rad(i))
       state_out%Cv(i)    = (kmp/mu) * ((12._x_precision * gammag -1._x_precision ) &
            * (1._x_precision - state_out%beta(i)) + state_out%beta(i)) &
            / ( state_out%beta(i) * (gammag - 1._x_precision))
       !------------limit condition to compute v-------------
       if (i .eq. 1) then
          state_out%v(i)  = 0._x_precision
       else if (i .eq. n_cell) then
          state_out%v(i)  = - 1._x_precision / S_in(i) / x(i) 
       else
          state_out%v(i)  = - 1._x_precision / S_in(i) / x(i) * ( ((state_out%nu(i+1) &
               * S_in(i+1)) - (state_out%nu(i-1) * S_in(i-1))) / (2._x_precision * (x(i+1) - x(i-1))))
       endif
       state_out%M_dot(i) = - state_out%v(i) * S_in(i) * x(i)
      !--------------Compute varaibles for Fz---------------
       kappa_ff(i) = 6.13e22_x_precision * state_0%rho_0 * state_out%rho(i) * &
            (state_0%T_0 * T_in(i))**(-7._x_precision/2._x_precision) 
       tau(i)      = 0.5_x_precision * sqrt(0.2_x_precision * (1._x_precision * para%X) * 6.13e22 * state_0%rho_0 &
            * state_out%rho(i) * (state_0%T_0 * T_in(i))**(-7._x_precision/2._x_precision)) * state_0%S_0 &
            * S_in(i) / x(i)
       epsil(i)    = 6.22e20_x_precision * (state_0%rho_0 * state_out%rho(i))**2  &
            * sqrt((state_0%T_0 * T_in(i))) 

       if (tau(i) .ge. 1.0) then
          state_out%Fz(i) = (4._x_precision * c * c * T_in(i)**4) / (27._x_precision * sqrt(3.0) * &
               (kappa_ff(i) + kappa_e) * (S_in(i)/x(i) * state_0%S_0)) 
       else
          state_out%Fz(i) = epsil(i) * state_0%H_0 * state_out%H(i)
       endif
    enddo
  end subroutine compute_variables
  !--------------------------------------------------------------
  subroutine dim_adim(mode,state_in)
    use mod_constants
    use mod_read_parameters
    implicit none
    !select 0 or 1 to adimension or dimension your state 
    integer         , intent(in)                 :: mode
    type(state)     , intent(inout)              :: state_in
	
    select case(mode)
    case(0)
       state_in%H     = state_in%H / state_0%H_0
       state_in%v     = state_in%v / state_0%v_0
       state_in%cs    = state_in%cs / state_0%cs_0
       state_in%S     = state_in%S / state_0%S_0
       state_in%T     = state_in%T / state_0%T_0
       state_in%rho   = state_in%rho / state_0%rho_0
       state_in%nu    = state_in%nu / state_0%nu_0
       state_in%M_dot = state_in%M_dot / state_0%M_dot_0
       state_in%P_rad = state_in%P_rad / state_0%P_rad_0
       state_in%P_gaz = state_in%P_gaz / state_0%P_gaz_0
       state_in%Fz    = state_in%Fz / state_0%Fz_0
    case(1)
       state_in%H     = state_in%H * state_0%H_0
       state_in%v     = state_in%v * state_0%v_0
       state_in%cs    = state_in%cs * state_0%cs_0
       state_in%S     = state_in%S * state_0%S_0
       state_in%T     = state_in%T * state_0%T_0
       state_in%rho   = state_in%rho * state_0%rho_0
       state_in%nu    = state_in%nu * state_0%nu_0
       state_in%M_dot = state_in%M_dot * state_0%M_dot_0
       state_in%P_rad = state_in%P_rad * state_0%P_rad_0
       state_in%P_gaz = state_in%P_gaz * state_0%P_gaz_0
       state_in%Fz    = state_in%Fz * state_0%Fz_0
    case default
       stop
    end select
  end subroutine dim_adim
 !--------------------------------------------------------------
  subroutine init_variable_0()
    use mod_constants
    use mod_read_parameters
    !Compute the initial adimention parameters
    real(kind = x_precision)         :: omega_max, rs, c2
    type(parameters)                 :: para

    c2 = c**2
    
    call get_parameters(para)
    rs            = 2._x_precision*G*para%M/c2
    omega_max     = sqrt( G*para%M/(rs)**3 )
    
    state_0%temps_0 = 2._x_precision / omega_max
    state_0%x_0     = sqrt(rs)
    state_0%H_0     = rs
    state_0%nu_0    = 2._x_precision/3._x_precision * omega_max * rs**2
    state_0%omega_0 = omega_max
    state_0%v_0     = omega_max * rs
    state_0%cs_0    = state_0%v_0
    state_0%S_0     = para%Mdot / (3._x_precision * pi * state_0%nu_0 )
    state_0%T_0     = (1._x_precision/sqrt(27.0) * 1._x_precision/48._x_precision * para%Mdot * c2 &
         / ( pi * rs * rs * stefan ))**0.25_x_precision
    state_0%rho_0   = state_0%S_0 / (2._x_precision*rs)
    
    state_0%M_dot_0 = para%Mdot
    state_0%P_rad_0 = 1._x_precision/3._x_precision * cst_rad * (state_0%T_0**4)
    state_0%P_gaz_0 = state_0%rho_0 * kmp * state_0%T_0 / mu
    state_0%Fz_0    = state_0%S_0 / 2._x_precision / state_0%temps_0
  end subroutine init_variable_0
  
end module mod_variables
