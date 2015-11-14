! Module that computes the variables, but also their dimension transformations

module mod_variables
  use mod_constants
  use mod_read_parameters

  implicit none

  type(adim_state) :: state_0

  private

  public           :: compute_variables, dim_adim, init_variable_0, state_0

contains

  ! Compute all variables for a given state, i.e. a new S and/or T
  subroutine compute_variables(state_out)
    implicit none

    integer                                     :: i
    real(kind = x_precision), dimension(n_cell) :: a1, b1, c1
    real(kind = x_precision), dimension(n_cell) :: Delta
    real(kind = x_precision), dimension(n_cell) :: kappa_ff, tau, epsilo
    real(kind = x_precision)                    :: kappa_e !cgs

    type(state), intent(inout)                  :: state_out

    ! Compute trinomial coefficients for H
    a1 = ((state_out%Omega**2) * (state_0%Omega_0**2) * state_out%S * state_0%S_0) / (2._x_precision * state_out%x)
    b1 = - (cst_rad * (state_out%T**4) * (state_0%T_0**4)) / (3._x_precision)
    c1 = - (params%RTM * state_out%T * state_out%S * state_0%S_0 ) / (2._x_precision * state_out%x)
    Delta = b1**2 - (4._x_precision * a1 * c1)

    ! Loop over all cells to update variables
    do i=1,n_cell
      state_out%H(i)     = - 0.5_x_precision * (b1(i) + sign(sqrt(Delta(i)),b1(i))) / a1(i) / state_0%H_0
      state_out%P_rad(i) = state_out%T(i)**4
      state_out%cs(i)    = state_out%Omega(i) * state_out%H(i)
      state_out%rho(i)   = state_out%S(i) / (state_out%H(i) * state_out%x(i))
      state_out%nu(i)    = params%alpha * state_out%cs(i) * state_out%H(i)
      state_out%P_gaz(i) = state_out%rho(i) * state_out%T(i)
      state_out%beta(i)  = state_out%P_gaz(i) / (state_out%P_gaz(i) + state_out%P_rad(i))
      state_out%Cv(i)    = params%RTM * ((12._x_precision * gammag - 1._x_precision ) * &
                           (1._x_precision - state_out%beta(i)) + state_out%beta(i)) / &
                           ( state_out%beta(i) * (gammag - 1._x_precision))

      ! Compute v while taking care of limit conditions
      if (i == 1) then
         state_out%v(i)  = 0._x_precision
      else if (i == n_cell) then
         state_out%v(i)  = - 1._x_precision / state_out%S(i) / state_out%x(i)
      else
         state_out%v(i)  = - 1._x_precision / state_out%S(i) / state_out%x(i) * ( ((state_out%nu(i) * &
                           state_out%S(i)) - (state_out%nu(i-1) * state_out%S(i-1))) /  &
                           ( (state_out%x(i) - state_out%x(i-1))))
      endif

      state_out%M_dot(i) = - state_out%v(i) * state_out%S(i) * state_out%x(i)

      ! Compute variables needed for Fz
      kappa_ff(i) = 6.13e22_x_precision * state_0%rho_0 * state_out%rho(i) * &
                    (state_0%T_0 * state_out%T(i))**(-7._x_precision/2._x_precision)

      tau(i)      = 0.5_x_precision * sqrt(params%kappa_e * kappa_ff(i)) * (state_out%S(i)/state_out%x(i) * state_0%S_0)

      epsilo(i)   = 6.22e20_x_precision * (state_0%rho_0 * state_out%rho(i))**2 * sqrt((state_0%T_0 * state_out%T(i)))

      ! Compute Fz depending on optical thickness
      if (tau(i) >= 1.0) then
         state_out%Fz(i) = (2._x_precision * c * c * state_out%T(i)**4) / (27._x_precision * &
                           sqrt(3.0) * (kappa_ff(i) + kappa_e) * (state_out%S(i)/state_out%x(i) * state_0%S_0))
      else
         state_out%Fz(i) = epsilo(i) * state_out%H(i) * state_0%temps_0 / state_0%rho_0
      endif
    enddo
  end subroutine compute_variables

  ! Transform dimensions of variables
  subroutine dim_adim(mode,state_in)
    implicit none

    integer,     intent(in)    :: mode     ! Select 0 or 1 to adimension or dimension your state
    type(state), intent(inout) :: state_in

    select case(mode)
    case(0)
      state_in%nu    = state_in%nu    / state_0%nu_0
      state_in%v     = state_in%v     / state_0%v_0
      state_in%cs    = state_in%cs    / state_0%cs_0
      state_in%S     = state_in%S     / state_0%S_0
      state_in%H     = state_in%H     / state_0%H_0
      state_in%M_dot = state_in%M_dot / state_0%M_dot_0
      state_in%rho   = state_in%rho   / state_0%rho_0
      state_in%T     = state_in%T     / state_0%T_0
      state_in%Fz    = state_in%Fz    / state_0%Fz_0
      state_in%Cv    = state_in%Cv    / state_0%Cv_0
      state_in%P_gaz = state_in%P_gaz / state_0%P_gaz_0
      state_in%P_rad = state_in%P_rad / state_0%P_rad_0
    case(1)
      state_in%nu    = state_in%nu    * state_0%nu_0
      state_in%v     = state_in%v     * state_0%v_0
      state_in%cs    = state_in%cs    * state_0%cs_0
      state_in%S     = state_in%S     * state_0%S_0
      state_in%H     = state_in%H     * state_0%H_0
      state_in%M_dot = state_in%M_dot * state_0%M_dot_0
      state_in%rho   = state_in%rho   * state_0%rho_0
      state_in%T     = state_in%T     * state_0%T_0
      state_in%Fz    = state_in%Fz    * state_0%Fz_0
      state_in%Cv    = state_in%Cv    * state_0%Cv_0
      state_in%P_gaz = state_in%P_gaz * state_0%P_gaz_0
      state_in%P_rad = state_in%P_rad * state_0%P_rad_0
    case default
       stop
    end select
  end subroutine dim_adim

  ! Compute the initial parameters needed for dimensionless variables
  subroutine init_variable_0()
    implicit none

    real(kind = x_precision)         :: rs, c2

    c2 = c**2
    rs = 2._x_precision*G*params%M/c2

    state_0%omega_0 = sqrt( G*params%M/(rs)**3 )
    state_0%temps_0 = 2._x_precision / state_0%omega_0
    state_0%nu_0    = 2._x_precision/3._x_precision * state_0%omega_0 * rs**2
    state_0%v_0     = state_0%omega_0 * rs
    state_0%cs_0    = state_0%v_0
    state_0%S_0     = params%Mdot / (3._x_precision * pi * state_0%nu_0 )
    state_0%H_0     = rs
    state_0%M_dot_0 = params%Mdot
    state_0%rho_0   = state_0%S_0 / (2._x_precision*rs)
    state_0%T_0     = (1._x_precision/sqrt(27.0) * 1._x_precision/48._x_precision * &
                      params%Mdot * c2 / ( pi * rs * rs * stefan ))**0.25_x_precision
    state_0%Fz_0    = state_0%S_0 / 2._x_precision / state_0%temps_0
    state_0%Cv_0    = 1._x_precision / state_0%T_0
    state_0%P_gaz_0 = state_0%rho_0 * params%RTM
    state_0%P_rad_0 = 1._x_precision/3._x_precision * cst_rad * (state_0%T_0**4)
  end subroutine init_variable_0

end module mod_variables
