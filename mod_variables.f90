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

    type(state), intent(inout)                  :: state_out

    ! Compute trinomial coefficients for H
    !a1 = (state_out%Omega * state_0%Omega_0)**2 * (state_out%S * state_0%S_0)
    !b1 = - 2._x_precision * cst_rad * (state_out%T * state_0%T_0)**4 * state_out%x / 3._x_precision
    !c1 = - (params%RTM * state_out%T * state_out%S * state_0%S_0)
    a1 = (state_out%Omega)**2 * state_0%Omega_0 * state_0%M_dot_0 * (state_out%S * state_0%S_0)
    b1 = - 4._x_precision * pi * cst_rad * (state_out%T * state_0%T_0)**4 * state_0%H_0 * state_out%x / 3._x_precision
    c1 = - 2._x_precision * pi * (params%RTM * state_out%T * state_out%S * state_0%S_0)
    Delta = b1**2 - 4._x_precision * a1 * c1

    ! Start computing variables
    !state_out%H     = - 0.5_x_precision * (b1 + sign(sqrt(Delta),b1)) / a1 / state_0%H_0
    state_out%H     = - 0.5_x_precision * (b1 + sign(sqrt(Delta),b1)) / a1
    state_out%P_rad = state_out%T**4
    state_out%cs    = state_out%Omega * state_out%H
    state_out%rho   = state_out%S / (state_out%H * state_out%x)
    state_out%nu    = params%alpha * state_out%cs * state_out%H
    state_out%P_gaz = state_out%rho * state_out%T
    state_out%beta  = state_out%P_gaz / (state_out%P_gaz + state_out%P_rad)
    state_out%Cv    = params%RTM * ((12._x_precision * gammag - 1._x_precision) * &
                      (1._x_precision - state_out%beta) + state_out%beta) / &
                      (state_out%beta * (gammag - 1._x_precision))

    ! Compute v while taking care of limit conditions
    state_out%v(1:n_cell-1) = - 1._x_precision / (state_out%S(1:n_cell-1) * state_out%x(1:n_cell-1)) * &
                              (state_out%nu(2:n_cell) * state_out%S(2:n_cell) - &
                               state_out%nu(1:n_cell-1) * state_out%S(1:n_cell-1)) / params%dx
    state_out%v(n_cell)     = - 1._x_precision / (state_out%S(n_cell) * state_out%x(n_cell))

    state_out%M_dot = - state_out%v * state_out%S * state_out%x

    ! Compute variables needed for Fz
    kappa_ff = 6.13e22_x_precision * state_0%rho_0 * state_out%rho * &
               (state_0%T_0 * state_out%T)**(-7._x_precision/2._x_precision)

    tau      = 0.5_x_precision * sqrt(params%kappa_e * kappa_ff) * (state_0%S_0 * state_out%S / state_out%x)

    epsilo   = 6.22e20_x_precision * (state_0%rho_0 * state_out%rho)**2 * sqrt(state_0%T_0 * state_out%T)

    ! Loop over all cells to update Fz
    do i=1,n_cell
      ! Compute Fz depending on optical thickness
      if (tau(i) >= 1.0) then
         state_out%Fz(i) = (2._x_precision * c * c * state_out%T(i)**4) / (27._x_precision * &
                           sqrt(3.0) * (kappa_ff(i) + params%kappa_e) * (state_out%S(i)/state_out%x(i) * state_0%S_0))
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

    real(kind = x_precision) :: rs, c2

    c2 = c**2
    rs = 2._x_precision*G*params%M/c2

    state_0%omega_0 = sqrt(G * params%M / rs**3)
    state_0%temps_0 = 2._x_precision / state_0%omega_0
    state_0%nu_0    = 2._x_precision * state_0%omega_0 * rs**2 / 3._x_precision
    state_0%v_0     = state_0%omega_0 * rs
    state_0%cs_0    = state_0%v_0
    state_0%S_0     = params%Mdot / (3._x_precision * pi * state_0%nu_0)
    state_0%H_0     = rs
    state_0%M_dot_0 = params%Mdot
    state_0%rho_0   = state_0%S_0 / (2._x_precision * rs)
    state_0%T_0     = (params%Mdot * c2 / (sqrt(27.0) * 48._x_precision * pi * rs**2 * stefan))**0.25_x_precision
    state_0%Fz_0    = state_0%S_0 / (2._x_precision * state_0%temps_0)
    state_0%Cv_0    = 1._x_precision / state_0%T_0
    state_0%P_gaz_0 = state_0%rho_0 * params%RTM
    state_0%P_rad_0 = cst_rad * state_0%T_0**4 / 3._x_precision
  end subroutine init_variable_0

end module mod_variables
