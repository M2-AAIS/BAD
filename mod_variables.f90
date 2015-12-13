! Module that computes the variables, but also their dimension transformations

module mod_variables
  use mod_constants
  use mod_read_parameters
  implicit none

  private

  public :: compute_variables, dim_adim, dS_dt

contains

  ! Compute all the variables for a given state, i.e. a new S and/or T
  subroutine compute_variables(state_out)
    implicit none

    type(state), intent(inout) :: state_out

    real(x_precision), dimension(n_cell)   :: coeff_a, coeff_b, coeff_c
    real(x_precision), dimension(n_cell)   :: Delta
    real(x_precision), dimension(n_cell)   :: kappa_ff, epsilo
    real(x_precision), dimension(n_cell)   :: Sx
    real(x_precision), dimension(n_cell)   :: dS_over_dt, dSigma_over_dx, dT_over_dx
    real(x_precision), dimension(n_cell)   :: Sigma
    real(x_precision), dimension(0:n_cell) :: nuS ! one more dimension to store nuS(0)

    !-------------------------------------------------------------
    ! Compute trinomial coefficients for H
    coeff_a = (x_state%Omega * state_0%Omega_0)**2 * (state_out%S * state_0%S_0)
    coeff_b = - 2._x_precision * cst_rad * (state_out%T * state_0%T_0)**4 * x_state%x / (3._x_precision * state_0%H_0)
    coeff_c = - params%RTM * state_out%T * state_out%S * state_0%S_0 / state_0%H_0**2
    Delta = coeff_b**2 - 4._x_precision * coeff_a * coeff_c

    !-------------------------------------------------------------
    ! Compute variable depending on S, H
    state_out%H    = - 0.5_x_precision * (coeff_b + sign(sqrt(Delta), coeff_b)) / coeff_a
    state_out%cs   = x_state%Omega * state_out%H
    state_out%nu   = params%alpha * state_out%cs * state_out%H
    state_out%rho  = state_out%S / (state_out%H * x_state%x)

    !-------------------------------------------------------------
    ! Compute v while taking care of limit conditions. This is not an advection equation, so use mean derivative.
    nuS(0)        = 0
    nuS(1:n_cell) = state_out%nu * state_out%S
    Sx            = state_out%S * x_state%x

    state_out%v(1:n_cell-1) = - 1._x_precision / Sx(1:n_cell-1) * (nuS(2:n_cell) - nuS(0:n_cell-2)) / (2 * params%dx)
    state_out%v(n_cell)     = - 1._x_precision / Sx(n_cell) * params%Mdot_kick_factor

    state_out%Mdot = - state_out%v * Sx

    ! Compute variables needed for Fz
    kappa_ff      = 6.13e22_x_precision * state_0%rho_0 * state_out%rho * &
                    (state_0%T_0 * state_out%T)**(-3.5_x_precision)


    state_out%tau = 0.5_x_precision * sqrt(params%kappa_e * kappa_ff) * (state_0%S_0 * state_out%S / x_state%x)

    epsilo        = 6.22e20_x_precision * (state_0%rho_0 * state_out%rho)**2 * sqrt(state_0%T_0 * state_out%T)

    ! Compute Fz depending on the optical thickness
    !FIXME: use code bellow when optically thick
    ! state_out%Fz = (2._x_precision * c**2 * x_state%x * state_out%T**4) / (27._x_precision * &
    !      sqrt(3._x_precision) * (kappa_ff + params%kappa_e) * state_out%S * state_0%S_0)

    !FIXME: We need to figure why it’s 0.006 and not 1
    where (state_out%tau >= 1)
       state_out%Fz = (2._x_precision * c**2 * x_state%x * state_out%T**4) / (27._x_precision * &
            sqrt(3._x_precision) * (kappa_ff + params%kappa_e) * state_out%S * state_0%S_0)
    elsewhere
       state_out%Fz = epsilo * state_out%H * state_0%temps_0 / state_0%rho_0
    end where

    ! Compute variables related to pressure
    state_out%Pgaz = state_out%rho * state_out%T
    state_out%Prad = state_out%T**4
    state_out%beta = 1._x_precision / (1._x_precision + state_0%beta_0 * state_out%Prad / state_out%Pgaz)
    state_out%Cv   = params%RTM * ((12._x_precision * (gammag - 1._x_precision) * &
                     (1._x_precision - state_out%beta)) + state_out%beta) / &
                     (state_out%beta * (gammag - 1._x_precision))

    ! Compute heating/cooling terms
    state_out%Qplus  = 3._x_precision * state_0%v_0**2 * state_out%nu * x_state%Omega**2
    state_out%Qminus = state_out%Fz * x_state%x / state_out%S

    !-----------------------
    ! Compute advecting term
    !-----------------------

    dS_over_dt = dS_dt(state_out)

    ! Compute dT/dx as the right spatial derivative
    dT_over_dx(1:n_cell-1) = (state_out%T(2:n_cell) - state_out%T(1:n_cell-1)) / params%dx
    dT_over_dx(n_cell)     = 0

    ! Compute d(S/x)/dx as the right spatial derivative
    Sigma = state_out%S / x_state%x

    dSigma_over_dx(1:n_cell-1) = (Sigma(2:n_cell) - Sigma(1:n_cell-1)) / params%dx
    dSigma_over_dx(n_cell)     = - dS_over_dt(n_cell) / state_out%v(n_cell)

    !Compute advecting term
    state_out%Qadv=params%RTM * ( 4._x_precision - 3._x_precision * state_out%beta ) / state_out%beta * &
                 state_out%T / state_out%S * (dS_over_dt + state_out%v * dSigma_over_dx) - &
                 state_out%Cv * state_out%v / x_state%x * dT_over_dx 

    !-----------------------

  end subroutine compute_variables

  ! Transform dimensions of variables
  subroutine dim_adim(mode, state_in)
    implicit none

    integer,     intent(in)    :: mode     ! Select 0 or 1 to adimension or dimension your state
    type(state), intent(inout) :: state_in

    select case(mode)
    case(0)
      state_in%nu   = state_in%nu   / state_0%nu_0
      state_in%v    = state_in%v    / state_0%v_0
      state_in%cs   = state_in%cs   / state_0%cs_0
      state_in%S    = state_in%S    / state_0%S_0
      state_in%H    = state_in%H    / state_0%H_0
      state_in%Mdot = state_in%Mdot / state_0%Mdot_0
      state_in%rho  = state_in%rho  / state_0%rho_0
      state_in%T    = state_in%T    / state_0%T_0
      state_in%Fz   = state_in%Fz   / state_0%Fz_0
      state_in%Cv   = state_in%Cv   / state_0%Cv_0
      state_in%Pgaz = state_in%Pgaz / state_0%Pgaz_0
      state_in%Prad = state_in%Prad / state_0%Prad_0
    case(1)
      state_in%nu   = state_in%nu   * state_0%nu_0
      state_in%v    = state_in%v    * state_0%v_0
      state_in%cs   = state_in%cs   * state_0%cs_0
      state_in%S    = state_in%S    * state_0%S_0
      state_in%H    = state_in%H    * state_0%H_0
      state_in%Mdot = state_in%Mdot * state_0%Mdot_0
      state_in%rho  = state_in%rho  * state_0%rho_0
      state_in%T    = state_in%T    * state_0%T_0
      state_in%Fz   = state_in%Fz   * state_0%Fz_0
      state_in%Cv   = state_in%Cv   * state_0%Cv_0
      state_in%Pgaz = state_in%Pgaz * state_0%Pgaz_0
      state_in%Prad = state_in%Prad * state_0%Prad_0
    case default
       stop
    end select
  end subroutine dim_adim

  
  function dS_dt(s)
    implicit none

    type(state), intent(in) :: s

    real(x_precision), dimension(n_cell)     :: dS_dt

    real(x_precision), dimension(0:n_cell+1) :: nuS ! two more for beginning / end

    ! see BAD-report for explanation of nuS(0) and nuS(n_cell+1)
    nuS(0)        = 0
    nuS(1:n_cell) = s%nu * s%S
    nuS(n_cell+1) = params%Mdot_kick_factor * params%dx + nuS(n_cell)

    dS_dt = 1._x_precision / x_state%x**2 * &
            (nuS(2:n_cell+1) - 2._x_precision * nuS(1:n_cell) + nuS(0:n_cell-1)) / &
            params%dx**2

  end function dS_dt

end module mod_variables
