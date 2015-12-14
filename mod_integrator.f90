! Module that do integration of S and T given the current state and timestep.

module mod_integrator
  use mod_constants
  use mod_timestep
  use mod_read_parameters
  use mod_variables

  implicit none
  private

  public :: do_timestep_T, do_timestep_S_imp, do_timestep_S_exp

contains

  ! Process the right term of dT*/dt*. Simplified case without convection.
  function dT_dt_imp(s)
    implicit none

    type(state), intent(in) :: s

    real(x_precision), dimension(n_cell) :: dT_dt_imp

    ! Second member of the dT/dt equation. Simplified case.
    dT_dt_imp = (s%Qplus - s%Qminus) / s%Cv

  end function dT_dt_imp

  ! Process the right term of dT*/dt*. Full expression with convection.
  function dT_dt_exp(s)
    implicit none

    type(state), intent(in) :: s

    real(x_precision), dimension(n_cell) :: dT_dt_exp

    ! Second member of the dT/dt equation.
    dT_dt_exp = dT_dt_imp(s) + s%Qadv / s%CV

  end function dT_dt_exp

  ! Process the temporal evolution of S with an implicit algorithm.
  ! For use on the stable part of the simulation.
  subroutine do_timestep_S_imp(s, dt)
    implicit none

    real(x_precision), intent(in)    :: dt
    type(state),       intent(inout) :: s

    real(x_precision), dimension(n_cell)   :: diag, x2
    real(x_precision), dimension(n_cell-1) :: diag_low, diag_up
    real(x_precision)                      :: overdx, overdx2, dtoverdx2
    integer                                :: info


    overdx    = 1._x_precision / params%dx
    overdx2   = overdx**2
    dtoverdx2 = dt*overdx2

    x2 = x_state%x**2

    ! Create the diagonals NOT VERIFIED YET
    diag     =  1._x_precision + 2._x_precision * s%nu / x2 * dtoverdx2
    diag_up  = -s%nu(2:n_cell  ) / x2(1:n_cell-1) * dtoverdx2
    diag_low = -s%nu(1:n_cell-1) / x2(2:  n_cell) * dtoverdx2

    diag(n_cell)       = 1._x_precision + s%nu(n_cell) / x2(n_cell) * dtoverdx2
    diag_low(n_cell-1) = - s%nu(n_cell-1) / x2(n_cell) * dtoverdx2

    s%S(n_cell) = s%S(n_cell) + dt * s%Mdot(n_cell) * overdx / x2(n_cell)

    ! Solving for S, this does modify s%S directly
    call dgtsv(n_cell, 1, diag_low, diag, diag_up, s%S, n_cell, info)

    if (info /= 0) then
      stop "Error in DGTSV call in subroutine do_timestep_S_imp"
    end if

  end subroutine do_timestep_S_imp

  ! Process the temporal evolution of S with an explicit algorithm.
  ! For use on the instable part of the simulation.
  subroutine do_timestep_S_exp(s, dt)
    implicit none

    real(x_precision), intent(in)    :: dt
    type(state),       intent(inout) :: s

    s%S = s%S + dt * dS_dt(s)

  end subroutine do_timestep_S_exp

  ! Process the temporal evolution of T with an explicit algorithm.
  subroutine do_timestep_T(s, dt)
    implicit none

    real(x_precision), intent(in)    :: dt
    type(state),       intent(inout) :: s


    real(x_precision), dimension(n_cell) :: dtemp, rhs, newT
    real(x_precision), dimension(n_cell) :: f0, fT
    type(state)                          :: s_deriv
    integer :: i

    dtemp     = s%T * 1.e-3_x_precision
    s_deriv   = s
    s_deriv%T = s%T + dtemp

    call compute_variables(s_deriv)
    call compute_variables(s)
    ! Let dT/dt = f0 and d²T/dt² = f1 so T(t+1) = T(t) + dt * f0 * (1 + dt * f1)

    f0 = dT_dt_imp(s)
    fT = (dT_dt_imp(s_deriv) - f0) / dtemp

    rhs = f0 / fT * (exp(fT*dt) - 1._x_precision)
    ! rhs = dt*f0 * (1._x_precision + 0.5_x_precision*dt*fT)

    newT = s%T + rhs

    ! T should not go negative, if it does we’re in trouble
    if (minval(newT) < 0) then
       print*, 'Convergence problem!'
       do i = 1, n_cell
           print*, i, newT(i)
       end do
       stop
    end if

    ! When the minimum value of T is ≤0, reduce the timestep
    ! do while (minval(newT) < 0)
    !    print*, 'Convergence problem. Reducing dt', dt
    !    dt = dt / 2
    !    rhs = f0 / fT * (exp(fT*dt) - 1)
    !    newT = s%T + rhs
    ! end do

    s%T = newT

  end subroutine do_timestep_T

end module mod_integrator
