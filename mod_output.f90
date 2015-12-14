! Module that provides a way to dump the state

module mod_output
  use mod_constants
  use mod_read_parameters
  use mod_variables

  implicit none

  private

  public :: snapshot

contains

  ! Save a snapshot of state (s) at iteration, time
  subroutine snapshot (s, iteration, time, unit)

    implicit none

    real(x_precision), intent(in) :: time
    integer,           intent(in) :: iteration
    integer,           intent(in) :: unit
    type(state),       intent(inout) :: s

    integer     :: i
    type(state) :: state_dim

    write(unit, fmt=*) '# ', iteration
    write(unit, fmt=*) '# ', time
    write(unit, fmt='(21(a16))') 'r', &
                                 !'Omega', &
                                 !'nu', &
                                 !'v', &
                                 'T', &
                                 'M_dot', &
                                 !'P_rad', &
                                 !'P_gas', &
                                 !'beta', &
                                 !'cs', &
                                 'H', &
                                 !'rho', &
                                 !'S', &
                                 'Sigma', &
                                 'tau', &
                                 !'Fz', &
                                 'Cv', &
                                 'Q_+', &
                                 'Q_-', &
                                 'f', &
                                 'Qadv'

    call compute_variables(s)

    state_dim = s
    call dim_adim(1, state_dim)
    do i = 1, n_cell
       write(unit, fmt='(21(e16.8e2))') r_state%r(i), &
                                        !r_state%Omega_r(i), &
                                        !state_dim%nu(i), &
                                        !state_dim%v(i), &
                                        state_dim%T(i), &
                                        state_dim%Mdot(i), &
                                        !state_dim%Prad(i), &
                                        !state_dim%Pgaz(i), &
                                        !state_dim%beta(i), &
                                        !state_dim%cs(i), &
                                        state_dim%H(i), &
                                        !state_dim%rho(i), &
                                        !state_dim%S(i), &
                                        state_dim%S(i) / x_state%x(i), &
                                        s%tau(i), &
                                        !state_dim%Fz(i), &
                                        state_dim%Cv(i), &
                                        s%Qplus(i), &
                                        s%Qminus(i), &
                                        s%Qplus(i) - s%Qminus(i), &
                                        s%Qadv(i)
    end do

  end subroutine snapshot

end module mod_output

