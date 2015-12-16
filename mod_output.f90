! Module that provides a way to dump the state

module mod_output
  use mod_constants
  use mod_read_parameters
  use mod_variables

  implicit none

  private

  logical :: write_headers = .true.

  public :: snapshot

contains

  ! Save a snapshot of state (s) at iteration, time
  subroutine snapshot (s, iteration, time, unit_full, unit_int)

    implicit none

    real(x_precision), intent(in) :: time
    integer,           intent(in) :: iteration
    integer,           intent(in) :: unit_full ! write the full state in this unit 
    integer,           intent(in) :: unit_int  ! write integrated quantities in this unit
    type(state),       intent(inout) :: s

    integer     :: i
    type(state) :: state_dim

    real(x_precision) :: dim_time, Ltot, Mtot, rs

    dim_time = time * state_0%temps_0

    write(unit_full, fmt=*) '# ', iteration
    write(unit_full, fmt=*) '# ', dim_time
    write(unit_full, fmt='(21(a16))') 'r', &
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
       write(unit_full, fmt='(21(e16.8e2))') r_state%r(i), &
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

    !-----------------------
    ! Compute some integrated variables on the disk
    !-----------------------
    rs = 2._x_precision * G * params%M / c**2
    Ltot = 8*pi*sum(x_state%x * state_dim%Fz)
    Mtot = 4*pi*rs**2 *sum(x_state%x**3 * state_dim%rho * state_dim%H)

    if (write_headers) then
       write (unit_int, fmt='(4(a16))') 'iteration', 'time', 'Ltot', 'Mtot'
       write_headers = .false.
    end if
    write (unit_int, fmt='(i16, 3(e16.8e2))') iteration, dim_time, Ltot, Mtot

  end subroutine snapshot

end module mod_output

