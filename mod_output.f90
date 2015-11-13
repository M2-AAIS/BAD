! Module that provides a way to dump the state

module mod_output
  use mod_constants

  implicit none

  private

  public :: snapshot

contains

  ! Save a snapshot of state (s) at iteration, time
  subroutine snapshot (s, iteration, time, unit)

    implicit none
    
    type (state), intent(in)               :: s
    real (kind=x_precision), intent(in)    :: time
    integer, intent(in)                    :: iteration, unit
    
    integer                                :: i 

    write(unit, fmt=*) '# ', iteration
    write(unit, fmt=*) '# ', time
    write(unit, fmt='(2a, 15(a16))') '# ', 'x', 'Omega', 'nu', 'v', 'T', 'P_rad',&
         'P_gaz', 'beta', 'cs', 'H', 'rho', 'S', 'Fz', 'M_dot', 'Cv'

    do i = 1, n_cell
       write(unit, fmt='(15(e16.8e2))') s%x(i), s%Omega(i), s%nu(i), s%v(i),&
            s%T(i), s%P_rad(i), s%P_gaz(i), s%beta(i), s%cs(i),&
            s%H(i), s%rho(i), S%S(i), S%Fz(i), S%M_dot(i), S%Cv(i)
    end do
  end subroutine snapshot
  
end module mod_output

