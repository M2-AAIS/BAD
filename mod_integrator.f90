! Module that exposes the integrator that transforms the state_in to the state_out

module mod_integrator
  
  private

  public :: do_timestep
  
contains
  
  subroutine do_timestep (state_in, state_out)
    use mod_constants
    use mod_timestep
    implicit none
    
    type (state), intent(in), dimension(:)        :: state_in
    type (state), intent(out), dimension(:)       :: state_out

    real(kind = x_precision)                      :: dt, dx2
    integer                                       :: i, info

    real(kind = x_precision), dimension(n_cell)   :: S, T, nu, diag
    real(kind = x_precision), dimension(n_cell-1) :: diag_low, diag_up
    
    ! Get the timestep
    call timestep(state_in, dt)

    dx2 = dx**2
    
    ! Create the diagonals and copy T and S
    do i = 1, n_cell
       T  = state_in(i)%T
       S  = state_in(i)%S

       diag(i) = 1/dt + 2/dx2 * (state_in(i)%nu / state_in(i)%x)
       
       if (i < n_step) then
          diag_up(i) = -1/dx2 * state_in(i+1)%nu / state_in(i)%x**2
       else ! i = n_cell
          diag(i)       = state_in(i)%nu
          diag_low(i-1) = state_in(i-1)%nu
       end if
       
       if (i > 1) then
          diag_low(i) = -1/dx2 * state_in(i-1)%nu / state_in(i)%x**2
       else ! i = 1
          diag(1) = 1
       end if
    end do
    
    ! Solving for S
    call dgtsv(n_cell, 1, diag_low, diag, diag_up, S, n_cell, info)
    
    if (info /= 0) then
       print *, * 'Ooops, something bad happened!'
    end if

    ! Solve for T
    ! TODO

    ! Copy the result
    do i = 1, n_cell
       state_out(i)   = state_in(i)
       state_out(i)%S = S(i)
       state_out(i)%T = T(i)
    end do
  end subroutine do_timestep
  
end module mod_integrator

