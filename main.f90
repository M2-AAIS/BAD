program black_hole_diffusion
  use mod_constants
  use mod_variables
  use mod_read_parameters
  use mod_s_curve
  use mod_integrator
  use mod_output

  implicit none

  integer                                     :: iteration, n, ios, j
  type(state)                                 :: s
  real(kind = x_precision)                    :: delta_S_max, delta_T_max, t
  real(kind = x_precision), dimension(n_cell) :: prev_S, prev_T, S_crit
  real(kind = x_precision)                    :: t_V, t_T  
  real(kind = x_precision)                    :: dt_V, dt_T

  ! FIXME
  S_crit = 1.e99_x_precision
  ! FIXME
  delta_S_max = 1e-2
  ! FIXME I love it
  delta_T_max = 1e-2

  ! Initial time = 0
  t = 0._x_precision

  ! Read the parameters 
  call get_parameters()
  
  ! Call init variable to create the adim state vector and
  ! generate the state_0
  call init_variable_0()

  ! Initiate the S_curve
  ! call s_curve(foo, bar)

  ! Copy the value of state_0 into state vector s
  do n = 1, n_cell
     s%T(n)       = 0.52_x_precision
     s%S(n)       = 1200._x_precision
  end do
  call compute_variables(s)

  ! Open unit 13 to write
  ! do it more safely
  open(13, file="output.dat", status="replace", iostat=ios)
  if (ios /= 0) then
     stop "Error while opening output file."
  end if

  ! Initialize prev_S and prev_T
  ! we multiply delta to prevent the code from thinking it converged
  prev_S = 1.2*s%S
  prev_T = 1.2*s%T
  
  ! FIXME
  ! write(*,*) state_0%temps_0
  t_T = params%t_T *state_0%temps_0 
  t_V = params%t_nu*state_0%temps_0
  write(*,*)'t_T, t_V:', t_T, t_V
  
  !  t_V = 1.3e4_x_precision / state_0%temps_0 
  !  t_T = 0.72_x_precision / state_0%temps_0

  dt_V = t_V / 2._x_precision
  dt_T = t_T / 2._x_precision
  
  ! call snapshot(s, iteration, t, 13)
  
  ! Start iterating
  do iteration = 1, 100!n_iterations
     ! Check that S is at a fixed point
     if (maxval(abs((prev_S - s%S)/s%S)) > delta_S_max) then
        ! Check here that S < S_crit
        if (maxval(s%S - S_crit) > 0) then
           ! Switch to explicit scheme
        else
           ! Switch to implicit scheme
           prev_S = s%S

           ! Integrate S
           print *, '#S integration'
           call snapshot(s, iteration, t, 13)
           call do_timestep_S(s, dt_V)
           ! Increment time
           t = t + dt_V

           print *, j, log10(s%T(50)*state_0%T_0), log10(s%S(50)*state_0%S_0)
           
           ! Recompute variables
           call compute_variables(s)

           ! Iterate while T hasn't converged
           j = 0
           do while (maxval(abs(prev_T - s%T)/s%T) > delta_T_max) 
              ! Check here that T < T_crit
              prev_T = s%T

              ! Integrate T
              print *, '# T integration'
              call do_timestep_T(s, dt_T)
              print *, j, log10(s%T(50)*state_0%T_0), log10(s%S(50)*state_0%S_0), t
              ! Increment time
              t = t + dt_T
              ! Recompute variables
              call compute_variables(s)
              if (mod(j, 50) == 0) then
                 call snapshot(s, iteration, t, 13)
              end if
              j = j+1

              ! FIXME : increment time
           end do
           ! Output things here
        end if
     else
        ! give a Mdot kick
     end if
     
     ! Increase M_0_dot if stalled
     ! TODO : increase M_dot_0
     ! state_0%M_dot_0 = state_0%M_dot_0 + 
  end do

  close(13)
  
end program black_hole_diffusion
