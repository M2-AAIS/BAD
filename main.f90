program black_hole_diffusion
  use mod_constants
  use mod_variables
  use mod_read_parameters
  use mod_s_curve
  use mod_integrator
  use mod_output

  implicit none

  integer                                     :: iteration, n, ios
  type(state)                                 :: s
  real(kind = x_precision)                    :: delta_S_max, delta_T_max, t
  real(kind = x_precision), dimension(n_cell) :: prev_S, prev_T, S_crit

  ! FIXME
  S_crit = 1.e99_x_precision
  ! FIXME
  delta_S_max = 1e-5
  ! FIXME I love it
  delta_T_max = 1e-5
  

  ! Initial time = 0
  t = 0

  ! Read the parameters 
  call get_parameters()
  
  ! Call init variable to create the adim state vector and
  ! generate the state_0
  call init_variable_0()

  ! Initiate the S_curve
  ! call s_curve(foo, bar)

  ! Copy the value of state_0 into state vector s
  do n = 1, n_cell
     s%x(n)     = (n-1) * params%dx + 3_x_precision ! x_min = 3
     s%Omega(n) = 1_x_precision/s%x(n)**3
     s%T(n)     = 1e-3_x_precision
     s%S(n)     = 1e-3_x_precision
  end do
  call compute_variables(s)

  ! Open unit 13 to write
  ! do it more safely
  open(13, file="output.dat", status="replace", iostat=ios)
  if (ios /= 0) then
     stop "Error while opening output file."
  end if

  ! Initialize prev_S and prev_T
  prev_S = s%S + 2*delta_S_max
  prev_T = s%T + 2*delta_T_max

  call snapshot(s, iteration, t, 13)
  
  ! Start iterating
  do iteration = 1, n_iterations
     ! Check that S is at a fixed point
     do while (maxval(abs(prev_S - s%S)) > delta_S_max)
        write(13, *) "While prev_S/s%S", iteration
        ! Check here that S < S_crit
        if (maxval(s%S - S_crit) > 0) then
           ! Switch to explicit scheme
        else
           ! Switch to implicit scheme
           prev_S = s%S

           ! Integrate S
           call do_timestep_S(s)

           write(13, *) "After S timestep"
           call snapshot(s, iteration, t, 13)
           
           ! Recompute variables
           call compute_variables(s)

           ! FIXME : increment time
 
           ! Iterate while T hasn't converged
           do while (maxval(abs(prev_T - s%T)) > delta_T_max)
              ! Check here that T < T_crit
              prev_T = s%T

              ! Integrate T
              call do_timestep_T(s)
              ! Recompute variables
              call compute_variables(s)
              write(13, *) "After integration of T integration"
              call snapshot(s, iteration, t, 13)

              ! FIXME : increment time
           end do
           ! Output things here
           write(13, *) "YOLO DUDE"
           call snapshot(s, iteration, t, 13)
        end if
     end do
     
     ! Increase M_0_dot if stalled
     ! TODO : increase M_dot_0
     ! state_0%M_dot_0 = state_0%M_dot_0 + 
  end do

  close(13)
  
end program black_hole_diffusion
