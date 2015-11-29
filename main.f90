program black_hole_diffusion
  use mod_constants
  use mod_variables
  use mod_read_parameters
  use mod_s_curve
  use mod_integrator
  use mod_output

  implicit none

  integer                                     :: iteration, ios, i
  type(state)                                 :: s
  real(kind = x_precision)                    :: delta_S_max, delta_T_max, t
  real(kind = x_precision), dimension(n_cell) :: prev_S, S_crit
  real(kind = x_precision)                    :: t_V, t_T
  real(kind = x_precision)                    :: dt_V, dt_T
  logical                                     :: T_converged

  ! FIXME
  S_crit = 1.e99_x_precision
  ! FIXME
  delta_S_max = 1e-5
  ! FIXME I love it
  delta_T_max = 1e-5

  ! Initial time = 0
  t = 0._x_precision

  ! Read the parameters, generate state_0 and create adim state
  call get_parameters()

  ! Initiate the S_curve
  ! call s_curve(foo, bar)

  ! Copy the value of state_0 into state vector s
  s%T = CI%T_ci / state_0%T_0
  s%S = CI%Sig_ci / state_0%S_0 * x_state%x
  
  ! H_over_r become H*
  CI%H_over_r = CI%H_over_r * r_state%r
  open(15, file="CI.dat", status="replace", iostat=ios)
  
  if (ios /= 0) then
     stop "Error while opening output file."
  end if

  ! Save initial conditions
  write(15, *)'r x T* S* T Sigma H* '
  do i= 1, n_cell
     write(15, *)r_state%r(i), x_state%x(i), s%T(i), s%S(i), CI%T_ci(i), CI%Sig_ci(i), CI%H_over_r(i)
  enddo
  
  close(15)

  call compute_variables(s)

  ! Open unit 13 to write
  ! do it more safely
  open(13, file="output.dat", status="replace", iostat=ios)
  if (ios /= 0) then
     stop "Error while opening output file."
  end if

  ! Initialize prev_S
  ! we multiply delta to prevent the code from thinking it converged
  prev_S = 1.2*s%S
  
  ! FIXME
  ! write(*,*) state_0%temps_0
  t_T = params%t_T *state_0%temps_0 
  t_V = params%t_nu*state_0%temps_0
  
  !  t_V = 1.3e4_x_precision / state_0%temps_0 
  !  t_T = 0.72_x_precision / state_0%temps_0

  dt_V = t_V / 2._x_precision
  dt_T = t_T / 2._x_precision

  write(*,*) 'dt_T, dt_V:', dt_T, dt_V

  iteration = 0
  ! Start iterating
  do while (iteration < n_iterations)
     ! Check that S is at a fixed point
     if (maxval(abs((prev_S - s%S)/s%S)) > delta_S_max) then
        ! Check here that S < S_crit
        if (maxval(s%S - S_crit) > 0) then
           print*, 'Exiting because of S'
           ! Switch to explicit scheme
        else
           ! Switch to implicit scheme
           prev_S = s%S

           ! Integrate S
           call do_timestep_S(s, dt_V)
           ! Increment time
           t = t + dt_V

           ! print *, j, log10(s%T(50)*state_0%T_0), log10(s%S(50)*state_0%S_0)
           
           ! Recompute variables
           call compute_variables(s)

           iteration = iteration + 1
           ! Iterate while T hasn't converged
           T_converged = .false.
           do while (.not. T_converged)
              ! Integrate T
              ! print *, '# T integration'
              call do_timestep_T(s, dt_T, T_converged, delta_T_max)
              ! Increment time
              t = t + dt_T
              
              ! Recompute variables
              call compute_variables(s)
              if (mod(iteration, 10000) == 0) then
                 call snapshot(s, iteration, t, 13)
                 print*,'snapshot', iteration, t
              end if
              ! print *, j, log10(s%T(50)*state_0%T_0), log10(s%S(50)*state_0%S_0), t
              iteration = iteration + 1
           end do
           ! Output things here
           if (mod(iteration, 10000) == 0) then
              call snapshot(s, iteration, t, 13)
              print*,'snapshot', iteration, t
           end if
        end if
     else
       ! print*, 'Mdot kick'
        ! give a Mdot kick
     end if
     
     ! Increase M_0_dot if stalled
     ! TODO : increase M_dot_0
     ! state_0%M_dot_0 = state_0%M_dot_0 + 
  end do

  close(13)
  
end program black_hole_diffusion
