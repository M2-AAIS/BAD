program black_hole_diffusion
  use mod_constants
  use mod_read_parameters
  use mod_variables
  !use mod_s_curve
  use mod_integrator
  use mod_output
  use mod_timestep
  use mod_s_curve
  use mod_distance

  implicit none

  !--------------------------- Parameters for s_curve-----------------------
  real(kind = x_precision), parameter         :: eps_in = 1.e-8_x_precision  ! Precision required for the dichotomy
  real(kind = x_precision), parameter         :: Tmin   = 2.5e-2_x_precision
  real(kind = x_precision), parameter         :: Tmax   = 4.49e0_x_precision
  real(kind = x_precision), parameter         :: Smin   = 2.36e1_x_precision
  real(kind = x_precision), parameter         :: Smax   = 2.36e3_x_precision

  real(kind = x_precision), dimension(n_cell) :: temperature_c
  real(kind = x_precision), dimension(n_cell) :: sigma_c
  !-------------------------------------------------------------------------

  integer                                     :: iteration, ios, i
  type(state)                                 :: s
  real(kind = x_precision)                    :: delta_S_max, delta_T_max, t
  real(kind = x_precision), dimension(n_cell) :: prev_S, S_crit
 ! real(kind = x_precision)                    :: t_nu, t_T
  real(kind = x_precision), dimension(n_cell) :: dt_nu, dt_T
  logical                                     :: T_converged
  real (kind = x_precision), dimension(n_cell):: dist
  

  ! FIXME
  S_crit = 1.e99_x_precision
  ! FIXME
  delta_S_max = 1e-4
  ! FIXME I love it
  delta_T_max = 1e-3

  ! Read the parameters, generate state_0 and create adim state
  call get_parameters()
  call make_temperature(Tmin, Tmax)
  call set_conditions(eps_in, Smin, Smax)
  call curve(temperature_c, sigma_c)
  ! Initiate the S_curve
  ! call s_curve(foo, bar)

  ! Copy the value of state_0 into state vector s
  s%T = CI%T_ci / state_0%T_0
  !s%T = 10._x_precision * CI%T_ci / state_0%T_0
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
    write(15, *)r_state%r(i), x_state%x(i), s%T(i), s%S(i), s%T(i)*state_0%T_0, s%S(i) * state_0%s_0 / x_state%x(i), CI%H_over_r(i)
  enddo

  close(15)

  call compute_variables(s)
  call timestep (s,dt_T, dt_nu)

  ! Open unit 13 to write
  ! do it more safely
  open(13, file="output.dat", status="replace", iostat=ios)
  if (ios /= 0) then
    stop "Error while opening output file."
  end if

  ! Initialize prev_S
  ! we multiply delta to prevent the code from thinking it converged
  prev_S = 1.2*s%S

 ! t_T  = params%t_T  !* state_0%temps_0
 ! t_nu = params%t_nu !*  state_0%temps_0

 ! dt_nu = t_nu / cst_dt
 ! dt_T  = t_T / cst_dt

  write(*,*) 'dt_T_min, dt_nu_min:', minval(dt_T), minval(dt_nu)
  write(*,*) 'dt_T_max, dt_nu_max:', maxval(dt_T), maxval(dt_nu)

  ! Initial time = 0
  t = 0._x_precision

  ! Start
  iteration = 0

  ! Write initial state
  call snapshot(s, iteration, t, 13)

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
          
          call do_timestep_S(s, minval(dt_nu))
          
          ! Increment time, number of iterations
          t = t + minval(dt_nu)
          iteration = iteration + 1
          ! Output things here
          if (mod(iteration, output_freq) == 0) then
             call snapshot(s, iteration, t, 13)
             print*,'snapshot', iteration, t
          end if

          ! Iterate while T hasn't converged
          T_converged = .false.
          do while (.not. T_converged)
             ! Integrate T
             call do_timestep_T(s, minval(dt_T), T_converged, delta_T_max)
                         
             ! Increment time, number of iterations
             t = t + minval(dt_T)
             iteration = iteration + 1
             
             if (mod(iteration, output_freq) == 0) then
                call snapshot(s, iteration, t, 13)
                print*,'snapshot', iteration, t
             end if
             
          end do
          ! Recompute variables when the system is stable
          call compute_variables(s)
          call timestep (s,dt_T, dt_nu)
          call distance(s, dist, temperature_c, sigma_c)
       end if
    else
       iteration = n_iterations
       !s%Mdot(n_cell) = s%Mdot(n_cell) * 1.01_x_precision
    end if
 end do
  
 close(13)

end program black_hole_diffusion
