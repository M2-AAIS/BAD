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
  real(kind = x_precision), parameter         :: eps_in = 1.e-7_x_precision  ! Precision required for the dichotomy
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
  real(kind = x_precision)                    :: t_nu, t_T
  real(kind = x_precision)                    :: dt_nu, dt_T, dt_tmp
  logical                                     :: T_converged
  real (kind = x_precision), dimension(n_cell):: dist, dist_crit
  

  ! FIXME
  S_crit = 1.e99_x_precision

  !----------------------------------------------
  ! Convergence criterion for S and T
  !----------------------------------------------
  delta_S_max = 1e-6
  delta_T_max = 1e-4

  ! Read the parameters, generate state_0 and create adim state
  call get_parameters()
  call make_temperature(Tmin, Tmax)
  call set_conditions(eps_in, Smin, Smax)
  call curve(temperature_c, sigma_c)

  ! Pass to T*_crit and S*_crit
  temperature_c  = temperature_c / state_0%T_0
  sigma_c = sigma_c / state_0%S_0 * x_state%x
  
  ! Initiate the S_curve
  ! call s_curve(foo, bar)

  !----------------------------------------------
  ! Set the initial conditions for S, T
  !----------------------------------------------
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
    write(15, *)r_state%r(i), x_state%x(i), s%T(i), s%S(i), s%T(i)*state_0%T_0, s%S(i) * state_0%s_0 / x_state%x(i), CI%H_over_r(i)
  enddo
  close(15)

  !----------------------------------------------
  ! Do an initial computation of the variables
  !----------------------------------------------
  call compute_variables(s)
  call timestep (s,dt_T, dt_nu)

  !----------------------------------------------
  ! Open file to write output to
  !----------------------------------------------
  open(13, file="output.dat", status="replace", iostat=ios)
  if (ios /= 0) then
    stop "Error while opening output file."
  end if

  !----------------------------------------------
  ! Initialize t, dt and iteration counter
  !----------------------------------------------
  t_T  = params%t_T  !* state_0%temps_0
  t_nu = params%t_nu !*  state_0%temps_0

  dt_nu = t_nu * 10._x_precision
  dt_T  = t_T * 10._x_precision
  
  ! Initial time = 0
  t = 0._x_precision

  ! Start
  iteration = 0

  ! Write initial state
  call snapshot(s, iteration, t, 13)

  ! Start iterations
  do while (iteration < n_iterations)
     ! Check that we are not close to S_critical
     if (maxval(s%S - S_crit) > 0) then
        print*, 'Exiting because of S'
        ! Switch to explicit scheme
     else
        ! Implicit schem
        prev_S = s%S

        ! Do a single S integration
        call do_timestep_S(s, dt_nu)

        ! Increment time, number of iterations
        t = t + dt_nu
        iteration = iteration + 1

        ! Do a snapshot
        if (mod(iteration, output_freq) == 0) then
           call snapshot(s, iteration, t, 13)
           print*,'snapshot', iteration, t
        end if

        ! Do T integrations
        T_converged = .false.
        do while (.not. T_converged)
           dt_tmp = dt_T
           ! Do a single T integration
           call do_timestep_T(s, dt_tmp, T_converged, delta_T_max)

           ! Increment time, number of iterations
           t = t + dt_tmp
           iteration = iteration + 1

           ! Do a snapshot
           if (mod(iteration, output_freq) == 0) then
              call snapshot(s, iteration, t, 13)
              print*,'snapshot', iteration, t
           end if
        end do

     end if
     if (maxval(abs((prev_S - s%S)/s%S)) < delta_S_max) then
        ! Mdot kick
        
        iteration = n_iterations
        !s%Mdot(n_cell) = s%Mdot(n_cell) * 1.01_x_precision
     end if

     ! Recompute variables and loop
     call compute_variables(s)
  end do

  close(13)

end program
