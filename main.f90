program black_hole_diffusion
  use mod_constants
  use mod_read_parameters
  use mod_variables
  use mod_integrator
  use mod_output
  use mod_timestep
  use mod_s_curve
  use mod_distance

  implicit none

  !--------------------------- Parameters for s_curve-----------------------
  real(x_precision), parameter         :: eps_in = 1.0e-7_x_precision ! Precision required for the dichotomy
  real(x_precision), parameter         :: Tmin   = 2.5e-2_x_precision / 2.
  real(x_precision), parameter         :: Tmax   = 4.49e0_x_precision * 2.
  real(x_precision), parameter         :: Smin   = 2.36e1_x_precision / 5.
  real(x_precision), parameter         :: Smax   = 2.36e3_x_precision * 5.

  real(x_precision), dimension(n_cell) :: temperature_c
  real(x_precision), dimension(n_cell) :: sigma_c
  real(x_precision), dimension(n_cell) :: S_c
  !-------------------------------------------------------------------------

  integer                              :: iteration, ios, i
  type(state)                          :: s
  real(x_precision)                    :: delta_S_max, delta_T_max, t
  real(x_precision), dimension(n_cell) :: prev_S
  real(x_precision), dimension(n_cell) :: dt_nu, dt_T
  real(x_precision)                    :: min_dt_nu, min_dt_T
  logical                              :: T_converged

  real(x_precision), dimension(n_cell) :: dist_crit
  real(x_precision)                    :: dt_pre_factor

  !----------------------------------------------
  ! Convergence criteria for S and T
  !----------------------------------------------
  delta_S_max = 1e-8
  delta_T_max = 1e-4

  ! Read the parameters, generate state_0 and create adim state
  call get_parameters()
  call make_temperature(Tmin, Tmax)
  call set_conditions(eps_in, Smin, Smax)
  call curve(temperature_c, sigma_c)

  ! Conversion into S*_crit
  S_c = sigma_c * x_state%x
  dist_crit = 400. * x_state%x

  ! Initiate the S_curve
  ! call s_curve(foo, bar)

  !----------------------------------------------
  ! Set the initial conditions for S, T
  !----------------------------------------------
  s%T = 8._x_precision*CI%T_ci / state_0%T_0
  s%S = CI%Sig_ci / state_0%S_0 * x_state%x

  ! H_over_r becoming H*
  CI%H_over_r = CI%H_over_r * r_state%r
  open(15, file="CI.dat", status="replace", iostat=ios)

  if (ios /= 0) then
     stop "Error while opening output file."
  end if

  ! Save initial conditions
  write(15, '(7(A16))')'r', 'x', 'T*', 'S*', 'T', 'Sigma', 'H*'
  do i= 1, n_cell
     write(15, '(7(e16.6e2))')r_state%r(i), x_state%x(i), s%T(i), s%S(i),&
                                            s%T(i)*state_0%T_0, s%S(i) * state_0%s_0 / x_state%x(i), CI%H_over_r(i)
  enddo
  close(15)

  !----------------------------------------------
  ! Do an initial computation of the variables
  !----------------------------------------------
  call compute_variables(s)
  call timestep (s, dt_T, dt_nu)
  min_dt_T = minval(dt_T)
  min_dt_nu = minval(dt_nu)

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

  ! Initialize prev_S
  ! we multiply delta to prevent the code from thinking it converged
  prev_S = 1.2*s%S

  ! t_T  = params%t_T  !* state_0%temps_0
  ! t_nu = params%t_nu !*  state_0%temps_0

  ! dt_nu = t_nu / cst_dt
  ! dt_T  = t_T / cst_dt

  write(*,*) 'dt_T_min, dt_nu_min:', min_dt_T, min_dt_nu
  write(*,*) 'dt_T_max, dt_nu_max:', maxval(dt_T), maxval(dt_nu)

  ! Initial time = 0
  t = 0._x_precision

  ! Start
  iteration = 0

  !----------------------------------------------
  ! Save initial snapshot
  !----------------------------------------------
  call snapshot(s, iteration, t, 13)

  !----------------------------------------------
  ! Start iterations
  !----------------------------------------------
  do while (iteration < n_iterations)
     prev_S = s%S
     
     ! Check that we are not close to S_critical
     if (1. - maxval(s%S/S_c) < 1.e-2) then
        call do_timestep_S_exp(s, min_dt_T*8.e-3)
        call do_timestep_T(s, min_dt_T*1.e-2, T_converged, delta_T_max)
        t = t + 2*min_dt_T*1e-2
        call compute_variables(s)
        call timestep (s, dt_T, dt_nu)

        if (mod(iteration, output_freq) == 0 .or. &
             (iteration > 830000 .and. mod(iteration, 10) == 0)) then
           call snapshot(s, iteration, t, 13)
           print*,'snapshot', iteration, t, 1 - maxval(s%S/S_c), dt_pre_factor
        end if
        iteration = iteration + 1


        ! Switch to explicit scheme
     else
        !----------------------------------------------
        ! S integration
        !----------------------------------------------

        ! Do a single S integration
        call do_timestep_S_imp(s, min_dt_nu)

        ! Update time, number of iterations
        t = t + min_dt_nu
        iteration = iteration + 1

        ! Do a snapshot
        if (mod(iteration, output_freq) == 0) then
           call snapshot(s, iteration, t, 13)
           print*,'snapshot', iteration, t
        end if

        !----------------------------------------------
        ! T integrations
        ! iterate while t hasn't converged
        !----------------------------------------------
        T_converged = .false.
        do while (.not. T_converged)
           ! Recompute the dt pre factor
           dt_pre_factor = pre_factor(s, S_c, dist_crit)

           ! Do a single T integration
           call do_timestep_T(s, min_dt_T, T_converged, delta_T_max)

           ! Increment time, number of iterations
           t = t + min_dt_T
           iteration = iteration + 1

           ! Do a snapshot
           if (mod(iteration, output_freq) == 0) then
              call snapshot(s, iteration, t, 13)
              print*,'snapshot', iteration, t
           end if
        end do

        ! Recompute variables when the system is stable
        call compute_variables(s)
        call timestep (s, dt_T, dt_nu)
        dt_pre_factor = pre_factor(s, S_c, dist_crit)
        min_dt_T = minval(dt_T) * dt_pre_factor
        min_dt_nu = minval(dt_nu) * dt_pre_factor


     end if

     !----------------------------------------------
     ! Mdot kick
     ! increase Mdot at the boundaries if S is stalled
     !----------------------------------------------
     if (maxval(abs((prev_S - s%S)/s%S)) < delta_S_max) then
        ! Mdot kick
        params%Mdot_kick_factor = params%Mdot_kick_factor * 1.5
        print*, 'Mdot kick!', params%Mdot_kick_factor
        iteration = iteration + 1
     end if
  end do

  close(13)

end program
