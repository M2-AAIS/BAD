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

  character(len=10)                      :: arg
  integer                                :: stp_value = 500000 ! select the stop iteration
  !  don't kill the job without stp_value
  real(x_precision), dimension(n_cell)   :: jnk1, jnk3
  character(len=100)                     :: line
  !--------------------------- Parameters for s_curve-----------------------
  real(x_precision), dimension(n_cell) :: T_c
  real(x_precision), dimension(n_cell) :: Sigma_c
  real(x_precision), dimension(n_cell) :: S_c
  !-------------------------------------------------------------------------

  integer                              :: start_iteration, iteration, ios, i, T_steps, S_steps
  type(state)                          :: s
  real(x_precision), dimension(n_cell) :: prev_S, prev_T, S_rel_diff
  real(x_precision)                    :: delta_S_max, delta_T_max, S_rel_diff_max
  real(x_precision)                    :: t
  real(x_precision)                    :: dt_nu, dt_T, pf
  integer :: p
  real(x_precision), dimension(n_cell) :: mean
  logical                              :: T_converged
  logical                              :: unstable, wasunstable

  real(x_precision), dimension(n_cell) :: dist_crit

  call getarg(1, arg)
  mean       = 0._x_precision
  !----------------------------------------------
  ! Convergence criteria for S and T
  !----------------------------------------------
  delta_S_max = 1e-8_x_precision
  delta_T_max = 1e-4_x_precision

  ! Read the parameters, generate state_0 and create adim state
  call get_parameters()
  call make_temperature()
  call set_conditions()
  call curve(T_c, Sigma_c)

  ! Conversion into S*_crit
  S_c = Sigma_c * x_state%x
  dist_crit = 200. * x_state%x

  !----------------------------------------------
  ! Set the initial conditions for S, T
  !----------------------------------------------
  open(15, file="CI.dat", status="replace", iostat=ios)
  if (ios /= 0) then
     stop "Error while opening output file."
  end if

  ! Save initial conditions
  write(15, '(7(A16))')'r', 'T', 'Sigma', 'H'

  do i = 1, n_cell
     write(15, '(7(e16.6e2))') r_state%r(i), IC%T(i), IC%Sigma(i), IC%H(i)
  enddo

  close(15)

  !-------------------------------------------------------------------
  !------------------------------Restart------------------------------
  !---------------------- Select n_iterations-------------------------
  !------Call ./simul load > let the job do until n_iterations--------
  !-----------------------Call ./simul restart------------------------
  !-------------------------------------------------------------------
  if (arg == 'restart' .or. arg == 'load') then
     if (arg == 'restart') then
        print*, 'Creating restart file and saving old output file to output.dat.prev'
        call system ("tail -n 259 output.dat > restart.dat")
        call system ("cp output.dat output.dat.prev")
     end if
     open(21, file="restart.dat", action='read', iostat=ios)

     if (ios /= 0) then
        stop "Error while opening output file."
     end if
     read(21,*)line, start_iteration
     read(21,*)line, t
     read(21,*)line
     do i = 1, n_cell
        read(21,*)jnk1(i), s%T(i), s%Mdot(i), jnk3(i), s%s(i)
     end do
     s%T = s%T / state_0%T_0
     s%S = s%S / state_0%S_0 * x_state%x
     s%Mdot = s%Mdot / state_0%Mdot_0
     close(21)

     if (maxval(s%T - T_c) < 0) then
       ! …or right from its critical S…
       unstable = maxval(s%S - S_c) > 0
     else
       ! …then we’re unstable.
       unstable = .true.
     end if

     iteration = start_iteration

  elseif (arg == 'start') then
     s%T = IC%T / state_0%T_0
     s%S = IC%Sigma / state_0%S_0 * x_state%x
     ! Initial time = 0
     t = 0._x_precision
     ! Start
     start_iteration = 0
     iteration = 0
     unstable = .false.
     s%Mdot(n_cell) = 1._x_precision
  else
     print*, 'Unsupported action "', arg, '". Call ./simul [start|load|restart].'
     stop
  end if

  !----------------------------------------------
  ! Do an initial computation of the variables
  !----------------------------------------------
  call compute_variables(s)

  !----------------------------------------------
  ! Open file to write output to
  !----------------------------------------------
  open(13, file="output.dat", status="replace", iostat=ios)
  if (ios /= 0) then
     stop "Error while opening output file."
  end if

  !----------------------------------------------
  ! Save initial snapshot
  !----------------------------------------------
  call snapshot(s, iteration, t, 13)

  !----------------------------------------------
  ! Compute initial timestep
  !----------------------------------------------
  call timestep_T(s, dt_T)
  call timestep_nu(s, dt_nu)
  write(*,*) 'dt_T_init, dt_nu_init:', dt_T, dt_nu

  !----------------------------------------------
  ! Start iterations
  !----------------------------------------------
  S_steps = 0
  do while (iteration < n_iterations)

     !--------------------------------------------
     ! Check if any point is in the unstable area
     !--------------------------------------------

     wasunstable = unstable
     ! If some point is above is critical temperature…
     if (maxval(s%T - T_c) < 0) then
       ! …or right from its critical S…
       unstable = maxval(s%S - S_c) > 0
     else
       ! …then we’re unstable.
       unstable = .true.
     end if

     if (unstable) then
        if (.not. wasunstable) then
           print*, 'Switched to explicit mode!'
        end if

        ! Do an explicit integration of both S and T over a thermic timestep
        call timestep_T(s, dt_T)
        dt_T = 0.01_x_precision * dt_T
        call do_timestep_S_exp(s, dt_T)
        call do_timestep_T(s, dt_T)

        ! Increase time, increment number of iterations
        t = t + dt_T
        iteration = iteration + 1

        if (mod(iteration, output_freq) == 0) then
           call snapshot(s, iteration, t, 13)
           if ((arg == 'load' .or. arg == 'restart') &
                .and. iteration-start_iteration > stp_value) then
              stop
           endif
           print*,'snapshot', iteration, t, 'exp', maxval(s%T - T_c), maxval(s%S - S_c), dt_T
        end if

     else
        S_steps = S_steps + 1 
        if (wasunstable) then
           print*, 'Switched back to implicit mode!'
           s%Mdot(n_cell) = s%Mdot(n_cell) / 4._x_precision
           print*, 'Mdot un-kick! un-YOLOOOOO', s%Mdot(n_cell)
        end if

        ! Integrate S, then integrate T until we’re back on the curve

        !----------------------------------------------
        ! S integration
        !----------------------------------------------

        call timestep_nu(s, dt_nu)

        ! Take into account the distance to the critical point
        pf =  pre_factor(s%S, S_c, dist_crit)
        dt_nu = dt_nu * pf

        ! Do a single S integration
        prev_S = s%S ! We need to keep the old value for comparison
        call do_timestep_S_imp(s, dt_nu)

        ! Increase time, increment number of iterations
        t = t + dt_nu
        iteration = iteration + 1

        if (mod(iteration, output_freq) == 0 ) then
           if ((arg == 'load' .or. arg == 'restart') &
                .and. iteration-start_iteration > stp_value) then
              stop
           endif
           call snapshot(s, iteration, t, 13)
           print*,'snapshot', iteration, t, 'S', dt_nu, dt_T, pf
        end if

        !----------------------------------------------
        ! T integrations
        ! Iterate while T hasn't converged
        !----------------------------------------------

        T_converged = .false.

        call compute_variables(s)
        
        T_steps = 0
        do while (.not. T_converged)
           T_steps = T_steps + 1
           if (T_steps >= 20 .and. mod(T_steps, 20) == 0) then
              print*, 'W: T is taking time to converge, already', T_steps, 'iterations'
           end if

           call compute_variables(s)
           call timestep_T(s, dt_T)

           ! Do a single T integration
           prev_T = s%T ! We need to keep the old value for comparison
           call do_timestep_T(s, dt_T)

           ! Increase time, increment number of iterations
           t = t + dt_T
           iteration = iteration + 1

           ! Do a snapshot
           if (mod(iteration, output_freq) == 0) then
              if ((arg == 'load' .or. arg == 'restart') &
                   .and. iteration-start_iteration > stp_value) then
                 stop
              endif
              call snapshot(s, iteration, t, 13)
              print*,'snapshot', iteration, t, 'T', dt_nu, dt_T, pf
           end if

           ! Check if T converged everywhere
           T_converged = maxval(abs(prev_T - s%T)/s%T) / dt_T < delta_T_max

        end do

        !----------------------------------------------
        ! Mdot kick
        ! Increase Mdot at r_max if S is stalled
        !
        ! we give an mdot kick when the relative variation is below delta_S_max
        ! and when the mean over last 10 iterations does not evolve anymore
        !----------------------------------------------
        S_rel_diff = (prev_S - s%S) / s%S / dt_nu
        S_rel_diff_max = maxval(dabs(S_rel_diff))

        p = mod(S_steps, mean_size)
        if (S_rel_diff_max < delta_S_max * mean_size) then
           ! - - - - - - - - - - - - - - - - - - - - - - -
           ! compute the rolling mean of the latest S when getting
           ! close to the convergence
           ! - - - - - - - - - - - - - - - - - - - - - - -
           if (S_rel_diff_max < delta_S_max) then
              s%Mdot(n_cell) = s%Mdot(n_cell) * 2._x_precision
              print*, 'Converged lately!'
              print*, 'Mdot kick! YOLOOOOO', s%Mdot(n_cell), S_steps
           else
              mean = S_rel_diff / mean_size + &
                   (mean_size - 1._x_precision) / mean_size * mean
              
              if (maxval(dabs(mean)) < delta_S_max) then
                 s%Mdot(n_cell) = s%Mdot(n_cell) * 2._x_precision
                 print*, 'Converged in mean!', maxval(dabs(mean)), S_rel_diff_max
                 print*, 'Mdot kick! YOLOOOOO', s%Mdot(n_cell)
              else
                 if (mod(iteration, output_freq) < 10) then
                    print*, maxval(dabs(mean)), S_rel_diff_max
                 end if
              end if
           end if
        end if
        
     end if

     ! Recompute variables when the system is stable
     call compute_variables(s)

  end do

  close(13)

end program
