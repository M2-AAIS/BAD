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
  integer                                :: stp_value = 50000 ! select the stop iteration 
  !  don't kill the job without stp_value
  real(x_precision), dimension(n_cell)   :: jnk1, jnk3
  character(len=100)                     :: line
  !--------------------------- Parameters for s_curve-----------------------
  real(x_precision), dimension(n_cell) :: T_c
  real(x_precision), dimension(n_cell) :: Sigma_c
  real(x_precision), dimension(n_cell) :: S_c
  !-------------------------------------------------------------------------

  integer                              :: iteration, ios, i
  type(state)                          :: s
  real(x_precision), dimension(n_cell) :: prev_S, prev_T
  real(x_precision)                    :: delta_S_max, delta_T_max
  real(x_precision)                    :: t
  real(x_precision)                    :: dt_nu, dt_T
  logical                              :: T_converged
  logical                              :: unstable

  real(x_precision), dimension(n_cell) :: dist_crit

  call getarg(1,arg)

  !----------------------------------------------
  ! Convergence criteria for S and T
  !----------------------------------------------
  delta_S_max = 1e-8
  delta_T_max = 1e-4

  ! Read the parameters, generate state_0 and create adim state
  call get_parameters()
  call make_temperature()
  call set_conditions()
  call curve(T_c, Sigma_c)

  ! Conversion into S*_crit
  S_c = Sigma_c * x_state%x
  dist_crit = 100. * x_state%x

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
  !--------------------------Process restart--------------------------
  !------------ Select stp_value to stop the firt execute-------------
  !--Call ./simul start > let the job do until stp_value for iteration
  !-----------------------Call ./simul restart------------------------
  !-------------------------------------------------------------------
  if (arg == 'restart') then  
     call system ("tail -n 259 output.dat > restart.dat" )
     open(21, file="restart.dat", action='read', iostat=ios)
     if (ios /= 0) then
        stop "Error while opening output file."
     end if
     read(21,*)line, iteration
     read(21,*)line, t
     read(21,*)line
     do i = 1, n_cell      
        read(21,*)jnk1(i), s%T(i), jnk3(i), s%s(i)
     end do
     s%T = s%T / state_0%T_0
     s%S = s%S / state_0%S_0 * x_state%x
     close(21)
     ! restart iteration FIXEME
  elseif (arg == 'load') then
     s%T = IC%T / state_0%T_0
     s%S = IC%Sigma / state_0%S_0 * x_state%x
     ! Initial time = 0
     t = 0._x_precision
     ! Start
     iteration = 0
  else
     print*, 'Unsupported action "', arg, '". Call ./simul [start|restart].'
     stop
  end if

  !----------------------------------------------
  ! Do an initial computation of the variables
  !----------------------------------------------
  s%Mdot(n_cell) = 1._x_precision
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
  do while (iteration < n_iterations)

     !--------------------------------------------
     ! Check if any point is in the unstable area
     !--------------------------------------------

     ! If some point is above is critical temperature…
     if (maxval(s%T - T_c) < 0) then
       ! …or right from is critical S…
       unstable = maxval(s%S - S_c) > 0
     else
       ! …then we’re unstable.
       unstable = .true.
     end if

     if (unstable) then

        ! Do an explicit integration of both S and T over a thermic timestep
        call timestep_T(s, dt_T)
        call do_timestep_S_exp(s, dt_T)
        call do_timestep_T(s, dt_T)

        ! Increase time, increment number of iterations
        t = t + dt_T
        iteration = iteration + 1

        if (mod(iteration, output_freq) == 0) then
           call snapshot(s, iteration, t, 13)
           if (arg == 'load' .and. iteration > stp_value) then
              stop
           endif
           print*,'snapshot', iteration, t, 'exp', maxval(s%T - T_c), maxval(s%S - S_c), dt_T
        end if

     else

        ! Integrate S, then integrate T until we’re back on the curve

        !----------------------------------------------
        ! S integration
        !----------------------------------------------

        call timestep_nu(s, dt_nu)

        ! Take into account the distance to the critical point
        dt_nu = dt_nu * pre_factor(s%S, S_c, dist_crit)

        ! Do a single S integration
        prev_S = s%S ! We need to keep the old value for comparison
        call do_timestep_S_imp(s, dt_nu)

        ! Increase time, increment number of iterations
        t = t + dt_nu
        iteration = iteration + 1

        if (mod(iteration, output_freq) == 0 ) then
           if (arg == 'load' .and. iteration > stp_value) then
              stop
           endif
           call snapshot(s, iteration, t, 13)
           print*,'snapshot', iteration, t, 'S', dt_nu
        end if

        !----------------------------------------------
        ! T integrations
        ! Iterate while T hasn't converged
        !----------------------------------------------

        T_converged = .false.

        do while (.not. T_converged)

           call timestep_T(s, dt_T)

           ! Do a single T integration
           prev_T = s%T ! We need to keep the old value for comparison
           call do_timestep_T(s, dt_T)

           ! Increase time, increment number of iterations
           t = t + dt_T
           iteration = iteration + 1

           ! Do a snapshot
           if (mod(iteration, output_freq) == 0) then
              if (arg == 'load' .and. iteration > stp_value) then
                 stop
              endif
              call snapshot(s, iteration, t, 13)
              print*,'snapshot', iteration, t, 'T', dt_T
           end if

           ! Check if T converged everywhere
           T_converged = maxval(abs(prev_T - s%T)/s%T) < delta_T_max

        end do

        !----------------------------------------------
        ! Mdot kick
        ! Increase Mdot at r_max if S is stalled
        !----------------------------------------------
        if (maxval(abs(prev_S - s%S)/s%S) / dt_nu < delta_S_max) then
           s%Mdot(n_cell) = s%Mdot(n_cell) * 2._x_precision
           print*, 'Mdot kick! YOLOOOOO', s%Mdot(n_cell)
        end if

     end if

     ! Recompute variables when the system is stable
     call compute_variables(s)

  end do

  close(13)

end program
