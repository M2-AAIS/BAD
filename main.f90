program black_hole_diffusion
  use mod_constants
  use mod_variables
  use mod_read_parameters
  use mod_s_curve

  implicit none

  integer     :: iteration, n, ios
  type(state) :: s
  real(kind = x_precision) :: delta_S_max, delta_T_max, t
  real(kind = x_precision), dimension(n_cell) :: prev_S, prev_T, S_crit

  ! FIXME
  S_crit = 1e15
  ! FIXME
  delta_S_max = 1e-5
  ! FIXME I love it
  delta_T_max = 1e-5
  

  ! Initial time = 0
  t = 0

  call get_parameters()

  ! Call init variable to create the adim state vector and
  ! generate the state_0
  call init_variable_0()

  ! Initiate the S_curve
  ! call s_curve(foo, bar)

  ! Copy the value of state_0 into state vector s
  do n = 1, n_cell
     s%Omega(n) = state_0%Omega_0
     s%x(n)     = state_0%x_0
     s%nu(n)    = state_0%nu_0
     s%v(n)     = state_0%v_0
     s%T(n)     = state_0%T_0
     s%P_rad(n) = state_0%P_rad_0
     s%P_gaz(n) = state_0%P_gaz_0
     ! FIXME
     ! s%beta(n)  = state_0%beta_0
     s%cs(n)    = state_0%cs_0
     s%H(n)     = state_0%H_0
     s%rho(n)   = state_0%rho_0
     s%S(n)     = state_0%S_0
     s%Fz(n)    = state_0%Fz_0
     s%M_dot(n) = state_0%M_dot_0
     s%Cv(n)    = state_0%Cv_0
  end do

  ! Open unit 13 to write
  ! do it more safely
  open(13, file="output.dat", status="replace", iostat=ios)
  if (ios /= 0) then
     stop "Error while opening output file."
  end if
  
  ! Start iterating
  do iteration = 1, n_iterations
     ! Check that S is at a fixed point
     do while (maxval(abs(prev_S - s%S)) > delta_S_max)
        ! Check here that S < S_crit
        if (maxval(s%S - S_crit) > 0) then
           ! Switch to explicit scheme
        else
           ! Switch to implicit scheme
           prev_S = s%S

           ! Integrate S
           call do_timestep_S(s)
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

              ! FIXME : increment time
           end do
           ! Output things here
           call snapshot(s, iteration, t, 13)
        end if
     end do
     
     ! Increase M_0_dot if stalled
     ! TODO : increase M_dot_0
     ! state_0%M_dot_0 = state_0%M_dot_0 + 
  end do

  close(13)
  
end program black_hole_diffusion
