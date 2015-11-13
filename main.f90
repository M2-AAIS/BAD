program black_hole_diffusion
  use mod_constants
  use mod_variables
  use mod_read_parameters
  use mod_s_curve

  implicit none

  integer     :: iteration, n
  type(state) :: s

  call get_parameters()

  ! Call init variable to create the adim state vector and
  ! generate the state_0
  call init_variable_0()

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

  ! Start iterating
  do iteration = 1, n_iterations

  end do

end program black_hole_diffusion
