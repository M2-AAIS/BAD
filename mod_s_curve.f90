module mod_s_curve
  use mod_constants
  use mod_read_parameters
  implicit none

  real(x_precision), dimension(nb_it) :: temperature ! Temperatures used for the dichotomy
  real(x_precision), dimension(nb_it) :: temp_real   ! Dimensioned version
  real(x_precision)                   :: eps         ! Precision required for the dichotomy
  real(x_precision)                   :: S_min       ! Minimum surface density limit
  real(x_precision)                   :: S_max       ! Maximum surface density limit

contains
  !-------------------------------------------------------------------------
  ! SUBROUTINES :
  !               make_temperature      : : Creates the required temperature array
  !                                     <-- input : Lower and upper limit for the temperature coordinate of the S curve

  !               set_conditions        : : Set conditions for the resolution
  !                                     <-- input : Resolution threshold, lower and upper limits for the Sigma coordinate of the S curve

  !               curve                 : : fake main
  !                                         --> output : Coordinates (Temperature, Sigma) of the first critical point

  !               first_critical_point  : :
  !                                     <-- input : Sigma for the optically thick medium
  !                                         --> output : (Temperature, Sigma) of the first critical point + index of the point

  !               second_critical_point : :
  !                                     <-- input : Sigma for the optically thin medium + number of points + index of the thick critical point
  !                                         --> output : (Temperature, Sigma) of the second critical point + index of the point

  !               build_s_curve         : : combining both thin and thick media
  !                                     <-- input : Sigma for both media + index of the thin critical point
  !                                         --> output : Sigma

  !               variables             : : Computing all the variables for Q+ and Q-
  !                                     <-- input : cell index, Temperature, Sigma, + optical depth indicator (0 or 1)
  !                                         --> output : f = Q+ - Q- , tau_eff

  !               dichotomy             : : In order to find the result of the equation Q+ - Q- = 0
  !                                     <-- input : cell index, temperature index + optical depth indicator (0 or 1)
  !                                         --> output : Sigma between a computed range of surface density [Smin, Smax] for a given temperature

  !-------------------------------------------------------------------------

  subroutine make_temperature(Tmin, Tmax)
    implicit none

    real(x_precision), intent(in) :: Tmin
    real(x_precision), intent(in) :: Tmax
    !------------------------------------------------------------------------

    real(x_precision) :: dt
    integer           :: i

    dt = (Tmax - Tmin) / (nb_it - 1)

    do i = 1, nb_it
      temperature(i) = dt * (i-1) + Tmin
    enddo

    temp_real  = temperature * state_0%T_0

  end subroutine make_temperature


  subroutine set_conditions(eps_in, Smin, Smax)
    implicit none

    real(x_precision), intent(in) :: eps_in
    real(x_precision), intent(in) :: Smin
    real(x_precision), intent(in) :: Smax
    !------------------------------------------------------------------------

    eps   = eps_in
    S_min = Smin
    S_max = Smax

  end subroutine set_conditions


  subroutine curve(temperature_c, sigma_c)
    implicit none

    real(x_precision), dimension(n_cell), intent(out) :: temperature_c  ! Temperature for the critical point
    real(x_precision), dimension(n_cell), intent(out) :: sigma_c        ! Surface density for the critical point
    !------------------------------------------------------------------------
    integer                              :: i,k            ! Iteration counters

    ! For variables
    real(x_precision)                    :: f              ! Q+ - Q-
    real(x_precision)                    :: tau_eff        ! Effective optical depth
    integer                                     :: optical_depth  ! Indicator for the optical thickness

    ! For the separate curves given the thickness
    real(x_precision), dimension(nb_it)  :: sigma_t_thick
    real(x_precision), dimension(nb_it)  :: sigma_t_thin
    real(x_precision), dimension(nb_it)  :: sigma          ! Surface density for the S curve

    ! For the critical points
    integer                              :: index_fcp
    real(x_precision)                    :: sigma_c_thick
    real(x_precision)                    :: temp_c_thick
    real(x_precision), dimension(n_cell) :: sigma_thick
    real(x_precision), dimension(n_cell) :: temp_thick
    real(x_precision), dimension(n_cell) :: tau_thick      ! Optical thickness for the first critical point
    integer                              :: index_scp
    real(x_precision)                    :: sigma_c_thin
    real(x_precision)                    :: temp_c_thin
    real(x_precision), dimension(n_cell) :: sigma_thin
    real(x_precision), dimension(n_cell) :: temp_thin
    real(x_precision), dimension(n_cell) :: tau_thin       ! Optical thickness for the second critical point

    !------------------------------------------------------------------------
    ! Loop over values of r/x (radius)
    !------------------------------------------------------------------------
    do k = 1, n_cell
      ! For each point of the S curve
      do i = 1, nb_it
        ! Optical thick case (tau >= 1)
        optical_depth = 1

        ! S found with the dichotomy approach
        sigma_t_thick(i) = dichotomy(k, i, optical_depth)

        ! Optical thin case (tau < 1)
        optical_depth = 0

        ! S found with the dichotomy approach
        sigma_t_thin(i)  = dichotomy(k, i, optical_depth)
      enddo

      ! For the first critical point (the one at the right of the S shape)
      call first_critical_point(sigma_t_thick, index_fcp, sigma_c_thick, temp_c_thick)

      ! For the second critical point (the one at the left of the S shape)
      call second_critical_point(sigma_t_thick, sigma_t_thin, index_fcp, index_scp, sigma_c_thin, temp_c_thin)

      ! Combining the two separate curves into one that forms the S shape
      call build_s_curve(sigma_t_thick, sigma_t_thin, index_scp, sigma)

      ! Display the two critical points
      !call display_critical_points(sigma_c_thin, temp_c_thin, sigma_c_thick, temp_c_thick, k)

      ! Saving the DIMENSIONED values of Sigma and T for the S curve in a file
      call save_data(k, sigma)

      ! ADIMENSIONED critical T and Sigma [OUTPUT]
      temperature_c(k) = temp_c_thick
      sigma_c(k)       = sigma_c_thick

      ! DIMENSIONED critical T and Sigma
      temp_thick(k)  =  temp_c_thick * state_0%T_0
      sigma_thick(k) =  sigma_c_thick * state_0%S_0
      temp_thin(k)   =  temp_c_thin * state_0%T_0
      sigma_thin(k)  =  sigma_c_thin * state_0%S_0

      ! Computing the tau_eff values for both critical points
      call variables(k, temp_c_thick, sigma_c_thick, f, optical_depth, tau_eff)
      tau_thick(k)   = tau_eff
      call variables(k, temp_c_thin, sigma_c_thin, f, optical_depth, tau_eff)
      tau_thin(k)    = tau_eff

    enddo

    call write_critical_points(sigma_thin, temp_thin, tau_thin, sigma_thick, temp_thick, tau_thick)

  end subroutine curve


  !-------------------------------------------------------------------------
  ! Subroutine in order to find the first critical point
  !-------------------------------------------------------------------------
  subroutine first_critical_point(sigma_real_thick, index_fcp, sigma_c_thick, temp_c_thick)
    implicit none

    real(x_precision), dimension(nb_it), intent(in)  :: sigma_real_thick

    integer,                             intent(out) :: index_fcp
    real(x_precision),                   intent(out) :: sigma_c_thick
    real(x_precision),                   intent(out) :: temp_c_thick
    !-----------------------------------------------------------------------

    integer :: i

    i = 1

    do while (sigma_real_thick(i) == 0)
      i = i + 1
    enddo

    do while (sigma_real_thick(i) < sigma_real_thick(i+1))
      index_fcp     = i
      sigma_c_thick = sigma_real_thick(i)
      i             = i + 1
    enddo
    temp_c_thick  = temperature(i)

  end subroutine first_critical_point


  !-------------------------------------------------------------------------
  ! Subroutine in order to find the second critical point
  !-------------------------------------------------------------------------
  subroutine second_critical_point(sigma_t_thick, sigma_t_thin, index_fcp, index_scp, sigma_c_thin, temp_c_thin)
    implicit none

    real(x_precision), dimension(nb_it), intent(in)  :: sigma_t_thick
    real(x_precision), dimension(nb_it), intent(in)  :: sigma_t_thin
    integer,                             intent(in)  :: index_fcp

    integer,                             intent(out) :: index_scp
    real(x_precision),                   intent(out) :: sigma_c_thin
    real(x_precision),                   intent(out) :: temp_c_thin
    !-----------------------------------------------------------------------

    integer :: i

    i = max(1, index_fcp)

    ! The change occurs when the difference of sigma changes its sign
    do while (sigma_t_thick(i) > sigma_t_thin(i))
      i = i + 1
    end do

    index_scp    = i - 1
    sigma_c_thin = sigma_t_thin(index_scp)
    temp_c_thin  = temperature(index_scp)

  endsubroutine second_critical_point


  !-------------------------------------------------------------------------
  ! Subroutine in order to build the S curve
  !-------------------------------------------------------------------------
  subroutine build_s_curve(sigma_real_thick, sigma_real_thin, index_scp, sigma)
    implicit none

    real(x_precision), dimension(nb_it), intent(in)  :: sigma_real_thick
    real(x_precision), dimension(nb_it), intent(in)  :: sigma_real_thin
    integer,                             intent(in)  :: index_scp

    real(x_precision), dimension(nb_it), intent(out) :: sigma
    !-----------------------------------------------------------------------

    integer :: i

    do i = 1, index_scp
      sigma(i) = sigma_real_thick(i)
    enddo
    do i = index_scp + 1, nb_it
      sigma(i) = sigma_real_thin(i)
    enddo

  end subroutine build_s_curve


  !-------------------------------------------------------------------------
  ! Subroutine in order to display the critical points
  !-------------------------------------------------------------------------
  subroutine display_critical_points(sigma_c_thin, temp_c_thin, sigma_c_thick, temp_c_thick, k)
    implicit none

    real(x_precision), intent(in) :: sigma_c_thin
    real(x_precision), intent(in) :: temp_c_thin
    real(x_precision), intent(in) :: sigma_c_thick
    real(x_precision), intent(in) :: temp_c_thick
    integer,           intent(in) :: k

    !------------------------------------------------------------------------

    write(*,*)'**** Critical Point',k,'********'
    write(*,*)'****************************************'
    write(*,"(' Optically thin (T, Sigma):',1p,E12.4,4x,1p,E12.4)") temp_c_thin, sigma_c_thin
    write(*,"(' Optically thick (T, Sigma):',1p,E12.4,4x,1p,E12.4)") temp_c_thick, sigma_c_thick
    write(*,*)'****************************************'

  end subroutine display_critical_points


  !-------------------------------------------------------------------------
  ! Subroutine to save in a file the critical points
  !-------------------------------------------------------------------------
  subroutine write_critical_points(sigma_c_thin, temp_c_thin, tau_thin, sigma_c_thick, temp_c_thick, tau_thick)

    implicit none

    real(x_precision), dimension(n_cell), intent(in) :: sigma_c_thin
    real(x_precision), dimension(n_cell), intent(in) :: temp_c_thin
    real(x_precision), dimension(n_cell), intent(in) :: tau_thin
    real(x_precision), dimension(n_cell), intent(in) :: sigma_c_thick
    real(x_precision), dimension(n_cell), intent(in) :: temp_c_thick
    real(x_precision), dimension(n_cell), intent(in) :: tau_thick
    !-----------------------------------------------------------------------
    integer             :: i
    integer             :: fid_c = 11
    character(len = 64) :: fname_c

    !-----------------------------------------------------------------------

    fname_c = 'critical_points/file.dat'

    open(fid_c,file  = fname_c, status='unknown',action='readwrite')
    write(fid_c,'(7(A16))') 'Radius','Temp_thin','Sigma_thin','Tau_thin',&
                                     'Temp_thick','Sigma_thick','Tau_thick'

    do i = 1, n_cell
      write(fid_c,'(7(e16.6e2))') r_state%r(i), temp_c_thin(i), sigma_c_thin(i), tau_thin(i),&
                                                temp_c_thick(i), sigma_c_thick(i), tau_thick(i)
    end do

    close(fid_c)

  end subroutine write_critical_points

  !-------------------------------------------------------------------------
  ! Subroutine to save in a file the S curve
  !-------------------------------------------------------------------------
  subroutine save_data(k, sigma)
    implicit none

    integer,                             intent(in) :: k
    real(x_precision), dimension(nb_it), intent(in) :: sigma

    !------------------------------------------------------------------------
    ! For the output file
    character(len = 8)                  :: number_of_cell
    character(len = 64)                 :: fname_tot
    integer                             :: fid_tot
    integer                             :: ios

    real(x_precision), dimension(nb_it) :: sigma_real ! Dimensioned surface density
    real(x_precision)                   :: f
    real(x_precision)                   :: tau_eff
    integer                             :: optical_depth
    integer                             :: j
    !------------------------------------------------------------------------

    write(number_of_cell,'(I5.5)') k
    fid_tot = 22 + k
    fname_tot = 's_curves/Temperature_Sigma_'//trim(number_of_cell)//'_tot.dat'

    open(fid_tot,file  = fname_tot, action='write', status = 'replace', iostat = ios)
    if (ios /= 0) then
      write(*,*)"Error while opening the ", fname_tot," file."
      stop
    endif

    write(fid_tot,'(3(A16))') 'Surface_density', 'Temperature', 'Optical_depth'

    do j = 1, nb_it
      if (sigma(j) /= 0) then
        call variables(k, temperature(j), sigma(j), f, optical_depth, tau_eff)

        sigma_real(j) = sigma(j) * state_0%S_0

        write(fid_tot,fmt = '(3(e16.6e2))') sigma_real(j), temp_real(j), tau_eff
      endif
    enddo

    close(fid_tot)

  end subroutine save_data


  !-------------------------------------------------------------------------
  !Subroutine in order to compute variables H, rho, cs, nu, Q_plus, Q_minus,
  !K_ff, K_e, tau_eff, E_ff,Fz given T, Sigma and k, the position
  !------------------------------------------------------------------------
  subroutine variables(k, T, Sigma, f, optical_depth, tau_eff)
    implicit none

    real(x_precision), intent(in)  :: T,Sigma
    integer,           intent(in)  :: optical_depth
    integer,           intent(in)  :: k

    real(x_precision), intent(out) :: f,tau_eff
    !-----------------------------------------------------------------------
    real(x_precision) :: coeff_a,coeff_b,coeff_c
    real(x_precision) :: delta

    real(x_precision) :: Omega
    real(x_precision) :: H
    real(x_precision) :: rho
    real(x_precision) :: cs
    real(x_precision) :: nu
    real(x_precision) :: K_ff
    real(x_precision) :: K_e
    real(x_precision) :: E_ff
    real(x_precision) :: Fz
    real(x_precision) :: Q_plus
    real(x_precision) :: Q_minus
    !------------------------------------------------------------------------

    Omega   = x_state%Omega(k)

    coeff_a = (Omega * state_0%Omega_0)**2 * Sigma * state_0%S_0 / 2._x_precision
    coeff_b = (-1._x_precision/3._x_precision) * cst_rad * (T * state_0%T_0)**4 / state_0%H_0
    coeff_c = - params%RTM * T * Sigma * state_0%S_0 / (2._x_precision * state_0%H_0**2)
    delta   = coeff_b**2 - 4._x_precision * coeff_a * coeff_c
    H       = -0.5_x_precision * (coeff_b + sign(sqrt(delta),coeff_b)) / coeff_a

    rho     = Sigma / H
    cs      = Omega * H
    nu      = params%alpha * cs * H
    K_ff    = 6.13e22_x_precision * state_0%rho_0 * rho * (state_0%T_0 * T)**(-3.5_x_precision)
    K_e     = params%kappa_e
    E_ff    = 6.22e20_x_precision * (state_0%rho_0 * rho)**2 * sqrt(state_0%T_0 * T)
    tau_eff = 0.5_x_precision * sqrt(K_e * K_ff) * Sigma * state_0%S_0

    !-------------------------------------------------------------------------
    !Select the case for the optical depth to compute Fz
    !-------------------------------------------------------------------------

    select case(optical_depth)

    case(1)

      Fz = 2._x_precision * c**2 * T**4 /(27._x_precision * sqrt(3._x_precision) * (K_ff + K_e) * Sigma * state_0%S_0)

    case (0)

      Fz = 4._x_precision * state_0%H_0 * E_ff * H / (state_0%Omega_0 * state_0%S_0)

    end select

    Q_plus  = 3._x_precision  * state_0%H_0**2 * nu * (Omega * state_0%Omega_0)**2
    Q_minus = Fz / Sigma

    f = Q_plus - Q_minus

  end subroutine variables


  !-------------------------------------------------------------------------
  ! Dichotomic function in order to determine the change of sign in a given
  ! interval [Smin,Smax] with an epsilon precision
  !-------------------------------------------------------------------------
  real(x_precision) function dichotomy(k, i, optical_depth)
    implicit none

    integer, intent(in) :: k ! Position in the disk, to get r/Omega
    integer, intent(in) :: i ! Iteration counter, to get the temperature
    integer, intent(in) :: optical_depth
    !-------------------------------------------------------------------------
    integer           :: j
    real(x_precision) :: T ! Temperature used for the dichotomy
    real(x_precision) :: tau_eff
    real(x_precision) :: Smin     ! Lowest point
    real(x_precision) :: Smax     ! Highest point
    real(x_precision) :: S_center ! Intermediate point
    real(x_precision) :: f_min    ! Q+ − Q− at the lowest point
    real(x_precision) :: f_max    ! Q+ − Q− at the highest point
    real(x_precision) :: f_center ! Q+ − Q− at the intermediate point
    !-------------------------------------------------------------------------

    T = temperature(i)

    Smin     = S_min
    Smax     = S_max

    S_center = (Smin + Smax) / 2._x_precision

    j = 0

    call variables(k, T, Smin, f_min, optical_depth, tau_eff)

    call variables(k, T, Smax, f_max, optical_depth, tau_eff)

    if ( f_max * f_min > 0.) then

      dichotomy = 0

    else if( f_max * f_min < 0.) then

      do while (dabs(Smax - Smin) >= eps .and. j < max_it)

        call variables(k, T, Smin, f_min, optical_depth, tau_eff)

        call variables(k, T, Smax, f_max, optical_depth, tau_eff)

        call variables(k, T, S_center, f_center, optical_depth, tau_eff)

        if (f_min * f_center > 0.) then

          Smin = S_center

        else if (f_max * f_center > 0.) then

          Smax = S_center

        endif

        S_center = (Smin + Smax) / 2._x_precision
        j = j + 1

      end do

      dichotomy = S_center

    endif

  end function dichotomy

end module mod_s_curve
