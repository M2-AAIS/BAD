module mod_s_curve
  use mod_constants
  use mod_read_parameters
  use mod_variables
  implicit none

contains
  !-------------------------------------------------------------------------
  ! SUBROUTINES :
  !               curve                 : : fake main
  !                                         --> output : 2 arrays for the coordinates (Temperature, S ) of the first critical point

  !               first_critical_point  : :
  !                                     <-- input : T and Sigma for the optically thick medium + number of points
  !                                         --> output : (Temperature, Sigma) of the first critical point + index of the point

  !               second_critical_point : :
  !                                     <-- input : T and Sigma for the optically thin medium + T for the thick one + number of points + index of the thick critical point + precision
  !                                         --> output : (Temperature, Sigma) of the second critical point + index of the point

  !               build_s_curve         : : combining both thin and thick media
  !                                     <-- input : T and Sigma for both media + index of the thin critical point
  !                                         --> output : array for (Temperature, Sigma)

  !               quadratic             : : For a positive real solution of a quadratic equation ax^2 + bx + c =0
  !                                     <-- input : factors a, b, c
  !                                         --> output : positive solution

  !               variables             : : Computing all the variables for Q+ and Q-
  !                                     <-- input : Temperature, Sigma, Omega + optical depth indicator (0 or 1) + more
  !                                         --> output : f = Q+ - Q-

  !               dichotomy             : : In order to find the result of the equation Q+ - Q- = 0
  !                                     <-- input : Temperature, Sigma, Omega + optical depth indicator (0 or 1) + precision
  !                                         --> output : Sigma between a computed range of surface density [Smin, Smax] for a given temperature

  !-------------------------------------------------------------------------

  subroutine curve(nb_it, max_it, eps, t_min, dt, s_min, s_max, temperature, s)
    implicit none

    integer                                   ,intent(in)  :: nb_it          ! Number of points for the S curve
    integer                                   ,intent(in)  :: max_it         ! Maximum number of dichotomy iterations
    real(kind = x_precision)                  ,intent(in)  :: eps            ! Precision required for the dichotomy
    real(kind = x_precision)                  ,intent(in)  :: t_min          ! Minimum temperature limit
    real(kind = x_precision)                  ,intent(in)  :: dt             ! Temperature step
    real(kind = x_precision)                  ,intent(in)  :: S_min          ! Minimum surface density limit
    real(kind = x_precision)                  ,intent(in)  :: S_max          ! Maximum surface density limit

    real(kind = x_precision),dimension(n_cell),intent(out) :: temperature    ! Temperature for the critical point
    real(kind = x_precision),dimension(n_cell),intent(out) :: s              ! Surface density for the critical point
    !------------------------------------------------------------------------
    integer                                                :: i,j,k

    ! For the temperature and s
    real(kind = x_precision)                               :: temp
    real(kind = x_precision)                               :: Smin
    real(kind = x_precision)                               :: Smax
    real(kind = x_precision)                               :: Sigma          ! Sigma

    real(kind = x_precision)                               :: omega          ! Angular velocity
    real(kind = x_precision)                               :: r              ! Radius
    real(kind = x_precision)                               :: f              ! Q+ - Q-
    real(kind = x_precision)                               :: tau_eff        ! Effective optical depth
    integer                                                :: optical_depth  ! Indicator for the optical thickness

    ! For the separate curves given the thickness
    real(kind = x_precision),dimension(nb_it)              :: temp_t_thick
    real(kind = x_precision),dimension(nb_it)              :: sigma_t_thick
    real(kind = x_precision),dimension(nb_it)              :: temp_t_thin
    real(kind = x_precision),dimension(nb_it)              :: sigma_t_thin
    real(kind = x_precision),dimension(nb_it)              :: temp_real      ! Temperature for the S curve
    real(kind = x_precision),dimension(nb_it)              :: sigma_real     ! Surface density for the S curve

    ! For the output file
    character(len = 8)                                     :: number_of_cell
    character(len = 64)                                    :: fname_tot
    integer                                                :: fid_tot
    integer                                                :: ios

    ! For the critical points
    integer                                                :: index_fcp
    real(kind = x_precision)                               :: sigma_c_thick
    real(kind = x_precision)                               :: temp_c_thick
    real(kind = x_precision),dimension(n_cell)             :: sigma_thick
    real(kind = x_precision),dimension(n_cell)             :: temp_thick
    real(kind = x_precision),dimension(nb_it)              :: tau_thick      ! Optical thickness for the first critical point
    integer                                                :: index_scp
    real(kind = x_precision)                               :: sigma_c_thin
    real(kind = x_precision)                               :: temp_c_thin
    real(kind = x_precision),dimension(n_cell)             :: sigma_thin
    real(kind = x_precision),dimension(n_cell)             :: temp_thin
    real(kind = x_precision),dimension(nb_it)              :: tau_thin       ! Optical thickness for the second critical point

    !------------------------------------------------------------------------
    ! Loop for n_cell values of x (radius)
    !------------------------------------------------------------------------

    do k = 1, n_cell
      r     = r_state%r(k)
      omega = x_state%Omega(k)

      ! For each point of the S curve
      do i = 1, nb_it
        temp          = dt * (i-1) + t_min   ! Chosing the temperature value

        ! Optical thick case (tau >= 1)
        optical_depth = 1

        ! S range
        Smin          = S_min
        Smax          = S_max

        ! S found with the dichotomy approach
        sigma         = dichotomy(Smin, Smax, max_it, eps, temp, omega, optical_depth)

        temp_t_thick(i)  = temp
        sigma_t_thick(i) = sigma


        ! Optical thin case (tau < 1)
        optical_depth = 0

        ! S range
        Smin          = S_min
        Smax          = S_max

        ! S found with the dichotomy approach
        sigma         = dichotomy(Smin, Smax, max_it, eps, temp, omega, optical_depth)

        temp_t_thin(i)   =  temp
        sigma_t_thin(i)  =  sigma

      enddo

      ! For the first critical point (the one at the right of the S shape)
      call first_critical_point(sigma_t_thick, temp_t_thick, index_fcp, sigma_c_thick, temp_c_thick, nb_it)

      ! For the second critical point (the one at the left of the S shape)
      call second_critical_point(sigma_t_thick, sigma_t_thin, temp_t_thick,&
       index_fcp, index_scp, sigma_c_thin, temp_c_thin, nb_it)

      ! Combining the two separate curves into one that forms the S shape
      call build_s_curve(sigma_t_thick, sigma_t_thin, temp_t_thick, nb_it, index_scp, sigma_real, temp_real)

      call display_critical_points(sigma_c_thin, temp_c_thin, sigma_c_thick, temp_c_thick, k)


      ! Saving the DIMENSIONED values of Sigma and T for the S curve in a file
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
        temp_real(j)  =  temp_real(j) 
        sigma_real(j) =  sigma_real(j) 

        call variables(temp_real(j), sigma_real(j), omega, f, optical_depth, tau_eff)

        temp_real(j)  =  temp_real(j) * state_0%T_0
        sigma_real(j) =  sigma_real(j) * state_0%S_0

        write(fid_tot,fmt = '(3(e16.6e2))') sigma_real(j), temp_real(j), tau_eff
      enddo
      
      close(fid_tot)

      ! ADIMENSIONED critical T and S [OUTPUT]

      temperature(k) = temp_c_thick
      s(k)           = sigma_c_thick 

      ! DIMENSIONED critical T and S
      temp_thick(k)  =  temp_c_thick * state_0%T_0
      sigma_thick(k) =  sigma_c_thick * state_0%S_0
      temp_thin(k)   =  temp_c_thin * state_0%T_0
      sigma_thin(k)  =  sigma_c_thin * state_0%S_0


      call variables(temp_c_thick, sigma_c_thick, omega, f, optical_depth, tau_eff)
      tau_thick(k)   = tau_eff
      call variables(temp_c_thin, sigma_c_thin, omega, f, optical_depth, tau_eff)
      tau_thin(k)    = tau_eff

    enddo

    call write_critical_points(n_cell, r_state%r, sigma_thin, temp_thin, sigma_thick, temp_thick,&
                                                  tau_thin, tau_thick)

  end subroutine curve


  !-------------------------------------------------------------------------
  ! Subroutine in order to find the first critical point
  !-------------------------------------------------------------------------
  subroutine first_critical_point(sigma_real_thick, temp_real_thick, index_fcp, sigma_c_thick, temp_c_thick, nb_it)
    implicit none

    integer                                  ,intent(in)  :: nb_it
    real(kind = x_precision),dimension(nb_it),intent(in)  :: sigma_real_thick
    real(kind = x_precision),dimension(nb_it),intent(in)  :: temp_real_thick

    integer                                  ,intent(out) :: index_fcp
    real(kind = x_precision)                 ,intent(out) :: sigma_c_thick
    real(kind = x_precision)                 ,intent(out) :: temp_c_thick
    !-----------------------------------------------------------------------

    integer                                               :: i

    i = 1
    do while (sigma_real_thick(i) < sigma_real_thick(i+1) .and. i < nb_it - 1)
      index_fcp     = i
      sigma_c_thick = sigma_real_thick(i)
      i             = i + 1
    enddo
    temp_c_thick  = temp_real_thick(i)

  end subroutine first_critical_point


  !-------------------------------------------------------------------------
  ! Subroutine in order to find the second critical point
  !-------------------------------------------------------------------------
  subroutine second_critical_point(sigma_real_thick, sigma_real_thin, temp_real_thick,&
         index_fcp, index_scp, sigma_c_thin, temp_c_thin, nb_it)
    implicit none

    integer                                  ,intent(in)  :: nb_it
    real(kind = x_precision),dimension(nb_it),intent(in)  :: sigma_real_thick
    real(kind = x_precision),dimension(nb_it),intent(in)  :: temp_real_thick
    real(kind = x_precision),dimension(nb_it),intent(in)  :: sigma_real_thin
    integer                                  ,intent(in)  :: index_fcp

    integer                                  ,intent(out) :: index_scp
    real(kind = x_precision)                 ,intent(out) :: sigma_c_thin
    real(kind = x_precision)                 ,intent(out) :: temp_c_thin
    !-----------------------------------------------------------------------

    integer                                               :: i


    i = max(1, index_fcp)

    ! The change occurs when the difference of sigma changes its sign
    do while ( (sigma_real_thick(i) > sigma_real_thin(i))  .and. i < nb_it)
      i = i + 1
    end do

    index_scp = i-1
    sigma_c_thin =  sigma_real_thin(index_scp)
    temp_c_thin = temp_real_thick(index_scp)

  endsubroutine second_critical_point


  !-------------------------------------------------------------------------
  ! Subroutine in order to build the S curve
  !-------------------------------------------------------------------------
  subroutine build_s_curve(sigma_real_thick, sigma_real_thin, temp_real_thick, nb_it, index_scp, sigma_real, temp_real)
    implicit none

    integer                                  ,intent(in)  :: nb_it
    integer                                  ,intent(in)  :: index_scp
    real(kind = x_precision),dimension(nb_it),intent(in)  :: sigma_real_thick
    real(kind = x_precision),dimension(nb_it),intent(in)  :: sigma_real_thin
    real(kind = x_precision),dimension(nb_it),intent(in)  :: temp_real_thick

    real(kind = x_precision),dimension(nb_it),intent(out) :: sigma_real
    real(kind = x_precision),dimension(nb_it),intent(out) :: temp_real
    !-----------------------------------------------------------------------

    integer                                               :: i

    do i = 1,index_scp
      sigma_real(i) = sigma_real_thick(i)
      temp_real(i) = temp_real_thick(i)
    enddo
    do i = index_scp + 1, nb_it
      sigma_real(i) = sigma_real_thin(i)
      temp_real(i) = temp_real_thick(i)
    enddo

  end subroutine build_s_curve


  !-------------------------------------------------------------------------
  ! Subroutine in order to display the critical points
  !-------------------------------------------------------------------------
  subroutine display_critical_points(sigma_c_thin, temp_c_thin,sigma_c_thick, temp_c_thick,k)
    implicit none

    real(kind = x_precision),intent(in) :: sigma_c_thin
    real(kind = x_precision),intent(in) :: temp_c_thin
    real(kind = x_precision),intent(in) :: sigma_c_thick
    real(kind = x_precision),intent(in) :: temp_c_thick
    integer                 ,intent(in) :: k
    !------------------------------------------------------------------------
    write(*,*)'**** Critical Point',k,'********'
    write(*,*)'****************************************'
    write(*,"(' Optically thin (T,sigma) :',1p,E12.4,4x,1p,E12.4)")temp_c_thin,sigma_c_thin
    write(*,"(' Optically thick (T,sigma):',1p,E12.4,4x,1p,E12.4)")temp_c_thick,sigma_c_thick

    write(*,*)'****************************************'

  end subroutine display_critical_points


  !-------------------------------------------------------------------------
  ! Subroutine to save in a file the critical points
  !-------------------------------------------------------------------------
  subroutine write_critical_points(n,radius, sigma_c_thin, temp_c_thin,sigma_c_thick, temp_c_thick,&
                                             tau_thin, tau_thick)

    implicit none

    integer                                  ,intent(in) :: n
    real(kind = x_precision),dimension(n)    ,intent(in) :: radius
    real(kind = x_precision),dimension(n)    ,intent(in) :: sigma_c_thin
    real(kind = x_precision),dimension(n)    ,intent(in) :: temp_c_thin
    real(kind = x_precision),dimension(n)    ,intent(in) :: sigma_c_thick
    real(kind = x_precision),dimension(n)    ,intent(in) :: temp_c_thick
    real(kind = x_precision),dimension(n)    ,intent(in) :: tau_thin
    real(kind = x_precision),dimension(n)    ,intent(in) :: tau_thick
    !-----------------------------------------------------------------------
    integer                                             :: i
    integer                                             :: fid_4 = 11
    character(len = 64)                                 :: fname_4

    !-----------------------------------------------------------------------

    fname_4 = 'critical_points/file.dat'

    open(fid_4,file  = fname_4, status='unknown',action='readwrite')
    write(fid_4,'(7(A16))') 'Radius','Temp_thin','Sigma_thin',&
                                     'Temp_thick','Sigma_thick','Tau_thin','Tau_thick'

      do i=1, n
        write(fid_4,'(7(e16.6e2))') radius(i), temp_c_thin(i),sigma_c_thin(i),&
                                                 temp_c_thick(i),sigma_c_thick(i),tau_thin(i),tau_thick(i)
      end do
    close(fid_4)

  end subroutine write_critical_points


  !-------------------------------------------------------------------------
  !Subroutine in order to compute variables H, rho, cs, nu, Q_plus, Q_minus,
  !K_ff, K_e, tau_eff, E_ff,Fz given T, Sigma and Omega
  !------------------------------------------------------------------------
  subroutine variables(T, Sigma, Omega, f, optical_depth, tau_eff)
    implicit none

    real(kind = x_precision),intent(in)  :: T,Sigma,Omega
    integer, intent(in)                  :: optical_depth
    real(kind = x_precision),intent(out) :: f,tau_eff
    real(kind = x_precision)             :: coeff_a,coeff_b,coeff_c
    real(kind = x_precision)             :: delta

    real(kind = x_precision)             :: H
    real(kind = x_precision)             :: rho
    real(kind = x_precision)             :: cs
    real(kind = x_precision)             :: nu
    real(kind = x_precision)             :: Q_plus
    real(kind = x_precision)             :: K_ff
    real(kind = x_precision)             :: K_e
    real(kind = x_precision)             :: E_ff
    real(kind = x_precision)             :: Fz
    real(kind = x_precision)             :: Q_minus
    !------------------------------------------------------------------------

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
  real(kind = x_precision) function dichotomy(S_min, S_max, max_it, eps, T, omega, optical_depth)
    implicit none

    real(kind = x_precision),intent(inout) :: S_min,S_max
    real(kind = x_precision),intent(in)    :: eps
    real(kind = x_precision),intent(in)    :: T
    real(kind = x_precision),intent(in)    :: omega
    integer,intent(in)                     :: max_it
    integer,intent(in)                     :: optical_depth

    integer                                :: j
    real(kind = x_precision)               :: tau_eff
    real(kind = x_precision)               :: Smin
    real(kind = x_precision)               :: Smax
    real(kind = x_precision)               :: f_min
    real(kind = x_precision)               :: f_max
    real(kind = x_precision)               :: f_center
    real(kind = x_precision)               :: S_center

    !-------------------------------------------------------------------------
    ! N-> Number of iterations for the dichotomy
    ! Smin, Smax -> Starting range
    ! eps -> Precision
    ! T-> Fixed variable
    !-------------------------------------------------------------------------

    Smin                 = S_min
    Smax                 = S_max

    S_center             = (Smin + Smax) / 2._x_precision

    j = 0

    call variables(T, Smin, Omega, f_min, optical_depth, tau_eff)

    call variables(T, Smax, Omega, f_max, optical_depth, tau_eff)

    !write(*,*)'fmin = ',f_min
    !write(*,*)'fmax = ',f_max

    if ( f_max * f_min > 0.) then
      !dichotomy = 0

    else if( f_max * f_min < 0.) then
      iteration:do while (dabs(Smax - Smin) >= eps .and. j < max_it)

        call variables(T, Smin, Omega, f_min, optical_depth, tau_eff)

        call variables(T, Smax, Omega, f_max, optical_depth, tau_eff)

        call variables(T, S_center, Omega, f_center, optical_depth, tau_eff)

        if (f_min * f_center > 0.) then

          Smin = S_center

        else if (f_max * f_center > 0.) then

          Smax = S_center

        endif

        S_center = (Smin + Smax) / 2._x_precision
        j = j + 1

      end do iteration

      dichotomy = S_center

    endif

  end function dichotomy

end module mod_s_curve
