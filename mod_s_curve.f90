module mod_s_curve
  use mod_constants
  use mod_read_parameters
  use mod_variables
  implicit none

! _____________________________________________________________________________________________

  ! BE SURE TO CREATE TWO DIRECTORIES IN YOUR WORKING DIRECTORY : s_curves ( & critical_points )

  !PROBLEMATIC FOR THE SECOND CRITICAL POINT.

!_____________________________________________________________________________________________

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

  !               intial_variables      : : Computing some variables
  !                                         --> output : rs, rmin, Mdot_0, Sigma_0, Omega_0, T_0 (+ nu_0 )

  !               variables             : : Computing all the variables for Q+ and Q-
  !                                     <-- input : Temperature, Sigma, Omega + optical depth indicator (0 or 1) + more
  !                                         --> output : ...


  !               dichotomy             : : In order to find the result of the equation Q+ - Q- = 0
  !                                     <-- input : Temperature, Sigma, Omega + optical depth indicator (0 or 1) + precision + more
  !                                         --> output : Sigma between a computed range of surface density [Smin, Smax] for a given temperature

  !
  !-------------------------------------------------------------------------

  subroutine curve(nb_it, eps, t_min, t_max, dt, s_min, s_max, temperature, s)
    implicit none

    integer                                   ,intent(in)  :: nb_it         
    real(kind = x_precision)                  ,intent(in)  :: eps          
    real(kind = x_precision)                  ,intent(in)  :: t_min        
    real(kind = x_precision)                  ,intent(in)  :: t_max         
    real(kind = x_precision)                  ,intent(in)  :: dt             
    real(kind = x_precision)                  ,intent(in)  :: S_min
    real(kind = x_precision)                  ,intent(in)  :: S_max

    real(kind = x_precision),dimension(n_cell),intent(out) :: temperature
    real(kind = x_precision),dimension(n_cell),intent(out) :: s
 !  __________________________________________________________________________  
    integer                                                :: i     = 0
    integer                                                :: j     = 0
    integer                                                :: k     = 0
    integer                                                :: file_length = 0

 ! For the temperature and s 
    real(kind = x_precision)                               :: temp  = 0.0d0
    real(kind = x_precision)                               :: Smin  = 0.0d0
    real(kind = x_precision)                               :: Smax  = 0.0d0
    real(kind = x_precision)                               :: Sigma = 0.0d0

    real(kind = x_precision)                               :: omega = 0.0d0
    real(kind = x_precision)                               :: r     = 0.0d0
    real(kind = x_precision)                               :: f     = 0.0d0
    integer                                                :: optical_depth = 0      ! Indicator for the optical thickness

 ! For the separate curves given the thickness
    real(kind = x_precision),dimension(nb_it)              :: temp_t_thick
    real(kind = x_precision),dimension(nb_it)              :: sigma_t_thick
    real(kind = x_precision),dimension(nb_it)              :: temp_t_thin
    real(kind = x_precision),dimension(nb_it)              :: sigma_t_thin
    real(kind = x_precision),dimension(nb_it)              :: temp_real             ! Temperature for the S curve
    real(kind = x_precision),dimension(nb_it)              :: sigma_real            ! Surface density for the S curve

 ! For the output file 
    character(len = 8)                                     :: number_of_cell
    character(len = 64)                                    :: fname_tot
    integer                                                :: fid_tot

 ! For the critical points 
    integer                                                :: index_fcp
    real(kind = x_precision)                               :: sigma_c_thick
    real(kind = x_precision)                               :: temp_c_thick
    real(kind = x_precision),dimension(n_cell)             :: sigma_thick
    real(kind = x_precision),dimension(n_cell)             :: temp_thick
    integer                                                :: index_scp
    real(kind = x_precision)                               :: sigma_c_thin
    real(kind = x_precision)                               :: temp_c_thin
    real(kind = x_precision),dimension(n_cell)             :: sigma_thin
    real(kind = x_precision),dimension(n_cell)             :: temp_thin

 !____________________________________________________________________________________

    do k = 1, n_cell
       file_length = 0

       r     = r_state%r(k)
       omega = x_state%Omega(k)

       write(number_of_cell,'(I5.5)') k
       fid_tot = 22 + k
       fname_tot = 's_curves/Temperature_Sigma_'//trim(number_of_cell)//'_tot.dat'

       do i = 1, nb_it
          temp          = dt * (i-1) + t_min

          ! Optical thick case
          optical_depth = 1

          Smin          = S_min
          Smax          = S_max

          sigma         = dichotomy(Smin, Smax, eps, temp, omega, optical_depth)

          call variables(temp, sigma, omega, f, optical_depth)

          temp_t_thick(i)  = temp
          sigma_t_thick(i) = sigma

          ! Optical thin case
          optical_depth = 0

          Smin          = S_min
          Smax          = S_max

          sigma         = dichotomy(Smin, Smax, eps, temp, omega, optical_depth)

          call variables(temp, sigma, omega, f, optical_depth)

          temp_t_thin(i)   =  temp
          sigma_t_thin(i)  =  sigma
 
        enddo


      call first_critical_point(sigma_t_thick, temp_t_thick, index_fcp,sigma_c_thick, temp_c_thick, nb_it)

      call second_critical_point(sigma_t_thick, sigma_t_thin, temp_t_thick,&
        index_fcp, index_scp, sigma_c_thin, temp_c_thin, nb_it)

      call build_s_curve(sigma_t_thick, sigma_t_thin, temp_t_thick, nb_it, index_scp, sigma_real, temp_real)

      call display_critical_points(sigma_c_thin, temp_c_thin,sigma_c_thick, temp_c_thick, k)


 ! Saving S and T for the S curve in a file
      do j = 1, nb_it
        temp_real(j)  =  temp_real(j) * state_0%T_0 
        sigma_real(j) =  sigma_real(j) * state_0%S_0 
        call save_data(fid_tot, fname_tot, nb_it, file_length, sigma_real(j), temp_real(j))
        file_length = file_length + 1
      enddo

      close(fid_tot)

 ! Critical T and S 
      temperature(k) = temp_c_thick * state_0%T_0
      s(k)           = sigma_c_thick * state_0%S_0

      temp_thick(k)  =  temp_c_thick * state_0%T_0 
      sigma_thick(k) =  sigma_c_thick * state_0%S_0 
      temp_thin(k)   =  temp_c_thin * state_0%T_0 
      sigma_thin(k)  =  sigma_c_thin * state_0%S_0 

    enddo

    call write_critical_points(n_cell, r_state%r, sigma_thin, temp_thin, sigma_thick, temp_thick)

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
 !____________________________________________________________________________________

    integer                                               :: i

    i = 1
    do while (sigma_real_thick(i) < sigma_real_thick(i+1) .and. i < nb_it - 1)
      index_fcp     = i
      sigma_c_thick = sigma_real_thick(i)
      i             = i + 1
    enddo
    temp_c_thick  = temp_real_thick(i)

  endsubroutine first_critical_point


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
 !____________________________________________________________________________________
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


  subroutine build_s_curve(sigma_real_thick, sigma_real_thin, temp_real_thick, nb_it, index_scp, sigma_real, temp_real)
    implicit none

    integer                                  ,intent(in)  :: nb_it
    integer                                  ,intent(in)  :: index_scp
    real(kind = x_precision),dimension(nb_it),intent(in)  :: sigma_real_thick
    real(kind = x_precision),dimension(nb_it),intent(in)  :: sigma_real_thin
    real(kind = x_precision),dimension(nb_it),intent(in)  :: temp_real_thick

    real(kind = x_precision),dimension(nb_it),intent(out) :: sigma_real
    real(kind = x_precision),dimension(nb_it),intent(out) :: temp_real
 !____________________________________________________________________________________
    integer::i = 0

    do i = 1,index_scp
      sigma_real(i) = sigma_real_thick(i)
      temp_real(i) = temp_real_thick(i)
    enddo
    do i = index_scp + 1, nb_it
      sigma_real(i) = sigma_real_thin(i)
      temp_real(i) = temp_real_thick(i)
    enddo

  end subroutine build_s_curve


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


  subroutine write_critical_points(n,radius, sigma_c_thin, temp_c_thin,sigma_c_thick, temp_c_thick)

    implicit none

    integer                                  ,intent(in):: n
    real(kind = x_precision),dimension(n)    ,intent(in):: radius
    real(kind = x_precision),dimension(n)    ,intent(in):: sigma_c_thin
    real(kind = x_precision),dimension(n)    ,intent(in):: temp_c_thin
    real(kind = x_precision),dimension(n)    ,intent(in):: sigma_c_thick
    real(kind = x_precision),dimension(n)    ,intent(in):: temp_c_thick
 !____________________________________________________________________________________
    integer                                             :: i
    integer                                             :: fid_4 = 11
    character(len = 64)                                 :: fname_4

    !------------------------------------------------------------------------

     fname_4 = 'critical_points/file.dat'

     open(fid_4,file  = fname_4, status='unknown',action='readwrite')
     write(fid_4,'(1p,A6,10x,1p,A9,8x,1p,A9,8x,1p,A10,6x,1p,A11)')'Radius',&
                   'Temp_thin','Sigma_thin','Temp_thick','Sigma_thick'

        do i=1, n
           write(fid_4,'(1p,E12.6,4x,1p,E12.6,4x,1p,E12.6,4x,1p,E12.6,4x,1p,E12.6)')radius(i),&
                      temp_c_thin(i),sigma_c_thin(i),temp_c_thick(i),sigma_c_thick(i)
        end do
     close(fid_4)

  end subroutine write_critical_points


   subroutine save_data(fid, filename, n, length, data_1, data_2)
    implicit none

    integer                                  ,intent(in):: fid
    character(len =64)                       ,intent(in):: filename
    integer                                  ,intent(in):: n
    real(kind = x_precision)                 ,intent(in):: data_1
    real(kind = x_precision)                 ,intent(in):: data_2

    integer                               ,intent(inout):: length 

    !------------------------------------------------------------------------

    if( length .eq. 0 ) then
       open(fid,file  = filename, status='unknown', position='append', action='readwrite')
    end if
    
    write(fid,'(1p,E12.6,4x,1p,E12.6)') data_1, data_2

    if( length .eq. n ) then
       close(fid)     
    end if

  end subroutine save_data

  !-------------------------------------------------------------------------
  !Subroutine for resolving a quadratic equation as long as the solutions are a set of real numbers
  !-------------------------------------------------------------------------
  subroutine quadratic(coeff_a, coeff_b, coeff_c, sol)
    implicit none

    real(kind=x_precision), intent(in)  :: coeff_a,coeff_b,coeff_c
    real(kind=x_precision), intent(out) :: sol
    real(kind=x_precision)              :: sol_1 = 0.0d0
    real(kind=x_precision)              :: sol_2 = 0.0d0
    real(kind=x_precision)              :: delta = 0.0d0
    !------------------------------------------------------------------------

    delta          = coeff_b**2 - 4_x_precision * coeff_a * coeff_c

    if (coeff_a == 0.0) then
      write(*,*)'Coefficient a in the quadratic equation is nought.'
      continue
    else

      if (delta < 0.) then
        write(*,*)'No solutions in the R field.'
      else
        sol_1         = -0.5_x_precision * (coeff_b + sign(sqrt(delta),coeff_b)) / coeff_a
        if (sol_1 < 0.d-12) then
          write(*,*)'Problem'
          stop
        end if
        sol_2 = (coeff_c / (coeff_a * sol_1))
        sol = max(sol_1,sol_2)
        if (sol /= sol_1) then
          write(*,*)'The solution was the other one !'
        end if
      end if

    end if
  end subroutine quadratic


  !-------------------------------------------------------------------------
  !Subroutine in order to compute variables H, rho, cs, nu, Q_plus, Q_minus,
  !K_ff, K_e, tau_eff, E_ff,Fz given T, Sigma and Omega
  !------------------------------------------------------------------------
  subroutine variables(T, Sigma, Omega, f, optical_depth)
    implicit none

    real(kind = x_precision),intent(in)  :: T,Sigma,Omega
    integer, intent(in)                  :: optical_depth
    real(kind = x_precision),intent(out) :: f
    real(kind = x_precision)             :: coeff_a=0.,coeff_b=0.,coeff_c=0.

    real(kind = x_precision)             :: H
    real(kind = x_precision)             :: rho
    real(kind = x_precision)             :: cs
    real(kind = x_precision)             :: nu
    real(kind = x_precision)             :: Q_plus
    real(kind = x_precision)             :: K_ff
    real(kind = x_precision)             :: K_e
    real(kind = x_precision)             :: E_ff
    real(kind = x_precision)             :: tau_eff
    real(kind = x_precision)             :: Fz
    real(kind = x_precision)             :: Q_minus
    !------------------------------------------------------------------------

    coeff_a = (Omega * state_0%Omega_0)**2 * Sigma * state_0%S_0 / 2._x_precision
    coeff_b = (-1._x_precision/3._x_precision) * cst_rad * (T * state_0%T_0)**4 / state_0%H_0
    coeff_c = - params%RTM * T * Sigma * state_0%S_0 / (2._x_precision * state_0%H_0**2)

    call quadratic(coeff_a , coeff_b , coeff_c , H)

    rho     = Sigma / H
    cs      = Omega * H
    nu      = params%alpha * cs * H
    K_ff    = 6.13d22 * state_0%rho_0 * rho * (state_0%T_0 * T)**(-3.5_x_precision)
    K_e     = params%kappa_e
    E_ff    = 6.22d20 * (state_0%rho_0 * rho)**2 * sqrt(state_0%T_0 * T)
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
  real(kind=x_precision) function dichotomy(S_min, S_max, eps, T, omega, optical_depth)
    implicit none

    real(kind=x_precision),intent(inout) :: S_min,S_max
    real(kind=x_precision),intent(in)    :: eps
    real(kind=x_precision),intent(in)    :: T
    real(kind=x_precision),intent(in)    :: omega
    integer,intent(in)                   :: optical_depth

    real(kind=x_precision)               :: Smin
    real(kind=x_precision)               :: Smax
    real(kind=x_precision)               :: f_min
    real(kind=x_precision)               :: f_max
    real(kind=x_precision)               :: f_center
    real(kind=x_precision)               :: S_center = 0._x_precision
    integer                              :: j = 0

    !-------------------------------------------------------------------------
    ! N-> Number of iterations for the dichotomy
    ! Smin, Smax -> Starting range
    ! eps -> Precision
    ! T-> Fixed variable
    !-------------------------------------------------------------------------

    Smin                 = S_min
    Smax                 = S_max

    S_center             = (Smin+Smax)/2.
    j = 0
    call variables(T, Smin, Omega, f_min, optical_depth)

    call variables(T, Smax, Omega, f_max, optical_depth)

    !write(*,*)'fmin = ',f_min
    !write(*,*)'fmax = ',f_max

    if ( f_max * f_min > 0.) then
      !write(*,*)'This function image does not switch its sign in this particular interval.'
      !dichotomy = 0

    else if( f_max * f_min < 0.) then
      iteration:do while ( dabs( Smax - Smin ) >= eps .and. j < 10000)

        call variables(T, Smin, Omega, f_min, optical_depth)

        call variables(T, Smax, Omega, f_max, optical_depth)

        call variables(T, S_center, Omega, f_center, optical_depth)

        if (f_min * f_center > 0.) then

          Smin = S_center

        else if (f_max * f_center > 0.) then

          Smax = S_center

        endif

        S_center = (Smin + Smax) * 1._x_precision / 2._x_precision
        j = j + 1

      end do iteration

      dichotomy = S_center

    endif

  end function dichotomy

end module mod_s_curve
