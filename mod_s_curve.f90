module mod_S_curve
  use mod_constants
  use mod_read_parameters
  use mod_variables
  implicit none

contains
  !-------------------------------------------------------------------------
  ! Subroutine in order to find the result of the equation Q+ - Q- = 0
  ! output : datafile with the temperature T and the surface density Sigma
  ! Range of surface density [Smin, Smax] and temperature [T_min, T_max]
  ! values should be changed in order to find more coherent results
  !-------------------------------------------------------------------------
  subroutine curve()
    implicit none
    integer                                                :: i     = 0
    integer                                                :: j     = 0
    integer                                                :: k     = 0
    integer                                                :: l     = 0

    real(kind = x_precision)                               :: temp  = 0.0d0
    real(kind = x_precision), parameter                    :: t_min = 4.0d-1
    real(kind = x_precision), parameter                    :: t_max = 3.49d0

    integer,                  parameter                    :: nb_it = 100
    real(kind = x_precision)                               :: eps   = 1.0d-4
    real(kind = x_precision)                               :: eps2   = 5.0d-1


    real(kind = x_precision)                               :: sigma = 0.0d0
    real(kind = x_precision)                               :: Smin  = 0.0d0
    real(kind = x_precision)                               :: Smax  = 0.0d0
    real(kind = x_precision)                               :: r     = 0.0d0
    real(kind = x_precision)                               :: omega = 0.0d0
    real(kind = x_precision)                               :: H     = 0.0d0
    real(kind = x_precision)                               :: rho   = 0.0d0
    real(kind = x_precision)                               :: cs    = 0.0d0
    real(kind = x_precision)                               :: nu    = 0.0d0
    real(kind = x_precision)                               :: Q_plus = 0.0d0
    real(kind = x_precision)                               :: Q_minus= 0.0d0
    real(kind = x_precision)                               :: K_ff   = 0.0d0
    real(kind = x_precision)                               :: K_e    = 0.0d0
    real(kind = x_precision)                               :: tau_eff= 0.0d0
    real(kind = x_precision)                               :: P_rad  = 0.0d0
    real(kind = x_precision)                               :: P_gaz  = 0.0d0
    real(kind = x_precision)                               :: f      = 0.0d0
    real(kind = x_precision)                               :: E_ff   = 0.0d0
    real(kind = x_precision)                               :: Fz     = 0.0d0
    integer                                                :: optical_depth =0

    real(kind = x_precision),dimension(nb_it)              :: temp_real_1 =0.0d0
    real(kind = x_precision),dimension(nb_it)              :: sigma_real_1=0.0d0
    real(kind = x_precision),dimension(nb_it)              :: temp_real_0 =0.0d0
    real(kind = x_precision),dimension(nb_it)              :: sigma_real_0=0.0d0
    real(kind = x_precision),dimension(nb_it)              :: temp_real =0.0d0
    real(kind = x_precision),dimension(nb_it)              :: sigma_real=0.0d0

    character(len = 8)                                     :: number_of_cell
    character(len = 64)                                    :: fname
    character(len = 64)                                    :: fname_2
    character(len = 64)                                    :: fname_3

    integer                                                :: fid
    integer                                                :: fid_2
    integer                                                :: fid_3

    integer                                                :: n_cell = 1


    integer                                                :: index_fcp
    real(kind = x_precision)                               :: sigma_c_thick
    real(kind = x_precision)                               :: temp_c_thick
    integer                                                :: index_scp
    real(kind = x_precision)                               :: sigma_c_thin
    real(kind = x_precision)                               :: temp_c_thin
    
    !------------------------------------------------------------------------

    !-------------------------------------------------------------------------
    ! Test for 1 value of r
    !-------------------------------------------------------------------------
     do k              = 1 , n_cell
         !omega          = x_state%Omega(k)
         r              = 10._x_precision*G*params%M/(c**2)
         omega          = sqrt(G*params%M/r**3) / state_0%Omega_0

         write(number_of_cell,'(I5.5)') k
         fid = 20 + k
         fid_2 = 21 + k
         fid_3 = 22 + k

         fname = 's_curves/Temperature_Sigma_'//trim(number_of_cell)//'_1.dat'
         fname_2 = 's_curves/Temperature_Sigma_'//trim(number_of_cell)//'_0.dat'
         fname_3 = 's_curves/Temperature_Sigma_'//trim(number_of_cell)//'_tot.dat'

         open(fid,file  = fname, status='unknown',action='readwrite')
         open(fid_2,file  = fname_2, status='unknown',action='readwrite')
         open(fid_3,file  = fname_3, status='unknown',action='readwrite')

         do i          = 1, nb_it
          Smin        = 1d-2
          Smax        = 1d54

          temp        = (t_max-t_min)/(nb_it-1)*(i-1) + t_min

          optical_depth = 1

          sigma       = dichotomy(Smin, Smax, eps, temp, omega, optical_depth)
          
        !  call variables(temp, sigma, omega, H, rho, cs, nu, Q_plus, Q_minus, K_ff,&
        !      K_e, tau_eff, P_rad, P_gaz, E_ff, Fz, f, optical_depth)

        !  call display_variables(temp,Omega,r, sigma, H, rho, cs, nu, Q_plus, Q_minus,&
        !      K_ff, K_e, tau_eff, P_rad, P_gaz, E_ff, Fz, f)

          temp_real_1(i)   = log10(temp * state_0%T_0)
          sigma_real_1(i)  = log10(sigma * state_0%S_0)

          do j          = 1, nb_it
             write(fid_2,'(1p,E12.6,4x,1p,E12.6,4x,1p,E12.6)') sigma_real_1(j), temp_real_1(j)
          enddo

          optical_depth = 0

          sigma       = dichotomy(Smin, Smax, eps, temp, omega, optical_depth)
          
        !  call variables(temp, sigma, omega, H, rho, cs, nu, Q_plus, Q_minus, K_ff,&
        !      K_e, tau_eff, P_rad, P_gaz, E_ff, Fz, f, optical_depth)

        !  call display_variables(temp, Omega, r, sigma, H, rho, cs, nu, Q_plus, Q_minus,&
        !      K_ff, K_e, tau_eff, P_rad, P_gaz, E_ff, Fz, f)

          temp_real_0(i)   = log10(temp * state_0%T_0)
          sigma_real_0(i)  = log10(sigma * state_0%S_0)

          do l          = 1, nb_it
             write(fid,'(1p,E12.6,4x,1p,E12.6,4x,1p,E12.6)') sigma_real_0(l), temp_real_0(l)
          enddo
          
       enddo


       call first_critical_point(sigma_real_1, temp_real_1, index_fcp,sigma_c_thick, temp_c_thick, nb_it)

       call second_critical_point(sigma_real_1, sigma_real_0, temp_real_1,&
         index_fcp, index_scp, sigma_c_thin, temp_c_thin, nb_it, eps2)

       call build_s_curve(sigma_real_1, sigma_real_0, temp_real_1,nb_it, index_scp, sigma_real, temp_real)


          do l          = 1, nb_it
             write(fid_3,'(1p,E12.6,4x,1p,E12.6,4x,1p,E12.6)') sigma_real(l), temp_real(l)
          enddo
       
       
       close(fid)
       close(fid_2)
       close(fid_3)

    enddo

  end subroutine curve


  
  !-------------------------------------------------------------------------
  ! Subroutine in order to find the first critical point
  !-------------------------------------------------------------------------
  subroutine first_critical_point(sigma_real_1,temp_real_1, index_fcp,sigma_c_thick, temp_c_thick, nb_it)
    implicit none
    integer,intent(in)                                   :: nb_it

    real(kind = x_precision),dimension(nb_it),intent(in) :: sigma_real_1
    real(kind = x_precision),dimension(nb_it),intent(in) :: temp_real_1

    integer                                  ,intent(out):: index_fcp
    real(kind = x_precision)                 ,intent(out):: sigma_c_thick
    real(kind = x_precision)                 ,intent(out):: temp_c_thick
    integer                                              :: i
    i = 1
    do while (sigma_real_1(i) .le. sigma_real_1(i+1) .and. i .lt. nb_it - 1)
       index_fcp = i
       sigma_c_thick = sigma_real_1(i)
       temp_c_thick = temp_real_1(i)
       i = i + 1
    enddo

    endsubroutine first_critical_point





 !-------------------------------------------------------------------------
  ! Subroutine in order to find the second critical point
  !-------------------------------------------------------------------------
    subroutine second_critical_point(sigma_real_1, sigma_real_0, temp_real_1,&
         index_fcp, index_scp, sigma_c_thin, temp_c_thin, nb_it, eps2)
    implicit none
    integer,intent(in)                                   :: nb_it

    real(kind = x_precision),dimension(nb_it),intent(in) :: sigma_real_1
    real(kind = x_precision),dimension(nb_it),intent(in) :: temp_real_1
    real(kind = x_precision),dimension(nb_it),intent(in) :: sigma_real_0

    integer                                  ,intent(in) :: index_fcp
    integer                                  ,intent(out):: index_scp
    real(kind = x_precision)                 ,intent(in):: eps2

    real(kind = x_precision)                 ,intent(out):: sigma_c_thin
    real(kind = x_precision)                 ,intent(out):: temp_c_thin
    integer                                              :: i = 0
            ! write(*,*)index_fcp

    do i = index_fcp, nb_it
      ! write(*,*)sigma_real_0(i), sigma_real_1(i)
     !  write(*,*)(dabs(sigma_real_1(i) - sigma_real_0(i+1)))
    end do
    
        do i = index_fcp, nb_it - 1

           if(dabs(sigma_real_1(i) - sigma_real_0(i+1)) .lt. eps2)then
           index_scp = i
            sigma_c_thin = (sigma_real_1(i) + sigma_real_0(i+1))/2._x_precision
           
           temp_c_thin = temp_real_1(i)
           end if
        enddo

    endsubroutine second_critical_point


    
    subroutine build_s_curve(sigma_real_1, sigma_real_0, temp_real_1,nb_it, index_scp, sigma_real, temp_real)
      implicit none
      integer,intent(in)                                   :: nb_it
      
      real(kind = x_precision),dimension(nb_it),intent(in) :: sigma_real_1
      real(kind = x_precision),dimension(nb_it),intent(in) :: sigma_real_0
      real(kind = x_precision),dimension(nb_it),intent(in) :: temp_real_1
      integer                                  ,intent(in) :: index_scp
      !real(kind = x_precision),dimension(:),allocatable,intent(out) :: sigma_real
      !real(kind = x_precision),dimension(:),allocatable,intent(out) :: temp_real
      real(kind = x_precision),dimension(nb_it),intent(out) :: sigma_real
      real(kind = x_precision),dimension(nb_it),intent(out) :: temp_real
      integer::i = 0

      do i = 1,index_scp
         sigma_real(i) = sigma_real_1(i)
         temp_real(i) = temp_real_1(i)
      enddo
      do i = index_scp + 1, nb_it
         sigma_real(i) = sigma_real_0(i)
         temp_real(i) = temp_real_1(i)
      enddo
      
    end subroutine build_s_curve

    
  !-------------------------------------------------------------------------
  !Subroutine for resolving a quadratic equation as long as the solutions are
  !a set of real numbers
  !-------------------------------------------------------------------------
  subroutine quadratic(coeff_a, coeff_b, coeff_c, sol)

    implicit none
    real(kind=x_precision), intent(in)                     :: coeff_a,coeff_b,coeff_c
    real(kind=x_precision), intent(out)                    :: sol
    real(kind=x_precision)                                 :: sol_1 = 0.0d0
    real(kind=x_precision)                                 :: sol_2 = 0.0d0
    real(kind=x_precision)                                 :: delta = 0.0d0
    !------------------------------------------------------------------------

    delta          = coeff_b**2 - 4_x_precision * coeff_a * coeff_c

    if (delta .lt. 0.) then
      write(*,*)'No solutions in the R field.'
   else
      sol_1         = -0.5_x_precision *(coeff_b + sign(sqrt(delta),coeff_b))/coeff_a
      if (sol_1 .lt. 0.d-12) then
         write(*,*)'Problem'
         stop
      endif
      sol_2 = (coeff_c / (coeff_a * sol_1))
      sol = max(sol_1,sol_2)

    end if
  end subroutine quadratic


  !-------------------------------------------------------------------------
  !Subroutine in order to compute variables H, rho, cs, nu, Q_plus, Q_minus,
  !K_ff, K_e, tau_eff, P_rad, P_gaz,E_ff,Fz for T, Sigma and Omega given
  !------------------------------------------------------------------------
  subroutine variables(T, Sigma, Omega, H, rho, cs, nu, Q_plus, Q_minus,&
       K_ff, K_e, tau_eff, P_rad, P_gaz, E_ff, Fz, f, optical_depth)
    implicit none

    real(kind = x_precision),intent(in)                      :: T,Sigma,Omega
    integer, intent(in)                                      :: optical_depth
    real(kind = x_precision)                                 :: coeff_a=0.,coeff_b=0.,coeff_c=0.

    real(kind = x_precision),intent(out)                     :: H
    real(kind = x_precision),intent(out)                     :: rho
    real(kind = x_precision),intent(out)                     :: cs
    real(kind = x_precision),intent(out)                     :: nu
    real(kind = x_precision),intent(out)                     :: Q_plus
    real(kind = x_precision),intent(out)                     :: K_ff
    real(kind = x_precision),intent(out)                     :: K_e
    real(kind = x_precision),intent(out)                     :: E_ff
    real(kind = x_precision),intent(out)                     :: tau_eff
    real(kind = x_precision),intent(out)                     :: Fz
    real(kind = x_precision),intent(out)                     :: Q_minus
    real(kind = x_precision),intent(out)                     :: P_rad
    real(kind = x_precision),intent(out)                     :: P_gaz
    real(kind = x_precision),intent(out)                     :: f
    !------------------------------------------------------------------------

    coeff_a              = (Omega**2 * state_0%Omega_0**2 * Sigma * state_0%S_0)/2._x_precision
    coeff_b              = (-1._x_precision/3._x_precision) * cst_rad*T**4 * state_0%T_0**4 / state_0%H_0
    coeff_c              = (-1._x_precision * params%RTM * T  *  Sigma * state_0%S_0)/(2._x_precision * state_0%H_0**2)

    call quadratic(coeff_a , coeff_b , coeff_c , H)

    rho                  = Sigma / H
    P_rad                = T**4
    P_gaz                = rho * T
    cs                   = Omega * H
    nu                   = params%alpha * cs * H
    K_ff                 = 6.13d22 * state_0%rho_0 * rho * (state_0%T_0 * T)**(-3.5_x_precision)
    K_e                  = params%kappa_e
    E_ff                 = 6.22d20 * (state_0%rho_0 * rho)**2 * sqrt(state_0%T_0 * T)
    tau_eff              = 0.5_x_precision * sqrt(K_e * K_ff) * Sigma * state_0%S_0


    !-------------------------------------------------------------------------
    !Select the case for the optical depth to compute Fz
    !-------------------------------------------------------------------------

   ! if (tau_eff .ge. 1.)  then
   !    optical_depth     = 1
   ! else
   !    optical_depth     = 0
   ! end if

    select case(optical_depth)

    case(1)

       Fz = 2._x_precision * c**2 * T**4 /(27._x_precision * sqrt(3._x_precision) &
            * (K_ff + K_e) * Sigma * state_0%S_0)
    case (0)

       Fz = 4._x_precision * state_0%H_0 * E_ff * H / (state_0%Omega_0 * state_0%S_0)

    end select

    Q_plus              = 3._x_precision  * state_0%H_0**2 * nu * Omega**2 * state_0%Omega_0**2
    Q_minus             = Fz  / Sigma

    f                   = Q_plus - Q_minus

  end subroutine variables


  !-------------------------------------------------------------------------
  ! Dichotomic function in order to determine the change of sign in a given
  ! interval [Smin,Smax] with an epsilon precision
  !-------------------------------------------------------------------------
  real(kind=x_precision) function dichotomy(Smin, Smax, eps, T, omega, optical_depth)
    use mod_read_parameters
    use mod_constants
    use mod_variables
    implicit none

    real(kind=x_precision),intent(inout)                     :: Smin,Smax
    real(kind=x_precision),intent(in)                        :: eps
    real(kind=x_precision),intent(in)                        :: T
    real(kind=x_precision),intent(in)                        :: omega
    integer,intent(in)                                       :: optical_depth

    real(kind=x_precision)                                   :: H
    real(kind=x_precision)                                   :: rho
    real(kind=x_precision)                                   :: cs
    real(kind=x_precision)                                   :: nu
    real(kind=x_precision)                                   :: Q_plus
    real(kind=x_precision)                                   :: Q_minus
    real(kind=x_precision)                                   :: K_ff
    real(kind=x_precision)                                   :: K_e
    real(kind=x_precision)                                   :: tau_eff
    real(kind=x_precision)                                   :: P_rad
    real(kind=x_precision)                                   :: P_gaz
    real(kind=x_precision)                                   :: f_min
    real(kind=x_precision)                                   :: f_max
    real(kind=x_precision)                                   :: f_center
    real(kind=x_precision)                                   :: E_ff
    real(kind=x_precision)                                   :: Fz
    real(kind=x_precision)                                   :: S_center = 0._x_precision
    integer                                                  :: j = 0


    !-------------------------------------------------------------------------
    ! N-> Number of iterations for the dichotomy
    ! Smin, Smax -> Starting range
    ! eps -> Precision
    ! T-> Fixed variable
    !-------------------------------------------------------------------------
    dichotomy             = (Smin+Smax)/2.
    j = 0
    call variables(T, Smin, Omega, H, rho, cs, nu, Q_plus, Q_minus, K_ff, K_e,&
         tau_eff, P_rad, P_gaz,E_ff,Fz,f_min,optical_depth)

    call variables(T, Smax, Omega, H, rho, cs, nu, Q_plus, Q_minus, K_ff, K_e, &
         tau_eff, P_rad, P_gaz,E_ff,Fz,f_max, optical_depth)

  !   write(*,*)'fmin = ',f_min
  !   write(*,*)'fmax = ',f_max

    if ( f_max * f_min .gt. 0.) then
   !    write(*,*)'This function image does not switch its sign in this particular interval.'
       dichotomy = 0

    else if( f_max * f_min .lt. 0.) then
       iteration:do while ( dabs( Smax - Smin ) .ge. eps .and. j .lt. 10000)


    call variables(T, Smin, Omega, H, rho, cs, nu, Q_plus, Q_minus, K_ff, K_e,&
         tau_eff, P_rad, P_gaz,E_ff,Fz,f_min, optical_depth)

    call variables(T, Smax, Omega, H, rho, cs, nu, Q_plus, Q_minus, K_ff, K_e, &
         tau_eff, P_rad, P_gaz,E_ff,Fz,f_max, optical_depth)

    call variables(T, S_center, Omega, H, rho, cs, nu, Q_plus, Q_minus, K_ff,&
         K_e, tau_eff, P_rad, P_gaz,E_ff,Fz,f_center, optical_depth)


          if(f_min * f_center .gt. 0.) then

             Smin         = S_center

          else if (f_max * f_center .gt. 0.) then
             Smax         = S_center

          endif

          S_center       = (Smin + Smax) * 1._x_precision / 2._x_precision
          j               = j + 1

       end do iteration
         dichotomy = S_center

    endif


  end function dichotomy


  !Subroutine in order to display variables
  !-------------------------------------------------------------------------
  subroutine display_variables(temp,Omega,r,sigma, H, rho, cs, nu, Q_plus, &
       Q_minus, K_ff, K_e, tau_eff, P_rad, P_gaz,E_ff,Fz,f)
    implicit none

    real(kind = x_precision), intent(in)                     :: temp
    real(kind = x_precision), intent(in)                     :: Omega
    real(kind = x_precision), intent(in)                     :: r
    real(kind = x_precision), intent(in)                     :: sigma
    real(kind = x_precision), intent(in)                     :: H
    real(kind = x_precision), intent(in)                     :: rho
    real(kind = x_precision), intent(in)                     :: cs
    real(kind = x_precision), intent(in)                     :: nu
    real(kind = x_precision), intent(in)                     :: Q_plus
    real(kind = x_precision), intent(in)                     :: Q_minus
    real(kind = x_precision), intent(in)                     :: K_ff
    real(kind = x_precision), intent(in)                     :: K_e
    real(kind = x_precision), intent(in)                     :: tau_eff
    real(kind = x_precision), intent(in)                     :: P_rad
    real(kind = x_precision), intent(in)                     :: P_gaz
    real(kind = x_precision), intent(in)                     :: f
    real(kind = x_precision), intent(in)                     :: E_ff
    real(kind = x_precision), intent(in)                     :: Fz
    !------------------------------------------------------------------------

    write(*,*)'               Variables                '
    write(*,*)'****************************************'
    write(*,"(' r           =',1p,E12.4)") r
    write(*,"(' Temp        =',1p,E12.4)") temp
    write(*,"(' Sigma       =',1p,E12.4)") sigma
    write(*,"(' Omega       =',1p,E12.4)") Omega
    write(*,"(' H           =',1p,E12.4)") H
    write(*,"(' rho         =',1p,E12.4)") rho
    write(*,"(' cs          =',1p,E12.4)") cs
    write(*,"(' nu          =',1p,E12.4)") nu
    write(*,"(' Q plus      =',1p,E12.4)") Q_plus
    write(*,"(' Q minus     =',1p,E12.4)") Q_minus
    write(*,"(' K_ff        =',1p,E12.4)") K_ff
    write(*,"(' K_e         =',1p,E12.4)") K_e
    write(*,"(' tau_eff     =',1p,E12.4)") tau_eff
    write(*,"(' P_gaz       =',1p,E12.4)") P_gaz
    write(*,"(' P_rad       =',1p,E12.4)") P_rad
    write(*,"(' E_ff        =',1p,E12.4)") E_ff
    write(*,"(' delta Q     =',1p,E12.4)") f
    write(*,"(' Fz          =',1p,E12.4)") Fz
    write(*,*)'****************************************'
  end subroutine display_variables


end module mod_S_curve
