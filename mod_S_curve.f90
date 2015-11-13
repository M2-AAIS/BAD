module mod_S_curve
  use mod_read_parameters
  use mod_constants
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
    real(kind = x_precision)                               :: temp  = 0.0d0
    real(kind = x_precision), parameter                    :: t_min = 5.715d-4
    real(kind = x_precision), parameter                    :: t_max = 1.815d-6
    integer,                  parameter                    :: nb_it = 200
    real(kind = x_precision)                               :: sigma = 0.0d0
    real(kind = x_precision)                               :: Smin  = 0.0d0
    real(kind = x_precision)                               :: Smax  = 0.0d0
    real(kind = x_precision)                               :: eps   = 1.0d-10
    real(kind = x_precision)                               :: omega = 0.0d0
    real(kind = x_precision)                               :: r     = 0.0d0
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
    real(kind = x_precision)                               :: rs
    real(kind = x_precision)                               :: rmin
    real(kind = x_precision)                               :: Mdot_0
    real(kind = x_precision)                               :: Sigma_0
    real(kind = x_precision)                               :: Omega_0
    real(kind = x_precision)                               :: rho_0
    real(kind = x_precision)                               :: T_0
    real(kind = x_precision)                               :: temp_reel=0.0d0
    real(kind = x_precision)                               :: sigma_reel=0.0d0
    character(len = 8)                                     :: number_of_cell
    character(len = 64)                                    :: fname
    integer                                                :: fid
    type(parameters)                                       :: param
    !------------------------------------------------------------------------
    call get_parameters(param)
    call display_parameters()
    call initial_variables(rs, rmin, Mdot_0, Sigma_0, Omega_0, T_0, rho_0)
    call display_initial_variables(rs, rmin, Mdot_0, Sigma_0, Omega_0, T_0, rho_0)
    !-------------------------------------------------------------------------
    ! Test for 1 value of r
    !-------------------------------------------------------------------------
    do j              = 1 , 1
       !do j              = 1 , n_cell
       !  r = (param%rmax-rmin)/(n_cell-1)*(j-1) + rmin
       r              = 10._x_precision*G*param%M/(c**2)     
       omega          = sqrt(G*param%M/r**3) / Omega_0 
       write(number_of_cell,'(I5.5)') j
       fid = 20 + j
       fname = 'Temperature_Sigma_'//trim(number_of_cell)//'.dat'
       open(fid,file  = fname, status='unknown',action='readwrite') 
       !  write(fid,*)"T     &         Sigma"
        do i          = 1 , nb_it
          Smin        = 5.0d-13
          Smax        = 3.0d13
          
          temp        = (t_max-t_min)/(nb_it-1)*(i-1) + t_min
          sigma       = dichotomy(Smin, Smax, eps, temp, omega, sigma_0, Omega_0,rs, T_0, rho_0)
          
          call variables(temp, sigma, omega, H, rho, cs, nu, Q_plus, Q_minus, K_ff,&
              K_e, tau_eff, P_rad, P_gaz,E_ff,Fz,f,Sigma_0, Omega_0,rs, T_0, rho_0)

          
          call display_variables(temp,Omega,r, sigma, H, rho, cs, nu, Q_plus, Q_minus,&
               K_ff, K_e, tau_eff, P_rad, P_gaz,E_ff,Fz,f)
          
          temp_reel   = temp * T_0
          sigma_reel  = sigma * sigma_0
          write(fid,'(1p,E12.6,4x,1p,E12.6)')temp_reel,sigma_reel
       enddo
       close(fid)
    enddo
    
  end subroutine curve


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
      sol_1=(-1.*coeff_b+sqrt(delta))/(2.*coeff_a)
      sol_2=(-1.*coeff_b-sqrt(delta))/(2.*coeff_a)
   end if
   if (sol_1 .gt. 0. ) then
      sol=sol_1
   else
      sol=sol_2
   end if
      !sol          = -0.5_x_precision *(coeff_b + sign(sqrt(delta),coeff_b))/coeff_a
    !end if
  end subroutine quadratic


  
  !-------------------------------------------------------------------------
  !Subroutine in order to compute initial variables rs, rmin, Mdot_0,
  !Sigma_0
  !-------------------------------------------------------------------------
  subroutine initial_variables(rs, rmin, Mdot_0, Sigma_0, Omega_0, T_0, rho_0)
    implicit none

    real(kind=x_precision),intent(out)                      :: rs
    real(kind=x_precision),intent(out)                      :: rmin
    real(kind=x_precision),intent(out)                      :: Mdot_0
    real(kind=x_precision),intent(out)                      :: Sigma_0
    real(kind=x_precision),intent(out)                      :: Omega_0
    real(kind=x_precision),intent(out)                      :: T_0
    real(kind=x_precision),intent(out)                      :: rho_0

    type(parameters)                                        :: param
    !------------------------------------------------------------------------

    call get_parameters(param)
    rs                 = 2._x_precision * G * param%M/(c**2)
    rmin               = 3._x_precision * rs
    Mdot_0             = param%Mdot
    Omega_0            = sqrt(G * param%M / rs**3 )
    Sigma_0            = Mdot_0 /(Omega_0 * rs**2 * 2_x_precision * pi)
    write(*,*)'Mdot_0',Mdot_0
    T_0                = (Mdot_0 * c**2 / (48._x_precision * pi * rs**2 * stefan * &
                         sqrt(27._x_precision) ) )**(1._x_precision/4._x_precision)
    rho_0              = Sigma_0 / (2._x_precision * rs)
    
  end subroutine initial_variables

  

  !-------------------------------------------------------------------------
  !Subroutine in order to compute variables H, rho, cs, nu, Q_plus, Q_minus,
  !K_ff, K_e, tau_eff, P_rad, P_gaz,E_ff,Fz for T, Sigma and Omega given
  !------------------------------------------------------------------------
  subroutine variables(T, Sigma, Omega, H, rho, cs, nu, Q_plus, Q_minus,&
       K_ff, K_e, tau_eff, P_rad, P_gaz,E_ff,Fz,f,Sigma_0, Omega_0,rs, T_0, rho_0)
    implicit none
 
    type(parameters)                                         :: param
    real(kind = x_precision),intent(in)                      :: T,Sigma,Omega
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
    real(kind = x_precision),intent(in)                      :: Sigma_0
    real(kind = x_precision),intent(in)                      :: Omega_0
    real(kind = x_precision),intent(in)                      :: T_0
    real(kind = x_precision),intent(in)                      :: rho_0
    real(kind = x_precision),intent(in)                      :: rs
    real(kind = x_precision)                                 :: nu_0

    real(kind = x_precision),intent(out)                     :: f 
    integer                                                  :: optical_depth = 0
    !------------------------------------------------------------------------

    call get_parameters(param)
   
    coeff_a              = (Omega**2 * Omega_0**2 * Sigma * Sigma_0)/2._x_precision
    coeff_b              = (-1._x_precision/3._x_precision) * cst_rad*T**4 * T_0**4 / rs
    coeff_c              = (-1._x_precision * param%RTM * T  *  Sigma * Sigma_0)/(2._x_precision * rs**2)

    call quadratic(coeff_a , coeff_b , coeff_c , H)

    rho                  = Sigma / H
    P_rad                = T**4
    P_gaz                = rho * T
    cs                   = Omega * H
    nu                   = param%alpha * cs * H
    K_ff                 = 6.13d22 * rho_0 * rho * (T_0 * T)**(-3.5_x_precision)
    K_e                  = 0.2_x_precision * (1._x_precision + param%X)
    E_ff                 = 6.22d20 * (rho_0 * rho)**2 * sqrt(T_0 * T)
    tau_eff              = 0.5_x_precision * sqrt(K_e * K_ff) * Sigma * Sigma_0
    nu_0                 = 2._x_precision * rs**2 * Omega_0/3._x_precision
    
    
    !-------------------------------------------------------------------------
    !Select the case for the opticaly depth to compute Fz
    !-------------------------------------------------------------------------

    optical_depth        = 1

  !  if (tau_eff .ge. 1.)  then
  !     optical_depth     = 1
  !  else
  !     optical_depth     = 0
  !  end if
    select case(optical_depth)
    case(1)
       Fz = 4._x_precision * c**2 * T**4 /(27._x_precision * sqrt(3._x_precision) &
            * (K_ff + K_e) * Sigma * Sigma_0)
    case (0)
    !   Fz = 4._x_precision * rs * E_ff * H / (Omega_0 * Sigma_0)

       Fz = 4._x_precision * rs * E_ff * H / (Omega_0 * Sigma_0)
    end select
    Q_plus              = 3._x_precision  * rs**2 * nu * Omega**2 * Omega_0**2
    !  Q_minus             = Fz  / Sigma * (1 / (3._x_precision * nu_0) * Omega_0 * rs**2)
    Q_minus             = Fz  / Sigma 

    f                   = Q_plus - Q_minus
    
  end subroutine variables
  
  
  !-------------------------------------------------------------------------
  ! Dichotomic function in order to determine the change of sign in a given
  ! interval [Smin,Smax] with an epsilon precision
  !-------------------------------------------------------------------------
  real(kind=x_precision) function dichotomy(Smin, Smax, eps, T, omega, sigma_0, Omega_0,rs, T_0, rho_0)
    use mod_read_parameters
    use mod_constants
    use mod_variables
    implicit none
    
    integer                                                  :: j = 0   
    real(kind=x_precision),intent(inout)                     :: Smin,Smax   
    real(kind=x_precision),intent(in)                        :: eps         
    real(kind=x_precision),intent(in)                        :: T           
    real(kind=x_precision),intent(in)                        :: omega
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
    real(kind=x_precision)                                   :: f_dichotomy
    real(kind=x_precision)                                   :: E_ff
    real(kind=x_precision)                                   :: Fz
    real(kind=x_precision),intent(in)                        :: Sigma_0
    real(kind=x_precision),intent(in)                        :: Omega_0
    real(kind=x_precision),intent(in)                        :: T_0
    real(kind=x_precision),intent(in)                        :: rho_0
    real(kind=x_precision),intent(in)                        :: rs

    !-------------------------------------------------------------------------
    ! N-> Number of iterations for the dichotomy
    ! Smin, Smax -> Starting range
    ! eps -> Precision
    ! T-> Fixed variable
    !-------------------------------------------------------------------------
    dichotomy             = (Smin+Smax)/2.
    j = 0
    call variables(T, Smin, Omega, H, rho, cs, nu, Q_plus, Q_minus, K_ff, K_e,&
         tau_eff, P_rad, P_gaz,E_ff,Fz,f_min,Sigma_0, Omega_0,rs,T_0,rho_0)
    
    call variables(T, Smax, Omega, H, rho, cs, nu, Q_plus, Q_minus, K_ff, K_e, &
         tau_eff, P_rad, P_gaz,E_ff,Fz,f_max,Sigma_0, Omega_0,rs,T_0, rho_0)
    
    call variables(T, dichotomy, Omega, H, rho, cs, nu, Q_plus, Q_minus, K_ff,&
         K_e, tau_eff, P_rad, P_gaz,E_ff,Fz,f_dichotomy,Sigma_0, Omega_0,rs,T_0, rho_0)
    
    write(*,*)'fmin = ',f_min
    write(*,*)'fmax = ',f_max
    
    if ( f_max * f_min .gt. 0.) then
       write(*,*)'This function image does not switch its sign in this particular interval.'
    endif

    if( f_max * f_min .lt. 0.) then
       iteration:do while ( dabs( Smax - Smin ) .ge. eps .and. j .lt. 1000)
          if(f_min * f_dichotomy .lt. 0.) then

 !   write(*,*)'fmin = ',f_min
 !   write(*,*)'fmax = ',f_max
 !   write(*,*)' '

             Smax         = dichotomy
             
             call variables(T, Smax, Omega, H, rho, cs, nu, Q_plus, Q_minus, K_ff, K_e, &
                  tau_eff, P_rad, P_gaz,E_ff,Fz,f_max,Sigma_0, Omega_0,rs,T_0, rho_0)
             
          else
             Smin         = dichotomy
             
             call variables(T, Smin, Omega, H, rho, cs, nu, Q_plus, Q_minus, K_ff, K_e,&
                  tau_eff, P_rad, P_gaz,E_ff,Fz,f_min,Sigma_0, Omega_0,rs,T_0, rho_0)
          endif
          
          dichotomy       = (Smin + Smax) * 1.0 / 2._x_precision

          call variables(T, dichotomy, Omega, H, rho, cs, nu, Q_plus, Q_minus, K_ff,&
               K_e, tau_eff, P_rad, P_gaz,E_ff,Fz,f_dichotomy,Sigma_0, Omega_0,rs,T_0, rho_0)
          j               = j + 1
    
       end do iteration
    endif
  
  end function dichotomy
  
  
  !-------------------------------------------------------------------------
  !Subroutine in order to display parameters
  !-------------------------------------------------------------------------
  subroutine display_parameters()
    type(parameters)                            :: param
    call get_parameters(param)
    write(*,*)'           Input Parameters             '
    write(*,*)'****************************************'
    write(*,"(' BH_mass     =',1p,E12.4)") param%M
    write(*,"(' Mdot        =',1p,E12.4)") param%Mdot
    write(*,"(' rmax        =',1p,E12.4)") param%rmax
    write(*,"(' alpha       =',1p,E12.4)") param%alpha
    write(*,"(' X           =',1p,E12.4)") param%X
    write(*,*)'****************************************'
    read(*,*)
  end subroutine display_parameters

  
  !-------------------------------------------------------------------------
  !Subroutine in order to display initial variables
  !-------------------------------------------------------------------------
  subroutine display_initial_variables(rs, rmin, Mdot_0, Sigma_0, Omega_0, T_0, rho_0)
    implicit none
    
    real(kind = x_precision), intent(in)                     :: rs
    real(kind = x_precision), intent(in)                     :: rmin
    real(kind = x_precision), intent(in)                     :: Mdot_0
    real(kind = x_precision), intent(in)                     :: Sigma_0
    real(kind = x_precision), intent(in)                     :: Omega_0
    real(kind = x_precision), intent(in)                     :: T_0
    real(kind = x_precision), intent(in)                     :: rho_0

    Type(parameters)                                         :: param
    call get_parameters(param)
    !-----------------------------------------------------------------------

    write(*,*)'           Initial Variables            '
    write(*,*)'****************************************'
    write(*,"(' Temp_0      =',1p,E12.4)") T_0
    write(*,"(' Sigma_0     =',1p,E12.4)") Sigma_0
    write(*,"(' Omega_0     =',1p,E12.4)") Omega_0
    write(*,"(' Omega_max   =',1p,E12.4)") sqrt(G*param%M / rmin**3)
    write(*,"(' H_0         =',1p,E12.4)") rs
    write(*,"(' Mdot_0      =',1p,E12.4)") Mdot_0
    write(*,"(' rho_0       =',1p,E12.4)") rho_0
    write(*,"(' rmin        =',1p,E12.4)") rmin
    write(*,"(' rs          =',1p,E12.4)") rs
    write(*,*)'****************************************'

    read(*,*)
  end subroutine display_initial_variables


  !-------------------------------------------------------------------------
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
