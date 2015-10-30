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
    integer::i = 0
    real(kind = x_precision)                               :: temp  = 0.0d0
    real(kind = x_precision), parameter                    :: t_min = 1.2d-2
    real(kind = x_precision), parameter                    :: t_max = 1.0d0
    integer,                  parameter                    :: nb_it = 2
    real(kind = x_precision)                               :: sigma = 0.0d0
    real(kind = x_precision)                               :: Smin  = 0.0d0
    real(kind = x_precision)                               :: Smax  = 0.0d0
    real(kind = x_precision)                               :: eps   = 1.0d-4
    real(kind = x_precision)                               :: omega = 0.0d0
    real(kind = x_precision)                               :: r     = 5.0d0
    type(parameters)                                       :: param
    !----------------------------------------------------------------------
    
    call get_parameters(param)
    call display_parameters()
    call initial_variables()
    !-------------------------------------------------------------------------
    ! Test for 1 value of r
    !-------------------------------------------------------------------------
    r              = 6*G*param%M/(c**2)     
    omega          = G*param%M/r**3

    open(18,file   = 'Temperature_Sigma.dat',status='unknown',action='readwrite') 
    write(18,*)"T Sigma"
    do i           = 1 , nb_it
       Smin        = 1.0d-5
       Smax        = 1.0d-0
       temp        = (t_max-t_min)/(nb_it-1)*(i-1) + t_min
       sigma       = dichotomy(Smin, Smax, eps, temp, omega)
       call display_variables(temp,Omega,r,sigma)
       write(18,*)temp," ",sigma
    enddo
    close(18)
  end subroutine curve


  !-------------------------------------------------------------------------
  !Subroutine for resolving a quadratic equation as long as the solutions are
  !a set of real numbers
  !-------------------------------------------------------------------------
  subroutine quadratic(coeff_a, coeff_b, coeff_c, sol)

    implicit none
    real(kind=x_precision), intent(in)                     :: coeff_a,coeff_b,coeff_c
    real(kind=x_precision), intent(inout)                  :: sol
    real(kind=x_precision)                                 :: delta=0.
    delta          = coeff_b**2 - 4 * coeff_a * coeff_c
    if (delta .lt. 0.) then
      write(*,*)'No solutions in the R field.'
    else
      sol          = -0.5d0 *(coeff_b+sign(sqrt(delta),coeff_b))/coeff_a
    end if
end subroutine quadratic



real(kind=x_precision) function f(T, Sigma, Omega)
  implicit none
 
  type(parameters)::param
  real(kind = x_precision),intent(in)                      :: T,Sigma,Omega
  real(kind = x_precision)                                 :: coeff_a=0.,coeff_b=0.,coeff_c=0.
  real(kind = x_precision)                                 :: H     = 0.
  real(kind = x_precision)                                 :: rho   = 0.
  real(kind = x_precision)                                 :: cs    = 0.
  real(kind = x_precision)                                 :: nu    = 0.
  real(kind = x_precision)                                 :: Q_plus= 0.
  real(kind = x_precision)                                 :: K_ff  = 0.
  real(kind = x_precision)                                 :: K_e   = 0.
  real(kind = x_precision)                                 :: E_ff  = 0.
  real(kind = x_precision)                                 :: tau_eff=0.
  real(kind = x_precision)                                 :: Fz     =0.
  real(kind = x_precision)                                 :: Q_minus=0.
  real(kind = x_precision)                                 :: P_rad  =0.
  real(kind = x_precision)                                 :: P_gaz = 0.
  real(kind = x_precision)                                 :: Mdot_0 =0.
  real(kind = x_precision)                                 :: Sigma_0=0.
  real(kind = x_precision)                                 :: rmin  = 0.
  real(kind = x_precision)                                 :: rs    = 0.
  integer                                                  :: optical_depth = 0
  call get_parameters(param)
  
  rs                   = 2._x_precision * G * param%M/(c**2)
  rmin                 = 3._x_precision * rs
  Mdot_0               = param%Mdot 
  Sigma_0              = Mdot_0 /(state_0%Omega_0 * rs**2 * 2 * pi)
  coeff_a              = (Omega**2 * state_0%Omega_0**2 * Sigma * Sigma_0)/2.
  coeff_b              = (-1._x_precision/3._x_precision) * cst_rad*T**4 * state_0%T_0**4
  coeff_c              = (-1._x_precision * kmp * T * state_0%T_0 * Sigma * Sigma_0)/(2._x_precision * mu)
  call quadratic(coeff_a,coeff_b,coeff_c,H)
  rho                  = Sigma/H
  cs                   = Omega * H
  nu                   = param%alpha * cs * H
  P_rad                = T**4
  P_gaz                = rho * T
  K_ff                 = 6.13d22 * state_0%rho_0 * rho * (state_0%T_0 * T)**(-3.5)
  K_e                  = 0.2_x_precision * (1_x_precision + param%X)
  E_ff                 = 6.22d20 * (state_0%rho_0 * rho)**2 * sqrt(state_0%T_0 * T)
  tau_eff              = 0.5_x_precision * sqrt(K_e * K_ff) * Sigma * Sigma_0    
  !-------------------------------------------------------------------------
  !Select the case for the opticaly depth
  !-------------------------------------------------------------------------
  if (tau_eff .ge. 1.)  then
     optical_depth     = 1
  else
     optical_depth     = 0
  end if
  !write(*,*)'tau',tau_eff
  select case(optical_depth)
  case(1)
     Fz = 4._x_precision * c**2 * T**4/(27. * sqrt(3._x_precision) * (K_ff + K_e) * Sigma * Sigma_0)
  case default
     Fz = 6.22d20 * 2._x_precision / state_0%Omega_0 * state_0%rho_0 *  H * rho**2 * sqrt(T*state_0%T_0)
  end select
  Q_minus             = 2._x_precision * Fz /Sigma
  Q_plus              = 9._x_precision /4._x_precision * nu * Omega**2
  f                   = Q_plus - Q_minus
  end function f


  !-------------------------------------------------------------------------
  ! Function in order to determine the change of sign in a given interval
  ! [Smin,Smax] with an epsilon precision
  !-------------------------------------------------------------------------
  real(kind=x_precision) function dichotomy(Smin, Smax, eps, T, omega)
    use mod_read_parameters
    use mod_constants
    use mod_variables
    implicit none
    
  integer                                                  :: j = 0   
  real(kind=x_precision),intent(inout)                     :: Smin,Smax   
  real(kind=x_precision),intent(in)                        :: eps         
  real(kind=x_precision),intent(in)                        :: T           
  real(kind=x_precision),intent(in)                        :: omega
  !-------------------------------------------------------------------------
  ! N-> Number of iterations for the dichotomy
  ! Smin, Smax -> Starting range
  ! eps -> Precision
  ! T-> Fixed variable
  !-------------------------------------------------------------------------
  dichotomy             = (Smin+Smax)/2.
  j = 0
  if ( f(T,Smin,omega) * f(T,Smax,omega) .gt. 0.) then
  endif
  if( f(T,Smin,omega) * f(T,Smax,omega) .lt. 0.) then
     iteration:do while ( dabs( Smax - Smin ) .ge. eps)
        if( f(T,Smin,omega) * f(T,dichotomy,omega) .lt. 0.) then
           Smax         = dichotomy
        else
           Smin         = dichotomy
        endif
        dichotomy       = (Smin+Smax)/2.0d0
        j               = j + 1
     end do iteration
  endif
end function dichotomy


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

subroutine initial_variables()
  write(*,*)'           Initial Variables            '
  write(*,*)'****************************************'
  write(*,"(' Temp_0      =',1p,E12.4)") state_0%T_0
  write(*,"(' Simga_0     =',1p,E12.4)") state_0%S_0/state_0%x_0
  write(*,"(' Omega_0     =',1p,E12.4)") state_0%Omega_0
  write(*,"(' H_0         =',1p,E12.4)") state_0%H_0
  write(*,"(' Mdot_0      =',1p,E12.4)") state_0%M_dot_0
  write(*,"(' rho_0       =',1p,E12.4)") state_0%rho_0

  write(*,*)'****************************************'
  read(*,*)
end subroutine initial_variables



subroutine display_variables(temp,Omega,r,sigma)
  implicit none
  real(kind = x_precision), intent(in)                     ::temp
  real(kind = x_precision), intent(in)                     ::Omega
  real(kind = x_precision), intent(in)                     ::r
  real(kind = x_precision), intent(in)                     ::sigma
  !real(kind = x_precision), intent(in)                     ::sol

  write(*,*)'               Variables                '
  write(*,*)'****************************************'
  write(*,"(' r           =',1p,E12.4)") r
  write(*,"(' Temp        =',1p,E12.4)") temp
  write(*,"(' Sigma       =',1p,E12.4)") sigma
  write(*,"(' Omega       =',1p,E12.4)") Omega
  !write(*,*)'H         = ', sol
  write(*,*)'****************************************'
  read(*,*)

end subroutine display_variables


end module mod_S_curve
