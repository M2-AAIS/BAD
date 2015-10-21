module mod_S_curve
  use mod_read_parameters
  use mod_constants
  use mod_variables
  implicit none
  
contains
  !-------------------------------------------------------------------------
  ! subroutine to find the result of equation Q+ + Q- = 0
  ! output : datafile with temperature and surface density
  ! Range of surface density [Smin, Smax] and temperature [t_min, t_max]
  ! values can be changed
  !-------------------------------------------------------------------------
  subroutine curve()  
    implicit none
    integer::i
    real(kind = x_precision)                               :: temp
    real(kind = x_precision), parameter                    :: t_min = 0.0d0
    real(kind = x_precision), parameter                    :: t_max = 10.0d0
    integer,                  parameter                    :: nb_it = 10
    real(kind = x_precision)                               :: sigma = 0.0d0
    real(kind = x_precision)                               :: Smin=0.0d0
    real(kind = x_precision)                               :: Smax=10.0d0
    real(kind = x_precision)                               :: eps=1.0d0
    real(kind = x_precision)                               :: omega = 0.0d0
    real(kind = x_precision)                               :: r = 0.0d0
    type(parameters)                                       :: param
    !----------------------------------------------------------------------
    
    call get_parameters(param)
  !-------------------------------------------------------------------------
  ! Test for 1 value of r
  !-------------------------------------------------------------------------
    r              =  6*G*param%M/(c**2)     
    omega          =  G*param%M/r**3

    open(18,file = 'Temperature_Sigma.dat',status='unknown',action='readwrite') 
    write(18,*)"T Sigma"
    do i           = 1 , nb_it
       temp        = (t_max-t_min)/(nb_it-1)*(i-1) + t_min
       sigma       = dichotomy(Smin, Smax, eps, temp, omega)
       
       write(18,*)temp," ",sigma
    enddo
    close(18)
  end subroutine curve


  !-------------------------------------------------------------------------
  !Subroutine for resolving a quadratic equation as long as the solutions are
  !a set of real numbers
  !-------------------------------------------------------------------------
  subroutine quadratic(coeff_a,coeff_b,coeff_c,sol)

    implicit none
    real(kind=x_precision), intent(in)                     ::coeff_a,coeff_b,coeff_c
    real(kind=x_precision), intent(out)                    ::sol
    real(kind=x_precision)                                 ::delta=0.
    
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
  real(kind=x_precision),intent(in)                        ::T,Sigma,Omega
  real(kind=x_precision)                                   ::coeff_a=0.,coeff_b=0.,coeff_c=0.
  real(kind=x_precision)                                   ::H=0.
  real(kind=x_precision)                                   ::rho=0.
  real(kind=x_precision)                                   ::cs=0.
  real(kind=x_precision)                                   ::nu=0.
  real(kind=x_precision)                                   ::Q_plus=0.
  real(kind=x_precision)                                   ::K_ff=0.
  real(kind=x_precision)                                   ::K_e=0.
  real(kind=x_precision)                                   ::E_ff=0.
  real(kind=x_precision)                                   ::tau_eff=0.
  real(kind=x_precision)                                   ::Fz=0.
  real(kind=x_precision)                                   ::Q_minus=0.
  real(kind=x_precision)                                   ::P_rad=0.
  real(kind=x_precision)                                   ::P_gaz=0.
  real(kind=x_precision)                                   ::Mdot_0=0.
  real(kind=x_precision)                                   ::Sigma_0=0.
  real(kind=x_precision)                                   ::rmin=0.
  real(kind=x_precision)                                   ::rs=0.
  integer                                                  ::optical_depth = 0
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

  select case(optical_depth)
  case(1:)
     Fz = 4 * c**2 * T**4/(27. * sqrt(3.) * (K_ff + K_e) * Sigma * Sigma_0)
  case default
     Fz = 6.22d20 * 2 / state_0%Omega_0 * state_0%rho_0 * state_0%T_0**2 * H * rho **2 * T**2
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
    
  integer                                                  ::j=0,N=100   
  real(kind=x_precision),intent(inout)                     ::Smin,Smax   
  real(kind=x_precision),intent(in)                        ::eps         
  real(kind=x_precision),intent(in)                        ::T           
  real(kind=x_precision),intent(in)                        ::omega
  !-------------------------------------------------------------------------
  ! N-> Number of iterations for the dichotomy
  ! Smin, Smax -> Starting range
  ! eps -> Precision
  ! T-> Fixed variable
  !-------------------------------------------------------------------------
  
  dichotomy             = (Smin+Smax)/2.
  if ( f(T,Smin,omega) * f(T,Smax,omega) .gt. 0.) then
     write(*,*)'This function image does not switch its sign in this particular interval.'
  endif
  if( f(T,Smin,omega) * f(T,Smax,omega) .lt. 0.) then
     iteration:do while ( abs( f(T,dichotomy,omega) ) .ge. eps .and. j .lt. N)
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

end module mod_S_curve
