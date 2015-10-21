module mod_S_curve
  use mod_read_parameters
  use mod_constants
  use mod_variables
  implicit none
  

contains
  subroutine curve()  
    implicit none
    integer::i
    real(kind=x_precision)::temp
    real(kind=x_precision), parameter::t_min = 0.
    real(kind=x_precision), parameter::t_max = 10.
    integer, parameter::nb_it=2
    real(kind=x_precision)::sigma=0.
    real(kind=x_precision)::Smin=0.
    real(kind=x_precision)::Smax=10.
    real(kind=x_precision)::eps=1.
    real(kind=x_precision)::omega
    real(kind=x_precision)::r = 1.
    type(parameters)::param
    call get_parameters(param)
    omega =  G*param%M/r**3

    open(18,file='Temperature_Sigma.dat',status='unknown',action='readwrite') 
    write(18,*)"T Sigma"
    
    do i=1,nb_it
       temp = (t_max-t_min)/(nb_it-1)*(i-1) + t_min
       sigma = dichotomy(Smin, Smax, eps, temp, omega)
       write(18,*)temp," ",sigma
    enddo
close(18)
end subroutine curve



!Subroutine for resolving a quadratic equation as long as the solutions are a set of real numbers
subroutine quadratic(coeff_a,coeff_b,coeff_c,sol)
implicit none
real(kind=x_precision),intent(in)::coeff_a,coeff_b,coeff_c
real(kind=x_precision),intent(out)::sol
real(kind=x_precision)::delta=0.
delta=coeff_b**2-4*coeff_a*coeff_c
if (delta .lt. 0.) then
   write(*,*)'No solutions in the R field.'
else
   sol=-0.5d0*(coeff_b+sign(sqrt(delta),coeff_b))/coeff_a
end if
end subroutine quadratic

!_________________________________________________________________________

!_________________________________________________________________________









real(kind=x_precision) function f(T, Sigma, Omega)

  implicit none
 
  type(parameters)::param
  real(kind=x_precision),intent(in)::T,Sigma,Omega
  real(kind=x_precision)::coeff_a=0.,coeff_b=0.,coeff_c=0.
  real(kind=x_precision)::A_coeff=1.
  real(kind=x_precision)::H=0.
  real(kind=x_precision)::rho=0.
  real(kind=x_precision)::cs=0.
  real(kind=x_precision)::nu=0.
  real(kind=x_precision)::Q_plus=0.
  real(kind=x_precision)::K_ff=0.
  real(kind=x_precision)::K_e=0.
  real(kind=x_precision)::E_ff=0.
  real(kind=x_precision)::tau_eff=0.
  real(kind=x_precision)::Fz=0.
  real(kind=x_precision)::Q_minus=0.
  real(kind=x_precision)::P_rad=0.
  real(kind=x_precision)::P_gaz=0.
  real(kind=x_precision)::Mdot_0=0.
  real(kind=x_precision)::Sigma_0=0.
  real(kind=x_precision)::rmin=0.
  real(kind=x_precision)::rs=0.
  integer::optical_depth=0 !Indicates the thickness of the medium
  call get_parameters(param)
  
  rs=2*G*param%M/(c**2)
  rmin=3*rs
  Mdot_0 = param%Mdot 
  Sigma_0 = Mdot_0 /(state_0%Omega_0*rs**2*2*pi)
  coeff_a = (Omega**2*state_0%Omega_0**2*Sigma*Sigma_0)/2.
  coeff_b = (-1./3.)*cst_rad*T**4*state_0%T_0**4
  coeff_c = (-1.* kmp * T*state_0%T_0 *Sigma * Sigma_0)/(2.*mu)
  call quadratic(coeff_a,coeff_b,coeff_c,H)
  !write(*,*)'H= ',H
  rho=Sigma/H
  cs = Omega*H
  nu = param%alpha*cs*H
  P_rad= T**4
  P_gaz = rho*T
  K_ff = 6.13e22 * state_0%rho_0 * rho * (state_0%T_0 * T)**(-3.5)
  K_e = 0.2*(1.+param%X)
  E_ff = 6.22e20 * (state_0%rho_0 * rho)**2 * sqrt(state_0%T_0 * T)
  tau_eff = 0.5 * sqrt(K_e * K_ff) * Sigma * Sigma_0    
  !if (tau_eff .ge. 1.) then
  !  optical_depth=1
  !else
  optical_depth=0
  !end if
  !select case(optical_depth)
  !  case(1:)
  !    Fz = (2.0 * cst_rad * c * (T_0 * T)**4) / (3.0 * (K_ff * K_e) * Sigma*Sigma_0
  ! case default
  Fz = E_ff*A_coeff*H*state_0%H_0
  !end select
  Q_minus = 2.*Fz /Sigma
  !write(*,*)'Fz= ',Fz
  !write(*,*)'Q_minus= ',Q_minus
  Q_plus = 9./4.* nu *Omega**2
  ! write(*,*)'Q_plus= ',Q_plus
  f = Q_plus - Q_minus
  !write(*,*)'f= ',f
  end function f


real(kind=x_precision) function dichotomy(Smin, Smax, eps, T, omega)
  use mod_read_parameters
  use mod_constants
  use mod_variables
  implicit none
  integer::j=0,N=100            ! Number of iterations for the dichotomy
  real(kind=x_precision),intent(inout)::Smin,Smax   ! Starting range
  real(kind=x_precision),intent(in)::eps            ! Precision
  real(kind=x_precision),intent(in)::T              ! Fixed variable
  real(kind=x_precision),intent(in)::omega

  dichotomy=(Smin+Smax)/2.
  if ( f(T,Smin,omega) * f(T,Smax,omega) .gt. 0.) then
     write(*,*)'This function image does not switch its sign in this particular interval.'
  endif
  if( f(T,Smin,omega) * f(T,Smax,omega) .lt. 0.) then
     iteration:do while ( abs( f(T,dichotomy,omega) ) .ge. eps .and. j .lt. N)
        if( f(T,Smin,omega) * f(T,dichotomy,omega) .lt. 0.) then
           Smax=dichotomy
        else
           Smin=dichotomy
      endif
      dichotomy=(Smin+Smax)/2.
      j=j+1
    end do iteration
  endif
end function dichotomy

end module mod_S_curve
