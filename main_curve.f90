program black_hole_s_curve
  use mod_constants
  use mod_read_parameters
  use mod_s_curve


! CAREFUL !!! The output arrays are for the NON dimensioned values of the surface density S (not Sigma) and the Temperature  whereas the produced data files used when ploting are dimensioned, in log10, for Sigma (not S) and the Temperature.

  implicit none

  integer                                    :: i

  integer                 ,parameter         :: nb_it = 1000
  real(kind = x_precision),parameter         :: eps   = 1.0d-4
  real(kind = x_precision),parameter         :: t_min = 2.0d-1
  real(kind = x_precision),parameter         :: t_max = 4.49d0
  real(kind = x_precision),parameter         :: dt    = (t_max-t_min) / (nb_it-1)
  real(kind = x_precision),parameter         :: S_min = 4.5d1 ! Greater values lead to FPE
  real(kind = x_precision),parameter         :: S_max = 2.5d3

  real(kind = x_precision),dimension(n_cell) :: temperature
  real(kind = x_precision),dimension(n_cell) :: s


  call get_parameters()

  call curve( nb_it, eps, t_min, t_max, dt, s_min, s_max, temperature, s)
  do i = 1, n_cell

     write(*,*)'**********************************'
     write(*,*)'T = ',temperature(i)
     write(*,*)'S = ',s(i)
     write(*,*)'**********************************'

  enddo

end program black_hole_s_curve
