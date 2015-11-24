program black_hole_s_curve
  use mod_constants
  use mod_variables
  use mod_read_parameters
  use mod_s_curve


! CAREFUL !!! The output arrays are for the NON dimensioned values of the surface density S (not Sigma) and the Temperature  whereas the produced data files used when ploting are dimensioned, in log10, for Sigma (not S) and the Temperature.

  implicit none

    integer                                                 :: i

    real(kind = x_precision),dimension(n_cell)              :: temperature
    real(kind = x_precision),dimension(n_cell)              :: s


  call get_parameters() 
  call init_variable_0()
  
     call curve( temperature,s)

  do i          = 1, n_cell

    
     write(*,*)'**********************************'
     write(*,*)'T = ',temperature(i)
     write(*,*)'S = ',s(i)
     write(*,*)'**********************************'

  enddo

  
end program black_hole_s_curve
