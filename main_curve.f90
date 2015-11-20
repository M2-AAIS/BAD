program black_hole_s_curve
  use mod_constants
  use mod_variables
  use mod_read_parameters
  use mod_s_curve

  implicit none

  call get_parameters() 
  call init_variable_0()
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
