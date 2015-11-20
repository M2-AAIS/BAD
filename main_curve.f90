program black_hole_s_curve
  use mod_constants
  use mod_variables
  use mod_read_parameters
  use mod_s_curve

  implicit none

  call get_parameters() 
  call init_variable_0()
  call curve()
  
end program black_hole_s_curve
