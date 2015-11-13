program black_hole_diffusion
  use mod_constants
  use mod_variables
  use mod_read_parameters
  use mod_S_curve

  implicit none

  call get_parameters() 
  call init_variable_0()
  call curve()
  
end program black_hole_diffusion
