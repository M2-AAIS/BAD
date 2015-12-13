module mod_maps
  use mod_constants
  use mod_read_parameters
  use mod_variables

  implicit none

  integer           :: nb_S  ! The number of points along Sigma axis
  integer           :: nb_T  ! The number of points along T axis
  real(x_precision) :: T_min ! Lowest value of Temperature range
  real(x_precision) :: T_max ! Highest value of Temperature range
  real(x_precision) :: S_min ! Lowest value of Sigma range
  real(x_precision) :: S_max ! Highest value of Sigma range


contains

  subroutine set_conditions(nbS, nbT, Tmin, Tmax, Smin, Smax)
    implicit none

    integer, intent(in) :: nbS
    integer, intent(in) :: nbT
    real(x_precision), intent(in) :: Tmin ! Lowest value of Temperature range
    real(x_precision), intent(in) :: Tmax ! Highest value of Temperature range
    real(x_precision), intent(in) :: Smin ! Lowest value of Sigma range
    real(x_precision), intent(in) :: Smax ! Highest value of Sigma range

    nb_S = nbS
    nb_T = nbT
    T_min = Tmin
    T_max = Tmax
    S_min = Smin
    S_max = Smax

  end subroutine set_conditions

  ! Build the T, Sigma grid given the number of points along each axis and the ranges
  subroutine build_grid(Q_res, tau_res)
    implicit none

    real(x_precision), dimension(n_cell,nb_T,nb_S), intent(out) :: Q_res   ! The Q+-Q- grids, one for each position in the disk
    real(x_precision), dimension(n_cell,nb_T,nb_S), intent(out) :: tau_res ! The tau grids, one for each position in the disk

    type(state)       :: s   ! A temporary state to compute the variables
    real(x_precision) :: dT  ! Temperature steps in the grid
    real(x_precision) :: dS  ! Sigma steps in the grid
    integer           :: i,j ! Loop counters

    ! Compute the Temperature and Sigma steps
    dT = (log(T_max) - log(T_min)) / (nb_T - 1)
    dS = (log(S_max) - log(S_min)) / (nb_S - 1)

    do i = 1, nb_T
       s%T = exp(dT * (i-1) + log(T_min))
       do j = 1, nb_S
          s%S = exp(dS * (j-1) + log(S_min)) * x_state%x

          call compute_variables(s) ! Compute the variables in each position of the state
          Q_res(:,i,j) = s%Qplus - s%Qminus
          tau_res(:,i,j) = s%tau
       end do
    end do

  end subroutine build_grid

  subroutine save_data(Q_res, tau_res)
    implicit none

    real(x_precision), dimension(n_cell,nb_T,nb_S), intent(in) :: Q_res   ! The Q+-Q- grids, one for each position in the disk
    real(x_precision), dimension(n_cell,nb_T,nb_S), intent(in) :: tau_res ! The tau grids, one for each position in the disk

    character(len = 64) :: fname       ! Name of the output file
    character(len = 16) :: line_fmt    ! Format of lines in the output file
    character(len = 8)  :: cell_number ! Cell number, to use in filename
    integer             :: fid         ! File descriptor
    integer             :: ios         ! I/O status
    integer             :: i,j,k       ! Loop counters

    write(line_fmt,'(A,I4,A)') '(',nb_S,'(e11.3e2))'

    do k = 1, n_cell

      write(cell_number,'(I5.5)') k
      fid = 30 + k
      fname = 'maps/map_'//trim(cell_number)//'.dat'

      open(fid, file = fname, action='write', status = 'replace', iostat = ios)
      if (ios /= 0) then
        write(*,*)"Error while opening the ", fname," file."
        stop
      endif

      write(fid, fmt = '(A)') '# nb_T'
      write(fid, fmt = '(I5.5)') nb_T

      write(fid, fmt = '(A)') '# nb_S'
      write(fid, fmt = '(I5.5)') nb_S

      write(fid, fmt = '(A)') '# T bounds'
      write(fid, fmt = '(2(e11.3e2))') T_min*state_0%T_0, T_max*state_0%T_0

      write(fid, fmt = '(A)') '# Sigma bounds'
      write(fid, fmt = '(2(e11.3e2))') S_min*state_0%S_0, S_max*state_0%S_0

      write(fid, fmt = '(A)') '# Q grid'
      do i = 1, nb_T
        write(fid, fmt = line_fmt) (Q_res(k,nb_T+1-i,j), j = 1, nb_S)
      end do

      write(fid, fmt = '(A)') '# tau grid'
      do i = 1, nb_T
        write(fid, fmt = line_fmt) (tau_res(k,nb_T+1-i,j), j = 1, nb_S)
      end do

      close(fid)

    enddo

  end subroutine save_data

end module mod_maps
