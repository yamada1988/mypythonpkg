subroutine calc(imol, iatm, k, A)
  implicit none
  integer(4), intent(in) :: iatm, imol, k
  real(8),    intent(in) :: A(k, iatm, imol)
  integer(4) :: iimol, iiatm, jjmol, jjatm, ik
 

  !write(*,"(a)")"** write from fortran"
  !write(*,"(a)")"molid atmid position "
  ! Memory-ordered access (row-major)
  do jjmol = 1, imol
      do jjatm = 1, iatm

      write(*, *) jjmol, jjatm, A(:, jjatm, jjmol) 
    end do
  end do


end subroutine calc
