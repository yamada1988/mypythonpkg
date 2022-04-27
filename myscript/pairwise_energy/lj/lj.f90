subroutine lj(imol, iatm, k, box, rcut, sigma, eps, A, e_lj)
  implicit none
  integer(4), intent(in) :: imol, iatm, k
  !convert column-major to row-major
  real(8),    intent(in) :: A(k, iatm, imol)
  real(8),    intent(in) :: box, rcut
  real(8),    intent(in) :: sigma(iatm, imol)
  real(8),    intent(in) :: eps(iatm, imol)
  real(8),    intent(inout) :: e_lj(imol, imol)
  integer(4) :: iimol, iiatm, jjmol, jjatm, ik
  real(8), dimension(3) :: pos
  real(8) :: e, s, r, elj_i_ii_j_jj
 

  !write(*,"(a)")"** write from fortran"
  !write(*,"(a)")"molid atmid position "
  ! Memory-ordered access (row-major)
  !$omp parallel private(pos, r, e, s, elj_i_ii_j_jj) 
  !$omp do schedule(static, 1)
  do jjmol = 1, imol
    do iimol = jjmol+1, imol
      do jjatm = 1, iatm
        do iiatm = 1, iatm

          pos = A(:, iiatm, iimol) - A(:, jjatm, jjmol)
          pos = pos - nint(pos/box)*box
          r = sqrt(sum(pos**2.0e0))
          if (r<=rcut) cycle
       
          e   = sqrt(eps(iiatm, iimol)*eps(jjatm, jjmol) )
          s = 0.50e0*(sigma(iiatm, iimol) + sigma(jjatm, jjmol))
          elj_i_ii_j_jj = 4.0*e*((s/r)**12-(s/r)**6)
          e_lj(jjmol, iimol) = e_lj(jjmol, iimol) + elj_i_ii_j_jj
        end do
      end do
      !write(*, *) iimol, jjmol, e_lj(jjmol, iimol) 
    end do
  end do
  !$omp end do
  !$omp end parallel


end subroutine lj
