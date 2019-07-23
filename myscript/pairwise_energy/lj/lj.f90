subroutine lj(k, iatm, imol, box, rcut, sigma, eps, A, e_lj)
  implicit none
  integer(4), intent(in) :: iatm, imol, k
  real(8),    intent(in) :: A(k, iatm, imol)
  real(8),    intent(in) :: box, rcut
  real(8),    intent(in) :: sigma(iatm, imol)
  real(8),    intent(in) :: eps(iatm, imol)
  real(8),    intent(inout) :: e_lj(imol, imol)
  integer(4) :: iimol, iiatm, jjmol, jjatm, ik
  real(8), dimension(3) :: pos
  real(8) :: e, s, r
 

  !write(*,"(a)")"** write from fortran"
  !write(*,"(a)")"molid atmid position "
  ! Memory-ordered access (row-major)
  do iimol = 1, imol
    do jjmol = iimol+1, imol
      do iiatm = 1, iatm
        do jjatm = 1, iatm

          pos = A(:, iiatm, iimol) - A(:, jjatm, jjmol)
          pos = pos - nint(pos/box)*box
          r = sqrt(sum(pos**2.0e0))
          if (r<=rcut) cycle
       
          e   = sqrt(eps(iiatm, iimol)*eps(jjatm, jjmol) )
          s = 0.50e0*(sigma(iiatm, iimol) + sigma(jjatm, jjmol))
          e_lj(jjmol, iimol) = e_lj(jjmol, iimol) + 4.0*e*((s/r)**12-(s/r)**6)
        end do
      end do
      !write(*, *) iimol, jjmol, e_lj(jjmol, iimol) 
    end do
  end do


end subroutine lj
