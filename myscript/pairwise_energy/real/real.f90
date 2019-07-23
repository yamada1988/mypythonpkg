subroutine ureal(k, iatm, imol, box, rcut, kappa, charge, A, e_real)
  implicit none
  integer(4), intent(in) :: iatm, imol, k
  real(8),    intent(in) :: A(k, iatm, imol)
  real(8),    intent(in) :: box, rcut
  real(8),    intent(in) :: kappa
  real(8),    intent(in) :: charge(iatm, imol)
  real(8),    intent(inout) :: e_real(imol, imol)
  real, parameter  :: FourPIeps0inv = 138.935456 
  integer(4) :: iimol, iiatm, jjmol, jjatm, ik
  real(8), dimension(3) :: pos
  real(8) :: chargeprod, r
 

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
          if ( r>rcut) cycle        

          chargeprod   =charge(iiatm, iimol)*charge(jjatm, jjmol)
          e_real(jjmol, iimol) = e_real(jjmol, iimol) + chargeprod * FourPIeps0inv  * (1.0e0 - erf(kappa*r))/r
        end do
      end do
      write(*, *) iimol, jjmol, e_real(jjmol, iimol) 
    end do
  end do


end subroutine ureal
