subroutine ureal(imol, iatm, k, box, rcut, kappa, charge, A, e_real)
  implicit none
  integer(4), intent(in) :: imol, iatm, k
  real(8),    intent(in) :: A(k, iatm, imol)
  real(8),    intent(in) :: box, rcut
  real(8),    intent(in) :: kappa
  real(8),    intent(in) :: charge(iatm, imol)
  real(8),    intent(inout) :: e_real(imol, imol)
  real, parameter  :: FourPIeps0inv = 138.935456 
  integer(4) :: iimol, iiatm, jjmol, jjatm, ik
  real(8), dimension(3) :: pos0, pos
  real(8) :: chargeprod, r, ereal_i_ii_j_jj
  integer(4) :: nx, ny, nz 


  !write(*,"(a)")"** write from fortran"
  !write(*,"(a)")"molid atmid position "
  ! Memory-ordered access (row-major)
  !$omp parallel private(pos0, pos, r, chargeprod, ereal_i_ii_j_jj) 
  !$omp do schedule(static, 1)
  do iimol = 1, imol
    do jjmol = iimol+1, imol
      do iiatm = 1, iatm
        do jjatm = 1, iatm

          pos0 = A(:, iiatm, iimol) - A(:, jjatm, jjmol)
          !pos0 = pos0 - nint(pos0/box)*box
          do nx = -1, 1
            do ny = -1, 1
              do nz = -1, 1
                pos =  pos0 + [box*nx, box*ny, box*nz]
                r = sqrt(sum(pos**2.0e0))
                !if ( r>rcut) cycle        

                chargeprod   =charge(iiatm, iimol)*charge(jjatm, jjmol)
                ereal_i_ii_j_jj = chargeprod * FourPIeps0inv  * (1.0e0 -erf(kappa*r))/r
                e_real(jjmol, iimol) = e_real(jjmol, iimol) + ereal_i_ii_j_jj
                !write(*, *) iimol, jjmol, iiatm, jjatm, nx, ny, nz, r, elj_i_ii_j_jj
              end do
            end do
          end do
        end do
      end do
      !write(*, *) iimol, jjmol, e_real(jjmol, iimol) 
    end do
  end do
  !$omp end do
  !$omp end parallel


end subroutine ureal
