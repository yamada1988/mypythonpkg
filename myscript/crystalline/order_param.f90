! This program calculates the ow-ow rdf
! for pure water in NVT ensemble with periodic cubic boundaries
program ordrprm
    use xdr, only: xtcfile
    implicit none
    type(xtcfile) :: xtc
    real(kind=8), allocatable :: pos_c3(:,:,:), com_c3(:,:,:), vec_c3(:,:,:), cos_theta(:,:), P2(:)
    real(kind=8) :: vec_zaxis(3)
    real(kind=8) :: box_dim, r_max, v
    integer :: n_chain, n_mon, i_chain, i_mon, n_pos, n_t, n_m, natm_chain
    integer :: i, k

    call xtc % init("PE_020_072_short.xtc")

    n_chain = 72
    n_mon = 20
    n_m = 6
    natm_chain = n_mon * n_m + 2
    n_pos = n_chain * n_mon

    vec_zaxis(:) = [0.0E0, 0.0E0, 1.0E0]
    allocate(pos_c3(3, n_mon*2, n_chain))
    allocate(com_c3(3, n_mon*2-1, n_chain))
    allocate(vec_c3(3, n_mon*2-2, n_chain))
    allocate(cos_theta(n_mon*2-2, n_chain))
    allocate(P2(10))


    n_t = 0
    do while ( xtc % STAT == 0 )    
      n_t = n_t + 1 
      call xtc % read
      box_dim = xtc % box(1,1)
      r_max = box_dim / 2d0

      ! read position of c3 into pos_c3
      do i_chain = 1, n_chain
        do i_mon = 1, n_mon*2
          if  ( i_mon < 3 ) then
              k = ((i_mon-1)*3 + 1)     + (i_chain-1)*natm_chain
          else if ( i_mon >= 3 ) then
              k = ((i_mon-1)*3 + 1) + 1 + (i_chain-1)*natm_chain ! gro
          end if

          pos_c3(:, i_mon, i_chain)= xtc % pos(:, k)  !get the position of OW (every 3rd atom)
          !print *, n_t, i_chain, i_mon, pos_c3(:, i_mon, i_chain), k
          i = i + 1
        end do
      end do

      ! calculate com of c3 for each monomers
      do i_chain = 1, n_chain
        do i_mon = 1, n_mon*2-1
          com_c3(:, i_mon, i_chain) = (pos_c3(:, i_mon, i_chain) + pos_c3(:, i_mon+1, i_chain))*0.50E0
          !print *, n_t, i_chain, i_mon, com_c3(:, i_mon, i_chain)
        end do
      end do


      ! calculate vectors of com_c3
      do i_chain = 1, n_chain
        do i_mon = 1, n_mon*2-2
          vec_c3(:, i_mon, i_chain) = (com_c3(:, i_mon+1, i_chain) - com_c3(:, i_mon, i_chain))
          v = sqrt(sum(vec_c3(:, i_mon, i_chain)**2))
          vec_c3(:, i_mon, i_chain) = vec_c3(:, i_mon, i_chain) /v
          !print *, n_t, i_chain, i_mon, vec_c3(:, i_mon, i_chain)
        end do
      end do


      ! calculate angles between z-axis and vec_c3
      do i_chain = 1, n_chain
        do i_mon = 1, n_mon*2-2
          cos_theta(i_mon, i_chain) = DOT_PRODUCT(vec_c3(:, i_mon, i_chain), vec_zaxis)
          !print *, n_t, i_chain, i_mon, cos_theta(i_mon, i_chain)
        end do
      end do     


      ! calculate P2(cos(theta))
      P2(n_t) = 0.0E0
      do i_chain = 1, n_chain
        do i_mon = 1, n_mon*2-2
          P2(n_t) = P2(n_t) + 1.0E0/2.0E0*(3.0E0*cos_theta(i_mon, i_chain)**2.0E0 - 1.0E0)
        end do
      end do
      P2(n_t) = P2(n_t) / (n_chain * (n_mon*2-2))
      print *, n_t, P2(n_t)

    end do

    ! 5. Close the file
    call xtc % close

end program ordrprm
