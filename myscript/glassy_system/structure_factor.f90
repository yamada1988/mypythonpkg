!  XDR Fortran Interface XTC Example Program with Wrapper
!  2014 (c) James W. Barnett <jbarnet4@tulane.edu>
!  https://github.com/wesbarnett/
!

program example

    ! 1. Use the xdr interface
    use xdr, only: xtcfile


    ! 2. Declare a variable of type xtcfile
    type(xtcfile) :: xtcf
    real(kind=8), parameter :: pi = 3.141592653589793238462643383
    real(kind=8), parameter :: sgm = 0.60 !nm
    real(kind=8), parameter :: dk = 0.050
    real(kind=8), allocatable :: pos(:,:), S(:)
    real(kind=8) :: box_dim, xr(3), k, sinkr, diff
    integer, parameter :: Nk = 40
    integer :: npos, i, j, it, nt, ik

    allocate(S(Nk))
    
    ! 3. Initialize it with the names of xtc files you want to read in and write out
    call xtcf % init("../../../MD/md/md.xtc")
    npos = xtcf % NATOMS / 4  !number of water molecules (NATOMS is obtained after calling init)
    allocate(pos(3, npos))



    ! 4. Read in each configuration. Everything is stored in the xtcfile type (precision, time,
    !    step, no of atoms, positions, etc.). Look in the xtc module for more details.
    !    You can save the positions in the loop for your calculations in another array, or 
    !    do your calculations after each read.

    call xtcf % read
    it = 0
    nt = 0
    Sk = 0.0E0
    do while ( xtcf % STAT == 0 )

        if ( it > 20000 ) then
            exit
        end if
        if ( mod(it, 5000) /= 0) then 
        call xtcf % read
        it = it + 1
        else
        ! Just an example to show what was read in
        write(*,'(a,f12.6,a,i0)') " Time (ps): ", xtcf % time, "  Num. Atoms: ", xtcf % NATOMS
        pos = xtcf % pos(:, 1:xtcf % NATOMS:4)  !get the position of OW (every 1st atom)
        box_dim = xtcf % box(1,1)

        do ik = 1, Nk
          k = ik*dk
          S(ik) = S(ik) + 1.0 / npos
          do i = 1, npos
            do j = i+1, npos
              xr = pos(:, i) - pos(:, j)
              xr = xr - box_dim * nint(xr/box_dim)  !wrap distance (periodic boundary condition)
              diff = sqrt(sum(xr**2.0))
              sinkr = sin(k*diff)/(k*diff)
              S(ik) = S(ik) + 2.0 * sinkr/(npos**2.0)
              ! write (*,*) i, j, sinkr
            end do
          end do 
        end do

        call xtcf % read
        it = it + 1
        nt = nt + 1
        end if
    end do

    do ik = 1, Nk
      k = ik * dk
      S(ik) = S(ik) / nt
      write (*,*) k, S(ik)
    end do
   
    ! 5. Close the file
    call xtcf % close

end program example
