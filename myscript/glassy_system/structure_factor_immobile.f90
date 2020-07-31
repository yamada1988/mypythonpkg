!  XDR Fortran Interface XTC Example Program with Wrapper
!  2014 (c) James W. Barnett <jbarnet4@tulane.edu>
!  https://github.com/wesbarnett/
!

program example
    !$ use omp_lib
    real(kind=8), parameter :: pi = 3.141592653589793238462643383
    real(kind=8), parameter :: sgm = 0.20 !nm
    real(kind=8), parameter :: dr = 0.000020
    real(kind=8), allocatable :: pos(:,:), S(:), S0(:), r(:), dr_(:), sinqr(:)
    real(kind=8) :: box_dim, xr(3), k, sinumkr, diff, r_max, rho, r_, q, sumq, numt, sqij, numij, k0(4), dk(4)
    integer :: i, j, it, nt, ik, jk, nhist, Nk, numk(5), lk
    character(50) :: fname, dummy
    integer, parameter :: npos = 10000
    integer, parameter :: Ntmax = 100000, dt = 2000, tau=taunum

    ! set k 
    Nk = 0
    k0(1) = 0.00010E0
    k0(2) = 0.10E0
    k0(3) = 1.0E0
    k0(4) = 10.0E0

    dk(1) = 1.0E0*10.0E0**(-2)
    dk(2) = 2.50E0*10.0E0**(-2)
    dk(3) = 1.0E0*10.0E0**(-1)
    dk(4) = 1.0E0*1.0E0
    numk(0) = 0
    numk(1) = 9
    numk(2) = 36
    numk(3) = 89
    numk(4) = 90
    do i = 1, 4 
      Nk = Nk + numk(i)
    end do

    q = pi / sgm

    allocate(S(Nk))
    allocate(S0(Nk))
    allocate(pos(3, npos))
    allocate(dr_(npos))
    allocate(sinqr(npos))
    S = 0.0E0
    S0 = 0.0E0 

    it = 10
    nt = 0

    write(fname, "('../../../DAT/drs/tau='i5.5'/dr'i5.5'.dat')") tau, tau+it
    open(10,file=fname,status='old')
    read (10, *) dummy, box_dim
    close(10)
    it = 0
    r_max = box_dim / 2.0E0
    nhist = ceiling(r_max / dr)
    rho = npos / (box_dim**3.0E0)
    allocate(r(nhist))
    r = [((i - 0.50E0) * dr, i = 1, nhist)]  !set r-scales as the middle points
    

    numt = 0.0
    sumq = 0.0E0
    sqij = 0.0E0
    numij = 0.0E0
    do it=10+tau,Ntmax,dt
        write (*,*) it
        write(fname, "('../../../DAT/drs/tau='i5.5'/dr'i5.5'.dat')") tau, it
        open(10,file=fname,status='old')
        read (10, *) dummy, box_dim
        do i = 1, npos
          read (10, *) pos(1, i), pos(2, i), pos(3, i), dr_(i)
          !if ( i == 1) then
          !  write(*,*) pos(:,i)
          !end if
        sinqr(i) = sin(q*dr_(i))/(q*dr_(i))
        sumq = sumq + sinqr(i)/npos
        end do
        close(10)
        numt = numt + 1.0E0
        write(*,*) sumq / numt

        do ik = 1, 4
          !$omp parallel private(xr, diff, sinumkr, k)
          !$omp do reduction(+:sqij, numij) 
          do jk = 1, numk(ik)
            k = k0(ik) + dk(ik)*(jk-1)
            lk = sum(numk(:ik-1)) + jk
            do i = 1, npos
              do j = 1, npos
                xr = pos(:, i) - pos(:, j)
                xr = xr - box_dim * nint(xr/box_dim)  !wrap distance (periodic boundary condition)
                diff = sqrt(sum(xr**2.0))
                if ( diff <= r_max ) then
                  if ( i == j ) then
                    sinumkr = 1.0E0
                  else
                    sinumkr = sin(k*diff)/(k*diff)
                  end if
                  !r_ = (ceiling(diff / dr) - 0.50E0) * dr
                  !sinumkr = sin(k*r_)/(k*r_)
                  S(lk) = S(lk) + 1.0E0 * sinumkr/npos * sinqr(i) * sinqr(j)
                  sqij = sqij + sinqr(i) * sinqr(j)
                  numij = numij + 1.0E0
                end if
              end do
            end do 
          end do
         !$omp end do
         !$omp end parallel
        end do
        write(*,*) sqij/numij
    end do

    sumq = sumq / dble(numt)
    sqij = sqij / numij
    write(*,*) sumq, sumq*sumq, sqij


    S = S / dble(numt)

    write (*,*) numt, nhist
    do ik=1, 4
      do jk=1, numk(ik)
        k = k0(ik) + dk(ik)*(jk-1)
        lk = sum(numk(:ik-1)) + jk
        S0(lk) = S0(lk) + 0.50E0*dr*rho*4.0*pi*r(1)**2.0E0*sin(k*r(1))/(k*r(1))
        S0(lk) = S0(lk) + 0.50E0*dr*rho*4.0*pi*r(nhist)**2.0E0*sin(k*r(nhist))/(k*r(nhist))
        do i=2, nhist-1
          r_ = r(i)
          S0(lk) = S0(lk) + dr*rho*4.0*pi*r_**2.0E0*sin(k*r_)/(k*r_)
        end do
      end do
    end do

    !S0 = S0 * sumq * sumq
    !S = S / (sumq*sumq)


    write(fname, "('./Sk'i5.5'_'f5.3'.dat')") tau, sgm
    open (17,file=fname)
    do ik = 1, 4
      do jk=1, numk(ik)
        k = k0(ik) + dk(ik)*(jk-1)
        lk = sum(numk(:ik-1)) + jk
        write(*,*) lk, k
        write (17, '(f10.6,2X,f10.4,2X,f10.4,2X,f10.4,2X,f10.4)') k, S(lk), S0(lk), S(lk)/(sumq*sumq)-S0(lk), S(lk)/(sqij)-S0(lk)
      end do
    end do
    close(17)
   

end program example
