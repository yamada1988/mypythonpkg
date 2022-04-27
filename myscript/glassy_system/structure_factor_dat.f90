!  XDR Fortran Interface XTC Example Program with Wrapper
!  2014 (c) James W. Barnett <jbarnet4@tulane.edu>
!  https://github.com/wesbarnett/
!

program example
    !$ use omp_lib
    real(kind=8), parameter :: pi = 3.141592653589793238462643383
    real(kind=8), parameter :: sgm = 0.60 !nm
    real(kind=8), parameter :: dk = 0.50, dr = 0.000020, k0 = 0.000010
    real(kind=8), allocatable :: pos(:,:), S(:), S0(:), r(:), dqr(:), sinqr(:)
    real(kind=8) :: box_dim, xr(3), k, sinkr, diff, r_max, rho, r_
    integer, parameter :: Nk = 400
    integer :: i, j, it, nt, ik, nhist
    character(50) :: fname, dummy
    integer, parameter :: npos = 10000
    integer, parameter :: Ntmax = 1000, dt = 10


    allocate(S(Nk))
    allocate(S0(Nk))
    allocate(pos(3, npos))
    allocate(dqr(npos))
    allocate(sinqr(npos))
    S = 0.0E0
    S0 = 0.0E0 

    it = 10
    nt = 0

    write(fname, "('../../../DAT/drs/tau=01000/dr'i5.5'.dat')") it
    open(10,file=fname,status='old')
    read (10, *) dummy, box_dim
    close(10)
    it = 0
    r_max = box_dim / 2.0E0
    nhist = ceiling(r_max / dr)
    rho = npos / (box_dim**3.0E0)
    allocate(r(nhist))
    r = [((i - 0.50E0) * dr, i = 1, nhist)]  !set r-scales as the middle points
    


    do it=10,Ntmax,dt
        write (*,*) it
        write(fname, "('../../../DAT/drs/tau=01000/dr'i5.5'.dat')") it
        open(10,file=fname,status='old')
        read (10, *) dummy, box_dim
        do i = 1, npos
          read (10, *) pos(1, i), pos(2, i), pos(3, i), dqr(i)
          !if ( i == 1) then
          !  write(*,*) pos(:,i)
          !end if
        end do
        close(10)

        !$omp parallel private(xr, diff, sinkr)
        !$omp do 
        do ik = 1, Nk
          k = (ik-1)*dk + k0
          do i = 1, npos
            do j = i+1, npos
              xr = pos(:, i) - pos(:, j)
              xr = xr - box_dim * nint(xr/box_dim)  !wrap distance (periodic boundary condition)
              diff = sqrt(sum(xr**2.0))
              if ( diff <= r_max ) then
                sinkr = sin(k*diff)/(k*diff)
                !r_ = (ceiling(diff / dr) - 0.50E0) * dr
                !sinkr = sin(k*r_)/(k*r_)
                S(ik) = S(ik) + 2.0E0 * sinkr/npos
              end if
            end do
          end do 
        end do
        !$omp end do
        !$omp end parallel
    end do


    nt = int(Ntmax/dt) 
    S = S / real(nt)
    S = S + 1.0E0

    write (*,*) nt, nhist
    do ik=1, Nk
      k = (ik-1)*dk + k0
      S0(ik) = S0(ik) + 0.50E0*dr*rho*4.0*pi*r(1)**2.0E0*sin(k*r(1))/(k*r(1))
      S0(ik) = S0(ik) + 0.50E0*dr*rho*4.0*pi*r(nhist)**2.0E0*sin(k*r(nhist))/(k*r(nhist))
      do i=2, nhist-1
        r_ = r(i)
        S0(ik) = S0(ik) + dr*rho*4.0*pi*r_**2.0E0*sin(k*r_)/(k*r_)
      end do
    end do
 

    open (17,file='Sk_dat.dat')
    do ik = 1, Nk
      k = (ik-1)*dk + k0
      write (17, '(f10.6,2X,f10.4,2X,f10.4,2X,f10.4)') k, S(ik), S0(ik), S(ik)-S0(ik)
    end do
    close(17)
   

end program example
