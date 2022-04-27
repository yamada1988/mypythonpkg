!  XDR Fortran Interface XTC Example Program with Wrapper
!  2014 (c) James W. Barnett <jbarnet4@tulane.edu>
!  https://github.com/wesbarnett/
!

program example
    !$ use omp_lib
    use xdr, only: xtcfile
    implicit none
    type(xtcfile) :: xtc
    real(kind=8), parameter :: pi = 3.141592653589793238462643383
    real(kind=8), allocatable :: pos(:, :, :), diff(:, :)
    real(kind=8) :: box_dim, xr(3), dxr(3)
    integer :: i, j, it, jt, nt
    character(50) :: fname, dummy
    integer :: npos
    integer, parameter :: Ntmax = 100000, dt = 10, dtau = taunum

    call xtc % init("../../../MD/md_long/short.xtc")
    npos = xtc % NATOMS   !number of water molecules (NATOMS is obtained after calling init)
    allocate(pos(3, npos, int(Ntmax/dt)))
    allocate(diff(npos, int(Ntmax/dt)))
    call xtc % read

    ! box information cannot be obtained until at least one read call
    box_dim = xtc % box(1,1)

    it = 0
    nt = 0 
    jt = 0
    do it=1,Ntmax
        call xtc % read
        if ( mod(it, dt) /= 0 ) cycle
        pos(:, :, int(it/dt)) = xtc % pos(:, :)  
        if ( mod(it, dtau) /= 0 ) cycle
        write(*,*) it, int(it/dt)
    end do
    call xtc % close

     write(*,*) "drs calculation..."
    !$omp parallel private(xr, dxr, it, jt, diff)
    !$omp do 
    do i = 1, npos
      !write(*,*) i
      do it = 1, int(Ntmax/dt)-int(dtau/dt)
        dxr = 0.0E0
        do jt = 0, int(dtau/dt)-1
          !write(*,*) jt
          xr = pos(:, i, it+jt+1) - pos(:, i, it+jt)
          xr = xr - box_dim * nint(xr/box_dim)  !wrap distance (periodic boundary condition)
          dxr = dxr + xr
        end do
        diff(i, it) = sqrt(sum(dxr**2.0E0))
      end do
    end do
    !$omp end do
    !$omp end parallel
       
     !do it = 1, int(Ntmax/dt)-int(dtau/dt)
     do it = 1+int(dtau/dt), int(Ntmax/dt)
       write(fname, "('../../../DAT/drs/tau='i5.5'/dr'i5.5'.dat')") int(dtau), int(it*dt)
       open (17,file=fname)
       if ( mod(it, 10) == 0) then
         write(*,*) int(it*dt)
       end if
       write(17, '(A, 2X, f8.5)') , '#', box_dim
       do i = 1, npos
         write (17, '(f8.5,2X,f8.5,2X,f8.5,2X,f7.5)') pos(1, i, it), pos(2, i, it), pos(3, i, it), diff(i, it-(int(dtau/dt)))
       end do
     end do
    close(17)
   

end program example

