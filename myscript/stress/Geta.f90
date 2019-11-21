program calc_sigma
  implicit none
  integer i, i0, j, k
  double precision, allocatable,dimension(:,:) :: P
  double precision, allocatable,dimension(:) :: t, Geta
  integer n, nt
  double precision, allocatable,dimension(:,:) :: G
  double precision :: V, kBT, alpha
  V = Lbox**3.0E0 !(nm^3)
  kBT = 2.479 !(kJ/mol)
  alpha = kBT/V ! (kJ/mol/nm^3)
  alpha = alpha * 1.604*10**(6.0E0) !(Pa)
  alpha = alpha * 10**(-5.0E0) !(bar)

  n = Frames
  nt = 51
  allocate( P(3,n) )
  allocate( t(n) )
  allocate( G(3,0:nt) )
  allocate( Geta(0:nt) )
  G = 0.0

  open(99, file='ifname' )
  do j = 1,1
  read (99, *)
  end do


  do i = 1, n
    read (99,*) t(i), P(1,i), P(2,i), P(3,i)
  end do
  close(99)


open(95, file='ofname')

  do i = 0, nt-1
    Geta(i) = 0.0
    do k = 1,3
      do i0 = 1, n-i
        G(k,i) = G(k,i) + P(k,i+i0)*P(k,i0)
      end do
      G(k,i) = G(k,i)/(alpha*(n-i+1))*10**(5.0E0)
    end do
    do k=1,3
      Geta(i) = Geta(i) + G(k,i)
    end do
    Geta(i) = Geta(i)/3.0E0
    write (*,*) t(i+1), Geta(i)/Geta(0), Geta(i), G(1,i), G(2,i), G(3,i)
    write(95, "(f12.5, 2X, f13.8, E14.3, E14.3, E14.3, E14.3)") t(i+1),
Geta(i)/Geta(0), Geta(i), G(1,i), G(2,i), G(3,i)

  end do

  close(95)
end
