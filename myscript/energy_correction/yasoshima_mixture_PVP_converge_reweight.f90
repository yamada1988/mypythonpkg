      program main
      implicit none
      save

!input data : set_num_soln; number of samples for soln
!input data : set_num_refs ; number of samples for refs
!           : set_line; file lines(weight_refs,engref) 

integer :: set_num_soln=ks
integer :: set_num_refs=kr
print *, "===== soln ====="
print *, "weight start"
call ave_weight_sln(set_num_soln)

contains
!=======================================================================
subroutine ave_weight_sln(set_num_soln)
implicit none
save

integer :: k,l,i
integer,intent(in) :: set_num_soln
integer :: set_line=5
integer, dimension(:),allocatable :: num
character(len=80) :: fname, fnamei
real*8,dimension(:),allocatable :: ave_w,sum_w,ave_w2, sum_w2, wU
real*8,dimension(:,:),allocatable :: w,w2
integer :: sample_num = 0

allocate(num(set_line),ave_w(set_line),sum_w(set_line),ave_w2(set_line),sum_w2(set_line))
allocate(wU(0:set_num_soln))
allocate(w(0:set_num_soln,set_line),w2(0:set_num_soln,set_line))

 ave_w = 0
 sum_w = 0
 ave_w2 = 0
 sum_w2 = 0
 wU = 0

 do l=0,set_num_soln-1

  write(fname,"('ERmod_',i4.4,'/soln/weight_soln')") l+1
  write(fnamei,"('ERmod_',i4.4,'/soln/ave_SysWght')") l+1
  open(21,file=fnamei)
  read(21,*) wU(l)
  close(21)

  if (access( fname, " ") /= 0) cycle
    sample_num = sample_num + 1
  open(10,file=fname)
  do i=1,set_line
   read(10,1200) num(i),w(l,i)
   w(l,i) = w(l,i)*wU(l)
   sum_w(i) = sum_w(i) + w(l,i) 
  end do
  close(10)
 end do

 open(20,file="aveERmod_PVP_conv/soln/weight_soln")
 do i=1,set_line
  ave_w(i) = sum_w(i)/dble(sample_num+1)
  write(20,1200) num(i),ave_w(i)
 end do
 close(20)
1200 format(I5,E20.8)

print *, "aveuv start"
call ave_uv_sln(set_num_soln,ave_w,w,sum_w, ave_w2,w2,sum_w2)

return
end subroutine ave_weight_sln
!=======================================================================

!=======================================================================
subroutine ave_uv_sln(set_num_soln,ave_w,w,sum_w,ave_w2,w2,sum_w2)
implicit none
save

integer :: k,l,i,j
integer,intent(in) :: set_num_soln
real*8,dimension(:),intent(in) :: ave_w,sum_w,ave_w2,sum_w2
real*8,dimension(0:,:),intent(in) :: w,w2
integer :: set_line=5
character(len=80) :: fname,fname2
character(len=25),dimension(:),allocatable :: line
real*8,dimension(:),allocatable :: ave_e, ave_e2
real*8,dimension(:,:),allocatable :: e_num, e_num2

allocate(ave_e(set_line), ave_e2(set_line))
allocate(e_num(0:set_num_soln,set_line), e_num2(0:set_num_soln,set_line))

  ave_e = 0.d0
  e_num = 0.0
  ave_e2 = 0.d0
  e_num2 = 0.0

  do l=0,set_num_soln-1
   write(fname,"('ERmod_',i4.4,'/soln/aveuv.tt')") l+1
   if (access( fname, " ") /= 0) cycle
   open(10,file=fname)
   do i=1,set_line
    read(10, *) j,e_num(l,i), e_num2(l,i)
    !print *, l+1, i, e_num(l,i), e_num2(l,i)
    ave_e(i) = ave_e(i) + e_num(l,i)*w(l,i)
    ave_e2(i) = ave_e2(i) + e_num2(l,i)*w(l,i)
   end do
   close(10)
  end do
  
  print *, ave_e
  print *, ave_e2
  write(fname2,"('aveERmod_PVP_conv/soln/aveuv.tt')") 
  open(20,file=fname2)
  do i=1,set_line
   ave_e(i) = ave_e(i)/sum_w(i)
   ave_e2(i) = ave_e2(i)/sum_w(i)
   !print *, i, ave_e(i), ave_e2(i)
   write(20,*) i, ave_e(i), ave_e2(i)
  end do
  close(20)
700 format(A12,f7.5, f7.5)

print *, "eng start"
call ave_eng_sln(set_num_soln,ave_w,w,sum_w, ave_w2,w2,sum_w2)

return
end subroutine ave_uv_sln
!=======================================================================


!=======================================================================
subroutine ave_eng_sln(set_num_soln,ave_w,w,sum_w, ave_w2,w2,sum_w2)
implicit none
save

integer :: k,l,i
integer,intent(in) :: set_num_soln
real*8,dimension(:),intent(in) :: ave_w,sum_w,ave_w2,sum_w2
real*8,dimension(0:,:),intent(in) :: w,w2
integer :: set_line=2596
character(len=80) :: fname,fname2,fname3
character(len=25),dimension(:),allocatable :: line
real*8,dimension(:),allocatable :: ave_e, ave_e2
real*8,dimension(:,:),allocatable :: e_num, e_num2

allocate(line(set_line),ave_e(set_line), ave_e2(set_line))
allocate(e_num(0:set_num_soln,set_line), e_num2(0:set_num_soln,set_line))

 do k=1,5
  ave_e = 0.d0
  e_num = 0.0

  do l=0,set_num_soln-1
   write(fname3,"('ERmod_',i4.4,'/soln/aveuv.tt')") l+1
   write(fname,"('ERmod_',i4.4,'/soln/engsln.',i2.2)") l+1,k
  if (access( fname3, " ") /= 0) cycle
   open(10,file=fname)
   do i=1,set_line
    read(10,660) line(i), e_num(l,i)
    ave_e(i) = ave_e(i) + e_num(l,i)*w(l,k)
   end do
   close(10)
  end do

  write(fname2,"('aveERmod_PVP_conv/soln/engsln.'i2.2)") k
  open(20,file=fname2)
  do i=1,set_line
   ave_e(i) = ave_e(i)/sum_w(k)
   write(20,660) line(i), ave_e(i)
  end do

  close(20)
 end do
600 format(A20,E25.15e2)
660 format(A20,g25.15)

print *, "===== refs ====="
print *, "weight start"
call ave_weight(set_num_refs)

return
end subroutine ave_eng_sln
!=======================================================================


!=======================================================================
subroutine ave_weight(set_num_refs)
implicit none
save

integer :: k,l,i
integer,intent(in) :: set_num_refs
integer :: set_line=5
integer, dimension(:),allocatable :: num
character(len=80) :: fname,fnamei
real*8,dimension(:),allocatable :: ave_w,sum_w,ave_s,wU
real*8,dimension(:,:),allocatable :: w,s
integer :: sample_num = 0

allocate(num(set_line),ave_w(set_line),sum_w(set_line),ave_s(set_line))
allocate(wU(0:set_num_refs))
allocate(w(0:set_num_refs,set_line),s(0:set_num_refs,set_line))

 ave_w = 0
 sum_w = 0
 ave_s = 0
 wU = 0

 do l=0,set_num_refs-1

  write(fname,"('ERmod_',i4.4,'/refs/weight_soln')") l+1
  write(fnamei,"('ERmod_',i4.4,'/refs/ave_SysWght')") l+1
  open(21,file=fnamei)
  read(21,*) wU(l)
  close(21)

  if (access( fname, " ") /= 0) cycle
  sample_num = sample_num+1
  open(10,file=fname)
  do i=1,set_line
   read(10,1000) num(i),w(l,i),s(l,i)
   w(l,i) = w(l,i)*wU(l)
   s(l,i) = s(l,i)*wU(l)
   sum_w(i) = sum_w(i) + w(l,i)
   ave_s(i) = ave_s(i) + s(l,i) * w(l,i)
  end do
  close(10)
 end do

 open(20,file="aveERmod_PVP_conv/refs/weight_refs")
 do i=1,set_line
  ave_w(i) = sum_w(i)/dble(sample_num+1)
  ave_s(i) = ave_s(i)/sum_w(i)
  write(20,1000) num(i),ave_w(i),ave_s(i)
 end do
 close(20)
1000 format(I5,E20.8,E20.7)

print *, "eng start"
call ave_eng(set_num_refs,ave_w,w,sum_w)
print *, "all jobs end"

return
end subroutine ave_weight
!=======================================================================


!=======================================================================
subroutine ave_eng(set_num_refs,ave_w,w,sum_w)
implicit none
save

integer :: k,l,i
integer,intent(in) :: set_num_refs
real*8,dimension(:),intent(in) :: ave_w,sum_w
real*8,dimension(0:,:),intent(in) :: w
integer :: set_line=2596
character(len=80) :: fname,fname2,fname3
character(len=25),dimension(:),allocatable :: line
real*8,dimension(:),allocatable :: ave_e
real*8,dimension(:,:),allocatable :: e_num

allocate(line(set_line),ave_e(set_line))
allocate(e_num(0:set_num_refs,set_line))

 do k=1,5
  ave_e = 0.d0
  e_num = 0.0

  do l=0,set_num_refs-1
   write(fname3,"('ERmod_',i4.4,'/refs/aveuv.tt')") l+1
   write(fname,"('ERmod_',i4.4,'/refs/engsln.',i2.2)") l+1,k
   if (access( fname3, " ") /= 0) cycle
   open(10,file=fname)
   do i=1,set_line
    read(10,660) line(i), e_num(l,i)
    ave_e(i) = ave_e(i) + e_num(l,i)*w(l,k)
   end do
   close(10)
  end do

  write(fname2,"('aveERmod_PVP_conv/refs/engref.'i2.2)") k
  open(20,file=fname2)
  do i=1,set_line
   ave_e(i) = ave_e(i)/sum_w(k)
   write(20,660) line(i), ave_e(i)
  end do
  close(20)
 end do
600 format(A20,E25.15e2)
660 format(A20,g25.15)

print *, "corref start"
call ave_cor(set_num_refs,ave_w,w,sum_w,set_line)

return
end subroutine ave_eng
!=======================================================================

!=======================================================================
subroutine ave_cor(set_num_refs,ave_w,w,sum_w,set_line)
implicit none
save

integer :: k,l,i,j
integer,intent(in) :: set_num_refs, set_line
real*8,dimension(:),intent(in) :: ave_w,sum_w
real*8,dimension(0:,:),intent(in) :: w
real*8,dimension(:,:), allocatable :: cormat_temp,rdcor
character(len=80) :: fname,fname2,fname3

allocate(cormat_temp(set_line,set_line),rdcor(set_line,set_line))

 do k=1,5
  cormat_temp = 0
  rdcor = 0

  do l=0,set_num_refs-1
   write(fname,"('ERmod_',i4.4,'/refs/corsln.',i2.2)") l+1,k
   write(fname3,"('ERmod_',i4.4,'/refs/aveuv.tt')") l+1
   if (access( fname3, " ") /= 0) cycle
   open(10,file=fname,status='old',form="UNFORMATTED")
   read(10) cormat_temp
   close(10)

   do i=1,set_line
    do j=1,set_line
     rdcor(i,j) = rdcor(i,j) + cormat_temp(i,j)*w(l,k)
    end do
   end do
  end do

  write(fname2,"('aveERmod_PVP_conv/refs/corref.'i2.2)") k
  open(20,file=fname2,form="UNFORMATTED")
  do i=1,set_line
   do j=1,set_line
    rdcor(i,j) = rdcor(i,j)/sum_w(k)
   end do
  end do
  write(20) rdcor
  close(20)
 end do

return
end subroutine ave_cor

!=======================================================================
end program main
