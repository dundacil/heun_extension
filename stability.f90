program stability
use mpi
implicit real*8 (a-h,o-z)

integer, parameter :: nd=150, nbin=2*nd+1, nnd=20, ndd=20, nalgo=17
! ialgo 1-nito      Ito schemes
! ialgo nstra-nalgo Stratonovich schemes
integer, parameter :: nito = 5, nstra = nito+1
integer :: ierr, rank, size, idum, i, jj, inner, ix, root
integer :: ninner
integer :: status(MPI_STATUS_SIZE)
real*8 :: dh, d1, x0, xlarge, dx, d0, dd, time, z0, w0
real*8 :: residui(nalgo),resiass(nalgo)
integer :: offset,nstable_local(nalgo),nstable_global(nalgo)
character*4 :: algo(nalgo)
character(len=120) :: fmt_string70
character(len=120) :: fmt_string71

external rnor

! Function definitions
!interface
!function rnor(idum) z
!integer :: idum
!real*8 :: rnor
!end function rnor
!end interface
! Code rewritten: ialgo 1-5  Ito schemes
!  ialgo 6-nalgo Stra

real*8 :: f0, f1, f2, g0, g1, g2
f0(x) = -x*(1.d0 + x*x)
f1(x) = -1.d0 - 3.d0*x*x
f2(x) = -6.d0*x
g0(x) = dd*(x*x + 1.d0)
g1(x) = dd*(2.d0*x)
g2(x) = dd*2.d0
g3(x) = 0.d0
frk(x) = f0(x) + 0.5d0*g0(x)*g1(x)

call MPI_Init(ierr)
call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)

root = 0
fmt_string70=""
write(fmt_string70, '(A,I0,A)') "(A,1x,A,1x,A,1x,", nalgo, "(A4,1x))"
fmt_string71=""
write(fmt_string71, '(A,I0,A)') "(A,1x,A,1x,", nalgo, "(A4,1x))"


print *,rank,size

do ialgo = 1,nalgo
   algo(ialgo) = "dumm"
   ide = int(ialgo/10)
   iun = mod(ialgo,10)
   write(algo(ialgo)(3:3),'(i1)') ide
   write(algo(ialgo)(4:4),'(i1)') iun
   if(rank == 0) then
      print *,ialgo, ide, iun
   endif
enddo
do idd = 1,ndd
do id = 1,nnd

!d0 = id*2.d-2
!dh = 11.d-2

d0 = idd*1.d-2
dh = 1.d-2*id
idum = -17 - 1000*rank   ! Different seed for each rank
tmax = 100.d0  !1000.d0
ninner = nint(tmax/dh)
nstable_local = 0
nstable_global = 0

d1 = dsqrt(dh)
dx = 1.d-2
dd = dsqrt(2.d0*d0)
xlarge = 100.d0

!!!!
ialgo=1   ! Ito 
algo(ialgo) = "Ito"
x0 = 0.d0
time = 0.d0
do i = 1, 1000
do inner = 1, ninner
do jj = 1, 10
time = time + dh
w0 = rnor(idum)
z0 = w0 * d1
x0 = x0 + dh * f0(x0) + g0(x0) * z0 
if (dabs(x0) .gt. xlarge) then
nstable_local(ialgo) = nstable_local(ialgo)+1
x0 = 0.d0
endif
enddo

enddo
enddo
!!!!

!!!!
ialgo=6  ! Stratonovich
algo(ialgo) = "Stra"
x0 = 0.d0
time = 0.d0
do i = 1, 1000
do inner = 1, ninner
do jj = 1, 10
time = time + dh
w0 = rnor(idum)
z0 = w0*d1  
x0 = x0+dh*f0(x0)+g0(x0)*z0+0.5d0*z0*z0*g0(x0)*g1(x0)
if (dabs(x0) .gt. xlarge) then
nstable_local(ialgo) = nstable_local(ialgo)+1
x0 = 0.d0
endif
enddo
enddo
enddo
!!!!

!!!!
ialgo=2  ! Heun to integrate a Ito dynamics using the correct first order ITO (not Euler)
algo(ialgo) = "HeIt"
x0 = 0.d0
time = 0.d0
do i = 1, 1000
do inner = 1, ninner
do jj = 1, 10
time = time + dh
w0 = rnor(idum)
z0 = w0*d1  
x1 = x0+dh*(f0(x0)-0.5d0*g0(x0)*g1(x0))+g0(x0)*z0
x2 = x0+dh*(f0(x1)-0.5d0*g0(x1)*g1(x1))+g0(x1)*z0
x0 = 0.5d0*(x1+x2)
if (dabs(x0) .gt. xlarge) then
nstable_local(ialgo) = nstable_local(ialgo)+1
x0 = 0.d0
endif
enddo
enddo
enddo
!!!!


!!!!
ialgo=4  ! Heun from Euler PECECCC for Ito
algo(ialgo) = "HeIC"
x0 = 0.d0
time = 0.d0
do i = 1, 1000
do inner = 1, ninner
do jj = 1, 10
time = time + dh
w0 = rnor(idum)
z0 = w0*d1  
x1 = x0+dh*(f0(x0)-0.5d0*g0(x0)*g1(x0))+g0(x0)*z0
x2 = x0+dh*(f0(x1)-0.5d0*g0(x1)*g1(x1))+g0(x1)*z0
do k = 1,3
x3 = 0.5d0*(x1+x2)
x2 = x0+dh*(f0(x3)-0.5d0*g0(x3)*g1(x3))+g0(x3)*z0
enddo
x0 = 0.5d0*(x1+x2)
if (dabs(x0) .gt. xlarge) then
nstable_local(ialgo) = nstable_local(ialgo)+1
x0 = 0.d0
endif
enddo
enddo
enddo
!!!!


!!!!
ialgo=7  ! Heun to integrate a Stra dynamics from Euler-Ito
algo(ialgo) = "Heun"
x0 = 0.d0
time = 0.d0
do i = 1, 1000
do inner = 1, ninner
do jj = 1, 10
time = time + dh
w0 = rnor(idum)
z0 = w0*d1  
x1 = x0+dh*f0(x0)+g0(x0)*z0
x2 = x0+dh*f0(x1)+g0(x1)*z0
x0 = 0.5d0*(x1+x2)
if (dabs(x0) .gt. xlarge) then
nstable_local(ialgo) = nstable_local(ialgo)+1
x0 = 0.d0
endif
enddo
enddo
enddo
!!!!

!!!!
ialgo=8  ! Heun to integrate a Stra dynamics from a Strat with - sign in the second step
algo(ialgo) = "HeSm"
x0 = 0.d0
time = 0.d0
do i = 1, 1000
do inner = 1, ninner
do jj = 1, 10
time = time + dh
w0 = rnor(idum)
z0 = w0*d1  
x1 = x0+dh*f0(x0)+g0(x0)*z0 + 0.5d0*z0*z0*g0(x0)*g1(x0)
x2 = x0+dh*f0(x1)+g0(x1)*z0 - 0.5d0*z0*z0*g0(x1)*g1(x1)
x0 = 0.5d0*(x1+x2)
if (dabs(x0) .gt. 
xlarge) then
nstable_local(ialgo) = nstable_local(ialgo)+1
x0 = 0.d0
endif
enddo
enddo
enddo
!!!!

!!!!
ialgo=17  ! Heun from a Strat with - sign in the second step PECECCC
algo(ialgo) = "HSmC"
x0 = 0.d0
time = 0.d0
do i = 1, 1000
do inner = 1, ninner
do jj = 1, 10
time = time + dh
w0 = rnor(idum)
z0 = w0*d1  
x1 = x0+dh*f0(x0)+g0(x0)*z0 + 0.5d0*z0*z0*g0(x0)*g1(x0)
x2 = x0+dh*f0(x1)+g0(x1)*z0 - 0.5d0*z0*z0*g0(x1)*g1(x1)
do k = 1,3
x3 = 0.5d0*(x1+x2)
x2 = x0+dh*f0(x3)+g0(x3)*z0 - 0.5d0*z0*z0*g0(x3)*g1(x3)
enddo
x0 = 0.5d0*(x1+x2)
if (dabs(x0) .gt. xlarge) then
nstable_local(ialgo) = nstable_local(ialgo)+1
x0 = 0.d0
endif
enddo
enddo
enddo
!!!!

!!!!
ialgo=9  ! Heun to integrate a Stra dynamics from Ito first order + Strat term
algo(ialgo) = "HeSp"
x0 = 0.d0
time = 0.d0
do i = 1, 1000
do inner = 1, ninner
do jj = 1, 10
time = time + dh
w0 = rnor(idum)
z0 = w0*d1  
x1 = x0+dh*(f0(x0)-0.5d0*g0(x0)*g1(x0))+g0(x0)*z0+0.5d0*z0*z0*g0(x0)*g1(x0)
x2 = x0+dh*(f0(x1)-0.5d0*g0(x1)*g1(x1))+g0(x1)*z0+0.5d0*z0*z0*g0(x1)*g1(x1)
x0 = 0.5d0*(x1+x2)
if (dabs(x0) .gt. xlarge) then
nstable_local(ialgo) = nstable_local(ialgo)+1
x0 = 0.d0
endif
enddo
enddo
enddo
!!!!

!!!!
ialgo=10  ! Heun to integrate a Stra dynamics from Strato secondo order (with a - sign on C)
algo(ialgo) = "HS2m"
x0 = 0.d0
time = 0.d0
do i = 1, 1000
do inner = 1, ninner
do jj = 1, 10
time = time + dh
w0 = rnor(idum)
z0 = w0*d1  
x1 = x0+dh*(f0(x0)-0.5d0*g0(x0)*g1(x0))+g0(x0)*z0+0.5d0*z0*z0*g0(x0)*g1(x0)
x2 = x0+dh*(f0(x1)+0.5d0*g0(x1)*g1(x1))+g0(x1)*z0-0.5d0*z0*z0*g0(x1)*g1(x1)
x0 = 0.5d0*(x1+x2)
if (dabs(x0) .gt. xlarge) then
nstable_local(ialgo) = nstable_local(ialgo)+1
x0 = 0.d0
endif
enddo
enddo
enddo
!!!!

!!!!
ialgo=11  ! Heun from Euler PECECCC
algo(ialgo) = "HeEC"
x0 = 0.d0
time = 0.d0
do i = 1, 1000
do inner = 1, ninner
do jj = 1, 10
time = time + dh
w0 = rnor(idum)
z0 = w0*d1  
x1 = x0+dh*f0(x0)+g0(x0)*z0
x2 = x0+dh*f0(x1)+g0(x1)*z0
do k = 1,3
x3 = 0.5d0*(x1+x2)
x2 = x0+dh*f0(x3)+g0(x3)*z0
enddo
x0 = 0.5d0*(x1+x2)
if (dabs(x0) .gt. xlarge) then
nstable_local(ialgo) = nstable_local(ialgo)+1
x0 = 0.d0
endif
enddo
enddo
enddo
!!!!


!!!!
ialgo=12  ! One step collocation to order h^3/2
algo(ialgo) = "T3/2"
x0 = 0.d0
time = 0.d0
zd2 = dsqrt(dh)/(dsqrt(3.d0))
do i = 1, 1000
do inner = 1, ninner
do jj = 1, 10
time = time + dh
w1 = rnor(idum)
w2 = rnor(idum)
z1 = w1*d1
z2 = dh*0.5d0*(z1+w2*zd2)
x0 = x0+dh*f0(x0)+g0(x0)*z1+0.5d0*z1*z1*g0(x0)*g1(x0)+&
   z2*(g0(x0)*f1(x0)-f0(x0)*g1(x0))+dh*z1*f0(x0)*g1(x0)+&
   z1*z1*z1*(g1(x0)*g1(x0)+g2(x0)*g0(x0))*g0(x0)/6.d0
if (dabs(x0) .gt. xlarge) then
nstable_local(ialgo) = nstable_local(ialgo)+1
x0 = 0.d0
endif
enddo
enddo
enddo
!!!!


!!!!
ialgo=13  ! One step collocation h^2 come da articolo con Vincenzo
algo(ialgo) = "Vinc"
x0 = 0.d0
time = 0.d0
zd2 = dsqrt(dh)/(dsqrt(3.d0))
do i = 1, 1000
do inner = 1, ninner
do jj = 1, 10
time = time + dh
w1 = rnor(idum)
w2 = rnor(idum)
w3 = rnor(idum)
z1 = w1*d1
z2 = dh*0.5d0*(z1+w2*zd2)
z3 = dh*(z1*z1+dh*(0.5d0+w3))/3.d0
g0x = g0(x0)
g1x = g1(x0)
g2x = g2(x0)
g3x = g3(x0)
f0x = f0(x0)
f1x = f1(x0)
f2x = f2(x0)
x0 = x0+dh*f0x+g0x*z1+0.5d0*z1*z1*g0x*g1x+&
   z2*(g0x*f1x-f0x*g1x)+dh*z1*f0x*g1x+&
   z1*z1*z1*(g1x*g1x+g2x*g0x)*g0x/6.d0 &
   +0.5d0*dh*dh*f1x*f0x+0.5d0*f1x*g1x*g0x*(z1*z2-z3) &
   +0.5d0*f2x*g0x*g0x*(z1*z2-z3) &
   +g1x*(g0x*f1x-f0x*g1x)*z3+0.5d0*(g1x*f0x*g1x+g2x*g0x*f0x*0.5d0)*(dh*z1*z1-z1*z2+z3) &
     + (g1x*(g1x*g1x*g0x+g2x*g0x*g0x)+g3x*g0x*g0x*g0x)/12.d0*z1*z1*z1*z1 
if (dabs(x0) .gt. xlarge) then
nstable_local(ialgo) = nstable_local(ialgo)+1
x0 = 0.d0
endif
enddo
enddo
enddo
!!!!

!!!!
ialgo=14  ! One step collocation h^2 come da CUP, ma corretto
algo(ialgo) = "Cupp"
x0 = 0.d0
time = 0.d0
zd2 = dsqrt(dh)/(dsqrt(3.d0))
zd3 = dh**2/6.d0
do i = 1, 1000
do inner = 1, ninner
do jj = 1, 10
time = time + dh
w1 = rnor(idum)
w2 = rnor(idum)
w3 = rnor(idum)
z1 = w1*d1
z2 = dh*0.5d0*(z1+w2*zd2)
z3 = dh*(z1*z1-dh)/6.d0+zd3*w3
g0x = g0(x0)
g1x = g1(x0)
g2x = g2(x0)
g3x = g3(x0)
f0x = f0(x0)
f1x = f1(x0)
f2x = f2(x0)
x0 = x0+dh*f0x+g0x*z1+0.5d0*z1*z1*g0x*g1x+&
   z2*(g0x*f1x-f0x*g1x)+dh*z1*f0x*g1x+&
   z1*z1*z1*(g1x*g1x+g2x*g0x)*g0x/6.d0 &
   +0.5d0*dh*dh*f1x*f0x+0.5d0*f1x*g1x*g0x*(z1*z2-z3) &
   +0.5d0*f2x*g0x*g0x*(z1*z2-z3) &
   +g1x*(g0x*f1x-f0x*g1x)*z3+0.5d0*(g1x*f0x*g1x+g2x*g0x*f0x*0.5d0)*(dh*z1*z1-z1*z2+z3) &
     + (g1x*(g1x*g1x*g0x+g2x*g0x*g0x)+g3x*g0x*g0x*g0x)/12.d0*z1*z1*z1*z1 
if (dabs(x0) .gt. xlarge) then
nstable_local(ialgo) = nstable_local(ialgo)+1
x0 = 0.d0
endif
enddo
enddo
enddo
!!!!

!!!!
ialgo=15  ! One step collocation calcolo gemini
algo(ialgo) = "Gemi"
x0 = 0.d0
time = 0.d0
zd2 = dsqrt(dh)/(dsqrt(3.d0))
zd3 = dsqrt(dh*dh*dh/6.d0)
do i = 1, 1000
do inner = 1, ninner
do jj = 1, 10
time = time + dh
w1 = rnor(idum)
w2 = rnor(idum)
w3 = rnor(idum)
z1 = w1*d1
z2 = dh*0.5d0*(z1+w2*zd2)
z3 = z1*z2 - 0.5d0*dh*dh - zd3*w3
g0x = g0(x0)
g1x = g1(x0)
g2x = g2(x0)
g3x = g3(x0)
f0x = f0(x0)
f1x = f1(x0)
f2x = f2(x0)
x0 = x0+dh*f0x+g0x*z1+0.5d0*z1*z1*g0x*g1x+&
   z2*(g0x*f1x-f0x*g1x)+dh*z1*f0x*g1x+&
   z1*z1*z1*(g1x*g1x+g2x*g0x)*g0x/6.d0 &
   +0.5d0*dh*dh*f1x*f0x+0.5d0*f1x*g1x*g0x*(z1*z2-z3) &
   +0.5d0*f2x*g0x*g0x*(z1*z2-z3) &
   +g1x*(g0x*f1x-f0x*g1x)*z3+0.5d0*(g1x*f0x*g1x+g2x*g0x*f0x*0.5d0)*(dh*z1*z1-z1*z2+z3) &
     + (g1x*(g1x*g1x*g0x+g2x*g0x*g0x)+g3x*g0x*g0x*g0x)/12.d0*z1*z1*z1*z1 
if (dabs(x0) .gt. xlarge) then
nstable_local(ialgo) = nstable_local(ialgo)+1
x0 = 0.d0
endif
enddo
enddo
enddo
!!!!

!!!!
ialgo=16  ! efficient RK rewritte for Strat
algo(ialgo) = "RKSt"
x0 = 0.d0
time = 0.d0
s3 = dsqrt(3.d0*dh)
zd2 = dsqrt(dh)/(dsqrt(3.d0))
do i = 1, 1000
do inner = 1, ninner
do jj = 1, 10
time = time + dh
w1 = rnor(idum)
z1 = w1*d1
f1x = frk(x0)
g1x = g0(x0)
x1 = x0 + dh*f1x+z1*g1x
f2x = frk(x1)
g2x = g0(x0-2.30/3*g1x*(z1+s3))
g3x = g0(x0+2.d0/9*g1x*(3.d0*z1+s3))
g4x = g0(x0-20.d0/27.d0*f1x*dh+10.d0/27.d0*(g2x-g1x)*z1-10.d0/27.d0*g2x*s3)
x0 = x0 +0.5d0*dh*(f1x+f2x)+z1/40.d0*(37.d0*g1x+30.d0*g3x-27.d0*g4x)+&
    s3/16.d0*(8*g1x+g2x-9*g3x)
   
if (dabs(x0) .gt. xlarge) then
nstable_local(ialgo) = nstable_local(ialgo)+1
x0 = 0.d0
endif
enddo
enddo
enddo

!!!!
ialgo=3  ! efficient RK Ito
algo(ialgo) = "RKIt"
x0 = 0.d0
time = 0.d0
s3 = dsqrt(3.d0*dh)
zd2 = dsqrt(dh)/(dsqrt(3.d0))
do i = 1, 1000
do inner = 1, ninner
do jj = 1, 10
time = time + dh
w1 = rnor(idum)
z1 = w1*d1
f1x = f0(x0)
g1x = g0(x0)
x1 = x0 + dh*f1x+z1*g1x
f2x = f0(x1)
g2x = g0(x0-2.30/3*g1x*(z1+s3))
g3x = g0(x0+2.d0/9*g1x*(3.d0*z1+s3))
g4x = g0(x0-20.d0/27.d0*f1x*dh+10.d0/27.d0*(g2x-g1x)*z1-10.d0/27.d0*g2x*s3)
x0 = x0 +0.5d0*dh*(f1x+f2x)+z1/40.d0*(37.d0*g1x+30.d0*g3x-27.d0*g4x)+&
    s3/16.d0*(8*g1x+g2x-9*g3x)
   
if (dabs(x0) .gt. xlarge) then
nstable_local(ialgo) = nstable_local(ialgo)+1
x0 = 0.d0
endif
enddo
enddo
enddo


!!!!
! Reduce all local rhos to global sum on root
call MPI_Reduce(nstable_local, nstable_global, nalgo, MPI_INTEGER, MPI_SUM, root, MPI_COMM_WORLD, ierr)


if (rank == root) then


if(id.eq.1.and.idd.eq.1) then
write(70,trim(fmt_string70)) "descriptor","dh","d0",(algo(ialgo),ialgo=1,nalgo)
print *,fmt_string
write(71,trim(fmt_string71)) "dh","d0",(algo(ialgo),ialgo=1,nalgo)
endif

write(70,*) dh,d0,(nstable_global(ialgo),ialgo=1,nalgo)
write(71,*) dh,d0,(nstable_global(ialgo),ialgo=1,nalgo)
 
end if


enddo
if (rank == root) then
write(71,*) " "
endif
enddo


call MPI_Finalize(ierr)
end program stability

