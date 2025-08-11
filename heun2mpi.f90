program mpi_rho_sum
use mpi
implicit real*8 (a-h,o-z)

integer, parameter :: nd=150, nbin=2*nd+1, nnd=20, nalgo=17
! ialgo 1-nito      Ito schemes
! ialgo nstra-nalgo Stratonovich schemes
integer, parameter :: nito = 5, nstra = nito+1
integer :: ierr, rank, size, idum, i, jj, inner, ix, root
integer :: ninner
integer :: status(MPI_STATUS_SIZE)
real*8 :: rho_local(nbin,nalgo), rho_global(nbin,nalgo)
real*8 :: dh, d1, x0, xlarge, dx, d0, dd, time, z0, w0
real*8 :: residui(nalgo),resiass(nalgo)
integer :: offset, nstable(nalgo)
character*4 :: algo(nalgo)
character(len=80) :: fmt_string

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

print *,rank,size

do ialgo = 1,nalgo
   algo(ialgo) = "dumm"
   ide = int(ialgo/10)
   iun = mod(ialgo,10)
   write(algo(3:3),'(i1)') ide
   write(algo(4:4),'(i1)') iun  
enddo

do id = 1,nnd
rho_local = 0.d0
rho_global = 0.d0
!d0 = id*2.d-2
!dh = 11.d-2

d0 = 21.d-2
dh = 1.d-2*id
idum = -17 - 1000*rank   ! Different seed for each rank
ninner = 10000


d1 = dsqrt(dh)
dx = 1.d-2
dd = dsqrt(2.d0*d0)
xlarge = 10.d0

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
if (dabs(x0) .gt. xlarge) x0 = 0.d0
enddo
ix = nint(x0 / dx)
if (iabs(ix) .le. nd) rho_local(ix + nd + 1,ialgo) = rho_local(ix + nd + 1,ialgo) + 1
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
if (dabs(x0) .gt. xlarge) x0 = 0.d0
enddo
ix = nint(x0 / dx)
if (iabs(ix) .le. nd) rho_local(ix + nd + 1,ialgo) = rho_local(ix + nd + 1,ialgo) + 1
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
if (dabs(x0) .gt. xlarge) x0 = 0.d0
enddo
ix = nint(x0 / dx)
if (iabs(ix) .le. nd) rho_local(ix + nd + 1,ialgo) = rho_local(ix + nd + 1,ialgo) + 1
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
if (dabs(x0) .gt. xlarge) x0 = 0.d0
enddo
ix = nint(x0 / dx)
if (iabs(ix) .le. nd) rho_local(ix + nd + 1,ialgo) = rho_local(ix + nd + 1,ialgo) + 1
enddo
enddo
!!!!

!!!!
ialgo=8  ! Heun to integrate a Stra dynamics from a Strat with - sign in the second step
algo(ialgo) = "HeS-"
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
if (dabs(x0) .gt. xlarge) x0 = 0.d0
enddo
ix = nint(x0 / dx)
if (iabs(ix) .le. nd) rho_local(ix + nd + 1,ialgo) = rho_local(ix + nd + 1,ialgo) + 1
enddo
enddo
!!!!

!!!!
ialgo=17  ! Heun from a Strat with - sign in the second step PECECCC
algo(ialgo) = "HS-C"
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
if (dabs(x0) .gt. xlarge) x0 = 0.d0
enddo
ix = nint(x0 / dx)
if (iabs(ix) .le. nd) rho_local(ix + nd + 1,ialgo) = rho_local(ix + nd + 1,ialgo) + 1
enddo
enddo
!!!!

!!!!
ialgo=9  ! Heun to integrate a Stra dynamics from Ito first order + Strat term
algo(ialgo) = "HeS+"
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
if (dabs(x0) .gt. xlarge) x0 = 0.d0
enddo
ix = nint(x0 / dx)
if (iabs(ix) .le. nd) rho_local(ix + nd + 1,ialgo) = rho_local(ix + nd + 1,ialgo) + 1
enddo
enddo
!!!!

!!!!
ialgo=10  ! Heun to integrate a Stra dynamics from Strato secondo order (with a - sign on C)
algo(ialgo) = "HS2-"
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
if (dabs(x0) .gt. xlarge) x0 = 0.d0
enddo
ix = nint(x0 / dx)
if (iabs(ix) .le. nd) rho_local(ix + nd + 1,ialgo) = rho_local(ix + nd + 1,ialgo) + 1
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
if (dabs(x0) .gt. xlarge) x0 = 0.d0
enddo
ix = nint(x0 / dx)
if (iabs(ix) .le. nd) rho_local(ix + nd + 1,ialgo) = rho_local(ix + nd + 1,ialgo) + 1
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
if (dabs(x0) .gt. xlarge) x0 = 0.d0
enddo
ix = nint(x0 / dx)
if (iabs(ix) .le. nd) rho_local(ix + nd + 1,ialgo) = rho_local(ix + nd + 1,ialgo) + 1
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
if (dabs(x0) .gt. xlarge) x0 = 0.d0
enddo
ix = nint(x0 / dx)
if (iabs(ix) .le. nd) rho_local(ix + nd + 1,ialgo) = rho_local(ix + nd + 1,ialgo) + 1
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
if (dabs(x0) .gt. xlarge) x0 = 0.d0
enddo
ix = nint(x0 / dx)
if (iabs(ix) .le. nd) rho_local(ix + nd + 1,ialgo) = rho_local(ix + nd + 1,ialgo) + 1
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
if (dabs(x0) .gt. xlarge) x0 = 0.d0
enddo
ix = nint(x0 / dx)
if (iabs(ix) .le. nd) rho_local(ix + nd + 1,ialgo) = rho_local(ix + nd + 1,ialgo) + 1
enddo
enddo
!!!!

!!!!
ialgo=16  ! efficient RK rewritte for Ito
algo(ialgo) = "RKef"
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
   
if (dabs(x0) .gt. xlarge) x0 = 0.d0
enddo
ix = nint(x0 / dx)
if (iabs(ix) .le. nd) rho_local(ix + nd + 1,ialgo) = rho_local(ix + nd + 1,ialgo) + 1
enddo
enddo

!!!!
ialgo=3  ! efficient RK 
algo(ialgo) = "RKef"
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
   
if (dabs(x0) .gt. xlarge) x0 = 0.d0
enddo
ix = nint(x0 / dx)
if (iabs(ix) .le. nd) rho_local(ix + nd + 1,ialgo) = rho_local(ix + nd + 1,ialgo) + 1
enddo
enddo


!!!!
! Reduce all local rhos to global sum on root
call MPI_Reduce(rho_local, rho_global, nbin*nalgo, MPI_DOUBLE_PRECISION, MPI_SUM, root, MPI_COMM_WORLD, ierr)


if (rank == root) then
rho_global=rho_global/size/ninner/1000/dx
pi = 4.d0*datan(1.d0)
residui = 0.d0
resiass = 0.d0
aa = 1.d0/(2.d0*d0)

! Ito like
anorm = gamma(2.d0+aa)/gamma(1.5d0+aa)/dsqrt(pi)
print *,anorm

fmt_string=""
write(fmt_string, '(A,I0,A)') "(A,1x,A,1x,", nito, "(A4,1x))"
write(200+id,fmt_string) "x ","theo ",(algo(ialgo),ialgo=1,nito)
fmt_string=""
write(fmt_string, '(A,I0,A)') "(A,1x,", nito, "(A4,1x))"
write(400+id,fmt_string) "x ",(algo(ialgo),ialgo=1,nito)

do i = -nd, nd
x = i*dx
theory =  anorm/((1+x**2)**(2.d0+aa))
!write(200+id,*) x, theory,rho_global(i+nd+1,1),&
!rho_global(i+nd+1,2)
!write(400+id,*) x,theory-rho_global(i+nd+1,1),&
!theory-rho_global(i+nd+1,2)
!do ialgo = 1,2
!residui(ialgo) = residui(ialgo) + (theory-rho_global(i+nd+1,ialgo))**2
!resiass(ialgo) = resiass(ialgo) + (1.d0-rho_global(i+nd+1,ialgo)/theory)**2
!enddo
write(200+id,*) x, theory,(rho_global(i+nd+1,ialgo),ialgo=1,nito)
write(400+id,*) x, (theory-rho_global(i+nd+1,ialgo),ialgo=1,nito)
do ialgo = 1,nito
residui(ialgo) = residui(ialgo) + (theory-rho_global(i+nd+1,ialgo))**2
resiass(ialgo) = resiass(ialgo) + (1.d0-rho_global(i+nd+1,ialgo)/theory)**2
enddo

enddo


! Strato like
anorm = gamma(1.d0+aa)/gamma(0.5d0+aa)/dsqrt(pi)
print *,anorm,gamma(2.d0)
fmt_string=""
write(fmt_string, '(A,I0,A)') "(A,1x,A,1x,", nalgo-nstra+1, "(A4,1x))"
write(300+id,fmt_string) "x ","theo ",(algo(ialgo),ialgo=nstra,nalgo)
fmt_string=""
write(fmt_string, '(A,I0,A)') "(A,1x,", nalgo-nstra+1, "(A4,1x))"
write(500+id,fmt_string) "x ",(algo(ialgo),ialgo=nstra,nalgo)

!write(300+id,'(A,A,' // char(ichar('0') + nalgo-nstra+1) // 'A)')  "x ","theo ",(algo(ialgo),ialgo=nstra,nalgo)
!write(500+id,'(A,' // char(ichar('0') + nalgo-nstra+1) // 'A)') "x ",(algo(ialgo),ialgo=nstra,nalgo)
do i = -nd, nd
x = i*dx
theory = anorm/((1+x**2)**(1.d0+aa))
!write(300+id,*) x, theory,rho_global(i+nd+1,3),&
!rho_global(i+nd+1,4),rho_global(i+nd+1,5),rho_global(i+nd+1,6),rho_global(i+nd+1,7),&
!rho_global(i+nd+1,8),rho_global(i+nd+1,9),rho_global(i+nd+1,10),rho_global(i+nd+1,11)
!write(500+id,*) x, theory-rho_global(i+nd+1,3),&
!theory-rho_global(i+nd+1,4),theory-rho_global(i+nd+1,5),&
!theory-rho_global(i+nd+1,6),theory-rho_global(i+nd+1,7),&
!theory-rho_global(i+nd+1,8),theory-rho_global(i+nd+1,9),&
!theory-rho_global(i+nd+1,10),theory-rho_global(i+nd+1,10)
write(300+id,*) x, theory,(rho_global(i+nd+1,ialgo),ialgo=nstra,nalgo)
write(500+id,*) x, (theory-rho_global(i+nd+1,ialgo),ialgo=nstra,nalgo)
do ialgo = nstra,nalgo
residui(ialgo) = residui(ialgo) + (theory-rho_global(i+nd+1,ialgo))**2
resiass(ialgo) = resiass(ialgo) + (1.d0-rho_global(i+nd+1,ialgo)/theory)**2
enddo

enddo
if(id.eq.1) then
fmt_string=""
write(fmt_string, '(A,I0,A)') "(A,1x,A,1x,", nalgo, "(A4,1x))"
write(80,fmt_string) "dh","d0",(algo(ialgo),ialgo=1,nalgo)
write(81,fmt_string) "dh","d0",(algo(ialgo),ialgo=1,nalgo)
endif
write(80,*) dh,d0,(residui(ialgo)/(2*nd+1),ialgo=1,nalgo)
write(81,*) dh,d0,(resiass(ialgo)/(2*nd+1),ialgo=1,nalgo)

 CLOSE(200+id)
 CLOSE(400+id)
 
 CLOSE(300+id)
 CLOSE(500+id)
 
end if


enddo



call MPI_Finalize(ierr)
end program mpi_rho_sum
!
!Notes:
!
!
!
! set key autotitle columnhead
!
!

!* This assumes rnor(idum) is a user-defined function generating normal-distributed random numbers using idum as the seed (e.g. Box–Muller method or similar).
!* We re-index rho from ρ(1, -nd\:nd) to a 1D rho\_local(nbin) with index shift: ix + nd + 1.
!* You can replace the print loop with file I/O if needed.

!To compile:

!mpif90 -O2 -o mpi\_rho mpi\_rho.f90

!To run:

!mpirun -np 4 ./mpi\_rho

!Let me know if you need help adapting rnor, adding MPI I/O, or improving performance.

