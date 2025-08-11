program stability
use mpi
implicit real*8 (a-h,o-z)

integer, parameter :: nd=150, nbin=2*nd+1, nnd=20, ndd=20, nalgo=17
! ialgo 1-nito      Ito schemes
! ialgo nstra-nalgo Stratonovich schemes
integer, parameter :: nito = 5, nstra = nito+1, ndh=13,  npoichk=16, nnoise=npoichk*2**ndh
integer :: ierr, rank, size, idum, i, jj, inner, ix, root
integer :: ninner
integer :: status(MPI_STATUS_SIZE)
real*8 :: dh, d1, x0, xlarge, dx, d0, dd, time, z0, w0
real*8 :: residui(nalgo),resiass(nalgo),res_loc(npoichk,ndh),res_glo(npoichk,ndh)
integer :: offset,nstable_local(nalgo),nstable_global(nalgo)
character*4 :: algo(nalgo)
real*8 :: noise(nnoise),xsave(npoichk,0:ndh),xexact(0:nnoise),noise2(nnoise)
character(len=120) :: fmt_string70
character(len=120) :: fmt_string71
logical :: got

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
!fmt_string70=""
!write(fmt_string70, '(A,I0,A)') "(A,1x,A,1x,A,1x,", nalgo, "(A4,1x))"
!fmt_string71=""
!write(fmt_string71, '(A,I0,A)') "(A,1x,A,1x,", nalgo, "(A4,1x))"


print *,rank,size

nave = 100000
idum = -17 - 1000*rank   ! Different seed for each rank
idum2 = -37 - 1000*rank   ! Different seed for each rank
dh0 = 1.d0/npoichk
ndh0 = 1 !  ndh-7
a=1.5d0 ! 2.d0
pi = 4.d0*datan(1.d0)

do ialgo = 1,nalgo
   algo(ialgo) = "dumm"
   ide = int(ialgo/10)
   iun = mod(ialgo,10)
   write(algo(ialgo)(3:3),'(i1)') ide
   write(algo(ialgo)(4:4),'(i1)') iun
   if(rank == 0) then
!      print *,ialgo, ide, iun
   endif
enddo

! Nella simulazione, dobbiamo anzitutto decidere quanti punti massimi faremo
! diciamo 16*2^ndh punti massimo, con ndh=10
! il passo di integrazione elementare per ndh=0 sara` 0.1d0
! quindi stiamo seguendo cose per un tempo di 1.6
! per ndh=10 il passo di integrazione sara` 0.1/2^10
! Serve un vettore del rumore pari a 2^(ndh+4) elementi
! Facciamo quindi un loop sulla media iave: calcoliamo il vettore del rumore e poi procediamo
! Potendo fare diverse intensita` del rumore, prima mettiamo intensita` del rumore

do id = 1,1 !5   ! ipotizziamo 5 intensita` del rumore
 d0 = 5.d-2 ! 1.d0 ! 2.d-2 * id   ! rumore intensita` d0
 dd = dsqrt(2*d0) ! 0.5d0 ! dsqrt(2.d0*d0)

 res_glo = 0.d0
 res_loc = 0.d0
 do iave=1, nave

! First generate the noise
  do inoise=1,nnoise
   noise(inoise) = rnor(idum)
   noise2(inoise) = rnor(idum2)
  enddo
! We have the noise, now we integrate

! Now we generate the initial point according to the Stratonovich equilibrium
  aa = 1.d0/(2.d0*d0)
  got = .false.

  do while (.not.got)
    xs = dtan(pi*(ran2(idum)-0.5d0)) ! Distributed according to a Lorentz distribution
    ratio = 1.d0/(1+xs**2)**aa ! rejection method, Lorentzian is larger than the one we need
    if(ratio.gt.ran2(idum)) got=.true.
  enddo
!  xs = 2.d0*ran2(idum) -1.d0
!  xexact(0)=xs
  xsave = 0.d0

!  theory = 1.d0/((1+x**2)**(1.d0+aa))
! we now pick the integration time step
! compute exact
  dhm = dh0/2**ndh
  ddgm = dd*dsqrt(dhm)
  d1 = dsqrt(dhm)
  zd2 = dsqrt(dhm)/(dsqrt(3.d0))
!  do ipoint=1,nnoise
!   xexact(ipoint) = xexact(ipoint-1)*dexp((a-0.5d0*dd*dd)*dhm+ddgm*noise(ipoint))
!  enddo
  do idh = ndh0,ndh
   dh = dh0/2**idh
!   d1 = dsqrt(dh)
! For each dh, we need to identify the ninner  and the nskip
   ninner = 2**idh
   nskip = 2**(ndh-idh)
   time = 0.d0
   wtot = 0.d0
   x0 = xs
   do ipoints = 1,npoichk  ! these are the 16 points which we use to evaluate things
! equivalent noise
      do inner=1,ninner
         w0 = 0.d0
         w1 = 0.d0
         do iw = 1,nskip
            w0 = w0 + noise(iw+(ipoints-1)*ninner*nskip+(inner-1)*nskip)
            w1 = w1 + noise2(iw+(ipoints-1)*ninner*nskip+(inner-1)*nskip)
         enddo
         wtot = wtot + w0
         z1 = w0 * d1
         z2 = dh*0.5d0*(z1+w1*zd2)
         x0 = x0+dh*f0(x0)+g0(x0)*z1+0.5d0*z1*z1*g0(x0)*g1(x0)+&
            z2*(g0(x0)*f1(x0)-f0(x0)*g1(x0))+dh*z1*f0(x0)*g1(x0)+&
            z1*z1*z1*(g1(x0)*g1(x0)+g2(x0)*g0(x0))*g0(x0)/6.d0
         time = time + dh
!         x2 = x0 + dh*f0(x1) + g0(x1) * z0
!         x0 = 0.5d0*(x1+x2)
      enddo
      if(idh.eq.ndh) then
!       xexact(ipoints) = xs*dexp((a-0.5d0*dd*dd)*time + d1*dd*wtot)
!       xexact(ipoints) = xs*dexp(a*time + d1*dd*wtot)
      endif
      xsave(ipoints,idh) = x0
!      if(rank==root)      print *,time,x0,dh,w0
   enddo
  enddo
  if (rank.eq.root.and.iave.eq.1) then
     do ipoints=1,npoichk
        write(100,8) ipoints*dh*2**ndh,(xsave(ipoints,ii),ii=1,ndh)
     enddo
     close(100)
  endif
  if (rank.eq.root.and.iave.eq.nave) then
     do ipoints=1,npoichk
        write(101,8) ipoints*dh*2**ndh,(xsave(ipoints,ii),ii=1,ndh)
     enddo
     close(101)
  endif  
8 format(15(1x,e10.4))
! here we have the points, compute things
  do ipoints = 1,npoichk
   do idh=ndh0,ndh
    res_loc(ipoints,idh) = res_loc(ipoints,idh) + dabs(xsave(ipoints,ndh)-xsave(ipoints,idh))
   enddo
  enddo
 enddo
 call MPI_Reduce(res_loc, res_glo, npoichk*ndh, MPI_DOUBLE_PRECISION, MPI_SUM, root, MPI_COMM_WORLD, ierr)
 if (rank == root) then
   res_glo=res_glo/(size*nave)
   do idh = ndh0,ndh
    dh = dh0/2**idh
    do ipoints = 1,npoichk
     write(10+ipoints,*) dh,res_glo(ipoints,idh)
    enddo
   enddo
 endif
      
enddo


call MPI_Finalize(ierr)
end program stability
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

!mpif90 -O2 -o stability.x stability.f90

!To run:

!mpirun -np 4 ./stability.x

!Let me know if you need help adapting rnor, adding MPI I/O, or improving performance.

