program stabilityarea
implicit real*8 (a-h,o-z)
integer, parameter :: nd=150, nbin=2*nd+1, nnd=20, ndd=20, nalgo=17
! ialgo 1-nito      Ito schemes
! ialgo nstra-nalgo Stratonovich schemes
integer, parameter :: nito = 5, nstra = nito+1
integer :: algo(nalgo)
character :: string*120

read(70,*) string
do i = 1,10000
   read(70,*,end=12) dh,dd,(algo(ialgo),ialgo=1,nalgo)
   do ialgo=1,nalgo
      if(algo(ialgo).eq.0) then
         write(700+ialgo,*) dh,dd
        else
         write(800+ialgo,*) dh,dd
      endif
   enddo
enddo
12 stop
end

   
