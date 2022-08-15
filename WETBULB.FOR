	program wetbulb
c---------------------------------------------------------------------
c Adds to the 20 year meteo file the calculated wet bulb temperature
c
c A. Argiriou, 26-3-99
c----------------------------------------------------------------------

	integer x(10)

	open (100, file='1977to1996')
	open (200, file='wetbul.dat')

 10	read (100,*,end=100) (x(i),i=1,10)
	T=real(x(8))/10.
	RH=real(x(9))
	PA=real(x(10))/10./760.
	call WETBUL(T,RH,PA,TWB)
	write (200,'(i4,3(2x,i2),2(2x,f4.1))')(x(i),i=1,4),T,TWB
	goto10
	
  100	close(100)
	close(200)
	stop 'Normal Termination'
	end
