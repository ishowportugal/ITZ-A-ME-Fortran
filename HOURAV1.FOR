      program hourav1
c**********************************************************************
c     Calculation of hourly averages from 10 min data of PASSYS test cell
c     experiments. The time is not written in the output file as with the 
c     code houravg.for
c     The 10 min values files including data of one sensor must be first 
c     joined with the time values, using the software join.for
c
c     Athanassios A. Argiriou, 13 March 1994
c************************************************************************
      character filin*10,filout*10
c
      write(*,*)' Enter the input file name '
      read(*,'(A)')filin
      open(100,file=filin)
      write(*,*)' Enter the output file name '
      read(*,'(A)')filout
      open(200,file=filout)
c
      icount=0
      x1=0.
      read(100,*,err=100)time,x
   10 call ticonv(time,iday,ihour,min)
      icount=icount+1
      min1=min
      iday1=iday
      ihour1=ihour
      icheck=ihour1+1
      if((ihour1+1).eq.24)icheck=0
      x1=x1+x
  20  read(100,*,err=100)time,x
      call ticonv(time,iday,ihour,min)
c
      if(ihour.eq.icheck.and.min.gt.30)then
      x1=x1/float(icount)
      write(200,*)x1
      icount=0
      x1=0
      goto 10
      else
      x1=x1+x
      icount=icount+1
      goto 20
      endif
100   continue
      close(100)
      close(200)
      end
      