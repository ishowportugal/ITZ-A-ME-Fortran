      PROGRAM AVGCHAN
C**********************************************************************
C     CALCULATES THE AVERAGE VALUE AT EACH TIME STEP OF THE VALUES
C     CONTAINED IN SEPARATE SINGLE COLUMN ASCII FILES
C     THE TOTAL NUMBER AND THE NAMES OF THE INPUT FILES MUST
C     BE LISTED IN THE FILE AVGCHAN
C     
C     ATHANASSIOS A. ARGIRIOU, 14 MARCH 1994
C**********************************************************************
      character*10 infil(10),filout
      dimension logun(10)
      open(100,file='avgchan.in')
      read(100,*)nfiles
      do i=1,nfiles
      read(100,'(A)')infil(i)
      write(*,'(a)')infil(i)
      logun(i)=101*i
      open(logun(i),file=infil(i))
      enddo
c
      write(*,*)' Enter output file name '
      read(*,'(A)')filout
      open(1111,file=filout)
c
  20  xtot=0.
      do i=1,nfiles
      read(logun(i),*,err=200)x
      xtot=xtot+x
      enddo
      xtot=xtot/float(nfiles)
      write(1111,*)xtot
      goto 20
 200  close(1111)
      do i=1,nfiles
      close(logun(i))
      enddo
      end
      