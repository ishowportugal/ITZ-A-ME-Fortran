      program makecas
C***********************************************************************
C     THIS CODE CREATES A TIME STEP CASUAL GAINS FILE COMPATIBLE TO ESP-R
C
C     ATHANASSIOS A. ARGIRIOU, 12-4-94
C***********************************************************************
      CHARACTER*12 INFIL,OUTFIL
C
      write(*,'(A)')' Enter the input file name'
      READ(*,'(A)')INFIL
C
      write(*,'(A)')' Enter the output file name'
      READ(*,'(A)')OUTFIL

      OPEN(100,FILE=INFIL)
      OPEN(200,FILE=OUTFIL)
C
      NX1=0
      NX2=0
  20  READ(100,*,ERR=30)X
      NX=NINT(X)     
      WRITE(200,'(3I7)')NX,NX1,NX2
      GOTO20
 30   CONTINUE
      CLOSE(100)
      CLOSE(200)
      END
