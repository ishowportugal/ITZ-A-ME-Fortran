      program join
c *********************************************************************
c     JOINS VARIOUS ASCII SINGLE COLUMN FILES IN PARALLEL
C     (AT THE SAME TIME STEP)
C
C     ATHANASSIOS A. ARGIRIOU, 13 MARCH 1994
C**********************************************************************
      character ifile(10)*10,IND*1,FILOUT*10
      DIMENSION LOGIC(10),X(10)
      write(*,*)' Enter the number of files to be joined'
      read(*,*)nfiles
  10  DO I=1,NFILES
      WRITE(*,*)' ENTER NAME OF FILE ',i
      READ(*,'(A)')IFILE(I)
      ENDDO
      WRITE(*,*)' THE FILES TO BE JOINED ARE '
      DO I=1,NFILES
      WRITE(*,*)IFILE(I)
      ENDDO
      WRITE(*,*)' IF OK THEN Y(ES)'
      READ(*,'(A)')IND
      IF(IND.NE.'Y')GOTO 10
C
      WRITE(*,*)' ENTER THE NAME OF THE OUTPUT FILE'
      READ(*,'(A)')FILOUT
      OPEN(999,FILE=FILOUT)
C
      DO I=1,NFILES
      LOGIC(I)=100*I
      OPEN(LOGIC(I),FILE=IFILE(I))
      ENDDO
   20 CONTINUE
      DO I=1,NFILES
      READ(LOGIC(I),*,ERR=100)X(I)
      ENDDO
      WRITE(999,*)(X(I),I=1,NFILES)
      GOTO 20
 100  CONTINUE
      DO I=1,NFILES
      CLOSE(LOGIC(I))
      ENDDO
      CLOSE(999)   
      END