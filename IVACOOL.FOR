      PROGRAM IVACOOL
C*******************************************************************
C This program allows the evaluation of the performance of an
C indirect evaporative cooler
C Athanassios Argiriou, 14 - 6 - 1990
C*****************************************************************
C
      CHARACTER LOCAT*10, INAM*6
      DIMENSION TA(24),RH(24)
       DATA PA/1./,VA/1./
C 
C Functions
C
      EFF(VA)=1./(1.+.47*VA**.3)
C
C INPUT PARAMETERS
C
      WRITE(*,9200)
 9200 FORMAT(' ENTER THE CLIMATE FILE NAME')
      READ(*,9201)INAM
 9201 FORMAT(A6)
      CALL FICHR(100,INAM,'OLD',1)
      READ(100,*)NCIT
      DO 20 L=1,NCIT
      READ(100,9000)LOCAT
 9000 FORMAT(A10)
      DO 20 LM=1,5
      READ(100,*)MONTH
      READ(100,*)(TA(J),J=1,24)
      READ(100,*)(RH(J),J=1,24)
C
C INITIALISATION
C
      Q=0.
      EFI=EFF(VA)
      WRITE(*,9001)LOCAT,MONTH
 9001 FORMAT(/,' Indirect evaporative cooling for ',A10,' Month ',
     &I2,/)
C
      WRITE(*,9002)VA,EFI
 9002 FORMAT(' Air velocity = ',F2.0,' m/s',' Efficiency = ',F3.2,/)
      WRITE(*,9003)
 9003 FORMAT('               INLET              OUTLET')
      WRITE(*,9005)
 9005 FORMAT(' HOUR     T     RH    TWB     T    RH    TWB        Q')
      WRITE(*,9006)
 9006 FORMAT('          C      %     C      C     %     C      Wh/m2'
     &)
C
C
      DO 10 I=1,24
C
      CALL WETBUL(TA(I),RH(I),PA,TWB)
      TO=TA(I)-EFI*(TA(I)-TWB)
      RH2=RH(I)*PWS(TA(I))/PWS(TO)
      CALL WETBUL(TO,RH2,PA,TWBO)
      TMOY=(TA(I)+TO)/2.
      Q0=VA*DAIR(TMOY)*CPAIR(TMOY)*(TA(I)-TO)
C
      WRITE(*,9004)I,TA(I),RH(I),TWB,TO,RH2,TWBO,Q0
 9004 FORMAT(2X,I2,6(3X,F4.1),F9.1)
C
      Q=Q+Q0
  10  CONTINUE
      WRITE(*,9100)Q
 9100 FORMAT(/,' DAILY OFFERED COOLING ENERGY = ',F9.1,' Wh/m2',/)
 20   CONTINUE
      CALL FICHR(100,INAM,'OLD',2)
      END
