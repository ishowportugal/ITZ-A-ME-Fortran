      PROGRAM EVACOOL
C*******************************************************************
C This program allows the evaluation of the performance of an evapo-
C rative cooling system
C Athanassios Argiriou, 3 - 6 - 1990
C*****************************************************************
C
      CHARACTER LOCAT*10, INAM*6
      DIMENSION TA(24),RH(24)
      DATA PA/1./,VFAN,QF/1500.,100./
C 
C Functions
C
      EFF(TA,TO,TWB)=(TA-TO)/(TA-TWB)
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
      WRITE(*,9001)LOCAT,MONTH
 9001 FORMAT(/,' Direct evaporative cooling for ',A10,'  Month ',I2,/)
C
      CALL QAIR(QF,VFAN,QA)
      WRITE(*,9002)VFAN,QF,QA
 9002 FORMAT(' Fan speed = ',F5.0,' r.p.m., Water flow rate = '
     *,F4.0,' kg/h, Air flow rate =',F5.1,' kg/h',/)
      WRITE(*,9005)
 9005 FORMAT(' HOUR   T IN    RH IN   TWB    T OUT    RH OUT      Q0')
      WRITE(*,9006)
 9006 FORMAT('          C       %      C       C         %         kJ')
C
C
      DO 10 I=1,24
C
      CALL WETBUL(TA(I),RH(I),PA,TWB)
      CALL TAIROUT(QF,TA(I),TWB,VFAN,TO)
      CALL RELOUT(TWB,PA,TO,RHO)
      CALL EVAWAT(QF,TA(I),TWB,VFAN,QE)
      Q0=QA*CPAIR((TA(I)+TO)/2.)*(TA(I)-TO)/1000.
C
      WRITE(*,9004)I,TA(I),RH(I),TWB,TO,RHO,Q0
 9004 FORMAT(2X,I2,4X,F4.1,5X,F4.1,4X,F4.1,4X,F4.1,5X,F4.1,6X,F9.1)
      Q=Q+Q0
  10  CONTINUE
      Q=Q/3600.
      WRITE(*,9100)Q
 9100 FORMAT(/,' DAILY OFFERED COOLING ENERGY = ',F9.1,' kWh',/)
 20   CONTINUE
      CALL FICHR(100,INAM,'OLD',2)
      END
C
C
      SUBROUTINE EVAWAT(QF,TA,TWB,VFAN,QE)
C****************************************************************
C Calculation of the evaporated water QE (kg/h) during the cooling
C process
C Athanassios Argiriou, 3 - 6 - 1990
C*****************************************************************
C
      DATA C1,C2,C3,C4,C5/.07,.14,.002E-04,1.4,.4/
C
      QFF=QF**C2
      DT=TA-TWB
      QE=C1*QFF*DT+C3*QFF*VFAN**C4*DT**C5
      RETURN
      END
C
C
      SUBROUTINE QAIR(QF,VFAN,QA)
C*******************************************************************
C Air mass flow rate QA (kg/h) at the output of the cooler
C Athanassios Argiriou, 3 - 6 - 1990
C*******************************************************************
C
      DATA C1,C2,C3/-39.7,1.46E-05,.2/
C
      QA=C1+C2*QF*VFAN+C3*VFAN
      RETURN
      END
C
C
      SUBROUTINE TAIROUT(QF,TA,TWB,VFAN,TO)
C*******************************************************************
C Calculation of the outlet air temparature TO (C), from the cooler
C Athanassios Argiriou, 3 - 6 - 1990
C*******************************************************************
C
      DATA C1,C2,C3,C4,C5,C6/1.1,.23,.09,1.18E-04,2.16,-.61/
C
      T1=QF**C3
      T2=TA-TWB
      TO=TA-C1-C2*T1*T2-C4*T1*T2**C5/VFAN**C6
      RETURN
      END
C
C
      SUBROUTINE RELOUT(TWB,PA,T,RHO)
C*******************************************************************
C Calculation of the relative humidity (%) of the outlet air
C Athanassios Argiriou, 10 - 6 - 1990
C********************************************************************
C
      DATA C1,C2,C3,C4/1093.,.444,.24,.556/
C
      TF(X)=X*1.8+32.
C
      PWSAT=PWS(TWB)
      CALL HUMRAT(PWSAT,TWB,PA,WS)
C
      TWBF=TF(TWB)
      TAF=TF(T)
      WO=((C1-C4*TWBF)*WS-C3*(TAF-TWBF))/(C1+C2*TAF-TWBF)
C relative humidity of the outlet air
      PWO=WO*PA/(.62198+WO)
      PWSATO=PWS(T)
      RHO=PWO/PWSATO*100.
C
      RETURN
      END
