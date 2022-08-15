      PROGRAM RADCOOL
C
C***************************************************************
C THE PROGRAM CALCULATES THE COOLING POTENTIAL FOR A GIVEN AREA
C THE COOLING POTENTIAL IS EVALUATED IN TERMS OF THE STAGNATION 
C TEMPERATURE AND THE NET RADIATIVE FLUX OF A RADIATOR
C ATHANASSIOS ARGIRIOU 7 - 5 - 1990
C***************************************************************
C
      CHARACTER LOCAT*10
      DIMENSION JFREQ(30)
C
C CONSTANTS
C
      DATA SBOLTZ/5.67E-08/
C
C RADIATOR CHARACTERISTICS
C
      DATA ERAD/.9/,XR,YR,ZR,DZR/10.,1.,.01,.002/,FV/2.5/
C
C WINDSCREEN CHARACTERISTICS
C
       DATA TRANS/.75/,DWSC/1.56E-02/,CWSC/.46/,Z1/.03/
C
C FUNCTIONS
C
      TFOUT(T,TS,D,V,CP,U,X,Z)=TS+(T-TS)*EXP(-U*X/D/V/CP/Z)
C
C DEFINITION OF WEATHER CONDITIONS FOR THE CALCULATION
C
 30   WRITE(*,9300)
 9300 FORMAT(' CALCULATIONS WITH AVERAGE WEATHER CONDITIONS (1)',/,
     &'           OR WITH OPTIMAL WEATHER CONDITIONS (2)?')
      READ(*,*)ICOND
      IF (ICOND.GT.2.OR.ICOND.LT.1) GO TO 30
C
C INPUT PARAMETERS
C
      OPEN(100,FILE='KLIMAT')
      READ(100,*)NCIT
      DO 20 L=1,NCIT
      READ(100,9000)LOCAT
 9000 FORMAT(A10)
      READ(100,*)SLAT
      DO 20 LL=1,5
      READ(100,*)DAY,MONTH,TMEAN,TAMPL,RH,WV,OSC
C
      IF (ICOND.EQ.2) OSC=0.0
C
      WRITE(*,9600)LOCAT,MONTH,FV
 9600 FORMAT(' Radiative coolers performance in ',A10,' MONTH '
     &,I2,' Fluid velocity ',F4.1,' m/s')
C
C DAILY COOLING PERIOD
C
      CALL COOLPER(DAY,MONTH,SLAT,ISTART,ISTOP)
C
C INITIALISATION
C
      WRITE(*,9100)
 9100 FORMAT(19X,'NO WINDSCREEN',20X,'WINDSCREEN')
      WRITE(*,9200)
 9200 FORMAT(1X,'HOUR  T AMB   DTSKY    T STAG   DT   T OUT',4X,
     *'T STAG   DT   T OUT')
      Q0=0.
      Q1=0.
      DO 40 NN=1,30
 40   JFREQ(NN)=0
C
      DO 10 I=ISTART,ISTOP
C
      IF (I.LE.24) THEN
      TIME=FLOAT(I)
      ELSE
      TIME=FLOAT(I-24)
      ENDIF
C
C AMBIENT TEMPERATURE CALCULATION
C
      TA=TSIN(TMEAN,TAMPL,TIME)
      TABS=TA+273.15
C
C CALCULATION OF THE THERMOPHYSICAL PARAMETERS
C
      ESKY=EBERD(TDEW(RH,TA),TIME)*CA(OSC)
      HE0=HW0(WV)+4.*ERAD*SBOLTZ*TABS**3
      HE1=HW1(WV)+4.*TRANS*ERAD*SBOLTZ*TABS**3
      CALL HFLU(TA,YR,ZR,FV,HF)
      U0=1./(1./HE0+1./HF+DZR/CSTEEL(TA))
      U1=1./(1./HE1+1./HF+DZR/CSTEEL(TA))
C
C SKY TEMPERATURE DEPRESSION
C
      DTSKY=(1.-ESKY**.25)*TABS
C HISTOGRAM
      IFREQ=INT(DTSKY)+1
      JFREQ(IFREQ)=JFREQ(IFREQ)+1
C
C STAGNATION TEMPERATURE
C
      TSTAG0=STAG(TA,ERAD,ESKY,HE0)
      TSTAG1=STAG(TA,ERAD*TRANS,ESKY,HE1)
C
C FLUID OUTLET TEMPERATURE
C
      TF0=TFOUT(TA,TSTAG0,DAIR(TA),FV,CPAIR(TA),U0,XR,ZR)
      TF1=TFOUT(TA,TSTAG1,DAIR(TA),FV,CPAIR(TA),U1,XR,ZR)
C
C COOLING POWER
C
      TM0=(TA+TF0)/2.
      TM1=(TA+TF1)/2.
      Q0=Q0+POWER(TA,TF0,FV,DAIR(TM0),CPAIR(TM0),XR,ZR)
      Q1=Q1+POWER(TA,TF1,FV,DAIR(TM1),CPAIR(TM1),XR,ZR)
C
C CONDENSATION CONTROL
C
      CALL DEWTEST(TA,TF1,FV,TM1,RH,HF,WV,Z1)
C
C RESULTS
C
      WRITE(*,9005)TIME,TA,DTSKY,TSTAG0,TA-TSTAG0,TF0,TSTAG1,
     *TA-TSTAG1,TF1
 9005 FORMAT(2X,F3.0,2(3X,F4.1),2(5X,F4.1,3(3X,F4.1)))
 10   CONTINUE
      WRITE(*,9500)Q0,Q1
 9500 FORMAT(/,1X,'DAILY COOLING LOAD WITHOUT WINDSCREEN  ',F7.2,
     *2X,'Wh/m2',/,1X,'DAILY COOLING LOAD WITH    WINDSCREEN  ',
     *F7.2,2X,'Wh/m2',/)
C SKY TEMPERATURE DEPRESSION FREQUENCIES
      WRITE(*,9990)
 9990 FORMAT(' INTERVAL   No EVENTS')
      DO 50 IFREQ=1,30
      IF (JFREQ(IFREQ).GT.0)THEN
      WRITE(*,9999)IFREQ,JFREQ(IFREQ)
 9999 FORMAT(5X,I2,17X,I2)
      ENDIF
 50   CONTINUE
      WRITE(*,9998)
 9998 FORMAT(/)
 20   CONTINUE
      CLOSE(100)
      END
C
      SUBROUTINE HFLU(TA,Y,Z,FV,HF)
C
      DH=2*Y*Z/(Y+Z)
      REY=RE(FV,Y,Z,DAIR(TA),DVAIR(TA))
      IF (REY .LE. 10000.)THEN
      HF=5.385*CAIR(TA)/DH
      ELSE
      HF=0.019068*(REY**.785884)*CAIR(TA)/DH
      ENDIF
      RETURN
      END
C
      SUBROUTINE COOLPER(DAY,MONTH,SLAT,ISTART,ISTOP)
C***************************************************************
C DETERMINATION OF THE NOCTURNAL COOLING PERIOD
C ATHANASSIOS ARGIRIOU, 30 - 5 - 1990
C***************************************************************
C
      CALL SDEC(DAY,MONTH,DEC)
      CALL RISET(SLAT,DEC,SR,SST,DLG)
      DAY1=DAY+1.
      CALL SDEC(DAY1,MONTH,DEC1)
      CALL RISET(SLAT,DEC1,SR1,SST1,DLG1)
      ISTART=IFIX(SST)+1
      ISTOP=IFIX(SR1)+24
      RETURN
      END
C
C
      FUNCTION STAG(T,ER,ES,HE)
C
C CALCULATION OF THE STAGNATION TEMPERATURE FOR AN UNCOVERED
C RADIATOR
C
      DATA SBOLTZ/5.67E-08/
      STAG=T-ER*SBOLTZ*(1.-ES)*(T+273.15)**4/HE
      RETURN
      END
C
C
      FUNCTION POWER(TA,TF,FV,D,CP,X,Z)
C
C CALCULATION OF THE COOLING LOAD (Wh/m2)
C
      POWER=D*FV*Z*CP*(TA-TF)/X
      RETURN
      END
C
C
      FUNCTION HRG(TRAD,Z1)
      HRG=CAIR(TRAD)/Z1
      RETURN 
      END
C
C
      SUBROUTINE DEWTEST(TA,TF1,FV,TM1,RH,HF,WV,Z1)
C******************************************************************
C ROUTINE FOR CHECKING THE POSSIBILITY OF DEW FORMATION ON THE
C WINDSCREEN OF A RADIATOR
C ATHANASSIOS ARGIRIOU, 30 - 5 - 1990
C******************************************************************
C
      TRAD=TM1-DAIR(TM1)*FV*(TA-TF1)/HF
      TWSC=(TA*HW0(WV)+TRAD*HRG(TRAD,Z1))/(HW0(WV)+HRG(TRAD,Z1))
      IF((TWSC-TDEW(RH,TA)).GT.2.)THEN
      GO TO 800
      ELSE
         IF((TWSC-TDEW(RH,TA)).GT.1.)THEN
         WRITE(*,9800)
 9800    FORMAT(1X,'DANGER FOR DEW FORMATION')
         ELSE
         WRITE(*,9810)
 9810    FORMAT(1X,'ATTENTION!!-DEW FORMATION')
         ENDIF
      ENDIF
 800  CONTINUE     
      RETURN
      END
C*****************************************************
      SUBROUTINE SDEC(DAY,MONTH,DEC)
C*****************************************************
C THIS ROUTINE CALCULATES SOLAR DECLINATION (DEC)
C ATHANASSIOS ARGIRIOU 2 - 5 - 1990
C
      PI=4.*ATAN(1.)
C
C CALCULATION OF DAY NO
C
      GOTO(1,2,3,4,5,6,7,8,9,10,11,12)MONTH
 1    DA=DAY
      GOTO 13
 2    DA=DAY+31.
      GOTO 13
 3    DA=DAY+59.
      GOTO 13
 4    DA=DAY+90.
      GOTO 13
 5    DA=DAY+120.
      GOTO 13
 6    DA=DAY+151.
      GOTO 13
 7    DA=DAY+181.
      GOTO 13
 8    DA=DAY+212.
      GOTO 13
 9    DA=DAY+243.
      GOTO 13
 10   DA=DAY+273.
      GOTO 13
 11   DA=DAY+304.
      GOTO 13
 12   DA=DAY+334.
C CALCULATION OF DECLINATION
 13   DEC=23.44*SIN(2.*PI*(DA-81.)/365.)
      RETURN
      END
C
C***********************************************************
      SUBROUTINE SANG(SHOUR,DAY,MONTH,SLAT,DEC,SALT,SAZM)
C***********************************************************
C
C THIS ROUTINE CALCULATES THE SOLAR ALTITUDE (SALT) & SOLAR AZIMUTH
C (SAZM) USING THE ROUTINE SDEC TO CALCULATE THE DECLINATION. THE SOL
C HOUR (SHOUR) AND THE LATITUDE OF THE PLACE (SLAT) IS USED AS INPUT PAR
C AMETER
C
C ATHANASSIOS ARGIRIOU 2 - 5 - 1990
C
      PI=4.*ATAN(1.)
      R=PI/180.
C
C DECLINATION CALCULATION
C
      CALL SDEC(DAY,MONTH,DEC)
C
C ALTITUDE CALCULATION
C
      DEC=DEC*R
      SLAT=SLAT*R
      HR=HRANG(SHOUR)*R
      SALT=ASIN(SIN(DEC)*SIN(SLAT)+COS(DEC)*COS(SLAT)*COS(HR))
C
C AZIMUTH CALCULATION
      SAZM=ASIN(-SIN(HR)*COS(DEC)/COS(SALT))
      IF (SHOUR.LE.12.)SAZM=PI-ABS(SAZM)
      IF (SHOUR.GT.12.)SAZM=PI+ABS(SAZM)
C
      RETURN
      END
C
C*********************************************************
      FUNCTION HRANG(X)
C*********************************************************
C THIS FUNCTION CALCULATES THE HOUR ANGLE CORRESPONDING TO THE SOLAR
C HOUR X.
C ATHANASSIOS ARGIRIOU, 2 - 5 - 1990
C
      IF(X.GE.12.)THEN
      HRANG=15.*(X-12.)
      ELSE
      HRANG=15.*(X+12.)
      ENDIF
      RETURN
      END
C
C
      SUBROUTINE RISET(SLAT,DEC,SR,SST,DLG)
C************************************************************
C CALCULATION OF THE SUNRISE AND SUNSET HOURS AND OF THE 
C DAY LENGTH
C SR:SUNRISE HOUR, SST:SUNSET HOUR, DLG:DAYLENGTH (hrs)
C ATHANASSIOS ARGIRIOU, 21 - 5 - 1990
C*************************************************************
C
      PI=4.*ATAN(1.)
      R=PI/180.
C
      SR=(360.-ACOS(-TAN(SLAT*R)*TAN(DEC*R))/R)/15.-12.
      DLG=(12.-SR)*2.
      SST=SR+DLG
C
      RETURN
      END
C
C

      FUNCTION TDEW(RH,TA)
C
C************************************************************
C CALCULATION OF THE DEW POINT TEMPERATURE AS A FUNCTION OF
C THE RELATIVE HUMIDITY (RH) AND THE AMBIENT TEMPERATURE (TA)
C RH %, TA deg C
C ATHANASSIOS ARGIRIOU, 8 - 5 - 1990
C************************************************************
C
      CTE=17.2693882
      RHD=RH/100.
      D1=CTE*TA/(TA+237.3)
      D2=ALOG(RHD)
      TDEW=237.3*(D2+D1)/(CTE-D2-D1)
      RETURN
      END
C
C
      FUNCTION EBERD(TD,TI)
C
C***************************************************************
C BERDAHL & MARTIN FORMULA
C CALCULATION OF THE SKY EMISSIVITY AS A FUNCTION OF THE DEW
C POINT TEMPERATURE (TD, deg C).
C ATHANASSIOS ARGIRIOU, 8 - 5 - 1990.
C***************************************************************
C
      EBERD=0.711+.56*TD/100.+.73*(TD/100.)**2+
     *.013*COS(2.*4.*ATAN(1.)*TI/24.)
      RETURN
      END
C
C
      FUNCTION HW0(X)
C
C***************************************************************
C THIS FUNCTION GIVES THE CONVECTIVE HEAT TRANSFER COEFFICIENT
C FOR A PLATE EXPOSED TO WIND.
C ATHANASSIOS ARGIRIOU, 10 - 5 - 1990.
C*****************************************************************
C
      IF(X.LE.4.)THEN
      HW0=5.7+3.8*X
      ELSE
      HW0=7.3*X**.8
      ENDIF
      RETURN
      END
C
C
      FUNCTION HW1(X)
C
C***************************************************************
C THIS FUNCTION GIVES THE CONVECTIVE HEAT TRANSFER COEFFICIENT
C FOR A SINGLE COVERED RADIATOR
C ATHANASSIOS ARGIRIOU, 30 - 5 - 1990.
C*****************************************************************
C
      HW1=.5+1.2*SQRT(X*3.6)
      RETURN
      END
C
C
      FUNCTION CA(X)
C*****************************************************************
C CORRECTION FOR THE CALCULATION OF THE SKY EMISSIVITY FOR A CLOU-
C DY SKY AS A FUNCTION OF THE TENTHS OF OPAQUE SKY COVER X
C 0<=X<=10
C ATHANASSIOS ARGIRIOU, 14 - 5 - 1990
C*****************************************************************
C
      IF(X.GE.0.AND.X.LE.10.)THEN
      CA=1.+.0224*X-.0035*X*X+.00028*X**3
      ELSE
      WRITE(*,9000)
 9000 FORMAT(1X,'TENTHS OF SKY COVER OUT OF RANGE')
      STOP
      ENDIF
      RETURN
      END
C
      FUNCTION DAIR(T)
C
C***********************************************************
C AIR DENSITY (kg/m3) AT 1,013 bar, -27 C<T<77 C
C ATHANASSIOS ARGIRIOU, 19 - 5 - 1990
C************************************************************
C
      DAIR=1.2997*EXP(-.35779E-02*T)
      RETURN
      END
C
C
      FUNCTION CPAIR(T)
C
C***********************************************************
C AIR SPECIFIC HEAT (J/kg/C) AT 1,013 bar, 0 C<T<77 C
C ATHANASSIOS ARGIRIOU, 19 - 5 - 1990
C************************************************************
C
      IF (T.GE.0.AND.T.LE.27.) THEN
      CPAIR=1003.2+.066796*T
      ELSE
      CPAIR=1003.6+.051484*T
      ENDIF
      RETURN
      END
C
C
      FUNCTION DVAIR(T)
C
C***********************************************************
C AIR DYNAMIC VISCOSITY (kg/m/s) AT 1,013 bar, 0 C<T<77 C
C ATHANASSIOS ARGIRIOU, 19 - 5 - 1990
C************************************************************
C
      DVAIR=.17175E-04+.4526E-07*T
      RETURN
      END
C
C
      FUNCTION CAIR(T)
C
C***********************************************************
C AIR THERMAL CONDUCTIVITY (W/m/C) AT 1,013 bar, 0 C<T<40 C
C ATHANASSIOS ARGIRIOU, 19 - 5 - 1990
C************************************************************
C
      CAIR=.02328+.8475E-04*T
      RETURN
      END
C
C
      FUNCTION TSIN(TMEAN,TAMPL,TIME)
C***********************************************************
C SINUSOIDAL TEMPERATURE VARIATION AS A FUNCTION OF TIME
C ATHANASSIOS ARGIRIOU, 22 - 5 - 1990
C************************************************************
C
      PI=4.*ATAN(1.)
      TSIN=TMEAN+TAMPL*COS(2.*PI*TIME/24.+7.*PI/12.)
      RETURN
      END
C
C
      FUNCTION CSTEEL(T)
C*****************************************************************
C THERMAL CONDUCTIVITY OF STEINLESS STEEL AS A FUNCTION OF TEMPERA
C TURE T (C)
C ATHANASSIOS ARGIRIOU, 22 - 5 - 1990
C******************************************************************
C
      CSTEEL=14.59+.0139*T
      RETURN
      END
C
C
      FUNCTION RE(V,Y,Z,DEN,VIS)
C**************************************************************
C CALCULATION OF THE REYNOLD'S NUMBER
C ATHANASSIOS ARGIRIOU, 22 - 5 - 1990
C***************************************************************
      DH=2*Z*Y/(Y+Z)
      RE=V*DH*DEN/VIS
      RETURN
      END
C
C
      FUNCTION PR(CP,VIS,COND)
C**************************************************************
C CALCULATION OF THE PRANDTL NUMBER
C ATHANASSIOS ARGIRIOU, 22 - 5 - 1990
C***************************************************************
      PR=CP*VIS/COND
      RETURN
      END
C
C
      FUNCTION PWS(T)
C*******************************************************************
C SATURATION PRESSURE OVER LIQUID WATER FOR THE TEMPERATURE RANGE
C OF -50 C TO 100 C
C SATURATION VAPOUR PRESSURE (atm)
C ATHANASSIOS ARGIRIOU, 3 - 6 - 1990
C*****************************************************************
C
      DATA C1,C2,C3,C4,C5,C6,C7/10.79586,5.02808,1.50474E-04,
     *-8.29692,.42873E-03,4.76955,2.219583/
      T1=273.16/(273.16+T)
      PWS=10.**(C1*(1.-T1)+C2*ALOG10(T1)+
     *C3*(1.-10.**(C4*((1./T1)-1.)))+C5*10.**(C6*(1.-T1)-1.)-C7)
      RETURN
      END
C
C
      SUBROUTINE HUMRAT(PW,PA,W)
C***************************************************************
C CALCULATION OF THE HUMIDITY RATIO (W) AS A FUNCTION OF THE PAR
C TIAL PRESSURE OF WATER VAPOUR (PW, atm), AMBIENT TEMPERATURE
C (T, C) AND ATMOSPHERIC PRESSURE (PA, atm)
C ATHANASSIOS ARGIRIOU, 3 - 6 - 1990
C***************************************************************
C
      W=.62198*PW/(PA-PW)
      RETURN
      END
C
C
      SUBROUTINE WETBUL(T,RH,PA,TWB)
C*******************************************************************
C CALCULATION OF THE WET BULB TEMPERATURE AS A FUNCTION OF THE DRY
C BULB TEMPERATURE (T, C), REL. HUMIDITY (%), PRESSURE (PA, atm)
C ATHANASSIOS ARGIRIOU, 3 - 6 - 1990
C*******************************************************************
C
C
      DATA C1,C2,C3,C4/1093.,.444,.24,.556/
C
      TF(X)=X*1.8+32.
      TC(X)=(X-32.)/1.8
C
      TWB=T
 10   PW=PWS(T)*RH/100.
      PWSAT=PWS(TWB)
      CALL HUMRAT(PW,PA,W)
      CALL HUMRAT(PWSAT,PA,WS)
C
      TWB1=TC((W*(C1+C2*TF(T))-C1*WS+C3*TF(T))/(W-C4*WS+C3))
C
      TE=ABS(TWB-TWB1)
      IF(TE.LE..1)THEN
      GO TO 20
      ELSE
      IF(TWB.GT.TWB1)THEN
      TWB=TWB-1.
      GO TO 10
      ELSE
      IF(TWB.LT.TWB1)THEN
      TWB=TWB+.01
      GO TO 10
      ENDIF
      ENDIF
      ENDIF
 20   CONTINUE
      RETURN
      END
C
C
