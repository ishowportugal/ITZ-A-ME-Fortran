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
      FUNCTION RE(V,XL,DEN,VIS)
C**************************************************************
C CALCULATION OF THE REYNOLD'S NUMBER
C ATHANASSIOS ARGIRIOU, 22 - 5 - 1990
C***************************************************************
      RE=V*XL*DEN/VIS
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
      SUBROUTINE HUMRAT(PW,T,PA,W)
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
      CALL HUMRAT(PW,T,PA,W)
      CALL HUMRAT(PWSAT,TWB,PA,WS)
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