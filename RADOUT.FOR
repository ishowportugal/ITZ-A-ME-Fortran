      PROGRAM RADOYT
C
C***************************************************************
C THE PROGRAM CALCULATES THE COOLING POTENTIAL FOR A GIVEN AREA
C THE COOLING POTENTIAL IS EVALUATED IN TERMS OF THE STAGNATION 
C TEMPERATURE AND THE NET RADIATIVE FLUX OF A RADIATOR
C ATHANASSIOS ARGIRIOU 7 - 7 - 1994
C CEC DG XII PASCOOL PROJECT
C***************************************************************
C
C CONSTANTS
C
      DATA SBOLTZ/5.67E-08/
C
C RADIATOR CHARACTERISTICS
C
      DATA ERAD/.9/,XR,ZR,DZR/2.,.01,.002/,CSTEEL/14.6/
C
C WINDSCREEN CHARACTERISTICS
C
       DATA TRANS/.75/
C
C FUNCTIONS
C
      TFOUT(T,TS,D,V,CP,U,X,Z)=TS+(T-TS)*EXP(-U*X/D/V/CP/Z)
C
C
      WRITE(*,*) ' ENTER THE AIR VELOCITY UNDER THE PLATE (m/s)'
      READ(*,*)FV
C
C DEFINITION OF WEATHER CONDITIONS FOR THE CALCULATION
C
      WRITE(*,*)' ENTER THE HOUR OF THE DAY'
      READ(*,*)TIME
      WRITE(*,*)'ENTER THE AMBIENT TEMPERATURE (C)'
      READ(*,*)TA
      WRITE(*,*)' RELATIVE HUMIDITY (%)'
      READ(*,*)RH
      WRITE(*,*)' WIND VELOCITY (m/s)'
      READ(*,*)WV
      WRITE(*,*)' ENTER THE CLOUD COVER (0 - 10)'      
      READ(*,*)OSC
C
      OSC=OSC/10.
C         
C CALCULATION OF THE THERMOPHYSICAL PARAMETERS
C
      TABS=TA+273.15
      ESKY=EBERD(TDEW(RH,TA),TIME)*CA(OSC)
      HE0=HW0(WV)+4.*ERAD*SBOLTZ*TABS**3
      HE1=HW1(WV)+4.*TRANS*ERAD*SBOLTZ*TABS**3
      CALL HFLU(TA,XR,FV,HF)
      U0=1./(1./HE0+1./HF+DZR/CSTEEL)
      U1=1./(1./HE1+1./HF+DZR/CSTEEL)
C
C SKY TEMPERATURE DEPRESSION
C
      DTSKY=(1.-ESKY**.25)*TABS
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
C RESULTS
C
      WRITE(*,1001)TSTAG0
 1001 FORMAT(' STAGNATION TEMPERATURE WITHOUT WIND SCREEN', F5.1,' C')

      WRITE(*,1002)TSTAG1
 1002 FORMAT(' STAGNATION TEMPERATURE WITH WIND SCREEN', F5.1,' C')

      WRITE(*,1003)TF0
 1003 FORMAT(' OUTLET TEMPERATURE WITHOUT WIND SCREEN', F5.1,' C')

      WRITE(*,1004)TF1
 1004 FORMAT(' OUTLET TEMPERATURE WITH WIND SCREEN', F5.1,' C')
C
      END
C
      SUBROUTINE HFLU(TA,X,FV,HF)
C***************************************************************
C CONVECTIVE HEAT TRANSFER COEFFICIENT FOR FORCED CONVECTION
C UNDER A HORIZONTAL PLATE
C ATHANASSIOS ARGIRIOU, 30 - 5 - 1990
C***************************************************************
C
      REY=RE(FV,X,DAIR(TA),DVAIR(TA))
      PRA=PR(CPAIR(TA),DVAIR(TA),CAIR(TA))
      IF (REY.LE.500000.)THEN
      HF=.664*REY**.5*PRA**(1./3.)*CAIR(TA)/X
      ELSE
      HF=PRA**(1./3.)*(.036*REY**.8-836.)*CAIR(TA)/X
      ENDIF
      RETURN
      END
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
C
      FUNCTION HRG(TRAD,Z1)
      HRG=CAIR(TRAD)/Z1
      RETURN 
      END
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
      