       PROGRAM SKYTEM
C
C***************************************************************
C THE PROGRAM CALCULATES THE SKY EMMISSIVITY AND THE SKY TEMPERA
C TURE DEPRESSION OF A GIVEN LOCATION
C ATHANASSIOS ARGIRIOU UNIVESRITY OF ATHENS 7 - 7 - 1994
C CEC DG XII - PASCOOL PROJECT
C***************************************************************
C
C DEFINITION OF WEATHER CONDITIONS FOR THE CALCULATION
C
      WRITE(*,*)' ENTER THE HOUR OF THE DAY'
      READ(*,*)TIME
      WRITE(*,*)'ENTER THE AMBIENT TEMPERATURE (C)'
      READ(*,*)TA
      WRITE(*,*)' RELATIVE HUMIDITY (%)'
      READ(*,*)RH
      WRITE(*,*)' ENTER THE CLOUD COVER (0 - 10)'      
      READ(*,*)OSC
C
      OSC=OSC/10.
C         
C CALCULATION OF THE THERMOPHYSICAL PARAMETERS
C
      TABS=TA+273.15
      ESKY=EBERD(TDEW(RH,TA),TIME)*CA(OSC)
C
C SKY TEMPERATURE DEPRESSION
C
      DTSKY=(1.-ESKY**.25)*TABS
C
      WRITE(*,1000)ESKY,DTSKY
 1000 FORMAT(' THE SKY EMISSIVITY IS:', F5.3,/,' THE SKY TEMPERATURE DEP
     &RESSION IS:', F5.1, ' C')
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
      TDEW=237.3*(D2+D1)/(CTE-D2+D1)
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
      