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
      ISTOP=IFIX(SR1)
      RETURN
      END

