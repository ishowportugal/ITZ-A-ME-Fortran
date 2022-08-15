       PROGRAM SHADIN
C******************************************************
C Design of shading devices using the one point method
C Athanassios Argiriou, 16.10.1990
C******************************************************
C
      WRITE(*,9000)
 9000 FORMAT(1X,' ENTER THE WIDTH OF THE WINDOW',/,
     &'          HEIGHT OF THE WINDOW',/,
     &'          AZIMUTH OF THE WINDOW (NORTH=0, CLOCKWISE)',/,
     &'          LATITUDE OF THE PLACE',/,
     &'          DAY,MONTH,HOUR OF MAX. SHADING')
      READ(*,*)WD,HE,WAZ,SLAT,D,M,HR
C
      CALL MPOINT(WD,HE,WAZ,SLAT,D,M,HR)
C
      END
C
      SUBROUTINE MPOINT(WD,HE,WAZ,SLAT,D,M,HR)
C
      PI=4.*ATAN(1.)
      R=PI/180.
      WAZ=WAZ*R
      CALL SANG(HR,D,M,SLAT,DEC,SALT,SAZM)
C
      IF (WAZ.GT.SAZM) THEN
         XM=HE*SIN(WAZ-SAZM)/TAN(SALT)-WD/2.
         YM=HE*COS(WAZ-SAZM)/TAN(SALT)
         ZM=0.
      ELSE
         XM=HE*SIN(WAZ-SAZM)/TAN(SALT)+WD/2.
         YM=HE*COS(SAZM-WAZ)/TAN(SALT)
         ZM=0.
      ENDIF
C
      IF (ABS(XM).GT.WD/2.) THEN
         IF (XM.GT.0.) THEN
            XM2=WD/2.
         ELSE
            XM2=-WD/2.
         ENDIF
         XM1=XM2
         YM2=WD*YM/(ABS(XM)+WD/2.)
         YM1=YM2
         ZM1=HE-YM2/COS(WAZ-SAZM)*TAN(SALT)
      WRITE(*,9100)XM2,YM2,XM1,YM1,ZM1
 9100 FORMAT(' X,Y FOR OVERHANG : ',2F7.2,/,' X,Y,Z  FOR SIDE FIN : ',
     &3F7.2)
      RETURN
      ELSE
      WRITE(*,9200)XM,YM
 9200 FORMAT(' X,Y FOR OVERHANG : ',2F7.2)
      ENDIF
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