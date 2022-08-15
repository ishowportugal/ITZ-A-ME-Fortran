      PROGRAM AIRHUM
C DETERMINATION OF THE MOIST AIR PROPERTIES USING ASHRAE CORRELATIONS.

      REAL PSYDAT(9),CEN(10)
      OPEN(11,FILE='METEO.MET')
      OPEN(12,FILE='RESULT.AIR')
C      DO J=1,9
C         PSYDAT(J)=0
C      END DO

 158  WRITE(*,15)'(METEO.MET) is the file which contains climatic data'
      WRITE(*,15)'(RESULT.AIR) is the file which contains results'
      WRITE(*,*)'ENTER THE SIMULATION TIME (Hours)'
      READ(*,*)KTP
      WRITE(*,*)'ENTER THE TOTAL SYSTEM PRESSURE (Atm)'
      READ(*,*) TSP
      WRITE(*,*)'ENTER WHICH UNITS DO YOU USE (1 For SI 2 For English)'
      READ(*,*)IUNT
      WRITE(*,*)'ENTER WHICH MODE DO YOU USE (1 to 6 )'
      WRITE (*,*) '1 for dry bulb and wet bulb'
      WRITE (*,*) '2 for dry bulb and relative humidity'
      WRITE (*,*) '3 for dry bulb and dew point'
      WRITE (*,*) '4 for dry bulb and humidity rate'
      WRITE (*,*) '5 for dry bulb and enthalpy'
      WRITE (*,*) '6 for humidity rate and enthalpy'
      READ(*,*)  MT
      WRITE(*,*) 'Enter the number of rows in climatic file '
      READ(*,*) NN
      WRITE(*,*) 'Enter the row number of the first climatic data '
      READ(*,*) NC1
      WRITE(*,*) 'Enter Multiplying coefficient'
      READ(*,*) AA
      WRITE(*,*) 'Enter the row number of the second climatic data '
      READ(*,*) NC2
      WRITE(*,*) 'Enter Multiplying coefficient'
      READ(*,*) BB
      PSYDAT(1)=TSP
      PSYDAT(8)=1.18
      PSYDAT(9)=1.21
      DT=1.
      WRITE(12,*)'  P(atm)  Tdry(øC)  Twb(øC)    RH     Tdp(øC)   W(e/s)
     & H(j/gas)  RauD     RauM'
      DO K=1,KTP
      TIME=K*DT
      READ(11,*)(CEN(I),I=1,NN)
        IF (MT.EQ.1) THEN
           PSYDAT(2)=CEN(NC1)*AA
           PSYDAT(3)=CEN(NC2)*BB

        ELSEIF (MT.EQ.2) THEN
           PSYDAT(2)=CEN(NC1)*AA
           PSYDAT(4)=CEN(NC2)*BB

        ELSEIF (MT.EQ.3) THEN
           PSYDAT(2)=CEN(NC1)*AA
           PSYDAT(5)=CEN(NC2)*BB

        ELSEIF (MT.EQ.4) THEN
           PSYDAT(2)=CEN(NC1)*AA
           PSYDAT(6)=CEN(NC2)*BB

        ELSEIF (MT.EQ.5) THEN
           PSYDAT(2)=CEN(NC1)*AA
           PSYDAT(7)=CEN(NC2)*BB

        ELSEIF (MT.EQ.6) THEN
           PSYDAT(6)=CEN(NC1)*AA
           PSYDAT(7)=CEN(NC2)*BB
        ELSE
        GO TO 158
        ENDIF
      CALL PSYCH(TIME,IUNT,MT,1,PSYDAT)
      WRITE(12,30)(PSYDAT(I),I=1,9)
      END DO
  30  FORMAT(1X,20(F8.4,1X))
  15  FORMAT(8X,'**',A,'**')
      END

      SUBROUTINE PSYCH(TIME,IUNITS,MODE,WBMD,PSYDAT)


      INTEGER WBMD,WBMODE
      DIMENSION PSYDAT(9)
      DATA RA/287.055/, PAOLD/0./
C
C  MODES 1-5:
C  THESE MODES TAKE AS INPUT PATM (IN ATMOSPHERES), A DRY BULB
C  TEMP., AND ONE OTHER PROPERTY: WET BULB TEMP., REL.HUMIDITY(FRACTION),       
C  DEW PT.TEMP.,HUMIDITY RATIO, OR ENTHALPY, DEPENDING ON MODE.
C  OUTPUTS ARE HUMIDITY RATIO (OR REL.HUMIDITY IN MODE 4), WET BULB TEMP.,      
C  ENTHALPY (OR REL.HUMIDITY IN MODE 5), MIXTURE DENSITY, AND DRY
C  AIR DENSITY.
C
C  MODE 6:
C  MODE 6 TAKES AS AN INPUT PATM, HUMIDITY RATIO AND ENTHALPY, AND
C  RETURNS ALL OF THE OTHER PROPERTIES.
C
C  THE WET BULB TEMPERATURE IS ONLY CALCULATED IF WBMODE EQUALS ONE.
C
C  TEMPERATURES ARE IN CELSIUS (IUNITS=1) OR FAHRENHEIT (IUNITS=2).
C  ENTHALPY IS IN KJ/KG (IUNITS=1) OR BTU/LBM (IUNITS=2), AND
C  DENSITY IS IN KG/M**3 (IUNITS=1) OR LBM/FT**3 (IUNITS=2).
C  THE REFERENCE STATES FOR ENTHALPIES ARE:
C         HAIR=0.0 AT 0. DEG C AND 0. DEG F
C         HW(LIQUID)=0.0 AT 0. DEG C AND 32. DEG F
C
C  THE PSYDAT ARRAY CONTAINS THE MOIST AIR PROPERTIES.
C
C  ASSIGN LOCAL VARIABLES:
      WBMODE = WBMD
      PATM = PSYDAT(1)
      TDB =  PSYDAT(2)
      TWB =  PSYDAT(3)
      TWBOLD=-9.999E20
      RH =   PSYDAT(4)
      TDP =  PSYDAT(5)
      W =    PSYDAT(6)
      H =    PSYDAT(7)

C
C  UNIT CONVERSIONS
C
      IF(IUNITS .EQ. 2) THEN
        TDB = (TDB - 32.)/1.8
        TWB = (TWB - 32.)/1.8
        TDP = (TDP - 32.)/1.8
        H   = (H - 7.687)/0.43002
      ENDIF
C-- CHECK THAT THE TOTAL PRESSURE IS WITHIN THE IDEAL GAS RANGE.
      IF (PATM .GT. 5.0) THEN
          STOP
      END IF
      IF (PATM .LE. 0.0) THEN

         STOP
      END IF
C
C   FOR MODE 6, CHECK THAT THE ENTHALPY IS GREATER THAN THE
C   SATURATION ENTHALPY (MINIMUM) FOR THE GIVEN HUMIDITY RATIO
C   AND THAT THE HUMIDITY RATIO IS GREATER THAN 0.
C
      IF (MODE .EQ. 6) THEN
            IF (W .LT. 0.) THEN
               W = 0.
               GO TO 99
            ELSE IF (W .EQ. 0.) THEN
               GO TO 99
            END IF
         CALL DEWPT (PATM,W,TDP)
         HDP = TDP + W*(2501. + 1.805*TDP)
         HMIN = HDP
           IF (H .LT. HMIN) THEN

               H = HMIN
           END IF
99       TDB = (H - 2501.*W)/(1. + 1.805*W)
       END IF
C
C  FIND SATURATION PRESSURE OF WATER AT WET BULB, DRY BULB, OR
C  DEW POINT TEMPERATURE.
C
      CALL SAT (TDB,PSATDB)
      GOTO (1,2,3,2,2,2) ,MODE
C-- CHECK FOR IMPOSSIBLE WET BULB TEMPERATURES AND CORRECT THEM
C   IF POSSIBLE.
1     IF (TWB .GT. TDB .AND. PSATDB .GE. PATM) THEN

         STOP
      ELSE IF (TWB .GT. TDB) THEN
         TWB = TDB
         PSAT = PSATDB
         GOTO 5
      END IF
      CALL SAT (TWB,PSAT)
C-- ERROR: IF PSATWB IS GREATER THAN PATM.
         IF (PSAT .GE. PATM) THEN

             STOP
         END IF
      GOTO 5
2     PSAT = PSATDB
      GOTO 5
C-- CHECK FOR IMPOSSIBLE DEW POINT TEMPERATURES AND
C   CORRECT THEM IF POSSIBLE.
3     IF (TDP .GT. TDB .AND. PSATDB .GE. PATM) THEN

         STOP
      ELSE IF (TDP .GT. TDB) THEN
         TDP=TDB
      END IF
      CALL SAT(TDP,PSAT)
C-- ERROR: IF PSATDP IS GREATER THAN PATM.
      IF (PSAT .GE. PATM) THEN

         STOP
      END IF
5     CONTINUE
C
C   CALCULATE HUMIDITY RATIO AND WET BULB TEMPERATURE
C
      GO TO (10,20,30,40,50,60), MODE
C
C  MODE 1 -- DRY BULB AND WET BULB SUPPLIED
C
10    IF (TWB .LE. 0.) THEN
        P = PSAT - 5.704E-4*(TDB-TWB)*PATM
        W = .62198 * P/(PATM-P)
      ELSE
        WSAT = .62198 * PSAT/(PATM-PSAT)
        W = WSAT - (TDB-TWB)*(0.24 + .441*WSAT)/(597.31
     .      + 0.441*TDB - TWB)
      ENDIF
      IF (W .LT. 0.0) THEN

         W = 0.0
         H = TDB
         TWBOLD = TWB
         WBMODE = 1
         PSAT = PSATDB
         GO TO 60
      END IF
      H = TDB + W*(2501. + 1.805*TDB)
      GO TO 100
C
C  MODE 2 -- DRY BULB AND RELATIVE HUMIDITY SUPPLIED
C
20    IF (RH .LT. 0.) THEN

          RH=0.0
      ELSE IF (PSAT .GE. PATM) THEN
          RHMAX = PATM/PSAT
          IF (RH .GE. (.99*RHMAX)) THEN

             STOP
          END IF
      ELSE IF (RH .GT. 1.) THEN

          RH=1.0
      END IF
      W = .62198 * PSAT*RH/(PATM-PSAT*RH)
      GO TO 40
C
C  MODE 3 -- DRY BULB AND DEW POINT SUPPLIED
C
30    W = .62198 * PSAT/(PATM-PSAT)
C
C  FIND ENTHALPY FOR MODES 2 - 4
C
40    IF (PSATDB .LT. PATM) THEN
         WMAX = .62198 * PSATDB/(PATM-PSATDB)
         IF (W .GT. WMAX) THEN

             W = WMAX
         END IF
      END IF
      IF (W .LT. 0.0) THEN

          W = 0.
      END IF
      H = TDB + W*(2501. + 1.805*TDB)
      GO TO 60
C
C  MODE 5 -- DRY BULB AND ENTHALPY SUPPLIED
C
50    IF (PSATDB .LT. PATM) THEN
         WMAX = .62198 * PSAT/(PATM-PSAT)
         HMAX = TDB + WMAX*(2501. + 1.805*TDB)
         IF (H .GT. HMAX) THEN

             H = HMAX
         END IF
      END IF
      HMIN = TDB
      IF (H .LT. HMIN) THEN

          H = HMIN
      END IF
        W = (H-TDB)/(2501.+1.805*TDB)
C
C  FIND WET BULB TEMPERATURE FOR MODES 2 - 6 IF WBMODE EQUALS 1.
C
60    IF (WBMODE .NE. 1) THEN
         TWB=TDB
         GO TO 100
      END IF
      DPRESS = ABS(1.-PATM)
C  THE FOLLOWING CORRELATION IS FOR 1 ATMOSPHERES TOTAL PRESSURE.
C  IF OUTSIDE THE CORRELATION RANGE, THE CORRELATION IS USED FOR
C  THE INITIAL GUESS IN THE ITERATIVE METHOD.
      IF (H .GT. 0. .AND. H .LT. 2000.) THEN
          Y = ALOG(H*.43002+7.687)
          IF (H .LE. 9.473) THEN
             TWB=-17.4422+1.9356*Y+.7556*Y**2+.5406*Y**3
          ELSE IF (H .GT. 9.473) THEN
             TWB=-.6008-22.04556*Y+11.4356*Y**2-.97667*Y**3
          END IF
       ELSE
          TWB = 9.99999E25
       END IF
C  USE A NEWTON'S ITERATIVE METHOD TO FIND THE WET BULB
C  TEMPERATURE.
       IF (DPRESS.GT..001 .OR. H.LE.0. .OR. H.GT.275.) THEN
          ITEST=0
          IF (ABS(PATM-PAOLD) .GT. 1.0E-10) THEN
          CALL BOIL(PATM,TBOIL)
          PAOLD = PATM
          END IF
C-- INITIAL GUESS
          TWBNEW=AMIN1(TWB,(TBOIL-0.1),TDB)
70        IF (TWBNEW .GE. (TBOIL-0.09)) TWBNEW=TBOIL-0.1
          CALL SAT(TWBNEW,PSAT)
          WSSTAR=.62198*PSAT/(PATM-PSAT)
          IF(MODE .EQ. 5) THEN
            W = (H-TDB)/(2501.+1.805*TDB)
          END IF
          WNEW=((2501.-2.381*TWBNEW)*WSSTAR-(TDB-TWBNEW))/
     .         (2501.+1.805*TDB-4.186*TWBNEW)
          ERR = W - WNEW
          IF (ABS(ERR) .LE. (.01*W)) GO TO 75
          IF (W .EQ. 0.) THEN
             IF (ABS(ERR) .LE. .0001) GO TO 75
          END IF
          ITEST = ITEST + 1
          IF (ITEST .GE. 25) GO TO 75
C--FIND THE SLOPE OF THE ERROR FUNCTION
          TSLOPE = 0.999*TWBNEW
          IF (TWBNEW .EQ. 0.) TSLOPE = -.005
          CALL SAT(TSLOPE,PSLOPE)
          WSSLP=.62198*PSLOPE/(PATM-PSLOPE)
          WSLOPE=((2501.-2.381*TSLOPE)*WSSLP-(TDB-TSLOPE))/
     .         (2501.+1.805*TDB-4.186*TSLOPE)
          ERRSLP = W - WSLOPE
          DERRDT = (ERRSLP - ERR)/(TSLOPE - TWBNEW)
          TWBNEW = TWBNEW - ERR/DERRDT
          GO TO 70
75        TWB = TWBNEW
      END IF
      IF (TWB .LT. TWBOLD) TWB = TWBOLD
100   CONTINUE
C
C  FIND RELATIVE HUMIDITY, DEW POINT, MIXTURE DENSITY,
C  AND DRY AIR DENSITY
C
      PV  = PATM*W/(.62198+W)
      IF(MODE .NE. 2) RH = PV/PSATDB
      IF(MODE .NE. 3 .AND. PV .GT. 0.) THEN
        CALL DEWPT(PATM,W,TDP)
      ELSE IF (MODE .NE. 3 .AND. PV .LE. 0) THEN
C  FOR DRY AIR, THERE IS NO DEW POINT TEMPERATURE
        TDP=-9.99999E25
      ENDIF
      SPCVOL = RA*(TDB+273.15)/(PATM*101325)*(1+1.6078*W)
      RHOWA = 1/SPCVOL
      RHOWM = RHOWA*(1+W)
C
C  CONVERT OUTPUTS TO APPROPRIATE UNITS
C
      IF(IUNITS .EQ. 2) THEN
        H   = H*0.43002+7.68
        TDB = 1.8*TDB + 32.
        TWB = 1.8*TWB + 32.
        IF (TDP .GT. -9.99999E24) THEN
         TDP = 1.8*TDP + 32.
        END IF
        RHOWM = RHOWM/16.02
        RHOWA = RHOWA/16.02
      ENDIF
C SET OUTPUTS
      IF(MODE .EQ. 6) PSYDAT(2) = TDB
      PSYDAT(3) = TWB
      PSYDAT(4) = RH
      PSYDAT(5) = TDP
      PSYDAT(6) = W
      PSYDAT(7) = H
      PSYDAT(8) = RHOWM
      PSYDAT(9) = RHOWA
C
      RETURN
C

      END
C
C**********************************************************************
C
C  SUBROUTINE FOR FINDING SATURATION PRESSURE OF WATER AT A GIVEN
C  TEMPERATURE
C
      SUBROUTINE SAT(TIN,PSAT)


C
C  THE FOLLOWING CORRELATION FOR THE SATURATION PRESSURE OF
C  WATER VAPOR (IN PASCALS) AS A FUNCTION OF TEMPERATURE IS
C  TAKEN FROM THE 1985 ASHRAE FUNDAMENTALS HANDBOOK (SI).
C
      DATA C1/-5674.5359/,C2/6.3925247/,C3/-0.9677843E-2/
      DATA C4/0.62215701E-6/,C5/0.20747825E-8/,C6/-0.9484024E-12/
      DATA C7/4.1635019/,C8/-5800.2206/,C9/1.3914993/
      DATA C10/-0.048640239/,C11/0.41764768E-4/,C12/-0.14452093E-7/
      DATA C13/6.5459673/
C
      T = TIN + 273.15
      IF (T .LE. 0.) THEN

         STOP
      END IF
C
C  SATURATION PRESSURE OVER ICE (-100 C TO 0 C)
C
      IF (T .LT. 273.15) THEN
         PSAT=EXP(C1/T+C2+C3*T+C4*T**2+C5*T**3+C6*T**4+C7*ALOG(T))
C
C  SATURATION PRESSURE OVER LIQUID WATER (0 C TO 200 C)
C
      ELSE IF (T .GE. 273.15) THEN
         PSAT=EXP(C8/T+C9+C10*T+C11*T**2+C12*T**3+C13*ALOG(T))
      END IF
C
C  TEMPERATURE OUT OF THE RANGE USED FOR THE CORRELATION
C
      IF (T .LT. 173.15 .OR. T .GT. 473.15) THEN
          STOP
      END IF
C  CONVERT PRESSURE FROM PASCALS TO ATMOSPHERES
      PSAT = PSAT/101325
      RETURN
C
      END
C
C******************************************************************
C
C  SUBROUTINE FOR FINDING THE DEW POINT TEMPERATURE GIVEN
C  THE HUMIDITY RATIO.  THE CORRELATION IS FROM THE 1981
C  ASHRAE FUNDAMETALS HANDBOOK.  THE DEW POINT TEMPERATURE
C  IS IN DEGREES C AND ATMOSPHERIC PRESSURE IN PASCALS.
C
      SUBROUTINE DEWPT(PATM,W,TDP)



C
      PV  = PATM*W/(.62198+W)
      Y = ALOG(1.013E05*PV)
      TDP = -35.957 - 1.8726*Y + 1.1689*Y*Y
      IF(TDP .LT. 0.) TDP = -60.45 + 7.0322*Y + 0.3700*Y*Y
C
      IF (TDP.GT.70. .OR. TDP.LT.-60.) THEN
         STOP
      END IF
      RETURN
      END
C
C******************************************************************
C
C SUBROUTINE FOR FINDING THE BOILING TEMPERATURE OF WATER GIVEN
C THE TOTAL PRESSURE.  A NEWTON'S METHOD IS USED WITH THE
C SATURATED WATER VAPOR PRESSURE CORRELATION FROM THE 1985 ASHRAE
C FUNDAMENTALS HANDBOOK.
C
      SUBROUTINE BOIL(PATM,TBOIL)
C
      DATA C1/-5674.5359/,C2/6.3925247/,C3/-0.9677843E-2/
      DATA C4/0.62215701E-6/,C5/0.20747825E-8/,C6/-0.9484024E-12/
      DATA C7/4.1635019/,C8/-5800.2206/,C9/1.3914993/
      DATA C10/-0.048640239/,C11/0.41764768E-4/,C12/-0.14452093E-7/
      DATA C13/6.5459673/
C
      PBOIL = PATM*101325
      ITEST = 0
      IF (PBOIL .LT. 611.21) GO TO 100
C--USING ASHRAE CORRELATION FOR 0 C TO 200 C.
C  FIRST GUESS TBOIL EQUALS 100 C.
      T1 = 100 + 273.15
10    ZZ=C8/T1+C9+C10*T1+C11*T1**2+C12*T1**3+C13*ALOG(T1)
      P1=EXP(ZZ)
      ERR = PBOIL - P1
      IF (ABS(ERR) .LE. (.01*PBOIL)) THEN
          TBOIL = T1 - 273.15
          RETURN
      END IF
      ITEST = ITEST +1
      IF (ITEST .GT. 100) STOP
      DERRDT=-(P1)*(-C8/T1**2+C10+2*C11*T1+3*C12*T1**2+C13/T1)
      TOLD = T1
      T1 = T1 - ERR/DERRDT
      GO TO 10
C
100   CONTINUE
C--USING ASHRAE CORRELATION FOR -100 C TO 0 C.
C  FIRST GUESS TBOIL EQUALS 0 C.
      T1 = 273.15
110   ZZ=C1/T1+C2+C3*T1+C4*T1**2+C5*T1**3+C6*T1**4+C7*ALOG(T1)
      P1=EXP(ZZ)
      ERR = PBOIL - P1
      IF (ABS(ERR) .LE. (.01*PBOIL)) THEN
          TBOIL = T1 - 273.15
          RETURN
      END IF
      ITEST = ITEST +1
      DERRDT=-P1*(-C1/T1**2+C3+2*C4*T1+3*C5*T1**2+4*C6*T1**3+C7/T1)
      TOLD = T1
      T1 = T1 - ERR/DERRDT
      GO TO 110
      END
