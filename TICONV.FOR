      SUBROUTINE ticonv(TIDEC,IDAY,IHOUR,MIN)
c******************************************************************
c     READS THE TIME IN DAYS AND DECIMALS OF DAYS AND CONVERTS IT
C     IN DAY HOUR MINUTE
C
C     INPUT PARAMETERS
C     TIDEC = TIME IN DAYS AND DECIMALS OF DAYS
C
C     OUTPUT PARAMETERS
C     IDAY : JULIAN DAY
C     IHOUR: HOUR
C     MIN  : MINUTE
C
C     ATHANASSIOS ARGIRIOU, 13 MARCH 1994
C********************************************************************
C
C CALCULATION OF THE DAY
      IDAY=INT(TIDEC)
C CALCULATION OF THE HOUR
      RES=(TIDEC-FLOAT(IDAY))*1440.
      IHOUR=INT(RES/60.)
C CALCULATION OF THE MINUTES
      MIN=INT(RES-FLOAT(IHOUR)*60.)
C
      RETURN
      END
      