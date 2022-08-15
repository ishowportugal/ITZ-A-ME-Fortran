      PROGRAM SOLTEM
C
C THIS PROGRAM CALCULATES THE SOIL TEMPERATURE FOR A GIVEN DEPTH
C AND FOR A GIVEN DAY OF THE YEAR
C AS : amplitude of ground surface temperature
C TXT: ground temperature at depth X and day T
C X  : depth in feet
C TG : average annual ground temperature
C TDS: ground diffusivity (ft2/day)
C T  : day of the year (1 - 365)
C PC : day of the year of the lowest average ground temperature
C 
      PI=4.*ATAN(1.)
      WRITE(*,9000)
 9000 FORMAT(' ENTER AS,X,TG,TDS,T,PC',/)
      READ(*,*)AS,X,TG,TDS,T,PC
      X=X/.3048
      TDS=TDS/.0929*86400.
      TXT=TG-(AS*EXP(-X*(PI/(365.*TDS))**.5)*COS((2.*PI/365.)*
     &(T-PC-(X/2.)*((365./(PI*TDS))**.5))))
      WRITE(*,9100)T,X,TXT
 9100 FORMAT(' THE SOIL TEMPERATURE DURING THE',F4.0,' DAY OF THE'
     &' YEAR, AT A DEPTH OF ',F3.1,' IS ',F4.1,' C',/)
      END
