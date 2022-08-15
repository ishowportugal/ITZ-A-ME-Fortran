      PROGRAM CURVE
      CHARACTER INAM*6,STA*3
      WRITE(*,9000)
9000  FORMAT(' Quel type de regression',/,
     *' 1 : lineaire y=a+b*x',/,' 2 : exponentiele y=a*exp(b*x)'
     *,/,' 3 : logarithmique y=a+b*ln(x)',/,
     *' 4 : puissance y=a*x**b ?')
      READ(*,*)I0
      WRITE(*,9001)
9001  FORMAT(' Le nom du fichier contenant les points?')
      READ(*,9002)INAM
9002  FORMAT(A6)
      STA='OLD'
      CALL FICHR(100,INAM,STA,1)
      READ(100,*)N
      DO 50 I=1,N
      READ(100,*)X,Y
      GOTO (40,10,20,30)I0
 10   Y=ALOG(Y)
      GOTO 40
 20   X=ALOG(X)
      GOTO 40
 30   X=ALOG(X)
      Y=ALOG(Y)
 40   SX1=SX1+X
      SX2=SX2+X*X
      SY1=SY1+Y
      SY2=SY2+Y*Y
 50   SXY=SXY+X*Y
      CALL FICHR(100,INAM,STA,2)
C
      DEL=FLOAT(N)*SX2-SX1*SX1
      A=(SY1*SX2-SXY*SX1)/DEL
      B=(FLOAT(N)*SXY-SX1*SY1)/DEL
      DEL=FLOAT(N)*SY2-SY1*SY1
      BP=(FLOAT(N)*SXY-SY1*SX1)/DEL
      R=B*BP
      GOTO(70,60,70,60)I0
 60   A=EXP(A)
 70   WRITE(*,9005)A,B,R
9005  FORMAT(/' a = ',E12.5,/,' b = ',E12.5,/,
     *'r**2 = ',E12.5)
      END
