      FUNCTION DICHO(XMAX,XMIN)
C
C**********************************************************
C CALCUL DE LA RACINE D'UNE FONCTION (FNEW) A L'AIDE DE LA
C METHODE DE DICHOTOMIE.
C ATHANASSIOS ARGIRIOU, 6 - 5 - 1990.
C***********************************************************
C
      EPS=XMAX*1.E-05
      X2=XMAX
      X1=XMIN
      Y1=FNEW(XMIN)
 5    X=(X1+X2)/2.
      Y=FNEW(X)
      R=Y*Y1
      IF(R)10,10,20
 10   X2=X
      Y2=Y
      GOTO 30
 20   X1=X
      Y1=Y
 30   EP=ABS(X2-X1)
      IF(EP-EPS)40,5,5
 40   DICHO=(X1+X2)/2.
      IF (DICHO-XMIN.LT.EPS.OR.XMAX-DICHO.LT.EPS)THEN
      WRITE(*,9000)
 9000 FORMAT(1X,'DICHO: ATTENTION X NEAR THE BOUNDARIES')
      ENDIF
      RETURN
      END
C
      FUNCTION RAPHNEW(X0)
C
C*********************************************************
C CALCUL DE LA RACINE D'UNE EQUATION A L'AIDE DE LA METHODE
C DE NEWTON - RAPHSON
C FNEW(X) LA FONCTION EN QUESTION, DNEW(X) 1ERE DERIVEE DE FNEW(X)
C ATHANASSIOS ARGIRIOU, 6 - 5 - 1990
C**********************************************************
C
      EPS=X0*1.E-05
 5    F=FNEW(X0)
      F1=DNEW(X0)
      RAPHNEW=X0-F/F1
      IF(ABS(X0-RAPHNEW).LT.EPS)GOTO 10
      X0=RAPHNEW
      GOTO 5
 10   RETURN
      END
C
C
      FUNCTION SIMPN(X1,X2,DX,F)
C***************************************************************
C INTEGRATION PAR LA METHODE DE SIMPSON
C X1,X2 : BORNES D'INTEGRATION
C DX    : PAS D'INTEGRATION
C F     : FONCTION INTEGREE
C ATHANASSIOS ARGIRIOU, 16 - 6 - 1990
C****************************************************************
C
      I=INT((X2-X1)/DX)
      IF(REAL(I)/2.-REAL(I/2).GT..1) I=I+1
      DX=(X2-X1)/REAL(I)
      X=X1
      I1=I+1
      DO 10 J=1,I1
      IF(XF.GT.3.) THEN
      XF=2.
      ELSE
      XF=4.
      ENDIF
      IF(J.EQ.1.OR.J.EQ.I1) XF=1.
      SIMPN=SIMPN+F(X)*XF
 10   X=X1+DX*REAL(J)
      SIMPN=DX*SIMPN/3.
      RETURN
      END
C
C
      SUBROUTINE BESJ(X,N,BJ,D,IER)
C******************************************************************
C CALCULE LA FONCTION DE BESSEL J POUR UN ORDRE ET UN ARGUMENT DONNES
C
C	X    - ARGUMENT DE LA FONCTION
C	N    - ORDRE DE LA FONCTION
C	BJ   - VALEUR DE LA FONCTION DE BESSEL
C	D    - PRECISION DEMANDEE
C	IER  - CODE D'ERREUR
C	        IER=0, PAS D'ERREUR
C		IER=1, N EST NEGATIF
C		IER=2, X EST NEGATIF OU ZERO
C		IER=3, LA PRECISION REQUISE N'EST PAS ATTEINTE
C		IER=4, VALEUR DE N NON COMPATIBLE AVEC CELLE DE X
C		       (VOIR REMARQUES)
C
C	REMARQUES
C       - N DOIT ETRE SURERIEUR OU EGAL A ZERO, MAIS IL DOIT ETRE
C            INFERIEUR DE:
C            20+10*X-X**2/3     POUR X INFERIEUR OU EGAL A 15
C            90+X/2             POUR X SUPERIEUR A 15
C
C ATHANASSIOS ARGIRIOU, 16 - 6 - 1990
C**********************************************************************
C
      BJ=0.
      IF (N) 10,20,20
 10   IER=1
      RETURN
 20   IF (X) 30,30,31
 30   IER=2
      RETURN
 31   IF(X-15.) 32,32,34
 32   NTEST=20.+10.*X-X**2/3.
      GO TO 36
 34   NTEST=90.+X/2.
 36   IF(N-NTEST) 40,38,38
 38   IER=4
      RETURN
 40   IER=0
      N1=N+1
      BPREV=.0
      IF(X-5.) 50,60,60
 50   MA=X+6.
      GO TO 70
 60   MA=1.4*X+60./X
 70   MB=N+IFIX(X)/4+2
      MZERO=MAX0(MA,MB)
      MMAX=NTEST
 100  DO 190 M=MZERO,MMAX,3
      FM1=1.E-28
      FM=.0
      ALPHA=.0
      IF(M-(M/2)*2) 120,110,120
 110  JT=-1
      GO TO 130
 120  JT=1
 130  M2=M-2
      DO 160 K=1,M2
      MK=M-K
      BMK=2.*FLOAT(MK)*FM1/X-FM
      FM=FM1
      FM1=BMK
      IF(MK-N-1) 150,140,150
 140  BJ=BMK
 150  JT=-JT
      S=1+JT
 160  ALPHA=ALPHA+BMK*S
      BMK=2.*FM1/X-FM
      IF(N) 180,170,180
 170  BJ=BMK
 180  ALPHA=ALPHA+BMK
      BJ=BJ/ALPHA
      IF(ABS(BJ-BPREV)-ABS(D*BJ)) 200,200,190
 190  BPREV=BJ
      IER=3
 200  RETURN
      END
C
C
      SUBROUTINE I0(X,RIC)
C*******************************************************************
C CALCULE LA FONCTION DE BESSEL MODIFIEE D'ORDRE 0 ET DE 1ER ESPECE
C X : ARGUMENT D'ENTREE - RIC : ARGUMENT DE SORTIE
C REMARQUES: UNE VALEUR DE X TROP GRANDE PEUT PROVOQUER UN OVERFLOW
C DE LA FONCTION EXP
C ATHANASSIOS ARGIRIOU, 18 - 6 - 1990
C*******************************************************************
C
      RIC=ABS(X)
      IF(RIC-3.75)1,1,2
 1    Z=(X/3.75)**2
      RIC=(((((4.5813E-3*Z+3.60768E-2)*Z+2.659732E-1)*Z+1.206749)
     $*Z+3.089942)*Z+3.515623)*Z+1.
      RETURN
 2    Z=3.75/RIC
      RIC=EXP(RIC)/SQRT(RIC)*((((((((3.92377E-3*Z-1.647633E-2)*Z
     $+2.635537E-2)*Z-2.057706E-2)*Z+9.16281E-3)*Z-1.57565E-3)*Z
     $+2.25319E-3)*Z-1.328592E-2)*Z+3.989423E-1)
      RETURN
      END
C
C
      SUBROUTINE INUE(X,N,XI,RI)
C*******************************************************************
C COMPUTE THE MODIFIED BESSEL FUNCTIONS I FOR ORDERS 1 TO N
C 	- X 	GIVEN ARGUMENT OF THE BESSEL FUNCTIONS I
C	- N	GIVEN MAXIMUM ORDER OF BESSEL FUNCTIONS I
C	- XI	GIVEN VALUE OF BESSEL FUNCTION I OF ORDER ZERO FOR 
C		ARGUMENT X
C	- RI	RESULTANT VECTOR OF DIMENSION N, CONTAINING THE VALUES 
C		OF THE FUNCTIONS I FOR ORDERS 1 TO N
C ATHANASSIOS ARGIRIOU, 18 - 6 - 1990
C**********************************************************************
      DIMENSION RI(N)
      IF (N) 10,10,1
    1 FN=N+N
      G1=X/FN
      IF(ABS(X)-5.E-04)4,4,2
   2  A0=1.
      A1=0.
      B0=0.
      B1=1.
      F1=FN
  3   F1=F1+2.
      AW=F1/ABS(X)
      A=AW*A1+A0
      B=AW*B1+B0
      A0=A1
      B0=B1
      A1=A
      B1=B
      G0=G1
      G1=A/B
      IF(ABS(G1-G0)/G1-1.E-06)4,4,3
 4    IF (X) 5,6,6
 5    G1=-G1
 6    K=N
 7    G1=X/(FN+X*G1)
      RI(K)=G1
      FN=FN-2.
      K=K-1
      IF (K) 8,8,7
 8    FI=XI
      DO 9 I=1,N
      FI=FI*RI(I)
 9    RI(I)=FI
 10   RETURN
      END
