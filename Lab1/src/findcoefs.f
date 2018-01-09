C -----------------------------------------------------------------------
C  Program to run a least squares minimization to standard star 
C  photometry to find transformation coefficients: 
C  instrumental --> standard
C
C  ISYA 2004                               (c) Ignasi Ribas (IEEC, Spain)
C -----------------------------------------------------------------------
      PROGRAM FINDCOEFS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAX=1000,N=2)
      DIMENSION YINST(MAX),BINST(MAX)
      DIMENSION XMA(MAX,N),XDEP(MAX)
      DIMENSION XV(MAX),XCV(N,N),XCR(N,N),XR(N),XERR(MAX)
      DIMENSION YMA(MAX,N),YDEP(MAX)
      DIMENSION YV(MAX),YCV(N,N),YCR(N,N),YR(N),YERR(MAX)

      Open(1,file='standards.lst',status='old')
      Do 1 k=1,MAX
       Read(1,*,end=2) yinst(k),xerr(k),binst(k),sigb,ystd,ydep(k)
       xma(k,1)=1.d0 
       yma(k,1)=1.d0
       xma(k,2)=binst(k)-yinst(k)
       yma(k,2)=binst(k)-yinst(k)
       yerr(k)=dsqrt(xerr(k)*xerr(k)+sigb*sigb)
       xdep(k)=ystd-yinst(k)
  1   Continue
  2   Close(1)
      ns=k-1

      CALL MINQ(XMA,XDEP,XERR,NS,XR,XCV,XCR,XV,XSIG,XR1)
      CALL MINQ(YMA,YDEP,YERR,NS,YR,YCV,YCR,YV,YSIG,YR1)

      Open(2,file='coefficients.dat')
      Write(2,50) 
  50  Format(20x,'FIT RESULTS',/,20x,'-----------')
      Write(2,*) 
      Write(2,100) XR(1),XR(2)
      Write(2,200) SQRT(XCV(1,1)),SQRT(XCV(2,2))
      Write(2,*) 
      Write(2,400) XSIG,XR1
      Write(2,*) 
      Write(2,*) 
      Write(2,300) YR(1),YR(2)
      Write(2,200) SQRT(YCV(1,1)),SQRT(YCV(2,2))
      Write(2,*) 
      Write(2,400) YSIG,YR1
 100  Format('y_std-y_inst =    ',f7.3,'  + ',f7.3,' * (b-y)_inst')
 200  Format('               +/-',f7.3,' +/-',f7.3)
 300  Format('   (b-y)_std =    ',f7.3,'  + ',f7.3,' * (b-y)_inst')
 400  Format(5x,'sigma = ',f7.3,5x,'correlation coeff = ',f7.3)
      Write(2,*)
      Write(2,500)
 500  Format('  y_inst','   O-C ','  sigma',' (b-y)_i','   O-C ',
     $ '  sigma')
      Do 40 i=1,ns
       calc1=xr(1)+xma(i,2)*xr(2)
       res1=xdep(i)-calc1
       calc2=yr(1)+yma(i,2)*yr(2)
       res2=ydep(i)-calc2
       Write(2,'(2(f8.3,2f7.3))') yinst(i),-res1,xerr(i),yma(i,2),res2,
     $  yerr(i)
 40   Continue
      Close(2)

C      WRITE(*,900)
C      Do 10 i=1,n
C        WRITE(*,800) (XCV(i,j),j=1,n)
C 10   Continue
C      WRITE(*,901)
C      Do 20 i=1,n
C        WRITE(*,800) (YCV(i,j),j=1,n)
C 20   Continue

 800  Format(8x,4f9.5)
 900  Format(/,/,9x,'Covariance matrix 1')
 901  Format(/,/,9x,'Covariance matrix 2')

      Stop
      End
*
*     -------------
*     LEAST SQUARES
*     -------------
*
      SUBROUTINE MINQ (A,ALTINCR,ERRY,NS,R,CV,CR,V,SIGMA,R2)
      PARAMETER (MAX=1000,N=2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(MAX,N),ALTINCR(MAX),ERRY(MAX)
      DIMENSION P(MAX,MAX),ATP(N,MAX),AT(N,MAX),C(N,N)
      DIMENSION V(MAX),CV(N,N),CR(N,N)
      DIMENSION D(N,MAX),E(MAX),F(MAX),G(MAX),R(N)
      DIMENSION LL(N),MM(N)
*
*     -------------------------------
*     INITIALIZATION OF WEIGHT MATRIX 
*     -------------------------------
*
      PT=0.
      DO 5 I=1,NS
         DO 6 J=1,NS
            IF (I.EQ.J) THEN
               P(I,J)=1.D0/ERRY(I)
            ELSE
               P(I,J)=0.D0
            ENDIF
6        CONTINUE
         PT=PT+P(I,I)
5     CONTINUE
      DO 7 I=1,NS
7     P(I,I)=NS*(P(I,I)/PT)
*
*     ------------------------------------
*     CALCULATION OF THE TRANSPOSED MATRIX
*     -------------------------------------
*
      DO 10 I=1,NS
         DO 15 J=1,N
           AT(J,I)=A(I,J)
15       CONTINUE
10    CONTINUE
*
*     ----------------------------------
*     PRODUCT OF TRANSPOSED TIMES WEIGHT
*     ----------------------------------
*
      DO 20 I=1,N
         DO 25 J=1,NS
            ATP(I,J)=0.D0
           DO 30 K=1,NS
               ATP(I,J)=ATP(I,J)+AT(I,K)*P(K,J)
30         CONTINUE
25       CONTINUE
20    CONTINUE
*
*     --------------------------
*     PRODUCT WITH SYSTEM MATRIX
*     --------------------------
*
      DO 35 I=1,N
         DO 40 J=1,N
           C(I,J)=0.D0
            DO 45 K=1,NS
              C(I,J)=C(I,J)+ATP(I,K)*A(K,J)
45          CONTINUE
40       CONTINUE
35    CONTINUE
*
*     -------
*     INVERSE
*     -------
*
      CALL INVERT (C,N,DET,LL,MM)
*
*     -----------------------------------
*     PRODUCT OF INVERSE TIMES TRANSPOSED
*     -----------------------------------
*
      DO 50 I=1,N
         DO 55 J=1,NS
            D(I,J)=0.D0
            DO 60 K=1,N
              D(I,J)=D(I,J)+C(I,K)*AT(K,J)
60          CONTINUE
55       CONTINUE
50    CONTINUE
*
*     -----------------------------------
*     PRODUCT OF WEIGHT TIMES DATA MATRIX
*     -----------------------------------
*
      em1=0.
      em2=0.
      DO 65 I=1,NS
         E(I)=0.D0
         em1=em1+ALTINCR(I)
         em2=em2+A(I,2)
         DO 70 J=1,NS
            E(I)=E(I)+P(I,J)*ALTINCR(J)
70    CONTINUE
65    CONTINUE
      em1=em1/NS
      em2=em2/NS
      dyt=0.
      dxt=0.
      dxy=0.
      Do 77 I=1,NS
      DXY=DXY+(P(I,I)/NS)*(A(I,2)-em2)*(ALTINCR(I)-em1)
      DXT=DXT+(P(I,I)/NS)*(A(I,2)-em2)*(A(I,2)-em2)
77    DYT=DYT+(P(I,I)/NS)*(ALTINCR(I)-em1)*(ALTINCR(I)-em1)
      R2=DXY/(SQRT(DXT*DYT))
*
*     -------------------------
*     SOLUTION BY LEAST SQUARES
*     -------------------------
*
      DO 75 I=1,N
         R(I)=0.D0
         DO 80 J=1,NS
            R(I)=R(I)+D(I,J)*E(J)
80       CONTINUE
75    CONTINUE
*
*     ------------------------------------
*     CALCULATION OF THE COVARIANCE MATRIX 
*     ------------------------------------
*
      DO 85 I=1,NS
         F(I)=0.D0
         DO 90 J=1,N
           F(I)=F(I)+A(I,J)*R(J)
90       CONTINUE
        V(I)=F(I)-ALTINCR(I)
85    CONTINUE
*
      XNUMSIG=0.D0
      DO 95 I=1,NS
         G(I)=0.D0
         DO 100 J=1,NS
            G(I)=G(I)+V(J)*P(J,I)
100      CONTINUE
         XNUMSIG=XNUMSIG+G(I)*V(I)
95    CONTINUE
      SIGMAQ=XNUMSIG/(NS-N)
      SIGMA=DSQRT(SIGMAQ)
*
      DO 105 I=1,N
         DO 110 J=1,N
            CV(I,J)=SIGMAQ*C(I,J)
110      CONTINUE
105   CONTINUE
*
*     ------------------
*     CORRELATION MATRIX
*     ------------------
*
      DO 115 I=1,N
         DO 120 J=1,N
            CR(I,J)=C(I,J)/DSQRT(C(I,I)*C(J,J))
120      CONTINUE
115   CONTINUE
*
      RETURN
      END
*
*     ----------------
*     MATRIX INVERSION
*     ----------------
*
      SUBROUTINE INVERT (A,N,D,LL,MM)
      DIMENSION A(2*N),LL(N),MM(N)
      DOUBLE PRECISION A,D,BIGA,HOLD,DABS
*     --------------------------
*     SEARCH FOR LARGEST ELEMENT
*     --------------------------
      D=1.0D0
      NK=-N
      DO 80 K=1,N
      NK=NK+N
      LL(K)=K
      MM(K)=K
      KK=NK+K
      BIGA=A(KK)
      DO 20 J=K,N
      IZ=N*(J-1)
      DO 20 I=K,N
      IJ=IZ+I
10    IF (DABS(BIGA)-DABS(A(IJ))) 15,20,20
15    BIGA=A(IJ)
      LL(K)=I
      MM(K)=J
20    CONTINUE
*     ----------------
*     INTERCHANGE ROWS
*     ----------------
      J=LL(K)
      IF(J-K) 35,35,25
25    KI=K-N
      DO 30 I=1,N
      KI=KI+N
      HOLD=-A(KI)
      JI=KI-K+J
      A(KI)=A(JI)
30    A(JI) =HOLD
*     -------------------
*     INTERCHANGE COLUMNS
*     -------------------
35    I=MM(K)
      IF(I-K) 45,45,38
38    JP=N*(I-1)
      DO 40 J=1,N
      JK=NK+J
      JI=JP+J
      HOLD=-A(JK)
      A(JK)=A(JI)
40    A(JI) =HOLD
*     -------------------------------------------------------
*     DIVIDE COLUMN BY MINUS PIVOT (VALUE OF PIVOT ELEMENT IS
*     CONTAINED IN BIGA)
*     -------------------------------------------------------
45    IF(BIGA) 48,46,48
46    D=0.0D0
      RETURN
48    DO 55 I=1,N
      IF(I-K) 50,55,50
50    IK=NK+I
      A(IK)=A(IK)/(-BIGA)
55    CONTINUE
*     -------------
*     REDUCE MATRIX
*     -------------
      DO 65 I=1,N
      IK=NK+I
      HOLD=A(IK)
      IJ=I-N
      DO 65 J=1,N
      IJ=IJ+N
      IF(I-K) 60,65,60
60    IF(J-K) 62,65,62
62    KJ=IJ-I+K
      A(IJ)=HOLD*A(KJ)+A(IJ)
65    CONTINUE
*     -------------------
*     DIVIDE ROW BY PIVOT
*     -------------------
      KJ=K-N
      DO 75 J=1,N
      KJ=KJ+N
      IF(J-K) 70,75,70
70    A(KJ)=A(KJ)/BIGA
75    CONTINUE
*     -----------------
*     PRODUCT OF PIVOTS
*     -----------------
      D=D*BIGA
*     ---------------------------
*     REPLACE PIVOT BY RECIPROCAL
*     ---------------------------
      A(KK)=1.0D0/BIGA
80    CONTINUE
*     --------------------------------
*     FINAL ROW AND COLUMN INTERCHANGE
*     --------------------------------
      K=N
100   K=(K-1)
      IF(K) 150,150,105
105   I=LL(K)
      IF(I-K) 120,120,108
108   JQ=N*(K-1)
      JR=N*(I-1)
      DO 110 J=1,N
      JK=JQ+J
      HOLD=A(JK)
      JI=JR+J
      A(JK)=-A(JI)
110   A(JI) =HOLD
120   J=MM(K)
      IF(J-K) 100,100,125
125   KI=K-N
      DO 130 I=1,N
      KI=KI+N
      HOLD=A(KI)
      JI=KI-K+J
      A(KI)=-A(JI)
130   A(JI) =HOLD
      GO TO 100
150   RETURN
      END
