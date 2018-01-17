      SUBROUTINE pymodels(feh,dist,eby,XA,nfit,imod,Vout,byout,YLout,
     & YTout,YGout,YMPout,YLMout,Rout)
!f2py double precision, intent(in) :: feh,dist,eby,XA
!f2py intent(in) :: nfit
!f2py integer, intent(in) :: imod
!     feh   metalicity
!     dist  distance (pc)
!     eby   reddening
!     XA    age (Myr)
!     nfit  name of the output file (30 characters maximum)
!f2py double precision, dimension(730),intent(out) :: Vout,byout,YLout,YTout,YGout,YMPout,YLMout,Rout
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION DLO0(35,1190,10),DLM0(35,10),DT0(35,1190,10)
      DIMENSION GG0(35,1190,10),XMAS0(35,1190,10),ZM(10)
      DIMENSION DDL0(35,1190,10),XMF(1190)
      DIMENSION BU1(35,1190),BU2(35),DDL(35,1190),DLM(35)
      DIMENSION DLMI(10),DTI(10),GGI(10),XMASI(10),DLOI(10),DDLI(10)
      DIMENSION DLAF(1190),XLMF(1190),GF(1190),DLF(1190),DTF(1190)
      DIMENSION XMAS(35,1190),GG(35,1190),DT(35,1190),DLO(35,1190)
      DIMENSION Zint(10)

      CHARACTER*30 NFIT

      xmbsun=4.75d0

!117   print*,'Which set of evolutionary models?'
!      print*,' 1.- Cassisi overshooting '
!      print*,' 2.- Cassisi standard '
!      print*,' 3.- Padova group '
!      read(*,'(i1)') imod
!      If(imod.gt.3.or.imod.lt.1) goto 117
!      print*,' '
!      print*,'Name of the output file'
!      print*,' '
!      read(*,'(a30)') nfit
 
 197  iso=2
C197  print*,' '
C     print*,'Compute evolutionary track (1) or isochrone (2)?'
C     print*,' '
C     read(*,'(i1)') iso
 
 198  Open(1,file=nfit,status='unknown')

! 199  print*,' '
!      print*,'Enter the metallicity [Fe/H]'
!      print*,' '
!      read(*,*) feh
      Zest=0.0198d0*10.d0**(feh)
      Zlog=DLOG10(Zest)

      If (iso.eq.1) then
        print*,' '
        print*,'Enter the value for the mass (ZAMS) in Msol'
        print*,' '
        read(*,*) XME
        XLME=DLOG10(XME)
        write(1,71) Zest,XME
 71     Format('# Z = ',f7.4,';  M = ',f7.3,' Msun')
        write(1,73)
 73     Format('# log L ','  log T ','  log g ','  age (Myr) ',
     1'  M/Mo  ',' log M ',' R/Ro')
      elseif (iso.eq.2) then

!      print*,' '
!      print*,'Enter the distance (pc)'
!      print*,' '
!      read(*,*) dist
      distmod=5.d0*dlog10(dist)-5.d0

!      print*,' '
!      print*,'Enter the reddening E(b-y)'
!      print*,' '
!      read(*,*) eby

!      print*,' '
!      print*,'Enter the the age (Myr)'
!      print*,' '
!      read(*,*) XA

      IF(XA.LT.1.d0) THEN
         XLA=0.d0
        ELSE
         XLA=DLOG10(XA*1.d6)
        ENDIF
      else
        close(1)
        goto 197
      endif

      If(imod.eq.3) then
        NZ=5
        il=1
328     continue
        if(il.eq.1) then
            open(10,status='old',file='data/TRACKS0004.ITA')
        else if(il.eq.2) then
            close(10)
            open(10,status='old',file='data/TRACKS0040.ITA')
        else if(il.eq.3) then
            close(10)
            open(10,status='old',file='data/TRACKS0080.ITA')
        else if(il.eq.4) then
            close(10)
            open(10,status='old',file='data/TRACKS0200.ITA')
        else if(il.eq.5) then
            close(10)
            open(10,status='old',file='data/TRACKS0500.ITA')
        else
            close(10)
            goto 122
        endif

        READ(10,*)
        READ(10,502) NMA,JMP,ZM(il)
502     FORMAT(I3,I4,f7.4)
        READ(10,*)
        DO 177 KI=1,NMA
         READ(10,503) XMASR
503      FORMAT(F7.3)
         DLM0(KI,il)=DLOG10(XMASR)
         READ(10,*)
         KM=1
178      READ(10,504) DDL0(KI,KM,il),DT0(KI,KM,il),DLO0(KI,KM,il),XM
     1  ,GG0(KI,KM,il)
504      FORMAT(F7.3,F6.3,F11.4,F8.3,F11.3)
         IF(XM.LT.0.01d0) GOTO 180
         XMAS0(KI,KM,il)=DLOG10(XM)
         KM=KM+1
         GOTO 178
180      DO 179 I=1,4
          READ(10,*,END=177)
179      CONTINUE
177     CONTINUE
                                                                       
        il=il+1
                                                                        
        GOTO 328
122     CONTINUE

      else
        NZ=9
        il=1
327     continue
        if(il.eq.1) then
          if(imod.eq.1) then
          open(10,status='old',file='data/TRACKS0001_o.CAS')
          else
          open(10,status='old',file='data/TRACKS0001_s.CAS')
          endif
        else if(il.eq.2) then
          close(10)
          if(imod.eq.1) then
          open(10,status='old',file='data/TRACKS0003_o.CAS')
          else
          open(10,status='old',file='data/TRACKS0003_s.CAS')
          endif
        else if(il.eq.3) then
          close(10)
          if(imod.eq.1) then
          open(10,status='old',file='data/TRACKS0010_o.CAS')
          else
          open(10,status='old',file='data/TRACKS0010_s.CAS')
          endif
        else if(il.eq.4) then
          close(10)
          if(imod.eq.1) then
          open(10,status='old',file='data/TRACKS0020_o.CAS')
          else
          open(10,status='old',file='data/TRACKS0020_s.CAS')
          endif
        else if(il.eq.5) then
          close(10)
          if(imod.eq.1) then
          open(10,status='old',file='data/TRACKS0040_o.CAS')
          else
          open(10,status='old',file='data/TRACKS0040_s.CAS')
          endif
        else if(il.eq.6) then
          close(10)
          if(imod.eq.1) then
          open(10,status='old',file='data/TRACKS0100_o.CAS')
          else
          open(10,status='old',file='data/TRACKS0100_s.CAS')
          endif
        else if(il.eq.7) then
          close(10)
          if(imod.eq.1) then
          open(10,status='old',file='data/TRACKS0200_o.CAS')
          else
          open(10,status='old',file='data/TRACKS0200_s.CAS')
          endif
        else if(il.eq.8) then
          close(10)
          if(imod.eq.1) then
          open(10,status='old',file='data/TRACKS0300_o.CAS')
          else
          open(10,status='old',file='data/TRACKS0300_s.CAS')
          endif
        else if(il.eq.9) then
          close(10)
          if(imod.eq.1) then
          open(10,status='old',file='data/TRACKS0400_o.CAS')
          else
          open(10,status='old',file='data/TRACKS0400_s.CAS')
          endif
        else
          close(10)
          goto 121
        endif

        READ(10,*)
        READ(10,402) NMA,JMP,ZM(il)
402     FORMAT(I3,I6,f7.4)
        READ(10,*)
        DO 77 KI=1,NMA
         READ(10,403) XMASR
403      FORMAT(F5.2)
         DLM0(KI,il)=DLOG10(XMASR)
         READ(10,*)
         KM=1
78       READ(10,404) DDL0(KI,KM,il),DT0(KI,KM,il),DLO0(KI,KM,il),XM
     1  ,GG0(KI,KM,il)
404      FORMAT(F7.3,F6.3,1X,F11.4,F7.3,f10.3)
         IF(XM.LT.0.01d0) GOTO 80
         XMAS0(KI,KM,il)=DLOG10(XM)
         KM=KM+1
         GOTO 78
80       DO 79 I=1,4
          READ(10,*,END=77)
79       CONTINUE
77      CONTINUE

        il=il+1

        GOTO 327
121     CONTINUE
      endif

      IF((Zest-ZM(NZ))*(Zest-ZM(1)).gt.0.) stop 'Fe/H' 

      write(1,72) feh,XA,dist,eby
 72   Format('# [Fe/H] = ',f7.3,';  age = ',f9.3,' Myr;  distance = ',
     1  f9.2,' pc;  E(b-y) = ',f7.3,' mag')
      write(1,74)
 74   Format('#   V   ',' (b-y) ','  log L ','  log T ','  log g ',
     1 '   M/Mo  ',' log M ',' R/Ro')
 
      DO 32 M=1,NMA
       DO 32 J=1,JMP
        DO 33 K=1,NZ
          If(Zest.gt.0.d0) then
           Zint(K)=DLOG10(ZM(K))
          else
            Zint(K)=ZM(K)
          endif
          DLMI(K)=DLM0(M,K)
          XMASI(K)=XMAS0(M,J,K)
          GGI(K)=GG0(M,J,K)
          DTI(K)=DT0(M,J,K)
          DLOI(K)=DLO0(M,J,K)
          DDLI(K)=DDL0(M,J,K)
 33     CONTINUE
        IF(J.EQ.1) THEN
         Call Polint (2,Zlog,DLM(M),Zint,DLMI,BU1,0,NZ,2)
        ENDIF
        Call Polint (2,Zlog,XMAS(M,J),Zint,XMASI,BU1,0,NZ,2)
        Call Polint (2,Zlog,GG(M,J),Zint,GGI,BU1,0,NZ,2)
        Call Polint (2,Zlog,DT(M,J),Zint,DTI,BU1,0,NZ,2)
        Call Polint (2,Zlog,DDL(M,J),Zint,DDLI,BU1,0,NZ,2)
        Call Polint (2,Zlog,DLO(M,J),Zint,DLOI,BU1,0,NZ,2)
 32   CONTINUE

      K=0
      XLMEI=DLM(NMA)
      CALL UVBY (0,feh,4.5d0,3.762d0,a1,a2,a3,a4)
 
499   IF (iso.eq.2) THEN
        XLME=DLOG10((10.D0**XLMEI)+0.001D0*K)
        IF((XLME-DLM(NMA))*(XLME-DLM(1)).gt.0.d0) GOTO 899
      ELSE
        IF((XLME-DLM(NMA))*(XLME-DLM(1)).gt.0.d0) THEN
          CLOSE(1)
          GOTO 198
        ENDIF
      ENDIF

      DO 898 J=1,JMP
        Call Polint (2,XLME,XLMF(J),DLM,BU2,XMAS,J,NMA,0)
        Call Polint (2,XLME,GF(J),DLM,BU2,GG,J,NMA,0)
        Call Polint (2,XLME,DTF(J),DLM,BU2,DT,J,NMA,0)
        Call Polint (2,XLME,DLAF(J),DLM,BU2,DLO,J,NMA,0)
        Call Polint (2,XLME,DLF(J),DLM,BU2,DDL,J,NMA,0)
898   CONTINUE

      If(iso.eq.1) THEN
        DO 897 J=1,JMP
          XMF(J)=10.D0**(XLMF(J))
          Rad=10.d0**(0.5d0*(4.43772d0-GF(J)+XLMF(J)))
          WRITE(1,976) DLF(J),DTF(J),GF(J),DLAF(J),
     1 XMF(J),XLMF(J),Rad
976       FORMAT(3f8.4,f11.4,f8.3,f7.3,f7.3)
897     CONTINUE
      else
        IF(XA.LT.DLAF(1)) THEN
          YMP=10.D0**XLMF(1)
          YLM=XLMF(1)
          YG=GF(1)
          YT=DTF(1)
          YL=DLF(1)
          GOTO 258
        ENDIF

        DO 289 J=2,JMP
          IF((XA-DLAF(J-1))*(XA-DLAF(J)).LE.0.D0) THEN
            JC=J
            GOTO 290
          ENDIF
289     CONTINUE
        GOTO 899
 
290     DGZM=GF(1)-GF(JC)
C       IF(DGZM.LT.0.15d0.AND.XLME.GT..301d0) THEN
C         aint=(XA-DLAF(JC-1))/(DLAF(JC)-DLAF(JC-1))
C       ELSE
          aint=(XLA-DLOG10(1.d6*DLAF(JC-1)))/
     1 (DLOG10(1.d6*DLAF(JC))-DLOG10(1.d6*DLAF(JC-1)))
C       ENDIF
        YLM=XLMF(JC-1)+aint*(XLMF(JC)-XLMF(JC-1))
        YMP=10.D0**YLM
        YG=GF(JC-1)+aint*(GF(JC)-GF(JC-1))
        YT=DTF(JC-1)+aint*(DTF(JC)-DTF(JC-1))
        YL=DLF(JC-1)+aint*(DLF(JC)-DLF(JC-1))
 
        Rad=10.d0**(0.5d0*(4.43772d0-YG+YLM))

        CALL UVBY(1,feh,yg,yt,by0,xm0,xc0,bcv)           

        by=by0+eby
        xmbol=xmbsun-2.5d0*yl
        vmag=xmbol-bcv+distmod+4.27d0*eby

258     WRITE(1,977) vmag,by,YL,YT,YG,YMP,YLM,Rad
977     FORMAT(f8.3,f7.3,3f8.4,f8.3,2f7.3)
 
        K=K+1
        GOTO 499

      endif

899   CLOSE(1)
!      STOP
      END subroutine pymodels
C********************************************************************
      Subroutine polint(n,x,y,xa,ya,yb,Kpu,j,ico)
        Implicit double precision(a-h,o-z)
        Dimension xa(35),ya(35),yb(35,1190),c(35),d(35)

        If(ico.eq.1) then
          noff=0
          goto 15
        elseif (ico.eq.0) then
          Do 99 kk=1,j
 99       ya(kk)=yb(kk,Kpu)
        endif

        Do 5 i=1,j
          If (ico.eq.2) then
            dif=x-xa(i)
          elseif (ico.eq.0) then
            dif=xa(i)-x
          endif
          If(dif.lt.0.d0) then
            noff=i-n/2-1
            If(noff.lt.0) noff=0
            goto 6
          endif
  5     Continue
        noff=j-n
  6     If(noff+n.gt.j) noff=j-n

C
C    Subrutina polint propiament dita
C

 15     ns=noff+1
        dif=dabs(x-xa(noff+1))
        do i=noff+1,noff+n
          dift=dabs(x-xa(i))
          if (dift.lt.dif) then
            ns=i
            dif=dift
          endif
          c(i)=ya(i)
          d(i)=ya(i)
        enddo
        y=ya(ns)
        ns=ns-1
        do m=1,n-1
          do i=noff+1,noff+n-m
            ho=xa(i)-x
            hp=xa(i+m)-x
            w=c(i+1)-d(i)
            den=ho-hp
            if (den.eq.0.d0) pause 'errada a polint'
            den=w/den
            d(i)=hp*den
            c(i)=ho*den
          enddo
          if (2*(ns-noff).lt.n-m) then
            dy=c(ns+1)
          else
            dy=d(ns)
            ns=ns-1
          endif
          y=y+dy
        enddo
        return
        END
********************************************************************
      SUBROUTINE UVBY(MODE,FE,GV,TEFF,SBY,SM1,SC1,BCV)           
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
C ----------------------------------------------------------------------
C *** THIS READS, AND INTERPOLATES IN, UVBYLO.DATA AND UVBYHI.DATA (SEE
C     CLEM, VANDENBERG, GRUNDAHL, & BELL, AJ, IN PRESS).            
C *** SET MODE=0 TO READ THE TABLES AND TO INTERPOLATE THEM TO THE
C                DESIRED VALUE OF [Fe/H]  (= FE IN THE ARGUMENT LIST)   
C *** SET MODE=1 TO INTERPOLATE FOR THE COLORS AT THE DESIRED VALUES OF
C                [Fe/H], log g, and log Teff (GV AND TEFF REPRESENT THE
C                LAST TWO OF THESE QUANTITIES IN THE ARGUMENT LIST).
C                (NOTE THAT THE [Fe/H] INTERPOLATION IS CARRIED OUT ONLY
C                WHEN THE SUBROUTINE IS CALLED WITH MODE=0 OR MODE=-1.)
C *** SET MODE=-1 IF THE INPUT TABLES, WHICH HAVE ALREADY BEEN READ
C                USING MODE=0, ARE TO BE RE-INTERPOLATED TO BE
C                CONSISTENT WITH A NEW VALUE OF [Fe/H].        
C *** THE OUTPUT CONSISTS OF THE b-y, m1, AND c1 INDICES, ALONG 
C     WITH THE BOLOMETRIC CORRECTION TO V (ON THE SCALE WHERE THE SUN 
C     HAS M_bol = 4.75 and M_V = 4.84).  THESE QUANTITIES ARE 
C     REPRESENTED, IN TURN, BY SBY, SM1, SC1, AND BCV IN THE 
C     ARGUMENT LIST.   NOTE THAT IT IS ASSUMED THE BOLOMETRIC 
C     CORRECTIONS TO STROMGREN y ARE IDENTICAL TO THOSE IN JOHNSON V.
C *** MODE MUST BE SET TO ZERO THE FIRST TIME THAT THIS SUBROUTINE IS
C     CALLED, AS IT IS ONLY WHEN MODE=0 THAT UVBYLO.DATA AND UVBYHI.DATA
C     ARE READ.  ONCE THE TABLES HAVE BEEN INTERPOLATED TO A PARTICULAR
C     VALUE OF [Fe/H] (USING MODE=0 OR MODE=-1), INTERPOLATIONS FOR THE
C     COLORS APPROPRIATE TO INPUT VALUES OF [Fe/H], log g, and log Teff
C     MAY BE CARRIED OUT *ANY NUMBER OF TIMES* USING MODE=1 (THE ONLY
C     VALUE OF MODE FOR WHICH COLOR INFORMATION IS OBTAINED).
C *** NOTE THAT THE UVBYLO.DATA AND UVBYHI.DATA FILES ARE "ATTACHED" TO
C     UNITS 12 AND 13, RESPECTIVELY.  WHENEVER AN [Fe/H] INTERPOLATION 
C     IS CARRIED OUT, A MESSAGE IS SENT TO UNIT 6 ("ATTACHED" TO THE  
C     TERMINAL) TO INDICATE THE VALUE OF [FE/H] THAT APPLIES TO
C     SUBSEQUENT INTERPOLATIONS FOR log g and log Teff. 
C *** CHKUVBY.FOR (ALSO ON DISK) MAY BE USED TO TEST THE INTERPOLATION
C     SUBROUTINE
C ----------------------------------------------------------------------
      REAL*8 A(51,7,8,4),B(51,7,4),BY(13,12,8),M1(13,12,8),C1(13,12,8),
     1  BC(13,12,8),C(13,12),D(13,12),E(13,12),F(13,12),H(4,4),TB(51),
     2  TT(51),GG(7),TA(13),T(13),G(12),FEH(8),O(4),P(4),Q(4),R(4),S(4),
     3  AG(4),AT(4),X(13),Y(13)
      CHARACTER*19 NAMCOL
      CHARACTER*10 NAMFE                                       
 1000 FORMAT(I3,13X,I3,13X,I3,13X,I3)                                   
 1001 FORMAT(13F6.0)                                                    
 1002 FORMAT(20F4.1)                                                    
 1003 FORMAT(13F6.3)                                                    
 1004 FORMAT('     Color grid interpolated to [FE/H] =',F6.2)           
 1005 FORMAT(1X,7HLOG G =,F6.3,11H LOG TEFF =,F7.4,8H OUTSIDE,          
     1   1X,11HCOLOR TABLE)                                             
 1006 FORMAT(A10,F5.2,A19)
 1007 FORMAT(1X,' **** INPUT DATA FILE DOES NOT EXIST **** ')      
      SAVE
      IF(MODE) 11,3,16                                                  
C *** WHEN MODE=0 THE COLOR TRANSFORMATION TABLES ARE READ
    3 TBND=DLOG10(5500.d0)                                                
      GBND=2.
      OPEN(UNIT=12,ERR=49,FILE='data/uvbylo.data',STATUS='OLD')
      OPEN(UNIT=13,ERR=49,FILE='data/uvbyhi.data',STATUS='OLD')      
      READ(12,1000) NT,NG,NFE,NDX                                       
      READ(12,1001) (TA(I),I=1,NT)                                      
      READ(12,1002) (G(I),I=1,NG)                                       
      DO 4 I=1,NT                                                       
    4 T(I)=DLOG10(TA(I))                                                
      DO 5 K=1,NFE                                                      
      READ(12,1006) NAMFE,FEH(K),NAMCOL                                    
      DO 5 J=1,NG                                                       
      READ(12,1003) (BY(I,J,K),I=1,NT)                                  
    5 CONTINUE                                                          
      DO 6 K=1,NFE                                                      
      READ(12,1006) NAMFE,FEH(K),NAMCOL                                 
      DO 6 J=1,NG                                                       
      READ(12,1003) (M1(I,J,K),I=1,NT)                                  
    6 CONTINUE                                                          
      DO 7 K=1,NFE                                                      
      READ(12,1006) NAMFE,FEH(K),NAMCOL                                   
      DO 7 J=1,NG                                                       
      READ(12,1003) (C1(I,J,K),I=1,NT)                                  
    7 CONTINUE                                                          
      DO 8 K=1,NFE                                                      
      READ(12,1006) NAMFE,FEH(K),NAMCOL                                   
      DO 8 J=1,NG
      READ(12,1003) (BC(I,J,K),I=1,NT)
    8 CONTINUE                                             
      READ(13,1000) NTT,NGG,NFE,NDXX                                    
      READ(13,1001) (TB(I),I=1,NTT)                                     
      READ(13,1002) (GG(I),I=1,NGG)                                     
      DO 9 I=1,NTT                                                      
    9 TT(I)=DLOG10(TB(I))                                               
      DO 10 L=1,NDXX                                                    
      DO 10 K=1,NFE                                                     
      READ(13,1006) NAMFE,FEH(K),NAMCOL
      DO 10 J=1,NGG                                                     
      READ(13,1003) (A(I,J,K,L),I=1,NTT)
   10 CONTINUE
      CLOSE(UNIT=12,STATUS='KEEP')
      CLOSE(UNIT=13,STATUS='KEEP')         
C *** WHEN MODE=0 OR MODE=-1, COLOR TRANSFORMATION TABLES ARE
C *** CREATED FOR THE INPUT [Fe/H] VALUE USING LINEAR INTERPOLATION
   11 DO 12 M=2,NFE                                                     
      K=M-1                                                             
      IF(FE.LE.FEH(M)) GO TO 13                                         
   12 CONTINUE                                                          
   13 M=K+1                                                             
      SLOPE=(FE-FEH(K))/(FEH(M)-FEH(K))                                 
      DO 14 J=1,NG                                                      
      DO 14 I=1,NT                                                      
      C(I,J)=BY(I,J,K)+SLOPE*(BY(I,J,M)-BY(I,J,K))
      D(I,J)=M1(I,J,K)+SLOPE*(M1(I,J,M)-M1(I,J,K))
      E(I,J)=C1(I,J,K)+SLOPE*(C1(I,J,M)-C1(I,J,K))
      F(I,J)=BC(I,J,K)+SLOPE*(BC(I,J,M)-BC(I,J,K))
   14 CONTINUE                                                          
      DO 15 L=1,NDXX                                                    
      DO 15 J=1,NGG                                                     
      DO 15 I=1,NTT                                                     
      B(I,J,L)=A(I,J,K,L)+SLOPE*(A(I,J,M,L)-A(I,J,K,L))                 
   15 CONTINUE                                                          
!      WRITE(*,1004) FE  
C *** WHEN MODE=0 OR MODE=-1, CONTROL RETURNS TO THE CALLING
C *** PROGRAM ONCE THE [Fe/H] INTERPOLATION IS CARRIED OUT (I.E., NO
C *** INTERPOLATIONS ARE PERFORMED FOR INPUT VALUES OF GV AND TEFF)
      GO TO 50
C *** WHEN MODE=1, INTERPOLATIONS ARE CARRIED OUT FOR THE INPUT VALUES
C *** OF [Fe/H], log g, AND log Teff (= FE, GV, AND TEFF IN THE ARGUMENT
C *** LIST)
   16 IF(TEFF.LE.TBND) GO TO 18                                         
      IF(GV.GE.GBND) GO TO 40                                           
      IF(TEFF.LE.T(NT)) GO TO 18                                        
C *** EXECUTION HALTS WITH A "STOP30" CODE IF THE INPUT TEMPERATURE IS
C *** OUTSIDE THE RANGE OF THE TABLES ON THE HIGH SIDE
      WRITE(*,1005) GV,TEFF                                             
      STOP30
C *** THE NEXT SECTION ASSUMES THAT THE LOW-TEMPERATURE TABLES ARE THE
C *** RELEVANT ONES TO USE
   18 NGM=NG-1                                                          
      DO 19 I=3,NGM                                                     
      MG=I-2                                                            
      IF(GV.LE.G(I)) GO TO 20                                           
   19 CONTINUE                                                          
      IF(GV.LE.G(NG)) GO TO 20                                          
C *** SOME EXTRAPOLATION OUTSIDE THE RANGE OF log g CONSIDERED IN THE 
C *** TABLES IS PERMITTED.  TO AVOID ANY EXTRAPOLATIONS, THE NEXT TWO
C *** LINES SHOULD BE ACTIVATED   
C     WRITE(*,1005) GV,TEFF
C     STOP31
   20 R(1)=G(MG)                                                        
      R(2)=G(MG+1)                                                      
      R(3)=G(MG+2)                                                      
      R(4)=G(MG+3)                                                      
      CALL LGRAN4(R,AG,GV)                                              
      NTM=NT-1                                                          
      DO 25 I=3,NTM                                                     
      MT=I-2                                                            
      IF(TEFF.LE.T(I)) GO TO 30                                         
   25 CONTINUE                                                          
      IF(TEFF.LE.T(NT)) GO TO 30                                        
C *** SOME EXTRAPOLATION OUTSIDE THE RANGE OF log Teff CONSIDERED IN
C *** THE TABLES IS PERMITTED.  TO AVOID ANY EXTRAPOLATIONS, THE NEXT
C *** TWO LINES SHOULD BE ACTIVATED.
C     WRITE(*,1005) GV,TEFF
C     STOP32
   30 R(1)=T(MT)                                                        
      R(2)=T(MT+1)                                                      
      R(3)=T(MT+2)                                                      
      R(4)=T(MT+3)                                                      
      CALL LGRAN4(R,AT,TEFF)                                            
      L=MG-1                                                            
      DO 35 I=1,4                                                       
      L=L+1                                                             
      O(I)=AT(1)*C(MT,L)+AT(2)*C(MT+1,L)+AT(3)*C(MT+2,L)+AT(4)*C(MT+3,L)
      P(I)=AT(1)*D(MT,L)+AT(2)*D(MT+1,L)+AT(3)*D(MT+2,L)+AT(4)*D(MT+3,L)
      Q(I)=AT(1)*E(MT,L)+AT(2)*E(MT+1,L)+AT(3)*E(MT+2,L)+AT(4)*E(MT+3,L)
   35 R(I)=AT(1)*F(MT,L)+AT(2)*F(MT+1,L)+AT(3)*F(MT+2,L)+AT(4)*F(MT+3,L)
C *** 4-POINT LAGRANGIAN INTERPOLATION IS USED TO DERIVE THE VALUES OF
C *** b-y, m1, c1, AND BC_V CORRESPONDING TO THE INPUT VALUES OF
C *** [Fe/H], log g, and log Teff.
      SBY=AG(1)*O(1)+AG(2)*O(2)+AG(3)*O(3)+AG(4)*O(4)
      SM1=AG(1)*P(1)+AG(2)*P(2)+AG(3)*P(3)+AG(4)*P(4)                   
      SC1=AG(1)*Q(1)+AG(2)*Q(2)+AG(3)*Q(3)+AG(4)*Q(4)                   
      BCV=AG(1)*R(1)+AG(2)*R(2)+AG(3)*R(3)+AG(4)*R(4)                  
C *** SPLINE INTERPOLATION IS USED TO FIND THE VALUE OF BC_V AT LOW
C *** TEMPERATURES.
      DO 38 I=1,NT                                                      
      X(I)=T(I)                                                         
   38 Y(I)=AG(1)*F(I,MG)+AG(2)*F(I,MG+1)+AG(3)*F(I,MG+2)+AG(4)*F(I,MG+3)
      CALL INTEP(TEFF,BCV,X,Y,NT,IER)                                  
      GO TO 50                                                          
C *** THE NEXT SECTION ASSUMES THAT THE HIGH-TEMPERATURE TABLES ARE THE
C *** RELEVANT ONES TO USE.
   40 IF(TEFF.LE.TT(NTT)) GO TO 42                                      
      WRITE(*,1005) GV,TEFF                                             
      STOP33                                                            
   42 NGM=NGG-1                                                         
      DO 43 I=3,NGM                                                     
      MG=I-2                                                            
      IF(GV.LE.GG(I)) GO TO 44                                          
   43 CONTINUE                                                          
      IF(GV.LE.GG(NGG)) GO TO 44                                        
C *** SOME EXTRAPOLATION OUTSIDE THE RANGE OF log g CONSIDERED IN THE
C *** TABLES IS PERMITTED.  TO AVOID ANY EXTRAPOLATIONS, THE NEXT TWO
C *** LINES SHOULD BE ACTIVATED.
C     WRITE(*,1005) GV,TEFF   
C     STOP34
   44 R(1)=GG(MG)                                                       
      R(2)=GG(MG+1)                                                     
      R(3)=GG(MG+2)                                                     
      R(4)=GG(MG+3)                                                     
      CALL LGRAN4(R,AG,GV)                                              
      NTM=NTT-1                                                         
      DO 45 I=3,NTM                                                     
      MT=I-2                                                            
      IF(TEFF.LE.TT(I)) GO TO 46                                        
   45 CONTINUE                                                          
      IF(TEFF.LE.TT(NTT)) GO TO 46                                      
C *** SOME EXTRAPOLATION OUTSIDE THE RANGE OF log Teff CONSIDERED IN
C *** THE TABLES IS PERMITTED.  TO AVOID ANY EXTRAPOLATIONS, THE NEXT
C *** TWO LINES SHOULD SHOULD BE ACTIVATED.
C     WRITE(*,1005) GV,TEFF
C     STOP35
   46 R(1)=TT(MT)                                                       
      R(2)=TT(MT+1)                                                     
      R(3)=TT(MT+2)                                                     
      R(4)=TT(MT+3)                                                     
      CALL LGRAN4(R,AT,TEFF)                                            
      L=MG-1                                                            
      DO 48 I=1,4                                                       
      L=L+1                                                             
      DO 48 K=1,NDXX                                                    
      H(I,K)=AT(1)*B(MT,L,K)+AT(2)*B(MT+1,L,K)+AT(3)*B(MT+2,L,K)+       
     1   AT(4)*B(MT+3,L,K)                                              
   48 CONTINUE                                                          
C *** 4-POINT LAGRANGIAN INTERPOLATION IS USED TO DERIVE THE VALUES OF
C *** b-y, m1, c1, AND BC_V CORRESPONDING TO THE INPUT VALUES OF
C *** [Fe/H], log g, and log Teff.
      SBY=AG(1)*H(1,1)+AG(2)*H(2,1)+AG(3)*H(3,1)+AG(4)*H(4,1)           
      SM1=AG(1)*H(1,2)+AG(2)*H(2,2)+AG(3)*H(3,2)+AG(4)*H(4,2)           
      SC1=AG(1)*H(1,3)+AG(2)*H(2,3)+AG(3)*H(3,3)+AG(4)*H(4,3)           
      BCV=AG(1)*H(1,4)+AG(2)*H(2,4)+AG(3)*H(3,4)+AG(4)*H(4,4)
      GO TO 50
   49 WRITE(*,1007)
      STOP          
   50 RETURN                                                            
      END                                                               
C ----------------------------------------------------------------------
      SUBROUTINE LGRAN4(X,A,XX)                                         
C *** A 4-POINT LAGRANGIAN INTERPOLATION SUBROUTINE
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      SAVE
      REAL*8 X(4),A(4)                                                  
      R1=(X(1)-X(2))*(X(1)-X(3))*(X(1)-X(4))                            
      R2=(X(2)-X(1))*(X(2)-X(3))*(X(2)-X(4))                            
      R3=(X(3)-X(1))*(X(3)-X(2))*(X(3)-X(4))                            
      R4=(X(4)-X(1))*(X(4)-X(2))*(X(4)-X(3))                            
      A(1)=((XX-X(2))*(XX-X(3))*(XX-X(4)))/R1                           
      A(2)=((XX-X(1))*(XX-X(3))*(XX-X(4)))/R2                           
      A(3)=((XX-X(1))*(XX-X(2))*(XX-X(4)))/R3                           
      A(4)=((XX-X(1))*(XX-X(2))*(XX-X(3)))/R4                           
      RETURN                                                            
      END                                                               
C ----------------------------------------------------------------------
      SUBROUTINE INTEP(XP,P,X,F,N,IER)                                  
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
C *** PURPOSE:  To interpolate a function value P for a given argument  
C     XP from a table of N values (X,F).  This is a spline interpolation
C     scheme based on Hermite polynomials.  The source is U.S. Airforce 
C     Surveys in Geophysics No 272.                                     
C *** USAGE:  For random values of XP                                   
C               CALL INTEP(XP,P,X,F,N,IER)                              
C     or after the first call to INTEP with monotonically increasing or 
C     decreasing values of XP consistent with the X vector              
C               CALL EINTEP(XP,P,X,F,N,IER)                             
C     DESCRIPTION OF PARAMETERS:                                        
C     XP  - the chosen argument value                                   
C     P   - the resultant interpolated value                            
C     X   - the vector of independent values                            
C     F   - the vector of dependent values                              
C     N   - the number of points in the (X,F) vectors                   
C     IER - the resultant error parameter (set = 2 if XP is beyond      
C           either extreme of X; in which case P is set to the value of 
C           F at that extremum)                                         
      REAL*8 LP1,LP2,L1,L2                                              
      DIMENSION F(13),X(13)                                               
      IER=1                                                             
      IO=1                                                              
      IUP=0                                                             
      IF(X(2).LT.X(1)) IUP=1                                            
      N1=N-1                                                            
    5 IF((XP.GE.X(N).AND.IUP.EQ.0).OR.(XP.LE.X(N).AND.IUP.EQ.1)) THEN   
       P=F(N)                                                           
       GO TO 6                                                          
      ELSE IF((XP.LE.X(1).AND.IUP.EQ.0).OR.                             
     1   (XP.GE.X(1).AND.IUP.EQ.1)) THEN                                
       P=F(1)                                                           
    6  IER=2                                                            
       RETURN                                                           
      ENDIF                                                             
      ENTRY EINTEP(XP,P,X,F,N,IER)                                      
    8 DO 1 I=IO,N                                                       
      IF(XP.LT.X(I).AND.IUP.EQ.0) GO TO 2                               
      IF(XP.GT.X(I).AND.IUP.EQ.1) GO TO 2                               
    1 CONTINUE                                                          
      GO TO 5                                                           
    2 I=I-1                                                             
      IF(I.EQ.IO-1) GO TO 4                                             
      IO=I+1                                                            
      LP1=1./(X(I)-X(I+1))                                              
      LP2=1./(X(I+1)-X(I))                                              
      IF(I.EQ.1) FP1=(F(2)-F(1))/(X(2)-X(1))                            
      IF(I.EQ.1) GO TO 3                                                
      FP1=(F(I+1)-F(I-1))/(X(I+1)-X(I-1))                               
    3 IF(I.GE.N1) FP2=(F(N)-F(N-1))/(X(N)-X(N-1))                       
      IF(I.EQ.N1) GO TO 4                                               
      FP2=(F(I+2)-F(I))/(X(I+2)-X(I))                                   
    4 XPI1=XP-X(I+1)                                                    
      XPI=XP-X(I)                                                       
      L1=XPI1*LP1                                                       
      L2=XPI*LP2                                                        
      P=F(I)*(1.-2.*LP1*XPI)*L1*L1+F(I+1)*(1.-2.*LP2*XPI1)*L2*L2+       
     1   FP2*XPI1*L2*L2+FP1*XPI*L1*L1                                   
      RETURN                                                            
      END                                                               
