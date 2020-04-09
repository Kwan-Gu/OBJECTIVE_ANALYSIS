$DEBUG
C*********************************************************************
C    THIS PROGRAM IS TO FILTER BY CARTWRIGHT'S FILTER AND
C          CALCULATE THE MEAN AND VARIANCE
C                           MODIFIED BY Lyu Sang Jin, 1998. 5. 5
C*********************************************************************
C
C          PARAMETERS;IFC,MM,IB1,IB2
C
C      IFC=0; CARTWRIGHT'S REAL FILTER
C          1; CARTWRIGHT'S REAL AND CONJUGATE FILTER
C
      PARAMETER(IX1=100000,M=10000)
      REAL*8 XM,XMM,VAR
      COMMON/A/NW,BK(20001,2),t
      COMMON/B/DAY(IX1)
      COMMON/C/X(IX1)
      CHARACTER*12 infile,outfile
      CHARACTER*50 inFMT,outfmt
      real   ib1,ib2
C
       NW=5
       PHI=3.1415926536
C
C     N ; NO. OF DATA IN EACH SET
C     MM; HALF NO. OF FILTER WEIGHT TO BE GENERATED
C         AND SHOULD BE EVEN NUMBER.
C     NDT; INTERVAL OF THE FILTERED DATA
C     T ; SAMPLE TIME INTERVAL
C     IB1; LOW CUTOFF FREQ., IB2; HIGH CUTOFF FREQ.
C
      OPEN(NW,FILE='FILT.LIS')
C
      OPEN(111,FILE='filter.con',status='old')
      READ(111,*)
      READ(111,'(40x,A12)') infile
      READ(111,'(40x,A12)') outfile
      READ(111,'(40x,i5)') N
      READ(111,'(40x,F10.4)') T
      READ(111,'(40x,i8)') MM
      READ(111,'(40x,f10.4)') IB1
      READ(111,'(40x,f10.4)') IB2
      READ(111,'(40x,i8)') NDT
      READ(111,'(40x,i8)') NCUT
      READ(111,'(40x,i8)') IFC
      READ(111,'(40x,a50)') inFMT
      READ(111,'(40x,a50)') outFMT
      
      NU=11
      M1=MM+1
C
C************DATA READ IN******************************************
C
      OPEN(NU,FILE=INFILE,STATUS='OLD')
C
      do j = 1,n 
c         write(*,*) j
         READ(NU,infmt) DAY(j),X(J)
      enddo
C
C*****MEAN AND VARIANCE OF THE FILTERED DATA
C
      VAR=0.0
      XM=0.0
      DO 70 II=1,N
   70 XM=XM+X(II)
      XMM=XM/N
      DO 80 II=1,N
   80 X(II)=X(II)-XMM
      DO 90 II=1,N
   90 VAR=VAR+X(II)*X(II)
      XVAR=VAR/N
C
      WRITE(NW,210) INFILE
      WRITE(NW,220) N,NCUT,MM,IFC,NDT
      WRITE(NW,230) IB1,IB2,T

cccccccccccccccccccccccccccccc
c     detrend
cccccccccccccccccccccccccccccc

c      call detrnd(n)

C
  210 FORMAT(2X,'FILTERING TO THE STATION;  ',A12)
  220 FORMAT(3X,'RAW DATA NO.; ',I10,2X,'NO.CUT',I5,' NO.WEIGHT; ',I10,/
     *,7X,'KINDS OF FILTER ; ',I10,3X,'INTERVAL OF FILTERED DATA; ',I10)
  230 FORMAT(7X,'LOW Fc ; ',f10.4,3X,'HIGH Fc ; ',f10.4,/,10X
     *,'SAMLPE DT ; ',F10.5)
C
C******CARTWRIGHT'S FILTER GENERATION*********************
C
C     LOWPASS FILTER WEIGHT GENERATION
C       IB1=0, Fc;(IB2+0.5)/MM (in Nyquist frequency)
C
C     HIGHPASS
C       IB2=MM, Fc;(IB1-0.5)/MM (in Nyquist frequency)
C
C     BANDPASS  (in Nyquist frequency)
C        0<IB1<IB2<MM, Fc;(IB1-0.5)/MM---(IB2+0.5)/MM
C
      CALL FCART(MM,IB1,IB2,IFC)
C
C******ARRANGE THE FILTER WEIGHT
C
      IBKC=IFC+1
      M2=2*MM+1
      DO 819 KL=1,IBKC
      TEM=BK(M1,KL)
      DO 51 J=1,M1
   51 BK(J+MM,KL)=BK(J,KL)
      BK(M2,KL)=TEM
      DO 52 J=1,MM
      J1=J-1
   52 BK(J,KL)=BK(M2-J1,KL)
  819 CONTINUE
C
C******FILTERING PROCEDURE
C
      NP=0
      M3=M1+NCUT
      DO 60 J=M3,N-MM,NDT
      XF=0.0
      K=J-M1
      NP=NP+1
      DO 61 II=1,M2
  61  XF=XF+X(K+II)*BK(II,IFC+1)
      X(NP)=XF
   60 CONTINUE
      N=NP
C
C************END OF CARTWRIGHT'S FILTERING******************
C
       WRITE(NW,240) N,NCUT
  240 FORMAT(5X,'DATA NO. AFTER FILTERING; ',I10,' NO. CUT; ',I5)
      WRITE(NW,301) XMM,XVAR
 301  FORMAT(5X,'MEAN = ',F15.8,5X,'VARIANCE = ',F15.8)
C
C*****WRITING THE FILTERED DATA**************************************
C
      NOP=21
      OPEN(NOP,FILE=OUTFILE)
      DO 10 J=1,N
        JJ=(J-1)*NDT+MM+1
        XX=X(J)+XMM 
        WRITE(NOP,outFMT) DAY(JJ),XX
   10 CONTINUE
C
           STOP
           END

ccccccccccccccccccccccccccccccccccccc
      SUBROUTINE DETRND(N)
ccccccccccccccccccccccccccccccccccccc
      PARAMETER(IX1=100000)
      COMMON/c/X(IX1)
C
C     DETRENDING BY THE L.S.M
C
      TBAR=FLOAT(N+1)/2.
      ZN=FLOAT(N)
      SMTT=ZN*(ZN*ZN-1.0)/12.0
      SUMX=0.0
      DO 40 I=1,N
   40 SUMX=SUMX+X(I)*(FLOAT(I)-TBAR)
      BETA=Sumx/SMTT
      DO 50 I=1,N
   50 X(I)=X(I)-BETA*(FLOAT(I)-TBAR)
      RETURN
      END


C*********************************************************************
      SUBROUTINE FCART(NNN,NLO,NHI,INDX)
C*********************************************************************
      PARAMETER(IDM=5)
      REAL*8 C(20000),S(20000),F(nnn),G(nnn),H(6)
      COMMON/A/NW,BK(20001,2),t
      real nlo,nhi  
C
C      GENERATE HIGH-, LOW-, OR BAND-PASS FILTER.
C      WITH CONJUGATE IF REQUIRED. AMPLITUDE RESPONSE LISTRED
C        SOURCE; D.E. CARTWRIGHT
C
      DO 1 N=1,2500
      X=N/10000.
      C(N)=COS(3.1415926536*X)
      S(N)=SIN(3.1415926536*X)
      M=5000-N
      C(M)=S(N)
    1 S(M)=C(N)
      C(5000)=0.
      S(5000)=1.0
      DO 2 N=5001,10000
      C(N)=-S(N-5000)
    2 S(N)=C(N-5000)
   3  CONTINUE
      IF(NNN) 99,99,103
  103 XNN=0.5/NNN
      XLO=(2.*NLO-1.)*XNN
      XHI=(2.*NHI+1.)*XNN
      IF(NLO)4,4,5
    4 WRITE(NW,501) NNN,NHI,XHI
      IGO=-1
      XLO=0.
      GO TO 8
   5  IF(NHI-NNN) 6,7,7
   6  WRITE(NW,502) NNN,NLO,NHI,XLO,XHI
      IGO=0
      GO TO 8
    7 WRITE(NW,503)NNN,NLO,XLO
      IGO=1
      XHI=1.
    8 WRITE(NW,504)
C
  500 FORMAT(4I4)
  501 FORMAT(//31H  LOWPASS FILTER OF HALF-LENGTH,I4/
     *15H CUT-OFF NUMBER,f10.4,16H  I.E. FREQUENCY,F7.4,9H NYQUISTS)
  502 FORMAT(//31H BANDPASS FILTER OF HALF-LENGTH,I4/
     *16H CUT-OFF NUMBERS,2f10.4,18H  I.E. FREQUENCIES,2F10.6,9H NYQUIST
     *S)
  503 FORMAT(//31H HIGHPASS FILETR OF HALF-LENGTH,I4/
     *15H CUT-OFF NUMBER,f10.4,16H  I.E. FREQUENCY,F7.4,9H NYQUISTS)
  504 FORMAT(/45H MULTIPLIERS ROUNDED AT 8 D.P., ALSO ON FILES//5H REAL)
C
C
      F0=XHI-XLO
      G0=0.0
      F(NNN)=0.0
      G(NNN)=0.0
      SUMF=F0
      NN=NNN-1
      I=1
      DO 22 NT=1,NN
      X=3.1415926536*NT
      Z=XNN*X
      C0=COS(Z)
      S0=SIN(Z)
      ZZ=XNN*C0*C0/S0
      IF(IGO) 12,12,13
   12 Z=XHI*X
      C2=COS(Z)
      S2=SIN(Z)
   13 IF(IGO) 15,14,14
   14 Z=XLO*X
      C1=COS(Z)
      S1=SIN(Z)
   15 IF(IGO) 16,17,18
   16 S1=0.
      G(NT)=(C0-C2)*ZZ
      GO TO 21
   17 G(NT)=(C1-C2)*ZZ
      GO TO 21
   18 S2=0.
      I=-I
      IF(I) 19,19,20
   19 G(NT)=(C1+C0)*ZZ
      GO TO 21
   20 G(NT)=(C1-C0)*ZZ
   21 F(NT)=(S2-S1)*ZZ
   22 SUMF=SUMF+2.*F(NT)
C
C
      IBK=1
      I=0
      J=-1
      K=1
      NT=0
      IGO=1
      H(K)=F0
   23 IF(H(K)) 24,125,25
   24 H(K)=H(K)-0.5E-8
      GO TO 125
   25 H(K)=H(K)+0.5E-8
  125 GO TO (28,27,27,27),IGO
  126 IGO=2
   26 H(K)=F(NT)
      GO TO 23
   27 IF(K-6) 29,28,28
   28 J=J+1
      WRITE(NW,505) I,J,(H(N),N=1,K)
C
C
      IF(J.GE.1) GO TO 887
      BK(1,IBK)=H(1)
      GO TO 889
  887 NOP=(J-1)*6+1
      DO 888 N=1,K
  888 BK(NOP+N,IBK)=H(N)
  889 CONTINUE
  505 FORMAT(2I2,6F12.8)
C
C
      K=0
   29 NT=NT+1
      IF(NT-NNN) 30,30,31
   30 K=K+1
      GO TO (126,26,136,36),IGO
   31 IF(K) 33,33,28
   33 GO TO (126,34,136,37),IGO
   34 IF(INDX) 37,37,35
   35 WRITE(NW,506)
  506 FORMAT(//10H CONJUGATE)
C
C
      IBK=2
      I=1
      J=-1
      K=1
      NT=0
      IGO=3
      H(K)=0.
      GO TO 28
  136 IGO=4
   36 H(K)=G(K)
      GO TO 23
   37 WRITE(NW,507)
  507 FORMAT(//49H AMPLITUDE RESPONSE AT INTERVAL OF 0.01 NYQUISTS//
     *11X,4HREAL,6X,9HCONJUGATE)
  508 FORMAT(F10.4,2F12.6)
      X=0.
      IF(INDX) 38,38,39
   38 WRITE(NW,508) X,SUMF
      GO TO 40
   39 SUMG=0.0
      WRITE(NW,508) X,SUMF,SUMG
   40 DO 47 N=1,10000
      X=N
      X=0.0001*X+0.00001
      SUMF=F0
      SUMG=0.
      Z=2.
      K=0
      DO 44 NT=1,NN
      K=K+N
      IF(K-10000) 42,42,41
   41 K=K-10000
      Z=-Z
   42 IF(INDX) 44,44,43
   43 SUMG=SUMG+Z*G(NT)*S(k)
   44 SUMF=SUMF+z*F(NT)*C(k)
      IF(INDX) 45,45,46
   45 WRITE(NW,508) X,SUMF
      GO TO 47
   46 WRITE(NW,508) X,SUMF,SUMG
   47 CONTINUE
C        GO TO 3
   99 CONTINUE
        RETURN
        END




