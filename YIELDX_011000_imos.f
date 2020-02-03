
!!**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
!! * YIELDX_011000_imos.f *                        galprop package * 4/14/2000 
!!**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

c### the code was changed to make it compatible with NAG f95 compiler (Linux).
c### Changes are marked by "imos"; I.Moskalenko 2/15/2000
c###
      SUBROUTINE YIELDX(IZ, IA, JZ, JA, EJ, QJ)                 ! Dec 1999
C*****************************************************************************                                          
C       IZ IA :  charge and mass of the nucleus
C       JZ JA :  charge and mass of spallation product (He-6 and heavier)
C       EJ    :  energy in MeV/n
C       QJ    :  cross section in mb
C
C       The code was developed at NRL by Rein Silberberg and C. H. Tsao.
C       It is designed for quick estimation of cross sections only.
C       The results are not suited for physical interpretation.
C
C       This latest version is undertaken with NASA's sponsorship
C       and the assistance of Dr. A.F. Barghouty at the Roanoke College
C       
C         E-mail:  tsao@crs2.nrl.navy.mil
C
C       The code is available at:  http//spdsch.phys.lsu.edu
C*****************************************************************************
      COMMON /FS/ QR,QE,QF,QH,FE,FF,FA,FZ,CJ,PN,GA,ANZJ,AA,AE,AC,EC
      QJ = 0
      IF (IZ*IA .EQ. JZ*JA)  RETURN
      A = IA
      Z = IZ
      IJ = IZ - JZ
      NN = (IA - IZ) - (JA - JZ)
      X = NN
      IF (IA.LT.JA .OR.  IZ.GT.92 .OR.  IZ.LE. 2)  RETURN
      IF (IJ.LT.-1 .OR.  NN.LT. 0 .OR.  JZ.LE. 1)  RETURN
      IF (IJ.EQ.-1 .AND. IA.GE.JA)  THEN
        CALL PXN(Z, A, X, EJ, QJ)
        RETURN
      ENDIF
      IF (IZ-JZ.GT.IA-JA .OR. IZ*IA.EQ.JZ*JA)  RETURN
      IF (JZ.LE. 4)   CALL YIELD1(IZ, IA, JZ, JA, EJ, QJ)
      IF (IZ.EQ. 6 .AND. EJ.GT.200. .AND. EJ.LE.400.)
     1QJ = QJ*(1.-0.002*(EJ-200. ))
      IF (IZ.EQ. 6 .AND. EJ.GT.400. .AND. EJ.LE.1000.)
     1QJ = QJ*0.6
      IF (IZ.EQ. 6 .AND. EJ.GT.1000..AND. EJ.LE.5000.)
     1QJ = QJ*(.6+.0001*(EJ-1000.))
      IF (JZ.GE. 5 .AND. IZ.GE. 5 .AND. IZ.LE.16)
     2CALL YIELD2 (IZ, IA, JZ, JA, EJ, QJ)
      IF (JZ.GE. 5 .AND. IZ.GE.17 .AND. IZ.LE.20)
     3CALL YIELD3 (IZ, IA, JZ, JA, EJ, QJ)
      IF (JZ.GE. 5 .AND. IZ.GE.21 .AND. IZ.LE.92)
     4CALL YIELD4 (IZ, IA, JZ, JA, EJ, QJ)
      IF (IZ.EQ.6.AND.IA.EQ.12.AND.JZ.EQ.4.AND.JA.EQ. 8)  QJ = QJ*1.8
      IF (IZ.EQ.7.AND.IA.EQ.14.AND.JZ.EQ.6.AND.JA.EQ.12)  QJ = QJ*1.8
      IF (IZ.EQ.7.AND.IA.EQ.14.AND.JZ.EQ.6.AND.JA.EQ.13)  QJ = QJ*0.5
      IF (IZ.EQ.8.AND.IA.EQ.16.AND.JZ.EQ.6.AND.JA.EQ.12)  QJ = QJ*1.8
      IF (IZ.EQ.8.AND.IA.EQ.16.AND.JZ.EQ.7.AND.JA.EQ.14)  QJ = QJ*1.8
      IF (IZ.EQ.8.AND.IA.EQ.16.AND.JZ.EQ.7.AND.JA.EQ.15)  QJ = QJ*1.5
        IF (IZ.EQ. 8 .OR.IZ.EQ.10) THEN                                 
          IF (IZ*2.EQ.IA.AND.JA-2*JZ.GE.2.AND.JZ.GE.5)    QJ=QJ*0.7     
        ENDIF                                                           
        IF (IZ.GE. 9.AND.IZ.LE.16.AND.2*JZ-JA.EQ.1.AND.JZ.GE.9)QJ=QJ*.7 
        IF (IZ.GE. 6.AND.IZ.LE.11.AND.2*JZ-JA.EQ.2.AND.JZ.GE.9)QJ=QJ*.4 
        IF (IZ.GE.10.AND.IZ.LE.13)  THEN                                
          IF ((JZ.EQ.6.OR.JZ.EQ.8).AND.JZ*2.EQ.JA)  QJ=QJ*2.            
        ENDIF                                                           
        IF((IZ.EQ.10.AND.IA.EQ.20).or.(IZ.EQ.12.AND.IA.EQ.24)) THEN     
          IF (JZ.EQ.7.AND.(JA.EQ.14.OR.JA.EQ.15))            QJ=QJ*1.5  
        ENDIF                                                           
        IF (IZ.GE.10.AND.IZ.LE.16.AND.JZ.EQ.9 )              QJ=QJ*0.8  
C
        IF ((IZ.EQ.12.OR.IZ.EQ.14).AND.(2*IZ.EQ.IA)) THEN               
          IF (IZ-JZ.EQ.2.AND.IA-JA.EQ.4)                     QJ=QJ*1.6  
          IF (IZ-JZ.EQ.1.AND.IA-JA.EQ.1)                     QJ=QJ*1.6  
        ENDIF                                                           
C
        IF ((IZ.EQ.16).AND.(IA.EQ.32)) THEN                             
          IF (IZ-JZ.EQ.2.AND.IA-JA.EQ.4)                     QJ=QJ*1.6  
          IF (IZ-JZ.EQ.1.AND.IA-JA.EQ.1)                     QJ=QJ*1.2  
          IF (IZ-JZ.EQ.1.AND.IA-JA.EQ.2)                     QJ=QJ*0.8  
          IF (IZ-JZ.EQ.3.AND.IA-JA.EQ.5)                     QJ=QJ*1.4  
          IF (IZ-JZ.EQ.4.AND.IA-JA.EQ.8)                     QJ=QJ*1.4  
          IF (IZ-JZ.EQ.8.AND.IA-JA.EQ.16)                    QJ=QJ*1.4  
        ENDIF                                                           
        IF ((IZ.EQ.18.OR.IZ.EQ.20).AND.2*IZ.EQ.IA)  THEN                
          IF (IZ-JZ.EQ.1.OR.IZ-JZ.EQ.3)                      QJ=QJ*0.7  
        ENDIF                                                           
        IF (IZ.EQ.20.AND.IA.EQ.40.AND.(JZ.EQ.12.OR.JZ.EQ.14))QJ=QJ*2.4  
        IF (IZ.EQ.20.AND.IA.EQ.40.AND.(JZ.EQ.18.OR.JZ.EQ.16))QJ=QJ*1.4  
        IF (IZ.GE.24.AND.IZ.LE.28) THEN                                 
          if (JZ.GE.20.AND.JZ.LE.23.and.JA-JZ*2.GE.6)        QJ=QJ*0.5  
        ENDIF                                                           
        IF((IZ.EQ.26.AND.IA.EQ.56).OR.(IZ.EQ.24.AND.IA.EQ.52)) THEN     
          IF((JZ.EQ.20).OR.(JZ.EQ.18).OR.(JZ.EQ.16))         QJ=QJ*1.3  
        ENDIF                                                           
        IF (IZ.EQ.28.AND.IA.EQ.58.AND.(JZ.ge.19.and.JZ.le.22))QJ=QJ*1.4
        F2A= 1.+.9*EXP(-((EJ-1230)/150)**2)*EXP(-(ABS(IZ-JZ-12)/5.)**2)
        F2B= 1.+.9*EXP(-((EJ-1230)/950)**2)*EXP(-(ABS(IZ-JZ-12)/5.)**2)
        IF (IZ.GE.30.AND.IZ-JZ.GE.6.AND.EJ.LT.1230.)  QJ = QJ*F2A       
        IF (IZ.GE.30.AND.IZ-JZ.GE.6.AND.EJ.GE.1230.)  QJ = QJ*F2B
C
        IF (IZ.EQ. 6.AND.JZ.GE. 3.AND.JZ.LE. 4) QJ=QJ*1.4               
C
      IF (IZ.EQ.26.AND.IA.EQ.56) THEN
          IF (JZ.EQ.23)          QJ=QJ*(1. - 0.6*EXP(-((52-JA)/2.6)**2))
          IF (JZ.EQ.24.AND.JA.EQ.54) QJ=0.7*QJ
          IF (JZ.EQ.25.AND.JA.GE.54.AND.JA.LE.55)
     1                           QJ=QJ*(1.7 - (JA-54)*0.45)
          IF (JZ.EQ.17)          QJ = QJ*0.9
        ENDIF
        IF (IZ.EQ.21.AND.JZ.GE.6)  THEN
            CALL YIELD3(IZ,IA,JZ,JA,EJ,Q3)
            QJ = SQRT(QJ*Q3)
        ENDIF
        IF (IZ.EQ.20.AND.JZ.GE.6)  THEN
            CALL YIELD4(IZ,IA,JZ,JA,EJ,Q4)
            QJ = SQRT(QJ*Q4)
        ENDIF
        IF (JZ .EQ. 5)  THEN
          CALL YIELD1(IZ, IA, JZ, JA, EJ, Q1)
          QJ = SQRT(Q1*QJ)
          IF (IZ.EQ.7.AND.IA.EQ.14.AND.JZ.EQ.5.AND.JA.EQ.10) QJ=QJ*1.8
          IF (IZ.EQ.6.AND.IA.EQ.12.AND.JZ.EQ.5.AND.JA.EQ.10) QJ=QJ*1.3
          IF (IZ.EQ.6.AND.IA.EQ.12.AND.JZ.EQ.5.AND.JA.EQ.11) QJ=QJ*1.3
          RETURN
        ENDIF
C
C       The followings were added just prior to 26th ICRC date
C
        IF (IZ.EQ.20.AND.IA.EQ.40.AND.JZ.EQ.16.AND.JA.EQ.31) QJ=QJ*0.7
        IF (IZ.EQ.20.AND.IA.EQ.40.AND.JZ.EQ.15.AND.JA.EQ.31) QJ=QJ*1.3
        IF (IZ.EQ.18.AND.IA.EQ.40.AND.JZ.EQ.16.AND.JA.EQ.31) QJ=QJ*1.2
        IF (IZ.EQ.18.AND.IA.EQ.40.AND.JZ.EQ.15.AND.JA.EQ.31) QJ=QJ*1.2
        IF (IZ.EQ.14.AND.IA.EQ.28.AND.JZ.EQ.12.AND.JA.EQ.26) QJ=QJ*1.3
        IF (IZ.EQ.12.AND.IA.EQ.24.AND.JZ.EQ.12.AND.JA.EQ.23) QJ=QJ*1.4
        IF (IZ.EQ.10.AND.IA.EQ.20.AND.JZ.EQ.10.AND.JA.EQ.19) QJ=QJ*1.3
        IF (IZ.EQ.10.AND.IA.EQ.20.AND.JZ.EQ. 9.AND.JA.EQ.19) QJ=QJ*1.2
        IF (IZ.EQ. 8.AND.IA.EQ.16.AND.JZ.EQ. 7.AND.JA.EQ.14) QJ=QJ*0.8
        IF (IZ.EQ. 8.AND.IA.EQ.16.AND.JZ.EQ. 6.AND.JA.EQ.14) QJ=QJ*1.2
        IF (IZ.EQ. 6.AND.IA.EQ.12.AND.JZ.EQ. 4.AND.JA.EQ.10) QJ=QJ*1.5
          RETURN
        END
      SUBROUTINE YIELD1 (IZ, IA, JZ, JA, EJ, QJ)
C     IZ .GE. 3  AND JZ = 2 TO 4
      COMMON /FS/ QR,QE,QF,QH,FE,FF,FA,FZ,PJ,PN,GA,ANZJ,AA,AE,AC,EC
      COMMON /QG/ QI,G1,G2,G3,G4
      COMMON /ST/ ST,SS,T
      DIMENSION  BA(92), CJ(56,26), C2(56), C3(56), C4(56), C5(56)
      DATA   BA / 16 * 0 ., 36., 38., 40., 44.
     1,           45., 48., 50., 52., 55.,55.7,58.5, 61., 64., 67.
     1,           70., 73., 75., 78., 80., 82., 84., 86., 89., 92.
     1,           93., 96., 98.,100.,103.,106.,108.,111.,114.,118.
     1,          122.,125.,127.,130.,133.,134.,137.,139.,141.,146.
     1,          146.,149.,153.,156.,159.,161.,165.,166.,169.,172.
     1,          175.,178.,181.,183.,186.,188.,192.,194.,197.,200.
     1,          204.,206.,209.,  0.,0.,0.,0. ,226.,227.,232.,231.,238./
      EQUIVALENCE(CJ(1, 2),C2),(CJ(1, 3),C3),(CJ(1, 4),C4),(CJ(1, 5),C5)
      DATA  C2/ 5*1.0,0.60,1.00,1.00,1.00,1.00,46*1. /          
      DATA  C3/ 5*1.0,1.00,1.80,0.70,3.00,1.00,46*1. /          
      DATA  C4/ 6*1.0,0.95,1.00,0.65,1.30,1.00,45*1. /
      DATA  C5/ 7*1.0,0.20,1.00,0.70,1.70,1.00,44*1. /          
      DATA  QC/ 13.0/, QL / 13.0/ , B / 1.15/
      DATA  PC/  0. /, PL / 0.16/
      DATA  RC/ 1.80/, RL / 10.7/ , RU/ 0.25/
      DATA  S / 0.54/, C  / 0.32/
      T  = .003
      QR = 0.
      QE = 0.
      QF = 0.
      QH = 0.
      AE = 0.
      AC = 0.
      FE = 1.
      FF = 1.
      FA = 1.
      FZ = 1.
      PN = 1.
      GA = 1.
      AI = IA
      ZI = IZ
      AJ = JA
      ZJ = JZ
      AA = AI - AJ
      AM = BA(IZ)
      IF (AM .EQ. 0)    AM = IA
      AN = AJ - ZJ
      ANZJ = AN/ZJ
      PJ = CJ(JA,JZ)
      CN = .3*(AI-AM)/ZI
      KN = (IA - IZ) - (JA - JZ)
      JP = IZ - JZ + 1
      JN = JA - JZ
!### imos      MN = JN .AND. 1 ! JN,JZ = odd, TRUE =1 ### imos 1/21/2000
!### imos      MZ = JZ .AND. 1 ! JN,JZ = even,FALSE=0 ### imos 1/21/2000
      MN = 1                          ! ### imos 1/21/2000
      if (JN/2*2 .eq. JN) MN = 0      ! ### imos 1/21/2000
      MZ = 1                          ! ### imos 1/21/2000
      if (JZ/2*2 .eq. JZ) MZ = 0      ! ### imos 1/21/2000
!
      PN = 1.15
      IF (MZ .EQ. 1.AND.MN.LE.1)   PN = 0.9 - 0.1*MN
      EC = 68.7*AI**.866                                        
      IF (EC.LT.1250.)  EC = 1250.
      EI = EJ
      IF (EI.GT.EC)   EI = EC
      H3 = 1.
      IF (EI.LT.80.)  H3 = 1. - EXP(-(EI/25.)**4)
      IF((JP.GT. 2).OR. (KN.GT. 1))     GO TO  3
      CX = 1.
      CT = 1.                                                   
      IF (JZ.EQ.4 .AND. JA.EQ.9)   CX = CJ(JA,JZ)/PN
      IF (AI/ZI.GT.2.)  CT = (ZI-2.)/((AI-ZI)-2.)
      IF (ZI.GT.5 .OR. CT.GT.1. .OR. AI-ZI.LE.2.)   CT = 1.
 1    IF (JP.GT. 1)      GO TO  2
      IF (KN.GT. 1)                     GO TO  3
      IF (AI .LE. 40.)   QN = 24.*(1. + .01*AI)
      IF (AI .GT. 40.)   QN = 1.02*(AI - 7.)
      IF (AI .GE. 63.)   QN = 57.
      IF (EI.LT.2500.)
     1FE = (1.0 + 2.1*EXP(-(EI/100.)**2) + 0.4*EXP(-EI/350.))*H3
      QH = QN * CX
      QJ = QH
      IF (EI.LT.EC)      QJ = QH*FE
      RETURN
 2    IF (KN.NE. 0)                     GO TO  3
      IF (EJ.LT.2500.)
     1FE = (1.0 - EXP(-(EJ/230.)**2) + 2.2*EXP(-EJ/75.)
     1   + .33*EXP(-((EJ-900.)/500.)**2))*H3
      QH = 21.*CX*CT
      QJ = QH*FE
      RETURN
 3    IF (IZ.LT.29)      GO TO  5
      IF (JZ.NE. 4)      FE = (EI/EC)**(.4*ZJ)
      IF (JZ.EQ. 4)      FE = (EI/EC)**1.8
      FA = EXP(.01*(AI - 56.)*((AN/ZJ - CN) - .45))
      IF (FA.LE.1.)      FA = 1.
      EI = EC
      ZI = 29.
      AI = 63.
 5    AT = IA
      IF (IA.GE.104)  AT = 104
      IF (IA.GE.64)   GO TO 7
      IF (IA.GT.34)   GO TO 6
      IF (IA.GE.14.AND.JA.EQ. 6.AND.JZ.EQ. 2)  FZ = 1. + .1*(IA-14)
      GO TO 8
 6    IF (JA.EQ. 6.AND.JZ.EQ. 2)  FZ = 3.0*(1. + .02*(IA-34))
      IF (JA.EQ. 6.AND.JZ.EQ. 3)  FZ = 1.0 + .02*(IA-34)
      IF (JA.EQ. 7.AND.JZ.LE. 4)  FZ = 1.0 + .01*(IA-34)
      GO TO 8
 7    IF (JA.EQ. 6.AND.JZ.EQ. 2)  FZ = 4.8 + 0.045*(AT-64.)
      IF (JA.EQ. 6.AND.JZ.EQ. 3)  FZ = 1.6 + 0.015*(AT-64.)
      IF (JA.EQ. 7.AND.JZ.EQ. 3)  FZ = 1.3 + 0.015*(AT-64.)
      IF (JA.EQ. 7.AND.JZ.EQ. 4)  FZ = 1.3 + .0105*(AT-64.)
      IF (JA.EQ. 8.AND.JZ.EQ. 3)  FZ = 1.0 + .0225*(AT-64.)
      IF (JA.EQ. 9.AND.JZ.EQ. 3)  FZ = 1.0 + .0125*(AT-64.)
      IF (JA.EQ. 9.AND.JZ.EQ. 4)  FZ = 1.0 + .01  *(AT-64.)
      IF (JA.GE.10.AND.JZ.LE. 5)  FZ = 1.0 + .0050*(AT-64.)
8     DA = AI - AJ
      AA = DA
      AE = 31.5 + .052*(AI - 36.0)*(ALOG(EI) - 3.17)
      IF (AA.GT.AE)      DA = AE
      Q1 = QC
      IF (EI .LT. 1250.)  Q1 = QL*EXP(B*(1. - .0008*EI))
      PE = PC
      IF (EI .LT. 1250.) PE = PL*(1. - .0008*EI)
      RE = RC
      IF (EI .LT. 1250.) RE = RL/EI**RU
      IF((IZ.GT.20).AND.(EI.LT.1250.))   RE = 1.8
      IF((EI.LT.1250.).AND.(IZ.GE.21)) Q1=Q1*2*EXP(-((EI-650)/720)**2)
      IF (AI/ZI .GE. 2.)   SS = S - C*(AI/ZI - 2.)**1.4
      IF (AI/ZI .LT. 2.)   SS = S + C*(2. - AI/ZI)**1.4
      ST = (ZJ - (SS - T*AJ)*AJ)
      ZA = (ZJ - (SS - T*AJ)*AJ)**2
      QR = Q1*EXP(-PE*DA - RE*ZA)*FA*FZ*CJ(JA,JZ)
      QH = QC*EXP(-RC*ZA)*FA*FZ*CJ(JA,JZ)
      QE = QH*FE
      IF (IZ.LT.29)      QJ = QR
      IF (IZ .GE. 29)   QJ = QE
      IF (JP.EQ.3 .AND. KN.EQ.0)   QJ = 1.5*QR
C
C     FOLLOWING REVISIONS ARE ADDED 12/27/78
      IF (IZ .LE. 20)   RETURN
      QI = QJ
      EX = 68.7*56.**.866*EJ/EC
      AC = 31.5 + .045*(AI-36.)*(ALOG(AI)+1.23)
      G1 = 1.-0.6*(1.-EXP(-(EJ/1000.)**2))*EXP(-(EJ/2000.)**2)
     1       +0.2*(1.-EXP(-(EJ/3000.)**2))
      GX = 1.-0.6*(1.-EXP(-(EX/1000.)**2))*EXP(-(EX/2000.)**2)
     1       +0.2*(1.-EXP(-(EX/3000.)**2))
      G2 = 1.-0.4*(1.-EXP(-(EX/2000.)**2))*EXP(-((EX-1800.)/1800.)**2)
     1       +.17*(1.-EXP(-(EX/2000.)**2))
      IF (EJ.GT.2500.)   G1=GX
      IF (EJ.GT.EC .AND. EJ.LT.2500.)   G1=SQRT(G1*GX)
      IF (IZ.GE.21 .AND. IZ.LE.28)   QJ = QR*G1
      IF (IZ.GT.28 .AND. AA.GT.AC)   QJ = QE*G2
C
      END
      SUBROUTINE YIELD2 (IZ, IA, JZ, JA, EJ, QJ)
C     5 .LE. IZ .LE. 16  AND  5 .LE. JZ .LE. IZ
      COMMON /FS/ QR,QE,QF,QH,FE,FF,FA,FZ,PJ,PN,GA,ANZJ,AA,AE,AC,EC
      COMMON /ST/ ST,SS,T
      DIMENSION  CJ(56,26)
      DIMENSION  C5(56), C6(56), C7(56), C8(56), C9(56), CD(56)
      DIMENSION  D1(56), D2(56), D3(56), D4(56), D5(56), D6(56)
      EQUIVALENCE(CJ(1, 5),C5),(CJ(1, 6),C6),(CJ(1, 7),C7),(CJ(1, 8),C8)
      EQUIVALENCE(CJ(1, 9),C9),(CJ(1,10),CD),(CJ(1,11),D1),(CJ(1,12),D2)
      EQUIVALENCE(CJ(1,13),D3),(CJ(1,14),D4),(CJ(1,15),D5),(CJ(1,16),D6)
      DATA  C5/ 7*1.0,0.20,1.00,0.70,1.70,1.00,44*1.0/          
      DATA  C6/ 9*1.0,0.60,46*1.0/                              
      DATA  C7/12*1.0,0.39,1.00,2.00,41*1.0/                    
      DATA  C8/14*1.0,1.20,1.00,1.00,39*1.0/, C9/17*1.0,1.00,38*1.0/
      DATA  CD/18*1.0,0.60,1.00,0.83,35*1.0/, D1/21*1.0,1.20,34*1.0/
      DATA  D2,D3,D4,D5,D6/56*1.0,56*1.0,56*1.0,56*1.0,56*1.0/
      DATA  PC,PL,PU/.075,2.60,0.50/, S,C/.502,.26/
      DATA  RC,RL,RU/1.60,10.2,0.26/
      T  = .0005
      QR = 0.
      QE = 0.
      QF = 0.
      QH = 0.
      AE = 0.
      AC = 0.
      FE = 1.
      FF = 1.
      FA = 1.
      FZ = 1.
      GA = 1.
      AI = IA
      ZI = IZ
      AJ = JA
      ZJ = JZ
      AA = AI - AJ
      AN = AJ - ZJ
      ANZJ = AN/ZJ
      PJ = CJ(JA,JZ)
      PN = 1.15
      JN = JA - JZ
!### imos      MN = JN .AND. 1 ! JN,JZ = odd, TRUE =1 ### imos 1/21/2000
!### imos      MZ = JZ .AND. 1 ! JN,JZ = even,FALSE=0 ### imos 1/21/2000
      MN = 1                          ! ### imos 1/21/2000
      if (JN/2*2 .eq. JN) MN = 0      ! ### imos 1/21/2000
      MZ = 1                          ! ### imos 1/21/2000
      if (JZ/2*2 .eq. JZ) MZ = 0      ! ### imos 1/21/2000
!
      CX = 1.
      IF((MZ.EQ. 1).AND.(MN.LE. 1))   PN = 0.9 - 0.1*MN
      EC = 68.7*AI**.866                                        
      IF (EC.LT.1250.)  EC = 1250.
      EI = EJ
      IF (EI.GT.EC)  EI = EC
      H3 = 1.
      IF (EI.LT.80.)  H3 = 1. - EXP(-(EI/25.)**4)
      KN = (IA - IZ) - (JA - JZ)
      JP = IZ - JZ + 1
      IF((JZ.EQ.7.AND.JA.EQ.13).OR.(JZ.EQ.10.AND.JA.EQ.19)) CX=CJ(JA,JZ)
      CT = (ZI - 2.)/((AI - ZI) - 2.)
      IF (ZI.GT.5. OR. CT.GT.1.)   CT = 1.
      IF (JP.GT. 1)      GO TO  2
 1    IF (KN.NE. 1)      GO TO  3
      IF (EI.LT.2500.)
     1FE =(1.0 + 2.1*EXP(-(EI/100.)**2) + 0.4*EXP(-EI/350.))*H3
      QH = 24.*(1.+ .01*AI)*CX
      QJ = QH
      IF (EI.LT.EC)      QJ = QH*FE
      GO TO 10                                                  
 2    IF((JP.GT. 2).OR. (KN.NE. 0))     GO TO  3
      IF (EJ.LT.2500.)
     1FE =(1.0 - EXP (-(EJ/230.)**2) + 2.2*EXP(-EJ/75.)
     1    + .33*EXP(-((EJ-900.)/500.)**2))*H3
      QH = 21.*CX*CT
      QJ = QH*FE
      GO TO 10                                                  
 3    F1 = 1. - .3*ALOG(AI/20.)
      PE = PC
      IF (EC .GT. 1250.)  PE = 0.77/AI**(2./3.)
      IF (EI.LT.EC)  PE = PL/EI**PU
      RE = RC
      IF (EI .LT. 1250.)  RE = RL/EI**RU
      IF (AI/ZI .GE. 2.)   SS = S - C*(AI/ZI - 2.)**1.4
      IF (AI/ZI .LT. 2.)   SS = S + C*(2. - AI/ZI)**1.4
      ST = (ZJ - (SS - T*AJ)*AJ)
      ZA = (ZJ - (SS - T*AJ)*AJ)**2
      Q1 = 27.6*(AI**(2./3.)-1.)*F1*PE*RE**.5/(1. - EXP (-PE*AI))
      QC = 27.6*(AI**(2./3.)-1.)*F1*PC*RC**.5/(1. - EXP (-PC*AI))
      QR = Q1*EXP(-PE*(AI - AJ))*EXP(-RE*ZA)*CJ(JA,JZ)*PN
      QJ = QR
      QH = QC*EXP (-PC*(AI - AJ))*EXP (-RC*ZA)*CJ(JA,JZ)*PN
      IF((JP.LT. 3).OR. (KN.NE. 0))     GO TO 10                
      QJ = QR*AMIN1(.0022*AJ*AJ,1.)                                     
      IF (JP.NE. 4)      GO TO 10                               
      IF (QJ.GT..5)      QJ = .5
 10   IF (IZ.GE.14 .AND. EI.GE.500.)    QJ=QJ*(1.+.12*(IZ-13))  
      IF (IZ.GE.14 .AND. EI.GE.200. .AND. EI.LT.500.)           
     &    QJ = QJ*(1.+.12*(IZ-13)*EXP(-((EI-500)/350)**2))      
      IF (JP.EQ. 2 .AND. KN.EQ. 1 .AND. IZ.GE.12)  QJ = QJ*1.7
!### imos      IF (AJ/ZJ .GT. 1.8 .OR. IZ.GT.10) RETURN        
      IF (AJ/ZJ-1.e-7 .GT. 1.8 .OR. IZ.GT.10) RETURN     ! ### imos 1/21/2000
      IF (JP.NE. 1 .OR.  KN.NE. 2)   QJ = QJ*0.3                
      IF (JP.EQ. 1 .AND. KN.EQ. 2)   QJ = QJ*0.5                
      END
      SUBROUTINE YIELD3(IZ, IA, JZ, JA, EJ, QJ)
C     IZ = 17 - 20  AND  JZ = 5 TO IZ
      COMMON /FS/ QR,QE,QF,QH,FE,FF,FA,FZ,CJ,PN,GA,ANZJ,AA,AE,AC,EC
      COMMON /ST/ ST,S,T
      DIMENSION   BA(92)
      DATA   BA / 16 * 0 ., 36., 38., 40., 44.
     1,           45., 48., 50., 52., 55.,55.7,58.5, 61., 64., 67.
     1,           70., 73., 75., 78., 80., 82., 84., 86., 89., 92.
     1,           93., 96., 98.,100.,103.,106.,108.,111.,114.,118.
     1,          122.,125.,127.,130.,133.,134.,137.,139.,141.,146.
     1,          146.,149.,153.,156.,159.,161.,165.,166.,169.,172.
     1,          175.,178.,181.,183.,186.,188.,192.,194.,197.,200.
     1,          204.,206.,209.,  0.,0.,0.,0. ,226.,227.,232.,231.,238./
        AI = IA
        ZI = IZ
        AJ = JA
        ZJ = JZ
        AM = BA(IZ)
        IF (AM .EQ. 0.)  AM = IA
        EC = 68.7*AI**.866                                      
        IF (EC.LT.1250.)   EC = 1250.
        EI = EJ
        IF (EI.GT.  EC )   EI = EC
        QR = 0.
        QE = 0.
        QF = 0.
        QH = 0.
        AE = 0.
        AC = 0.
        DX=0
        FE = 1.
        FF = 1.
        FA = 1.
        FZ = 1.
        GA = 1.
        CJ = 1.
        AA = AI - AJ
        H3 = 1. - EXP(-(AMIN1(EI,80.)/25.)**4)
        H4 = AMIN1(AMAX1((400./EI-2.2),1.),2.)
        AN = AJ - ZJ
        ANZJ = AN/ZJ
        JN = JA - JZ
!### imos      MN = JN .AND. 1 ! JN,JZ = odd, TRUE =1 ### imos 1/21/2000
!### imos      MZ = JZ .AND. 1 ! JN,JZ = even,FALSE=0 ### imos 1/21/2000
      MN = 1                          ! ### imos 1/21/2000
      if (JN/2*2 .eq. JN) MN = 0      ! ### imos 1/21/2000
      MZ = 1                          ! ### imos 1/21/2000
      if (JZ/2*2 .eq. JZ) MZ = 0      ! ### imos 1/21/2000
!
        PN = 1.
        IF (MN+MZ .EQ. 2)   PN = 0.85
        IF (MN+MZ .EQ. 0)   PN = 1.25
        IF (MN-MZ .EQ. 1)   PN = 0.90
        JP = IZ - JZ+ 1
        KN = IA - IZ - JN
        IF((JP.GT. 2).OR. (KN.GT. 3).OR. (IA.LT.35))   GO TO 14
        IF (JP.EQ. 2)   GO TO  2
 1    IF (EI.LT.2500.)
     1FE =(1.0 + 2.1*EXP(-(EI/100.)**2) + 0.4*EXP(-EI/350.))*H3
        QH = 24.*(1. + .01*AI)
        QJ = QH
        DX = 3.
        IF (KN .GT. 1)   DX = 15.
        IF (EI.LT.EC)   QJ = QH*FE
        IF (KN.EQ. 1)   GO TO  5
        XN = KN
        XA = 1.17                                                       
        IF (IZ.LE.30)  XA=1.6                                           
        QH = QH*EXP(1. - XN**(XA-.0048*AI))
        QJ = QH
        IF (EI.GE.EC)   GO TO  5
        GO TO  4
 2    IF (EJ.LT.2500.)
     1FE =(1.0 - EXP (-(EJ/230.)**2) + 2.2/EXP (EJ/75.)
     1    + .33/EXP (((EJ-900.)/500.)**2))*H3
        QH = 21.
        QJ = QH
        QJ = QH * FE
        IF (KN.EQ. 0)   DX = -3.
        IF (KN.EQ. 1)   DX = -1.
        IF (KN.EQ. 0)   GO TO  5
        IF (KN.GT. 2)   GO TO 14
        IF (IA.GE.35.AND.IA.LE.40.AND.KN.GT.1)  GO TO 14
        QH = 17.
        QJ = QH
        IF (EI.GE.EC)   GO TO  5
 4      KA = IA - JA
        IF (KA.GE. 8)   GO TO  5
        FE = (1. + 1.9/EXP((AA/7.9)**2 + (EI/420.)**1.4))*H3*H4
        QJ = QH*FE
 5      DD = DX*(AI - AM)/ZI
        YA = EXP(DD)
        IF (DD.GT. 0)   YA = 2. - 1./YA
        IF((IA.LT.35).OR. (IA.GT.209))   YA = 1.
        IF((IA.LT.70).AND.(JP.EQ. 3 ))   YA = 1.
        QH = YA*QH
        QJ = YA*QJ
        RETURN
 14     DA = AI - AJ
        IF((JZ.EQ. 7).AND.(JA.EQ.13))     CJ = 0.39
        F1 = 1. - .3*ALOG(AI/20.)
        IF (EI.GE.600.)   F2 = 1.
        IF (EI.LT.600.)   F2 = EXP (.90 - .0015*EI)
        IF (EI.GE.EC)   GO TO 15
        PE = 20./EI**.77
 15     IF (EI.LT.EC)   GO TO 16
        F2 = 1.00
        PE = 1.98/AI**.92                                       
 16     IF (EI.LT.1250.)   RE = 10.2/EI**0.26
      IF (EI.GE.1250.)   RE = 1.60
      Q1 = 27.6*(AI**(2./3.)-1)*F1*F2*PE*RE**.5/(1. - EXP (-PE*AI))
      S  = 0.502 - 0.08*ABS (AI/ZI - 2.)
      T  = 0.0005
      ST = (ZJ - (S - 0.0005*AJ)*AJ)
      ZA = ABS(ZJ - (S - 0.0005*AJ)*AJ)**2
      QJ = Q1*EXP (-PE*DA-RE*ZA)
      QR = QJ*PN*CJ
      QJ = QR
      PH = .77/AI**.667
      RH = 1.60
      QC = 27.6*(AI**(2./3.)-1)* F1*PH*RH**.5/(1. - EXP (-PH*AI))
      QH = QC*EXP (-PH*DA - RH*ZA)*PN*CJ
      IF((JP.NE. 3).OR. (KN.GT. 0).OR. (IA.LT.35))   GO TO 20   
      FF = 1.
      FE = 1.
      QJ = QH
 20   IF (IZ.LT.14 .OR. IZ .GT.19 .OR. EI.LT.200.)   RETURN     
      IF (EI.GT.500.)  QJ = QJ*(1.+.12*(IZ-13))
      IF (EI.LE.500.)  QJ = QJ*(1.+.12*(IZ-13)*EXP(-((EI-500)/350)**2)) 
      END
      SUBROUTINE YIELD4(IZ, IA, JZ, JA, EJ, QJ)                 
C     IZ .GE. 21  AND  JZ .GE. 5                        
      COMMON /FS/ QR,QE,QF,QH,FE,FF,FA,FZ,CJ,PN,GA,ANZJ,AA,AE,AC,EC
      COMMON /Q4/ C0, F1, F2, F3, PE, DA, RJ, ZA, PH, YA, Q1, QM, HM
      COMMON /QG/ QI, G1, G2, G3, G4
      COMMON /ST/ ST,S,T2
      DIMENSION   BA(92)
      DATA   BA / 16 * 0 ., 36., 38., 40., 44.
     1,           45., 48., 50., 52., 55.,55.7,58.5, 61., 64., 67.
     1,           70., 73., 75., 78., 80., 82., 84., 86., 89., 92.
     1,           93., 96., 98.,100.,103.,106.,108.,111.,114.,118.
     1,          122.,125.,127.,130.,133.,134.,137.,139.,141.,146.
     1,          146.,149.,153.,156.,159.,161.,165.,166.,169.,172.
     1,          175.,178.,181.,183.,186.,188.,192.,194.,197.,200.
     1,          204.,206.,209.,  0.,0.,0.,0. ,226.,227.,232.,231.,238./
      DATA    S0,S1,T3/.482,.07,3.E-7/,          R0,R1/11.8,0.45 /
      DATA    P0,P1   /1.98, 0.92/,              E0,E1/20.3,1.169/
      DATA    C1,C2,C3,C4/144.,.367,0.3,0.7/,    D1,D2/.0365,1.23/
      T2 = 2.8E-4
      AI = IA
      ZI = IZ
      AJ = JA
      ZJ = JZ
      AM = BA(IZ) + IA*(300/(INT(BA(IZ))+300))
      KM = 0
      FP = 1.                                                           
      QJ = 0.
      QR = 0.
      QE = 0.
      QF = 0.
      QM = 0.
      HM = 0.
      QH = 0.
      AC = 0.
      AE = 0.
      DM = 0.
      DX = 0.
      KFLAG = 0
      FE = 1.
      FF = 1.
      FM = 1.
      FA = 1.
      FZ = 1.
      GA = 1.
      CJ = 1.
      PN = 1.
      YA = 1.
      IF((IA.LE.JA).OR.(IZ.LT.JZ).OR.(IA-JA.LT.IZ-JZ))   RETURN
      AA = AI - AJ
      DA = AA
      AN = AJ - ZJ
      ANZJ = AN/ZJ
      ANZI = AI/ZI - 1.
      JN = JA - JZ
      JP = IZ - JZ + 1
      CN = .3*(AI - AM)/ZI
      KN = (IA - IZ) - JN
      XN = KN
!### imos      MN = JN .AND. 1 ! JN,JZ = odd, TRUE =1 ### imos 1/21/2000
!### imos      MZ = JZ .AND. 1 ! JN,JZ = even,FALSE=0 ### imos 1/21/2000
      MN = 1                          ! ### imos 1/21/2000
      if (JN/2*2 .eq. JN) MN = 0      ! ### imos 1/21/2000
      MZ = 1                          ! ### imos 1/21/2000
      if (JZ/2*2 .eq. JZ) MZ = 0      ! ### imos 1/21/2000
!
      IF (MN+MZ .EQ. 2)   PN = 0.85
      IF (MN+MZ .EQ. 0)   PN = 1.25
      IF (MN-MZ .EQ. 1)   PN = 0.90
      PN = PN + (1.-PN)*(1.-EXP(-((IA-100)/35.)**2))                    
      AC = 31.5 + .045*(AI-36.)*(ALOG(AI)+1.23)
      EI = EJ
      EC = E0*AI**E1                                                    
      IF (EC.GT.4000.)   EC = 4000.
      IF (EC.LT.1250.)   EC = 1250.
      IF (EI.GT.  EC )   EI = EC
      FMAX = 1800./EI
      IF (IZ.GE.84)   FMAX = (1800./EI)**(6.56-.067*ZI)
      IF (FMAX .GT. 4.)   FMAX = 4.
      IF (EI.GT.1800.)   FMAX = 1.
      IF (IZ.LE.28 .AND. JZ.EQ.20 .AND. JA.EQ.19)   CJ = 0.6
      IF((JZ.EQ. 7).AND.(JA.EQ.13))   CJ = 0.39
      IF((JZ.EQ.53).AND.(IZ.GE.71).AND.(IZ.LE.83))   CJ = 1. 
      P500 = 20./( 500.**.77)
      P1000= 20./(1000.**.77)
      P2000= 20./(2000.**.77)                                           
      P3000= 20./(3000.**.77)
      IF (IZ.EQ.92.AND.JZ.EQ.88)  QH = 1.20*EXP(-.70*ABS(JA-224.))
      IF (IZ.EQ.92.AND.JZ.EQ.89)  QH = 1.60*EXP(-.15*ABS(JA-224.))
      IF (IZ.EQ.92.AND.JZ.EQ.90)  QH = 8.00*EXP(-.25*ABS(JA-233.))
      IF (IZ.EQ.92.AND.JZ.EQ.91)  QH = 18.5*EXP(-.55*ABS(JA-234.))      
      IF (IZ.EQ.92.AND.JZ.EQ.92)  QH = 55.0*EXP(-.80*ABS(JA-237.))
      IF (JP .GE. 6 .OR. IA .LT. 35)  GO TO 100
      if (JP.GE.4.AND.JP.LE.5.AND.KN.GE.KM)  GO TO 100                  
C
C     PERIPHERAL
C
      H3 = 1. - EXP(-(EI/(15. + 10.*AA))**4)
      H4 = 400./EI - 2.2
      IF (H4 .LT. 1.)   H4 = 1.
      IF (H4 .GT. 2.)   H4 = 2.
      HE = (1. + 1.9*EXP(-(DA/7.9)**2 - (EI/420.)**1.4)
     1         -(1.-EXP(-(DA/12.)**8))*EXP(-(EI/420.)**3))*H3*H4
      EC = 2500.
      IF((JP.GE. 3).AND.(JP.LE. 5))   GO TO  3
      IF (JP.EQ. 2)      GO TO  2
 1    IF (EI.LT.2500.)
     1FE =(1.0 + 2.1*EXP(-(EI/100.)**2) + 0.4*EXP(-EI/350.))*H3
      IF (KN .GT. 1)   FE = HE
      KM = AI/20. + 1.5*(ABS(238.-AI)/167.)**2.5 + .8
      IF (IA.LE.70)   KM = 3 + AI/66.
      IF (KN.GT.KM)      GO TO 100                                      
      IF (IA .LE. 40)  QA = 24.0*(1. + .01*AI)
      IF (IA .GT. 40)  QA = 1.02*(AI - 7.)
      IF (IA .GE. 63)  QA = 57.
      QH = QA
      XA = 1.17                                                         
      IF (IZ .LE. 30)  XA = 1.5                                         
      IF (KN .NE. 1 )  QH = QA*EXP(1.-XN**(XA-.0048*AI))
      IF (KN .EQ. 1 )  DX = 3.
      IF (KN .GT. 1 )  DX = 15.
      QH = QH*(1. + 0.15*(IZ/80.)**2)                                   
      GO TO 6
 2    FE =(1.0 - EXP(-(EJ/230.)**2) + 2.2*EXP(-EJ/75.)
     1    + .33*EXP(-((EJ-900.)/500.)**2))*H3
      IF (KN .GT. 0 )   FE = HE
      KM = 10.*(1. - EXP(-((AI-39.)/54.5)**2))
      IF (KM.LT. 2)      KM = 2
      IF (IA.GE.35.AND.IA.LE.40.AND.KM.GT.1)  KM = 1                    
      IF (KN.GT.KM)      GO TO 100                                      
      QH = 17.
        IF (KN.GT. 0 .AND. IZ.LE.50)  QH = 27
        IF (KN.GT. 0 .AND. IZ.GT.50)  QH = 21
        IF (KN.EQ. 0 .AND. IZ.LE.28)  QH = 33
        IF (KN.EQ. 0 .AND. IZ.LE.23)  QH = 28
        IF (KN.EQ. 0 .AND. IZ.LE. 5)  QH = 21
      QJ = QH
      IF (KN.EQ. 0)   DX = -3.
      IF (KN.EQ. 1)   DX = -1.
      GO TO 6
 3    FN = AMIN1((EI/420.)**(.72 - .18*XN)*H3,1.)
      IF (IA .LE. 70)   FE = 1. - EXP(-(EI/35.)**4)
      IF (IA .GT. 70 .AND. KN .LE. 4)   FE = FN
      IF (IA .GT. 70 .AND. KN .GT. 4)   FE = H3
      IF (EI .LT. 200. )   FE = FE*EI/200.
      KM = AI/25. + 0.5
      XM = KM
      IF (KN .GT. KM)   KFLAG = 1
      FM=FE
      IF (KN .LE. 4 )  FM=(EI/420.)**(.72-.18*XM)*H3
      IF (FM .GT. 1.)   FM = 1.
      IF (EI .LT. 200. )   FM = FM*EI/200.
      QH = .2 + 60.*(EXP(-((AI - 89.)/25.)**2-((XN-4.6)/2.)**2)
     1   + (1.-EXP(-(AI/135.)**3))*EXP(-((XN-AI**.46)/AI**.27)**2))
      HM = .2 + 60.*(EXP(-((AI-89.)/25.)**2-((XM-4.6)/2.)**2)
     1   + (1.-EXP(-(AI/135.)**3))*EXP(-((XM-AI**.46)/AI**.27)**2))
      IF (KN.EQ. 0)      DX = -10.
      IF (KN.EQ. 1)      DX = -3.
      IF (IA.LE.70)  GO TO 100
      IF (KN.LE. 2)  QH = (1.2E5/AJ**2.2)*0.12*EXP(XN/0.85)             
      IF (KN.GE. 3 .AND. IZ.LE.26)  QH = (1.2E5/AJ**2.2)*EXP((XN-2)/2.)
 4    IF ((JZ.LT.82.AND.KN.GT.KM) .OR. (IZ.GE.88.AND.JP.GT.3))   KFLAG=1
      QH = 0.1**(JP-3)*QH
      HM = 0.1**(JP-3)*HM
 6    IF (EI .GE. EC)   FE = 1.
      IF (DX*(AI-AM) .LE. 0)  YA = EXP(DX*(AI-AM)/ZI)
      IF (DX*(AI-AM) .GT. 0)   YA = 2. - EXP(-DX*(AI-AM)/ZI)
      IF (JP.GE.3 .AND. IA.LE.70)   YA = 1.
 10   IF (IA.GT.157.AND.EI.GT.500.AND.JP*KN.NE.1)
     1YA = YA*(1. - .012*(AI-157.)*(1. - EXP(-((EI-500.)/290.)**2)))
     1       *(1. - EXP(-0.6*(1.*JP)**.8-((92.-ZI)/2.7)**2))
      QH = QH*YA
      IF (IZ.EQ.92.AND.JZ.EQ.88)  QH = 1.20*EXP(-.70*ABS(JA-224.))
      IF (IZ.EQ.92.AND.JZ.EQ.89)  QH = 1.60*EXP(-.15*ABS(JA-224.))
      IF (IZ.EQ.92.AND.JZ.EQ.90)  QH = 8.00*EXP(-.25*ABS(JA-233.))
      IF (IZ.EQ.92.AND.JZ.EQ.91)  QH = 18.5*EXP(-.55*ABS(JA-234.))      
      IF (IZ.EQ.92.AND.JZ.EQ.92)  QH = 55.0*EXP(-.80*ABS(JA-237.))
      IF (JP.GT.5 .OR. (JP.GT.3.AND.IZ.GE.88))  GO TO 100               
      QJ = QH*FE
      HM = HM*YA
      QM = HM*FM
        IF (KFLAG.EQ.0)  GO TO 200                                      
 100  EC = E0*AI**E1                                            
      IF (EC.GT.4000.)  EC=4000.
      IF (IA.GT.100)  FP = 1. - 0.20*(IA-100)/100.                      
      IF (IA.GT.180)  FP = 0.84                                         
      ANZX = AMIN1(ANZJ,ANZI)
      ANZC = ANZJ - CN
      IF (ANZC .GT. 1.5)  ANZC = 1.5
      AZ = ANZJ
      EE = 450./EI
      IF (JA.GT.36 .AND. JZ.LT.88 .AND. IZ.EQ.92)
     1  FF = EXP(33.*((ANZJ - 1.36)+.00006*EI-1200./(EI+80.)**2))
      IF (JA.GT.36 .AND. JZ.LT.88 .AND. IZ.EQ.90)
     1  FF = EXP(16.*EE**.6*((ANZJ-CN)+.00006*EI - 1200./(EI+80)**2))
      IF (JA.GE.36 .AND. JZ.LT.IZ .AND. IZ.GE.84 .AND. IZ.LT.90)
     1  FF = EXP(19.*EE**.6*((ANZX-CN)-ALOG(5.5/EI**.07)))
      GE = EE
      H1 = AJ - 105. + .68*(207.-AI)
      IF (EI .LE. 200.)   GE = 2.25
      IF (EI .GT. 450. .AND. EI .LT. 600.)  H1 = .0006*H1*(600. - EI)
      IF (EI .LE. 450.)                     H1 = .0900*H1
      IF (EI .GE. 600. .OR.  AI .LT. 180. .OR . H1 .LT. 0)   H1 = 0
      AE=.0065*(207-IA)*EXP(-ABS(EI-700.)/700.)
      IF (IA.GE.110 .AND. IZ.LE.83 .AND. JZ.LT.IZ)
     1FF = EXP(19.*EE**.6*((ANZX-CN)-ALOG(5.5/EI**.07) - AE)
     1         -(((AJ/AI-.46)/.15)**2*GE+H1))
      IF (IA.GE.110 .AND. IZ.LE.83 .AND. JZ.LT.IZ .AND. ANZX-CN.LT.1.3)
     1    FF = FF*EXP((1.3-ANZX+CN)*(350./EI)**4/(1. + (130./EI)**4))
      IF (FF.GT.FMAX)    FF = FMAX
      IF (JZ.GE.88)      GO TO 110
      EE = EI/EC
      A4 = (AI/63.)**4
      IF (EE.GE..3)      FE = EE**(2.2 - .01*A4)
      IF (EE.LT..3)      FE = EE**(2.2 + .01*A4)*.3**(-.02*A4)
      IF ((JZ.GE. 5).AND.(JZ.LE. 8))   FE = EE**1.8
      IF ((JZ.GE. 9).AND.(JZ.LE.10))   FE = EE**2.0
      GF = 1.
      EX = EI
      IF (EX .LT. 500.)  EX = 500.
      IF((JA.GE.18).AND.(JA.LE.30))
     1GF = EXP(2600.*((ANZX-CN)+.01*AI-3.1)/EX)
      IF((JA.GE.31).AND.(JA.LE.56))
     1GF = EXP(5600.*((ANZX-CN)+.0035*AI-1.72)/EX)
      IF(GF .LT. 1.)  GF = 1.
      FE = FE*GF
      IF((JZ.LT.30).OR.  (IZ.LE.83))      GO TO 140
C
C     URANIUM GROUP, ZP .GT. CU
C
 110  FU=1.
      IF (IZ.GE.88 .AND. JZ.GE.82 .AND. (KN.GT.KM .OR. JP.GT.3) )
     1FU=(1.-EXP(-.6*(1.*JP)**.8-((92.-ZI)/2.7)**2))
     1  *(1.-0.8*(EXP(-((ZI-ZJ)/5.4)**8)))
      IF (IZ.GE.88 .AND. JZ.GE.57 .AND. JZ.LT.82)
     1FU=1.0-0.36*EXP(-((82.-ZJ)/27.)**2-((92.-ZI)/2.7)**2)
      CJ=CJ*FU
 120  IF((JZ.GE.30).AND.(JZ.LE.35))
     1  QH =  6.5*EXP(- 80.*(AN/ZJ - 1.26)**2)*PN
      IF((JZ.GE.36).AND.(JZ.LE.43))
     1  QH = 10.0*EXP(- 55.*(AN/ZJ - 1.33)**2)*PN
      IF((JZ.GE.44).AND.(JZ.LE.50))
     1  QH = 10.0*EXP(- 70.*(AN/ZJ - 1.36)**2)*PN
      IF (JZ.LE.50)    GO TO 130
      EE = EI - 100
      IF (EE .LT. 150.)  EE = 150.
      HE = .03*(ALOG(EE))**2
      IF (HE.GT. 1.)   HE = 1.
      ZZ = JZ
      IF (JZ .LT. 55)  ZZ = 55.
      HS = (1.5*EXP(-1.0/ZJ*(ZZ-55.)**2) + 3.0) * FU * HE
      HF =  5.0*EXP(-10./ZJ*(ZZ-55.)**1.5)
      CON = 80. + 220.*(1. - EXP(-(ZJ-55.)**2/4.))*(JZ/55)
      IF((JZ.GE.51).AND.(JA.LE.155))
     1  QH = HF*EXP(-CON *(AN/ZJ - 1.48)**2)*PN
     1     + HS*EXP(-260.*(AN/ZJ - 1.29)**2)*PN
      IF((JZ.GE.66).AND.(JZ.LE.87))   GO TO 140
      IF (IZ.LT.88)  GOTO 140
      QJ = QH
      IF (IZ .LE. 90)   QH = 0.6*QH
 130  QF = QH*FF
      QJ = QF
      GO TO 200
C
 140  IF (IZ.GE.29 .AND. JZ.EQ. 5)  FA = EXP(.02*(AI-56)*(ANZC-0.6))
      IF((IZ.GE.29).AND.(JZ.GE. 6).AND.(JZ.LE. 8))
     1    FA = EXP (.020*(AI - 56.)*(ANZC-0.7))
      IF((IZ.GE.29).AND.(JZ.GE. 9)             )
     1    FA = EXP (.020*(AI - 56.)*(ANZC-0.9))
      IF (JZ.LE.11)      FA = FA*EXP(-2.5*(ANZC-1.))            
      IF (FA.LT.1.)      FA = 1.0
      IF (IA.GE.110)     FM = FF
      IF (IA.LT.110)     FM = FE
      IF((IA.GE.110).AND.(JA.LT.57).AND.(AJ.GT.0.23*AI)) FM=AMAX1(FE,FF)
      IF((IA.GE.110).AND.(JA.LT.57).AND.(AJ.LE.0.23*AI)) FM=FE
      DMAX = 31.5 + .045*(AI - 36.0)*(ALOG(AI) + 1.23)
      AC = DMAX
      AE = AC
      AH = DA
      IF (AH .GT. DMAX)   AH = DMAX
      IF (EI .LT. EC)  DMAX = 31.5 + D1*(AI-36)*(ALOG(EI) - D2) 
      AE = DMAX
      IF (DA.GT.DMAX)    DA = DMAX
 150  F1 = EXP(-.25 + .0074*AI)
      IF (EI.GT.750.)    F1 = 1.00 + (F1 - 1.00)*((EC-EI)/(EC-750.))**2
      IF (EI.GE. EC )    F1 = 1.00
      F2 = EXP(1.73 - .0071*EI)
      IF (F2 .LT. 1.)  F2 = 1.
      PE = 20./EI**.77
      PH = P0 /AI**P1
      IF((IZ.GE.20).AND.(IZ.LE.30))
     1PE = PE*(1. - .32*EXP (-((EI - 100.)/100.)**2))
      IF (IA.GT.100)
     1PE = PE*(1. - .000015*(AI-100.)*(EC+150.)/(EI+150.))
      IF (IA.LE.71 .AND. EJ.GE.  EC )  PE = PH
      IF (IA.GT.72 .AND. EJ.Ge.3000.)  PE = PH                          
      IF (EJ.GE.EC .AND. EJ.GE.3000.)  PE = PH
      PA = PE * AI
      HA = PH * AI
      RJ = R0*AJ**(-R1)*(1.-0.4*(IZ/88))
      IF (AJ .LT. 40.)  RJ = 1.29*AJ**.15*(1.-0.4*(IZ/88))
      IF((AA.LT.AC).AND.(JZ.EQ.40))  DM =-1.1
      IF((AA.LT.AC).AND.(JZ.EQ.42))  DM = 0.8
      IF((AA.LT.AC).AND.(JZ.EQ.53))  DM =-0.5           
      IF((AA.LT.AC).AND.(JZ.EQ.55))  DM = 1.7
      IF((AA.LT.AC).AND.(JZ.EQ.60))  DM =-1.1
      IF((AA.LT.AC).AND.(JZ.EQ.61))  DM =-1.3
      IF((AA.LT.AC).AND.(JZ.EQ.63))  DM = 0.9
      IF((AA.LT.AC).AND.(JZ.EQ.64))  DM = 0.7
      IF((IZ.GT.76).AND.(JZ.GT.62))  DM = DM + 1.
      AJ = AJ + DM
      S  = S0 - S1*(AI - AM)/ZI
      ST = ZJ - (S - T2*AJ - T3*AJ*AJ)*AJ
      IF (ST .LT.-1.)  ZA = ABS(ST)**1.30               
      IF (ST .GE.-1.)  ZA = ABS(ST)**1.50                       
      IF (ST .GT. 1.)  ZA = ABS(ST)**1.75                       
      Q0 = C1*F1*F2*PE*AI**C2/(1.-C3/PA-(C4-C3/PA)*EXP(-PA))    
      H0 = C1*1.00 *PH*AI**C2/(1.-C3/HA-(C4-C3/HA)*EXP(-HA))    
      QR = CJ*PN*EXP(-PE*DA-RJ*ZA) * Q0
      QH = CJ*PN*EXP(-PH*AH-RJ*ZA) * H0
      AX = 0.5*(S - SQRT(S*S - 4.*T2*ZJ))/T2 - DM               
      KX = (AI - AX) - (JP-1)
      XP = KX
      IF (IA.GT.70.AND.KN.GT.KM.AND.KN.LT.KX.AND.KFLAG.EQ.1) GOTO 155
      KFLAG=0
      IF (EJ.GT.EC)   QR = QH
      IF (EI.GE.2500. .AND. QR.GT.QH)  QR=QH                            
      F500=EXP(.0074*AI-.25)
      Q500 =CJ*PN*C1*F500 *AI**C2*P500 *EXP(-P500 *DA-RJ*ZA)
     &     /(1.-C3/P500 /AI-(C4-C3/P500 /AI)*EXP(-P500 *AI))
      F1000= 1.+(EXP(.0074*AI-.25)-1.)*((EC-1.E3)/(EC-750.))**2
      Q1000=CJ*PN*C1*F1000*AI**C2*P1000*EXP(-P1000*DA-RJ*ZA)
     &     /(1.-C3/P1000/AI-(C4-C3/P1000/AI)*EXP(-P1000*AI))
      F2000= 1.+(EXP(.0074*AI-.25)-1.)*((EC-2.E3)/(EC-750.))**2         
      Q2000=CJ*PN*C1*F2000*AI**C2*P2000*EXP(-P2000*DA-RJ*ZA)            
     &     /(1.-C3/P2000/AI-(C4-C3/P2000/AI)*EXP(-P2000*AI))            
      F3000= 1.+(EXP(.0074*AI-.25)-1.)*((EC-3.E3)/(EC-750.))**2
      Q3000=CJ*PN*C1*F3000*AI**C2*P3000*EXP(-P3000*DA-RJ*ZA)
     &     /(1.-C3/P3000/AI-(C4-C3/P3000/AI)*EXP(-P3000*AI))
      FDAE = AMIN1((AI-AJ)/(.14*AI)*(EI/1000.)**(-2./3.),2.)
        FDAE1 = (AI-AJ)/(.14*AI)                                        
      IF (IZ.GE.29.AND.IZ.LE.83.AND.EI.LE.1.E3)   QR = QR*FDAE
      IF (IZ.GE.29.AND.IZ.LE.83.AND.EI.GE.1.E3.AND.EI.LE.3.E3)          
     &  QR = Q1000*FDAE1 + (QH-Q1000*FDAE1)*((EI-1000)/2000.)           
      QJ = QR
      QI = QR
      IF (IZ.LT.88 .OR. JP.GT.3)  GOTO 160
      QJ = QR*FU
      GO TO 200
 155  F1000=1.00+(EXP(.0074*AI-.25)-1.00)*((EC-1000.)/(EC-750.))**2
      Q1000=CJ*PN*C1*F1000*AI**C2*P1000*EXP(-P1000*DA-RJ*ZA)
     1     /(1.-C3/P1000/AI-(C4-C3/P1000/AI)*EXP(-P1000*AI))
      X1000=Q1000*EXP(-P1000*(AI-AX))/EXP(-P1000*DA-RJ*ZA)
      HR=QH*EXP(-PH*(AI-AX))/EXP(-PH*AH-RJ*ZA)
      IF (X1000.GT.HR)  X1000=HR
      F500 =EXP(.0074*AI-.25)
      Q500 =CJ*PN*C1*F500*AI**C2*P500*EXP(-P500*DA-RJ*ZA)
     1     /(1.-C3/P500/AI-(C4-C3/P500/AI)*EXP(-P500*AI))
      X500=Q500*EXP(-P500*(AI-AX))/EXP(-P500*DA-RJ*ZA)
      QR=QR*EXP(-PE*(AI-AX))/EXP(-PE*DA-RJ*ZA)
      IF (EI.GT.EC)                QR=HR
      IF (EI.GE.1000. .AND. QR.GT.HR) QR=HR
      IF (EI.LT.1000. .AND. EI.GT.500.)
     1QR=X500+(X1000-X500)*(EI-500.)/500.
      QH=HM+(HR-HM)*(XN-XM)/(XP-XM)
      QJ=QM+(QR-QM)*(XN-XM)/(XP-XM)
      KFLAG=0
      GO TO 200
 160  DD = AC - AA
      XZ = 1.29 + 0.005*DD + CN
      IF (QH .GT. 0.)  RH = QR/QH
      IF (GF .GT. 0.)  FR = FE/GF
      IF((IA.GE. 69).AND.(RH.LT.FR).AND.(EI.LT.300.))  QR = QH*FR
      IF (IZ.GE.76.AND.IZ.LE.80) FF=FF*AMIN1(1050./EI+EI/EC,6.)
      IF((IA.GE.110).AND.(IA.LE.238))   QF = QH*FF*FA
      IF((IA.GE. 69).AND.(JA.LE. 57))   QE = QH*FE*FA
      IF (IA.LT.110. OR. IA.GT.209 .OR. JA.LE.56)  GO TO 190
      IF (DD.GE. 20)                               GO TO 170
      IF (IA.LE.125)                               GO TO 180
      IF (AZ.LE.XZ .OR. JZ.GT.57)                  GO TO 170
      KZ=51+IZ/76+IZ/80
      IF (AZ.GT.XZ .AND.JZ.LT.KZ)                  GO TO 180
      ER = EI**(2./3.)
      IF (ANZC.LT.1.56)  GA = (.03 + (ZJ-KZ)*.007)*(1.56-ANZC)*ER
      IF (ANZC.GE.1.56)  GA = 0
      IF (GA  .GT.  1.)  GA = 1.
      IF((JZ.GE. KZ).AND.(JZ.LE.57).AND.(DD.LT.20.).AND.(AZ.GT.XZ))
     1    QJ = QR**GA*QF**(1.-GA)
      GO TO 200
 170  QJ = QR
      GO TO 200
 180  QJ = AMAX1(QR, QF)
      GO TO 200
 190  IF((IA.GE.210).AND.(JA.GE.57))                         QJ=QF
      IF((IA.GE.110).AND.(JA.LE.56).AND.(AJ.GT.0.23*AI))     QJ=QH*FA*FM
      IF((IA.GE.110).AND.(JA.LE.56).AND.(AJ.LE.0.23*AI))     QJ=QE
      IF((IA.LT.110).AND.(IA.GE.69).AND.(DD.GE. 0 ))         QJ = QR
      IF((IA.LT.110).AND.(IA.GE.69).AND.(DD.LT. 0 ))         QJ = QE
      IF (IA.LT. 69)                                         QJ = QR
      IF (QJ.EQ.QR .AND. QI.GT.0.)   QH = QH*QR/QI + .000001
      IF((IZ.GE.90).AND.(JZ.GE.66))
     1QJ = HF*EXP(-CON*(AN/JZ - 1.48)**2)*PN*(0.6+0.2*(IZ-90))
     1      + HS*EXP(- RJ*ZA )*PN*(0.6+0.2*(IZ-90))
      IF((QJ.EQ.QE).OR. (QJ.EQ.QF))   GA = 0.
      IF (QJ.EQ.QR)                   GA = 1.
      IF (JP.LT.3 .OR. JP.GT.5)   GO TO 200
      FE = 1. - EXP(-(EI/35.)**4)
      IF (EI .LT. 200.)   FE = FE*EI/200.
      IF (JP.GE.3 .AND. KN.GE.1 .AND. IA.LE.70)   GO TO 200
        IF (IZ.GT.30.AND.JP.GE.3.AND.KN.GE.KM)  GO TO 200               
      QJ = QH*FE*YA
 200  EX = E0*56.**E1*EJ/EC
      G1 = 1.-0.6*(1.-EXP(-(EX/1000.)**2))*EXP(-(EX/2000.)**2)
     1       +0.2*(1.-EXP(-(EX/3000.)**2))
      G2 = 1.-0.4*(1.-EXP(-(EX/2000.)**2))*EXP(-((EX-1800.)/1800.)**2)
     1       +.17*(1.-EXP(-(EX/2000.)**2))
      G3 = 1.+.25*(1.-EXP(-(EX/1500.)**2))*EXP(-((EX-1500.)/1800.)**2)
     1       -.05*(1.-EXP(-(EX/2000.)**2))
      G4 = 1.-0.1*(1.-EXP(-(EX/4000.)**2))
      IF (IZ.LE.28 .AND. JZ.LE.4 .AND. JA.LE.12)  QJ=QJ*G1
      IF (IZ.LE.28 .AND. AA.GT.AC)                QJ=QJ*G2
      IF (IZ.GT.28 .AND. JA.LE.56 .AND. AA.GT.AC) QJ=QJ*G2
      IF (AA.GE.7. .AND. AA.LE.AC-13.)            QJ=QJ*G3
      IF (IZ.GE.90 .AND. JA.GT.56 .AND. AA.GE.7.) QJ=QJ*G4
      MX=10
      IF (IA.LT.150.AND. JP.LT. 3)  MX=9
      IF (IA.GE.150.AND. JP.EQ. 3)  MX=11
      IF (JP.EQ. 1 .AND. KN.GT.MX)  QJ=QJ*(.1+.9*EXP(-(XN-MX)**2/4.))
      IF (JP.EQ. 2 .AND. KN.GE.MX)  QJ=QJ*(.1+.9*EXP(-(XN-MX)**2/4.))
      IF (JP.EQ. 3 .AND. KN.GE.MX)  QJ=QJ*(.5+.5*EXP(-(XN-MX)**2/4.))
      QJ = QJ*FP                                                        
      END
      SUBROUTINE PXN(Z, A, X, E, QJ)
      COMMON /A/ Q1, Q2, QI, F1, F2, F3
      EI = 500. + 300.*X
      Q1 = 3.3*EXP(-ABS(6.9-X)**2.8/67.45-90*X**2.35*ABS(2.5-A/Z)**5)
      QA = 3.3*EXP(-ABS(6.9-X)**2.8/67.45 - .04638*X**2.35)
      QB = 3.3*EXP(-ABS(6.9-X)**2.8/67.45)
      IF (A/Z.GE.2.28 .AND. A/Z.LE.2.50)
     1Q1 = QA**((2.5-A/Z)/0.22) * QB**((A/Z-2.28)/0.22)
      IF (A/Z .GT. 2.5)   Q1 = QB
      QI = Q1
      XZ = ((A - Z) - X + 1)/Z
      Q2 = 0.5*EXP(-90.*ABS(1.5-XZ)**5)
      IF (XZ .LE. 1.2)   QI = AMIN1(Q1, Q2)
      FI = 1.0                                          !### imos 1/21/2000
      IF (Z  .GE. 81.)   FI = 1.5*EXP(-X*((Z-80.)/12.)**5)
      IF (FI .GT. 1.0)   FI = 1.
      IF (Z  .GE. 81.)   QI = QI*FI
      IZ = Z +.1
      NX = X +.1
      IA = A +.1
      NT = IA-IZ
      N  = NT-NX
      ED = E
      IF (NX.GE.3 .AND. IZ.GE.39 .AND. NT.GE.28 .AND. N .LT.50)
     1    ED = E - 10.
      IF (NX.GE.3 .AND. IZ.GE.59 .AND. NT.GE.50 .AND. N .LT.82)
     1    ED = E - 3.
      B  =  A
      IF (B .LT. 35.)   B = 35.
      FX = (12.+.1*X)*X - (1.-1./X**.5)*B**(2./3.)
      ALOGY=(1.+1.5*X)*ALOG(ABS(ED/FX))
      IF (ALOGY .GT. 10.)  ALOGY=10.
      F1 = 1. - EXP(-EXP(ALOGY))
      C  = 1. + 0.03*X*(A-208.)
      IF (C .LE. 1.)   C = 1.
      D  =  27.5 - 0.1*(A + (200.-A)/X**.5)
      IF ((1.-X)*(208.-A) .GT. 0)   D = D - 0.03*(1.-X)*(208.-A)
      F  = E
      IF (NX.GE.3 .AND. IZ.GE.39 .AND. NT.GE.28 .AND. N .LT.50)
     1    F = E - 3.
      IF (NX.GE.3 .AND. IZ.GE.59 .AND. NT.GE.50 .AND. N .LT.82)
     1    F = E - 3.
      IF (NX.EQ.1 .AND. IZ.GE.79 .AND. IZ.LE.83)  F = E + 5.
      F2 = 3500/C * EXP(-.6*X-.5*(((F-D)-5.*X**1.34)/(6.-2.5/X**4))**2)
     1             * (1. - EXP(-(0.03*A/(2.-1./X**4))**3))
      G  =  0.01*(A-208.)
      IF (A .LT. 208.)   G = 0.
      F3 = (1300./(E + 20.*X**1.5/E))**(1.3 - G)
      FH = (1300./(EI+20.*X**1.5/EI))**(1.3-G)
      FE = F1*(F2 + F3)
      IF (E .GE. EI)   QJ = FH*QI
      IF (E .LT. EI)   QJ = FE*QI
      IF (Z.EQ.6. .AND. A-X+1..EQ.13.)   QJ = QJ*0.4
      IF (Z.EQ.3. .AND. A-X+1..EQ. 9.)   QJ = QJ*0.65
      END
