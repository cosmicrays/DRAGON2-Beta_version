c!     **.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
c!      * crn6.f *                                      galprop package * 2001/05/11
c!     **"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
c!     
c!     **.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
c*     Joint Institute for Nuclear Research
c*     Laboratory of Computing Techniques and Automation
c*     V.S.Barashenkov,A.Polanski
c*     
c*     Electronic Guide
c*     for Nuclear Cross Sections
c*     
c*     A short description of a fortran CROSEC code providing the integral cross
cc*     -sections for pion-nucleus,nucleon-nucleus and nucleus-nucleus interactions
c*     (total,nonelastic,elastic).The hadron-nucleus cross-sections are obtained,by
c*     means of interpolation between evaluated experimental data at target mass 
c*     numbers A>=4 and energies from 14(20)MeV up to 1 TeV.The nucleus-nucleus cross
c*     sections are calculated with the help of approximation formula with fitted
c*     coefficient at energies above several Mev/nucleon.The CROSEC code uses 168 KB
c*     memory and can be used in interactive way to generate separate values or ta-
c*     bles of cross sections on a display screen or file.Fractions of the CROSEC
c*code can be used as subroutines employed by other codes.A set of all available
c*experimental integral hadron-nucleus cross-sections at energies exceeding 14
c*MeV and plots of evaluated total,nonelastic and elastic cross sections for pio
c*and nucleon interactions are presented in/3-5/.Two methods have been employed
c*to calculate dependence of cross-sections vs energy.At high energies,where
c*the projectile de Brogle wave length is significantly smaller then the size
c*of the target nucleus quasioptical model is used.The parameters of the
c*model have been fitted to obtain best agreement of calculated and experimental
c*data.The high-energy region has been divided into separate intervals with
c*characteristic behavior of cross sections.(For example,the region near the
c*minimum of nucleon cross-sections at energy about 200 MeV - the resonance
c*region.In case of pion-nucleus cross sections,the interval of smooth cross
c*section alterations at energy above 1 GeV).A set of parameters has been
c*defined for each interval.Phenomenological approximation of cross sections
c*was used at lower energies/3-5/.
c*  The known experimental information on nucleus-nucleus cross-sections is
c*insufficient to compile a detailed plots curves,especially if one considers
c*the great number of particular interesting nuclei pairs.In this case approxi-
c*mation relations can be used with coefficients fitted by means of comparison
c*with the known experimental data/5/.Of course,the accuracy of the results is
c*lower then that for hadron-nucleus interactions.
c*   The CROSEC code includes following modules:
c*c* main module-reading input and printing separate values or tables of hadron
c*  -nucleus and nucleus-nucleus cross sections(mbarns)
c*c* FUNCTION SIGHAD-calculation of pion-nucleus and nucleon-nucleus cross
c*  -sections.
c*c* FUNCTION SIGION-calculation of nucleus-nucleus cross-sections
c*c* FUNCTION FHS-   calculation of high-energy nucleus-nucleus cross
c*  -sections
c*c* FUNCTION FC- calculation of low-energy parameters
c*c* BLOCK DATA A-target nucleus mass numbers
c*c* BLOCK DATA B-projectile kinetic energies
c*c* BLOCK DATA C-nucleus-nucleus cross section parameters
c*c* BLOCK DATA S1-neutron-nucleus cross-sections
c*c* BLOCK DATA S2-proton-nucleus cross-sections
c*c* BLOCK DATA S3-pi-meson-nucleus cross sections
c*c* BLOCK DATA S4-pi+meson-nucleus cross-sections
c*c* Integer BLOCKS B,S1,S2,S3,S4 there are on the BARPOL.FOR or on the BARPOL.DAT 
c*After starting the program you must insert the eight numbers in free format:
c*ITYPE
c*PA_PZ_TA_TZ_E1_E2_ES
c*representing respectively:
c*ITYPE=1 OR ITYPE=2 - Calculation of total or nonelastic cross-sections,
c*ITYPE=3            - Calculation at once of total,nonelastic and elastic cross
c*                     sections,
c* ITYPE=-1            -Read cross sections from BARPOL.DAT
c* ITYPE=0             -Exit
c*PA,PZ-projectile mass and charge numbers(for pions PA<0.2)
c*TA,TZ-the same for target nucleus,
c*E1,E2-the lowest and the highest energies in the considered interval,
c*ES   -energy step,
c*   For example if one needs to obtain cross sections for pion(pi-)+nucleus
c*Pb-207.19 in the energy range from 20 MeV to 200 MeV with a step 10 Mev one
c*must enter(in free format):
c*3 
c*0.15 -1 207.19 82 20 200 10
c*
c*The SIGHAD and SIGION functions can be used separately with others codes
c*to calculate a current value of hadron-nucleus or nucleus-nucleus cross
c*-section:
c*CS=SIGHAD(IS,PA,PZ,TA,TZ,E)
c*CS=SIGION(IS,PA,PZ,TA,TZ,E)
c*where:E-projectile energy in MeV (for hadrons) or in MeV/nucleon(for nucleus)
c*IS=1 or IS=2 -calculation of total or nonelastic cross-sections,
c*PA,PZ,TA,TZ-explained above.
c*   The CROSEC code is written in FORTRAN.
c*   Further information can be requested at e-mail addresses:
c*   barashen@lcta30.jinr.dubna.su
c*   polanski@cyf.gov.pl
c*   olek@neutron.kth.se
c*References:
c*1.Barashenkov V.S.,Polanski A.,Comm. JINR E2-94-417,Dubna, 1994.
c*2.Barashenkov V.S.Interaction cross-sections of particles and nuclei
c*  JINR.Dubna,1993.(In Russian)
c*3.Barashenkov V.S.Gareeva G.F.,Polanski A.,Comm. JINR 10-92-214,Dubna 1992.(In Russ.)
c*4.Barashenkov V.S.,Polanski A.,Comm. JINR B2-90-489,Dubna, 1990.(In Russ.)
c*5.V.S.Barashenkov, A. Polanski, A.N. Sosnin. Comm. JINR P2- 90-159,Dubna 1990.(In Russ.)
c*6.V.S.Barashenkov, A. Polanski, A.N. Sosnin.Comm. JINR P2-89-753, Dubna 1989.(In Russ.)
!
C     PROGRAM CROSEC(INPUT,OUTPUT)
      subroutine CROSEC

      IMPLICIT REAL*8 (A-H,O-Z)
 5000 FORMAT(' ****************************************************' 
     *  /,   ' CODE FOR CALCULATION OF NUCLEON-NUCLEUS,PION-NUCLEUS'
     *  /,   ' AND NUCLEUS-NUCLEUS TOTAL,NONELASTIC AND ELASTIC    '
     *  /,   ' ++++++++++++++++ CROSS-SECTIONS(MBARNS)+++++++++++++'
     *  /,   ' ****************************************************'
     *  /,   ' WRITTEN BY VLADILEN S.BARASHENKOV-JINR-DUBNA(RUSSIA)' 
     *  /,   ' AND ALEKSANDER POLANSKI-SINS-SWIERK(POLAND)         '
     *  /,   ' ****************************************************'
     *  /,   ' *********** e mail:polanski@ipj.gov.pl ************')

C     Input parameters:
* ITYPE
* PA_PZ_TA_TZ_E1_E2_ES
* representing respectively:
* ITYPE=1 OR ITYPE=2 - Calculation of total or inelastic cross-sections,
* ITYPE=3            - Calculation at once of total,inelastic and elastic cross
*                      sections,
* ITYPE=-1            -Read cross sections from BARPOL.DAT
* ITYPE=0             -Exit
* PA,PZ-projectile mass and charge numbers(for pions PA<0.2)
* TA,TZ-the same for target nucleus,
* E1,E2-the lowest and the highest energies in the considered interval,
* ES   -energy step,
* To break the procedure one must insert ITYPE=0 or PA=0 
*    For example if one needs to obtain cross sections for proton+nucleus
* Pb-207.19 in the energy range from 20 MeV to 40 MeV with a step 2.5 MeV one
* must enter(in free format):
* 3 
* 1 1 207.19 82 20 40 2.5 
C     REGION OF APPLICABILITY OF THIS CODE:
C     FROM 14 MEV UP TO 1 TEV FOR NUCLEON-NUCLEUS COLISIONS
C     FROM 20 MEV UP TO 1 TEV FOR PION-NUCLEUS COLISIONS
C     FROM 1.0 MEV/NUCLEON  UP TO 1 TEV/NUCLEON FOR NUCLEUS-NUCLEUS
C     COLISIONS 
C     A1,A2<240
c     A2>=4
 6000 FORMAT(' ENTER:'/' ITYPE (ITYPE=0 exit,ITYPE>0 continue)')
     
 6001 FORMAT(' ENTER:'/' PA_PZ_TA_TZ_E1_E2_ES'/)
      ITYPE=-1
      CALL SIGTAP2(ITYPE)
      PRINT 5000
 100  PRINT 6000
      IN=5
      READ(IN,*) ITYPE
      IF(ITYPE.LT.0) CALL SIGTAP2(ITYPE)
      IF(ITYPE.EQ.0) GO TO 200
      IF(ITYPE.LT.0) GO TO 100
      PRINT 6001
      READ(IN,*) PA,PZ,TA,TZ,E1,E2,ES
c      WRITE(6,*) PA,PZ,TA,TZ,E1,E2,ES
      NE1=0
      ND=1
      ISS=ITYPE
      IF(ES.GT.0.0.AND.E2.GT.0.0) NE1=(E2-E1)/ES
      NE=1+ABS(NE1)
      IF(E1.GT.E2) ES=-ES
      IF(ITYPE.GT.2) ND=2
      IF(PA.EQ.1.) WRITE(6,1000)
      IF(PA.LT.1.) WRITE(6,2000)
      IF(PA.GT.1.) WRITE(6,3000)
      DO 10 IE=1,NE
      CR1=0.
      CR2=0.
      CR3=0.
      DO 20 I=1,ND
      IF(ITYPE.GT.2) ISS=I
      E=E1+(IE-1)*ES
      IF(PA-1.0) 11,11,12
   11 CS=SIGHAD(ISS,PA,PZ,TA,TZ,E)
      GO TO 13
   12 CS=SIGION(ISS,PA,PZ,TA,TZ,E)
   13 CONTINUE
      IF(ISS.EQ.1) CR1=CS
      IF(ISS.EQ.2) CR2=CS
   20 CONTINUE
      IF(ITYPE.GT.2) CR3=CR1-CR2
      WRITE(6,7000) E,CR1,CR2,CR3
   10 CONTINUE
      GO TO 100 
 1000 FORMAT('     NUCLEON NUCLEUS CROSS-SECTIONS (MBARNS) '/
     *        ' ENERGY(MEV)  TOTAL   NONELASTIC     ELASTIC ')
 2000 FORMAT('         PION-NUCLEUS CROSS-SECTIONS(MBARNS) '/
     *       ' ENERGY(MEV)   TOTAL   NONELASTIC     ELASTIC ')
 3000 FORMAT('         NUCLEUS-NUCLEUS CROSS-SECTIONS(MBARNS)'/
     *       ' ENERGY(MEV/NUC) TOTAL   NONELASTIC     ELASTIC ')
 7000 FORMAT(4(F8.1,4X))
 200  CONTINUE
      END


      REAL*8 FUNCTION BINT(U,E,F,N,IS)
      IMPLICIT REAL*8 (A-H,O-Z)
C     LINEAR INTERPOLATION  IS=1
C     QUADRATIC INTERPOLATION  IS=2
      DIMENSION E(N),F(N)
      IF(IS.LE.0.or.n.eq.1) BINT=U
      IF(IS.LE.0.or.n.eq.1) RETURN
      IF(N.GT.2) GO TO 10
      X1=E(1)
      Y1=F(1)
      X2=E(2)
      Y2=F(2)
      GO TO 8
 10   CONTINUE     
      IF(U-E(1))1,1,2
  1   X1=E(1)
      Y1=F(1)
      X2=E(2)
      Y2=F(2)
      X3=E(3)      
      Y3=F(3)  
      GO TO 7
  2   IF(U-E(N-1)) 3,4,4
 4    X2=E(N-1)
      X3=E(N)      
      Y2=F(N-1)
      Y3=F(N)
      X1=E(N-2)
      Y1=F(N-2)  
      GO TO 7
 3    CONTINUE
      IF(N.LE.2) GO TO 7
      N1=N-1
      DO 5 J=2,N1
      IF(U-E(J)) 6,5,5
  6   X1=E(J-1)
      X2=E(J)
      X3=E(J+1)
      Y1=F(J-1)
      Y2=F(J)
      Y3= F(J+1)
      GO TO 7
 5    CONTINUE
 7    CONTINUE
      IF(IS.NE.2.OR.N.EQ.2)  GOTO 8
       BINT=Y1*(((U-X2)*(U-X3))/((X1-X2)*(X1-X3)))+
     *      Y2*(((U-X1)*(U-X3))/((X2-X1)*(X2-X3)))+
     *      Y3*(((U-X1)*(U-X2))/((X3-X1)*(X3-X2)))
  8   CONTINUE
      IF(IS.EQ.1.OR.N.EQ.2) BINT=Y1+(U-X1)*(Y2-Y1)/(X2-X1)
      RETURN
      END


      REAL*8 FUNCTION SIGION(ISS,A1,Z1,A2,Z2,T)
      IMPLICIT REAL*8 (A-H,O-Z)
C     FOR CALCULATION OF NUCLEUS-NUCLEUS TOTAL (ISS=1)
C     AND INELASTIC (ISS=2) CROSS SECTIONS
C     A1,Z1 - PROJECTILE MASS AND CHARGE NUMBERS (A1>1)
C     A2,Z2 - THE SAME FOR TARGET NUCLEUS (3<A2<240)
C     T - LAB. KINETIC ENERGY OF PROGECTALE (1 MEV/NUCLEON< T <1 TEV/NUCLEON)

      COMMON/CX/CX(38)
      COMMON /FH/AMP,AMT,AP,AT,B0,R0
C      WRITE(*,*)' SIGION',ISS,A1,Z1,A2,Z2,T
      IS=3-ISS
      IF(A1.LT.1.0D0.OR.A1.GT.240.0D0.OR.A2.LT.3.0D0.OR.A2.GT.240.0D0)
     *GO TO 101
      IF(DABS(Z1).LT.1.0D0) GO TO 101
      IF(T.LT.1.0D0) GO TO 101
      SIGION=0.0D0
      TP=T/A1
      AP=A1**0.333333
      AT=A2**0.333333
      AMP=A1*930.63D0
      AMT=A2*930.63D0
C     PARAMETER FOR CALCULATION OF NUCLEAR RADIUS
      R0=1.4D0
      IF(DABS(A1-4.0D0) .LT. 0.1D0) R0=1.3D0
      B0=1.44D0*Z1*Z2
      I=1
       IF(IS.EQ.2) I=20
C     SELECTION OF PROGECTALES
C     HEVY ION
      IF(A1.GT.4.1D0) N=I
C     ALFA,HELION,TRITON
      IF(A1.GT.2.1D0 .AND. A1.LT.4.1D0) N=I+6
C     DEUTRON
      IF(A1.LT.2.1D0) N=I+12
C     HIGH-ENERGY CROSS-SECTION
C     SELECTION OF PROGECTALE ENERGY
      IF(TP.LT.CX(N+1)) K=2
      IF(TP.LT.CX(N+4)) K=5
C     CROSS-SACTION PARAMRTERS
      C=CX(I)
      IF(TP.LT.CX(N+1)) C=CX(N+K)+CX(N+K+1)*DLOG10(TP)
      CP=CX(N+5)+CX(N+6)
      IF(TP.LT.10.) GO TO 1
C     HIGH-ENERGY CROSS-SECTION
      SIGION=FHS(IS,T,C)
      RETURN
C     CALCULATION OF LOW-ENERGY CROSS-SECTION
C     NORMALUSED HIGH-ENERGY CROSS-SECTION
 1     SH10=FHS(IS,10.*A1,CP)
      R0=1.45D0
      IF(DABS(A1-4.0D0).LT. 0.1D0) R0=1.4D0
C     RENORMALUSED COULOMB BARRIER
      B=B0/R0/(AP+AT)
C     LOW-ENERGY CROSS-SECTION
      SIGION=SH10*FC(T,B)/FC(10.*A1,B)
      IF(SIGION) 101,100,100
  101 CONTINUE
C      WRITE(*,1001)
C      WRITE(*,*)' ISS',ISS,' A1',A1,' Z1',Z1,' A2',A2,' Z2',Z2,' T',T
C      PAUSE
      SIGION=1.0D-07
C 1001 FORMAT(' ERROR IN INPUT OF PARAMETERS OF FUNCTION SIGION')
  100 CONTINUE
      RETURN
      END
      REAL*8 FUNCTION FC(T,B)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /FH/AMP,AMT,AP,AT,B0,R0
C     CMS ENERGY
      TC=T*AMT/(AMP+AMT)
      X=(TC-B)/1.2
      IF(X.GT.5) GO TO 1
      D=1.+EXP(X)
      FC=DLOG(D)/TC
      RETURN
  1   FC=X/TC
      RETURN
      END


      REAL*8 FUNCTION FHS(IS,E,C)
      IMPLICIT REAL*8 (A-H,O-Z)
C     CALCULATION OF HIGH-ENERGY TOTAL (IS=2) AND
C     INELASTIC (IS=1) CROSS-SECTIONS

C      E - LAB. KINETIC ENERGY OF PROGECTALE(MEV)
      COMMON/FH/AMP,AMT,AP,AT,B0,R0
C     SQUED PROGECTALE CMS MOMENTUM
      PPC=AMT*AMT*E*(E+2.*AMP)/((AMP+AMT)**2+2.*AMT*E)
C     DE BROGLE WAVE LANGTH
      AL=1.41*140./SQRT(PPC)
      EC=SQRT(PPC+AMP*AMP)-AMP
C     COULOMB BARRIER
      B=B0/R0/(AP+AT+AL)
      FHS=31.416*1.21*(1.-B/EC)*(AP+AT+1.85*AP*AT/(AP+AT)
     *+AL-C)**2*IS
      RETURN
      END


      BLOCK DATA C
      IMPLICIT REAL*8 (A-H,O-Z)
C     NUCLEUS-NUCLEUS CROSS-SECTION PARAMETERS

      COMMON/CX/CX(38)
      DATA  CX/2.07,560.,0.8,0.426,
     *           100.,-2.05,1.9,
     *           200.,0.07,0.87,
     *           20.,-1.55,2.1,
     *           700.,-1.01,1.08,
     *           400.,-0.59,0.94,
     *        2.45,225.,-2.25,2.,
     *           100.,-4.61,3.18,
     *           185.,-3.,2.4,
     *           185.,-3.,2.4,
     *           185.,-4.77,3.18,
     *           185.,-4.77,3.18/
      END


      REAL*8 FUNCTION SIGHAD(IS,A1,Z1,A2,Z2,T)
      IMPLICIT REAL*8 (A-H,O-Z)
C     CODE FOR CALCULATION OF NUCLEON-NUCLEUS AND PION-NUCLEUS
C     TOTAL(IS=1) AND NONELASTIC(IS=2)CROSS-SECTIONS(MBARNS) 
C     IS=0 EXIT
C     A1=1.0 FOR NUCLEON OR 0<A1<0.2 FOR PIONS
C     Z1-PROJECTILE CHARGE NUMBER
C     A2,Z2- TARGET NUCLEUS MASS AND CHARGE NUMBERS (4.0<=A2<=239.0)
C     T- PROJECTILE PARTICLE KINETIC ENERGY(MEV;14(20)MEV<T<1TEV)
C     AAI1(K),AAI2(K)- TARGET MASS NUMBERS for nucleons and pions
C     IENER(L,JS)-PROJECTILE ENERGY(MEV)
C     CROSS-SECTIONS ARE STORED IN ISIG1(K,L),ISIG2(K,L),ISIG3(K,L),ISIG4(K,L) 
C      ISIG(JS,K,L)- JS=1 NEUTRON CROSS SECTIONS  
C                    JS=2 PROTON CROSS SECTIONS 
C                    JS=3PI- MESON CROSS SECTIONS  
C                    JS=4PI+ MESON CROSS SECTIONS 
C      K- VS. TARGET MASS NUMBERS
C      L- VS. ENERGY
      COMMON /BARPO1/ AAI1(24),AAI2(20)
      COMMON/BARPO/ NEL(2),NE(4),IENER(62,4),ISIG(4,24,62)
      DIMENSION SS(24),EE(3),FF(3),SIG(2),JSIG(2)
      SIGHAD=0.0
      IF(IS.LE.0) RETURN
      IF(A1.LE.0.OR.A1.GT.1.0) GO TO 101
      IF(A1.LT.1.0.AND.A1.GT.0.2) GO TO 101  
      IF(A2.LT.1.0.OR.A2.GT.250.) GO TO 103 
      IF(ABS(Z1).GT.1.0) GO TO 101
      IF(Z1.EQ.0.0.and.T.LT.10.0) GO TO 102
      IF(T.LT.1.0) go to 102
      NS=4
C     NE(JS) - NUMBER OF ENERGY POINTS (62 for protons,53 for neutrons, 43 FOR PIONS-+)
C     NEL(IN) - NUMBER OF TARGETS (23 FOR NUCLEONS AND 19 FOR PIONS)
      IF(Z1.EQ.0.0.AND.A1.EQ.1.) JS=1
      IF(Z1.EQ.1.0.AND.A1.EQ.1.) JS=2
      IF(Z1.LT.0.0.AND.A1.LE.0.2) JS=3
      IF(Z1.EQ.1.0.AND.A1.LE.0.2) JS=4 
      IN=1
      IF(JS.GT.2)IN=2
C     FLAG JS=1  CROSS SECTIONS FOR NEUTRONS 
C     FLAG JS=2  CROSS SECTIONS FOR PROTONS 
C     FLAG JS=3  CROSS SECTIONS FOR MESONS PI- 
C     FLAG JS=4  CROSS SECTIONS FOR MESONS PI+
      IS1=1
      IS2=2
c      write(6,*) NEL(1),NEL(2),NE(1),NE(2),NE(3),NE(4),JS
      IF(T.GE.1.E6) T=1.E6
C     SELECTION OF ENERGY IENER(LP,JS)<T<IENER(LP+1,JS)
      DO 80 L=1,NE(JS)
      TI=IENER(L,JS)
      IF(TI-T) 80,90,90
  80  CONTINUE
  90  LP=L-1
      IF(L.LE.2) LP=1
      IF(L.GE.NE(JS))LP=NE(JS)-1
      NPAR=1
      NEE=3
C     NPAR=1 FOR NEUTRON,PROTON,PI-,P+.
C     NPAR=2 FOR MESON PI-0
      IF(Z1.EQ.0.0.AND.A1.LE.0.2) NPAR=2
      DO 55 IPAR=1,NPAR
      ZP1=Z1
      IF(IPAR.EQ.1.AND.NPAR.EQ.2) ZP1=-1.
      IF(IPAR.EQ.2.AND.NPAR.EQ.2) ZP1=1.
      DO 50 LI=1,NEE
      L1=LP-1+LI 
      EE(LI)=IENER(L1,JS)
      DO 30 K=1,NEL(IN)
C     JSIG(1)-TOTAL CROSS SECTIONS
C     JSIG(2)-NONELASTIC CROSS SECTIONS
      JSIG(1)=ISIG(JS,K,L1)/10000
      JSIG(2)=ISIG(JS,K,L1)-10000*JSIG(1)
 30   SS(K)=JSIG(IS)      
C     QUADRATIC INTERPOLATION OF CROSS SECTONS VS TARGET MASS
      IF(IN.EQ.1) FF(LI)=BINT(A2,AAI1,SS,NEL(IN),IS2)
      IF(IN.EQ.2) FF(LI)=BINT(A2,AAI2,SS,NEL(IN),IS2)
c      IF(IN.EQ.1) write(6,*)A2,AAI1,SS,NEL(IN),IS2 
c      IF(IN.EQ.2) write(6,*)A2,AAI2,SS,NEL(IN),IS2
  50  CONTINUE  
C     LINEAR INTERPOLATION OF CROSS SECTIONS VS ENERGY
 55   SIG(IPAR)=BINT(T,EE,FF,NEE,IS1)
c      write(6,*) T,EE,FF,NEE,IS1 
      IF(NPAR.EQ.1) SIGHAD=SIG(1)
      IF(NPAR.EQ.2) SIGHAD=(SIG(1)+SIG(2))/2
c      write(6,*) ' sig(1),sig(2) ',sig(1),sig(2)
      IF(SIGHAD) 101,100,100
  101 SIGHAD=0.0
c      WRITE(6,1001)A1,Z1,A2,T,sig(1),sig(2)
      GO TO 100
  102 WRITE(6,1002)A1,Z1,A2,T,sig(1),sig(2)
      GO TO 100
  103 WRITE(6,1003)A1,Z1,A2,T,sig(1),sig(2)
 1001 FORMAT(3F5.1,3E12.4,' ERROR IN INPUT (A1,Z1) OF FUNCTION SIGHAD')
 1002 FORMAT(3F5.1,3E12.4,' ERROR IN INPUT(Proj.ene)OF FUNCTION SIGHAD')
 1003 FORMAT(3F5.1,3E12.4,' ERROR IN INPUT (A2,Z2) OF FUNCTION SIGHAD')
  100 CONTINUE
      RETURN
      END


      BLOCK DATA ASIG
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /BARPO1/ AAI1(24),AAI2(20)     
C     AAI1- for nucleus, AAI2- for pions TARGET NUCLEUS MASS NUMBERS

C    1   H, D, He,  Li,  Be,  C,   N,   O,     Na,  Al,   S,   Ca, 
c       Ti,   Fe,   Cu,   Br,   Mo,  Cd,   Sn,   Ba,     W,     Pb, 
c        U,  Cf      
c
c    3   He,  Be,  C,   N,    O,   Na,  Al,  Ca,  Ti,  Fe,  Cu,  Br, 
c         Mo,  Cd,    Sn,    Ba,     W,    Pb,     U,     Cf          
      DATA AAI1/
     *1.0,2.0,4.0,6.90,9.01,12.00,14.00,16.00,23.00,26.98,32.08,40.08,  
     *47.90,55.85,63.50,79.90,95.94,112.40,118.7,137.34,183.85, 207.19,
     *238.03,250.0/
      DATA AAI2/
     *4.0,9.01,12.0,14.0,16.0,23.0,26.98,40.08,47.90,55.85,63.55,79.90,
     *95.94, 112.40, 118.69, 137.34,183.85,207.19,238.03,250.0/
      END


      SUBROUTINE SIGTAP2 (ITYPE)
      IMPLICIT REAL*8 (A-H,O-Z)
C     READ CROSS SECTIONS FROM BARPOL.DAT
C 2  4 55 62 43 43 24 20 Barashenkov and Polanski 
C (1)Neutrons,(2)Protons (3)Pi- (4)Pi+ Cross Sections
      COMMON /BARPO1/ AAI1(24),AAI2(20)
      COMMON/BARPO/ NEL(2),NE(4),IENER(62,4),ISIG(4,24,62)
      CHARACTER TITLE(60)
      DIMENSION ISIGN(24),ISIGE(24)
      IF(ITYPE.GE.0) RETURN
C     N=2 - LOGICAL UNIT OF TAPE
C     NS=4
C     NE(J) - NUMBER OF ENERGY POINTS 
C     NEL(IN) - NUMBER OF TARGETS
C     IN=1 FOR NUCLEONS
C     IN=2 FOR PIONS 
C     IENER(L,J) THE ENERGY BEAMS IN MEV
C     ISIG(J,K,L)- CROSS SECTIONS
C     (J=1) CROSS SECTIONS FOR NEUTRONS 
C     (J=2) CROSS SECTIONS FOR PROTONS
C     (J=3) CROSS SECTIONS FOR PI- 
C     (J=4) CROSS SECTIONS FOR PI+
C     K- VS. TARGET MASS NUMBERS
C     L- VS. ENERGY
C
      
      OPEN(2,FILE='data/galprop_barpol.dat',STATUS='OLD',
     1 FORM='FORMATTED')
c      OPEN(8,FILE='BARPOL8.DAT',STATUS='NEW',FORM='FORMATTED')
C------------------------------------------------------------------------------------
      INP1=2
      INP2=8
      REWIND INP1       
      IX0=0
      READ(INP1,3001) N,NS,(NE(JS),JS=1,4),(NEL(K),K=1,2),
     1 (TITLE(I),I=1,60)
c      WRITE(INP2,3001)N,NS,(NE(JS),JS=1,4),(NEL(K),K=1,2),
c     1 (TITLE(I),I=1,60) 
      DO 170 JS=1,NS
      IN=1
      IF(JS.GT.2)IN=2
      READ(INP1,3000) IDUM 
c      WRITE(INP2,3000) JS
      DO 170 L=1,NE(JS)
      IF(IN.EQ.1) READ(INP1,5100)IENE,(ISIGE(K),K=1,NEL(IN)),IDUM,IEN,
     1 (ISIGN(K),K=1,NEL(IN)) 
      IF(IN.EQ.2) READ(INP1,4300)IENE,(ISIGE(K),K=1,NEL(IN)),IDUM,IEN,
     1 (ISIGN(K),K=1,NEL(IN))
c      WRITE(INP2,3030)IENE,(ISIGE(K),K=1,NEL(IN)),JS,IEN,
c     1 (ISIGN(K),K=1,NEL(IN))
      IENER(L,JS)=IENE
      IF(L.GT.(NE(JS)-3)) IENER(L,JS)=IENE*1000
      DO 170 K=1,NEL(IN)
      ISIG(JS,K,L)=(ISIGN(K)+ ISIGE(K))*10000+ISIGN(K)
  170 CONTINUE
C---------------------------------------------------------------------
      IF(ITYPE.NE.-2) GO TO 200
      WRITE(INP2,3001)N,NS,(NE(JS),JS=1,4),(NEL(K),K=1,2),
     1 (TITLE(I),I=1,60)     
      DO 180 JS=1,NS
      WRITE(INP2,1015) (IENER(L,JS),L=1,NE(JS))
      IN=1
      IF(JS.GT.2)IN=2
      DO 180 K=1,NEL(IN)
      WRITE(INP2,1015)(ISIG(JS,K,L),L=1,NE(JS)),K,JS 
  180 CONTINUE 
  200 CONTINUE
 3000 FORMAT(20I4)
 5100 FORMAT(51I5)  
 4300 FORMAT(43I5) 
 3030 FORMAT(52I5) 
 3020 FORMAT(A1,2I8)
 3001 FORMAT(8I3,60A1)
 1015 FORMAT(8I9)
      RETURN
      END 
