!!**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
!! * WNEWTR_FUNC_aws.f *                           galprop package * 4/14/2000 
!!**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

C
C Webber Cross section program  - WNEWTR version on floppy dated 8/19/93
C      
C Here the original Webber code has been slightly modified and put into
C form suitable for calling from an arbitrary program.  In fact, only the
C WNEWTR main was altered and split into two subroutines: SET_SIGMA which
C is used to initialize the parameter arrays and WSIGMA which provides the
C special case handling for Beryllium isotopes.  
C
C To use these routines you should first call SET_SIGMA with a logical file
C number and the parameter database filename as arguments.  After the
C parameters are initialized then WSIGMA can be called to return a cross
C section for a particular incident beam (ZI, AI, ENERGY) going to a
C fragment (ZF, AF).  Please note that WSIGMA is a double precision
C function, ENERGY is a double precision argument and ZI, AI, ZF and AF are
C integer arguments.  An example of how these routine are used can be found
C WNEWTR_MAIN.FOR  
C
C---------------------------------------------------------------------
C Modified by A.W.Strong to run correctly on SunSPARC and under fortran90
C 3 changes, labelled 'A.W. Strong Dec 1997' 
c Modified by I.V.Moskalenko to run correctly on PIII under Linux with
c fortran95 NAG compiler, 1/27/2000;
c 1 change, labeled 'imos'
C---------------------------------------------------------------------
C
      SUBROUTINE SET_SIGMA(CDR)                          ! IMOS20020502
C
C  Read in Webber Cross section program initialization parameters
C
C  CDR: Input Device Number (Integer)
C  FILENAME : Webber parameter file name
C
C  NOTE: For this version of the Webber code the parameter file name
C        should be galprop_WNEWTR_082693.CDR.dat
C
        IMPLICIT DOUBLE PRECISION (A-Z) 

        CHARACTER*50 FILENAME
        INTEGER CDR
        INTEGER I,J,K ! A.W. Strong 15 Dec 1997 (required by fortran90)
        DOUBLE PRECISION ODD,SIGMAZF,N0,DELTAZF
        DOUBLE PRECISION ENERGY
        DOUBLE PRECISION EM,DEM,M,EN,DEN,N,EO,DEO,O,B,BEEN,BEDEN,BEN,BEB
!                                                                    ||| A.W. Strong 15 Dec 1997
        COMMON/DATA/ODD,SIGMAZF(100),N0(100),DELTAZF(100)
        COMMON/EDATA/ EM(0:20),DEM(0:20),M(0:20),
     +                EN(0:20),DEN(0:20),N(0:20),
     +                EO(0:20),DEO(0:20),O(0:20),
     +                B(0:20),
     +                BEEN(0:13),BEDEN(0:13),BEN(0:13),BEB(0:13)

C
C Read in the control file
C
        FILENAME="data/galprop_WNEWTR_082693.CDR.dat"                       ! IMOS20020502
        OPEN(UNIT=CDR,FILE=FILENAME,STATUS='OLD',ERR=8)
        do 444 k=1,5   ! IMOS20020325 /comments added to WNEWTR_082693.CDR.dat
           read(CDR,*) !
 444    continue       !
C
C Read in energy (in MeV/Nuc)
C
10      READ(CDR,*) ENERGY
C
C Read in SIGMAZF,N0,DELTAZF table
C
        READ(CDR,*) J
        Z1 = 0
        DO 21 K=1,J
        READ(CDR,*) I, SIGMAZF(I), N0(I), DELTAZF(I)
21      IF( Z1 .EQ. 0 ) Z1 = I
        Z2 = 27
C
C Read in EM,DEM,M,EN,DEN,N,B,EO,DEO,O table
C
        READ(CDR,*) J
        DO 32 I=0,J-1
        READ(CDR,*)EM(I),DEM(I),M(I),EN(I),
     +             DEN(I),N(I),B(I),EO(I),DEO(I),O(I)
   32   CONTINUE
C
C Read in BEEN,BEDEN,BEN table ( special case for BE )
C
        READ(CDR,*) J
        DO 42 I=0,J-1
        READ(CDR,*)BEEN(I),BEDEN(I),BEN(I),BEB(I)
   42   CONTINUE
C
      close(CDR)                          !IMOS20010101
      RETURN
C
8     WRITE(*,9) FILENAME
9     FORMAT(' ?FATAL - Unable to open file ',A20)
      STOP
      END  ! of SET_SIGMA
C
C******************************************************************************
C
        DOUBLE PRECISION FUNCTION WSIGMA(ZI,AI,ZF,AF,ENERGY)
C
C Interface between Webber's main routine, SIGMA, and user applications
C
        DOUBLE PRECISION ENERGY,SIGMA
	INTEGER	ZI,ZF,AI,AF

        WSIGMA = 0.0
	IF (ZI .EQ. 4 .AND. AI .EQ. 8) THEN
          WSIGMA = 0.0
        ELSE IF (ZF .EQ. 4 .AND. AF .EQ. 8) THEN
          WSIGMA = 0.0
	ELSE IF (ZF .EQ. 4 .AND. AF .EQ. 7) THEN
          WSIGMA = SIGMA(ZI,AI,ZF,AF,ENERGY) +
     *             SIGMA(ZI,AI,ZF,AF+1,ENERGY)
        ELSE
	  WSIGMA = SIGMA(ZI,AI,ZF,AF,ENERGY)
        ENDIF
	RETURN
	END	! of FUNCTION WSIGMA
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C Subroutine to calculate total cross section
C
        SUBROUTINE SECTOTAL( IT,ER,SIGTOTA )
        REAL*4 ER(10), SIGTOTA(10)
        REAL AA,AT
        INTEGER IT
        AT = IT
        AA = 5.3 - 2.63*ALOG(AT)
        SIGA = 45.0*AT**0.7*(1.0+0.016*SIN(AA) )
        DO 480 K=1,10
        EE = ER(K)/200.0
        EX = EXP(EE)
        EE = 10.9/ER(K)**0.28
        SIGTOTA(K) = SIGA*(1.0-0.62*EX*SIN(EE) )
480     CONTINUE
        RETURN
        END
C
C Subroutine to write decay constants
C
C
        SUBROUTINE WDECAY(OUT,AI,ZI)

        REAL CONST
        INTEGER OUT,AI,ZI,T,NUL,M
        DATA M/0/
        DATA NUL/0/
        M = M + 1
C
C T = Number of disintegration product
C N = Number of radioactive element(be10,al26,cl36,mn54)
C
        IF( AI .EQ. 10 .AND. ZI .EQ. 4 ) THEN
                CONST = 0.4332E+00
                T = 6
                N = 1
        ELSE IF( AI .EQ. 26 .AND. ZI .EQ. 13 ) THEN
                CONST = 0.7940E+00
                T = 22
                N = 2
        ELSE IF( AI .EQ. 36 .AND. ZI .EQ. 17 ) THEN
                CONST = 0.2265E+01
                T = 36
                N = 3
        ELSE IF( AI .EQ. 54 .AND. ZI .EQ. 25 ) THEN
                CONST = 0.6931E+00
                T = 67
                N = 4
        ELSE
                CONST = 0.0
                T = 0
                N = 0
        END IF

        WRITE(OUT,10) CONST,T,NUL,NUL,N,M
10      FORMAT(E10.4,I5,I5,I10,2I5)
        RETURN
        END

C
C Function SIGMA
C
        DOUBLE PRECISION FUNCTION SIGMA(ZI,AI,ZF,AF,E)
        IMPLICIT DOUBLE PRECISION (A-Z) 

        DOUBLE PRECISION ODD,SIGMAZF,N0,DELTAZF,EDEP
        DOUBLE PRECISION SIGMA0,F,NBAR,DZF,F2,F3,E,E1,E2
        DOUBLE PRECISION TEMP,TEMP2
        INTEGER ZI,ZF,AI,AF,NZI,NZF
        INTEGER I,J,K

        COMMON/DATA/ODD,SIGMAZF(100),N0(100),DELTAZF(100)

C
C Variables:
C               NZI     Neutron excess of initial Z
C               NZF     Neutron excess of final Z
C               ZI,AI   Z and A of Initial particle
C               AF,AF   Z and A of Final particle
C


C
C Figure out if this case if ILLEGAL based on AI,ZI,AF,ZF and return
C sigma = 0 if non valid 
C
        SIGMA = 0.0
        IF( ZF .LT. 4 ) RETURN
        IF( ZF .EQ. 0 ) RETURN
        IF( ZI .LT. ZF )RETURN 
        IF( (ZI .EQ. ZF) .AND. (AI .LE. AF) ) RETURN
        IF( AF .GE. AI ) RETURN
        IF( AI-AF .LT. ZI-ZF ) RETURN
C
C  Calcualte ODD term for SIGMAZF
C
        ODD = 1.0
        IF( ZI .EQ.  9 ) GOTO 1
        GOTO 2
1       IF( (ZI-ZF) .EQ. 1 ) ODD = 0.78
        IF( (ZI-ZF) .EQ. 2 ) ODD = 0.65
2       CONTINUE
C
        IF ( ZI .EQ. 26 ) GOTO 3
        GOTO 4
3       IF( (ZI-ZF) .EQ. 8 ) ODD = 0.93
        IF( (ZI-ZF) .EQ. 10) ODD = 0.90
        IF( (ZI-ZF) .EQ. 12) ODD = 0.60
4       CONTINUE
C
        IF ( ZI .EQ. 20 ) GOTO 5
        GOTO 6
5       IF( (ZI-ZF) .EQ. 5) ODD = 1.20
        IF( (ZI-ZF) .EQ. 6) ODD = 1.20
        IF( (ZI-ZF) .EQ. 7) ODD = 1.15
6       CONTINUE
C
        IF ( ZI .EQ. 18 ) GOTO 7
        GOTO 8
7       IF( (ZI-ZF) .EQ. 2) ODD = 0.90
        IF( (ZI-ZF) .EQ. 4) ODD = 1.10
8       CONTINUE
C
        IF ( ZI .EQ. 16 ) GOTO 9
        GOTO 10
9       IF( (ZI-ZF) .EQ. 4) ODD = 1.18
        IF( (ZI-ZF) .EQ. 7) ODD = 1.50
10      CONTINUE
C
        IF (ZI .EQ. 14) GOTO 13
        GOTO 14
13      IF( (ZI-ZF) .EQ. 3) ODD = 0.70
        IF( (ZI-ZF) .EQ. 4) ODD = 0.84
        IF( (ZI-ZF) .EQ. 7) ODD = 0.80
14      CONTINUE
C
        IF (ZI .EQ. 13) GOTO 15
        GOTO 16
15      IF( (ZI-ZF) .EQ. 1) ODD = 0.83
        IF( (ZI-ZF) .EQ. 2) ODD = 0.60
16      CONTINUE
C
        IF (ZI .EQ. 12) GOTO 17
        GOTO 18
17      IF( (ZI-ZF) .EQ. 1) ODD = 1.10
        IF( (ZI-ZF) .EQ. 3) ODD = 1.18
        IF( (ZI-ZF) .EQ. 7) ODD = 1.45
18      CONTINUE
C       
        IF (ZI .EQ. 11) GOTO 19
        GOTO 20
19      IF( (ZI-ZF) .EQ. 1) ODD = 0.81
        IF( (ZI-ZF) .EQ. 2) ODD = 0.78
20      CONTINUE
C
        IF (ZI .EQ. 10) GOTO 21
        GOTO 22
21      IF( (ZI-ZF) .EQ. 3) ODD = 1.25
22      CONTINUE
C
        IF (ZI .EQ. 8) GOTO 25
        GOTO 26       
25      IF( (ZI-ZF) .EQ. 1) ODD = 1.08
        IF( (ZI-ZF) .EQ. 2) ODD = 0.93
26      CONTINUE

        IF (ZI .EQ. 7) GOTO 27
        GOTO 28
27      IF( (ZI-ZF) .EQ. 1) ODD = 1.07
        IF( (ZI-ZF) .EQ. 2) ODD = 0.80
28      CONTINUE
C
C
        SIGMAZF(4) = 16.40
        IF( ZF .EQ. 4 .AND. ZI .EQ. 5 ) SIGMAZF(4) = 36.0
        NZI = AI-2*ZI
        NZF = AF-2*ZF
        DZF = 0.32*(DBLE(ZF)**0.390)
        IF(ZF.EQ.4) DZF = DZF + 0.20
        F = 1.0

C
C See if we need do Neutron stripping
C
        IF(ZI .NE. ZF) GOTO 40
        IF( AI-AF .EQ. 1 ) GOTO 30
        IF( AI-AF .EQ. 2 ) GOTO 32
        IF( AI-AF .EQ. 3 ) GOTO 34
        IF( AI-AF .ge. 4 ) return    !### imos 1/27/2000 to avoid uncertainty
        GOTO 70  
C
C       Neutron stripping ( 1 )
C
C       SIGMA = (4+0.6Nzi)(5+0.25Zi)( 1 - (Zi - (15 + 2Nzi))/15 )
C                                               ^ This term only used
C                                                 when >= 0 and Zi >= 15
C       SIGMA = SIGMA*( Zi - 26 )
C                          ^for Zi-26 >= 1.0
C
C       for Zi = 7, SIGMA = SIGMA * ( 0.35 + (Nzi/3) )
C
30      TEMP = ( DBLE(ZI) - (15.0 + 2.0*DBLE(NZI)) )/15.0
        IF( (TEMP .LT. 0.0) .OR. (ZI .LT. 15.0) ) TEMP = 0.0
        SIGMA = (1.0-TEMP )*(4.0+0.6*DBLE(NZI))*(5.0+0.25*DBLE(ZI))
        TEMP = DBLE(ZI) - 26.0
        IF(TEMP .LT. 1)TEMP = 1
        SIGMA = SIGMA*TEMP
        IF(ZI .EQ. 7 ) SIGMA = SIGMA * (0.35 + (NZI/3.0) )
        GOTO 90
C
C       Neutron stripping ( 2 )
C
C       SIGMA = 1.8(1 + ( 11.4 - Zi**0.7 )Nzi )
C
C       SIGMA = SIGMA * ( Zi - 26 )
C                           ^ for Zi-26 >= 1.0
C
32      SIGMA = 1.8 * ( 1.0 + (11.4 - DBLE(ZI)**0.7)*DBLE(NZI) )
        TEMP = DBLE(ZI) - 26.0
        IF(TEMP .LT. 1)TEMP = 1
        SIGMA = SIGMA*TEMP
        GOTO 90
C
C       Neutron stripping ( 3 )
C
C       SIGMA = 0.3( 1 + ( 7.5 - Zi**0.5 )Nzi )
C
34      SIGMA = 0.3 * ( 1.0 + (7.5 - DBLE(ZI)**0.5)*DBLE(NZI) )
        GOTO 90

C
C Check for possible Proton stripping
C
40      IF( ZI-ZF .EQ. 1 .AND. AI-AF .EQ. 1) GOTO 50
        IF( ZI-ZF .EQ. 2 .AND. AI-AF .EQ. 2) GOTO 60
        GOTO 70
C
C       Proton stripping ( 1 )
C
C       SIGMA = 0.50 * SIGMAzf * ( 0.7 - 0.05Nzi ) for ZF <> 5
C       SIGMA = 0.78 * SIGMAzf * ( 0.7 - 0.05Nzi ) for ZF = 5
C
C       FOR Zi = 7, SIGMA = SIGMA * 0.4
C
50      TEMP = 0.50
        IF(ZF .EQ. 5)TEMP = 0.78
        SIGMA = TEMP*ODD*SIGMAZF(ZF)*(0.70 - 0.05*DBLE(NZI))
        IF(ZI .EQ. 7) SIGMA = SIGMA*0.4
        GOTO 90
C
C       Proton stripping ( 2 )
C
C       SIGMA = ( 1.85 + (2.5(Zf+2) - 23.0)) * ( 1 - 0.23Nzi )
C                       ^ this term used
C                         if >= 0
C
60      TEMP = 2.5*(DBLE(ZF)+2.0) - 23.0
        IF( TEMP .LT. 0.0 ) TEMP = 0.0
        SIGMA = ( 1.85 + TEMP ) * ( 1.0 - 0.23*DBLE(NZI) )
        GOTO 90
C
C Normal case for SIGMA
C
C
C SIGMA = SIGMA0(ZF,ZI) * F(AF,AI) * Edependance(E,ZF,ZI)
C
C Calculate SIGMA0
C
C       SIGMA0 = SIGMAzf * EXP( -1(ABS(Nzi-N0zf)/8.5)
C                        * EXP( -1(Zi - Zf)/DELTAzf )
C
70      IF ( DELTAZF(ZF) .EQ. 0.0) GOTO 1111  
        SIGMA0 = ODD*SIGMAZF(ZF)*
     +          DEXP(-1.0*( DABS(DBLE(NZI)-N0(ZF)) )/8.5)*
     +          DEXP(-1.0*( DBLE(ZI)-DBLE(ZF) )/DELTAZF(ZF))
        
C
C Calculate F(AF,AI)
C
1111    F2 = DBLE(ZF)-DBLE(ZI)+10.0
        IF(F2.LT.0) F2 = 0
        F3 = DBLE(ZF)-18.0
        IF(F3.LT.0) F3 = 0

        IF( MOD(ZF,2) .EQ. 0 ) GOTO 110
C ZF(odd)
        IF( ZF .EQ. 11 ) GOTO 81
        IF( ZF .EQ. 13 ) GOTO 81
        IF( ZF .EQ. 15 ) GOTO 81
        TEMP = 0.0
        GOTO 82
81      TEMP = -0.07*DBLE(NZI)*(28.0-DBLE(ZI))/10.0
82      F2 = F2**1.4
        NBAR = (0.0400+0.0065*DBLE(NZI))*(DBLE(ZF)**1.2) +
     +          0.06*(F3**0.3) +
     +          DBLE(NZI)*((28.0-DBLE(ZI))/820.0)*F2 + TEMP
        GOTO 111 
C ZF(even)
110     TEMP = 0.0
        IF( ZI .LT. 18 ) GOTO 113
        IF( ZF .EQ. 12 ) GOTO 112
        IF( ZF .EQ. 14 ) GOTO 112
        IF( ZF .EQ. 16 ) GOTO 112
        GOTO 113
112     TEMP = 0.06*((DBLE(NZI)+36.0)/10.0)*(28.0-DBLE(ZI))/10.0
113     TEMP2 = -0.08
        IF( ZF .EQ. 4 .OR. ZF .EQ. 10 ) TEMP2 = 0.0
        NBAR = TEMP2 + (0.0040+0.0009*DBLE(NZI))*(DBLE(ZF)**2  ) -
     +          0.22*(F3**0.5) +
     +          DBLE(NZI)*((28.0-DBLE(ZI))/820.0)*F2 + TEMP
111     CONTINUE
C
C Do Mean mass adjustments to NBAR
C
        IF(ZF .EQ.  4) NBAR = NBAR + 0.27
        IF(ZF .EQ.  5) NBAR = NBAR + 0.36
        IF(ZF .EQ.  6) NBAR = NBAR + 0.10
        IF(ZF .EQ.  7) NBAR = NBAR + 0.12
        IF(ZF .EQ.  8) NBAR = NBAR - 0.15
        IF(ZF .EQ.  9) NBAR = NBAR - 0.05
        IF(ZF .EQ. 10) NBAR = NBAR + 0.27
        IF(ZF .EQ. 11) NBAR = NBAR + 0.05
        IF(ZF .EQ. 12) NBAR = NBAR + 0.02
C       IF(ZF .EQ. 14) NBAR = NBAR + 0.05
        IF(ZF .EQ. 14) NBAR = NBAR - 0.08
        IF(ZF .EQ. 16) NBAR = NBAR - 0.05
        IF(ZF .EQ. 17) NBAR = NBAR - 0.03
        IF(ZF .EQ. 19) NBAR = NBAR + 0.23
        IF(ZF .EQ. 24) NBAR = NBAR - 0.24

        F = (DEXP((-1.0*((DBLE(NZF)-NBAR)**2)) /
     +          (2.0*(DZF**2.0)) ) ) / 
     +          ( DZF * ((2.0*3.1415)**0.5) )
C
C SIGMA = SIGMA0 * F * Energy dependance term
C
C Note: Special case for Neutron/Proton stripping
C       SIGMA = SIGMA * Energy dependance term
C
        SIGMA = SIGMA0 * F
90      SIGMA = SIGMA * EDEP( E, ZI, ZF )
C
C Re-Normalize Energy dependance for E dep at 2000 and 600
C
        E1 = 2000.0
        E2 = 600.0
        SIGMA = SIGMA * EDEP( E1,ZI,ZF ) / EDEP( E2,ZI,ZF )
C
C Special case for 1p + 1n,
C
C       SIGMA = SIGMA + 1.6*(Zi-9)*(1-Ni/4)
C                                    ^ this term used for Ni <= 4
C                                      and Zi <= 9
C
        TEMP = 0.0
        IF( ZI .LE. 9 .OR. ZF .EQ. 4) GOTO 91
        IF(NZI .LE. 4) TEMP = 1.0 - (DBLE(NZI)/4.0)
91      IF(AI-AF .EQ. 2 .AND. ZI-ZF .EQ. 1)
     +                  SIGMA = SIGMA + 1.6*(DBLE(ZI)-9.0)*TEMP

        IF(AF.EQ.10 .AND. ZF.EQ.4 .AND. ZI.GE.7)SIGMA = SIGMA*3.3
C
C ****** Uncomment next line to return ISM cross sections
C
C       SIGMA = 0.93*SIGMA + 0.07*SIGMA*SIGHE( ZI, ZF, E )

9999    RETURN
        END
C
C SIGHE - SIGMA CORRECTION FOR HE ( 4/89 )
C
        DOUBLE PRECISION FUNCTION SIGHE( ZI, ZF, E )

        INTEGER ZI, ZF
        DOUBLE PRECISION E, V, U, F, D, DFE

        V = 1.230

        U = 0.030 * ( 1.0 + 600.0 / E )
        IF( E .LT. 300.0 ) U = 0.090
        IF( E .GT. 1500.0 ) U = 0.042

        F = 1.00 - 0.096 * ( 26 - ZI )

        DFE = 4.40 * ( 1.0 + 0.33 * ( E - 1500.0 ) / 1000.0 )
        IF( E .LT. 300.0 ) DFE = 2.66
        IF( E .GT. 1500.0 ) DFE = 4.40

        D = F * DFE

        SIGHE = DEXP( U * ( DABS((ZI-ZF) - D) )**V )

        RETURN
        END

C
C Energy Dependance term
C
C
C       F(E,ZI,ZI-ZF) = 1 + M_term + N_term + B_term + O_term
C
C
        DOUBLE PRECISION FUNCTION EDEP(E,ZI,ZF)

        DOUBLE PRECISION E,G,H,HO,TEMPA,TEMPB,TEMPC
        DOUBLE PRECISION MTERM,NTERM,OTERM,BTERM
        INTEGER ZI,ZF,DZ
        DOUBLE PRECISION EM,DEM,M,EN,DEN,N,EO,DEO,O,B,BEEN,BEDEN,BEN,BEB
!                                                                    ||| A.W. Strong 15 Dec 1997
        COMMON/EDATA/EM(0:20),DEM(0:20),M(0:20),
     +               EN(0:20),DEN(0:20),N(0:20),
     +               EO(0:20),DEO(0:20),O(0:20),
     +               B(0:20),
     +               BEEN(0:13),BEDEN(0:13),BEN(0:13),BEB(0:13)
C
C Get index into tables ( DZ = ZI - ZF, IF DZ > 18, then DZ = 18 )
C
C Note: Sliding scale for DZ, DZ = (Zi-Zf)*((26/Zi)**0.5)
C
C Also get G = (ZI/26)**2.5, and H = (ZI/26)**1.4, HO = (ZI/26)**0.4
C

        DZ = ZI - ZF
        DZ = INT( FLOAT(DZ)*( (26.0/FLOAT(ZI))**0.5 ) )
        IF( DZ .GT. 18 ) DZ = 18
        G = (DBLE(ZI)/26.0)**2.5
        H = (DBLE(ZI)/26.0)**1.4
        HO = (DBLE(ZI)/26.0)**0.4
        IF(ZI .LE. 26) GOTO 1
        G = 1.0
        H = 1.0
        HO = 1.0
1       CONTINUE
C
C Do special case for ZF = 4 (Be)
C
C Check for special case for BE ( ZF = 4 )
C
C   F = 1.0 + Ndz * EXP( -(E-ENdz)**2 )/DENdz**2 ) + B ( E - 2000 ) / 2000
C
C where EN and DEN are specific for BE ( EN = BEEB, DEN = BEDEN , B = BEB )
C
        IF( ZF .NE. 4 ) GOTO 50
        DZ = ZI - ZF
        IF( DZ .GT. 12 ) DZ = 12
   
        IF (BEDEN(DZ) .NE. 0.0)  GOTO 5958
        TEMPA = -700.0
        GOTO 5959

5958    TEMPA = -1.0*(E - BEEN(DZ))*(E - BEEN(DZ))/(BEDEN(DZ)*BEDEN(DZ))
        IF(TEMPA .LT. -700.0) TEMPA = -700.0

5959    NTERM = 0.0
        IF( BEN(DZ) .NE. 0 ) NTERM = BEN(DZ)*DEXP(TEMPA)

        BTERM = 0.0
        IF((E .LT. 4000.0) .AND. (E .GT. 2000.0))
     +          BTERM = BEB(DZ)*(E-2000.0)/2000.0
        
        EDEP = 1 + NTERM + BTERM

!        print*,'zf E dz beb(dz) edep,nterm,bterm ', zf,E,dz,beb(dz),edep,nterm,bterm
        RETURN
C
C Here for everything else
C
C Calculate TEMPa,TEMPb,TEMPc which are used in Mterm,Nterm,...etc
C
C TEMPA = -1.0 * (E - EMdz)**2 / DEMdz**2
C TEMPB = -1.0 * (E - ENdz)**2 / DENdz**2
C TEMPC = -1.0 * (E - EOdz)**2 / DEOdz**2
C
50      TEMPA =  -1.0 * (E - EM(DZ)) * (E - EM(DZ))/(DEM(DZ)*DEM(DZ))
        TEMPB =  -1.0 * (E - EN(DZ)) * (E - EN(DZ))/(DEN(DZ)*DEN(DZ))

        IF ( DEO(DZ) .NE. 0) GOTO 7979
        TEMPC = -700
        GOTO 7989
7979    TEMPC =  -1.0 * (E - EO(DZ)) * (E - EO(DZ))/(DEO(DZ)*DEO(DZ))
C
C The following is a FUDGE for underflow of TEMPx. If FORTRAN takes
C DEXP(x) where x is less than around -700, the result would be less than
C 10D-38 (pretty small), so our error from this FUDGE can be ignored.
C
7989    IF(TEMPA .LT. -700.0) TEMPA = -700.0
        IF(TEMPB .LT. -700.0) TEMPB = -700.0
        IF(TEMPC .LT. -700.0) TEMPC = -700.0
C
C M_term = Mdz * G * EXP( -(E - Em)**2 / DEMdz**2 )
C
        MTERM = M(DZ)*G*DEXP(TEMPA)
C
C N_term = Ndz * H * EXP( -(E - En)**2 / DENdz**2 )
C
C Note: If Ndz = 0, then N_term = 0
C       if E <= EN(Dz) then EXP(~~~) == 1
C
        NTERM = 0.0
        IF( N(DZ) .EQ. 0 ) GOTO 9087
        IF( E .LE. EN(DZ) ) THEN
            NTERM = N(DZ) * H
        ELSE
            NTERM = N(DZ) * H * DEXP( TEMPB )
        ENDIF
9087    CONTINUE
C
C O_term = Odz * HO * EXP( -(E - Eo)**2 / DEOdz**2 )
C
C Note: if Odz = 0, then O_term = 0
C
        OTERM = 0.0
        IF(O(DZ) .NE. 0 ) OTERM = O(DZ)*HO*DEXP(TEMPC)
C
C B_term = (Zi/26) * Bdz * (E-2000)/2000
C
C        BTERM = 0.0
C        IF( E .GE. 2000 ) THEN
C            TEMP = (E-2000.0)/2000.0
C            IF( TEMP .GT. 1) TEMP = 1
C            BTERM = (DBLE(ZI)/26.0) * B(DZ) * TEMP
C        ENDIF
        TEMP = (E-2000.0)/2000.0
        IF( TEMP .GT. 1.0 ) TEMP = 1.0
        IF( ZI .LT. 26 ) THEN
            BTERM = (DBLE(ZI)/26.0) * B(DZ) * TEMP
        ELSE
            BTERM = B(DZ) * TEMP
        ENDIF
C
C Finally, put it all together
C
        EDEP = 1.0 + MTERM + NTERM + BTERM + OTERM

C        IF( EDEP .GT. 0 ) RETURN
C        WRITE(*,100)ZI,ZF,DZ,E,MTERM,NTERM,OTERM,BTERM,EDEP
100     FORMAT(' ',3I5,2X,6F10.5)
        RETURN
        END
C
C
C
