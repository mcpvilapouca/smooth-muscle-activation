C********************************************************************
C Record of revisions:                                              |
C        Date        Programmer        Description of change        |
C        ====        ==========        =====================        |
C     03/10/2017  Maria Vila Pouca        Mechanical Model          |
C     03/01/2018  Maria Vila Pouca        Chemical Model            |
C--------------------------------------------------------------------
C     Description:
C     UEXTERNALDB: mount database; fibers directions;.inc, .inp files
C     UMAT: NH matrix + Holzapfel fibers + activation
C--------------------------------------------------------------------
      SUBROUTINE UEXTERNALDB(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
      INCLUDE 'PARAM_UMAT.F'
C       this subroutine get the fibers directions resorting the
C      .inc files and read the fibers direction
C
C     UEXTERNAL just called once; work in parallel computing
      COMMON /KFIB/FIBORI
      COMMON /KFIB/FIBORI2
      DIMENSION TIME(2)
      REAL*8 FIBORI(NELEM,4),FIBORI2(NELEM,4)
      CHARACTER(256) FILENAME
      CHARACTER(256) JOBDIR
C     LOP=0 --> START OF THE ANALYSIS
      IF(LOP.EQ.0.OR.LOP.EQ.4) THEN
C
        CALL GETOUTDIR(JOBDIR,LENJOBDIR)
C
          FILENAME=JOBDIR(:LENJOBDIR)//'/'//DIR1
          OPEN(15,FILE=FILENAME)
          DO I=1,NELEM
             READ(15,*) (FIBORI(I,J),J=1,4)
          END DO
           CLOSE(15)
C
          FILENAME=JOBDIR(:LENJOBDIR)//'/'//DIR2
          OPEN(15,FILE=FILENAME)
          DO I=1,NELEM
             READ(15,*) (FIBORI2(I,J),J=1,4)
          END DO
           CLOSE(15)
      END IF
C
      RETURN
C
      END
C
C---------------------------------------------------------------------
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATEV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
      IMPLICIT NONE
      INCLUDE 'PARAM_UMAT.F'
C
      COMMON /KFIB/FIBORI
      COMMON /KFIB/FIBORI2
C
      CHARACTER*8 CMNAME
      DOUBLE PRECISION STRESS(NTENS),STATEV(NSTATEV),
     1 DDSDDT(NTENS),DRPLDE(NTENS),STRAN(NTENS),DSTRAN(NTENS),
     2 TIME(2),PREDEF(1),DPRED(1),PROPS(NPROPS),COORDS(3),
     3 DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),DDSDDE(NTENS,NTENS),
     4 FIBORI(NELEM,4),FIBORI2(NELEM,4),CCMAT1(NTENS,NTENS)
C
      DOUBLE PRECISION SSE, SPD, SCD, RPL, DRPLDT, DTIME, TEMP,
     1   DTEMP,PNEWDT,CELENT
C
      INTEGER NDI, NSHR, NTENS, NSTATEV, NPROPS, NOEL, NPT,
     1 LAYER, KSPT, KSTEP, KINC
C
C---------------------------------------------------------------------
C     LOCAL ARRAYS
C---------------------------------------------------------------------
C     STATEV - STATE VARIABLES ARRAY
C     UNIT   - IDENTITY TENSOR
C     II1    - POINTERS VECTOR
C     DISTGR - DEVIATORIC DEFORMATION GRADIENT (DISTORTION TENSOR)
C     C      -  RIGHT CAUCHY-GREEN TENSOR
C     B      -  LEFT CAUCHY-GREEN TENSOR
C     CBAR   - DEVIATORIC RIGHT CAUCHY-GREEN TENSOR
C     BBAR   - DEVIATORIC LEFT CAUCHY-GREEN TENSOR
C     VORIF  - FIBER ORIENTATION IN UNDERFORMED CONFIGURATION (FAMILY 1)
C     VORIF2 - FIBER ORIENTATION IN UNDERFORMED CONFIGURATION (FAMILY 2)
C     M0,N0  - STRUCTURAL FIBERS TENSOR IN UNDEFORMED  CONFIGURATION
C     M,MM,N,NN - STRUCTURAL FIBERS TENSOR IN DEFORMED  CONFIGURATION
C     MM0,NN0 -FIBERS 4ORDER TENSOR IN UNDEFORMED CONFIGURATION
C     VD     - FIBER ORIENTATION IN DERFORMED CONFIGURATION (FAMILY 1)
C     VD2    - FIBER ORIENTATION IN DERFORMED CONFIGURATION (FAMILY 2)
C     SSX    - AUXILIAR TENSORS
C     SX     - CAUCHY STRESSES OF EACH CONTRIBUTION: SVOL, SM, SF, SDEV
C     SIGMA  - TOTAL CAUCHY STRESS
C     DDSDDEXX - JACOBIAN MATERIAL CONTRIBUTIONS
C---------------------------------------------------------------------
C     LOCAL VARIABLES
C---------------------------------------------------------------------
C     DET       - DETERMINANT OF THE DEFORMATION GRADIENT
C     I1BAR     - 1ST INVARIANT OF CBAR
C     DNORM     - NORM OF THE FIBER DIRECTION VECTOR
C     I4BAR     - SQUARE OF THE FIBER1 STRETCH
C     STRETCH4B - DEVIATORIC STRETCH OF FIBER 1
C     STRETCH4  - STRETCH OF FIBER 1
C     I6BAR     - SQUARE OF THE FIBER2 STRETCH
C     STRETCH6B - DEVIATORIC STRETCH OF FIBER 2
C     STRETCH6  - STRETCH OF FIBER 2
C     DUDXX     - 1ST DERIVATIVE OF PHI RESPECTIVE TO XX
C     D2UDXXX - 2ND DERIVATIVES OF PHI RESPECTIVE TO XXX
C---------------------------------------------------------------------
C                COUNTERS AND STATE INDICATORS
      INTEGER  II1(6),II2(6),ISTAT,INOEL,I,J,K,I1,J1,L,
     1         FLAG1,FLAG2
C
C                KINEMATIC ARRAY VARIABLES
      DOUBLE PRECISION  UNIT(3,3),C(3,3),B(3,3),CBAR(3,3),BBAR(3,3),
     1        DISTGR(3,3),DFGRD1T(3,3),DISTGRT(3,3),CINV(3,3),DETC,
     1        ISTATC,CBARINV(3,3),DETCBAR,ISTATCBAR,
C                FIBERS STRUCTURE VARIABLES
     2        VORIF(3),VD(3),M(3,3),M0(3,3),MM(3,3),
     3        VORIF2(3),VD2(3),N(3,3),N0(3,3),NN(3,3),MM0(3,3,3,3),
     4        NN0(3,3,3,3),S(3,3),
C                STRESS VARIABLES
     5        SBAR(3,3),SBARM(3,3),SBARF4(3,3),SBARF6(3,3),SISO(3,3),
     6        SISOM(3,3),SISOF4(3,3),SISOF6(3,3),SVOL(3,3),
     7        SPK(3,3),SIGMA(3,3), JAC, stressF4(3,3),stressF6(3,3),
C                MATERIAL TANGENT VARIABLES
     8        PP(3,3,3,3),PPMOD(3,3,3,3),PPT(3,3,3,3),CCISO(3,3,3,3),
     9        CCISO1(3,3,3,3),CCISO2(3,3,3,3),CCISO3(3,3,3,3),
     9        CCVOL(3,3,3,3), CCMAT(3,3,3,3), CCSPATIAL(3,3,3,3),
     9        DDSIGDDE(3,3,3,3),ESC,AAUX(3,3,3,3),
     9        DDSDDEJR(3,3,3,3), CCBAR(3,3,3,3)
C
C                MATERIAL PARAMETERS VARIABLES
C                MECHANICAL MODEL
      DOUBLE PRECISION  D1,C10,K11,K12,K21,K22,EOPT,GAMMA,ETA,VC,KXB,
     1                  KC, ZET1, ZET6,
C                CHEMICAL MODEL
     2           AXB, ALPHAA,ALPHAB,ALPHAC,ALPHAD,K1,K2,K3,K7,
     3           ALPHAAn1,ALPHABn1,ALPHACn1,ALPHADn1,CAMLCK,CA2_SCALE,
     4           NK1,CA2,
C                KINEMATIC SCALAR VARIABLES
     2        SCALE,VAR0,VAR1,VAR2,VAR3,VAR4,VAR5,
     3        INV1,I1BAR,PP1,PP2,A,STRETCH4B,STRETCH4,STRETCH6B,
     3        STRETCH6, STRETCH4B2, STRETCH6B2, SUM1, SUM2,
     5        DNORM,DNORM2,ST4INV, ST4INV2, ST6INV, ST6INV2,J23,J43,
     6        EA4, EA6,I4BAR,I6BAR,DTDEA4,DTDEA6,EAIn4,EAIn6,
     7        DTDEAn4,DTDEAn6,Ee4,Ee6,DUDEAI4,DUDEEI4,
C                STRAIN-ENERGY DERIVATIVES VARIABLES
     6        DUDI1, DUDJ, DUDST4, DUDST6, D2UDST4, D2UDST6, D2DJ
C----------------------------------------------------------------------
C     MATERIAL CONSTANTS
C----------------------------------------------------------------------
      D1 = PROPS(1)
      C10 = PROPS(2)
      K11  = PROPS(3)
      K12  = PROPS(4)
      K21  = PROPS(5)
      K22  = PROPS(6)
      EOPT  = PROPS(7)
      GAMMA = PROPS(8)
      ETA = PROPS(9)
      VC = PROPS(10)
      KXB = PROPS(11)
      KC = PROPS(12)
      ZET1 = PROPS(13)
      K2=PROPS(14)
      K3=PROPS(15)
      K7=PROPS(16)
      CAMLCK=PROPS(17)
      NK1=PROPS(18)
      CA2_SCALE=PROPS(19)
C----------------------------------------------------------------------
C     INITIALIZATIONS
C----------------------------------------------------------------------
C     AUXILIAR VARIABLES
      JAC=ZERO
      SCALE=ZERO
      I1BAR=ZERO
      C=ZERO
      CBAR=ZERO
      B=ZERO
      BBAR=ZERO
      CINV=ZERO
      CBARINV=ZERO
      INV1=ZERO
      DUDJ=ZERO
      D2DJ=ZERO
      DUDI1=ZERO
C
      I4BAR=ZERO
      I6BAR=ZERO
      STRETCH4=ONE
      STRETCH4B=ZERO
      STRETCH4B2=ZERO
      STRETCH6=ONE
      STRETCH6B=ZERO
      STRETCH6B2=ZERO
      VAR0=ZERO
      VAR1=ZERO
      VAR2=ZERO
      DUDST4=ZERO
      D2UDST4=ZERO
      VAR3=ZERO
      VAR4=ZERO
      VAR5=ZERO
      DUDST6=ZERO
      D2UDST6=ZERO
      ST4INV=ZERO
      ST6INV=ZERO
      J23=ZERO
      J43=ZERO
      ST4INV2=ZERO
      ST6INV2=ZERO
C
      ISTAT=1
      PP1=ZERO
      PP2=ZERO
      FLAG1=0
      FLAG2=0
C
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C     IDENTITY TENSOR DEFINITION                                       *
C----------------------------------------------------------------------
C
      CAll ONEM(UNIT)
C     INICIALIZAÇAO STATEV
      IF (TIME(1).EQ.0.D0.AND.KSTEP.EQ.1) THEN
c
      STATEV(5)=ZERO
      STATEV(6)=ZERO
      STATEV(7)=ZERO
      STATEV(8)=ZERO
      STATEV(9)=ONE
      STATEV(10)=ZERO
      STATEV(11)=ZERO
      STATEV(12)=ZERO
c
      STATEV(14)=ZERO
      STATEV(15)=ZERO
      STATEV(16)=ZERO
      STATEV(17)=ZERO
      STATEV(18)=ZERO
      STATEV(19)=ZERO
      END IF
C
C
      EA4=STATEV(5)
      EA6=STATEV(6)
      DTDEA4=STATEV(7)
      DTDEA6=STATEV(8)
      ALPHAA=STATEV(9)
      ALPHAB=STATEV(10)
      ALPHAC=STATEV(11)
      ALPHAD=STATEV(12)
C
C----------------------------------------------------------------------
C     POINTERS DEFINITION
C----------------------------------------------------------------------
C     POINTERS TO STORAGE STRESS AND DDSDDE ARRAYS
      II1(1)=1
      II1(2)=2
      II1(3)=3
      II1(4)=1
      II1(5)=1
      II1(6)=2
C
      II2(1)=1
      II2(2)=2
      II2(3)=3
      II2(4)=2
      II2(5)=3
      II2(6)=3
C
C----------------------------------------------------------------------
C        CHEMICAL ACTIVATION (CA2+ CONCENTRATION)
C----------------------------------------------------------------------
          CA2=(PREDEF(1)+DPRED(1))*CA2_SCALE
C       NOTES ABOUT FIELD VARIABLES
c       BEFORE THE STEP, ADD:
c       *initial condition, type=field
c        SET OF NODES, INITIAL VALUE
c        INSIDE THE STEP
C          *field, variable=1, amplitude=ACTIVATION CURVE
C*         *node_set. INICIAL VALUE/FINAL VALUE
           K1=((CA2**(NK1)))/((CA2**(NK1))+(CAMLCK**(NK1)))
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C     CHEMICAL MODEL UPDATE
C----------------------------------------------------------------------
      CALL CHEMICAL(ALPHAA,ALPHAB,ALPHAC,ALPHAD,K1,K2,K3,K7,
     1    ALPHAAn1,ALPHABn1,ALPHACn1,ALPHADn1,DTIME)
C
      ALPHAA=ALPHAAn1
      ALPHAB=ALPHABn1
      ALPHAC=ALPHACn1
      ALPHAD=ALPHADn1
c
C
      AXB=ALPHAC+ALPHAD
C
C----------------------------------------------------------------------
C     JACOBIAN AND DISTORTION TENSOR
C----------------------------------------------------------------------
C     JACOBIAN
      JAC=ZERO
      JAC = DFGRD1(1,1) * DFGRD1(2,2) * DFGRD1(3,3)
     1    - DFGRD1(1,2) * DFGRD1(2,1) * DFGRD1(3,3)
C
      IF (NSHR .EQ. 3) THEN
          JAC = JAC + DFGRD1(1,2) * DFGRD1(2,3) * DFGRD1(3,1)
     1              + DFGRD1(1,3) * DFGRD1(3,2) * DFGRD1(2,1)
     2              - DFGRD1(1,3) * DFGRD1(3,1) * DFGRD1(2,2)
     3              - DFGRD1(2,3) * DFGRD1(3,2) * DFGRD1(1,1)
      END IF
C
C     AUXILIAR VARIABLES
      SCALE = JAC**(-ONE/THREE)
C
      DISTGR=SCALE*DFGRD1
C
C----------------------------------------------------------------------
C     RIGHT AND LEFT CAUCHY-GREEN TENSOR AND INVARIANTS
C----------------------------------------------------------------------
C
C     TRANSPOSE DFGRD1 AND DISTGR
      DFGRD1T=TRANSPOSE(DFGRD1)
      DISTGRT=TRANSPOSE(DISTGR)
C
C     RIGHT AND LEFT CAUCHY-GREEN TENSORS
      CALL M3MULT(DFGRD1T,DFGRD1,C)
      CALL M3MULT(DISTGRT,DISTGR,CBAR)
C
      CALL M3MULT(DFGRD1,DFGRD1T,B)
      CALL M3MULT(DISTGR,DISTGRT,BBAR)
C
C     INVERSE OF LEFT CAUCHY-GREEN TENSOR
      CALL MATINV3D(C,CINV,DETC,ISTATC)
      CALL MATINV3D(CBAR,CBARINV,DETCBAR,ISTATCBAR)
C
C     INVARIANTS OF C AND CBAR
      CALL TRACE(C,INV1)
      I1BAR=JAC**(-TWO/THREE)*INV1
C----------------------------------------------------------------------
C     UNIT VECTOR IN THE DIRECTION OF THE UNDEFORMED FIBERS
C----------------------------------------------------------------------
C
        INOEL=0
        I=0
        DO I=1,NELEM
C               ELEMENT IDENTIFICATION
            IF(NOEL.EQ.INT(FIBORI(I,1))) THEN
                INOEL=I
            ENDIF
        ENDDO
C
C     FIBORI - FIBER ORIENTATION - FAMILY 1
             DNORM=DSQRT(FIBORI(INOEL,2)*FIBORI(INOEL,2)+
     1                   FIBORI(INOEL,3)*FIBORI(INOEL,3)+
     2                   FIBORI(INOEL,4)*FIBORI(INOEL,4))
C
C      FIBORI2 - FIBER ORIENTATION - FAMILY 2
             DNORM2=DSQRT(FIBORI2(INOEL,2)*FIBORI2(INOEL,2)+
     1                   FIBORI2(INOEL,3)*FIBORI2(INOEL,3)+
     2                   FIBORI2(INOEL,4)*FIBORI2(INOEL,4))
C
C       UNDERFORMED FIBER ORIENTATION TENSOR
C
        DO I=1,NDI
        J=I+1
C       FIBER ORIENTATION NORMALIZED VECTOR - FAMILY 1
        VORIF(I)=FIBORI(INOEL,J)/DNORM
C       FIBER ORIENTATION NORMALIZED VECTOR - FAMILY 2
        VORIF2(I)=FIBORI2(INOEL,J)/DNORM2
        END DO
C
      DO I=1,NDI
       DO J=1,NDI
C       STRUCTURAL TENSOR - FAMILY 1
       M0(I,J)=VORIF(I)*VORIF(J)
C       STRUCTURAL TENSOR - FAMILY 2
       N0(I,J)=VORIF2(I)*VORIF2(J)
       END DO
      END DO
C----------------------------------------------------------------------
C     STRETCH RATIO IN THE FIBER DIRECTION
C----------------------------------------------------------------------
C
        I4BAR=ZERO
        I6BAR=ZERO
      DO I=1,NDI
        DO J=1, NDI
            I4BAR=I4BAR+CBAR(I,J)*M0(I,J)
            I6BAR=I6BAR+CBAR(I,J)*N0(I,J)
        ENDDO
      ENDDO
C     FAMILY 1
      STRETCH4B=DSQRT(I4BAR)
      STRETCH4=STRETCH4B/SCALE
C
C     FAMILY 2
      STRETCH6B=DSQRT(I6BAR)
      STRETCH6=STRETCH6B/SCALE
C
C----------------------------------------------------------------------
C     UNIT VECTOR IN THE DIRECTION OF THE DEFORMED FIBERS
C----------------------------------------------------------------------
C
      DO I=1,NDI
         SUM1=ZERO
         SUM2=ZERO
         DO J=1,NDI
          SUM1=SUM1+DFGRD1(I,J)*VORIF(J)
          SUM2=SUM2+DFGRD1(I,J)*VORIF2(J)
         ENDDO
C     FIBER DIRECTIONS IN THE DEFORMED CONFIGURATION
C               -FAMILY 1
         VD(I)=SUM1
C
C               -FAMILY 2
         VD2(I)=SUM2
      ENDDO
      DNORM=DSQRT(VD(1)*VD(1)+
     1             VD(2)*VD(2)+
     2             VD(3)*VD(3))
      DNORM2=DSQRT(VD2(1)*VD2(1)+
     1             VD2(2)*VD2(2)+
     2             VD2(3)*VD2(3))
C           COSINE OF THE ANGLE BETWEEN FIBERS
      DO I=1,NDI
       VD(I)=VD(I)/DNORM
       VD2(I)=VD2(I)/DNORM2
       A=A+VD(I)*VD2(I)
      END DO
C
C
C       STRUCTURE TENSOR IN THE DEFORMED CONFIGURATION - FAMILY 1
      M=MATMUL(M0,DFGRD1T)
      M=MATMUL(DFGRD1,M)
      MM=MATMUL(M0,DISTGRT)
      MM=MATMUL(DISTGR,MM)
C
C       STRUCTURE TENSOR IN THE DEFORMED CONFIGURATION - FAMILY 2
      N=MATMUL(N0,DFGRD1T)
      N=MATMUL(DFGRD1,N)
      NN=MATMUL(N0,DISTGRT)
      NN=MATMUL(DISTGR,NN)
C
C----------------------------------------------------------------------
C     STRAIN-ENERGY DERIVATIVES WITH RESPECT TO INVARIANTS
C----------------------------------------------------------------------
C     FIRST AND SECOND DERIVATIVES PHI - VOLUME-PRESERVING PART
      DUDJ=(ONE/D1)*(JAC-ONE/JAC)
      D2DJ=(ONE/D1)*(ONE+ONE/(JAC**TWO))
C----------------------------------------------------------------------
C     FIRST AND SECOND DERIVATIVES PHI
C     MATRIX CONTRIBUTION
      DUDI1=C10
C     FIBERS CONTRIBUTION
C     FAMILY 1
      STRETCH4B2=STRETCH4B**TWO
      VAR0=(STRETCH4B2-ONE)
      VAR1=VAR0**TWO
      VAR2=EXP((-((EOPT-EA4)**TWO))/(TWO*(GAMMA**TWO)))
      IF(I4BAR.GE.ONE) THEN
C
      DUDST4=TWO*ZET1*K11*VAR0*STRETCH4B*EXP(K12*VAR1)+
     +       ZET1*VAR2*AXB*KXB*(STRETCH4B-ONE-EA4)
C
      D2UDST4=FOUR*ZET1*K11*STRETCH4B2*EXP(K12*VAR1)+
     +       TWO*ZET1*K11*VAR0*EXP(K12*VAR1)+
     +       (8.D0)*ZET1*K11*VAR1*STRETCH4B2*K12*EXP(K12*VAR1)+
     +       ZET1*VAR2*AXB*KXB
      ELSE
      DUDST4=ZET1*VAR2*AXB*KXB*(STRETCH4B-ONE-EA4)
      D2UDST4=ZET1*VAR2*AXB*KXB
      END IF
c
C     FAMILY 2
      STRETCH6B2=STRETCH6B**TWO
      VAR3=(STRETCH6B2-ONE)
      VAR4=VAR3**TWO
      VAR5=EXP((-((EOPT-EA6)**TWO))/(TWO*(GAMMA**TWO)))
      IF(I6BAR.GE.ONE) THEN
C
      DUDST6=TWO*(ONE-ZET1)*K21*VAR3*STRETCH6B*EXP(K22*VAR4)+
     +       (ONE-ZET1)*VAR5*AXB*KXB*(STRETCH6B-ONE-EA6)
C
      D2UDST6=FOUR*(ONE-ZET1)*K21*STRETCH6B2*EXP(K22*VAR4)+
     +       TWO*(ONE-ZET1)*K21*VAR3*EXP(K22*VAR4)+
     +       (8.D0)*(ONE-ZET1)*K21*VAR4*STRETCH6B2*K22*EXP(K22*VAR4)+
     +       (ONE-ZET1)*VAR5*AXB*KXB
      ELSE
      DUDST6=(ONE-ZET1)*VAR5*AXB*KXB*(STRETCH6B-ONE-EA6)
      D2UDST6=(ONE-ZET1)*VAR5*AXB*KXB
      END IF
c
C
C----------------------------------------------------------------------
C     SECOND PIOLA-KIRCHHOFF STRESS TENSOR
C----------------------------------------------------------------------
C     AUXILIAR VARIABLES
      ST4INV=(STRETCH4B)**(-ONE)
      ST6INV=(STRETCH6B)**(-ONE)
      J23=JAC**(-TWO/THREE)
C
C     ISOTROPIC CONTRIBUTION (SBAR)
      SBARM=TWO*DUDI1*UNIT
      SBARF4=DUDST4*(ST4INV)*M0
      SBARF6=DUDST6*((STRETCH6B)**(-ONE))*N0
      SBAR=SBARM+SBARF4+SBARF6
C     ISOTROPIC CONTRIBUTION (SISO)
      SISOM=DUDI1*(TWO*J23*UNIT-(TWO/THREE)*I1BAR*CINV)
      SISOF4=DUDST4*(J23*ST4INV*M0-(ONE/THREE)*STRETCH4B*CINV)
      SISOF6=DUDST6*(J23*ST6INV*N0-(ONE/THREE)*STRETCH6B*CINV)
      SISO=SISOM+SISOF4+SISOF6
C     VOLUMETRIC CONTRIBUTION (SVOL)
      SVOL=JAC*DUDJ*CINV
c
C
C     SECOND PIOLA-KIRCHHOFF STRESS TENSOR (SPK)
      S=SISO+SVOL
C----------------------------------------------------------------------
C     CAUCHY STRESS TENSOR
C----------------------------------------------------------------------
      CALL PUSH2(S,DFGRD1,JAC,SIGMA,NDI)
      CALL PUSH2(SISOF4,DFGRD1,JAC,stressF4,NDI)
      CALL PUSH2(SISOF6,DFGRD1,JAC,stressF6,NDI)
C----------------------------------------------------------------------
C     JACOBIAN MATERIAL MATRIX
C----------------------------------------------------------------------
C     AUXILIAR VARIABLES
      J43=JAC**(-FOUR/THREE)
      ST4INV2=STRETCH4B**(-TWO)
      ST6INV2=STRETCH6B**(-TWO)
C
C     MATERIAL MATRIX(CCBAR)
      CALL TENSORPROD2(M0,M0,MM0,NDI)
      CALL TENSORPROD2(N0,N0,NN0,NDI)
      DO I=1,NDI
       DO J=1,NDI
         DO K=1,NDI
          DO L=1,NDI
          CCBAR(I,J,K,L)=J43*ST4INV2*(D2UDST4-ST4INV*DUDST4)*
     *    MM0(I,J,K,L)+J43*ST6INV2*(D2UDST6-ST6INV*DUDST6)*NN0(I,J,K,L)
          END DO
         END DO
       END DO
      END DO
C     MATERIAL MATRIX(CCISO)
      CALL PROJECTION(C,CINV,NDI,PP)
      CALL PROJECTIONMOD(CINV,NDI,PPMOD)
      CALL PROJECTIONT(C,CINV,NDI,PPT)
C
      CALL DOUBLEDOT4(PP,CCBAR,AAUX,NDI)
      CALL DOUBLEDOT4(AAUX,PPT,CCISO1,NDI)
C
      CALL DOUBLEDOT2(SBAR,C,ESC,NDI)
C
      CALL TENSORPROD2(CINV,SISO,CCISO2,NDI)
      CALL TENSORPROD2(SISO,CINV,CCISO3,NDI)
C
      CCISO=CCISO1+(TWO/THREE)*J23*ESC*PPMOD-(TWO/THREE)*(CCISO2+CCISO3)
C
C     VOLUMETRIC CONTRIBUTION
      DO I=1,NDI
       DO J=1,NDI
         DO K=1,NDI
          DO L=1,NDI
          CCVOL(I,J,K,L)=(TWO/D1)*(JAC**TWO)*(CINV(I,J)*CINV(K,L))-
     -     (ONE/D1)*((JAC**TWO)-ONE)*(CINV(I,K)*CINV(J,L)+CINV(I,L)*
     *     CINV(J,K))
          END DO
         END DO
       END DO
      END DO

C     MATERIAL DDSDDE
      CCMAT=CCISO+CCVOL
c
      CALL PUSH4(CCMAT,DFGRD1,JAC,CCSPATIAL,NDI)
C
C     SPATIAL TANGENT MODULI BASED ON THE JAUMANN
C                        RATE OF THE KIRCHHOFF STRESS TENSOR
      DO I=1,NDI
         DO J=1,NDI
            DO K=1,NDI
               DO L=1,NDI
C         ------JAUMMAN RATE PART------
                   DDSDDEJR(I,J,K,L)=
     +             (ONE/TWO)*(UNIT(I,K)*SIGMA(J,L)
     +             +SIGMA(I,K)*UNIT(J,L)+UNIT(I,L)*SIGMA(J,K)
     +             +SIGMA(I,L)*UNIT(J,K))
C        -----SPATIAL TANGENT MODULI------
                   DDSIGDDE(I,J,K,L)=DDSDDEJR(I,J,K,L)+
     +                               CCSPATIAL(I,J,K,L)
               END DO
             END DO
          END DO
       END DO
C----------------------------------------------------------------------
C     STRESS AND JACOBIAN MATRIX STORAGE
C----------------------------------------------------------------------

      DO I1=1,NTENS
C       STRESS VECTOR
         STRESS(I1)=SIGMA(II1(I1),II2(I1))
         DO J1=1,NTENS
C       DDSDDE - FULLY SIMMETRY IMPOSED
            PP1=DDSIGDDE(II1(I1),II2(I1),II1(J1),II2(J1))
            PP2=DDSIGDDE(II1(I1),II2(I1),II2(J1),II1(J1))
            DDSDDE(I1,J1)=(ONE/TWO)*(PP1+PP2)
         END DO
      END DO
            DO I1=1,NTENS
         DO J1=1,NTENS
C       DDSDDE - FULLY SIMMETRY IMPOSED
            PP1=CCMAT(II1(I1),II2(I1),II1(J1),II2(J1))
            PP2=CCMAT(II1(I1),II2(I1),II2(J1),II1(J1))
            CCMAT1(I1,J1)=(ONE/TWO)*(PP1+PP2)
         END DO
      END DO
C
C----------------------------------------------------------------------
C     EAI UPDATE
C----------------------------------------------------------------------
C     FAMILY 1
      CALL EAI_UPDATE(EA4,DTDEA4,STRETCH4B,VC,ETA,GAMMA,EOPT,KC,
     1             ALPHAC,DTIME,EAIn4,ZET1,AXB,KXB,DTDEAn4)
c
      DTDEA4=DTDEAn4
      EA4=EAIn4
      Ee4=STRETCH4B-ONE-EA4
c     FAMILY 2
      ZET6=ONE-ZET1
      CALL EAI_UPDATE(EA6,DTDEA6,STRETCH6B,VC,ETA,GAMMA,EOPT,KC,
     1             ALPHAC,DTIME,EAIn6,ZET6,AXB,KXB,DTDEAn6)

      DTDEA6=DTDEAn6
      EA6=EAIn6
      Ee6=STRETCH6B-ONE-EA6
c
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C     STATE VARIABLES
C----------------------------------------------------------------------
C
        STATEV(1)=  JAC
        STATEV(2)=  A
        STATEV(3)=  STRETCH4B
        STATEV(4)=  STRETCH6B
        STATEV(5)=  EA4
        STATEV(6)=  EA6
        STATEV(7)=  DTDEA4
        STATEV(8)=  DTDEA6
        STATEV(9)=  ALPHAA
        STATEV(10)= ALPHAB
        STATEV(11)= ALPHAC
        STATEV(12)= ALPHAD
        STATEV(13)= CA2
c
C     DIRECTIONS OF FIBERS IN THE DEFORMED CONFIGURATION
        STATEV(14)=VD(1)
        STATEV(15)=VD(2)
        STATEV(16)=VD(3)
        STATEV(17)=VD2(1)
        STATEV(18)=VD2(2)
        STATEV(19)=VD2(3)
C      END DO
C----------------------------------------------------------------------
      RETURN
      END
C----------------------------------------------------------------------
C END OF MAIN UMAT ROUTINE
C----------------------------------------------------------------------
C****************************************************************************
C     UTILITY SUBROUTINES
C****************************************************************************
C----------------------------------------------------------------------
       SUBROUTINE CHEMICAL(ALPHAA,ALPHAB,ALPHAC,ALPHAD,K1,K2,K3,K7,
     1    ALPHAAn1,ALPHABn1,ALPHACn1,ALPHADn1,DTIME)
c
       Implicit None
c
       DOUBLE PRECISION ALPHAA,ALPHAB,ALPHAC,ALPHAD,K1,K2,K3,K7,K4,
     1   ALPHAAn,ALPHABn,ALPHACn,ALPHADn,ALPHAAp,ALPHABp,ALPHACp,
     2   ALPHADp,FAn,FBn,FCn,FDn,FAp,FBp,FCp,FDp,DTIME,
     3   ALPHAAn1,ALPHABn1,ALPHACn1,ALPHADn1
C
      K4=K3/(4.D0)
      ALPHAAn=ALPHAA
      ALPHABn=ALPHAB
      ALPHACn=ALPHAC
      ALPHADn=ALPHAD
C
      FAn=(-K1)*(ALPHAAn)+K2*(ALPHABn)+K7*(ALPHADn)
      FBN=K1*(ALPHAAn)-(K2+K3)*(ALPHABn)+K4*(ALPHACn)
      FCN=K3*(ALPHABn)-(K4+K2)*(ALPHACn)+K1*(ALPHADn)
      FDN=(K2*ALPHACn)-(K1+K7)*(ALPHADn)
C
C     PREDICTOR
      ALPHAAp=ALPHAAn+DTIME*FAn
      ALPHABp=ALPHABn+DTIME*FBn
      ALPHACp=ALPHACn+DTIME*FCn
      ALPHADp=ALPHADn+DTIME*FDn
C
      FAp=(-K1)*(ALPHAAp)+K2*(ALPHABp)+K7*(ALPHADp)
      FBp=K1*(ALPHAAp)-(K2+K3)*(ALPHABp)+K4*(ALPHACp)
      FCp=K3*(ALPHABp)-(K4+K2)*(ALPHACp)+K1*(ALPHADp)
      FDp=K2*(ALPHACp)-(K1+K7)*(ALPHADp)
C
C     NEXT INCREMENT N+1
      ALPHAAn1=((1.D0)/(2.D0))*(ALPHAAn+ALPHAAp+DTIME*FAp)
      ALPHABn1=((1.D0)/(2.D0))*(ALPHABn+ALPHABp+DTIME*FBp)
      ALPHACn1=((1.D0)/(2.D0))*(ALPHACn+ALPHACp+DTIME*FCp)
      ALPHADn1=((1.D0)/(2.D0))*(ALPHADn+ALPHADp+DTIME*FDp)
C
      RETURN
C
      end SUBROUTINE CHEMICAL
C----------------------------------------------------------------------
       SUBROUTINE DOUBLEDOT2(A,B,C,NDI)
C
       Implicit None
C
       INTEGER I,J,K,L,NDI
C
       DOUBLE PRECISION A(3,3),B(3,3),C
C
      C=0.D0
      DO I=1,NDI
       DO J=1,NDI
         C=C+A(I,J)*B(J,I)
       END DO
      END DO
C
      RETURN
C
      end SUBROUTINE DOUBLEDOT2
C----------------------------------------------------------------------
       SUBROUTINE DOUBLEDOT4(A,B,C,NDI)
C
       Implicit None
C
       INTEGER I,J,K,L,NDI, M, N
C
       DOUBLE PRECISION A(3,3,3,3),B(3,3,3,3),C(3,3,3,3),SUM
C
      DO I=1,NDI
       DO J=1,NDI
         DO K=1,NDI
          DO L=1,NDI
           SUM=0.D0
           DO M=1,NDI
            DO N=1,NDI
             SUM=SUM+A(I,J,M,N)*B(M,N,K,L)
            END DO
           END DO
           C(I,J,K,L)=SUM
          END DO
         END DO
       END DO
      END DO
C
      RETURN

      end SUBROUTINE DOUBLEDOT4
C----------------------------------------------------------------------
       SUBROUTINE EAI_UPDATE(EAI,DTDEAI,STRETCH,VC,ETA,GAMMA,EOPT,KC,
     1             ALPHAC,DTIME,EAIn1,ZETI,AXB,KXB,DTDEAIn1)
c
       Implicit None
c
       INTEGER I,J,K,L,NDI
c
       DOUBLE PRECISION EAI,VC,ETA,GAMMA,EOPT,KC,ALPHAC, STRETCH, A1n,
     1  A2,A3n,DTIME,A4n,DUDEAIn,Nin,EEIn,EEIp,DUDEEIn,A1p,A3p,A4p,
     2  DUDEAIp,DUDEEIp,Nip,Fn,Fp,ZETI,DTDEAI,EAIn1,EAIn,AXB,KXB,EAIp,
     3  DTDEAIn1,EEIn1,Nin1,DUDEAIn1,DUDEEIn1,A3n1,A4n1
c
       EAIn=EAI
       EEIn=STRETCH-(1.D0)-EAIn
c
C
      IF (DTDEAI.GT.-VC) THEN
C      AUXILIAR VARIABLES
       A2=(2.D0)*((GAMMA)**(2.D0))
       A3n=(-(EOPT-EAIn)**(2.D0))
        Nin=EXP(A3n/A2)
        A4n=Nin*KC*ALPHAC
        DUDEAIn=(ZETI*(EOPT-EAIn)*Nin*
     *           AXB*KXB*((EEIn)**(2.d0)))/A2
        DUDEEIn=ZETI*Nin*AXB*KXB*EEIn
c
C
        Fn=((ETA+A4n)**(-1.D0))*(-(A4n*VC)-DUDEAIn+DUDEEIn)
C     PREDICTOR
        EAIp=EAIn+DTIME*Fn
        EEIp=STRETCH-(1.D0)-EAIp
c
        A3p=(-(EOPT-EAIp)**(2.D0))
        Nip=EXP(A3p/A2)
        A4p=Nip*KC*ALPHAC
        DUDEAIp=(ZETI*(EOPT-EAIp)*Nip*
     *           AXB*KXB*((EEIp)**(2.d0)))/A2
        DUDEEIp=ZETI*Nip*AXB*KXB*EEIp
c
        Fp=((ETA+A4p)**(-1.D0))*(-(A4p*VC)-DUDEAIp+DUDEEIp)
c
       EAIn1=(1.D0/2.D0)*(EAIn+EAIp+DTIME*Fp)
       EEIn1=STRETCH-(1.D0)-EAIn1
c
        A3n1=(-(EOPT-EAIn1)**(2.D0))
        Nin1=EXP(A3n1/A2)
        A4n1=Nin1*KC*ALPHAC
        DUDEAIn1=(ZETI*(EOPT-EAIn1)*Nin1*
     *           AXB*KXB*((EEIn1)**(2.d0)))/A2
        DUDEEIn1=ZETI*Nin1*AXB*KXB*EEIn1
c
C
        DTDEAIn1=((ETA+A4n1)**(-1.D0))*(-(A4n1*VC)-DUDEAIn1+DUDEEIn1)
c
       ELSE
C      AUXILIAR VARIABLES
       A2=(2.D0)*((GAMMA)**(2.D0))
       A3n=-(EOPT-EAIn)**(2.D0)
C
        DUDEAIn=(ZETI*(EOPT-EAIn)*EXP(A3n/A2)*
     *           AXB*KXB*((EEIn)**(2.d0)))/A2
        DUDEEIn=ZETI*EXP(A3n/A2)*AXB*KXB*EEIn
C
        Fn=((ETA)**(-1.D0))*(-DUDEAIn+DUDEEIn)
C
C     PREDICTOR
        EAIp=EAIn+DTIME*Fn
        EEIp=STRETCH-(1.D0)-EAIp
c
       A3p=-(EOPT-EAIp)**(2.D0)
c
       DUDEAIp=(ZETI*(EOPT-EAIp)*EXP(A3p/A2)*
     *           AXB*KXB*((EEIp)**(2.d0)))/A2
       DUDEEIp=ZETI*EXP(A3p/A2)*AXB*KXB*EEIp
c
       Fp=((ETA)**(-1.D0))*(-DUDEAIp+DUDEEIp)
C
       EAIn1=(1.D0/2.D0)*(EAIn+EAIp+DTIME*Fp)
       EEIn1=STRETCH-(1.D0)-EAIn1
       A3n1=(-(EOPT-EAIn1)**(2.D0))
        Nin1=EXP(A3n1/A2)
        DUDEAIn1=(ZETI*(EOPT-EAIn1)*Nin1*
     *           AXB*KXB*((EEIn1)**(2.d0)))/A2
        DUDEEIn1=ZETI*Nin1*AXB*KXB*EEIn1
c
        DTDEAIn1=((ETA)**(-1.D0))*(-DUDEAIn1+DUDEEIn1)
       END IF
C
      RETURN
C
      end SUBROUTINE EAI_UPDATE
C----------------------------------------------------------------------
      SUBROUTINE M3MULT(A,B,C)
C
      IMPLICIT NONE
C
      DOUBLE PRECISION A(3,3),B(3,3), C(3,3), SUM
      INTEGER I,J,K
C
      DO I=1,3
       DO J=1,3
        SUM=0.D0
          DO K=1,3
           SUM=SUM+A(I,K)*B(K,J)
          END DO
        C(I,J)=SUM
       END DO
      END DO
C
      RETURN
C
      END SUBROUTINE M3MULT
C----------------------------------------------------------------------
       SUBROUTINE PROJECTION(C,CINV,NDI,PP)
C
       Implicit None
C
       INTEGER I,J,K,L,NDI
C
       DOUBLE PRECISION C(3,3),CINV(3,3),PP(3,3,3,3)
C
      DO I=1,NDI
       DO J=1,NDI
         DO K=1,NDI
          DO L=1,NDI
          IF ((I.EQ.K).AND.(J.EQ.L)) THEN
          PP(I,J,K,L)=1.D0-((1.D0)/(3.D0))*CINV(I,J)*C(K,L)
          ELSE
          PP(I,J,K,L)=-((1.D0)/(3.D0))*CINV(I,J)*C(K,L)
          END IF
          END DO
         END DO
       END DO
      END DO
C
      RETURN
C
      end SUBROUTINE PROJECTION
C----------------------------------------------------------------------
       SUBROUTINE PROJECTIONMOD(CINV,NDI,PPMOD)
C
       Implicit None
C
       INTEGER I,J,K,L,NDI
C
       DOUBLE PRECISION CINV(3,3),PPMOD(3,3,3,3)
C
      DO I=1,NDI
       DO J=1,NDI
         DO K=1,NDI
          DO L=1,NDI
          PPMOD(I,J,K,L)=((1.D0)/(2.D0))*(CINV(I,K)*CINV(J,L)+
     +      CINV(I,L)*CINV(J,K))-((1.D0)/(3.D0))*(CINV(I,J)*CINV(K,L))
          END DO
         END DO
       END DO
      END DO
C
      RETURN
C
      end SUBROUTINE PROJECTIONMOD
C----------------------------------------------------------------------
       SUBROUTINE PROJECTIONT(C,CINV,NDI,PPT)
C
       Implicit None
C
       INTEGER I,J,K,L,NDI
C
       DOUBLE PRECISION C(3,3),CINV(3,3),PPT(3,3,3,3)
C
      DO I=1,NDI
       DO J=1,NDI
         DO K=1,NDI
          DO L=1,NDI
          IF ((I.EQ.K).AND.(J.EQ.L)) THEN
          PPT(I,J,K,L)=1.D0-((1.D0)/(3.D0))*C(I,J)*CINV(K,L)
          ELSE
          PPT(I,J,K,L)=-((1.D0)/(3.D0))*C(I,J)*CINV(K,L)
          END IF
          END DO
         END DO
       END DO
      END DO
C
      RETURN
C
      end SUBROUTINE PROJECTIONT
C----------------------------------------------------------------------
       SUBROUTINE PULL2(A,DFGRD1,DET,PULL2A,NDI)
C
       Implicit None
C
       INTEGER NDI,I,J,K,L,ISTAT
C
       DOUBLE PRECISION A(NDI,NDI),DFGRD1(NDI,NDI),PULL2A(NDI,NDI),DET,
     1                  DFGRD1_INV(NDI,NDI),DET_DFGRD1,CONT
C
       CAll MATINV3D(DFGRD1,DFGRD1_INV,DET_DFGRD1,ISTAT)
C
      DO I=1,NDI
       DO J=1,NDI
       CONT=0.D0
         DO K=1,NDI
          DO L=1,NDI
          CONT=CONT+DET*(DFGRD1_INV(I,K)*A(K,L)*
     *      DFGRD1_INV(J,L))
          END DO
         END DO
         PULL2A(I,J)=CONT
       END DO
      END DO
C
      RETURN
C
      end SUBROUTINE
C----------------------------------------------------------------------
      SUBROUTINE PULL4(CC,DFGRD1,DET,PULL4CC,NDI)
C
       Implicit None
C
       INTEGER NDI,I,J,K,L,A,B,C,D,ISTAT
C
       DOUBLE PRECISION CC(NDI,NDI,NDI,NDI),DFGRD1(NDI,NDI),
     1                  PULL4CC(NDI,NDI,NDI,NDI),DET,
     1                  DFGRD1_INV(NDI,NDI),DET_DFGRD1,CONT
C
       CAll MATINV3D(DFGRD1,DFGRD1_INV,DET_DFGRD1,ISTAT)
C
      DO I=1,NDI
       DO J=1,NDI
         DO K=1,NDI
          DO L=1,NDI
          CONT=0.D0
          DO A=1,NDI
           DO B=1,NDI
            DO C=1,NDI
             DO D=1,NDI
             CONT=CONT+DET*(DFGRD1_INV(I,A)*DFGRD1_INV(J,B)*
     *       DFGRD1_INV(K,C)*DFGRD1_INV(L,D)*CC(A,B,C,D))
             END DO
            END DO
           END DO
          END DO
         PULL4CC(I,J,K,L)=CONT
         END DO
        END DO
       END DO
      END DO
C
      RETURN
C
      end SUBROUTINE
C----------------------------------------------------------------------
       SUBROUTINE PUSH2(A,DFGRD1,DET,PUSH2A,NDI)
C
       Implicit None
C
       INTEGER NDI,I,J,K,L
C
       DOUBLE PRECISION A(NDI,NDI),DFGRD1(NDI,NDI),PUSH2A(NDI,NDI),DET,
     1                  CONT
C
      DO I=1,NDI
       DO J=1,NDI
       CONT=0.D0
         DO K=1,NDI
          DO L=1,NDI
          CONT=CONT+(1.d0/DET)*(DFGRD1(I,K)*A(K,L)*
     *      DFGRD1(J,L))
          END DO
         END DO
         PUSH2A(I,J)=CONT
       END DO
      END DO
C
      RETURN
C
      end SUBROUTINE
C----------------------------------------------------------------------
      SUBROUTINE PUSH4(CC,DFGRD1,DET,PUSH4CC,NDI)
C
       Implicit None
C
       INTEGER NDI,I,J,K,L,A,B,C,D
C
       DOUBLE PRECISION CC(NDI,NDI,NDI,NDI),DFGRD1(NDI,NDI),
     1                  PUSH4CC(NDI,NDI,NDI,NDI),DET,CONT
C
      DO I=1,NDI
       DO J=1,NDI
         DO K=1,NDI
          DO L=1,NDI
          CONT=0.D0
          DO A=1,NDI
           DO B=1,NDI
            DO C=1,NDI
             DO D=1,NDI
             CONT=CONT+(1.D0/DET)*(DFGRD1(I,A)*DFGRD1(J,B)*
     *       DFGRD1(K,C)*DFGRD1(L,D)*CC(A,B,C,D))
             END DO
            END DO
           END DO
          END DO
         PUSH4CC(I,J,K,L)=CONT
         END DO
        END DO
       END DO
      END DO
C
      RETURN
C
      end SUBROUTINE
C----------------------------------------------------------------------
       SUBROUTINE TENSORPROD2(A,B,C,NDI)
C
       Implicit None
C
       INTEGER I,J,K,L,NDI
C
       DOUBLE PRECISION A(3,3),B(3,3),C(3,3,3,3)
C
      DO I=1,NDI
       DO J=1,NDI
         DO K=1,NDI
          DO L=1,NDI
          C(I,J,K,L)=A(I,J)*B(K,L)
          END DO
         END DO
       END DO
      END DO
C
      RETURN
C
      end SUBROUTINE TENSORPROD2
C----------------------------------------------------------------------
      SUBROUTINE TRACE(A,SUM)
C
      IMPLICIT NONE
C
      DOUBLE PRECISION A(3,3),SUM
      INTEGER I,J
C
      SUM=0.D0
      DO I=1,3
       DO J=1,3
       IF (I.EQ.J) THEN
        SUM=SUM+A(I,J)
       END IF
       END DO
      END DO
C
      RETURN
C
      END SUBROUTINE TRACE
C----------------------------------------------------------------------
      SUBROUTINE MATINV3D(A,A_INV,DET_A,ISTAT)
C
C      RETURNS A_INV, THE INVERSE AND DET_A, THE DETERMINANT
C     DET OF THE ORIGINAL MATRIX, NOT OF THE INVERSE
C      RECURSIVE EXPRESSION OBTAINED FROM MATHEMATICA
      IMPLICIT NONE
C
      INTEGER ISTAT
C
      REAL*8 A(3,3),A_INV(3,3),DET_A,DET_A_INV
C
C
      ISTAT = 1
C
      DET_A = A(1,1)*(A(2,2)*A(3,3) - A(3,2)*A(2,3)) -
     +        A(2,1)*(A(1,2)*A(3,3) - A(3,2)*A(1,3)) +
     +        A(3,1)*(A(1,2)*A(2,3) - A(2,2)*A(1,3))
C
      IF (DET_A .LE. 0.D0) THEN
        WRITE(*,*) 'WARNING: SUBROUTINE MATINV3D:'
        WRITE(*,*) 'WARNING: DET OF MAT=',DET_A
        ISTAT = 0
        RETURN
      END IF
C
      DET_A_INV = 1.D0/DET_A
C
      A_INV(1,1) = DET_A_INV*(A(2,2)*A(3,3)-A(3,2)*A(2,3))
      A_INV(1,2) = DET_A_INV*(A(3,2)*A(1,3)-A(1,2)*A(3,3))
      A_INV(1,3) = DET_A_INV*(A(1,2)*A(2,3)-A(2,2)*A(1,3))
      A_INV(2,1) = DET_A_INV*(A(3,1)*A(2,3)-A(2,1)*A(3,3))
      A_INV(2,2) = DET_A_INV*(A(1,1)*A(3,3)-A(3,1)*A(1,3))
      A_INV(2,3) = DET_A_INV*(A(2,1)*A(1,3)-A(1,1)*A(2,3))
      A_INV(3,1) = DET_A_INV*(A(2,1)*A(3,2)-A(3,1)*A(2,2))
      A_INV(3,2) = DET_A_INV*(A(3,1)*A(1,2)-A(1,1)*A(3,2))
      A_INV(3,3) = DET_A_INV*(A(1,1)*A(2,2)-A(2,1)*A(1,2))
C
C
      RETURN
      END SUBROUTINE MATINV3D
C----------------------------------------------------------------------
      SUBROUTINE MATINV2D(A,A_INV,DET_A,ISTAT)
C
C     RETURNS A_INV, THE INVERSE, AND DET_A, THE DETERMINANT
C     NOTE THAT THE DET IS OF THE ORIGINAL MATRIX, NOT THE
C     INVERSE
C
      IMPLICIT NONE
C
      INTEGER ISTAT
C
      REAL*8 A(2,2),A_INV(2,2),DET_A,DET_A_INV
C
C
      ISTAT = 1.D0
C
      DET_A = A(1,1)*A(2,2) - A(1,2)*A(2,1)
C
      IF (DET_A .LE. 0.D0) THEN
        WRITE(*,*) 'WARNING: SUBROUTINE MATINV2D:'
        WRITE(*,*) 'WARNING: DET OF MAT=',DET_A
        ISTAT = 0
        RETURN
      END IF
C
      DET_A_INV = 1.D0/DET_A
C
      A_INV(1,1) =  DET_A_INV*A(2,2)
      A_INV(1,2) = -DET_A_INV*A(1,2)
      A_INV(2,1) = -DET_A_INV*A(2,1)
      A_INV(2,2) =  DET_A_INV*A(1,1)
C
C
      RETURN
      END SUBROUTINE MATINV2D
C----------------------------------------------------------------------
C
      SUBROUTINE MDET(A,DET)
C
C      THIS SUBROUTINE CALCULATES THE DETERMINANT
C      OF A 3 BY 3 MATRIX [A]
C
      IMPLICIT NONE
C
      REAL*8  A(3,3),DET
C
C
      DET = A(1,1)*A(2,2)*A(3,3)
     +	  + A(1,2)*A(2,3)*A(3,1)
     +	  + A(1,3)*A(2,1)*A(3,2)
     +	  - A(3,1)*A(2,2)*A(1,3)
     +	  - A(3,2)*A(2,3)*A(1,1)
     +	  - A(3,3)*A(2,1)*A(1,2)
C
C
      RETURN
      END SUBROUTINE MDET
C
C----------------------------------------------------------------------
C
      SUBROUTINE ONEM(A)
C
C      THIS SUBROUTINE STORES THE IDENTITY MATRIX IN THE
C      3 BY 3 MATRIX [A]
C
      IMPLICIT NONE
C
      INTEGER I,J
C
      DOUBLE PRECISION A(3,3)
C
C
      DO I=1,3
         DO J=1,3
	           IF (I .EQ. J) THEN
              A(I,J) = 1.D0
            ELSE
              A(I,J) = 0.D0
            END IF
         END DO
      END DO
C
C
      RETURN
      END SUBROUTINE ONEM
C
!****************************************************************************
