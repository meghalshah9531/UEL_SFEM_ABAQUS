C ------------------------------------------------------------
C     ABAQUS USER ELEMENT ROUTINE 
C     FOR: 4-Noded Qaudrilateral Element 
C          Linear Static Analysis 
C     Author: MEGHAL SHAH
C     INSTITUTION: TU BERGAKADEMIE FREIBERG
C ------------------------------------------------------------
      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1 PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2 KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,NPREDF,
     3 LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,PERIOD)
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     1 SVARS(*),ENERGY(8),COORDS(MCRD,NNODE),U(NDOFEL),
     2 DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),PARAMS(*),
     3 JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     4 PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)
C     
C     VARIABLES INITIALIZATION
      INTEGER :: I,J,K
      REAL :: POINT1,POINT2,E,NU,CONST_C
      DOUBLE PRECISION, DIMENSION(2,2) :: JACOBIAN
      DOUBLE PRECISION, DIMENSION(3,3) :: CMATRIX
      DOUBLE PRECISION, DIMENSION(2,4) :: GAUSSPTS
      DOUBLE PRECISION, DIMENSION(2,4) :: DSHAPE
      DOUBLE PRECISION :: DET_J 
      DOUBLE PRECISION, DIMENSION(2,2) :: INV_JACOBIAN
      DOUBLE PRECISION, DIMENSION(2) :: A1
      DOUBLE PRECISION, DIMENSION(2) :: A2
      DOUBLE PRECISION, DIMENSION(2) :: A3
      DOUBLE PRECISION, DIMENSION(2) :: A4
      DOUBLE PRECISION, DIMENSION(3,8) :: BMATRIX 
      DOUBLE PRECISION, DIMENSION(8,3) :: BTRANSPOSE
      DOUBLE PRECISION, DIMENSION(3,8) :: TEMP 
      DOUBLE PRECISION, DIMENSION(3) :: STRESS
      DOUBLE PRECISION, DIMENSION(3) :: STRAIN
C
C     ZERO INITIALIZATION 
      DO I = 1,NDOFEL
          DO J = 1,NDOFEL
              AMATRX(I,J) = 0.0
          END DO
      END DO 
C
      DO I = 1,NDOFEL
          RHS(I,1) = 0.0
      END DO 
C 
C     INTEGER PROPERTIES 
      E = PROPS(1)
      NU = PROPS(2)
C           
C     FOR PLANE STRESS ANALYSIS 
      CONST_C = E / (1-NU*NU)
      CMATRIX(1,1) = 1.0
      CMATRIX(1,2) = NU
      CMATRIX(1,3) = 0.0
      CMATRIX(2,1) = NU
      CMATRIX(2,2) = 1.0
      CMATRIX(2,3) = 0.0
      CMATRIX(3,1) = 0.0
      CMATRIX(3,2) = 0.0
      CMATRIX(3,3) = 0.5*(1-NU)
      CMATRIX = CMATRIX * CONST_C
C
C     GAUSS POINTS:
      GAUSSPTS(1,1) = -0.5773
      GAUSSPTS(2,1) = -0.5773
      GAUSSPTS(1,2) = 0.5773
      GAUSSPTS(2,2) = -0.5773
      GAUSSPTS(1,3) = 0.5773
      GAUSSPTS(2,3) = 0.5773
      GAUSSPTS(1,4) = -0.5773
      GAUSSPTS(2,4) = 0.5773
C           
C     GAUSS POINT INTEGRATION
      DO IP = 1,4
          POINT1 = GAUSSPTS(1,IP)
          POINT2 = GAUSSPTS(2,IP)
C     DETERMINATION OF DERIVATIVES OF THE SHAPE FUNCTIONS
          DSHAPE(1,1) = -0.25 * (1 - POINT2)
          DSHAPE(1,2) =  0.25 * (1 - POINT2)
          DSHAPE(1,3) =  0.25 * (1 + POINT2)
          DSHAPE(1,4) = -0.25 * (1 + POINT2)
          DSHAPE(2,1) = -0.25 * (1 - POINT1)
          DSHAPE(2,2) = -0.25 * (1 + POINT1)
          DSHAPE(2,3) =  0.25 * (1 + POINT1)
          DSHAPE(2,4) =  0.25 * (1 - POINT1)
C     DETERMINE JACOBIAN OF THE TRANSFORMATION
          COORDINATES = TRANSPOSE(COORDS)
          JACOBIAN = MATMUL(DSHAPE,COORDINATES)
C     DETERMINANT OF THE JACOBIAN MATRIX 
          DET_J = ABS(JACOBIAN(1,1)*JACOBIAN(2,2) - JACOBIAN(1,2)
     1                *JACOBIAN(2,1))
C     INVERSE OF THE JACOBIAN 
          INV_JACOBIAN(1,1) = JACOBIAN(2,2)
          INV_JACOBIAN(1,2) = -JACOBIAN(1,2)
          INV_JACOBIAN(2,1) = -JACOBIAN(2,1)
          INV_JACOBIAN(2,2) = JACOBIAN(1,1)
          INV_JACOBIAN = INV_JACOBIAN/DET_J
C     B_MATRIX
          A1 = MATMUL(INV_JACOBIAN,DSHAPE(:,1))
          A2 = MATMUL(INV_JACOBIAN,DSHAPE(:,2))
          A3 = MATMUL(INV_JACOBIAN,DSHAPE(:,3))
          A4 = MATMUL(INV_JACOBIAN,DSHAPE(:,4))
C
          BMATRIX(1,1) = A1(1)
          BMATRIX(1,2) = 0.0
          BMATRIX(1,3) = A2(1)
          BMATRIX(1,4) = 0.0
          BMATRIX(1,5) = A3(1)
          BMATRIX(1,6) = 0.0
          BMATRIX(1,7) = A4(1)
          BMATRIX(1,8) = 0.0
C
          BMATRIX(2,1) = 0.0
          BMATRIX(2,2) = A1(2)
          BMATRIX(2,3) = 0.0
          BMATRIX(2,4) = A2(2)
          BMATRIX(2,5) = 0.0
          BMATRIX(2,6) = A3(2)
          BMATRIX(2,7) = 0.0
          BMATRIX(2,8) = A4(2)
C
          BMATRIX(3,1) = A1(2)
          BMATRIX(3,2) = A1(1)
          BMATRIX(3,3) = A2(2)
          BMATRIX(3,4) = A2(1)
          BMATRIX(3,5) = A3(2)
          BMATRIX(3,6) = A3(1)
          BMATRIX(3,7) = A4(2)
          BMATRIX(3,8) = A4(1)
C
          BTRANSPOSE = TRANSPOSE(BMATRIX)
C     AMATRIX 
          DO I=1,3
              DO J = 1,8
                  TEMP(I,J) = 0.0
              END DO 
          END DO
          TEMP = MATMUL(CMATRIX,BMATRIX)
          AMATRX = AMATRX + DET_J *MATMUL(BTRANSPOSE,TEMP)
          DO I = 1,3
              DO J = 1,8
                  TEMP(I,J) = 0.0
              END DO 
          END DO 
C     CALCULATION OF STRESS AND STRAIN 
          STRAIN = MATMUL(BMATRIX,U)
          STRESS = MATMUL(CMATRIX,STRAIN)
C     RHS MATRIX
          DO I = 1,NDOFEL
              DO J = 1,3
                  RHS(I,1) = RHS(I,1) - BMATRIX(J,I)*STRESS(J)*DET_J
              END DO
          END DO       
C     SVARS 
          SVARS(1) = STRAIN(1)
          SVARS(2) = STRAIN(2) 
          SVARS(3) = STRAIN(3)
          SVARS(4) = STRESS(1)
          SVARS(5) = STRESS(2) 
          SVARS(6) = STRESS(3)
          SVARS(7) = PROPS(2)*(STRESS(1)+STRESS(2))
C     PRINTING 
          WRITE(6,*) '************SVARS*****************'
          DO I = 1,7
              WRITE(6,*) SVARS(I)
          END DO
      END DO
      RETURN 
      END 