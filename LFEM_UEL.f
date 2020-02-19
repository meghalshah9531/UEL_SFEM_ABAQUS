C     ----------------------------------------------------------------
C     USER ELEMENT SUBROUTINE 
C     FOR LINEAR STATIC ANALYSIS UNDER PLANE STRESS CONDITION
C     ELEMENT: BILINIEAR QAUD ELEMENT
C     SHAPE FUNCITON: BILINIEAR SHAPE FUNCTION
C     BY: MEGHAL SHAH 
C     AT: TU BERGAKADEMIE FREIBERG 
C     ----------------------------------------------------------------
      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1 PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2 KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,NPREDF,
     3 LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,PERIOD)
C
      INCLUDE 'ABA_PARAM.INC'
	PARAMETER (ZERO = 0.D0, ONE = 1.D0, HALF = 0.5D0)
C
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     1 SVARS(*),ENERGY(8),COORDS(MCRD,NNODE),U(NDOFEL),
     2 DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),PARAMS(*),
     3 JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     4 PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)
	DIMENSION COORDINATES(4,2)
C     DECLARATION OF THE VARIABLES 
	INTEGER :: I,J
	REAL :: CONST_C, E, NU, POINT1, POINT2 
	DOUBLE PRECISION, DIMENSION(3,3) :: CMATRIX
	DOUBLE PRECISION, DIMENSION(2,4) :: GAUSSPTS
	DOUBLE PRECISION, DIMENSION(2,4) :: DSHAPE
	
C     INITIALIZATION OF MATRICES 
	IF (NRHS.EQ.1) THEN  ! RHS VECTOR 
		DO I = 1,MLVARX
			RHS(I,1) = ZERO
		END DO 
		DO I = 1, NDOFEL
			DO J = 1, NDOFEL
				AMATRX(I,J) = ZERO
			END DO 
          END DO 
      END IF 
C     FOR PLANE STRESS LINEAR STATIC ANALYSIS
      E = PROPS(1)
	NU = PROPS(2)
	CONST_C = E/(1-NU*NU)
	CMATRIX(1,1) = ONE 
	CMATRIX(1,2) = NU 
	CMATRIX(1,3) = ZERO
	CMATRIX(2,1) = NU
	CMATRIX(2,2) = ONE 
	CMATRIX(2,3) = ZERO 
	CMATRIX(3,1) = ZERO 
	CMATRIX(3,2) = ZERO
	CMATRIX(3,3) = HALF*(1-NU)
	CMATRIX = CMATRIX*CONST_C
C     COORDINATES (TRANSPOSE OF THE COORDS)
	COORDINATES(1,1) = COORDS(1,1)
	COORDINATES(1,2) = COORDS(2,1)
      COORDINATES(2,1) = COORDS(1,2)
      COORDINATES(2,2) = COORDS(2,2)
      COORDINATES(3,1) = COORDS(1,3)
      COORDINATES(3,2) = COORDS(2,3)
      COORDINATES(4,1) = COORDS(1,4)
      COORDINATES(4,2) = COORDS(2,4)
C     GAUSS POINTS 
	GAUSSPTS(1,1) = -0.5773
	GAUSSPTS(2,1) = -0.5773
	GAUSSPTS(1,2) = 0.5773
	GAUSSPTS(2,2) = -0.5773
	GAUSSPTS(1,3) = 0.5773
	GAUSSPTS(2,3) = 0.5773
	GAUSSPTS(1,4) = -0.5773
	GAUSSPTS(2,4) = 0.5773 
C     INTEGRATION OVER THE GAUSS POINTS 
	DO IP = 1,4
		POINT1 = GAUSSPTS(1,IP)
		POINT2 = GAUSSPTS(2,IP)
C             DETERMINATION OF THE SHAPE FUNCTIONS 
		DSHAPE(1,1) = -0.25 * (1-POINT2)
          DSHAPE(1,2) = 0.25 * (1-POINT2)
          DSHAPE(1,3) = 0.25 * (1+POINT2)
          DSHAPE(1,4) = -0.25 * (1+POINT2)
          DSHAPE(2,1) = -0.25 * (1-POINT1)
          DSHAPE(2,2) = -0.25 * (1+POINT1)
          DSHAPE(2,3) = 0.25 * (1+POINT1)
          DSHAPE(2,4) = 0.25 * (1-POINT1)
      RETURN
      END