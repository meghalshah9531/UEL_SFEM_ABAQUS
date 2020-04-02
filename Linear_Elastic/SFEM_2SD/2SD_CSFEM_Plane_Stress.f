CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     USER ELEMENT SUBROUTINE (UEL)                                  C
C     FOR LINEAR ELASTICITY                                          C
C     WITH NEW METHOD OF SMOOTHED FINITE ELEMENT METHOD              C
C     AUTHOR: MEGHAL SHAH                                            C
C     INSTITUTION: TU BERGAKADEMIE FREIBERG                          C
C     PROJECT: PERSONAL PROGRAMMING PROJECT                          C
C     ELEMENT HAVING 2-DOMAINS                                       C
C                      N4 --x-- N6 --x-- N3                          C
C                       !       !        !                           C
C                       x  1SD  x  2SD   x                           C
C                       !       !        !                           C
C                      N1 --x-- N5 --x-- N2                          C
C      WHERE:   N1-N4 = Nodes of the QUAD element                    C
C               N5-N6 = Adding nodes to make 2-Domains of QUAD       C
C               x = Integration points at each edge                  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
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
C     ----------------------------------------------------------------
C     INFORMATION ABOUT VARIABLES USED IN UEL
C     CMATRIX : CONSTITUTIVE MATRIX 
C     E: YOUNG'S MODULUS OF THE MATERIAL 
C     NU: POISSION'S RATIO OF THE MATERIAL 
C     N1, N2, N3, N4 : FOUR NODES OF QUAD ELEMENT      
C     AREA_CELL: AREA OF THE CELL
C     ----------------------------------------------------------------
C     DECLARATION OF THE VARIABLES
      REAL(8) NU, E
      INTEGER I, J, K, M, N
      DOUBLE PRECISION, DIMENSION(2) :: N1
      DOUBLE PRECISION, DIMENSION(2) :: N2
      DOUBLE PRECISION, DIMENSION(2) :: N3
      DOUBLE PRECISION, DIMENSION(2) :: N4
      DOUBLE PRECISION, DIMENSION(2) :: N5
      DOUBLE PRECISION, DIMENSION(2) :: N6
      DOUBLE PRECISION, DIMENSION(7,4) :: SHAPE_V
      DOUBLE PRECISION, DIMENSION(8,3) :: TEMP
      DOUBLE PRECISION, DIMENSION(8,8) :: TEMP2
      DOUBLE PRECISION, DIMENSION(3,3) :: CMATRIX
      DOUBLE PRECISION, DIMENSION(2) :: AREA_CELL
      DOUBLE PRECISION, DIMENSION(12,4) :: B
      DOUBLE PRECISION, DIMENSION(8,6) :: BT
      DOUBLE PRECISION, DIMENSION(2,2) :: AKMAT 
C --------------------------------------------------------------------
C     SHAPE_V: LINEAR SHAPE FUNCTION VALUES AT INTEGRATION POITNTS 
C     B: STRAIN-DISPLACEMENT MATRIX WITH MODIFIED STRAIN FIELD 
C     BT: TRANSPOSE OF STRAIN-DISPLACEMENT
C     CMATRIX: MATRIAL MATRIX 
C     AKMAT: 2X2 MATRIX FOR ASSEMBLY OF K-MATRIX 
C --------------------------------------------------------------------
C     INITALIZATION OF MATRICES 
      CALL KASET2(AMATRX,8,8)
      CALL KASET2(B1,3,8)
      CALL KASET2(B2,3,8)
      CALL KASET2(TEMP,8,3)
      CALL KASET2(TEMP2,8,8)
      CALL KASET2(SHAPE_V,7,4)
      AREA_CELL(1) = 0.0
      AREA_CELL(2) = 0.0
C     PRINTING OUT ALL THE VARIABLES FOR THE ANALYSIS
      WRITE(7,*) '***********************************'
      WRITE(7,*) 'CURRENT ELEMENT NUMBER::::', JELEM
      WRITE(7,*) 'NUMBER OF NODES:::', NNODE
      WRITE(7,*) 'DOFs:::', MCRD
      WRITE(7,*) 'TOTAL DOF OF ELEMENT::: ', NDOFEL
      WRITE(7,*) '************************************'
C     EXTRACT COORDS OF FOUR NODES QUAD ELEMENT
      N1(1) = COORDS(1,1)
      N1(2) = COORDS(2,1)
      N2(1) = COORDS(1,2)
      N2(2) = COORDS(2,2)
      N3(1) = COORDS(1,3)
      N3(2) = COORDS(2,3)
      N4(1) = COORDS(1,4)
      N4(2) = COORDS(2,4)
C      WRITE(7,*) 'COORDINATES OF nODES :::'
C      WRITE(7,*) N1,N2,N3,N4
C     
C   COORDINATES OF THE ADDED POINTS  
      N5(1) = (COORDS(1,1)+COORDS(1,2))/2
      N5(2) = (COORDS(2,1)+COORDS(2,2))/2
      N6(1) = (COORDS(1,3)+COORDS(1,4))/2
      N6(2) = (COORDS(2,3)+COORDS(2,4))/2
C      WRITE(7,*) N5,N6  
C   VALUES OF SHAPE FUNCTIONs AT INTEGRATION POINTS
      CALL SHAPE(SHAPE_V)
C      PRINT *, "SHAPE VALUES AT INTEGRATION POINTS"
C      PRINT *, SHAPE_V      
      RETURN 
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     SUBROUITNES                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     SUBROUTINE FOR INITIALIZATION OF MATRIX WITH ZEROs
      SUBROUTINE KASET2(DMATRIX, IDIMX, IDIMY)
      INCLUDE 'ABA_PARAM.INC'
      PARAMETER (ZERO = 0.0D0)
      DIMENSION DMATRIX(IDIMX, IDIMY)
      DO I = 1, IDIMX
        DO J = 1, IDIMY
            DMATRIX(I,J) = ZERO
        END DO
      END DO
      RETURN
      END
C --------------------------------------------------------------------
C     FUNCTION TO STORE SHAPE FUNCTIONs VALUES
      SUBROUTINE K_SHAPE(SHAPE_V)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION SHAPE_V(7,4)
      SHAPE_V(1,1) = 0.75D0
      SHAPE_V(1,2) = 0.25D0
      SHAPE_V(1,3) = 0.0D0
      SHAPE_V(1,4) = 0.0D0
      SHAPE_V(2,1) = 0.25D0
      SHAPE_V(2,2) = 0.25D0
      SHAPE_V(2,3) = 0.25D0
      SHAPE_V(2,4) = 0.25D0
      SHAPE_V(3,1) = 0.0D0
      SHAPE_V(3,2) = 0.0D0
      SHAPE_V(3,3) = 0.25D0
      SHAPE_V(3,4) = 0.75D0
      SHAPE_V(4,1) = 0.5D0
      SHAPE_V(4,2) = 0.0D0
      SHAPE_V(4,3) = 0.0D0
      SHAPE_V(4,4) = 0.5D0
      SHAPE_V(5,1) = 0.25D0
      SHAPE_V(5,2) = 0.75D0
      SHAPE_V(5,3) = 0.0D0
      SHAPE_V(5,4) = 0.0D0
      SHAPE_V(6,1) = 0.0D0
      SHAPE_V(6,2) = 0.5D0
      SHAPE_V(6,3) = 0.5D0
      SHAPE_V(6,4) = 0.0D0
      SHAPE_V(7,1) = 0.0D0
      SHAPE_V(7,2) = 0.0D0
      SHAPE_V(7,3) = 0.75D0
      SHAPE_V(7,4) = 0.25D0
      RETURN
      END
C --------------------------------------------------------------------
