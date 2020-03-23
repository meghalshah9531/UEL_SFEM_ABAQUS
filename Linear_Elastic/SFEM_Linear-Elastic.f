CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C           USER ELEMENT SUBROUTINE (UEL)                            C
C           LINEAR ELASTIC MATERIAL                                  C
C           SMOOTHED FINITE ELEMENT METHOD WITH 4-NODED QUAD ELE     C
C           AUTHOR: MEGHAL SHAH                                      C
C           INTITUTION: TU BERGAKADEMIE FREIBERG                     C
C           PROJECT: PERSONAL PROGRAMMING PROJECT                    C
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
C
C   Variables ...
C   E: Young's Modulus
C   ENU: Poission's Ratio
C   N1, N2, N3, N4: Four Nodes of QUAD element
C   AREA_CELL: Area of the Cell 
C   
C   Declaration of the Variables
      DOUBLE PRECISION, DIMENSION(2) :: N1
      DOUBLE PRECISION, DIMENSION(2) :: N2
      DOUBLE PRECISION, DIMENSION(2) :: N3
      DOUBLE PRECISION, DIMENSION(2) :: N4
      DOUBLE PRECISION, DIMENSION(3,8) :: B_MATRIX
      DOUBLE PRECISION, DIMENSION(8,3) :: B_TRANSPOSE
      DOUBLE PRECISION, DIMENSION(8,3) :: TEMP
      DOUBLE PRECISION, DIMENSION(8,8) :: TEMP2
      DOUBLE PRECISION, DIMENSION(3,3) :: C_MATRIX
      DOUBLE PRECISION, DIMENSION(3,1) :: STRESS
      DOUBLE PRECISION, DIMENSION(3,1) :: STRAIN
      DOUBLE PRECISION, DIMENSION(4,4) :: GAUSVAL
C   Initiallization of Variables 
      CALL KZEROs(AMATRX,8,8)
      CALL KZEROs(B_MATRIX,3,8)
      CALL KZEROs(B_TRANSPOSE,8,3)
      CALL KZEROs(TEMP,8,3)
      CALL KZEROs(TEMP2,8,8)
C   Writing out Varibles for Check 
      WRITE(7,*) '*******************************'
      WRITE(7,*) 'ELEMENT No:', JELEM
      WRITE(7,*) 'NNODEs:', NNODE
      WRITE(7,*) 'DOF:', MCRD
      WRITE(7,*) 'TOTAL DOFs:', NDOFEL
      WRITE(7,*) '*******************************'
C   Extract COORDs of Four Nodes 
      N1(1) = COORDS(1,1)
      N1(2) = COORDS(2,1)
      N2(1) = COORDS(1,2)
      N2(2) = COORDS(2,2)
      N3(1) = COORDS(1,3)
      N3(2) = COORDS(2,3) 
      N4(1) = COORDS(1,4)
      N4(2) = COORDS(2,4)
C   Storing Values of Shape functions @ Integration Points
      GAUSVAL(1,1) = 0.5D0
      GAUSVAL(1,2) = 0.5D0
      GAUSVAL(1,3) = 0.0D0
      GAUSVAL(1,4) = 0.0D0
      GAUSVAL(2,1) = 0.0D0
      GAUSVAL(2,2) = 0.5D0
      GAUSVAL(2,3) = 0.5D0
      GAUSVAL(2,4) = 0.0D0
      GAUSVAL(3,1) = 0.0D0
      GAUSVAL(3,2) = 0.0D0
      GAUSVAL(3,3) = 0.5D0
      GAUSVAL(3,4) = 0.5D0
      GAUSVAL(4,1) = 0.5D0
      GAUSVAL(4,2) = 0.0D0
      GAUSVAL(4,3) = 0.0D0
      GAUSVAL(4,4) = 0.5D0
C   Material Properties 
      E = PROPS(1)
      ENU = PROPS(2)
C   Plane Strain Condition
      C_MATRIX(1,1) = E / (1-NU*NU)
      C_MATRIX(1,2) = (E*NU)/(1-NU*NU)
      C_MATRIX(1,3) = 0.0
      C_MATRIX(2,1) = (E*NU)/(1-NU*NU)
      C_MATRIX(2,2) = E / (1-NU*NU)
      C_MATRIX(2,3) = 0.0
      C_MATRIX(3,1) = 0.0
      C_MATRIX(3,2) = 0.0
      C_MATRIX(3,3) = (E*(1-NU))/(1-NU*NU)
      WRITE(7,*) 'C_MATRIX', C_MATRIX
C   Strain-Displacement Matrix 
      

















        RETURN 
        END
C   ----------------
C   Subroutines 
C   ----------------
      SUBROUTINE KZEROs(DMAT, IX, IY)
      INCLUDE 'ABA_PARAM.INC'
      PARAMETER (ZERO = 0.0D0)
      DIMENSION DMAT(IX, IY)
      DO I = 1, IX 
        DO J = 1, IY
            DMAT(I,J) = ZERO
        END DO 
      END DO 
      RETURN
      END