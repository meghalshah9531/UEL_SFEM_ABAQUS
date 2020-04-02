CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     USER ELEMENT SUBROUTINE (UEL)                                  C
C     --- LINEAR ELASTICITY ---                                      C
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
C    INFORMATION ABOUT VARIABLES USED IN UEL
C    CMATRIX : CONSTITUTIVE MATRIX 
C    E: YOUNG'S MODULUS OF THE MATERIAL 
C    NU: POISSION'S RATIO OF THE MATERIAL 
C    N1, N2, N3, N4 : FOUR NODES OF QUAD ELEMENT      
C    AREA_CELL: AREA OF THE CELL
C    SHAPE_V: LINEAR SHAPE FUNCTION VALUES AT INTEGRATION POITNTS 
C    B: STRAIN-DISPLACEMENT MATRIX WITH MODIFIED STRAIN FIELD 
C    BT: TRANSPOSE OF STRAIN-DISPLACEMENT
C    CMATRIX: MATRIAL MATRIX 
C    KMAT: 2X2 MATRIX FOR ASSEMBLY OF K-MATRIX 
C     ----------------------------------------------------------------
C    DECLARATION OF THE VARIABLES
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
C    INITALIZATION OF MATRICES 
      CALL KASET2(AMATRX,8,8)
      CALL KASET2(B1,3,8)
      CALL KASET2(B2,3,8)
      CALL KASET2(TEMP,8,3)
      CALL KASET2(TEMP2,8,8)
      CALL KASET2(SHAPE_V,7,4)
      AREA_CELL(1) = 0.0
      AREA_CELL(2) = 0.0
C    PRINTING OUT ALL THE VARIABLES FOR THE ANALYSIS
      WRITE(7,*) '***********************************'
      WRITE(7,*) 'CURRENT ELEMENT NUMBER::::', JELEM
      WRITE(7,*) 'NUMBER OF NODES:::', NNODE
      WRITE(7,*) 'DOFs:::', MCRD
      WRITE(7,*) 'TOTAL DOF OF ELEMENT::: ', NDOFEL
      WRITE(7,*) '************************************'
C    EXTRACT COORDS OF FOUR NODES QUAD ELEMENT
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
C   MATERIAL MATRIX FOR PLANE STRESS ANALYSIS 
C   MATERIAL PROPERTIES
      E = PROPS(1)
      NU = PROPS(2)
      CALL MATERIAL_MATRIX(CMATRIX, E, NU)
      RETURN 
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     SUBROUITNES                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   Initialization of Matrix with Zeros
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
C   Function to store Shape function Values at Integration Points 
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
C   Material Matrix 
      SUBROUTINE MATERIAL_MATRIX(CMATRIX,E,NU)
      INCLUDE 'ABA_PARAM.INC'
      DOUBLE PRECISION, DIMENSION(3,3) :: CMATRIX
      REAL(8) E, NU
      CMATRIX(1,1) = E / (1-NU*NU)
      CMATRIX(1,2) = (E*NU)/(1-NU*NU)
      CMATRIX(1,3) = 0.0
      CMATRIX(2,1) = (E*NU)/(1-NU*NU)
      CMATRIX(2,2) = E / (1-NU*NU)
      CMATRIX(2,3) = 0.0
      CMATRIX(3,1) = 0.0
      CMATRIX(3,2) = 0.0
      CMATRIX(3,3) = (E*(1-NU))/(1-NU*NU)
      RETURN
      END
C --------------------------------------------------------------------
C   Strain-Displacement Matrix Calculation
      SUBROUTINE STRAIN_D(N1,N2,N3,N4,N5,N6,SHAPE_V,B,BT,AREA_CELL)
      INCLUDE 'ABA_PARAM.INC'
      DOUBLE PRECISION, DIMENSION(2) :: N1
      DOUBLE PRECISION, DIMENSION(2) :: N2
      DOUBLE PRECISION, DIMENSION(2) :: N3
      DOUBLE PRECISION, DIMENSION(2) :: N4
      DOUBLE PRECISION, DIMENSION(2) :: N5
      DOUBLE PRECISION, DIMENSION(2) :: N6
      DOUBLE PRECISION, DIMENSION(2) :: AREA_CELL
      DOUBLE PRECISION, DIMENSION(7,4) :: SHAPE_V
      DOUBLE PRECISION, DIMENSION(12,4) :: B
      DOUBLE PRECISION, DIMENSION(8,6) :: BT
C   @ NODE 1
C             
      CALL B_COMPONENT(N1,N5,N6,N4,SHAPE_V(1,1),SHAPE_V(2,1),
     1  SHAPE_V(3,1),SHAPE_V(4,1),B(1:3,1:2),BT(1:2,1:3),AREA_CELL(1))
      CALL B_COMPONENT(N5,N2,N3,N6,SHAPE_V(5,1),SHAPE_V(6,1),
     1  SHAPE_V(7,1),SHAPE_V(2,1),B(1:3,3:4),BT(1:2,4:6),AREA_CELL(2))
C   @ NODE 2
      CALL B_COMPONENT(N1,N5,N6,N4,SHAPE_V(1,2),SHAPE_V(2,2),
     1  SHAPE_V(3,2),SHAPE_V(4,2),B(4:6,1:2),BT(3:4,1:3),AREA_CELL(1))
      CALL B_COMPONENT(N5,N2,N3,N6,SHAPE_V(5,2),SHAPE_V(6,2),
     1  SHAPE_V(7,2),SHAPE_V(2,2),B(4:6,3:4),BT(3:4,4:6),AREA_CELL(2))
C
      CALL B_COMPONENT(N1,N5,N6,N4,SHAPE_V(1,3),SHAPE_V(2,3),
     1  SHAPE_V(3,3),SHAPE_V(4,3),B(7:9,1:2),BT(5:6,1:3),AREA_CELL(1))
      CALL B_COMPONENT(N5,N2,N3,N6,SHAPE_V(5,3),SHAPE_V(6,3),
     1  SHAPE_V(7,3),SHAPE_V(2,3),B(7:9,3:4),BT(5:6,4:6),AREA_CELL(2))
C   
      CALL B_COMPONENT(N1,N5,N6,N4,SHAPE_V(1,4),SHAPE_V(2,4),
     1SHAPE_V(3,4),SHAPE_V(4,4),B(10:12,1:2),BT(7:8,1:3),AREA_CELL(1))
      CALL B_COMPONENT(N5,N2,N3,N6,SHAPE_V(5,4),SHAPE_V(6,4),
     1SHAPE_V(7,4),SHAPE_V(2,4),B(10:12,3:4),BT(7:8,4:6),AREA_CELL(2))
      RETURN
      END 
C --------------------------------------------------------------------
C   STRAIN-DISPLACEMENT COMPONENT FOR EACH NODE
      SUBROUTINE B_COMPONENT(P,Q,R,S,INT_PT1,INT_PT2,INT_PT3,INT_PT4,
     1   Bi,BiT,AREA_CELL)
      INCLUDE 'ABA_PARAM.INC'
      REAL(8) INT_PT1,INT_PT2,INT_PT3,INT_PT4,AREA_CELL
      REAL(8) X(4),Y(4),L(4)
      REAL(8) NORMALS(2,4),Bi(3,2),BiT(2,3)
      REAL(8) P(2),Q(2),R(2),S(2)
      REAL(8) BIx,BIy
C   L: LENGTH OF FOUR BOUNDARIES 
C   X(4),Y(4): X & Y COORDINATES OF FOUR NODES 
C   NORMALS(4): OUTWARD NORMALS AT INTEGRATION PTs
C   INT_PT: INTEGRATION POINTs
C   AREA_CELL: AREA_CELL OF THE CELL
      X(1) = P(1)
      Y(1) = P(2)
      X(2) = Q(1)
      Y(2) = Q(2)
      X(3) = R(1)
      Y(3) = R(2)
      X(4) = S(1)
      Y(4) = S(2)
C   Area of the cell 
      CALL AREA_QUAD(X,Y,AREA_CELL,L)
C   Calculate Normals at Integration Points 
      NORMALS(1,1) = 0.5D0 * (Y(2)-Y(1))/(L(1)/2)
      NORMALS(2,1) = 0.5D0 * (X(2)-X(1))/(L(1)/2)
      NORMALS(1,2) = 0.5D0 * (Y(3)-Y(2))/(L(2)/2)
      NORMALS(2,2) = 0.5D0 * (X(3)-X(2))/(L(2)/2)
      NORMALS(1,3) = 0.5D0 * (Y(4)-Y(3))/(L(3)/2)
      NORMALS(2,3) = 0.5D0 * (X(4)-X(3))/(L(3)/2)
      NORMALS(1,4) = 0.5D0 * (Y(1)-Y(4))/(L(4)/2)
      NORMALS(2,4) = 0.5D0 * (X(1)-X(4))/(L(4)/2)
C   COMPONENT OF STRAIN-DISPLACEMENT FOR NODE      
      BIx = (NORMALS(1,1)*INT_PT1*L(1) + NORMALS(1,2)*INT_PT2*L(2) +
     1       NORMALS(1,3)*INT_PT3*L(3) + NORMALS(1,4)*INT_PT4*L(4))
      BIx = BIx / AREA_CELL
      BIy = (NORMALS(2,1)*INT_PT1*L(1) + NORMALS(2,2)*INT_PT2*L(2) +
     1       NORMALS(2,3)*INT_PT3*L(3) + NORMALS(2,4)*INT_PT4*L(4))
      BIy = BIy / AREA_CELL
      CALL KASET2(Bi,3,2)
      Bi(1,1) = Bi(1,1) + BIx
      Bi(2,2) = Bi(2,2) + BIy
      Bi(3,1) = Bi(3,1) + BIy
      Bi(3,2) = Bi(3,2) + BIx
      CALL TRANSPOSE_MATRIX(Bi,BiT,3,2)
      RETURN
      END