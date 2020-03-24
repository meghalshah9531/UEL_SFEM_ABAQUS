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
      CALL B_I(N1,N2,N3,N4,GAUSVAL(1,1),GAUSVAL(1,2),GAUSVAL(1,3),
     1  GAUSVAL(1,4), B_MATRIX(1:3,1:2), B_TRANSPOSE(1:2,1:3),AREA_CELL)
      CALL B_I(N1,N2,N3,N4,GAUSVAL(2,1),GAUSVAL(2,2),GAUSVAL(2,3),
     1  GAUSVAL(2,4), B_MATRIX(1:3,3:4), B_TRANSPOSE(3:4,1:3),AREA_CELL)
      CALL B_I(N1,N2,N3,N4,GAUSVAL(3,1),GAUSVAL(3,2),GAUSVAL(3,3),
     1  GAUSVAL(3,4), B_MATRIX(1:3,5:6), B_TRANSPOSE(5:6,1:3),AREA_CELL)
      CALL B_I(N1,N2,N3,N4,GAUSVAL(4,1),GAUSVAL(4,2),GAUSVAL(4,3),
     1  GAUSVAL(4,4), B_MATRIX(1:3,7:8), B_TRANSPOSE(7:8,1:3),AREA_CELL)
      WRITE(7,*) 'Strain-Displacement Matrix'
      WRITE(7,*) B_MATRIX
C      
C   Stiffness Matrix 
      TEMP = MATMUL(B_TRANSPOSE, C_MATRIX)
      WRITE(7,*) 'TEMP', TEMP
      TEMP2 = MATMUL(TEMP, B_MATRIX)
      TEMP2 = TEMP2*AREA_CELL
      DO I = 1, NDOFEL
        DO J = 1, NDOFEL
            AMATRX(I,J) = TEMP2(I,J)
        END DO 
      END DO 
      WRITE(7,*) 'AMATRIX'
      WRITE(7,*) AMATRX
C   
C     STRAIN :::
      DO I = 1,3
        DO J = 1, NDOFEL
            STRAIN(I,1) = B_MATRIX(I,J)*U(J)
        END DO 
      END DO
C   
C     STRESS :::
      DO I = 1,3
        DO J = 1,3
            STRESS(I,1) = C_MATRIX(I,J)*STRAIN(I,1)
        END DO 
      END DO 
C   
C     RHS :::
      DO I = 1,NDOFEL
        DO J = 1,NRHS
            RHS(I,1) = RHS(I,1)-AMATRX(I,J)*U(J)
        END DO 
      END DO
C   
      SVARS(1) = STRAIN(1,1)
      SVARS(2) = STRAIN(2,1)
      SVARS(3) = STRAIN(3,1)
      SVARS(4) = STRESS(1,1)
      SVARS(5) = STRESS(2,1)
      SVARS(6) = STRESS(3,1)
      SVARS(7) = AREA_CELL
C   
      DO I = 1,7 
        WRITE(7,*) 'sVARS', SVARS(I)
      END DO 
      RETURN 
      END
C   ----------------
C   Subroutines 
C   ----------------
C   Initiallization of Matrices with Zeros
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
C   Strain-Displacement Matrix
      SUBROUTINE B_I(N1,N2,N3,N4,INT_1,INT_2,INT_3,INT_4,Bi,BiT,AREA)
      INCLUDE 'ABA_PARAM.INC'
      REAL :: INT_1, INT_2, INT_3, INT_4, AREA
      DOUBLE PRECISION, DIMENSION(4) :: X
      DOUBLE PRECISION, DIMENSION(4) :: Y
      DOUBLE PRECISION, DIMENSION(4) :: L
      DOUBLE PRECISION, DIMENSION(2,4) :: NORMALS(2,4)
      DOUBLE PRECISION, DIMENSION(3,2) :: Bi
      DOUBLE PRECISION, DIMENSION(2,3) :: BiT
      REAL :: BIx, BIy 
      X(1) = N1(1)
      Y(1) = N1(2) 
      X(2) = N2(1)
      Y(2) = N2(2)
      X(3) = N3(1)
      Y(3) = N3(2)
      X(4) = N4(1)
      Y(4) = N4(2)	
C       Area of the cell 
      CALL AREA_QUAD(X,Y,AREA,L)
      WRITE(7,*) '****AREA****'
      WRITE(7,*) AREA
C       Normals
      CALL FIND_NORMALs(X,Y,L,NORMALS)
      BIx = (NORMALS(1,1)*INT_1*L(1) + NORMALS(1,2)*INT_2*L(2) +
     1       NORMALS(1,3)*INT_3*L(3) + NORMALS(1,4)*INT_4*L(4))/AREA
      BIy = (NORMALS(2,1)*INT_1*L(1) + NORMALS(2,2)*INT_2*L(2) +
     1       NORMALS(2,3)*INT_3*L(3) + NORMALS(2,4)*INT_4*L(4))/AREA
      CALL KZEROs(Bi,3,2)
      Bi(1,1) = Bi(1,1) + BIx
      Bi(2,2) = Bi(2,2) + BIy
      Bi(3,1) = Bi(3,1) + BIy
      Bi(3,2) = Bi(3,2) + BIx
      WRITE(7,*) 'After Replacement'
      WRITE(7,*) Bi
      CALL TRANSPOSE_MAT(Bi,BiT,3,2)
      RETURN
      END
C   Area of the cell 
      SUBROUTINE AREA_QUAD(X,Y,AREA,L)
      INCLUDE 'ABA_PARAM.INC'
      REAL(8) X(4),Y(4),AREA,L(4)
      DO I = 1,3
        CALL EDGE_LENGTH(X(I),Y(I),X(I+1),Y(I+1),L(1))
      END DO 
      CALL EDGE_LENGTH(X(4),Y(4),X(1),Y(1),L(4))
      AREA = L(1)*L(2)
      RETURN
      END
C   Length of the edge of each cell 
      SUBROUTINE EDGE_LENGTH(a,b,c,d,l)
      INCLUDE 'ABA_PARAM.INC'
      REAL(8) a,b,c,d,l 
      l = SQRT((a-c)**2 + (b-d)**2)
      RETURN
      END
C   Outward Normals @ Integration points
      SUBROUTINE FIND_NORMALs(X,Y,L,NORMALS)
      INCLUDE 'ABA_PARAM.INC'
      REAL(8) X(4), Y(4), L(4), NORMALS(2,4)
      DO I = 1,3
        NORMALS(1,I)=0.5D0*(Y(I+1)-Y(I))/(LENGTH(I)/2)
        NORMALS(2,I)=-0.5D0*(X(I+1)-X(I))/(LENGTH(I)/2)
      END DO
      NORMALS(1,4)=0.5D0*(Y(1)-Y(4))/(LENGTH(4)/2)
      NORMALS(2,4)=-0.5D0*(X(1)-X(4))/(LENGTH(4)/2)
      RETURN
      END
C   Transpose Matrix function
      SUBROUTINE TRANSPOSE_MAT(DMAT,D_TRANS,I,J)
      INCLUDE 'ABA_PARAM.INC'
      REAL(8) DMAT(I,J), D_TRANS(J,I)
      DO II = 1,I
        DO JJ = 1,J
            D_TRANS(JJ,II) = DMAT(I,J)
        END DO
      END DO
      RETURN
      END
            