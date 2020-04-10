CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   USER ELEMENT SUBROUTINE (UEL)                                    C
C   --- LINEAR PIEZOELECTRICITY ---                                  C
C   SMOOTH FINITE ELEMENT METHOD (SFEM) -- 2-SMOOTHED DOMAINS        C
C   AUTHOR: MEGHAL SHAH                                              C
C   INSTITUTION: TU BERGAKADEMIE FREIBERG                            C
C   PROJECT: PERSONAL PROGRAMMING PROJECT                            C
C                      N4 --x-- N6 --x-- N3                          C
C                       !       !        !                           C
C                       x  1SD  x  2SD   x                           C
C                       !       !        !                           C
C                      N1 --x-- N5 --x-- N2                          C
C                                                                    C
C   N1-N4 = Nodes of Quadrilateral Element                           C
C   N5-N6 = Added Nodes to make 2-domains                            C
C   x = Integration Points at each edge of Element                   C
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
        DIMENSION CMATRIX(3,3), SHAPE_V(7,4), eMATRIX(2,3), 
     1  gMATRIX(2,2)
C   INFORMATION ABOUT VARIABLES 
        REAL :: c11,c13,c33,c55,e31,e33,e15,g11,g33
        INTEGER :: I,J,K
        DOUBLE PRECISION, DIMENSION(2) :: N1
        DOUBLE PRECISION, DIMENSION(2) :: N2
        DOUBLE PRECISION, DIMENSION(2) :: N3
        DOUBLE PRECISION, DIMENSION(2) :: N4
        DOUBLE PRECISION, DIMENSION(2) :: N5
        DOUBLE PRECISION, DIMENSION(2) :: N6
        DOUBLE PRECISION, DIMENSION(3,2) :: eTMATRIX
        DOUBLE PRECISION, DIMENSION(2,2) :: gTMATRIX
        DOUBLE PRECISION, DIMENSION(6,8) :: B_u
        DOUBLE PRECISION, DIMENSION(8,6) :: B_uT
        DOUBLE PRECISION, DIMENSION(4,4) :: Bphi
        DOUBLE PRECISION, DIMENSION(4,4) :: BphiT
        DOUBLE PRECISION, DIMENSION(8,8) :: K_u
        DOUBLE PRECISION, DIMENSION(8,4) :: K_uphi
        DOUBLE PRECISION, DIMENSION(4,8) :: K_phiu
        DOUBLE PRECISION, DIMENSION(4,4) :: K_phi
        DOUBLE PRECISION, DIMENSION(2,2) :: TEMP
        DOUBLE PRECISION, DIMENSION(2,1) :: TEMP2
        DOUBLE PRECISION, DIMENSION(1,1) :: TEMP3
        DOUBLE PRECISION, DIMENSION(2) :: AREA_CELL
C   INITIALIZATION OF MATRICES 
        CALL KASET2(SHAPE_V,7,4)
        CALL KASET2(TEMP,2,2)
        CALL KASET2(TEMP2,2,1)
        CALL KASET2(TEMP3,1,1)
        CALL KASET2(CMATRIX,3,3)
        CALL KASET2(eMATRIX,2,3)
        CALL KASET2(eTMATRIX,3,2)
        CALL KASET2(gMATRIX,2,2)
        CALL KASET2(gTMATRIX,2,2)
        CALL KASET2(AMATRX,12,12)
        CALL KASET2(B_u,6,8)
        CALL KASET2(Bphi,4,4)
        CALL KASET2(K_u,8,8)
        CALL KASET2(K_uphi,8,4)
        CALL KASET2(K_phi,4,4)
        AREA_CELL(1) = 0.0
        AREA_CELL(2) = 0.0
C   PRINTING OUT VARIABLES 
        WRITE(7,*) '***********************************'
        WRITE(7,*) 'CURRENT ELEMENT NUMBER::::', JELEM
        WRITE(7,*) 'NUMBER OF NODES:::', NNODE
        WRITE(7,*) 'DOFs:::', MCRD
        WRITE(7,*) 'TOTAL DOF OF ELEMENT::: ', NDOFEL
        WRITE(7,*) '************************************'
C   EXCTRACT COORDS OF FOUR NODES QUAD ELEMENT
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
C   VALUES OF SHAPE FUNCTION AT INTEGRATION POINTS 
        CALL K_SHAPE(SHAPE_V)
C   MATERIAL PROPERTIES 
        c11 = PROPS(1)
        c13 = PROPS(2)
        c33 = PROPS(3)
        c55 = PROPS(4)
        e31 = PROPS(5)
        e33 = PROPS(6)
        e15 = PROPS(7)
        g11 = PROPS(8)
        g33 = PROPS(9)
C   ELASTIC STIFFNESS MATRIX
        CALL ELASTIC_MATRIX(CMATRIX,c11,c13,c33,c55)
C   PIEZOELECTRIC MATRIX
        CALL PIEZOELECTRIC_MATRIX(eMATRIX,e31,e33,e15)
c      WRITE(7,*) 'eMATRIX', eMATRIX
        CALL TRANSPOSE_MATRIX(eMATRIX,eTMATRIX,2,3)
C   DIELECTRIC CONSTANT MATRIX 
        CALL DIELECTRIC_MATRIX(gMATRIX,g11,g33)
c      WRITE(7,*) 'gMATRIX', gMATRIX
        CALL TRANSPOSE_MATRIX(gMATRIX,gTMATRIX,2,2)
C   STRAIN-DISPLACEMENT MATRICES 
        CALL STRAIN_MECH(N1,N2,N3,N4,N5,N6,SHAPE_V,B_u,B_uT,AREA_CELL)
C   STRIAN-ELECTRICAL POTENTIAL MATRICES 
        CALL STRAIN_ELE(N1,N2,N3,N4,N5,N6,SHAPE_V,Bphi,BphiT,
     1   AREA_CELL)
        RETURN 
        END SUBROUTINE UEL 
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
C   Function to store Shape Function Values at Integration points 
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
        SUBROUTINE ELASTIC_MATRIX(CMATRIX,c11,c13,c33,c55)
        INCLUDE 'ABA_PARAM.INC'
        DOUBLE PRECISION, DIMENSION(3,3) :: CMATRIX
        REAL(8) c11,c13,c33,c55
        CMATRIX(1,1) = c11
        CMATRIX(1,2) = c13
        CMATRIX(1,3) = 0.0
        CMATRIX(2,1) = c13
        CMATRIX(2,2) = c33
        CMATRIX(2,3) = 0.0
        CMATRIX(3,1) = 0.0
        CMATRIX(3,2) = 0.0
        CMATRIX(3,3) = c55
        RETURN
        END
C --------------------------------------------------------------------
C   Transpose Of the Matrix 
        SUBROUTINE TRANSPOSE_MATRIX(MAT,TRANS_MAT,I,J)
        INCLUDE 'ABA_PARAM.INC'
        REAL(8) MAT(I,J), TRANS_MAT(J,I)
        INTEGER I,J,M,N
        DO M = 1,I
            DO N = 1,J
                TRANS_MAT(N,M) = MAT(M,N)
            END DO 
        END DO 
        RETURN
        END
C --------------------------------------------------------------------
C   Dielectric Constant Matrix 
        SUBROUTINE DIELECTRIC_MATRIX(gMATRIX,g11,g33)
        INCLUDE 'ABA_PARAM.INC'
        REAL(8) g11,g33
        DOUBLE PRECISION, DIMENSION(2,2) :: gMATRIX
        gMATRIX(1,1) = g11
        gMATRIX(1,2) = 0.0
        gMATRIX(2,1) = 0.0
        gMATRIX(2,2) = g33
        RETURN
        END
C --------------------------------------------------------------------
C   Strain-Displacement Matrix Calculation
        SUBROUTINE STRAIN_MECH(N1,N2,N3,N4,N5,N6,SHAPE_V,B,BT,AREA_CELL)
        INCLUDE 'ABA_PARAM.INC'
        DOUBLE PRECISION, DIMENSION(2) :: N1
        DOUBLE PRECISION, DIMENSION(2) :: N2
        DOUBLE PRECISION, DIMENSION(2) :: N3
        DOUBLE PRECISION, DIMENSION(2) :: N4
        DOUBLE PRECISION, DIMENSION(2) :: N5
        DOUBLE PRECISION, DIMENSION(2) :: N6
        DOUBLE PRECISION, DIMENSION(2) :: AREA_CELL
        DOUBLE PRECISION, DIMENSION(7,4) :: SHAPE_V
        DOUBLE PRECISION, DIMENSION(6,8) :: B_u
        DOUBLE PRECISION, DIMENSION(8,6) :: B_uT
C   @ NODE 1
C             
        CALL B_ui(N1,N5,N6,N4,SHAPE_V(1,1),SHAPE_V(2,1),
     1  SHAPE_V(3,1),SHAPE_V(4,1),B_u(1:3,1:2),B_uT(1:2,1:3),
     2  AREA_CELL(1))
        CALL B_ui(N5,N2,N3,N6,SHAPE_V(5,1),SHAPE_V(6,1),
     1  SHAPE_V(7,1),SHAPE_V(2,1),B_u(4:6,1:2),B_uT(1:2,4:6),
     2  AREA_CELL(2))
C   @ NODE 2
        CALL B_ui(N1,N5,N6,N4,SHAPE_V(1,2),SHAPE_V(2,2),
     1  SHAPE_V(3,2),SHAPE_V(4,2),B_u(1:3,3:4),B_uT(3:4,1:3),
     2  AREA_CELL(1))
        CALL B_ui(N5,N2,N3,N6,SHAPE_V(5,2),SHAPE_V(6,2),
     1  SHAPE_V(7,2),SHAPE_V(2,2),B_u(4:6,3:4),B_uT(3:4,4:6),
     2  AREA_CELL(2))
C   @ NODE 3
        CALL B_ui(N1,N5,N6,N4,SHAPE_V(1,3),SHAPE_V(2,3),
     1  SHAPE_V(3,3),SHAPE_V(4,3),B_u(1:3,5:6),B_uT(5:6,1:3),
     2  AREA_CELL(1))
        CALL B_ui(N5,N2,N3,N6,SHAPE_V(5,3),SHAPE_V(6,3),
     1  SHAPE_V(7,3),SHAPE_V(2,3),B_u(4:6,5:6),B_uT(5:6,4:6),
     2  AREA_CELL(2))
C   NODE 4
        CALL B_ui(N1,N5,N6,N4,SHAPE_V(1,4),SHAPE_V(2,4),
     1  SHAPE_V(3,4),SHAPE_V(4,4),B_u(1:3,7:8),B_uT(7:8,1:3),
     2  AREA_CELL(1))
        CALL B_ui(N5,N2,N3,N6,SHAPE_V(5,4),SHAPE_V(6,4),
     1  SHAPE_V(7,4),SHAPE_V(2,4),B_u(4:6,7:8),B_uT(7:8,4:6),
     2  AREA_CELL(2))
        RETURN
        END
C --------------------------------------------------------------------
C   STRAIN-DISPLACEMENT COMPONENT FOR EACH NODE
        SUBROUTINE B_ui(P,Q,R,S,INT_PT1,INT_PT2,INT_PT3,INT_PT4,
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
C        WRITE(7,*) '******AREA_CELL*******'
C        WRITE(7,*) AREA_CELL
C   Calculate Normals at Integration Points 
        NORMALS(1,1) = 0.5D0 * (Y(2)-Y(1))/(L(1)/2)
        NORMALS(2,1) = 0.5D0 * (X(2)-X(1))/(L(1)/2)
        NORMALS(1,2) = 0.5D0 * (Y(3)-Y(2))/(L(2)/2)
        NORMALS(2,2) = 0.5D0 * (X(3)-X(2))/(L(2)/2)
        NORMALS(1,3) = 0.5D0 * (Y(4)-Y(3))/(L(3)/2)
        NORMALS(2,3) = 0.5D0 * (X(4)-X(3))/(L(3)/2)
        NORMALS(1,4) = 0.5D0 * (Y(1)-Y(4))/(L(4)/2)
        NORMALS(2,4) = 0.5D0 * (X(1)-X(4))/(L(4)/2)
C        WRITE(7,*) '*******NORMALS******'
C        WRITE(7,*) NORMALS      
C   COMPONENT OF STRAIN-DISPLACEMENT FOR NODE      
        BIx = (NORMALS(1,1)*INT_PT1*L(1) + NORMALS(1,2)*INT_PT2*L(2) +
     1   NORMALS(1,3)*INT_PT3*L(3) + NORMALS(1,4)*INT_PT4*L(4))
        BIx = BIx / AREA_CELL
        BIy = (NORMALS(2,1)*INT_PT1*L(1) + NORMALS(2,2)*INT_PT2*L(2) +
     1       NORMALS(2,3)*INT_PT3*L(3) + NORMALS(2,4)*INT_PT4*L(4))
        BIy = BIy / AREA_CELL
        CALL KASET2(Bi,3,2)
        Bi(1,1) = Bi(1,1) + BIx
        Bi(2,2) = Bi(2,2) + BIy
        Bi(3,1) = Bi(3,1) + BIy
        Bi(3,2) = Bi(3,2) + BIx
C        WRITE(7,*) 'AFTER REPLACEMENT'
C        WRITE(7,*) Bi
C        WRITE(7,*) 'bi(1,1)',Bi(1,1)
C        WRITE(7,*) 'BI(1,2)',Bi(1,2)
C        WRITE(7,*) 'bi(2,1)',Bi(2,1)
C        WRITE(7,*) 'Bi(2,2)',Bi(2,2)
C        WRITE(7,*) 'Bi', Bi      
        CALL TRANSPOSE_MATRIX(Bi,BiT,3,2)
        RETURN
        END
C --------------------------------------------------------------------
C   Find Area of the element for Regular and Irregular Meshes 
         SUBROUTINE AREA_QUAD(X,Y,AREA_CELL,LENGTH)
        INCLUDE 'ABA_PARAM.INC'
        REAL(8) X(4),Y(4),AREA_CELL,AREA1,AREA2,LENGTH(4)
        REAL(8) ANGLE(4), DIA(2)
C   Length of all edges       
        LENGTH(1) = SQRT((X(1)-X(2))**2 + (Y(1)-Y(2))**2)
        LENGTH(2) = SQRT((X(2)-X(3))**2 + (Y(2)-Y(3))**2)
        LENGTH(3) = SQRT((X(3)-X(4))**2 + (Y(3)-Y(4))**2)
        LENGTH(4) = SQRT((X(4)-X(1))**2 + (Y(4)-Y(1))**2)
C   length of Diagonal
        DIA(1) = SQRT((X(1)-X(3))**2 + (Y(1)-Y(3))**2)
        DIA(2) = SQRT((X(2)-X(4))**2 + (Y(2)-Y(4))**2)
C   
        CALL ANGLE_QUAD(X,Y,ANGLE,LENGTH,DIA)
        IF (ANGLE(1) >= 3.14 .OR. ANGLE(3) >= 3.14) THEN 
        AREA1 = 0.5*(X(1)*Y(2)+X(2)*Y(3)+X(3)*Y(1)-
     1        X(2)*Y(1)-X(3)*Y(2)-X(1)*Y(3))
        AREA2 = 0.5*(X(1)*Y(3)+X(3)*Y(4)+X(4)*Y(1)-
     1        X(1)*Y(4)-X(3)*Y(1)-X(4)*Y(3))
        AREA_CELL = ABS(AREA1) + ABS(AREA2)
        ELSE IF (ANGLE(2) >= 3.14 .OR. ANGLE(4) >= 3.14) THEN
        AREA1 = 0.5*(X(1)*Y(2)+X(2)*Y(4)+X(4)*Y(1)-
     1        X(1)*Y(4)-X(2)*Y(1)-X(4)*Y(2))
        AREA2 = 0.5*(X(2)*Y(3)+X(3)*Y(4)+X(4)*Y(2)-
     1        X(2)*Y(4)-X(3)*Y(2)-X(4)*Y(3))
        AREA_CELL = ABS(AREA1) + ABS(AREA2)
        ELSE IF (ANGLE(1)<3.14 .AND. ANGLE(2)<3.14 .AND. ANGLE(3)<3.14
     1          .AND. ANGLE(4)<3.14 .AND. DIA(1)<=DIA(2)) THEN
        AREA1 = 0.5*(X(1)*Y(2)+X(2)*Y(3)+X(3)*Y(1)-
     1        X(1)*Y(3)-X(2)*Y(1)-X(3)*Y(2))
        AREA2 = 0.5*(X(1)*Y(3)+X(3)*Y(4)+X(4)*Y(1)-
     1        X(1)*Y(4)-X(3)*Y(1)-X(4)*Y(3))
        AREA_CELL = ABS(AREA1) + ABS(AREA2)                        
        ELSE IF (ANGLE(1)<3.14 .AND. ANGLE(2)<3.14 .AND. ANGLE(3)<3.14
     1        .AND. ANGLE(4)<3.14 .AND. DIA(1)>DIA(2)) THEN
        AREA1 = 0.5*(X(1)*Y(2)+X(2)*Y(4)+X(4)*Y(1)-
     1        X(1)*Y(4)-X(2)*Y(1)-X(4)*Y(2))
        AREA2 = 0.5*(X(2)*Y(3)+X(3)*Y(4)+X(4)*Y(2)-
     1        X(2)*Y(4)-X(3)*Y(2)-X(4)*Y(3))
        AREA_CELL = ABS(AREA1) + ABS(AREA2)
        END IF
        RETURN
        END
C --------------------------------------------------------------------
C   ANGLE OF QUAD
        SUBROUTINE ANGLE_QUAD(X,Y,ANGLE,EDGE,DIA)
        INCLUDE 'ABA_PARAM.INC'
        REAL(8) X(4),Y(4),ANGLE(4),EDGE(4),DIA(2)
        REAL(8) ADJ_A1,ADJ_A2,ADJ_B1,ADJ_B2,ADJ_C1,ADJ_C2,ADJ_D1,ADJ_D2
C               
        ADJ_A1=(DIA(1)**2+edge(4)**2-edge(3)**2)/(2*DIA(1)*edge(4))
        ADJ_A2=(DIA(1)**2+edge(1)**2-edge(2)**2)/(2*DIA(1)*edge(1))
        ADJ_B1=(DIA(2)**2+edge(1)**2-edge(4)**2)/(2*DIA(2)*edge(1))
        ADJ_B2=(DIA(2)**2+edge(2)**2-edge(3)**2)/(2*DIA(2)*edge(2))
        ADJ_C1=(DIA(1)**2+edge(2)**2-edge(1)**2)/(2*DIA(1)*edge(2))
        ADJ_C2=(DIA(1)**2+edge(3)**2-edge(4)**2)/(2*DIA(1)*edge(3))
        ADJ_D1=(DIA(2)**2+edge(3)**2-edge(2)**2)/(2*DIA(2)*edge(3))
        ADJ_D2=(DIA(2)**2+edge(4)**2-edge(1)**2)/(2*DIA(2)*edge(4))
C           
        ANGLE(1) = acos(ADJ_A1)+acos(ADJ_A2)        
        ANGLE(2) = acos(ADJ_B1)+acos(ADJ_B2)
        ANGLE(3) = acos(ADJ_C1)+acos(ADJ_C2)
        ANGLE(4) = acos(ADJ_D1)+acos(ADJ_D2)
C        PRINT *, ANGLE
        RETURN
        END                   
