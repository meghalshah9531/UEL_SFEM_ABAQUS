CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     USER ELEMENT SUBROUTINE (UEL)                                  C
C     FOR LINEAR ELASTICITY                                          C
C     WITH NEW METHOD OF SMOOTHED FINITE ELEMENT METHOD              C
C     AUTHOR: MEGHAL SHAH                                            C
C     INSTITUTION: TU BERGAKADEMIE FREIBERG                          C
C     PROJECT: PERSONAL PROGRAMMING PROJECT                          C
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
      REAL(8) AREA_CELL, NU, E
      INTEGER I, J, K, M, N
      DOUBLE PRECISION, DIMENSION(2) :: N1
      DOUBLE PRECISION, DIMENSION(2) :: N2
      DOUBLE PRECISION, DIMENSION(2) :: N3
      DOUBLE PRECISION, DIMENSION(2) :: N4 
      DOUBLE PRECISION, DIMENSION(3,8) :: B
      DOUBLE PRECISION, DIMENSION(8,3) :: B_TRANSPOSE
      DOUBLE PRECISION, DIMENSION(4,4) :: GAUSVAL
      DOUBLE PRECISION, DIMENSION(8,3) :: TEMP
      DOUBLE PRECISION, DIMENSION(8,8) :: TEMP2
      DOUBLE PRECISION, DIMENSION(3,3) :: CMATRIX
      DOUBLE PRECISION, DIMENSION(3,1) :: STRESS
      DOUBLE PRECISION, DIMENSION(3,1) :: STRAIN
C     INITALIZATION OF MATRICES 
      CALL KASET2(AMATRX,8,8)
      CALL KASET2(B,3,8)
      CALL KASET2(TEMP,8,3)
      CALL KASET2(TEMP2,8,8)
      DO I = 1,NDOFEL
        DO J = 1,NRHS
            RHS(I,J) = 0.d0
        END DO 
      END DO
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
C     STORE OF INTEGRATION POINT VALUES       
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
C      write(7,*) 'gaussval'
C      WRITE(7,*) GAUSVAL
C     STRAIN-DISPLACEMENT MATRIX  
      CALL B_I(N1,N2,N3,N4,GAUSVAL(1,1),GAUSVAL(2,1),GAUSVAL(3,1),
     1     GAUSVAL(4,1),B(1:3,1:2),B_TRANSPOSE(1:2,1:3),AREA_CELL)
      CALL B_I(N1,N2,N3,N4,GAUSVAL(1,2),GAUSVAL(2,2),GAUSVAL(3,2),
     1     GAUSVAL(4,2),B(1:3,3:4),B_TRANSPOSE(3:4,1:3),AREA_CELL)
      CALL B_I(N1,N2,N3,N4,GAUSVAL(1,3),GAUSVAL(2,3),GAUSVAL(3,3),
     1     GAUSVAL(4,3),B(1:3,5:6),B_TRANSPOSE(5:6,1:3),AREA_CELL)
      CALL B_I(N1,N2,N3,N4,GAUSVAL(1,4),GAUSVAL(2,4),GAUSVAL(3,4),
     1     GAUSVAL(4,4),B(1:3,7:8),B_TRANSPOSE(7:8,1:3),AREA_CELL)
      PRINT *, 'AREA OF THE CELL :::', AREA_CELL
C      WRITE(7,*) 'STRAIN-DISPLACEMENT MATRIX ***********'
C      WRITE(7,*) B
C     MATERIAL PROPERTIES 
      E = PROPS(1)
      NU = PROPS(2)
C      WRITE(7,*) 'E', E 
C      WRITE(7,*) 'NU', NU 
C     PLANE STRAIN CONDITION 
C      CMATRIX(1,1)=E*(1.0-NU)/((1.0+NU)*(1.0-2*NU))
C      CMATRIX(1,2)=E*NU/((1.0+NU)*(1.0-2*NU))
C      CMATRIX(1,3)=0.0
C      CMATRIX(2,1)=E*NU/((1.0+NU)*(1.0-2*NU))
C      CMATRIX(2,2)=E*(1.0-NU)/((1.0+NU)*(1.0-2*NU))
C      CMATRIX(2,3)=0.0
C      CMATRIX(3,1)=0.0
C      CMATRIX(3,2)=0.0
C      CMATRIX(3,3)=(1.0-2*NU)/(2*(1.0+NU)*(1.0-2*NU))
C     PLANE STRESS ANALYSIS 
      CMATRIX(1,1) = E / (1-NU*NU)
      CMATRIX(1,2) = (E*NU)/(1-NU*NU)
      CMATRIX(1,3) = 0.0
      CMATRIX(2,1) = (E*NU)/(1-NU*NU)
      CMATRIX(2,2) = E / (1-NU*NU)
      CMATRIX(2,3) = 0.0
      CMATRIX(3,1) = 0.0
      CMATRIX(3,2) = 0.0
      CMATRIX(3,3) = (E*(1-NU))/(1-NU*NU)
C      WRITE(7,*) 'CMATRIX',CMATRIX
C
C     STIFFNESS MATRIX
      TEMP = MATMUL(B_TRANSPOSE,CMATRIX)
      WRITE(7,*) 'TEMP', TEMP
      TEMP2 = MATMUL(TEMP,B)
      TEMP2 = TEMP2*AREA_CELL
      DO I = 1,NDOFEL
        DO J = 1, NDOFEL
            AMATRX(I,J) = TEMP2(I,J)
        END DO 
      END DO 
      WRITE(7,*) 'AMATRIX:::'
      WRITE(7,*) AMATRX
C     SAVING STRESS,STRAIN,STRAIN-DISPLACEMENT AND AREA_CELL FOR OUTPUT
      DO I = 1,3
        DO J = 1,NDOFEL
            STRAIN(I,1) = B(I,J)*U(J)
        END DO 
      END DO 
      PRINT *, 'STRAIN', STRAIN
      DO I = 1,3
        DO J = 1,3
            STRESS(I,1) = CMATRIX(I,J)*STRAIN(I,1)
        END DO
      END DO 
      PRINT *, 'STRESS', STRESS 
C     SAVE SVARS
      SVARS(1) = STRAIN(1,1)
      SVARS(2) = STRAIN(2,1)
      SVARS(3) = STRAIN(3,1)
      SVARS(4) = STRESS(1,1)
      SVARS(5) = STRESS(2,1)
      SVARS(6) = STRESS(3,1)
      SVARS(7) = AREA_CELL
C     RESIDUAL VECTOR - RHS 
      DO I = 1,NDOFEL
        DO J = 1,NDOFEL
            RHS(I,1) = RHS(I,1)-AMATRX(I,J)*U(J)
        END DO 
      END DO 
C     WRITING OUTPUT ***** 
      WRITE(*,*) JELEM
      DO I = 1,7
        WRITE(7,*) 'SVARS',SVARS(I)
      END DO
      DO I = 1,8
            WRITE(7,*) 'u', U(I)
      END DO                      
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
C     ----------------------------------------------------------------
C     SUBROUTINE FOR STRAIN-DISPLACEMENT MATRIX 
      SUBROUTINE B_I(P,Q,R,S,INT_PT1,INT_PT2,INT_PT3,INT_PT4,Bi,BiT,
     1   AREA_CELL)
      INCLUDE 'ABA_PARAM.INC'
      REAL(8) INT_PT1,INT_PT2,INT_PT3,INT_PT4,AREA_CELL
      REAL(8) X(4),Y(4),L(4)
      REAL(8) NORMALS(2,4),Bi(3,2),BiT(2,3)
      REAL(8) P(2),Q(2),R(2),S(2)
      REAL(8) BIx,BIy
C       L: LENGTH OF FOUR BOUNDARIES 
C       X(4),Y(4): X & Y COORDINATES OF FOUR NODES 
C       NORMALS(4): OUTWARD NORMALS AT INTEGRATION PTs
C       INT_PT: INTEGRATION POINTs
C       AREA_CELL: AREA_CELL OF THE CELL 
        X(1) = P(1)
        Y(1) = P(2)
        X(2) = Q(1)
        Y(2) = Q(2)
        X(3) = R(1)
        Y(3) = R(2)
        X(4) = S(1)
        Y(4) = S(2)   
C       GET AREA_CELL OF THE CELL
        CALL AREA_QUAD(X,Y,AREA_CELL,L)
C        WRITE(7,*) '******AREA_CELL*******'
C        WRITE(7,*) AREA_CELL
C       NORMALs OUTWARD AT INTEGRATION POINTS
        CALL FIND_NORMALs(X,Y,L,NORMALS)
C        WRITE(7,*) '*******NORMALS******'
C        WRITE(7,*) NORMALS
C       STRAIN-DISPLACEMENT COMPONENT         
        BIx = (NORMALS(1,1)*INT_PT1*L(1) + NORMALS(1,2)*INT_PT2*L(2) +
     1         NORMALS(1,3)*INT_PT3*L(3) + NORMALS(1,4)*INT_PT4*L(4))
        BIx = BIx / AREA_CELL
        BIy = (NORMALS(2,1)*INT_PT1*L(1) + NORMALS(2,2)*INT_PT2*L(2) +
     1         NORMALS(2,3)*INT_PT3*L(3) + NORMALS(2,4)*INT_PT4*L(4))
        BIy = BIy / AREA_CELL
        CALL KASET2(Bi,3,2)
C        WRITE(7,*) 'KASET2-bi',Bi
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
C       TRANSPOSE OF STRAIN-DISPLACEMENT COMPONENT 
        CALL TRANSPOSE_MATRIX(Bi,BiT,3,2)
        RETURN
        END 
C       --------------------------------------------------------------
C       SUBROUTINE FOR CALCULATING AREA_CELL OF CELL
        SUBROUTINE AREA_QUAD(X,Y,AREA_CELL,LENGTH)
        INCLUDE 'ABA_PARAM.INC'
        REAL(8) X(4),Y(4),AREA_CELL,LENGTH(4)
        DO I = 1,3
            CALL EDGE_LENGTH(X(I),Y(I),X(I+1),Y(I+1),LENGTH(I))
        END DO 
        CALL EDGE_LENGTH(X(4),Y(4),X(1),Y(1),LENGTH(4))
        AREA_CELL = LENGTH(1)*LENGTH(2)
        RETURN
        END             
C       --------------------------------------------------------------
C       SUBROUTINE TO FIND ANGLE OF QUADRATURE 
C        SUBROUTINE ANGLE_QUAD(X,Y,ANGLE,EDGE,DIA)
C            INCLUDE 'ABA_PARAM.INC'
C            REAL(8) X(4),Y(4),ANGLE(4),EDGE(4),DIA(2)
C            REAL(8) coA1,coA2,coB1,coB2,coC1,coC2,coD1,coD2
C            CALL EDGE_LENGTH(X(1),Y(1),X(3),Y(3),DIA(1))
C            CALL EDGE_LENGTH(X(2),Y(2),X(4),Y(4),DIA(2))
CC               
C            DO I = 1,3
C                CALL EDGE_LENGTH(X(I),Y(I),X(I+1),Y(I+1),EDGE(I))
C            END DO 
C            CALL EDGE_LENGTH(X(4),Y(4),X(1),Y(1),EDGE(4))
C            coA1=(d1**2+edge(4)**2-edge(3)**2)/(2*d1*edge(4))
C            coA2=(d1**2+edge(1)**2-edge(2)**2)/(2*d1*edge(1))
C            coB1=(d2**2+edge(1)**2-edge(4)**2)/(2*d2*edge(1))
C            coB2=(d2**2+edge(2)**2-edge(3)**2)/(2*d2*edge(2))
C            coC1=(d1**2+edge(2)**2-edge(1)**2)/(2*d1*edge(2))
C            coC2=(d1**2+edge(3)**2-edge(4)**2)/(2*d1*edge(3))
C            coD1=(d2**2+edge(3)**2-edge(2)**2)/(2*d2*edge(3))
C            coD2=(d2**2+edge(4)**2-edge(1)**2)/(2*d2*edge(4))
C           
C            ANGLE(1) = acos(coA1)+acos(coA2)        
C            ANGLE(2) = acos(coB1)+acos(coB2)
C            ANGLE(3) = acos(coC1)+acos(coC2)
C            ANGLE(4) = acos(coD1)+acos(coD2)
C            RETURN
C        END
C       --------------------------------------------------------------
C       CALCULATE LENGTH OF EDGE OF CELL 
        SUBROUTINE EDGE_LENGTH(a,b,c,d,l)
        INCLUDE 'ABA_PARAM.INC'
        REAL(8) a,b,c,d,l 
        l = SQRT((a-c)**2 + (b-d)**2)
        RETURN
        END
C       --------------------------------------------------------------
C       FIND OUT NORMALs AT INTEGRATION POINTs
        SUBROUTINE FIND_NORMALs(X,Y,LENGTH,NORMALS)
        INCLUDE 'ABA_PARAM.INC'
        REAL(8) X(4),Y(4),LENGTH(4),NORMALS(2,4)
        DO I = 1,3 
            NORMALS(1,I)=0.5D0*(Y(I+1)-Y(I))/(LENGTH(I)/2)
            NORMALS(2,I)=-0.5D0*(X(I+1)-X(I))/(LENGTH(I)/2)
        END DO 
        NORMALS(1,4)=0.5D0*(Y(1)-Y(4))/(LENGTH(4)/2)
        NORMALS(2,4)=-0.5D0*(X(1)-X(4))/(LENGTH(4)/2)
        RETURN
        END
C       --------------------------------------------------------------
C       TRANSPOSE OF THE MATRIX 
        SUBROUTINE TRANSPOSE_MATRIX(MAT,TRANS_MAT,I,J)
        INCLUDE 'ABA_PARAM.INC'
        REAL(8) MAT(I,J), TRANS_MAT(J,I)
        DO M = 1,I
            DO N = 1,J
                TRANS_MAT(N,M) = MAT(M,N)
            END DO 
        END DO
        RETURN
        END