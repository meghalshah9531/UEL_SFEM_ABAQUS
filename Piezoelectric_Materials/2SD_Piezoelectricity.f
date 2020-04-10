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
C                             
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