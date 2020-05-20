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
        DIMENSION SHAPE_V(7,4),CMATRIX(3,3),eMATRIX(2,3),
     1  eTMATRIX(3,2),gMATRIX(2,2),gTMATRIX(2,2), AREA_CELL(2),
     2 U_MECH(8), U_ELE(4), U_NEW(12)
C     
      REAL(8) c11,c13,c33,c55,e31,e33,e15,g11,g33
      INTEGER I,J,K,M,N
      DOUBLE PRECISION, DIMENSION(2) :: N1
      DOUBLE PRECISION, DIMENSION(2) :: N2
      DOUBLE PRECISION, DIMENSION(2) :: N3
      DOUBLE PRECISION, DIMENSION(2) :: N4
      DOUBLE PRECISION, DIMENSION(2) :: N5
      DOUBLE PRECISION, DIMENSION(2) :: N6
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
      DOUBLE PRECISION, DIMENSION(3,8) :: Bu_SD1
      DOUBLE PRECISION, DIMENSION(3,8) :: Bu_SD2
      DOUBLE PRECISION, DIMENSION(3,1) :: STRAIN_SD_MECH
      DOUBLE PRECISION, DIMENSION(3,1) :: STRESS_SD_MECH
      DOUBLE PRECISION, DIMENSION(2,4) :: Bphi_1SD
      DOUBLE PRECISION, DIMENSION(2,4) :: Bphi_2SD
      DOUBLE PRECISION, DIMENSION(2,4) :: Bphi_SD
      DOUBLE PRECISION, DIMENSION(2,1) :: ELE_FIELD_SD
      DOUBLE PRECISION, DIMENSION(2,1) :: D1
      DOUBLE PRECISION, DIMENSION(2,1) :: D2
      DOUBLE PRECISION, DIMENSION(2,1) :: ELE_DISP
      DOUBLE PRECISION, DIMENSION(3,8) :: Bu_SD
C   INITIALIZATION OF MATRICES 
      CALL KASET2(SHAPE_V,7,4)
      CALL KASET2(CMATRIX,3,3)
      CALL KASET2(TEMP,2,2)
      CALL KASET2(eMATRIX,2,3)
      CALL KASET2(gMATRIX,2,2)
      CALL KASET2(B_u,6,8)
      CALL KASET2(Bphi,4,4)
      CALL KASET2(K_u,8,8)
      CALL KASET2(K_uphi,8,4)
      CALL KASET2(K_phi,4,4)
      CALL KASET2(TEMP2,2,1)
      CALL KASET2(AMATRX,12,12)
      CALL KASET2(STRAIN_SD_MECH,3,1)
      CALL KASET2(STRESS_SD_MECH,3,1)
      CALL KASET2(ELE_FIELD_SD,2,1)
      CALL KASET2(D1,2,1)
      CALL KASET2(D2,2,1)
      AREA_CELL(1) = 0.0
      AREA_CELL(2) = 0.0
C   PRINTING OUT ALL THE VARIABLES
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
C      WRITE(7,*) N5,N6
C   VALUES OF SHAPE FUNCTIONs AT INTEGRATION POINTS
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
C      PRINT *, 'CMATRIX'
C      PRINT *, CMATRIX
C   PIEZOELECTRIC MATRIX
      CALL PIEZOELECTRIC_MATRIX(eMATRIX,e31,e33,e15)
C      PRINT *, eMATRIX
      CALL TRANSPOSE_MATRIX(eMATRIX,eTMATRIX,2,3)
C      PRINT *, eTMATRIX
C   DIELECTRIC CONSTANT MATRIX 
      CALL DIELECTRIC_MATRIX(gMATRIX,g11,g33)
C      PRINT *, gMATRIX
      CALL TRANSPOSE_MATRIX(gMATRIX,gTMATRIX,2,2)
C      WRITE(7,*) gTMATRIX
C   STRAIN-DISPLACEMENT MATRICES
      CALL STRAIN_MECH(N1,N2,N3,N4,N5,N6,SHAPE_V,B_u,B_uT,AREA_CELL)
C      PRINT *, 'B_U**********'
C      PRINT *, B_u
C   STRIAN-ELECTRICAL POTENTIAL MATRICES 
      CALL STRAIN_ELE(N1,N2,N3,N4,N5,N6,SHAPE_V,Bphi,BphiT,AREA_CELL)
C   Kuu (MECHANICAL)-MATRIX CALCULATION 
      DO I = 1,4 ! <-- NO. OF NODES IN QUAD ELEMENT
        DO J = 1,4 ! <-- NO. OF NODES IN QUAD ELEMENT
            DO K = 1,2 ! <-- REPRESENTS NO. OF CELLS 
                CALL KMECH(B_uT((2*I-1):2*I,(3*K-2):3*K),CMATRIX,
     1  B_u((3*K-2):3*K,(2*J-1):2*J),AREA_CELL(K),TEMP)
                K_u((2*I-1):2*I,(2*J-1):2*J) = K_u((2*I-1):2*I,
     1  (2*J-1):2*J) + TEMP
            END DO 
        END DO 
      END DO 
C   Kuphi (MECH-ELE)- MATRIX 
      DO I = 1,4
        DO J = 1,4
            DO K = 1,2
                CALL KMECH_ELE(B_uT((2*I-1):2*I,(3*K-2):3*K),eTMATRIX,
     1  Bphi((2*K-1):2*K,J:J),AREA_CELL(K),TEMP2)
                K_uphi((2*I-1):2*I,J:J) = K_uphi((2*I-1):2*I,J:J)
     1          +TEMP2
            END DO 
        END DO 
      END DO
C   Kphiu 
      CALL TRANSPOSE_MATRIX(K_uphi,K_phiu,8,4)
C   Kphi
      DO I = 1,4
        DO J = 1,4
            DO K = 1,2
                CALL KELE(BphiT(I,(2*K-1):2*K),gMATRIX,
     1           Bphi((2*K-1):2*K,J),AREA_CELL(K),TEMP3)
                K_phi(I,J) = K_phi(I,J) + TEMP3(1,1)
            END DO 
        END DO 
      END DO
C   K - ASSEMBLY 
      AMATRX(1,1) = K_u(1,1)
      AMATRX(1,2) = K_u(1,2)
      AMATRX(1,3) = K_uphi(1,1)
      AMATRX(1,4) = K_u(1,3)
      AMATRX(1,5) = K_u(1,4)
      AMATRX(1,6) = K_uphi(1,2)
      AMATRX(1,7) = K_u(1,5)
      AMATRX(1,8) = K_u(1,6)
      AMATRX(1,9) = K_uphi(1,3)
      AMATRX(1,10) = K_u(1,7)
      AMATRX(1,11) = K_u(1,8)
      AMATRX(1,12) = K_uphi(1,4)
C
      AMATRX(2,1) = K_u(2,1)
      AMATRX(2,2) = K_u(2,2)
      AMATRX(2,3) = K_uphi(2,1)
      AMATRX(2,4) = K_u(2,3)
      AMATRX(2,5) = K_u(2,4)
      AMATRX(2,6) = K_uphi(2,2)
      AMATRX(2,7) = K_u(2,5)
      AMATRX(2,8) = K_u(2,6)
      AMATRX(2,9) = K_uphi(2,3)
      AMATRX(2,10) = K_u(2,7)
      AMATRX(2,11) = K_u(2,8)
      AMATRX(2,12) = K_uphi(2,4)
C   
      AMATRX(3,1) = K_phiu(1,1)
      AMATRX(3,2) = K_phiu(1,2)
      AMATRX(3,3) = K_phi(1,1)
      AMATRX(3,4) = K_phiu(1,3)
      AMATRX(3,5) = K_phiu(1,4)
      AMATRX(3,6) = K_phi(1,2)
      AMATRX(3,7) = K_phiu(1,5)
      AMATRX(3,8) = K_phiu(1,6)
      AMATRX(3,9) = K_phi(1,3)
      AMATRX(3,10) = K_phiu(1,7)
      AMATRX(3,11) = K_phiu(1,8)
      AMATRX(3,12) = K_phi(1,4)
C                  
      AMATRX(4,1) = K_u(3,1)
      AMATRX(4,2) = K_u(3,2)
      AMATRX(4,3) = K_uphi(3,1)
      AMATRX(4,4) = K_u(3,3)
      AMATRX(4,5) = K_u(3,4)
      AMATRX(4,6) = K_uphi(3,2)
      AMATRX(4,7) = K_u(3,5)
      AMATRX(4,8) = K_u(3,6)
      AMATRX(4,9) = K_uphi(3,3)
      AMATRX(4,10) = K_u(3,7)
      AMATRX(4,11) = K_u(3,8)
      AMATRX(4,12) = K_uphi(3,4)
C
      AMATRX(5,1) = K_u(4,1)
      AMATRX(5,2) = K_u(4,2)
      AMATRX(5,3) = K_uphi(4,1)
      AMATRX(5,4) = K_u(4,3)
      AMATRX(5,5) = K_u(4,4)
      AMATRX(5,6) = K_uphi(4,2)
      AMATRX(5,7) = K_u(4,5)
      AMATRX(5,8) = K_u(4,6)
      AMATRX(5,9) = K_uphi(4,3)
      AMATRX(5,10) = K_u(4,7)
      AMATRX(5,11) = K_u(4,8)
      AMATRX(5,12) = K_uphi(4,4)
C   
      AMATRX(6,1) = K_phiu(2,1)
      AMATRX(6,2) = K_phiu(2,2)
      AMATRX(6,3) = K_phi(2,1)
      AMATRX(6,4) = K_phiu(2,3)
      AMATRX(6,5) = K_phiu(2,4)
      AMATRX(6,6) = K_phi(2,2)
      AMATRX(6,7) = K_phiu(2,5)
      AMATRX(6,8) = K_phiu(2,6)
      AMATRX(6,9) = K_phi(2,3)
      AMATRX(6,10) = K_phiu(2,7)
      AMATRX(6,11) = K_phiu(2,8)
      AMATRX(6,12) = K_phi(2,4)
C                  
      AMATRX(7,1) = K_u(5,1)
      AMATRX(7,2) = K_u(5,2)
      AMATRX(7,3) = K_uphi(5,1)
      AMATRX(7,4) = K_u(5,3)
      AMATRX(7,5) = K_u(5,4)
      AMATRX(7,6) = K_uphi(5,2)
      AMATRX(7,7) = K_u(5,5)
      AMATRX(7,8) = K_u(5,6)
      AMATRX(7,9) = K_uphi(5,3)
      AMATRX(7,10) = K_u(5,7)
      AMATRX(7,11) = K_u(5,8)
      AMATRX(7,12) = K_uphi(5,4)
C
      AMATRX(8,1) = K_u(6,1)
      AMATRX(8,2) = K_u(6,2)
      AMATRX(8,3) = K_uphi(6,1)
      AMATRX(8,4) = K_u(6,3)
      AMATRX(8,5) = K_u(6,4)
      AMATRX(8,6) = K_uphi(6,2)
      AMATRX(8,7) = K_u(6,5)
      AMATRX(8,8) = K_u(6,6)
      AMATRX(8,9) = K_uphi(6,3)
      AMATRX(8,10) = K_u(6,7)
      AMATRX(8,11) = K_u(6,8)
      AMATRX(8,12) = K_uphi(6,4)
C   
      AMATRX(9,1) = K_phiu(3,1)
      AMATRX(9,2) = K_phiu(3,2)
      AMATRX(9,3) = K_phi(3,1)
      AMATRX(9,4) = K_phiu(3,3)
      AMATRX(9,5) = K_phiu(3,4)
      AMATRX(9,6) = K_phi(3,2)
      AMATRX(9,7) = K_phiu(3,5)
      AMATRX(9,8) = K_phiu(3,6)
      AMATRX(9,9) = K_phi(3,3)
      AMATRX(9,10) = K_phiu(3,7)
      AMATRX(9,11) = K_phiu(3,8)
      AMATRX(9,12) = K_phi(3,4)
C
      AMATRX(10,1) = K_u(7,1)
      AMATRX(10,2) = K_u(7,2)
      AMATRX(10,3) = K_uphi(7,1)
      AMATRX(10,4) = K_u(7,3)
      AMATRX(10,5) = K_u(7,4)
      AMATRX(10,6) = K_uphi(7,2)
      AMATRX(10,7) = K_u(7,5)
      AMATRX(10,8) = K_u(7,6)
      AMATRX(10,9) = K_uphi(7,3)
      AMATRX(10,10) = K_u(7,7)
      AMATRX(10,11) = K_u(7,8)
      AMATRX(10,12) = K_uphi(7,4)
C
      AMATRX(11,1) = K_u(8,1)
      AMATRX(11,2) = K_u(8,2)
      AMATRX(11,3) = K_uphi(8,1)
      AMATRX(11,4) = K_u(8,3)
      AMATRX(11,5) = K_u(8,4)
      AMATRX(11,6) = K_uphi(8,2)
      AMATRX(11,7) = K_u(8,5)
      AMATRX(11,8) = K_u(8,6)
      AMATRX(11,9) = K_uphi(8,3)
      AMATRX(11,10) = K_u(8,7)
      AMATRX(11,11) = K_u(8,8)
      AMATRX(11,12) = K_uphi(8,4)
C   
      AMATRX(12,1) = K_phiu(4,1)
      AMATRX(12,2) = K_phiu(4,2)
      AMATRX(12,3) = K_phi(4,1)
      AMATRX(12,4) = K_phiu(4,3)
      AMATRX(12,5) = K_phiu(4,4)
      AMATRX(12,6) = K_phi(4,2)
      AMATRX(12,7) = K_phiu(4,5)
      AMATRX(12,8) = K_phiu(4,6)
      AMATRX(12,9) = K_phi(4,3)
      AMATRX(12,10) = K_phiu(4,7)
      AMATRX(12,11) = K_phiu(4,8)
      AMATRX(12,12) = K_phi(4,4) 
C   MECHANICAL DISPLACEMENT 
      U_MECH(1) = U(1)
      U_MECH(2) = U(2)
      U_MECH(3) = U(4)
      U_MECH(4) = U(5)
      U_MECH(5) = U(7)
      U_MECH(6) = U(8)
      U_MECH(7) = U(10)
      U_MECH(8) = U(11)
C   ELECTRICAL DISPLACEMENT
      U_ELE(1) = U(3)
      U_ELE(2) = U(6)
      U_ELE(3) = U(9)
      U_ELE(4) = U(12)
C   STRAIN AND STRESS CALCULATION
      Bu_SD1 = B_u(1:3,1:8)
      Bu_SD2 = B_u(4:6,1:8)
      Bu_SD = Bu_SD1 + Bu_SD2
c      DO I = 1,3 
c        DO J = 1,8 
c            STRAIN_SD_MECH(I,1) = STRAIN_SD_MECH(I,1) + 
c     1      Bu_SD(I,J)*U_MECH(J)
c        END DO 
c      END DO 
C       STRAIN FIELD OF SD1       
      DO I = 1,3 
        DO J = 1,8
            STRAIN_SD_MECH(I,1)=STRAIN_SD_MECH(I,1)+
     1  0.5*Bu_SD1(I,J)*U_MECH(J)
        END DO 
      END DO
C       STRAIN FIELD OF SD2
      DO M = 1,3
        DO N = 1,8
            STRAIN_SD_MECH(M,1)=STRAIN_SD_MECH(M,1)+
     1  0.5*Bu_SD2(M,N)*U_MECH(N)
        END DO
      END DO
C   STRESS FIELD 
c      DO I = 1,3 
c        DO J = 1,3
c            STRESS_SD_MECH(I,J)=STRESS_SD_MECH(I,J)+
c     1  CMATRIX(I,J)*STRAIN_SD_MECH(J,1)
c        END DO
c      END DO
C   ELECTRICAL FIELD CALCULATION 
      Bphi_1SD = Bphi(1:2,1:4)
      Bphi_2SD = Bphi(3:4,1:4)
      Bphi_SD = Bphi_1SD + Bphi_2SD 
      DO I = 1,2
        DO J = 1,4 
            ELE_FIELD_SD(I,1) = ELE_FIELD_SD(I,1) - 
     1          0.5*Bphi_SD(I,J) * U_ELE(J)
        END DO 
      END DO 
      STRESS_SD_MECH = MATMUL(CMATRIX,STRAIN_SD_MECH) - 
     1 MATMUL(eTMATRIX,ELE_FIELD_SD)
C      ELE_FIELD_SD = -MATMUL(Bphi_SD,U_ELE)
C   ELE. FIELD OF SD1 
c      DO I = 1,2 
c        DO J = 1,4 
c            ELE_FIELD_SD(I,1) = ELE_FIELD_SD(I,1)-
c     1  0.5*Bphi_1SD(I,J)*U_ELE(J)
c        END DO
c      END DO
C   ELE. FIELD OF SD2 c
c      DO K = 1,2 
c        DO M = 1,4 
c            ELE_FIELD_SD(I,1) = ELE_FIELD_SD(I,1)-
c     1  0.5*Bphi_2SD(K,M)*U_ELE(M)
c        END DO
c      END DO
c      ELE_FIELD_SD = -ELE_FIELD_SD
C   DIELECTRIC DISPLACEMENT VECTOR 
C   D = (eMATRIX) * (STRAIN) + (gMATRIX) * (ELE_FIELD)
      D1 = MATMUL(eMATRIX, STRAIN_SD_MECH)
c      DO I = 1,2 
c        DO J = 1,3
c            D1(I,1) =eMATRIX(I,J)*STRAIN_SD_MECH(J,1)
c        END DO
c      END DO
C   
      D2 = MATMUL(gMATRIX, ELE_FIELD_SD)
c      DO I = 1,2
c        DO J = 1,2
c            D2(I,1)= D2(I,1) +gMATRIX(I,J)*ELE_FIELD_SD(J,1)
c       END DO 
c      END DO
C    
      ELE_DISP = D1 + D2       
C   S-VARIABLES 
      SVARS(1) = STRESS_SD_MECH(1,1)
      SVARS(2) = STRESS_SD_MECH(2,1)
      SVARS(3) = STRESS_SD_MECH(3,1)
      SVARS(4) = STRAIN_SD_MECH(1,1)
      SVARS(5) = STRAIN_SD_MECH(2,1)
      SVARS(6) = STRAIN_SD_MECH(3,1)
      SVARS(7) = ELE_FIELD_SD(1,1)
      SVARS(8) = ELE_FIELD_SD(2,1)
      SVARS(9) = ELE_DISP(1,1)
      SVARS(10) = ELE_DISP(2,1)
C --------------------------------------------------------------------
      DO I = 1,NDOFEL
        DO J = 1,NDOFEL
           RHS(I,1) = RHS(I,1) - AMATRX(I,J)*U(J)
        END DO 
      END DO
C   WRITING OUT OUTPUT 
      WRITE(7,*) JELEM
      WRITE(7,*) 'AREA1',SVARS(1)
      WRITE(7,*) 'AREA2',SVARS(2)
      WRITE(7,*) 'U1x  =', U(1)
      WRITE(7,*) 'U1y  =', U(2)
      WRITE(7,*) 'Phi1 =', U(3)
      WRITE(7,*) 'U2x  =', U(4)
      WRITE(7,*) 'U2y  =', U(5)
      WRITE(7,*) 'Phi2 =', U(6)
      WRITE(7,*) 'U3x  =', U(7)
      WRITE(7,*) 'U3y  =', U(8)
      WRITE(7,*) 'Phi3 =', U(9)
      WRITE(7,*) 'U4x  =', U(10)
      WRITE(7,*) 'U4y  =', U(11)
      WRITE(7,*) 'Phi4 =', U(12)
C   PRINTING OUT -- OUTPUT MATRIX 
c      WRITE(7,*) 'OUTPUT MATRIX :::'
c      WRITE(7,*) OUTPUT_MATRIX      
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
C   PIEZOELECTRIC MATRIX 
        SUBROUTINE PIEZOELECTRIC_MATRIX(eMATRIX,e31,e33,e15)
        INCLUDE 'ABA_PARAM.INC'
        REAL(8) e15,e31,e33
        DOUBLE PRECISION, DIMENSION(2,3) :: eMATRIX
        eMATRIX(1,1) = 0.0
        eMATRIX(1,2) = 0.0
        eMATRIX(1,3) = e15
        eMATRIX(2,1) = e31
        eMATRIX(2,2) = e33
        eMATRIX(2,3) = 0.0
        RETURN
        END
C --------------------------------------------------------------------
C   Transpose Of the Matrix 
        SUBROUTINE TRANSPOSE_MATRIX(MAT,TRANS_MAT,I,J)
        INCLUDE 'ABA_PARAM.INC'
        INTEGER I,J,M,N
        DOUBLE PRECISION, DIMENSION(I,J) :: MAT
        DOUBLE PRECISION, DIMENSION(J,I) :: TRANS_MAT
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
        SUBROUTINE STRAIN_MECH(N1,N2,N3,N4,N5,N6,SHAPE_V,B_u,B_uT,
     1       AREA_CELL)
        INCLUDE 'ABA_PARAM.INC'
        DOUBLE PRECISION, DIMENSION(2) :: N1
        DOUBLE PRECISION, DIMENSION(2) :: N2
        DOUBLE PRECISION, DIMENSION(2) :: N3
        DOUBLE PRECISION, DIMENSION(2) :: N4
        DOUBLE PRECISION, DIMENSION(2) :: N5
        DOUBLE PRECISION, DIMENSION(2) :: N6
        DOUBLE PRECISION, DIMENSION(7,4) :: SHAPE_V
        DOUBLE PRECISION, DIMENSION(2) :: AREA_CELL
        DOUBLE PRECISION, DIMENSION(6,8) :: B_u
        DOUBLE PRECISION, DIMENSION(8,6) :: B_uT
C                -                  -
C               |    STRAIN-DISP    |
C               |        SD1        |
C               |      (3 x 8)      |
C                  ----------------
C               |                   |
C               |       SD2         |
C               |     (3 x 8)       |
C                -                 - (6 x 8)
C   @ NODE 1 -- CELL 1 & 2 
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
        NORMALS(2,1) = -0.5D0 * (X(2)-X(1))/(L(1)/2)
        NORMALS(1,2) = 0.5D0 * (Y(3)-Y(2))/(L(2)/2)
        NORMALS(2,2) = -0.5D0 * (X(3)-X(2))/(L(2)/2)
        NORMALS(1,3) = 0.5D0 * (Y(4)-Y(3))/(L(3)/2)
        NORMALS(2,3) = -0.5D0 * (X(4)-X(3))/(L(3)/2)
        NORMALS(1,4) = 0.5D0 * (Y(1)-Y(4))/(L(4)/2)
        NORMALS(2,4) = -0.5D0 * (X(1)-X(4))/(L(4)/2)
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
c        PRINT *, 'bI****'
c        PRINT *, Bi(1,1)
c        PRINT *, Bi(2,2) 
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
c        PRINT *, LENGTH
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
c        PRINT *, AREA_CELL
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
        RETURN
        END
C --------------------------------------------------------------------
C   MATRIX-MULTIPLICATION 
        SUBROUTINE MATRIX_MULTIPLICATION(A,B,C,L,N,M)
        INCLUDE 'ABA_PARAM.INC'
        INTEGER I,J,L,N,M,K
        DOUBLE PRECISION, DIMENSION(L,N) :: A
        DOUBLE PRECISION, DIMENSION(N,M) :: B
        DOUBLE PRECISION, DIMENSION(L,M) :: C
        CALL KASET2(C,L,M)
        DO I = 1,L
            DO J = 1,M 
                DO K = 1,N
                    C(I,J) = C(I,J) + A(I,K)*B(K,J)
                END DO 
            END DO 
        END DO
        RETURN
        END
C --------------------------------------------------------------------
        SUBROUTINE STRAIN_ELE(N1,N2,N3,N4,N5,N6,SHAPE_V,Bphi,BphiT,
     1   AREA_CELL)
        INCLUDE 'ABA_PARAM.INC'
        DOUBLE PRECISION, DIMENSION(2) :: N1
        DOUBLE PRECISION, DIMENSION(2) :: N2
        DOUBLE PRECISION, DIMENSION(2) :: N3
        DOUBLE PRECISION, DIMENSION(2) :: N4
        DOUBLE PRECISION, DIMENSION(2) :: N5
        DOUBLE PRECISION, DIMENSION(2) :: N6
        DOUBLE PRECISION, DIMENSION(2) :: AREA_CELL
        DOUBLE PRECISION, DIMENSION(7,4) :: SHAPE_V
        DOUBLE PRECISION, DIMENSION(4,4) :: Bphi
        DOUBLE PRECISION, DIMENSION(4,4) :: BphiT
C               -                  -
C               |    STRAIN-DISP    |
C               |        SD1        |
C               |      (2 x 4)      |
C                 ----------------
C               |                   |
C               |       SD2         |
C               |     (2 x 4)       |
C               -                  - (4 x 4)        
C   @ NODE 1     
        CALL B_phiI(N1,N5,N6,N4,SHAPE_V(1,1),SHAPE_V(2,1),
     1  SHAPE_V(3,1),SHAPE_V(4,1),Bphi(1:2,1),BphiT(1,1:2),
     2  AREA_CELL(1))
        CALL B_phiI(N5,N2,N3,N6,SHAPE_V(5,1),SHAPE_V(6,1),
     1  SHAPE_V(7,1),SHAPE_V(2,1),Bphi(3:4,1),BphiT(1,3:4),
     2  AREA_CELL(2))
C   @ NODE 2
        CALL B_phiI(N1,N5,N6,N4,SHAPE_V(1,2),SHAPE_V(2,2),
     1  SHAPE_V(3,2),SHAPE_V(4,2),Bphi(1:2,2),BphiT(2,1:2),
     2  AREA_CELL(1))
        CALL B_phiI(N5,N2,N3,N6,SHAPE_V(5,2),SHAPE_V(6,2),
     1  SHAPE_V(7,2),SHAPE_V(2,2),Bphi(3:4,2),BphiT(2,3:4),
     2  AREA_CELL(2))
C   @ NODE 3
        CALL B_phiI(N1,N5,N6,N4,SHAPE_V(1,3),SHAPE_V(2,3),
     1  SHAPE_V(3,3),SHAPE_V(4,3),Bphi(1:2,3),BphiT(3,1:2),
     2  AREA_CELL(1))
        CALL B_phiI(N5,N2,N3,N6,SHAPE_V(5,3),SHAPE_V(6,3),
     1  SHAPE_V(7,3),SHAPE_V(2,3),Bphi(3:4,3),BphiT(3,3:4),
     2  AREA_CELL(2))
C   @ NODE 4 
        CALL B_phiI(N1,N5,N6,N4,SHAPE_V(1,4),SHAPE_V(2,4),
     1  SHAPE_V(3,4),SHAPE_V(4,4),Bphi(1:2,4),BphiT(4,1:2),
     2  AREA_CELL(1))
        CALL B_phiI(N5,N2,N3,N6,SHAPE_V(5,4),SHAPE_V(6,4),
     1  SHAPE_V(7,4),SHAPE_V(2,4),Bphi(3:4,4),BphiT(4,3:4),
     2  AREA_CELL(2))
        RETURN
        END
C --------------------------------------------------------------------
C   STRAIN-ELECTRIC POTENTIAL COMPONENT FOR EACH NODE
        SUBROUTINE B_phiI(P,Q,R,S,INT_PT1,INT_PT2,INT_PT3,INT_PT4,
     1   Bi,BiT,AREA_CELL)
        INCLUDE 'ABA_PARAM.INC'
        REAL(8) INT_PT1,INT_PT2,INT_PT3,INT_PT4,AREA_CELL
        REAL(8) X(4),Y(4),L(4)
        REAL(8) NORMALS(2,4),Bi(2,1),BiT(1,2)
        REAL(8) P(2),Q(2),R(2),S(2)
        REAL(8) BIx,BIy
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
        NORMALS(2,1) = -0.5D0 * (X(2)-X(1))/(L(1)/2)
        NORMALS(1,2) = 0.5D0 * (Y(3)-Y(2))/(L(2)/2)
        NORMALS(2,2) = -0.5D0 * (X(3)-X(2))/(L(2)/2)
        NORMALS(1,3) = 0.5D0 * (Y(4)-Y(3))/(L(3)/2)
        NORMALS(2,3) = -0.5D0 * (X(4)-X(3))/(L(3)/2)
        NORMALS(1,4) = 0.5D0 * (Y(1)-Y(4))/(L(4)/2)
        NORMALS(2,4) = -0.5D0 * (X(1)-X(4))/(L(4)/2)  
C   COMPONENT OF STRAIN-DISPLACEMENT FOR NODE      
        BIx = (NORMALS(1,1)*INT_PT1*L(1) + NORMALS(1,2)*INT_PT2*L(2) +
     1       NORMALS(1,3)*INT_PT3*L(3) + NORMALS(1,4)*INT_PT4*L(4))
        BIx = BIx / AREA_CELL
        BIy = (NORMALS(2,1)*INT_PT1*L(1) + NORMALS(2,2)*INT_PT2*L(2) +
     1       NORMALS(2,3)*INT_PT3*L(3) + NORMALS(2,4)*INT_PT4*L(4))
        BIy = BIy / AREA_CELL
        CALL KASET2(Bi,2,1)
        Bi(1,1) = Bi(1,1) + BIx
        Bi(2,1) = Bi(2,1) + BIy
c        PRINT *, '*******bi_Electrical***********'
c        PRINT *, Bi
        CALL TRANSPOSE_MATRIX(Bi,BiT,2,1)
        RETURN
        END
C --------------------------------------------------------------------
C   Kuu MATRIX 
        SUBROUTINE KMECH(KA,CMATRIX,KB,AREA,TEMP)
        INCLUDE 'ABA_PARAM.INC'
        DOUBLE PRECISION, DIMENSION(2,3) :: KA 
        DOUBLE PRECISION, DIMENSION(3,3) :: CMATRIX
        DOUBLE PRECISION, DIMENSION(3,2) :: KB 
        DOUBLE PRECISION, DIMENSION(2,2) :: TEMP
        DOUBLE PRECISION, DIMENSION(2,3) :: KDUMMY
        REAL(8) AREA
        CALL KASET2(KDUMMY,2,3)
        CALL MATRIX_MULTIPLICATION(KA,CMATRIX,KDUMMY,2,3,3)
        CALL MATRIX_MULTIPLICATION(KDUMMY,KB,TEMP,2,3,2)
        DO I = 1,2
            DO J = 1,2 
                TEMP(I,J) = TEMP(I,J)*AREA
            END DO 
        END DO
        RETURN
        END
C --------------------------------------------------------------------
C   K-MECH-ELE MATRIX 
        SUBROUTINE KMECH_ELE(KA,eTMATRIX,KB,AREA,TEMP2)
        INCLUDE 'ABA_PARAM.INC'
        DOUBLE PRECISION, DIMENSION(2,3) :: KA
        DOUBLE PRECISION, DIMENSION(3,2) :: eTMATRIX
        DOUBLE PRECISION, DIMENSION(2,1) :: KB
        DOUBLE PRECISION, DIMENSION(2,2) :: KDUMMY2
        DOUBLE PRECISION, DIMENSION(2,1) :: TEMP2
        REAL(8) AREA
        INTEGER I
        CALL KASET2(TEMP2,2,1)
        CALL KASET2(KDUMMY2,2,2)
        CALL MATRIX_MULTIPLICATION(KA,eTMATRIX,KDUMMY2,2,3,2)
        CALL MATRIX_MULTIPLICATION(KDUMMY2,KB,TEMP2,2,2,1)
        DO I = 1,2 
            TEMP2(I,1) = TEMP2(I,1)*AREA
        END DO 
        RETURN
        END
C --------------------------------------------------------------------
C KELE MATRIX 
        SUBROUTINE KELE(KA,gTMATRIX,KB,AREA,TEMP3)
        INCLUDE 'ABA_PARAM.INC'
        DOUBLE PRECISION, DIMENSION(1,2) :: KA
        DOUBLE PRECISION, DIMENSION(2,2) :: gTMATRIX
        DOUBLE PRECISION, DIMENSION(2,1) :: KB
        DOUBLE PRECISION, DIMENSION(1,2) :: KDUMMY3
        DOUBLE PRECISION, DIMENSION(1,1) :: TEMP3
        REAL(8) AREA
        CALL KASET2(TEMP3,1,1)
        CALL MATRIX_MULTIPLICATION(KA,gTMATRIX,KDUMMY3,1,2,2)
        CALL MATRIX_MULTIPLICATION(KDUMMY3,KB,TEMP3,1,2,1)
        TEMP3 = -TEMP3 * AREA
        RETURN
        END
C --------------------------------------------------------------------        
