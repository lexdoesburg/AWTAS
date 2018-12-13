      module MatrixSolvers

	contains

!----------------------------------------------------------------------------------

cc        
      SUBROUTINE SEVEN1(C1,C2,C3,C4,C5,C6,C7,F,XXX,N)
      IMPLICIT REAL*8 (A-H,O-Z)
!     AC: Altered dimensions from 400 to N (=2M):
      DIMENSION C1(N),C2(N),C3(N),C4(N)
      DIMENSION C5(N),C6(N),C7(N)
      DIMENSION XXX(N),F(N)
      M=N/2
      J=1
      J1=2
      IF(C7(J).NE.0.0) THEN
        FX=C7(J)/C6(J1)
        C4(J)= C4(J)-C3(J1)*FX
        C5(J)= C5(J)-C4(J1)*FX
        C6(J)= C6(J)-C5(J1)*FX
        C7(J)= C7(J)-C6(J1)*FX
        F(J)= F(J)-F(J1)*FX
	ENDIF
      DO 30 I=2,M-1
         J=2*I-1
         J1=J+1
      IF(C7(J).NE.0.0) THEN
        FX=C7(J)/C6(J1)
        C2(J)= C2(J)-C1(J1)*FX
        C3(J)= C3(J)-C2(J1)*FX
        C4(J)= C4(J)-C3(J1)*FX
        C5(J)= C5(J)-C4(J1)*FX
        C6(J)= C6(J)-C5(J1)*FX
        C7(J)= C7(J)-C6(J1)*FX
        F(J)= F(J)-F(J1)*FX
      ENDIF
      IF(C1(J1).NE.0.0) THEN
        FX=C1(J1)/C2(J)
        C1(J1)= C1(J1)-C2(J)*FX
        C2(J1)= C2(J1)-C3(J)*FX
        C3(J1)= C3(J1)-C4(J)*FX
        C4(J1)= C4(J1)-C5(J)*FX
        C5(J1)= C5(J1)-C6(J)*FX
        C6(J1)= C6(J1)-C7(J)*FX
        F(J1)= F(J1)-F(J)*FX
      ENDIF
30    CONTINUE
      J=2*M-1
      J1=J+1
      IF(C1(J1).NE.0.0) THEN
        FX=C1(J1)/C2(J)
        C1(J1)= C1(J1)-C2(J)*FX
        C2(J1)= C2(J1)-C3(J)*FX
        C3(J1)= C3(J1)-C4(J)*FX
        C4(J1)= C4(J1)-C5(J)*FX
        F(J1)= F(J1)-F(J)*FX
      ENDIF
cc
      CALL FIVED(C2,C3,C4,C5,C6,F,XXX,N)
cc
      RETURN
      END SUBROUTINE SEVEN1
cc       
cc----7--1--------2---------3---------4---------5---------6---------7---------8 
cc        
!      SUBROUTINE SEVEN2(C1,C2,C3,C4,C5,C6,C7,F,XXX,N)
cc
cc  This is an extra solver which uses the "Numerical Recipes" LU stuff
cc  It can be all removed
!   Commented out- AC (7/7/99)
cc
!      IMPLICIT REAL*8 (A-H,O-Z)
!      DIMENSION C1(400),C2(400),C3(400),C4(400)
!      DIMENSION C5(400),C6(400),C7(400)
!      DIMENSION XXX(400),F(400)
!      DIMENSION Z(400,400)
!      DIMENSION INDX(400)
!      COMMON/SEV2/Z
!      COMMON/SEV1/INDX
!      NP=400
!      M=N/2
!      DO 31 I=1,N
!         DO 32 J=1,N
!            Z(I,J)=0.0
!32       CONTINUE
!31    CONTINUE
!      J=1
!      J1=2
!      K4=J
!      K5=J+1
!      K6=J+2
!      K7=J+3
!      Z(J,K4)=C4(J)
!      Z(J,K5)=C5(J)
!      Z(J,K6)=C6(J)
!      Z(J,K7)=C7(J)
!      Z(J1,K4)=C3(J1)
!      Z(J1,K5)=C4(J1)
!      Z(J1,K6)=C5(J1)
!      Z(J1,K7)=C6(J1)
!      DO 33 I=2,M-1
!         J=2*I-1
!         J1=J+1
!	   K0=J-3
!         K1=J-2
!         K2=J-1
!         K3=J
!         K4=J+1
!         K5=J+2
!         K6=J+3
!	   K7=J+4
!	IF(I.GT.2) THEN
!	  Z(J,K0)=C1(J)
!	ENDIF
!        Z(J,K1)=C2(J)
!        Z(J,K2)=C3(J)
!        Z(J,K3)=C4(J)
!        Z(J,K4)=C5(J)
!        Z(J,K5)=C6(J)
!        Z(J,K6)=C7(J)
!        Z(J1,K1)=C1(J1)
!        Z(J1,K2)=C2(J1)
!        Z(J1,K3)=C3(J1)
!        Z(J1,K4)=C4(J1)
!        Z(J1,K5)=C5(J1)
!        Z(J1,K6)=C6(J1)
!	IF(I.LT.M-1) THEN
!	  Z(J1,K7)=C7(J1)
!	ENDIF
!33    CONTINUE
!      J=2*M-1
!      J1=J+1
!	K0=J-3
!	K1=J-2
!      K2=J-1
!      K3=J
!      K4=J+1
!      K5=J+2
!	Z(J,K0)=C1(J)
!      Z(J,K1)=C2(J)
!      Z(J,K2)=C3(J)
!      Z(J,K3)=C4(J)
!      Z(J,K4)=C5(J)
!      Z(J1,K1)=C1(J1)
!      Z(J1,K2)=C2(J1)
!      Z(J1,K3)=C3(J1)
!      Z(J1,K4)=C4(J1)
cc
!      CALL SOLVEX(Z,N,NP,INDX,F,XXX)   
!      RETURN
!      END
cc       
cc----7--1--------2---------3---------4---------5---------6---------7---------8 
cc        
      SUBROUTINE FIVED(A,B,C,D,E,F,V,L)
      IMPLICIT REAL*8 (A-H,O-Z)
!     AC: Altered dimensions from 400 to L:
      DIMENSION A(L),B(L),C(L),D(L),E(L),F(L),V(L)
      DIMENSION OME(L),BET(L),GAM(L),DEL(L),H(L)
cc
cc---This subroutine solves the 5-diagonal system
cc---
cc---   C(1)*V(1)+D(1)*V(2)+E(1)*V(3)=F(1)
cc---
cc---   B(2)*V(1)+C(2)*V(2)+D(2)*V(3)+E(2)*V(4)=F(2)
cc---
cc---   A(3)*V(1)+B(3)*V(2)+C(3)*V(3)+D(3)*V(4)+E(3)*V(5)=F(3)
cc---              ...
cc---   A(I)*V(I-2)+B(I)*V(I-1)+C(I)*V(I)+D(I)*V(I+1)+E(I)*V(I+2)=F(I)
cc---              ...
cc---   A(L-1)*V(L-3)+B(L-1)*V(L-2)+C(L-1)*V(L-1)+D(L-1)*V(L)=F(L-1)
cc---
cc---        A(L)*V(L-2)+B(L)*V(L-1)+C(L)*V(L)=F(L)
cc---
      LL=L-1
      OME1=C(1)
      OME(1)=OME1
      BET1=D(1)/OME1
      BET(1)=BET1
      GAM1=E(1)/OME1
      GAM(1)=GAM1
      DEL2=B(2)
      DEL(2)=DEL2
      OME2=C(2)-DEL2*BET1
      OME(2)=OME2
      BET(2)=(D(2)-DEL2*GAM1)/OME2
      GAM(2)=E(2)/OME2
      do N=3,L-2
         DELN=B(N)-A(N)*BET(N-2)
         DEL(N)=DELN
         OMEN=C(N)-A(N)*GAM(N-2)-DELN*BET(N-1)
         OME(N)=OMEN
         BET(N)=(D(N)-DELN*GAM(N-1))/OMEN
         GAM(N)=E(N)/OMEN
      end do
      N=L-1
      AN=A(N)
      DELN=B(N)-AN*BET(N-2)
      DEL(N)=DELN
      OMEN=C(N)-AN*GAM(N-2)-DELN*BET(N-1)
      OME(N)=OMEN
      BET(N)=(D(N)-DELN*GAM(N-1))/OMEN
      GAM(N)=0.0
      N=L
      AN=A(N)
      DELN=B(N)-AN*BET(N-2)
      DEL(N)=DELN
      OMEN=C(N)-AN*GAM(N-2)-DELN*BET(N-1)
      OME(N)=OMEN
      BET(N)=0.0
      GAM(N)=0.0
      H(1)=F(1)/OME(1)
      H(2)=(F(2)-DEL(2)*H(1))/OME(2)
      do N=3,L
         H(N)=(F(N)-A(N)*H(N-2)-DEL(N)*H(N-1))/OME(N)
      end do
      V(L)=H(L)
      V(LL)=H(LL)-BET(LL)*V(L)
      do N=2,LL
         NN=LL+1-N
         V(NN)=H(NN)-BET(NN)*V(NN+1)-GAM(NN)*V(NN+2)
      end do
      RETURN
      END SUBROUTINE FIVED
cc       

cc----7--1--------2---------3---------4---------5---------6---------7---------8 
cc
      SUBROUTINE LUDCMP(A,N,NP,INDX,D)   
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NMAX=100, TINY=1.0D-20)
      DIMENSION A(NP,NP),VV(NMAX)
      DIMENSION INDX(N)
      D=1.0
      DO 12 I=1,N
         AAMAX=0.0
         DO 11 J=1,N
            IF(ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
11          CONTINUE
         IF(AAMAX.EQ.0.0) PAUSE 'singular matrix'
         VV(I)=1.0/AAMAX
12       CONTINUE
      DO 19 J=1,N
         DO 14 I=1,J-1
            SUM=A(I,J)
            DO 13 K=1,I-1
               SUM=SUM-A(I,K)*A(K,J)
13             CONTINUE
            A(I,J)=SUM
14          CONTINUE
         AAMAX=0.0
         DO 16 I=J,N
            SUM=A(I,J)
            DO 15 K=1,J-1
               SUM=SUM-A(I,K)*A(K,J)
15             CONTINUE
            A(I,J)=SUM
            DUM=VV(I)*ABS(SUM)
            IF(DUM.GE.AAMAX) THEN
               IMAX=I
               AAMAX=DUM
            ENDIF
16          CONTINUE
         IF(J.NE.IMAX) THEN
            DO 17 K=1,N
               DUM=A(IMAX,K)
               A(IMAX,K)=A(J,K)
               A(J,K)=DUM
17             CONTINUE
            D=-D
            VV(IMAX)=VV(J)
         ENDIF
         INDX(J)=IMAX
         IF(A(J,J).EQ.0.0) A(J,J)=TINY
         IF(J.NE.N) THEN
            DUM=1.0/A(J,J)
            DO 18 I=J+1,N
               A(I,J)=A(I,J)*DUM
18             CONTINUE
         ENDIF
19       CONTINUE
      RETURN
      END SUBROUTINE LUDCMP
cc
cc ==================================================================
cc
      SUBROUTINE LUBKSB(A,N,NP,INDX,B)   
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(NP,NP),B(N)
      DIMENSION INDX(N)
      II=0
      DO 12 I=1,N
         LL=INDX(I)
         SUM=B(LL)
         B(LL)=B(I)
         IF(II.NE.0) THEN
             DO 11 J=II,I-1
                SUM=SUM-A(I,J)*B(J)
11              CONTINUE
         ELSE IF(SUM.NE.0.0) THEN
             II=I
         ENDIF
         B(I)=SUM
12       CONTINUE
      DO 14 I=N,1,-1
         SUM=B(I)
         IF(I.LT.N) THEN
             DO 13 J=I+1,N
                SUM=SUM-A(I,J)*B(J)
13              CONTINUE
         ENDIF
         B(I)=SUM/A(I,I)
14       CONTINUE
      RETURN
      END SUBROUTINE LUBKSB
cc
cc ==================================================================
cc
      SUBROUTINE SOLVEX(A,N,NP,INDX,B,X)   
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(NP,NP),B(N),X(N)
	DIMENSION INDX(N)
      CALL LUDCMP(A,N,NP,INDX,D)
      CALL LUBKSB(A,N,NP,INDX,B)
      DO 1 J=1,N
1     X(J)=B(J)
      RETURN
      END SUBROUTINE SOLVEX
  
!----------------------------------------------------------------------------------

	end module MatrixSolvers