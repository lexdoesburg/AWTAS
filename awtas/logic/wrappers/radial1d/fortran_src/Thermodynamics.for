      module Thermodynamics

      use variable_types

!     IGOOD is used as a flag for reporting errors in the thermodynamic routines.
      integer(I4B):: IGOOD

      contains

cc----7--1--------2---------3---------4---------5---------6---------7---------8 
cc        
      SUBROUTINE ACCUM(BM,BE,IPH,P,T,SV,RHOL,RHOV,HL,HV,POR,RHOR,CR,
     2                  P0,T0,COMP,COMT,AAA,PERFAC)
      implicit none
	real(DP):: BM,BE,P,T,SV,RHOL,RHOV,HL,HV,POR,RHOR,CR,P0,T0
	real(DP):: COMP,COMT,AAA,PERFAC,PORN,SL
	integer(I4B):: IPH
cc
cc  This subroutine calculates the mass and energy accumulation terms BM and BE
cc  It requires values for IPH,P,T,SV,RHOL,RHOV,HL,HV,POR,RHOR,CR
cc      
cc  Calculate modified porosity and permeability
cc
      PORN=POR*(1.0D+00+COMP*(P-P0)+COMT*(T-T0))
      CALL PERF(POR,PORN,AAA,PERFAC)
cc 
      select case (IPH)
cc  single phase liquid
      case(0)
          BM=PORN*RHOL
          BE=PORN*(RHOL*HL-P)+(1.0D0 - PORN)*RHOR*CR*T
cc  two-phase
      case(1)
	    SL=1.0D0-SV
          BM=PORN*(SL*RHOL+SV*RHOV)
          BE=PORN*(SL*RHOL*HL+SV*RHOV*HV-P)+(1.0D0 - PORN)*RHOR*CR*T
cc  single phase steam
      case(2)
          BM=PORN*RHOV
          BE=PORN*(RHOV*HV-P)+(1.0D0 - PORN)*RHOR*CR*T
      end select
      RETURN
      end subroutine ACCUM
cc       
cc----7--1--------2---------3---------4---------5---------6---------7---------8 
cc        
      SUBROUTINE TRANS(TM,TE,IPH,SV,HL,HV,VISL,VISV,PER,PERFAC)

      implicit none       
	real(DP):: TM,TE,SV,HL,VISL,VISV,PER,PERFAC
	integer(I4B):: IPH
	real(DP):: HV,XKRL,XKRV
	integer(I4B):: IRELP
cc
cc  This subroutine calculates the mass and energy transmissibility terms TM and TE
cc  It requires values of SV,HL,HV,VISL,VISV. It calculates relative permeabilities if necessary
cc      
      select case (IPH)
cc  single phase liquid
      case(0)
          TM=PER*PERFAC*1.0D0/VISL
          TE=HL*TM
cc  two-phase
      case(1)
          CALL RELP(SV,IRELP,XKRL,XKRV)
          TM=PER*PERFAC*(XKRL/VISL + XKRV/VISV)
          TE=PER*PERFAC*(HL*XKRL/VISL + HV*XKRV/VISV)
cc  single phase steam
      case(2)
           TM=PER*PERFAC*1.0D0/VISV
           TE=HV*TM
      end select
      RETURN
      end subroutine TRANS
cc       
cc----7--1--------2---------3---------4---------5---------6---------7---------8 
cc        
      SUBROUTINE TERMS(PX,TX,SVX,IPHX,PX1,TX1,TX2,
     1             BMX,BEX,TMX,TEX,BMX1,BEX1,TMX1,TEX1,
     2             BMX2,BEX2,TMX2,TEX2,DELPX,DELXX,
     3             PORX,CR,RHOR,FACP,FACT,FACS,
     4             PERX,P0X,T0X,COMP,COMT,AAA)

      implicit none
	real(DP):: PX,TX,SVX,PX1,TX1,TX2,BMX,BEX,TMX,TEX,BMX1,BEX1
	real(DP):: TMX1,TEX1,BMX2,BEX2,TMX2,TEX2,DELPX,DELXX
      real(DP):: PORX,CR,RHOR,FACP,FACT,FACS,PERX,P0X,T0X,COMP,COMT,AAA
	integer(I4B):: IPHX
      real(DP):: RHOL,RHOV,HL,HV,VISL,VISV,PERFAC,RHOL1,RHOV1,HL1,HV1
	real(DP):: VISL1,VISV1,PERFAC1,RHOL2,RHOV2,HL2,HV2,VISL2,VISV2
	real(DP):: PERFAC2,T0,SVX2
cc
cc  This subroutine calculates accumulation and transmissibility terms at each block
cc  It also perturbs the primary variables and recalculates them
cc
cc  First calculate the basic terms
cc
      CALL THERMO(PX,TX,SVX,IPHX,RHOL,RHOV,HL,HV,VISL,VISV)
	if (IGOOD>0) then
	!   print *, 'IGOOD 1'
	  return
	end if
      CALL ACCUM(BMX,BEX,IPHX,PX,TX,SVX,
     1	       RHOL,RHOV,HL,HV,PORX,RHOR,CR,
     2         P0X,T0X,COMP,COMT,AAA,PERFAC)
      CALL TRANS(TMX,TEX,IPHX,SVX,
     1           HL,HV,VISL,VISV,PERX,PERFAC)
cc
cc  now perturb the primary variables to calculate derivatives
cc
      select case (IPHX)
cc  -------------------------------------------------------------------------------
cc  single phase liquid
cc
      case (0)
cc
cc  make sure it stays liquid - increase pressure
cc
        DELPX=FACP*PX
        PX1=PX+DELPX
	  TX1=TX
        CALL THERMO(PX1,TX,SVX,IPHX,RHOL1,RHOV1,HL1,HV1,VISL1,VISV1)
	  if (IGOOD>0) then
	!     print *, 'IGOOD 2'
	    return
	  end if
        CALL ACCUM(BMX1,BEX1,IPHX,PX1,TX,SVX,
     1             RHOL1,RHOV1,HL1,HV1,PORX,RHOR,CR,
     2             P0X,T0X,COMP,COMT,AAA,PERFAC1)
        CALL TRANS(TMX1,TEX1,IPHX,SVX,
     1             HL1,HV1,VISL1,VISV1,PERX,PERFAC1)
cc
cc  make sure it stays liquid - decrease temperature
        DELXX=-FACT*TX
        TX2=TX+DELXX
        CALL THERMO(PX,TX2,SVX,IPHX,RHOL2,RHOV2,HL2,HV2,VISL2,VISV2)
        if (IGOOD>0) then
      !     print *, 'IGOOD 3'
          return
        end if
        CALL ACCUM(BMX2,BEX2,IPHX,PX,TX2,SVX,
     1             RHOL2,RHOV2,HL2,HV2,PORX,RHOR,CR,
     2             P0X,T0X,COMP,COMT,AAA,PERFAC2)
        CALL TRANS(TMX2,TEX2,IPHX,SVX,
     1             HL2,HV2,VISL2,VISV2,PERX,PERFAC2)
cc
cc  -------------------------------------------------------------------------------
cc  2-phase conditions
cc
      case (1)
cc  
cc  increase pressure
        DELPX=FACP*PX
        PX1=PX+DELPX
        TX2=TX
	  T0=0.0D0
        CALL TSAT(PX1,T0,TX1)
	! Added in case of problems in TSAT (AC 8/00):
	  if (IGOOD>0) then
	!     print *, 'IGOOD 4'
	    return
	  end if
        CALL THERMO(PX1,TX1,SVX,IPHX,RHOL1,RHOV1,HL1,HV1,VISL1,VISV1)
	  if (IGOOD>0) then
	!     print *, 'IGOOD 5'
	    return
	  end if
        CALL ACCUM(BMX1,BEX1,IPHX,PX1,TX1,SVX,
     1             RHOL1,RHOV1,HL1,HV1,PORX,RHOR,CR,
     2             P0X,T0X,COMP,COMT,AAA,PERFAC1)
        CALL TRANS(TMX1,TEX1,IPHX,SVX,
     1             HL1,HV1,VISL1,VISV1,PERX,PERFAC1)
cc  increase saturation
        DELXX=FACS*SVX+1.0D-10
        IF(SVX+DELXX>=1.0D0) DELXX=-DELXX
        SVX2=SVX+DELXX
        CALL ACCUM(BMX2,BEX2,IPHX,PX,TX,SVX2,
     1             RHOL,RHOV,HL,HV,PORX,RHOR,CR,
     2             P0X,T0X,COMP,COMT,AAA,PERFAC2)
        CALL TRANS(TMX2,TEX2,IPHX,SVX2,
     1             HL,HV,VISL,VISV,PERX,PERFAC2)
cc
cc  -------------------------------------------------------------------------------
cc  single phase gas
cc
      case (2)
cc
cc  make sure it stays gas - decrease pressure
cc
        DELPX=-FACP*PX
        PX1=PX+DELPX
	  TX1=TX
        CALL THERMO(PX1,TX,SVX,IPHX,RHOL1,RHOV1,HL1,HV1,VISL1,VISV1)
	  if (IGOOD>0) return
        CALL ACCUM(BMX1,BEX1,IPHX,PX1,TX,SVX,
     1             RHOL1,RHOV1,HL1,HV1,PORX,RHOR,CR,
     2             P0X,T0X,COMP,COMT,AAA,PERFAC1)
        CALL TRANS(TMX1,TEX1,IPHX,SVX,
     1             HL1,HV1,VISL1,VISV1,PERX,PERFAC1)
cc
cc  make sure it stays gas - increase temperature
        DELXX=FACT*TX
        TX2=TX+DELXX
        CALL THERMO(PX,TX2,SVX,IPHX,RHOL2,RHOV2,HL2,HV2,VISL2,VISV2)
	  if (IGOOD>0) return
        CALL ACCUM(BMX2,BEX2,IPHX,PX,TX2,SVX,
     1             RHOL2,RHOV2,HL2,HV2,PORX,RHOR,CR,
     2             P0X,T0X,COMP,COMT,AAA,PERFAC2)
        CALL TRANS(TMX2,TEX2,IPHX,SVX,
     1             HL2,HV2,VISL2,VISV2,PERX,PERFAC2)
      end select
      return
      end subroutine TERMS
cc       
!----------------------------------------------------------------------------------	 
cc        
      SUBROUTINE COWAT(TF,PP,D,U)
C--------- Fast COWAT M.J.O'Sullivan - 17 SEPT 1990 ---------
C     20 September 1990.  VAX needs double precision, CRAY does not.
cc 
cc  This is the fast version of COWAT. Karsten has used it in version 2.0 of TOUGH2

!     AC: this routine is too horrendous for proper variable typing.
cc
      IMPLICIT REAL(DP) (A-H,O-Z)
      DATA A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12 /
     16.824687741E3,-5.422063673E2,-2.096666205E4,3.941286787E4,
     2-13.466555478E4,29.707143084E4,-4.375647096E5,42.954208335E4,
     3-27.067012452E4,9.926972482E4,-16.138168904E3,7.982692717E0/
      DATA A13,A14,A15,A16,A17,A18,A19,A20,A21,A22,A23 /
     4-2.616571843E-2,1.522411790E-3,2.284279054E-2,2.421647003E2,
     51.269716088E-10,2.074838328E-7,2.174020350E-8,1.105710498E-9,
     61.293441934E1,1.308119072E-5,6.047626338E-14/
      DATA SA1,SA2,SA3,SA4,SA5,SA6,SA7,SA8,SA9,SA10,SA11,SA12 /
     18.438375405E-1,5.362162162E-4,1.720000000E0,7.342278489E-2,
     24.975858870E-2,6.537154300E-1,1.150E-6,1.51080E-5,
     31.41880E-1,7.002753165E0,2.995284926E-4,2.040E-1   /
C
      TKR=(TF+273.15)/647.3
      TKR2=TKR*TKR
      TKR3=TKR*TKR2
      TKR4=TKR2*TKR2
      TKR5=TKR2*TKR3
      TKR6=TKR4*TKR2
      TKR7=TKR4*TKR3
      TKR8=TKR4*TKR4
      TKR9=TKR4*TKR5
      TKR10=TKR4*TKR6
      TKR11=TKR*TKR10
      TKR19=TKR8*TKR11
      TKR18=TKR8*TKR10
      TKR20=TKR10*TKR10
      PNMR=PP/2.212E7
      PNMR2=PNMR*PNMR
      PNMR3=PNMR*PNMR2
      PNMR4=PNMR*PNMR3
      Y=1.-SA1*TKR2-SA2/TKR6
      ZP=SA3*Y*Y-2.*SA4*TKR+2.*SA5*PNMR
      if (ZP<0.D0) then
      ! print *, 'IGOOD 6'
	  IGOOD=2
	  return
	end if
      Z=Y+ DSQRT(ZP)
      if (Z<0.D0) then ! Added AC (8/00)
      ! print *, 'IGOOD 7'
	  IGOOD=2
	  return
	end if
      CZ=Z**(5./17.)
      PAR1=A12*SA5/CZ
      CC1=SA6-TKR
      CC2=CC1*CC1
      CC4=CC2*CC2
      CC8=CC4*CC4
      CC10=CC2*CC8
      AA1=SA7+TKR19
      PAR2=A13+A14*TKR+A15*TKR2+A16*CC10+A17/AA1
      PAR3=(A18+2.*A19*PNMR+3.*A20*PNMR2)/(SA8+TKR11)
      DD1=SA10+PNMR
      DD2=DD1*DD1
      DD4=DD2*DD2
      PAR4=A21*TKR18*(SA9+TKR2)*(-3./DD4+SA11)
      PAR5=3.*A22*(SA12-TKR)*PNMR2+4.*A23/TKR20*PNMR3
      VMKR=PAR1+PAR2-PAR3-PAR4+PAR5
      V=VMKR*3.17E-3
      D=1./V
      YD=-2.*SA1*TKR+6.*SA2/TKR7
      SNUM= A10+A11*TKR
      SNUM=SNUM*TKR + A9
      SNUM=SNUM*TKR + A8
      SNUM=SNUM*TKR + A7
      SNUM=SNUM*TKR + A6
      SNUM=SNUM*TKR + A5
      SNUM=SNUM*TKR + A4
      SNUM=SNUM*TKR2 - A2
      PRT1=A12*(Z*(17.*(Z/29.-Y/12.)+5.*TKR*YD/12.)+SA4*TKR-
     1(SA3-1.)*TKR*Y*YD)/CZ
      PRT2=PNMR*(A13-A15*TKR2+A16*(9.*TKR+SA6)*CC8*CC1
     2+A17*(19.*TKR19+AA1)/(AA1*AA1))
      BB1=SA8+TKR11
      BB2=BB1*BB1
      PRT3=(11.*TKR11+BB1)/BB2*(A18*PNMR+A19*
     3PNMR2+A20*PNMR3)
      EE1=SA10+PNMR
      EE3=EE1*EE1*EE1
      PRT4=A21*TKR18*(17.*SA9+19.*TKR2)*(1./EE3+SA11*PNMR)
      PRT5=A22*SA12*PNMR3+21.*A23/TKR20*PNMR4
      ENTR= A1*TKR - SNUM +PRT1+PRT2-PRT3+PRT4+PRT5
      H=ENTR*70120.4
      U=H-PP*V
      return
      end subroutine COWAT
cc       
cc----7--1--------2---------3---------4---------5---------6---------7---------8 
cc        
      SUBROUTINE SUPST(T,P,D,U)
C--------- Fast SUPST M.J.O'Sullivan - 17 SEPT 1990 ---------
C     20 September 1990.  VAX needs double precision, CRAY does not.
      IMPLICIT real(DP) (A-H,O-Z)
!     AC: this routine is too horrendous for proper variable typing.
      REAL(DP):: I1

      DATA B0,B01,B02,B03,B04,B05/
     116.83599274,28.56067796,-54.38923329,0.4330662834,-0.6547711697,
     28.565182058E-2/
      DATA B11,B12,B21,B22,B23,B31,B32,B41,B42/
     16.670375918E-2,1.388983801,8.390104328E-2,2.614670893E-2,-3.373439
     2453E-2,4.520918904E-1,1.069036614E-1,-5.975336707E-1,-8.847535804E
     3-2/
      DATA B51,B52,B53,B61,B62,B71,B72,B81,B82/
     15.958051609E-1,-5.159303373E-1,2.075021122E-1,1.190610271E-1,-9.86
     27174132E-2,1.683998803E-1,-5.809438001E-2,6.552390126E-3,5.7102186
     349E-4/
      DATA B90,B91,B92,B93,B94,B95,B96/
     11.936587558E2,-1.388522425E3,4.126607219E3,-6.508211677E3,
     25.745984054E3,-2.693088365E3,5.235718623E2/
      DATA SB,SB61,SB71,SB81,SB82/
     17.633333333E-1,4.006073948E-1,8.636081627E-2,-8.532322921E-1,
     23.460208861E-1/
C
      THETA=(T+273.15)/647.3
      BETA=P/2.212E7
      I1=4.260321148
      X=EXP(SB*(1.-THETA))
C
      X2=X*X
      X3=X2*X
      X4=X3*X
      X5=X4*X
      X6=X5*X
      X8=X6*X2
      X10=X6*X4
      X11=X10*X
      X14=X10*X4
      X18=X14*X4
      X19=X18*X
      X24=X18*X6
      X27=X24*X3
C
      THETA2=THETA*THETA
      THETA3=THETA2*THETA
      THETA4=THETA3*THETA
C
      BETA2=BETA*BETA
      BETA3=BETA2*BETA
      BETA4=BETA3*BETA
      BETA5=BETA4*BETA
      BETA6=BETA5*BETA
      BETA7=BETA6*BETA
C
      BETAL=15.74373327-34.17061978*THETA+19.31380707*THETA2
      DBETAL=-34.17061978+38.62761414*THETA
      R=BETA/BETAL
      R2=R*R
      R4=R2*R2
      R6=R4*R2
      R10=R6*R4
C
      CHI2=I1*THETA/BETA
      SC=(B11*X10+B12)*X3
      CHI2=CHI2-SC
      SC=B21*X18+B22*X2+B23*X
      CHI2=CHI2-2*BETA*SC
      SC=(B31*X8+B32)*X10
      CHI2=CHI2-3*BETA2*SC
      SC=(B41*X11+B42)*X14
      CHI2=CHI2-4*BETA3*SC
      SC=(B51*X8+B52*X4+B53)*X24
      CHI2=CHI2-5*BETA4*SC
C
      SD1=1./BETA4+SB61*X14
      SD2=1./BETA5+SB71*X19
      SD3=1./BETA6+(SB81*X27+SB82)*X27
C
      SN=(B61*X+B62)*X11
cover CHI2=CHI2-SN/SD12*4/BETA5
      chi2=chi2-(sn/sd1*4/beta5)/sd1
      SN=(B71*X6+B72)*X18
cover CHI2=CHI2-SN/SD22*5/BETA6
      chi2=chi2-(sn/sd2*5/beta6)/sd2
      SN=(B81*X10+B82)*X14
cover CHI2=CHI2-SN/SD32*6/BETA7
      chi2=chi2-(sn/sd3*6/beta7)/sd3
      SC=B96
      SC=SC*X+B95
      SC=SC*X+B94
      SC=SC*X+B93
      SC=SC*X+B92
      SC=SC*X+B91
      SC=SC*X+B90
      CHI2=CHI2+11.*R10*SC
      V=CHI2*0.00317
      D=1./V
C
      OS1=SB*THETA
      EPS2=0.0+B0*THETA-(-B01+B03*THETA2+2*B04*THETA3+3*B05*THETA4)
      SC=(B11*(1.+13.*OS1)*X10+B12*(1.+3.*OS1))*X3
      EPS2=EPS2-BETA*SC
      SC=B21*(1.+18.*OS1)*X18+B22*(1.+2.*OS1)*X2+B23*(1.+OS1)*X
      EPS2=EPS2-BETA2*SC
      SC=(B31*(1.+18.*OS1)*X8+B32*(1.+10.*OS1))*X10
      EPS2=EPS2-BETA3*SC
      SC=(B41*(1.+25.*OS1)*X11+B42*(1.+14.*OS1))*X14
      EPS2=EPS2-BETA4*SC
      SC=(B51*(1.+32.*OS1)*X8+B52*(1.+28.*OS1)*X4+
     1 B53*(1.+24.*OS1))*X24
      EPS2=EPS2-BETA5*SC
C
      SN6=14.*SB61*X14
      SN7=19.*SB71*X19
      SN8=(54.*SB81*X27+27.*SB82)*X27
      OS5= 1+11.*OS1-OS1*SN6/SD1
      SC=(B61*X*(OS1+OS5)+B62*OS5)*(X11/SD1)
      EPS2=EPS2-SC
      OS6= 1.+24.*OS1-OS1*SN7/SD2
      SC=(B71*X6*OS6+B72*(OS6-6.*OS1))*(X18/SD2)
      EPS2=EPS2-SC
      OS7= 1.+24.*OS1-OS1*SN8/SD3
      SC=(B81*X10*OS7+B82*(OS7-10.* OS1))*(X14/SD3)
      EPS2=EPS2-SC
      OS2=1+THETA*10.0*DBETAL/BETAL
      SC= (OS2+6*OS1)*B96
      SC=SC*X + (OS2+5*OS1)*B95
      SC=SC*X + (OS2+4*OS1)*B94
      SC=SC*X + (OS2+3*OS1)*B93
      SC=SC*X + (OS2+2*OS1)*B92
      SC=SC*X + (OS2+OS1)*B91
      SC=SC*X + OS2*B90
      EPS2=EPS2+BETA*R10*SC
      H=EPS2*70120.4
      U=H-P*V
      RETURN
      END SUBROUTINE SUPST
cc       
cc----7--1--------2---------3---------4---------5---------6---------7---------8 
cc
!     New version of SAT from Mike 14/10/99.

      subroutine SAT(T,P)

      implicit none
      real(DP):: T,P
	real(DP):: A1,A2,A3,A4,A5,A6,A7,A8,A9
	real(DP):: B0,B1,B2,B3,C1,C2
	real(DP):: X,X1,X2,X3,SC,PC,TC,EX

      DATA A1,A2,A3,A4,A5,A6,A7,A8,A9/
     1-7.691234564,-2.608023696E1,-1.681706546E2,6.423285504E1,
     2-1.189646225E2,4.167117320,2.097506760E1,1.E9,6./
	DATA B0,B1,B2,B3,C1,C2/610.768,44.4063,1.4080,3.0667E-02,
     1-1.310279E-4,2.667390E-02/

      if (T<1.d0) then
	  X=T
	  X2=X*X
	  X3=X2*X
	  P=B0+B1*X+B2*X2+B3*X3
	else if ((T>=1.d0).and.(T<=350.d0)) then
        TC=(T+273.15d0)/647.3d0
        X1=1.-TC
        X2=X1*X1
        SC=A5*X1+A4
        SC=SC*X1+A3
        SC=SC*X1+A2
        SC=SC*X1+A1
        SC=SC*X1
        PC=EXP(SC/(TC*(1.d0+A6*X1+A7*X2))-X1/(A8*X2+A9))
        P=PC*2.212d7
	else if (T>350.d0) then
	  X=(T-350.d0)*0.2
	  EX=7.218407d0+C1*X+C2*X*X
	  P=10.d0**EX
	end if

      return
      end subroutine SAT
       
cc----7--1--------2---------3---------4---------5---------6---------7---------8 
cc
!     Old TSAT

!      SUBROUTINE TSAT(PX,TX00,TS)
!      IMPLICIT REAL*8 (A-H,O-Z)       
C
C-----FIND SATURATION TEMPERATURE TS AT PRESSURE PX.
C
C     TX0 IS STARTING TEMPERATURE FOR ITERATION.
C
!      TX0=TX00
!      IF(TX0.NE.0.) GOTO 2
C
C-----COME HERE TO OBTAIN ROUGH STARTING VALUE FOR ITERATION.
!      IF(PX.GE.656.613) THEN
!         TX0=4606./(24.02-LOG(PX)) - 273.15
!         TX0=MAX(TX0,5.0D0)
!	ELSEIF(PX.LT.656.613) THEN
!	   CC=PX-610.768
!	   ACC=ABS(CC)
!	   SIGN=CC/ACC
!	   TX0=(ACC/3.0667E-02)**(1.0/3.0)
!	   TX0=SIGN*TX0
!	ENDIF
C
!    2 CONTINUE
!      TS=TX0
!      DT=TS*1.E-8
!      TSD=TS+DT
C
!    1 CONTINUE
C
!      CALL SAT(TS,PS)
!      IF(IGOOD.NE.0) RETURN
!      IF(ABS((PX-PS)/PX).LE.1.E-10) RETURN
C
!      TSD=TS+DT
!      CALL SAT(TSD,PSD)
!      TS=TS+(PX-PS)*DT/(PSD-PS)
!      GOTO 1
C
!      END SUBROUTINE TSAT
cc       

!--------------------------------------------------------------------------------------

      subroutine TSAT(PX,TX00,TS)

      implicit none
    !     Arguments:
      real(DP),intent(in):: PX,TX00
      real(DP),intent(out)::TS
    !     Local variables:
      real(DP):: TX0,DT,PSD,CC,ACC,SIGN,PS,RelTol
      integer(I4B):: i
      real(DP):: Residual
      real(DP),parameter:: PTol=1.0d-10 ! Pressure tolerance for Newton iteration
      integer(I4B),parameter:: MaxIt=15  ! Max iterations for Newton
!      real(DP),parameter:: MinP=0.5d5 ! MinP,MaxP are the pressure bounds
!      real(DP),parameter:: MaxP=1000.d5
      real(DP),parameter:: MinP=1.d3 ! MinP,MaxP are the pressure bounds
      real(DP),parameter:: MaxP=2212.d4

!     FIND SATURATION TEMPERATURE TS AT PRESSURE PX.
!     TX0 IS STARTING TEMPERATURE FOR ITERATION.
!     Modified (AC 8/00)- screens P values and checks Newton iteration count.

!     Screen input pressure:
      if ((PX<MinP).or.(PX>MaxP)) then
      ! print *, 'IGOOD 8'
      ! print *, PX, MinP, MaxP
	  IGOOD=2
	else
        TX0=TX00
        if (TX0==0.0d0) then
!         Estimate starting temperature for iteration:
          if (PX>=656.613d0) then
            TX0=4606.d0/(24.02d0-LOG(PX))-273.15d0
            TX0=MAX(TX0,5.0d0)
	    else
	      CC=PX-610.768d0
	      ACC=ABS(CC)
	      SIGN=CC/ACC
	      TX0=(ACC/3.0667d-02)**(1.d0/3.d0)
	      TX0=SIGN*TX0
	    end if
	  end if
!       Initialise:
        TS=TX0
        DT=TS*1.d-8
	  RelTol=PTol*PX
	  i=0
	  Residual=1.d0
!       Newton loop: 
        do while ((Residual>RelTol).and.(IGOOD==0))
          call SAT(TS,PS)
	    call SAT(TS+DT,PSD)
          TS=TS+(PX-PS)*DT/(PSD-PS)
	    Residual=abs(PX-PS)
	    i=i+1
	    if (i>MaxIt) then
	      ! print *, 'IGOOD 9'
	      IGOOD=2 ! Too many iterations
	    end if
        end do

	end if

      end subroutine TSAT
!       
cc----7--1--------2---------3---------4---------5---------6---------7---------8 
        
      subroutine SIGMA(T,ST)

	implicit none
	real(DP),intent(in):: T
	real(DP),intent(out):: ST

C-----COMPUTE SURFACE TENSION OF WATER, USING THE
C     "INTERNATIONAL REPRESENTATION OF THE SURFACE TENSION OF
C      WATER SUBSTANCE" (1975).

	if (T<374.15) then
        ST=1.-0.625*(374.15-T)/647.3
        ST=ST*.2358*((374.15-T)/647.3)**1.256
	else
	  ST=0.0
	end if
      return
      end subroutine SIGMA
       
cc----7--1--------2---------3---------4---------5---------6---------7---------8 
cc        
      SUBROUTINE VIS(T,P,D,VW,VS,PS)
C
      implicit none
	real(DP):: T,P,D,VW,VS,PS
	real(DP):: EX,PHI,AM,V1       
C
      EX=247.8/(T+133.15)
      PHI=1.0467*(T-31.85)
      AM=1.+PHI*(P-PS)*1.E-11
      VW=1.E-7*AM*241.4*10.**EX
C
      V1=.407*T+80.4
      if (T<=350.0) then
	  VS=1.E-7*(V1-D*(1858.-5.9*T)*1.E-3)
	else
        VS=1.E-7*(V1+.353*D+676.5E-6*D**2+102.1E-9*D**3)
	end if

      RETURN
      END SUBROUTINE VIS
cc       
cc----7--1--------2---------3---------4---------5---------6---------7---------8 
cc        
      SUBROUTINE VISW(T,P,PS,VW)
C
      implicit none
	real(DP):: T,P,PS,VW
	real(DP):: EX,PHI,AM
	      
      EX=247.8/(T+133.15)
      PHI=1.0467*(T-31.85)
      AM=1.+PHI*(P-PS)*1.E-11
      VW=1.E-7*AM*241.4*10.**EX
C
      RETURN
      END SUBROUTINE VISW
cc       
cc----7--1--------2---------3---------4---------5---------6---------7---------8 
cc        
      SUBROUTINE VISS(T,P,D,VS)
C
      implicit none
	real(DP):: T,P,D,VS
	real(DP):: V1      
C
      V1=.407*T+80.4
      if (T<=350.0) then
	  VS=1.E-7*(V1-D*(1858.-5.9*T)*1.E-3)
      else
        VS=1.E-7*(V1+.353*D+676.5E-6*D**2+102.1E-9*D**3)
	end if
C
      RETURN
      END SUBROUTINE VISS
cc       
cc----7--1--------2---------3---------4---------5---------6---------7---------8 
cc        
      SUBROUTINE THERC(T,P,D,CONW,CONS,PS)
C
      implicit none
	real(DP):: T,P,D,CONW,CONS,PS
	real(DP):: A0,A1,A2,A3,A4,B0,B1,B2,B3,C0,C1,C2,C3,T0       
      real(DP):: T1,T2,T3,T4,CON1,CON2,CON3,CONS1

!      SAVE ICALL,A0,A1,A2,A3,A4,B0,B1,B2,B3,C0,C1,C2,C3,T0
      DATA A0,A1,A2,A3,A4/-922.47,2839.5,-1800.7,525.77,-73.440/
      DATA B0,B1,B2,B3/-.94730,2.5186,-2.0012,.51536/
      DATA C0,C1,C2,C3/1.6563E-3,-3.8929E-3,2.9323E-3,-7.1693E-4/
      DATA T0/273.15/
C
!      DATA ICALL/0/
!      ICALL=ICALL+1
!      AC: removed write statement for use of routine in DLL:
!      IF(ICALL.EQ.1) WRITE(11,899)
!  899 FORMAT(6X,'THERC    1.0       4 MARCH     1991',6X,
!     X'THERMAL CONDUCTIVITY OF WATER AND VAPOR AS FUNCTION OF',
!     X' TEMPERATURE AND PRESSURE')
C
      T1=(T+273.15)/T0
      T2=T1*T1
      T3=T2*T1
      T4=T3*T1
C
C     IF(P-PS.LT.0.) GOTO1
      CON1=A0+A1*T1+A2*T2+A3*T3+A4*T4
      CON2=(P-PS)*(B0+B1*T1+B2*T2+B3*T3)*1.E-5
      CON3=(P-PS)*(P-PS)*(C0+C1*T1+C2*T2+C3*T3)*1.E-10
      CONW=(CON1+CON2+CON3)*1.E-3
      CON1=17.6+5.87E-2*T
      CON2=1.04E-4*T*T
      CON3=4.51E-8*T**3
      CONS1=1.E-3*(CON1+CON2-CON3)
      CONS=CONS1+1.E-6*(103.51+.4198*T-2.771E-5*T*T)*D
     A+1.E-9*D*D*2.1482E14/T**4.2
C
C     PRINT 10,T,P,PS,CON
C   10 FORMAT(5H T = ,E12.6,5H P = ,E12.6,6H PS = ,E12.6,7H CON = ,E12.6)
C
!      RETURN
!    1 CONTINUE
C     PRINT 2,T,P,PS
!    2 FORMAT(8H AT T = ,E12.6,5H P = ,E12.6,19H IS LESS THAN PS = ,E12.6
!     A)
      RETURN
      END SUBROUTINE THERC
cc       
cc----7--1--------2---------3---------4---------5---------6---------7---------8 
cc
      SUBROUTINE THERMO(PX,TX,SVX,IPHX,RHOL,RHOV,HL,HV,VISL,VISV)

      implicit none
	real(DP):: PX,TX,SVX,RHOL,RHOV,HL,HV,VISL,VISV
	integer(I4B):: IPHX
	real(DP):: PS,DW,UW,DVISL,DS,US,DVISV
	 
	select case (IPHX)  
	case(0)    
cc       single phase liquid water
         CALL SAT(TX,PS)
         CALL COWAT(TX,PX,DW,UW)
	   if (IGOOD>0) then
	!     print *, 'IGOOD 10'
	    return
	   end if
         RHOL=DW
         HL=UW+PX/RHOL
         RHOV=0.0
         HV=0.0
         CALL VISW(TX,PX,PS,DVISL)
	   VISL=DVISL/RHOL
         VISV=0.0
	case(1)
cc       2-phase 
         CALL COWAT(TX,PX,DW,UW)
	   if (IGOOD>0) then
	!     print *, 'IGOOD 11'
	    return
	   end if
         CALL SUPST(TX,PX,DS,US)
         RHOL=DW
         HL=UW+PX/RHOL
         RHOV=DS
         HV=US+PX/RHOV
         CALL VIS(TX,PX,DS,DVISL,DVISV,PX)
	   VISL=DVISL/RHOL
	   VISV=DVISV/RHOV
	case(2)
cc       single phase dry steam
         CALL SUPST(TX,PX,DS,US)
         RHOL=0.0
         HL=0.0
         RHOV=DS
         HV=US+PX/RHOV
         CALL VISS(TX,PX,DS,DVISV)
         VISL=0.0
	   VISV=DVISV/RHOV
      end select
      RETURN
      END SUBROUTINE THERMO      
cc       

cc----7--1--------2---------3---------4---------5---------6---------7---------8 
cc
      SUBROUTINE RELP(SV,IRELP,XKRL,XKRV)

      implicit none
	real(DP)::  SV,XKRL,XKRV
	integer(I4B):: IRELP
	real(DP):: SL,SLX    
	real(DP), parameter:: SLR0=0.55, SVR0=0.45
	
	SL=1.0D0-SV
      SLX=SL-SLR0
	if (SLX<=0.0D0) then
	   XKRL=0.0D0
	else
	   XKRL=SLX/(1.0D0-SLR0)
	end if
	if (SV<SVR0) THEN
	   XKRV=SV/SVR0
	else
	   XKRV=1.0D0
	end if

	RETURN
	END SUBROUTINE RELP
cc       

!------------------------------------------------------------------------------

      SUBROUTINE PERF(POR0X,PORX,AAA,PERFAC)  
!     Calculates permeability update factor PERFAC, based on the change
!     in porosity, for temperature- and pressure-dependent permeability.
      
	implicit none
	real(DP):: POR0X,PORX,AAA,PERFAC
	real(DP):: R1,CC1,R2,CC2,CC3

      R1=PORX/POR0X
      CC1=R1*R1*R1
      R2=(1.0D0-POR0X)/(1.0D0-PORX)
      CC2=R2*R2
      CC3=DEXP(AAA*(PORX-POR0X))
      PERFAC=CC1*CC2*CC3

      RETURN
      end subroutine PERF

!----------------------------------------------------------------------------------

	end module Thermodynamics
