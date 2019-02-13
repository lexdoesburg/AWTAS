
      module NumericalSimulator1D_Routines

	use Thermodynamics
	use MatrixSolvers
	use problem_data

	contains

cc----7--1--------2---------3---------4---------5---------6---------7---------8 
cc        
      SUBROUTINE SOLVE(P,T,SV,X,IPH,HF,M,POR,PER,A,V,CR,RHOR,COND,
     1             QQMM,HIN,XX,RR,DT,DELR,BMOLD,BEOLD,RMAX,
     2             P0,T0,COMP,COMT,AAA,QMM,XLAM,PRECH,HRECH,THICK)
!       New version for temp/pressure-dependent porosity & permeability (11/99).
!       Deliverability added 8/00.
!       Variable compressibility added 4/02.
 
      IMPLICIT REAL*8 (A-H,O-Z)       
      DIMENSION P(M),T(M),SV(M),X(M)
      DIMENSION A(M),V(M),DELR(M),COMP(M)
	DIMENSION C1(2*M),C2(2*M),C3(2*M),C4(2*M)   
	DIMENSION C5(2*M),C6(2*M),C7(2*M)
      DIMENSION XX(2*M),R(2*M),RR(2*M)
      DIMENSION IPH(M)
	DIMENSION POR(M),PER(M)
	DIMENSION BMOLD(M),BEOLD(M)
      DIMENSION P0(M),T0(M)
      ! COMMON/WRIT/IWRIT,IT
	! CHARACTER*6 NAME
      DATA FACP,FACT,FACS,EnergyScale /1.0E-8,1.0E-8,1.0E-8,1.0E-6/
	real*8, parameter:: small=1.0D-5

cc
cc  SOLVE calculates the increments to the solution for one step of the Newton-Raphson method
cc
cc  The labelling syatem below is used for the blocks and block boundaries
cc
cc             I-1             |            I            |            I+1
cc                             |                         |
cc              *              |            *            |             *
cc                             |                         |
cc              A              |            B            |             C
cc                             M                         P
cc  
cc  pressure            PA,PB,PC
cc  temperature         TA,TB,TC 
cc  2nd variable        XA,XB,XC  (X = T or SV)
cc  phase indicator     IPHA,IPHB,IPHC  (IPH = 0 water, = 1 2-phase, = 2 steam)
cc  saturation          SVA,SVB,SVC
cc  mass accumulation   BMB,BMC
cc  energy accumulation BEB,BEC 
cc  mass flux           FMM,FMP
cc  energy flux         FEM,FEP
cc
cc  Now start the construction of the Jacobian matrix and the residuals                                                
cc
cc  The Jacobian matrix J has a maximum of 6 non-zero elements in each row.
cc  The equations and unknowns are ordered mass-energy alternately
cc  For each mass equation, index J, the coefficients are
cc  C1(J) = 0.0
cc  C2(J) = DRM/DPA
cc  C3(J) = DRM/DXA
cc  C4(J) = DRM/DPB
cc  C5(J) = DRM/DXB
cc  C6(J) = DRM/DPC
cc  C7(J) = DRM/DXC
cc  For each energy equation, index J1, the coefficients are
cc  C1(J1) = DRE/DPA
cc  C2(J1) = DRE/DXA
cc  C3(J1) = DRE/DPB
cc  C4(J1) = DRE/DXB
cc  C5(J1) = DRE/DPC
cc  C6(J1) = DRE/DXC
cc  C7(J1) = 0.0
cc  derivatives for mass equation
cc  Note that the equations for the 1st and last blocks are special
cc  
cc  block 1 (the well block)
cc
      RMAX=0.0D0
      PB=P(1)
      TB=T(1)
      XB=X(1)
      SVB=SV(1)
      SLB=1.0-SVB
      IPHB=IPH(1)
      VB=V(1)
      PORB=POR(1)
      PERB=PER(1)
      P0B=P0(1)
      T0B=T0(1)
      COMPB=COMP(1)
cc
      CALL TERMS(PB,TB,SVB,IPHB,PB1,TB1,TB2,
     1             BMB,BEB,TMB,TEB,BMB1,BEB1,TMB1,TEB1,
     2             BMB2,BEB2,TMB2,TEB2,DELPB,DELXB,
     3             PORB,CR,RHOR,FACP,FACT,FACS,
     4             PERB,P0B,T0B,COMPB,COMT,AAA)
cc 
cc  TERMS calculates accumulation terms, and transmissibility terms for PB,TB,SVB and 
cc  incremented values of variables. 
cc  "1" indicates values for PB1=PB+DELPB, TB or SVB are not changed, DELPB is calculated
cc  by TERMS.
cc  "2" indicates values for TB2=TB+DELXB (or SVB2=SVB+DELXB), PB is not changed, DELXB is 
cc  calculated by TERMS.
cc 
	 IF(IGOOD.EQ.2) RETURN
cc  now calculate derivatives
cc  accumulation terms
      BMB_PB=(BMB1-BMB)/DELPB
      BMB_XB=(BMB2-BMB)/DELXB
      BEB_PB=(BEB1-BEB)/DELPB
      BEB_XB=(BEB2-BEB)/DELXB
cc  production/injection term
      if (OnDeliv) then ! Pump on deliverability:
        if(QQMM>=0.0D0) then
          QMM=QQMM
          QEM=QMM*HIN
          QMM_PB=0.0
          QMM_XB=0.0
          QEM_PB=0.0
          QEM_XB=0.0
        else ! Production:
          if(PB<=PCutoff) then
            QMM=0.0
            QEM=0.0
            QMM_PB=0.0
            QMM_XB=0.0
            QEM_PB=0.0
            QEM_XB=0.0
          else
            QQMMX=-(TMB/PERB)*ProdIndex*(PB-PCutoff)
	      AQQMMX=ABS(QQMMX)
	      AQQMM=ABS(QQMM)
            if(AQQMMX<AQQMM) then
               QMM=QQMMX
               QEM=QMM*TEB/TMB
               QMM1=-(TMB1/PERB)*ProdIndex*(PB1-PCutoff)
               QMM2=-(TMB2/PERB)*ProdIndex*(PB-PCutoff)
               QMM_PB=(QMM1-QMM)/DELPB
               QMM_XB=(QMM2-QMM)/DELXB
               QEM1=QMM1*TEB1/TMB1
               QEM2=QMM2*TEB2/TMB2
               QEM_PB=(QEM1 - QEM)/DELPB
               QEM_XB=(QEM2 - QEM)/DELXB
            else 
               QMM=QQMM          
               QMM_PB=0.0
               QMM_XB=0.0
               QEM=QMM*TEB/TMB
               QEM1=QMM*TEB1/TMB1
               QEM2=QMM*TEB2/TMB2
               QEM_PB=(QEM1 - QEM)/DELPB
               QEM_XB=(QEM2 - QEM)/DELXB
            end if
          end if
        end if
      else ! No deliverability:
	  QMM=QQMM
        if(QMM>=0.0) then
          QEM=QMM*HIN
          QEM_PB=0.0
          QEM_XB=0.0
        else
          QEM=QMM*TEB/TMB
          QEM_PB=(QMM*TEB1/TMB1 - QEM)/DELPB
          QEM_XB=(QMM*TEB2/TMB2 - QEM)/DELXB
        end if
	end if
c     Flowing enthalpy calculation (AC 20/9/99):
      if (dabs(QMM)>small) then
	  HF=QEM/QMM
	else
	  HF=0.0D0
	end if
cc
cc  recharge terms (added AC 15/12/00):
cc
      QMRB=XLAM*(PB-PRECH)
      QMRB1=XLAM*(PB1-PRECH)
      QMRB_PB=XLAM
      QMRB_XB=0.0
      if (QMRB.LT.0.0) then
         QERB=QMRB*HRECH
         QERB_PB=XLAM*HRECH
         QERB_XB=0.0
      else
         QERB=QMRB*TEB/TMB
         QERB_PB=(QMRB*TEB1/TMB1 - QERB)/DELPB
         QERB_XB=(QMRB*TEB2/TMB2 - QERB)/DELXB
      end if  
cc
cc  block 2
cc
      PC=P(2)
      TC=T(2)
      SVC=SV(2)
      XC=X(2)
      SLC=1.0-SVC
      IPHC=IPH(2)
      VC=V(2)
      PORC=POR(2)
      PERC=PER(2)
      P0C=P0(2)
      T0C=T0(2)
      AP=A(1)
      DRP=DELR(1)
	COMPC=COMP(2)
cc
      CALL TERMS(PC,TC,SVC,IPHC,PC1,TC1,TC2,
     1             BMC,BEC,TMC,TEC,BMC1,BEC1,TMC1,TEC1,
     2             BMC2,BEC2,TMC2,TEC2,DELPC,DELXC,
     3             PORC,CR,RHOR,FACP,FACT,FACS,
     4             PERC,P0C,T0C,COMPC,COMT,AAA)
cc
	 IF(IGOOD.EQ.2) RETURN
      BMC_PC=(BMC1-BMC)/DELPC
      BMC_XC=(BMC2-BMC)/DELXC
      BEC_PC=(BEC1-BEC)/DELPC
      BEC_XC=(BEC2-BEC)/DELXC
cc

cc  recharge terms at C block (added AC 15/12/00):
cc
      QMRC=XLAM*(PC-PRECH)
      QMRC1=XLAM*(PC1-PRECH)
      QMRC_PC=XLAM
      QMRC_XC=0.0
      if (QMRC.LT.0.0) then
         QERC=QMRC*HRECH
         QERC_PC=XLAM*HRECH
         QERC_XC=0.0
      else
         QERC=QMRC*TEC/TMC
         QERC_PC=(QMRC*TEC1/TMC1 - QERC)/DELPC
         QERC_XC=(QMRC*TEC2/TMC2 - QERC)/DELXC
      end if  
cc
cc  calculate flux terms using upstream weighting
cc  1st PC > PB
cc
      IF(PC>PB) THEN
cc  basic flows
         FMP=- TMC*(PC-PB)/DRP
         FEP=- TEC*(PC-PB)/DRP - COND*(TC-TB)/DRP
cc  flows with incremented pressure in block C
         FMP1C=- TMC1*(PC1-PB)/DRP
         FEP1C=- TEC1*(PC1-PB)/DRP - COND*(TC1-TB)/DRP
cc  flows with incremented temperature or saturation in block C
         FMP2C=- TMC2*(PC-PB)/DRP
         FEP2C=- TEC2*(PC-PB)/DRP - COND*(TC2-TB)/DRP
cc  flows with incremented pressure in block B
	 FMP1B=- TMC*(PC-PB1)/DRP
	 FEP1B=- TEC*(PC-PB1)/DRP - COND*(TC-TB1)/DRP
cc  flows with incremented temperature or saturation in block B
	 FEP2B=- TEC*(PC-PB)/DRP - COND*(TC-TB2)/DRP
cc  mass flux derivatives
        FMP_PC=(FMP1C - FMP)/DELPC
        FMP_XC=(FMP2C - FMP)/DELXC
        FMP_PB=(FMP1B - FMP)/DELPB
        FMP_XB=0.0D0
cc  energy flux derivatives
        FEP_PC=(FEP1C - FEP)/DELPC
        FEP_XC=(FEP2C - FEP)/DELXC
        FEP_PB=(FEP1B - FEP)/DELPB
        FEP_XB=(FEP2B - FEP)/DELXB
	HFFF=TEC/TMC
cc
cc  or PC < PB or = PB
cc
      else
cc  fluxes
cc  basic flows
        FMP=- TMB*(PC-PB)/DRP
        FEP=- TEB*(PC-PB)/DRP - COND*(TC-TB)/DRP
cc  flows with incremented pressure in block C
	  FMP1C=- TMB*(PC1-PB)/DRP
	  FEP1C=- TEB*(PC1-PB)/DRP - COND*(TC1-TB)/DRP
cc  flows with incremented temperature or saturation in block C
	  FEP2C=- TEB*(PC-PB)/DRP - COND*(TC2-TB)/DRP
cc  flows with incremented pressure in block B
	  FMP1B=- TMB1*(PC-PB1)/DRP
	  FEP1B=- TEB1*(PC-PB1)/DRP - COND*(TC-TB1)/DRP 
cc  flows with incremented temperature or saturation in block B
	  FMP2B=- TMB2*(PC-PB)/DRP
	  FEP2B=- TEB2*(PC-PB)/DRP - COND*(TC-TB2)/DRP
cc  mass flux derivatives
        FMP_PC=(FMP1C - FMP)/DELPC
        FMP_XC=0.0D0
        FMP_PB=(FMP1B - FMP)/DELPB
        FMP_XB=(FMP2B - FMP)/DELXB
cc  energy flux derivatives
        FEP_PC=(FEP1C - FEP)/DELPC
        FEP_XC=(FEP2C - FEP)/DELXC
        FEP_PB=(FEP1B - FEP)/DELPB
        FEP_XB=(FEP2B - FEP)/DELXB
	HFFF=TEB/TMB
      end if
	I=1
	FFMP=FMP*AP
	FFEP=FEP*AP
	FCON=- AP*COND*(TC-TB)/DRP
cc	WRITE(16,333) I,FFMP,FFEP,FCON,HFFF
cc333   FORMAT(1H ,I5,4(2X,E12.5))
cc
cc  set up coefficients in Jacobian
cc
      J=1
      J1=2
cc  mass residual
cc
cc  mass residual is given by
cc  R(I)=(BMB-BMOLD(I)) + (FMP*AP - QMM)*DT/VB
cc  Note that the -ve sign is include below. We solve [J]XX=-R
cc
      CCB=DT/VB
	CCD=DT
      ! CCD=DT/THICK ! (Added AC 15/12/00, + added terms to next 2 lines:)	
      R(J)=-(BMB-BMOLD(1)+(FMP*AP-QMM)*CCB)-QMRB*CCD
      RR(J)=-(BMB-BMOLD(1)+(FMP*AP-QMM)*CCB)-QMRB*CCD
      RRR=ABS(R(J))
      RMAX=MAX(RMAX,RRR)
      C1(J)=0.0D0
      C2(J)=0.0D0
      C3(J)=0.0D0
      C4(J)=BMB_PB
      C5(J)=BMB_XB
	! Added terms below for deliverability (QMM_PB, QMM_XB) and recharge:
      C4(J)=C4(J)+(AP*FMP_PB-QMM_PB)*CCB+QMRB_PB*CCD 
      C5(J)=C5(J)+(AP*FMP_XB-QMM_XB)*CCB+QMRB_XB*CCD
      C6(J)=AP*FMP_PC*CCB
      C7(J)=AP*FMP_XC*CCB
cc
cc  energy residual
      ! Added QERB*CCD terms to next 2 lines for recharge:
      R(J1)=-EnergyScale*(BEB-BEOLD(1)+(FEP*AP-QEM)*CCB+QERB*CCD)
      RR(J1)=-EnergyScale*(BEB-BEOLD(1)+(FEP*AP-QEM)*CCB+QERB*CCD)
      RRR=ABS(R(J1))
      RMAX=MAX(RMAX,RRR)
cc  derivatives for energy equation
      C1(J1)=0.0D0
      C2(J1)=0.0D0
      C7(J1)=0.0D0
      C3(J1)=BEB_PB*EnergyScale
      C4(J1)=BEB_XB*EnergyScale
      C3(J1)=C3(J1)+(AP*FEP_PB*CCB)*EnergyScale
      C4(J1)=C4(J1)+(AP*FEP_XB*CCB)*EnergyScale
      C5(J1)=(AP*FEP_PC*CCB)*EnergyScale
      C6(J1)=(AP*FEP_XC*CCB)*EnergyScale
      C3(J1)=C3(J1)-(QEM_PB*CCB)*EnergyScale+QERB_PB*CCD*EnergyScale
      C4(J1)=C4(J1)-(QEM_XB*CCB)*EnergyScale+QERB_XB*CCD*EnergyScale
cc
cc  iterate over internal blocks
      do I=2,M-1
cc  shuffle variables along
      IC=I+1
      PA=PB
      TA=TB
      XA=XB
      SVA=SVB
      SLA=SLB
      IPHA=IPHB
      VA=VB
      PB=PC
      PB1=PC1
	TB1=TC1
      TB=TC
      TB2=TC2
      XB=XC
      SVB=SVC
      SLB=SLC
      IPHB=IPHC
      VB=VC
      BMB=BMC
      BEB=BEC
	! Added for recharge (AC 15/12/00):
      QMRB=QMRC
      QMRB_PB=QMRC_PC
      QMRB_XB=QMRC_XC
      QERB=QMRC
      QERB_PB=QERC_PC
      QERB_XB=QERC_XC
      ! End of recharge section

      DELPB=DELPC
      DELXB=DELXC
      AM=AP
      DRM=DRP
	TMB=TMC
	TEB=TEC
	TMB1=TMC1
	TEB1=TEC1
	TMB2=TMC2
	TEB2=TEC2

	COMPA=COMPB
	COMPB=COMPC

cc  shuffle along derived quantities
      BMB_PB=BMC_PC 
      BMB_XB=BMC_XC     
      BEB_PB=BEC_PC     
      BEB_XB=BEC_XC
      FMM=FMP
      FEM=FEP
      FMM_PA=FMP_PB
      FMM_XA=FMP_XB     
      FMM_PB=FMP_PC
      FMM_XB=FMP_XC     
      FEM_PA=FEP_PB
      FEM_XA=FEP_XB     
      FEM_PB=FEP_PC
      FEM_XB=FEP_XC     
      PC=P(IC)
      TC=T(IC)
      XC=X(IC)
      SVC=SV(IC)
      SLC=1.0-SVC
      IPHC=IPH(IC)
      VC=V(IC)
      PORC=POR(IC)
      PERC=PER(IC)
      P0C=P0(IC)
      T0C=T0(IC)
      AP=A(I)
      DRP=DELR(I)

	COMPC=COMP(IC)
cc
       CALL TERMS(PC,TC,SVC,IPHC,PC1,TC1,TC2,
     1            BMC,BEC,TMC,TEC,BMC1,BEC1,TMC1,TEC1,
     2            BMC2,BEC2,TMC2,TEC2,DELPC,DELXC,
     3            PORC,CR,RHOR,FACP,FACT,FACS,
     4            PERC,P0C,T0C,COMPC,COMT,AAA)
cc
	 IF(IGOOD.EQ.2) RETURN
cc  now calculate derivatives
cc  accumulation terms
      BMC_PC=(BMC1-BMC)/DELPC
      BMC_XC=(BMC2-BMC)/DELXC
      BEC_PC=(BEC1-BEC)/DELPC
      BEC_XC=(BEC2-BEC)/DELXC

cc  recharge terms at C block (added AC 15/12/00):
cc
      QMRC=XLAM*(PC-PRECH)
      QMRC1=XLAM*(PC1-PRECH)
      QMRC_PC=XLAM
      QMRC_XC=0.0
      if (QMRC.LT.0.0) then
         QERC=QMRC*HRECH
         QERC_PC=XLAM*HRECH
         QERC_XC=0.0
      else
         QERC=QMRC*TEC/TMC
         QERC_PC=(QMRC*TEC1/TMC1 - QERC)/DELPC
         QERC_XC=(QMRC*TEC2/TMC2 - QERC)/DELXC
      end if  
cc

cc
cc  calculate flux terms using upstream weighting
cc  1st PC > PB
cc
      IF(PC>PB) THEN
cc  basic flows
        FMP=- TMC*(PC-PB)/DRP
        FEP=- TEC*(PC-PB)/DRP - COND*(TC-TB)/DRP
cc  flows with incremented pressure in block C
        FMP1C=- TMC1*(PC1-PB)/DRP
        FEP1C=- TEC1*(PC1-PB)/DRP - COND*(TC1-TB)/DRP
cc  flows with incremented temperature or saturation in block C
        FMP2C=- TMC2*(PC-PB)/DRP
        FEP2C=- TEC2*(PC-PB)/DRP - COND*(TC2-TB)/DRP
cc  flows with incremented pressure in block B
	FMP1B=- TMC*(PC-PB1)/DRP
	FEP1B=- TEC*(PC-PB1)/DRP - COND*(TC-TB1)/DRP
cc  flows with incremented temperature or saturation in block B
	FEP2B=- TEC*(PC-PB)/DRP - COND*(TC-TB2)/DRP
cc  mass flux derivatives
        FMP_PC=(FMP1C - FMP)/DELPC
        FMP_XC=(FMP2C - FMP)/DELXC
        FMP_PB=(FMP1B - FMP)/DELPB
        FMP_XB=0.0D0
cc  energy flux derivatives
        FEP_PC=(FEP1C - FEP)/DELPC
        FEP_XC=(FEP2C - FEP)/DELXC
        FEP_PB=(FEP1B - FEP)/DELPB
        FEP_XB=(FEP2B - FEP)/DELXB
	HFFF=TEC/TMC
      else
cc
cc  basic flows
        FMP=- TMB*(PC-PB)/DRP
        FEP=- TEB*(PC-PB)/DRP - COND*(TC-TB)/DRP
cc  flows with incremented pressure in block C
	  FMP1C=- TMB*(PC1-PB)/DRP
	  FEP1C=- TEB*(PC1-PB)/DRP - COND*(TC1-TB)/DRP
cc  flows with incremented temperature or saturation in block C
	  FEP2C=- TEB*(PC-PB)/DRP - COND*(TC2-TB)/DRP
cc  flows with incremented pressure in block B
	  FMP1B=- TMB1*(PC-PB1)/DRP
	  FEP1B=- TEB1*(PC-PB1)/DRP - COND*(TC-TB1)/DRP
cc  flows with incremented temperature or saturation in block B
	  FMP2B=- TMB2*(PC-PB)/DRP
	  FEP2B=- TEB2*(PC-PB)/DRP - COND*(TC-TB2)/DRP
cc  mass flux derivatives
        FMP_PC=(FMP1C - FMP)/DELPC
        FMP_XC=0.0D0
        FMP_PB=(FMP1B - FMP)/DELPB
        FMP_XB=(FMP2B - FMP)/DELXB
cc  energy flux derivatives
        FEP_PC=(FEP1C - FEP)/DELPC
        FEP_XC=(FEP2C - FEP)/DELXC
        FEP_PB=(FEP1B - FEP)/DELPB
        FEP_XB=(FEP2B - FEP)/DELXB
	HFFF=TEB/TMB
      end if
	FFMP=FMP*AP
	FFEP=FEP*AP
	FCON=- AP*COND*(TC-TB)/DRP
cc	WRITE(16,333) I,FFMP,FFEP,FCON,HFFF
cc
cc  set up coefficients in Jacobian
cc
      J=2*I-1
      J1=J+1
cc  mass residual
      CCB=DT/VB
	! Recharge terms added (AC 15/12/00):
      R(J)=-(BMB-BMOLD(I)+(FMP*AP-FMM*AM)*CCB)-QMRB*CCD
      RR(J)=-(BMB-BMOLD(I)+(FMP*AP-FMM*AM)*CCB)-QMRB*CCD
	RRR=ABS(R(J))
	RMAX=MAX(RMAX,RRR)
cc  derivatives for mass equation
      C1(J)=0.0D0
      C4(J)=BMB_PB
      C5(J)=BMB_XB
      C4(J)=C4(J)+AP*FMP_PB*CCB
      C5(J)=C5(J)+AP*FMP_XB*CCB
      C6(J)=AP*FMP_PC*CCB
      C7(J)=AP*FMP_XC*CCB
      C2(J)=-AM*FMM_PA*CCB
      C3(J)=-AM*FMM_XA*CCB
	! Recharge terms added (AC 15/12/00):
      C4(J)=C4(J)-AM*FMM_PB*CCB+CCD*QMRB_PB
      C5(J)=C5(J)-AM*FMM_XB*CCB+CCD*QMRB_XB
cc  energy residual (recharge terms added AC 15/12/00): 
      R(J1)=-(BEB-BEOLD(I)+(FEP*AP-FEM*AM)*CCB+QERB*CCD)*EnergyScale
      RR(J1)=-(BEB-BEOLD(I)+(FEP*AP-FEM*AM)*CCB+QERB*CCD)*EnergyScale
	RRR=ABS(R(J1))
	RMAX=MAX(RMAX,RRR)
cc  derivatives for energy equation
      C7(J1)=0.0D0
      C3(J1)=BEB_PB*EnergyScale
      C4(J1)=BEB_XB*EnergyScale
	! Recharge terms added (AC 15/12/00):
      C3(J1)=C3(J1)+(AP*FEP_PB*CCB+CCD*QERB_PB)*EnergyScale
      C4(J1)=C4(J1)+(AP*FEP_XB*CCB+CCD*QERB_XB)*EnergyScale
      C5(J1)=(AP*FEP_PC*CCB)*EnergyScale
      C6(J1)=(AP*FEP_XC*CCB)*EnergyScale
      C1(J1)=(-AM*FEM_PA*CCB)*EnergyScale
      C2(J1)=(-AM*FEM_XA*CCB)*EnergyScale
      C3(J1)=C3(J1)-(AM*FEM_PB*CCB)*EnergyScale
      C4(J1)=C4(J1)-(AM*FEM_XB*CCB)*EnergyScale
      end do
cc
cc set up the i=m equation
cc
cc  shuffle variables along
      PA=PB
      TA=TB
	XA=XB
      SVA=SVB
      SLA=SLB
      IPHA=IPHB
      VA=VB
      PB=PC
	PB1=PC1
	TB1=TC1
      TB=TC
	TB2=TC2
	XB=XC
      SVB=SVC
      SLB=SLC
      IPHB=IPHC
      VB=VC
	BMB=BMC
	BEB=BEC
	! Recharge terms added (AC 15/12/00):
      QMRB=QMRC
      QMRB_PB=QMRC_PC
      QMRB_XB=QMRC_XC
      QERB=QERC
      QERB_PB=QERC_PC
      QERB_XB=QERC_XC
      ! End of recharge terms
      FMM=FMP
      FEM=FEP
      AM=AP

cc  shuffle along derived quanitities
      BMB_PB=BMC_PC 
      BMB_XB=BMC_XC     
      BEB_PB=BEC_PC     
      BEB_XB=BEC_XC
      FMM=FMP
      FEM=FEP
      FMM_PA=FMP_PB
      FMM_XA=FMP_XB     
      FMM_PB=FMP_PC
      FMM_XB=FMP_XC     
      FEM_PA=FEP_PB
      FEM_XA=FEP_XB     
      FEM_PB=FEP_PC
      FEM_XB=FEP_XC     
cc  set up coefficients in Jacobian
cc
      CCB=DT/VB
      J=2*M-1
      J1=J+1
cc  mass residual (recharge terms added AC 15/12/00):
      R(J)=-(BMB-BMOLD(M) - FMM*AM*CCB)-CCD*QMRB
      RR(J)=-(BMB-BMOLD(M) - FMM*AM*CCB)-CCD*QMRB
	RRR=ABS(R(J))
	RMAX=MAX(RMAX,RRR)
cc  derivatives for mass equation
      C1(J)=0.0D0
      C2(J)=0.0D0
      C3(J)=0.0D0
      C4(J)=0.0D0
      C5(J)=0.0D0
      C6(J)=0.0D0
      C7(J)=0.0D0
      C4(J)=BMB_PB
      C5(J)=BMB_XB
      C2(J)=-AM*FMM_PA*CCB
      C3(J)=-AM*FMM_XA*CCB
	! Recharge terms added (AC 15/12/00):
      C4(J)=C4(J)-AM*FMM_PB*CCB+CCD*QMRB_PB
      C5(J)=C5(J)-AM*FMM_XB*CCB+CCD*QMRB_XB
cc  energy residual (recharge terms added AC 15/12/00):
      R(J1)=-(BEB-BEOLD(M) - FEM*AM*CCB+CCD*QERB)*EnergyScale
      RR(J1)=-(BEB-BEOLD(M) - FEM*AM*CCB+CCD*QERB)*EnergyScale
	RRR=ABS(R(J1))
	RMAX=MAX(RMAX,RRR)
cc  derivatives for energy equation
      C1(J1)=0.0D0
      C2(J1)=0.0D0
      C3(J1)=0.0D0
      C4(J1)=0.0D0
      C5(J1)=0.0D0
      C6(J1)=0.0D0
      C7(J1)=0.0D0
	! Recharge terms added (AC 15/12/00): 
      C3(J1)=(BEB_PB+CCD*QERB_PB)*EnergyScale
      C4(J1)=(BEB_XB+CCD*QERB_XB)*EnergyScale
      C1(J1)=(-AM*FEM_PA*CCB)*EnergyScale
      C2(J1)=(-AM*FEM_XA*CCB)*EnergyScale
      C3(J1)=C3(J1)-(AM*FEM_PB*CCB)*EnergyScale
      C4(J1)=C4(J1)-(AM*FEM_XB*CCB)*EnergyScale
      M2=2*M

      CALL SEVEN1(C1,C2,C3,C4,C5,C6,C7,R,XX,M2)
      RETURN
      end subroutine SOLVE
cc       
cc----7--1--------2---------3---------4---------5---------6---------7---------8 
cc        
      SUBROUTINE UPDATE(M,IPH,P,T,SV,X,XX,EPS)
      IMPLICIT REAL*8 (A-H,O-Z)       
      DIMENSION P(M),T(M),X(M),SV(M),XX(2*M)
      DIMENSION IPH(M)
cc
cc  This subroutine adds the increments on to the primary variables and carries out
cc  phase changes
cc
      do I=1,M
         J=2*I-1
         J1=J+1
         PP=P(I)
         TP=T(I)
         IPHP=IPH(I)
         PT=PP+XX(J)
         XT=X(I)+XX(J1)
cc
cc  -----------------------------------------------------------------------------
         select case (IPHP)
	   case(0)
cc
cc  currently single phase liquid water
            TT=XT
            CALL SAT(TT,PS)
            if (PT>PS) then
cc  remain single phase liquid water
               P(I)=PT
               T(I)=TT
               X(I)=TT
               SV(I)=0.0
               IPH(I)=0
            else
cc  change to 2-phase
               P(I)=PS
               T(I)=TT
               SV(I)=EPS
               X(I)=EPS
               IPH(I)=1
            end if
cc
cc  -----------------------------------------------------------------------------
         case(1)
cc
cc  currently 2-phase
            SVT=XT
            CALL SAT(TP,PS)
            IF((SVT>0.0).AND.(SVT<1.0)) THEN
cc  remain 2-phase
               P(I)=PT
               SV(I)=SVT
               X(I)=SVT
               IPH(I)=1
	         T0=0.0
               CALL TSAT(PT,T0,TT)
	         if (IGOOD>0) return
               T(I)=TT
            ELSEIF(SVT<=0.0) THEN
cc  change to single phase liquid water
               P(I)=1.000001*PS
               T(I)=TP
               X(I)=TP
               SV(I)=0.0
               IPH(I)=0
            ELSEIF(SVT>=1.0) THEN
cc  change to single phase dry steam
               P(I)=0.999999*PS
               T(I)=TP
               X(I)=TP
               SV(I)=1.0
               IPH(I)=2
            ENDIF
cc
cc  -----------------------------------------------------------------------------
         case(2)
cc
cc  currently single phase dry steam
            TT=XT
            CALL SAT(TT,PS)
            if (PT<PS) then
cc  remain single phase dry steam
               P(I)=PT
               T(I)=TT
               X(I)=TT
               SV(I)=1.0
               IPH(I)=2
            else
cc  change to 2-phase
               P(I)=PS
               T(I)=TT
               SV(I)=1.0-EPS
               X(I)=1.0-EPS
               IPH(I)=1
            end if
	   end select
      end do
      RETURN
      end subroutine UPDATE
cc       
cc----7--1--------2---------3---------4---------5---------6---------7---------8 
cc
      SUBROUTINE INIT1(M,POLD,XOLD,TOLD,SVOLD,IPHOLD)
      IMPLICIT REAL*8 (A-H,O-Z)       
      DIMENSION POLD(M),XOLD(M),TOLD(M),SVOLD(M),IPHOLD(M)
      DIMENSION P(M),T(M),SV(M)
      do I=1,M
         PP=POLD(I)
         XP=XOLD(I)
         IF(I==1) THEN ! Added to fix issue of TZERO being unintialised
           TZERO=0.0 
         END IF
	   CALL TSAT(PP,TZERO,TS)
         IF(XP<1.0) THEN
cc  2-phase conditions
            TZERO=0.0
            TP=TS 
            TOLD(I)=TS
            SVOLD(I)=XP
            SVP=XP
            IPHOLD(I)=1
            IPHP=1
         ELSEIF(XP>1.5) THEN
            IF(XP>=TS) THEN
cc   single phase dry steam
               TP=XP
               TOLD(I)=XP
               SVOLD(I)=1.0
               SVP=1.0
               IPHOLD(I)=2
               IPHP=2
            else
cc   single phase liquid water
               TP=XP
               TOLD(I)=XP
               SVOLD(I)=0.0
               SVP=0.0
               IPHOLD(I)=0
               IPHP=0
            end if
         ENDIF
      end do
      RETURN
      end subroutine INIT1
cc       
cc----7--1--------2---------3---------4---------5---------6---------7---------8 
cc
      SUBROUTINE INIT2(M,POLD,XOLD,TOLD,SVOLD,IPHOLD,
     1                   BMOLD,BEOLD,POR,RHOR,CR,P0,T0,COMP,COMT,AAA)
      IMPLICIT REAL*8 (A-H,O-Z)       
      DIMENSION POLD(M),XOLD(M),TOLD(M),SVOLD(M),IPHOLD(M)
      DIMENSION BMOLD(M),BEOLD(M)
      DIMENSION POR(M),P0(M),T0(M),COMP(M)
      do I=1,M
         PP=POLD(I)
         TP=TOLD(I)
         P0P=P0(I)
         T0P=T0(I)
	       SVP=SVOLD(I)
	       PORP=POR(I)
         IPHP=IPHOLD(I)
	       COMPP=COMP(I)
         CALL THERMO(PP,TP,SVP,IPHP,RHOL,RHOV,HL,HV,VISL,VISV)
         CALL ACCUM(BM,BE,IPHP,PP,TP,SVP,RHOL,RHOV,HL,HV,PORP,RHOR,CR,
     1              P0P,T0P,COMPP,COMT,AAA,PERFAC)
         BMOLD(I)=BM
         BEOLD(I)=BE
      end do
      RETURN
      end subroutine INIT2
cc       
! -----------------------------------------------------------------------------------
      
	end module NumericalSimulator1D_Routines