      module NumericalSimulatorMultiLayer_Routines

	use Thermodynamics
	use MatrixSolvers

	contains

cc----7--1--------2---------3---------4---------5---------6---------7---------8 
cc        
      SUBROUTINE SOLVE(P,T,SV,X,IPH,HF,XX,BMOLD,BEOLD,MR,NZ,PORFR,PORMA,
     1           PERFR,PERMA,AR,AZ,V,CRF,CRM,RHORF,RHORM,CONDF,CONDM,
     2           QMM,HIN,DT,DELR,DELZ,RMAXMM,RMAXME,RMAXFM,RMAXFE,
     3           P0,T0,COMPF,COMTF,AAAF,COMPM,COMTM,AAAM,PRAT)
      IMPLICIT REAL*8 (A-H,O-Z)      
      real*8:: P(MR,NZ),T(MR,NZ),SV(MR,NZ),X(MR,NZ)
      real*8:: P0(MR,NZ),T0(MR,NZ)
      real*8:: AR(MR),V(MR,NZ),DELR(MR)
      real*8:: AZ(MR),DELZ(NZ-1)
      real*8:: C1(2*MR),C2(2*MR),C3(2*MR),C4(2*MR)
      real*8:: C5(2*MR),C6(2*MR),C7(2*MR)
      real*8:: XX(2*MR,NZ),R(2*MR),RR(2*MR)
      integer*4:: IPH(MR,NZ)
      real*8:: PORFR(MR),PERFR(MR),PRAT(MR)
      real*8:: PORMA(MR),PERMA(MR)
      real*8:: BMOLD(MR,NZ),BEOLD(MR,NZ)
      real*8:: ALM1(MR,NZ),ALM2(MR,NZ),ALE1(MR,NZ)
      real*8:: ALE2(MR,NZ),ROM(MR,NZ),ROE(MR,NZ)
	real*8:: XXX(2*MR)
      CHARACTER*6:: NAME
	CHARACTER*5:: NXXX,NAM1,NAM2
	real*8, parameter:: FACP=1.0D-8
      real*8, parameter:: FACT=1.0D-8
	real*8, parameter:: FACS=1.0D-8
	real*8, parameter:: SCALE=1.0D-6
	real*8, parameter:: small=1.0D-5
cc
cc  SOLVE calculates the increments to the solution for one step of the Newton-Raphson method
cc
cc
cc  The labelling system below is used for the matrix blocks and block boundaries
cc
cc         -------------------------------------------------------------------
cc             I,J-1           |           I,J           |            I,J+1
cc                             |                         |
cc              *              |            *            |             *
cc                             |                         |
cc              E              |            F            |             G
cc         -------------------------------------------------------------------
cc                             ZM                        ZP
cc  
cc  pressure            PE,PF,PG
cc  temperature         TE,TF,TG 
cc  2nd variable        XE,XF,XG  (X = T or SV)
cc  phase indicator     IPHE,IPHF,IPHG  (IPH = 0 water, = 1 2-phase, = 2 steam)
cc  saturation          SVE,SVF,SVG
cc  mass accumulation   BMF,BMG
cc  energy accumulation BEF,BEG 
cc  mass flux           FMZM,FMZP
cc  energy flux         FEZM,FEZP
cc
cc  The coefficients in the mass equation for each block are defined by
cc  CME1 = DRM/DPE
cc  CME2 = DRM/DXE
cc  CMF1 = DRM/DPF
cc  CMF2 = DRM/DXF
cc  CMG1 = DRM/DPG
cc  CMG2 = DRM/DXG
cc  For each energy equation the coefficients are
cc  CEE1 = DRE/DPE
cc  CEE2 = DRE/DXE
cc  CEF1 = DRE/DPF
cc  CEF2 = DRE/DXF
cc  CEG1 = DRE/DPG
cc  CEG2 = DRE/DXG
cc  
cc  Note that the equations for the 1st and last blocks are special
cc  Most of these quantities do not have to be stored and therefore 
cc  are not indexed 
cc
cc  Set max residuals to zero
cc
      RMAXMM=0.0D0
      RMAXME=0.0D0
      RMAXFM=0.0D0
      RMAXFE=0.0D0
cc
cc  Loop over all columns
cc
      DO 400 I=1,MR
cc  
cc  Start at the innermost matrix block JE=NZ and set up quantities
cc
         JE=NZ
	   P0E=P0(I,JE)
         T0E=T0(I,JE)
	   PE=P(I,JE)
         TE=T(I,JE)
         XE=X(I,JE)
         SVE=SV(I,JE)
         SLE=1.0-SVE
         IPHE=IPH(I,JE)
         VE=V(I,JE)
	   CRE=CRM
	   RHORE=RHORM
         PORE=PORMA(I)
         PERE=PERMA(I)
         COMP=COMPM
         COMT=COMTM
         AAA=AAAM
cc
         CALL TERMS(PE,TE,SVE,IPHE,PE1,TE1,TE2,
     1             BME,BEE,TME,TEE,BME1,BEE1,TME1,TEE1,
     2             BME2,BEE2,TME2,TEE2,DELPE,DELXE,
     3             PORE,CRE,RHORE,FACP,FACT,FACS,
     4             PERE,P0E,T0E,COMP,COMT,AAA)
cc 
cc  TERMS calculates accumulation terms, and transmissibility terms for PE,TE,SVE and 
cc  incremented values of variables. 
cc  "1" indicates values for PE1=PE+DELPE, TE or SVE are not changed, DELPE is calculated
cc  by TERMS.
cc  "2" indicates values for TE2=TE+DELXE (or SVE2=SVE+DELXE), PE is not changed, DELXE is 
cc  calculated by TERMS.
cc 
	 IF(IGOOD.EQ.2) RETURN
cc  now calculate derivatives
cc  accumulation terms
         BME_PE=(BME1-BME)/DELPE
         BME_XE=(BME2-BME)/DELXE
         BEE_PE=(BEE1-BEE)/DELPE
         BEE_XE=(BEE2-BEE)/DELXE
cc
cc
cc  work down a matrix column 
cc
         DO 401 J=1,NZ-1
cc  shuffle variables along
            JF=NZ+1-J
            JE=JF-1
            JG=JF+1
cc
cc   shuffle terms along
cc
            IF(JF.LT.NZ) THEN
               PG=PF
               TG=TF
               XG=XF
               SVG=SVF
               SLG=SLF
               IPHG=IPHF
               VG=VF
            ENDIF
               PF=PE
               PF1=PE1
               TF1=TE1
               TF=TE
               TF2=TE2
               XF=XE
               SVF=SVE
               SLF=SLE
               IPHF=IPHE
               VF=VE
               BMF=BME
               BEF=BEE
               DELPF=DELPE
               DELXF=DELXE
	       TMF=TME
	       TEF=TEE
	       TMF1=TME1
	       TEF1=TEE1
	       TMF2=TME2
	       TEF2=TEE2         
            IF(JF.LT.NZ) THEN
               AZP=AZM
               DZP=DZM     
	      ENDIF
cc
cc  shuffle along derived quanitities
cc
               BMF_PF=BME_PE 
               BMF_XF=BME_XE     
               BEF_PF=BEE_PE     
               BEF_XF=BEE_XE
            IF(JF.LT.NZ) THEN
               FMZP=FMZM
               FEZP=FEZM
               FMZP_PF=FMZM_PE
               FMZP_XF=FMZM_XE     
               FMZP_PG=FMZM_PF
               FMZP_XG=FMZM_XF     
               FEZP_PF=FEZM_PE
               FEZP_XF=FEZM_XE     
               FEZP_PG=FEZM_PF
               FEZP_XG=FEZM_XF 
            ENDIF    
            PE=P(I,JE)
            TE=T(I,JE)
            P0E=P0(I,JE)
            T0E=T0(I,JE)
            XE=X(I,JE)
            SVE=SV(I,JE)
            SLE=1.0-SVE
            IPHE=IPH(I,JE)
            VE=V(I,JE)
            PORE=PORMA(I)
            PERE=PERMA(I)
            CRE=CRM
            RHORE=RHORM
            AZM=AZ(I)
            DZM=DELZ(JE)
  	    COND=CONDM
            COMP=COMPM
            COMT=COMTM
            AAA=AAAM
cc
cc  permeability for the fracture matrix interface (JE=1) is a special case
cc  and the porosity and other properties for the fracture should be used
cc
            IF(JE.EQ.1) THEN
               PORE=PORFR(I)
               CRE=CRF
               RHORE=RHORF
               COMP=COMPF
               COMT=COMTF
               AAA=AAAF
            ENDIF
            CALL TERMS(PE,TE,SVE,IPHE,PE1,TE1,TE2,
     1            BME,BEE,TME,TEE,BME1,BEE1,TME1,TEE1,
     2            BME2,BEE2,TME2,TEE2,DELPE,DELXE,
     3            PORE,CRE,RHORE,FACP,FACT,FACS,
     4            PERE,P0E,T0E,COMPM,COMTM,AAAM)
cc
	    IF(IGOOD.EQ.2) RETURN
cc  now calculate derivatives
cc  accumulation terms
            BME_PE=(BME1-BME)/DELPE
            BME_XE=(BME2-BME)/DELXE
            BEE_PE=(BEE1-BEE)/DELPE
            BEE_XE=(BEE2-BEE)/DELXE
cc
cc  calculate flux terms using upstream weighting
cc  1st PF > PE
cc
            IF(PF.GT.PE) THEN
cc  basic flows
               FMZM=- TMF*(PF-PE)/DZM
               FEZM=- TEF*(PF-PE)/DZM - COND*(TF-TE)/DZM
cc  flows with incremented pressure in block F
               FMZM1F=- TMF1*(PF1-PE)/DZM
               FEZM1F=- TEF1*(PF1-PE)/DZM - COND*(TF1-TE)/DZM
cc  flows with incremented temperature or saturation in block F
               FMZM2F=- TMF2*(PF-PE)/DZM
               FEZM2F=- TEF2*(PF-PE)/DZM - COND*(TF2-TE)/DZM
cc  flows with incremented pressure in block E
	       FMZM1E=- TMF*(PF-PE1)/DZM
	       FEZM1E=- TEF*(PF-PE1)/DZM - COND*(TF-TE1)/DZM
cc  flows with incremented temperature or saturation in block E
	       FEZM2E=- TEF*(PF-PE)/DZM - COND*(TF-TE2)/DZM
cc  mass flux derivatives
               FMZM_PF=(FMZM1F - FMZM)/DELPF
               FMZM_XF=(FMZM2F - FMZM)/DELXF
               FMZM_PE=(FMZM1E - FMZM)/DELPE
               FMZM_XE=0.0D0
cc  energy flux derivatives
               FEZM_PF=(FEZM1F - FEZM)/DELPF
               FEZM_XF=(FEZM2F - FEZM)/DELXF
               FEZM_PE=(FEZM1E - FEZM)/DELPE
               FEZM_XE=(FEZM2E - FEZM)/DELXE
            ELSEIF(PF.LE.PE) THEN
cc
cc  basic flows
               FMZM=- TME*(PF-PE)/DZM
               FEZM=- TEE*(PF-PE)/DZM - COND*(TF-TE)/DZM
cc  flows with incremented pressure in block F
	       FMZM1F=- TME*(PF1-PE)/DZM
	       FEZM1F=- TEE*(PF1-PE)/DZM - COND*(TF1-TE)/DZM
cc  flows with incremented temperature or saturation in block F
	       FEZM2F=- TEE*(PF-PE)/DZM - COND*(TF2-TE)/DZM
cc  flows with incremented pressure in block E
	       FMZM1E=- TME1*(PF-PE1)/DZM
	       FEZM1E=- TEE1*(PF-PE1)/DZM - COND*(TF-TE1)/DZM
cc  flows with incremented temperature or saturation in block E
	       FMZM2E=- TME2*(PF-PE)/DZM
	       FEZM2E=- TEE2*(PF-PE)/DZM - COND*(TF-TE2)/DZM
cc  mass flux derivatives
               FMZM_PF=(FMZM1F - FMZM)/DELPF
               FMZM_XF=0.0D0
               FMZM_PE=(FMZM1E - FMZM)/DELPE
               FMZM_XE=(FMZM2E - FMZM)/DELXE
cc  energy flux derivatives
               FEZM_PF=(FEZM1F - FEZM)/DELPF
               FEZM_XF=(FEZM2F - FEZM)/DELXF
               FEZM_PE=(FEZM1E - FEZM)/DELPE
               FEZM_XE=(FEZM2E - FEZM)/DELXE
            ENDIF
cc  
cc    set up coefficients 
cc  
cc  
cc  mass terms
cc
            CCF=DT/VF
            CME1=-AZM*FMZM_PE*CCF
            CME2=-AZM*FMZM_XE*CCF
            IF(JF.EQ.NZ) THEN
               RM=-(BMF-BMOLD(I,JF)-FMZM*AZM*CCF)
               CMF1=BMF_PF-AZM*FMZM_PF*CCF
               CMF2=BMF_XF-AZM*FMZM_XF*CCF
               CMG1=0.0D0
               CMG2=0.0D0
            ELSEIF(JF.LT.NZ) THEN
               RM=-(BMF-BMOLD(I,JF)+(FMZP*AZP-FMZM*AZM)*CCF)
               CMF1=BMF_PF+(AZP*FMZP_PF-AZM*FMZM_PF)*CCF
               CMF2=BMF_XF+(AZP*FMZP_XF-AZM*FMZM_XF)*CCF
               CMG1=AZP*FMZP_PG*CCF
               CMG2=AZP*FMZP_XG*CCF
            ENDIF
            RRR=ABS(RM)
            RMAXMM=MAX(RMAXMM,RRR)
cc
cc  energy terms
cc
            RE=-SCALE*(BEF-BEOLD(I,JF)+(FEZP*AZP-FEZM*AZM)*CCF)
            CEE1=-SCALE*AZM*FEZM_PE*CCF
            CEE2=-SCALE*AZM*FEZM_XE*CCF
            IF(JF.EQ.NZ) THEN
               RE=-SCALE*(BEF-BEOLD(I,JF)-FEZM*AZM*CCF)
               CEF1=SCALE*(BEF_PF-AZM*FEZM_PF*CCF)
               CEF2=SCALE*(BEF_XF-AZM*FEZM_XF*CCF)
               CEG1=0.0D0
               CEG2=0.0D0
            ELSEIF(JF.LT.NZ) THEN
               RE=-SCALE*(BEF-BEOLD(I,JF)+(FEZP*AZP-FEZM*AZM)*CCF)
               CEF1=SCALE*(BEF_PF+(AZP*FEZP_PF-AZM*FEZM_PF)*CCF)
               CEF2=SCALE*(BEF_XF+(AZP*FEZP_XF-AZM*FEZM_XF)*CCF)
               CEG1=SCALE*AZP*FEZP_PG*CCF
               CEG2=SCALE*AZP*FEZP_XG*CCF
            ENDIF
            RRR=ABS(RE)
            RMAXME=MAX(RMAXME,RRR)
cc  
cc  now calculate amended coefficients
cc
            IF(JF.EQ.NZ) THEN
               RMS=RM
               CME1S=CME1
               CME2S=CME2
               CMF1S=CMF1
               CMF2S=CMF2
cc
cc  energy terms
cc
               RES=RE
               CEE1S=CEE1
               CEE2S=CEE2
               CEF1S=CEF1
               CEF2S=CEF2
            ELSEIF(JF.LT.NZ) THEN
               RMS=RM-CMG1*ROM(I,JG)-CMG2*ROE(I,JG)
               CME1S=CME1
               CME2S=CME2
               CMF1S=CMF1-CMG1*ALM1(I,JF)-CMG2*ALE1(I,JF)
               CMF2S=CMF2-CMG1*ALM2(I,JF)-CMG2*ALE2(I,JF)
cc
cc  energy terms
cc
               RES=RE-CEG1*ROM(I,JG)-CEG2*ROE(I,JG)
               CEE1S=CEE1
               CEE2S=CEE2
               CEF1S=CEF1-CEG1*ALM1(I,JF)-CEG2*ALE1(I,JF)
               CEF2S=CEF2-CEG1*ALM2(I,JF)-CEG2*ALE2(I,JF)
            ENDIF
cc
cc  now calculate the terms required for later use
cc
            DELTA=CEF2S*CMF1S-CMF2S*CEF1S
            ALM1(I,JE)=(CEF2S*CME1S-CMF2S*CEE1S)/DELTA
            ALM2(I,JE)=(CEF2S*CME2S-CMF2S*CEE2S)/DELTA
            ROM(I,JF)=(CEF2S*RMS-CMF2S*RES)/DELTA
            ALE1(I,JE)=-(CEF1S*CME1S-CMF1S*CEE1S)/DELTA
            ALE2(I,JE)=-(CEF1S*CME2S-CMF1S*CEE2S)/DELTA
            ROE(I,JF)=-(CEF1S*RMS-CMF1S*RES)/DELTA
401      CONTINUE
cc
cc  fracture equation for block (I,1)
cc  The labelling syatem below is used for the blocks and block boundaries
cc
cc          -----------------------------------------------------------------
cc                             |                         |
cc                             |                         |
cc                             |          I,2            |
cc                             |           *             |
cc                             |           D             |
cc                             |                         |
cc                             |                         |
cc          ------------------------------ Q ---------------------------------- 
cc                             |                         |
cc             I-1,1           |            I,1          |            I+1,1
cc                             |                         |
cc              *              |            *            |             *
cc                             |                         |
cc              A              |            B            |             C
cc                             |                         |
cc           ------------------------------------------------------------------
cc                             M                         P
cc  
cc
cc
cc  first rename final matrix variables using fracture notation
cc
            IB=I-1
            PRATC=PRAT(I)
            PC=PE
            TC=TE
            XC=XE
            SVC=SVE
            SLC=1.0-SVC
            IPHC=IPHE
            VC=VE
            PORC=PORE
            BMC=BME
            BEC=BEE
            BMC_PC=BME_PE
            BEC_PC=BEE_PE
            BMC_XC=BME_XE
            BEC_XC=BEE_XE
            TMC=TME*PRATC
            TEC=TEE*PRATC
            TMC1=TME1*PRATC
            TEC1=TEE1*PRATC
            PC1=PE1
            TC1=TE1
            TMC2=TME2*PRATC
            TEC2=TEE2*PRATC
            TC2=TE2 
	      DELPC=DELPE
	      DELXC=DELXE
cc  
       IF(IB.GT.0) THEN
            AP=AR(IB)
            DRP=DELR(IB)
       ENDIF
cc
		  COND=CONDF
cc
cc  only calculate fluxes and set up fracture equations for I>1
cc
		  IF(I.GE.2) THEN
cc
cc  calculate flux terms using upstream weighting
cc  1st PC > PB
cc
               IF(PC.GT.PB) THEN
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
cc
cc  or PC < PB or = PB
cc
      ELSEIF(PC.LE.PB) THEN
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
      ENDIF
cc
cc  set up mass coefficients and mass residual for (I-1)th equation
cc
      CCB=DT/VB
      IF(IB.EQ.1) THEN
         RMB=-(BMB-BMOLD(IB,1)+(FMP*AP+FMQ*AQ-QMM)*CCB)
         CMA1=0.0D0
         CMA2=0.0D0
         CMB1=BMB_PB+(FMP_PB*AP+FMQ_PB*AQ)*CCB
         CMB2=BMB_XB+(FMP_XB*AP+FMQ_XB*AQ)*CCB
      ELSEIF(IB.GT.1) THEN
         RMB=-(BMB-BMOLD(IB,1)+(FMP*AP+FMQ*AQ-FMM*AM)*CCB)
         CMA1=-AM*FMM_PA*CCB
         CMA2=-AM*FMM_XA*CCB
         CMB1=BMB_PB+(FMP_PB*AP+FMQ_PB*AQ-FMM_PB*AM)*CCB
         CMB2=BMB_XB+(FMP_XB*AP+FMQ_XB*AQ-FMM_XB*AM)*CCB
      ENDIF
      CMC1=FMP_PC*AP*CCB
      CMC2=FMP_XC*AP*CCB
      CMD1=FMQ_PD*AQ*CCB
      CMD2=FMQ_XD*AQ*CCB
      L=2*IB-1
      R(L)=RMB-CMD1*ROM(IB,2)-CMD2*ROE(IB,2)
      RR(L)=R(L)
      RRR=ABS(RMB)
      RMAXFM=MAX(RMAXFM,RRR)
      C1(L)=0.0D0
      C2(L)=CMA1
      C3(L)=CMA2
      C4(L)=CMB1-CMD1*ALM1(IB,1)-CMD2*ALE1(IB,1)
      C5(L)=CMB2-CMD1*ALM2(IB,1)-CMD2*ALE2(IB,1)
      C6(L)=CMC1
      C7(L)=CMC2
cc
cc  energy equation coefficients and energy residual
      IF(IB.EQ.1) THEN
cc  production/injection term
      IF(QMM.GE.0.0) THEN
         QEM=QMM*HIN
         QEM_PB=0.0D0
         QEM_XB=0.0D0
      ELSEIF(QMM.LT.0.0) THEN
         QEM=QMM*TEB/TMB
         QEM_PB=(QMM*TEB1/TMB1 - QEM)/DELPB
         QEM_XB=(QMM*TEB2/TMB2 - QEM)/DELXB
      ENDIF
cc
c     Flowing enthalpy calculation (AC 13/10/99):
      if (dabs(QMM)>small) then
	  HF=QEM/QMM
	else
	  HF=0.0D0
	end if

         REB=-SCALE*(BEB-BEOLD(IB,1)+(FEP*AP+FEQ*AQ-QEM)*CCB)
         CEA1=0.0D0
         CEA2=0.0D0
         CEB1=SCALE*(BEB_PB+(FEP_PB*AP+FEQ_PB*AQ-QEM_PB)*CCB)
         CEB2=SCALE*(BEB_XB+(FEP_XB*AP+FEQ_XB*AQ-QEM_XB)*CCB)
      ELSEIF(IB.GT.1) THEN
         REB=-SCALE*(BEB-BEOLD(IB,1)+(FEP*AP+FEQ*AQ-FEM*AM)*CCB)
         CEA1=-SCALE*AM*FEM_PA*CCB
         CEA2=-SCALE*AM*FEM_XA*CCB
         CEB1=SCALE*(BEB_PB+(FEP_PB*AP+FEQ_PB*AQ-FEM_PB*AM)*CCB)
         CEB2=SCALE*(BEB_XB+(FEP_XB*AP+FEQ_XB*AQ-FEM_XB*AM)*CCB)
      ENDIF
      CEC1=SCALE*(FEP_PC*AP*CCB)
      CEC2=SCALE*(FEP_XC*AP*CCB)
      CED1=SCALE*(FEQ_PD*AQ*CCB)
      CED2=SCALE*(FEQ_XD*AQ*CCB)
      L1=2*IB
      R(L1)=REB-CED1*ROM(IB,2)-CED2*ROE(IB,2)
      RR(L1)=R(L1)
      RRR=ABS(REB)
      RMAXFE=MAX(RMAXFE,RRR)
      C1(L1)=CEA1
      C2(L1)=CEA2
      C3(L1)=CEB1-CED1*ALM1(IB,1)-CED2*ALE1(IB,1)
      C4(L1)=CEB2-CED1*ALM2(IB,1)-CED2*ALE2(IB,1)
      C5(L1)=CEC1
      C6(L1)=CEC2
      C7(L1)=0.0D0
      ENDIF
cc
cc  now shuffle terms along ready for next calculation at (I+1)th block
cc
      IF(I.GT.1) THEN
      PA=PB
      TA=TB
      XA=XB
      SVA=SVB
      SLA=SLB
      IPHA=IPHB
      VA=VB
      AM=AP
      DRM=DRP
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
	ENDIF
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
      DELPB=DELPC
      DELXB=DELXC
	TMB=TMC
	TEB=TEC
	TMB1=TMC1
	TEB1=TEC1
	TMB2=TMC2
	TEB2=TEC2
cc  shuffle along derived quanitities
      BMB_PB=BMC_PC 
      BMB_XB=BMC_XC     
      BEB_PB=BEC_PC     
      BEB_XB=BEC_XC
 
cc  
cc  calculate vertical flow and other terms required for I+1th equation
cc
            AQ=AZM
            FMQ=FMZM
            FEQ=FEZM
            FMQ_PD=FMZM_PF
            FMQ_PB=FMZM_PE
            FMQ_XD=FMZM_XF
            FMQ_XB=FMZM_XE
            FEQ_PD=FEZM_PF
            FEQ_PB=FEZM_PE
            FEQ_XD=FEZM_XF
            FEQ_XB=FEZM_XE
      
cc
cc  now complete the equation for the end block if I=MR
cc  
         IF(I.EQ.MR) THEN
            IB=MR
            CCB=DT/VB
            RMB=-(BMB-BMOLD(IB,1)+(FMQ*AQ-FMM*AM)*CCB)
            CMA1=-AM*FMM_PA*CCB
            CMA2=-AM*FMM_XA*CCB
            CMB1=BMB_PB+(FMQ_PB*AQ-FMM_PB*AM)*CCB
            CMB2=BMB_XB+(FMQ_XB*AQ-FMM_XB*AM)*CCB
            CMC1=0.0D0
            CMC2=0.0D0
            CMD1=FMQ_PD*AQ*CCB
            CMD2=FMQ_XD*AQ*CCB
            L=2*IB-1
            R(L)=RMB-CMD1*ROM(IB,2)-CMD2*ROE(IB,2)
            RR(L)=R(L)
            RRR=ABS(RMB)
            RMAXFM=MAX(RMAXFM,RRR)
            C1(L)=0.0D0
            C2(L)=CMA1
            C3(L)=CMA2
            C4(L)=CMB1-CMD1*ALM1(IB,1)-CMD2*ALE1(IB,1)
            C5(L)=CMB2-CMD1*ALM2(IB,1)-CMD2*ALE2(IB,1)
            C6(L)=CMC1
            C7(L)=CMC2
cc
cc  energy equation coefficients and energy residual
cc
            REB=-SCALE*(BEB-BEOLD(IB,1)+(FEQ*AQ-FEM*AM)*CCB)
            CEA1=-SCALE*AM*FEM_PA*CCB
            CEA2=-SCALE*AM*FEM_XA*CCB
            CEB1=SCALE*(BEB_PB+(FEQ_PB*AQ-FEM_PB*AM)*CCB)
            CEB2=SCALE*(BEB_XB+(FEQ_XB*AQ-FEM_XB*AM)*CCB)
            CEC1=0.0D0
            CEC2=0.0D0
            CED1=SCALE*(FEQ_PD*AQ*CCB)
            CED2=SCALE*(FEQ_XD*AQ*CCB)
            L1=2*IB
            R(L1)=REB-CED1*ROM(IB,2)-CED2*ROE(IB,2)
            RR(L1)=R(L1)
            RRR=ABS(REB)
            RMAXFE=MAX(RMAXFE,RRR)
            C1(L1)=CEA1
            C2(L1)=CEA2
            C3(L1)=CEB1-CED1*ALM1(IB,1)-CED2*ALE1(IB,1)
            C4(L1)=CEB2-CED1*ALM2(IB,1)-CED2*ALE2(IB,1)
            C5(L1)=CEC1
            C6(L1)=CEC2
            C7(L1)=0.0D0
         ENDIF
400   CONTINUE
cc
cc  now solve the amended fracture equations
cc
      M2=2*MR
      CALL SEVEN1(C1,C2,C3,C4,C5,C6,C7,R,XXX,M2)
      DO 402 I=1,M2
         XX(I,1)=XXX(I)
402   CONTINUE
406   FORMAT(1H ,I5,2X,E12.5)
cc
cc  now back substitute to calculate the matrix unknowns
cc
      DO 403 I=1,MR
         DO 404 J=2,NZ
            L=2*I-1
            L1=2*I
            JA=J-1
            XX(L,J)=ROM(I,J)-ALM1(I,JA)*XX(L,JA)-ALM2(I,JA)*XX(L1,JA)
            XX(L1,J)=ROE(I,J)-ALE1(I,JA)*XX(L,JA)-ALE2(I,JA)*XX(L1,JA)
404      CONTINUE
403   CONTINUE
      RETURN
      end subroutine SOLVE
cc       
cc----7--1--------2---------3---------4---------5---------6---------7---------8 
cc        
      SUBROUTINE UPDATE(MR,NZ,IPH,P,T,SV,X,XX,EPS)
      IMPLICIT REAL*8 (A-H,O-Z)       
      DIMENSION P(MR,NZ),T(MR,NZ),X(MR,NZ),SV(MR,NZ),XX(2*MR,NZ)
      DIMENSION IPH(MR,NZ)
      COMMON/WRIT/IWRIT,IT,IW
cc
cc  This subroutine adds the increments on to the primary variables and carries out
cc  phase changes
cc
      do J=1,NZ
         do I=1,MR
            L=2*I-1
            L1=L+1
            PP=P(I,J)
            TP=T(I,J)
            IPHP=IPH(I,J)
            PT=PP+XX(L,J)
            XT=X(I,J)+XX(L1,J)
cc
            select case (IPHP)
cc  -----------------------------------------------------------------------------
            case (0)
cc
cc  currently single phase liquid water
               TT=XT
               CALL SAT(TT,PS)
               if (PT>PS) then
cc  remain single phase liquid water
                  P(I,J)=PT
                  T(I,J)=TT
                  X(I,J)=TT
                  SV(I,J)=0.0
                  IPH(I,J)=0
               else
cc  change to 2-phase
                  P(I,J)=PS
                  T(I,J)=TT
                  SV(I,J)=EPS
                  X(I,J)=EPS
                  IPH(I,J)=1
               end if
cc
cc  -----------------------------------------------------------------------------
            case (1)
cc
cc  currently 2-phase
               SVT=XT
               CALL SAT(TP,PS)
               IF((SVT.GT.0.0).AND.(SVT.LT.1.0)) THEN
cc  remain 2-phase
                  P(I,J)=PT
                  SV(I,J)=SVT
                  X(I,J)=SVT
                  IPH(I,J)=1
	          T0=0.0
                  CALL TSAT(PT,T0,TT)
	            if (IGOOD>0) return 
                  T(I,J)=TT
               ELSEIF(SVT.LE.0.0) THEN
cc  change to single phase liquid water
                  P(I,J)=1.000001*PS
                  T(I,J)=TP
                  X(I,J)=TP
                  SV(I,J)=0.0
                  IPH(I,J)=0
               ELSEIF(SVT.GE.1.0) THEN
cc  change to single phase dry steam
                  P(I,J)=0.999999*PS
                  T(I,J)=TP
                  X(I,J)=TP
                  SV(I,J)=1.0
                  IPH(I,J)=2
               ENDIF
cc
cc  -----------------------------------------------------------------------------
            case (2)
cc
cc  currently single phase dry steam
               TT=XT
               CALL SAT(TT,PS)
               if (PT<PS) then
cc  remain single phase dry steam
                  P(I,J)=PT
                  T(I,J)=TT
                  X(I,J)=TT
                  SV(I,J)=1.0
                  IPH(I,J)=2
               else
cc  change to 2-phase
                  P(I,J)=PS
                  T(I,J)=TT
                  SV(I,J)=1.0-EPS
                  X(I,J)=1.0-EPS
                  IPH(I,J)=1
               end if
	     end select
        end do
      end do
      RETURN
      end subroutine UPDATE
cc       
cc----7--1--------2---------3---------4---------5---------6---------7---------8 
cc
      SUBROUTINE INIT1(MR,NZ,POLD,XOLD,TOLD,SVOLD,IPHOLD)
      IMPLICIT REAL*8 (A-H,O-Z)       
      DIMENSION POLD(MR,NZ),XOLD(MR,NZ),TOLD(MR,NZ)
	DIMENSION SVOLD(MR,NZ),IPHOLD(MR,NZ)
      DIMENSION P(MR,NZ),T(MR,NZ),SV(MR,NZ)
      do J=1,NZ
         do I=1,MR
            PP=POLD(I,J)
            XP=XOLD(I,J)
	      CALL TSAT(PP,TZERO,TS)
            IF(XP.LT.1.0) THEN
cc  2-phase conditions
               TZERO=0.0
               TP=TS 
               TOLD(I,J)=TS
               SVOLD(I,J)=XP
               SVP=XP
               IPHOLD(I,J)=1
               IPHP=1
            ELSEIF(XP.GT.1.5) THEN
               if (XP>=TS) then
cc   single phase dry steam
                  TP=XP
                  TOLD(I,J)=XP
                  SVOLD(I,J)=1.0
                  SVP=1.0
                  IPHOLD(I,J)=2
                  IPHP=2
               else
cc   single phase liquid water
                  TP=XP
                  TOLD(I,J)=XP
                  SVOLD(I,J)=0.0
                  SVP=0.0
                  IPHOLD(I,J)=0
                  IPHP=0
               end if
             ENDIF
          end do
       end do
       RETURN
       end subroutine INIT1
cc       
cc----7--1--------2---------3---------4---------5---------6---------7---------8 
cc
      SUBROUTINE INIT2(MR,NZ,POLD,XOLD,TOLD,SVOLD,IPHOLD,
     1                  BMOLD,BEOLD,PORFR,RHORF,CRF,PORMA,RHORM,CRM,
     2                  P0,T0,COMPF,COMTF,AAAF,COMPM,COMTM,AAAM)
      IMPLICIT REAL*8 (A-H,O-Z)       
      DIMENSION POLD(MR,NZ),XOLD(MR,NZ),TOLD(MR,NZ)
      DIMENSION P0(MR,NZ),T0(MR,NZ)
	DIMENSION SVOLD(MR,NZ),IPHOLD(MR,NZ)
      DIMENSION BMOLD(MR,NZ),BEOLD(MR,NZ)
      DIMENSION PORFR(MR),PORMA(MR)
      do J=1,NZ
         do I=1,MR
            PP=POLD(I,J)
            TP=TOLD(I,J)
            P0P=P0(I,J)
            T0P=T0(I,J)
	      SVP=SVOLD(I,J)
	      if (J==1) then
               PORP=PORFR(I)
	         RHOR=RHORF
	         CR=CRF
               COMP=COMPF
               COMT=COMTF
               AAA=AAAF
            else
               PORP=PORMA(I)
	         RHOR=RHORM
	         CR=CRM
               COMP=COMPM
               COMT=COMTM
               AAA=AAAM
            end if
            IPHP=IPHOLD(I,J)
            CALL THERMO(PP,TP,SVP,IPHP,RHOL,RHOV,HL,HV,VISL,VISV)
            CALL ACCUM(BM,BE,IPHP,PP,TP,SVP,RHOL,RHOV,HL,HV,PORP,
     1                 RHOR,CR,P0P,T0P,COMP,COMT,AAA,PERFAC)
            BMOLD(I,J)=BM
            BEOLD(I,J)=BE
          end do
       end do
       RETURN
       end subroutine INIT2
cc       
cc----7--1--------2---------3---------4---------5---------6---------7---------8 

      end module NumericalSimulatorMultiLayer_Routines