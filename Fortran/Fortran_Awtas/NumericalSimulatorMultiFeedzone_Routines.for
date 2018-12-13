      module NumericalSimulatorMultiFeedzone_Routines

	use MultiFeedzoneThermodynamics
	use MatrixSolvers
	use variable_types

	contains

cc----7--1--------2---------3---------4---------5---------6---------7---------8 

      SUBROUTINE SOLVEMF(P,T,SV,X,IPH,XX,BMOLD,BEOLD,NR,MZ,PORWB,PORMA,
     1      PERWB,PERMA,AR,AZ,V,CRWB,CRMA,RHORWB,RHORMA,CONDWB,CONDMA,
     2      QMM,HIN,DT,DELR,DELZ,RMAXMM,RMAXME,RMAXFM,RMAXFE,
     3      P0,T0,COMPWB,COMTWB,AAAWB,COMPMA,COMTMA,AAAMA,PRAT,IACT)

      implicit none
	integer(I4B):: MZ,NR
      real(DP):: P(MZ,NR),T(MZ,NR),SV(MZ,NR),X(MZ,NR)
      integer(I4B):: IPH(MZ,NR)
 	real(DP):: CRWB,RHORWB,CONDWB,QMM,HIN,DT
	real(DP):: RMAXM,RMAXE,RMAXFM,RMAXFE
	real(DP):: RMAXMM,RMAXME,P0E,T0E,PE,TE,XE,SVE,SLE,VE,CRE,RHORE
	real(DP):: PORE,PERE,COMP,COMT,AAA,PE1,TE1,TE2,BME,BEE,TMLE,TMVE
	real(DP):: TELE,TEVE,BME1,BEE1,TMLE1,TMVE1,TELE1,TEVE1,BME2,BEE2
	real(DP):: TMLE2,TMVE2,TELE2,TEVE2,DELPE,DELXE,RHOLE,RHOVE,RHOLE1
	real(DP):: RHOVE1,RHOLE2,RHOVE2,PERFAC,TME,TEE,TME1,TEE1,TME2,TEE2
	real(DP):: BME_PE, BME_XE,BEE_PE,BEE_XE,PG,PF,TG,TF,XG,XF,SVG,SVF
	real(DP):: SLG,SLF,VG,VF,PF1,TF1,TF2,BMF,BEF,DELPF,DELXF,TMF,TEF
	real(DP):: TMF1,TEF1,TMF2,TEF2,ARP,ARM,DRP,DRM,BMF_PF,BMF_XF,BEF_PF
	real(DP):: BEF_XF,FMRP,FMRM,FERP,FERM,FMRP_PF,FMRP_PE,FMRP_XF
	real(DP):: FMRP_XE,FMRP_PG,FMRP_XG,FERP_PF,FERP_FERM_PE
	real(DP):: FERP_XF,FERM_XE,FERP_PG,FMRM_PE,FMRM_XE,FMRM_PF,FMRM_XF
	real(DP):: FERM_PE,FERM_PF,FERP_XG,FERM_XF,COND
	real(DP):: FMRM1F,FERM1F,FMRM2F,FERM2F,FMRM1E,FERM1E,FERM2E,FMRM2E
	real(DP):: CCF,CME1,CME2,RM,CMF1,CMF2,CMG1,CMG2,RRR,RE,CEE1,CEE2
	real(DP):: CEF1,CEF2,CEG1,CEG2,RMS,CME1S,CME2S,CMF1S,CMF2S,RES
	real(DP):: CEE1S,CEE2S,CEF1S,CEF2S,DELTA,PRATC,PC,TC,XC,SVC,SLC,VC
	real(DP):: PORC,BMC,BEC,BMC_PC,BEC_PC,BMC_XC,BEC_XC,TMLC,TMVC,TELC
	real(DP):: TEVC,TMLC1,TMVC1,TELC1,TEVC1,TMLC2,TMVC2,TELC2,TEVC2
	real(DP):: RHOLC,RHOVC,RHOLC1,RHOVC1,RHOLC2,RHOVC2,PC1,TC1,TC2,DELPC
	real(DP):: DELXC,AP,DZP,GRADL,PB,RHOLB,GRADV,RHOVB,GRADL1C,GRADV1C
	real(DP):: GRADL2C,GRADV2C,GRADL1B,PB1,RHOLB1,GRADV1B,RHOVB1,GRADL2B
	real(DP):: RHOLB2,GRADV2B,RHOVB2,FMPL,FEPL,FMPL1C,FEPL1C,FMPL2C
	real(DP):: FEPL2C,FMPL1B,FMPL2B,FEPL2B,TMLB,TELB,TMLB1,TELB1,TMLB2
	real(DP):: TELB2,FMPV,FEPV,FMPV1C,FEPV1C,FEPL1B,FMPV2C,FEPV2C,FMPV1B
	real(DP):: FEPV1B,FMPV2B,TMVB,TEVB,TMVB1,TEVB1,TMVB2,TEVB2,FMP,FEP
	real(DP):: TB,FMP1C,FEP1C,FMP2C,FEP2C,FMP1B,FEP1B,TB1,FMP2B,TB2
	real(DP):: FMP_PC,FMP_XC,FMP_PB,DELPB,FEPV2B,FEP2B,FMP_XB,DELXB
	real(DP):: FEP_PC,FEP_XC,FEP_PB,FEP_XB,FMQ,AQ,CCB,VB,RMB,BMB,CMB1
	real(DP):: BMB_PB,FMQ_PB,CMB2,BMB_XB,FMQ_XB,CMA1,CMA2,DELBM,DELFM
	real(DP):: FMM,AM,FMM_PB,FMM_XB,FMM_PA,FMM_XA,CMC1,CMC2,CMD1,FMQ_PD
	real(DP):: CMD2,FMQ_XD,QEM,QEM_PB,QEM_XB,QEM1,QEM2,REB,BEB,FEQ,CEB1
	real(DP):: BEB_PB,FEQ_PB,CEB2,BEB_XB,FEQ_XB,CEA1,CEA2,FEM,FEM_XB
	real(DP):: FEM_PA,FEM_XA,CEC1,CEC2,FEM_PB,CED1,FEQ_PD,CED2,FEQ_XD
	real(DP):: PA,TA,XA,XB,SVA,SVB,SLA,SLB,VA,DZM,Z1,Z2,Z3,Z4,Z5,Z6
	real(DP):: COMPWB,COMTWB,AAAWB
      real(DP):: P0(MZ,NR),T0(MZ,NR)
      real(DP):: AR(MZ,NR),V(MZ,NR),DELR(NR)
      real(DP):: AZ(NR),DELZ(MZ-1)
      real(DP):: C1(2*NR),C2(2*NR),C3(2*NR),C4(2*NR)
      real(DP):: C5(2*NR),C6(2*NR),C7(2*NR)
      real(DP):: XX(2*MZ,NR),R(2*NR),RR(2*NR)
      real(DP):: PORWB(MZ),PERWB(MZ),PRAT(MZ)
      real(DP):: PORMA(MZ,NR),PERMA(MZ,NR)
      real(DP):: BMOLD(MZ,NR),BEOLD(MZ,NR)
      real(DP):: ALM1(MZ,NR),ALM2(MZ,NR),ALE1(MZ,NR)
      real(DP):: ALE2(MZ,NR),ROM(MZ,NR),ROE(MZ,NR)
      real(DP):: XXX(2*NR)
      real(DP):: FMR(MZ,NR),FER(MZ,NR)
      real(DP):: FMZ(MZ),FEZ(MZ),FMX(MZ),CONDMA(MZ)
      real(DP):: RHORMA(MZ),COMPMA(MZ),COMTMA(MZ),AAAMA(MZ),CRMA(MZ)
      integer(I4B):: IACT(MZ)
	integer(I4B):: I,J,JE,JF,JG,IPHE,IPHG,IPHF,IB,IPHC,L,L1,IPHA,IPHB,M2
	integer(I4B):: JA
      CHARACTER*6 NAME
      CHARACTER*5 NXXX,NAM1,NAM2
	real(DP), parameter:: FACP=1.0D-8, FACT=1.0D-8, FACS=1.0D-8
	real(DP), parameter:: SCALE=1.0D-6, G=9.81D0
cc
cc  SOLVE calculates the increments to the solution for one step of the Newton-Raphson method
cc
cc
cc  The labelling system below is used for the matrix blocks 
cc  and block boundaries in each layer
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
cc  mass flux           FMRM,FMRP
cc  energy flux         FERM,FERP
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
      RMAXMM=0.0
      RMAXME=0.0
      RMAXFM=0.0
      RMAXFE=0.0
cc
cc  Loop over all layers
cc
      do I=1,MZ
cc  
cc  Start at the innermost matrix block JE=NR and set up quantities
cc  Include only "active layers"
cc
         IF(IACT(I)==1) THEN
            JE=NR
	      P0E=P0(I,JE)
            T0E=T0(I,JE)
	      PE=P(I,JE)
            TE=T(I,JE)
            XE=X(I,JE)
            SVE=SV(I,JE)
            SLE=1.0-SVE
            IPHE=IPH(I,JE)
            VE=V(I,JE)
	      CRE=CRMA(I)
	      RHORE=RHORMA(I)
            PORE=PORMA(I,JE)
            PERE=PERMA(I,JE)
            COMP=COMPMA(I)
            COMT=COMTMA(I)
            AAA=AAAMA(I)
cc
         CALL TERMSMF(PE,TE,SVE,IPHE,PE1,TE1,TE2,BME,BEE,
     1        TMLE,TMVE,TELE,TEVE,BME1,BEE1,TMLE1,TMVE1,TELE1,TEVE1,
     2        BME2,BEE2,TMLE2,TMVE2,TELE2,TEVE2,DELPE,DELXE,
     3        RHOLE,RHOVE,RHOLE1,RHOVE1,RHOLE2,RHOVE2,
     4        PORE,CRE,RHORE,FACP,FACT,FACS,PERE,P0E,T0E,COMP,COMT,AAA,
     5        PERFAC)
	      TME=TMLE+TMVE
            TEE=TELE+TEVE
            TME1=TMLE1+TMVE1
            TEE1=TELE1+TEVE1
            TME2=TMLE2+TMVE2
            TEE2=TELE2+TEVE2
cc 
cc  TERMS calculates accumulation terms, and transmissibility terms for PE,TE,SVE and 
cc  incremented values of variables. 
cc  "1" indicates values for PE1=PE+DELPE, TE or SVE are not changed, DELPE is calculated
cc  by TERMS.
cc  "2" indicates values for TE2=TE+DELXE (or SVE2=SVE+DELXE), PE is not changed, DELXE is 
cc  calculated by TERMS.
cc 
	 IF(IGOOD==2) RETURN
cc  now calculate derivatives
cc  accumulation terms
         BME_PE=(BME1-BME)/DELPE
         BME_XE=(BME2-BME)/DELXE
         BEE_PE=(BEE1-BEE)/DELPE
         BEE_XE=(BEE2-BEE)/DELXE
cc
cc
cc  work back along a layer 
cc
         do J=1,NR-1
cc  shuffle variables along
            JF=NR+1-J
            JE=JF-1
            JG=JF+1
cc
cc   shuffle terms along
cc
            IF(JF<NR) THEN
               PG=PF
               TG=TF
               XG=XF
               SVG=SVF
               SLG=SLF
               IPHG=IPHF
               VG=VF
            END IF
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
            IF(JF<NR) THEN
              ARP=ARM
              DRP=DRM     
	      END IF
cc
cc  shuffle along derived quanitities
cc
            BMF_PF=BME_PE 
            BMF_XF=BME_XE     
            BEF_PF=BEE_PE     
            BEF_XF=BEE_XE
            IF(JF<NR) THEN
               FMRP=FMRM
               FERP=FERM
               FMRP_PF=FMRM_PE
               FMRP_XF=FMRM_XE     
               FMRP_PG=FMRM_PF
               FMRP_XG=FMRM_XF     
               FERP_PF=FERM_PE
               FERP_XF=FERM_XE     
               FERP_PG=FERM_PF
               FERP_XG=FERM_XF 
            END IF    
            PE=P(I,JE)
            TE=T(I,JE)
            P0E=P0(I,JE)
            T0E=T0(I,JE)
            XE=X(I,JE)
            SVE=SV(I,JE)
            SLE=1.0-SVE
            IPHE=IPH(I,JE)
            VE=V(I,JE)
            PORE=PORMA(I,JE)
            PERE=PERMA(I,JE)
            CRE=CRMA(I)
            RHORE=RHORMA(I)
            ARM=AR(I,JE)
            DRM=DELR(JE)
  	      COND=CONDMA(I)
            COMP=COMPMA(I)
            COMT=COMTMA(I)
            AAA=AAAMA(I)
cc
cc  permeability for the wellbore/matrix interface (JE=1) is a special case
cc  and the porosity and other properties for the wellbore should be used
cc
            IF(JE==1) THEN
               PORE=PORWB(I)
               CRE=CRWB
               RHORE=RHORWB
               COMP=COMPWB
               COMT=COMTWB
               AAA=AAAWB
            END IF
            CALL TERMSMF(PE,TE,SVE,IPHE,PE1,TE1,TE2,BME,BEE,
     1         TMLE,TMVE,TELE,TEVE,BME1,BEE1,TMLE1,TMVE1,TELE1,TEVE1,
     2         BME2,BEE2,TMLE2,TMVE2,TELE2,TEVE2,DELPE,DELXE,
     3         RHOLE,RHOVE,RHOLE1,RHOVE1,RHOLE2,RHOVE2,
     4         PORE,CRE,RHORE,FACP,FACT,FACS,PERE,P0E,T0E,COMP,COMT,AAA,
     5         PERFAC)
            TME=TMLE+TMVE
            TEE=TELE+TEVE
            TME1=TMLE1+TMVE1
            TEE1=TELE1+TEVE1
            TME2=TMLE2+TMVE2
            TEE2=TELE2+TEVE2
cc
	    IF (IGOOD==2) RETURN
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
            IF(PF>PE) THEN
cc  basic flows
               FMRM=- TMF*(PF-PE)/DRM
               FERM=- TEF*(PF-PE)/DRM - COND*(TF-TE)/DRM
cc  flows with incremented pressure in block F
               FMRM1F=- TMF1*(PF1-PE)/DRM
               FERM1F=- TEF1*(PF1-PE)/DRM - COND*(TF1-TE)/DRM
cc  flows with incremented temperature or saturation in block F
               FMRM2F=- TMF2*(PF-PE)/DRM
               FERM2F=- TEF2*(PF-PE)/DRM - COND*(TF2-TE)/DRM
cc  flows with incremented pressure in block E
	         FMRM1E=- TMF*(PF-PE1)/DRM
	         FERM1E=- TEF*(PF-PE1)/DRM - COND*(TF-TE1)/DRM
cc  flows with incremented temperature or saturation in block E
	         FERM2E=- TEF*(PF-PE)/DRM - COND*(TF-TE2)/DRM
cc  mass flux derivatives
               FMRM_PF=(FMRM1F - FMRM)/DELPF
               FMRM_XF=(FMRM2F - FMRM)/DELXF
               FMRM_PE=(FMRM1E - FMRM)/DELPE
               FMRM_XE=0.0
cc  energy flux derivatives
               FERM_PF=(FERM1F - FERM)/DELPF
               FERM_XF=(FERM2F - FERM)/DELXF
               FERM_PE=(FERM1E - FERM)/DELPE
               FERM_XE=(FERM2E - FERM)/DELXE
            ELSE
cc
cc  basic flows
               FMRM=- TME*(PF-PE)/DRM
               FERM=- TEE*(PF-PE)/DRM - COND*(TF-TE)/DRM
cc  flows with incremented pressure in block F
	       FMRM1F=- TME*(PF1-PE)/DRM
	       FERM1F=- TEE*(PF1-PE)/DRM - COND*(TF1-TE)/DRM
cc  flows with incremented temperature or saturation in block F
	       FERM2F=- TEE*(PF-PE)/DRM - COND*(TF2-TE)/DRM
cc  flows with incremented pressure in block E
	       FMRM1E=- TME1*(PF-PE1)/DRM
	       FERM1E=- TEE1*(PF-PE1)/DRM - COND*(TF-TE1)/DRM
cc  flows with incremented temperature or saturation in block E
	       FMRM2E=- TME2*(PF-PE)/DRM
	       FERM2E=- TEE2*(PF-PE)/DRM - COND*(TF-TE2)/DRM
cc  mass flux derivatives
               FMRM_PF=(FMRM1F - FMRM)/DELPF
               FMRM_XF=0.0
               FMRM_PE=(FMRM1E - FMRM)/DELPE
               FMRM_XE=(FMRM2E - FMRM)/DELXE
cc  energy flux derivatives
               FERM_PF=(FERM1F - FERM)/DELPF
               FERM_XF=(FERM2F - FERM)/DELXF
               FERM_PE=(FERM1E - FERM)/DELPE
               FERM_XE=(FERM2E - FERM)/DELXE
            END IF
            FMR(I,JE)=FMRM*ARM
            FER(I,JE)=FERM*ARM

cc  set up coefficients 
   
cc  mass terms

            CCF=DT/VF
            CME1=-ARM*FMRM_PE*CCF
            CME2=-ARM*FMRM_XE*CCF
            IF(JF==NR) THEN
               RM=-(BMF-BMOLD(I,JF)-FMRM*ARM*CCF)
               CMF1=BMF_PF-ARM*FMRM_PF*CCF
               CMF2=BMF_XF-ARM*FMRM_XF*CCF
               CMG1=0.0
               CMG2=0.0
            ELSE IF(JF<NR) THEN
               RM=-(BMF-BMOLD(I,JF)+(FMRP*ARP-FMRM*ARM)*CCF)
               CMF1=BMF_PF+(ARP*FMRP_PF-ARM*FMRM_PF)*CCF
               CMF2=BMF_XF+(ARP*FMRP_XF-ARM*FMRM_XF)*CCF
               CMG1=ARP*FMRP_PG*CCF
               CMG2=ARP*FMRP_XG*CCF
            END IF
            RRR=ABS(RM)
            RMAXMM=MAX(RMAXMM,RRR)
cc
cc  energy terms
cc
            RE=-SCALE*(BEF-BEOLD(I,JF)+(FERP*ARP-FERM*ARM)*CCF)
            CEE1=-SCALE*ARM*FERM_PE*CCF
            CEE2=-SCALE*ARM*FERM_XE*CCF
            IF(JF==NR) THEN
               RE=-SCALE*(BEF-BEOLD(I,JF)-FERM*ARM*CCF)
               CEF1=SCALE*(BEF_PF-ARM*FERM_PF*CCF)
               CEF2=SCALE*(BEF_XF-ARM*FERM_XF*CCF)
               CEG1=0.0
               CEG2=0.0
            ELSE IF(JF<NR) THEN
               RE=-SCALE*(BEF-BEOLD(I,JF)+(FERP*ARP-FERM*ARM)*CCF)
               CEF1=SCALE*(BEF_PF+(ARP*FERP_PF-ARM*FERM_PF)*CCF)
               CEF2=SCALE*(BEF_XF+(ARP*FERP_XF-ARM*FERM_XF)*CCF)
               CEG1=SCALE*ARP*FERP_PG*CCF
               CEG2=SCALE*ARP*FERP_XG*CCF
            END IF
            RRR=ABS(RE)
            RMAXME=MAX(RMAXME,RRR)

cc  now calculate amended coefficients

            IF(JF==NR) THEN
               RMS=RM
               CME1S=CME1
               CME2S=CME2
               CMF1S=CMF1
               CMF2S=CMF2

cc  energy terms

               RES=RE
               CEE1S=CEE1
               CEE2S=CEE2
               CEF1S=CEF1
               CEF2S=CEF2
            ELSE IF(JF<NR) THEN
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
            END IF
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

        end do
      END IF
cc
cc  Calculate thermodynamic quantities for wellbore block when layer is inactive
cc  Store as  xxE
cc
      IF(IACT(I)==0) THEN
            PE=P(I,1)
            TE=T(I,1)
            P0E=P0(I,1)
            T0E=T0(I,1)
            XE=X(I,1)
            SVE=SV(I,1)
            SLE=1.0-SVE
            IPHE=IPH(I,1)
            VE=V(I,1)
	      PERE=PERMA(I,1)
            PORE=PORWB(I)
            CRE=CRWB
            RHORE=RHORWB
            COMP=COMPWB
            COMT=COMTWB
            AAA=AAAWB
            CALL TERMSMF(PE,TE,SVE,IPHE,PE1,TE1,TE2,BME,BEE,
     1         TMLE,TMVE,TELE,TEVE,BME1,BEE1,TMLE1,TMVE1,TELE1,TEVE1,
     2         BME2,BEE2,TMLE2,TMVE2,TELE2,TEVE2,DELPE,DELXE,
     3         RHOLE,RHOVE,RHOLE1,RHOVE1,RHOLE2,RHOVE2,
     4         PORE,CRE,RHORE,FACP,FACT,FACS,PERE,P0E,T0E,COMP,COMT,AAA,
     5         PERFAC)
	      TME=TMLE+TMVE
            TEE=TELE+TEVE
            TME1=TMLE1+TMVE1
            TEE1=TELE1+TEVE1
            TME2=TMLE2+TMVE2
            TEE2=TELE2+TEVE2
            BME_PE=(BME1-BME)/DELPE
            BME_XE=(BME2-BME)/DELXE
            BEE_PE=(BEE1-BEE)/DELPE
            BEE_XE=(BEE2-BEE)/DELXE
       END IF
cc
cc  wellbore equation for block (I,1)
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
cc  first rename final matrix variables using wellbore notation
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
         TMLC=TMLE*PRATC
         TMVC=TMVE*PRATC
         TELC=TELE*PRATC
         TEVC=TEVE*PRATC
         TMLC1=TMLE1*PRATC
         TMVC1=TMVE1*PRATC
         TELC1=TELE1*PRATC
         TEVC1=TEVE1*PRATC
         TMLC2=TMLE2*PRATC
         TMVC2=TMVE2*PRATC
         TELC2=TELE2*PRATC
         TEVC2=TEVE2*PRATC
         RHOLC=RHOLE
         RHOVC=RHOVE
         RHOLC1=RHOLE1
         RHOVC1=RHOVE1
         RHOLC2=RHOLE2
         RHOVC2=RHOVE2
         PC1=PE1
         TC1=TE1
         TC2=TE2 
	 DELPC=DELPE
	 DELXC=DELXE
cc  
         IF(IB>0) THEN
            AP=AZ(1)
            DZP=DELZ(IB)
         END IF
cc
         COND=CONDWB
cc
cc  only calculate fluxes and set up fracture equations for I>1
cc
         IF(I>=2) THEN
cc
cc  calculate liquid and vapour gradients
cc
            GRADL=(PB-PC)/DZP + 0.5*(RHOLB+RHOLC)*G
            GRADV=(PB-PC)/DZP + 0.5*(RHOVB+RHOVC)*G
            GRADL1C=(PB-PC1)/DZP + 0.5*(RHOLB+RHOLC1)*G
            GRADV1C=(PB-PC1)/DZP + 0.5*(RHOVB+RHOVC1)*G
            GRADL2C=(PB-PC)/DZP + 0.5*(RHOLB+RHOLC2)*G
            GRADV2C=(PB-PC)/DZP + 0.5*(RHOVB+RHOVC2)*G
            GRADL1B=(PB1-PC)/DZP + 0.5*(RHOLB1+RHOLC)*G
            GRADV1B=(PB1-PC)/DZP + 0.5*(RHOVB1+RHOVC)*G
            GRADL2B=(PB-PC)/DZP + 0.5*(RHOLB2+RHOLC)*G
            GRADV2B=(PB-PC)/DZP + 0.5*(RHOVB2+RHOVC)*G
cc
cc  calculate flux terms using upstream weighting
cc  1st liquid
cc
            IF(GRADL<0.0) THEN
cc  basic flows
               FMPL=TMLC*GRADL
               FEPL=TELC*GRADL 
cc  flows with incremented pressure in block C
               FMPL1C=TMLC1*GRADL1C
               FEPL1C=TELC1*GRADL1C 
cc  flows with incremented temperature or saturation in block C
               FMPL2C=TMLC2*GRADL2C
               FEPL2C=TELC2*GRADL2C 
cc  flows with incremented pressure in block B
	         FMPL1B=TMLC*GRADL1B
	         FEPL1B=TELC*GRADL1B 
cc  flows with incremented temperature or saturation in block B
	         FMPL2B=TMLC*GRADL2B
	         FEPL2B=TELC*GRADL2B 
            ELSE
cc  basic flows
               FMPL=TMLB*GRADL
               FEPL=TELB*GRADL 
cc  flows with incremented pressure in block C
               FMPL1C=TMLB*GRADL1C
               FEPL1C=TELB*GRADL1C 
cc  flows with incremented temperature or saturation in block C
               FMPL2C=TMLB*GRADL2C
               FEPL2C=TELB*GRADL2C 
cc  flows with incremented pressure in block B
	         FMPL1B=TMLB1*GRADL1B
	         FEPL1B=TELB1*GRADL1B 
cc  flows with incremented temperature or saturation in block B
	         FMPL2B=TMLB2*GRADL2B
	         FEPL2B=TELB2*GRADL2B
            END IF
cc
cc  2nd vapour
cc
            IF(GRADV<0.0) THEN
cc  basic flows
               FMPV=TMVC*GRADV
               FEPV=TEVC*GRADV 
cc  flows with incremented pressure in block C
               FMPV1C=TMVC1*GRADV1C
               FEPV1C=TEVC1*GRADV1C 
cc  flows with incremented temperature or saturation in block C
               FMPV2C=TMVC2*GRADV2C
               FEPV2C=TEVC2*GRADV2C 
cc  flows with incremented pressure in block B
	         FMPV1B=TMVC*GRADV1B
	         FEPV1B=TEVC*GRADV1B 
cc  flows with incremented temperature or saturation in block B
	         FMPV2B=TMVC*GRADV2B
	         FEPV2B=TEVC*GRADV2B 
            ELSE
cc  basic flows
               FMPV=TMVB*GRADV
               FEPV=TEVB*GRADV 
cc  flows with incremented pressure in block C
               FMPV1C=TMVB*GRADV1C
               FEPV1C=TEVB*GRADV1C 
cc  flows with incremented temperature or saturation in block C
               FMPV2C=TMVB*GRADV2C
               FEPV2C=TEVB*GRADV2C 
cc  flows with incremented pressure in block B
	         FMPV1B=TMVB1*GRADV1B
	         FEPV1B=TEVB1*GRADV1B 
cc  flows with incremented temperature or saturation in block B
	         FMPV2B=TMVB2*GRADV2B
	         FEPV2B=TEVB2*GRADV2B 
            END IF

               FMP=FMPL+FMPV
               FEP=FEPL+FEPV-COND*(TC-TB)/DZP
               FMP1C=FMPL1C+FMPV1C
               FEP1C=FEPL1C+FEPV1C-COND*(TC1-TB)/DZP
               FMP2C=FMPL2C+FMPV2C
               FEP2C=FEPL2C+FEPV2C-COND*(TC2-TB)/DZP
               FMP1B=FMPL1B+FMPV1B
               FEP1B=FEPL1B+FEPV1B-COND*(TC-TB1)/DZP
               FMP2B=FMPL2B+FMPV2B
               FEP2B=FEPL2B+FEPV2B-COND*(TC-TB2)/DZP
cc  mass flux derivatives
               FMP_PC=(FMP1C - FMP)/DELPC
               FMP_XC=(FMP2C - FMP)/DELXC
               FMP_PB=(FMP1B - FMP)/DELPB
               FMP_XB=(FMP2B - FMP)/DELXB
cc  energy flux derivatives
               FEP_PC=(FEP1C - FEP)/DELPC
               FEP_XC=(FEP2C - FEP)/DELXC
               FEP_PB=(FEP1B - FEP)/DELPB
               FEP_XB=(FEP2B - FEP)/DELXB
               FMZ(IB)=FMP*AP
               FEZ(IB)=FEP*AP
               FMX(IB)=FMQ*AQ

cc
cc  set up mass coefficients and mass residual for (I-1)th equation
cc
            CCB=DT/VB
            IF(IB==1) THEN
               IF(IACT(IB)==1) THEN
                  RMB=-(BMB-BMOLD(IB,1)+(FMP*AP+FMQ*AQ-QMM)*CCB)
                  CMB1=BMB_PB+(FMP_PB*AP+FMQ_PB*AQ)*CCB
                  CMB2=BMB_XB+(FMP_XB*AP+FMQ_XB*AQ)*CCB
               ELSE
                  RMB=-(BMB-BMOLD(IB,1)+(FMP*AP-QMM)*CCB)
                  CMB1=BMB_PB+(FMP_PB*AP)*CCB
                  CMB2=BMB_XB+(FMP_XB*AP)*CCB
               END IF
               CMA1=0.0
               CMA2=0.0
		     DELBM=BMB-BMOLD(IB,1)
	         DELFM=FMP*AP-FMM*AM
            ELSE
               IF(IACT(IB)==1) THEN
                  RMB=-(BMB-BMOLD(IB,1)+(FMP*AP+FMQ*AQ-FMM*AM)*CCB)
                  CMB1=BMB_PB+(FMP_PB*AP+FMQ_PB*AQ-FMM_PB*AM)*CCB
                  CMB2=BMB_XB+(FMP_XB*AP+FMQ_XB*AQ-FMM_XB*AM)*CCB
               ELSE
                  RMB=-(BMB-BMOLD(IB,1)+(FMP*AP-FMM*AM)*CCB)
                  CMB1=BMB_PB+(FMP_PB*AP-FMM_PB*AM)*CCB
                  CMB2=BMB_XB+(FMP_XB*AP-FMM_XB*AM)*CCB
               END IF
               CMA1=-AM*FMM_PA*CCB
               CMA2=-AM*FMM_XA*CCB
	         DELBM=BMB-BMOLD(IB,1)
	         DELFM=FMP*AP-FMM*AM
            ENDIF
            CMC1=FMP_PC*AP*CCB
            CMC2=FMP_XC*AP*CCB
            L=2*IB-1
            IF(IACT(IB)==1) THEN
               CMD1=FMQ_PD*AQ*CCB
               CMD2=FMQ_XD*AQ*CCB
               R(L)=RMB-CMD1*ROM(IB,2)-CMD2*ROE(IB,2)
               C4(L)=CMB1-CMD1*ALM1(IB,1)-CMD2*ALE1(IB,1)
               C5(L)=CMB2-CMD1*ALM2(IB,1)-CMD2*ALE2(IB,1)
            ELSE
               CMD1=0.0
               CMD2=0.0
               R(L)=RMB
               C4(L)=CMB1
               C5(L)=CMB2
            ENDIF
            RR(L)=R(L)
            RRR=ABS(RMB)
            RMAXFM=MAX(RMAXFM,RRR)
            C1(L)=0.0
            C2(L)=CMA1
            C3(L)=CMA2
            C6(L)=CMC1
            C7(L)=CMC2
cc
cc  energy equation coefficients and energy residual
            IF(IB==1) THEN
cc  production/injection term
               IF(QMM>=0.0) THEN
                  QEM=QMM*HIN
                  QEM_PB=0.0
                  QEM_XB=0.0
               ELSE
                  QEM=QMM*(TEVB+TELB)/(TMVB+TMLB)
                  QEM1=QMM*(TEVB1+TELB1)/(TMVB1+TMLB1)
                  QEM2=QMM*(TEVB2+TELB2)/(TMVB2+TMLB2)
                  QEM_PB=(QEM1 - QEM)/DELPB
                  QEM_XB=(QEM2 - QEM)/DELXB
               END IF
cc
               IF(IACT(IB)==1) THEN
                REB=-SCALE*(BEB-BEOLD(IB,1)+(FEP*AP+FEQ*AQ-QEM)*CCB)
                CEB1=SCALE*(BEB_PB+(FEP_PB*AP+FEQ_PB*AQ-QEM_PB)*CCB)
                CEB2=SCALE*(BEB_XB+(FEP_XB*AP+FEQ_XB*AQ-QEM_XB)*CCB)
               ELSE
                REB=-SCALE*(BEB-BEOLD(IB,1)+(FEP*AP-QEM)*CCB)
                CEB1=SCALE*(BEB_PB+(FEP_PB*AP-QEM_PB)*CCB)
                CEB2=SCALE*(BEB_XB+(FEP_XB*AP-QEM_XB)*CCB)
               ENDIF
               CEA1=0.0
               CEA2=0.0
            ELSE IF (IB>1) THEN
               IF(IACT(IB)==1) THEN
                 REB=-SCALE*(BEB-BEOLD(IB,1)+(FEP*AP+FEQ*AQ-FEM*AM)*CCB)
                 CEB1=SCALE*(BEB_PB+(FEP_PB*AP+FEQ_PB*AQ-FEM_PB*AM)*CCB)
                 CEB2=SCALE*(BEB_XB+(FEP_XB*AP+FEQ_XB*AQ-FEM_XB*AM)*CCB)
               ELSE
                 REB=-SCALE*(BEB-BEOLD(IB,1)+(FEP*AP-FEM*AM)*CCB)
                 CEB1=SCALE*(BEB_PB+(FEP_PB*AP-FEM_PB*AM)*CCB)
                 CEB2=SCALE*(BEB_XB+(FEP_XB*AP-FEM_XB*AM)*CCB)
               END IF
               CEA1=-SCALE*AM*FEM_PA*CCB
               CEA2=-SCALE*AM*FEM_XA*CCB
            END IF
            CEC1=SCALE*(FEP_PC*AP*CCB)
            CEC2=SCALE*(FEP_XC*AP*CCB)
            L1=2*IB
            IF(IACT(IB)==1) THEN
               CED1=SCALE*(FEQ_PD*AQ*CCB)
               CED2=SCALE*(FEQ_XD*AQ*CCB)
               R(L1)=REB-CED1*ROM(IB,2)-CED2*ROE(IB,2)
               C3(L1)=CEB1-CED1*ALM1(IB,1)-CED2*ALE1(IB,1)
               C4(L1)=CEB2-CED1*ALM2(IB,1)-CED2*ALE2(IB,1)
            ELSE
               CED1=0.0
               CED2=0.0
               R(L1)=REB
               C3(L1)=CEB1
               C4(L1)=CEB2
            END IF
            R(L1)=REB-CED1*ROM(IB,2)-CED2*ROE(IB,2)
            RR(L1)=R(L1)
            RRR=ABS(REB)
            RMAXFE=MAX(RMAXFE,RRR)
            C1(L1)=CEA1
            C2(L1)=CEA2
            C5(L1)=CEC1
            C6(L1)=CEC2
            C7(L1)=0.0
         END IF
cc
cc  now shuffle terms along ready for next calculation at (I+1)th block
cc
         IF(I>1) THEN
            PA=PB
            TA=TB
            XA=XB
            SVA=SVB
            SLA=SLB
            IPHA=IPHB
            VA=VB
            AM=AP
            DZM=DZP
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
         END IF
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
	   TMLB=TMLC
	   TELB=TELC
	   TMVB=TMVC
	   TEVB=TEVC
	   TMLB1=TMLC1
	   TELB1=TELC1
	   TMVB1=TMVC1
	   TEVB1=TEVC1
	   TMLB2=TMLC2
	   TELB2=TELC2
	   TMVB2=TMVC2
	   TEVB2=TEVC2
	   RHOLB=RHOLC
         RHOVB=RHOVC
	   RHOLB1=RHOLC1
         RHOVB1=RHOVC1
	   RHOLB2=RHOLC2
         RHOVB2=RHOVC2
cc  shuffle along derived quanitities
         BMB_PB=BMC_PC 
         BMB_XB=BMC_XC     
         BEB_PB=BEC_PC     
         BEB_XB=BEC_XC
cc  
cc  calculate radial flow and other terms required for Ith mass and
cc  energy equations which are set up at next (I+1th) step
cc
         IF(IACT(I)==1) THEN
            AQ=ARM
            FMQ=FMRM
            FEQ=FERM
            FMQ_PD=FMRM_PF
            FMQ_PB=FMRM_PE
            FMQ_XD=FMRM_XF
            FMQ_XB=FMRM_XE
            FEQ_PD=FERM_PF
            FEQ_PB=FERM_PE
            FEQ_XD=FERM_XF
            FEQ_XB=FERM_XE
         ELSE
            AQ=0.0
            FMQ=0.0
            FEQ=0.0
            FMQ_PD=0.0
            FMQ_PB=0.0
            FMQ_XD=0.0
            FMQ_XB=0.0
            FEQ_PD=0.0
            FEQ_PB=0.0
            FEQ_XD=0.0
            FEQ_XB=0.0
         ENDIF
cc
cc  now complete the equation for the end block if I=MZ
cc  
         IF(I==MZ) THEN
            IB=MZ
            CCB=DT/VB
            L=2*IB-1
            IF(IACT(IB)==1) THEN
               RMB=-(BMB-BMOLD(IB,1)+(FMQ*AQ-FMM*AM)*CCB)
               CMB1=BMB_PB+(FMQ_PB*AQ-FMM_PB*AM)*CCB
               CMB2=BMB_XB+(FMQ_XB*AQ-FMM_XB*AM)*CCB
               CMD1=FMQ_PD*AQ*CCB
               CMD2=FMQ_XD*AQ*CCB
               R(L)=RMB-CMD1*ROM(IB,2)-CMD2*ROE(IB,2)
               C4(L)=CMB1-CMD1*ALM1(IB,1)-CMD2*ALE1(IB,1)
               C5(L)=CMB2-CMD1*ALM2(IB,1)-CMD2*ALE2(IB,1)
            ELSE
               RMB=-(BMB-BMOLD(IB,1)+(-FMM*AM)*CCB)
               CMB1=BMB_PB+(-FMM_PB*AM)*CCB
               CMB2=BMB_XB+(-FMM_XB*AM)*CCB
               CMD1=0.0
               CMD2=0.0
               R(L)=RMB
               C4(L)=CMB1
               C5(L)=CMB2
            END IF
            CMA1=-AM*FMM_PA*CCB
            CMA2=-AM*FMM_XA*CCB
            CMC1=0.0
            CMC2=0.0
            RR(L)=R(L)
            RRR=ABS(RMB)
            RMAXFM=MAX(RMAXFM,RRR)
            C1(L)=0.0
            C2(L)=CMA1
            C3(L)=CMA2
            C6(L)=CMC1
            C7(L)=CMC2
cc
cc  energy equation coefficients and energy residual
cc
            L1=2*IB
            IF(IACT(IB)==1) THEN
               REB=-SCALE*(BEB-BEOLD(IB,1)+(FEQ*AQ-FEM*AM)*CCB)
               CEB1=SCALE*(BEB_PB+(FEQ_PB*AQ-FEM_PB*AM)*CCB)
               CEB2=SCALE*(BEB_XB+(FEQ_XB*AQ-FEM_XB*AM)*CCB)
               CED1=SCALE*(FEQ_PD*AQ*CCB)
               CED2=SCALE*(FEQ_XD*AQ*CCB)
               R(L1)=REB-CED1*ROM(IB,2)-CED2*ROE(IB,2)
               C3(L1)=CEB1-CED1*ALM1(IB,1)-CED2*ALE1(IB,1)
               C4(L1)=CEB2-CED1*ALM2(IB,1)-CED2*ALE2(IB,1)
            ELSE
               REB=-SCALE*(BEB-BEOLD(IB,1)+(-FEM*AM)*CCB)
               CEB1=SCALE*(BEB_PB+(-FEM_PB*AM)*CCB)
               CEB2=SCALE*(BEB_XB+(-FEM_XB*AM)*CCB)
               CED1=0.0
               CED2=0.0
               R(L1)=REB
               C3(L1)=CEB1
               C4(L1)=CEB2
            END IF
            CEA1=-SCALE*AM*FEM_PA*CCB
            CEA2=-SCALE*AM*FEM_XA*CCB
            CEC1=0.0
            CEC2=0.0
            RR(L1)=R(L1)
            RRR=ABS(REB)
            RMAXFE=MAX(RMAXFE,RRR)
            C1(L1)=CEA1
            C2(L1)=CEA2
            C5(L1)=CEC1
            C6(L1)=CEC2
            C7(L1)=0.0
	      Z1=BMB
	      Z2=BMOLD(IB,1)
	      Z3=FMM
	      Z4=AM
	      Z5=VB
	      Z6=DZM
         ENDIF
      end do

cc
cc  now solve the amended wellbore equations
cc
      M2=2*MZ
      CALL SEVEN1(C1,C2,C3,C4,C5,C6,C7,R,XXX,M2)
!      do I=1,M2
!         XX(I,1)=XXX(I)
!      end do
	XX(1:M2,1)=XXX(1:M2)
cc
cc  now back substitute to calculate the matrix unknowns
cc
      do I=1,MZ
         if (IACT(I)==1) then
	      L=2*I-1
            L1=2*I
            do J=2,NR
              JA=J-1
              XX(L,J)=ROM(I,J)-ALM1(I,JA)*XX(L,JA)-ALM2(I,JA)*XX(L1,JA)
              XX(L1,J)=ROE(I,J)-ALE1(I,JA)*XX(L,JA)-ALE2(I,JA)*XX(L1,JA)
           end do
         end if
      end do

!     Added trap for NaN values in XX array (AC 3/03):
	if (any(isnan(XX))) then
	  IGOOD=1
	end if

      RETURN
      END SUBROUTINE SOLVEMF
cc       
cc----7--1--------2---------3---------4---------5---------6---------7---------8 

      subroutine UPDATE(NR,MZ,IPH,P,T,SV,X,XX,EPS)

      implicit none
	integer(I4B):: NR,MZ     
      real(DP):: P(MZ,NR),T(MZ,NR),X(MZ,NR),SV(MZ,NR),XX(2*MZ,NR)
	real(DP):: EPS,PP,TP,PT,XT,TT,PS,SVT,T0
      integer(I4B):: IPH(MZ,NR)
	integer(I4B):: I,J,L,L1,IPHP

cc  This subroutine adds the increments on to the primary variables and carries out
cc  phase changes

      do I=1,MZ
         do J=1,NR
            L=2*I-1
            L1=L+1
            PP=P(I,J)
            TP=T(I,J)
            IPHP=IPH(I,J)
            PT=PP+XX(L,J)
            XT=X(I,J)+XX(L1,J)

cc  -----------------------------------------------------------------------------
            select case (IPHP)
            case (0)  
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
cc  currently 2-phase
               SVT=XT
               CALL SAT(TP,PS)
               IF((SVT>0.0d0).AND.(SVT<1.0d0)) THEN
cc  remain 2-phase
                  P(I,J)=PT
                  SV(I,J)=SVT
                  X(I,J)=SVT
                  IPH(I,J)=1
	            T0=0.0d0
                  CALL TSAT(PT,T0,TT)
	            if (IGOOD>0) return 
                  T(I,J)=TT
               ELSE IF(SVT<=0.0d0) THEN
cc  change to single phase liquid water
                  P(I,J)=1.000001*PS
                  T(I,J)=TP
                  X(I,J)=TP
                  SV(I,J)=0.0
                  IPH(I,J)=0
               ELSE IF(SVT>=1.0) THEN
cc  change to single phase dry steam
                  P(I,J)=0.999999*PS
                  T(I,J)=TP
                  X(I,J)=TP
                  SV(I,J)=1.0
                  IPH(I,J)=2
               END IF
cc
cc  -----------------------------------------------------------------------------
            case (2)
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
      return
      end subroutine UPDATE
cc       
cc----7--1--------2---------3---------4---------5---------6---------7---------8 

      subroutine INIT1(NR,MZ,POLD,XOLD,TOLD,SVOLD,IPHOLD)

      implicit none       
	integer(I4B):: NR,MZ
      real(DP):: POLD(MZ,NR),XOLD(MZ,NR),TOLD(MZ,NR)
      real(DP):: SVOLD(MZ,NR)
	integer(I4B):: IPHOLD(MZ,NR)
      real(DP):: P(MZ,NR),T(MZ,NR),SV(MZ,NR)
	real(DP):: PP,XP,TZERO,TS,TP,SVP
	integer(I4B):: I,J,IPHP

      do I=1,MZ
         do J=1,NR
            PP=POLD(I,J)
            XP=XOLD(I,J)
	      CALL TSAT(PP,TZERO,TS)
            IF(XP<1.0d0) THEN
cc  2-phase conditions
               TZERO=0.0d0
               TP=TS 
               TOLD(I,J)=TS
               SVOLD(I,J)=XP
               SVP=XP
               IPHOLD(I,J)=1
               IPHP=1
            ELSEIF(XP>1.5d0) THEN
               IF(XP>=TS) THEN
cc   single phase dry steam
                  TP=XP
                  TOLD(I,J)=XP
                  SVOLD(I,J)=1.0d0
                  SVP=1.0d0
                  IPHOLD(I,J)=2
                  IPHP=2
               ELSE
cc   single phase liquid water
                  TP=XP
                  TOLD(I,J)=XP
                  SVOLD(I,J)=0.0d0
                  SVP=0.0d0
                  IPHOLD(I,J)=0
                  IPHP=0
               ENDIF
            ENDIF
         end do
      end do
      return
      end subroutine INIT1
cc       
cc----7--1--------2---------3---------4---------5---------6---------7---------8 
cc
      subroutine INIT2(NR,MZ,POLD,XOLD,TOLD,SVOLD,IPHOLD,
     1                  BMOLD,BEOLD,PORWB,RHORWB,CRWB,PORMA,RHORMA,CRMA,
     2                  P0,T0,COMPWB,COMTWB,AAAWB,COMPMA,COMTMA,AAAMA)

      implicit none
	integer(I4B):: NR,MZ    
      real(DP):: POLD(MZ,NR),XOLD(MZ,NR),TOLD(MZ,NR)
      real(DP):: P0(MZ,NR),T0(MZ,NR)
      real(DP):: SVOLD(MZ,NR)
	integer(I4B):: IPHOLD(MZ,NR)
      real(DP):: BMOLD(MZ,NR),BEOLD(MZ,NR)
      real(DP):: PORWB(MZ),PORMA(MZ,NR)
      real(DP):: RHORMA(MZ),COMPMA(MZ),COMTMA(MZ),AAAMA(MZ)
	real(DP):: CRMA(MZ),CONDMA(MZ)
	real(DP):: RHORWB,CRWB,COMPWB,COMTWB,AAAWB,PP,TP,P0P,T0P
	real(DP):: SVP,PORP,RHOR,CR,COMP,COMT,AAA,RHOL,RHOV,HL,HV
	real(DP):: VISL,VISV,BM,BE,PERFAC
	integer(I4B):: I,J,IPHP
      
	do I=1,MZ
         do J=1,NR
            PP=POLD(I,J)
            TP=TOLD(I,J)
            P0P=P0(I,J)
            T0P=T0(I,J)
	      SVP=SVOLD(I,J)
	      IF(J==1) THEN
               PORP=PORWB(I)
	         RHOR=RHORWB
	         CR=CRWB
               COMP=COMPWB
               COMT=COMTWB
               AAA=AAAWB
            ELSE
               PORP=PORMA(I,J)
     	         RHOR=RHORMA(I)
	         CR=CRMA(I)
               COMP=COMPMA(I)
               COMT=COMTMA(I)
               AAA=AAAMA(I)
            END IF
            IPHP=IPHOLD(I,J)
            CALL THERMO(PP,TP,SVP,IPHP,RHOL,RHOV,HL,HV,VISL,VISV)
            CALL ACCUM(BM,BE,IPHP,PP,TP,SVP,RHOL,RHOV,HL,HV,PORP,
     1                 RHOR,CR,P0P,T0P,COMP,COMT,AAA,PERFAC)
            BMOLD(I,J)=BM
            BEOLD(I,J)=BE
         end do
      end do
      return
      end subroutine INIT2
cc       
cc----7--1--------2---------3---------4---------5---------6---------7---------8 

      end module NumericalSimulatorMultiFeedzone_Routines