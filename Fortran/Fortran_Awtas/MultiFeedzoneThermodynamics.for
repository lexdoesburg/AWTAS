      module MultiFeedzoneThermodynamics

	use Thermodynamics

      contains

cc----7--1--------2---------3---------4---------5---------6---------7---------8 
cc        
      SUBROUTINE TRANSMF(TML,TMV,TEL,TEV,IPH,SV,
     1                 HL,HV,VISL,VISV,PER,PERFAC)
      IMPLICIT REAL*8 (A-H,O-Z)       
cc
cc  This subroutine calculates the mass and energy transmissibility terms TM and TE
cc  It requires values of SV,HL,HV,VISL,VISV. It calculates relative permeabilities if necessary
cc      
      select case (IPH)
cc  single phase liquid
      case (0)
         TML=PER*PERFAC*1.0D0/VISL
         TMV=0.0
         TEL=HL*TML
         TEV=0.0
cc  two-phase
      case (1)
         CALL RELP(SV,IRELP,XKRL,XKRV)
         TML=PER*PERFAC*XKRL/VISL 
         TEL=HL*TML 
         TMV=PER*PERFAC*XKRV/VISV
         TEV=HV*TMV
cc  single phase steam
      case (2)
         TML=0.0D0
         TMV=PER*PERFAC*1.0D0/VISV
         TEL=0.0D0
         TEV=HV*TMV
      end select
      RETURN
      END SUBROUTINE TRANSMF
cc       
cc----7--1--------2---------3---------4---------5---------6---------7---------8 
cc        
      SUBROUTINE TERMSMF(PX,TX,SVX,IPHX,PX1,TX1,TX2,BMX,BEX,
     1      TMLX,TMVX,TELX,TEVX,BMX1,BEX1,TMLX1,TMVX1,TELX1,TEVX1,
     2      BMX2,BEX2,TMLX2,TMVX2,TELX2,TEVX2,DELPX,DELXX,
     3      RHOLX,RHOVX,RHOLX1,RHOVX1,RHOLX2,RHOVX2,PORX,CRX,
     4      RHORX,FACP,FACT,FACS,PERX,P0X,T0X,COMP,COMT,AAA,PERFAC)
      IMPLICIT REAL*8 (A-H,O-Z)
cc
cc  This subroutine calculates accumulation and transmissibility terms at each block
cc  It also perturbs the primary variables and recalculates them
cc
cc  First calculate the basic terms
cc
      CALL THERMO(PX,TX,SVX,IPHX,RHOL,RHOV,HL,HV,VISL,VISV)
	if (IGOOD>0) return
      RHOLX=RHOL
      RHOVX=RHOV
      CALL ACCUM(BMX,BEX,IPHX,PX,TX,SVX,
     1	       RHOL,RHOV,HL,HV,PORX,RHORX,CRX,
     2         P0X,T0X,COMP,COMT,AAA,PERFAC)
      CALL TRANSMF(TMLX,TMVX,TELX,TEVX,IPHX,SVX,
     1           HL,HV,VISL,VISV,PERX,PERFAC)
cc
cc  now perturb the primary variables to calculate derivatives
cc
cc  -------------------------------------------------------------------------------
cc  single phase liquid
cc
      select case (IPHX)

	case (0)
cc
cc  make sure it stays liquid - increase pressure
cc
         DELPX=FACP*PX
         PX1=PX+DELPX
	   TX1=TX
         CALL THERMO(PX1,TX,SVX,IPHX,RHOL1,RHOV1,HL1,HV1,VISL1,VISV1)
	   if (IGOOD>0) return
         RHOLX1=RHOL1
         RHOVX1=RHOV1
         CALL ACCUM(BMX1,BEX1,IPHX,PX1,TX,SVX,
     1             RHOL1,RHOV1,HL1,HV1,PORX,RHORX,CRX,
     2             P0X,T0X,COMP,COMT,AAA,PERFAC1)
         CALL TRANSMF(TMLX1,TMVX1,TELX1,TEVX1,IPHX,SVX,
     1             HL1,HV1,VISL1,VISV1,PERX,PERFAC1)
cc
cc  make sure it stays liquid - decrease temperature
         DELXX=-FACT*TX
         TX2=TX+DELXX
         CALL THERMO(PX,TX2,SVX,IPHX,RHOL2,RHOV2,HL2,HV2,VISL2,VISV2)
         if (IGOOD>0) return
         RHOLX2=RHOL2
         RHOVX2=RHOV2
         CALL ACCUM(BMX2,BEX2,IPHX,PX,TX2,SVX,
     1             RHOL2,RHOV2,HL2,HV2,PORX,RHORX,CRX,
     2             P0X,T0X,COMP,COMT,AAA,PERFAC2)
         CALL TRANSMF(TMLX2,TMVX2,TELX2,TEVX2,IPHX,SVX,
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
	   if (IGOOD>0) return
         CALL THERMO(PX1,TX1,SVX,IPHX,RHOL1,RHOV1,HL1,HV1,VISL1,VISV1)
	   if (IGOOD>0) return
         RHOLX1=RHOL1
         RHOVX1=RHOV1
         CALL ACCUM(BMX1,BEX1,IPHX,PX1,TX1,SVX,
     1             RHOL1,RHOV1,HL1,HV1,PORX,RHORX,CRX,
     2             P0X,T0X,COMP,COMT,AAA,PERFAC1)
         CALL TRANSMF(TMLX1,TMVX1,TELX1,TEVX1,IPHX,SVX,
     1             HL1,HV1,VISL1,VISV1,PERX,PERFAC1)
cc  increase saturation
         DELXX=FACS*SVX+1.0D-10
         IF(SVX+DELXX>=1.0D0) DELXX=-DELXX
         SVX2=SVX+DELXX
         RHOLX2=RHOL
         RHOVX2=RHOV
         CALL ACCUM(BMX2,BEX2,IPHX,PX,TX,SVX2,
     1             RHOL,RHOV,HL,HV,PORX,RHORX,CRX,
     2             P0X,T0X,COMP,COMT,AAA,PERFAC2)
         CALL TRANSMF(TMLX2,TMVX2,TELX2,TEVX2,IPHX,SVX2,
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
         RHOLX1=RHOL1
         RHOVX1=RHOV1
         CALL ACCUM(BMX1,BEX1,IPHX,PX1,TX,SVX,
     1             RHOL1,RHOV1,HL1,HV1,PORX,RHORX,CRX,
     2             P0X,T0X,COMP,COMT,AAA,PERFAC1)
         CALL TRANSMF(TMLX1,TMVX1,TELX1,TEVX1,IPHX,SVX,
     1             HL1,HV1,VISL1,VISV1,PERX,PERFAC1)
cc
cc  make sure it stays gas - increase temperature
         DELXX=FACT*TX
         TX2=TX+DELXX
         CALL THERMO(PX,TX2,SVX,IPHX,RHOL2,RHOV2,HL2,HV2,VISL2,VISV2)
         RHOLX2=RHOL2
         RHOVX2=RHOV2
         CALL ACCUM(BMX2,BEX2,IPHX,PX,TX2,SVX,
     1             RHOL2,RHOV2,HL2,HV2,PORX,RHORX,CRX,
     2             P0X,T0X,COMP,COMT,AAA,PERFAC2)
         CALL TRANSMF(TMLX2,TMVX2,TELX2,TEVX2,IPHX,SVX,
     1             HL2,HV2,VISL2,VISV2,PERX,PERFAC2)
      end select
      RETURN
      END SUBROUTINE TERMSMF
cc       
cc----7--1--------2---------3---------4---------5---------6---------7---------8 


      end module MultiFeedzoneThermodynamics