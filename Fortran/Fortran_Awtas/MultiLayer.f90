module MultiLayerSimulator

  use variable_types
  use problem_data
  implicit none

  contains

! ------------------------------------------------------------------------------------

  function MultiLayer(variable,updatemodelprogress)

!   Numerical simulator for a radially-symmetric multi-layer model. 

!   Fixed parameter values for this model:
!   1: Layer elevation (not used)
!   2: Layer thickness
!   3: Action well radius
!   4: Fracture thickness
!   5: Fracture spacing
!   6: Fractured rock specific heat
!   7: Matrix rock specific heat
!   8: Fractured rock heat conductivity
!   9: Matrix rock heat conductivity
!  10: Fractured rock density
!  11: Matrix rock density
!  12: Fractured rock compressibility
!  13: Matrix rock compressibility 
!  14: Fractured rock thermal expansion coefficient
!  15: Matrix rock thermal expansion coefficient
!  16: Fractured rock permeability/porosity correlation coefficient
!  17: Matrix rock permeability/porosity correlation coefficient 

!   Reservoir conditions for this model:
!   1: P (pressure)
!   2: X (temperature/ vapour saturation)

!   Variables for this model:
!   1: kFracture (fracture permeability)
!   2: phiFracture (fracture porosity)
!   3: kMatrix (matrix permeability)
!   4: phiMatrix (matrix porosity)

	use NumericalSimulatorMultiLayer

!   Argument Variables:
    real(DP),intent(in) :: variable(:)
    real(DP) :: MultiLayer(TotalNData)
!   Fixed parameters:	                    
	real(DP)    :: LayerThickness,ActionWellRadius
	real(DP)    :: FractureThickness,FractureSpacing
	real(DP)    :: CRF,CRM,CONDF,CONDM,RHORF,RHORM
	real(DP)    :: COMPF,COMPM,COMTF,COMTM,AAAF,AAAM 
!   Reservoir conditions:
    real(DP)    :: InitialPressure,InitialX
!   Variable parameters:
    real(DP)    :: kFracture,kMatrix,phiFracture,phiMatrix
!   Locals:
    integer(I4B)::NFractures

    external updatemodelprogress

!   Unpack fixed parameters:
    LayerThickness=FixedParameter(2)
	ActionWellRadius=FixedParameter(3) 
	FractureThickness=FixedParameter(4)
    FractureSpacing=FixedParameter(5)
    CRF=FixedParameter(6)     
    CRM=FixedParameter(7)     
	CONDF=FixedParameter(8)   
	CONDM=FixedParameter(9)   
	RHORF=FixedParameter(10)   
	RHORM=FixedParameter(11)  
	COMPF=FixedParameter(12) 
	COMPM=FixedParameter(13) 
	! Variable k/phi has been implemented in the back end, but not the front:
	COMTF=0.0_dp ! FixedParameter(14)
    COMTM=0.0_dp ! FixedParameter(15)
	AAAF=0.0_dp  ! FixedParameter(16)
	AAAM=0.0_dp  ! FixedParameter(17)
!   Unpack reservoir conditions:
    InitialPressure=ReservoirCondition(1) 
	InitialX=ReservoirCondition(2) 
!   Unpack variable parameters:
    kFracture=variable(1)   
    phiFracture=variable(2) 
    kMatrix=variable(3)   
    phiMatrix=variable(4) 

    NFractures=round(LayerThickness/FractureSpacing)

	call SetupGrid(LayerThickness,ActionWellRadius,FractureThickness,&
	  FractureSpacing)
	  
    MultiLayer=NumericalSolutionMultiLayer(CRF,CRM,CONDF,CONDM,RHORF,RHORM,&
      COMPF,COMPM,COMTF,COMTM,AAAF,AAAM,NFractures,InitialPressure,InitialX,&
	  kFracture,kMatrix,phiFracture,phiMatrix,updatemodelprogress)

    return    
  end function MultiLayer 

! ------------------------------------------------------------------------

  function round(x)
!   Rounding to nearest integer.
!   Argument:
    real(DP),intent(in)::x
	integer(I4B)::round
	round=int(x+0.5_dp)
	return
  end function round

! ------------------------------------------------------------------------

end module MultiLayerSimulator