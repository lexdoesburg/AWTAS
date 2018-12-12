module MINCSimulator

  use variable_types
  use problem_data
  implicit none

  contains

! ------------------------------------------------------------------------------------

  function MINC(variable,updatemodelprogress)

!   Numerical simulator for a radially-symmetric MINC model. 

!   Fixed parameter values for this model:
!   1: Layer elevation (not used)
!   2: Layer thickness
!   3: Action well radius
!   4: Fracture volume fraction
!   5: Fracture spacing
!   6: # fracture planes (1..3)
!   7: Fractured rock specific heat
!   8: Matrix rock specific heat
!   9: Fractured rock heat conductivity
!  10: Matrix rock heat conductivity
!  11: Fractured rock density
!  12: Matrix rock density
!  13: Fractured rock compressibility
!  14: Matrix rock compressibility

!  Not yet implemented: 
!  15: Fractured rock thermal expansion coefficient
!  16: Matrix rock thermal expansion coefficient
!  17: Fractured rock permeability/porosity correlation coefficient
!  18: Matrix rock permeability/porosity correlation coefficient 

!   Reservoir conditions for this model:
!   1: P (pressure)
!   2: X (temperature/ vapour saturation)

!   Variables for this model:
!   1: kFracture (fracture permeability)
!   2: phiFracture (fracture porosity)
!   3: kMatrix (matrix permeability)
!   4: phiMatrix (matrix porosity)

	use NumericalSimulatorMINC

!   Argument Variables:
    real(DP),intent(in) :: variable(:)
    real(DP) :: MINC(TotalNData)
!   Fixed parameters:	                    
	real(DP)    :: LayerThickness,ActionWellRadius
	real(DP)    :: FractureVolFraction,FractureSpacing
    ! NB: this is generally an integer, but it's easier to 
    ! leave it real, & allows for generalisation to 'fractional dimension' MINC:
    real(DP)    :: NumFracturePlanes
	real(DP)    :: CRF,CRM,CONDF,CONDM,RHORF,RHORM
	real(DP)    :: COMPF,COMPM,COMTF,COMTM,AAAF,AAAM 
!   Reservoir conditions:
    real(DP)    :: InitialPressure,InitialX
!   Variable parameters:
    real(DP)    :: kFracture,kMatrix,phiFracture,phiMatrix

    external updatemodelprogress

!   Unpack fixed parameters:
    LayerThickness=FixedParameter(2)
	ActionWellRadius=FixedParameter(3)
	FractureVolFraction=FixedParameter(4)
    FractureSpacing=FixedParameter(5)
    NumFracturePlanes=FixedParameter(6)
    CRF=FixedParameter(7)     
    CRM=FixedParameter(8)     
	CONDF=FixedParameter(9)   
	CONDM=FixedParameter(10)   
	RHORF=FixedParameter(11)   
	RHORM=FixedParameter(12)  
	COMPF=FixedParameter(13) 
	COMPM=FixedParameter(14) 
	! Variable k/phi has been implemented in the back end, but not the front:
	COMTF=0.0_dp ! FixedParameter(15)
    COMTM=0.0_dp ! FixedParameter(16)
	AAAF=0.0_dp  ! FixedParameter(17)
	AAAM=0.0_dp  ! FixedParameter(18)
!   Unpack reservoir conditions:
    InitialPressure=ReservoirCondition(1) 
	InitialX=ReservoirCondition(2) 
!   Unpack variable parameters:
    kFracture=variable(1)   
    phiFracture=variable(2) 
    kMatrix=variable(3)   
    phiMatrix=variable(4) 

	call SetupGrid(LayerThickness,ActionWellRadius,FractureVolFraction,&
	  FractureSpacing,NumFracturePlanes)
	  
    MINC=NumericalSolutionMINC(CRF,CRM,CONDF,CONDM,RHORF,RHORM,&
      COMPF,COMPM,COMTF,COMTM,AAAF,AAAM,InitialPressure,InitialX,&
	  kFracture,kMatrix,phiFracture,phiMatrix,updatemodelprogress)

    return    
  end function MINC 

! ------------------------------------------------------------------------

end module MINCSimulator