module LeakySimulator

  use variable_types
  use problem_data
  implicit none

  contains

! ------------------------------------------------------------------------------------

  function Leaky(variable,updatemodelprogress)

!   Numerical simulator for a radially-symmetric single layer of 
!   leaky aquifer.  Rock compressibility effects are included,
!   but not other temperature/ pressure effects on permeability or porosity.

!   Fixed parameter values for this model:
!   1: Layer elevation (not used)
!   2: Layer thickness
!   3: Action well radius
!   4: Rock specific heat
!   5: Rock heat conductivity
!   6: Rock density
!   7: Rock compressibility

!   Reservoir conditions for this model:
!   1: P (pressure)
!   2: X (temperature/ vapour saturation)

!   Variables for this model:
!   1: k (permeability)
!   2: phi (porosity)
!   3: RechargeCoef (leaky aquifer recharge coefficient) 

	use NumericalSimulator1D

!   Argument Variables:
    real(DP),intent(in) :: variable(:)
    real(DP) :: Leaky(TotalNData)
!   Fixed parameters:	                    
	real(DP)    :: LayerThickness,ActionWellRadius,CR,COND,RHOR,Compressibility
	real(DP)    :: AAA,Skin 
!   Reservoir conditions:
    real(DP)    :: InitialPressure,InitialX
!   Variable parameters:
    real(DP)    :: k,phi,RechargeCoef

    external updatemodelprogress

!   Unpack fixed parameters:
    LayerThickness=FixedParameter(2)
	ActionWellRadius=FixedParameter(3) 
    CR=FixedParameter(4)     
	COND=FixedParameter(5)   
	RHOR=FixedParameter(6)
	Compressibility=FixedParameter(7) 
	COMT=0.0_dp
	AAA=0.0_dp
	Skin=0.0_dp
!   Unpack reservoir conditions:
    InitialPressure=ReservoirCondition(1) 
	InitialX=ReservoirCondition(2) 
!   Unpack variable parameters:
    k=variable(1)   
    phi=variable(2) 
	RechargeCoef=variable(3)

	call SetupGrid(LayerThickness,ActionWellRadius) 
    call AssignBlockProperties(k,phi,Skin,ActionWellRadius,RechargeCoef,Compressibility)
	 
    Leaky=NumericalSolution1D(CR,COND,RHOR,AAA,&
	   InitialPressure,InitialX,updatemodelprogress,LayerThickness)

    return    
  end function Leaky 

! ------------------------------------------------------------------------------------

end module LeakySimulator