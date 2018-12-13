module FractionalVariableKPhiSimulator

  use variable_types
  use problem_data
  implicit none

  contains

! ------------------------------------------------------------------------------------

  function FractionalVariableKPhi(variable,updatemodelprogress)

!   Numerical simulator for a radially-symmetric single layer of 
!   fractional-dimension medium, with permeability & porosity dependent
!   on temperature and pressure.

!   Fixed parameter values for this model:
!   1: Layer elevation (not used)
!   2: Layer thickness
!   3: Action well radius
!   4: Rock specific heat
!   5: Rock heat conductivity
!   6: Rock density
!   7: Rock compressibility
!   8: Rock expansivity

!   Reservoir conditions for this model:
!   1: P (pressure)
!   2: X (temperature/ vapour saturation)

!   Variables for this model:
!   1: k (permeability)
!   2: phi (porosity)
!   3: Rock k/phi parameter (sensitivity of k to phi) 
!   4: fractional dimension

	use NumericalSimulator1D

!   Argument Variables:
    real(DP),intent(in) :: variable(:)
    real(DP) :: FractionalVariableKPhi(TotalNData)
!   Fixed parameters:	                    
	real(DP)    :: LayerThickness,ActionWellRadius,CR,COND,RHOR,Compressibility
	real(DP)    :: AAA,Skin,RechargeCoef
!   Reservoir conditions:
    real(DP)    :: InitialPressure,InitialX
!   Variable parameters:
    real(DP)    :: k,phi,FractionalDimension

    external updatemodelprogress

!   Unpack fixed parameters:
    LayerThickness=FixedParameter(2)
	ActionWellRadius=FixedParameter(3) 
    CR=FixedParameter(4)     
	COND=FixedParameter(5)   
	RHOR=FixedParameter(6) 
	Compressibility=FixedParameter(7)  
	COMT=FixedParameter(8)
	Skin=0.0_dp
	RechargeCoef=0.0_dp
!   Unpack reservoir conditions:
    InitialPressure=ReservoirCondition(1) 
	InitialX=ReservoirCondition(2) 
!   Unpack variable parameters:
    k=variable(1)   
    phi=variable(2)
	AAA=variable(3) 
	FractionalDimension=variable(4)

	call SetupFractionalGrid(LayerThickness,ActionWellRadius,FractionalDimension)  
    call AssignBlockProperties(k,phi,Skin,ActionWellRadius,RechargeCoef,Compressibility)

    FractionalVariableKPhi=NumericalSolution1D(CR,COND,RHOR,AAA,&
	    InitialPressure,InitialX,updatemodelprogress,LayerThickness)

    return    
  end function FractionalVariableKPhi 

! ------------------------------------------------------------------------------------

end module FractionalVariableKPhiSimulator