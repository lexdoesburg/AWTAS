module FractionalSimulator

  use variable_types
  use problem_data
  implicit none

  contains

! ------------------------------------------------------------------------------------

  function Fractional(variable,updatemodelprogress)

!   Numerical simulator for a radially-symmetric single layer of 
!   fractional-dimension medium.

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
!   3: fractional dimension

	use NumericalSimulator1D

!   Argument Variables:
    real(DP),intent(in) :: variable(:)
    real(DP) :: Fractional(TotalNData)
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
	COMT=0.0_dp
	AAA=0.0_dp 
	Skin=0.0_dp
	RechargeCoef=0.0_dp
!   Unpack reservoir conditions:
    InitialPressure=ReservoirCondition(1) 
	InitialX=ReservoirCondition(2) 
!   Unpack variable parameters:
    k=variable(1)   
    phi=variable(2)
	FractionalDimension=variable(3) 

	call SetupFractionalGrid(LayerThickness,ActionWellRadius,FractionalDimension)  
    call AssignBlockProperties(k,phi,Skin,ActionWellRadius,RechargeCoef,Compressibility)

    Fractional=NumericalSolution1D(CR,COND,RHOR,AAA,&
	    InitialPressure,InitialX,updatemodelprogress,LayerThickness)

    return    
  end function Fractional 

! ------------------------------------------------------------------------------------

end module FractionalSimulator