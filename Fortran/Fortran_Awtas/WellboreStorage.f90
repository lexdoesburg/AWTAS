module WellboreStorageSimulator

  use variable_types
  use problem_data
  implicit none

  contains

! ------------------------------------------------------------------------------------

  function WellboreStorage(variable,updatemodelprogress)

!   Numerical simulator for homogeneous porous model with simple wellbore storage model,
!   containing a single large well block with high compressibility and volume representing
!   the wellbore volume.

!   Fixed parameter values for this model:
!   1: Layer elevation (not used)
!   2: Layer thickness
!   3: Action well radius
!   4: Rock specific heat
!   5: Rock heat conductivity
!   6: Rock density
!   7: Rock compressibility
!   8: Wellbore volume

!   Reservoir conditions for this model:
!   1: P (pressure)
!   2: X (temperature/ vapour saturation)

!   Variables for this model:
!   1: k (permeability)
!   2: phi (porosity)
!   3: wellcomp (wellbore compressibility)

	use NumericalSimulator1D

!   Argument Variables:
    real(DP),intent(in) :: variable(:)
    real(DP) :: WellboreStorage(TotalNData)
!   Fixed parameters:	                    
	real(DP)    :: LayerThickness,ActionWellRadius,CR,COND,RHOR,Compressibility
	real(DP)    :: AAA,Skin,RechargeCoef,WellVol
!   Reservoir conditions:
    real(DP)    :: InitialPressure,InitialX
!   Variable parameters:
    real(DP)    :: k,phi,wellcomp

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
	WellVol=FixedParameter(8)
!   Unpack reservoir conditions:
    InitialPressure=ReservoirCondition(1) 
	InitialX=ReservoirCondition(2) 
!   Unpack variable parameters:
    k=variable(1)   
    phi=variable(2)
	wellcomp=variable(3)

	call SetupGrid(LayerThickness,ActionWellRadius)
	 
    call AssignBlockProperties(k,phi,Skin,ActionWellRadius,RechargeCoef,Compressibility)

	call SetupWellboreParameters(WellVol,wellcomp)
	 
    WellboreStorage=NumericalSolution1D(CR,COND,RHOR,AAA,&
	   InitialPressure,InitialX,updatemodelprogress,LayerThickness)

    return    
  end function WellboreStorage 

! ------------------------------------------------------------------------------------

end module WellboreStorageSimulator