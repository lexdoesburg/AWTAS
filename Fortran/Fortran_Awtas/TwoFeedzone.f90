module TwoFeedzoneSimulator

  use variable_types
  use problem_data
  use MultiFeedzoneSimulator
  implicit none

  contains

! ------------------------------------------------------------------------------------

  function TwoFeedzone(variable,updatemodelprogress)

!   Numerical simulator for a radially-symmetric two-feedzone model- 
!   a special case of the general multi-feedzone model.

!   Variable k/phi parameters have been omitted.  The wellbore permeability
!   and porosity are constant with depth. 

!   Fixed parameter values for this model:
!   1: Action well radius
!   2: Wellbore specific heat
!   3: Wellbore conductivity
!   4: Wellbore rock density
!   5: Wellbore compressibility
!   6: Wellhead elevation
!   7: Upper feedzone elevation
!   8: Lower feedzone elevation
!   9: Wellbore permeability 
!  10: Wellbore porosity
!  11: Upper feedzone thickness
!  12: Upper feedzone specific heat
!  13: Upper feedzone conductivity
!  14: Upper feedzone rock density
!  15: Upper feedzone compressibility
!  16: Lower feedzone thickness
!  17: Lower feedzone specific heat
!  18: Lower feedzone conductivity
!  19: Lower feedzone rock density
!  20: Lower feedzone compressibility

!   Reservoir conditions for this model:
!   1: Wellhead pressure
!   2...: (Depth(i),Temperature(i),i=1,NumTemperatures)

!   Variables for this model:
!   1: Upper feedzone permeability
!   2: Upper feedzone porosity
!   3: Lower feedzone permeability
!   4: Lower feedzone porosity

	use NumericalSimulatorMultiFeedzone
	use ModelProgress

    implicit none
!   Argument Variables:
    real(DP),intent(in) :: variable(:)
    real(DP) :: TwoFeedzone(TotalNData)
!   Fixed parameters:	                    
	real(DP)    :: ActionWellRadius
	real(DP)    :: CRWB,CONDWB,RHORWB,COMPWB,COMTWB,AAAWB
!   Locals:
    integer(I4B)::i,NumTopInactiveLayers,NumMiddleInactiveLayers
	real(DP):: DepthToUpperFeedzoneTop,DistBetweenFeedzones
	real(DP):: DepthToLowerFeedzoneTop,DepthToUpperFeedzoneBottom
	real(DP):: TopInactiveBlockSize,MiddleInactiveBlockSize
	real(DP), parameter::  TargetVerticalBlockSize=100.0_dp
	real(DP), parameter::  DefaultCR=1000.0_dp
	real(DP), parameter::  DefaultCond=2.5_dp
	real(DP), parameter::  DefaultCOMP=0.0_dp
	real(DP), parameter::  DefaultCOMT=0.0_dp
	real(DP), parameter::  DefaultAAA=0.0_dp
	real(DP), parameter::  DefaultRhoR=2500.0_dp
	real(DP), parameter::  DefaultPermeability=1.0E-18_dp
	real(DP), parameter::  DefaultPorosity=0.1_dp

    external updatemodelprogress

!   Default number of radial blocks:
    NR=99

!   Determine vertical block structure:
    DepthToUpperFeedzoneTop=FixedParameter(6)-(FixedParameter(7)+0.5_dp*FixedParameter(11))
	if (DepthToUpperFeedzoneTop>0.0_dp) then
      NumTopInactiveLayers=max(NINT(DepthToUpperFeedzoneTop/TargetVerticalBlockSize+0.5_dp),1)
      TopInactiveBlockSize=DepthToUpperFeedzoneTop/NumTopInactiveLayers
	else
	  NumTopInactiveLayers=0
	end if

    DepthToLowerFeedzoneTop=FixedParameter(6)-(FixedParameter(8)+0.5_dp*FixedParameter(16))
	DepthToUpperFeedzoneBottom=FixedParameter(6)-(FixedParameter(7)-0.5_dp*FixedParameter(11))
    DistBetweenFeedzones=DepthToLowerFeedzoneTop-DepthToUpperFeedzoneBottom
	if (DistBetweenFeedzones>0.0_dp) then
      NumMiddleInactiveLayers=max(NINT(DistBetweenFeedzones/TargetVerticalBlockSize+0.5_dp),1)
      MiddleInactiveBlockSize=DistBetweenFeedzones/NumMiddleInactiveLayers
	else
	  NumMiddleInactiveLayers=0
	end if

    MZ=NumTopInactiveLayers+NumMiddleInactiveLayers+2

    allocate(DZ(MZ),IACT(MZ))
	allocate(PERWB(MZ),PORWB(MZ))
    allocate(CRMA(MZ),CONDMA(MZ),RHORMA(MZ))
	allocate(COMPMA(MZ),COMTMA(MZ),AAAMA(MZ))
    allocate(PERMA(MZ,NR),PORMA(MZ,NR))

!   Unpack parameters:
	ActionWellRadius=FixedParameter(1) 
    CRWB=FixedParameter(2)     
	CONDWB=FixedParameter(3)   
	RHORWB=FixedParameter(4)   
	COMPWB=FixedParameter(5) 
	COMTWB=0.0_dp
	AAAWB=0.0_dp
	PERWB=FixedParameter(9)
	PORWB=FixedParameter(10)

    do i=1,NumTopInactiveLayers
      DZ(i)=TopInactiveBlockSize
	  IACT(i)=0
      CRMA(i)=DefaultCR ! Irrelevant values here...
	  CONDMA(i)=DefaultCond
	  RHORMA(i)=DefaultRhoR
	  COMPMA(i)=DefaultCOMP
	  COMTMA(i)=DefaultCOMT
	  AAAMA(i)=DefaultAAA
	  PERMA(i,:)=DefaultPermeability
	  PORMA(i,:)=DefaultPorosity
	end do

    i=NumTopInactiveLayers+1  ! Upper feedzone:
    DZ(i)=FixedParameter(11)
	IACT(i)=1
    CRMA(i)=FixedParameter(12) 
	CONDMA(i)=FixedParameter(13)
	RHORMA(i)=FixedParameter(14)
	COMPMA(i)=FixedParameter(15)
	COMTMA(i)=DefaultCOMT
	AAAMA(i)=DefaultAAA
	PERMA(i,:)=variable(1)
	PORMA(i,:)=variable(2)

	do i=NumTopInactiveLayers+2,MZ-1 ! Middle inactive layers:
      DZ(i)=MiddleInactiveBlockSize
	  IACT(i)=0
      CRMA(i)=DefaultCR ! Irrelevant values here...
	  CONDMA(i)=DefaultCond
	  RHORMA(i)=DefaultRhoR
	  COMPMA(i)=DefaultCOMP
	  COMTMA(i)=DefaultCOMT
	  AAAMA(i)=DefaultAAA
	  PERMA(i,:)=DefaultPermeability
	  PORMA(i,:)=DefaultPorosity
	end do

	i=MZ  ! Lower feedzone:
    DZ(i)=FixedParameter(16)
	IACT(i)=1
    CRMA(i)=FixedParameter(17) 
	CONDMA(i)=FixedParameter(18)
	RHORMA(i)=FixedParameter(19)
	COMPMA(i)=FixedParameter(20)
	COMTMA(i)=DefaultCOMT
	AAAMA(i)=DefaultAAA
	PERMA(i,:)=variable(3)
	PORMA(i,:)=variable(4)

	call SetupGrid(ActionWellRadius)

	call CalculateInitialConditions(CRWB,CONDWB,RHORWB,COMPWB,&
	  COMTWB,AAAWB,updatemodelprogress)

    TwoFeedzone=NumericalSolutionMultiFeedzone(CRWB,CONDWB,RHORWB,COMPWB,&
	  COMTWB,AAAWB,.false.,updatemodelprogress)

	call DestroyGridArrays
    call DestroyParameterArrays
    deallocate(DoneDataPoints)

    return    
  end function TwoFeedzone 

! ------------------------------------------------------------------------

end module TwoFeedzoneSimulator