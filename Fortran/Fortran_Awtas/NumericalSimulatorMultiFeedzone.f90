module NumericalSimulatorMultiFeedzone

  use utility_functions
  use variable_types
  use problem_data
  use NumericalSimulatorMultiFeedzone_Routines
  implicit none

! nb: NR is the number of grid blocks in the radial direction,
! MZ is the number in the vertical.

  integer(I4B)              :: MZ,NR
  real(DP), allocatable     :: DZ(:)
  real(DP), allocatable     :: PRAT(:)
  real(DP), allocatable     :: V(:,:)       
  real(DP), allocatable     :: P(:,:),T(:,:),SV(:,:),X(:,:)
  real(DP), allocatable     :: POLD(:,:),TOLD(:,:),SVOLD(:,:)
  real(DP), allocatable     :: XOLD(:,:)
  real(DP), allocatable     :: BMOLD(:,:),BEOLD(:,:)
  integer(I4B), allocatable :: IPH(:,:),IPHOLD(:,:)
  real(DP), allocatable     :: DR(:),AR(:,:),DELR(:)
  real(DP), allocatable     :: AZ(:),DELZ(:)
  real(DP), allocatable     :: XX(:,:),RR(:)
! Parameter arrays:
  integer(I4B), allocatable:: IACT(:)
  real(DP), allocatable:: PERWB(:),PORWB(:)
  real(DP), allocatable:: CRMA(:),CONDMA(:),RHORMA(:)
  real(DP), allocatable:: COMPMA(:),COMTMA(:),AAAMA(:)
  real(DP),allocatable:: PERMA(:,:),PORMA(:,:)

! This stores the grid block each obs. point is in:
  integer(I4B), allocatable :: ObsPointGridBlock(:,:)
! Flowing enthalpy variables:
  real(DP)                  :: HF, HFOLD

  contains

! ------------------------------------------------------------------------------------

  function NumericalSolutionMultiFeedzone(CRWB,CONDWB,RHORWB,COMPWB,&
	  COMTWB,AAAWB,RunToSS,updatemodelprogress)

    use ModelProgress
	implicit none

!   Argument Variables:
    real(DP),intent(in) :: CRWB,CONDWB,RHORWB,COMPWB,COMTWB,AAAWB
    real(DP) :: NumericalSolutionMultiFeedzone(TotalNData)
	logical(LGT),intent(in):: RunToSS
!   Local variables:
    real(DP)    :: P0(MZ,NR),T0(MZ,NR)
	integer(I4B):: I,IT,J
	real(DP)    :: dt,HIN,QMM,RMAX,RMAXMM,RMAXME,RMAXFM,RMAXFE
	real(DP)    :: time,TestEndTime
	integer(I4B):: NTimeSteps, FlowIndex
	logical(LGT):: ResetTimeStepSize, StepFlows, Finished
	integer(I4B):: UpdateCounter
	integer(I4B), parameter :: UpdateInterval=20 ! #iterations before updating progress
    real(DP), parameter     :: RTOL=1.E-7_dp  ! Newton-Raphson iteration tolerance
	integer(I4B), parameter :: IMAX=10        ! Max permitted no. of iterations
	real(DP), parameter     :: EPS=1.0E-6_dp  ! A smallish number
    real(DP),parameter      :: BigTime=1.0E15_dp ! Max. end time for SS runs
	integer(I4B), parameter :: OutFile=15
    external updatemodelprogress

    open(unit=OutFile,name='AWTASLog_MF.dat',status='replace',dispose='delete')

	DoneDataPoints=0

    if (RunToSS) then
	  TestEndTime=BigTime
	else
      TestEndTime=maxval(TestData(:)%time) ! Last observation time
	  NumModelRuns=NumModelRuns+1
	end if

    call AssignBlockProperties
	call GetPumpingData(TestEndTime,StepFlows,RunToSS)

!   Initialise result:
    NumericalSolutionMultiFeedzone=TestData%ModelledValue

    IGOOD=0     ! 'ok' parameter passed from thermo subroutines

    P0=POLD
	T0=TOLD

	time=0.0_dp
	ResetTimeStepSize=.true.
	FlowIndex=1
	NTimeSteps=0
	call UpdateTimeStepSize(time,TestEndTime,0,FlowIndex,ResetTimeStepSize,&
	  StepFlows,RunToSS,dt)
    UpdateCounter=0
    Finished=.false.

!   Main time-stepping loop:
    do while (not(Finished))

!     Calculate the accumulation terms BMOLD(I),BEOLD(I) at the old time:
	  call INIT2(NR,MZ,POLD,XOLD,TOLD,SVOLD,IPHOLD,&
               BMOLD,BEOLD,PORWB,RHORWB,CRWB,PORMA,RHORMA,CRMA,&
               P0,T0,COMPWB,COMTWB,AAAWB,COMPMA,COMTMA,AAAMA)

      time=time+dt
!     Calculate the mass flow rate and enthalpy (injection only) for the given time:
      call GetFlows(time,FlowIndex,QMM,HIN,.true.,ResetTimeStepSize)
  
!     Start of Newton iterations:
10    IT=0

      if (UpdateCounter<UpdateInterval) then
	    UpdateCounter=UpdateCounter+1
	  else
!       Display progress & reset counter:
        call updatemodelprogress(NumModelRuns,dt,int(100.*time/TestEndTime), analysisstopped)
		UpdateCounter=0
	  end if

!     Terminate the model run if analysis has been stopped:
	  if (analysisstopped) exit

!     Initialise:
      P=POLD         
      T=TOLD         
      SV=SVOLD         
      X=XOLD         
      IPH=IPHOLD 
	  HF=HFOLD

!     Notation: P,T, ... are new time values, POLD, TOLD, ... are old time values
!     at the start of each time step, or after a Newton-Raphson failure, new values are set 
!     equal to old values.

11    CONTINUE

!     Carry out most of one step of the Newton-Raphson process:
         CALL SOLVEMF(P,T,SV,X,IPH,XX,BMOLD,BEOLD,NR,MZ,PORWB,PORMA,&
             PERWB,PERMA,AR,AZ,V,CRWB,CRMA,RHORWB,RHORMA,CONDWB,CONDMA,&
             QMM,HIN,DT,DELR,DELZ,RMAXMM,RMAXME,RMAXFM,RMAXFE,&
             P0,T0,COMPWB,COMTWB,AAAWB,COMPMA,COMTMA,AAAMA,PRAT,IACT)

!      write (Outfile,'(i4,1x,i2,1x,i2,1x,E10.4,1x,E10.4)') NTimeSteps,IT,IGOOD,time,dt

  	  if (IGOOD>0) then  ! Problems... reduce timestep:
	    call ReduceTimestep(time,dt,FlowIndex,QMM,HIN,ResetTimeStepSize)
	    IGOOD=0
        GO TO 10 
      end if	  
	  	 
!     Add on the solution increments XX and make any necessary phase changes:
      call UPDATE(NR,MZ,IPH,P,T,SV,X,XX,EPS)

!     Check for convergence:
      RMAX=max(RMAXMM,RMAXME,RMAXFM,RMAXFE)
      if (RMAX>=RTOL) then

        IT=IT+1
        if (IT>IMAX) then
		  call ReduceTimestep(time,dt,FlowIndex,QMM,HIN,ResetTimeStepSize)
          GO TO 10 
        end if 

        GO TO 11

      end if

!     Interpolate model response at observation point locations & times:
      call InterpolateResponse(time,dt,TestEndTime,NumericalSolutionMultiFeedzone)

!     Update time step size:
      call UpdateTimeStepSize(time,TestEndTime,IT,FlowIndex,ResetTimeStepSize,&
	  StepFlows,RunToSS,dt)

      if (RunToSS) then
	    Finished=SSConverged(dt)
	  else
	    Finished=(time>=TestEndTime)
	  end if

!     Update 'old' variables:
      POLD=P         
      TOLD=T         
      SVOLD=SV         
      XOLD=X         
      IPHOLD=IPH
	  HFOLD=HF

	  NTimeSteps=NTimeSteps+1

    end do  !  End of main time-stepping loop.

    call updatemodelprogress(NumModelRuns,dt,100, analysisstopped)

    call DestroyInterpolationArrays

	close(OutFile)

    return
  end function NumericalSolutionMultiFeedzone  

! ------------------------------------------------------------------------------------

  subroutine InterpolateResponse(NewTime,dt,TestEndTime,Response)
!   Interpolates response in space & time.  Pressure, temperature and enthalpy
!   values are taken from the fracture layer (layer 1).

!   Argument Variables:
    real(DP), intent(in)    :: NewTime, dt,TestEndTime
	real(DP), intent(inout) :: Response(TotalNData)

!   Local variables:
    real(DP)     :: OldTime
    integer(I4B) :: DataIndex,i,j
	integer(I4B) :: FirstDataPointNotDone
	real(DP)     :: ObsTime, theta, OneMinusTheta
	real(DP)     :: OldValue, NewValue
	integer(I4B) :: Blocki,Blockj
	real(DP)     :: SimulatedValue, DataOffset

    OldTime=NewTime-dt

!   Loop over observation points:
    do i=1,NObsPoints

      Blocki=ObsPointGridBlock(i,2) ! Vertical position
      Blockj=ObsPointGridBlock(i,1) ! Horizontal position
	  DataOffset=ObsPoint(i)%DataOffset

!     Loop over remaining datapoints:
      FirstDataPointNotDone=DoneDataPoints(i)+1

      do j=FirstDataPointNotDone,ObsPoint(i)%NData

!       Index of datapoint in the TestData array:
	    DataIndex=ObsPoint(i)%DataIndex+j-1

		ObsTime=TestData(DataIndex)%time

!       Scaled position of ObsTime within the time step:
		theta=(ObsTime-OldTime)/dt

!       If theta<0, reset to 0:
        theta=max(theta,0.0_dp)  
!       (This shouldn't really happen, but could if there 
!       are datapoints for t<0, for example.)

	    if ((theta<1.0_dp).or.(NewTime==TestEndTime)) then

!         Interpolate modelled value.  Since the underlying model is
!         a finite-volume model, spatial interpolation is assumed 
!         constant over each grid block.

		  OneMinusTheta=1.0_dp-theta

          select case (ObsPoint(i)%Property)
!         case(0) (flow) isn't simulated by this model.
		  case(1) ! Pressure:
		    OldValue=InterpolateBlockValue(POLD,Blocki,Blockj)
			NewValue=InterpolateBlockValue(P,Blocki,Blockj)
	      case(2) ! Temperature:
		    OldValue=InterpolateBlockValue(TOLD,Blocki,Blockj)
			NewValue=InterpolateBlockValue(T,Blocki,Blockj)
		  case(3) ! Enthalpy:
		    OldValue=HFOLD 
			NewValue=HF
		  end select

		  SimulatedValue=OneMinusTheta*OldValue+theta*NewValue
		  Response(DataIndex)=SimulatedValue-DataOffset

          DoneDataPoints(i)=DoneDataPoints(i)+1

        else ! theta>=1:

!         Finished for this observation point.
		  exit		  

		end if

	  end do ! End datapoint loop

	end do ! End observation point loop

	return
  end subroutine InterpolateResponse

!-----------------------------------------------------------------------------------

  logical(LGT) function SSConverged(dt)
!   Determines whether steady-state run has converged.
!   Argument:
    real(DP),intent(in):: dt
!   Locals:
    real(DP):: MaxPChange,MaxTChange
	integer(I4B):: i,j
	real(DP),parameter:: PTol=1.0E-3_dp  ! 0.01% tolerances
	real(DP),parameter:: TTol=1.0E-3_dp
    real(DP),parameter:: MinConvergedDT=1.E6_dp

    MaxPChange=0.0_dp
	MaxTChange=0.0_dp
    do i=1,MZ
	  do j=1,NR
        MaxPChange=max(MaxPChange,(P(i,j)-POLD(i,j))/POLD(i,j))
        MaxTChange=max(MaxTChange,(T(i,j)-TOLD(i,j))/TOLD(i,j))
      end do
	end do

    SSConverged=((MaxPChange<=PTol).and.(MaxTChange<=TTol).and.(dt>=MinConvergedDT))
    
    return
  end function SSConverged

!-----------------------------------------------------------------------------------

  function InterpolateBlockValue(Values,Blocki,Blockj)
!   Interpolates array Values at block (Blocki,Blockj)

!   Argument variables:
    real(DP), intent(in):: Values(MZ,NR)
	integer(I4B), intent(in):: Blocki,Blockj
	real(DP):: InterpolateBlockValue

	InterpolateBlockValue=Values(Blocki,Blockj)

    return
  end function InterpolateBlockValue

! ----------------------------------------------------------------------------------
!  Gridding routines
! ----------------------------------------------------------------------------------

  subroutine SetupGrid(ActionWellRadius)
!   Sets up grid for the simulation.

!   Argument variables:
    real(DP), intent(in) :: ActionWellRadius

	call SetupBlockSizes(ActionWellRadius)
    call CalculateGridGeometryParameters
!   Work out observation point positions on the grid:
    call GetObsPointGridPositions
	
	return
  end subroutine SetupGrid
  
! ------------------------------------------------------------------------

  subroutine SetupBlockSizes(ActionWellRadius)
!   Argument variables:
    real(DP), intent(in)   :: ActionWellRadius
!   Local variables:
    real(DP),parameter     :: ConstDR=0.10_dp
	real(DP),parameter     :: GrowthFactor=1.2_dp
    integer(I4B),parameter :: NConstBlocks=10
	integer(I4B)           :: i
	real(DP)               :: alpha

!   Allocate main local arrays:
    allocate (PRAT(MZ))
    allocate (V(MZ,NR),P(MZ,NR),T(MZ,NR),SV(MZ,NR),X(MZ,NR))
    allocate (POLD(MZ,NR),TOLD(MZ,NR),SVOLD(MZ,NR))
    allocate (XOLD(MZ,NR))
    allocate (BMOLD(MZ,NR),BEOLD(MZ,NR))
    allocate (IPH(MZ,NR),IPHOLD(MZ,NR))
    allocate (DR(NR),AR(MZ,NR),DELR(NR))
    allocate (AZ(NR),DELZ(MZ-1))
    allocate (XX(2*MZ,NR),RR(2*NR))

!   Determine the radial block sizes:
    DR(1)=ActionWellRadius  ! Well block
	do i=2,NConstBlocks+1
	  DR(i)=ConstDR
	end do
	do i=NConstBlocks+2,NR
	  DR(i)=GrowthFactor*DR(i-1)
	end do

    return
  end subroutine SetupBlockSizes

! ------------------------------------------------------------------------

  subroutine CalculateGridGeometryParameters
!   Calculates block volumes and interface areas.
!   Local variables:
    real(DP):: R0,R1,thick
	integer(I4B):: i,j

    R1=0.0_dp
    do j=1,NR
      R0=R1
      R1=R0+DR(j)
      do i=1,MZ
        thick=DZ(i)
        V(i,j)=PI_D*(R0+R1)*(R1-R0)*thick
        if (i<MZ) then
          DELZ(i)=0.5_dp*(DZ(i)+DZ(i+1))
        end if
        if (j<NR) then
          AR(i,j)=2.0_dp*PI_D*R1*DZ(i)
        end if       
      end do
      AZ(j)=PI_D*(R0+R1)*(R1-R0)
      if (j<NR) then
        DELR(j)=0.5_dp*(DR(j)+DR(j+1))
      end if       
    end do

    return
  end subroutine CalculateGridGeometryParameters

! ------------------------------------------------------------------------

  subroutine GetObsPointGridPositions

!   Local variables:
    integer(I4B):: i,block
    real(DP)    :: ObsPointR,ObsPointZ,R0,R1,Z0,Z1
	real(DP)    :: WellheadElev

    allocate(ObsPointGridBlock(NObsPoints,2))

	WellheadElev=FixedParameter(6)

    do i=1,NObsPoints

!     Radial position of observation point:
      ObsPointR=ObsPoint(i)%Position%x(1)

	  R1=0.0_dp
	  do block=1,NR
        R0=R1
		R1=R0+DR(block)
		if ((R0<=ObsPointR).and.(ObsPointR<R1)) then
		  ObsPointGridBlock(i,1)=block
		  exit
	    end if
      end do

!     Vertical position of observation point:
      ObsPointZ=ObsPoint(i)%Position%x(3)
	  if (ObsPointZ>WellheadElev) then
	    ObsPointGridBlock(i,2)=1 ! Tool above wellhead
	  else
	    Z1=WellheadElev
	    do block=1,MZ
          Z0=Z1
		  Z1=Z0-DZ(block)
		  if ((Z0>=ObsPointZ).and.(ObsPointZ>Z1)) then
		    ObsPointGridBlock(i,2)=block
		    exit
	      end if
	    end do
		if (ObsPointZ<Z1) then ! Tool below lower feedzone:
		  ObsPointGridBlock(i,2)=MZ
		end if
      end if

	end do

!   Also allocate and initialise the array holding the number of
!   datapoints that have been modelled so far in each observation point:
    allocate(DoneDataPoints(NObsPoints))

    return
  end subroutine GetObsPointGridPositions

! ------------------------------------------------------------------------

  subroutine AssignBlockProperties
!   Local:
    integer(I4B):: i

    do i=1,MZ
	  PRAT(i)=PERWB(i)/PERMA(i,1)
    end do

	return
  end subroutine AssignBlockProperties

! ------------------------------------------------------------------------
      
  subroutine CalculateInitialConditions(CRWB,CONDWB,RHORWB,COMPWB,&
	  COMTWB,AAAWB,updatemodelprogress)

!   Takes specified initial temperature profile in the well and calculates
!   approximate pressure profile.  This is done by making all layers inactive
!   (so only the well is simulated) and running the numerical simulator until 
!   an approximate steady state is reached.
 
    implicit none 
!   Argument Variables:
    real(DP),intent(inout) :: CRWB
	real(DP),intent(in):: CONDWB,RHORWB,COMPWB,COMTWB,AAAWB
!   Locals:
    integer(I4B):: i,j
	integer(I4B):: IACTValue(MZ)
	real(DP):: soln(TotalNData)
	real(DP),parameter:: BigCR=1.0E3_dp

    external updatemodelprogress

    call InterpolateTemperatures
	call AssignHydrostaticPressures
	call AssignInitialX

    call INIT1(NR,MZ,POLD,XOLD,TOLD,SVOLD,IPHOLD)
	HFOLD=0.0_dp

!   Make all layers inactive:
    IACTValue=IACT ! (save original IACT array)
	IACT=0

!   Run to steady state to recalculate pressures:
	soln=NumericalSolutionMultiFeedzone(CRWB,CONDWB,RHORWB,COMPWB,&
	  COMTWB,AAAWB,.true.,updatemodelprogress)

	POLD=P
	XOLD=X
	TOLD=T
	SVOLD=SV
	HFOLD=HF
	IPH=IPHOLD

!   Restore original IACT array:
	IACT=IACTValue

	return
  end subroutine CalculateInitialConditions

! ------------------------------------------------------------------------

  subroutine AssignHydrostaticPressures
!   Assigns hydrostatic pressures to the POLD array.

    use Thermodynamics

    implicit none
!   Locals:
    integer(I4B):: i
	real(DP):: Ps,Ts,rho,U
	real(DP), parameter:: g=9.8_dp

!   Extrapolate surface temperature from top two block values:
    Ts=TOLD(1,1)-0.5_dp*DZ(1)*(TOLD(2,1)-TOLD(1,1))/DELZ(1)
    Ps=ReservoirCondition(1)

    IGOOD=0
!   Get density rho at surface:
	call COWAT(Ts,Ps,rho,U)

!   Pressure in first layer:
 	POLD(1,:)=Ps+rho*g*0.5_dp*DZ(1)

    do i=1,MZ-1
	  call COWAT(TOLD(i,1),POLD(i,1),rho,U)
	  POLD(i+1,:)=POLD(i,1)+rho*g*DELZ(i)
	end do

    return
  end subroutine AssignHydrostaticPressures

! ------------------------------------------------------------------------

  subroutine InterpolateTemperatures

    implicit none
!   Interpolates downhole temperature profile to give initial temperature
!   conditions at each layer in the wellbore.

!   It's assumed that the temperature data are in depth order.

!   Locals:
    real(DP),allocatable:: Depth(:),Temp(:)
	integer(I4B):: NTempData,i,j
    real(DP):: LayerDepth(MZ),d0,d1,theta,z,t0,t1
	logical(LGT):: OutOfData

	NTempData=(NReservoirConditions-1)/2  ! Should be an even number

	if (NTempData>1) then
	  allocate(Depth(NTempData),Temp(NTempData))
      do i=1,NTempData
	    Depth(i)=ReservoirCondition(2*i)
	    Temp(i)=ReservoirCondition(2*i+1)
	  end do

      LayerDepth(1)=0.5_dp*DZ(1)
	  do i=2,MZ
	    LayerDepth(i)=LayerDepth(i-1)+DELZ(i-1)
	  end do

      j=1
	  d0=Depth(1)
	  d1=Depth(2)
	  t0=Temp(1)
	  t1=Temp(2)
	  OutOfData=.false.

      do i=1,MZ
	    if (OutOfData) then
	      theta=1.0_dp
	    else
          z=LayerDepth(i)
	      theta=(z-d0)/(d1-d0)
	      if (theta<0.0_dp) then
	        theta=0.0_dp
	      else
	        do while ((theta>1.0_dp).and.(not(OutOfData))) 
		      j=j+1
		      if (j>=NTempData) then
		        OutOfData=.true.
			    theta=1.0_dp
		      else
		        d0=d1
			    d1=Depth(j+1)
			    theta=(z-d0)/(d1-d0)
			    t0=t1
			    t1=Temp(j+1)
		      end if
		    end do
	      end if
	    end if
        TOLD(i,:)=(1._dp-theta)*t0+theta*t1
      end do
      deallocate(Depth,Temp)

	else ! Only one temperature data point:
	  TOLD(:,:)=ReservoirCondition(3)
	end if

    return
  end subroutine InterpolateTemperatures

! ------------------------------------------------------------------------

  subroutine AssignInitialX
!   Calculates approximate starting values for the X variable.
    implicit none
!   Locals:
    integer(I4B):: i
	real(DP):: PSat
	real(DP),parameter:: DefaultTwoPhaseSat=0.1_dp

	do i=1,MZ
	  call SAT(TOLD(i,1),PSat)
	  if (POLD(i,1)<=PSat) then ! 2-phase
	    POLD(i,:)=PSat
		XOLD(i,:)=DefaultTwoPhaseSat
	  else ! 1-phase
	    XOLD(i,:)=TOLD(i,1)
	  end if
	end do

    return
  end subroutine AssignInitialX

! ------------------------------------------------------------------------

  subroutine DestroyGridArrays
    implicit none
    deallocate (DZ,PRAT)
    deallocate (V,P,T,SV,X)
    deallocate (POLD,TOLD,SVOLD,XOLD)
    deallocate (BMOLD,BEOLD,IPH,IPHOLD)
    deallocate (DR,AR,DELR)
    deallocate (AZ,DELZ)
    deallocate (XX,RR)
	deallocate(ObsPointGridBlock)
    return
  end subroutine DestroyGridArrays

! ------------------------------------------------------------------------------------

  subroutine DestroyParameterArrays
    implicit none
    deallocate(IACT)
	deallocate(PERWB,PORWB)
	deallocate(CRMA,CONDMA,RHORMA)
	deallocate(COMPMA,COMTMA,AAAMA)
	deallocate(PERMA,PORMA)
    return
  end subroutine DestroyParameterArrays

! ------------------------------------------------------------------------

end module NumericalSimulatorMultiFeedzone