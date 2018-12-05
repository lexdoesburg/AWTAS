module NumericalSimulator1D

  use utility_functions
  use variable_types
  use problem_data
  use NumericalSimulator1D_Routines
  implicit none

! nb: M is the number of grid blocks.

  integer(I4B)              :: M
  real(DP), allocatable     :: DR(:),PER(:),POR(:),V(:)     
  real(DP), allocatable     :: P(:),T(:),SV(:),X(:)
  real(DP), allocatable     :: POLD(:),TOLD(:),SVOLD(:),XOLD(:)
  real(DP), allocatable     :: BMOLD(:),BEOLD(:)
  integer(I4B), allocatable :: IPH(:),IPHOLD(:)
  real(DP), allocatable     :: A(:),DELR(:),XK(:)
  real(DP), allocatable     :: XX(:),RR(:)
  real(DP), allocatable     :: COMP(:)
  real(DP)                  :: COMT
! Flowing enthalpy variables:
  real(DP)                  :: HF,HFOLD
! Recharge (leakage) coefficient:
  real(DP)                  :: XLAM

! This stores the grid block each obs. point is in:
  integer(I4B), allocatable :: ObsPointGridBlock(:)
! Weighting factors for extrapolating values at the sandface:
  real(DP):: SandFaceWeight1,SandFaceWeight2,SandFaceWeight3
! QQMM is the modelled value of flow rate:
  real(DP):: QQMM,ModFlowRate,ModFlowRateOld

  contains

! ------------------------------------------------------------------------------------

  function NumericalSolution1D(CR,COND,RHOR,AAA,&
    InitialPressure,InitialX,updatemodelprogress,LayerThickness)

!   Fast numerical simulator to calculate 1D radial solution.
!   Simulator originally written by Mike O'Sullivan 1999.
!   Temperature- & pressure-dependent permeability/porosity added (AC 11/99).
!   Skin effect included (AC 7/00).
!   Deliverability option added (AC 8/00).
!   Variable compressibility added (AC 4/02).

    use ModelProgress

!   Argument Variables:
    real(DP),intent(in) :: CR,COND,RHOR,AAA
	real(DP),intent(in) :: InitialPressure,InitialX,LayerThickness
    real(DP) :: NumericalSolution1D(TotalNData)
!   Local variables:
    real(DP)    :: P0(M),T0(M)
	integer(I4B):: I,IT,J
	real(DP)    :: dt,HIN,QMM,RMAX,time
	real(DP)    :: TestEndTime
	integer(I4B):: NTimeSteps, FlowIndex
	logical(LGT):: ResetTimeStepSize, StepFlows
	integer(I4B):: UpdateCounter
	real(DP)    :: PRECH,HRECH
	integer(I4B), parameter :: UpdateInterval=20 ! #iterations before updating progress
    real(DP), parameter     :: RTOL=1.E-6_dp  ! Newton-Raphson iteration tolerance
	integer(I4B), parameter :: IMAX=8         ! Max permitted no. of iterations
	real(DP), parameter     :: EPS=1.0E-6_dp  ! A smallish number

    external updatemodelprogress

	NumModelRuns=NumModelRuns+1

    TestEndTime=maxval(TestData(:)%time) ! Last observation time

    call AssignInitialConditions(InitialPressure,InitialX)
	call GetPumpingData(TestEndTime,StepFlows,.false.)

!   Initialise result:
    NumericalSolution1D=TestData%ModelledValue

    IGOOD=0     ! 'ok' parameter passed from thermo subroutines

!   Read in initial conditions and decide on the phase conditions:
    call INIT1(M,POLD,XOLD,TOLD,SVOLD,IPHOLD)
	P0=POLD
	T0=TOLD
	ModFlowRateOld=FLO(1)

!   Calculate recharge parameters (PRECH,TRECH):
    call GetLeakyRechargeParameters(POLD(1),TOLD(1),SVOLD(1),IPHOLD(1),PRECH,HRECH)

	time=0.0_dp
	ResetTimeStepSize=.true.
	FlowIndex=1
	NTimeSteps=0
	call UpdateTimeStepSize(time,TestEndTime,0,FlowIndex,ResetTimeStepSize,&
	  StepFlows,.false.,dt)
    UpdateCounter=0

!   Main time-stepping loop:
    do while (time<TestEndTime)

!     Calculate the accumulation terms BMOLD(I),BEOLD(I) at the old time:
	  call INIT2(M,POLD,XOLD,TOLD,SVOLD,IPHOLD,&
                 BMOLD,BEOLD,POR,RHOR,CR,P0,T0,COMP,COMT,AAA)

      time=time+dt
!     Calculate the mass flow rate and enthalpy (injection only) for the given time:
      call GetFlows(time,FlowIndex,QQMM,HIN,.true.,ResetTimeStepSize)
  
!     Start of Newton iterations:
10    IT=0

      if (UpdateCounter<UpdateInterval) then
	    UpdateCounter=UpdateCounter+1
	  else
!       Display progress & reset counter:
        call updatemodelprogress(NumModelRuns,dt,int(100.0*time/TestEndTime), &
		  analysisstopped)
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
      call SOLVE(P,T,SV,X,IPH,HF,M,POR,PER,A,V,CR,RHOR,COND,&
                 QQMM,HIN,XX,RR,DT,DELR,BMOLD,BEOLD,RMAX,&
                 P0,T0,COMP,COMT,AAA,QMM,XLAM,PRECH,HRECH,LayerThickness)

  	  if (IGOOD>0) then  ! Problems... reduce timestep:
	    call ReduceTimestep(time,dt,FlowIndex,QQMM,HIN,ResetTimeStepSize)
	    IGOOD=0
        GO TO 10 
      end if	  
	  	 
!     Add on the solution increments XX and make any necessary phase changes:
      call UPDATE(M,IPH,P,T,SV,X,XX,EPS)

!     Check for convergence:
      if (RMAX>=RTOL) then

        IT=IT+1
        if (IT.GT.IMAX) then
		  call ReduceTimestep(time,dt,FlowIndex,QQMM,HIN,ResetTimeStepSize)
          GO TO 10 
        end if 

        GO TO 11

      end if

	  ModFlowRate=QMM

!     Interpolate model response at observation point locations & times:
      call InterpolateResponse(time,dt,TestEndTime,NumericalSolution1D)

!     Update 'old' variables:
      POLD=P         
      TOLD=T         
      SVOLD=SV         
      XOLD=X         
      IPHOLD=IPH
	  HFOLD=HF
	  ModFlowRateOld=ModFlowRate

	  NTimeSteps=NTimeSteps+1

!     Update time step size:
      call UpdateTimeStepSize(time,TestEndTime,IT,FlowIndex,ResetTimeStepSize,&
	  StepFlows,.false.,dt)

    end do  !  End of main time-stepping loop.

    call updatemodelprogress(NumModelRuns,dt,100, analysisstopped)

!   Deallocate the main local arrays:
	call DestroyGridArrays
	call DestroyInterpolationArrays
	deallocate(DoneDataPoints)

    return
    
  end function NumericalSolution1D 

! ------------------------------------------------------------------------------------

  subroutine InterpolateResponse(NewTime,dt,TestEndTime,Response)
!   Interpolates model Response at the observation points and times,
!   from the results generated at the model block centres and modelled times.
!   Modelled results are calculated only for observations within the 
!   current model timestep [NewTime-dt,NewTime).
!   It's assumed that the data are in chronological order. 

!   Argument Variables:
    real(DP), intent(in)    :: NewTime, dt,TestEndTime
	real(DP), intent(inout) :: Response(TotalNData)
!   Local variables:
    real(DP)     :: OldTime
    integer(I4B) :: DataIndex,i,j
	integer(I4B) :: FirstDataPointNotDone
	real(DP)     :: ObsTime, theta, OneMinusTheta
	real(DP)     :: OldValue, NewValue
	integer(I4B) :: BlockNo
	real(DP)     :: SimulatedValue, DataOffset

    OldTime=NewTime-dt

!   Loop over observation points:
    do i=1,NObsPoints

      BlockNo=ObsPointGridBlock(i)
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
          case(0) ! Flows on deliverability:
		    if (ObsPoint(i)%IsPumpObsPoint) then
			  OldValue=ModFlowRateOld
			  NewValue=ModFlowRate
			end if 
		  case(1) ! Pressure:
		    OldValue=InterpolateBlockValue(POLD,BlockNo)
			NewValue=InterpolateBlockValue(P,BlockNo)
	      case(2) ! Temperature:
		    OldValue=InterpolateBlockValue(TOLD,BlockNo)
			NewValue=InterpolateBlockValue(T,BlockNo)
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

  function InterpolateBlockValue(Values,BlockNo)
!   Interpolates array Values at block BlockNo.  If BlockNo=0, extrapolates
!   quadratically to the sandface (ie. to the well radius).

!   Argument variables:
    real(DP), intent(in):: Values(M)
	integer(I4B), intent(in):: BlockNo
	real(DP):: InterpolateBlockValue

    if (BlockNo==0) then 
      InterpolateBlockValue=SandFaceWeight1*Values(1)+ &
	    SandFaceWeight2*Values(2)+SandFaceWeight3*Values(3)
	else
	  InterpolateBlockValue=Values(BlockNo)
	end if

    return
  end function InterpolateBlockValue

! ----------------------------------------------------------------------------------

  subroutine CalculateSandfaceWeights
!   Calculates weighting factors for sandface extrapolation.  
!   This is done by fitting a quadratic through the 
!   first 3 block centres to fit the value at the sandface.  

!   Locals:
    real(DP):: d1,x2,x3

    d1=DR(1)
	x2=0.5_dp*(d1+DR(2))
	x3=x2+0.5_dp*(DR(2)+DR(3))
    SandFaceWeight1=(4._dp*x2*x3+d1*(2._dp*(x2+x3)+d1))/(4._dp*x2*x3)
	SandFaceWeight2=d1*(2._dp*x3+d1)/(4._dp*x2*(x2-x3))
	SandFaceWeight3=d1*(2._dp*x2+d1)/(4._dp*x3*(x3-x2))
    
	return
  end subroutine CalculateSandfaceWeights

! ----------------------------------------------------------------------------------
!  Gridding routines
! ----------------------------------------------------------------------------------

  subroutine SetupGrid(LayerThickness,ActionWellRadius)
!   Sets up radial grid for the simulation, of thickness 'LayerThickness'.

!   Argument variables:
    real(DP), intent(in) :: LayerThickness, ActionWellRadius

	call SetupBlockSizes(ActionWellRadius) 
    call CalculateGridGeometryParameters(ActionWellRadius,LayerThickness)
!   Work out observation point positions on the grid:
    call GetObsPointGridPositions(ActionWellRadius)
!   Calculate weighting factors for sandface interpolation:
    if (not(WellBlockIncl)) call CalculateSandfaceWeights
	
	return
  end subroutine SetupGrid
  
! ------------------------------------------------------------------------

  subroutine SetupFractionalGrid(LayerThickness,ActionWellRadius,FractionalDimension)
!   Sets up fractional-dimension radial grid for the simulation, of thickness 
!   'LayerThickness' and dimension FractionalDimension.

!   Argument variables:
    real(DP), intent(in) :: LayerThickness, ActionWellRadius, FractionalDimension

	call SetupBlockSizes(ActionWellRadius) 
    call CalculateFractionalGridGeometryParameters(ActionWellRadius,&
	  LayerThickness,FractionalDimension)
!   Work out observation point positions on the grid:
    call GetObsPointGridPositions(ActionWellRadius)
!   Calculate weighting factors for sandface interpolation:
    call CalculateSandfaceWeights
	
	return
  end subroutine SetupFractionalGrid
    
! ------------------------------------------------------------------------

  subroutine SetupBlockSizes(ActionWellRadius)

!   The provisional grid used consists of:
!   - NConstBlocks blocks of constant radial thickness DR(i)=ConstDR
!   - blocks of exponentially increasing radial thickness (up to block M),
!     where DR(i)=GrowthFactor*DR(i-1).

    implicit none
	real(DP), intent(in)   :: ActionWellRadius
!   Local variables:
    real(DP),parameter     :: ConstDR=0.05_dp
	real(DP),parameter     :: GrowthFactor=1.2_dp
    integer(I4B),parameter :: NConstBlocks=10 
	integer(I4B)           :: i
!   No. of blocks (hard-coded for now):
	M=100

!   Allocate main local arrays:
	allocate(DR(M),PER(M),POR(M),V(M))
	allocate(P(M),T(M),SV(M),X(M))
	allocate(POLD(M),TOLD(M),SVOLD(M),XOLD(M))
	allocate(BMOLD(M),BEOLD(M))
	allocate(IPH(M),IPHOLD(M))
	allocate(A(M),DELR(M),XK(M),COMP(M))
	allocate(XX(2*M),RR(2*M)) ! Solving for 2 variables per block

!   Determine the block sizes, DR:
    if (WellBlockIncl) then
      DR(1)=ActionWellRadius
	  do i=2,NConstBlocks
	    DR(i)=ConstDR
	  end do     
	else ! No well block:     
	  do i=1,NConstBlocks
	    DR(i)=ConstDR
	  end do
	end if

	do i=1+NConstBlocks,M
	  DR(i)=GrowthFactor*DR(i-1)
	end do

    return
  end subroutine SetupBlockSizes

! ------------------------------------------------------------------------

  subroutine CalculateGridGeometryParameters(rw,b)

!   Calculates block interface areas, volumes and connection distances
!   for simple grid.  rw is the well radius; b is the layer thickness.

!   Argument variable:
    real(DP), intent(in) :: rw,b
!   Local variables:
    real(DP):: R0,R1
	integer(I4B):: i

    if (WellBlockIncl) then
	  R1=0._dp
	else
      R1=rw
	end if

    do i=1,M
      R0=R1
      R1=R0+DR(i)
      V(i)=PI_D*(R1+R0)*(R1-R0)*b
      if (i<M) then
        A(i)=2.0_dp*PI_D*R1*b
        DELR(i)=0.5_dp*(DR(i)+DR(i+1))
      end if
    end do

    return
  end subroutine CalculateGridGeometryParameters

! ------------------------------------------------------------------------

  subroutine CalculateFractionalGridGeometryParameters(rw,b,N)

!   Calculates block interface areas, volumes and connection distances
!   for fractional dimension grid.  rw is the well radius;
!   b is the layer thickness; N is the fractional dimension.

    use GammaFunction
!   Argument variables:
    real(DP), intent(in) :: rw,b,N
!   Local variables:
    real(DP):: R0,R1
	integer(I4B):: i
	real(DP):: HalfN,NMinus1,AlphaN
	real(DP):: V0,V1
	real(DP):: AreaCoefficient

    HalfN=0.5*N
	NMinus1=N-1.0_dp
    AlphaN=2.0_dp/DGAMMA(HalfN)*PI_D**HalfN
    AreaCoefficient=AlphaN*b**(3.0_dp-N)

    R1=rw
	V1=AreaCoefficient/N*R1**N

    do i=1,M
      R0=R1
	  V0=V1
      R1=R0+DR(i)
	  V1=AreaCoefficient/N*R1**N
      V(i)=V1-V0
      if (i<M) then
        A(i)=AreaCoefficient*R1**NMinus1
        DELR(i)=0.5*(DR(i)+DR(i+1))
      end if
    end do

    return
  end subroutine CalculateFractionalGridGeometryParameters

! ------------------------------------------------------------------------

  subroutine GetObsPointGridPositions(rw)
! Work out grid position of each observation point- this is used to 
! interpolate the model response from the grid values to the actual
! observation point position.  Here a piecewise-constant representation
! of pressure & temperatures over the grid is assumed, so we need only
! work out which block each observation point is in.

!   Argument:
    real(DP), intent(in):: rw
!   Local variables:
    integer(I4B):: i,block,StartBlock
    real(DP)    :: ObsPointR,R0,R1

    allocate(ObsPointGridBlock(NObsPoints))

    do i=1,NObsPoints

!     Radial position of observation point:
      ObsPointR=ObsPoint(i)%Position%x(1)

	  if ((rw>small_d).and.(ObsPointR<=rw)) then ! In well:
	    if (WellBlockIncl) then
		  ObsPointGridBlock(i)=1
		else
	      ObsPointGridBlock(i)=0 ! interpolate to sandface.
		end if
      else ! Search grid:
	    R1=rw
		if (WellBlockIncl) then
		  StartBlock=2
		else
		  StartBlock=1
		end if
	    do block=StartBlock,M
          R0=R1
		  R1=R0+DR(block)
		  if ((R0<=ObsPointR).and.(ObsPointR<R1)) then
		    ObsPointGridBlock(i)=block
		    exit
	      end if
        end do

	  end if

	end do

!   Also allocate and initialise the array holding the number of
!   datapoints that have been modelled so far in each observation point:
    allocate(DoneDataPoints(NObsPoints))
	DoneDataPoints=0

    return
  end subroutine GetObsPointGridPositions

! ------------------------------------------------------------------------

  subroutine AssignBlockProperties(k,phi,SkinFactor,ActionWellRadius,RechargeCoef,Compressibility)
!   Assigns permeability k and porosity phi to all blocks in the grid,
!   except for an altered-permeability 'skin zone' around the well.
!   The 'Safety' parameter keeps SkinFactor away from the asymptote in Ks
!   that occurs at SkinFactor=log(ActionWellRadius/Rs).

!   Argument Variables:
	real(DP), intent(in) :: k,phi,SkinFactor,ActionWellRadius,RechargeCoef,Compressibility
!   Local variables:
    real(DP):: Ks,TargetRs,Rs,R0,R1
	integer(I4B):: i,NumSkinBlocks
	real(DP), parameter :: Safety=0.1  

	TargetRs=max(2.0_dp,dexp(-SkinFactor+Safety))*ActionWellRadius

	! Line skin zone up with a block boundary:
	if (WellBlockIncl) then
	  R1=0._dp
	else
	  R1=ActionWellRadius
	end if
	Rs=R1
	do i=1,M
	  R1=R1+DR(i)
      if (R1>=TargetRs) then
	    Rs=R1
		NumSkinBlocks=i
	    exit
	  end if
	end do
	Ks=k/(1.0_dp+SkinFactor/dlog(Rs/ActionWellRadius))

!   Assign permeabilities:
	do i=1,NumSkinBlocks
	  PER(i)=Ks
	end do
	do i=NumSkinBlocks+1,M
	  PER(i)=k
	end do

	POR=phi

	XLAM=RechargeCoef
	COMP=Compressibility

	return
  end subroutine AssignBlockProperties

! ------------------------------------------------------------------------
      
  subroutine AssignInitialConditions(InitialPressure,InitialX)
!   Assigns initial conditions (InitialPressure,InitialX) to all 
!   blocks in the grid.

!   Argument Variables:
    real(DP),intent(in) :: InitialPressure,InitialX

	POLD=InitialPressure
	XOLD=InitialX
	HFOLD=0.0_dp

	return
  end subroutine AssignInitialConditions

! ------------------------------------------------------------------------

  subroutine GetLeakyRechargeParameters(PINIT,TINIT,SVINIT,IPHINIT,PRECH,HRECH)

!   Argument variables:
    real(DP),intent(in) :: PINIT,TINIT,SVINIT
	integer(I4B),intent(in):: IPHINIT
    real(DP),intent(out):: PRECH,HRECH
!   Locals:
    real(DP):: PP,TP,SVP,PORP,TE,TM
	real(DP):: RHOL,RHOV,HL,HV,VISL,VISV,XKRL,XKRV
	integer(I4B):: IPHP,IRELP

    PP=PINIT
    TP=TINIT
	SVP=SVINIT
    IPHP=IPHINIT
    call THERMO(PP,TP,SVP,IPHP,RHOL,RHOV,HL,HV,VISL,VISV)
    call RELP(SVINIT,IRELP,XKRL,XKRV)
    PRECH=PINIT
	select case (IPHP)
	case(0) ! liquid
	  HRECH=HL
	case(1) ! 2-phase
      TM=XKRL/VISL + XKRV/VISV
      TE=HL*XKRL/VISL + HV*XKRV/VISV
      HRECH=TE/TM
	case(2) ! steam
	  HRECH=HV
	end select

    return
  end subroutine GetLeakyRechargeParameters       

! ------------------------------------------------------------------------

  subroutine SetupWellboreParameters(WellVol,WellComp)
    implicit none
	real(DP), intent(in):: WellVol,WellComp

	V(1)=WellVol
    COMP(1)=WellComp

	return
  end subroutine SetupWellboreParameters

! ------------------------------------------------------------------------

  subroutine DestroyGridArrays
    deallocate(DR,PER,POR,V)
	deallocate(P,T,SV,X)
	deallocate(POLD,TOLD,SVOLD,XOLD)
	deallocate(BMOLD,BEOLD) 
	deallocate(IPH,IPHOLD) 
	deallocate(A,DELR,XK) 
	deallocate(XX,RR,COMP)
	deallocate(ObsPointGridBlock)
    return
  end subroutine DestroyGridArrays

! ------------------------------------------------------------------------------------

end module NumericalSimulator1D  


