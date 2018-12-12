module problem_data

	! This module stores all the model-independent problem data, together with 
	! routines for operating on it.
	
		use variable_types
		implicit none
	
		type ThreeVector
			real(DP) x(3)
		end type ThreeVector
	
		type pump_type
			integer(I4B) :: Scheme, NData
		real(DP) :: Enthalpy
		logical(LGT) :: StepFlows
		logical(LGT) :: OnDeliverability
		real(DP)     :: ProdIndex,CutoffPressure
		end type pump_type
	
		type pump_data_type
			real(DP) :: time, rate
		end type pump_data_type
		
		type ObsPoint_type
			integer(I4B) :: Property, NData
			real(DP)     :: Error, Weight
		real(DP)     :: DataOffset
			integer(I4B) :: DataIndex
			type(ThreeVector) :: Position
		logical(LGT) :: IsPumpObsPoint  ! For pumps on deliverability
		integer(I4B) :: PumpNo          ! Index of corresponding pump on deliv.
		end type ObsPoint_type
	
		type datapoint
			 real(DP) :: time,value,ModelledValue,error,weight
		end type datapoint
		
		integer(I4B) :: NPumps, NObsPoints 
		integer(I4B) :: MaxNPumpTimes, TotalNData, MaxNPumpSchemeParams
	
		type(pump_type),allocatable      :: Pump(:)
		type(ObsPoint_type),allocatable  :: ObsPoint(:)
		type(datapoint),allocatable      :: TestData(:)
		type(pump_data_type),allocatable :: PumpData(:,:)
		real(DP),allocatable             :: PumpSchemeParams(:,:)
		
		integer(I4B) :: NReservoirConditions
		integer(I4B) :: NFixedParameters
		real(DP),allocatable :: ReservoirCondition(:)
		real(DP),allocatable :: FixedParameter(:) 
	
		integer(I4B) :: ModelType
		logical(LGT) :: ErrorsKnown  
		logical(LGT) :: WellBlockIncl
	
	! This stores the number of data points so far modelled for each obs. point:
		integer(I4B), allocatable :: DoneDataPoints(:)
	
	! Flow data structures for f77 routines:
		integer(I4B):: NFLOWS
		real(DP), allocatable     :: TIM(:),FLO(:),ENT(:)
		logical(LGT) :: OnDeliv
		real(DP) :: ProdIndex,PCutoff
	
		contains
	
	!--------------------------------------------------------------------------------
	
		subroutine GetPumpingData(TestEndTime,StepFlows,RunToSS)
	!   Gets pumping data from the Problem_data module & fills the local
	!   arrays TIM,FLO & ENT.
	!   For pre-specified flows, data are generated from the pumping parameters.
	
	!   Argument variables:
			real(DP) :: TestEndTime
		logical(LGT):: StepFlows,RunToSS
	!   Local variables:
			integer(I4B) :: IntervalsPerCycle,i
		real(DP)     :: NumCycles,Amplitude,Period,Offset,Phase
		real(DP)     :: ConstantFlow,Duration,SamplingInterval,time
	
			if (RunToSS) then
			NFLOWS=1
		else
	!     Determine total number of pumping times:
			select case (Pump(1)%Scheme)
			case(0) ! Measured flows:
					NFLOWS=Pump(1)%NData	 
			case(1) ! Constant rate:
				ConstantFlow=PumpSchemeParams(1,1)
					NFLOWS=1 
				Pump(1)%StepFlows=.true.
			case(2) ! Sinusoidal:
					IntervalsPerCycle=20 ! Number of points generated per cycle of sinusoidal data
				Amplitude=PumpSchemeParams(1,1)
				Period=PumpSchemeParams(1,2)
				Offset=0.0_dp
				Phase=0.0_dp
				NumCycles=TestEndTime/Period
					NFLOWS=IntervalsPerCycle*NumCycles+1
				Pump(1)%StepFlows=.false.
				case(3) ! Offset sinusoidal:
					IntervalsPerCycle=20 
				Amplitude=PumpSchemeParams(1,1)
				Period=PumpSchemeParams(1,2)
				Offset=PumpSchemeParams(1,3)
				Phase=PumpSchemeParams(1,4)
				NumCycles=TestEndTime/Period
					NFLOWS=IntervalsPerCycle*NumCycles+1
				Pump(1)%StepFlows=.false.
			case(4) ! Step flow:
				ConstantFlow=PumpSchemeParams(1,1)
				Duration=PumpSchemeParams(1,2)
				NFLOWS=2
				Pump(1)%StepFlows=.true.
			end select
		end if
	
		allocate(TIM(NFLOWS),FLO(NFLOWS),ENT(NFLOWS))
	
			if (RunToSS) then ! Turn off pump for steady state run
			TIM(1)=0.0_dp
			FLO(1)=0.0_dp
		else
	!     Fill the local arrays:
			select case (Pump(1)%Scheme)
			case(0) ! Measured flows:
					TIM=PumpData(1,1:NFLOWS)%time
				FLO=PumpData(1,1:NFLOWS)%rate    	 
			case(1) ! Constant rate:
					TIM(1)=0.0_dp
				FLO(1)=ConstantFlow
			case(2,3) ! Sinusoidal:
					SamplingInterval=Period/IntervalsPerCycle
					time=0.0_dp
				do i=1,NFLOWS
					TIM(i)=time
				FLO(i)=Offset+Amplitude*dsin(2.0_dp*PI_d*time/Period+Phase)
				time=time+SamplingInterval
				end do
				case(4) ! Step:
				TIM(1)=0.0_dp
				FLO(1)=ConstantFlow
				TIM(2)=Duration
				FLO(2)=0.0_dp
			end select
			end if
	
	!   Fill enthalpy array- just one value:
		ENT=Pump(1)%Enthalpy
	
	!   Flag for flow interpolation type:
		StepFlows=Pump(1)%StepFlows
	
	!   Deliverability parameters:
			OnDeliv=Pump(1)%OnDeliverability
			ProdIndex=Pump(1)%ProdIndex
		PCutoff=Pump(1)%CutoffPressure
	
		return
		end subroutine GetPumpingData
	
	! ------------------------------------------------------------------------
	
		subroutine UpdateTimeStepSize(time,TestEndTime,NumIts,FlowIndex,&
		ResetTimeStepSize,StepFlows,RunToSS,dt)
	!   Updates time step size.
	
	!   Argument variables:
		real(DP), intent(in)    :: time,TestEndTime
		integer(I4B),intent(in) :: NumIts, FlowIndex
		logical(LGT),intent(inout):: ResetTimeStepSize
		logical(LGT),intent(in) :: StepFlows
		logical(LGT),intent(in) :: RunToSS
			real(DP), intent(inout) :: dt
	!   Local variables:
			real(DP) :: MinObsDT,MaxTimeStepSize
			integer(I4B),parameter :: MinNumTimeSteps=100
		real(DP),parameter     :: GrowthFactor=1.2_dp
		real(DP),parameter     :: SafetyFactor=1.05_dp
		integer(I4B),parameter :: SmallNumIts=4
		real(DP), parameter    :: StartTimeStepSize=1.E-2_dp !1.0_dp    
	
			MaxTimeStepSize=TestEndTime/MinNumTimeSteps
	
			if (time<TestEndTime) then
	
			if (ResetTimeStepSize) then ! Reset back to starting value
				dt=StartTimeStepSize
				ResetTimeStepSize=.false.
				else
	
	!       Increase dt if no. of iterations is small enough:
					if (NumIts<=SmallNumIts) then
						dt=min(GrowthFactor*dt,MaxTimeStepSize)
					end if
	
					if (.NOT.RunToSS) then
	!         Cut down to minimum observation data time interval:
						call FindMinObsTimeInterval(MinObsDT)
					dt=min(dt,MinObsDT)
			end if
	
	!       If passing a stepped flow rate, reduce timestep to hit it exactly,
	!       and reset following time step size to starting value:
					if ((StepFlows).and.(FlowIndex<NFLOWS).and.(time+dt>TIM(FlowIndex+1))) then
					if (TIM(FlowIndex+1)-time<StartTimeStepSize) then
					dt=TIM(FlowIndex+1)-time
					else
					dt=TIM(FlowIndex+1)-time
				end if
				ResetTimeStepSize=.true.
					end if
	
	!       Make sure TestEndTime isn't overshot:
				dt=min(dt,SafetyFactor*(TestEndTime-time))
	
				end if
	
		end if
	
			return
		end subroutine UpdateTimeStepSize
	
	! ------------------------------------------------------------------------
	
		subroutine FindMinObsTimeInterval(MinObsDT)
	!   Finds the current minimum time interval (over all observation points) between
	!   observations.
	
	!   Argument:
		real(DP), intent(inout):: MinObsDT
	!   Local variables:
			integer(I4B):: i, DataIndex
		real(DP) :: ObsTimeInterval
		real(DP),parameter:: big=1.E9
	
			MinObsDT=big
		do i=1,NObsPoints
			if ((DoneDataPoints(i)>0).and.(DoneDataPoints(i)<ObsPoint(i)%NData)) then
				DataIndex=ObsPoint(i)%DataIndex+DoneDataPoints(i)-1
					ObsTimeInterval=TestData(DataIndex+1)%time-TestData(DataIndex)%time
			MinObsDT=min(ObsTimeInterval,MinObsDT)	  
			end if
		end do
	
			return
		end subroutine FindMinObsTimeInterval
	
	! ------------------------------------------------------------------------
	
		subroutine ReduceTimestep(time,dt,FlowIndex,Q,H,ResetTimeStepSize)
	!   Discards current time step, reduces timestep size and resets simulation
	!   back to the start of the time step.
	
	!   Argument Variables:
			real(DP), intent(inout)    :: time,dt
		integer(I4B),intent(inout) :: FlowIndex
		real(DP), intent(inout)    :: Q,H
		logical(LGT), intent(in)   :: ResetTimeStepSize
			real(DP), parameter :: ReductionFactor=0.2_dp
	
			time=time-dt
			dt=ReductionFactor*dt
			time=time+dt
			call GetFlows(time,FlowIndex,Q,H,.false.,ResetTimeStepSize)
	
			return
		end subroutine ReduceTimestep
	
	! ------------------------------------------------------------------------
	
		subroutine GetFlows(time,FlowIndex,Q,H,SearchForwards,ResetTimeStepSize)
	
	!   Gets flow Q & enthalpy H from pumping record at 'time'. If 'time' is within
	!   the interval [TIM(FlowIndex),TIM(FlowIndex+1)), then FlowIndex is unchanged.
	!   Otherwise, a new FlowIndex is searched for, starting from the currrent value 
	!   (in the forwards direction if SearchForwards is true, and backwards otherwise).
	!   If this search fails, a search is made in the opposite direction.
	
	!   If ResetTimeStepSize is true, the current flow rate must be kept in use for 
	!   the current time step, so FlowIndex is not altered. 
	
	!   Argument Variables:
			real(DP), intent(in) :: time
		integer(I4B), intent(inout) :: FlowIndex
		real(DP), intent(out) :: Q,H
		logical(LGT), intent(in) :: SearchForwards, ResetTimeStepSize
	!   Locals:
			integer(I4B) :: Index, Increment
		logical(LGT) :: Found
			real(DP) :: theta,OneMinusTheta
		integer(I4B) :: StartFlowIndex
	
	!   Check to see if new flows need to be found:
	
			if (.NOT.ResetTimeStepSize) then
	
			 if ((time<TIM(FlowIndex)).or.(time>=TIM(min(FlowIndex+1,NFLOWS)))) then
	
			if (time<TIM(1)) then ! Pump is starting after test begins.
					FlowIndex=0
	
				else if (time>=TIM(NFLOWS)) then ! After last pump time.
					FlowIndex=NFLOWS
	
				else
				if (FlowIndex<1) then
				StartFlowIndex=1
				else if (FlowIndex>=NFLOWS) then
				StartFlowIndex=NFLOWS-1
			else
				StartFlowIndex=FlowIndex
			end if
	
				Index=StartFlowIndex
				if (SearchForwards) then 
					Increment=1
					else
					Increment=-1
				end if 
	
				Found=.false.
					do while (.NOT.Found)
						if((time>=TIM(Index)).and.(time<TIM(Index+1))) then
					Found=.true.
					else
					Index=Index+Increment
					! Try other direction if going outside range of pumping times:
					if ((Index<1).or.(Index>=NFLOWS)) then
						Index=StartFlowIndex
						Increment=-Increment
					end if		     
						end if
					end do
	
				FlowIndex=Index
	
			end if
	
			 end if ! Found new flow index.
	
	!    Interpolate flow and enthalpy values:
	
			 if (FlowIndex<1) then
				Q=0.0_dp
				H=FLO(1)
			 else if (FlowIndex>=NFLOWS) then
				Q=FLO(NFLOWS)
				H=ENT(NFLOWS)
		 else
			if (Pump(1)%StepFlows) then
				! If step flows then use value at FlowIndex:
				theta=0.0_dp
				else
					theta=(time-TIM(FlowIndex))/(TIM(FlowIndex+1)-TIM(FlowIndex))
			end if
			OneMinusTheta=1.0_dp-theta
				Q=FLO(FlowIndex)*OneMinusTheta + FLO(FlowIndex+1)*theta
				H=ENT(FlowIndex)*OneMinusTheta + ENT(FlowIndex+1)*theta
			 end if
	
		end if
	
		return	  
		end subroutine GetFlows
			
	!--------------------------------------------------------------------------------
	
		subroutine DestroyProblemDataArrays
			deallocate (Pump,ObsPoint,TestData,PumpData)
		deallocate (PumpSchemeParams)
		deallocate (ReservoirCondition,FixedParameter)
			return
		end subroutine DestroyProblemDataArrays
	
	!--------------------------------------------------------------------------------
	
		subroutine DestroyInterpolationArrays
		deallocate(TIM,FLO,ENT)
			return
		end subroutine DestroyInterpolationArrays
	
	!--------------------------------------------------------------------------------
	
	end module problem_data
	