module convert_input

contains

  subroutine convert_input_variables(InputNPumps, PumpingScheme, &
        NumPumpTimes, InputMaxNPumpTimes, PumpTime, PumpRate, PumpEnthalpy, &
		PumpStepFlows, PumpOnDeliv, PumpObsPoint, PumpProdIndex, PumpCutoffPressure, &
		InputMaxNPumpSchemeParams, InputPumpSchemeParams, &
		InputNObsPoints, ObsPointProperty, &
		ObsPointNumData, ObsPointPosition, ObsPointError, ObsPointWeight,&
		ObsPointDataOffset, InputTotalNData, ObservationTime, &
		ObservationValue, InputModelType, InputNReservoirConditions, &
		InputNFixedParameters, InputNVariables, InputReservoirCondition, &
		InputFixedParameter, ParameterBound, ParameterScale, &
		ParameterLogScale)

!   This routine takes the intermediary input arrays, passed to the fortran90 back end
!   by the Delphi front end, and uses them to fill the data structures in the 'problem_data'
!   and 'variable_parameters' modules.  The models then use these data structures; the 
!   intermediary arrays are not used further.
        
    use variable_types
    use problem_data
    use variable_parameters
	use models
    implicit none

!   Argument variables:
    integer(I4B), intent(IN) :: InputNPumps
    integer(I4B), intent(IN) :: PumpingScheme(1:InputNPumps)
    integer(I4B), intent(IN) :: NumPumpTimes(1:InputNPumps)
    integer, intent(IN)      :: InputMaxNPumpTimes
    real(DP), intent(IN)     :: PumpTime(1:InputMaxNPumpTimes,1:InputNPumps)
    real(DP), intent(IN)     :: PumpRate(1:InputMaxNPumpTimes,1:InputNPumps)
    real(DP), intent(IN)     :: PumpEnthalpy(1:InputNPumps)
	integer(I4B), intent(IN) :: PumpStepFlows(1:InputNPumps)
	integer(I4B), intent(IN) :: PumpOnDeliv(1:InputNPumps)
    integer(I4B), intent(IN) :: PumpObsPoint(1:InputNPumps)
    real(DP), intent(IN)     :: PumpProdIndex(1:InputNPumps)
    real(DP), intent(IN)     :: PumpCutoffPressure(1:InputNPumps)
	integer(I4B), intent(IN) :: InputMaxNPumpSchemeParams
	real(DP), intent(IN)     :: InputPumpSchemeParams(1:InputMaxNPumpSchemeParams,1:InputNPumps)
    integer(I4B), intent(IN) :: InputNObsPoints
    integer(I4B), intent(IN) :: ObsPointProperty(1:InputNObsPoints)
    integer(I4B), intent(IN) :: ObsPointNumData(1:InputNObsPoints)
    real(DP), intent(IN)     :: ObsPointPosition(1:InputNObsPoints,1:3)   
    real(DP), intent(IN)     :: ObsPointError(1:InputNObsPoints)
    real(DP), intent(IN)     :: ObsPointWeight(1:InputNObsPoints)
    real(DP), intent(IN)     :: ObsPointDataOffset(1:InputNObsPoints)
    integer(I4B), intent(IN) :: InputTotalNData
    real(DP), intent(IN)     :: ObservationTime(1:InputTotalNData)
    real(DP), intent(IN)     :: ObservationValue(1:InputTotalNData)
    integer(I4B), intent(IN) :: InputModelType
    integer(I4B), intent(IN) :: InputNReservoirConditions
    integer(I4B), intent(IN) :: InputNFixedParameters
    integer(I4B), intent(IN) :: InputNVariables
    real(DP), intent(IN)     :: InputReservoirCondition(1:InputNReservoirConditions)
    real(DP), intent(IN)     :: InputFixedParameter(1:InputNFixedParameters)
    real(DP), intent(IN)     :: ParameterBound(1:2,1:InputNVariables)
    real(DP), intent(IN)     :: ParameterScale(1:InputNVariables)
    integer(I4B), intent(IN) :: ParameterLogScale(1:InputNVariables)

!   Locals:
    integer(I4B) :: PumpNo,ObsPointNo,i,lower,upper
    real(DP)     :: DefaultError

!   Arbitrary default observation error (in case errors are not specified):
    DefaultError=1.0E5_DP

!   Copy problem dimensions from 'dummy' intermediary variables:
    NPumps=InputNPumps
    NObsPoints=InputNObsPoints 
    TotalNData=InputTotalNData 
    MaxNPumpTimes=InputMaxNPumpTimes
	MaxNPumpSchemeParams=InputMaxNPumpSchemeParams
    ModelType=InputModelType
	call SetWellBlockIncl
    NReservoirConditions=InputNReservoirConditions
    NFixedParameters=InputNFixedParameters
    NVariables=InputNVariables
    
!   Allocate main arrays:
    allocate(Pump(NPumps),ObsPoint(NObsPoints),TestData(TotalNData))
    allocate(PumpData(NPumps,MaxNPumpTimes))
	allocate(PumpSchemeParams(NPumps,MaxNPumpSchemeParams))
    allocate(ReservoirCondition(NReservoirConditions),FixedParameter(NFixedParameters))
    allocate(scale_factor(NVariables),lowerbounds(NVariables))
    allocate(upperbounds(NVariables),logscale(NVariables))

!   Copy fixed model parameters & reservoir conditions:
    FixedParameter=InputFixedParameter
  	ReservoirCondition=InputReservoirCondition

!   Fill in pump data structures:
    Pump(:)%Scheme=PumpingScheme(:)
	Pump(:)%Enthalpy=PumpEnthalpy(:)
    Pump(:)%NData=NumPumpTimes(:)
	Pump(:)%OnDeliverability=PumpOnDeliv(:)
	Pump(:)%ProdIndex=PumpProdIndex(:)
	Pump(:)%CutoffPressure=PumpCutoffPressure
	ObsPoint(:)%IsPumpObsPoint=.false. ! (initialise)
	ObsPoint(:)%PumpNo=0
    do PumpNo=1, NPumps
	  if (PumpStepFlows(PumpNo)==1) then
	    Pump(PumpNo)%StepFlows=.true.
      else
	    Pump(PumpNo)%StepFlows=.false.
      end if	      
      do i=1,Pump(PumpNo)%NData
        PumpData(PumpNo,i)%time=PumpTime(i,PumpNo)
        PumpData(PumpNo,i)%rate=PumpRate(i,PumpNo)
      end do
	  if (Pump(PumpNo)%OnDeliverability) then
	    ObsPoint(PumpObsPoint(PumpNo))%IsPumpObsPoint=.true.
        ObsPoint(PumpObsPoint(PumpNo))%PumpNo=PumpNo
	  end if
    end do
    PumpSchemeParams=transpose(InputPumpSchemeParams)
    
!   Fill in ObsPoint data structures:
    ObsPoint%Property=ObsPointProperty
    ObsPoint%NData=ObsPointNumData
    ObsPoint%Error=ObsPointError
	ObsPoint%Weight=ObsPointWeight
	ObsPoint%DataOffset=ObsPointDataOffset
    do ObsPointNo=1,NObsPoints
      ObsPoint(ObsPointNo)%Position%x=ObsPointPosition(ObsPointNo,:)
	end do      
    
!   Construct ObsPoint%DataIndex array- this stores the starting index in the
!   TestData array of the observation data for each ObsPoint:
    ObsPoint(:)%DataIndex=0
    ObsPoint(1)%DataIndex=1
    if (NObsPoints>1) then
      do ObsPointNo=2,NObsPoints
        ObsPoint(ObsPointNo)%DataIndex= &
		  ObsPoint(ObsPointNo-1)%DataIndex+ObsPoint(ObsPointNo-1)%NData
      end do
    end if
    
!   Fill in main observation data structure:
    TestData%time=ObservationTime
    TestData%value=ObservationValue
    
!  Work out whether or not observation errors have been specified for all
!  data points:
    if (minval(ObsPoint%Error)>small_D) then
      ErrorsKnown=.true.
    else
      ErrorsKnown=.false.
    end if

!   Fill in observation errors and weights:
    do ObsPointNo=1,NObsPoints
      lower=ObsPoint(ObsPointNo)%DataIndex
      upper=lower+ObsPoint(ObsPointNo)%NData-1
	  if (ErrorsKnown) then
        TestData(lower:upper)%error=ObsPoint(ObsPointNo)%error
	  else
		TestData(lower:upper)%error=DefaultError
	  end if
	  TestData(lower:upper)%weight=ObsPoint(ObsPointNo)%Weight 
    end do
    
!  Fill in 'variable_parameters' data structures:
    lowerbounds=ParameterBound(1,:)
    upperbounds=ParameterBound(2,:)
    
!  Variable scaling parameters:
    scale_factor=ParameterScale
    do i=1, NVariables
      if (ParameterLogScale(i)==1) then
        logscale(i)=.true.
      else
        logscale(i)=.false.
      end if
    end do
    
!   Factor for scaling down objective function:
!    ObjectiveScale=real(TotalNData)
    ObjectiveScale=1.0_dp
!   (Took this scaling out as it is not critical, and makes output
!   of chisquare values from nl2sno difficult.)

  return  
  end subroutine convert_input_variables

end module convert_input
