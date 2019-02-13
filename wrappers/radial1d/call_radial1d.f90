module call_radial1d

    use variable_types
    use problem_data
    ! use variable_parameters
    ! use models
  
    implicit none
  
    contains
    

    !------------------------------------------------------------------------------------------------------------------------------------------------------------------
    ! Calling simulator
    !------------------------------------------------------------------------------------------------------------------------------------------------------------------

    subroutine Radial1D(Porosity, Permeability, LayerThickness, ActionWellRadius, RockSpecificHeat, RockHeatConductivity, RockDensity, InitialPressure, InitialX, &
                        InjectionWell, InjectionEnthalpy, NumPumpTimes, NumObservationPoints, TotalNumData, PumpingScheme, MassFlowrate, FlowDuration, PumpTime, PumpRate, &
                        Time, ObsPointRadialLocation, ObsPointNumData, ObsPointProperty, Deliverability, ProductionIndex, CutoffPressure, updatemodelprogress, &
                        ExecutionFlag, ModelledValue)

        ! Modules
        use HomogeneousPorousSimulator

        ! Input arguments
        real(DP), intent(in) :: Porosity, Permeability ! Variables
        real(DP), intent(in) :: LayerThickness, ActionWellRadius, RockSpecificHeat, RockHeatConductivity, RockDensity, RockCompressibility ! Fixed parameters
        real(DP), intent(in) :: InitialPressure, InitialX ! Initial reservoir conditions
        logical(LGT), intent(in) :: InjectionWell ! True if the well is an injection well, false if production
        real(DP), intent(in) :: InjectionEnthalpy
        integer(I4B), intent(in) :: NumPumpTimes, NumObservationPoints, TotalNumData, PumpingScheme
        real(DP), intent(in) :: MassFlowrate, FlowDuration ! Flow parameters for step flow and constant flow
        real(DP), intent(in), dimension(NumPumpTimes) :: PumpTime, PumpRate ! Flow parameters for measured flows
        real(DP), intent(in), dimension(TotalNumData) :: Time
        integer(I4B), intent(in), dimension(NumObservationPoints) :: ObsPointRadialLocation, ObsPointNumData, ObsPointProperty ! Observation point info arrays
        logical(LGT), intent(in) :: Deliverability
        real(DP), intent(in) :: ProductionIndex, CutoffPressure

        ! Outputs
        integer(I4B), intent(out) :: ExecutionFlag ! Flag which signals whether the program executed properly or not. Flag = 1 if everything ran fine, 2 if missing input, 3 if simulator did not run as expected
        real(DP), intent(out), dimension(TotalNumData) :: ModelledValue

        ! Local variables
        integer(I4B) :: i

        ! External subroutines
        external updatemodelprogress

        ! Assign key variables in problem_data module
        NPumps = 1 ! Only one pump (ask Mike if you can have more in this model) Could change to an input of number of pumps
        MaxNPumpTimes = NumPumpTimes
        NObsPoints = NumObservationPoints
        TotalNData = TotalNumData

        if (PumpingScheme == 4) then
            MaxNPumpSchemeParams = 2
        else  ! PumpingScheme == 1 (don't really need it if not)
            MaxNPumpSchemeParams = 1
        end if

        ModelType = 1
        call SetWellBlockIncl

        select case(ModelType)
        case(1) ! Homogeneous porous layer model
            NVariables = 2
            NFixedParameters = 7
            NReservoirConditions = 2
        end select

        ! Allocate arrays declared in problem_data
        allocate(Pump(NPumps))
        allocate(ObsPoint(NObsPoints))
        allocate(PumpData(NPumps,MaxNPumpTimes))
        allocate(PumpSchemeParams(NPumps,MaxNPumpSchemeParams))
        allocate(ReservoirCondition(NReservoirConditions))
        allocate(FixedParameter(NFixedParameters))
        allocate(TestData(TotalNData))
        allocate(variable(NVariables))

        ! Fill variable, parameter and reservoir condition arrays
        call FillParametersArray(LayerThickness, ActionWellRadius, RockSpecificHeat, RockHeatConductivity, RockDensity, RockCompressibility)
        call FillReservoirConditionsArray(InitialPressure, InitialX)
        call FillVariablesArray(Porosity, Permeability)

        ! Fill pump array
        Pump(1)%Scheme = PumpingScheme
        Pump(1)%NData = NumPumpTimes
        select case(PumpingScheme)
        case default
            ! Should probably check inputs are fine before array allocation
            ExecutionFlag = 2
            call DeallocateArrays
            return
        case(0) ! Measured flows (flowrates differ at different times)
            Pump(1)%StepFlows = .True.
            do i = 1, Pump(1)%NData
                PumpData(1,i)%time = PumpTime(i)
                PumpData(1,i)%rate = PumpRate(i)
            end do
        case(1) ! Constant rate
            PumpSchemeParams(1,1) = MassFlowrate
            Pump(1)%StepFlows = .True. ! could be false tbh - since its not actually a step
        case(4) ! Step flow (flow switched off after some time)
            PumpSchemeParams(1,1) = MassFlowrate
            PumpSchemeParams(1,2) = FlowDuration
            Pump(1)%StepFlows = .True.
        end select

        if (InjectionWell == .True.) then
            Pump(1)%Enthalpy = InjectionEnthalpy
        else
            Pump(1)%Enthalpy = 0.0_dp          
        end if

        if (Deliverability == .True.) then
            Pump(1)%OnDeliverability = .True.
            Pump(1)%ProdIndex = ProductionIndex
            Pump(1)%CutoffPressure = CutoffPressure

            ! Assign variables in problem_data module
            OnDeliv = .True.
            ProdIndex = ProductionIndex
            PCutoff = CutoffPressure
        else
            Pump(1)%OnDeliverability = .False.
            Pump(1)%ProdIndex = 0.0_dp ! check reasonable value
            Pump(1)%CutoffPressure = 0.0_dp ! check reasonable value

            ! Assign variables in problem_data module
            OnDeliv = .False.
            ProdIndex = 0.0_dp
            PCutoff = 0.0_dp          
        end if

        ! Observation points
        ObsPoint(:)%Error = 0.0_dp ! could remove from problem data
        ObsPoint(:)%Weight = 1.0_dp ! could remove from problem data
        ObsPoint(:)%DataOffset = 0.0_dp ! could remove from problem data
        ObsPoint(:)%Position%x(2) = 0.0_dp
        ObsPoint(:)%Position%x(3) = 0.0_dp
        ObsPoint(:)%DataIndex = 1
        ObsPoint(:)%IsPumpObsPoint = .False.
        ObsPoint(:)%PumpNo = 0
        do i = 1, NObsPoints
            ObsPoint(i)%Position%x(1) = ObsPointRadialLocation(i)
            ObsPoint(i)%NData = ObsPointNumData(i)
            if (NObsPoints > 1 .and. i > 1) then
                ObsPoint(i)%DataIndex = ObsPoint(i-1)%DataIndex + ObsPoint(i-1)%NData
            end if

            ObsPoint(i)%Property = ObsPointProperty(i) ! 0 = flows on deliverability, 1 = pressure, 2 = temperature, 3 = enthalpy
            if (ObsPoint(i)%Property == 0) then ! Deliverability
                ObsPoint(i)%IsPumpObsPoint = .True. ! Not sure
                ObsPoint(i)%PumpNo = 1 ! We only have one pump currently
            end if
        end do

        ! Testdata array
        TestData(:)%error = 0.0_dp ! could remove from problem data
        TestData(:)%weight = 1.0_dp ! could remove from problem data
        TestData(:)%value = 0.0_dp ! could remove from problem data
        TestData(:)%ModelledValue = 0.0_dp ! could remove from problem data
        TestData%time = Time ! this is in the structure of observation points stacked on top of each other ! could be improved if observation points themselves just had time arrays possibly

        ! Calculate modelled value
        ModelledValue = HomogeneousPorous(variable,updatemodelprogress)

        ! Deallocate arrays
        deallocate(Pump)
        deallocate(ObsPoint)
        deallocate(PumpData)
        deallocate(PumpSchemeParams)
        deallocate(ReservoirCondition)
        deallocate(FixedParameter)
        deallocate(TestData)
        deallocate(variable)
        
        return
    end subroutine Radial1D
  
    !------------------------------------------------------------------------------------------------------------------------------------------------------------------
    ! Helper routines
    !------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    subroutine FillParametersArray(LayerThickness, ActionWellRadius, RockSpecificHeat, RockHeatConductivity, RockDensity, RockCompressibility)
        real(DP), intent(in) :: LayerThickness, ActionWellRadius, RockSpecificHeat, RockHeatConductivity, RockDensity, RockCompressibility
        FixedParameter(1) = 0.0_dp
        FixedParameter(2) = LayerThickness
        FixedParameter(3) = ActionWellRadius
        FixedParameter(4) = RockSpecificHeat
        FixedParameter(5) = RockHeatConductivity
        FixedParameter(6) = RockDensity
        FixedParameter(7) = RockCompressibility
        return
    end subroutine FillParametersArray

    subroutine FillReservoirConditionsArray(InitialPressure, InitialX)
        real(DP), intent(in) :: InitialPressure, InitialX
        ReservoirCondition(1) = InitialPressure
        ReservoirCondition(2) = InitialX
        return
    end subroutine FillReservoirConditionsArray

    subroutine FillVariablesArray(Porosity, Permeability)
        real(DP), intent(in) :: Porosity, Permeability
        variable(1) = Permeability
        variable(2) = Porosity
        return
    end subroutine FillVariablesArray

    !------------------------------------------------------------------------------------------------------------------------------------------------------------------
    ! Dummy routine for model progress
    !------------------------------------------------------------------------------------------------------------------------------------------------------------------
    subroutine updatemodelprogress(nummodelruns,timestepsize, &
        progressfraction, analysisstopped)
        ! This is a dummy routine to pass into the fast solver.
        use variable_types
        integer(I4B), intent(in) :: nummodelruns
        real(DP), intent(in) :: timestepsize,progressfraction
        logical(LGT), intent(inout) :: analysisstopped
        return
    end subroutine updatemodelprogress

  end module call_radial1d
  