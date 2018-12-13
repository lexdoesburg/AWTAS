module NumericalSimulatorMINC

  use utility_functions
  use variable_types
  use problem_data
  use NumericalSimulatorMINC_Routines
  implicit none

! nb: MR is the number of grid blocks in the radial direction,
! NZ is the number in the vertical direction.

  integer(I4B)              :: MR,NZ
  real(DP), allocatable     :: PERFR(:),PORFR(:)
  real(DP), allocatable     :: PERMA(:),PORMA(:),PRAT(:)
  real(DP), allocatable     :: V(:,:)       
  real(DP), allocatable     :: P(:,:),T(:,:),SV(:,:),X(:,:)
  real(DP), allocatable     :: POLD(:,:),TOLD(:,:),SVOLD(:,:)
  real(DP), allocatable     :: XOLD(:,:)
  real(DP), allocatable     :: BMOLD(:,:),BEOLD(:,:)
  integer(I4B), allocatable :: IPH(:,:),IPHOLD(:,:)
  real(DP), allocatable     :: DR(:),AR(:),DELR(:)
  real(DP), allocatable     :: DZ(:),AZ(:,:),DELZ(:)
  real(DP), allocatable     :: XX(:,:),RR(:)
  real(DP), allocatable     :: Alpha(:),Beta(:)

! This stores the radial grid block each obs. point is in:
  integer(I4B), allocatable :: ObsPointGridBlock(:)
! Flowing enthalpy variables:
  real(DP)                  :: HF, HFOLD
! Weighting factors for extrapolating values at the sandface:
  real(DP):: SandFaceWeight1,SandFaceWeight2,SandFaceWeight3

  contains

! ------------------------------------------------------------------------------------

  function NumericalSolutionMINC(CRF,CRM,CONDF,CONDM,RHORF,RHORM,&
      COMPF,COMPM,COMTF,COMTM,AAAF,AAAM,InitialPressure,InitialX,&
      kFracture,kMatrix,phiFracture,phiMatrix,updatemodelprogress)

    use ModelProgress

!   Argument Variables:
    real(DP),intent(in) :: CRF,CRM,CONDF,CONDM,RHORF,RHORM
    real(DP),intent(in) :: COMPF,COMPM,COMTF,COMTM,AAAF,AAAM
    real(DP),intent(in) :: InitialPressure,InitialX
    real(DP),intent(in) :: kFracture,kMatrix,phiFracture,phiMatrix
    real(DP) :: NumericalSolutionMINC(TotalNData)

!   Local variables:
    real(DP)    :: P0(MR,NZ),T0(MR,NZ)
    integer(I4B):: I,IT,J
    real(DP)    :: dt,HIN,QMM,RMAX,RMAXMM,RMAXME,RMAXFM,RMAXFE
    real(DP)    :: time,TestEndTime
    integer(I4B):: NTimeSteps, FlowIndex
    logical(LGT):: ResetTimeStepSize, StepFlows
    integer(I4B):: UpdateCounter
    integer(I4B), parameter :: UpdateInterval=20 ! #iterations before updating progress
    real(DP), parameter     :: RTOL=1.E-5_dp  ! Newton-Raphson iteration tolerance
    integer(I4B), parameter :: IMAX=10        ! Max permitted no. of iterations
    real(DP), parameter     :: EPS=1.0E-6_dp  ! A smallish number

    external updatemodelprogress

    NumModelRuns=NumModelRuns+1

    TestEndTime=maxval(TestData(:)%time) ! Last observation time

    call AssignBlockProperties(kFracture,kMatrix,phiFracture,phiMatrix)
    call AssignInitialConditions(InitialPressure,InitialX)
    call GetPumpingData(TestEndTime,StepFlows,.false.)

!   Initialise result:
    NumericalSolutionMINC=TestData%ModelledValue

    IGOOD=0     ! 'ok' parameter passed from thermo subroutines

!   Read in initial conditions and decide on the phase conditions:
    call INIT1(MR,NZ,POLD,XOLD,TOLD,SVOLD,IPHOLD)
    P0=POLD
    T0=TOLD

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
      call INIT2(MR,NZ,POLD,XOLD,TOLD,SVOLD,IPHOLD,&
               BMOLD,BEOLD,PORFR,RHORF,CRF,PORMA,RHORM,CRM,&
               P0,T0,COMPF,COMTF,AAAF,COMPM,COMTM,AAAM)

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
      call SOLVE(P,T,SV,X,IPH,HF,XX,BMOLD,BEOLD,MR,NZ,PORFR,PORMA,&
               PERFR,PERMA,AR,AZ,V,CRF,CRM,RHORF,RHORM,CONDF,CONDM,&
               QMM,HIN,DT,DELR,DELZ,RMAXMM,RMAXME,RMAXFM,RMAXFE,&
               P0,T0,COMPF,COMTF,AAAF,COMPM,COMTM,AAAM,PRAT)

        if (IGOOD>0) then  ! Problems... reduce timestep:
        call ReduceTimestep(time,dt,FlowIndex,QMM,HIN,ResetTimeStepSize)
        IGOOD=0
        GO TO 10 
      end if      
           
!     Add on the solution increments XX and make any necessary phase changes:
      call UPDATE(MR,NZ,IPH,P,T,SV,X,XX,EPS)

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
      call InterpolateResponse(time,dt,TestEndTime,NumericalSolutionMINC)

!     Update 'old' variables:
      POLD=P         
      TOLD=T         
      SVOLD=SV         
      XOLD=X         
      IPHOLD=IPH
      HFOLD=HF

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
  end function NumericalSolutionMINC  

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
!         case(0) (flow) isn't simulated by this model.
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
    real(DP), intent(in):: Values(MR,NZ)
    integer(I4B), intent(in):: BlockNo
    real(DP):: InterpolateBlockValue

    if (BlockNo==0) then ! Sandface extrapolation:
      InterpolateBlockValue=SandFaceWeight1*Values(1,1)+ &
        SandFaceWeight2*Values(2,1)+SandFaceWeight3*Values(3,1)
    else
      InterpolateBlockValue=Values(BlockNo,1)
    end if

    return
  end function InterpolateBlockValue

! ----------------------------------------------------------------------------------

  subroutine CalculateSandfaceWeights
!   Calculates weighting factors for sandface extrapolation.  This is done by
!   fitting a quadratic through the first 3 block centres to fit the value at the
!   sandface.

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

  subroutine SetupGrid(LayerThickness,ActionWellRadius,FractureVolFraction,&
      FractureSpacing,NumFracturePlanes)
!   Sets up radial MINC grid for the simulation.

!   Argument variables:
    real(DP), intent(in) :: LayerThickness, ActionWellRadius
    real(DP), intent(in) :: FractureVolFraction,FractureSpacing
    real(DP), intent(in) :: NumFracturePlanes

    call SetupGridGeometry(ActionWellRadius,LayerThickness,&
      FractureVolFraction,FractureSpacing,NumFracturePlanes) 
!   Work out observation point positions on the grid:
    call GetObsPointGridPositions(ActionWellRadius)
!   Calculate weights for sandface extrapolation:
    call CalculateSandfaceWeights
    
    return
  end subroutine SetupGrid
  
! ------------------------------------------------------------------------

  subroutine SetupGridGeometry(ActionWellRadius,LayerThickness,&
      FractureVolFraction,FractureSpacing,NFP) 

!   Sets up block sizes & geometry parameters (areas & volumes) for
!   MINC grid with NFP fracture planes.

!   Argument variable:
    real(DP), intent(in)   :: ActionWellRadius,LayerThickness
    real(DP), intent(in)   :: FractureVolFraction,FractureSpacing
    real(DP), intent(in)   :: NFP
!   Local variables:
    real(DP),parameter     :: ConstDR=0.10_dp
    real(DP),parameter     :: GrowthFactor=1.2_dp
    integer(I4B),parameter :: NConstBlocks=10
    integer(I4B)           :: i,j
    real(DP),allocatable   :: Z(:)
    real(DP)               :: R0,R1,VTot,ARTot,NFPInv

!   No. of blocks:
    MR=100
    NZ=8

!   Allocate main local arrays:
    allocate (PERFR(MR),PERMA(MR))
    allocate (PORFR(MR),PORMA(MR),PRAT(MR))
    allocate (V(MR,NZ),P(MR,NZ),T(MR,NZ),SV(MR,NZ),X(MR,NZ))
    allocate (POLD(MR,NZ),TOLD(MR,NZ),SVOLD(MR,NZ))
    allocate (XOLD(MR,NZ))
    allocate (BMOLD(MR,NZ),BEOLD(MR,NZ))
    allocate (IPH(MR,NZ),IPHOLD(MR,NZ))
    allocate (DR(MR),AR(MR),DELR(MR))
    allocate (DZ(NZ),AZ(MR,NZ-1),DELZ(NZ-1))
    allocate (XX(2*MR,NZ),RR(2*MR))
    allocate (Alpha(NZ),Beta(NZ))

!   Determine the radial block sizes:
    do i=1,NConstBlocks
      DR(i)=ConstDR
    end do
    do i=1+NConstBlocks,MR
      DR(i)=GrowthFactor*DR(i-1)
    end do

!   Calculate volume fractions for matrix blocks:
    call GetMatrixVolumeFractions(NZ,FractureVolFraction,Alpha,Beta)

!   Determine vertical block sizes:
    allocate (Z(NZ+1))
    NFPInv=1._dp/NFP
    Z(1)=0._dp
    do j=2,NZ+1
      Z(j)=0.5_dp*FractureSpacing*(1._dp-(1._dp-Beta(j-1))**NFPInv)
    end do
    do j=1,NZ
      DZ(j)=Z(j+1)-Z(j)
    end do

!   Vertical block centre spacings:    
    do j=1,NZ-2
      DELZ(j)=0.5_dp*(DZ(j)+DZ(j+1))
    end do
    DELZ(NZ-1)=0.5_dp*DZ(NZ-1)+DZ(NZ)/(NFP+2._dp)

    R1=ActionWellRadius
    do i=1,MR

      R0=R1
      R1=R0+DR(i)
      if (i<MR) then
        DELR(i)=0.5_dp*(DR(i)+DR(i+1))
      end if
      VTot=PI_D*(R1+R0)*(R1-R0)*LayerThickness
      ARTot=2._dp*PI_D*R1*LayerThickness
      AR(i)=ARTot

      do j=1,NZ
        V(i,j)=Alpha(j)*VTot
        if (j<NZ) then
          AZ(i,j)=(1._dp-Alpha(1))*VTot*2._dp*NFP/FractureSpacing*&
            (1._dp-Beta(j))**(1._dp-NFPInv)
        end if
      end do

    end do

    deallocate (Z)

    return
  end subroutine SetupGridGeometry

! ------------------------------------------------------------------------

  subroutine GetMatrixVolumeFractions(P,Alpha1,Alpha,Beta)

  ! Sets up volume fractions for matrix blocks (alpha) and cumulative
  ! volume fractions (beta).  Alpha(2) is set equal to Alpha1, the fracture
  ! volume fraction.  Volume fractions increase geometrically thereafter, 
  ! according to the growth factor w.

!   Argument variables:
    integer(I4B),intent(in):: P
    real(DP), intent(in)   :: Alpha1
    real(DP)               :: Alpha(P),Beta(P)
!   Locals:
    real(DP) :: w,delw
    real(DP) :: rhs,f,fdash,wPm2,Onemw,OnemwPm2,AlphaSum
    integer(I4B) :: i
    logical(LGT) :: Converged
    real(DP),parameter     :: Startingw=2.0_dp
    integer(I4B),parameter :: MaxIts=20
    real(DP),parameter     :: tolerance=1.E-8_dp

    ! Calculate growth factor w to fit specified # blocks:
    w=Startingw
    Converged=.false.
    do i=1,MaxIts
      wPm2=w**(P-2)
      Onemw=1._dp-w
      OnemwPm2=1._dp-wPm2
      f=Alpha1*(2._dp+w*OnemwPm2/Onemw)-1._dp
      fdash=Alpha1*(OnemwPm2/Onemw+(w*OnemwPm2-(P-2)*wPm2*Onemw)/(Onemw*Onemw))
      delw=-f/fdash
      w=w+delw
      if (dabs(delw)<=tolerance) then
        Converged=.true.
        exit
      end if 
    end do

    ! Fill in Alpha array:
    Alpha(1)=Alpha1
    Alpha(2)=Alpha(1)
    do i=3,P
      Alpha(i)=w*Alpha(i-1)
    end do

    ! Fill in Beta array:
    Beta(1)=0.0_dp
    AlphaSum=0.0_dp
    do i=2,P
      AlphaSum=AlphaSum+Alpha(i)
      Beta(i)=AlphaSum/(1._dp-Alpha1)
    end do

    return
  end subroutine GetMatrixVolumeFractions

! ------------------------------------------------------------------------

  subroutine GetObsPointGridPositions(rw)

!   Argument:
    real(DP), intent(in):: rw
!   Local variables:
    integer(I4B):: i,block
    real(DP)    :: ObsPointR,R0,R1

    allocate(ObsPointGridBlock(NObsPoints))

    do i=1,NObsPoints

!     Radial position of observation point:
      ObsPointR=ObsPoint(i)%Position%x(1)

      if ((rw>small_d).and.(ObsPointR<=rw)) then
        ObsPointGridBlock(i)=0
      else
!       Search grid:
        R1=rw
        do block=1,MR
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

  subroutine AssignBlockProperties(kFracture,kMatrix,phiFracture,phiMatrix)
!   Argument Variables:
    real(DP), intent(in) :: kFracture,kMatrix,phiFracture,phiMatrix

    PERFR=kFracture
    PERMA=kMatrix

    PORFR=phiFracture
    PORMA=phiMatrix

    PRAT=PERFR/PERMA

    return
  end subroutine AssignBlockProperties

! ------------------------------------------------------------------------
      
  subroutine AssignInitialConditions(InitialPressure,InitialX)
!   Argument Variables:
    real(DP),intent(in) :: InitialPressure,InitialX

    POLD=InitialPressure
    XOLD=InitialX
    HFOLD=0.0_dp

    return
  end subroutine AssignInitialConditions

! ------------------------------------------------------------------------

  subroutine DestroyGridArrays
    deallocate (PERFR,PERMA)
    deallocate (PORFR,PORMA,PRAT)
    deallocate (V,P,T,SV,X)
    deallocate (POLD,TOLD,SVOLD,XOLD)
    deallocate (BMOLD,BEOLD,IPH,IPHOLD)
    deallocate (DR,AR,DELR)
    deallocate (DZ,AZ,DELZ)
    deallocate (XX,RR)
    deallocate (Alpha,Beta)
    deallocate(ObsPointGridBlock)
    return
  end subroutine DestroyGridArrays

! ------------------------------------------------------------------------------------

end module NumericalSimulatorMINC
