! For theis solution
! gfortran -o pseudodata pseudodata.f90 variable_types.o problem_data.o variable_parameters.o models.o noise.o theis_solution.o utility_functions.o

! For homogeneous porous simulator
! gfortran -o pseudodata pseudodata.f90 variable_types.o problem_data.o variable_parameters.o models.o noise.o theis_solution.o utility_functions.o homogeneousporous.o thermodynamics.o numericalsimulator1d_routines.o NumericalSimulator1D.o gammafunction.o matrixsolvers.o modelprogress.o
program PseudoData

! Generates artificial test data by adding random noise to the results
! from an AWTAS model.

! Assume 1 pump, and 1 observation point per observation well.

use variable_types
use problem_data
use variable_parameters
use models
use noise

implicit none

! Global variables:
integer(I4B) :: out,dat
character(len=60) :: outfile
integer(I4B) :: MaxTotalNData
integer(I4B) :: ObsPointNo,i,DataNo
real(DP) :: t,t0,dt,t1
real(DP), allocatable :: variable(:)
real(DP) :: k,nu,phi,rho,c,b,Q0,P0,rw,period
real(DP) :: Pressure0,X0,thick,CR,COND,RHOR,COMP,COMT,AAA
real(DP) :: CRF,CONDF,RHORF,COMPF,COMTF,AAAF
real(DP) :: CRM,CONDM,RHORM,COMPM,COMTM,AAAM
real(DP) :: fracth,fracsp
type(datapoint),allocatable      :: ReadData(:)
integer(I4B) :: lower,upper
character    :: DataFileName*80
integer(I4B) :: StepRates
real(DP) :: WellRadius

external updatemodelprogress

dat=10
out=11

MaxTotalNData=5000

write (*,*) ' Output file name:'
read (*,'(a)') outfile

open (unit=out,file=outfile,status='replace')
write (*,*) ' Number of observation points:'
read(*,*) NObsPoints

! Get model type:
write(*,*) ' Model type:'
write(*,*) '   0: analytical'
write(*,*) '   1: numerical homogeneous'
write(*,*) '   2: numerical fractional'
write(*,*) '   3: numerical homogeneous (variable k/phi)'
write(*,*) '   4: numerical fractional (variable k/phi)'
write(*,*) '   5: multi-layer:'
read(*,*) ModelType

! Problem dimensions:
NPumps=1

MaxNPumpTimes=1000

select case (ModelType)
case(0)
  NVariables=2
  NFixedParameters=7
  NReservoirConditions=2
case(1)
  NVariables=2
  NFixedParameters=7
  NReservoirConditions=2
case(2)
  NVariables=3
  NFixedParameters=7
  NReservoirConditions=2
case(3)
  NVariables=2
  NFixedParameters=8
  NReservoirConditions=2
case(4)
  NVariables=4
  NFixedParameters=8
  NReservoirConditions=2
case(5)
  NVariables=4
  NFixedParameters=17
  NReservoirConditions=2
end select

allocate(Pump(NPumps),ObsPoint(NObsPoints))
allocate(ReadData(MaxTotalNData))
allocate(PumpData(NPumps,MaxNPumpTimes))
allocate(FixedParameter(NFixedParameters),ReservoirCondition(NReservoirConditions))
allocate(variable(NVariables))

Pump(1)%NData=1

! write (*,*) ' Enter model parameters:'
! read (*,*) (variable(i),i=1,NVariables)

select case (ModelType)

  case(0) ! Analytical:

    write(*,*) ' Pumping type (1:constant, 2:sinusoidal):'
    read(*,*) Pump(1)%Scheme

    select case (Pump(1)%Scheme)
	case (1)
      open(unit=dat,file='wellfit_theis.dat',status='old')
      read(dat,*) k,nu,phi,rho,c,b,Q0,P0
      variable(1)=k
      variable(2)=phi
      ReservoirCondition(1)=P0
      FixedParameter(1)=0.0_DP
      FixedParameter(2)=b
	  FixedParameter(3)=0.1
      FixedParameter(4)=nu
      FixedParameter(5)=rho
      FixedParameter(6)=c
      PumpData(1,1)%rate=Q0
	  close(dat)

	case (2)
      open(unit=dat,file='wellfit_sin.dat',status='old')
      read(dat,*) k,nu,phi,rho,c,b,Q0,P0
      read(dat,*) rw,period
      ReservoirCondition(1)=P0
      FixedParameter(1)=0.0_DP
      FixedParameter(2)=b
	  FixedParameter(3)=rw
      FixedParameter(4)=nu
      FixedParameter(5)=rho
      FixedParameter(6)=c
      PumpData(1,1)%rate=Q0
      PumpData(1,1)%time=period
	  close(dat)

	end select


  case(1,2) ! Numerical models:

    write(*,*) 'Data file:'
    read(*,'(a)') DataFileName
    open(unit=dat,file=DataFileName,status='old')
	read(dat,*) Pressure0,X0
	read(dat,*) rw,thick,CR,COND,RHOR,COMP
    ReservoirCondition(1)=Pressure0
	ReservoirCondition(2)=X0
    FixedParameter(1)=0.0_DP
	FixedParameter(2)=thick
	FixedParameter(3)=rw
	FixedParameter(4)=CR
	FixedParameter(5)=COND
	FixedParameter(6)=RHOR
    FixedParameter(7)=COMP
	read(dat,*) Pump(1)%Enthalpy, StepRates
	if (StepRates==0) then
	  Pump(1)%StepFlows=.false.
	else
	  Pump(1)%StepFlows=.true.
	end if
    read(dat,*) Pump(1)%Scheme,Pump(1)%NData
	do i=1,Pump(1)%NData
	  read(dat,*) PumpData(1,i)%time,PumpData(1,i)%rate
	end do
	close(dat)

  case(3,4) ! Numerical models with variable k/phi:

    write(*,*) 'Data file:'
    read(*,'(a)') DataFileName
    open(unit=dat,file=DataFileName,status='old')
	read(dat,*) Pressure0,X0
	read(dat,*) rw,thick,CR,COND,RHOR,COMP,COMT
    ReservoirCondition(1)=Pressure0
	ReservoirCondition(2)=X0
    FixedParameter(1)=0.0_DP
	FixedParameter(2)=thick
	FixedParameter(3)=rw
	FixedParameter(4)=CR
	FixedParameter(5)=COND
	FixedParameter(6)=RHOR
    FixedParameter(7)=COMP
    FixedParameter(8)=COMT
	read(dat,*) Pump(1)%Enthalpy, StepRates
	if (StepRates==0) then
	  Pump(1)%StepFlows=.false.
	else
	  Pump(1)%StepFlows=.true.
	end if
    read(dat,*) Pump(1)%Scheme,Pump(1)%NData
	do i=1,Pump(1)%NData
	  read(dat,*) PumpData(1,i)%time,PumpData(1,i)%rate
	end do
	close(dat)

  case(5) ! Multi-layer model:

    write(*,*) 'Data file:'
    read(*,'(a)') DataFileName
    open(unit=dat,file=DataFileName,status='old')
	read(dat,*) Pressure0,X0
	read(dat,*) rw,thick
	read(dat,*) fracth,fracsp
	read(dat,*) CRF,CONDF,RHORF,COMPF
	read(dat,*) CRM,CONDM,RHORM,COMPM
    ReservoirCondition(1)=Pressure0
	ReservoirCondition(2)=X0
    FixedParameter(1)=0.0_DP
	FixedParameter(2)=thick
	FixedParameter(3)=rw
    FixedParameter(4)=fracth
	FixedParameter(5)=fracsp
	FixedParameter(6)=CRF
	FixedParameter(7)=CRM
	FixedParameter(8)=CONDF
	FixedParameter(9)=CONDM
	FixedParameter(10)=RHORF
	FixedParameter(11)=RHORM
    FixedParameter(12)=COMPF
    FixedParameter(13)=COMPM
    FixedParameter(14)=0.0_dp
    FixedParameter(15)=0.0_dp
	FixedParameter(16)=0.0_dp
    FixedParameter(17)=0.0_dp
	read(dat,*) Pump(1)%Enthalpy, StepRates
	if (StepRates==0) then
	  Pump(1)%StepFlows=.false.
	else
	  Pump(1)%StepFlows=.true.
	end if
    read(dat,*) Pump(1)%Scheme,Pump(1)%NData
	do i=1,Pump(1)%NData
	  read(dat,*) PumpData(1,i)%time,PumpData(1,i)%rate
	end do
	close(dat)

end select

DataNo=0
TotalNData=0

do ObsPointNo=1,NObsPoints

  write (*,'(a,i2,a)') 'Observation Point ',ObsPointNo,'...'

  select case (ModelType)
  case(0) ! Analytical:
    write (*,*) 'Property (0: flow, 1:pressure):'
  case(1:5) ! Numerical:
    write (*,*) 'Property (1:pressure, 2:temperature, 3: enthalpy):'
  end select

  read (*,*) ObsPoint(ObsPointNo)%Property

!  write (*,*) 'Data offset:'
!  read (*,*) ObsPoint(ObsPointNo)%DataOffset
   ObsPoint(ObsPointNo)%DataOffset=0.0

  write (*,*) 'Distance from action well:'
  read (*,*) ObsPoint(ObsPointNo)%Position%x(1)

  ObsPoint(ObsPointNo)%Position%x(2)=0.0  ! Not used.
  ObsPoint(ObsPointNo)%Position%x(3)=0.0  ! Not used.

  write (*,*) 'Error:'
  read (*,*) ObsPoint(ObsPointNo)%Error
!   ObsPoint(ObsPointNo)%Error=0.0

  if (ObsPointNo>1) then
    ObsPoint(ObsPointNo)%DataIndex=ObsPoint(ObsPointNo-1)%DataIndex+ObsPoint(ObsPointNo-1)%NData
  else
    ObsPoint(ObsPointNo)%DataIndex=1
  end if

  write (*,*) 'Time range (start t, dt, end t) :'
  write (*,*) 'Attempting to read'
  read (*,*) t0,dt,t1
  write (*,*) 'Attempt 1'
  ObsPoint(ObsPointNo)%NData=0
  write (*,*) 'Attempt 2'
  do t=t0,t1,dt
    write (*,*) 'Attempt 3'
    DataNo=DataNo+1
    write (*,*) 'Attempt 4'
    if (DataNo<=MaxTotalNData) then
      write (*,*) 'Attempt 5'
      ReadData(DataNo)%time=t
      write (*,*) 'Attempt 6'
      ObsPoint(ObsPointNo)%NData=ObsPoint(ObsPointNo)%NData+1
    else
      write (*,*) 'Attempt 7'
      write (*,'(a,i5)') 'Too many data- can handle only',MaxTotalNData
      stop
    end if

  end do
  write (*,*) 'Attempt 8'
  TotalNData=TotalNData+ObsPoint(ObsPointNo)%NData

end do  ! ObsPoint loop.
write (*,*) 'Attempt 9'
allocate(TestData(TotalNData))
write (*,*) 'Attempt 10'
TestData=ReadData(1:TotalNData)
write (*,*) 'Attempt 11'
! Assign ObsPoint errors to TestData:
ObsPoint%Weight=1.0_dp
write (*,*) 'Attempt 12'
TestData%Weight=1.0_dp
write (*,*) 'Attempt 13'
do ObsPointNo=1,NObsPoints
  write (*,*) 'Attempt 14'
  lower=ObsPoint(ObsPointNo)%DataIndex
  write (*,*) 'Attempt 15'
  upper=lower+ObsPoint(ObsPointNo)%NData-1
  write (*,*) 'Attempt 16'
  TestData(lower:upper)%error=ObsPoint(ObsPointNo)%error
end do
write (*,*) 'Attempt 17'
! Generate modelled values:
! TestData%ModelledValue=1
TestData%ModelledValue=model(variable,updatemodelprogress)
write (*,*) 'Attempt 18'
! Add random noise:
 call AddNoise

! Write file:
write (*,*) 'Attempt 19'
write (out,'(10(e15.5,1x))') (variable(i),i=1,NVariables)
write (*,*) 'Attempt 20'
write (out,'(i15)') NObsPoints
write (*,*) 'Attempt 21'
write (out,*)
write (*,*) 'Attempt 22'
DataNo=0
write (*,*) 'Attempt 23'
do ObsPointNo=1,NObsPoints
  write (*,*) 'Attempt 24'
  write (out,'(i15,1x,i15)') ObsPointNo,ObsPoint(ObsPointNo)%Property
  write (*,*) 'Attempt 25'
  write (out,'(e15.5,1x,e15.5)') ObsPoint(ObsPointNo)%Position%x(1),ObsPoint(ObsPointNo)%DataOffset
  write (*,*) 'Attempt 26'
  write (out,'(i15,1x,e15.5)') ObsPoint(ObsPointNo)%NData,ObsPoint(ObsPointNo)%Error
  write (*,*) 'Attempt 27'
  do i=1,ObsPoint(ObsPointNo)%NData
    write (*,*) 'Attempt 28'
    DataNo=DataNo+1
    write (*,*) 'Attempt 29'
    write (out,'(e30.25,1x,e30.25)') TestData(DataNo)%time,TestData(DataNo)%ModelledValue
  end do

end do
write (*,*) 'Attempt 30'
close(out)

stop
end program PseudoData

!.....................................................................................

subroutine updatemodelprogress(nummodelruns,timestepsize, &
  progressfraction, analysisstopped)
! This is a dummy routine to pass into the fast solver.
  use variable_types
  integer(I4B), intent(in) :: nummodelruns
  real(DP), intent(in) :: timestepsize,progressfraction
  logical(LGT), intent(inout) :: analysisstopped

  return
end subroutine updatemodelprogress
