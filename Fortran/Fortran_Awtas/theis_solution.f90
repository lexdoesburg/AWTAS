module theis_solution
  use utility_functions
  use variable_types
  implicit none

  contains

! -----------------------------------------------------------------------------
  function AnalyticalTheis(variable)

!   This uses the Theis curve to calculate pressure at each (r,t) using the
!   given variables.  Analytic solution due to Mike O'Sullivan.

!   FixedParameter values for this model:
!   1: layer elevation (not used)
!   2: B   (layer thickness)
!   3: action well radius (not used)
!   4: nu (viscosity)
!   5: rho (fluid density)
!   6: C   (compressibility)

!   Reservoir conditions for this model:
!   1: Pressure

!   Variables for this model:
!   1: k (permeability)
!   2: phi (porosity)

    use problem_data

!   Argument Variables:
    real(DP),intent(in) :: variable(:)
    real(DP) :: AnalyticalTheis(TotalNData)

!   Local Variables:
!   Fixed parameters:
    real(DP) :: nu,rho,B,C
!   Reservoir conditions:
	real(DP) :: P0
!   Variable parameters:
    real(DP) :: k,phi
!   Local variables:
    real(DP) :: CC,D,x ! Intermediate quantities
    real(DP) :: r,t,pressure,flow,Q0
    integer(I4B) :: ObsPointNo,i,ObsPointOffset,DataNo
	real(DP) :: DataOffset
  write (*,*) 'Inside theis 1'
!   Unpack fixed parameters:
    B=FixedParameter(2)
    write (*,*) 'Inside theis 2'
    nu=FixedParameter(4)
    write (*,*) 'Inside theis 3'
    rho=FixedParameter(5)
    write (*,*) 'Inside theis 4'
    C=FixedParameter(6)
    write (*,*) 'Inside theis 5'
!   Unpack reservoir conditions:
    P0=ReservoirCondition(1)
    write (*,*) 'Inside theis 6'
!   Unpack variable parameters:
    k=variable(1)
    write (*,*) 'Inside theis 7'
    phi=variable(2)
    write (*,*) 'Inside theis 8'

!   Get pumping rate (only one constant-rate pump):
    ! Q0=PumpSchemeParams(1,1)
    Q0=PumpData(1,1)%rate
    write (*,*) 'Inside theis 9'

    CC=Q0*nu/(4.0_dp*pi_d*B*k)
    write (*,*) 'Inside theis 10'
    D=k/(nu*phi*rho*C)
    write (*,*) 'Inside theis 11'

!   Loop over ObsPoints:
    do ObsPointNo=1,NObsPoints
      write (*,*) 'Inside theis 12'

      ObsPointOffset=ObsPoint(ObsPointNo)%DataIndex-1
      write (*,*) 'Inside theis 13'
      r=ObsPoint(ObsPointNo)%Position%x(1)
      write (*,*) 'Inside theis 14'
	  DataOffset=ObsPoint(ObsPointNo)%DataOffset
    write (*,*) 'Inside theis 15'

      select case (ObsPoint(ObsPointNo)%Property)

        case(0)  ! Flow:
          write (*,*) 'Inside theis 16'
!         Loop over datapoints:
          do i=1,ObsPoint(ObsPointNo)%NData
            write (*,*) 'Inside theis 17'
            DataNo=ObsPointOffset+i
            write (*,*) 'Inside theis 18'
            t=TestData(DataNo)%time
            write (*,*) 'Inside theis 19'
            if (t >= small_D) then
              write (*,*) 'Inside theis 20'
              x=0.25_dp*r**2.0_dp/(D*t)
              write (*,*) 'Inside theis 21'
              flow=Q0*dexp(-x)
              write (*,*) 'Inside theis 22'
            else  ! t<=0
              write (*,*) 'Inside theis 23'
              flow=0.0_DP
              write (*,*) 'Inside theis 24'
            end if
            write (*,*) 'Inside theis 25'
            AnalyticalTheis(DataNo)=flow-DataOffset

          end do  ! Datapoint loop

        case(1)  ! Pressure:
          write (*,*) 'Inside theis 26'
!         Loop over datapoints:
          do i=1,ObsPoint(ObsPointNo)%NData
            write (*,*) 'Inside theis 27'
            DataNo=ObsPointOffset+i
            write (*,*) 'Inside theis 28'
            t=TestData(DataNo)%time
            write (*,*) 'Inside theis 29'
            if (t >= small_D) then
              write (*,*) 'Inside theis 30'
              x=0.25_dp*r**2.0_dp/(D*t)
              write (*,*) 'Inside theis 31'
              pressure=P0+CC*E1(x)
              write (*,*) 'Inside theis 32'
            else  ! t<=0
              write (*,*) 'Inside theis 33'
              pressure=P0
              write (*,*) 'Inside theis 34'
            end if
            write (*,*) 'Inside theis 35'
            AnalyticalTheis(DataNo)=pressure-DataOffset
            write (*,*) 'Inside theis 36'

          end do  ! Datapoint loop

      end select

    end do  ! ObsPoint loop

    return

  end function AnalyticalTheis
! -----------------------------------------------------------------------------

end module theis_solution
