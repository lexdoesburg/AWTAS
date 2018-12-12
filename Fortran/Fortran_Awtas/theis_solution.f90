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

!   Unpack fixed parameters:
    B=FixedParameter(2)
    nu=FixedParameter(4)
    rho=FixedParameter(5)
    C=FixedParameter(6)
!   Unpack reservoir conditions:
    P0=ReservoirCondition(1)
!   Unpack variable parameters:
    k=variable(1)
    phi=variable(2)
      
!   Get pumping rate (only one constant-rate pump):
    Q0=PumpSchemeParams(1,1)
    
    CC=Q0*nu/(4.0_dp*pi_d*B*k)
    D=k/(nu*phi*rho*C)

!   Loop over ObsPoints:
    do ObsPointNo=1,NObsPoints
    
      ObsPointOffset=ObsPoint(ObsPointNo)%DataIndex-1
      r=ObsPoint(ObsPointNo)%Position%x(1)
	  DataOffset=ObsPoint(ObsPointNo)%DataOffset
      
      select case (ObsPoint(ObsPointNo)%Property)
      
        case(0)  ! Flow:
        
!         Loop over datapoints:
          do i=1,ObsPoint(ObsPointNo)%NData
      
            DataNo=ObsPointOffset+i
            t=TestData(DataNo)%time
        
            if (t >= small_D) then        
              x=0.25_dp*r**2.0_dp/(D*t)          
              flow=Q0*dexp(-x)
            else  ! t<=0
              flow=0.0_DP
            end if
            
            AnalyticalTheis(DataNo)=flow-DataOffset

          end do  ! Datapoint loop
 
        case(1)  ! Pressure:

!         Loop over datapoints:
          do i=1,ObsPoint(ObsPointNo)%NData
      
            DataNo=ObsPointOffset+i
            t=TestData(DataNo)%time
        
            if (t >= small_D) then        
              x=0.25_dp*r**2.0_dp/(D*t)          
              pressure=P0+CC*E1(x)
            else  ! t<=0
              pressure=P0
            end if
            
            AnalyticalTheis(DataNo)=pressure-DataOffset

          end do  ! Datapoint loop

      end select          
      
    end do  ! ObsPoint loop      
    
    return
    
  end function AnalyticalTheis
! -----------------------------------------------------------------------------

end module theis_solution
