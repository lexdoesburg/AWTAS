module sinusoidal_isothermal

  use utility_functions
  use variable_types
  implicit none

  contains

! -------------------------------------------------------------------
  function AnalyticalSinusoidal(variable)

!   Generates analytical solution for sinusoidal pumping of a single-layer
!   model.

!   FixedParameter values for this model:
!   1: layer elevation (not used)
!   2: B   (layer thickness)
!   3: rw  (action well radius)
!   4: nu  (viscosity)
!   5: rho (fluid density)
!   6: C   (compressibility)

!   Reservoir conditions for this model:
!   1: Pressure

!   Variables for this model:
!   1: k (permeability)
!   2: phi (porosity)

    use problem_data
    
!   Argument Variables:
    real(DP),intent(IN) :: variable(:)
    real(DP) :: AnalyticalSinusoidal(TotalNData)
!   Fixed parameters:
    real(DP) :: nu,rho,B,C,rw 
!   Reservoir conditions:
	real(DP) :: P0
!   Variable parameters:
    real(DP) :: k,phi         
!   Local Variables:
    real(DP) :: amplitude,phase
    real(DP) :: alpha,beta,XN0W,PHI0W,XN1W,PHI1W
    real(DP) :: XN0,PHI0,XN1,PHI1
    real(DP) :: CC,D,omega,rstar,X,factor
    real(DP) :: ForcingAmplitude,ForcingPeriod
    integer(I4B) :: ObsPointOffset,ObsPointNo,i,DataNo
    real(DP) :: r,t,pressure,DataOffset
    
!   Unpack fixed parameters:
    B=FixedParameter(2)
	rw=FixedParameter(3)
    nu=FixedParameter(4)
    rho=FixedParameter(5)
    C=FixedParameter(6)
!   Unpack reservoir conditions:
    P0=ReservoirCondition(1)
!   Unpack variable parameters:
    k=variable(1)
    phi=variable(2)
    
!   Get pumping parameters (only one pump):   
    ForcingAmplitude=PumpSchemeParams(1,1)    
    ForcingPeriod=PumpSchemeParams(1,2)
    
!   Intermediate quantities:
    CC=ForcingAmplitude/(2.0_dp*pi_D*B*k/nu)
    D=(k/nu)/(phi*rho*C)
    omega= 2.0_dp*pi_D/ForcingPeriod  
    alpha=dsqrt(omega/D)
    beta=alpha*rw

!   This calculates useful numbers (XN1W,PHI1W) at the well:
    call kelvin(beta,XN0W,PHI0W,XN1W,PHI1W)
    
!   Loop over ObsPoints:
    do ObsPointNo=1,NObsPoints
    
      ObsPointOffset=ObsPoint(ObsPointNo)%DataIndex-1
      r=ObsPoint(ObsPointNo)%Position%x(1)
	  DataOffset=ObsPoint(ObsPointNo)%DataOffset
      
      rstar=r/rw
      X=beta*rstar
      factor=dexp(-(X-beta)/sqrt2_D)
      call kelvin(X,XN0,PHI0,XN1,PHI1)
      
      select case (ObsPoint(ObsPointNo)%Property)
      
        case(0)  ! Flow:

!         Calculate flow amplitude & phase at the ObsPoint position:                
          amplitude=ForcingAmplitude*dsqrt(rstar)*factor*XN1/XN1W
          phase=PHI1W-PHI1
          
!         Loop over datapoints:
          do i=1,ObsPoint(ObsPointNo)%NData
      
            DataNo=ObsPointOffset+i
            t=TestData(DataNo)%time
        
            AnalyticalSinusoidal(DataNo)=amplitude*dsin(omega*t-phase)-DataOffset
           
          end do  ! Datapoint loop
 
        case(1)  ! Pressure:
                        
!         Calculate pressure amplitude & phase at the ObsPoint position:                
          amplitude=CC/dsqrt(rstar)*factor*XN0/(beta*XN1W)
          phase=PHI1W-PHI0-0.25_DP*pi_D
           
!         Loop over datapoints:
          do i=1,ObsPoint(ObsPointNo)%NData
      
            DataNo=ObsPointOffset+i
            t=TestData(DataNo)%time
            
            pressure=P0-amplitude*dsin(omega*t-phase)
                   
            AnalyticalSinusoidal(DataNo)=pressure-DataOffset
            
          end do  ! Datapoint loop

      end select          
      
    end do  ! ObsPoint loop      
    
    return
    
  end function AnalyticalSinusoidal

!------------------------------------------------------------------------------
  
end module sinusoidal_isothermal
