module models
! Routines related to well testing models...

  use problem_data
  use variable_parameters
  use variable_types
  use theis_solution
  ! use sinusoidal_isothermal
  use HomogeneousPorousSimulator
  ! use FractionalSimulator
  ! use VariableKPhiSimulator
  ! use FractionalVariableKPhiSimulator
  ! use MultiLayerSimulator
  ! use SkinSimulator
  ! use FractionalSkinSimulator
  ! use VariableKPhiSkinSimulator
  ! use LeakySimulator
  ! use MINCSimulator
  ! use WellboreStorageSimulator
  ! use TwoFeedzoneSimulator

  implicit none

  contains

! -----------------------------------------------------------------------------
    function model(variable,updatemodelprogress)
! Calculates model response, given the parameter values, and returns
! the values in the array model_result.  Values are calculated for all
! radii and times in the 'TestData' array.

! Argument Variables:
      implicit none
      real(DP),intent(in)::variable(NVariables)
      real(DP)::model(TotalNData)
	  external updatemodelprogress

      select case(ModelType)
      case (0) ! Homogeneous porous layer- analytical:
        write (*,*) 'Inside model 1'
        select case(Pump(1)%Scheme)
        case (1)
          write (*,*) 'Inside model 2'
          model=AnalyticalTheis(variable)
          write (*,*) 'Inside model 3'
        ! case (2)
        !   model=AnalyticalSinusoidal(variable)
        end select
      case (1)  ! Homogeneous porous layer- numerical:
          write (*,*) 'Calling homogeneousporous from model'
		      model=HomogeneousPorous(variable,updatemodelprogress)
	  ! case (2)  ! Numerical fractional dimension model:
		! model=Fractional(variable,updatemodelprogress)
	  ! case (3)  ! Homog. porous model with variable k/phi:
		! model=VariableKPhi(variable,updatemodelprogress)
	  ! case (4)  ! Fractional dimension model with variable k/phi:
		! model=FractionalVariableKPhi(variable,updatemodelprogress)
    !   case(5)   ! Multi-layer model:
		!  model=MultiLayer(variable,updatemodelprogress)
	  ! case(6)   ! Skin model:
		!  model=Skin(variable,updatemodelprogress)
	  ! case(7)   ! Fractional dimension skin model:
		!  model=FractionalSkin(variable,updatemodelprogress)
	  ! case(8)   ! Variable k/phi with skin:
	  !    model=VariableKPhiSkin(variable,updatemodelprogress)
	  ! case(9)   ! Leaky aquifer model:
	  !    model=Leaky(variable,updatemodelprogress)
	  ! case(10)   ! MINC model:
	  !    model=MINC(variable,updatemodelprogress)
	  ! case(11)   ! Wellbore storage model:
	  !    model=WellboreStorage(variable,updatemodelprogress)
      ! case(12)   ! Two-feedzone model:
	    !  model=TwoFeedzone(variable,updatemodelprogress)
      end select

      return
    end function model
! -----------------------------------------------------------------------------

  subroutine SetWellBlockIncl
    implicit none
	select case(ModelType)
	case default
	  WellBlockIncl=.false.
	case(11)
	  WellBlockIncl=.true.
	case(12)
	  WellBlockIncl=.true.
	end select
    return
  end subroutine SetWellBlockIncl

end module models
