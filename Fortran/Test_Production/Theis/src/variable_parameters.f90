module variable_parameters

  use variable_types
  implicit none
! Module containing variable bounds & scaling parameters.

  integer(I4B)             :: NVariables
  real(DP),allocatable     :: scale_factor(:)
  real(DP),allocatable     :: lowerbounds(:)
  real(DP),allocatable     :: upperbounds(:)
  logical(LGT),allocatable :: logscale(:)
  real(DP)                 :: ObjectiveScale

  contains

! -----------------------------------------------------------------------------

  function scale_variable(variable)
! "Scale Variable" transforms physical parameters to scaled variables xi. Those
! for which logscale(i) is true are scaled logarithmically.

! Argument variables:
    real(DP),intent(IN)::variable(NVariables)
    real(DP)::scale_variable(NVariables)

! Function Unit:
    where(logscale)
      scale_variable=dlog(variable/scale_factor)
    elsewhere
      scale_variable=variable/scale_factor
    end where

    return
  end function scale_variable
  
! -----------------------------------------------------------------------------

  function unscale_variable(xi)
! "Unscale Variable" reverses the scale transformation.

! Argument Variables:
    real(DP),intent(IN)::xi(NVariables)
    real(DP)::unscale_variable(NVariables)

! Function Unit:
    where (logscale)
      unscale_variable=scale_factor*dexp(xi)
    elsewhere
      unscale_variable=xi*scale_factor
    end where

    return
  end function unscale_variable

!----------------------------------------------------------------------------------

  subroutine DestroyVariableParameterArrays
    deallocate(scale_factor,lowerbounds,upperbounds,logscale)
	return
  end subroutine DestroyVariableParameterArrays

!----------------------------------------------------------------------------------

end module variable_parameters
