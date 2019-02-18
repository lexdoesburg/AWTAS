MODULE variable_types

! =============================================================================

! Module: 
!    variable_types

! Description:
!    Based on the nrtype module listed in "NUMERICAL RECIPES in 
!   Fortran 90" (NR) (W. H. Press, S. A. Teukolsky, W. T. Vetterling, 
!   B. P. Flannery). Published 1996 by Cambridge University Press. 
!   Pages 1361-1362. 
!    This module contains definitions for a number of named parameters,
!    and some elementary derived data types. The most important are
!   thost that define the KIND types for most of the variables used.
!    I4B,I2B,I1B for integer variables; SP and DP for real/real*8; 
!    SPC and DPC for the corresponding complex cases; LGT is the default 
!    for logical types.

! Created: 
!    Aaron Richard Hochwimmer - 25 April 1998
!   ------------------------
!      Entered in the nrtype.f90 file from NR.

! Modified:
!    Aaron Richard Hochwimmer - 23 June 1998
!   ------------------------
!      Converting routines to higher precision. Real*8 DP (1.0D0 is 
!      the highest precision available under standard Digital Visual
!     FORTRAN). Followed suggestion NR pg 1362. 

! =============================================================================

! Integer Declarations:
! ---------------------
! Symbolic names for KIND types of 4-, 2-, and 1-byte integers:
  INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
  INTEGER, PARAMETER :: I2B = SELECTED_INT_KIND(4)
  INTEGER, PARAMETER :: I1B = SELECTED_INT_KIND(2)

! Floating Point Declarations:
! ----------------------------
! Symbolic names for KIND types of single- and double-precision reals:
! ARH 23/6/98: Redefine SP to be the same precision as DP.
  INTEGER, PARAMETER :: SP = KIND(1.0)
  INTEGER, PARAMETER :: DP = KIND(1.0D0)

! Complex Declarations:
! ---------------------
! Symbolic names for KIND types of single- and double-precision complex:
  INTEGER, PARAMETER :: SPC = KIND((1.0,1.0))
  INTEGER, PARAMETER :: DPC = KIND((1.0D0,1.0D0))

! Logical Declarations:
! ---------------------
! Symbolic names for KIND type of default logical:
  INTEGER, PARAMETER :: LGT = KIND(.TRUE.)

! Global Parameters:
! ------------------
! Frequently used mathematical constants (with precision to spare):
!    Single Precision: (SP Versions)
  REAL(SP), PARAMETER :: PI=3.141592653589793238462643383279502884197_sp
  REAL(SP), PARAMETER :: PIO2=1.57079632679489661923132169163975144209858_sp
  REAL(SP), PARAMETER :: TWOPI=6.283185307179586476925286766559005768394_sp
  REAL(SP), PARAMETER :: SQRT2=1.41421356237309504880168872420969807856967_sp
  REAL(SP), PARAMETER :: EULER=0.5772156649015328606065120900824024310422_sp
!    Double Precision: (DP Versions)
  REAL(DP), PARAMETER :: PI_D=3.141592653589793238462643383279502884197_dp
  REAL(DP), PARAMETER :: PIO2_D=1.57079632679489661923132169163975144209858_dp
  REAL(DP), PARAMETER :: TWOPI_D=6.283185307179586476925286766559005768394_dp
  real(DP), parameter :: SQRT2_D=1.41421356237309504880168872420969807856967_dp
  real(DP), parameter :: small_D=1.e-5_dp

! Global Derived Types:
! ---------------------
! Derived data types for sparse matrices, single and 
! double precision versions:

  TYPE sprs2_sp
    INTEGER(I4B) :: n,len
    REAL(SP), DIMENSION(:), POINTER :: val
    INTEGER(I4B), DIMENSION(:), POINTER :: irow
    INTEGER(I4B), DIMENSION(:), POINTER :: jcol
  END TYPE sprs2_sp

  TYPE sprs2_dp
    INTEGER(I4B) :: n,len
    REAL(DP), DIMENSION(:), POINTER :: val
    INTEGER(I4B), DIMENSION(:), POINTER :: irow
    INTEGER(I4B), DIMENSION(:), POINTER :: jcol
  END TYPE sprs2_dp

END MODULE variable_types
! =============================================================================
