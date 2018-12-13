module noise

contains

!-----------------------------------------------------------------------

  subroutine AddNoise

    ! Adds noise to the TestData%ModelledValue results, according to the 
    ! size of the Tool%Error array.

    use variable_types
    use problem_data
    use variable_parameters

    implicit none
    ! Locals:
    integer(I4B) :: i,DataNo,StartInt,idum
    real(DP) :: gas

    write (*,*) ' Enter a random starting integer:'
    read (*,*) StartInt

    ! Generate some initial numbers, to avoid repeats:
    idum=StartInt
    do i=1,StartInt
      gas=gasdev(idum)
    end do   

    do DataNo=1,TotalNData
      TestData(DataNo)%ModelledValue=TestData(DataNo)%ModelledValue+ &
      gasdev(idum)*TestData(DataNo)%Error
    end do

    return

  end subroutine AddNoise

!--------------------------------------------------------------------------

    function gasdev(idum)
! Gaussian random number generator from Numerical Recipes (2nd ed.)

      use variable_types

      implicit none
      INTEGER(I4B),INTENT(INOUT)::idum
      REAL(DP)::gasdev
      
! Local Variables:
      INTEGER(I4B),SAVE::iset
      REAL(DP),SAVE::gset
      REAL(DP)::fac,rsq,v1,v2

! Function Unit:
      DATA iset /0/
      if (iset == 0) then
! ARH: Init rsq
        rsq=1.0_dp
        do while ((rsq >=1.0_dp) .or. (rsq == 0.0_dp))
          v1=2.*ran1(idum)-1.
          v2=2.*ran1(idum)-1.
          rsq=v1**2+v2**2
        end do
        fac=sqrt(-2.0_dp*log(rsq)/rsq)
        gset=v1*fac
        gasdev=v2*fac
        iset=1
      else
        gasdev=gset
        iset=0
      end if

      return
    end function gasdev
! -------------------------------------------------------------------

    function ran1(IDUM)
! Random number generator from Numerical Recipes (2nd ed.)

      use variable_types

      implicit none
      INTEGER(I4B),INTENT(INOUT)::idum
      REAL(DP)::ran1

! Local Variables:
      INTEGER(I4B),PARAMETER::IA=16807, IM=2147483647, IQ=127773,&
        IR=2836,NTAB=32,NDIV=1+(IM-1)/NTAB
      INTEGER(I4B)::j,k
      INTEGER,SAVE::iv(NTAB),iy
      REAL(DP),PARAMETER::AM=1.0_dp/IM,EPS=1.2e-7_dp,RNMX=1.0_dp-EPS
      
! Function Unit:
      DATA iv /NTAB*0/, iy /0/

      if ((idum <= 0) .or. (iy == 0)) then
        idum=max(-idum,1)
        do j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
        end do
        iy=iv(1)
      end if
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)

      return
    end function ran1

! -------------------------------------------------------------------

end module noise