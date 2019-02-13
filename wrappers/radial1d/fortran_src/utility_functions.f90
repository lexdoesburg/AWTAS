module utility_functions
! Various utility functions - Exponential Integral, Kelvin, etc

  use variable_types
  implicit none

  contains

! -----------------------------------------------------------------------------
    function poly(a,x,n)
! Evaluates the polynomial a0+a1.x+a2.x^2+a3.x^3+...+an.x^n, using nested
! multiplication.

! Argument Variables:
      integer(I4B),intent(IN)::n
      real(DP),intent(IN)::a(0:n),x
      real(DP)::poly

! Locals:
      integer(I4B)::k

! Function Unit:
      poly=a(n)
      
      k_loop: do k=n-1,0,-1
        poly=x*poly+a(k)
      end do k_loop

    return
    end function poly
! -----------------------------------------------------------------------------
    subroutine kelvin(X,XN0,PHI0,XN1,PHI1) 

! "Kelvin" routine
! Last Modified: 26 Jan 1998 by ARH
!    replaced implicit variables by real*8 declarations. Also use defns
!   of PI_D and SQRT2_D that exist in variable_types.f90
! Description: Returns scaled versions of the amplitudes
! XN0=sqrt(X)*exp(X/SQRT2_D)*N0, XN1=sqrt(X)*exp(X/SQRT2_D)*N1 

! Argument Variables:
      real(DP),intent(IN)::X
      real(DP),intent(OUT)::XN0,PHI0,XN1,PHI1
  
! Local Variables: (replaced implict declarations).
      real(DP)::XX,RX,XX2,XX3,XX4,XX6,XX8,XX10,XX12,XX14,XX16,XX18,&
        XX20,XX22,XX24,XX26,XX28,CCC,BER,BEI,XKER,XKEI,BERP,BEIP,XKERP,&
        BER1,BEI1,XKEIP,XKER1,XKEI1,RX2,RX3,RX4,RX5,RX6,RETHET,XMTHET,&
        RECC,XMCC,CC,REF,XMF,REPHI,XMPHI,PP,MM,SIGN,PPA

! Subroutine Unit:
      XX=0.125_dp*X
      RX=1.0_dp/XX

      XX_Check: if (XX <= 1.0_dp) then
        XX2=XX*XX
        XX3=XX2*XX
        XX4=XX3*XX
        XX6=XX4*XX2
        XX8=XX4*XX4
        XX10=XX8*XX2
        XX12=XX8*XX4
        XX14=XX12*XX2
        XX16=XX12*XX4
        XX18=XX16*XX2
        XX20=XX16*XX4
        XX22=XX20*XX2
        XX24=XX20*XX4
        XX26=XX24*XX2
        XX28=XX24*XX4
        CCC=1.0_dp-64.0_dp*XX4+113.77777774_dp*XX8
        CCC=CCC-32.36345652_dp*XX12+2.64191397_dp*XX16
        CCC=CCC-0.08349609_dp*XX20+0.00122552_dp*XX24
        BER=CCC-0.00000901_dp*XX28
        CCC=16.0_dp*XX2-113.77777774_dp*XX6
        CCC=CCC+72.81777742_dp*XX10-10.56765779_dp*XX14
        CCC=CCC+0.52185615_dp*XX18-0.01103667_dp*XX22
        BEI=CCC+0.00011346_dp*XX26
        CCC=-LOG(0.5_dp*X)*BER+0.25_dp*PI_D*BEI    
        CCC=CCC-0.57721566_dp-59.05819744_dp*XX4
        CCC=CCC+171.36272133_dp*XX8-60.60977451_dp*XX12
        CCC=CCC+5.65539121_dp*XX16-0.19636347_dp*XX20
        XKER=CCC+0.00309699_dp*XX24-0.00002458_dp*XX28
        CCC=-LOG(0.5_dp*X)*BEI-0.25_dp*PI_D*BER    
        CCC=CCC+6.76454936_dp*XX2-142.91827687_dp*XX6
        CCC=CCC+124.2356965_dp*XX10-21.30060904_dp*XX14
        CCC=CCC+1.17509064_dp*XX18-0.02695875_dp*XX22
        XKEI=CCC+0.00029532_dp*XX26
        CCC=-4.0_dp*XX2+14.22222222_dp*XX6
        CCC=CCC-6.06814810_dp*XX10+0.66047849_dp*XX14
        CCC=CCC-0.02609253_dp*XX18+0.00045957_dp*XX22
        CCC=CCC-0.00000394_dp*XX26
        BERP=CCC*X
        CCC=0.5_dp-10.66666666_dp*XX4+11.37777772_dp*XX8
        CCC=CCC-2.31167514_dp*XX12+0.14677204_dp*XX16
        CCC=CCC-0.00379386_dp*XX20+0.00004609_dp*XX24
        BEIP=CCC*X
        CCC=-3.69113734_dp*XX2+21.42034017_dp*XX6
        CCC=CCC-11.36433272_dp*XX10+1.41384780_dp*XX14
        CCC=CCC-0.06136358_dp*XX18+0.00116137_dp*XX22
        CCC=CCC-0.00001075_dp*XX26
        XKERP=-LOG(0.5_dp*X)*BERP-BER/X+0.25_dp*PI_D*BEIP+X*CCC   
        CCC=0.21139217_dp-13.39858846_dp*XX4
        CCC=CCC+19.41182758_dp*XX8-4.65950823_dp*XX12
        CCC=CCC+0.33049424_dp*XX16-0.00926707_dp*XX20
        CCC=CCC+0.00011997_dp*XX24
        XKEIP=-LOG(0.5_dp*X)*BEIP-BEI/X-0.25_dp*PI_D*BERP+X*CCC   
        BER1=(BERP-BEIP)/SQRT2_D     
        BEI1=(BERP+BEIP)/SQRT2_D     
        XKER1=(XKERP-XKEIP)/SQRT2_D    
        XKEI1=(XKERP+XKEIP)/SQRT2_D       
        XN0=SQRT(XKER*XKER+XKEI*XKEI)
        XN0=SQRT(X)*EXP(X/SQRT2_D)*XN0
        PHI0=ATAN(XKEI/XKER)
        XN1=SQRT(XKERP*XKERP+XKEIP*XKEIP)
        XN1=SQRT(X)*EXP(X/SQRT2_D)*XN1
        PHI1=0.25_dp*PI + ATAN(XKEIP/XKERP)-PI_D  
        if (X > 1.71854306_dp) PHI0=PHI0-PI_D
        if (X > 6.1273188_dp) PHI0=PHI0-PI_D
        if (X > 2.66583954_dp) PHI1=PHI1-PI_D 
        if (X > 7.1719675_dp) PHI1=PHI1-PI_D 
      else
        RX2=RX*RX
        RX3=RX2*RX
        RX4=RX3*RX
        RX5=RX4*RX
        RX6=RX5*RX  
        CCC=-0.0110486_dp*RX+0.0000906_dp*RX3
        RETHET=CCC-0.0000252_dp*RX4+0.0000034_dp*RX5+0.0000006_dp*RX6
        CCC=-0.3926991_dp+0.0110485_dp*RX-0.0009765_dp*RX2
        CCC=CCC+0.0000901_dp*RX3-0.0000051_dp*RX5
        XMTHET=CCC+0.0000019_dp*RX6
        RECC= RETHET
        XMCC=-X/SQRT2_D + XMTHET
        CC=PI_D/2.0_dp
        REF=sqrt(CC)*exp(RECC)*cos(XMCC)
        XMF=sqrt(CC)*exp(RECC)*sin(XMCC)
        XKER=REF
        XKEI=XMF
        CCC=0.7071068_dp+0.0625001_dp-0.0013813_dp*RX2
        CCC=CCC-0.0000005_dp*RX3
        REPHI=CCC+0.0000346_dp*RX4-0.0000117_dp*RX5+0.0000016_dp*RX6
        CCC=0.7071068_dp+0.0000001_dp*RX+0.0013811_dp*RX2
        CCC=CCC-0.0002452_dp*RX3+0.0000338_dp*RX4
        XMPHI=CCC+0.0000024_dp*RX5-0.0000032_dp*RX6
        XKERP=-(REF*REPHI-XMF*XMPHI)
        XKEIP=-(REF*XMPHI+XMF*REPHI)
        XN0=sqrt(XKER*XKER+XKEI*XKEI)
        PHI0=atan(XKEI/XKER)
        XN1=sqrt(XKERP*XKERP+XKEIP*XKEIP)
        PHI1=0.25_dp*PI_D+atan(XKEIP/XKERP)
        PP=PHI0+X/SQRT2_D+0.39720_dp
        MM=int(PP/PI_D)
        PHI0=PHI0-MM*PI_D
        PP=PHI1+X/SQRT2_D+1.96350_dp
        PPA=abs(PP)
        SIGN=PP/PPA
        MM=int(PPA/PI_D)
        PHI1=PHI1-(MM+1)*PI_D
      end if XX_Check

      return
      end subroutine kelvin
! -----------------------------------------------------------------------------
      function E1(x)
! Exponential integral function.  Used for calculating Theis solution.

! Argument Variables:
        real(DP),intent(IN)::x
        real(DP)::E1

! Local Variables:
        real(DP)::a(0:5),b(0:5)
        integer(I4B)::i

! Function Unit:
        x_check: if (x <= 1.0_dp) then
          a=(/-0.57721566_dp,0.99999193_dp,-0.24991055_dp,0.05519968_dp,&
            -0.00976004_dp,0.00107857_dp/)
          E1=-dlog(x)+poly(a,x,5)
        else
          a=(/0.2677737343_dp,8.6347608925_dp,18.0590169730_dp,8.5733287401_dp,&
            1.0_dp,0.0_dp/)
          b=(/3.9584969228_dp,21.0996530827_dp,25.6329561486_dp,9.5733223454_dp,&
            1.0_dp,0.0_dp/)
          E1=dexp(-x)/x*poly(a,x,4)/poly(b,x,4)
        end if x_check

      return
    end function E1
! -----------------------------------------------------------------------------

end module utility_functions
