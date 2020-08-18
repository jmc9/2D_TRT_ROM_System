MODULE TEMP_FUNCTIONS
    USE INTERPOLATORS
    IMPLICIT NONE

CONTAINS

!============================================================================================================!
!
!============================================================================================================!

FUNCTION kapB_calc(Temp,nu_low,nu_high)

  IMPLICIT NONE
  REAL*8,INTENT(IN):: Temp, nu_low, nu_high

  REAL*8:: kapB_calc, sig, z1, z2
  REAL*8,PARAMETER:: sig0=27d0*1d9

  z1=nu_high/Temp
  z2=nu_low/Temp

  IF (z2 .GT. 10d0) THEN
    sig=-EXP(-(z1-z2))*(z1**3 + 3d0*z1**2 + 6d0*z1 + 7.28d0) + (z2**3 + 3d0*z2**2 + 6d0*z2 + 7.28d0)
    kapB_calc=sig0*(1d0-EXP(-(nu_high-nu_low)/Temp))/(Temp**3*sig)

  ELSE
    sig = sigint(z1) - sigint(z2)
    kapB_calc=sig0*(EXP(-nu_low/Temp) - EXP(-nu_high/Temp))/(Temp**3*sig)

  END IF

END FUNCTION

!============================================================================================================!
!
!============================================================================================================!

FUNCTION kapB_dT_calc(Temp,nu_low,nu_high)

  IMPLICIT NONE
  REAL*8,INTENT(IN):: Temp, nu_low, nu_high
  REAL*8::kapB_dT_calc, sig, z1, z2
  REAL*8,PARAMETER:: sig0=27d0*1d9

  z1=nu_high/Temp
  z2=nu_low/Temp

  IF (nu_low .LE. 1d-15) THEN !nu_g(1)=0

    IF (z2 .GT. 10d0) THEN
      sig=-EXP(-(z1-z2))*(z1**3 + 3d0*z1**2 + 6d0*z1 + 7.28d0) + (z2**3 + 3d0*z2**2 + 6d0*z2 + 7.28d0)
      kapB_dT_calc = (kapB_calc(Temp,nu_low,nu_high)/Temp)*&
                     ( z2 - 3d0 + ( ( (nu_high**4*EXP(-(nu_high-nu_low)/Temp))/(1d0 - EXP(-z1)) ) )&
                     /( Temp**4*sig ) ) -&
                     (sig0*(nu_high - nu_low)*EXP(-(nu_high-nu_low)/Temp))/(Temp**5*sig)

    ELSE
      sig = sigint(z1) - sigint(z2)
      kapB_dT_calc = (kapB_calc(Temp,nu_low,nu_high)/Temp)*( z2 - 3d0 + ( ( (nu_high**4*EXP(-z1))/(1d0 - EXP(-z1)) ) )&
                     /( Temp**4*sig ) ) -&
                     (sig0*(nu_high - nu_low)*EXP(-z1))/(Temp**5*sig)

    END IF

  ELSE

    IF (z2 .GT. 10d0) THEN
      sig=-EXP(-(z1-z2))*(z1**3 + 3d0*z1**2 + 6d0*z1 + 7.28d0) + (z2**3 + 3d0*z2**2 + 6d0*z2 + 7.28d0)
      kapB_dT_calc = (kapB_calc(Temp,nu_low,nu_high)/Temp)*( z2 - 3d0 +&
                     ( ( (nu_high**4*EXP(-(nu_high-nu_low)/Temp))/(1d0 - EXP(-z1)) ) -&
                     ( (nu_low**4)/(1d0 - EXP(-z2)) ) )/( Temp**4*sig ) ) -&
                     (sig0*(nu_high - nu_low)*EXP(-(nu_high-nu_low)/Temp))/(Temp**5*sig)

    ELSE
      sig = sigint(z1) - sigint(z2)
      kapB_dT_calc = (kapB_calc(Temp,nu_low,nu_high)/Temp)*( z2 - 3d0 + ( ( (nu_high**4*EXP(-z1))/(1d0 - EXP(-z1)) ) -&
                     ( (nu_low**4*EXP(-z2))/(1d0 - EXP(-z2)) ) )/( Temp**4*sig ) ) -&
                     (sig0*(nu_high - nu_low)*EXP(-z1))/(Temp**5*sig)

    END IF

  END IF

END FUNCTION

!============================================================================================================!
!
!============================================================================================================!

FUNCTION Bg_planck_calc(Tin,nu_low,nu_high,comp_unit)

  IMPLICIT NONE
  REAL*8::Bg_planck_calc
  REAL*8::Tin, sig, nu_low, nu_high, comp_unit
  REAL*8,PARAMETER:: c=2.99792458d2
  REAL*8,PARAMETER:: erg=6.24150934d11
  REAL*8,PARAMETER:: h=6.62613d-19

  IF (Tin .LE. 1d-15) THEN !essentially Tin=0
    Bg_planck_calc = 0d0
  ELSE
    sig = (sigint(nu_high/Tin) - sigint(nu_low/Tin))
    sig = sig/(h**3*c**2*erg**4*comp_unit)

    Bg_planck_calc = 2d0*Tin**4*sig
  END IF

END FUNCTION

!============================================================================================================!
!
!============================================================================================================!

FUNCTION sigint(z)

  IMPLICIT NONE
  REAL*8:: sigint, q, theta
  REAL*8,INTENT(IN):: z
  REAL*8,PARAMETER:: pi=4d0*ATAN(1d0)

  ! IF (z .LE. 2d0) THEN
  !   sigint = z**3*(1d0/3d0 - z/8d0 + z**2/62.4d0)
  ! ELSE
  !   sigint = pi**4/15d0 - EXP(-z)*(z**3 + 3d0*z**2 + 6d0*z + 7.28d0)
  ! END IF

  IF (z .LE. 1.95d0) THEN
    sigint = z**3*(1d0/3d0 - z/8d0 + z**2/62.4d0)
  ELSE IF (z .GT. 2.05d0) THEN
    sigint = pi**4/15d0 - EXP(-z)*(z**3 + 3d0*z**2 + 6d0*z + 7.28d0)
  ELSE
    sigint = 1.95d0**3*(1d0/3d0 - 1.95d0/8d0 + 1.95d0**2/62.4d0)
    q = pi**4/15d0 - EXP(-2.05d0)*(2.05d0**3 + 3d0*2.05d0**2 + 6d0*2.05d0 + 7.28d0)
    CALL Interpolator(sigint,q,1.95d0,2.05d0,z)
  END IF

  ! theta = 1d0/(1d0-EXP(-(z-2d0)/(1d-3)))
  ! sigint = (1d0-theta)*(z**3*(1d0/3d0 - z/8d0 + z**2/64d0)) +&
  !   theta*(pi**4/15d0 - EXP(-z)*(z**3 + 3d0*z**2 + 6d0*z + 7.28d0))

END FUNCTION

!============================================================================================================!
!
!============================================================================================================!

SUBROUTINE kapR_hat_calc_all(kapR_hat,kapE,delx)

  IMPLICIT NONE
  REAL*8,INTENT(IN):: kapE(:,:), delx(:)
  REAL*8,INTENT(OUT):: kapR_hat(:,:)
  INTEGER:: cells, groups
  INTEGER:: i, g

  cells = SIZE(kapE,1)
  groups = SIZE(kapE,2)

  kapR_hat(1,:) = kapE(1,:)
  kapR_hat(cells+1,:) = kapE(cells,:)
  DO g=1,groups
    DO i=2,cells
      kapR_hat(i,g) = kapR_hat_calc(kapE(i-1,g),kapE(i,g),delx(i-1),delx(i))
    END DO
  END DO

END SUBROUTINE

FUNCTION kapR_hat_calc(kapE_1,kapE_2,delx_1,delx_2)

  IMPLICIT NONE
  REAL*8,INTENT(IN):: kapE_1, kapE_2, delx_1, delx_2
  REAL*8:: kapR_hat_calc

  kapR_hat_calc=(kapE_1*delx_1 + kapE_2*delx_2)/(delx_1 + delx_2)

END FUNCTION

!============================================================================================================!
!
!============================================================================================================!

FUNCTION GREY_kapE_dT_calc(kapE_dT_mod,kapE_dT_flag,enrgy_strc,nu_g,MG_E,Temp,Temp_mold,GREY_kapE,GREY_kapE_mold,&
  GREY_kapE_dT,MG_it)

  CHARACTER(100),INTENT(IN):: kapE_dT_mod, kapE_dT_flag
  CHARACTER(100),INTENT(IN):: enrgy_strc
  REAL*8,INTENT(IN):: nu_g(:), MG_E(:), Temp, Temp_mold, GREY_kapE, GREY_kapE_mold, GREY_kapE_dT
  INTEGER,INTENT(IN):: MG_it
  REAL*8:: GREY_kapE_dT_calc
  REAL*8:: sum1, sum2
  REAL*8,ALLOCATABLE:: kapB_dT(:)
  INTEGER:: groups, g

  IF (kapE_dT_mod .EQ. 'exact') THEN

    IF (kapE_dT_flag .EQ. 'on') THEN

      groups = SIZE(nu_g,1)-1
      ALLOCATE(kapB_dT(groups))
      DO g=1,groups
        IF (enrgy_strc .EQ. 'ONE_GROUP') THEN
          kapB_dT(g) = kapB_dT_calc(Temp,nu_g(1),nu_g(2))
        ELSE
          kapB_dT(g) = kapB_dT_calc(Temp,nu_g(g),nu_g(g+1))
        END IF
      END DO

      sum1=0d0
      sum2=0d0
      DO g=1,groups
        sum1 = sum1 + kapB_dT(g)*MG_E(g)
        sum2 = sum2 + MG_E(g)
      END DO
      GREY_kapE_dT_calc = sum1/sum2

      DEALLOCATE(kapb_dT)

    ELSE IF (kapE_dT_flag .EQ. 'off') THEN
      GREY_kapE_dT_calc = 0d0

    ELSE
      STOP 'unknown value for "kapE_dT_flag"'

    END IF

  ELSE

    IF (kapE_dT_flag .EQ. 'on') THEN
      IF (MG_it .LE. 1) THEN
        GREY_kapE_dT_calc = 0d0
      ELSE
        IF (ABS(Temp-Temp_mold) .LT. 1d-15) THEN
          GREY_kapE_dT_calc = GREY_kapE_dT
        ELSE
          GREY_kapE_dT_calc = ( GREY_kapE - GREY_kapE_mold )/( Temp - Temp_mold )
        END IF
      END IF

    ELSE IF (kapE_dT_flag .EQ. 'off') THEN
      GREY_kapE_dT_calc = 0d0

    ELSE
      STOP 'unknown value for "kapE_dT_flag"'

    END IF

  END IF

END FUNCTION

!============================================================================================================!
!
!============================================================================================================!


END MODULE TEMP_FUNCTIONS
