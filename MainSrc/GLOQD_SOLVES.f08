MODULE GLOQD_SOLVES

  USE QD_SOLVERS

  IMPLICIT NONE

CONTAINS

!==================================================================================================================================!
!
!==================================================================================================================================!
SUBROUTINE OLD_GREY_COEFS(c, Delt, Theta,Gold,Temp_old,KapE_bar_old,E_avg_old,GQD_Src_old,Gold_Hat,Rhat_old)

  REAL*8,INTENT(IN):: c, Delt, Theta
  REAL*8,INTENT(IN):: Gold(:,:,:), Temp_old(:,:), KapE_bar_old(:,:), E_avg_old(:,:), GQD_Src_old(:,:)
  REAL*8,INTENT(OUT):: Gold_Hat(:,:), Rhat_old(:,:)

  INTEGER:: N_x, N_y, N_g, i, j, g

  N_x = SIZE(Gold,1)
  N_y = SIZE(Gold,2)
  N_g = SIZE(Gold,3)

  Gold_Hat = 0d0
  DO g=1,N_g
    DO j=1,N_y
      DO i=1,N_x
        Gold_Hat(i,j) = Gold_Hat(i,j) + Gold(i,j,g)
      END DO
    END DO
  END DO

  Rhat_old = 0d0
  DO j=1,N_y
    DO i=1,N_x
      Gold_Hat(i,j) = Gold_Hat(i,j) + A(i,j)*E_avg_old(i,j)/(Theta*Delt)
      Rhat_old(i,j) = cv/(Theta*Delt)*Temp_old(i,j) + ((1d0-Theta)/Theta)*(c*KapE_bar_old(i,j)*E_avg_old(i,j) - GQD_Src_old(i,j))
    END DO
  END DO


END SUBROUTINE OLD_GREY_COEFS

!==================================================================================================================================!
!
!==================================================================================================================================!
SUBROUTINE GREY_COEFS(Eg_avg,Eg_edgV,Eg_edgH,Fxg_edgV,Fyg_edgH,fg_avg_xx,fg_avg_yy,fg_edgV_xx,fg_edgV_xy,fg_edgH_yy,&
  fg_edgH_xy,KapE,KapR,Delx,Dely,A,c,Delt,Theta,KapE_Bar,DC_xx,DL_xx,DR_xx,DC_yy,DB_yy,DT_yy,DL_xy,DR_xy,DB_xy,DT_xy,&
  PL,PR,PB,PT)

  REAL*8,INTENT(IN):: Eg_avg(:,:,:), Eg_edgV(:,:,:), Eg_edgH(:,:,:)
  REAL*8,INTENT(IN):: Fxg_edgV(:,:,:), Fyg_edgH(:,:,:)
  REAL*8,INTENT(IN):: fg_avg_xx(:,:,:), fg_avg_yy(:,:,:)
  REAL*8,INTENT(IN):: fg_edgV_xx(:,:,:), fg_edgV_xy(:,:,:)
  REAL*8,INTENT(IN):: fg_edgH_yy(:,:,:), fg_edgH_xy(:,:,:)
  REAL*8,INTENT(IN):: KapE(:,:,:), KapR(:,:,:)
  REAL*8,INTENT(IN):: Delx(:), Dely(:), A(:,:)
  REAL*8,INTENT(IN):: c, Delt, Theta

  REAL*8,INTENT(OUT):: KapE_Bar(:,:)
  REAL*8,INTENT(OUT):: DC_xx(:,:), DL_xx(:,:), DR_xx(:,:)
  REAL*8,INTENT(OUT):: DC_yy(:,:), DB_yy(:,:), DT_yy(:,:)
  REAL*8,INTENT(OUT):: DL_xy(:,:), DR_xy(:,:), DB_xy(:,:), DT_xy(:,:)
  REAL*8,INTENT(OUT):: PL(:,:), PR(:,:), PB(:,:), PT(:,:)

  REAL*8:: Tilde_KapR, sum1, sum2, sum3, sum4, sum5
  INTEGER:: N_x, N_y, N_g, i, j, g

  N_x = SIZE(Eg_avg,1)
  N_y = SIZE(Eg_avg,2)
  N_g = SIZE(Eg_avg,3)

  KapE_Bar = 0d0
  DC_xx = 0d0
  DL_xx = 0d0
  DR_xx = 0d0
  DC_yy = 0d0
  DB_yy = 0d0
  DT_yy = 0d0
  DL_xy = 0d0
  DR_xy = 0d0
  DB_xy = 0d0
  DT_xy = 0d0
  PL = 0d0
  PR = 0d0
  PB = 0d0
  PT = 0d0
  DO g=1,N_g
    DO j=1,N_y
      DO i=1,N_x
        Tilde_KapR = KapR(i,j,g) + 1d0/(Theta*c*Delt)
        KapE_Bar(i,j) = KapE_Bar(i,j) + KapE(i,j,g)*Eg_avg(i,j,g)
        DC_xx(i,j) = DC_xx(i,j) + fg_avg_xx(i,j,g)*Eg_avg(i,j,g)/Tilde_KapR

        !D parameters
        DL_xx(i,j) = DL_xx(i,j) + fg_edgV_xx(i,j,g)*Eg_edgV(i,j,g)/Tilde_KapR
        DR_xx(i,j) = DR_xx(i,j) + fg_edgV_xx(i+1,j,g)*Eg_edgV(i+1,j,g)/Tilde_KapR
        DC_yy(i,j) = DC_yy(i,j) + fg_avg_yy(i,j,g)*Eg_avg(i,j,g)/Tilde_KapR
        DB_yy(i,j) = DB_yy(i,j) + fg_edgH_yy(i,j,g)*Eg_edgH(i,j,g)/Tilde_KapR
        DT_yy(i,j) = DT_yy(i,j) + fg_edgH_yy(i,j+1,g)*Eg_edgH(i,j+1,g)/Tilde_KapR
        DL_xy(i,j) = DL_xy(i,j) + fg_edgV_xy(i,j,g)*Eg_edgV(i,j,g)/Tilde_KapR
        DR_xy(i,j) = DR_xy(i,j) + fg_edgV_xy(i+1,j,g)*Eg_edgV(i+1,j,g)/Tilde_KapR
        DB_xy(i,j) = DB_xy(i,j) + fg_edgH_xy(i,j,g)*Eg_edgH(i,j,g)/Tilde_KapR
        DT_xy(i,j) = DT_xy(i,j) + fg_edgH_xy(i,j+1,g)*Eg_edgH(i,j+1,g)/Tilde_KapR

        !P parameters
        PL(i,j) = PL(i,j) + (A(i,j)/(2d0*Theta*c*Delt)*Fxg_edgV(i,j,g) + Pold_L(i,j,g))/Tilde_KapR
        PR(i,j) = PR(i,j) + (A(i,j)/(2d0*Theta*c*Delt)*Fxg_edgV(i+1,j,g) + Pold_R(i,j,g))/Tilde_KapR
        PB(i,j) = PB(i,j) + (A(i,j)/(2d0*Theta*c*Delt)*Fyg_edgH(i,j,g) + Pold_B(i,j,g))/Tilde_KapR
        PT(i,j) = PT(i,j) + (A(i,j)/(2d0*Theta*c*Delt)*Fyg_edgH(i,j+1,g) + Pold_T(i,j,g))/Tilde_KapR
      END DO
    END DO
  END DO

  DO j=1,N_y
    DO i=1,N_x
      !precalculating weighting sums
      sum1 = SUM(Eg_avg(i,j,:))
      sum2 = SUM(Eg_edgV(i,j,:))
      sum3 = SUM(Eg_edgV(i+1,j,:))
      sum4 = SUM(Eg_edgH(i,j,:))
      sum5 = SUM(Eg_edgH(i,j+1,:))

      !dividing weighted coefficients by respective weights
      KapE_Bar(i,j) = KapE_Bar(i,j)/sum1
      DC_xx(i,j) = DC_xx(i,j)/sum1
      DL_xx(i,j) = DL_xx(i,j)/sum2
      DR_xx(i,j) = DR_xx(i,j)/sum3
      DC_yy(i,j) = DC_yy(i,j)/sum1
      DB_yy(i,j) = DB_yy(i,j)/sum4
      DT_yy(i,j) = DT_yy(i,j)/sum5
      DL_xy(i,j) = DL_xy(i,j)/sum2
      DR_xy(i,j) = DR_xy(i,j)/sum3
      DB_xy(i,j) = DB_xy(i,j)/sum4
      DT_xy(i,j) = DT_xy(i,j)/sum5
    END DO
  END DO


END SUBROUTINE GREY_COEFS

!==================================================================================================================================!
!Subroutine Cbar_Calc
!
! Calculates all C_bar terms for grey QD boundary conditions
!
! NOTES::
!
! WARNINGS::
!
! OUTPUTS::
!   Cb_L - "C bar" coefficient for the left-boundary GQD condition; dim(N_y)
!   Cb_B - "C bar" coefficient for the bottom-boundary GQD condition; dim(N_x)
!   Cb_R - "C bar" coefficient for the right-boundary GQD condition; dim(N_y)
!   Cb_T - "C bar" coefficient for the top-boundary GQD condition; dim(N_x)
!
! INPUTS::
!   Cb_L - "Cg" coefficient for the left-boundary MGQD condition; dim(N_y,N_g)
!   Cb_B - "Cg" coefficient for the bottom-boundary MGQD condition; dim(N_x,N_g)
!   Cb_R - "Cg" coefficient for the right-boundary MGQD condition; dim(N_y,N_g)
!   Cb_T - "Cg" coefficient for the top-boundary MGQD condition; dim(N_x,N_g)
!   Eg_avg - cell-averaged multigroup radiation energy densities; dim(N_x,N_y,N_g)
!   Eg_edgV - cell-edge multigroup radiation energy densities on 'verticle' cell edges (x = constant); dim(N_x+1,N_y,N_g)
!   Eg_edgH - cell-edge multigroup radiation energy densities on 'horizontal' cell edges (y = constant); dim(N_x,N_y+1,N_g)
!   Eg_in_L - multigroup energy density of radiation moving into the domain through the left boundary; dim(N_y,N_g)
!   Eg_in_B - multigroup energy density of radiation moving into the domain through the bottom boundary; dim(N_x,N_g)
!   Eg_in_R - multigroup energy density of radiation moving into the domain through the right boundary; dim(N_y,N_g)
!   Eg_in_T - multigroup energy density of radiation moving into the domain through the top boundary; dim(N_x,N_g)
!
!==================================================================================================================================!
SUBROUTINE Cbar_Calc(Cb_L,Cb_B,Cb_R,Cb_T,Cg_L,Cg_B,Cg_R,Cg_T,Eg_avg,Eg_edgV,Eg_edgH,Eg_in_L,Eg_in_B,Eg_in_R,Eg_in_T)

  REAL*8,INTENT(OUT):: Cb_L(:), Cb_B(:), Cb_R(:), Cb_T(:)
  REAL*8,INTENT(IN):: Cg_L(:,:), Cg_B(:,:), Cg_R(:,:), Cg_T(:,:)
  REAL*8,INTENT(IN):: Eg_avg(:,:,:), Eg_edgV(:,:,:), Eg_edgH(:,:,:)
  REAL*8,INTENT(IN):: Eg_in_L(:,:), Eg_in_B(:,:), Eg_in_R(:,:), Eg_in_T(:,:)

  REAL*8:: sum1, sum2, q
  INTEGER:: N_x, N_y, N_g, i, j, g

  N_x = SIZE(Eg_avg,1)
  N_y = SIZE(Eg_avg,2)
  N_g = SIZE(Eg_avg,3)

  Cb_L = 0d0
  Cb_R = 0d0
  DO j=1,N_y
    sum1 = 0d0
    sum2 = 0d0
    DO g=1,N_g
      !Left-boundary grey condition
      q = Eg_edgV(1,j,g) - Eg_in_L(j,g)
      Cb_L(j) = Cb_L(j) + Cg_L(j,g)*q
      sum1 = sum1 + q
      !Right-boundary grey condition
      q = Eg_edgV(N_x+1,j,g) - Eg_in_R(j,g)
      Cb_R(j) = Cb_R(j) + Cg_R(j,g)*q
      sum2 = sum2 + q
    END DO
    !dividing by respective weights
    Cb_L(j) = Cb_L(j)/sum1
    Cb_R(j) = Cb_R(j)/sum2
  END DO

  Cb_B = 0d0
  Cb_T = 0d0
  DO i=1,N_x
    sum1 = 0d0
    sum2 = 0d0
    DO g=1,N_g
      !Bottom-boundary grey condition
      q = Eg_edgH(i,1,g) - Eg_in_B(i,g)
      Cb_B(i) = Cb_B(i) + Cg_B(i,g)*q
      sum1 = sum1 + q
      !Top-boundary grey condition
      q = Eg_edgH(i,N_y+1,g) - Eg_in_T(i,g)
      Cb_T(i) = Cb_T(i) + Cg_R(i,g)*q
      sum2 = sum2 + q
    END DO
    !dividing by respective weights
    Cb_B(i) = Cb_B(i)/sum1
    Cb_T(i) = Cb_T(i)/sum2
  END DO

END SUBROUTINE Cbar_Calc

!==================================================================================================================================!
!Subroutine GQD_In_Calc
!
! Calculates the total energy density of radiation moving into the domain through each boundary
!
! NOTES::
!
! WARNINGS::
!
! OUTPUTS::
!   E_in_L - total energy density of radiation moving into the domain through the left boundary; dim(N_y)
!   E_in_B - total energy density of radiation moving into the domain through the bottom boundary; dim(N_x)
!   E_in_R - total energy density of radiation moving into the domain through the right boundary; dim(N_y)
!   E_in_T - total energy density of radiation moving into the domain through the top boundary; dim(N_x)
!   F_in_L - total flux of radiation moving into the domain through the left boundary; dim(N_y)
!   F_in_B - total flux of radiation moving into the domain through the bottom boundary; dim(N_x)
!   F_in_R - total flux of radiation moving into the domain through the right boundary; dim(N_y)
!   F_in_T - total flux of radiation moving into the domain through the top boundary; dim(N_x)
!
! INPUTS::
!   Eg_in_L - multigroup energy density of radiation moving into the domain through the left boundary; dim(N_y,N_g)
!   Eg_in_B - multigroup energy density of radiation moving into the domain through the bottom boundary; dim(N_x,N_g)
!   Eg_in_R - multigroup energy density of radiation moving into the domain through the right boundary; dim(N_y,N_g)
!   Eg_in_T - multigroup energy density of radiation moving into the domain through the top boundary; dim(N_x,N_g)
!   Fg_in_L - multigroup flux of radiation moving into the domain through the left boundary; dim(N_y,N_g)
!   Fg_in_B - multigroup flux of radiation moving into the domain through the bottom boundary; dim(N_x,N_g)
!   Fg_in_R - multigroup flux of radiation moving into the domain through the right boundary; dim(N_y,N_g)
!   Fg_in_T - multigroup flux of radiation moving into the domain through the top boundary; dim(N_x,N_g)
!
!==================================================================================================================================!
SUBROUTINE GQD_In_Calc(E_in_L,E_in_B,E_in_R,E_in_T,F_in_L,F_in_B,F_in_R,F_in_T,Eg_in_L,Eg_in_B,Eg_in_R,Eg_in_T,&
  Fg_in_L,Fg_in_B,Fg_in_R,Fg_in_T)

  REAL*8,INTENT(OUT):: E_in_L(:), E_in_B(:), E_in_R(:), E_in_T(:)
  REAL*8,INTENT(OUT):: F_in_L(:), F_in_B(:), F_in_R(:), F_in_T(:)
  REAL*8,INTENT(IN):: Eg_in_L(:,:), Eg_in_B(:,:), Eg_in_R(:,:), Eg_in_T(:,:)
  REAL*8,INTENT(IN):: Fg_in_L(:,:), Fg_in_B(:,:), Fg_in_R(:,:), Fg_in_T(:,:)

  INTEGER:: N_x, N_y, i

  N_x = SIZE(E_in_B,1)
  N_y = SIZE(E_in_L,1)

  DO i=1,N_y
    E_in_L(i) = SUM(Eg_in_L(i,:))
    E_in_R(i) = SUM(Eg_in_R(i,:))
    F_in_L(i) = SUM(Fg_in_L(i,:))
    F_in_R(i) = SUM(Fg_in_R(i,:))
  END DO

  DO i=1,N_x
    E_in_B(i) = SUM(Eg_in_B(i,:))
    E_in_T(i) = SUM(Eg_in_T(i,:))
    F_in_B(i) = SUM(Fg_in_B(i,:))
    F_in_T(i) = SUM(Fg_in_T(i,:))
  END DO

END SUBROUTINE GQD_In_Calc

!==================================================================================================================================!
!
!==================================================================================================================================!
SUBROUTINE EGP_FV_NEWT()

  !OUTPUTS
  REAL*8,INTENT(INOUT):: E_avg(:,:), E_edgV(:,:), E_edgH(:,:), Temp(:,:), KapE_Bar(:,:)
  REAL*8,INTENT(OUT):: Fx_edgV(:,:), Fy_edgH(:,:)
  REAL*8,INTENT(OUT):: Q_bar(:,:), Q_bar_dT(:,:,:), KapE_Bar_dT(:,:,:)
  REAL*8,INTENT(OUT):: Del_E_avg(:,:,:), Del_E_edgV(:,:,:), Del_E_edgH(:,:,:), Del_T(:,:,:)
  REAL*8,INTENT(OUT):: Del_Fx_edgV(:,:,:), Del_Fy_edgH(:,:,:)
  INTEGER,INTENT(OUT):: Its

  !INPUTS
  REAL*8,INTENT(IN):: Temp_Mold(:,:), KapE_bar_Mold(:,:)
  REAL*8,INTENT(IN):: Theta, c, Delt
  REAL*8,INTENT(IN):: Delx(:), Dely(:), A(:,:)
  REAL*8,INTENT(IN):: Cb_L(:), Cb_B(:), Cb_R(:), Cb_T(:)
  REAL*8,INTENT(IN):: E_in_L(:), E_in_B(:), E_in_R(:), E_in_T(:)
  REAL*8,INTENT(IN):: F_in_L(:), F_in_B(:), F_in_R(:), F_in_T(:)
  REAL*8,INTENT(IN):: DC_xx(:,:), DL_xx(:,:), DR_xx(:,:)
  REAL*8,INTENT(IN):: DC_yy(:,:), DB_yy(:,:), DT_yy(:,:)
  REAL*8,INTENT(IN):: DL_xy(:,:), DR_xy(:,:), DB_xy(:,:), DT_xy(:,:)
  REAL*8,INTENT(IN):: PL(:,:), PR(:,:), PB(:,:), PT(:,:)
  REAL*8,INTENT(IN):: Gold_Hat(:,:), Rhat_old(:,:)
  REAL*8,INTENT(IN):: Chi, line_src
  REAL*8,INTENT(IN):: Kap0, cv, Comp_Unit

  !INTERNALS(ARRAYS)
  REAL*8,ALLOCATABLE:: EB_L(:,:), EB_B(:,:), EB_C(:,:), EB_R(:,:), EB_T(:,:)
  REAL*8,ALLOCATABLE:: MBx_C(:,:), MBx_R(:,:), MBx_B(:,:), MBx_T(:,:)
  REAL*8,ALLOCATABLE:: MBy_C(:,:), MBy_T(:,:), MBy_L(:,:), MBy_R(:,:)
  REAL*8,ALLOCATABLE:: Cp_L(:), Cp_B(:), Cp_R(:), Cp_T(:)
  REAL*8,ALLOCATABLE:: BC_L(:), BC_B(:), BC_R(:), BC_T(:)
  REAL*8,ALLOCATABLE:: MBx_RHS(:,:), MBy_RHS(:,:), W(:,:)
  REAL*8,ALLOCATABLE:: dr_B(:,:), dr_T(:,:), dr_ML(:,:), dr_MB(:,:), dr_MR(:,:), dr_MT(:,:), dr_G(:,:)
  REAL*8,ALLOCATABLE:: KapE_dT_Org(:,:)

  !INTERNALS(NON-ARRAYS)
  REAL*8,PARAMETER:: h=6.62613d-19    !erg*sh
  REAL*8,PARAMETER:: c=2.99792458d2   !cm/sh
  REAL*8,PARAMETER:: pi=4d0*ATAN(1d0)
  REAL*8,PARAMETER:: erg=6.24150934d11  !eV/erg -- 1 erg is this many ev's
  LOGICAL:: Converged
  INTEGER:: N_x, N_y, i, j
  REAL*8:: KapE_Bound


  !===========================================================================!
  !                                                                           !
  !     SETUP OF PARAMETERS BEFORE NEWTON ITERATIONS                          !
  !                                                                           !
  !===========================================================================!
  !calculating constants
  sig_R=2d0*pi**5/(15d0*c**2*h**3*erg**4*comp_unit)
  aR=4d0*sig_R/c

  !determining domain size
  N_x = SIZE(Temp,1)
  N_y = SIZE(Temp,2)

  !allocating all internal arrays
  ALLOCATE(EB_L(N_x,N_y),EB_B(N_x,N_y),EB_C(N_x,N_y),EB_R(N_x,N_y),EB_T(N_x,N_y))
  ALLOCATE(MBx_C(N_x,N_y),MBx_R(N_x-1,N_y),MBx_B(N_x,N_y),MBx_T(N_x,N_y))
  ALLOCATE(MBy_C(N_x,N_y),MBy_T(N_x,N_y-1),MBy_L(N_x,N_y),MBy_R(N_x,N_y))
  ALLOCATE(Cp_L(N_y),Cp_B(N_x),Cp_R(N_y),Cp_T(N_x))
  ALLOCATE(BC_L(N_y),BC_B(N_x),BC_R(N_y),BC_T(N_x))
  ALLOCATE(MBx_RHS(N_x,N_y),MBy_RHS(N_x,N_y),W(N_x,N_y))
  ALLOCATE(dr_B(N_x,N_y),dr_T(N_x,N_y),dr_ML(N_x,N_y),dr_MB(N_x,N_y),dr_MR(N_x,N_y),dr_MT(N_x,N_y),dr_G(N_x,N_y))
  ALLOCATE(KapE_dT_Org(N_x,N_y))

  !calculating KapE_Bar derivative with respect to temperature
  DO j=1,N_y
    DO i=1,N_x
      KapE_Bar_dT(i,j) = FC_KapE_Bar_dT(Temp(i,j),Temp_Mold(i,j),KapE_Bar(i,j),KapE_Bar_Mold(i,j))
    END DO
  END DO
  KapE_dT_Org = KapE_Bar_dT !storing KapE_Bar_dT for restoration at each Newton iteration

  !===========================================================================!
  !                                                                           !
  !     START OF NEWTON ITERATIONS                                            !
  !                                                                           !
  !===========================================================================!
  Converged = .FALSE.
  Its = 0
  DO WHILE (.NOT.Converged)
    Its = Its + 1

    !calculating source and source T-derivative
    DO j=1,N_y
      DO i=1,N_x
        Q_bar(i,j) = 15d0*kap0*c*aR*Temp(i,j)/pi**4
        Q_bar_dT(i,j) = 15d0*kap0*c*aR/pi**4
      END DO
    END DO

    !damping KapE_Bar_dT when/where too large
    KapE_Bar_dT = KapE_dT_Org !restoring original KapE_Bar derivative
    DO j=1,N_y
      DO i=1,N_x
        !if KapE_Bar_dT exceeds the given bound in cell (i,j), reduce for this iteration
        KapE_Bound = (cv/Delt + Q_bar_dT(i,j))/(c*E_avg(i,j))
        IF (KapE_Bar_dT(i,j) .GT. chi*KapE_Bound) KapE_Bar_dT(i,j) = chi*KapE_Bound
      END DO
    END DO

    !===========================================================================!
    !                                                                           !
    !     CALCULATING ALL COEFFICIENTS OF THE LINEAR SYSTEM                     !
    !                                                                           !
    !===========================================================================!
    DO j=1,N_y
      DO i=1,N_x
        !'delta r' quantities (momentum bal eqs)
        dr_ML(i,j) = A(i,j)*Fx_edgV(i,j)/2d0 + c*Dely(j)*(DC_xx(i,j)*E_avg(i,j) - DL_xx(i,j)*E_edgV(i,j)) +&
         c*Delx(i)*(DT_xy(i,j)*E_edgH(i,j+1) - DB_xy(i,j)*E_edgH(i,j))/2d0 - PL(i,j)
        dr_MR(i,j) = A(i,j)*Fx_edgV(i,j+1)/2d0 + c*Dely(j)*(DR_xx(i,j)*E_edgV(i+1,j) - DC_xx(i,j)*E_avg(i,j)) +&
         c*Delx(i)*(DT_xy(i,j)*E_edgH(i,j+1) - DB_xy(i,j)*E_edgH(i,j))/2d0 - PR(i,j)
        dr_MB(i,j) = A(i,j)*Fy_edgH(i,j)/2d0 + c*Delx(i)*(DC_yy(i,j)*E_avg(i,j) - DB_yy(i,j)*E_edgH(i,j)) +&
         c*Delx(i)*(DR_xy(i,j)*E_edgV(i+1,j) - DL_xy(i,j)*E_edgV(i,j))/2d0 - PB(i,j)
        dr_MT(i,j) = A(i,j)*Fy_edgH(i+1,j)/2d0 + c*Delx(i)*(DT_yy(i,j)*E_edgH(i,j+1) - DC_yy(i,j)*E_avg(i,j)) +&
         c*Delx(i)*(DR_xy(i,j)*E_edgV(i+1,j) - DL_xy(i,j)*E_edgV(i,j))/2d0 - PT(i,j)

        !'delta r' quantities (rad energy bal/ material energy bal)
        dr_B(i,j) = A(i,j)*(1d0/(Theta*Delt) + c**KapE_Bar(i,j))*E_avg(i,j) + Dely(j)*(Fx_edgV(i+1,j) - Fx_edgV(i,j)) +&
         Delx(i)*(Fy_edgH(i,j+1) - Fy_edgH(i,j)) - A(i,j)*Q_bar(i,j)
        dr_T(i,j) = cv/(Theta*Delt)*Temp(i,j) + Q_bar(i,j) - c*KapE_Bar(i,j)*E_avg(i,j)

        !'delta r_G' is right hand side for linearized ebal
        dr_G(i,j) = Gold_Hat(i,j) - dr_B(i,j) -&
         (A(i,j)/W(i,j))*( c*KapE_Bar_dT(i,j,Its)*E_avg(i,j) - Q_bar_dT(i,j,Its) )*(Rhat_old(i,j)-dr_T(i,j)) +&
         2d0*Dely(j)*(dr_MR(i,j) - dr_ML(i,j))/A(i,j) + 2d0*Delx(i)*(dr_MT(i,j) - dr_MB(i,j))/A(i,j)

        !'omega' in delta T eq
        W(i,j) = Cv/(Theta*Delt) + Q_bar_dT(i,j,Its) - c*KapE_Bar_dT(i,j,Its)*E_avg(i,j)

        !coefficients for the left hand side of the energy balance eq
        EB_C(i,j) = A(i,j)*( 1d0/(Theta*Delt) +&
         c*KapE_Bar(i,j)*( 1d0 + (c*KapE_Bar_dT(i,j,Its)*E_avg(i,j) - Q_bar_dT(i,j,Its))/W(i,j) ) ) +&
         4d0*c*Dely(j)**2*DC_xx(i,j)/A(i,j) + 4d0*c*Delx(i)**2*DC_yy(i,j)/A(i,j)
        EB_L(i,j) = -2d0*c*Dely(j)**2*DL_xx(i,j)/A(i,j)
        EB_R(i,j) = -2d0*c*Dely(j)**2*DR_xx(i,j)/A(i,j)
        EB_B(i,j) = -2d0*c*Delx(i)**2*DB_yy(i,j)/A(i,j)
        EB_T(i,j) = -2d0*c*Delx(i)**2*DT_yy(i,j)/A(i,j)

        !coefficients for the left hand side of the x-momentum balance eq
        MBx_C(i,j) = 2d0*c*Dely(j)*DC_xx(i,j)/A(i,j)
        MBx_T(i,j) = -c*Delx(i)*DT_xy(i,j)/A(i,j)
        MBx_B(i,j) = c*Delx(i)*DB_xy(i,j)/A(i,j)

        !coefficients for the left hand side of the y-momentum balance eq
        MBy_C(i,j) = 2d0*c*Delx(i)*DC_yy(i,j)/A(i,j)
        MBy_R(i,j) = -c*Dely(j)*DR_xy(i,j)/A(i,j)
        MBy_L(i,j) = c*Dely(j)*DL_xy(i,j)/A(i,j)
      END DO
    END DO

    DO j=1,N_y-1
      DO i=1,N_x
        !center cell face coefficient for the left hand side of the y-momentum balance eq
        MBy_T(i,j) = -2d0*c*Delx(i)*( DT_yy(i,j)/A(i,j) + DB_yy(i,j+1)/A(i,j+1) )

        !right hand side of y-momentum balance eq
        MBy_RHS(i,j) = dr_MT(i,j)/A(i,j) - dr_MB(i,j+1)/A(i,j+1)
      END DO
    END DO

    DO j=1,N_y
      DO i=1,N_x-1
        !center cell face coefficient for the left hand side of the x-momentum balance eq
        MBx_R(i,j) = -2d0*c*Dely(j)*( DR_yy(i,j)/A(i,j) + DL_yy(i+1,j)/A(i+1,j) )

        !right hand side of x-momentum balance eq
        MBx_RHS(i,j) = dr_MR(i,j)/A(i,j) - dr_ML(i+1,j)/A(i+1,j)
      END DO
    END DO

    !left/right boundary condition coefficients and right hand sides
    DO j=1,N_y
      Cp_L(j) = c*Cb_L(j) - 2d0*c*Dely(j)*DL_xx(1,j)/A(1,j)
      Cp_R(j) = c*Cb_R(j) + 2d0*c*Dely(j)*DR_xx(N_x,j)/A(N_x,j)
      BC_L(j) = c*Cb_L(j)*(E_in_L(j) - E_edgV(1,j)) + Fx_edgV(1,j) - F_in_L(j) - 2d0*dr_ML(1,j)/A(1,j)
      BC_R(j) = c*Cb_R(j)*(E_in_R(j) - E_edgV(N_x+1,j)) + Fx_edgV(N_x+1,j) - F_in_R(j) - 2d0*dr_MR(N_x,j)/A(N_x,j)
    END DO

    !bottom/top boundary condition coefficients and right hand sides
    DO i=1,N_x
      Cp_B(j) = c*Cb_B(i) - 2d0*c*Delx(i)*DB_yy(i,1)/A(i,1)
      Cp_T(j) = c*Cb_T(i) + 2d0*c*Delx(i)*DT_yy(i,N_y)/A(i,N_y)
      BC_B(j) = c*Cb_B(i)*(E_in_B(i) - E_edgH(i,1)) + Fy_edgH(i,1) - F_in_B(i) - 2d0*dr_MB(i,1)/A(i,1)
      BC_T(j) = c*Cb_T(i)*(E_in_T(i) - E_edgH(i,N_y+1)) + Fy_edgH(i,N_y+1) - F_in_T(i) - 2d0*dr_MT(i,N_y)/A(i,N_y)
    END DO

    !===========================================================================!
    !                                                                           !
    !     BUILDING AND INVERTING LINEAR SYSTEM                                  !
    !     (solving for Delta E's)                                               !
    !                                                                           !
    !===========================================================================!
    CALL QD_FV(Del_E_avg(:,:,Its),Del_E_edgV(:,:,Its),Del_E_edgH(:,:,Its),EB_L,EB_B,EB_C,EB_R,EB_T,MBx_C,MBx_R,MBx_B,MBx_T,MBy_C,&
      MBy_T,MBy_L,MBy_R,EB_RHS,MBx_RHS,MBy_RHS,Cp_L,Cp_B,Cp_R,Cp_T,BC_L,BC_B,BC_R,BC_T)

    !===========================================================================!
    !                                                                           !
    !     CALCULATING Delta T's & Delta F's from Delta E's                      !
    !                                                                           !
    !===========================================================================!
    DO j=1,N_y
      DO i=1,N_x
        !Delta T
        Del_T(i,j) = ( c*KapE_Bar(i,j)*Del_E_avg(i,j,Its) + Rhat_old(i,j) - dr_T(i,j) )/W(i,j)

        !Delta Fx on left(i-1/2) edges
        Del_Fx_edgV(i,j) = -2d0*c*Dely(j)*(DC_xx(i,j)*Del_E_avg(i,j,Its) - DL_xx(i,j)*Del_E_edgV(i,j,Its))/A(i,j) -&
         c*Delx(i)*(DT_xy(i,j)*Del_E_edgH(i,j+1,Its) - DT_xy(i,j)*Del_E_edgH(i,j,Its))/A(i,j) - 2d0*dr_ML(i,j)/A(i,j)

        !Delta Fy on bottom(j-1/2) edges
        Del_Fy_edgH(i,j) = -2d0*c*Delx(i)*(DC_yy(i,j)*Del_E_avg(i,j,Its) - DB_yy(i,j)*Del_E_edgH(i,j,Its))/A(i,j) -&
         c*Dely(j)*(DR_xy(i,j)*Del_E_edgH(i+1,j,Its) - DB_xy(i,j)*Del_E_edgH(i,j,Its))/A(i,j) - 2d0*dr_MB(i,j)/A(i,j)
      END DO
      !Delta Fx on right-most edge
      i = N_x
      Del_Fx_edgV(i+1,j) = -2d0*c*Dely(j)*(DR_xx(i,j)*Del_E_edgV(i+1,j,Its) - DC_xx(i,j)*Del_E_avg(i,j,Its))/A(i,j) -&
       c*Delx(i)*(DT_xy(i,j)*Del_E_edgH(i,j+1,Its) - DT_xy(i,j)*Del_E_edgH(i,j,Its))/A(i,j) - 2d0*dr_ML(i,j)/A(i,j)
    END DO

    j = N_y
    DO i=1,N_x
      !Delta Fy onTOP-MOST edge
      Del_Fy_edgH(i,j+1) = -2d0*c*Delx(i)*(Dt_yy(i,j)*Del_E_edgH(i,j+1,Its) - Dc_yy(i,j)*Del_E_avg(i,j,Its))/A(i,j) -&
       c*Dely(j)*(DR_xy(i,j)*Del_E_edgH(i+1,j,Its) - DB_xy(i,j)*Del_E_edgH(i,j,Its))/A(i,j) - 2d0*dr_MB(i,j)/A(i,j)
    END DO

    !===========================================================================!
    !                                                                           !
    !     PERFORMING LINE SEARCH                                                !
    !                                                                           !
    !===========================================================================!

  END DO

END SUBROUTINE EGP_FV_NEWT

!==================================================================================================================================!
!
!==================================================================================================================================!
FUNCTION FC_KapE_Bar_dT(Temp_l,Temp_lold,KapE_Bar_l,KapE_Bar_lold)
  REAL*8:: FC_KapE_Bar_dT
  REAL*8,INTENT(IN):: Temp_l, Temp_lold, KapE_Bar_l, KapE_Bar_lold

  FC_KapE_Bar_dT = (KapE_Bar_l - KapE_Bar_lold)/(Temp_l - Temp_lold)

END FUNCTION

!==================================================================================================================================!
!
!==================================================================================================================================!

END MODULE GLOQD_SOLVES
