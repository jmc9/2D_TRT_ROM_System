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

  REAL*8,INTENT(OUT):: E_avg(:,:), E_edgV(:,:), E_edgH(:,:), Temp(:,:)

  REAL*8,INTENT(IN):: Delx(:), Dely(:), A(:,:)
  REAL*8,INTENT(IN):: Cb_L(:), Cb_B(:), Cb_R(:), Cb_T(:)
  REAL*8,INTENT(IN):: E_in_L(:), E_in_B(:), E_in_R(:), E_in_T(:)
  REAL*8,INTENT(IN):: F_in_L(:), F_in_B(:), F_in_R(:), F_in_T(:)
  REAL*8,INTENT(IN):: KapE_Bar(:,:)
  REAL*8,INTENT(IN):: DC_xx(:,:), DL_xx(:,:), DR_xx(:,:)
  REAL*8,INTENT(IN):: DC_yy(:,:), DB_yy(:,:), DT_yy(:,:)
  REAL*8,INTENT(IN):: DL_xy(:,:), DR_xy(:,:), DB_xy(:,:), DT_xy(:,:)
  REAL*8,INTENT(IN):: PL(:,:), PR(:,:), PB(:,:), PT(:,:)
  REAL*8,INTENT(IN):: Gold_Hat(:,:), Rhat_old(:,:)

  REAL*8,ALLOCATABLE:: Del_E_avg(:,:), Del_E_edgV(:,:), Del_E_edgH(:,:), Del_Temp(:,:)
  REAL*8,ALLOCATABLE:: EB_L(:,:), EB_B(:,:), EB_C(:,:), EB_R(:,:), EB_T(:,:)
  REAL*8,ALLOCATABLE:: MBx_C(:,:), MBx_R(:,:), MBx_B(:,:), MBx_T(:,:)
  REAL*8,ALLOCATABLE:: MBy_C(:,:), MBy_T(:,:), MBy_L(:,:), MBy_R(:,:)
  REAL*8,ALLOCATABLE:: Cp_L(:), Cp_B(:), Cp_R(:), Cp_T(:)
  REAL*8,ALLOCATABLE:: BC_L(:), BC_B(:), BC_R(:), BC_T(:)
  REAL*8,ALLOCATABLE:: MBx_RHS(:,:), MBy_RHS(:,:)

END SUBROUTINE EGP_FV_NEWT

!==================================================================================================================================!
!
!==================================================================================================================================!

END MODULE GLOQD_SOLVES
