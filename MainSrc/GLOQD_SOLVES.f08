MODULE GLOQD_SOLVES

  USE QD_SOLVERS

  IMPLICIT NONE

CONTAINS

!============================================================================================================!
!
!============================================================================================================!
SUBROUTINE OLD_GREY_COEFS()

  REAL*8,INTENT(IN):: c, Delt, Theta
  REAL*8,INTENT(IN):: Gold(:,:,:), Temp_old(:,:), KapE_bar_old(:,:), E_avg_old(:,:), GQD_Src_old(:,:)
  REAL*8,INTENT(OUT):: Gold_Bar(:,:), Rhat_old(:,:)

  INTEGER:: N_x, N_y, N_g, i, j, g

  N_x = SIZE(Gold,1)
  N_y = SIZE(Gold,2)
  N_g = SIZE(Gold,3)

  Gold_Bar = 0d0
  DO g=1,N_g
    DO j=1,N_y
      DO i=1,N_x
        Gold_Bar(i,j) = Gold_Bar(i,j) + Gold(i,j,g)
        Rhat_old(i,j) = cv/(Theta*Delt)*Temp_old(i,j) + ((1d0-Theta)/Theta)*(c*KapE_bar_old(i,j)*E_avg_old(i,j) - GQD_Src_old(i,j))
      END DO
    END DO
  END DO


END SUBROUTINE OLD_GREY_COEFS

!============================================================================================================!
!
!============================================================================================================!
SUBROUTINE GREY_COEFS()

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
        DL_xx(i,j) = DL_xx(i,j) + fg_edgV_xx(i,j,g)*Eg_edgV(i,j,g)/Tilde_KapR
        DR_xx(i,j) = DR_xx(i,j) + fg_edgV_xx(i+1,j,g)*Eg_edgV(i+1,j,g)/Tilde_KapR
        DC_yy(i,j) = DC_yy(i,j) + fg_avg_yy(i,j,g)*Eg_avg(i,j,g)/Tilde_KapR
        DB_yy(i,j) = DB_yy(i,j) + fg_edgH_yy(i,j,g)*Eg_edgH(i,j,g)/Tilde_KapR
        DT_yy(i,j) = DT_yy(i,j) + fg_edgH_yy(i,j+1,g)*Eg_edgH(i,j+1,g)/Tilde_KapR
        DL_xy(i,j) = DL_xy(i,j) + fg_edgV_xy(i,j,g)*Eg_edgV(i,j,g)/Tilde_KapR
        DR_xy(i,j) = DR_xy(i,j) + fg_edgV_xy(i+1,j,g)*Eg_edgV(i+1,j,g)/Tilde_KapR
        DB_xy(i,j) = DB_xy(i,j) + fg_edgH_xy(i,j,g)*Eg_edgH(i,j,g)/Tilde_KapR
        DT_xy(i,j) = DT_xy(i,j) + fg_edgH_xy(i,j+1,g)*Eg_edgH(i,j+1,g)/Tilde_KapR
        PL(i,j) = PL(i,j) + (A(i,j)/(2d0*Theta*c*Delt)*Fxg_edgV(i,j,g) + Pold_L(i,j,g))/Tilde_KapR
        PR(i,j) = PR(i,j) + (A(i,j)/(2d0*Theta*c*Delt)*Fxg_edgV(i+1,j,g) + Pold_R(i,j,g))/Tilde_KapR
        PB(i,j) = PB(i,j) + (A(i,j)/(2d0*Theta*c*Delt)*Fyg_edgH(i,j,g) + Pold_B(i,j,g))/Tilde_KapR
        PT(i,j) = PT(i,j) + (A(i,j)/(2d0*Theta*c*Delt)*Fyg_edgH(i,j+1,g) + Pold_T(i,j,g))/Tilde_KapR
      END DO
    END DO
  END DO

  DO j=1,N_y
    DO i=1,N_x
      sum1 = SUM(Eg_avg(i,j,:))
      sum2 = SUM(Eg_edgV(i,j,:))
      sum3 = SUM(Eg_edgV(i+1,j,:))
      sum4 = SUM(Eg_edgH(i,j,:))
      sum5 = SUM(Eg_edgH(i,j+1,:))
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

!============================================================================================================!
!
!============================================================================================================!
SUBROUTINE EGP_FV_NEWT()

  REAL*8,INTENT(OUT):: E_avg(:,:), E_edgV(:,:), E_edgH(:,:), Temp(:,:)

  REAL*8,INTENT(IN):: Delx(:), Dely(:), A(:,:)

  REAL*8,ALLOCATABLE:: Del_E_avg(:,:), Del_E_edgV(:,:), Del_E_edgH(:,:), Del_Temp(:,:)
  REAL*8,ALLOCATABLE:: EB_L(:,:), EB_B(:,:), EB_C(:,:), EB_R(:,:), EB_T(:,:)
  REAL*8,ALLOCATABLE:: MBx_C(:,:), MBx_R(:,:), MBx_B(:,:), MBx_T(:,:)
  REAL*8,ALLOCATABLE:: MBy_C(:,:), MBy_T(:,:), MBy_L(:,:), MBy_R(:,:)
  REAL*8,ALLOCATABLE:: Cp_L(:), Cp_B(:), Cp_R(:), Cp_T(:)
  REAL*8,ALLOCATABLE:: BC_L(:), BC_B(:), BC_R(:), BC_T(:)
  REAL*8,ALLOCATABLE:: MBx_RHS(:,:), MBy_RHS(:,:)

END SUBROUTINE EGP_FV_NEWT

!============================================================================================================!
!
!============================================================================================================!

END MODULE GLOQD_SOLVES
