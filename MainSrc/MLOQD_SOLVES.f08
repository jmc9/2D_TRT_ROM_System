MODULE MLOQD_SOLVES

  USE QD_SOLVERS

  IMPLICIT NONE

CONTAINS

!==================================================================================================================================!
!
!==================================================================================================================================!
SUBROUTINE OLD_MGQD_COEFS(Eg_avg_old,Eg_edgV_old,Eg_edgH_old,Fxg_edgV_old,Fyg_edgH_old,fg_avg_xx_old,fg_avg_yy_old,&
  fg_edgV_xx_old,fg_edgV_xy_old,fg_edgH_yy_old,fg_edgH_xy_old,KapE_old,KapR_old,Src_old,Delx,Dely,A,c,Delt,Theta,&
  G_old,Pold_L,Pold_B,Pold_R,Pold_T)

  REAL*8,INTENT(IN):: Eg_avg_old(:,:,:), Eg_edgV_old(:,:,:), Eg_edgH_old(:,:,:)
  REAL*8,INTENT(IN):: Fxg_edgV_old(:,:,:), Fyg_edgH_old(:,:,:)
  REAL*8,INTENT(IN):: fg_avg_xx_old(:,:,:), fg_avg_yy_old(:,:,:)
  REAL*8,INTENT(IN):: fg_edgV_xx_old(:,:,:), fg_edgV_xy_old(:,:,:)
  REAL*8,INTENT(IN):: fg_edgH_yy_old(:,:,:), fg_edgH_xy_old(:,:,:)
  REAL*8,INTENT(IN):: KapE_old(:,:,:), KapR_old(:,:,:), Src_old(:,:,:)
  REAL*8,INTENT(IN):: Delx(:), Dely(:), A(:,:)
  REAL*8,INTENT(IN):: c, Delt, Theta

  REAL*8,INTENT(OUT):: G_old(:,:,:), Pold_L(:,:,:), Pold_B(:,:,:), Pold_R(:,:,:), Pold_T(:,:,:)

  INTEGER:: N_x, N_y, N_g, i, j, g

  N_x = SIZE(Eg_avg_old,1)
  N_y = SIZE(Eg_avg_old,2)
  N_g = SIZE(Eg_avg_old,3)

  DO g=1,N_g
    DO j=1,N_y
      DO i=1,N_x

        G_old(i,j,g) = ((Theta-1d0)/Theta)*( &
        Dely(j)*(Fxg_edgV_old(i+1,j,g)-Fxg_edgV_old(i,j,g)) &
        + Delx(i)*(Fyg_edgH_old(i,j+1,g)-Fyg_edgH_old(i,j,g)) &
        + A(i,j)*(c*KapE_old(i,j,g)*Eg_avg_old(i,j,g)-Src_old(i,j,g)) )

        Pold_L(i,j,g) = ((Theta-1d0)/Theta)*( &
        c*Dely(j)*(fg_avg_xx_old(i,j,g)*Eg_avg_old(i,j,g)-fg_edgV_xx_old(i,j,g)*Eg_edgV_old(i,j,g)) &
        + c*Delx(i)*(fg_edgH_xy_old(i,j+1,g)*Eg_edgH_old(i,j+1,g)-fg_edgH_xy_old(i,j,g)*Eg_edgH_old(i,j,g))/2d0 &
        + A(i,j)*(KapR_old(i,j,g)*Fxg_edgV_old(i,j,g))/2d0 )

        Pold_B(i,j,g) = ((Theta-1d0)/Theta)*( &
        c*Delx(i)*(fg_avg_yy_old(i,j,g)*Eg_avg_old(i,j,g)-fg_edgH_yy_old(i,j,g)*Eg_edgH_old(i,j,g)) &
        + c*Dely(j)*(fg_edgV_xy_old(i+1,j,g)*Eg_edgV_old(i+1,j,g)-fg_edgV_xy_old(i,j,g)*Eg_edgV_old(i,j,g))/2d0 &
        + A(i,j)*(KapR_old(i,j,g)*Fyg_edgH_old(i,j,g))/2d0 )

        Pold_R(i,j,g) = ((Theta-1d0)/Theta)*( &
        c*Dely(j)*(fg_edgV_xx_old(i+1,j,g)*Eg_edgV_old(i+1,j,g)-fg_avg_xx_old(i,j,g)*Eg_avg_old(i,j,g)) &
        + c*Delx(i)*(fg_edgH_xy_old(i,j+1,g)*Eg_edgH_old(i,j+1,g)-fg_edgH_xy_old(i,j,g)*Eg_edgH_old(i,j,g))/2d0 &
        + A(i,j)*(KapR_old(i,j,g)*Fxg_edgV_old(i+1,j,g))/2d0 )

        Pold_T(i,j,g) = ((Theta-1d0)/Theta)*( &
        c*Delx(i)*(fg_edgH_yy_old(i,j+1,g)*Eg_edgH_old(i,j+1,g)-fg_avg_yy_old(i,j,g)*Eg_avg_old(i,j,g)) &
        + c*Dely(j)*(fg_edgV_xy_old(i+1,j,g)*Eg_edgV_old(i+1,j,g)-fg_edgV_xy_old(i,j,g)*Eg_edgV_old(i,j,g))/2d0 &
        + A(i,j)*(KapR_old(i,j,g)*Fyg_edgH_old(i,j+1,g))/2d0 )

      END DO
    END DO
  END DO


END SUBROUTINE OLD_MGQD_COEFS

!==================================================================================================================================!
!
!==================================================================================================================================!
SUBROUTINE MLOQD_FV(Eg_avg,Eg_edgV,Eg_edgH,Fxg_edgV,Fyg_edgH,fg_avg_xx,fg_avg_xy,fg_avg_yy,fg_edgV_xx,fg_edgV_xy,fg_edgH_yy,&
  fg_edgH_xy,Cg_L,Cg_B,Cg_R,Cg_T,Eg_in_L,Eg_in_B,Eg_in_R,Eg_in_T,Fg_in_L,Fg_in_B,Fg_in_R,Fg_in_T,Src,KapE,KapR,Delx,Dely,A,c,&
  Delt,Theta,Open_Threads,Res_Calc,MGQD_Residual,G_old,Pold_L,Pold_B,Pold_R,Pold_T,Eg_avg_old,Eg_edgV_old,Eg_edgH_old,&
  Fxg_edgV_old,Fyg_edgH_old)

  REAL*8,INTENT(OUT):: Eg_avg(:,:,:), Eg_edgV(:,:,:), Eg_edgH(:,:,:)
  REAL*8,INTENT(OUT):: Fxg_edgV(:,:,:), Fyg_edgH(:,:,:)
  REAL*8,INTENT(OUT):: MGQD_Residual(:,:,:,:,:)

  REAL*8,INTENT(IN):: G_old(:,:,:), Pold_L(:,:,:), Pold_B(:,:,:), Pold_R(:,:,:), Pold_T(:,:,:)
  REAL*8,INTENT(IN):: fg_avg_xx(:,:,:), fg_avg_xy(:,:,:), fg_avg_yy(:,:,:)
  REAL*8,INTENT(IN):: fg_edgV_xx(:,:,:), fg_edgV_xy(:,:,:)
  REAL*8,INTENT(IN):: fg_edgH_yy(:,:,:), fg_edgH_xy(:,:,:)
  REAL*8,INTENT(IN):: Cg_L(:,:), Cg_B(:,:), Cg_R(:,:), Cg_T(:,:)
  REAL*8,INTENT(IN):: Eg_in_L(:,:), Eg_in_B(:,:), Eg_in_R(:,:), Eg_in_T(:,:)
  REAL*8,INTENT(IN):: Fg_in_L(:,:), Fg_in_B(:,:), Fg_in_R(:,:), Fg_in_T(:,:)
  REAL*8,INTENT(IN):: Src(:,:,:)
  REAL*8,INTENT(IN):: Eg_avg_old(:,:,:), Eg_edgV_old(:,:,:), Eg_edgH_old(:,:,:)
  REAL*8,INTENT(IN):: Fxg_edgV_old(:,:,:), Fyg_edgH_old(:,:,:)
  REAL*8,INTENT(IN):: KapE(:,:,:), KapR(:,:,:)
  REAL*8,INTENT(IN):: Delx(:), Dely(:), A(:,:)
  REAL*8,INTENT(IN):: c, Delt
  REAL*8,INTENT(IN):: Theta
  INTEGER,INTENT(IN):: Open_Threads
  LOGICAL,INTENT(IN):: Res_Calc

  REAL*8,ALLOCATABLE:: Xi(:,:)
  REAL*8,ALLOCATABLE:: MBx_RHS(:,:), MBy_RHS(:,:)
  REAL*8,ALLOCATABLE:: Ghat(:,:)
  REAL*8,ALLOCATABLE:: Phat_L(:,:), Phat_B(:,:), Phat_R(:,:), Phat_T(:,:)

  REAL*8,ALLOCATABLE:: EB_L(:,:), EB_B(:,:), EB_C(:,:), EB_R(:,:), EB_T(:,:)
  REAL*8,ALLOCATABLE:: MBx_C(:,:), MBx_R(:,:), MBx_B(:,:), MBx_T(:,:)
  REAL*8,ALLOCATABLE:: MBy_C(:,:), MBy_T(:,:), MBy_L(:,:), MBy_R(:,:)
  REAL*8,ALLOCATABLE:: Cp_L(:), Cp_B(:), Cp_R(:), Cp_T(:)
  REAL*8,ALLOCATABLE:: BC_L(:), BC_B(:), BC_R(:), BC_T(:)
  REAL*8,ALLOCATABLE:: E1(:,:),E2(:,:),E3(:,:)

  INTEGER:: N_g, N_x, N_y, Threads
  INTEGER:: i, j, g

  N_g = SIZE(Eg_avg,3)
  N_y = SIZE(Eg_avg,2)
  N_x = SIZE(Eg_avg,1)

  ALLOCATE(Xi(N_x,N_y))
  ALLOCATE(Ghat(N_x,N_y))
  ALLOCATE(Phat_L(N_x,N_y),Phat_B(N_x,N_y),Phat_R(N_x,N_y),Phat_T(N_x,N_y))

  ALLOCATE(EB_L(N_x,N_y),EB_B(N_x,N_y),EB_C(N_x,N_y),EB_R(N_x,N_y),EB_T(N_x,N_y))
  ALLOCATE(MBx_C(N_x,N_y),MBx_R(N_x-1,N_y),MBx_B(N_x,N_y),MBx_T(N_x,N_y))
  ALLOCATE(MBy_C(N_x,N_y),MBy_T(N_x,N_y-1),MBy_L(N_x,N_y),MBy_R(N_x,N_y))
  ALLOCATE(MBx_RHS(N_x,N_y),MBy_RHS(N_x,N_y))

  ALLOCATE(Cp_L(N_y), Cp_B(N_x), Cp_R(N_y), Cp_T(N_x))
  ALLOCATE(BC_L(N_y), BC_B(N_x), BC_R(N_y), BC_T(N_x))

  ALLOCATE(E1(N_x,N_y),E2(N_x+1,N_y),E3(N_x,N_y+1))

  DO g=1,N_g

    DO j=1,N_y
      DO i=1,N_x
        Xi(i,j) = A(i,j)*(1d0/(Theta*c*Delt)+KapR(i,j,g))/2d0

        Phat_L(i,j) = ( A(i,j)*Fxg_edgV_old(i,j,g)/(2d0*Theta*c*Delt) + Pold_L(i,j,g) )/Xi(i,j)
        Phat_B(i,j) = ( A(i,j)*Fyg_edgH_old(i,j,g)/(2d0*Theta*c*Delt) + Pold_B(i,j,g) )/Xi(i,j)
        Phat_R(i,j) = ( A(i,j)*Fxg_edgV_old(i+1,j,g)/(2d0*Theta*c*Delt) + Pold_R(i,j,g) )/Xi(i,j)
        Phat_T(i,j) = ( A(i,j)*Fyg_edgH_old(i,j+1,g)/(2d0*Theta*c*Delt) + Pold_T(i,j,g) )/Xi(i,j)

        Ghat(i,j) = A(i,j)*Src(i,j,g) + A(i,j)*Eg_avg_old(i,j,g)/(Theta*Delt) + G_old(i,j,g) &
        + Dely(j)*(Phat_L(i,j)-Phat_R(i,j)) + Delx(i)*(Phat_B(i,j)-Phat_T(i,j))
      END DO
    END DO

    DO j=1,N_y
      DO i=1,N_x
        EB_C(i,j) = ( A(i,j)*(1d0/(Theta*Delt)+c*KapE(i,j,g)) &
        + 2d0*c*Dely(j)**2*fg_avg_xx(i,j,g)/Xi(i,j) + 2d0*c*Delx(i)**2*fg_avg_yy(i,j,g)/Xi(i,j) )
        EB_R(i,j) = -c*Dely(j)**2*fg_edgV_xx(i+1,j,g)/Xi(i,j)
        EB_L(i,j) = -c*Dely(j)**2*fg_edgV_xx(i,j,g)/Xi(i,j)
        EB_T(i,j) = -c*Delx(i)**2*fg_edgH_yy(i,j+1,g)/Xi(i,j)
        EB_B(i,j) = -c*Delx(i)**2*fg_edgH_yy(i,j,g)/Xi(i,j)

        MBx_C(i,j) = c*Dely(j)*fg_avg_xx(i,j,g)/Xi(i,j)
        MBx_T(i,j) = -c*Delx(i)*fg_edgH_xy(i,j+1,g)/(2d0*Xi(i,j))
        MBx_B(i,j) = c*Delx(i)*fg_edgH_xy(i,j,g)/(2d0*Xi(i,j))

        MBy_C(i,j) = c*Delx(i)*fg_avg_yy(i,j,g)/Xi(i,j)
        MBy_R(i,j) = -c*Dely(j)*fg_edgV_xy(i+1,j,g)/(2d0*Xi(i,j))
        MBy_L(i,j) = c*Dely(j)*fg_edgV_xy(i,j,g)/(2d0*Xi(i,j))
      END DO
    END DO

    DO j=1,N_y
      DO i=1,N_x-1
        MBx_R(i,j) = -(c*Dely(j)/Xi(i,j) + c*Dely(j)/Xi(i+1,j))*fg_edgV_xx(i+1,j,g)
      END DO
    END DO

    DO j=1,N_y-1
      DO i=1,N_x
        MBy_T(i,j) = -(c*Delx(i)/Xi(i,j) + c*Delx(i)/Xi(i,j+1))*fg_edgH_yy(i,j+1,g)
      END DO
    END DO

    DO j=1,N_y
      Cp_L(j) = c*Cg_L(j,g) - c*Dely(j)*fg_edgV_xx(1,j,g)/Xi(1,j)
      Cp_R(j) = c*Cg_R(j,g) + c*Dely(j)*fg_edgV_xx(N_x+1,j,g)/Xi(N_x,j)
      BC_L(j) = c*Cg_L(j,g)*Eg_in_L(j,g) - Fg_in_L(j,g) + Phat_L(1,j)
      BC_R(j) = c*Cg_R(j,g)*Eg_in_R(j,g) - Fg_in_R(j,g) + Phat_R(N_x,j)
    END DO

    DO i=1,N_x
      Cp_B(i) = c*Cg_B(i,g) - c*Delx(i)*fg_edgH_yy(i,1,g)/Xi(i,1)
      Cp_T(i) = c*Cg_T(i,g) + c*Delx(i)*fg_edgH_yy(i,N_y+1,g)/Xi(i,N_y)
      BC_B(i) = c*Cg_B(i,g)*Eg_in_B(i,g) - Fg_in_B(i,g) + Phat_B(i,1)
      BC_T(i) = c*Cg_T(i,g)*Eg_in_T(i,g) - Fg_in_T(i,g) + Phat_T(i,N_y)
    END DO

    DO j=1,N_y
      DO i=1,N_x
        IF (i .NE. N_x) MBx_RHS(i,j) = Phat_L(i+1,j) - Phat_R(i,j)
        IF (j .NE. N_y) MBy_RHS(i,j) = Phat_B(i,j+1) - Phat_T(i,j)
      END DO
    END DO

    !===========================================================================!
    !                                                                           !
    !     Finding Eg's by inverting the reduced linear system                   !
    !                                                                           !
    !===========================================================================!
    ! CALL QD_FV(Eg_avg(:,:,g),Eg_edgV(:,:,g),Eg_edgH(:,:,g),EB_L,EB_B,EB_C,EB_R,EB_T,MBx_C,MBx_R,MBx_B,MBx_T,MBy_C,&
    !   MBy_T,MBy_L,MBy_R,Ghat,MBx_RHS,MBy_RHS,Cp_L,Cp_B,Cp_R,Cp_T,BC_L,BC_B,BC_R,BC_T)
    CALL QD_FV(E1,E2,E3,EB_L,EB_B,EB_C,EB_R,EB_T,MBx_C,MBx_R,MBx_B,MBx_T,MBy_C,&
      MBy_T,MBy_L,MBy_R,Ghat,MBx_RHS,MBy_RHS,Cp_L,Cp_B,Cp_R,Cp_T,BC_L,BC_B,BC_R,BC_T)
    Eg_avg(:,:,g)=E1
    Eg_edgV(:,:,g)=E2
    Eg_edgH(:,:,g)=E3

    !===========================================================================!
    !                                                                           !
    !     Finding Fg's with auxilliary equations                                !
    !                                                                           !
    !===========================================================================!
    DO j=1,N_y
      DO i=1,N_x
        Fxg_edgV(i,j,g) = -(c*Dely(j)/(xi(i,j)))*(fg_avg_xx(i,j,g)*Eg_avg(i,j,g) - fg_edgV_xx(i,j,g)*Eg_edgV(i,j,g)) - &
        (c*Delx(i)/(2d0*xi(i,j)))*(fg_edgH_xy(i,j+1,g)*Eg_edgH(i,j+1,g) - fg_edgH_xy(i,j,g)*Eg_edgH(i,j,g)) + Phat_L(i,j)

        Fyg_edgH(i,j,g) = -(c*Delx(i)/(xi(i,j)))*(fg_avg_yy(i,j,g)*Eg_avg(i,j,g) - fg_edgH_yy(i,j,g)*Eg_edgH(i,j,g)) - &
        (c*Dely(j)/(2d0*xi(i,j)))*(fg_edgV_xy(i+1,j,g)*Eg_edgV(i+1,j,g) - fg_edgV_xy(i,j,g)*Eg_edgV(i,j,g)) + Phat_B(i,j)
      END DO

      Fxg_edgV(N_x+1,j,g) = &
      -(c*Dely(j)/(xi(N_x,j)))*(fg_edgV_xx(N_x+1,j,g)*Eg_edgV(N_x+1,j,g) - fg_avg_xx(N_x,j,g)*Eg_avg(N_x,j,g)) - &
      (c*Delx(N_x)/(2d0*xi(N_x,j)))*(fg_edgH_xy(N_x,j+1,g)*Eg_edgH(N_x,j+1,g) - fg_edgH_xy(N_x,j,g)*Eg_edgH(N_x,j,g)) + &
      Phat_R(N_x,j)
    END DO
    DO i=1,N_x
      Fyg_edgH(i,N_y+1,g) = &
      -(c*Delx(i)/(xi(i,N_y)))*(fg_edgH_yy(i,N_y+1,g)*Eg_edgH(i,N_y+1,g) - fg_avg_yy(i,N_y,g)*Eg_avg(i,N_y,g)) - &
      (c*Dely(N_y)/(2d0*xi(i,N_y)))*(fg_edgV_xy(i+1,N_y,g)*Eg_edgV(i+1,N_y,g) - fg_edgV_xy(i,N_y,g)*Eg_edgV(i,N_y,g)) + &
      Phat_T(i,N_y)
    END DO

  END DO !End loop over N_g

  !--------------------------------------------------!
  !              Residual Calculations               !
  !--------------------------------------------------!
  IF (Res_Calc) THEN
  MGQD_Residual = 0d0
  DO g=1,N_g
    DO j=1,N_y
      DO i=1,N_x

        !Energy Balance
        MGQD_Residual(i,j,g,1,1) = A(i,j)*(1d0/(Theta*Delt)+c*KapE(i,j,g))*Eg_avg(i,j,g) + &
        Dely(j)*(Fxg_edgV(i+1,j,g)-Fxg_edgV(i,j,g)) + Delx(i)*(Fyg_edgH(i,j+1,g)-Fyg_edgH(i,j,g)) - &
        (A(i,j)*Src(i,j,g) + A(i,j)*Eg_avg_old(i,j,g)/(Theta*Delt) + G_old(i,j,g))

        MGQD_Residual(i,j,g,1,2) = MGQD_Residual(i,j,g,1,1)/Eg_avg(i,j,g)

        !Momentum Balance (left half-cell)
        MGQD_Residual(i,j,g,2,1) = A(i,j)*(1d0/(Theta*c*Delt)+KapR(i,j,g))*Fxg_edgV(i,j,g)/2d0 + &
        c*Dely(j)*(fg_avg_xx(i,j,g)*Eg_avg(i,j,g) - fg_edgV_xx(i,j,g)*Eg_edgV(i,j,g)) + &
        c*Delx(i)*(fg_edgH_xy(i,j+1,g)*Eg_edgH(i,j+1,g) - fg_edgH_xy(i,j,g)*Eg_edgH(i,j,g))/2d0 - &
        (A(i,j)*Fxg_edgV_old(i,j,g)/(2d0*Theta*c*Delt)+Pold_L(i,j,g))

        MGQD_Residual(i,j,g,2,2) = MGQD_Residual(i,j,g,2,1)/Fxg_edgV(i,j,g)

        !Momentum Balance (bottom half-cell)
        MGQD_Residual(i,j,g,3,1) = A(i,j)*(1d0/(Theta*c*Delt)+KapR(i,j,g))*Fyg_edgH(i,j,g)/2d0 + &
        c*Delx(i)*(fg_avg_yy(i,j,g)*Eg_avg(i,j,g) - fg_edgH_yy(i,j,g)*Eg_edgH(i,j,g)) + &
        c*Dely(j)*(fg_edgV_xy(i+1,j,g)*Eg_edgV(i+1,j,g) - fg_edgV_xy(i,j,g)*Eg_edgV(i,j,g))/2d0 - &
        (A(i,j)*Fyg_edgH_old(i,j,g)/(2d0*Theta*c*Delt)+Pold_B(i,j,g))

        MGQD_Residual(i,j,g,3,2) = MGQD_Residual(i,j,g,3,1)/Fyg_edgH(i,j,g)

        !Momentum Balance (right half-cell)
        MGQD_Residual(i,j,g,4,1) = A(i,j)*(1d0/(Theta*c*Delt)+KapR(i,j,g))*Fxg_edgV(i+1,j,g)/2d0 + &
        c*Dely(j)*(fg_edgV_xx(i+1,j,g)*Eg_edgV(i+1,j,g) - fg_avg_xx(i,j,g)*Eg_avg(i,j,g)) + &
        c*Delx(i)*(fg_edgH_xy(i,j+1,g)*Eg_edgH(i,j+1,g) - fg_edgH_xy(i,j,g)*Eg_edgH(i,j,g))/2d0 - &
        (A(i,j)*Fxg_edgV_old(i+1,j,g)/(2d0*Theta*c*Delt)+Pold_R(i,j,g))

        MGQD_Residual(i,j,g,4,2) = MGQD_Residual(i,j,g,4,1)/Fxg_edgV(i+1,j,g)

        !Momentum Balance (top half-cell)
        MGQD_Residual(i,j,g,5,1) = A(i,j)*(1d0/(Theta*c*Delt)+KapR(i,j,g))*Fyg_edgH(i,j+1,g)/2d0 + &
        c*Delx(i)*(fg_edgH_yy(i,j+1,g)*Eg_edgH(i,j+1,g) - fg_avg_yy(i,j,g)*Eg_avg(i,j,g)) + &
        c*Dely(j)*(fg_edgV_xy(i+1,j,g)*Eg_edgV(i+1,j,g) - fg_edgV_xy(i,j,g)*Eg_edgV(i,j,g))/2d0 - &
        (A(i,j)*Fyg_edgH_old(i,j+1,g)/(2d0*Theta*c*Delt)+Pold_T(i,j,g))

        MGQD_Residual(i,j,g,5,2) = MGQD_Residual(i,j,g,5,1)/Fyg_edgH(i,j+1,g)

      END DO
    END DO
  END DO
  END IF

  DEALLOCATE(Xi)
  DEALLOCATE(Ghat)
  DEALLOCATE(Phat_L,Phat_B,Phat_R,Phat_T)

  DEALLOCATE(EB_L,EB_B,EB_C,EB_R,EB_T)
  DEALLOCATE(MBx_C,MBx_R,MBx_B,MBx_T)
  DEALLOCATE(MBy_C,MBy_T,MBy_L,MBy_R)
  DEALLOCATE(MBx_RHS,MBy_RHS)

  DEALLOCATE(Cp_L,Cp_B,Cp_R,Cp_T)
  DEALLOCATE(BC_L,BC_B,BC_R,BC_T)

  DEALLOCATE(E1,E2,E3)


END SUBROUTINE MLOQD_FV

!==================================================================================================================================!
!
!==================================================================================================================================!
SUBROUTINE Cg_Calc(Cg_L,Cg_B,Cg_R,Cg_T,I_edgV,I_edgH,Omega_x,Omega_y,quad_weight,Comp_Unit)
  REAL*8,INTENT(OUT):: Cg_L(:,:), Cg_B(:,:), Cg_R(:,:), Cg_T(:,:)
  REAL*8,INTENT(IN):: I_edgV(:,:,:,:), I_edgH(:,:,:,:)
  REAL*8,INTENT(IN):: Omega_x(:), Omega_y(:), quad_weight(:)
  REAL*8,INTENT(IN):: Comp_Unit

  REAL*8:: sum, eps
  INTEGER:: N_x, N_y, N_m, N_g
  INTEGER:: i, m, g, m1, m2

  !Determining array sizes
  N_x = SIZE(Cg_B,1)
  N_y = SIZE(Cg_L,1)
  N_m = SIZE(quad_weight,1)
  N_g = SIZE(Cg_L,2)

  !setting base eps value and scaling this based on the computational units
  eps=1d-25
  eps=eps/comp_unit

  Cg_L = 0d0
  Cg_B = 0d0
  Cg_R = 0d0
  Cg_T = 0d0
  DO g=1,N_g

    m1 = N_m/4+1
    m2 = 3*N_m/4
    DO i=1,N_y
      sum=0d0
      DO m=m1,m2
        Cg_L(i,g) = Cg_L(i,g) + Omega_x(m)*quad_weight(m)*(I_edgV(1,i,m,g) + eps)
        sum = sum + quad_weight(m)*(I_edgV(1,i,m,g) + eps)
      END DO
      Cg_L(i,g) = Cg_L(i,g)/sum
    END DO

    m1 = N_m/2+1
    m2 = N_m
    DO i=1,N_x
      sum=0d0
      DO m=m1,m2
        Cg_B(i,g) = Cg_B(i,g) + Omega_y(m)*quad_weight(m)*(I_edgH(i,1,m,g) + eps)
        sum = sum + quad_weight(m)*(I_edgH(i,1,m,g) + eps)
      END DO
      Cg_B(i,g) = Cg_B(i,g)/sum
    END DO

    DO i=1,N_y
      sum=0d0
      m1 = 1
      m2 = N_m/4
      DO m=m1,m2
        Cg_R(i,g) = Cg_R(i,g) + Omega_x(m)*quad_weight(m)*(I_edgV(SIZE(I_edgV,1),i,m,g) + eps)
        sum = sum + quad_weight(m)*(I_edgV(SIZE(I_edgV,1),i,m,g) + eps)
      END DO
      m1 = 3*N_m/4+1
      m2 = N_m
      DO m=m1,m2
        Cg_R(i,g) = Cg_R(i,g) + Omega_x(m)*quad_weight(m)*(I_edgV(SIZE(I_edgV,1),i,m,g) + eps)
        sum = sum + quad_weight(m)*(I_edgV(SIZE(I_edgV,1),i,m,g) + eps)
      END DO
      Cg_R(i,g) = Cg_R(i,g)/sum
    END DO

    m1 = 1
    m2 = N_m/2
    DO i=1,N_x
      sum=0d0
      DO m=m1,m2
        Cg_T(i,g) = Cg_T(i,g) + Omega_y(m)*quad_weight(m)*(I_edgH(i,SIZE(I_edgH,2),m,g) + eps)
        sum = sum + quad_weight(m)*(I_edgH(i,SIZE(I_edgH,2),m,g) + eps)
      END DO
      Cg_T(i,g) = Cg_T(i,g)/sum
    END DO

  END DO

END SUBROUTINE Cg_Calc

!==================================================================================================================================!
!
!==================================================================================================================================!
SUBROUTINE MGQD_In_Calc(Eg_in_L,Eg_in_B,Eg_in_R,Eg_in_T,Fg_in_L,Fg_in_B,Fg_in_R,Fg_in_T,I_edgV,I_edgH,Omega_x,Omega_y,&
  quad_weight,c,Comp_Unit)

  REAL*8,INTENT(OUT):: Eg_in_L(:,:), Eg_in_B(:,:), Eg_in_R(:,:), Eg_in_T(:,:)
  REAL*8,INTENT(OUT):: Fg_in_L(:,:), Fg_in_B(:,:), Fg_in_R(:,:), Fg_in_T(:,:)

  REAL*8,INTENT(IN):: I_edgV(:,:,:,:), I_edgH(:,:,:,:)
  REAL*8,INTENT(IN):: Omega_x(:), Omega_y(:), quad_weight(:)
  REAL*8,INTENT(IN):: c, Comp_Unit

  REAL*8:: eps
  INTEGER:: N_x, N_y, N_m, N_g
  INTEGER:: i, m, g, m1, m2

  !Determining array sizes
  N_x = SIZE(Eg_in_B,1)
  N_y = SIZE(Eg_in_L,1)
  N_m = SIZE(quad_weight,1)
  N_g = SIZE(Eg_in_L,2)

  !setting base eps value and scaling this based on the computational units
  eps=1d-25
  eps=eps/comp_unit

  Eg_in_L = 0d0
  Fg_in_L = 0d0
  Eg_in_B = 0d0
  Fg_in_B = 0d0
  Eg_in_R = 0d0
  Fg_in_R = 0d0
  Eg_in_T = 0d0
  Fg_in_T = 0d0
  DO g=1,N_g

    DO i=1,N_y
      m1 = 1
      m2 = N_m/4
      DO m=m1,m2
        Fg_in_L(i,g) = Fg_in_L(i,g) + Omega_x(m)*quad_weight(m)*(I_edgV(1,i,m,g) + eps)
        Eg_in_L(i,g) = Eg_in_L(i,g) + quad_weight(m)*(I_edgV(1,i,m,g) + eps)
      END DO
      m1 = 3*N_m/4+1
      m2 = N_m
      DO m=m1,m2
        Fg_in_L(i,g) = Fg_in_L(i,g) + Omega_x(m)*quad_weight(m)*(I_edgV(1,i,m,g) + eps)
        Eg_in_L(i,g) = Eg_in_L(i,g) + quad_weight(m)*(I_edgV(1,i,m,g) + eps)
      END DO
    END DO

    m1 = 1
    m2 = N_m/2
    DO i=1,N_x
      DO m=m1,m2
        Fg_in_B(i,g) = Fg_in_B(i,g) + Omega_y(m)*quad_weight(m)*(I_edgH(i,1,m,g) + eps)
        Eg_in_B(i,g) = Eg_in_B(i,g) + quad_weight(m)*(I_edgH(i,1,m,g) + eps)
      END DO
    END DO

    m1 = N_m/4+1
    m2 = 3*N_m/4
    DO i=1,N_y
      DO m=m1,m2
        Fg_in_R(i,g) = Fg_in_R(i,g) + Omega_x(m)*quad_weight(m)*(I_edgV(SIZE(I_edgV,1),i,m,g) + eps)
        Eg_in_R(i,g) = Eg_in_R(i,g) + quad_weight(m)*(I_edgV(SIZE(I_edgV,1),i,m,g) + eps)
      END DO
    END DO

    m1 = N_m/2+1
    m2 = N_m
    DO i=1,N_x
      DO m=m1,m2
        Fg_in_T(i,g) = Fg_in_T(i,g) + Omega_y(m)*quad_weight(m)*(I_edgH(i,SIZE(I_edgH,2),m,g) + eps)
        Eg_in_T(i,g) = Eg_in_T(i,g) + quad_weight(m)*(I_edgH(i,SIZE(I_edgH,2),m,g) + eps)
      END DO
    END DO

  END DO
  Eg_in_L = Eg_in_L/c
  Eg_in_B = Eg_in_B/c
  Eg_in_R = Eg_in_R/c
  Eg_in_T = Eg_in_T/c

END SUBROUTINE MGQD_In_Calc

!==================================================================================================================================!
!
!==================================================================================================================================!
SUBROUTINE fg_Calc(fg_avg_xx,fg_avg_xy,fg_avg_yy,fg_edgV_xx,fg_edgV_xy,fg_edgH_yy,fg_edgH_xy,Hg_avg_xx,Hg_avg_xy,Hg_avg_yy,&
  Hg_edgV_xx,Hg_edgV_xy,Hg_edgH_yy,Hg_edgH_xy,HO_Eg_avg,HO_Eg_edgV,HO_Eg_edgH,c,Comp_Unit)

  REAL*8,INTENT(OUT):: fg_avg_xx(:,:,:), fg_avg_xy(:,:,:), fg_avg_yy(:,:,:)
  REAL*8,INTENT(OUT):: fg_edgV_xx(:,:,:), fg_edgV_xy(:,:,:)
  REAL*8,INTENT(OUT):: fg_edgH_yy(:,:,:), fg_edgH_xy(:,:,:)

  REAL*8,INTENT(IN):: Hg_avg_xx(:,:,:), Hg_avg_xy(:,:,:), Hg_avg_yy(:,:,:)
  REAL*8,INTENT(IN):: Hg_edgV_xx(:,:,:), Hg_edgV_xy(:,:,:)
  REAL*8,INTENT(IN):: Hg_edgH_yy(:,:,:), Hg_edgH_xy(:,:,:)
  REAL*8,INTENT(IN):: HO_Eg_avg(:,:,:), HO_Eg_edgV(:,:,:), HO_Eg_edgH(:,:,:)
  REAL*8,INTENT(IN):: c, Comp_Unit

  REAL*8:: eps
  INTEGER:: N_x, N_y, N_g
  INTEGER:: i, j, g

  !Determining array sizes
  N_x = SIZE(fg_avg_xx,1)
  N_y = SIZE(fg_avg_xx,2)
  N_g = SIZE(fg_avg_xx,3)

  !setting base eps value and scaling this based on the computational units
  eps=1d-25
  eps=eps/comp_unit

  DO g=1,N_g
    DO j=1,N_y
      DO i=1,N_x

        fg_avg_xx(i,j,g) = Hg_avg_xx(i,j,g)/(c*HO_Eg_avg(i,j,g))
        fg_avg_xy(i,j,g) = Hg_avg_xy(i,j,g)/(c*HO_Eg_avg(i,j,g))
        fg_avg_yy(i,j,g) = Hg_avg_yy(i,j,g)/(c*HO_Eg_avg(i,j,g))

      END DO
    END DO

    DO j=1,N_y
      DO i=1,N_x+1

        fg_edgV_xx(i,j,g) = Hg_edgV_xx(i,j,g)/(c*HO_Eg_edgV(i,j,g))
        fg_edgV_xy(i,j,g) = Hg_edgV_xy(i,j,g)/(c*HO_Eg_edgV(i,j,g))

      END DO
    END DO

    DO j=1,N_y+1
      DO i=1,N_x

        fg_edgH_yy(i,j,g) = Hg_edgH_yy(i,j,g)/(c*HO_Eg_edgH(i,j,g))
        fg_edgH_xy(i,j,g) = Hg_edgH_xy(i,j,g)/(c*HO_Eg_edgH(i,j,g))

      END DO
    END DO
  END DO

END SUBROUTINE fg_Calc

!==================================================================================================================================!
!Subroutine COLLAPSE_MG_EF
!
! 'collapses' multigroup radiation energy densities and fluxes into low-order grey quantities:
! --> The total radiation energy density (E)
! --> The total radiation flux (F)
!==================================================================================================================================!
SUBROUTINE COLLAPSE_MG_EF(Eg_avg,Eg_edgV,Eg_edgH,Fxg_edgV,Fyg_edgH,Open_Threads,E_avg,E_edgV,E_edgH,Fx_edgV,Fy_edgH)

  REAL*8,INTENT(IN):: Eg_avg(:,:,:), Eg_edgV(:,:,:), Eg_edgH(:,:,:)
  REAL*8,INTENT(IN):: Fxg_edgV(:,:,:), Fyg_edgH(:,:,:)
  INTEGER,INTENT(IN):: Open_Threads

  REAL*8,INTENT(OUT):: E_avg(:,:), E_edgV(:,:), E_edgH(:,:)
  REAL*8,INTENT(OUT):: Fx_edgV(:,:), Fy_edgH(:,:)

  INTEGER:: i, j
  INTEGER:: N_x, N_y
  INTEGER:: Threads

  !Determining array sizes
  N_x = SIZE(Eg_avg,1)
  N_y = SIZE(Eg_avg,2)

  !Total E is the sum of all Eg's
  !Total F is the sum of all Fg's

  !cell-averaged values
  DO j=1,N_y
    DO i=1,N_x
      E_avg(i,j) = SUM(Eg_avg(i,j,:))
    END DO
  END DO

  !verticle cell-edge values
  DO j=1,N_y
    DO i=1,N_x+1
      E_edgV(i,j) = SUM(Eg_edgV(i,j,:))
      Fx_edgV(i,j) = SUM(Fxg_edgV(i,j,:))
    END DO
  END DO

  !horizontal cell-edge values
  DO j=1,N_y+1
    DO i=1,N_x
      E_edgH(i,j) = SUM(Eg_edgH(i,j,:))
      Fy_edgH(i,j) = SUM(Fyg_edgH(i,j,:))
    END DO
  END DO

END SUBROUTINE COLLAPSE_MG_EF

!==================================================================================================================================!
!
!==================================================================================================================================!

END MODULE MLOQD_SOLVES
