MODULE GLOQD_SOLVES

  USE QD_SOLVERS

  IMPLICIT NONE

CONTAINS

!==================================================================================================================================!
!
!==================================================================================================================================!
SUBROUTINE OLD_GREY_COEFS(c,Delt,Theta,cv,A,Gold,Temp_old,KapE_bar_old,E_avg_old,GQD_Src_old,Gold_Hat,Rhat_old)

  REAL*8,INTENT(IN):: c, Delt, Theta, cv, A(:,:)
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
SUBROUTINE GREY_COEFS(Eg_avg,Eg_edgV,Eg_edgH,Fxg_edgV_old,Fyg_edgH_old,fg_avg_xx,fg_avg_yy,fg_edgV_xx,fg_edgV_xy,fg_edgH_yy,&
  fg_edgH_xy,KapE,KapR,A,c,Delt,Theta,Pold_L,Pold_R,Pold_B,Pold_T,KapE_Bar,DC_xx,DL_xx,DR_xx,DC_yy,DB_yy,DT_yy,DL_xy,&
  DR_xy,DB_xy,DT_xy,PL,PR,PB,PT)

  REAL*8,INTENT(IN):: Eg_avg(:,:,:), Eg_edgV(:,:,:), Eg_edgH(:,:,:)
  REAL*8,INTENT(IN):: Fxg_edgV_old(:,:,:), Fyg_edgH_old(:,:,:)
  REAL*8,INTENT(IN):: fg_avg_xx(:,:,:), fg_avg_yy(:,:,:)
  REAL*8,INTENT(IN):: fg_edgV_xx(:,:,:), fg_edgV_xy(:,:,:)
  REAL*8,INTENT(IN):: fg_edgH_yy(:,:,:), fg_edgH_xy(:,:,:)
  REAL*8,INTENT(IN):: KapE(:,:,:), KapR(:,:,:)
  REAL*8,INTENT(IN):: A(:,:)
  REAL*8,INTENT(IN):: c, Delt, Theta
  REAL*8,INTENT(IN):: Pold_L(:,:,:), Pold_R(:,:,:), Pold_B(:,:,:), Pold_T(:,:,:)

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

        !D parameters
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

        !P parameters
        PL(i,j) = PL(i,j) + (A(i,j)/(2d0*Theta*c*Delt)*Fxg_edgV_old(i,j,g) + Pold_L(i,j,g))/Tilde_KapR
        PR(i,j) = PR(i,j) + (A(i,j)/(2d0*Theta*c*Delt)*Fxg_edgV_old(i+1,j,g) + Pold_R(i,j,g))/Tilde_KapR
        PB(i,j) = PB(i,j) + (A(i,j)/(2d0*Theta*c*Delt)*Fyg_edgH_old(i,j,g) + Pold_B(i,j,g))/Tilde_KapR
        PT(i,j) = PT(i,j) + (A(i,j)/(2d0*Theta*c*Delt)*Fyg_edgH_old(i,j+1,g) + Pold_T(i,j,g))/Tilde_KapR
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
      Cb_T(i) = Cb_T(i) + Cg_T(i,g)*q
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
!Subroutine EGP_FV_NEWT
!
! Solves the 2D effective grey problem (EGP) discretized with a second-order finite volumes (FV) scheme with Newton's method
!
! NOTES::
!   The effective grey problem couples together the grey quasidiffusion system and the material energy balance equation
!   A reduced form of the linear system is used to solve ONLY for change in ENERGY DENSITIES (Delta E) via matrix inversion
!   Auxilliary equations are used to solve for other Delta terms (for other unknowns)
!
!   The Newton iterations employ a line search and safety bounds
!
!==================================================================================================================================!
SUBROUTINE EGP_FV_NEWT(E_avg,E_edgV,E_edgH,Temp,KapE_Bar,Fx_edgV,Fy_edgH,Q_bar,KapE_Bar_dT,Its,Deltas,dresiduals,Temp_Mold,&
  KapE_bar_Mold,Theta,Delt,Delx,Dely,A,Cb_L,Cb_B,Cb_R,Cb_T,E_in_L,E_in_B,E_in_R,E_in_T,F_in_L,F_in_B,F_in_R,F_in_T,&
  DC_xx,DL_xx,DR_xx,DC_yy,DB_yy,DT_yy,DL_xy,DR_xy,DB_xy,DT_xy,PL,PR,PB,PT,Gold_Hat,Rhat_old,Kap0,cv,Comp_Unit,Chi,&
  line_src,E_Bound_Low,T_Bound_Low,Eps1,Eps2,Maxits,MGQD_It,Use_Line_Search,Use_Safety_Search,Res_Calc,kapE_dT_flag,GQD_Kits)

  !OUTPUTS
  REAL*8,INTENT(INOUT):: E_avg(:,:), E_edgV(:,:), E_edgH(:,:)
  REAL*8,INTENT(INOUT):: Fx_edgV(:,:), Fy_edgH(:,:)
  REAL*8,INTENT(INOUT):: Temp(:,:), KapE_Bar(:,:)
  REAL*8,INTENT(OUT):: Q_bar(:,:), KapE_Bar_dT(:,:)
  INTEGER,INTENT(OUT):: Its, GQD_Kits(*)
  REAL*8,ALLOCATABLE,INTENT(OUT):: Deltas(:,:), dresiduals(:,:)

  !INPUTS
  REAL*8,INTENT(IN):: Temp_Mold(:,:), KapE_bar_Mold(:,:)
  REAL*8,INTENT(IN):: Theta, Delt
  REAL*8,INTENT(IN):: Delx(:), Dely(:), A(:,:)
  REAL*8,INTENT(IN):: Cb_L(:), Cb_B(:), Cb_R(:), Cb_T(:)
  REAL*8,INTENT(IN):: E_in_L(:), E_in_B(:), E_in_R(:), E_in_T(:)
  REAL*8,INTENT(IN):: F_in_L(:), F_in_B(:), F_in_R(:), F_in_T(:)
  REAL*8,INTENT(IN):: DC_xx(:,:), DL_xx(:,:), DR_xx(:,:)
  REAL*8,INTENT(IN):: DC_yy(:,:), DB_yy(:,:), DT_yy(:,:)
  REAL*8,INTENT(IN):: DL_xy(:,:), DR_xy(:,:), DB_xy(:,:), DT_xy(:,:)
  REAL*8,INTENT(IN):: PL(:,:), PR(:,:), PB(:,:), PT(:,:)
  REAL*8,INTENT(IN):: Gold_Hat(:,:), Rhat_old(:,:)
  REAL*8,INTENT(IN):: Kap0, cv, Comp_Unit
  REAL*8,INTENT(IN):: Chi, line_src, E_Bound_Low, T_Bound_Low, Eps1, Eps2
  INTEGER,INTENT(IN):: Maxits, MGQD_It
  LOGICAL,INTENT(IN):: Use_Line_Search, Use_Safety_Search, Res_Calc, kapE_dT_flag

  !INTERNALS(ARRAYS)
  REAL*8,ALLOCATABLE:: Del_T(:,:), Del_E_avg(:,:), Del_E_edgV(:,:), Del_E_edgH(:,:)
  REAL*8,ALLOCATABLE:: Del_Fx_edgV(:,:), Del_Fy_edgH(:,:)
  REAL*8,ALLOCATABLE:: Q_bar_dT(:,:)!, KapE_Bar_dT(:,:)
  REAL*8,ALLOCATABLE:: EB_L(:,:), EB_B(:,:), EB_C(:,:), EB_R(:,:), EB_T(:,:)
  REAL*8,ALLOCATABLE:: MBx_C(:,:), MBx_R(:,:), MBx_B(:,:), MBx_T(:,:)
  REAL*8,ALLOCATABLE:: MBy_C(:,:), MBy_T(:,:), MBy_L(:,:), MBy_R(:,:)
  REAL*8,ALLOCATABLE:: Cp_L(:), Cp_B(:), Cp_R(:), Cp_T(:)
  REAL*8,ALLOCATABLE:: BC_L(:), BC_B(:), BC_R(:), BC_T(:)
  REAL*8,ALLOCATABLE:: MBx_RHS(:,:), MBy_RHS(:,:), W(:,:)
  REAL*8,ALLOCATABLE:: dr_B(:,:), dr_T(:,:), dr_ML(:,:), dr_MB(:,:), dr_MR(:,:), dr_MT(:,:), dr_G(:,:)
  REAL*8,ALLOCATABLE:: KapE_dT_Org(:,:)
  REAL*8,ALLOCATABLE:: Deltas2(:,:)
  REAL*8,ALLOCATABLE:: Del_T_old(:,:), Del_E_avg_old(:,:), Del_E_edgV_old(:,:), Del_E_edgH_old(:,:)
  REAL*8,ALLOCATABLE:: Del_Fx_edgV_old(:,:), Del_Fy_edgH_old(:,:)

  !INTERNALS(NON-ARRAYS)
  REAL*8,PARAMETER:: h=6.62613d-19    !erg*sh
  REAL*8,PARAMETER:: c=2.99792458d2   !cm/sh
  REAL*8,PARAMETER:: pi=4d0*ATAN(1d0)
  REAL*8,PARAMETER:: erg=6.24150934d11  !eV/erg -- 1 erg is this many ev's
  LOGICAL:: Converged
  INTEGER:: N_x, N_y, i, j, LS_Maxit=10, LS_Its
  REAL*8:: KapE_Bound, dres, LS, aR, sig_R, LSbnd, norm, norm_old
  REAL*8:: T_Norm, T_Eps, Fx_Norm, Fx_Eps, Fy_Norm, Fy_Eps
  REAL*8:: Ea_Norm, Ea_Eps, Ev_Norm, Ev_Eps, Eh_Norm, Eh_Eps


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
  ALLOCATE(Del_T(N_x,N_y),Del_E_avg(N_x,N_y),Del_E_edgV(N_x+1,N_y),Del_E_edgH(N_x,N_y+1))
  ALLOCATE(Del_Fx_edgV(N_x+1,N_y),Del_Fy_edgH(N_x,N_y+1))
  ALLOCATE(Q_bar_dT(N_x,N_y))
  ALLOCATE(EB_L(N_x,N_y),EB_B(N_x,N_y),EB_C(N_x,N_y),EB_R(N_x,N_y),EB_T(N_x,N_y))
  ALLOCATE(MBx_C(N_x,N_y),MBx_R(N_x-1,N_y),MBx_B(N_x,N_y),MBx_T(N_x,N_y))
  ALLOCATE(MBy_C(N_x,N_y),MBy_T(N_x,N_y-1),MBy_L(N_x,N_y),MBy_R(N_x,N_y))
  ALLOCATE(Cp_L(N_y),Cp_B(N_x),Cp_R(N_y),Cp_T(N_x))
  ALLOCATE(BC_L(N_y),BC_B(N_x),BC_R(N_y),BC_T(N_x))
  ALLOCATE(MBx_RHS(N_x,N_y),MBy_RHS(N_x,N_y),W(N_x,N_y))
  ALLOCATE(dr_B(N_x,N_y),dr_T(N_x,N_y),dr_ML(N_x,N_y),dr_MB(N_x,N_y),dr_MR(N_x,N_y),dr_MT(N_x,N_y),dr_G(N_x,N_y))
  ALLOCATE(KapE_dT_Org(N_x,N_y))
  IF (Use_Line_Search) THEN
    ALLOCATE(Del_T_old(N_x,N_y),Del_E_avg_old(N_x,N_y),Del_E_edgV_old(N_x+1,N_y),Del_E_edgH_old(N_x,N_y+1))
    ALLOCATE(Del_Fx_edgV_old(N_x+1,N_y),Del_Fy_edgH_old(N_x,N_y+1))
  END IF
  Del_T = 0d0
  Del_E_avg = 0d0
  Del_E_edgV = 0d0
  Del_E_edgH = 0d0
  Del_Fx_edgV = 0d0
  Del_Fy_edgH = 0d0

  !calculating KapE_Bar derivative with respect to temperature
  IF ((MGQD_It .GT. 1).AND.(kapE_dT_flag)) THEN
    DO j=1,N_y
      DO i=1,N_x
        IF (ABS(Temp(i,j)-Temp_Mold(i,j)) .LT. 1d-15) THEN
          KapE_Bar_dT(i,j) = 0d0
        ELSE
          KapE_Bar_dT(i,j) = FC_KapE_Bar_dT(Temp(i,j),Temp_Mold(i,j),KapE_Bar(i,j),KapE_Bar_Mold(i,j))
        END IF
      END DO
    END DO
    KapE_dT_Org = KapE_Bar_dT !storing KapE_Bar_dT for restoration at each Newton iteration

  ELSE
    KapE_Bar_dT = 0d0

  END IF

  !===========================================================================!
  !                                                                           !
  !     START OF NEWTON ITERATIONS                                            !
  !                                                                           !
  !===========================================================================!
  Converged = .FALSE.
  Its = 0
  Newton_Its: DO WHILE ((.NOT.Converged).AND.(Its .LT. Maxits))
    Its = Its + 1

    Del_T_old = Del_T
    Del_E_avg_old = Del_E_avg
    Del_E_edgV_old = Del_E_edgV
    Del_E_edgH_old = Del_E_edgH
    Del_Fx_edgV_old = Del_Fx_edgV
    Del_Fy_edgH_old = Del_Fy_edgH

    !calculating source and source T-derivative
    DO j=1,N_y
      DO i=1,N_x
        Q_bar(i,j) = 15d0*kap0*c*aR*Temp(i,j)/pi**4
        Q_bar_dT(i,j) = 15d0*kap0*c*aR/pi**4
      END DO
    END DO

    IF ((MGQD_It .GT. 1).AND.(kapE_dT_flag)) THEN
      !damping KapE_Bar_dT when/where too large
      KapE_Bar_dT = KapE_dT_Org !restoring original KapE_Bar derivative
      DO j=1,N_y
        DO i=1,N_x
          !if KapE_Bar_dT exceeds the given bound in cell (i,j), reduce for this iteration
          KapE_Bound = (cv/(Theta*Delt) + Q_bar_dT(i,j))/(c*E_avg(i,j))
          IF (KapE_Bar_dT(i,j) .GT. chi*KapE_Bound) KapE_Bar_dT(i,j) = chi*KapE_Bound
        END DO
      END DO
    END IF

    !===========================================================================!
    !                                                                           !
    !     CALCULATING ALL COEFFICIENTS OF THE LINEAR SYSTEM                     !
    !                                                                           !
    !===========================================================================!
    DO j=1,N_y
      DO i=1,N_x
        !residual for material energy balance eq
        dr_T(i,j) = dr_T_Calc(Theta,c,cv,Delt,Temp(i,j),E_avg(i,j),KapE_Bar(i,j),Q_Bar(i,j),Rhat_old(i,j))

        !residual for rad energy balance eq
        dr_B(i,j) = dr_B_Calc(Theta,c,Delt,Delx(i),Dely(j),A(i,j),KapE_Bar(i,j),Q_Bar(i,j),E_avg(i,j),&
         Fx_edgV(i+1,j),Fx_edgV(i,j),Fy_edgH(i,j+1),Fy_edgH(i,j),Gold_hat(i,j))

        !residual for left-half-cell rad x-momentum balance eq
        dr_ML(i,j) = dr_ML_Calc(c,A(i,j),Delx(i),Dely(j),DC_xx(i,j),DL_xx(i,j),DT_xy(i,j),DB_xy(i,j),PL(i,j),&
         Fx_edgV(i,j),E_avg(i,j),E_edgV(i,j),E_edgH(i,j+1),E_edgH(i,j))

        !residual for right-half-cell rad x-momentum balance eq
        dr_MR(i,j) = dr_MR_Calc(c,A(i,j),Delx(i),Dely(j),DR_xx(i,j),DC_xx(i,j),DT_xy(i,j),DB_xy(i,j),PR(i,j),&
         Fx_edgV(i+1,j),E_edgV(i+1,j),E_avg(i,j),E_edgH(i,j+1),E_edgH(i,j))

        !residual for bottom-half-cell rad y-momentum balance eq
        dr_MB(i,j) = dr_MB_Calc(c,A(i,j),Delx(i),Dely(j),DC_yy(i,j),DB_yy(i,j),DL_xy(i,j),DR_xy(i,j),PB(i,j),&
         Fy_edgH(i,j),E_avg(i,j),E_edgH(i,j),E_edgV(i+1,j),E_edgV(i,j))

        !residual for top-half-cell rad y-momentum balance eq
        dr_MT(i,j) = dr_MT_Calc(c,A(i,j),Delx(i),Dely(j),DT_yy(i,j),DC_yy(i,j),DL_xy(i,j),DR_xy(i,j),PT(i,j),&
         Fy_edgH(i,j+1),E_edgH(i,j+1),E_avg(i,j),E_edgV(i+1,j),E_edgV(i,j))

        !'omega' in delta T eq
        W(i,j) = Cv/(Theta*Delt) + Q_bar_dT(i,j) - c*KapE_Bar_dT(i,j)*E_avg(i,j)

        !'delta r_G' is right hand side for linearized ebal
        dr_G(i,j) = - dr_B(i,j) +&
         dr_T(i,j)*(A(i,j)/W(i,j))*( c*KapE_Bar_dT(i,j)*E_avg(i,j) - Q_bar_dT(i,j) ) +&
         2d0*Dely(j)*(dr_MR(i,j) - dr_ML(i,j))/A(i,j) + 2d0*Delx(i)*(dr_MT(i,j) - dr_MB(i,j))/A(i,j)

        !coefficients for the left hand side of the energy balance eq
        EB_C(i,j) = A(i,j)*( 1d0/(Theta*Delt) +&
         c*KapE_Bar(i,j)*( 1d0 + (c*KapE_Bar_dT(i,j)*E_avg(i,j) - Q_bar_dT(i,j))/W(i,j) ) ) +&
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
        MBy_RHS(i,j) = 2d0*dr_MT(i,j)/A(i,j) - 2d0*dr_MB(i,j+1)/A(i,j+1)
      END DO
    END DO

    DO j=1,N_y
      DO i=1,N_x-1
        !center cell face coefficient for the left hand side of the x-momentum balance eq
        MBx_R(i,j) = -2d0*c*Dely(j)*( DR_xx(i,j)/A(i,j) + DL_xx(i+1,j)/A(i+1,j) )

        !right hand side of x-momentum balance eq
        MBx_RHS(i,j) = 2d0*dr_MR(i,j)/A(i,j) - 2d0*dr_ML(i+1,j)/A(i+1,j)
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
      Cp_B(i) = c*Cb_B(i) - 2d0*c*Delx(i)*DB_yy(i,1)/A(i,1)
      Cp_T(i) = c*Cb_T(i) + 2d0*c*Delx(i)*DT_yy(i,N_y)/A(i,N_y)
      BC_B(i) = c*Cb_B(i)*(E_in_B(i) - E_edgH(i,1)) + Fy_edgH(i,1) - F_in_B(i) - 2d0*dr_MB(i,1)/A(i,1)
      BC_T(i) = c*Cb_T(i)*(E_in_T(i) - E_edgH(i,N_y+1)) + Fy_edgH(i,N_y+1) - F_in_T(i) - 2d0*dr_MT(i,N_y)/A(i,N_y)
    END DO

    !===========================================================================!
    !                                                                           !
    !     BUILDING AND INVERTING LINEAR SYSTEM                                  !
    !     (solving for Delta E's)                                               !
    !                                                                           !
    !===========================================================================!
    CALL QD_FV(Del_E_avg,Del_E_edgV,Del_E_edgH,GQD_Kits(Its),EB_L,EB_B,EB_C,EB_R,EB_T,MBx_C,MBx_R,MBx_B,MBx_T,MBy_C,&
      MBy_T,MBy_L,MBy_R,dr_G,MBx_RHS,MBy_RHS,Cp_L,Cp_B,Cp_R,Cp_T,BC_L,BC_B,BC_R,BC_T)

    !===========================================================================!
    !                                                                           !
    !     CALCULATING Delta T's & Delta F's from Delta E's                      !
    !                                                                           !
    !===========================================================================!
    DO j=1,N_y
      DO i=1,N_x
        !Delta T
        Del_T(i,j) = ( c*KapE_Bar(i,j)*Del_E_avg(i,j) - dr_T(i,j) )/W(i,j)

        !Delta Fx on left(i-1/2) edges
        Del_Fx_edgV(i,j) = -2d0*c*Dely(j)*(DC_xx(i,j)*Del_E_avg(i,j) - DL_xx(i,j)*Del_E_edgV(i,j))/A(i,j) -&
         c*Delx(i)*(DT_xy(i,j)*Del_E_edgH(i,j+1) - DT_xy(i,j)*Del_E_edgH(i,j))/A(i,j) - 2d0*dr_ML(i,j)/A(i,j)

        !Delta Fy on bottom(j-1/2) edges
        Del_Fy_edgH(i,j) = -2d0*c*Delx(i)*(DC_yy(i,j)*Del_E_avg(i,j) - DB_yy(i,j)*Del_E_edgH(i,j))/A(i,j) -&
         c*Dely(j)*(DR_xy(i,j)*Del_E_edgV(i+1,j) - DL_xy(i,j)*Del_E_edgV(i,j))/A(i,j) - 2d0*dr_MB(i,j)/A(i,j)
      END DO
      !Delta Fx on right-most edge
      i = N_x
      Del_Fx_edgV(i+1,j) = -2d0*c*Dely(j)*(DR_xx(i,j)*Del_E_edgV(i+1,j) - DC_xx(i,j)*Del_E_avg(i,j))/A(i,j) -&
       c*Delx(i)*(DT_xy(i,j)*Del_E_edgH(i,j+1) - DT_xy(i,j)*Del_E_edgH(i,j))/A(i,j) - 2d0*dr_MR(i,j)/A(i,j)
    END DO

    j = N_y
    DO i=1,N_x
      !Delta Fy on top-most edge
      Del_Fy_edgH(i,j+1) = -2d0*c*Delx(i)*(DT_yy(i,j)*Del_E_edgH(i,j+1) - DC_yy(i,j)*Del_E_avg(i,j))/A(i,j) -&
       c*Dely(j)*(DR_xy(i,j)*Del_E_edgV(i+1,j) - DL_xy(i,j)*Del_E_edgV(i,j))/A(i,j) - 2d0*dr_MT(i,j)/A(i,j)
    END DO

    !===========================================================================!
    !                                                                           !
    !     PERFORMING LINE SEARCH                                                !
    !                                                                           !
    !===========================================================================!
    LS_Switch: IF ((Use_Line_Search).AND.(Its .GT. 1)) THEN

      !Delta T search
      norm = NORM2(Del_T)
      norm_old = NORM2(Del_T_old)
      CALL LINE_SEARCH(norm,norm_old,line_src,LS_Maxit,Del_T,Del_E_avg,Del_E_edgV,Del_E_edgH,Del_Fx_edgV,Del_Fy_edgH,LS_Its)

      !Delta E_avg search
      norm = NORM2(Del_E_avg)
      norm_old = NORM2(Del_E_avg_old)
      CALL LINE_SEARCH(norm,norm_old,line_src,LS_Maxit,Del_T,Del_E_avg,Del_E_edgV,Del_E_edgH,Del_Fx_edgV,Del_Fy_edgH,LS_Its)

      !Delta E_edgV search
      norm = NORM2(Del_E_edgV)
      norm_old = NORM2(Del_E_edgV_old)
      CALL LINE_SEARCH(norm,norm_old,line_src,LS_Maxit,Del_T,Del_E_avg,Del_E_edgV,Del_E_edgH,Del_Fx_edgV,Del_Fy_edgH,LS_Its)

      !Delta E_edgH search
      norm = NORM2(Del_E_edgH)
      norm_old = NORM2(Del_E_edgH_old)
      CALL LINE_SEARCH(norm,norm_old,line_src,LS_Maxit,Del_T,Del_E_avg,Del_E_edgV,Del_E_edgH,Del_Fx_edgV,Del_Fy_edgH,LS_Its)

      !Delta Fx_edgV search
      norm = NORM2(Del_Fx_edgV)
      norm_old = NORM2(Del_Fx_edgV_old)
      CALL LINE_SEARCH(norm,norm_old,line_src,LS_Maxit,Del_T,Del_E_avg,Del_E_edgV,Del_E_edgH,Del_Fx_edgV,Del_Fy_edgH,LS_Its)

      !Delta Fy_edgH search
      norm = NORM2(Del_Fy_edgH)
      norm_old = NORM2(Del_Fy_edgH_old)
      CALL LINE_SEARCH(norm,norm_old,line_src,LS_Maxit,Del_T,Del_E_avg,Del_E_edgV,Del_E_edgH,Del_Fx_edgV,Del_Fy_edgH,LS_Its)

    END IF LS_Switch

    !===========================================================================!
    !                                                                           !
    !     PERFORMING SAFETY LINE SEARCH                                         !
    !                                                                           !
    !===========================================================================!
    Safe_Switch: IF (Use_Safety_Search) THEN

    LS = 1d0
    Safe_Search: DO

      DO j=1,N_y
        DO i=1,N_x

          !--------------------------------------------------!
          !                                                  !
          !     Temperature                                  !
          !                                                  !
          !--------------------------------------------------!
          Safe_T: DO
            !if temperature found below lower bound, perform line search
            IF (Temp(i,j)+Del_T(i,j)/LS .LT. T_Bound_Low) THEN
              LS = LS*line_src !increase total line search step
              CYCLE Safe_T !cycle loop to check if residual still bigger

            ELSE !if/when temperature above lower bound, move out of loop
              EXIT Safe_T

            END IF
          END DO Safe_T
          IF (LS .GT. 1d0) THEN
            CALL LS_Reset(LS,Del_T,Del_E_avg,Del_E_edgV,Del_E_edgH,Del_Fx_edgV,Del_Fy_edgH)
            LS = 1d0
          END IF

          !--------------------------------------------------!
          !                                                  !
          !     cell-averaged Energy Density                 !
          !                                                  !
          !--------------------------------------------------!
          Safe_EA: DO
            !if E found below lower bound, perform line search
            IF (E_avg(i,j)+Del_E_avg(i,j)/LS .LT. E_Bound_Low) THEN
              LS = LS*line_src !increase total line search step
              CYCLE Safe_EA !cycle loop to check if residual still bigger

            ELSE !if/when E above lower bound, move out of loop
              EXIT Safe_EA

            END IF
          END DO Safe_EA
          IF (LS .GT. 1d0) THEN
            CALL LS_Reset(LS,Del_T,Del_E_avg,Del_E_edgV,Del_E_edgH,Del_Fx_edgV,Del_Fy_edgH)
            LS = 1d0
          END IF

        END DO
      END DO

      DO j=1,N_y
        DO i=1,N_x+1

          !--------------------------------------------------!
          !                                                  !
          !     cell-edge (vertical) Energy Density          !
          !                                                  !
          !--------------------------------------------------!
          Safe_EV: DO
            !if E found below lower bound, perform line search
            IF (E_edgV(i,j)+Del_E_edgV(i,j)/LS .LT. E_Bound_Low) THEN
              LS = LS*line_src !increase total line search step
              CYCLE Safe_EV !cycle loop to check if residual still bigger

            ELSE !if/when E above lower bound, move out of loop
              EXIT Safe_EV

            END IF
          END DO Safe_EV
          IF (LS .GT. 1d0) THEN
            CALL LS_Reset(LS,Del_T,Del_E_avg,Del_E_edgV,Del_E_edgH,Del_Fx_edgV,Del_Fy_edgH)
            LS = 1d0
          END IF

        END DO
      END DO

      DO j=1,N_y+1
        DO i=1,N_x

          !--------------------------------------------------!
          !                                                  !
          !     cell-edge (horizontal) Energy Density        !
          !                                                  !
          !--------------------------------------------------!
          Safe_EH: DO
            !if E found below lower bound, perform line search
            IF (E_edgH(i,j)+Del_E_edgH(i,j)/LS .LT. E_Bound_Low) THEN
              LS = LS*line_src !increase total line search step
              CYCLE Safe_EH !cycle loop to check if residual still bigger

            ELSE !if/when E above lower bound, move out of loop
              EXIT Safe_EH

            END IF
          END DO Safe_EH
          IF (LS .GT. 1d0) THEN
            CALL LS_Reset(LS,Del_T,Del_E_avg,Del_E_edgV,Del_E_edgH,Del_Fx_edgV,Del_Fy_edgH)
            LS = 1d0
          END IF

        END DO
      END DO

    END DO Safe_Search

    END IF Safe_Switch

    !===========================================================================!
    !                                                                           !
    !     EVALUATING NEW SOLUTION VECTOR                                        !
    !                                                                           !
    !===========================================================================!
    Temp = Temp + Del_T
    E_avg = E_avg + Del_E_avg
    E_edgH = E_edgH + Del_E_edgH
    E_edgV = E_edgV + Del_E_edgV
    Fx_edgV = Fx_edgV + Del_Fx_edgV
    Fy_edgH = Fy_edgH + Del_Fy_edgH
    KapE_Bar = KapE_Bar + KapE_Bar_dT*Del_T
    Q_bar = 15d0*kap0*c*aR*Temp/pi**4

    !===========================================================================!
    !                                                                           !
    !     CHECKING CONVERGENCE                                                  !
    !                                                                           !
    !===========================================================================!
    T_Eps = Eps1*MAXVAL(ABS(Temp)) + Eps2
    Ea_Eps = Eps1*MAXVAL(ABS(E_avg)) + Eps2
    Ev_Eps = Eps1*MAXVAL(ABS(E_edgV)) + Eps2
    Eh_Eps = Eps1*MAXVAL(ABS(E_edgH)) + Eps2
    Fx_Eps = Eps1*MAXVAL(ABS(Fx_edgV)) + Eps2
    Fy_Eps = Eps1*MAXVAL(ABS(Fy_edgH)) + Eps2

    T_Norm = MAXVAL(ABS(Del_T))
    Ea_Norm = MAXVAL(ABS(Del_E_avg))
    Ev_Norm = MAXVAL(ABS(Del_E_edgV))
    Eh_Norm = MAXVAL(ABS(Del_E_edgH))
    Fx_Norm = MAXVAL(ABS(Del_Fx_edgV))
    Fy_Norm = MAXVAL(ABS(Del_Fy_edgH))

    Converged = .TRUE.
    IF (T_Norm .GT. T_Eps) Converged = .FALSE.
    IF (Ea_Norm .GT. Ea_Eps) Converged = .FALSE.
    IF (Ev_Norm .GT. Ev_Eps) Converged = .FALSE.
    IF (Eh_Norm .GT. Eh_Eps) Converged = .FALSE.
    IF (Fx_Norm .GT. Fx_Eps) Converged = .FALSE.
    IF (Fy_Norm .GT. Fy_Eps) Converged = .FALSE.

    IF (Res_Calc) THEN
      IF (Its .EQ. 1) THEN
        ALLOCATE(Deltas(1,6))
        Deltas(Its,1) = T_Norm
        Deltas(Its,2) = Ea_Norm
        Deltas(Its,3) = Ev_Norm
        Deltas(Its,4) = Eh_Norm
        Deltas(Its,5) = Fx_Norm
        Deltas(Its,6) = Fy_Norm

        ALLOCATE(dresiduals(1,6))
        dresiduals(Its,1) = MAXVAL(ABS(dr_T))
        dresiduals(Its,2) = MAXVAL(ABS(dr_B))
        dresiduals(Its,3) = MAXVAL(ABS(dr_ML))
        dresiduals(Its,4) = MAXVAL(ABS(dr_MB))
        dresiduals(Its,5) = MAXVAL(ABS(dr_MR))
        dresiduals(Its,6) = MAXVAL(ABS(dr_MT))

      ELSE
        ALLOCATE(Deltas2(Its-1,6))
        Deltas2 = Deltas
        DEALLOCATE(Deltas)
        ALLOCATE(Deltas(Its,6))
        Deltas(1:Its-1,:) = Deltas2
        ! DEALLOCATE(Deltas2)
        Deltas(Its,1) = T_Norm
        Deltas(Its,2) = Ea_Norm
        Deltas(Its,3) = Ev_Norm
        Deltas(Its,4) = Eh_Norm
        Deltas(Its,5) = Fx_Norm
        Deltas(Its,6) = Fy_Norm

        Deltas2 = dresiduals
        DEALLOCATE(dresiduals)
        ALLOCATE(dresiduals(Its,6))
        dresiduals(1:Its-1,:) = Deltas2
        DEALLOCATE(Deltas2)
        dresiduals(Its,1) = MAXVAL(ABS(dr_T))
        dresiduals(Its,2) = MAXVAL(ABS(dr_B))
        dresiduals(Its,3) = MAXVAL(ABS(dr_ML))
        dresiduals(Its,4) = MAXVAL(ABS(dr_MB))
        dresiduals(Its,5) = MAXVAL(ABS(dr_MR))
        dresiduals(Its,6) = MAXVAL(ABS(dr_MT))

      END IF
    END IF

  END DO Newton_Its

  DEALLOCATE(Del_T,Del_E_avg,Del_E_edgV,Del_E_edgH)
  DEALLOCATE(Del_Fx_edgV,Del_Fy_edgH)
  DEALLOCATE(Q_bar_dT)
  DEALLOCATE(EB_L,EB_B,EB_C,EB_R,EB_T)
  DEALLOCATE(MBx_C,MBx_R,MBx_B,MBx_T)
  DEALLOCATE(MBy_C,MBy_T,MBy_L,MBy_R)
  DEALLOCATE(Cp_L,Cp_B,Cp_R,Cp_T)
  DEALLOCATE(BC_L,BC_B,BC_R,BC_T)
  DEALLOCATE(MBx_RHS,MBy_RHS,W)
  DEALLOCATE(dr_B,dr_T,dr_ML,dr_MB,dr_MR,dr_MT,dr_G)
  DEALLOCATE(KapE_dT_Org)
  IF (Use_Line_Search) THEN
    DEALLOCATE(Del_T_old,Del_E_avg_old,Del_E_edgV_old,Del_E_edgH_old)
    DEALLOCATE(Del_Fx_edgV_old,Del_Fy_edgH_old)
  END IF

END SUBROUTINE EGP_FV_NEWT

!==================================================================================================================================!
!
!==================================================================================================================================!
SUBROUTINE LINE_SEARCH(norm,norm_old,line_src,maxit,Del_T,Del_E_avg,Del_E_edgV,Del_E_edgH,Del_Fx_edgV,Del_Fy_edgH,Its)
  REAL*8,INTENT(IN):: norm, norm_old, line_src
  INTEGER,INTENT(IN):: maxit
  REAL*8,INTENT(INOUT):: Del_T(:,:), Del_E_avg(:,:), Del_E_edgV(:,:), Del_E_edgH(:,:), Del_Fx_edgV(:,:), Del_Fy_edgH(:,:)
  INTEGER,INTENT(OUT):: Its
  REAL*8:: LS

  Its = 0
  LS = 1d0
  SRC: DO WHILE (Its .LT. maxit)
    Its = Its + 1
    IF (norm/LS .GT. norm_old) THEN
      LS = LS*line_src !increase total line search step
      CYCLE SRC !cycle loop to check if residual still bigger

    ELSE !if/when residual detected smaller than last iterate, move out of loop
      EXIT SRC

    END IF
  END DO SRC
  IF (LS .GT. 1d0) CALL LS_Reset(LS,Del_T,Del_E_avg,Del_E_edgV,Del_E_edgH,Del_Fx_edgV,Del_Fy_edgH)

END SUBROUTINE LINE_SEARCH

!==================================================================================================================================!
!
!==================================================================================================================================!
SUBROUTINE LS_Reset(LS,Del_T,Del_E_avg,Del_E_edgV,Del_E_edgH,Del_Fx_edgV,Del_Fy_edgH)

  REAL*8,INTENT(INOUT):: Del_T(:,:), Del_E_avg(:,:), Del_E_edgV(:,:), Del_E_edgH(:,:), Del_Fx_edgV(:,:), Del_Fy_edgH(:,:)
  REAL*8,INTENT(IN):: LS

  Del_T = Del_T/LS
  Del_E_avg = Del_E_avg/LS
  Del_E_edgV = Del_E_edgV/LS
  Del_E_edgH = Del_E_edgH/LS
  Del_Fx_edgV = Del_Fx_edgV/LS
  Del_Fy_edgH = Del_Fy_edgH/LS

END SUBROUTINE LS_Reset

!==================================================================================================================================!
!
!==================================================================================================================================!
FUNCTION FC_KapE_Bar_dT(Temp_l,Temp_lold,KapE_Bar_l,KapE_Bar_lold)
  REAL*8:: FC_KapE_Bar_dT
  REAL*8,INTENT(IN):: Temp_l, Temp_lold, KapE_Bar_l, KapE_Bar_lold

  FC_KapE_Bar_dT = (KapE_Bar_l - KapE_Bar_lold)/(Temp_l - Temp_lold)

END FUNCTION FC_KapE_Bar_dT

!==================================================================================================================================!
!
!==================================================================================================================================!
FUNCTION dr_T_Calc(Theta,c,cv,Delt,Temp,E_avg,KapE_Bar,Q_Bar,Rhat)
  REAL*8:: dr_T_Calc
  REAL*8,INTENT(IN):: Theta, c, cv, Delt
  REAL*8,INTENT(IN):: Temp, E_avg, KapE_Bar, Q_Bar, Rhat

  dr_T_Calc = cv/(Theta*Delt)*Temp + Q_bar - c*KapE_Bar*E_avg - Rhat

END FUNCTION dr_T_Calc

!==================================================================================================================================!
!
!==================================================================================================================================!
FUNCTION dr_B_Calc(Theta,c,Delt,Delx,Dely,A,KapE_Bar,Q_Bar,E_avg,Fx_edgV_R,Fx_edgV_L,Fy_edgH_T,Fy_edgH_B,Gold_hat)
  REAL*8:: dr_B_Calc
  REAL*8,INTENT(IN):: Theta, c
  REAL*8,INTENT(IN):: Delt, Delx, Dely, A
  REAL*8,INTENT(IN):: KapE_Bar, Q_Bar, Gold_hat
  REAL*8,INTENT(IN):: E_avg, Fx_edgV_R, Fx_edgV_L, Fy_edgH_T, Fy_edgH_B

  dr_B_Calc = A*(1d0/(Theta*Delt) + c*KapE_Bar)*E_avg + Dely*(Fx_edgV_R - Fx_edgV_L) + Delx*(Fy_edgH_T - Fy_edgH_B) -&
   A*Q_bar - Gold_hat

END FUNCTION dr_B_Calc

!==================================================================================================================================!
!
!==================================================================================================================================!
FUNCTION dr_ML_Calc(c,A,Delx,Dely,DC_xx,DL_xx,DT_xy,DB_xy,PL,Fx_edgV_L,E_avg,E_edgV_L,E_edgH_T,E_edgH_B)
  REAL*8:: dr_ML_Calc
  REAL*8,INTENT(IN):: c, A, Delx, Dely
  REAL*8,INTENT(IN):: DC_xx, DL_xx, DT_xy, DB_xy, PL
  REAL*8,INTENT(IN):: Fx_edgV_L, E_avg, E_edgV_L, E_edgH_T, E_edgH_B

  dr_ML_Calc = A*Fx_edgV_L/2d0 + c*Dely*(DC_xx*E_avg - DL_xx*E_edgV_L) + c*Delx*(DT_xy*E_edgH_T - DB_xy*E_edgH_B)/2d0 - PL

END FUNCTION dr_ML_Calc

!==================================================================================================================================!
!
!==================================================================================================================================!
FUNCTION dr_MR_Calc(c,A,Delx,Dely,DR_xx,DC_xx,DT_xy,DB_xy,PR,Fx_edgV_R,E_edgV_R,E_avg,E_edgH_T,E_edgH_B)
  REAL*8:: dr_MR_Calc
  REAL*8,INTENT(IN):: c, A, Delx, Dely
  REAL*8,INTENT(IN):: DR_xx, DC_xx, DT_xy, DB_xy, PR
  REAL*8,INTENT(IN):: Fx_edgV_R, E_edgV_R, E_avg, E_edgH_T, E_edgH_B

  dr_MR_Calc = A*Fx_edgV_R/2d0 + c*Dely*(DR_xx*E_edgV_R - DC_xx*E_avg) + c*Delx*(DT_xy*E_edgH_T - DB_xy*E_edgH_B)/2d0 - PR

END FUNCTION dr_MR_Calc

!==================================================================================================================================!
!
!==================================================================================================================================!
FUNCTION dr_MB_Calc(c,A,Delx,Dely,DC_yy,DB_yy,DL_xy,DR_xy,PB,Fy_edgH_B,E_avg,E_edgH_B,E_edgV_R,E_edgV_L)
  REAL*8:: dr_MB_Calc
  REAL*8,INTENT(IN):: c, A, Delx, Dely
  REAL*8,INTENT(IN):: DC_yy, DB_yy, DL_xy, DR_xy, PB
  REAL*8,INTENT(IN):: Fy_edgH_B, E_avg, E_edgH_B, E_edgV_R, E_edgV_L

  dr_MB_Calc = A*Fy_edgH_B/2d0 + c*Delx*(DC_yy*E_avg - DB_yy*E_edgH_B) + c*Dely*(DR_xy*E_edgV_R - DL_xy*E_edgV_L)/2d0 - PB

END FUNCTION dr_MB_Calc

!==================================================================================================================================!
!
!==================================================================================================================================!
FUNCTION dr_MT_Calc(c,A,Delx,Dely,DT_yy,DC_yy,DL_xy,DR_xy,PT,Fy_edgH_T,E_edgH_T,E_avg,E_edgV_R,E_edgV_L)
  REAL*8:: dr_MT_Calc
  REAL*8,INTENT(IN):: c, A, Delx, Dely
  REAL*8,INTENT(IN):: DT_yy, DC_yy, DL_xy, DR_xy, PT
  REAL*8,INTENT(IN):: Fy_edgH_T, E_edgH_T, E_avg, E_edgV_R, E_edgV_L

  dr_MT_Calc = A*Fy_edgH_T/2d0 + c*Delx*(DT_yy*E_edgH_T - DC_yy*E_avg) + c*Dely*(DR_xy*E_edgV_R - DL_xy*E_edgV_L)/2d0 - PT

END FUNCTION dr_MT_Calc

END MODULE GLOQD_SOLVES
