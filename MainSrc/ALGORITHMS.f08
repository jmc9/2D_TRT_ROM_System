MODULE ALGORITHMS

  USE TRANSPORT_SOLVES
  USE UPDATES
  USE CONVERGENCE_CHECKS
  USE MLOQD_SOLVES
  USE OUTPUTS

  IMPLICIT NONE

CONTAINS

!============================================================================================================!
!
!============================================================================================================!
SUBROUTINE TRT_MLQD_ALGORITHM(Omega_x,Omega_y,quad_weight,Delx,Dely,A,Delt,Final_Time,c,cV,Kap0,Comp_Unit,Bg,KapE,&
  kapB,RT_Src,Temp,Temp_old,I_crn_old,I_avg,I_edgV,I_edgH,I_crn,Ic_edgV,Ic_edgH,RT_Residual,Hg_avg_xx,Hg_avg_xy,&
  Hg_avg_yy,Hg_edgV_xx,Hg_edgV_xy,Hg_edgH_yy,Hg_edgH_xy,HO_Eg_avg,HO_Eg_edgV,HO_Eg_edgH,HO_Fxg_edgV,HO_Fyg_edgH,&
  HO_E_avg,HO_E_edgV,HO_E_edgH,HO_Fx_edgV,HO_Fy_edgH,Temp_Times,HO_E_avg_Times,Nu_g,Conv_HO,Maxit_RTE,Conv_Type,&
  Start_Time,Threads,BC_Type,bcT_left,bcT_lower,bcT_right,bcT_upper,MGQD_Src,KapR,KapE_old,KapR_old,Theta,Eg_avg,&
  Eg_edgV,Eg_edgH,Fxg_edgV,Fyg_edgH,fg_avg_xx,fg_avg_xy,fg_avg_yy,fg_edgV_xx,fg_edgV_xy,fg_edgH_yy,fg_edgH_xy,&
  fg_avg_xx_old,fg_avg_xy_old,fg_avg_yy_old,fg_edgV_xx_old,fg_edgV_xy_old,fg_edgH_yy_old,fg_edgH_xy_old,Cg_L,Cg_B,&
  Cg_R,Cg_T,Eg_in_L,Eg_in_B,Eg_in_R,Eg_in_T,Fg_in_L,Fg_in_B,Fg_in_R,Fg_in_T,Eg_avg_old,Eg_edgV_old,Eg_edgH_old,&
  Fxg_edgV_old,Fyg_edgH_old,MGQD_Src_old,Res_Calc,MGQD_Residual,MGQD_E_avg,MGQD_E_edgV,MGQD_E_edgH,MGQD_Fx_edgV,&
  MGQD_Fy_edgH,MGQD_E_avg_Times,Maxit_MLOQD,run_type,Conv_LO,G_old,Pold_L,Pold_B,Pold_R,Pold_T)

  REAL*8,INTENT(IN):: Omega_x(:), Omega_y(:), quad_weight(:), Nu_g(:)
  REAL*8,INTENT(IN):: Delx(:), Dely(:), A(:,:), Delt
  REAL*8,INTENT(IN):: Final_Time, Start_Time
  REAL*8,INTENT(IN):: c, cV, Kap0
  REAL*8,INTENT(IN):: Comp_Unit, Conv_HO, Conv_LO
  REAL*8,INTENT(IN):: bcT_left, bcT_lower, bcT_right, bcT_upper
  INTEGER,INTENT(IN):: Conv_Type, Maxit_RTE, Threads, BC_Type(:), Maxit_MLOQD
  REAL*8,INTENT(IN):: Theta
  LOGICAL,INTENT(IN):: Res_Calc
  CHARACTER(100),INTENT(IN):: run_type

  REAL*8,INTENT(INOUT):: Bg(:,:,:), KapE(:,:,:), KapB(:,:,:), KapR(:,:,:)
  REAL*8,INTENT(INOUT):: KapE_old(:,:,:), KapR_old(:,:,:)
  REAL*8,INTENT(INOUT):: MGQD_Src(:,:,:), MGQD_Src_old(:,:,:), RT_Src(:,:,:,:)
  REAL*8,INTENT(INOUT):: Temp(:,:), Temp_old(:,:)

  REAL*8,INTENT(INOUT):: I_crn_old(:,:,:,:)
  REAL*8,INTENT(INOUT):: I_avg(:,:,:,:), I_edgV(:,:,:,:), I_edgH(:,:,:,:)
  REAL*8,INTENT(INOUT):: I_crn(:,:,:,:), Ic_edgV(:,:,:,:), Ic_edgH(:,:,:,:)
  REAL*8,INTENT(INOUT):: RT_Residual(:,:,:,:,:,:)

  REAL*8,INTENT(OUT):: Hg_avg_xx(:,:,:), Hg_avg_xy(:,:,:), Hg_avg_yy(:,:,:)
  REAL*8,INTENT(OUT):: Hg_edgV_xx(:,:,:), Hg_edgV_xy(:,:,:)
  REAL*8,INTENT(OUT):: Hg_edgH_yy(:,:,:), Hg_edgH_xy(:,:,:)
  REAL*8,INTENT(OUT):: HO_Eg_avg(:,:,:), HO_Eg_edgV(:,:,:), HO_Eg_edgH(:,:,:)
  REAL*8,INTENT(OUT):: HO_Fxg_edgV(:,:,:), HO_Fyg_edgH(:,:,:)
  REAL*8,INTENT(OUT):: HO_E_avg(:,:), HO_E_edgV(:,:), HO_E_edgH(:,:)
  REAL*8,INTENT(OUT):: HO_Fx_edgV(:,:), HO_Fy_edgH(:,:)

  REAL*8,INTENT(OUT):: Eg_avg(:,:,:), Eg_edgV(:,:,:), Eg_edgH(:,:,:)
  REAL*8,INTENT(OUT):: Fxg_edgV(:,:,:), Fyg_edgH(:,:,:)
  REAL*8,INTENT(OUT):: MGQD_Residual(:,:,:,:,:,:,:)
  REAL*8,INTENT(OUT):: MGQD_E_avg(:,:), MGQD_E_edgV(:,:), MGQD_E_edgH(:,:)
  REAL*8,INTENT(OUT):: MGQD_Fx_edgV(:,:), MGQD_Fy_edgH(:,:)

  REAL*8,INTENT(INOUT):: fg_avg_xx(:,:,:), fg_avg_xy(:,:,:), fg_avg_yy(:,:,:)
  REAL*8,INTENT(INOUT):: fg_edgV_xx(:,:,:), fg_edgV_xy(:,:,:)
  REAL*8,INTENT(INOUT):: fg_edgH_yy(:,:,:), fg_edgH_xy(:,:,:)
  REAL*8,INTENT(INOUT):: fg_avg_xx_old(:,:,:), fg_avg_xy_old(:,:,:), fg_avg_yy_old(:,:,:)
  REAL*8,INTENT(INOUT):: fg_edgV_xx_old(:,:,:), fg_edgV_xy_old(:,:,:)
  REAL*8,INTENT(INOUT):: fg_edgH_yy_old(:,:,:), fg_edgH_xy_old(:,:,:)
  REAL*8,INTENT(INOUT):: Cg_L(:,:), Cg_B(:,:), Cg_R(:,:), Cg_T(:,:)
  REAL*8,INTENT(INOUT):: Eg_in_L(:,:), Eg_in_B(:,:), Eg_in_R(:,:), Eg_in_T(:,:)
  REAL*8,INTENT(INOUT):: Fg_in_L(:,:), Fg_in_B(:,:), Fg_in_R(:,:), Fg_in_T(:,:)
  REAL*8,INTENT(INOUT):: Eg_avg_old(:,:,:), Eg_edgV_old(:,:,:), Eg_edgH_old(:,:,:)
  REAL*8,INTENT(INOUT):: Fxg_edgV_old(:,:,:), Fyg_edgH_old(:,:,:)
  REAL*8,INTENT(INOUT):: G_old(:,:,:), Pold_L(:,:,:), Pold_B(:,:,:), Pold_R(:,:,:), Pold_T(:,:,:)

  REAL*8,INTENT(OUT):: Temp_Times(:,:,:), HO_E_avg_Times(:,:,:), MGQD_E_avg_Times(:,:,:)

  REAL*8,ALLOCATABLE:: Temp_RTold(:,:), Temp_RTold2(:,:)
  REAL*8,ALLOCATABLE:: Temp_MGQDold(:,:), Temp_MGQDold2(:,:)
  REAL*8,ALLOCATABLE:: HO_E_avg_RTold(:,:), HO_E_avg_RTold2(:,:)
  REAL*8,ALLOCATABLE:: MGQD_E_avg_MGQDold(:,:), MGQD_E_avg_MGQDold2(:,:)
  REAL*8,ALLOCATABLE:: TR_Tnorm(:), TR_Enorm(:), TR_Trho(:), TR_Erho(:)
  REAL*8,ALLOCATABLE:: MGQD_Tnorm(:,:), MGQD_Enorm(:,:), MGQD_Trho(:,:), MGQD_Erho(:,:)
  REAL*8:: Time
  INTEGER,ALLOCATABLE:: MGQD_Its(:)
  INTEGER:: N_x, N_y, N_g
  INTEGER:: RT_Its, RT_start_Its, t
  LOGICAL:: RT_Conv, MGQD_conv, Tconv, Econv

  N_g = SIZE(Eg_avg,3)
  N_y = SIZE(Eg_avg,2)
  N_x = SIZE(Eg_avg,1)

  ALLOCATE(Temp_RTold(SIZE(TEMP,1),SIZE(TEMP,2)))
  ALLOCATE(Temp_RTold2(SIZE(TEMP,1),SIZE(TEMP,2)))
  ALLOCATE(Temp_MGQDold(SIZE(TEMP,1),SIZE(TEMP,2)))
  ALLOCATE(Temp_MGQDold2(SIZE(TEMP,1),SIZE(TEMP,2)))
  ALLOCATE(HO_E_avg_RTold(SIZE(HO_E_avg,1),SIZE(HO_E_avg,2)))
  ALLOCATE(HO_E_avg_RTold2(SIZE(HO_E_avg,1),SIZE(HO_E_avg,2)))
  ALLOCATE(MGQD_E_avg_MGQDold(SIZE(MGQD_E_avg,1),SIZE(MGQD_E_avg,2)))
  ALLOCATE(MGQD_E_avg_MGQDold2(SIZE(MGQD_E_avg,1),SIZE(MGQD_E_avg,2)))

  ALLOCATE(MGQD_Its(Maxit_RTE))
  ALLOCATE(TR_Tnorm(Maxit_RTE),TR_Enorm(Maxit_RTE),TR_Trho(Maxit_RTE),TR_Erho(Maxit_RTE))
  ALLOCATE(MGQD_Tnorm(Maxit_RTE,Maxit_MLOQD),MGQD_Enorm(Maxit_RTE,Maxit_MLOQD))
  ALLOCATE(MGQD_Trho(Maxit_RTE,Maxit_MLOQD),MGQD_Erho(Maxit_RTE,Maxit_MLOQD))

  IF ( run_type .EQ. 'tr_no_qd' ) THEN
    RT_start_Its = 1
  ELSE IF ( run_type .EQ. 'mlqd' ) THEN
    RT_start_Its = 2
  ELSE
    RT_start_Its = Maxit_RTE + 1
  END IF

  Time = Start_Time
  t = 0
  DO
    HO_E_avg_RTold2 = 0d0
    HO_E_avg_RTold = 0d0
    Temp_RTold2 = 0d0
    Temp_RTold = 0d0
    RT_Its = 0
    t = t + 1
    Time = Time + Delt
    RT_Conv = .FALSE.
    DO WHILE ((.NOT. RT_Conv).AND.(RT_Its .LT. Maxit_RTE))
      RT_Its = RT_Its + 1

      IF (RT_Its .GE. RT_start_Its) THEN

        CALL TRANSPORT_SCB(I_avg,I_edgV,I_edgH,I_crn,Ic_edgV,Ic_edgH,RT_Residual(:,:,:,:,:,RT_Its),Omega_x,Omega_y,Delx,Dely,A,&
          KapE,RT_Src,I_crn_old,c,Delt,Threads,Res_Calc)

        CALL COLLAPSE_INTENSITIES(Threads,I_avg,I_edgV,I_edgH,Omega_x,Omega_y,quad_weight,Comp_Unit,Hg_avg_xx,Hg_avg_xy,Hg_avg_yy,&
          Hg_edgV_xx,Hg_edgV_xy,Hg_edgH_yy,Hg_edgH_xy,HO_Eg_edgV,HO_Eg_edgH,HO_Eg_avg,HO_Fxg_edgV,HO_Fyg_edgH,HO_E_edgV,HO_E_edgH,&
          HO_E_avg,HO_Fx_edgV,HO_Fy_edgH)

        IF (MAXVAL(BC_Type) .GT. 0) THEN
          CALL RT_BC_UPDATE(BC_Type,bcT_left,bcT_lower,bcT_right,bcT_upper,Comp_Unit,Nu_g,I_edgV,I_edgH,Ic_edgV,Ic_edgH)
          CALL MGQD_In_Calc(Eg_in_L,Eg_in_B,Eg_in_R,Eg_in_T,Fg_in_L,Fg_in_B,Fg_in_R,Fg_in_T,I_edgV,I_edgH,Omega_x,Omega_y,&
            quad_weight,c,Comp_Unit)
        END IF

        CALL fg_Calc(fg_avg_xx,fg_avg_xy,fg_avg_yy,fg_edgV_xx,fg_edgV_xy,fg_edgH_yy,fg_edgH_xy,Hg_avg_xx,Hg_avg_xy,Hg_avg_yy,&
          Hg_edgV_xx,Hg_edgV_xy,Hg_edgH_yy,Hg_edgH_xy,HO_Eg_avg,HO_Eg_edgV,HO_Eg_edgH,c,Comp_Unit)

        CALL Cg_Calc(Cg_L,Cg_B,Cg_R,Cg_T,I_edgV,I_edgH,Omega_x,Omega_y,quad_weight,Comp_Unit)

      END IF

      IF ( run_type .EQ. 'tr_no_qd' ) THEN
        CALL MEB_SOLVE(Temp,HO_Eg_avg,Bg,KapE,Temp_old,Delt,cV,Kap0,Comp_Unit)
        CALL MATERIAL_UPDATE(RT_Src,MGQD_Src,KapE,KapB,KapR,Bg,Temp,Comp_Unit,Nu_g)

      ELSE
        MGQD_E_avg_MGQDold2 = 0d0
        MGQD_E_avg_MGQDold = HO_E_avg
        Temp_MGQDold2 = 0d0
        Temp_MGQDold = Temp
        MGQD_Its(RT_Its) = 0
        MGQD_conv = .FALSE.
        DO WHILE ((.NOT. MGQD_conv).AND.(MGQD_Its(RT_Its) .LT. Maxit_MLOQD))
          MGQD_Its(RT_Its) = MGQD_Its(RT_Its) + 1

          CALL MLOQD_FV(Eg_avg,Eg_edgV,Eg_edgH,Fxg_edgV,Fyg_edgH,fg_avg_xx,fg_avg_xy,fg_avg_yy,fg_edgV_xx,fg_edgV_xy,fg_edgH_yy,&
            fg_edgH_xy,Cg_L,Cg_B,Cg_R,Cg_T,Eg_in_L,Eg_in_B,Eg_in_R,Eg_in_T,Fg_in_L,Fg_in_B,Fg_in_R,Fg_in_T,MGQD_Src,KapE,KapR,&
            Delx,Dely,A,c,Delt,Theta,Threads,Res_Calc,MGQD_Residual(:,:,:,:,:,MGQD_Its(RT_Its),RT_Its),G_old,Pold_L,Pold_B,&
            Pold_R,Pold_T,Eg_avg_old,Eg_edgV_old,Eg_edgH_old,Fxg_edgV_old,Fyg_edgH_old)
          CALL COLLAPSE_MG_EF(Eg_avg,Eg_edgV,Eg_edgH,Fxg_edgV,Fyg_edgH,Threads,MGQD_E_avg,MGQD_E_edgV,MGQD_E_edgH,MGQD_Fx_edgV,&
            MGQD_Fy_edgH)

          CALL MEB_SOLVE(Temp,Eg_avg,Bg,KapE,Temp_old,Delt,cV,Kap0,Comp_Unit)
          CALL MATERIAL_UPDATE(RT_Src,MGQD_Src,KapE,KapB,KapR,Bg,Temp,Comp_Unit,Nu_g)

          Tconv = CONVERGENCE(Conv_Type,conv_lo,Temp,Temp_MGQDold,Temp_MGQDold2,MGQD_Its(RT_Its),&
            MGQD_Tnorm(RT_Its,MGQD_Its(RT_Its)),MGQD_Trho(RT_Its,MGQD_Its(RT_Its)))
          Econv = CONVERGENCE(Conv_Type,conv_lo,MGQD_E_avg,MGQD_E_avg_MGQDold,MGQD_E_avg_MGQDold2,MGQD_Its(RT_Its),&
            MGQD_Enorm(RT_Its,MGQD_Its(RT_Its)),MGQD_Erho(RT_Its,MGQD_Its(RT_Its)))
          MGQD_conv = Tconv.AND.Econv

          MGQD_E_avg_MGQDold2 = MGQD_E_avg_MGQDold
          MGQD_E_avg_MGQDold = MGQD_E_avg
          Temp_MGQDold2 = Temp_MGQDold
          Temp_MGQDold = Temp

        END DO

      END IF

      Tconv = CONVERGENCE(Conv_Type,conv_ho,Temp,Temp_RTold,Temp_RTold2,RT_Its,TR_Tnorm(RT_Its),TR_Trho(RT_Its))
      Econv = CONVERGENCE(Conv_Type,conv_ho,HO_E_avg,HO_E_avg_RTold,HO_E_avg_RTold2,RT_Its,TR_Enorm(RT_Its),TR_Erho(RT_Its))
      RT_conv = Tconv.AND.Econv

      HO_E_avg_RTold2 = HO_E_avg_RTold
      HO_E_avg_RTold = HO_E_avg
      Temp_RTold2 = Temp_RTold
      Temp_RTold = Temp

      write(*,*) RT_Its, MGQD_Its(RT_Its), TR_Tnorm(RT_Its), TR_Enorm(RT_Its)
      IF (RT_Its .GE. Maxit_RTE) EXIT

    END DO
    write(*,*) Time, RT_Its, TR_Tnorm(RT_Its), TR_Enorm(RT_Its)

    CALL ITERATIVE_OUTPUT(103,t,RT_its,MGQD_its,0,Time,c,TR_Tnorm,TR_Enorm,TR_Trho,TR_Erho,MGQD_Tnorm,&
      MGQD_Enorm,MGQD_Trho,MGQD_Erho,RT_Residual,MGQD_Residual)

    Temp_Times(:,:,t) = Temp
    HO_E_avg_Times(:,:,t) = HO_E_avg
    MGQD_E_avg_Times(:,:,t) = MGQD_E_avg

    Temp_Old = Temp
    I_crn_old = I_crn

    fg_avg_xx_old = fg_avg_xx
    fg_avg_xy_old = fg_avg_xy
    fg_avg_yy_old = fg_avg_yy
    fg_edgV_xx_old = fg_edgV_xx
    fg_edgV_xy_old = fg_edgV_xy
    fg_edgH_yy_old = fg_edgH_yy
    fg_edgH_xy_old = fg_edgH_xy
    MGQD_Src_old = MGQD_Src
    Eg_avg_old = Eg_avg
    Eg_edgV_old = Eg_edgV
    Eg_edgH_old = Eg_edgH
    Fxg_edgV_old = Fxg_edgV
    Fyg_edgH_old = Fyg_edgH
    KapE_old = KapE
    KapR_old = KapR

    CALL OLD_MGQD_COEFS(Eg_avg_old,Eg_edgV_old,Eg_edgH_old,Fxg_edgV_old,Fyg_edgH_old,fg_avg_xx_old,fg_avg_yy_old,&
      fg_edgV_xx_old,fg_edgV_xy_old,fg_edgH_yy_old,fg_edgH_xy_old,KapE_old,KapR_old,MGQD_Src_old,Delx,Dely,A,c,Delt,&
      Theta,G_old,Pold_L,Pold_B,Pold_R,Pold_T)

    IF (Time .GE. Final_Time) EXIT
  END DO

END SUBROUTINE TRT_MLQD_ALGORITHM

!============================================================================================================!
!
!============================================================================================================!

!============================================================================================================!
!
!============================================================================================================!

END MODULE ALGORITHMS
