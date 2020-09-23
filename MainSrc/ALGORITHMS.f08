MODULE ALGORITHMS

  USE TRANSPORT_SOLVES
  USE UPDATES
  USE CONVERGENCE_CHECKS
  USE MLOQD_SOLVES
  ! USE GLOQD_SOLVES
  USE OUTPUTS
  USE INITIALIZERS
  USE netcdf
  USE NCDF_IO

  IMPLICIT NONE

CONTAINS

!==================================================================================================================================!
!
!==================================================================================================================================!
SUBROUTINE TRT_MLQD_ALGORITHM(Omega_x,Omega_y,quad_weight,Nu_g,Delx,Dely,Delt,tlen,Theta,Start_Time,c,cV,h,pi,Kap0,&
  erg,Comp_Unit,Conv_ho,Conv_lo,Conv_gr1,Conv_gr2,bcT_left,bcT_bottom,bcT_right,bcT_top,Tini,chi,line_src,E_Bound_Low,&
  T_Bound_Low,database_gen,use_grey,Conv_Type,Maxit_RTE,Threads,BC_Type,Maxit_MLOQD,Maxit_GLOQD,N_x,N_y,N_m,N_g,&
  N_t,Res_Calc,Use_Line_Search,Use_Safety_Search,run_type,kapE_dT_flag,outID,N_x_ID,N_y_ID,N_m_ID,N_g_ID,N_t_ID,&
  N_edgV_ID,N_edgH_ID,N_xc_ID,N_yc_ID,Quads_ID,RT_Its_ID,MGQD_Its_ID,GQD_Its_ID,Norm_Types_ID,MGQD_ResTypes_ID,&
  Boundaries_ID)

  !---------------Solution Parameters----------------!
  REAL*8,INTENT(IN):: Omega_x(:), Omega_y(:), quad_weight(:), Nu_g(:)
  REAL*8,INTENT(IN):: Delx(:), Dely(:), Delt, Theta, tlen
  REAL*8,INTENT(IN):: Start_Time
  REAL*8,INTENT(IN):: c, cV, h, pi, Kap0, erg
  REAL*8,INTENT(IN):: Comp_Unit, Conv_ho, Conv_lo, Conv_gr1, Conv_gr2
  REAL*8,INTENT(IN):: bcT_left, bcT_bottom, bcT_right, bcT_top, Tini
  REAL*8,INTENT(IN):: chi, line_src, E_Bound_Low, T_Bound_Low
  INTEGER,INTENT(IN):: N_x, N_y, N_m, N_g, N_t
  INTEGER,INTENT(IN):: database_gen, use_grey, Conv_Type, Threads, BC_Type(:), Maxit_RTE, Maxit_MLOQD, Maxit_GLOQD
  LOGICAL,INTENT(IN):: Res_Calc, Use_Line_Search, Use_Safety_Search, kapE_dT_flag
  CHARACTER(*),INTENT(IN):: run_type

  !----------------Output File ID's------------------!
  INTEGER,INTENT(IN):: outID
  INTEGER,INTENT(IN):: N_x_ID, N_y_ID, N_m_ID, N_g_ID, N_t_ID, N_edgV_ID, N_edgH_ID, N_xc_ID, N_yc_ID, Quads_ID
  INTEGER,INTENT(IN):: RT_Its_ID, MGQD_Its_ID, GQD_Its_ID, Norm_Types_ID, MGQD_ResTypes_ID, Boundaries_ID

  !---------------Material Properties----------------!
  REAL*8,ALLOCATABLE:: Bg(:,:,:), KapE(:,:,:), KapB(:,:,:), KapR(:,:,:), A(:,:)
  REAL*8,ALLOCATABLE:: KapE_old(:,:,:), KapR_old(:,:,:)
  REAL*8,ALLOCATABLE:: MGQD_Src(:,:,:), MGQD_Src_old(:,:,:), RT_Src(:,:,:,:)

  !--------------Radiation Intensities---------------!
  REAL*8,ALLOCATABLE:: I_crn_old(:,:,:,:)
  REAL*8,ALLOCATABLE:: I_avg(:,:,:,:), I_edgV(:,:,:,:), I_edgH(:,:,:,:)
  REAL*8,ALLOCATABLE:: I_crn(:,:,:,:), Ic_edgV(:,:,:,:), Ic_edgH(:,:,:,:)

  !--------------High-Order Quantities---------------!
  REAL*8,ALLOCATABLE:: Hg_avg_xx(:,:,:), Hg_avg_yy(:,:,:)
  REAL*8,ALLOCATABLE:: Hg_edgV_xx(:,:,:), Hg_edgV_xy(:,:,:)
  REAL*8,ALLOCATABLE:: Hg_edgH_yy(:,:,:), Hg_edgH_xy(:,:,:)
  REAL*8,ALLOCATABLE:: HO_Eg_avg(:,:,:), HO_Eg_edgV(:,:,:), HO_Eg_edgH(:,:,:)
  REAL*8,ALLOCATABLE:: HO_Fxg_edgV(:,:,:), HO_Fyg_edgH(:,:,:)
  REAL*8,ALLOCATABLE:: HO_E_avg(:,:), HO_E_edgV(:,:), HO_E_edgH(:,:)
  REAL*8,ALLOCATABLE:: HO_Fx_edgV(:,:), HO_Fy_edgH(:,:)

  !------------------MGQD Solution-------------------!
  REAL*8,ALLOCATABLE:: Eg_avg(:,:,:), Eg_edgV(:,:,:), Eg_edgH(:,:,:)
  REAL*8,ALLOCATABLE:: Fxg_edgV(:,:,:), Fyg_edgH(:,:,:)
  REAL*8,ALLOCATABLE:: MGQD_E_avg(:,:), MGQD_E_edgV(:,:), MGQD_E_edgH(:,:)
  REAL*8,ALLOCATABLE:: MGQD_Fx_edgV(:,:), MGQD_Fy_edgH(:,:)

  !-----------------MGQD Parameters------------------!
  REAL*8,ALLOCATABLE:: fg_avg_xx(:,:,:), fg_avg_yy(:,:,:)
  REAL*8,ALLOCATABLE:: fg_edgV_xx(:,:,:), fg_edgV_xy(:,:,:)
  REAL*8,ALLOCATABLE:: fg_edgH_yy(:,:,:), fg_edgH_xy(:,:,:)
  REAL*8,ALLOCATABLE:: fg_avg_xx_old(:,:,:), fg_avg_yy_old(:,:,:)
  REAL*8,ALLOCATABLE:: fg_edgV_xx_old(:,:,:), fg_edgV_xy_old(:,:,:)
  REAL*8,ALLOCATABLE:: fg_edgH_yy_old(:,:,:), fg_edgH_xy_old(:,:,:)
  REAL*8,ALLOCATABLE:: Cg_L(:,:), Cg_B(:,:), Cg_R(:,:), Cg_T(:,:)
  REAL*8,ALLOCATABLE:: Eg_in_L(:,:), Eg_in_B(:,:), Eg_in_R(:,:), Eg_in_T(:,:)
  REAL*8,ALLOCATABLE:: Fg_in_L(:,:), Fg_in_B(:,:), Fg_in_R(:,:), Fg_in_T(:,:)
  REAL*8,ALLOCATABLE:: Eg_avg_old(:,:,:), Eg_edgV_old(:,:,:), Eg_edgH_old(:,:,:)
  REAL*8,ALLOCATABLE:: Fxg_edgV_old(:,:,:), Fyg_edgH_old(:,:,:)
  REAL*8,ALLOCATABLE:: G_old(:,:,:), Pold_L(:,:,:), Pold_B(:,:,:), Pold_R(:,:,:), Pold_T(:,:,:)

  !-------------------EGP Solution-------------------!
  REAL*8,ALLOCATABLE:: Temp(:,:), Temp_old(:,:)
  REAL*8,ALLOCATABLE:: E_avg(:,:), E_edgV(:,:), E_edgH(:,:)
  REAL*8,ALLOCATABLE:: Fx_edgV(:,:), Fy_edgH(:,:)
  REAL*8,ALLOCATABLE:: E_avg_old(:,:), KapE_bar_old(:,:),  GQD_Src_old(:,:)
  REAL*8,ALLOCATABLE:: Fx_edgV_old(:,:), Fy_edgH_old(:,:)
  REAL*8,ALLOCATABLE:: Gold_Hat(:,:), Rhat_old(:,:)
  REAL*8,ALLOCATABLE:: KapE_Bar(:,:), KapE_Bar_dT(:,:), GQD_Src(:,:)
  REAL*8,ALLOCATABLE:: DC_xx(:,:), DL_xx(:,:), DR_xx(:,:)
  REAL*8,ALLOCATABLE:: DC_yy(:,:), DB_yy(:,:), DT_yy(:,:)
  REAL*8,ALLOCATABLE:: DL_xy(:,:), DR_xy(:,:), DB_xy(:,:), DT_xy(:,:)
  REAL*8,ALLOCATABLE:: PL(:,:), PR(:,:), PB(:,:), PT(:,:)
  REAL*8,ALLOCATABLE:: Cb_L(:), Cb_B(:), Cb_R(:), Cb_T(:)
  REAL*8,ALLOCATABLE:: E_in_L(:), E_in_B(:), E_in_R(:), E_in_T(:)
  REAL*8,ALLOCATABLE:: F_in_L(:), F_in_B(:), F_in_R(:), F_in_T(:)
  REAL*8,ALLOCATABLE:: KapE_bar_MGQDold(:,:)

  !--------------------------------------------------!
  REAL*8,ALLOCATABLE:: RT_Residual(:,:,:), MGQD_Residual(:,:,:), MGQD_BC_Residual(:,:), Deltas(:,:), dres(:,:)
  INTEGER,ALLOCATABLE:: RT_ResLoc_x(:,:,:), RT_ResLoc_y(:,:,:)
  REAL*8,ALLOCATABLE:: Temp_RTold(:,:), Temp_RTold2(:,:)
  REAL*8,ALLOCATABLE:: Temp_MGQDold(:,:), Temp_MGQDold2(:,:)
  REAL*8,ALLOCATABLE:: HO_E_avg_RTold(:,:), HO_E_avg_RTold2(:,:)
  REAL*8,ALLOCATABLE:: MGQD_E_avg_MGQDold(:,:), MGQD_E_avg_MGQDold2(:,:)
  REAL*8:: TR_Tnorm, TR_Enorm, TR_Trho, TR_Erho
  REAL*8:: MGQD_Tnorm, MGQD_Enorm, MGQD_Trho, MGQD_Erho
  REAL*8:: Time
  INTEGER:: MGQD_Its, EGP_Its, Status
  INTEGER:: RT_Its, RT_start_Its, t
  LOGICAL:: RT_Conv, MGQD_conv, Tconv, Econv

  !-------------------Variable ID's------------------!
  INTEGER:: Temp_ID, E_avg_ID, E_edgV_ID, E_edgH_ID, MGQD_E_avg_ID, MGQD_E_edgV_ID, MGQD_E_edgH_ID, HO_E_avg_ID
  INTEGER:: HO_E_edgV_ID, HO_E_edgH_ID, Fx_edgV_ID, Fy_edgH_ID, MGQD_Fx_edgV_ID, MGQD_Fy_edgH_ID, HO_Fx_edgV_ID
  INTEGER:: HO_Fy_edgH_ID, Eg_avg_ID, Eg_edgV_ID, Eg_edgH_ID, HO_Eg_avg_ID, HO_Eg_edgV_ID, HO_Eg_edgH_ID
  INTEGER:: Fxg_edgV_ID, Fyg_edgH_ID, HO_Fxg_edgV_ID, HO_Fyg_edgH_ID, I_avg_ID, I_edgV_ID, I_edgH_ID
  INTEGER:: KapE_Bar_ID, KapB_ID, KapE_ID, KapR_ID, Bg_ID
  INTEGER:: RT_Residual_ID, MGQD_Residual_ID, MGQD_BC_Residual_ID
  INTEGER:: Del_T_ID, Del_E_avg_ID, Del_E_edgV_ID, Del_E_edgH_ID, Del_Fx_edgV_ID, Del_Fy_edgH_ID
  INTEGER:: RT_ItCount_ID, MGQD_ItCount_ID, GQD_ItCount_ID
  INTEGER:: RT_Tnorm_ID, RT_Enorm_ID, MGQD_Tnorm_ID, MGQD_Enorm_ID
  INTEGER:: RT_Trho_ID, RT_Erho_ID, MGQD_Trho_ID, MGQD_Erho_ID
  INTEGER:: Cg_L_ID, Cg_B_ID, Cg_R_ID, Cg_T_ID, Eg_in_L_ID, Eg_in_B_ID, Eg_in_R_ID, Eg_in_T_ID
  INTEGER:: Fg_in_L_ID, Fg_in_B_ID, Fg_in_R_ID, Fg_in_T_ID
  INTEGER:: Cb_L_ID, Cb_B_ID, Cb_R_ID, Cb_T_ID, E_in_L_ID, E_in_B_ID, E_in_R_ID, E_in_T_ID
  INTEGER:: F_in_L_ID, F_in_B_ID, F_in_R_ID, F_in_T_ID
  INTEGER:: fg_avg_xx_ID, fg_avg_yy_ID, fg_edgV_xx_ID, fg_edgV_xy_ID, fg_edgH_yy_ID, fg_edgH_xy_ID
  INTEGER:: DC_xx_ID, DL_xx_ID, DR_xx_ID, DC_yy_ID, DB_yy_ID, DT_yy_ID, DL_xy_ID, DB_xy_ID, DR_xy_ID, DT_xy_ID
  INTEGER:: G_old_ID, Pold_L_ID, Pold_B_ID, Pold_R_ID, Pold_T_ID
  INTEGER:: Gold_hat_ID, Rhat_old_ID, PL_ID, PB_ID, PR_ID, PT_ID
  INTEGER:: dr_T_ID, dr_B_ID, dr_ML_ID, dr_MB_ID, dr_MR_ID, dr_MT_ID

  !===========================================================================!
  !                                                                           !
  !     INITIALIZING ARRAYS                                                   !
  !                                                                           !
  !===========================================================================!
  CALL MISC_INIT(Delx,Dely,A)
  CALL RT_INIT(I_avg,I_edgV,I_edgH,I_crn,I_crn_old,Ic_edgV,Ic_edgH,Hg_avg_xx,Hg_avg_yy,Hg_edgV_xx,&
    Hg_edgV_xy,Hg_edgH_yy,Hg_edgH_xy,HO_Eg_avg,HO_Eg_edgV,HO_Eg_edgH,HO_Fxg_edgV,HO_Fyg_edgH,HO_E_avg,&
    HO_E_edgV,HO_E_edgH,HO_Fx_edgV,HO_Fy_edgH,N_y,N_x,N_m,N_g,Tini,comp_unit,nu_g,bcT_left,bcT_right,&
    bcT_top,bcT_bottom,BC_Type,pi,c)
  CALL MGQD_INIT(Eg_avg,Eg_edgV,Eg_edgH,Fxg_edgV,Fyg_edgH,Eg_avg_old,Eg_edgV_old,Eg_edgH_old,Fxg_edgV_old,Fyg_edgH_old,&
    fg_avg_xx,fg_avg_yy,fg_edgV_xx,fg_edgV_xy,fg_edgH_yy,fg_edgH_xy,fg_avg_xx_old,fg_avg_yy_old,fg_edgV_xx_old,fg_edgV_xy_old,&
    fg_edgH_yy_old,fg_edgH_xy_old,Cg_L,Cg_B,Cg_R,Cg_T,Eg_in_L,Eg_in_B,Eg_in_R,Eg_in_T,Fg_in_L,Fg_in_B,Fg_in_R,Fg_in_T,I_edgV,&
    I_edgH,Omega_x,Omega_y,quad_weight,c,Comp_Unit,N_y,N_x,N_g,MGQD_E_avg,MGQD_E_edgV,MGQD_E_edgH,MGQD_Fx_edgV,MGQD_Fy_edgH,&
    G_old,Pold_L,Pold_B,Pold_R,Pold_T,BC_Type,Tini,nu_g,pi,Threads)
  CALL GQD_INIT(E_avg,E_edgV,E_edgH,Fx_edgV,Fy_edgH,E_avg_old,KapE_bar_old,GQD_Src_old,Gold_Hat,Rhat_old,&
    KapE_Bar,KapE_Bar_dT,GQD_Src,DC_xx,DL_xx,DR_xx,DC_yy,DB_yy,DT_yy,DL_xy,DB_xy,DR_xy,DT_xy,PL,PB,PR,PT,Cb_L,&
    Cb_B,Cb_R,Cb_T,E_in_L,E_in_B,E_in_R,E_in_T,F_in_L,F_in_B,F_in_R,F_in_T,KapE_bar_MGQDold,Eg_in_L,Eg_in_B,&
    Eg_in_R,Eg_in_T,Fg_in_L,Fg_in_B,Fg_in_R,Fg_in_T,Tini,nu_g,c,Comp_Unit,pi,N_y,N_x,N_g,BC_Type)
  CALL TEMP_INIT(Temp,RT_Src,MGQD_Src,MGQD_Src_old,KapE,KapB,KapR,KapE_old,KapR_old,Bg,N_y,N_x,N_m,N_g,Tini,&
    Comp_Unit,Nu_g,Temp_Old,Threads)

  ALLOCATE(Temp_RTold(SIZE(TEMP,1),SIZE(TEMP,2)))
  ALLOCATE(Temp_RTold2(SIZE(TEMP,1),SIZE(TEMP,2)))
  ALLOCATE(Temp_MGQDold(SIZE(TEMP,1),SIZE(TEMP,2)))
  ALLOCATE(Temp_MGQDold2(SIZE(TEMP,1),SIZE(TEMP,2)))
  ALLOCATE(HO_E_avg_RTold(SIZE(HO_E_avg,1),SIZE(HO_E_avg,2)))
  ALLOCATE(HO_E_avg_RTold2(SIZE(HO_E_avg,1),SIZE(HO_E_avg,2)))
  ALLOCATE(MGQD_E_avg_MGQDold(SIZE(MGQD_E_avg,1),SIZE(MGQD_E_avg,2)))
  ALLOCATE(MGQD_E_avg_MGQDold2(SIZE(MGQD_E_avg,1),SIZE(MGQD_E_avg,2)))

  IF (Res_Calc) THEN
    ALLOCATE(RT_Residual(4,N_g,5),RT_ResLoc_x(4,N_g,5),RT_ResLoc_y(4,N_g,5))
    ALLOCATE(MGQD_Residual(N_g,5,5),MGQD_BC_Residual(N_g,4))
  END IF

  !===========================================================================!
  !                                                                           !
  !     INITIALIZING OUTPUT FILE                                              !
  !                                                                           !
  !===========================================================================!
  CALL OUTFILE_VARDEFS(outID,Res_Calc,N_x_ID,N_y_ID,N_m_ID,N_g_ID,N_t_ID,N_edgV_ID,N_edgH_ID,N_xc_ID,N_yc_ID,Quads_ID,&
    RT_Its_ID,MGQD_Its_ID,GQD_Its_ID,Norm_Types_ID,MGQD_ResTypes_ID,Boundaries_ID,Temp_ID,E_avg_ID,E_edgV_ID,E_edgH_ID,&
    MGQD_E_avg_ID,MGQD_E_edgV_ID,MGQD_E_edgH_ID,HO_E_avg_ID,HO_E_edgV_ID,HO_E_edgH_ID,Fx_edgV_ID,Fy_edgH_ID,&
    MGQD_Fx_edgV_ID,MGQD_Fy_edgH_ID,HO_Fx_edgV_ID,HO_Fy_edgH_ID,Eg_avg_ID,Eg_edgV_ID,Eg_edgH_ID,HO_Eg_avg_ID,&
    HO_Eg_edgV_ID,HO_Eg_edgH_ID,Fxg_edgV_ID,Fyg_edgH_ID,HO_Fxg_edgV_ID,HO_Fyg_edgH_ID,I_avg_ID,I_edgV_ID,I_edgH_ID,&
    KapE_Bar_ID,KapB_ID,KapE_ID,KapR_ID,Bg_ID,RT_Residual_ID,MGQD_Residual_ID,MGQD_BC_Residual_ID,Del_T_ID,Del_E_avg_ID,&
    Del_E_edgV_ID,Del_E_edgH_ID,Del_Fx_edgV_ID,Del_Fy_edgH_ID,RT_ItCount_ID,MGQD_ItCount_ID,GQD_ItCount_ID,&
    RT_Tnorm_ID,RT_Enorm_ID,MGQD_Tnorm_ID,MGQD_Enorm_ID,RT_Trho_ID,RT_Erho_ID,&
    MGQD_Trho_ID,MGQD_Erho_ID,Cg_L_ID,Cg_B_ID,Cg_R_ID,Cg_T_ID,Eg_in_L_ID,Eg_in_B_ID,&
    Eg_in_R_ID,Eg_in_T_ID,Fg_in_L_ID,Fg_in_B_ID,Fg_in_R_ID,Fg_in_T_ID,Cb_L_ID,Cb_B_ID,Cb_R_ID,Cb_T_ID,E_in_L_ID,&
    E_in_B_ID,E_in_R_ID,E_in_T_ID,F_in_L_ID,F_in_B_ID,F_in_R_ID,F_in_T_ID,fg_avg_xx_ID,fg_avg_yy_ID,fg_edgV_xx_ID,&
    fg_edgV_xy_ID,fg_edgH_yy_ID,fg_edgH_xy_ID,DC_xx_ID,DL_xx_ID,DR_xx_ID,DC_yy_ID,DB_yy_ID,DT_yy_ID,DL_xy_ID,DB_xy_ID,&
    DR_xy_ID,DT_xy_ID,G_old_ID,Pold_L_ID,Pold_B_ID,Pold_R_ID,Pold_T_ID,Gold_hat_ID,Rhat_old_ID,PL_ID,&
    PB_ID,PR_ID,PT_ID,dr_T_ID,dr_B_ID,dr_ML_ID,dr_MB_ID,dr_MR_ID,dr_MT_ID)

  IF ( run_type .EQ. 'tr_no_qd' ) THEN
    RT_start_Its = 1
  ELSE IF ( run_type .EQ. 'mlqd' ) THEN
    RT_start_Its = 2
  ELSE
    RT_start_Its = Maxit_RTE + 1
  END IF

  !===========================================================================!
  !                                                                           !
  !     PROBLEM SOLVE (BEGIN TIME STEP LOOP)                                  !
  !                                                                           !
  !===========================================================================!
  Time = Start_Time
  t = 0
  DO

    !setting time and t to next time step
    t = t + 1
    Time = Time + Delt

    !moving last time step data to _old arrays
    Temp_Old = Temp
    I_crn_old = I_crn
    fg_avg_xx_old = fg_avg_xx
    fg_avg_yy_old = fg_avg_yy
    fg_edgV_xx_old = fg_edgV_xx
    fg_edgV_xy_old = fg_edgV_xy
    fg_edgH_yy_old = fg_edgH_yy
    fg_edgH_xy_old = fg_edgH_xy
    MGQD_Src_old = MGQD_Src
    GQD_Src_old = GQD_Src
    E_avg_old = E_avg
    Eg_avg_old = Eg_avg
    Eg_edgV_old = Eg_edgV
    Eg_edgH_old = Eg_edgH
    Fxg_edgV_old = Fxg_edgV
    Fyg_edgH_old = Fyg_edgH
    KapE_Bar_old = KapE_Bar
    KapE_old = KapE
    KapR_old = KapR

    !calculating coefficients for MGQD solve that only depend on the last time step solution
    CALL OLD_MGQD_COEFS(Eg_avg_old,Eg_edgV_old,Eg_edgH_old,Fxg_edgV_old,Fyg_edgH_old,fg_avg_xx_old,fg_avg_yy_old,&
      fg_edgV_xx_old,fg_edgV_xy_old,fg_edgH_yy_old,fg_edgH_xy_old,KapE_old,KapR_old,MGQD_Src_old,Delx,Dely,A,c,Delt,&
      Theta,G_old,Pold_L,Pold_B,Pold_R,Pold_T,Threads)
    CALL OLD_GREY_COEFS(c,Delt,Theta,cv,A,G_old,Temp_old,KapE_bar_old,E_avg_old,GQD_Src_old,Gold_Hat,Rhat_old)

    !===========================================================================!
    !                                                                           !
    !     OUTER ITERATION LOOP (RTE ITERATIONS)                                 !
    !                                                                           !
    !===========================================================================!
    HO_E_avg_RTold2 = 0d0
    HO_E_avg_RTold = 0d0
    Temp_RTold2 = 0d0
    Temp_RTold = 0d0
    RT_Its = 0
    MGQD_Its = 0
    RT_Conv = .FALSE.
    DO WHILE ((.NOT. RT_Conv).AND.(RT_Its .LT. Maxit_RTE))
      RT_Its = RT_Its + 1

      !the RTE is not necessarily solved all the time
      !only solving the RTE on and after the 2nd outer iteration for the mlqd algorithm
      !for methods that don't use the RTE it is never solved
      IF (RT_Its .GE. RT_start_Its) THEN

        !solving the RTE
        CALL TRANSPORT_SCB(I_avg,I_edgV,I_edgH,I_crn,Ic_edgV,Ic_edgH,Omega_x,Omega_y,Delx,Dely,A,&
          KapE,RT_Src,I_crn_old,c,Delt,Threads,RT_Residual,RT_ResLoc_x,RT_ResLoc_y,Res_Calc)

        !writing norms of the RTE residual to the output file
        IF (Res_Calc) THEN
          Status = nf90_put_var(outID,RT_Residual_ID,RT_Residual,(/1,1,1,RT_Its,t/),(/4,N_g,5,1,1/))
          CALL HANDLE_ERR(Status)
        END IF

        !calculating low-order quantities from the high-order intensities
        CALL COLLAPSE_INTENSITIES(Threads,I_avg,I_edgV,I_edgH,Omega_x,Omega_y,quad_weight,Comp_Unit,Hg_avg_xx,Hg_avg_yy,&
          Hg_edgV_xx,Hg_edgV_xy,Hg_edgH_yy,Hg_edgH_xy,HO_Eg_edgV,HO_Eg_edgH,HO_Eg_avg,HO_Fxg_edgV,HO_Fyg_edgH,HO_E_edgV,&
          HO_E_edgH,HO_E_avg,HO_Fx_edgV,HO_Fy_edgH)

        !if ANY boundary is non-vaccuum (reflective) then the boundary intensities must be updated
        IF (MAXVAL(BC_Type) .GT. 0) THEN
          !updating 'incoming' radiation intensities
          CALL RT_BC_UPDATE(BC_Type,bcT_left,bcT_bottom,bcT_right,bcT_top,Comp_Unit,Nu_g,I_edgV,I_edgH,Ic_edgV,Ic_edgH)
          !updating 'incoming' energy densities and fluxes accordingly
          CALL MGQD_In_Calc(Eg_in_L,Eg_in_B,Eg_in_R,Eg_in_T,Fg_in_L,Fg_in_B,Fg_in_R,Fg_in_T,I_edgV,I_edgH,Omega_x,Omega_y,&
            quad_weight,c,Comp_Unit,BC_Type,Threads)

          IF (use_grey .EQ. 1) CALL GQD_In_Calc(E_in_L,E_in_B,E_in_R,E_in_T,F_in_L,F_in_B,F_in_R,F_in_T,Eg_in_L,&
            Eg_in_B,Eg_in_R,Eg_in_T,Fg_in_L,Fg_in_B,Fg_in_R,Fg_in_T)
        END IF

        !calculating QD factors from the low-order quantities calculated from RTE solution
        CALL fg_Calc(fg_avg_xx,fg_avg_yy,fg_edgV_xx,fg_edgV_xy,fg_edgH_yy,fg_edgH_xy,Hg_avg_xx,Hg_avg_yy,&
          Hg_edgV_xx,Hg_edgV_xy,Hg_edgH_yy,Hg_edgH_xy,HO_Eg_avg,HO_Eg_edgV,HO_Eg_edgH,c,Threads)

        !calculating new MGQD boundary factors with the current iterate's intensities
        CALL Cg_Calc(Cg_L,Cg_B,Cg_R,Cg_T,I_edgV,I_edgH,Omega_x,Omega_y,quad_weight,Comp_Unit,BC_Type,Threads)

      END IF

      IF ( run_type .EQ. 'tr_no_qd' ) THEN
        !solve the MEB equation with the MGQD solution to find new Temp
        CALL MEB_SOLVE(Temp,HO_Eg_avg,Bg,KapE,Temp_old,Delt,cV,Kap0,Comp_Unit)
        !update material properties with new Temp
        CALL MATERIAL_UPDATE(RT_Src,MGQD_Src,KapE,KapB,KapR,Bg,Temp,Comp_Unit,Nu_g,Threads)

      ELSE
        !===========================================================================!
        !                                                                           !
        !     MGQD ITERATION LOOP                                                   !
        !                                                                           !
        !===========================================================================!
        MGQD_E_avg_MGQDold2 = 0d0
        MGQD_E_avg_MGQDold = HO_E_avg
        Temp_MGQDold2 = 0d0
        Temp_MGQDold = Temp
        MGQD_Its = 0
        MGQD_conv = .FALSE.
        DO WHILE ((.NOT. MGQD_conv).AND.(MGQD_Its .LT. Maxit_MLOQD))
          MGQD_Its = MGQD_Its + 1

          !solve the MGQD linear system
          CALL MLOQD_FV(Eg_avg,Eg_edgV,Eg_edgH,Fxg_edgV,Fyg_edgH,fg_avg_xx,fg_avg_yy,fg_edgV_xx,fg_edgV_xy,fg_edgH_yy,&
            fg_edgH_xy,Cg_L,Cg_B,Cg_R,Cg_T,Eg_in_L,Eg_in_B,Eg_in_R,Eg_in_T,Fg_in_L,Fg_in_B,Fg_in_R,Fg_in_T,MGQD_Src,KapE,&
            KapR,Delx,Dely,A,c,Delt,Theta,Threads,Res_Calc,MGQD_Residual,MGQD_BC_Residual,G_old,Pold_L,Pold_B,Pold_R,&
            Pold_T,Eg_avg_old,Fxg_edgV_old,Fyg_edgH_old)

          !calculate grey solution from MGQD solution
          CALL COLLAPSE_MG_EF(Eg_avg,Eg_edgV,Eg_edgH,Fxg_edgV,Fyg_edgH,Threads,MGQD_E_avg,MGQD_E_edgV,MGQD_E_edgH,MGQD_Fx_edgV,&
            MGQD_Fy_edgH)

          !writing MGQD residuals to output file
          IF (Res_Calc) THEN
            Status = nf90_put_var(outID,MGQD_BC_Residual_ID,MGQD_BC_Residual,(/1,1,MGQD_Its,RT_Its,t/),(/N_g,4,1,1,1/))
            CALL HANDLE_ERR(Status)
            Status = nf90_put_var(outID,MGQD_Residual_ID,MGQD_Residual,(/1,1,1,MGQD_Its,RT_Its,t/),(/N_g,5,5,1,1,1/))
            CALL HANDLE_ERR(Status)
          END IF

          IF (use_grey .EQ. 1) THEN
            !===========================================================================!
            !                                                                           !
            !     EFFECTIVE GREY PROBLEM                                                !
            !                                                                           !
            !===========================================================================!
            CALL Cbar_Calc(Cb_L,Cb_B,Cb_R,Cb_T,Cg_L,Cg_B,Cg_R,Cg_T,Eg_avg,Eg_edgV,Eg_edgH,Eg_in_L,Eg_in_B,Eg_in_R,Eg_in_T)

            !calculating all coefficients needed for the EGP solve
            KapE_bar_MGQDold = KapE_bar
            CALL GREY_COEFS(Eg_avg,Eg_edgV,Eg_edgH,Fxg_edgV_old,Fyg_edgH_old,fg_avg_xx,fg_avg_yy,fg_edgV_xx,fg_edgV_xy,fg_edgH_yy,&
              fg_edgH_xy,KapE,KapR,A,c,Delt,Theta,Pold_L,Pold_R,Pold_B,Pold_T,KapE_Bar,DC_xx,DL_xx,DR_xx,DC_yy,DB_yy,DT_yy,DL_xy,&
              DR_xy,DB_xy,DT_xy,PL,PR,PB,PT)

            !solve the nonlinear EGP with newton iterations
            CALL EGP_FV_NEWT(E_avg,E_edgV,E_edgH,Temp,KapE_Bar,Fx_edgV,Fy_edgH,GQD_Src,KapE_Bar_dT,EGP_Its,Deltas,dres,&
              Temp_MGQDold2,KapE_bar_MGQDold,Theta,Delt,Delx,Dely,A,Cb_L,Cb_B,Cb_R,Cb_T,E_in_L,E_in_B,E_in_R,E_in_T,F_in_L,&
              F_in_B,F_in_R,F_in_T,DC_xx,DL_xx,DR_xx,DC_yy,DB_yy,DT_yy,DL_xy,DR_xy,DB_xy,DT_xy,PL,PR,PB,PT,Gold_Hat,&
              Rhat_old,Kap0,cv,Comp_Unit,Chi,line_src,E_Bound_Low,T_Bound_Low,Conv_gr1,Conv_gr2,Maxit_GLOQD,MGQD_Its,&
              Use_Line_Search,Use_Safety_Search,Res_Calc,kapE_dT_flag)

            !writing the count of EGP iterations to output file
            IF (Res_Calc) THEN
              Status = nf90_put_var(outID,Del_T_ID,Deltas(:,1),(/1,MGQD_Its,RT_Its,t/),(/EGP_Its,1,1,1/))
              CALL HANDLE_ERR(Status)
              Status = nf90_put_var(outID,Del_E_avg_ID,Deltas(:,2),(/1,MGQD_Its,RT_Its,t/),(/EGP_Its,1,1,1/))
              CALL HANDLE_ERR(Status)
              Status = nf90_put_var(outID,Del_E_EdgV_ID,Deltas(:,3),(/1,MGQD_Its,RT_Its,t/),(/EGP_Its,1,1,1/))
              CALL HANDLE_ERR(Status)
              Status = nf90_put_var(outID,Del_E_EdgH_ID,Deltas(:,4),(/1,MGQD_Its,RT_Its,t/),(/EGP_Its,1,1,1/))
              CALL HANDLE_ERR(Status)
              Status = nf90_put_var(outID,Del_Fx_EdgV_ID,Deltas(:,5),(/1,MGQD_Its,RT_Its,t/),(/EGP_Its,1,1,1/))
              CALL HANDLE_ERR(Status)
              Status = nf90_put_var(outID,Del_Fy_EdgH_ID,Deltas(:,6),(/1,MGQD_Its,RT_Its,t/),(/EGP_Its,1,1,1/))
              CALL HANDLE_ERR(Status)

              Status = nf90_put_var(outID,dr_T_ID,dres(:,1),(/1,MGQD_Its,RT_Its,t/),(/EGP_Its,1,1,1/))
              CALL HANDLE_ERR(Status)
              Status = nf90_put_var(outID,dr_B_ID,dres(:,2),(/1,MGQD_Its,RT_Its,t/),(/EGP_Its,1,1,1/))
              CALL HANDLE_ERR(Status)
              Status = nf90_put_var(outID,dr_ML_ID,dres(:,3),(/1,MGQD_Its,RT_Its,t/),(/EGP_Its,1,1,1/))
              CALL HANDLE_ERR(Status)
              Status = nf90_put_var(outID,dr_MB_ID,dres(:,4),(/1,MGQD_Its,RT_Its,t/),(/EGP_Its,1,1,1/))
              CALL HANDLE_ERR(Status)
              Status = nf90_put_var(outID,dr_MR_ID,dres(:,5),(/1,MGQD_Its,RT_Its,t/),(/EGP_Its,1,1,1/))
              CALL HANDLE_ERR(Status)
              Status = nf90_put_var(outID,dr_MT_ID,dres(:,6),(/1,MGQD_Its,RT_Its,t/),(/EGP_Its,1,1,1/))
              CALL HANDLE_ERR(Status)
            END IF

            !writing the count of EGP iterations to output file
            Status = nf90_put_var(outID,GQD_ItCount_ID,(/EGP_Its/),(/MGQD_Its,RT_Its,t/),(/1,1,1/))
            CALL HANDLE_ERR(Status)

          ELSE
            !solve the MEB equation with the MGQD solution to find new Temp
            CALL MEB_SOLVE(Temp,Eg_avg,Bg,KapE,Temp_old,Delt,cV,Kap0,Comp_Unit)

          END IF

          !update material properties with new Temp
          CALL MATERIAL_UPDATE(RT_Src,MGQD_Src,KapE,KapB,KapR,Bg,Temp,Comp_Unit,Nu_g,Threads)

          !check convergence of MGQD solution (in grey form) for E, T
          Tconv = CONVERGENCE(Conv_Type,conv_lo,Temp,Temp_MGQDold,Temp_MGQDold2,MGQD_Its,&
            MGQD_Tnorm,MGQD_Trho)
          Econv = CONVERGENCE(Conv_Type,conv_lo,MGQD_E_avg,MGQD_E_avg_MGQDold,MGQD_E_avg_MGQDold2,MGQD_Its,&
            MGQD_Enorm,MGQD_Erho)
          MGQD_conv = Tconv.AND.Econv !if both T and E are converged, the MGQD iterations have successfully converged

          !writing MGQD convergence status (norms and spectral radii) to output file
          Status = nf90_put_var(outID,MGQD_Tnorm_ID,(/MGQD_Tnorm/),(/MGQD_Its,RT_Its,t/),(/1,1,1/))
          CALL HANDLE_ERR(Status)
          Status = nf90_put_var(outID,MGQD_Enorm_ID,(/MGQD_Enorm/),(/MGQD_Its,RT_Its,t/),(/1,1,1/))
          CALL HANDLE_ERR(Status)
          Status = nf90_put_var(outID,MGQD_Trho_ID,(/MGQD_Trho/),(/MGQD_Its,RT_Its,t/),(/1,1,1/))
          CALL HANDLE_ERR(Status)
          Status = nf90_put_var(outID,MGQD_Erho_ID,(/MGQD_Erho/),(/MGQD_Its,RT_Its,t/),(/1,1,1/))
          CALL HANDLE_ERR(Status)

          !preparing for next iteration by moving current solution -> last iterate solution
          MGQD_E_avg_MGQDold2 = MGQD_E_avg_MGQDold
          MGQD_E_avg_MGQDold = MGQD_E_avg
          Temp_MGQDold2 = Temp_MGQDold
          Temp_MGQDold = Temp

          !writing current iterate to terminal
          write(*,*) 'MGQD:     ',MGQD_Its, EGP_Its, MGQD_Tnorm, MGQD_Enorm

        END DO

        !writing the count of MGQD iterations to output file
        Status = nf90_put_var(outID,MGQD_ItCount_ID,(/MGQD_Its/),(/RT_Its,t/),(/1,1/))
        CALL HANDLE_ERR(Status)

      END IF

      Tconv = CONVERGENCE(Conv_Type,conv_ho,Temp,Temp_RTold,Temp_RTold2,RT_Its,TR_Tnorm,TR_Trho)
      Econv = CONVERGENCE(Conv_Type,conv_ho,HO_E_avg,HO_E_avg_RTold,HO_E_avg_RTold2,RT_Its,TR_Enorm,TR_Erho)
      RT_conv = Tconv.AND.Econv !if both T and E are converged, the outer/RTE iterations have successfully converged

      !writing RTE convergence status (norms and spectral radii) to output file
      Status = nf90_put_var(outID,RT_Tnorm_ID,(/TR_Tnorm/),(/RT_Its,t/),(/1,1/))
      CALL HANDLE_ERR(Status)
      Status = nf90_put_var(outID,RT_Enorm_ID,(/TR_Enorm/),(/RT_Its,t/),(/1,1/))
      CALL HANDLE_ERR(Status)
      Status = nf90_put_var(outID,RT_Trho_ID,(/TR_Trho/),(/RT_Its,t/),(/1,1/))
      CALL HANDLE_ERR(Status)
      Status = nf90_put_var(outID,RT_Erho_ID,(/TR_Erho/),(/RT_Its,t/),(/1,1/))
      CALL HANDLE_ERR(Status)

      !preparing for next iteration by moving current solution -> last iterate solution
      HO_E_avg_RTold2 = HO_E_avg_RTold
      HO_E_avg_RTold = HO_E_avg
      Temp_RTold2 = Temp_RTold
      Temp_RTold = Temp

      !writing current iterate to terminal
      write(*,*) RT_Its, MGQD_Its, TR_Tnorm, TR_Enorm

    END DO

    !writing finished time step to terminal
    write(*,*) Time, RT_Its, TR_Tnorm, TR_Enorm

    !writing the count of outer/RTE iterations to output file
    Status = nf90_put_var(outID,RT_ItCount_ID,(/RT_Its/),(/t/),(/1/))
    CALL HANDLE_ERR(Status)

    CALL NF_PUT_t_VAR(outID,Temp_ID,Temp,t)
    CALL NF_PUT_t_VAR(outID,MGQD_E_avg_ID,MGQD_E_avg,t)
    CALL NF_PUT_t_VAR(outID,MGQD_E_edgV_ID,MGQD_E_edgV,t)
    CALL NF_PUT_t_VAR(outID,MGQD_E_edgH_ID,MGQD_E_edgH,t)
    CALL NF_PUT_t_VAR(outID,MGQD_Fx_edgV_ID,MGQD_Fx_edgV,t)
    CALL NF_PUT_t_VAR(outID,MGQD_Fy_edgH_ID,MGQD_Fy_edgH,t)
    CALL NF_PUT_t_VAR(outID,HO_E_avg_ID,HO_E_avg,t)
    CALL NF_PUT_t_VAR(outID,HO_E_edgV_ID,HO_E_edgV,t)
    CALL NF_PUT_t_VAR(outID,HO_E_edgH_ID,HO_E_edgH,t)
    CALL NF_PUT_t_VAR(outID,HO_Fx_edgV_ID,HO_Fx_edgV,t)
    CALL NF_PUT_t_VAR(outID,HO_Fy_edgH_ID,HO_Fy_edgH,t)

    CALL NF_PUT_t_VAR(outID,Eg_avg_ID,Eg_avg,t)
    CALL NF_PUT_t_VAR(outID,Eg_edgV_ID,Eg_edgV,t)
    CALL NF_PUT_t_VAR(outID,Eg_edgH_ID,Eg_edgH,t)
    CALL NF_PUT_t_VAR(outID,Fxg_edgV_ID,Fxg_edgV,t)
    CALL NF_PUT_t_VAR(outID,Fyg_edgH_ID,Fyg_edgH,t)
    CALL NF_PUT_t_VAR(outID,HO_Eg_avg_ID,HO_Eg_avg,t)
    CALL NF_PUT_t_VAR(outID,HO_Eg_edgV_ID,HO_Eg_edgV,t)
    CALL NF_PUT_t_VAR(outID,HO_Eg_edgH_ID,HO_Eg_edgH,t)
    CALL NF_PUT_t_VAR(outID,HO_Fxg_edgV_ID,HO_Fxg_edgV,t)
    CALL NF_PUT_t_VAR(outID,HO_Fyg_edgH_ID,HO_Fyg_edgH,t)

    CALL NF_PUT_t_VAR(outID,I_avg_ID,I_avg,t)
    CALL NF_PUT_t_VAR(outID,I_edgV_ID,I_edgV,t)
    CALL NF_PUT_t_VAR(outID,I_edgH_ID,I_edgH,t)

    CALL NF_PUT_t_VAR(outID,KapE_Bar_ID,KapE_Bar,t)
    CALL NF_PUT_t_VAR(outID,KapB_ID,KapB,t)
    CALL NF_PUT_t_VAR(outID,KapE_ID,KapE,t)
    CALL NF_PUT_t_VAR(outID,KapR_ID,KapR,t)
    CALL NF_PUT_t_VAR(outID,Bg_ID,Bg,t)

    CALL NF_PUT_t_VAR(outID,Cg_L_ID,Cg_L,t)
    CALL NF_PUT_t_VAR(outID,Cg_B_ID,Cg_B,t)
    CALL NF_PUT_t_VAR(outID,Cg_R_ID,Cg_R,t)
    CALL NF_PUT_t_VAR(outID,Cg_T_ID,Cg_T,t)
    CALL NF_PUT_t_VAR(outID,Eg_in_L_ID,Eg_in_L,t)
    CALL NF_PUT_t_VAR(outID,Eg_in_B_ID,Eg_in_B,t)
    CALL NF_PUT_t_VAR(outID,Eg_in_R_ID,Eg_in_R,t)
    CALL NF_PUT_t_VAR(outID,Eg_in_T_ID,Eg_in_T,t)
    CALL NF_PUT_t_VAR(outID,Fg_in_L_ID,Fg_in_L,t)
    CALL NF_PUT_t_VAR(outID,Fg_in_B_ID,Fg_in_B,t)
    CALL NF_PUT_t_VAR(outID,Fg_in_R_ID,Fg_in_R,t)
    CALL NF_PUT_t_VAR(outID,Fg_in_T_ID,Fg_in_T,t)

    CALL NF_PUT_t_VAR(outID,fg_avg_xx_ID,fg_avg_xx,t)
    CALL NF_PUT_t_VAR(outID,fg_avg_yy_ID,fg_avg_yy,t)
    CALL NF_PUT_t_VAR(outID,fg_edgV_xx_ID,fg_edgV_xx,t)
    CALL NF_PUT_t_VAR(outID,fg_edgV_xy_ID,fg_edgV_xy,t)
    CALL NF_PUT_t_VAR(outID,fg_edgH_yy_ID,fg_edgH_yy,t)
    CALL NF_PUT_t_VAR(outID,fg_edgH_xy_ID,fg_edgH_xy,t)

    CALL NF_PUT_t_VAR(outID,G_old_ID,G_old,t)
    CALL NF_PUT_t_VAR(outID,Pold_L_ID,Pold_L,t)
    CALL NF_PUT_t_VAR(outID,Pold_B_ID,Pold_B,t)
    CALL NF_PUT_t_VAR(outID,Pold_R_ID,Pold_R,t)
    CALL NF_PUT_t_VAR(outID,Pold_T_ID,Pold_T,t)

    IF (use_grey .EQ. 1) THEN
      CALL NF_PUT_t_VAR(outID,E_avg_ID,E_avg,t)
      CALL NF_PUT_t_VAR(outID,E_edgV_ID,E_edgV,t)
      CALL NF_PUT_t_VAR(outID,E_edgH_ID,E_edgH,t)
      CALL NF_PUT_t_VAR(outID,Fx_edgV_ID,Fx_edgV,t)
      CALL NF_PUT_t_VAR(outID,Fy_edgH_ID,Fy_edgH,t)

      CALL NF_PUT_t_VAR(outID,Cb_L_ID,Cb_L,t)
      CALL NF_PUT_t_VAR(outID,Cb_B_ID,Cb_B,t)
      CALL NF_PUT_t_VAR(outID,Cb_R_ID,Cb_R,t)
      CALL NF_PUT_t_VAR(outID,Cb_T_ID,Cb_T,t)
      CALL NF_PUT_t_VAR(outID,E_in_L_ID,E_in_L,t)
      CALL NF_PUT_t_VAR(outID,E_in_B_ID,E_in_B,t)
      CALL NF_PUT_t_VAR(outID,E_in_R_ID,E_in_R,t)
      CALL NF_PUT_t_VAR(outID,E_in_T_ID,E_in_T,t)
      CALL NF_PUT_t_VAR(outID,F_in_L_ID,F_in_L,t)
      CALL NF_PUT_t_VAR(outID,F_in_B_ID,F_in_B,t)
      CALL NF_PUT_t_VAR(outID,F_in_R_ID,F_in_R,t)
      CALL NF_PUT_t_VAR(outID,F_in_T_ID,F_in_T,t)

      CALL NF_PUT_t_VAR(outID,DC_xx_ID,DC_xx,t)
      CALL NF_PUT_t_VAR(outID,DL_xx_ID,DL_xx,t)
      CALL NF_PUT_t_VAR(outID,DR_xx_ID,DR_xx,t)
      CALL NF_PUT_t_VAR(outID,DC_yy_ID,DC_yy,t)
      CALL NF_PUT_t_VAR(outID,DB_yy_ID,DB_yy,t)
      CALL NF_PUT_t_VAR(outID,DT_yy_ID,DT_yy,t)
      CALL NF_PUT_t_VAR(outID,DL_xy_ID,DL_xy,t)
      CALL NF_PUT_t_VAR(outID,DB_xy_ID,DB_xy,t)
      CALL NF_PUT_t_VAR(outID,DR_xy_ID,DR_xy,t)
      CALL NF_PUT_t_VAR(outID,DT_xy_ID,DT_xy,t)

      CALL NF_PUT_t_VAR(outID,Gold_hat_ID,Gold_hat,t)
      CALL NF_PUT_t_VAR(outID,Rhat_old_ID,Rhat_old,t)
      CALL NF_PUT_t_VAR(outID,PL_ID,PL,t)
      CALL NF_PUT_t_VAR(outID,PB_ID,PB,t)
      CALL NF_PUT_t_VAR(outID,PR_ID,PR,t)
      CALL NF_PUT_t_VAR(outID,PT_ID,PT,t)
    END IF

    IF (Time .GE. tlen) EXIT
  END DO

END SUBROUTINE TRT_MLQD_ALGORITHM

!==================================================================================================================================!
!
!==================================================================================================================================!

END MODULE ALGORITHMS
