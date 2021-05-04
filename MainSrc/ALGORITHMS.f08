MODULE ALGORITHMS

  USE TRANSPORT_SOLVES
  USE UPDATES
  USE CONVERGENCE_CHECKS
  USE MLOQD_SOLVES
  USE GLOQD_SOLVES
  USE OUTPUTS
  USE INPUTS
  USE INITIALIZERS
  USE netcdf
  USE NCDF_IO
  USE POD_ROUTINES
  USE DMD_ROUTINES
  USE GRID_FUNCTIONS
  USE WAVEPROP_TOOLS

  IMPLICIT NONE

CONTAINS

!==================================================================================================================================!
!
!==================================================================================================================================!
SUBROUTINE TRT_MLQD_ALGORITHM(Omega_x,Omega_y,quad_weight,Nu_g,Delx,Dely,Delt,tlen,Theta,Start_Time,c,cV,h,pi,Kap0,&
  erg,Comp_Unit,Conv_ho,Conv_lo,Conv_gr1,Conv_gr2,bcT_left,bcT_bottom,bcT_right,bcT_top,Tini,chi,line_src,E_Bound_Low,&
  T_Bound_Low,use_grey,Conv_Type,Maxit_RTE,Threads,BC_Type,Maxit_MLOQD,Maxit_GLOQD,N_x,N_y,N_m,N_g,&
  N_t,Res_Calc,Use_Line_Search,Use_Safety_Search,run_type,kapE_dT_flag,outID,N_x_ID,N_y_ID,N_m_ID,N_g_ID,N_t_ID,&
  N_edgV_ID,N_edgH_ID,N_xc_ID,N_yc_ID,Quads_ID,RT_Its_ID,MGQD_Its_ID,GQD_Its_ID,Norm_Types_ID,MGQD_ResTypes_ID,&
  Boundaries_ID,out_freq,I_out,HO_Eg_out,HO_Fg_out,HO_E_out,HO_F_out,Eg_out,Fg_out,MGQD_E_out,MGQD_F_out,QDfg_out,&
  E_out,F_out,D_out,old_parms_out,its_out,conv_out,kap_out,Src_out,POD_dsets,POD_err,POD_Type,xlen,ylen,&
  Direc_Diff,Mat,Kappa_Mult,restart_outfile,restart_freq,restart_infile,N_dsets,DMD_dsets,DMD_Type,dset_times,&
  N_dsets_ID,Init_out,qdf_infile, qdfin_fxy_flag, qdfin_bctype)

  !---------------Solution Parameters----------------!
  REAL*8,INTENT(IN):: Omega_x(:), Omega_y(:), quad_weight(:), Nu_g(:), Kappa_Mult(:)
  REAL*8,INTENT(IN):: Delx(:), Dely(:), Delt, Theta, tlen, xlen, ylen
  REAL*8,INTENT(IN):: Start_Time, dset_times(:)
  REAL*8,INTENT(IN):: c, cV, h, pi, Kap0, erg
  REAL*8,INTENT(IN):: Comp_Unit, Conv_ho, Conv_lo, Conv_gr1, Conv_gr2
  REAL*8,INTENT(IN):: bcT_left, bcT_bottom, bcT_right, bcT_top, Tini
  REAL*8,INTENT(IN):: chi, line_src, E_Bound_Low, T_Bound_Low, POD_err(:)
  INTEGER,INTENT(IN):: N_x, N_y, N_m, N_g, N_t, Direc_Diff, Mat(:,:), N_dsets, qdfin_fxy_flag
  INTEGER,INTENT(IN):: use_grey, Conv_Type, Threads, BC_Type(:), Maxit_RTE, Maxit_MLOQD, Maxit_GLOQD
  LOGICAL,INTENT(IN):: Res_Calc, Use_Line_Search, Use_Safety_Search, kapE_dT_flag
  CHARACTER(*),INTENT(IN):: run_type, POD_Type, POD_dsets(:), restart_outfile, restart_infile
  CHARACTER(*),INTENT(IN):: DMD_Type, DMD_dsets(:), qdf_infile, qdfin_bctype

  !----------------Output File ID's------------------!
  INTEGER,INTENT(IN):: outID
  INTEGER,INTENT(IN):: N_x_ID, N_y_ID, N_m_ID, N_g_ID, N_t_ID, N_edgV_ID, N_edgH_ID, N_xc_ID, N_yc_ID, Quads_ID
  INTEGER,INTENT(IN):: RT_Its_ID, MGQD_Its_ID, GQD_Its_ID, Norm_Types_ID, MGQD_ResTypes_ID, Boundaries_ID, N_dsets_ID
  INTEGER,INTENT(IN):: out_freq, I_out, HO_Eg_out, HO_Fg_out, HO_E_out, HO_F_out, restart_freq
  INTEGER,INTENT(IN):: Eg_out, Fg_out, MGQD_E_out, MGQD_F_out, QDfg_out
  INTEGER,INTENT(IN):: E_out, F_out, D_out
  INTEGER,INTENT(IN):: old_parms_out, its_out, conv_out, kap_out, Src_out, Init_out

  !---------------Material Properties----------------!
  REAL*8,ALLOCATABLE:: Bg(:,:,:), KapE(:,:,:), KapB(:,:,:), KapR(:,:,:), A(:,:)
  REAL*8,ALLOCATABLE:: KapE_old(:,:,:), KapR_old(:,:,:)
  REAL*8,ALLOCATABLE:: MGQD_Src(:,:,:), MGQD_Src_old(:,:,:), RT_Src(:,:,:,:)

  !--------------Radiation Intensities---------------!
  REAL*8,ALLOCATABLE:: I_crn_old(:,:,:,:)
  REAL*8,ALLOCATABLE:: I_avg(:,:,:,:), I_edgV(:,:,:,:), I_edgH(:,:,:,:)
  REAL*8,ALLOCATABLE:: I_crn(:,:,:,:), Ic_edgV(:,:,:,:), Ic_edgH(:,:,:,:)

  !--------------High-Order Quantities---------------!
  REAL*8,ALLOCATABLE:: Hg_avg_xx(:,:,:), Hg_avg_yy(:,:,:), Hg_avg_xy(:,:,:)
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
  REAL*8,ALLOCATABLE:: fg_avg_xx(:,:,:), fg_avg_yy(:,:,:), fg_avg_xy(:,:,:)
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
  REAL*8:: Time, Fdt_Weight, Fin_Weight, Ein_Weight
  INTEGER:: MGQD_Its, EGP_Its, Status
  INTEGER:: RT_Its, RT_start_Its, t, HO_Form
  INTEGER,ALLOCATABLE:: MGQD_Kits(:), GQD_Kits(:)
  LOGICAL:: RT_Conv, MGQD_conv, Tconv, Econv

  !-------------------Variable ID's------------------!
  INTEGER:: Temp_ID, E_avg_ID, E_edgV_ID, E_edgH_ID, MGQD_E_avg_ID, MGQD_E_edgV_ID, MGQD_E_edgH_ID, HO_E_avg_ID
  INTEGER:: HO_E_edgV_ID, HO_E_edgH_ID, Fx_edgV_ID, Fy_edgH_ID, MGQD_Fx_edgV_ID, MGQD_Fy_edgH_ID, HO_Fx_edgV_ID
  INTEGER:: HO_Fy_edgH_ID, Eg_avg_ID, Eg_edgV_ID, Eg_edgH_ID, HO_Eg_avg_ID, HO_Eg_edgV_ID, HO_Eg_edgH_ID
  INTEGER:: Fxg_edgV_ID, Fyg_edgH_ID, HO_Fxg_edgV_ID, HO_Fyg_edgH_ID, I_avg_ID, I_edgV_ID, I_edgH_ID
  INTEGER:: KapE_Bar_ID, KapB_ID, KapE_ID, KapR_ID, Bg_ID
  INTEGER:: RT_Residual_ID, MGQD_Residual_ID, MGQD_BC_Residual_ID
  INTEGER:: Del_T_ID, Del_E_avg_ID, Del_E_edgV_ID, Del_E_edgH_ID, Del_Fx_edgV_ID, Del_Fy_edgH_ID
  INTEGER:: RT_ItCount_ID, MGQD_ItCount_ID, GQD_ItCount_ID, MGQD_KIts_ID, GQD_KIts_ID
  INTEGER:: RT_Tnorm_ID, RT_Enorm_ID, MGQD_Tnorm_ID, MGQD_Enorm_ID
  INTEGER:: RT_Trho_ID, RT_Erho_ID, MGQD_Trho_ID, MGQD_Erho_ID
  INTEGER:: Cg_L_ID, Cg_B_ID, Cg_R_ID, Cg_T_ID, Eg_in_L_ID, Eg_in_B_ID, Eg_in_R_ID, Eg_in_T_ID
  INTEGER:: Fg_in_L_ID, Fg_in_B_ID, Fg_in_R_ID, Fg_in_T_ID
  INTEGER:: Cb_L_ID, Cb_B_ID, Cb_R_ID, Cb_T_ID, E_in_L_ID, E_in_B_ID, E_in_R_ID, E_in_T_ID
  INTEGER:: F_in_L_ID, F_in_B_ID, F_in_R_ID, F_in_T_ID
  INTEGER:: fg_avg_xx_ID, fg_avg_yy_ID, fg_avg_xy_ID, fg_edgV_xx_ID, fg_edgV_xy_ID, fg_edgH_yy_ID, fg_edgH_xy_ID
  INTEGER:: DC_xx_ID, DL_xx_ID, DR_xx_ID, DC_yy_ID, DB_yy_ID, DT_yy_ID, DL_xy_ID, DB_xy_ID, DR_xy_ID, DT_xy_ID
  INTEGER:: G_old_ID, Pold_L_ID, Pold_B_ID, Pold_R_ID, Pold_T_ID
  INTEGER:: Gold_hat_ID, Rhat_old_ID, PL_ID, PB_ID, PR_ID, PT_ID
  INTEGER:: dr_T_ID, dr_B_ID, dr_ML_ID, dr_MB_ID, dr_MR_ID, dr_MT_ID
  INTEGER:: rrank_BCg_ID, rrank_fg_avg_xx_ID, rrank_fg_edgV_xx_ID, rrank_fg_avg_yy_ID, rrank_fg_edgH_yy_ID
  INTEGER:: rrank_fg_edgV_xy_ID, rrank_fg_edgH_xy_ID
  INTEGER:: qdf_infileID
  INTEGER:: Temp_XWvSpeed_ID, E_XWvSpeed_ID

  !--------------------------------------------------!
  REAL*8,ALLOCATABLE:: C_BCg(:), S_BCg(:), U_BCg(:), V_BCg(:)
  COMPLEX*16,ALLOCATABLE:: L_BCg(:), W_BCg(:), B_BCg(:)
  INTEGER,ALLOCATABLE:: rrank_BCg(:)

  REAL*8,ALLOCATABLE:: C_fg_avg_xx(:), S_fg_avg_xx(:), U_fg_avg_xx(:), V_fg_avg_xx(:)
  REAL*8,ALLOCATABLE:: C_fg_edgV_xx(:), S_fg_edgV_xx(:), U_fg_edgV_xx(:), V_fg_edgV_xx(:)
  REAL*8,ALLOCATABLE:: C_fg_avg_yy(:), S_fg_avg_yy(:), U_fg_avg_yy(:), V_fg_avg_yy(:)
  REAL*8,ALLOCATABLE:: C_fg_edgH_yy(:), S_fg_edgH_yy(:), U_fg_edgH_yy(:), V_fg_edgH_yy(:)
  REAL*8,ALLOCATABLE:: C_fg_edgV_xy(:), S_fg_edgV_xy(:), U_fg_edgV_xy(:), V_fg_edgV_xy(:)
  REAL*8,ALLOCATABLE:: C_fg_edgH_xy(:), S_fg_edgH_xy(:), U_fg_edgH_xy(:), V_fg_edgH_xy(:)
  !
  COMPLEX*16,ALLOCATABLE:: L_fg_avg_xx(:),  W_fg_avg_xx(:),  B_fg_avg_xx(:)
  COMPLEX*16,ALLOCATABLE:: L_fg_edgV_xx(:), W_fg_edgV_xx(:), B_fg_edgV_xx(:)
  COMPLEX*16,ALLOCATABLE:: L_fg_avg_yy(:),  W_fg_avg_yy(:),  B_fg_avg_yy(:)
  COMPLEX*16,ALLOCATABLE:: L_fg_edgH_yy(:), W_fg_edgH_yy(:), B_fg_edgH_yy(:)
  COMPLEX*16,ALLOCATABLE:: L_fg_edgV_xy(:), W_fg_edgV_xy(:), B_fg_edgV_xy(:)
  COMPLEX*16,ALLOCATABLE:: L_fg_edgH_xy(:), W_fg_edgH_xy(:), B_fg_edgH_xy(:)
  !
  INTEGER,ALLOCATABLE:: rrank_fg_avg_xx(:), rrank_fg_edgV_xx(:), rrank_fg_avg_yy(:), rrank_fg_edgH_yy(:)
  INTEGER,ALLOCATABLE:: rrank_fg_edgV_xy(:), rrank_fg_edgH_xy(:)

  REAL*8,ALLOCATABLE:: C_I_avg(:), S_I_avg(:), U_I_avg(:), V_I_avg(:)
  REAL*8,ALLOCATABLE:: C_I_edgV(:), S_I_edgV(:), U_I_edgV(:), V_I_edgV(:)
  REAL*8,ALLOCATABLE:: C_I_edgH(:), S_I_edgH(:), U_I_edgH(:), V_I_edgH(:)
  INTEGER,ALLOCATABLE:: rrank_I_avg(:), rrank_I_edgV(:), rrank_I_edgH(:)

  REAL*8:: tlen_d, xlen_d, ylen_d, Tini_d, Delt_d, Start_Time_d, DMD_Time
  REAL*8,ALLOCATABLE:: Delx_d(:), Dely_d(:), bcT_d(:)
  INTEGER:: dN_x, dN_y, dN_m, dN_g, dN_t, fg_pod_out, fg_dmd_out, DMDgsum, PODgsum
  INTEGER:: Current_Database, DMDgsum_ID, PODgsum_ID, PODerr_ID, gsum
  INTEGER,ALLOCATABLE:: BC_Type_d(:)

  REAL*8,ALLOCATABLE:: Sim_Grid_Avg(:), Sim_Grid_EdgV(:), Sim_Grid_EdgH(:), Sim_Grid_Bnds(:)
  REAL*8,ALLOCATABLE:: Dat_Grid_Avg(:), Dat_Grid_EdgV(:), Dat_Grid_EdgH(:), Dat_Grid_Bnds(:)
  REAL*8,ALLOCATABLE:: Sim_TGrid(:), Dat_TGrid(:)
  INTEGER,ALLOCATABLE:: GMap_xyAvg(:), GMap_xyEdgV(:), GMap_xyEdgH(:), GMap_xyBnds(:)
  INTEGER,ALLOCATABLE:: VMap_xyAvg(:), VMap_xyEdgV(:), VMap_xyEdgH(:), VMap_xyBnds(:)
  INTEGER,ALLOCATABLE:: TMap(:)

  REAL*8:: E_xWvSpd, T_xWvSpd

  INTEGER:: resf_unit = 307

  !===========================================================================!
  !                                                                           !
  !     INITIALIZING ARRAYS                                                   !
  !                                                                           !
  !===========================================================================!
  CALL MISC_INIT(Delx,Dely,A)
  CALL RT_INIT(I_avg,I_edgV,I_edgH,I_crn,I_crn_old,Ic_edgV,Ic_edgH,Hg_avg_xx,Hg_avg_yy,Hg_avg_xy,Hg_edgV_xx,&
    Hg_edgV_xy,Hg_edgH_yy,Hg_edgH_xy,HO_Eg_avg,HO_Eg_edgV,HO_Eg_edgH,HO_Fxg_edgV,HO_Fyg_edgH,HO_E_avg,&
    HO_E_edgV,HO_E_edgH,HO_Fx_edgV,HO_Fy_edgH,N_y,N_x,N_m,N_g,Tini,comp_unit,nu_g,bcT_left,bcT_right,&
    bcT_top,bcT_bottom,BC_Type,pi,c)
  CALL MGQD_INIT(Eg_avg,Eg_edgV,Eg_edgH,Fxg_edgV,Fyg_edgH,Eg_avg_old,Eg_edgV_old,Eg_edgH_old,Fxg_edgV_old,Fyg_edgH_old,&
    fg_avg_xx,fg_avg_yy,fg_avg_xy,fg_edgV_xx,fg_edgV_xy,fg_edgH_yy,fg_edgH_xy,fg_avg_xx_old,fg_avg_yy_old,fg_edgV_xx_old,&
    fg_edgV_xy_old,fg_edgH_yy_old,fg_edgH_xy_old,Cg_L,Cg_B,Cg_R,Cg_T,Eg_in_L,Eg_in_B,Eg_in_R,Eg_in_T,Fg_in_L,Fg_in_B,Fg_in_R,&
    Fg_in_T,I_edgV,I_edgH,Omega_x,Omega_y,quad_weight,c,Comp_Unit,N_y,N_x,N_g,MGQD_E_avg,MGQD_E_edgV,MGQD_E_edgH,MGQD_Fx_edgV,&
    MGQD_Fy_edgH,G_old,Pold_L,Pold_B,Pold_R,Pold_T,BC_Type,Tini,nu_g,pi,Threads,Maxit_GLOQD,MGQD_Kits,GQD_Kits)
  CALL GQD_INIT(E_avg,E_edgV,E_edgH,Fx_edgV,Fy_edgH,E_avg_old,KapE_bar_old,GQD_Src_old,Gold_Hat,Rhat_old,&
    KapE_Bar,KapE_Bar_dT,GQD_Src,DC_xx,DL_xx,DR_xx,DC_yy,DB_yy,DT_yy,DL_xy,DB_xy,DR_xy,DT_xy,PL,PB,PR,PT,Cb_L,&
    Cb_B,Cb_R,Cb_T,E_in_L,E_in_B,E_in_R,E_in_T,F_in_L,F_in_B,F_in_R,F_in_T,KapE_bar_MGQDold,Eg_in_L,Eg_in_B,&
    Eg_in_R,Eg_in_T,Fg_in_L,Fg_in_B,Fg_in_R,Fg_in_T,Tini,nu_g,c,Comp_Unit,pi,N_y,N_x,N_g,BC_Type)
  CALL TEMP_INIT(Temp,RT_Src,MGQD_Src,MGQD_Src_old,KapE,KapB,KapR,KapE_old,KapR_old,Bg,N_y,N_x,N_m,N_g,Tini,&
    Comp_Unit,Nu_g,Temp_Old,Threads,Mat,Kappa_Mult)

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

  IF (restart_infile .NE. '') THEN
    CALL RESTART_IN(Time, Temp, I_crn, fg_avg_xx, fg_avg_yy, fg_edgV_xx, fg_edgV_xy, fg_edgH_yy, fg_edgH_xy, RT_Src, MGQD_Src,&
      GQD_Src, E_avg, Eg_avg, Eg_edgV, Eg_edgH, Fxg_edgV, Fyg_edgH, KapE_Bar, KapE, KapR, Cg_L, Cg_B, Cg_R, Cg_T, restart_infile)
  END IF

  !===========================================================================!
  !                                                                           !
  !     INITIALIZING OUTPUT FILE                                              !
  !                                                                           !
  !===========================================================================!
  fg_pod_out = 0
  fg_dmd_out = 0
  IF ((run_type .EQ. 'mg_pod').AND.(POD_Type .EQ. 'fg')) THEN
    fg_pod_out = 1
  ELSE IF ((run_type .EQ. 'mg_dmd').AND.(DMD_Type .EQ. 'fg')) THEN
    fg_dmd_out = 1
  END IF
  CALL OUTFILE_VARDEFS(outID,Res_Calc,out_freq,I_out,HO_Eg_out,HO_Fg_out,HO_E_out,HO_F_out,Eg_out,Fg_out,MGQD_E_out,&
    MGQD_F_out,QDfg_out,fg_pod_out,E_out,F_out,D_out,old_parms_out,its_out,conv_out,kap_out,Src_out,N_x_ID,N_y_ID,N_m_ID,N_g_ID,&
    N_t_ID,N_edgV_ID,N_edgH_ID,N_xc_ID,N_yc_ID,Quads_ID,RT_Its_ID,MGQD_Its_ID,GQD_Its_ID,Norm_Types_ID,MGQD_ResTypes_ID,&
    Boundaries_ID,Temp_ID,E_avg_ID,E_edgV_ID,E_edgH_ID,MGQD_E_avg_ID,MGQD_E_edgV_ID,MGQD_E_edgH_ID,HO_E_avg_ID,HO_E_edgV_ID,&
    HO_E_edgH_ID,Fx_edgV_ID,Fy_edgH_ID,MGQD_Fx_edgV_ID,MGQD_Fy_edgH_ID,HO_Fx_edgV_ID,HO_Fy_edgH_ID,Eg_avg_ID,Eg_edgV_ID,&
    Eg_edgH_ID,HO_Eg_avg_ID,HO_Eg_edgV_ID,HO_Eg_edgH_ID,Fxg_edgV_ID,Fyg_edgH_ID,HO_Fxg_edgV_ID,HO_Fyg_edgH_ID,I_avg_ID,&
    I_edgV_ID,I_edgH_ID,KapE_Bar_ID,KapB_ID,KapE_ID,KapR_ID,Bg_ID,RT_Residual_ID,MGQD_Residual_ID,MGQD_BC_Residual_ID,&
    Del_T_ID,Del_E_avg_ID,Del_E_edgV_ID,Del_E_edgH_ID,Del_Fx_edgV_ID,Del_Fy_edgH_ID,RT_ItCount_ID,MGQD_ItCount_ID,&
    GQD_ItCount_ID,MGQD_KIts_ID,GQD_KIts_ID,RT_Tnorm_ID,RT_Enorm_ID,MGQD_Tnorm_ID,MGQD_Enorm_ID,RT_Trho_ID,RT_Erho_ID,&
    MGQD_Trho_ID,MGQD_Erho_ID,Cg_L_ID,Cg_B_ID,Cg_R_ID,Cg_T_ID,Eg_in_L_ID,Eg_in_B_ID,Eg_in_R_ID,Eg_in_T_ID,Fg_in_L_ID,&
    Fg_in_B_ID,Fg_in_R_ID,Fg_in_T_ID,Cb_L_ID,Cb_B_ID,Cb_R_ID,Cb_T_ID,E_in_L_ID,E_in_B_ID,E_in_R_ID,E_in_T_ID,F_in_L_ID,&
    F_in_B_ID,F_in_R_ID,F_in_T_ID,fg_avg_xx_ID,fg_avg_yy_ID,fg_avg_xy_ID,fg_edgV_xx_ID,fg_edgV_xy_ID,fg_edgH_yy_ID,&
    fg_edgH_xy_ID,DC_xx_ID,DL_xx_ID,DR_xx_ID,DC_yy_ID,DB_yy_ID,DT_yy_ID,DL_xy_ID,DB_xy_ID,DR_xy_ID,DT_xy_ID,G_old_ID,&
    Pold_L_ID,Pold_B_ID,Pold_R_ID,Pold_T_ID,Gold_hat_ID,Rhat_old_ID,PL_ID,PB_ID,PR_ID,PT_ID,dr_T_ID,dr_B_ID,dr_ML_ID,&
    dr_MB_ID,dr_MR_ID,dr_MT_ID,rrank_BCg_ID,rrank_fg_avg_xx_ID,rrank_fg_edgV_xx_ID,rrank_fg_avg_yy_ID,rrank_fg_edgH_yy_ID,&
    rrank_fg_edgV_xy_ID,rrank_fg_edgH_xy_ID,fg_dmd_out,N_dsets_ID,DMDgsum_ID,PODgsum_ID,PODerr_ID,Temp_XWvSpeed_ID,E_XWvSpeed_ID)

  IF (Init_out .EQ. 1) THEN
    CALL TIMESTEP_OUTS(outID,1,out_freq,I_out,HO_Eg_out,HO_Fg_out,HO_E_out,HO_F_out,Eg_out,Fg_out,MGQD_E_out,MGQD_F_out,&
      QDfg_out,E_out,F_out,D_out,old_parms_out,its_out,conv_out,kap_out,Src_out,use_grey,Temp_ID,E_avg_ID,E_edgV_ID,E_edgH_ID,&
      MGQD_E_avg_ID,MGQD_E_edgV_ID,MGQD_E_edgH_ID,HO_E_avg_ID,HO_E_edgV_ID,HO_E_edgH_ID,Fx_edgV_ID,Fy_edgH_ID,&
      MGQD_Fx_edgV_ID,MGQD_Fy_edgH_ID,HO_Fx_edgV_ID,HO_Fy_edgH_ID,Eg_avg_ID,Eg_edgV_ID,Eg_edgH_ID,HO_Eg_avg_ID,&
      HO_Eg_edgV_ID,HO_Eg_edgH_ID,Fxg_edgV_ID,Fyg_edgH_ID,HO_Fxg_edgV_ID,HO_Fyg_edgH_ID,I_avg_ID,I_edgV_ID,I_edgH_ID,&
      KapE_Bar_ID,KapB_ID,KapE_ID,KapR_ID,Bg_ID,Cg_L_ID,Cg_B_ID,Cg_R_ID,Cg_T_ID,Eg_in_L_ID,Eg_in_B_ID,Eg_in_R_ID,&
      Eg_in_T_ID,Fg_in_L_ID,Fg_in_B_ID,Fg_in_R_ID,Fg_in_T_ID,Cb_L_ID,Cb_B_ID,Cb_R_ID,Cb_T_ID,E_in_L_ID,E_in_B_ID,&
      E_in_R_ID,E_in_T_ID,F_in_L_ID,F_in_B_ID,F_in_R_ID,F_in_T_ID,fg_avg_xx_ID,fg_avg_yy_ID,fg_avg_xy_ID,fg_edgV_xx_ID,&
      fg_edgV_xy_ID,fg_edgH_yy_ID,fg_edgH_xy_ID,DC_xx_ID,DL_xx_ID,DR_xx_ID,DC_yy_ID,DB_yy_ID,DT_yy_ID,DL_xy_ID,&
      DB_xy_ID,DR_xy_ID,DT_xy_ID,G_old_ID,Pold_L_ID,Pold_B_ID,Pold_R_ID,Pold_T_ID,Gold_hat_ID,Rhat_old_ID,PL_ID,&
      PB_ID,PR_ID,PT_ID, Temp,E_avg,E_edgV,E_edgH,MGQD_E_avg,MGQD_E_edgV,MGQD_E_edgH,HO_E_avg,HO_E_edgV,HO_E_edgH,&
      Fx_edgV,Fy_edgH,MGQD_Fx_edgV,MGQD_Fy_edgH,HO_Fx_edgV,HO_Fy_edgH,Eg_avg,Eg_edgV,Eg_edgH,HO_Eg_avg,HO_Eg_edgV,&
      HO_Eg_edgH,Fxg_edgV,Fyg_edgH,HO_Fxg_edgV,HO_Fyg_edgH,I_avg,I_edgV,I_edgH,KapE_Bar,KapB,KapE,KapR,Bg,Cg_L,Cg_B,&
      Cg_R,Cg_T,Eg_in_L,Eg_in_B,Eg_in_R,Eg_in_T,Fg_in_L,Fg_in_B,Fg_in_R,Fg_in_T,Cb_L,Cb_B,Cb_R,Cb_T,E_in_L,E_in_B,&
      E_in_R,E_in_T,F_in_L,F_in_B,F_in_R,F_in_T,fg_avg_xx,fg_avg_yy,fg_avg_xy,fg_edgV_xx,fg_edgV_xy,fg_edgH_yy,fg_edgH_xy,DC_xx,&
      DL_xx,DR_xx,DC_yy,DB_yy,DT_yy,DL_xy,DB_xy,DR_xy,DT_xy,G_old,Pold_L,Pold_B,Pold_R,Pold_T,Gold_hat,Rhat_old,PL,PB,PR,PT)
  END IF

  !===========================================================================!
  !                                                                           !
  !     CONFIGURING ALGORITHM BASED ON 'run_type'                             !
  !                                                                           !
  !===========================================================================!
  IF ( run_type .EQ. 'tr_no_qd' ) THEN
    RT_start_Its = 1
    HO_Form = 0 !Use Boltzmann solver
    ! Use QD (p1) equations for low-order
    Fdt_Weight = 1d0
    ! Quasidiffusion BC
    Fin_Weight = 1d0
    Ein_Weight = 1d0

  ELSE IF ( run_type .EQ. 'mlqd' ) THEN
    RT_start_Its = 2
    HO_Form = 0 !Use Boltzmann solver
    ! Use QD (p1) equations for low-order
    Fdt_Weight = 1d0
    ! Quasidiffusion BC
    Fin_Weight = 1d0
    Ein_Weight = 1d0

  ELSE IF (( run_type .EQ. 'diff' ).OR.( run_type .EQ. 'fld' )) THEN
    RT_start_Its = 1
    HO_Form = 2 !Turn off boltzmann solve
    ! Use diffusion equations for low-order (turn off dF/dt)
    Fdt_Weight = 0d0
    ! Marshak BC
    Fin_Weight = 2d0
    Ein_Weight = 0d0

  ELSE IF ( run_type .EQ. 'p1' ) THEN
    RT_start_Its = 1
    HO_Form = 2 !Turn off boltzmann solve
    ! Use QD (p1) equations for low-order
    Fdt_Weight = 1d0
    ! Marshak BC
    Fin_Weight = 2d0
    Ein_Weight = 0d0

  ELSE IF ( run_type .EQ. 'p13' ) THEN
    RT_start_Its = 1
    HO_Form = 2 !Turn off boltzmann solve
    ! Use QD (p1) equations for low-order with 1/3 weight on dF/dt
    Fdt_Weight = 1d0/3d0
    ! Marshak BC
    Fin_Weight = 2d0
    Ein_Weight = 0d0

  ELSE IF ((run_type .EQ. 'mg_pod').AND.(POD_Type .EQ. 'fg')) THEN
    RT_start_Its = 1
    HO_Form = 1 !Turn off boltzmann solve, replace with POD approximation of QD tensor
    ! Use QD (p1) equations for low-order
    Fdt_Weight = 1d0
    ! Quasidiffusion BC
    Fin_Weight = 1d0
    Ein_Weight = 1d0

  ELSE IF ( run_type .EQ. 'mg_dmd' ) THEN
    RT_start_Its = 1
    HO_Form = 3 !Turn off boltzmann solve, replace with DMD approximation of QD tensor
    ! Use QD (p1) equations for low-order
    Fdt_Weight = 1d0
    ! Quasidiffusion BC
    Fin_Weight = 1d0
    Ein_Weight = 1d0

  ELSE IF ( run_type .EQ. 'qdf_in' ) THEN
    RT_start_Its = 1
    HO_Form = 4 !Direct input of some QD tensor
    ! Use QD (p1) equations for low-order
    Fdt_Weight = 1d0
    !Choosing boundary condition
    IF (qdfin_bctype .EQ. 'qd') THEN! Quasidiffusion BC
      Fin_Weight = 1d0
      Ein_Weight = 1d0
    ELSE IF (qdfin_bctype .EQ. 'marshak') THEN! Marshak BC
      Fin_Weight = 2d0
      Ein_Weight = 0d0
    END IF

  ELSE
    WRITE(*,*)
    WRITE(*,*)' ERROR IN ALGORITHMS: Unknown run_type/POD_Type detected'
    STOP
    RT_start_Its = Maxit_RTE + 1
  END IF

  !===========================================================================!
  !                                                                           !
  !     PROBLEM SOLVE (BEGIN TIME STEP LOOP)                                  !
  !                                                                           !
  !===========================================================================!
  IF (Init_out .EQ. 0) THEN
    t = 0
  ELSE
    t = 1
  END IF
  Time = Start_Time
  ! t = 0
  Current_Database = 0
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
    !     CHECKING DMD/POD DATABASE TIME RANGES                                 !
    !     -> LOADING NEW DATABASE WHEN ENTERING SPECIFIED TEMPORAL RANGE        !
    !                                                                           !
    !===========================================================================!
    IF ( ANY( HO_Form .EQ. (/1,3/) ) ) THEN
      IF (time .GE. dset_times(Current_Database + 1)) THEN
        Current_Database = Current_Database + 1

        !--------------------------------------------------!
        !                     Using DMD                    !
        !--------------------------------------------------!
        IF (HO_Form .EQ. 3) THEN
          CALL INPUT_fg_DMD(DMD_dsets(Current_Database), DMDgsum, dN_x, dN_y, dN_g, L_BCg, W_BCg, B_BCg, C_BCg, rrank_BCg,&
            L_fg_avg_xx, W_fg_avg_xx, B_fg_avg_xx, C_fg_avg_xx, rrank_fg_avg_xx, L_fg_edgV_xx, W_fg_edgV_xx, B_fg_edgV_xx,&
            C_fg_edgV_xx, rrank_fg_edgV_xx, L_fg_avg_yy, W_fg_avg_yy, B_fg_avg_yy, C_fg_avg_yy, rrank_fg_avg_yy, L_fg_edgH_yy,&
            W_fg_edgH_yy, B_fg_edgH_yy, C_fg_edgH_yy, rrank_fg_edgH_yy, L_fg_edgV_xy, W_fg_edgV_xy, B_fg_edgV_xy, C_fg_edgV_xy,&
            rrank_fg_edgV_xy, L_fg_edgH_xy, W_fg_edgH_xy, B_fg_edgH_xy, C_fg_edgH_xy, rrank_fg_edgH_xy, xlen_d, ylen_d,&
            Tini_d, Delx_d, Dely_d, bcT_d, BC_Type_d, Start_Time_d)

          CALL GENERATE_GRIDS(Delx_d, Dely_d, Delx, Dely, xlen_d, ylen_d, xlen, ylen, dN_x, dN_y, N_x, N_y,&
            Sim_Grid_Avg, Sim_Grid_EdgV, Sim_Grid_EdgH, Sim_Grid_Bnds, Dat_Grid_Avg, Dat_Grid_EdgV, Dat_Grid_EdgH,&
            Dat_Grid_Bnds, GMap_xyAvg, GMap_xyEdgV, GMap_xyEdgH, GMap_xyBnds, VMap_xyAvg, VMap_xyEdgV, VMap_xyEdgH,&
            VMap_xyBnds)

          gsum = DMDgsum
          Status = nf90_put_var(outID,DMDgsum_ID,(/DMDgsum/),(/Current_Database/),(/1/))
          CALL HANDLE_ERR(Status)

        !--------------------------------------------------!
        !                     Using POD                    !
        !--------------------------------------------------!
        ELSE IF (HO_Form .EQ. 1) THEN
          CALL INPUT_fg_POD(POD_dsets(Current_Database), PODgsum, POD_err(Current_Database), dN_x, dN_y, dN_m, dN_g, dN_t,&
            C_BCg, S_BCg, U_BCg, V_BCg, rrank_BCg, C_fg_avg_xx, S_fg_avg_xx, U_fg_avg_xx, V_fg_avg_xx, rrank_fg_avg_xx,&
            C_fg_edgV_xx, S_fg_edgV_xx, U_fg_edgV_xx, V_fg_edgV_xx, rrank_fg_edgV_xx, C_fg_avg_yy, S_fg_avg_yy, U_fg_avg_yy,&
            V_fg_avg_yy, rrank_fg_avg_yy, C_fg_edgH_yy, S_fg_edgH_yy, U_fg_edgH_yy, V_fg_edgH_yy, rrank_fg_edgH_yy, C_fg_edgV_xy,&
            S_fg_edgV_xy, U_fg_edgV_xy, V_fg_edgV_xy, rrank_fg_edgV_xy, C_fg_edgH_xy, S_fg_edgH_xy, U_fg_edgH_xy, V_fg_edgH_xy,&
            rrank_fg_edgH_xy, tlen_d, xlen_d, ylen_d, Tini_d, Delt_d, Delx_d, Dely_d, bcT_d, BC_Type_d, Start_Time_d)

          CALL GENERATE_GRIDS(Delx_d, Dely_d, Delt_d, Delx, Dely, Delt, xlen_d, ylen_d, Start_Time_d, xlen, ylen, Start_Time,&
            dN_x, dN_y, dN_t, N_x, N_y, N_t, Sim_Grid_Avg, Sim_Grid_EdgV, Sim_Grid_EdgH, Sim_Grid_Bnds, Dat_Grid_Avg,&
            Dat_Grid_EdgV, Dat_Grid_EdgH,Dat_Grid_Bnds, Sim_TGrid, Dat_TGrid, GMap_xyAvg, GMap_xyEdgV, GMap_xyEdgH, GMap_xyBnds,&
            VMap_xyAvg, VMap_xyEdgV, VMap_xyEdgH, VMap_xyBnds, TMap)

          gsum = PODgsum
          Status = nf90_put_var(outID,PODgsum_ID,(/PODgsum/),(/Current_Database/),(/1/))
          CALL HANDLE_ERR(Status)

        END IF
        ! WRITE(*,*) 'Successfully loaded database ',Current_Database

        !--------------------------------------------------!
        !   Writing ranks of each dataset to output file   !
        !--------------------------------------------------!
        IF (gsum .EQ. 1) THEN
          Status = nf90_put_var(outID,rrank_BCg_ID,rrank_BCg,(/1,Current_Database/),(/1,1/))
          CALL HANDLE_ERR(Status)
          Status = nf90_put_var(outID,rrank_fg_avg_xx_ID,rrank_fg_avg_xx,(/1,Current_Database/),(/1,1/))
          CALL HANDLE_ERR(Status)
          Status = nf90_put_var(outID,rrank_fg_EdgV_xx_ID,rrank_fg_EdgV_xx,(/1,Current_Database/),(/1,1/))
          CALL HANDLE_ERR(Status)
          Status = nf90_put_var(outID,rrank_fg_avg_yy_ID,rrank_fg_avg_yy,(/1,Current_Database/),(/1,1/))
          CALL HANDLE_ERR(Status)
          Status = nf90_put_var(outID,rrank_fg_EdgH_yy_ID,rrank_fg_EdgH_yy,(/1,Current_Database/),(/1,1/))
          CALL HANDLE_ERR(Status)
          Status = nf90_put_var(outID,rrank_fg_EdgV_xy_ID,rrank_fg_EdgV_xy,(/1,Current_Database/),(/1,1/))
          CALL HANDLE_ERR(Status)
          Status = nf90_put_var(outID,rrank_fg_EdgH_xy_ID,rrank_fg_EdgH_xy,(/1,Current_Database/),(/1,1/))
          CALL HANDLE_ERR(Status)
        ELSE
          Status = nf90_put_var(outID,rrank_BCg_ID,rrank_BCg,(/1,Current_Database/),(/N_g,1/))
          CALL HANDLE_ERR(Status)
          Status = nf90_put_var(outID,rrank_fg_avg_xx_ID,rrank_fg_avg_xx,(/1,Current_Database/),(/N_g,1/))
          CALL HANDLE_ERR(Status)
          Status = nf90_put_var(outID,rrank_fg_EdgV_xx_ID,rrank_fg_EdgV_xx,(/1,Current_Database/),(/N_g,1/))
          CALL HANDLE_ERR(Status)
          Status = nf90_put_var(outID,rrank_fg_avg_yy_ID,rrank_fg_avg_yy,(/1,Current_Database/),(/N_g,1/))
          CALL HANDLE_ERR(Status)
          Status = nf90_put_var(outID,rrank_fg_EdgH_yy_ID,rrank_fg_EdgH_yy,(/1,Current_Database/),(/N_g,1/))
          CALL HANDLE_ERR(Status)
          Status = nf90_put_var(outID,rrank_fg_EdgV_xy_ID,rrank_fg_EdgV_xy,(/1,Current_Database/),(/N_g,1/))
          CALL HANDLE_ERR(Status)
          Status = nf90_put_var(outID,rrank_fg_EdgH_xy_ID,rrank_fg_EdgH_xy,(/1,Current_Database/),(/N_g,1/))
          CALL HANDLE_ERR(Status)
        END IF

      END IF
    END IF

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
      !--------------------------------------------------!
      !                   HO_FORM = 0                    !
      !                                                  !
      !            SOLVING THE FULL-ORDER RTE            !
      !                                                  !
      ! (Only solving for iterations after RT_start_Its) !
      !--------------------------------------------------!
      IF ((HO_Form .EQ. 0).AND.(RT_Its .GE. RT_start_Its)) THEN

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
          Hg_avg_xy,Hg_edgV_xx,Hg_edgV_xy,Hg_edgH_yy,Hg_edgH_xy,HO_Eg_edgV,HO_Eg_edgH,HO_Eg_avg,HO_Fxg_edgV,HO_Fyg_edgH,&
          HO_E_edgV,HO_E_edgH,HO_E_avg,HO_Fx_edgV,HO_Fy_edgH)

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
        CALL fg_Calc(fg_avg_xx,fg_avg_yy,fg_avg_xy,fg_edgV_xx,fg_edgV_xy,fg_edgH_yy,fg_edgH_xy,Hg_avg_xx,Hg_avg_yy,Hg_avg_xy,&
          Hg_edgV_xx,Hg_edgV_xy,Hg_edgH_yy,Hg_edgH_xy,HO_Eg_avg,HO_Eg_edgV,HO_Eg_edgH,c,Threads)

        !calculating new MGQD boundary factors with the current iterate's intensities
        CALL Cg_Calc(Cg_L,Cg_B,Cg_R,Cg_T,I_edgV,I_edgH,Omega_x,Omega_y,quad_weight,Comp_Unit,BC_Type,Threads)

      !--------------------------------------------------!
      !                   HO_FORM = 1                    !
      !                                                  !
      !    CONSTRUCTING MULTIGROUP QD FACTORS FROM -     !
      !    - TSVD REPRESENTATION                         !
      !--------------------------------------------------!
      ELSE IF (HO_Form .EQ. 1) THEN

        CALL POD_RECONSTRUCT_fg(fg_avg_xx,fg_avg_yy,fg_edgV_xx,fg_edgV_xy,fg_edgH_yy,fg_edgH_xy,C_fg_avg_xx,&
          S_fg_avg_xx,U_fg_avg_xx,V_fg_avg_xx,C_fg_edgV_xx,S_fg_edgV_xx,U_fg_edgV_xx,V_fg_edgV_xx,C_fg_avg_yy,S_fg_avg_yy,&
          U_fg_avg_yy,V_fg_avg_yy,C_fg_edgH_yy,S_fg_edgH_yy,U_fg_edgH_yy,V_fg_edgH_yy,C_fg_edgV_xy,S_fg_edgV_xy,U_fg_edgV_xy,&
          V_fg_edgV_xy,C_fg_edgH_xy,S_fg_edgH_xy,U_fg_edgH_xy,V_fg_edgH_xy,rrank_fg_avg_xx,rrank_fg_edgV_xx,rrank_fg_avg_yy,&
          rrank_fg_edgH_yy,rrank_fg_edgV_xy,rrank_fg_edgH_xy,dN_x,dN_y,dN_g,dN_t,N_x,N_y,N_g,N_t,t,PODgsum,Sim_Grid_Avg,&
          Sim_Grid_EdgV,Sim_Grid_EdgH,Dat_Grid_Avg,Dat_Grid_EdgV,Dat_Grid_EdgH,Sim_TGrid,Dat_TGrid,GMap_xyAvg,GMap_xyEdgV,&
          GMap_xyEdgH,VMap_xyAvg,VMap_xyEdgV,VMap_xyEdgH,TMap,Threads)

        CALL POD_RECONSTRUCT_BCg(Cg_L,Cg_B,Cg_R,Cg_T,C_BCg,S_BCg,U_BCg,V_BCg,rrank_BCg,dN_x,dN_y,dN_g,dN_t,N_x,N_y,N_g,N_t,&
          t,PODgsum,Sim_Grid_Bnds,Dat_Grid_Bnds,Sim_TGrid,Dat_TGrid,GMap_xyBnds,VMap_xyBnds,TMap,Threads)

        IF (Direc_Diff .EQ. 1) THEN
          fg_edgV_xy = 0d0
          fg_edgH_xy = 0d0
        END IF

      !--------------------------------------------------!
      !                   HO_FORM = 2                    !
      !                                                  !
      !          SKIP CALCULATION OF QD FACTORS          !
      !                                                  !
      ! (QD factors initialized to 1/3)                  !
      !--------------------------------------------------!
      ELSE IF (HO_Form .EQ. 2) THEN

        ! CONTINUE !do nothing
        ! write(*,*) maxit_RTE

      !--------------------------------------------------!
      !                   HO_FORM = 3                    !
      !                                                  !
      !    CONSTRUCTING MULTIGROUP QD FACTORS FROM -     !
      !    - DMD EXPANSION                               !
      !--------------------------------------------------!
      ELSE IF (HO_Form .EQ. 3) THEN
        DMD_Time = Time - Start_Time_d

        CALL DMD_RECONSTRUCT_fg(fg_avg_xx, fg_avg_yy, fg_edgV_xx, fg_edgV_xy, fg_edgH_yy, fg_edgH_xy,&
          L_fg_avg_xx, B_fg_avg_xx, W_fg_avg_xx, C_fg_avg_xx, L_fg_edgV_xx, B_fg_edgV_xx, W_fg_edgV_xx, C_fg_edgV_xx,&
          L_fg_avg_yy, B_fg_avg_yy, W_fg_avg_yy, C_fg_avg_yy, L_fg_edgH_yy, B_fg_edgH_yy, W_fg_edgH_yy, C_fg_edgH_yy,&
          L_fg_edgV_xy, B_fg_edgV_xy, W_fg_edgV_xy, C_fg_edgV_xy, L_fg_edgH_xy, B_fg_edgH_xy, W_fg_edgH_xy, C_fg_edgH_xy,&
          rrank_fg_avg_xx, rrank_fg_edgV_xx, rrank_fg_avg_yy, rrank_fg_edgH_yy, rrank_fg_edgV_xy, rrank_fg_edgH_xy,&
          dN_x, dN_y, dN_g, N_x, N_y, N_g, DMD_Time, DMDgsum, Sim_Grid_Avg, Sim_Grid_EdgV, Sim_Grid_EdgH, Dat_Grid_Avg,&
          Dat_Grid_EdgV, Dat_Grid_EdgH, GMap_xyAvg, GMap_xyEdgV, GMap_xyEdgH, VMap_xyAvg, VMap_xyEdgV, VMap_xyEdgH, Threads)

        CALL DMD_RECONSTRUCT_BCg(Cg_L, Cg_B, Cg_R, Cg_T, L_BCg, B_BCg, W_BCg, C_BCg, rrank_BCg, dN_x, dN_y, dN_g, N_x, N_y, N_g,&
          DMD_Time, DMDgsum, Sim_Grid_Bnds, Dat_Grid_Bnds, GMap_xyBnds, VMap_xyBnds, Threads)

          ! write(*,*) fg_avg_xx
          ! stop

      !--------------------------------------------------!
      !                   HO_FORM = 4                    !
      !                                                  !
      !    READING IN MULTIGROUP QD TENSOR GIVEN BY -    !
      !    - SOME PREVIOUS PROBLEM SOLUTION              !
      !--------------------------------------------------!
      ELSE IF (HO_Form .EQ. 4) THEN

        CALL NF_OPEN_FILE(qdf_infileID,qdf_infile,'old','r')

        CALL NF_INQ_VAR_3D(qdf_infileID, 'fg_avg_xx', fg_avg_xx, (/1,1,1,t/), (/N_x, N_y, N_g, 1/))
        CALL NF_INQ_VAR_3D(qdf_infileID, 'fg_avg_yy', fg_avg_yy, (/1,1,1,t/), (/N_x, N_y, N_g, 1/))
        CALL NF_INQ_VAR_3D(qdf_infileID, 'fg_edgV_xx', fg_edgV_xx, (/1,1,1,t/), (/N_x+1, N_y, N_g, 1/))
        CALL NF_INQ_VAR_3D(qdf_infileID, 'fg_edgH_yy', fg_edgH_yy, (/1,1,1,t/), (/N_x, N_y+1, N_g, 1/))

        IF (qdfin_fxy_flag .EQ. 1) THEN
          CALL NF_INQ_VAR_3D(qdf_infileID, 'fg_edgV_xy', fg_edgV_xy, (/1,1,1,t/), (/N_x+1, N_y, N_g, 1/))
          CALL NF_INQ_VAR_3D(qdf_infileID, 'fg_edgH_xy', fg_edgH_xy, (/1,1,1,t/), (/N_x, N_y+1, N_g, 1/))
        END IF

        IF (qdfin_bctype .EQ. 'qd') THEN
          CALL NF_INQ_VAR_2D(qdf_infileID, 'Cg_L', Cg_L, (/1,1,t/), (/N_y, N_g, 1/))
          CALL NF_INQ_VAR_2D(qdf_infileID, 'Cg_B', Cg_B, (/1,1,t/), (/N_x, N_g, 1/))
          CALL NF_INQ_VAR_2D(qdf_infileID, 'Cg_R', Cg_R, (/1,1,t/), (/N_y, N_g, 1/))
          CALL NF_INQ_VAR_2D(qdf_infileID, 'Cg_T', Cg_T, (/1,1,t/), (/N_x, N_g, 1/))
        END IF

        CALL NF_CLOSE_FILE(qdf_infileID)

      END IF

      IF ( run_type .EQ. 'tr_no_qd' ) THEN
        !solve the MEB equation with the MGQD solution to find new Temp
        CALL MEB_SOLVE(Temp,HO_Eg_avg,Bg,KapE,Temp_old,Delt,cV,Kap0,Comp_Unit)
        !update material properties with new Temp
        CALL MATERIAL_UPDATE(RT_Src,MGQD_Src,KapE,KapB,KapR,Bg,Temp,Comp_Unit,Nu_g,Threads,Mat,Kappa_Mult)

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
            Pold_T,Eg_avg_old,Fxg_edgV_old,Fyg_edgH_old,MGQD_Kits,Fdt_Weight,Fin_Weight,Ein_Weight)

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
              DR_xy,DB_xy,DT_xy,PL,PR,PB,PT,Fdt_Weight)

            E_avg = MGQD_E_avg
            E_edgV = MGQD_E_edgV
            E_edgH = MGQD_E_edgH
            Fx_edgV = MGQD_Fx_edgV
            Fy_edgH = MGQD_Fy_edgH

            !solve the nonlinear EGP with newton iterations
            CALL EGP_FV_NEWT(E_avg,E_edgV,E_edgH,Temp,KapE_Bar,Fx_edgV,Fy_edgH,GQD_Src,KapE_Bar_dT,EGP_Its,Deltas,dres,&
              Temp_MGQDold2,KapE_bar_MGQDold,Theta,Delt,Delx,Dely,A,Cb_L,Cb_B,Cb_R,Cb_T,E_in_L,E_in_B,E_in_R,E_in_T,F_in_L,&
              F_in_B,F_in_R,F_in_T,DC_xx,DL_xx,DR_xx,DC_yy,DB_yy,DT_yy,DL_xy,DR_xy,DB_xy,DT_xy,PL,PR,PB,PT,Gold_Hat,&
              Rhat_old,Kap0,cv,Comp_Unit,Chi,line_src,E_Bound_Low,T_Bound_Low,Conv_gr1,Conv_gr2,Maxit_GLOQD,MGQD_Its,&
              Use_Line_Search,Use_Safety_Search,Res_Calc,kapE_dT_flag,GQD_Kits,Fin_Weight,Ein_Weight)

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
            IF (its_out .EQ. 1) THEN !checking if iteration counts are to be output
              Status = nf90_put_var(outID,GQD_ItCount_ID,(/EGP_Its/),(/MGQD_Its,RT_Its,t/),(/1,1,1/))
              CALL HANDLE_ERR(Status)
              Status = nf90_put_var(outID,GQD_Kits_ID,GQD_Kits,(/1,MGQD_Its,RT_Its,t/),(/EGP_Its,1,1,1/))
              CALL HANDLE_ERR(Status)
            END IF

          ELSE
            !solve the MEB equation with the MGQD solution to find new Temp
            CALL MEB_SOLVE(Temp,Eg_avg,Bg,KapE,Temp_old,Delt,cV,Kap0,Comp_Unit)

          END IF

          !update material properties with new Temp
          CALL MATERIAL_UPDATE(RT_Src,MGQD_Src,KapE,KapB,KapR,Bg,Temp,Comp_Unit,Nu_g,Threads,Mat,Kappa_Mult)
          IF (run_type .EQ. 'fld') CALL FLD_COEFFICIENT_CALC(KapR,Eg_avg,Eg_edgV,Eg_edgH,Delx,Dely,Threads)

          !check convergence of MGQD solution (in grey form) for E, T
          Tconv = CONVERGENCE(Conv_Type,conv_lo,Temp,Temp_MGQDold,Temp_MGQDold2,MGQD_Its,&
            MGQD_Tnorm,MGQD_Trho)
          Econv = CONVERGENCE(Conv_Type,conv_lo,MGQD_E_avg,MGQD_E_avg_MGQDold,MGQD_E_avg_MGQDold2,MGQD_Its,&
            MGQD_Enorm,MGQD_Erho)
          MGQD_conv = Tconv.AND.Econv !if both T and E are converged, the MGQD iterations have successfully converged

          !writing MGQD convergence status (norms and spectral radii) to output file
          IF (conv_out .EQ. 1) THEN !checking if convergence info is to be output
            Status = nf90_put_var(outID,MGQD_Tnorm_ID,(/MGQD_Tnorm/),(/MGQD_Its,RT_Its,t/),(/1,1,1/))
            CALL HANDLE_ERR(Status)
            Status = nf90_put_var(outID,MGQD_Enorm_ID,(/MGQD_Enorm/),(/MGQD_Its,RT_Its,t/),(/1,1,1/))
            CALL HANDLE_ERR(Status)
            Status = nf90_put_var(outID,MGQD_Trho_ID,(/MGQD_Trho/),(/MGQD_Its,RT_Its,t/),(/1,1,1/))
            CALL HANDLE_ERR(Status)
            Status = nf90_put_var(outID,MGQD_Erho_ID,(/MGQD_Erho/),(/MGQD_Its,RT_Its,t/),(/1,1,1/))
            CALL HANDLE_ERR(Status)
          END IF

          !preparing for next iteration by moving current solution -> last iterate solution
          MGQD_E_avg_MGQDold2 = MGQD_E_avg_MGQDold
          MGQD_E_avg_MGQDold = MGQD_E_avg
          Temp_MGQDold2 = Temp_MGQDold
          Temp_MGQDold = Temp

          !writing current iterate to terminal
          ! write(*,*) 'MGQD:     ',MGQD_Its, EGP_Its, MGQD_Tnorm, MGQD_Enorm

        END DO

        !writing the count of MGQD iterations to output file
        IF (its_out .EQ. 1) THEN !checking if iteration counts are to be output
          Status = nf90_put_var(outID,MGQD_ItCount_ID,(/MGQD_Its/),(/RT_Its,t/),(/1,1/))
          CALL HANDLE_ERR(Status)
          Status = nf90_put_var(outID,MGQD_Kits_ID,MGQD_Kits,(/1,RT_Its,t/),(/N_g,1,1/))
          CALL HANDLE_ERR(Status)
        END IF

      END IF

      Tconv = CONVERGENCE(Conv_Type,conv_ho,Temp,Temp_RTold,Temp_RTold2,RT_Its,TR_Tnorm,TR_Trho)
      Econv = CONVERGENCE(Conv_Type,conv_ho,HO_E_avg,HO_E_avg_RTold,HO_E_avg_RTold2,RT_Its,TR_Enorm,TR_Erho)
      RT_conv = Tconv.AND.Econv !if both T and E are converged, the outer/RTE iterations have successfully converged

      !writing RTE convergence status (norms and spectral radii) to output file
      IF (conv_out .EQ. 1) THEN !checking if convergence info is to be output
        Status = nf90_put_var(outID,RT_Tnorm_ID,(/TR_Tnorm/),(/RT_Its,t/),(/1,1/))
        CALL HANDLE_ERR(Status)
        Status = nf90_put_var(outID,RT_Enorm_ID,(/TR_Enorm/),(/RT_Its,t/),(/1,1/))
        CALL HANDLE_ERR(Status)
        Status = nf90_put_var(outID,RT_Trho_ID,(/TR_Trho/),(/RT_Its,t/),(/1,1/))
        CALL HANDLE_ERR(Status)
        Status = nf90_put_var(outID,RT_Erho_ID,(/TR_Erho/),(/RT_Its,t/),(/1,1/))
        CALL HANDLE_ERR(Status)
      END IF

      !preparing for next iteration by moving current solution -> last iterate solution
      HO_E_avg_RTold2 = HO_E_avg_RTold
      HO_E_avg_RTold = HO_E_avg
      Temp_RTold2 = Temp_RTold
      Temp_RTold = Temp

      !writing current iterate to terminal
      ! write(*,*) RT_Its, MGQD_Its, TR_Tnorm, TR_Enorm

    END DO

    !writing finished time step to terminal
    IF ( run_type .EQ. 'tr_no_qd' ) THEN
      write(*,*) Time, RT_Its, TR_Tnorm, TR_Enorm
    ELSE IF ( run_type .EQ. 'mlqd' ) THEN
      write(*,*) Time, RT_Its, TR_Tnorm, TR_Enorm
    ELSE IF ((run_type .EQ. 'mg_pod').AND.(POD_Type .EQ. 'fg')) THEN
      write(*,*) Time, MGQD_Its, MGQD_Tnorm, MGQD_Enorm
    ELSE IF ((run_type .EQ. 'mg_dmd').AND.(DMD_Type .EQ. 'fg')) THEN
      write(*,*) Time, MGQD_Its, MGQD_Tnorm, MGQD_Enorm
    ELSE IF (run_type .EQ. 'qdf_in') THEN
      write(*,*) Time, MGQD_Its, MGQD_Tnorm, MGQD_Enorm
    ELSE IF ((run_type .EQ. 'diff').OR.(run_type .EQ. 'fld').OR.(run_type .EQ. 'p1').OR.(run_type .EQ. 'p13')) THEN
      write(*,*) Time, MGQD_Its, MGQD_Tnorm, MGQD_Enorm
    ELSE
    END IF

    !calculating, writing speed of radiation & temperature wave fronts
    E_xWvSpd = XWAVE_SPEED_edg(E_avg, E_edgV, E_avg_old, Dely, Delx, Delt, ylen, N_y, N_x)
    CALL NF_PUT_t_VAR(outID, E_XWvSpeed_ID, E_xWvSpd, t)
    T_xWvSpd = XWAVE_SPEED_cnt(Temp, Temp_old, Dely, Delx, Delt, ylen, N_y, N_x)
    CALL NF_PUT_t_VAR(outID, Temp_XWvSpeed_ID, T_xWvSpd, t)

    !writing the count of outer/RTE iterations to output file
    IF (its_out .EQ. 1) THEN !checking if iteration counts are to be output
      Status = nf90_put_var(outID,RT_ItCount_ID,(/RT_Its/),(/t/),(/1/))
      CALL HANDLE_ERR(Status)
    END IF

    CALL TIMESTEP_OUTS(outID,t,out_freq,I_out,HO_Eg_out,HO_Fg_out,HO_E_out,HO_F_out,Eg_out,Fg_out,MGQD_E_out,MGQD_F_out,&
      QDfg_out,E_out,F_out,D_out,old_parms_out,its_out,conv_out,kap_out,Src_out,use_grey,Temp_ID,E_avg_ID,E_edgV_ID,E_edgH_ID,&
      MGQD_E_avg_ID,MGQD_E_edgV_ID,MGQD_E_edgH_ID,HO_E_avg_ID,HO_E_edgV_ID,HO_E_edgH_ID,Fx_edgV_ID,Fy_edgH_ID,&
      MGQD_Fx_edgV_ID,MGQD_Fy_edgH_ID,HO_Fx_edgV_ID,HO_Fy_edgH_ID,Eg_avg_ID,Eg_edgV_ID,Eg_edgH_ID,HO_Eg_avg_ID,&
      HO_Eg_edgV_ID,HO_Eg_edgH_ID,Fxg_edgV_ID,Fyg_edgH_ID,HO_Fxg_edgV_ID,HO_Fyg_edgH_ID,I_avg_ID,I_edgV_ID,I_edgH_ID,&
      KapE_Bar_ID,KapB_ID,KapE_ID,KapR_ID,Bg_ID,Cg_L_ID,Cg_B_ID,Cg_R_ID,Cg_T_ID,Eg_in_L_ID,Eg_in_B_ID,Eg_in_R_ID,&
      Eg_in_T_ID,Fg_in_L_ID,Fg_in_B_ID,Fg_in_R_ID,Fg_in_T_ID,Cb_L_ID,Cb_B_ID,Cb_R_ID,Cb_T_ID,E_in_L_ID,E_in_B_ID,&
      E_in_R_ID,E_in_T_ID,F_in_L_ID,F_in_B_ID,F_in_R_ID,F_in_T_ID,fg_avg_xx_ID,fg_avg_yy_ID,fg_avg_xy_ID,fg_edgV_xx_ID,&
      fg_edgV_xy_ID,fg_edgH_yy_ID,fg_edgH_xy_ID,DC_xx_ID,DL_xx_ID,DR_xx_ID,DC_yy_ID,DB_yy_ID,DT_yy_ID,DL_xy_ID,&
      DB_xy_ID,DR_xy_ID,DT_xy_ID,G_old_ID,Pold_L_ID,Pold_B_ID,Pold_R_ID,Pold_T_ID,Gold_hat_ID,Rhat_old_ID,PL_ID,&
      PB_ID,PR_ID,PT_ID, Temp,E_avg,E_edgV,E_edgH,MGQD_E_avg,MGQD_E_edgV,MGQD_E_edgH,HO_E_avg,HO_E_edgV,HO_E_edgH,&
      Fx_edgV,Fy_edgH,MGQD_Fx_edgV,MGQD_Fy_edgH,HO_Fx_edgV,HO_Fy_edgH,Eg_avg,Eg_edgV,Eg_edgH,HO_Eg_avg,HO_Eg_edgV,&
      HO_Eg_edgH,Fxg_edgV,Fyg_edgH,HO_Fxg_edgV,HO_Fyg_edgH,I_avg,I_edgV,I_edgH,KapE_Bar,KapB,KapE,KapR,Bg,Cg_L,Cg_B,&
      Cg_R,Cg_T,Eg_in_L,Eg_in_B,Eg_in_R,Eg_in_T,Fg_in_L,Fg_in_B,Fg_in_R,Fg_in_T,Cb_L,Cb_B,Cb_R,Cb_T,E_in_L,E_in_B,&
      E_in_R,E_in_T,F_in_L,F_in_B,F_in_R,F_in_T,fg_avg_xx,fg_avg_yy,fg_avg_xy,fg_edgV_xx,fg_edgV_xy,fg_edgH_yy,fg_edgH_xy,DC_xx,&
      DL_xx,DR_xx,DC_yy,DB_yy,DT_yy,DL_xy,DB_xy,DR_xy,DT_xy,G_old,Pold_L,Pold_B,Pold_R,Pold_T,Gold_hat,Rhat_old,PL,PB,PR,PT)

    !Creating, writing restart file on each 'restart_freq'-th time step
    IF (restart_freq .GT. 0) THEN
      IF (MOD(t, restart_freq) .EQ. 0) CALL RESTART_OUT(Time, Temp, I_crn, fg_avg_xx, fg_avg_yy, fg_edgV_xx, fg_edgV_xy,&
        fg_edgH_yy, fg_edgH_xy, RT_Src, MGQD_Src, GQD_Src, E_avg, Eg_avg, Eg_edgV, Eg_edgH, Fxg_edgV, Fyg_edgH, KapE_Bar,&
        KapE, KapR,Cg_L, Cg_B, Cg_R, Cg_T,restart_outfile, resf_unit)
    END IF

    IF (Time .GE. tlen) EXIT
  END DO

END SUBROUTINE TRT_MLQD_ALGORITHM

!==================================================================================================================================!
!
!==================================================================================================================================!
SUBROUTINE Tgen_QDf(Omega_x, Omega_y, quad_weight, Nu_g, Kappa_Mult, c, cV, h, pi, Kap0, erg, Comp_Unit,&
  N_m, Mat, threads, BC_Type, Res_Calc, TfileName, outfile, enrgy_strc,&
  I_out, HO_Eg_out, HO_Fg_out, HO_E_out, HO_F_out, QDfg_out, quadrature, run_type)

  !---------------Solution Parameters----------------!
  REAL*8,INTENT(IN):: Omega_x(:), Omega_y(:), quad_weight(:), Nu_g(:), Kappa_Mult(:)
  REAL*8,INTENT(IN):: c, cV, h, pi, Kap0, erg
  REAL*8,INTENT(IN):: Comp_Unit
  INTEGER,INTENT(IN):: N_m, Mat(:,:)
  INTEGER,INTENT(IN):: Threads, BC_Type(:)
  LOGICAL,INTENT(IN):: Res_Calc
  CHARACTER(*),INTENT(IN):: TfileName, outfile, enrgy_strc, quadrature, run_type
  INTEGER,INTENT(IN):: I_out, HO_Eg_out, HO_Fg_out, HO_E_out, HO_F_out, QDfg_out

  REAL*8,ALLOCATABLE:: Temp(:,:)
  REAL*8,ALLOCATABLE:: fg_avg_xx(:,:,:), fg_avg_yy(:,:,:), fg_avg_xy(:,:,:)
  REAL*8,ALLOCATABLE:: fg_edgV_xx(:,:,:), fg_edgV_xy(:,:,:)
  REAL*8,ALLOCATABLE:: fg_edgH_yy(:,:,:), fg_edgH_xy(:,:,:)
  REAL*8,ALLOCATABLE:: Cg_L(:,:), Cg_B(:,:), Cg_R(:,:), Cg_T(:,:)

  !---------------Material Properties----------------!
  REAL*8:: Delt, xlen, ylen, tlen
  REAL*8,ALLOCATABLE:: xpts_avg(:), xpts_edgV(:), ypts_avg(:), ypts_edgH(:), tpts(:)
  REAL*8,ALLOCATABLE:: A(:,:), Delx(:), Dely(:)
  REAL*8,ALLOCATABLE:: Bg(:,:,:), KapE(:,:,:), KapB(:,:,:), KapR(:,:,:)
  REAL*8,ALLOCATABLE:: KapE_old(:,:,:), KapR_old(:,:,:), Temp_old(:,:)
  REAL*8,ALLOCATABLE:: MGQD_Src(:,:,:), MGQD_Src_old(:,:,:), RT_Src(:,:,:,:)

  !--------------Radiation Intensities---------------!
  REAL*8,ALLOCATABLE:: I_crn_old(:,:,:,:)
  REAL*8,ALLOCATABLE:: I_avg(:,:,:,:), I_edgV(:,:,:,:), I_edgH(:,:,:,:)
  REAL*8,ALLOCATABLE:: I_crn(:,:,:,:), Ic_edgV(:,:,:,:), Ic_edgH(:,:,:,:)

  !--------------High-Order Quantities---------------!
  REAL*8,ALLOCATABLE:: Hg_avg_xx(:,:,:), Hg_avg_yy(:,:,:), Hg_avg_xy(:,:,:)
  REAL*8,ALLOCATABLE:: Hg_edgV_xx(:,:,:), Hg_edgV_xy(:,:,:)
  REAL*8,ALLOCATABLE:: Hg_edgH_yy(:,:,:), Hg_edgH_xy(:,:,:)
  REAL*8,ALLOCATABLE:: HO_Eg_avg(:,:,:), HO_Eg_edgV(:,:,:), HO_Eg_edgH(:,:,:)
  REAL*8,ALLOCATABLE:: HO_Fxg_edgV(:,:,:), HO_Fyg_edgH(:,:,:)
  REAL*8,ALLOCATABLE:: HO_E_avg(:,:), HO_E_edgV(:,:), HO_E_edgH(:,:)
  REAL*8,ALLOCATABLE:: HO_Fx_edgV(:,:), HO_Fy_edgH(:,:)

  !-----------------Misc. Variables------------------!
  REAL*8:: bcT_left, bcT_bottom, bcT_right, bcT_top, Tini
  REAL*8,ALLOCATABLE:: RT_Residual(:,:,:)
  INTEGER,ALLOCATABLE:: RT_ResLoc_x(:,:,:), RT_ResLoc_y(:,:,:)
  INTEGER:: TfileID, outID
  INTEGER:: N_x_ID, N_y_ID, N_m_ID, N_g_ID, N_t_ID, N_edgV_ID, N_edgH_ID, N_xc_ID, N_yc_ID, Quads_ID
  INTEGER:: Boundaries_ID, c_ID, h_ID, pi_ID, erg_ID, Comp_Unit_ID, cv_ID
  INTEGER:: Temp_ID, HO_E_avg_ID
  INTEGER:: HO_E_edgV_ID, HO_E_edgH_ID, HO_Fx_edgV_ID
  INTEGER:: HO_Fy_edgH_ID, HO_Eg_avg_ID, HO_Eg_edgV_ID, HO_Eg_edgH_ID
  INTEGER:: HO_Fxg_edgV_ID, HO_Fyg_edgH_ID, I_avg_ID, I_edgV_ID, I_edgH_ID
  INTEGER:: Cg_L_ID, Cg_B_ID, Cg_R_ID, Cg_T_ID
  INTEGER:: fg_avg_xx_ID, fg_avg_yy_ID, fg_avg_xy_ID, fg_edgV_xx_ID, fg_edgV_xy_ID, fg_edgH_yy_ID, fg_edgH_xy_ID
  INTEGER:: N_x, N_y, N_g, N_t
  INTEGER:: t

  !===========================================================================!
  !                                                                           !
  !     INITIALIZATION                                                        !
  !                                                                           !
  !===========================================================================!
  write(*,*) TfileName
  CALL NF_OPEN_FILE(TfileID,TfileName,'old','r')
  write(*,*) TfileName
  CALL NF_INQ_DIM(TfileID,"N_x",N_x)
  CALL NF_INQ_DIM(TfileID,"N_y",N_y)
  CALL NF_INQ_DIM(TfileID,"N_g",N_g)
  CALL NF_INQ_DIM(TfileID,"N_t",N_t)

  ALLOCATE(Delx(N_x), Dely(N_y))
  ALLOCATE(xpts_avg(N_x), xpts_edgV(N_x+1))
  ALLOCATE(ypts_avg(N_y), ypts_edgH(N_y+1))
  ALLOCATE(tpts(N_t))

  CALL NF_INQ_VAR_0D(TfileID,"tlen",tlen)
  CALL NF_INQ_VAR_0D(TfileID,"xlen",xlen)
  CALL NF_INQ_VAR_0D(TfileID,"ylen",ylen)
  CALL NF_INQ_VAR_0D(TfileID,"Delt",Delt)
  CALL NF_INQ_VAR_1D(TfileID,"Delx",Delx,(/1/),(/N_x/))
  CALL NF_INQ_VAR_1D(TfileID,"Dely",Dely,(/1/),(/N_y/))
  CALL NF_INQ_VAR_0D(TfileID,"bcT_left",bcT_left)
  CALL NF_INQ_VAR_0D(TfileID,"bcT_bottom",bcT_bottom)
  CALL NF_INQ_VAR_0D(TfileID,"bcT_right",bcT_right)
  CALL NF_INQ_VAR_0D(TfileID,"bcT_top",bcT_top)
  CALL NF_INQ_VAR_0D(TfileID,"Tini",Tini)

  CALL NF_INQ_VAR_1D(TfileID,"tpts",tpts,(/1/),(/N_t/))
  CALL NF_INQ_VAR_1D(TfileID,"xpts_avg",xpts_avg,(/1/),(/N_x/))
  CALL NF_INQ_VAR_1D(TfileID,"xpts_edgV",xpts_edgV,(/1/),(/N_x+1/))
  CALL NF_INQ_VAR_1D(TfileID,"ypts_avg",ypts_avg,(/1/),(/N_y/))
  CALL NF_INQ_VAR_1D(TfileID,"ypts_edgH",ypts_edgH,(/1/),(/N_y+1/))

  !calling 'init' routines to allocate space for all needed arrays and initialize values
  CALL MISC_INIT(Delx,Dely,A)
  CALL RT_INIT(I_avg,I_edgV,I_edgH,I_crn,I_crn_old,Ic_edgV,Ic_edgH,Hg_avg_xx,Hg_avg_yy,Hg_avg_xy,Hg_edgV_xx,&
    Hg_edgV_xy,Hg_edgH_yy,Hg_edgH_xy,HO_Eg_avg,HO_Eg_edgV,HO_Eg_edgH,HO_Fxg_edgV,HO_Fyg_edgH,HO_E_avg,&
    HO_E_edgV,HO_E_edgH,HO_Fx_edgV,HO_Fy_edgH,N_y,N_x,N_m,N_g,Tini,comp_unit,nu_g,bcT_left,bcT_right,&
    bcT_top,bcT_bottom,BC_Type,pi,c)
  CALL TEMP_INIT(Temp,RT_Src,MGQD_Src,MGQD_Src_old,KapE,KapB,KapR,KapE_old,KapR_old,Bg,N_y,N_x,N_m,N_g,Tini,&
    Comp_Unit,Nu_g,Temp_Old,Threads,Mat,Kappa_Mult)

  !Allocate Multigroup quasidiffusion tensors
  ALLOCATE(fg_avg_xx(N_x,N_y,N_g))
  ALLOCATE(fg_avg_yy(N_x,N_y,N_g))
  ALLOCATE(fg_avg_xy(N_x,N_y,N_g))
  ALLOCATE(fg_edgV_xx(N_x+1,N_y,N_g))
  ALLOCATE(fg_edgV_xy(N_x+1,N_y,N_g))
  ALLOCATE(fg_edgH_yy(N_x,N_y+1,N_g))
  ALLOCATE(fg_edgH_xy(N_x,N_y+1,N_g))
  !Allocate MGQD boundary factors
  ALLOCATE(Cg_L(N_y,N_g))
  ALLOCATE(Cg_B(N_x,N_g))
  ALLOCATE(Cg_R(N_y,N_g))
  ALLOCATE(Cg_T(N_x,N_g))

  !Create output file and initialize all arrays/ values
  CALL OUTFILE_INIT_Tgen(outID, N_x_ID, N_y_ID, N_m_ID, N_g_ID, N_t_ID, N_edgV_ID, N_edgH_ID, N_xc_ID, N_yc_ID, Quads_ID,&
    Boundaries_ID, c_ID, h_ID, pi_ID, erg_ID, Comp_Unit_ID, cv_ID, c, h, pi, erg, Comp_Unit, cv, xlen, ylen, Delx, Dely, tlen,&
    Delt,bcT_left, bcT_bottom, bcT_right, bcT_top, Tini, N_x, N_y, N_m, N_g, N_t, threads, BC_type, outfile, quadrature,&
    enrgy_strc, I_out, HO_Eg_out, HO_Fg_out, HO_E_out, HO_F_out, QDfg_out, xpts_avg, xpts_edgV, ypts_avg, ypts_edgH, tpts,&
    run_type, Temp_ID, HO_E_avg_ID, HO_E_edgV_ID, HO_E_edgH_ID, HO_Fx_edgV_ID, HO_Fy_edgH_ID, HO_Eg_avg_ID, HO_Eg_edgV_ID,&
    HO_Eg_edgH_ID, HO_Fxg_edgV_ID, HO_Fyg_edgH_ID, I_avg_ID, I_edgV_ID, I_edgH_ID, Cg_L_ID, Cg_B_ID, Cg_R_ID, Cg_T_ID,&
    fg_avg_xx_ID, fg_avg_yy_ID, fg_avg_xy_ID, fg_edgV_xx_ID, fg_edgV_xy_ID, fg_edgH_yy_ID, fg_edgH_xy_ID)

  !===========================================================================!
  !                                                                           !
  !     PROBLEM SOLVE (BEGIN TIME STEP LOOP)                                  !
  !                                                                           !
  !===========================================================================!
  DO t=1,N_t
    write(*,*) 'Time: ',tpts(t)

    CALL NF_INQ_VAR_2D(TfileID, 'Temperature', Temp, (/1,1,t/), (/N_x, N_y, 1/)) !reading Temp at time t from input file
    CALL NF_PUT_t_VAR(outID,Temp_ID,Temp,t) !writing same Temp to output file

    !calculating material properties at T=Temp
    CALL MATERIAL_UPDATE(RT_Src,MGQD_Src,KapE,KapB,KapR,Bg,Temp,Comp_Unit,Nu_g,Threads,Mat,Kappa_Mult)

    !solving the RTE
    I_crn_old = I_crn
    CALL TRANSPORT_SCB(I_avg,I_edgV,I_edgH,I_crn,Ic_edgV,Ic_edgH,Omega_x,Omega_y,Delx,Dely,A,&
      KapE,RT_Src,I_crn_old,c,Delt,Threads,RT_Residual,RT_ResLoc_x,RT_ResLoc_y,Res_Calc)
    !Output intensities
    IF (I_out .EQ. 1) THEN
      CALL NF_PUT_t_VAR(outID,I_avg_ID,I_avg,t)
      CALL NF_PUT_t_VAR(outID,I_edgV_ID,I_edgV,t)
      CALL NF_PUT_t_VAR(outID,I_edgH_ID,I_edgH,t)
    END IF

    !calculating low-order quantities from the high-order intensities
    CALL COLLAPSE_INTENSITIES(Threads,I_avg,I_edgV,I_edgH,Omega_x,Omega_y,quad_weight,Comp_Unit,Hg_avg_xx,Hg_avg_yy,&
      Hg_avg_xy,Hg_edgV_xx,Hg_edgV_xy,Hg_edgH_yy,Hg_edgH_xy,HO_Eg_edgV,HO_Eg_edgH,HO_Eg_avg,HO_Fxg_edgV,HO_Fyg_edgH,&
      HO_E_edgV,HO_E_edgH,HO_E_avg,HO_Fx_edgV,HO_Fy_edgH)
    !Output total rad energy densities
    IF (HO_E_out .EQ. 1) THEN
      CALL NF_PUT_t_VAR(outID,HO_E_avg_ID,HO_E_avg,t)
      CALL NF_PUT_t_VAR(outID,HO_E_edgV_ID,HO_E_edgV,t)
      CALL NF_PUT_t_VAR(outID,HO_E_edgH_ID,HO_E_edgH,t)
    END IF
    !Output total rad fluxes
    IF (HO_F_out .EQ. 1) THEN
      CALL NF_PUT_t_VAR(outID,HO_Fx_edgV_ID,HO_Fx_edgV,t)
      CALL NF_PUT_t_VAR(outID,HO_Fy_edgH_ID,HO_Fy_edgH,t)
    END IF
    !Output multigroup rad energy densities
    IF (HO_Eg_out .EQ. 1) THEN
      CALL NF_PUT_t_VAR(outID,HO_Eg_avg_ID,HO_Eg_avg,t)
      CALL NF_PUT_t_VAR(outID,HO_Eg_edgV_ID,HO_Eg_edgV,t)
      CALL NF_PUT_t_VAR(outID,HO_Eg_edgH_ID,HO_Eg_edgH,t)
    END IF
    !Output total rad fluxes
    IF (HO_Fg_out .EQ. 1) THEN
      CALL NF_PUT_t_VAR(outID,HO_Fxg_edgV_ID,HO_Fxg_edgV,t)
      CALL NF_PUT_t_VAR(outID,HO_Fyg_edgH_ID,HO_Fyg_edgH,t)
    END IF

    !calculating QD factors from the low-order quantities calculated from RTE solution
    CALL fg_Calc(fg_avg_xx,fg_avg_yy,fg_avg_xy,fg_edgV_xx,fg_edgV_xy,fg_edgH_yy,fg_edgH_xy,Hg_avg_xx,Hg_avg_yy,Hg_avg_xy,&
      Hg_edgV_xx,Hg_edgV_xy,Hg_edgH_yy,Hg_edgH_xy,HO_Eg_avg,HO_Eg_edgV,HO_Eg_edgH,c,Threads)
    !Output QD factors
    CALL NF_PUT_t_VAR(outID,fg_avg_xx_ID,fg_avg_xx,t)
    CALL NF_PUT_t_VAR(outID,fg_avg_yy_ID,fg_avg_yy,t)
    CALL NF_PUT_t_VAR(outID,fg_avg_xy_ID,fg_avg_xy,t)
    CALL NF_PUT_t_VAR(outID,fg_edgV_xx_ID,fg_edgV_xx,t)
    CALL NF_PUT_t_VAR(outID,fg_edgV_xy_ID,fg_edgV_xy,t)
    CALL NF_PUT_t_VAR(outID,fg_edgH_yy_ID,fg_edgH_yy,t)
    CALL NF_PUT_t_VAR(outID,fg_edgH_xy_ID,fg_edgH_xy,t)

    !calculating new MGQD boundary factors with the current iterate's intensities
    CALL Cg_Calc(Cg_L,Cg_B,Cg_R,Cg_T,I_edgV,I_edgH,Omega_x,Omega_y,quad_weight,Comp_Unit,BC_Type,Threads)
    !Output boundary factors
    CALL NF_PUT_t_VAR(outID,Cg_L_ID,Cg_L,t)
    CALL NF_PUT_t_VAR(outID,Cg_B_ID,Cg_B,t)
    CALL NF_PUT_t_VAR(outID,Cg_R_ID,Cg_R,t)
    CALL NF_PUT_t_VAR(outID,Cg_T_ID,Cg_T,t)

  END DO

  !close files
  CALL NF_CLOSE_FILE(TfileID)
  CALL NF_CLOSE_FILE(outID)

END SUBROUTINE Tgen_QDf

END MODULE ALGORITHMS
