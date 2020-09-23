MODULE OUTPUTS

  USE NCDF_IO

  IMPLICIT NONE

CONTAINS

!==================================================================================================================================!
!
!==================================================================================================================================!
SUBROUTINE OUTFILE_INIT(outID,N_x_ID,N_y_ID,N_m_ID,N_g_ID,N_t_ID,N_edgV_ID,N_edgH_ID,N_xc_ID,N_yc_ID,Quads_ID,RT_Its_ID,&
  MGQD_Its_ID,GQD_Its_ID,Norm_Types_ID,MGQD_ResTypes_ID,Boundaries_ID,c_ID,h_ID,pi_ID,erg_ID,Comp_Unit_ID,cv_ID,&
  c,h,pi,erg,Comp_Unit,cv,chi,conv_ho,conv_lo,conv_gr1,conv_gr2,line_src,xlen,ylen,Delx,Dely,tlen,Delt,bcT_left,&
  bcT_bottom,bcT_right,bcT_top,Tini,E_Bound_Low,T_Bound_Low,N_x,N_y,N_m,N_g,N_t,database_gen,use_grey,maxit_RTE,&
  maxit_MLOQD,maxit_GLOQD,conv_type,threads,BC_type,outfile,run_type,kapE_dT_flag,quadrature,enrgy_strc,Theta,&
  Use_Line_Search,Use_Safety_Search)

  INTEGER,INTENT(OUT):: outID
  INTEGER,INTENT(OUT):: N_x_ID, N_y_ID, N_m_ID, N_g_ID, N_t_ID, N_edgV_ID, N_edgH_ID, N_xc_ID, N_yc_ID, Quads_ID
  INTEGER,INTENT(OUT):: RT_Its_ID, MGQD_Its_ID, GQD_Its_ID, Norm_Types_ID, MGQD_ResTypes_ID, Boundaries_ID
  INTEGER,INTENT(OUT):: c_ID, h_ID, pi_ID, erg_ID, Comp_Unit_ID, cv_ID

  REAL*8,INTENT(IN):: c, h, pi, erg, Comp_Unit, cv
  REAL*8,INTENT(IN):: chi, conv_ho, conv_lo, conv_gr1, conv_gr2, line_src, xlen, ylen, Delx(:), Dely(:), tlen, Delt, Theta
  REAL*8,INTENT(IN):: bcT_left, bcT_bottom, bcT_right, bcT_top, Tini, E_Bound_Low, T_Bound_Low
  INTEGER,INTENT(IN):: N_x, N_y, N_m, N_g, N_t
  INTEGER,INTENT(IN):: database_gen, use_grey, maxit_RTE, maxit_MLOQD, maxit_GLOQD, conv_type, threads, BC_type(:)
  CHARACTER(*),INTENT(IN):: outfile, run_type, quadrature, enrgy_strc
  LOGICAL,INTENT(IN):: Use_Line_Search, Use_Safety_Search, kapE_dT_flag

  INTEGER:: Status
  INTEGER:: xlen_ID, ylen_ID
  INTEGER:: Delx_ID, Dely_ID, tlen_ID, Delt_ID, BC_type_ID, bcT_left_ID, bcT_bottom_ID, bcT_right_ID, bcT_top_ID, Tini_ID
  INTEGER:: Z(0)

  !Opening output file
  CALL NF_OPEN_FILE(outID,outfile,'New')

  !===========================================================================!
  !                                                                           !
  !     DEFINING DIMENSIONS                                                   !
  !                                                                           !
  !===========================================================================!
  CALL NF_DEF_DIM(N_x_ID,outID,'N_x',N_x)
  CALL NF_DEF_DIM(N_y_ID,outID,'N_y',N_y)
  CALL NF_DEF_DIM(N_m_ID,outID,'N_m',N_m)
  CALL NF_DEF_DIM(N_g_ID,outID,'N_g',N_g)
  CALL NF_DEF_DIM(N_t_ID,outID,'N_t',N_t)
  CALL NF_DEF_DIM(N_edgV_ID,outID,'N_edgV',N_x+1)
  CALL NF_DEF_DIM(N_edgH_ID,outID,'N_edgH',N_y+1)
  CALL NF_DEF_DIM(N_xc_ID,outID,'N_xc',N_x*2)
  CALL NF_DEF_DIM(N_yc_ID,outID,'N_yc',N_y*2)
  CALL NF_DEF_DIM(Quads_ID,outID,'Quadrants',4)
  CALL NF_DEF_DIM(RT_Its_ID,outID,'RT_Its-dim',0)
  CALL NF_DEF_DIM(MGQD_Its_ID,outID,'MGQD_Its-dim',0)
  CALL NF_DEF_DIM(GQD_Its_ID,outID,'GQD_Its-dim',0)
  CALL NF_DEF_DIM(Norm_Types_ID,outID,'Norm_Types',5) !inf, 1, 2, l1, l2
  CALL NF_DEF_DIM(MGQD_ResTypes_ID,outID,'MGQD_ResTypes',5)
  CALL NF_DEF_DIM(Boundaries_ID,outID,'Boundaries',4)

  !===========================================================================!
  !                                                                           !
  !     DEFINING CONSTANTS                                                    !
  !                                                                           !
  !===========================================================================!
  CALL NF_DEF_VAR(c_ID,outID,Z,'c--lightspeed','Double')
  CALL NF_DEF_UNIT(outID,c_ID,'cm/sh')

  CALL NF_DEF_VAR(h_ID,outID,Z,'h--Boltzmann_const','Double')
  CALL NF_DEF_UNIT(outID,h_ID,'erg*sh')

  CALL NF_DEF_VAR(cv_ID,outID,Z,'cv--specific_heat','Double')
  CALL NF_DEF_UNIT(outID,cv_ID,'erg/(ev**4 cm**2 sh)')

  CALL NF_DEF_VAR(erg_ID,outID,Z,'ev<->erg_conversion','Double')
  CALL NF_DEF_UNIT(outID,erg_ID,'eV/erg')

  CALL NF_DEF_VAR(pi_ID,outID,Z,'pi','Double')

  CALL NF_DEF_VAR(Comp_Unit_ID,outID,Z,'Computational_Unit','Double')

  !===========================================================================!
  !                                                                           !
  !     DEFINING PROBLEM PARAMETERS & USER INPUTS                             !
  !                                                                           !
  !===========================================================================!
  Status = nf90_put_att(outID,NF90_GLOBAL,'run_type',run_type)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_att(outID,NF90_GLOBAL,'database_gen',database_gen)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_att(outID,NF90_GLOBAL,'use_grey',use_grey)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_att(outID,NF90_GLOBAL,'maxit_RTE',maxit_RTE)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_att(outID,NF90_GLOBAL,'maxit_MLOQD',maxit_MLOQD)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_att(outID,NF90_GLOBAL,'maxit_GLOQD',maxit_GLOQD)
  CALL HANDLE_ERR(Status)

  IF (kapE_dT_flag) THEN
    Status = nf90_put_att(outID,NF90_GLOBAL,'kapE_dT_flag','on')
  ELSE
    Status = nf90_put_att(outID,NF90_GLOBAL,'kapE_dT_flag','off')
  END IF
  CALL HANDLE_ERR(Status)

  Status = nf90_put_att(outID,NF90_GLOBAL,'chi',chi)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_att(outID,NF90_GLOBAL,'conv_ho',conv_ho)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_att(outID,NF90_GLOBAL,'conv_lo',conv_lo)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_att(outID,NF90_GLOBAL,'conv_gr',(/conv_gr1,conv_gr2/))
  CALL HANDLE_ERR(Status)

  Status = nf90_put_att(outID,NF90_GLOBAL,'conv_type',conv_type)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_att(outID,NF90_GLOBAL,'quadrature',quadrature)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_att(outID,NF90_GLOBAL,'threads',threads)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_att(outID,NF90_GLOBAL,'enrgy_strc',enrgy_strc)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_att(outID,NF90_GLOBAL,'BC_type',BC_type)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_att(outID,NF90_GLOBAL,'Theta',Theta)
  CALL HANDLE_ERR(Status)

  IF (Use_Line_Search) THEN
    Status = nf90_put_att(outID,NF90_GLOBAL,'Use_Line_Search','True')
    CALL HANDLE_ERR(Status)
    Status = nf90_put_att(outID,NF90_GLOBAL,'line_src',line_src)
    CALL HANDLE_ERR(Status)

  ELSE
    Status = nf90_put_att(outID,NF90_GLOBAL,'Use_Line_Search','False')
    CALL HANDLE_ERR(Status)

  END IF

  IF (Use_Safety_Search) THEN
    Status = nf90_put_att(outID,NF90_GLOBAL,'Use_Safety_Search','True')
    CALL HANDLE_ERR(Status)
    Status = nf90_put_att(outID,NF90_GLOBAL,'E_Bound_Low',E_Bound_Low)
    CALL HANDLE_ERR(Status)
    Status = nf90_put_att(outID,NF90_GLOBAL,'T_Bound_Low',T_Bound_Low)
    CALL HANDLE_ERR(Status)

  ELSE
    Status = nf90_put_att(outID,NF90_GLOBAL,'Use_Safety_Search','False')
    CALL HANDLE_ERR(Status)

  END IF

  CALL NF_DEF_VAR(xlen_ID,outID,Z,'xlen','Double')
  CALL NF_DEF_UNIT(outID,xlen_ID,'cm')
  CALL NF_DEF_VAR(ylen_ID,outID,Z,'ylen','Double')
  CALL NF_DEF_UNIT(outID,ylen_ID,'cm')
  CALL NF_DEF_VAR(Delx_ID,outID,(/N_x_ID/),'Delx','Double')
  CALL NF_DEF_UNIT(outID,Delx_ID,'cm')
  CALL NF_DEF_VAR(Dely_ID,outID,(/N_y_ID/),'Dely','Double')
  CALL NF_DEF_UNIT(outID,Dely_ID,'cm')
  CALL NF_DEF_VAR(tlen_ID,outID,Z,'tlen','Double')
  CALL NF_DEF_UNIT(outID,tlen_ID,'sh')
  CALL NF_DEF_VAR(Delt_ID,outID,Z,'Delt','Double')
  CALL NF_DEF_UNIT(outID,Delt_ID,'sh')
  CALL NF_DEF_VAR(bcT_left_ID,outID,Z,'bcT_left','Double')
  CALL NF_DEF_UNIT(outID,bcT_left_ID,'ev')
  CALL NF_DEF_VAR(bcT_bottom_ID,outID,Z,'bcT_bottom','Double')
  CALL NF_DEF_UNIT(outID,bcT_bottom_ID,'ev')
  CALL NF_DEF_VAR(bcT_right_ID,outID,Z,'bcT_right','Double')
  CALL NF_DEF_UNIT(outID,bcT_right_ID,'ev')
  CALL NF_DEF_VAR(bcT_top_ID,outID,Z,'bcT_top','Double')
  CALL NF_DEF_UNIT(outID,bcT_top_ID,'ev')
  CALL NF_DEF_VAR(Tini_ID,outID,Z,'Tini','Double')
  CALL NF_DEF_UNIT(outID,Tini_ID,'ev')

  !===========================================================================!
  !                                                                           !
  !     TAKING FILE OUT OF DEFINE MODE                                        !
  !                                                                           !
  !===========================================================================!
  CALL NF_ENDDEF(outID)

  !===========================================================================!
  !                                                                           !
  !     FILLING CONSTANTS & PARAMETERS                                        !
  !                                                                           !
  !===========================================================================!
  Status = nf90_put_var(outID,c_ID,c)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_var(outID,h_ID,h)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_var(outID,pi_ID,pi)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_var(outID,erg_ID,erg)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_var(outID,Comp_Unit_ID,Comp_Unit)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_var(outID,cv_ID,cv)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_var(outID,xlen_ID,xlen)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_var(outID,ylen_ID,ylen)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_var(outID,Delx_ID,Delx,(/1/),(/N_x/))
  CALL HANDLE_ERR(Status)

  Status = nf90_put_var(outID,Dely_ID,Dely,(/1/),(/N_y/))
  CALL HANDLE_ERR(Status)

  Status = nf90_put_var(outID,tlen_ID,tlen)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_var(outID,Delt_ID,Delt)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_var(outID,bcT_left_ID,bcT_left)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_var(outID,bcT_bottom_ID,bcT_bottom)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_var(outID,bcT_right_ID,bcT_right)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_var(outID,bcT_top_ID,bcT_top)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_var(outID,Tini_ID,Tini)
  CALL HANDLE_ERR(Status)

END SUBROUTINE OUTFILE_INIT

!==================================================================================================================================!
!
!==================================================================================================================================!
SUBROUTINE OUTFILE_VARDEFS(outID,Res_Calc,N_x_ID,N_y_ID,N_m_ID,N_g_ID,N_t_ID,N_edgV_ID,N_edgH_ID,N_xc_ID,N_yc_ID,Quads_ID,&
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
  fg_edgV_xy_ID,fg_edgH_yy_ID,fg_edgH_xy_ID,DC_xx_ID,DL_xx_ID,DR_xx_ID,DC_yy_ID,DB_yy_ID,DT_yy_ID,DL_xy_ID,&
  DB_xy_ID,DR_xy_ID,DT_xy_ID,G_old_ID,Pold_L_ID,Pold_B_ID,Pold_R_ID,Pold_T_ID,Gold_hat_ID,Rhat_old_ID,PL_ID,&
  PB_ID,PR_ID,PT_ID,dr_T_ID,dr_B_ID,dr_ML_ID,dr_MB_ID,dr_MR_ID,dr_MT_ID)

  INTEGER,INTENT(IN):: outID
  LOGICAL,INTENT(IN):: Res_Calc
  INTEGER,INTENT(IN):: N_x_ID, N_y_ID, N_m_ID, N_g_ID, N_t_ID, N_edgV_ID, N_edgH_ID, N_xc_ID, N_yc_ID, Quads_ID
  INTEGER,INTENT(IN):: RT_Its_ID, MGQD_Its_ID, GQD_Its_ID, Norm_Types_ID, MGQD_ResTypes_ID, Boundaries_ID
  INTEGER,INTENT(OUT):: Temp_ID, E_avg_ID, E_edgV_ID, E_edgH_ID, MGQD_E_avg_ID, MGQD_E_edgV_ID, MGQD_E_edgH_ID, HO_E_avg_ID
  INTEGER,INTENT(OUT):: HO_E_edgV_ID, HO_E_edgH_ID, Fx_edgV_ID, Fy_edgH_ID, MGQD_Fx_edgV_ID, MGQD_Fy_edgH_ID, HO_Fx_edgV_ID
  INTEGER,INTENT(OUT):: HO_Fy_edgH_ID, Eg_avg_ID, Eg_edgV_ID, Eg_edgH_ID, HO_Eg_avg_ID, HO_Eg_edgV_ID, HO_Eg_edgH_ID
  INTEGER,INTENT(OUT):: Fxg_edgV_ID, Fyg_edgH_ID, HO_Fxg_edgV_ID, HO_Fyg_edgH_ID, I_avg_ID, I_edgV_ID, I_edgH_ID
  INTEGER,INTENT(OUT):: KapE_Bar_ID, KapB_ID, KapE_ID, KapR_ID, Bg_ID
  INTEGER,INTENT(OUT):: RT_Residual_ID, MGQD_Residual_ID, MGQD_BC_Residual_ID
  INTEGER,INTENT(OUT):: Del_T_ID, Del_E_avg_ID, Del_E_edgV_ID, Del_E_edgH_ID, Del_Fx_edgV_ID, Del_Fy_edgH_ID
  INTEGER,INTENT(OUT):: RT_ItCount_ID, MGQD_ItCount_ID, GQD_ItCount_ID
  INTEGER,INTENT(OUT):: RT_Tnorm_ID, RT_Enorm_ID, MGQD_Tnorm_ID, MGQD_Enorm_ID
  INTEGER,INTENT(OUT):: RT_Trho_ID, RT_Erho_ID, MGQD_Trho_ID, MGQD_Erho_ID
  INTEGER,INTENT(OUT):: Cg_L_ID, Cg_B_ID, Cg_R_ID, Cg_T_ID, Eg_in_L_ID, Eg_in_B_ID, Eg_in_R_ID, Eg_in_T_ID
  INTEGER,INTENT(OUT):: Fg_in_L_ID, Fg_in_B_ID, Fg_in_R_ID, Fg_in_T_ID
  INTEGER,INTENT(OUT):: Cb_L_ID, Cb_B_ID, Cb_R_ID, Cb_T_ID, E_in_L_ID, E_in_B_ID, E_in_R_ID, E_in_T_ID
  INTEGER,INTENT(OUT):: F_in_L_ID, F_in_B_ID, F_in_R_ID, F_in_T_ID
  INTEGER,INTENT(OUT):: fg_avg_xx_ID, fg_avg_yy_ID, fg_edgV_xx_ID, fg_edgV_xy_ID, fg_edgH_yy_ID, fg_edgH_xy_ID
  INTEGER,INTENT(OUT):: DC_xx_ID, DL_xx_ID, DR_xx_ID, DC_yy_ID, DB_yy_ID, DT_yy_ID, DL_xy_ID, DB_xy_ID, DR_xy_ID, DT_xy_ID
  INTEGER,INTENT(OUT):: G_old_ID, Pold_L_ID, Pold_B_ID, Pold_R_ID, Pold_T_ID
  INTEGER,INTENT(OUT):: Gold_hat_ID, Rhat_old_ID, PL_ID, PB_ID, PR_ID, PT_ID
  INTEGER,INTENT(OUT):: dr_T_ID, dr_B_ID, dr_ML_ID, dr_MB_ID, dr_MR_ID, dr_MT_ID

  INTEGER:: Status

  !===========================================================================!
  !                                                                           !
  !     ENTERING FILE INTO OF DEFINE MODE                                     !
  !                                                                           !
  !===========================================================================!
  CALL NF_REDEF(outID)

  !===========================================================================!
  !                                                                           !
  !     DEFINING ITERATION COUNTS & CONVERGENCE NORMS                         !
  !                                                                           !
  !===========================================================================!
  CALL NF_DEF_VAR(RT_ItCount_ID,outID,(/N_t_ID/),'RT_Its-count','Int')
  CALL NF_DEF_VAR(MGQD_ItCount_ID,outID,(/RT_Its_ID,N_t_ID/),'MGQD_Its-count','Int')
  CALL NF_DEF_VAR(GQD_ItCount_ID,outID,(/MGQD_Its_ID,RT_Its_ID,N_t_ID/),'GQD_Its-count','Int')

  CALL NF_DEF_VAR(RT_Tnorm_ID,outID,(/RT_Its_ID,N_t_ID/),'RT_Tnorm','Real')
  CALL NF_DEF_VAR(RT_Enorm_ID,outID,(/RT_Its_ID,N_t_ID/),'RT_Enorm','Real')
  CALL NF_DEF_VAR(MGQD_Tnorm_ID,outID,(/MGQD_Its_ID,RT_Its_ID,N_t_ID/),'MGQD_Tnorm','Real')
  CALL NF_DEF_VAR(MGQD_Enorm_ID,outID,(/MGQD_Its_ID,RT_Its_ID,N_t_ID/),'MGQD_Enorm','Real')

  CALL NF_DEF_VAR(RT_Trho_ID,outID,(/RT_Its_ID,N_t_ID/),'RT_Trho','Real')
  CALL NF_DEF_VAR(RT_Erho_ID,outID,(/RT_Its_ID,N_t_ID/),'RT_Erho','Real')
  CALL NF_DEF_VAR(MGQD_Trho_ID,outID,(/MGQD_Its_ID,RT_Its_ID,N_t_ID/),'MGQD_Trho','Real')
  CALL NF_DEF_VAR(MGQD_Erho_ID,outID,(/MGQD_Its_ID,RT_Its_ID,N_t_ID/),'MGQD_Erho','Real')

  !===========================================================================!
  !                                                                           !
  !     DEFINING VARIABLES & UNITS                                            !
  !                                                                           !
  !===========================================================================!
  !--------------------------------------------------!
  !            Boundary Coefficients                 !
  !--------------------------------------------------!
  CALL NF_DEF_VAR(Cb_L_ID,outID,(/N_y_ID,N_t_ID/),'Cb_L','Double')
  CALL NF_DEF_VAR(Cb_B_ID,outID,(/N_x_ID,N_t_ID/),'Cb_B','Double')
  CALL NF_DEF_VAR(Cb_R_ID,outID,(/N_y_ID,N_t_ID/),'Cb_R','Double')
  CALL NF_DEF_VAR(Cb_T_ID,outID,(/N_x_ID,N_t_ID/),'Cb_T','Double')

  CALL NF_DEF_VAR(E_in_L_ID,outID,(/N_y_ID,N_t_ID/),'E_in_L','Double')
  CALL NF_DEF_VAR(E_in_B_ID,outID,(/N_x_ID,N_t_ID/),'E_in_B','Double')
  CALL NF_DEF_VAR(E_in_R_ID,outID,(/N_y_ID,N_t_ID/),'E_in_R','Double')
  CALL NF_DEF_VAR(E_in_T_ID,outID,(/N_x_ID,N_t_ID/),'E_in_T','Double')

  CALL NF_DEF_VAR(F_in_L_ID,outID,(/N_y_ID,N_t_ID/),'F_in_L','Double')
  CALL NF_DEF_VAR(F_in_B_ID,outID,(/N_x_ID,N_t_ID/),'F_in_B','Double')
  CALL NF_DEF_VAR(F_in_R_ID,outID,(/N_y_ID,N_t_ID/),'F_in_R','Double')
  CALL NF_DEF_VAR(F_in_T_ID,outID,(/N_x_ID,N_t_ID/),'F_in_T','Double')

  CALL NF_DEF_VAR(Cg_L_ID,outID,(/N_y_ID,N_g_ID,N_t_ID/),'Cg_L','Double')
  CALL NF_DEF_VAR(Cg_B_ID,outID,(/N_x_ID,N_g_ID,N_t_ID/),'Cg_B','Double')
  CALL NF_DEF_VAR(Cg_R_ID,outID,(/N_y_ID,N_g_ID,N_t_ID/),'Cg_R','Double')
  CALL NF_DEF_VAR(Cg_T_ID,outID,(/N_x_ID,N_g_ID,N_t_ID/),'Cg_T','Double')

  CALL NF_DEF_VAR(Eg_in_L_ID,outID,(/N_y_ID,N_g_ID,N_t_ID/),'Eg_in_L','Double')
  CALL NF_DEF_VAR(Eg_in_B_ID,outID,(/N_x_ID,N_g_ID,N_t_ID/),'Eg_in_B','Double')
  CALL NF_DEF_VAR(Eg_in_R_ID,outID,(/N_y_ID,N_g_ID,N_t_ID/),'Eg_in_R','Double')
  CALL NF_DEF_VAR(Eg_in_T_ID,outID,(/N_x_ID,N_g_ID,N_t_ID/),'Eg_in_T','Double')

  CALL NF_DEF_VAR(Fg_in_L_ID,outID,(/N_y_ID,N_g_ID,N_t_ID/),'Fg_in_L','Double')
  CALL NF_DEF_VAR(Fg_in_B_ID,outID,(/N_x_ID,N_g_ID,N_t_ID/),'Fg_in_B','Double')
  CALL NF_DEF_VAR(Fg_in_R_ID,outID,(/N_y_ID,N_g_ID,N_t_ID/),'Fg_in_R','Double')
  CALL NF_DEF_VAR(Fg_in_T_ID,outID,(/N_x_ID,N_g_ID,N_t_ID/),'Fg_in_T','Double')

  !--------------------------------------------------!
  !                 Temperature                      !
  !--------------------------------------------------!
  CALL NF_DEF_VAR(Temp_ID,outID,(/N_x_ID,N_y_ID,N_t_ID/),'Temperature','Double')
  CALL NF_DEF_UNIT(outID,Temp_ID,'ev/cm^3')

  !===========================================================================!
  !--------------------------------------------------!
  !       Total Radiation Energy Densities           !
  !--------------------------------------------------!
  !----------------(cell-averaged)-------------------!
  CALL NF_DEF_VAR(E_avg_ID,outID,(/N_x_ID,N_y_ID,N_t_ID/),'E_avg_Grey','Double') !from GQD eqs
  CALL NF_DEF_UNIT(outID,E_avg_ID,'ev/cm^3')

  CALL NF_DEF_VAR(MGQD_E_avg_ID,outID,(/N_x_ID,N_y_ID,N_t_ID/),'E_avg_MGQD','Double') !from MGQD eqs
  CALL NF_DEF_UNIT(outID,MGQD_E_avg_ID,'ev/cm^3')

  CALL NF_DEF_VAR(HO_E_avg_ID,outID,(/N_x_ID,N_y_ID,N_t_ID/),'E_avg_HO','Double') !from High-Order eqs
  CALL NF_DEF_UNIT(outID,HO_E_avg_ID,'ev/cm^3')

  !--------('vertical' cell edges, x-const)----------!
  CALL NF_DEF_VAR(E_edgV_ID,outID,(/N_edgV_ID,N_y_ID,N_t_ID/),'E_edgV_Grey','Double') !from GQD eqs
  CALL NF_DEF_UNIT(outID,E_edgV_ID,'ev/cm^3')

  CALL NF_DEF_VAR(MGQD_E_edgV_ID,outID,(/N_edgV_ID,N_y_ID,N_t_ID/),'E_edgV_MGQD','Double') !from MGQD eqs
  CALL NF_DEF_UNIT(outID,MGQD_E_edgV_ID,'ev/cm^3')

  CALL NF_DEF_VAR(HO_E_edgV_ID,outID,(/N_edgV_ID,N_y_ID,N_t_ID/),'E_edgV_HO','Double') !from High-Order eqs
  CALL NF_DEF_UNIT(outID,HO_E_edgV_ID,'ev/cm^3')

  !-------('horizontal' cell edges, x-const)---------!
  CALL NF_DEF_VAR(E_edgH_ID,outID,(/N_x_ID,N_edgH_ID,N_t_ID/),'E_edgH_Grey','Double') !from GQD eqs
  CALL NF_DEF_UNIT(outID,E_edgH_ID,'ev/cm^3')

  CALL NF_DEF_VAR(MGQD_E_edgH_ID,outID,(/N_x_ID,N_edgH_ID,N_t_ID/),'E_edgH_MGQD','Double') !from MGQD eqs
  CALL NF_DEF_UNIT(outID,MGQD_E_edgH_ID,'ev/cm^3')

  CALL NF_DEF_VAR(HO_E_edgH_ID,outID,(/N_x_ID,N_edgH_ID,N_t_ID/),'E_edgH_HO','Double') !from High-Order eqs
  CALL NF_DEF_UNIT(outID,HO_E_edgH_ID,'ev/cm^3')

  !===========================================================================!
  !--------------------------------------------------!
  !            Total Radiation Fluxes                !
  !--------------------------------------------------!
  !--------('vertical' cell edges, x-const)----------!
  CALL NF_DEF_VAR(Fx_edgV_ID,outID,(/N_edgV_ID,N_y_ID,N_t_ID/),'Fx_edgV_Grey','Double') !from GQD eqs
  CALL NF_DEF_UNIT(outID,Fx_edgV_ID,'ev/cm^3')

  CALL NF_DEF_VAR(MGQD_Fx_edgV_ID,outID,(/N_edgV_ID,N_y_ID,N_t_ID/),'Fx_edgV_MGQD','Double') !from MGQD eqs
  CALL NF_DEF_UNIT(outID,MGQD_Fx_edgV_ID,'ev/cm^3')

  CALL NF_DEF_VAR(HO_Fx_edgV_ID,outID,(/N_edgV_ID,N_y_ID,N_t_ID/),'Fx_edgV_HO','Double') !from High-Order eqs
  CALL NF_DEF_UNIT(outID,HO_Fx_edgV_ID,'ev/cm^3')

  !-------('horizontal' cell edges, x-const)---------!
  CALL NF_DEF_VAR(Fy_edgH_ID,outID,(/N_x_ID,N_edgH_ID,N_t_ID/),'Grey_Fy_edgH_Grey','Double') !from GQD eqs
  CALL NF_DEF_UNIT(outID,Fy_edgH_ID,'ev/cm^3')

  CALL NF_DEF_VAR(MGQD_Fy_edgH_ID,outID,(/N_x_ID,N_edgH_ID,N_t_ID/),'Fy_edgH_MGQD','Double') !from MGQD eqs
  CALL NF_DEF_UNIT(outID,MGQD_Fy_edgH_ID,'ev/cm^3')

  CALL NF_DEF_VAR(HO_Fy_edgH_ID,outID,(/N_x_ID,N_edgH_ID,N_t_ID/),'Fy_edgH_HO','Double') !from High-Order eqs
  CALL NF_DEF_UNIT(outID,HO_Fy_edgH_ID,'ev/cm^3')

  !===========================================================================!
  !--------------------------------------------------!
  !     Multigroup Radiation Energy Densities        !
  !--------------------------------------------------!
  !----------------(cell-averaged)-------------------!
  CALL NF_DEF_VAR(Eg_avg_ID,outID,(/N_x_ID,N_y_ID,N_g_ID,N_t_ID/),'Eg_avg_MGQD','Double') !from MGQD eqs
  CALL NF_DEF_UNIT(outID,Eg_avg_ID,'ev/cm^3')

  CALL NF_DEF_VAR(HO_Eg_avg_ID,outID,(/N_x_ID,N_y_ID,N_g_ID,N_t_ID/),'Eg_avg_HO','Double') !from High-Order eqs
  CALL NF_DEF_UNIT(outID,HO_Eg_avg_ID,'ev/cm^3')

  !--------('vertical' cell edges, x-const)----------!
  CALL NF_DEF_VAR(Eg_edgV_ID,outID,(/N_edgV_ID,N_y_ID,N_g_ID,N_t_ID/),'Eg_edgV_MGQD','Double') !from MGQD eqs
  CALL NF_DEF_UNIT(outID,Eg_edgV_ID,'ev/cm^3')

  CALL NF_DEF_VAR(HO_Eg_edgV_ID,outID,(/N_edgV_ID,N_y_ID,N_g_ID,N_t_ID/),'Eg_edgV_HO','Double') !from High-Order eqs
  CALL NF_DEF_UNIT(outID,HO_Eg_edgV_ID,'ev/cm^3')

  !-------('horizontal' cell edges, x-const)---------!
  CALL NF_DEF_VAR(Eg_edgH_ID,outID,(/N_x_ID,N_edgH_ID,N_g_ID,N_t_ID/),'Eg_edgH_MGQD','Double') !from MGQD eqs
  CALL NF_DEF_UNIT(outID,Eg_edgH_ID,'ev/cm^3')

  CALL NF_DEF_VAR(HO_Eg_edgH_ID,outID,(/N_x_ID,N_edgH_ID,N_g_ID,N_t_ID/),'Eg_edgH_HO','Double') !from High-Order eqs
  CALL NF_DEF_UNIT(outID,HO_Eg_edgH_ID,'ev/cm^3')

  !===========================================================================!
  !--------------------------------------------------!
  !          Multigroup Radiation Fluxes             !
  !--------------------------------------------------!
  !--------('vertical' cell edges, x-const)----------!
  CALL NF_DEF_VAR(Fxg_edgV_ID,outID,(/N_edgV_ID,N_y_ID,N_g_ID,N_t_ID/),'Fxg_edgV_MGQD','Double') !from MGQD eqs
  CALL NF_DEF_UNIT(outID,Fxg_edgV_ID,'ev/cm^3')

  CALL NF_DEF_VAR(HO_Fxg_edgV_ID,outID,(/N_edgV_ID,N_y_ID,N_g_ID,N_t_ID/),'Fxg_edgV_HO','Double') !from High-Order eqs
  CALL NF_DEF_UNIT(outID,HO_Fxg_edgV_ID,'ev/cm^3')

  !-------('horizontal' cell edges, x-const)---------!
  CALL NF_DEF_VAR(Fyg_edgH_ID,outID,(/N_x_ID,N_edgH_ID,N_g_ID,N_t_ID/),'Fyg_edgH_MGQD','Double') !from MGQD eqs
  CALL NF_DEF_UNIT(outID,Fyg_edgH_ID,'ev/cm^3')

  CALL NF_DEF_VAR(HO_Fyg_edgH_ID,outID,(/N_x_ID,N_edgH_ID,N_g_ID,N_t_ID/),'Fyg_edgH_HO','Double') !from High-Order eqs
  CALL NF_DEF_UNIT(outID,HO_Fyg_edgH_ID,'ev/cm^3')

  !===========================================================================!
  !--------------------------------------------------!
  !             Radiation Intensities                !
  !--------------------------------------------------!
  !----------------(cell-averaged)-------------------!
  CALL NF_DEF_VAR(I_avg_ID,outID,(/N_x_ID,N_y_ID,N_m_ID,N_g_ID,N_t_ID/),'I_avg','Double')
  CALL NF_DEF_UNIT(outID,I_avg_ID,'erg/(Ster*cm^3)')

  !--------('vertical' cell edges, x-const)----------!
  CALL NF_DEF_VAR(I_edgV_ID,outID,(/N_edgV_ID,N_y_ID,N_m_ID,N_g_ID,N_t_ID/),'I_edgV','Double')
  CALL NF_DEF_UNIT(outID,I_edgV_ID,'erg/(Ster*cm^3)')

  !-------('horizontal' cell edges, x-const)---------!
  CALL NF_DEF_VAR(I_edgH_ID,outID,(/N_x_ID,N_edgH_ID,N_m_ID,N_g_ID,N_t_ID/),'I_edgH','Double')
  CALL NF_DEF_UNIT(outID,I_edgH_ID,'erg/(Ster*cm^3)')

  !===========================================================================!
  !--------------------------------------------------!
  !                   Opacities                      !
  !--------------------------------------------------!
  CALL NF_DEF_VAR(KapE_Bar_ID,outID,(/N_x_ID,N_y_ID,N_t_ID/),'KapE_Bar','Double')
  CALL NF_DEF_UNIT(outID,KapE_Bar_ID,'1/cm')

  CALL NF_DEF_VAR(KapB_ID,outID,(/N_x_ID,N_y_ID,N_g_ID,N_t_ID/),'KapB','Double')
  CALL NF_DEF_UNIT(outID,KapB_ID,'1/cm')

  CALL NF_DEF_VAR(KapE_ID,outID,(/N_x_ID,N_y_ID,N_g_ID,N_t_ID/),'KapE','Double')
  CALL NF_DEF_UNIT(outID,KapB_ID,'1/cm')

  CALL NF_DEF_VAR(KapR_ID,outID,(/N_x_ID,N_y_ID,N_g_ID,N_t_ID/),'KapR','Double')
  CALL NF_DEF_UNIT(outID,KapB_ID,'1/cm')

  CALL NF_DEF_VAR(Bg_ID,outID,(/N_x_ID,N_y_ID,N_g_ID,N_t_ID/),'Bg','Double')

  !===========================================================================!
  !--------------------------------------------------!
  !                   QD factors                     !
  !--------------------------------------------------!
  CALL NF_DEF_VAR(fg_avg_xx_ID,outID,(/N_x_ID,N_y_ID,N_g_ID,N_t_ID/),'fg_avg_xx','Double')
  CALL NF_DEF_VAR(fg_avg_yy_ID,outID,(/N_x_ID,N_y_ID,N_g_ID,N_t_ID/),'fg_avg_yy','Double')
  CALL NF_DEF_VAR(fg_edgV_xx_ID,outID,(/N_edgV_ID,N_y_ID,N_g_ID,N_t_ID/),'fg_edgV_xx','Double')
  CALL NF_DEF_VAR(fg_edgV_xy_ID,outID,(/N_edgV_ID,N_y_ID,N_g_ID,N_t_ID/),'fg_edgV_xy','Double')
  CALL NF_DEF_VAR(fg_edgH_yy_ID,outID,(/N_x_ID,N_edgH_ID,N_g_ID,N_t_ID/),'fg_edgH_yy','Double')
  CALL NF_DEF_VAR(fg_edgH_xy_ID,outID,(/N_x_ID,N_edgH_ID,N_g_ID,N_t_ID/),'fg_edgH_xy','Double')

  CALL NF_DEF_VAR(G_old_ID,outID,(/N_x_ID,N_y_ID,N_g_ID,N_t_ID/),'Gold_L','Double')
  CALL NF_DEF_VAR(Pold_L_ID,outID,(/N_x_ID,N_y_ID,N_g_ID,N_t_ID/),'Pold_L','Double')
  CALL NF_DEF_VAR(Pold_B_ID,outID,(/N_x_ID,N_y_ID,N_g_ID,N_t_ID/),'Pold_B','Double')
  CALL NF_DEF_VAR(Pold_R_ID,outID,(/N_x_ID,N_y_ID,N_g_ID,N_t_ID/),'Pold_R','Double')
  CALL NF_DEF_VAR(Pold_T_ID,outID,(/N_x_ID,N_y_ID,N_g_ID,N_t_ID/),'Pold_T','Double')

  CALL NF_DEF_VAR(DC_xx_ID,outID,(/N_x_ID,N_y_ID,N_t_ID/),'DC_xx','Double')
  CALL NF_DEF_VAR(DL_xx_ID,outID,(/N_x_ID,N_y_ID,N_t_ID/),'DL_xx','Double')
  CALL NF_DEF_VAR(DR_xx_ID,outID,(/N_x_ID,N_y_ID,N_t_ID/),'DR_xx','Double')
  CALL NF_DEF_VAR(DC_yy_ID,outID,(/N_x_ID,N_y_ID,N_t_ID/),'DC_yy','Double')
  CALL NF_DEF_VAR(DB_yy_ID,outID,(/N_x_ID,N_y_ID,N_t_ID/),'DB_yy','Double')
  CALL NF_DEF_VAR(DT_yy_ID,outID,(/N_x_ID,N_y_ID,N_t_ID/),'DT_yy','Double')
  CALL NF_DEF_VAR(DL_xy_ID,outID,(/N_x_ID,N_y_ID,N_t_ID/),'DL_xy','Double')
  CALL NF_DEF_VAR(DB_xy_ID,outID,(/N_x_ID,N_y_ID,N_t_ID/),'DB_xy','Double')
  CALL NF_DEF_VAR(DR_xy_ID,outID,(/N_x_ID,N_y_ID,N_t_ID/),'DR_xy','Double')
  CALL NF_DEF_VAR(DT_xy_ID,outID,(/N_x_ID,N_y_ID,N_t_ID/),'DT_xy','Double')

  CALL NF_DEF_VAR(Gold_hat_ID,outID,(/N_x_ID,N_y_ID,N_t_ID/),'Gold_hat','Double')
  CALL NF_DEF_VAR(Rhat_old_ID,outID,(/N_x_ID,N_y_ID,N_t_ID/),'Rhat_old','Double')
  CALL NF_DEF_VAR(PL_ID,outID,(/N_x_ID,N_y_ID,N_t_ID/),'PL','Double')
  CALL NF_DEF_VAR(PB_ID,outID,(/N_x_ID,N_y_ID,N_t_ID/),'PB','Double')
  CALL NF_DEF_VAR(PR_ID,outID,(/N_x_ID,N_y_ID,N_t_ID/),'PR','Double')
  CALL NF_DEF_VAR(PT_ID,outID,(/N_x_ID,N_y_ID,N_t_ID/),'PT','Double')

  !===========================================================================!
  !                                                                           !
  !     DEFINING RESIDUALS                                                    !
  !                                                                           !
  !===========================================================================!
  IF (Res_Calc) THEN
    CALL NF_DEF_VAR(RT_Residual_ID,outID,(/Quads_ID,N_g_ID,Norm_Types_ID,RT_Its_ID,N_t_ID/),'RT_Residual','Real')

    CALL NF_DEF_VAR(MGQD_Residual_ID,outID,&
     (/N_g_ID,MGQD_ResTypes_ID,Norm_Types_ID,MGQD_Its_ID,RT_Its_ID,N_t_ID/),'MGQD_Residual','Real')
    CALL NF_DEF_VAR(MGQD_BC_Residual_ID,outID,&
     (/N_g_ID,Boundaries_ID,MGQD_Its_ID,RT_Its_ID,N_t_ID/),'MGQD_BC_Residual','Real')

    CALL NF_DEF_VAR(Del_T_ID,outID,(/GQD_Its_ID,MGQD_Its_ID,RT_Its_ID,N_t_ID/),'Delta_T','Real')
    CALL NF_DEF_VAR(Del_E_avg_ID,outID,(/GQD_Its_ID,MGQD_Its_ID,RT_Its_ID,N_t_ID/),'Delta_E_avg','Real')
    CALL NF_DEF_VAR(Del_E_edgV_ID,outID,(/GQD_Its_ID,MGQD_Its_ID,RT_Its_ID,N_t_ID/),'Delta_E_edgV','Real')
    CALL NF_DEF_VAR(Del_E_edgH_ID,outID,(/GQD_Its_ID,MGQD_Its_ID,RT_Its_ID,N_t_ID/),'Delta_E_edgH','Real')
    CALL NF_DEF_VAR(Del_Fx_edgV_ID,outID,(/GQD_Its_ID,MGQD_Its_ID,RT_Its_ID,N_t_ID/),'Delta_Fx_edgV','Real')
    CALL NF_DEF_VAR(Del_Fy_edgH_ID,outID,(/GQD_Its_ID,MGQD_Its_ID,RT_Its_ID,N_t_ID/),'Delta_Fy_edgH','Real')

    CALL NF_DEF_VAR(dr_T_ID,outID,(/GQD_Its_ID,MGQD_Its_ID,RT_Its_ID,N_t_ID/),'dr_T','Real')
    CALL NF_DEF_VAR(dr_B_ID,outID,(/GQD_Its_ID,MGQD_Its_ID,RT_Its_ID,N_t_ID/),'dr_B','Real')
    CALL NF_DEF_VAR(dr_ML_ID,outID,(/GQD_Its_ID,MGQD_Its_ID,RT_Its_ID,N_t_ID/),'dr_ML','Real')
    CALL NF_DEF_VAR(dr_MR_ID,outID,(/GQD_Its_ID,MGQD_Its_ID,RT_Its_ID,N_t_ID/),'dr_MR','Real')
    CALL NF_DEF_VAR(dr_MB_ID,outID,(/GQD_Its_ID,MGQD_Its_ID,RT_Its_ID,N_t_ID/),'dr_MB','Real')
    CALL NF_DEF_VAR(dr_MT_ID,outID,(/GQD_Its_ID,MGQD_Its_ID,RT_Its_ID,N_t_ID/),'dr_MT','Real')
  END IF

  !===========================================================================!
  !                                                                           !
  !     TAKING FILE OUT OF DEFINE MODE                                        !
  !                                                                           !
  !===========================================================================!
  CALL NF_ENDDEF(outID)

END SUBROUTINE OUTFILE_VARDEFS

!==================================================================================================================================!
!
!==================================================================================================================================!

END MODULE OUTPUTS
