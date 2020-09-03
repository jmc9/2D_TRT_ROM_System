MODULE OUTPUTS

  USE NCDF_IO

  IMPLICIT NONE

CONTAINS

!==================================================================================================================================!
!
!==================================================================================================================================!
SUBROUTINE OUTFILE_INIT(outID,N_x_ID,N_y_ID,N_m_ID,N_g_ID,N_t_ID,N_edgV_ID,N_edgH_ID,N_xc_ID,N_yc_ID,Quads_ID,RT_Its_ID,&
  MGQD_Its_ID,GQD_Its_ID,Norm_Types_ID,MGQD_ResTypes_ID,Boundaries_ID,c_ID,h_ID,pi_ID,erg_ID,Comp_Unit_ID,cv_ID,&
  c,h,pi,erg,Comp_Unit,cv,chi,conv_ho,conv_lo,conv_gr,line_src,xlen,ylen,Delx,Dely,tlen,Delt,bcT_left,bcT_bottom,&
  bcT_right,bcT_top,Tini,N_x,N_y,N_m,N_g,N_t,database_gen,use_grey,maxit_RTE,maxit_MLOQD,maxit_GLOQD,conv_type,&
  threads,BC_type,outfile,run_type,kapE_dT,quad,enrgy_strc)

  INTEGER,INTENT(OUT):: outID
  INTEGER,INTENT(OUT):: N_x_ID, N_y_ID, N_m_ID, N_g_ID, N_t_ID, N_edgV_ID, N_edgH_ID, N_xc_ID, N_yc_ID, Quads_ID
  INTEGER,INTENT(OUT):: RT_Its_ID, MGQD_Its_ID, GQD_Its_ID, Norm_Types_ID, MGQD_ResTypes_ID, Boundaries_ID
  INTEGER,INTENT(OUT):: c_ID, h_ID, pi_ID, erg_ID, Comp_Unit_ID, cv_ID

  REAL*8,INTENT(IN):: c, h, pi, erg, Comp_Unit, cv
  REAL*8,INTENT(IN):: chi, conv_ho, conv_lo, conv_gr, line_src, xlen, ylen, Delx(:), Dely(:), tlen, Delt
  REAL*8,INTENT(IN):: bcT_left, bcT_bottom, bcT_right, bcT_top, Tini
  INTEGER,INTENT(IN):: N_x, N_y, N_m, N_g, N_t
  INTEGER,INTENT(IN):: database_gen, use_grey, maxit_RTE, maxit_MLOQD, maxit_GLOQD, conv_type, threads, BC_type(:)
  CHARACTER(*),INTENT(IN):: outfile, run_type, kapE_dT, quad, enrgy_strc

  INTEGER:: Status
  INTEGER:: run_type_ID, database_gen_ID, use_grey_ID, maxit_RTE_ID, maxit_MLOQD_ID, maxit_GLOQD_ID, kapE_dT_ID, chi_ID
  INTEGER:: conv_ho_ID, conv_lo_ID, conv_gr_ID, conv_type_ID, quad_ID, threads_ID, enrgy_strc_ID, line_src_ID, xlen_ID, ylen_ID
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
  CALL NF_DEF_DIM(RT_Its_ID,outID,'RT_Its (dim)',0)
  CALL NF_DEF_DIM(MGQD_Its_ID,outID,'MGQD_Its (dim)',0)
  CALL NF_DEF_DIM(GQD_Its_ID,outID,'GQD_Its (dim)',0)
  CALL NF_DEF_DIM(Norm_Types_ID,outID,'Norm_Types',5) !inf, 1, 2, l1, l2
  CALL NF_DEF_DIM(MGQD_ResTypes_ID,outID,'MGQD_ResTypes',5)
  CALL NF_DEF_DIM(Boundaries_ID,outID,'Boundaries',4)

  !===========================================================================!
  !                                                                           !
  !     DEFINING CONSTANTS                                                    !
  !                                                                           !
  !===========================================================================!
  CALL NF_DEF_VAR(c_ID,outID,Z,'c (lightspeed)','Double')
  CALL NF_DEF_UNIT(outID,c_ID,'cm/sh')

  CALL NF_DEF_VAR(h_ID,outID,Z,'h (Boltzmann const)','Double')
  CALL NF_DEF_UNIT(outID,h_ID,'erg*sh')

  CALL NF_DEF_VAR(cv_ID,outID,Z,'cv (specific heat)','Double')
  CALL NF_DEF_UNIT(outID,cv_ID,'erg/(ev**4 cm**2 sh)')

  CALL NF_DEF_VAR(erg_ID,outID,Z,'ev <-> erg conversion','Double')
  CALL NF_DEF_UNIT(outID,erg_ID,'eV/erg')

  CALL NF_DEF_VAR(pi_ID,outID,Z,'pi','Double')

  CALL NF_DEF_VAR(Comp_Unit_ID,outID,Z,'Computational Unit','Double')

  !===========================================================================!
  !                                                                           !
  !     DEFINING PROBLEM PARAMETERS & USER INPUTS                             !
  !                                                                           !
  !===========================================================================!
  CALL NF_DEF_VAR(run_type_ID,outID,Z,'run_type','String')
  CALL NF_DEF_VAR(database_gen_ID,outID,Z,'database_gen','Int')
  CALL NF_DEF_VAR(use_grey_ID,outID,Z,'use_grey','Int')
  CALL NF_DEF_VAR(maxit_RTE_ID,outID,Z,'maxit_RTE','Int')
  CALL NF_DEF_VAR(maxit_MLOQD_ID,outID,Z,'maxit_MLOQD','Int')
  CALL NF_DEF_VAR(maxit_GLOQD_ID,outID,Z,'maxit_GLOQD','Int')
  CALL NF_DEF_VAR(kapE_dT_ID,outID,Z,'kapE_dT','String')
  CALL NF_DEF_VAR(chi_ID,outID,Z,'chi','Double')
  CALL NF_DEF_VAR(conv_ho_ID,outID,Z,'conv_ho','Double')
  CALL NF_DEF_VAR(conv_lo_ID,outID,Z,'conv_lo','Double')
  CALL NF_DEF_VAR(conv_gr_ID,outID,Z,'conv_gr','Double')
  CALL NF_DEF_VAR(conv_type_ID,outID,Z,'conv_type','Int')
  CALL NF_DEF_VAR(quad_ID,outID,Z,'quad','String')
  CALL NF_DEF_VAR(threads_ID,outID,Z,'threads','Int')
  CALL NF_DEF_VAR(enrgy_strc_ID,outID,Z,'enrgy_strc','String')
  CALL NF_DEF_VAR(line_src_ID,outID,Z,'line_src','Double')

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
  CALL NF_DEF_VAR(BC_type_ID,outID,(/Boundaries_ID/),'BC_type','Int')
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

  Status = nf90_put_var(outID,run_type_ID,run_type)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_var(outID,database_gen_ID,database_gen)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_var(outID,use_grey_ID,use_grey)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_var(outID,maxit_RTE_ID,maxit_RTE)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_var(outID,maxit_MLOQD_ID,maxit_MLOQD)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_var(outID,maxit_GLOQD_ID,maxit_GLOQD)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_var(outID,kapE_dT_ID,kapE_dT)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_var(outID,chi_ID,chi)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_var(outID,conv_ho_ID,conv_ho)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_var(outID,conv_lo_ID,conv_lo)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_var(outID,conv_gr_ID,conv_gr)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_var(outID,conv_type_ID,conv_type)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_var(outID,quad_ID,quad)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_var(outID,threads_ID,threads)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_var(outID,enrgy_strc_ID,enrgy_strc)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_var(outID,line_src_ID,line_src)
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

  Status = nf90_put_var(outID,BC_type_ID,BC_type,(/1/),(/4/))
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
SUBROUTINE OUTFILE_VARDEFS(outID,N_x_ID,N_y_ID,N_m_ID,N_g_ID,N_t_ID,N_edgV_ID,N_edgH_ID,N_xc_ID,N_yc_ID,Quads_ID,&
  RT_Its_ID,MGQD_Its_ID,GQD_Its_ID,Norm_Types_ID,MGQD_ResTypes_ID,Boundaries_ID,Temp_ID,E_avg_ID,E_edgV_ID,E_edgH_ID,&
  MGQD_E_avg_ID,MGQD_E_edgV_ID,MGQD_E_edgH_ID,HO_E_avg_ID,HO_E_edgV_ID,HO_E_edgH_ID,Fx_edgV_ID,Fy_edgH_ID,&
  MGQD_Fx_edgV_ID,MGQD_Fy_edgH_ID,HO_Fx_edgV_ID,HO_Fy_edgH_ID,Eg_avg_ID,Eg_edgV_ID,Eg_edgH_ID,HO_Eg_avg_ID,&
  HO_Eg_edgV_ID,HO_Eg_edgH_ID,Fxg_edgV_ID,Fyg_edgH_ID,HO_Fxg_edgV_ID,HO_Fyg_edgH_ID,I_avg_ID,I_edgV_ID,I_edgH_ID,&
  RT_Residual_ID,MGQD_Residual_ID,MGQD_BC_Residual_ID,Del_T_ID,Del_E_avg_ID,Del_E_edgV_ID,Del_E_edgH_ID,&
  Del_Fx_edgV_ID,Del_Fy_edgH_ID,Del_T_Norms_ID,Del_E_avg_Norms_ID,Del_E_edgV_Norms_ID,Del_E_edgH_Norms_ID,&
  Del_Fx_edgV_Norms_ID,Del_Fy_edgH_Norms_ID,RT_ItCount_ID,MGQD_ItCount_ID,GQD_ItCount_ID)

  INTEGER,INTENT(IN):: outID
  INTEGER,INTENT(IN):: N_x_ID, N_y_ID, N_m_ID, N_g_ID, N_t_ID, N_edgV_ID, N_edgH_ID, N_xc_ID, N_yc_ID, Quads_ID
  INTEGER,INTENT(IN):: RT_Its_ID, MGQD_Its_ID, GQD_Its_ID, Norm_Types_ID, MGQD_ResTypes_ID, Boundaries_ID
  INTEGER,INTENT(OUT):: Temp_ID, E_avg_ID, E_edgV_ID, E_edgH_ID, MGQD_E_avg_ID, MGQD_E_edgV_ID, MGQD_E_edgH_ID, HO_E_avg_ID
  INTEGER,INTENT(OUT):: HO_E_edgV_ID, HO_E_edgH_ID, Fx_edgV_ID, Fy_edgH_ID, MGQD_Fx_edgV_ID, MGQD_Fy_edgH_ID, HO_Fx_edgV_ID
  INTEGER,INTENT(OUT):: HO_Fy_edgH_ID, Eg_avg_ID, Eg_edgV_ID, Eg_edgH_ID, HO_Eg_avg_ID, HO_Eg_edgV_ID, HO_Eg_edgH_ID
  INTEGER,INTENT(OUT):: Fxg_edgV_ID, Fyg_edgH_ID, HO_Fxg_edgV_ID, HO_Fyg_edgH_ID, I_avg_ID, I_edgV_ID, I_edgH_ID
  INTEGER,INTENT(OUT):: RT_Residual_ID, MGQD_Residual_ID, MGQD_BC_Residual_ID
  INTEGER,INTENT(OUT):: Del_T_ID, Del_E_avg_ID, Del_E_edgV_ID, Del_E_edgH_ID, Del_Fx_edgV_ID, Del_Fy_edgH_ID, Del_T_Norms_ID
  INTEGER,INTENT(OUT):: Del_E_avg_Norms_ID, Del_E_edgV_Norms_ID, Del_E_edgH_Norms_ID, Del_Fx_edgV_Norms_ID, Del_Fy_edgH_Norms_ID
  INTEGER,INTENT(OUT):: RT_ItCount_ID, MGQD_ItCount_ID, GQD_ItCount_ID

  INTEGER:: Status

  !===========================================================================!
  !                                                                           !
  !     ENTERING FILE INTO OF DEFINE MODE                                     !
  !                                                                           !
  !===========================================================================!
  CALL NF_REDEF(outID)

  !===========================================================================!
  !                                                                           !
  !     DEFINING ITERATION COUNTS                                             !
  !                                                                           !
  !===========================================================================!
  CALL NF_DEF_VAR(RT_ItCount_ID,outID,(/N_t_ID/),'RT_Its (count)','Double')
  CALL NF_DEF_VAR(MGQD_ItCount_ID,outID,(/RT_Its_ID,N_t_ID/),'MGQD_Its (count)','Double')
  CALL NF_DEF_VAR(GQD_ItCount_ID,outID,(/MGQD_Its_ID,RT_Its_ID,N_t_ID/),'GQD_Its (count)','Double')

  !===========================================================================!
  !                                                                           !
  !     DEFINING VARIABLES & UNITS                                            !
  !                                                                           !
  !===========================================================================!
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
  CALL NF_DEF_VAR(E_avg_ID,outID,(/N_x_ID,N_y_ID,N_t_ID/),'E_avg (Grey)','Double') !from GQD eqs
  CALL NF_DEF_UNIT(outID,E_avg_ID,'ev/cm^3')

  CALL NF_DEF_VAR(MGQD_E_avg_ID,outID,(/N_x_ID,N_y_ID,N_t_ID/),'E_avg (MGQD)','Double') !from MGQD eqs
  CALL NF_DEF_UNIT(outID,MGQD_E_avg_ID,'ev/cm^3')

  CALL NF_DEF_VAR(HO_E_avg_ID,outID,(/N_x_ID,N_y_ID,N_t_ID/),'E_avg (HO)','Double') !from High-Order eqs
  CALL NF_DEF_UNIT(outID,HO_E_avg_ID,'ev/cm^3')

  !--------('vertical' cell edges, x-const)----------!
  CALL NF_DEF_VAR(E_edgV_ID,outID,(/N_edgV_ID,N_y_ID,N_t_ID/),'E_edgV (Grey)','Double') !from GQD eqs
  CALL NF_DEF_UNIT(outID,E_edgV_ID,'ev/cm^3')

  CALL NF_DEF_VAR(MGQD_E_edgV_ID,outID,(/N_edgV_ID,N_y_ID,N_t_ID/),'E_edgV (MGQD)','Double') !from MGQD eqs
  CALL NF_DEF_UNIT(outID,MGQD_E_edgV_ID,'ev/cm^3')

  CALL NF_DEF_VAR(HO_E_edgV_ID,outID,(/N_edgV_ID,N_y_ID,N_t_ID/),'E_edgV (HO)','Double') !from High-Order eqs
  CALL NF_DEF_UNIT(outID,HO_E_edgV_ID,'ev/cm^3')

  !-------('horizontal' cell edges, x-const)---------!
  CALL NF_DEF_VAR(E_edgH_ID,outID,(/N_x_ID,N_edgH_ID,N_t_ID/),'E_edgH (Grey)','Double') !from GQD eqs
  CALL NF_DEF_UNIT(outID,E_edgH_ID,'ev/cm^3')

  CALL NF_DEF_VAR(MGQD_E_edgH_ID,outID,(/N_x_ID,N_edgH_ID,N_t_ID/),'E_edgH (MGQD)','Double') !from MGQD eqs
  CALL NF_DEF_UNIT(outID,MGQD_E_edgH_ID,'ev/cm^3')

  CALL NF_DEF_VAR(HO_E_edgH_ID,outID,(/N_x_ID,N_edgH_ID,N_t_ID/),'E_edgH (HO)','Double') !from High-Order eqs
  CALL NF_DEF_UNIT(outID,HO_E_edgH_ID,'ev/cm^3')

  !===========================================================================!
  !--------------------------------------------------!
  !            Total Radiation Fluxes                !
  !--------------------------------------------------!
  !--------('vertical' cell edges, x-const)----------!
  CALL NF_DEF_VAR(Fx_edgV_ID,outID,(/N_edgV_ID,N_y_ID,N_t_ID/),'Fx_edgV (Grey)','Double') !from GQD eqs
  CALL NF_DEF_UNIT(outID,Fx_edgV_ID,'ev/cm^3')

  CALL NF_DEF_VAR(MGQD_Fx_edgV_ID,outID,(/N_edgV_ID,N_y_ID,N_t_ID/),'Fx_edgV (MGQD)','Double') !from MGQD eqs
  CALL NF_DEF_UNIT(outID,MGQD_Fx_edgV_ID,'ev/cm^3')

  CALL NF_DEF_VAR(HO_Fx_edgV_ID,outID,(/N_edgV_ID,N_y_ID,N_t_ID/),'Fx_edgV (HO)','Double') !from High-Order eqs
  CALL NF_DEF_UNIT(outID,HO_Fx_edgV_ID,'ev/cm^3')

  !-------('horizontal' cell edges, x-const)---------!
  CALL NF_DEF_VAR(Fy_edgH_ID,outID,(/N_x_ID,N_edgH_ID,N_t_ID/),'Fy_edgH (Grey)','Double') !from GQD eqs
  CALL NF_DEF_UNIT(outID,Fy_edgH_ID,'ev/cm^3')

  CALL NF_DEF_VAR(MGQD_Fy_edgH_ID,outID,(/N_x_ID,N_edgH_ID,N_t_ID/),'Fy_edgH (MGQD)','Double') !from MGQD eqs
  CALL NF_DEF_UNIT(outID,MGQD_Fy_edgH_ID,'ev/cm^3')

  CALL NF_DEF_VAR(HO_Fy_edgH_ID,outID,(/N_x_ID,N_edgH_ID,N_t_ID/),'Fy_edgH (HO)','Double') !from High-Order eqs
  CALL NF_DEF_UNIT(outID,HO_Fy_edgH_ID,'ev/cm^3')

  !===========================================================================!
  !--------------------------------------------------!
  !     Multigroup Radiation Energy Densities        !
  !--------------------------------------------------!
  !----------------(cell-averaged)-------------------!
  CALL NF_DEF_VAR(Eg_avg_ID,outID,(/N_x_ID,N_y_ID,N_g_ID,N_t_ID/),'Eg_avg (MGQD)','Double') !from MGQD eqs
  CALL NF_DEF_UNIT(outID,Eg_avg_ID,'ev/cm^3')

  CALL NF_DEF_VAR(HO_Eg_avg_ID,outID,(/N_x_ID,N_y_ID,N_g_ID,N_t_ID/),'Eg_avg (HO)','Double') !from High-Order eqs
  CALL NF_DEF_UNIT(outID,HO_Eg_avg_ID,'ev/cm^3')

  !--------('vertical' cell edges, x-const)----------!
  CALL NF_DEF_VAR(Eg_edgV_ID,outID,(/N_edgV_ID,N_y_ID,N_g_ID,N_t_ID/),'Eg_edgV (MGQD)','Double') !from MGQD eqs
  CALL NF_DEF_UNIT(outID,Eg_edgV_ID,'ev/cm^3')

  CALL NF_DEF_VAR(HO_Eg_edgV_ID,outID,(/N_edgV_ID,N_y_ID,N_g_ID,N_t_ID/),'Eg_edgV (HO)','Double') !from High-Order eqs
  CALL NF_DEF_UNIT(outID,HO_Eg_edgV_ID,'ev/cm^3')

  !-------('horizontal' cell edges, x-const)---------!
  CALL NF_DEF_VAR(Eg_edgH_ID,outID,(/N_x_ID,N_edgH_ID,N_g_ID,N_t_ID/),'Eg_edgH (MGQD)','Double') !from MGQD eqs
  CALL NF_DEF_UNIT(outID,Eg_edgH_ID,'ev/cm^3')

  CALL NF_DEF_VAR(HO_Eg_edgH_ID,outID,(/N_x_ID,N_edgH_ID,N_g_ID,N_t_ID/),'Eg_edgH (HO)','Double') !from High-Order eqs
  CALL NF_DEF_UNIT(outID,HO_Eg_edgH_ID,'ev/cm^3')

  !===========================================================================!
  !--------------------------------------------------!
  !          Multigroup Radiation Fluxes             !
  !--------------------------------------------------!
  !--------('vertical' cell edges, x-const)----------!
  CALL NF_DEF_VAR(Fxg_edgV_ID,outID,(/N_edgV_ID,N_y_ID,N_g_ID,N_t_ID/),'Fxg_edgV (MGQD)','Double') !from MGQD eqs
  CALL NF_DEF_UNIT(outID,Fxg_edgV_ID,'ev/cm^3')

  CALL NF_DEF_VAR(HO_Fxg_edgV_ID,outID,(/N_edgV_ID,N_y_ID,N_g_ID,N_t_ID/),'Fxg_edgV (HO)','Double') !from High-Order eqs
  CALL NF_DEF_UNIT(outID,HO_Fxg_edgV_ID,'ev/cm^3')

  !-------('horizontal' cell edges, x-const)---------!
  CALL NF_DEF_VAR(Fyg_edgH_ID,outID,(/N_x_ID,N_edgH_ID,N_g_ID,N_t_ID/),'Fyg_edgH (MGQD)','Double') !from MGQD eqs
  CALL NF_DEF_UNIT(outID,Fyg_edgH_ID,'ev/cm^3')

  CALL NF_DEF_VAR(HO_Fyg_edgH_ID,outID,(/N_x_ID,N_edgH_ID,N_g_ID,N_t_ID/),'Fyg_edgH (HO)','Double') !from High-Order eqs
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
  !                                                                           !
  !     DEFINING RESIDUALS                                                    !
  !                                                                           !
  !===========================================================================!
  CALL NF_DEF_VAR(RT_Residual_ID,outID,(/Quads_ID,N_g_ID,Norm_Types_ID,RT_Its_ID,N_t_ID/),'RT_Residual','Real')

  CALL NF_DEF_VAR(MGQD_Residual_ID,outID,&
   (/N_g_ID,MGQD_ResTypes_ID,Norm_Types_ID,MGQD_Its_ID,RT_Its_ID,N_t_ID/),'MGQD_Residual','Real')
  CALL NF_DEF_VAR(MGQD_BC_Residual_ID,outID,&
   (/N_g_ID,Boundaries_ID,MGQD_Its_ID,RT_Its_ID,N_t_ID/),'MGQD_BC_Residual','Real')

  CALL NF_DEF_VAR(Del_T_ID,outID,(/N_x_ID,N_y_ID,GQD_Its_ID,MGQD_Its_ID,RT_Its_ID,N_t_ID/),'Delta T','Double')
  CALL NF_DEF_VAR(Del_E_avg_ID,outID,(/N_x_ID,N_y_ID,GQD_Its_ID,MGQD_Its_ID,RT_Its_ID,N_t_ID/),'Delta E_avg','Double')
  CALL NF_DEF_VAR(Del_E_edgV_ID,outID,(/N_edgV_ID,N_y_ID,GQD_Its_ID,MGQD_Its_ID,RT_Its_ID,N_t_ID/),'Delta E_edgV','Double')
  CALL NF_DEF_VAR(Del_E_edgH_ID,outID,(/N_x_ID,N_edgH_ID,GQD_Its_ID,MGQD_Its_ID,RT_Its_ID,N_t_ID/),'Delta E_edgH','Double')
  CALL NF_DEF_VAR(Del_Fx_edgV_ID,outID,(/N_edgV_ID,N_y_ID,GQD_Its_ID,MGQD_Its_ID,RT_Its_ID,N_t_ID/),'Delta Fx_edgV','Double')
  CALL NF_DEF_VAR(Del_Fy_edgH_ID,outID,(/N_x_ID,N_edgH_ID,GQD_Its_ID,MGQD_Its_ID,RT_Its_ID,N_t_ID/),'Delta Fy_edgH','Double')

  CALL NF_DEF_VAR(Del_T_Norms_ID,outID,&
   (/Norm_Types_ID,GQD_Its_ID,MGQD_Its_ID,RT_Its_ID,N_t_ID/),'Delta T Norms','Double')
  CALL NF_DEF_VAR(Del_E_avg_Norms_ID,outID,&
   (/Norm_Types_ID,GQD_Its_ID,MGQD_Its_ID,RT_Its_ID,N_t_ID/),'Delta E_avg Norms','Double')
  CALL NF_DEF_VAR(Del_E_edgV_Norms_ID,outID,&
   (/Norm_Types_ID,GQD_Its_ID,MGQD_Its_ID,RT_Its_ID,N_t_ID/),'Delta E_edgV Norms','Double')
  CALL NF_DEF_VAR(Del_E_edgH_Norms_ID,outID,&
   (/Norm_Types_ID,GQD_Its_ID,MGQD_Its_ID,RT_Its_ID,N_t_ID/),'Delta E_edgH Norms','Double')
  CALL NF_DEF_VAR(Del_Fx_edgV_Norms_ID,outID,&
   (/Norm_Types_ID,GQD_Its_ID,MGQD_Its_ID,RT_Its_ID,N_t_ID/),'Delta Fx_edgV Norms','Double')
  CALL NF_DEF_VAR(Del_Fy_edgH_Norms_ID,outID,&
   (/Norm_Types_ID,GQD_Its_ID,MGQD_Its_ID,RT_Its_ID,N_t_ID/),'Delta Fy_edgH Norms','Double')

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
SUBROUTINE OUTFILE_HEADER(f_unit,datevalues,database_gen,database_add,run_type,restart_infile,use_grey,chi,&
  conv_ho,conv_lo,conv_gr,comp_unit,line_src,maxit_RTE,maxit_MLOQD,maxit_GLOQD,conv_type,N_m,threads,&
  kapE_dT_flag,enrgy_strc,Nwt_upbnd,xlen,ylen,N_x,N_y,tlen,delt,bc_type,bcT_left,bcT_right,bcT_upper,&
  bcT_lower,Tini,out_times,out_time_steps,n_out_times,restart_outlen,TEMP_out,GREY_E_out,GREY_F_out,&
  GREY_kap_out,GREY_fsmall_out,MG_fsmall_out,res_history_out,outfile,restart_outfile,decomp_outfile,TEMP_outfile,&
  GREY_E_outfile,GREY_F_outfile,GREY_kap_outfile,GREY_fsmall_outfile,MG_fsmall_outfile,res_history_outfile,&
  HO_E_out,HO_E_outfile,delx,dely,N_g,nu_g,N_t,sig_r,ar,sig0,pi,c,h,cv,Omega_x,Omega_y,quad_weight,quadrature)

  IMPLICIT NONE

  !INPUT VARIABLES
  INTEGER,INTENT(IN):: f_unit, datevalues(:)
  REAL*8,INTENT(IN):: chi, conv_ho, conv_lo, conv_gr, line_src, Nwt_upbnd, comp_unit
  REAL*8,INTENT(IN):: xlen, ylen, tlen, delt, bcT_left, bcT_right, bcT_upper, bcT_lower, Tini
  REAL*8,INTENT(IN):: Delx(:), Dely(:)
  REAL*8,INTENT(IN):: sig_r, ar, sig0, pi, c, h, cv
  REAL*8,INTENT(IN):: out_times(:), nu_g(:), Omega_x(:), Omega_y(:), quad_weight(:)
  INTEGER,INTENT(IN):: database_gen, database_add, use_grey, bc_type(:)
  INTEGER,INTENT(IN):: maxit_RTE, maxit_MLOQD, maxit_GLOQD, conv_type, N_m, threads, N_g, N_t
  INTEGER,INTENT(IN):: N_x, N_y
  INTEGER,INTENT(IN):: n_out_times, restart_outlen, TEMP_out, GREY_E_out, GREY_F_out, HO_E_out
  INTEGER,INTENT(IN):: GREY_kap_out, GREY_fsmall_out, MG_fsmall_out, res_history_out
  INTEGER,INTENT(IN):: out_time_steps(:)
  CHARACTER(100),INTENT(IN):: outfile, restart_outfile, decomp_outfile, TEMP_outfile, GREY_E_outfile, HO_E_outfile
  CHARACTER(100),INTENT(IN):: GREY_F_outfile, GREY_kap_outfile, GREY_fsmall_outfile, MG_fsmall_outfile, res_history_outfile
  CHARACTER(100),INTENT(IN):: run_type, restart_infile
  CHARACTER(100),INTENT(IN):: kapE_dT_flag, enrgy_strc, quadrature

  !LOCAL VARIABLES
  INTEGER:: g, m

  !CODE VERSION AND RUN DATE DATA
  WRITE(f_unit,'(A)') '-----------------------------------------------------------'
  WRITE(f_unit,'(A)') 'JOSEPH COALE 2D DISSERTATION CODE'
  WRITE(f_unit,'(2(A,I2.2),A,I4.4)') 'RUN DATE = ',datevalues(2),'/',datevalues(3),'/',datevalues(1)
  WRITE(f_unit,'(3(A,I2.2))') 'RUN TIME = ',datevalues(5),':',datevalues(6),':',datevalues(7)

  !INPUT ECHO
  WRITE(f_unit,'(A)') '-----------------------------------------------------------'
  WRITE(f_unit,'(A)') 'INPUT ECHO:'
  WRITE(f_unit,*)
  WRITE(f_unit,'(A)') '[RUN_STATE]'
  WRITE(f_unit,'(2A)')   '   run_type       ',run_type
  WRITE(f_unit,'(A,I1)') '   database_gen   ',database_gen
  WRITE(f_unit,'(A,I1)') '   database_add   ',database_add
  WRITE(f_unit,'(A,I1)') '   use_grey       ',use_grey
  IF (restart_infile .NE. '') WRITE(f_unit,'(2A)') '   restart_infile ',restart_infile
  WRITE(f_unit,*)

  WRITE(f_unit,'(A)') '[SOLVER_OPTS]'
  WRITE(f_unit,'(A,I3.3)')  '   maxit_RTE     ',maxit_RTE
  WRITE(f_unit,'(A,I4.4)')  '   maxit_MLOQD   ',maxit_MLOQD
  WRITE(f_unit,'(A,I4.4)')  '   maxit_GLOQD   ',maxit_GLOQD
  WRITE(f_unit,'(2A)')      '   kapE_dT       ',kapE_dT_flag
  WRITE(f_unit,'(A,ES7.1)') '   chi           ',chi
  WRITE(f_unit,'(A,ES7.1)') '   conv_ho       ',conv_ho
  WRITE(f_unit,'(A,ES7.1)') '   conv_lo       ',conv_lo
  WRITE(f_unit,'(A,ES7.1)') '   conv_gr       ',conv_gr
  WRITE(f_unit,'(A,I1)')    '   conv_type     ',conv_type
  WRITE(f_unit,'(2A)')      '   quad          ',quadrature
  WRITE(f_unit,'(A,I2.2)')  '   threads       ',threads
  WRITE(f_unit,'(2A)')      '   enrgy_strc    ',enrgy_strc
  WRITE(f_unit,'(A,ES7.1)') '   comp_unit     ',comp_unit
  WRITE(f_unit,'(A,ES7.1)') '   line_src      ',line_src
  WRITE(f_unit,*)

  WRITE(f_unit,'(A)') '[PARAMETERS]'
  WRITE(f_unit,'(A,ES9.3,A)') '   xlen      ',xlen,' (cm)'
  WRITE(f_unit,'(A,ES9.3,A)') '   ylen      ',ylen,' (cm)'
  WRITE(f_unit,'(A,I5)')      '   N_x       ',N_x
  WRITE(f_unit,'(A,I5)')      '   N_y       ',N_y
  WRITE(f_unit,'(A,ES9.3,A)') '   tlen      ',tlen,' (sh)'
  WRITE(f_unit,'(A,ES9.3,A)') '   delt      ',delt,' (sh)'
  WRITE(f_unit,'(A,4I2)')     '   bc_type   ',bc_type
  WRITE(f_unit,'(A,ES9.3,A)') '   bcT_left  ',bcT_left,' (ev)'
  WRITE(f_unit,'(A,ES9.3,A)') '   bcT_upper ',bcT_upper,' (ev)'
  WRITE(f_unit,'(A,ES9.3,A)') '   bcT_right ',bcT_right,' (ev)'
  WRITE(f_unit,'(A,ES9.3,A)') '   bcT_lower ',bcT_lower,' (ev)'
  WRITE(f_unit,'(A,ES9.3,A)') '   Tini      ',Tini,' (ev)'
  WRITE(f_unit,*)

  WRITE(f_unit,'(A)') '[OUTPUT_OPTS]'
  WRITE(f_unit,'(A,I1)') '   TEMP_out              ',TEMP_out
  WRITE(f_unit,'(A,I1)') '   GREY_E_out            ',GREY_E_out
  WRITE(f_unit,'(A,I1)') '   HO_E_out              ',HO_E_out
  WRITE(f_unit,'(A,I1)') '   GREY_F_out            ',GREY_F_out
  WRITE(f_unit,'(A,I1)') '   GREY_kap_out          ',GREY_kap_out
  WRITE(f_unit,'(A,I1)') '   GREY_fsmall_out       ',GREY_fsmall_out
  WRITE(f_unit,'(A,I1)') '   MG_fsmall_out         ',MG_fsmall_out
  WRITE(f_unit,'(A,I1)') '   res_history_out       ',res_history_out
  IF (restart_outlen .GT. 0) THEN
     WRITE(f_unit,'(A,I1)') '   restart_outlen        ',restart_outlen
     WRITE(f_unit,'(2A)')   '   restart_outfile       ',restart_outfile
   END IF
  WRITE(f_unit,'(2A)')   '   outfile               ',outfile
  IF (database_gen .EQ. 1)   WRITE(f_unit,'(2A)') '   decomp_outfile        ',decomp_outfile
  IF (TEMP_out .EQ. 1)     WRITE(f_unit,'(2A)') '   TEMP_outfile          ',TEMP_outfile
  IF (GREY_E_out .EQ. 1)   WRITE(f_unit,'(2A)') '   GREY_E_outfile        ',GREY_E_outfile
  IF (HO_E_out .EQ. 1)   WRITE(f_unit,'(2A)')   '   HO_E_outfile          ',HO_E_outfile
  IF (GREY_F_out .EQ. 1)   WRITE(f_unit,'(2A)') '   GREY_F_outfile        ',GREY_F_outfile
  IF (GREY_kap_out .EQ. 1) WRITE(f_unit,'(2A)') '   GREY_kap_outfile      ',GREY_kap_outfile
  IF (GREY_fsmall_out .EQ. 1) WRITE(f_unit,'(2A)') '   GREY_fsmall_outfile   ',GREY_fsmall_outfile
  IF (MG_fsmall_out .EQ. 1)   WRITE(f_unit,'(2A)') '   MG_fsmall_outfile     ',MG_fsmall_outfile
  IF (res_history_out .EQ. 1)   WRITE(f_unit,'(2A)') '   res_history_outfile   ',res_history_outfile
  WRITE(f_unit,*)

  !DISCRETIZATION DATA (PROBLEM DOMAIN AND ENERGY GROUPS)
  WRITE(f_unit,'(A)') '-----------------------------------------------------------'
  WRITE(f_unit,'(A)') 'DISCRETIZATION:'
  WRITE(f_unit,*)
  WRITE(f_unit,'(A,I5.5)') 'time_steps ',N_t
  WRITE(f_unit,'(A,I5.5)') 'N_x ',N_x
  WRITE(f_unit,'(A,I5.5)') 'N_y ',N_y
  WRITE(f_unit,'(A,ES9.3,A)') 'delx    ',Delx(1),' (cm)'
  WRITE(f_unit,'(A,ES9.3,A)') 'dely    ',Dely(1),' (cm)'
  WRITE(f_unit,'(A,I5.5)') 'N_g  ',N_g
  WRITE(f_unit,*)
  WRITE(f_unit,'(A)') 'group     energy interval (ev)'
  DO g=1,N_g
    WRITE(f_unit,'(I2.2,A,ES12.5,A,ES12.5,A)') g,'  [',nu_g(g),' -> ',nu_g(g+1),' ]'
  END DO
  WRITE(f_unit,*)
  WRITE(f_unit,'(A)') 'quadrature     angle (Omega_x, Omega_y, weight)'
  DO m=1,N_m
    WRITE(f_unit,'(I3.3,A,ES12.5,A,ES12.5,A,ES12.5,A)') m,'  [',Omega_x(m),' , ',Omega_y(m),' , ',quad_weight(m),' ]'
  END DO
  WRITE(f_unit,*)
  WRITE(f_unit,'(A,ES12.5)') 'SUM(Omega_x(:)) ',SUM(Omega_x(:))
  WRITE(f_unit,'(A,ES12.5)') 'SUM(Omega_y(:)) ',SUM(Omega_y(:))
  WRITE(f_unit,'(A,ES12.5)') 'SUM(quad_weight(:)) ',SUM(quad_weight(:))

  !CONSTANT VALUES ECHO
  WRITE(f_unit,'(A)') '-----------------------------------------------------------'
  WRITE(f_unit,'(A)') 'CONSTANTS:'
  WRITE(f_unit,*)
  WRITE(f_unit,'(A,ES22.15)') 'pi =    ',pi
  WRITE(f_unit,'(A,ES22.15,A)') 'c =     ',c,' (cm/sh)'
  WRITE(f_unit,'(A,ES22.15,A)') 'h =     ',h,' (erg*sh)'
  WRITE(f_unit,'(A,ES22.15,A)') 'Sig_R = ',sig_R*comp_unit,' (erg/(ev**4 cm**2 sh))'
  WRITE(f_unit,'(A,ES22.15,A)') 'a_R =   ',aR*comp_unit,' (erg/(ev**4 cm**3))'
  WRITE(f_unit,'(A,ES22.15,A)') 'Sig_0 = ',sig0,' (ev**3/cm)'
  WRITE(f_unit,'(A,ES22.15,A)') 'Cv =    ',cv*comp_unit,' (erg/(ev cm**3))'
  WRITE(f_unit,'(A,ES22.15,A)') 'comp_unit = ',comp_unit,' (erg/unit)'

END SUBROUTINE OUTFILE_HEADER

!==================================================================================================================================!
!
!==================================================================================================================================!
SUBROUTINE ITERATIVE_OUTPUT(f_unit,Tstep,TR_its,MGQD_its,GQD_its,Time,c,TR_Tnorm,TR_Enorm,TR_Trho,TR_Erho,MGQD_Tnorm,&
  MGQD_Enorm,MGQD_Trho,MGQD_Erho,TR_Residual,MGQD_Residual)

  INTEGER,INTENT(IN):: f_unit, Tstep, TR_its, MGQD_its(:), GQD_its
  REAL*8,INTENT(IN):: Time, c
  REAL*8,INTENT(IN):: TR_Residual(:,:,:,:,:,:), MGQD_Residual(:,:,:,:,:,:,:)
  REAL*8,INTENT(IN):: TR_Tnorm(:), TR_Enorm(:), TR_Trho(:), TR_Erho(:)
  REAL*8,INTENT(IN):: MGQD_Tnorm(:,:), MGQD_Enorm(:,:), MGQD_Trho(:,:), MGQD_Erho(:,:)
  INTEGER:: ti, mi, N_m, N_g, g, m1, m2, m

  N_m = SIZE(TR_Residual,3)
  N_g = SIZE(TR_Residual,4)

  WRITE(f_unit,'(A)') '*========================='
  WRITE(f_unit,'(A,I4)')    'TIME STEP ',Tstep
  WRITE(f_unit,'(A,ES9.3)') 'TIME      ',Time
  WRITE(f_unit,'(A,ES9.3)') 'ct        ',Time*c
  WRITE(f_unit,'(A,I4)')    'TR_its    ',TR_its
  WRITE(f_unit,'(A,I4)')    'MGQD_its  ',SUM(MGQD_its)
  WRITE(f_unit,'(A,I4)')    'GQD_its   ',GQD_its

  WRITE(f_unit,'(A)') '*'
  WRITE(f_unit,'(A)') 'Transport Iteration history'
  WRITE(f_unit,'(A)') 'TRit  Enorm     Erho      Tnorm     Trho'
  DO ti=1,TR_its
    WRITE(f_unit,'(I4.4,4ES10.3)') ti, TR_Tnorm(ti), TR_Trho(ti), TR_Enorm(ti), TR_Erho(ti)
  END DO

  WRITE(f_unit,'(A)') '*'
  WRITE(f_unit,'(A)') 'Transport Residuals (Summary)'
  WRITE(f_unit,'(A)') 'TRit  quad1     quad2     quad3     quad4'
  DO ti=1,TR_its
    WRITE(f_unit,'(I4.4,4ES10.3)') ti, MAXVAL(ABS(TR_Residual(:,:,1:N_m/4,:,1,ti))),&
                                       MAXVAL(ABS(TR_Residual(:,:,N_m/4+1:N_m/2,:,1,ti))),&
                                       MAXVAL(ABS(TR_Residual(:,:,N_m/2+1:3*N_m/4,:,1,ti))),&
                                       MAXVAL(ABS(TR_Residual(:,:,3*N_m/4+1:N_m,:,1,ti)))
  END DO

  WRITE(f_unit,'(A)') '*'
  WRITE(f_unit,'(A)') 'Transport Residuals (In Depth)'
  WRITE(f_unit,'(A)',ADVANCE='NO') 'TRit  quad    g01'
  DO g=2,N_g
    WRITE(f_unit,'(A,I2.2)',ADVANCE='NO') '       g',g
  END DO
  WRITE(f_unit,*)
  DO ti=1,TR_its
    m1=1
    m2=N_m/4
    DO m=1,4
      WRITE(f_unit,'(I4.4,I4,A)',ADVANCE='NO') ti, m,'   '
      DO g=1,N_g
        WRITE(f_unit,'(ES10.3)',ADVANCE='NO') MAXVAL(ABS(TR_Residual(:,:,m1:m2,g,1,ti)))
      END DO
      WRITE(f_unit,*)
      m1=m1+N_m/4
      m2=m2+N_m/4
    END DO
  END DO

  WRITE(f_unit,'(A)') '*'
  WRITE(f_unit,'(A)') 'MGQD Iteration history'
  WRITE(f_unit,'(A)') 'TRit MGQDit   Enorm     Erho      Tnorm     Trho'
  DO ti=1,TR_its
    DO mi=1,MGQD_its(ti)
      WRITE(f_unit,'(I4.4,I6.4,A,4ES10.3)') ti, mi, '  ', MGQD_Tnorm(ti,mi), MGQD_Trho(ti,mi),&
        MGQD_Enorm(ti,mi), MGQD_Erho(ti,mi)
    END DO
  END DO

  WRITE(f_unit,'(A)') '*'
  WRITE(f_unit,'(A)') 'MGQD Residuals (In Depth)'
  WRITE(f_unit,'(A)',ADVANCE='NO') 'TRit MGQDit  EQ      g01'
  DO g=2,N_g
    WRITE(f_unit,'(A,I2.2)',ADVANCE='NO') '       g',g
  END DO
  WRITE(f_unit,*)
  DO ti=1,TR_its
    DO mi=1,MGQD_its(ti)

      WRITE(f_unit,'(I4.4,I4,A)',ADVANCE='NO') ti, mi,'     EB   '
      DO g=1,N_g
        WRITE(f_unit,'(ES10.3)',ADVANCE='NO') MAXVAL(ABS(MGQD_Residual(:,:,g,1,1,mi,ti)))
      END DO
      WRITE(f_unit,*)

      WRITE(f_unit,'(I4.4,I4,A)',ADVANCE='NO') ti, mi,'     MBxL '
      DO g=1,N_g
        WRITE(f_unit,'(ES10.3)',ADVANCE='NO') MAXVAL(ABS(MGQD_Residual(:,:,g,2,1,mi,ti)))
      END DO
      WRITE(f_unit,*)

      WRITE(f_unit,'(I4.4,I4,A)',ADVANCE='NO') ti, mi,'     MBxR '
      DO g=1,N_g
        WRITE(f_unit,'(ES10.3)',ADVANCE='NO') MAXVAL(ABS(MGQD_Residual(:,:,g,3,1,mi,ti)))
      END DO
      WRITE(f_unit,*)

      WRITE(f_unit,'(I4.4,I4,A)',ADVANCE='NO') ti, mi,'     MByB '
      DO g=1,N_g
        WRITE(f_unit,'(ES10.3)',ADVANCE='NO') MAXVAL(ABS(MGQD_Residual(:,:,g,4,1,mi,ti)))
      END DO
      WRITE(f_unit,*)

      WRITE(f_unit,'(I4.4,I4,A)',ADVANCE='NO') ti, mi,'     MByT '
      DO g=1,N_g
        WRITE(f_unit,'(ES10.3)',ADVANCE='NO') MAXVAL(ABS(MGQD_Residual(:,:,g,5,1,mi,ti)))
      END DO
      WRITE(f_unit,*)
    END DO
  END DO
  WRITE(f_unit,'(A)') '*========================='

END SUBROUTINE ITERATIVE_OUTPUT

!==================================================================================================================================!
!
!==================================================================================================================================!
SUBROUTINE SPECIFIC_OUTPUTS(Temp_Times,HO_E_avg_Times,GREY_E_avg_Times,&
  Times,f_unit,datevalues,database_gen,database_add,run_type,restart_infile,&
  use_grey,chi,conv_ho,conv_lo,conv_gr,comp_unit,line_src,maxit_RTE,maxit_MLOQD,maxit_GLOQD,conv_type,N_m,threads,&
  kapE_dT_flag,enrgy_strc,Nwt_upbnd,xlen,ylen,N_x,N_y,tlen,delt,bc_type,bcT_left,bcT_right,bcT_upper,bcT_lower,Tini,&
  out_times,out_time_steps,n_out_times,restart_outlen,TEMP_out,GREY_E_out,GREY_F_out,GREY_kap_out,GREY_fsmall_out,&
  MG_fsmall_out,res_history_out,outfile,restart_outfile,decomp_outfile,TEMP_outfile,GREY_E_outfile,GREY_F_outfile,&
  GREY_kap_outfile,GREY_fsmall_outfile,MG_fsmall_outfile,res_history_outfile,HO_E_out,HO_E_outfile,delx,dely,N_g,&
  nu_g,N_t,sig_r,ar,sig0,pi,c,h,cv,Omega_x,Omega_y,quad_weight,quadrature,MGQD_E_avg_Times)

  REAL*8,INTENT(IN):: Temp_Times(:,:,:), HO_E_avg_Times(:,:,:), GREY_E_avg_Times(:,:,:), MGQD_E_avg_Times(:,:,:)
  REAL*8,INTENT(IN):: Times(:)

  !INPUT VARIABLES
  INTEGER,INTENT(IN):: f_unit, datevalues(:)
  REAL*8,INTENT(IN):: chi, conv_ho, conv_lo, conv_gr, line_src, Nwt_upbnd, comp_unit
  REAL*8,INTENT(IN):: xlen, ylen, tlen, delt, bcT_left, bcT_right, bcT_upper, bcT_lower, Tini
  REAL*8,INTENT(IN):: Delx(:), Dely(:)
  REAL*8,INTENT(IN):: sig_r, ar, sig0, pi, c, h, cv
  REAL*8,INTENT(IN):: out_times(:), nu_g(:), Omega_x(:), Omega_y(:), quad_weight(:)
  INTEGER,INTENT(IN):: database_gen, database_add, use_grey, bc_type(:)
  INTEGER,INTENT(IN):: maxit_RTE, maxit_MLOQD, maxit_GLOQD, conv_type, N_m, threads, N_g, N_t
  INTEGER,INTENT(IN):: N_x, N_y
  INTEGER,INTENT(IN):: n_out_times, restart_outlen, TEMP_out, GREY_E_out, GREY_F_out, HO_E_out
  INTEGER,INTENT(IN):: GREY_kap_out, GREY_fsmall_out, MG_fsmall_out, res_history_out
  INTEGER,INTENT(IN):: out_time_steps(:)
  CHARACTER(100),INTENT(IN):: outfile, restart_outfile, decomp_outfile, TEMP_outfile, GREY_E_outfile, HO_E_outfile
  CHARACTER(100),INTENT(IN):: GREY_F_outfile, GREY_kap_outfile, GREY_fsmall_outfile, MG_fsmall_outfile, res_history_outfile
  CHARACTER(100),INTENT(IN):: run_type, restart_infile
  CHARACTER(100),INTENT(IN):: kapE_dT_flag, enrgy_strc, quadrature

  CHARACTER(100):: Data_Name, file

  IF (TEMP_out .EQ. 1) THEN
    Data_Name = 'Temperature'
    CALL SPCF_OUTPUT(Temp_Times,Times,TEMP_outfile,Data_Name,f_unit,datevalues,database_gen,database_add,run_type,restart_infile,&
      use_grey,chi,conv_ho,conv_lo,conv_gr,comp_unit,line_src,maxit_RTE,maxit_MLOQD,maxit_GLOQD,conv_type,N_m,threads,&
      kapE_dT_flag,enrgy_strc,Nwt_upbnd,xlen,ylen,N_x,N_y,tlen,delt,bcT_left,bcT_right,bcT_upper,bcT_lower,Tini,&
      out_times,out_time_steps,n_out_times,restart_outlen,TEMP_out,GREY_E_out,GREY_F_out,GREY_kap_out,GREY_fsmall_out,&
      MG_fsmall_out,res_history_out,outfile,restart_outfile,decomp_outfile,TEMP_outfile,GREY_E_outfile,GREY_F_outfile,&
      GREY_kap_outfile,GREY_fsmall_outfile,MG_fsmall_outfile,res_history_outfile,HO_E_out,HO_E_outfile,delx,dely,N_g,&
      nu_g,N_t,sig_r,ar,sig0,pi,c,h,cv,Omega_x,Omega_y,quad_weight,quadrature,bc_type)
  END IF

  IF (HO_E_out .EQ. 1) THEN
    Data_Name = 'High-Order Radiation Energy Density'
    CALL SPCF_OUTPUT(HO_E_avg_Times,Times,HO_E_outfile,Data_Name,f_unit,datevalues,database_gen,database_add,run_type,&
      restart_infile,use_grey,chi,conv_ho,conv_lo,conv_gr,comp_unit,line_src,maxit_RTE,maxit_MLOQD,maxit_GLOQD,conv_type,N_m,&
      threads,kapE_dT_flag,enrgy_strc,Nwt_upbnd,xlen,ylen,N_x,N_y,tlen,delt,bcT_left,bcT_right,bcT_upper,bcT_lower,Tini,&
      out_times,out_time_steps,n_out_times,restart_outlen,TEMP_out,GREY_E_out,GREY_F_out,GREY_kap_out,GREY_fsmall_out,&
      MG_fsmall_out,res_history_out,outfile,restart_outfile,decomp_outfile,TEMP_outfile,GREY_E_outfile,GREY_F_outfile,&
      GREY_kap_outfile,GREY_fsmall_outfile,MG_fsmall_outfile,res_history_outfile,HO_E_out,HO_E_outfile,delx,dely,N_g,&
      nu_g,N_t,sig_r,ar,sig0,pi,c,h,cv,Omega_x,Omega_y,quad_weight,quadrature,bc_type)
  END IF

  IF (HO_E_out .EQ. 1) THEN
    Data_Name = 'MGQD Radiation Energy Density'
    file = 'E_MGQD.out'
    CALL SPCF_OUTPUT(MGQD_E_avg_Times,Times,file,Data_Name,f_unit,datevalues,database_gen,database_add,run_type,&
      restart_infile,use_grey,chi,conv_ho,conv_lo,conv_gr,comp_unit,line_src,maxit_RTE,maxit_MLOQD,maxit_GLOQD,conv_type,N_m,&
      threads,kapE_dT_flag,enrgy_strc,Nwt_upbnd,xlen,ylen,N_x,N_y,tlen,delt,bcT_left,bcT_right,bcT_upper,bcT_lower,Tini,&
      out_times,out_time_steps,n_out_times,restart_outlen,TEMP_out,GREY_E_out,GREY_F_out,GREY_kap_out,GREY_fsmall_out,&
      MG_fsmall_out,res_history_out,outfile,restart_outfile,decomp_outfile,TEMP_outfile,GREY_E_outfile,GREY_F_outfile,&
      GREY_kap_outfile,GREY_fsmall_outfile,MG_fsmall_outfile,res_history_outfile,HO_E_out,HO_E_outfile,delx,dely,N_g,&
      nu_g,N_t,sig_r,ar,sig0,pi,c,h,cv,Omega_x,Omega_y,quad_weight,quadrature,bc_type)
  END IF

  IF (GREY_E_out .EQ. 1) THEN
    Data_Name = 'Grey Radiation Energy Density'
    CALL SPCF_OUTPUT(GREY_E_avg_Times,Times,GREY_E_outfile,Data_Name,f_unit,datevalues,database_gen,database_add,run_type,&
      restart_infile,use_grey,chi,conv_ho,conv_lo,conv_gr,comp_unit,line_src,maxit_RTE,maxit_MLOQD,maxit_GLOQD,conv_type,N_m,&
      threads,kapE_dT_flag,enrgy_strc,Nwt_upbnd,xlen,ylen,N_x,N_y,tlen,delt,bcT_left,bcT_right,bcT_upper,bcT_lower,Tini,&
      out_times,out_time_steps,n_out_times,restart_outlen,TEMP_out,GREY_E_out,GREY_F_out,GREY_kap_out,GREY_fsmall_out,&
      MG_fsmall_out,res_history_out,outfile,restart_outfile,decomp_outfile,TEMP_outfile,GREY_E_outfile,GREY_F_outfile,&
      GREY_kap_outfile,GREY_fsmall_outfile,MG_fsmall_outfile,res_history_outfile,HO_E_out,HO_E_outfile,delx,dely,N_g,&
      nu_g,N_t,sig_r,ar,sig0,pi,c,h,cv,Omega_x,Omega_y,quad_weight,quadrature,bc_type)
  END IF

END SUBROUTINE SPECIFIC_OUTPUTS

!============================================================================================================!
!
!============================================================================================================!
SUBROUTINE SPCF_OUTPUT(Data,Times,File,Data_Name,f_unit,datevalues,database_gen,database_add,run_type,restart_infile,&
  use_grey,chi,conv_ho,conv_lo,conv_gr,comp_unit,line_src,maxit_RTE,maxit_MLOQD,maxit_GLOQD,conv_type,N_m,threads,&
  kapE_dT_flag,enrgy_strc,Nwt_upbnd,xlen,ylen,N_x,N_y,tlen,delt,bcT_left,bcT_right,bcT_upper,bcT_lower,Tini,&
  out_times,out_time_steps,n_out_times,restart_outlen,TEMP_out,GREY_E_out,GREY_F_out,GREY_kap_out,GREY_fsmall_out,&
  MG_fsmall_out,res_history_out,outfile,restart_outfile,decomp_outfile,TEMP_outfile,GREY_E_outfile,GREY_F_outfile,&
  GREY_kap_outfile,GREY_fsmall_outfile,MG_fsmall_outfile,res_history_outfile,HO_E_out,HO_E_outfile,delx,dely,N_g,&
  nu_g,N_t,sig_r,ar,sig0,pi,c,h,cv,Omega_x,Omega_y,quad_weight,quadrature,bc_type)

  REAL*8,INTENT(IN):: Data(:,:,:)
  REAL*8,INTENT(IN):: Times(:)
  CHARACTER(100),INTENT(IN):: File, Data_Name

  !INPUT VARIABLES
  INTEGER,INTENT(IN):: f_unit, datevalues(:)
  REAL*8,INTENT(IN):: chi, conv_ho, conv_lo, conv_gr, line_src, Nwt_upbnd, comp_unit
  REAL*8,INTENT(IN):: xlen, ylen, tlen, delt, bcT_left, bcT_right, bcT_upper, bcT_lower, Tini
  REAL*8,INTENT(IN):: Delx(:), Dely(:)
  REAL*8,INTENT(IN):: sig_r, ar, sig0, pi, c, h, cv
  REAL*8,INTENT(IN):: out_times(:), nu_g(:), Omega_x(:), Omega_y(:), quad_weight(:)
  INTEGER,INTENT(IN):: database_gen, database_add, use_grey, bc_type(:)
  INTEGER,INTENT(IN):: maxit_RTE, maxit_MLOQD, maxit_GLOQD, conv_type, N_m, threads, N_g, N_t
  INTEGER,INTENT(IN):: N_x, N_y
  INTEGER,INTENT(IN):: n_out_times, restart_outlen, TEMP_out, GREY_E_out, GREY_F_out, HO_E_out
  INTEGER,INTENT(IN):: GREY_kap_out, GREY_fsmall_out, MG_fsmall_out, res_history_out
  INTEGER,INTENT(IN):: out_time_steps(:)
  CHARACTER(100),INTENT(IN):: outfile, restart_outfile, decomp_outfile, TEMP_outfile, GREY_E_outfile, HO_E_outfile
  CHARACTER(100),INTENT(IN):: GREY_F_outfile, GREY_kap_outfile, GREY_fsmall_outfile, MG_fsmall_outfile, res_history_outfile
  CHARACTER(100),INTENT(IN):: run_type, restart_infile
  CHARACTER(100),INTENT(IN):: kapE_dT_flag, enrgy_strc, quadrature

  INTEGER:: N_t2, t, i, j
  INTEGER:: err

  N_t2 = SIZE(Times,1)

  OPEN(UNIT=f_unit,FILE=File,STATUS='REPLACE',ACTION='WRITE',IOSTAT=err)
!   making sure file exists/opens, if not tells user
  IF (err .NE. 0) THEN
      WRITE(*,'(3A)') 'The file, ',File,', could not open properly.'
      STOP
  END IF

  CALL OUTFILE_HEADER(f_unit,datevalues,database_gen,database_add,run_type,restart_infile,use_grey,chi,&
    conv_ho,conv_lo,conv_gr,comp_unit,line_src,maxit_RTE,maxit_MLOQD,maxit_GLOQD,conv_type,N_m,threads,&
    kapE_dT_flag,enrgy_strc,Nwt_upbnd,xlen,ylen,N_x,N_y,tlen,delt,bc_type,bcT_left,bcT_right,bcT_upper,&
    bcT_lower,Tini,out_times,out_time_steps,n_out_times,restart_outlen,TEMP_out,GREY_E_out,GREY_F_out,&
    GREY_kap_out,GREY_fsmall_out,MG_fsmall_out,res_history_out,outfile,restart_outfile,decomp_outfile,TEMP_outfile,&
    GREY_E_outfile,GREY_F_outfile,GREY_kap_outfile,GREY_fsmall_outfile,MG_fsmall_outfile,res_history_outfile,&
    HO_E_out,HO_E_outfile,delx,dely,N_g,nu_g,N_t,sig_r,ar,sig0,pi,c,h,cv,Omega_x,Omega_y,quad_weight,quadrature)

  WRITE(f_unit,'(A)') '-----------------------------------------------------------'
  WRITE(f_unit,'(A)') 'OUT_TIMES:'
  WRITE(f_unit,*)
  WRITE(f_unit,'(A)') 'time_step  |  time (ns)  |  c*t'
  DO t=1,N_t2
    WRITE(f_unit,'(I5.5,A,ES12.5,A,ES12.5)') t,' , ',Times(t),' , ',c*Times(t)
  END DO

  WRITE(f_unit,'(A)') '-----------------------------------------------------------'
  WRITE(f_unit,'(A)') trim(Data_Name)
  WRITE(f_unit,*)

  WRITE(f_unit,'(A)',ADVANCE='NO') ' i   j        x             y               '
  DO t=1,N_t2
    WRITE(f_unit,'(I5.5,A)',ADVANCE='NO') t,'                  '
  END DO
  WRITE(f_unit,*)
  DO j=1,N_y
    DO i=1,N_x
      WRITE(f_unit,'(I3.3,A,I3.3,2ES14.6)',ADVANCE='NO') i,' ', j, (DBLE(i)-.5d0)*Delx(i), (DBLE(j)-.5d0)*Dely(j)
      DO t=1,N_t2
        WRITE(f_unit,'(ES23.15)',ADVANCE='NO') Data(i,j,t)
      END DO
      WRITE(f_unit,*)
    END DO
  END DO

  CLOSE ( f_unit, STATUS='KEEP')

END SUBROUTINE SPCF_OUTPUT

!==================================================================================================================================!
!
!==================================================================================================================================!

END MODULE OUTPUTS
