MODULE OUTPUTS

  USE NCDF_IO
  USE netcdf

  IMPLICIT NONE

CONTAINS

!==================================================================================================================================!
!
!==================================================================================================================================!
SUBROUTINE OUTFILE_INIT(outID,N_x_ID,N_y_ID,N_m_ID,N_g_ID,N_t_ID,N_edgV_ID,N_edgH_ID,N_xc_ID,N_yc_ID,Quads_ID,RT_Its_ID,&
  MGQD_Its_ID,GQD_Its_ID,Norm_Types_ID,MGQD_ResTypes_ID,Boundaries_ID,c_ID,h_ID,pi_ID,erg_ID,Comp_Unit_ID,cv_ID,&
  c,h,pi,erg,Comp_Unit,cv,chi,conv_ho,conv_lo,conv_gr1,conv_gr2,line_src,xlen,ylen,Delx,Dely,tlen,Delt,bcT_left,&
  bcT_bottom,bcT_right,bcT_top,Tini,E_Bound_Low,T_Bound_Low,N_x,N_y,N_m,N_g,N_t,use_grey,maxit_RTE,&
  maxit_MLOQD,maxit_GLOQD,conv_type,threads,BC_type,outfile,run_type,kapE_dT_flag,quadrature,enrgy_strc,Theta,&
  Use_Line_Search,Use_Safety_Search,I_out,HO_Eg_out,HO_Fg_out,HO_E_out,HO_F_out,Eg_out,Fg_out,MGQD_E_out,MGQD_F_out,&
  QDfg_out,E_out,F_out,D_out,old_parms_out,its_out,conv_out,kap_out,Src_out,POD_Type,Direc_Diff,&
  xpts_avg,xpts_edgV,ypts_avg,ypts_edgH,tpts,N_dsets,DMD_dsets,DMD_Type,dset_times,N_dsets_ID,POD_dsets)

  INTEGER,INTENT(OUT):: outID
  INTEGER,INTENT(OUT):: N_x_ID, N_y_ID, N_m_ID, N_g_ID, N_t_ID, N_edgV_ID, N_edgH_ID, N_xc_ID, N_yc_ID, Quads_ID
  INTEGER,INTENT(OUT):: RT_Its_ID, MGQD_Its_ID, GQD_Its_ID, Norm_Types_ID, MGQD_ResTypes_ID, Boundaries_ID
  INTEGER,INTENT(OUT):: c_ID, h_ID, pi_ID, erg_ID, Comp_Unit_ID, cv_ID, N_dsets_ID

  REAL*8,INTENT(IN):: c, h, pi, erg, Comp_Unit, cv
  REAL*8,INTENT(IN):: chi, conv_ho, conv_lo, conv_gr1, conv_gr2, line_src, xlen, ylen, Delx(:), Dely(:), tlen, Delt, Theta
  REAL*8,INTENT(IN):: xpts_avg(:), xpts_edgV(:), ypts_avg(:), ypts_edgH(:), tpts(:), dset_times(:)
  REAL*8,INTENT(IN):: bcT_left, bcT_bottom, bcT_right, bcT_top, Tini, E_Bound_Low, T_Bound_Low
  INTEGER,INTENT(IN):: N_x, N_y, N_m, N_g, N_t
  INTEGER,INTENT(IN):: use_grey, maxit_RTE, maxit_MLOQD, maxit_GLOQD, conv_type, threads, BC_type(:)
  INTEGER,INTENT(IN):: I_out, HO_Eg_out, HO_Fg_out, HO_E_out, HO_F_out
  INTEGER,INTENT(IN):: Eg_out, Fg_out, MGQD_E_out, MGQD_F_out, QDfg_out, E_out, F_out, D_out
  INTEGER,INTENT(IN):: old_parms_out, its_out, conv_out, kap_out, Src_out, Direc_Diff, N_dsets
  CHARACTER(*),INTENT(IN):: outfile, run_type, quadrature, enrgy_strc, POD_Type, DMD_Type, DMD_dsets(:), POD_dsets(:)
  LOGICAL,INTENT(IN):: Use_Line_Search, Use_Safety_Search, kapE_dT_flag

  INTEGER:: Status
  INTEGER:: xlen_ID, ylen_ID
  INTEGER:: Delx_ID, Dely_ID, tlen_ID, Delt_ID, BC_type_ID, bcT_left_ID, bcT_bottom_ID, bcT_right_ID, bcT_top_ID, Tini_ID
  INTEGER:: xpts_avg_ID, xpts_edgV_ID, ypts_avg_ID, ypts_edgH_ID, tpts_ID
  INTEGER:: Z(0), i
  CHARACTER(50):: buf

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

  IF (run_type .EQ. 'mg_pod') THEN
    Status = nf90_put_att(outID,NF90_GLOBAL,'POD_Type',POD_Type)
    CALL HANDLE_ERR(Status)

    ! Status = nf90_put_att(outID,NF90_GLOBAL,'PODgsum',PODgsum)
    ! CALL HANDLE_ERR(Status)

    Status = nf90_put_att(outID,NF90_GLOBAL,'Direc_Diff',Direc_Diff)
    CALL HANDLE_ERR(Status)

    ! Status = nf90_put_att(outID,NF90_GLOBAL,'POD_err',POD_err)
    ! CALL HANDLE_ERR(Status)

    CALL NF_DEF_DIM(N_dsets_ID,outID,'N_dsets',N_dsets)

    Status = nf90_put_att(outID,NF90_GLOBAL,'dset_times',dset_times)
    CALL HANDLE_ERR(Status)

    DO i=1,N_dsets
      write(buf,'(A,I1)') 'POD_dset',i
      Status = nf90_put_att(outID,NF90_GLOBAL,buf,POD_dsets(i))
      CALL HANDLE_ERR(Status)
    END DO

  ELSE IF (run_type .EQ. 'mg_dmd') THEN
    Status = nf90_put_att(outID,NF90_GLOBAL,'DMD_Type',DMD_Type)
    CALL HANDLE_ERR(Status)

    CALL NF_DEF_DIM(N_dsets_ID,outID,'N_dsets',N_dsets)

    Status = nf90_put_att(outID,NF90_GLOBAL,'dset_times',dset_times)
    CALL HANDLE_ERR(Status)

    DO i=1,N_dsets
      write(buf,'(A,I1)') 'DMD_dset',i
      Status = nf90_put_att(outID,NF90_GLOBAL,buf,DMD_dsets(i))
      CALL HANDLE_ERR(Status)
    END DO

  END IF

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

  Status = nf90_put_att(outID,NF90_GLOBAL,'I_out',I_out)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_att(outID,NF90_GLOBAL,'D_out',D_out)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_att(outID,NF90_GLOBAL,'fg_out',QDfg_out)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_att(outID,NF90_GLOBAL,'old_parms_out',old_parms_out)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_att(outID,NF90_GLOBAL,'HO_Eg_out',HO_Eg_out)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_att(outID,NF90_GLOBAL,'HO_E_out',HO_E_out)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_att(outID,NF90_GLOBAL,'HO_Fg_out',HO_Fg_out)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_att(outID,NF90_GLOBAL,'HO_F_out',HO_F_out)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_att(outID,NF90_GLOBAL,'Eg_out',Eg_out)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_att(outID,NF90_GLOBAL,'MGQD_E_out',MGQD_E_out)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_att(outID,NF90_GLOBAL,'Fg_out',Fg_out)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_att(outID,NF90_GLOBAL,'MGQD_F_out',MGQD_F_out)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_att(outID,NF90_GLOBAL,'E_out',E_out)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_att(outID,NF90_GLOBAL,'F_out',F_out)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_att(outID,NF90_GLOBAL,'its_out',its_out)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_att(outID,NF90_GLOBAL,'conv_out',conv_out)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_att(outID,NF90_GLOBAL,'kap_out',kap_out)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_att(outID,NF90_GLOBAL,'Src_out',Src_out)
  CALL HANDLE_ERR(Status)

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

  CALL NF_DEF_VAR(xpts_avg_ID,outID,(/N_x_ID/),'xpts_avg','Double')
  Status = nf90_put_att(outID,xpts_avg_ID,'bnds',(/0d0,xlen/))
  CALL NF_DEF_VAR(xpts_edgV_ID,outID,(/N_edgV_ID/),'xpts_edgV','Double')
  Status = nf90_put_att(outID,xpts_edgV_ID,'bnds',(/0d0,xlen/))
  CALL NF_DEF_VAR(ypts_avg_ID,outID,(/N_y_ID/),'ypts_avg','Double')
  Status = nf90_put_att(outID,ypts_avg_ID,'bnds',(/0d0,ylen/))
  CALL NF_DEF_VAR(ypts_edgH_ID,outID,(/N_edgH_ID/),'ypts_edgH','Double')
  Status = nf90_put_att(outID,ypts_edgH_ID,'bnds',(/0d0,ylen/))
  CALL NF_DEF_VAR(tpts_ID,outID,(/N_t_ID/),'tpts','Double')
  Status = nf90_put_att(outID,tpts_ID,'bnds',(/tpts(1),tlen/))

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

  !
  !

  Status = nf90_put_var(outID,xpts_avg_ID,xpts_avg,(/1/),(/N_x/))
  CALL HANDLE_ERR(Status)

  Status = nf90_put_var(outID,xpts_edgV_ID,xpts_edgV,(/1/),(/N_x+1/))
  CALL HANDLE_ERR(Status)

  Status = nf90_put_var(outID,ypts_avg_ID,ypts_avg,(/1/),(/N_y/))
  CALL HANDLE_ERR(Status)

  Status = nf90_put_var(outID,ypts_edgH_ID,ypts_edgH,(/1/),(/N_y+1/))
  CALL HANDLE_ERR(Status)

  Status = nf90_put_var(outID,tpts_ID,tpts,(/1/),(/N_t/))
  CALL HANDLE_ERR(Status)

END SUBROUTINE OUTFILE_INIT

!==================================================================================================================================!
!
!==================================================================================================================================!
SUBROUTINE OUTFILE_INIT_Tgen(outID, N_x_ID, N_y_ID, N_m_ID, N_g_ID, N_t_ID, N_edgV_ID, N_edgH_ID, N_xc_ID, N_yc_ID, Quads_ID,&
  Boundaries_ID, c_ID, h_ID, pi_ID, erg_ID, Comp_Unit_ID, cv_ID, c, h, pi, erg, Comp_Unit, cv, xlen, ylen, Delx, Dely, tlen, Delt,&
  bcT_left, bcT_bottom, bcT_right, bcT_top, Tini, N_x, N_y, N_m, N_g, N_t, threads, BC_type, outfile, quadrature, enrgy_strc,&
  I_out, HO_Eg_out, HO_Fg_out, HO_E_out, HO_F_out, QDfg_out, xpts_avg, xpts_edgV, ypts_avg, ypts_edgH, tpts, run_type, Temp_ID,&
  HO_E_avg_ID, HO_E_edgV_ID, HO_E_edgH_ID, HO_Fx_edgV_ID, HO_Fy_edgH_ID, HO_Eg_avg_ID, HO_Eg_edgV_ID, HO_Eg_edgH_ID,&
  HO_Fxg_edgV_ID, HO_Fyg_edgH_ID, I_avg_ID, I_edgV_ID, I_edgH_ID, Cg_L_ID, Cg_B_ID, Cg_R_ID, Cg_T_ID, fg_avg_xx_ID, fg_avg_yy_ID,&
  fg_avg_xy_ID, fg_edgV_xx_ID, fg_edgV_xy_ID, fg_edgH_yy_ID, fg_edgH_xy_ID)

  INTEGER,INTENT(OUT):: outID
  INTEGER,INTENT(OUT):: N_x_ID, N_y_ID, N_m_ID, N_g_ID, N_t_ID, N_edgV_ID, N_edgH_ID, N_xc_ID, N_yc_ID, Quads_ID
  INTEGER,INTENT(OUT):: Boundaries_ID, c_ID, h_ID, pi_ID, erg_ID, Comp_Unit_ID, cv_ID
  INTEGER,INTENT(OUT):: Temp_ID, HO_E_avg_ID
  INTEGER,INTENT(OUT):: HO_E_edgV_ID, HO_E_edgH_ID, HO_Fx_edgV_ID
  INTEGER,INTENT(OUT):: HO_Fy_edgH_ID, HO_Eg_avg_ID, HO_Eg_edgV_ID, HO_Eg_edgH_ID
  INTEGER,INTENT(OUT):: HO_Fxg_edgV_ID, HO_Fyg_edgH_ID, I_avg_ID, I_edgV_ID, I_edgH_ID
  INTEGER,INTENT(OUT):: Cg_L_ID, Cg_B_ID, Cg_R_ID, Cg_T_ID
  INTEGER,INTENT(OUT):: fg_avg_xx_ID, fg_avg_yy_ID, fg_avg_xy_ID, fg_edgV_xx_ID, fg_edgV_xy_ID, fg_edgH_yy_ID, fg_edgH_xy_ID

  REAL*8,INTENT(IN):: c, h, pi, erg, Comp_Unit, cv
  REAL*8,INTENT(IN):: xlen, ylen, Delx(:), Dely(:), tlen, Delt
  REAL*8,INTENT(IN):: xpts_avg(:), xpts_edgV(:), ypts_avg(:), ypts_edgH(:), tpts(:)
  REAL*8,INTENT(IN):: bcT_left, bcT_bottom, bcT_right, bcT_top, Tini
  INTEGER,INTENT(IN):: N_x, N_y, N_m, N_g, N_t
  INTEGER,INTENT(IN):: threads, BC_type(:)
  INTEGER,INTENT(IN):: I_out, HO_Eg_out, HO_Fg_out, HO_E_out, HO_F_out
  INTEGER,INTENT(IN):: QDfg_out
  CHARACTER(*),INTENT(IN):: outfile, quadrature, enrgy_strc, run_type

  INTEGER:: Status
  INTEGER:: xlen_ID, ylen_ID
  INTEGER:: Delx_ID, Dely_ID, tlen_ID, Delt_ID, BC_type_ID, bcT_left_ID, bcT_bottom_ID, bcT_right_ID, bcT_top_ID, Tini_ID
  INTEGER:: xpts_avg_ID, xpts_edgV_ID, ypts_avg_ID, ypts_edgH_ID, tpts_ID
  INTEGER:: Z(0), i
  CHARACTER(50):: buf

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
  Status = nf90_put_att(outID,NF90_GLOBAL,'T_run_type',run_type)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_att(outID,NF90_GLOBAL,'quadrature',quadrature)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_att(outID,NF90_GLOBAL,'threads',threads)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_att(outID,NF90_GLOBAL,'enrgy_strc',enrgy_strc)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_att(outID,NF90_GLOBAL,'BC_type',BC_type)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_att(outID,NF90_GLOBAL,'I_out',I_out)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_att(outID,NF90_GLOBAL,'fg_out',QDfg_out)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_att(outID,NF90_GLOBAL,'HO_Eg_out',HO_Eg_out)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_att(outID,NF90_GLOBAL,'HO_E_out',HO_E_out)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_att(outID,NF90_GLOBAL,'HO_Fg_out',HO_Fg_out)
  CALL HANDLE_ERR(Status)

  Status = nf90_put_att(outID,NF90_GLOBAL,'HO_F_out',HO_F_out)
  CALL HANDLE_ERR(Status)

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

  CALL NF_DEF_VAR(xpts_avg_ID,outID,(/N_x_ID/),'xpts_avg','Double')
  Status = nf90_put_att(outID,xpts_avg_ID,'bnds',(/0d0,xlen/))
  CALL NF_DEF_VAR(xpts_edgV_ID,outID,(/N_edgV_ID/),'xpts_edgV','Double')
  Status = nf90_put_att(outID,xpts_edgV_ID,'bnds',(/0d0,xlen/))
  CALL NF_DEF_VAR(ypts_avg_ID,outID,(/N_y_ID/),'ypts_avg','Double')
  Status = nf90_put_att(outID,ypts_avg_ID,'bnds',(/0d0,ylen/))
  CALL NF_DEF_VAR(ypts_edgH_ID,outID,(/N_edgH_ID/),'ypts_edgH','Double')
  Status = nf90_put_att(outID,ypts_edgH_ID,'bnds',(/0d0,ylen/))
  CALL NF_DEF_VAR(tpts_ID,outID,(/N_t_ID/),'tpts','Double')
  Status = nf90_put_att(outID,tpts_ID,'bnds',(/tpts(1),tlen/))

  !===========================================================================!
  !                                                                           !
  !     DEFINING VARIABLES & UNITS                                            !
  !                                                                           !
  !===========================================================================!
  !--------------------------------------------------!
  !            Boundary Coefficients                 !
  !--------------------------------------------------!
  CALL NF_DEF_VAR(Cg_L_ID,outID,(/N_y_ID,N_g_ID,N_t_ID/),'Cg_L','Double')
  Status = nf90_put_att(outID,Cg_L_ID,'N_grids',1)
  Status = nf90_put_att(outID,Cg_L_ID,'grid0','tpts')
  CALL NF_DEF_VAR(Cg_B_ID,outID,(/N_x_ID,N_g_ID,N_t_ID/),'Cg_B','Double')
  Status = nf90_put_att(outID,Cg_B_ID,'N_grids',1)
  Status = nf90_put_att(outID,Cg_B_ID,'grid0','tpts')
  CALL NF_DEF_VAR(Cg_R_ID,outID,(/N_y_ID,N_g_ID,N_t_ID/),'Cg_R','Double')
  Status = nf90_put_att(outID,Cg_R_ID,'N_grids',1)
  Status = nf90_put_att(outID,Cg_R_ID,'grid0','tpts')
  CALL NF_DEF_VAR(Cg_T_ID,outID,(/N_x_ID,N_g_ID,N_t_ID/),'Cg_T','Double')
  Status = nf90_put_att(outID,Cg_T_ID,'N_grids',1)
  Status = nf90_put_att(outID,Cg_T_ID,'grid0','tpts')

  !--------------------------------------------------!
  !                 Temperature                      !
  !--------------------------------------------------!
  CALL NF_DEF_VAR(Temp_ID,outID,(/N_x_ID,N_y_ID,N_t_ID/),'Temperature','Double')
  CALL NF_DEF_UNIT(outID,Temp_ID,'ev/cm^3')

  !===========================================================================!
  !--------------------------------------------------!
  !         Total Radiation Energy Densities         !
  !--------------------------------------------------!
  IF (HO_E_out .EQ. 1) THEN
    !cell-averaged
    CALL NF_DEF_VAR(HO_E_avg_ID,outID,(/N_x_ID,N_y_ID,N_t_ID/),'E_avg_HO','Double') !from High-Order eqs
    CALL NF_DEF_UNIT(outID,HO_E_avg_ID,'ev/cm^3')
    Status = nf90_put_att(outID,HO_E_avg_ID,'N_grids',3)
    Status = nf90_put_att(outID,HO_E_avg_ID,'grid0','tpts')
    Status = nf90_put_att(outID,HO_E_avg_ID,'grid1','ypts_avg')
    Status = nf90_put_att(outID,HO_E_avg_ID,'grid2','xpts_avg')

    !'vertical' cell edges, x-const
    CALL NF_DEF_VAR(HO_E_edgV_ID,outID,(/N_edgV_ID,N_y_ID,N_t_ID/),'E_edgV_HO','Double') !from High-Order eqs
    CALL NF_DEF_UNIT(outID,HO_E_edgV_ID,'ev/cm^3')
    Status = nf90_put_att(outID,HO_E_edgV_ID,'N_grids',3)
    Status = nf90_put_att(outID,HO_E_edgV_ID,'grid0','tpts')
    Status = nf90_put_att(outID,HO_E_edgV_ID,'grid1','ypts_avg')
    Status = nf90_put_att(outID,HO_E_edgV_ID,'grid2','xpts_edgV')

    !'horizontal' cell edges, x-const
    CALL NF_DEF_VAR(HO_E_edgH_ID,outID,(/N_x_ID,N_edgH_ID,N_t_ID/),'E_edgH_HO','Double') !from High-Order eqs
    CALL NF_DEF_UNIT(outID,HO_E_edgH_ID,'ev/cm^3')
    Status = nf90_put_att(outID,HO_E_edgH_ID,'N_grids',3)
    Status = nf90_put_att(outID,HO_E_edgH_ID,'grid0','tpts')
    Status = nf90_put_att(outID,HO_E_edgH_ID,'grid1','ypts_edgH')
    Status = nf90_put_att(outID,HO_E_edgH_ID,'grid2','xpts_avg')

  END IF

  !===========================================================================!
  !--------------------------------------------------!
  !            Total Radiation Fluxes                !
  !--------------------------------------------------!
  IF (HO_F_out .EQ. 1) THEN
    !'vertical' cell edges, x-const
    CALL NF_DEF_VAR(HO_Fx_edgV_ID,outID,(/N_edgV_ID,N_y_ID,N_t_ID/),'Fx_edgV_HO','Double') !from High-Order eqs
    CALL NF_DEF_UNIT(outID,HO_Fx_edgV_ID,'ev/cm^3')
    Status = nf90_put_att(outID,HO_Fx_edgV_ID,'N_grids',3)
    Status = nf90_put_att(outID,HO_Fx_edgV_ID,'grid0','tpts')
    Status = nf90_put_att(outID,HO_Fx_edgV_ID,'grid1','ypts_avg')
    Status = nf90_put_att(outID,HO_Fx_edgV_ID,'grid2','xpts_edgV')

    !'horizontal' cell edges, x-const
    CALL NF_DEF_VAR(HO_Fy_edgH_ID,outID,(/N_x_ID,N_edgH_ID,N_t_ID/),'Fy_edgH_HO','Double') !from High-Order eqs
    CALL NF_DEF_UNIT(outID,HO_Fy_edgH_ID,'ev/cm^3')
    Status = nf90_put_att(outID,HO_Fy_edgH_ID,'N_grids',3)
    Status = nf90_put_att(outID,HO_Fy_edgH_ID,'grid0','tpts')
    Status = nf90_put_att(outID,HO_Fy_edgH_ID,'grid1','ypts_edgH')
    Status = nf90_put_att(outID,HO_Fy_edgH_ID,'grid2','xpts_avg')

  END IF

  !===========================================================================!
  !--------------------------------------------------!
  !     Multigroup Radiation Energy Densities        !
  !--------------------------------------------------!
  IF (HO_Eg_out .EQ. 1) THEN
    !cell-averaged
    CALL NF_DEF_VAR(HO_Eg_avg_ID,outID,(/N_x_ID,N_y_ID,N_g_ID,N_t_ID/),'Eg_avg_HO','Double') !from High-Order eqs
    CALL NF_DEF_UNIT(outID,HO_Eg_avg_ID,'ev/cm^3')
    Status = nf90_put_att(outID,HO_Eg_avg_ID,'N_grids',3)
    Status = nf90_put_att(outID,HO_Eg_avg_ID,'grid0','tpts')
    Status = nf90_put_att(outID,HO_Eg_avg_ID,'grid1','ypts_avg')
    Status = nf90_put_att(outID,HO_Eg_avg_ID,'grid2','xpts_avg')

    !'vertical' cell edges, x-const
    CALL NF_DEF_VAR(HO_Eg_edgV_ID,outID,(/N_edgV_ID,N_y_ID,N_g_ID,N_t_ID/),'Eg_edgV_HO','Double') !from High-Order eqs
    CALL NF_DEF_UNIT(outID,HO_Eg_edgV_ID,'ev/cm^3')
    Status = nf90_put_att(outID,HO_Eg_edgV_ID,'N_grids',3)
    Status = nf90_put_att(outID,HO_Eg_edgV_ID,'grid0','tpts')
    Status = nf90_put_att(outID,HO_Eg_edgV_ID,'grid1','ypts_avg')
    Status = nf90_put_att(outID,HO_Eg_edgV_ID,'grid2','xpts_edgV')

    !'horizontal' cell edges, x-const
    CALL NF_DEF_VAR(HO_Eg_edgH_ID,outID,(/N_x_ID,N_edgH_ID,N_g_ID,N_t_ID/),'Eg_edgH_HO','Double') !from High-Order eqs
    CALL NF_DEF_UNIT(outID,HO_Eg_edgH_ID,'ev/cm^3')
    Status = nf90_put_att(outID,HO_Eg_edgH_ID,'N_grids',3)
    Status = nf90_put_att(outID,HO_Eg_edgH_ID,'grid0','tpts')
    Status = nf90_put_att(outID,HO_Eg_edgH_ID,'grid1','ypts_edgH')
    Status = nf90_put_att(outID,HO_Eg_edgH_ID,'grid2','xpts_avg')

  END IF

  !===========================================================================!
  !--------------------------------------------------!
  !          Multigroup Radiation Fluxes             !
  !--------------------------------------------------!
  IF (HO_Fg_out .EQ. 1) THEN
    !'vertical' cell edges, x-const
    CALL NF_DEF_VAR(HO_Fxg_edgV_ID,outID,(/N_edgV_ID,N_y_ID,N_g_ID,N_t_ID/),'Fxg_edgV_HO','Double') !from High-Order eqs
    CALL NF_DEF_UNIT(outID,HO_Fxg_edgV_ID,'ev/cm^3')
    Status = nf90_put_att(outID,HO_Fxg_edgV_ID,'N_grids',3)
    Status = nf90_put_att(outID,HO_Fxg_edgV_ID,'grid0','tpts')
    Status = nf90_put_att(outID,HO_Fxg_edgV_ID,'grid1','ypts_avg')
    Status = nf90_put_att(outID,HO_Fxg_edgV_ID,'grid2','xpts_edgV')

    !'horizontal' cell edges, x-const
    CALL NF_DEF_VAR(HO_Fyg_edgH_ID,outID,(/N_x_ID,N_edgH_ID,N_g_ID,N_t_ID/),'Fyg_edgH_HO','Double') !from High-Order eqs
    CALL NF_DEF_UNIT(outID,HO_Fyg_edgH_ID,'ev/cm^3')
    Status = nf90_put_att(outID,HO_Fyg_edgH_ID,'N_grids',3)
    Status = nf90_put_att(outID,HO_Fyg_edgH_ID,'grid0','tpts')
    Status = nf90_put_att(outID,HO_Fyg_edgH_ID,'grid1','ypts_edgH')
    Status = nf90_put_att(outID,HO_Fyg_edgH_ID,'grid2','xpts_avg')

  END IF

  !===========================================================================!
  !--------------------------------------------------!
  !             Radiation Intensities                !
  !--------------------------------------------------!
  IF (I_out .EQ. 1) THEN
    !cell-averaged
    CALL NF_DEF_VAR(I_avg_ID,outID,(/N_x_ID,N_y_ID,N_m_ID,N_g_ID,N_t_ID/),'I_avg','Double')
    CALL NF_DEF_UNIT(outID,I_avg_ID,'erg/(Ster*cm^3)')
    Status = nf90_put_att(outID,I_avg_ID,'N_grids',3)
    Status = nf90_put_att(outID,I_avg_ID,'grid0','tpts')
    Status = nf90_put_att(outID,I_avg_ID,'grid1','ypts_avg')
    Status = nf90_put_att(outID,I_avg_ID,'grid2','xpts_avg')

    !'vertical' cell edges, x-const
    CALL NF_DEF_VAR(I_edgV_ID,outID,(/N_edgV_ID,N_y_ID,N_m_ID,N_g_ID,N_t_ID/),'I_edgV','Double')
    CALL NF_DEF_UNIT(outID,I_edgV_ID,'erg/(Ster*cm^3)')
    Status = nf90_put_att(outID,I_edgV_ID,'N_grids',3)
    Status = nf90_put_att(outID,I_edgV_ID,'grid0','tpts')
    Status = nf90_put_att(outID,I_edgV_ID,'grid1','ypts_avg')
    Status = nf90_put_att(outID,I_edgV_ID,'grid2','xpts_edgV')

    !'horizontal' cell edges, x-const
    CALL NF_DEF_VAR(I_edgH_ID,outID,(/N_x_ID,N_edgH_ID,N_m_ID,N_g_ID,N_t_ID/),'I_edgH','Double')
    CALL NF_DEF_UNIT(outID,I_edgH_ID,'erg/(Ster*cm^3)')
    Status = nf90_put_att(outID,I_edgH_ID,'N_grids',3)
    Status = nf90_put_att(outID,I_edgH_ID,'grid0','tpts')
    Status = nf90_put_att(outID,I_edgH_ID,'grid1','ypts_edgH')
    Status = nf90_put_att(outID,I_edgH_ID,'grid2','xpts_avg')

  END IF

  !===========================================================================!
  !--------------------------------------------------!
  !                   QD factors                     !
  !--------------------------------------------------!
  IF (QDfg_out .EQ. 1) THEN
    CALL NF_DEF_VAR(fg_avg_xx_ID,outID,(/N_x_ID,N_y_ID,N_g_ID,N_t_ID/),'fg_avg_xx','Double')
    Status = nf90_put_att(outID,fg_avg_xx_ID,'N_grids',3)
    Status = nf90_put_att(outID,fg_avg_xx_ID,'grid0','tpts')
    Status = nf90_put_att(outID,fg_avg_xx_ID,'grid1','ypts_avg')
    Status = nf90_put_att(outID,fg_avg_xx_ID,'grid2','xpts_avg')

    CALL NF_DEF_VAR(fg_avg_yy_ID,outID,(/N_x_ID,N_y_ID,N_g_ID,N_t_ID/),'fg_avg_yy','Double')
    Status = nf90_put_att(outID,fg_avg_yy_ID,'N_grids',3)
    Status = nf90_put_att(outID,fg_avg_yy_ID,'grid0','tpts')
    Status = nf90_put_att(outID,fg_avg_yy_ID,'grid1','ypts_avg')
    Status = nf90_put_att(outID,fg_avg_yy_ID,'grid2','xpts_avg')

    CALL NF_DEF_VAR(fg_avg_xy_ID,outID,(/N_x_ID,N_y_ID,N_g_ID,N_t_ID/),'fg_avg_xy','Double')
    Status = nf90_put_att(outID,fg_avg_xy_ID,'N_grids',3)
    Status = nf90_put_att(outID,fg_avg_xy_ID,'grid0','tpts')
    Status = nf90_put_att(outID,fg_avg_xy_ID,'grid1','ypts_avg')
    Status = nf90_put_att(outID,fg_avg_xy_ID,'grid2','xpts_avg')

    CALL NF_DEF_VAR(fg_edgV_xx_ID,outID,(/N_edgV_ID,N_y_ID,N_g_ID,N_t_ID/),'fg_edgV_xx','Double')
    Status = nf90_put_att(outID,fg_edgV_xx_ID,'N_grids',3)
    Status = nf90_put_att(outID,fg_edgV_xx_ID,'grid0','tpts')
    Status = nf90_put_att(outID,fg_edgV_xx_ID,'grid1','ypts_avg')
    Status = nf90_put_att(outID,fg_edgV_xx_ID,'grid2','xpts_edgV')

    CALL NF_DEF_VAR(fg_edgV_xy_ID,outID,(/N_edgV_ID,N_y_ID,N_g_ID,N_t_ID/),'fg_edgV_xy','Double')
    Status = nf90_put_att(outID,fg_edgV_xy_ID,'N_grids',3)
    Status = nf90_put_att(outID,fg_edgV_xy_ID,'grid0','tpts')
    Status = nf90_put_att(outID,fg_edgV_xy_ID,'grid1','ypts_avg')
    Status = nf90_put_att(outID,fg_edgV_xy_ID,'grid2','xpts_edgV')

    CALL NF_DEF_VAR(fg_edgH_yy_ID,outID,(/N_x_ID,N_edgH_ID,N_g_ID,N_t_ID/),'fg_edgH_yy','Double')
    Status = nf90_put_att(outID,fg_edgH_yy_ID,'N_grids',3)
    Status = nf90_put_att(outID,fg_edgH_yy_ID,'grid0','tpts')
    Status = nf90_put_att(outID,fg_edgH_yy_ID,'grid1','ypts_edgH')
    Status = nf90_put_att(outID,fg_edgH_yy_ID,'grid2','xpts_avg')

    CALL NF_DEF_VAR(fg_edgH_xy_ID,outID,(/N_x_ID,N_edgH_ID,N_g_ID,N_t_ID/),'fg_edgH_xy','Double')
    Status = nf90_put_att(outID,fg_edgH_xy_ID,'N_grids',3)
    Status = nf90_put_att(outID,fg_edgH_xy_ID,'grid0','tpts')
    Status = nf90_put_att(outID,fg_edgH_xy_ID,'grid1','ypts_edgH')
    Status = nf90_put_att(outID,fg_edgH_xy_ID,'grid2','xpts_avg')

  END IF

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

  !
  !

  Status = nf90_put_var(outID,xpts_avg_ID,xpts_avg,(/1/),(/N_x/))
  CALL HANDLE_ERR(Status)

  Status = nf90_put_var(outID,xpts_edgV_ID,xpts_edgV,(/1/),(/N_x+1/))
  CALL HANDLE_ERR(Status)

  Status = nf90_put_var(outID,ypts_avg_ID,ypts_avg,(/1/),(/N_y/))
  CALL HANDLE_ERR(Status)

  Status = nf90_put_var(outID,ypts_edgH_ID,ypts_edgH,(/1/),(/N_y+1/))
  CALL HANDLE_ERR(Status)

  Status = nf90_put_var(outID,tpts_ID,tpts,(/1/),(/N_t/))
  CALL HANDLE_ERR(Status)

END SUBROUTINE OUTFILE_INIT_Tgen

!==================================================================================================================================!
!
!==================================================================================================================================!
SUBROUTINE OUTFILE_VARDEFS(outID,Res_Calc,out_freq,I_out,HO_Eg_out,HO_Fg_out,HO_E_out,HO_F_out,Eg_out,Fg_out,MGQD_E_out,&
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

  INTEGER,INTENT(IN):: outID
  LOGICAL,INTENT(IN):: Res_Calc
  INTEGER,INTENT(IN):: out_freq, I_out, HO_Eg_out, HO_Fg_out, HO_E_out, HO_F_out
  INTEGER,INTENT(IN):: Eg_out, Fg_out, MGQD_E_out, MGQD_F_out, QDfg_out, fg_pod_out, fg_dmd_out, N_dsets_ID
  INTEGER,INTENT(IN):: E_out, F_out, D_out
  INTEGER,INTENT(IN):: old_parms_out, its_out, conv_out, kap_out, Src_out
  INTEGER,INTENT(IN):: N_x_ID, N_y_ID, N_m_ID, N_g_ID, N_t_ID, N_edgV_ID, N_edgH_ID, N_xc_ID, N_yc_ID, Quads_ID
  INTEGER,INTENT(IN):: RT_Its_ID, MGQD_Its_ID, GQD_Its_ID, Norm_Types_ID, MGQD_ResTypes_ID, Boundaries_ID
  INTEGER,INTENT(OUT):: Temp_ID, E_avg_ID, E_edgV_ID, E_edgH_ID, MGQD_E_avg_ID, MGQD_E_edgV_ID, MGQD_E_edgH_ID, HO_E_avg_ID
  INTEGER,INTENT(OUT):: HO_E_edgV_ID, HO_E_edgH_ID, Fx_edgV_ID, Fy_edgH_ID, MGQD_Fx_edgV_ID, MGQD_Fy_edgH_ID, HO_Fx_edgV_ID
  INTEGER,INTENT(OUT):: HO_Fy_edgH_ID, Eg_avg_ID, Eg_edgV_ID, Eg_edgH_ID, HO_Eg_avg_ID, HO_Eg_edgV_ID, HO_Eg_edgH_ID
  INTEGER,INTENT(OUT):: Fxg_edgV_ID, Fyg_edgH_ID, HO_Fxg_edgV_ID, HO_Fyg_edgH_ID, I_avg_ID, I_edgV_ID, I_edgH_ID
  INTEGER,INTENT(OUT):: KapE_Bar_ID, KapB_ID, KapE_ID, KapR_ID, Bg_ID
  INTEGER,INTENT(OUT):: RT_Residual_ID, MGQD_Residual_ID, MGQD_BC_Residual_ID
  INTEGER,INTENT(OUT):: Del_T_ID, Del_E_avg_ID, Del_E_edgV_ID, Del_E_edgH_ID, Del_Fx_edgV_ID, Del_Fy_edgH_ID
  INTEGER,INTENT(OUT):: RT_ItCount_ID, MGQD_ItCount_ID, GQD_ItCount_ID, MGQD_KIts_ID, GQD_KIts_ID
  INTEGER,INTENT(OUT):: RT_Tnorm_ID, RT_Enorm_ID, MGQD_Tnorm_ID, MGQD_Enorm_ID
  INTEGER,INTENT(OUT):: RT_Trho_ID, RT_Erho_ID, MGQD_Trho_ID, MGQD_Erho_ID
  INTEGER,INTENT(OUT):: Cg_L_ID, Cg_B_ID, Cg_R_ID, Cg_T_ID, Eg_in_L_ID, Eg_in_B_ID, Eg_in_R_ID, Eg_in_T_ID
  INTEGER,INTENT(OUT):: Fg_in_L_ID, Fg_in_B_ID, Fg_in_R_ID, Fg_in_T_ID
  INTEGER,INTENT(OUT):: Cb_L_ID, Cb_B_ID, Cb_R_ID, Cb_T_ID, E_in_L_ID, E_in_B_ID, E_in_R_ID, E_in_T_ID
  INTEGER,INTENT(OUT):: F_in_L_ID, F_in_B_ID, F_in_R_ID, F_in_T_ID
  INTEGER,INTENT(OUT):: fg_avg_xx_ID, fg_avg_yy_ID, fg_avg_xy_ID, fg_edgV_xx_ID, fg_edgV_xy_ID, fg_edgH_yy_ID, fg_edgH_xy_ID
  INTEGER,INTENT(OUT):: DC_xx_ID, DL_xx_ID, DR_xx_ID, DC_yy_ID, DB_yy_ID, DT_yy_ID, DL_xy_ID, DB_xy_ID, DR_xy_ID, DT_xy_ID
  INTEGER,INTENT(OUT):: G_old_ID, Pold_L_ID, Pold_B_ID, Pold_R_ID, Pold_T_ID
  INTEGER,INTENT(OUT):: Gold_hat_ID, Rhat_old_ID, PL_ID, PB_ID, PR_ID, PT_ID
  INTEGER,INTENT(OUT):: dr_T_ID, dr_B_ID, dr_ML_ID, dr_MB_ID, dr_MR_ID, dr_MT_ID
  INTEGER,INTENT(OUT):: rrank_BCg_ID, rrank_fg_avg_xx_ID, rrank_fg_edgV_xx_ID, rrank_fg_avg_yy_ID, rrank_fg_edgH_yy_ID
  INTEGER,INTENT(OUT):: rrank_fg_edgV_xy_ID, rrank_fg_edgH_xy_ID, DMDgsum_ID, PODgsum_ID, PODerr_ID
  INTEGER,INTENT(OUT):: Temp_XWvSpeed_ID, E_XWvSpeed_ID

  INTEGER:: Status
  CHARACTER(50):: Location = 'MODULE: OUTPUTS / SUBROUTINE: OUTFILE_VARDEFS'

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
  IF (its_out .EQ. 1) THEN
    CALL NF_DEF_VAR(RT_ItCount_ID,outID,(/N_t_ID/),'RT_Its-count','Int')
    CALL NF_DEF_VAR(MGQD_ItCount_ID,outID,(/RT_Its_ID,N_t_ID/),'MGQD_Its-count','Int')
    CALL NF_DEF_VAR(GQD_ItCount_ID,outID,(/MGQD_Its_ID,RT_Its_ID,N_t_ID/),'GQD_Its-count','Int')

    CALL NF_DEF_VAR(MGQD_KIts_ID,outID,(/N_g_ID,RT_Its_ID,N_t_ID/),'MGQD_MatVec_prds','Int')
    CALL NF_DEF_VAR(GQD_KIts_ID,outID,(/GQD_Its_ID,MGQD_Its_ID,RT_Its_ID,N_t_ID/),'GQD_MatVec_prds','Int')
  END IF

  IF (conv_out .EQ. 1) THEN
    CALL NF_DEF_VAR(RT_Tnorm_ID,outID,(/RT_Its_ID,N_t_ID/),'RT_Tnorm','Real')
    CALL NF_DEF_VAR(RT_Enorm_ID,outID,(/RT_Its_ID,N_t_ID/),'RT_Enorm','Real')
    CALL NF_DEF_VAR(MGQD_Tnorm_ID,outID,(/MGQD_Its_ID,RT_Its_ID,N_t_ID/),'MGQD_Tnorm','Real')
    CALL NF_DEF_VAR(MGQD_Enorm_ID,outID,(/MGQD_Its_ID,RT_Its_ID,N_t_ID/),'MGQD_Enorm','Real')

    CALL NF_DEF_VAR(RT_Trho_ID,outID,(/RT_Its_ID,N_t_ID/),'RT_Trho','Real')
    CALL NF_DEF_VAR(RT_Erho_ID,outID,(/RT_Its_ID,N_t_ID/),'RT_Erho','Real')
    CALL NF_DEF_VAR(MGQD_Trho_ID,outID,(/MGQD_Its_ID,RT_Its_ID,N_t_ID/),'MGQD_Trho','Real')
    CALL NF_DEF_VAR(MGQD_Erho_ID,outID,(/MGQD_Its_ID,RT_Its_ID,N_t_ID/),'MGQD_Erho','Real')
  END IF

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
  Status = nf90_put_att(outID,Cg_L_ID,'N_grids',1)
  Status = nf90_put_att(outID,Cg_L_ID,'grid0','tpts')
  CALL NF_DEF_VAR(Cg_B_ID,outID,(/N_x_ID,N_g_ID,N_t_ID/),'Cg_B','Double')
  Status = nf90_put_att(outID,Cg_B_ID,'N_grids',1)
  Status = nf90_put_att(outID,Cg_B_ID,'grid0','tpts')
  CALL NF_DEF_VAR(Cg_R_ID,outID,(/N_y_ID,N_g_ID,N_t_ID/),'Cg_R','Double')
  Status = nf90_put_att(outID,Cg_R_ID,'N_grids',1)
  Status = nf90_put_att(outID,Cg_R_ID,'grid0','tpts')
  CALL NF_DEF_VAR(Cg_T_ID,outID,(/N_x_ID,N_g_ID,N_t_ID/),'Cg_T','Double')
  Status = nf90_put_att(outID,Cg_T_ID,'N_grids',1)
  Status = nf90_put_att(outID,Cg_T_ID,'grid0','tpts')

  CALL NF_DEF_VAR(Eg_in_L_ID,outID,(/N_y_ID,N_g_ID,N_t_ID/),'Eg_in_L','Double')
  CALL NF_DEF_VAR(Eg_in_B_ID,outID,(/N_x_ID,N_g_ID,N_t_ID/),'Eg_in_B','Double')
  CALL NF_DEF_VAR(Eg_in_R_ID,outID,(/N_y_ID,N_g_ID,N_t_ID/),'Eg_in_R','Double')
  CALL NF_DEF_VAR(Eg_in_T_ID,outID,(/N_x_ID,N_g_ID,N_t_ID/),'Eg_in_T','Double')

  CALL NF_DEF_VAR(Fg_in_L_ID,outID,(/N_y_ID,N_g_ID,N_t_ID/),'Fg_in_L','Double')
  CALL NF_DEF_VAR(Fg_in_B_ID,outID,(/N_x_ID,N_g_ID,N_t_ID/),'Fg_in_B','Double')
  CALL NF_DEF_VAR(Fg_in_R_ID,outID,(/N_y_ID,N_g_ID,N_t_ID/),'Fg_in_R','Double')
  CALL NF_DEF_VAR(Fg_in_T_ID,outID,(/N_x_ID,N_g_ID,N_t_ID/),'Fg_in_T','Double')

  !===========================================================================!
  !--------------------------------------------------!
  !                 Wave Speeds                      !
  !--------------------------------------------------!
  CALL NF_DEF_VAR(Temp_XWvSpeed_ID,outID,(/N_t_ID/),'Temp_xWave_Speed','Double')
  CALL NF_DEF_VAR(E_XWvSpeed_ID,outID,(/N_t_ID/),'E_xWave_Speed','Double')

  !===========================================================================!
  !--------------------------------------------------!
  !                 Temperature                      !
  !--------------------------------------------------!
  CALL NF_DEF_VAR(Temp_ID,outID,(/N_x_ID,N_y_ID,N_t_ID/),'Temperature','Double')
  CALL NF_DEF_UNIT(outID,Temp_ID,'ev/cm^3')

  !===========================================================================!
  !--------------------------------------------------!
  !         Total Radiation Energy Densities         !
  !--------------------------------------------------!
  !------------------(grey loqd)---------------------!
  IF (E_out .EQ. 1) THEN
    !cell-averaged
    CALL NF_DEF_VAR(E_avg_ID,outID,(/N_x_ID,N_y_ID,N_t_ID/),'E_avg_Grey','Double') !from GQD eqs
    CALL NF_DEF_UNIT(outID,E_avg_ID,'ev/cm^3')
    Status = nf90_put_att(outID,E_avg_ID,'N_grids',3)
    Status = nf90_put_att(outID,E_avg_ID,'grid0','tpts')
    Status = nf90_put_att(outID,E_avg_ID,'grid1','ypts_avg')
    Status = nf90_put_att(outID,E_avg_ID,'grid2','xpts_avg')

    !'vertical' cell edges, x-const
    CALL NF_DEF_VAR(E_edgV_ID,outID,(/N_edgV_ID,N_y_ID,N_t_ID/),'E_edgV_Grey','Double') !from GQD eqs
    CALL NF_DEF_UNIT(outID,E_edgV_ID,'ev/cm^3')
    Status = nf90_put_att(outID,E_edgV_ID,'N_grids',3)
    Status = nf90_put_att(outID,E_edgV_ID,'grid0','tpts')
    Status = nf90_put_att(outID,E_edgV_ID,'grid1','ypts_avg')
    Status = nf90_put_att(outID,E_edgV_ID,'grid2','xpts_edgV')

    !'horizontal' cell edges, x-const
    CALL NF_DEF_VAR(E_edgH_ID,outID,(/N_x_ID,N_edgH_ID,N_t_ID/),'E_edgH_Grey','Double') !from GQD eqs
    CALL NF_DEF_UNIT(outID,E_edgH_ID,'ev/cm^3')
    Status = nf90_put_att(outID,E_edgH_ID,'N_grids',3)
    Status = nf90_put_att(outID,E_edgH_ID,'grid0','tpts')
    Status = nf90_put_att(outID,E_edgH_ID,'grid1','ypts_edgH')
    Status = nf90_put_att(outID,E_edgH_ID,'grid2','xpts_avg')

  END IF

  !---------------(multigroup loqd)------------------!
  IF (MGQD_E_out .EQ. 1) THEN
    !cell-averaged
    CALL NF_DEF_VAR(MGQD_E_avg_ID,outID,(/N_x_ID,N_y_ID,N_t_ID/),'E_avg_MGQD','Double') !from MGQD eqs
    CALL NF_DEF_UNIT(outID,MGQD_E_avg_ID,'ev/cm^3')
    Status = nf90_put_att(outID,MGQD_E_avg_ID,'N_grids',3)
    Status = nf90_put_att(outID,MGQD_E_avg_ID,'grid0','tpts')
    Status = nf90_put_att(outID,MGQD_E_avg_ID,'grid1','ypts_avg')
    Status = nf90_put_att(outID,MGQD_E_avg_ID,'grid2','xpts_avg')

    !'vertical' cell edges, x-const
    CALL NF_DEF_VAR(MGQD_E_edgV_ID,outID,(/N_edgV_ID,N_y_ID,N_t_ID/),'E_edgV_MGQD','Double') !from MGQD eqs
    CALL NF_DEF_UNIT(outID,MGQD_E_edgV_ID,'ev/cm^3')
    Status = nf90_put_att(outID,MGQD_E_edgV_ID,'N_grids',3)
    Status = nf90_put_att(outID,MGQD_E_edgV_ID,'grid0','tpts')
    Status = nf90_put_att(outID,MGQD_E_edgV_ID,'grid1','ypts_avg')
    Status = nf90_put_att(outID,MGQD_E_edgV_ID,'grid2','xpts_edgV')

    !'horizontal' cell edges, x-const
    CALL NF_DEF_VAR(MGQD_E_edgH_ID,outID,(/N_x_ID,N_edgH_ID,N_t_ID/),'E_edgH_MGQD','Double') !from MGQD eqs
    CALL NF_DEF_UNIT(outID,MGQD_E_edgH_ID,'ev/cm^3')
    Status = nf90_put_att(outID,MGQD_E_edgH_ID,'N_grids',3)
    Status = nf90_put_att(outID,MGQD_E_edgH_ID,'grid0','tpts')
    Status = nf90_put_att(outID,MGQD_E_edgH_ID,'grid1','ypts_edgH')
    Status = nf90_put_att(outID,MGQD_E_edgH_ID,'grid2','xpts_avg')

  END IF

  !------------------(high-order)--------------------!
  IF (HO_E_out .EQ. 1) THEN
    !cell-averaged
    CALL NF_DEF_VAR(HO_E_avg_ID,outID,(/N_x_ID,N_y_ID,N_t_ID/),'E_avg_HO','Double') !from High-Order eqs
    CALL NF_DEF_UNIT(outID,HO_E_avg_ID,'ev/cm^3')
    Status = nf90_put_att(outID,HO_E_avg_ID,'N_grids',3)
    Status = nf90_put_att(outID,HO_E_avg_ID,'grid0','tpts')
    Status = nf90_put_att(outID,HO_E_avg_ID,'grid1','ypts_avg')
    Status = nf90_put_att(outID,HO_E_avg_ID,'grid2','xpts_avg')

    !'vertical' cell edges, x-const
    CALL NF_DEF_VAR(HO_E_edgV_ID,outID,(/N_edgV_ID,N_y_ID,N_t_ID/),'E_edgV_HO','Double') !from High-Order eqs
    CALL NF_DEF_UNIT(outID,HO_E_edgV_ID,'ev/cm^3')
    Status = nf90_put_att(outID,HO_E_edgV_ID,'N_grids',3)
    Status = nf90_put_att(outID,HO_E_edgV_ID,'grid0','tpts')
    Status = nf90_put_att(outID,HO_E_edgV_ID,'grid1','ypts_avg')
    Status = nf90_put_att(outID,HO_E_edgV_ID,'grid2','xpts_edgV')

    !'horizontal' cell edges, x-const
    CALL NF_DEF_VAR(HO_E_edgH_ID,outID,(/N_x_ID,N_edgH_ID,N_t_ID/),'E_edgH_HO','Double') !from High-Order eqs
    CALL NF_DEF_UNIT(outID,HO_E_edgH_ID,'ev/cm^3')
    Status = nf90_put_att(outID,HO_E_edgH_ID,'N_grids',3)
    Status = nf90_put_att(outID,HO_E_edgH_ID,'grid0','tpts')
    Status = nf90_put_att(outID,HO_E_edgH_ID,'grid1','ypts_edgH')
    Status = nf90_put_att(outID,HO_E_edgH_ID,'grid2','xpts_avg')

  END IF

  !===========================================================================!
  !--------------------------------------------------!
  !            Total Radiation Fluxes                !
  !--------------------------------------------------!
  !------------------(grey loqd)---------------------!
  IF (F_out .EQ. 1) THEN
    !'vertical' cell edges, x-const
    CALL NF_DEF_VAR(Fx_edgV_ID,outID,(/N_edgV_ID,N_y_ID,N_t_ID/),'Fx_edgV_Grey','Double') !from GQD eqs
    CALL NF_DEF_UNIT(outID,Fx_edgV_ID,'ev/cm^3')
    Status = nf90_put_att(outID,Fx_edgV_ID,'N_grids',3)
    Status = nf90_put_att(outID,Fx_edgV_ID,'grid0','tpts')
    Status = nf90_put_att(outID,Fx_edgV_ID,'grid1','ypts_avg')
    Status = nf90_put_att(outID,Fx_edgV_ID,'grid2','xpts_edgV')

    !'horizontal' cell edges, x-const
    CALL NF_DEF_VAR(Fy_edgH_ID,outID,(/N_x_ID,N_edgH_ID,N_t_ID/),'Fy_edgH_Grey','Double') !from GQD eqs
    CALL NF_DEF_UNIT(outID,Fy_edgH_ID,'ev/cm^3')
    Status = nf90_put_att(outID,Fy_edgH_ID,'N_grids',3)
    Status = nf90_put_att(outID,Fy_edgH_ID,'grid0','tpts')
    Status = nf90_put_att(outID,Fy_edgH_ID,'grid1','ypts_edgH')
    Status = nf90_put_att(outID,Fy_edgH_ID,'grid2','xpts_avg')

  END IF

  !---------------(multigroup loqd)------------------!
  IF (MGQD_F_out .EQ. 1) THEN
    !'vertical' cell edges, x-const
    CALL NF_DEF_VAR(MGQD_Fx_edgV_ID,outID,(/N_edgV_ID,N_y_ID,N_t_ID/),'Fx_edgV_MGQD','Double') !from MGQD eqs
    CALL NF_DEF_UNIT(outID,MGQD_Fx_edgV_ID,'ev/cm^3')
    Status = nf90_put_att(outID,MGQD_Fx_edgV_ID,'N_grids',3)
    Status = nf90_put_att(outID,MGQD_Fx_edgV_ID,'grid0','tpts')
    Status = nf90_put_att(outID,MGQD_Fx_edgV_ID,'grid1','ypts_avg')
    Status = nf90_put_att(outID,MGQD_Fx_edgV_ID,'grid2','xpts_edgV')

    !'horizontal' cell edges, x-const
    CALL NF_DEF_VAR(MGQD_Fy_edgH_ID,outID,(/N_x_ID,N_edgH_ID,N_t_ID/),'Fy_edgH_MGQD','Double') !from MGQD eqs
    CALL NF_DEF_UNIT(outID,MGQD_Fy_edgH_ID,'ev/cm^3')
    Status = nf90_put_att(outID,MGQD_Fy_edgH_ID,'N_grids',3)
    Status = nf90_put_att(outID,MGQD_Fy_edgH_ID,'grid0','tpts')
    Status = nf90_put_att(outID,MGQD_Fy_edgH_ID,'grid1','ypts_edgH')
    Status = nf90_put_att(outID,MGQD_Fy_edgH_ID,'grid2','xpts_avg')


  END IF

  !------------------(high-order)--------------------!
  IF (HO_F_out .EQ. 1) THEN
    !'vertical' cell edges, x-const
    CALL NF_DEF_VAR(HO_Fx_edgV_ID,outID,(/N_edgV_ID,N_y_ID,N_t_ID/),'Fx_edgV_HO','Double') !from High-Order eqs
    CALL NF_DEF_UNIT(outID,HO_Fx_edgV_ID,'ev/cm^3')
    Status = nf90_put_att(outID,HO_Fx_edgV_ID,'N_grids',3)
    Status = nf90_put_att(outID,HO_Fx_edgV_ID,'grid0','tpts')
    Status = nf90_put_att(outID,HO_Fx_edgV_ID,'grid1','ypts_avg')
    Status = nf90_put_att(outID,HO_Fx_edgV_ID,'grid2','xpts_edgV')

    !'horizontal' cell edges, x-const
    CALL NF_DEF_VAR(HO_Fy_edgH_ID,outID,(/N_x_ID,N_edgH_ID,N_t_ID/),'Fy_edgH_HO','Double') !from High-Order eqs
    CALL NF_DEF_UNIT(outID,HO_Fy_edgH_ID,'ev/cm^3')
    Status = nf90_put_att(outID,HO_Fy_edgH_ID,'N_grids',3)
    Status = nf90_put_att(outID,HO_Fy_edgH_ID,'grid0','tpts')
    Status = nf90_put_att(outID,HO_Fy_edgH_ID,'grid1','ypts_edgH')
    Status = nf90_put_att(outID,HO_Fy_edgH_ID,'grid2','xpts_avg')

  END IF

  !===========================================================================!
  !--------------------------------------------------!
  !     Multigroup Radiation Energy Densities        !
  !--------------------------------------------------!
  !---------------(multigroup loqd)------------------!
  IF (Eg_out .EQ. 1) THEN
    !cell-averaged
    CALL NF_DEF_VAR(Eg_avg_ID,outID,(/N_x_ID,N_y_ID,N_g_ID,N_t_ID/),'Eg_avg_MGQD','Double') !from MGQD eqs
    CALL NF_DEF_UNIT(outID,Eg_avg_ID,'ev/cm^3')
    Status = nf90_put_att(outID,Eg_avg_ID,'N_grids',3)
    Status = nf90_put_att(outID,Eg_avg_ID,'grid0','tpts')
    Status = nf90_put_att(outID,Eg_avg_ID,'grid1','ypts_avg')
    Status = nf90_put_att(outID,Eg_avg_ID,'grid2','xpts_avg')

    !'vertical' cell edges, x-const
    CALL NF_DEF_VAR(Eg_edgV_ID,outID,(/N_edgV_ID,N_y_ID,N_g_ID,N_t_ID/),'Eg_edgV_MGQD','Double') !from MGQD eqs
    CALL NF_DEF_UNIT(outID,Eg_edgV_ID,'ev/cm^3')
    Status = nf90_put_att(outID,Eg_edgV_ID,'N_grids',3)
    Status = nf90_put_att(outID,Eg_edgV_ID,'grid0','tpts')
    Status = nf90_put_att(outID,Eg_edgV_ID,'grid1','ypts_avg')
    Status = nf90_put_att(outID,Eg_edgV_ID,'grid2','xpts_edgV')

    !'horizontal' cell edges, x-const
    CALL NF_DEF_VAR(Eg_edgH_ID,outID,(/N_x_ID,N_edgH_ID,N_g_ID,N_t_ID/),'Eg_edgH_MGQD','Double') !from MGQD eqs
    CALL NF_DEF_UNIT(outID,Eg_edgH_ID,'ev/cm^3')
    Status = nf90_put_att(outID,Eg_edgH_ID,'N_grids',3)
    Status = nf90_put_att(outID,Eg_edgH_ID,'grid0','tpts')
    Status = nf90_put_att(outID,Eg_edgH_ID,'grid1','ypts_edgH')
    Status = nf90_put_att(outID,Eg_edgH_ID,'grid2','xpts_avg')

  END IF

  !------------------(high-order)--------------------!
  IF (HO_Eg_out .EQ. 1) THEN
    !cell-averaged
    CALL NF_DEF_VAR(HO_Eg_avg_ID,outID,(/N_x_ID,N_y_ID,N_g_ID,N_t_ID/),'Eg_avg_HO','Double') !from High-Order eqs
    CALL NF_DEF_UNIT(outID,HO_Eg_avg_ID,'ev/cm^3')
    Status = nf90_put_att(outID,HO_Eg_avg_ID,'N_grids',3)
    Status = nf90_put_att(outID,HO_Eg_avg_ID,'grid0','tpts')
    Status = nf90_put_att(outID,HO_Eg_avg_ID,'grid1','ypts_avg')
    Status = nf90_put_att(outID,HO_Eg_avg_ID,'grid2','xpts_avg')

    !'vertical' cell edges, x-const
    CALL NF_DEF_VAR(HO_Eg_edgV_ID,outID,(/N_edgV_ID,N_y_ID,N_g_ID,N_t_ID/),'Eg_edgV_HO','Double') !from High-Order eqs
    CALL NF_DEF_UNIT(outID,HO_Eg_edgV_ID,'ev/cm^3')
    Status = nf90_put_att(outID,HO_Eg_edgV_ID,'N_grids',3)
    Status = nf90_put_att(outID,HO_Eg_edgV_ID,'grid0','tpts')
    Status = nf90_put_att(outID,HO_Eg_edgV_ID,'grid1','ypts_avg')
    Status = nf90_put_att(outID,HO_Eg_edgV_ID,'grid2','xpts_edgV')

    !'horizontal' cell edges, x-const
    CALL NF_DEF_VAR(HO_Eg_edgH_ID,outID,(/N_x_ID,N_edgH_ID,N_g_ID,N_t_ID/),'Eg_edgH_HO','Double') !from High-Order eqs
    CALL NF_DEF_UNIT(outID,HO_Eg_edgH_ID,'ev/cm^3')
    Status = nf90_put_att(outID,HO_Eg_edgH_ID,'N_grids',3)
    Status = nf90_put_att(outID,HO_Eg_edgH_ID,'grid0','tpts')
    Status = nf90_put_att(outID,HO_Eg_edgH_ID,'grid1','ypts_edgH')
    Status = nf90_put_att(outID,HO_Eg_edgH_ID,'grid2','xpts_avg')

  END IF

  !===========================================================================!
  !--------------------------------------------------!
  !          Multigroup Radiation Fluxes             !
  !--------------------------------------------------!
  !---------------(multigroup loqd)------------------!
  IF (Fg_out .EQ. 1) THEN
    !'vertical' cell edges, x-const
    CALL NF_DEF_VAR(Fxg_edgV_ID,outID,(/N_edgV_ID,N_y_ID,N_g_ID,N_t_ID/),'Fxg_edgV_MGQD','Double') !from MGQD eqs
    CALL NF_DEF_UNIT(outID,Fxg_edgV_ID,'ev/cm^3')
    Status = nf90_put_att(outID,Fxg_edgV_ID,'N_grids',3)
    Status = nf90_put_att(outID,Fxg_edgV_ID,'grid0','tpts')
    Status = nf90_put_att(outID,Fxg_edgV_ID,'grid1','ypts_avg')
    Status = nf90_put_att(outID,Fxg_edgV_ID,'grid2','xpts_edgV')

    !'horizontal' cell edges, x-const
    CALL NF_DEF_VAR(Fyg_edgH_ID,outID,(/N_x_ID,N_edgH_ID,N_g_ID,N_t_ID/),'Fyg_edgH_MGQD','Double') !from MGQD eqs
    CALL NF_DEF_UNIT(outID,Fyg_edgH_ID,'ev/cm^3')
    Status = nf90_put_att(outID,Fyg_edgH_ID,'N_grids',3)
    Status = nf90_put_att(outID,Fyg_edgH_ID,'grid0','tpts')
    Status = nf90_put_att(outID,Fyg_edgH_ID,'grid1','ypts_edgH')
    Status = nf90_put_att(outID,Fyg_edgH_ID,'grid2','xpts_avg')

  END IF

  !------------------(high-order)--------------------!
  IF (HO_Fg_out .EQ. 1) THEN
    !'vertical' cell edges, x-const
    CALL NF_DEF_VAR(HO_Fxg_edgV_ID,outID,(/N_edgV_ID,N_y_ID,N_g_ID,N_t_ID/),'Fxg_edgV_HO','Double') !from High-Order eqs
    CALL NF_DEF_UNIT(outID,HO_Fxg_edgV_ID,'ev/cm^3')
    Status = nf90_put_att(outID,HO_Fxg_edgV_ID,'N_grids',3)
    Status = nf90_put_att(outID,HO_Fxg_edgV_ID,'grid0','tpts')
    Status = nf90_put_att(outID,HO_Fxg_edgV_ID,'grid1','ypts_avg')
    Status = nf90_put_att(outID,HO_Fxg_edgV_ID,'grid2','xpts_edgV')

    !'horizontal' cell edges, x-const
    CALL NF_DEF_VAR(HO_Fyg_edgH_ID,outID,(/N_x_ID,N_edgH_ID,N_g_ID,N_t_ID/),'Fyg_edgH_HO','Double') !from High-Order eqs
    CALL NF_DEF_UNIT(outID,HO_Fyg_edgH_ID,'ev/cm^3')
    Status = nf90_put_att(outID,HO_Fyg_edgH_ID,'N_grids',3)
    Status = nf90_put_att(outID,HO_Fyg_edgH_ID,'grid0','tpts')
    Status = nf90_put_att(outID,HO_Fyg_edgH_ID,'grid1','ypts_edgH')
    Status = nf90_put_att(outID,HO_Fyg_edgH_ID,'grid2','xpts_avg')

  END IF

  !===========================================================================!
  !--------------------------------------------------!
  !             Radiation Intensities                !
  !--------------------------------------------------!
  IF (I_out .EQ. 1) THEN
    !cell-averaged
    CALL NF_DEF_VAR(I_avg_ID,outID,(/N_x_ID,N_y_ID,N_m_ID,N_g_ID,N_t_ID/),'I_avg','Double')
    CALL NF_DEF_UNIT(outID,I_avg_ID,'erg/(Ster*cm^3)')
    Status = nf90_put_att(outID,I_avg_ID,'N_grids',3)
    Status = nf90_put_att(outID,I_avg_ID,'grid0','tpts')
    Status = nf90_put_att(outID,I_avg_ID,'grid1','ypts_avg')
    Status = nf90_put_att(outID,I_avg_ID,'grid2','xpts_avg')

    !'vertical' cell edges, x-const
    CALL NF_DEF_VAR(I_edgV_ID,outID,(/N_edgV_ID,N_y_ID,N_m_ID,N_g_ID,N_t_ID/),'I_edgV','Double')
    CALL NF_DEF_UNIT(outID,I_edgV_ID,'erg/(Ster*cm^3)')
    Status = nf90_put_att(outID,I_edgV_ID,'N_grids',3)
    Status = nf90_put_att(outID,I_edgV_ID,'grid0','tpts')
    Status = nf90_put_att(outID,I_edgV_ID,'grid1','ypts_avg')
    Status = nf90_put_att(outID,I_edgV_ID,'grid2','xpts_edgV')

    !'horizontal' cell edges, x-const
    CALL NF_DEF_VAR(I_edgH_ID,outID,(/N_x_ID,N_edgH_ID,N_m_ID,N_g_ID,N_t_ID/),'I_edgH','Double')
    CALL NF_DEF_UNIT(outID,I_edgH_ID,'erg/(Ster*cm^3)')
    Status = nf90_put_att(outID,I_edgH_ID,'N_grids',3)
    Status = nf90_put_att(outID,I_edgH_ID,'grid0','tpts')
    Status = nf90_put_att(outID,I_edgH_ID,'grid1','ypts_edgH')
    Status = nf90_put_att(outID,I_edgH_ID,'grid2','xpts_avg')

  END IF

  !===========================================================================!
  !--------------------------------------------------!
  !                   Opacities                      !
  !--------------------------------------------------!
  IF (kap_out .EQ. 1) THEN
    CALL NF_DEF_VAR(KapE_Bar_ID,outID,(/N_x_ID,N_y_ID,N_t_ID/),'KapE_Bar','Double')
    CALL NF_DEF_UNIT(outID,KapE_Bar_ID,'1/cm')

    CALL NF_DEF_VAR(KapB_ID,outID,(/N_x_ID,N_y_ID,N_g_ID,N_t_ID/),'KapB','Double')
    CALL NF_DEF_UNIT(outID,KapB_ID,'1/cm')

    CALL NF_DEF_VAR(KapE_ID,outID,(/N_x_ID,N_y_ID,N_g_ID,N_t_ID/),'KapE','Double')
    CALL NF_DEF_UNIT(outID,KapB_ID,'1/cm')

    CALL NF_DEF_VAR(KapR_ID,outID,(/N_x_ID,N_y_ID,N_g_ID,N_t_ID/),'KapR','Double')
    CALL NF_DEF_UNIT(outID,KapB_ID,'1/cm')
  END IF

  IF (Src_out .EQ. 1) CALL NF_DEF_VAR(Bg_ID,outID,(/N_x_ID,N_y_ID,N_g_ID,N_t_ID/),'Bg','Double')

  !===========================================================================!
  !--------------------------------------------------!
  !                   QD factors                     !
  !--------------------------------------------------!
  IF (QDfg_out .EQ. 1) THEN
    CALL NF_DEF_VAR(fg_avg_xx_ID,outID,(/N_x_ID,N_y_ID,N_g_ID,N_t_ID/),'fg_avg_xx','Double')
    Status = nf90_put_att(outID,fg_avg_xx_ID,'N_grids',3)
    Status = nf90_put_att(outID,fg_avg_xx_ID,'grid0','tpts')
    Status = nf90_put_att(outID,fg_avg_xx_ID,'grid1','ypts_avg')
    Status = nf90_put_att(outID,fg_avg_xx_ID,'grid2','xpts_avg')

    CALL NF_DEF_VAR(fg_avg_yy_ID,outID,(/N_x_ID,N_y_ID,N_g_ID,N_t_ID/),'fg_avg_yy','Double')
    Status = nf90_put_att(outID,fg_avg_yy_ID,'N_grids',3)
    Status = nf90_put_att(outID,fg_avg_yy_ID,'grid0','tpts')
    Status = nf90_put_att(outID,fg_avg_yy_ID,'grid1','ypts_avg')
    Status = nf90_put_att(outID,fg_avg_yy_ID,'grid2','xpts_avg')

    CALL NF_DEF_VAR(fg_avg_xy_ID,outID,(/N_x_ID,N_y_ID,N_g_ID,N_t_ID/),'fg_avg_xy','Double')
    Status = nf90_put_att(outID,fg_avg_xy_ID,'N_grids',3)
    Status = nf90_put_att(outID,fg_avg_xy_ID,'grid0','tpts')
    Status = nf90_put_att(outID,fg_avg_xy_ID,'grid1','ypts_avg')
    Status = nf90_put_att(outID,fg_avg_xy_ID,'grid2','xpts_avg')

    CALL NF_DEF_VAR(fg_edgV_xx_ID,outID,(/N_edgV_ID,N_y_ID,N_g_ID,N_t_ID/),'fg_edgV_xx','Double')
    Status = nf90_put_att(outID,fg_edgV_xx_ID,'N_grids',3)
    Status = nf90_put_att(outID,fg_edgV_xx_ID,'grid0','tpts')
    Status = nf90_put_att(outID,fg_edgV_xx_ID,'grid1','ypts_avg')
    Status = nf90_put_att(outID,fg_edgV_xx_ID,'grid2','xpts_edgV')

    CALL NF_DEF_VAR(fg_edgV_xy_ID,outID,(/N_edgV_ID,N_y_ID,N_g_ID,N_t_ID/),'fg_edgV_xy','Double')
    Status = nf90_put_att(outID,fg_edgV_xy_ID,'N_grids',3)
    Status = nf90_put_att(outID,fg_edgV_xy_ID,'grid0','tpts')
    Status = nf90_put_att(outID,fg_edgV_xy_ID,'grid1','ypts_avg')
    Status = nf90_put_att(outID,fg_edgV_xy_ID,'grid2','xpts_edgV')

    CALL NF_DEF_VAR(fg_edgH_yy_ID,outID,(/N_x_ID,N_edgH_ID,N_g_ID,N_t_ID/),'fg_edgH_yy','Double')
    Status = nf90_put_att(outID,fg_edgH_yy_ID,'N_grids',3)
    Status = nf90_put_att(outID,fg_edgH_yy_ID,'grid0','tpts')
    Status = nf90_put_att(outID,fg_edgH_yy_ID,'grid1','ypts_edgH')
    Status = nf90_put_att(outID,fg_edgH_yy_ID,'grid2','xpts_avg')

    CALL NF_DEF_VAR(fg_edgH_xy_ID,outID,(/N_x_ID,N_edgH_ID,N_g_ID,N_t_ID/),'fg_edgH_xy','Double')
    Status = nf90_put_att(outID,fg_edgH_xy_ID,'N_grids',3)
    Status = nf90_put_att(outID,fg_edgH_xy_ID,'grid0','tpts')
    Status = nf90_put_att(outID,fg_edgH_xy_ID,'grid1','ypts_edgH')
    Status = nf90_put_att(outID,fg_edgH_xy_ID,'grid2','xpts_avg')

  END IF

  IF (D_out .EQ. 1) THEN
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
  END IF

  IF (old_parms_out .EQ. 1) THEN
    CALL NF_DEF_VAR(G_old_ID,outID,(/N_x_ID,N_y_ID,N_g_ID,N_t_ID/),'Gold_L','Double')
    CALL NF_DEF_VAR(Pold_L_ID,outID,(/N_x_ID,N_y_ID,N_g_ID,N_t_ID/),'Pold_L','Double')
    CALL NF_DEF_VAR(Pold_B_ID,outID,(/N_x_ID,N_y_ID,N_g_ID,N_t_ID/),'Pold_B','Double')
    CALL NF_DEF_VAR(Pold_R_ID,outID,(/N_x_ID,N_y_ID,N_g_ID,N_t_ID/),'Pold_R','Double')
    CALL NF_DEF_VAR(Pold_T_ID,outID,(/N_x_ID,N_y_ID,N_g_ID,N_t_ID/),'Pold_T','Double')

    CALL NF_DEF_VAR(Gold_hat_ID,outID,(/N_x_ID,N_y_ID,N_t_ID/),'Gold_hat','Double')
    CALL NF_DEF_VAR(Rhat_old_ID,outID,(/N_x_ID,N_y_ID,N_t_ID/),'Rhat_old','Double')
    CALL NF_DEF_VAR(PL_ID,outID,(/N_x_ID,N_y_ID,N_t_ID/),'PL','Double')
    CALL NF_DEF_VAR(PB_ID,outID,(/N_x_ID,N_y_ID,N_t_ID/),'PB','Double')
    CALL NF_DEF_VAR(PR_ID,outID,(/N_x_ID,N_y_ID,N_t_ID/),'PR','Double')
    CALL NF_DEF_VAR(PT_ID,outID,(/N_x_ID,N_y_ID,N_t_ID/),'PT','Double')
  END IF

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

  IF (fg_pod_out .EQ. 1) THEN
    CALL NF_DEF_VAR(PODgsum_ID,outID,(/N_dsets_ID/),'PODgsum','Int')
    CALL NF_DEF_VAR(PODerr_ID,outID,(/N_dsets_ID/),'PODerr','Double')
    CALL NF_DEF_VAR(rrank_BCg_ID,outID,(/N_g_ID,N_dsets_ID/),'rrank_BCg','Int')
    CALL NF_DEF_VAR(rrank_fg_avg_xx_ID,outID,(/N_g_ID,N_dsets_ID/),'rrank_fg_avg_xx','Int')
    CALL NF_DEF_VAR(rrank_fg_edgV_xx_ID,outID,(/N_g_ID,N_dsets_ID/),'rrank_fg_edgV_xx','Int')
    CALL NF_DEF_VAR(rrank_fg_avg_yy_ID,outID,(/N_g_ID,N_dsets_ID/),'rrank_fg_avg_yy','Int')
    CALL NF_DEF_VAR(rrank_fg_edgH_yy_ID,outID,(/N_g_ID,N_dsets_ID/),'rrank_fg_edgH_yy','Int')
    CALL NF_DEF_VAR(rrank_fg_edgV_xy_ID,outID,(/N_g_ID,N_dsets_ID/),'rrank_fg_edgV_xy','Int')
    CALL NF_DEF_VAR(rrank_fg_edgH_xy_ID,outID,(/N_g_ID,N_dsets_ID/),'rrank_fg_edgH_xy','Int')
  ELSE IF (fg_dmd_out .EQ. 1) THEN
    CALL NF_DEF_VAR(DMDgsum_ID,outID,(/N_dsets_ID/),'DMDgsum','Int')
    CALL NF_DEF_VAR(rrank_BCg_ID,outID,(/N_g_ID,N_dsets_ID/),'rrank_BCg','Int')
    CALL NF_DEF_VAR(rrank_fg_avg_xx_ID,outID,(/N_g_ID,N_dsets_ID/),'rrank_fg_avg_xx','Int')
    CALL NF_DEF_VAR(rrank_fg_edgV_xx_ID,outID,(/N_g_ID,N_dsets_ID/),'rrank_fg_edgV_xx','Int')
    CALL NF_DEF_VAR(rrank_fg_avg_yy_ID,outID,(/N_g_ID,N_dsets_ID/),'rrank_fg_avg_yy','Int')
    CALL NF_DEF_VAR(rrank_fg_edgH_yy_ID,outID,(/N_g_ID,N_dsets_ID/),'rrank_fg_edgH_yy','Int')
    CALL NF_DEF_VAR(rrank_fg_edgV_xy_ID,outID,(/N_g_ID,N_dsets_ID/),'rrank_fg_edgV_xy','Int')
    CALL NF_DEF_VAR(rrank_fg_edgH_xy_ID,outID,(/N_g_ID,N_dsets_ID/),'rrank_fg_edgH_xy','Int')
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
SUBROUTINE TIMESTEP_OUTS(outID,t,out_freq,I_out,HO_Eg_out,HO_Fg_out,HO_E_out,HO_F_out,Eg_out,Fg_out,MGQD_E_out,MGQD_F_out,&
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
  DL_xx,DR_xx,DC_yy,DB_yy,DT_yy,DL_xy,DB_xy,DR_xy,DT_xy,G_old, Pold_L,Pold_B,Pold_R,Pold_T,Gold_hat,Rhat_old,PL,PB,PR,PT)

  INTEGER,INTENT(IN):: outID, t
  INTEGER,INTENT(IN):: out_freq, I_out, HO_Eg_out, HO_Fg_out, HO_E_out, HO_F_out
  INTEGER,INTENT(IN):: Eg_out, Fg_out, MGQD_E_out, MGQD_F_out, QDfg_out
  INTEGER,INTENT(IN):: E_out, F_out, D_out
  INTEGER,INTENT(IN):: old_parms_out, its_out, conv_out, kap_out, Src_out, use_grey
  INTEGER,INTENT(IN):: Temp_ID, E_avg_ID, E_edgV_ID, E_edgH_ID, MGQD_E_avg_ID, MGQD_E_edgV_ID, MGQD_E_edgH_ID, HO_E_avg_ID
  INTEGER,INTENT(IN):: HO_E_edgV_ID, HO_E_edgH_ID, Fx_edgV_ID, Fy_edgH_ID, MGQD_Fx_edgV_ID, MGQD_Fy_edgH_ID, HO_Fx_edgV_ID
  INTEGER,INTENT(IN):: HO_Fy_edgH_ID, Eg_avg_ID, Eg_edgV_ID, Eg_edgH_ID, HO_Eg_avg_ID, HO_Eg_edgV_ID, HO_Eg_edgH_ID
  INTEGER,INTENT(IN):: Fxg_edgV_ID, Fyg_edgH_ID, HO_Fxg_edgV_ID, HO_Fyg_edgH_ID, I_avg_ID, I_edgV_ID, I_edgH_ID
  INTEGER,INTENT(IN):: KapE_Bar_ID, KapB_ID, KapE_ID, KapR_ID, Bg_ID
  INTEGER,INTENT(IN):: Cg_L_ID, Cg_B_ID, Cg_R_ID, Cg_T_ID, Eg_in_L_ID, Eg_in_B_ID, Eg_in_R_ID, Eg_in_T_ID
  INTEGER,INTENT(IN):: Fg_in_L_ID, Fg_in_B_ID, Fg_in_R_ID, Fg_in_T_ID
  INTEGER,INTENT(IN):: Cb_L_ID, Cb_B_ID, Cb_R_ID, Cb_T_ID, E_in_L_ID, E_in_B_ID, E_in_R_ID, E_in_T_ID
  INTEGER,INTENT(IN):: F_in_L_ID, F_in_B_ID, F_in_R_ID, F_in_T_ID
  INTEGER,INTENT(IN):: fg_avg_xx_ID, fg_avg_yy_ID, fg_avg_xy_ID, fg_edgV_xx_ID, fg_edgV_xy_ID, fg_edgH_yy_ID, fg_edgH_xy_ID
  INTEGER,INTENT(IN):: DC_xx_ID, DL_xx_ID, DR_xx_ID, DC_yy_ID, DB_yy_ID, DT_yy_ID, DL_xy_ID, DB_xy_ID, DR_xy_ID, DT_xy_ID
  INTEGER,INTENT(IN):: G_old_ID, Pold_L_ID, Pold_B_ID, Pold_R_ID, Pold_T_ID
  INTEGER,INTENT(IN):: Gold_hat_ID, Rhat_old_ID, PL_ID, PB_ID, PR_ID, PT_ID

  REAL*8,INTENT(IN):: Temp(:,:), E_avg(:,:), E_edgV(:,:), E_edgH(:,:), MGQD_E_avg(:,:), MGQD_E_edgV(:,:), MGQD_E_edgH(:,:)
  REAL*8,INTENT(IN):: HO_E_avg(:,:), HO_E_edgV(:,:), HO_E_edgH(:,:), Fx_edgV(:,:), Fy_edgH(:,:), MGQD_Fx_edgV(:,:)
  REAL*8,INTENT(IN):: MGQD_Fy_edgH(:,:), HO_Fx_edgV(:,:), HO_Fy_edgH(:,:), Eg_avg(:,:,:), Eg_edgV(:,:,:), Eg_edgH(:,:,:)
  REAL*8,INTENT(IN):: HO_Eg_avg(:,:,:), HO_Eg_edgV(:,:,:), HO_Eg_edgH(:,:,:), Fxg_edgV(:,:,:), Fyg_edgH(:,:,:), HO_Fxg_edgV(:,:,:)
  REAL*8,INTENT(IN):: HO_Fyg_edgH(:,:,:), I_avg(:,:,:,:), I_edgV(:,:,:,:), I_edgH(:,:,:,:), KapE_Bar(:,:), KapB(:,:,:)
  REAL*8,INTENT(IN):: KapE(:,:,:), KapR(:,:,:), Bg(:,:,:), Cg_L(:,:), Cg_B(:,:), Cg_R(:,:), Cg_T(:,:), Eg_in_L(:,:), Eg_in_B(:,:)
  REAL*8,INTENT(IN):: Eg_in_R(:,:), Eg_in_T(:,:), Fg_in_L(:,:), Fg_in_B(:,:), Fg_in_R(:,:), Fg_in_T(:,:)
  REAL*8,INTENT(IN):: Cb_L(:), Cb_B(:), Cb_R(:), Cb_T(:), E_in_L(:), E_in_B(:), E_in_R(:), E_in_T(:)
  REAL*8,INTENT(IN):: F_in_L(:), F_in_B(:), F_in_R(:), F_in_T(:)
  REAL*8,INTENT(IN):: fg_avg_xx(:,:,:), fg_avg_yy(:,:,:), fg_avg_xy(:,:,:), fg_edgV_xx(:,:,:), fg_edgV_xy(:,:,:), fg_edgH_yy(:,:,:)
  REAL*8,INTENT(IN):: fg_edgH_xy(:,:,:), DC_xx(:,:), DL_xx(:,:), DR_xx(:,:), DC_yy(:,:), DB_yy(:,:), DT_yy(:,:), DL_xy(:,:)
  REAL*8,INTENT(IN):: DB_xy(:,:), DR_xy(:,:), DT_xy(:,:), G_old(:,:,:), Pold_L(:,:,:), Pold_B(:,:,:), Pold_R(:,:,:), Pold_T(:,:,:)
  REAL*8,INTENT(IN):: Gold_hat(:,:), Rhat_old(:,:), PL(:,:), PB(:,:), PR(:,:), PT(:,:)

  CALL NF_PUT_t_VAR(outID,Temp_ID,Temp,t)

  IF (use_grey .EQ. 1) THEN
    IF (E_out .EQ. 1) THEN
      CALL NF_PUT_t_VAR(outID,E_avg_ID,E_avg,t)
      CALL NF_PUT_t_VAR(outID,E_edgV_ID,E_edgV,t)
      CALL NF_PUT_t_VAR(outID,E_edgH_ID,E_edgH,t)
    END IF

    IF (F_out .EQ. 1) THEN
      CALL NF_PUT_t_VAR(outID,Fx_edgV_ID,Fx_edgV,t)
      CALL NF_PUT_t_VAR(outID,Fy_edgH_ID,Fy_edgH,t)
    END IF

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

    IF (D_out .EQ. 1) THEN
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
    END IF

    IF (old_parms_out .EQ. 1) THEN
      CALL NF_PUT_t_VAR(outID,Gold_hat_ID,Gold_hat,t)
      CALL NF_PUT_t_VAR(outID,Rhat_old_ID,Rhat_old,t)
      CALL NF_PUT_t_VAR(outID,PL_ID,PL,t)
      CALL NF_PUT_t_VAR(outID,PB_ID,PB,t)
      CALL NF_PUT_t_VAR(outID,PR_ID,PR,t)
      CALL NF_PUT_t_VAR(outID,PT_ID,PT,t)
    END IF
  END IF

  IF (MGQD_E_out .EQ. 1) THEN
    CALL NF_PUT_t_VAR(outID,MGQD_E_avg_ID,MGQD_E_avg,t)
    CALL NF_PUT_t_VAR(outID,MGQD_E_edgV_ID,MGQD_E_edgV,t)
    CALL NF_PUT_t_VAR(outID,MGQD_E_edgH_ID,MGQD_E_edgH,t)
  END IF

  IF (MGQD_F_out .EQ. 1) THEN
    CALL NF_PUT_t_VAR(outID,MGQD_Fx_edgV_ID,MGQD_Fx_edgV,t)
    CALL NF_PUT_t_VAR(outID,MGQD_Fy_edgH_ID,MGQD_Fy_edgH,t)
  END IF

  IF (Eg_out .EQ. 1) THEN
    CALL NF_PUT_t_VAR(outID,Eg_avg_ID,Eg_avg,t)
    CALL NF_PUT_t_VAR(outID,Eg_edgV_ID,Eg_edgV,t)
    CALL NF_PUT_t_VAR(outID,Eg_edgH_ID,Eg_edgH,t)
  END IF

  IF (Fg_out .EQ. 1) THEN
    CALL NF_PUT_t_VAR(outID,Fxg_edgV_ID,Fxg_edgV,t)
    CALL NF_PUT_t_VAR(outID,Fyg_edgH_ID,Fyg_edgH,t)
  END IF

  IF (HO_E_out .EQ. 1) THEN
    CALL NF_PUT_t_VAR(outID,HO_E_avg_ID,HO_E_avg,t)
    CALL NF_PUT_t_VAR(outID,HO_E_edgV_ID,HO_E_edgV,t)
    CALL NF_PUT_t_VAR(outID,HO_E_edgH_ID,HO_E_edgH,t)
  END IF

  IF (HO_F_out .EQ. 1) THEN
    CALL NF_PUT_t_VAR(outID,HO_Fx_edgV_ID,HO_Fx_edgV,t)
    CALL NF_PUT_t_VAR(outID,HO_Fy_edgH_ID,HO_Fy_edgH,t)
  END IF

  IF (HO_Eg_out .EQ. 1) THEN
    CALL NF_PUT_t_VAR(outID,HO_Eg_avg_ID,HO_Eg_avg,t)
    CALL NF_PUT_t_VAR(outID,HO_Eg_edgV_ID,HO_Eg_edgV,t)
    CALL NF_PUT_t_VAR(outID,HO_Eg_edgH_ID,HO_Eg_edgH,t)
  END IF

  IF (HO_Fg_out .EQ. 1) THEN
    CALL NF_PUT_t_VAR(outID,HO_Fxg_edgV_ID,HO_Fxg_edgV,t)
    CALL NF_PUT_t_VAR(outID,HO_Fyg_edgH_ID,HO_Fyg_edgH,t)
  END IF

  IF (I_out .EQ. 1) THEN
    CALL NF_PUT_t_VAR(outID,I_avg_ID,I_avg,t)
    CALL NF_PUT_t_VAR(outID,I_edgV_ID,I_edgV,t)
    CALL NF_PUT_t_VAR(outID,I_edgH_ID,I_edgH,t)
  END IF

  IF (QDfg_out .EQ. 1) THEN
    CALL NF_PUT_t_VAR(outID,fg_avg_xx_ID,fg_avg_xx,t)
    CALL NF_PUT_t_VAR(outID,fg_avg_yy_ID,fg_avg_yy,t)
    CALL NF_PUT_t_VAR(outID,fg_avg_xy_ID,fg_avg_xy,t)
    CALL NF_PUT_t_VAR(outID,fg_edgV_xx_ID,fg_edgV_xx,t)
    CALL NF_PUT_t_VAR(outID,fg_edgV_xy_ID,fg_edgV_xy,t)
    CALL NF_PUT_t_VAR(outID,fg_edgH_yy_ID,fg_edgH_yy,t)
    CALL NF_PUT_t_VAR(outID,fg_edgH_xy_ID,fg_edgH_xy,t)
  END IF

  IF (kap_out .EQ. 1) THEN
    CALL NF_PUT_t_VAR(outID,KapE_Bar_ID,KapE_Bar,t)
    CALL NF_PUT_t_VAR(outID,KapB_ID,KapB,t)
    CALL NF_PUT_t_VAR(outID,KapE_ID,KapE,t)
    CALL NF_PUT_t_VAR(outID,KapR_ID,KapR,t)
  END IF

  IF (Src_out .EQ. 1) CALL NF_PUT_t_VAR(outID,Bg_ID,Bg,t)

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

  IF (old_parms_out .EQ. 1) THEN
    CALL NF_PUT_t_VAR(outID,G_old_ID,G_old,t)
    CALL NF_PUT_t_VAR(outID,Pold_L_ID,Pold_L,t)
    CALL NF_PUT_t_VAR(outID,Pold_B_ID,Pold_B,t)
    CALL NF_PUT_t_VAR(outID,Pold_R_ID,Pold_R,t)
    CALL NF_PUT_t_VAR(outID,Pold_T_ID,Pold_T,t)
  END IF

END SUBROUTINE

!==================================================================================================================================!
!
!==================================================================================================================================!
SUBROUTINE RESTART_OUT(Time, Temp, I_crn, fg_avg_xx, fg_avg_yy, fg_edgV_xx, fg_edgV_xy, fg_edgH_yy, fg_edgH_xy, RT_Src,&
  MGQD_Src,GQD_Src, E_avg, Eg_avg, Eg_edgV, Eg_edgH, Fxg_edgV, Fyg_edgH, KapE_Bar, KapE, KapR, Cg_L, Cg_B, Cg_R, Cg_T,&
  restart_outfile, resf_unit)

  REAL*8,INTENT(IN):: Time, Temp(:,:), I_crn(:,:,:,:), RT_Src(:,:,:,:)
  REAL*8,INTENT(IN):: fg_avg_xx(:,:,:), fg_avg_yy(:,:,:), fg_edgV_xx(:,:,:)
  REAL*8,INTENT(IN):: fg_edgV_xy(:,:,:), fg_edgH_yy(:,:,:), fg_edgH_xy(:,:,:)
  REAL*8,INTENT(IN):: MGQD_Src(:,:,:), GQD_Src(:,:), E_avg(:,:), Eg_avg(:,:,:)
  REAL*8,INTENT(IN):: Eg_edgV(:,:,:), Eg_edgH(:,:,:), Fxg_edgV(:,:,:)
  REAL*8,INTENT(IN):: Fyg_edgH(:,:,:), KapE_Bar(:,:), KapE(:,:,:), KapR(:,:,:)
  REAL*8,INTENT(IN):: Cg_L(:,:), Cg_B(:,:), Cg_R(:,:), Cg_T(:,:)
  CHARACTER(100),INTENT(IN):: restart_outfile
  INTEGER,INTENT(IN):: resf_unit
  INTEGER:: err

  OPEN(UNIT=resf_unit,FILE=restart_outfile,FORM='UNFORMATTED',STATUS='REPLACE',ACTION='WRITE',IOSTAT=err)
  !   making sure file exists/opens, if not tells user
  IF (err .NE. 0) THEN
      WRITE(*,'(3A)') 'The file, ',restart_outfile,', could not open properly.'
      STOP
  END IF

  WRITE(resf_unit) Time
  WRITE(resf_unit) Temp
  WRITE(resf_unit) I_crn
  WRITE(resf_unit) fg_avg_xx
  WRITE(resf_unit) fg_avg_yy
  WRITE(resf_unit) fg_edgV_xx
  WRITE(resf_unit) fg_edgV_xy
  WRITE(resf_unit) fg_edgH_yy
  WRITE(resf_unit) fg_edgH_xy
  WRITE(resf_unit) RT_Src
  WRITE(resf_unit) MGQD_Src
  WRITE(resf_unit) GQD_Src
  WRITE(resf_unit) E_avg
  WRITE(resf_unit) Eg_avg
  WRITE(resf_unit) Eg_edgV
  WRITE(resf_unit) Eg_edgH
  WRITE(resf_unit) Fxg_edgV
  WRITE(resf_unit) Fyg_edgH
  WRITE(resf_unit) KapE_Bar
  WRITE(resf_unit) KapE
  WRITE(resf_unit) KapR
  WRITE(resf_unit) Cg_L
  WRITE(resf_unit) Cg_B
  WRITE(resf_unit) Cg_R
  WRITE(resf_unit) Cg_T

  CLOSE ( resf_unit, STATUS='KEEP')

END SUBROUTINE RESTART_OUT

!==================================================================================================================================!
!
!==================================================================================================================================!

END MODULE OUTPUTS
