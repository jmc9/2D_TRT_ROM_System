!==================================================================================================================================!
!
!==================================================================================================================================!
program main

    USE INPUTS
    USE OUTPUTS
    USE INITIALIZERS
    USE TRANSPORT_SOLVES
    USE LA_TOOLS
    USE ALGORITHMS
    ! USE NCDF_IO

    IMPLICIT NONE

    !--------------------------------------------------!
    !               PHYSICAL CONSTANTS                 !
    !--------------------------------------------------!
    REAL*8,PARAMETER:: h=6.62613d-19    !erg*sh
    REAL*8,PARAMETER:: c=2.99792458d2   !cm/sh
    REAL*8,PARAMETER:: pi=4d0*ATAN(1d0)
    REAL*8,PARAMETER:: erg=6.24150934d11  !eV/erg -- 1 erg is this many ev's
    REAL*8,PARAMETER:: Kap0=27d0*1d9   !(ev**3/cm)
    REAL*8:: comp_unit   !erg/unit
    REAL*8:: sig_R !(erg/(ev**4 cm**2 sh))
    REAL*8:: aR !(erg/(ev**4 cm**3))
    REAL*8:: Cv !=0.5917d0*ar*(T_rad)**3 (erg/(ev cm**3))

    !--------------------------------------------------!
    !                   CODE INPUTS                    !
    !--------------------------------------------------!
    REAL*8:: chi, conv_ho, conv_lo, conv_gr, line_src, Nwt_upbnd
    REAL*8:: xlen, ylen, Tlen, delt, bcT_left, bcT_right, bcT_top, bcT_bottom, Tini
    REAL*8,ALLOCATABLE:: Delx(:), Dely(:)
    REAL*8,ALLOCATABLE:: out_times(:)
    REAL*8,ALLOCATABLE:: nu_g(:)
    INTEGER:: database_gen, database_add, use_grey
    INTEGER:: maxit_RTE, maxit_MLOQD, maxit_GLOQD, conv_type, N_m, threads
    INTEGER:: N_x, N_y, N_t, N_g
    INTEGER:: n_out_times, restart_outlen, TEMP_out, GREY_E_out, GREY_F_out, HO_E_out
    INTEGER:: GREY_kap_out, GREY_fsmall_out, MG_fsmall_out, res_history_out
    INTEGER:: BC_Type(4)
    INTEGER,ALLOCATABLE:: out_time_steps(:)
    CHARACTER(100):: outfile, restart_outfile, decomp_outfile, TEMP_outfile, GREY_E_outfile, HO_E_outfile
    CHARACTER(100):: GREY_F_outfile, GREY_kap_outfile, GREY_fsmall_outfile, MG_fsmall_outfile, res_history_outfile
    CHARACTER(100):: run_type, restart_infile
    CHARACTER(100):: kapE_dT_flag, enrgy_strc, quadrature
    REAL*8,ALLOCATABLE:: Omega_x(:), Omega_y(:), quad_weight(:)
    REAL*8:: Theta=1d0
    LOGICAL:: Res_Calc = .TRUE.

    !--------------------------------------------------!
    !                                                  !
    !--------------------------------------------------!
    REAL*8,ALLOCATABLE:: I_avg(:,:,:,:), I_avg_old(:,:,:,:), I_edgV(:,:,:,:), I_edgH(:,:,:,:)
    REAL*8,ALLOCATABLE:: I_crn(:,:,:,:), I_crn_old(:,:,:,:), Ic_edgV(:,:,:,:), Ic_edgH(:,:,:,:)
    REAL*8,ALLOCATABLE:: RT_Src(:,:,:,:)
    REAL*8,ALLOCATABLE:: Hg_avg_xx(:,:,:), Hg_avg_xy(:,:,:), Hg_avg_yy(:,:,:)
    REAL*8,ALLOCATABLE:: Hg_edgV_xx(:,:,:), Hg_edgV_xy(:,:,:)
    REAL*8,ALLOCATABLE:: Hg_edgH_yy(:,:,:), Hg_edgH_xy(:,:,:)
    REAL*8,ALLOCATABLE:: HO_Eg_avg(:,:,:), HO_Eg_edgV(:,:,:), HO_Eg_edgH(:,:,:)
    REAL*8,ALLOCATABLE:: HO_Fxg_edgV(:,:,:), HO_Fyg_edgH(:,:,:)
    REAL*8,ALLOCATABLE:: HO_E_avg(:,:), HO_E_edgV(:,:), HO_E_edgH(:,:)
    REAL*8,ALLOCATABLE:: HO_Fx_edgV(:,:), HO_Fy_edgH(:,:)

    !--------------------------------------------------!
    !                                                  !
    !--------------------------------------------------!
    REAL*8,ALLOCATABLE:: Eg_avg(:,:,:), Eg_edgV(:,:,:), Eg_edgH(:,:,:)
    REAL*8,ALLOCATABLE:: Fxg_edgV(:,:,:), Fyg_edgH(:,:,:)
    REAL*8,ALLOCATABLE:: fg_avg_xx(:,:,:), fg_avg_xy(:,:,:), fg_avg_yy(:,:,:)
    REAL*8,ALLOCATABLE:: fg_edgV_xx(:,:,:), fg_edgV_xy(:,:,:)
    REAL*8,ALLOCATABLE:: fg_edgH_yy(:,:,:), fg_edgH_xy(:,:,:)
    REAL*8,ALLOCATABLE:: fg_avg_xx_old(:,:,:), fg_avg_xy_old(:,:,:), fg_avg_yy_old(:,:,:)
    REAL*8,ALLOCATABLE:: fg_edgV_xx_old(:,:,:), fg_edgV_xy_old(:,:,:)
    REAL*8,ALLOCATABLE:: fg_edgH_yy_old(:,:,:), fg_edgH_xy_old(:,:,:)
    REAL*8,ALLOCATABLE:: Cg_L(:,:), Cg_B(:,:), Cg_R(:,:), Cg_T(:,:)
    REAL*8,ALLOCATABLE:: Eg_in_L(:,:), Eg_in_B(:,:), Eg_in_R(:,:), Eg_in_T(:,:)
    REAL*8,ALLOCATABLE:: Fg_in_L(:,:), Fg_in_B(:,:), Fg_in_R(:,:), Fg_in_T(:,:)
    REAL*8,ALLOCATABLE:: MGQD_Src(:,:,:), MGQD_Src_old(:,:,:)
    REAL*8,ALLOCATABLE:: Eg_avg_old(:,:,:), Eg_edgV_old(:,:,:), Eg_edgH_old(:,:,:)
    REAL*8,ALLOCATABLE:: Fxg_edgV_old(:,:,:), Fyg_edgH_old(:,:,:)
    REAL*8,ALLOCATABLE:: MGQD_Residual(:,:,:,:,:,:,:)
    REAL*8,ALLOCATABLE:: MGQD_E_avg(:,:), MGQD_E_edgV(:,:), MGQD_E_edgH(:,:)
    REAL*8,ALLOCATABLE:: MGQD_Fx_edgV(:,:), MGQD_Fy_edgH(:,:)
    REAL*8,ALLOCATABLE:: G_old(:,:,:), Pold_L(:,:,:), Pold_B(:,:,:), Pold_R(:,:,:), Pold_T(:,:,:)

    !--------------------------------------------------!
    !                                                  !
    !--------------------------------------------------!
    !!!!!!!!!!
    !NONE OF THESE ARE ALLOCATED YET
    REAL*8,ALLOCATABLE:: E_avg(:,:), E_edgV(:,:), E_edgH(:,:)
    REAL*8,ALLOCATABLE:: Fx_edgV(:,:), Fy_edgH(:,:)
    !!!!!!!!!!

    !--------------------------------------------------!
    !                                                  !
    !--------------------------------------------------!
    INTEGER:: outID
    INTEGER:: N_x_ID, N_y_ID, N_m_ID, N_g_ID, N_t_ID, N_edgV_ID, N_edgH_ID
    INTEGER:: RT_Its_ID, MGQD_Its_ID, GQD_Its_ID, Norm_Types_ID, Boundaries_ID
    INTEGER:: c_ID, h_ID, pi_ID, erg_ID, Comp_Unit_ID, cv_ID
    INTEGER:: Temp_ID, E_avg_ID, E_edgV_ID, E_edgH_ID, MGQD_E_avg_ID, MGQD_E_edgV_ID, MGQD_E_edgH_ID, HO_E_avg_ID
    INTEGER:: HO_E_edgV_ID, HO_E_edgH_ID, Fx_edgV_ID, Fy_edgH_ID, MGQD_Fx_edgV_ID, MGQD_Fy_edgH_ID, HO_Fx_edgV_ID
    INTEGER:: HO_Fy_edgH_ID, Eg_avg_ID, Eg_edgV_ID, Eg_edgH_ID, HO_Eg_avg_ID, HO_Eg_edgV_ID, HO_Eg_edgH_ID
    INTEGER:: Fxg_edgV_ID, Fyg_edgH_ID, HO_Fxg_edgV_ID, HO_Fyg_edgH_ID, I_avg_ID, I_edgV_ID, I_edgH_ID
    INTEGER:: RT_Residual_ID, RT_ResNorms_ID, MGQD_Residual_ID, MGQD_ResNorms_ID
    INTEGER:: Del_T_ID, Del_E_avg_ID, Del_E_edgV_ID, Del_E_edgH_ID, Del_Fx_edgV_ID, Del_Fy_edgH_ID, Del_T_Norms_ID
    INTEGER:: Del_E_avg_Norms_ID, Del_E_edgV_Norms_ID, Del_E_edgH_Norms_ID, Del_Fx_edgV_Norms_ID, Del_Fy_edgH_Norms_ID
    INTEGER:: RT_ItCount_ID, MGQD_ItCount_ID, GQD_ItCount_ID

    !--------------------------------------------------!
    !                                                  !
    !--------------------------------------------------!
    REAL*8,ALLOCATABLE:: KapE(:,:,:), KapB(:,:,:), KapR(:,:,:), Bg(:,:,:), Temp(:,:)
    REAL*8,ALLOCATABLE:: KapE_old(:,:,:), KapR_old(:,:,:)
    REAL*8,ALLOCATABLE:: Temp_Old(:,:)
    REAL*8,ALLOCATABLE:: Temp_Times(:,:,:)
    REAL*8,ALLOCATABLE:: HO_E_avg_Times(:,:,:), GREY_E_avg_Times(:,:,:), MGQD_E_avg_Times(:,:,:)
    REAL*8,ALLOCATABLE:: RT_Residual(:,:,:,:,:,:)
    REAL*8,ALLOCATABLE:: A(:,:), Times(:)

    !--------------------------------------------------!
    !                 MISCELLANIOUS                    !
    !--------------------------------------------------!
    INTEGER:: datevalues(8), err
    INTEGER:: solve_time1, solve_time2, clockrate
    REAL*8:: solve_time, Start_Time

    CALL DATE_AND_TIME(VALUES=datevalues)

    WRITE(*,*) ' ___  _____    _______ _____ _______   _____   ____  __  __                 _      '
    WRITE(*,*) '|__ \|  __ \  |__   __|  __ \__   __| |  __ \ / __ \|  \/  |               | |     '
    WRITE(*,*) '   ) | |  | |    | |  | |__) | | |    | |__) | |  | | \  / |   ___ ___   __| | ___ '
    WRITE(*,*) '  / /| |  | |    | |  |  _  /  | |    |  _  /| |  | | |\/| |  / __/ _ \ / _` |/ _ \'
    WRITE(*,*) ' / /_| |__| |    | |  | | \ \  | |    | | \ \| |__| | |  | | | (_| (_) | (_| |  __/'
    WRITE(*,*) '|____|_____/     |_|  |_|  \_\ |_|    |_|  \_\\____/|_|  |_|  \___\___/ \__,_|\___|'
    WRITE(*,*) 'Joseph M. Coale - NCSU'

    CALL SYSTEM_CLOCK(solve_time1)

    CALL INPUT(database_gen,database_add,run_type,restart_infile,use_grey,chi,conv_ho,conv_lo,conv_gr,&
      comp_unit,line_src,maxit_RTE,maxit_MLOQD,maxit_GLOQD,conv_type,N_m,threads,kapE_dT_flag,enrgy_strc,&
      Nwt_upbnd,erg,xlen,ylen,N_x,N_y,tlen,delt,bcT_left,bcT_right,bcT_top,bcT_bottom,Tini,sig_R,ar,&
      pi,c,h,delx,dely,cv,out_times,out_time_steps,n_out_times,restart_outlen,TEMP_out,GREY_E_out,GREY_F_out,&
      GREY_kap_out,GREY_fsmall_out,MG_fsmall_out,res_history_out,outfile,restart_outfile,decomp_outfile,TEMP_outfile,&
      GREY_E_outfile,GREY_F_outfile,GREY_kap_outfile,GREY_fsmall_outfile,MG_fsmall_outfile,res_history_outfile,&
      HO_E_out,HO_E_outfile,nu_g,N_g,Omega_x,Omega_y,quad_weight,N_t,quadrature,BC_Type)

    CALL OUTFILE_INIT(outID,N_x_ID,N_y_ID,N_m_ID,N_g_ID,N_t_ID,N_edgV_ID,N_edgH_ID,RT_Its_ID,MGQD_Its_ID,GQD_Its_ID,&
      Norm_Types_ID,Boundaries_ID,c_ID,h_ID,pi_ID,erg_ID,Comp_Unit_ID,cv_ID,c,h,pi,erg,Comp_Unit,cv,chi,conv_ho,conv_lo,&
      conv_gr,line_src,xlen,ylen,Delx,Dely,tlen,Delt,bcT_left,bcT_bottom,bcT_right,bcT_top,Tini,N_x,N_y,N_m,N_g,N_t,&
      database_gen,use_grey,maxit_RTE,maxit_MLOQD,maxit_GLOQD,conv_type,threads,BC_type,outfile,run_type,kapE_dT_flag,&
      quadrature,enrgy_strc)

    ! Start_Time = 0d0
    ! CALL MISC_INIT(Delx,Dely,Delt,Start_Time,N_t,A,Times)
    ! CALL RT_INIT(I_avg,I_avg_old,I_edgV,I_edgH,I_crn,I_crn_old,Ic_edgV,Ic_edgH,Hg_avg_xx,Hg_avg_xy,Hg_avg_yy,&
    !   Hg_edgV_xx,Hg_edgV_xy,Hg_edgH_yy,Hg_edgH_xy,HO_Eg_avg,HO_Eg_edgV,HO_Eg_edgH,HO_Fxg_edgV,HO_Fyg_edgH,HO_E_avg,&
    !   HO_E_edgV,HO_E_edgH,HO_Fx_edgV,HO_Fy_edgH,RT_Residual,N_y,N_x,N_m,N_g,&
    !   Tini,comp_unit,nu_g,bcT_left,bcT_right,bcT_top,bcT_bottom,BC_Type,maxit_RTE)
    ! CALL MGQD_INIT(Eg_avg,Eg_edgV,Eg_edgH,Fxg_edgV,Fyg_edgH,Eg_avg_old,Eg_edgV_old,Eg_edgH_old,Fxg_edgV_old,Fyg_edgH_old,&
    !   fg_avg_xx,fg_avg_xy,fg_avg_yy,fg_edgV_xx,fg_edgV_xy,fg_edgH_yy,fg_edgH_xy,fg_avg_xx_old,fg_avg_xy_old,fg_avg_yy_old,&
    !   fg_edgV_xx_old,fg_edgV_xy_old,fg_edgH_yy_old,fg_edgH_xy_old,Cg_L,Cg_B,Cg_R,Cg_T,Eg_in_L,Eg_in_B,Eg_in_R,Eg_in_T,Fg_in_L,&
    !   Fg_in_B,Fg_in_R,Fg_in_T,MGQD_Residual,I_edgV,I_edgH,Omega_x,Omega_y,quad_weight,c,Comp_Unit,N_y,N_x,N_g,MGQD_E_avg,&
    !   MGQD_E_edgV,MGQD_E_edgH,MGQD_Fx_edgV,MGQD_Fy_edgH,G_old,Pold_L,Pold_B,Pold_R,Pold_T,maxit_MLOQD,maxit_RTE,BC_Type)
    ! CALL COLLAPSE_INTENSITIES(Threads,I_avg,I_edgV,I_edgH,Omega_x,Omega_y,quad_weight,Comp_Unit,Hg_avg_xx,Hg_avg_xy,Hg_avg_yy,&
    !   Hg_edgV_xx,Hg_edgV_xy,Hg_edgH_yy,Hg_edgH_xy,Eg_edgV,Eg_edgH,Eg_avg,Fxg_edgV,Fyg_edgH,HO_E_edgV,&
    !   HO_E_edgH,HO_E_avg,HO_Fx_edgV,HO_Fy_edgH)
    ! CALL TEMP_INIT(Temp,RT_Src,MGQD_Src,MGQD_Src_old,KapE,KapB,KapR,KapE_old,KapR_old,Bg,N_y,N_x,N_m,N_g,N_t,Tini,&
    !   Comp_Unit,Nu_g,Temp_Times,HO_E_avg_Times,GREY_E_avg_Times,Temp_Old,MGQD_E_avg_Times)

  !   OPEN(UNIT=103,FILE=res_history_outfile,STATUS='REPLACE',ACTION='WRITE',IOSTAT=err)
  ! !   making sure file exists/opens, if not tells user
  !   IF (err .NE. 0) THEN
  !       WRITE(*,'(3A)') 'The file, ',res_history_outfile,', could not open properly.'
  !       STOP
  !   END IF

! CALL OUTFILE_INIT(outID,N_x_ID,N_y_ID,N_m_ID,N_g_ID,N_t_ID,N_edgV_ID,N_edgH_ID,RT_Its_ID,MGQD_Its_ID,GQD_Its_ID,&
!   Norm_Types_ID,Boundaries_ID,c_ID,h_ID,pi_ID,erg_ID,Comp_Unit_ID,cv_ID,Temp_ID,E_avg_ID,E_edgV_ID,E_edgH_ID,MGQD_E_avg_ID,&
!   MGQD_E_edgV_ID,MGQD_E_edgH_ID,HO_E_avg_ID,HO_E_edgV_ID,HO_E_edgH_ID,Fx_edgV_ID,Fy_edgH_ID,MGQD_Fx_edgV_ID,MGQD_Fy_edgH_ID,&
!   HO_Fx_edgV_ID,HO_Fy_edgH_ID,Eg_avg_ID,Eg_edgV_ID,Eg_edgH_ID,HO_Eg_avg_ID,HO_Eg_edgV_ID,HO_Eg_edgH_ID,Fxg_edgV_ID,Fyg_edgH_ID,&
!   HO_Fxg_edgV_ID,HO_Fyg_edgH_ID,I_avg_ID,I_edgV_ID,I_edgH_ID,RT_Residual_ID,RT_ResNorms_ID,MGQD_Residual_ID,MGQD_ResNorms_ID,&
!   Del_T_ID,Del_E_avg_ID,Del_E_edgV_ID,Del_E_edgH_ID,Del_Fx_edgV_ID,Del_Fy_edgH_ID,Del_T_Norms_ID,Del_E_avg_Norms_ID,&
!   Del_E_edgV_Norms_ID,Del_E_edgH_Norms_ID,Del_Fx_edgV_Norms_ID,Del_Fy_edgH_Norms_ID,RT_ItCount_ID,MGQD_ItCount_ID,GQD_ItCount_ID,&
!   c,h,pi,erg,Comp_Unit,cv,chi,conv_ho,conv_lo,conv_gr,line_src,xlen,ylen,Delx,Dely,tlen,Delt,bcT_left,bcT_bottom,bcT_right,&
!   bcT_top,Tini,N_x,N_y,N_m,N_g,N_t,database_gen,use_grey,maxit_RTE,maxit_MLOQD,maxit_GLOQD,conv_type,threads,BC_type,outfile,&
!   run_type,kapE_dT_flag,quadrature,enrgy_strc)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! CALL TRT_MLQD_ALGORITHM(Omega_x,Omega_y,quad_weight,Delx,Dely,A,Delt,Tlen,c,cV,Kap0,Comp_Unit,Bg,KapE,&
    !   kapB,RT_Src,Temp,Temp_old,I_crn_old,I_avg,I_edgV,I_edgH,I_crn,Ic_edgV,Ic_edgH,Hg_avg_xx,Hg_avg_xy,&
    !   Hg_avg_yy,Hg_edgV_xx,Hg_edgV_xy,Hg_edgH_yy,Hg_edgH_xy,HO_Eg_avg,HO_Eg_edgV,HO_Eg_edgH,HO_Fxg_edgV,HO_Fyg_edgH,&
    !   HO_E_avg,HO_E_edgV,HO_E_edgH,HO_Fx_edgV,HO_Fy_edgH,Temp_Times,HO_E_avg_Times,Nu_g,conv_ho,maxit_RTE,conv_type,&
    !   Start_Time,Threads,BC_Type,bcT_left,bcT_bottom,bcT_right,bcT_top,MGQD_Src,KapR,KapE_old,KapR_old,Theta,Eg_avg,&
    !   Eg_edgV,Eg_edgH,Fxg_edgV,Fyg_edgH,fg_avg_xx,fg_avg_xy,fg_avg_yy,fg_edgV_xx,fg_edgV_xy,fg_edgH_yy,fg_edgH_xy,&
    !   fg_avg_xx_old,fg_avg_xy_old,fg_avg_yy_old,fg_edgV_xx_old,fg_edgV_xy_old,fg_edgH_yy_old,fg_edgH_xy_old,Cg_L,Cg_B,&
    !   Cg_R,Cg_T,Eg_in_L,Eg_in_B,Eg_in_R,Eg_in_T,Fg_in_L,Fg_in_B,Fg_in_R,Fg_in_T,Eg_avg_old,Eg_edgV_old,Eg_edgH_old,&
    !   Fxg_edgV_old,Fyg_edgH_old,MGQD_Src_old,Res_Calc,MGQD_E_avg,MGQD_E_edgV,MGQD_E_edgH,MGQD_Fx_edgV,&
    !   MGQD_Fy_edgH,MGQD_E_avg_Times,Maxit_MLOQD,run_type,Conv_LO,G_old,Pold_L,Pold_B,Pold_R,Pold_T)

    CALL TRT_MLQD_ALGORITHM(Omega_x,Omega_y,quad_weight,Nu_g,Delx,Dely,A,Delt,Theta,Start_Time,c,cV,h,pi,Kap0,&
      erg,Comp_Unit,Conv_ho,Conv_lo,Conv_gr,bcT_left,bcT_bottom,bcT_right,bcT_top,Tini,chi,line_src,database_gen,&
      use_grey,Conv_Type,Maxit_RTE,Threads,BC_Type,Maxit_MLOQD,Maxit_GLOQD,N_x,N_y,N_m,N_g,N_t,Res_Calc,run_type,&
      kapE_dT_flag,outID,N_x_ID,N_y_ID,N_m_ID,N_g_ID,N_t_ID,N_edgV_ID,N_edgH_ID,RT_Its_ID,MGQD_Its_ID,GQD_Its_ID,&
      Norm_Types_ID)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! CLOSE ( 103, STATUS='KEEP')

    ! OPEN(UNIT=10,FILE="output.out",STATUS='REPLACE',ACTION='WRITE',IOSTAT=err)
    ! !   making sure file exists/opens, if not tells user
    ! IF (err .NE. 0) THEN
    !     WRITE(*,'(A)') 'The file, output.out, could not open properly.'
    !     STOP
    ! END IF
    !
    ! CALL OUTFILE_HEADER(10,datevalues,database_gen,database_add,run_type,restart_infile,use_grey,chi,&
    !   conv_ho,conv_lo,conv_gr,comp_unit,line_src,maxit_RTE,maxit_MLOQD,maxit_GLOQD,conv_type,N_m,threads,&
    !   kapE_dT_flag,enrgy_strc,Nwt_upbnd,xlen,ylen,N_x,N_y,tlen,delt,bc_type,bcT_left,bcT_right,bcT_top,&
    !   bcT_bottom,Tini,out_times,out_time_steps,n_out_times,restart_outlen,TEMP_out,GREY_E_out,GREY_F_out,&
    !   GREY_kap_out,GREY_fsmall_out,MG_fsmall_out,res_history_out,outfile,restart_outfile,decomp_outfile,TEMP_outfile,&
    !   GREY_E_outfile,GREY_F_outfile,GREY_kap_outfile,GREY_fsmall_outfile,MG_fsmall_outfile,res_history_outfile,&
    !   HO_E_out,HO_E_outfile,Delx,Dely,N_g,nu_g,N_t,sig_r,ar,Kap0,pi,c,h,cv,Omega_x,Omega_y,quad_weight,quadrature)
    !
    !   CLOSE ( 10, STATUS='KEEP')
    !
    ! CALL SPECIFIC_OUTPUTS(Temp_Times,HO_E_avg_Times,GREY_E_avg_Times,&
    !   Times,10,datevalues,database_gen,database_add,run_type,restart_infile,&
    !   use_grey,chi,conv_ho,conv_lo,conv_gr,comp_unit,line_src,maxit_RTE,maxit_MLOQD,maxit_GLOQD,conv_type,N_m,threads,&
    !   kapE_dT_flag,enrgy_strc,Nwt_upbnd,xlen,ylen,N_x,N_y,tlen,delt,bc_type,bcT_left,bcT_right,bcT_top,bcT_bottom,Tini,&
    !   out_times,out_time_steps,n_out_times,restart_outlen,TEMP_out,GREY_E_out,GREY_F_out,GREY_kap_out,GREY_fsmall_out,&
    !   MG_fsmall_out,res_history_out,outfile,restart_outfile,decomp_outfile,TEMP_outfile,GREY_E_outfile,GREY_F_outfile,&
    !   GREY_kap_outfile,GREY_fsmall_outfile,MG_fsmall_outfile,res_history_outfile,HO_E_out,HO_E_outfile,delx,dely,N_g,&
    !   nu_g,N_t,sig_r,ar,Kap0,pi,c,h,cv,Omega_x,Omega_y,quad_weight,quadrature,MGQD_E_avg_Times)

    CALL SYSTEM_CLOCK(solve_time2,clockrate)
    solve_time=(DBLE(solve_time2)-DBLE(solve_time1))/DBLE(clockrate)

    WRITE(*,*)
    WRITE(*,'(A,ES16.8,A)') "Solve took ",solve_time," seconds."

end
