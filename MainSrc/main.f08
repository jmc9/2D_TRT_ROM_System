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
    REAL*8:: chi, conv_ho, conv_lo, conv_gr1, conv_gr2, line_src, E_Bound_Low, T_Bound_Low
    REAL*8:: xlen, ylen, Tlen, delt, bcT_left, bcT_right, bcT_top, bcT_bottom, Tini
    REAL*8,ALLOCATABLE:: Delx(:), Dely(:), xpts_avg(:), xpts_edgV(:), ypts_avg(:), ypts_edgH(:), tpts(:)
    REAL*8,ALLOCATABLE:: out_times(:), Kappa_Mult(:), dset_times(:), POD_err(:)
    REAL*8,ALLOCATABLE:: nu_g(:)
    INTEGER,ALLOCATABLE:: Mat(:,:)
    INTEGER:: use_grey, Direc_Diff, N_dsets, qdfin_fxy_flag
    INTEGER:: maxit_RTE, maxit_MLOQD, maxit_GLOQD, conv_type, N_m, threads
    INTEGER:: N_x, N_y, N_t, N_g
    INTEGER:: BC_Type(4)
    INTEGER:: out_freq, I_out, HO_Eg_out, HO_Fg_out, HO_E_out, HO_F_out, restart_freq
    INTEGER:: Eg_out, Fg_out, MGQD_E_out, MGQD_F_out, QDfg_out
    INTEGER:: E_out, F_out, D_out
    INTEGER:: old_parms_out, its_out, conv_out, kap_out, Src_out, Init_out
    CHARACTER(100):: run_type, restart_infile, outfile, Test, restart_outfile
    CHARACTER(100):: enrgy_strc, quadrature
    CHARACTER(100):: POD_Type, DMD_Type, Tgen_file, qdf_infile, qdfin_bctype
    CHARACTER(100),ALLOCATABLE:: DMD_dsets(:), POD_dsets(:)
    REAL*8,ALLOCATABLE:: Omega_x(:), Omega_y(:), quad_weight(:)
    REAL*8:: Theta
    LOGICAL:: Res_Calc, kapE_dT_flag
    LOGICAL:: Use_Line_Search, Use_Safety_Search

    !--------------------------------------------------!
    !                                                  !
    !--------------------------------------------------!
    INTEGER:: outID
    INTEGER:: N_x_ID, N_y_ID, N_m_ID, N_g_ID, N_t_ID, N_edgV_ID, N_edgH_ID, N_xc_ID, N_yc_ID, Quads_ID
    INTEGER:: RT_Its_ID, MGQD_Its_ID, GQD_Its_ID, Norm_Types_ID, MGQD_ResTypes_ID, Boundaries_ID, N_dsets_ID
    INTEGER:: c_ID, h_ID, pi_ID, erg_ID, Comp_Unit_ID, cv_ID

    !--------------------------------------------------!
    !                 MISCELLANIOUS                    !
    !--------------------------------------------------!
    INTEGER:: datevalues(8)
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

    CALL INPUT(run_type,restart_infile,use_grey,Test,Mat,Kappa_Mult,chi,conv_ho,conv_lo,conv_gr1,&
      conv_gr2,comp_unit,line_src,E_Bound_Low,T_Bound_Low,Theta,maxit_RTE,maxit_MLOQD,maxit_GLOQD,conv_type,N_m,&
      threads,kapE_dT_flag,enrgy_strc,erg,xlen,ylen,N_x,N_y,tlen,delt,bcT_left,bcT_right,bcT_top,bcT_bottom,&
      Tini,sig_R,ar,pi,c,h,delx,dely,cv,outfile,out_freq,I_out,HO_Eg_out,HO_Fg_out,HO_E_out,HO_F_out,Eg_out,Fg_out,MGQD_E_out,&
      MGQD_F_out,QDfg_out,E_out,F_out,D_out,old_parms_out,its_out,conv_out,kap_out,Src_out,nu_g,N_g,Omega_x,Omega_y,&
      quad_weight,N_t,quadrature,BC_Type,Use_Line_Search,Use_Safety_Search,Res_Calc,POD_err,POD_Type,POD_dsets,&
      Direc_Diff,xpts_avg,xpts_edgV,ypts_avg,ypts_edgH,tpts,N_dsets,DMD_dsets,DMD_Type,restart_outfile,restart_freq,&
      Start_Time,dset_times,Init_out, Tgen_file, qdf_infile, qdfin_fxy_flag, qdfin_bctype)

    IF (run_type .EQ. 'Tgen_qdf') THEN
      CALL Tgen_QDf(Omega_x, Omega_y, quad_weight, Nu_g, Kappa_Mult, c, cV, h, pi, Kap0, erg, Comp_Unit,&
        N_m, Mat, threads, BC_Type, Res_Calc, Tgen_file, outfile, enrgy_strc,&
        I_out, HO_Eg_out, HO_Fg_out, HO_E_out, HO_F_out, QDfg_out, quadrature, run_type)

    ELSE
      CALL OUTFILE_INIT(outID,N_x_ID,N_y_ID,N_m_ID,N_g_ID,N_t_ID,N_edgV_ID,N_edgH_ID,N_xc_ID,N_yc_ID,Quads_ID,RT_Its_ID,&
        MGQD_Its_ID,GQD_Its_ID,Norm_Types_ID,MGQD_ResTypes_ID,Boundaries_ID,c_ID,h_ID,pi_ID,erg_ID,Comp_Unit_ID,cv_ID,&
        c,h,pi,erg,Comp_Unit,cv,chi,conv_ho,conv_lo,conv_gr1,conv_gr2,line_src,xlen,ylen,Delx,Dely,tlen,Delt,bcT_left,&
        bcT_bottom,bcT_right,bcT_top,Tini,E_Bound_Low,T_Bound_Low,N_x,N_y,N_m,N_g,N_t,use_grey,maxit_RTE,&
        maxit_MLOQD,maxit_GLOQD,conv_type,threads,BC_type,outfile,run_type,kapE_dT_flag,quadrature,enrgy_strc,Theta,&
        Use_Line_Search,Use_Safety_Search,I_out,HO_Eg_out,HO_Fg_out,HO_E_out,HO_F_out,Eg_out,Fg_out,MGQD_E_out,MGQD_F_out,&
        QDfg_out,E_out,F_out,D_out,old_parms_out,its_out,conv_out,kap_out,Src_out,POD_Type,Direc_Diff,&
        xpts_avg,xpts_edgV,ypts_avg,ypts_edgH,tpts,N_dsets,DMD_dsets,DMD_Type,dset_times,N_dsets_ID,POD_dsets)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CALL TRT_MLQD_ALGORITHM(Omega_x,Omega_y,quad_weight,Nu_g,Delx,Dely,Delt,tlen,Theta,Start_Time,c,cV,h,pi,Kap0,&
        erg,Comp_Unit,Conv_ho,Conv_lo,Conv_gr1,Conv_gr2,bcT_left,bcT_bottom,bcT_right,bcT_top,Tini,chi,line_src,E_Bound_Low,&
        T_Bound_Low,use_grey,Conv_Type,Maxit_RTE,Threads,BC_Type,Maxit_MLOQD,Maxit_GLOQD,N_x,N_y,N_m,N_g,&
        N_t,Res_Calc,Use_Line_Search,Use_Safety_Search,run_type,kapE_dT_flag,outID,N_x_ID,N_y_ID,N_m_ID,N_g_ID,N_t_ID,&
        N_edgV_ID,N_edgH_ID,N_xc_ID,N_yc_ID,Quads_ID,RT_Its_ID,MGQD_Its_ID,GQD_Its_ID,Norm_Types_ID,MGQD_ResTypes_ID,&
        Boundaries_ID,out_freq,I_out,HO_Eg_out,HO_Fg_out,HO_E_out,HO_F_out,Eg_out,Fg_out,MGQD_E_out,MGQD_F_out,QDfg_out,&
        E_out,F_out,D_out,old_parms_out,its_out,conv_out,kap_out,Src_out,POD_dsets,POD_err,POD_Type,xlen,ylen,&
        Direc_Diff,Mat,Kappa_Mult,restart_outfile,restart_freq,restart_infile,N_dsets,DMD_dsets,DMD_Type,dset_times,&
        N_dsets_ID,Init_out,qdf_infile, qdfin_fxy_flag, qdfin_bctype)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      CALL NF_CLOSE_FILE(outID)

    END IF

    CALL SYSTEM_CLOCK(solve_time2,clockrate)
    solve_time=(DBLE(solve_time2)-DBLE(solve_time1))/DBLE(clockrate)

    WRITE(*,*)
    WRITE(*,'(A,ES16.8,A)') "Solve took ",solve_time," seconds."

end
