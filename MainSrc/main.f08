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
    INTEGER:: outID
    INTEGER:: N_x_ID, N_y_ID, N_m_ID, N_g_ID, N_t_ID, N_edgV_ID, N_edgH_ID, N_xc_ID, N_yc_ID, Quads_ID
    INTEGER:: RT_Its_ID, MGQD_Its_ID, GQD_Its_ID, Norm_Types_ID, MGQD_ResTypes_ID, Boundaries_ID
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

    CALL INPUT(database_gen,database_add,run_type,restart_infile,use_grey,chi,conv_ho,conv_lo,conv_gr,&
      comp_unit,line_src,maxit_RTE,maxit_MLOQD,maxit_GLOQD,conv_type,N_m,threads,kapE_dT_flag,enrgy_strc,&
      Nwt_upbnd,erg,xlen,ylen,N_x,N_y,tlen,delt,bcT_left,bcT_right,bcT_top,bcT_bottom,Tini,sig_R,ar,&
      pi,c,h,delx,dely,cv,out_times,out_time_steps,n_out_times,restart_outlen,TEMP_out,GREY_E_out,GREY_F_out,&
      GREY_kap_out,GREY_fsmall_out,MG_fsmall_out,res_history_out,outfile,restart_outfile,decomp_outfile,TEMP_outfile,&
      GREY_E_outfile,GREY_F_outfile,GREY_kap_outfile,GREY_fsmall_outfile,MG_fsmall_outfile,res_history_outfile,&
      HO_E_out,HO_E_outfile,nu_g,N_g,Omega_x,Omega_y,quad_weight,N_t,quadrature,BC_Type)

    CALL OUTFILE_INIT(outID,N_x_ID,N_y_ID,N_m_ID,N_g_ID,N_t_ID,N_edgV_ID,N_edgH_ID,N_xc_ID,N_yc_ID,Quads_ID,RT_Its_ID,&
      MGQD_Its_ID,GQD_Its_ID,Norm_Types_ID,MGQD_ResTypes_ID,Boundaries_ID,c_ID,h_ID,pi_ID,erg_ID,Comp_Unit_ID,cv_ID,&
      c,h,pi,erg,Comp_Unit,cv,chi,conv_ho,conv_lo,conv_gr,line_src,xlen,ylen,Delx,Dely,tlen,Delt,bcT_left,bcT_bottom,&
      bcT_right,bcT_top,Tini,N_x,N_y,N_m,N_g,N_t,database_gen,use_grey,maxit_RTE,maxit_MLOQD,maxit_GLOQD,conv_type,&
      threads,BC_type,outfile,run_type,kapE_dT_flag,quadrature,enrgy_strc)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Start_Time=0d0
    CALL TRT_MLQD_ALGORITHM(Omega_x,Omega_y,quad_weight,Nu_g,Delx,Dely,Delt,tlen,Theta,Start_Time,c,cV,h,pi,Kap0,&
      erg,Comp_Unit,Conv_ho,Conv_lo,Conv_gr,bcT_left,bcT_bottom,bcT_right,bcT_top,Tini,chi,line_src,database_gen,&
      use_grey,Conv_Type,Maxit_RTE,Threads,BC_Type,Maxit_MLOQD,Maxit_GLOQD,N_x,N_y,N_m,N_g,N_t,Res_Calc,run_type,&
      kapE_dT_flag,outID,N_x_ID,N_y_ID,N_m_ID,N_g_ID,N_t_ID,N_edgV_ID,N_edgH_ID,N_xc_ID,N_yc_ID,Quads_ID,RT_Its_ID,&
      MGQD_Its_ID,GQD_Its_ID,Norm_Types_ID,MGQD_ResTypes_ID,Boundaries_ID)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    CALL SYSTEM_CLOCK(solve_time2,clockrate)
    solve_time=(DBLE(solve_time2)-DBLE(solve_time1))/DBLE(clockrate)

    WRITE(*,*)
    WRITE(*,'(A,ES16.8,A)') "Solve took ",solve_time," seconds."

end