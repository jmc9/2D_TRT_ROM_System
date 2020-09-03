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
  erg,Comp_Unit,Conv_ho,Conv_lo,Conv_gr,bcT_left,bcT_bottom,bcT_right,bcT_top,Tini,chi,line_src,database_gen,&
  use_grey,Conv_Type,Maxit_RTE,Threads,BC_Type,Maxit_MLOQD,Maxit_GLOQD,N_x,N_y,N_m,N_g,N_t,Res_Calc,run_type,&
  kapE_dT_flag,outID,N_x_ID,N_y_ID,N_m_ID,N_g_ID,N_t_ID,N_edgV_ID,N_edgH_ID,N_xc_ID,N_yc_ID,Quads_ID,RT_Its_ID,&
  MGQD_Its_ID,GQD_Its_ID,Norm_Types_ID,MGQD_ResTypes_ID,Boundaries_ID)

  !---------------Solution Parameters----------------!
  REAL*8,INTENT(IN):: Omega_x(:), Omega_y(:), quad_weight(:), Nu_g(:)
  REAL*8,INTENT(IN):: Delx(:), Dely(:), Delt, Theta, tlen
  REAL*8,INTENT(IN):: Start_Time
  REAL*8,INTENT(IN):: c, cV, h, pi, Kap0, erg
  REAL*8,INTENT(IN):: Comp_Unit, Conv_ho, Conv_lo, Conv_gr
  REAL*8,INTENT(IN):: bcT_left, bcT_bottom, bcT_right, bcT_top, Tini
  REAL*8,INTENT(IN):: chi, line_src
  INTEGER,INTENT(IN):: N_x, N_y, N_m, N_g, N_t
  INTEGER,INTENT(IN):: database_gen, use_grey, Conv_Type, Threads, BC_Type(:), Maxit_RTE, Maxit_MLOQD, Maxit_GLOQD
  LOGICAL,INTENT(IN):: Res_Calc
  CHARACTER(*),INTENT(IN):: run_type, kapE_dT_flag

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
  ! REAL*8,ALLOCATABLE:: E_avg(:,:), E_edgV(:,:), E_edgH(:,:)
  ! REAL*8,ALLOCATABLE:: Fx_edgV(:,:), Fy_edgH(:,:)

  !--------------------------------------------------!
  REAL*8,ALLOCATABLE:: RT_Residual(:,:,:), MGQD_Residual(:,:,:), MGQD_BC_Residual(:,:)
  INTEGER,ALLOCATABLE:: RT_ResLoc_x(:,:,:), RT_ResLoc_y(:,:,:), MGQD_ResLoc_x(:,:,:), MGQD_ResLoc_y(:,:,:)
  REAL*8,ALLOCATABLE:: Temp_RTold(:,:), Temp_RTold2(:,:)
  REAL*8,ALLOCATABLE:: Temp_MGQDold(:,:), Temp_MGQDold2(:,:)
  REAL*8,ALLOCATABLE:: HO_E_avg_RTold(:,:), HO_E_avg_RTold2(:,:)
  REAL*8,ALLOCATABLE:: MGQD_E_avg_MGQDold(:,:), MGQD_E_avg_MGQDold2(:,:)
  REAL*8:: TR_Tnorm, TR_Enorm, TR_Trho, TR_Erho
  REAL*8:: MGQD_Tnorm, MGQD_Enorm, MGQD_Trho, MGQD_Erho
  REAL*8:: Time
  INTEGER:: MGQD_Its, Status
  INTEGER:: RT_Its, RT_start_Its, t
  LOGICAL:: RT_Conv, MGQD_conv, Tconv, Econv

  !-------------------Variable ID's------------------!
  INTEGER:: Temp_ID, E_avg_ID, E_edgV_ID, E_edgH_ID, MGQD_E_avg_ID, MGQD_E_edgV_ID, MGQD_E_edgH_ID, HO_E_avg_ID
  INTEGER:: HO_E_edgV_ID, HO_E_edgH_ID, Fx_edgV_ID, Fy_edgH_ID, MGQD_Fx_edgV_ID, MGQD_Fy_edgH_ID, HO_Fx_edgV_ID
  INTEGER:: HO_Fy_edgH_ID, Eg_avg_ID, Eg_edgV_ID, Eg_edgH_ID, HO_Eg_avg_ID, HO_Eg_edgV_ID, HO_Eg_edgH_ID
  INTEGER:: Fxg_edgV_ID, Fyg_edgH_ID, HO_Fxg_edgV_ID, HO_Fyg_edgH_ID, I_avg_ID, I_edgV_ID, I_edgH_ID
  INTEGER:: RT_Residual_ID, MGQD_Residual_ID, MGQD_BC_Residual_ID
  INTEGER:: Del_T_ID, Del_E_avg_ID, Del_E_edgV_ID, Del_E_edgH_ID, Del_Fx_edgV_ID, Del_Fy_edgH_ID, Del_T_Norms_ID
  INTEGER:: Del_E_avg_Norms_ID, Del_E_edgV_Norms_ID, Del_E_edgH_Norms_ID, Del_Fx_edgV_Norms_ID, Del_Fy_edgH_Norms_ID
  INTEGER:: RT_ItCount_ID, MGQD_ItCount_ID, GQD_ItCount_ID

  !===========================================================================!
  !                                                                           !
  !     INITIALIZING ARRAYS                                                   !
  !                                                                           !
  !===========================================================================!
  CALL MISC_INIT(Delx,Dely,A)
  CALL RT_INIT(I_avg,I_edgV,I_edgH,I_crn,I_crn_old,Ic_edgV,Ic_edgH,Hg_avg_xx,Hg_avg_yy,Hg_edgV_xx,&
    Hg_edgV_xy,Hg_edgH_yy,Hg_edgH_xy,HO_Eg_avg,HO_Eg_edgV,HO_Eg_edgH,HO_Fxg_edgV,HO_Fyg_edgH,HO_E_avg,&
    HO_E_edgV,HO_E_edgH,HO_Fx_edgV,HO_Fy_edgH,N_y,N_x,N_m,N_g,Tini,comp_unit,nu_g,bcT_left,bcT_right,&
    bcT_top,bcT_bottom,BC_Type)
  CALL MGQD_INIT(Eg_avg,Eg_edgV,Eg_edgH,Fxg_edgV,Fyg_edgH,Eg_avg_old,Eg_edgV_old,Eg_edgH_old,Fxg_edgV_old,Fyg_edgH_old,&
    fg_avg_xx,fg_avg_yy,fg_edgV_xx,fg_edgV_xy,fg_edgH_yy,fg_edgH_xy,fg_avg_xx_old,fg_avg_yy_old,fg_edgV_xx_old,fg_edgV_xy_old,&
    fg_edgH_yy_old,fg_edgH_xy_old,Cg_L,Cg_B,Cg_R,Cg_T,Eg_in_L,Eg_in_B,Eg_in_R,Eg_in_T,Fg_in_L,Fg_in_B,Fg_in_R,Fg_in_T,I_edgV,&
    I_edgH,Omega_x,Omega_y,quad_weight,c,Comp_Unit,N_y,N_x,N_g,MGQD_E_avg,MGQD_E_edgV,MGQD_E_edgH,MGQD_Fx_edgV,MGQD_Fy_edgH,&
    G_old,Pold_L,Pold_B,Pold_R,Pold_T,BC_Type)
  CALL COLLAPSE_INTENSITIES(Threads,I_avg,I_edgV,I_edgH,Omega_x,Omega_y,quad_weight,Comp_Unit,Hg_avg_xx,Hg_avg_yy,&
    Hg_edgV_xx,Hg_edgV_xy,Hg_edgH_yy,Hg_edgH_xy,Eg_edgV,Eg_edgH,Eg_avg,Fxg_edgV,Fyg_edgH,HO_E_edgV,&
    HO_E_edgH,HO_E_avg,HO_Fx_edgV,HO_Fy_edgH)
  CALL TEMP_INIT(Temp,RT_Src,MGQD_Src,MGQD_Src_old,KapE,KapB,KapR,KapE_old,KapR_old,Bg,N_y,N_x,N_m,N_g,Tini,&
    Comp_Unit,Nu_g,Temp_Old)

  !===========================================================================!
  !                                                                           !
  !     INITIALIZING OUTPUT FILE                                              !
  !                                                                           !
  !===========================================================================!
  CALL OUTFILE_VARDEFS(outID,N_x_ID,N_y_ID,N_m_ID,N_g_ID,N_t_ID,N_edgV_ID,N_edgH_ID,N_xc_ID,N_yc_ID,Quads_ID,&
    RT_Its_ID,MGQD_Its_ID,GQD_Its_ID,Norm_Types_ID,MGQD_ResTypes_ID,Boundaries_ID,Temp_ID,E_avg_ID,E_edgV_ID,E_edgH_ID,&
    MGQD_E_avg_ID,MGQD_E_edgV_ID,MGQD_E_edgH_ID,HO_E_avg_ID,HO_E_edgV_ID,HO_E_edgH_ID,Fx_edgV_ID,Fy_edgH_ID,&
    MGQD_Fx_edgV_ID,MGQD_Fy_edgH_ID,HO_Fx_edgV_ID,HO_Fy_edgH_ID,Eg_avg_ID,Eg_edgV_ID,Eg_edgH_ID,HO_Eg_avg_ID,&
    HO_Eg_edgV_ID,HO_Eg_edgH_ID,Fxg_edgV_ID,Fyg_edgH_ID,HO_Fxg_edgV_ID,HO_Fyg_edgH_ID,I_avg_ID,I_edgV_ID,I_edgH_ID,&
    RT_Residual_ID,MGQD_Residual_ID,MGQD_BC_Residual_ID,Del_T_ID,Del_E_avg_ID,Del_E_edgV_ID,Del_E_edgH_ID,&
    Del_Fx_edgV_ID,Del_Fy_edgH_ID,Del_T_Norms_ID,Del_E_avg_Norms_ID,Del_E_edgV_Norms_ID,Del_E_edgH_Norms_ID,&
    Del_Fx_edgV_Norms_ID,Del_Fy_edgH_Norms_ID,RT_ItCount_ID,MGQD_ItCount_ID,GQD_ItCount_ID)

  ALLOCATE(Temp_RTold(SIZE(TEMP,1),SIZE(TEMP,2)))
  ALLOCATE(Temp_RTold2(SIZE(TEMP,1),SIZE(TEMP,2)))
  ALLOCATE(Temp_MGQDold(SIZE(TEMP,1),SIZE(TEMP,2)))
  ALLOCATE(Temp_MGQDold2(SIZE(TEMP,1),SIZE(TEMP,2)))
  ALLOCATE(HO_E_avg_RTold(SIZE(HO_E_avg,1),SIZE(HO_E_avg,2)))
  ALLOCATE(HO_E_avg_RTold2(SIZE(HO_E_avg,1),SIZE(HO_E_avg,2)))
  ALLOCATE(MGQD_E_avg_MGQDold(SIZE(MGQD_E_avg,1),SIZE(MGQD_E_avg,2)))
  ALLOCATE(MGQD_E_avg_MGQDold2(SIZE(MGQD_E_avg,1),SIZE(MGQD_E_avg,2)))

  ALLOCATE(RT_Residual(4,N_g,5),RT_ResLoc_x(4,N_g,5),RT_ResLoc_y(4,N_g,5))
  ALLOCATE(MGQD_Residual(N_g,5,5),MGQD_ResLoc_x(N_g,5,5),MGQD_ResLoc_y(N_g,5,5),MGQD_BC_Residual(N_g,4))

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

    t = t + 1
    Time = Time + Delt

    Temp_Old = Temp
    I_crn_old = I_crn

    fg_avg_xx_old = fg_avg_xx
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
    ! CALL OLD_GREY_COEFS(c,Delt,Theta,cv,A,Gold,Temp_old,KapE_bar_old,E_avg_old,GQD_Src_old,Gold_Hat,Rhat_old)

    HO_E_avg_RTold2 = 0d0
    HO_E_avg_RTold = 0d0
    Temp_RTold2 = 0d0
    Temp_RTold = 0d0
    RT_Its = 0
    MGQD_Its = 0
    RT_Conv = .FALSE.
    DO WHILE ((.NOT. RT_Conv).AND.(RT_Its .LT. Maxit_RTE))
      RT_Its = RT_Its + 1

      IF (RT_Its .GE. RT_start_Its) THEN

        CALL TRANSPORT_SCB(I_avg,I_edgV,I_edgH,I_crn,Ic_edgV,Ic_edgH,Omega_x,Omega_y,Delx,Dely,A,&
          KapE,RT_Src,I_crn_old,c,Delt,Threads,RT_Residual,RT_ResLoc_x,RT_ResLoc_y,Res_Calc)

        Status = nf90_put_var(outID,RT_Residual_ID,RT_Residual,(/1,1,1,RT_Its,t/),(/4,N_g,5,1,1/))
        CALL HANDLE_ERR(Status)

        CALL COLLAPSE_INTENSITIES(Threads,I_avg,I_edgV,I_edgH,Omega_x,Omega_y,quad_weight,Comp_Unit,Hg_avg_xx,Hg_avg_yy,&
          Hg_edgV_xx,Hg_edgV_xy,Hg_edgH_yy,Hg_edgH_xy,HO_Eg_edgV,HO_Eg_edgH,HO_Eg_avg,HO_Fxg_edgV,HO_Fyg_edgH,HO_E_edgV,&
          HO_E_edgH,HO_E_avg,HO_Fx_edgV,HO_Fy_edgH)

        IF (MAXVAL(BC_Type) .GT. 0) THEN
          CALL RT_BC_UPDATE(BC_Type,bcT_left,bcT_bottom,bcT_right,bcT_top,Comp_Unit,Nu_g,I_edgV,I_edgH,Ic_edgV,Ic_edgH)
          CALL MGQD_In_Calc(Eg_in_L,Eg_in_B,Eg_in_R,Eg_in_T,Fg_in_L,Fg_in_B,Fg_in_R,Fg_in_T,I_edgV,I_edgH,Omega_x,Omega_y,&
            quad_weight,c,Comp_Unit)
        END IF

        CALL fg_Calc(fg_avg_xx,fg_avg_yy,fg_edgV_xx,fg_edgV_xy,fg_edgH_yy,fg_edgH_xy,Hg_avg_xx,Hg_avg_yy,&
          Hg_edgV_xx,Hg_edgV_xy,Hg_edgH_yy,Hg_edgH_xy,HO_Eg_avg,HO_Eg_edgV,HO_Eg_edgH,c,Comp_Unit)

        CALL Cg_Calc(Cg_L,Cg_B,Cg_R,Cg_T,I_edgV,I_edgH,Omega_x,Omega_y,quad_weight,Comp_Unit,BC_Type)

      END IF

      IF ( run_type .EQ. 'tr_no_qd' ) THEN
        CALL MEB_SOLVE(Temp,HO_Eg_avg,Bg,KapE,Temp_old,Delt,cV,Kap0,Comp_Unit)
        CALL MATERIAL_UPDATE(RT_Src,MGQD_Src,KapE,KapB,KapR,Bg,Temp,Comp_Unit,Nu_g)

      ELSE
        MGQD_E_avg_MGQDold2 = 0d0
        MGQD_E_avg_MGQDold = HO_E_avg
        Temp_MGQDold2 = 0d0
        Temp_MGQDold = Temp
        MGQD_Its = 0
        MGQD_conv = .FALSE.
        DO WHILE ((.NOT. MGQD_conv).AND.(MGQD_Its .LT. Maxit_MLOQD))
          MGQD_Its = MGQD_Its + 1

          CALL MLOQD_FV(Eg_avg,Eg_edgV,Eg_edgH,Fxg_edgV,Fyg_edgH,fg_avg_xx,fg_avg_yy,fg_edgV_xx,fg_edgV_xy,fg_edgH_yy,&
            fg_edgH_xy,Cg_L,Cg_B,Cg_R,Cg_T,Eg_in_L,Eg_in_B,Eg_in_R,Eg_in_T,Fg_in_L,Fg_in_B,Fg_in_R,Fg_in_T,MGQD_Src,KapE,KapR,&
            Delx,Dely,A,c,Delt,Theta,Threads,Res_Calc,MGQD_Residual,MGQD_ResLoc_x,MGQD_ResLoc_y,MGQD_BC_Residual,G_old,Pold_L,&
            Pold_B,Pold_R,Pold_T,Eg_avg_old,Fxg_edgV_old,Fyg_edgH_old)
          CALL COLLAPSE_MG_EF(Eg_avg,Eg_edgV,Eg_edgH,Fxg_edgV,Fyg_edgH,Threads,MGQD_E_avg,MGQD_E_edgV,MGQD_E_edgH,MGQD_Fx_edgV,&
            MGQD_Fy_edgH)

          !Causes errors sometimes (bug somewhere)
          Status = nf90_put_var(outID,MGQD_BC_Residual_ID,MGQD_BC_Residual,(/1,1,MGQD_Its,RT_Its,t/),(/N_g,4,1,1,1/))
          CALL HANDLE_ERR(Status)
          Status = nf90_put_var(outID,MGQD_Residual_ID,MGQD_Residual,(/1,1,1,MGQD_Its,RT_Its,t/),(/N_g,5,5,1,1,1/))
          CALL HANDLE_ERR(Status)

          CALL MEB_SOLVE(Temp,Eg_avg,Bg,KapE,Temp_old,Delt,cV,Kap0,Comp_Unit)
          CALL MATERIAL_UPDATE(RT_Src,MGQD_Src,KapE,KapB,KapR,Bg,Temp,Comp_Unit,Nu_g)

          Tconv = CONVERGENCE(Conv_Type,conv_lo,Temp,Temp_MGQDold,Temp_MGQDold2,MGQD_Its,&
            MGQD_Tnorm,MGQD_Trho)
          Econv = CONVERGENCE(Conv_Type,conv_lo,MGQD_E_avg,MGQD_E_avg_MGQDold,MGQD_E_avg_MGQDold2,MGQD_Its,&
            MGQD_Enorm,MGQD_Erho)
          MGQD_conv = Tconv.AND.Econv

          MGQD_E_avg_MGQDold2 = MGQD_E_avg_MGQDold
          MGQD_E_avg_MGQDold = MGQD_E_avg
          Temp_MGQDold2 = Temp_MGQDold
          Temp_MGQDold = Temp

          write(*,*) 'MGQD:     ',MGQD_Its, MGQD_Tnorm, MGQD_Enorm

        END DO

      END IF

      Tconv = CONVERGENCE(Conv_Type,conv_ho,Temp,Temp_RTold,Temp_RTold2,RT_Its,TR_Tnorm,TR_Trho)
      Econv = CONVERGENCE(Conv_Type,conv_ho,HO_E_avg,HO_E_avg_RTold,HO_E_avg_RTold2,RT_Its,TR_Enorm,TR_Erho)
      RT_conv = Tconv.AND.Econv

      HO_E_avg_RTold2 = HO_E_avg_RTold
      HO_E_avg_RTold = HO_E_avg
      Temp_RTold2 = Temp_RTold
      Temp_RTold = Temp

      write(*,*) RT_Its, MGQD_Its, TR_Tnorm, TR_Enorm
      IF (RT_Its .GE. Maxit_RTE) EXIT

    END DO
    write(*,*) Time, RT_Its, TR_Tnorm, TR_Enorm

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

    IF (Time .GE. tlen) EXIT
  END DO

END SUBROUTINE TRT_MLQD_ALGORITHM

!==================================================================================================================================!
!
!==================================================================================================================================!

END MODULE ALGORITHMS
