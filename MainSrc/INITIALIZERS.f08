MODULE INITIALIZERS

  USE UPDATES
  USE LA_TOOLS
  USE MLOQD_SOLVES

  IMPLICIT NONE

CONTAINS

!==================================================================================================================================!
!
!==================================================================================================================================!
SUBROUTINE RT_INIT(I_avg,I_edgV,I_edgH,I_crn,I_crn_old,Ic_edgV,Ic_edgH,Hg_avg_xx,Hg_avg_yy,Hg_edgV_xx,&
  Hg_edgV_xy,Hg_edgH_yy,Hg_edgH_xy,HO_Eg_avg,HO_Eg_edgV,HO_Eg_edgH,HO_Fxg_edgV,HO_Fyg_edgH,HO_E_avg,&
  HO_E_edgV,HO_E_edgH,HO_Fx_edgV,HO_Fy_edgH,N_y,N_x,N_m,N_g,Tini,comp_unit,nu_g,bcT_left,bcT_right,&
  bcT_upper,bcT_lower,BC_Type,pi,c)
  REAL*8,ALLOCATABLE,INTENT(OUT):: I_avg(:,:,:,:), I_edgV(:,:,:,:), I_edgH(:,:,:,:)
  REAL*8,ALLOCATABLE,INTENT(OUT):: I_crn(:,:,:,:), I_crn_old(:,:,:,:), Ic_edgV(:,:,:,:), Ic_edgH(:,:,:,:)
  REAL*8,ALLOCATABLE,INTENT(OUT):: Hg_avg_xx(:,:,:), Hg_avg_yy(:,:,:)
  REAL*8,ALLOCATABLE,INTENT(OUT):: Hg_edgV_xx(:,:,:), Hg_edgV_xy(:,:,:)
  REAL*8,ALLOCATABLE,INTENT(OUT):: Hg_edgH_yy(:,:,:), Hg_edgH_xy(:,:,:)
  REAL*8,ALLOCATABLE,INTENT(OUT):: HO_Eg_avg(:,:,:), HO_Eg_edgV(:,:,:), HO_Eg_edgH(:,:,:)
  REAL*8,ALLOCATABLE,INTENT(OUT):: HO_Fxg_edgV(:,:,:), HO_Fyg_edgH(:,:,:)
  REAL*8,ALLOCATABLE,INTENT(OUT):: HO_E_avg(:,:), HO_E_edgV(:,:), HO_E_edgH(:,:)
  REAL*8,ALLOCATABLE,INTENT(OUT):: HO_Fx_edgV(:,:), HO_Fy_edgH(:,:)

  REAL*8,INTENT(IN):: Tini, bcT_left, bcT_right, bcT_upper, bcT_lower
  REAL*8,INTENT(IN):: comp_unit, nu_g(:), pi, c
  INTEGER,INTENT(IN):: N_y, N_x, N_m, N_g, BC_Type(:)

  REAL*8:: bg
  INTEGER:: g

  !--------------------------------------------------!
  !             Allocating array sizes               !
  !--------------------------------------------------!

  !Intensities along the spatial grid
  ALLOCATE(I_avg(N_x,N_y,N_m,N_g))
  ALLOCATE(I_edgV(N_x+1,N_y,N_m,N_g))
  ALLOCATE(I_edgH(N_x,N_y+1,N_m,N_g))

  !Intensities along the SCB subcell grid
  ALLOCATE(I_crn(N_x*2,N_y*2,N_m,N_g))
  ALLOCATE(I_crn_old(N_x*2,N_y*2,N_m,N_g))
  ALLOCATE(Ic_edgV(N_x+1,N_y*2,N_m,N_g))
  ALLOCATE(Ic_edgH(N_x*2,N_y+1,N_m,N_g))

  !Multigroup second angular moment tensors
  ALLOCATE(Hg_avg_xx(N_x,N_y,N_g))
  ALLOCATE(Hg_avg_yy(N_x,N_y,N_g))
  ALLOCATE(Hg_edgV_xx(N_x+1,N_y,N_g))
  ALLOCATE(Hg_edgV_xy(N_x+1,N_y,N_g))
  ALLOCATE(Hg_edgH_yy(N_x,N_y+1,N_g))
  ALLOCATE(Hg_edgH_xy(N_x,N_y+1,N_g))

  !Multigroup rad fluxes
  ALLOCATE(HO_Fxg_edgV(N_x+1,N_y,N_g))
  ALLOCATE(HO_Fyg_edgH(N_x,N_y+1,N_g))

  !Multigroup rad energy densities
  ALLOCATE(HO_Eg_avg(N_x,N_y,N_g))
  ALLOCATE(HO_Eg_edgV(N_x+1,N_y,N_g))
  ALLOCATE(HO_Eg_edgH(N_x,N_y+1,N_g))

  !Total rad fluxes
  ALLOCATE(HO_Fx_edgV(N_x+1,N_y))
  ALLOCATE(HO_Fy_edgH(N_x,N_y+1))

  !Total rad energy densities
  ALLOCATE(HO_E_avg(N_x,N_y))
  ALLOCATE(HO_E_edgV(N_x+1,N_y))
  ALLOCATE(HO_E_edgH(N_x,N_y+1))

  !Residuals for the transport solve
  ! ALLOCATE(RT_Residual(N_x*2,N_y*2,N_m,N_g,2,maxit_RTE))

  HO_E_avg = 0d0
  HO_E_edgV = 0d0
  HO_E_edgH = 0d0
  DO g=1,N_g
    !black body radiation distribution at initial Temp
    bg = Bg_planck_calc(Tini,nu_g(g),nu_g(g+1),comp_unit)
    I_avg(:,:,:,g) = bg
    I_edgV(:,:,:,g) = bg
    I_edgH(:,:,:,g) = bg
    I_crn(:,:,:,g) = bg
    Ic_edgV(:,:,:,g) = bg
    Ic_edgH(:,:,:,g) = bg
    HO_Eg_avg(:,:,g) = 4d0*pi*bg/c
    HO_Eg_edgV(:,:,g) = 4d0*pi*bg/c
    HO_Eg_edgH(:,:,g) = 4d0*pi*bg/c
    HO_E_avg = HO_E_avg + 4d0*pi*bg/c
    HO_E_edgV = HO_E_edgV + 4d0*pi*bg/c
    HO_E_edgH = HO_E_edgH + 4d0*pi*bg/c
  END DO
  HO_Fxg_edgV = 0d0
  HO_Fyg_edgH = 0d0
  HO_Fx_edgV = 0d0
  HO_Fy_edgH = 0d0

  ! CALL RT_BC_UPDATE(BC_Type,bcT_left,bcT_lower,bcT_right,bcT_upper,Comp_Unit,Nu_g,I_edgV,I_edgH,Ic_edgV,Ic_edgH)
  CALL RT_BC_UPDATE((/0,0,0,0/),bcT_left,bcT_lower,bcT_right,bcT_upper,Comp_Unit,Nu_g,I_edgV,I_edgH,Ic_edgV,Ic_edgH)
  I_crn_old = I_crn

END SUBROUTINE RT_INIT

!==================================================================================================================================!
!
!==================================================================================================================================!
SUBROUTINE TEMP_INIT(Temp,RT_Src,MGQD_Src,MGQD_Src_old,KapE,KapB,KapR,KapE_old,KapR_old,Bg,N_y,N_x,N_m,N_g,Tini,&
  Comp_Unit,Nu_g,Temp_Old,Threads)
  REAL*8,ALLOCATABLE,INTENT(INOUT):: Temp(:,:), Temp_Old(:,:)
  REAL*8,ALLOCATABLE,INTENT(OUT):: RT_Src(:,:,:,:), MGQD_Src(:,:,:), MGQD_Src_old(:,:,:)
  REAL*8,ALLOCATABLE,INTENT(OUT):: KapE(:,:,:), KapB(:,:,:), KapR(:,:,:), Bg(:,:,:)
  REAL*8,ALLOCATABLE,INTENT(OUT):: KapE_old(:,:,:), KapR_old(:,:,:)

  REAL*8,INTENT(IN):: Tini, comp_unit, nu_g(:)
  INTEGER,INTENT(IN):: N_y, N_x, N_m, N_g, Threads

  ALLOCATE(RT_Src(N_x*2,N_y*2,N_m,N_g))
  ALLOCATE(MGQD_Src(N_x,N_y,N_g))
  ALLOCATE(MGQD_Src_old(N_x,N_y,N_g))
  ALLOCATE(KapE(N_x,N_y,N_g))
  ALLOCATE(KapB(N_x,N_y,N_g))
  ALLOCATE(KapR(N_x,N_y,N_g))
  ALLOCATE(KapE_old(N_x,N_y,N_g))
  ALLOCATE(KapR_old(N_x,N_y,N_g))
  ALLOCATE(Bg(N_x,N_y,N_g))
  ALLOCATE(Temp(N_x,N_y))
  ALLOCATE(Temp_Old(N_x,N_y))

  Temp = Tini
  ! Temp_Old = Temp
  CALL MATERIAL_UPDATE(RT_Src,MGQD_Src,KapE,KapB,KapR,Bg,Temp,Comp_Unit,Nu_g,Threads)

  ! KapE_old = KapE
  ! KapR_old = KapR
  ! MGQD_Src_old = MGQD_Src

END SUBROUTINE TEMP_INIT

!==================================================================================================================================!
!
!==================================================================================================================================!
SUBROUTINE MISC_INIT(Delx,Dely,A)
  REAL*8,INTENT(IN):: Delx(:), Dely(:)
  REAL*8,ALLOCATABLE,INTENT(OUT):: A(:,:)
  INTEGER:: N_x, N_y

  N_x = SIZE(Delx,1)
  N_y = SIZE(Dely,1)
  ALLOCATE(A(N_x,N_y))
  CALL OUTER_PROD(A,Delx,Dely)

END SUBROUTINE MISC_INIT

!==================================================================================================================================!
!
!==================================================================================================================================!
SUBROUTINE MGQD_INIT(Eg_avg,Eg_edgV,Eg_edgH,Fxg_edgV,Fyg_edgH,Eg_avg_old,Eg_edgV_old,Eg_edgH_old,Fxg_edgV_old,Fyg_edgH_old,&
  fg_avg_xx,fg_avg_yy,fg_edgV_xx,fg_edgV_xy,fg_edgH_yy,fg_edgH_xy,fg_avg_xx_old,fg_avg_yy_old,fg_edgV_xx_old,fg_edgV_xy_old,&
  fg_edgH_yy_old,fg_edgH_xy_old,Cg_L,Cg_B,Cg_R,Cg_T,Eg_in_L,Eg_in_B,Eg_in_R,Eg_in_T,Fg_in_L,Fg_in_B,Fg_in_R,Fg_in_T,I_edgV,&
  I_edgH,Omega_x,Omega_y,quad_weight,c,Comp_Unit,N_y,N_x,N_g,MGQD_E_avg,MGQD_E_edgV,MGQD_E_edgH,MGQD_Fx_edgV,MGQD_Fy_edgH,&
  G_old,Pold_L,Pold_B,Pold_R,Pold_T,BC_Type,Tini,nu_g,pi,Open_Threads)

  REAL*8,ALLOCATABLE,INTENT(OUT):: Eg_avg(:,:,:), Eg_edgV(:,:,:), Eg_edgH(:,:,:)
  REAL*8,ALLOCATABLE,INTENT(OUT):: Fxg_edgV(:,:,:), Fyg_edgH(:,:,:)
  REAL*8,ALLOCATABLE,INTENT(OUT):: Eg_avg_old(:,:,:), Eg_edgV_old(:,:,:), Eg_edgH_old(:,:,:)
  REAL*8,ALLOCATABLE,INTENT(OUT):: Fxg_edgV_old(:,:,:), Fyg_edgH_old(:,:,:)
  REAL*8,ALLOCATABLE,INTENT(OUT):: fg_avg_xx(:,:,:), fg_avg_yy(:,:,:)
  REAL*8,ALLOCATABLE,INTENT(OUT):: fg_edgV_xx(:,:,:), fg_edgV_xy(:,:,:)
  REAL*8,ALLOCATABLE,INTENT(OUT):: fg_edgH_yy(:,:,:), fg_edgH_xy(:,:,:)
  REAL*8,ALLOCATABLE,INTENT(OUT):: fg_avg_xx_old(:,:,:), fg_avg_yy_old(:,:,:)
  REAL*8,ALLOCATABLE,INTENT(OUT):: fg_edgV_xx_old(:,:,:), fg_edgV_xy_old(:,:,:)
  REAL*8,ALLOCATABLE,INTENT(OUT):: fg_edgH_yy_old(:,:,:), fg_edgH_xy_old(:,:,:)
  REAL*8,ALLOCATABLE,INTENT(OUT):: Cg_L(:,:), Cg_B(:,:), Cg_R(:,:), Cg_T(:,:)
  REAL*8,ALLOCATABLE,INTENT(OUT):: Eg_in_L(:,:), Eg_in_B(:,:), Eg_in_R(:,:), Eg_in_T(:,:)
  REAL*8,ALLOCATABLE,INTENT(OUT):: Fg_in_L(:,:), Fg_in_B(:,:), Fg_in_R(:,:), Fg_in_T(:,:)
  REAL*8,ALLOCATABLE,INTENT(OUT):: MGQD_E_avg(:,:), MGQD_E_edgV(:,:), MGQD_E_edgH(:,:)
  REAL*8,ALLOCATABLE,INTENT(OUT):: MGQD_Fx_edgV(:,:), MGQD_Fy_edgH(:,:)
  REAL*8,ALLOCATABLE,INTENT(OUT):: G_old(:,:,:), Pold_L(:,:,:), Pold_B(:,:,:), Pold_R(:,:,:), Pold_T(:,:,:)

  REAL*8,INTENT(IN):: I_edgV(:,:,:,:), I_edgH(:,:,:,:)
  REAL*8,INTENT(IN):: Omega_x(:), Omega_y(:), quad_weight(:)
  REAL*8,INTENT(IN):: Tini, nu_g(:)
  REAL*8,INTENT(IN):: c, Comp_Unit, pi
  INTEGER,INTENT(IN):: N_y, N_x, N_g
  INTEGER,INTENT(IN):: BC_Type(:), Open_Threads

  REAL*8:: bg
  INTEGER:: g

  !Multigroup quasidiffusion tensors
  ALLOCATE(fg_avg_xx(N_x,N_y,N_g))
  ALLOCATE(fg_avg_yy(N_x,N_y,N_g))
  ALLOCATE(fg_edgV_xx(N_x+1,N_y,N_g))
  ALLOCATE(fg_edgV_xy(N_x+1,N_y,N_g))
  ALLOCATE(fg_edgH_yy(N_x,N_y+1,N_g))
  ALLOCATE(fg_edgH_xy(N_x,N_y+1,N_g))

  !Multigroup quasidiffusion tensors (last time step)
  ALLOCATE(fg_avg_xx_old(N_x,N_y,N_g))
  ALLOCATE(fg_avg_yy_old(N_x,N_y,N_g))
  ALLOCATE(fg_edgV_xx_old(N_x+1,N_y,N_g))
  ALLOCATE(fg_edgV_xy_old(N_x+1,N_y,N_g))
  ALLOCATE(fg_edgH_yy_old(N_x,N_y+1,N_g))
  ALLOCATE(fg_edgH_xy_old(N_x,N_y+1,N_g))

  !Multigroup rad fluxes
  ALLOCATE(Fxg_edgV(N_x+1,N_y,N_g))
  ALLOCATE(Fyg_edgH(N_x,N_y+1,N_g))

  !Multigroup rad fluxes (last time step)
  ALLOCATE(Fxg_edgV_old(N_x+1,N_y,N_g))
  ALLOCATE(Fyg_edgH_old(N_x,N_y+1,N_g))

  !Multigroup rad energy densities
  ALLOCATE(Eg_avg(N_x,N_y,N_g))
  ALLOCATE(Eg_edgV(N_x+1,N_y,N_g))
  ALLOCATE(Eg_edgH(N_x,N_y+1,N_g))

  !Multigroup rad energy densities (last time step)
  ALLOCATE(Eg_avg_old(N_x,N_y,N_g))
  ALLOCATE(Eg_edgV_old(N_x+1,N_y,N_g))
  ALLOCATE(Eg_edgH_old(N_x,N_y+1,N_g))

  !MGQD boundary factors
  ALLOCATE(Cg_L(N_y,N_g))
  ALLOCATE(Cg_B(N_x,N_g))
  ALLOCATE(Cg_R(N_y,N_g))
  ALLOCATE(Cg_T(N_x,N_g))

  !Incoming multigroup rad energy densities
  ALLOCATE(Eg_in_L(N_y,N_g))
  ALLOCATE(Eg_in_B(N_x,N_g))
  ALLOCATE(Eg_in_R(N_y,N_g))
  ALLOCATE(Eg_in_T(N_x,N_g))

  !Incoming multigroup rad fluxes
  ALLOCATE(Fg_in_L(N_y,N_g))
  ALLOCATE(Fg_in_B(N_x,N_g))
  ALLOCATE(Fg_in_R(N_y,N_g))
  ALLOCATE(Fg_in_T(N_x,N_g))

  !Total rad fluxes
  ALLOCATE(MGQD_Fx_edgV(N_x+1,N_y))
  ALLOCATE(MGQD_Fy_edgH(N_x,N_y+1))

  !Total rad energy densities
  ALLOCATE(MGQD_E_avg(N_x,N_y))
  ALLOCATE(MGQD_E_edgV(N_x+1,N_y))
  ALLOCATE(MGQD_E_edgH(N_x,N_y+1))

  !Residuals for the transport solve
  ! ALLOCATE(MGQD_Residual(N_x,N_y,N_g,5,2,maxit_MLOQD,maxit_RTE))

  ALLOCATE(G_old(N_x,N_y,N_g),Pold_L(N_x,N_y,N_g),Pold_B(N_x,N_y,N_g),Pold_R(N_x,N_y,N_g),Pold_T(N_x,N_y,N_g))

  MGQD_E_avg = 0d0
  MGQD_E_edgV = 0d0
  MGQD_E_edgH = 0d0
  DO g=1,N_g
    !black body radiation distribution at initial Temp
    bg = Bg_planck_calc(Tini,nu_g(g),nu_g(g+1),comp_unit)
    Eg_avg(:,:,g) = 4d0*pi*bg/c
    Eg_edgV(:,:,g) = 4d0*pi*bg/c
    Eg_edgH(:,:,g) = 4d0*pi*bg/c
    MGQD_E_avg = MGQD_E_avg + 4d0*pi*bg/c
    MGQD_E_edgV = MGQD_E_edgV + 4d0*pi*bg/c
    MGQD_E_edgH = MGQD_E_edgH + 4d0*pi*bg/c
  END DO
  Fxg_edgV = 0d0
  Fyg_edgH = 0d0
  MGQD_Fx_edgV = 0d0
  MGQD_Fy_edgH = 0d0

  fg_avg_xx = 1d0/3d0
  fg_avg_yy = 1d0/3d0
  fg_edgV_xx = 1d0/3d0
  fg_edgV_xy = 0d0
  fg_edgH_yy = 1d0/3d0
  fg_edgH_xy = 0d0
  fg_avg_xx_old = 1d0/3d0
  fg_avg_yy_old = 1d0/3d0
  fg_edgV_xx_old = 1d0/3d0
  fg_edgV_xy_old = 0d0
  fg_edgH_yy_old = 1d0/3d0
  fg_edgH_xy_old = 0d0

  IF (BC_Type(1) .EQ. 0) THEN
    Cg_L=-0.5d0
  ELSE IF (BC_Type(1) .EQ. 1) THEN
    Cg_L=0d0
  END IF

  IF (BC_Type(2) .EQ. 0) THEN
    Cg_B=-0.5d0
  ELSE IF (BC_Type(2) .EQ. 1) THEN
    Cg_B=0d0
  END IF

  IF (BC_Type(3) .EQ. 0) THEN
    Cg_R=0.5d0
  ELSE IF (BC_Type(3) .EQ. 1) THEN
    Cg_R=0d0
  END IF

  IF (BC_Type(4) .EQ. 0) THEN
    Cg_T=0.5d0
  ELSE IF (BC_Type(4) .EQ. 1) THEN
    Cg_T=0d0
  END IF

  CALL MGQD_In_Calc(Eg_in_L,Eg_in_B,Eg_in_R,Eg_in_T,Fg_in_L,Fg_in_B,Fg_in_R,Fg_in_T,I_edgV,I_edgH,Omega_x,Omega_y,&
    quad_weight,c,Comp_Unit,BC_Type,Open_Threads)

END SUBROUTINE MGQD_INIT

!==================================================================================================================================!
!
!==================================================================================================================================!

END MODULE INITIALIZERS
