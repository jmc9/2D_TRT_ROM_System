MODULE UPDATES

  USE TEMP_FUNCTIONS

  IMPLICIT NONE

CONTAINS

!==================================================================================================================================!
!
!==================================================================================================================================!
SUBROUTINE RT_BC_UPDATE(BC_Type,bcT_left,bcT_lower,bcT_right,bcT_upper,Comp_Unit,Nu_g,I_edgV,I_edgH,Ic_edgV,Ic_edgH)
  REAL*8,INTENT(IN):: bcT_left, bcT_lower, bcT_right, bcT_upper
  REAL*8,INTENT(IN):: Comp_Unit, Nu_g(:)
  INTEGER,INTENT(IN):: BC_Type(:)

  REAL*8,INTENT(OUT):: I_edgV(:,:,:,:), I_edgH(:,:,:,:), Ic_edgV(:,:,:,:), Ic_edgH(:,:,:,:)

  INTEGER:: m, g, N_x, N_y, N_m, N_g

  N_x = SIZE(I_edgH,1)
  N_y = SIZE(I_edgV,2)
  N_m = SIZE(I_edgV,3)
  N_g = SIZE(I_edgV,4)

  !left boundary condition for Omega_x>0
  IF (BC_Type(1) .EQ. 0) THEN
    DO g=1,N_g
      I_edgV(1,:,1:N_m/4,g) = Bg_planck_calc(bcT_left,nu_g(g),nu_g(g+1),comp_unit)
      I_edgV(1,:,3*N_m/4+1:N_m,g) = Bg_planck_calc(bcT_left,nu_g(g),nu_g(g+1),comp_unit)
      Ic_edgV(1,:,1:N_m/4,g) = Bg_planck_calc(bcT_left,nu_g(g),nu_g(g+1),comp_unit)
      Ic_edgV(1,:,3*N_m/4+1:N_m,g) = Bg_planck_calc(bcT_left,nu_g(g),nu_g(g+1),comp_unit)
    END DO
  ELSE IF (BC_Type(1) .EQ. 1) THEN
    DO g=1,N_g
      DO m=1,N_m/4
        I_edgV(1,:,m,g) = I_edgV(1,:,m+N_m/4,g)
        Ic_edgV(1,:,m,g) = Ic_edgV(1,:,m+N_m/4,g)
      END DO
      DO m=3*N_m/4+1,N_m
        I_edgV(1,:,m,g) = I_edgV(1,:,m-N_m/4,g)
        Ic_edgV(1,:,m,g) = Ic_edgV(1,:,m-N_m/4,g)
      END DO
    END DO
  END IF

  !bottom boundary condition for Omega_y>0
  IF (BC_Type(2) .EQ. 0) THEN
    DO g=1,N_g
      I_edgH(:,1,1:N_m/2,g) = Bg_planck_calc(bcT_lower,nu_g(g),nu_g(g+1),comp_unit)
      Ic_edgH(:,1,1:N_m/2,g) = Bg_planck_calc(bcT_lower,nu_g(g),nu_g(g+1),comp_unit)
    END DO
  ELSE IF (BC_Type(2) .EQ. 1) THEN
    DO g=1,N_g
      DO m=1,N_m/4
        I_edgH(:,1,m,g) = I_edgH(:,1,m+3*N_m/4,g)
        Ic_edgH(:,1,m,g) = Ic_edgH(:,1,m+3*N_m/4,g)
      END DO
      DO m=N_m/4+1,N_m/2
        I_edgH(:,1,m,g) = I_edgH(:,1,m+N_m/4,g)
        Ic_edgH(:,1,m,g) = Ic_edgH(:,1,m+N_m/4,g)
      END DO
    END DO
  END IF

  !right boundary condition for Omega_x<0
  IF (BC_Type(3) .EQ. 0) THEN
    DO g=1,N_g
      I_edgV(N_x+1,:,N_m/4+1:3*N_m/4,g) = Bg_planck_calc(bcT_right,nu_g(g),nu_g(g+1),comp_unit)
      Ic_edgV(N_x+1,:,N_m/4+1:3*N_m/4,g) = Bg_planck_calc(bcT_right,nu_g(g),nu_g(g+1),comp_unit)
    END DO
  ELSE IF (BC_Type(3) .EQ. 1) THEN
    DO g=1,N_g
      DO m=N_m/4+1,N_m/2
        I_edgV(N_x+1,:,m,g) = I_edgV(N_x+1,:,m-N_m/4,g)
        Ic_edgV(N_x+1,:,m,g) = Ic_edgV(N_x+1,:,m-N_m/4,g)
      END DO
      DO m=N_m/2+1,3*N_m/4
        I_edgV(N_x+1,:,m,g) = I_edgV(N_x+1,:,m+N_m/4,g)
        Ic_edgV(N_x+1,:,m,g) = Ic_edgV(N_x+1,:,m+N_m/4,g)
      END DO
    END DO
  END IF

  !top boundary condition for Omega_y<0
  IF (BC_Type(4) .EQ. 0) THEN
    DO g=1,N_g
      I_edgH(:,N_y+1,N_m/2+1:N_m,g) = Bg_planck_calc(bcT_upper,nu_g(g),nu_g(g+1),comp_unit)
      Ic_edgH(:,N_y+1,N_m/2+1:N_m,g) = Bg_planck_calc(bcT_upper,nu_g(g),nu_g(g+1),comp_unit)
    END DO
  ELSE IF (BC_Type(4) .EQ. 1) THEN
    DO g=1,N_g
      DO m=N_m/2+1,3*N_m/4
        I_edgH(:,N_y+1,m,g) = I_edgH(:,N_y+1,m-N_m/4,g)
        Ic_edgH(:,N_y+1,m,g) = Ic_edgH(:,N_y+1,m-N_m/4,g)
      END DO
      DO m=3*N_m/4+1,N_m
        I_edgH(:,N_y+1,m,g) = I_edgH(:,N_y+1,m-3*N_m/4,g)
        Ic_edgH(:,N_y+1,m,g) = Ic_edgH(:,N_y+1,m-3*N_m/4,g)
      END DO
    END DO
  END IF


END SUBROUTINE RT_BC_UPDATE

!==================================================================================================================================!
!
!==================================================================================================================================!
SUBROUTINE MATERIAL_UPDATE(RT_Src,MGQD_Src,KapE,KapB,KapR,Bg,Temp,Comp_Unit,Nu_g,Open_Threads)
  REAL*8,INTENT(OUT):: RT_Src(:,:,:,:), MGQD_Src(:,:,:)
  REAL*8,INTENT(OUT):: KapE(:,:,:), KapB(:,:,:), KapR(:,:,:), Bg(:,:,:)

  REAL*8,INTENT(IN):: Temp(:,:), comp_unit, nu_g(:)
  INTEGER,INTENT(IN):: Open_Threads

  REAL*8,PARAMETER:: pi=4d0*ATAN(1d0)
  INTEGER:: N_y, N_x, N_m, N_g, Threads
  INTEGER:: i, j, g, m

  N_x = SIZE(Temp,1)
  N_y = SIZE(Temp,2)
  N_m = SIZE(RT_Src,3)
  N_g = SIZE(Nu_g,1)-1

  !$ Threads = Open_Threads
  !$ IF (Threads .GT. N_g) Threads = N_g
  !$OMP PARALLEL DEFAULT(SHARED) NUM_THREADS(Threads) PRIVATE(g,j,i)
  !$OMP DO
  DO g=1,N_g
    DO j=1,N_y
      DO i=1,N_x
        Bg(i,j,g) = Bg_planck_calc(Temp(i,j),nu_g(g),nu_g(g+1),comp_unit)
        KapB(i,j,g) = kapB_calc(Temp(i,j),nu_g(g),nu_g(g+1))
        KapE(i,j,g) = KapB(i,j,g)
        KapR(i,j,g) = KapB(i,j,g)
        MGQD_Src(i,j,g) = 4d0*pi*KapB(i,j,g)*Bg(i,j,g)
      END DO
    END DO

    DO m=1,N_m
      DO j=1,N_y
        DO i=1,N_x
          RT_Src(i*2-1,j*2-1,m,g) = KapB(i,j,g)*Bg(i,j,g)
          RT_Src(i*2,j*2-1,m,g) = KapB(i,j,g)*Bg(i,j,g)
          RT_Src(i*2-1,j*2,m,g) = KapB(i,j,g)*Bg(i,j,g)
          RT_Src(i*2,j*2,m,g) = KapB(i,j,g)*Bg(i,j,g)
        END DO
      END DO
    END DO

  END DO
  !$OMP END DO
  !$OMP END PARALLEL

END SUBROUTINE

!==================================================================================================================================!
!
!==================================================================================================================================!
SUBROUTINE MEB_SOLVE(Temp,Eg_HO,Bg,KapE,Temp_old,Delt,cV,Kap0,Comp_Unit)
  REAL*8,INTENT(INOUT):: Temp(:,:)
  REAL*8,INTENT(IN):: Eg_HO(:,:,:), Bg(:,:,:), KapE(:,:,:)
  REAL*8,INTENT(IN):: Temp_old(:,:), Delt, cV, Kap0, Comp_Unit

  REAL*8:: ar, sig_r
  REAL*8,ALLOCATABLE:: E(:,:), KapE_Bar(:,:)
  REAL*8,PARAMETER:: h=6.62613d-19    !erg*sh
  REAL*8,PARAMETER:: c=2.99792458d2   !cm/sh
  REAL*8,PARAMETER:: pi=4d0*ATAN(1d0)
  REAL*8,PARAMETER:: erg=6.24150934d11
  INTEGER:: i, j, g
  INTEGER:: N_x, N_y, N_g

  !DETERMINING ARRAY SIZES
  N_x = SIZE(Temp,1)
  N_y = SIZE(Temp,2)
  N_g = SIZE(Bg,3)
  sig_R=2d0*pi**5/(15d0*c**2*h**3*erg**4*Comp_Unit)
  aR=4d0*sig_R/c !1/(erg**3 cm**3)

  ALLOCATE(E(N_x,N_y),KapE_Bar(N_x,N_y))

  E = 0d0
  KapE_Bar = 0d0
  DO j=1,N_y
    DO i=1,N_x
      DO g=1,N_g
        E(i,j) = E(i,j) + Eg_HO(i,j,g)
        KapE_Bar(i,j) = KapE_Bar(i,j) + KapE(i,j,g)*Eg_HO(i,j,g)
      END DO
      KapE_Bar(i,j) = KapE_Bar(i,j)/E(i,j)
    END DO
  END DO

  DO j=1,N_y
    DO i=1,N_x
      Temp(i,j) = (c**2*h**3)/(c**2*h**3*cv*Comp_Unit+8d0*pi*Delt*Kap0/erg**4)*(Delt*c*KapE_Bar(i,j)*E(i,j) &
      + cv*Temp_old(i,j))*Comp_Unit
    END DO
  END DO

END SUBROUTINE

!==================================================================================================================================!
!
!==================================================================================================================================!

END MODULE UPDATES
