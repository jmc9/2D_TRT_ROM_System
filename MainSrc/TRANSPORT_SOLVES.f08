MODULE TRANSPORT_SOLVES

  USE LA_TOOLS

  IMPLICIT NONE

CONTAINS

!==================================================================================================================================!
!Subroutine TRANSPORT_SCB
!
! Solves the 2D radiative transfer equation discretized with simple corner balance
!==================================================================================================================================!
SUBROUTINE TRANSPORT_SCB(I_avg,I_edgV,I_edgH,I_crn,Ic_edgV,Ic_edgH,RT_Residual,Omega_x,Omega_y,Delx,Dely,A,KapE_in,Src_in,&
  I_crn_old,c,Delt,Open_Threads,Res_Calc)
  REAL*8,INTENT(IN):: Omega_x(:), Omega_y(:)
  REAL*8,INTENT(IN):: Delx(:), Dely(:), A(:,:)
  REAL*8,INTENT(IN):: KapE_in(:,:,:), Src_in(:,:,:,:)
  REAL*8,INTENT(IN):: I_crn_old(:,:,:,:)
  REAL*8,INTENT(IN):: c, Delt
  INTEGER,INTENT(IN):: Open_Threads
  LOGICAL,INTENT(IN):: Res_Calc

  REAL*8,INTENT(OUT):: I_avg(:,:,:,:), I_edgV(:,:,:,:), I_edgH(:,:,:,:)
  REAL*8,INTENT(OUT):: I_crn(:,:,:,:), Ic_edgV(:,:,:,:), Ic_edgH(:,:,:,:)
  REAL*8,INTENT(OUT):: RT_Residual(:,:,:,:,:)

  REAL*8,ALLOCATABLE:: Src(:,:,:,:), KapE(:,:,:)
  REAL*8,ALLOCATABLE:: B(:,:), Q(:)
  INTEGER:: N_g, N_m, N_x, N_y, Threads
  INTEGER:: left, right, bot, top
  INTEGER:: i, j, g, m
  INTEGER:: m1, m2

  ALLOCATE(B(4,4))
  ALLOCATE(Q(4))

  N_m = SIZE(Omega_x,1)
  N_g = SIZE(I_avg,4)
  N_x = SIZE(Delx,1)
  N_y = SIZE(Dely,1)

  ALLOCATE(Src(N_x*2,N_y*2,N_m,N_g))
  ALLOCATE(KapE(N_x,N_y,N_g))
  KapE = KapE_in + 1d0/(c*Delt)
  Src = Src_in + I_crn_old/(c*Delt)

  !$ threads = open_threads
  !$ IF (threads .GT. N_g) threads = N_g
  !$OMP PARALLEL DEFAULT(SHARED) NUM_THREADS(threads) &
  !$OMP& PRIVATE(m,i,j,m1,m2,left,right,bot,top,Q,B)
  !$OMP DO
  DO g=1,N_g
    !--------------------------------------------------!
    !            Omega_x > 0, Omega_y > 0              !
    !--------------------------------------------------!
    m1 = 1
    m2 = N_m/4
    DO m=m1,m2
      DO j=1,N_y
        DO i=1,N_x

          !corner positions
          left = 2*i-1
          right = 2*i
          bot = 2*j-1
          top = 2*j

          !build linear system
          CALL SCB_BGEN(B,Omega_x(m),Omega_y(m),Delx(i),Dely(j),KapE(i,j,g))
          Q(1) = A(i,j)*Src(left,bot,m,g)/4d0 + Omega_x(m)*Dely(j)*Ic_edgV(i,bot,m,g)/2d0 &
          + Omega_y(m)*Delx(i)*Ic_edgH(left,j,m,g)/2d0
          Q(2) = A(i,j)*Src(right,bot,m,g)/4d0 + Omega_y(m)*Delx(i)*Ic_edgH(right,j,m,g)/2d0
          Q(3) = A(i,j)*Src(right,top,m,g)/4d0
          Q(4) = A(i,j)*Src(left,top,m,g)/4d0 + Omega_x(m)*Dely(j)*Ic_edgV(i,top,m,g)/2d0

          !solve subcell system for corner intensities
          CALL LP_MATSOLVE(B,Q)
          I_crn(left,bot,m,g)  = Q(1)
          I_crn(right,bot,m,g) = Q(2)
          I_crn(right,top,m,g) = Q(3)
          I_crn(left,top,m,g)  = Q(4)

          !upwinding corner intensities to corner edges
          Ic_edgV(i+1,bot,m,g)   = I_crn(right,bot,m,g)
          Ic_edgV(i+1,top,m,g)   = I_crn(right,top,m,g)
          Ic_edgH(right,j+1,m,g) = I_crn(right,top,m,g)
          Ic_edgH(left,j+1,m,g)  = I_crn(left,top,m,g)

          !calculating forward cell edge intensities from corner edge intensities
          I_edgV(i+1,j,m,g) = (Ic_edgV(i+1,bot,m,g) + Ic_edgV(i+1,top,m,g))/2d0
          I_edgH(i,j+1,m,g) = (Ic_edgH(left,j+1,m,g) + Ic_edgH(right,j+1,m,g))/2d0

          !calculating cell average intenity from corner intensities
          I_avg(i,j,m,g) = (I_crn(left,bot,m,g) + I_crn(right,bot,m,g) + I_crn(right,top,m,g) + I_crn(left,top,m,g))/4d0

        END DO !end of x cell loop
      END DO !end of y cell loop
    END DO !end of angle loop

    !--------------------------------------------------!
    !            Omega_x < 0, Omega_y > 0              !
    !--------------------------------------------------!
    m1 = N_m/4 + 1
    m2 = N_m/2
    DO m=m1,m2
      DO j=1,N_y
        DO i=N_x,1,-1

          !corner positions
          left = 2*i-1
          right = 2*i
          bot = 2*j-1
          top = 2*j

          !build linear system
          CALL SCB_BGEN(B,Omega_x(m),Omega_y(m),Delx(i),Dely(j),KapE(i,j,g))
          Q(1) = A(i,j)*Src(left,bot,m,g)/4d0 + Omega_y(m)*Delx(i)*Ic_edgH(left,j,m,g)/2d0
          Q(2) = A(i,j)*Src(right,bot,m,g)/4d0 - Omega_x(m)*Dely(j)*Ic_edgV(i+1,bot,m,g)/2d0 &
          + Omega_y(m)*Delx(i)*Ic_edgH(right,j,m,g)/2d0
          Q(3) = A(i,j)*Src(right,top,m,g)/4d0 - Omega_x(m)*Dely(j)*Ic_edgV(i+1,top,m,g)/2d0
          Q(4) = A(i,j)*Src(left,top,m,g)/4d0

          !solve subcell system for corner intensities
          CALL LP_MATSOLVE(B,Q)
          I_crn(left,bot,m,g)  = Q(1)
          I_crn(right,bot,m,g) = Q(2)
          I_crn(right,top,m,g) = Q(3)
          I_crn(left,top,m,g)  = Q(4)

          !upwinding corner intensities to corner edges
          Ic_edgV(i,bot,m,g)   = I_crn(left,bot,m,g)
          Ic_edgV(i,top,m,g)   = I_crn(left,top,m,g)
          Ic_edgH(left,j+1,m,g)  = I_crn(left,top,m,g)
          Ic_edgH(right,j+1,m,g) = I_crn(right,top,m,g)

          !calculating forward cell edge intensities from corner edge intensities
          I_edgV(i,j,m,g) = (Ic_edgV(i,bot,m,g) + Ic_edgV(i,top,m,g))/2d0
          I_edgH(i,j+1,m,g) = (Ic_edgH(left,j+1,m,g) + Ic_edgH(right,j+1,m,g))/2d0

          !calculating cell average intenity from corner intensities
          I_avg(i,j,m,g) = (I_crn(left,bot,m,g) + I_crn(right,bot,m,g) + I_crn(right,top,m,g) + I_crn(left,top,m,g))/4d0

        END DO !end of x cell loop
      END DO !end of y cell loop
    END DO !end of angle loop

    !--------------------------------------------------!
    !            Omega_x < 0, Omega_y < 0              !
    !--------------------------------------------------!
    m1 = N_m/2 + 1
    m2 = 3*N_m/4
    DO m=m1,m2
      DO j=N_y,1,-1
        DO i=N_x,1,-1

          !corner positions
          left = 2*i-1
          right = 2*i
          bot = 2*j-1
          top = 2*j

          !build linear system
          CALL SCB_BGEN(B,Omega_x(m),Omega_y(m),Delx(i),Dely(j),KapE(i,j,g))
          Q(1) = A(i,j)*Src(left,bot,m,g)/4d0
          Q(2) = A(i,j)*Src(right,bot,m,g)/4d0 - Omega_x(m)*Dely(j)*Ic_edgV(i+1,bot,m,g)/2d0
          Q(3) = A(i,j)*Src(right,top,m,g)/4d0 - Omega_x(m)*Dely(j)*Ic_edgV(i+1,top,m,g)/2d0 &
          - Omega_y(m)*Delx(i)*Ic_edgH(right,j+1,m,g)/2d0
          Q(4) = A(i,j)*Src(left,top,m,g)/4d0 - Omega_y(m)*Delx(i)*Ic_edgH(left,j+1,m,g)/2d0

          !solve subcell system for corner intensities
          CALL LP_MATSOLVE(B,Q)
          I_crn(left,bot,m,g)  = Q(1)
          I_crn(right,bot,m,g) = Q(2)
          I_crn(right,top,m,g) = Q(3)
          I_crn(left,top,m,g)  = Q(4)

          !upwinding corner intensities to corner edges
          Ic_edgV(i,top,m,g)   = I_crn(left,top,m,g)
          Ic_edgV(i,bot,m,g)   = I_crn(left,bot,m,g)
          Ic_edgH(left,j,m,g)  = I_crn(left,bot,m,g)
          Ic_edgH(right,j,m,g) = I_crn(right,bot,m,g)

          !calculating forward cell edge intensities from corner edge intensities
          I_edgV(i,j,m,g) = (Ic_edgV(i,bot,m,g) + Ic_edgV(i,top,m,g))/2d0
          I_edgH(i,j,m,g) = (Ic_edgH(left,j,m,g) + Ic_edgH(right,j,m,g))/2d0

          !calculating cell average intenity from corner intensities
          I_avg(i,j,m,g) = (I_crn(left,bot,m,g) + I_crn(right,bot,m,g) + I_crn(right,top,m,g) + I_crn(left,top,m,g))/4d0

        END DO !end of x cell loop
      END DO !end of y cell loop
    END DO !end of angle loop

    !--------------------------------------------------!
    !            Omega_x > 0, Omega_y < 0              !
    !--------------------------------------------------!
    m1 = 3*N_m/4 + 1
    m2 = N_m
    DO m=m1,m2
      DO j=N_y,1,-1
        DO i=1,N_x

          !corner positions
          left = 2*i-1
          right = 2*i
          bot = 2*j-1
          top = 2*j

          !build linear system
          CALL SCB_BGEN(B,Omega_x(m),Omega_y(m),Delx(i),Dely(j),KapE(i,j,g))
          Q(1) = A(i,j)*Src(left,bot,m,g)/4d0 + Omega_x(m)*Dely(j)*Ic_edgV(i,bot,m,g)/2d0
          Q(2) = A(i,j)*Src(right,bot,m,g)/4d0
          Q(3) = A(i,j)*Src(right,top,m,g)/4d0 - Omega_y(m)*Delx(i)*Ic_edgH(right,j+1,m,g)/2d0
          Q(4) = A(i,j)*Src(left,top,m,g)/4d0 + Omega_x(m)*Dely(j)*Ic_edgV(i,top,m,g)/2d0 &
          - Omega_y(m)*Delx(i)*Ic_edgH(left,j+1,m,g)/2d0

          !solve subcell system for corner intensities
          CALL LP_MATSOLVE(B,Q)
          I_crn(left,bot,m,g)  = Q(1)
          I_crn(right,bot,m,g) = Q(2)
          I_crn(right,top,m,g) = Q(3)
          I_crn(left,top,m,g)  = Q(4)

          !upwinding corner intensities to corner edges
          Ic_edgV(i+1,top,m,g) = I_crn(right,top,m,g)
          Ic_edgV(i+1,bot,m,g) = I_crn(right,bot,m,g)
          Ic_edgH(right,j,m,g) = I_crn(right,bot,m,g)
          Ic_edgH(left,j,m,g)  = I_crn(left,bot,m,g)

          !calculating forward cell edge intensities from corner edge intensities
          I_edgV(i+1,j,m,g) = (Ic_edgV(i+1,bot,m,g) + Ic_edgV(i+1,top,m,g))/2d0
          I_edgH(i,j,m,g) = (Ic_edgH(left,j,m,g) + Ic_edgH(right,j,m,g))/2d0

          !calculating cell average intenity from corner intensities
          I_avg(i,j,m,g) = (I_crn(left,bot,m,g) + I_crn(right,bot,m,g) + I_crn(right,top,m,g) + I_crn(left,top,m,g))/4d0

        END DO !end of x cell loop
      END DO !end of y cell loop
    END DO !end of angle loop
  END DO
  !$OMP END DO

  !--------------------------------------------------!
  !              Residual Calculations               !
  !--------------------------------------------------!
  IF (Res_Calc) THEN
  !$OMP DO
  DO g=1,N_g
    DO m=1,N_m
      DO j=1,N_y
        DO i=1,N_x

          !corner positions
          left = 2*i-1
          right = 2*i
          bot = 2*j-1
          top = 2*j

          !bottom left
          RT_Residual(left,bot,m,g,1) = (Omega_x(m)*Dely(j) + Omega_y(m)*Delx(i) + KapE(i,j,g)*A(i,j))*I_crn(left,bot,m,g)/4d0 &
          + Omega_x(m)*Dely(j)*I_crn(right,bot,m,g)/4d0 + Omega_y(m)*Delx(i)*I_crn(left,top,m,g)/4d0 &
          - Omega_x(m)*Dely(j)*Ic_edgV(i,bot,m,g)/2d0 - Omega_y(m)*Delx(i)*Ic_edgH(left,j,m,g)/2d0 &
          - A(i,j)*Src(left,bot,m,g)/4d0

          RT_Residual(left,bot,m,g,2) = RT_Residual(left,bot,m,g,1)/I_crn(left,bot,m,g)

          !bottom right
          RT_Residual(right,bot,m,g,1) = (-Omega_x(m)*Dely(j) + Omega_y(m)*Delx(i) + KapE(i,j,g)*A(i,j))*I_crn(right,bot,m,g)/4d0 &
          - Omega_x(m)*Dely(j)*I_crn(left,bot,m,g)/4d0 + Omega_y(m)*Delx(i)*I_crn(right,top,m,g)/4d0 &
          + Omega_x(m)*Dely(j)*Ic_edgV(i+1,bot,m,g)/2d0 - Omega_y(m)*Delx(i)*Ic_edgH(right,j,m,g)/2d0 &
          - A(i,j)*Src(right,bot,m,g)/4d0

          RT_Residual(right,bot,m,g,2) = RT_Residual(right,bot,m,g,1)/I_crn(right,bot,m,g)

          !top right
          RT_Residual(right,top,m,g,1) = (-Omega_x(m)*Dely(j) - Omega_y(m)*Delx(i) + KapE(i,j,g)*A(i,j))*I_crn(right,top,m,g)/4d0 &
          - Omega_x(m)*Dely(j)*I_crn(left,top,m,g)/4d0 - Omega_y(m)*Delx(i)*I_crn(right,bot,m,g)/4d0 &
          + Omega_x(m)*Dely(j)*Ic_edgV(i+1,top,m,g)/2d0 + Omega_y(m)*Delx(i)*Ic_edgH(right,j+1,m,g)/2d0 &
          - A(i,j)*Src(right,top,m,g)/4d0

          RT_Residual(right,top,m,g,2) = RT_Residual(right,top,m,g,1)/I_crn(right,top,m,g)

          !top left
          RT_Residual(left,top,m,g,1) = (Omega_x(m)*Dely(j) - Omega_y(m)*Delx(i) + KapE(i,j,g)*A(i,j))*I_crn(left,top,m,g)/4d0 &
          + Omega_x(m)*Dely(j)*I_crn(right,top,m,g)/4d0 - Omega_y(m)*Delx(i)*I_crn(left,bot,m,g)/4d0 &
          - Omega_x(m)*Dely(j)*Ic_edgV(i,top,m,g)/2d0 + Omega_y(m)*Delx(i)*Ic_edgH(left,j+1,m,g)/2d0 &
          - A(i,j)*Src(left,top,m,g)/4d0

          RT_Residual(left,top,m,g,2) = RT_Residual(left,top,m,g,1)/I_crn(left,top,m,g)

        END DO
      END DO
    END DO
  END DO
  !$OMP END DO
  END IF
  !$OMP END PARALLEL

  DEALLOCATE(B,Q)

END SUBROUTINE TRANSPORT_SCB

!==================================================================================================================================!
!
!==================================================================================================================================!
SUBROUTINE SCB_BGEN(B,Omega_x,Omega_y,Delx,Dely,KapE)
  REAL*8,INTENT(IN):: Omega_x, Omega_y
  REAL*8,INTENT(IN):: Delx, Dely
  REAL*8,INTENT(IN):: KapE

  REAL*8,INTENT(OUT):: B(:,:)

  REAL*8:: A
  REAL*8:: gamma, Del_Omx, Del_Omy

  A = Delx*Dely
  gamma = (ABS(Omega_x)*Dely + ABS(Omega_y)*Delx + KapE*A)/4d0
  Del_Omx = Omega_x*Dely/4d0
  Del_Omy = Omega_y*Delx/4d0

  B = 0d0

  !Diagonal
  B(1,1) = gamma
  B(2,2) = gamma
  B(3,3) = gamma
  B(4,4) = gamma

  !Upper Diagonal
  B(1,2) = Del_Omx
  B(2,3) = Del_Omy
  B(3,4) = -Del_Omx

  !Upper Corner
  B(1,4) = Del_Omy

  !Lower Diagonal
  B(2,1) = -Del_Omx
  B(3,2) = -Del_Omy
  B(4,3) = Del_Omx

  !Lower Corner
  B(4,1) = -Del_Omy

END SUBROUTINE SCB_BGEN

!==================================================================================================================================!
!Subroutine COLLAPSE_INTENSITIES
!
! 'collapses' radiation intensities into low-order multigroup and grey quantities:
! --> H is the second angular moment of I (tensor)
! --> F is the first angular moment of I (vector)
! --> E is the zeroth angular moment of I (scalar)
!==================================================================================================================================!
SUBROUTINE COLLAPSE_INTENSITIES(Open_Threads,I_avg,I_edgV,I_edgH,Omega_x,Omega_y,quad_weight,Comp_Unit,Hg_avg_xx,Hg_avg_xy,&
  Hg_avg_yy,Hg_edgV_xx,Hg_edgV_xy,Hg_edgH_yy,Hg_edgH_xy,Eg_edgV,Eg_edgH,Eg_avg,Fxg_edgV,Fyg_edgH,E_edgV,E_edgH,E_avg,Fx_edgV,&
  Fy_edgH)

  REAL*8,INTENT(IN):: I_avg(:,:,:,:)
  REAL*8,INTENT(IN):: I_edgV(:,:,:,:), I_edgH(:,:,:,:)
  REAL*8,INTENT(IN):: Omega_x(:), Omega_y(:), quad_weight(:)
  REAL*8,INTENT(IN):: Comp_Unit
  INTEGER,INTENT(IN):: Open_Threads

  REAL*8,INTENT(OUT):: Hg_avg_xx(:,:,:), Hg_avg_xy(:,:,:), Hg_avg_yy(:,:,:)
  REAL*8,INTENT(OUT):: Hg_edgV_xx(:,:,:), Hg_edgV_xy(:,:,:)
  REAL*8,INTENT(OUT):: Hg_edgH_yy(:,:,:), Hg_edgH_xy(:,:,:)
  REAL*8,INTENT(OUT):: Eg_avg(:,:,:), Eg_edgV(:,:,:), Eg_edgH(:,:,:)
  REAL*8,INTENT(OUT):: Fxg_edgV(:,:,:), Fyg_edgH(:,:,:)
  REAL*8,INTENT(OUT):: E_avg(:,:), E_edgV(:,:), E_edgH(:,:)
  REAL*8,INTENT(OUT):: Fx_edgV(:,:), Fy_edgH(:,:)

  INTEGER:: i, j, m, g
  INTEGER:: N_x, N_y, N_m, N_g
  INTEGER:: Threads
  REAL*8:: eps
  REAL*8,PARAMETER:: c=2.99792458d2   !cm/sh

  !setting base eps value and scaling this based on the computational units
  eps=1d-25
  eps=eps/comp_unit

  !Determining array sizes
  N_x = SIZE(I_avg,1)
  N_y = SIZE(I_avg,2)
  N_m = SIZE(I_avg,3)
  N_g = SIZE(I_avg,4)

  !Pre-setting arrays to zero
  Hg_avg_xx = 0d0
  Hg_avg_xy = 0d0
  Hg_avg_yy = 0d0
  Hg_edgV_xx = 0d0
  Hg_edgV_xy = 0d0
  Hg_edgH_yy = 0d0
  Hg_edgH_xy = 0d0
  Eg_avg = 0d0
  Eg_edgV = 0d0
  Eg_edgH = 0d0
  Fxg_edgV = 0d0
  Fyg_edgH = 0d0
  E_avg = 0d0
  E_edgV = 0d0
  E_edgH = 0d0
  Fx_edgV = 0d0
  Fy_edgH = 0d0

  !--------------------------------------------------!
  !             Multigroup E's and F's               !
  !--------------------------------------------------!
  !$ threads = open_threads
  !$ IF (threads .GT. N_g) threads = N_m
  !$OMP PARALLEL DEFAULT(SHARED) NUM_THREADS(threads) &
  !$OMP& PRIVATE(g,m,i,j)
  !$OMP DO
  DO g=1,N_g
    DO m=1,N_m

      !cell-averaged values
      DO j=1,N_y
        DO i=1,N_x
          Hg_avg_xx(i,j,g) = Hg_avg_xx(i,j,g) + Omega_x(m)**2*quad_weight(m)*(I_avg(i,j,m,g) + eps)
          Hg_avg_yy(i,j,g) = Hg_avg_yy(i,j,g) + Omega_y(m)**2*quad_weight(m)*(I_avg(i,j,m,g) + eps)
          Hg_avg_xy(i,j,g) = Hg_avg_xy(i,j,g) + Omega_x(m)*Omega_y(m)*quad_weight(m)*(I_avg(i,j,m,g) + eps)
          Eg_avg(i,j,g) = Eg_avg(i,j,g) + quad_weight(m)*(I_avg(i,j,m,g) + eps)
        END DO
      END DO

      !verticle cell-edge values
      DO j=1,N_y
        DO i=1,N_x+1
          Hg_edgV_xx(i,j,g) = Hg_edgV_xx(i,j,g) + Omega_x(m)**2*quad_weight(m)*(I_edgV(i,j,m,g) + eps)
          Hg_edgV_xy(i,j,g) = Hg_edgV_xy(i,j,g) + Omega_x(m)*Omega_y(m)*quad_weight(m)*(I_edgV(i,j,m,g) + eps)
          Fxg_edgV(i,j,g) = Fxg_edgV(i,j,g) + Omega_x(m)*quad_weight(m)*(I_edgV(i,j,m,g) + eps)
          Eg_edgV(i,j,g) = Eg_edgV(i,j,g) + quad_weight(m)*(I_edgV(i,j,m,g) + eps)
        END DO
      END DO

      !horizontal cell-edge values
      DO j=1,N_y+1
        DO i=1,N_x
          Hg_edgH_yy(i,j,g) = Hg_edgH_yy(i,j,g) + Omega_y(m)**2*quad_weight(m)*(I_edgH(i,j,m,g) + eps)
          Hg_edgH_xy(i,j,g) = Hg_edgH_xy(i,j,g) + Omega_x(m)*Omega_y(m)*quad_weight(m)*(I_edgH(i,j,m,g) + eps)
          Fyg_edgH(i,j,g) = Fyg_edgH(i,j,g) + Omega_x(m)*quad_weight(m)*(I_edgH(i,j,m,g) + eps)
          Eg_edgH(i,j,g) = Eg_edgH(i,j,g) + quad_weight(m)*(I_edgH(i,j,m,g) + eps)
        END DO
      END DO

    END DO
  END DO
  !$OMP END DO
  !$OMP END PARALLEL
  Eg_avg = Eg_avg/c !Dividing all values by c to find radiation energy density (instead of scalar rad flux)
  Eg_edgV = Eg_edgV/c
  Eg_edgH = Eg_edgH/c

  !--------------------------------------------------!
  !             Multigroup E's and F's               !
  !--------------------------------------------------!
  !Total E is the sum of all Eg's
  !Total F is the sum of all Fg's

  !cell-averaged values
  DO j=1,N_y
    DO i=1,N_x
      E_avg(i,j) = SUM(Eg_avg(i,j,:))
    END DO
  END DO

  !verticle cell-edge values
  DO j=1,N_y
    DO i=1,N_x+1
      E_edgV(i,j) = SUM(Eg_edgV(i,j,:))
      Fx_edgV(i,j) = SUM(Fxg_edgV(i,j,:))
    END DO
  END DO

  !horizontal cell-edge values
  DO j=1,N_y+1
    DO i=1,N_x
      E_edgH(i,j) = SUM(Eg_edgH(i,j,:))
      Fy_edgH(i,j) = SUM(Fyg_edgH(i,j,:))
    END DO
  END DO

END SUBROUTINE COLLAPSE_INTENSITIES
!==================================================================================================================================!
!
!==================================================================================================================================!


END MODULE TRANSPORT_SOLVES
