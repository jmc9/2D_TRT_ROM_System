MODULE WAVEPROP_TOOLS

  IMPLICIT NONE

! INTERFACE XWAVE_SPEED
!   MODULE PROCEDURE XWAVE_SPEED_edg, XWAVE_SPEED_cnt
! END INTERFACE

CONTAINS

!==================================================================================================================================!
!
!==================================================================================================================================!
FUNCTION YMEAN(x, dy, ylen, N_y, N_x, xp)
  REAL*8:: YMEAN
  REAL*8,INTENT(IN):: x(*), dy(*), ylen
  INTEGER,INTENT(IN):: N_y, N_x, xp
  INTEGER:: j, p

  YMEAN = 0d0
  p = xp
  DO j=1,N_y
    YMEAN = YMEAN + x(p) * dy(j)
    p = p + N_x
  END DO
  YMEAN = YMEAN / ylen
END FUNCTION YMEAN

!==================================================================================================================================!
! first-order finite difference derivative
!==================================================================================================================================!
FUNCTION FD_DV1(xl, xr, dx)
  REAL*8:: FD_DV1
  REAL*8,INTENT(IN):: xl, xr, dx

  FD_DV1 = (xr-xl)/dx
END FUNCTION FD_DV1

!==================================================================================================================================!
! uses cell-face values of x to determine x-derivatives at cell-centers
!==================================================================================================================================!
FUNCTION XWAVE_SPEED_edg(xc, xe, xold, dy, dx, dt, ylen, N_y, N_x)
  REAL*8:: XWAVE_SPEED_edg
  REAL*8,INTENT(IN):: xc(*), xe(*), xold(*)
  REAL*8,INTENT(IN):: dy(*), dx(*), dt, ylen
  INTEGER,INTENT(IN):: N_y, N_x
  REAL*8:: tsum, xsum
  REAL*8:: xl, xr, tl, tr
  INTEGER:: i

  tsum = 0d0
  xsum = 0d0

  xr = YMEAN(xe, dy, ylen, N_y, N_x, 1)
  DO i=1,N_x
    xl = xr
    xr = YMEAN(xe, dy, ylen, N_y, N_x, i+1)
    xsum = xsum + ABS( FD_DV1(xl, xr, dx(i)) )

    tl = YMEAN(xold, dy, ylen, N_y, N_x, i)
    tr = YMEAN(xc, dy, ylen, N_y, N_x, i)
    tsum = tsum + ABS( FD_DV1(tl, tr, dt) )
  END DO

  XWAVE_SPEED_edg = tsum/xsum
END FUNCTION XWAVE_SPEED_edg

!==================================================================================================================================!
! no cell-face values available, instead finds derivatives using central differences on cell centered values
!==================================================================================================================================!
FUNCTION XWAVE_SPEED_cnt(xc, xold, dy, dx, dt, ylen, N_y, N_x)
  REAL*8:: XWAVE_SPEED_cnt
  REAL*8,INTENT(IN):: xc(*), xold(*)
  REAL*8,INTENT(IN):: dy(*), dx(*), dt, ylen
  INTEGER,INTENT(IN):: N_y, N_x
  REAL*8:: tsum, xsum
  REAL*8:: xl, xr, tl, tr
  INTEGER:: i

  tsum = 0d0
  xsum = 0d0

  !left boundary - have to use forward difference formula
  xl = YMEAN(xc, dy, ylen, N_y, N_x, 1)
  xr = YMEAN(xc, dy, ylen, N_y, N_x, 2)
  xsum = xsum + ABS( FD_DV1(xl, xr, dx(1)) )
  !center cells can use central difference
  DO i=2,N_x-1
    xl = YMEAN(xc, dy, ylen, N_y, N_x, i-1)
    xr = YMEAN(xc, dy, ylen, N_y, N_x, i+1)
    xsum = xsum + ABS( FD_DV1(xl, xr, dx(i))/2d0 )
  END DO
  !right boundary - have to use backward difference
  xl = YMEAN(xc, dy, ylen, N_y, N_x, N_x-1)
  xr = YMEAN(xc, dy, ylen, N_y, N_x, N_x)
  xsum = xsum + ABS( FD_DV1(xl, xr, dx(N_x)) )

  DO i=1,N_x
    tl = YMEAN(xold, dy, ylen, N_y, N_x, i)
    tr = YMEAN(xc, dy, ylen, N_y, N_x, i)
    tsum = tsum + ABS( FD_DV1(tl, tr, dt) )
  END DO

  XWAVE_SPEED_cnt = tsum/xsum
END FUNCTION XWAVE_SPEED_cnt

!==================================================================================================================================!
!
!==================================================================================================================================!


END MODULE WAVEPROP_TOOLS
