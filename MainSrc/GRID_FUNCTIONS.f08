MODULE GRID_FUNCTIONS

  IMPLICIT NONE

CONTAINS

!==================================================================================================================================!
!
!==================================================================================================================================!
SUBROUTINE GRIDMAP_GEN_AVG(GridMap,Delx,Dely,N_x,N_y)
  REAL*8,INTENT(OUT):: GridMap(*)
  REAL*8,INTENT(IN):: Delx(*), Dely(*)
  INTEGER,INTENT(IN):: N_x, N_y
  REAL*8:: xloc(N_x), yloc
  INTEGER:: i, j, p

  xloc(1) = Delx(1)/2d0
  DO i=2,N_x
    xloc(i) = xloc(i-1) + (Delx(i-1) + Delx(i))/2d0
  END DO

  p = 0
  yloc = 0d0
  DO j=1,N_y
    yloc = yloc + Dely(j)/2d0
    DO i=1,N_x
      p = p + 1
      GridMap(p) = xloc(i)
      p = p + 1
      GridMap(p) = yloc
    END DO
    yloc = yloc + Dely(j)/2d0
  END DO

END SUBROUTINE GRIDMAP_GEN_AVG

!==================================================================================================================================!
!
!==================================================================================================================================!
SUBROUTINE GRIDMAP_GEN_EDGV(GridMap,Delx,Dely,N_x,N_y)
  REAL*8,INTENT(OUT):: GridMap(*)
  REAL*8,INTENT(IN):: Delx(*), Dely(*)
  INTEGER,INTENT(IN):: N_x, N_y
  REAL*8:: xloc(N_x), yloc
  INTEGER:: i, j, p

  xloc(1) = 0d0
  DO i=2,N_x
    xloc(i) = xloc(i-1) + Delx(i-1)
  END DO

  p = 0
  yloc = 0d0
  DO j=1,N_y
    yloc = yloc + Dely(j)/2d0
    DO i=1,N_x
      p = p + 1
      GridMap(p) = xloc(i)
      p = p + 1
      GridMap(p) = yloc
    END DO
    yloc = yloc + Dely(j)/2d0
  END DO

END SUBROUTINE GRIDMAP_GEN_EDGV

!==================================================================================================================================!
!
!==================================================================================================================================!
SUBROUTINE GRIDMAP_GEN_EDGH(GridMap,Delx,Dely,N_x,N_y)
  REAL*8,INTENT(OUT):: GridMap(*)
  REAL*8,INTENT(IN):: Delx(*), Dely(*)
  INTEGER,INTENT(IN):: N_x, N_y
  REAL*8:: xloc(N_x), yloc
  INTEGER:: i, j, p

  xloc(1) = Delx(1)/2d0
  DO i=2,N_x
    xloc(i) = xloc(i-1) + (Delx(i-1) + Delx(i))/2d0
  END DO

  p = 0
  yloc = 0d0
  DO j=1,N_y
    DO i=1,N_x
      p = p + 1
      GridMap(p) = xloc(i)
      p = p + 1
      GridMap(p) = yloc
    END DO
    yloc = yloc + Dely(j)
  END DO

END SUBROUTINE GRIDMAP_GEN_EDGH

!==================================================================================================================================!
! SUBROUTINE MAP_GRIDS
!
!DEF:
!  Finds a map between two grids on the x,y plane
!  Maps Grid1 -> Grid2
!
!  Let 1pts be the number of points on Grid1 and 2pts the number of points on Grid2
!  Grid1 has dimension= 2*(1pts) . . . Grid2 has dimension= 2*(2pts)
!  Map has dimension= 4*(1pts)
!  For i=1,2,...,1pts: Grid1(2*i-1 : 2*i) = (x1_i, y1_i) where x1_i and y1_i are the x,y coordinates of the i-th point on Grid1
!  For j=1,2,...,2pts: Grid2(2*j-1 : 2*j) = (x2_j, y2_j) where x2_j and y2_j are the x,y coordinates of the j-th point on Grid2
!  For k=1,2,...,2pts: Map(4*k-3 : 4*k) = (p1, p2, p3, p4) where p1-4 are the locations(i) in Grid1 of the x-coordinates of the closest four points from Grid1 to point k in Grid2
!  -- Example: the closest point in Grid1 to point k in Grid2 has coordinates (x,y) = (Grid1(Map(4*k-3)), Grid1(Map(4*k-3)+1))
!
!INPUT:
!
!OUTPUT:
!
!==================================================================================================================================!
SUBROUTINE MAP_GRIDS(Map,Grid1,Grid2,len1,len2,Vec_Locs)
  INTEGER,INTENT(OUT):: Map(*)
  INTEGER,INTENT(OUT),OPTIONAL:: Vec_Locs(*)
  REAL*8,INTENT(IN):: Grid1(*), Grid2(*)
  INTEGER,INTENT(IN):: len1, len2
  REAL*8:: d(4), r
  INTEGER:: l(4), i, j, mp

  mp = 0
  G2: DO j=1,len2,2 !Moving through the 'mapped-to' grid
    d = HUGE(0d0) !Presetting distance to point j as the largest number possible

    G1: DO i=1,len1,2 !Moving through the 'mapped-from' grid

      r = DIST(Grid1(j),Grid1(j+1),Grid2(i),Grid2(i+1)) !Calculating distance between points j & i
      IF (r .LT. d(1)) THEN !Found closest i to j (so far)
        !Shuffling last closest points down a level
        d(4) = d(3)
        d(3) = d(2)
        d(2) = d(1)
        !Recording location (i) and distance (r)
        d(1) = r
        l(1) = i
      ELSE IF (r .LT. d(2)) THEN !Found second closest i to j (so far)
        !Shuffling last closest points down a level
        d(4) = d(3)
        d(3) = d(2)
        !Recording location (i) and distance (r)
        d(2) = r
        l(2) = i
      ELSE IF (r .LT. d(3)) THEN !Found third closest i to j (so far)
        !Shuffling last closest points down a level
        d(4) = d(3)
        !Recording location (i) and distance (r)
        d(3) = r
        l(3) = i
      ELSE IF (r .LT. d(4)) THEN !Found fourth closest i to j (so far)
        !Recording location (i) and distance (r)
        d(4) = r
        l(4) = i
      END IF

    END DO G1

    !Recording closest four points from grid2 to point i on grid1
    Map(mp+1) = l(1)
    Map(mp+2) = l(2)
    Map(mp+3) = l(3)
    Map(mp+4) = l(4)
    mp = mp + 4 !Moving to next map position

  END DO G2

  IF (PRESENT(Vec_Locs)) THEN
    i=len2*2
    DO j=1,i
      Vec_Locs(j) = (Map(j)+1)/2
    END DO
  END IF

END SUBROUTINE MAP_GRIDS

!==================================================================================================================================!
! FUNCTION DIST(xp1, yp1, xp2, yp2)
!
!DEF:
!  Calculates the Euclidian distance between two points on the x,y plane
!
!INPUT:
!  xp1 (REAL*8) - x coordinate of point 1
!  yp1 (REAL*8) - y coordinate of point 1
!  xp2 (REAL*8) - x coordinate of point 2
!  xp2 (REAL*8) - y coordinate of point 2
!
!OUTPUT:
!  DIST (REAL*8) - Euclidian distance between points (xp1,yp1) and (xp2,yp2)
!
!==================================================================================================================================!
FUNCTION DIST(xp1, yp1, xp2, yp2)
  REAL*8:: DIST
  REAL*8,INTENT(IN):: xp1, yp1, xp2, yp2

  DIST = ABS( SQRT(xp1**2 + yp1**2) - SQRT(xp2**2 + yp2**2) )

END FUNCTION DIST

!==================================================================================================================================!
! FUNCTION BILINEAR_INTERPOLATE
!
!DEF:
!  Performs Bilinear Interpolation to find the value of a function f(x,y) given: f_11 = f(x1,y1), f_21 = f(x2,y1)
!                                                                                f_12 = f(x1,y2), f_22 = f(x2,y2)
!
!INPUT:
!
!OUTPUT:
!  BILINEAR_INTERPOLATE (REAL*8) - the function value f(x,y)
!
!==================================================================================================================================!
FUNCTION BILINEAR_INTERPOLATE(f_11,f_12,f_21,f_22,x,x1,x2,y,y1,y2)
  REAL*8:: BILINEAR_INTERPOLATE
  REAL*8,INTENT(IN):: x, x1, x2, y, y1, y2
  REAL*8,INTENT(IN):: f_11, f_12, f_21, f_22
  REAL*8:: xx1, xx2, yy1, yy2

  xx1 = x - x1
  xx2 = x2 - x
  yy1 = y - y1
  yy2 = y2 - y

  BILINEAR_INTERPOLATE = ( f_11*xx2*yy2 + f_21*xx1*yy2 + f_12*xx2*yy1 + f_22*xx1*yy1 )/( (x2 - x1)*(y2 - y1) )

END FUNCTION BILINEAR_INTERPOLATE

!==================================================================================================================================!
!
!==================================================================================================================================!

END MODULE GRID_FUNCTIONS
