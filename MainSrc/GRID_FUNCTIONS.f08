MODULE GRID_FUNCTIONS

  IMPLICIT NONE

CONTAINS

!==================================================================================================================================!
! Map f1 onto f2
!==================================================================================================================================!
SUBROUTINE FMAP(f1,f2,Grid1,Grid2,GMap,VMap,Glen1,Glen2,flen1,flen2)
  REAL*8,INTENT(OUT):: f2(*)
  REAL*8,INTENT(IN):: f1(*), Grid1(*), Grid2(*)
  INTEGER,INTENT(IN):: GMap(*), VMap(*), Glen1, Glen2, flen1, flen2
  INTEGER:: olen, i, j, p1, p2

  olen = flen2/Glen2

  p2 = 0
  DO i=1,olen
    p1 = (i-1)*Glen1
    DO j=1,Glen2

      p2 = p2 + 1
      f2(p2) = BILINEAR_INTERPOLATE(f1(p1 + Vmap(4*j-3)), f1(p1 + Vmap(4*j-1)), f1(p1 + Vmap(4*j-2)), f1(p1 + Vmap(4*j)), &
      Grid2(2*j-1), Grid1(GMap(4*j-3)), Grid1(GMap(4*j)), Grid2(2*j), Grid1(GMap(4*j-3)+1), Grid1(GMap(4*j)+1))

    END DO
  END DO

END SUBROUTINE FMAP

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
!
!==================================================================================================================================!
! SUBROUTINE GRIDMAP_GEN_BNDS(GridMap,Delx,Dely,N_x,N_y,xlen,ylen)
!   REAL*8,INTENT(OUT):: GridMap(*)
!   REAL*8,INTENT(IN):: Delx(*), Dely(*), xlen, ylen
!   INTEGER,INTENT(IN):: N_x, N_y
!   INTEGER:: p
!
!   p = 1
!   CALL GRIDMAP_GEN_LR_BND(GridMap(p),0d0,Dely,N_y)
!
!   p = p + 2*N_y
!   CALL GRIDMAP_GEN_BT_BND(GridMap(p),Delx,0d0,N_x)
!
!   p = p + 2*N_x
!   CALL GRIDMAP_GEN_LR_BND(GridMap(p),xlen,Dely,N_y)
!
!   p = p + 2*N_y
!   CALL GRIDMAP_GEN_BT_BND(GridMap(p),Delx,ylen,N_x)
!
! END SUBROUTINE GRIDMAP_GEN_BNDS

SUBROUTINE GRIDMAP_GEN_LR_BND(GridMap,x,Dely,N_y)
  REAL*8,INTENT(OUT):: GridMap(*)
  REAL*8,INTENT(IN):: x, Dely(*)
  INTEGER,INTENT(IN):: N_y
  REAL*8:: yloc
  INTEGER:: i, p

  p = 0
  yloc = 0d0
  DO i=1,N_y
    yloc = yloc + Dely(i)/2d0

    p = p + 1
    GridMap(p) = x
    p = p + 1
    GridMap(p) = yloc

    yloc = yloc + Dely(i)/2d0
  END DO

END SUBROUTINE GRIDMAP_GEN_LR_BND

SUBROUTINE GRIDMAP_GEN_BT_BND(GridMap,Delx,y,N_x)
  REAL*8,INTENT(OUT):: GridMap(*)
  REAL*8,INTENT(IN):: Delx(*), y
  INTEGER,INTENT(IN):: N_x
  REAL*8:: xloc
  INTEGER:: i, p

  p = 0
  xloc = 0d0
  DO i=1,N_x
    xloc = xloc + Delx(i)/2d0

    p = p + 1
    GridMap(p) = xloc
    p = p + 1
    GridMap(p) = y

    xloc = xloc + Delx(i)/2d0
  END DO

END SUBROUTINE GRIDMAP_GEN_BT_BND

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
SUBROUTINE MAP_GRIDS(Map,Grid1,Grid2,len1,len2,xlen1,Vec_Locs)
  INTEGER,INTENT(OUT):: Map(*)
  INTEGER,INTENT(OUT),OPTIONAL:: Vec_Locs(*)
  REAL*8,INTENT(IN):: Grid1(*), Grid2(*)
  INTEGER,INTENT(IN):: len1, len2, xlen1
  REAL*8:: d(4), r
  INTEGER:: l(4), i, j, mp, k

  mp = 0
  G2: DO j=1,len2,2 !Moving through the 'mapped-to' grid
    d = HUGE(0d0) !Presetting distance to point j as the largest number possible

    G1: DO i=1,len1,2 !Moving through the 'mapped-from' grid

      r = DIST(Grid1(i),Grid1(i+1),Grid2(j),Grid2(j+1)) !Calculating distance between points j & i
      IF (r .LT. d(1)) THEN !Found closest i to j (so far)
        !Recording location (i) and distance (r)
        d(1) = r
        l(1) = i
      END IF

    END DO G1

    !l(1) = bottom left
    !l(2) = bottom right
    !l(3) = top left
    !l(4) = top right

    Xcheck: IF ((Grid1(l(1)) - Grid2(j)) .LT. -1d-10) THEN
      LYcheck: IF ((Grid1(l(1)+1) - Grid2(j+1)) .LT. -1d-10) THEN !bottom left
        l(2) = l(1) + 2
        l(3) = l(1) + xlen1
        l(4) = l(3) + 2

      ELSE IF ((Grid1(l(1)+1) - Grid2(j+1)) .GT. 1d-10) THEN LYcheck !top left
        l(3) = l(1)
        l(1) = l(3) - xlen1
        l(2) = l(1) + 2
        l(4) = l(3) + 2

      ELSE LYcheck !exact left
        l(2) = l(1) + 2
        l(3) = l(1)
        l(4) = l(2)

      END IF LYcheck

    ELSE IF ((Grid1(l(1)) - Grid2(j)) .GT. 1d-10) THEN Xcheck
      RYcheck: IF ((Grid1(l(1)+1) - Grid2(j+1)) .LT. -1d-10) THEN !bottom right
        l(2) = l(1)
        l(1) = l(2) - 2
        l(3) = l(1) + xlen1
        l(4) = l(3) + 2

      ELSE IF ((Grid1(l(1)+1) - Grid2(j+1)) .GT. 1d-10) THEN RYcheck !top right
        l(4) = l(1)
        l(2) = l(4) - xlen1
        l(1) = l(2) - 2
        l(3) = l(4) - 2

      ELSE RYcheck !exact right
        l(2) = l(1)
        l(1) = l(2) - 2
        l(3) = l(1)
        l(4) = l(2)

      END IF RYcheck

    ELSE Xcheck
      Ycheck: IF ((Grid1(l(1)+1) - Grid2(j+1)) .LT. -1d-10) THEN !exact bottom
        l(2) = l(1)
        l(3) = l(1) + xlen1
        l(4) = l(3)

      ELSE IF ((Grid1(l(1)+1) - Grid2(j+1)) .GT. 1d-10) THEN Ycheck !exact top
        l(3) = l(1)
        l(1) = l(3) - xlen1
        l(2) = l(1)
        l(4) = l(3)

      ELSE Ycheck !exact point
        l(2) = l(1)
        l(3) = l(1)
        l(4) = l(1)

      END IF Ycheck

    END IF Xcheck

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

  DIST = ABS( SQRT((xp1-xp2)**2 + (yp1-yp2)**2) )

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


  IF ((x2 - x1) .EQ. 0) THEN !x-coordinates of all pts are the same
    IF ((y2 - y1) .EQ. 0) THEN !already have exact point, return this value
      BILINEAR_INTERPOLATE = f_11

    ELSE !only interpolate along y direction
      yy1 = y - y1
      yy2 = y2 - y

      BILINEAR_INTERPOLATE = ( f_11*yy2 + f_12*yy1 )/( (y2 - y1) )

    END IF

  ELSE IF ((y2 - y1) .EQ. 0) THEN !y-coordinates of all pts are the same, interpolate only along x-direction
    xx1 = x - x1
    xx2 = x2 - x

    BILINEAR_INTERPOLATE = ( f_11*xx2 + f_21*xx1)/( (x2 - x1) )

  ELSE !4 unique points, can perform bilinear interpolation across x and y directions
    xx1 = x - x1
    xx2 = x2 - x
    yy1 = y - y1
    yy2 = y2 - y

    BILINEAR_INTERPOLATE = ( f_11*xx2*yy2 + f_21*xx1*yy2 + f_12*xx2*yy1 + f_22*xx1*yy1 )/( (x2 - x1)*(y2 - y1) )

  END IF

END FUNCTION BILINEAR_INTERPOLATE

!==================================================================================================================================!
!
!==================================================================================================================================!
SUBROUTINE AVG_APND_BND(f_out,f_in,LR_Bnd,BT_Bnd,Crn_Bnd,N_g,N_y,N_x)
  REAL*8,INTENT(OUT):: f_out(*)
  REAL*8,INTENT(IN):: f_in(*), LR_Bnd(*), BT_Bnd(*), Crn_Bnd(*)
  INTEGER,INTENT(IN):: N_g, N_y, N_x

  INTEGER:: g, j, i
  INTEGER:: p_out, p_in, p_LR, p_BT, p_C

  p_out=0
  p_in=0
  p_LR=0
  p_BT=0
  p_C=0
  DO g=1,N_g

    !Bottom-Left boundary corner
    p_out = p_out + 1
    p_C = p_C + 1
    f_out(p_out) = Crn_Bnd(p_C)

    !Bottom boundary
    DO i=1,N_x

      p_out = p_out + 1
      p_BT = p_BT + 1

      f_out(p_out) = BT_Bnd(p_BT)

    END DO

    !Bottom-right boundary corner
    p_out = p_out + 1
    p_C = p_C + 1
    f_out(p_out) = Crn_Bnd(p_C)

    !
    DO j=1,N_y

      !Left boundary
      p_out = p_out + 1
      p_LR = p_LR + 1

      f_out(p_out) = LR_Bnd(p_LR)

      !Cell-centers
      DO i=1,N_x

        p_out = p_out + 1
        p_in = p_in + 1

        f_out(p_out) = f_in(p_in)

      END DO

      !Right boundary
      p_out = p_out + 1
      p_LR = p_LR + 1

      f_out(p_out) = LR_Bnd(p_LR)

    END DO

    !Top-left boundary corner
    p_out = p_out + 1
    p_C = p_C + 1
    f_out(p_out) = Crn_Bnd(p_C)

    !Top boundary
    DO i=1,N_x

      p_out = p_out + 1
      p_BT = p_BT + 1

      f_out(p_out) = BT_Bnd(p_BT)

    END DO

    !Top-right boundary corner
    p_out = p_out + 1
    p_C = p_C + 1
    f_out(p_out) = Crn_Bnd(p_C)
  END DO

END SUBROUTINE AVG_APND_BND

!==================================================================================================================================!
!
!==================================================================================================================================!
SUBROUTINE COLLECT_LR_BND(bnd,f,N_g,N_y,N_x)
  REAL*8,INTENT(OUT):: bnd(*)
  REAL*8,INTENT(IN):: f(*)
  INTEGER,INTENT(IN):: N_g, N_y, N_x
  INTEGER:: g, j, i, p1, p2

  p1=0
  p2=0
  DO g=1,N_g

    DO j=1,N_y

      p1 = p1 + 1
      p2 = p2 + 1
      bnd(p2) = f(p1)

      p1 = p1 + (N_x-1)
      p2 = p2 + 1
      bnd(p2) = f(p1)

    END DO

  END DO

END SUBROUTINE COLLECT_LR_BND

!==================================================================================================================================!
!
!==================================================================================================================================!
SUBROUTINE COLLECT_BT_BND(bnd,f,N_g,N_y,N_x)
  REAL*8,INTENT(OUT):: bnd(*)
  REAL*8,INTENT(IN):: f(*)
  INTEGER,INTENT(IN):: N_g, N_y, N_x
  INTEGER:: g, j, i, p1, p2

  p1=0
  p2=0
  DO g=1,N_g
    DO j=1,N_x

      p1 = p1 + 1
      p2 = p2 + 1

      bnd(p2) = f(p1)

    END DO

    p1 = p1 + (N_y-2)*N_x
    DO j=1,N_x

      p1 = p1 + 1
      p2 = p2 + 1

      bnd(p2) = f(p1)

    END DO
  END DO

END SUBROUTINE COLLECT_BT_BND

!==================================================================================================================================!
!
!==================================================================================================================================!
SUBROUTINE COLLECT_CRN_BND(bnd,f,N_g,N_y,N_x)
  REAL*8,INTENT(OUT):: bnd(*)
  REAL*8,INTENT(IN):: f(*)
  INTEGER,INTENT(IN):: N_g, N_y, N_x
  INTEGER:: g, j, i, p1, p2

  p1=0
  p2=0
  DO g=1,N_g
    p1 = p1 + 1
    p2 = p2 + 1
    bnd(p2) = f(p1)

    p1 = p1 + (N_x-1)
    p2 = p2 + 1
    bnd(p2) = f(p1)

    p1 = p1 + (N_y-2)*N_x + 1
    p2 = p2 + 1
    bnd(p2) = f(p1)

    p1 = p1 + (N_x-1)
    p2 = p2 + 1
    bnd(p2) = f(p1)
  END DO

END SUBROUTINE COLLECT_CRN_BND

!==================================================================================================================================!
!
!==================================================================================================================================!
SUBROUTINE EXTND_BT_BND(f_out,f_in,N_g,N_y,N_x)
  REAL*8,INTENT(OUT):: f_out(*)
  REAL*8,INTENT(IN):: f_in(*)
  INTEGER,INTENT(IN):: N_g, N_y, N_x
  INTEGER:: g, j, i, p1, p2

  p1=0
  p2=0
  DO g=1,N_g

    DO i=1,N_x
      p1 = p1 + 1
      p2 = p2 + 1
      f_out(p2) = f_in(p1)
    END DO

    p1 = p1 - N_x
    DO j=1,N_y

      DO i=1,N_x
        p1 = p1 + 1
        p2 = p2 + 1
        f_out(p2) = f_in(p1)
      END DO

    END DO

    p1 = p1 - N_x
    DO i=1,N_x
      p1 = p1 + 1
      p2 = p2 + 1
      f_out(p2) = f_in(p1)
    END DO

  END DO

END SUBROUTINE EXTND_BT_BND

!==================================================================================================================================!
!
!==================================================================================================================================!
SUBROUTINE EXTND_LR_BND(f_out,f_in,N_g,N_y,N_x)
  REAL*8,INTENT(OUT):: f_out(*)
  REAL*8,INTENT(IN):: f_in(*)
  INTEGER,INTENT(IN):: N_g, N_y, N_x
  INTEGER:: g, j, i, p1, p2

  p1=0
  p2=0
  DO g=1,N_g

    DO j=1,N_y
      p2 = p2 + 1
      f_out(p2) = f_in(p1+1)

      DO i=1,N_x
        p1 = p1 + 1
        p2 = p2 + 1
        f_out(p2) = f_in(p1)
      END DO

      p2 = p2 + 1
      f_out(p2) = f_in(p1)

    END DO

  END DO

END SUBROUTINE EXTND_LR_BND

END MODULE GRID_FUNCTIONS
