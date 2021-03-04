MODULE DMD_ROUTINES

  USE netcdf
  USE NCDF_IO
  USE GRID_FUNCTIONS

  IMPLICIT NONE

CONTAINS

!==================================================================================================================================!
! SUBROUTINE DMD_RECONSTRUCT
!
!DEF:
!  Reconstructs data at a specific instant of time from a DMD expansion
!
!INPUT:
!  N_g, N_y, N_x - size of the grid in energy/ y-length/ x-length, respectively
!                - Note that these are the grid sizes from the original decomposition, *NOT* the grid for the current simulation
!
!  DMDgsum - flag to denote whether the given decomposition was performed over the entire phase space or groupwise for each qd factor
!
!  Time - NOTe this is the Time RELATIVE TO THE START TIME OF THE DECOMPOSITION FIT
!
!OUTPUT:
!
!==================================================================================================================================!
SUBROUTINE DMD_RECONSTRUCT(dat, L, B, W, C, rank, len, N_g, DMDgsum, Time, Open_Threads)
  REAL*8,INTENT(OUT):: dat(*)
  REAL*8,INTENT(IN):: Time
  REAL*8,INTENT(IN):: C(*)
  COMPLEX*16,INTENT(IN):: L(*), B(*), W(*)
  INTEGER,INTENT(IN):: DMDgsum, len, N_g, rank(*), Open_Threads
  INTEGER:: Threads
  COMPLEX*16:: dmd_exp
  COMPLEX*16,ALLOCATABLE:: dat_(:)
  INTEGER:: i, g, r, p1, p2, p3, glen

  !initializing data array to zero
  ALLOCATE(dat_(len))
  DO i=1,len
    dat_(i) = 0d0
  END DO

  !
  ! Case 1: groupwise decomposition of data
  !
  IF (DMDgsum .EQ. 0) THEN
    glen = len/N_g

    p2 = 0
    p3 = 0
    !$ Threads = Open_Threads
    !$OMP PARALLEL DEFAULT(SHARED) NUM_THREADS(Threads) PRIVATE(g, r, i, dmd_exp, p1, p2, p3)
    !$OMP DO
    DO g=1,N_g
      p3 = SUM(rank(1:g))-rank(1)
      p2 = p3*glen
      DO r=1,rank(g)
        p1 = (g-1)*glen
        dmd_exp = EXP(L(p3+r) * Time)

        DO i=1,glen
          p1 = p1 + 1
          p2 = p2 + 1

          dat_(p1) = dat_(p1) + B(p3+r) * W(p2) * dmd_exp

        END DO
      END DO

    END DO
    !$OMP END DO
    !$OMP END PARALLEL

  !
  ! Case 2: decomposition of data over entire phase space
  !
  ELSE IF (DMDgsum .EQ. 1) THEN

    !$ Threads = Open_Threads
    p1 = 0
    DO r=1,rank(1)
      dmd_exp = EXP(L(r) * Time)
      !$OMP PARALLEL DEFAULT(SHARED) NUM_THREADS(Threads) PRIVATE(i)
      !$OMP DO
      DO i=1,len

        dat_(i) = dat_(i) + B(r) * W(p1+i) * dmd_exp

      END DO
      !$OMP END DO
      !$OMP END PARALLEL
      p1 = p1 + len
    END DO

  END IF

  !Taking the real component of the calculated expansion as the solution (assuming the solution to be a real entitity)
  !$OMP PARALLEL DEFAULT(SHARED) NUM_THREADS(Threads) PRIVATE(i)
  !$OMP DO
  DO i=1,len
    dat(i) = REALPART(dat_(i)) + C(i)
  END DO
  !$OMP END DO
  !$OMP END PARALLEL
  DEALLOCATE(dat_)

END SUBROUTINE DMD_RECONSTRUCT

!==================================================================================================================================!
! SUBROUTINE DMD_RECONSTRUCT_fg
!
!DEF:
!  Reconstructs multigroup qd-factors at a specific time instant (fg(t)) as a DMD expansion
!
!INPUT:
!  dN_g, dN_y, dN_x - size of the decomposotion grid in energy/ y-length/ x-length, respectively
!                   - Note that these are the grid sizes from the original decomposition, *NOT* the grid for the current simulation
!
!  N_g, N_y, N_x - size of the simulation grid in energy/ y-length/ x-length, respectively
!                - Note that these are the grid sizes defined for the current simulation
!
!  DMDgsum - flag to denote whether the given decomposition was performed over the entire phase space or groupwise for each boundary factor
!
!  Time - NOTe this is the Time RELATIVE TO THE START TIME OF THE DECOMPOSITION FIT
!
!OUTPUT:
!
!==================================================================================================================================!
SUBROUTINE DMD_RECONSTRUCT_fg(fg_avg_xx, fg_avg_yy, fg_edgV_xx, fg_edgV_xy, fg_edgH_yy, fg_edgH_xy,&
  L_fg_avg_xx, B_fg_avg_xx, W_fg_avg_xx, C_fg_avg_xx, L_fg_edgV_xx, B_fg_edgV_xx, W_fg_edgV_xx, C_fg_edgV_xx,&
  L_fg_avg_yy, B_fg_avg_yy, W_fg_avg_yy, C_fg_avg_yy, L_fg_edgH_yy, B_fg_edgH_yy, W_fg_edgH_yy, C_fg_edgH_yy,&
  L_fg_edgV_xy, B_fg_edgV_xy, W_fg_edgV_xy, C_fg_edgV_xy, L_fg_edgH_xy, B_fg_edgH_xy, W_fg_edgH_xy, C_fg_edgH_xy,&
  rrank_fg_avg_xx, rrank_fg_edgV_xx, rrank_fg_avg_yy, rrank_fg_edgH_yy, rrank_fg_edgV_xy, rrank_fg_edgH_xy,&
  dN_x, dN_y, dN_g, N_x, N_y, N_g, Time, DMDgsum, Sim_Grid_Avg, Sim_Grid_EdgV, Sim_Grid_EdgH, Dat_Grid_Avg,&
  Dat_Grid_EdgV, Dat_Grid_EdgH, GMap_xyAvg, GMap_xyEdgV, GMap_xyEdgH, VMap_xyAvg, VMap_xyEdgV, VMap_xyEdgH, Open_Threads)

  REAL*8,INTENT(OUT):: fg_avg_xx(*), fg_avg_yy(*)
  REAL*8,INTENT(OUT):: fg_edgV_xx(*), fg_edgV_xy(*)
  REAL*8,INTENT(OUT):: fg_edgH_yy(*), fg_edgH_xy(*)

  REAL*8,INTENT(IN):: C_fg_avg_xx(*), C_fg_avg_yy(*)
  REAL*8,INTENT(IN):: C_fg_edgV_xx(*), C_fg_edgH_yy(*), C_fg_edgV_xy(*), C_fg_edgH_xy(*)
  COMPLEX*16,INTENT(IN):: L_fg_avg_xx(*),  B_fg_avg_xx(*),  W_fg_avg_xx(*)
  COMPLEX*16,INTENT(IN):: L_fg_edgV_xx(*), B_fg_edgV_xx(*), W_fg_edgV_xx(*)
  COMPLEX*16,INTENT(IN):: L_fg_avg_yy(*),  B_fg_avg_yy(*),  W_fg_avg_yy(*)
  COMPLEX*16,INTENT(IN):: L_fg_edgH_yy(*), B_fg_edgH_yy(*), W_fg_edgH_yy(*)
  COMPLEX*16,INTENT(IN):: L_fg_edgV_xy(*), B_fg_edgV_xy(*), W_fg_edgV_xy(*)
  COMPLEX*16,INTENT(IN):: L_fg_edgH_xy(*), B_fg_edgH_xy(*), W_fg_edgH_xy(*)
  INTEGER,INTENT(IN):: rrank_fg_avg_xx(*), rrank_fg_edgV_xx(*), rrank_fg_avg_yy(*), rrank_fg_edgH_yy(*)
  INTEGER,INTENT(IN):: rrank_fg_edgV_xy(*), rrank_fg_edgH_xy(*)
  INTEGER,INTENT(IN):: dN_x, dN_y, dN_g
  INTEGER,INTENT(IN):: N_x, N_y, N_g
  INTEGER,INTENT(IN):: DMDgsum, Open_Threads
  REAL*8,INTENT(IN):: Sim_Grid_Avg(*), Sim_Grid_EdgV(*), Sim_Grid_EdgH(*)
  REAL*8,INTENT(IN):: Dat_Grid_Avg(*), Dat_Grid_EdgV(*), Dat_Grid_EdgH(*)
  REAL*8,INTENT(IN):: Time
  INTEGER,INTENT(IN):: GMap_xyAvg(*), GMap_xyEdgV(*), GMap_xyEdgH(*)
  INTEGER,INTENT(IN):: VMap_xyAvg(*), VMap_xyEdgV(*), VMap_xyEdgH(*)

  REAL*8,ALLOCATABLE:: f(:), f2(:), xxbnd(:), yybnd(:), crnbnd(:), bnd2(:)
  INTEGER:: len, g, i, j, p1, p2, p3

  ALLOCATE(xxbnd((2*dN_y)*dN_g), yybnd((2*dN_x)*dN_g), crnbnd(4*dN_g))

  !----------------------------------------
  !
  ! VERTICAL (y=const) CELL FACE DATA
  !
  !----------------------------------------

  !----------vector lengths----------!

  len = (dN_x+1)*(dN_y+2)*dN_g !length of EdgV data vector with boundary cells above/below the grid (for mapping purposes)
  ALLOCATE(f2(len)) !allocating space to hold reconstituted DMD data, with appended boundary cells

  len = (dN_x+1)*dN_y*dN_g !length of EdgV-type data vector
  ALLOCATE(f(len)) !allocating space to hold reconstituted DMD data

  !----------fg_edgV_xx----------!

  !reconstructing fg_edgV_xx from decomposition (on the same grid it was decomposed on)
  CALL DMD_RECONSTRUCT(f, L_fg_edgV_xx, B_fg_edgV_xx, W_fg_edgV_xx, C_fg_edgV_xx, rrank_fg_edgV_xx,&
  len, N_g, DMDgsum, Time, Open_Threads)

  !storing the fg_xx boundary data (for the left/right boundaries only)
  CALL COLLECT_LR_BND(xxbnd, f, dN_g, dN_y, dN_x+1)

  !Appending on extra 'boundary data' above/below the grid for mapping purposes
  CALL EXTND_BT_BND(f2, f, dN_g, dN_y, dN_x+1)

  !mapping the vector of fg_edgV_xx data from the decomposition grid to the simulation grid
  CALL FMAP(f2, fg_edgV_xx, Dat_Grid_EdgV, Sim_Grid_EdgV, GMap_xyEdgV, VMap_xyEdgV, (dN_x+1)*(dN_y+2), (N_x+1)*N_y, &
  (dN_x+1)*(dN_y+2)*N_g, (N_x+1)*N_y*N_g, Open_Threads)

  !----------fg_edgV_xy----------!

  !reconstructing fg_edgV_xy from decomposition (on the same grid it was decomposed on)
  CALL DMD_RECONSTRUCT(f, L_fg_edgV_xy, B_fg_edgV_xy, W_fg_edgV_xy, C_fg_edgV_xy, rrank_fg_edgV_xy,&
  len, N_g, DMDgsum, Time, Open_Threads)

  !Appending on extra 'boundary data' above/below the grid for mapping purposes
  CALL EXTND_BT_BND(f2, f, dN_g, dN_y, dN_x+1)

  !mapping the vector of fg_edgV_xy data from the decomposition grid to the simulation grid
  CALL FMAP(f2, fg_edgV_xy, Dat_Grid_EdgV, Sim_Grid_EdgV, GMap_xyEdgV, VMap_xyEdgV, (dN_x+1)*(dN_y+2), (N_x+1)*N_y, &
  (dN_x+1)*(dN_y+2)*N_g, (N_x+1)*N_y*N_g, Open_Threads)

  DEALLOCATE(f,f2) !deallocating DMD data vectors

  !----------------------------------------
  !
  ! HORIZONTAL (x=const) CELL FACE DATA
  !
  !----------------------------------------

  !----------vector lengths----------!

  len = (dN_x+2)*(dN_y+1)*dN_g !length of EdgH data vector with boundary cells left/right of the grid (for mapping purposes)
  ALLOCATE(f2(len)) !allocating space to hold reconstituted DMD data, with appended boundary cells

  len = dN_x*(dN_y+1)*dN_g !length of EdgH-type data vector
  ALLOCATE(f(len)) !allocating space to hold reconstituted DMD data

  !----------fg_edgH_yy----------!

  !reconstructing fg_edgH_yy from decomposition (on the same grid it was decomposed on)
  CALL DMD_RECONSTRUCT(f, L_fg_edgH_yy, B_fg_edgH_yy, W_fg_edgH_yy, C_fg_edgH_yy, rrank_fg_edgH_yy,&
  len, N_g, DMDgsum, Time, Open_Threads)

  !storing the fg_yy boundary data (for the top/bottom boundaries only)
  CALL COLLECT_BT_BND(yybnd, f, dN_g, dN_y+1, dN_x)

  !Appending on extra 'boundary data' left/right of the grid for mapping purposes
  CALL EXTND_LR_BND(f2, f, dN_g, dN_y+1, dN_x)

  !mapping the vector of fg_edgH_yy data from the decomposition grid to the simulation grid
  CALL FMAP(f2, fg_edgH_yy, Dat_Grid_EdgH, Sim_Grid_EdgH, GMap_xyEdgH, VMap_xyEdgH, (dN_x+2)*(dN_y+1), N_x*(N_y+1), &
  (dN_x+2)*(dN_y+1)*dN_g, N_x*(N_y+1)*N_g, Open_Threads)

  !----------fg_edgH_xy----------!

  !reconstructing fg_edgH_xy from decomposition (on the same grid it was decomposed on)
  CALL DMD_RECONSTRUCT(f, L_fg_edgH_xy, B_fg_edgH_xy, W_fg_edgH_xy, C_fg_edgH_xy, rrank_fg_edgH_xy,&
  len, N_g, DMDgsum, Time, Open_Threads)

  !Appending on extra 'boundary data' left/right of the grid for mapping purposes
  CALL EXTND_LR_BND(f2, f, dN_g, dN_y+1, dN_x)

  !mapping the vector of fg_edgH_xy data from the decomposition grid to the simulation grid
  CALL FMAP(f2, fg_edgH_xy, Dat_Grid_EdgH, Sim_Grid_EdgH, GMap_xyEdgH, VMap_xyEdgH, (dN_x+2)*(dN_y+1), N_x*(N_y+1), &
  (dN_x+2)*(dN_y+1)*dN_g, N_x*(N_y+1)*N_g, Open_Threads)

  DEALLOCATE(f,f2) !deallocating DMD data vector

  !----------------------------------------
  !
  ! CELL AVERAGED DATA
  !
  !----------------------------------------

  len = ((dN_x+2)*(dN_y+2))*dN_g !length of cell-averaged data vector with boundary cells (for mapping purposes)
  ALLOCATE(f2(len)) !allocating space to hold reconstituted DMD data, with appended boundary cells

  len = dN_x*dN_y*dN_g !length of cell-averaged data vector
  ALLOCATE(f(len)) !allocating space to hold reconstituted DMD data

  !reconstructing fg_avg_xx from decomposition (on the same grid it was decomposed on)
  CALL DMD_RECONSTRUCT(f, L_fg_avg_xx, B_fg_avg_xx, W_fg_avg_xx, C_fg_avg_xx, rrank_fg_avg_xx,&
  len, N_g, DMDgsum, Time, Open_Threads)

  !appending boundary data to the vector of cell-averaged fg_xx (on the decomposition grid)
  !Note that this is only done for the purposes of mapping to the simulation grid
  ALLOCATE(bnd2((2*dN_x)*dN_g))
  CALL COLLECT_BT_BND(bnd2, f, dN_g, dN_y, dN_x)
  CALL COLLECT_CRN_BND(crnbnd, f, dN_g, dN_y, dN_x)

  CALL AVG_APND_BND(f2, f, xxbnd, bnd2, crnbnd, dN_g, dN_y, dN_x)
  DEALLOCATE(bnd2,xxbnd)

  !mapping the vector of fg_avg_xx data from the decomposition grid to the simulation grid
  CALL FMAP(f2, fg_avg_xx, Dat_Grid_Avg, Sim_Grid_Avg, GMap_xyAvg, VMap_xyAvg, (dN_x+2)*(dN_y+2), N_x*N_y, &
  (dN_x+2)*(dN_y+2)*dN_g, N_x*N_y*N_g, Open_Threads)

  !reconstructing fg_avg_yy from decomposition (on the same grid it was decomposed on)
  CALL DMD_RECONSTRUCT(f, L_fg_avg_yy, B_fg_avg_yy, W_fg_avg_yy, C_fg_avg_yy, rrank_fg_avg_yy,&
  len, N_g, DMDgsum, Time, Open_Threads)

  !appending boundary data to the vector of cell-averaged fg_xx (on the decomposition grid)
  !Note that this is only done for the purposes of mapping to the simulation grid
  ALLOCATE(bnd2((2*dN_y)*dN_g))
  CALL COLLECT_LR_BND(bnd2, f, dN_g, dN_y, dN_x)
  CALL COLLECT_CRN_BND(crnbnd, f, dN_g, dN_y, dN_x)

  CALL AVG_APND_BND(f2, f, bnd2, yybnd, crnbnd, dN_g, dN_y, dN_x)
  DEALLOCATE(bnd2, yybnd, crnbnd)

  !mapping the vector of fg_avg_yy data from the decomposition grid to the simulation grid
  CALL FMAP(f2, fg_avg_yy, Dat_Grid_Avg, Sim_Grid_Avg, GMap_xyAvg, VMap_xyAvg, (dN_x+2)*(dN_y+2), N_x*N_y, &
  (dN_x+2)*(dN_y+2)*dN_g, N_x*N_y*N_g, Open_Threads)

  DEALLOCATE(f,f2) !deallocating DMD data vector

END SUBROUTINE DMD_RECONSTRUCT_fg

!==================================================================================================================================!
! SUBROUTINE DMD_RECONSTRUCT_BCg
!
!DEF:
!  Reconstructs multigroup boundary factors at a specific time instant as a DMD expansion
!
!INPUT:
!  dN_g, dN_y, dN_x - size of the decomposotion grid in energy/ y-length/ x-length, respectively
!                   - Note that these are the grid sizes from the original decomposition, *NOT* the grid for the current simulation
!
!  N_g, N_y, N_x - size of the simulation grid in energy/ y-length/ x-length, respectively
!                - Note that these are the grid sizes defined for the current simulation
!
!  DMDgsum - flag to denote whether the given decomposition was performed over the entire phase space or groupwise for each boundary factor
!
!  Time - NOTe this is the Time RELATIVE TO THE START TIME OF THE DECOMPOSITION FIT
!
!OUTPUT:
!
!==================================================================================================================================!
SUBROUTINE DMD_RECONSTRUCT_BCg(Cg_L, Cg_B, Cg_R, Cg_T, L_BCg, B_BCg, W_BCg, C_BCg, rrank_BCg, dN_x, dN_y, dN_g, N_x, N_y, N_g,&
  Time, DMDgsum, Sim_Grid_Bnds, Dat_Grid_Bnds, GMap_xyBnds, VMap_xyBnds, Open_Threads)

  REAL*8,INTENT(OUT):: Cg_L(*), Cg_B(*), Cg_R(*), Cg_T(*)

  REAL*8,INTENT(IN):: C_BCg(*)
  COMPLEX*16,INTENT(IN):: L_BCg(*), B_BCg(*), W_BCg(*)
  REAL*8,INTENT(IN):: Sim_Grid_Bnds(*), Dat_Grid_Bnds(*)
  INTEGER,INTENT(IN):: rrank_BCg(*), GMap_xyBnds(*), VMap_xyBnds(*)
  INTEGER,INTENT(IN):: dN_x, dN_y, dN_g, N_x, N_y, N_g, DMDgsum, Open_Threads
  REAL*8,INTENT(IN):: Time

  REAL*8,ALLOCATABLE:: BC(:), BC2(:), CL2(:), CB2(:), CR2(:), CT2(:)
  INTEGER:: len, g, i, j, p, pl, pb, pr, pt

  len = dN_g*2*(dN_x+dN_y) !length of boundary-factor data vector
  ALLOCATE(BC(len)) !allocating space to hold reconstituted DMD data

  !Calculating the DMD expansion at t = Time
  CALL DMD_RECONSTRUCT(BC, L_BCg, B_BCg, W_BCg, C_BCg, rrank_BCg, len, N_g, DMDgsum, Time, Open_Threads)

  !temporary arrays to hold each boundary factor on the decomposition grid
  ALLOCATE(CL2(dN_g*(dN_y+2)), CB2(dN_g*(dN_x+2)), CR2(dN_g*(dN_y+2)), CT2(dN_g*(dN_x+2)))

  !partitioning the vector BC into the four seperate boundary factors
  p = 0
  pl = 0
  pb = 0
  pr = 0
  pt = 0
  DO g=1,dN_g

    !left boundary factor
    pl = pl + 1
    CL2(pl) = BC(p+1)
    DO i=1,dN_y
      p = p + 1
      pl = pl + 1
      CL2(pl) = BC(p)
    END DO
    pl = pl + 1
    CL2(pl) = BC(p)

    !bottom boundary factor
    pb = pb + 1
    CB2(pb) = BC(p+1)
    DO i=1,dN_x
      p = p + 1
      pb = pb + 1
      CB2(pb) = BC(p)
    END DO
    pb = pb + 1
    CB2(pb) = BC(p)

    !right boundary factor
    pr = pr + 1
    CR2(pr) = BC(p+1)
    DO i=1,dN_y
      p = p + 1
      pr = pr + 1
      CR2(pr) = BC(p)
    END DO
    pr = pr + 1
    CR2(pr) = BC(p)

    !top boundary factor
    pt = pt + 1
    CT2(pt) = BC(p+1)
    DO i=1,dN_x
      p = p + 1
      pt = pt + 1
      CT2(pt) = BC(p)
    END DO
    pt = pt + 1
    CT2(pt) = BC(p)

  END DO

  !mapping the vector of ... data from the decomposition grid to the simulation grid
  CALL FMAP(CL2, Cg_L, Dat_Grid_Bnds(1), Sim_Grid_Bnds(1), GMap_xyBnds(1), VMap_xyBnds(1),&
  dN_y+2, N_y, dN_g*(dN_y+2), N_g*N_y, Open_Threads)

  !mapping the vector of ... data from the decomposition grid to the simulation grid
  CALL FMAP(CB2, Cg_B, Dat_Grid_Bnds(2*dN_y+5), Sim_Grid_Bnds(2*N_y+1), GMap_xyBnds(4*N_y+1), VMap_xyBnds(4*N_y+1), &
  dN_x+2, N_x, dN_g*(dN_x+2), N_g*N_x, Open_Threads)

  !mapping the vector of ... data from the decomposition grid to the simulation grid
  CALL FMAP(CR2, Cg_R, Dat_Grid_Bnds(2*(dN_y+dN_x)+9), Sim_Grid_Bnds(2*(N_y+N_x)+1), GMap_xyBnds(4*(N_y+N_x)+1), &
  VMap_xyBnds(4*(N_y+N_x)+1), dN_y+2, N_y, dN_g*(dN_y+2), N_g*N_y, Open_Threads)

  !mapping the vector of ... data from the decomposition grid to the simulation grid
  CALL FMAP(CT2, Cg_T, Dat_Grid_Bnds(4*dN_y+2*dN_x+13), Sim_Grid_Bnds(4*N_y+2*N_x+1), GMap_xyBnds(8*N_y+4*N_x+1), &
  VMap_xyBnds(8*N_y+4*N_x+1), dN_x+2, N_x, dN_g*(dN_x+2), N_g*N_x, Open_Threads)

  DEALLOCATE(BC,CL2,CB2,CR2,CT2)


END SUBROUTINE DMD_RECONSTRUCT_BCg


!==================================================================================================================================!
!
!==================================================================================================================================!
SUBROUTINE INPUT_fg_DMD(Fname, DMDgsum, N_x, N_y, N_g, L_BCg, W_BCg, B_BCg, C_BCg, rrank_BCg, L_fg_avg_xx, W_fg_avg_xx,&
  B_fg_avg_xx, C_fg_avg_xx, rrank_fg_avg_xx, L_fg_edgV_xx, W_fg_edgV_xx, B_fg_edgV_xx, C_fg_edgV_xx, rrank_fg_edgV_xx,&
  L_fg_avg_yy, W_fg_avg_yy, B_fg_avg_yy, C_fg_avg_yy, rrank_fg_avg_yy, L_fg_edgH_yy, W_fg_edgH_yy, B_fg_edgH_yy, C_fg_edgH_yy,&
  rrank_fg_edgH_yy, L_fg_edgV_xy, W_fg_edgV_xy, B_fg_edgV_xy, C_fg_edgV_xy, rrank_fg_edgV_xy, L_fg_edgH_xy, W_fg_edgH_xy,&
  B_fg_edgH_xy, C_fg_edgH_xy, rrank_fg_edgH_xy, xlen, ylen, Tini, Delx, Dely, bcT, BC_Type, Start_Time)

  CHARACTER(*),INTENT(IN):: Fname
  INTEGER,INTENT(OUT):: DMDgsum
  REAL*8,INTENT(OUT):: xlen, ylen, Tini, Start_Time
  REAL*8,ALLOCATABLE,INTENT(OUT):: Delx(:), Dely(:), bcT(:)
  REAL*8,ALLOCATABLE,INTENT(OUT):: C_BCg(:), C_fg_avg_xx(:), C_fg_avg_yy(:)
  REAL*8,ALLOCATABLE,INTENT(OUT):: C_fg_edgV_xx(:), C_fg_edgH_yy(:), C_fg_edgV_xy(:), C_fg_edgH_xy(:)
  COMPLEX*16,ALLOCATABLE,INTENT(OUT):: L_BCg(:),        W_BCg(:),        B_BCg(:)
  COMPLEX*16,ALLOCATABLE,INTENT(OUT):: L_fg_avg_xx(:),  W_fg_avg_xx(:),  B_fg_avg_xx(:)
  COMPLEX*16,ALLOCATABLE,INTENT(OUT):: L_fg_edgV_xx(:), W_fg_edgV_xx(:), B_fg_edgV_xx(:)
  COMPLEX*16,ALLOCATABLE,INTENT(OUT):: L_fg_avg_yy(:),  W_fg_avg_yy(:),  B_fg_avg_yy(:)
  COMPLEX*16,ALLOCATABLE,INTENT(OUT):: L_fg_edgH_yy(:), W_fg_edgH_yy(:), B_fg_edgH_yy(:)
  COMPLEX*16,ALLOCATABLE,INTENT(OUT):: L_fg_edgV_xy(:), W_fg_edgV_xy(:), B_fg_edgV_xy(:)
  COMPLEX*16,ALLOCATABLE,INTENT(OUT):: L_fg_edgH_xy(:), W_fg_edgH_xy(:), B_fg_edgH_xy(:)
  INTEGER,INTENT(OUT):: N_x, N_y, N_g
  INTEGER,ALLOCATABLE,INTENT(OUT):: rrank_fg_avg_xx(:), rrank_fg_edgV_xx(:), rrank_fg_avg_yy(:), rrank_fg_edgH_yy(:)
  INTEGER,ALLOCATABLE,INTENT(OUT):: rrank_fg_edgV_xy(:), rrank_fg_edgH_xy(:), rrank_BCg(:)
  INTEGER,ALLOCATABLE,INTENT(OUT):: BC_Type(:)

  CHARACTER(10):: dcmp_type
  INTEGER:: ncID, clen, g, Status, ID2
  REAL*8:: bnds(2)
  CHARACTER(50):: Location = 'MODULE: DMD_ROUTINES / SUBROUTINE: INPUT_fg_DMD'

  CALL NF_OPEN_FILE(ncID,Fname,'old','r')

  CALL NF_INQ_DIM(ncID,"N_x",N_x)
  CALL NF_INQ_DIM(ncID,"N_y",N_y)
  CALL NF_INQ_DIM(ncID,"N_g",N_g)

  ALLOCATE(Delx(N_x),Dely(N_y),BC_Type(4),bcT(4))

  CALL NF_INQ_VAR_0D(ncID,"xlen",xlen)
  CALL NF_INQ_VAR_0D(ncID,"ylen",ylen)
  CALL NF_INQ_VAR_1D(ncID,"Delx",Delx,(/1/),(/N_x/))
  CALL NF_INQ_VAR_1D(ncID,"Dely",Dely,(/1/),(/N_y/))
  CALL NF_INQ_VAR_0D(ncID,"bcT_left",bcT(1))
  CALL NF_INQ_VAR_0D(ncID,"bcT_bottom",bcT(2))
  CALL NF_INQ_VAR_0D(ncID,"bcT_right",bcT(3))
  CALL NF_INQ_VAR_0D(ncID,"bcT_top",bcT(4))
  CALL NF_INQ_VAR_0D(ncID,"Tini",Tini)

  Status = nf90_inq_varid(ncID, "tpts", ID2)
  CALL HANDLE_ERR(Status, Location)
  Status = nf90_get_att(ncID, ID2, "bnds", bnds)
  CALL HANDLE_ERR(Status, Location)
  Start_Time = bnds(1)

  Status = nf90_get_att(ncID, NF90_GLOBAL, "dcmp_type", dcmp_type)
  CALL HANDLE_ERR(Status,Location)
  IF (dcmp_type(1:4) .EQ. "DMDg") THEN
    DMDgsum = 0
  ELSE IF (dcmp_type(1:3) .EQ. "DMD") THEN
    DMDgsum = 1
  ELSE
    WRITE(*,'(5A)') "DMD_ROUTINES :: INPUT_fg_DMD - Unrecognized dcmp_type (",trim(dcmp_type),") in dataset (",trim(Fname),")"
    STOP
  END IF

  clen = 2*(N_x+N_y)
  CALL READ_DMD(ncID, 'BCg', DMDgsum, clen, N_g, L_BCg, W_BCg, B_BCg, C_BCg, rrank_BCg)

  clen = N_x*N_y
  CALL READ_DMD(ncID, 'fg_avg_xx', DMDgsum, clen, N_g, L_fg_avg_xx, W_fg_avg_xx, B_fg_avg_xx, C_fg_avg_xx, rrank_fg_avg_xx)

  clen = (N_x+1)*N_y
  CALL READ_DMD(ncID, 'fg_edgV_xx', DMDgsum, clen, N_g, L_fg_edgV_xx, W_fg_edgV_xx, B_fg_edgV_xx, C_fg_edgV_xx, rrank_fg_edgV_xx)

  clen = N_x*N_y
  CALL READ_DMD(ncID, 'fg_avg_yy', DMDgsum, clen, N_g, L_fg_avg_yy, W_fg_avg_yy, B_fg_avg_yy, C_fg_avg_yy, rrank_fg_avg_yy)

  clen = N_x*(N_y+1)
  CALL READ_DMD(ncID, 'fg_edgH_yy', DMDgsum, clen, N_g, L_fg_edgH_yy, W_fg_edgH_yy, B_fg_edgH_yy, C_fg_edgH_yy, rrank_fg_edgH_yy)

  clen = (N_x+1)*N_y
  CALL READ_DMD(ncID, 'fg_edgV_xy', DMDgsum, clen, N_g, L_fg_edgV_xy, W_fg_edgV_xy, B_fg_edgV_xy, C_fg_edgV_xy, rrank_fg_edgV_xy)

  clen = N_x*(N_y+1)
  CALL READ_DMD(ncID, 'fg_edgH_xy', DMDgsum, clen, N_g, L_fg_edgH_xy, W_fg_edgH_xy, B_fg_edgH_xy, C_fg_edgH_xy, rrank_fg_edgH_xy)

  CALL NF_CLOSE_FILE(ncID)

END SUBROUTINE INPUT_fg_DMD

!==================================================================================================================================!
!Subroutine READ_DMD
!
! NOTES::
!
! WARNINGS::
!
! OUTPUTS::
!
! INPUTS::
!
!==================================================================================================================================!
SUBROUTINE READ_DMD(ncID, name, DMDgsum, clen, N_g, L, W, B, Base, N_modes)
  INTEGER,INTENT(IN):: ncID, clen, DMDgsum, N_g
  CHARACTER(*),INTENT(IN):: name
  INTEGER,ALLOCATABLE,INTENT(OUT):: N_modes(:)
  COMPLEX*16,ALLOCATABLE,INTENT(OUT):: L(:), W(:), B(:)
  REAL*8,ALLOCATABLE,INTENT(OUT):: Base(:)
  REAL*8,ALLOCATABLE:: r(:), c(:)
  CHARACTER(100):: buf1, buf2
  CHARACTER(50):: Location = 'MODULE: DMD_ROUTINES / SUBROUTINE: READ_DMD'
  INTEGER:: ModeCount, ModeSum, g, i, j
  INTEGER:: Status

  !Checking if using groupwise DMD or DMD over full phase space
  IF (DMDgsum .EQ. 1) THEN
    ALLOCATE(N_modes(1))
  ELSE
    ALLOCATE(N_modes(N_g))
  END IF

  !Reading in the number of DMD modes (N_modes)
  write(buf1, '(2A)') "N_modes_",name
  Status = nf90_get_att(ncID, NF90_GLOBAL, buf1, N_modes)
  CALL HANDLE_ERR(Status,Location)

  ModeCount = SUM(N_modes) !collecting total number of DMD modes

  !--------------------------------------------------!
  !                 DMD Eigenvalues                  !
  !--------------------------------------------------!
  IF (ALLOCATED(L)) DEALLOCATE(L)
  ALLOCATE(L(ModeCount), r(ModeCount), c(ModeCount))
  WRITE(buf1,'(2A)') 'eL_real_',name
  WRITE(buf2,'(2A)') 'eL_imag_',name
  IF (DMDgsum .EQ. 1) THEN
    CALL NF_INQ_VAR_1D(ncID, buf1, r, (/1/), (/N_modes(1)/))
    CALL NF_INQ_VAR_1D(ncID, buf2, c, (/1/), (/N_modes(1)/))
  ELSE
    ModeSum = 1
    DO g=1,N_g
      CALL NF_INQ_VAR_1D(ncID, buf1, r( ModeSum : ModeSum + N_modes(g) - 1 ), (/1,g/), (/N_modes(g),1/))
      CALL NF_INQ_VAR_1D(ncID, buf2, c( ModeSum : ModeSum + N_modes(g) - 1 ), (/1,g/), (/N_modes(g),1/))
      ModeSum = ModeSum + N_modes(g)
    END DO
  END IF

  DO i=1,ModeCount
    L(i) = COMPLEX(r(i), c(i))
  END DO

  !--------------------------------------------------!
  !           DMD Expansion Coefficients             !
  !--------------------------------------------------!
  IF (ALLOCATED(B)) DEALLOCATE(B)
  ALLOCATE(B(ModeCount))
  WRITE(buf1,'(2A)') 'b_real_',name
  WRITE(buf2,'(2A)') 'b_imag_',name
  IF (DMDgsum .EQ. 1) THEN
    CALL NF_INQ_VAR_1D(ncID, buf1, r, (/1/), (/N_modes(1)/))
    CALL NF_INQ_VAR_1D(ncID, buf2, c, (/1/), (/N_modes(1)/))
  ELSE
    ModeSum = 1
    DO g=1,N_g
      CALL NF_INQ_VAR_1D(ncID, buf1, r( ModeSum : ModeSum + N_modes(g) - 1 ), (/1,g/), (/N_modes(g),1/))
      CALL NF_INQ_VAR_1D(ncID, buf2, c( ModeSum : ModeSum + N_modes(g) - 1 ), (/1,g/), (/N_modes(g),1/))
      ModeSum = ModeSum + N_modes(g)
    END DO
  END IF

  DO i=1,ModeCount
    B(i) = COMPLEX(r(i), c(i))
  END DO
  DEALLOCATE(r, c)

  !--------------------------------------------------!
  !                 DMD 'Base Flow'                  !
  !--------------------------------------------------!
  IF (ALLOCATED(Base)) DEALLOCATE(Base)
  ALLOCATE(Base(clen*N_g))
  WRITE(buf1,'(2A)') 'C_',name
  IF (DMDgsum .EQ. 1) THEN
    CALL NF_INQ_VAR_1D(ncID, buf1, Base, (/1/), (/clen*N_g/))
  ELSE
    ModeSum = 1
    DO g=1,N_g
      CALL NF_INQ_VAR_1D(ncID, buf1, Base( ModeSum : ModeSum + N_modes(g) - 1 ), (/1,g/), (/clen,1/))
      ModeSum = ModeSum + clen
    END DO
  END IF

  !--------------------------------------------------!
  !                    DMD Modes                     !
  !--------------------------------------------------!
  IF (ALLOCATED(W)) DEALLOCATE(W)
  ModeCount = ModeCount * clen * N_g
  ALLOCATE(W(ModeCount), r(ModeCount), c(ModeCount))
  WRITE(buf1,'(2A)') 'Wproj_real_',name
  WRITE(buf2,'(2A)') 'Wproj_imag_',name
  IF (DMDgsum .EQ. 1) THEN
    CALL NF_INQ_VAR_1D(ncID, buf1, r, (/1,1/), (/clen*N_g, N_modes(1)/))
    CALL NF_INQ_VAR_1D(ncID, buf2, c, (/1,1/), (/clen*N_g, N_modes(1)/))
  ELSE
    ModeSum = 1
    DO g=1,N_g
      CALL NF_INQ_VAR_1D(ncID, buf1, r( ModeSum : ModeSum + N_modes(g)*clen - 1 ), (/1,1,g/), (/clen, N_modes(g), 1/))
      CALL NF_INQ_VAR_1D(ncID, buf2, c( ModeSum : ModeSum + N_modes(g)*clen - 1 ), (/1,1,g/), (/clen, N_modes(g), 1/))
      ModeSum = ModeSum + N_modes(g)*clen
    END DO
  END IF

  ! ModeCount = ModeCount * clen * N_g
  DO i=1,ModeCount
    W(i) = COMPLEX(r(i), c(i))
  END DO
  DEALLOCATE(r, c)

END SUBROUTINE READ_DMD

!==================================================================================================================================!
!
!==================================================================================================================================!

END MODULE DMD_ROUTINES
