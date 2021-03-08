MODULE POD_ROUTINES

  USE NCDF_IO
  USE GRID_FUNCTIONS

  IMPLICIT NONE

INTERFACE READ_POD
  MODULE PROCEDURE READ_POD_FULL, READ_POD_GPART
END INTERFACE

CONTAINS

!==================================================================================================================================!
! SUBROUTINE POD_RECONSTRUCT_fg
!
!DEF:
!  Reconstructs multigroup qd-factors at a specific time instant (fg(t)) from a truncated SVD representation
!
!INPUT:
!  N_t, N_g, N_y, N_x - size of the grid in time/ energy/ y-length/ x-length, respectively
!                     - Note that these are the grid sizes from the original decomposition, *NOT* the grid for the current simulation
!  PODgsum - flag to denote whether the given decomposition was performed over the entire phase space or groupwise for each qd factor
!
!OUTPUT:
!
!==================================================================================================================================!
SUBROUTINE POD_RECONSTRUCT_fg(fg_avg_xx,fg_avg_yy,fg_edgV_xx,fg_edgV_xy,fg_edgH_yy,fg_edgH_xy,C_fg_avg_xx,&
  S_fg_avg_xx,U_fg_avg_xx,V_fg_avg_xx,C_fg_edgV_xx,S_fg_edgV_xx,U_fg_edgV_xx,V_fg_edgV_xx,C_fg_avg_yy,S_fg_avg_yy,&
  U_fg_avg_yy,V_fg_avg_yy,C_fg_edgH_yy,S_fg_edgH_yy,U_fg_edgH_yy,V_fg_edgH_yy,C_fg_edgV_xy,S_fg_edgV_xy,U_fg_edgV_xy,&
  V_fg_edgV_xy,C_fg_edgH_xy,S_fg_edgH_xy,U_fg_edgH_xy,V_fg_edgH_xy,rrank_fg_avg_xx,rrank_fg_edgV_xx,rrank_fg_avg_yy,&
  rrank_fg_edgH_yy,rrank_fg_edgV_xy,rrank_fg_edgH_xy,dN_x,dN_y,dN_g,dN_t,N_x,N_y,N_g,N_t,t,PODgsum,Sim_Grid_Avg,&
  Sim_Grid_EdgV,Sim_Grid_EdgH,Dat_Grid_Avg,Dat_Grid_EdgV,Dat_Grid_EdgH,Sim_TGrid,Dat_TGrid,GMap_xyAvg,GMap_xyEdgV,&
  GMap_xyEdgH,VMap_xyAvg,VMap_xyEdgV,VMap_xyEdgH,TMap,Open_Threads)

  REAL*8,INTENT(OUT):: fg_avg_xx(*), fg_avg_yy(*)
  REAL*8,INTENT(OUT):: fg_edgV_xx(*), fg_edgV_xy(*)
  REAL*8,INTENT(OUT):: fg_edgH_yy(*), fg_edgH_xy(*)

  REAL*8,INTENT(IN):: C_fg_avg_xx(*), S_fg_avg_xx(*), U_fg_avg_xx(*), V_fg_avg_xx(*)
  REAL*8,INTENT(IN):: C_fg_edgV_xx(*), S_fg_edgV_xx(*), U_fg_edgV_xx(*), V_fg_edgV_xx(*)
  REAL*8,INTENT(IN):: C_fg_avg_yy(*), S_fg_avg_yy(*), U_fg_avg_yy(*), V_fg_avg_yy(*)
  REAL*8,INTENT(IN):: C_fg_edgH_yy(*), S_fg_edgH_yy(*), U_fg_edgH_yy(*), V_fg_edgH_yy(*)
  REAL*8,INTENT(IN):: C_fg_edgV_xy(*), S_fg_edgV_xy(*), U_fg_edgV_xy(*), V_fg_edgV_xy(*)
  REAL*8,INTENT(IN):: C_fg_edgH_xy(*), S_fg_edgH_xy(*), U_fg_edgH_xy(*), V_fg_edgH_xy(*)
  INTEGER,INTENT(IN):: rrank_fg_avg_xx(*), rrank_fg_edgV_xx(*), rrank_fg_avg_yy(*), rrank_fg_edgH_yy(*)
  INTEGER,INTENT(IN):: rrank_fg_edgV_xy(*), rrank_fg_edgH_xy(*)
  INTEGER,INTENT(IN):: dN_x, dN_y, dN_g, dN_t, N_x, N_y, N_g, N_t, t, PODgsum
  REAL*8,INTENT(IN):: Sim_Grid_Avg(*), Sim_Grid_EdgV(*), Sim_Grid_EdgH(*)
  REAL*8,INTENT(IN):: Dat_Grid_Avg(*), Dat_Grid_EdgV(*), Dat_Grid_EdgH(*)
  REAL*8,INTENT(IN):: Sim_TGrid(*), Dat_TGrid(*)
  INTEGER,INTENT(IN):: GMap_xyAvg(*), GMap_xyEdgV(*), GMap_xyEdgH(*)
  INTEGER,INTENT(IN):: VMap_xyAvg(*), VMap_xyEdgV(*), VMap_xyEdgH(*)
  INTEGER,INTENT(IN):: TMap(*)
  INTEGER,INTENT(IN):: Open_Threads

  REAL*8,ALLOCATABLE:: f(:), f2(:), f3(:), xxbnd(:), yybnd(:), crnbnd(:), bnd2(:)
  INTEGER:: len, g, i, j, p1, p2, p3

  ALLOCATE(xxbnd((2*dN_y)*dN_g), yybnd((2*dN_x)*dN_g), crnbnd(4*dN_g))

  !----------------------------------------
  !
  ! VERTICAL (y=const) CELL FACE DATA
  !
  !----------------------------------------

  !----------vector lengths----------!

  len = (dN_x+1)*(dN_y+2)*dN_g !length of EdgV data vector with boundary cells above/below the grid
  ALLOCATE(f2(len)) !allocating space to hold reconstituted POD data, with appended boundary cells

  len = (dN_x+1)*dN_y*dN_g !length of EdgV-type data vector
  ALLOCATE(f(len),f3(len)) !allocating space to hold reconstituted POD data

  !----------fg_edgV_xx----------!

  !reconstructing fg_edgV_xx from decomposition (on the same grid it was decomposed on)
  ! CALL POD_RECONSTRUCT(f,C_fg_edgV_xx,S_fg_edgV_xx,U_fg_edgV_xx,V_fg_edgV_xx,rrank_fg_edgV_xx,len,dN_g,dN_t,t,PODgsum)
  IF (TMap(2*t-1) .EQ. TMap(2*t)) THEN
    CALL POD_RECONSTRUCT(f,C_fg_edgV_xx,S_fg_edgV_xx,U_fg_edgV_xx,V_fg_edgV_xx,rrank_fg_edgV_xx,len,dN_g,dN_t,TMap(2*t),PODgsum)
  ELSE
    CALL POD_RECONSTRUCT(f,C_fg_edgV_xx,S_fg_edgV_xx,U_fg_edgV_xx,V_fg_edgV_xx,rrank_fg_edgV_xx,len,dN_g,dN_t,TMap(2*t-1),PODgsum)
    CALL POD_RECONSTRUCT(f3,C_fg_edgV_xx,S_fg_edgV_xx,U_fg_edgV_xx,V_fg_edgV_xx,rrank_fg_edgV_xx,len,dN_g,dN_t,TMap(2*t),PODgsum)
    CALL LINEAR_INTERPOLATION(f,f3,Sim_TGrid(t),Dat_TGrid(TMap(2*t-1)),Dat_TGrid(TMap(2*t)),len)
  END IF

  !storing the fg_xx boundary data (for the left/right boundaries only)
  CALL COLLECT_LR_BND(xxbnd,f,dN_g,dN_y,dN_x+1)

  !Appending on extra 'boundary data' above/below the grid for mapping purposes
  CALL EXTND_BT_BND(f2,f,dN_g,dN_y,dN_x+1)

  !mapping the vector of fg_edgV_xx data from the decomposition grid to the simulation grid
  CALL FMAP(f2,fg_edgV_xx,Dat_Grid_EdgV,Sim_Grid_EdgV,GMap_xyEdgV,VMap_xyEdgV, (dN_x+1)*(dN_y+2), (N_x+1)*N_y, &
  (dN_x+1)*(dN_y+2)*N_g, (N_x+1)*N_y*N_g, Open_Threads)

  !----------fg_edgV_xy----------!

  !reconstructing fg_edgV_xy from decomposition (on the same grid it was decomposed on)
  ! CALL POD_RECONSTRUCT(f,C_fg_edgV_xy,S_fg_edgV_xy,U_fg_edgV_xy,V_fg_edgV_xy,rrank_fg_edgV_xy,len,dN_g,dN_t,t,PODgsum)
  IF (TMap(2*t-1) .EQ. TMap(2*t)) THEN
    CALL POD_RECONSTRUCT(f,C_fg_edgV_xy,S_fg_edgV_xy,U_fg_edgV_xy,V_fg_edgV_xy,rrank_fg_edgV_xy,len,dN_g,dN_t,TMap(2*t),PODgsum)
  ELSE
    CALL POD_RECONSTRUCT(f,C_fg_edgV_xy,S_fg_edgV_xy,U_fg_edgV_xy,V_fg_edgV_xy,rrank_fg_edgV_xy,len,dN_g,dN_t,TMap(2*t-1),PODgsum)
    CALL POD_RECONSTRUCT(f3,C_fg_edgV_xy,S_fg_edgV_xy,U_fg_edgV_xy,V_fg_edgV_xy,rrank_fg_edgV_xy,len,dN_g,dN_t,TMap(2*t),PODgsum)
    CALL LINEAR_INTERPOLATION(f,f3,Sim_TGrid(t),Dat_TGrid(TMap(2*t-1)),Dat_TGrid(TMap(2*t)),len)
  END IF

  !Appending on extra 'boundary data' above/below the grid for mapping purposes
  CALL EXTND_BT_BND(f2,f,dN_g,dN_y,dN_x+1)

  !mapping the vector of fg_edgV_xy data from the decomposition grid to the simulation grid
  CALL FMAP(f2,fg_edgV_xy,Dat_Grid_EdgV,Sim_Grid_EdgV,GMap_xyEdgV,VMap_xyEdgV, (dN_x+1)*(dN_y+2), (N_x+1)*N_y, &
  (dN_x+1)*(dN_y+2)*N_g, (N_x+1)*N_y*N_g, Open_Threads)

  DEALLOCATE(f,f2,f3) !deallocating POD data vector


  !----------------------------------------
  !
  ! HORIZONTAL (x=const) CELL FACE DATA
  !
  !----------------------------------------

  !----------vector lengths----------!

  len = (dN_x+2)*(dN_y+1)*dN_g !length of EdgH data vector with boundary cells left/right of the grid
  ALLOCATE(f2(len)) !allocating space to hold reconstituted POD data, with appended boundary cells

  len = dN_x*(dN_y+1)*dN_g !length of EdgH-type data vector
  ALLOCATE(f(len),f3(len)) !allocating space to hold reconstituted POD data

  !----------fg_edgH_yy----------!

  !reconstructing fg_edgH_yy from decomposition (on the same grid it was decomposed on)
  ! CALL POD_RECONSTRUCT(f,C_fg_edgH_yy,S_fg_edgH_yy,U_fg_edgH_yy,V_fg_edgH_yy,rrank_fg_edgH_yy,len,dN_g,dN_t,t,PODgsum)
  IF (TMap(2*t-1) .EQ. TMap(2*t)) THEN
    CALL POD_RECONSTRUCT(f,C_fg_edgH_yy,S_fg_edgH_yy,U_fg_edgH_yy,V_fg_edgH_yy,rrank_fg_edgH_yy,len,dN_g,dN_t,TMap(2*t),PODgsum)
  ELSE
    CALL POD_RECONSTRUCT(f,C_fg_edgH_yy,S_fg_edgH_yy,U_fg_edgH_yy,V_fg_edgH_yy,rrank_fg_edgH_yy,len,dN_g,dN_t,TMap(2*t-1),PODgsum)
    CALL POD_RECONSTRUCT(f3,C_fg_edgH_yy,S_fg_edgH_yy,U_fg_edgH_yy,V_fg_edgH_yy,rrank_fg_edgH_yy,len,dN_g,dN_t,TMap(2*t),PODgsum)
    CALL LINEAR_INTERPOLATION(f,f3,Sim_TGrid(t),Dat_TGrid(TMap(2*t-1)),Dat_TGrid(TMap(2*t)),len)
  END IF

  !storing the fg_yy boundary data (for the top/bottom boundaries only)
  CALL COLLECT_BT_BND(yybnd,f,dN_g,dN_y+1,dN_x)

  !Appending on extra 'boundary data' left/right of the grid for mapping purposes
  CALL EXTND_LR_BND(f2,f,dN_g,dN_y+1,dN_x)

  !mapping the vector of fg_edgH_yy data from the decomposition grid to the simulation grid
  CALL FMAP(f2,fg_edgH_yy,Dat_Grid_EdgH,Sim_Grid_EdgH,GMap_xyEdgH,VMap_xyEdgH, (dN_x+2)*(dN_y+1), N_x*(N_y+1), &
  (dN_x+2)*(dN_y+1)*dN_g, N_x*(N_y+1)*N_g, Open_Threads)

  !----------fg_edgH_xy----------!

  !reconstructing fg_edgH_xy from decomposition (on the same grid it was decomposed on)
  ! CALL POD_RECONSTRUCT(f,C_fg_edgH_xy,S_fg_edgH_xy,U_fg_edgH_xy,V_fg_edgH_xy,rrank_fg_edgH_xy,len,dN_g,dN_t,t,PODgsum)
  IF (TMap(2*t-1) .EQ. TMap(2*t)) THEN
    CALL POD_RECONSTRUCT(f,C_fg_edgH_xy,S_fg_edgH_xy,U_fg_edgH_xy,V_fg_edgH_xy,rrank_fg_edgH_xy,len,dN_g,dN_t,TMap(2*t),PODgsum)
  ELSE
    CALL POD_RECONSTRUCT(f,C_fg_edgH_xy,S_fg_edgH_xy,U_fg_edgH_xy,V_fg_edgH_xy,rrank_fg_edgH_xy,len,dN_g,dN_t,TMap(2*t-1),PODgsum)
    CALL POD_RECONSTRUCT(f3,C_fg_edgH_xy,S_fg_edgH_xy,U_fg_edgH_xy,V_fg_edgH_xy,rrank_fg_edgH_xy,len,dN_g,dN_t,TMap(2*t),PODgsum)
    CALL LINEAR_INTERPOLATION(f,f3,Sim_TGrid(t),Dat_TGrid(TMap(2*t-1)),Dat_TGrid(TMap(2*t)),len)
  END IF

  !Appending on extra 'boundary data' left/right of the grid for mapping purposes
  CALL EXTND_LR_BND(f2,f,dN_g,dN_y+1,dN_x)

  !mapping the vector of fg_edgH_xy data from the decomposition grid to the simulation grid
  CALL FMAP(f2,fg_edgH_xy,Dat_Grid_EdgH,Sim_Grid_EdgH,GMap_xyEdgH,VMap_xyEdgH, (dN_x+2)*(dN_y+1), N_x*(N_y+1), &
  (dN_x+2)*(dN_y+1)*dN_g, N_x*(N_y+1)*N_g, Open_Threads)

  DEALLOCATE(f,f2,f3) !deallocating POD data vector

  !----------------------------------------
  !
  ! CELL AVERAGED DATA
  !
  !----------------------------------------

  len = ((dN_x+2)*(dN_y+2))*dN_g !length of cell-averaged data vector with boundary cells
  ALLOCATE(f2(len)) !allocating space to hold reconstituted POD data, with appended boundary cells

  len = dN_x*dN_y*dN_g !length of cell-averaged data vector
  ALLOCATE(f(len),f3(len)) !allocating space to hold reconstituted POD data

  !reconstructing fg_avg_xx from decomposition (on the same grid it was decomposed on)
  ! CALL POD_RECONSTRUCT(f,C_fg_avg_xx,S_fg_avg_xx,U_fg_avg_xx,V_fg_avg_xx,rrank_fg_avg_xx,len,dN_g,dN_t,t,PODgsum)
  IF (TMap(2*t-1) .EQ. TMap(2*t)) THEN
    CALL POD_RECONSTRUCT(f,C_fg_avg_xx,S_fg_avg_xx,U_fg_avg_xx,V_fg_avg_xx,rrank_fg_avg_xx,len,dN_g,dN_t,TMap(2*t),PODgsum)
  ELSE
    CALL POD_RECONSTRUCT(f,C_fg_avg_xx,S_fg_avg_xx,U_fg_avg_xx,V_fg_avg_xx,rrank_fg_avg_xx,len,dN_g,dN_t,TMap(2*t-1),PODgsum)
    CALL POD_RECONSTRUCT(f3,C_fg_avg_xx,S_fg_avg_xx,U_fg_avg_xx,V_fg_avg_xx,rrank_fg_avg_xx,len,dN_g,dN_t,TMap(2*t),PODgsum)
    CALL LINEAR_INTERPOLATION(f,f3,Sim_TGrid(t),Dat_TGrid(TMap(2*t-1)),Dat_TGrid(TMap(2*t)),len)
  END IF

  !appending boundary data to the vector of cell-averaged fg_xx (on the decomposition grid)
  !Note that this is only done for the purposes of mapping to the simulation grid
  ALLOCATE(bnd2((2*dN_x)*dN_g))
  CALL COLLECT_BT_BND(bnd2,f,dN_g,dN_y,dN_x)
  CALL COLLECT_CRN_BND(crnbnd,f,dN_g,dN_y,dN_x)

  CALL AVG_APND_BND(f2,f,xxbnd,bnd2,crnbnd,dN_g,dN_y,dN_x)
  DEALLOCATE(bnd2,xxbnd)

  !mapping the vector of fg_avg_xx data from the decomposition grid to the simulation grid
  CALL FMAP(f2,fg_avg_xx,Dat_Grid_Avg,Sim_Grid_Avg,GMap_xyAvg,VMap_xyAvg, (dN_x+2)*(dN_y+2), N_x*N_y, &
  (dN_x+2)*(dN_y+2)*dN_g, N_x*N_y*N_g, Open_Threads)

  !reconstructing fg_avg_yy from decomposition (on the same grid it was decomposed on)
  ! CALL POD_RECONSTRUCT(f,C_fg_avg_yy,S_fg_avg_yy,U_fg_avg_yy,V_fg_avg_yy,rrank_fg_avg_yy,len,dN_g,dN_t,t,PODgsum)
  IF (TMap(2*t-1) .EQ. TMap(2*t)) THEN
    CALL POD_RECONSTRUCT(f,C_fg_avg_yy,S_fg_avg_yy,U_fg_avg_yy,V_fg_avg_yy,rrank_fg_avg_yy,len,dN_g,dN_t,TMap(2*t),PODgsum)
  ELSE
    CALL POD_RECONSTRUCT(f,C_fg_avg_yy,S_fg_avg_yy,U_fg_avg_yy,V_fg_avg_yy,rrank_fg_avg_yy,len,dN_g,dN_t,TMap(2*t-1),PODgsum)
    CALL POD_RECONSTRUCT(f3,C_fg_avg_yy,S_fg_avg_yy,U_fg_avg_yy,V_fg_avg_yy,rrank_fg_avg_yy,len,dN_g,dN_t,TMap(2*t),PODgsum)
    CALL LINEAR_INTERPOLATION(f,f3,Sim_TGrid(t),Dat_TGrid(TMap(2*t-1)),Dat_TGrid(TMap(2*t)),len)
  END IF

  !appending boundary data to the vector of cell-averaged fg_xx (on the decomposition grid)
  !Note that this is only done for the purposes of mapping to the simulation grid
  ALLOCATE(bnd2((2*dN_y)*dN_g))
  CALL COLLECT_LR_BND(bnd2,f,dN_g,dN_y,dN_x)
  CALL COLLECT_CRN_BND(crnbnd,f,dN_g,dN_y,dN_x)

  CALL AVG_APND_BND(f2,f,bnd2,yybnd,crnbnd,dN_g,dN_y,dN_x)
  DEALLOCATE(bnd2,yybnd,crnbnd)

  !mapping the vector of fg_avg_yy data from the decomposition grid to the simulation grid
  CALL FMAP(f2,fg_avg_yy,Dat_Grid_Avg,Sim_Grid_Avg,GMap_xyAvg,VMap_xyAvg, (dN_x+2)*(dN_y+2), N_x*N_y, &
  (dN_x+2)*(dN_y+2)*dN_g, N_x*N_y*N_g, Open_Threads)

  DEALLOCATE(f,f2,f3) !deallocating POD data vector

END SUBROUTINE POD_RECONSTRUCT_fg

!==================================================================================================================================!
! SUBROUTINE POD_RECONSTRUCT_fg
!
!DEF:
!  Reconstructs multigroup qd-factors at a specific time instant (fg(t)) from a truncated SVD representation
!
!INPUT:
!  N_t, N_g, N_y, N_x - size of the grid in time/ energy/ y-length/ x-length, respectively
!                     - Note that these are the grid sizes from the original decomposition, *NOT* the grid for the current simulation
!  PODgsum - flag to denote whether the given decomposition was performed over the entire phase space or groupwise for each qd factor
!
!OUTPUT:
!
!==================================================================================================================================!
SUBROUTINE POD_RECONSTRUCT_BCg(Cg_L,Cg_B,Cg_R,Cg_T,C_BCg,S_BCg,U_BCg,V_BCg,rrank_BCg,dN_x,dN_y,dN_g,dN_t,N_x,N_y,N_g,N_t,&
  t,PODgsum,Sim_Grid_Bnds,Dat_Grid_Bnds,Sim_TGrid,Dat_TGrid,GMap_xyBnds,VMap_xyBnds,TMap,Open_Threads)

  REAL*8,INTENT(OUT):: Cg_L(*), Cg_B(*), Cg_R(*), Cg_T(*)

  REAL*8,INTENT(IN):: C_BCg(*), S_BCg(*), U_BCg(*), V_BCg(*)
  REAL*8,INTENT(IN):: Sim_Grid_Bnds(*), Dat_Grid_Bnds(*), Sim_TGrid(*), Dat_TGrid(*)
  INTEGER,INTENT(IN):: rrank_BCg(*), GMap_xyBnds(*), VMap_xyBnds(*), TMap(*)
  INTEGER,INTENT(IN):: dN_x, dN_y, dN_g, dN_t, N_x, N_y, N_g, N_t, t, PODgsum
  INTEGER,INTENT(IN):: Open_Threads

  REAL*8,ALLOCATABLE:: BC(:), BC2(:), CL2(:), CB2(:), CR2(:), CT2(:)
  INTEGER:: len, g, i, j, p, pl, pb, pr, pt

  len = dN_g*2*(dN_x+dN_y) !length of cell-averaged data vector
  ALLOCATE(BC(len),BC2(len)) !allocating space to hold reconstituted POD data

  ! CALL POD_RECONSTRUCT(BC,C_BCg,S_BCg,U_BCg,V_BCg,rrank_BCg,len,dN_g,dN_t,t,PODgsum)
  IF (TMap(2*t-1) .EQ. TMap(2*t)) THEN
    CALL POD_RECONSTRUCT(BC,C_BCg,S_BCg,U_BCg,V_BCg,rrank_BCg,len,dN_g,dN_t,TMap(2*t),PODgsum)
  ELSE
    CALL POD_RECONSTRUCT(BC,C_BCg,S_BCg,U_BCg,V_BCg,rrank_BCg,len,dN_g,dN_t,TMap(2*t-1),PODgsum)
    CALL POD_RECONSTRUCT(BC2,C_BCg,S_BCg,U_BCg,V_BCg,rrank_BCg,len,dN_g,dN_t,TMap(2*t),PODgsum)
    CALL LINEAR_INTERPOLATION(BC,BC2,Sim_TGrid(t),Dat_TGrid(TMap(2*t-1)),Dat_TGrid(TMap(2*t)),len)
  END IF

  ALLOCATE(CL2(dN_g*(dN_y+2)), CB2(dN_g*(dN_x+2)), CR2(dN_g*(dN_y+2)), CT2(dN_g*(dN_x+2)))

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
  CALL FMAP(CB2, Cg_B, Dat_Grid_Bnds(2*dN_y+5), Sim_Grid_Bnds(2*N_y+1), GMap_xyBnds(4*N_y+1), VMap_xyBnds(4*N_y+1),&
  dN_x+2, N_x, dN_g*(dN_x+2), N_g*N_x, Open_Threads)

  !mapping the vector of ... data from the decomposition grid to the simulation grid
  CALL FMAP(CR2, Cg_R, Dat_Grid_Bnds(2*(dN_y+dN_x)+9), Sim_Grid_Bnds(2*(N_y+N_x)+1), GMap_xyBnds(4*(N_y+N_x)+1),&
  VMap_xyBnds(4*(N_y+N_x)+1), dN_y+2, N_y, dN_g*(dN_y+2), N_g*N_y, Open_Threads)

  !mapping the vector of ... data from the decomposition grid to the simulation grid
  CALL FMAP(CT2, Cg_T, Dat_Grid_Bnds(4*dN_y+2*dN_x+13), Sim_Grid_Bnds(4*N_y+2*N_x+1), GMap_xyBnds(8*N_y+4*N_x+1),&
  VMap_xyBnds(8*N_y+4*N_x+1), dN_x+2, N_x, dN_g*(dN_x+2), N_g*N_x, Open_Threads)

  DEALLOCATE(BC,BC2,CL2,CB2,CR2,CT2)


END SUBROUTINE POD_RECONSTRUCT_BCg

!==================================================================================================================================!
! SUBROUTINE POD_RECONSTRUCT
!
!DEF:
!  Reconstructs data at a specific instant of time from a truncated SVD representation
!
!INPUT:
!  N_t, N_g, N_y, N_x - size of the grid in time/ energy/ y-length/ x-length, respectively
!                     - Note that these are the grid sizes from the original decomposition, *NOT* the grid for the current simulation
!
!  PODgsum - flag to denote whether the given decomposition was performed over the entire phase space or groupwise for each qd factor
!
!  t - time step number to locate the desired instant of time to reconstruct data from (= Time/Delt)
!
!OUTPUT:
!
!==================================================================================================================================!
SUBROUTINE POD_RECONSTRUCT(dat,C,S,U,V,rank,len,N_g,N_t,t,PODgsum)
  REAL*8,INTENT(OUT):: dat(*)
  REAL*8,INTENT(IN):: C(*), S(*), U(*), V(*)
  INTEGER,INTENT(IN):: PODgsum, len, N_g, N_t, t, rank(*)
  INTEGER:: i, g, r, p1, p2, p3, p4, glen

  !initializing data array to zero
  DO i=1,len
    dat(i) = 0d0
  END DO

  !
  ! Case 1: groupwise decomposition of data
  !
  IF (PODgsum .EQ. 0) THEN
    glen = len/N_g

    p2 = 0
    p3 = 0
    DO g=1,N_g
      p4 = N_t*p3 + (t-1)*rank(g)
      DO r=1,rank(g)
        p1 = (g-1)*glen
        p4 = p4 + 1

        DO i=1,glen
          p1 = p1 + 1
          p2 = p2 + 1

          dat(p1) = dat(p1) + S(p3+r)*U(p2)*V(p4)

        END DO
      END DO

      p3 = p3 + rank(g)
    END DO

  !
  ! Case 2: decomposition of data over entire phase space
  !
  ELSE IF (PODgsum .EQ. 1) THEN

    p2 = 0
    DO r=1,rank(1)
      p3 = r + (t-1)*rank(1)
      DO i=1,len
        p2 = p2 + 1

        dat(i) = dat(i) + S(r)*U(p2)*V(p3)

      END DO
    END DO

  END IF

  !adding centering vector back to data
  DO i=1,len
    dat(i) = dat(i) + C(i)
  END DO

END SUBROUTINE POD_RECONSTRUCT

!==================================================================================================================================!
!
!==================================================================================================================================!
SUBROUTINE INPUT_fg_POD(Fname,PODgsum,eps,N_x,N_y,N_m,N_g,N_t,C_BCg,S_BCg,U_BCg,V_BCg,rrank_BCg,C_fg_avg_xx,S_fg_avg_xx,&
  U_fg_avg_xx,V_fg_avg_xx,rrank_fg_avg_xx,C_fg_edgV_xx,S_fg_edgV_xx,U_fg_edgV_xx,V_fg_edgV_xx,rrank_fg_edgV_xx,C_fg_avg_yy,&
  S_fg_avg_yy,U_fg_avg_yy,V_fg_avg_yy,rrank_fg_avg_yy,C_fg_edgH_yy,S_fg_edgH_yy,U_fg_edgH_yy,V_fg_edgH_yy,rrank_fg_edgH_yy,&
  C_fg_edgV_xy,S_fg_edgV_xy,U_fg_edgV_xy,V_fg_edgV_xy,rrank_fg_edgV_xy,C_fg_edgH_xy,S_fg_edgH_xy,U_fg_edgH_xy,V_fg_edgH_xy,&
  rrank_fg_edgH_xy,tlen,xlen,ylen,Tini,Delt,Delx,Dely,bcT,BC_Type,Start_Time)

  CHARACTER(*),INTENT(IN):: Fname
  INTEGER,INTENT(OUT):: PODgsum
  REAL*8,INTENT(IN):: eps
  REAL*8,INTENT(OUT):: tlen, xlen, ylen, Tini, Delt, Start_Time
  REAL*8,ALLOCATABLE,INTENT(OUT):: Delx(:), Dely(:), bcT(:)
  REAL*8,ALLOCATABLE,INTENT(OUT):: C_BCg(:), S_BCg(:), U_BCg(:), V_BCg(:)
  REAL*8,ALLOCATABLE,INTENT(OUT):: C_fg_avg_xx(:), S_fg_avg_xx(:), U_fg_avg_xx(:), V_fg_avg_xx(:)
  REAL*8,ALLOCATABLE,INTENT(OUT):: C_fg_edgV_xx(:), S_fg_edgV_xx(:), U_fg_edgV_xx(:), V_fg_edgV_xx(:)
  REAL*8,ALLOCATABLE,INTENT(OUT):: C_fg_avg_yy(:), S_fg_avg_yy(:), U_fg_avg_yy(:), V_fg_avg_yy(:)
  REAL*8,ALLOCATABLE,INTENT(OUT):: C_fg_edgH_yy(:), S_fg_edgH_yy(:), U_fg_edgH_yy(:), V_fg_edgH_yy(:)
  REAL*8,ALLOCATABLE,INTENT(OUT):: C_fg_edgV_xy(:), S_fg_edgV_xy(:), U_fg_edgV_xy(:), V_fg_edgV_xy(:)
  REAL*8,ALLOCATABLE,INTENT(OUT):: C_fg_edgH_xy(:), S_fg_edgH_xy(:), U_fg_edgH_xy(:), V_fg_edgH_xy(:)
  INTEGER,INTENT(OUT):: N_x, N_y, N_m, N_g, N_t
  INTEGER,ALLOCATABLE,INTENT(OUT):: rrank_fg_avg_xx(:), rrank_fg_edgV_xx(:), rrank_fg_avg_yy(:), rrank_fg_edgH_yy(:)
  INTEGER,ALLOCATABLE,INTENT(OUT):: rrank_fg_edgV_xy(:), rrank_fg_edgH_xy(:), rrank_BCg(:)
  INTEGER,ALLOCATABLE,INTENT(OUT):: BC_Type(:)

  CHARACTER(10):: dcmp_type
  CHARACTER(50):: Location = 'MODULE: POD_ROUTINES / SUBROUTINE: INPUT_fg_POD'
  INTEGER:: ncID, clen, g, Status, ID2
  REAL*8:: bnds(2)

  CALL NF_OPEN_FILE(ncID,Fname,'old','r')

  CALL NF_INQ_DIM(ncID,"N_x",N_x)
  CALL NF_INQ_DIM(ncID,"N_y",N_y)
  CALL NF_INQ_DIM(ncID,"N_g",N_g)
  CALL NF_INQ_DIM(ncID,"N_t",N_t)

  ALLOCATE(Delx(N_x),Dely(N_y),BC_Type(4),bcT(4))

  CALL NF_INQ_VAR_0D(ncID,"tlen",tlen)
  CALL NF_INQ_VAR_0D(ncID,"xlen",xlen)
  CALL NF_INQ_VAR_0D(ncID,"ylen",ylen)
  CALL NF_INQ_VAR_0D(ncID,"Delt",Delt)
  CALL NF_INQ_VAR_1D(ncID,"Delx",Delx,(/1/),(/N_x/))
  CALL NF_INQ_VAR_1D(ncID,"Dely",Dely,(/1/),(/N_y/))
  ! CALL NF_INQ_intVAR_1D(ncID,"BC_Type",BC_Type,(/1/),(/4/))
  ! CALL NF_INQ_VAR_1D(ncID,"bcT",bcT,(/1/),(/4/))
  CALL NF_INQ_VAR_0D(ncID,"bcT_left",bcT(1))
  CALL NF_INQ_VAR_0D(ncID,"bcT_bottom",bcT(2))
  CALL NF_INQ_VAR_0D(ncID,"bcT_right",bcT(3))
  CALL NF_INQ_VAR_0D(ncID,"bcT_top",bcT(4))
  CALL NF_INQ_VAR_0D(ncID,"Tini",Tini)

  Status = nf90_inq_varid(ncID, "tpts", ID2)
  CALL HANDLE_ERR(Status, Location)
  Status = nf90_get_att(ncID, ID2, "bnds", bnds)
  CALL HANDLE_ERR(Status, Location)
  Start_Time = bnds(1) - Delt

  Status = nf90_get_att(ncID, NF90_GLOBAL, "dcmp_type", dcmp_type)
  CALL HANDLE_ERR(Status,Location)
  IF (dcmp_type(1:4) .EQ. "PODg") THEN
    PODgsum = 0
  ELSE IF (dcmp_type(1:3) .EQ. "POD") THEN
    PODgsum = 1
  ELSE
    WRITE(*,'(5A)') "POD_ROUTINES :: INPUT_fg_POD - Unrecognized dcmp_type (",trim(dcmp_type),") in dataset (",trim(Fname),")"
    STOP
  END IF

  IF (PODgsum .EQ. 1) THEN
    ALLOCATE(rrank_fg_avg_xx(1), rrank_fg_edgV_xx(1), rrank_fg_avg_yy(1), rrank_fg_edgH_yy(1))
    ALLOCATE(rrank_fg_edgV_xy(1), rrank_fg_edgH_xy(1), rrank_BCg(1))

    clen = N_g*2*(N_x+N_y)
    CALL READ_POD(ncID,'BCg',clen,N_t,eps,C_BCg,S_BCg,U_BCg,V_BCg,rrank_BCg(1))

    clen = N_x*N_y*N_g
    CALL READ_POD(ncID,'fg_avg_xx',clen,N_t,eps,C_fg_avg_xx,S_fg_avg_xx,U_fg_avg_xx,V_fg_avg_xx,rrank_fg_avg_xx(1))

    clen = (N_x+1)*N_y*N_g
    CALL READ_POD(ncID,'fg_edgV_xx',clen,N_t,eps,C_fg_edgV_xx,S_fg_edgV_xx,U_fg_edgV_xx,V_fg_edgV_xx,rrank_fg_edgV_xx(1))

    clen = N_x*N_y*N_g
    CALL READ_POD(ncID,'fg_avg_yy',clen,N_t,eps,C_fg_avg_yy,S_fg_avg_yy,U_fg_avg_yy,V_fg_avg_yy,rrank_fg_avg_yy(1))

    clen = N_x*(N_y+1)*N_g
    CALL READ_POD(ncID,'fg_edgH_yy',clen,N_t,eps,C_fg_edgH_yy,S_fg_edgH_yy,U_fg_edgH_yy,V_fg_edgH_yy,rrank_fg_edgH_yy(1))

    clen = (N_x+1)*N_y*N_g
    CALL READ_POD(ncID,'fg_edgV_xy',clen,N_t,eps,C_fg_edgV_xy,S_fg_edgV_xy,U_fg_edgV_xy,V_fg_edgV_xy,rrank_fg_edgV_xy(1))

    clen = N_x*(N_y+1)*N_g
    CALL READ_POD(ncID,'fg_edgH_xy',clen,N_t,eps,C_fg_edgH_xy,S_fg_edgH_xy,U_fg_edgH_xy,V_fg_edgH_xy,rrank_fg_edgH_xy(1))

  ELSE
    ALLOCATE(rrank_fg_avg_xx(N_g), rrank_fg_edgV_xx(N_g), rrank_fg_avg_yy(N_g), rrank_fg_edgH_yy(N_g))
    ALLOCATE(rrank_fg_edgV_xy(N_g), rrank_fg_edgH_xy(N_g), rrank_BCg(N_g))

    DO g=1,N_g
      clen = 2*(N_x+N_y)
      CALL READ_POD(ncID,'BCg',clen,N_t,N_g,g,eps,C_BCg,S_BCg,U_BCg,V_BCg,rrank_BCg(g))

      clen = N_x*N_y
      CALL READ_POD(ncID,'fg_avg_xx',clen,N_t,N_g,g,eps,C_fg_avg_xx,S_fg_avg_xx,U_fg_avg_xx,V_fg_avg_xx,rrank_fg_avg_xx(g))

      clen = (N_x+1)*N_y
      CALL READ_POD(ncID,'fg_edgV_xx',clen,N_t,N_g,g,eps,C_fg_edgV_xx,S_fg_edgV_xx,U_fg_edgV_xx,V_fg_edgV_xx,rrank_fg_edgV_xx(g))

      clen = N_x*N_y
      CALL READ_POD(ncID,'fg_avg_yy',clen,N_t,N_g,g,eps,C_fg_avg_yy,S_fg_avg_yy,U_fg_avg_yy,V_fg_avg_yy,rrank_fg_avg_yy(g))

      clen = N_x*(N_y+1)
      CALL READ_POD(ncID,'fg_edgH_yy',clen,N_t,N_g,g,eps,C_fg_edgH_yy,S_fg_edgH_yy,U_fg_edgH_yy,V_fg_edgH_yy,rrank_fg_edgH_yy(g))

      clen = (N_x+1)*N_y
      CALL READ_POD(ncID,'fg_edgV_xy',clen,N_t,N_g,g,eps,C_fg_edgV_xy,S_fg_edgV_xy,U_fg_edgV_xy,V_fg_edgV_xy,rrank_fg_edgV_xy(g))

      clen = N_x*(N_y+1)
      CALL READ_POD(ncID,'fg_edgH_xy',clen,N_t,N_g,g,eps,C_fg_edgH_xy,S_fg_edgH_xy,U_fg_edgH_xy,V_fg_edgH_xy,rrank_fg_edgH_xy(g))

    END DO

  END IF

  CALL NF_CLOSE_FILE(ncID)

END SUBROUTINE INPUT_fg_POD

!==================================================================================================================================!
!
!==================================================================================================================================!
SUBROUTINE INPUT_Ig_POD(Fname,PODgsum,N_x,N_y,N_m,N_g,N_t,eps,C_I_avg,S_I_avg,U_I_avg,V_I_avg,rrank_I_avg,&
  C_I_edgV,S_I_edgV,U_I_edgV,V_I_edgV,rrank_I_edgV,C_I_edgH,S_I_edgH,U_I_edgH,V_I_edgH,rrank_I_edgH)

  CHARACTER(*),INTENT(IN):: Fname
  INTEGER,INTENT(IN):: N_x, N_y, N_m, N_g, N_t, PODgsum
  REAL*8,INTENT(IN):: eps
  REAL*8,ALLOCATABLE,INTENT(OUT):: C_I_avg(:), S_I_avg(:), U_I_avg(:), V_I_avg(:)
  REAL*8,ALLOCATABLE,INTENT(OUT):: C_I_edgV(:), S_I_edgV(:), U_I_edgV(:), V_I_edgV(:)
  REAL*8,ALLOCATABLE,INTENT(OUT):: C_I_edgH(:), S_I_edgH(:), U_I_edgH(:), V_I_edgH(:)
  INTEGER,ALLOCATABLE,INTENT(OUT):: rrank_I_avg(:), rrank_I_edgV(:), rrank_I_edgH(:)
  INTEGER:: ncID, dN_x, dN_y, dN_m, dN_g, dN_t, clen, g

  CALL NF_OPEN_FILE(ncID,Fname,'old','r')

  CALL NF_INQ_DIM(ncID,"N_x",dN_x)
  CALL NF_INQ_DIM(ncID,"N_y",dN_y)
  CALL NF_INQ_DIM(ncID,"N_m",dN_m)
  CALL NF_INQ_DIM(ncID,"N_g",dN_g)
  CALL NF_INQ_DIM(ncID,"N_t",dN_t)

  IF (dN_x .NE. N_x) THEN
    WRITE(*,'(A)') 'N_x from dataset does not match N_x for problem'
    WRITE(*,'(A,I3.3)') 'N_x (dataset) = ',dN_x
    WRITE(*,'(A,I3.3)') 'N_x (problem) = ',N_x
    STOP
  ELSE IF (dN_y .NE. N_y) THEN
    WRITE(*,'(A)') 'N_y from dataset does not match N_y for problem'
    WRITE(*,'(A,I3.3)') 'N_y (dataset) = ',dN_y
    WRITE(*,'(A,I3.3)') 'N_y (problem) = ',N_y
    STOP
  ELSE IF (dN_g .NE. N_g) THEN
    WRITE(*,'(A)') 'N_g from dataset does not match N_g for problem'
    WRITE(*,'(A,I3.3)') 'N_g (dataset) = ',dN_g
    WRITE(*,'(A,I3.3)') 'N_g (problem) = ',N_g
    STOP
  ELSE IF (dN_m .NE. N_m) THEN
    WRITE(*,'(A)') 'N_m from dataset does not match N_m for problem'
    WRITE(*,'(A,I3.3)') 'N_m (dataset) = ',dN_m
    WRITE(*,'(A,I3.3)') 'N_m (problem) = ',N_m
    STOP
  END IF

  IF (PODgsum .EQ. 1) THEN
    ALLOCATE(rrank_I_avg(1), rrank_I_edgV(1), rrank_I_edgH(1))

    clen = dN_x*dN_y*dN_m*dN_g
    CALL READ_POD(ncID,'I_avg',clen,dN_t,eps,C_I_avg,S_I_avg,U_I_avg,V_I_avg,rrank_I_avg(1))

    clen = (dN_x+1)*dN_y*dN_m*dN_g
    CALL READ_POD(ncID,'I_edgV',clen,dN_t,eps,C_I_edgV,S_I_edgV,U_I_edgV,V_I_edgV,rrank_I_edgV(1))

    clen = dN_x*(dN_y+1)*dN_m*dN_g
    CALL READ_POD(ncID,'I_edgH',clen,dN_t,eps,C_I_edgH,S_I_edgH,U_I_edgH,V_I_edgH,rrank_I_edgH(1))

  ELSE
    ALLOCATE(rrank_I_avg(dN_g), rrank_I_edgV(dN_g), rrank_I_edgH(dN_g))

    DO g=1,dN_g
      clen = dN_x*dN_y*dN_m
      CALL READ_POD(ncID,'I_avg',clen,dN_t,dN_g,g,eps,C_I_avg,S_I_avg,U_I_avg,V_I_avg,rrank_I_avg(g))

      clen = (dN_x+1)*dN_y*dN_m
      CALL READ_POD(ncID,'I_edgV',clen,dN_t,dN_g,g,eps,C_I_edgV,S_I_edgV,U_I_edgV,V_I_edgV,rrank_I_edgV(g))

      clen = dN_x*(dN_y+1)*dN_m
      CALL READ_POD(ncID,'I_edgH',clen,dN_t,dN_g,g,eps,C_I_edgH,S_I_edgH,U_I_edgH,V_I_edgH,rrank_I_edgH(g))

    END DO

  END IF

  CALL NF_CLOSE_FILE(ncID)

END SUBROUTINE INPUT_Ig_POD

!==================================================================================================================================!
!Subroutine READ_POD
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
SUBROUTINE READ_POD_FULL(ncID,name,clen,N_t,eps,C,S,U,V,rrank)
  REAL*8,INTENT(IN):: eps
  INTEGER,INTENT(IN):: ncID, clen, N_t
  CHARACTER(*),INTENT(IN):: name
  REAL*8,ALLOCATABLE,INTENT(OUT):: C(:), S(:), U(:), V(:)
  INTEGER,INTENT(OUT):: rrank
  REAL*8,ALLOCATABLE:: Sn(:)
  INTEGER:: rank
  CHARACTER(100):: name2

  CALL CLEAR_POD(C,S,U,V)

  !Determining rank of data matrix to be read in
  rank = MIN(clen,N_t)

  !Reading the centering vector
  ALLOCATE(C(clen))
  WRITE(name2,'(2A)') 'C_',name
  CALL NF_INQ_VAR_1D(ncID,name2,C,(/1/),(/clen/))

  !Reading the vector of singular values
  ALLOCATE(S(rank),Sn(rank))
  WRITE(name2,'(2A)') 'S_',name
  CALL NF_INQ_VAR_1D(ncID,name2,S,(/1/),(/rank/))

  !Normalizing singular values; calculating reduced rank
  CALL RRANK_CALC(S,eps,rank,rrank)

  !Resizing vector of singular values to hold only the first 'rrank' values
  DEALLOCATE(Sn)
  ALLOCATE(Sn(rrank))
  Sn = S(1:rrank)
  DEALLOCATE(S)
  ALLOCATE(S(rrank))
  S = Sn
  DEALLOCATE(Sn)

  !Reading singular vector matrices U and V (only first rrank vectors in each matrix)
  ALLOCATE(U(clen*rrank),V(rrank*N_t))
  WRITE(name2,'(2A)') 'U_',name
  CALL NF_INQ_VAR_1D(ncID,name2,U,(/1,1/),(/clen,rrank/))
  WRITE(name2,'(2A)') 'Vt_',name
  CALL NF_INQ_VAR_1D(ncID,name2,V,(/1,1/),(/rrank,N_t/))

END SUBROUTINE READ_POD_FULL

SUBROUTINE READ_POD_GPART(ncID,name,clen,N_t,N_g,g,eps,C,S,U,V,rrank)
  REAL*8,INTENT(IN):: eps
  INTEGER,INTENT(IN):: ncID, clen, N_t, N_g, g
  CHARACTER(*),INTENT(IN):: name
  REAL*8,ALLOCATABLE,INTENT(INOUT):: C(:), S(:), U(:), V(:)
  INTEGER,INTENT(OUT):: rrank
  REAL*8,ALLOCATABLE:: Sn(:), S2(:), U2(:), V2(:)
  INTEGER:: rank
  CHARACTER(100):: name2

  IF (g.EQ.1) CALL CLEAR_POD(C,S,U,V)

  !Determining rank of data matrix to be read in
  rank = MIN(clen,N_t)

  !Reading the centering vector
  IF (.NOT. ALLOCATED(C)) ALLOCATE(C(N_g*clen))
  WRITE(name2,'(2A)') 'C_',name
  CALL NF_INQ_VAR_1D(ncID,name2,C((g-1)*clen+1:g*clen),(/1,g/),(/clen,1/))

  !Reading the vector of singular values
  ALLOCATE(S2(rank),Sn(rank))
  rrank = rank
  WRITE(name2,'(2A)') 'S_',name
  CALL NF_INQ_VAR_1D(ncID,name2,S2,(/1,g/),(/rank,1/))

  !Normalizing singular values; calculating reduced rank
  CALL RRANK_CALC(S2,eps,rank,rrank)

  !Resizing vector of singular values to hold only the first 'rrank' values
  DEALLOCATE(Sn)
  IF (ALLOCATED(S)) THEN
    ALLOCATE(Sn(SIZE(S)))
    Sn = S
    DEALLOCATE(S)
    ALLOCATE(S(SIZE(Sn)+rrank))
    S(1:SIZE(Sn)) = Sn
    S(SIZE(Sn)+1:SIZE(S)) = S2(1:rrank)
    DEALLOCATE(S2, Sn)
  ELSE
    ALLOCATE(S(rrank))
    S(1:rrank) = S2(1:rrank)
    DEALLOCATE(S2)
  END IF

  !Resizing matrices of singular vectors (U, V) to add on new vectors (number of new vectors = rrank)
  IF (ALLOCATED(U)) THEN
    ALLOCATE(U2(SIZE(U)))
    U2 = U
    DEALLOCATE(U)
    ALLOCATE(U(SIZE(U2)+clen*rrank))
    U(1:SIZE(U2)) = U2
    DEALLOCATE(U2)
  ELSE
    ALLOCATE(U(clen*rrank))
  END IF

  IF (ALLOCATED(V)) THEN
    ALLOCATE(V2(SIZE(V)))
    V2 = V
    DEALLOCATE(V)
    ALLOCATE(V(SIZE(V2)+N_t*rrank))
    V(1:SIZE(V2)) = V2
    DEALLOCATE(V2)
  ELSE
    ALLOCATE(V(N_t*rrank))
  END IF

  !Reading singular vector matrices U and V (only first rrank vectors in each matrix)
  WRITE(name2,'(2A)') 'U_',name
  CALL NF_INQ_VAR_1D(ncID,name2,U(SIZE(U)-clen*rrank+1:SIZE(U)),(/1,1,g/),(/clen,rrank,1/))
  WRITE(name2,'(2A)') 'Vt_',name
  CALL NF_INQ_VAR_1D(ncID,name2,V(SIZE(V)-N_t*rrank+1:SIZE(V)),(/1,1,g/),(/rrank,N_t,1/))

END SUBROUTINE READ_POD_GPART

!==================================================================================================================================!
!Subroutine RRANK_CALC
!
! Finds the reduced rank satisfying a desired POD error level
!
! NOTES:: If one wants the reduced rank to satisfy the 'relative' POD error, simply normalize the vector of singular vlues prior to
!         calling this subroutine such that sig(1)=1
!
! WARNINGS::
!
! OUTPUTS::
!
! INPUTS::
!
!==================================================================================================================================!
SUBROUTINE RRANK_CALC(sig,eps,len,rank)
  REAL*8,INTENT(IN):: sig(*), eps
  INTEGER,INTENT(IN):: len
  INTEGER,INTENT(OUT):: rank
  REAL*8:: sum1, sum2, eps2
  INTEGER:: i, j

  sum1 = 0d0
  DO i=1,len
    sum1 = sum1 + sig(i)**2
  END DO
  eps2 = eps*SQRT(sum1)

  rank = len
  DO i=1,len
    sum2 = 0d0
    DO j=i+1,len
      sum2 = sum2 + sig(j)**2
    END DO
    sum2 = SQRT(sum2)

    IF (sum2 .LE. eps2) THEN
      rank = i
      EXIT
    END IF

  END DO

END SUBROUTINE RRANK_CALC

!==================================================================================================================================!
!
!==================================================================================================================================!
SUBROUTINE CLEAR_POD(C,S,U,V)
  REAL*8,ALLOCATABLE,INTENT(INOUT):: C(:), S(:), U(:), V(:)
  IF (ALLOCATED(C)) DEALLOCATE(C)
  IF (ALLOCATED(S)) DEALLOCATE(S)
  IF (ALLOCATED(U)) DEALLOCATE(U)
  IF (ALLOCATED(V)) DEALLOCATE(V)
END SUBROUTINE CLEAR_POD
!==================================================================================================================================!
!
!==================================================================================================================================!

END MODULE POD_ROUTINES
