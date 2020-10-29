MODULE POD_ROUTINES

  USE NCDF_IO

  IMPLICIT NONE

INTERFACE READ_POD
  MODULE PROCEDURE READ_POD_FULL, READ_POD_GPART
END INTERFACE

CONTAINS

!==================================================================================================================================!
!
!==================================================================================================================================!
SUBROUTINE POD_RECONSTRUCT_fg(fg_avg_xx,fg_avg_yy,fg_edgV_xx,fg_edgV_xy,fg_edgH_yy,fg_edgH_xy,C_fg_avg_xx,&
  S_fg_avg_xx,U_fg_avg_xx,V_fg_avg_xx,C_fg_edgV_xx,S_fg_edgV_xx,U_fg_edgV_xx,V_fg_edgV_xx,C_fg_avg_yy,S_fg_avg_yy,&
  U_fg_avg_yy,V_fg_avg_yy,C_fg_edgH_yy,S_fg_edgH_yy,U_fg_edgH_yy,V_fg_edgH_yy,C_fg_edgV_xy,S_fg_edgV_xy,U_fg_edgV_xy,&
  V_fg_edgV_xy,C_fg_edgH_xy,S_fg_edgH_xy,U_fg_edgH_xy,V_fg_edgH_xy,rrank_fg_avg_xx,rrank_fg_edgV_xx,rrank_fg_avg_yy,&
  rrank_fg_edgH_yy,rrank_fg_edgV_xy,rrank_fg_edgH_xy,N_x,N_y,N_g,N_t,t,PODgsum)

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
  INTEGER,INTENT(IN):: N_x, N_y, N_g, N_t, t, PODgsum

  INTEGER:: len

  len = N_x*N_y*N_g
  CALL POD_RECONSTRUCT(fg_avg_xx,C_fg_avg_xx,S_fg_avg_xx,U_fg_avg_xx,V_fg_avg_xx,rrank_fg_avg_xx,len,N_g,N_t,t,PODgsum)
  CALL POD_RECONSTRUCT(fg_avg_yy,C_fg_avg_yy,S_fg_avg_yy,U_fg_avg_yy,V_fg_avg_yy,rrank_fg_avg_yy,len,N_g,N_t,t,PODgsum)

  len = (N_x+1)*N_y*N_g
  CALL POD_RECONSTRUCT(fg_edgV_xx,C_fg_edgV_xx,S_fg_edgV_xx,U_fg_edgV_xx,V_fg_edgV_xx,rrank_fg_edgV_xx,len,N_g,N_t,t,PODgsum)
  CALL POD_RECONSTRUCT(fg_edgV_xy,C_fg_edgV_xy,S_fg_edgV_xy,U_fg_edgV_xy,V_fg_edgV_xy,rrank_fg_edgV_xy,len,N_g,N_t,t,PODgsum)

  len = N_x*(N_y+1)*N_g
  CALL POD_RECONSTRUCT(fg_edgH_yy,C_fg_edgH_yy,S_fg_edgH_yy,U_fg_edgH_yy,V_fg_edgH_yy,rrank_fg_edgH_yy,len,N_g,N_t,t,PODgsum)
  CALL POD_RECONSTRUCT(fg_edgH_xy,C_fg_edgH_xy,S_fg_edgH_xy,U_fg_edgH_xy,V_fg_edgH_xy,rrank_fg_edgH_xy,len,N_g,N_t,t,PODgsum)

END SUBROUTINE POD_RECONSTRUCT_fg

!==================================================================================================================================!
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
SUBROUTINE INPUT_fg_POD(Fname,PODgsum,eps,N_x,N_y,N_m,N_g,N_t,C_fg_avg_xx,S_fg_avg_xx,U_fg_avg_xx,V_fg_avg_xx,rrank_fg_avg_xx,&
  C_fg_edgV_xx,S_fg_edgV_xx,U_fg_edgV_xx,V_fg_edgV_xx,rrank_fg_edgV_xx,C_fg_avg_yy,S_fg_avg_yy,U_fg_avg_yy,V_fg_avg_yy,&
  rrank_fg_avg_yy,C_fg_edgH_yy,S_fg_edgH_yy,U_fg_edgH_yy,V_fg_edgH_yy,rrank_fg_edgH_yy,C_fg_edgV_xy,S_fg_edgV_xy,U_fg_edgV_xy,&
  V_fg_edgV_xy,rrank_fg_edgV_xy,C_fg_edgH_xy,S_fg_edgH_xy,U_fg_edgH_xy,V_fg_edgH_xy,rrank_fg_edgH_xy,tlen,xlen,ylen,Tini,Delt,&
  Delx,Dely,bcT,BC_Type)

  CHARACTER(*),INTENT(IN):: Fname
  INTEGER,INTENT(IN):: PODgsum
  REAL*8,INTENT(IN):: eps
  REAL*8,INTENT(OUT):: tlen, xlen, ylen, Tini, Delt
  REAL*8,ALLOCATABLE,INTENT(OUT):: Delx(:), Dely(:), bcT(:)
  REAL*8,ALLOCATABLE,INTENT(OUT):: C_fg_avg_xx(:), S_fg_avg_xx(:), U_fg_avg_xx(:), V_fg_avg_xx(:)
  REAL*8,ALLOCATABLE,INTENT(OUT):: C_fg_edgV_xx(:), S_fg_edgV_xx(:), U_fg_edgV_xx(:), V_fg_edgV_xx(:)
  REAL*8,ALLOCATABLE,INTENT(OUT):: C_fg_avg_yy(:), S_fg_avg_yy(:), U_fg_avg_yy(:), V_fg_avg_yy(:)
  REAL*8,ALLOCATABLE,INTENT(OUT):: C_fg_edgH_yy(:), S_fg_edgH_yy(:), U_fg_edgH_yy(:), V_fg_edgH_yy(:)
  REAL*8,ALLOCATABLE,INTENT(OUT):: C_fg_edgV_xy(:), S_fg_edgV_xy(:), U_fg_edgV_xy(:), V_fg_edgV_xy(:)
  REAL*8,ALLOCATABLE,INTENT(OUT):: C_fg_edgH_xy(:), S_fg_edgH_xy(:), U_fg_edgH_xy(:), V_fg_edgH_xy(:)
  INTEGER,INTENT(OUT):: N_x, N_y, N_m, N_g, N_t
  INTEGER,ALLOCATABLE,INTENT(OUT):: rrank_fg_avg_xx(:), rrank_fg_edgV_xx(:), rrank_fg_avg_yy(:), rrank_fg_edgH_yy(:)
  INTEGER,ALLOCATABLE,INTENT(OUT):: rrank_fg_edgV_xy(:), rrank_fg_edgH_xy(:)
  INTEGER,ALLOCATABLE,INTENT(OUT):: BC_Type(:)

  INTEGER:: ncID, clen, g

  CALL NF_OPEN_FILE(ncID,Fname,'old','r')

  CALL NF_INQ_DIM(ncID,"N_x",N_x)
  CALL NF_INQ_DIM(ncID,"N_y",N_y)
  CALL NF_INQ_DIM(ncID,"N_m",N_m)
  CALL NF_INQ_DIM(ncID,"N_g",N_g)
  CALL NF_INQ_DIM(ncID,"N_t",N_t)

  ALLOCATE(Delx(N_x),Dely(N_y),BC_Type(4),bcT(4))

  CALL NF_INQ_VAR_0D(ncID,"tlen",tlen)
  CALL NF_INQ_VAR_0D(ncID,"xlen",xlen)
  CALL NF_INQ_VAR_0D(ncID,"ylen",ylen)
  CALL NF_INQ_VAR_0D(ncID,"Delt",Delt)
  CALL NF_INQ_VAR_1D(ncID,"Delx",Delx,(/1/),(/N_x/))
  CALL NF_INQ_VAR_1D(ncID,"Dely",Dely,(/1/),(/N_y/))
  CALL NF_INQ_intVAR_1D(ncID,"BC_Type",BC_Type,(/1/),(/4/))
  CALL NF_INQ_VAR_1D(ncID,"bcT",bcT,(/1/),(/4/))
  CALL NF_INQ_VAR_0D(ncID,"Tini",Tini)

  IF (PODgsum .EQ. 1) THEN
    ALLOCATE(rrank_fg_avg_xx(1), rrank_fg_edgV_xx(1), rrank_fg_avg_yy(1), rrank_fg_edgH_yy(1))
    ALLOCATE(rrank_fg_edgV_xy(1), rrank_fg_edgH_xy(1))

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
    ALLOCATE(rrank_fg_edgV_xy(N_g), rrank_fg_edgH_xy(N_g))

    DO g=1,N_g
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
    CALL READ_POD(ncID,'Ig_avg',clen,dN_t,eps,C_I_avg,S_I_avg,U_I_avg,V_I_avg,rrank_I_avg(1))

    clen = (dN_x+1)*dN_y*dN_m*dN_g
    CALL READ_POD(ncID,'Ig_edgV',clen,dN_t,eps,C_I_edgV,S_I_edgV,U_I_edgV,V_I_edgV,rrank_I_edgV(1))

    clen = dN_x*(dN_y+1)*dN_m*dN_g
    CALL READ_POD(ncID,'Ig_edgH',clen,dN_t,eps,C_I_edgH,S_I_edgH,U_I_edgH,V_I_edgH,rrank_I_edgH(1))

  ELSE
    ALLOCATE(rrank_I_avg(dN_g), rrank_I_edgV(dN_g), rrank_I_edgH(dN_g))

    DO g=1,dN_g
      clen = dN_x*dN_y*dN_m
      CALL READ_POD(ncID,'Ig_avg',clen,dN_t,dN_g,g,eps,C_I_avg,S_I_avg,U_I_avg,V_I_avg,rrank_I_avg(g))

      clen = (dN_x+1)*dN_y*dN_m
      CALL READ_POD(ncID,'Ig_edgV',clen,dN_t,dN_g,g,eps,C_I_edgV,S_I_edgV,U_I_edgV,V_I_edgV,rrank_I_edgV(g))

      clen = dN_x*(dN_y+1)*dN_m
      CALL READ_POD(ncID,'Ig_edgH',clen,dN_t,dN_g,g,eps,C_I_edgH,S_I_edgH,U_I_edgH,V_I_edgH,rrank_I_edgH(g))

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
  Sn = S/S(1)
  CALL RRANK_CALC(Sn,eps,rank,rrank)

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
  Sn = S2/S2(1)
  CALL RRANK_CALC(Sn,eps,rank,rrank)

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
  REAL*8:: sum
  INTEGER:: i, j

  rank = len
  DO i=1,len
    sum = 0d0
    DO j=i+1,len
      sum = sum + sig(j)**2
    END DO
    sum = SQRT(sum)

    IF (sum .LE. eps) THEN
      rank = i
      EXIT
    END IF

  END DO

END SUBROUTINE RRANK_CALC

!==================================================================================================================================!
!
!==================================================================================================================================!

END MODULE POD_ROUTINES
