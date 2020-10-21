MODULE POD_ROUTINES

  USE NCDF_IO

  IMPLICIT NONE

CONTAINS

!==================================================================================================================================!
!
!==================================================================================================================================!
SUBROUTINE INPUT_fg_POD(Fname,N_x,N_y,N_g,N_t,eps,C_fg_avg_xx,S_fg_avg_xx,U_fg_avg_xx,V_fg_avg_xx,rrank_fg_avg_xx)
  CHARACTER(*),INTENT(IN):: Fname
  INTEGER,INTENT(IN):: N_x, N_y, N_g, N_t
  REAL*8,INTENT(IN):: eps
  REAL*8,ALLOCATABLE,INTENT(OUT):: C_fg_avg_xx(:), S_fg_avg_xx(:), U_fg_avg_xx(:,:), V_fg_avg_xx(:,:)
  INTEGER,INTENT(OUT):: rrank_fg_avg_xx
  INTEGER:: ncID, dN_x, dN_y, dN_g, dN_t, clen

  CALL NF_OPEN_FILE(ncID,Fname,'old','r')

  CALL NF_INQ_DIM(ncID,"N_x",dN_x)
  CALL NF_INQ_DIM(ncID,"N_y",dN_y)
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
  END IF

  clen = dN_x*dN_y*dN_g
  CALL READ_POD(ncID,'fg_avg_xx',clen,dN_t,eps,C_fg_avg_xx,S_fg_avg_xx,U_fg_avg_xx,V_fg_avg_xx,rrank_fg_avg_xx)

  CALL NF_CLOSE_FILE(ncID)

END SUBROUTINE INPUT_fg_POD

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
SUBROUTINE READ_POD(ncID,name,clen,N_t,eps,C,S,U,V,rrank)
  REAL*8,INTENT(IN):: eps
  INTEGER,INTENT(IN):: ncID, clen, N_t
  CHARACTER(*),INTENT(IN):: name
  REAL*8,ALLOCATABLE,INTENT(OUT):: C(:), S(:), U(:,:), V(:,:)
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
  ALLOCATE(U(clen,rrank),V(rrank,N_t))
  WRITE(name2,'(2A)') 'U_',name
  CALL NF_INQ_VAR_2D(ncID,name2,U,(/1,1/),(/clen,rrank/))
  WRITE(name2,'(2A)') 'Vt_',name
  CALL NF_INQ_VAR_2D(ncID,name2,V,(/1,1/),(/rrank,N_t/))

END SUBROUTINE READ_POD

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
