MODULE POD_Routines

  USE ISO_C_BINDING
  IMPLICIT NONE

CONTAINS

!==================================================================================================================================!
!
!==================================================================================================================================!
SUBROUTINE SVD_CALC(dat,N_t,N_y,N_x,center,umat,sig,vtmat) BIND(c, name="SVD_CALC")
  INTEGER,INTENT(IN):: N_t, N_y, N_x
  REAL(c_double),INTENT(IN):: dat(*)
  REAL(c_double),INTENT(OUT):: center(*), umat(*), sig(*), vtmat(*)

  REAL*8,ALLOCATABLE:: mat(:,:), u(:,:), s(:), vt(:,:)
  REAL*8:: rank
  INTEGER:: t, i, j, r

  rank = MIN(N_x*N_y,N_t)
  ALLOCATE(mat(N_x*N_y,N_t))
  ALLOCATE(u(N_x*N_y,N_x*N_y))
  ALLOCATE(vt(N_t,N_t))
  ALLOCATE(s(rank))

  DO t=1,N_t
    DO j=1,N_y
      DO i=1,N_x
        center(i+(j-1)*N_x) = dat(i+(j-1)*N_x+(t-1)*N_x*N_y)
      END DO
    END DO
  END DO

  DO j=1,N_y
    DO i=1,N_x
      center(i+(j-1)*N_x) = center(i+(j-1)*N_x)/DBLE(N_t)
    END DO
  END DO

  DO t=1,N_t
    DO j=1,N_y
      DO i=1,N_x
        mat(i+(j-1)*N_x,t) = dat(i+(j-1)*N_x+(t-1)*N_x*N_y) - center(i+(j-1)*N_x)
      END DO
    END DO
  END DO

  CALL SVD_DECOMP(mat,u,s,vt)

  DO r=1,rank
    DO j=1,N_y
      DO i=1,N_x
        umat(i+(j-1)*N_x+(r-1)*N_x*N_y) = u(i+(j-1)*N_x,r)
      END DO
    END DO
  END DO

  DO t=1,N_t
    DO r=1,rank
      vtmat(r+(t-1)*rank) = vt(r,t)
    END DO
  END DO

  DO r=1,rank
    sig(r) = s(r)
  END DO

  DEALLOCATE(mat,u,vt,s)

END SUBROUTINE SVD_CALC

!==================================================================================================================================!
!Subroutine svd_decomp
!
! this subroutine will take any matrix and decompose it with svd such that A = U * Sig * V^T
! where A is the matrix to be decomposed
! !!! REQUIRES LAPACK !!!
! WARNINGS:: WILL DESTROY THE INPUT MATRIX datamat
!
! INPUTS:: datamat, dim1, dim2
!   datamat (REAL) -- matrix to be decomposed by SVD, dimensions (dim1 x dim2)
!
! OUTPUTS:: umat, sig, vtmat
!   umat (REAL) -- matrix of dimension(dim1 x dim1), the 'left singular vectors' of datamat are the columns
!   sig (REAL) -- vector of dimension(MIN(dim1,dim2)), the 'singular values' of datamat
!   vtmat (REAL) -- matrix of dimension(dim2 x dim2), the 'right singular vectors' of datamat are the rows
!
! INTERNALS:: dim1, dim2, LWORK, info, WORK
!   dim1 (INT) -- first dimension of 'datamat' (rows)
!   dim2 (INT) -- second dimension of 'datamat' (columns)
!
!==================================================================================================================================!
SUBROUTINE SVD_DECOMP(datamat,umat,sig,vtmat)
  INTEGER:: dim1, dim2
  INTEGER:: LWORK, info
  REAL*8,INTENT(INOUT):: datamat(:,:)
  REAL*8,INTENT(OUT):: umat(:,:), sig(:), vtmat(:,:)
  REAL*8,ALLOCATABLE:: WORK(:)

  dim1 = SIZE(datamat,1)
  dim2 = SIZE(datamat,2)
  LWORK=MAX(1,3*MIN(dim1,dim2) + MAX(dim1,dim2),5*MIN(dim1,dim2)) !honestly still dont know what the fuck this is but it works
  ALLOCATE(WORK(LWORK))

  CALL DGESVD('A','A',dim1,dim2,datamat,dim1,sig,umat,dim1,vtmat,dim2,WORK,LWORK,info)  !lapack SVD function
  IF (info .NE. 0) STOP "SVD failed"

  DEALLOCATE(WORK)

END SUBROUTINE SVD_DECOMP

END MODULE POD_Routines
