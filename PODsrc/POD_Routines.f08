MODULE POD_Routines

  USE ISO_C_BINDING
  IMPLICIT NONE

CONTAINS

!==================================================================================================================================!
!
!==================================================================================================================================!
SUBROUTINE SVD_CALC(dat,N_t,clen,center,umat,sig,vtmat) BIND(c, name="SVD_CALC")
  INTEGER(c_size_t),INTENT(IN):: N_t, clen
  REAL(c_double),INTENT(IN):: dat(*)
  REAL(c_double),INTENT(OUT):: center(*), umat(*), sig(*), vtmat(*)

  REAL*8,ALLOCATABLE:: mat(:,:), u(:,:), s(:), vt(:,:)
  INTEGER(c_size_t):: rank
  INTEGER(c_size_t):: t, i, j, r, p1

  rank = MIN(clen,N_t)
  ALLOCATE(mat(clen,N_t))
  ALLOCATE(u(clen,rank))
  ALLOCATE(vt(rank,N_t))
  ALLOCATE(s(rank))

  !filling in centering vector with zeros
  DO i=1,clen
    center(i) = 0d0
  END DO

  !calculating centering vector as column sum of data matrix
  p1 = 0
  DO t=1,N_t
    DO i=1,clen
      p1 = p1 + 1
      center(i) = center(i) + dat(p1)
    END DO
  END DO

  !averaging centering vector
  DO i=1,clen
    center(i) = center(i)/DBLE(N_t)
  END DO

  !centering the data matrix
  p1 = 0
  DO t=1,N_t
    DO i=1,clen
      p1 = p1 + 1
      mat(i,t) = dat(p1) - center(i)
    END DO
  END DO

  !calculating the SVD of the datamatrix
  CALL SVD_DECOMP(mat,u,s,vt)

  !vectorizing left singular vector matrix
  p1 = 0
  DO r=1,rank
    DO i=1,clen
      p1 = p1 + 1
      umat(p1) = u(i,r)
    END DO
  END DO

  !vectorizing right singular vector matrix
  p1 = 0
  DO t=1,N_t
    DO r=1,rank
      p1 = p1 + 1
      vtmat(p1) = vt(r,t)
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
!   umat (REAL) -- matrix of dimension(dim1 x MIN(dim1,dim2)), the 'left singular vectors' of datamat are the columns
!   sig (REAL) -- vector of dimension(MIN(dim1,dim2)), the 'singular values' of datamat
!   vtmat (REAL) -- matrix of dimension(MIN(dim1,dim2) x dim2), the 'right singular vectors' of datamat are the rows
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

  CALL DGESVD('S','O',dim1,dim2,datamat,dim1,sig,umat,dim1,vtmat,dim2,WORK,LWORK,info)  !lapack SVD function
  IF (info .NE. 0) STOP "SVD failed"

  vtmat = datamat(1:MIN(dim1,dim2),1:dim2)

  DEALLOCATE(WORK)

END SUBROUTINE SVD_DECOMP

END MODULE POD_Routines
