MODULE CONVERGENCE_CHECKS

  IMPLICIT NONE

CONTAINS

FUNCTION CONVERGENCE(Conv_Type,Eps,x,x_old,x_old2,Its,Norm,Rho)

  IMPLICIT NONE
  REAL*8,INTENT(IN):: Eps
  REAL*8,INTENT(IN):: x(:,:), x_old(:,:), x_old2(:,:)
  INTEGER,INTENT(IN):: Conv_Type, Its

  REAL*8,INTENT(OUT):: Norm
  REAL*8,INTENT(OUT):: Rho

  REAL*8:: New_Eps
  LOGICAL:: CONVERGENCE

  IF (Its .GT. 1) THEN
    Rho = MAXVAL(ABS(x - x_old))/MAXVAL(ABS(x_old - x_old2))
  ELSE
    Rho = 2d0
  END IF

  IF (conv_type .EQ. 1) THEN
    Norm = MAXVAL(ABS(x - x_old))/MAXVAL(ABS(x))

  ELSE IF (conv_type .EQ. 2) THEN
    Norm = MAXVAL(ABS(1d0-x_old/x))

  ELSE IF (conv_type .EQ. 3) THEN
    Norm = MAXVAL(ABS(x - x_old))/MINVAL(ABS(x))

  END IF

  New_Eps = ABS(1d0/Rho-1d0)*Eps
  IF (Norm .GT. New_Eps) THEN
    CONVERGENCE = .FALSE.
  ELSE
    CONVERGENCE = .TRUE.
  END IF

END FUNCTION

END MODULE CONVERGENCE_CHECKS
