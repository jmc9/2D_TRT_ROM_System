MODULE QD_SOLVERS

  USE LA_TOOLS

  IMPLICIT NONE

CONTAINS

!==================================================================================================================================!
!Subroutine QD_FV
!
! Solves the 2D radiative transfer quasidiffusion equations discretized with a finite volumes scheme
!
! NOTES::
!   This subroutine solves the one-group QD equations, multigroup problems need only wrap this in a loop over groups
!   This subroutine ONLY SOLVES FOR ENERGY DENSITIES using a reduced linear system (obtained from algebraic manipulation)
!   As such, attempting to use coefficients from the non-reduced QD system will not yield the correct solution
!   If the user would like to find radiation fluxes, calculate those using the solution from this subroutine
!   The 'reduced linear system' mentioned above is derived in the documentation (my notes)
!
! WARNINGS::
!
! OUTPUTS::
!   E_avg - cell-averaged radiation energy densities
!   E_edgV - cell-edge radiation energy densities on 'verticle' cell edges (x = constant)
!   E_edgH - cell-edge radiation energy densities on 'horizontal' cell edges (y = constant)
!
! INPUTS::
!   (Note that L,B,C,R,T corresponds to the cell left-edge, bottom-edge, center, right-edge, top-edge, respectively)
!   EB_[L,B,C,R,T] - coefficients attached to radiation energy densities present in the radiation energy balance equation
!   MBx_[C,R,B,T] - coefficients attached to radiation energy densities present in the radiation x-momentum balance equation
!   MBy_[C,R,B,T] - coefficients attached to radiation energy densities present in the radiation y-momentum balance equation
!   EB_RHS - right hand side of the radiation energy balance equation
!   MBx_RHS - right hand side of the radiation x-momentum balance equation
!   MBy_RHS - right hand side of the radiation y-momentum balance equation
!   Cp_[L,B,R,T] - coefficients attached to boundary cell edge radiation energy densities present in the boundary conditions
!   BC_[L,B,R,T] - right hand sides of the boundary conditions
!
!==================================================================================================================================!
RECURSIVE SUBROUTINE QD_FV(E_avg,E_edgV,E_edgH,Its,EB_L,EB_B,EB_C,EB_R,EB_T,MBx_C,MBx_R,MBx_B,MBx_T,MBy_C,MBy_T,MBy_L,MBy_R,&
  EB_RHS,MBx_RHS,MBy_RHS,Cp_L,Cp_B,Cp_R,Cp_T,BC_L,BC_B,BC_R,BC_T)

  REAL*8,INTENT(INOUT):: E_avg(:,:), E_edgV(:,:), E_edgH(:,:)
  INTEGER,INTENT(OUT):: Its

  REAL*8,INTENT(IN):: EB_L(:,:), EB_B(:,:), EB_C(:,:), EB_R(:,:), EB_T(:,:)
  REAL*8,INTENT(IN):: MBx_C(:,:), MBx_R(:,:), MBx_B(:,:), MBx_T(:,:)
  REAL*8,INTENT(IN):: MBy_C(:,:), MBy_T(:,:), MBy_L(:,:), MBy_R(:,:)
  REAL*8,INTENT(IN):: EB_RHS(:,:), MBx_RHS(:,:), MBy_RHS(:,:)
  REAL*8,INTENT(IN):: Cp_L(:), Cp_B(:), Cp_R(:), Cp_T(:)
  REAL*8,INTENT(IN):: BC_L(:), BC_B(:), BC_R(:), BC_T(:)

  REAL*8,ALLOCATABLE:: mat(:,:), sol(:), sol2(:)
  REAL*8,ALLOCATABLE:: mat_sp(:), mat_Usp(:), w(:), fpar(:), w2(:)
  INTEGER,ALLOCATABLE:: ja(:), ia(:), jlu(:), ju(:), jw(:), ipar(:)
  INTEGER:: N_x, N_y, i, j, r, t, k, row_n, row_size, b, br, err, pj
  INTEGER:: Nrows, Nelm

  N_y = SIZE(E_avg,2)
  N_x = SIZE(E_avg,1)

  Nrows = N_x+N_y+3*N_x*N_y
  Nelm = 19*N_x*N_y+N_x+N_y

  ALLOCATE(mat_sp(Nelm), ja(Nelm), ia(Nrows+1))
  ALLOCATE(sol(Nrows))

  pj = 1

  !===========================================================================!
  !                                                                           !
  !     Building the linear system (part 1)                                   !
  !     First set of 'Blocks' for the bottom row along the bottom boundary    !
  !                                                                           !
  !     First unknown is the left boundary face                               !
  !     Cell unknowns ordered as: bottom, center, right, top                  !
  !     First cell unknowns: left, bottom, center, right, top                 !
  !                                                                           !
  !     Total unknowns in this row: 4*N_x + 1                                 !
  !                                                                           !
  !===========================================================================!
  !--------------------------------------------------!
  !               Bottom Left Corner                 !
  !--------------------------------------------------!
  !1 = left
  !2 = bottom
  !3 = center
  !4 = right
  !5 = top

  !left edge BC
  ia(1) = pj
  sol(1) = BC_L(1)
  CALL SPARSE_IN(pj, mat_sp, Cp_L(1),     ja, 1)
  CALL SPARSE_IN(pj, mat_sp, -MBx_B(1,1), ja, 2)
  CALL SPARSE_IN(pj, mat_sp, MBx_C(1,1),  ja, 3)
  CALL SPARSE_IN(pj, mat_sp, -MBx_T(1,1), ja, 5)

  !bottom edge BC
  ia(2) = pj
  sol(2) = BC_B(1)
  CALL SPARSE_IN(pj, mat_sp, -MBy_L(1,1), ja, 1)
  CALL SPARSE_IN(pj, mat_sp, Cp_B(1),     ja, 2)
  CALL SPARSE_IN(pj, mat_sp, MBy_C(1,1),  ja, 3)
  CALL SPARSE_IN(pj, mat_sp, -MBy_R(1,1), ja, 4)

  !rad energy balance
  ia(3) = pj
  sol(3) = EB_RHS(1,1)
  CALL SPARSE_IN(pj, mat_sp, EB_L(1,1), ja, 1)
  CALL SPARSE_IN(pj, mat_sp, EB_B(1,1), ja, 2)
  CALL SPARSE_IN(pj, mat_sp, EB_C(1,1), ja, 3)
  CALL SPARSE_IN(pj, mat_sp, EB_R(1,1), ja, 4)
  CALL SPARSE_IN(pj, mat_sp, EB_T(1,1), ja, 5)

  !rad (x) momentum balance
  ia(4) = pj
  sol(4) = MBx_RHS(1,1)
  CALL SPARSE_IN(pj, mat_sp, MBx_B(1,1),  ja, 2)
  CALL SPARSE_IN(pj, mat_sp, MBx_C(1,1),  ja, 3)
  CALL SPARSE_IN(pj, mat_sp, MBx_R(1,1),  ja, 4)
  CALL SPARSE_IN(pj, mat_sp, MBx_T(1,1),  ja, 5)
  CALL SPARSE_IN(pj, mat_sp, -MBx_B(2,1), ja, 6)
  CALL SPARSE_IN(pj, mat_sp, MBx_C(2,1),  ja, 7)
  CALL SPARSE_IN(pj, mat_sp, -MBx_T(2,1), ja, 9)

  !rad (y) momentum balance
  t = 4*N_x + 1 !points to the top face of the cell(N_x,1)
  ia(5) = pj
  sol(5) = MBy_RHS(1,1)
  CALL SPARSE_IN(pj, mat_sp, MBy_L(1,1),  ja, 1)
  CALL SPARSE_IN(pj, mat_sp, MBy_C(1,1),  ja, 3)
  CALL SPARSE_IN(pj, mat_sp, MBy_R(1,1),  ja, 4)
  CALL SPARSE_IN(pj, mat_sp, MBy_T(1,1),  ja, 5)
  CALL SPARSE_IN(pj, mat_sp, -MBy_L(1,2), ja, t+1)
  CALL SPARSE_IN(pj, mat_sp, MBy_C(1,2),  ja, t+2)
  CALL SPARSE_IN(pj, mat_sp, -MBy_R(1,2), ja, t+3)

  !--------------------------------------------------!
  !                   Bottom Edge                    !
  !--------------------------------------------------!
  ! pj = 28
  DO i=2,N_x
    k = 4*i + 1 !points to top edge of cell(i,1)
    r = 4*(i+1) + 1 !points to top edge of cell(i+1,1)
    t = 4*N_x + 1 + 3*i + 1 !points to top edge of cell(i,2)

    !k = top
    !k-1 = right
    !k-2 = center
    !k-3 = bottom
    !k-4 = top(i-1)
    !k-5 = right(i-1) = left(i)

    !r = top(i+1)
    !r-1 = right(i+1)
    !r-2 = center(i+1)
    !r-3 = bottom(i+1)
    !r-4 = top(i) = k

    !t = top(j+1)
    !t-1 = right(j+1)
    !t-2 = center(j+1)
    !t-3 = top(i-1,j+1)
    !t-4 = right(i-1,j+1) = left(j+1)

    !bottom edge BC
    ia(k-3) = pj
    sol(k-3) = BC_B(i)
    CALL SPARSE_IN(pj, mat_sp, -MBy_L(i,1), ja, k-5)
    CALL SPARSE_IN(pj, mat_sp, Cp_B(i),     ja, k-3)
    CALL SPARSE_IN(pj, mat_sp, MBy_C(i,1),  ja, k-2)
    CALL SPARSE_IN(pj, mat_sp, -MBy_R(i,1), ja, k-1)

    !rad energy balance
    ia(k-2) = pj
    sol(k-2) = EB_RHS(i,1)
    CALL SPARSE_IN(pj, mat_sp, EB_L(i,1), ja, k-5)
    CALL SPARSE_IN(pj, mat_sp, EB_B(i,1), ja, k-3)
    CALL SPARSE_IN(pj, mat_sp, EB_C(i,1), ja, k-2)
    CALL SPARSE_IN(pj, mat_sp, EB_R(i,1), ja, k-1)
    CALL SPARSE_IN(pj, mat_sp, EB_T(i,1), ja, k)

    !rad (x) momentum balance
    IF (i .NE. N_x) THEN
      ia(k-1) = pj
      sol(k-1) = MBx_RHS(i,1)
      CALL SPARSE_IN(pj, mat_sp, MBx_B(i,1),    ja, k-3)
      CALL SPARSE_IN(pj, mat_sp, MBx_C(i,1),    ja, k-2)
      CALL SPARSE_IN(pj, mat_sp, MBx_R(i,1),    ja, k-1)
      CALL SPARSE_IN(pj, mat_sp, MBx_T(i,1),    ja, k)
      CALL SPARSE_IN(pj, mat_sp, -MBx_B(i+1,1), ja, r-3)
      CALL SPARSE_IN(pj, mat_sp, MBx_C(i+1,1),  ja, r-2)
      CALL SPARSE_IN(pj, mat_sp, -MBx_T(i+1,1), ja, r)

    ELSE
      !Bottom right corner right face boundary condition
      ia(k-1) = pj
      sol(k-1) = BC_R(1)
      CALL SPARSE_IN(pj, mat_sp, -MBx_B(N_x,1), ja, k-3)
      CALL SPARSE_IN(pj, mat_sp, -MBx_C(N_x,1), ja, k-2)
      CALL SPARSE_IN(pj, mat_sp, Cp_R(1),       ja, k-1)
      CALL SPARSE_IN(pj, mat_sp, -MBx_T(N_x,1), ja, k)

    END IF

    !rad (y) momentum balance
    ia(k) = pj
    sol(k) = MBy_RHS(i,1)
    CALL SPARSE_IN(pj, mat_sp, MBy_L(i,1),  ja, k-5)
    CALL SPARSE_IN(pj, mat_sp, MBy_C(i,1),  ja, k-2)
    CALL SPARSE_IN(pj, mat_sp, MBy_R(i,1),  ja, k-1)
    CALL SPARSE_IN(pj, mat_sp, MBy_T(i,1),  ja, k)
    CALL SPARSE_IN(pj, mat_sp, -MBy_L(i,2), ja, t-4)
    CALL SPARSE_IN(pj, mat_sp, MBy_C(i,2),  ja, t-2)
    CALL SPARSE_IN(pj, mat_sp, -MBy_R(i,2), ja, t-1)


  END DO


  !===========================================================================!
  !                                                                           !
  !     Building the linear system (part 2)                                   !
  !     Second set of 'Blocks' for the center rows with no top/bottom -       !
  !     - boundary faces                                                      !
  !                                                                           !
  !     First unknown is the left boundary face                               !
  !     Cell unknowns ordered as: center, right, top                          !
  !     First cell unknowns: left, center, right, top                         !
  !                                                                           !
  !     Total unknowns in each row: 3*N_x + 1                                 !
  !                                                                           !
  !===========================================================================!
  row_n = 4*N_x + 1 !row_n points to the right boundary cell edge of the last (j-1) row
  row_size = 3*N_x + 1 !row_size contains the number of unknowns contained in the current (j) row
  DO j=2,N_y
    !--------------------------------------------------!
    !                  Left Boundary                   !
    !--------------------------------------------------!
    i = 1
    k = row_n + 4 !points to top edge of cell(i,j)
    IF (j .EQ. 2) THEN !j=2 has special pointer to bottom face since j=1 has an extra unknown per cell
      b = 5 !points to bottom edge of cell(i,j)
      br = b+4 !points to bottom edge of cell(i+1,j)
    ELSE
      b = k - row_size !points to bottom edge of cell(i,j)
      br = b+3 !points to bottom edge of cell(i+1,j)
    END IF
    r = k + 3 !points to top edge of cell(i+1,j)
    t = k + row_size !points to top edge of cell(i,j+1)

    !k = top
    !k-1 = right
    !k-2 = center
    !k-3 = left

    !b = top(j-1) = bottom
    !b+4 = top(i+1,j-1) = bottom(i+1)

    !left edge BC
    ia(k-3) = pj
    sol(k-3) = BC_L(j)
    CALL SPARSE_IN(pj, mat_sp, Cp_L(j),     ja, k-3)
    CALL SPARSE_IN(pj, mat_sp, -MBx_B(1,j), ja, b)
    CALL SPARSE_IN(pj, mat_sp, MBx_C(1,j),  ja, k-2)
    CALL SPARSE_IN(pj, mat_sp, -MBx_T(1,j), ja, k)

    !rad energy balance
    ia(k-2) = pj
    sol(k-2) = EB_RHS(i,j)
    CALL SPARSE_IN(pj, mat_sp, EB_L(i,j), ja, k-3)
    CALL SPARSE_IN(pj, mat_sp, EB_B(i,j), ja, b)
    CALL SPARSE_IN(pj, mat_sp, EB_C(i,j), ja, k-2)
    CALL SPARSE_IN(pj, mat_sp, EB_R(i,j), ja, k-1)
    CALL SPARSE_IN(pj, mat_sp, EB_T(i,j), ja, k)

    !rad (x) momentum balance
    ia(k-1) = pj
    sol(k-1) = MBx_RHS(i,j)
    CALL SPARSE_IN(pj, mat_sp, MBx_B(i,j),    ja, b)
    CALL SPARSE_IN(pj, mat_sp, MBx_C(i,j),    ja, k-2)
    CALL SPARSE_IN(pj, mat_sp, MBx_R(i,j),    ja, k-1)
    CALL SPARSE_IN(pj, mat_sp, MBx_T(i,j),    ja, k)
    CALL SPARSE_IN(pj, mat_sp, -MBx_B(i+1,j), ja, br)
    CALL SPARSE_IN(pj, mat_sp, MBx_C(i+1,j),  ja, r-2)
    CALL SPARSE_IN(pj, mat_sp, -MBx_T(i+1,j), ja, r)

    !if on the top cell row, must invoke the top edge boundary condition in place of the y-momentum balance
    IF (j .NE. N_y) THEN
      !rad (y) momentum balance
      ia(k) = pj
      sol(k) = MBy_RHS(i,j)
      CALL SPARSE_IN(pj, mat_sp, MBy_L(i,j),    ja, k-3)
      CALL SPARSE_IN(pj, mat_sp, MBy_C(i,j),    ja, k-2)
      CALL SPARSE_IN(pj, mat_sp, MBy_R(i,j),    ja, k-1)
      CALL SPARSE_IN(pj, mat_sp, MBy_T(i,j),    ja, k)
      CALL SPARSE_IN(pj, mat_sp, -MBy_L(i,j+1), ja, t-3)
      CALL SPARSE_IN(pj, mat_sp, MBy_C(i,j+1),  ja, t-2)
      CALL SPARSE_IN(pj, mat_sp, -MBy_R(i,j+1), ja, t-1)

    ELSE
      !Top edge BC
      ia(k) = pj
      sol(k) = BC_T(i)
      CALL SPARSE_IN(pj, mat_sp, -MBy_L(i,j), ja, k-3)
      CALL SPARSE_IN(pj, mat_sp, -MBy_C(i,j), ja, k-2)
      CALL SPARSE_IN(pj, mat_sp, -MBy_R(i,j), ja, k-1)
      CALL SPARSE_IN(pj, mat_sp, Cp_T(i),     ja, k)

    END IF

    !--------------------------------------------------!
    !                      Center                      !
    !--------------------------------------------------!
    DO i=2,N_x

      k = row_n + 3*i + 1
      IF (j .EQ. 2) THEN
        b = 4*i + 1
        br = b+4
      ELSE
        b = k - row_size
        br = b+3
      END IF
      r = k + 3
      t = k + row_size

      !k = top
      !k-1 = right
      !k-2 = center
      !k-3 = top(i-1)
      !k-4 = right(i-1) = left
      !above rules apply to r for (i+1) and t for (j+1)

      !b = top(j-1) = bottom
      !b+4 = top(i+1,j-1) = bottom(i+1) --> ONLY FOR j=2
      !b+3 = top(i+1,j-1) = bottom(i+1) --> FOR ALL j>2

      !rad energy balance
      ia(k-2) = pj
      sol(k-2) = EB_RHS(i,j)
      CALL SPARSE_IN(pj, mat_sp, EB_L(i,j), ja, k-4)
      CALL SPARSE_IN(pj, mat_sp, EB_B(i,j), ja, b)
      CALL SPARSE_IN(pj, mat_sp, EB_C(i,j), ja, k-2)
      CALL SPARSE_IN(pj, mat_sp, EB_R(i,j), ja, k-1)
      CALL SPARSE_IN(pj, mat_sp, EB_T(i,j), ja, k)

      !if on the right-most cell of each row, must invoke the right edge boundary condition in place of the x-momentum balance
      IF (i .NE. N_x) THEN
        !rad (x) momentum balance
        ia(k-1) = pj
        sol(k-1) = MBx_RHS(i,j)
        CALL SPARSE_IN(pj, mat_sp, MBx_B(i,j),    ja, b)
        CALL SPARSE_IN(pj, mat_sp, MBx_C(i,j),    ja, k-2)
        CALL SPARSE_IN(pj, mat_sp, MBx_R(i,j),    ja, k-1)
        CALL SPARSE_IN(pj, mat_sp, MBx_T(i,j),    ja, k)
        CALL SPARSE_IN(pj, mat_sp, -MBx_B(i+1,j), ja, br)
        CALL SPARSE_IN(pj, mat_sp, MBx_C(i+1,j),  ja, r-2)
        CALL SPARSE_IN(pj, mat_sp, -MBx_T(i+1,j), ja, r)

      ELSE
        !right edge BC
        ia(k-1) = pj
        sol(k-1) = BC_R(j)
        CALL SPARSE_IN(pj, mat_sp, -MBx_B(N_x,j), ja, b)
        CALL SPARSE_IN(pj, mat_sp, -MBx_C(N_x,j), ja, k-2)
        CALL SPARSE_IN(pj, mat_sp, Cp_R(j),       ja, k-1)
        CALL SPARSE_IN(pj, mat_sp, -MBx_T(N_x,j), ja, k)

      END IF

      !if on the top cell row, must invoke the top edge boundary condition in place of the y-momentum balance
      IF (j .NE. N_y) THEN
        !rad (y) momentum balance
        ia(k) = pj
        sol(k) = MBy_RHS(i,j)
        CALL SPARSE_IN(pj, mat_sp, MBy_L(i,j),    ja, k-4)
        CALL SPARSE_IN(pj, mat_sp, MBy_C(i,j),    ja, k-2)
        CALL SPARSE_IN(pj, mat_sp, MBy_R(i,j),    ja, k-1)
        CALL SPARSE_IN(pj, mat_sp, MBy_T(i,j),    ja, k)
        CALL SPARSE_IN(pj, mat_sp, -MBy_L(i,j+1), ja, t-4)
        CALL SPARSE_IN(pj, mat_sp, MBy_C(i,j+1),  ja, t-2)
        CALL SPARSE_IN(pj, mat_sp, -MBy_R(i,j+1), ja, t-1)

      ELSE
        !Top edge BC
        ia(k) = pj
        sol(k) = BC_T(i)
        CALL SPARSE_IN(pj, mat_sp, -MBy_L(i,j), ja, k-4)
        CALL SPARSE_IN(pj, mat_sp, -MBy_C(i,j), ja, k-2)
        CALL SPARSE_IN(pj, mat_sp, -MBy_R(i,j), ja, k-1)
        CALL SPARSE_IN(pj, mat_sp, Cp_T(i),     ja, k)

      END IF

    END DO

    row_n = row_n + row_size !moving row_n to the right boundary cell edge of the current row
  END DO !End loop over N_y

  ia(Nrows+1) = pj

  !===========================================================================!
  !                                                                           !
  !     Solving the now built linear system                                   !
  !                                                                           !
  !     Linear solver: BI-CGSTAB preconditioned with ILUT                     !
  !     -> Using tools from the SPARSKIT2 library                             !
  !                                                                           !
  !===========================================================================!
  !allocating all arrays for the linear solve
  ALLOCATE(mat_Usp(Nelm), jlu(Nelm), ju(Nrows))
  ALLOCATE(w(Nrows+1), jw(2*(Nrows)))
  ALLOCATE(sol2(Nrows), w2(Nrows*8), ipar(16), fpar(16))

  !initializing all new arrays to zero
  mat_Usp = 0d0
  jlu = 0
  ju = 0
  w = 0d0
  jw = 0
  ipar = 0
  fpar = 0d0
  w2 = 0d0
  sol2 = 0d0

  !Performing preconditioning with ILUT
  !using a fill level of 3 and 'drop tolerance' of 1d-6
  CALL ILUT(Nrows ,mat_sp, ja, ia, 3, 1d-6, mat_Usp, jlu, ju, Nelm, w, jw, err)
  IF (err .NE. 0) THEN
    WRITE(*,*) err
    STOP 'BAD ILUT'
  END IF

  !setting options for BI-CGSTAB
  ipar(1) = 0  !tell solver to initialize new linear solve
  ipar(2) = 2  !use right preconditioning
  ipar(3) = -1 !termination condition
  ipar(4) = Nrows*8 !length of the work space
  ipar(5) = 10   !useless for BI-CGSTAB, default is 10
  ipar(6) = 2000 !max number of mat-vec products

  fpar(1) = 1d-13 !relative tolerance
  fpar(2) = 1d-15 !absolute tolerance

  !putting initial guess for solution in the vector sol2
  sol2(1) = E_edgV(1,1)
  DO i=1,N_x
    k = 4*i + 1
    sol2(k-3) = E_edgH(i,1)
    sol2(k-2) = E_avg(i,1)
    sol2(k-1) = E_edgV(i+1,1)
    sol2(k) = E_edgH(i,2)
  END DO

  row_n = 4*N_x + 1
  row_size = 3*N_x + 1

  DO j=2,N_y
    k = row_n + 4
    sol2(k-3) = E_edgV(1,j)
    DO i=1,N_x
      k = row_n + 3*i + 1
      sol2(k-2) = E_avg(i,j)
      sol2(k-1) = E_edgV(i+1,j)
      sol2(k) = E_edgH(i,j+1)
    END DO
    row_n = row_n + row_size
  END DO

  !--------------------------------------------------!
  !      Performing linear solve via BI-CGSTAB       !
  !--------------------------------------------------!
  KLOOP: DO
    CALL BCGSTAB(Nrows, sol, sol2, ipar, fpar, w2)

    IF (ipar(1) .EQ. 0) THEN !successful termination
      EXIT KLOOP
    ELSE IF (ipar(1) .EQ. 1) THEN !require matrix vector product with A
      CALL AMUX(Nrows, w2(ipar(8)), w2(ipar(9)), mat_sp, ja, ia)
      CYCLE KLOOP
    ELSE IF (ipar(1) .EQ. 2) THEN !require matrix vector product with A^T
      CALL ATMUX(Nrows, w2(ipar(8)), w2(ipar(9)), mat_sp, ja, ia)
      CYCLE KLOOP
    ELSE IF ((ipar(1) .EQ. 3).OR.(ipar(1) .EQ. 5)) THEN !require matrix vector product with M
      CALL LUSOL(Nrows, w2(ipar(8)), w2(ipar(9)), mat_Usp, jlu, ju)
      CYCLE KLOOP
    ELSE IF ((ipar(1) .EQ. 4).OR.(ipar(1) .EQ. 6)) THEN !require matrix vector product with M^T
      CALL LUTSOL(Nrows, w2(ipar(8)), w2(ipar(9)), mat_Usp, jlu, ju)
      CYCLE KLOOP
    ELSE IF (ipar(1) .LT. 0) THEN !UNsuccessful termination
      ALLOCATE(mat(Nrows, Nrows))
      CALL CSRDNS(Nrows, Nrows, mat_sp, ja, ia, mat, Nrows, err)
      IF (err .NE. 0) STOP 'CSRDNS TERMINATED TOO EARLY'
      CALL GAUSS_COMPLETE(mat, sol) !if BI-CGSTAB fails to solve the system, resort to direct solution by gaussian elimination
      DEALLOCATE(mat)
      sol2 = sol
      ipar(7) = -1 !flag that a direct solve was used
      EXIT KLOOP
    END IF

  END DO KLOOP
  Its = ipar(7) !store total matrix-vector product count

  !deallocate all arrays except for solution vector
  DEALLOCATE(mat_sp, mat_Usp, sol, ja, ia, jlu, ju, w, jw, w2, ipar, fpar)

  !--------------------------------------------------!
  !               Collecting Solution                !
  !--------------------------------------------------!
  E_edgV(1,1) = sol2(1) !left
  DO i=1,N_x
    k = 4*i + 1
    E_edgH(i,1)   = sol2(k-3) !bottom
    E_avg(i,1)    = sol2(k-2) !center
    E_edgV(i+1,1) = sol2(k-1) !right
    E_edgH(i,2)   = sol2(k)   !top
  END DO

  row_n = 4*N_x + 1
  row_size = 3*N_x + 1

  DO j=2,N_y
    k = row_n + 4
    E_edgV(1,j) = sol2(k-3) !left
    DO i=1,N_x
      k = row_n + 3*i + 1
      E_avg(i,j)    = sol2(k-2) !center
      E_edgV(i+1,j) = sol2(k-1) !right
      E_edgH(i,j+1) = sol2(k)   !top
    END DO
    row_n = row_n + row_size
  END DO

  DEALLOCATE(sol2)

END SUBROUTINE QD_FV

!==================================================================================================================================!
!
!==================================================================================================================================!
SUBROUTINE SPARSE_IN(p, mat_sp, val, ja, col)
  INTEGER,INTENT(INOUT):: p, ja(*)
  REAL*8,INTENT(INOUT):: mat_sp(*)
  INTEGER,INTENT(IN):: col
  REAL*8,INTENT(IN):: val

  IF (val .NE. 0d0) THEN
    mat_sp(p) = val
    ja(p) = col
    p = p + 1
  END IF

END SUBROUTINE

!==================================================================================================================================!
!
!==================================================================================================================================!


END MODULE QD_SOLVERS
