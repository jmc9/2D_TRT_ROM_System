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
SUBROUTINE QD_FV(E_avg,E_edgV,E_edgH,EB_L,EB_B,EB_C,EB_R,EB_T,MBx_C,MBx_R,MBx_B,MBx_T,MBy_C,MBy_T,MBy_L,MBy_R,&
  EB_RHS,MBx_RHS,MBy_RHS,Cp_L,Cp_B,Cp_R,Cp_T,BC_L,BC_B,BC_R,BC_T)

  REAL*8,INTENT(OUT):: E_avg(:,:), E_edgV(:,:), E_edgH(:,:)

  REAL*8,INTENT(IN):: EB_L(:,:), EB_B(:,:), EB_C(:,:), EB_R(:,:), EB_T(:,:)
  REAL*8,INTENT(IN):: MBx_C(:,:), MBx_R(:,:), MBx_B(:,:), MBx_T(:,:)
  REAL*8,INTENT(IN):: MBy_C(:,:), MBy_T(:,:), MBy_L(:,:), MBy_R(:,:)
  REAL*8,INTENT(IN):: EB_RHS(:,:), MBx_RHS(:,:), MBy_RHS(:,:)
  REAL*8,INTENT(IN):: Cp_L(:), Cp_B(:), Cp_R(:), Cp_T(:)
  REAL*8,INTENT(IN):: BC_L(:), BC_B(:), BC_R(:), BC_T(:)

  REAL*8,ALLOCATABLE:: mat(:,:), sol(:)
  INTEGER:: N_x, N_y, i, j, r, t, k, row_n, row_size, b, br

  N_y = SIZE(E_avg,2)
  N_x = SIZE(E_avg,1)

  ALLOCATE(mat(N_x+N_y+3*N_x*N_y,N_x+N_y+3*N_x*N_y))
  ALLOCATE(sol(N_x+N_y+3*N_x*N_y))

  mat=0d0
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
  mat(1,1) = Cp_L(1) !left
  mat(1,2) = -MBx_B(1,1) !bottom
  mat(1,3) = MBx_C(1,1) !center
  mat(1,5) = -MBx_T(1,1) !top
  sol(1) = BC_L(1)

  !bottom edge BC
  mat(2,1) = -MBy_L(1,1)
  mat(2,2) = Cp_B(1)
  mat(2,3) = MBy_C(1,1)
  mat(2,4) = -MBy_R(1,1)
  sol(2) = BC_B(1)

  !rad energy balance
  mat(3,1) = EB_L(1,1)
  mat(3,2) = EB_B(1,1)
  mat(3,3) = EB_C(1,1)
  mat(3,4) = EB_R(1,1)
  mat(3,5) = EB_T(1,1)
  sol(3) = EB_RHS(1,1)

  !rad (x) momentum balance
  mat(4,2) = MBx_B(1,1)
  mat(4,3) = MBx_C(1,1)
  mat(4,4) = MBx_R(1,1)
  mat(4,5) = MBx_T(1,1)
  mat(4,6) = -MBx_B(2,1)
  mat(4,7) = MBx_C(2,1)
  mat(4,9) = -MBx_T(2,1)
  sol(4) = MBx_RHS(1,1)

  !rad (y) momentum balance
  t = 4*N_x + 1 !points to the top face of the cell(N_x,1)
  mat(5,1) = MBy_L(1,1)
  mat(5,3) = MBy_C(1,1)
  mat(5,4) = MBy_R(1,1)
  mat(5,5) = MBy_T(1,1)
  mat(5,t+1) = -MBy_L(1,2)
  mat(5,t+2) = MBy_C(1,2)
  mat(5,t+3) = -MBy_R(1,2)
  sol(5) = MBy_RHS(1,1)

  !--------------------------------------------------!
  !                   Bottom Edge                    !
  !--------------------------------------------------!
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
    mat(k-3,k-5) = -MBy_L(i,1)
    mat(k-3,k-3) = Cp_B(i)
    mat(k-3,k-2) = MBy_C(i,1)
    mat(k-3,k-1) = -MBy_R(i,1)
    sol(k-3) = BC_B(i)

    !rad energy balance
    mat(k-2,k-5) = EB_L(i,1)
    mat(k-2,k-3) = EB_B(i,1)
    mat(k-2,k-2) = EB_C(i,1)
    mat(k-2,k-1) = EB_R(i,1)
    mat(k-2,k)   = EB_T(i,1)
    sol(k-2) = EB_RHS(i,1)

    !rad (x) momentum balance
    IF (i .NE. N_x) THEN
      mat(k-1,k-3) = MBx_B(i,1)
      mat(k-1,k-2) = MBx_C(i,1)
      mat(k-1,k-1) = MBx_R(i,1)
      mat(k-1,k)   = MBx_T(i,1)
      mat(k-1,r-3) = -MBx_B(i+1,1)
      mat(k-1,r-2) = MBx_C(i+1,1)
      mat(k-1,r)   = -MBx_T(i+1,1)
      sol(k-1) = MBx_RHS(i,1)
    END IF

    !rad (y) momentum balance
    mat(k,k-5) = MBy_L(i,1)
    mat(k,k-2) = MBy_C(i,1)
    mat(k,k-1) = MBy_R(i,1)
    mat(k,k)   = MBy_T(i,1)
    mat(k,t-4) = -MBy_L(i,2)
    mat(k,t-2) = MBy_C(i,2)
    mat(k,t-1) = -MBy_R(i,2)
    sol(k) = MBy_RHS(i,1)

  END DO

  !--------------------------------------------------!
  !              Bottom right Corner                 !
  !--------------------------------------------------!
  i = N_x
  k = 4*N_x + 1

  !k = top
  !k-1 = right
  !k-2 = center
  !k-3 = bottom
  !k-4 = top(N_x-1)
  !k-5 = right(N_x-1) = left

  !right edge BC
  mat(k-1,k-3) = -MBx_B(N_x,1)
  mat(k-1,k-2) = -MBx_C(N_x,1)
  mat(k-1,k-1) = Cp_R(1)
  mat(k-1,k)   = -MBx_T(N_x,1)

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
    mat(k-3,k-3) = Cp_L(j) !left
    mat(k-3,b)   = -MBx_B(1,j) !bottom
    mat(k-3,k-2) = MBx_C(1,j) !center
    mat(k-3,k)   = -MBx_T(1,j) !top
    sol(k-3) = BC_L(j)

    !rad energy balance
    mat(k-2,k-3) = EB_L(i,j)
    mat(k-2,b)   = EB_B(i,j)
    mat(k-2,k-2) = EB_C(i,j)
    mat(k-2,k-1) = EB_R(i,j)
    mat(k-2,k)   = EB_T(i,j)
    sol(k-2) = EB_RHS(i,j)

    !rad (x) momentum balance
    mat(k-1,b)   = MBx_B(i,j)
    mat(k-1,k-2) = MBx_C(i,j)
    mat(k-1,k-1) = MBx_R(i,j)
    mat(k-1,k)   = MBx_T(i,j)
    mat(k-1,br) = -MBx_B(i+1,j)
    mat(k-1,r-2) = MBx_C(i+1,j)
    mat(k-1,r)   = -MBx_T(i+1,j)
    sol(k-1) = MBx_RHS(i,j)

    !if on the top cell row, must invoke the top edge boundary condition in place of the y-momentum balance
    IF (j .NE. N_y) THEN
      !rad (y) momentum balance
      mat(k,k-3) = MBy_L(i,j)
      mat(k,k-2) = MBy_C(i,j)
      mat(k,k-1) = MBy_R(i,j)
      mat(k,k)   = MBy_T(i,j)
      mat(k,t-3) = -MBy_L(i,j+1)
      mat(k,t-2) = MBy_C(i,j+1)
      mat(k,t-1) = -MBy_R(i,j+1)
      sol(k) = MBy_RHS(i,j)
    ELSE
      !Top edge BC
      mat(k,k-3) = -MBy_L(i,j)
      mat(k,k-2) = -MBy_C(i,j)
      mat(k,k-1) = -MBy_R(i,j)
      mat(k,k)   = Cp_T(i)
      sol(k) = BC_T(i)
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
      mat(k-2,k-4) = EB_L(i,j)
      mat(k-2,b)   = EB_B(i,j)
      mat(k-2,k-2) = EB_C(i,j)
      mat(k-2,k-1) = EB_R(i,j)
      mat(k-2,k)   = EB_T(i,j)
      sol(k-2) = EB_RHS(i,j)

      !if on the right-most cell of each row, must invoke the right edge boundary condition in place of the x-momentum balance
      IF (i .NE. N_x) THEN
        !rad (x) momentum balance
        mat(k-1,b)   = MBx_B(i,j)
        mat(k-1,k-2) = MBx_C(i,j)
        mat(k-1,k-1) = MBx_R(i,j)
        mat(k-1,k)   = MBx_T(i,j)
        mat(k-1,br) = -MBx_B(i+1,j)
        mat(k-1,r-2) = MBx_C(i+1,j)
        mat(k-1,r)   = -MBx_T(i+1,j)
        sol(k-1) = MBx_RHS(i,j)
      ELSE
        !right edge BC
        mat(k-1,b)   = -MBx_B(N_x,j)
        mat(k-1,k-2) = -MBx_C(N_x,j)
        mat(k-1,k-1) = Cp_R(j)
        mat(k-1,k)   = -MBx_T(N_x,j)
        sol(k-1) = BC_R(j)
      END IF

      !if on the top cell row, must invoke the top edge boundary condition in place of the y-momentum balance
      IF (j .NE. N_y) THEN
        !rad (y) momentum balance
        mat(k,k-4) = MBy_L(i,j)
        mat(k,k-2) = MBy_C(i,j)
        mat(k,k-1) = MBy_R(i,j)
        mat(k,k)   = MBy_T(i,j)
        mat(k,t-4) = -MBy_L(i,j+1)
        mat(k,t-2) = MBy_C(i,j+1)
        mat(k,t-1) = -MBy_R(i,j+1)
        sol(k) = MBy_RHS(i,j)
      ELSE
        !Top edge BC
        mat(k,k-4) = -MBy_L(i,j)
        mat(k,k-2) = -MBy_C(i,j)
        mat(k,k-1) = -MBy_R(i,j)
        mat(k,k)   = Cp_T(i)
        sol(k) = BC_T(i)
      END IF

    END DO

    row_n = row_n + row_size !moving row_n to the right boundary cell edge of the current row
  END DO !End loop over N_y

  !===========================================================================!
  !                                                                           !
  !     Solving the now built linear system                                   !
  !                                                                           !
  !     Just using a direct solve right now to test                           !
  !                                                                           !
  !===========================================================================!
  CALL LP_MATSOLVE(mat,sol) !GAUSS_COMPLETE(mat,sol)

  E_edgV(1,1) = sol(1) !left
  DO i=1,N_x
    k = 4*i + 1
    E_edgH(i,1)   = sol(k-3) !bottom
    E_avg(i,1)    = sol(k-2) !center
    E_edgV(i+1,1) = sol(k-1) !right
    E_edgH(i,2)   = sol(k)   !top
  END DO

  row_n = 4*N_x + 1
  row_size = 3*N_x + 1

  DO j=2,N_y
    k = row_n + 4
    E_edgV(1,j) = sol(k-3) !left
    DO i=1,N_x
      k = row_n + 3*i + 1
      E_avg(i,j)    = sol(k-2) !center
      E_edgV(i+1,j) = sol(k-1) !right
      E_edgH(i,j+1) = sol(k)   !top
    END DO
    row_n = row_n + row_size
  END DO

  DEALLOCATE(mat,sol)

END SUBROUTINE QD_FV

!==================================================================================================================================!
!
!==================================================================================================================================!

END MODULE QD_SOLVERS
