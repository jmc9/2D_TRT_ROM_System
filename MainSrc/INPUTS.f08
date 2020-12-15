MODULE INPUTS

  IMPLICIT NONE

CONTAINS

!==================================================================================================================================!
!
!==================================================================================================================================!
SUBROUTINE INPUT(run_type,restart_infile,use_grey,Test,Mat,Kappa_Mult,chi,conv_ho,conv_lo,conv_gr1,&
  conv_gr2,comp_unit,line_src,E_Bound_Low,T_Bound_Low,Theta,maxit_RTE,maxit_MLOQD,maxit_GLOQD,conv_type,N_m,&
  threads,kapE_dT_flag,enrgy_strc,erg,xlen,ylen,N_x,N_y,tlen,delt,bcT_left,bcT_right,bcT_upper,bcT_lower,&
  Tini,sig_R,ar,pi,c,h,delx,dely,cv,outfile,out_freq,I_out,HO_Eg_out,HO_Fg_out,HO_E_out,HO_F_out,Eg_out,Fg_out,MGQD_E_out,&
  MGQD_F_out,QDfg_out,E_out,F_out,D_out,old_parms_out,its_out,conv_out,kap_out,Src_out,nu_g,N_g,Omega_x,Omega_y,&
  quad_weight,N_t,quadrature,BC_Type,Use_Line_Search,Use_Safety_Search,Res_Calc,POD_err,PODgsum,POD_Type,POD_dset,&
  Direc_Diff,xpts_avg,xpts_edgV,ypts_avg,ypts_edgH,tpts)

  IMPLICIT NONE

  !INPUT VARIABLES
  REAL*8,INTENT(IN):: erg, pi, c, h

  !OUTPUT VARIABLES
  REAL*8,INTENT(INOUT):: comp_unit

  INTEGER,INTENT(OUT):: use_grey
  CHARACTER(100),INTENT(OUT):: run_type, restart_infile, Test
  INTEGER,ALLOCATABLE,INTENT(OUT):: Mat(:,:)
  REAL*8,ALLOCATABLE,INTENT(OUT):: Kappa_Mult(:)

  REAL*8,INTENT(OUT):: chi, conv_ho, conv_lo, conv_gr1, conv_gr2, line_src, E_Bound_Low, T_Bound_Low, Theta
  INTEGER,INTENT(OUT):: maxit_RTE, maxit_MLOQD, maxit_GLOQD, conv_type, threads
  CHARACTER(100),INTENT(OUT):: enrgy_strc, quadrature
  LOGICAL,INTENT(OUT):: Use_Line_Search, Use_Safety_Search, Res_Calc, kapE_dT_flag

  REAL*8,INTENT(OUT):: POD_err
  INTEGER,INTENT(OUT):: PODgsum, Direc_Diff
  CHARACTER(100),INTENT(OUT):: POD_Type, POD_dset

  REAL*8,INTENT(OUT):: xlen, ylen, tlen, delt, bcT_left, bcT_right, bcT_upper, bcT_lower, Tini
  REAL*8,ALLOCATABLE,INTENT(OUT):: Delx(:), Dely(:), xpts_avg(:), xpts_edgV(:), ypts_avg(:), ypts_edgH(:), tpts(:)
  INTEGER,INTENT(OUT):: N_x, N_y, N_t, BC_Type(:)

  CHARACTER(100),INTENT(OUT):: outfile
  INTEGER,INTENT(OUT):: out_freq, I_out, HO_Eg_out, HO_Fg_out, HO_E_out, HO_F_out
  INTEGER,INTENT(OUT):: Eg_out, Fg_out, MGQD_E_out, MGQD_F_out, QDfg_out
  INTEGER,INTENT(OUT):: E_out, F_out, D_out
  INTEGER,INTENT(OUT):: old_parms_out, its_out, conv_out, kap_out, Src_out

  REAL*8,ALLOCATABLE,INTENT(OUT):: nu_g(:)
  INTEGER,INTENT(OUT):: N_g

  INTEGER,INTENT(INOUT):: N_m
  REAL*8,ALLOCATABLE,INTENT(OUT):: Omega_x(:), Omega_y(:), quad_weight(:)

  REAL*8,INTENT(OUT):: sig_R, ar, cv

  !LOCAL VARIABLES
  INTEGER:: inpunit = 10
  INTEGER:: err, i, j

  OPEN(UNIT=inpunit,FILE="input/input.inp",STATUS='OLD',ACTION='READ',IOSTAT=err)
  !   making sure file exists/opens, if not tells user
  IF (err .NE. 0) THEN
      WRITE(*,'(A)') 'The file, input.inp, could not open properly.'
      STOP
  END IF

  CALL INPUT_RUN_STATE(inpunit,run_type,restart_infile,use_grey,Res_Calc,Test)

  CALL INPUT_SOLVER_OPTS(inpunit,chi,conv_ho,conv_lo,conv_gr1,conv_gr2,comp_unit,line_src,maxit_RTE,maxit_MLOQD,&
    maxit_GLOQD,conv_type,threads,kapE_dT_flag,enrgy_strc,quadrature,E_Bound_Low,T_Bound_Low,Use_Line_Search,&
    Use_Safety_Search,Theta)

  sig_R=2d0*pi**5/(15d0*c**2*h**3*erg**4*comp_unit) !(erg/(ev**4 cm**2 sh))
  aR=4d0*sig_R/c

  IF (run_type .EQ. 'mg_pod') THEN
    CALL INPUT_POD_OPTS(inpunit,POD_err,PODgsum,POD_Type,POD_dset,Direc_Diff)

    IF (POD_Type .EQ. 'fg') THEN
      maxit_RTE = 1
    END IF

  ELSE IF ((run_type .EQ. 'diff').OR.(run_type .EQ. 'p1').OR.(run_type .EQ. 'p1/3')) THEN
    maxit_RTE = 1
  END IF

  CALL INPUT_PARAMETERS(inpunit,erg,xlen,ylen,N_x,N_y,tlen,delt,BC_Type,bcT_left,bcT_right,bcT_upper,&
    bcT_lower,Tini)

  ALLOCATE(Delx(N_x))
  ALLOCATE(Dely(N_y))
  Delx = xlen/REAL(N_x,8)
  Dely = ylen/REAL(N_y,8)
  N_t = NINT(tlen/delt)
  Cv=0.5917d0*ar*(bcT_left)**3*((1.38d-16)**4*(11600d0**4)*erg**4)

  ALLOCATE(xpts_avg(N_x), xpts_edgV(N_x+1))
  ALLOCATE(ypts_avg(N_y), ypts_edgH(N_y+1))
  ALLOCATE(tpts(N_t))
  xpts_edgV(1) = 0d0
  DO i=1,N_x
    xpts_edgv(i+1) = xpts_edgv(i) + Delx(i)
  END DO
  xpts_avg(1) = Delx(1)/2d0
  DO i=2,N_x
    xpts_avg(i) = xpts_avg(i-1) + (Delx(i) + Delx(i-1))/2d0
  END DO
  ypts_edgH(1) = 0d0
  DO i=1,N_y
    ypts_edgH(i+1) = ypts_edgH(i) + Dely(i)
  END DO
  ypts_avg(1) = Dely(1)/2d0
  DO i=2,N_y
    ypts_avg(i) = ypts_avg(i-1) + (Dely(i) + Dely(i-1))/2d0
  END DO
  tpts(1) = Delt
  DO i=2,N_t
    tpts(i) = tpts(i-1) + Delt
  END DO

  CALL INPUT_TEST_TYPE(Test,N_x,N_y,Delx,Dely,xlen,ylen,Mat,Kappa_Mult)

  CALL INPUT_ENERGY(inpunit,enrgy_strc,N_g,nu_g)

  CALL INPUT_QUAD(quadrature,N_m,Omega_x,Omega_y,quad_weight)

  CALL INPUT_OUTPUT_OPTS(inpunit,outfile,out_freq,I_out,HO_Eg_out,HO_Fg_out,HO_E_out,HO_F_out,Eg_out,Fg_out,MGQD_E_out,&
    MGQD_F_out,QDfg_out,E_out,F_out,D_out,old_parms_out,its_out,conv_out,kap_out,Src_out)

  CLOSE ( inpunit, STATUS='KEEP')

END SUBROUTINE INPUT

!==================================================================================================================================!
!
!==================================================================================================================================!
SUBROUTINE LOCATE_BLOCK(file,block,block_found)

  IMPLICIT NONE

  !INPUT VARIABLES
  INTEGER,INTENT(IN):: file
  CHARACTER(100),INTENT(IN):: block

  !OUTPUT VARIABLES
  INTEGER,INTENT(OUT):: block_found

  !LOCAL VARIABLES
  CHARACTER(1000):: line
  CHARACTER(100):: key
  INTEGER:: io, io2

  REWIND(file) !go back to start of file
  block_found = 0
  DO
    READ(file,'(A)',IOSTAT=io) line

    IF (io .GT. 0) THEN !io > 0 means bad read
      STOP 'Something went wrong reading file in LOCATE_BLOCK'

    ELSE IF (io .LT. 0) THEN !io < 0 signals end of file
      EXIT

    ELSE !checking which key is specified and putting argument in correct variable
      READ(line,*,IOSTAT=io2) key !reading in key/argument pairs

      IF ((key(1:1) .EQ. '*').OR.(key(1:1) .EQ. '!').OR.(key(1:1) .EQ. '#').OR.(key(1:1) .EQ. ' ')) THEN !detecting commented lines
        CYCLE !if commented line, go to next line

      ELSE IF (trim(key) .EQ. block) THEN !found start of input block
        block_found = 1 !set flag
        EXIT

      ELSE !line not commented, but does not contain block name
        CYCLE

      END IF

    END IF

  END DO

END SUBROUTINE LOCATE_BLOCK

!==================================================================================================================================!
!
!==================================================================================================================================!

SUBROUTINE INPUT_RUN_STATE(inpunit,run_type,restart_infile,use_grey,Res_Calc,Test)

  IMPLICIT NONE

  !INPUT VARIABLES
  INTEGER,INTENT(IN):: inpunit

  !OUTPUT VARIABLES
  INTEGER,INTENT(OUT):: use_grey
  CHARACTER(100),INTENT(OUT):: run_type, restart_infile, Test
  LOGICAL,INTENT(OUT):: Res_Calc

  !LOCAL VARIABLES
  CHARACTER(1000):: line
  CHARACTER(100):: key
  CHARACTER(100):: block
  CHARACTER(100),DIMENSION(3):: args
  INTEGER:: io, io2, block_found

  !DEFAULT VALUES
  run_type = 'mlqd'
  Res_Calc = .TRUE.
  restart_infile = ''
  use_grey = 1
  Test = 'FC'

  block = '[RUN_STATE]'
  CALL LOCATE_BLOCK(inpunit,block,block_found)

  IF (block_found .EQ. 0) STOP '[RUN_STATE] block was not located in input file'

  DO
    READ(inpunit,'(A)',IOSTAT=io) line

    IF (io .GT. 0) THEN !io > 0 means bad read
      STOP 'Something went wrong reading general.inp'

    ELSE IF (io .LT. 0) THEN !io < 0 signals end of file
      EXIT

    ELSE !checking which key is specified and putting argument in correct variable
      READ(line,*,IOSTAT=io2) key, args !reading in key/argument pairs

      IF (key(1:1) .EQ. '[') THEN
        EXIT

      ELSE IF (trim(key) .EQ. 'run_type') THEN
        READ(args(1),*) run_type
        IF ( ALL(run_type .NE. &
        (/'mlqd     ','tr_no_qd ','mg_pod   ','gr_pod   ','p1       ','p13      ','diff     ','fld      '/)) ) THEN
          STOP 'unrecognized run_type (source - subroutine INPUT_RUN_STATE :: module INPUTS)'
        END IF

      ELSE IF (trim(key) .EQ. 'Test') THEN
        READ(args(1),*) Test

      ELSE IF (trim(key) .EQ. 'restart') THEN
        READ(args(1),*) restart_infile

      ELSE IF (trim(key) .EQ. 'use_grey') THEN
        READ(args(1),*) use_grey

      ELSE IF (trim(key) .EQ. 'res_calc') THEN
        READ(args(1),*) io
        IF (io .EQ. 1) THEN
          Res_Calc = .TRUE.
        ELSE
          Res_Calc = .FALSE.
        END IF

      END IF

    END IF

  END DO

  IF (ANY(run_type .EQ. (/'rte_no_qd'/)) ) THEN
    use_grey = 0
  END IF

END SUBROUTINE INPUT_RUN_STATE

!==================================================================================================================================!
!
!==================================================================================================================================!

SUBROUTINE INPUT_SOLVER_OPTS(inpunit,chi,conv_ho,conv_lo,conv_gr1,conv_gr2,comp_unit,line_src,maxit_RTE,maxit_MLOQD,&
  maxit_GLOQD,conv_type,threads,kapE_dT_flag,enrgy_strc,quadrature,E_Bound_Low,T_Bound_Low,Use_Line_Search,&
  Use_Safety_Search,Theta)

  IMPLICIT NONE

  !INPUT VARIABLES
  INTEGER,INTENT(IN):: inpunit

  !OUTPUT VARIABLES
  REAL*8,INTENT(OUT):: chi, conv_ho, conv_lo, conv_gr1, conv_gr2, comp_unit, line_src, Theta
  REAL*8,INTENT(OUT):: E_Bound_Low, T_Bound_Low
  INTEGER,INTENT(OUT):: maxit_RTE, maxit_MLOQD, maxit_GLOQD, conv_type, threads
  CHARACTER(100),INTENT(OUT):: enrgy_strc, quadrature
  LOGICAL:: Use_Line_Search, Use_Safety_Search, kapE_dT_flag

  !LOCAL VARIABLES
  CHARACTER(1000):: line
  CHARACTER(100):: key
  CHARACTER(100):: block
  CHARACTER(100),DIMENSION(3):: args
  INTEGER:: io, io2, block_found
  INTEGER:: conv_ho_flag, conv_lo_flag, conv_gr_flag, bnd_flag

  !DEFAULT VALUES
  maxit_RTE = 25
  maxit_MLOQD = 500
  maxit_GLOQD = 100
  kapE_dT_flag = .TRUE.
  chi = 0.7d0
  Theta = 1d0
  conv_ho = 1d-8
  conv_lo = 1d-9
  conv_gr1 = 1d-10
  conv_gr2 = 1d-15
  conv_type = 1
  quadrature = 'abu36'
  threads = 1
  enrgy_strc = 'JCP'
  comp_unit = 1d13
  Use_Line_Search = .TRUE.
  line_src = 2d0
  Use_Safety_Search = .FALSE.
  E_Bound_Low = -HUGE(1d0)
  T_Bound_Low = -HUGE(1d0)

  conv_ho_flag = 0
  conv_lo_flag = 0
  conv_gr_flag = 0
  bnd_flag = 0

  block = '[SOLVER_OPTS]'
  CALL LOCATE_BLOCK(inpunit,block,block_found)

  IF (block_found .EQ. 0) STOP '[SOLVER_OPTS] block was not located in input file'

  DO
    READ(inpunit,'(A)',IOSTAT=io) line

    IF (io .GT. 0) THEN !io > 0 means bad read
      STOP 'Something went wrong reading general.inp'

    ELSE IF (io .LT. 0) THEN !io < 0 signals end of file
      EXIT

    ELSE !checking which key is specified and putting argument in correct variable
      READ(line,*,IOSTAT=io2) key, args !reading in key/argument pairs

      IF (key(1:1) .EQ. '[') THEN
        EXIT

      ELSE IF (trim(key) .EQ. 'maxit_RTE') THEN
        READ(args(1),*) maxit_RTE
        IF (maxit_RTE .LE. 0) THEN
          WRITE(*,*)
          WRITE(*,'(A)') '*** WARNING ***'
          WRITE(*,'(A)') 'Cannot have zero or less max RTE iterations'
          WRITE(*,'(A)') 'Setting maxit_RTE to default'
          maxit_RTE = 15
        END IF

      ELSE IF (trim(key) .EQ. 'maxit_MLOQD') THEN
        READ(args(1),*) maxit_MLOQD
        IF (maxit_MLOQD .LE. 0) THEN
          WRITE(*,*)
          WRITE(*,'(A)') '*** WARNING ***'
          WRITE(*,'(A)') 'Cannot have zero or less max MLOQD iterations'
          WRITE(*,'(A)') 'Setting maxit_MLOQD to default'
          maxit_MLOQD = 100
        END IF

      ELSE IF (trim(key) .EQ. 'maxit_GLOQD') THEN
        READ(args(1),*) maxit_GLOQD
        IF (maxit_GLOQD .LE. 0) THEN
          WRITE(*,*)
          WRITE(*,'(A)') '*** WARNING ***'
          WRITE(*,'(A)') 'Cannot have zero or less max GLOQD iterations'
          WRITE(*,'(A)') 'Setting maxit_GLOQD to default'
          maxit_GLOQD = 100
        END IF

      ELSE IF (trim(key) .EQ. 'chi') THEN
        READ(args(1),*) chi
        IF (chi .LT. 0) THEN
          WRITE(*,*)
          WRITE(*,'(A)') '*** WARNING ***'
          WRITE(*,'(A)') '"chi" cannot have a negative value'
          WRITE(*,'(A)') 'Setting chi to default (0.7)'
          chi = 0.7d0
        ELSE IF (chi .GE. 1) THEN
          WRITE(*,*)
          WRITE(*,'(A)') '*** WARNING ***'
          WRITE(*,'(A)') '"chi" cannot have a value greater than or equal to 1'
          WRITE(*,'(A)') 'Setting chi to default (0.7)'
          chi = 0.7d0
        END IF

      ELSE IF (trim(key) .EQ. 'theta') THEN
        READ(args(1),*) Theta
        IF (Theta .LT. 0) THEN
          WRITE(*,*)
          WRITE(*,'(A)') '*** WARNING ***'
          WRITE(*,'(A)') '"Theta" cannot have a negative value'
          WRITE(*,'(A)') 'Setting Theta to default (1)'
          Theta = 1d0
        ELSE IF (Theta .GT. 1) THEN
          WRITE(*,*)
          WRITE(*,'(A)') '*** WARNING ***'
          WRITE(*,'(A)') '"Theta" cannot have a value greater than 1'
          WRITE(*,'(A)') 'Setting Theta to default (1)'
          Theta = 1d0
        ELSE IF (Theta .EQ. 0) THEN
          WRITE(*,*)
          WRITE(*,'(A)') '*** WARNING ***'
          WRITE(*,'(A)') '"Theta" cannot be =0 due to 1/Theta terms included in the discretization'
          WRITE(*,'(A)') 'Setting Theta to 0.00001, but continue with caution - this may cause numerical issues'
          Theta = 1d-5
        END IF

      ELSE IF (trim(key) .EQ. 'conv_ho') THEN
        READ(args(1),*) conv_ho
        conv_ho_flag = 1
        IF (conv_ho .LE. 0) THEN
          WRITE(*,*)
          WRITE(*,'(A)') '*** WARNING ***'
          WRITE(*,'(A)') '"conv_ho" cannot have a negative value'
          WRITE(*,'(A)') 'Setting conv_ho to default'
          conv_ho = 1d-8
        ELSE IF (conv_ho .GE. 1) THEN
          WRITE(*,*)
          WRITE(*,'(A)') '*** WARNING ***'
          WRITE(*,'(A)') '"conv_ho" cannot have a value greater than 1'
          WRITE(*,'(A)') 'Setting conv_ho to default'
          conv_ho = 1d-8
        END IF

      ELSE IF (trim(key) .EQ. 'conv_lo') THEN
        READ(args(1),*) conv_lo
        conv_lo_flag = 1
        IF (conv_lo .LE. 0) THEN
          WRITE(*,*)
          WRITE(*,'(A)') '*** WARNING ***'
          WRITE(*,'(A)') '"conv_lo" cannot have a negative value'
          WRITE(*,'(A)') 'Setting conv_lo to default'
          conv_lo = 1d-9
        ELSE IF (conv_lo .GE. 1) THEN
          WRITE(*,*)
          WRITE(*,'(A)') '*** WARNING ***'
          WRITE(*,'(A)') '"conv_lo" cannot have a value greater than 1'
          WRITE(*,'(A)') 'Setting conv_lo to default'
          conv_lo = 1d-9
        END IF

      ELSE IF (trim(key) .EQ. 'conv_gr') THEN
        READ(args(1),*) conv_gr1
        READ(args(2),*) conv_gr2
        conv_gr_flag = 1
        IF ((conv_gr1 .LE. 0).OR.(conv_gr2 .LE. 0)) THEN
          WRITE(*,*)
          WRITE(*,'(A)') '*** WARNING ***'
          WRITE(*,'(A)') '"conv_gr" cannot have a negative value'
          WRITE(*,'(A)') 'Setting conv_gr to default'
          conv_gr1 = 1d-10
          conv_gr2 = 1d-15
        ELSE IF ((conv_gr1 .GE. 1).OR.(conv_gr2 .GE. 1)) THEN
          WRITE(*,*)
          WRITE(*,'(A)') '*** WARNING ***'
          WRITE(*,'(A)') '"conv_gr" cannot have a value greater than 1'
          WRITE(*,'(A)') 'Setting conv_gr to default'
          conv_gr1 = 1d-10
          conv_gr2 = 1d-15
        END IF

      ELSE IF (trim(key) .EQ. 'conv_type') THEN
        READ(args(1),*) conv_type
        IF ( ALL(conv_type .NE. (/1,2,3/)) ) THEN
          WRITE(*,*)
          WRITE(*,'(A)') '*** WARNING ***'
          WRITE(*,'(A)') 'Wrong value set for "conv_type" input, Omega_xst be contained in [1,2,3]'
          WRITE(*,'(A)') 'Setting conv_type to default'
          conv_type = 1
        END IF

      ELSE IF (trim(key) .EQ. 'quad') THEN
        READ(args(1),*) quadrature

      ELSE IF (trim(key) .EQ. 'threads') THEN
        READ(args(1),*) threads
        IF (threads .LE. 0) THEN
          WRITE(*,*)
          WRITE(*,'(A)') '*** WARNING ***'
          WRITE(*,'(A)') '"threads" cannot have a value less than 1'
          WRITE(*,'(A)') 'Setting threads to default'
          threads = 1
        END IF

      ELSE IF (trim(key) .EQ. 'enrgy_strc') THEN
        READ(args(1),*) enrgy_strc

      ELSE IF (trim(key) .EQ. 'kapE_dT') THEN
        ! READ(args(1),*) kapE_dT_flag
        IF ( ALL(args(1) .NE. (/ 'on ' , 'off' /)) ) THEN
          WRITE(*,*)
          WRITE(*,'(A)') '*** WARNING ***'
          WRITE(*,'(A)') 'unknown flag detected for "kapE_dT" input'
          WRITE(*,'(A)') 'acceptable "kapE_dT" inputs are [on, off]'
          WRITE(*,'(A)') 'Setting kapE_dT_flag to default (off)'
          kapE_dT_flag = .FALSE.
        ELSE IF (args(1) .EQ. 'on') THEN
          kapE_dT_flag = .TRUE.
        ELSE IF (args(1) .EQ. 'off') THEN
          kapE_dT_flag = .FALSE.
        END IF


      ELSE IF (trim(key) .EQ. 'comp_unit') THEN
        READ(args(1),*) comp_unit
        IF (comp_unit .LT. 1) THEN
          WRITE(*,*)
          WRITE(*,'(A)') '*** WARNING ***'
          WRITE(*,'(A)') '"comp_unit" cannot have a value less than 1'
          WRITE(*,'(A)') 'Setting comp_unit to default (1d13)'
          comp_unit = 1d13
        END IF

      ELSE IF (trim(key) .EQ. 'line_src') THEN
        IF ( ANY( trim(args(1)) .EQ. (/'off','Off','OFF'/) ) ) THEN
          Use_Line_Search = .FALSE.
          line_src = 1d0
        ELSE
          Use_Line_Search = .TRUE.
          READ(args(1),*) line_src
        END IF

      ELSE IF (trim(key) .EQ. 'E_low_bnd') THEN
        bnd_flag = bnd_flag + 1
        IF ( ANY( trim(args(1)) .EQ. (/'off','Off','OFF'/) ) ) THEN
          IF (bnd_flag .EQ. 1) Use_Safety_Search = .FALSE.
          E_Bound_Low = -HUGE(1d0)
        ELSE
          Use_Safety_Search = .TRUE.
          READ(args(1),*) E_Bound_Low
        END IF

      ELSE IF (trim(key) .EQ. 'T_low_bnd') THEN
        bnd_flag = bnd_flag + 1
        IF ( ANY( trim(args(1)) .EQ. (/'off','Off','OFF'/) ) ) THEN
          IF (bnd_flag .EQ. 1) Use_Safety_Search = .FALSE.
          E_Bound_Low = -HUGE(1d0)
        ELSE
          Use_Safety_Search = .TRUE.
          READ(args(1),*) E_Bound_Low
        END IF

      END IF

    END IF

  END DO


  IF (conv_ho_flag .EQ. 1) THEN
    IF (conv_lo_flag .EQ. 1) THEN
      IF (conv_gr_flag .EQ. 0) conv_gr1 = conv_lo/10d0
    ELSE
      conv_lo = conv_ho/10d0
      IF (conv_gr_flag .EQ. 0) conv_gr1 = conv_ho/100d0
    END IF
  ELSE
    IF (conv_lo_flag .EQ. 1) THEN
      conv_ho = conv_lo*10d0
      IF (conv_gr_flag .EQ. 0) conv_gr1 = conv_lo/10d0
    ELSE
      IF (conv_gr_flag .EQ. 1) THEN
        conv_ho = conv_gr1*100d0
        conv_lo = conv_gr1*10d0
      ELSE
        CONTINUE
      END IF
    END IF
  END IF

END SUBROUTINE INPUT_SOLVER_OPTS

!==================================================================================================================================!
!
!==================================================================================================================================!
SUBROUTINE INPUT_POD_OPTS(inpunit,POD_err,PODgsum,POD_Type,POD_dset,Direc_Diff)

  IMPLICIT NONE

  !INPUT VARIABLES
  INTEGER,INTENT(IN):: inpunit

  !OUTPUT VARIABLES
  REAL*8,INTENT(OUT):: POD_err
  INTEGER,INTENT(OUT):: PODgsum,Direc_Diff
  CHARACTER(100),INTENT(OUT):: POD_Type, POD_dset

  !LOCAL VARIABLES
  CHARACTER(1000):: line
  CHARACTER(100):: key
  CHARACTER(100):: block
  CHARACTER(100),DIMENSION(3):: args
  INTEGER:: io, io2, block_found

  !DEFAULT VALUES
  POD_dset = 'PODout.h5'
  POD_err = 1d-5
  PODgsum = 1
  POD_Type = 'fg'
  Direc_Diff = 0

  block = '[POD_OPTS]'
  CALL LOCATE_BLOCK(inpunit,block,block_found)

  IF (block_found .EQ. 0) STOP '[POD_OPTS] block was not located in input file'

  DO
    READ(inpunit,'(A)',IOSTAT=io) line

    IF (io .GT. 0) THEN !io > 0 means bad read
      STOP 'Something went wrong reading general.inp'

    ELSE IF (io .LT. 0) THEN !io < 0 signals end of file
      EXIT

    ELSE !checking which key is specified and putting argument in correct variable
      READ(line,*,IOSTAT=io2) key, args !reading in key/argument pairs

      IF (key(1:1) .EQ. '[') THEN
        EXIT

      ELSE IF (trim(key) .EQ. 'dataset') THEN
        READ(args(1),*) POD_dset

      ELSE IF (trim(key) .EQ. 'POD_type') THEN
        READ(args(1),*) POD_type
        IF ( ALL(POD_type .NE. (/'fg','Ig'/)) ) THEN
          STOP 'unrecognized POD_type (source - subroutine INPUT_POD_OPTS :: module INPUTS)'
        END IF

      ELSE IF (trim(key) .EQ. 'gsum') THEN
        READ(args(1),*) PODgsum
        IF (ALL(PODgsum .NE. (/0,1/))) STOP 'POD_Type must be 0 or 1'

      ELSE IF (trim(key) .EQ. 'POD_err') THEN
        READ(args(1),*) POD_err
        IF ((POD_err .LE. 0).OR.(POD_err .GT. 1)) STOP 'POD_err must be between 0 and 1'

      ELSE IF (trim(key) .EQ. 'direc_diff') THEN
        READ(args(1),*) Direc_Diff
        IF (ALL(Direc_Diff .NE. (/0,1/))) STOP 'direc_diff must be 0 or 1'

      END IF

    END IF

  END DO

END SUBROUTINE INPUT_POD_OPTS

!==================================================================================================================================!
!
!==================================================================================================================================!

SUBROUTINE INPUT_PARAMETERS(inpunit,erg,xlen,ylen,N_x,N_y,tlen,delt,BC_Type,bcT_left,bcT_right,bcT_upper,&
  bcT_lower,Tini)

  IMPLICIT NONE

  !INPUT VARIABLES
  REAL*8,INTENT(IN):: erg
  INTEGER,INTENT(IN):: inpunit

  !OUTPUT VARIABLES
  REAL*8,INTENT(OUT):: xlen, ylen, tlen, delt, bcT_left, bcT_right, bcT_upper, bcT_lower, Tini
  INTEGER,INTENT(OUT):: N_x, N_y, BC_Type(:)

  !LOCAL VARIABLES
  CHARACTER(1000):: line
  CHARACTER(100):: key, unit
  CHARACTER(100):: block
  CHARACTER(100),DIMENSION(4):: args
  INTEGER:: io, io2, block_found

  !DEFAULT VALUES
  xlen = 6d0
  ylen = 6d0
  N_x = 10
  N_y = 10
  tlen = 6d-1
  delt = 2d-3
  bcT_left = 1d3
  bcT_right = 0d0
  bcT_upper = 0d0
  bcT_lower = 0d0
  Tini = 1d0
  BC_Type = (/0,0,0,0/)

  block = '[PARAMETERS]'
  CALL LOCATE_BLOCK(inpunit,block,block_found)

  IF (block_found .EQ. 0) STOP '[PARAMETERS] block was not located in input file'

  DO
    READ(inpunit,'(A)',IOSTAT=io) line

    IF (io .GT. 0) THEN !io > 0 means bad read
      STOP 'Something went wrong reading general.inp'

    ELSE IF (io .LT. 0) THEN !io < 0 signals end of file
      EXIT

    ELSE !checking which key is specified and putting argument in correct variable
      READ(line,*,IOSTAT=io2) key, args !reading in key/argument pairs

      IF (key(1:1) .EQ. '[') THEN
        EXIT

      ELSE IF (trim(key) .EQ. 'xlen') THEN
        READ(args(1),*) xlen
        READ(args(2),*,IOSTAT=io2) unit
        IF (unit .EQ. 'cm') THEN
          xlen = xlen
        ELSE IF (unit .EQ. 'm') THEN
          xlen = xlen*100d0
        ELSE IF (unit .EQ. 'mm') THEN
          xlen = xlen/10d0
        ELSE IF (unit .EQ. '!singular') THEN
          WRITE(*,*)
          WRITE(*,'(A)') '*** WARNING ***'
          WRITE(*,'(A)') 'NO UNITS GIVEN FOR "xlen", ASSUMING cm'
          WRITE(*,*)
        ELSE
          WRITE(*,'(3A)') 'UNKNOWN UNITS (',trim(unit),') GIVEN FOR "xlen"'
          WRITE(*,'(A)') 'ACCEPTABLE UNITS ARE [cm, m, mm]'
          STOP
        END IF

      ELSE IF (trim(key) .EQ. 'ylen') THEN
        READ(args(1),*) ylen
        READ(args(2),*,IOSTAT=io2) unit
        IF (unit .EQ. 'cm') THEN
          ylen = ylen
        ELSE IF (unit .EQ. 'm') THEN
          ylen = ylen*100d0
        ELSE IF (unit .EQ. 'mm') THEN
          ylen = ylen/10d0
        ELSE IF (unit .EQ. '!singular') THEN
          WRITE(*,*)
          WRITE(*,'(A)') '*** WARNING ***'
          WRITE(*,'(A)') 'NO UNITS GIVEN FOR "ylen", ASSUMING cm'
          WRITE(*,*)
        ELSE
          WRITE(*,'(3A)') 'UNKNOWN UNITS (',trim(unit),') GIVEN FOR "ylen"'
          WRITE(*,'(A)') 'ACCEPTABLE UNITS ARE [cm, m, mm]'
          STOP
        END IF

      ELSE IF (trim(key) .EQ. 'x_cells') THEN
        READ(args(1),*) N_x
        IF (N_x .LE. 0) THEN
          WRITE(*,*)
          WRITE(*,*) 'x_cells be a positive integer (>0)'
          WRITE(*,*) 'resetting x_cells = 10'
          N_x = 10
        END IF

      ELSE IF (trim(key) .EQ. 'y_cells') THEN
        READ(args(1),*) N_y
        IF (N_y .LE. 0) THEN
          WRITE(*,*)
          WRITE(*,*) 'y_cells be a positive integer (>0)'
          WRITE(*,*) 'resetting y_cells = 10'
          N_y = 10
        END IF

      ELSE IF (trim(key) .EQ. 'tlen') THEN
        READ(args(1),*) tlen
        READ(args(2),*,IOSTAT=io2) unit
        IF (unit .EQ. 'sh') THEN
          tlen = tlen
        ELSE IF (unit .EQ. 'ns') THEN
          tlen = tlen*10d0
        ELSE IF (unit .EQ. 'sec') THEN
          tlen = tlen/1d8
        ELSE IF (unit .EQ. '!singular') THEN
          WRITE(*,*)
          WRITE(*,'(A)') '*** WARNING ***'
          WRITE(*,'(A)') 'NO UNITS GIVEN FOR "tlen", ASSUMING sh'
          WRITE(*,*)
        ELSE
          WRITE(*,'(3A)') 'UNKNOWN UNITS (',trim(unit),') GIVEN FOR "tlen"'
          WRITE(*,'(A)') 'ACCEPTABLE UNITS ARE [sh, ns, sec]'
          STOP
        END IF

      ELSE IF (trim(key) .EQ. 'delt') THEN
        READ(args(1),*) delt
        READ(args(2),*,IOSTAT=io2) unit
        IF (unit .EQ. 'sh') THEN
          delt = delt
        ELSE IF (unit .EQ. 'ns') THEN
          delt = delt*10d0
        ELSE IF (unit .EQ. 'sec') THEN
          delt = delt/1d8
        ELSE IF (unit .EQ. '!singular') THEN
          WRITE(*,*)
          WRITE(*,'(A)') '*** WARNING ***'
          WRITE(*,'(A)') 'NO UNITS GIVEN FOR "delt", ASSUMING sh'
          WRITE(*,*)
        ELSE
          WRITE(*,'(3A)') 'UNKNOWN UNITS (',trim(unit),') GIVEN FOR "delt"'
          WRITE(*,'(A)') 'ACCEPTABLE UNITS ARE [sh, ns, sec]'
          STOP
        END IF

      ELSE IF (trim(key) .EQ. 'bc_type') THEN
        READ(args(1:4),*) BC_Type
        IF (MAXVAL(BC_Type) .LT. 0) THEN
          STOP 'bc_type cannot be less than 0'
        ELSE IF (MAXVAL(BC_Type) .GT. 2) THEN
            STOP 'bc_type cannot be greater than 2'
        END IF

      ELSE IF (trim(key) .EQ. 'bcT_left') THEN
        READ(args(1),*) bcT_left
        READ(args(2),*,IOSTAT=io2) unit
        IF (unit .EQ. 'erg') THEN
          bcT_left = bcT_left*erg
        ELSE IF (unit .EQ. 'ev') THEN
          bcT_left = bcT_left
        ELSE IF (unit .EQ. '!singular') THEN
          WRITE(*,*)
          WRITE(*,'(A)') '*** WARNING ***'
          WRITE(*,'(A)') 'NO UNITS GIVEN FOR "bcT_left", ASSUMING ev'
          WRITE(*,*)
          bcT_left = bcT_left
        ELSE
          WRITE(*,'(3A)') 'UNKNOWN UNITS (',trim(unit),') GIVEN FOR "bcT_left"'
          WRITE(*,'(A)') 'ACCEPTABLE UNITS ARE [erg, ev]'
          STOP
        END IF
        ! Cv=0.5917d0*ar*(bcl)**3*((1.38d-16)**4*(11600d0**4)*erg**4) !technically the ((1.38d-16)**4*(11600d0**4)*erg**4) should =1 but it doesnt
        ! cv=8.0901395d10                                           !taking it out would be more correct, but it gets consistent results with dima so fuck it

      ELSE IF (trim(key) .EQ. 'bcT_right') THEN
        READ(args(1),*) bcT_right
        READ(args(2),*,IOSTAT=io2) unit
        IF (unit .EQ. 'erg') THEN
          bcT_right = bcT_right*erg
        ELSE IF (unit .EQ. 'ev') THEN
          bcT_right = bcT_right
        ELSE IF (unit .EQ. '!singular') THEN
          WRITE(*,*)
          WRITE(*,'(A)') '*** WARNING ***'
          WRITE(*,'(A)') 'NO UNITS GIVEN FOR "bcT_right", ASSUMING ev'
          WRITE(*,*)
          bcT_right = bcT_right
        ELSE
          WRITE(*,'(3A)') 'UNKNOWN UNITS (',trim(unit),') GIVEN FOR "bcT_right"'
          WRITE(*,'(A)') 'ACCEPTABLE UNITS ARE [erg, ev]'
          STOP
        END IF

      ELSE IF (trim(key) .EQ. 'bcT_upper') THEN
        READ(args(1),*) bcT_upper
        READ(args(2),*,IOSTAT=io2) unit
        IF (unit .EQ. 'erg') THEN
          bcT_upper = bcT_upper*erg
        ELSE IF (unit .EQ. 'ev') THEN
          bcT_upper = bcT_upper
        ELSE IF (unit .EQ. '!singular') THEN
          WRITE(*,*)
          WRITE(*,'(A)') '*** WARNING ***'
          WRITE(*,'(A)') 'NO UNITS GIVEN FOR "bcT_upper", ASSUMING ev'
          WRITE(*,*)
          bcT_upper = bcT_upper
        ELSE
          WRITE(*,'(3A)') 'UNKNOWN UNITS (',trim(unit),') GIVEN FOR "bcT_upper"'
          WRITE(*,'(A)') 'ACCEPTABLE UNITS ARE [erg, ev]'
          STOP
        END IF

      ELSE IF (trim(key) .EQ. 'bcT_lower') THEN
        READ(args(1),*) bcT_lower
        READ(args(2),*,IOSTAT=io2) unit
        IF (unit .EQ. 'erg') THEN
          bcT_lower = bcT_lower*erg
        ELSE IF (unit .EQ. 'ev') THEN
          bcT_lower = bcT_lower
        ELSE IF (unit .EQ. '!singular') THEN
          WRITE(*,*)
          WRITE(*,'(A)') '*** WARNING ***'
          WRITE(*,'(A)') 'NO UNITS GIVEN FOR "bcT_lower", ASSUMING ev'
          WRITE(*,*)
          bcT_lower = bcT_lower
        ELSE
          WRITE(*,'(3A)') 'UNKNOWN UNITS (',trim(unit),') GIVEN FOR "bcT_lower"'
          WRITE(*,'(A)') 'ACCEPTABLE UNITS ARE [erg, ev]'
          STOP
        END IF

      ELSE IF (trim(key) .EQ. 'Tini') THEN
        READ(args(1),*) Tini
        READ(args(2),*,IOSTAT=io2) unit
        IF (unit .EQ. 'erg') THEN
          Tini = Tini*erg
        ELSE IF (unit .EQ. 'ev') THEN
          Tini = Tini
        ELSE IF (unit .EQ. '!singular') THEN
          WRITE(*,*)
          WRITE(*,'(A)') '*** WARNING ***'
          WRITE(*,'(A)') 'NO UNITS GIVEN FOR "Tini", ASSUMING ev'
          WRITE(*,*)
          Tini = Tini
        ELSE
          WRITE(*,'(3A)') 'UNKNOWN UNITS (',trim(unit),') GIVEN FOR "Tini"'
          WRITE(*,'(A)') 'ACCEPTABLE UNITS ARE [erg, ev]'
          STOP
        END IF

      END IF

    END IF

  END DO

  ! IF (xlen .LE. 0d0) STOP 'Negative or zero value for "xlen" was detected in general.inp, aborting program'
  ! IF (cell_xl .LE. 0d0) STOP 'Negative or zero value for "delx" was detected in general.inp, aborting program'
  ! IF (tlen .LE. 0d0) STOP 'Negative or zero value for "tlen" was detected in general.inp, aborting program'
  ! IF (delt .LE. 0d0) STOP 'Negative or zero value for "delt" was detected in general.inp, aborting program'
  ! IF (Tini .LE. 0d0) STOP 'Negative or zero value for "Tini" was detected in general.inp, aborting program'
  ! IF (bcl .LT. 0d0) STOP 'Negative value for "bcl" was detected in general.inp, aborting program'
  ! IF (bcr .LT. 0d0) STOP 'Negative value for "bcr" was detected in general.inp, aborting program'

END SUBROUTINE INPUT_PARAMETERS

!==================================================================================================================================!
!
!==================================================================================================================================!

SUBROUTINE INPUT_OUTPUT_OPTS(inpunit,outfile,out_freq,I_out,HO_Eg_out,HO_Fg_out,HO_E_out,HO_F_out,Eg_out,Fg_out,MGQD_E_out,&
  MGQD_F_out,QDfg_out,E_out,F_out,D_out,old_parms_out,its_out,conv_out,kap_out,Src_out)

  IMPLICIT NONE

  !INPUT VARIABLES
  INTEGER,INTENT(IN):: inpunit

  !OUTPUT VARIABLES
  CHARACTER(100),INTENT(OUT):: outfile
  INTEGER,INTENT(OUT):: out_freq, I_out, HO_Eg_out, HO_Fg_out, HO_E_out, HO_F_out
  INTEGER,INTENT(OUT):: Eg_out, Fg_out, MGQD_E_out, MGQD_F_out, QDfg_out
  INTEGER,INTENT(OUT):: E_out, F_out, D_out
  INTEGER,INTENT(OUT):: old_parms_out, its_out, conv_out, kap_out, Src_out

  !LOCAL VARIABLES
  CHARACTER(1000):: line
  CHARACTER(100):: key
  CHARACTER(100):: block
  CHARACTER(100),DIMENSION(3):: args
  INTEGER:: io, io2, block_found

  !DEFAULT VALUES
  outfile = 'output.h5'
  out_freq  = 1
  I_out     = 1
  HO_Eg_out = 1
  HO_Fg_out = 1
  HO_E_out  = 1
  HO_F_out  = 1
  Eg_out    = 1
  Fg_out    = 1
  MGQD_E_out = 1
  MGQD_F_out = 1
  QDfg_out    = 1
  E_out     = 1
  F_out     = 1
  D_out     = 1
  old_parms_out = 1
  its_out  = 1
  conv_out = 1
  kap_out  = 1
  Src_out  = 1

  block = '[OUTPUT_OPTS]'
  CALL LOCATE_BLOCK(inpunit,block,block_found)

  IF (block_found .EQ. 0) STOP '[OUTPUT_OPTS] block was not located in input file'

  DO
    READ(inpunit,'(A)',IOSTAT=io) line

    IF (io .GT. 0) THEN !io > 0 means bad read
      STOP 'Something went wrong reading general.inp'

    ELSE IF (io .LT. 0) THEN !io < 0 signals end of file
      EXIT

    ELSE !checking which key is specified and putting argument in correct variable
      READ(line,*,IOSTAT=io2) key, args !reading in key/argument pairs

      IF (key(1:1) .EQ. '[') THEN
        EXIT

      ELSE IF (trim(key) .EQ. 'outfile') THEN
        READ(args(1),*) outfile

      ELSE IF (trim(key) .EQ. 'out_freq') THEN
        READ(args(1),*) out_freq
        IF (ALL(out_freq .NE. (/0,1/))) THEN
          WRITE(*,*)
          WRITE(*,'(A)') '*** WARNING ***'
          WRITE(*,'(A)') 'out_freq must be 0 or 1'
          WRITE(*,'(A)') 'Setting out_freq to 0'
          out_freq = 0
        END IF

      ELSE IF (trim(key) .EQ. 'I_out') THEN
        READ(args(1),*) I_out
        IF (ALL(I_out .NE. (/0,1/))) THEN
          WRITE(*,*)
          WRITE(*,'(A)') '*** WARNING ***'
          WRITE(*,'(A)') 'I_out must be 0 or 1'
          WRITE(*,'(A)') 'Setting I_out to 0'
          I_out = 0
        END IF

      ELSE IF (trim(key) .EQ. 'HO_Eg_out') THEN
        READ(args(1),*) HO_Eg_out
        IF (ALL(HO_Eg_out .NE. (/0,1/))) THEN
          WRITE(*,*)
          WRITE(*,'(A)') '*** WARNING ***'
          WRITE(*,'(A)') 'HO_Eg_out must be 0 or 1'
          WRITE(*,'(A)') 'Setting HO_Eg_out to 0'
          HO_Eg_out = 0
        END IF

      ELSE IF (trim(key) .EQ. 'HO_Fg_out') THEN
        READ(args(1),*) HO_Fg_out
        IF (ALL(HO_Fg_out .NE. (/0,1/))) THEN
          WRITE(*,*)
          WRITE(*,'(A)') '*** WARNING ***'
          WRITE(*,'(A)') 'HO_Fg_out must be 0 or 1'
          WRITE(*,'(A)') 'Setting HO_Fg_out to 0'
          HO_Fg_out = 0
        END IF

      ELSE IF (trim(key) .EQ. 'HO_E_out') THEN
        READ(args(1),*) HO_E_out
        IF (ALL(HO_E_out .NE. (/0,1/))) THEN
          WRITE(*,*)
          WRITE(*,'(A)') '*** WARNING ***'
          WRITE(*,'(A)') 'HO_E_out must be 0 or 1'
          WRITE(*,'(A)') 'Setting HO_E_out to 0'
          HO_E_out = 0
        END IF

      ELSE IF (trim(key) .EQ. 'HO_F_out') THEN
        READ(args(1),*) HO_F_out
        IF (ALL(HO_F_out .NE. (/0,1/))) THEN
          WRITE(*,*)
          WRITE(*,'(A)') '*** WARNING ***'
          WRITE(*,'(A)') 'HO_F_out must be 0 or 1'
          WRITE(*,'(A)') 'Setting HO_F_out to 0'
          HO_F_out = 0
        END IF

      ELSE IF (trim(key) .EQ. 'Eg_out') THEN
        READ(args(1),*) Eg_out
        IF (ALL(HO_F_out .NE. (/0,1/))) THEN
          WRITE(*,*)
          WRITE(*,'(A)') '*** WARNING ***'
          WRITE(*,'(A)') 'Eg_out must be 0 or 1'
          WRITE(*,'(A)') 'Setting Eg_out to 0'
          Eg_out = 0
        END IF

      ELSE IF (trim(key) .EQ. 'Fg_out') THEN
        READ(args(1),*) Fg_out
        IF (ALL(HO_F_out .NE. (/0,1/))) THEN
          WRITE(*,*)
          WRITE(*,'(A)') '*** WARNING ***'
          WRITE(*,'(A)') 'Fg_out must be 0 or 1'
          WRITE(*,'(A)') 'Setting Fg_out to 0'
          Fg_out = 0
        END IF

      ELSE IF (trim(key) .EQ. 'MGQD_E_out') THEN
        READ(args(1),*) MGQD_E_out
        IF (ALL(MGQD_E_out .NE. (/0,1/))) THEN
          WRITE(*,*)
          WRITE(*,'(A)') '*** WARNING ***'
          WRITE(*,'(A)') 'MGQD_E_out must be 0 or 1'
          WRITE(*,'(A)') 'Setting MGQD_E_out to 0'
          MGQD_E_out = 0
        END IF

      ELSE IF (trim(key) .EQ. 'MGQD_F_out') THEN
        READ(args(1),*) MGQD_F_out
        IF (ALL(MGQD_F_out .NE. (/0,1/))) THEN
          WRITE(*,*)
          WRITE(*,'(A)') '*** WARNING ***'
          WRITE(*,'(A)') 'MGQD_F_out must be 0 or 1'
          WRITE(*,'(A)') 'Setting MGQD_F_out to 0'
          MGQD_F_out = 0
        END IF

      ELSE IF (trim(key) .EQ. 'fg_out') THEN
        READ(args(1),*) QDfg_out
        IF (ALL(QDfg_out .NE. (/0,1/))) THEN
          WRITE(*,*)
          WRITE(*,'(A)') '*** WARNING ***'
          WRITE(*,'(A)') 'fg_out must be 0 or 1'
          WRITE(*,'(A)') 'Setting fg_out to 0'
          QDfg_out = 0
        END IF

      ELSE IF (trim(key) .EQ. 'E_out') THEN
        READ(args(1),*) E_out
        IF (ALL(E_out .NE. (/0,1/))) THEN
          WRITE(*,*)
          WRITE(*,'(A)') '*** WARNING ***'
          WRITE(*,'(A)') 'E_out must be 0 or 1'
          WRITE(*,'(A)') 'Setting E_out to 0'
          E_out = 0
        END IF

      ELSE IF (trim(key) .EQ. 'F_out') THEN
        READ(args(1),*) F_out
        IF (ALL(F_out .NE. (/0,1/))) THEN
          WRITE(*,*)
          WRITE(*,'(A)') '*** WARNING ***'
          WRITE(*,'(A)') 'F_out must be 0 or 1'
          WRITE(*,'(A)') 'Setting F_out to 0'
          F_out = 0
        END IF

      ELSE IF (trim(key) .EQ. 'D_out') THEN
        READ(args(1),*) D_out
        IF (ALL(D_out .NE. (/0,1/))) THEN
          WRITE(*,*)
          WRITE(*,'(A)') '*** WARNING ***'
          WRITE(*,'(A)') 'D_out must be 0 or 1'
          WRITE(*,'(A)') 'Setting D_out to 0'
          D_out = 0
        END IF

      ELSE IF (trim(key) .EQ. 'old_parms_out') THEN
        READ(args(1),*) old_parms_out
        IF (ALL(old_parms_out .NE. (/0,1/))) THEN
          WRITE(*,*)
          WRITE(*,'(A)') '*** WARNING ***'
          WRITE(*,'(A)') 'old_parms_out must be 0 or 1'
          WRITE(*,'(A)') 'Setting old_parms_out to 0'
          old_parms_out = 0
        END IF

      ELSE IF (trim(key) .EQ. 'its_out') THEN
        READ(args(1),*) its_out
        IF (ALL(its_out .NE. (/0,1/))) THEN
          WRITE(*,*)
          WRITE(*,'(A)') '*** WARNING ***'
          WRITE(*,'(A)') 'its_out must be 0 or 1'
          WRITE(*,'(A)') 'Setting its_out to 0'
          its_out = 0
        END IF

      ELSE IF (trim(key) .EQ. 'conv_out') THEN
        READ(args(1),*) conv_out
        IF (ALL(conv_out .NE. (/0,1/))) THEN
          WRITE(*,*)
          WRITE(*,'(A)') '*** WARNING ***'
          WRITE(*,'(A)') 'conv_out must be 0 or 1'
          WRITE(*,'(A)') 'Setting conv_out to 0'
          conv_out = 0
        END IF

      ELSE IF (trim(key) .EQ. 'kap_out') THEN
        READ(args(1),*) kap_out
        IF (ALL(conv_out .NE. (/0,1/))) THEN
          WRITE(*,*)
          WRITE(*,'(A)') '*** WARNING ***'
          WRITE(*,'(A)') 'kap_out must be 0 or 1'
          WRITE(*,'(A)') 'Setting kap_out to 0'
          kap_out = 0
        END IF

      ELSE IF (trim(key) .EQ. 'Src_out') THEN
        READ(args(1),*) Src_out
        IF (ALL(Src_out .NE. (/0,1/))) THEN
          WRITE(*,*)
          WRITE(*,'(A)') '*** WARNING ***'
          WRITE(*,'(A)') 'Src_out must be 0 or 1'
          WRITE(*,'(A)') 'Setting Src_out to 0'
          Src_out = 0
        END IF

      END IF

    END IF

  END DO

END SUBROUTINE INPUT_OUTPUT_OPTS

!==================================================================================================================================!
!
!==================================================================================================================================!

SUBROUTINE INPUT_ENERGY(inpunit,enrgy_strc,N_g,nu_g)

  IMPLICIT NONE
  INTEGER,INTENT(IN):: inpunit
  CHARACTER(100),INTENT(INOUT):: enrgy_strc

  REAL*8,ALLOCATABLE,INTENT(OUT):: nu_g(:)
  INTEGER,INTENT(OUT):: N_g

  CHARACTER(100):: block
  INTEGER:: g, block_found

  IF (trim(enrgy_strc) .EQ. 'custom') THEN

    block = '[ENERGY_STRUCTURE]'
    CALL LOCATE_BLOCK(inpunit,block,block_found)

    IF (block_found .EQ. 0) STOP '[ENERGY_STRUCTURE] block was not located in input file'

    READ(inpunit,*) N_g
    ALLOCATE(nu_g(N_g+1))

    nu_g(1)=0d0
    DO g=1,N_g
      READ(10,*) nu_g(g+1)
    END DO

  ELSE IF ((trim(enrgy_strc) .EQ. 'JCP').OR.(trim(enrgy_strc) .EQ. 'jcp')) THEN

    N_g = 17
    ALLOCATE(nu_g(N_g+1))
    nu_g = (/ 0d0, 7.075d2, 1.415d3, 2.123d3, 2.830d3, 3.538d3, 4.245d3, 5.129d3, &
              6.014d3, 6.898d3, 7.783d3, 8.667d3, 9.551d3, 1.044d4, 1.132d4, 1.220d4, &
              1.309d4, 1d10/)

  ELSE IF ((trim(enrgy_strc) .EQ. 'NE795').OR.(trim(enrgy_strc) .EQ. 'ne795')) THEN

    N_g = 17
    ALLOCATE(nu_g(N_g+1))
    nu_g = (/ 0d0, 1.9047d0, 3.6278d0, 6.9097d0, 13.161d0, 25.067d0, 47.744d0, 90.937d0, &
              173.21d0, 329.90d0, 628.35d0, 1196.8d0, 2279.5d0, 4341.7d0, 8269.5d0, 15751d0, &
              3d4, 1d10/)

  ELSE IF ((trim(enrgy_strc) .EQ. 'ONE_GROUP').OR.(trim(enrgy_strc) .EQ. 'one_group')) THEN

    N_g = 1
    ALLOCATE(nu_g(N_g+1))
    nu_g = (/ 0d0, 1d10/) !Should never reference nu_g(3)
    enrgy_strc = 'ONE_GROUP' !This will be used as a flag in the code

  ELSE

    STOP 'Invalid "enrgy_strc"'

  END IF

END SUBROUTINE INPUT_ENERGY

!==================================================================================================================================!
!
!==================================================================================================================================!

SUBROUTINE INPUT_QUAD(quadrature,N_m,Omega_x,Omega_y,quad_weight)

  IMPLICIT NONE
  INTEGER,INTENT(INOUT):: N_m
  REAL*8,ALLOCATABLE,INTENT(OUT):: Omega_x(:), Omega_y(:), quad_weight(:)
  CHARACTER(100),INTENT(OUT):: quadrature
  REAL*8,PARAMETER:: pi=4d0*ATAN(1d0)

  IF (quadrature .EQ. 'abu36') THEN
    N_m = 36*4
  ELSE IF(quadrature .EQ. 'abu20') THEN
    N_m = 20*4
  ELSE IF(quadrature .EQ. '1d_test') THEN
    N_m = 8
  ELSE
    STOP 'unrecognized input for quad'
  END IF

  ALLOCATE(Omega_x(N_m))
  ALLOCATE(Omega_y(N_m))
  ALLOCATE(quad_weight(N_m))

  IF (quadrature .EQ. 'abu36') THEN
    Omega_x(1:N_m/4) = (/9.717784813336d-1, 7.656319455497d-1, 4.445439440056d-1, 9.701698603928d-1,&
      7.633693960835d-1, 1.483114568272d-1, 9.622473153642d-1, 7.524467626583d-1, 9.410672772109d-1,&
      4.288508824476d-1, 7.241384940891d-1, 8.997294996538d-1, 6.711819639118d-1, 1.293388490485d-1,&
      8.335743322378d-1, 3.670788892962d-1, 5.909368760506d-1, 7.417637460141d-1, 6.279014865859d-1,&
      4.869502395267d-1, 2.519357740235d-1, 4.996274255819d-1, 7.447663982495d-2, 3.685670882907d-1,&
      3.673697853806d-1, 2.488983098158d-1, 1.196054590036d-1, 2.426233944222d-1, 1.419535016004d-1,&
      1.362124657777d-1, 1.670387759191d-2, 5.995074957044d-2, 5.695764868253d-2, 2.447911451942d-2,&
      1.160393058611d-2, 1.096881837272d-2/)

    Omega_y(1:N_m/4) = (/1.096881837272d-2, 1.160393058611d-2, 2.447911451942d-2, 5.695764868253d-2,&
      5.995074957044d-2, 1.670387759191d-2, 1.362124657777d-1, 1.419535016004d-1, 2.426233944222d-1,&
      1.196054590036d-1, 2.488983098158d-1, 3.673697853806d-1, 3.685670882907d-1, 7.447663982495d-2,&
      4.996274255819d-1, 2.519357740235d-1, 4.869502395267d-1, 6.279014865859d-1, 7.417637460141d-1,&
      5.909368760506d-1, 3.670788892962d-1, 8.335743322378d-1, 1.293388490485d-1, 6.711819639118d-1,&
      8.997294996538d-1, 7.241384940891d-1, 4.288508824476d-1, 9.410672772109d-1, 7.524467626583d-1,&
      9.622473153642d-1, 1.483114568272d-1, 7.633693960835d-1, 9.701698603928d-1, 4.445439440056d-1,&
      7.656319455497d-1, 9.717784813336d-1/)

    quad_weight(1:N_m/4) = (/8.454511187252d-3, 8.352354145856d-3, 1.460888798152d-2, 1.913728513580d-2,&
      1.873220073879d-2, 6.404244616724d-3, 2.863542971348d-2, 2.759429759588d-2, 3.648716160597d-2,&
      2.995376809966d-2, 3.442681426024d-2, 4.244873302980d-2, 3.901232700510d-2, 1.162080754372d-2,&
      4.642823955812d-2, 3.798783310581d-2, 4.130171453748d-2, 4.841339013884d-2, 4.841339013884d-2,&
      4.130171453748d-2, 3.798783310581d-2, 4.642823955812d-2, 1.162080754372d-2, 3.901232700510d-2,&
      4.244873302980d-2, 3.442681426024d-2, 2.995376809966d-2, 3.648716160597d-2, 2.759429759588d-2,&
      2.863542971348d-2, 6.404244616724d-3, 1.873220073879d-2, 1.913728513580d-2, 1.460888798152d-2,&
      8.352354145856d-3, 8.454511187252d-3/)

  ELSE IF(quadrature .EQ. 'abu20') THEN
    Omega_x(1:N_m/4) = (/9.713274064903d-1, 7.645615896150d-1, 4.424202396002d-1, 9.586898685237d-1,&
      7.375714298063d-1, 1.409476441875d-1, 9.028558915298d-1, 3.858240341629d-1, 6.313311043797d-1,&
      7.770210099715d-1, 5.837054752370d-1, 4.332989313333d-1, 2.221674140412d-1, 3.596178122512d-1,&
      4.908227124734d-2, 2.057068622698d-1, 1.593344524838d-1, 4.982847370367d-2, 4.210110375297d-2,&
      3.157215799340d-2/)

    Omega_y(1:N_m/4) = (/3.157215799340d-2, 4.210110375297d-2, 4.982847370367d-2, 1.593344524838d-1,&
      2.057068622698d-1, 4.908227124734d-2, 3.596178122512d-1, 2.221674140412d-1, 4.332989313333d-1,&
      5.837054752370d-1, 7.770210099715d-1, 6.313311043797d-1, 3.858240341629d-1, 9.028558915298d-1,&
      1.409476441875d-1, 7.375714298063d-1, 9.586898685237d-1, 4.424202396002d-1, 7.645615896150d-1,&
      9.713274064903d-1/)

    quad_weight(1:N_m/4) = (/2.419260514149d-2, 2.998205782366d-2, 2.932993043666d-2, 5.213067212540d-2,&
      6.147460425028d-2, 1.802505216045d-2, 7.185542471164d-2, 5.322055875020d-2, 7.796304620960d-2,&
      8.182604839076d-2, 8.182604839076d-2, 7.796304620960d-2, 5.322055875020d-2, 7.185542471164d-2,&
      1.802505216045d-2, 6.147460425028d-2, 5.213067212540d-2, 2.932993043666d-2, 2.998205782366d-2,&
      2.419260514149d-2/)

  ELSE IF(quadrature .EQ. '1d_test') THEN
    Omega_x(1:N_m/4) = (/0.92388d0, 0.382683d0/)

    Omega_y(1:N_m/4) = (/0.382683d0, 0.92388d0/)

    ! quad_weight(1:N_m/4) = (/1.5708d0, 1.5708d0/)
    quad_weight(1:N_m/4) = (/.5d0, .5d0/)

  END IF

  Omega_x(N_m/4+1:N_m*2/4) = -Omega_x(1:N_m/4)
  Omega_x(N_m*2/4+1:N_m*3/4) = -Omega_x(1:N_m/4)
  Omega_x(N_m*3/4+1:N_m*4/4) = Omega_x(1:N_m/4)

  Omega_y(N_m/4+1:N_m*2/4) = Omega_y(1:N_m/4)
  Omega_y(N_m*2/4+1:N_m*3/4) = -Omega_y(1:N_m/4)
  Omega_y(N_m*3/4+1:N_m*4/4) = -Omega_y(1:N_m/4)

  quad_weight(N_m/4+1:N_m*2/4) = quad_weight(1:N_m/4)
  quad_weight(N_m*2/4+1:N_m*3/4) = quad_weight(1:N_m/4)
  quad_weight(N_m*3/4+1:N_m*4/4) = quad_weight(1:N_m/4)
  quad_weight = 4d0*pi*quad_weight/SUM(quad_weight(:))

END SUBROUTINE INPUT_QUAD

!==================================================================================================================================!
!
!==================================================================================================================================!

SUBROUTINE INPUT_TEST_TYPE(Test,N_x,N_y,Delx,Dely,xlen,ylen,Mat,Kappa_Mult)

  IMPLICIT NONE
  CHARACTER(*),INTENT(IN):: Test
  INTEGER,INTENT(IN):: N_x, N_y
  REAL*8,INTENT(IN):: Delx(*), Dely(*), xlen, ylen
  INTEGER,ALLOCATABLE,INTENT(OUT):: Mat(:,:)
  REAL*8,ALLOCATABLE,INTENT(OUT):: Kappa_Mult(:)
  INTEGER:: Nmat, xm, ym, i, j
  REAL*8:: xsum, ysum

  IF (ANY(TEST .EQ.(/ 'FC ', 'F-C', 'fc ', 'f-c' /))) THEN
    Nmat = 1
    ALLOCATE(Mat(N_x,N_y), Kappa_Mult(Nmat))
    Mat = 1
    Kappa_Mult = 1d0

  ELSE IF (ANY(TEST .EQ.(/ 'BLOCK', 'BLK  ', 'block', 'blk  ' /))) THEN
    Nmat = 2
    ALLOCATE(Mat(N_x,N_y), Kappa_Mult(Nmat))
    Kappa_Mult(1) = 1d0
    Kappa_Mult(2) = 1d1

    !block is 2.4 x 2.4 cm^2
    ysum = 0d0
   DO j=1,N_y
     ysum = ysum + Dely(j)/2d0
     xsum = 0d0
      DO i=1,N_x
        xsum = xsum + Delx(i)/2d0
        IF ((xsum .LT. xlen/2d0 - 1.2d0).OR.(xsum .GT. xlen/2d0 + 1.2d0)) THEN
          Mat(i,j) = 1

        ELSE
          IF ((ysum .LT. ylen/2d0 - 1.2d0).OR.(ysum .GT. ylen/2d0 + 1.2d0)) THEN
            Mat(i,j) = 1

          ELSE
            Mat(i,j) = 2

          END IF

        END IF
        xsum = xsum + Delx(i)/2d0
      END DO
      ysum = ysum + Dely(j)/2d0
    END DO

  ELSE
    STOP 'Unrecognized Test type in INPUT_TEST_TYPE'

  END IF

END SUBROUTINE INPUT_TEST_TYPE

END MODULE INPUTS
