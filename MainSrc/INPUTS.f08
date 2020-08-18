MODULE INPUTS

  IMPLICIT NONE

CONTAINS

!============================================================================================================!
!
!============================================================================================================!
SUBROUTINE INPUT(database_gen,database_add,run_type,restart_infile,use_grey,chi,conv_ho,conv_lo,conv_gr,&
  comp_unit,line_src,maxit_RTE,maxit_MLOQD,maxit_GLOQD,conv_type,N_m,threads,kapE_dT_flag,enrgy_strc,&
  Nwt_upbnd,erg,xlen,ylen,x_cells,y_cells,tlen,delt,bcT_left,bcT_right,bcT_upper,bcT_lower,Tini,sig_R,ar,&
  pi,c,h,delx,dely,cv,out_times,out_time_steps,n_out_times,restart_outlen,TEMP_out,GREY_E_out,GREY_F_out,&
  GREY_kap_out,GREY_fsmall_out,MG_fsmall_out,res_history_out,outfile,restart_outfile,decomp_outfile,TEMP_outfile,&
  GREY_E_outfile,GREY_F_outfile,GREY_kap_outfile,GREY_fsmall_outfile,MG_fsmall_outfile,res_history_outfile,&
  E_ho_out,E_ho_outfile,nu_g,N_g,Omega_x,Omega_y,quad_weight,tpoints,quadrature,BC_Type)

  IMPLICIT NONE

  !INPUT VARIABLES
  REAL*8,INTENT(IN):: erg, pi, c, h

  !OUTPUT VARIABLES
  REAL*8,INTENT(INOUT):: comp_unit

  INTEGER,INTENT(OUT):: database_gen, database_add, use_grey
  CHARACTER(100),INTENT(OUT):: run_type, restart_infile

  REAL*8,INTENT(OUT):: chi, conv_ho, conv_lo, conv_gr, line_src, Nwt_upbnd
  INTEGER,INTENT(OUT):: maxit_RTE, maxit_MLOQD, maxit_GLOQD, conv_type, threads
  CHARACTER(100),INTENT(OUT):: kapE_dT_flag, enrgy_strc, quadrature

  REAL*8,INTENT(OUT):: xlen, ylen, tlen, delt, bcT_left, bcT_right, bcT_upper, bcT_lower, Tini
  REAL*8,ALLOCATABLE,INTENT(OUT):: Delx(:), Dely(:)
  INTEGER,INTENT(OUT):: x_cells, y_cells, tpoints, BC_Type(:)

  REAL*8,ALLOCATABLE,INTENT(OUT):: out_times(:)
  INTEGER,ALLOCATABLE,INTENT(OUT):: out_time_steps(:)
  INTEGER,INTENT(OUT):: n_out_times, restart_outlen, TEMP_out, GREY_E_out, GREY_F_out, E_ho_out
  INTEGER,INTENT(OUT):: GREY_kap_out, GREY_fsmall_out, MG_fsmall_out, res_history_out
  CHARACTER(100),INTENT(OUT):: outfile, restart_outfile, decomp_outfile, TEMP_outfile, GREY_E_outfile, E_ho_outfile
  CHARACTER(100),INTENT(OUT):: GREY_F_outfile, GREY_kap_outfile, GREY_fsmall_outfile, MG_fsmall_outfile, res_history_outfile

  REAL*8,ALLOCATABLE,INTENT(OUT):: nu_g(:)
  INTEGER,INTENT(OUT):: N_g

  INTEGER,INTENT(INOUT):: N_m
  REAL*8,ALLOCATABLE,INTENT(OUT):: Omega_x(:), Omega_y(:), quad_weight(:)

  REAL*8,INTENT(OUT):: sig_R, ar, cv

  !LOCAL VARIABLES
  INTEGER:: inpunit = 10
  INTEGER:: err

  OPEN(UNIT=inpunit,FILE="input/input.inp",STATUS='OLD',ACTION='READ',IOSTAT=err)
  !   making sure file exists/opens, if not tells user
  IF (err .NE. 0) THEN
      WRITE(*,'(A)') 'The file, input.inp, could not open properly.'
      STOP
  END IF

  CALL INPUT_RUN_STATE(inpunit,database_gen,database_add,run_type,restart_infile,use_grey)

  CALL INPUT_SOLVER_OPTS(inpunit,chi,conv_ho,conv_lo,conv_gr,comp_unit,line_src,maxit_RTE,maxit_MLOQD,maxit_GLOQD,&
    conv_type,threads,kapE_dT_flag,enrgy_strc,quadrature,Nwt_upbnd)

  sig_R=2d0*pi**5/(15d0*c**2*h**3*erg**4*comp_unit) !(erg/(ev**4 cm**2 sh))
  aR=4d0*sig_R/c

  CALL INPUT_PARAMETERS(inpunit,erg,xlen,ylen,x_cells,y_cells,tlen,delt,BC_Type,bcT_left,bcT_right,bcT_upper,&
    bcT_lower,Tini)

  ALLOCATE(Delx(x_cells))
  ALLOCATE(Dely(y_cells))
  Delx = xlen/REAL(x_cells,8)
  Dely = ylen/REAL(y_cells,8)
  tpoints = NINT(tlen/delt)
  Cv=0.5917d0*ar*(bcT_left)**3*((1.38d-16)**4*(11600d0**4)*erg**4)

  CALL INPUT_ENERGY(inpunit,enrgy_strc,N_g,nu_g)

  CALL INPUT_QUAD(quadrature,N_m,Omega_x,Omega_y,quad_weight)

  CALL INPUT_OUTPUT_OPTS(inpunit,out_times,out_time_steps,n_out_times,restart_outlen,TEMP_out,GREY_E_out,GREY_F_out,&
    GREY_kap_out,GREY_fsmall_out,MG_fsmall_out,res_history_out,outfile,restart_outfile,decomp_outfile,TEMP_outfile,&
    GREY_E_outfile,GREY_F_outfile,GREY_kap_outfile,GREY_fsmall_outfile,MG_fsmall_outfile,res_history_outfile,&
    E_ho_out,E_ho_outfile)

  CLOSE ( inpunit, STATUS='KEEP')

END SUBROUTINE INPUT

!============================================================================================================!
!
!============================================================================================================!
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

!============================================================================================================!
!
!============================================================================================================!

SUBROUTINE INPUT_RUN_STATE(inpunit,database_gen,database_add,run_type,restart_infile,use_grey)

  IMPLICIT NONE

  !INPUT VARIABLES
  INTEGER,INTENT(IN):: inpunit

  !OUTPUT VARIABLES
  INTEGER,INTENT(OUT):: database_gen, database_add, use_grey
  CHARACTER(100),INTENT(OUT):: run_type, restart_infile

  !LOCAL VARIABLES
  CHARACTER(1000):: line
  CHARACTER(100):: key
  CHARACTER(100):: block
  CHARACTER(100),DIMENSION(3):: args
  INTEGER:: io, io2, block_found

  !DEFAULT VALUES
  run_type = 'mlqd'
  restart_infile = ''
  database_gen = 0
  database_add = 0
  use_grey = 1

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

      ELSE IF (trim(key) .EQ. 'database_gen') THEN
        READ(args(1),*) database_gen

      ELSE IF (trim(key) .EQ. 'database_add') THEN
        READ(args(1),*) database_add

      ELSE IF (trim(key) .EQ. 'restart') THEN
        READ(args(1),*) restart_infile

      ELSE IF (trim(key) .EQ. 'use_grey') THEN
        READ(args(1),*) use_grey

      END IF

    END IF

  END DO

  IF ((run_type .NE. 'mlqd').AND.(database_gen .EQ. 1)) THEN
    WRITE(*,*) '*** WARNING ***'
    WRITE(*,*) 'database_gen = 1 IN INPUT DECK WHILE run_type = ', trim(run_type)
    WRITE(*,*) 'RESETTING database_gen = 0'
    database_gen = 0
    database_add = 0
  END IF

  IF (ANY(run_type .EQ. (/'rte_no_qd'/)) ) THEN
    use_grey = 0
  END IF

END SUBROUTINE INPUT_RUN_STATE

!============================================================================================================!
!
!============================================================================================================!

SUBROUTINE INPUT_SOLVER_OPTS(inpunit,chi,conv_ho,conv_lo,conv_gr,comp_unit,line_src,maxit_RTE,maxit_MLOQD,maxit_GLOQD,&
  conv_type,threads,kapE_dT_flag,enrgy_strc,quadrature,Nwt_upbnd)

  IMPLICIT NONE

  !INPUT VARIABLES
  INTEGER,INTENT(IN):: inpunit

  !OUTPUT VARIABLES
  REAL*8,INTENT(OUT):: chi, conv_ho, conv_lo, conv_gr, comp_unit, line_src, Nwt_upbnd
  INTEGER,INTENT(OUT):: maxit_RTE, maxit_MLOQD, maxit_GLOQD, conv_type, threads
  CHARACTER(100),INTENT(OUT):: kapE_dT_flag, enrgy_strc, quadrature

  !LOCAL VARIABLES
  CHARACTER(1000):: line
  CHARACTER(100):: key
  CHARACTER(100):: block
  CHARACTER(100),DIMENSION(3):: args
  INTEGER:: io, io2, block_found
  INTEGER:: conv_ho_flag, conv_lo_flag, conv_gr_flag

  !DEFAULT VALUES
  maxit_RTE = 25
  maxit_MLOQD = 500
  maxit_GLOQD = 100
  kapE_dT_flag = 'off'
  chi = 0.7d0
  conv_ho = 1d-8
  conv_lo = 1d-9
  conv_type = 1
  quadrature = 'abu36'
  threads = 1
  enrgy_strc = 'JCP'
  comp_unit = 1d13
  line_src = 5d0
  Nwt_upbnd = 1d10

  conv_ho_flag = 0
  conv_lo_flag = 0
  conv_gr_flag = 0

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
        READ(args(1),*) conv_gr
        conv_gr_flag = 1
        IF (conv_gr .LE. 0) THEN
          WRITE(*,*)
          WRITE(*,'(A)') '*** WARNING ***'
          WRITE(*,'(A)') '"conv_lo" cannot have a negative value'
          WRITE(*,'(A)') 'Setting conv_lo to default'
          conv_gr = 1d-10
        ELSE IF (conv_gr .GE. 1) THEN
          WRITE(*,*)
          WRITE(*,'(A)') '*** WARNING ***'
          WRITE(*,'(A)') '"conv_lo" cannot have a value greater than 1'
          WRITE(*,'(A)') 'Setting conv_lo to default'
          conv_gr = 1d-10
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
        READ(args(1),*) kapE_dT_flag
        IF ( ALL(kapE_dT_flag .NE. (/ 'on ' , 'off' /)) ) THEN
          WRITE(*,*)
          WRITE(*,'(A)') '*** WARNING ***'
          WRITE(*,'(A)') 'unknown flag detected for "kapE_dT" input'
          WRITE(*,'(A)') 'acceptable "kapE_dT" inputs are [on, off]'
          WRITE(*,'(A)') 'Setting kapE_dT_flag to default (off)'
          kapE_dT_flag = 'off'
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
        READ(args(1),*) line_src

      ELSE IF (trim(key) .EQ. 'Nwt_upbnd') THEN
        READ(args(1),*) Nwt_upbnd

      END IF

    END IF

  END DO


  IF (conv_ho_flag .EQ. 1) THEN
    IF (conv_lo_flag .EQ. 1) THEN
      IF (conv_gr_flag .EQ. 0) conv_gr = conv_lo/10d0
    ELSE
      conv_lo = conv_ho/10d0
      IF (conv_gr_flag .EQ. 0) conv_gr = conv_ho/100d0
    END IF
  ELSE
    IF (conv_lo_flag .EQ. 1) THEN
      conv_ho = conv_lo*10d0
      IF (conv_gr_flag .EQ. 0) conv_gr = conv_lo/10d0
    ELSE
      IF (conv_gr_flag .EQ. 1) THEN
        conv_ho = conv_gr*100d0
        conv_lo = conv_gr*10d0
      ELSE
        CONTINUE
      END IF
    END IF
  END IF

END SUBROUTINE INPUT_SOLVER_OPTS

!============================================================================================================!
!
!============================================================================================================!

SUBROUTINE INPUT_PARAMETERS(inpunit,erg,xlen,ylen,x_cells,y_cells,tlen,delt,BC_Type,bcT_left,bcT_right,bcT_upper,&
  bcT_lower,Tini)

  IMPLICIT NONE

  !INPUT VARIABLES
  REAL*8,INTENT(IN):: erg
  INTEGER,INTENT(IN):: inpunit

  !OUTPUT VARIABLES
  REAL*8,INTENT(OUT):: xlen, ylen, tlen, delt, bcT_left, bcT_right, bcT_upper, bcT_lower, Tini
  INTEGER,INTENT(OUT):: x_cells, y_cells, BC_Type(:)

  !LOCAL VARIABLES
  CHARACTER(1000):: line
  CHARACTER(100):: key, unit
  CHARACTER(100):: block
  CHARACTER(100),DIMENSION(4):: args
  INTEGER:: io, io2, block_found

  !DEFAULT VALUES
  xlen = 6d0
  ylen = 6d0
  x_cells = 10
  y_cells = 10
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
        READ(args(1),*) x_cells
        IF (x_cells .LE. 0) THEN
          WRITE(*,*)
          WRITE(*,*) 'x_cells Omega_xst be a positive integer (>0)'
          WRITE(*,*) 'resetting x_cells = 10'
          x_cells = 10
        END IF

      ELSE IF (trim(key) .EQ. 'y_cells') THEN
        READ(args(1),*) y_cells
        IF (y_cells .LE. 0) THEN
          WRITE(*,*)
          WRITE(*,*) 'y_cells Omega_xst be a positive integer (>0)'
          WRITE(*,*) 'resetting y_cells = 10'
          y_cells = 10
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

!============================================================================================================!
!
!============================================================================================================!

SUBROUTINE INPUT_OUTPUT_OPTS(inpunit,out_times,out_time_steps,n_out_times,restart_outlen,TEMP_out,GREY_E_out,GREY_F_out,&
  GREY_kap_out,GREY_fsmall_out,MG_fsmall_out,res_history_out,outfile,restart_outfile,decomp_outfile,TEMP_outfile,&
  GREY_E_outfile,GREY_F_outfile,GREY_kap_outfile,GREY_fsmall_outfile,MG_fsmall_outfile,res_history_outfile,&
  E_ho_out,E_ho_outfile)

  IMPLICIT NONE

  !OUTPUT VARIABLES
  REAL*8,ALLOCATABLE,INTENT(OUT):: out_times(:)
  INTEGER,ALLOCATABLE,INTENT(OUT):: out_time_steps(:)
  INTEGER,INTENT(OUT):: inpunit, n_out_times, restart_outlen, TEMP_out, GREY_E_out, GREY_F_out, E_ho_out
  INTEGER,INTENT(OUT):: GREY_kap_out, GREY_fsmall_out, MG_fsmall_out, res_history_out
  CHARACTER(100),INTENT(OUT):: outfile, restart_outfile, decomp_outfile, TEMP_outfile, GREY_E_outfile, E_ho_outfile
  CHARACTER(100),INTENT(OUT):: GREY_F_outfile, GREY_kap_outfile, GREY_fsmall_outfile, MG_fsmall_outfile, res_history_outfile

  !LOCAL VARIABLES
  CHARACTER(1000):: line
  CHARACTER(100):: key
  CHARACTER(100):: block
  CHARACTER(100),DIMENSION(3):: args
  INTEGER:: io, io2, block_found

  !DEFAULT VALUES
  n_out_times = 0
  restart_outlen = 0
  TEMP_out = 0
  GREY_E_out = 0
  E_ho_out = 0
  GREY_F_out = 0
  GREY_kap_out = 0
  GREY_fsmall_out = 0
  MG_fsmall_out = 0
  res_history_out = 0
  outfile = 'output.out'
  restart_outfile = 'restart.out'
  decomp_outfile = 'MG_QD_decomp.out'
  TEMP_outfile = 'TEMP.out'
  GREY_E_outfile = 'GREY_E.out'
  E_ho_outfile = 'E_ho.out'
  GREY_F_outfile = 'GREY_F.out'
  GREY_kap_outfile = 'GREY_kaps.out'
  GREY_fsmall_outfile = 'QD_factors.out'
  MG_fsmall_outfile = 'MG_QD_factors.out'
  res_history_outfile = 'iteration_history.out'

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

      ELSE IF (trim(key) .EQ. 'restart_outlen') THEN
        READ(args(1),*) restart_outlen
        IF (TEMP_out .LT. 0) THEN
          WRITE(*,*)
          WRITE(*,'(A)') '*** WARNING ***'
          WRITE(*,'(A)') 'Wrong value set for "TEMP_out" input, Omega_xst be >= 0'
          WRITE(*,'(A)') 'Setting restart_outlen to default'
          restart_outlen = 0
        END IF

        ELSE IF (trim(key) .EQ. 'TEMP_out') THEN
          READ(args(1),*) TEMP_out
          IF ((TEMP_out .NE. 0) .AND. (TEMP_out .NE. 1)) THEN
            WRITE(*,*)
            WRITE(*,'(A)') '*** WARNING ***'
            WRITE(*,'(A)') 'Wrong value set for "TEMP_out" input, Omega_xst be contained in [0,1]'
            WRITE(*,'(A)') 'Setting TEMP_out to default'
            TEMP_out = 0
          END IF

        ELSE IF (trim(key) .EQ. 'GREY_E_out') THEN
          READ(args(1),*) GREY_E_out
          IF ((GREY_E_out .NE. 0) .AND. (GREY_E_out .NE. 1)) THEN
            WRITE(*,*)
            WRITE(*,'(A)') '*** WARNING ***'
            WRITE(*,'(A)') 'Wrong value set for "GREY_E_out" input, Omega_xst be contained in [0,1]'
            WRITE(*,'(A)') 'Setting GREY_E_out to default'
            GREY_E_out = 0
          END IF

        ELSE IF (trim(key) .EQ. 'E_ho_out') THEN
          READ(args(1),*) E_ho_out
          IF ((E_ho_out .NE. 0) .AND. (E_ho_out .NE. 1)) THEN
            WRITE(*,*)
            WRITE(*,'(A)') '*** WARNING ***'
            WRITE(*,'(A)') 'Wrong value set for "E_ho_out" input, Omega_xst be contained in [0,1]'
            WRITE(*,'(A)') 'Setting E_ho_out to default'
            E_ho_out = 0
          END IF

        ELSE IF (trim(key) .EQ. 'GREY_F_out') THEN
          READ(args(1),*) GREY_F_out
          IF ((GREY_F_out .NE. 0) .AND. (GREY_F_out .NE. 1)) THEN
            WRITE(*,*)
            WRITE(*,'(A)') '*** WARNING ***'
            WRITE(*,'(A)') 'Wrong value set for "GREY_F_out" input, Omega_xst be contained in [0,1]'
            WRITE(*,'(A)') 'Setting GREY_F_out to default'
            GREY_F_out = 0
          END IF

        ELSE IF (trim(key) .EQ. 'GREY_kap_out') THEN
          READ(args(1),*) GREY_kap_out
          IF ((GREY_kap_out .NE. 0) .AND. (GREY_kap_out .NE. 1)) THEN
            WRITE(*,*)
            WRITE(*,'(A)') '*** WARNING ***'
            WRITE(*,'(A)') 'Wrong value set for "GREY_kap_out" input, Omega_xst be contained in [0,1]'
            WRITE(*,'(A)') 'Setting GREY_kap_out to default'
            GREY_kap_out = 0
          END IF

        ELSE IF (trim(key) .EQ. 'GREY_fsmall_out') THEN
          READ(args(1),*) GREY_fsmall_out
          IF ((GREY_fsmall_out .NE. 0) .AND. (GREY_fsmall_out .NE. 1)) THEN
            WRITE(*,*)
            WRITE(*,'(A)') '*** WARNING ***'
            WRITE(*,'(A)') 'Wrong value set for "GREY_fsmall_out" input, Omega_xst be contained in [0,1]'
            WRITE(*,'(A)') 'Setting GREY_fsmall_out to default'
            GREY_fsmall_out = 0
          END IF

        ELSE IF (trim(key) .EQ. 'MG_fsmall_out') THEN
          READ(args(1),*) MG_fsmall_out
          IF ((MG_fsmall_out .NE. 0) .AND. (MG_fsmall_out .NE. 1)) THEN
            WRITE(*,*)
            WRITE(*,'(A)') '*** WARNING ***'
            WRITE(*,'(A)') 'Wrong value set for "MG_fsmall_out" input, Omega_xst be contained in [0,1]'
            WRITE(*,'(A)') 'Setting MG_fsmall_out to default'
            MG_fsmall_out = 0
          END IF

        ELSE IF (trim(key) .EQ. 'res_history_out') THEN
          READ(args(1),*) res_history_out
          IF ((res_history_out .NE. 0) .AND. (res_history_out .NE. 1)) THEN
            WRITE(*,*)
            WRITE(*,'(A)') '*** WARNING ***'
            WRITE(*,'(A)') 'Wrong value set for "res_history_out" input, Omega_xst be contained in [0,1]'
            WRITE(*,'(A)') 'Setting res_history_out to default'
            res_history_out = 0
          END IF

        ELSE IF (trim(key) .EQ. 'outfile') THEN
          READ(args(1),*) outfile

        ELSE IF (trim(key) .EQ. 'restart_outfile') THEN
          READ(args(1),*) restart_outfile

        ELSE IF (trim(key) .EQ. 'decomp_outfile') THEN
          READ(args(1),*) decomp_outfile

        ELSE IF (trim(key) .EQ. 'TEMP_outfile') THEN
          READ(args(1),*) TEMP_outfile

        ELSE IF (trim(key) .EQ. 'GREY_E_outfile') THEN
          READ(args(1),*) GREY_E_outfile

        ELSE IF (trim(key) .EQ. 'E_ho_outfile') THEN
          READ(args(1),*) E_ho_outfile

        ELSE IF (trim(key) .EQ. 'GREY_F_outfile') THEN
          READ(args(1),*) GREY_F_outfile

        ELSE IF (trim(key) .EQ. 'GREY_kap_outfile') THEN
          READ(args(1),*) GREY_kap_outfile

        ELSE IF (trim(key) .EQ. 'GREY_fsmall_outfile') THEN
          READ(args(1),*) GREY_fsmall_outfile

        ELSE IF (trim(key) .EQ. 'MG_fsmall_outfile') THEN
          READ(args(1),*) MG_fsmall_outfile

        ELSE IF (trim(key) .EQ. 'res_history_outfile') THEN
          READ(args(1),*) res_history_outfile

        ELSE IF (trim(key) .EQ. 'out_times') THEN
          READ(args(1),*) n_out_times
          IF ((n_out_times .LE. 0)) THEN
            WRITE(*,*)
            WRITE(*,'(A)') '*** WARNING ***'
            WRITE(*,'(A)') 'You cannot have a zero or less out times, n_out_times </= 0'
            WRITE(*,'(A)') 'Aborting run'
            STOP
          END IF
          ALLOCATE(out_times(n_out_times))
          ALLOCATE(out_time_steps(n_out_times))
          READ(line,*,IOSTAT=io2) key, n_out_times, out_times(:)

      END IF

    END IF

  END DO

  IF (((TEMP_out .EQ. 1).OR.(GREY_E_out .EQ. 1).OR.(GREY_F_out .EQ. 1).OR.(GREY_kap_out .EQ. 1).OR.(GREY_fsmall_out .EQ. 1))&
    &.AND.(n_out_times .EQ. 0)) STOP 'No out times were specified, aborting program'

END SUBROUTINE INPUT_OUTPUT_OPTS

!============================================================================================================!
!
!============================================================================================================!

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

!============================================================================================================!
!
!============================================================================================================!

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

END MODULE INPUTS
