MODULE OUTPUTS

  IMPLICIT NONE

CONTAINS

!==================================================================================================================================!
!
!==================================================================================================================================!
SUBROUTINE OUTFILE_HEADER(f_unit,datevalues,database_gen,database_add,run_type,restart_infile,use_grey,chi,&
  conv_ho,conv_lo,conv_gr,comp_unit,line_src,maxit_RTE,maxit_MLOQD,maxit_GLOQD,conv_type,N_m,threads,&
  kapE_dT_flag,enrgy_strc,Nwt_upbnd,xlen,ylen,N_x,N_y,tlen,delt,bc_type,bcT_left,bcT_right,bcT_upper,&
  bcT_lower,Tini,out_times,out_time_steps,n_out_times,restart_outlen,TEMP_out,GREY_E_out,GREY_F_out,&
  GREY_kap_out,GREY_fsmall_out,MG_fsmall_out,res_history_out,outfile,restart_outfile,decomp_outfile,TEMP_outfile,&
  GREY_E_outfile,GREY_F_outfile,GREY_kap_outfile,GREY_fsmall_outfile,MG_fsmall_outfile,res_history_outfile,&
  HO_E_out,HO_E_outfile,delx,dely,N_g,nu_g,N_t,sig_r,ar,sig0,pi,c,h,cv,Omega_x,Omega_y,quad_weight,quadrature)

  IMPLICIT NONE

  !INPUT VARIABLES
  INTEGER,INTENT(IN):: f_unit, datevalues(:)
  REAL*8,INTENT(IN):: chi, conv_ho, conv_lo, conv_gr, line_src, Nwt_upbnd, comp_unit
  REAL*8,INTENT(IN):: xlen, ylen, tlen, delt, bcT_left, bcT_right, bcT_upper, bcT_lower, Tini
  REAL*8,INTENT(IN):: Delx(:), Dely(:)
  REAL*8,INTENT(IN):: sig_r, ar, sig0, pi, c, h, cv
  REAL*8,INTENT(IN):: out_times(:), nu_g(:), Omega_x(:), Omega_y(:), quad_weight(:)
  INTEGER,INTENT(IN):: database_gen, database_add, use_grey, bc_type(:)
  INTEGER,INTENT(IN):: maxit_RTE, maxit_MLOQD, maxit_GLOQD, conv_type, N_m, threads, N_g, N_t
  INTEGER,INTENT(IN):: N_x, N_y
  INTEGER,INTENT(IN):: n_out_times, restart_outlen, TEMP_out, GREY_E_out, GREY_F_out, HO_E_out
  INTEGER,INTENT(IN):: GREY_kap_out, GREY_fsmall_out, MG_fsmall_out, res_history_out
  INTEGER,INTENT(IN):: out_time_steps(:)
  CHARACTER(100),INTENT(IN):: outfile, restart_outfile, decomp_outfile, TEMP_outfile, GREY_E_outfile, HO_E_outfile
  CHARACTER(100),INTENT(IN):: GREY_F_outfile, GREY_kap_outfile, GREY_fsmall_outfile, MG_fsmall_outfile, res_history_outfile
  CHARACTER(100),INTENT(IN):: run_type, restart_infile
  CHARACTER(100),INTENT(IN):: kapE_dT_flag, enrgy_strc, quadrature

  !LOCAL VARIABLES
  INTEGER:: g, m

  !CODE VERSION AND RUN DATE DATA
  WRITE(f_unit,'(A)') '-----------------------------------------------------------'
  WRITE(f_unit,'(A)') 'JOSEPH COALE 2D DISSERTATION CODE'
  WRITE(f_unit,'(2(A,I2.2),A,I4.4)') 'RUN DATE = ',datevalues(2),'/',datevalues(3),'/',datevalues(1)
  WRITE(f_unit,'(3(A,I2.2))') 'RUN TIME = ',datevalues(5),':',datevalues(6),':',datevalues(7)

  !INPUT ECHO
  WRITE(f_unit,'(A)') '-----------------------------------------------------------'
  WRITE(f_unit,'(A)') 'INPUT ECHO:'
  WRITE(f_unit,*)
  WRITE(f_unit,'(A)') '[RUN_STATE]'
  WRITE(f_unit,'(2A)')   '   run_type       ',run_type
  WRITE(f_unit,'(A,I1)') '   database_gen   ',database_gen
  WRITE(f_unit,'(A,I1)') '   database_add   ',database_add
  WRITE(f_unit,'(A,I1)') '   use_grey       ',use_grey
  IF (restart_infile .NE. '') WRITE(f_unit,'(2A)') '   restart_infile ',restart_infile
  WRITE(f_unit,*)

  WRITE(f_unit,'(A)') '[SOLVER_OPTS]'
  WRITE(f_unit,'(A,I3.3)')  '   maxit_RTE     ',maxit_RTE
  WRITE(f_unit,'(A,I4.4)')  '   maxit_MLOQD   ',maxit_MLOQD
  WRITE(f_unit,'(A,I4.4)')  '   maxit_GLOQD   ',maxit_GLOQD
  WRITE(f_unit,'(2A)')      '   kapE_dT       ',kapE_dT_flag
  WRITE(f_unit,'(A,ES7.1)') '   chi           ',chi
  WRITE(f_unit,'(A,ES7.1)') '   conv_ho       ',conv_ho
  WRITE(f_unit,'(A,ES7.1)') '   conv_lo       ',conv_lo
  WRITE(f_unit,'(A,ES7.1)') '   conv_gr       ',conv_gr
  WRITE(f_unit,'(A,I1)')    '   conv_type     ',conv_type
  WRITE(f_unit,'(2A)')      '   quad          ',quadrature
  WRITE(f_unit,'(A,I2.2)')  '   threads       ',threads
  WRITE(f_unit,'(2A)')      '   enrgy_strc    ',enrgy_strc
  WRITE(f_unit,'(A,ES7.1)') '   comp_unit     ',comp_unit
  WRITE(f_unit,'(A,ES7.1)') '   line_src      ',line_src
  WRITE(f_unit,*)

  WRITE(f_unit,'(A)') '[PARAMETERS]'
  WRITE(f_unit,'(A,ES9.3,A)') '   xlen      ',xlen,' (cm)'
  WRITE(f_unit,'(A,ES9.3,A)') '   ylen      ',ylen,' (cm)'
  WRITE(f_unit,'(A,I5)')      '   N_x       ',N_x
  WRITE(f_unit,'(A,I5)')      '   N_y       ',N_y
  WRITE(f_unit,'(A,ES9.3,A)') '   tlen      ',tlen,' (sh)'
  WRITE(f_unit,'(A,ES9.3,A)') '   delt      ',delt,' (sh)'
  WRITE(f_unit,'(A,4I2)')     '   bc_type   ',bc_type
  WRITE(f_unit,'(A,ES9.3,A)') '   bcT_left  ',bcT_left,' (ev)'
  WRITE(f_unit,'(A,ES9.3,A)') '   bcT_upper ',bcT_upper,' (ev)'
  WRITE(f_unit,'(A,ES9.3,A)') '   bcT_right ',bcT_right,' (ev)'
  WRITE(f_unit,'(A,ES9.3,A)') '   bcT_lower ',bcT_lower,' (ev)'
  WRITE(f_unit,'(A,ES9.3,A)') '   Tini      ',Tini,' (ev)'
  WRITE(f_unit,*)

  WRITE(f_unit,'(A)') '[OUTPUT_OPTS]'
  WRITE(f_unit,'(A,I1)') '   TEMP_out              ',TEMP_out
  WRITE(f_unit,'(A,I1)') '   GREY_E_out            ',GREY_E_out
  WRITE(f_unit,'(A,I1)') '   HO_E_out              ',HO_E_out
  WRITE(f_unit,'(A,I1)') '   GREY_F_out            ',GREY_F_out
  WRITE(f_unit,'(A,I1)') '   GREY_kap_out          ',GREY_kap_out
  WRITE(f_unit,'(A,I1)') '   GREY_fsmall_out       ',GREY_fsmall_out
  WRITE(f_unit,'(A,I1)') '   MG_fsmall_out         ',MG_fsmall_out
  WRITE(f_unit,'(A,I1)') '   res_history_out       ',res_history_out
  IF (restart_outlen .GT. 0) THEN
     WRITE(f_unit,'(A,I1)') '   restart_outlen        ',restart_outlen
     WRITE(f_unit,'(2A)')   '   restart_outfile       ',restart_outfile
   END IF
  WRITE(f_unit,'(2A)')   '   outfile               ',outfile
  IF (database_gen .EQ. 1)   WRITE(f_unit,'(2A)') '   decomp_outfile        ',decomp_outfile
  IF (TEMP_out .EQ. 1)     WRITE(f_unit,'(2A)') '   TEMP_outfile          ',TEMP_outfile
  IF (GREY_E_out .EQ. 1)   WRITE(f_unit,'(2A)') '   GREY_E_outfile        ',GREY_E_outfile
  IF (HO_E_out .EQ. 1)   WRITE(f_unit,'(2A)')   '   HO_E_outfile          ',HO_E_outfile
  IF (GREY_F_out .EQ. 1)   WRITE(f_unit,'(2A)') '   GREY_F_outfile        ',GREY_F_outfile
  IF (GREY_kap_out .EQ. 1) WRITE(f_unit,'(2A)') '   GREY_kap_outfile      ',GREY_kap_outfile
  IF (GREY_fsmall_out .EQ. 1) WRITE(f_unit,'(2A)') '   GREY_fsmall_outfile   ',GREY_fsmall_outfile
  IF (MG_fsmall_out .EQ. 1)   WRITE(f_unit,'(2A)') '   MG_fsmall_outfile     ',MG_fsmall_outfile
  IF (res_history_out .EQ. 1)   WRITE(f_unit,'(2A)') '   res_history_outfile   ',res_history_outfile
  WRITE(f_unit,*)

  !DISCRETIZATION DATA (PROBLEM DOMAIN AND ENERGY GROUPS)
  WRITE(f_unit,'(A)') '-----------------------------------------------------------'
  WRITE(f_unit,'(A)') 'DISCRETIZATION:'
  WRITE(f_unit,*)
  WRITE(f_unit,'(A,I5.5)') 'time_steps ',N_t
  WRITE(f_unit,'(A,I5.5)') 'N_x ',N_x
  WRITE(f_unit,'(A,I5.5)') 'N_y ',N_y
  WRITE(f_unit,'(A,ES9.3,A)') 'delx    ',Delx(1),' (cm)'
  WRITE(f_unit,'(A,ES9.3,A)') 'dely    ',Dely(1),' (cm)'
  WRITE(f_unit,'(A,I5.5)') 'N_g  ',N_g
  WRITE(f_unit,*)
  WRITE(f_unit,'(A)') 'group     energy interval (ev)'
  DO g=1,N_g
    WRITE(f_unit,'(I2.2,A,ES12.5,A,ES12.5,A)') g,'  [',nu_g(g),' -> ',nu_g(g+1),' ]'
  END DO
  WRITE(f_unit,*)
  WRITE(f_unit,'(A)') 'quadrature     angle (Omega_x, Omega_y, weight)'
  DO m=1,N_m
    WRITE(f_unit,'(I3.3,A,ES12.5,A,ES12.5,A,ES12.5,A)') m,'  [',Omega_x(m),' , ',Omega_y(m),' , ',quad_weight(m),' ]'
  END DO
  WRITE(f_unit,*)
  WRITE(f_unit,'(A,ES12.5)') 'SUM(Omega_x(:)) ',SUM(Omega_x(:))
  WRITE(f_unit,'(A,ES12.5)') 'SUM(Omega_y(:)) ',SUM(Omega_y(:))
  WRITE(f_unit,'(A,ES12.5)') 'SUM(quad_weight(:)) ',SUM(quad_weight(:))

  !CONSTANT VALUES ECHO
  WRITE(f_unit,'(A)') '-----------------------------------------------------------'
  WRITE(f_unit,'(A)') 'CONSTANTS:'
  WRITE(f_unit,*)
  WRITE(f_unit,'(A,ES22.15)') 'pi =    ',pi
  WRITE(f_unit,'(A,ES22.15,A)') 'c =     ',c,' (cm/sh)'
  WRITE(f_unit,'(A,ES22.15,A)') 'h =     ',h,' (erg*sh)'
  WRITE(f_unit,'(A,ES22.15,A)') 'Sig_R = ',sig_R*comp_unit,' (erg/(ev**4 cm**2 sh))'
  WRITE(f_unit,'(A,ES22.15,A)') 'a_R =   ',aR*comp_unit,' (erg/(ev**4 cm**3))'
  WRITE(f_unit,'(A,ES22.15,A)') 'Sig_0 = ',sig0,' (ev**3/cm)'
  WRITE(f_unit,'(A,ES22.15,A)') 'Cv =    ',cv*comp_unit,' (erg/(ev cm**3))'
  WRITE(f_unit,'(A,ES22.15,A)') 'comp_unit = ',comp_unit,' (erg/unit)'

END SUBROUTINE OUTFILE_HEADER

!==================================================================================================================================!
!
!==================================================================================================================================!
SUBROUTINE ITERATIVE_OUTPUT(f_unit,Tstep,TR_its,MGQD_its,GQD_its,Time,c,TR_Tnorm,TR_Enorm,TR_Trho,TR_Erho,MGQD_Tnorm,&
  MGQD_Enorm,MGQD_Trho,MGQD_Erho,TR_Residual,MGQD_Residual)

  INTEGER,INTENT(IN):: f_unit, Tstep, TR_its, MGQD_its(:), GQD_its
  REAL*8,INTENT(IN):: Time, c
  REAL*8,INTENT(IN):: TR_Residual(:,:,:,:,:,:), MGQD_Residual(:,:,:,:,:,:,:)
  REAL*8,INTENT(IN):: TR_Tnorm(:), TR_Enorm(:), TR_Trho(:), TR_Erho(:)
  REAL*8,INTENT(IN):: MGQD_Tnorm(:,:), MGQD_Enorm(:,:), MGQD_Trho(:,:), MGQD_Erho(:,:)
  INTEGER:: ti, mi, N_m, N_g, g, m1, m2, m

  N_m = SIZE(TR_Residual,3)
  N_g = SIZE(TR_Residual,4)

  WRITE(f_unit,'(A)') '*========================='
  WRITE(f_unit,'(A,I4)')    'TIME STEP ',Tstep
  WRITE(f_unit,'(A,ES9.3)') 'TIME      ',Time
  WRITE(f_unit,'(A,ES9.3)') 'ct        ',Time*c
  WRITE(f_unit,'(A,I4)')    'TR_its    ',TR_its
  WRITE(f_unit,'(A,I4)')    'MGQD_its  ',SUM(MGQD_its)
  WRITE(f_unit,'(A,I4)')    'GQD_its   ',GQD_its

  WRITE(f_unit,'(A)') '*'
  WRITE(f_unit,'(A)') 'Transport Iteration history'
  WRITE(f_unit,'(A)') 'TRit  Enorm     Erho      Tnorm     Trho'
  DO ti=1,TR_its
    WRITE(f_unit,'(I4.4,4ES10.3)') ti, TR_Tnorm(ti), TR_Trho(ti), TR_Enorm(ti), TR_Erho(ti)
  END DO

  WRITE(f_unit,'(A)') '*'
  WRITE(f_unit,'(A)') 'Transport Residuals (Summary)'
  WRITE(f_unit,'(A)') 'TRit  quad1     quad2     quad3     quad4'
  DO ti=1,TR_its
    WRITE(f_unit,'(I4.4,4ES10.3)') ti, MAXVAL(ABS(TR_Residual(:,:,1:N_m/4,:,1,ti))),&
                                       MAXVAL(ABS(TR_Residual(:,:,N_m/4+1:N_m/2,:,1,ti))),&
                                       MAXVAL(ABS(TR_Residual(:,:,N_m/2+1:3*N_m/4,:,1,ti))),&
                                       MAXVAL(ABS(TR_Residual(:,:,3*N_m/4+1:N_m,:,1,ti)))
  END DO

  WRITE(f_unit,'(A)') '*'
  WRITE(f_unit,'(A)') 'Transport Residuals (In Depth)'
  WRITE(f_unit,'(A)',ADVANCE='NO') 'TRit  quad    g01'
  DO g=2,N_g
    WRITE(f_unit,'(A,I2.2)',ADVANCE='NO') '       g',g
  END DO
  WRITE(f_unit,*)
  DO ti=1,TR_its
    m1=1
    m2=N_m/4
    DO m=1,4
      WRITE(f_unit,'(I4.4,I4,A)',ADVANCE='NO') ti, m,'   '
      DO g=1,N_g
        WRITE(f_unit,'(ES10.3)',ADVANCE='NO') MAXVAL(ABS(TR_Residual(:,:,m1:m2,g,1,ti)))
      END DO
      WRITE(f_unit,*)
      m1=m1+N_m/4
      m2=m2+N_m/4
    END DO
  END DO

  WRITE(f_unit,'(A)') '*'
  WRITE(f_unit,'(A)') 'MGQD Iteration history'
  WRITE(f_unit,'(A)') 'TRit MGQDit   Enorm     Erho      Tnorm     Trho'
  DO ti=1,TR_its
    DO mi=1,MGQD_its(ti)
      WRITE(f_unit,'(I4.4,I6.4,A,4ES10.3)') ti, mi, '  ', MGQD_Tnorm(ti,mi), MGQD_Trho(ti,mi),&
        MGQD_Enorm(ti,mi), MGQD_Erho(ti,mi)
    END DO
  END DO

  WRITE(f_unit,'(A)') '*'
  WRITE(f_unit,'(A)') 'MGQD Residuals (In Depth)'
  WRITE(f_unit,'(A)',ADVANCE='NO') 'TRit MGQDit  EQ      g01'
  DO g=2,N_g
    WRITE(f_unit,'(A,I2.2)',ADVANCE='NO') '       g',g
  END DO
  WRITE(f_unit,*)
  DO ti=1,TR_its
    DO mi=1,MGQD_its(ti)

      WRITE(f_unit,'(I4.4,I4,A)',ADVANCE='NO') ti, mi,'     EB   '
      DO g=1,N_g
        WRITE(f_unit,'(ES10.3)',ADVANCE='NO') MAXVAL(ABS(MGQD_Residual(:,:,g,1,1,mi,ti)))
      END DO
      WRITE(f_unit,*)

      WRITE(f_unit,'(I4.4,I4,A)',ADVANCE='NO') ti, mi,'     MBxL '
      DO g=1,N_g
        WRITE(f_unit,'(ES10.3)',ADVANCE='NO') MAXVAL(ABS(MGQD_Residual(:,:,g,2,1,mi,ti)))
      END DO
      WRITE(f_unit,*)

      WRITE(f_unit,'(I4.4,I4,A)',ADVANCE='NO') ti, mi,'     MBxR '
      DO g=1,N_g
        WRITE(f_unit,'(ES10.3)',ADVANCE='NO') MAXVAL(ABS(MGQD_Residual(:,:,g,3,1,mi,ti)))
      END DO
      WRITE(f_unit,*)

      WRITE(f_unit,'(I4.4,I4,A)',ADVANCE='NO') ti, mi,'     MByB '
      DO g=1,N_g
        WRITE(f_unit,'(ES10.3)',ADVANCE='NO') MAXVAL(ABS(MGQD_Residual(:,:,g,4,1,mi,ti)))
      END DO
      WRITE(f_unit,*)

      WRITE(f_unit,'(I4.4,I4,A)',ADVANCE='NO') ti, mi,'     MByT '
      DO g=1,N_g
        WRITE(f_unit,'(ES10.3)',ADVANCE='NO') MAXVAL(ABS(MGQD_Residual(:,:,g,5,1,mi,ti)))
      END DO
      WRITE(f_unit,*)
    END DO
  END DO
  WRITE(f_unit,'(A)') '*========================='

END SUBROUTINE ITERATIVE_OUTPUT

!==================================================================================================================================!
!
!==================================================================================================================================!
SUBROUTINE SPECIFIC_OUTPUTS(Temp_Times,HO_E_avg_Times,GREY_E_avg_Times,&
  Times,f_unit,datevalues,database_gen,database_add,run_type,restart_infile,&
  use_grey,chi,conv_ho,conv_lo,conv_gr,comp_unit,line_src,maxit_RTE,maxit_MLOQD,maxit_GLOQD,conv_type,N_m,threads,&
  kapE_dT_flag,enrgy_strc,Nwt_upbnd,xlen,ylen,N_x,N_y,tlen,delt,bc_type,bcT_left,bcT_right,bcT_upper,bcT_lower,Tini,&
  out_times,out_time_steps,n_out_times,restart_outlen,TEMP_out,GREY_E_out,GREY_F_out,GREY_kap_out,GREY_fsmall_out,&
  MG_fsmall_out,res_history_out,outfile,restart_outfile,decomp_outfile,TEMP_outfile,GREY_E_outfile,GREY_F_outfile,&
  GREY_kap_outfile,GREY_fsmall_outfile,MG_fsmall_outfile,res_history_outfile,HO_E_out,HO_E_outfile,delx,dely,N_g,&
  nu_g,N_t,sig_r,ar,sig0,pi,c,h,cv,Omega_x,Omega_y,quad_weight,quadrature,MGQD_E_avg_Times)

  REAL*8,INTENT(IN):: Temp_Times(:,:,:), HO_E_avg_Times(:,:,:), GREY_E_avg_Times(:,:,:), MGQD_E_avg_Times(:,:,:)
  REAL*8,INTENT(IN):: Times(:)

  !INPUT VARIABLES
  INTEGER,INTENT(IN):: f_unit, datevalues(:)
  REAL*8,INTENT(IN):: chi, conv_ho, conv_lo, conv_gr, line_src, Nwt_upbnd, comp_unit
  REAL*8,INTENT(IN):: xlen, ylen, tlen, delt, bcT_left, bcT_right, bcT_upper, bcT_lower, Tini
  REAL*8,INTENT(IN):: Delx(:), Dely(:)
  REAL*8,INTENT(IN):: sig_r, ar, sig0, pi, c, h, cv
  REAL*8,INTENT(IN):: out_times(:), nu_g(:), Omega_x(:), Omega_y(:), quad_weight(:)
  INTEGER,INTENT(IN):: database_gen, database_add, use_grey, bc_type(:)
  INTEGER,INTENT(IN):: maxit_RTE, maxit_MLOQD, maxit_GLOQD, conv_type, N_m, threads, N_g, N_t
  INTEGER,INTENT(IN):: N_x, N_y
  INTEGER,INTENT(IN):: n_out_times, restart_outlen, TEMP_out, GREY_E_out, GREY_F_out, HO_E_out
  INTEGER,INTENT(IN):: GREY_kap_out, GREY_fsmall_out, MG_fsmall_out, res_history_out
  INTEGER,INTENT(IN):: out_time_steps(:)
  CHARACTER(100),INTENT(IN):: outfile, restart_outfile, decomp_outfile, TEMP_outfile, GREY_E_outfile, HO_E_outfile
  CHARACTER(100),INTENT(IN):: GREY_F_outfile, GREY_kap_outfile, GREY_fsmall_outfile, MG_fsmall_outfile, res_history_outfile
  CHARACTER(100),INTENT(IN):: run_type, restart_infile
  CHARACTER(100),INTENT(IN):: kapE_dT_flag, enrgy_strc, quadrature

  CHARACTER(100):: Data_Name, file

  IF (TEMP_out .EQ. 1) THEN
    Data_Name = 'Temperature'
    CALL SPCF_OUTPUT(Temp_Times,Times,TEMP_outfile,Data_Name,f_unit,datevalues,database_gen,database_add,run_type,restart_infile,&
      use_grey,chi,conv_ho,conv_lo,conv_gr,comp_unit,line_src,maxit_RTE,maxit_MLOQD,maxit_GLOQD,conv_type,N_m,threads,&
      kapE_dT_flag,enrgy_strc,Nwt_upbnd,xlen,ylen,N_x,N_y,tlen,delt,bcT_left,bcT_right,bcT_upper,bcT_lower,Tini,&
      out_times,out_time_steps,n_out_times,restart_outlen,TEMP_out,GREY_E_out,GREY_F_out,GREY_kap_out,GREY_fsmall_out,&
      MG_fsmall_out,res_history_out,outfile,restart_outfile,decomp_outfile,TEMP_outfile,GREY_E_outfile,GREY_F_outfile,&
      GREY_kap_outfile,GREY_fsmall_outfile,MG_fsmall_outfile,res_history_outfile,HO_E_out,HO_E_outfile,delx,dely,N_g,&
      nu_g,N_t,sig_r,ar,sig0,pi,c,h,cv,Omega_x,Omega_y,quad_weight,quadrature,bc_type)
  END IF

  IF (HO_E_out .EQ. 1) THEN
    Data_Name = 'High-Order Radiation Energy Density'
    CALL SPCF_OUTPUT(HO_E_avg_Times,Times,HO_E_outfile,Data_Name,f_unit,datevalues,database_gen,database_add,run_type,&
      restart_infile,use_grey,chi,conv_ho,conv_lo,conv_gr,comp_unit,line_src,maxit_RTE,maxit_MLOQD,maxit_GLOQD,conv_type,N_m,&
      threads,kapE_dT_flag,enrgy_strc,Nwt_upbnd,xlen,ylen,N_x,N_y,tlen,delt,bcT_left,bcT_right,bcT_upper,bcT_lower,Tini,&
      out_times,out_time_steps,n_out_times,restart_outlen,TEMP_out,GREY_E_out,GREY_F_out,GREY_kap_out,GREY_fsmall_out,&
      MG_fsmall_out,res_history_out,outfile,restart_outfile,decomp_outfile,TEMP_outfile,GREY_E_outfile,GREY_F_outfile,&
      GREY_kap_outfile,GREY_fsmall_outfile,MG_fsmall_outfile,res_history_outfile,HO_E_out,HO_E_outfile,delx,dely,N_g,&
      nu_g,N_t,sig_r,ar,sig0,pi,c,h,cv,Omega_x,Omega_y,quad_weight,quadrature,bc_type)
  END IF

  IF (HO_E_out .EQ. 1) THEN
    Data_Name = 'MGQD Radiation Energy Density'
    file = 'E_MGQD.out'
    CALL SPCF_OUTPUT(MGQD_E_avg_Times,Times,file,Data_Name,f_unit,datevalues,database_gen,database_add,run_type,&
      restart_infile,use_grey,chi,conv_ho,conv_lo,conv_gr,comp_unit,line_src,maxit_RTE,maxit_MLOQD,maxit_GLOQD,conv_type,N_m,&
      threads,kapE_dT_flag,enrgy_strc,Nwt_upbnd,xlen,ylen,N_x,N_y,tlen,delt,bcT_left,bcT_right,bcT_upper,bcT_lower,Tini,&
      out_times,out_time_steps,n_out_times,restart_outlen,TEMP_out,GREY_E_out,GREY_F_out,GREY_kap_out,GREY_fsmall_out,&
      MG_fsmall_out,res_history_out,outfile,restart_outfile,decomp_outfile,TEMP_outfile,GREY_E_outfile,GREY_F_outfile,&
      GREY_kap_outfile,GREY_fsmall_outfile,MG_fsmall_outfile,res_history_outfile,HO_E_out,HO_E_outfile,delx,dely,N_g,&
      nu_g,N_t,sig_r,ar,sig0,pi,c,h,cv,Omega_x,Omega_y,quad_weight,quadrature,bc_type)
  END IF

  IF (GREY_E_out .EQ. 1) THEN
    Data_Name = 'Grey Radiation Energy Density'
    CALL SPCF_OUTPUT(GREY_E_avg_Times,Times,GREY_E_outfile,Data_Name,f_unit,datevalues,database_gen,database_add,run_type,&
      restart_infile,use_grey,chi,conv_ho,conv_lo,conv_gr,comp_unit,line_src,maxit_RTE,maxit_MLOQD,maxit_GLOQD,conv_type,N_m,&
      threads,kapE_dT_flag,enrgy_strc,Nwt_upbnd,xlen,ylen,N_x,N_y,tlen,delt,bcT_left,bcT_right,bcT_upper,bcT_lower,Tini,&
      out_times,out_time_steps,n_out_times,restart_outlen,TEMP_out,GREY_E_out,GREY_F_out,GREY_kap_out,GREY_fsmall_out,&
      MG_fsmall_out,res_history_out,outfile,restart_outfile,decomp_outfile,TEMP_outfile,GREY_E_outfile,GREY_F_outfile,&
      GREY_kap_outfile,GREY_fsmall_outfile,MG_fsmall_outfile,res_history_outfile,HO_E_out,HO_E_outfile,delx,dely,N_g,&
      nu_g,N_t,sig_r,ar,sig0,pi,c,h,cv,Omega_x,Omega_y,quad_weight,quadrature,bc_type)
  END IF

END SUBROUTINE SPECIFIC_OUTPUTS

!============================================================================================================!
!
!============================================================================================================!
SUBROUTINE SPCF_OUTPUT(Data,Times,File,Data_Name,f_unit,datevalues,database_gen,database_add,run_type,restart_infile,&
  use_grey,chi,conv_ho,conv_lo,conv_gr,comp_unit,line_src,maxit_RTE,maxit_MLOQD,maxit_GLOQD,conv_type,N_m,threads,&
  kapE_dT_flag,enrgy_strc,Nwt_upbnd,xlen,ylen,N_x,N_y,tlen,delt,bcT_left,bcT_right,bcT_upper,bcT_lower,Tini,&
  out_times,out_time_steps,n_out_times,restart_outlen,TEMP_out,GREY_E_out,GREY_F_out,GREY_kap_out,GREY_fsmall_out,&
  MG_fsmall_out,res_history_out,outfile,restart_outfile,decomp_outfile,TEMP_outfile,GREY_E_outfile,GREY_F_outfile,&
  GREY_kap_outfile,GREY_fsmall_outfile,MG_fsmall_outfile,res_history_outfile,HO_E_out,HO_E_outfile,delx,dely,N_g,&
  nu_g,N_t,sig_r,ar,sig0,pi,c,h,cv,Omega_x,Omega_y,quad_weight,quadrature,bc_type)

  REAL*8,INTENT(IN):: Data(:,:,:)
  REAL*8,INTENT(IN):: Times(:)
  CHARACTER(100),INTENT(IN):: File, Data_Name

  !INPUT VARIABLES
  INTEGER,INTENT(IN):: f_unit, datevalues(:)
  REAL*8,INTENT(IN):: chi, conv_ho, conv_lo, conv_gr, line_src, Nwt_upbnd, comp_unit
  REAL*8,INTENT(IN):: xlen, ylen, tlen, delt, bcT_left, bcT_right, bcT_upper, bcT_lower, Tini
  REAL*8,INTENT(IN):: Delx(:), Dely(:)
  REAL*8,INTENT(IN):: sig_r, ar, sig0, pi, c, h, cv
  REAL*8,INTENT(IN):: out_times(:), nu_g(:), Omega_x(:), Omega_y(:), quad_weight(:)
  INTEGER,INTENT(IN):: database_gen, database_add, use_grey, bc_type(:)
  INTEGER,INTENT(IN):: maxit_RTE, maxit_MLOQD, maxit_GLOQD, conv_type, N_m, threads, N_g, N_t
  INTEGER,INTENT(IN):: N_x, N_y
  INTEGER,INTENT(IN):: n_out_times, restart_outlen, TEMP_out, GREY_E_out, GREY_F_out, HO_E_out
  INTEGER,INTENT(IN):: GREY_kap_out, GREY_fsmall_out, MG_fsmall_out, res_history_out
  INTEGER,INTENT(IN):: out_time_steps(:)
  CHARACTER(100),INTENT(IN):: outfile, restart_outfile, decomp_outfile, TEMP_outfile, GREY_E_outfile, HO_E_outfile
  CHARACTER(100),INTENT(IN):: GREY_F_outfile, GREY_kap_outfile, GREY_fsmall_outfile, MG_fsmall_outfile, res_history_outfile
  CHARACTER(100),INTENT(IN):: run_type, restart_infile
  CHARACTER(100),INTENT(IN):: kapE_dT_flag, enrgy_strc, quadrature

  INTEGER:: N_t2, t, i, j
  INTEGER:: err

  N_t2 = SIZE(Times,1)

  OPEN(UNIT=f_unit,FILE=File,STATUS='REPLACE',ACTION='WRITE',IOSTAT=err)
!   making sure file exists/opens, if not tells user
  IF (err .NE. 0) THEN
      WRITE(*,'(3A)') 'The file, ',File,', could not open properly.'
      STOP
  END IF

  CALL OUTFILE_HEADER(f_unit,datevalues,database_gen,database_add,run_type,restart_infile,use_grey,chi,&
    conv_ho,conv_lo,conv_gr,comp_unit,line_src,maxit_RTE,maxit_MLOQD,maxit_GLOQD,conv_type,N_m,threads,&
    kapE_dT_flag,enrgy_strc,Nwt_upbnd,xlen,ylen,N_x,N_y,tlen,delt,bc_type,bcT_left,bcT_right,bcT_upper,&
    bcT_lower,Tini,out_times,out_time_steps,n_out_times,restart_outlen,TEMP_out,GREY_E_out,GREY_F_out,&
    GREY_kap_out,GREY_fsmall_out,MG_fsmall_out,res_history_out,outfile,restart_outfile,decomp_outfile,TEMP_outfile,&
    GREY_E_outfile,GREY_F_outfile,GREY_kap_outfile,GREY_fsmall_outfile,MG_fsmall_outfile,res_history_outfile,&
    HO_E_out,HO_E_outfile,delx,dely,N_g,nu_g,N_t,sig_r,ar,sig0,pi,c,h,cv,Omega_x,Omega_y,quad_weight,quadrature)

  WRITE(f_unit,'(A)') '-----------------------------------------------------------'
  WRITE(f_unit,'(A)') 'OUT_TIMES:'
  WRITE(f_unit,*)
  WRITE(f_unit,'(A)') 'time_step  |  time (ns)  |  c*t'
  DO t=1,N_t2
    WRITE(f_unit,'(I5.5,A,ES12.5,A,ES12.5)') t,' , ',Times(t),' , ',c*Times(t)
  END DO

  WRITE(f_unit,'(A)') '-----------------------------------------------------------'
  WRITE(f_unit,'(A)') trim(Data_Name)
  WRITE(f_unit,*)

  WRITE(f_unit,'(A)',ADVANCE='NO') ' i   j        x             y               '
  DO t=1,N_t2
    WRITE(f_unit,'(I5.5,A)',ADVANCE='NO') t,'                  '
  END DO
  WRITE(f_unit,*)
  DO j=1,N_y
    DO i=1,N_x
      WRITE(f_unit,'(I3.3,A,I3.3,2ES14.6)',ADVANCE='NO') i,' ', j, (DBLE(i)-.5d0)*Delx(i), (DBLE(j)-.5d0)*Dely(j)
      DO t=1,N_t2
        WRITE(f_unit,'(ES23.15)',ADVANCE='NO') Data(i,j,t)
      END DO
      WRITE(f_unit,*)
    END DO
  END DO

  CLOSE ( f_unit, STATUS='KEEP')

END SUBROUTINE SPCF_OUTPUT

!==================================================================================================================================!
!
!==================================================================================================================================!

END MODULE OUTPUTS
