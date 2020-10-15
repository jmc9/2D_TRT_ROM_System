#include <stdlib.h>
#include <stdio.h>
#include <netcdf.h>
#include <string.h>
#include <math.h>

void INPUT(char *infile, char *dsfile, char *outfile);
void GET_DIMS(int ncid, size_t *N_t, size_t *N_g, size_t *N_y, size_t *N_x);
void HANDLE_ERR(int Status, char Location[]);
void OPEN_NCFILE(char *fname, int *ncid);

int GENERATE_POD(int ncid_in, int ncid_out, char *dname, size_t len, size_t N_t, size_t N_g, size_t N_y, size_t N_x, int gsum,
  int Sid, int Uid, int Vtid);

int DEF_DIMS(int ncid_out, size_t N_t, size_t N_g, size_t N_y, size_t N_x, size_t rank_avg, size_t rank_edgV, size_t rank_edgH,
  size_t clen_avg, size_t clen_edgV, size_t clen_edgH, int *N_t_ID, int *N_g_ID, int *N_y_ID, int *N_x_ID, int *rank_avg_ID,
  int *rank_edgV_ID, int *rank_edgH_ID, int *clen_avg_ID, int *clen_edgV_ID, int *clen_edgH_ID);

int DEF_VARS(int ncid, int N_g_ID, int rank_avg_ID, int rank_edgV_ID, int rank_edgH_ID, int clen_avg_ID, int clen_edgV_ID,
   int clen_edgH_ID, int N_t_ID, int *S_fg_avg_xx_ID, int *U_fg_avg_xx_ID, int *Vt_fg_avg_xx_ID, int *S_fg_edgV_xx_ID,
   int *U_fg_edgV_xx_ID, int *Vt_fg_edgV_xx_ID, int *S_fg_avg_yy_ID, int *U_fg_avg_yy_ID, int *Vt_fg_avg_yy_ID,
   int *S_fg_edgH_yy_ID, int *U_fg_edgH_yy_ID, int *Vt_fg_edgH_yy_ID, int *S_fg_edgV_xy_ID, int *U_fg_edgV_xy_ID,
   int *Vt_fg_edgV_xy_ID, int *S_fg_edgH_xy_ID, int *U_fg_edgH_xy_ID, int *Vt_fg_edgH_xy_ID);

int OUTPUT_fg_POD(int ncid_in, int ncid_out, size_t N_t, size_t N_g, size_t N_y, size_t N_x, int gsum, int S_fg_avg_xx_ID,
  int U_fg_avg_xx_ID, int Vt_fg_avg_xx_ID, int S_fg_edgV_xx_ID, int U_fg_edgV_xx_ID, int Vt_fg_edgV_xx_ID,
  int S_fg_avg_yy_ID, int U_fg_avg_yy_ID, int Vt_fg_avg_yy_ID, int S_fg_edgH_yy_ID, int U_fg_edgH_yy_ID,
  int Vt_fg_edgH_yy_ID, int S_fg_edgV_xy_ID, int U_fg_edgV_xy_ID, int Vt_fg_edgV_xy_ID, int S_fg_edgH_xy_ID,
  int U_fg_edgH_xy_ID, int Vt_fg_edgH_xy_ID);

#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })

//================================================================================================================================//
//
//================================================================================================================================//
int main()
{
  // comment
  int ncid_in, ncid_out, err;
  int N_t_ID, N_g_ID, N_y_ID, N_x_ID;
  int rank_avg_ID, rank_edgV_ID, rank_edgH_ID, clen_avg_ID, clen_edgV_ID, clen_edgH_ID;
  int ndims, dimids[4];
  int S_fg_avg_xx_ID, U_fg_avg_xx_ID, Vt_fg_avg_xx_ID, S_fg_edgV_xx_ID, U_fg_edgV_xx_ID, Vt_fg_edgV_xx_ID;
  int S_fg_avg_yy_ID, U_fg_avg_yy_ID, Vt_fg_avg_yy_ID, S_fg_edgH_yy_ID, U_fg_edgH_yy_ID, Vt_fg_edgH_yy_ID;
  int S_fg_edgV_xy_ID, U_fg_edgV_xy_ID, Vt_fg_edgV_xy_ID, S_fg_edgH_xy_ID, U_fg_edgH_xy_ID, Vt_fg_edgH_xy_ID;
  char infile[25], dsfile[100], dname[25], outfile[25];
  size_t N_t, N_g, N_y, N_x, glen, glen_edgV, glen_edgH;
  size_t rank_avg, rank_edgV, rank_edgH, clen_avg, clen_edgV, clen_edgH;
  char loc[5] = "main";

  strcpy(infile,"input.inp"); //setting input file name (default)

  INPUT(infile,dsfile,outfile); //reading input file
  printf("%s %s\n",dsfile,outfile);

  err = nc_open(dsfile,NC_NOWRITE,&ncid_in); HANDLE_ERR(err,loc); //opening NetCDF dataset
  GET_DIMS(ncid_in,&N_t,&N_g,&N_y,&N_x); //finding dimensions of problem domain
  glen = N_t*N_g*N_y*N_x;
  glen_edgV = N_t*N_g*N_y*(N_x+1);
  glen_edgH = N_t*N_g*(N_y+1)*N_x;

  rank_avg = min(N_x*N_y,N_t); rank_edgV = min((N_x+1)*N_y,N_t); rank_edgH = min(N_x*(N_y+1),N_t);
  clen_avg = N_x*N_y; clen_edgV = (N_x+1)*N_y; clen_edgH = N_x*(N_y+1);

  err = nc_create(outfile,NC_CLOBBER,&ncid_out); HANDLE_ERR(err,loc); //opening NetCDF dataset

  err = DEF_DIMS(ncid_out,N_t,N_g,N_y,N_x,rank_avg,rank_edgV,rank_edgH,clen_avg,clen_edgV,clen_edgH,&N_t_ID,&N_g_ID,&N_y_ID,
    &N_x_ID,&rank_avg_ID,&rank_edgV_ID,&rank_edgH_ID,&clen_avg_ID,&clen_edgV_ID,&clen_edgH_ID);

  err = DEF_VARS(ncid_out,N_g_ID,rank_avg_ID,rank_edgV_ID,rank_edgH_ID,clen_avg_ID,clen_edgV_ID,clen_edgH_ID,N_t_ID,
    &S_fg_avg_xx_ID,&U_fg_avg_xx_ID,&Vt_fg_avg_xx_ID,&S_fg_edgV_xx_ID,&U_fg_edgV_xx_ID,&Vt_fg_edgV_xx_ID,&S_fg_avg_yy_ID,
    &U_fg_avg_yy_ID,&Vt_fg_avg_yy_ID,&S_fg_edgH_yy_ID,&U_fg_edgH_yy_ID,&Vt_fg_edgH_yy_ID,&S_fg_edgV_xy_ID,&U_fg_edgV_xy_ID,
    &Vt_fg_edgV_xy_ID,&S_fg_edgH_xy_ID,&U_fg_edgH_xy_ID,&Vt_fg_edgH_xy_ID);

  err = nc_enddef(ncid_out); HANDLE_ERR(err,loc);

  err = OUTPUT_fg_POD(ncid_in,ncid_out,N_t,N_g,N_y,N_x,0,
     S_fg_avg_xx_ID,U_fg_avg_xx_ID,Vt_fg_avg_xx_ID,S_fg_edgV_xx_ID,
     U_fg_edgV_xx_ID,Vt_fg_edgV_xx_ID,S_fg_avg_yy_ID,U_fg_avg_yy_ID,Vt_fg_avg_yy_ID,
     S_fg_edgH_yy_ID,U_fg_edgH_yy_ID,Vt_fg_edgH_yy_ID,S_fg_edgV_xy_ID,U_fg_edgV_xy_ID,
     Vt_fg_edgV_xy_ID,S_fg_edgH_xy_ID,U_fg_edgH_xy_ID,Vt_fg_edgH_xy_ID);

  err = nc_close(ncid_in); HANDLE_ERR(err,loc); //closing NetCDF dataset
  err = nc_close(ncid_out); HANDLE_ERR(err,loc); //closing NetCDF dataset

  printf("yuh\n");
}

//================================================================================================================================//
//
//================================================================================================================================//
