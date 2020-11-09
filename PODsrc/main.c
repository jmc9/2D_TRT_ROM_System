#include <stdlib.h>
#include <stdio.h>
#include <netcdf.h>
#include <string.h>
#include <math.h>

//FROM NCDF_IO.c
void HANDLE_ERR(const int Status, const char *Location);

//FROM INPUTS.c
void INPUT(const char *infile, char *dsfile, char *outfile, int *gsum, int *fg_pod, int *Ig_pod, int *Mean_Ig_pod);
void GET_DIMS(const int ncid, size_t *N_t, size_t *N_g, size_t *N_m, size_t *N_y, size_t *N_x, double *tlen, double *Delt,
  double *xlen, double *ylen, double **Delx, double **Dely, int *BC_Type, double *bcT, double *Tini);

//FROM OUTPUTS.c
int DEF_DIMS(const int ncid_out, const size_t N_t, const size_t N_g, const size_t N_m, const size_t N_y, const size_t N_x,
  const size_t rank_BC, const size_t rank_avg, const size_t rank_edgV,const size_t rank_edgH, const size_t BClen,
  const size_t clen_avg, const size_t clen_edgV, const size_t clen_edgH, const double tlen, const double Delt, const double xlen,
  const double ylen, double *Delx, double *Dely, int *BC_Type, double *bcT, const double Tini, int *N_t_ID, int *N_g_ID,
  int *N_m_ID, int *N_y_ID, int *N_x_ID, int *rank_BC_ID, int *rank_avg_ID, int *rank_edgV_ID, int *rank_edgH_ID, int *BClen_ID,
  int *clen_avg_ID, int *clen_edgV_ID, int *clen_edgH_ID);

int DEF_GBCs(const int ncid, const int gsum, const int N_g_ID, const int rank_BC_ID, const int BClen_ID, const int N_t_ID,
   int *C_BCg_ID, int *S_BCg_ID, int *U_BCg_ID, int *Vt_BCg_ID);

int DEF_fg_VARS(const int ncid, const int gsum, const int N_g_ID, const int rank_avg_ID, const int rank_edgV_ID,
   const int rank_edgH_ID, const int clen_avg_ID,const int clen_edgV_ID, const int clen_edgH_ID, const int N_t_ID,
   int *C_fg_avg_xx_ID, int *S_fg_avg_xx_ID, int *U_fg_avg_xx_ID,int *Vt_fg_avg_xx_ID, int *C_fg_edgV_xx_ID,
   int *S_fg_edgV_xx_ID, int *U_fg_edgV_xx_ID, int *Vt_fg_edgV_xx_ID, int *C_fg_avg_yy_ID, int *S_fg_avg_yy_ID,
   int *U_fg_avg_yy_ID, int *Vt_fg_avg_yy_ID,int *C_fg_edgH_yy_ID, int *S_fg_edgH_yy_ID, int *U_fg_edgH_yy_ID,
   int *Vt_fg_edgH_yy_ID, int *C_fg_edgV_xy_ID, int *S_fg_edgV_xy_ID, int *U_fg_edgV_xy_ID, int *Vt_fg_edgV_xy_ID,
   int *C_fg_edgH_xy_ID, int *S_fg_edgH_xy_ID, int *U_fg_edgH_xy_ID, int *Vt_fg_edgH_xy_ID);

int DEF_Ig_VARS(const int ncid, const int gsum, const int N_g_ID, const int rank_avg_ID, const int rank_edgV_ID,
  const int rank_edgH_ID, const int clen_avg_ID,const int clen_edgV_ID, const int clen_edgH_ID, const int N_t_ID,
  int *C_Ig_avg_ID, int *S_Ig_avg_ID, int *U_Ig_avg_ID, int *Vt_Ig_avg_ID,int *C_Ig_edgV_ID, int *S_Ig_edgV_ID,
  int *U_Ig_edgV_ID, int *Vt_Ig_edgV_ID, int *C_Ig_edgH_ID, int *S_Ig_edgH_ID, int *U_Ig_edgH_ID, int *Vt_Ig_edgH_ID);

int DEF_meanIg_VARS(const int ncid, const int gsum, const int N_g_ID, const int rank_avg_ID, const int rank_edgV_ID,
  const int rank_edgH_ID, const int clen_avg_ID,const int clen_edgV_ID, const int clen_edgH_ID, const int N_t_ID,
  int *C_Ig_avg_ID, int *S_Ig_avg_ID, int *U_Ig_avg_ID, int *Vt_Ig_avg_ID,int *C_Ig_edgV_ID, int *S_Ig_edgV_ID,
  int *U_Ig_edgV_ID, int *Vt_Ig_edgV_ID, int *C_Ig_edgH_ID, int *S_Ig_edgH_ID, int *U_Ig_edgH_ID, int *Vt_Ig_edgH_ID);

int OUTPUT_BCg_POD(const int ncid_in, const int ncid_out, const size_t N_t, const size_t N_g, const size_t N_y, const size_t N_x,
  const int gsum, const int C_BCg_ID,const int S_BCg_ID, const int U_BCg_ID, const int Vt_BCg_ID);

int OUTPUT_fg_POD(const int ncid_in, const int ncid_out, const size_t N_t, const size_t N_g, const size_t N_y, const size_t N_x,
  const int gsum, const int C_fg_avg_xx_ID,const int S_fg_avg_xx_ID, const int U_fg_avg_xx_ID, const int Vt_fg_avg_xx_ID,
  const int C_fg_edgV_xx_ID, const int S_fg_edgV_xx_ID, const int U_fg_edgV_xx_ID,const int Vt_fg_edgV_xx_ID,
  const int C_fg_avg_yy_ID, const int S_fg_avg_yy_ID, const int U_fg_avg_yy_ID, const int Vt_fg_avg_yy_ID,
  const int C_fg_edgH_yy_ID, const int S_fg_edgH_yy_ID, const int U_fg_edgH_yy_ID, const int Vt_fg_edgH_yy_ID,
  const int C_fg_edgV_xy_ID, const int S_fg_edgV_xy_ID, const int U_fg_edgV_xy_ID, const int Vt_fg_edgV_xy_ID,
  const int C_fg_edgH_xy_ID, const int S_fg_edgH_xy_ID, const int U_fg_edgH_xy_ID, const int Vt_fg_edgH_xy_ID);

int OUTPUT_Ig_POD(const int ncid_in, const int ncid_out, const size_t N_t, const size_t N_g, const size_t N_m, const size_t N_y,
  const size_t N_x, const int gsum, const int C_Ig_avg_ID, const int S_Ig_avg_ID, const int U_Ig_avg_ID, const int Vt_Ig_avg_ID,
  const int C_Ig_edgV_ID, const int S_Ig_edgV_ID, const int U_Ig_edgV_ID, const int Vt_Ig_edgV_ID, const int C_Ig_edgH_ID,
  const int S_Ig_edgH_ID, const int U_Ig_edgH_ID, const int Vt_Ig_edgH_ID);

int OUTPUT_meanIg_POD(const int ncid_in, const int ncid_out, const size_t N_t, const size_t N_g, const size_t N_m, const size_t N_y,
  const size_t N_x, const int gsum, const int C_Ig_avg_ID, const int S_Ig_avg_ID, const int U_Ig_avg_ID, const int Vt_Ig_avg_ID,
  const int C_Ig_edgV_ID, const int S_Ig_edgV_ID, const int U_Ig_edgV_ID, const int Vt_Ig_edgV_ID, const int C_Ig_edgH_ID,
  const int S_Ig_edgH_ID, const int U_Ig_edgH_ID, const int Vt_Ig_edgH_ID);

//LOCAL DEFINITIONS
#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })

//================================================================================================================================//
//
//================================================================================================================================//
int main()
{
  //
  int ncid_in, ncid_out, err, gsum, fg_pod, Ig_pod, Mean_Ig_pod;
  int N_t_ID, N_g_ID, N_m_ID, N_y_ID, N_x_ID;
  int rank_BC_ID, rank_avg_ID, rank_edgV_ID, rank_edgH_ID, BClen_ID, clen_avg_ID, clen_edgV_ID, clen_edgH_ID;
  //
  int C_BCg_ID, S_BCg_ID, U_BCg_ID, Vt_BCg_ID;
  int C_fg_avg_xx_ID, S_fg_avg_xx_ID, U_fg_avg_xx_ID, Vt_fg_avg_xx_ID;
  int C_fg_edgV_xx_ID, S_fg_edgV_xx_ID, U_fg_edgV_xx_ID, Vt_fg_edgV_xx_ID;
  int C_fg_avg_yy_ID, S_fg_avg_yy_ID, U_fg_avg_yy_ID, Vt_fg_avg_yy_ID;
  int C_fg_edgH_yy_ID, S_fg_edgH_yy_ID, U_fg_edgH_yy_ID, Vt_fg_edgH_yy_ID;
  int C_fg_edgV_xy_ID, S_fg_edgV_xy_ID, U_fg_edgV_xy_ID, Vt_fg_edgV_xy_ID;
  int C_fg_edgH_xy_ID, S_fg_edgH_xy_ID, U_fg_edgH_xy_ID, Vt_fg_edgH_xy_ID;
  //
  int C_Ig_avg_ID, S_Ig_avg_ID, U_Ig_avg_ID, Vt_Ig_avg_ID;
  int C_Ig_edgV_ID, S_Ig_edgV_ID, U_Ig_edgV_ID, Vt_Ig_edgV_ID;
  int C_Ig_edgH_ID, S_Ig_edgH_ID, U_Ig_edgH_ID, Vt_Ig_edgH_ID;
  //
  double tlen, Delt, xlen, ylen, *Delx, *Dely, bcT[4], Tini;
  int BC_Type[4];
  size_t N_t, N_g, N_m, N_y, N_x, gscale, mscale;
  size_t rank_BC, rank_avg, rank_edgV, rank_edgH, BClen, clen_avg, clen_edgV, clen_edgH;
  char infile[25], dsfile[100], outfile[25];
  char loc[5] = "main";

  //===========================================================================//
  //                                                                           //
  //     READING INPUTS                                                        //
  //                                                                           //
  //===========================================================================//
  printf("Program start\nReading inputs\n");
  strcpy(infile,"input.inp"); //setting input file name (default)

  INPUT(infile,dsfile,outfile,&gsum,&fg_pod,&Ig_pod,&Mean_Ig_pod); //reading input file

  err = nc_open(dsfile,NC_NOWRITE,&ncid_in); HANDLE_ERR(err,loc); //opening NetCDF dataset
  GET_DIMS(ncid_in,&N_t,&N_g,&N_m,&N_y,&N_x,&tlen,&Delt,&xlen,&ylen,&Delx,&Dely,BC_Type,bcT,&Tini); //finding dimensions of problem domain

  if (gsum == 1){ gscale = N_g; }
  else{ gscale = 1; }
  if ((Ig_pod == 1) || (Mean_Ig_pod == 1)){ mscale = N_m; }
  else{ mscale = 1; }
  clen_avg = gscale*mscale*N_x*N_y; clen_edgV = gscale*mscale*(N_x+1)*N_y; clen_edgH = gscale*mscale*N_x*(N_y+1);
  rank_avg = min(clen_avg,N_t); rank_edgV = min(clen_edgV,N_t); rank_edgH = min(clen_edgH,N_t);

  if(fg_pod == 1){
    BClen = gscale*2*(N_x+N_y);
    rank_BC = min(BClen,N_t);
  }
  else{
    BClen = 0;
    rank_BC = 0;
  }

  //===========================================================================//
  //                                                                           //
  //     INITIALIZING OUTPUT FILE                                              //
  //                                                                           //
  //===========================================================================//
  printf("Initializing output file\n");
  err = nc_create(outfile,NC_NETCDF4,&ncid_out); HANDLE_ERR(err,loc); //opening NetCDF dataset

  err = DEF_DIMS(ncid_out,N_t,N_g,N_m,N_y,N_x,rank_BC,rank_avg,rank_edgV,rank_edgH,BClen,clen_avg,clen_edgV,clen_edgH,tlen,Delt,
    xlen,ylen,Delx,Dely,BC_Type,bcT,Tini,&N_t_ID,&N_g_ID,&N_m_ID,&N_y_ID,&N_x_ID,&rank_BC_ID,&rank_avg_ID,&rank_edgV_ID,
    &rank_edgH_ID,&BClen_ID,&clen_avg_ID,&clen_edgV_ID,&clen_edgH_ID);
  if (err != 0){ printf("Error occured while defining dimensions in %s\n",outfile); exit(1); }

  if (fg_pod == 1){
    err = DEF_GBCs(ncid_out,gsum,N_g_ID,rank_BC_ID,BClen_ID,N_t_ID,&C_BCg_ID,&S_BCg_ID,&U_BCg_ID,&Vt_BCg_ID);

    err = DEF_fg_VARS(ncid_out,gsum,N_g_ID,rank_avg_ID,rank_edgV_ID,rank_edgH_ID,clen_avg_ID,clen_edgV_ID,clen_edgH_ID,N_t_ID,
      &C_fg_avg_xx_ID,&S_fg_avg_xx_ID,&U_fg_avg_xx_ID,&Vt_fg_avg_xx_ID,&C_fg_edgV_xx_ID,&S_fg_edgV_xx_ID,&U_fg_edgV_xx_ID,
      &Vt_fg_edgV_xx_ID,&C_fg_avg_yy_ID,&S_fg_avg_yy_ID,&U_fg_avg_yy_ID,&Vt_fg_avg_yy_ID,&C_fg_edgH_yy_ID,&S_fg_edgH_yy_ID,
      &U_fg_edgH_yy_ID,&Vt_fg_edgH_yy_ID,&C_fg_edgV_xy_ID,&S_fg_edgV_xy_ID,&U_fg_edgV_xy_ID,&Vt_fg_edgV_xy_ID,&C_fg_edgH_xy_ID,
      &S_fg_edgH_xy_ID,&U_fg_edgH_xy_ID,&Vt_fg_edgH_xy_ID);
  }
  else if (Ig_pod == 1){
    err = DEF_Ig_VARS(ncid_out,gsum,N_g_ID,rank_avg_ID,rank_edgV_ID,rank_edgH_ID,clen_avg_ID,clen_edgV_ID,clen_edgH_ID,
       N_t_ID,&C_Ig_avg_ID,&S_Ig_avg_ID,&U_Ig_avg_ID,&Vt_Ig_avg_ID,&C_Ig_edgV_ID,&S_Ig_edgV_ID,&U_Ig_edgV_ID,&Vt_Ig_edgV_ID,
       &C_Ig_edgH_ID,&S_Ig_edgH_ID,&U_Ig_edgH_ID,&Vt_Ig_edgH_ID);
  }
  else if (Mean_Ig_pod == 1){
    err = DEF_meanIg_VARS(ncid_out,gsum,N_g_ID,rank_avg_ID,rank_edgV_ID,rank_edgH_ID,clen_avg_ID,clen_edgV_ID,clen_edgH_ID,
       N_t_ID,&C_Ig_avg_ID,&S_Ig_avg_ID,&U_Ig_avg_ID,&Vt_Ig_avg_ID,&C_Ig_edgV_ID,&S_Ig_edgV_ID,&U_Ig_edgV_ID,&Vt_Ig_edgV_ID,
       &C_Ig_edgH_ID,&S_Ig_edgH_ID,&U_Ig_edgH_ID,&Vt_Ig_edgH_ID);
  }

  if (err != 0){ printf("Error occured while defining variables in %s\n",outfile); exit(1); }

  err = nc_enddef(ncid_out); HANDLE_ERR(err,loc); //changing output file to data mode

  //===========================================================================//
  //                                                                           //
  //     PERFORMING THE POD                                                    //
  //                                                                           //
  //===========================================================================//
  printf("Beginning POD calculations...\n");

  if (fg_pod == 1){
    err = OUTPUT_BCg_POD(ncid_in,ncid_out,N_t,N_g,N_y,N_x,gsum,C_BCg_ID,S_BCg_ID,U_BCg_ID,Vt_BCg_ID);

    err = OUTPUT_fg_POD(ncid_in,ncid_out,N_t,N_g,N_y,N_x,gsum,C_fg_avg_xx_ID,S_fg_avg_xx_ID,U_fg_avg_xx_ID,Vt_fg_avg_xx_ID,
       C_fg_edgV_xx_ID,S_fg_edgV_xx_ID,U_fg_edgV_xx_ID,Vt_fg_edgV_xx_ID,C_fg_avg_yy_ID,S_fg_avg_yy_ID,U_fg_avg_yy_ID,
       Vt_fg_avg_yy_ID,C_fg_edgH_yy_ID,S_fg_edgH_yy_ID,U_fg_edgH_yy_ID,Vt_fg_edgH_yy_ID,C_fg_edgV_xy_ID,S_fg_edgV_xy_ID,
       U_fg_edgV_xy_ID,Vt_fg_edgV_xy_ID,C_fg_edgH_xy_ID,S_fg_edgH_xy_ID,U_fg_edgH_xy_ID,Vt_fg_edgH_xy_ID);
  }
  else if (Ig_pod == 1){
    err = OUTPUT_Ig_POD(ncid_in,ncid_out,N_t,N_g,N_m,N_y,N_x,gsum,C_Ig_avg_ID,S_Ig_avg_ID,U_Ig_avg_ID,Vt_Ig_avg_ID,
       C_Ig_edgV_ID,S_Ig_edgV_ID,U_Ig_edgV_ID,Vt_Ig_edgV_ID,C_Ig_edgH_ID,S_Ig_edgH_ID,U_Ig_edgH_ID,Vt_Ig_edgH_ID);
  }
  else if (Mean_Ig_pod == 1){
    err = OUTPUT_meanIg_POD(ncid_in,ncid_out,N_t,N_g,N_m,N_y,N_x,gsum,C_Ig_avg_ID,S_Ig_avg_ID,U_Ig_avg_ID,Vt_Ig_avg_ID,
       C_Ig_edgV_ID,S_Ig_edgV_ID,U_Ig_edgV_ID,Vt_Ig_edgV_ID,C_Ig_edgH_ID,S_Ig_edgH_ID,U_Ig_edgH_ID,Vt_Ig_edgH_ID);
  }

  //===========================================================================//
  //                                                                           //
  //     END OF PROGRAM CLEAN UP                                               //
  //                                                                           //
  //===========================================================================//
  err = nc_close(ncid_in); HANDLE_ERR(err,loc); //closing NetCDF dataset
  err = nc_close(ncid_out); HANDLE_ERR(err,loc); //closing NetCDF dataset

  printf("Program completed\n");
}

//================================================================================================================================//
//
//================================================================================================================================//
