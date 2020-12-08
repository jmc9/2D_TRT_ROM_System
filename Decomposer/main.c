/*================================================================================================================================*/
/*
  Decomposer.c
  Author: Joseph M. Coale
  jmcoale@ncsu.edu - josephcoale912@gmail.com
*/
/*================================================================================================================================*/
/* Importing Packages */
/*================================================================================================================================*/
#include <stdlib.h>
#include <stdio.h>
#include <netcdf.h>
#include <string.h>
#include <math.h>

/*================================================================================================================================*/
/* Importing Functions */
/*================================================================================================================================*/
/* ----- FROM NCDF_IO.c ----- */
void Handle_Err(const int Status, const char *Location);

/* ----- FROM INPUTS.c ----- */
int Input(const char *infile, char *dsfile, char *outfile, int *dcmp_type, int *dcmp_data, int *gsum);

void Get_Dims(const int ncid, size_t *N_t, size_t *N_g, size_t *N_m, size_t *N_y, size_t *N_x, double *tlen, double *Delt,
  double *xlen, double *ylen, double **Delx, double **Dely, int *BC_Type, double *bcT, double *Tini);

/* ----- FROM OUTPUTS.c ----- */
int Def_Dims(const int ncid_out, const size_t N_t, const size_t N_g, const size_t N_m, const size_t N_y, const size_t N_x,
  const size_t rank_BC, const size_t rank_avg, const size_t rank_edgV,const size_t rank_edgH, const size_t BClen,
  const size_t clen_avg, const size_t clen_edgV, const size_t clen_edgH, const double tlen, const double Delt, const double xlen,
  const double ylen, double *Delx, double *Dely, int *BC_Type, double *bcT, const double Tini, int *N_t_ID, int *N_g_ID,
  int *N_m_ID, int *N_y_ID, int *N_x_ID, int *rank_BC_ID, int *rank_avg_ID, int *rank_edgV_ID, int *rank_edgH_ID, int *BClen_ID,
  int *clen_avg_ID, int *clen_edgV_ID, int *clen_edgH_ID);

int Def_DCMP_Vars(const int ncid, const int dcmp_type, const int dcmp_data, const int gsum, const int N_g_ID, const int rank_BC_ID,
  const int rank_avg_ID,const int rank_edgV_ID, const int rank_edgH_ID, const int BClen_ID, const int clen_avg_ID, const int clen_edgV_ID,
  const int clen_edgH_ID, const int N_t_ID, int **BCg_IDs, int **fg_avg_xx_IDs, int **fg_edgV_xx_IDs, int **fg_avg_yy_IDs, int **fg_edgH_yy_IDs,
  int **fg_edgV_xy_IDs, int **fg_edgH_xy_IDs, int **Ig_avg_IDs, int **Ig_edgV_IDs, int **Ig_edgH_IDs);

// int Def_Pod_Vars(const int ncid, const char *vname, const int rank_ID, const int clen_ID, const int N_t_ID, int *DCMP_IDs);
int Decompose_Data(const int ncid_in, const int ncid_out, const int dcmp_type, const int dcmp_data, const int gsum, const size_t N_t,
  const size_t N_g, const size_t N_m, const size_t N_x, const size_t N_y, const int *BCg_IDs, const int *fg_avg_xx_IDs, const int *fg_edgV_xx_IDs,
  const int *fg_avg_yy_IDs, const int *fg_edgH_yy_IDs,const int *fg_edgV_xy_IDs, const int *fg_edgH_xy_IDs, const int *Ig_avg_IDs,
  const int *Ig_edgV_IDs, const int *Ig_edgH_IDs);

/* ----- LOCAL DEFINITIONS ----- */
#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })

/*================================================================================================================================*/
/* Main Program */
/*================================================================================================================================*/
int main()
{
  //basic (local) program variables
  int err;
  char loc[5] = "main";

  //input variables
  int dcmp_data, dcmp_type, gsum;
  int ncid_in, ncid_out;
  char infile[25], outfile[25], dsfile[100];

  //problem parameters
  double tlen, Delt, xlen, ylen, *Delx, *Dely, bcT[4], Tini;
  int BC_Type[4];
  size_t N_t, N_g, N_m, N_y, N_x;
  size_t mscale, gscale;
  size_t rank_BC, rank_avg, rank_edgV, rank_edgH, BClen, clen_avg, clen_edgV, clen_edgH;
  int N_t_ID, N_g_ID, N_m_ID, N_y_ID, N_x_ID;
  int rank_BC_ID, rank_avg_ID, rank_edgV_ID, rank_edgH_ID, BClen_ID, clen_avg_ID, clen_edgV_ID, clen_edgH_ID;

  //QDf decomposition IDs
  int *BCg_IDs, *fg_avg_xx_IDs, *fg_edgV_xx_IDs, *fg_avg_yy_IDs, *fg_edgH_yy_IDs, *fg_edgV_xy_IDs, *fg_edgH_xy_IDs;

  //I decomposition IDs
  int *Ig_avg_IDs, *Ig_edgV_IDs, *Ig_edgH_IDs;


  //big text header output to terminal on program execution
  printf("  ____                                                           \n");
  printf(" |  _ \\  ___  ___ ___  _ __ ___  _ __   ___  ___  ___ _ __      \n");
  printf(" | | | |/ _ \\/ __/ _ \\| '_ ` _ \\| '_ \\ / _ \\/ __|/ _ \\ '__|\n");
  printf(" | |_| |  __/ (_| (_) | | | | | | |_) | (_) \\__ \\  __/ |       \n");
  printf(" |____/ \\___|\\___\\___/|_| |_| |_| .__/ \\___/|___/\\___|_|    \n");
  printf("                                |_|                              \n");
  printf("\n");

  /*===========================================================================*/
  /*                                                                            /
  /                               READING INPUTS                                /
  /                                                                            */
  /*===========================================================================*/
  printf("Program start\nReading inputs\n");
  strcpy(infile,"input.inp"); //setting input file name (default)

  //opening, reading input file
  err = Input(infile,dsfile,outfile,&dcmp_type,&dcmp_data,&gsum);
  if (err != 0){
    printf("Error detected upon user input, aborting program\n");
    exit(1);
  }

  //opening, parsing datafile
  printf("Reading datafile\n");
  err = nc_open(dsfile,NC_NOWRITE,&ncid_in); Handle_Err(err,loc); //opening NetCDF dataset
  Get_Dims(ncid_in,&N_t,&N_g,&N_m,&N_y,&N_x,&tlen,&Delt,&xlen,&ylen,&Delx,&Dely,BC_Type,bcT,&Tini); //finding dimensions of problem domain

  //calculating ranks and vector lengths of data
  if (gsum == 1){ gscale = N_g; } //if decomposing data over all groups, must include length of groups
  else{ gscale = 1; }
  if ((dcmp_data == 1) || (dcmp_data == 2)){ mscale = N_m; } //if decomposing intensities, must include length of angular discretization
  else{ mscale = 1; }
  clen_avg = N_x*N_y*mscale*gscale; clen_edgV = (N_x+1)*N_y*mscale*gscale; clen_edgH = N_x*(N_y+1)*mscale*gscale; //clen='column length'
  rank_avg = min(clen_avg,N_t); rank_edgV = min(clen_edgV,N_t); rank_edgH = min(clen_edgH,N_t); //rank of each datatype

  if(dcmp_data == 0){ //if decomposing QD factors, must also include decomposition of boundary factors
    BClen = 2*(N_x+N_y)*N_g;
    rank_BC = min(BClen,N_t);
  }
  else{
    BClen = 0;
    rank_BC = 0;
  }

  /*===========================================================================*/
  /*                                                                            /
  /                          INITIALIZING OUTPUT FILE                           /
  /                                                                            */
  /*===========================================================================*/
  printf("Initializing output file\n");
  err = nc_create(outfile,NC_NETCDF4,&ncid_out); Handle_Err(err,loc); //opening NetCDF dataset

  err = Def_Dims(ncid_out,N_t,N_g,N_m,N_y,N_x,rank_BC,rank_avg,rank_edgV,rank_edgH,BClen,clen_avg,clen_edgV,clen_edgH,tlen,Delt,
    xlen,ylen,Delx,Dely,BC_Type,bcT,Tini,&N_t_ID,&N_g_ID,&N_m_ID,&N_y_ID,&N_x_ID,&rank_BC_ID,&rank_avg_ID,&rank_edgV_ID,
    &rank_edgH_ID,&BClen_ID,&clen_avg_ID,&clen_edgV_ID,&clen_edgH_ID);
  if (err != 0){ printf("Error occured while defining dimensions in %s\n",outfile); exit(1); }

  err = Def_DCMP_Vars(ncid_out,dcmp_type,dcmp_data,gsum,N_g_ID,rank_BC_ID,rank_avg_ID,rank_edgV_ID,rank_edgH_ID,BClen_ID,clen_avg_ID,
    clen_edgV_ID,clen_edgH_ID,N_t_ID,&BCg_IDs,&fg_avg_xx_IDs,&fg_edgV_xx_IDs,&fg_avg_yy_IDs,&fg_edgH_yy_IDs,&fg_edgV_xy_IDs,
    &fg_edgH_xy_IDs,&Ig_avg_IDs,&Ig_edgV_IDs,&Ig_edgH_IDs);

  if (err != 0){ printf("Error occured while defining variables in %s\n",outfile); exit(1); }

  err = nc_enddef(ncid_out); Handle_Err(err,loc); //changing output file to data mode

  /*===========================================================================*/
  /*                                                                            /
  /                          PERFORMING DECOMPOSITION                           /
  /                                                                            */
  /*===========================================================================*/
  err = Decompose_Data(ncid_in,ncid_out,dcmp_type,dcmp_data,gsum,N_t,N_g,N_m,N_x,N_y,BCg_IDs,fg_avg_xx_IDs,fg_edgV_xx_IDs,
    fg_avg_yy_IDs,fg_edgH_yy_IDs,fg_edgV_xy_IDs,fg_edgH_xy_IDs,Ig_avg_IDs,Ig_edgV_IDs,Ig_edgH_IDs);

  /*===========================================================================*/
  /*                                                                            /
  /                               END OF PROGRAM                                /
  /                                                                            */
  /*===========================================================================*/
  err = nc_close(ncid_in); Handle_Err(err,loc); //closing NetCDF dataset
  err = nc_close(ncid_out); Handle_Err(err,loc); //closing NetCDF dataset

  //successfull completion of program
  printf("Program finish\n");
  exit(0);

}
