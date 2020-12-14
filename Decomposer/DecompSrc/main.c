/*================================================================================================================================*/
/*
  Decomposer.c
  Author: Joseph M. Coale
  jmcoale@ncsu.edu - josephcoale912@gmail.com
*/
/*================================================================================================================================*/
/* Importing Packages */
/*================================================================================================================================*/
/* ----- EXTERNAL ----- */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/* ----- LOCAL ----- */
#include "Data_Handling.h"

/*================================================================================================================================*/
/* Importing Functions */
/*================================================================================================================================*/
/* ----- FROM NCDF_IO.c ----- */
void Handle_Err(const int Status, const char *Location);

/* ----- FROM INPUTS.c ----- */
int Input(const char *infile, char *dsfile, char *outfile, int *dcmp_type, int *gsum, double *svd_eps, Spec **Prb_specs, size_t *N_specs,
  Data **Dcmp_data, size_t *N_data, Data **Disc_Wts, size_t *N_wts);

void Get_Dims(const int ncid, int *BC_Type, Spec *Prb_specs, const size_t N_specs,
  Data *Dcmp_data, const size_t N_data, ncdim **dims, size_t *N_dims, Data *Disc_Wts, const size_t N_wts);

int Get_Disc(const int ncid, Data *Disc_Wts, const size_t N_wts, ncdim *dims);

/* ----- FROM OUTPUTS.c ----- */
int Def_Dims(const int ncid, Data *Dcmp_data, const size_t N_data, ncdim *dims, const size_t N_dims, ncdim **clen, ncdim **rank, const int gsum);

int Def_Disc(const int ncid, Data *Disc_Wts, const size_t N_wts, ncdim *dims);

int Def_Specs(const int ncid, Spec *Prb_specs, const size_t N_specs);

int Def_DCMP_Vars(const int ncid, const int dcmp_type, const int gsum, Data *Dcmp_data, const size_t N_data, ncdim *clen,
  ncdim *rank, ncdim *dims, Data **Decomp);

int Decompose_Data(const int ncid_in, const int ncid_out, const int dcmp_type, const int dcmp_data, const int gsum,
  const size_t N_t, const size_t N_g, const size_t N_m, const size_t N_x, const size_t N_y, const size_t rank_BC,
  const size_t rank_avg, const size_t rank_edgV, const size_t rank_edgH, const double svd_eps, const int *BCg_IDs,
  const int *fg_avg_xx_IDs, const int *fg_edgV_xx_IDs, const int *fg_avg_yy_IDs, const int *fg_edgH_yy_IDs,
  const int *fg_edgV_xy_IDs, const int *fg_edgH_xy_IDs, const int *Ig_avg_IDs,const int *Ig_edgV_IDs, const int *Ig_edgH_IDs);

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
  double svd_eps;

  //problem parameters
  // double Delt, *Delx, *Dely;
  int BC_Type[4];
  size_t N_t, N_g, N_m, N_y, N_x;
  size_t mscale, gscale, tscale;
  size_t rank_BC, rank_avg, rank_edgV, rank_edgH, BClen, clen_avg, clen_edgV, clen_edgH;
  int N_t_ID, N_g_ID;
  int rank_BC_ID, rank_avg_ID, rank_edgV_ID, rank_edgH_ID, BClen_ID, clen_avg_ID, clen_edgV_ID, clen_edgH_ID;

  //QDf decomposition IDs
  int *BCg_IDs, *fg_avg_xx_IDs, *fg_edgV_xx_IDs, *fg_avg_yy_IDs, *fg_edgH_yy_IDs, *fg_edgV_xy_IDs, *fg_edgH_xy_IDs;

  //I decomposition IDs
  int *Ig_avg_IDs, *Ig_edgV_IDs, *Ig_edgH_IDs;

  //
  size_t N_specs, N_data, N_dims, N_Wts;
  Spec *Prb_specs;
  Data *Dcmp_data, *Disc_Wts, *Decomp;
  ncdim *dims, *clen, *rank;

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
  err = Input(infile,dsfile,outfile,&dcmp_type,&gsum,&svd_eps,&Prb_specs,&N_specs,&Dcmp_data,&N_data,&Disc_Wts,&N_Wts);
  if (err != 0){
    printf("Error detected upon user input, aborting program\n");
    exit(1);
  }

  //opening, parsing datafile
  printf("Reading datafile\n");
  err = nc_open(dsfile,NC_NOWRITE,&ncid_in); Handle_Err(err,loc); //opening NetCDF dataset
  Get_Dims(ncid_in,BC_Type,Prb_specs,N_specs,Dcmp_data,N_data,&dims,&N_dims,Disc_Wts,N_Wts); //finding dimensions of problem domain
  err = Get_Disc(ncid_in,Disc_Wts,N_Wts,dims); //reading in discretization of domain

  /*===========================================================================*/
  /*                                                                            /
  /                          INITIALIZING OUTPUT FILE                           /
  /                                                                            */
  /*===========================================================================*/
  printf("Initializing output file\n");
  err = nc_create(outfile,NC_NETCDF4,&ncid_out); Handle_Err(err,loc); //opening NetCDF dataset

  err = Def_Dims(ncid_out,Dcmp_data,N_data,dims,N_dims,&clen,&rank,gsum); //defining dimensions in output file
  if (err != 0){ printf("Error occured while defining dimensions in %s\n",outfile); exit(1); }

  Def_DCMP_Vars(ncid_out, dcmp_type, gsum, Dcmp_data, N_data, clen, rank, dims, &Decomp);
  if (err != 0){ printf("Error occured while defining variables in %s\n",outfile); exit(1); }

  err = nc_enddef(ncid_out); Handle_Err(err,loc); //changing output file to data mode

  err =  Def_Disc(ncid_out,Disc_Wts,N_Wts,dims); //writing discretization to output file
  err =  Def_Specs(ncid_out,Prb_specs,N_specs); //writing problem specs to output file

  err = nc_close(ncid_in); Handle_Err(err,loc); //closing NetCDF dataset
  err = nc_close(ncid_out); Handle_Err(err,loc); //closing NetCDF dataset

  exit(0);

  /*===========================================================================*/
  /*===========================================================================*/
  //                  END OF RE-CONSTRUCTION SO FAR
  /*===========================================================================*/
  /*===========================================================================*/

  /*===========================================================================*/
  /*                                                                            /
  /                          PERFORMING DECOMPOSITION                           /
  /                                                                            */
  /*===========================================================================*/
  err = Decompose_Data(ncid_in,ncid_out,dcmp_type,dcmp_data,gsum,N_t,N_g,N_m,N_x,N_y,rank_BC,rank_avg,rank_edgV,rank_edgH,svd_eps,BCg_IDs,
    fg_avg_xx_IDs,fg_edgV_xx_IDs,fg_avg_yy_IDs,fg_edgH_yy_IDs,fg_edgV_xy_IDs,fg_edgH_xy_IDs,Ig_avg_IDs,Ig_edgV_IDs,Ig_edgH_IDs);

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
