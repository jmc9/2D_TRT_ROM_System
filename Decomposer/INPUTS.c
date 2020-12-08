/*================================================================================================================================*/
/*
  INPUTS.c
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

/*================================================================================================================================*/
/* Importing Functions */
/*================================================================================================================================*/
/* ----- FROM NCDF_IO.c ----- */
void Handle_Err(const int Status, const char *Location);

/* ----- FROM MISC_PROCS.c ----- */
int Delimit(const char *Line, const char Del, char **Parts, const int Nparts);

/*================================================================================================================================*/
/**/
/*================================================================================================================================*/
void Get_Dims(const int ncid, size_t *N_t, size_t *N_g, size_t *N_m, size_t *N_y, size_t *N_x, double *tlen, double *Delt,
  double *xlen, double *ylen, double **Delx, double **Dely, int *BC_Type, double *bcT, double *Tini)
{
  int err;
  int N_t_ID, N_g_ID, N_m_ID, N_y_ID, N_x_ID, tlen_ID, Delt_ID, xlen_ID, ylen_ID, Delx_ID, Dely_ID, bcT_ID[4], Tini_ID;
  char loc[8] = "Get_Dims";

  //reading in dimension ID's
  err = nc_inq_dimid(ncid,"N_t",&N_t_ID); Handle_Err(err,loc);
  err = nc_inq_dimid(ncid,"N_g",&N_g_ID); Handle_Err(err,loc);
  err = nc_inq_dimid(ncid,"N_m",&N_m_ID); Handle_Err(err,loc);
  err = nc_inq_dimid(ncid,"N_y",&N_y_ID); Handle_Err(err,loc);
  err = nc_inq_dimid(ncid,"N_x",&N_x_ID); Handle_Err(err,loc);

  err = nc_inq_varid(ncid,"tlen",&tlen_ID); Handle_Err(err,loc);
  err = nc_inq_varid(ncid,"Delt",&Delt_ID); Handle_Err(err,loc);
  err = nc_inq_varid(ncid,"xlen",&xlen_ID); Handle_Err(err,loc);
  err = nc_inq_varid(ncid,"ylen",&ylen_ID); Handle_Err(err,loc);
  err = nc_inq_varid(ncid,"Delx",&Delx_ID); Handle_Err(err,loc);
  err = nc_inq_varid(ncid,"Dely",&Dely_ID); Handle_Err(err,loc);
  err = nc_inq_varid(ncid,"bcT_left",&bcT_ID[0]); Handle_Err(err,loc);
  err = nc_inq_varid(ncid,"bcT_bottom",&bcT_ID[1]); Handle_Err(err,loc);
  err = nc_inq_varid(ncid,"bcT_right",&bcT_ID[2]); Handle_Err(err,loc);
  err = nc_inq_varid(ncid,"bcT_top",&bcT_ID[3]); Handle_Err(err,loc);
  err = nc_inq_varid(ncid,"Tini",&Tini_ID); Handle_Err(err,loc);

  //reading in dimension values
  err = nc_inq_dimlen(ncid,N_t_ID,N_t); Handle_Err(err,loc);
  err = nc_inq_dimlen(ncid,N_g_ID,N_g); Handle_Err(err,loc);
  err = nc_inq_dimlen(ncid,N_m_ID,N_m); Handle_Err(err,loc);
  err = nc_inq_dimlen(ncid,N_y_ID,N_y); Handle_Err(err,loc);
  err = nc_inq_dimlen(ncid,N_x_ID,N_x); Handle_Err(err,loc);

  err = nc_get_var_double(ncid,tlen_ID,tlen); Handle_Err(err,loc);
  err = nc_get_var_double(ncid,Delt_ID,Delt); Handle_Err(err,loc);
  err = nc_get_var_double(ncid,xlen_ID,xlen); Handle_Err(err,loc);
  err = nc_get_var_double(ncid,ylen_ID,ylen); Handle_Err(err,loc);

  *Delx = (double *)malloc(sizeof(double)*(*N_x));
  *Dely = (double *)malloc(sizeof(double)*(*N_y));
  err = nc_get_var_double(ncid,Delx_ID,Delx[0]); Handle_Err(err,loc);
  err = nc_get_var_double(ncid,Dely_ID,Dely[0]); Handle_Err(err,loc);

  err = nc_get_att_int(ncid,NC_GLOBAL,"BC_type",BC_Type); Handle_Err(err,loc);
  err = nc_get_var_double(ncid,bcT_ID[0],&bcT[0]); Handle_Err(err,loc);
  err = nc_get_var_double(ncid,bcT_ID[1],&bcT[1]); Handle_Err(err,loc);
  err = nc_get_var_double(ncid,bcT_ID[2],&bcT[2]); Handle_Err(err,loc);
  err = nc_get_var_double(ncid,bcT_ID[3],&bcT[3]); Handle_Err(err,loc);

  err = nc_get_var_double(ncid,Tini_ID,Tini); Handle_Err(err,loc);

}

/*================================================================================================================================*/
/*Input is used to load user input from a specified input file
  infile - name of input file
  dsfile - name of datafile holding data to decompose
  outfile - name of file to output decompositon to
  dcmp_type - type of decomposition to use
  dcmp_data - type of data to decompose */
/*================================================================================================================================*/
int Input(const char *infile, char *dsfile, char *outfile, int *dcmp_type, int *dcmp_data, int *gsum)
{
  FILE *inpf;
  char line[256];
  char **inps, delim = ' ';
  char dcmp_type_in[10], dcmp_data_in[10];
  int err;
  size_t nargs = 2;
  size_t inp_size = 25;

  inpf = fopen(infile,"r"); //opening input file
  //check if input file exists, if not terminate program
  if (inpf == NULL){
    printf("Error occured whilst opening file: %s\n",infile);
    return 1;
  }

  //allocating 2D array that holds input arguments
  inps = (char **)malloc(sizeof(char *)*nargs);
  for (size_t n=0; n < nargs; n++){
    inps[n] = (char *)malloc(sizeof(char)*inp_size);
  }

  //setting defaults
  strcpy(dsfile,"output.h5");
  strcpy(outfile,"DCMPout.h5");
  strcpy(dcmp_type_in,"POD");
  strcpy(dcmp_data_in,"QDf");
  *gsum = 1;

  //moving line by line through file to grab inputs
  while ( fgets(line, 255, inpf) != NULL ){ //placing current line of input file into 'line' character array

    for (size_t n=0; n < nargs; n++){
      for (size_t i=0; i < inp_size; i++){
        inps[n][i] = 0;
      }
    }

    //delimiting input line into seperate arguments
    err = Delimit(line,delim,inps,(int)nargs);

    //checking for valid input flags, storing associated arguments
    if (strcmp(inps[0],"dataset") == 0){
      memset(dsfile,0,10);
      strcpy(dsfile,inps[1]);
    }
    else if (strcmp(inps[0],"outfile") == 0){
      memset(outfile,0,10);
      strcpy(outfile,inps[1]);
    }
    else if (strcmp(inps[0],"dcmp_type") == 0){
      memset(dcmp_type_in,0,10);
      strcpy(dcmp_type_in,inps[1]);
    }
    else if (strcmp(inps[0],"dcmp_data") == 0){
      memset(dcmp_data_in,0,10);
      strcpy(dcmp_data_in,inps[1]);
    }

  }

  //deallocating input array
  free(inps);

  //closing input file
  err = fclose(inpf);
  if (err != 0){
    printf("Error occured whilst closing file: %s",infile);
    return 1;
  }

  //translating input dcmp_type to an integer value
  if (strcmp(dcmp_type_in,"POD") == 0){
    *dcmp_type = 0;
  }
  else if (strcmp(dcmp_type_in,"PODg") == 0){
    *dcmp_type = 0;
    *gsum = 0;
  }
  else if (strcmp(dcmp_type_in,"DMD") == 0){
    *dcmp_type = 1;
  }
  else{
    printf("Error! [location INPUTS.c/Input]\n");
    printf("Unrecognized dcmp_type in input file! (%s)\n",dcmp_type_in);
    printf("Valid dcmp_type include: POD, PODg, DMD\n");
    return 1;
  }

  //translating input dcmp_data to an integer value
  if (strcmp(dcmp_data_in,"QDf") == 0){
    (*dcmp_data) = 0;
  }
  else if (strcmp(dcmp_data_in,"I") == 0){
    (*dcmp_data) = 1;
  }
  else if (strcmp(dcmp_data_in,"Mean_I") == 0){
    (*dcmp_data) = 2;
  }
  else{
    printf("Error! [location INPUTS.c/Input]\n");
    printf("Unrecognized dcmp_data in input file! (%s)\n",dcmp_data_in);
    printf("Valid dcmp_data include: QDf, I, Mean_I\n");
    return 1;
  }

  //successfull termination
  return 0;

}
