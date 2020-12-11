/*================================================================================================================================*/
/*
  INPUTS.c
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

/* ----- LOCAL ----- */
#include "Data_Handling.h"

/*================================================================================================================================*/
/* Importing Functions */
/*================================================================================================================================*/
/* ----- FROM NCDF_IO.c ----- */
void Handle_Err(const int Status, const char *Location);

int Get_Spec(const int ncid, Spec *spec);

/* ----- FROM MISC_PROCS.c ----- */
int Delimit(const char *Line, const char Del, char **Parts, const int Nparts);

int Load_Specs(char ***spec_names, size_t *N_specs);

/*================================================================================================================================*/
/**/
/*================================================================================================================================*/
void Get_Dims(const int ncid, size_t *N_t, size_t *N_g, size_t *N_m, size_t *N_y, size_t *N_x, double *Delt,
  double **Delx, double **Dely, int *BC_Type, Spec *Prb_specs, const size_t N_specs)
{
  int err;
  int N_t_ID, N_g_ID, N_m_ID, N_y_ID, N_x_ID, Delt_ID, Delx_ID, Dely_ID;
  char loc[8] = "Get_Dims";

  //reading in dimension ID's
  err = nc_inq_dimid(ncid,"N_t",&N_t_ID); Handle_Err(err,loc);
  err = nc_inq_dimid(ncid,"N_g",&N_g_ID); Handle_Err(err,loc);
  err = nc_inq_dimid(ncid,"N_m",&N_m_ID); Handle_Err(err,loc);
  err = nc_inq_dimid(ncid,"N_y",&N_y_ID); Handle_Err(err,loc);
  err = nc_inq_dimid(ncid,"N_x",&N_x_ID); Handle_Err(err,loc);

  for (size_t i=0; i<N_specs; i++){
    err = Get_Spec(ncid,&Prb_specs[i]);
  }

  err = nc_inq_varid(ncid,"Delt",&Delt_ID); Handle_Err(err,loc);
  err = nc_inq_varid(ncid,"Delx",&Delx_ID); Handle_Err(err,loc);
  err = nc_inq_varid(ncid,"Dely",&Dely_ID); Handle_Err(err,loc);

  //reading in dimension values
  err = nc_inq_dimlen(ncid,N_t_ID,N_t); Handle_Err(err,loc);
  err = nc_inq_dimlen(ncid,N_g_ID,N_g); Handle_Err(err,loc);
  err = nc_inq_dimlen(ncid,N_m_ID,N_m); Handle_Err(err,loc);
  err = nc_inq_dimlen(ncid,N_y_ID,N_y); Handle_Err(err,loc);
  err = nc_inq_dimlen(ncid,N_x_ID,N_x); Handle_Err(err,loc);

  *Delx = (double *)malloc(sizeof(double)*(*N_x));
  *Dely = (double *)malloc(sizeof(double)*(*N_y));
  err = nc_get_var_double(ncid,Delx_ID,Delx[0]); Handle_Err(err,loc);
  err = nc_get_var_double(ncid,Dely_ID,Dely[0]); Handle_Err(err,loc);
  err = nc_get_var_double(ncid,Delt_ID,Delt); Handle_Err(err,loc);

  err = nc_get_att_int(ncid,NC_GLOBAL,"BC_type",BC_Type); Handle_Err(err,loc);

}

/*================================================================================================================================*/
/*Input is used to load user input from a specified input file
  infile - name of input file
  dsfile - name of datafile holding data to decompose
  outfile - name of file to output decompositon to
  dcmp_type - type of decomposition to use
  dcmp_data - type of data to decompose */
/*================================================================================================================================*/
int Input(const char *infile, char *dsfile, char *outfile, int *dcmp_type, int *dcmp_data, int *gsum, double *svd_eps, char ***spec_names, size_t *N_specs)
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
  *svd_eps = -1.;
  *N_specs = 1;

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
    else if (strcmp(inps[0],"svd_eps") == 0){
      sscanf(inps[1], "%le", svd_eps);
    }
    else if (strcmp(inps[0],"Prb_spec") == 0){

      //finding number of problem specs to read in
      sscanf(inps[1], "%ld", N_specs);
      if (*N_specs <= 0){ //checking for valid N_specs
        printf("Invalid Prb_spec = %ld! must be > 0\n",*N_specs);
        return 1;
      }

      //Prb_spec will hold entire input line for Prb_spec input
      char **Prb_spec;
      Prb_spec = (char **)malloc(sizeof(char*)*(*N_specs+2));
      for (size_t i=0; i<(*N_specs+2); i++){
        Prb_spec[i] = (char *)malloc(sizeof(char)*20);
        memset(Prb_spec[i],0,20);
      }

      //allocating array of spec names to output
      *spec_names = (char **)malloc(sizeof(char*)*(*N_specs));
      for (size_t i=0; i<*N_specs; i++){
        (*spec_names)[i] = (char *)malloc(sizeof(char)*20);
        memset((*spec_names)[i],0,20);
      }

      //copying names of problem specs into spec_names array
      err = Delimit(line,delim,Prb_spec,*(int*)N_specs+2);
      for (size_t i=0; i<*N_specs; i++){
        strcpy((*spec_names)[i],Prb_spec[i+2]);
      }

      for (size_t i=0; i<(*N_specs+2); i++){
        free(Prb_spec[i]);
      }
      free(Prb_spec);
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
  else if (strcmp(dcmp_type_in,"DMDg") == 0){
    *dcmp_type = 1;
    *gsum = 0;
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

  //loading spec_names (checking for flags to 'default' sets of spec_names)
  err = Load_Specs(spec_names,N_specs);

  //successfull termination
  return 0;

}

/*================================================================================================================================*/
/**/
/*================================================================================================================================*/
int Load_Specs(char ***spec_names, size_t *N_specs)
{
  int err = 0;

  char **temp; //temporary array
  temp = (char **)malloc(sizeof(char*)*(*N_specs));

  //copying spec_names array to temporary array and freeing spec_names
  for (size_t i=0; i<(*N_specs); i++){
    temp[i] = (char *)malloc(sizeof(char)*20); //allocate string
    memset(temp[i],0,20); //scrub memory
    strcpy(temp[i],(*spec_names)[i]); //copy data
    free((*spec_names)[i]); //free string
  }
  free(*spec_names); //free pointer array

  //holding input number of specs, resetting N_specs
  size_t ns = *N_specs;
  *N_specs = 0;

  //finding updated number of specs if any default flags encountered
  for (size_t n=0; n<ns; n++){
    if (strcmp(temp[n],"TRT") == 0){ //TRT = preset flags for TRT problem
      *N_specs = *N_specs + 8;
    }
    else{
      *N_specs = *N_specs + 1;
    }
  }

  //re-allocating spec_names array according to updated N_specs
  *spec_names = (char **)malloc(sizeof(char*)*(*N_specs));
  for (size_t i=0; i<(*N_specs); i++){
    (*spec_names)[i] = (char *)malloc(sizeof(char)*20);
    memset((*spec_names)[i],0,20);
  }

  //loading in all spec_names
  size_t p=0;
  for (size_t n=0; n<ns; n++){
    if (strcmp(temp[n],"TRT") == 0){ //TRT flag detected, loading default specs
      strcpy((*spec_names)[p],"tlen"); p++;
      strcpy((*spec_names)[p],"xlen"); p++;
      strcpy((*spec_names)[p],"ylen"); p++;
      strcpy((*spec_names)[p],"bcT_left"); p++;
      strcpy((*spec_names)[p],"bcT_bottom"); p++;
      strcpy((*spec_names)[p],"bcT_right"); p++;
      strcpy((*spec_names)[p],"bcT_top"); p++;
      strcpy((*spec_names)[p],"Tini"); p++;
    }
    else{
      strcpy((*spec_names)[p],temp[n]); p++;
    }
  }

  //freeing temporary array
  for (size_t n=0; n<ns; n++){
    free(temp[n]);
  }
  free(temp);

  //termination
  return err;
}
