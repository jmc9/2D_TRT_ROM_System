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

void Sort_Uniq_int(int *list, const size_t len, size_t *ulen);

/*================================================================================================================================*/
/**/
/*================================================================================================================================*/
int Read_List(char *line, char ***list_out, const size_t llen, const size_t lstrt)
{
  int err = 0;
  char delim = ' ';

  //list will hold entire input line for Prb_spec input
  char **list;
  list = (char **)malloc(sizeof(char*)*(llen+lstrt));
  for (size_t i=0; i<(llen+lstrt); i++){
    list[i] = (char *)malloc(sizeof(char)*20);
    memset(list[i],0,20);
  }

  //allocating array of list to output
  *list_out = (char **)malloc(sizeof(char*)*llen);
  for (size_t i=0; i<llen; i++){
    (*list_out)[i] = (char *)malloc(sizeof(char)*20);
    memset((*list_out)[i],0,20);
  }

  //copying list without leading elements into list_out array
  err = Delimit(line,delim,list,(int)(llen+lstrt));
  for (size_t i=0; i<llen; i++){
    strcpy((*list_out)[i],list[i+lstrt]);
  }

  for (size_t i=0; i<(llen+lstrt); i++){
    free(list[i]);
  }
  free(list);

  return err;
}

/*================================================================================================================================*/
/**/
/*================================================================================================================================*/
int Load_Specs(char **spec_names, size_t *N_specs, Spec **Prb_specs)
{
  int err = 0;

  //holding input number of specs, resetting N_specs
  size_t ns = *N_specs;
  *N_specs = 0;

  //finding updated number of specs if any default flags encountered
  for (size_t n=0; n<ns; n++){
    if (strcmp(spec_names[n],"TRT") == 0){ //TRT = preset flags for TRT problem
      *N_specs = *N_specs + 8;
    }
    else{
      *N_specs = *N_specs + 1;
    }
  }

  *Prb_specs = (Spec *)malloc(sizeof(Spec)*(*N_specs));

  //loading in all spec_names
  size_t p=0;
  for (size_t n=0; n<ns; n++){
    if (strcmp(spec_names[n],"TRT") == 0){ //TRT flag detected, loading default specs
      strcpy((*Prb_specs)[p].name,"tlen"); p++;
      strcpy((*Prb_specs)[p].name,"xlen"); p++;
      strcpy((*Prb_specs)[p].name,"ylen"); p++;
      strcpy((*Prb_specs)[p].name,"bcT_left"); p++;
      strcpy((*Prb_specs)[p].name,"bcT_bottom"); p++;
      strcpy((*Prb_specs)[p].name,"bcT_right"); p++;
      strcpy((*Prb_specs)[p].name,"bcT_top"); p++;
      strcpy((*Prb_specs)[p].name,"Tini"); p++;
    }
    else{
      strcpy((*Prb_specs)[p].name,spec_names[n]); p++;
    }
  }

  //freeing temporary array
  for (size_t n=0; n<ns; n++){
    free(spec_names[n]);
  }
  free(spec_names);

  //termination
  return err;
}

/*================================================================================================================================*/
/**/
/*================================================================================================================================*/
int Load_Data(char **data_names, size_t *N_data, Data **Dcmp_data)
{
  int err = 0;

  //holding input number of specs, resetting N_data
  size_t ns = *N_data;
  *N_data = 0;

  //finding updated number of specs if any default flags encountered
  for (size_t n=0; n<ns; n++){
    if (strcmp(data_names[n],"QDf") == 0){ //QDf = preset flag for TRT problem - QD factor decomposition
      *N_data = *N_data + 7;
    }
    else if (strcmp(data_names[n],"I") == 0){ //QDf = preset flag for TRT problem - QD factor decomposition
      *N_data = *N_data + 3;
    }
    else if (strcmp(data_names[n],"Mean_I") == 0){ //QDf = preset flag for TRT problem - QD factor decomposition
      *N_data = *N_data + 3;
    }
    else{
      *N_data = *N_data + 1;
    }
  }

  *Dcmp_data = (Data *)malloc(sizeof(Data)*(*N_data));

  //loading in all data_names
  size_t p=0;
  for (size_t n=0; n<ns; n++){
    if (strcmp(data_names[n],"QDf") == 0){ //TRT flag detected, loading default specs
      strcpy((*Dcmp_data)[p].name,"BCg");
      (*Dcmp_data)[p].opt[0] = 1;
      (*Dcmp_data)[p].opt[1] = 4;
      strcpy((*Dcmp_data)[p].cdat,"Cg_L,Cg_B,Cg_R,Cg_T");
      p++;

      strcpy((*Dcmp_data)[p].name,"fg_avg_xx");  (*Dcmp_data)[p].opt[0] = 0; p++;
      strcpy((*Dcmp_data)[p].name,"fg_edgV_xx"); (*Dcmp_data)[p].opt[0] = 0; p++;
      strcpy((*Dcmp_data)[p].name,"fg_avg_yy");  (*Dcmp_data)[p].opt[0] = 0; p++;
      strcpy((*Dcmp_data)[p].name,"fg_edgH_yy"); (*Dcmp_data)[p].opt[0] = 0; p++;
      strcpy((*Dcmp_data)[p].name,"fg_edgV_xy"); (*Dcmp_data)[p].opt[0] = 0; p++;
      strcpy((*Dcmp_data)[p].name,"fg_edgH_xy"); (*Dcmp_data)[p].opt[0] = 0; p++;
    }
    else if (strcmp(data_names[n],"I") == 0){ //TRT flag detected, loading default specs
      strcpy((*Dcmp_data)[p].name,"I_avg"); (*Dcmp_data)[p].opt[0] = 0; p++;
      strcpy((*Dcmp_data)[p].name,"I_edgV"); (*Dcmp_data)[p].opt[0] = 0; p++;
      strcpy((*Dcmp_data)[p].name,"I_edgH"); (*Dcmp_data)[p].opt[0] = 0; p++;
    }
    else{
      strcpy((*Dcmp_data)[p].name,data_names[n]); (*Dcmp_data)[p].opt[0] = 0; p++;
    }
  }

  //freeing temporary array
  for (size_t n=0; n<ns; n++){
    free(data_names[n]);
  }
  free(data_names);

  //termination
  return err;
}

/*================================================================================================================================*/
/**/
/*================================================================================================================================*/
int Load_Wts(char **wt_names, size_t *N_wts, Data **Disc_Wts)
{
  int err = 0;

  //holding input number of specs, resetting N_data
  size_t ns = *N_wts;
  *N_wts = 0;

  //finding updated number of specs if any default flags encountered
  for (size_t n=0; n<ns; n++){
    if (strcmp(wt_names[n],"default") == 0){ //default = preset flag for TRT problem
      *N_wts = *N_wts + 3;
    }
    else{
      *N_wts = *N_wts + 1;
    }
  }

  *Disc_Wts = (Data *)malloc(sizeof(Data)*(*N_wts));

  //loading in all data_names
  size_t p=0;
  for (size_t n=0; n<ns; n++){
    if (strcmp(wt_names[n],"default") == 0){ //TRT flag detected, loading default specs
      strcpy((*Disc_Wts)[p].name,"Delt"); (*Disc_Wts)[p].opt[0] = 0; p++;
      strcpy((*Disc_Wts)[p].name,"Delx"); (*Disc_Wts)[p].opt[0] = 0; p++;
      strcpy((*Disc_Wts)[p].name,"Dely"); (*Disc_Wts)[p].opt[0] = 0; p++;
    }
    else{
      strcpy((*Disc_Wts)[p].name,wt_names[n]); p++;
    }
  }

  //freeing temporary array
  for (size_t n=0; n<ns; n++){
    free(wt_names[n]);
  }
  free(wt_names);

  //termination
  return err;
}

/*================================================================================================================================*/
/**/
/*================================================================================================================================*/
int Read_Data_Dims(const int ncid, Data *data)
{
  int err, id;
  char loc[15] = "Read_Data_Dims";

  if (data->opt[0] == 0){ //opt = 0 --> decomposing the data as is
    err = nc_inq_varid(ncid, data->name, &data->id); Handle_Err(err,loc); //find data ID
    err = nc_inq_varndims(ncid, data->id, (int*)&data->ndims); Handle_Err(err,loc); //find number of dimensions
    data->dimids = malloc(sizeof(int)*(size_t)data->ndims); //allocate dimids as ndims
    err = nc_inq_vardimid(ncid, data->id, data->dimids); Handle_Err(err,loc); //read in dimension id's for each dimension

  }
  else if (data->opt[0] == 1){ //opt = 1 --> stack data into new constructed datamatrix to decompose
    size_t ndat = (size_t)data->opt[1]; //opt[1] holds the number of data matrices to to stack

    //allocating dats to hold names of construction data matrices
    char **dats = (char**)malloc(sizeof(char*)*ndat);
    for (size_t j=0; j<ndat; j++){
      dats[j] = (char*)malloc(sizeof(char)*50);
      memset(dats[j],0,50);
    }

    char del=',';
    err = Delimit(data->cdat,del,dats,(int)ndat); //finding names of construction data, assumed to be comma delimited

    err = nc_inq_varid(ncid, dats[0], &data->id); Handle_Err(err,loc); /*setting the data ID to the ID of the *first* construction data
                                                                       -- this is for the purposes of grid manipulation later on,
                                                                       -- since construction data for a dataset with opt[0]=1 are
                                                                       -- assumed to all reside on at least the same grid in time */

    int *nd = malloc(sizeof(int)*ndat); //nd holds the number of dimensions for each construction data
    int **ds = malloc(sizeof(int*)*ndat); //ds holds the dimensions of each construction data
    for (size_t j=0; j<ndat; j++){
      err = nc_inq_varid(ncid, dats[j], &id); Handle_Err(err,loc); //find data ID
      err = nc_inq_varndims(ncid, id, &nd[j]); Handle_Err(err,loc); //find number of dimensions
      ds[j] = (int*)malloc(sizeof(int)*(size_t)nd[j]); //allocate dimids as ndims
      err = nc_inq_vardimid(ncid, id, ds[j]); Handle_Err(err,loc); //read in dimension id's for each dimension

    }

    //making sure each construction data has the same dimensionality
    for (size_t j=1; j<ndat; j++){
      if (nd[j] != nd[j-1]){
        printf("Dimension counts of data to be stacked are not the same!");
        exit(1);
      }
    }

    //storing collective dimension id's for all construction matrices into data.dimids
    data->ndims = (size_t)nd[0] + ndat - 1;
    data->dimids = malloc(sizeof(int)*data->ndims);
    for (size_t j=0; j<data->ndims-ndat; j++){
      data->dimids[j] = ds[0][j];
    }
    size_t p = 0;
    for (size_t j=data->ndims-ndat; j<data->ndims; j++){
      data->dimids[j] = ds[p][nd[p]-1]; p++;
    }

    //freeing all allocated memory
    for (size_t j=0; j<ndat; j++){
      free(ds[j]);
    }
    free(ds);
    free(nd);

  }

  return err;
}

/*================================================================================================================================*/
/**/
/*================================================================================================================================*/
void Get_Dims(const int ncid, int *BC_Type, Spec *Prb_specs, const size_t N_specs,
  Data *Dcmp_data, const size_t N_data, ncdim **dims, size_t *N_dims, Data *Disc_Wts, const size_t N_wts)
{
  int err;
  char loc[9] = "Get_Dims";

  *N_dims = 0;

  //collecting dimensions of each decomp data
  for (size_t i=0; i<N_data; i++){

    err = Read_Data_Dims(ncid,&Dcmp_data[i]);
    *N_dims = *N_dims + Dcmp_data[i].ndims;

  }
  //collecting dimensions of each discretization weight
  for (size_t i=0; i<N_wts; i++){

    err = Read_Data_Dims(ncid,&Disc_Wts[i]);
    *N_dims = *N_dims + Disc_Wts[i].ndims;

  }

  int *dimlist = malloc(sizeof(int)*(*N_dims)); //allocating array to hold each dimension ID
  size_t p = 0; //counting variable

  //putting the dimensions of each Dcmp_data into dimlist
  for (size_t i=0; i<N_data; i++){
    for (size_t j=0; j<Dcmp_data[i].ndims; j++){
      dimlist[p] = Dcmp_data[i].dimids[j];
      p++;
    }
  }

  //putting the dimensions of each Disc_Wt into dimlist
  for (size_t i=0; i<N_wts; i++){
    for (size_t j=0; j<Disc_Wts[i].ndims; j++){
      dimlist[p] = Disc_Wts[i].dimids[j];
      p++;
    }
  }

  //sorting dimlist to have all unique dimensions in the first indexes
  //counting number of unique dimensions
  size_t N_dims_uniq; //number of unique dimensions
  Sort_Uniq_int(dimlist,*N_dims,&N_dims_uniq);

  //allocating dims as number of unique dimensions, collecting names and lengths
  *N_dims = N_dims_uniq;
  *dims = (ncdim*)malloc(sizeof(ncdim)*(*N_dims));
  for (size_t i=0; i<*N_dims; i++){
    (*dims)[i].id = dimlist[i];
    memset((*dims)[i].name,0,10);
    err = nc_inq_dimname(ncid,(*dims)[i].id,(*dims)[i].name); Handle_Err(err,loc);
    err = nc_inq_dimlen(ncid,(*dims)[i].id,&(*dims)[i].len); Handle_Err(err,loc);

  }
  free(dimlist);

  //replacing Dcmp_data.dimids with the location in the dims array corresponding to each dimid
  for (size_t j=0; j<N_data; j++){
    for (size_t k=0; k<Dcmp_data[j].ndims; k++){
      for (size_t i=0; i<*N_dims; i++){
        if (Dcmp_data[j].dimids[k] == (*dims)[i].id){
          Dcmp_data[j].dimids[k] = (int)i;
          break;
        }
      }
    }
  }

  //replacing Disc_Wts.dimids with the location in the dims array corresponding to each dimid
  for (size_t j=0; j<N_wts; j++){
    for (size_t k=0; k<Disc_Wts[j].ndims; k++){
      for (size_t i=0; i<*N_dims; i++){
        if (Disc_Wts[j].dimids[k] == (*dims)[i].id){
          Disc_Wts[j].dimids[k] = (int)i;
          break;
        }
      }
    }
  }

  //reading in problem specs
  for (size_t i=0; i<N_specs; i++){
    err = Get_Spec(ncid,&Prb_specs[i]);
  }

  // err = nc_get_att_int(ncid,NC_GLOBAL,"BC_type",BC_Type); Handle_Err(err,loc);

}

/*================================================================================================================================*/
/**/
/*================================================================================================================================*/
int Get_Grids(const int ncid, Data *Dcmp_data, const size_t N_data, Data **Grids, size_t *N_grids, ncdim *dims, const size_t N_dims)
{
  int err;
  char loc[10] = "Get_Grids";

  size_t N_grids_ = 0;
  //collecting the grids over which each each dataset is defined over
  for (size_t i=0; i<N_data; i++){

      err = nc_get_att_int(ncid, Dcmp_data[i].id, "N_grids", (int*)&Dcmp_data[i].ngrids); //finding the number of grids defined for dataset-i
      Dcmp_data[i].gridids = (int*)malloc(sizeof(int)*Dcmp_data[i].ngrids); //allocating array to hold the nc-ID's of each grid
      N_grids_ = N_grids_ + Dcmp_data[i].ngrids; //adding the number of grids defined for dataset-i to the number of all grids
      for (size_t j=0; j<Dcmp_data[i].ngrids; j++){
        char buf[10], gname[50];
        memset(buf,0,10); memset(gname,0,50);
        sprintf(buf,"grid%ld",j);
        err = nc_get_att_text(ncid, Dcmp_data[i].id, buf, gname); Handle_Err(err,loc); //getting the name of the jth grid for dataset-i
        err = nc_inq_varid(ncid, gname, &Dcmp_data[i].gridids[j]); Handle_Err(err,loc); //getting the nc-ID of the jth grid given its name
      }

  }

  //putting the nc_ID's of each grid collected above into one list
  int *gridlist = malloc(sizeof(Data)*N_grids_);
  size_t p = 0;
  for (size_t i=0; i<N_data; i++){
    for (size_t j=0; j<Dcmp_data[i].ngrids; j++){
      gridlist[p] = Dcmp_data[i].gridids[j];
      p++;
    }
  }

  //finding the number of *unique* grids by nc-ID
  Sort_Uniq_int(gridlist,N_grids_,N_grids);

  //creating an array of grid structs (Data type) to carry information of each unique grid
  *Grids = (Data*)malloc(sizeof(Data)*(*N_grids));
  for (size_t i=0; i<(*N_grids); i++){
    (*Grids)[i].id = gridlist[i];
    memset((*Grids)[i].name,0,10);
    err = nc_inq_varname(ncid, (*Grids)[i].id, (*Grids)[i].name); Handle_Err(err,loc);
    err = nc_inq_varndims(ncid, (*Grids)[i].id, (int*)&(*Grids)[i].ndims); Handle_Err(err,loc);
    (*Grids)[i].dimids = (int*)malloc(sizeof(int)*(*Grids)[i].ndims);
    err = nc_inq_vardimid(ncid, (*Grids)[i].id, (*Grids)[i].dimids); Handle_Err(err,loc);

    for (size_t j=0; j<(*Grids)[i].ndims; j++){
      for (size_t k=0; k<N_dims; k++){
        if ((*Grids)[i].dimids[j] == dims[k].id){
          (*Grids)[i].dimids[j] = (int)k;
          break;
        }
      }
    }

    size_t dlen = 1.;
    for (size_t j=0; j<(*Grids)[i].ndims; j++){
      dlen = dlen*dims[(*Grids)[i].dimids[j]].len;
    }
    (*Grids)[i].dat = (double*)malloc(sizeof(double)*dlen);
    err = nc_get_var_double(ncid, (*Grids)[i].id, &(*Grids)[i].dat[0]); Handle_Err(err,loc);

    (*Grids)[i].bnds = (double*)malloc(sizeof(double)*2);
    err = nc_get_att_double(ncid, (*Grids)[i].id, "bnds", (*Grids)[i].bnds); Handle_Err(err,loc);
  }
  free(gridlist);

  //for each dataset, replacing the nc-ID's of its grids with location of each grid in the 'Grids' array
  for (size_t i=0; i<N_data; i++){
    for (size_t j=0; j<Dcmp_data[i].ngrids; j++){
      for (size_t k=0; k<(*N_grids); k++){
        if (Dcmp_data[i].gridids[j] == (*Grids)[k].id){
          Dcmp_data[i].gridids[j] = (int)k;
          break;
        }
      }
    }
  }

  return err;
}

/*================================================================================================================================*/
/**/
/*================================================================================================================================*/
int Get_Disc(const int ncid, Data *Disc_Wts, const size_t N_wts, ncdim *dims)
{
  int err;
  char loc[9] = "Get_Disc";
  for (size_t i=0; i<N_wts; i++){
    if (Disc_Wts[i].ndims == 0){
      Disc_Wts[i].dat = (double*)malloc(sizeof(double));
      err = nc_get_var_double(ncid,Disc_Wts[i].id,Disc_Wts[i].dat); Handle_Err(err,loc);
    }
    else{
      size_t dsum = 0.;
      for (size_t j=0; j<Disc_Wts[i].ndims; j++){
        dsum = dsum + dims[Disc_Wts[i].dimids[j]].len;
      }
      Disc_Wts[i].dat = (double*)malloc(sizeof(double)*dsum);
      err = nc_get_var_double(ncid,Disc_Wts[i].id,Disc_Wts[i].dat); Handle_Err(err,loc);
    }
  }
  return err;
}

/*================================================================================================================================*/
/*Input is used to load user input from a specified input file
  infile - name of input file
  dsfile - name of datafile holding data to decompose
  outfile - name of file to output decompositon to
  dcmp_type - type of decomposition to use
  dcmp_data - type of data to decompose */
/*================================================================================================================================*/
int Input(const char *infile, char *dsfile, char *outfile, int *dcmp_type, int *gsum, double *svd_eps, Spec **Prb_specs, size_t *N_specs,
  Data **Dcmp_data, size_t *N_data, Data **Disc_Wts, size_t *N_wts, int *base_subtract)
{
  FILE *inpf;
  char line[256];
  char **inps, delim = ' ';
  char dcmp_type_in[10], **spec_names, **data_names, **wt_names;
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
  *gsum = 1;
  *svd_eps = -1.;
  *N_specs = 0;
  *base_subtract = 0;
  int data_check = 0;
  int wts_check = 0;


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
    else if (strcmp(inps[0],"svd_eps") == 0){
      sscanf(inps[1], "%le", svd_eps);
    }
    else if (strcmp(inps[0],"dcmp_data") == 0){

      //finding number of problem specs to read in
      int n_data;
      sscanf(inps[1], "%d", &n_data);
      if (n_data <= 0){ //checking for valid N_data
        printf("Invalid N_data = %ld! must be > 0\n",*N_data);
        return 1;
      }
      *N_data = (size_t)n_data;

      //reading in list of data_names
      if (*N_data > 0){
        err = Read_List(line,&data_names,*N_data,2);
      }

      //setting flag that data was found on input
      data_check = 1;

    }
    else if (strcmp(inps[0],"Prb_spec") == 0){

      //finding number of problem specs to read in
      int n_specs;
      sscanf(inps[1], "%d", &n_specs);
      if (n_specs < 0){ //checking for valid N_specs
        printf("Invalid Prb_spec = %ld! must be >= 0\n",*N_specs);
        return 1;
      }
      *N_specs = (size_t)n_specs;

      //reading in list of spec_names
      if (*N_specs > 0){
        err = Read_List(line,&spec_names,*N_specs,2);
      }

    }
    else if (strcmp(inps[0],"Disc_wts") == 0){

      //finding number of problem specs to read in
      int n_wts;
      sscanf(inps[1], "%d", &n_wts);
      if (n_wts < 0){ //checking for valid N_specs
        printf("Invalid Disc_wts = %ld! must be >= 0\n",*N_wts);
        return 1;
      }
      *N_wts = (size_t)n_wts;

      //reading in list of wt_names
      if (*N_wts > 0){
        err = Read_List(line,&wt_names,*N_wts,2);
      }

      //setting flag that Disc_wts was found on input
      wts_check = 1;

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
  else if (strcmp(dcmp_type_in,"DMDb") == 0){
    *dcmp_type = 1;
    *base_subtract = 1;
  }
  else{
    printf("Error! [location INPUTS.c/Input]\n");
    printf("Unrecognized dcmp_type in input file! (%s)\n",dcmp_type_in);
    printf("Valid dcmp_type include: POD, PODg, DMD\n");
    return 1;
  }

  if (data_check == 0){
    *N_data = 1;
    data_names = (char **)malloc(sizeof(char *));
    data_names[0] = (char *)malloc(sizeof(char)*20);
    memset(data_names[0],0,20);
    strcpy(data_names[0],"QDf");
  }
  if (wts_check == 0){
    *N_wts = 1;
    wt_names = (char **)malloc(sizeof(char *));
    wt_names[0] = (char *)malloc(sizeof(char)*20);
    memset(wt_names[0],0,20);
    strcpy(wt_names[0],"default");
  }

  //loading Prb_specs (checking for flags to 'default' sets of spec_names)
  if (*N_specs > 0){
    err = Load_Specs(spec_names,N_specs,Prb_specs);
  }

  //loading Dcmp_data
  if (*N_data > 0){
    err = Load_Data(data_names,N_data,Dcmp_data);
  }

  //loading Disc_Wts
  if (*N_wts > 0){
    err = Load_Wts(wt_names,N_wts,Disc_Wts);
  }

  //successfull termination
  return 0;

}
