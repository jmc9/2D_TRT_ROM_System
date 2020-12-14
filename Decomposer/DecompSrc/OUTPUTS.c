/*================================================================================================================================*/
/*
  OUTPUTS.c
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

void Get_Var_Double(const int ncid, const char *name, double **var, const size_t size);

int Write_Spec(const int ncid, Spec *spec);

/* ----- FROM CPOD_ROUTINES.c ----- */
int Generate_POD(const double *data, const int ncid_out, const char *dname, const size_t N_t, const size_t N_g, const size_t clen,
  const size_t rank, const int *DCMP_IDs);

int Generate_DMD(const double *data, const int ncid_out, const char *dname, const size_t N_t, const size_t N_g, const size_t clen,
  const size_t rank, const double svd_eps, const int *DCMP_IDs);

/* ----- FROM MISC_PROCS.c ----- */
void Sort_Uniq_sizet(size_t *list, const size_t len, size_t *ulen);

/* ----- LOCAL DEFINITIONS ----- */
#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })

/*================================================================================================================================*/
/**/
/*================================================================================================================================*/
int Def_Dims(const int ncid, Data *Dcmp_data, const size_t N_data, ncdim *dims, const size_t N_dims, ncdim **clen, ncdim **rank, const int gsum)
{
  int err;
  char loc[9] = "Def_Dims";

  size_t *clen_ = malloc(sizeof(size_t)*N_data); //array holding clen of each variable

  //calculating clen of each variable
  size_t k;
  for (size_t i=0; i<N_data; i++){

    if (Dcmp_data[i].opt[0] == 0){ //if decomposing the data as is, clen is simply the product of its dimensions (excluding time)
      if (gsum == 1){ k = 1; } //full phase-space decomposition: all dimensions except time go into clen
      else{ k = 2; } //groupwise decomposition: all dimensions except time and leading dimension go into clen
      clen_[i] = 1.; //initializing clen_ to 1

      for (size_t j=k; j<Dcmp_data[i].ndims; j++){ //taking product of dimension lengths
        clen_[i] = clen_[i]*dims[Dcmp_data[i].dimids[j]].len;
      }

    }
    else if (Dcmp_data[i].opt[0] == 1){ //if decomposing stacked data, clen depends on the summation of their dimensions
      clen_[i] = 0.; //initializing clen_ to 0
      for (size_t j=2; j<Dcmp_data[i].ndims; j++){ //summing the dimensions of stacked data (excluding time and leading dimensions)
        clen_[i] = clen_[i] + dims[Dcmp_data[i].dimids[j]].len;
      }

      if (gsum == 1){ //for full phase space decomposition, multiply stacked length with that of the leading dimension
        clen_[i] = clen_[i]*dims[Dcmp_data[i].dimids[1]].len;
      }

    }

    //prepping to replace variable dimensions with clen/rank
    free(Dcmp_data[i].dimids); //freeing previously malloc'd dimid arrays
    Dcmp_data[i].ndims = 2; //only 2 dimensions per dataset: clen, rank
    Dcmp_data[i].dimids = (int*)malloc(sizeof(int)*2);
    Dcmp_data[i].dimids[0] = (int)clen_[i];
  }

  size_t N_clen; //N_clen = number of unique 'clen' dimensions
  Sort_Uniq_sizet(clen_,N_data,&N_clen); //sort clen_ to have unique values as the first N_clen elements

  //allocating memory for clen, rank dimension lists
  *clen = (ncdim*)malloc(sizeof(ncdim)*N_clen);
  *rank = (ncdim*)malloc(sizeof(ncdim)*N_clen);

  //filling clen, rank dimension descriptions
  for (size_t i=0; i<N_clen; i++){
    (*clen)[i].len = clen_[i];
    (*rank)[i].len = min((*clen)[i].len,dims[0].len);

    char buf[10];
    memset(buf,0,10);
    sprintf(buf,"clen%ld",i);
    strcpy((*clen)[i].name,buf);

    memset(buf,0,10);
    sprintf(buf,"rank%ld",i);
    strcpy((*rank)[i].name,buf);

  }
  free(clen_);

  //mapping each dataset to its respective clen, rank dimensions
  for (size_t i=0; i<N_data; i++){
    for (size_t j=0; j<N_clen; j++){
      if (Dcmp_data[i].dimids[0] == (int)(*clen)[j].len){
        Dcmp_data[i].dimids[0] = (int)j;
        Dcmp_data[i].dimids[1] = (int)j;
        
        break;
      }
    }
  }

  //writing original dimensions (from input datafile) to the output ncdf file and replacing id's from input file
  for (size_t i=0; i<N_dims; i++){
    err = nc_def_dim(ncid,dims[i].name,dims[i].len,&dims[i].id); Handle_Err(err,loc);
  }

  //writing clen, rank dimensions
  for (size_t i=0; i<N_clen; i++){
    err = nc_def_dim(ncid,(*clen)[i].name,(*clen)[i].len,&(*clen)[i].id); Handle_Err(err,loc);
    err = nc_def_dim(ncid,(*rank)[i].name,(*rank)[i].len,&(*rank)[i].id); Handle_Err(err,loc);
  }

  //termination
  return err;
}

/*================================================================================================================================*/
/**/
/*================================================================================================================================*/
int Def_Disc(const int ncid, Data *Disc_Wts, const size_t N_wts, ncdim *dims)
{
  int err;
  char loc[9] = "Def_Disc";

  for (size_t i=0; i<N_wts; i++){

    if (Disc_Wts[i].ndims == 0){
      err = nc_def_var(ncid,Disc_Wts[i].name,NC_DOUBLE,0,0,&Disc_Wts[i].id); Handle_Err(err,loc);
    }
    else{
      int *d = malloc(sizeof(int)*Disc_Wts[i].ndims);
      for (size_t j=0; j<Disc_Wts[i].ndims; j++){
        d[j] = dims[Disc_Wts[i].dimids[j]].id;
      }
      err = nc_def_var(ncid,Disc_Wts[i].name,NC_DOUBLE,(int)Disc_Wts[i].ndims,d,&Disc_Wts[i].id); Handle_Err(err,loc);
      free(d);
    }

    err = nc_put_var_double(ncid,Disc_Wts[i].id,Disc_Wts[i].dat); Handle_Err(err,loc);

  }

  return err;
}

/*================================================================================================================================*/
/**/
/*================================================================================================================================*/
int Def_Specs(const int ncid, Spec *Prb_specs, const size_t N_specs)
{
  int err;

  for (size_t i=0; i<N_specs; i++){
    err = Write_Spec(ncid,&Prb_specs[i]);
  }

  return err;
}

/*================================================================================================================================*/
/**/
/*================================================================================================================================*/
int Def_POD_Vars(const int ncid, char *dname, const int rank_ID, const int clen_ID, const int N_t_ID, const int N_g_ID,
  const int gsum, Data *Decomp)
{
  int err;
  char loc[13] = "Def_POD_Vars";
  size_t ndims;
  int *dimids, p1, p2;

  if (gsum == 1){ ndims = 2; p1 = 0; p2 = 1;}
  else{ ndims = 3; p1 = 1; p2 = 2;}

  dimids = (int *)malloc(sizeof(int)*ndims);
  dimids[0] = N_g_ID;

  dimids[p1] = clen_ID;
  memset(Decomp[0].name,0,50);
  strcpy(Decomp[0].name,"C_");
  strcat(Decomp[0].name,dname);
  err = nc_def_var(ncid,Decomp[0].name,NC_DOUBLE,(int)ndims-1,dimids,&Decomp[0].id); Handle_Err(err,loc);

  dimids[p1] = rank_ID;
  memset(Decomp[1].name,0,50);
  strcpy(Decomp[1].name,"S_");
  strcat(Decomp[1].name,dname);
  err = nc_def_var(ncid,Decomp[1].name,NC_DOUBLE,(int)ndims-1,dimids,&Decomp[1].id); Handle_Err(err,loc);

  dimids[p1] = rank_ID; dimids[p2] = clen_ID;
  memset(Decomp[2].name,0,50);
  strcpy(Decomp[2].name,"U_");
  strcat(Decomp[2].name,dname);
  err = nc_def_var(ncid,Decomp[2].name,NC_DOUBLE,(int)ndims,dimids,&Decomp[2].id); Handle_Err(err,loc);

  dimids[p1] = N_t_ID; dimids[p2] = rank_ID;
  memset(Decomp[3].name,0,50);
  strcpy(Decomp[3].name,"Vt_");
  strcat(Decomp[3].name,dname);
  err = nc_def_var(ncid,Decomp[3].name,NC_DOUBLE,(int)ndims,dimids,&Decomp[3].id); Handle_Err(err,loc);

  free(dimids);

  return 0;

}

/*================================================================================================================================*/
/**/
/*================================================================================================================================*/
int Def_DMD_Vars(const int ncid, char *dname, const int rank_ID, const int clen_ID, const int N_t_ID, const int N_g_ID,
  const int gsum, Data *Decomp)
{
  int err;
  char loc[13] = "Def_DMD_Vars";
  size_t ndims = 2;
  int *dimids, p1, p2;

  if (gsum == 1){ ndims = 2; p1 = 0; p2 = 1;}
  else{ ndims = 3; p1 = 1; p2 = 2;}

  dimids = (int *)malloc(sizeof(int)*ndims);
  dimids[0] = N_g_ID;

  dimids[p1] = rank_ID;
  memset(Decomp[0].name,0,50);
  strcpy(Decomp[0].name,"L_real_");
  strcat(Decomp[0].name,dname);
  err = nc_def_var(ncid,Decomp[0].name,NC_DOUBLE,(int)ndims-1,dimids,&Decomp[0].id); Handle_Err(err,loc);

  memset(Decomp[1].name,0,50);
  strcpy(Decomp[1].name,"L_imag_");
  strcat(Decomp[1].name,dname);
  err = nc_def_var(ncid,Decomp[1].name,NC_DOUBLE,(int)ndims-1,dimids,&Decomp[1].id); Handle_Err(err,loc);

  dimids[p1] = rank_ID; dimids[p2] = clen_ID;
  memset(Decomp[2].name,0,50);
  strcpy(Decomp[2].name,"W_real_");
  strcat(Decomp[2].name,dname);
  err = nc_def_var(ncid,Decomp[2].name,NC_DOUBLE,(int)ndims,dimids,&Decomp[2].id); Handle_Err(err,loc);

  memset(Decomp[3].name,0,50);
  strcpy(Decomp[3].name,"W_imag_");
  strcat(Decomp[3].name,dname);
  err = nc_def_var(ncid,Decomp[3].name,NC_DOUBLE,(int)ndims,dimids,&Decomp[3].id); Handle_Err(err,loc);

  free(dimids);

  return 0;

}

/*================================================================================================================================*/
/**/
/*================================================================================================================================*/
int Def_DCMP_Vars(const int ncid, const int dcmp_type, const int gsum, Data *Dcmp_data, const size_t N_data, ncdim *clen,
  ncdim *rank, ncdim *dims, Data **Decomp)
{

  int err;

  if (dcmp_type == 0 && gsum == 1){ err = nc_put_att_text(ncid,NC_GLOBAL,"dcmp_type",4,"POD"); }
  else if (dcmp_type == 0 && gsum == 0){ err = nc_put_att_text(ncid,NC_GLOBAL,"dcmp_type",5,"PODg"); }
  else if (dcmp_type == 1 && gsum == 1){ err = nc_put_att_text(ncid,NC_GLOBAL,"dcmp_type",4,"DMD"); }
  else if (dcmp_type == 1 && gsum == 0){ err = nc_put_att_text(ncid,NC_GLOBAL,"dcmp_type",5,"DMDg"); }

  /*------------------------------------------------------------/
  /                             POD                             /
  /------------------------------------------------------------*/
  /*POD uses 4 variables: ID[0] = C (centering vector),
                          ID[1] = S (singular values)
                          ID[2] = U (left singular vectors)
                          ID[3] = Vt (right singular vectors)
  */
  if (dcmp_type == 0){ //decompose with the POD

    *Decomp = (Data*)malloc(sizeof(Data)*N_data*4);
    size_t p = 0;
    for (size_t i=0; i<N_data; i++){
      err = Def_POD_Vars(ncid, Dcmp_data[i].name, rank[Dcmp_data[i].dimids[1]].id, clen[Dcmp_data[i].dimids[0]].id,
        dims[0].id, dims[1].id, gsum, &(*Decomp)[p]);
      p = p + 4;

    }
  }
  /*------------------------------------------------------------/
  /                             DMD                             /
  /------------------------------------------------------------*/
  /*DMD uses 4 variables: ID[0] = L_real (DMD eigenvalues, real component),
                          ID[1] = L_imag (DMD eigenvalues, imaginary component),
                          ID[2] = W_real (DMD modes/ eigenvectors, real component)
                          ID[3] = W_imag (DMD modes/ eigenvectors, imaginary component)
  */
  else if (dcmp_type == 1){ //decompose with the DMD

    *Decomp = (Data*)malloc(sizeof(Data)*N_data*4);
    size_t p = 0;
    for (size_t i=0; i<N_data; i++){
      err = Def_DMD_Vars(ncid, Dcmp_data[i].name, rank[Dcmp_data[i].dimids[1]].id, clen[Dcmp_data[i].dimids[0]].id,
        dims[0].id, dims[1].id, gsum, &(*Decomp)[p]);
      p = p + 4;

    }

  }

  return err;

}

/*================================================================================================================================*/
/**/
/*================================================================================================================================*/
int Generate_DCMP(const double *data, const int ncid_out, const char *dname, const size_t N_t, const size_t N_g, const size_t clen,
  const size_t rank, const double svd_eps, const int *DCMP_IDs, const int dcmp_type)
{
  int err;

  if (dcmp_type == 0){ //decompose with the POD
    err = Generate_POD(data,ncid_out,dname,N_t,N_g,clen,rank,DCMP_IDs);
  }
  else if (dcmp_type == 1){ //decompose with the DMD
    err = Generate_DMD(data,ncid_out,dname,N_t,N_g,clen,rank,svd_eps,DCMP_IDs);
  }

  return err;
}

/*================================================================================================================================*/
/**/
/*================================================================================================================================*/
int Output_BC_DCMP(const int ncid_in, const int ncid_out, const size_t N_t, const size_t N_g, const size_t N_y, const size_t N_x,
  const size_t rank_BC, const double svd_eps, const int gsum, const int *BCg_IDs, const int dcmp_type)
{
  int err;
  char dname[25];
  size_t clen, n_g, gscale, len, p1, p2;
  double *data, *data2;

  if (gsum == 1){ //if decomposing data over all groups, must include length of groups
    n_g = 0; gscale = N_g;
  }
  else{ //decomposing data by group
    n_g = N_g; gscale = 1;
  }

  //calculating vector lengths and data size
  clen = gscale*2*(N_x+N_y);
  len = N_t*N_g*2*(N_x+N_y);
  data = (double *)malloc(sizeof(double)*len);

  printf("... BCg start\n");

  /*------------------------------------------------------------/
  /               Read Cg_L, store in data vector               /
  /------------------------------------------------------------*/
  memset(dname,0,25); strcpy(dname,"Cg_L");
  len = N_t*N_g*N_y; Get_Var_Double(ncid_in,dname,&data2,len);

  len = N_t*N_g; p1=0; p2=0;
  for(size_t i=0; i<len; i++){
    for(size_t j=0; j<N_y; j++){
      data[p1] = data2[p2];
      p1 = p1 + 1;
      p2 = p2 + 1;
    }
    p1 = p1 + N_y + 2*N_x;
  }
  free(data2);

  /*------------------------------------------------------------/
  /               Read Cg_B, store in data vector               /
  /------------------------------------------------------------*/
  memset(dname,0,25); strcpy(dname,"Cg_B");
  len = N_t*N_g*N_x; Get_Var_Double(ncid_in,dname,&data2,len);

  len = N_t*N_g; p1=0; p2=0;
  for(size_t i=0; i<len; i++){
    p1 = p1 + N_y;
    for(size_t j=0; j<N_x; j++){
      data[p1] = data2[p2];
      p1 = p1 + 1;
      p2 = p2 + 1;
    }
    p1 = p1 + N_y + N_x;
  }
  free(data2);

  /*------------------------------------------------------------/
  /               Read Cg_R, store in data vector               /
  /------------------------------------------------------------*/
  memset(dname,0,25); strcpy(dname,"Cg_R");
  len = N_t*N_g*N_y; Get_Var_Double(ncid_in,dname,&data2,len);

  len = N_t*N_g; p1=0; p2=0;
  for(size_t i=0; i<len; i++){
    p1 = p1 + N_y + N_x;
    for(size_t j=0; j<N_y; j++){
      data[p1] = data2[p2];
      p1 = p1 + 1;
      p2 = p2 + 1;
    }
    p1 = p1 + N_x;
  }
  free(data2);

  /*------------------------------------------------------------/
  /               Read Cg_T, store in data vector               /
  /------------------------------------------------------------*/
  memset(dname,0,25); strcpy(dname,"Cg_T");
  len = N_t*N_g*N_x; Get_Var_Double(ncid_in,dname,&data2,len);

  len = N_t*N_g; p1=0; p2=0;
  for(size_t i=0; i<len; i++){
    p1 = p1 + 2*N_y + N_x;
    for(size_t j=0; j<N_x; j++){
      data[p1] = data2[p2];
      p1 = p1 + 1;
      p2 = p2 + 1;
    }
  }
  free(data2);

  /*------------------------------------------------------------/
  /     Now perform Decomposition on the vector containing      /
  /     all boundary factors                                    /
  /------------------------------------------------------------*/
  memset(dname,0,25); strcpy(dname,"BCg");
  err = Generate_DCMP(data,ncid_out,dname,N_t,n_g,clen,rank_BC,svd_eps,BCg_IDs,dcmp_type);
  free(data);

  return err;
}

/*================================================================================================================================*/
/**/
/*================================================================================================================================*/
int Output_QDf_DCMP(const int ncid_in, const int ncid_out, const size_t N_t, const size_t N_g, const size_t N_y, const size_t N_x,
  const size_t rank_avg, const size_t rank_edgV, const size_t rank_edgH, const double svd_eps, const int gsum,
  const int *fg_avg_xx_IDs, const int *fg_edgV_xx_IDs, const int *fg_avg_yy_IDs, const int *fg_edgH_yy_IDs, const int *fg_edgV_xy_IDs,
  const int *fg_edgH_xy_IDs, const int dcmp_type)
{
  int err;
  char dname[25];
  size_t clen, n_g, gscale, len;
  double *data;

  if (gsum == 1){
    n_g = 0; gscale = N_g;
  }
  else{
    n_g = N_g; gscale = 1;
  }

  printf("... fg_avg_xx start\n");
  memset(dname,0,25); strcpy(dname,"fg_avg_xx"); clen = gscale*N_y*N_x;
  len = N_t*N_g*N_y*N_x; Get_Var_Double(ncid_in,dname,&data,len);
  err = Generate_DCMP(data,ncid_out,dname,N_t,n_g,clen,rank_avg,svd_eps,fg_avg_xx_IDs,dcmp_type);
  free(data);

  printf("... fg_edgV_xx start\n");
  memset(dname,0,25); strcpy(dname,"fg_edgV_xx"); clen = gscale*N_y*(N_x+1);
  len = N_t*N_g*N_y*(N_x+1); Get_Var_Double(ncid_in,dname,&data,len);
  err = Generate_DCMP(data,ncid_out,dname,N_t,n_g,clen,rank_edgV,svd_eps,fg_edgV_xx_IDs,dcmp_type);
  free(data);

  printf("... fg_avg_yy start\n");
  memset(dname,0,25); strcpy(dname,"fg_avg_yy"); clen = gscale*N_y*N_x;
  len = N_t*N_g*N_y*N_x; Get_Var_Double(ncid_in,dname,&data,len);
  err = Generate_DCMP(data,ncid_out,dname,N_t,n_g,clen,rank_avg,svd_eps,fg_avg_yy_IDs,dcmp_type);
  free(data);

  printf("... fg_edgH_yy start\n");
  memset(dname,0,25); strcpy(dname,"fg_edgH_yy"); clen = gscale*(N_y+1)*N_x;
  len = N_t*N_g*(N_y+1)*N_x; Get_Var_Double(ncid_in,dname,&data,len);
  err = Generate_DCMP(data,ncid_out,dname,N_t,n_g,clen,rank_edgH,svd_eps,fg_edgH_yy_IDs,dcmp_type);
  free(data);

  printf("... fg_edgV_xy start\n");
  memset(dname,0,25); strcpy(dname,"fg_edgV_xy"); clen = gscale*N_y*(N_x+1);
  len = N_t*N_g*N_y*(N_x+1); Get_Var_Double(ncid_in,dname,&data,len);
  err = Generate_DCMP(data,ncid_out,dname,N_t,n_g,clen,rank_edgV,svd_eps,fg_edgV_xy_IDs,dcmp_type);
  free(data);

  printf("... fg_edgH_xy start\n");
  memset(dname,0,25); strcpy(dname,"fg_edgH_xy"); clen = gscale*(N_y+1)*N_x;
  len = N_t*N_g*(N_y+1)*N_x; Get_Var_Double(ncid_in,dname,&data,len);
  err = Generate_DCMP(data,ncid_out,dname,N_t,n_g,clen,rank_edgH,svd_eps,fg_edgH_xy_IDs,dcmp_type);
  free(data);

  return err;
}

/*================================================================================================================================*/
/**/
/*================================================================================================================================*/
int Output_I_DCMP(const int ncid_in, const int ncid_out, const size_t N_t, const size_t N_g, const size_t N_m, const size_t N_y,
  const size_t N_x, const size_t rank_avg, const size_t rank_edgV, const size_t rank_edgH, const double svd_eps, const int gsum,
  const int *Ig_avg_IDs, const int *Ig_edgV_IDs, const int *Ig_edgH_IDs, const int dcmp_type)
{
  int err;
  char dname[25];
  size_t clen, n_g, gscale, len;
  double *data;

  if (gsum == 1){
    n_g = 0; gscale = N_g;
  }
  else{
    n_g = N_g; gscale = 1;
  }

  printf("... I_avg start\n");
  memset(dname,0,25); strcpy(dname,"I_avg"); clen = gscale*N_m*N_y*N_x;
  len = N_t*N_g*N_m*N_y*N_x; Get_Var_Double(ncid_in,dname,&data,len);
  err = Generate_DCMP(data,ncid_out,dname,N_t,n_g,clen,rank_avg,svd_eps,Ig_avg_IDs,dcmp_type);
  free(data);

  printf("... I_edgV start\n");
  memset(dname,0,25); strcpy(dname,"I_edgV"); clen = gscale*N_m*N_y*(N_x+1);
  len = N_t*N_g*N_m*N_y*(N_x+1); Get_Var_Double(ncid_in,dname,&data,len);
  err = Generate_DCMP(data,ncid_out,dname,N_t,n_g,clen,rank_edgV,svd_eps,Ig_edgV_IDs,dcmp_type);
  free(data);

  printf("... I_edgH start\n");
  memset(dname,0,25); strcpy(dname,"I_edgH"); clen = gscale*N_m*(N_y+1)*N_x;
  len = N_t*N_g*N_m*(N_y+1)*N_x; Get_Var_Double(ncid_in,dname,&data,len);
  err = Generate_DCMP(data,ncid_out,dname,N_t,n_g,clen,rank_edgH,svd_eps,Ig_edgH_IDs,dcmp_type);
  free(data);

  return err;
}

/*================================================================================================================================*/
/**/
/*================================================================================================================================*/
int Output_meanI_DCMP(const int ncid_in, const int ncid_out, const size_t N_t, const size_t N_g, const size_t N_m, const size_t N_y,
  const size_t N_x, const size_t rank_avg, const size_t rank_edgV, const size_t rank_edgH, const double svd_eps, const int gsum,
  const int *Ig_avg_IDs, const int *Ig_edgV_IDs, const int *Ig_edgH_IDs, const int dcmp_type)
{
  int err;
  char dname[25];
  size_t clen, n_g, gscale, len, len2, p1, p2;
  double *data, *data2;
  double c;

  c = 299.792458;
  if (gsum == 1){
    n_g = 0; gscale = N_g;
  }
  else{
    n_g = N_g; gscale = 1;
  }

  /*------------------------------------------------------------/
  /       Average I_Avg with scalar intensities, decompose      /
  /------------------------------------------------------------*/
  printf("... mean I_avg start\n");
  memset(dname,0,25); strcpy(dname,"Eg_avg_HO");
  len = N_t*N_g*N_y*N_x; Get_Var_Double(ncid_in,dname,&data2,len);

  memset(dname,0,25); strcpy(dname,"I_avg");
  len = N_t*N_g*N_m*N_y*N_x; Get_Var_Double(ncid_in,dname,&data,len);

  len = N_t*N_g;
  len2 = N_y*N_x;
  p1 = 0;
  for(size_t i=0; i<len; i++){
    for(size_t m=0; m<N_m; m++){
      p2 = i*len2;
      for(size_t j=0; j<len2; j++){

        data[p1] = data[p1]/(data2[p2]*c);
        p1 = p1 + 1;
        p2 = p2 + 1;

      }
    }
  }
  free(data2);

  memset(dname,0,25); strcpy(dname,"Mean_I_avg"); clen = gscale*N_m*N_y*N_x;
  err = Generate_DCMP(data,ncid_out,dname,N_t,n_g,clen,rank_avg,svd_eps,Ig_avg_IDs,dcmp_type);
  free(data);


  /*------------------------------------------------------------/
  /      Average I_EdgV with scalar intensities, decompose      /
  /------------------------------------------------------------*/
  printf("... mean I_edgV start\n");
  memset(dname,0,25); strcpy(dname,"Eg_edgV_HO");
  len = N_t*N_g*N_y*(N_x+1); Get_Var_Double(ncid_in,dname,&data2,len);

  memset(dname,0,25); strcpy(dname,"I_edgV");
  len = N_t*N_g*N_m*N_y*(N_x+1); Get_Var_Double(ncid_in,dname,&data,len);

  len = N_t*N_g;
  len2 = N_y*(N_x+1);
  p1 = 0;
  for(size_t i=0; i<len; i++){
    for(size_t m=0; m<N_m; m++){
      p2 = i*len2;
      for(size_t j=0; j<len2; j++){

        data[p1] = data[p1]/(data2[p2]*c);
        p1 = p1 + 1;
        p2 = p2 + 1;

      }
    }
  }
  free(data2);

  memset(dname,0,25); strcpy(dname,"Mean_I_edgV"); clen = gscale*N_m*N_y*(N_x+1);
  err = Generate_DCMP(data,ncid_out,dname,N_t,n_g,clen,rank_edgV,svd_eps,Ig_edgV_IDs,dcmp_type);
  free(data);

  /*------------------------------------------------------------/
  /      Average I_EdgH with scalar intensities, decompose      /
  /------------------------------------------------------------*/
  printf("... mean I_edgH start\n");
  memset(dname,0,25); strcpy(dname,"Eg_edgH_HO");
  len = N_t*N_g*(N_y+1)*N_x; Get_Var_Double(ncid_in,dname,&data2,len);

  memset(dname,0,25); strcpy(dname,"I_edgH");
  len = N_t*N_g*N_m*(N_y+1)*N_x; Get_Var_Double(ncid_in,dname,&data,len);

  len = N_t*N_g;
  len2 = (N_y+1)*N_x;
  p1 = 0;
  for(size_t i=0; i<len; i++){
    for(size_t m=0; m<N_m; m++){
      p2 = i*len2;
      for(size_t j=0; j<len2; j++){

        data[p1] = data[p1]/(data2[p2]*c);
        p1 = p1 + 1;
        p2 = p2 + 1;

      }
    }
  }
  free(data2);

  memset(dname,0,25); strcpy(dname,"Mean_I_edgH"); clen = gscale*N_m*(N_y+1)*N_x;
  err = Generate_DCMP(data,ncid_out,dname,N_t,n_g,clen,rank_edgH,svd_eps,Ig_edgH_IDs,dcmp_type);
  free(data);

  return err;
}

/*================================================================================================================================*/
/**/
/*================================================================================================================================*/
int Decompose_Data(const int ncid_in, const int ncid_out, const int dcmp_type, const int dcmp_data, const int gsum,
  const size_t N_t, const size_t N_g, const size_t N_m, const size_t N_x, const size_t N_y, const size_t rank_BC,
  const size_t rank_avg, const size_t rank_edgV, const size_t rank_edgH, const double svd_eps, const int *BCg_IDs,
  const int *fg_avg_xx_IDs, const int *fg_edgV_xx_IDs, const int *fg_avg_yy_IDs, const int *fg_edgH_yy_IDs,
  const int *fg_edgV_xy_IDs, const int *fg_edgH_xy_IDs, const int *Ig_avg_IDs,const int *Ig_edgV_IDs, const int *Ig_edgH_IDs)
{
  int err;

  if (dcmp_data == 0){ //decompose QD factors (and boundary factors)

    err = Output_BC_DCMP(ncid_in,ncid_out,N_t,N_g,N_y,N_x,rank_BC,svd_eps,gsum,BCg_IDs,dcmp_type);

    err = Output_QDf_DCMP(ncid_in,ncid_out,N_t,N_g,N_y,N_x,rank_avg,rank_edgV,rank_edgH,svd_eps,gsum,fg_avg_xx_IDs,
       fg_edgV_xx_IDs,fg_avg_yy_IDs,fg_edgH_yy_IDs,fg_edgV_xy_IDs,fg_edgH_xy_IDs,dcmp_type);

  }
  else if (dcmp_data == 1){ //decompose Intensities

    err = Output_I_DCMP(ncid_in,ncid_out,N_t,N_g,N_m,N_y,N_x,rank_avg,rank_edgV,rank_edgH,svd_eps,gsum,Ig_avg_IDs,
      Ig_edgV_IDs,Ig_edgH_IDs,dcmp_type);

  }
  else if (dcmp_data == 2){ //decompose mean Intensities

    err = Output_meanI_DCMP(ncid_in,ncid_out,N_t,N_g,N_m,N_y,N_x,rank_avg,rank_edgV,rank_edgH,svd_eps,gsum,Ig_avg_IDs,
      Ig_edgV_IDs,Ig_edgH_IDs,dcmp_type);

  }

  return err;
}
