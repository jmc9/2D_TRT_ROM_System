/*================================================================================================================================*/
/*
  OUTPUTS.c
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

void Get_Var_Double(const int ncid, const char *name, double **var, const size_t size);

/* ----- FROM CPOD_ROUTINES.c ----- */
int Generate_POD(const double *data, const int ncid_out, const char *dname, const size_t N_t, const size_t N_g, const size_t clen,
  const size_t rank, const int *DCMP_IDs);

int Generate_DMD(const double *data, const int ncid_out, const char *dname, const size_t N_t, const size_t N_g, const size_t clen,
  const size_t rank, const int *DCMP_IDs);

/* ----- LOCAL DEFINITIONS ----- */
#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })

/*================================================================================================================================*/
/**/
/*================================================================================================================================*/
int Def_Dims(const int ncid_out, const size_t N_t, const size_t N_g, const size_t N_m, const size_t N_y, const size_t N_x,
  const size_t rank_BC, const size_t rank_avg, const size_t rank_edgV,const size_t rank_edgH, const size_t BClen,
  const size_t clen_avg, const size_t clen_edgV, const size_t clen_edgH, const double tlen, const double Delt, const double xlen,
  const double ylen, double *Delx, double *Dely, int *BC_Type, double *bcT, const double Tini, int *N_t_ID, int *N_g_ID,
  int *N_m_ID, int *N_y_ID, int *N_x_ID, int *rank_BC_ID, int *rank_avg_ID, int *rank_edgV_ID, int *rank_edgH_ID, int *BClen_ID,
  int *clen_avg_ID, int *clen_edgV_ID, int *clen_edgH_ID)
{
  int err, id, Bnds_ID, dims[1];
  char loc[9] = "OUT_DIMS";

  err = nc_def_dim(ncid_out,"N_x",N_x,N_x_ID); Handle_Err(err,loc);
  err = nc_def_dim(ncid_out,"N_y",N_y,N_y_ID); Handle_Err(err,loc);
  err = nc_def_dim(ncid_out,"N_m",N_m,N_m_ID); Handle_Err(err,loc);
  err = nc_def_dim(ncid_out,"N_g",N_g,N_g_ID); Handle_Err(err,loc);
  err = nc_def_dim(ncid_out,"N_t",N_t,N_t_ID); Handle_Err(err,loc);
  err = nc_def_dim(ncid_out,"Boundaries",4,&Bnds_ID); Handle_Err(err,loc);

  if (rank_BC != 0){
    err = nc_def_dim(ncid_out,"rank_BC",rank_BC,rank_BC_ID); Handle_Err(err,loc);
  }
  err = nc_def_dim(ncid_out,"rank_avg",rank_avg,rank_avg_ID); Handle_Err(err,loc);
  err = nc_def_dim(ncid_out,"rank_edgV",rank_edgV,rank_edgV_ID); Handle_Err(err,loc);
  err = nc_def_dim(ncid_out,"rank_edgH",rank_edgH,rank_edgH_ID); Handle_Err(err,loc);

  if (BClen != 0){
    err = nc_def_dim(ncid_out,"BClen",BClen,BClen_ID); Handle_Err(err,loc);
  }
  err = nc_def_dim(ncid_out,"clen_avg",clen_avg,clen_avg_ID); Handle_Err(err,loc);
  err = nc_def_dim(ncid_out,"clen_edgV",clen_edgV,clen_edgV_ID); Handle_Err(err,loc);
  err = nc_def_dim(ncid_out,"clen_edgH",clen_edgH,clen_edgH_ID); Handle_Err(err,loc);

  err = nc_def_var(ncid_out,"tlen",NC_DOUBLE,0,0,&id); Handle_Err(err,loc);
  err = nc_put_var_double(ncid_out,id,&tlen); Handle_Err(err,loc);

  err = nc_def_var(ncid_out,"Delt",NC_DOUBLE,0,0,&id); Handle_Err(err,loc);
  err = nc_put_var_double(ncid_out,id,&Delt); Handle_Err(err,loc);

  err = nc_def_var(ncid_out,"xlen",NC_DOUBLE,0,0,&id); Handle_Err(err,loc);
  err = nc_put_var_double(ncid_out,id,&xlen); Handle_Err(err,loc);

  err = nc_def_var(ncid_out,"ylen",NC_DOUBLE,0,0,&id); Handle_Err(err,loc);
  err = nc_put_var_double(ncid_out,id,&ylen); Handle_Err(err,loc);

  dims[0] = (*N_x_ID);
  err = nc_def_var(ncid_out,"Delx",NC_DOUBLE,1,dims,&id); Handle_Err(err,loc);
  err = nc_put_var_double(ncid_out,id,Delx); Handle_Err(err,loc);

  dims[0] = (*N_y_ID);
  err = nc_def_var(ncid_out,"Dely",NC_DOUBLE,1,dims,&id); Handle_Err(err,loc);
  err = nc_put_var_double(ncid_out,id,Dely); Handle_Err(err,loc);

  dims[0] = Bnds_ID;
  err = nc_def_var(ncid_out,"BC_Type",NC_DOUBLE,1,dims,&id); Handle_Err(err,loc);
  err = nc_put_var_int(ncid_out,id,BC_Type); Handle_Err(err,loc);

  dims[0] = Bnds_ID;
  err = nc_def_var(ncid_out,"bcT",NC_DOUBLE,1,dims,&id); Handle_Err(err,loc);
  err = nc_put_var_double(ncid_out,id,bcT); Handle_Err(err,loc);

  err = nc_def_var(ncid_out,"Tini",NC_DOUBLE,0,0,&id); Handle_Err(err,loc);
  err = nc_put_var_double(ncid_out,id,&Tini); Handle_Err(err,loc);

  return err;
}

/*================================================================================================================================*/
/**/
/*================================================================================================================================*/
int Def_POD_Vars(const int ncid, const char *vname, const int rank_ID, const int clen_ID, const int N_t_ID, const int N_g_ID,
  const int gsum, int *DCMP_IDs)
{
  int err;
  char loc[13] = "Def_POD_Vars";
  size_t ndims = 2;
  int *dimids, p1, p2;
  char name[128];

  if (gsum == 1){ ndims = 2; p1 = 0; p2 = 1;}
  else{ ndims = 3; p1 = 1; p2 = 2;}

  dimids = (int *)malloc(sizeof(int)*ndims);
  dimids[0] = N_g_ID;

  dimids[p1] = clen_ID;
  memset(name,0,128); strcpy(name,"C_"); strcat(name,vname);
  err = nc_def_var(ncid,name,NC_DOUBLE,(int)ndims-1,dimids,&DCMP_IDs[0]); Handle_Err(err,loc);

  dimids[p1] = rank_ID;
  memset(name,0,128); strcpy(name,"S_"); strcat(name,vname);
  err = nc_def_var(ncid,name,NC_DOUBLE,(int)ndims-1,dimids,&DCMP_IDs[1]); Handle_Err(err,loc);

  dimids[p1] = rank_ID; dimids[p2] = clen_ID;
  memset(name,0,128); strcpy(name,"U_"); strcat(name,vname);
  err = nc_def_var(ncid,name,NC_DOUBLE,(int)ndims,dimids,&DCMP_IDs[2]); Handle_Err(err,loc);

  dimids[p1] = N_t_ID; dimids[p2] = rank_ID;
  memset(name,0,128); strcpy(name,"Vt_"); strcat(name,vname);
  err = nc_def_var(ncid,name,NC_DOUBLE,(int)ndims,dimids,&DCMP_IDs[3]); Handle_Err(err,loc);

  free(dimids);

  return 0;

}

/*================================================================================================================================*/
/**/
/*================================================================================================================================*/
int Def_DMD_Vars(const int ncid, const char *vname, const int rank_ID, const int clen_ID, const int N_t_ID, const int N_g_ID,
  const int gsum, int *DCMP_IDs)
{
  int err;
  char loc[13] = "Def_DMD_Vars";
  size_t ndims = 2;
  int *dimids, p1, p2;
  char name[128];

  if (gsum == 1){ ndims = 2; p1 = 0; p2 = 1;}
  else{ ndims = 3; p1 = 1; p2 = 2;}

  dimids = (int *)malloc(sizeof(int)*ndims);
  dimids[0] = N_g_ID;

  dimids[p1] = rank_ID;
  memset(name,0,128); strcpy(name,"L_real_"); strcat(name,vname);
  err = nc_def_var(ncid,name,NC_DOUBLE,(int)ndims-1,dimids,&DCMP_IDs[0]); Handle_Err(err,loc);

  memset(name,0,128); strcpy(name,"L_imag_"); strcat(name,vname);
  err = nc_def_var(ncid,name,NC_DOUBLE,(int)ndims-1,dimids,&DCMP_IDs[1]); Handle_Err(err,loc);

  dimids[p1] = rank_ID; dimids[p2] = clen_ID;
  memset(name,0,128); strcpy(name,"W_real_"); strcat(name,vname);
  err = nc_def_var(ncid,name,NC_DOUBLE,(int)ndims,dimids,&DCMP_IDs[2]); Handle_Err(err,loc);

  memset(name,0,128); strcpy(name,"W_imag_"); strcat(name,vname);
  err = nc_def_var(ncid,name,NC_DOUBLE,(int)ndims,dimids,&DCMP_IDs[3]); Handle_Err(err,loc);

  free(dimids);

  return 0;

}

/*================================================================================================================================*/
/**/
/*================================================================================================================================*/
int Def_DCMP_Vars(const int ncid, const int dcmp_type, const int dcmp_data, const int gsum, const int N_g_ID, const int rank_BC_ID,
  const int rank_avg_ID,const int rank_edgV_ID, const int rank_edgH_ID, const int BClen_ID, const int clen_avg_ID, const int clen_edgV_ID,
  const int clen_edgH_ID, const int N_t_ID, int **BCg_IDs, int **fg_avg_xx_IDs, int **fg_edgV_xx_IDs, int **fg_avg_yy_IDs, int **fg_edgH_yy_IDs,
  int **fg_edgV_xy_IDs, int **fg_edgH_xy_IDs, int **Ig_avg_IDs, int **Ig_edgV_IDs, int **Ig_edgH_IDs)
{

  int err;

  /*------------------------------------------------------------/
  /                             POD                             /
  /------------------------------------------------------------*/
  /*POD uses 4 variables: ID[0] = C (centering vector),
                          ID[1] = S (singular values)
                          ID[2] = U (left singular vectors)
                          ID[3] = Vt (right singular vectors)
  */
  if (dcmp_type == 0){ //decompose with the POD
    if (dcmp_data == 0){ //decompose QD factors (and boundary factors)

      *BCg_IDs = (int *)malloc(sizeof(int)*4);
      *fg_avg_xx_IDs = (int *)malloc(sizeof(int)*4);
      *fg_edgV_xx_IDs = (int *)malloc(sizeof(int)*4);
      *fg_avg_yy_IDs = (int *)malloc(sizeof(int)*4);
      *fg_edgH_yy_IDs = (int *)malloc(sizeof(int)*4);
      *fg_edgV_xy_IDs = (int *)malloc(sizeof(int)*4);
      *fg_edgH_xy_IDs = (int *)malloc(sizeof(int)*4);

      err = Def_POD_Vars(ncid,"BCg",rank_BC_ID,BClen_ID,N_t_ID,N_g_ID,gsum,*BCg_IDs);
      err = Def_POD_Vars(ncid,"fg_avg_xx",rank_avg_ID,clen_avg_ID,N_t_ID,N_g_ID,gsum,*fg_avg_xx_IDs);
      err = Def_POD_Vars(ncid,"fg_edgV_xx",rank_edgV_ID,clen_edgV_ID,N_t_ID,N_g_ID,gsum,*fg_edgV_xx_IDs);
      err = Def_POD_Vars(ncid,"fg_avg_yy",rank_avg_ID,clen_avg_ID,N_t_ID,N_g_ID,gsum,*fg_avg_yy_IDs);
      err = Def_POD_Vars(ncid,"fg_edgH_yy",rank_edgH_ID,clen_edgH_ID,N_t_ID,N_g_ID,gsum,*fg_edgH_yy_IDs);
      err = Def_POD_Vars(ncid,"fg_edgV_xy",rank_edgV_ID,clen_edgV_ID,N_t_ID,N_g_ID,gsum,*fg_edgV_xy_IDs);
      err = Def_POD_Vars(ncid,"fg_edgH_xy",rank_edgH_ID,clen_edgH_ID,N_t_ID,N_g_ID,gsum,*fg_edgH_xy_IDs);

    }
    else if (dcmp_data == 1){ //decompose Intensities

      *Ig_avg_IDs = (int *)malloc(sizeof(int)*4);
      *Ig_edgV_IDs = (int *)malloc(sizeof(int)*4);
      *Ig_edgH_IDs = (int *)malloc(sizeof(int)*4);

      err = Def_POD_Vars(ncid,"Ig_avg",rank_avg_ID,clen_avg_ID,N_t_ID,N_g_ID,gsum,*Ig_avg_IDs);
      err = Def_POD_Vars(ncid,"Ig_edgV",rank_edgV_ID,clen_edgV_ID,N_t_ID,N_g_ID,gsum,*Ig_edgV_IDs);
      err = Def_POD_Vars(ncid,"Ig_edgH",rank_edgH_ID,clen_edgH_ID,N_t_ID,N_g_ID,gsum,*Ig_edgH_IDs);

    }
    else if (dcmp_data == 2){ //decompose mean Intensities

      *Ig_avg_IDs = (int *)malloc(sizeof(int)*4);
      *Ig_edgV_IDs = (int *)malloc(sizeof(int)*4);
      *Ig_edgH_IDs = (int *)malloc(sizeof(int)*4);

      err = Def_POD_Vars(ncid,"Mean_Ig_avg",rank_avg_ID,clen_avg_ID,N_t_ID,N_g_ID,gsum,*Ig_avg_IDs);
      err = Def_POD_Vars(ncid,"Mean_Ig_edgV",rank_edgV_ID,clen_edgV_ID,N_t_ID,N_g_ID,gsum,*Ig_edgV_IDs);
      err = Def_POD_Vars(ncid,"Mean_Ig_edgH",rank_edgH_ID,clen_edgH_ID,N_t_ID,N_g_ID,gsum,*Ig_edgH_IDs);

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
    if (dcmp_data == 0){ //decompose QD factors (and boundary factors)

      *BCg_IDs = (int *)malloc(sizeof(int)*4);
      *fg_avg_xx_IDs = (int *)malloc(sizeof(int)*4);
      *fg_edgV_xx_IDs = (int *)malloc(sizeof(int)*4);
      *fg_avg_yy_IDs = (int *)malloc(sizeof(int)*4);
      *fg_edgH_yy_IDs = (int *)malloc(sizeof(int)*4);
      *fg_edgV_xy_IDs = (int *)malloc(sizeof(int)*4);
      *fg_edgH_xy_IDs = (int *)malloc(sizeof(int)*4);

      err = Def_DMD_Vars(ncid,"BCg",rank_BC_ID,BClen_ID,N_t_ID,N_g_ID,gsum,*BCg_IDs);
      err = Def_DMD_Vars(ncid,"fg_avg_xx",rank_avg_ID,clen_avg_ID,N_t_ID,N_g_ID,gsum,*fg_avg_xx_IDs);
      err = Def_DMD_Vars(ncid,"fg_edgV_xx",rank_edgV_ID,clen_edgV_ID,N_t_ID,N_g_ID,gsum,*fg_edgV_xx_IDs);
      err = Def_DMD_Vars(ncid,"fg_avg_yy",rank_avg_ID,clen_avg_ID,N_t_ID,N_g_ID,gsum,*fg_avg_yy_IDs);
      err = Def_DMD_Vars(ncid,"fg_edgH_yy",rank_edgH_ID,clen_edgH_ID,N_t_ID,N_g_ID,gsum,*fg_edgH_yy_IDs);
      err = Def_DMD_Vars(ncid,"fg_edgV_xy",rank_edgV_ID,clen_edgV_ID,N_t_ID,N_g_ID,gsum,*fg_edgV_xy_IDs);
      err = Def_DMD_Vars(ncid,"fg_edgH_xy",rank_edgH_ID,clen_edgH_ID,N_t_ID,N_g_ID,gsum,*fg_edgH_xy_IDs);

    }
    else if (dcmp_data == 1){ //decompose Intensities

      *Ig_avg_IDs = (int *)malloc(sizeof(int)*4);
      *Ig_edgV_IDs = (int *)malloc(sizeof(int)*4);
      *Ig_edgH_IDs = (int *)malloc(sizeof(int)*4);

      err = Def_DMD_Vars(ncid,"Ig_avg",rank_avg_ID,clen_avg_ID,N_t_ID,N_g_ID,gsum,*Ig_avg_IDs);
      err = Def_DMD_Vars(ncid,"Ig_edgV",rank_edgV_ID,clen_edgV_ID,N_t_ID,N_g_ID,gsum,*Ig_edgV_IDs);
      err = Def_DMD_Vars(ncid,"Ig_edgH",rank_edgH_ID,clen_edgH_ID,N_t_ID,N_g_ID,gsum,*Ig_edgH_IDs);

    }
    else if (dcmp_data == 2){ //decompose mean Intensities

      *Ig_avg_IDs = (int *)malloc(sizeof(int)*4);
      *Ig_edgV_IDs = (int *)malloc(sizeof(int)*4);
      *Ig_edgH_IDs = (int *)malloc(sizeof(int)*4);

      err = Def_DMD_Vars(ncid,"Mean_Ig_avg",rank_avg_ID,clen_avg_ID,N_t_ID,N_g_ID,gsum,*Ig_avg_IDs);
      err = Def_DMD_Vars(ncid,"Mean_Ig_edgV",rank_edgV_ID,clen_edgV_ID,N_t_ID,N_g_ID,gsum,*Ig_edgV_IDs);
      err = Def_DMD_Vars(ncid,"Mean_Ig_edgH",rank_edgH_ID,clen_edgH_ID,N_t_ID,N_g_ID,gsum,*Ig_edgH_IDs);

    }
  }

  return err;

}

/*================================================================================================================================*/
/**/
/*================================================================================================================================*/
int Generate_DCMP(const double *data, const int ncid_out, const char *dname, const size_t N_t, const size_t N_g, const size_t clen,
  const size_t rank, const int *DCMP_IDs, const int dcmp_type)
{
  int err;

  if (dcmp_type == 0){ //decompose with the POD
    err = Generate_POD(data,ncid_out,dname,N_t,N_g,clen,rank,DCMP_IDs);
  }
  else if (dcmp_type == 1){ //decompose with the DMD
    err = Generate_DMD(data,ncid_out,dname,N_t,N_g,clen,rank,DCMP_IDs);
  }

  return err;
}

/*================================================================================================================================*/
/**/
/*================================================================================================================================*/
int Output_BC_DCMP(const int ncid_in, const int ncid_out, const size_t N_t, const size_t N_g, const size_t N_y, const size_t N_x,
  const size_t rank_BC, const int gsum, const int *BCg_IDs, const int dcmp_type)
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
  memset(dname,0,25); strcpy(dname,"Cg");
  err = Generate_DCMP(data,ncid_out,dname,N_t,n_g,clen,rank_BC,BCg_IDs,dcmp_type);
  free(data);

  return err;
}

/*================================================================================================================================*/
/**/
/*================================================================================================================================*/
int Output_QDf_DCMP(const int ncid_in, const int ncid_out, const size_t N_t, const size_t N_g, const size_t N_y, const size_t N_x,
  const size_t rank_avg, const size_t rank_edgV, const size_t rank_edgH, const int gsum, const int *fg_avg_xx_IDs, const int *fg_edgV_xx_IDs, const int *fg_avg_yy_IDs, const int *fg_edgH_yy_IDs,
  const int *fg_edgV_xy_IDs, const int *fg_edgH_xy_IDs, const int dcmp_type)
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
  err = Generate_DCMP(data,ncid_out,dname,N_t,n_g,clen,rank_avg,fg_avg_xx_IDs,dcmp_type);
  free(data);

  printf("... fg_edgV_xx start\n");
  memset(dname,0,25); strcpy(dname,"fg_edgV_xx"); clen = gscale*N_y*(N_x+1);
  len = N_t*N_g*N_y*(N_x+1); Get_Var_Double(ncid_in,dname,&data,len);
  err = Generate_DCMP(data,ncid_out,dname,N_t,n_g,clen,rank_edgV,fg_edgV_xx_IDs,dcmp_type);
  free(data);

  printf("... fg_avg_yy start\n");
  memset(dname,0,25); strcpy(dname,"fg_avg_yy"); clen = gscale*N_y*N_x;
  len = N_t*N_g*N_y*N_x; Get_Var_Double(ncid_in,dname,&data,len);
  err = Generate_DCMP(data,ncid_out,dname,N_t,n_g,clen,rank_avg,fg_avg_yy_IDs,dcmp_type);
  free(data);

  printf("... fg_edgH_yy start\n");
  memset(dname,0,25); strcpy(dname,"fg_edgH_yy"); clen = gscale*(N_y+1)*N_x;
  len = N_t*N_g*(N_y+1)*N_x; Get_Var_Double(ncid_in,dname,&data,len);
  err = Generate_DCMP(data,ncid_out,dname,N_t,n_g,clen,rank_edgH,fg_edgH_yy_IDs,dcmp_type);
  free(data);

  printf("... fg_edgV_xy start\n");
  memset(dname,0,25); strcpy(dname,"fg_edgV_xy"); clen = gscale*N_y*(N_x+1);
  len = N_t*N_g*N_y*(N_x+1); Get_Var_Double(ncid_in,dname,&data,len);
  err = Generate_DCMP(data,ncid_out,dname,N_t,n_g,clen,rank_edgV,fg_edgV_xy_IDs,dcmp_type);
  free(data);

  printf("... fg_edgH_xy start\n");
  memset(dname,0,25); strcpy(dname,"fg_edgH_xy"); clen = gscale*(N_y+1)*N_x;
  len = N_t*N_g*(N_y+1)*N_x; Get_Var_Double(ncid_in,dname,&data,len);
  err = Generate_DCMP(data,ncid_out,dname,N_t,n_g,clen,rank_edgH,fg_edgH_xy_IDs,dcmp_type);
  free(data);

  return err;
}

/*================================================================================================================================*/
/**/
/*================================================================================================================================*/
int Output_I_DCMP(const int ncid_in, const int ncid_out, const size_t N_t, const size_t N_g, const size_t N_m, const size_t N_y,
  const size_t N_x, const size_t rank_avg, const size_t rank_edgV, const size_t rank_edgH, const int gsum, const int *Ig_avg_IDs, const int *Ig_edgV_IDs, const int *Ig_edgH_IDs, const int dcmp_type)
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
  err = Generate_DCMP(data,ncid_out,dname,N_t,n_g,clen,rank_avg,Ig_avg_IDs,dcmp_type);
  free(data);

  printf("... I_edgV start\n");
  memset(dname,0,25); strcpy(dname,"I_edgV"); clen = gscale*N_m*N_y*(N_x+1);
  len = N_t*N_g*N_m*N_y*(N_x+1); Get_Var_Double(ncid_in,dname,&data,len);
  err = Generate_DCMP(data,ncid_out,dname,N_t,n_g,clen,rank_edgV,Ig_edgV_IDs,dcmp_type);
  free(data);

  printf("... I_edgH start\n");
  memset(dname,0,25); strcpy(dname,"I_edgH"); clen = gscale*N_m*(N_y+1)*N_x;
  len = N_t*N_g*N_m*(N_y+1)*N_x; Get_Var_Double(ncid_in,dname,&data,len);
  err = Generate_DCMP(data,ncid_out,dname,N_t,n_g,clen,rank_edgH,Ig_edgH_IDs,dcmp_type);
  free(data);

  return err;
}

/*================================================================================================================================*/
/**/
/*================================================================================================================================*/
int Output_meanI_DCMP(const int ncid_in, const int ncid_out, const size_t N_t, const size_t N_g, const size_t N_m, const size_t N_y,
  const size_t N_x, const size_t rank_avg, const size_t rank_edgV, const size_t rank_edgH, const int gsum, const int *Ig_avg_IDs, const int *Ig_edgV_IDs, const int *Ig_edgH_IDs, const int dcmp_type)
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
  err = Generate_DCMP(data,ncid_out,dname,N_t,n_g,clen,rank_avg,Ig_avg_IDs,dcmp_type);
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
  err = Generate_DCMP(data,ncid_out,dname,N_t,n_g,clen,rank_edgV,Ig_edgV_IDs,dcmp_type);
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
  err = Generate_DCMP(data,ncid_out,dname,N_t,n_g,clen,rank_edgH,Ig_edgH_IDs,dcmp_type);
  free(data);

  return err;
}

/*================================================================================================================================*/
/**/
/*================================================================================================================================*/
int Decompose_Data(const int ncid_in, const int ncid_out, const int dcmp_type, const int dcmp_data, const int gsum,
  const size_t N_t, const size_t N_g, const size_t N_m, const size_t N_x, const size_t N_y, const size_t rank_BC,
  const size_t rank_avg, const size_t rank_edgV, const size_t rank_edgH, const int *BCg_IDs, const int *fg_avg_xx_IDs,
  const int *fg_edgV_xx_IDs, const int *fg_avg_yy_IDs, const int *fg_edgH_yy_IDs,const int *fg_edgV_xy_IDs,
  const int *fg_edgH_xy_IDs, const int *Ig_avg_IDs,const int *Ig_edgV_IDs, const int *Ig_edgH_IDs)
{
  int err;

  if (dcmp_data == 0){ //decompose QD factors (and boundary factors)

    err = Output_BC_DCMP(ncid_in,ncid_out,N_t,N_g,N_y,N_x,rank_BC,gsum,BCg_IDs,dcmp_type);

    err = Output_QDf_DCMP(ncid_in,ncid_out,N_t,N_g,N_y,N_x,rank_avg,rank_edgV,rank_edgH,gsum,fg_avg_xx_IDs,
       fg_edgV_xx_IDs,fg_avg_yy_IDs,fg_edgH_yy_IDs,fg_edgV_xy_IDs,fg_edgH_xy_IDs,dcmp_type);

  }
  else if (dcmp_data == 1){ //decompose Intensities

    err = Output_I_DCMP(ncid_in,ncid_out,N_t,N_g,N_m,N_y,N_x,rank_avg,rank_edgV,rank_edgH,gsum,Ig_avg_IDs,
      Ig_edgV_IDs,Ig_edgH_IDs,dcmp_type);

  }
  else if (dcmp_data == 2){ //decompose mean Intensities

    err = Output_meanI_DCMP(ncid_in,ncid_out,N_t,N_g,N_m,N_y,N_x,rank_avg,rank_edgV,rank_edgH,gsum,Ig_avg_IDs,
      Ig_edgV_IDs,Ig_edgH_IDs,dcmp_type);

  }

  return err;
}
