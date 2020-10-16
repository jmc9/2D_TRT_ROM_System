#include <stdlib.h>
#include <stdio.h>
#include <netcdf.h>
#include <string.h>

void HANDLE_ERR(int Status, char Location[]);

int GENERATE_POD(int ncid_in, int ncid_out, char *dname, size_t N_t, size_t N_g, size_t clen, size_t rank, int gsum,
  int Sid, int Uid, int Vtid);

#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })

//================================================================================================================================//
//
//================================================================================================================================//
int DEF_DIMS(int ncid_out, size_t N_t, size_t N_g, size_t N_y, size_t N_x, size_t rank_avg, size_t rank_edgV, size_t rank_edgH,
   size_t clen_avg, size_t clen_edgV, size_t clen_edgH, int *N_t_ID, int *N_g_ID, int *N_y_ID, int *N_x_ID, int *rank_avg_ID,
   int *rank_edgV_ID, int *rank_edgH_ID, int *clen_avg_ID, int *clen_edgV_ID, int *clen_edgH_ID)
{
  int err;
  char loc[9] = "OUT_DIMS";

  err = nc_def_dim(ncid_out,"N_x",N_x,N_x_ID); HANDLE_ERR(err,loc);
  err = nc_def_dim(ncid_out,"N_y",N_y,N_y_ID); HANDLE_ERR(err,loc);
  err = nc_def_dim(ncid_out,"N_g",N_g,N_g_ID); HANDLE_ERR(err,loc);
  err = nc_def_dim(ncid_out,"N_t",N_t,N_t_ID); HANDLE_ERR(err,loc);

  err = nc_def_dim(ncid_out,"rank_avg",rank_avg,rank_avg_ID); HANDLE_ERR(err,loc);
  err = nc_def_dim(ncid_out,"rank_edgV",rank_edgV,rank_edgV_ID); HANDLE_ERR(err,loc);
  err = nc_def_dim(ncid_out,"rank_edgH",rank_edgH,rank_edgH_ID); HANDLE_ERR(err,loc);

  err = nc_def_dim(ncid_out,"clen_avg",clen_avg,clen_avg_ID); HANDLE_ERR(err,loc);
  err = nc_def_dim(ncid_out,"clen_edgV",clen_edgV,clen_edgV_ID); HANDLE_ERR(err,loc);
  err = nc_def_dim(ncid_out,"clen_edgH",clen_edgH,clen_edgH_ID); HANDLE_ERR(err,loc);

  return err;
}

//================================================================================================================================//
//
//================================================================================================================================//
int DEF_POD_VARS(int ncid, int gsum, char *vname, int N_g_ID, int rank_ID, int clen_ID, int N_t_ID, int *Sid, int *Uid, int *Vtid)
{
  int err;
  char loc[13] = "DEF_POD_VARS";
  size_t ndims;
  int *dimids, p1, p2;
  char name[128];

  if (gsum == 1){ ndims = 2; p1 = 0; p2 = 1;}
  else{ ndims = 3; p1 = 1; p2 = 2;}

  dimids = (int *)malloc(sizeof(int)*ndims);
  dimids[0] = N_g_ID;

  dimids[p1] = rank_ID;
  memset(name,0,128); strcpy(name,"S_"); strcat(name,vname);
  err = nc_def_var(ncid,name,NC_DOUBLE,(int)ndims-1,dimids,Sid); HANDLE_ERR(err,loc);

  dimids[p1] = rank_ID; dimids[p2] = clen_ID;
  memset(name,0,128); strcpy(name,"U_"); strcat(name,vname);
  err = nc_def_var(ncid,name,NC_DOUBLE,(int)ndims,dimids,Uid); HANDLE_ERR(err,loc);

  dimids[p1] = N_t_ID; dimids[p2] = rank_ID;
  memset(name,0,128); strcpy(name,"Vt_"); strcat(name,vname);
  err = nc_def_var(ncid,name,NC_DOUBLE,(int)ndims,dimids,Vtid); HANDLE_ERR(err,loc);

  free(dimids);

  return 0;

}

//================================================================================================================================//
//
//================================================================================================================================//
int DEF_VARS(int ncid, int gsum, int N_g_ID, int rank_avg_ID, int rank_edgV_ID, int rank_edgH_ID, int clen_avg_ID, int clen_edgV_ID,
   int clen_edgH_ID, int N_t_ID, int *S_fg_avg_xx_ID, int *U_fg_avg_xx_ID, int *Vt_fg_avg_xx_ID, int *S_fg_edgV_xx_ID,
   int *U_fg_edgV_xx_ID, int *Vt_fg_edgV_xx_ID, int *S_fg_avg_yy_ID, int *U_fg_avg_yy_ID, int *Vt_fg_avg_yy_ID,
   int *S_fg_edgH_yy_ID, int *U_fg_edgH_yy_ID, int *Vt_fg_edgH_yy_ID, int *S_fg_edgV_xy_ID, int *U_fg_edgV_xy_ID,
   int *Vt_fg_edgV_xy_ID, int *S_fg_edgH_xy_ID, int *U_fg_edgH_xy_ID, int *Vt_fg_edgH_xy_ID)
{
  int err;

  err = DEF_POD_VARS(ncid,gsum,"fg_avg_xx",N_g_ID,rank_avg_ID,clen_avg_ID,N_t_ID,&(*S_fg_avg_xx_ID),&(*U_fg_avg_xx_ID),&(*Vt_fg_avg_xx_ID));
  err = DEF_POD_VARS(ncid,gsum,"fg_edgV_xx",N_g_ID,rank_edgV_ID,clen_edgV_ID,N_t_ID,&(*S_fg_edgV_xx_ID),&(*U_fg_edgV_xx_ID),&(*Vt_fg_edgV_xx_ID));
  err = DEF_POD_VARS(ncid,gsum,"fg_avg_yy",N_g_ID,rank_avg_ID,clen_avg_ID,N_t_ID,&(*S_fg_avg_yy_ID),&(*U_fg_avg_yy_ID),&(*Vt_fg_avg_yy_ID));
  err = DEF_POD_VARS(ncid,gsum,"fg_edgH_yy",N_g_ID,rank_edgH_ID,clen_edgH_ID,N_t_ID,&(*S_fg_edgH_yy_ID),&(*U_fg_edgH_yy_ID),&(*Vt_fg_edgH_yy_ID));
  err = DEF_POD_VARS(ncid,gsum,"fg_edgV_xy",N_g_ID,rank_edgV_ID,clen_edgV_ID,N_t_ID,&(*S_fg_edgV_xy_ID),&(*U_fg_edgV_xy_ID),&(*Vt_fg_edgV_xy_ID));
  err = DEF_POD_VARS(ncid,gsum,"fg_edgH_xy",N_g_ID,rank_edgH_ID,clen_edgH_ID,N_t_ID,&(*S_fg_edgH_xy_ID),&(*U_fg_edgH_xy_ID),&(*Vt_fg_edgH_xy_ID));

}

//================================================================================================================================//
//
//================================================================================================================================//
int OUTPUT_fg_POD(int ncid_in, int ncid_out, size_t N_t, size_t N_g, size_t N_y, size_t N_x, int gsum, int S_fg_avg_xx_ID,
  int U_fg_avg_xx_ID, int Vt_fg_avg_xx_ID, int S_fg_edgV_xx_ID, int U_fg_edgV_xx_ID, int Vt_fg_edgV_xx_ID,
  int S_fg_avg_yy_ID, int U_fg_avg_yy_ID, int Vt_fg_avg_yy_ID, int S_fg_edgH_yy_ID, int U_fg_edgH_yy_ID,
  int Vt_fg_edgH_yy_ID, int S_fg_edgV_xy_ID, int U_fg_edgV_xy_ID, int Vt_fg_edgV_xy_ID, int S_fg_edgH_xy_ID,
  int U_fg_edgH_xy_ID, int Vt_fg_edgH_xy_ID)
{
  int err;
  char dname[25];
  size_t clen, rank, n_g, gscale;

  if (gsum == 1){
    n_g = 0; gscale = N_g;
  }
  else{
    n_g = N_g; gscale = 1;
  }

  memset(dname,0,25); strcpy(dname,"fg_avg_xx"); clen = gscale*N_y*N_x; rank = min(clen,N_t);
  err = GENERATE_POD(ncid_in,ncid_out,dname,N_t,n_g,clen,rank,gsum,S_fg_avg_xx_ID,U_fg_avg_xx_ID,Vt_fg_avg_xx_ID);

  memset(dname,0,25); strcpy(dname,"fg_edgV_xx"); clen = gscale*N_y*(N_x+1); rank = min(clen,N_t);
  err = GENERATE_POD(ncid_in,ncid_out,dname,N_t,n_g,clen,rank,gsum,S_fg_edgV_xx_ID,U_fg_edgV_xx_ID,Vt_fg_edgV_xx_ID);

  memset(dname,0,25); strcpy(dname,"fg_avg_yy"); clen = gscale*N_y*N_x; rank = min(clen,N_t);
  err = GENERATE_POD(ncid_in,ncid_out,dname,N_t,n_g,clen,rank,gsum,S_fg_avg_yy_ID,U_fg_avg_yy_ID,Vt_fg_avg_yy_ID);

  memset(dname,0,25); strcpy(dname,"fg_edgH_yy\0"); clen = gscale*(N_y+1)*N_x; rank = min(clen,N_t);
  err = GENERATE_POD(ncid_in,ncid_out,dname,N_t,n_g,clen,rank,gsum,S_fg_edgH_yy_ID,U_fg_edgH_yy_ID,Vt_fg_edgH_yy_ID);

  memset(dname,0,25); strcpy(dname,"fg_edgV_xy\0"); clen = gscale*N_y*(N_x+1); rank = min(clen,N_t);
  err = GENERATE_POD(ncid_in,ncid_out,dname,N_t,n_g,clen,rank,gsum,S_fg_edgV_xy_ID,U_fg_edgV_xy_ID,Vt_fg_edgV_xy_ID);

  memset(dname,0,25); strcpy(dname,"fg_edgH_xy\0"); clen = gscale*(N_y+1)*N_x; rank = min(clen,N_t);
  err = GENERATE_POD(ncid_in,ncid_out,dname,N_t,n_g,clen,rank,gsum,S_fg_edgH_xy_ID,U_fg_edgH_xy_ID,Vt_fg_edgH_xy_ID);

  return err;
}
