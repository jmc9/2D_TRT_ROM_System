#include <stdlib.h>
#include <stdio.h>
#include <netcdf.h>
#include <string.h>

//FROM NCDF_IO.c
void HANDLE_ERR(const int Status, const char *Location);

//FROM CPOD_ROUTINES.c
int GENERATE_POD(const int ncid_in, const int ncid_out, const char *dname, const size_t N_t, const size_t N_g, const size_t clen,
  size_t rank, int Cid, int Sid, int Uid, int Vtid);

//LOCAL FUNCTIONS
#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })

//================================================================================================================================//
//
//================================================================================================================================//
int DEF_DIMS(const int ncid_out, const size_t N_t, const size_t N_g, const size_t N_m, const size_t N_y, const size_t N_x,
   const size_t rank_avg, const size_t rank_edgV,const size_t rank_edgH, const size_t clen_avg, const size_t clen_edgV,
   const size_t clen_edgH, int *N_t_ID, int *N_g_ID, int *N_m_ID, int *N_y_ID,int *N_x_ID, int *rank_avg_ID, int *rank_edgV_ID,
   int *rank_edgH_ID, int *clen_avg_ID, int *clen_edgV_ID, int *clen_edgH_ID)
{
  int err;
  char loc[9] = "OUT_DIMS";

  err = nc_def_dim(ncid_out,"N_x",N_x,N_x_ID); HANDLE_ERR(err,loc);
  err = nc_def_dim(ncid_out,"N_y",N_y,N_y_ID); HANDLE_ERR(err,loc);
  err = nc_def_dim(ncid_out,"N_m",N_m,N_m_ID); HANDLE_ERR(err,loc);
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
int DEF_POD_VARS(const int ncid, const int gsum, const char *vname, const int N_g_ID, const int rank_ID, const int clen_ID,
  const int N_t_ID, int *Cid, int *Sid, int *Uid, int *Vtid)
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

  dimids[p1] = clen_ID;
  memset(name,0,128); strcpy(name,"C_"); strcat(name,vname);
  err = nc_def_var(ncid,name,NC_DOUBLE,(int)ndims-1,dimids,Cid); HANDLE_ERR(err,loc);

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
int DEF_fg_VARS(const int ncid, const int gsum, const int N_g_ID, const int rank_avg_ID, const int rank_edgV_ID,
   const int rank_edgH_ID, const int clen_avg_ID,const int clen_edgV_ID, const int clen_edgH_ID, const int N_t_ID,
   int *C_fg_avg_xx_ID, int *S_fg_avg_xx_ID, int *U_fg_avg_xx_ID,int *Vt_fg_avg_xx_ID, int *C_fg_edgV_xx_ID,
   int *S_fg_edgV_xx_ID, int *U_fg_edgV_xx_ID, int *Vt_fg_edgV_xx_ID, int *C_fg_avg_yy_ID, int *S_fg_avg_yy_ID,
   int *U_fg_avg_yy_ID, int *Vt_fg_avg_yy_ID,int *C_fg_edgH_yy_ID, int *S_fg_edgH_yy_ID, int *U_fg_edgH_yy_ID,
   int *Vt_fg_edgH_yy_ID, int *C_fg_edgV_xy_ID, int *S_fg_edgV_xy_ID, int *U_fg_edgV_xy_ID, int *Vt_fg_edgV_xy_ID,
   int *C_fg_edgH_xy_ID, int *S_fg_edgH_xy_ID, int *U_fg_edgH_xy_ID, int *Vt_fg_edgH_xy_ID)
{
  int err;

  err = DEF_POD_VARS(ncid,gsum,"fg_avg_xx",N_g_ID,rank_avg_ID,clen_avg_ID,N_t_ID,&(*C_fg_avg_xx_ID),&(*S_fg_avg_xx_ID),&(*U_fg_avg_xx_ID),&(*Vt_fg_avg_xx_ID));
  err = DEF_POD_VARS(ncid,gsum,"fg_edgV_xx",N_g_ID,rank_edgV_ID,clen_edgV_ID,N_t_ID,&(*C_fg_edgV_xx_ID),&(*S_fg_edgV_xx_ID),&(*U_fg_edgV_xx_ID),&(*Vt_fg_edgV_xx_ID));
  err = DEF_POD_VARS(ncid,gsum,"fg_avg_yy",N_g_ID,rank_avg_ID,clen_avg_ID,N_t_ID,&(*C_fg_avg_yy_ID),&(*S_fg_avg_yy_ID),&(*U_fg_avg_yy_ID),&(*Vt_fg_avg_yy_ID));
  err = DEF_POD_VARS(ncid,gsum,"fg_edgH_yy",N_g_ID,rank_edgH_ID,clen_edgH_ID,N_t_ID,&(*C_fg_edgH_yy_ID),&(*S_fg_edgH_yy_ID),&(*U_fg_edgH_yy_ID),&(*Vt_fg_edgH_yy_ID));
  err = DEF_POD_VARS(ncid,gsum,"fg_edgV_xy",N_g_ID,rank_edgV_ID,clen_edgV_ID,N_t_ID,&(*C_fg_edgV_xy_ID),&(*S_fg_edgV_xy_ID),&(*U_fg_edgV_xy_ID),&(*Vt_fg_edgV_xy_ID));
  err = DEF_POD_VARS(ncid,gsum,"fg_edgH_xy",N_g_ID,rank_edgH_ID,clen_edgH_ID,N_t_ID,&(*C_fg_edgH_xy_ID) ,&(*S_fg_edgH_xy_ID),&(*U_fg_edgH_xy_ID),&(*Vt_fg_edgH_xy_ID));

  return err;

}

//================================================================================================================================//
//
//================================================================================================================================//
int DEF_Ig_VARS(const int ncid, const int gsum, const int N_g_ID, const int rank_avg_ID, const int rank_edgV_ID,
   const int rank_edgH_ID, const int clen_avg_ID,const int clen_edgV_ID, const int clen_edgH_ID, const int N_t_ID,
   int *C_Ig_avg_ID, int *S_Ig_avg_ID, int *U_Ig_avg_ID, int *Vt_Ig_avg_ID,int *C_Ig_edgV_ID, int *S_Ig_edgV_ID,
   int *U_Ig_edgV_ID, int *Vt_Ig_edgV_ID, int *C_Ig_edgH_ID, int *S_Ig_edgH_ID, int *U_Ig_edgH_ID, int *Vt_Ig_edgH_ID)
{
  int err;

  err = DEF_POD_VARS(ncid,gsum,"Ig_avg",N_g_ID,rank_avg_ID,clen_avg_ID,N_t_ID,&(*C_Ig_avg_ID),&(*S_Ig_avg_ID),&(*U_Ig_avg_ID),&(*Vt_Ig_avg_ID));
  err = DEF_POD_VARS(ncid,gsum,"Ig_edgV",N_g_ID,rank_edgV_ID,clen_edgV_ID,N_t_ID,&(*C_Ig_edgV_ID),&(*S_Ig_edgV_ID),&(*U_Ig_edgV_ID),&(*Vt_Ig_edgV_ID));
  err = DEF_POD_VARS(ncid,gsum,"Ig_edgH",N_g_ID,rank_edgH_ID,clen_edgH_ID,N_t_ID,&(*C_Ig_edgH_ID),&(*S_Ig_edgH_ID),&(*U_Ig_edgH_ID),&(*Vt_Ig_edgH_ID));

  return err;

}

//================================================================================================================================//
//
//================================================================================================================================//
int OUTPUT_fg_POD(const int ncid_in, const int ncid_out, const size_t N_t, const size_t N_g, const size_t N_y, const size_t N_x,
  const int gsum, const int C_fg_avg_xx_ID,const int S_fg_avg_xx_ID, const int U_fg_avg_xx_ID, const int Vt_fg_avg_xx_ID,
  const int C_fg_edgV_xx_ID, const int S_fg_edgV_xx_ID, const int U_fg_edgV_xx_ID,const int Vt_fg_edgV_xx_ID,
  const int C_fg_avg_yy_ID, const int S_fg_avg_yy_ID, const int U_fg_avg_yy_ID, const int Vt_fg_avg_yy_ID,
  const int C_fg_edgH_yy_ID, const int S_fg_edgH_yy_ID, const int U_fg_edgH_yy_ID, const int Vt_fg_edgH_yy_ID,
  const int C_fg_edgV_xy_ID, const int S_fg_edgV_xy_ID, const int U_fg_edgV_xy_ID, const int Vt_fg_edgV_xy_ID,
  const int C_fg_edgH_xy_ID, const int S_fg_edgH_xy_ID, const int U_fg_edgH_xy_ID, const int Vt_fg_edgH_xy_ID)
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

  printf("... fg_avg_xx start\n");
  memset(dname,0,25); strcpy(dname,"fg_avg_xx"); clen = gscale*N_y*N_x; rank = min(clen,N_t);
  err = GENERATE_POD(ncid_in,ncid_out,dname,N_t,n_g,clen,rank,C_fg_avg_xx_ID,S_fg_avg_xx_ID,U_fg_avg_xx_ID,Vt_fg_avg_xx_ID);

  printf("... fg_edgV_xx start\n");
  memset(dname,0,25); strcpy(dname,"fg_edgV_xx"); clen = gscale*N_y*(N_x+1); rank = min(clen,N_t);
  err = GENERATE_POD(ncid_in,ncid_out,dname,N_t,n_g,clen,rank,C_fg_edgV_xx_ID,S_fg_edgV_xx_ID,U_fg_edgV_xx_ID,Vt_fg_edgV_xx_ID);

  printf("... fg_avg_yy start\n");
  memset(dname,0,25); strcpy(dname,"fg_avg_yy"); clen = gscale*N_y*N_x; rank = min(clen,N_t);
  err = GENERATE_POD(ncid_in,ncid_out,dname,N_t,n_g,clen,rank,C_fg_avg_yy_ID,S_fg_avg_yy_ID,U_fg_avg_yy_ID,Vt_fg_avg_yy_ID);

  printf("... fg_edgH_yy start\n");
  memset(dname,0,25); strcpy(dname,"fg_edgH_yy"); clen = gscale*(N_y+1)*N_x; rank = min(clen,N_t);
  err = GENERATE_POD(ncid_in,ncid_out,dname,N_t,n_g,clen,rank,C_fg_edgH_yy_ID,S_fg_edgH_yy_ID,U_fg_edgH_yy_ID,Vt_fg_edgH_yy_ID);

  printf("... fg_edgV_xy start\n");
  memset(dname,0,25); strcpy(dname,"fg_edgV_xy"); clen = gscale*N_y*(N_x+1); rank = min(clen,N_t);
  err = GENERATE_POD(ncid_in,ncid_out,dname,N_t,n_g,clen,rank,C_fg_edgV_xy_ID,S_fg_edgV_xy_ID,U_fg_edgV_xy_ID,Vt_fg_edgV_xy_ID);

  printf("... fg_edgH_xy start\n");
  memset(dname,0,25); strcpy(dname,"fg_edgH_xy"); clen = gscale*(N_y+1)*N_x; rank = min(clen,N_t);
  err = GENERATE_POD(ncid_in,ncid_out,dname,N_t,n_g,clen,rank,C_fg_edgH_xy_ID,S_fg_edgH_xy_ID,U_fg_edgH_xy_ID,Vt_fg_edgH_xy_ID);

  return err;
}

//================================================================================================================================//
//
//================================================================================================================================//
int OUTPUT_Ig_POD(const int ncid_in, const int ncid_out, const size_t N_t, const size_t N_g, const size_t N_m, const size_t N_y,
  const size_t N_x, const int gsum, const int C_Ig_avg_ID, const int S_Ig_avg_ID, const int U_Ig_avg_ID, const int Vt_Ig_avg_ID,
  const int C_Ig_edgV_ID, const int S_Ig_edgV_ID, const int U_Ig_edgV_ID, const int Vt_Ig_edgV_ID, const int C_Ig_edgH_ID,
  const int S_Ig_edgH_ID, const int U_Ig_edgH_ID, const int Vt_Ig_edgH_ID)
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

  printf("... I_avg start\n");
  memset(dname,0,25); strcpy(dname,"I_avg"); clen = gscale*N_m*N_y*N_x; rank = min(clen,N_t);
  err = GENERATE_POD(ncid_in,ncid_out,dname,N_t,n_g,clen,rank,C_Ig_avg_ID,S_Ig_avg_ID,U_Ig_avg_ID,Vt_Ig_avg_ID);

  printf("... I_edgV start\n");
  memset(dname,0,25); strcpy(dname,"I_edgV"); clen = gscale*N_m*N_y*(N_x+1); rank = min(clen,N_t);
  err = GENERATE_POD(ncid_in,ncid_out,dname,N_t,n_g,clen,rank,C_Ig_edgV_ID,S_Ig_edgV_ID,U_Ig_edgV_ID,Vt_Ig_edgV_ID);

  printf("... I_edgH start\n");
  memset(dname,0,25); strcpy(dname,"I_edgH"); clen = gscale*N_m*(N_y+1)*N_x; rank = min(clen,N_t);
  err = GENERATE_POD(ncid_in,ncid_out,dname,N_t,n_g,clen,rank,C_Ig_edgH_ID,S_Ig_edgH_ID,U_Ig_edgH_ID,Vt_Ig_edgH_ID);

  return err;
}
