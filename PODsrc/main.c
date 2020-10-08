#include <stdlib.h>
#include <stdio.h>
#include <netcdf.h>
#include <string.h>

// int TEST(int a, int b);

void SVD_CALC(const double *, const int *, const int *, const int *, double *umat, double *sig, double *vtmat);

void INPUT(char *infile, char *dsfile);
void copy(char to[], char from[]);
void GET_DIMS(int ncid, size_t *N_t, size_t *N_g, size_t *N_y, size_t *N_x);
void GET_VAR_DOUBLE(int ncid, char name[], double **var, size_t size);

/* Handle NetCDF errors by printing an error message and exiting with a
 * non-zero status. */
// #define ERRCODE 2
#define NC_ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(2);}
#define NC_NOERR   0

//================================================================================================================================//
//
//================================================================================================================================//
int main()
{
  // comment
  int ncid, err;
  int fg_avg_xx_ID, fg_avg_xy_ID, fg_avg_yy_ID, fg_edgV_xx_ID, fg_edgV_xy_ID, fg_edgH_yy_ID, fg_edgH_xy_ID;
  char infile[9], dsfile[100];
  size_t N_t, N_g, N_y, N_x, len, rank;
  int n_t, n_g, n_y, n_x, i, j, g, t;
  double *fg_avg_xx, *temp, *umat, *sig, *vtmat;
  // int t, g, j, i;

  copy(infile,"input.inp"); //setting input file name (default)

  INPUT(infile,dsfile); //reading input file
  printf("%s\n",dsfile);

  err = nc_open(dsfile,NC_NOWRITE,&ncid); //opening NetCDF dataset
  if(err != NC_NOERR) NC_ERR(err)

  GET_DIMS(ncid,&N_t,&N_g,&N_y,&N_x);
  printf("%ld %ld %ld %ld \n",N_t,N_g,N_y,N_x);
  n_t = (int)N_t; n_g = (int)N_g; n_y = (int)N_y; n_x = (int)N_x;

  //reading in multigroup qd factor ID's
  // err = nc_inq_varid(ncid,"fg_avg_xy",&fg_avg_xy_ID); if(err != NC_NOERR) NC_ERR(err);
  // err = nc_inq_varid(ncid,"fg_avg_yy",&fg_avg_yy_ID); if(err != NC_NOERR) NC_ERR(err);
  // err = nc_inq_varid(ncid,"fg_edgV_xx",&fg_edgV_xx_ID); if(err != NC_NOERR) NC_ERR(err);
  // err = nc_inq_varid(ncid,"fg_edgV_xy",&fg_edgV_xy_ID); if(err != NC_NOERR) NC_ERR(err);
  // err = nc_inq_varid(ncid,"fg_edgH_yy",&fg_edgH_yy_ID); if(err != NC_NOERR) NC_ERR(err);
  // err = nc_inq_varid(ncid,"fg_edgH_xy",&fg_edgH_xy_ID); if(err != NC_NOERR) NC_ERR(err);

  len = N_t*N_g*N_y*N_x;
  //reading in multigroup qd factor data
  GET_VAR_DOUBLE(ncid,"fg_avg_xx",&fg_avg_xx,len);

  rank = N_t; //min(N_t,N_y*N_x);
  umat = (double *)malloc(sizeof(double)*N_y*N_x*rank);
  vtmat = (double *)malloc(sizeof(double)*N_t*rank);
  sig = (double *)malloc(sizeof(double)*rank);

  temp = (double *)malloc(sizeof(double)*N_t*N_y*N_x);
  for(t=0;t<n_t;t++){
    for(j=0;j<n_y;j++){
      for(i=0;i<n_x;i++){
        temp[i+j*n_x+t*n_x*n_y] = fg_avg_xx[i+j*n_x+t*n_x*n_y*n_g];
      }
    }
  }
  // printf("yuh\n");
  SVD_CALC(temp,&n_t,&n_y,&n_x,&umat[0],&sig[0],&vtmat[0]);
  printf("%e\n",fg_avg_xx[100]);

  err = nc_close(ncid); //closing NetCDF dataset
  if(err != NC_NOERR) NC_ERR(err);

  // free(fg_avg_xx);

  printf("yuh\n");
}
