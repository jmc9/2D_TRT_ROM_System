#include <stdlib.h>
#include <stdio.h>
#include <netcdf.h>
#include <string.h>

// int TEST(int a, int b);

void SVD_CALC(const double *, const int *n_t, const int *n_y, const int *n_x, double *center, double *umat, double *sig, double *vtmat);

void INPUT(char *infile, char *dsfile);
void copy(char to[], char from[]);
void GET_DIMS(int ncid, size_t *N_t, size_t *N_g, size_t *N_y, size_t *N_x);
void GET_VAR_DOUBLE(int ncid, char name[], double **var, size_t size);

int gnuplot_1d(char *title, double *data, double *crd, int dim, char *plttyp, int logscale, int pt, char *lc);

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
  char infile[9], dsfile[100], lc[10], plttyp[1], title[10];
  size_t N_t, N_g, N_y, N_x, len, rank;
  int n_t, n_g, n_y, n_x, i, j, g, t;
  double *fg_avg_xx, *temp, *center, *umat, *sig, *vtmat, *sigp;
  // int t, g, j, i;
  // FILE *gnuplot_Pipe;

  copy(infile,"input.inp"); //setting input file name (default)

  INPUT(infile,dsfile); //reading input file
  printf("%s\n",dsfile);

  err = nc_open(dsfile,NC_NOWRITE,&ncid); //opening NetCDF dataset
  if(err != NC_NOERR) NC_ERR(err);

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
  center = (double *)malloc(sizeof(double)*N_y*N_x);
  umat = (double *)malloc(sizeof(double)*N_y*N_x*rank);
  vtmat = (double *)malloc(sizeof(double)*N_t*rank);
  sig = (double *)malloc(sizeof(double)*rank);

  sigp = (double *)malloc(sizeof(double)*rank);
  for(i=0;i<rank;i++){
    sigp[i] = i+1;
  }

  temp = (double *)malloc(sizeof(double)*N_t*N_y*N_x);
  for(t=0;t<n_t;t++){
    for(j=0;j<n_y;j++){
      for(i=0;i<n_x;i++){
        temp[i+j*n_x+t*n_x*n_y] = fg_avg_xx[i+j*n_x+t*n_x*n_y*n_g];
      }
    }
  }
  SVD_CALC(temp,&n_t,&n_y,&n_x,&center[0],&umat[0],&sig[0],&vtmat[0]);
  printf("%e\n",fg_avg_xx[100]);

  err = nc_close(ncid); //closing NetCDF dataset
  if(err != NC_NOERR) NC_ERR(err);

  copy(lc,"black"); copy(plttyp,"p"); copy(title,"Sigs");
  err = gnuplot_1d(title,&sig[0],&sigp[0],rank,&plttyp,10,7,&lc);
  if(err != 0){
    if(err == 1){ printf("Failed to open gnuplot_Pipe"); exit(2); }
    if(err == 2){ printf("Failed to close gnuplot_Pipe"); exit(2); }
  }

  free(fg_avg_xx);

  printf("yuh\n");
}
