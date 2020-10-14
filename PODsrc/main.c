#include <stdlib.h>
#include <stdio.h>
#include <netcdf.h>
#include <string.h>
#include <math.h>

// int TEST(int a, int b);

void SVD_CALC(const double *, const int *n_t, const int *n_y, const int *n_x, double *center, double *umat, double *sig, double *vtmat);

void INPUT(char *infile, char *dsfile);
void copy(char to[], char from[]);
void GET_DIMS(int ncid, size_t *N_t, size_t *N_g, size_t *N_y, size_t *N_x);
void GET_VAR_DOUBLE(int ncid, char name[], double **var, size_t size);
void HANDLE_ERR(int Status, char Location[]);

int gnuplot_1d(char *title, double *data, double *crd, int dim, char *plttyp, int logscale, int pt, char *lc);
void GNUP_ERR(int err);

int gdat_reform(int n_t, int n_g, int n_y, int n_x, int g, double *gdat, double *vec);

int POD_CALC(int ncid, char *dname, size_t glen, size_t N_t, size_t N_g, size_t N_y, size_t N_x, int g, double **center,
  double **umat, double **sig, double **vtmat, double **sigp, size_t *rank);
int GENERATE_POD(int ncid, char *dname, size_t len, size_t N_t, size_t N_g, size_t N_y, size_t N_x);

#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

 #define min(a,b) \
    ({ __typeof__ (a) _a = (a); \
        __typeof__ (b) _b = (b); \
      _a < _b ? _a : _b; })

//================================================================================================================================//
//
//================================================================================================================================//
int main()
{
  // comment
  int ncid, err;
  char infile[9], dsfile[100], lc[10], plttyp[1], title[10], dname[25];
  size_t N_t, N_g, N_y, N_x, glen, glen_edgV, glen_edgH, rank;
  int n_t, n_g, n_y, n_x, i, j, g, t;
  double *fg_avg_xx, *temp, *center, *umat, *sig, *vtmat, *sigp;
  char loc[4] = "main";

  copy(infile,"input.inp"); //setting input file name (default)

  INPUT(infile,dsfile); //reading input file
  printf("%s\n",dsfile);

  err = nc_open(dsfile,NC_NOWRITE,&ncid); HANDLE_ERR(err,loc); //opening NetCDF dataset

  GET_DIMS(ncid,&N_t,&N_g,&N_y,&N_x);
  n_t = (int)N_t; n_g = (int)N_g; n_y = (int)N_y; n_x = (int)N_x;
  glen = N_t*N_g*N_y*N_x;
  glen_edgV = N_t*N_g*N_y*(N_x+1);
  glen_edgH = N_t*N_g*(N_y+1)*N_x;

  copy(dname,"fg_avg_xx\0"); err = GENERATE_POD(ncid,dname,glen,N_t,N_g,N_y,N_x);
  copy(dname,"fg_edgV_xx\0"); err = GENERATE_POD(ncid,dname,glen_edgV,N_t,N_g,N_y,(N_x+1));
  copy(dname,"fg_avg_yy\0"); err = GENERATE_POD(ncid,dname,glen,N_t,N_g,N_y,N_x);

  err = nc_close(ncid); HANDLE_ERR(err,loc); //closing NetCDF dataset

  printf("yuh\n");
}

//================================================================================================================================//
//
//================================================================================================================================//
// int

//================================================================================================================================//
//
//================================================================================================================================//
// int gdat_reform(int n_t, int n_g, int n_y, int n_x, int g, double *gdat, double *vec)
// {
//   int i, j, t;
//
//   for(t=0;t<n_t;t++){
//     for(j=0;j<n_y;j++){
//       for(i=0;i<n_x;i++){
//         vec[i+j*n_x+t*n_x*n_y] = gdat[i+j*n_x+g*n_x*n_y+t*n_x*n_y*n_g];
//       }
//     }
//   }
//
//   return 0;
// }
