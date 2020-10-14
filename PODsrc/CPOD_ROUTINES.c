#include <stdlib.h>
#include <stdio.h>
#include <string.h>

//FROM NCDF_IO.c
void HANDLE_ERR(int Status, char Location[]);
void GET_VAR_DOUBLE(int ncid, char name[], double **var, size_t size);

//FROM POD_ROUTINES.f08
void SVD_CALC(const double *, const int *n_t, const int *n_y, const int *n_x, double *center, double *umat, double *sig, double *vtmat);

//FROM GNUPLOT_ROUTINES.c
int gnuplot_1d(char *title, double *data, double *crd, int dim, char *plttyp, int logscale, int pt, char *lc, char *saveas);
void GNUP_ERR(int err);

//LOCAL FUNCTIONS
int POD_CALC(int ncid, char *dname, size_t len, size_t N_t, size_t N_g, size_t N_y, size_t N_x, int g, double **center,
  double **umat, double **sig, double **vtmat, double **sigp, size_t *rank);
int gdat_reform(int n_t, int n_g, int n_y, int n_x, int g, double *gdat, double *vec);

#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })

//================================================================================================================================//
// int GENERATE_POD
//
// GENERATE_POD
//================================================================================================================================//
int GENERATE_POD(int ncid, char *dname, size_t len, size_t N_t, size_t N_g, size_t N_y, size_t N_x)
{
  int err, i;
  size_t rank;
  char lc[10], plttyp[3], title[100], saveas[100];
  double *center, *umat, *sig, *vtmat, *sigp;

  err = POD_CALC(ncid,dname,len,N_t,N_g,N_y,N_x,0,&center,&umat,&sig,&vtmat,&sigp,&rank);

  strcpy(lc,"black"); //color of the plotted data
  strcpy(plttyp,"p"); //type of plot
  strcpy(title,dname); strcat(title,"_Singular_Values\0"); //create plot title
  strcpy(saveas,title); strcat(saveas,".png\0"); //create saved name of plot

  //replace all underscores in title as spaces
  i=0;
  while(title[i] != '\0'){
    if(title[i] == '_'){
      title[i] = ' ';
    }
    i++;
  }

  err = gnuplot_1d(title,sig,sigp,(int)rank,plttyp,10,7,lc,saveas); GNUP_ERR(err);

  //deallocating arrays
  free(center); free(umat); free(sig); free(vtmat); free(sigp);

  return 0;
}

//================================================================================================================================//
// int POD_CALC
//
// POD_CALC calculates the POD of a data-matrix read from a given NetCDF dataset
//================================================================================================================================//
int POD_CALC(int ncid, char *dname, size_t len, size_t N_t, size_t N_g, size_t N_y, size_t N_x, int g, double **center,
  double **umat, double **sig, double **vtmat, double **sigp, size_t *rank)
{
  double *data, *temp;
  int n_t, n_g, n_y, n_x;
  int i, err;

  //finding integer versions of the given dimensions
  n_t = (int)N_t; n_g = (int)N_g; n_y = (int)N_y; n_x = (int)N_x;

  //reading in datamatrix from NetCDF dataset
  GET_VAR_DOUBLE(ncid,dname,&data,len);

  //allocating POD arrays
  *rank = min(N_x*N_y,N_t); //calculating rank
  *center = (double *)malloc(sizeof(double)*N_y*N_x);
  *umat = (double *)malloc(sizeof(double)*N_y*N_x*(*rank));
  *vtmat = (double *)malloc(sizeof(double)*N_t*(*rank));
  *sig = (double *)malloc(sizeof(double)*(*rank));

  //forming array of singular value indices
  *sigp = (double *)malloc(sizeof(double)*(*rank));
  for(i=0;i<(int)(*rank);i++){
    (*sigp)[i] = (double)(i+1);
  }

  //checking whether a multigroup data-matrix is being used or not
  if(n_g>0){ //if the datamatrix is multigroup, must either isolate one group or form a large multigroup matrix
    //reforming the datamatrix to isolate a single group
    temp = (double *)malloc(sizeof(double)*N_t*N_y*N_x);
    err = gdat_reform(n_t,n_g,n_y,n_x,g,data,temp);

    //calculating the SVD of the datamatrix
    SVD_CALC(temp,&n_t,&n_y,&n_x,&*center[0],&*umat[0],&*sig[0],&*vtmat[0]);
  }
  else{ //if the datamatrix is not multigroup, can immediately procede with finding the SVD
    //calculating the SVD of the datamatrix
    SVD_CALC(data,&n_t,&n_y,&n_x,&*center[0],&*umat[0],&*sig[0],&*vtmat[0]);
  }

  //deallocating arrays
  free(temp); free(data);

  return 0;
}

//================================================================================================================================//
//
//================================================================================================================================//
int gdat_reform(int n_t, int n_g, int n_y, int n_x, int g, double *gdat, double *vec)
{
  int i, j, t;

  for(t=0;t<n_t;t++){
    for(j=0;j<n_y;j++){
      for(i=0;i<n_x;i++){
        vec[i+j*n_x+t*n_x*n_y] = gdat[i+j*n_x+g*n_x*n_y+t*n_x*n_y*n_g];
      }
    }
  }

  return 0;

}
