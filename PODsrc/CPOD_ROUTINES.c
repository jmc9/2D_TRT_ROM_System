#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <netcdf.h>

//FROM NCDF_IO.c
void HANDLE_ERR(int Status, char Location[]);
void GET_VAR_DOUBLE(int ncid, char name[], double **var, size_t size);

//FROM POD_ROUTINES.f08
void SVD_CALC(const double *, const int *n_t, const int *clen, double *center, double *umat, double *sig, double *vtmat);

//FROM GNUPLOT_ROUTINES.c
int gnuplot_1d(char *title, double *data, double *crd, int dim, char *plttyp, int logscale, int pt, char *lc, char *saveas);
void GNUP_ERR(int err);

//FROM MISC_PROCS.c
int make_dir(char *dir);

//LOCAL FUNCTIONS
int POD_CALC(double *data, size_t N_t, size_t N_g, size_t N_y, size_t N_x, int g, double **center, double **umat,
  double **sig, double **vtmat, double **sigp, size_t *rank);
int gdat_reform(int n_t, int n_g, int n_y, int n_x, int g, double *gdat, double *vec);

int SIG_PLOT(char *pname, double *sig, double *sigp, int rank, char *drop);



#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })

//================================================================================================================================//
// int GENERATE_POD
//
// GENERATE_POD
//================================================================================================================================//
int GENERATE_POD(int ncid_in, int ncid_out, char *dname, size_t len, size_t N_t, size_t N_g, size_t N_y, size_t N_x, int gsum,
  int Sid, int Uid, int Vtid)
{
  int err, n_g, g;
  size_t rank;
  char buf[25], pname[25], drop[25];
  double *data, *center, *umat, *sig, *vtmat, *sigp;

  size_t startp[2], countp[2];
  ptrdiff_t stridep[2];

  //creating directory to put plots
  err = make_dir(dname); //setting up directory
  strcpy(drop,dname); strcat(drop,"/"); //creating path to directory

  //reading in datamatrix from NetCDF dataset
  GET_VAR_DOUBLE(ncid_in,dname,&data,len);

  n_g = (int)N_g; //integer version of group dimension (number of energy groups)

  //checking type of dataset to perform POD on
  if(n_g > 0){ //if n_g>0, then a multigroup dataset has been detected


    if(gsum == 1){ //gsum == 1 means to perform POD on the entire multigroup dataset at once
      //large singule multigroup matrix will go here
    }
    else{ //gsum != 1 means to perform POD on each group of the dataset individually

      //loop over energy groups
      for(g=0;g<n_g;g++){

        //find the POD modes and singular values of a given groupwise datamatrix
        err = POD_CALC(data,N_t,N_g,N_y,N_x,g,&center,&umat,&sig,&vtmat,&sigp,&rank);
        // printf("%lu\n",rank) ;

        startp[0] = (size_t)g; startp[1] = 0;
        countp[0] = 1; countp[1] = rank;
        stridep[0] = 1; stridep[1] = 1;
        err = nc_put_vars(ncid_out,Sid,startp,countp,stridep,sig);

        //plot the singular values
        strcpy(pname,dname); sprintf(buf,"_g%d",g+1); strcat(pname,buf);
        err = SIG_PLOT(pname,sig,sigp,(int)rank,drop);

      } //end g loop
    }
  }
  else{
    err = POD_CALC(data,N_t,N_g,N_y,N_x,0,&center,&umat,&sig,&vtmat,&sigp,&rank);
  }

  //deallocating arrays
  free(data); free(center); free(umat); free(sig); free(vtmat); free(sigp);

  return 0;
}

//================================================================================================================================//
// int SIG_PLOT
//
// SIG_PLOT
//================================================================================================================================//
int SIG_PLOT(char *pname, double *sig, double *sigp, int rank, char *drop)
{
  int i, err;
  char lc[10], plttyp[3], title[100], saveas[100];

  strcpy(lc,"black"); //color of the plotted data
  strcpy(plttyp,"p"); //type of plot
  strcpy(title,pname); strcat(title,"_Singular_Values\0"); //create plot title
  strcpy(saveas,drop); strcat(saveas,title); strcat(saveas,".png\0"); //create saved name of plot

  //replace all underscores in title as spaces
  i=0;
  while(title[i] != '\0'){
    if(title[i] == '_'){
      title[i] = ' ';
    }
    i++;
  }

  err = gnuplot_1d(title,sig,sigp,rank,plttyp,10,7,lc,saveas); GNUP_ERR(err);

  return err;
}

//================================================================================================================================//
// int POD_CALC
//
// POD_CALC calculates the POD of a data-matrix read from a given NetCDF dataset
//================================================================================================================================//
int POD_CALC(double *data, size_t N_t, size_t N_g, size_t N_y, size_t N_x, int g, double **center,
  double **umat, double **sig, double **vtmat, double **sigp, size_t *rank)
{
  double *temp;
  int n_t, n_g, n_y, n_x, len_xy;
  int i, err;

  //finding integer versions of the given dimensions
  n_t = (int)N_t; n_g = (int)N_g; n_y = (int)N_y; n_x = (int)N_x;
  len_xy = n_x*n_y;

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
    SVD_CALC(temp,&n_t,&len_xy,&*center[0],&*umat[0],&*sig[0],&*vtmat[0]);
  }
  else{ //if the datamatrix is not multigroup, can immediately procede with finding the SVD
    //calculating the SVD of the datamatrix
    SVD_CALC(data,&n_t,&len_xy,&*center[0],&*umat[0],&*sig[0],&*vtmat[0]);
  }

  //deallocating arrays
  free(temp);

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
