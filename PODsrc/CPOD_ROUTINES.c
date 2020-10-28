#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <netcdf.h>

//FROM NCDF_IO.c
void HANDLE_ERR(const int Status, const char *Location);
void GET_VAR_DOUBLE(const int ncid, const char *name, double **var, size_t size);

//FROM POD_ROUTINES.f08
void SVD_CALC(const double *, const size_t *n_t, const size_t *clen, double *center, double *umat, double *sig, double *vtmat);

//FROM GNUPLOT_ROUTINES.c
int gnuplot_1d(const char *title, const double *data, const double *crd, const int dim, const char *plttyp, const int logscale,
  const int pt, const char *lc, const char *saveas);
void GNUP_ERR(const int err);

//FROM MISC_PROCS.c
int make_dir(const char *dir);

//LOCAL FUNCTIONS
int POD_CALC(const double *data, const size_t N_t, const size_t N_g, const size_t clen, const size_t rank, const size_t g, double **center,
  double **umat, double **sig, double **vtmat, double **sigp);
int gdat_reform(const size_t N_t, const size_t N_g, const size_t clen, const size_t g, const double *gdat, double *vec);

int SIG_PLOT(const char *pname, const double *sig, const double *sigp, const int rank, const char *drop);

#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })

//================================================================================================================================//
// int GENERATE_POD
//
// GENERATE_POD
//================================================================================================================================//
int GENERATE_POD(const int ncid_in, const int ncid_out, const char *dname, const size_t N_t, const size_t N_g, const size_t clen,
  size_t rank, int Cid, int Sid, int Uid, int Vtid)
{
  int err;
  char loc[13] = "GENERATE_POD";
  size_t len;
  char buf[25], pname[25], drop[25];
  double *data, *center, *umat, *sig, *vtmat, *sigp;

  size_t startp[3], countp[3];
  ptrdiff_t stridep[3];

  //creating directory to put plots
  err = make_dir(dname); //setting up directory
  strcpy(drop,dname); strcat(drop,"/"); //creating path to directory

  if (N_g > 0){
    len = N_g*N_t*clen;
  }
  else{
    len = N_t*clen;
  }

  //reading in datamatrix from NetCDF dataset
  GET_VAR_DOUBLE(ncid_in,dname,&data,len);

  //checking type of dataset to perform POD on
  if(N_g > 0){ //if N_g>0, then a multigroup dataset has been detected
    printf("    -- groupwise decomposition detected\n");

    //loop over energy groups
    for(size_t g=0; g<N_g; g++){
      printf("    -- Start POD on group %lu\n",g+1);

      //find the POD modes and singular values of a given groupwise datamatrix
      err = POD_CALC(data,N_t,N_g,clen,rank,g,&center,&umat,&sig,&vtmat,&sigp);

      //write centering (column-avg'd data) vector to outfile
      startp[0] = (size_t)g; startp[1] = 0; startp[2] = 0;
      countp[0] = 1; countp[1] = clen; countp[2] = 0;
      stridep[0] = 1; stridep[1] = 1; stridep[2] = 0;
      err = nc_put_vars(ncid_out,Cid,startp,countp,stridep,center); HANDLE_ERR(err,loc);

      //write singular value vector to outfile
      startp[0] = (size_t)g; startp[1] = 0; startp[2] = 0;
      countp[0] = 1; countp[1] = rank; countp[2] = 0;
      stridep[0] = 1; stridep[1] = 1; stridep[2] = 0;
      err = nc_put_vars(ncid_out,Sid,startp,countp,stridep,sig); HANDLE_ERR(err,loc);

      //write left singular vector matrix to outfile
      startp[0] = (size_t)g; startp[1] = 0; startp[2] = 0;
      countp[0] = 1; countp[1] = rank; countp[2] = clen;
      stridep[0] = 1; stridep[1] = 1; stridep[2] = 1;
      err = nc_put_vars(ncid_out,Uid,startp,countp,stridep,umat); HANDLE_ERR(err,loc);

      //write right singular vector matrix to outfile
      startp[0] = (size_t)g; startp[1] = 0; startp[2] = 0;
      countp[0] = 1; countp[1] = N_t; countp[2] = rank;
      stridep[0] = 1; stridep[1] = 1; stridep[2] = 1;
      err = nc_put_vars(ncid_out,Vtid,startp,countp,stridep,vtmat); HANDLE_ERR(err,loc);

      //plot the singular values
      strcpy(pname,dname); sprintf(buf,"_g%lu",g+1); strcat(pname,buf);
      err = SIG_PLOT(pname,sig,sigp,(int)rank,drop);

      //deallocating arrays
      free(center); free(umat); free(sig); free(vtmat); free(sigp);

    } //end g loop

  }
  else{
    printf("    -- Start POD on full phase space\n");

    err = POD_CALC(data,N_t,N_g,clen,rank,0,&center,&umat,&sig,&vtmat,&sigp);

    //write centering (column-avg'd data) vector to outfile
    startp[0] = 0; startp[1] = 0; startp[2] = 0;
    countp[0] = clen; countp[1] = 0; countp[2] = 0;
    stridep[0] = 1; stridep[1] = 0; stridep[2] = 0;
    err = nc_put_vars(ncid_out,Cid,startp,countp,stridep,center); HANDLE_ERR(err,loc);

    //write singular value vector to outfile
    startp[0] = 0; startp[1] = 0; startp[2] = 0;
    countp[0] = rank; countp[1] = 0; countp[2] = 0;
    stridep[0] = 1; stridep[1] = 0; stridep[2] = 0;
    err = nc_put_vars(ncid_out,Sid,startp,countp,stridep,sig); HANDLE_ERR(err,loc);

    //write left singular vector matrix to outfile
    startp[0] = 0; startp[1] = 0; startp[2] = 0;
    countp[0] = rank; countp[1] = clen; countp[2] = 0;
    stridep[0] = 1; stridep[1] = 1; stridep[2] = 0;
    err = nc_put_vars(ncid_out,Uid,startp,countp,stridep,umat); HANDLE_ERR(err,loc);

    //write right singular vector matrix to outfile
    startp[0] = 0; startp[1] = 0; startp[2] = 0;
    countp[0] = N_t; countp[1] = rank; countp[2] = 0;
    stridep[0] = 1; stridep[1] = 1; stridep[2] = 0;
    err = nc_put_vars(ncid_out,Vtid,startp,countp,stridep,vtmat); HANDLE_ERR(err,loc);

    //plot the singular values
    strcpy(pname,dname); err = SIG_PLOT(pname,sig,sigp,(int)rank,drop);

    //deallocating arrays
    free(center); free(umat); free(sig); free(vtmat); free(sigp);

  }

  free(data);

  return 0;
}

//================================================================================================================================//
// int SIG_PLOT
//
// SIG_PLOT
//================================================================================================================================//
int SIG_PLOT(const char *pname, const double *sig, const double *sigp, const int rank, const char *drop)
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
int POD_CALC(const double *data, const size_t N_t, const size_t N_g, const size_t clen, const size_t rank, const size_t g, double **center,
  double **umat, double **sig, double **vtmat, double **sigp)
{
  double *temp;
  int i, err;

  //allocating POD arrays
  *center = (double *)malloc(sizeof(double)*clen);
  *umat = (double *)malloc(sizeof(double)*clen*rank);
  *vtmat = (double *)malloc(sizeof(double)*N_t*rank);
  *sig = (double *)malloc(sizeof(double)*rank);

  //forming array of singular value indices
  *sigp = (double *)malloc(sizeof(double)*rank);
  for(i=0;i<(int)rank;i++){
    (*sigp)[i] = (double)(i+1);
  }

  //checking whether a multigroup data-matrix is being used or not
  if(N_g>0){ //if the datamatrix is multigroup, must either isolate one group or form a large multigroup matrix
    //reforming the datamatrix to isolate a single group
    temp = (double *)malloc(sizeof(double)*N_t*clen);
    err = gdat_reform(N_t,N_g,clen,g,data,temp);

    //calculating the SVD of the datamatrix
    SVD_CALC(temp,&N_t,&clen,&*center[0],&*umat[0],&*sig[0],&*vtmat[0]);

    //deallocating arrays
    free(temp);
  }
  else{ //if the datamatrix is not multigroup, can immediately procede with finding the SVD
    //calculating the SVD of the datamatrix
    SVD_CALC(data,&N_t,&clen,&*center[0],&*umat[0],&*sig[0],&*vtmat[0]);
  }

  return 0;
}

//================================================================================================================================//
//
//================================================================================================================================//
int gdat_reform(const size_t n_t, const size_t n_g, const size_t clen, const size_t g, const double *gdat, double *vec)
{

  size_t p1 = 0;
  size_t p2 = 0;
  for(size_t t=0; t<n_t; t++){
    p2 = g*clen+t*n_g*clen;
    for(size_t i=0; i<clen; i++){
      vec[p1] = gdat[i+p2];
      p1 = p1 + 1;
    }
  }

  return 0;

}
