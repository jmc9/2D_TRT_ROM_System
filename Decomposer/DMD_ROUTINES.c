/*================================================================================================================================*/
/*
  DMD_ROUTINES.c
  Author: Joseph M. Coale
  jmcoale@ncsu.edu - josephcoale912@gmail.com
*/
/*================================================================================================================================*/
/* Importing Packages */
/*================================================================================================================================*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <netcdf.h>
#include <complex.h>

/*================================================================================================================================*/
/* Importing Functions */
/*================================================================================================================================*/
/* ----- FROM NCDF_IO.c ----- */
void Handle_Err(const int Status, const char *Location);
void Get_Var_Double(const int ncid, const char *name, double **var, size_t size);

/* ----- FROM GNUPLOT_ROUTINES.c ----- */
int gnuplot_1d(const char *title, const double *data, const double *crd, const int dim, const char *plttyp, const int logscale,
  const int pt, const char *lc, const char *saveas);
void GNUp_Err(const int err);

/* ----- FROM MISC_PROCS.c ----- */
int Make_Dir(const char *dir);

/* ----- FROM LA_ROUTINES.c ----- */
int SVD_Calc(double *dat, const size_t N_t, const size_t clen, double *umat, double *sig, double *vtmat);

int gdat_reform(const size_t N_t, const size_t N_g, const size_t clen, const size_t g, const double *gdat, double *vec);

int Transpose_Double(double *a, double *b, const size_t a_rows, const size_t a_cols);

int MatMul_Double(double *c, double *a, double *b, const size_t rows, const size_t cols, const size_t mdim);

int MatMul_cDouble(double complex *c, double complex *a, double complex *b, const size_t rows, const size_t cols, const size_t mdim);

int EIG_Calc(double *a, size_t len, double complex *wmat, double complex *lambda);

/* ----- LOCAL DEFINITIONS ----- */
#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })

/*================================================================================================================================*/
/**/
/*================================================================================================================================*/
int DMD_Calc(const double *data, const size_t N_t, const size_t N_g, const size_t clen, const size_t rank, const size_t g,
  double complex **wmat, double complex **lambda, double **lamp)
{
  double *x, *y;
  double *xu, *xs, *xvt, *xv, *a;
  double complex *w, *x2;
  int err;
  size_t px, py, pa;

  //allocating DMD arrays
  *wmat = (double complex *)malloc(sizeof(double complex)*clen*rank);
  *lambda = (double complex *)malloc(sizeof(double complex)*rank);
  w = (double complex *)malloc(sizeof(double complex)*rank*rank);

  //forming array of singular value indices
  *lamp = (double *)malloc(sizeof(double)*rank);
  for(size_t i=0; i<rank; i++){
    (*lamp)[i] = (double)(i+1);
  }

  //checking whether a multigroup data-matrix is being used or not
  if(N_g>0){ //if the datamatrix is multigroup, must either isolate one group or form a large multigroup matrix
    // //reforming the datamatrix to isolate a single group
    // temp = (double *)malloc(sizeof(double)*N_t*clen);
    // err = gdat_reform(N_t,N_g,clen,g,data,temp);
    //
    // //calculating the SVD of the datamatrix
    // err = SVD_Calc(temp,N_t,clen,*center,*umat,*sig,*vtmat);
    //
    // //deallocating arrays
    // free(temp);
  }
  else{ //if the datamatrix is not multigroup, can immediately procede with finding the SVD

    x = (double *)malloc(sizeof(double)*(N_t-1)*clen);
    y = (double *)malloc(sizeof(double)*(N_t-1)*clen);

    px = 0; py = 0; pa = 0;
    //copying first column of A to X
    for(size_t i=0; i<clen; i++){
      x[px] = data[pa];
      px++;
      pa++;
    }

    //copying columns 2 through N_t-1 of A to X and Y
    for(size_t t=1; t<(N_t-1); t++){
      for(size_t i=0; i<clen; i++){
        x[px] = data[pa];
        y[py] = data[pa];
        px++;
        py++;
        pa++;
      }
    }

    //copying last column of A to Y
    for(size_t i=0; i<clen; i++){
      y[py] = data[pa];
      py++;
      pa++;
    }

    //allocating space for SVD X=U*S*V^T
    xu = (double *)malloc(sizeof(double)*clen*rank);     //xu = U
    xvt = (double *)malloc(sizeof(double)*(N_t-1)*rank); //xvt = V^T
    xs = (double *)malloc(sizeof(double)*rank);          //xs = S

    //calculating the SVD of X
    err = SVD_Calc(x,N_t-1,clen,xu,xs,xvt);
    free(x);

    //transposing V^T (finding V = (V^T)^T)
    xv = (double *)malloc(sizeof(double)*(N_t-1)*rank);
    err = Transpose_Double(xvt,xv,rank,(N_t-1));
    free(xvt);

    //multiplying V*S^{-1}
    px = 0;
    for(size_t i=0; i<rank; i++){
      for(size_t j=0; j<(N_t-1); j++){
        xv[px] = xv[px]/xs[i];
      }
    }
    free(xs);

    //multiplying Y*(V*S^{-1}) (i.e y*xv), storing in x
    x = (double *)malloc(sizeof(double)*clen*rank);
    err = MatMul_Double(x,y,xv,clen,rank,N_t-1); //finding x = y*xv
    free(y); free(xv);

    //transposing U, using y as temporary matrix
    y = (double *)malloc(sizeof(double)*clen*rank);
    for(size_t i=0; i<clen*rank; i++){ //copying U into y
      y[i] = xu[i];
    }
    err = Transpose_Double(y,xu,clen,rank); //finding U^T, writing into xu
    free(y);

    //calculating A = U^T * Y * V * S^{-1} (i.e. A = U^T * B where B = Y * V * S^{-1})
    a = (double *)malloc(sizeof(double)*rank*rank);
    err = MatMul_Double(a,xu,x,rank,rank,clen); //finding a = xu*x where xu=U^T and x=Y * V * S^{-1}
    free(xu);

    //calculating the eigendecomposition of a (a*w = l*w)
    err = EIG_Calc(a,rank,w,*lambda);

    x2 = (double complex *)malloc(sizeof(double complex)*clen*rank);
    for (size_t i=0; i<clen*rank; i++){
      x2[i] = x[i] + 0.*I;
    }
    free(x);

    err = MatMul_cDouble(*wmat,x2,w,clen,rank,rank);
    free(x2); free(w);

    px = 0;
    for (size_t j=0; j<rank; j++){
      for (size_t i=0; i<clen; i++){
        (*wmat)[px] = (*wmat)[px]/ (*lambda)[j];
        px++;
      }
    }

  }

  return err;
}

/*================================================================================================================================*/
/**/
/*================================================================================================================================*/
int Generate_DMD(const double *data, const int ncid_out, const char *dname, const size_t N_t, const size_t N_g, const size_t clen,
  const size_t rank, const int *DCMP_IDs)
{
  int err;
  char loc[13] = "Generate_DMD";
  char buf[25], pname[25], drop[25];
  double complex *wmat, *lambda;
  double *lamp, *temp;

  size_t startp[3], countp[3];
  ptrdiff_t stridep[3];

  //creating directory to put plots
  err = Make_Dir(dname); //setting up directory
  strcpy(drop,dname); strcat(drop,"/"); //creating path to directory

  //checking type of dataset to perform POD on
  if(N_g > 0){ //if N_g>0, then a multigroup dataset has been detected
    printf("    -- groupwise decomposition detected\n");

    //loop over energy groups
    // for(size_t g=0; g<N_g; g++){
    //   printf("    -- Start POD on group %lu\n",g+1);
    //
    //   //find the POD modes and singular values of a given groupwise datamatrix
    //   err = POD_Calc(data,N_t,N_g,clen,rank,g,&center,&umat,&sig,&vtmat,&sigp);
    //
    //   //write centering (column-avg'd data) vector to outfile
    //   startp[0] = (size_t)g; startp[1] = 0; startp[2] = 0;
    //   countp[0] = 1; countp[1] = clen; countp[2] = 0;
    //   stridep[0] = 1; stridep[1] = 1; stridep[2] = 0;
    //   err = nc_put_vars(ncid_out,DCMP_IDs[0],startp,countp,stridep,center); Handle_Err(err,loc);
    //
    //   //write singular value vector to outfile
    //   startp[0] = (size_t)g; startp[1] = 0; startp[2] = 0;
    //   countp[0] = 1; countp[1] = rank; countp[2] = 0;
    //   stridep[0] = 1; stridep[1] = 1; stridep[2] = 0;
    //   err = nc_put_vars(ncid_out,DCMP_IDs[1],startp,countp,stridep,sig); Handle_Err(err,loc);
    //
    //   //write left singular vector matrix to outfile
    //   startp[0] = (size_t)g; startp[1] = 0; startp[2] = 0;
    //   countp[0] = 1; countp[1] = rank; countp[2] = clen;
    //   stridep[0] = 1; stridep[1] = 1; stridep[2] = 1;
    //   err = nc_put_vars(ncid_out,DCMP_IDs[2],startp,countp,stridep,umat); Handle_Err(err,loc);
    //
    //   //write right singular vector matrix to outfile
    //   startp[0] = (size_t)g; startp[1] = 0; startp[2] = 0;
    //   countp[0] = 1; countp[1] = N_t; countp[2] = rank;
    //   stridep[0] = 1; stridep[1] = 1; stridep[2] = 1;
    //   err = nc_put_vars(ncid_out,DCMP_IDs[3],startp,countp,stridep,vtmat); Handle_Err(err,loc);
    //
    //   //plot the singular values
    //   strcpy(pname,dname); sprintf(buf,"_g%lu",g+1); strcat(pname,buf);
    //   err = Sig_Plot(pname,sig,sigp,(int)rank,drop);
    //
    //   //deallocating arrays
    //   free(center); free(umat); free(sig); free(vtmat); free(sigp);
    //
    // } //end g loop

  }
  else{
    printf("    -- Start DMD on full phase space\n");

    err = DMD_Calc(data,N_t,N_g,clen,rank,0,&wmat,&lambda,&lamp);

    temp = (double *)malloc(sizeof(double)*rank);
    for (size_t i=0; i<rank; i++){
      temp[i] = creal(lambda[i]);
    }

    //write DMD eigenvalue (real component) vector to outfile
    startp[0] = 0; startp[1] = 0; startp[2] = 0;
    countp[0] = rank; countp[1] = 0; countp[2] = 0;
    stridep[0] = 1; stridep[1] = 0; stridep[2] = 0;
    err = nc_put_vars(ncid_out,DCMP_IDs[0],startp,countp,stridep,temp); Handle_Err(err,loc);

    for (size_t i=0; i<rank; i++){
      temp[i] = cimag(lambda[i]);
    }

    //write DMD eigenvalue (imaginary component) vector to outfile
    err = nc_put_vars(ncid_out,DCMP_IDs[1],startp,countp,stridep,temp); Handle_Err(err,loc);

    free(temp);
    temp = (double *)malloc(sizeof(double)*rank*clen);
    for (size_t i=0; i<rank*clen; i++){
      temp[i] = creal(wmat[i]);
    }

    //write DMD mode matrix (real component) to outfile
    startp[0] = 0; startp[1] = 0; startp[2] = 0;
    countp[0] = rank; countp[1] = clen; countp[2] = 0;
    stridep[0] = 1; stridep[1] = 1; stridep[2] = 0;
    err = nc_put_vars(ncid_out,DCMP_IDs[2],startp,countp,stridep,temp); Handle_Err(err,loc);

    for (size_t i=0; i<rank*clen; i++){
      temp[i] = cimag(wmat[i]);
    }

    //write DMD mode matrix (imaginary component) to outfile
    err = nc_put_vars(ncid_out,DCMP_IDs[3],startp,countp,stridep,temp); Handle_Err(err,loc);

    //plot the singular values
    // strcpy(pname,dname); err = Sig_Plot(pname,sig,sigp,(int)rank,drop);

    //deallocating arrays
    free(wmat); free(lambda); free(lamp); free(temp);

  }

  return err;
}
