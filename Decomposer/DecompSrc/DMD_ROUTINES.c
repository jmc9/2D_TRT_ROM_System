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
// #include <netcdf.h>
#include <complex.h>
#include <math.h>

/* ----- LOCAL ----- */
#include "Data_Handling.h"

/*================================================================================================================================*/
/* Importing Functions */
/*================================================================================================================================*/
/* ----- FROM NCDF_IO.c ----- */
void Handle_Err(const int Status, const char *Location);
void Get_Var_Double(const int ncid, const char *name, double **var, size_t size);

/* ----- FROM GNUPLOT_ROUTINES.c ----- */
int gnuplot_1d(const char *title, const double *data, const double *crd, const int dim, const char *plttyp, const int logscale,
  const int pt, const char *lc, const char *saveas, const double *xbnd);
void GNUp_Err(const int err);

/* ----- FROM MISC_PROCS.c ----- */
int Make_Dir(const char *dir);

/* ----- FROM LA_ROUTINES.c ----- */
int SVD_Calc(double *dat, const size_t N_t, const size_t clen, double *umat, double *sig, double *vtmat);

int SVD_Calc_cmplx(double complex *dat, const size_t N_t, const size_t clen, double complex *umat, double *sig, double complex *vtmat);

size_t SVD_Rank_Calc(const size_t frank, const double svd_eps, const int opt, const double *s);

int gdat_reform(const size_t N_t, const size_t N_g, const size_t clen, const size_t g, const double *gdat, double *vec);

int Transpose_Double(double *a, double *b, const size_t a_rows, const size_t a_cols);

int Transpose_Double_cmplx(double complex *a, double complex *b, const size_t a_rows, const size_t a_cols);

int MatMul_Double(double *c, double *a, double *b, const size_t rows, const size_t cols, const size_t mdim);

int MatMul_cDouble(double complex *c, double complex *a, double complex *b, const size_t rows, const size_t cols, const size_t mdim);

int EIG_Calc(double *a, size_t len, double complex *wmat, double complex *lambda);

int cGaussElim(const size_t n, double complex *a, double complex *b);

/* ----- LOCAL DEFINITIONS ----- */
#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })

/*================================================================================================================================*/
/**/
/*================================================================================================================================*/
int Lambda_Plot(const char *pname, const double *lambda_r, const double *lambda_i, const int rank, const char *drop)
{
 int i, err;
 char lc[10], plttyp[3], title[100], saveas[100];

 strcpy(lc,"black"); //color of the plotted data
 strcpy(plttyp,"p"); //type of plot
 strcpy(title,pname); strcat(title,"_DMD_eigenvalues\0"); //create plot title
 strcpy(saveas,drop); strcat(saveas,title); strcat(saveas,".png\0"); //create saved name of plot

 //replace all underscores in title as spaces
 i=0;
 while(title[i] != '\0'){
   if(title[i] == '_'){
     title[i] = ' ';
   }
   i++;
 }

 double xbnd[2]={-1,1};
 err = gnuplot_1d(title,lambda_i,lambda_r,rank,plttyp,1,7,lc,saveas,xbnd); GNUp_Err(err);

 return err;
}

/*================================================================================================================================*/
/* Omega_Calc finds:
                        omega = ln( lambda ) / ( Delta t )
   where lambda is assumed imaginary (i.e. lambda = a + ib, a and b may be negative)

   omega is found by transforming ln( a + ib ) = ln( r * exp( i * theta ) ) = ln(r) + i * theta
   where r = sqrt( a^2 + b^2 )
         theta = arctan( |b| / |a| )
*/
/*================================================================================================================================*/
int Omega_Calc(double complex **omega, double complex *lambda, const size_t xr, const double delt)
{
  int err = 0;

  double pi = 4.*atan(1.);
  *omega = (double complex *)malloc(sizeof(double complex)*xr);

  for (size_t i=0; i<xr; i++){
    double clam = cimag(lambda[i]);
    double rlam = creal(lambda[i]);

    //finding trivial cases (real or imaginary component are zero)
    if (clam == 0.){
      if (rlam>0.){
        (*omega)[i] = log(fabs(rlam))/delt + 0.*I;
      }
      else if (rlam<0.){
        (*omega)[i] = (log(fabs(rlam)) + pi*I)/delt;
      }
      continue;
    }
    else if (rlam == 0.){
      if (clam > 0.){
        (*omega)[i] = ( log(fabs(clam)) + (pi / 2.)*I )/delt;
      }
      else{
        (*omega)[i] = ( log(fabs(clam)) + (3.*pi / 2.)*I )/delt;
      }
      continue;
    }

    //both components nonzero, calculate omega as normal
    double theta = atan( fabs(clam) / fabs(rlam) );
    double radius = sqrt( pow(rlam, 2) + pow(clam, 2) );
    if ((rlam > 0.)&&(clam > 0.)){ //first quadrant
      theta = theta;
    }
    else if ((rlam < 0.)&&(clam > 0.)){ //second quadrant
      theta = pi - theta;
    }
    else if ((rlam < 0.)&&(clam < 0.)){ //third quadrant
      theta += pi;
    }
    else if ((rlam > 0.)&&(clam < 0.)){ //fourth quadrant
      theta = 2.*pi - theta;
    }
    (*omega)[i] = ( log(radius) + theta*I )/delt;
  }

  return err;
}

/*================================================================================================================================*/
/**/
/*================================================================================================================================*/
int DMD_Calc(const double *data, const size_t N_t, const size_t N_g, const size_t clen, const size_t rank, const size_t g,
   const double svd_eps, double complex **wmat, double complex **lambda, double complex **wmat_tild, double **umat, size_t *N_modes,
   double **cvec, double complex **wmat_proj)
{
  double *x, *y;
  double *xu, *xs, *xvt, *xv, *a;
  double complex *w, *x2, *x3;
  int err;
  size_t px, py, pa, xr;

  *cvec = (double *)malloc(sizeof(double)*clen);

  //checking whether a multigroup data-matrix is being used or not
  if(N_g>0){ //if the datamatrix is multigroup, must either isolate one group or form a large multigroup matrix
    //reforming the datamatrix to isolate a single group
    a = (double *)malloc(sizeof(double)*N_t*clen);
    err = gdat_reform(N_t,N_g,clen,g,data,a);

    for(size_t i=0; i<clen; i++){
      (*cvec)[i] = 0.;
    }

    x = (double *)malloc(sizeof(double)*(N_t-1)*clen);
    y = (double *)malloc(sizeof(double)*(N_t-1)*clen);

    px = 0; py = 0; pa = 0;
    //copying first column of A to X
    for(size_t i=0; i<clen; i++){
      x[px] = a[pa] - (*cvec)[i];
      px++;
      pa++;
    }

    //copying columns 2 through N_t-1 of A to X and Y
    for(size_t t=1; t<(N_t-1); t++){
      for(size_t i=0; i<clen; i++){
        x[px] = a[pa] - (*cvec)[i];
        y[py] = a[pa] - (*cvec)[i];
        px++;
        py++;
        pa++;
      }
    }

    //copying last column of A to Y
    for(size_t i=0; i<clen; i++){
      y[py] = a[pa] - (*cvec)[i];
      py++;
      pa++;
    }

    //deallocating arrays
    free(a);
  }
  else{ //if the datamatrix is not multigroup, can immediately procede with finding the SVD

    for(size_t i=0; i<clen; i++){
      (*cvec)[i] = 0.;
    }

    x = (double *)malloc(sizeof(double)*(N_t-1)*clen);
    y = (double *)malloc(sizeof(double)*(N_t-1)*clen);

    px = 0; py = 0; pa = 0;
    //copying first column of A to X
    for(size_t i=0; i<clen; i++){
      x[px] = data[pa] - (*cvec)[i];
      px++;
      pa++;
    }

    //copying columns 2 through N_t-1 of A to X and Y
    for(size_t t=1; t<(N_t-1); t++){
      for(size_t i=0; i<clen; i++){
        x[px] = data[pa] - (*cvec)[i];
        y[py] = data[pa] - (*cvec)[i];
        px++;
        py++;
        pa++;
      }
    }

    //copying last column of A to Y
    for(size_t i=0; i<clen; i++){
      y[py] = data[pa] - (*cvec)[i];
      py++;
      pa++;
    }

  }

  //allocating space for SVD X=U*S*V^T
  xu = (double *)malloc(sizeof(double)*clen*rank);     //xu = U
  xvt = (double *)malloc(sizeof(double)*(N_t-1)*rank); //xvt = V^T
  xs = (double *)malloc(sizeof(double)*rank);          //xs = S

  //calculating the SVD of X
  err = SVD_Calc(x,N_t-1,clen,xu,xs,xvt);
  free(x);

  //checking whether to use reduced- or full-rank SVD of x in calculations
  if (svd_eps < 0.){ //negative svd_eps flags the use of full-rank SVD
    xr = rank;
    printf("    -- Using full-rank SVD\n");
  }
  else{ //positive svd_eps flags the use of reduced-rank SVD
    xr = SVD_Rank_Calc(rank,svd_eps,1,xs);
    printf("    -- Truncating SVD of data with rank %ld\n",xr);

    //reallocating V^T with reduced row count
    double *z = malloc(sizeof(double)*(N_t-1)*xr); //temporary array
    size_t p = 0;
    size_t p2 = 0;
    for (size_t j=0; j<(N_t-1); j++){
      for (size_t i=0; i<xr; i++){
        z[p] = xvt[p2]; //copying first xr rows of xvt into z
        p++;
        p2++;
      }
      p2 = p2 + rank - xr;
    }
    free(xvt);

    xvt = (double *)malloc(sizeof(double)*(N_t-1)*xr); //allocating xvt with rank = xr
    p = 0;
    for (size_t j=0; j<(N_t-1); j++){
      for (size_t i=0; i<xr; i++){
        xvt[p] = z[p]; //copying data back into xvt
        p++;
      }
    }
    free(z);

    //reallocating U with reduced column count
    z = (double *)malloc(sizeof(double)*clen*xr); //temporary array
    p = 0;
    for (size_t j=0; j<xr; j++){
      for (size_t i=0; i<clen; i++){
        z[p] = xu[p]; //copying first xr columns of xu into z
        p++;
      }
    }
    free(xu);

    xu = (double *)malloc(sizeof(double)*clen*xr); //allocating xu with rank = xr
    p = 0;
    for (size_t j=0; j<xr; j++){
      for (size_t i=0; i<clen; i++){
        xu[p] = z[p]; //copying data back into xu
        p++;
      }
    }
    free(z);

  }

  //allocating DMD arrays
  *wmat = (double complex *)malloc(sizeof(double complex)*clen*xr);
  *wmat_proj = (double complex *)malloc(sizeof(double complex)*clen*xr);
  *lambda = (double complex *)malloc(sizeof(double complex)*xr);
  *wmat_tild = (double complex *)malloc(sizeof(double complex)*xr*xr);
  *umat = (double *)malloc(sizeof(double)*clen*xr);

  //transposing V^T (finding V = (V^T)^T)
  xv = (double *)malloc(sizeof(double)*(N_t-1)*xr);
  err = Transpose_Double(xvt,xv,xr,(N_t-1));
  free(xvt);

  //multiplying V*S^{-1}
  px = 0;
  for(size_t i=0; i<xr; i++){
    for(size_t j=0; j<(N_t-1); j++){
      xv[px] = xv[px]/xs[i];
      px++;
    }
  }
  free(xs);

  //multiplying Y*(V*S^{-1}) (i.e y*xv), storing in x
  x = (double *)malloc(sizeof(double)*clen*xr);
  err = MatMul_Double(x,y,xv,clen,xr,N_t-1); //finding x = y*xv
  free(y); free(xv);

  //transposing U, using y as temporary matrix
  y = (double *)malloc(sizeof(double)*clen*xr);
  for(size_t i=0; i<clen*xr; i++){ //copying U into y
    y[i] = xu[i];
    (*umat)[i] = xu[i];
  }
  err = Transpose_Double(y,xu,clen,xr); //finding U^T, writing into xu
  free(y);

  //calculating A = U^T * Y * V * S^{-1} (i.e. A = U^T * B where B = Y * V * S^{-1})
  a = (double *)malloc(sizeof(double)*xr*xr);
  err = MatMul_Double(a,xu,x,xr,xr,clen); //finding a = xu*x where xu=U^T and x=Y * V * S^{-1}
  free(xu);

  //calculating the eigendecomposition of a (a*w = l*w)
  err = EIG_Calc(a,xr,*wmat_tild,*lambda);

  x2 = (double complex *)malloc(sizeof(double complex)*clen*xr);
  for (size_t i=0; i<clen*xr; i++){
    x2[i] = x[i] + 0.*I;
  }
  free(x);

  err = MatMul_cDouble(*wmat,x2,*wmat_tild,clen,xr,xr);
  free(x2);

  px = 0;
  for (size_t j=0; j<xr; j++){
    for (size_t i=0; i<clen; i++){
      (*wmat)[px] = (*wmat)[px]/ (*lambda)[j];
      px++;
    }
  }

  size_t ulen = clen*xr;
  double complex *uc = (double complex*)malloc(sizeof(double complex)*ulen);
  for (size_t i=0; i<ulen;){
    uc[i] = (*umat)[i] + 0.*I;
    ++i;
  }
  err = MatMul_cDouble(*wmat_proj,uc,*wmat_tild,clen,xr,xr);
  free(uc);

  *N_modes = xr;

  return err;
}

/*================================================================================================================================*/
/**/
/*================================================================================================================================*/
int Coef_Calc(double complex **b, const double *data, double complex *wmat_tild, double *umat, const size_t N_t, const size_t N_g, const size_t clen,
  const size_t rank, const size_t g, const size_t xr, double complex *lambda, double *cvec, double complex *wmat)
{
  int err;
  int p;

  //a0 holds the fist data-vector
  double complex *a0 = (double complex *)malloc(sizeof(double complex)*clen);

  if(N_g>0){
    //reforming the datamatrix to isolate a single group
    double *a = (double *)malloc(sizeof(double)*N_t*clen);
    err = gdat_reform(N_t,N_g,clen,g,data,a);

    //copying first column of A to a0
    for(size_t i=0; i<clen; i++){
      a0[i] = a[i] + 0.*I;
    }
    free(a);

  }
  else{
    //copying first column of A to a0
    for(size_t i=0; i<clen; i++){
      a0[i] = data[i] - cvec[i] + 0.*I;
    }

  }

  //transposing umat -> umat_T
  size_t ulen = clen*xr;
  double *umat_T= (double *)malloc(sizeof(double)*ulen);
  err = Transpose_Double(umat, umat_T, clen, xr);

  //cumat_T holds umat_T in complex space
  double complex *cumat_T = (double complex *)malloc(sizeof(double complex)*ulen);
  for(size_t i=0; i<ulen; i++){
    cumat_T[i] = umat_T[i] + 0.*I;
  }

  //calculating U^T * a0, storing in b
  *b = (double complex *)malloc(sizeof(double complex)*xr);
  err = MatMul_cDouble(*b,cumat_T,a0,xr,1,clen);
  free(cumat_T);

  //making a copy of wmat_tild since cGaussElim will destroy the matrix it is given
  size_t wlen = xr*xr;
  double complex *w = (double complex *)malloc(sizeof(double complex)*wlen);

  for(size_t i=0; i<wlen; i++){
    w[i] = wmat_tild[i];
  }

  //solving W * x = U^T * a0 for x, storing solution in b
  err = cGaussElim(xr, w, *b);

  free(w);

  return err;
}

/*================================================================================================================================*/
/**/
/*================================================================================================================================*/
int Generate_DMD(const double *data, const int ncid_out, const char *dname, const size_t N_t, const size_t N_g, const size_t clen,
  const size_t rank, const double svd_eps, Data *Decomp, const double delt)
{
  int err;
  int *N_modes;
  char loc[13] = "Generate_DMD";
  char buf[25], pname[25], drop[25];
  double complex *wmat, *lambda, *wmat_tild, *b, *omega, *wmat_proj;
  double *temp, *temp2, *umat, *cvec;
  size_t nc_p1, nc_p2;

  size_t startp[3], countp[3], xr, gcount;
  ptrdiff_t stridep[3];

  //creating directory to put plots
  err = Make_Dir(dname); //setting up directory
  strcpy(drop,dname); strcat(drop,"/"); //creating path to directory

  //checking type of dataset to perform POD on
  if(N_g > 0){ //if N_g>0, then a multigroup dataset has been detected
    printf("    -- groupwise decomposition detected\n");

    N_modes = (int *)malloc(sizeof(int)*N_g);

    gcount = N_g;

    nc_p1=1; nc_p2=2;
    startp[0] = 0; startp[1] = 0; startp[2] = 0;
    countp[0] = 1; countp[1] = 0; countp[2] = 0;
    stridep[0] = 1; stridep[1] = 0; stridep[2] = 0;

  }
  else{
    printf("    -- Start DMD on full phase space\n");

    N_modes = (int *)malloc(sizeof(int)*1);

    gcount = 1;

    nc_p1=0; nc_p2=1;
    startp[0] = 0; startp[1] = 0; startp[2] = 0;
    countp[0] = 0; countp[1] = 0; countp[2] = 0;
    stridep[0] = 0; stridep[1] = 0; stridep[2] = 0;

  }

  for (size_t g=0; g<gcount; g++){
    startp[0] = g;

    //Perform the DMD algorithm, calculate DMD eigenpairs, reduced eigenvectors and left singular values of first orbital data matrix (X)
    err = DMD_Calc(data,N_t,N_g,clen,rank,g,svd_eps,&wmat,&lambda,&wmat_tild,&umat,&xr,&cvec,&wmat_proj);

    //Calculate coefficients of DMD expansion fit to the training data
    err = Coef_Calc(&b, data, wmat_tild, umat, N_t, N_g, clen, rank, g, xr, lambda, cvec, wmat);

    //calculating exponential eigenvalues (omegas)
    err = Omega_Calc(&omega, lambda, xr, delt);

    /*--------------------------------------
         OUTPUT DMD eigenvalues
    --------------------------------------*/
    //seperating real/imaginary components of DMD eigenvalues for output
    temp = (double *)malloc(sizeof(double)*xr);
    temp2 = (double *)malloc(sizeof(double)*xr);
    for (size_t i=0; i<xr; i++){
      temp[i] = creal(lambda[i]);
      temp2[i] = cimag(lambda[i]);
    }
    free(lambda);
    //setting parameters for NCDF file write
    startp[nc_p1] = 0; startp[nc_p2] = 0;
    countp[nc_p1] = xr; countp[nc_p2] = 0;
    stridep[nc_p1] = 1; stridep[nc_p2] = 0;
    err = nc_put_vars(ncid_out,Decomp[0].id,startp,countp,stridep,temp); Handle_Err(err,loc);
    err = nc_put_vars(ncid_out,Decomp[1].id,startp,countp,stridep,temp2); Handle_Err(err,loc);

    /*--------------------------------------
         OUTPUT DMD eigenvalues (exponential form)
    --------------------------------------*/
    //seperating real/imaginary components of exponential eigenvalues (omegas) for output
    for (size_t i=0; i<xr; i++){
      temp[i] = creal(omega[i]);
      temp2[i] = cimag(omega[i]);
    }
    free(omega);
    //using the same NCDF write parameters defined above
    err = nc_put_vars(ncid_out,Decomp[4].id,startp,countp,stridep,temp); Handle_Err(err,loc);
    err = nc_put_vars(ncid_out,Decomp[5].id,startp,countp,stridep,temp2); Handle_Err(err,loc);

    /*--------------------------------------
         OUTPUT DMD expansion coefficients
    --------------------------------------*/
    //seperating real/imaginary components of DMD expansion coefficients for output
    for (size_t i=0; i<xr; i++){
      temp[i] = creal(b[i]);
      temp2[i] = cimag(b[i]);
    }
    free(b);
    //using the same NCDF write parameters defined above
    err = nc_put_vars(ncid_out,Decomp[9].id,startp,countp,stridep,temp); Handle_Err(err,loc);
    err = nc_put_vars(ncid_out,Decomp[10].id,startp,countp,stridep,temp2); Handle_Err(err,loc);

    /*--------------------------------------
         OUTPUT DMD modes
    --------------------------------------*/
    free(temp); free(temp2); //need to re-allocate temp arrays for modes
    temp = (double *)malloc(sizeof(double)*xr*clen);
    temp2 = (double *)malloc(sizeof(double)*xr*clen);
    //seperating real/imaginary components
    for (size_t i=0; i<xr*clen; i++){
      temp[i] = creal(wmat[i]);
      temp2[i] = cimag(wmat[i]);
    }
    free(wmat);
    //setting parameters for NCDF file write
    startp[nc_p1] = 0; startp[nc_p2] = 0;
    countp[nc_p1] = xr; countp[nc_p2] = clen;
    stridep[nc_p1] = 1; stridep[nc_p2] = 1;
    err = nc_put_vars(ncid_out,Decomp[2].id,startp,countp,stridep,temp); Handle_Err(err,loc);
    err = nc_put_vars(ncid_out,Decomp[3].id,startp,countp,stridep,temp2); Handle_Err(err,loc);

    /*--------------------------------------
         OUTPUT DMD projected modes
    --------------------------------------*/
    //seperating real/imaginary components for output
    for (size_t i=0; i<xr*clen; i++){
      temp[i] = creal(wmat_proj[i]);
      temp2[i] = cimag(wmat_proj[i]);
    }
    free(wmat_proj);
    //using the same NCDF write parameters defined above
    err = nc_put_vars(ncid_out,Decomp[12].id,startp,countp,stridep,temp); Handle_Err(err,loc);
    err = nc_put_vars(ncid_out,Decomp[13].id,startp,countp,stridep,temp2); Handle_Err(err,loc);

    /*--------------------------------------
         OUTPUT left singular vector matrix of X
    --------------------------------------*/
    //using the same NCDF write parameters defined above
    err = nc_put_vars(ncid_out,Decomp[8].id,startp,countp,stridep,umat); Handle_Err(err,loc);
    free(umat);

    /*--------------------------------------
         OUTPUT DMD reduced modes
    --------------------------------------*/
    free(temp); free(temp2); //need to re-allocate temp arrays for reduced modes
    temp = (double *)malloc(sizeof(double)*xr*xr);
    temp2 = (double *)malloc(sizeof(double)*xr*xr);
    //seperating real/imaginary components for output
    for (size_t i=0; i<xr*xr; i++){
      temp[i] = creal(wmat_tild[i]);
      temp2[i] = cimag(wmat_tild[i]);
    }
    free(wmat_tild);
    //setting parameters for NCDF file write
    startp[nc_p1] = 0; startp[nc_p2] = 0;
    countp[nc_p1] = xr; countp[nc_p2] = xr;
    stridep[nc_p1] = 1; stridep[nc_p2] = 1;
    err = nc_put_vars(ncid_out,Decomp[6].id,startp,countp,stridep,temp); Handle_Err(err,loc);
    err = nc_put_vars(ncid_out,Decomp[7].id,startp,countp,stridep,temp2); Handle_Err(err,loc);

    //final deallocation of temp arrays
    free(temp); free(temp2);

    /*--------------------------------------
         OUTPUT DMD 'centering' vector
    --------------------------------------*/
    //setting parameters for NCDF file write
    startp[nc_p1] = 0; startp[nc_p2] = 0;
    countp[nc_p1] = clen; countp[nc_p2] = 0;
    stridep[nc_p1] = 1; stridep[nc_p2] = 0;
    err = nc_put_vars(ncid_out,Decomp[11].id,startp,countp,stridep,cvec); Handle_Err(err,loc);
    free(cvec);

    N_modes[g] = (int)xr;

  }

  /*--------------------------------------
       OUTPUT number of calculated DMD modes
  --------------------------------------*/
  memset(pname,0,25); strcpy(pname,"N_modes_"); strcat(pname,dname);
  err = nc_put_att_int(ncid_out,NC_GLOBAL,pname,NC_INT,gcount,N_modes); Handle_Err(err,loc);

  return err;
}
