/*================================================================================================================================*/
/*
  LA_ROUTINES.c
  Author: Joseph M. Coale
  jmcoale@ncsu.edu - josephcoale912@gmail.com
*/
/*================================================================================================================================*/
/* Importing Packages */
/*================================================================================================================================*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <lapacke.h>
#include <cblas.h>
#include <complex.h>

/*================================================================================================================================*/
/* Importing Functions */
/*================================================================================================================================*/
/* ----- LOCAL DEFINITIONS ----- */
#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })

#define max(a,b) \
  ({ __typeof__ (a) _a = (a); \
      __typeof__ (b) _b = (b); \
    _a > _b ? _a : _b; })

/*================================================================================================================================*/
/**/
/*================================================================================================================================*/
int SVD_Calc(double *dat, const size_t N_t, const size_t clen, double *umat, double *sig, double *vtmat)
{
  int err;
  double *work;
  char jobu[2], jobvt[2];
  size_t Lwork;

  //allocating workspace for dgesvd
  Lwork = max(3*min(clen,N_t)+max(clen,N_t),5*min(clen,N_t));
  Lwork = max((size_t)1,Lwork);
  work = (double *)malloc(sizeof(double)*Lwork);

  //options for dgesvd
  strcpy(jobu,"S"); strcpy(jobvt,"S");

  //finding the SVD of the data matrix (dat) with the LAPACKE package
  err = LAPACKE_dgesvd(LAPACK_COL_MAJOR,*jobu,*jobvt,(int)clen,(int)N_t,dat,(int)clen,sig,umat,(int)clen,vtmat,(int)N_t,work);
  if (err!=0){printf("SVD failed!");}

  //deallocating local arrays
  free(work);

  return err;
}

/*================================================================================================================================*/
/**/
/*================================================================================================================================*/
int SVD_Calc_Cntr(double *dat, const size_t N_t, const size_t clen, double *center, double *umat, double *sig, double *vtmat)
{
  int p, err;
  double nt;

  // temp = (double *)malloc(sizeof(double)*clen*N_t);

  //filling in centering vector with zeros
  for(size_t i=0; i<clen; i++){
    center[i] = 0.;
  }

  //calculating centering vector as column sum of data matrix
  p = 0;
  for(size_t t=0; t<N_t; t++){
    for(size_t i=0; i<clen; i++){
      center[i] = center[i] + dat[p];
      // temp[p] = dat[p];
      p++;
    }
  }

  //averaging centering vector
  nt = (double)N_t;
  for(size_t i=0; i<clen; i++){
    center[i] = center[i]/nt;
  }

  //centering the data matrix
  p = 0;
  for(size_t t=0; t<N_t; t++){
    for(size_t i=0; i<clen; i++){
      dat[p] = dat[p] - center[i];
      p++;
    }
  }

  err = SVD_Calc(dat,N_t,clen,umat,sig,vtmat);

  return err;
}

/*================================================================================================================================*/
/**/
/*================================================================================================================================*/
int EIG_Calc(double *a, size_t len, double complex *wmat, double complex *lambda)
{
  int err;
  double *wmat_l, *wmat_r, *wmat_i, *lambda_r, *lambda_i;
  char jobvl[2], jobvr[2];

  wmat_r = (double *)malloc(sizeof(double)*len*len);
  wmat_i = (double *)malloc(sizeof(double)*len*len);
  lambda_r = (double *)malloc(sizeof(double)*len);
  lambda_i = (double *)malloc(sizeof(double)*len);

  //options for dgeev
  strcpy(jobvl,"N"); strcpy(jobvr,"V");

  err = LAPACKE_dgeev(LAPACK_COL_MAJOR,*jobvl,*jobvr,(int)len,a,(int)len,lambda_r,lambda_i,wmat_l,1,wmat_r,(int)len);

  size_t p = 0;
  for (size_t i=1; i<len; i++){
    if (lambda_r[i] == lambda_r[i-1]){
      p = i*len;
      for (size_t j=0; j<len; j++){
        wmat_i[p] = wmat_r[p];
        wmat_i[p-len] = -wmat_r[p];
        wmat_r[p] = wmat_r[p-len];
        p++;
      }
    }
  }

  p = 0;
  for (size_t i=0; i<len; i++){

    lambda[i] = lambda_r[i] + lambda_i[i]*I;

    for (size_t j=0; j<len; j++){

      wmat[p] = wmat_r[p] + wmat_i[p]*I;
      p++;

    }
  }

  return err;
}

/*================================================================================================================================*/
/**/
/*================================================================================================================================*/
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

/*================================================================================================================================*/
/**/
/*================================================================================================================================*/
int Transpose_Double(double *a, double *b, const size_t a_rows, const size_t a_cols)
{
  int err = 0;
  size_t pa = 0;
  size_t pb = 0;
  for(size_t j=0; j<a_rows; j++){
    pa = j;
    for(size_t i=0; i<a_cols; i++){
      b[pb] = a[pa];
      pb++;
      pa = pa + a_rows;
    }
  }

  return err;
}

/*================================================================================================================================*/
/* C = A*B */
/*================================================================================================================================*/
int MatMul_Double(double *c, double *a, double *b, const size_t rows, const size_t cols, const size_t mdim)
{
  int err = 0;

  size_t pa = 0;
  size_t pb = 0;
  size_t pc = 0;
  for(size_t j=0; j<cols; j++){
    for(size_t i=0; i<rows; i++){

      c[pc] = 0.;
      pa = i;
      pb = j*mdim;
      for(size_t k=0; k<mdim; k++){
        c[pc] = c[pc] + a[pa]*b[pb];
        pa = pa + rows;
        pb++;
      }
      pc++;

    }
  }

  return err;
}

/*================================================================================================================================*/
/* C = A*B */
/*================================================================================================================================*/
int MatMul_cDouble(double complex *c, double complex *a, double complex *b, const size_t rows, const size_t cols, const size_t mdim)
{
  int err = 0;

  size_t pa = 0;
  size_t pb = 0;
  size_t pc = 0;
  for(size_t j=0; j<cols; j++){
    for(size_t i=0; i<rows; i++){

      c[pc] = 0.;
      pa = i;
      pb = j*mdim;
      for(size_t k=0; k<mdim; k++){
        c[pc] = c[pc] + a[pa]*b[pb];
        pa = pa + rows;
        pb++;
      }
      pc++;

    }
  }

  return err;
}
