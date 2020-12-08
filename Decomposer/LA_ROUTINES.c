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
int SVD_Calc(const double *dat, const size_t N_t, const size_t clen, double *center, double *umat, double *sig, double *vtmat)
{
  int p, err;
  double nt, *work, *temp;
  char jobu[2], jobvt[2];
  size_t Lwork;

  temp = (double *)malloc(sizeof(double)*clen*N_t);

  //filling in centering vector with zeros
  for(size_t i=0; i<clen; i++){
    center[i] = 0.;
  }

  //calculating centering vector as column sum of data matrix
  p = 0;
  for(size_t t=0; t<N_t; t++){
    for(size_t i=0; i<clen; i++){
      center[i] = center[i] + dat[p];
      temp[p] = dat[p];
      p = p + 1;
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
      temp[p] = temp[p] - center[i];
      p = p + 1;
    }
  }

  //allocating workspace for dgesvd
  Lwork = max(3*min(clen,N_t)+max(clen,N_t),5*min(clen,N_t));
  Lwork = max((size_t)1,Lwork);
  work = (double *)malloc(sizeof(double)*Lwork);

  //options for dgesvd
  strcpy(jobu,"S"); strcpy(jobvt,"S");

  //finding the SVD of the data matrix (dat) with the LAPACKE package
  err = LAPACKE_dgesvd(LAPACK_COL_MAJOR,*jobu,*jobvt,(int)clen,(int)N_t,temp,(int)clen,sig,umat,(int)clen,vtmat,(int)N_t,work);
  if (err!=0){printf("SVD failed!");}

  //deallocating local arrays
  free(work); free(temp);

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