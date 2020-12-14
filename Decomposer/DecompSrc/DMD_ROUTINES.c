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

size_t SVD_Rank_Calc(const size_t frank, const double svd_eps, const int opt, const double *s);

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
/**/
/*================================================================================================================================*/
int DMD_Calc(const double *data, const size_t N_t, const size_t N_g, const size_t clen, const size_t rank, const size_t g,
   const double svd_eps, double complex **wmat, double complex **lambda, size_t *N_modes)
{
  double *x, *y;
  double *xu, *xs, *xvt, *xv, *a;
  double complex *w, *x2;
  int err;
  size_t px, py, pa, xr;

  //allocating DMD arrays
  *wmat = (double complex *)malloc(sizeof(double complex)*clen*rank);
  *lambda = (double complex *)malloc(sizeof(double complex)*rank);
  w = (double complex *)malloc(sizeof(double complex)*rank*rank);

  //checking whether a multigroup data-matrix is being used or not
  if(N_g>0){ //if the datamatrix is multigroup, must either isolate one group or form a large multigroup matrix
    //reforming the datamatrix to isolate a single group
    a = (double *)malloc(sizeof(double)*N_t*clen);
    err = gdat_reform(N_t,N_g,clen,g,data,a);

    x = (double *)malloc(sizeof(double)*(N_t-1)*clen);
    y = (double *)malloc(sizeof(double)*(N_t-1)*clen);

    px = 0; py = 0; pa = 0;
    //copying first column of A to X
    for(size_t i=0; i<clen; i++){
      x[px] = a[pa];
      px++;
      pa++;
    }

    //copying columns 2 through N_t-1 of A to X and Y
    for(size_t t=1; t<(N_t-1); t++){
      for(size_t i=0; i<clen; i++){
        x[px] = a[pa];
        y[py] = a[pa];
        px++;
        py++;
        pa++;
      }
    }

    //copying last column of A to Y
    for(size_t i=0; i<clen; i++){
      y[py] = a[pa];
      py++;
      pa++;
    }

    //deallocating arrays
    free(a);
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
  else{ //negative svd_eps flags the use of reduced-rank SVD
    xr = SVD_Rank_Calc(rank,svd_eps,0,xs);
    printf("    -- Truncating SVD of data with rank %ld\n",xr);
  }

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
  }
  err = Transpose_Double(y,xu,clen,xr); //finding U^T, writing into xu
  free(y);

  //calculating A = U^T * Y * V * S^{-1} (i.e. A = U^T * B where B = Y * V * S^{-1})
  a = (double *)malloc(sizeof(double)*xr*xr);
  err = MatMul_Double(a,xu,x,xr,xr,clen); //finding a = xu*x where xu=U^T and x=Y * V * S^{-1}
  free(xu);

  //calculating the eigendecomposition of a (a*w = l*w)
  err = EIG_Calc(a,xr,w,*lambda);

  x2 = (double complex *)malloc(sizeof(double complex)*clen*xr);
  for (size_t i=0; i<clen*xr; i++){
    x2[i] = x[i] + 0.*I;
  }
  free(x);

  err = MatMul_cDouble(*wmat,x2,w,clen,xr,xr);
  free(x2); free(w);

  px = 0;
  for (size_t j=0; j<xr; j++){
    for (size_t i=0; i<clen; i++){
      (*wmat)[px] = (*wmat)[px]/ (*lambda)[j];
      px++;
    }
  }

  *N_modes = xr;

  return err;
}

/*================================================================================================================================*/
/**/
/*================================================================================================================================*/
int Generate_DMD(const double *data, const int ncid_out, const char *dname, const size_t N_t, const size_t N_g, const size_t clen,
  const size_t rank, const double svd_eps, Data *Decomp)
{
  int err;
  int *N_modes;
  char loc[13] = "Generate_DMD";
  char buf[25], pname[25], drop[25];
  double complex *wmat, *lambda;
  double *temp, *temp2;

  size_t startp[3], countp[3], xr;
  ptrdiff_t stridep[3];

  //creating directory to put plots
  err = Make_Dir(dname); //setting up directory
  strcpy(drop,dname); strcat(drop,"/"); //creating path to directory

  //checking type of dataset to perform POD on
  if(N_g > 0){ //if N_g>0, then a multigroup dataset has been detected
    printf("    -- groupwise decomposition detected\n");

    N_modes = (int *)malloc(sizeof(int)*N_g);

    //loop over energy groups
    for(size_t g=0; g<N_g; g++){
      printf("    -- Start DMD on group %lu\n",g+1);

      //find the POD modes and singular values of a given groupwise datamatrix
      err = DMD_Calc(data,N_t,N_g,clen,rank,g,svd_eps,&wmat,&lambda,&xr);

      temp = (double *)malloc(sizeof(double)*xr);
      for (size_t i=0; i<xr; i++){
        temp[i] = creal(lambda[i]);
      }

      //write singular value vector to outfile
      startp[0] = (size_t)g; startp[1] = 0; startp[2] = 0;
      countp[0] = 1; countp[1] = xr; countp[2] = 0;
      stridep[0] = 1; stridep[1] = 1; stridep[2] = 0;
      err = nc_put_vars(ncid_out,Decomp[0].id,startp,countp,stridep,temp); Handle_Err(err,loc);

      temp2 = (double *)malloc(sizeof(double)*xr);
      for (size_t i=0; i<xr; i++){
        temp2[i] = cimag(lambda[i]);
      }

      //write DMD eigenvalue (imaginary component) vector to outfile
      err = nc_put_vars(ncid_out,Decomp[1].id,startp,countp,stridep,temp2); Handle_Err(err,loc);

      //plot the singular values
      strcpy(pname,dname); sprintf(buf,"_g%lu",g+1); strcat(pname,buf);
      err = Lambda_Plot(pname,temp,temp2,(int)xr,drop);

      free(temp);
      temp = (double *)malloc(sizeof(double)*xr*clen);
      for (size_t i=0; i<xr*clen; i++){
        temp[i] = creal(wmat[i]);
      }

      //write left singular vector matrix to outfile
      startp[0] = (size_t)g; startp[1] = 0; startp[2] = 0;
      countp[0] = 1; countp[1] = xr; countp[2] = clen;
      stridep[0] = 1; stridep[1] = 1; stridep[2] = 1;
      err = nc_put_vars(ncid_out,Decomp[2].id,startp,countp,stridep,temp); Handle_Err(err,loc);

      for (size_t i=0; i<xr*clen; i++){
        temp[i] = cimag(wmat[i]);
      }

      //write DMD mode matrix (imaginary component) to outfile
      err = nc_put_vars(ncid_out,Decomp[3].id,startp,countp,stridep,temp); Handle_Err(err,loc);

      //deallocating arrays
      free(wmat); free(lambda); free(temp); free(temp2);

      N_modes[g] = (int)xr;

    } //end g loop

    memset(pname,0,25); strcpy(pname,"N_modes_"); strcat(pname,dname);
    err = nc_put_att_int(ncid_out,NC_GLOBAL,pname,NC_INT,N_g,N_modes); Handle_Err(err,loc);
    free(N_modes);

  }
  else{
    printf("    -- Start DMD on full phase space\n");

    err = DMD_Calc(data,N_t,N_g,clen,rank,0,svd_eps,&wmat,&lambda,&xr);

    temp = (double *)malloc(sizeof(double)*xr);
    for (size_t i=0; i<xr; i++){
      temp[i] = creal(lambda[i]);
    }

    //write DMD eigenvalue (real component) vector to outfile
    startp[0] = 0; startp[1] = 0; startp[2] = 0;
    countp[0] = xr; countp[1] = 0; countp[2] = 0;
    stridep[0] = 1; stridep[1] = 0; stridep[2] = 0;
    err = nc_put_vars(ncid_out,Decomp[0].id,startp,countp,stridep,temp); Handle_Err(err,loc);

    temp2 = (double *)malloc(sizeof(double)*xr);
    for (size_t i=0; i<xr; i++){
      temp2[i] = cimag(lambda[i]);
    }

    //write DMD eigenvalue (imaginary component) vector to outfile
    err = nc_put_vars(ncid_out,Decomp[1].id,startp,countp,stridep,temp2); Handle_Err(err,loc);

    //plot the singular values
    strcpy(pname,dname); err = Lambda_Plot(pname,temp,temp2,(int)xr,drop);

    free(temp);
    temp = (double *)malloc(sizeof(double)*xr*clen);
    for (size_t i=0; i<xr*clen; i++){
      temp[i] = creal(wmat[i]);
    }

    //write DMD mode matrix (real component) to outfile
    startp[0] = 0; startp[1] = 0; startp[2] = 0;
    countp[0] = xr; countp[1] = clen; countp[2] = 0;
    stridep[0] = 1; stridep[1] = 1; stridep[2] = 0;
    err = nc_put_vars(ncid_out,Decomp[2].id,startp,countp,stridep,temp); Handle_Err(err,loc);

    for (size_t i=0; i<xr*clen; i++){
      temp[i] = cimag(wmat[i]);
    }

    //write DMD mode matrix (imaginary component) to outfile
    err = nc_put_vars(ncid_out,Decomp[3].id,startp,countp,stridep,temp); Handle_Err(err,loc);

    //deallocating arrays
    free(wmat); free(lambda); free(temp); free(temp2);

    N_modes = (int *)malloc(sizeof(int)*1);
    N_modes[0] = (int)xr;
    memset(pname,0,25); strcpy(pname,"N_modes_"); strcat(pname,dname);
    err = nc_put_att_int(ncid_out,NC_GLOBAL,pname,NC_INT,1,N_modes); Handle_Err(err,loc);

  }

  return err;
}
// int Generate_DMD(const double *data, const int ncid_out, const char *dname, const size_t N_t, const size_t N_g, const size_t clen,
//   const size_t rank, const double svd_eps, const int *DCMP_IDs)
// {
//   int err;
//   int *N_modes;
//   char loc[13] = "Generate_DMD";
//   char buf[25], pname[25], drop[25];
//   double complex *wmat, *lambda;
//   double *temp, *temp2;
//
//   size_t startp[3], countp[3], xr;
//   ptrdiff_t stridep[3];
//
//   //creating directory to put plots
//   err = Make_Dir(dname); //setting up directory
//   strcpy(drop,dname); strcat(drop,"/"); //creating path to directory
//
//   //checking type of dataset to perform POD on
//   if(N_g > 0){ //if N_g>0, then a multigroup dataset has been detected
//     printf("    -- groupwise decomposition detected\n");
//
//     N_modes = (int *)malloc(sizeof(int)*N_g);
//
//     //loop over energy groups
//     for(size_t g=0; g<N_g; g++){
//       printf("    -- Start DMD on group %lu\n",g+1);
//
//       //find the POD modes and singular values of a given groupwise datamatrix
//       err = DMD_Calc(data,N_t,N_g,clen,rank,g,svd_eps,&wmat,&lambda,&xr);
//
//       temp = (double *)malloc(sizeof(double)*xr);
//       for (size_t i=0; i<xr; i++){
//         temp[i] = creal(lambda[i]);
//       }
//
//       //write singular value vector to outfile
//       startp[0] = (size_t)g; startp[1] = 0; startp[2] = 0;
//       countp[0] = 1; countp[1] = xr; countp[2] = 0;
//       stridep[0] = 1; stridep[1] = 1; stridep[2] = 0;
//       err = nc_put_vars(ncid_out,DCMP_IDs[0],startp,countp,stridep,temp); Handle_Err(err,loc);
//
//       temp2 = (double *)malloc(sizeof(double)*xr);
//       for (size_t i=0; i<xr; i++){
//         temp2[i] = cimag(lambda[i]);
//       }
//
//       //write DMD eigenvalue (imaginary component) vector to outfile
//       err = nc_put_vars(ncid_out,DCMP_IDs[1],startp,countp,stridep,temp2); Handle_Err(err,loc);
//
//       //plot the singular values
//       strcpy(pname,dname); sprintf(buf,"_g%lu",g+1); strcat(pname,buf);
//       err = Lambda_Plot(pname,temp,temp2,(int)xr,drop);
//
//       free(temp);
//       temp = (double *)malloc(sizeof(double)*xr*clen);
//       for (size_t i=0; i<xr*clen; i++){
//         temp[i] = creal(wmat[i]);
//       }
//
//       //write left singular vector matrix to outfile
//       startp[0] = (size_t)g; startp[1] = 0; startp[2] = 0;
//       countp[0] = 1; countp[1] = xr; countp[2] = clen;
//       stridep[0] = 1; stridep[1] = 1; stridep[2] = 1;
//       err = nc_put_vars(ncid_out,DCMP_IDs[2],startp,countp,stridep,temp); Handle_Err(err,loc);
//
//       for (size_t i=0; i<xr*clen; i++){
//         temp[i] = cimag(wmat[i]);
//       }
//
//       //write DMD mode matrix (imaginary component) to outfile
//       err = nc_put_vars(ncid_out,DCMP_IDs[3],startp,countp,stridep,temp); Handle_Err(err,loc);
//
//       //deallocating arrays
//       free(wmat); free(lambda); free(temp); free(temp2);
//
//       N_modes[g] = (int)xr;
//
//     } //end g loop
//
//     memset(pname,0,25); strcpy(pname,"N_modes_"); strcat(pname,dname);
//     err = nc_put_att_int(ncid_out,NC_GLOBAL,pname,NC_INT,N_g,N_modes); Handle_Err(err,loc);
//     free(N_modes);
//
//   }
//   else{
//     printf("    -- Start DMD on full phase space\n");
//
//     err = DMD_Calc(data,N_t,N_g,clen,rank,0,svd_eps,&wmat,&lambda,&xr);
//
//     temp = (double *)malloc(sizeof(double)*xr);
//     for (size_t i=0; i<xr; i++){
//       temp[i] = creal(lambda[i]);
//     }
//
//     //write DMD eigenvalue (real component) vector to outfile
//     startp[0] = 0; startp[1] = 0; startp[2] = 0;
//     countp[0] = xr; countp[1] = 0; countp[2] = 0;
//     stridep[0] = 1; stridep[1] = 0; stridep[2] = 0;
//     err = nc_put_vars(ncid_out,DCMP_IDs[0],startp,countp,stridep,temp); Handle_Err(err,loc);
//
//     temp2 = (double *)malloc(sizeof(double)*xr);
//     for (size_t i=0; i<xr; i++){
//       temp2[i] = cimag(lambda[i]);
//     }
//
//     //write DMD eigenvalue (imaginary component) vector to outfile
//     err = nc_put_vars(ncid_out,DCMP_IDs[1],startp,countp,stridep,temp2); Handle_Err(err,loc);
//
//     //plot the singular values
//     strcpy(pname,dname); err = Lambda_Plot(pname,temp,temp2,(int)xr,drop);
//
//     free(temp);
//     temp = (double *)malloc(sizeof(double)*xr*clen);
//     for (size_t i=0; i<xr*clen; i++){
//       temp[i] = creal(wmat[i]);
//     }
//
//     //write DMD mode matrix (real component) to outfile
//     startp[0] = 0; startp[1] = 0; startp[2] = 0;
//     countp[0] = xr; countp[1] = clen; countp[2] = 0;
//     stridep[0] = 1; stridep[1] = 1; stridep[2] = 0;
//     err = nc_put_vars(ncid_out,DCMP_IDs[2],startp,countp,stridep,temp); Handle_Err(err,loc);
//
//     for (size_t i=0; i<xr*clen; i++){
//       temp[i] = cimag(wmat[i]);
//     }
//
//     //write DMD mode matrix (imaginary component) to outfile
//     err = nc_put_vars(ncid_out,DCMP_IDs[3],startp,countp,stridep,temp); Handle_Err(err,loc);
//
//     //deallocating arrays
//     free(wmat); free(lambda); free(temp); free(temp2);
//
//     N_modes = (int *)malloc(sizeof(int)*1);
//     N_modes[0] = (int)xr;
//     memset(pname,0,25); strcpy(pname,"N_modes_"); strcat(pname,dname);
//     err = nc_put_att_int(ncid_out,NC_GLOBAL,pname,NC_INT,1,N_modes); Handle_Err(err,loc);
//
//   }
//
//   return err;
// }
