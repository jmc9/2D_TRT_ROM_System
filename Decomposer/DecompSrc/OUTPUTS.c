/*================================================================================================================================*/
/*
  OUTPUTS.c
  Author: Joseph M. Coale
  jmcoale@ncsu.edu - josephcoale912@gmail.com
*/
/*================================================================================================================================*/
/* Importing Packages */
/*================================================================================================================================*/
/* ----- EXTERNAL ----- */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/* ----- LOCAL ----- */
#include "Data_Handling.h"

/*================================================================================================================================*/
/* Importing Functions */
/*================================================================================================================================*/
/* ----- FROM NCDF_IO.c ----- */
void Handle_Err(const int Status, const char *Location);

void Get_Var_Double(const int ncid, const char *name, double **var, const size_t size);

int Write_Spec(const int ncid, Spec *spec);

/* ----- FROM CPOD_ROUTINES.c ----- */
int Generate_POD(const double *data, const int ncid_out, const char *dname, const size_t N_t, const size_t N_g, const size_t clen,
  const size_t rank, Data *Decomp);

int Generate_DMD(const double *data, const int ncid_out, const char *dname, const size_t N_t, const size_t N_g, const size_t clen,
  const size_t rank, const double svd_eps, Data *Decomp, const double delt);

/* ----- FROM MISC_PROCS.c ----- */
void Sort_Uniq_sizet(size_t *list, const size_t len, size_t *ulen);
void Sort_Uniq_int(int *list, const size_t len, size_t *ulen);

int Delimit(const char *line, const char del, char **parts, const int nparts);

/* ----- LOCAL DEFINITIONS ----- */
#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })

/*================================================================================================================================*/
/**/
/*================================================================================================================================*/
int Def_Dims(const int ncid, Data *Dcmp_data, const size_t N_data, ncdim *dims, const size_t N_dims, ncdim **clen,
  ncdim **rank, const int gsum, const int dcmp_type, int ***dcdims)
{
  int err;
  char loc[9] = "Def_Dims";

  size_t *clen_ = malloc(sizeof(size_t)*N_data); //array holding clen of each variable
  *dcdims = (int**)malloc(sizeof(int*)*N_data);

  //calculating clen of each variable
  size_t k;
  for (size_t i=0; i<N_data; i++){

    if (Dcmp_data[i].opt[0] == 0){ //if decomposing the data as is, clen is simply the product of its dimensions (excluding time)
      if (gsum == 1){ k = 1; } //full phase-space decomposition: all dimensions except time go into clen
      else{ k = 2; } //groupwise decomposition: all dimensions except time and leading dimension go into clen
      clen_[i] = 1.; //initializing clen_ to 1

      for (size_t j=k; j<Dcmp_data[i].ndims; j++){ //taking product of dimension lengths
        clen_[i] = clen_[i]*dims[Dcmp_data[i].dimids[j]].len;
      }

    }
    else if (Dcmp_data[i].opt[0] == 1){ //if decomposing stacked data, clen depends on the summation of their dimensions
      clen_[i] = 0.; //initializing clen_ to 0
      //taking summation of stacking dimensions
      for (size_t j=0; j<(size_t)Dcmp_data[i].opt[1]; j++){
        clen_[i] = clen_[i] + dims[Dcmp_data[i].dimids[Dcmp_data[i].ndims-(size_t)Dcmp_data[i].opt[1]+j]].len;
      }

      if ((gsum == 1) && ((Dcmp_data[i].ndims-(size_t)Dcmp_data[i].opt[1]) != 1)){ //for full phase space decomposition, multiply stacked length with that of the leading dimension
        clen_[i] = clen_[i]*dims[Dcmp_data[i].dimids[1]].len;
      }

    }

    (*dcdims)[i] = (int*)malloc(sizeof(int)*2);
    (*dcdims)[i][0] = (int)clen_[i];
  }

  size_t N_clen; //N_clen = number of unique 'clen' dimensions
  Sort_Uniq_sizet(clen_,N_data,&N_clen); //sort clen_ to have unique values as the first N_clen elements

  //allocating memory for clen, rank dimension lists
  *clen = (ncdim*)malloc(sizeof(ncdim)*N_clen);
  *rank = (ncdim*)malloc(sizeof(ncdim)*N_clen);

  //filling clen, rank dimension descriptions
  for (size_t i=0; i<N_clen; i++){
    (*clen)[i].len = clen_[i];

    if (dcmp_type == 0){
      (*rank)[i].len = min((*clen)[i].len,dims[0].len);
    }
    else if (dcmp_type == 1){
      (*rank)[i].len = min((*clen)[i].len,dims[0].len-1);
    }
    else{
      printf("Invalid dcmp_type = %d found in Def_Dims!",dcmp_type);
      exit(1);
    }


    char buf[10];
    memset(buf,0,10);
    sprintf(buf,"clen%ld",i);
    strcpy((*clen)[i].name,buf);

    memset(buf,0,10);
    sprintf(buf,"rank%ld",i);
    strcpy((*rank)[i].name,buf);

  }
  free(clen_);

  //mapping each dataset to its respective clen, rank dimensions
  for (size_t i=0; i<N_data; i++){
    for (size_t j=0; j<N_clen; j++){
      if ((*dcdims)[i][0] == (int)(*clen)[j].len){
        (*dcdims)[i][0] = (int)j;
        (*dcdims)[i][1] = (int)j;

        break;
      }
    }
  }

  //writing original dimensions (from input datafile) to the output ncdf file and replacing id's from input file
  for (size_t i=0; i<N_dims; i++){
    err = nc_def_dim(ncid,dims[i].name,dims[i].len,&dims[i].id); Handle_Err(err,loc);
  }

  //writing clen, rank dimensions
  for (size_t i=0; i<N_clen; i++){
    err = nc_def_dim(ncid,(*clen)[i].name,(*clen)[i].len,&(*clen)[i].id); Handle_Err(err,loc);
    err = nc_def_dim(ncid,(*rank)[i].name,(*rank)[i].len,&(*rank)[i].id); Handle_Err(err,loc);
  }

  //termination
  return err;
}

/*================================================================================================================================*/
/**/
/*================================================================================================================================*/
int Def_Grids(const int ncid, Data *Dcmp_data, const size_t N_data, ncdim *dims, const size_t N_dims, Data *Grids, const size_t N_grids)
{
  int err;
  char loc[10] = "Def_Grids";

  for (size_t i=0; i<N_grids; i++){
    int *d = malloc(sizeof(int)*Grids[i].ndims);
    for (size_t j=0; j<Grids[i].ndims; j++){
      d[j] = dims[Grids[i].dimids[j]].id;
    }
    err = nc_def_var(ncid,Grids[i].name,NC_DOUBLE,(int)Grids[i].ndims,d,&Grids[i].id); Handle_Err(err,loc);
    free(d);

    err = nc_put_var_double(ncid,Grids[i].id,Grids[i].dat); Handle_Err(err,loc);
    err = nc_put_att_double(ncid,Grids[i].id,"bnds",NC_DOUBLE,2,Grids[i].bnds); Handle_Err(err,loc);
  }

  size_t *dcounts = malloc(sizeof(size_t)*N_data);
  for (size_t i=0; i<N_data; i++){
    dcounts[i] = Dcmp_data[i].ndims;
  }

  for (size_t i=0; i<N_data; i++){

    int id, id2;
    err = nc_def_var(ncid, Dcmp_data[i].name, NC_DOUBLE, 0, &id2, &id); Handle_Err(err,loc);
    err = nc_put_att_int(ncid, id, "N_dims", NC_INT, 1, (int*)&Dcmp_data[i].ndims); Handle_Err(err,loc);
    for (size_t j=0; j<Dcmp_data[i].ndims; j++){
      char buf[10];
      memset(buf,0,10);
      sprintf(buf,"dim%ld",j);
      err = nc_put_att_text(ncid, id, buf, 10, dims[Dcmp_data[i].dimids[j]].name); Handle_Err(err,loc);
    }

    if (Dcmp_data[i].opt[0] == 0){

      err = nc_put_att_int(ncid, id, "N_grids", NC_INT, 1, (int*)&Dcmp_data[i].ngrids); Handle_Err(err,loc);
      for (size_t j=0; j<Dcmp_data[i].ngrids; j++){
        char buf[10];
        memset(buf,0,10);
        sprintf(buf,"grid%ld",j);
        err = nc_put_att_text(ncid, id, buf, 10, Grids[Dcmp_data[i].gridids[j]].name); Handle_Err(err,loc);
      }

    }

  }

  //termination
  return err;
}

/*================================================================================================================================*/
/**/
/*================================================================================================================================*/
int Def_Disc(const int ncid, Data *Disc_Wts, const size_t N_wts, ncdim *dims)
{
  int err;
  char loc[9] = "Def_Disc";

  for (size_t i=0; i<N_wts; i++){

    if (Disc_Wts[i].ndims == 0){
      err = nc_def_var(ncid,Disc_Wts[i].name,NC_DOUBLE,0,0,&Disc_Wts[i].id); Handle_Err(err,loc);
    }
    else{
      int *d = malloc(sizeof(int)*Disc_Wts[i].ndims);
      for (size_t j=0; j<Disc_Wts[i].ndims; j++){
        d[j] = dims[Disc_Wts[i].dimids[j]].id;
      }
      err = nc_def_var(ncid,Disc_Wts[i].name,NC_DOUBLE,(int)Disc_Wts[i].ndims,d,&Disc_Wts[i].id); Handle_Err(err,loc);
      free(d);
    }

    err = nc_put_var_double(ncid,Disc_Wts[i].id,Disc_Wts[i].dat); Handle_Err(err,loc);

  }

  return err;
}

/*================================================================================================================================*/
/**/
/*================================================================================================================================*/
int Def_Specs(const int ncid, Spec *Prb_specs, const size_t N_specs)
{
  int err;

  for (size_t i=0; i<N_specs; i++){
    err = Write_Spec(ncid,&Prb_specs[i]);
  }

  return err;
}

/*================================================================================================================================*/
/**/
/*================================================================================================================================*/
int Def_POD_Vars(const int ncid, char *dname, const int rank_ID, const int clen_ID, const int N_t_ID, const int N_g_ID,
  const int gsum, Data *Decomp)
{
  int err;
  char loc[13] = "Def_POD_Vars";
  size_t ndims;
  int p1, p2;

  if (gsum == 1){ ndims = 2; p1 = 0; p2 = 1;}
  else{ ndims = 3; p1 = 1; p2 = 2;}

  for (size_t i=0; i<4; i++){
    Decomp[i].ndims = ndims;
    Decomp[i].dimids = (int *)malloc(sizeof(int)*ndims);
    Decomp[i].dimids[0] = N_g_ID;
  }

  Decomp[0].dimids[p1] = clen_ID;
  memset(Decomp[0].name,0,50);
  strcpy(Decomp[0].name,"C_");
  strcat(Decomp[0].name,dname);
  err = nc_def_var(ncid,Decomp[0].name,NC_DOUBLE,(int)ndims-1,Decomp[0].dimids,&Decomp[0].id); Handle_Err(err,loc);

  Decomp[1].dimids[p1] = rank_ID;
  memset(Decomp[1].name,0,50);
  strcpy(Decomp[1].name,"S_");
  strcat(Decomp[1].name,dname);
  err = nc_def_var(ncid,Decomp[1].name,NC_DOUBLE,(int)ndims-1,Decomp[1].dimids,&Decomp[1].id); Handle_Err(err,loc);

  Decomp[2].dimids[p1] = rank_ID; Decomp[2].dimids[p2] = clen_ID;
  memset(Decomp[2].name,0,50);
  strcpy(Decomp[2].name,"U_");
  strcat(Decomp[2].name,dname);
  err = nc_def_var(ncid,Decomp[2].name,NC_DOUBLE,(int)ndims,Decomp[2].dimids,&Decomp[2].id); Handle_Err(err,loc);

  Decomp[3].dimids[p1] = N_t_ID; Decomp[3].dimids[p2] = rank_ID;
  memset(Decomp[3].name,0,50);
  strcpy(Decomp[3].name,"Vt_");
  strcat(Decomp[3].name,dname);
  err = nc_def_var(ncid,Decomp[3].name,NC_DOUBLE,(int)ndims,Decomp[3].dimids,&Decomp[3].id); Handle_Err(err,loc);

  return 0;

}

/*================================================================================================================================*/
/**/
/*================================================================================================================================*/
int Def_DMD_Vars(const int ncid, char *dname, const int rank_ID, const int clen_ID, const int N_t_ID, const int N_g_ID,
  const int gsum, Data *Decomp)
{
  int err;
  char loc[13] = "Def_DMD_Vars";
  size_t ndims = 2;
  int p1, p2;

  if (gsum == 1){ ndims = 2; p1 = 0; p2 = 1;}
  else{ ndims = 3; p1 = 1; p2 = 2;}

  for (size_t i=0; i<6; i++){
    Decomp[i].ndims = ndims;
    Decomp[i].dimids = (int *)malloc(sizeof(int)*ndims);
    Decomp[i].dimids[0] = N_g_ID;
  }

  Decomp[0].dimids[p1] = rank_ID;
  memset(Decomp[0].name,0,50);
  strcpy(Decomp[0].name,"L_real_");
  strcat(Decomp[0].name,dname);
  err = nc_def_var(ncid,Decomp[0].name,NC_DOUBLE,(int)ndims-1,Decomp[0].dimids,&Decomp[0].id); Handle_Err(err,loc);

  Decomp[1].dimids[p1] = rank_ID;
  memset(Decomp[1].name,0,50);
  strcpy(Decomp[1].name,"L_imag_");
  strcat(Decomp[1].name,dname);
  err = nc_def_var(ncid,Decomp[1].name,NC_DOUBLE,(int)ndims-1,Decomp[1].dimids,&Decomp[1].id); Handle_Err(err,loc);

  Decomp[2].dimids[p1] = rank_ID; Decomp[2].dimids[p2] = clen_ID;
  memset(Decomp[2].name,0,50);
  strcpy(Decomp[2].name,"W_real_");
  strcat(Decomp[2].name,dname);
  err = nc_def_var(ncid,Decomp[2].name,NC_DOUBLE,(int)ndims,Decomp[2].dimids,&Decomp[2].id); Handle_Err(err,loc);

  Decomp[3].dimids[p1] = rank_ID; Decomp[3].dimids[p2] = clen_ID;
  memset(Decomp[3].name,0,50);
  strcpy(Decomp[3].name,"W_imag_");
  strcat(Decomp[3].name,dname);
  err = nc_def_var(ncid,Decomp[3].name,NC_DOUBLE,(int)ndims,Decomp[3].dimids,&Decomp[3].id); Handle_Err(err,loc);

  Decomp[4].dimids[p1] = rank_ID;
  memset(Decomp[4].name,0,50);
  strcpy(Decomp[4].name,"eL_real_");
  strcat(Decomp[4].name,dname);
  err = nc_def_var(ncid,Decomp[4].name,NC_DOUBLE,(int)ndims-1,Decomp[4].dimids,&Decomp[4].id); Handle_Err(err,loc);

  Decomp[5].dimids[p1] = rank_ID;
  memset(Decomp[5].name,0,50);
  strcpy(Decomp[5].name,"eL_imag_");
  strcat(Decomp[5].name,dname);
  err = nc_def_var(ncid,Decomp[5].name,NC_DOUBLE,(int)ndims-1,Decomp[5].dimids,&Decomp[5].id); Handle_Err(err,loc);

  return 0;

}

/*================================================================================================================================*/
/**/
/*================================================================================================================================*/
int Def_DCMP_Vars(const int ncid, const int dcmp_type, const int gsum, Data *Dcmp_data, const size_t N_data, ncdim *clen,
  ncdim *rank, ncdim *dims, Data **Decomp, int **dcdims)
{

  int err;

  if (dcmp_type == 0 && gsum == 1){ err = nc_put_att_text(ncid,NC_GLOBAL,"dcmp_type",4,"POD"); }
  else if (dcmp_type == 0 && gsum == 0){ err = nc_put_att_text(ncid,NC_GLOBAL,"dcmp_type",5,"PODg"); }
  else if (dcmp_type == 1 && gsum == 1){ err = nc_put_att_text(ncid,NC_GLOBAL,"dcmp_type",4,"DMD"); }
  else if (dcmp_type == 1 && gsum == 0){ err = nc_put_att_text(ncid,NC_GLOBAL,"dcmp_type",5,"DMDg"); }

  /*------------------------------------------------------------/
  /                             POD                             /
  /------------------------------------------------------------*/
  /*POD uses 4 variables: ID[0] = C (centering vector),
                          ID[1] = S (singular values)
                          ID[2] = U (left singular vectors)
                          ID[3] = Vt (right singular vectors)
  */
  if (dcmp_type == 0){ //decompose with the POD

    *Decomp = (Data*)malloc(sizeof(Data)*N_data*4);
    size_t p = 0;
    for (size_t i=0; i<N_data; i++){
      err = Def_POD_Vars(ncid, Dcmp_data[i].name, rank[dcdims[i][1]].id, clen[dcdims[i][0]].id,
        dims[0].id, dims[1].id, gsum, &(*Decomp)[p]);
      p = p + 4;

    }
  }
  /*------------------------------------------------------------/
  /                             DMD                             /
  /------------------------------------------------------------*/
  /*DMD uses 4 variables: ID[0] = L_real (DMD eigenvalues, real component),
                          ID[1] = L_imag (DMD eigenvalues, imaginary component),
                          ID[2] = W_real (DMD modes/ eigenvectors, real component)
                          ID[3] = W_imag (DMD modes/ eigenvectors, imaginary component)
  */
  else if (dcmp_type == 1){ //decompose with the DMD

    *Decomp = (Data*)malloc(sizeof(Data)*N_data*6);
    size_t p = 0;
    for (size_t i=0; i<N_data; i++){
      err = Def_DMD_Vars(ncid, Dcmp_data[i].name, rank[dcdims[i][1]].id, clen[dcdims[i][0]].id,
        dims[0].id, dims[1].id, gsum, &(*Decomp)[p]);
      p = p + 6;

    }

  }

  return err;

}

/*================================================================================================================================*/
/**/
/*================================================================================================================================*/
int Generate_DCMP(const double *data, const int ncid_out, const char *dname, const size_t N_t, const size_t N_g, const size_t clen,
  const size_t rank, const double svd_eps, Data *Decomp, const int dcmp_type, const double delt)
{
  int err;

  if (dcmp_type == 0){ //decompose with the POD
    err = Generate_POD(data,ncid_out,dname,N_t,N_g,clen,rank,Decomp);
  }
  else if (dcmp_type == 1){ //decompose with the DMD
    err = Generate_DMD(data,ncid_out,dname,N_t,N_g,clen,rank,svd_eps,Decomp,delt);
  }

  return err;
}

/*================================================================================================================================*/
/**/
/*================================================================================================================================*/
// int Output_meanI_DCMP(const int ncid_in, const int ncid_out, const size_t N_t, const size_t N_g, const size_t N_m, const size_t N_y,
//   const size_t N_x, const size_t rank_avg, const size_t rank_edgV, const size_t rank_edgH, const double svd_eps, const int gsum,
//   const int *Ig_avg_IDs, const int *Ig_edgV_IDs, const int *Ig_edgH_IDs, const int dcmp_type)
// {
//   int err;
//   char dname[25];
//   size_t clen, n_g, gscale, len, len2, p1, p2;
//   double *data, *data2;
//   double c;
//
//   c = 299.792458;
//   if (gsum == 1){
//     n_g = 0; gscale = N_g;
//   }
//   else{
//     n_g = N_g; gscale = 1;
//   }
//
//   /*------------------------------------------------------------/
//   /       Average I_Avg with scalar intensities, decompose      /
//   /------------------------------------------------------------*/
//   printf("... mean I_avg start\n");
//   memset(dname,0,25); strcpy(dname,"Eg_avg_HO");
//   len = N_t*N_g*N_y*N_x; Get_Var_Double(ncid_in,dname,&data2,len);
//
//   memset(dname,0,25); strcpy(dname,"I_avg");
//   len = N_t*N_g*N_m*N_y*N_x; Get_Var_Double(ncid_in,dname,&data,len);
//
//   len = N_t*N_g;
//   len2 = N_y*N_x;
//   p1 = 0;
//   for(size_t i=0; i<len; i++){
//     for(size_t m=0; m<N_m; m++){
//       p2 = i*len2;
//       for(size_t j=0; j<len2; j++){
//
//         data[p1] = data[p1]/(data2[p2]*c);
//         p1 = p1 + 1;
//         p2 = p2 + 1;
//
//       }
//     }
//   }
//   free(data2);
//
//   memset(dname,0,25); strcpy(dname,"Mean_I_avg"); clen = gscale*N_m*N_y*N_x;
//   err = Generate_DCMP(data,ncid_out,dname,N_t,n_g,clen,rank_avg,svd_eps,Ig_avg_IDs,dcmp_type);
//   free(data);
//
//
//   /*------------------------------------------------------------/
//   /      Average I_EdgV with scalar intensities, decompose      /
//   /------------------------------------------------------------*/
//   printf("... mean I_edgV start\n");
//   memset(dname,0,25); strcpy(dname,"Eg_edgV_HO");
//   len = N_t*N_g*N_y*(N_x+1); Get_Var_Double(ncid_in,dname,&data2,len);
//
//   memset(dname,0,25); strcpy(dname,"I_edgV");
//   len = N_t*N_g*N_m*N_y*(N_x+1); Get_Var_Double(ncid_in,dname,&data,len);
//
//   len = N_t*N_g;
//   len2 = N_y*(N_x+1);
//   p1 = 0;
//   for(size_t i=0; i<len; i++){
//     for(size_t m=0; m<N_m; m++){
//       p2 = i*len2;
//       for(size_t j=0; j<len2; j++){
//
//         data[p1] = data[p1]/(data2[p2]*c);
//         p1 = p1 + 1;
//         p2 = p2 + 1;
//
//       }
//     }
//   }
//   free(data2);
//
//   memset(dname,0,25); strcpy(dname,"Mean_I_edgV"); clen = gscale*N_m*N_y*(N_x+1);
//   err = Generate_DCMP(data,ncid_out,dname,N_t,n_g,clen,rank_edgV,svd_eps,Ig_edgV_IDs,dcmp_type);
//   free(data);
//
//   /*------------------------------------------------------------/
//   /      Average I_EdgH with scalar intensities, decompose      /
//   /------------------------------------------------------------*/
//   printf("... mean I_edgH start\n");
//   memset(dname,0,25); strcpy(dname,"Eg_edgH_HO");
//   len = N_t*N_g*(N_y+1)*N_x; Get_Var_Double(ncid_in,dname,&data2,len);
//
//   memset(dname,0,25); strcpy(dname,"I_edgH");
//   len = N_t*N_g*N_m*(N_y+1)*N_x; Get_Var_Double(ncid_in,dname,&data,len);
//
//   len = N_t*N_g;
//   len2 = (N_y+1)*N_x;
//   p1 = 0;
//   for(size_t i=0; i<len; i++){
//     for(size_t m=0; m<N_m; m++){
//       p2 = i*len2;
//       for(size_t j=0; j<len2; j++){
//
//         data[p1] = data[p1]/(data2[p2]*c);
//         p1 = p1 + 1;
//         p2 = p2 + 1;
//
//       }
//     }
//   }
//   free(data2);
//
//   memset(dname,0,25); strcpy(dname,"Mean_I_edgH"); clen = gscale*N_m*(N_y+1)*N_x;
//   err = Generate_DCMP(data,ncid_out,dname,N_t,n_g,clen,rank_edgH,svd_eps,Ig_edgH_IDs,dcmp_type);
//   free(data);
//
//   return err;
// }

/*================================================================================================================================*/
/**/
/*================================================================================================================================*/
int Decompose_Data(const int ncid_in, const int ncid_out, const int dcmp_type, const int gsum, Data *Dcmp_data,
  const size_t N_data, Data *Grids, const size_t N_grids, Data *Decomp, ncdim *clen, ncdim *rank, ncdim *dims,
  int **dcdims, const double svd_eps)
{
  int err;
  char loc[15] = "Decompose_Data";
  size_t n_g;
  double delt;

  if (gsum == 1){ n_g = 0; }
  else{ n_g = dims[1].len; }

  size_t p = 0;
  for (size_t i=0; i<N_data; i++){
    printf("Decomposing %s\n",Dcmp_data[i].name);

    //
    //
    //
    delt = Grids[Dcmp_data[i].gridids[0]].dat[1] - Grids[Dcmp_data[i].gridids[0]].dat[0];
    // exit(0);
    //
    //
    //

    /*------------------------------------------------------------/
    /                     As-is decomposition                     /
    /------------------------------------------------------------*/
    if (Dcmp_data[i].opt[0] == 0){
      printf("... Performing decomposition of data as input\n");

      size_t dlen = 1;
      for (size_t j=0; j<Dcmp_data[i].ndims; j++){
        dlen = dlen*dims[Dcmp_data[i].dimids[j]].len;
      }
      nc_inq_varname(ncid_in,Dcmp_data[i].id,Dcmp_data[i].name);

      Dcmp_data[i].dat = (double*)malloc(sizeof(double)*dlen);
      err = nc_get_var_double(ncid_in,Dcmp_data[i].id,&Dcmp_data[i].dat[0]); Handle_Err(err,loc);

      err = Generate_DCMP(Dcmp_data[i].dat, ncid_out, Dcmp_data[i].name, dims[0].len, n_g, clen[dcdims[i][0]].len,
        rank[dcdims[i][1]].len, svd_eps, &Decomp[p], dcmp_type, delt);

      free(Dcmp_data[i].dat);

    }
    /*------------------------------------------------------------/
    /               Decomposition of 'stacked' data               /
    /------------------------------------------------------------*/
    else if (Dcmp_data[i].opt[0] == 1){
      printf("... Performing decomposition of stacked data\n");

      size_t dlen = 0.; //initializing dlen to 0
      //adding 'stacking' dimensions together
      for (size_t j=0; j<(size_t)Dcmp_data[i].opt[1]; j++){
        dlen = dlen + dims[Dcmp_data[i].dimids[Dcmp_data[i].ndims-(size_t)Dcmp_data[i].opt[1]+j]].len;
      }
      //multiplying stacked dimension by leading dimensions
      for (size_t j=0; j<(Dcmp_data[i].ndims-(size_t)Dcmp_data[i].opt[1]); j++){
        dlen = dlen*dims[Dcmp_data[i].dimids[j]].len;
      }

      Dcmp_data[i].dat = (double*)malloc(sizeof(double)*dlen); //allocating space for data matrix

      //allocating 2D character array to hold names of variables to stack together
      char **stk_names = (char**)malloc(sizeof(char*)*(size_t)Dcmp_data[i].opt[1]);
      for (size_t j=0; j<(size_t)Dcmp_data[i].opt[1]; j++){
        stk_names[j] = (char*)malloc(sizeof(char)*50);
      }

      err = Delimit(Dcmp_data[i].cdat,',',stk_names,Dcmp_data[i].opt[1]); //finding names of stacking variables from comma delimited list

      //calculating length of leading dimensions (product)
      size_t lead_len = 1.;
      for (size_t j=0; j<(Dcmp_data[i].ndims-(size_t)Dcmp_data[i].opt[1]); j++){
        lead_len = lead_len*dims[Dcmp_data[i].dimids[j]].len;
      }

      //loop over number of stacking variables
      for (size_t j=0; j<(size_t)Dcmp_data[i].opt[1]; j++){

        //stk_len is the length of the current stacking dimension
        size_t stk_len = dims[Dcmp_data[i].dimids[Dcmp_data[i].ndims-(size_t)Dcmp_data[i].opt[1]+j]].len;
        dlen = stk_len*lead_len; //dlen = total dimensional length of stacking variable

        double *data_ = malloc(sizeof(double)*dlen); //allocating space to hold stacking data
        int id;
        err = nc_inq_varid(ncid_in,stk_names[j],&id); Handle_Err(err,loc); //finding variable
        err = nc_get_var_double(ncid_in,id,&data_[0]); Handle_Err(err,loc); //reading in data

        //calculating len_t = 'top length' (length of stacking dimensions for variables stacked above current variable)
        size_t len_t = 0.;
        for (size_t k=0; k<j; k++){
          len_t = len_t + dims[Dcmp_data[i].dimids[Dcmp_data[i].ndims-(size_t)Dcmp_data[i].opt[1]+k]].len;
        }

        //calculating len_b = 'bottom length' (length of stacking dimensions for variables stacked below current variable)
        size_t len_b = 0.;
        for (size_t k=j+1; k<(size_t)Dcmp_data[i].opt[1]; k++){
          len_b = len_b + dims[Dcmp_data[i].dimids[Dcmp_data[i].ndims-(size_t)Dcmp_data[i].opt[1]+k]].len;
        }

        //putting data into stacked data matrix
        size_t p1 = 0;
        size_t p2 = 0;
        for(size_t k=0; k<lead_len; k++){
          p1 = p1 + len_t;
          for(size_t l=0; l<stk_len; l++){
            Dcmp_data[i].dat[p1] = data_[p2];
            p1 = p1 + 1;
            p2 = p2 + 1;
          }
          p1 = p1 + len_b;
        }
        free(data_);

      }//end loop over number of stacking variables

      //decomposing stacked data matrix
      err = Generate_DCMP(Dcmp_data[i].dat, ncid_out, Dcmp_data[i].name, dims[0].len, n_g, clen[dcdims[i][0]].len,
        rank[dcdims[i][1]].len, svd_eps, &Decomp[p], dcmp_type, delt);

      //freeing all allocated space
      free(Dcmp_data[i].dat);
      for (size_t j=0; j<(size_t)Dcmp_data[i].opt[1]; j++){
        free(stk_names[j]);
      }
      free(stk_names);

    }
    /*------------------------------------------------------------/
    /           Error - Unrecognized decomposition type           /
    /------------------------------------------------------------*/
    else{
      printf("Invalid Dcmp_data.opt = %d detected in Decompose_Data!\n",Dcmp_data[i].opt[0]);
      exit(1);
    }

    p = p + 4;

  }

  return err;
}
