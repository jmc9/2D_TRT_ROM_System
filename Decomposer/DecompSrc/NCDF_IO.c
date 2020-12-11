/*================================================================================================================================*/
/*
  NCDF_IO.c
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
/**/
/*================================================================================================================================*/
void Handle_Err(const int Status, const char *Location)
{
  if (Status != NC_NOERR){
    printf("Error! [location %s]\n",Location);
    printf("NCDF error message: %s\n",nc_strerror(Status));
  }
}

/*================================================================================================================================*/
/**/
/*================================================================================================================================*/
void Get_Var_Double(const int ncid, const char *name, double **var, const size_t size)
{
  int err, vID;
  char loc[15] = "GET_VAR_DOUBLE";

  err = nc_inq_varid(ncid,name,&vID); Handle_Err(err,loc);

  *var = (double *)malloc(sizeof(double)*size);

  err = nc_get_var_double(ncid,vID,var[0]); Handle_Err(err,loc);

}

/*================================================================================================================================*/
/**/
/*================================================================================================================================*/
void Open_NCfile(char *fname, int *ncid)
{
  int err;
  char loc[12] = "OPEN_NCFILE";

  err = nc_create(fname, NC_CLOBBER, &(*ncid)); Handle_Err(err,loc);

}

/*================================================================================================================================*/
/*  */
/*================================================================================================================================*/
int Get_Spec(const int ncid, Spec *spec)
{
  int err, vID;
  char loc[8] = "Get_Spec";

  union
  {
    int i;
    double d;
  } val;

  err = nc_inq_varid(ncid,spec->name,&vID); Handle_Err(err,loc);
  err = nc_inq_vartype(ncid,vID,&spec->type); Handle_Err(err,loc);

  if (spec->type == sp_dbl){
    spec->data = malloc(sizeof(double));
    err = nc_get_var_double(ncid,vID,&val.d); Handle_Err(err,loc);
    *(double*)spec->data = val.d;
  }
  else if (spec->type == sp_int){
    spec->data = malloc(sizeof(int));
    err = nc_get_var_int(ncid,vID,&val.i); Handle_Err(err,loc);
    *(int*)spec->data = val.i;
  }
  else{
    printf("Get_Spec only handles doubles, ints right now\n");
    return 1;
  }

  return err;
}
