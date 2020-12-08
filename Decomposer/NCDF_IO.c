/*================================================================================================================================*/
/*
  NCDF_IO.c
  Author: Joseph M. Coale
  jmcoale@ncsu.edu - josephcoale912@gmail.com
*/
/*================================================================================================================================*/
/* Importing Packages */
/*================================================================================================================================*/
#include <stdlib.h>
#include <stdio.h>
#include <netcdf.h>

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
