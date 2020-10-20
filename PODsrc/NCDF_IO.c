#include <stdlib.h>
#include <stdio.h>
#include <netcdf.h>

//================================================================================================================================//
//
//================================================================================================================================//
void HANDLE_ERR(const int Status, const char *Location)
{
  if (Status != NC_NOERR){
    printf("***   NETCDF ERROR ENCOUNTERED   ***\n");
    printf("***   Error occured in %s\n",Location);
    printf("***   %s\n",nc_strerror(Status));
    exit(2);
  }
}

//================================================================================================================================//
//
//================================================================================================================================//
void GET_VAR_DOUBLE(const int ncid, const char *name, double **var, const size_t size)
{
  int err, vID;
  char loc[15] = "GET_VAR_DOUBLE";

  err = nc_inq_varid(ncid,name,&vID); HANDLE_ERR(err,loc);

  *var = (double *)malloc(sizeof(double)*size);

  err = nc_get_var_double(ncid,vID,var[0]); HANDLE_ERR(err,loc);

}

//================================================================================================================================//
//
//================================================================================================================================//
void OPEN_NCFILE(char *fname, int *ncid)
{
  int err;
  char loc[12] = "OPEN_NCFILE";
  err = nc_create(fname, NC_CLOBBER, &(*ncid)); HANDLE_ERR(err,loc);
}
