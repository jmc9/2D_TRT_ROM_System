#include <stdlib.h>
#include <stdio.h>
#include <netcdf.h>

// #define NC_NOERR   0;

//================================================================================================================================//
//
//================================================================================================================================//
void HANDLE_ERR(int Status, char Location[])
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
void GET_VAR_DOUBLE(int ncid, char name[], double **var, size_t size)
{
  int err, vID;
  char loc[14] = "GET_VAR_DOUBLE";

  err = nc_inq_varid(ncid,name,&vID); HANDLE_ERR(err,loc);

  *var = (double *)malloc(sizeof(double)*size);

  err = nc_get_var_double(ncid,vID,var[0]); HANDLE_ERR(err,loc);

}
