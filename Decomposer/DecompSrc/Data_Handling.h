/*================================================================================================================================*/
/*
  Data_Handling.h
  Author: Joseph M. Coale
  jmcoale@ncsu.edu - josephcoale912@gmail.com
*/
/*================================================================================================================================*/
/*  */
/*================================================================================================================================*/
#include <netcdf.h>

typedef struct Dataset
{
  char name[20];
  int ncID;
  int *dimIDs;
  size_t *dims;
  double *data;
} Dataset;

typedef struct Spec
{
  char name[20];
  int type;
  void *data;
} Spec;

typedef struct ncdim
{
  char name[10];
  size_t len;
  int id;
} ncdim;

#define sp_int NC_INT
#define sp_dbl NC_DOUBLE
