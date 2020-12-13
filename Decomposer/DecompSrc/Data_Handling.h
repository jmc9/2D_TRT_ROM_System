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

typedef struct Data
{
  char name[50];
  size_t ndims;
  int *dimids;
  int opt[3];
  char cdat[50];
  double *dat;
  int id;
} Data;

#define sp_int NC_INT
#define sp_dbl NC_DOUBLE
