#include <stdlib.h>
#include <stdio.h>
#include <netcdf.h>
#include <string.h>

void HANDLE_ERR(int Status, char Location[]);

int delimit(char *line, char del, char **parts, int nparts);

//================================================================================================================================//
//
//================================================================================================================================//
void GET_DIMS(int ncid, size_t *N_t, size_t *N_g, size_t *N_y, size_t *N_x)
{
  int err;
  int N_t_ID, N_g_ID, N_y_ID, N_x_ID;
  char loc[8] = "GET_DIMS";

  //reading in dimension ID's
  err = nc_inq_dimid(ncid,"N_t",&N_t_ID); HANDLE_ERR(err,loc);
  err = nc_inq_dimid(ncid,"N_g",&N_g_ID); HANDLE_ERR(err,loc);
  err = nc_inq_dimid(ncid,"N_y",&N_y_ID); HANDLE_ERR(err,loc);
  err = nc_inq_dimid(ncid,"N_x",&N_x_ID); HANDLE_ERR(err,loc);

  //reading in dimension values
  err = nc_inq_dimlen(ncid,N_t_ID,N_t); HANDLE_ERR(err,loc);
  err = nc_inq_dimlen(ncid,N_g_ID,N_g); HANDLE_ERR(err,loc);
  err = nc_inq_dimlen(ncid,N_y_ID,N_y); HANDLE_ERR(err,loc);
  err = nc_inq_dimlen(ncid,N_x_ID,N_x); HANDLE_ERR(err,loc);

}

//================================================================================================================================//
//
//================================================================================================================================//
void INPUT(char *infile, char *dsfile, char *outfile)
{
  FILE *inpf;
  char line[256];
  char **inps, delim = ' ';
  int n, i, err;
  size_t nargs = 2;

  inpf = fopen(infile,"r"); //opening input file
  //check if input file exists, if not terminate program
  if (inpf == NULL){
    printf("Error occured whilst opening file: %s",infile);
    exit(1);
  }

  inps = (char **)malloc(sizeof(char *)*nargs);
  for (n=0; n < (int)nargs; n++){
    inps[n] = (char *)malloc(sizeof(char)*25);
  }

  //moving line by line through file to grab inputs
  while ( fgets(line, 255, inpf) != NULL ){ //placing current line of input file into 'line' character array

    for (n=0; n < (int)nargs; n++){
      for (i=0; i<25; i++){
        inps[n][i] = 0;
      }
    }
    
    err = delimit(line,delim,inps,(int)nargs);

    if (strcmp(inps[0],"dataset") == 0){
      strcpy(dsfile,inps[1]);
    }
    else if (strcmp(inps[0],"outfile") == 0){
      strcpy(outfile,inps[1]);
    }

  }

  fclose(inpf); //closing input file

}
