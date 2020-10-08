#include <stdlib.h>
#include <stdio.h>
#include <netcdf.h>
#include <string.h>

void INPUT(char *infile, char *dsfile);
void copy(char to[], char from[]);
void HANDLE_ERR(int Status, char Location[]);

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
void INPUT(char infile[], char dsfile[])
{
  FILE *inpf;
  char line[256];
  char *inps;

  inpf = fopen(infile,"r"); //opening input file
  //check if input file exists, if not terminate program
  if (inpf == NULL){
    printf("Error occured whilst opening file: %s",infile);
    exit(1);
  }

  //moving line by line through file to grab inputs
  while ( fgets(line, 255, inpf) != NULL ){ //placing current line of input file into 'line' character array

    inps = strtok(line," "); //splitting line into seperate strings, delimited by a space
    //looping through each delimited string
    while (inps != NULL){
      if(strcmp(inps,"dataset")){
        copy(dsfile,inps);
      }
      inps = strtok(NULL," "); //moving 'inps' to the next delimited part of line
    }

  }

  fclose(inpf); //closing input file
}

//================================================================================================================================//
/* copy: copy 'from' into 'to'; assume to is big enough */
//================================================================================================================================//
void copy(char to[], char from[])
{
  int i;
  i = 0;
  while ((to[i] = from[i]) != '\0')
  ++i;
}
