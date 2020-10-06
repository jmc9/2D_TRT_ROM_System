#include <stdlib.h>
#include <stdio.h>
#include <netcdf.h>
#include <string.h>

// int TEST(int a, int b);
void input(char *infile, char *dsfile);
void copy(char to[], char from[]);

/* Handle NetCDF errors by printing an error message and exiting with a
 * non-zero status. */
// #define ERRCODE 2
#define NC_ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(2);}
#define NC_NOERR   0

//================================================================================================================================//
//
//================================================================================================================================//
int main()
{
  // comment
  int ncid, err;
  int fg_avg_xx_ID, fg_avg_xy_ID, fg_avg_yy_ID, fg_edgV_xx_ID, fg_edgV_xy_ID, fg_edgH_yy_ID, fg_edgH_xy_ID;
  char infile[9], dsfile[100];
  // double fg_avg_xx;

  copy(infile,"input.inp"); //setting input file name (default)

  input(infile,dsfile); //reading input file
  printf("%s\n",dsfile);

  err = nc_open(dsfile,NC_NOWRITE,&ncid); //opening NetCDF dataset
  if(err != NC_NOERR) NC_ERR(err)

  //reading in multigroup qd factors from the dataset
  err = nc_inq_varid(ncid,"fg_avg_xx",&fg_avg_xx_ID); if(err != NC_NOERR) NC_ERR(err);
  err = nc_inq_varid(ncid,"fg_avg_xy",&fg_avg_xy_ID); if(err != NC_NOERR) NC_ERR(err);
  err = nc_inq_varid(ncid,"fg_avg_yy",&fg_avg_yy_ID); if(err != NC_NOERR) NC_ERR(err);
  err = nc_inq_varid(ncid,"fg_edgV_xx",&fg_edgV_xx_ID); if(err != NC_NOERR) NC_ERR(err);
  err = nc_inq_varid(ncid,"fg_edgV_xy",&fg_edgV_xy_ID); if(err != NC_NOERR) NC_ERR(err);
  err = nc_inq_varid(ncid,"fg_edgH_yy",&fg_edgH_yy_ID); if(err != NC_NOERR) NC_ERR(err);
  err = nc_inq_varid(ncid,"fg_edgH_xy",&fg_edgH_xy_ID); if(err != NC_NOERR) NC_ERR(err);

  //
  // err = nc_get_var_int(ncid,fg_avg_xx_ID,&fg_avg_xx[0][0][0]); if(err != NC_NOERR) NC_ERR(err);

  err = nc_close(ncid);
  if(err != NC_NOERR) NC_ERR(err) //closing NetCDF dataset

  printf("yuh\n");
}

//================================================================================================================================//
//
//================================================================================================================================//
void input(char infile[], char dsfile[])
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
