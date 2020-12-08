/*================================================================================================================================*/
/*
  MISC_PROCS.c
  Author: Joseph M. Coale
  jmcoale@ncsu.edu - josephcoale912@gmail.com
*/
/*================================================================================================================================*/
/* Importing Packages */
/*================================================================================================================================*/
#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <string.h>

/*================================================================================================================================*/
/**/
/*================================================================================================================================*/
int Make_Dir(const char *dir)
{
  int err = 0;
  char buf[100];
  struct stat st;

  //first, check if directory already exists
  if (stat(dir,&st) == -1){ //if directory does not exist, create
    err = mkdir(dir,0777);
  }
  else{ //if directory does exist, delete everything in it
    sprintf(buf,"rm -rf %s/*",dir);
    system(buf);
  }

  return err;
}

/*================================================================================================================================*/
/**/
/*================================================================================================================================*/
int Delimit(const char *line, const char del, char **parts, const int nparts)
{
  int err = 0;
  int i, j, n;

  i=0; j=0; n=0;
  while (line[i] != '\0' && line[i] != '\n' && line[i] != (char)13){

    if (line[i] != del){
      parts[n][j] = line[i];
    }
    else{
      parts[n][j] = '\0';
      n++; j=-1;
      if (n == nparts){
        return 0;
      }
    }

    i++; j++;
  }

  return err;
}
