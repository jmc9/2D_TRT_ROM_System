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
  while ((line[i] != '\0') && (line[i] != '\n') && (line[i] != (char)13)){

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

/*================================================================================================================================*/
/**/
/*================================================================================================================================*/
void Sort_Uniq_int(int *list, const size_t len, size_t *ulen)
{
  *ulen = 1;
  for (size_t i=0; i<len; i++){
    for (size_t j=i+1; j<len; j++){

      int new = 1;
      for (size_t k=0; k<i+1; k++){
        if (list[j] == list[k]){
          new = 0;
          break;
        }
      }

      if (new == 1){
        int temp = list[i+1];
        list[i+1] = list[j];
        list[j] = temp;
        *ulen = *ulen + 1;
        break;
      }

    }
  }
}

void Sort_Uniq_sizet(size_t *list, const size_t len, size_t *ulen)
{
  *ulen = 1;
  for (size_t i=0; i<len; i++){
    for (size_t j=i+1; j<len; j++){

      int new = 1;
      for (size_t k=0; k<i+1; k++){
        if (list[j] == list[k]){
          new = 0;
          break;
        }
      }

      if (new == 1){
        size_t temp = list[i+1];
        list[i+1] = list[j];
        list[j] = temp;
        *ulen = *ulen + 1;
        break;
      }

    }
  }
}
