/*================================================================================================================================*/
/*
  GNUPLOT_ROUTINES.c
  Author: Joseph M. Coale
  jmcoale@ncsu.edu - josephcoale912@gmail.com
*/
/*================================================================================================================================*/
/* Importing Packages */
/*================================================================================================================================*/
#include <stdlib.h>
#include <stdio.h>

/*================================================================================================================================*/
/*
  function: gnuplot_1d

  this function opens a pipe to gnuplot and graphs a 1D vector as a line or scatter plot

  returns an integer error code:
    0 = successful
    1 = gnuplot failed to open
    2 = gnuplot failed to close
    3 = incorrect input to function detected

  inputs:
  char *title  = title of the plot
    double *data = vector of data to be plotted (y coordinates)
    double *crd  = vector of coordinates of each data point (x coordinates)
    int dim      = dimension/length of data and crd vectors
    char *plttyp = character array containing the desired plot type ("p","l","lp")
    int logscale = base of the log to generate a semilog y-axis with (if logscale = 1, linear y axis is used)
    int pt       = gnuplot point type
    char *lc     = character array containing gnuplot line color
*/
/*================================================================================================================================*/
int gnuplot_1d(const char *title, const double *data, const double *crd, const int dim, const char *plttyp, const int logscale,
  const int pt, const char *lc, const char *saveas, const double *xbnd)
{
  FILE *gnuplot_Pipe;
  int err, i;

  gnuplot_Pipe = popen("gnuplot","w"); //opening pipe to gnuplot
  //check if gnuplot opens correctly
  if(gnuplot_Pipe == NULL){
    return 1;
  }

  fprintf(gnuplot_Pipe, "set term png\n"); //setting gnu to png mode
  fprintf(gnuplot_Pipe, "set output \"%s\"\n",saveas); //telling gnu to output plot in png format

  fprintf(gnuplot_Pipe, "set title \"%s\"\n",title); //setting plot title

  //check if using semilog format
  if(logscale > 1){ //using semilog format - set logscale and reformat y axis ticks
    fprintf(gnuplot_Pipe, "set logscale y %d\n", logscale);
    fprintf(gnuplot_Pipe, "set format y \"\%.0s10^{\%T}\"\n");
  }
  else if(logscale < 1){ //negative logscale impossible, terminate with errorcode 3
    return 3;
  }

  fprintf(gnuplot_Pipe, "set xrange [%le:%le]\n",xbnd[0],xbnd[1]); //moving x-axis bounds to just outside data domain

  fprintf(gnuplot_Pipe, "set grid y\n"); //turning on y-grid

  /* beginning plot
     parameter with x:
       x -> p = points (scatter plot); l = lines (line plot); lp = line points (line plot with points)
     parameter pt x: (pt = pointtype)
       x -> integer to specify point type (7 = filled circle)
     parameter lc x: (lc = linecolor)
       x -> string specifying color, i.e. "black"  */
  fprintf(gnuplot_Pipe, "plot '-' with p pt %d lc \"%s\"\n", pt, lc);

  //looping through data points
  for(i=0;i<(int)dim;i++){
    fprintf(gnuplot_Pipe, "%le %le\n", crd[i], data[i]);
  }

  fprintf(gnuplot_Pipe, "e"); //finishing plot, outputting

  err = pclose(gnuplot_Pipe); //closing pipe to gnuplot
  //check if gnuplot closes correctly
  if(err == -1){
    return 2;
  }

  return 0; //function completed successfully
}

/*================================================================================================================================*/
/**/
/*================================================================================================================================*/
void GNUp_Err(const int err)
{
  if(err != 0){
    if(err == 1){ printf("Failed to open gnuplot_Pipe"); exit(2); }
    if(err == 2){ printf("Failed to close gnuplot_Pipe"); exit(2); }
  }
}
