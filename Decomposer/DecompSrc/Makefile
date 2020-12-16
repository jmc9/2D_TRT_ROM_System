# Makefile

INC_DIR = /usr/include

#--------------------------------------------

CC = gcc
CFLAGS = -O0 -g
CWARNINGS =  -Wall -Wextra -Wconversion
CRUNTIMES = -fbounds-check #-fsanitize=bounds

NCDF = -lnetcdf
LAPACK = -llapacke -lblas

#--------------------------------------------

SRC_OBJ = main.o \
					INPUTS.o \
					NCDF_IO.o \
					MISC_PROCS.o \
					OUTPUTS.o \
					POD_ROUTINES.o \
					DMD_ROUTINES.o \
					GNUPLOT_ROUTINES.o \
					LA_ROUTINES.o 

EXEC = decomposer.exe

.SUFFIXES :
.SUFFIXES : .c .o

.PHONY : default all clean

default : all

all :
	$(MAKE) $(EXEC)

$(EXEC) : $(SRC_OBJ)
	$(FC) $^ -o $@ $(NCDF) $(LAPACK)
	ls -fl $@

.c.o :
	$(CC) $(CFLAGS) $(CWARNINGS) $(CRUNTIMES) -c -I$(INC_DIR) $<

clean :
	rm -f *.o

main.o :  NCDF_IO.o INPUTS.o OUTPUTS.o Data_Handling.h
INPUTS.o : NCDF_IO.o Data_Handling.h
NCDF_IO.o : Data_Handling.h
MISC_PROCS.o :
OUTPUTS.o : POD_ROUTINES.o DMD_ROUTINES.o Data_Handling.h
POD_ROUTINES.o : GNUPLOT_ROUTINES.o MISC_PROCS.o LA_ROUTINES.o Data_Handling.h
GNUPLOT_ROUTINES.o :
LA_ROUTINES.o :
DMD_ROUTINES.o : GNUPLOT_ROUTINES.o MISC_PROCS.o LA_ROUTINES.o Data_Handling.h
