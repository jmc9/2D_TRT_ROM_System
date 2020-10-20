# Makefile

INC_DIR = /usr/include

#--------------------------------------------

CC = gcc
CFLAGS = -O0 -g
CWARNINGS =  -Wall -Wextra -Wconversion
CRUNTIMES = -fbounds-check #-fsanitize=bounds

FC = gfortran
FFLAGS = -O0 -g
FWARNINGS =  -Wall -Wextra -Wunused-dummy-argument -Wconversion -Wconversion-extra -Wintrinsic-shadow -Wsurprising
FRUNTIMES = -fbacktrace -fbounds-check -fcheck=all

NCDF = -lnetcdf
LAPACK = -llapack -lblas

#--------------------------------------------

SRC_OBJ = main.o \
					INPUTS.o \
					NCDF_IO.o \
					GNUPLOT_ROUTINES.o \
          POD_Routines.o \
					CPOD_ROUTINES.o \
					MISC_PROCS.o \
					OUTPUTS.o

EXEC = podgen.exe

.SUFFIXES :
.SUFFIXES : .f08 .c .o

.PHONY : default all clean

default : all

all :
	$(MAKE) $(EXEC)

$(EXEC) : $(SRC_OBJ)
	$(FC) $^ -o $@ $(NCDF) $(LAPACK)
	ls -fl $@

.c.o :
	$(CC) $(CFLAGS) $(CWARNINGS) $(CRUNTIMES) -c -I$(INC_DIR) $<

.f08.o :
	$(FC) $(FFLAGS) $(FWARNINGS) $(FRUNTIMES) -c -I$(INC_DIR) $<

clean :
	rm -f *.o

main.o :  NCDF_IO.o INPUTS.o POD_Routines.o GNUPLOT_ROUTINES.o CPOD_ROUTINES.o OUTPUTS.o
INPUTS.o : NCDF_IO.o
NCDF_IO.o :
GNUPLOT_ROUTINES.o :
POD_Routines.o :
CPOD_ROUTINES.o : MISC_PROCS.o
MISC_PROCS.o :
OUTPUTS.o :
