# Makefile

INC_DIR = /usr/include

#--------------------------------------------

CC = gcc
CFLAGS = -O0
CWARNINGS =  -Wall -Wextra -Wconversion
CRUNTIMES = -fbounds-check #-fsanitize=bounds

FC = gfortran
FFLAGS = -O0
FWARNINGS =  -Wall -Wextra -Wunused-dummy-argument -Wconversion -Wconversion-extra -Wintrinsic-shadow -Wsurprising
FRUNTIMES = -fbacktrace -fbounds-check -fcheck=all

NCDF = -lnetcdf

#--------------------------------------------

SRC_OBJ = main.o \
          POD_Routines.o

EXEC = podgen.exe

.SUFFIXES :
.SUFFIXES : .f08 .c .o

.PHONY : default all clean

default : all

all :
	$(MAKE) $(EXEC)

$(EXEC) : $(SRC_OBJ)
	$(FC) $^ -o $@ $(NCDF)
	ls -fl $@

.c.o :
	$(CC) $(CFLAGS) $(CWARNINGS) $(CRUNTIMES) -c -I$(INC_DIR) $<

.f08.o :
	$(FC) $(FFLAGS) $(FWARNINGS) $(FRUNTIMES) -c -I$(INC_DIR) $<

clean :
	rm -f *.o

main.o : POD_Routines.o
POD_Routines.o :
