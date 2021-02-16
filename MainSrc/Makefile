# Makefile

MOD_DIR = ./mods
OBJ_DIR = ./obj
MATH_DIR = /mnt/c/math_routines
SKIT_DIR = /mnt/c/SPARSKIT2
INC_DIR = /usr/include

#--------------------------------------------

FC = gfortran
FFLAGS = -O0 -fopenmp -I $(MATH_DIR) -I $(SKIT_DIR)
WARNINGS =  -Wall -Wextra -Wunused-dummy-argument -Wconversion -Wconversion-extra -Wintrinsic-shadow -Wsurprising
RUNTIMES = -fbacktrace -fbounds-check -fcheck=all
LAPACK = -llapack -lblas
NCDF = -lnetcdff

#--------------------------------------------

SRC_OBJ = main.o \
					INPUTS.o \
					OUTPUTS.o \
					TRANSPORT_SOLVES.o \
					INITIALIZERS.o \
					UPDATES.o \
					ALGORITHMS.o \
					CONVERGENCE_CHECKS.o \
					MLOQD_SOLVES.o \
					QD_SOLVERS.o \
					TEMP_FUNCTIONS.o \
					GLOQD_SOLVES.o \
					NCDF_IO.o \
					POD_ROUTINES.o \
					DISTDOT.o \
					GRID_FUNCTIONS.o

EXEC = 2dtrt.exe

.SUFFIXES :
.SUFFIXES : .f08 .o .f

.PHONY : default all clean

default : all

all :
	$(MAKE) $(EXEC)
	@mkdir -p $(MOD_DIR)
	mv *mod $(MOD_DIR)

$(EXEC) : $(SRC_OBJ)
	$(FC) $(FFLAGS) $(WARNINGS) $(RUNTIMES) $^ -o $@ -L$(MATH_DIR) -ljmcmath -lskit -L$(SKIT_DIR) -Wl,-rpath=$(MATH_DIR) $(LAPACK) $(NCDF)
	ls -fl $@

.f08.o :
	$(FC) $(FFLAGS) $(WARNINGS) $(RUNTIMES) -c -I$(MOD_DIR) -I$(INC_DIR) $<

clean :
	rm -f *.o
	rm -rf $(MOD_DIR)
	rm -rf $(OBJ_DIR)

main.o : INPUTS.o OUTPUTS.o ALGORITHMS.o
INPUTS.o :
OUTPUTS.o : NCDF_IO.o
TRANSPORT_SOLVES.o :
INITIALIZERS.o : TEMP_FUNCTIONS.o UPDATES.o MLOQD_SOLVES.o GLOQD_SOLVES.o
UPDATES.o : TEMP_FUNCTIONS.o
ALGORITHMS.o : TRANSPORT_SOLVES.o MLOQD_SOLVES.o GLOQD_SOLVES.o UPDATES.o CONVERGENCE_CHECKS.o OUTPUTS.o INITIALIZERS.o  POD_ROUTINES.o GRID_FUNCTIONS.o INPUTS.o
CONVERGENCE_CHECKS.o :
MLOQD_SOLVES.o : QD_SOLVERS.o
GLOQD_SOLVES.o : QD_SOLVERS.o
QD_SOLVERS.o : DISTDOT.o
NCDF_IO.o :
POD_ROUTINES.o : NCDF_IO.o  GRID_FUNCTIONS.o
DISTDOT.o :
GRID_FUNCTIONS.o :
