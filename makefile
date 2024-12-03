.SUFFIXES:

MAIN= kolmogorov_moran
PROGRAM = kolmogorov_moran
SRC=.f90
OBJDIR = obj
MODDIR = mod
FC = gfortran
FLAGS = -fmax-errors=3 -ffree-line-length-0 -J$(MODDIR)

all: $(PROGRAM)

$(PROGRAM): ./$(MAIN)$(SRC)
	$(FC) $(FLAGS) ./$(MAIN)$(SRC) -o $(PROGRAM)