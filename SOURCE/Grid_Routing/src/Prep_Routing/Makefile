# Makefile jroundy Jan 2012
FC = ifort
TFLAGS = -g -CB -check all -warn all -O0
PFLAGS = -CB -check all -O3 -openmp
FLAGS = $(PFLAGS)

OBJ = 	dist.o \
	area_m.o \
	point2grid.o \
	Prep_Routing.o
EXE = Prep_Routing.exe
# Makefile
$(EXE): $(OBJ)
	$(FC) -o $(EXE) -warn all $(FLAGS) $(OBJ) $(LIB)
%.o: %.f90
	$(FC) -c $(FLAGS) $<
# Cleaning everything
clean:
	rm -rf *.o *.mod *.exe
# End of the makefile
