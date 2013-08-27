# Makefile jroundy Jan 2012
FC = /opt/intel/bin/ifort
TFLAGS = -g -CB -check all -warn all -O0
PFLAGS = -CB -check all -O3 -openmp
FLAGS = $(PFLAGS)

OBJ = 	isleap.o \
	daysMonth.o \
	zeroMonth.o \
	Rout_Grid.o \
	Grid_Routing.o
EXE = Grid_Routing.exe
# Makefile
$(EXE): $(OBJ)
	$(FC) -o $(EXE) -warn all $(FLAGS) $(OBJ)
%.o: %.f90
	$(FC) -c $(FLAGS) $<
# Cleaning everything
clean:
	rm -rf *.o *.mod *.exe
# End of the makefile
