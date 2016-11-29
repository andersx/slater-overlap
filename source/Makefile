FC = ifort 
FFLAGS = -O3 -xHost

COMPILE_MOD = 


all: main

main: ini_reader.o coordinates.o printer.o
	$(FC) $(FFLAGS) \
		main.f90 \
		-o ../bin/ppqm ini_reader.o coordinates.o printer.o

coordinates.o: coordinates.f90 printer.o
	$(FC) $(FFLAGS) -c \
		coordinates.f90 -o coordinates.o

printer.o: 
	$(FC) $(FFLAGS) -c \
	printer.f90 -o printer.o

ini_reader.o: 
	$(FC) $(FFLAGS) -c \
	ini_reader.f90 -o ini_reader.o

test:
	cd ../test && ./test_all


overlap_integrals.o: integrals/overlap_integrals.f90 overlap_backend.o
	$(FC) $(FFLAGS) -c \
	integrals/overlap_integrals.f90 -o overlap_integrals.o

overlap_backend.o: integrals/overlap_backend.f90
	$(FC) $(FFLAGS) -c \
	integrals/overlap_backend.f90 -o overlap_backend.o

clean:
	rm -f *.mod *.o
