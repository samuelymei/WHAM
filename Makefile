.SUFFIXES: .f90

FC = ifort
FFLAGS = -CB
INCLUDE = 
LIBS = 

EXE = WHAM.x

OBJS = precision_m.o constant.o snapshot.o bin.o react_coord_bin.o simulation.o WHAM.o WHAM_caller.o

all:	${EXE}


$(EXE):$(OBJS)
	$(FC) -o $@ $(FFLAGS) $(OBJS) $(LIBS)

%.o %mod: %.f90
	$(FC) -c $(FFLAGS) $(INCLUDE) $<

clean:
	/bin/rm -f $(EXE) $(OBJS) *.mod
