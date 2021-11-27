UNAMES = $(shell uname -s)

EXEC = dipole_coordinates.exe

FORTRAN = gfortran -std=legacy
FORTRAN = gfortran

FFLAGS = -g -fcheck=all
FFLAGS = -O2

OBJS = mod_precision.o mod_const.o dipole_mu.o dipole_psi.o \
		 dipole2spherical.o dipole_coordinates.o mod_igrf.o

MODS = mod_precision.mod mod_const.mod mod_igrf.mod

.PHONY: dipole clean 

$(EXEC): $(MODS) $(OBJS)
	$(FORTRAN) $(FFLAGS) $(OBJS) -o $@

%.o : %.f90
	$(FORTRAN) $(FFLAGS) -c $< 

%.mod: %.f90
	$(FORTRAN) $(FFLAGS) -c $<

clean:
	rm -f *.o *.exe *.mod

