UNAMES = $(shell uname -s)

EXEC = general_coordinates.exe

FORTRAN = gfortran -std=legacy
FORTRAN = gfortran

LFLAGS =
LIBS =

FFLAGS = -g -fcheck=all
FFLAGS = -O2
 
OBJS = mod_precision.o mod_const.o PlmSchmidt_d1.o PlmIndex.o \
		 mod_igrf.o rksuite_90.o mod_function.o zero_rc.o \
		 mod_fieldline.o general_coordinates.o \

MODS = mod_precision.mod mod_const.mod mod_igrf.mod rksuite_90.mod \
		 mod_function.mod mod_fieldline.mod

$(EXEC):	$(MODS) $(OBJS)
	$(FORTRAN) -o $@ $(OBJS) $(LFLAGS) $(LIBS)
 
%.o %.mod: %.f90
	$(FORTRAN) -c $(FFLAGS) $<

%.o %.mod: %.f95
	$(FORTRAN) -c $(FFLAGS) $<

.f90.o:
	$(FORTRAN) -c $(FFLAGS) $<

.f95.o:
	$(FORTRAN) -c $(FFLAGS) $<

clean:
	rm -f *.o *.mod *.exe

