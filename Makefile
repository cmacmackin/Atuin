# The compiler
FC = gfortran
CC = gcc
# flags for debugging or for maximum performance, comment as necessary
FCFLAGS = -g -fbounds-check -ffpe-trap=i,z,o
#FCFLAGS = -O3
# flags forall (e.g. look for system .mod files, required in gfortran)
#FCFLAGS += -I/usr/include

# libraries needed for linking, unused in the examples
#LDFLAGS = -li_need_this_lib

# List of executables to be built within the package
PROGRAMS = randgen.o physvals.o particle.o pdstrbtn.o forces.o tree2d.o leap2d.o 

# "make" builds all
all: $(PROGRAMS) atuin.o
	 $(FC) -o atuin $(PROGRAMS) atuin.o

%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.f
	$(FC) $(FCFLAGS) -c $<

%.o: %.F
	$(FC) $(FCFLAGS) -c $<

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

%.o: %.F90
	$(FC) $(FCFLAGS) -c $<

%.o: %.f95
	$(FC) $(FCFLAGS) -c $<

%.o: %.F95
	$(FC) $(FCFLAGS) -c $<

posixwrapper.o : posixwrapper.c
				 $(CC) -c posixwrapper.c
	
video: posixwrapper.o fortranposix.o gnuplot_fortran95.o physvals.o randgen.o particle.o pdstrbtn.o mkvideo.o
	   $(FC) -o mkvideo posixwrapper.o fortranposix.o gnuplot_fortran95.o physvals.o randgen.o particle.o pdstrbtn.o mkvideo.o

hist: posixwrapper.o fortranposix.o gnuplot_fortran95.o physvals.o randgen.o particle.o pdstrbtn.o mkhist.o
	  $(FC) -o mkhist posixwrapper.o fortranposix.o gnuplot_fortran95.o physvals.o randgen.o particle.o pdstrbtn.o mkhist.o

plothist: posixwrapper.o fortranposix.o gnuplot_fortran95.o physvals.o plothist.o
	  $(FC) -o plothist posixwrapper.o fortranposix.o gnuplot_fortran95.o physvals.o plothist.o

histdat: physvals.o randgen.o particle.o pdstrbtn.o mkhistdat.o
	     $(FC) -o histdat physvals.o randgen.o particle.o pdstrbtn.o mkhistdat.o

clean:
	/bin/rm  -rf *.o ./tmp *.mod
