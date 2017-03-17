################################################################################
# Makefile for building a HepMC executable with the gcc compiler
#
# Matt.Dobbs@CERN.CH  1.2000
#
# This makefiles also works to compile the example programs.
# I.E.: syntax for compiling the example_MyPythia.cc example:
#       gmake example_MyPythia.exe
# or simply   gmake all     to compile all examples.
#
################################################################################ Define directory paths 
#    You may have to change GENSERdir and/or other variables
#
  HepMCdir             = /usr/lib64
  HepMClib             = -L$(HepMCdir)/lib -lHepMC
  HepMCfiolib          = -L$(HepMCdir)/lib -lHepMCfio
  GENSERdir            = 
  CLHEPdir             = 
  FASTJETCONFIG        := /home/rbertens/Documents/CERN/alice/JET_SOFTWARE/fastjet-3.0.6_BUILD/bin/fastjet-config

################################################################################ Compiler options
#
  CXX           = g++
  F77		= gfortran
  INCLUDES 	= -I/usr/include -I$(CLHEPdir)/include
#  CXXFLAGS      =   -O -ansi -pedantic -Wall -D_GNU_SOURCE -O2 -g $(INCLUDES)
  CXXFLAGS      =   -O -D_GNU_SOURCE -O2 -g $(INCLUDES)
  FLAGS 	= $(DFLG) $(INCDIR)
ifeq "$(CXX)" "g++"
   FLAGS 	+= -fno-second-underscore
endif

  EXAMPLES	= example_BuildEventFromScratch.exe	\
		  example_EventSelection.exe	\
		  example_UsingIterators.exe
  LINK_LIBS     =  
  UNAME = $(shell uname)

ifneq (,$(strip $(CLHEPdir)))
  EXAMPLES	+= example_VectorConversion.exe
ifeq "$(UNAME)" "Darwin"
  CLHEP_LIB     = -L$(CLHEPdir)/lib -lCLHEP 
else
  CLHEP_LIB     = -L$(CLHEPdir)/lib -lCLHEP -Wl,-rpath -Wl,$(CLHEPdir)/lib
endif  
endif

ifeq "$(UNAME)" "Darwin"
else
  LINK_LIBS     += -Wl,-rpath -Wl,$(HepMCdir)/lib
endif

  HDRS          = $(HepMCdir)/include/HepMC/*.h *.h


FASTJETLIBS := $(shell $(FASTJETCONFIG) --libs)
FASTJETFLAGS := $(shell $(FASTJETCONFIG) --cxxflags)

ROOTLIBS := $(shell root-config  --libs)
ROOTFLAGS := $(shell root-config --cflags)

CXXFLAGS += $(FASTJETFLAGS) $(ROOTFLAGS)
FLAGS += $(FASTJETLIBS) $(ROOTLIBS)

################################################################################ 

PLATFORM=$(shell uname)

ifeq "$(PLATFORM)" "Linux"
    LINK_LIBS     += -lgfortran
endif

################################################################################ definition of the compiler options
#	-I location of directory containing include files
#	-L location of directory containing libraries
#       -lname include the library from -L location called libname.a
#	   -lg2c is the library containing info on converting fortran to C
#          -lf   is the library containing the intrinsic for HPUX only.
#	-shared make a shared library as output
#	-fPIC produce position independent code
#        necessary on some platforms (including HPUX) for -shared
# 	-fpic ^^ same(?)
#	-O optimizes
#	-g produces output for the debugger
#       -pg produces output for gprof profiler
#       note: if you want to see all warnings and ensure ansi standard 
#             compatibility, use:
#             -pipe -ansi -pedantic -fnonnull-objects \
#             -W -Wall -Wwrite-strings -Wpointer-arith -Wnested-externs \
#             -Woverloaded-virtual -Wbad-function-cast -fnonnull-objects
#       The proper order for cernlib libraries is:
#       -lpawlib -lgraflib -lgrafX11 -lmathlib -lkernlib -lpacklib -ljetset74
#
# makefile syntax:
#        for target thedir/target.suf from source anotherdir/source.suf2
#        ${*D}  = thedir
#        ${*F}  = target
#        $*     = thedir/target
#        $@     = thedir/target.suf
#        $<     = anotherdir/source.suf2
#  

###############################################################################
#
.SUFFIXES:      .o .cxx .exe
all:	analyze_hepmc_jets

analyze_hepmc_jets: analyze_hepmc_jets.o
	@echo "Building $@ ..."
	$(CXX) $(FLAGS) analyze_hepmc_jets.o \
		$(HepMClib) \
	        $(LINK_LIBS) -o $@

analyze_hepmc_jet_frag: analyze_hepmc_jet_frag.o
	@echo "Building $@ ..."
	$(CXX) $(FLAGS) analyze_hepmc_jet_frag.o \
		$(HepMClib) \
	        $(LINK_LIBS) -o $@

analyze_hepmc_jet_shapes: analyze_hepmc_jet_shapes.o ExampleShapes.o
	@echo "Building $@ ..."
	$(CXX) $(FLAGS) $^ \
		$(HepMClib) \
	        $(LINK_LIBS) -lfastjetcontribfragile -o $@

analyze_hepmc_jet_shapes_constsub: analyze_hepmc_jet_shapes_constsub.o ExampleShapes.o
	@echo "Building $@ ..."
	$(CXX) $(FLAGS) $^ \
		$(HepMClib) \
	        $(LINK_LIBS) -lfastjetcontribfragile -o $@

analyze_hepmc_hjet: analyze_hepmc_hjet.o
	@echo "Building $@ ..."
	$(CXX) $(FLAGS) analyze_hepmc_hjet.o \
		$(HepMClib) \
	        $(LINK_LIBS) -o $@

analyze_hepmc_dihadrons: analyze_hepmc_dihadrons.o
	@echo "Building $@ ..."
	$(CXX) $(FLAGS) analyze_hepmc_dihadrons.o \
		$(HepMClib) \
	        $(LINK_LIBS) -o $@

analyze_hepmc_prim_sec: analyze_hepmc_prim_sec.o
	@echo "Building $@ ..."
	$(CXX) $(FLAGS) analyze_hepmc_prim_sec.o \
		$(HepMClib) \
	        $(LINK_LIBS) -o $@
###############################################################################
# instructions for building a .o file from a .cc file
#
.cxx.o:         $(HDRS) $<
	@echo "Compiling $< with $(CXX) ..."
	@$(CXX) $(CXXFLAGS) -c $< -o $@

###############################################################################
# gmake clean       removes all garbage from HepMC directories.
#
clean:
	rm -f *.o

###############################################################################
# gmake distclean       removes all compiled libraries, executables, +garbage
#                       to prepare the package for distribution
distclean: 
	$(MAKE) clean --no-print-directory
	rm -f *.exe
	rm -f *.dat

