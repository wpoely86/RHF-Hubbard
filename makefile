###############################################################################
#
#  makefile template for the sources
#
###############################################################################

# -----------------------------------------------------------------------------
#   Sources for all modules
# -----------------------------------------------------------------------------
BINNAME = HF
CPPSRC	= HF.cpp\
            Matrix.cpp\
            Vector.cpp\
            LibInt.cpp\
            R.cpp\
            Gauss.cpp\
            input.cpp\
            CI_SPM.cpp\
            CI_SPPM.cpp\
            CI_SPPM_m.cpp\
            CI_TPM.cpp\
            CI_TPPM.cpp\
            CartInt.cpp\
            Transform.cpp\
            SI_SPM.cpp\
            SI_TPM.cpp\
            SphInt.cpp\
            Tools.cpp\
            RSPM.cpp\
            DIIS.cpp

OBJ	= $(CPPSRC:.cpp=.o)

# -----------------------------------------------------------------------------
#   These are the standard libraries, include paths and compiler settings
# -----------------------------------------------------------------------------

BRIGHT_ROOT= .

INCLUDE = ./include -I/usr/local/libint/2.0.0-stable/include -I/usr/local/libint/2.0.0-stable/include/libint2 

LIBS= -L/usr/local/libint/2.0.0-stable/lib -llapack -lblas -lgmp -lgmpxx -lgsl -lhdf5 -lint2

CC	= gcc
CXX	= g++

# -----------------------------------------------------------------------------
#   Compiler & Linker flags
# -----------------------------------------------------------------------------
CFLAGS	= -I$(INCLUDE) -g -Wall
LDFLAGS	= -g -Wall


# =============================================================================
#   Targets & Rules
# =============================================================================
all:
	@echo
	@echo '  +++ Building $(BINNAME)...'
	@echo	
	$(MAKE) $(BRIGHT_ROOT)/$(BINNAME) 
	@if test $?; then \
	   echo; echo '*************** FAILED! ***************' ; echo; \
	 else \
	   echo; echo '  +++ $(BINNAME) has been built successfully!'; \
	   echo; \
	 fi

# -----------------------------------------------------------------------------
#   The default way to compile all source modules
# -----------------------------------------------------------------------------
%.o:	%.for makefile
	@echo; echo "Compiling $(@:.o=.for) ..."
	$(FF) -c $(FFLAGS) $(SFLAGS) $(@:.o=.for) -o $@

%.o:	%.c makefile
	@echo; echo "Compiling $(@:.o=.c) ..."
	$(CC) -c $(CFLAGS) $(SFLAGS) $(@:.o=.c) -o $@

%.o:	%.cpp makefile
	@echo; echo "Compiling $(@:.o=.cpp) ..."
	$(CXX) -c $(CFLAGS) $(SFLAGS) $(DEFS) $(@:.o=.cpp) -o $@


# -----------------------------------------------------------------------------
#   Link everything together
# -----------------------------------------------------------------------------
$(BRIGHT_ROOT)/$(BINNAME):	makefile $(OBJ) 
	@echo; echo "Linker: creating $(BRIGHT_ROOT)/$(BINNAME) ..."
	$(CXX) $(LDFLAGS) $(SFLAGS) -o $(BRIGHT_ROOT)/$(BINNAME) $(OBJ) $(LIBS)

# -----------------------------------------------------------------------------
#   Create everything newly from scratch
# -----------------------------------------------------------------------------
new:	clean all

# -----------------------------------------------------------------------------
#   Clean up all object files
# -----------------------------------------------------------------------------
clean:
	@echo -n '  +++ Cleaning all object files ... '
	@echo -n $(OBJ)
	@rm -f $(OBJ) read.o
	@echo 'Done.'

# -----------------------------------------------------------------------------
#   Make new documentation using doxygen
# -----------------------------------------------------------------------------
doc:
	@doxygen doc-config

read: read.cpp $(CPPSRC:bp_sdp.cpp=)
	$(MAKE) read.o $(OBJ:bp_sdp.o=) DEFS="-DPQG"
	$(CXX) $(LDFLAGS) $(SFLAGS) -o read read.o $(OBJ:bp_sdp.o=) $(LIBS)

# ====================== End of file 'makefile.in' ========================== #
