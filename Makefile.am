#Name of executable
check_PROGRAMS=muenster_lubri_beam

#Sources for the executable
muenster_lubri_beam_SOURCES = muenster_lubri_beam.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
muenster_lubri_beam_LDADD = -L@libdir@ -lbeam -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


# Include path for library headers: All library headers live in 
# the include directory which we specify with -I
# Automake will replace the variable @includedir@ with the actual
# include directory.
AM_CPPFLAGS += -I@includedir@  
