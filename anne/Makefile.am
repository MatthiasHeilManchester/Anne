
# Names of executable
check_PROGRAMS=anne

#-----------------------------------------------------------------

# Sources for executable
anne_SOURCES = anne.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
anne_LDADD = -L@libdir@ -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


#-----------------------------------------------------------------

# Include path for library headers: All library headers live in 
# the include directory which we specify with -I
INCLUDES = -I@includedir@ 


# Files that aren't mentioned as sources above but are still supposed
# to be distributed
EXTRA_DIST = run.bash
