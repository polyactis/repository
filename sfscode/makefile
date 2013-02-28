# Set installation path and src/util paths
EXEDIR = bin
UTILDIR = src/util

# Set optimization flags
# for gcc -O3 should work
# on MacOSX can use -fast instead
#  to debug, use -g then run using gdb
#  to run with profiler gprof, use -g -pg
CFLAGS = -Wall -O3 -funroll-loops

# Set the compiler
CC = gcc $(CFLAGS)

# Link the following libraries
LIBS = -lm 
GSLLIBS = -lgsl -lgslcblas

#list all programs to be installed by make
PROGS = sfs_code convertSFS_CODE 

all: introduce $(PROGS) install
	@echo "Finished compiling"
	@echo ""

introduce:
	@echo "Building SFS_CODE program files"


install: 
	@echo "Installing SFS_CODE binaries in $(EXEDIR)"
	@mkdir -p $(EXEDIR)
	@mv $(PROGS) $(EXEDIR)
	@echo "Finished installation."
	@echo ""

clean: 
	@echo "Removing object and executable files"
	@rm -f */*.o
	@rm -f $(UTILDIR)/*.o
	@rm -f $(EXEDIR)/*
	@echo "Finished cleanup."
	@echo ""

SF = src/sfs_code_DIR
sfs_code : $(SF)/sfs_code.o $(SF)/btfuncs.c $(SF)/gencontextrate.h $(SF)/sfs_code.h
	$(CC) $(SF)/sfs_code.o $(LIBS) -o sfs_code
$(SF)/sfs_code.o : $(SF)/btfuncs.c $(SF)/gencontextrate.h $(SF)/sfs_code.h $(UTILDIR)/nrutil.c $(UTILDIR)/nrutil.h $(UTILDIR)/myrand1.c makefile

convertSFS_CODE : $(SF)/convertSFS_CODE.o $(SF)/sfs_code.h
	$(CC) $(SF)/convertSFS_CODE.o $(LIBS) -o convertSFS_CODE
$(SF)/convertSFS_CODE.o : $(SF)/sfs_code.h $(UTILDIR)/nrutil.c $(UTILDIR)/nrutil.h $(UTILDIR)/myrand1.c
