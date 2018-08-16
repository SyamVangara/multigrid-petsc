LIBFILES   =
IDIR	   = include
SDIR	   = src
BDIR	   = obj
TARGET    = poisson
_OBJ	   = poisson.o problem.o solver.o mesh.o array.o matbuild.o
OBJ	   = $(patsubst %,$(BDIR)/%,$(_OBJ))

CFLAGS     = -I$(IDIR) -std=c99 -fopenmp 
CPPFLAGS   =
LIBS	   = -lm
CLEANFILES = $(TARGET)

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

all: $(TARGET)

$(BDIR)/array.o: $(SDIR)/array.c $(IDIR)/array.h 
	-$(PETSC_COMPILE) -c -o $@ $< $(CFLAGS) 

$(BDIR)/matbuild.o: $(SDIR)/matbuild.c $(IDIR)/solver.h $(IDIR)/array.h $(IDIR)/mesh.h
	-$(PETSC_COMPILE) -c -o $@ $< $(CFLAGS) 

$(BDIR)/mesh.o: $(SDIR)/mesh.c  $(IDIR)/mesh.h $(IDIR)/array.h $(IDIR)/problem.h
	-$(PETSC_COMPILE) -c -o $@ $< $(CFLAGS) 

$(BDIR)/poisson.o: $(SDIR)/poisson.c  $(IDIR)/header.h $(IDIR)/array.h $(IDIR)/mesh.h $(IDIR)/solver.h
	-$(PETSC_COMPILE) -c -o $@ $< $(CFLAGS) 

$(BDIR)/problem.o: $(SDIR)/problem.c  $(IDIR)/problem.h $(IDIR)/array.h
	-$(PETSC_COMPILE) -c -o $@ $< $(CFLAGS) 

$(BDIR)/solver.o: $(SDIR)/solver.c  $(IDIR)/solver.h $(IDIR)/array.h $(IDIR)/mesh.h
	-$(PETSC_COMPILE) -c -o $@ $< $(CFLAGS) 

$(TARGET) : $(OBJ) 
	-${CLINKER} -o $(TARGET) $(OBJ) ${PETSC_KSP_LIB}
#	${RM} *.o

.PHONY: print_vars

#clean:
#	rm -f $(BDIR)/*.o *~ core $(INCDIR)/*~

print_vars:
	-@echo "PETSC_DIR: $(PETSC_DIR)"
	-@echo "CLINKER: $(CLINKER)"
	-@echo "CXX_COMPILE: $(PETSC_CXXCOMPILE)"
	-@echo "C_COMPILE: $(PETSC_COMPILE)"

tags : $(TARGET)
	ctags -R --c++-kinds=+p --fields=+iaS --extra=+q .

