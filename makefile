
CFLAGS     = -std=c99 -g #-O3 
CPPFLAGS   =
LIBFILES   =
TARGET     = poisson
OBJFILES   = poisson.o problem.o solver.o mesh.o array.o matbuild.o
CLEANFILES = $(TARGET)

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

all: $(TARGET)

$(TARGET) : $(OBJFILES)
	${CLINKER} -o $(TARGET) $(OBJFILES) ${PETSC_KSP_LIB}
	${RM} *.o
