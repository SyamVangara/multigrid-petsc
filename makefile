#  -g    adds debugging information to the executable file
#  -Wall turns on most, but not all, compiler warnings
#
# for C++ define  CC = g++
#include ${PETSC_DIR}/lib/petsc/conf/variables
#include ${PETSC_DIR}/lib/petsc/conf/rules
#
#CC = mpicc -O3
#CFLAGS  = -std=c99 -Wall
#
## External libraries:
#LIBS = -lm
#
## Pre-defined macros for conditional compilation
#DEFS = #-DDEBUG_FLAG -DEXPERIMENTAL=0
#
#BIN = poisson
#
#DEP = header.h 
#
#OBJS = solver.o mesh.o array.o matbuild.o
#
#$(BIN): $(BIN).o $(OBJS) $(DEP) chkopts
#	#$(CC) $(CFLAGS) $(DEFS) $(OBJS) $(BIN).c -o $(BIN) $(LIBS)
#	-${CLINKER} -o $(BIN) $(OBJS) $(BIN).o ${PETSC_LIB}
#	${RM} *.o
#
#%.o: %.c %.h
#	#$(CC) -c $(CFLAGS) $(DEFS) $< -o $@
#
##clean: 
##	$(RM) count *.o *~
#
#depend:
#	makedepend -Y -- $(CFLAGS) $(DEFS) -- *.c
#########################################################
#CFLAGS   = -O3 -std=c99 -Wall 
#FFLAGS   =
#CPPFLAGS =
#FPPFLAGS =
#
#include ${PETSC_DIR}/lib/petsc/conf/variables
#include ${PETSC_DIR}/lib/petsc/conf/rules
#
#BIN = poisson
#DEP = header.h 
#OBJS = solver.o mesh.o array.o matbuild.o
#
#CLEANFILES = $(BIN)
#
#
#$(BIN): $(BIN).o $(OBJS) $(DEP) chkopts
#	#$(CC) $(CFLAGS) $(DEFS) $(OBJS) $(BIN).c -o $(BIN) $(LIBS)
#	-${CLINKER} -o $(BIN) $(OBJS) $(BIN).o ${PETSC_LIB}
#	#${RM} *.o
#
#%.o: %.c %.h chkopts
#	#$(CC) -c $(CFLAGS) $(DEFS) $< -o $@
#	-${CLINKER} -c $< ${PETSC_LIB} -o $@
#
##########################################################

CFLAGS     = -O3 

CPPFLAGS   =

LIBFILES   =

TARGET     = poisson

OBJFILES   = poisson.o solver.o mesh.o array.o matbuild.o

CLEANFILES = $(TARGET)



include ${PETSC_DIR}/lib/petsc/conf/variables

include ${PETSC_DIR}/lib/petsc/conf/rules



all: $(TARGET)



$(TARGET) : $(OBJFILES)
	${CLINKER} -o $(TARGET) $(OBJFILES) ${PETSC_KSP_LIB}
	${RM} *.o
