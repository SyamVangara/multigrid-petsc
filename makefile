#  -g    adds debugging information to the executable file
#  -Wall turns on most, but not all, compiler warnings
#
# for C++ define  CC = g++
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

CC = mpicc #-O
CFLAGS  = -std=c99 -Wall

# External libraries:
LIBS = -lm

# Pre-defined macros for conditional compilation
DEFS = #-DDEBUG_FLAG -DEXPERIMENTAL=0

BIN = poisson

DEP = header.h 

OBJS = solver.o mesh.o array.o

$(BIN): $(BIN).o $(OBJS) $(DEP) chkopts
	#$(CC) $(CFLAGS) $(DEFS) $(OBJS) $(BIN).c -o $(BIN) $(LIBS)
	-${CLINKER} -o $(BIN) $(OBJS) $(BIN).o ${PETSC_LIB}

%.o: %.c %.h
	$(CC) -c $(CFLAGS) $(DEFS) $< -o $@

#clean: 
#	$(RM) count *.o *~

depend:
	makedepend -Y -- $(CFLAGS) $(DEFS) -- *.c
