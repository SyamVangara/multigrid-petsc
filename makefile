include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

main: main.o chkopts
	-${CLINKER} -o main main.o ${PETSC_LIB}
	${RM} main.o

