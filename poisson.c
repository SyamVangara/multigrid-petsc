#include "header.h"

#define ERROR_MSG(message) (fprintf(stderr,"Error:%s:%d: %s\n",__FILE__,__LINE__,(message)))
#define ERROR_RETURN(message) {ERROR_MSG(message);return ierr;}
#define CHKERR_PRNT(message) {if(ierr==1) {ERROR_MSG(message);}}
#define CHKERR_RETURN(message) {if(ierr==1) {ERROR_RETURN(message);}}

#define PI 3.14159265358979323846
#define DIMENSION 2
#define FUNC(i,j) (-2*PI*PI*sin(PI*coord[0][(j)])*sin(PI*coord[1][(i)]))
#define SOL(i,j) (sin(PI*coord[0][(j)])*sin(PI*coord[1][(i)]))

static int ipow(int base, int exp);
int totalUnknowns(int *n, int totalGrids);

void PrintInfo(Problem prob, Mesh mesh, Indices indices, Operator op, Solver solver, PostProcess pp, int cyc, int meshflag, int mappingStyleflag);
void ViewMeshInfo(Mesh mesh);
void ViewGridsInfo(Indices indices);
void ViewIndexMapsInfoLevel(Level level, int l);
void ViewIndexMapsInfo(Indices indices);
void ViewRangesInfo(Indices indices);
void ViewSolverInfo(Indices indices, Solver solver);
void ViewOperatorInfo(Operator op);
void ViewLinSysMatsInfo(Assembly assem, int view);
void ViewGridTransferMatsInfo(Assembly assem, int view, int cyc);

int main(int argc, char *argv[]) {

	PetscInitialize(&argc, &argv, "poisson.in", 0);

	Problem		prob;
	Mesh		mesh;
	Indices		indices;
	Operator	op;
	Solver		solver;
	PostProcess	pp;	
	
	int		cyc;
	int		meshflag;
	int		mappingStyleflag;
	int		vmax = 2;

	int	ierr=0;
	
	SetUpProblem(&prob);
	
//	freopen("poisson.in", "r", stdin);
//	freopen("poisson.out", "w", stdout);
//	freopen("poisson.err", "w", stderr);
	
	PetscOptionsGetInt(NULL, NULL, "-npts", mesh.n, NULL);
	PetscOptionsGetInt(NULL, NULL, "-mesh", &meshflag, NULL);
	PetscOptionsGetInt(NULL, NULL, "-iter", &(solver.numIter), NULL);
	PetscOptionsGetInt(NULL, NULL, "-grids", &(indices.totalGrids), NULL);
	PetscOptionsGetInt(NULL, NULL, "-levels", &(indices.levels), NULL);
	PetscOptionsGetInt(NULL, NULL, "-cycle", &(cyc), NULL);
	PetscOptionsGetInt(NULL, NULL, "-map", &(mappingStyleflag), NULL);
	PetscOptionsGetIntArray(NULL, NULL, "-v", solver.v, &vmax, NULL);
	PetscOptionsGetInt(NULL, NULL, "-moreNorm", &(solver.moreInfo), NULL);
	
	if (indices.levels>1 && cyc==3) {
		PetscPrintf(PETSC_COMM_WORLD, "For now only one level is allowed for delayed cycling"); 
		PetscFinalize();
		return 0;
	}

	for (int i=1;i<DIMENSION;i++) { 
		mesh.n[i]  = mesh.n[0];      // No. of points in each dimension
	}
	for (int i=0;i<DIMENSION;i++) {
		mesh.bounds[i*2] = 0.0;    // Lower bound in each dimension
		mesh.bounds[i*2+1] = 1.0;  // Upper bound in each dimension
	}
	

	if (meshflag == 0) SetUpMesh(&mesh, UNIFORM);
	if (meshflag == 1) SetUpMesh(&mesh, NONUNIFORM);

//	ViewMeshInfo(mesh);
	
	indices.coarseningFactor = 2;
	SetUpIndices(&mesh, &indices);

//	ViewGridsInfo(indices);

	mapping(&indices, mappingStyleflag);
	
//	ViewIndexMapsInfo(indices);
//	ViewRangesInfo(indices);
	
	SetUpOperator(&indices, &op);
	GridTransferOperators(op, indices);

//	ViewOperatorInfo(op);
	
	if (cyc == 0) SetUpSolver(&indices, &solver, VCYCLE);
	if (cyc == 1) SetUpSolver(&indices, &solver, ICYCLE);
	if (cyc == 2) SetUpSolver(&indices, &solver, ECYCLE);
	if (cyc == 3) SetUpSolver(&indices, &solver, D1CYCLE);

//	ViewSolverInfo(indices, solver);

	Assemble(&prob, &mesh, &indices, &op, &solver);

//	ViewLinSysMatsInfo(*(solver.assem), 0);
//	ViewGridTransferMatsInfo(*(solver.assem), 0, cyc);
	
	Solve(&solver);
	
	SetUpPostProcess(&pp);
	Postprocessing(&prob, &mesh, &indices, &solver, &pp);
	
	PrintInfo(prob, mesh, indices, op, solver, pp, cyc, meshflag, mappingStyleflag);
	
	DestroyPostProcess(&pp);
	DestroySolver(&solver);
	DestroyOperator(&op);
	DestroyIndices(&indices);
	DestroyMesh(&mesh);
	PetscFinalize();

	return 0;
}

int ipow(int base, int exp) {

	int result = 1;
	while (exp) {
		if (exp & 1)
			result *= base;
		exp >>= 1;
		base *= base;
	}
	return result;
}

int totalUnknowns(int *n, int totalGrids) {
		
	int length, n0;

	n0 = n[0]-2;
	length=n0*n0;
	for (int i=1;i<totalGrids;i++) {
		n0 = (n0-1)/2;
		length = length + n0*n0;
	}
	return length;
}

void PrintInfo(Problem prob, Mesh mesh, Indices indices, Operator op, Solver solver, PostProcess pp, int cyc, int meshflag, int mappingStyleflag) {
	// Prints complete some info of problem, grids, solver
	
	int	procs, rank;

	MPI_Comm_size(PETSC_COMM_WORLD, &procs);
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	
	if (rank==0) {
	printf("=============================================================\n");
	printf("Size:				%d x %d\n", mesh.n[0], mesh.n[1]);
	if (meshflag==0) printf("Mesh Type:			Uniform\n");
	if (meshflag==1) printf("Mesh Type:			Non Uniform\n");
	printf("Number of grids:		%d\n",op.totalGrids);
	printf("Number of levels:		%d\n",solver.assem->levels);
	printf("Number of grids per level:	");
	for (int l=0;l<indices.levels;l++) {
		printf("%d	", indices.level[l].grids);
	}
	printf("\n");
	printf("Number of unknowns per level:	");
	for (int l=0;l<indices.levels;l++) {
		printf("%d	", indices.level[l].global.ni);
	}
	printf("\n");
	if (mappingStyleflag == 0) printf("Mapping style :			Grid after grid\n");
	if (mappingStyleflag == 1) printf("Mapping style :			Through the grids\n");
	if (mappingStyleflag == 2) printf("Mapping style :			Local grid after grid\n");
	if (cyc == 3) printf("Cycle :				D1-Cycle\n");
	if (cyc == 2) printf("Cycle :				E-Cycle\n");
	if (cyc == 1) printf("Cycle :				I-Cycle\n");
	if (cyc == 0) printf("Cycle :				V-Cycle\n");
	
	if (cyc == 3) printf("Number of smoothing steps :	%d per iteration \n", solver.v[0]);
	if (cyc == 2) printf("Number of smoothing steps :	%d per RHS update \n", solver.v[0]);
	if (cyc == 0) printf("Number of smoothing steps :	%d(fine) %d(coarsest)\n", solver.v[0], solver.v[1]);
	printf("Number of processes:		%d\n",procs);
	printf("Number of iterations:		%d\n",solver.numIter);
	printf("=============================================================\n");
	}

}

void ViewMeshInfo(Mesh mesh) {
	// Prints the info in mesh data structure

	int	procs, rank;
	
	MPI_Comm_size(PETSC_COMM_WORLD, &procs);
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	
	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = %d: n[0] = %d, n[1] = %d\n", rank, mesh.n[0], mesh.n[1]);
	PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
	
	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = %d: bounds[0:1] = %f, %f; bounds[2:3] = %f, %f\n", rank, mesh.bounds[0], mesh.bounds[1], mesh.bounds[2], mesh.bounds[3]);
	PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
	
	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = %d: h = %f\n", rank, mesh.h);
	PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
      
	for (int i=0;i<DIMENSION;i++) {
		PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = %d: coord[%d]:",rank,i);
		for (int j=0;j<mesh.n[i];j++) {
			PetscSynchronizedPrintf(PETSC_COMM_WORLD," %f ",mesh.coord[i][j]);
		}
		PetscSynchronizedPrintf(PETSC_COMM_WORLD,"\n");
	}
	PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
}

void ViewGridsInfo(Indices indices) {
	// Prints the info of GridId in each level
	
	int	procs, rank;
	
	MPI_Comm_size(PETSC_COMM_WORLD, &procs);
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

	for (int l=0;l<indices.levels;l++) {
		for (int lg=0;lg<indices.level[l].grids;lg++) {
			PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = %d; level: %d; gridId: %d; (ni, hi): (%d, %f); (nj, hj): (%d, %f)\n",rank,l,indices.level[l].gridId[lg], indices.level[l].grid[lg].ni, indices.level[l].h[lg][0], indices.level[l].grid[lg].nj, indices.level[l].h[lg][1]);
		}
	}
	PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
}

void ViewIndexMapsInfoLevel(Level level, int l) {
/********************************************************************************
 *
 * Print the Index maps of the given level
 *
 ********************************************************************************/ 
	int	procs, rank;
	
	MPI_Comm_size(PETSC_COMM_WORLD, &procs);
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	
	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = %d; level: %d; ranges: ", rank, l);
	for (int p=0;p<procs+1;p++) {
		PetscSynchronizedPrintf(PETSC_COMM_WORLD,"%d ",level.ranges[p]);
	}
	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"\n");
	PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
	
	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = %d: level = %d: Map from grid indices to global index:\n",rank,l);
	for (int lg=0; lg<level.grids; lg++) {
		for (int i=0; i<level.grid[lg].ni; i++) {
			for (int j=0; j<level.grid[lg].nj; j++) {
				PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = %d: level = %d: level_gridId = %d: gridId = %d: (row = %d, col = %d) -> (global_index = %d)\n",rank,l,lg,level.gridId[lg],i,j,level.grid[lg].data[i*level.grid[lg].nj+j]);
			}
		}
	}
	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"\n");
	PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
	
	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = %d: level = %d: Map from global index to grid indices:\n",rank,l);
	for (int i=0;i<level.global.ni;i++) {
		PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = %d: level = %d: (global_index = %d) -> (row = %d, col = %d, gridId = %d)\n",rank,l,i,level.global.data[i*level.global.nj+0], level.global.data[i*level.global.nj+1], level.global.data[i*level.global.nj+2]);
	}
	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"\n");
	PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
}

void ViewIndexMapsInfo(Indices indices) {
	// Prints the info of index maps between global and grids
	
	int	procs, rank;
	
	MPI_Comm_size(PETSC_COMM_WORLD, &procs);
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	
	for (int l=0;l<indices.levels;l++) {
		ViewIndexMapsInfoLevel(indices.level[l], l);
	}
}

void ViewRangesInfo(Indices indices) {
	// Prints the info of index maps between global and grids
	
	int	procs, rank;
	
	MPI_Comm_size(PETSC_COMM_WORLD, &procs);
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	
	for (int l=0;l<indices.levels;l++) {
		PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = %d; level: %d; ranges: ", rank, l);
		for (int p=0;p<procs+1;p++) {
			PetscSynchronizedPrintf(PETSC_COMM_WORLD,"%d ",indices.level[l].ranges[p]);
		}
		PetscSynchronizedPrintf(PETSC_COMM_WORLD,"\n");
	}
	PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
}

void ViewSolverInfo(Indices indices, Solver solver) {
	// Prints the info in Solver struct
	
	int	procs, rank;
	
	MPI_Comm_size(PETSC_COMM_WORLD, &procs);
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	
	if (solver.cycle==VCYCLE) {
		PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = %d; numIter = %d; Cycle = VCYCLE\n",rank,solver.numIter);
	} else {
		PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = %d; numIter = %d; Cycle = ICYCLE\n",rank,solver.numIter);
	}
	PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
}

void ViewOperatorInfo(Operator op) {
	// Prints the info in Operator struct
	
	int	procs, rank;
	
	MPI_Comm_size(PETSC_COMM_WORLD, &procs);
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	
	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = %d; Total num of grids = %d:\n",rank,op.totalGrids);
	PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
	for (int l=0;l<op.totalGrids-1;l++) {
		PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = %d; res[%d]:\n",rank,l);
		for (int i=0;i<op.res[l].ni;i++) {
			for (int j=0;j<op.res[l].nj;j++) {
				PetscSynchronizedPrintf(PETSC_COMM_WORLD,"%f ",op.res[l].data[i*op.res[l].nj+j]);
			}
			PetscSynchronizedPrintf(PETSC_COMM_WORLD,"\n");
		}
		PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = %d; pro[%d]:\n",rank,l);
		for (int i=0;i<op.pro[l].ni;i++) {
			for (int j=0;j<op.pro[l].nj;j++) {
				PetscSynchronizedPrintf(PETSC_COMM_WORLD,"%f ",op.pro[l].data[i*op.pro[l].nj+j]);
			}
			PetscSynchronizedPrintf(PETSC_COMM_WORLD,"\n");
		}
	}
	PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
}

void ViewLinSysMatsInfo(Assembly assem, int view) {
	// Prints the info of Assembly struct
	
	int	procs, rank;
	
	MPI_Comm_size(PETSC_COMM_WORLD, &procs);
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	
	for (int l=0;l<assem.levels;l++) {
		PetscPrintf(PETSC_COMM_WORLD,"A[%d]:\n",l);
		if (view == 0) {
			MatView(assem.A[l],PETSC_VIEWER_STDOUT_WORLD);
			MatView(assem.A2[l],PETSC_VIEWER_STDOUT_WORLD);
		}
		if (view == 1) {
			MatView(assem.A[l],PETSC_VIEWER_DRAW_WORLD);
			MatView(assem.A2[l],PETSC_VIEWER_DRAW_WORLD);
		}
		VecView(assem.b[l],PETSC_VIEWER_STDOUT_WORLD);
	}
	PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
}

void ViewGridTransferMatsInfo(Assembly assem, int view, int cyc) {
	// Prints the info of Assembly struct
	
	int	procs, rank;
	
	MPI_Comm_size(PETSC_COMM_WORLD, &procs);
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

	if (cyc == 3)
	{
		PetscPrintf(PETSC_COMM_WORLD,"\nres:\n");
		if (view == 0) MatView(assem.res[0],PETSC_VIEWER_STDOUT_WORLD);
		if (view == 1) MatView(assem.res[0],PETSC_VIEWER_DRAW_WORLD);
		PetscPrintf(PETSC_COMM_WORLD,"\npro:\n");
		if (view == 0) MatView(assem.pro[0],PETSC_VIEWER_STDOUT_WORLD);
		if (view == 1) MatView(assem.pro[0],PETSC_VIEWER_DRAW_WORLD);
		PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
		ISView(*(assem.topIS),PETSC_VIEWER_STDOUT_WORLD);
		ISView(*(assem.bottomIS),PETSC_VIEWER_STDOUT_WORLD);
	}

	for (int l=0;l<assem.levels-1;l++) {
		PetscPrintf(PETSC_COMM_WORLD,"res[%d]:\n",l);
		if (view == 0) MatView(assem.res[l],PETSC_VIEWER_STDOUT_WORLD);
		if (view == 1) MatView(assem.res[l],PETSC_VIEWER_DRAW_WORLD);
		PetscPrintf(PETSC_COMM_WORLD,"pro[%d]:\n",l);
		if (view == 0) MatView(assem.pro[l],PETSC_VIEWER_STDOUT_WORLD);
		if (view == 1) MatView(assem.pro[l],PETSC_VIEWER_DRAW_WORLD);
	}
	PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
}
