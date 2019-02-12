#include "header.h"

#define ERROR_MSG(message) (fprintf(stderr,"ERROR: %s:%d: %s\n",__FILE__,__LINE__,(message)))
#define ERROR_RETURN(message) {ERROR_MSG(message);return ierr;}
#define CHKERR_PRNT(message) {if(ierr != 0) {ERROR_MSG(message);}}
#define CHKERR_RETURN(message) {if(ierr != 0) {ERROR_RETURN(message);}}

#define pERROR_MSG(message) (PetscPrintf(PETSC_COMM_WORLD,"ERROR: %s:%d: %s\n",__FILE__,__LINE__,(message)))
#define pERROR_RETURN(message) {pERROR_MSG(message);return ierr;}
#define pCHKERR_PRNT(message) {if(ierr != 0) {pERROR_MSG(message);}}
#define pCHKERR_RETURN(message) {if(ierr != 0) {pERROR_RETURN(message);}}

#define PI 3.14159265358979323846
#define MAX_DIMENSION 3
#define FUNC2D(i,j) (-2*PI*PI*sin(PI*coord[0][(j)])*sin(PI*coord[1][(i)]))
#define SOL2D(i,j) (sin(PI*coord[0][(j)])*sin(PI*coord[1][(i)]))

//static int ipow(int base, int exp);
int totalUnknowns(int *n, int totalGrids);

void PrintInfo(Grids *grids, Solver *solver);
void ViewGridsInfo(Grids grids, int verbose);
//void ViewIndicesInfo(Indices indices);
void ViewMatResInfo(Solver solver);
void ViewMatAInfo(Solver solver);
void ViewVecbInfo(Solver solver);
void ViewLevelsInfo(Solver solver, int verbose);
//void ViewIndexMapsInfoLevel(Level level, int l);
//void ViewIndexMapsInfo(Indices indices);
//void ViewRangesInfo(Indices indices);
//void ViewSolverInfo(Indices indices, Solver solver);
//void ViewOperatorInfo(Operator op);
//void ViewLinSysMatsInfo(Assembly assem, int view);
//void ViewGridTransferMatsInfo(Assembly assem, int view, int cyc);

int main(int argc, char *argv[]) {
	
	int	procs, rank;
	int 	provided;
	
	MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
	
	MPI_Comm_size(MPI_COMM_WORLD, &procs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (provided == MPI_THREAD_MULTIPLE) {
		if (rank == 0) printf("\nMPI_THREAD_MULTIPLE is provided!\n");
	} else {
		if (rank == 0) printf("\nMPI_THREAD_MULTIPLE is not provided!\n");
	}
	PetscInitialize(&argc, &argv, "poisson.in", 0);

	Problem		prob;
	Grids		grids;
//	Indices		indices;
//	Operator	op;
	Solver		solver;
//	PostProcess	pp;	

	int	ierr=0;

	SetUpProblem(&prob);
	
//	freopen("poisson.in", "r", stdin);
//	freopen("poisson.out", "w", stdout);
//	freopen("poisson.err", "w", stderr);
	
	ierr = CreateGrids(&grids); pCHKERR_PRNT("Grids creation failed");
	if (ierr == 1) {
		DestroyGrids(&grids);
		PetscFinalize();
		MPI_Finalize();
		return 0;
	}
	ViewGridsInfo(grids, 2);
	PetscBarrier(PETSC_NULL);
	ierr = CreateSolver(&grids, &solver); pCHKERR_PRNT("Solver creation failed");
	if (ierr == 1) {
		DestroySolver(&solver);
		DestroyGrids(&grids);
		PetscFinalize();
		MPI_Finalize();
		return 0;
	}
	ViewLevelsInfo(solver, 2);
	ViewMatResInfo(solver);
//	ViewMatAInfo(solver);
//	ViewVecbInfo(solver);
	
	ierr = Solve(&solver); pCHKERR_PRNT("Solver failed");
	if (ierr == 1) {
		DestroySolver(&solver);
		DestroyGrids(&grids);
		PetscFinalize();
		MPI_Finalize();
		return 0;
	}
	
//	SetUpPostProcess(&pp);
	PostProcessing(&grids, &solver);
	
	PrintInfo(&grids, &solver);
	
//	DestroyPostProcess(&pp);

	DestroySolver(&solver);
	DestroyGrids(&grids);
	PetscFinalize();
	MPI_Finalize();
	return 0;
	
	//PetscOptionsGetInt(NULL, NULL, "-map", &(mappingStyleflag), NULL);
	//
	//if (indices.levels>1 && (cyc==3 || cyc==4 || cyc==7)) {
	//	PetscPrintf(PETSC_COMM_WORLD, "For now only one level is allowed for delayed cycling\n"); 
	//	PetscFinalize();
	//	MPI_Finalize();
	//	return 0;
	//}

//	//SetUpMesh(&mesh, meshflag);
//	//if (meshflag == 0) SetUpMesh(&mesh, UNIFORM);
//	//if (meshflag == 1) SetUpMesh(&mesh, NONUNIFORM1);
//	//if (meshflag == 2) SetUpMesh(&mesh, NONUNIFORM2);

//	//ViewMeshInfo(mesh);
	//
	//indices.coarseningFactor = 2;
	//SetUpIndices(&mesh, &indices);

//	//ViewIndicesInfo(indices);

	//mapping(&indices, mappingStyleflag);
	//
//	//ViewIndexMapsInfo(indices);
//	//ViewRangesInfo(indices);
	//
	//SetUpOperator(&indices, &op);
	//GridTransferOperators(op, indices);

//	//ViewOperatorInfo(op);
	//
	//if (cyc == 0) SetUpSolver(&indices, &solver, VCYCLE);
	//if (cyc == 1) SetUpSolver(&indices, &solver, ICYCLE);
	//if (cyc == 2) SetUpSolver(&indices, &solver, ECYCLE);
	//if (cyc == 3) SetUpSolver(&indices, &solver, D1CYCLE);
	//if (cyc == 4) SetUpSolver(&indices, &solver, D2CYCLE);
	//if (cyc == 7) SetUpSolver(&indices, &solver, D1PSCYCLE);
	//if (cyc == 8) SetUpSolver(&indices, &solver, PetscPCMG);
	//if (cyc == 9) SetUpSolver(&indices, &solver, FILTER);
	//if (cyc == 10) SetUpSolver(&indices, &solver, VFILTER);
	//if (cyc == 11) SetUpSolver(&indices, &solver, ADDITIVE);
	//if (cyc == 12) SetUpSolver(&indices, &solver, ADDITIVEScaled);

//	//ViewSolverInfo(indices, solver);

	//Assemble(&prob, &mesh, &indices, &op, &solver);

//	//ViewLinSysMatsInfo(*(solver.assem), 0);
//	//ViewGridTransferMatsInfo(*(solver.assem), 0, cyc);
	//
	//Solve(&solver);
	//
	//SetUpPostProcess(&pp);
	//Postprocessing(&prob, &mesh, &indices, &solver, &pp);
	//
	//PrintInfo(prob, mesh, indices, op, solver, pp, cyc, meshflag, mappingStyleflag);
	//
	//DestroyPostProcess(&pp);
	//DestroySolver(&solver);
	//DestroyOperator(&op);
	//DestroyIndices(&indices);
	//DestroyMesh(&mesh);
	//PetscFinalize();
	//MPI_Finalize();

	//return 0;
}

//int ipow(int base, int exp) {
//
//	int result = 1;
//	while (exp) {
//		if (exp & 1)
//			result *= base;
//		exp >>= 1;
//		base *= base;
//	}
//	return result;
//}

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

void PrintInfo(Grids *grids, Solver *solver) {
	// Prints complete some info of problem, grids, solver
	
	int	procs, rank;

	MPI_Comm_size(PETSC_COMM_WORLD, &procs);
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	
	if (rank==0) {
	Grid	*fgrid = grids->grid;
	Levels	*levels = solver->levels;
	Level	*level = levels->level;
	int	dimension = grids->topo->dimension;
	int	*gridtype = grids->topo->gridtype;

	printf("\n=============================================================\n");
	printf("Dimension:			%d\n", dimension);
	printf("Grid size:			(%d", fgrid->n[0]);
	for (int i=1; i<dimension; i++) {
		printf(", %d", fgrid->n[i]);
	}
	printf(")\n");
	printf("Characteristic length:		%lf\n", fgrid->h);
	if (gridtype[0] == 0) printf("Grid Type:			(Uniform");
	if (gridtype[0] == 1) printf("Grid Type:			(Non-Uniform-Cos");
	if (gridtype[0] == 2) printf("Grid Type:			(Non-Uniform-Exp");
	for (int i=1; i<dimension; i++) {
		if (gridtype[i] == 0) printf(", Uniform");
		if (gridtype[i] == 1) printf(", Non-Uniform-Cos");
		if (gridtype[i] == 2) printf(", Non-Uniform-Exp");
	}
	printf(")\n");
//	double	commTocomp = 1.0;
//	double	*para = fgrid->para;
//	for (int i=0; i<dimension; i++) {
//		commTocomp *= fgrid->un[i];
//	}
//	commTocomp = para[0]/commTocomp;
//	printf("TotalCommCost =			%d\n", (int)para[0]);
//	printf("MaxLoad =			%d\n", (int)para[1]);
//	printf("Comm-to-Comp =			%lf\n", commTocomp);
//	printf("LoadFactor =			%lf\n", para[2]);
//	printf("nInterfaces =			%d\n", (int)para[3]);
	printf("Number of grids:		%d\n", grids->ngrids);
	printf("Number of levels:		%d\n", levels->nlevels);
	printf("Number of grids per level:	(%d", level[0].ngrids);
	for (int l=1;l<levels->nlevels;l++) {
		printf(", %d", level[l].ngrids);
	}
	printf(")\n");
	if (solver->cycle == 1) printf("Cycle :				V-Cycle\n");
	if (solver->cycle == 0) printf("Cycle :				No MG\n");
	
	printf("Number of processes:		%d\n",procs);
	printf("Number of iterations:		%d\n",solver->numIter);
	printf("=============================================================\n");
	}

}

void ViewTopoInfo(Topo *topo) {
	// Prints the info in Topo data structure
	
	int	dimension = topo->dimension;
	
	PetscPrintf(PETSC_COMM_WORLD,"Topo:\n");
	PetscPrintf(PETSC_COMM_WORLD,"dimension = %d\n", dimension);

	for (int dim = 0; dim<dimension; dim++) {
		PetscPrintf(PETSC_COMM_WORLD,"gridtype[%d] = %d  ", dim, topo->gridtype[dim]);
	}
	PetscPrintf(PETSC_COMM_WORLD,"\n");

	for (int dim = 0; dim<dimension; dim++) {
		PetscPrintf(PETSC_COMM_WORLD,"procs_%d = %d  ", dim, topo->dimProcs[dim]);
	}
	PetscPrintf(PETSC_COMM_WORLD,"\n");

	for (int dim = 0; dim<dimension; dim++) {
		for (int i=0; i<2; i++) {
			PetscPrintf(PETSC_COMM_WORLD,"bounds[%d][%d] = %f  ", dim, i, topo->bounds[dim][i]);
		}
	}
	PetscPrintf(PETSC_COMM_WORLD,"\n\n");
}

void ViewEblockInfo(Nblock (*eblock)[2][2], int dimension) {
	// Prints the info in edge blocks
	
	int	procs, rank;
	
	MPI_Comm_size(MPI_COMM_WORLD, &procs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	PetscPrintf(PETSC_COMM_WORLD,"Neighbor-Edge blocks: \n");
	PetscSynchronizedPrintf(PETSC_COMM_WORLD, "Rank = %d: \n", rank);
	for (int i=0; i<dimension; i++) {
		for (int j=0; j<2; j++) {
			for (int k=0; k<2; k++) {
				PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Direction(%d, %d, %d):	rank = %d", i, 2*j-1, 2*k-1, eblock[i][j][k].rank);
				PetscSynchronizedPrintf(PETSC_COMM_WORLD, ";	blockID = ( ");
				for (int dim=0; dim<dimension; dim++)
					PetscSynchronizedPrintf(PETSC_COMM_WORLD,"%d ", eblock[i][j][k].blockID[dim]);
				PetscSynchronizedPrintf(PETSC_COMM_WORLD, ");	ln =");
				for (int dim=0; dim<dimension; dim++)
					PetscSynchronizedPrintf(PETSC_COMM_WORLD," %d", eblock[i][j][k].ln[dim]);
				PetscSynchronizedPrintf(PETSC_COMM_WORLD, "\n");
			}
		}
	}
	PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
}

void ViewCblockInfo(Nblock (*cblock)[2][2], int dimension) {
	// Prints the info in corner blocks
	
	int	procs, rank;
	
	MPI_Comm_size(MPI_COMM_WORLD, &procs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	PetscPrintf(PETSC_COMM_WORLD,"Neighbor-Corner blocks: \n");
	PetscSynchronizedPrintf(PETSC_COMM_WORLD, "Rank = %d: \n", rank);
	for (int i=0; i<2; i++) {
		for (int j=0; j<2; j++) {
			for (int k=0; k<2; k++) {
				PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Direction(%d, %d, %d):	rank = %d", 2*i-1, 2*j-1, 2*k-1, cblock[i][j][k].rank);
				PetscSynchronizedPrintf(PETSC_COMM_WORLD, ";	blockID = ( ");
				for (int dim=0; dim<dimension; dim++)
					PetscSynchronizedPrintf(PETSC_COMM_WORLD,"%d ", cblock[i][j][k].blockID[dim]);
				PetscSynchronizedPrintf(PETSC_COMM_WORLD, ");	ln =");
				for (int dim=0; dim<dimension; dim++)
					PetscSynchronizedPrintf(PETSC_COMM_WORLD," %d", cblock[i][j][k].ln[dim]);
				PetscSynchronizedPrintf(PETSC_COMM_WORLD, "\n");
			}
		}
	}
	PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
}

void ViewNblockInfo(Nblock (*nblock)[2], int dimension) {
	// Prints the info in face blocks
	
	int	procs, rank;
	
	MPI_Comm_size(MPI_COMM_WORLD, &procs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	PetscPrintf(PETSC_COMM_WORLD,"Neighbor-blocks: \n");
	PetscSynchronizedPrintf(PETSC_COMM_WORLD, "Rank = %d: \n", rank);
	for (int i=0; i<dimension; i++) {
		for (int j=0; j<2; j++) {
			PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Direction(%d, %d):	rank = %d", i, 2*j-1, nblock[i][j].rank);
			PetscSynchronizedPrintf(PETSC_COMM_WORLD, ";	blockID = ( ");
			for (int k=0; k<dimension; k++)
				PetscSynchronizedPrintf(PETSC_COMM_WORLD,"%d ", nblock[i][j].blockID[k]);
			PetscSynchronizedPrintf(PETSC_COMM_WORLD, ");	ln =");
			for (int k=0; k<dimension; k++)
				PetscSynchronizedPrintf(PETSC_COMM_WORLD," %d", nblock[i][j].ln[k]);
			PetscSynchronizedPrintf(PETSC_COMM_WORLD, "\n");
		}
	}
	PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
}

void ViewGridInfo(Grid grid, int verbose) {
	// Prints the info in Grid data structure

	int	dimension = grid.topo->dimension;
	int	*blockID = grid.topo->blockID;
	int	**range = grid.range;
//	int	*inc = grid.inc;
//	for (int i=0; i<dimension; i++)
//		range[i] = grid.range[i];
	double	*para = grid.para;
	double	**coord = grid.coord;
	double	**dx = grid.dx;
	int	procs, rank;
	
	MPI_Comm_size(MPI_COMM_WORLD, &procs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	PetscPrintf(PETSC_COMM_WORLD,"Grid-%d:	\n", grid.id);
	PetscPrintf(PETSC_COMM_WORLD,"n =");
	for (int dim = 0; dim<dimension; dim++) {
		PetscPrintf(PETSC_COMM_WORLD," %d", grid.n[dim]);
	}
	PetscPrintf(PETSC_COMM_WORLD,";	un =");
	for (int dim = 0; dim<dimension; dim++) {
		PetscPrintf(PETSC_COMM_WORLD," %d", grid.un[dim]);
	}
	PetscPrintf(PETSC_COMM_WORLD,";	TotalN = %d (%d)", grid.tn, grid.tun);
	PetscPrintf(PETSC_COMM_WORLD,"\n");
	
//	double	commTocomp = para[0]/(grid.un[0]*grid.un[1]*grid.un[2]);
	double	commTocomp = 1.0;
	for (int i=0; i<dimension; i++) {
		commTocomp *= grid.un[i];
	}
	commTocomp = para[0]/commTocomp;
	PetscPrintf(PETSC_COMM_WORLD,"TotalCommCost = %d, MaxLoad = %d, Comm-to-Comp = %lf, LoadFactor = %lf, nInterfaces = %d\n", (int)para[0], (int)para[1], commTocomp, para[2], (int)para[3]);
	
	PetscPrintf(PETSC_COMM_WORLD,"h = %f\n", grid.h);
      	
	if (verbose > 0) {
		PetscSynchronizedPrintf(PETSC_COMM_WORLD, "Rank = %d:	blockID = ( ", rank);
		for (int i=0; i<dimension; i++)
			PetscSynchronizedPrintf(PETSC_COMM_WORLD,"%d ", blockID[i]);
		PetscSynchronizedPrintf(PETSC_COMM_WORLD, ");	range =");
		for (int i=0; i<dimension; i++)
			PetscSynchronizedPrintf(PETSC_COMM_WORLD," (%d-%d)", range[i][blockID[i]], range[i][blockID[i]+1]);
		PetscSynchronizedPrintf(PETSC_COMM_WORLD, ";	ln =");
		for (int i=0; i<dimension; i++)
			PetscSynchronizedPrintf(PETSC_COMM_WORLD," %d", grid.ln[i]);
		PetscSynchronizedPrintf(PETSC_COMM_WORLD,";	LocalTotalN = %d", grid.tln);
		PetscSynchronizedPrintf(PETSC_COMM_WORLD, "\n");
		PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
	}
	if (verbose > 1) {
		ViewNblockInfo(grid.nblock, dimension);
		if (dimension == 3) ViewEblockInfo(grid.eblock, dimension);
		ViewCblockInfo(grid.cblock, dimension);
		for (int i=0;i<dimension;i++) {
			PetscPrintf(PETSC_COMM_WORLD,"coord[%d]:",i);
			for (int j=0;j<grid.n[i];j++) {
				PetscPrintf(PETSC_COMM_WORLD," %f ",coord[i][j]);
			}
			PetscPrintf(PETSC_COMM_WORLD,"\n");
		}
		for (int i=0;i<dimension;i++) {
			PetscPrintf(PETSC_COMM_WORLD,"dx[%d]:",i);
			for (int j=0;j<grid.n[i]-1;j++) {
				PetscPrintf(PETSC_COMM_WORLD," %f ",dx[i][j]);
			}
			PetscPrintf(PETSC_COMM_WORLD,"\n");
		}
	}
	PetscPrintf(PETSC_COMM_WORLD,"\n");
}

void ViewGridsInfo(Grids grids, int verbose) {
	// Prints the info of Grids data structure
	
	int ngrids = grids.ngrids;
	PetscPrintf(PETSC_COMM_WORLD,"Grids: \n");
	PetscPrintf(PETSC_COMM_WORLD,"Total no. of grids = %d\n", ngrids);
	for (int i=0; i<ngrids-1; i++) {
		PetscPrintf(PETSC_COMM_WORLD,"coarsening factors[%d] = ", i);
		for (int dim=0; dim<grids.topo->dimension; dim++) {
			PetscPrintf(PETSC_COMM_WORLD,"%d  ", grids.cfactor[i][dim]);
		}
		PetscPrintf(PETSC_COMM_WORLD,"\n");
	}
	PetscPrintf(PETSC_COMM_WORLD,"\n");
	
	ViewTopoInfo(grids.topo);
	for (int i=0; i<ngrids; i++) {
		ViewGridInfo(grids.grid[i], verbose);
	}
}

void ViewBCindicesInfo(BCindices (*bcindices)[2], int dimension) {
	// Prints the info in Nblock
	
	int	procs, rank;
	
	MPI_Comm_size(MPI_COMM_WORLD, &procs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	PetscSynchronizedPrintf(PETSC_COMM_WORLD, "Rank = %d: \n", rank);
	for (int i=0; i<dimension; i++) {
		for (int j=0; j<2; j++) {
			PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Dim-%d, Dir(%d):	rank = %d", i, 2*j-1, bcindices[i][j].rank);
			PetscSynchronizedPrintf(PETSC_COMM_WORLD, ";	blockID = ( ");
			for (int k=0; k<dimension; k++)
				PetscSynchronizedPrintf(PETSC_COMM_WORLD,"%d ", bcindices[i][j].blockID[k]);
			PetscSynchronizedPrintf(PETSC_COMM_WORLD, ");	start = %ld", bcindices[i][j].bcStartIndex);
			PetscSynchronizedPrintf(PETSC_COMM_WORLD, "	gstart = %ld", bcindices[i][j].bcGStartIndex);
			PetscSynchronizedPrintf(PETSC_COMM_WORLD, "	Inc = ( ");
			for (int k=0; k<dimension; k++)
				PetscSynchronizedPrintf(PETSC_COMM_WORLD,"%ld ", bcindices[i][j].bcInc[k]);
			PetscSynchronizedPrintf(PETSC_COMM_WORLD, ")\n");
		}
	}
	PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);

}

void ViewEBCindicesInfo(BCindices (*bcindices)[2][2], int dimension) {
	// Prints the info in Nblock
	
	int	procs, rank;
	
	MPI_Comm_size(MPI_COMM_WORLD, &procs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	PetscSynchronizedPrintf(PETSC_COMM_WORLD, "Rank = %d: \n", rank);
	for (int i=0; i<dimension; i++) {
		for (int j=0; j<2; j++) {
			for (int k=0; k<2; k++) {
				PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Dim-%d, Dir(%d, %d):	rank = %d", i, 2*j-1, 2*k-1, bcindices[i][j][k].rank);
				PetscSynchronizedPrintf(PETSC_COMM_WORLD, ";	blockID = ( ");
				for (int dim=0; dim<dimension; dim++)
					PetscSynchronizedPrintf(PETSC_COMM_WORLD,"%d ", bcindices[i][j][k].blockID[dim]);
				PetscSynchronizedPrintf(PETSC_COMM_WORLD, ");	start = %ld", bcindices[i][j][k].bcStartIndex);
				PetscSynchronizedPrintf(PETSC_COMM_WORLD, "	gstart = %ld", bcindices[i][j][k].bcGStartIndex);
				PetscSynchronizedPrintf(PETSC_COMM_WORLD, "	Inc = ( ");
				for (int dim=0; dim<dimension; dim++)
					PetscSynchronizedPrintf(PETSC_COMM_WORLD,"%ld ", bcindices[i][j][k].bcInc[dim]);
				PetscSynchronizedPrintf(PETSC_COMM_WORLD, ")\n");
			}
		}
	}
	PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
}

void ViewCBCindicesInfo(BCindices (*bcindices)[2][2], int dimension) {
	// Prints the info in Nblock
	
	int	procs, rank;
	
	MPI_Comm_size(MPI_COMM_WORLD, &procs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	PetscSynchronizedPrintf(PETSC_COMM_WORLD, "Rank = %d: \n", rank);
	for (int i=0; i<2; i++) {
		for (int j=0; j<2; j++) {
			for (int k=0; k<dimension-1; k++) {
				PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Dir(%d,%d,%d):	rank = %d", 2*i-1, 2*j-1, 2*k-1, bcindices[i][j][k].rank);
				PetscSynchronizedPrintf(PETSC_COMM_WORLD, ";	blockID = ( ");
				for (int dim=0; dim<dimension; dim++)
					PetscSynchronizedPrintf(PETSC_COMM_WORLD,"%d ", bcindices[i][j][k].blockID[dim]);
				PetscSynchronizedPrintf(PETSC_COMM_WORLD, ");	start = %ld", bcindices[i][j][k].bcStartIndex);
				PetscSynchronizedPrintf(PETSC_COMM_WORLD, "	gstart = %ld", bcindices[i][j][k].bcGStartIndex);
				PetscSynchronizedPrintf(PETSC_COMM_WORLD, "	Inc = ( ");
				for (int dim=0; dim<dimension; dim++)
					PetscSynchronizedPrintf(PETSC_COMM_WORLD,"%ld ", bcindices[i][j][k].bcInc[dim]);
				PetscSynchronizedPrintf(PETSC_COMM_WORLD, ")\n");
			}
		}
	}
	PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
}

void ViewLevelInfo(Level level, int dimension, int verbose) {
	// Prints the info of Level data structure
	
	int	procs, rank;
	
	MPI_Comm_size(MPI_COMM_WORLD, &procs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	PetscPrintf(PETSC_COMM_WORLD,"No. of grids = %d\n", level.ngrids);
	
	PetscPrintf(PETSC_COMM_WORLD,"gridIDs = ");
	for (int lg=0;lg<level.ngrids;lg++) {
		PetscPrintf(PETSC_COMM_WORLD,"%d ",level.gridId[lg]);
	}
	PetscPrintf(PETSC_COMM_WORLD,"\n");
	
	if (verbose > 0) {
		PetscPrintf(PETSC_COMM_WORLD, "Global Level index ranges = \n");
		PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank: %d; ",rank);
		for (int lg=0;lg<level.ngrids;lg++) {
			PetscSynchronizedPrintf(PETSC_COMM_WORLD,"(%ld-%ld: %d %d %d) ",level.ranges[lg], level.ranges[lg+1], level.inc[lg][0], level.inc[lg][1], level.inc[lg][2]);
		}
		PetscSynchronizedPrintf(PETSC_COMM_WORLD, "\n");
		PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
		PetscPrintf(PETSC_COMM_WORLD, "Global Grid index ranges = \n");
		PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank: %d; ",rank);
		for (int lg=0;lg<level.ngrids;lg++) {
			PetscSynchronizedPrintf(PETSC_COMM_WORLD,"(%ld-%ld) ",level.granges[lg][0], level.granges[lg][1]);
		}
		PetscSynchronizedPrintf(PETSC_COMM_WORLD, "\n");
		PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
	} 
	
	if (verbose > 1) {
		PetscPrintf(PETSC_COMM_WORLD,"\nNeighbor block indices: \n");
		for (int lg=0; lg<level.ngrids; lg++) {
			PetscPrintf(PETSC_COMM_WORLD,"\nGrid-%d: \n", level.gridId[lg]);
			PetscPrintf(PETSC_COMM_WORLD,"Face BC: \n");
			ViewBCindicesInfo(level.bcindices[lg], dimension);
			if (dimension == 3) {
				PetscPrintf(PETSC_COMM_WORLD,"\nEdge BC: \n");
				ViewEBCindicesInfo(level.ebcindices[lg], dimension);
			}
			PetscPrintf(PETSC_COMM_WORLD,"\nCorner BC: \n");
			ViewCBCindicesInfo(level.cbcindices[lg], dimension);
		}
	}
	PetscPrintf(PETSC_COMM_WORLD,"\n");
}

void ViewLevelsInfo(Solver solver, int verbose) {
	// Prints the info of Indices data structure
	
	int	procs, rank;
	
	MPI_Comm_size(PETSC_COMM_WORLD, &procs);
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

	Levels *levels = solver.levels;
	
	int dimension = levels->dimension;
	PetscPrintf(PETSC_COMM_WORLD,"Levels: \n");
	PetscPrintf(PETSC_COMM_WORLD,"Total no. of levels = %d\n", levels->nlevels);
	for (int l=0;l<levels->nlevels;l++) {
		PetscPrintf(PETSC_COMM_WORLD,"Level-%d:\n", l);
		ViewLevelInfo(levels->level[l], dimension, verbose);
	}
	PetscPrintf(PETSC_COMM_WORLD,"\n");
}

void ViewMatResInfo(Solver solver) {
	// Prints the Mat A info of each level
	
	int nlevels = solver.levels->nlevels;
	int lngrids = solver.levels->level[nlevels-1].ngrids;
	int nres = solver.levels->level[nlevels-1].gridId[lngrids-1];
	Mat *res = solver.levels->res;
	
	for (int l=0; l<nres; l++) MatView(res[l], PETSC_VIEWER_STDOUT_WORLD);
}

void ViewMatAInfo(Solver solver) {
	// Prints the Mat A info of each level
	
	int nlevels = solver.levels->nlevels;
	Mat *A = solver.levels->A;
	
	for (int l=0; l<nlevels; l++) MatView(A[l], PETSC_VIEWER_STDOUT_WORLD);
}

void ViewVecbInfo(Solver solver) {
	// Prints the Vec b info of first of first level
	
	PetscPrintf(PETSC_COMM_WORLD,"Vec b of Level 0: \n");
	VecView(solver.levels->b[0], PETSC_VIEWER_STDOUT_WORLD);
}

//void ViewIndicesInfo(Indices indices) {
//	// Prints the info of Indices data structure
//	
//	int	procs, rank;
//	
//	MPI_Comm_size(PETSC_COMM_WORLD, &procs);
//	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
//
//	PetscPrintf(PETSC_COMM_WORLD,"Indices: \n");
//	PetscPrintf(PETSC_COMM_WORLD,"Total no. of levels = %d\n", indices.levels);
//	for (int l=0;l<indices.levels;l++) {
//		PetscPrintf(PETSC_COMM_WORLD,"Level-%d:\n", l);
//		ViewLevelInfo(indices.level[l]);
//	}
//	PetscPrintf(PETSC_COMM_WORLD,"\n");
//}

//void ViewIndexMapsInfoLevel(Level level, int l) {
///********************************************************************************
// *
// * Print the Index maps of the given level
// *
// ********************************************************************************/ 
//	int	procs, rank;
//	
//	MPI_Comm_size(PETSC_COMM_WORLD, &procs);
//	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
//	
//	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = %d; level: %d; ranges: ", rank, l);
//	for (int p=0;p<procs+1;p++) {
//		PetscSynchronizedPrintf(PETSC_COMM_WORLD,"%d ",level.ranges[p]);
//	}
//	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"\n");
//	PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
//	
//	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = %d: level = %d: Map from grid indices to global index:\n",rank,l);
//	for (int lg=0; lg<level.grids; lg++) {
//		for (int i=0; i<level.grid[lg].ni; i++) {
//			for (int j=0; j<level.grid[lg].nj; j++) {
//				PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = %d: level = %d: level_gridId = %d: gridId = %d: (row = %d, col = %d) -> (global_index = %d)\n",rank,l,lg,level.gridId[lg],i,j,level.grid[lg].data[i*level.grid[lg].nj+j]);
//			}
//		}
//	}
//	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"\n");
//	PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
//	
//	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = %d: level = %d: Map from global index to grid indices:\n",rank,l);
//	for (int i=0;i<level.global.ni;i++) {
//		PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = %d: level = %d: (global_index = %d) -> (row = %d, col = %d, gridId = %d)\n",rank,l,i,level.global.data[i*level.global.nj+0], level.global.data[i*level.global.nj+1], level.global.data[i*level.global.nj+2]);
//	}
//	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"\n");
//	PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
//}
//
//void ViewIndexMapsInfo(Indices indices) {
//	// Prints the info of index maps between global and grids
//	
//	int	procs, rank;
//	
//	MPI_Comm_size(PETSC_COMM_WORLD, &procs);
//	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
//	
//	for (int l=0;l<indices.levels;l++) {
//		ViewIndexMapsInfoLevel(indices.level[l], l);
//	}
//}
//
//void ViewRangesInfo(Indices indices) {
//	// Prints the info of index maps between global and grids
//	
//	int	procs, rank;
//	
//	MPI_Comm_size(PETSC_COMM_WORLD, &procs);
//	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
//	
//	for (int l=0;l<indices.levels;l++) {
//		PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = %d; level: %d; ranges: ", rank, l);
//		for (int p=0;p<procs+1;p++) {
//			PetscSynchronizedPrintf(PETSC_COMM_WORLD,"%d ",indices.level[l].ranges[p]);
//		}
//		PetscSynchronizedPrintf(PETSC_COMM_WORLD,"\n");
//	}
//	PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
//}
//
//void ViewSolverInfo(Indices indices, Solver solver) {
//	// Prints the info in Solver struct
//	
//	int	procs, rank;
//	
//	MPI_Comm_size(PETSC_COMM_WORLD, &procs);
//	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
//	
//	if (solver.cycle==VCYCLE) {
//		PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = %d; numIter = %d; Cycle = VCYCLE\n",rank,solver.numIter);
//	} else {
//		PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = %d; numIter = %d; Cycle = ICYCLE\n",rank,solver.numIter);
//	}
//	PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
//}
//
//void ViewOperatorInfo(Operator op) {
//	// Prints the info in Operator struct
//	
//	int	procs, rank;
//	
//	MPI_Comm_size(PETSC_COMM_WORLD, &procs);
//	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
//	
//	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = %d; Total num of grids = %d:\n",rank,op.totalGrids);
//	PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
//	for (int l=0;l<op.totalGrids-1;l++) {
//		PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = %d; res[%d]:\n",rank,l);
//		for (int i=0;i<op.res[l].ni;i++) {
//			for (int j=0;j<op.res[l].nj;j++) {
//				PetscSynchronizedPrintf(PETSC_COMM_WORLD,"%f ",op.res[l].data[i*op.res[l].nj+j]);
//			}
//			PetscSynchronizedPrintf(PETSC_COMM_WORLD,"\n");
//		}
//		PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = %d; pro[%d]:\n",rank,l);
//		for (int i=0;i<op.pro[l].ni;i++) {
//			for (int j=0;j<op.pro[l].nj;j++) {
//				PetscSynchronizedPrintf(PETSC_COMM_WORLD,"%f ",op.pro[l].data[i*op.pro[l].nj+j]);
//			}
//			PetscSynchronizedPrintf(PETSC_COMM_WORLD,"\n");
//		}
//	}
//	PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
//}
//
//void ViewLinSysMatsInfo(Assembly assem, int view) {
//	// Prints the info of Assembly struct
//	
//	int	procs, rank;
//	
//	MPI_Comm_size(PETSC_COMM_WORLD, &procs);
//	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
//	
//	for (int l=0;l<assem.levels;l++) {
//		PetscPrintf(PETSC_COMM_WORLD,"A[%d]:\n",l);
//		if (view == 0) {
//			MatView(assem.A[l],PETSC_VIEWER_STDOUT_WORLD);
//			MatView(assem.A2[l],PETSC_VIEWER_STDOUT_WORLD);
//		}
//		if (view == 1) {
//			MatView(assem.A[l],PETSC_VIEWER_DRAW_WORLD);
//			MatView(assem.A2[l],PETSC_VIEWER_DRAW_WORLD);
//		}
//		VecView(assem.b[l],PETSC_VIEWER_STDOUT_WORLD);
//	}
//	PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
//}
//
//void ViewGridTransferMatsInfo(Assembly assem, int view, int cyc) {
//	// Prints the info of Assembly struct
//	
//	int	procs, rank;
//	
//	MPI_Comm_size(PETSC_COMM_WORLD, &procs);
//	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
//
//	if (cyc == 3)
//	{
//		PetscPrintf(PETSC_COMM_WORLD,"\nres:\n");
//		if (view == 0) MatView(assem.res[0],PETSC_VIEWER_STDOUT_WORLD);
//		if (view == 1) MatView(assem.res[0],PETSC_VIEWER_DRAW_WORLD);
//		PetscPrintf(PETSC_COMM_WORLD,"\npro:\n");
//		if (view == 0) MatView(assem.pro[0],PETSC_VIEWER_STDOUT_WORLD);
//		if (view == 1) MatView(assem.pro[0],PETSC_VIEWER_DRAW_WORLD);
//		PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
//		ISView(*(assem.topIS),PETSC_VIEWER_STDOUT_WORLD);
//		ISView(*(assem.bottomIS),PETSC_VIEWER_STDOUT_WORLD);
//	}
//
//	for (int l=0;l<assem.levels-1;l++) {
//		PetscPrintf(PETSC_COMM_WORLD,"res[%d]:\n",l);
//		if (view == 0) MatView(assem.res[l],PETSC_VIEWER_STDOUT_WORLD);
//		if (view == 1) MatView(assem.res[l],PETSC_VIEWER_DRAW_WORLD);
//		PetscPrintf(PETSC_COMM_WORLD,"pro[%d]:\n",l);
//		if (view == 0) MatView(assem.pro[l],PETSC_VIEWER_STDOUT_WORLD);
//		if (view == 1) MatView(assem.pro[l],PETSC_VIEWER_DRAW_WORLD);
//	}
//	PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
//}
