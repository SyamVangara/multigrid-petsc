#include "header.h"
#include <time.h>
#include <string.h>
#include <petscksp.h>

#define ERROR_MSG(message) (fprintf(stderr,"Error:%s:%d: %s\n",__FILE__,__LINE__,(message)))
#define ERROR_RETURN(message) {ERROR_MSG(message);return ierr;}
#define CHKERR_PRNT(message) {if(ierr==1) {ERROR_MSG(message);}}
#define CHKERR_RETURN(message) {if(ierr==1) {ERROR_RETURN(message);}}

#define PI 3.14159265358979323846
#define DIMENSION 2
#define FUNC(i,j) (-2*PI*PI*sin(PI*coord[0][(j)])*sin(PI*coord[1][(i)]))
#define SOL(i,j) (sin(PI*coord[0][(j)])*sin(PI*coord[1][(i)]))

#define METRICS(i,j) (metrics.data[(i)*metrics.nj+(j)])
#define F(i,j) (f.data[((i)*f.nj+(j))])
#define U(i,j) (u.data[((i)*u.nj+(j))])
#define isGRIDtoGLOBAL(l,i,j) (IsGridToGlobal[l].data[((i)*IsGridToGlobal[l].nj+(j))])
#define isGLOBALtoGRID(l,i,j) (IsGlobalToGrid[l].data[((i)*IsGlobalToGrid[l].nj+(j))])
#define isSTENCIL(l,i,j) (IsStencil[l].data[((i)*IsStencil[l].nj+(j))])
#define isRESSTENCIL(l,i,j) (IsResStencil[l].data[((i)*IsResStencil[l].nj+(j))])
#define isPROSTENCIL(l,i,j) (IsProStencil[l].data[((i)*IsProStencil[l].nj+(j))])

void GetFuncValues2d(double **coord, ArrayInt2d *IsGlobalToGrid, double *f, IsRange *range);
void GetError(double **coord, int *n, Array2d u, double *error);
void UpdateBC(double **coord, double *u, int *n);
static int ipow(int base, int exp);
void CreateArrayOfIS(int n, int levels, IS *idx);
void prolongStencil2D(double ***IH2h, int m, int n);
void restrictStencil2D(double ***Ih2H, int m, int n);
int totalUnknowns(int *n, int totalGrids);
static void GetSol(double *u, double *px, int *n, int levels, const int *ranges, int numProcs, int rank);
PetscErrorCode myMonitor(KSP ksp, PetscInt n, PetscReal rnormAtn, double *rnorm);
void stencilIndices(ArrayInt2d *IsGlobalToGrid, ArrayInt2d *IsGridToGlobal, ArrayInt2d *IsStencil, IsRange *range, int levels);
void restrictionStencilIndices(ArrayInt2d *IsGlobalToGrid, ArrayInt2d *IsGridToGlobal, ArrayInt2d *IsResStencil, IsRange *range, int levels);

void ViewMeshInfo(Mesh mesh);
void ViewGridsInfo(Indices indices);
void ViewIndexMapsInfo(Indices indices);
void ViewRangesInfo(Indices indices);
void ViewSolverInfo(Indices indices, Solver solver);
void ViewOperatorInfo(Operator op);
void ViewLinSysMatsInfo(Assembly assem, int view);
void ViewGridTransferMatsInfo(Assembly assem, int view);

int main(int argc, char *argv[]) {
	
	PetscInitialize(&argc, &argv, 0, 0);

	Problem		prob;
	Mesh		mesh;
	Indices		indices;
	Operator	op;
	Assembly	assem;
	Solver		solver;
	PostProcess	pp;	
	
	int	cyc;
	int	mappingStyleflag;

	int	ierr=0;
	int	procs, rank;
	
	MPI_Comm_size(PETSC_COMM_WORLD, &procs);
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	
	SetUpProblem(&prob);
	
	freopen("poisson.in", "r", stdin);
//	freopen("poisson.out", "w", stdout);
//	freopen("poisson.err", "w", stderr);
	
	MPI_Barrier(PETSC_COMM_WORLD);
	scanf("%d",mesh.n);	
	scanf("%d",&(solver.numIter));
	scanf("%d",&(indices.totalGrids));
	scanf("%d",&(indices.levels));
	scanf("%d",&(cyc));
	scanf("%d",&(mappingStyleflag));
//	indices.levels = 2;
	
	for (int i=1;i<DIMENSION;i++) { 
		mesh.n[i]  = mesh.n[0];      // No. of points in each dimension
	}
	for (int i=0;i<DIMENSION;i++) {
		mesh.bounds[i*2] = 0.0;    // Lower bound in each dimension
		mesh.bounds[i*2+1] = 1.0;  // Upper bound in each dimension
	}
	

//	SetUpMesh(&mesh, UNIFORM);
	SetUpMesh(&mesh, NONUNIFORM);

//	ViewMeshInfo(mesh);
	
	// Indices maps; number of local unknowns	
	indices.coarseningFactor = 2;
	SetUpIndices(&mesh, &indices);

//	ViewGridsInfo(indices);

	mapping(&indices, mappingStyleflag);

//	ViewIndexMapsInfo(indices);
//	ViewRangesInfo(indices);
	
	SetUpOperator(&indices, &op);
	GridTransferOperators(op, indices);

//	ViewOperatorInfo(op);
	
	SetUpAssembly(&indices, &assem);
	Assemble(&prob, &mesh, &indices, &op, &assem);

//	ViewLinSysMatsInfo(assem, 0);
//	ViewGridTransferMatsInfo(assem, 0);

	if (cyc == 0) SetUpSolver(&indices, &solver, VCYCLE);
	if (cyc == 1) SetUpSolver(&indices, &solver, ICYCLE);

//	ViewSolverInfo(indices, solver);

	Solve(&assem, &solver);
	SetPostProcess(&pp);
	Postprocessing(&prob, &mesh, &indices, &assem, &solver, &pp);
	
	if (rank==0) {
	
	printf("=============================================================\n");
	printf("Size:				%d x %d\n", mesh.n[0], mesh.n[1]);
	printf("Number of grids:		%d\n",op.totalGrids);
	printf("Number of levels:		%d\n",assem.levels);
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
	printf("Number of processes:		%d\n",procs);
	printf("Number of iterations:		%d\n",solver.numIter);
	printf("=============================================================\n");
	}

	DestroyPostProcess(&pp);
	DestroySolver(&solver);
	DestroyAssembly(&assem);
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

void prolongStencil2D(double ***IH2h, int m, int n){
	// Builds prolongation 2D stencilwise operator (*IH2h)
	// Stencil size: m x n
	
	//double	**IH2h;
	int	ierr;

	ierr = malloc2d(IH2h,m,n); CHKERR_PRNT("malloc failed");
	for (int lj=0;lj<3;lj++) {
 		(*IH2h)[0][lj]= 0.5 - 0.25*fabs(1-lj);
 		(*IH2h)[1][lj]= 1.0 - 0.5*fabs(1-lj);
 		(*IH2h)[2][lj]= 0.5 - 0.25*fabs(1-lj);
	}
	//return IH2h;
}

void restrictStencil2D(double ***Ih2H, int m, int n){
	// Builds prolongation 2D stencilwise operator
	// Stencil size: m x n
	
	//double **Ih2H;
	int	ierr;

	ierr = malloc2d(Ih2H,m,n); CHKERR_PRNT("malloc failed");
	for (int lj=0;lj<3;lj++) {
 		(*Ih2H)[0][lj]= 0.0;
 		(*Ih2H)[1][lj]= 0.0;
 		(*Ih2H)[2][lj]= 0.0;
	}
	(*Ih2H)[1][1] = 1.0;
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

void GetSol(double *u, double *px, int *n, int levels, const int *ranges, int numProcs, int rank) {
	
	int	r;
	
	if (rank!=0) {
		MPI_Send(px, ranges[rank+1]-ranges[rank], MPI_DOUBLE, 0, rank, PETSC_COMM_WORLD);
	}
	else if (rank==0) {
	
		int	length, n0;
		double	*x;
	
//		n0 = n[0]-2;
//		length=n0*n0;
//		for (int i=1;i<levels;i++) {
//			n0 = (n0-1)/2;
//			length = length + n0*n0;
//		}
		length = totalUnknowns(n, levels);
//		length = ((n0+1)*(n0+1)*(ipow(4,levels)-1))/(3*ipow(4,levels-1))-(2*(n0+1)*(ipow(2,levels)-1))/(ipow(2,levels-1))+levels;
		x = (double *)malloc(length*sizeof(double)); 
		
		for (int i=0;i<ranges[1];i++) x[i] = px[i];
		
		for (int i=1;i<numProcs;i++) {
			MPI_Recv(&(x[ranges[i]]), ranges[i+1]-ranges[i], MPI_DOUBLE, i, i, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	
		r = 0;
		for (int i=1;i<n[1]-1;i++) {
			for (int j=1;j<n[0]-1;j++) {
				u[i*n[0]+j] = x[r];
				r = r+1;
			}
		}
		
		free(x);
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
      
	for (int i;i<DIMENSION;i++) {
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

void ViewIndexMapsInfo(Indices indices) {
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
		for (int lg=0;lg<indices.level[l].grids;lg++) {
			for (int i=0;i<indices.level[l].grid[lg].ni;i++) {
				for (int j=0;j<indices.level[l].grid[lg].nj;j++) {
					PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = %d: level: %d: grid[%d:%d]: (%d,%d): %d\n",rank,l,lg,indices.level[l].gridId[lg],i,j,indices.level[l].grid[lg].data[i*indices.level[l].grid[lg].nj+j]);
				}
			}
		}
		PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
	
		for (int i=0;i<indices.level[l].global.ni;i++) {
			PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = %d: level: %d: global[%d] = (row = %d, col = %d, grid = %d)\n",rank,l,i,indices.level[l].global.data[i*indices.level[l].global.nj+0], indices.level[l].global.data[i*indices.level[l].global.nj+1], indices.level[l].global.data[i*indices.level[l].global.nj+2]);
		}
	}
	PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
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
		if (view == 0) MatView(assem.A[l],PETSC_VIEWER_STDOUT_WORLD);
		if (view == 1) MatView(assem.A[l],PETSC_VIEWER_DRAW_WORLD);
		VecView(assem.b[l],PETSC_VIEWER_STDOUT_WORLD);
	}
	PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
}

void ViewGridTransferMatsInfo(Assembly assem, int view) {
	// Prints the info of Assembly struct
	
	int	procs, rank;
	
	MPI_Comm_size(PETSC_COMM_WORLD, &procs);
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	
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
