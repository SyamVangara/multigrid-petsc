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
int totalUnknowns(int *n, int levels);
static void GetSol(double *u, double *px, int *n, int levels, const int *ranges, int numProcs, int rank);
PetscErrorCode myMonitor(KSP ksp, PetscInt n, PetscReal rnormAtn, double *rnorm);
void stencilIndices(ArrayInt2d *IsGlobalToGrid, ArrayInt2d *IsGridToGlobal, ArrayInt2d *IsStencil, IsRange *range, int levels);
void restrictionStencilIndices(ArrayInt2d *IsGlobalToGrid, ArrayInt2d *IsGridToGlobal, ArrayInt2d *IsResStencil, IsRange *range, int levels);

void ViewMeshInfo(Mesh mesh);
void ViewGridsInfo(Indices indices);
void ViewIndexMapsInfo(Indices indices);
void ViewSolverInfo(Indices indices, Solver solver);
void ViewOperatorInfo(Operator op);
void ViewLinSysInfo(Assembly assem, int view);
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

	int	ierr=0;
	int	procs, rank;
	
	MPI_Comm_size(PETSC_COMM_WORLD, &procs);
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	
	SetUpProblem(&prob);
	
/*
	//double	weight=(2.0/3.0);
	int	ierr=0;
	Array2d	metrics;
	double	error[3], As[5], *px, *rnorm;
	Array2d	u;
	double	**opIH2h, **opIh2H;
	FILE	*solData, *errData, *resData;
	
	int	ln; //number of local unknowns

	int	rowStart, rowEnd;
	
		IsRange	*range;
		IsRange	*gridId;
	StencilIndices	indices;
//	ArrayInt2d	*IsStencil, *IsResStencil, *IsProStencil;
	
//	KSP	solver;
//	PC	pc;
//	Mat	A;
//	Vec	b, x;

*/	
	freopen("poisson.in", "r", stdin);
//	freopen("poisson.out", "w", stdout);
//	freopen("poisson.err", "w", stderr);
	
	MPI_Barrier(PETSC_COMM_WORLD);
	scanf("%d",mesh.n);	
	scanf("%d",&(solver.numIter));
	scanf("%d",&(indices.totalGrids));
	indices.levels = indices.totalGrids;
	
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
	
//	IsStencil = malloc(levels*sizeof(ArrayInt2d));
//	IsResStencil = malloc(levels*sizeof(ArrayInt2d));
	
	// Indices maps; number of local unknowns	
	indices.coarseningFactor = 2;
	SetUpIndices(&mesh, &indices);
//	ViewGridsInfo(indices);
	mapping(&indices);
//	ViewIndexMapsInfo(indices);
	
	SetUpOperator(&indices, &op);
	GridTransferOperators(op, indices);
//	ViewOperatorInfo(op);
	
	SetUpAssembly(&indices, &assem);
	Assemble(&prob, &mesh, &indices, &op, &assem);
//	ViewLinSysInfo(assem, 0);
//	printf("I am here\n");
	ViewGridTransferMatsInfo(assem, 0);

	SetUpSolver(&indices, &solver, VCYCLE);
//	ViewSolverInfo(indices, solver);
//	ln = range[1]-range[0];
//	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = %d: from %d to %d; local unknowns = %d\n",rank,range[0],range[1]-1,range[1]-range[0]);
//	PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
	
	
//	stencilIndices(IsGlobalToGrid, IsGridToGlobal, IsStencil, range, levels);

//	for (int l=0;l<levels;l++) {
//		for (int i=range[l].start;i<range[l].end;i++) {
//			PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = %d: level = %d: (%d,%d): isStencil[%d]: %d %d %d %d %d\n",rank,l,isGLOBALtoGRID(l,i,0),isGLOBALtoGRID(l,i,1),i-range[l].start,isSTENCIL(l,i-range[l].start,0),isSTENCIL(l,i-range[l].start,1),isSTENCIL(l,i-range[l].start,2),isSTENCIL(l,i-range[l].start,3),isSTENCIL(l,i-range[l].start,4));
//		}
//	PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
//	}

//	restrictionStencilIndices(IsGlobalToGrid, IsGridToGlobal, IsResStencil, range, levels);

//	for (int l=1;l<levels;l++) {
//		for (int i=range[l].start;i<range[l].end;i++) {
//			PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = %d: level = %d: (%d,%d): isResStencil[%d]:",rank,l,isGLOBALtoGRID(l,i,0),isGLOBALtoGRID(l,i,1),i-range[l].start);
//			for (int j=0;j<9;j++) { 
//				PetscSynchronizedPrintf(PETSC_COMM_WORLD," %d ",isRESSTENCIL(l-1,i-range[l].start,j));
//			}
//			PetscSynchronizedPrintf(PETSC_COMM_WORLD,"\n");
//		}
//	PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
//	}
	
//	IsProStencil = IsResStencil; // Prolongation stencil is same as restriction stencil
	
//FD Lab 
//		for (int i=range[l].start;i<range[l].end;i++) {
//			PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = %d: level = %d: (%d,%d): isPROStencil[%d]:",rank,l,isGLOBALtoGRID(l,i,0),isGLOBALtoGRID(l,i,1),i-range[l].start);
//			for (int j=0;j<9;j++) { 
//				PetscSynchronizedPrintf(PETSC_COMM_WORLD," %d ",isPROSTENCIL(l-1,i-range[l].start,j));
//			}
//			PetscSynchronizedPrintf(PETSC_COMM_WORLD,"\n");
//		}
//	PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
//	}
	
//	ierr = NonUniformMeshY(&coord,n,bounds,&h,DIMENSION,&TransformFunc); CHKERR_PRNT("meshing failed");
//	ierr = MetricCoefficients2D(&metrics,coord,map.level[0].,range,bounds,DIMENSION,&MetricCoefficientsFunc2D); CHKERR_PRNT("Metrics computation failed");
	
//	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = %d: metrics.ni = %d, metrics.nj = %d\n",rank,metrics.ni,metrics.nj);
//	for (int i=range[0];i<range[1];i++) {
//		PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = %d: (%d,%d): metrics[%d]: %f %f %f %f %f\n",rank,isGLOBALtoGRID(0,i,0),isGLOBALtoGRID(0,i,1),i,METRICS(i-range[0],0),METRICS(i-range[0],1),METRICS(i-range[0],2),METRICS(i-range[0],3),METRICS(i-range[0],4));
//	}
//	PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
	
//	f = malloc((range[0].end-range[0].start)*sizeof(double));if (f==NULL) ERROR_MSG("malloc failed");
//	GetFuncValues2d(coord, IsGlobalToGrid, f, range);

//	for (int i=range[0];i<range[1];i++) {
//		PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = %d: (%d,%d): f[%d]: %f\n",rank,isGLOBALtoGRID(0,i,0),isGLOBALtoGRID(0,i,1),i-range[0],f[i-range[0]]);
//	}
//	PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
/*	
	if (rank==0) {	
	
	// Memory allocation of RHS, solution and residual
//	f.ni = n[1]-2;
//	f.nj = n[0]-2;
//	f.data = malloc(f.ni*f.nj*sizeof(double));if (f.data==NULL) ERROR_MSG("malloc failed");
	u.ni = n[1]-2;
	u.nj = n[0]-2;
	u.data = malloc(u.ni*u.nj*sizeof(double));if (u.data==NULL) ERROR_MSG("malloc failed");
	rnorm = malloc((numIter+1)*sizeof(double));if (rnorm==NULL) ERROR_MSG("malloc failed");
//	GetFuncValues2d(coord,n,f);
//	UpdateBC(coord,u,n);
	
	}	
	
	// Solver
	prolongStencil2D(&opIH2h, 3, 3);
	restrictStencil2D(&opIh2H, 3, 3);

	MPI_Barrier(PETSC_COMM_WORLD);
	MultigridPetsc(u, metrics, f, opIH2h, opIh2H, IsStencil, IsResStencil, IsProStencil, range, rnorm, levels, n, &numIter);
*/	
/**********************************************************************************/	
/*
	PetscLogStage	stage, stageSolve;
	
//	if (rank==0) printf("Matrix and vector constructions: ");
	MPI_Barrier(PETSC_COMM_WORLD);

	double initAWallTime = MPI_Wtime();
	clock_t initAT = clock();
	PetscLogStageRegister("Setup A", &stage);
	PetscLogStagePush(stage);
	
	A = matrixA(metrics, opIH2h, opIh2H, (n[0]-2), levels);
	MatGetOwnershipRange(A, &rowStart, &rowEnd);
	
	PetscLogStagePop();
	clock_t endAT = clock();
	double endAWallTime = MPI_Wtime();

//	MatView(A, PETSC_VIEWER_STDOUT_WORLD);

	MatCreateVecs(A,&x,&b);
	vecb(&b, f, opIh2H, (n[0]-2), levels);
//	VecView(b, PETSC_VIEWER_STDOUT_WORLD);
	MPI_Barrier(PETSC_COMM_WORLD);
//	if (rank==0) printf("done\n");

//	if (rank==0) printf("Solving...\n");
	MPI_Barrier(PETSC_COMM_WORLD);
	
	KSPCreate(PETSC_COMM_WORLD, &solver);
	KSPSetOperators(solver, A, A);
	KSPGetPC(solver,&pc);
	PCSetType(pc,PCASM);
	KSPMonitorSet(solver, myMonitor, rnorm, NULL);
	KSPSetTolerances(solver, 1.e-7, PETSC_DEFAULT, PETSC_DEFAULT, numIter);
	KSPSetFromOptions(solver);

	double initWallTime = MPI_Wtime();
	clock_t solverInitT = clock();
	PetscLogStageRegister("Solver", &stageSolve);
	PetscLogStagePush(stageSolve);
	
	KSPSolve(solver, b, x);
	
	PetscLogStagePop();
	clock_t solverT = clock();
	double endWallTime = MPI_Wtime();
	
	MPI_Barrier(PETSC_COMM_WORLD);
//	if (rank==0) printf("\nSolver done\n");
	
//	VecView(x,PETSC_VIEWER_STDOUT_WORLD);
	KSPGetIterationNumber(solver, &numIter);
	VecGetArray(x,&px);
	//VecGetOwnershipRange(x, &rowStart, &rowEnd);
	VecGetOwnershipRanges(x,&ranges);
	GetSol(u,px,n,levels,ranges,size,rank);
	VecRestoreArray(x,&px);
	
	MatDestroy(&A); VecDestroy(&b); VecDestroy(&x);
	KSPDestroy(&solver);

	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = [%d]; A construction cputime:        %lf\n",rank,(double)(endAT-initAT)/CLOCKS_PER_SEC);
	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = [%d]; A construction walltime:       %lf\n",rank,endAWallTime-initAWallTime);
	PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = [%d]; Solver cputime:                %lf\n",rank,(double)(solverT-solverInitT)/CLOCKS_PER_SEC);
	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = [%d]; Solver walltime:               %lf\n",rank,endWallTime-initWallTime);
	PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
*/
/**********************************************************************************/	
/*
	if (rank==0) {	
	// Error computation
//	printf("Post-processing: ");

	GetError(coord,n,u,error);
	
	// Output
	solData = fopen("uData.dat","w");
	resData = fopen("rData.dat","w");
	errData = fopen("eData.dat","w");
	
	for(int i=0;i<3;i++){
		printf("\nerror[%d] = %.16e\n",i,error[i]);
		fprintf(errData,"%.16e\n",error[i]);
	}

	for (int i=0;i<u.ni;i++) {
		for (int j=0;j<u.nj;j++) {
			fprintf(solData,"%.16e ",U(i,j));
		}
		fprintf(solData,"\n");
	}
		
	for (int i=0;i<numIter;i++) {
		fprintf(resData,"%.16e ",rnorm[i]);
	}
	fprintf(resData,"\n");

//	printf("done\n");
	}
	
	free(metrics.data);
	if (rank==0) {
	fclose(solData);
	fclose(resData);
	fclose(errData);
//	free2dArray(&coord);
//	free(f);
//	free(u);
//	free(metrics);
//	free(f.data);
	free(u.data);
	free(rnorm);
	}
//	
//	free(ln);
	for (int i=0;i<levels-1;i++) {
		free(IsResStencil[i].data);
	}
	for (int i=0;i<levels;i++) {
		free(IsStencil[i].data);
		free(IsGlobalToGrid[i].data);
		free(IsGridToGlobal[i].data);
	}
	free(range);
	free(IsGlobalToGrid);
	free(IsGridToGlobal);
	free(IsResStencil);
	free(IsStencil);
	free2dArray(&coord);
	free(f);
	free2dArray(&opIH2h);	
	free2dArray(&opIh2H);
*/
//	DeleteIndexMaps(&map);
//	free(range);
//	free(gridId);
	DestroySolver(&solver);
	DestroyAssembly(&assem);
	DestroyOperator(&op);
	DestroyIndices(&indices);
	DestroyMesh(&mesh);
	PetscFinalize();
/*
	if (rank==0) {
	int temp;
	temp = totalUnknowns(n,levels);
	
	printf("=============================================================\n");
	printf("Size:			%d^2\n",n[0]);
	printf("Number of unknowns:	%d\n",temp);
	printf("Number of levels:	%d\n",levels);
	printf("Number of processes:	%d\n",procs);
	printf("Number of iterations:	%d\n",numIter);
	printf("=============================================================\n");
	}
*/
	return 0;
}

PetscErrorCode  myMonitor(KSP ksp, PetscInt n, PetscReal rnormAtn, double *rnorm) {
	//Writes the l2-norm of residual at each iteration to an array
	//
	//rnorm[n] = rnormInstant
	
	int	rank;

	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	if (rank==0) rnorm[n] = rnormAtn;
	return 0;
}


//void CreateStencilIndices(int *n, int levels, IsRange *range, StencilIndices &indices) {
//	// Allocate memory to grid-to-global indices struct
//	
//	int	temp, lg;
//	int	totaln, n0, n1;
//	int	stencilSize = 5;
//
//	(*indices).levels = levels;
//	(*indices).level = malloc((*indices).levels*sizeof(StencilLevel));
//	for (int l=0;l<((*indices).levels);l++) {
//		CreateArrayInt2d(range[i].end-range[i].start, stencilSize, &((*indices).level[l].stencil));
//	}
//	
//		IsResStencil[i-1].ni = range[i].end-range[i].start;
//		IsResStencil[i-1].nj = stencilSize;
//		IsResStencil[i-1].data = malloc(IsResStencil[i-1].ni*IsResStencil[i-1].nj*sizeof(int));
//}
//
//void DeleteStencilIndices( *map) {
//	// Free memory of grid-to-global indices struct
//	
//	for (int l=0;l<(*map).levels;l++) {
//		for (int g=0;g<(*map).level[l].grids;g++) {
//			DeleteArrayInt2d((*map).level[l].grid+g);
////			free((*IsGridToGlobal)[l].grid[g].data);
//		}
//		free((*map).level[l].grid);
//		DeleteArrayInt2d(&((*map).level[l].global));
//	}
//	free((*map).level);
//}

//void Range(int *n, IsRange *gridId, int levels, IsRange *range) {
//	// Computes the range of global indices in this process for all levels	
//	//
//	// range[level].start = Starting global index in level
//	// range[level].end = 1+(ending global index) in level
//	
//	int	remainder, quotient;
//	int	procs, rank;
//	int	totaln, temp, n0, n1;
//	
//	MPI_Comm_size(PETSC_COMM_WORLD, &procs);
//	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
//	
//	for (int l=0;l<levels;l++) {
//		totaln = 0;
//		for (int g=gridId[l].start;g<gridId[l].end;g++) {
//			temp = ipow(2,g);
//			n0 = (n[0]-1)/temp - 1;
//			n1 = (n[0]-1)/temp - 1;
//			totaln = totaln + n0*n1; 
//		}
//		remainder = (totaln)%procs;
//		quotient  = (totaln)/procs;
//		if (rank<remainder) {
//			range[l].start = rank*(quotient + 1);
//			range[l].end = range[l].start + (quotient + 1);
//		}
//		else {
//			range[l].start = rank*quotient + remainder;
//			range[l].end = range[l].start + quotient;
//		}
//	}
//}

//void stencilIndices(IndexMaps map, ArrayInt2d *IsStencil, IsRange *range, IsRange *gridId, int levels) {
//	// Maps a global index of point in a level to global indices of points in its stencil at that level
//	//
//	// IsGridToGlobal[level][i][j] = globalIndex
//	// IsGlobalToGrid[level][globalIndex][0/1] = i/j
//	//
//	// levels - num of multigrid levels
//	// range[level].start - level global index start
//	// range[level].end   - (level global index end + 1) 
//	// IsStencil[level].data[i][j] - global indices of points (j) in stencil at point (i) in a given level
//	
//	int	i0, j0, g0, itemp;
//	int	stencilSize = 5; //Stencil size
//	double	*a, *b;
//	int	aj, bj;
//
//	for (int i=0;i<levels;i++) {
//		IsStencil[i].ni = range[i].end-range[i].start;
//		IsStencil[i].nj = stencilSize;
//		IsStencil[i].data = malloc(IsStencil[i].ni*IsStencil[i].nj*sizeof(int));
//	}
//	
//	for (int l=0;l<map.levels;l++) {
//		a  = map.level[l].global.data;
//		aj = map.level[l].global.nj;
//		for (int i=range[l].start;i<range[l].end;i++) {
//			//i0 - row    - y coord
//			//j0 - column - x coord
//			//A[0]*u(i0-1,j0) + A[1]*u(i0,j0-1) + A[2]*u(i0,j0) + A[3]*u(i0,j0+1) + A[4]*u(i0+1,j0) = f(i0,j0)
////			i0 = isGLOBALtoGRID(l,i,0);
////			j0 = isGLOBALtoGRID(l,i,1);
//			i0 = a[i*aj];
//			j0 = a[i*aj+1];
//			g0 = a[i*aj+2]-gridId[l].start;
//
//			b  = map.level[l].grid[g0].data;
//			bj = map.level[l].grid[g0].nj;
//			itemp = i-range[l].start;
//			if (i0-1<0) {
//				isSTENCIL(l,itemp,0) = -1; 
//			}
//			else {
//				isSTENCIL(l,itemp,0) = b[(i0-1)*bj+j0];//isGRIDtoGLOBAL(l,i0-1,j0); 
//			}
//
//			if (j0-1<0) {
//				isSTENCIL(l,itemp,1) = -1; 
//			}
//			else {
//				isSTENCIL(l,itemp,1) = b[(i0)*bj+j0-1];//isGRIDtoGLOBAL(l,i0,j0-1); 
//			}
//
//			if (j0+1>IsGridToGlobal[l].nj-1) {
//				isSTENCIL(l,itemp,3) = -1; 
//			}
//			else {
//				isSTENCIL(l,itemp,3) = b[(i0)*bj+j0+1];//isGRIDtoGLOBAL(l,i0,j0+1); 
//			}
//
//			if (i0+1>IsGridToGlobal[l].ni-1) {
//				isSTENCIL(l,itemp,4) = -1;
//			}
//			else {
//				isSTENCIL(l,itemp,4) = b[(i0+1)*bj+j0];//isGRIDtoGLOBAL(l,i0+1,j0);
//			}
//			
//			i0 = ipow(2,l)*(i0+1)-1; // fine grid index // knowledge of coarsening strategy used
//			j0 = ipow(2,l)*(j0+1)-1; // fine grid index // knowledge of coarsening strategy used
//			isSTENCIL(l,itemp,2) = b[(i0)*bj+j0];//isGRIDtoGLOBAL(0,i0,j0); // Inserting fine grid global index instead of current level global index
//		}
//	}
//}
//
//void restrictionStencilIndices(ArrayInt2d *IsGlobalToGrid, ArrayInt2d *IsGridToGlobal, ArrayInt2d *IsResStencil, IsRange *range, int levels) {
//	// Maps a global index of point in a coarse level to global indices of points in its restriction stencil (likely from a fine grid level)
//	//
//	// IsGridToGlobal[level][i][j] = globalIndex
//	// IsGlobalToGrid[level][globalIndex][0/1] = i/j
//	//
//	// levels - num of multigrid levels
//	// range[level].start - level global index start
//	// range[level].end   - (level global index end + 1) 
//	// IsResStencil[level].data[i][j] - global indices of points (j) in restriction stencil at point (i) in a given level
//	
//	int	i0, j0, itemp, count;
//	int	stencilSize = 9; //Stencil size
//	int	n = 3; //stencil size per dimension
//
//	for (int i=1;i<levels;i++) {
//		IsResStencil[i-1].ni = range[i].end-range[i].start;
//		IsResStencil[i-1].nj = stencilSize;
//		IsResStencil[i-1].data = malloc(IsResStencil[i-1].ni*IsResStencil[i-1].nj*sizeof(int));
//	}
//	
//	for (int l=1;l<levels;l++) {
//		for (int i=range[l].start;i<range[l].end;i++) {
//			
//			//A[0]*u(i0-1,j0) + A[1]*u(i0,j0-1) + A[2]*u(i0,j0) + A[3]*u(i0,j0+1) + A[4]*u(i0+1,j0) = f(i0,j0)
//			i0 = isGLOBALtoGRID(l,i,0); // l-level grid x-index
//			j0 = isGLOBALtoGRID(l,i,1); // l-level grid y-index
//
//			i0 = 2*(i0+1)-1; // l-level grid x-index // knowledge of coarsening strategy used
//			j0 = 2*(j0+1)-1; // l-level grid y-index // knowledge of coarsening strategy used
//			itemp = i-range[l].start;
//			count = 0;
//			for (int id=-1;id<2;id++) {
//				for (int jd=-1;jd<2;jd++) {
//					isRESSTENCIL(l-1,itemp,count) = isGRIDtoGLOBAL(l-1,i0+id,j0+jd);
//					count = count + 1;
//				}
//			}
//		}
//	}
//}

//void prolongationStencilIndices(ArrayInt2d *IsGlobalToGrid, ArrayInt2d *IsGridToGlobal, ArrayInt2d *IsProStencil, IsRange *range, int levels) {
//	// Maps a global index of point in a coarse level to global indices of points in its prolongation stencil (likely from a fine grid level)
//	//
//	// IsGridToGlobal[level][i][j] = globalIndex
//	// IsGlobalToGrid[level][globalIndex][0/1] = i/j
//	//
//	// levels - num of multigrid levels
//	// range[level].start - level global index start
//	// range[level].end   - (level global index end + 1) 
//	// IsProStencil[level].data[i][j] - global indices of points (j) in prolongation stencil at point (i) in a given level
//	
//	int	i0, j0, itemp, count;
//	int	stencilSize = 9; //Stencil size
//	int	n = 3; //stencil size per dimension
//
//	for (int i=1;i<levels;i++) {
//		IsProStencil[i-1].ni = range[i].end-range[i].start;
//		IsProStencil[i-1].nj = stencilSize;
//		IsProStencil[i-1].data = malloc(IsProStencil[i-1].ni*IsProStencil[i-1].nj*sizeof(int));
//	}
//	
//	for (int l=1;l<levels;l++) {
//		for (int i=range[l].start;i<range[l].end;i++) {
//			
//			//A[0]*u(i0-1,j0) + A[1]*u(i0,j0-1) + A[2]*u(i0,j0) + A[3]*u(i0,j0+1) + A[4]*u(i0+1,j0) = f(i0,j0)
//			i0 = isGLOBALtoGRID(l,i,0); // l-level grid x-index
//			j0 = isGLOBALtoGRID(l,i,1); // l-level grid y-index
//
//			i0 = 2*(i0+1)-1; // l-level grid x-index // knowledge of coarsening strategy used
//			j0 = 2*(j0+1)-1; // l-level grid y-index // knowledge of coarsening strategy used
//			itemp = i-range[l].start;
//			count = 0;
//			for (int id=-1;id<2;id++) {
//				for (int jd=-1;jd<2;jd++) {
//					isRESSTENCIL(l-1,itemp,count) = isGRIDtoGLOBAL(l-1,i0+id,j0+jd);
//					count = count + 1;
//				}
//			}
//		}
//	}
//}

void GetFuncValues2d(double **coord, ArrayInt2d *IsGlobalToGrid, double *f, IsRange *range) {

	// f(x,y) = -2*PI^2*sin(Pi*x)*sin(pi*y)	
//	for (int i=0;i<f.ni;i++) {
//		for (int j=0;j<f.nj;j++) {
//			F(i,j) = FUNC(i+1,j+1);
//		}
//	}
	int	i, j;
	for (int ig=range[0].start;ig<range[0].end;ig++) {
		i = isGLOBALtoGRID(0,ig,0);
		j = isGLOBALtoGRID(0,ig,1);
		f[ig-range[0].start] = FUNC(i+1,j+1);
	}

}

void UpdateBC(double **coord, double *u, int *n) {

	int iend;
	
	for (int j=0;j<n[0];j++) {
		u[j] = SOL(0,j);
	}
	
	iend = n[0]-1;
	for (int i=1;i<n[1]-1;i++) {
		u[i*n[0]] = SOL(i,0);
		u[i*n[0]+iend] = SOL(i,iend);
	}

	iend = n[1]-1;
	for (int j=0;j<n[0];j++) {
		u[iend*n[0]+j] = SOL(iend,j);
	}
	
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

/*
double func(double x, double y) {
	return -2.0*PI*PI*sin(PI*x)*sin(PI*y);
}
*/
void CreateArrayOfIS(int n, int levels, IS *idx) {
	//Gives block indices range for 2D case, where block size is n*n
	
	int first;

	first = 0;
	for (int i=0;i<levels;i++) {
		ISCreateStride(PETSC_COMM_SELF,n*n,first,1,&(idx[i]));
		first	= n*n;
		n	= (n-1)/2;
	}
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

int totalUnknowns(int *n, int levels) {
		
	int length, n0;

	n0 = n[0]-2;
	length=n0*n0;
	for (int i=1;i<levels;i++) {
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
	for (int l=0;l<indices.levels;l++) {
		PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = %d; level: %d: range start = %d, range end = %d\n",rank,l,solver.range[l][0],solver.range[l][1]);
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

void ViewLinSysInfo(Assembly assem, int view) {
	// Prints the info of Assembly struct
	
	int	procs, rank;
	
	MPI_Comm_size(PETSC_COMM_WORLD, &procs);
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	
	for (int l=0;l<assem.levels;l++) {
		PetscPrintf(PETSC_COMM_WORLD,"A[%d]:\n",l);
		if (view == 0) MatView(assem.level[l].A,PETSC_VIEWER_STDOUT_WORLD);
		if (view == 1) MatView(assem.level[l].A,PETSC_VIEWER_DRAW_WORLD);
		VecView(assem.level[l].b,PETSC_VIEWER_STDOUT_WORLD);
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
