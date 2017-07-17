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

void GetFuncValues2d(double **coord, int *n, double **f);
void GetError(double **coord, int *n, double **u, double *error);
void UpdateBC(double **coord, double **u, int *n);
//void OpA(double *A, double *metrics, double *h);
static int ipow(int base, int exp);
int JacobiMalloc(double ***f, double ***u, double ***r, int *n);
int MultigridMalloc(double ***f, double ***u, double ***r, int *n, int levels);
int AsyncMultigridMalloc(double ***f, double ***u, double ***r,int *n, int levels);
void CreateArrayOfIS(int n, int levels, IS *idx);
//void insertSubMatValues(Mat *subA, int nrows, Mat *A, int i, int j);
void prolongStencil2D(double ***IH2h, int m, int n);
void restrictStencil2D(double ***Ih2H, int m, int n);
//void insertSubVecValues(Vec *subV, Vec *V, int i0);
int totalUnknowns(int *n, int levels);
static void GetSol(double **u, double *px, int *n, int levels, const int *ranges, int numProcs, int rank);
//void GetSol(double **u, double *px, int *n);
double TransformFunc(double *bounds, double length, double xi);
void MetricCoefficientsFunc2D(double *metrics, double *bounds, double *lengths, double x, double y);
PetscErrorCode myMonitor(KSP ksp, PetscInt n, PetscReal rnormAtn, double *rnorm);

int main(int argc, char *argv[]) {
	
	//double	weight=(2.0/3.0);
	int	n[DIMENSION], ierr=0, levels, numIter;
	double	**coord, h, bounds[DIMENSION*2], ***metrics;
	double	**f, **u, **r, error[3], As[5], *px, *rnorm;
	double	**opIH2h, **opIh2H;
	FILE	*solData, *errData, *resData;
	
	int	size, rank;
	int	rowStart, rowEnd;
	
	const 	int	*ranges;
	
	KSP	solver;
	PC	pc;
	Mat	A;
	Vec	b, x;
	
	PetscInitialize(&argc, &argv, 0, 0);
	
	MPI_Comm_size(PETSC_COMM_WORLD, &size);
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	
	freopen("poisson.in", "r", stdin);
	freopen("poisson.out", "w", stdout);
	freopen("poisson.err", "w", stderr);
	
//	if (rank==0) printf("Inputs reading and memory allocation: ");
	MPI_Barrier(PETSC_COMM_WORLD);
	//printf("Enter the no .of points in each dimension = ");
	scanf("%d",n);	// unTotal is used temporarily
	//printf("Enter the no .of iterations = ");
	scanf("%d",&numIter);
	//printf("Enter the no .of Multigrid levels = ");
	scanf("%d",&levels);
	
//	clock_t begin = clock();

	for (int i=1;i<DIMENSION;i++) { 
		n[i]  = n[0];      // No. of points in each dimension
	}
	for (int i=0;i<DIMENSION;i++) {
		bounds[i*2] = 0.0;    // Lower bound in each dimension
		bounds[i*2+1] = 1.0;  // Upper bound in each dimension
	}
	
	if (rank==0) {	
	// Memory allocation of RHS, solution and residual
	ierr = JacobiMalloc(&f,&u,&r,n); CHKERR_PRNT("malloc failed");
	//ierr = MultigridMalloc(&f,&u,&r,n,levels); CHKERR_PRNT("malloc failed");
	//ierr = AsyncMultigridMalloc(&f,&u,&r,n,levels); CHKERR_PRNT("malloc failed");
	rnorm = (double *)malloc((numIter+1)*sizeof(double));if (rnorm==NULL) ERROR_MSG("malloc failed");
//	px = (double *)malloc((n[0]-2)*(n[1]-2)*sizeof(double));if (px==NULL) ERROR_MSG("malloc failed");
	
//	clock_t memT = clock();
	
//	printf("done\n");
	
//	printf("Meshing and metrics computation: ");
	// Meshing
//	ierr = UniformMesh(&coord,n,bounds,h,DIMENSION); CHKERR_PRNT("meshing failed");
	ierr = NonUniformMeshY(&coord,n,bounds,&h,DIMENSION,&TransformFunc); CHKERR_PRNT("meshing failed");
	ierr = MetricCoefficients2D(&metrics,coord,n,bounds,DIMENSION,&MetricCoefficientsFunc2D); CHKERR_PRNT("Metrics computation failed");
	
//	clock_t meshT = clock();
//	printf("done\n");
	
//	printf("RHS function values and BCs: ");
	// f values
	GetFuncValues2d(coord,n,f);
	
	// Stencil operator coefficients
//	OpA(As,h);
	
	// Update 'u' with boundary conditions
	UpdateBC(coord,u,n);
	}	
	// Solver
	//Jacobi(u,f,r,As,weight,rnorm,numIter,n); // Weighted Jacobi
	//Multigrid(u,f,r,As,weight,rnorm,levels,n,numIter); // Multigrid V-cycle
	//PMultigrid(u,f,r,As,weight,rnorm,levels,n,numIter);
	//AsyncMultigrid(u,f,r,As,weight,rnorm,n,numIter);
	prolongStencil2D(&opIH2h, 3, 3);
	restrictStencil2D(&opIh2H, 3, 3);

	MPI_Barrier(PETSC_COMM_WORLD);
//	if (rank==0) printf("done\n");
	MultigridPetsc(u, metrics, f, opIH2h, opIh2H, rnorm, levels, n, &numIter);
	
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

	for (int i=0;i<n[1];i++) {
		for (int j=0;j<n[0];j++) {
			fprintf(solData,"%.16e ",u[i][j]);
		}
		fprintf(solData,"\n");
	}
		
	for (int i=0;i<numIter;i++) {
		fprintf(resData,"%.16e ",rnorm[i]);
	}
	fprintf(resData,"\n");

//	printf("done\n");
	}
	
	if (rank==0) {
	fclose(solData);
	fclose(resData);
	fclose(errData);
	free2dArray(&coord);
	free3dArray(&metrics);
	free2dArray(&f);
	free2dArray(&u);
	free2dArray(&r);
	free(rnorm);
//	free(px);
	}
	
	free2dArray(&opIH2h);	
	free2dArray(&opIh2H);
	PetscFinalize();
	
	if (rank==0) {
	int temp;
	temp = totalUnknowns(n,levels);
	
	printf("=============================================================\n");
	printf("Size:			%d^2\n",n[0]);
	printf("Number of unknowns:	%d\n",temp);
	printf("Number of levels:	%d\n",levels);
	printf("Number of processes:	%d\n",size);
	printf("Number of iterations:	%d\n",numIter);
	printf("=============================================================\n");
	}
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

double TransformFunc(double *bounds, double length, double xi) {
	//Transformation function from computational to physical space
	//
	//bounds - lower and upper bounds of physical coordinate
	//length = (bounds[1]-bounds[0])
	//
	//x or y = T(xi)
	
	double val;
	val = bounds[1]-length*(cos(PI*0.5*xi));
//	val = xi;
	return val;
}

void MetricCoefficientsFunc2D(double *metrics, double *bounds, double *lengths, double x, double y) {
	//Computes following metrics at (x,y)
	//
	//metrics[0] = (xi_x)^2 + (xi_y)^2
	//metrics[1] = (eta_x)^2 + (eta_y)^2
	//metrics[2] = (xi_xx) + (xi_yy)
	//metrics[3] = (eta_xx) + (eta_yy)
	//metrics[4] = (xi_x)(eta_x) + (xi_y)(eta_y)
	//
	//bounds[0] - Lower bound of "x"
	//bounds[1] - Upper bound of "x"
	//bounds[2] - Lower bound of "y"
	//bounds[3] - Upper bound of "y"
	//
	//lengths[0] = bounds[1] - bounds[0]
	//lengths[1] = bounds[3] - bounds[2]
	double temp;

	temp = (lengths[1]*lengths[1]-(bounds[3]-y)*(bounds[3]-y));
	metrics[0] = 1.0;
	metrics[1] = 4.0/(PI*PI*temp);
	metrics[2] = 0.0;
	metrics[3] = (-2.0*(bounds[3]-y))/(PI*sqrt(temp*temp*temp)); 
	metrics[4] = 0.0;
//	metrics[0] = 1.0;
//	metrics[1] = 1.0;
//	metrics[2] = 0.0;
//	metrics[3] = 0.0; 
//	metrics[4] = 0.0;
}

void GetFuncValues2d(double **coord, int *n, double **f) {

	// f(x,y) = -2*PI^2*sin(Pi*x)*sin(pi*y)	
	for (int i=0;i<n[1];i++) {
		for (int j=0;j<n[0];j++) {
			f[i][j] = FUNC(i,j);
		}
	}

}

void GetError(double **coord, int *n, double **u, double *error) {
	
	// u(x,y) = sin(Pi*x)*sin(pi*y)	
	double diff;
	error[0] = 0.0;
	error[1] = 0.0;
	error[2] = 0.0;
	for (int i=1;i<n[1]-1;i++) {
		for (int j=1;j<n[0]-1;j++) {
			diff = fabs(u[i][j]-SOL(i,j));
			error[0] = fmax(diff,error[0]);
			error[1] = error[1] + diff;
			error[2] = error[2] + diff*diff;
		}
	}
	error[2] = sqrt(error[2]);
}

void UpdateBC(double **coord, double **u, int *n) {

	int iend;
	
	for (int j=0;j<n[0];j++) {
		u[0][j] = SOL(0,j);
	}
	
	iend = n[0]-1;
	for (int i=1;i<n[1]-1;i++) {
		u[i][0] = SOL(i,0);
		u[i][iend] = SOL(i,iend);
	}

	iend = n[1]-1;
	for (int j=0;j<n[0];j++) {
		u[iend][j] = SOL(iend,j);
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

int JacobiMalloc(double ***f, double ***u, double ***r, int *n) {
	
	int ierr = 0;

	ierr = malloc2d(f,n[1],n[0]); CHKERR_RETURN("malloc failed");
	ierr = malloc2d(u,n[1],n[0]); CHKERR_RETURN("malloc failed");
	ierr = malloc2d(r,n[1],n[0]); CHKERR_RETURN("malloc failed");
	
	for (int i=0;i<n[1];i++) {
		for (int j=0;j<n[0];j++) {
			(*u)[i][j] = 0.0;
			(*f)[i][j] = 0.0;	
			(*r)[i][j] = 0.0;	
		}
	}

	return ierr;
}

int MultigridMalloc(double ***f, double ***u, double ***r, int *n, int levels) {
	
	int TotalRows, n1, n0, *m, k, ierr = 0;

	TotalRows = (2*(n[1]-1)*(ipow(2,levels)-1))/(ipow(2,levels))+levels;
	m = (int *)malloc(TotalRows*sizeof(int)); if (m==NULL) ierr=1;ERROR_RETURN("malloc failed"); 
	k = 0;
	for (int i=0;i<levels;i++) {
		n1 = (n[1]+ipow(2,i)-1)/(ipow(2,i));
		n0 = (n[0]+ipow(2,i)-1)/(ipow(2,i));
		for (int j=0;j<n1;j++) {
			m[k] = n0;
			k = k+1;
		}
	}
	
	ierr = malloc2dY(f,TotalRows,m); CHKERR_RETURN("malloc failed");
	ierr = malloc2dY(u,TotalRows,m); CHKERR_RETURN("malloc failed");
	ierr = malloc2dY(r,TotalRows,m); CHKERR_RETURN("malloc failed");
	
	for (int i=0;i<TotalRows;i++) {
		for (int j=0;j<m[i];j++) {
			(*u)[i][j] = 0.0;
			(*f)[i][j] = 0.0;	
			(*r)[i][j] = 0.0;	
		}
	}
	free(m);
	return ierr;
}

int AsyncMultigridMalloc(double ***f, double ***u, double ***r,int *n, int levels) {
	
	int TotalRows, n1, n0, *m, k, ierr = 0;

	TotalRows = (2*(n[1]-1)*(ipow(2,levels)-1))/(ipow(2,levels))+levels;
	m = (int *)malloc(TotalRows*sizeof(int)); if (m==NULL) ierr=1;ERROR_RETURN("malloc failed"); 
	k = 0;
	for (int i=0;i<levels;i++) {
		n1 = (n[1]+ipow(2,i)-1)/(ipow(2,i));
		n0 = (n[0]+ipow(2,i)-1)/(ipow(2,i));
		for (int j=0;j<n1;j++) {
			m[k] = n0;
			k = k+1;
		}
	}
	
	ierr = malloc2dY(f,TotalRows,m); CHKERR_RETURN("malloc failed");
	ierr = malloc2dY(u,TotalRows,m); CHKERR_RETURN("malloc failed");
	ierr = malloc2dY(r,TotalRows,m); CHKERR_RETURN("malloc failed");
	
	for (int i=0;i<TotalRows;i++) {
		for (int j=0;j<m[i];j++) {
			(*u)[i][j] = 0.0;
			(*f)[i][j] = 0.0;	
			(*r)[i][j] = 0.0;	
		}
	}
	free(m);
	return ierr;
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

void GetSol(double **u, double *px, int *n, int levels, const int *ranges, int numProcs, int rank) {
	
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
				u[i][j] = x[r];
				r = r+1;
			}
		}
		
		free(x);
	}

}

