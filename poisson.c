#include "header.h"
#include <time.h>
#include <string.h>
#include "petscksp.h"

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
void OpA(double *A, double *metrics, double *h);
int ipow(int base, int exp);
int JacobiMalloc(double ***f, double ***u, double ***r, int *n);
int MultigridMalloc(double ***f, double ***u, double ***r, int *n, int levels);
int AsyncMultigridMalloc(double ***f, double ***u, double ***r,int *n, int levels);
void CreateArrayOfIS(int n, int levels, IS *idx);
void insertSubMatValues(Mat *subA, int nrows, Mat *A, int i, int j);
void prolongStencil2D(double ***IH2h, int m, int n);
void restrictStencil2D(double ***Ih2H, int m, int n);
void insertSubVecValues(Vec *subV, Vec *V, int i0);
void GetSol(double **u, double *px, int *n, int levels, const int *ranges, int numProcs, int rank);
//void GetSol(double **u, double *px, int *n);
double TransformFunc(double *bounds, double length, double xi);
void MetricCoefficientsFunc2D(double *metrics, double *bounds, double *lengths, double x, double y);

int main(int argc, char *argv[]) {
	
	//double	weight=(2.0/3.0);
	int	n[DIMENSION], ierr=0, levels, numIter;
	double	**coord, h, bounds[DIMENSION*2], ***metrics;
	double	**f, **u, **r, error[3], As[5], *px;//, *rnorm;
	double	**opIH2h, **opIh2H;
	FILE	*solData, *errData;//, *resData;
	
	int	size, rank;
	int	rowStart, rowEnd;
	
	const 	int	*ranges;
	
	PetscLogStage	stage, stageSolve;
	//PetscViewer	viewer;

	KSP	solver;
	PC	pc;
	Mat	A;
	Vec	b, x;
	int	iters;
	
	PetscInitialize(&argc, &argv, 0, 0);
	
	MPI_Comm_size(PETSC_COMM_WORLD, &size);
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
//	if (rank==0) {
	freopen("poisson.in", "r", stdin);
	//freopen("poisson.out", "w", stdout);
//	freopen("petsc.dat", "w", stdout);
//	freopen("poisson.err", "w", stderr);
//	}
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
	//rnorm = (double *)malloc((numIter+1)*sizeof(double));if (rnorm==NULL) ERROR_MSG("malloc failed");
	//px = (double *)malloc((n[0]-2)*(n[1]-2)*sizeof(double));if (px==NULL) ERROR_MSG("malloc failed");
	
//	clock_t memT = clock();
	// Meshing
//	ierr = UniformMesh(&coord,n,bounds,h,DIMENSION); CHKERR_PRNT("meshing failed");
	ierr = NonUniformMeshY(&coord,n,bounds,&h,DIMENSION,&TransformFunc); CHKERR_PRNT("meshing failed");
	ierr = MetricCoefficients2D(&metrics,coord,n,bounds,DIMENSION,&MetricCoefficientsFunc2D); CHKERR_PRNT("Metrics computation failed");
	
//	clock_t meshT = clock();
	
	// f values
	GetFuncValues2d(coord,n,f);
	
	// Stencil operator coefficients
//	OpA(As,h);
	
	// Update 'u' with boundary conditions
	UpdateBC(coord,u,n);
	}	
	clock_t initT = clock();
	// Solver
	//Jacobi(u,f,r,As,weight,rnorm,numIter,n); // Weighted Jacobi
	//Multigrid(u,f,r,As,weight,rnorm,levels,n,numIter); // Multigrid V-cycle
	//PMultigrid(u,f,r,As,weight,rnorm,levels,n,numIter);
	//AsyncMultigrid(u,f,r,As,weight,rnorm,n,numIter);
	prolongStencil2D(&opIH2h, 3, 3);
	restrictStencil2D(&opIh2H, 3, 3);

//	PetscPrintf(PETSC_COMM_SELF,"rank = %d, n = %d, numIter = %d, levels = %d\n",rank,n[0],numIter,levels);
	
	PetscLogStageRegister("Setup A", &stage);
	PetscLogStagePush(stage);
	A = matrixA(metrics, opIH2h, opIh2H, (n[0]-2), levels);
	MatGetOwnershipRange(A, &rowStart, &rowEnd);
//	PetscPrintf(PETSC_COMM_SELF,"rank = %d:A: rowStart = %d, rowEnd = %d\n",rank,rowStart,rowEnd);
	PetscLogStagePop();
//	MatView(A, PETSC_VIEWER_STDOUT_WORLD);
	//A = matrixA(As,(n[0]-2),levels);

	clock_t constrA = clock();
	MatCreateVecs(A,&x,&b);
	//b = vecb(f, opIh2H, (n[0]-2), levels);
	vecb(&b, f, opIh2H, (n[0]-2), levels);
//	VecView(b, PETSC_VIEWER_STDOUT_WORLD);
	//b = vecb(f,(n[0]-2),levels);

	VecGetOwnershipRange(b, &rowStart, &rowEnd);
//	PetscPrintf(PETSC_COMM_SELF,"rank = %d:b: rowStart = %d, rowEnd = %d\n",rank,rowStart,rowEnd);
//	clock_t constrb = clock();

	//VecDuplicate(b, &x);
	KSPCreate(PETSC_COMM_WORLD, &solver);
	KSPSetOperators(solver, A, A);
	KSPGetPC(solver,&pc);
	PCSetType(pc,PCASM);
	KSPSetTolerances(solver, 1.e-7, PETSC_DEFAULT, PETSC_DEFAULT, numIter);
	KSPSetFromOptions(solver);
	//PetscViewerASCIIOpen(PETSC_COMM_WORLD, "petsc.data", &viewer);
	//KSPSetResidualHistory(solver, rnorm, numIter, PETSC_TRUE);

	clock_t solverInitT = clock();

	PetscLogStageRegister("Solver", &stageSolve);
	PetscLogStagePush(stageSolve);
	KSPSolve(solver, b, x);
	PetscLogStagePop();
	
	clock_t solverT = clock();
	//KSPGetResidualHistory(solver, &rnorm, &numIter);
	//KSPView(solver, viewer);
//	VecView(x,PETSC_VIEWER_STDOUT_WORLD);
	KSPGetIterationNumber(solver, &iters);
	VecGetArray(x,&px);
	//VecGetOwnershipRange(x, &rowStart, &rowEnd);
	VecGetOwnershipRanges(x,&ranges);
//	for (int i=0;i<size+1;i++) PetscPrintf(PETSC_COMM_SELF,"ranges[%d] = %d\n",i,ranges[i]);
//	for (int i=0;i<ranges[rank+1]-ranges[rank];i++) PetscPrintf(PETSC_COMM_SELF,"rank = %d; px[%d] = %f\n",rank,i,px[i]);
	GetSol(u,px,n,levels,ranges,size,rank);
	//GetSol(u,px,n,rowStart,rowEnd);
	VecRestoreArray(x,&px);
	
	MatDestroy(&A); VecDestroy(&b); VecDestroy(&x);
	KSPDestroy(&solver);
	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = [%d]; A construction time:        %lf\n",rank,(double)(constrA-initT)/CLOCKS_PER_SEC);
	PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = [%d]; Solver time:                %lf\n",rank,(double)(solverT-solverInitT)/CLOCKS_PER_SEC);
	PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
//	return 0;
//	clock_t solverFinalizeT = clock();
	if (rank==0) {	
	// Error computation
	GetError(coord,n,u,error);
	
	// Output
	solData = fopen("uData.dat","w");
	//resData = fopen("rData.dat","w");
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
	}
/*	}
	for (int i=0;i<numIter+1;i++) {
		fprintf(resData,"%.16e ",rnorm[i]);
	}
	fprintf(resData,"\n");
*/
//	clock_t ppT = clock();
	
//	printf("Total time:                 %lf\n",(double)(ppT-begin)/CLOCKS_PER_SEC);
//	printf("Memory allocation time:     %lf\n",(double)(memT-begin)/CLOCKS_PER_SEC);
//	printf("Meshing time:               %lf\n",(double)(meshT-memT)/CLOCKS_PER_SEC);
//	printf("Initialization time:        %lf\n",(double)(initT-meshT)/CLOCKS_PER_SEC);
//	printf("A construction time:        %lf\n",(double)(constrA-initT)/CLOCKS_PER_SEC);
//	printf("b construction time:        %lf\n",(double)(constrb-constrA)/CLOCKS_PER_SEC);
//	printf("Solver Initialization time: %lf\n",(double)(solverInitT-constrb)/CLOCKS_PER_SEC);
//	printf("Solver time:                %lf\n",(double)(solverT-solverInitT)/CLOCKS_PER_SEC);
//	printf("Solver Finalization time:   %lf\n",(double)(solverFinalizeT-solverT)/CLOCKS_PER_SEC);
//	printf("Post processing time:       %lf\n",(double)(ppT-solverT)/CLOCKS_PER_SEC);
	if (rank==0) {
//	printf("=============================================================\n");
//	printf("Size:			%d^2\n",n[0]);
//	printf("Number of levels:	%d\n",levels);
//	printf("=============================================================\n");
	
	fclose(solData);
	//fclose(resData);
	fclose(errData);
	free2dArray(&coord);
	free3dArray(&metrics);
	free2dArray(&f);
	free2dArray(&u);
	}
	//free(rnorm);
	//free(px);
	
	free2dArray(&opIH2h);	
	free2dArray(&opIh2H);
	PetscFinalize();
	if (rank==0) {
	printf("=============================================================\n");
	printf("Size:			%d^2\n",n[0]);
	printf("Number of unknowns:	%d\n",(((n[0]-2+1)*(n[0]-2+1)*(ipow(4,levels)-1))/(3*ipow(4,levels-1))-(2*(n[0]-2+1)*(ipow(2,levels)-1))/(ipow(2,levels-1))+levels));
	printf("Number of levels:	%d\n",levels);
	printf("Number of processes:	%d\n",size);
	printf("=============================================================\n");
	}
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



void OpA(double *A, double *metrics, double *h) {
	//Computes the coefficients
	//A[0]*u(i,j-1) + A[1]*u(i-1,j) + A[2]*u(i,j) + A[3]*u(i+1,j) + A[4]*u(i,j+1) = f(i,j)
	//
	//metrics[5]	- metrics at a point
	//h[2]		- mesh width in computational domain in each direction
	
	double hy2, hx2;
	
	hx2 = h[0]*h[0];
	hy2 = h[1]*h[1];
	A[0] = (metrics[1]/hy2) - (metrics[3]/(2*h[1]));
	A[1] = (metrics[0]/hx2) - (metrics[2]/(2*h[0]));
	A[2] = -2.0*((metrics[0]/hx2) + (metrics[1]/hy2));
	A[3] = (metrics[0]/hx2) + (metrics[2]/(2*h[0]));
	A[4] = (metrics[1]/hy2) + (metrics[3]/(2*h[1]));
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

void insertSubMatValues(Mat *subA, int nrows, Mat *A, int i0, int j0) {
	//Insert values of sub matrix "subA" in "A"
	//
	//nrows	- number of rows in "subA"
	//ncols - number of columns in "subA"
	//
	//A(i0+subi, j0+subj) = subA(subi,subj)
	
	const	int	*subj;
	const	double	*vals;
		int	ncols;
	
	for (int subi=0; subi<nrows; subi++) {
		MatGetRow(*subA, subi, &ncols, &subj, &vals);
		for (int jcol=0; jcol<ncols; jcol++) {
			if (vals[jcol]!=0.0)	MatSetValue(*A, i0+subi, j0+subj[jcol], vals[jcol], INSERT_VALUES);
		}
		MatRestoreRow(*subA, subi, &ncols, &subj, &vals);
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

void insertSubVecValues(Vec *subV, Vec *V, int i0) {
	//Insert values of sub vector "subV" in a vector "V"
	//
	//V(i0+i) = subV(i)
	
	double	*vals;
	int	n;
	
	VecGetLocalSize(*subV, &n);
	VecGetArray(*subV, &vals);
	for (int i=0; i<n; i++) {
		if (vals[i]!=0.0) VecSetValue(*V, i0+i, vals[i], INSERT_VALUES);
	}
	VecRestoreArray(*subV, &vals);
}

void GetSol(double **u, double *px, int *n, int levels, const int *ranges, int numProcs, int rank) {
	
	int	r;
	
	if (rank!=0) {
		MPI_Send(px, ranges[rank+1]-ranges[rank], MPI_DOUBLE, 0, rank, PETSC_COMM_WORLD);
	}
	else if (rank==0) {
	
		int	length, n0;
		double	*x;
	
		n0 = n[0]-2;
		length = ((n0+1)*(n0+1)*(ipow(4,levels)-1))/(3*ipow(4,levels-1))-(2*(n0+1)*(ipow(2,levels)-1))/(ipow(2,levels-1))+levels;
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

