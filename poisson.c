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
void OpA(double *A, double *h);
int ipow(int base, int exp);
int JacobiMalloc(double ***f, double ***u, double ***r, int *n);
int MultigridMalloc(double ***f, double ***u, double ***r, int *n, int levels);
int AsyncMultigridMalloc(double ***f, double ***u, double ***r,int *n, int levels);
void CreateArrayOfIS(int n, int levels, IS *idx);
void insertSubMatValues(Mat *subA, int nrows, Mat *A, int i, int j);
Mat GridTransferMatrix(double **Is, int m, int nh, int nH, char *type);
Mat restrictionMatrix(double **Is, int m, int nh, int nH);
Mat prolongationMatrix(double **Is, int m, int nh, int nH);
void prolongStencil2D(double ***IH2h, int m, int n);
void restrictStencil2D(double ***Ih2H, int m, int n);
Mat matrixA(double *As, double **opIH2h, double **opIh2H, int n0, int levels);
void insertSubVecValues(Vec *subV, Vec *V, int i0);
Vec vecb(double **f, double **opIh2H, int n0, int levels);
void GetSol(double **u, double *px, int *n);

int main(int argc, char *argv[]) {
	
	//double	weight=(2.0/3.0);
	int	n[DIMENSION], ierr=0, levels, numIter;
	double	**coord, h[DIMENSION], bounds[DIMENSION*2];
	double	**f, **u, **r, error[3], As[5], *px;//, *rnorm;
	double	**opIH2h, **opIh2H;
	FILE	*solData, *errData;//, *resData;
	//PetscViewer	viewer;

	KSP	solver;
	Mat	A;
	Vec	b, x;
	int	iters;
	
	freopen("poisson.in", "r", stdin);
	//freopen("poisson.out", "w", stdout);
	//freopen("petsc.dat", "w", stdout);
	//freopen("poisson.err", "w", stderr);
	
	PetscInitialize(&argc, &argv, 0, 0);
	//printf("Enter the no .of points in each dimension = ");
	scanf("%d",n);	// unTotal is used temporarily
	//printf("Enter the no .of iterations = ");
	scanf("%d",&numIter);
	//printf("Enter the no .of Multigrid levels = ");
	scanf("%d",&levels);
	
	clock_t begin = clock();

	for (int i=1;i<DIMENSION;i++) { 
		n[i]  = n[0];      // No. of points in each dimension
	}
	for (int i=0;i<DIMENSION;i++) {
		bounds[i*2] = 0.0;    // Lower bound in each dimension
		bounds[i*2+1] = 1.0;  // Upper bound in each dimension
	}
	
	// Memory allocation of RHS, solution and residual
	ierr = JacobiMalloc(&f,&u,&r,n); CHKERR_PRNT("malloc failed");
	//ierr = MultigridMalloc(&f,&u,&r,n,levels); CHKERR_PRNT("malloc failed");
	//ierr = AsyncMultigridMalloc(&f,&u,&r,n,levels); CHKERR_PRNT("malloc failed");
	//rnorm = (double *)malloc((numIter+1)*sizeof(double));if (rnorm==NULL) ERROR_MSG("malloc failed");
	//px = (double *)malloc((n[0]-2)*(n[1]-2)*sizeof(double));if (px==NULL) ERROR_MSG("malloc failed");

	clock_t memT = clock();
	// Meshing
	ierr = UniformMesh(&coord,n,bounds,h,DIMENSION); CHKERR_PRNT("meshing failed");
	
	clock_t meshT = clock();
	
	// f values
	GetFuncValues2d(coord,n,f);
	
	// Stencil operator coefficients
	OpA(As,h);
	
	// Update 'u' with boundary conditions
	UpdateBC(coord,u,n);
	
	clock_t initT = clock();
	// Solver
	//Jacobi(u,f,r,As,weight,rnorm,numIter,n); // Weighted Jacobi
	//Multigrid(u,f,r,As,weight,rnorm,levels,n,numIter); // Multigrid V-cycle
	//PMultigrid(u,f,r,As,weight,rnorm,levels,n,numIter);
	//AsyncMultigrid(u,f,r,As,weight,rnorm,n,numIter);
	prolongStencil2D(&opIH2h, 3, 3);
	restrictStencil2D(&opIh2H, 3, 3);
	A = matrixA(As, opIH2h, opIh2H, (n[0]-2), levels);
	//A = matrixA(As,(n[0]-2),levels);

	clock_t constrA = clock();

	b = vecb(f, opIh2H, (n[0]-2), levels);
	//b = vecb(f,(n[0]-2),levels);

	clock_t constrb = clock();

	VecDuplicate(b, &x);
	KSPCreate(PETSC_COMM_WORLD, &solver);
	KSPSetOperators(solver, A, A);
	KSPSetTolerances(solver, 1.e-16, PETSC_DEFAULT, PETSC_DEFAULT, numIter);
	KSPSetFromOptions(solver);
	//PetscViewerASCIIOpen(PETSC_COMM_WORLD, "petsc.data", &viewer);
	//KSPSetResidualHistory(solver, rnorm, numIter, PETSC_TRUE);

	clock_t solverInitT = clock();

	KSPSolve(solver, b, x);
	
	clock_t solverT = clock();
	
	//KSPGetResidualHistory(solver, &rnorm, &numIter);
	//KSPView(solver, viewer);
	//VecView(x,PETSC_VIEWER_STDOUT_WORLD);
	KSPGetIterationNumber(solver, &iters);
	VecGetArray(x,&px);
	GetSol(u,px,n);
	
	MatDestroy(&A); VecDestroy(&b); VecDestroy(&x);
	KSPDestroy(&solver);
	PetscFinalize();

	clock_t solverFinalizeT = clock();
	
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
/*
	for (int i=0;i<numIter+1;i++) {
		fprintf(resData,"%.16e ",rnorm[i]);
	}
	fprintf(resData,"\n");
*/
	clock_t ppT = clock();
	
	printf("Total time:                 %lf\n",(double)(ppT-begin)/CLOCKS_PER_SEC);
	printf("Memory allocation time:     %lf\n",(double)(memT-begin)/CLOCKS_PER_SEC);
	printf("Meshing time:               %lf\n",(double)(meshT-memT)/CLOCKS_PER_SEC);
	printf("Initialization time:        %lf\n",(double)(initT-meshT)/CLOCKS_PER_SEC);
	printf("A construction time:        %lf\n",(double)(constrA-initT)/CLOCKS_PER_SEC);
	printf("b construction time:        %lf\n",(double)(constrb-constrA)/CLOCKS_PER_SEC);
	printf("Solver Initialization time: %lf\n",(double)(solverInitT-constrb)/CLOCKS_PER_SEC);
	printf("Solver time:                %lf\n",(double)(solverT-solverInitT)/CLOCKS_PER_SEC);
	printf("Solver Finalization time:   %lf\n",(double)(solverFinalizeT-solverT)/CLOCKS_PER_SEC);
	printf("Post processing time:       %lf\n",(double)(ppT-solverT)/CLOCKS_PER_SEC);

	printf("=============================================================\n");
	printf("Size:			%d^2\n",n[0]);
	printf("Number of levels:	%d\n",levels);
	printf("=============================================================\n");
	fclose(solData);
	//fclose(resData);
	fclose(errData);
	free2dArray(&coord);
	free2dArray(&f);
	free2dArray(&u);
	//free(rnorm);
	//free(px);
	
	free2dArray(&opIH2h);	
	free2dArray(&opIh2H);	
	return 0;
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

void OpA(double *A, double *h) {
	
	double hy2, hx2;
	
	hx2 = h[0]*h[0];
	hy2 = h[1]*h[1];
	A[0] = 1/hy2;
	A[1] = 1/hx2;
	A[2] = -2*((1/hx2)+(1/hy2));
	A[3] = 1/hx2;
	A[4] = 1/hy2;

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
	m = (int *)malloc(TotalRows*sizeof(int)); if (m==NULL) ERROR_RETURN("malloc failed"); 
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
	m = (int *)malloc(TotalRows*sizeof(int)); if (m==NULL) ERROR_RETURN("malloc failed"); 
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

Mat matrixA(double *As, double **opIH2h, double **opIh2H, int n0, int levels) {
	// Builds matrix "A" for implicit multigrid correction method
	// As		- Stencilwise differenctial operator
	// opIH2h	- Stencilwise prolongation operator
	// opIh2H	- Stencilwise restriction operator
	// n0		- Number of unknowns per dimension
	// levels	- Number of levels
	
	Mat	A, subA[levels], prolongMatrix[levels-1], restrictMatrix[levels-1];
	Mat	UB[levels-1], LB[levels-1];
	int	n[levels];
	int	rows[levels], cols[levels], scale;//, ncols;
	int	rowStart, rowEnd, blockRowStart[levels], blockColStart[levels];
	
	n[0] = n0;
	blockRowStart[0] = 0;
	blockColStart[0] = 0;
	rows[0] = n[0]*n[0];
	cols[0] = rows[0];
	for (int l=0;l<levels-1;l++) {
		blockRowStart[l+1] = blockRowStart[l] + rows[l];
		blockColStart[l+1] = blockRowStart[l+1];
		n[l+1] = (n[l]-1)/2;
		rows[l+1] = n[l+1]*n[l+1];
		cols[l+1] = rows[l+1];
		restrictMatrix[l] =  restrictionMatrix(opIh2H, 3, n[l], n[l+1]);
		prolongMatrix[l] =  prolongationMatrix(opIH2h, 3, n[l], n[l+1]);
	}
	
	MatCreate(PETSC_COMM_WORLD, &A);
	MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, blockRowStart[levels-1]+rows[levels-1], blockColStart[levels-1]+cols[levels-1]);
	MatSetFromOptions(A);
	MatSetUp(A);
	scale = 1;
	for (int l=0;l<levels;l++) {
		MatCreate(PETSC_COMM_WORLD, &(subA[l]));
		MatSetSizes(subA[l], PETSC_DECIDE, PETSC_DECIDE, rows[l], cols[l]);
		MatSetFromOptions(subA[l]);
		MatSetUp(subA[l]);
		MatGetOwnershipRange(subA[l], &rowStart, &rowEnd);
		for (int i=rowStart; i<rowEnd; i++) {
			if (i-n[l]>=0) {
				MatSetValue(subA[l], i, i-n[l], As[0]/scale, INSERT_VALUES);
			}
			if (i-1>=0 && i%n[l]!=0) {
				MatSetValue(subA[l], i, i-1, As[1]/scale, INSERT_VALUES); 
			}
			MatSetValue(subA[l], i, i, As[2]/scale, INSERT_VALUES);
			if (i+1<=rows[l]-1 && (i+1)%n[l]!=0) {
				MatSetValue(subA[l], i, i+1, As[3]/scale, INSERT_VALUES);
			}
			if (i+n[l]<=rows[l]-1) {
				MatSetValue(subA[l], i, i+n[l], As[4]/scale, INSERT_VALUES);
			}
		}
		MatAssemblyBegin(subA[l],MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(subA[l],MAT_FINAL_ASSEMBLY);
		insertSubMatValues(&(subA[l]), rows[l], &A, blockRowStart[l], blockColStart[l]);
		
		if (l!=levels-1) {
			MatMatMult(subA[l], prolongMatrix[l], MAT_INITIAL_MATRIX, PETSC_DEFAULT, &(UB[l]));
			MatMatMult(restrictMatrix[l], subA[l], MAT_INITIAL_MATRIX, PETSC_DEFAULT, &(LB[l]));
			
			insertSubMatValues(&(UB[l]), rows[l], &A, blockRowStart[l], blockColStart[l+1]);
			insertSubMatValues(&(LB[l]), rows[l+1], &A, blockRowStart[l+1], blockColStart[l]);
			
			for (int b=l+1;b<levels-1;b++) {
				MatMatMult(UB[b-1], prolongMatrix[b], MAT_INITIAL_MATRIX, PETSC_DEFAULT, &(UB[b]));
				MatMatMult(restrictMatrix[b], LB[b-1], MAT_INITIAL_MATRIX, PETSC_DEFAULT, &(LB[b]));
				
				insertSubMatValues(&(UB[b]), rows[l], &A, blockRowStart[l], blockColStart[b+1]);
				insertSubMatValues(&(LB[b]), rows[b+1], &A, blockRowStart[b+1], blockColStart[l]);
				
				MatDestroy(&(UB[b-1]));
				MatDestroy(&(LB[b-1]));
			}
			MatDestroy(&(UB[levels-2]));
			MatDestroy(&(LB[levels-2]));
			MatDestroy(&(prolongMatrix[l]));
			MatDestroy(&(restrictMatrix[l]));
		}
		MatDestroy(&(subA[l]));
		scale = scale*4;
	}
	MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

	return A;
}

Mat GridTransferMatrix(double **Is, int m, int nh, int nH, char *type) {
	// Is	- stencil wise grid transfer operator of size m*m
	// nh	- number of unknowns per dimension in fine grid "h"
	// nH	- number of unknowns per dimension in coarse grid "H"
	// type	- "Restriction" or "Prolongation"
	
	Mat	matI;
	int	rowStart, rowEnd, colStart, colEnd;
	char	res[15] = "Restriction", pro[15] = "Prolongation";
	int	flag;

	MatCreate(PETSC_COMM_WORLD, &matI);
	if (strcmp(type, res)) {
		flag = 0;
		MatSetSizes(matI, PETSC_DECIDE, PETSC_DECIDE, nH*nH, nh*nh);
	} else if (strcmp(type, pro)) {
		flag = 1;
		MatSetSizes(matI, PETSC_DECIDE, PETSC_DECIDE, nh*nh, nH*nH);
	} else {
		printf("'%s' is not a valid grid transfer operation; use '%s' or '%s'\n",type,res,pro);
		return matI;
	}
	MatSetFromOptions(matI);
	MatSetUp(matI);
	for (int bj=0;bj<nH;bj++) {
		colStart = bj*nH;
		colEnd   = colStart+nH;
		rowStart = (bj*nh)*((m+1)/2);
		for (int bi=0;bi<m;bi++) {
			for (int j=colStart;j<colEnd;j++) {
				rowEnd  = rowStart + m;
				for (int i=rowStart;i<rowEnd;i++) {
					if (flag == 0 && Is[bi][i-rowStart]!=0.0) {
						MatSetValue(matI, j, i, Is[bi][i-rowStart], INSERT_VALUES);
					} else if (Is[bi][i-rowStart]!=0.0) {
						MatSetValue(matI, i, j, Is[bi][i-rowStart], INSERT_VALUES);
					}
				}
				rowStart = rowStart + ((m+1)/2);
			}
			rowStart = rowEnd;
		}
	}
	
	MatAssemblyBegin(matI, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(matI, MAT_FINAL_ASSEMBLY);
	//MatView(matI, PETSC_VIEWER_STDOUT_WORLD);
	return matI;
}

Mat restrictionMatrix(double **Is, int m, int nh, int nH) {
	// Is	- stencil wise grid transfer operator of size m*m
	// nh	- number of unknowns per dimension in fine grid "h"
	// nH	- number of unknowns per dimension in coarse grid "H"
	
	Mat	matI;
	int	rowStart, rowEnd, colStart, colEnd;

	MatCreate(PETSC_COMM_WORLD, &matI);
	MatSetSizes(matI, PETSC_DECIDE, PETSC_DECIDE, nH*nH, nh*nh);
	MatSetFromOptions(matI);
	MatSetUp(matI);
	for (int bj=0;bj<nH;bj++) {
		colStart = bj*nH;
		colEnd   = colStart+nH;
		rowStart = (bj*nh)*((m+1)/2);
		for (int bi=0;bi<m;bi++) {
			for (int j=colStart;j<colEnd;j++) {
				rowEnd  = rowStart + m;
				for (int i=rowStart;i<rowEnd;i++) {
					if (Is[bi][i-rowStart]!=0.0) MatSetValue(matI, j, i, Is[bi][i-rowStart], INSERT_VALUES);
				}
				rowStart = rowStart + ((m+1)/2);
			}
			rowStart = rowEnd;
		}
	}
	
	MatAssemblyBegin(matI, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(matI, MAT_FINAL_ASSEMBLY);
	//MatView(matI, PETSC_VIEWER_STDOUT_WORLD);
	return matI;
}

Mat prolongationMatrix(double **Is, int m, int nh, int nH) {
	// Is	- stencil wise grid transfer operator of size m*m
	// nh	- number of unknowns per dimension in fine grid "h"
	// nH	- number of unknowns per dimension in coarse grid "H"
	
	Mat	matI;
	int	rowStart, rowEnd, colStart, colEnd;

	MatCreate(PETSC_COMM_WORLD, &matI);
	MatSetSizes(matI, PETSC_DECIDE, PETSC_DECIDE, nh*nh, nH*nH);
	MatSetFromOptions(matI);
	MatSetUp(matI);
	for (int bj=0;bj<nH;bj++) {
		colStart = bj*nH;
		colEnd   = colStart+nH;
		rowStart = (bj*nh)*((m+1)/2);
		for (int bi=0;bi<m;bi++) {
			for (int j=colStart;j<colEnd;j++) {
				rowEnd  = rowStart + m;
				for (int i=rowStart;i<rowEnd;i++) {
					if (Is[bi][i-rowStart]!=0.0) MatSetValue(matI, i, j, Is[bi][i-rowStart], INSERT_VALUES);
				}
				rowStart = rowStart + ((m+1)/2);
			}
			rowStart = rowEnd;
		}
	}
	
	MatAssemblyBegin(matI, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(matI, MAT_FINAL_ASSEMBLY);
	//MatView(matI, PETSC_VIEWER_STDOUT_WORLD);
	return matI;
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

Vec vecb(double **f, double **opIh2H, int n0, int levels) {
	// Build vector "b" for the implicit multigrid correction method
	// f		- 2D array containing right hand side values at each grid point
	// opIh2H 	- Stencilwise restriction operator
	// n		- Number of unknowns per dimension on finest grid
	// levels	- Number of levels
	
	Vec	b, subb[levels]; 
	Mat	restrictMatrix[levels-1];
	int	r, rowStart, rowEnd, TotalRows;//, i, j;
	int	n[levels];

	TotalRows = ((n0+1)*(n0+1)*(ipow(4,levels)-1))/(3*ipow(4,levels-1))-(2*(n0+1)*(ipow(2,levels)-1))/(ipow(2,levels-1))+levels;
	
	n[0] = n0;
	for (int l=0;l<levels-1;l++) {
		n[l+1] = (n[l]-1)/2;
		restrictMatrix[l] = restrictionMatrix(opIh2H, 3, n[l], n[l+1]);
	}	

	VecCreate(PETSC_COMM_WORLD, &b);
	VecSetSizes(b, PETSC_DECIDE, TotalRows);
	VecSetFromOptions(b);
	VecGetOwnershipRange(b, &rowStart, &rowEnd);
	
	VecCreate(PETSC_COMM_WORLD, &(subb[0]));
	VecSetSizes(subb[0], PETSC_DECIDE, n[0]*n[0]);
	VecSetFromOptions(subb[0]);
	VecGetOwnershipRange(subb[0], &rowStart, &rowEnd);
	r=0;
	for (int i=1;i<n[0]+1;i++) {
		for (int j=1;j<n[0]+1;j++) {
			VecSetValue(subb[0], r, f[i][j], INSERT_VALUES);
			r = r+1;
		}
	}
	VecAssemblyBegin(subb[0]);
	VecAssemblyEnd(subb[0]);
	
	insertSubVecValues(&(subb[0]), &(b), 0);
	
	r=0;
	for (int l=0;l<levels-1;l++) {
		r = r+n[l]*n[l];
		VecCreate(PETSC_COMM_WORLD, &(subb[l+1]));
		VecSetSizes(subb[l+1], PETSC_DECIDE, n[l+1]*n[l+1]);
		VecSetFromOptions(subb[l+1]);
		MatMult(restrictMatrix[l],subb[l],subb[l+1]);
		insertSubVecValues(&(subb[l+1]), &(b), r);
		MatDestroy(&(restrictMatrix[l]));
		VecDestroy(&(subb[l]));
	}
	VecDestroy(&(subb[levels-1]));

	VecAssemblyBegin(b);
	VecAssemblyEnd(b);

	return b;
}

void GetSol(double **u, double *px, int *n) {
	
	int	r;
	r = 0;
	for (int i=1;i<n[1]-1;i++) {
		for (int j=1;j<n[0]-1;j++) {
			u[i][j] = px[r];
			r = r+1;
		}
	}

}

