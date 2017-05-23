#include "header.h"
#include <time.h>
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
Mat GridTransferMatrix(double **Is, int m, int nh, int nH, int flag);
Mat matrixA(double *As, int n, int levels);
Vec vecb(double **f, int n, int levels);
void GetSol(double **u, double *px, int *n);

int main(int argc, char *argv[]) {
	
	//double	weight=(2.0/3.0);
	int	n[DIMENSION], ierr=0, levels, numIter;
	double	**coord, h[DIMENSION], bounds[DIMENSION*2];
	double	**f, **u, **r, error[3], As[5], *px;//,*rnorm ,**r;
	FILE	*solData, *errData;//, *resData;

	KSP	solver;
	Mat	A;
	Vec	b, x;
	int	iters;
	
	freopen("poisson.in", "r", stdin);
	//freopen("poisson.out", "w", stdout);
	freopen("poisson.err", "w", stderr);
	
	PetscInitialize(&argc, &argv, 0, 0);
	//printf("Enter the no .of points in each dimension = ");
	scanf("%d",n); // unTotal is used temporarily
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
	A = matrixA(As,(n[0]-2),levels);
	b = vecb(f,(n[0]-2),levels);
	VecDuplicate(b, &x);
	KSPCreate(PETSC_COMM_WORLD, &solver);
	KSPSetOperators(solver, A, A);
	KSPSetTolerances(solver, 1.e-15, 1.e-50, PETSC_DEFAULT, PETSC_DEFAULT);
	KSPSetFromOptions(solver);
	KSPSolve(solver, b, x);
	KSPGetIterationNumber(solver, &iters);
	VecGetArray(x,&px);
	GetSol(u,px,n);
/*
	r=0;
	for (int i=1;i<n[1]-1;i++) {
		for (int j=1;j<n[0]-1;j++) {
			VecGetValues(x,1,r,u[i][j])
			r = r+1;
		}
	}
*/
	MatDestroy(&A); VecDestroy(&b); VecDestroy(&x);
	KSPDestroy(&solver);
	PetscFinalize();

	clock_t solverT = clock();
	
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
	
	printf("Total time:             %lf\n",(double)(ppT-begin)/CLOCKS_PER_SEC);
	printf("Memory allocation time: %lf\n",(double)(memT-begin)/CLOCKS_PER_SEC);
	printf("Meshing time:           %lf\n",(double)(meshT-memT)/CLOCKS_PER_SEC);
	printf("Initialization time:    %lf\n",(double)(initT-meshT)/CLOCKS_PER_SEC);
	printf("Solver time:            %lf\n",(double)(solverT-initT)/CLOCKS_PER_SEC);
	printf("Post processing time:   %lf\n",(double)(ppT-solverT)/CLOCKS_PER_SEC);
	
	fclose(solData);
	//fclose(resData);
	fclose(errData);
	free2dArray(&coord);
	free2dArray(&f);
	free2dArray(&u);
	//free(rnorm);
	//free(px);
	
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
Mat matrixA(double *As, int n, int levels) {
	Mat	A, IH2h, Ih2H;
	int	r, localr, rowStart, rowEnd, TotalRows, scale;
	double	**opIH2h, **opIh2H;
	int	m = 3, ierr;

	//double	h, invh2;

	//h	= 1.0/(n+1);
	//invh2	= 1.0/(h*h);
	ierr = malloc2d(&opIH2h,m,m); CHKERR_PRNT("malloc failed");
	ierr = malloc2d(&opIh2H,m,m); CHKERR_PRNT("malloc failed");
	for (int lj=0;lj<3;lj++) {
 		opIH2h[0][lj]= 0.5 - 0.25*fabs(1-lj);
 		opIH2h[1][lj]= 1.0 - 0.5*fabs(1-lj);
 		opIH2h[2][lj]= 0.5 - 0.25*fabs(1-lj);
	}


	for (int lj=0;lj<3;lj++) {
 		opIh2H[0][lj]= 0.0;
 		opIh2H[1][lj]= 0.0;
 		opIh2H[2][lj]= 0.0;
	}
	opIh2H[1][1] = 1.0;

	TotalRows = ((n+1)*(n+1)*(ipow(4,levels)-1))/(3*ipow(4,levels-1))-(2*(n+1)*(ipow(2,levels)-1))/(ipow(2,levels-1))+levels;
	printf("TotalRows = %d\n",TotalRows);

	Ih2H =  GridTransferMatrix(opIh2H, 3, n, (n-1)/2, 0);
	IH2h =  GridTransferMatrix(opIH2h, 3, n, (n-1)/2, 1);
	MatCreate(PETSC_COMM_WORLD, &A);
	MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, TotalRows, TotalRows);
	MatSetFromOptions(A);
	MatSetUp(A);
	//MatMPIAIJSetPreallocation(A,5,NULL,5,NULL);
	MatGetOwnershipRange(A, &rowStart, &rowEnd);
	scale = 1;
	r = 0;
	for (int l=0;l<levels;l++) {
		rowEnd = n*n;
		for (localr=0; localr<rowEnd; localr++) {
			//i = r%n; j = r/n;
			if (localr-n>=0) {
				MatSetValue(A, r, r-n, As[0]/scale, INSERT_VALUES);
			}
			if (localr-1>=0 && localr%n!=0) {
				MatSetValue(A, r, r-1, As[1]/scale, INSERT_VALUES); 
			}
			MatSetValue(A, r, r, As[2]/scale, INSERT_VALUES);
			if (localr+1<=n*n-1 && (localr+1)%n!=0) {
				MatSetValue(A, r, r+1, As[3]/scale, INSERT_VALUES);
			}
			if (localr+n<=n*n-1) {
				MatSetValue(A, r, r+n, As[4]/scale, INSERT_VALUES);
			}
			r = r+1;
		}
		//rowStart = rowEnd;
		n = (n-1)/2;
		scale = scale*4;
	}
	MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
	MatDestroy(&IH2h);
	MatDestroy(&Ih2H);

	free2dArray(&opIH2h);	
	free2dArray(&opIh2H);	
	//MatView(A,PETSC_VIEWER_STDOUT_WORLD);
	return A;
}

Mat GridTransferMatrix(double **Is, int m, int nh, int nH, int flag) {
//	Is	- stencil wise grid transfer operator of size m*m
//	nh	- number of unknowns per dimension in fine grid "h"
//	nH	- number of unknowns per dimension in coarse grid "H"
//	flag	- "0" for "Restriction"; "non-zero" for "Interpolation/Prolongation"
	
	Mat	matI;
	int	rowStart, rowEnd, colStart, colEnd;

	MatCreate(PETSC_COMM_WORLD, &matI);
	if (flag == 0) {
		MatSetSizes(matI, PETSC_DECIDE, PETSC_DECIDE, nH*nH, nh*nh);
	} else {
		MatSetSizes(matI, PETSC_DECIDE, PETSC_DECIDE, nh*nh, nH*nH);
	}
	MatSetFromOptions(matI);
	MatSetUp(matI);
	//MatMPIAIJSetPreallocation(A,5,NULL,5,NULL);
	//MatGetOwnershipRange(A, &rowStart, &rowEnd);
	for (int bj=0;bj<nH;bj++) {
		colStart = bj*nH;
		colEnd   = colStart+nH;
		rowStart = (bj*nh)*((m+1)/2);
		for (int bi=0;bi<m;bi++) {
			for (int j=colStart;j<colEnd;j++) {
				rowEnd  = rowStart + m;
				for (int i=rowStart;i<rowEnd;i++) {
					if (flag==0 && Is[bi][i-rowStart]!=0.0) {
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
	MatView(matI, PETSC_VIEWER_STDOUT_WORLD);
	return matI;
}

Vec vecb(double **f, int n, int levels) {
	Vec	V;
	//int	r;
	int	r, rowStart, rowEnd, TotalRows;//, i, j;
	int	skip;
	//double	h;

	TotalRows = ((n+1)*(n+1)*(ipow(4,levels)-1))/(3*ipow(4,levels-1))-(2*(n+1)*(ipow(2,levels)-1))/(ipow(2,levels-1))+levels;
	
	VecCreate(PETSC_COMM_WORLD, &V);
	VecSetSizes(V, PETSC_DECIDE, TotalRows);
	VecSetFromOptions(V);
	VecGetOwnershipRange(V, &rowStart, &rowEnd);
	r=0;
	skip = 1;
	for (int l=0;l<levels;l++) {
		for (int i=skip;i<n+1;i=i+skip) {
			for (int j=skip;j<n+1;j=j+skip) {
				VecSetValue(V, r, f[i][j], INSERT_VALUES);
				r = r+1;
			}
		}
		skip = skip*2;
	}
	VecAssemblyBegin(V);
	VecAssemblyEnd(V);

	return V;
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

