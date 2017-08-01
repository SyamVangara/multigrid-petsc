#include "matbuild.h"

#define METRICS(i,j,k) (metrics.data[metrics.nk*((i)*metrics.nj+(j))+(k)])
#define PMETRICS(i,j) (metrics.data+metrics.nk*((i)*metrics.nj+(j)))
#define F(i,j) (f.data[((i)*f.nj+(j))])
#define U(i,j) (u.data[((i)*u.nj+(j))])

static int ipow(int base, int exp) {

	int result = 1;
	while (exp) {
		if (exp & 1)
			result *= base;
		exp >>= 1;
		base *= base;
	}
	return result;
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

void OpA(double *A, double *metrics, double *h) {
	//Computes the coefficients
	//
	//i - row    - y coord
	//j - column - x coord
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

Mat levelMatrixA(Array3d metrics, ArrayInt2d IsGlobalToGrid, ArrayInt2d IsGridToGlobal, int n, int l) {
	// Builds matrix "A" at a given multigrid level
	// metrics	- metric terms
	// n		- number of unknowns per dimension
	// l		- level
	
	int	rows, cols, idummy, jdummy;
	double	As[5], h[2];
	Mat	A;
	
	int 	procs, rank;
	int	ln;
	int	range[2];

	rows = n*n;
	cols = rows;
	
	MPI_Comm_size(PETSC_COMM_WORLD, &procs);
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	
	if (rank<(n*n)%procs) {
		ln = (n*n)/procs + 1;
	}
	else {
		ln = (n*n)/procs;
	}

	MatCreateAIJ(PETSC_COMM_WORLD, ln, ln, PETSC_DETERMINE, PETSC_DETERMINE, 5, PETSC_NULL, 4, PETSC_NULL, &A);
//	MatCreateSeqAIJ(PETSC_COMM_SELF, rows, cols, 5, NULL, A);
//	MatGetOwnershipRange(subA[l], &rowStart, &rowEnd);
//	printf("level: %d\n",l);
	MatGetOwnershipRange(A, range, range+1);
//	if (rank==0) {

	h[0] = 1.0/(n+1);
	h[1] = h[0];
	for (int i=range[0]; i<range[1]; i++) {
//		printf("\ni = %d, im = %d, jm = %d\n",i,ipow(2,l)*((i/n[l])+1)-1,ipow(2,l)*((i%n[l])+1)-1);	
//		OpA(As,metrics[ipow(2,l)*((i/n)+1)-1][ipow(2,l)*((i%n)+1)-1],h);
		idummy = ipow(2,l)*((i/n)+1)-1;
		jdummy = ipow(2,l)*((i%n)+1)-1;
		OpA(As,PMETRICS(idummy, jdummy),h);
	//	printf("\nrow = %d; As[0] = %f\n",i,As[0]);
		if (i-n>=0) {
			MatSetValue(A, i, i-n, As[0], INSERT_VALUES);
		}
		if (i-1>=0 && i%n!=0) {
			MatSetValue(A, i, i-1, As[1], INSERT_VALUES); 
		}
		MatSetValue(A, i, i, As[2], INSERT_VALUES);
		if (i+1<=rows-1 && (i+1)%n!=0) {
			MatSetValue(A, i, i+1, As[3], INSERT_VALUES);
		}
		if (i+n<=rows-1) {
			MatSetValue(A, i, i+n, As[4], INSERT_VALUES);
		}
	}
//	for (int i=0; i<rows; i++) {
////		printf("\ni = %d, im = %d, jm = %d\n",i,ipow(2,l)*((i/n[l])+1)-1,ipow(2,l)*((i%n[l])+1)-1);	
////		OpA(As,metrics[ipow(2,l)*((i/n)+1)-1][ipow(2,l)*((i%n)+1)-1],h);
//		idummy = ipow(2,l)*((i/n)+1)-1;
//		jdummy = ipow(2,l)*((i%n)+1)-1;
//		OpA(As,PMETRICS(idummy, jdummy),h);
//	//	printf("\nrow = %d; As[0] = %f\n",i,As[0]);
//		if (i-n>=0) {
//			MatSetValue(A, i, i-n, As[0], INSERT_VALUES);
//		}
//		if (i-1>=0 && i%n!=0) {
//			MatSetValue(A, i, i-1, As[1], INSERT_VALUES); 
//		}
//		MatSetValue(A, i, i, As[2], INSERT_VALUES);
//		if (i+1<=rows-1 && (i+1)%n!=0) {
//			MatSetValue(A, i, i+1, As[3], INSERT_VALUES);
//		}
//		if (i+n<=rows-1) {
//			MatSetValue(A, i, i+n, As[4], INSERT_VALUES);
//		}
//	}

//	}
	MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);

	return A;
}

Mat matrixA(double *metrics, double **opIH2h, double **opIh2H, int n0, int levels) {
	// Builds matrix "A" for implicit multigrid correction method
	// metrics	- metric terms
	// opIH2h	- Stencilwise prolongation operator
	// opIh2H	- Stencilwise restriction operator
	// n0		- Number of unknowns per dimension
	// levels	- Number of levels
	
	Mat	A, subA[levels], prolongMatrix[levels-1], restrictMatrix[levels-1];
	Mat	UB[levels-1], LB[levels-1];
	int	n[levels];
	int	rows[levels], cols[levels];//, ncols;
	int	rowStart, rowEnd, blockRowStart[levels], blockColStart[levels];
	double	As[5], h[2];

	int	rank;

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
//		restrictMatrix[l] =  restrictionMatrix(opIh2H, 3, n[l], n[l+1]);
//		prolongMatrix[l] =  prolongationMatrix(opIH2h, 3, n[l], n[l+1]);
	}

	MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, blockRowStart[levels-1]+rows[levels-1], blockColStart[levels-1]+cols[levels-1], 6*levels, PETSC_NULL, 6*levels, PETSC_NULL,&A);
//	MatCreate(PETSC_COMM_WORLD, &A);
//	MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, blockRowStart[levels-1]+rows[levels-1], blockColStart[levels-1]+cols[levels-1]);
//	MatSetFromOptions(A);
//	MatSetUp(A);
	
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	if (rank==0) {

	for (int l=0;l<levels-1;l++) {
		restrictMatrix[l] =  restrictionMatrix(opIh2H, 3, n[l], n[l+1]);
		prolongMatrix[l] =  prolongationMatrix(opIH2h, 3, n[l], n[l+1]);
	}

	for (int l=0;l<levels;l++) {

		h[0] = 1.0/(n[l]+1);
		h[1] = h[0];
		
		MatCreateSeqAIJ(PETSC_COMM_SELF, rows[l], cols[l], 5, NULL, &(subA[l]));
//		MatCreate(PETSC_COMM_WORLD, &(subA[l]));
//		MatSetSizes(subA[l], PETSC_DECIDE, PETSC_DECIDE, rows[l], cols[l]);
//		MatSetFromOptions(subA[l]);
//		MatSetUp(subA[l]);
		MatGetOwnershipRange(subA[l], &rowStart, &rowEnd);
	//	printf("level: %d\n",l);
		for (int i=rowStart; i<rowEnd; i++) {
	//		printf("\ni = %d, im = %d, jm = %d\n",i,ipow(2,l)*((i/n[l])+1)-1,ipow(2,l)*((i%n[l])+1)-1);	
			OpA(As,(metrics+5*((ipow(2,l)*((i/n[l])+1)-1)*(n0)+(ipow(2,l)*((i%n[l])+1)-1))),h);
		//	printf("\nrow = %d; As[0] = %f\n",i,As[0]);
			if (i-n[l]>=0) {
				MatSetValue(subA[l], i, i-n[l], As[0], INSERT_VALUES);
			}
			if (i-1>=0 && i%n[l]!=0) {
				MatSetValue(subA[l], i, i-1, As[1], INSERT_VALUES); 
			}
			MatSetValue(subA[l], i, i, As[2], INSERT_VALUES);
			if (i+1<=rows[l]-1 && (i+1)%n[l]!=0) {
				MatSetValue(subA[l], i, i+1, As[3], INSERT_VALUES);
			}
			if (i+n[l]<=rows[l]-1) {
				MatSetValue(subA[l], i, i+n[l], As[4], INSERT_VALUES);
			}
		}
		MatAssemblyBegin(subA[l],MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(subA[l],MAT_FINAL_ASSEMBLY);
		
		//MatView(subA[l], PETSC_VIEWER_STDOUT_WORLD);
		insertSubMatValues(&(subA[l]), rows[l], &A, blockRowStart[l], blockColStart[l]);
		
		if (l!=levels-1) {
			MatMatMult(subA[l], prolongMatrix[l], MAT_INITIAL_MATRIX, 1.0, &(UB[l]));
			MatMatMult(restrictMatrix[l], subA[l], MAT_INITIAL_MATRIX, 1.0, &(LB[l]));
			
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

	}
	
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

	//MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, nH*nH, nh*nh, &matI);
	MatCreateSeqAIJ(PETSC_COMM_SELF, nH*nH, nh*nh, 1, NULL, &matI);
//	MatCreate(PETSC_COMM_WORLD, &matI);
//	MatSetType(matI,MATMPIAIJ);
//	MatSetSizes(matI, PETSC_DECIDE, PETSC_DECIDE, nH*nH, nh*nh);
//	MatSetFromOptions(matI);
//	MatSetUp(matI);
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

	MatCreateSeqAIJ(PETSC_COMM_SELF, nh*nh, nH*nH, 4, NULL, &matI);
//	MatCreate(PETSC_COMM_WORLD, &matI);
//	MatSetSizes(matI, PETSC_DECIDE, PETSC_DECIDE, nh*nh, nH*nH);
//	MatSetFromOptions(matI);
//	MatSetUp(matI);
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

Mat restrictionMatrixMPI(double **Is, int m, int nh, int nH) {
	// Is	- stencil wise grid transfer operator of size m*m
	// nh	- number of unknowns per dimension in fine grid "h"
	// nH	- number of unknowns per dimension in coarse grid "H"
	
	Mat	matI;
	int	rowStart, rowEnd, colStart, colEnd;
	int	rank;

	//MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, nH*nH, nh*nh, &matI);
	MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, nH*nH, nh*nh, 1, PETSC_NULL, 1, PETSC_NULL, &matI);
//	MatCreate(PETSC_COMM_WORLD, &matI);
//	MatSetType(matI,MATMPIAIJ);
//	MatSetSizes(matI, PETSC_DECIDE, PETSC_DECIDE, nH*nH, nh*nh);
//	MatSetFromOptions(matI);
//	MatSetUp(matI);
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	if (rank==0) {

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
	
	}	
	MatAssemblyBegin(matI, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(matI, MAT_FINAL_ASSEMBLY);
	//MatView(matI, PETSC_VIEWER_STDOUT_WORLD);
	return matI;
}

Mat prolongationMatrixMPI(double **Is, int m, int nh, int nH) {
	// Is	- stencil wise grid transfer operator of size m*m
	// nh	- number of unknowns per dimension in fine grid "h"
	// nH	- number of unknowns per dimension in coarse grid "H"
	
	Mat	matI;
	int	rowStart, rowEnd, colStart, colEnd;
	int	rank;

	MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, nh*nh, nH*nH, 4, PETSC_NULL, 4, PETSC_NULL, &matI);
//	MatCreateSeqAIJ(PETSC_COMM_SELF, nh*nh, nH*nH, 4, NULL, &matI);
//	MatCreate(PETSC_COMM_WORLD, &matI);
//	MatSetSizes(matI, PETSC_DECIDE, PETSC_DECIDE, nh*nh, nH*nH);
//	MatSetFromOptions(matI);
//	MatSetUp(matI);
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	if (rank==0) {

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
	
	}	
	MatAssemblyBegin(matI, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(matI, MAT_FINAL_ASSEMBLY);
	//MatView(matI, PETSC_VIEWER_STDOUT_WORLD);
	return matI;
}

void levelvecb(Vec *b, double *f) {
	// Build vector "b" for fine grid
	// f - logically 2D array containing right hand side values at each grid point
	
	int	rowStart, rowEnd;

	VecGetOwnershipRange(*b, &rowStart, &rowEnd);
	for (int i=rowStart;i<rowEnd;i++) {
		VecSetValue(*b, i, f[i-rowStart], INSERT_VALUES);
	}
	VecAssemblyBegin(*b);
	VecAssemblyEnd(*b);

	//return b;
}


void vecb(Vec *b, Array2d f, double **opIh2H, int n0, int levels) {
	// Build vector "b" for the implicit multigrid correction method
	// f		- logically 2D array containing right hand side values at each grid point
	// opIh2H 	- Stencilwise restriction operator
	// n		- Number of unknowns per dimension on finest grid
	// levels	- Number of levels
	
	//Vec	b, subb[levels]; 
	Vec	subb[levels]; 
	Mat	restrictMatrix[levels-1];
	int	r, rowStart, rowEnd, TotalRows;//, i, j;
	int	n[levels];

	int	rank;

	TotalRows = ((n0+1)*(n0+1)*(ipow(4,levels)-1))/(3*ipow(4,levels-1))-(2*(n0+1)*(ipow(2,levels)-1))/(ipow(2,levels-1))+levels;
	
	n[0] = n0;
	for (int l=0;l<levels-1;l++) {
		n[l+1] = (n[l]-1)/2;
		//restrictMatrix[l] = restrictionMatrix(opIh2H, 3, n[l], n[l+1]);
	}	

//	VecCreate(PETSC_COMM_WORLD, &b);
//	VecSetSizes(b, PETSC_DECIDE, TotalRows);
//	VecSetFromOptions(b);
//	VecGetOwnershipRange(b, &rowStart, &rowEnd);

	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

	if (rank == 0) {

	for (int l=0;l<levels-1;l++) {
		restrictMatrix[l] = restrictionMatrix(opIh2H, 3, n[l], n[l+1]);
	}	
	VecCreate(PETSC_COMM_SELF, &(subb[0]));
	VecSetSizes(subb[0], PETSC_DECIDE, n[0]*n[0]);
	VecSetFromOptions(subb[0]);
	VecGetOwnershipRange(subb[0], &rowStart, &rowEnd);
	r=0;
	for (int i=0;i<f.ni;i++) {
		for (int j=0;j<f.nj;j++) {
			VecSetValue(subb[0], r, F(i,j) , INSERT_VALUES);
			r = r+1;
		}
	}
	VecAssemblyBegin(subb[0]);
	VecAssemblyEnd(subb[0]);
	
	//insertSubVecValues(&(subb[0]), &(b), 0);
	insertSubVecValues(&(subb[0]), b, 0);
	
	r=0;
	for (int l=0;l<levels-1;l++) {
		r = r+n[l]*n[l];
		VecCreate(PETSC_COMM_SELF, &(subb[l+1]));
		VecSetSizes(subb[l+1], PETSC_DECIDE, n[l+1]*n[l+1]);
		VecSetFromOptions(subb[l+1]);
		MatMult(restrictMatrix[l],subb[l],subb[l+1]);
		//insertSubVecValues(&(subb[l+1]), &(b), r);
		insertSubVecValues(&(subb[l+1]), b, r);
		MatDestroy(&(restrictMatrix[l]));
		VecDestroy(&(subb[l]));
	}
	VecDestroy(&(subb[levels-1]));
	
	}

	VecAssemblyBegin(*b);
	VecAssemblyEnd(*b);

	//return b;
}

