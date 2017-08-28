void GetFuncValues2d(double **coord, ArrayInt2d *IsGlobalToGrid, double *f, IsRange *range);
void GetError(double **coord, int *n, Array2d u, double *error);
void UpdateBC(double **coord, double *u, int *n);
void CreateArrayOfIS(int n, int levels, IS *idx);
void prolongStencil2D(double ***IH2h, int m, int n);
void restrictStencil2D(double ***Ih2H, int m, int n);
static void GetSol(double *u, double *px, int *n, int levels, const int *ranges, int numProcs, int rank);
PetscErrorCode myMonitor(KSP ksp, PetscInt n, PetscReal rnormAtn, double *rnorm);
void stencilIndices(ArrayInt2d *IsGlobalToGrid, ArrayInt2d *IsGridToGlobal, ArrayInt2d *IsStencil, IsRange *range, int levels);
void restrictionStencilIndices(ArrayInt2d *IsGlobalToGrid, ArrayInt2d *IsGridToGlobal, ArrayInt2d *IsResStencil, IsRange *range, int levels);

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

void MultigridPetsc(Array2d u, Array2d metrics, double *f, double **opIH2h, double **opIh2H, ArrayInt2d *IsStencil, ArrayInt2d *IsResStencil, ArrayInt2d *IsProStencil, IsRange *range, double *rnorm, int levels, int *fulln, int *m) {

	int	v[2], n[levels];
	Mat	A[levels], prolongMatrix[levels-1], restrictMatrix[levels-1];
	int	iter;
	double	rnormchk, bnorm;
	
	double	*px;
	const	int	*ranges;
	int	size, rank;
	
	KSP	solver[levels];
	PC	pc[levels];
	Vec	r[levels], rv[levels], x[levels], b[levels];//, xbuf[levels];
//	Vec	dummy1, dummy2;

	PetscLogStage	stage;
	
	MPI_Comm_size(PETSC_COMM_WORLD, &size);
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	//printf("Enter the number of fine grid sweeps = ");
	
	scanf("%d",v);
	//printf("Enter the number of coarse grid sweeps = ");
	scanf("%d",v+1);

//	printf("rank = %d; v[0] = %d\n",rank,v[0]);
//	printf("rank = %d; v[1] = %d\n",rank,v[1]);

	n[0] = fulln[0]-2;
	for (int i=1;i<levels;i++) {
		n[i] = (n[i-1]-1)/2;
	}

	for (int l=0;l<levels-1;l++) {
		restrictMatrix[l] =  restrictionMatrixMPI(opIh2H, 3, range[l], range[l+1], IsResStencil[l]);
//		restrictMatrix[l] =  restrictionMatrixMPI(opIh2H, 3, n[l], n[l+1]);
		prolongMatrix[l] =  prolongationMatrixMPI(opIH2h, 3, range[l], range[l+1], IsProStencil[l]);
//		MatView(restrictMatrix[l], PETSC_VIEWER_STDOUT_WORLD);
//		MatView(prolongMatrix[l], PETSC_VIEWER_STDOUT_WORLD);
	}

	for (int i=0;i<levels;i++) {
		A[i] = levelMatrixA(metrics, IsStencil[i], n[i], i);
//		A[i] = levelMatrixA(metrics, n[i], i);
//		MatView(A[i], PETSC_VIEWER_STDOUT_WORLD);
		MatCreateVecs(A[i],&(x[i]),&(rv[i]));
		VecDuplicate(rv[i],&(b[i]));
//		VecDuplicate(x[i],&(xbuf[i]));
	}
	levelvecb(&(b[0]),f);
	VecNorm(b[0], NORM_2, &bnorm);
//	printf("rank = %d, bnorm = %f\n", rank, bnorm);
//	VecView(r[0], PETSC_VIEWER_STDOUT_WORLD);
	
	KSPCreate(PETSC_COMM_WORLD, &(solver[0]));
//	KSPSetType(solver[0],KSPGMRES);
	KSPSetType(solver[0],KSPRICHARDSON);
	KSPRichardsonSetScale(solver[0],2.0/3.0);
	KSPSetOperators(solver[0], A[0], A[0]);
	KSPGetPC(solver[0],&(pc[0]));
//	PCSetType(pc[0],PCASM);
//	PCASMSetType(pc[0],PC_ASM_BASIC);
//	PCASMSetOverlap(pc[0],3);
//	PCASMSetTotalSubdomains(pc[0], 32, NULL, NULL);
	PCSetType(pc[0],PCJACOBI);
	KSPSetNormType(solver[0],KSP_NORM_NONE);
	KSPSetTolerances(solver[0], 1.e-7, PETSC_DEFAULT, PETSC_DEFAULT, v[0]);
	
	for (int i=1;i<levels-1;i++) {
		KSPCreate(PETSC_COMM_WORLD, &(solver[i]));
//		KSPSetType(solver[i],KSPGMRES);
		KSPSetType(solver[i],KSPRICHARDSON);
		KSPRichardsonSetScale(solver[i],2.0/3.0);
		KSPSetOperators(solver[i], A[i], A[i]);
		KSPGetPC(solver[i],&(pc[i]));
//		PCSetType(pc[i],PCASM);
//		PCASMSetType(pc[i],PC_ASM_BASIC);
//		PCASMSetOverlap(pc[i],3);
//		PCASMSetTotalSubdomains(pc[i], 32, NULL, NULL);
		PCSetType(pc[i],PCJACOBI);
		KSPSetNormType(solver[i],KSP_NORM_NONE);
		KSPSetTolerances(solver[i], 1.e-7, PETSC_DEFAULT, PETSC_DEFAULT, v[0]);
	}

	if (levels>1) {
		KSPCreate(PETSC_COMM_WORLD, &(solver[levels-1]));
//		KSPSetType(solver[levels-1],KSPGMRES);
		KSPSetType(solver[levels-1],KSPRICHARDSON);
		KSPRichardsonSetScale(solver[levels-1],2.0/3.0);
		KSPSetOperators(solver[levels-1], A[levels-1], A[levels-1]);
		KSPGetPC(solver[levels-1],&(pc[levels-1]));
//		PCSetType(pc[levels-1],PCASM);
//		PCASMSetType(pc[levels-1],PC_ASM_BASIC);
//		PCASMSetOverlap(pc[levels-1],3);
//		PCASMSetTotalSubdomains(pc[levels-1], 32, NULL, NULL);
		PCSetType(pc[levels-1],PCJACOBI);
		KSPSetNormType(solver[levels-1],KSP_NORM_NONE);
		KSPSetTolerances(solver[levels-1], 1.e-7, PETSC_DEFAULT, PETSC_DEFAULT, v[1]);
	}

	double initWallTime = MPI_Wtime();
	clock_t solverInitT = clock();
	PetscLogStageRegister("Solver", &stage);
	PetscLogStagePush(stage);
	iter = 0;
	rnormchk = bnorm;
	if (rank==0) rnorm[0] = 1.0;
	while (iter<*m && 100000000*bnorm > rnormchk && rnormchk > (1.e-7)*bnorm) {
		KSPSolve(solver[0], b[0], x[0]);
		if (iter==0) KSPSetInitialGuessNonzero(solver[0],PETSC_TRUE);
//		KSPBuildResidual(solver[0],NULL,NULL,&(r[0]));
//		VecView(r[0], PETSC_VIEWER_STDOUT_WORLD);
//		MatMult(restrictMatrix[0],r[0],b[1]);
//		VecView(b[1], PETSC_VIEWER_STDOUT_WORLD);
//		for (int l=1;l<levels-1;l++) {
		for (int l=1;l<levels;l++) {
			KSPBuildResidual(solver[l-1],NULL,rv[l-1],&(r[l-1]));
			MatMult(restrictMatrix[l-1],r[l-1],b[l]);
			KSPSolve(solver[l], b[l], x[l]);
			if (l!=levels-1) KSPSetInitialGuessNonzero(solver[l],PETSC_TRUE);
//			KSPBuildResidual(solver[l],NULL,NULL,&(r[l]));
//			MatMult(restrictMatrix[l],r[l],b[l+1]);
		}
//		KSPSolve(solver[levels-1], b[levels-1], x[levels-1]);
//		VecView(x[levels-1], PETSC_VIEWER_STDOUT_WORLD);
		for (int l=levels-2;l>=0;l=l-1) {
//			MatMult(prolongMatrix[l],x[l+1],xbuf[l]);
//			VecAXPY(x[l],1.0,xbuf[l]);
			MatMult(prolongMatrix[l],x[l+1],rv[l]);
			VecAXPY(x[l],1.0,rv[l]);
			KSPSolve(solver[l], b[l], x[l]);
//			KSPBuildResidual(solver[l],NULL,NULL,&(r[l]));
			if (l!=0) KSPSetInitialGuessNonzero(solver[l],PETSC_FALSE);
		}
//		MatMult(prolongMatrix[0],x[1],xbuf[0]);
//		VecView(xbuf[0], PETSC_VIEWER_STDOUT_WORLD);
//		VecView(x[0], PETSC_VIEWER_STDOUT_WORLD);
//		VecAXPY(x[0],1.0,xbuf[0]);
//		VecView(x[0], PETSC_VIEWER_STDOUT_WORLD);
//		KSPSolve(solver[0], b[0], x[0]);

		KSPBuildResidual(solver[0],NULL,rv[0],&(r[0]));
		VecNorm(r[0], NORM_2, &rnormchk);	
//		KSPGetResidualNorm(solver[0],&(rnormchk));
		
//		printf("rank = %d; iter: %d, ul = %f, rnorm = %f, ll = %f\n", rank, iter, 100000000*bnorm, rnormchk, (1.e-7)*bnorm);
		iter = iter + 1;
		if (rank==0) rnorm[iter] = rnormchk/bnorm;
	}
	*m = iter;
	PetscLogStagePop();
	clock_t solverT = clock();
	double endWallTime = MPI_Wtime();

	for (int i=0;i<levels;i++) {
		PetscPrintf(PETSC_COMM_WORLD,"---------------------------| level = %d |------------------------\n",i);
		KSPView(solver[i],PETSC_VIEWER_STDOUT_WORLD);
		PetscPrintf(PETSC_COMM_WORLD,"-----------------------------------------------------------------\n");
	}
//	VecView(x[0], PETSC_VIEWER_STDOUT_WORLD);

//	rnorm[0] = Residual(u,f,r,As,n);
//	for (int i=1;i<m+1;i++) {
//		Vcycle(u,f,r,As,w,v,levels,n);
//		rnorm[i] = Residual(u,f,r,As,n); 
//	}
//	printf("residual = %.16e\n",rnorm[m]);
	VecGetArray(x[0],&px);
	VecGetOwnershipRanges(x[0],&ranges);
	GetSol(u,px,fulln,levels,ranges,size,rank);
	VecRestoreArray(x[0],&px);
	
	for (int i=0;i<levels;i++) {
		MatDestroy(&(A[i]));
		VecDestroy(&(rv[i]));
		VecDestroy(&(x[i]));
		VecDestroy(&(b[i]));
	}
	for (int l=0;l<levels-1;l++) {
		MatDestroy(&(restrictMatrix[l]));
		MatDestroy(&(prolongMatrix[l]));
	}
	
	for (int i=0;i<levels;i++) {
		KSPDestroy(&(solver[i]));
	}
	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = [%d]; Solver cputime:                %lf\n",rank,(double)(solverT-solverInitT)/CLOCKS_PER_SEC);
	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = [%d]; Solver walltime:               %lf\n",rank,endWallTime-initWallTime);
	PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);

}

//Mat matrixA(double *metrics, double **opIH2h, double **opIh2H, int n0, int levels) {
//	// Builds matrix "A" for implicit multigrid correction method
//	// metrics	- metric terms
//	// opIH2h	- Stencilwise prolongation operator
//	// opIh2H	- Stencilwise restriction operator
//	// n0		- Number of unknowns per dimension
//	// levels	- Number of levels
//	
//	Mat	A, subA[levels], prolongMatrix[levels-1], restrictMatrix[levels-1];
//	Mat	UB[levels-1], LB[levels-1];
//	int	n[levels];
//	int	rows[levels], cols[levels];//, ncols;
//	int	rowStart, rowEnd, blockRowStart[levels], blockColStart[levels];
//	double	As[5], h[2];
//
//	int	rank;
//
//	n[0] = n0;
//	blockRowStart[0] = 0;
//	blockColStart[0] = 0;
//	rows[0] = n[0]*n[0];
//	cols[0] = rows[0];
//
//	for (int l=0;l<levels-1;l++) {
//		blockRowStart[l+1] = blockRowStart[l] + rows[l];
//		blockColStart[l+1] = blockRowStart[l+1];
//		n[l+1] = (n[l]-1)/2;
//		rows[l+1] = n[l+1]*n[l+1];
//		cols[l+1] = rows[l+1];
////		restrictMatrix[l] =  restrictionMatrix(opIh2H, 3, n[l], n[l+1]);
////		prolongMatrix[l] =  prolongationMatrix(opIH2h, 3, n[l], n[l+1]);
//	}
//
//	MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, blockRowStart[levels-1]+rows[levels-1], blockColStart[levels-1]+cols[levels-1], 6*levels, PETSC_NULL, 6*levels, PETSC_NULL,&A);
////	MatCreate(PETSC_COMM_WORLD, &A);
////	MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, blockRowStart[levels-1]+rows[levels-1], blockColStart[levels-1]+cols[levels-1]);
////	MatSetFromOptions(A);
////	MatSetUp(A);
//	
//	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
//	if (rank==0) {
//
//	for (int l=0;l<levels-1;l++) {
//		restrictMatrix[l] =  restrictionMatrix(opIh2H, 3, n[l], n[l+1]);
//		prolongMatrix[l] =  prolongationMatrix(opIH2h, 3, n[l], n[l+1]);
//	}
//
//	for (int l=0;l<levels;l++) {
//
//		h[0] = 1.0/(n[l]+1);
//		h[1] = h[0];
//		
//		MatCreateSeqAIJ(PETSC_COMM_SELF, rows[l], cols[l], 5, NULL, &(subA[l]));
////		MatCreate(PETSC_COMM_WORLD, &(subA[l]));
////		MatSetSizes(subA[l], PETSC_DECIDE, PETSC_DECIDE, rows[l], cols[l]);
////		MatSetFromOptions(subA[l]);
////		MatSetUp(subA[l]);
//		MatGetOwnershipRange(subA[l], &rowStart, &rowEnd);
//	//	printf("level: %d\n",l);
//		for (int i=rowStart; i<rowEnd; i++) {
//	//		printf("\ni = %d, im = %d, jm = %d\n",i,ipow(2,l)*((i/n[l])+1)-1,ipow(2,l)*((i%n[l])+1)-1);	
//			OpA(As,(metrics+5*((ipow(2,l)*((i/n[l])+1)-1)*(n0)+(ipow(2,l)*((i%n[l])+1)-1))),h);
//		//	printf("\nrow = %d; As[0] = %f\n",i,As[0]);
//			if (i-n[l]>=0) {
//				MatSetValue(subA[l], i, i-n[l], As[0], INSERT_VALUES);
//			}
//			if (i-1>=0 && i%n[l]!=0) {
//				MatSetValue(subA[l], i, i-1, As[1], INSERT_VALUES); 
//			}
//			MatSetValue(subA[l], i, i, As[2], INSERT_VALUES);
//			if (i+1<=rows[l]-1 && (i+1)%n[l]!=0) {
//				MatSetValue(subA[l], i, i+1, As[3], INSERT_VALUES);
//			}
//			if (i+n[l]<=rows[l]-1) {
//				MatSetValue(subA[l], i, i+n[l], As[4], INSERT_VALUES);
//			}
//		}
//		MatAssemblyBegin(subA[l],MAT_FINAL_ASSEMBLY);
//		MatAssemblyEnd(subA[l],MAT_FINAL_ASSEMBLY);
//		
//		//MatView(subA[l], PETSC_VIEWER_STDOUT_WORLD);
//		insertSubMatValues(&(subA[l]), rows[l], &A, blockRowStart[l], blockColStart[l]);
//		
//		if (l!=levels-1) {
//			MatMatMult(subA[l], prolongMatrix[l], MAT_INITIAL_MATRIX, 1.0, &(UB[l]));
//			MatMatMult(restrictMatrix[l], subA[l], MAT_INITIAL_MATRIX, 1.0, &(LB[l]));
//			
//			insertSubMatValues(&(UB[l]), rows[l], &A, blockRowStart[l], blockColStart[l+1]);
//			insertSubMatValues(&(LB[l]), rows[l+1], &A, blockRowStart[l+1], blockColStart[l]);
//			
//			for (int b=l+1;b<levels-1;b++) {
//				MatMatMult(UB[b-1], prolongMatrix[b], MAT_INITIAL_MATRIX, PETSC_DEFAULT, &(UB[b]));
//				MatMatMult(restrictMatrix[b], LB[b-1], MAT_INITIAL_MATRIX, PETSC_DEFAULT, &(LB[b]));
//				
//				insertSubMatValues(&(UB[b]), rows[l], &A, blockRowStart[l], blockColStart[b+1]);
//				insertSubMatValues(&(LB[b]), rows[b+1], &A, blockRowStart[b+1], blockColStart[l]);
//				
//				MatDestroy(&(UB[b-1]));
//				MatDestroy(&(LB[b-1]));
//			}
//			MatDestroy(&(UB[levels-2]));
//			MatDestroy(&(LB[levels-2]));
//			MatDestroy(&(prolongMatrix[l]));
//			MatDestroy(&(restrictMatrix[l]));
//		}
//		MatDestroy(&(subA[l]));
//
//	}
//	
//	}
//	MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
//	MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
//	return A;
//}


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

Mat restrictionMatrixMPI(double **Is, int m, IsRange rangeh, IsRange rangeH, ArrayInt2d IsResStencil) {
	// Is	- stencil wise grid transfer operator of size m*m
	// lnh	- number of unknowns in the immediate higher level in this process
	//
	// IsResStencil[i][j=0:8]: i-global index at current level; 
	// 			   j-global index of a point in restriction stencil
	
	Mat	matI;
//	int	range[2];
//	int	rank;
	int	*col;

	//MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, nH*nH, nh*nh, &matI);
	MatCreateAIJ(PETSC_COMM_WORLD, IsResStencil.ni, rangeh.end-rangeh.start, PETSC_DETERMINE, PETSC_DETERMINE, 1, PETSC_NULL, 1, PETSC_NULL, &matI);
	
//	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
//	if (rank==0) {
//	MatGetOwnershipRange(matI, range, range+1);
	for (int i=rangeH.start;i<rangeH.end;i++) {
		col = IsResStencil.data+((i-rangeH.start)*IsResStencil.nj);
		for (int j=0;j<IsResStencil.nj;j++) {
			if (Is[j/3][j%3]!=0.0) MatSetValue(matI, i, col[j], Is[j/3][j%3], INSERT_VALUES);
		}
	}
//	for (int bj=0;bj<nH;bj++) {
//		colStart = bj*nH;
//		colEnd   = colStart+nH;
//		rowStart = (bj*nh)*((m+1)/2);
//		for (int bi=0;bi<m;bi++) {
//			for (int j=colStart;j<colEnd;j++) {
//				rowEnd  = rowStart + m;
//				for (int i=rowStart;i<rowEnd;i++) {
//					if (Is[bi][i-rowStart]!=0.0) MatSetValue(matI, j, i, Is[bi][i-rowStart], INSERT_VALUES);
//				}
//				rowStart = rowStart + ((m+1)/2);
//			}
//			rowStart = rowEnd;
//		}
//	}
	
//	}	
	MatAssemblyBegin(matI, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(matI, MAT_FINAL_ASSEMBLY);
	//MatView(matI, PETSC_VIEWER_STDOUT_WORLD);
	return matI;
}

Mat prolongationMatrixMPI(double **Is, int m, IsRange rangeh, IsRange rangeH, ArrayInt2d IsProStencil) {
	// Is	- stencil wise grid transfer operator of size m*m
	// lnh	- number of unknowns in the immediate higher level in this process
	//
	// IsProStencil[i][j=0:8]: i-global index at current level; 
	// 			   j-global index of a point in prolongation stencil
	
	Mat	matI;
//	int	range[2];
//	int	rank;
	int	*row;

	MatCreateAIJ(PETSC_COMM_WORLD, rangeh.end-rangeh.start, IsProStencil.ni, PETSC_DETERMINE, PETSC_DETERMINE, 4, PETSC_NULL, 4, PETSC_NULL, &matI);
	
//	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
//	if (rank==0) {

//	MatGetOwnershipRange(matI, range, range+1);
	for (int i=rangeH.start;i<rangeH.end;i++) {
		row = IsProStencil.data+((i-rangeH.start)*IsProStencil.nj);
		for (int j=0;j<IsProStencil.nj;j++) {
			if (Is[j/3][j%3]!=0.0) MatSetValue(matI, row[j], i, Is[j/3][j%3], INSERT_VALUES);
		}
	}
//	for (int bj=0;bj<nH;bj++) {
//		colStart = bj*nH;
//		colEnd   = colStart+nH;
//		rowStart = (bj*nh)*((m+1)/2);
//		for (int bi=0;bi<m;bi++) {
//			for (int j=colStart;j<colEnd;j++) {
//				rowEnd  = rowStart + m;
//				for (int i=rowStart;i<rowEnd;i++) {
//					if (Is[bi][i-rowStart]!=0.0) MatSetValue(matI, i, j, Is[bi][i-rowStart], INSERT_VALUES);
//				}
//				rowStart = rowStart + ((m+1)/2);
//			}
//			rowStart = rowEnd;
//		}
//	}
	
//	}	
	MatAssemblyBegin(matI, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(matI, MAT_FINAL_ASSEMBLY);
	//MatView(matI, PETSC_VIEWER_STDOUT_WORLD);
	return matI;
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

void Res(Indices *indices, Operator *op, int factor, Assembly *assem) {
	// Assembles the restriction matrices
	// Restriction is only from primary grid of one level to all grids of the next level
	
	int	levels;
	
	levels = assem->levels;
	for (int l=0;l<levels-1;l++) {	
		if (indices->level[l].grids>1) {
			PetscPrintf(PETSC_COMM_WORLD, "For now, only 1 grid per level on all levels except last level is allowed in std multigrid\n");
			return;
		}
	}

	ArrayInt2d	grid0, grid1;
	ArrayInt2d	global0, global1;
	int		g0, g1;
	int		i0, j0, i1, j1;
	int		range0[2], range1[2];
	Mat		*res;
	int		opResni, opResnj;
	double		*opRes;
	double		weight;
	
	res = assem->res;
	for (int l=0;l<levels-1;l++) {
		g0 = indices->level[l].gridId[0];
		g1 = indices->level[l+1].gridId[0];

		opResni = op->res[g1-g0-1].ni;
		opResnj = op->res[g1-g0-1].nj;
		opRes = op->res[g1-g0-1].data;

		global0 = indices->level[l].global;
		global1 = indices->level[l+1].global;
		
		grid0 = indices->level[l].grid[0];
		grid1 = indices->level[l+1].grid[0];
		
		VecGetOwnershipRange(assem->b[l], range0, range0+1);	
		VecGetOwnershipRange(assem->b[l+1], range1, range1+1);	
		
		MatCreateAIJ(PETSC_COMM_WORLD, range1[1]-range1[0], range0[1]-range0[0], PETSC_DETERMINE, PETSC_DETERMINE, 1, PETSC_NULL, 1, PETSC_NULL, res+l);
		for (int row=range1[0];row<range1[1];row++) {
			i1 = global1.data[row*global1.nj];
			j1 = global1.data[row*global1.nj+1];

			i0 = ipow(factor,(g1-g0))*(i1+1)-1-(opResni)/2;
			j0 = ipow(factor,(g1-g0))*(j1+1)-1-(opResnj)/2;
			for (int i=i0;i<i0+opResni;i++) {
				for (int j=j0;j<j0+opResnj;j++) {
					weight = opRes[(i-i0)*opResnj+(j-j0)];
					if (weight != 0.0) MatSetValue(res[l], row, grid0.data[i*grid0.nj+j], weight, ADD_VALUES);
				}
			}
		}
	
		MatAssemblyBegin(res[l], MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(res[l], MAT_FINAL_ASSEMBLY);
	}
}

void Pro(Indices *indices, Operator *op, int factor, Assembly *assem) {
	// Assembles the prolongation matrix for level 1 to 0
	
	int	levels;
	
	levels = assem->levels;
	for (int l=0;l<levels-1;l++) {	
		if (indices->level[l].grids>1) {
			PetscPrintf(PETSC_COMM_WORLD, "For now, only 1 grid per level on all levels except last level is allowed in std multigrid\n");
			return;
		}
	}

	ArrayInt2d	grid0, grid1;
	ArrayInt2d	global0, global1;
	int		g0, g1;
	int		i0, j0, i1, j1;
	int		range0[2], range1[2];
	Mat		*pro;
	int		opProni, opPronj;
	double		*opPro;
	double		weight;
	
	pro = assem->pro;
	for (int l=0;l<levels-1;l++) {
		g0 = indices->level[l].gridId[0];
		g1 = indices->level[l+1].gridId[0];

		opProni = op->pro[g1-g0-1].ni;
		opPronj = op->pro[g1-g0-1].nj;
		opPro = op->pro[g1-g0-1].data;

		global0 = indices->level[l].global;
		global1 = indices->level[l+1].global;
		
		grid0 = indices->level[l].grid[0];
		grid1 = indices->level[l+1].grid[0];
		
		VecGetOwnershipRange(assem->b[l], range0, range0+1);	
		VecGetOwnershipRange(assem->b[l+1], range1, range1+1);	
		
		MatCreateAIJ(PETSC_COMM_WORLD, range0[1]-range0[0], range1[1]-range1[0], PETSC_DETERMINE, PETSC_DETERMINE, 4, PETSC_NULL, 4, PETSC_NULL, pro+l);
		for (int col=range1[0];col<range1[1];col++) {
			i1 = global1.data[col*global1.nj];
			j1 = global1.data[col*global1.nj+1];

			i0 = ipow(factor,(g1-g0))*(i1+1)-1-(opProni)/2;
			j0 = ipow(factor,(g1-g0))*(j1+1)-1-(opPronj)/2;
			for (int i=i0;i<i0+opProni;i++) {
				for (int j=j0;j<j0+opPronj;j++) {
					weight = opPro[(i-i0)*opPronj+(j-j0)];
					if (weight != 0.0) MatSetValue(pro[l], grid0.data[i*grid0.nj+j], col, weight, ADD_VALUES);
				}
			}
		}
	
		MatAssemblyBegin(pro[l], MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(pro[l], MAT_FINAL_ASSEMBLY);
	}
}

static void GetSol(Array2d u, double *px, int *n, int levels, const int *ranges, int numProcs, int rank) {
	
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
		for (int i=0;i<u.ni;i++) {
			for (int j=0;j<u.nj;j++) {
				U(i,j) = x[r];
				r = r+1;
			}
		}
		
		free(x);
	}

}

void UpdateRHS(double *A, double **u, double **r, int *n) {
	
	int iend;
	
	for (int j=1;j<n[0]-1;j++) {
		r[1][j] = r[1][j] - (u[0][j])*A[0];
		//f[0][j-1] = f[0][j-1] - (SOL(0,j))*A[0];
	}
	
	iend = n[0]-3;
	for (int i=1;i<n[1]-1;i++) {
		r[i][1]    = r[i][1]    - (u[i][0])*A[1];
		r[i][iend] = r[i][iend] - (u[i][iend+1])*A[1];
		//f[i-1][0]    = f[i-1][0]    - (SOL(i,0))*A[1];
		//f[i-1][iend] = f[i-1][iend] - (SOL(i,iend+2))*A[1];
	}

	iend = n[1]-3;
	for (int j=1;j<n[0]-1;j++) {
		r[iend][j] = r[iend][j] - (u[iend+1][j])*A[0];
		//f[iend][j-1] = f[iend][j-1] - (SOL(iend+2,j))*A[0];
	}
	

}

double Residual(double **u, double **f, double **r, double *As, int *n) {

	double resLinf;
//	*rnorm = 0.0;
	//*(rnorm+1) = 0.0;
	//*(rnorm+2) = 0.0;
	resLinf = 0.0;
	for (int i=1;i<n[1]-1;i++) {
		for (int j=1;j<n[0]-1;j++) {
			r[i][j] = f[i][j] - (As[0]*u[i-1][j]+As[1]*u[i][j-1]+As[2]*u[i][j]+As[3]*u[i][j+1]+As[4]*u[i+1][j]);
			//res = fabs(r[i][j]);
/*
			*rnorm = fmax(res,*rnorm);
			*(rnorm+1) = *(rnorm+1) + res;
			*(rnorm+2) = *(rnorm+2) + res*res;
*/
			//rnorm = rnorm + res*res;
			resLinf = fmax(resLinf,fabs(r[i][j]));
		}
	}
	//rnorm = sqrt(rnorm);
	return resLinf;
}

void JacobiStep(double **u, double **f, double *As, double w, int *n) {
	
	double temp;
	
	for (int i=1;i<n[1]-1;i++) {
		for (int j=1;j<n[0]-1;j++) {
			
			temp = f[i][j] - (As[0]*u[i-1][j]+As[1]*u[i][j-1]+As[3]*u[i][j+1]+As[4]*u[i+1][j]);
			u[i][j] = (1-w)*u[i][j] + (w/As[2])*temp;
			
			//u[i][j] = u[i][j] + (w/As[2])*r[i][j];
		}
	}
}

void Jacobi(double **u, double **f, double **r, double *As, double w, double *rnorm, int v,int *n) {
	
	int i=0;
	rnorm[0] = Residual(u,f,r,As,n);
	while (i<v && (1.0+0.5*rnorm[i])!=1.0) {
		i = i+1;
		JacobiStep(u,f,As,w,n);
		rnorm[i] = Residual(u,f,r,As,n);
		//GetResidual(*u,*f,As,shift,*r,nt);
		//res = norm(*r,nt);
	}
	printf("residual = %.16e\n",rnorm[i]);

}

void ResidualRestriction(double **f, double **r, int *n) {
	
	for (int i=2;i<n[1]-1;i=i+2) {
		for (int j=2;j<n[0]-1;j=j+2) {
			/*	
			f[n[1]+i/2][j/2] = f[i][j]-(As[0]*u[i-1][j]+As[1]*u[i][j-1]+As[2]*u[i][j]+As[3]*u[i][j+1]+As[4]*u[i+1][j]);
			*/
			f[n[1]+i/2][j/2] = r[i][j];
		}
	}
}

void ErrorCorrection(double **u, int *n, int flag) {
	
	double Iop[3][3];
	int    im, jm;

	for (int lj=0;lj<3;lj++) {
 		Iop[0][lj]= 0.5 - 0.25*fabs(1-lj);
 		Iop[1][lj]= 1.0 - 0.5*fabs(1-lj);
 		Iop[2][lj]= 0.5 - 0.25*fabs(1-lj);
	}
	for (int i=2;i<n[1]-1;i=i+2) {
		for (int j=2;j<n[0]-1;j=j+2) {
			im = n[1]+i/2;
			jm = j/2;
			for (int li=0;li<3;li++) {
				for (int lj=0;lj<3;lj++) {
			 		u[i+li-1][j+lj-1] = u[i+li-1][j+lj-1]+Iop[li][lj]*u[im][jm];
				}
			}
			if (flag==0) u[im][jm] = 0.0;
		}
	}
}

void SweepAndRestrict(double **u, double **f, double **r, double *As, double w, int v,int *n) {
	
	double tempRes;
	for (int i=0;i<v;i++) {
		//JacobiStep(u,f,As,w,n);
		//tempRes = Residual(u,f,r,As,n);
		JacobiStep(u,f,As,w,n);
	}
	tempRes = Residual(u,f,r,As,n);
	ResidualRestriction(f,r,n);
}

void CorrectAndSweep(double **u, double **f, double *As, double w, int v,int *n) {
	
	//double tempRes;
	ErrorCorrection(u,n,0);
	for (int i=0;i<v;i++) {
		//tempRes = Residual(u,f,r,As,n);
		JacobiStep(u,f,As,w,n);
	}

}

void Vcycle(double **u, double **f, double **r, double *As, double w, int *v,int levels,int *n) {
	
	double AsH[levels][5];//, res;
	int    nH[levels][2], nid[levels];
	
	for (int j=0;j<5;j++) {
		AsH[0][j] = As[j];
	}
	
	nH[0][0] = n[0];
	nH[0][1] = n[1];
	nid[0] = 0;
	for (int i=1;i<levels;i++) {
		for (int j=0;j<5;j++) {
			AsH[i][j] = 0.25*AsH[i-1][j];
		}
		nH[i][0] = (nH[i-1][0]+1)/2;
		nH[i][1] = (nH[i-1][1]+1)/2;
		nid[i] = nid[i-1] + nH[i-1][1];
	}
	
	for (int i=0;i<levels-1;i++) {
		SweepAndRestrict((u+nid[i]),(f+nid[i]),(r+nid[i]),AsH[i],w,v[0],nH[i]);
	}
	
	for (int i=0;i<v[1];i++) {
		//res = Residual((u+nid[levels-1]),(f+nid[levels-1]),(r+nid[levels-1]),AsH[levels-1],nH[levels-1]);
		JacobiStep((u+nid[levels-1]),(f+nid[levels-1]),AsH[levels-1],w,nH[levels-1]);
	}
	for (int i=levels-2;i>=0;i=i-1) {
		CorrectAndSweep((u+nid[i]),(f+nid[i]),AsH[i],w,v[0],nH[i]);
	}
}

void Multigrid(double **u, double **f, double **r, double *As, double w, double *rnorm, int levels, int *n,int m) {

	int v[2];

	//printf("Enter the number of fine grid sweeps = ");
	scanf("%d",v);
	//printf("Enter the number of coarse grid sweeps = ");
	scanf("%d",v+1);
	
	rnorm[0] = Residual(u,f,r,As,n);
	for (int i=1;i<m+1;i++) {
		Vcycle(u,f,r,As,w,v,levels,n);
		rnorm[i] = Residual(u,f,r,As,n); 
	}
	printf("residual = %.16e\n",rnorm[m]);

}

void Copy(double **u, double **r, int *n) {
	
	//double temp;
	//int im, jm;
	
	for (int i=1;i<n[1]-1;i++) {
		for (int j=1;j<n[0]-1;j++) {
			
			r[i][j] = u[i][j];
/*
			if ((i%2 == 0) && (j%2 == 0)) {
				im = n[1]+i/2;
				jm = j/2;
				r[im][jm] = u[im][jm];
			}
			//u[i][j] = (1-w)*u[i][j] + (w/As[2])*temp;
*/			
			//u[i][j] = u[i][j] + (w/As[2])*r[i][j];
		}
	}
}

void Subtract(double **u, double **r, int *n) {
	
	//double temp;
	//int im, jm;
	
	for (int i=1;i<n[1]-1;i++) {
		for (int j=1;j<n[0]-1;j++) {
			
			u[i][j] = u[i][j] - r[i][j];
/*
			if ((i%2 == 0) && (j%2 == 0)) {
				im = n[1]+i/2;
				jm = j/2;
				r[im][jm] = u[im][jm];
			}
			//u[i][j] = (1-w)*u[i][j] + (w/As[2])*temp;
*/			
			//u[i][j] = u[i][j] + (w/As[2])*r[i][j];
		}
	}
}

void Pcycle(double **u, double **f, double **r, double *As, double w,int levels,int *n) {
	
	double AsH[levels][5], res;
	int    nH[levels][2], nid[levels], flag;
	
	for (int j=0;j<5;j++) {
		AsH[0][j] = As[j];
	}
	
	nH[0][0] = n[0];
	nH[0][1] = n[1];
	nid[0] = 0;
	for (int i=1;i<levels;i++) {
		for (int j=0;j<5;j++) {
			AsH[i][j] = 0.25*AsH[i-1][j];
		}
		nH[i][0] = (nH[i-1][0]+1)/2;
		nH[i][1] = (nH[i-1][1]+1)/2;
		nid[i] = nid[i-1] + nH[i-1][1];
	}
	
	res = Residual(u,f,r,As,n);
	ResidualRestriction(f,r,n);
	for (int i=1;i<levels-1;i++) {
		res = Residual((u+nid[i]),(f+nid[i]),(r+nid[i]),AsH[i],nH[i]);
		ResidualRestriction((f+nid[i]),(r+nid[i]),nH[i]);
		//SweepAndRestrict((u+nid[i]),(f+nid[i]),(r+nid[i]),AsH[i],w,v[0],nH[i]);
	}
	for (int i=1;i<levels;i++) {
		Copy((u+nid[i]),(r+nid[i]),nH[i]);
	}
	//ErrorCorrection(u+nid[levels-2],nH[levels-2],0);
	flag = 1;
	for (int i=levels-2;i>=0;i=i-1) {
		ErrorCorrection(u+nid[i],nH[i],flag);
		//CorrectAndSweep((u+nid[i]),(f+nid[i]),AsH[i],w,v[0],nH[i]);
	}
	for (int i=1;i<levels;i++) {
		Subtract((u+nid[i]),(r+nid[i]),nH[i]);
	}
	for (int i=0;i<levels;i++) {
		//res = Residual((u+nid[levels-1]),(f+nid[levels-1]),(r+nid[levels-1]),AsH[levels-1],nH[levels-1]);
		JacobiStep((u+nid[i]),(f+nid[i]),AsH[i],w,nH[i]);
	}
}
void PMultigrid(double **u, double **f, double **r, double *As, double w, double *rnorm, int levels, int*n, int m) {
	
	int i, flag;
	double AI[5];	
/*
	ResidualRestriction(f,f,n); // Building f-tilda
	OpAI(As,AI);

	AsyncRres(u,f,r,As,n);
	AsyncCorrection(u,r,AI,n,flag);
	AsyncRestriction(u,r,As,n);
	rnorm[0] = AsyncResNorm(u,r,As,n);
*/
	i = 0;
	rnorm[i] = Residual(u,f,r,As,n);
	while (i<m && (1.0+0.5*rnorm[i])!=1.0) {
		i = i+1;
		Pcycle(u,f,r,As,w,levels,n);
		rnorm[i] = Residual(u,f,r,As,n);
		//GetResidual(*u,*f,As,shift,*r,nt);
		//res = norm(*r,nt);
	}
	printf("residual = %.16e\n",rnorm[i]);
	
}

double L2norm(double *a, int n) {
	
	double result;
	result = a[0]*a[0];
	for (int i=1;i<n;i++) {
		result = result + a[i]*a[i];
	}
	return sqrt(result);
}

double L1Norm(double *a, int n) {
	
	double result;
	result = fabs(a[0]);
	for (int i=1;i<n;i++) {
		result = result + fabs(a[i]);
	}
	return result;
}

double LiNorm(double *a, int n) {
	
	double result;
	result = fabs(a[0]);
	for (int i=1;i<n;i++) {
		result = fmax(result,fabs(a[i]));
	}
	return result;
}

void Initialization(double **u, int *n) {
	
	for (int i=0;i<n[1];i++) {
		for (int j=0;j<n[0];j++) {
			u[i][j] = 0.0;
		}
	}
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

int JacobiMalloc(double ***f, double ***u, double ***r, int *n) {
	
	int ierr = 0;

	ierr = malloc2d(f,n[1],n[0]); CHKERR_RETURN("malloc failed");
	ierr = malloc2d(u,n[1],n[0]); CHKERR_RETURN("malloc failed");
//	ierr = malloc2d(r,n[1],n[0]); CHKERR_RETURN("malloc failed");
	
	for (int i=0;i<n[1];i++) {
		for (int j=0;j<n[0];j++) {
			(*u)[i][j] = 0.0;
			(*f)[i][j] = 0.0;	
//			(*r)[i][j] = 0.0;	
		}
	}

	return ierr;
}


