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

