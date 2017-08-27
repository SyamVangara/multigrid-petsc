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


