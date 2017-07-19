#include "solver.h"

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

static void GetSol(double **u, double *px, int *n, int levels, const int *ranges, int numProcs, int rank) {
	
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

void MultigridPetsc(double **u, double ***metrics, double **f, double **opIH2h, double **opIh2H, double *rnorm, int levels, int *fulln, int *m) {

	int	v[2], n[levels];
	Mat	A[levels], prolongMatrix[levels-1], restrictMatrix[levels-1];
	int	iter;
	double	rnormchk, bnorm;
	
	double	*px;
	const	int	*ranges;
	int	size, rank;
	
	KSP	solver[levels];
	PC	pc[levels];
	Vec	r[levels], x[levels], b[levels];//, xbuf[levels];
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
		restrictMatrix[l] =  restrictionMatrixMPI(opIh2H, 3, n[l], n[l+1]);
		prolongMatrix[l] =  prolongationMatrixMPI(opIH2h, 3, n[l], n[l+1]);
//		MatView(restrictMatrix[l], PETSC_VIEWER_STDOUT_WORLD);
//		MatView(prolongMatrix[l], PETSC_VIEWER_STDOUT_WORLD);
	}

	for (int i=0;i<levels;i++) {
		A[i] = levelMatrixA(metrics, n[i], i);
//		MatView(A[i], PETSC_VIEWER_STDOUT_WORLD);
		MatCreateVecs(A[i],&(x[i]),&(r[i]));
		VecDuplicate(r[i],&(b[i]));
//		VecDuplicate(x[i],&(xbuf[i]));
	}
	vecb(&(b[0]),f,opIh2H,n[0],1);
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
//	PetscPrintf(PETSC_COMM_WORLD,"---------------------------| level = %d |------------------------\n",0);
//	KSPView(solver[0],PETSC_VIEWER_STDOUT_WORLD);
//	PetscPrintf(PETSC_COMM_WORLD,"-----------------------------------------------------------------\n");
	
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
//		PetscPrintf(PETSC_COMM_WORLD,"---------------------------| level = %d |------------------------\n",i);
//		KSPView(solver[i],PETSC_VIEWER_STDOUT_WORLD);
//		PetscPrintf(PETSC_COMM_WORLD,"-----------------------------------------------------------------\n");
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
//		PetscPrintf(PETSC_COMM_WORLD,"---------------------------| level = %d |------------------------\n",levels-1);
//		KSPView(solver[levels-1],PETSC_VIEWER_STDOUT_WORLD);
//		PetscPrintf(PETSC_COMM_WORLD,"-----------------------------------------------------------------\n");
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
			KSPBuildResidual(solver[l-1],NULL,NULL,&(r[l-1]));
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
			MatMult(prolongMatrix[l],x[l+1],r[l]);
			VecAXPY(x[l],1.0,r[l]);
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

		KSPBuildResidual(solver[0],NULL,NULL,&(r[0]));
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
		VecDestroy(&(r[i]));
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

