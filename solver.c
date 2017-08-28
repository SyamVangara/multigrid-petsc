#include "solver.h"

#define ERROR_MSG(message) (fprintf(stderr,"Error:%s:%d: %s\n",__FILE__,__LINE__,(message)))
#define ERROR_RETURN(message) {ERROR_MSG(message);return ierr;}
#define CHKERR_PRNT(message) {if(ierr==1) {ERROR_MSG(message);}}
#define CHKERR_RETURN(message) {if(ierr==1) {ERROR_RETURN(message);}}

#define METRICS(i,j,k) (metrics.data[metrics.nk*((i)*metrics.nj+(j))+(k)])
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

void SetUpSolver(Indices *indices, Solver *solver, Cycle cyc) {
	// Allocates memory to Solver struct
		
	solver->rnorm = malloc((solver->numIter+1)*sizeof(double));
	solver->cycle = cyc;
}

void DestroySolver(Solver *solver) {
	// Free the memory in Solver struct
	
	free(solver->rnorm);
}

void SetUpPostProcess(PostProcess *pp) {
	// Allocates memory to PostProcess struct
		
	int	procs, rank;
	
	MPI_Comm_size(PETSC_COMM_WORLD, &procs);
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	
	if (rank == 0) {
		pp->solData = fopen("uData.dat","w");
		pp->resData = fopen("rData.dat","w");
		pp->errData = fopen("eData.dat","w");
	}
}

void DestroyPostProcess(PostProcess *pp) {
	// Free the memory in PostProcess struct
	
	int	procs, rank;
	
	MPI_Comm_size(PETSC_COMM_WORLD, &procs);
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	
	if (rank == 0) {
		fclose(pp->solData);
		fclose(pp->resData);
		fclose(pp->errData);
	}
}

void GetError(Problem *prob, Mesh *mesh, Array2d u1, double *error) {
	
	// u(x,y) = sin(Pi*x)*sin(pi*y)	
	double	diff;
	double	**coord;
	double	sol;
	int	uni, unj;
	double	*u;

	coord = mesh->coord;
	uni   = u1.ni;
	unj   = u1.nj;
	u     = u1.data;
	error[0] = 0.0;
	error[1] = 0.0;
	error[2] = 0.0;
	for (int i=0;i<uni;i++) {
		for (int j=0;j<unj;j++) {
			sol = prob->SOLfunc(coord[0][j+1], coord[1][i+1]);
			diff = fabs(u[i*unj+j]-sol);
			error[0] = fmax(diff,error[0]);
			error[1] = error[1] + diff;
			error[2] = error[2] + diff*diff;
		}
	}
	error[2] = sqrt(error[2]);
}

void GetSol(Indices *indices, Assembly *assem, Array2d u) {
	
	int		r;
	double		*px;
	const	int	*ranges;
	int		gridId;
	int		globalni, globalnj, *global;
	int		i, j, g;
	int		count, localcount;
	double		*buffer;

	int	procs, rank;
	
	MPI_Comm_size(PETSC_COMM_WORLD, &procs);
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	
	VecGetArray(assem->u[0], &px);
	VecGetOwnershipRanges(assem->u[0], &ranges);
	
	globalni = indices->level[0].global.ni;
	globalnj = indices->level[0].global.nj;
	global   = indices->level[0].global.data;
	gridId   = indices->level[0].gridId[0];
	
	if (rank!=0) {
		buffer = malloc((ranges[rank+1]-ranges[rank])*sizeof(double));
		localcount = 0;
		for (int row=ranges[rank];row<ranges[rank+1];row++) {
			g = global[row*globalnj + 2];
			if (g != gridId) continue;
			buffer[localcount] = px[row-ranges[rank]];
			localcount += 1;
		}

		MPI_Send(&localcount, 1, MPI_DOUBLE, 0, rank, PETSC_COMM_WORLD);
		MPI_Send(buffer, localcount, MPI_DOUBLE, 0, rank*(localcount+1), PETSC_COMM_WORLD);
		free(buffer);
	}
	else if (rank==0) {
		int	totalN;
		int	gridni, gridnj;
		int	*grid;
		double		*buffer;
		
		gridni = indices->level[0].grid[0].ni;
		gridnj = indices->level[0].grid[0].nj;
		grid   = indices->level[0].grid[0].data;
		totalN = (gridni*gridnj);
		buffer = malloc(totalN*sizeof(double));
		
		count = 0;
		for (int row=ranges[0];row<ranges[1];row++) {
			g = global[row*globalnj + 2];
			if (g != gridId) continue;
			buffer[count] = px[row-ranges[0]];
			count += 1;
		}

		for (int i=1;i<procs;i++) {
			MPI_Recv(&localcount, 1, MPI_DOUBLE, i, i, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Recv(buffer+count, localcount, MPI_DOUBLE, i, i*(localcount+1), PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
			count += localcount;
		}
		count = 0;
		for (int row=ranges[0];row<ranges[procs];row++) {
			i = global[row*globalnj    ];
			j = global[row*globalnj + 1];
			g = global[row*globalnj + 2];
			if (g != gridId) continue;
			u.data[i*u.nj+j] = buffer[count];
			count += 1;
		}
		free(buffer);
	}
	VecRestoreArray(assem->u[0], &px);

}

void Postprocessing(Problem *prob, Mesh *mesh, Indices *indices, Assembly *assem, Solver *solver, PostProcess *pp) {
	// Computes error and writes data to files
	
	int		rank;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

	Array2d		u;
	
	if (rank==0) CreateArray2d(indices->level[0].grid[0].ni, indices->level[0].grid[0].nj, &u);
	GetSol(indices, assem, u);
//	GetError(prob, mesh, indices, assem, pp->error);
	
	if (rank==0) {	
	
	GetError(prob, mesh, u, pp->error);
	for(int i=0;i<3;i++){
		printf("\nerror[%d] = %.16e\n", i, pp->error[i]);
		fprintf(pp->errData,"%.16e\n",pp->error[i]);
	}

	for (int i=0;i<u.ni;i++) {
		for (int j=0;j<u.nj;j++) {
			fprintf(pp->solData,"%.16e ", u.data[i*u.nj+j]);
		}
		fprintf(pp->solData,"\n");
	}
		
	for (int i=0;i<solver->numIter;i++) {
		fprintf(pp->resData,"%.16e ",solver->rnorm[i]);
	}
	fprintf(pp->resData,"\n");

	DeleteArray2d(&u);
	}
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

void MultigridVcycle(Assembly *assem, Solver *solver) {

	int	iter;
	double	rnormchk, bnorm;
	
	double	*rnorm;
	int	maxIter;

	int	*v;
	int	levels;
	Mat 	*res;
	Mat 	*pro;
	Mat	*A;
	Vec	*b;
	Vec	*u;
	
	int	size, rank;
	
	MPI_Comm_size(PETSC_COMM_WORLD, &size);
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	
	maxIter = solver->numIter;	
	rnorm	= solver->rnorm;
	v	= solver->v;

	levels	= assem->levels;
	res	= assem->res;
	pro	= assem->pro;
	A	= assem->A;
	b	= assem->b;
	u	= assem->u;
	
	KSP	ksp[levels];
	PC	pc[levels];
	Vec	r[levels], rv[levels];//, xbuf[levels];
	
	PetscLogStage	stage;
	
	//printf("Enter the number of fine grid sweeps = ");
//	scanf("%d",v);
	//printf("Enter the number of coarse grid sweeps = ");
//	scanf("%d",v+1);

//	PetscOptionsGetIntArray(NULL, NULL, "-v", v, &vmax, NULL);
	
	for (int i=0;i<levels;i++) {
		VecDuplicate(b[i],&(rv[i]));
	}
	VecNorm(b[0], NORM_2, &bnorm);
	
	KSPCreate(PETSC_COMM_WORLD, &(ksp[0]));
//	KSPSetType(ksp[0],KSPGMRES);
	KSPSetType(ksp[0],KSPRICHARDSON);
//	KSPRichardsonSetScale(ksp[0],2.0/3.0);
	KSPSetOperators(ksp[0], A[0], A[0]);
//	KSPGetPC(ksp[0],&(pc[0]));
//	PCSetType(pc[0],PCASM);
//	PCASMSetType(pc[0],PC_ASM_BASIC);
//	PCASMSetOverlap(pc[0],3);
//	PCASMSetTotalSubdomains(pc[0], 32, NULL, NULL);
//	PCSetType(pc[0],PCJACOBI);
	KSPSetNormType(ksp[0],KSP_NORM_NONE);
	KSPSetTolerances(ksp[0], 1.e-7, PETSC_DEFAULT, PETSC_DEFAULT, v[0]);
	KSPSetFromOptions(ksp[0]);
	
	for (int i=1;i<levels-1;i++) {
		KSPCreate(PETSC_COMM_WORLD, &(ksp[i]));
//		KSPSetType(ksp[i],KSPGMRES);
		KSPSetType(ksp[i],KSPRICHARDSON);
//		KSPRichardsonSetScale(ksp[i],2.0/3.0);
		KSPSetOperators(ksp[i], A[i], A[i]);
//		KSPGetPC(ksp[i],&(pc[i]));
//		PCSetType(pc[i],PCASM);
//		PCASMSetType(pc[i],PC_ASM_BASIC);
//		PCASMSetOverlap(pc[i],3);
//		PCASMSetTotalSubdomains(pc[i], 32, NULL, NULL);
//		PCSetType(pc[i],PCJACOBI);
		KSPSetNormType(ksp[i],KSP_NORM_NONE);
		KSPSetTolerances(ksp[i], 1.e-7, PETSC_DEFAULT, PETSC_DEFAULT, v[0]);
		KSPSetFromOptions(ksp[i]);
	}

	if (levels>1) {
		KSPCreate(PETSC_COMM_WORLD, &(ksp[levels-1]));
//		KSPSetType(ksp[levels-1],KSPGMRES);
		KSPSetType(ksp[levels-1],KSPRICHARDSON);
//		KSPRichardsonSetScale(ksp[levels-1],2.0/3.0);
		KSPSetOperators(ksp[levels-1], A[levels-1], A[levels-1]);
//		KSPGetPC(ksp[levels-1],&(pc[levels-1]));
//		PCSetType(pc[levels-1],PCASM);
//		PCASMSetType(pc[levels-1],PC_ASM_BASIC);
//		PCASMSetOverlap(pc[levels-1],3);
//		PCASMSetTotalSubdomains(pc[levels-1], 32, NULL, NULL);
//		PCSetType(pc[levels-1],PCJACOBI);
		KSPSetNormType(ksp[levels-1],KSP_NORM_NONE);
		KSPSetTolerances(ksp[levels-1], 1.e-7, PETSC_DEFAULT, PETSC_DEFAULT, v[1]);
		KSPSetFromOptions(ksp[levels-1]);
	}

	double initWallTime = MPI_Wtime();
	clock_t solverInitT = clock();
	PetscLogStageRegister("Solver", &stage);
	PetscLogStagePush(stage);
	iter = 0;
	rnormchk = bnorm;
	if (rank==0) rnorm[0] = 1.0;
	while (iter<maxIter && 100000000*bnorm > rnormchk && rnormchk > (1.e-7)*bnorm) {
		KSPSolve(ksp[0], b[0], u[0]);
		if (iter==0) KSPSetInitialGuessNonzero(ksp[0],PETSC_TRUE);
		for (int l=1;l<levels;l++) {
			KSPBuildResidual(ksp[l-1],NULL,rv[l-1],&(r[l-1]));
			MatMult(res[l-1],r[l-1],b[l]);
			KSPSolve(ksp[l], b[l], u[l]);
			if (l!=levels-1) KSPSetInitialGuessNonzero(ksp[l],PETSC_TRUE);
		}
		for (int l=levels-2;l>=0;l=l-1) {
			MatMult(pro[l],u[l+1],rv[l]);
			VecAXPY(u[l],1.0,rv[l]);
			KSPSolve(ksp[l], b[l], u[l]);
			if (l!=0) KSPSetInitialGuessNonzero(ksp[l],PETSC_FALSE);
		}
		KSPBuildResidual(ksp[0],NULL,rv[0],&(r[0]));
		VecNorm(r[0], NORM_2, &rnormchk);	
		iter = iter + 1;
		if (rank==0) rnorm[iter] = rnormchk/bnorm;
	}
	PetscLogStagePop();
	clock_t solverT = clock();
	double endWallTime = MPI_Wtime();
	solver->numIter = iter;

	for (int i=0;i<levels;i++) {
		PetscPrintf(PETSC_COMM_WORLD,"---------------------------| level = %d |------------------------\n",i);
		KSPView(ksp[i],PETSC_VIEWER_STDOUT_WORLD);
		PetscPrintf(PETSC_COMM_WORLD,"-----------------------------------------------------------------\n");
	}
	for (int i=0;i<levels;i++) {
		VecDestroy(&(rv[i]));
	}
	for (int i=0;i<levels;i++) {
		KSPDestroy(&(ksp[i]));
	}
	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = [%d]; Solver cputime:                %lf\n",rank,(double)(solverT-solverInitT)/CLOCKS_PER_SEC);
	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = [%d]; Solver walltime:               %lf\n",rank,endWallTime-initWallTime);
	PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);

}

void MultigridIcycle(Assembly *assem, Solver *solver) {

	Mat	*A;
	Vec	*b;
	Vec	*u;
	
	int	size, rank;
	
	MPI_Comm_size(PETSC_COMM_WORLD, &size);
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	
	A	= assem->A;
	b	= assem->b;
	u	= assem->u;
	
	KSP	ksp;
	PC	pc;
	
	PetscLogStage	stage, stageSolve;
	
	KSPCreate(PETSC_COMM_WORLD, &ksp);
	KSPSetType(ksp,KSPRICHARDSON);
	KSPSetOperators(ksp, *A, *A);
	KSPGetPC(ksp,&pc);
//	PCSetType(pc,PCASM);
	KSPMonitorSet(ksp, myMonitor, solver->rnorm, NULL);
	KSPSetTolerances(ksp, 1.e-7, PETSC_DEFAULT, PETSC_DEFAULT, solver->numIter);
	KSPSetFromOptions(ksp);

	double initWallTime = MPI_Wtime();
	clock_t solverInitT = clock();
	PetscLogStageRegister("Solver", &stageSolve);
	PetscLogStagePush(stageSolve);
	
	KSPSolve(ksp, *b, *u);
	
	PetscLogStagePop();
	clock_t solverT = clock();
	double endWallTime = MPI_Wtime();
	KSPGetIterationNumber(ksp, &(solver->numIter));

//	VecView(*u,PETSC_VIEWER_STDOUT_WORLD);
	PetscPrintf(PETSC_COMM_WORLD,"---------------------------| level = 0 |------------------------\n");
	KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
	PetscPrintf(PETSC_COMM_WORLD,"----------------------------------------------------------------\n");
	KSPDestroy(&ksp);

	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = [%d]; Solver cputime:                %lf\n",rank,(double)(solverT-solverInitT)/CLOCKS_PER_SEC);
	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = [%d]; Solver walltime:               %lf\n",rank,endWallTime-initWallTime);
	PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
}

void Solve(Assembly *assem, Solver *solver){
	// Solves the problem with chosen multigrid cycle
	
	if (solver->cycle == VCYCLE) MultigridVcycle(assem, solver);
	if (solver->cycle == ICYCLE) MultigridIcycle(assem, solver);
}


