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

void Range(Indices *indices, Solver *solver) {
	// Computes the range of global indices in this process for all levels	
	//
	// range[level][0] = Starting global index in level
	// range[level][1] = 1+(ending global index) in level
	
	int	remainder, quotient, totaln;
	int	procs, rank;
	int	(*range)[2];
	
	MPI_Comm_size(PETSC_COMM_WORLD, &procs);
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	
	range = solver->range;
	for (int l=0;l<indices->levels;l++) {
		totaln = indices->level[l].global.ni;
		remainder = (totaln)%procs;
		quotient  = (totaln)/procs;
		if (rank<remainder) {
			range[l][0] = rank*(quotient + 1);
			range[l][1] = range[l][0] + (quotient + 1);
		}
		else {
			range[l][0] = rank*quotient + remainder;
			range[l][1] = range[l][0] + quotient;
		}
	}
}

void SetUpSolver(Indices *indices, Solver *solver, Cycle cyc) {
	// Allocates memory to Solver struct
		
	solver->rnorm = malloc((solver->numIter+1)*sizeof(double));
	solver->range = malloc(indices->levels*sizeof(int[2]));
	solver->cycle = cyc;
	Range(indices, solver);
}

void DestroySolver(Solver *solver) {
	// Free the memory in Solver struct
	
	free(solver->rnorm);
	free(solver->range);
}

void SetPostProcess(PostProcess *pp) {
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

//void GetError(Problem *prob, Mesh *mesh, Indices *indices, Assembly *assem, double *error) {
//	
//	// u(x,y) = sin(Pi*x)*sin(pi*y)	
//	double		diff, sol;
//	double		*u;
//	int		range[2];
//	ArrayInt2d	*grid;
//	int		gridId;
//	int		globalni, globalnj, *global;
//	int		ifine, jfine;
//	int		i, j, g;
//	double		**coord;
//	double		localError[3];
//
//	int		rank;
//	
//	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
//	
//	VecGetOwnershipRange(assem->level[0].u, range, range+1);
//	VecGetArray(assem->level[0].u, &u);
//	
//	globalni = indices->level[0].global.ni;
//	globalnj = indices->level[0].global.nj;
//	global   = indices->level[0].global.data;
//	gridId   = indices->level[0].gridId[0];
//
//	coord = mesh->coord;
//
//	localError[0] = 0.0;
//	localError[1] = 0.0;
//	localError[2] = 0.0;
//	for (int row=range[0];row<range[1];row++) {
//		i = global[row*globalnj    ];
//		j = global[row*globalnj + 1];
//		g = global[row*globalnj + 2];
//		if (g != gridId) continue;
//		sol = prob->SOLfunc(coord[0][j+1], coord[1][i+1]);
//		diff = fabs(u[row-range[0]]-sol);
//		localError[0] = fmax(diff,localError[0]);
//		localError[1] += diff;
//		localError[2] += diff*diff;
//	}
//	MPI_Reduce(localError, error, 1, MPI_DOUBLE, MPI_MAX, 0, PETSC_COMM_WORLD);
//	MPI_Reduce(localError+1, error+1, 2, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
//	if (rank == 0) error[2] = sqrt(error[2]);
//	
//	VecRestoreArray(assem->level[0].u, &u);
//}

static void GetSol1(Indices *indices, Assembly *assem, Array2d u) {
	
	int		r;
	double		*px;
	const	int	*ranges;
	int		gridId;
	int		globalni, globalnj, *global;
	int		i, j, g;
	int		count;
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
		count = 0;
		for (int row=ranges[rank];row<ranges[rank+1];row++) {
			g = global[row*globalnj + 2];
			if (g != gridId) continue;
			buffer[count] = px[row-ranges[rank]];
			count += 1;
		}

		MPI_Send(buffer, ranges[rank+1]-ranges[rank], MPI_DOUBLE, 0, rank, PETSC_COMM_WORLD);
		free(buffer);
		printf("rank: %d; I am here\n", rank);
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
			MPI_Recv(buffer+ranges[i], ranges[i+1]-ranges[i], MPI_DOUBLE, i, i, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	
		for (int row=ranges[0];row<ranges[procs];row++) {
			i = global[row*globalnj    ];
			j = global[row*globalnj + 1];
			g = global[row*globalnj + 2];
			if (g != gridId) continue;
			u.data[i*u.nj+j] = buffer[row];
		}
		free(buffer);
		printf("rank: %d; I am here\n", rank);
	}
	VecRestoreArray(assem->u[0], &px);

}

void Postprocessing(Problem *prob, Mesh *mesh, Indices *indices, Assembly *assem, Solver *solver, PostProcess *pp) {
	// Computes error and writes data to files
	
	int		rank;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

	Array2d		u;
	
	if (rank==0) CreateArray2d(indices->level[0].grid[0].ni, indices->level[0].grid[0].nj, &u);
	GetSol1(indices, assem, u);
//	GetError(prob, mesh, indices, assem, pp->error);
	
	if (rank==0) {	
	
//	pp->solData = fopen("uData.dat","w");
//	pp->resData = fopen("rData.dat","w");
//	pp->errData = fopen("eData.dat","w");
	
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

void ProlongationOperator(Array2d pro) {
	// Builds prolongation 2D stencilwise operator (pro)
	// Stencil size: pro.ni x pro.nj
	if(pro.ni != 3 || pro.nj != 3) {printf("Error in ProlongationOperator"); return;}
	for (int i=0;i<pro.ni;i++) {
 		pro.data[i*pro.nj]= 0.5 - 0.25*fabs(1-i);
 		pro.data[i*pro.nj+1]= 1.0 - 0.5*fabs(1-i);
 		pro.data[i*pro.nj+2]= 0.5 - 0.25*fabs(1-i);
	}
}

void RestrictionOperator(Array2d res) {
	// Builds Restriction 2D stencilwise operator (res)
	// Stencil size: res.ni x res.nj
	
	if(res.ni != 3 || res.nj != 3) {printf("Error in RestrictionOperator"); return;}
	for (int i=0;i<res.nj;i++) {
 		res.data[i*res.nj]= 0.0;
 		res.data[i*res.nj+1]= 0.0;
 		res.data[i*res.nj+2]= 0.0;
	}
	res.data[res.nj+1] = 1.0;
}

void GridTransferOperator(Array2d *Iop, int factor, int totalGrids) {
	// Builds stencilwise grid transfer operator between any two grids that have "x" grids inbetween
	// where x = {0,...,totalGrids-2}
	
	int 	ni0, nj0;
	double	*weight;
	int 	nil, njl;
	double	*datal;
	int 	niu, nju;
	double	*datau;
	int 	iu, ju;
	
	ni0 = Iop[0].ni;
	nj0 = Iop[0].nj;
	weight = Iop[0].data;
	for (int l=0;l<totalGrids-2;l++) {
		nil = Iop[l].ni;
		njl = Iop[l].nj;
		datal = Iop[l].data;
		
		niu = Iop[l+1].ni;
		nju = Iop[l+1].nj;
		datau = Iop[l+1].data;
		// Initialization
		for (int i=0;i<niu*nju;i++) {
			datau[i] = 0.0;
		}
		// Building the operator
		for (int il=0;il<nil;il++) {
			for (int jl=0;jl<njl;jl++) {
				iu = factor*(il+1)-1-ni0/2;
				ju = factor*(jl+1)-1-nj0/2;
				for (int i0=0;i0<ni0;i0++) {
					for (int j0=0;j0<nj0;j0++) {
						datau[(iu+i0)*nju+(ju+j0)] += weight[i0*nj0+j0]*datal[il*nil+jl];
					}
				}
			}
		}
	}

}

void GridTransferOperators(Operator op, Indices indices) {
	// Builds stencilwise grid transfer operators between any two grids that have "x" grids inbetween
	// where x = {0,...,levels-2}
	
	RestrictionOperator(op.res[0]);
	ProlongationOperator(op.pro[0]);
	
	GridTransferOperator(op.res, indices.coarseningFactor, op.totalGrids);
	GridTransferOperator(op.pro, indices.coarseningFactor, op.totalGrids);
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

//void MultigridVcycle(Array2d u, Array2d metrics, double *f, double **opIH2h, double **opIh2H, ArrayInt2d *IsStencil, ArrayInt2d *IsResStencil, ArrayInt2d *IsProStencil, IsRange *range, double *rnorm, int levels, int *fulln, int *m) {
//
//	int	v[2], n[levels];
//	Mat	A[levels], prolongMatrix[levels-1], restrictMatrix[levels-1];
//	int	iter;
//	double	rnormchk, bnorm;
//	
//	double	*px;
//	const	int	*ranges;
//
//	int	size, rank;
//	
//	KSP	solver[levels];
//	PC	pc[levels];
//	Vec	r[levels], rv[levels];//, xbuf[levels];
////	Vec	dummy1, dummy2;
//
//	PetscLogStage	stage;
//	
//	MPI_Comm_size(PETSC_COMM_WORLD, &size);
//	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
//	//printf("Enter the number of fine grid sweeps = ");
//	
//	scanf("%d",v);
//	//printf("Enter the number of coarse grid sweeps = ");
//	scanf("%d",v+1);
//
////	printf("rank = %d; v[0] = %d\n",rank,v[0]);
////	printf("rank = %d; v[1] = %d\n",rank,v[1]);
//
//	n[0] = fulln[0]-2;
//	for (int i=1;i<levels;i++) {
//		n[i] = (n[i-1]-1)/2;
//	}
//
//	for (int i=0;i<levels;i++) {
//		VecDuplicate(b[i],&(rv[i]));
//	}
//	VecNorm(level[0].b, NORM_2, &bnorm);
////	printf("rank = %d, bnorm = %f\n", rank, bnorm);
////	VecView(r[0], PETSC_VIEWER_STDOUT_WORLD);
//	
//	KSPCreate(PETSC_COMM_WORLD, &(solver[0]));
////	KSPSetType(solver[0],KSPGMRES);
//	KSPSetType(solver[0],KSPRICHARDSON);
//	KSPRichardsonSetScale(solver[0],2.0/3.0);
//	KSPSetOperators(solver[0], A[0], A[0]);
//	KSPGetPC(solver[0],&(pc[0]));
////	PCSetType(pc[0],PCASM);
////	PCASMSetType(pc[0],PC_ASM_BASIC);
////	PCASMSetOverlap(pc[0],3);
////	PCASMSetTotalSubdomains(pc[0], 32, NULL, NULL);
//	PCSetType(pc[0],PCJACOBI);
//	KSPSetNormType(solver[0],KSP_NORM_NONE);
//	KSPSetTolerances(solver[0], 1.e-7, PETSC_DEFAULT, PETSC_DEFAULT, v[0]);
//	
//	for (int i=1;i<levels-1;i++) {
//		KSPCreate(PETSC_COMM_WORLD, &(solver[i]));
////		KSPSetType(solver[i],KSPGMRES);
//		KSPSetType(solver[i],KSPRICHARDSON);
//		KSPRichardsonSetScale(solver[i],2.0/3.0);
//		KSPSetOperators(solver[i], A[i], A[i]);
//		KSPGetPC(solver[i],&(pc[i]));
////		PCSetType(pc[i],PCASM);
////		PCASMSetType(pc[i],PC_ASM_BASIC);
////		PCASMSetOverlap(pc[i],3);
////		PCASMSetTotalSubdomains(pc[i], 32, NULL, NULL);
//		PCSetType(pc[i],PCJACOBI);
//		KSPSetNormType(solver[i],KSP_NORM_NONE);
//		KSPSetTolerances(solver[i], 1.e-7, PETSC_DEFAULT, PETSC_DEFAULT, v[0]);
//	}
//
//	if (levels>1) {
//		KSPCreate(PETSC_COMM_WORLD, &(solver[levels-1]));
////		KSPSetType(solver[levels-1],KSPGMRES);
//		KSPSetType(solver[levels-1],KSPRICHARDSON);
//		KSPRichardsonSetScale(solver[levels-1],2.0/3.0);
//		KSPSetOperators(solver[levels-1], A[levels-1], A[levels-1]);
//		KSPGetPC(solver[levels-1],&(pc[levels-1]));
////		PCSetType(pc[levels-1],PCASM);
////		PCASMSetType(pc[levels-1],PC_ASM_BASIC);
////		PCASMSetOverlap(pc[levels-1],3);
////		PCASMSetTotalSubdomains(pc[levels-1], 32, NULL, NULL);
//		PCSetType(pc[levels-1],PCJACOBI);
//		KSPSetNormType(solver[levels-1],KSP_NORM_NONE);
//		KSPSetTolerances(solver[levels-1], 1.e-7, PETSC_DEFAULT, PETSC_DEFAULT, v[1]);
//	}
//
//	double initWallTime = MPI_Wtime();
//	clock_t solverInitT = clock();
//	PetscLogStageRegister("Solver", &stage);
//	PetscLogStagePush(stage);
//	iter = 0;
//	rnormchk = bnorm;
//	if (rank==0) rnorm[0] = 1.0;
//	while (iter<*m && 100000000*bnorm > rnormchk && rnormchk > (1.e-7)*bnorm) {
//		KSPSolve(solver[0], b[0], x[0]);
//		if (iter==0) KSPSetInitialGuessNonzero(solver[0],PETSC_TRUE);
////		KSPBuildResidual(solver[0],NULL,NULL,&(r[0]));
////		VecView(r[0], PETSC_VIEWER_STDOUT_WORLD);
////		MatMult(restrictMatrix[0],r[0],b[1]);
////		VecView(b[1], PETSC_VIEWER_STDOUT_WORLD);
////		for (int l=1;l<levels-1;l++) {
//		for (int l=1;l<levels;l++) {
//			KSPBuildResidual(solver[l-1],NULL,rv[l-1],&(r[l-1]));
//			MatMult(restrictMatrix[l-1],r[l-1],b[l]);
//			KSPSolve(solver[l], b[l], x[l]);
//			if (l!=levels-1) KSPSetInitialGuessNonzero(solver[l],PETSC_TRUE);
////			KSPBuildResidual(solver[l],NULL,NULL,&(r[l]));
////			MatMult(restrictMatrix[l],r[l],b[l+1]);
//		}
////		KSPSolve(solver[levels-1], b[levels-1], x[levels-1]);
////		VecView(x[levels-1], PETSC_VIEWER_STDOUT_WORLD);
//		for (int l=levels-2;l>=0;l=l-1) {
////			MatMult(prolongMatrix[l],x[l+1],xbuf[l]);
////			VecAXPY(x[l],1.0,xbuf[l]);
//			MatMult(prolongMatrix[l],x[l+1],rv[l]);
//			VecAXPY(x[l],1.0,rv[l]);
//			KSPSolve(solver[l], b[l], x[l]);
////			KSPBuildResidual(solver[l],NULL,NULL,&(r[l]));
//			if (l!=0) KSPSetInitialGuessNonzero(solver[l],PETSC_FALSE);
//		}
////		MatMult(prolongMatrix[0],x[1],xbuf[0]);
////		VecView(xbuf[0], PETSC_VIEWER_STDOUT_WORLD);
////		VecView(x[0], PETSC_VIEWER_STDOUT_WORLD);
////		VecAXPY(x[0],1.0,xbuf[0]);
////		VecView(x[0], PETSC_VIEWER_STDOUT_WORLD);
////		KSPSolve(solver[0], b[0], x[0]);
//
//		KSPBuildResidual(solver[0],NULL,rv[0],&(r[0]));
//		VecNorm(r[0], NORM_2, &rnormchk);	
////		KSPGetResidualNorm(solver[0],&(rnormchk));
//		
////		printf("rank = %d; iter: %d, ul = %f, rnorm = %f, ll = %f\n", rank, iter, 100000000*bnorm, rnormchk, (1.e-7)*bnorm);
//		iter = iter + 1;
//		if (rank==0) rnorm[iter] = rnormchk/bnorm;
//	}
//	*m = iter;
//	PetscLogStagePop();
//	clock_t solverT = clock();
//	double endWallTime = MPI_Wtime();
//
//	for (int i=0;i<levels;i++) {
//		PetscPrintf(PETSC_COMM_WORLD,"---------------------------| level = %d |------------------------\n",i);
//		KSPView(solver[i],PETSC_VIEWER_STDOUT_WORLD);
//		PetscPrintf(PETSC_COMM_WORLD,"-----------------------------------------------------------------\n");
//	}
////	VecView(x[0], PETSC_VIEWER_STDOUT_WORLD);
//
////	rnorm[0] = Residual(u,f,r,As,n);
////	for (int i=1;i<m+1;i++) {
////		Vcycle(u,f,r,As,w,v,levels,n);
////		rnorm[i] = Residual(u,f,r,As,n); 
////	}
////	printf("residual = %.16e\n",rnorm[m]);
//	VecGetArray(x[0],&px);
//	VecGetOwnershipRanges(x[0],&ranges);
//	GetSol(u,px,fulln,levels,ranges,size,rank);
//	VecRestoreArray(x[0],&px);
//	
//	for (int i=0;i<levels;i++) {
//		MatDestroy(&(A[i]));
//		VecDestroy(&(rv[i]));
//		VecDestroy(&(x[i]));
//		VecDestroy(&(b[i]));
//	}
//	for (int l=0;l<levels-1;l++) {
//		MatDestroy(&(restrictMatrix[l]));
//		MatDestroy(&(prolongMatrix[l]));
//	}
//	
//	for (int i=0;i<levels;i++) {
//		KSPDestroy(&(solver[i]));
//	}
//	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = [%d]; Solver cputime:                %lf\n",rank,(double)(solverT-solverInitT)/CLOCKS_PER_SEC);
//	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = [%d]; Solver walltime:               %lf\n",rank,endWallTime-initWallTime);
//	PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
//
//}

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


