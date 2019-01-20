#include "solver.h"

#define ERROR_MSG(message) (fprintf(stderr,"ERROR: %s:%d: %s\n",__FILE__,__LINE__,(message)))
#define ERROR_RETURN(message) {ERROR_MSG(message);return ierr;}
#define CHKERR_PRNT(message) {if(ierr != 0) {ERROR_MSG(message);}}
#define CHKERR_RETURN(message) {if(ierr != 0) {ERROR_RETURN(message);}}
#define PI 3.14159265358979323846

#define pERROR_MSG(message) (PetscPrintf(PETSC_COMM_WORLD,"ERROR: %s:%d: %s\n",__FILE__,__LINE__,(message)))
#define pERROR_RETURN(message) {pERROR_MSG(message);return ierr;}
#define pCHKERR_PRNT(message) {if(ierr != 0) {pERROR_MSG(message);}}
#define pCHKERR_RETURN(message) {if(ierr != 0) {pERROR_RETURN(message);}}

#define METRICS(i,j,k) (metrics.data[metrics.nk*((i)*metrics.nj+(j))+(k)])
#define F(i,j) (f.data[((i)*f.nj+(j))])
#define U(i,j) (u.data[((i)*u.nj+(j))])

typedef struct {
	Vec	rInner;
	Vec	residualInner; // dummy for the sake of kspbuildresidual
	Vec	*rGrid;
	int	grids;
	int	innerCount; // Inner sweep counter 
	double	**rNormGrid;
} D1cntx;

//static int ipow(int base, int exp) {
//
//	int result = 1;
//	while (exp) {
//		if (exp & 1)
//			result *= base;
//		exp >>= 1;
//		base *= base;
//	}
//	return result;
//}

void SetUpAssembly(Levels *levels, Assembly *assem, int cycle) {
	// Allocate memory for Assembly struct
	
	assem->levels = levels->nlevels;
	assem->A = malloc((assem->levels)*sizeof(Mat));
//	if (cycle == ECYCLE) assem->A2 = malloc((assem->levels)*sizeof(Mat)); 
	assem->u = malloc((assem->levels)*sizeof(Vec));
	assem->b = malloc((assem->levels)*sizeof(Vec));
//	if (cycle != D1CYCLE && cycle != D2CYCLE && cycle != D1PSCYCLE) {
		assem->res = malloc((assem->levels-1)*sizeof(Mat));
		assem->pro = malloc((assem->levels-1)*sizeof(Mat));
//	} else if (assem->moreInfo == 0) {
//		assem->bottomIS	= malloc((assem->levels)*sizeof(IS));
//		assem->topIS	= malloc((assem->levels)*sizeof(IS));
////		assem->subFineIS= malloc((assem->levels)*sizeof(IS));
//		assem->res = malloc((assem->levels)*sizeof(Mat));
//		assem->pro = malloc((assem->levels)*sizeof(Mat));
//	} else {
//		assem->bottomIS	= malloc((assem->levels)*sizeof(IS));
//		assem->topIS	= malloc((assem->levels)*sizeof(IS));
////		assem->subFineIS= malloc((assem->levels)*sizeof(IS));
//		assem->gridIS	= malloc((assem->levels)*sizeof(IS*));
//		for (int i=0; i<assem->levels; i++) {
//			assem->gridIS[i] = malloc((levels->level[i].grids)*sizeof(IS));
//		}
//		assem->res = malloc((assem->levels)*sizeof(Mat));
//		assem->pro = malloc((assem->levels)*sizeof(Mat));
//	}
}

void DestroyAssembly(Assembly *assem, int cycle) 
{
	// Free the memory in Assembly struct
	
	for (int l=0;l<assem->levels;l++) {
		MatDestroy(assem->A+l);
//		if (cycle == ECYCLE) MatDestroy(assem->A2+l);
		VecDestroy(assem->b+l);
		VecDestroy(assem->u+l);
//		if (cycle == D1CYCLE || cycle == D2CYCLE || cycle == D1PSCYCLE) {
//			ISDestroy(assem->bottomIS+l);
//			ISDestroy(assem->topIS+l);
//		}
	}
//	if (cycle == D1CYCLE || cycle == D2CYCLE || cycle == D1PSCYCLE) {
//		MatDestroy(assem->res);
//		MatDestroy(assem->pro);
//	}
	for (int l=0;l<assem->levels-1;l++) {
		MatDestroy(assem->res+l);
		MatDestroy(assem->pro+l);
	}
	free(assem->res);
	free(assem->pro);
	free(assem->A);
//	if (cycle == ECYCLE) free(assem->A2); 
	free(assem->b);
	free(assem->u);
//	if (cycle != D1CYCLE && cycle != D2CYCLE && cycle != D1PSCYCLE) return;
//	if (assem->moreInfo == 0) {
//		free(assem->bottomIS);
//		free(assem->topIS);
////		free(assem->subFineIS);
//	} else {
//		free(assem->bottomIS);
//		free(assem->topIS);
//		for (int i=0; i<assem->levels; i++) {
//			free(assem->gridIS[i]);
//		}
//		free(assem->gridIS);
////		free(assem->subFineIS);
//	}
}

void InitializeSolver(Solver *solver) {
	solver->rnorm	= NULL;
	solver->levels	= NULL;
}

int CreateSolver(Grids *grids, Solver *solver) {
	// Allocates memory to Solver struct
	
	int	ierr = 0;
	
	if (!grids || !solver) {
		pERROR_MSG("NULL pointers encountered");
		return 1;
	}
	InitializeSolver(solver);

	PetscBool	set;	
	ierr = PetscOptionsGetInt(NULL, NULL, "-cycle", &(solver->cycle), &set);
	if (!set || ierr) {
		PetscBarrier(PETSC_NULL);
		pERROR_MSG("Type of MG cycle for solver not set");
		pERROR_MSG("Set '-cycle n' for n-th type cycle");
		return 1;
	} else if (solver->cycle > 4) {
		PetscBarrier(PETSC_NULL);
		pERROR_MSG("Selected MG cycle doesn't exist");
		return 1;
	}
	solver->moreInfo = 0;
//	PetscOptionsGetInt(NULL, NULL, "-moreNorm", &(solver->moreInfo), NULL);
	ierr = PetscOptionsGetInt(NULL, NULL, "-iter", &(solver->numIter), &set);
	if (!set || ierr) {
		PetscBarrier(PETSC_NULL);
		pERROR_MSG("Number of iterations not set");
		pERROR_MSG("Set '-iter n' for n iterations");
		return 1;
	}
	int	vmax = 2;
	ierr = PetscOptionsGetIntArray(NULL, NULL, "-v", solver->v, &vmax, &set);
	if (!set || vmax < 2 || ierr) {
		PetscBarrier(PETSC_NULL);
		pERROR_MSG("No. of smoothing steps and coarse solver iterations not set properly");
		pERROR_MSG("Set '-v v1,v2' for v1 smoothing steps and v2 coarse solve iterations");
		return 1;
	}
//	solver->cycle = cyc;
//	solver->moreInfo = moreInfo;
	solver->rnorm = malloc((solver->numIter+1)*sizeof(double));
//	solver->assem = malloc(sizeof(Assembly));
	solver->levels = malloc(sizeof(Levels));
//	solver->assem->moreInfo = solver->moreInfo;	
//	SetUpAssembly(levels, solver->assem, solver->cycle);
	ierr = CreateLevels(grids, solver->levels); pCHKERR_RETURN("Levels creation failed");
//	if (solver->moreInfo == 0) return 0;
//	if (cyc == D1PSCYCLE) {
//		int	v	= solver->v[0];
//		solver->grids	= levels->level->ngrids;
//		solver->rNormGrid = malloc(solver->grids*sizeof(double *));
//		for (int i=0; i<solver->grids; i++) {
//			solver->rNormGrid[i] = malloc(numIter*2*(v+1)*sizeof(double));
//		}
//		solver->rNormGlobal = malloc(numIter*2*(v+1)*sizeof(double));
//	} else if (cyc == D1CYCLE || cyc == D2CYCLE) {
//		int	v	= solver->v[0];
//		solver->grids	= levels->level->ngrids;
//		solver->rNormGrid = malloc(solver->grids*sizeof(double *));
//		for (int i=0; i<solver->grids; i++) {
//			solver->rNormGrid[i] = malloc(numIter*(v+1)*sizeof(double));
//		}
//		solver->rNormGlobal = malloc(numIter*(v+1)*sizeof(double));
//	}
	return 0;
}

void DestroySolver(Solver *solver) {
	// Free the memory in Solver struct
	
	if (!solver) return;	
//	DestroyAssembly(solver->assem, solver->cycle);
//	free(solver->assem);
	DestroyLevels(solver->levels);
	if (solver->levels) free(solver->levels);
	if (solver->rnorm) free(solver->rnorm);
//	if (solver->moreInfo == 0) return;
//	for (int i=0; i<solver->grids; i++) {
//		free(solver->rNormGrid[i]);
//	}
//	free(solver->rNormGrid);
//	free(solver->rNormGlobal);
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
		pp->XgridData = fopen("XgridData.dat","w");
		pp->YgridData = fopen("YgridData.dat","w");
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
		fclose(pp->XgridData);
		fclose(pp->YgridData);
	}
}

//void fillJacobians(Problem *prob, Grids *grids, Level *level, Mat *A) {
//	// Fills Mat A with the Jacobian or Discretized PDE coefficients of all grids contained in this level
//
//	int		*a, *b;
//	int		ai, aj, bi, bj;
//	int		lg;
//	
//	int		grids, *gridId;
//	int		*ranges;
//	double		As[5];
//
//	int		i0, j0, g0;
//	int		ifine, jfine;
//
//	double		metrics[5], **coord;
//	
//	int	procs, rank;
//	
//	MPI_Comm_size(PETSC_COMM_WORLD, &procs);
//	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
//	
//	coord = mesh->coord;
//
//	ai = level->global.ni;
//	aj = level->global.nj;
//	a  = level->global.data;
//
//	grids = level->ngrids;
//	gridId = level->gridId;
//
//	ranges = level->ranges;	
//	
//	// Row-based fill:
//	for (int row=ranges[rank];row<ranges[rank+1];row++) {
//		//i0 - row    - y coord
//		//j0 - column - x coord
//		//A[0]*u(i0-1,j0) + A[1]*u(i0,j0-1) + A[2]*u(i0,j0) + A[3]*u(i0,j0+1) + A[4]*u(i0+1,j0) = f(i0,j0)
//		i0 = a[row*aj];
//		j0 = a[row*aj+1];
//		g0 = a[row*aj+2]; 
//		for (lg=0;lg<grids;lg++) {if (g0 == gridId[lg]) break;} 
//		
//		bi = level->grid[lg].ni;
//		bj = level->grid[lg].nj;
//		b  = level->grid[lg].data;
//		// fine grid point corresponding to (i0, j0)
//		ifine = ipow(factor,g0)*(i0+1)-1;
//		jfine = ipow(factor,g0)*(j0+1)-1;
//		
//		// Compute metrics (analytically) at physical point corresponding to fine grid point
//		mesh->MetricCoefficients(mesh, coord[0][jfine+1], coord[1][ifine+1], metrics);
//		prob->OpA(As, metrics, level->h[lg]); // Get coefficients
//
//		// Fill the matrix
//		if (i0-1>=0) {
//			MatSetValue(*A, row, b[(i0-1)*bj+j0], As[0], ADD_VALUES);
//		}
//		if (j0-1>=0) {
//			MatSetValue(*A, row, b[i0*bj+j0-1], As[1], ADD_VALUES);
//		}
//		MatSetValue(*A, row, row, As[2], ADD_VALUES);
//		if (j0+1<bj) {
//			MatSetValue(*A, row, b[i0*bj+j0+1], As[3], ADD_VALUES);
//		}
//		if (i0+1<bi) {
//			MatSetValue(*A, row, b[(i0+1)*bj+j0], As[4], ADD_VALUES);
//		}
//	}
//}

//void fillJacobians(Problem *prob, Mesh *mesh, Level *level, int factor, Mat *A) {
//	// Fills Mat A with the Jacobian or Discretized PDE coefficients of all grids this level possesses
//
//	int		*a, *b;
//	int		ai, aj, bi, bj;
//	int		lg;
//	
//	int		grids, *gridId;
//	int		*ranges;
//	double		As[5];
//
//	int		i0, j0, g0;
//	int		ifine, jfine;
//
//	double		metrics[5], **coord;
//	
//	int	procs, rank;
//	
//	MPI_Comm_size(PETSC_COMM_WORLD, &procs);
//	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
//	
//	coord = mesh->coord;
//
//	ai = level->global.ni;
//	aj = level->global.nj;
//	a  = level->global.data;
//
//	grids = level->ngrids;
//	gridId = level->gridId;
//
//	ranges = level->ranges;	
//	
//	// Row-based fill:
//	for (int row=ranges[rank];row<ranges[rank+1];row++) {
//		//i0 - row    - y coord
//		//j0 - column - x coord
//		//A[0]*u(i0-1,j0) + A[1]*u(i0,j0-1) + A[2]*u(i0,j0) + A[3]*u(i0,j0+1) + A[4]*u(i0+1,j0) = f(i0,j0)
//		i0 = a[row*aj];
//		j0 = a[row*aj+1];
//		g0 = a[row*aj+2]; 
//		for (lg=0;lg<grids;lg++) {if (g0 == gridId[lg]) break;} 
//		
//		bi = level->grid[lg].ni;
//		bj = level->grid[lg].nj;
//		b  = level->grid[lg].data;
//		// fine grid point corresponding to (i0, j0)
//		ifine = ipow(factor,g0)*(i0+1)-1;
//		jfine = ipow(factor,g0)*(j0+1)-1;
//		
//		// Compute metrics (analytically) at physical point corresponding to fine grid point
//		mesh->MetricCoefficients(mesh, coord[0][jfine+1], coord[1][ifine+1], metrics);
//		prob->OpA(As, metrics, level->h[lg]); // Get coefficients
//
//		// Fill the matrix
//		if (i0-1>=0) {
//			MatSetValue(*A, row, b[(i0-1)*bj+j0], As[0], ADD_VALUES);
//		}
//		if (j0-1>=0) {
//			MatSetValue(*A, row, b[i0*bj+j0-1], As[1], ADD_VALUES);
//		}
//		MatSetValue(*A, row, row, As[2], ADD_VALUES);
//		if (j0+1<bj) {
//			MatSetValue(*A, row, b[i0*bj+j0+1], As[3], ADD_VALUES);
//		}
//		if (i0+1<bi) {
//			MatSetValue(*A, row, b[(i0+1)*bj+j0], As[4], ADD_VALUES);
//		}
//	}
//}
//
//void fillRestrictionPortion(Problem *prob, Mesh *mesh, Operator *op, Level *level, int factor, Mat *A) {
//	// Fills restriction portions (such as I_h^H A_h) of level Mat A
//	
//	int		*a, *b;
//	int		ai, aj, bi, bj;
//	double		*res, *pro;
//	double		weight;
//	int		resni, resnj;
//	
//	int		grids, *gridId;
//	int		*ranges;
//	double		As[5];
//
//	int		i0, j0, g0;
//	int		ifine, jfine;
//	int		i1, j1, g1;
//
//	double		metrics[5], **coord;
//	
//	int	procs, rank;
//	
//	MPI_Comm_size(PETSC_COMM_WORLD, &procs);
//	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
//	
//	coord = mesh->coord;
//
//	ai = level->global.ni;
//	aj = level->global.nj;
//	a  = level->global.data;
//
//	grids = level->ngrids;
//	gridId = level->gridId;
//
//	ranges = level->ranges;	
//
//	// Row-based fill:
//	for (int row=ranges[rank];row<ranges[rank+1];row++) {
//		//i0 - row    - y coord
//		//j0 - column - x coord
//		//A[0]*u(i0-1,j0) + A[1]*u(i0,j0-1) + A[2]*u(i0,j0) + A[3]*u(i0,j0+1) + A[4]*u(i0+1,j0) = f(i0,j0)
//		i0 = a[row*aj];
//		j0 = a[row*aj+1];
//		g0 = a[row*aj+2]; 
//		for (int lg=0;lg<grids;lg++) {
//			g1 = gridId[lg];
//			if (g1 < g0) {	
//				// grid-to-global index map
//				bi = level->grid[lg].ni;
//				bj = level->grid[lg].nj;
//				b  = level->grid[lg].data;
//
//				// Restriction operator from g1 to g0
//				resni = op->res[g0-g1-1].ni;
//				resnj = op->res[g0-g1-1].nj;
//				res = op->res[g0-g1-1].data;
//				
//				// (i1, j1) in grid g1 corresponding to (i0, j0) in grid g0
//				i1 = ipow(factor,(g0-g1))*(i0+1)-1 - (resni)/2;
//				j1 = ipow(factor,(g0-g1))*(j0+1)-1 - (resnj)/2;	
//				for (int i=i1;i<i1 + resni;i++) {
//					for (int j=j1;j<j1 + resnj;j++) {
//						// fine grid point corresponding to (i, j) in grid g1
//						ifine = ipow(factor,(g1))*(i+1)-1;
//						jfine = ipow(factor,(g1))*(j+1)-1;
//
//						// Compute metrics at physical point corresponding to fine grid point
//						mesh->MetricCoefficients(mesh, coord[0][jfine+1], coord[1][ifine+1], metrics);
//						prob->OpA(As, metrics, level->h[lg]);// Get the coefficients 
//						
//						// Fill the values corresponding to restriction from (i, j) to (i0, j0)
//						weight = res[(i-i1)*resnj+(j-j1)];
//						if (weight == 0.0) continue;
//						if (i-1>=0) {
//							MatSetValue(*A, row, b[(i-1)*bj+j], weight*As[0], ADD_VALUES);
//						}
//						if (j-1>=0) {
//							MatSetValue(*A, row, b[i*bj+j-1], weight*As[1], ADD_VALUES);
//						}
//						MatSetValue(*A, row, b[i*bj+j], weight*As[2], ADD_VALUES);
//						if (j+1<bj) {
//							MatSetValue(*A, row, b[i*bj+j+1], weight*As[3], ADD_VALUES);
//						}
//						if (i+1<bi) {
//							MatSetValue(*A, row, b[(i+1)*bj+j], weight*As[4], ADD_VALUES);
//						}
//					}
//				}
//			}
//		}
//	}
//}
//
//void fillProlongationPortion(Problem *prob, Mesh *mesh, Operator *op, Level *level, int factor, Mat *A) {
//	// Fills prolongation portions (such as A_h I_H^h) of level Mat A
//	
//	int		*a, *b;
//	int		ai, aj, bi, bj;
//	double		*res, *pro;
//	double		weight;
//	int		proni, pronj;
//	
//	int		grids, *gridId;
//	int		*ranges;
//	double		As[5];
//
//	int		i0, j0, g0;
//	int		ifine, jfine;
//	int		i1, j1, g1;
//
//	double		metrics[5], **coord;
//	
//	int	procs, rank;
//	
//	MPI_Comm_size(PETSC_COMM_WORLD, &procs);
//	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
//	
//	coord = mesh->coord;
//
//	// Global-to-grid index map
//	ai = level->global.ni;
//	aj = level->global.nj;
//	a  = level->global.data;
//
//	grids = level->ngrids;
//	gridId = level->gridId;
//
//	ranges = level->ranges;	
//		
//	// Column based fill
//	for (int col=ranges[rank];col<ranges[rank+1];col++) {
//		i0 = a[col*aj];
//		j0 = a[col*aj+1];
//		g0 = a[col*aj+2];
//
//		for (int lg=0;lg<grids;lg++) {
//			g1 = gridId[lg];
//			if (g1 < g0) {
//				// Grid-to-global index map
//				bi = level->grid[lg].ni;
//				bj = level->grid[lg].nj;
//				b  = level->grid[lg].data;
//				
//				// Prolongation operator from g0 to g1
//				proni = op->pro[g0-g1-1].ni;
//				pronj = op->pro[g0-g1-1].nj;
//				pro = op->pro[g0-g1-1].data;
//				
//				// (i1, j1) on g1 corresponding to (i0, j0) on g0
//				i1 = ipow(factor,(g0-g1))*(i0+1)-1 - (proni)/2;
//				j1 = ipow(factor,(g0-g1))*(j0+1)-1 - (pronj)/2;
//
//				// Fill the values associated to interior points in 2d stencil 
//				for (int i=1;i<proni-1;i++) {
//					for (int j=1;j<pronj-1;j++) {
//						// fine grid point corresponding to (i+i1, j+j1) on g1
//						ifine = ipow(factor,(g1))*(i+i1+1)-1;
//						jfine = ipow(factor,(g1))*(j+j1+1)-1;
//
//						// Compute metrics at physical point correspoding to fine grid point
//						mesh->MetricCoefficients(mesh, coord[0][jfine+1], coord[1][ifine+1], metrics);
//						prob->OpA(As, metrics, level->h[lg]); // Get coefficients
//
//						// Fill the values
//						weight = As[0]*pro[(i-1)*pronj+(j)] + As[1]*pro[(i)*pronj+(j-1)]+ As[2]*pro[(i)*pronj+(j)]+ As[3]*pro[(i)*pronj+(j+1)]+ As[4]*pro[(i+1)*pronj+(j)];
//						if (weight != 0.0) MatSetValue(*A, b[(i+i1)*bj+j+j1], col, weight, ADD_VALUES);
//					}
//				}
//
//				// Fill the values associated to edges in 2d stencil except corners
//				for (int j=1;j<pronj-1;j++) {
//					ifine = ipow(factor,(g1))*(i1+1)-1;
//					jfine = ipow(factor,(g1))*(j+j1+1)-1;
//					mesh->MetricCoefficients(mesh, coord[0][jfine+1], coord[1][ifine+1], metrics);
//					prob->OpA(As, metrics, level->h[lg]);
//					weight = As[1]*pro[(j-1)]+ As[2]*pro[(j)]+ As[3]*pro[(j+1)]+ As[4]*pro[pronj+(j)];
//					if (weight != 0.0) MatSetValue(*A, b[(i1)*bj+j+j1], col, weight, ADD_VALUES);
//					
//					ifine = ipow(factor,(g1))*(proni+i1)-1;
//					jfine = ipow(factor,(g1))*(j+j1+1)-1;
//					mesh->MetricCoefficients(mesh, coord[0][jfine+1], coord[1][ifine+1], metrics);
//					prob->OpA(As, metrics, level->h[lg]);
//					weight = As[0]*pro[(proni-2)*pronj+(j)] + As[1]*pro[(proni-1)*pronj+(j-1)]+ As[2]*pro[(proni-1)*pronj+(j)]+ As[3]*pro[(proni-1)*pronj+(j+1)];
//					if (weight != 0.0) MatSetValue(*A, b[(proni-1+i1)*bj+j+j1], col, weight, ADD_VALUES);
//				}
//				
//				for (int i=1;i<proni-1;i++) {
//					ifine = ipow(factor,(g1))*(i+i1+1)-1;
//					jfine = ipow(factor,(g1))*(j1+1)-1;
//					mesh->MetricCoefficients(mesh, coord[0][jfine+1], coord[1][ifine+1], metrics);
//					prob->OpA(As, metrics, level->h[lg]);
//					weight = As[0]*pro[(i-1)*pronj] + As[2]*pro[(i)*pronj]+ As[3]*pro[(i)*pronj+(1)]+ As[4]*pro[(i+1)*pronj];
//					if (weight != 0.0) MatSetValue(*A, b[(i+i1)*bj+j1], col, weight, ADD_VALUES);
//						
//					ifine = ipow(factor,(g1))*(i+i1+1)-1;
//					jfine = ipow(factor,(g1))*(pronj+j1)-1;
//					mesh->MetricCoefficients(mesh, coord[0][jfine+1], coord[1][ifine+1], metrics);
//					prob->OpA(As, metrics, level->h[lg]);
//					weight = As[0]*pro[(i-1)*pronj+(pronj-1)] + As[1]*pro[(i)*pronj+(pronj-2)]+ As[2]*pro[(i)*pronj+(pronj-1)]+ As[4]*pro[(i+1)*pronj+(pronj-1)];
//					if (weight != 0.0) MatSetValue(*A, b[(i+i1)*bj+pronj-1+j1], col, weight, ADD_VALUES);
//				}
//				
//				// Fill the values associated to corners in 2d stencil
//				ifine = ipow(factor,(g1))*(i1+1)-1;
//				jfine = ipow(factor,(g1))*(j1+1)-1;
//				mesh->MetricCoefficients(mesh, coord[0][jfine+1], coord[1][ifine+1], metrics);
//				prob->OpA(As, metrics, level->h[lg]);
//				weight = As[2]*pro[0]+ As[3]*pro[1]+ As[4]*pro[pronj];
//				if (weight != 0.0) MatSetValue(*A, b[(i1)*bj+j1], col, weight, ADD_VALUES);
//				
//				ifine = ipow(factor,(g1))*(i1+1)-1;
//				jfine = ipow(factor,(g1))*(pronj+j1)-1;
//				mesh->MetricCoefficients(mesh, coord[0][jfine+1], coord[1][ifine+1], metrics);
//				prob->OpA(As, metrics, level->h[lg]);
//				weight = As[1]*pro[(pronj-2)]+ As[2]*pro[(pronj-1)]+ As[4]*pro[pronj+(pronj-1)];
//				if (weight != 0.0) MatSetValue(*A, b[(i1)*bj+pronj-1+j1], col, weight, ADD_VALUES);
//				
//				ifine = ipow(factor,(g1))*(proni+i1)-1;
//				jfine = ipow(factor,(g1))*(j1+1)-1;
//				mesh->MetricCoefficients(mesh, coord[0][jfine+1], coord[1][ifine+1], metrics);
//				prob->OpA(As, metrics, level->h[lg]);
//				weight = As[0]*pro[(proni-2)*pronj] + As[2]*pro[(proni-1)*pronj]+ As[3]*pro[(proni-1)*pronj+(1)];
//				if (weight != 0.0) MatSetValue(*A, b[(proni-1+i1)*bj+j1], col, weight, ADD_VALUES);
//				
//				ifine = ipow(factor,(g1))*(proni+i1)-1;
//				jfine = ipow(factor,(g1))*(pronj+j1)-1;
//				mesh->MetricCoefficients(mesh, coord[0][jfine+1], coord[1][ifine+1], metrics);
//				prob->OpA(As, metrics, level->h[lg]);
//				weight = As[0]*pro[(proni-2)*pronj+(pronj-1)] + As[1]*pro[(proni-1)*pronj+(pronj-2)]+ As[2]*pro[(proni-1)*pronj+(pronj-1)];
//				if (weight != 0.0) MatSetValue(*A, b[(proni-1+i1)*bj+pronj-1+j1], col, weight, ADD_VALUES);
//			}
//		}
//	}
//}
//
//void levelMatrixA(Problem *prob, Mesh *mesh, Operator *op, Level *level, int factor, Mat *A) {
//	// Build matrix "A" for a given level
//	// level - contains global-to-grid, grid-to-global index maps
//	// factor - coarsening factor
//	
//	int	rank;
//	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
//	
//	int		grids;
//	int		*ranges;
//	grids = level->ngrids;
//	ranges = level->ranges;	
//	
//	MatCreateAIJ(PETSC_COMM_WORLD, ranges[rank+1]-ranges[rank], ranges[rank+1]-ranges[rank], PETSC_DETERMINE, PETSC_DETERMINE, 6*grids, PETSC_NULL, 6*grids, PETSC_NULL, A);
//
//	fillJacobians(prob, mesh, level, factor, A);
//	fillRestrictionPortion(prob, mesh, op, level, factor, A);
//	fillProlongationPortion(prob, mesh, op, level, factor, A);
//	
//	MatAssemblyBegin(*A,MAT_FINAL_ASSEMBLY);
//	MatAssemblyEnd(*A,MAT_FINAL_ASSEMBLY);
//}

void levelMatrixA1(Problem *prob, Grids *grids, Solver *solver) {
	// Build matrix "A" for a given level
	
	int	rank;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	
	int	nnz	= 2*(grids->topo->dimension)+1; // Number of non-zeros per row
	int	nlevels	= solver->levels->nlevels;
	Level	*level	= solver->levels->level;
	Mat	*A	= solver->levels->A;

	for (int l=0; l<nlevels; l++) {
		int	ngrids = level[l].ngrids;
		int	size = (int) (level[l].ranges[ngrids] - level[l].ranges[0]);
		MatCreateAIJ(PETSC_COMM_WORLD, size, size, PETSC_DETERMINE, PETSC_DETERMINE, nnz, PETSC_NULL, nnz, PETSC_NULL, A+l);

//		fillJacobians(prob, grids, level+l, A+l);
//		fillRestrictionPortion(prob, mesh, op, level, factor, A);
//		fillProlongationPortion(prob, mesh, op, level, factor, A);
		
		MatAssemblyBegin(A[l], MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(A[l], MAT_FINAL_ASSEMBLY);
	}
}
//
//void levelMatrixA1(Problem *prob, Mesh *mesh, Operator *op, Level *level, int factor, Mat *A) {
//	// Build matrix "A" for a given level
//	// level - contains global-to-grid, grid-to-global index maps
//	// factor - coarsening factor
//	
//	int	rank;
//	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
//	
//	int		grids;
//	int		*ranges;
//	grids = level->ngrids;
//	ranges = level->ranges;	
//	
//	MatCreateAIJ(PETSC_COMM_WORLD, ranges[rank+1]-ranges[rank], ranges[rank+1]-ranges[rank], PETSC_DETERMINE, PETSC_DETERMINE, 6, PETSC_NULL, 6, PETSC_NULL, A);
//
//	fillJacobians(prob, mesh, level, factor, A);
////	fillRestrictionPortion(prob, mesh, op, level, factor, A);
////	fillProlongationPortion(prob, mesh, op, level, factor, A);
//	
//	MatAssemblyBegin(*A,MAT_FINAL_ASSEMBLY);
//	MatAssemblyEnd(*A,MAT_FINAL_ASSEMBLY);
//}
//
//void levelMatrixA2(Problem *prob, Mesh *mesh, Operator *op, Level *level, int factor, Mat *A) {
//	// Build matrix "A" for a given level
//	// level - contains global-to-grid, grid-to-global index maps
//	// factor - coarsening factor
//	
//	int	rank;
//	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
//	
//	int		grids;
//	int		*ranges;
//	grids = level->ngrids;
//	ranges = level->ranges;	
//	
//	MatCreateAIJ(PETSC_COMM_WORLD, ranges[rank+1]-ranges[rank], ranges[rank+1]-ranges[rank], PETSC_DETERMINE, PETSC_DETERMINE, 6*(grids-1), PETSC_NULL, 6*(grids-1), PETSC_NULL, A);
//
////	fillJacobians(prob, mesh, level, factor, A);
//	fillRestrictionPortion(prob, mesh, op, level, factor, A);
//	fillProlongationPortion(prob, mesh, op, level, factor, A);
//	
//	MatAssemblyBegin(*A,MAT_FINAL_ASSEMBLY);
//	MatAssemblyEnd(*A,MAT_FINAL_ASSEMBLY);
//}
//
//void levelvecb(Problem *prob, Mesh *mesh, Operator *op, Level *level, int factor, Vec *b) {
//	// Build vector "b" for a given level
//	// f - logically 2D array containing right hand side values at each grid point
//	
//	int		*a;
//	int		ai, aj;
//	double		*res;
//	int		resni, resnj;
//	int		grids, *gridId;
//	
//	int		range[2];
//	double  	**coord;
//
//	int		i0, j0, g0;
//	int		ifine, jfine;
//	int		i1, j1, g1;
//
//	double		value;
//	
//	coord = mesh->coord;
//
//	ai = level->global.ni;
//	aj = level->global.nj;
//	a  = level->global.data;
//
//	grids = level->ngrids;
//	gridId = level->gridId;
//		
//	g1 = gridId[0];
//	
//	VecGetOwnershipRange(*b, range, range+1);
//	for (int row=range[0];row<range[1];row++) {
//		i0 = a[row*aj];
//		j0 = a[row*aj+1];
//		g0 = a[row*aj+2];
//		if (g0==g1) {
//			ifine = ipow(factor,g1)*(i0+1)-1; 
//			jfine = ipow(factor,g1)*(j0+1)-1;
//			value = prob->Ffunc(coord[0][jfine+1], coord[1][ifine+1]);
//			VecSetValue(*b, row, value, INSERT_VALUES);
//		} else {
//			resni = op->res[g0-g1-1].ni;
//			resnj = op->res[g0-g1-1].nj;
//			res = op->res[g0-g1-1].data;
//			
//			i1 = ipow(factor,(g0-g1))*(i0+1)-1 - (resni)/2;
//			j1 = ipow(factor,(g0-g1))*(j0+1)-1 - (resnj)/2;	
//			value = 0.0;
//			for (int i=i1;i<i1+resni;i++) {
//				for (int j=j1;j<j1+resnj;j++) {
//					ifine = ipow(factor,g1)*(i+1)-1; 
//					jfine = ipow(factor,g1)*(j+1)-1;
//					value += (prob->Ffunc(coord[0][jfine+1], coord[1][ifine+1]))*res[(i-i1)*resnj+(j-j1)];
//				}
//			}
//			VecSetValue(*b, row, value, INSERT_VALUES);
//		}
//	}
//	VecAssemblyBegin(*b);
//	VecAssemblyEnd(*b);
//
//	//return b;
//}
//
//void CreateSubLevel(Level *level, Level *sublevel, int flag) {
///********************************************************************************
// *
// * Allocate memory to the sub-level
// * 
// * Inputs:
// * 	level - source level to create sub-level
// * 	flag  - 0  : sub-level consists of all grids except last one
// * 		1  : sub-level consists of only last one
// * 		2  : sub-level consists of all grids except first one
// * 		3  : sub-level consists of only frist one
// * Output:
// * 	sublevel - sub-level of given level
// *
// ********************************************************************************/ 	
//	int	procs;
//	MPI_Comm_size(PETSC_COMM_WORLD, &procs);
//	
//	const	int	M = 3; // = cardinality{i,j,g}
//		int	totaln;
//	
//	int	*subGridId, *gridId;
//	int	(*subh)[2], (*h)[2];
//	int	subGlobalni, subGlobalnj;
//	int	subGridni, subGridnj;
//	int	grids;
//
//	gridId	= level->gridId;
//	h	= level->h;
//	grids	= level->ngrids;
//
//	if (flag == 1 || flag == 3) {
//		sublevel->ngrids	 = 1;
//	} else {
//		sublevel->ngrids	 = grids-1;
//	}
//
//	sublevel->gridId = malloc(sublevel->ngrids*sizeof(int));
//	sublevel->h	 = malloc(sublevel->ngrids*sizeof(double[2]));
//	sublevel->ranges = malloc((procs+1)*sizeof(int));
//	sublevel->grid	 = malloc(sublevel->ngrids*sizeof(ArrayInt2d));
//	
//	subGridId = sublevel->gridId;
//	subh	  = sublevel->h;
//
//	if (flag == 0) {
//		for (int lg=0; lg<sublevel->ngrids; lg++) {
//			subGridId[lg] = gridId[lg];
//			subh[lg][0] = h[lg][0];
//			subh[lg][1] = h[lg][1];
//			CreateArrayInt2d(level->grid[lg].ni, level->grid[lg].nj,&(sublevel->grid[lg]));
//		}
//		totaln = level->global.ni - (level->grid[grids-1].ni)*(level->grid[grids-1].nj);
//	} else if (flag == 1) {
//		subGridId[0] = gridId[grids-1];
//		subh[0][0] = h[grids-1][0];
//		subh[0][1] = h[grids-1][1];
//		CreateArrayInt2d(level->grid[grids-1].ni, level->grid[grids-1].nj,&(sublevel->grid[0]));
//		totaln = (level->grid[grids-1].ni)*(level->grid[grids-1].nj);
//	} else if (flag == 2) {
//		for (int lg=0; lg<sublevel->ngrids; lg++) {
//			subGridId[lg] = gridId[lg+1];
//			subh[lg][0] = h[lg+1][0];
//			subh[lg][1] = h[lg+1][1];
//			CreateArrayInt2d(level->grid[lg+1].ni, level->grid[lg+1].nj,&(sublevel->grid[lg]));
//		}
//		totaln = level->global.ni - (level->grid[0].ni)*(level->grid[0].nj);
//	} else if (flag == 3) {
//		subGridId[0] = gridId[0];
//		subh[0][0] = h[0][0];
//		subh[0][1] = h[0][1];
//		CreateArrayInt2d(level->grid[0].ni, level->grid[0].nj,&(sublevel->grid[0]));
//		totaln = (level->grid[0].ni)*(level->grid[0].nj);
//	}
//	CreateArrayInt2d(totaln, M, &(sublevel->global));
//}
//
//void DestroySubLevel(Level *sublevel) {
//	// Free the memory of sub-level
//	
//	for (int g=0;g<sublevel->ngrids;g++) {
//		DeleteArrayInt2d(&(sublevel->grid[g]));
//	}
//	DeleteArrayInt2d(&(sublevel->global));
//	free(sublevel->grid);
//	free(sublevel->h);
//	free(sublevel->gridId);
//	free(sublevel->ranges);
//}
//
//void subIS_based_on_grids(Level *level, int length, int *idg, IS *indexSet) {
///*********************************************************************************
// *
// * Creates a Petsc index set to hold a subset of global levels 
// * that belong to the grids given by an array
// *
// * Input:
// * 	level  - level info (mainly required for global index map and ranges)
// * 	length - length of "idg" array
// * 	idg    - Array containing the grid Ids
// * Output:
// * 	indexSet - Index set of global levels belonging to grids given by "idg"
// * 
// *********************************************************************************/
//	if (length > level->ngrids) {
//		ERROR_MSG("Error in extraction of sub global index sets");
//		return;
//	}
//	
//	int	rank;
//	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
//	
//	int	*ranges;
//	int	*global, nj;
//
//	ranges	= level->ranges;
//	global	= level->global.data;
//	nj	= level->global.nj;
//
//	int	*idx;   // Temporary sub levels holding array
//	idx	= calloc(ranges[rank+1]-ranges[rank], sizeof(int));  // Allocate maximum possible
//
//	int	g0, subi = 0;
//	for (int i=ranges[rank];i<ranges[rank+1];i++) {
//		g0 = global[i*nj+2]; // Get the grid Id of i^th global index
//		for (int g=0;g<length;g++) {
//			if (g0 == idg[g]) {
//				idx[subi] = i;
//				subi += 1;
//				break;
//			}
//		}
//	}
//	ISCreateGeneral(PETSC_COMM_WORLD, subi, idx, PETSC_COPY_VALUES, indexSet);
//	free(idx);
//}
//
//void getSubIS(Level *level, Level *sublevel, IS *indexSet) {
///*********************************************************************************
// *
// * Creates a Petsc index set that serves as a map from 'sub-level' global levels 
// * to 'level' global levels.
// *
// * Input:
// * 	level 	  - level info (mainly required for grid index map)
// * 	sublevel  - sub-level info (mainly required for its global index map and ranges)
// * Output:
// * 	indexSet - Petsc index set 
// * 
// *********************************************************************************/
//	int	rank;
//	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
//	
//	int	subRangeStart, subRangeEnd;
//	int	*subGlobal, subGlobalnj;
//
//	subGlobal	= sublevel->global.data;
//	subGlobalnj	= sublevel->global.nj;
//	subRangeStart	= sublevel->ranges[rank];
//	subRangeEnd	= sublevel->ranges[rank+1];
//
//	int		*gridId, grids;
//	ArrayInt2d	*grid;
//
//	grids	= level->ngrids;
//	gridId	= level->gridId;
//	grid	= level->grid;
//	
//	int	*idx;   // Temporary array holding sub-levels
//	idx	= calloc(subRangeEnd-subRangeStart, sizeof(int));
//
//	int	i, j;	
//	int	gId;
//	int	check = 0;
//	for (int row=subRangeStart; row<subRangeEnd; row++) {
//		i   = subGlobal[row*subGlobalnj  ];
//		j   = subGlobal[row*subGlobalnj+1];
//		gId = subGlobal[row*subGlobalnj+2];
//		for (int lg=0;lg<grids;lg++) {
//			if (gId == gridId[lg]) {
//				idx[row-subRangeStart] = grid[lg].data[i*grid[lg].nj+j];
//				check += 1;
//				break;
//			}
//		}
//	}
//	
//	ISCreateGeneral(PETSC_COMM_WORLD, subRangeEnd-subRangeStart, idx, PETSC_COPY_VALUES, indexSet);
//	free(idx);
//}
//
//void ComputeSubMaps(Level *level, Level *sublevel) {
///********************************************************************************
// *
// * Create global-to-grid and grid-to-global index maps for sub vectors, using
// * maps of actual vector.
// *
// * Input:
// * 	level	 - source level to use to create index maps in sub-level
// * 	sublevel - sub-level of the given "level"
// * Output:
// *	sublevel->global - Global to grid index map
// *	sublevel->grid[] - Grids to global index map
// *
// ********************************************************************************/
//	int	procs;
//	MPI_Comm_size(PETSC_COMM_WORLD, &procs);
//	
//	int	*global;
//	int	globalni, globalnj;
//	int	*ranges;
//
//	global		= level->global.data;
//	globalni	= level->global.ni;
//	globalnj	= level->global.nj;
//	ranges		= level->ranges;
//
//	int		*subGlobalData;
//	int		subGlobalni, subGlobalnj;
//	int		*subGridId;
//	int		subGrids;
//	int		*subRanges;
//	ArrayInt2d	*subGrid;
//	
//	subGlobalData	= sublevel->global.data;
//	subGlobalni	= sublevel->global.ni;
//	subGlobalnj	= sublevel->global.nj;
//	subGridId	= sublevel->gridId;	
//	subGrids	= sublevel->ngrids;
//	subGrid		= sublevel->grid;
//	subRanges	= sublevel->ranges;
//
//	int	i, j, g, subi = 0;
//	int	rank = 0;
//	subRanges[rank] = 0;
//	for (int ig=0; ig<globalni; ig++) {
//		i = global[ig*globalnj  ]; 
//		j = global[ig*globalnj+1]; 
//		g = global[ig*globalnj+2]; 
//		for (int slg=0; slg<subGrids; slg++) { // slg: sub_local_gridId
//			if (g == subGridId[slg]) {
//				subGlobalData[subi*subGlobalnj  ] = i;
//				subGlobalData[subi*subGlobalnj+1] = j;
//				subGlobalData[subi*subGlobalnj+2] = g;
//				
//				subGrid[slg].data[i*subGrid[slg].nj+j] = subi;
//				subi += 1;
//				break;
//			}
//		}
//		if (ig+1 == ranges[rank+1]) {
//			subRanges[rank+1] = subi;
//			rank += 1;
//		}
//	}
//}
//
//void Res_delayed(Levels *levels, Operator *op, Level *bottomlevel, Level *toplevel, Assembly *assem)
//{
///*********************************************************************************
// *
// * Builds restriction matrix for delayed cycling v1 within the level.!Only works
// * for one level for now.
// * 
// *********************************************************************************/ 	
//	// Case of more than 1 levels should handled outside, perhaps in poisson.c
//	int	procs, rank;
//	MPI_Comm_size(PETSC_COMM_WORLD, &procs);
//	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
//	
//	int	*topRanges, *botRanges;
//	topRanges = toplevel->ranges;
//	botRanges = bottomlevel->ranges;
//	
//	Mat	*res;
//	res = assem->res;
//	MatCreateAIJ(PETSC_COMM_WORLD, botRanges[rank+1]-botRanges[rank], topRanges[rank+1]-topRanges[rank], PETSC_DETERMINE, PETSC_DETERMINE, 9, PETSC_NULL, 9, PETSC_NULL, res);
//	
//	int	opResni, opResnj;
//	double	*opRes;
//	opResni = op->res[0].ni;
//	opResnj = op->res[0].nj;
//	opRes	= op->res[0].data;
//	
//	int	*botGlobal, botGlobalni, botGlobalnj;
//	botGlobal   = bottomlevel->global.data;
//	botGlobalni = bottomlevel->global.ni;
//	botGlobalnj = bottomlevel->global.nj;
//
//	ArrayInt2d	*topGrid;
//	topGrid = toplevel->grid;
//
//	int	topGrids, *topGridId;
//	topGrids = toplevel->ngrids;
//	topGridId = toplevel->gridId;
//	
//	int	factor;
//	factor = levels->coarseningFactor;
//	
//	int	*topGridData, topGridnj;
//	int	i1, j1, g1; // bottom grid stuff
//	int	i0, j0;// top grid stuff
//	double	weight;
//	for (int row=botRanges[rank]; row<botRanges[rank+1]; row++) {
//		i1 = botGlobal[row*botGlobalnj  ];
//		j1 = botGlobal[row*botGlobalnj+1];
//		g1 = botGlobal[row*botGlobalnj+2];
//		
//		// Identify the reference point (i0, j0) of restriction stencil on (g1-1)^th grid
//		i0 = factor*(i1+1)-1-(opResni)/2;
//		j0 = factor*(j1+1)-1-(opResnj)/2;
//		
//		// Searching for the local_gridId of "g1-1" in toplevel; Need a better way!
//		for (int lg=0; lg<topGrids; lg++) {
//			if (g1-1 == topGridId[lg]) {
//				topGridData = topGrid[lg].data;
//				topGridnj = topGrid[lg].nj;
//				break;
//			}
//		}
//		for (int i=i0; i<i0+opResni; i++) {
//			for (int j=j0; j<j0+opResnj; j++) {
//				weight = opRes[(i-i0)*opResnj+(j-j0)];
//				if (weight != 0.0) MatSetValue(*res, row, topGridData[i*topGridnj+j], weight, ADD_VALUES);
//			}
//		}
//	}
//	MatAssemblyBegin(*res, MAT_FINAL_ASSEMBLY);
//	MatAssemblyEnd(*res, MAT_FINAL_ASSEMBLY);
////	MatView(*res,PETSC_VIEWER_STDOUT_WORLD);
//
//}
//
//void Pro_delayed(Levels *levels, Operator *op, Level *bottomlevel, Level *toplevel, Assembly *assem) 
//{
///*********************************************************************************
// *
// * Builds prolongation matrix for delayed cycling v1 within the level.!Only works
// * for one level for now.
// * 
// *********************************************************************************/ 	
//	// Case of more than 1 levels should handled outside, perhaps in poisson.c
//	int	procs, rank;
//	MPI_Comm_size(PETSC_COMM_WORLD, &procs);
//	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
//	
//	int	*topRanges, *ranges;
//	topRanges = toplevel->ranges;
//	ranges = levels->level->ranges;
//	
//	Mat	*pro;
//	pro = assem->pro;
//	MatCreateAIJ(PETSC_COMM_WORLD, topRanges[rank+1]-topRanges[rank], ranges[rank+1]-ranges[rank], PETSC_DETERMINE, PETSC_DETERMINE, 4, PETSC_NULL, 4, PETSC_NULL, pro);
//	
//	int	opProni, opPronj;
//	double	*opPro;
//	opProni = op->pro[0].ni;
//	opPronj = op->pro[0].nj;
//	opPro	= op->pro[0].data;
//	
//	int	*global, globalni, globalnj;
//	global   = levels->level->global.data;
//	globalni = levels->level->global.ni;
//	globalnj = levels->level->global.nj;
//
//	ArrayInt2d	*topGrid;
//	topGrid = toplevel->grid;
//
//	int	*topGridId;
//	int	topGrids;
//	topGridId = toplevel->gridId;
//	topGrids  = toplevel->ngrids;
//	
//	
//	int	factor;
//	factor = levels->coarseningFactor;
//	
//	int	gfine;
//	gfine = levels->level->gridId[0];
//
//	int	*topGridData, topGridnj;
//	int	i1, j1, g1; // bottom level stuff
//	int	i0, j0;     // top level stuff
//	double	weight;
//	for (int col=ranges[rank]; col<ranges[rank+1]; col++) {
//		i1 = global[col*globalnj  ];
//		j1 = global[col*globalnj+1];
//		g1 = global[col*globalnj+2];
//		if (g1 == gfine) continue;	
//		// Identify the reference point (i0, j0) of restriction stencil on (g1-1)^th grid
//		i0 = factor*(i1+1)-1-(opProni)/2;
//		j0 = factor*(j1+1)-1-(opPronj)/2;
//		
//		// Searching for the local_gridId of "g1-1" in toplevel; Need a better way!
//		for (int lg=0; lg<topGrids; lg++) {
//			if (g1-1 == topGridId[lg]) {
//				topGridData = topGrid[lg].data;
//				topGridnj = topGrid[lg].nj;
//				break;
//			}
//		}
//		for (int i=i0; i<i0+opProni; i++) {
//			for (int j=j0; j<j0+opPronj; j++) {
//				weight = opPro[(i-i0)*opPronj+(j-j0)];
//				if (weight != 0.0) MatSetValue(*pro, topGridData[i*topGridnj+j], col, weight, ADD_VALUES);
//			}
//		}
//	}
//	MatAssemblyBegin(*pro, MAT_FINAL_ASSEMBLY);
//	MatAssemblyEnd(*pro, MAT_FINAL_ASSEMBLY);
////	MatView(*pro,PETSC_VIEWER_STDOUT_WORLD);
//}

//void Res(Levels *levels, Operator *op, int factor, Assembly *assem) {
//	// Assembles the restriction matrices
//	// Restriction is only from primary grid of one level to all grids of the next level
//	
//	int	levels;
//	
//	levels = assem->levels;
//	for (int l=0;l<levels-1;l++) {	
//		if (levels->level[l].grids>1) {
//			PetscPrintf(PETSC_COMM_WORLD, "For now, only 1 grid per level on all levels except last level is allowed in std multigrid\n");
//			return;
//		}
//	}
//
//	ArrayInt2d	grid0, grid1;
//	ArrayInt2d	global0, global1;
//	int		g0, g1;
//	int		i0, j0, i1, j1;
//	int		range0[2], range1[2];
//	Mat		*res;
//	int		opResni, opResnj;
//	double		*opRes;
//	double		weight;
//	
//	res = assem->res;
//	for (int l=0;l<levels-1;l++) {
//		g0 = levels->level[l].gridId[0];
//
//		global0 = levels->level[l].global;
//		global1 = levels->level[l+1].global;
//		
//		grid0 = levels->level[l].grid[0];
//		
//		VecGetOwnershipRange(assem->b[l], range0, range0+1);	
//		VecGetOwnershipRange(assem->b[l+1], range1, range1+1);	
//		
//		MatCreateAIJ(PETSC_COMM_WORLD, range1[1]-range1[0], range0[1]-range0[0], PETSC_DETERMINE, PETSC_DETERMINE, 9, PETSC_NULL, 9, PETSC_NULL, res+l);
//		for (int row=range1[0];row<range1[1];row++) {
//			i1 = global1.data[row*global1.nj];
//			j1 = global1.data[row*global1.nj+1];
//			g1 = global1.data[row*global1.nj+2];
//
//			opResni = op->res[g1-g0-1].ni;
//			opResnj = op->res[g1-g0-1].nj;
//			opRes = op->res[g1-g0-1].data;
//			
//			i0 = ipow(factor,(g1-g0))*(i1+1)-1-(opResni)/2;
//			j0 = ipow(factor,(g1-g0))*(j1+1)-1-(opResnj)/2;
//			for (int i=i0;i<i0+opResni;i++) {
//				for (int j=j0;j<j0+opResnj;j++) {
//					weight = opRes[(i-i0)*opResnj+(j-j0)];
//					if (weight != 0.0) MatSetValue(res[l], row, grid0.data[i*grid0.nj+j], weight, ADD_VALUES);
//				}
//			}
//		}
//	
//		MatAssemblyBegin(res[l], MAT_FINAL_ASSEMBLY);
//		MatAssemblyEnd(res[l], MAT_FINAL_ASSEMBLY);
//	}
//}
//
//void Pro(Levels *levels, Operator *op, int factor, Assembly *assem) {
//	// Assembles the prolongation matrix for level 1 to 0
//	
//	int	levels;
//	
//	levels = assem->levels;
//	for (int l=0;l<levels-1;l++) {	
//		if (levels->level[l].grids>1) {
//			PetscPrintf(PETSC_COMM_WORLD, "For now, only 1 grid per level on all levels except last level is allowed in std multigrid\n");
//			return;
//		}
//	}
//
//	ArrayInt2d	grid0, grid1;
//	ArrayInt2d	global0, global1;
//	int		g0, g1;
//	int		i0, j0, i1, j1;
//	int		range0[2], range1[2];
//	Mat		*pro;
//	int		opProni, opPronj;
//	double		*opPro;
//	double		weight;
//	
//	pro = assem->pro;
//	for (int l=0;l<levels-1;l++) {
//		g0 = levels->level[l].gridId[0];
//
//		global0 = levels->level[l].global;
//		global1 = levels->level[l+1].global;
//		
//		grid0 = levels->level[l].grid[0];
//		
//		VecGetOwnershipRange(assem->b[l], range0, range0+1);	
//		VecGetOwnershipRange(assem->b[l+1], range1, range1+1);	
//		
//		MatCreateAIJ(PETSC_COMM_WORLD, range0[1]-range0[0], range1[1]-range1[0], PETSC_DETERMINE, PETSC_DETERMINE, 4*(levels->level[l+1].grids), PETSC_NULL, 4*(levels->level[l+1].grids), PETSC_NULL, pro+l);
//		for (int col=range1[0];col<range1[1];col++) {
//			i1 = global1.data[col*global1.nj];
//			j1 = global1.data[col*global1.nj+1];
//			g1 = global1.data[col*global1.nj+2];
//
//			opProni = op->pro[g1-g0-1].ni;
//			opPronj = op->pro[g1-g0-1].nj;
//			opPro = op->pro[g1-g0-1].data;
//
//			i0 = ipow(factor,(g1-g0))*(i1+1)-1-(opProni)/2;
//			j0 = ipow(factor,(g1-g0))*(j1+1)-1-(opPronj)/2;
//			for (int i=i0;i<i0+opProni;i++) {
//				for (int j=j0;j<j0+opPronj;j++) {
//					weight = opPro[(i-i0)*opPronj+(j-j0)];
//					if (weight != 0.0) MatSetValue(pro[l], grid0.data[i*grid0.nj+j], col, weight, ADD_VALUES);
//				}
//			}
//		}
//	
//		MatAssemblyBegin(pro[l], MAT_FINAL_ASSEMBLY);
//		MatAssemblyEnd(pro[l], MAT_FINAL_ASSEMBLY);
//	}
//}
//
//void Assemble(Problem *prob, Grids *grids, Levels *levels, Operator *op, Solver *solver) {
//	// Assembles matrices, vectors and index sets in all levels
//	int		factor;
//	Assembly	*assem;
//	
//	assem = solver->assem;
//	factor = levels->coarseningFactor;
//	for (int l=0;l<assem->levels;l++) {
//		levelMatrixA(prob, mesh, op, &(levels->level[l]), factor, assem->A+l);
//		MatCreateVecs(assem->A[l], assem->u+l, assem->b+l);
//	}
//	// Only the zeroth level vec b is created
//	levelvecb(prob, mesh, op, levels->level, factor, assem->b);
//
//	if (assem->levels > 1) { 
//		Res(levels, op, factor, assem);
//		Pro(levels, op, factor, assem);
//	}
//}

//void Assemble(Problem *prob, Mesh *mesh, Levels *levels, Operator *op, Solver *solver) {
//	// Assembles matrices, vectors and index sets in all levels
//	int		factor;
//	Assembly	*assem;
//	
//	assem = solver->assem;
//	factor = levels->coarseningFactor;
//	for (int l=0;l<assem->levels;l++) {
//		if (solver->cycle == ECYCLE) {
//			levelMatrixA1(prob, mesh, op, &(levels->level[l]), factor, assem->A+l);
//			levelMatrixA2(prob, mesh, op, &(levels->level[l]), factor, assem->A2+l);
//		} else if (solver->cycle == D1CYCLE || solver->cycle == D2CYCLE || solver->cycle == D1PSCYCLE){
//			levelMatrixA1(prob, mesh, op, &(levels->level[l]), factor, assem->A+l);
//		} else {
//			levelMatrixA(prob, mesh, op, &(levels->level[l]), factor, assem->A+l);
//		}
//		MatCreateVecs(assem->A[l], assem->u+l, assem->b+l);
//	}
//	// Only the zeroth level vec b is created
//	levelvecb(prob, mesh, op, levels->level, factor, assem->b);
//
//	if (solver->cycle == D1CYCLE || solver->cycle == D2CYCLE || solver->cycle == D1PSCYCLE) {
//		Level	toplevel, bottomlevel;
////		Level	subFineLevel; // To be extracted from toplevel
//		
//		CreateSubLevel(levels->level, &toplevel, 0);
//		ComputeSubMaps(levels->level, &toplevel);
//		getSubIS(levels->level, &toplevel, assem->topIS);
//
////		CreateSubLevel(&toplevel, &subFineLevel, 3);
////		ComputeSubMaps(&toplevel, &subFineLevel);
////		getSubIS(&toplevel, &subFineLevel, assem->subFineIS);
//		
//		CreateSubLevel(levels->level, &bottomlevel, 2);
//		ComputeSubMaps(levels->level, &bottomlevel);
//		getSubIS(levels->level, &bottomlevel, assem->bottomIS);
//		
//		Res_delayed(levels, op, &bottomlevel, &toplevel, assem);
//		Pro_delayed(levels, op, &bottomlevel, &toplevel, assem);
//		
//		DestroySubLevel(&toplevel);
////		DestroySubLevel(&subFineLevel);
//		DestroySubLevel(&bottomlevel);
//		if (solver->moreInfo != 0) {
//			for (int i=0; i<solver->grids; i++) {
//				subIS_based_on_grids(levels->level, 1, levels->level->gridId+i, assem->gridIS[0]+i);
//			}
//		}
//
//	} else if (assem->levels > 1) { 
//		Res(levels, op, factor, assem);
//		Pro(levels, op, factor, assem);
//	}
//}

//void GetError(Problem *prob, Mesh *mesh, Array2d u1, double *error) {
//	
//	// u(x,y) = sin(Pi*x)*sin(pi*y)	
//	double	diff;
//	double	**coord;
//	double	sol;
//	int	uni, unj;
//	double	*u;
//
//	coord = mesh->coord;
//	uni   = u1.ni;
//	unj   = u1.nj;
//	u     = u1.data;
//	error[0] = 0.0;
//	error[1] = 0.0;
//	error[2] = 0.0;
//	for (int i=0;i<uni;i++) {
//		for (int j=0;j<unj;j++) {
//			sol = prob->SOLfunc(coord[0][j+1], coord[1][i+1]);
//			diff = fabs(u[i*unj+j]-sol);
//			error[0] = fmax(diff,error[0]);
//			error[1] = error[1] + diff;
//			error[2] = error[2] + diff*diff;
//		}
//	}
//	error[2] = sqrt(error[2]);
//}
//
//void GetSol(Levels *levels, Assembly *assem, Array2d u) {
//	
//	int		r;
//	double		*px;
//	const	int	*ranges;
//	int		gridId;
//	int		globalni, globalnj, *global;
//	int		i, j, g;
//	int		count, localcount;
//	double		*buffer;
//
//	int	procs, rank;
//	
//	MPI_Comm_size(PETSC_COMM_WORLD, &procs);
//	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
//	
//	VecGetArray(assem->u[0], &px);
//	VecGetOwnershipRanges(assem->u[0], &ranges);
//	
//	globalni = levels->level[0].global.ni;
//	globalnj = levels->level[0].global.nj;
//	global   = levels->level[0].global.data;
//	gridId   = levels->level[0].gridId[0];
//	
//	if (rank!=0) {
//		buffer = malloc((ranges[rank+1]-ranges[rank])*sizeof(double));
//		localcount = 0;
//		for (int row=ranges[rank];row<ranges[rank+1];row++) {
//			g = global[row*globalnj + 2];
//			if (g != gridId) continue;
//			buffer[localcount] = px[row-ranges[rank]];
//			localcount += 1;
//		}
//
//		MPI_Send(&localcount, 1, MPI_DOUBLE, 0, rank, PETSC_COMM_WORLD);
//		MPI_Send(buffer, localcount, MPI_DOUBLE, 0, rank*(localcount+1), PETSC_COMM_WORLD);
//		free(buffer);
//	}
//	else if (rank==0) {
//		int	totalN;
//		int	gridni, gridnj;
//		int	*grid;
//		double	*buffer;
//		
//		gridni = levels->level[0].grid[0].ni;
//		gridnj = levels->level[0].grid[0].nj;
//		grid   = levels->level[0].grid[0].data;
//		totalN = (gridni*gridnj);
//		buffer = malloc(totalN*sizeof(double));
//		
//		count = 0;
//		for (int row=ranges[0];row<ranges[1];row++) {
//			g = global[row*globalnj + 2];
//			if (g != gridId) continue;
//			buffer[count] = px[row-ranges[0]];
//			count += 1;
//		}
//
//		for (int i=1;i<procs;i++) {
//			MPI_Recv(&localcount, 1, MPI_DOUBLE, i, i, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
//			MPI_Recv(buffer+count, localcount, MPI_DOUBLE, i, i*(localcount+1), PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
//			count += localcount;
//		}
//		count = 0;
//		for (int row=ranges[0];row<ranges[procs];row++) {
//			i = global[row*globalnj    ];
//			j = global[row*globalnj + 1];
//			g = global[row*globalnj + 2];
//			if (g != gridId) continue;
//			u.data[i*u.nj+j] = buffer[count];
//			count += 1;
//		}
//		free(buffer);
//	}
//	VecRestoreArray(assem->u[0], &px);
//
//}
//
//void Postprocessing(Problem *prob, Mesh *mesh, Levels *levels, Solver *solver, PostProcess *pp) {
//	// Computes error and writes data to files
//	
//	int		rank;
//	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
//
//	Array2d		u;
//	
//	if (rank==0) CreateArray2d(levels->level[0].grid[0].ni, levels->level[0].grid[0].nj, &u);
//	GetSol(levels, solver->assem, u);
////	GetError(prob, mesh, levels, assem, pp->error);
//	
//	if (rank==0) {	
//	
//	GetError(prob, mesh, u, pp->error);
//	for (int i=0;i<3;i++) {
//		printf("\nerror[%d] = %.16e\n", i, pp->error[i]);
//		fprintf(pp->errData,"%.16e\n",pp->error[i]);
//	}
//	
//	double **coord;
//	coord = mesh->coord;
//	for (int i=0;i<u.ni;i++) {
//		for (int j=0;j<u.nj;j++) {
//			fprintf(pp->XgridData,"%lf    ", coord[0][j]);
//			fprintf(pp->YgridData,"%lf    ", coord[1][i]);
//			fprintf(pp->solData,"%.16e    ", u.data[i*u.nj+j]);
//		}
//		fprintf(pp->XgridData,"\n");
//		fprintf(pp->YgridData,"\n");
//		fprintf(pp->solData,"\n");
//	}
//		
//	for (int i=0;i<solver->numIter+1;i++) {
//		fprintf(pp->resData,"%.16e ",solver->rnorm[i]);
//	}
//	fprintf(pp->resData,"\n");
//	printf("\nRelative residual = %.16e\n", solver->rnorm[solver->numIter]);
//	
//	if (solver->moreInfo != 0) {	
//		FILE	*rGlobalData;
//		rGlobalData = fopen("rGlobal.dat","w");
//		for (int i=0;i<solver->numIter*(solver->v[0]+1);i++) {
//			fprintf(rGlobalData,"%.16e ",solver->rNormGlobal[i]);
//		}
//		fprintf(rGlobalData,"\n");
//		fclose(rGlobalData);
//
//		char	fileName[12];
//		for (int i=0; i<solver->grids; i++) {
//			FILE	*rGridData;
//			sprintf(fileName, "rGrid%d.dat",i);
//			rGridData = fopen(fileName, "w");
//			for (int j=0; j<solver->numIter*(solver->v[0]+1); j++) {
//				fprintf(rGridData, "%.16e ", (solver->rNormGrid)[i][j]);
//			}
//			fprintf(rGridData, "\n");
//			fclose(rGridData);
//		}
//	}
//
//	DeleteArray2d(&u);
//	}
//}
//
//PetscErrorCode  rNormGridMonitor(KSP ksp, PetscInt n, PetscReal rnormAtn, void *toCastInfo) {
//	//Writes the l2-norm of residual at each iteration to an array
//	
//	D1cntx	*info;
//	info = (D1cntx*)toCastInfo;
//	KSPBuildResidual(ksp, NULL, info->rInner, &(info->residualInner));
////	VecView(info->rInner, PETSC_VIEWER_STDOUT_WORLD);
////	VecView(info->rGrid[2], PETSC_VIEWER_STDOUT_WORLD);
//	double	temp;
//	for (int i=0; i<info->grids; i++) {
////		VecNorm(info->rGrid[i], NORM_2, (info->rNormGrid[i])+(info->innerCount));
//		VecDot(info->rGrid[i], info->rGrid[i], &temp);
//		info->rNormGrid[i][info->innerCount] = sqrt(temp);
////		PetscPrintf(PETSC_COMM_WORLD,"count: %d; rNormGrid[%d] = %lf\n", info->innerCount, i, temp);
//	}
//	(info->innerCount)++;
//	return 0;
//}
//
//PetscErrorCode  myMonitor(KSP ksp, PetscInt n, PetscReal rnormAtn, double *rnorm) {
//	//Writes the l2-norm of residual at each iteration to an array
//	//
//	//rnorm[n] = rnormInstant
//	
////	int	rank;
//
////	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
////	if (rank==0) rnorm[n] = rnormAtn;
//	rnorm[n] = rnormAtn;
//	return 0;
//}

int MultigridVcycle(Solver *solver) {

	int	size, rank;
	
	MPI_Comm_size(PETSC_COMM_WORLD, &size);
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	
	Mat 	*res;
	Mat 	*pro;
	Mat	*A;
	Vec	*b;
	Vec	*u;
	
	int	iter;
	double	rnormchk, bnorm;
	
	double	*rnorm	= solver->rnorm;
	int	maxIter	= solver->numIter;
	int	*v	= solver->v;
	int	levels	= solver->levels->nlevels;

//	KSP	ksp[levels];
////	PC	pc[levels];
//	Vec	r[levels], rv[levels];//, xbuf[levels];
//	
//	PetscLogStage	stage;
//	
//	for (int i=0;i<levels;i++) {
//		VecDuplicate(b[i],&(rv[i]));
//	}
//	
//	KSPCreate(PETSC_COMM_WORLD, &(ksp[0]));
////	KSPSetType(ksp[0],KSPGMRES);
//	KSPSetType(ksp[0],KSPRICHARDSON);
////	KSPRichardsonSetScale(ksp[0],2.0/3.0);
//	KSPSetOperators(ksp[0], A[0], A[0]);
////	KSPGetPC(ksp[0],&(pc[0]));
////	PCSetType(pc[0],PCASM);
////	PCASMSetType(pc[0],PC_ASM_BASIC);
////	PCASMSetOverlap(pc[0],3);
////	PCASMSetTotalSubdomains(pc[0], 32, NULL, NULL);
////	PCSetType(pc[0],PCJACOBI);
//	KSPSetNormType(ksp[0],KSP_NORM_NONE);
//	KSPSetTolerances(ksp[0], 1.e-7, PETSC_DEFAULT, PETSC_DEFAULT, v[0]);
//	KSPSetFromOptions(ksp[0]);
//	
//	for (int i=1;i<levels-1;i++) {
//		KSPCreate(PETSC_COMM_WORLD, &(ksp[i]));
////		KSPSetType(ksp[i],KSPGMRES);
//		KSPSetType(ksp[i],KSPRICHARDSON);
////		KSPRichardsonSetScale(ksp[i],2.0/3.0);
//		KSPSetOperators(ksp[i], A[i], A[i]);
////		KSPGetPC(ksp[i],&(pc[i]));
////		PCSetType(pc[i],PCASM);
////		PCASMSetType(pc[i],PC_ASM_BASIC);
////		PCASMSetOverlap(pc[i],3);
////		PCASMSetTotalSubdomains(pc[i], 32, NULL, NULL);
////		PCSetType(pc[i],PCJACOBI);
//		KSPSetNormType(ksp[i],KSP_NORM_NONE);
//		KSPSetTolerances(ksp[i], 1.e-7, PETSC_DEFAULT, PETSC_DEFAULT, v[0]);
//		KSPSetFromOptions(ksp[i]);
//	}
//
//	if (levels>1) {
//		KSPCreate(PETSC_COMM_WORLD, &(ksp[levels-1]));
////		KSPSetType(ksp[levels-1],KSPGMRES);
//		KSPSetType(ksp[levels-1],KSPRICHARDSON);
////		KSPRichardsonSetScale(ksp[levels-1],2.0/3.0);
//		KSPSetOperators(ksp[levels-1], A[levels-1], A[levels-1]);
////		KSPGetPC(ksp[levels-1],&(pc[levels-1]));
////		PCSetType(pc[levels-1],PCASM);
////		PCASMSetType(pc[levels-1],PC_ASM_BASIC);
////		PCASMSetOverlap(pc[levels-1],3);
////		PCASMSetTotalSubdomains(pc[levels-1], 32, NULL, NULL);
////		PCSetType(pc[levels-1],PCJACOBI);
//		KSPSetNormType(ksp[levels-1],KSP_NORM_NONE);
//		KSPSetTolerances(ksp[levels-1], 1.e-7, PETSC_DEFAULT, PETSC_DEFAULT, v[1]);
//		KSPSetFromOptions(ksp[levels-1]);
//	}
//
//	VecNorm(b[0], NORM_2, &bnorm);
//	
//	VecSet(u[0], 0.0); // Note: This should be moved out of this function?
//	VecDuplicate(b[0],&(rv[0]));
//	MatMult(A[0], u[0], rv[0]);
//	VecAXPY(rv[0], -1.0, b[0]);
//	VecNorm(rv[0], NORM_2, &rnormchk);
////	if (rank==0) rnorm[0] = rnormchk;
//	rnorm[0] = rnormchk;
//
//	iter = 0;
////	rnormchk = bnorm;
////	if (rank==0) rnorm[0] = 1.0;
//	
//	double initWallTime = MPI_Wtime();
//	clock_t solverInitT = clock();
//	PetscLogStageRegister("Solver", &stage);
//	PetscLogStagePush(stage);
//	while (iter<maxIter && 100000000*bnorm > rnormchk && rnormchk > (1.e-7)*bnorm) {
//		KSPSolve(ksp[0], b[0], u[0]);
//		if (iter==0) KSPSetInitialGuessNonzero(ksp[0],PETSC_TRUE);
//		for (int l=1;l<levels;l++) {
//			KSPBuildResidual(ksp[l-1],NULL,rv[l-1],&(r[l-1]));
//			MatMult(res[l-1],r[l-1],b[l]);
//			KSPSolve(ksp[l], b[l], u[l]);
//			if (l!=levels-1) KSPSetInitialGuessNonzero(ksp[l],PETSC_TRUE);
//		}
//		for (int l=levels-2;l>=0;l=l-1) {
//			MatMult(pro[l],u[l+1],rv[l]);
//			VecAXPY(u[l],1.0,rv[l]);
//			KSPSolve(ksp[l], b[l], u[l]);
//			if (l!=0) KSPSetInitialGuessNonzero(ksp[l],PETSC_FALSE);
//		}
//		KSPBuildResidual(ksp[0],NULL,rv[0],&(r[0]));
//		VecNorm(r[0], NORM_2, &rnormchk);	
//		iter = iter + 1;
////		if (rank==0) rnorm[iter] = rnormchk/bnorm;
//		rnorm[iter] = rnormchk;
//	}
//	PetscLogStagePop();
//	clock_t solverT = clock();
//	double endWallTime = MPI_Wtime();
//	rnormchk = rnorm[0];
//	for (int i=0;i<(maxIter+1);i++) {
//		rnorm[i] = rnorm[i]/rnormchk;
//	}
//	solver->numIter = iter;
//
//	for (int i=0;i<levels;i++) {
//		PetscPrintf(PETSC_COMM_WORLD,"---------------------------| level = %d |------------------------\n",i);
//		KSPView(ksp[i],PETSC_VIEWER_STDOUT_WORLD);
//		PetscPrintf(PETSC_COMM_WORLD,"-----------------------------------------------------------------\n");
//	}
//	for (int i=0;i<levels;i++) {
//		VecDestroy(&(rv[i]));
//	}
//	for (int i=0;i<levels;i++) {
//		KSPDestroy(&(ksp[i]));
//	}
//	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = [%d]; Solver cputime:                %lf\n",rank,(double)(solverT-solverInitT)/CLOCKS_PER_SEC);
//	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = [%d]; Solver walltime:               %lf\n",rank,endWallTime-initWallTime);
//	PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
	
	return 0;
}

//void MultigridVcycle(Solver *solver) {
//
//	int	iter;
//	double	rnormchk, bnorm;
//	
//	double	*rnorm;
//	int	maxIter;
//
//	int	*v;
//	int	levels;
//	Mat 	*res;
//	Mat 	*pro;
//	Mat	*A;
//	Vec	*b;
//	Vec	*u;
//	
//	int	size, rank;
//	
//	MPI_Comm_size(PETSC_COMM_WORLD, &size);
//	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
//	
//	maxIter = solver->numIter;	
//	rnorm	= solver->rnorm;
//	v	= solver->v;
//
//	levels	= solver->assem->levels;
//	res	= solver->assem->res;
//	pro	= solver->assem->pro;
//	A	= solver->assem->A;
//	b	= solver->assem->b;
//	u	= solver->assem->u;
//	
//	KSP	ksp[levels];
////	PC	pc[levels];
//	Vec	r[levels], rv[levels];//, xbuf[levels];
//	
//	PetscLogStage	stage;
//	
//	//printf("Enter the number of fine grid sweeps = ");
////	scanf("%d",v);
//	//printf("Enter the number of coarse grid sweeps = ");
////	scanf("%d",v+1);
//
////	PetscOptionsGetIntArray(NULL, NULL, "-v", v, &vmax, NULL);
//	
//	for (int i=0;i<levels;i++) {
//		VecDuplicate(b[i],&(rv[i]));
//	}
//	
//	KSPCreate(PETSC_COMM_WORLD, &(ksp[0]));
////	KSPSetType(ksp[0],KSPGMRES);
//	KSPSetType(ksp[0],KSPRICHARDSON);
////	KSPRichardsonSetScale(ksp[0],2.0/3.0);
//	KSPSetOperators(ksp[0], A[0], A[0]);
////	KSPGetPC(ksp[0],&(pc[0]));
////	PCSetType(pc[0],PCASM);
////	PCASMSetType(pc[0],PC_ASM_BASIC);
////	PCASMSetOverlap(pc[0],3);
////	PCASMSetTotalSubdomains(pc[0], 32, NULL, NULL);
////	PCSetType(pc[0],PCJACOBI);
//	KSPSetNormType(ksp[0],KSP_NORM_NONE);
//	KSPSetTolerances(ksp[0], 1.e-7, PETSC_DEFAULT, PETSC_DEFAULT, v[0]);
//	KSPSetFromOptions(ksp[0]);
//	
//	for (int i=1;i<levels-1;i++) {
//		KSPCreate(PETSC_COMM_WORLD, &(ksp[i]));
////		KSPSetType(ksp[i],KSPGMRES);
//		KSPSetType(ksp[i],KSPRICHARDSON);
////		KSPRichardsonSetScale(ksp[i],2.0/3.0);
//		KSPSetOperators(ksp[i], A[i], A[i]);
////		KSPGetPC(ksp[i],&(pc[i]));
////		PCSetType(pc[i],PCASM);
////		PCASMSetType(pc[i],PC_ASM_BASIC);
////		PCASMSetOverlap(pc[i],3);
////		PCASMSetTotalSubdomains(pc[i], 32, NULL, NULL);
////		PCSetType(pc[i],PCJACOBI);
//		KSPSetNormType(ksp[i],KSP_NORM_NONE);
//		KSPSetTolerances(ksp[i], 1.e-7, PETSC_DEFAULT, PETSC_DEFAULT, v[0]);
//		KSPSetFromOptions(ksp[i]);
//	}
//
//	if (levels>1) {
//		KSPCreate(PETSC_COMM_WORLD, &(ksp[levels-1]));
////		KSPSetType(ksp[levels-1],KSPGMRES);
//		KSPSetType(ksp[levels-1],KSPRICHARDSON);
////		KSPRichardsonSetScale(ksp[levels-1],2.0/3.0);
//		KSPSetOperators(ksp[levels-1], A[levels-1], A[levels-1]);
////		KSPGetPC(ksp[levels-1],&(pc[levels-1]));
////		PCSetType(pc[levels-1],PCASM);
////		PCASMSetType(pc[levels-1],PC_ASM_BASIC);
////		PCASMSetOverlap(pc[levels-1],3);
////		PCASMSetTotalSubdomains(pc[levels-1], 32, NULL, NULL);
////		PCSetType(pc[levels-1],PCJACOBI);
//		KSPSetNormType(ksp[levels-1],KSP_NORM_NONE);
//		KSPSetTolerances(ksp[levels-1], 1.e-7, PETSC_DEFAULT, PETSC_DEFAULT, v[1]);
//		KSPSetFromOptions(ksp[levels-1]);
//	}
//
//	VecNorm(b[0], NORM_2, &bnorm);
//	
//	VecSet(u[0], 0.0); // Note: This should be moved out of this function?
//	VecDuplicate(b[0],&(rv[0]));
//	MatMult(A[0], u[0], rv[0]);
//	VecAXPY(rv[0], -1.0, b[0]);
//	VecNorm(rv[0], NORM_2, &rnormchk);
////	if (rank==0) rnorm[0] = rnormchk;
//	rnorm[0] = rnormchk;
//
//	iter = 0;
////	rnormchk = bnorm;
////	if (rank==0) rnorm[0] = 1.0;
//	
//	double initWallTime = MPI_Wtime();
//	clock_t solverInitT = clock();
//	PetscLogStageRegister("Solver", &stage);
//	PetscLogStagePush(stage);
//	while (iter<maxIter && 100000000*bnorm > rnormchk && rnormchk > (1.e-7)*bnorm) {
//		KSPSolve(ksp[0], b[0], u[0]);
//		if (iter==0) KSPSetInitialGuessNonzero(ksp[0],PETSC_TRUE);
//		for (int l=1;l<levels;l++) {
//			KSPBuildResidual(ksp[l-1],NULL,rv[l-1],&(r[l-1]));
//			MatMult(res[l-1],r[l-1],b[l]);
//			KSPSolve(ksp[l], b[l], u[l]);
//			if (l!=levels-1) KSPSetInitialGuessNonzero(ksp[l],PETSC_TRUE);
//		}
//		for (int l=levels-2;l>=0;l=l-1) {
//			MatMult(pro[l],u[l+1],rv[l]);
//			VecAXPY(u[l],1.0,rv[l]);
//			KSPSolve(ksp[l], b[l], u[l]);
//			if (l!=0) KSPSetInitialGuessNonzero(ksp[l],PETSC_FALSE);
//		}
//		KSPBuildResidual(ksp[0],NULL,rv[0],&(r[0]));
//		VecNorm(r[0], NORM_2, &rnormchk);	
//		iter = iter + 1;
////		if (rank==0) rnorm[iter] = rnormchk/bnorm;
//		rnorm[iter] = rnormchk;
//	}
//	PetscLogStagePop();
//	clock_t solverT = clock();
//	double endWallTime = MPI_Wtime();
//	rnormchk = rnorm[0];
//	for (int i=0;i<(maxIter+1);i++) {
//		rnorm[i] = rnorm[i]/rnormchk;
//	}
//	solver->numIter = iter;
//
//	for (int i=0;i<levels;i++) {
//		PetscPrintf(PETSC_COMM_WORLD,"---------------------------| level = %d |------------------------\n",i);
//		KSPView(ksp[i],PETSC_VIEWER_STDOUT_WORLD);
//		PetscPrintf(PETSC_COMM_WORLD,"-----------------------------------------------------------------\n");
//	}
//	for (int i=0;i<levels;i++) {
//		VecDestroy(&(rv[i]));
//	}
//	for (int i=0;i<levels;i++) {
//		KSPDestroy(&(ksp[i]));
//	}
//	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = [%d]; Solver cputime:                %lf\n",rank,(double)(solverT-solverInitT)/CLOCKS_PER_SEC);
//	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = [%d]; Solver walltime:               %lf\n",rank,endWallTime-initWallTime);
//	PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
//
//}
//
//void MultigridAdditive(Solver *solver) {
//
//	int	iter;
//	double	rnormchk, bnorm;
//	
//	double	*rnorm;
//	int	maxIter;
//
//	int	*v;
//	int	levels;
//	Mat 	*res;
//	Mat 	*pro;
//	Mat	*A;
//	Vec	*b;
//	Vec	*u;
//	
//	int	size, rank;
//	
//	MPI_Comm_size(PETSC_COMM_WORLD, &size);
//	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
//	
//	maxIter = solver->numIter;	
//	rnorm	= solver->rnorm;
//	v	= solver->v;
//
//	levels	= solver->assem->levels;
//	res	= solver->assem->res;
//	pro	= solver->assem->pro;
//	A	= solver->assem->A;
//	b	= solver->assem->b;
//	u	= solver->assem->u;
//	
//	if (levels < 2) {ERROR_MSG("Cannot use Additive cycle for levels < 2; use I-cycle for levels = 1"); return;}
//
//	KSP	ksp[levels];
////	PC	pc[levels];
//	Vec	r;//, xbuf[levels];
//	
//	PetscLogStage	stage;
//		
//	VecDuplicate(b[0], &r);
//	
//	KSPCreate(PETSC_COMM_WORLD, &(ksp[0]));
//	PetscObjectSetOptionsPrefix(ksp[0], "fine_");
//	KSPSetType(ksp[0],KSPRICHARDSON);
//	KSPSetOperators(ksp[0], A[0], A[0]);
//	KSPSetNormType(ksp[0],KSP_NORM_NONE);
//	KSPSetTolerances(ksp[0], PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, v[0]);
//	KSPSetInitialGuessNonzero(ksp[0],PETSC_TRUE);
//	KSPSetFromOptions(ksp[0]);
//	
//	for (int i=1;i<levels-1;i++) {
//		KSPCreate(PETSC_COMM_WORLD, &(ksp[i]));
//		PetscObjectSetOptionsPrefix(ksp[i], "levels_");
//		KSPSetType(ksp[i],KSPRICHARDSON);
//		KSPSetOperators(ksp[i], A[i], A[i]);
//		KSPSetNormType(ksp[i],KSP_NORM_NONE);
//		KSPSetTolerances(ksp[i], PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, v[0]);
//		KSPSetFromOptions(ksp[i]);
//	}
//
//	KSPCreate(PETSC_COMM_WORLD, &(ksp[levels-1]));
//	PetscObjectSetOptionsPrefix(ksp[levels-1], "coarse_");
//	KSPSetType(ksp[levels-1],KSPRICHARDSON);
//	KSPSetOperators(ksp[levels-1], A[levels-1], A[levels-1]);
//	KSPSetNormType(ksp[levels-1],KSP_NORM_NONE);
//	KSPSetTolerances(ksp[levels-1], PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, v[1]);
//	KSPSetFromOptions(ksp[levels-1]);
//
//	VecNorm(b[0], NORM_2, &bnorm);
//	
////	VecSet(u[0], 0.0); // Note: Should this be moved out of this function?
//	for (int l=0; l<levels; l++) {
//		VecSet(u[l], 0.0); // Note: Should this be moved out of this function?
//	}
////	VecDuplicate(b[0],&(r[0]));
//	MatResidual(A[0], b[0], u[0], r);
////	MatMult(A[0], u[0], rv[0]);
////	VecAXPY(rv[0], -1.0, b[0]);
//	VecNorm(r, NORM_2, &rnormchk);
//	rnorm[0] = rnormchk;
//
//	iter = 0;
//	
//	double initWallTime = MPI_Wtime();
//	clock_t solverInitT = clock();
//	PetscLogStageRegister("Solver", &stage);
//	PetscLogStagePush(stage);
//	// Note: Following algorithm is written assuming zer initial solution on coarse grids
//	while (iter<maxIter && 100000000*bnorm > rnormchk && rnormchk > (1.e-7)*bnorm) {
//		MatMult(res[0], r, b[1]);
//		for (int l=1;l<levels-1;l++) {
//			MatMult(res[l], b[l], b[l+1]);
//		}
//		for (int l=0; l<levels; l++) {
//			KSPSolve(ksp[l], b[l], u[l]);
//		}
//		for (int l=levels-2;l>0;l=l-1) {
//			MatMult(pro[l], u[l+1], b[l]);
//			VecAXPY(u[l], 1.0, b[l]);
//		}
//		MatMult(pro[0], u[1], r);
//		VecAXPY(u[0], 1.0, r);
//		MatResidual(A[0], b[0], u[0], r);
//		VecNorm(r, NORM_2, &rnormchk);	
//		iter = iter + 1;
//		rnorm[iter] = rnormchk;
//	}
//	PetscLogStagePop();
//	clock_t solverT = clock();
//	double endWallTime = MPI_Wtime();
//	rnormchk = rnorm[0];
//	for (int i=0;i<(maxIter+1);i++) {
//		rnorm[i] = rnorm[i]/rnormchk;
//	}
//	solver->numIter = iter;
//
//	for (int i=0;i<levels;i++) {
//		PetscPrintf(PETSC_COMM_WORLD,"---------------------------| level = %d |------------------------\n",i);
//		KSPView(ksp[i],PETSC_VIEWER_STDOUT_WORLD);
//		PetscPrintf(PETSC_COMM_WORLD,"-----------------------------------------------------------------\n");
//	}
//	VecDestroy(&r);
//	for (int i=0;i<levels;i++) {
//		KSPDestroy(&(ksp[i]));
//	}
//	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = [%d]; Solver cputime:                %lf\n",rank,(double)(solverT-solverInitT)/CLOCKS_PER_SEC);
//	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = [%d]; Solver walltime:               %lf\n",rank,endWallTime-initWallTime);
//	PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
//
//}
//
//void MultigridAdditiveScaled(Solver *solver) {
//
//	int	iter;
//	double	rnormchk, bnorm;
//	
//	double	*rnorm;
//	int	maxIter;
//
//	int	*v;
//	int	levels;
//	Mat 	*res;
//	Mat 	*pro;
//	Mat	*A;
//	Vec	*b;
//	Vec	*u;
//	
//	int	size, rank;
//	
//	MPI_Comm_size(PETSC_COMM_WORLD, &size);
//	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
//	
//	maxIter = solver->numIter;	
//	rnorm	= solver->rnorm;
//	v	= solver->v;
//
//	levels	= solver->assem->levels;
//	res	= solver->assem->res;
//	pro	= solver->assem->pro;
//	A	= solver->assem->A;
//	b	= solver->assem->b;
//	u	= solver->assem->u;
//	
//	if (levels < 2) {ERROR_MSG("Cannot use Additive cycle for levels < 2; use I-cycle for levels = 1"); return;}
//
//	KSP	ksp[levels];
////	PC	pc[levels];
//	Vec	r0, r1[levels];
//	
//	PetscLogStage	stage;
//		
//	VecDuplicate(b[0], &(r0));
//	for (int i=0;i<levels;i++) {
//		VecDuplicate(b[i], &(r1[i]));
//	}
//	
//	KSPCreate(PETSC_COMM_WORLD, &(ksp[0]));
//	PetscObjectSetOptionsPrefix(ksp[0], "fine_");
//	KSPSetType(ksp[0],KSPRICHARDSON);
//	KSPSetOperators(ksp[0], A[0], A[0]);
//	KSPSetNormType(ksp[0],KSP_NORM_NONE);
//	KSPSetTolerances(ksp[0], PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, v[0]);
//	KSPSetInitialGuessNonzero(ksp[0],PETSC_TRUE);
//	KSPSetFromOptions(ksp[0]);
//	
//	for (int i=1;i<levels-1;i++) {
//		KSPCreate(PETSC_COMM_WORLD, &(ksp[i]));
//		PetscObjectSetOptionsPrefix(ksp[i], "levels_");
//		KSPSetType(ksp[i],KSPRICHARDSON);
//		KSPSetOperators(ksp[i], A[i], A[i]);
//		KSPSetNormType(ksp[i],KSP_NORM_NONE);
//		KSPSetTolerances(ksp[i], PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, v[0]);
//		KSPSetFromOptions(ksp[i]);
//	}
//
//	KSPCreate(PETSC_COMM_WORLD, &(ksp[levels-1]));
//	PetscObjectSetOptionsPrefix(ksp[levels-1], "coarse_");
//	KSPSetType(ksp[levels-1],KSPRICHARDSON);
//	KSPSetOperators(ksp[levels-1], A[levels-1], A[levels-1]);
//	KSPSetNormType(ksp[levels-1],KSP_NORM_NONE);
//	KSPSetTolerances(ksp[levels-1], PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, v[1]);
//	KSPSetFromOptions(ksp[levels-1]);
//
//	VecNorm(b[0], NORM_2, &bnorm);
//	
//	VecSet(u[0], 0.0); // Note: Should this be moved out of this function?
//	for (int l=1; l<levels; l++) {
//		VecSet(u[l], 0.0); // Note: Should this be moved out of this function?
//	}
////	VecDuplicate(b[0],&(r[0]));
//	MatResidual(A[0], b[0], u[0], r0);
////	MatMult(A[0], u[0], rv[0]);
////	VecAXPY(rv[0], -1.0, b[0]);
//	VecNorm(r0, NORM_2, &rnormchk);
//	rnorm[0] = rnormchk;
//
//	iter = 0;
//	double lambda[levels], r0Dot[levels];
//	
//	double initWallTime = MPI_Wtime();
//	clock_t solverInitT = clock();
//	PetscLogStageRegister("Solver", &stage);
//	PetscLogStagePush(stage);
//	while (iter<maxIter && 100000000*bnorm > rnormchk && rnormchk > (1.e-7)*bnorm) {
//		MatMult(res[0], r0, b[1]);
//		for (int l=1;l<levels-1;l++) {
//			MatMult(res[l], b[l], b[l+1]);
//		}
////		VecTDot(r0, r0, r0Dot);
//		for (int l=1; l<levels-1; l++) {
//			VecTDot(b[l], b[l], r0Dot+l);
//		}
//		#pragma omp parallel for
//		for (int l=0; l<levels; l++) {
////			int t_total = omp_get_num_threads();
////			int t_id = omp_get_thread_num();
////			printf("Thread %d out of %d is solving level %d\n", t_id, t_total, l);
//			KSPSolve(ksp[l], b[l], u[l]);
//		}
//		for (int l=0; l<levels-1; l++) {
//			MatResidual(A[l], b[l], u[l], r1[l]);
//		}
//		VecTDot(r0, r1[0], lambda);
//		lambda[0] = lambda[0]/(rnormchk*rnormchk);
//		for (int l=1; l<levels-1; l++) {
//			VecTDot(b[l], r1[l], lambda+l);
//			lambda[l] = lambda[l]/r0Dot[l];
//		}
//		for (int l=levels-2;l>=0;l=l-1) {
//			VecScale(u[l+1], lambda[l]);
//			MatMult(pro[l], u[l+1], r1[l]);
//			VecAXPY(u[l], 1.0, r1[l]);
//		}
//		MatResidual(A[0], b[0], u[0], r0);
//		VecNorm(r0, NORM_2, &rnormchk);	
//		iter = iter + 1;
//		rnorm[iter] = rnormchk;
//	}
//	PetscLogStagePop();
//	clock_t solverT = clock();
//	double endWallTime = MPI_Wtime();
//	rnormchk = rnorm[0];
//	for (int i=0;i<(maxIter+1);i++) {
//		rnorm[i] = rnorm[i]/rnormchk;
//	}
//	solver->numIter = iter;
//
//	for (int i=0;i<levels;i++) {
//		PetscPrintf(PETSC_COMM_WORLD,"---------------------------| level = %d |------------------------\n",i);
//		KSPView(ksp[i],PETSC_VIEWER_STDOUT_WORLD);
//		PetscPrintf(PETSC_COMM_WORLD,"-----------------------------------------------------------------\n");
//	}
//	VecDestroy(&(r0));
//	for (int i=0;i<levels;i++) {
//		VecDestroy(&(r1[i]));
//	}
//	for (int i=0;i<levels;i++) {
//		KSPDestroy(&(ksp[i]));
//	}
//	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = [%d]; Solver cputime:                %lf\n",rank,(double)(solverT-solverInitT)/CLOCKS_PER_SEC);
//	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = [%d]; Solver walltime:               %lf\n",rank,endWallTime-initWallTime);
//	PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
//
//}
//
//void MultigridFilter(Solver *solver) {
//
//	int	iter;
//	double	rnormchk, bnorm;
//	
//	double	*rnorm;
//	int	maxIter;
//
//	int	*v;
//	int	levels;
//	Mat 	*res;
//	Mat 	*pro;
//	Mat	*A;
//	Vec	*b;
//	Vec	*u;
//	
//	int	size, rank;
//	
//	MPI_Comm_size(PETSC_COMM_WORLD, &size);
//	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
//	
//	maxIter = solver->numIter;	
//	rnorm	= solver->rnorm;
//	v	= solver->v;
//
//	levels	= solver->assem->levels;
//	res	= solver->assem->res;
//	pro	= solver->assem->pro;
//	A	= solver->assem->A;
//	b	= solver->assem->b;
//	u	= solver->assem->u;
//	
//	if (levels < 2) {ERROR_MSG("Cannot use Filtered cycle for levels < 2; use I-cycle for levels = 1"); return;}
//
//	Mat	filter[levels-1];
//	
//	for (int l=0; l<levels-1; l++) {
//		MatMatMult(pro[l], res[l], MAT_INITIAL_MATRIX, PETSC_DEFAULT, filter+l);
////		MatView(filter[l], PETSC_VIEWER_STDOUT_WORLD);
//	}
//
//	KSP	ksp[levels];
////	PC	pc[levels];
//	Vec	r, e;//, xbuf[levels];
//	
//	PetscLogStage	stage;
//		
//	VecDuplicate(b[0], &r);
//	VecDuplicate(u[0], &e);
//	
//	KSPCreate(PETSC_COMM_WORLD, &(ksp[0]));
//	PetscObjectSetOptionsPrefix(ksp[0], "fine_");
//	KSPSetType(ksp[0],KSPRICHARDSON);
//	KSPSetOperators(ksp[0], A[0], A[0]);
//	KSPSetNormType(ksp[0],KSP_NORM_NONE);
//	KSPSetTolerances(ksp[0], PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, v[0]);
////	KSPSetInitialGuessNonzero(ksp[0],PETSC_TRUE);
//	KSPSetFromOptions(ksp[0]);
//	
//	for (int i=1;i<levels-1;i++) {
//		KSPCreate(PETSC_COMM_WORLD, &(ksp[i]));
//		PetscObjectSetOptionsPrefix(ksp[i], "levels_");
//		KSPSetType(ksp[i],KSPRICHARDSON);
//		KSPSetOperators(ksp[i], A[i], A[i]);
//		KSPSetNormType(ksp[i],KSP_NORM_NONE);
//		KSPSetTolerances(ksp[i], PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, v[0]);
//		KSPSetFromOptions(ksp[i]);
//	}
//
//	KSPCreate(PETSC_COMM_WORLD, &(ksp[levels-1]));
//	PetscObjectSetOptionsPrefix(ksp[levels-1], "coarse_");
//	KSPSetType(ksp[levels-1],KSPRICHARDSON);
//	KSPSetOperators(ksp[levels-1], A[levels-1], A[levels-1]);
//	KSPSetNormType(ksp[levels-1],KSP_NORM_NONE);
//	KSPSetTolerances(ksp[levels-1], PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, v[1]);
//	KSPSetFromOptions(ksp[levels-1]);
//
//	VecNorm(b[0], NORM_2, &bnorm);
//	
//	for (int i=0; i<levels; i++) {
//		VecSet(u[i], 0.0);
//	}
//	VecSet(r, 0.0);
//	VecSet(e, 0.0);
////	VecDuplicate(b[0],&(r[0]));
//	MatResidual(A[0], b[0], u[0], r);
////	MatMult(A[0], u[0], rv[0]);
////	VecAXPY(rv[0], -1.0, b[0]);
//	VecNorm(r, NORM_2, &rnormchk);
//	rnorm[0] = rnormchk;
//
//	iter = 0;
//	
//	double initWallTime = MPI_Wtime();
//	clock_t solverInitT = clock();
//	PetscLogStageRegister("Solver", &stage);
//	PetscLogStagePush(stage);
//	// Note: Following algorithm is written assuming initial solution 
//	// on coarse grids is always zero
//	while (iter<maxIter && 100000000*bnorm > rnormchk && rnormchk > (1.e-7)*bnorm) {
////		MatResidual(A[0], b[0], u[0], r);
//		MatMult(filter[0], r, e);
//		VecAXPY(r, -1.0, e);
//		MatMult(res[0], e, b[1]);
//		for (int l=1;l<levels-1;l++) {
//			MatMult(filter[l], b[l], u[l]);
//			VecAXPY(b[l], -1.0, u[l]);
//			MatMult(res[l], u[l], b[l+1]);
//		}
//		KSPSolve(ksp[0], r, e);
//		for (int l=1;l<levels;l++) {
//			KSPSolve(ksp[l], b[l], u[l]);
//		}
//		for (int l=levels-2;l>0;l=l-1) {
//			MatMult(pro[l], u[l+1], b[l]);
//			VecAXPY(u[l], 1.0, b[l]);
//		}
//		MatMult(pro[0],u[1],r);
//		VecAXPBYPCZ(u[0], 1.0, 1.0, 1.0, e, r);
//		MatResidual(A[0], b[0], u[0], r);
//		VecNorm(r, NORM_2, &rnormchk);	
//		iter = iter + 1;
//		rnorm[iter] = rnormchk;
//	}
//	PetscLogStagePop();
//	clock_t solverT = clock();
//	double endWallTime = MPI_Wtime();
//	rnormchk = rnorm[0];
//	for (int i=0;i<(maxIter+1);i++) {
//		rnorm[i] = rnorm[i]/rnormchk;
//	}
//	solver->numIter = iter;
//
//	for (int i=0;i<levels;i++) {
//		PetscPrintf(PETSC_COMM_WORLD,"---------------------------| level = %d |------------------------\n",i);
//		KSPView(ksp[i],PETSC_VIEWER_STDOUT_WORLD);
//		PetscPrintf(PETSC_COMM_WORLD,"-----------------------------------------------------------------\n");
//	}
//	for (int l=0; l<levels-1; l++) {
//		MatDestroy(filter+l);
//	}
//	VecDestroy(&r);
//	VecDestroy(&e);
//	for (int i=0;i<levels;i++) {
//		KSPDestroy(&(ksp[i]));
//	}
//	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = [%d]; Solver cputime:                %lf\n",rank,(double)(solverT-solverInitT)/CLOCKS_PER_SEC);
//	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = [%d]; Solver walltime:               %lf\n",rank,endWallTime-initWallTime);
//	PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
//
//}
//
//void MultigridVFilter(Solver *solver) {
//
//	int	iter;
//	double	rnormchk, bnorm;
//	
//	double	*rnorm;
//	int	maxIter;
//
//	int	*v;
//	int	levels;
//	Mat 	*res;
//	Mat 	*pro;
//	Mat	*A;
//	Vec	*b;
//	Vec	*u;
//	
//	int	size, rank;
//	
//	MPI_Comm_size(PETSC_COMM_WORLD, &size);
//	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
//	
//	maxIter = solver->numIter;	
//	rnorm	= solver->rnorm;
//	v	= solver->v;
//
//	levels	= solver->assem->levels;
//	res	= solver->assem->res;
//	pro	= solver->assem->pro;
//	A	= solver->assem->A;
//	b	= solver->assem->b;
//	u	= solver->assem->u;
//	
//	if (levels < 2) {ERROR_MSG("Cannot use V-Filter cycle for levels < 2; use I-cycle for levels = 1"); return;}
//
//	Mat	filter[levels-1];
//	
//	for (int l=0; l<levels-1; l++) {
//		MatMatMult(pro[l], res[l], MAT_INITIAL_MATRIX, PETSC_DEFAULT, filter+l);
////		MatView(filter[l], PETSC_VIEWER_STDOUT_WORLD);
//	}
//
//	KSP	ksp[levels];
////	PC	pc[levels];
//	Vec	r[levels], e[levels];//, xbuf[levels];
//	
//	PetscLogStage	stage;
//		
//	for (int i=0;i<levels;i++) {
//		VecDuplicate(b[i], &(r[i]));
//		VecDuplicate(u[i], &(e[i]));
//	}
//	
//	KSPCreate(PETSC_COMM_WORLD, &(ksp[0]));
//	PetscObjectSetOptionsPrefix(ksp[0], "fine_");
//	KSPSetType(ksp[0],KSPRICHARDSON);
//	KSPSetOperators(ksp[0], A[0], A[0]);
//	KSPSetNormType(ksp[0],KSP_NORM_NONE);
//	KSPSetTolerances(ksp[0], PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, v[0]);
//	KSPSetInitialGuessNonzero(ksp[0],PETSC_TRUE);
//	KSPSetFromOptions(ksp[0]);
//	
//	for (int i=1;i<levels-1;i++) {
//		KSPCreate(PETSC_COMM_WORLD, &(ksp[i]));
//		PetscObjectSetOptionsPrefix(ksp[i], "levels_");
//		KSPSetType(ksp[i],KSPRICHARDSON);
//		KSPSetOperators(ksp[i], A[i], A[i]);
//		KSPSetNormType(ksp[i],KSP_NORM_NONE);
//		KSPSetTolerances(ksp[i], PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, v[0]);
//		KSPSetFromOptions(ksp[i]);
//	}
//
//	KSPCreate(PETSC_COMM_WORLD, &(ksp[levels-1]));
//	PetscObjectSetOptionsPrefix(ksp[levels-1], "coarse_");
//	KSPSetType(ksp[levels-1],KSPRICHARDSON);
//	KSPSetOperators(ksp[levels-1], A[levels-1], A[levels-1]);
//	KSPSetNormType(ksp[levels-1],KSP_NORM_NONE);
//	KSPSetTolerances(ksp[levels-1], PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, v[1]);
//	KSPSetFromOptions(ksp[levels-1]);
//
//	VecNorm(b[0], NORM_2, &bnorm);
//	
//	VecSet(u[0], 0.0); // Note: Should this be moved out of this function?
////	VecDuplicate(b[0],&(r[0]));
//	MatResidual(A[0], b[0], u[0], r[0]);
////	MatMult(A[0], u[0], rv[0]);
////	VecAXPY(rv[0], -1.0, b[0]);
//	VecNorm(r[0], NORM_2, &rnormchk);
//	rnorm[0] = rnormchk;
//
//	iter = 0;
//	
//	double initWallTime = MPI_Wtime();
//	clock_t solverInitT = clock();
//	PetscLogStageRegister("Solver", &stage);
//	PetscLogStagePush(stage);
//	while (iter<maxIter && 100000000*bnorm > rnormchk && rnormchk > (1.e-7)*bnorm) {
//		KSPSolve(ksp[0], b[0], u[0]);
//		for (int l=1;l<levels-1;l++) {
////			KSPBuildResidual(ksp[l-1],NULL,rv[l-1],&(r[l-1]));
//			MatResidual(A[l-1], b[l-1], u[l-1], r[l-1]);
//			MatMult(filter[l-1], r[l-1], e[l-1]);
//			VecAXPY(r[l-1], -1.0, e[l-1]);
//			MatMult(res[l-1], e[l-1], b[l]);
//			VecSet(e[l-1], 0.0);
//			KSPSolve(ksp[l-1], r[l-1], e[l-1]);
//			KSPSolve(ksp[l], b[l], u[l]);
//			KSPSetInitialGuessNonzero(ksp[l], PETSC_TRUE);
//		}
//		MatResidual(A[levels-2], b[levels-2], u[levels-2], r[levels-2]);
//		MatMult(filter[levels-2], r[levels-2], e[levels-2]);
//		VecAXPY(r[levels-2], -1.0, e[levels-2]);
//		MatMult(res[levels-2], e[levels-2], b[levels-1]);
//		VecSet(e[levels-2], 0.0);
//		KSPSolve(ksp[levels-2], r[levels-2], e[levels-2]);
//		KSPSolve(ksp[levels-1], b[levels-1], u[levels-1]);
//		for (int l=levels-2;l>0;l=l-1) {
//			MatMult(pro[l], u[l+1], r[l]);
//			VecAXPBYPCZ(u[l], 1.0, 1.0, 1.0, e[l], r[l]);
////			VecAXPY(u[l], 1.0, r[l]);
//			KSPSolve(ksp[l], b[l], u[l]);
//			KSPSetInitialGuessNonzero(ksp[l], PETSC_FALSE);
//		}
//		MatMult(pro[0],u[1],r[0]);
//		VecAXPBYPCZ(u[0], 1.0, 1.0, 1.0, e[0], r[0]);
////		VecAXPY(u[0],1.0,r[0]);
//		KSPSolve(ksp[0], b[0], u[0]);
////		KSPBuildResidual(ksp[0],NULL,rv[0],&(r[0]));
//		MatResidual(A[0], b[0], u[0], r[0]);
//		VecNorm(r[0], NORM_2, &rnormchk);	
//		iter = iter + 1;
//		rnorm[iter] = rnormchk;
//	}
//	PetscLogStagePop();
//	clock_t solverT = clock();
//	double endWallTime = MPI_Wtime();
//	rnormchk = rnorm[0];
//	for (int i=0;i<(maxIter+1);i++) {
//		rnorm[i] = rnorm[i]/rnormchk;
//	}
//	solver->numIter = iter;
//
//	for (int i=0;i<levels;i++) {
//		PetscPrintf(PETSC_COMM_WORLD,"---------------------------| level = %d |------------------------\n",i);
//		KSPView(ksp[i],PETSC_VIEWER_STDOUT_WORLD);
//		PetscPrintf(PETSC_COMM_WORLD,"-----------------------------------------------------------------\n");
//	}
//	for (int l=0; l<levels-1; l++) {
//		MatDestroy(filter+l);
//	}
//	for (int i=0;i<levels;i++) {
//		VecDestroy(&(r[i]));
//		VecDestroy(&(e[i]));
//	}
//	for (int i=0;i<levels;i++) {
//		KSPDestroy(&(ksp[i]));
//	}
//	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = [%d]; Solver cputime:                %lf\n",rank,(double)(solverT-solverInitT)/CLOCKS_PER_SEC);
//	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = [%d]; Solver walltime:               %lf\n",rank,endWallTime-initWallTime);
//	PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
//
//}
//
//void MultigridPetscPCMG(Solver *solver) {
//
//	int	iter;
//	double	rnormchk, bnorm;
//	
//	double	*rnorm;
//	int	maxIter;
//
//	int	*v;
//	int	levels;
//	Mat 	*res;
//	Mat 	*pro;
//	Mat	*A;
//	Vec	*b;
//	Vec	*u;
//	
//	int	size, rank;
//	
//	MPI_Comm_size(PETSC_COMM_WORLD, &size);
//	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
//	
//	maxIter = solver->numIter;	
//	rnorm	= solver->rnorm;
//	v	= solver->v;
//
//	levels	= solver->assem->levels;
//	res	= solver->assem->res;
//	pro	= solver->assem->pro;
//	A	= solver->assem->A;
//	b	= solver->assem->b;
//	u	= solver->assem->u;
//
//	KSP	ksp, kspTemp;
//	PC	pc;
//	
//	KSPCreate(PETSC_COMM_WORLD, &ksp);
//	KSPSetType(ksp, KSPRICHARDSON);
//	KSPSetOperators(ksp, A[0], A[0]);
//	KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED);
//	KSPSetResidualHistory(ksp, solver->rnorm, solver->numIter, PETSC_FALSE);
//	KSPSetTolerances(ksp, 1.e-7, PETSC_DEFAULT, PETSC_DEFAULT, solver->numIter);
//	KSPGetPC(ksp, &pc);
//	PCSetType(pc, PCMG);
//	PCMGSetLevels(pc, levels, NULL);
//	PCMGSetNumberSmooth(pc, v[0]);
//
//	PCMGGetCoarseSolve(pc, &kspTemp);
//	KSPSetType(kspTemp, KSPRICHARDSON);
//	KSPSetOperators(kspTemp, A[levels-1], A[levels-1]);
//	KSPSetTolerances(kspTemp, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, v[1]);
//
//	for (int i=1; i<levels; i++) {
//		PCMGGetSmoother(pc, i, &kspTemp);
//		KSPSetType(kspTemp, KSPRICHARDSON);
//		KSPSetOperators(kspTemp, A[levels-i-1], A[levels-i-1]);
//		PCMGSetInterpolation(pc, i, pro[levels-i-1]);
//		PCMGSetRestriction(pc, i, res[levels-i-1]);
//	}
//
//	Vec	r[levels];
//	
//	for (int i=0;i<levels;i++) {
//		VecDuplicate(b[i],&(r[i]));
//	}
//
//	PCMGSetR(pc, levels-1, r[0]);
//	for (int i=1; i<levels-1; i++) {
//		PCMGSetRhs(pc, i, b[levels-i-1]);
//		PCMGSetX(pc, i, u[levels-i-1]);
//		PCMGSetR(pc, i, r[levels-i-1]);
//	}
//	PCMGSetRhs(pc, 0, b[levels-1]);
//	PCMGSetX(pc, 0, u[levels-1]);
//	
//	KSPSetFromOptions(ksp);
//	
//	PetscLogStage	stage;
//	
//	double initWallTime = MPI_Wtime();
//	clock_t solverInitT = clock();
//	PetscLogStageRegister("Solver", &stage);
//	PetscLogStagePush(stage);
//	
//	KSPSolve(ksp, b[0], u[0]);
//	
//	PetscLogStagePop();
//	clock_t solverT = clock();
//	double endWallTime = MPI_Wtime();
//
//	KSPGetIterationNumber(ksp, &(solver->numIter));
//
//	double	rnorm0 = solver->rnorm[0];
//	for (int i=0;i<(solver->numIter+1);i++) {
//		solver->rnorm[i] = solver->rnorm[i]/rnorm0;
//	}
//
//	KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
//	
//	for (int i=0;i<levels;i++) {
//		VecDestroy(&(r[i]));
//	}
//	KSPDestroy(&ksp);
//
//	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = [%d]; Solver cputime:                %lf\n",rank,(double)(solverT-solverInitT)/CLOCKS_PER_SEC);
//	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = [%d]; Solver walltime:               %lf\n",rank,endWallTime-initWallTime);
//	PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
//
//}
//
//void MultigridIcycle(Solver *solver) {
//
//	Mat	*A;
//	Vec	*b;
//	Vec	*u;
//
//	int	size, rank;
//
//	MPI_Comm_size(PETSC_COMM_WORLD, &size);
//	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
//	
//	A	= solver->assem->A;
//	b	= solver->assem->b;
//	u	= solver->assem->u;
//	
//	KSP	ksp;
//	PC	pc;
//	
//	PetscLogStage	stage, stageSolve;
//	
//	KSPCreate(PETSC_COMM_WORLD, &ksp);
//	KSPSetType(ksp,KSPRICHARDSON);
//	KSPSetOperators(ksp, *A, *A);
//	KSPGetPC(ksp,&pc);
////	PCSetType(pc,PCASM);
//	KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED);
//	KSPSetResidualHistory(ksp, solver->rnorm, solver->numIter, PETSC_FALSE);
////	KSPMonitorSet(ksp, myMonitor, solver->rnorm, NULL);
//	KSPSetTolerances(ksp, 1.e-7, PETSC_DEFAULT, PETSC_DEFAULT, solver->numIter);
//	KSPSetFromOptions(ksp);
//	
//	// Compute initial residual and its norm
////	Vec	r;
////	
////	VecSet(*u, 0.0); // Note: This should be moved out of this function?
////	VecDuplicate(*b,&r);
////	MatMult(*A, *u, r);
////	VecAXPY(r, -1.0, *b);
////	VecNorm(r, NORM_2, &rnorm0);
////	VecDestroy(&r);
//	
//	// Solve the system
//	double initWallTime = MPI_Wtime();
//	clock_t solverInitT = clock();
//	PetscLogStageRegister("Solver", &stageSolve);
//	PetscLogStagePush(stageSolve);
//	
//	KSPSolve(ksp, *b, *u);
//	
//	PetscLogStagePop();
//	clock_t solverT = clock();
//	double endWallTime = MPI_Wtime();
////	KSPGetResidualHistory(ksp, NULL, &(solver->numIter));
//	KSPGetIterationNumber(ksp, &(solver->numIter));
//
//	double	rnorm0 = solver->rnorm[0];
//	for (int i=0;i<(solver->numIter+1);i++) {
//		solver->rnorm[i] = solver->rnorm[i]/rnorm0;
//	}
//
////	VecView(*u,PETSC_VIEWER_STDOUT_WORLD);
//	PetscPrintf(PETSC_COMM_WORLD,"---------------------------| level = 0 |------------------------\n");
//	KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
//	PetscPrintf(PETSC_COMM_WORLD,"----------------------------------------------------------------\n");
//	KSPDestroy(&ksp);
//
//	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = [%d]; Solver cputime:                %lf\n",rank,(double)(solverT-solverInitT)/CLOCKS_PER_SEC);
//	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = [%d]; Solver walltime:               %lf\n",rank,endWallTime-initWallTime);
//	PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
//}
//
//void MultigridEcycle(Solver *solver) {
//
//	Mat	*A1;
//	Mat	*A2;
//	Vec	r;
//	Vec	*b;
//	Vec	*u;
//	
//	double	*norm, chkNorm, bnorm;
//	int	v, maxIter, iter;
//	
//	int	size, rank;
//	
//	MPI_Comm_size(PETSC_COMM_WORLD, &size);
//	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
//	
//	A1	= solver->assem->A;
//	A2	= solver->assem->A2;
//	b	= solver->assem->b;
//	u	= solver->assem->u;
//	maxIter	= solver->numIter;
//	norm	= solver->rnorm;
//	v	= solver->v[0];
////	v	= 1;
//
//	KSP	ksp;
//	PC	pc;
//	
//	PetscLogStage	stage, stageSolve;
//	
//	KSPCreate(PETSC_COMM_WORLD, &ksp);
//	KSPSetType(ksp,KSPRICHARDSON);
//	KSPSetOperators(ksp, *A1, *A1);
//	KSPGetPC(ksp,&pc);
////	PCSetType(pc,PCASM);
////	KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED);
////	KSPMonitorSet(ksp, myMonitor, norm, NULL);
//	KSPSetNormType(ksp,KSP_NORM_NONE);
//	KSPSetTolerances(ksp, 1.e-7, PETSC_DEFAULT, PETSC_DEFAULT, v);
//	KSPSetFromOptions(ksp);
//	KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);
//	
//	MatScale(*A2, -1);
//	VecDuplicate(*b, &r);
////	VecCopy(*b, *r);
//	VecSet(*u, 0); // Note: This should be moved out of this function?
//	VecNorm(*b, NORM_2, &bnorm);
//	
//	MatMult(*A1, *u, r);
//	VecAXPY(r, -1.0, *b);
//	VecNorm(r, NORM_2, &chkNorm);
//	norm[0] = chkNorm;
//	
//	iter = 0;
////	chkNorm = bnorm;
//	double initWallTime = MPI_Wtime();
//	clock_t solverInitT = clock();
//	PetscLogStageRegister("Solver", &stageSolve);
//	PetscLogStagePush(stageSolve);
//	while (iter<maxIter && 100000000*bnorm > chkNorm && chkNorm > (1.e-7)*bnorm) {
//		MatMultAdd(*A2, *u, *b, r);
//		KSPSolve(ksp, r, *u);
//		iter += 1;
//		
//		MatMult(*A1, *u, r);
//		VecAXPY(r, -1.0, *b);
//		VecNorm(r, NORM_2, &chkNorm);
//		norm[iter] = chkNorm;
//	}
//	
//	PetscLogStagePop();
//	clock_t solverT = clock();
//	double endWallTime = MPI_Wtime();
//	chkNorm = norm[0];
//	for (int i=0;i<(maxIter+1);i++) {
//		norm[i] = norm[i]/chkNorm;
//	}
//	solver->numIter = iter;
////	KSPGetIterationNumber(ksp, &(solver->numIter));
//
////	VecView(*u,PETSC_VIEWER_STDOUT_WORLD);
//	PetscPrintf(PETSC_COMM_WORLD,"---------------------------| level = 0 |------------------------\n");
//	KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
//	PetscPrintf(PETSC_COMM_WORLD,"----------------------------------------------------------------\n");
//	VecDestroy(&r);
//	KSPDestroy(&ksp);
//
//	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = [%d]; Solver cputime:                %lf\n",rank,(double)(solverT-solverInitT)/CLOCKS_PER_SEC);
//	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = [%d]; Solver walltime:               %lf\n",rank,endWallTime-initWallTime);
//	PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
//}
//
//void MultigridD2cycle(Solver *solver) {
//	
//	Mat	*res;
//	Mat	*pro;
//	Mat	*A;
//	Vec	*b;
//	Vec	*u;
//	IS	*botIS;
//	IS	*topIS;
////	IS	*subFineIS;
//	IS	**gridIS;
//
//	Vec	r;
//	Vec	bBot;
//	Vec	rTop;
//	Vec	uTop;
//	
//	double	*norm, chkNorm, bnorm;
//	int	v, maxIter, iter;
//	
//	int	size, rank;
//	
//	MPI_Comm_size(PETSC_COMM_WORLD, &size);
//	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
//	
//	res		= solver->assem->res;
//	pro		= solver->assem->pro;
//	A		= solver->assem->A;
//	b		= solver->assem->b;
//	u		= solver->assem->u;
//	botIS		= solver->assem->bottomIS;
//	topIS		= solver->assem->topIS;
////	subFineIS	= solver->assem->subFineIS;
//	gridIS		= solver->assem->gridIS;
//	maxIter		= solver->numIter;
//	norm		= solver->rnorm;
//	v		= solver->v[0];
//
//	KSP	ksp;
//	PC	pc;
//	
////	ISView(*topIS, PETSC_VIEWER_STDOUT_WORLD);
////	ISView(*subFineIS, PETSC_VIEWER_STDOUT_WORLD);
//	
////	for (int i=0; i<solver->grids; i++) {
////		ISView(gridIS[0][i], PETSC_VIEWER_STDOUT_WORLD);
////	}
//	
//	D1cntx	info;
//	if (solver->moreInfo != 0) {
//		info.innerCount = 0;
//		VecDuplicate(*b, &(info.rInner));
//		info.grids = solver->grids;
//		info.rGrid = malloc(info.grids*sizeof(Vec));
//		info.rNormGrid = solver->rNormGrid;
//		for (int i=0; i<solver->grids; i++) {
//			VecGetSubVector(info.rInner, gridIS[0][i], info.rGrid+i);
//		}
//	}
//
//	PetscLogStage	stage, stageSolve;
//	
//	KSPCreate(PETSC_COMM_WORLD, &ksp);
//	KSPSetType(ksp,KSPRICHARDSON);
//	KSPSetOperators(ksp, *A, *A);
//	KSPGetPC(ksp,&pc);
////	PCSetType(pc,PCASM);
//	if (solver->moreInfo == 0) {
//		KSPSetNormType(ksp,KSP_NORM_NONE);
//	} else {
//		KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED);
//		KSPSetResidualHistory(ksp, solver->rNormGlobal, maxIter*(v+1), PETSC_FALSE);
//		KSPMonitorSet(ksp, rNormGridMonitor, &info, NULL);
//	}
//	KSPSetTolerances(ksp, 1.e-16, PETSC_DEFAULT, PETSC_DEFAULT, v);
//	KSPSetFromOptions(ksp);
//	KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);
//	
//	VecDuplicate(*b, &r);
//	VecGetSubVector(*b, *botIS, &bBot);
//	VecGetSubVector(r, *topIS, &rTop);
//	VecGetSubVector(*u, *topIS, &uTop);
//
//	VecNorm(*b, NORM_2, &bnorm);
//	VecSet(*u, 0); // Note: This should be moved out of this function?
//	
//	Vec	residual;
//
//	MatMult(*A, *u, r);
//	VecAYPX(r, -1.0, *b);
//	VecNorm(r, NORM_2, &chkNorm);
//	norm[0] = chkNorm;
//	
//	iter = 0;
//	double initWallTime = MPI_Wtime();
//	clock_t solverInitT = clock();
//	PetscLogStageRegister("Solver", &stageSolve);
//	PetscLogStagePush(stageSolve);
//	while (iter<maxIter && 100000000*bnorm > chkNorm && chkNorm > (1.e-7)*bnorm) {
//		MatMult(*res, rTop, bBot);
//		KSPSolve(ksp, *b, *u);
//		MatMult(*pro, *u, rTop);
//		VecAXPY(uTop, 1.0, rTop);
//		KSPBuildResidual(ksp, NULL, r, &residual);
//		VecNorm(r, NORM_2, &chkNorm);
//		iter += 1;
//		norm[iter] = chkNorm;
//	}
//	
//	PetscLogStagePop();
//	clock_t solverT = clock();
//	double endWallTime = MPI_Wtime();
//	VecRestoreSubVector(*b, *botIS, &bBot);
//	VecRestoreSubVector(r, *topIS, &rTop);
//	VecRestoreSubVector(*u, *topIS, &uTop);
//	if (solver->moreInfo != 0) {
//		for (int i=0; i<info.grids; i++) {
////			VecView(uGrid[i],PETSC_VIEWER_STDOUT_WORLD);
//			VecRestoreSubVector(info.rInner, gridIS[0][i], info.rGrid+i);
//		}
//		free(info.rGrid);
//		VecDestroy(&(info.rInner));
//	}
//
//	solver->numIter = iter;
//	chkNorm = norm[0];
//	for (int i=0;i<solver->numIter+1;i++) {
//		norm[i] = norm[i]/chkNorm;
//	}
//	if (solver->moreInfo != 0) {
//		chkNorm = solver->rNormGlobal[0];
//		for (int i=0;i<solver->numIter*(v+1);i++) {
//			solver->rNormGlobal[i] = solver->rNormGlobal[i]/chkNorm;
//		}
//		for (int g=0; g<solver->grids; g++) {
//			chkNorm = solver->rNormGrid[g][0];
//			for (int i=0;i<solver->numIter*(v+1);i++) {
//				solver->rNormGrid[g][i] = solver->rNormGrid[g][i]/chkNorm;
//			}
//		}
//	}
//
//	PetscPrintf(PETSC_COMM_WORLD,"---------------------------| level = 0 |------------------------\n");
//	KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
//	PetscPrintf(PETSC_COMM_WORLD,"----------------------------------------------------------------\n");
//	VecDestroy(&r);
//	KSPDestroy(&ksp);
//
//	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = [%d]; Solver cputime:                %lf\n",rank,(double)(solverT-solverInitT)/CLOCKS_PER_SEC);
//	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = [%d]; Solver walltime:               %lf\n",rank,endWallTime-initWallTime);
//	PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
//}
//
//void MultigridD1PScycle(Solver *solver) {
//	
//	Mat	*res;
//	Mat	*pro;
//	Mat	*A;
//	Vec	*b;
//	Vec	*u;
//	IS	*botIS;
//	IS	*topIS;
////	IS	*subFineIS;
//	IS	**gridIS;
//
//	Vec	r;
//	Vec	bBot;
//	Vec	rTop;
//	Vec	uTop;
//	Vec	cTop;
//	
//	double	*norm, chkNorm, bnorm;
//	int	v, maxIter, iter;
//	
//	int	size, rank;
//	
//	MPI_Comm_size(PETSC_COMM_WORLD, &size);
//	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
//	
//	res		= solver->assem->res;
//	pro		= solver->assem->pro;
//	A		= solver->assem->A;
//	b		= solver->assem->b;
//	u		= solver->assem->u;
//	botIS		= solver->assem->bottomIS;
//	topIS		= solver->assem->topIS;
////	subFineIS	= solver->assem->subFineIS;
//	gridIS		= solver->assem->gridIS;
//	maxIter		= solver->numIter;
//	norm		= solver->rnorm;
//	v		= solver->v[0];
//
//	KSP	ksp;
//	PC	pc;
//	
////	ISView(*topIS, PETSC_VIEWER_STDOUT_WORLD);
////	ISView(*subFineIS, PETSC_VIEWER_STDOUT_WORLD);
//	
////	for (int i=0; i<solver->grids; i++) {
////		ISView(gridIS[0][i], PETSC_VIEWER_STDOUT_WORLD);
////	}
//	
//	D1cntx	info;
//	if (solver->moreInfo != 0) {
//		info.innerCount = 0;
//		VecDuplicate(*b, &(info.rInner));
//		info.grids = solver->grids;
//		info.rGrid = malloc(info.grids*sizeof(Vec));
//		info.rNormGrid = solver->rNormGrid;
//		for (int i=0; i<solver->grids; i++) {
//			VecGetSubVector(info.rInner, gridIS[0][i], info.rGrid+i);
//		}
//	}
//
//	PetscLogStage	stage, stageSolve;
//	
//	KSPCreate(PETSC_COMM_WORLD, &ksp);
//	KSPSetType(ksp,KSPRICHARDSON);
//	KSPSetOperators(ksp, *A, *A);
//	KSPGetPC(ksp,&pc);
////	PCSetType(pc,PCASM);
//	if (solver->moreInfo == 0) {
//		KSPSetNormType(ksp,KSP_NORM_NONE);
//	} else {
//		KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED);
//		KSPSetResidualHistory(ksp, solver->rNormGlobal, maxIter*2*(v+1), PETSC_FALSE);
//		KSPMonitorSet(ksp, rNormGridMonitor, &info, NULL);
//	}
//	KSPSetTolerances(ksp, 1.e-16, PETSC_DEFAULT, PETSC_DEFAULT, v);
//	KSPSetFromOptions(ksp);
//	KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);
//	
//	VecDuplicate(*b, &r);
//	VecGetSubVector(*b, *botIS, &bBot);
//	VecGetSubVector(r, *topIS, &rTop);
//	VecGetSubVector(*u, *topIS, &uTop);
//	VecDuplicate(uTop, &cTop);	
//
//	VecNorm(*b, NORM_2, &bnorm);
//	VecSet(*u, 0); // Note: This should be moved out of this function?
//	
//	Vec	residual;
//
//	MatMult(*A, *u, r);
//	VecAYPX(r, -1.0, *b);
//	VecNorm(r, NORM_2, &chkNorm);
//	norm[0] = chkNorm;
//	
//	iter = 0;
//	double initWallTime = MPI_Wtime();
//	clock_t solverInitT = clock();
//	PetscLogStageRegister("Solver", &stageSolve);
//	PetscLogStagePush(stageSolve);
//	while (iter<maxIter && 100000000*bnorm > chkNorm && chkNorm > (1.e-7)*bnorm) {
//		MatMult(*pro, *u, cTop);
//		VecAXPY(uTop, 1.0, cTop);
//		KSPSolve(ksp, *b, *u);
//		MatMult(*res, rTop, bBot);
//		KSPSolve(ksp, *b, *u);
//		KSPBuildResidual(ksp, NULL, r, &residual);
//		VecNorm(r, NORM_2, &chkNorm);
//		iter += 1;
//		norm[iter] = chkNorm;
//	}
//	
//	PetscLogStagePop();
//	clock_t solverT = clock();
//	double endWallTime = MPI_Wtime();
//	VecRestoreSubVector(*b, *botIS, &bBot);
//	VecRestoreSubVector(r, *topIS, &rTop);
//	VecRestoreSubVector(*u, *topIS, &uTop);
//	if (solver->moreInfo != 0) {
//		for (int i=0; i<info.grids; i++) {
////			VecView(uGrid[i],PETSC_VIEWER_STDOUT_WORLD);
//			VecRestoreSubVector(info.rInner, gridIS[0][i], info.rGrid+i);
//		}
//		free(info.rGrid);
//		VecDestroy(&(info.rInner));
//	}
//
//	solver->numIter = iter;
//	chkNorm = norm[0];
//	for (int i=0;i<solver->numIter+1;i++) {
//		norm[i] = norm[i]/chkNorm;
//	}
//	if (solver->moreInfo != 0) {
//		chkNorm = solver->rNormGlobal[0];
//		for (int i=0;i<solver->numIter*(v+1);i++) {
//			solver->rNormGlobal[i] = solver->rNormGlobal[i]/chkNorm;
//		}
//		for (int g=0; g<solver->grids; g++) {
//			chkNorm = solver->rNormGrid[g][0];
//			for (int i=0;i<solver->numIter*(v+1);i++) {
//				solver->rNormGrid[g][i] = solver->rNormGrid[g][i]/chkNorm;
//			}
//		}
//	}
//
//	PetscPrintf(PETSC_COMM_WORLD,"---------------------------| level = 0 |------------------------\n");
//	KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
//	PetscPrintf(PETSC_COMM_WORLD,"----------------------------------------------------------------\n");
//	VecDestroy(&cTop);
//	VecDestroy(&r);
//	KSPDestroy(&ksp);
//
//	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = [%d]; Solver cputime:                %lf\n",rank,(double)(solverT-solverInitT)/CLOCKS_PER_SEC);
//	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = [%d]; Solver walltime:               %lf\n",rank,endWallTime-initWallTime);
//	PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
//}
//
//void MultigridD1cycle(Solver *solver) {
//	
//	Mat	*res;
//	Mat	*pro;
//	Mat	*A;
//	Vec	*b;
//	Vec	*u;
//	IS	*botIS;
//	IS	*topIS;
////	IS	*subFineIS;
//	IS	**gridIS;
//
//	Vec	r;
//	Vec	bBot;
//	Vec	rTop;
//	Vec	uTop;
//	
//	double	*norm, chkNorm, bnorm;
//	int	v, maxIter, iter;
//	
//	int	size, rank;
//	
//	MPI_Comm_size(PETSC_COMM_WORLD, &size);
//	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
//	
//	res		= solver->assem->res;
//	pro		= solver->assem->pro;
//	A		= solver->assem->A;
//	b		= solver->assem->b;
//	u		= solver->assem->u;
//	botIS		= solver->assem->bottomIS;
//	topIS		= solver->assem->topIS;
////	subFineIS	= solver->assem->subFineIS;
//	gridIS		= solver->assem->gridIS;
//	maxIter		= solver->numIter;
//	norm		= solver->rnorm;
//	v		= solver->v[0];
//
//	KSP	ksp;
//	PC	pc;
//	
////	ISView(*topIS, PETSC_VIEWER_STDOUT_WORLD);
////	ISView(*subFineIS, PETSC_VIEWER_STDOUT_WORLD);
//	
////	for (int i=0; i<solver->grids; i++) {
////		ISView(gridIS[0][i], PETSC_VIEWER_STDOUT_WORLD);
////	}
//	
//	D1cntx	info;
//	if (solver->moreInfo != 0) {
//		info.innerCount = 0;
//		VecDuplicate(*b, &(info.rInner));
//		info.grids = solver->grids;
//		info.rGrid = malloc(info.grids*sizeof(Vec));
//		info.rNormGrid = solver->rNormGrid;
//		for (int i=0; i<solver->grids; i++) {
//			VecGetSubVector(info.rInner, gridIS[0][i], info.rGrid+i);
//		}
//	}
//
//	PetscLogStage	stage, stageSolve;
//	
//	KSPCreate(PETSC_COMM_WORLD, &ksp);
//	KSPSetType(ksp,KSPRICHARDSON);
//	KSPSetOperators(ksp, *A, *A);
//	KSPGetPC(ksp,&pc);
////	PCSetType(pc,PCASM);
//	if (solver->moreInfo == 0) {
//		KSPSetNormType(ksp,KSP_NORM_NONE);
//	} else {
//		KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED);
//		KSPSetResidualHistory(ksp, solver->rNormGlobal, maxIter*(v+1), PETSC_FALSE);
//		KSPMonitorSet(ksp, rNormGridMonitor, &info, NULL);
//	}
//	KSPSetTolerances(ksp, 1.e-16, PETSC_DEFAULT, PETSC_DEFAULT, v);
//	KSPSetFromOptions(ksp);
//	KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);
//	
//	VecDuplicate(*b, &r);
//	VecGetSubVector(*b, *botIS, &bBot);
//	VecGetSubVector(r, *topIS, &rTop);
//	VecGetSubVector(*u, *topIS, &uTop);
//
//	VecNorm(*b, NORM_2, &bnorm);
//	VecSet(*u, 0); // Note: This should be moved out of this function?
//	
//	Vec	residual;
//
//	MatMult(*A, *u, r);
//	VecAYPX(r, -1.0, *b);
//	VecNorm(r, NORM_2, &chkNorm);
//	norm[0] = chkNorm;
//	
//	iter = 0;
//	double initWallTime = MPI_Wtime();
//	clock_t solverInitT = clock();
//	PetscLogStageRegister("Solver", &stageSolve);
//	PetscLogStagePush(stageSolve);
//	while (iter<maxIter && 100000000*bnorm > chkNorm && chkNorm > (1.e-7)*bnorm) {
//		MatMult(*res, rTop, bBot);
//		MatMult(*pro, *u, rTop);
//		VecAXPY(uTop, 1.0, rTop);
//		KSPSolve(ksp, *b, *u);
//		KSPBuildResidual(ksp, NULL, r, &residual);
//		VecNorm(r, NORM_2, &chkNorm);
//		iter += 1;
//		norm[iter] = chkNorm;
//	}
//	
//	PetscLogStagePop();
//	clock_t solverT = clock();
//	double endWallTime = MPI_Wtime();
//	VecRestoreSubVector(*b, *botIS, &bBot);
//	VecRestoreSubVector(r, *topIS, &rTop);
//	VecRestoreSubVector(*u, *topIS, &uTop);
//	if (solver->moreInfo != 0) {
//		for (int i=0; i<info.grids; i++) {
////			VecView(uGrid[i],PETSC_VIEWER_STDOUT_WORLD);
//			VecRestoreSubVector(info.rInner, gridIS[0][i], info.rGrid+i);
//		}
//		free(info.rGrid);
//		VecDestroy(&(info.rInner));
//	}
//
//	solver->numIter = iter;
//	chkNorm = norm[0];
//	for (int i=0;i<solver->numIter+1;i++) {
//		norm[i] = norm[i]/chkNorm;
//	}
//	if (solver->moreInfo != 0) {
//		chkNorm = solver->rNormGlobal[0];
//		for (int i=0;i<solver->numIter*(v+1);i++) {
//			solver->rNormGlobal[i] = solver->rNormGlobal[i]/chkNorm;
//		}
//		for (int g=0; g<solver->grids; g++) {
//			chkNorm = solver->rNormGrid[g][0];
//			for (int i=0;i<solver->numIter*(v+1);i++) {
//				solver->rNormGrid[g][i] = solver->rNormGrid[g][i]/chkNorm;
//			}
//		}
//	}
//
//	PetscPrintf(PETSC_COMM_WORLD,"---------------------------| level = 0 |------------------------\n");
//	KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
//	PetscPrintf(PETSC_COMM_WORLD,"----------------------------------------------------------------\n");
//	VecDestroy(&r);
//	KSPDestroy(&ksp);
//
//	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = [%d]; Solver cputime:                %lf\n",rank,(double)(solverT-solverInitT)/CLOCKS_PER_SEC);
//	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = [%d]; Solver walltime:               %lf\n",rank,endWallTime-initWallTime);
//	PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
//}

int Solve(Solver *solver){
	// Solves the problem with chosen multigrid cycle
	
	//Assemble(Problem *prob, Mesh *mesh, Levels *levels, Operator *op, Solver *solver);
	int	ierr=0;
	if (solver->cycle == 0) {
		ierr = MultigridVcycle(solver); pCHKERR_RETURN("Multigrid V-cycle solver failed");
	}
//	if (solver->cycle == ICYCLE) MultigridIcycle(solver);
//	if (solver->cycle == ECYCLE) MultigridEcycle(solver);
//	if (solver->cycle == D1CYCLE) MultigridD1cycle(solver);
//	if (solver->cycle == D2CYCLE) MultigridD2cycle(solver);
//	if (solver->cycle == D1PSCYCLE) MultigridD1PScycle(solver);
//	if (solver->cycle == PetscPCMG) MultigridPetscPCMG(solver);
//	if (solver->cycle == FILTER) MultigridFilter(solver);
//	if (solver->cycle == VFILTER) MultigridVFilter(solver);
//	if (solver->cycle == ADDITIVE) MultigridAdditive(solver);
//	if (solver->cycle == ADDITIVEScaled) MultigridAdditiveScaled(solver);
	
	return 0;	
}


