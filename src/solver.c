#include "solver.h"

#define ERROR_MSG(message) {PetscBarrier(PETSC_NULL);(fprintf(stderr,"ERROR: %s:%d: %s\n",__FILE__,__LINE__,(message)));}
//#define ERROR_MSG(message) (fprintf(stderr,"ERROR: %s:%d: %s\n",__FILE__,__LINE__,(message)))
#define ERROR_RETURN(message) {ERROR_MSG(message);return ierr;}
#define CHKERR_PRNT(message) {if(ierr != 0) {ERROR_MSG(message);}}
#define CHKERR_RETURN(message) {if(ierr != 0) {ERROR_RETURN(message);}}
#define PI 3.14159265358979323846

#define pERROR_MSG(message) {PetscBarrier(PETSC_NULL);(PetscPrintf(PETSC_COMM_WORLD,"ERROR: %s:%d: %s\n",__FILE__,__LINE__,(message)));}
//#define pERROR_MSG(message) (PetscPrintf(PETSC_COMM_WORLD,"ERROR: %s:%d: %s\n",__FILE__,__LINE__,(message)))
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

void InitializeSolver(Solver *solver) {
	solver->rnorm	= NULL;
	solver->levels	= NULL;
}

void DestroySolver(Solver *solver) {
	// Free the memory in Solver struct
	
	if (!solver) return;	
	DestroyLevels(solver->levels);
	if (solver->levels) free(solver->levels);
	if (solver->rnorm) free(solver->rnorm);
}

//void CreatePostProcess(PostProcess *pp) {
//	// Allocates memory to PostProcess struct
//		
//	int	procs, rank;
//	
//	MPI_Comm_size(PETSC_COMM_WORLD, &procs);
//	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
//	
//	if (rank == 0) {
//		pp->solData = fopen("uData.dat","w");
//		pp->resData = fopen("rData.dat","w");
//		pp->errData = fopen("eData.dat","w");
//		pp->XgridData = fopen("XgridData.dat","w");
//		pp->YgridData = fopen("YgridData.dat","w");
//	}
//}
//
//void DestroyPostProcess(PostProcess *pp) {
//	// Free the memory in PostProcess struct
//	
//	int	procs, rank;
//	
//	MPI_Comm_size(PETSC_COMM_WORLD, &procs);
//	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
//	
//	if (rank == 0) {
//		fclose(pp->solData);
//		fclose(pp->resData);
//		fclose(pp->errData);
//		fclose(pp->XgridData);
//		fclose(pp->YgridData);
//	}
//}

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

inline double FD2Der2OrderSide(double dx, double del) {
	// Gives 2nd order FD right point coefficient for 2nd derivative
	
	return 2.0/(dx*del);
}

inline double FD2Der2OrderMid(double *dx) {
	// Gives 2nd order FD mid point coefficient for 2nd derivative
	
	return -2.0/(dx[0]*dx[1]);
}

inline void FillDirMatA(int igstart, int istart,
			long int inc0, long int inc1, long int inc2,
			int ln0, int ln1, int ln2,
			double *xcoord, double *dx, Mat *A) {
	// Fill left, mid and right along a direction
	
	// Fill "i-center" for all points
	for (int i=0; i<ln0; i++) {
		double val = FD2Der2OrderMid(dx+i);
		for (int k=0; k<ln2; k++) {
			for (int j=0; j<ln1; j++) {
				int row = igstart + i*inc0 + j*inc1 + k*inc2;
				MatSetValues(*A, 1, &row, 1, &row, &val, ADD_VALUES);
			}
		}
	}
	
	double *del = malloc(ln0*sizeof(double));
	
	for (int i=0; i<ln0; i++) {
		int ic = istart + i;
		del[i] = xcoord[ic+1]-xcoord[ic-1];
	}

	// Fill "i-right" for all points
	for (int i=0; i<ln0-1; i++) {
		double val = FD2Der2OrderSide(dx[i+1], del[i]);
		for (int k=0; k<ln2; k++) {
			for (int j=0; j<ln1; j++) {
				int row = igstart + i*inc0 + j*inc1 + k*inc2;
				int col = row+inc0;
				MatSetValues(*A, 1, &row, 1, &col, &val, ADD_VALUES);
			}
		}
	}
	// Fill "i-left" for all points
	for (int i=1; i<ln0; i++) {
		double val = FD2Der2OrderSide(dx[i-1], del[i]);
		for (int k=0; k<ln2; k++) {
			for (int j=0; j<ln1; j++) {
				int row = igstart + i*inc0 + j*inc1 + k*inc2;
				int col = row-inc0;
				MatSetValues(*A, 1, &row, 1, &col, &val, ADD_VALUES);
			}
		}
	}
	free(del);
}

inline void FillBCDirMatA(long int igstart, long int istart, 
			long int bcStartIndex0, long int bcinc01, long int bcinc02, 
			long int bcStartIndex1, long int bcinc11, long int bcinc12,
			long int inc0, long int inc1, long int inc2,
			int ln0, int ln1, int ln2,
			double *dx, double *coord, Mat *A) {
	// Fill contribution from adjacent blocks at BCs of this block in 
	// a single direction
	
	
	double del = coord[istart+1]-coord[istart-1];
	double val = FD2Der2OrderSide(dx[istart-1], del);
	if (bcStartIndex0 >= 0) {
	for (int k=0; k<ln2; k++) {
		for (int j=0; j<ln1; j++) {
			int row = igstart + j*inc1 + k*inc2;
			int col = bcStartIndex0 + j*bcinc01 + k*bcinc02;
			MatSetValues(*A, 1, &row, 1, &col, &val, ADD_VALUES);
		}
	}
	}
	if (bcStartIndex1 >= 0) {
	del = coord[istart+ln0]-coord[istart+ln0-2];
	val = FD2Der2OrderSide(dx[istart+ln0-1], del);
	for (int k=0; k<ln2; k++) {
		for (int j=0; j<ln1; j++) {
			int row = igstart + (ln0-1)*inc0 + j*inc1 + k*inc2;
			int col = bcStartIndex1 + j*bcinc11 + k*bcinc12;
			MatSetValues(*A, 1, &row, 1, &col, &val, ADD_VALUES);
		}
	}
	}
}

static long int GetFGridIndex(long int i, int cfactor) {
	// Gives Fine grid index corresponding to the coarse grid index
	return (i+1)*cfactor-1;
}

static void GetFGridPoint2D(long int i1, long int j1, 
			long int *i0, long int *j0, 
			int cfactor) {
	// Gives fine grid point corresponding to the coarse grid point in 2D
	
	*i0 = GetFGridIndex(i1, cfactor);
	*j0 = GetFGridIndex(j1, cfactor);
}

static void GetFGridPoint3D(long int i1, long int j1, long int k1,
			long int *i0, long int *j0, long int *k0,
			int cfactor) {
	// Gives fine grid point corresponding to the coarse grid point in 3D
	
	*i0 = GetFGridIndex(i1, cfactor);
	*j0 = GetFGridIndex(j1, cfactor);
	*k0 = GetFGridIndex(k1, cfactor);
}

//int GetIndexIfCorner3D(int p, int q, int r, int *ln, BCindices (*cbcindices)[2][2]) {
//	// Check if (p, q, r) is a corner. If yes, give appropriate index
//	// If not, give negative number
//	
//	if ((p >= 0 && p < ln[0]) || (q >= 0 && q < ln[1]) || (r >= 0 && r < ln[2]))
//		return -1;
//
//	int i = (p<0)? 0: 1;
//	int j = (q<0)? 0: 1;
//	int k = (r<0)? 0: 1;
//
//	return cbcindices[i][j][k].bcGStartIndex;
//}

int GetColForBCFace3D(int dim, int *id, int i, int j, int k, BCindices (*bcindices)[2]) {
	// Compute Col for i,j,k in BC face
	
	BCindices *bc = &(bcindices[dim][id[dim]]);
	int col = bc->bcGStartIndex;
	long int *inc = bc->bcInc;

	if (dim == 0)
		col += j*inc[1]+k*inc[2];
	else if (dim == 1)
		col += i*inc[0]+k*inc[2];
	else if (dim == 2)
		col += j*inc[1]+i*inc[0];

	return col;
}

int GetColForBCEdge3D(int dim, int *id, int i, int j, int k, BCindices (*ebcindices)[2][2]) {
	// Compute Col for i,j,k in BC edge
	
	int col = -1;	
	if (dim == 0) {
		BCindices *bc = &(ebcindices[0][id[1]][id[2]]);
		col = bc->bcGStartIndex + i*bc->bcInc[0]; 
	} else if (dim == 1) {
		BCindices *bc = &(ebcindices[1][id[2]][id[0]]);
		col = bc->bcGStartIndex + j*bc->bcInc[1]; 
	} else if (dim == 2) {
		BCindices *bc = &(ebcindices[2][id[0]][id[1]]);
		col = bc->bcGStartIndex + k*bc->bcInc[2];
	}
	return col;
}

void GetBCType3D(int p, int q, int r, int *ln0, int *id) {
	// Compute BC type and give it back in "id"
	
	id[0] = (p<0)? 0: ((p>ln0[0]-1)? 1: 6);
	id[1] = (q<0)? 0: ((q>ln0[1]-1)? 1: 9);
	id[2] = (r<0)? 0: ((r>ln0[2]-1)? 1: 12);
}

void FillMatResBC3D(Grids *grids, int lg0, Level *level0, int lg1, Level *level1, 
			double *weights, Mat *res) {
	// Fills Mat res for block BC points in 3D
	
	int	*blockID = grids->topo->blockID;

	int	g0 = level0->gridId[lg0];
	int	g1 = level1->gridId[lg1];

	Grid 	*grid0 = grids->grid+g0;
	Grid 	*grid1 = grids->grid+g1;

	int	*ln0 = grid0->ln;
	int	*ln1 = grid1->ln;

	long int *inc0 = level0->inc[lg0];
	long int *inc1 = level1->inc[lg1];

	long int gstart0 = level0->granges[lg0][0];
	long int gstart1 = level1->granges[lg1][0];

	long int istart0 = grid0->range[0][blockID[0]];
	long int jstart0 = grid0->range[1][blockID[1]];
	long int kstart0 = grid0->range[2][blockID[2]];

	long int istart1 = grid1->range[0][blockID[0]];
	long int jstart1 = grid1->range[1][blockID[1]];
	long int kstart1 = grid1->range[2][blockID[2]];

	int	cols[27];
	
	BCindices	(*bcindices)[2] = level0->bcindices[lg0];
	BCindices	(*ebcindices)[2][2] = level0->ebcindices[lg0];
	BCindices	(*cbcindices)[2][2] = level0->cbcindices[lg0];
	for (int k=0; k<ln1[2]; k++) {
	for (int j=0; j<ln1[1]; j++) {
	for (int i=0; i<ln1[0]; i++) {
		if (i != 0 && i != ln1[0]-1 
			&& j != 0 && j != ln1[1]-1
			&& k != 0 && k != ln1[2]-1) continue;
		int row = gstart1 + i*inc1[0] + j*inc1[1] + k*inc1[2];
		long int i0, j0, k0;
		GetFGridPoint3D(istart1+i-1, jstart1+j-1, kstart1+k-1, &i0, &j0, &k0, 2);
		i0 = i0-istart0;
		j0 = j0-jstart0;
		k0 = k0-kstart0;
		for (int r=0; r<3; r++) {
		for (int q=0; q<3; q++) {
		for (int p=0; p<3; p++) {
			int id[3], total;
			GetBCType3D(i0+p, j0+q, k0+r, ln0, id);
			total = id[0]+id[1]+id[2];
			if (total == 27) { // Interior
				cols[r*9+q*3+p] = gstart0
						+(i0+p)*inc0[0]
						+(j0+q)*inc0[1]
						+(k0+r)*inc0[2];
			} else if (total > 14) { // Faces
				int dim = (22-total)/3;	
				cols[r*9+q*3+p] = 
					GetColForBCFace3D(dim, id, i0+p, j0+q, k0+r, bcindices); 
			} else if (total > 5) { // Edges
				int dim = (total-6)/3;
				cols[r*9+q*3+p] = 
					GetColForBCEdge3D(dim, id, i0+p, j0+q, k0+r, ebcindices); 
			} else { // Corners
				cols[r*9+q*3+p] = cbcindices[id[0]][id[1]][id[2]].bcGStartIndex;
			}

//			col = GetIndexIfCorner3D(i0+p, j0+q, k0+r, ln0, cbcindices);
//			if (i0+p < 0 && j0+q < 0 && k0+r < 0)	// [Corners
//				cols[r*9+q*3+p] = cbcindices[0][0][0].bcGStartIndex; 
//			else if (i0+p < 0 && j0+q < 0 && k0+r > ln0[2]-1) 
//				cols[r*9+q*3+p] = cbcindices[0][0][1].bcGStartIndex;
//			else if (i0+p < 0 && j0+q > ln0[1]-1 && k0+r < 0) 
//				cols[r*9+q*3+p] = cbcindices[0][1][0].bcGStartIndex;
//			else if (i0+p < 0 && j0+q > ln0[1]-1 && k0+r > ln0[2]-1) 
//				cols[r*9+q*3+p] = cbcindices[0][1][1].bcGStartIndex;
//			else if (i0+p > ln0[0]-1 && j0+q < 0 && k0+r < 0) 
//				cols[r*9+q*3+p] = cbcindices[1][0][0].bcGStartIndex;
//			else if (i0+p > ln0[0]-1 && j0+q < 0 && k0+r > ln0[2]-1) 
//				cols[r*9+q*3+p] = cbcindices[1][0][1].bcGStartIndex;
//			else if (i0+p > ln0[0]-1 && j0+q > ln0[1]-1 && k0+r < 0) 
//				cols[r*9+q*3+p] = cbcindices[1][1][0].bcGStartIndex;
//			else if (i0+p > ln0[0]-1 && j0+q > ln0[1]-1 && k0+r > ln0[2]-1) 
//				cols[r*9+q*3+p] = cbcindices[1][1][1].bcGStartIndex; // Corners]
//			else if (j0+q < 0 && k0+r < 0)	// [Edges i
//				cols[r*9+q*3+p] = ebcindices[0][0][0].bcGStartIndex
//						+ (i0+p)*ebcindices[0][0][0].bcInc[0]; 
//			else if (j0+q < 0 && k0+r > ln0[2]-1)
//				cols[r*9+q*3+p] = ebcindices[0][0][1].bcGStartIndex
//						+ (i0+p)*ebcindices[0][0][1].bcInc[0]; 
//			else if (j0+q > ln0[1]-1 && k0+r < 0)
//				cols[r*9+q*3+p] = ebcindices[0][1][0].bcGStartIndex
//						+ (i0+p)*ebcindices[0][1][0].bcInc[0]; 
//			else if (j0+q > ln0[1]-1 && k0+r > ln0[2]-1)
//				cols[r*9+q*3+p] = ebcindices[0][1][1].bcGStartIndex
//						+ (i0+p)*ebcindices[0][1][1].bcInc[0]; // Edges i]
//			else if (i0+p < 0 && k0+r < 0)	// [Edges j
//				cols[r*9+q*3+p] = ebcindices[1][0][0].bcGStartIndex
//						+ (j0+q)*ebcindices[1][0][0].bcInc[1]; 
//			else if (i0+p < 0 && k0+r > ln0[2]-1)
//				cols[r*9+q*3+p] = ebcindices[1][1][0].bcGStartIndex
//						+ (j0+q)*ebcindices[1][1][0].bcInc[1]; 
//			else if (i0+p > ln0[0]-1 && k0+r < 0)
//				cols[r*9+q*3+p] = ebcindices[1][0][1].bcGStartIndex
//						+ (j0+q)*ebcindices[1][0][1].bcInc[1]; 
//			else if (i0+p > ln0[0]-1 && k0+r > ln0[2]-1)
//				cols[r*9+q*3+p] = ebcindices[1][1][1].bcGStartIndex
//						+ (j0+q)*ebcindices[1][1][1].bcInc[1]; // Edges j]
//			else if (i0+p < 0 && j0+q < 0)	// [Edges k
//				cols[r*9+q*3+p] = ebcindices[2][0][0].bcGStartIndex
//						+ (k0+r)*ebcindices[2][0][0].bcInc[2]; 
//			else if (i0+p < 0 && j0+q > ln0[1]-1)
//				cols[r*9+q*3+p] = ebcindices[2][0][1].bcGStartIndex
//						+ (j0+q)*ebcindices[2][0][1].bcInc[2]; 
//			else if (i0+p > ln0[0]-1 && k0+r < 0)
//				cols[r*9+q*3+p] = ebcindices[1][0][1].bcGStartIndex
//						+ (j0+q)*ebcindices[1][0][1].bcInc[1]; 
//			else if (i0+p > ln0[0]-1 && k0+r > ln0[2]-1)
//				cols[r*9+q*3+p] = ebcindices[1][1][1].bcGStartIndex
//						+ (j0+q)*ebcindices[1][1][1].bcInc[1]; // Edges k]
//			else if (i0+p < 0)
//				cols[r*9+q*3+p] = bcindices[0][0].bcGStartIndex
//					+(j0+q)*bcindices[0][0].bcInc[1];
//			else if (i0+p > ln0[0]-1)
//				cols[r*9+q*3+p] = bcindices[0][1].bcGStartIndex
//					+(j0+q)*bcindices[0][1].bcInc[1];
//			else if (j0+q < 0)
//				cols[r*9+q*3+p] = bcindices[1][0].bcGStartIndex
//					+(i0+p)*bcindices[1][0].bcInc[0];
//			else if (j0+q > ln0[1]-1)
//				cols[r*9+q*3+p] = bcindices[1][1].bcGStartIndex
//					+(i0+p)*bcindices[1][1].bcInc[0];
//			else
//				cols[r*9+q*3+p] = gstart0
//					+(i0+p)*inc0[0]
//					+(j0+q)*inc0[1]
//					+(k0+r)*inc0[2];
		}
		}
		}
		MatSetValues(*res, 1, &row, 27, cols, weights, INSERT_VALUES);
	}
	}
	}
}

void FillMatResBC2D(Grids *grids, int lg0, Level *level0, int lg1, Level *level1, 
			double *weights, Mat *res) {
	// Fills Mat res for block BC points in 2D
	
	int	*blockID = grids->topo->blockID;

	int	g0 = level0->gridId[lg0];
	int	g1 = level1->gridId[lg1];

	Grid 	*grid0 = grids->grid+g0;
	Grid 	*grid1 = grids->grid+g1;

	int	*ln0 = grid0->ln;
	int	*ln1 = grid1->ln;

	long int *inc0 = level0->inc[lg0];
	long int *inc1 = level1->inc[lg1];

	long int gstart0 = level0->granges[lg0][0];
	long int gstart1 = level1->granges[lg1][0];

	long int istart0 = grid0->range[0][blockID[0]];
	long int jstart0 = grid0->range[1][blockID[1]];

	long int istart1 = grid1->range[0][blockID[0]];
	long int jstart1 = grid1->range[1][blockID[1]];

	int	cols[9];
	
	BCindices	(*bcindices)[2] = level0->bcindices[lg0];
	BCindices	(*cbcindices)[2][2] = level0->cbcindices[lg0];
	for (int i=0; i<ln1[0]; i++) {
		for (int j=0; j<ln1[1]; j++) {
			if (i != 0 && i != ln1[0]-1 && j != 0 && j != ln1[1]-1) continue;
			int row = gstart1 + i*inc1[0] + j*inc1[1];
			long int i0, j0;
			GetFGridPoint2D(istart1+i-1, jstart1+j-1, &i0, &j0, 2);
			i0 = i0-istart0;
			j0 = j0-jstart0;
			for (int q=0; q<3; q++) {
				for (int p=0; p<3; p++) {
					if (i0+p < 0 && j0+q < 0) 
						cols[q*3+p] = cbcindices[0][0][0].bcGStartIndex;
					else if (i0+p < 0 && j0+q > ln0[1]-1) 
						cols[q*3+p] = cbcindices[0][1][0].bcGStartIndex;
					else if (i0+p > ln0[0]-1 && j0+q < 0) 
						cols[q*3+p] = cbcindices[1][0][0].bcGStartIndex;
					else if (i0+p > ln0[0]-1 && j0+q > ln0[1]-1) 
						cols[q*3+p] = cbcindices[1][1][0].bcGStartIndex;
					else if (i0+p < 0)
						cols[q*3+p] = bcindices[0][0].bcGStartIndex
							+(j0+q)*bcindices[0][0].bcInc[1];
					else if (i0+p > ln0[0]-1)
						cols[q*3+p] = bcindices[0][1].bcGStartIndex
							+(j0+q)*bcindices[0][1].bcInc[1];
					else if (j0+q < 0)
						cols[q*3+p] = bcindices[1][0].bcGStartIndex
							+(i0+p)*bcindices[1][0].bcInc[0];
					else if (j0+q > ln0[1]-1)
						cols[q*3+p] = bcindices[1][1].bcGStartIndex
							+(i0+p)*bcindices[1][1].bcInc[0];
					else
						cols[q*3+p] = gstart0
							+(i0+p)*inc0[0]+(j0+q)*inc0[1];
				}
			}
			MatSetValues(*res, 1, &row, 9, cols, weights, INSERT_VALUES);
		}
	}
}

void FillMatResInterior3D(Grids *grids, int lg0, Level *level0, int lg1, Level *level1, 
			double *weights, Mat *res) {
	// Fills Mat res for interior coarse points
	
	int	*blockID = grids->topo->blockID;

	int	g0 = level0->gridId[lg0];
	int	g1 = level1->gridId[lg1];

	Grid 	*grid0 = grids->grid+g0;
	Grid 	*grid1 = grids->grid+g1;

	int	*ln1 = grid1->ln;

	long int *inc0 = level0->inc[lg0];
	long int *inc1 = level1->inc[lg1];

	long int gstart0 = level0->granges[lg0][0];
	long int gstart1 = level1->granges[lg1][0];

	long int istart0 = grid0->range[0][blockID[0]];
	long int jstart0 = grid0->range[1][blockID[1]];
	long int kstart0 = grid0->range[2][blockID[2]];

	long int istart1 = grid1->range[0][blockID[0]];
	long int jstart1 = grid1->range[1][blockID[1]];
	long int kstart1 = grid1->range[2][blockID[2]];

	int	cols[27];
	for (int k=1; k<ln1[2]-1; k++) {
		for (int j=1; j<ln1[1]-1; j++) {
			for (int i=1; i<ln1[0]-1; i++) {
				int row = gstart1 + i*inc1[0] + j*inc1[1] + k*inc1[2];
				long int i0, j0, k0;
				GetFGridPoint3D(istart1+i-1, jstart1+j-1, kstart1+k-1, &i0, &j0, &k0, 2);
				i0 = i0-istart0;
				j0 = j0-jstart0;
				k0 = k0-kstart0;
				for (int r=0; r<3; r++) {
					for (int q=0; q<3; q++) {
						for (int p=0; p<3; p++) {
							cols[r*9+q*3+p] = gstart0+
								(i0+p)*inc0[0]+
								(j0+q)*inc0[1]+
								(k0+r)*inc0[2];
						}
					}
				}
				MatSetValues(*res, 1, &row, 27, cols, weights, INSERT_VALUES);
			}
		}
	}
}

void FillMatResInterior2D(Grids *grids, int lg0, Level *level0, int lg1, Level *level1, 
			double *weights, Mat *res) {
	// Fills Mat res for interior coarse points
	
	int	*blockID = grids->topo->blockID;

	int	g0 = level0->gridId[lg0];
	int	g1 = level1->gridId[lg1];

	Grid 	*grid0 = grids->grid+g0;
	Grid 	*grid1 = grids->grid+g1;

	int	*ln1 = grid1->ln;

	long int *inc0 = level0->inc[lg0];
	long int *inc1 = level1->inc[lg1];

	long int gstart0 = level0->granges[lg0][0];
	long int gstart1 = level1->granges[lg1][0];

	long int istart0 = grid0->range[0][blockID[0]];
	long int jstart0 = grid0->range[1][blockID[1]];

	long int istart1 = grid1->range[0][blockID[0]];
	long int jstart1 = grid1->range[1][blockID[1]];

	int	cols[9];
	for (int i=1; i<ln1[0]-1; i++) {
		for (int j=1; j<ln1[1]-1; j++) {
			int row = gstart1 + i*inc1[0] + j*inc1[1];
			long int i0, j0;
			GetFGridPoint2D(istart1+i-1, jstart1+j-1, &i0, &j0, 2);
			i0 = i0-istart0;
			j0 = j0-jstart0;
			for (int q=0; q<3; q++) {
				for (int p=0; p<3; p++) {
					cols[q*3+p] = gstart0+
						(i0+p)*inc0[0]+
						(j0+q)*inc0[1];
				}
			}
			MatSetValues(*res, 1, &row, 9, cols, weights, INSERT_VALUES);
		}
	}
}

void FillMatRes(Grids *grids, int lg0, Level *level0, int lg1, Level *level1, Mat *res){
	// Fills Mat res with restriction weights from grid lg0 in leve0 to grid lg1 in level1
	
	int	dimension = grids->topo->dimension;
	
	// Assuming 3x3x3 restriction stencil
	double	opx[3] = {0.25, 0.5, 0.25};
	double	opy[3], opz[3];
	for (int i=0; i<3; i++) {
		opy[i] = opx[i];
		opz[i] = opx[i];
	}

	if (dimension == 2) {
		double	weights[9];
		for (int j=0; j<3; j++) {
			for (int i=0; i<3; i++) {
				weights[j*3+i] = opx[i]*opy[j];
			}
		}
		FillMatResInterior2D(grids, lg0, level0, lg1, level1, weights, res);
		FillMatResBC2D(grids, lg0, level0, lg1, level1, weights, res);
	} else if (dimension == 3) {
		double	weights[27];
		for (int k=0; k<3; k++) {
			for (int j=0; j<3; j++) {
				for (int i=0; i<3; i++) {
					weights[k*9+j*3+i] = opx[i]*opy[j]*opz[k];
				}
			}
		}
		FillMatResInterior3D(grids, lg0, level0, lg1, level1, weights, res);
		FillMatResBC3D(grids, lg0, level0, lg1, level1, weights, res);
	}
}

void FillMatA(Grids *grids, Level *level, Mat *A) {
	// Fills Mat A with the Jacobian or Discretized PDE coefficients of all grids this level possesses

	long int	*ranges = level->ranges;
	int		ngrids = level->ngrids;
	int		*gridId = level->gridId;
	Grid		*grid = grids->grid;
	int		*l = grids->topo->blockID;
	int		dimension = grids->topo->dimension;
	
	if (dimension == 2) {
		for (int lg=0; lg<ngrids; lg++) {
			int igstart = ranges[lg];
			long int *inc = level->inc[lg];
			int g = gridId[lg];
			int *ln = grid[g].ln;
			double *dx = grid[g].dx[0];
			double *dy = grid[g].dx[1];
			double *xcoord = grid[g].coord[0];
			double *ycoord = grid[g].coord[1];
			
			if (ln[0] == 0 || ln[1] == 0) continue;
			int istart = grid[g].range[0][l[0]];
			FillDirMatA(igstart, istart, inc[0], inc[1], 0, ln[0], ln[1], 1,
					xcoord, dx+(istart-1), A); // Fill along i^th direction
			
			int jstart = grid[g].range[1][l[1]];
			FillDirMatA(igstart, jstart, inc[1], inc[0], 0, ln[1], ln[0], 1,
					ycoord, dy+(jstart-1), A); // Fill along j^th direction
			long int bcStartIndex0	= level->bcindices[lg][0][0].bcStartIndex;
			long int bcStartIndex1 = level->bcindices[lg][0][1].bcStartIndex;
			long int *bcinc0	= level->bcindices[lg][0][0].bcInc;
			long int *bcinc1	= level->bcindices[lg][0][1].bcInc;
			FillBCDirMatA(igstart, istart, bcStartIndex0, bcinc0[1], 0, 
				bcStartIndex1, bcinc1[1], 0, 
				inc[0], inc[1], 0, ln[0], ln[1], 1, dx, xcoord, A);
			
			bcStartIndex0 = level->bcindices[lg][1][0].bcStartIndex;
			bcStartIndex1 = level->bcindices[lg][1][1].bcStartIndex;
			bcinc0	= level->bcindices[lg][1][0].bcInc;
			bcinc1	= level->bcindices[lg][1][1].bcInc;
			FillBCDirMatA(igstart, jstart, bcStartIndex0,  bcinc0[0], 0,
				bcStartIndex1, bcinc1[0], 0,
				inc[1], inc[0], 0, ln[1], ln[0], 1, dy, ycoord, A);
		}
	} else if (dimension == 3) {
		for (int lg=0; lg<ngrids; lg++) {
			int igstart = ranges[lg];
			long int *inc = level->inc[lg];
			int g = gridId[lg];
			int *ln = grid[g].ln;
			double *dx = grid[g].dx[0];
			double *dy = grid[g].dx[1];
			double *dz = grid[g].dx[2];
			double *xcoord = grid[g].coord[0];
			double *ycoord = grid[g].coord[1];
			double *zcoord = grid[g].coord[2];
			
			if (ln[0] == 0 || ln[1] == 0 || ln[2] == 0) continue;
			// Contributions from within block
			int istart = grid[g].range[0][l[0]];
			FillDirMatA(igstart, istart, inc[0], inc[1], inc[2], 
					ln[0], ln[1], ln[2], xcoord, 
					dx+(istart-1), A); // Fill along i^th direction
			
			int jstart = grid[g].range[1][l[1]];
			FillDirMatA(igstart, jstart, inc[1], inc[0], inc[2], 
					ln[1], ln[0], ln[2], ycoord, 
					dy+(jstart-1), A); // Fill along j^th direction
			
			int kstart = grid[g].range[2][l[2]];
			FillDirMatA(igstart, kstart, inc[2], inc[0], inc[1], 
					ln[2], ln[0], ln[1], zcoord, 
					dz+(kstart-1), A); // Fill along z^th direction
			
			// Contributions from adjacent blocks
			long int bcStartIndex0	= level->bcindices[lg][0][0].bcStartIndex;
			long int bcStartIndex1	= level->bcindices[lg][0][1].bcStartIndex;
			long int *bcinc0	= level->bcindices[lg][0][0].bcInc;
			long int *bcinc1	= level->bcindices[lg][0][1].bcInc;
			FillBCDirMatA(igstart, istart, bcStartIndex0, bcinc0[1], bcinc0[2],
				bcStartIndex1,  bcinc1[1], bcinc1[2],
				inc[0], inc[1], inc[2], ln[0], ln[1], ln[2], dx, xcoord, A);
			
			bcStartIndex0 = level->bcindices[lg][1][0].bcStartIndex;
			bcStartIndex1 = level->bcindices[lg][1][1].bcStartIndex;
			bcinc0	= level->bcindices[lg][1][0].bcInc;
			bcinc1	= level->bcindices[lg][1][1].bcInc;
			FillBCDirMatA(igstart, jstart, bcStartIndex0, bcinc0[0], bcinc0[2], 
				bcStartIndex1, bcinc1[0], bcinc1[2],
				inc[1], inc[0], inc[2], ln[1], ln[0], ln[2], dy, ycoord, A);
			
			bcStartIndex0 = level->bcindices[lg][2][0].bcStartIndex;
			bcStartIndex1 = level->bcindices[lg][2][1].bcStartIndex;
			bcinc0	= level->bcindices[lg][2][0].bcInc;
			bcinc1	= level->bcindices[lg][2][1].bcInc;
			FillBCDirMatA(igstart, kstart, bcStartIndex0, bcinc0[0], bcinc0[1], 
				bcStartIndex1, bcinc1[0], bcinc1[1],
				inc[2], inc[0], inc[1], ln[2], ln[0], ln[1], dz, zcoord, A);
		}
	}
}

void AssembleMatRes(Grids *grids, int lg0, Level *level0, int lg1, Level *level1, Mat *res) {
	// Build restriction matrix "res" from "lg0" in level0 to "lg1" in level1
	
	int	g0 = level0->gridId[lg0];	
	int	g1 = level1->gridId[lg1];
	int	dimension = grids->topo->dimension;
	if (g1-g0 != 1) return; // Only successive grids are allowed!

	int	dnz = 3;
	for (int i=1; i<dimension; i++) dnz *= 3;
	int	onz = dnz-1;
	Grid	*grid = grids->grid;
	int	nrows = grid[g1].tln;
	int	ncols = grid[g0].tln;
	MatCreateAIJ(PETSC_COMM_WORLD, nrows, ncols, PETSC_DETERMINE, PETSC_DETERMINE, dnz, PETSC_NULL, onz, PETSC_NULL, res);

	FillMatRes(grids, lg0, level0, lg1, level1, res);
	
	MatAssemblyBegin(*res, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(*res, MAT_FINAL_ASSEMBLY);
}

void AssembleLevelMatA(Grids *grids, Level *level, Mat *A) {
	// Build matrix "A" for a given level
	
	int	nnz	= 2*(grids->topo->dimension)+1; // Number of non-zeros per row
	int	ngrids = level->ngrids;
	int	size = (int) (level->ranges[ngrids] - level->ranges[0]);
	MatCreateAIJ(PETSC_COMM_WORLD, size, size, PETSC_DETERMINE, PETSC_DETERMINE, nnz, PETSC_NULL, nnz, PETSC_NULL, A);

	FillMatA(grids, level, A);
	
	MatAssemblyBegin(*A, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(*A, MAT_FINAL_ASSEMBLY);
}

double RHSfunc2D(double x, double y) {
	// Gives the f(x,y)
	
	return -2*PI*PI*sin(PI*x)*sin(PI*y);
}

double RHSfunc3D(double x, double y, double z) {
	// Gives the f(x,y)
	
	return -3*PI*PI*sin(PI*x)*sin(PI*y)*sin(PI*z);
}

double Sol2D(double x, double y) {
	// Gives the f(x,y)
	
	return sin(PI*x)*sin(PI*y);
}

double Sol3D(double x, double y, double z) {
	// Gives the f(x,y)
	
	return sin(PI*x)*sin(PI*y)*sin(PI*z);
}

void ApplyBCLevelVecb3D(Grid *grid, Level *level,  Vec *b) {
	
	double	*xcoord = grid->coord[0];
	double	*ycoord = grid->coord[1];
	double	*zcoord = grid->coord[2];
	double	*dx = grid->dx[0];
	double	*dy = grid->dx[1];
	double	*dz = grid->dx[2];
	int	*blockID = grid->topo->blockID;
	long int istart = grid->range[0][blockID[0]];
	long int jstart = grid->range[1][blockID[1]];
	long int kstart = grid->range[2][blockID[2]];
	int	*ln = grid->ln;

	long int gstart = level->ranges[0];
	long int *inc = level->inc[0];

	// Along i^th-direction 
	double *val = malloc(ln[2]*ln[1]*sizeof(double));
	int *row = malloc(ln[2]*ln[1]*sizeof(int));
	
	long int bcrank0 = level->bcindices[0][0][0].rank;
	long int bcrank1 = level->bcindices[0][0][1].rank;
	double del = xcoord[istart+1]-xcoord[istart-1];
	double coeff = FD2Der2OrderSide(dx[istart-1], del);
	int count = 0;
	if (bcrank0 < 0) {
	for (int k=0; k<ln[2]; k++) {
		for (int j=0; j<ln[1]; j++) {
			row[count] = gstart + j*inc[1] + k*inc[2];
			val[count] = -1*coeff*Sol3D(xcoord[istart-1], ycoord[jstart+j], zcoord[kstart+k]);
			count++;
		}
	}
	VecSetValues(*b, ln[2]*ln[1], row, val, ADD_VALUES);
	}
	
	if (bcrank1 < 0) {
	del = xcoord[istart+ln[0]]-xcoord[istart+ln[0]-2];
	coeff = FD2Der2OrderSide(dx[istart+ln[0]-1], del);
	long int igstart = gstart + (ln[0]-1)*inc[0];
	count = 0;
	for (int k=0; k<ln[2]; k++) {
		for (int j=0; j<ln[1]; j++) {
			row[count] = igstart + j*inc[1] + k*inc[2];
			val[count] = -1*coeff*Sol3D(xcoord[istart+ln[0]], ycoord[jstart+j], zcoord[kstart+k]);
			count++;
		}
	}
	VecSetValues(*b, ln[2]*ln[1], row, val, ADD_VALUES);
	}
	
	free(row);	
	free(val);	
	
	// Along j^th-direction 
	val = malloc(ln[2]*ln[0]*sizeof(double));
	row = malloc(ln[2]*ln[0]*sizeof(int));
	
	bcrank0 = level->bcindices[0][1][0].rank;
	bcrank1 = level->bcindices[0][1][1].rank;
	del = ycoord[jstart+1]-ycoord[jstart-1];
	coeff = FD2Der2OrderSide(dy[jstart-1], del);
	count = 0;
	if (bcrank0 < 0) {
	for (int k=0; k<ln[2]; k++) {
		for (int i=0; i<ln[0]; i++) {
			row[count] = gstart + i*inc[0] + k*inc[2];
			val[count] = -1*coeff*Sol3D(xcoord[istart+i], ycoord[jstart-1], zcoord[kstart+k]);
			count++;
		}
	}
	VecSetValues(*b, ln[2]*ln[0], row, val, ADD_VALUES);
	}
	
	if (bcrank1 < 0) {
	del = ycoord[jstart+ln[1]]-ycoord[jstart+ln[1]-2];
	coeff = FD2Der2OrderSide(dy[jstart+ln[1]-1], del);
	long int jgstart = gstart + (ln[1]-1)*inc[1];
	count = 0;
	for (int k=0; k<ln[2]; k++) {
		for (int i=0; i<ln[0]; i++) {
			row[count] = jgstart + i*inc[0] + k*inc[2];
			val[count] = -1*coeff*Sol3D(xcoord[istart+i], ycoord[jstart+ln[1]], zcoord[kstart+k]);
			count++;
		}
	}
	VecSetValues(*b, ln[2]*ln[0], row, val, ADD_VALUES);
	}
	
	free(row);	
	free(val);	
	
	// Along k^th-direction 
	val = malloc(ln[0]*ln[1]*sizeof(double));
	row = malloc(ln[0]*ln[1]*sizeof(int));
	
	bcrank0 = level->bcindices[0][2][0].rank;
	bcrank1 = level->bcindices[0][2][1].rank;
	del = zcoord[kstart+1]-zcoord[kstart-1];
	coeff = FD2Der2OrderSide(dz[kstart-1], del);
	count = 0;
	if (bcrank0 < 0) {
	for (int j=0; j<ln[1]; j++) {
		for (int i=0; i<ln[0]; i++) {
			row[count] = gstart + i*inc[0] + j*inc[1];
			val[count] = -1*coeff*Sol3D(xcoord[istart+i], ycoord[jstart+j], zcoord[kstart-1]);
			count++;
		}
	}
	VecSetValues(*b, ln[0]*ln[1], row, val, ADD_VALUES);
	}
	
	if (bcrank1 < 0) {
	del = zcoord[kstart+ln[2]]-zcoord[kstart+ln[2]-2];
	coeff = FD2Der2OrderSide(dz[kstart+ln[2]-1], del);
	long int kgstart = gstart + (ln[2]-1)*inc[2];
	count = 0;
	for (int j=0; j<ln[1]; j++) {
		for (int i=0; i<ln[0]; i++) {
			row[count] = kgstart + i*inc[0] + j*inc[1];
			val[count] = -1*coeff*Sol3D(xcoord[istart+i], ycoord[jstart+j], zcoord[kstart+ln[2]]);
			count++;
		}
	}
	VecSetValues(*b, ln[0]*ln[1], row, val, ADD_VALUES);
	}
	
	free(row);	
	free(val);	
}

void ApplyBCLevelVecb2D(Grid *grid, Level *level,  Vec *b) {
	
	double	*xcoord = grid->coord[0];
	double	*ycoord = grid->coord[1];
	double	*dx = grid->dx[0];
	double	*dy = grid->dx[1];
	int	*blockID = grid->topo->blockID;
	long int istart = grid->range[0][blockID[0]];
	long int jstart = grid->range[1][blockID[1]];
	int	*ln = grid->ln;

	long int gstart = level->ranges[0];
	long int *inc = level->inc[0];

	// Along i^th-direction 
	double *val = malloc(ln[1]*sizeof(double));
	int *row = malloc(ln[1]*sizeof(int));
	
	long int bcrank0 = level->bcindices[0][0][0].rank;
	long int bcrank1 = level->bcindices[0][0][1].rank;
	double del = xcoord[istart+1]-xcoord[istart-1];
	double coeff = FD2Der2OrderSide(dx[istart-1], del);
	int count = 0;
	if (bcrank0 < 0) {
	for (int j=0; j<ln[1]; j++) {
		row[count] = gstart + j*inc[1];
		val[count] = -1*coeff*Sol2D(xcoord[istart-1], ycoord[jstart+j]);
		count++;
	}
	VecSetValues(*b, ln[1], row, val, ADD_VALUES);
	}
	
	if (bcrank1 < 0) {
	del = xcoord[istart+ln[0]]-xcoord[istart+ln[0]-2];
	coeff = FD2Der2OrderSide(dx[istart+ln[0]-1], del);
	long int igstart = gstart + (ln[0]-1)*inc[0];
	count = 0;
	for (int j=0; j<ln[1]; j++) {
		row[count] = igstart + j*inc[1];
		val[count] = -1*coeff*Sol2D(xcoord[istart+ln[0]], ycoord[jstart+j]);
		count++;
	}
	VecSetValues(*b, ln[1], row, val, ADD_VALUES);
	}
	
	free(row);	
	free(val);	
	
	// Along j^th-direction 
	val = malloc(ln[0]*sizeof(double));
	row = malloc(ln[0]*sizeof(int));
	
	bcrank0 = level->bcindices[0][1][0].rank;
	bcrank1 = level->bcindices[0][1][1].rank;
	del = ycoord[jstart+1]-ycoord[jstart-1];
	coeff = FD2Der2OrderSide(dy[jstart-1], del);
	count = 0;
	if (bcrank0 < 0) {
	for (int i=0; i<ln[0]; i++) {
		row[count] = gstart + i*inc[0];
		val[count] = -1*coeff*Sol2D(xcoord[istart+i], ycoord[jstart-1]);
		count++;
	}
	VecSetValues(*b, ln[0], row, val, ADD_VALUES);
	}
	
	if (bcrank1 < 0) {
	del = ycoord[jstart+ln[1]]-ycoord[jstart+ln[1]-2];
	coeff = FD2Der2OrderSide(dy[jstart+ln[1]-1], del);
	long int jgstart = gstart + (ln[1]-1)*inc[1];
	count = 0;
	for (int i=0; i<ln[0]; i++) {
		row[count] = jgstart + i*inc[0];
		val[count] = -1*coeff*Sol2D(xcoord[istart+i], ycoord[jstart+ln[1]]);
		count++;
	}
	VecSetValues(*b, ln[0], row, val, ADD_VALUES);
	}
	
	free(row);	
	free(val);	
}

void ApplyBCLevelVecb(Grids *grids, Level *level, Vec *b) {
	// Apply BC to RHS for first grid in the given level
	
	Grid	*grid = grids->grid;
	int	g = level->gridId[0];
	int	tln = grid[g].tln;
	int	dimension = grids->topo->dimension;
	
	if (tln == 0) return;

	if (dimension == 2) {
		ApplyBCLevelVecb2D(grid, level,  b);
	} else if (dimension == 3) {
		ApplyBCLevelVecb3D(grid, level,  b);
	}
}

void FillLevelVecb(int lg, Grid *grid, Level *level, Vec *b) {
	// Builds vector "b" for a given grid "lg" in a given level
	
	int *blockId = grid->topo->blockID;
	int dimension = grid->topo->dimension;

	int g	= level->gridId[lg];
	long int *inc = level->inc[lg];
	int igstart = level->ranges[lg];
	
	int *ln	= grid[g].ln;
	int tln = grid[g].tln;
	
	double *val = malloc(tln*sizeof(double));
	int *row = malloc(tln*sizeof(int));
	
	VecSet(*b, 0.0);
	if (dimension == 3) {
		int istart = grid[g].range[0][blockId[0]];
		int jstart = grid[g].range[1][blockId[1]];
		int kstart = grid[g].range[2][blockId[2]];

		double *xcoord = grid[g].coord[0];
		double *ycoord = grid[g].coord[1];
		double *zcoord = grid[g].coord[2];
		
		int count = 0;
		for (int k=0; k<ln[2]; k++) {
			for (int j=0; j<ln[1]; j++) {
				for (int i=0; i<ln[0]; i++) {
					row[count] = igstart + i*inc[0] + j*inc[1] + k*inc[2];
					val[count] = RHSfunc3D(xcoord[istart+i], ycoord[jstart+j], zcoord[kstart+k]);
					count++;
				}
			}
		}
	} else if (dimension == 2) {
		int istart = grid[g].range[0][blockId[0]];
		int jstart = grid[g].range[1][blockId[1]];

		double *xcoord = grid[g].coord[0];
		double *ycoord = grid[g].coord[1];
		
		int count = 0;
		for (int j=0; j<ln[1]; j++) {
			for (int i=0; i<ln[0]; i++) {
				row[count] = igstart + i*inc[0] + j*inc[1];
				val[count] = RHSfunc2D(xcoord[istart+i], ycoord[jstart+j]);
				count++;
			}
		}
	}
	VecSetValues(*b, tln, row, val, ADD_VALUES);
	free(row);	
	free(val);	
}

void AssembleLevels(Grids *grids, Levels *levels) {
	// Build matrix "A" for a given level
	
	int	nlevels	= levels->nlevels;
	Level	*level	= levels->level;
	Mat	*res	= levels->res;
	Mat	*A	= levels->A;
	Vec	*b	= levels->b;
	Vec	*u	= levels->u;

	for (int l=0; l<nlevels; l++) {
		AssembleLevelMatA(grids, level+l, A+l);
		MatCreateVecs(A[l], b+l, u+l);
	}
	FillLevelVecb(0, grids->grid, level, b);
	ApplyBCLevelVecb(grids, level, b);
	VecAssemblyBegin(*b);
	VecAssemblyEnd(*b);
	
	for (int l=0; l<nlevels; l++) {
		int ngrids = level[l].ngrids;
		for (int lg=0; lg<ngrids-1; lg++) {
			int g = level[l].gridId[lg];
			AssembleMatRes(grids, lg, level+l, lg+1, level+l, res+g);
		}
	}

	for (int l=0; l<nlevels-1; l++) {
		int ngrids = level[l].ngrids;
		int g = level[l].gridId[ngrids-1];
		AssembleMatRes(grids, ngrids-1, level+l, 0, level+l+1, res+g);
	}
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
	ierr = PetscOptionsGetReal(NULL, NULL, "-rtol", &(solver->rtol), &set);
	if (!set) {
		PetscBarrier(PETSC_NULL);
		pERROR_MSG("Relative tolerance for solver is not set!");
		pERROR_MSG("Set '-rtol value'");
		return 1;
	}
	solver->rnorm = malloc((solver->numIter+1)*sizeof(double));
	solver->levels = malloc(sizeof(Levels));
	ierr = CreateLevels(grids, solver->levels); pCHKERR_RETURN("Levels creation failed");
	AssembleLevels(grids, solver->levels);
	return 0;
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

void GetSubIS(int lg, Grid *grid, Level *level, IS *indexSet) {
	// Get global index set for a grid in given level
	
	int dimension = grid->topo->dimension;

	int g	= level->gridId[lg];
	long int *inc = level->inc[lg];
	int igstart = level->ranges[lg];

	int *ln	= grid[g].ln;
	int tln = grid[g].tln;
	
	int *row = malloc(tln*sizeof(int));
	
	if (dimension == 3) {
		int count = 0;
		for (int k=0; k<ln[2]; k++) {
			for (int j=0; j<ln[1]; j++) {
				for (int i=0; i<ln[0]; i++) {
					row[count] = igstart + i*inc[0] + j*inc[1] + k*inc[2];
					count++;
				}
			}
		}
	} else if (dimension == 2) {
		int count = 0;
		for (int j=0; j<ln[1]; j++) {
			for (int i=0; i<ln[0]; i++) {
				row[count] = igstart + i*inc[0] + j*inc[1];
				count++;
			}
		}
	}
	
	ISCreateGeneral(PETSC_COMM_WORLD, tln, row, PETSC_COPY_VALUES, indexSet);
	free(row);
}


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

void GetSol(int lg, Grid *grid, Level *level, Vec *usol) {
	// Computes solution for a grid in a given level
	
	int *blockId = grid->topo->blockID;
	int dimension = grid->topo->dimension;

	int g	= level->gridId[lg];
	long int *inc = level->inc[lg];
//	int igstart = level->ranges[lg];
	int igstart, igend;
	VecGetOwnershipRange(*usol, &igstart, &igend);

	int *ln	= grid[g].ln;
	int tln = grid[g].tln;
	
	double *val = malloc(tln*sizeof(double));
	int *row = malloc(tln*sizeof(int));
	
	VecSet(*usol, 0.0);
	if (dimension == 3) {
		int istart = grid[g].range[0][blockId[0]];
		int jstart = grid[g].range[1][blockId[1]];
		int kstart = grid[g].range[2][blockId[2]];

		double *xcoord = grid[g].coord[0];
		double *ycoord = grid[g].coord[1];
		double *zcoord = grid[g].coord[2];
		
		int count = 0;
		for (int k=0; k<ln[2]; k++) {
			for (int j=0; j<ln[1]; j++) {
				for (int i=0; i<ln[0]; i++) {
					row[count] = igstart + i*inc[0] + j*inc[1] + k*inc[2];
					val[count] = Sol3D(xcoord[istart+i], ycoord[jstart+j], zcoord[kstart+k]);
					count++;
				}
			}
		}
	} else if (dimension == 2) {
		int istart = grid[g].range[0][blockId[0]];
		int jstart = grid[g].range[1][blockId[1]];

		double *xcoord = grid[g].coord[0];
		double *ycoord = grid[g].coord[1];
		
		int count = 0;
		for (int j=0; j<ln[1]; j++) {
			for (int i=0; i<ln[0]; i++) {
				row[count] = igstart + i*inc[0] + j*inc[1];
				val[count] = Sol2D(xcoord[istart+i], ycoord[jstart+j]);
				count++;
			}
		}
	}
	VecSetValues(*usol, tln, row, val, ADD_VALUES);
	free(row);	
	free(val);	
	VecAssemblyBegin(*usol);
	VecAssemblyEnd(*usol);
}

void GetError(Vec *u, Vec *usol, double *error) {
	// Get L_inf, L1 and L2 error norms
	// Note: usol is reused to store error!
	// Note: usol looses its original values!
	
	VecAXPY(*usol, -1.0, *u);
	VecNorm(*usol, NORM_INFINITY, error);
	VecNorm(*usol, NORM_1, error+1);
	VecNorm(*usol, NORM_2, error+2);
}

void GetGridVecFromLevelVec(int lg, Grid *grid, Level *level, Vec *u, Vec *ulg, IS *islg) {
	// Create a sub-vector for grid "lg" in the given level using level vector
	GetSubIS(lg, grid, level, islg);
	VecGetSubVector(*u, *islg, ulg);
}

void RestoreGridVecToLevelVec(Vec *u, Vec *ulg, IS *islg) {
	// Restore the sub-vector (grid vector) to level
	VecRestoreSubVector(*u, *islg, ulg);
	ISDestroy(islg);
}

int GetSubVecs(int n, IS *is, Vec *v, Vec **subv) {
	// Extract "n" sub-vectors from "v" using index sets "is"
	
	*subv = malloc(n*sizeof(Vec));
	int ierr = 0;
	for (int i=0; i<n; i++) {
//		ierr = VecCreate(PetscObjectComm((PetscObject)(*v)), (*subv)+i); pCHKERR_RETURN("Vec creation failed");
		ierr = VecGetSubVector(*v, is[i], (*subv)+i); pCHKERR_RETURN("Getting sub-vector failed");
	}
	return 0;
}

void RestoreSubVecs(int n, IS *is, Vec *v, Vec **subv) {
	// Restore "n" sub-vectors to "v" using index sets "is"
	
	for (int i=0; i<n; i++) {
		VecRestoreSubVector(*v, is[i], (*subv)+i);
	}
	free(*subv);
}

int CreateSubMats(int n, IS *is, Mat *A, Mat **subA) {
	// Extract "n" sub-vectors from "v" using index sets "is"
	
	*subA = malloc(n*sizeof(Mat));
	int ierr = 0;
	for (int i=0; i<n; i++) {
//		ierr = VecCreate(PetscObjectComm((PetscObject)(*v)), (*subv)+i); pCHKERR_RETURN("Vec creation failed");
		ierr = MatCreateSubMatrix(*A, is[i], is[i], MAT_INITIAL_MATRIX, (*subA)+i); pCHKERR_RETURN("Sub-matrix creation failed");
//		ierr = VecGetSubVector(*v, is[i], (*subv)+i); pCHKERR_RETURN("Getting sub-vector failed");
	}
	return 0;
}

void DestroySubMats(int n, Mat **subA) {
	// Restore "n" sub-vectors to "v" using index sets "is"
	
	for (int i=0; i<n; i++) {
		MatDestroy((*subA)+i);
	}
	free(*subA);
}

void WriteToFiles(Solver *solver) {
	
	FILE	*errData = fopen("eData.dat","w");
	double	*error = solver->error;
	for (int i=0;i<3;i++) {
		printf("error[%d] = %.16e\n", i, error[i]);
		fprintf(errData,"%.16e\n", error[i]);
	}
	printf("\n");
	fclose(errData);
	
	FILE	*resData = fopen("rData.dat","w");
	double	*rnorm = solver->rnorm;
	int	numIter = solver->numIter;	
	for (int i=0;i<numIter+1;i++) {
		fprintf(resData,"%.16e\n", rnorm[i]);
	}
	printf("Relative Residual norm = %.16e after %d iteration(s)\n\n", rnorm[numIter], numIter);
	fclose(resData);
}

int PostProcessing(Grids *grids, Solver *solver) {
	// Computes error and writes data to files
	
	Grid	*grid = grids->grid;
	Levels	*levels = solver->levels;
	Level	*level = levels->level; // First level

//	IS is0;
	int ierr = 0;
	IS *is = level->is;
	Vec usol;
	Vec *ugrid = NULL;
	Vec *u = levels->u;
	if (level->ngrids == 1) {
		ugrid = u;
	} else {
		// Extract finest grid computed solution from level vec
//		GetGridVecFromLevelVec(0, grid, level, u, ugrid, &is0);
		ugrid = malloc(sizeof(Vec));
		ierr = VecGetSubVector(*u, is[0], ugrid); pCHKERR_RETURN("Getting sub-vector failed");
	}
	VecDuplicate(*ugrid, &usol);
	GetSol(0, grid, level, &usol); // Get exact solution on finest grid in first level

	double *error = solver->error;
	GetError(ugrid, &usol, error);
//	if (level->ngrids > 1) RestoreGridVecToLevelVec(u, ugrid, &is0);
	if (level->ngrids > 1) {
		VecRestoreSubVector(*u, is[0], ugrid);
		free(ugrid);
	}
	VecDestroy(&usol);
	
	int	rank;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	
	if (rank == 0)	WriteToFiles(solver);
	return 0;
}

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

int NoMultigrid(Solver *solver) {
	// Solve using a solver without multigrid
	
	Levels	*levels = solver->levels;
	int	nlevels = levels->nlevels;
	int	ngrids = levels->level->ngrids;

	if (nlevels != 1 || ngrids != 1) {
		PetscBarrier(PETSC_NULL);
		pERROR_MSG("No more than 1 level or grid is allowed for No-MG solver");
		return 1;
	}
	
	int	rank;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	
	double	rtol	= solver->rtol;	
	double	*rnorm	= solver->rnorm;
	int	numIter	= solver->numIter;

	Mat	*A = levels->A;
	Vec	*b = levels->b;
	Vec	*u = levels->u;

	KSP	ksp;
	
	PetscLogStage	stageSolve;
	
	KSPCreate(PETSC_COMM_WORLD, &ksp);
	KSPSetType(ksp,KSPRICHARDSON);
	KSPSetOperators(ksp, *A, *A);
	KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED);
	KSPSetResidualHistory(ksp, rnorm, numIter, PETSC_FALSE);
	KSPSetTolerances(ksp, rtol, PETSC_DEFAULT, PETSC_DEFAULT, numIter);
	KSPSetFromOptions(ksp);
	
	// Solve the system
	PetscBarrier(PETSC_NULL);
	double initWallTime = MPI_Wtime();
	clock_t solverInitT = clock();
	PetscLogStageRegister("Solver", &stageSolve);
	PetscLogStagePush(stageSolve);
	
	KSPSolve(ksp, *b, *u);
	
	PetscLogStagePop();
	clock_t solverT = clock();
	double endWallTime = MPI_Wtime();
	PetscBarrier(PETSC_NULL);
	KSPGetIterationNumber(ksp, &(solver->numIter));

	double	rnorm0 = rnorm[0];
	for (int i=0;i<(numIter+1);i++) rnorm[i] = rnorm[i]/rnorm0;

	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = [%d]; Solver cputime:                %lf\n",rank,(double)(solverT-solverInitT)/CLOCKS_PER_SEC);
	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = [%d]; Solver walltime:               %lf\n",rank,endWallTime-initWallTime);
	PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
	
	PetscPrintf(PETSC_COMM_WORLD,"\n");

//	VecView(*u,PETSC_VIEWER_STDOUT_WORLD);
//	PetscPrintf(PETSC_COMM_WORLD,"---------------------------| level = 0 |------------------------\n");
//	KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
//	PetscPrintf(PETSC_COMM_WORLD,"----------------------------------------------------------------\n");
	KSPDestroy(&ksp);

	return 0;
}

static void CreateProMatsFromResMats(int dimension, int length, Mat *res, Mat **pro) {
	// Creates pro mats from res mats
	// Must be free the memory outside
	
	*pro = malloc(length*sizeof(Mat));
	double	scalar = 2.0;
	for (int i=1; i<dimension; i++) scalar *= 2.0;
	for (int l=0; l<length; l++) {
		MatTranspose(res[l], MAT_INITIAL_MATRIX, (*pro)+l);
		MatScale((*pro)[l], scalar);
	}
}

static void DestroyProMats(int length, Mat **pro) {
	// Destroy pro mats
	
	for (int l=0; l<length; l++) {
		MatDestroy((*pro)+l);
	}
	free(*pro);
}

int MultigridVcycle(Solver *solver) {
	// Solve using MG-V-Cyle solver
	
	Levels	*levels = solver->levels;
	int	nlevels = levels->nlevels;
	int	dimension = levels->dimension;

	int flag=0;
	for (int l=0; l<nlevels; l++)
		flag += levels->level[l].ngrids-1;

	if (nlevels == 1) {
		PetscBarrier(PETSC_NULL);
		pERROR_MSG("V-cycle MG solver requires minimum of 2 levels");
		return 1;
	}
	if (flag != 0) {
		PetscBarrier(PETSC_NULL);
		pERROR_MSG("Only one grid per level is allowed in V-cycle MG solver");
		return 1;
	}
	
	int	rank;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	
	double	rtol	= solver->rtol;
	double	*rnorm	= solver->rnorm;
	int	numIter	= solver->numIter;
	int	*v	= solver->v;

	Mat 	*res = levels->res;
	Mat	*A = levels->A;
	Vec	*b = levels->b;
	Vec	*u = levels->u;
	Mat 	*pro;
	
	CreateProMatsFromResMats(dimension, nlevels-1, res, &pro);

	KSP	ksp[nlevels];
	PC	temp;
//	PC	pc[levels];
	Vec	r[nlevels], rv[nlevels];
	
	PetscLogStage	stage;
	
	for (int i=0; i<nlevels; i++) {
		VecDuplicate(b[i], &(rv[i]));
	}
	
	for (int i=0; i<nlevels-1; i++) {
		KSPCreate(PETSC_COMM_WORLD, &(ksp[i]));
		KSPSetType(ksp[i],KSPRICHARDSON);
		KSPSetOperators(ksp[i], A[i], A[i]);
		PetscObjectSetOptionsPrefix((PetscObject)ksp[i], "levels_");
		KSPGetPC(ksp[i], &temp);
		PetscObjectSetOptionsPrefix((PetscObject)temp, "levels_");
		KSPSetNormType(ksp[i],KSP_NORM_NONE);
		KSPSetTolerances(ksp[i], rtol, PETSC_DEFAULT, PETSC_DEFAULT, v[0]);
		KSPSetFromOptions(ksp[i]);
	}

	KSPCreate(PETSC_COMM_WORLD, &(ksp[nlevels-1]));
	PetscObjectSetOptionsPrefix((PetscObject)ksp[nlevels-1], "coarse_");
	KSPGetPC(ksp[nlevels-1], &temp);
	PetscObjectSetOptionsPrefix((PetscObject)temp, "coarse_");
	KSPSetType(ksp[nlevels-1],KSPRICHARDSON);
	KSPSetOperators(ksp[nlevels-1], A[nlevels-1], A[nlevels-1]);
	KSPSetNormType(ksp[nlevels-1],KSP_NORM_NONE);
	KSPSetTolerances(ksp[nlevels-1], rtol, PETSC_DEFAULT, PETSC_DEFAULT, v[1]);
	KSPSetFromOptions(ksp[nlevels-1]);
	
	double bnorm, rnormmin, rnormmax, rnormchk;
	VecNorm(b[0], NORM_2, &bnorm);
	rnormmax = 100000000*bnorm;
	rnormmin = rtol*bnorm;

	VecSet(u[0], 0.0); // Note: This should be moved out of this function?
	MatMult(A[0], u[0], rv[0]);
	VecAXPY(rv[0], -1.0, b[0]);
	VecNorm(rv[0], NORM_2, &rnormchk);
	rnorm[0] = rnormchk;

	int iter = 0;
	
	// Solve the system
	PetscBarrier(PETSC_NULL);
	double initWallTime = MPI_Wtime();
	clock_t solverInitT = clock();
	PetscLogStageRegister("Solver", &stage);
	PetscLogStagePush(stage);
	while (iter<numIter && rnormmax > rnormchk && rnormchk > rnormmin) {
		KSPSolve(ksp[0], b[0], u[0]);
		if (iter==0) KSPSetInitialGuessNonzero(ksp[0],PETSC_TRUE);
		for (int l=1;l<nlevels;l++) {
			KSPBuildResidual(ksp[l-1],NULL,rv[l-1],&(r[l-1]));
			MatMult(res[l-1],r[l-1],b[l]);
			KSPSolve(ksp[l], b[l], u[l]);
			if (l!=nlevels-1) KSPSetInitialGuessNonzero(ksp[l],PETSC_TRUE);
		}
		for (int l=nlevels-2;l>=0;l=l-1) {
			MatMult(pro[l],u[l+1],rv[l]);
			VecAXPY(u[l],1.0,rv[l]);
			KSPSolve(ksp[l], b[l], u[l]);
			if (l!=0) KSPSetInitialGuessNonzero(ksp[l],PETSC_FALSE);
		}
		KSPBuildResidual(ksp[0],NULL,rv[0],&(r[0]));
		VecNorm(r[0], NORM_2, &rnormchk);	
		iter = iter + 1;
		rnorm[iter] = rnormchk;
	}
	PetscLogStagePop();
	clock_t solverT = clock();
	double endWallTime = MPI_Wtime();
	PetscBarrier(PETSC_NULL);
	rnormchk = rnorm[0];
	for (int i=0;i<(iter+1);i++) {
		rnorm[i] = rnorm[i]/rnormchk;
	}
	solver->numIter = iter;

	for (int i=0;i<nlevels;i++) {
		PetscPrintf(PETSC_COMM_WORLD,"---------------------------| level = %d |------------------------\n",i);
		KSPView(ksp[i],PETSC_VIEWER_STDOUT_WORLD);
		PetscPrintf(PETSC_COMM_WORLD,"-----------------------------------------------------------------\n");
	}
	
	DestroyProMats(nlevels-1, &pro);

	for (int i=0;i<nlevels;i++) {
		VecDestroy(&(rv[i]));
	}
	for (int i=0;i<nlevels;i++) {
		KSPDestroy(&(ksp[i]));
	}
	PetscPrintf(PETSC_COMM_WORLD,"\n");
	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = [%d]; Solver cputime:                %lf\n",rank,(double)(solverT-solverInitT)/CLOCKS_PER_SEC);
	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = [%d]; Solver walltime:               %lf\n",rank,endWallTime-initWallTime);
	PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
	
	return 0;
}

int MultigridAdditive(Solver *solver) {

	Levels	*levels = solver->levels;
	int	nlevels = levels->nlevels;
	int	dimension = levels->dimension;

	if (nlevels != 1) {
		PetscBarrier(PETSC_NULL);
		pERROR_MSG("Only one level is allowed in Additive MG solver");
		return 1;
	}
	
	int	ngrids = levels->level->ngrids;
	if (ngrids == 1) {
		PetscBarrier(PETSC_NULL);
		pERROR_MSG("Additive MG solver requires minimum of 2 grids");
		return 1;
	}

	int	rank;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

	double	rtol	= solver->rtol;
	double	*rnorm	= solver->rnorm;
	int	numIter	= solver->numIter;
	int	*v	= solver->v;

	Mat 	*res = levels->res;
	Mat	*A = levels->A;
	Mat	Afine;

	Vec	*b = levels->b;
	Vec	*u = levels->u;
	Vec	*subb, *subu;
	Vec	rfine;
	IS	*is = levels->level->is;
	Mat 	*pro;
	
	int	ierr=0;
	CreateProMatsFromResMats(dimension, ngrids-1, res, &pro);
	ierr = GetSubVecs(ngrids, is, b, &subb); pCHKERR_RETURN("Getting sub-vectors failed");
	ierr = GetSubVecs(ngrids, is, u, &subu); pCHKERR_RETURN("Getting sub-vectors failed");
	ierr = VecDuplicate(subb[0], &rfine); pCHKERR_RETURN("Vector duplication failed");
	ierr = MatCreateSubMatrix(*A, is[0], is[0], MAT_INITIAL_MATRIX, &Afine); pCHKERR_RETURN("Sub-matrix creation failed");

	KSP	ksp;
	
	PetscLogStage	stage;
		
	KSPCreate(PETSC_COMM_WORLD, &ksp);
	KSPSetType(ksp, KSPRICHARDSON);
	KSPSetOperators(ksp, *A, *A);
	KSPSetNormType(ksp, KSP_NORM_NONE);
	KSPSetTolerances(ksp, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, v[0]);
	KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
	KSPSetFromOptions(ksp);
	
	double	bnorm, rnormmin, rnormmax, rnormchk;
	VecNorm(subb[0], NORM_2, &bnorm);
	rnormmax = 100000000*bnorm;
	rnormmin = rtol*bnorm;

	VecSet(*u, 0.0); // Note: Should this be moved out of this function?
	MatResidual(Afine, subb[0], subu[0], rfine);
	VecNorm(rfine, NORM_2, &rnormchk);
	rnorm[0] = rnormchk;

	int iter = 0;
	
	PetscBarrier(PETSC_NULL);
	double initWallTime = MPI_Wtime();
	clock_t solverInitT = clock();
	PetscLogStageRegister("Solver", &stage);
	PetscLogStagePush(stage);
	// Note: Following algorithm is written assuming zer initial solution on coarse grids
	while (iter<numIter && rnormmax > rnormchk && rnormchk > rnormmin) {
		MatMult(res[0], rfine, subb[1]);
		for (int l=1;l<ngrids-1;l++) {
			MatMult(res[l], subb[l], subb[l+1]);
		}
		KSPSolve(ksp, *b, *u);
		for (int l=ngrids-2;l>0;l=l-1) {
//			MatMultAdd(pro[l], subu[l+1], subu[l], subu[l]);
			MatMult(pro[l], subu[l+1], subb[l]);
			VecAXPY(subu[l], 1.0, subb[l]);
			VecSet(subu[l+1], 0.0);
		}
//		MatMultAdd(pro[0], subu[1], subu[0], subu[0]);
		MatMult(pro[0], subu[1], rfine);
		VecAXPY(subu[0], 1.0, rfine);
		VecSet(subu[1], 0.0);
		MatResidual(Afine, subb[0], subu[0], rfine);
		VecNorm(rfine, NORM_2, &rnormchk);	
		iter = iter + 1;
		rnorm[iter] = rnormchk;
	}
	PetscLogStagePop();
	clock_t solverT = clock();
	double endWallTime = MPI_Wtime();
	rnormchk = rnorm[0];
	for (int i=0;i<(iter+1);i++) {
		rnorm[i] = rnorm[i]/rnormchk;
	}
	solver->numIter = iter;

	PetscPrintf(PETSC_COMM_WORLD,"---------------------------| level = 0 |------------------------\n");
	KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
	PetscPrintf(PETSC_COMM_WORLD,"----------------------------------------------------------------\n");

	MatDestroy(&Afine);
	VecDestroy(&rfine);
	RestoreSubVecs(ngrids, is, u, &subu);
	RestoreSubVecs(ngrids, is, b, &subb);
	DestroyProMats(ngrids-1, &pro);
	KSPDestroy(&ksp);

	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = [%d]; Solver cputime:                %lf\n",rank,(double)(solverT-solverInitT)/CLOCKS_PER_SEC);
	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = [%d]; Solver walltime:               %lf\n",rank,endWallTime-initWallTime);
	PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);

	return 0;
}

int MultigridAdditiveScaled(Solver *solver) {

	Levels	*levels = solver->levels;
	int	nlevels = levels->nlevels;
	int	dimension = levels->dimension;

	if (nlevels != 1) {
		PetscBarrier(PETSC_NULL);
		pERROR_MSG("Only one level is allowed in Additive-scaled MG solver");
		return 1;
	}
	
	int	ngrids = levels->level->ngrids;
	if (ngrids == 1) {
		PetscBarrier(PETSC_NULL);
		pERROR_MSG("Additive-scaled MG solver requires minimum of 2 grids");
		return 1;
	}

	int	rank;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	
	double	rtol	= solver->rtol;
	double	*rnorm	= solver->rnorm;
	int	numIter	= solver->numIter;
	int	*v	= solver->v;

	Mat 	*res = levels->res;
	Mat	*A = levels->A;
	Mat	*subA;

	Vec	*b = levels->b;
	Vec	*u = levels->u;
	Vec	*subb, *subu, *subr;
	Vec	rfine;
	Vec	r;
	IS	*is = levels->level->is;
	Mat 	*pro;
	
	int	ierr=0;
	CreateProMatsFromResMats(dimension, ngrids-1, res, &pro);
	ierr = VecDuplicate(*b, &r); pCHKERR_RETURN("Vector duplication failed");
	ierr = GetSubVecs(ngrids, is, &r, &subr); pCHKERR_RETURN("Getting sub-vectors failed");
	ierr = GetSubVecs(ngrids, is, b, &subb); pCHKERR_RETURN("Getting sub-vectors failed");
	ierr = GetSubVecs(ngrids, is, u, &subu); pCHKERR_RETURN("Getting sub-vectors failed");
	ierr = CreateSubMats(ngrids, is, A, &subA); pCHKERR_RETURN("Sub-matrices creation failed");
	ierr = VecDuplicate(subb[0], &rfine); pCHKERR_RETURN("Vector duplication failed");
//	ierr = MatCreateSubMatrix(*A, is[0], is[0], MAT_INITIAL_MATRIX, &Afine); pCHKERR_RETURN("Sub-matrix creation failed");
	
	KSP	ksp;
	
	PetscLogStage	stage;
		
	KSPCreate(PETSC_COMM_WORLD, &ksp);
	KSPSetType(ksp, KSPRICHARDSON);
	KSPSetOperators(ksp, *A, *A);
	KSPSetNormType(ksp, KSP_NORM_NONE);
	KSPSetTolerances(ksp, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, v[0]);
	KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
	KSPSetFromOptions(ksp);
	
	double	bnorm, rnormmin, rnormmax, rnormchk;
	double	lambda[ngrids], r0Dot[ngrids];
	VecNorm(subb[0], NORM_2, &bnorm);
	rnormmax = 100000000*bnorm;
	rnormmin = rtol*bnorm;

	VecSet(*u, 0.0); // Note: Should this be moved out of this function?
	MatResidual(subA[0], subb[0], subu[0], rfine);
	VecNorm(rfine, NORM_2, &rnormchk);
	rnorm[0] = rnormchk;
	r0Dot[0] = (rnormchk*rnormchk);

	int iter = 0;
	MPI_Comm comm;
	PetscObjectGetComm((PetscObject)(*b), &comm);
	
	PetscBarrier(PETSC_NULL);
	double initWallTime = MPI_Wtime();
	clock_t solverInitT = clock();
	PetscLogStageRegister("Solver", &stage);
	PetscLogStagePush(stage);
	while (iter<numIter && rnormmax > rnormchk && rnormchk > rnormmin) {
		MatMult(res[0], rfine, subb[1]);
		for (int l=1; l<ngrids-1; l++) {
			VecTDotBegin(subb[l], subb[l], r0Dot+l);
			MatMult(res[l], subb[l], subb[l+1]);
		}
		if(ngrids > 2) PetscCommSplitReductionBegin(comm);
		KSPSolve(ksp, *b, *u);
		for (int l=1; l<ngrids-1; l++) {
			VecTDotEnd(subb[l], subb[l], r0Dot+l);
		}
//		MatResidual(*A, *b, *u, r);
		if (ngrids > 2) {
			MatResidual(subA[ngrids-2], subb[ngrids-2], subu[ngrids-2], subr[ngrids-2]);
			VecTDotBegin(subb[ngrids-2], subr[ngrids-2], lambda+ngrids-2);
			PetscCommSplitReductionBegin(comm);
		}
		for (int l=ngrids-3; l>0; l=l-1) {
			MatResidual(subA[l], subb[l], subu[l], subr[l]);
			VecTDotEnd(subb[l+1], subr[l+1], lambda+l+1);
			VecTDotBegin(subb[l], subr[l], lambda+l);
			PetscCommSplitReductionBegin(comm);

			lambda[l+1] = lambda[l+1]/r0Dot[l+1];
			VecScale(subu[l+2], lambda[l+1]);
//			MatMultAdd(pro[l+1], subu[l+2], subu[l+1], subu[l+1]);
			MatMult(pro[l+1], subu[l+2], subr[l+1]);
			VecAXPY(subu[l+1], 1.0, subr[l+1]);
			VecSet(subu[l+2], 0.0);
		}
		MatResidual(subA[0], subb[0], subu[0], subr[0]);
		if (ngrids > 2) VecTDotEnd(subb[1], subr[1], lambda+1);
		VecTDotBegin(rfine, subr[0], lambda);
		PetscCommSplitReductionBegin(comm);

		if (ngrids > 2) {
			lambda[1] = lambda[1]/r0Dot[1];
			VecScale(subu[2], lambda[1]);
//			MatMultAdd(pro[1], subu[2], subu[1], subu[1]);
			MatMult(pro[1], subu[2], subr[1]);
			VecAXPY(subu[1], 1.0, subr[1]);
			VecSet(subu[2], 0.0);
		}
		VecTDotEnd(rfine, subr[0], lambda);

		lambda[0] = lambda[0]/r0Dot[0];
		VecScale(subu[1], lambda[0]);
//		MatMultAdd(pro[0], subu[1], subu[0], subu[0]);
		MatMult(pro[0], subu[1], subr[0]);
		VecAXPY(subu[0], 1.0, subr[0]);
		VecSet(subu[1], 0.0);
		
//		for (int l=ngrids-2;l>=0;l=l-1) {
//			lambda[l] = lambda[l]/r0Dot[l];
//			VecScale(subu[l+1], lambda[l]);
////			MatMultAdd(pro[l], subu[l+1], subu[l], subu[l]);
//			MatMult(pro[l], subu[l+1], subr[l]);
//			VecAXPY(subu[l], 1.0, subr[l]);
//			VecSet(subu[l+1], 0.0);
//		}

		MatResidual(subA[0], subb[0], subu[0], rfine);
		VecTDot(rfine, rfine, r0Dot);
		rnormchk = sqrt(r0Dot[0]);
//		VecNorm(rfine, NORM_2, &rnormchk);
		iter = iter + 1;
		rnorm[iter] = rnormchk;
	}
	PetscLogStagePop();
	clock_t solverT = clock();
	double endWallTime = MPI_Wtime();
	rnormchk = rnorm[0];
	for (int i=0;i<(iter+1);i++) {
		rnorm[i] = rnorm[i]/rnormchk;
	}
	solver->numIter = iter;

	PetscPrintf(PETSC_COMM_WORLD,"---------------------------| level = 0 |------------------------\n");
	KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
	PetscPrintf(PETSC_COMM_WORLD,"-----------------------------------------------------------------\n");
	
//	MatDestroy(&Afine);
	DestroySubMats(ngrids, &subA);
	RestoreSubVecs(ngrids, is, &r, &subr);
	RestoreSubVecs(ngrids, is, u, &subu);
	RestoreSubVecs(ngrids, is, b, &subb);
	DestroyProMats(ngrids-1, &pro);
	VecDestroy(&rfine);
	VecDestroy(&r);
	KSPDestroy(&ksp);

	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = [%d]; Solver cputime:                %lf\n",rank,(double)(solverT-solverInitT)/CLOCKS_PER_SEC);
	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank = [%d]; Solver walltime:               %lf\n",rank,endWallTime-initWallTime);
	PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
	
	return 0;
}

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
	
	int	ierr=0;
	PetscBarrier(PETSC_NULL);
	switch(solver->cycle) {
		case 0:
			ierr = NoMultigrid(solver); pCHKERR_RETURN("Iterative solver (No MG) failed");
			break;
		case 1:
			ierr = MultigridVcycle(solver); pCHKERR_RETURN("Multigrid V-cycle solver failed");
			break;
		case 2:
			ierr = MultigridAdditive(solver); pCHKERR_RETURN("Additive MG solver failed");
			break;
		case 3:
			ierr = MultigridAdditiveScaled(solver); pCHKERR_RETURN("Additive-scaled MG solver failed");
			break;
	}
	PetscBarrier(PETSC_NULL);
//	if (solver->cycle == ECYCLE) MultigridEcycle(solver);
//	if (solver->cycle == D1CYCLE) MultigridD1cycle(solver);
//	if (solver->cycle == D2CYCLE) MultigridD2cycle(solver);
//	if (solver->cycle == D1PSCYCLE) MultigridD1PScycle(solver);
//	if (solver->cycle == PetscPCMG) MultigridPetscPCMG(solver);
//	if (solver->cycle == FILTER) MultigridFilter(solver);
	
	return 0;	
}


