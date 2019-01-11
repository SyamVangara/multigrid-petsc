//#include "matbuild.h"
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

#define METRICS(i,j) (metrics.data[(i)*metrics.nj+(j)])
//#define METRICS(i,j,k) (metrics.data[metrics.nk*((i)*metrics.nj+(j))+(k)])
#define PMETRICS(i) ((metrics.data)+((i)*metrics.nj))
#define F(i,j) (f.data[((i)*f.nj+(j))])
#define U(i,j) (u.data[((i)*u.nj+(j))])

//static int ipow(int base, int exp) {
//	// Computes integer power of a real or integer
//	//
//	// base - real or int
//	// exp - int
//	// output: base^(exp)
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

void AssignGridID(Levels *levels, int ngrids) {
	// Assigns gridIDs for all the grids in each level
	
	int	nlevels = levels->nlevels;
	Level	*level = levels->level;
	
	int q = 0; 
	for (int l=0;l<nlevels;l++) {
		level[l].ngrids = 1;
		q += 1;
	}
	level[nlevels-1].ngrids += (ngrids-q);
	
	int count = 0;
	for (int l=0;l<nlevels;l++) {
		level[l].gridId = malloc(level[l].ngrids*sizeof(int));
		for (int g=0;g<level[l].ngrids;g++) {
			level[l].gridId[g] = count;
			count += 1;
		}
	}
}

//void CreateIndexMaps(Mesh *mesh, Indices *indices) {
//	// Allocate memory to global-to-grid and grid-to-global maps
//	
//	int	temp, g;
//	int	totaln, n0, n1, fn0, fn1;
//	int	factor;
//	const	int	M = 3; // i,j,g
//	
//	factor = indices->coarseningFactor;
//	fn0 = mesh->n[0];
//	fn1 = mesh->n[1];
//	for (int l=0;l<indices->levels;l++) {
//		totaln = 0;
//		for (int lg=0;lg<indices->level[l].grids;lg++) {
//			g = indices->level[l].gridId[lg];
//			temp = ipow(factor,g);
//			n0 = (fn0-1)/temp - 1;
//			n1 = (fn1-1)/temp - 1;
//			totaln = totaln + n0*n1; 
//			CreateArrayInt2d(n1, n0,&(indices->level[l].grid[lg]));
//		}
//		CreateArrayInt2d(totaln, M, &(indices->level[l].global));
//	}
//}
//
//void DestroyIndexMaps(Indices *indices) {
//	// Free the memory of global-to-grid and grid-to-global maps	
//	
//	for (int l=0;l<indices->levels;l++) {
//		for (int g=0;g<indices->level[l].grids;g++) {
//			DeleteArrayInt2d(&(indices->level[l].grid[g]));
//		}
//		DeleteArrayInt2d(&(indices->level[l].global));
//	}
//}

long int GetBlockStart(Grids *grids, Level *level, int *blockID) {
	// Compute the starting global index for given block
	
	int	dimension = grids->topo->dimension;
	Grid	*grid = grids->grid;
	int	ngrids = level->ngrids;
	int	*gridID = level->gridID;

	long int start = 0;
	for (int i=0; i<dimension; i++)	start = (blockID[i]<0 ? -1: 0);
	
	if (start<0) return start;

	int	crange[MAX_DIMENSION][2]; // Coordinate ranges
	for (int i=0; i<MAX_DIMENSION; i++)
		for (int j=0; j<2; j++)
			crange[i][j] = 0;
	
	for (int lg=0; lg<ngrids; lg++) {
		Grid	*grid = grid+gridID[lg];
		int	**range = grid->range;
		int	*un = grid->un;
		
		for (int dim=0; dim<dimension; dim++)
			for (int i=0; i<2; i++)
				crange[dim][i] = range[dim][blockID[dim]+i];
		
		long int lg0 = (crange[1][0]-1)*(un[0]) 
			+ (crange[0][0]-1)*(crange[1][1]-crange[1][0]);
		if (dimension == 3) {
			lg0 = lg0*(crange[2][1]-crange[2][0]) + (crange[2][0]-1)*un[1]*un[0];
		}
		start += lg0;
	}

	return start;
}

long int GetBlockGridStart(Grids *grids, Level *level, int *blockID, int targetlg) {
	// Compute the starting global index for given block and grid (targetlg: local grid ID)
	
	Grid	*grid = grids->grid;
	int	ngrids = level->ngrids;
	int	*gridID = level->gridId;
	int	dimension = grids->topo->dimension;

	long int start = GetBlockStart(grids, level, blockID);
	if (start<0) return start;
	for (int lg=0; lg<targetlg; lg++) {
		int **range =grid[gridID[lg]].range;
		long int temp = range[0][blockID[0]+1]-range[0][blockID[0]];
		for (int i=1; i<dimension; i++) {
			temp *= range[i][blockID[i]+1]-range[i][blockID[i]];
		}
		start += temp;
	}

	return start;
}

void GetRanges(Grids *grids, Level *level) {
//void GetBlockRanges(int dimension, Grids *grids, int ngrids, int *gridID, long int *ranges) {
	// Get global index range in this block for each grid in the list of gridID
	
	Grid	*grid = grids->grid;
	int	*blockid = grids->topo->blockID;
	int	ngrids = level->ngrids;
	int	*gridID = level->gridId;
	int	dimension = grids->topo->dimension;
	
	long int *ranges = level->ranges;

	ranges[0] = GetBlockStart(grids, level, blockid);
	for (int lg=1; lg<ngrids+1; lg++) {
		ranges[lg] = ranges[lg-1] + grid[gridID[lg-1]].tln;
	}	
//	int	crange[MAX_DIMENSION][2]; // Coordinate ranges
//	for (int i=0; i<MAX_DIMENSION; i++)
//		for (int j=0; j<2; j++)
//			crange[i][j] = 0;
//	const	int	nbc = 2; // Assuming 2 BC points per direction
//	
//	long int	g0 = 0; // Lower global index for this level's first grid
//	long int	*deltag = malloc(ngrids*sizeof(long int)); // No. of unknown points for each grid in this level
//
//	for (int lg=0; lg<ngrids; lg++) {
//		Grid	*grid = grids->grid+gridID[lg];
//		int	*n = grid->n;
//		
//		for (int dim=0; dim<dimension; dim++)
//			for (int i=0; i<2; i++)
//				crange[dim][i] = grid->range[dim][blockid[dim]+i];
//		
//		long int lg0 = (crange[1][0]-1)*(n[0]-nbc) 
//			+ (crange[0][0]-1)*(crange[1][1]-crange[1][0]);
//		deltag[lg] = (crange[0][1]-crange[0][0])*(crange[1][1]-crange[1][0]);
//		if (dimension == 3) {
//			lg0 = lg0*(crange[2][1]-crange[2][0]) 
//				+ (crange[2][0]-1)*(n[1]-nbc)*(n[0]-nbc);
//			deltag[lg] *= (crange[2][1]-crange[2][0]);
//		}
//		g0 += lg0;
//	}
//
//	ranges[0] = g0;
//	for (int lg=1; lg<ngrids; lg++)	ranges[lg] = ranges[lg-1] + deltag[lg-1];
//	free(deltag);
}

//void GetRanges(Grids *grids, Level *level) {
//	// Get range of global indices for each grid in a given level
//	// This assumes lexicographical ordering of grid points
//	
//	int	dimension = grids->topo->dimension;
//	int	*p = grids->topo->blockID;
//	int	lngrids = level->ngrids;
//	int	*gridID = level->gridId;
//	long int	(*grange)[2] = level->ranges; // Global index ranges
//	int	crange[MAX_DIMENSION][2]; // Coordinate ranges
//	for (int i=0; i<MAX_DIMENSION; i++)
//		for (int j=0; j<2; j++)
//			crange[i][j] = 0;
//	const	int	nbc = 2; // Assuming 2 BC points per direction
//	
//	long int	g0 = 0; // Lower global index for this level's first grid
//	long int	*deltag = malloc(lngrids*sizeof(long int)); // No. of unknown points for each grid in this level
//
//	for (int lg=0; lg<lngrids; lg++) {
//		Grid	*grid = grids->grid+gridID[lg];
//		int	*n = grid->n;
//		
//		for (int dim=0; dim<dimension; dim++)
//			for (int i=0; i<2; i++)
//				crange[dim][i] = grid->range[dim][p[dim]+i];
//		
//		long int lg0 = (crange[1][0]-1)*(n[0]-nbc) 
//			+ (crange[0][0]-1)*(crange[1][1]-crange[1][0]);
//		deltag[lg] = (crange[0][1]-crange[0][0])*(crange[1][1]-crange[1][0]);
//		if (dimension == 3) {
//			lg0 = lg0*(crange[2][1]-crange[2][0]) 
//				+ (crange[2][0]-1)*(n[1]-nbc)*(n[0]-nbc);
//			deltag[lg] *= (crange[2][1]-crange[2][0]);
//		}
//		g0 += lg0;
//	}
//
//	grange[0][0] = g0;
//	grange[0][1] = grange[0][0] + deltag[0];
//	for (int lg=1; lg<lngrids; lg++) {
//		grange[lg][0] = grange[lg-1][1];
//		grange[lg][1] = grange[lg][0] + deltag[lg];
//	}
//	free(deltag);
//}

void GetGridIncrements(int dimension, int *ln, long int *inc) {
	// Computes increments in each direction 
	inc[0] = 1;
	for (int i=1; i<dimension; i++)
		inc[i] = inc[i-1]*ln[i-1];
}

void GetBCindicesFromNblock(BCindices *bcindices, Nblock *nblock) {
	// Get global indices of neighboring block from the block's grid details
	

}

//void GetIncrements(Grids *grids, Level *level) {
//	
//	long int (*inc)[MAX_DIMENSION] = level->inc;
//	int *gridId = level->gridId;
//	int *blockID = grids->topo->blockID;
//	int *l = grids->topo->dimProcs;
//	int dimension = grids->topo->dimension;	
//	int ngrids = level->ngrids;
//
//	for (int lg=0; lg<ngrids; lg++) {
//		int g = gridId[lg]; // local to global grid index
////		int **range = grids->grid[g].range;
//		int *ln = grids->grid[g].ln;
//
//		GetGridIncrements(dimension, ln, inc[lg]);
////		inc[lg][0] = 1;
////		for (int i=1; i<dimension; i++) {
////			inc[lg][i] = inc[lg][i-1]*(range[i-1][blockID[i-1]+1]-range[i-1][blockID[i-1]]);
////		}
//	}
//	
////	long int *bcStartIndex = level->bcStartIndex;
////	long int *bcStart;
////	bcStart = malloc(ngrids*sizeof(long int));
////	
////	int tblockID[MAX_DIMENSION];
////	for (int i=0; i<MAX_DIMENSION; i++) tblockID[i] = blockID[i];
////
////	for (int i=0; i<dimension; i++)	{
////		for (int j=-1; j<2; j=j+2) {
////			tblockID[i] = tblockID[i] + j;
////			if (tblockID[i] >= 0 && tblockID[i] < l[i]) {
////				GetBlockRanges(grids, ngrids, gridId, dimension, tblockID, bcStart);
////			} else {
////				for (int lg=0; lg<ngrids; lg++) bcStart[lg] = -1;
////			}
////			long int bcinc[MAX_DIMENSION];
////			for (int lg=0; lg<ngrids; lg++) {
////				int g = gridId[lg];
////				
////				GetGridIncrements(dimension, tblockID, range, bcinc);
////				bcIncIndex[lg][i][(j+1)/2][]
////			}
////			tblockID[i] = tblockID[i] - j;
////			
////			for (int lg=0; lg<ngrids; lg++) 
////				bcStartIndex[lg][i][(j+1)/2] = bcStart[lg];
////		}
////	}
////	free(bcStart);
//}

void InitializeLevel(Level *level) {
	level->gridId	= NULL;
	level->ranges	= NULL;
	level->inc	= NULL;
	level->bcindices= NULL;
}

void InitializeLevels(Levels *levels) {
	levels->level = NULL;
}

void CreateBCindices(Grids *grids, Level *level, int dim, int dir, int targetlg) {
	
	Grid	*grid = grids->grid;
	int	g = level->gridID[targetlg];
	
	level->bcindices[targetlg][dim][dir].rank = grid[g].nblock[dim][dir].rank;
	int *iblockID = level->bcindices[targetlg][dim][dir].blockID;
	int *gblockID = grid[g].nblock[dim][dir].blockID;

	for (int i=0; i<MAX_DIMENSION; i++) iblockID[i] = gblockID[i]; 
	level->bcindices[targetlg][dim][dir].bcStartIndex = GetBlockGridStart(grids, level, iblockID, targetlg);
	int tinc[MAX_DIMENSION];
	int dimension = grids->topo->dimension;
	GetGridIncrements(dimension, grid[g].nblock[dim][dir].ln, tinc);
	
	int count = 0;
	for (int i=0; i<dimension; i++) {
		if (i != dim) {
			level->bcindices[targetlg][dim][dir].bcInc[count] = tinc[i];
			count++;
		} else {
			level->bcindices[targetlg][dim][dir].bcStartIndex += ((1-dir)*inc[i]);
		}
	}
}

void fillLevel(Grids *grids, Level *level) {
	// Fill the Level struct
	
	int	ngrids = level->ngrids;
	int	*gridId = level->gridId;
	int	dimension = grids->topo->dimension;
	Grid	*grid = grids->grid;

	level->ranges = malloc((ngrids+1)*sizeof(long int));
	GetRanges(grids, level);

	level->inc = malloc(ngrids*sizeof(long int[MAX_DIMENSION]));
	for (int lg=0; lg<ngrids; lg++) {
		GetGridIncrements(dimension, grid[gridId[lg]].ln, level->inc[lg]);
	}
//	GetIncrements(grids, level);
	
	level->bcindices = malloc((levels->level[i].ngrids)*sizeof(BCindices[MAX_DIMENSION][2]));
	for (int lg=0; lg<ngrids; lg++) {
		for (int i=0; i<dimension; i++) {
			for (int j=0; j<2; j++) {
				CreateBCindices(grids, level, i, j, lg);
			}
		}
	}
}

int CreateLevels(Grids *grids, Levels *levels) {
	// Get the GridId range, then allocate memory
	
	int	procs;
	MPI_Comm_size(PETSC_COMM_WORLD, &procs);
	
	if (!grids || !levels) {
		pERROR_MSG("NULL pointer encountered");
		return 1;
	}
	InitializeLevels(levels);
	
	int		ierr = 0;
	PetscBool	set;	
	ierr = PetscOptionsGetInt(NULL, NULL, "-nlevels", &(levels->nlevels), &set);
	if (!set || ierr) {
		PetscBarrier(PETSC_NULL);
		pERROR_MSG("Number of levels for MG solver not set");
		pERROR_MSG("Set '-nlevels n' for n levels");
		return 1;
	} else if (grids->ngrids < levels->nlevels) {
		PetscBarrier(PETSC_NULL);
		pERROR_MSG("Number of grids cannot be less than no. of levels");
		return 1;
	}
	
	levels->level = malloc(levels->nlevels*sizeof(Level));
	for (int i=0; i<levels->nlevels; i++) InitializeLevel(levels->level+i);

	AssignGridID(levels, grids->ngrids);
	for (int i=0;i<levels->nlevels;i++) {
		fillLevel(grids, levels->level+i);
//		levels->level[i].ranges = malloc((levels->level[i].ngrids)*sizeof(long int[2]));
//		GetRanges(grids, levels->level+i);
//		levels->level[i].inc = malloc((levels->level[i].ngrids)*sizeof(long int[MAX_DIMENSION]));
//		levels->level[i].bcIndex = malloc((levels->level[i].ngrids)*sizeof(long int[MAX_DIMENSION][2][3]));
//		levels->level[i].bcindices = malloc((levels->level[i].ngrids)*sizeof(BCindices[MAX_DIMENSION][2]));
//		GetIncrements(grids, levels->level+i);
	}
	return 0;
}

//void SetUpIndices(Mesh *mesh, Indices *indices) {
//	// Get the GridId range, then allocate memory
//	
//	int	procs;
//	MPI_Comm_size(PETSC_COMM_WORLD, &procs);
//	
//	indices->level = malloc(indices->levels*sizeof(Level));
//	GridId(indices);	
//	for (int i=0;i<indices->levels;i++) {
//		indices->level[i].ranges = malloc((procs+1)*sizeof(int));
//		indices->level[i].grid = malloc(indices->level[i].grids*sizeof(ArrayInt2d));
//		indices->level[i].h = malloc(indices->level[i].grids*sizeof(double[2]));
//	}
//	CreateIndexMaps(mesh, indices);
//	for (int l=0;l<indices->levels;l++) {
//		for (int lg=0;lg<indices->level[l].grids;lg++) {
//			indices->level[l].h[lg][0] = 1.0/(indices->level[l].grid[lg].ni+1);
//			indices->level[l].h[lg][1] = 1.0/(indices->level[l].grid[lg].nj+1);
//		}
//	}
//}

void DestroyLevels(Levels *levels) {
	// Free the memory allocated to indices
	
	if (!levels) return;
	if (levels->level) {
		for (int l=0;l<levels->nlevels;l++) {
			if (levels->level[l].gridId) free(levels->level[l].gridId);
			if (levels->level[l].ranges) free(levels->level[l].ranges);
			if (levels->level[l].inc) free(levels->level[l].inc);
			if (levels->level[l].bcindices) free(levels->level[l].bcindices);
		}
		free(levels->level);
	}
}

//void DestroyIndices(Indices *indices) {
//	// Free the memory allocated to indices
//	
////	DestroyIndexMaps(indices);
//	for (int l=0;l<indices->levels;l++) {
////		free(indices->level[l].grid);
////		free(indices->level[l].h);
//		free(indices->level[l].gridId);
//		free(indices->level[l].ranges);
//	}
//	free(indices->level);
//}

//void GetRanges(int *ranges, int totaln) {
//	// Computes the ranges of indices for all processes for a given total number of indices	
//	//
//	// length of ranges = procs + 1
//	// range[p] - index start in process "p"
//	// range[p+1] - (index end + 1) in process "p"
//	
//	int	remainder, quotient;
//	int	procs, rank;
//	
//	MPI_Comm_size(PETSC_COMM_WORLD, &procs);
//	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
//	
//	remainder = (totaln)%procs;
//	quotient  = (totaln)/procs;
//	ranges[0] = 0;
//	for (int p=0;p<procs;p++) {
//		if (p<remainder) {
//			ranges[p+1] = ranges[p] + (quotient + 1);
//		}
//		else {
//			ranges[p+1] = ranges[p] + quotient;
//		}
//	}
//}
//
//void mappingLocalGridAfterGrid(Indices *indices) {
//	// Maps global indices to grid unknowns and vice-versa
//	// Mapping style: Within the process points are arranged in a Grid-After-After pattern
//	
//	int	count;
//	int	grids, factor, gfactor;
//	int 	g;
//	int	*ranges, *gridranges;
//	ArrayInt2d	a, b, fine;
//	
//	int	procs, rank;
//	
//	MPI_Comm_size(PETSC_COMM_WORLD, &procs);
//	
//	factor = indices->coarseningFactor;
//	for (int l=0;l<indices->levels;l++) {
//		ranges 	= indices->level[l].ranges;
//		fine 	= indices->level[l].grid[0]; // Fine grid of level "l"
//		grids	= indices->level[l].grids;
//		GetRanges(ranges, fine.ni*fine.nj); // Domain decomposition on fine grid gives ranges of global indices for fine grid in each rank
//
//		// Compute the number of points of each grid in each rank and store them temporarily in "gridranges"
//		// Logically: gridranges[rank][g] - num of points of grid "g" in "rank"
//		gridranges = calloc((procs*grids)+1, sizeof(int));
//		count	= 0; // counts number of fine grid points
//		rank	= 0; // Keeps track of the rank
//		for (int i=0;i<fine.ni;i++) {
//			for (int j=0;j<fine.nj;j++) {
//				gfactor = 1;
//				for (int lg=0;lg<grids;lg++) {
//					// Check if the (i,j) has a corresponding point on grid "lg"
//					// If yes then count it towards num of points that grid has in the "rank"
//					if ((i+1)%gfactor!=0 || (j+1)%gfactor!=0) continue; 
//					gfactor =  gfactor*factor;
//					gridranges[rank*(grids)+lg] += 1;
//				}
//				count += 1;
//				if (count == ranges[rank+1]) {
//					// update the ending global index of ranges in "rank"
//					ranges[rank+1] = ranges[rank] + gridranges[rank*(grids)]; 
//					for (int lg=1;lg<grids;lg++) {
//						ranges[rank+1] += gridranges[rank*(grids)+lg];
//					}
//					rank += 1;
//				}
//			}
//		}
//	
//		// Then compute ranges for grids in each rank and use "gridranges" to store them
//		// Logically:	gridranges[rank][g] - starting index of grid "g" in rank
//		// 		gridranges[rank][g+1] - (ending index + 1) of grid "g" in rank	
//		gridranges[procs*grids] = ranges[procs];
//		for (int i=(procs*grids)-1;i>=0;i=i-1) {
//			gridranges[i] = gridranges[i+1]-gridranges[i];
//		}
//		
//		// Mapping
//		a = indices->level[l].global;
//		for (int lg=0;lg<indices->level[l].grids;lg++) {
//			g = indices->level[l].gridId[lg];
//			b = indices->level[l].grid[lg];
//			rank  = 0;
//			count = gridranges[lg];
//			for (int i=0;i<b.ni;i++) {
//				for (int j=0;j<b.nj;j++) {
//					while (count == gridranges[rank*grids+lg+1]) {
//						rank += 1;
//						count = gridranges[rank*grids+lg];
//					}
//					b.data[i*b.nj+j] = count;
//					a.data[count*a.nj+0] = i;
//					a.data[count*a.nj+1] = j;
//					a.data[count*a.nj+2] = g;
//					count = count + 1;
//				}
//			}
//		}
//		free(gridranges);
//	}
//}
//
//void mappingThroughGrids(Indices *indices) {
//	// Maps global indices to grid unknowns and vice-versa
//	// Mapping style: places points of any grid corresponding to same (x, y) physical coordinates next to each other
//	
//	int	count, fcount;
//	int	factor, gfactor;
//	int	ig, jg;
//	int	gfine, g;
//	int	*ranges;
//	ArrayInt2d	a, b, fine;
//	
//	int	procs, rank;
//	
//	MPI_Comm_size(PETSC_COMM_WORLD, &procs);
//	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
//	
//	factor = indices->coarseningFactor;
//	for (int l=0;l<indices->levels;l++) {
//		count = 0; // counts number of points mapped
//		a 	= indices->level[l].global;
//		ranges 	= indices->level[l].ranges;
//		gfine 	= indices->level[l].gridId[0];
//		fine 	= indices->level[l].grid[0];
//
//		// Compute the ranges of fine grid indices for processes
//		GetRanges(ranges, fine.ni*fine.nj);
//		fcount = 0; // counts number of fine grid points mapped
//		rank = 0; // Keeps track of the rank
//		for (int i=0;i<fine.ni;i++) {
//			for (int j=0;j<fine.nj;j++) {
//				for (int lg=0;lg<indices->level[l].grids;lg++) {
//					g = indices->level[l].gridId[lg];
//					gfactor = ipow(factor, g-gfine);
//					if ((i+1)%gfactor!=0 || (j+1)%gfactor!=0) continue;
//					ig = (i+1)/gfactor-1;
//					jg = (j+1)/gfactor-1;
//					b = indices->level[l].grid[lg];
//					b.data[ig*b.nj+jg] = count;
//					a.data[count*a.nj+0] = ig;
//					a.data[count*a.nj+1] = jg;
//					a.data[count*a.nj+2] = g;
//					count = count + 1;
//				}
//				fcount += 1;
//				if (fcount == ranges[rank+1]) {
//					ranges[rank+1] = count; // update the number of points in "rank"
//					rank += 1;
//				}
//			}
//		}
//	}
//}
//
//void mappingGridAfterGrid(Indices *indices) {
//	// Maps global indices to grid unknowns and vice-versa
//	// Mapping style: places all the points of a grid together
//	
//	int	count;
//	int	g;
//	int	*ranges;
//	ArrayInt2d	a, b;
//	
//	for (int l=0;l<indices->levels;l++) {
//		count = 0;
//		a = indices->level[l].global;
//		ranges 	= indices->level[l].ranges;
//		for (int lg=0;lg<indices->level[l].grids;lg++) {
//			g = indices->level[l].gridId[lg];
//			b = indices->level[l].grid[lg];
//			for (int i=0;i<b.ni;i++) {
//				for (int j=0;j<b.nj;j++) {
//					b.data[i*b.nj+j] = count;
//					a.data[count*a.nj+0] = i;
//					a.data[count*a.nj+1] = j;
//					a.data[count*a.nj+2] = g;
//					count = count + 1;
//				}
//			}
//		}
//		GetRanges(ranges, count);
//	}
//
//}
//
//void mapping(Indices *indices, int mappingStyleflag) {
//	// Maps global indices to grid unknowns and vice-versa
//	
//       	if (mappingStyleflag == 0) {
//		mappingGridAfterGrid(indices);
//	} else if (mappingStyleflag == 1) {
//		mappingThroughGrids(indices);
//	} else if (mappingStyleflag == 2) {
//		mappingLocalGridAfterGrid(indices);
//	} else {
//		printf("Unknown indices mapping style\n");
//	}
//}


//int CreateOperator(Grids *grids, Operator *op) {
//	// Allocates memory to Operator struct
//	int	order; // order of grid transfer operators
//	int	ngrids = grids->ngrids;
//	int	stencilSize;
//	
//	order = grids->cfactor; // order of grid transfer operators is same as the coarsening factor
//	op->res = malloc((ngrids-1)*sizeof(ArrayInt2d));
//	op->pro = malloc((ngrids-1)*sizeof(ArrayInt2d));
//	stencilSize = 1;
//	for (int i=0;i<ngrids-1;i++) {
//		stencilSize = (stencilSize+1)*order-1;
//		CreateArray2d(stencilSize, stencilSize, op->res+i);
//		CreateArray2d(stencilSize, stencilSize, op->pro+i);
//	}
//	
//	return 0;
//}

//int CreateOperator(Indices *indices, Operator *op) {
//	// Allocates memory to Operator struct
//	int	order; // order of grid transfer operators
//	int	stencilSize;
//	
//	order = indices->coarseningFactor; // order of grid transfer operators is same as the coarsening factor
//	op->totalGrids = indices->totalGrids;
//	op->res = malloc((op->totalGrids-1)*sizeof(ArrayInt2d));
//	op->pro = malloc((op->totalGrids-1)*sizeof(ArrayInt2d));
//	stencilSize = 1;
//	for (int i=0;i<op->totalGrids-1;i++) {
//		stencilSize = (stencilSize+1)*order-1;
//		CreateArray2d(stencilSize, stencilSize, op->res+i);
//		CreateArray2d(stencilSize, stencilSize, op->pro+i);
//	}
//
//}

//void DestroyOperator(Operator *op) {
//	// Free the memory in Operator struct
//	
//	for (int i=0;i<op->totalGrids-1;i++) {
//		DeleteArray2d(op->res+i);
//		DeleteArray2d(op->pro+i);
//	}
//	free(op->res);
//	free(op->pro);
//}

//void DestroyOperator(Operator *op) {
//	// Free the memory in Operator struct
//	
//	for (int i=0;i<op->totalGrids-1;i++) {
//		DeleteArray2d(op->res+i);
//		DeleteArray2d(op->pro+i);
//	}
//	free(op->res);
//	free(op->pro);
//}

//void GridTransferOperator(Array2d *Iop, int factor, int totalGrids) {
//	// Builds stencilwise grid transfer operator between any two grids that have "x" grids inbetween
//	// where x = {0,...,totalGrids-2}
//	
//	int 	ni0, nj0;
//	double	*weight;
//	int 	nil, njl;
//	double	*datal;
//	int 	niu, nju;
//	double	*datau;
//	int 	iu, ju;
//	
//	ni0 = Iop[0].ni;
//	nj0 = Iop[0].nj;
//	weight = Iop[0].data;
//	for (int l=0;l<totalGrids-2;l++) {
//		nil = Iop[l].ni;
//		njl = Iop[l].nj;
//		datal = Iop[l].data;
//		
//		niu = Iop[l+1].ni;
//		nju = Iop[l+1].nj;
//		datau = Iop[l+1].data;
//		// Initialization
//		for (int i=0;i<niu*nju;i++) {
//			datau[i] = 0.0;
//		}
//		// Building the operator
//		for (int il=0;il<nil;il++) {
//			for (int jl=0;jl<njl;jl++) {
//				iu = factor*(il+1)-1-ni0/2;
//				ju = factor*(jl+1)-1-nj0/2;
//				for (int i0=0;i0<ni0;i0++) {
//					for (int j0=0;j0<nj0;j0++) {
//						datau[(iu+i0)*nju+(ju+j0)] += weight[i0*nj0+j0]*datal[il*nil+jl];
//					}
//				}
//			}
//		}
//	}
//
//}
//
//void ProlongationOperator(Array2d pro) {
//	// Builds prolongation 2D stencilwise operator (pro)
//	// Stencil size: pro.ni x pro.nj
//	if(pro.ni != 3 || pro.nj != 3) {printf("Error in ProlongationOperator\n"); return;}
//	for (int i=0;i<pro.ni;i++) {
// 		pro.data[i*pro.nj]= 0.5 - 0.25*fabs(1-i);
// 		pro.data[i*pro.nj+1]= 1.0 - 0.5*fabs(1-i);
// 		pro.data[i*pro.nj+2]= 0.5 - 0.25*fabs(1-i);
//	}
//}

//void RestrictionOperator(Array2d res) {
//	// Builds Restriction 2D stencilwise operator (res)
//	// Stencil size: res.ni x res.nj
//	
//	if(res.ni != 3 || res.nj != 3) {printf("Error in RestrictionOperator\n"); return;}
//	for (int i=0;i<res.nj;i++) {
// 		res.data[i*res.nj]= 0.0;
// 		res.data[i*res.nj+1]= 0.0;
// 		res.data[i*res.nj+2]= 0.0;
//	}
//	res.data[res.nj+1] = 1.0;
//}

//void RestrictionOperator(Array2d res) {
//	// Builds Restriction 2D stencilwise operator (res)
//	// Stencil size: res.ni x res.nj
//	if(res.ni != 3 || res.nj != 3) {printf("Error in RestrictionOperator\n"); return;}
//	for (int i=0;i<res.nj;i++) {
// 		res.data[i*res.nj]= 0.125 - 0.0625*fabs(1-i);
// 		res.data[i*res.nj+1]= 0.25 - 0.125*fabs(1-i);
// 		res.data[i*res.nj+2]= 0.125 - 0.0625*fabs(1-i);
//	}
//}
//
//void GridTransferOperators(Operator op, Levels levels) {
//	// Builds stencilwise grid transfer operators between any two grids that have "x" grids inbetween
//	// where x = {0,...,levels-2}
//	if (op.totalGrids < 2) return;
//	RestrictionOperator(op.res[0]);
//	ProlongationOperator(op.pro[0]);
//	
//	GridTransferOperator(op.res, levels.coarseningFactor, op.totalGrids);
//	GridTransferOperator(op.pro, levels.coarseningFactor, op.totalGrids);
//}

//void Assemble(Problem *prob, Mesh *mesh, Indices *indices, Operator *op, Assembly *assem) {
//	// Assembles matrices and vectors in all levels
//	int	factor;
//
//	factor = indices->coarseningFactor;
//	for (int l=0;l<assem->levels;l++) {
//		levelMatrixA(prob, mesh, op, &(indices->level[l]), factor, assem->A+l);
//		MatCreateVecs(assem->A[l], assem->u+l, assem->b+l);
//	}
//	// Only the zeroth level vec b is created
//	levelvecb(prob, mesh, op, indices->level, factor, assem->b);
//
//	if (assem->levels < 2) return; 
//	Res(indices, op, factor, assem);
//	Pro(indices, op, factor, assem);
//}

