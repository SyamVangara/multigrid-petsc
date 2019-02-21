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

int LevelsGridIDsCheck(Levels *levels) {
	// Checks if levels have asecnding and non-ovelapping
	// grids in them. This is needed untill the Levels and 
	// Grids data structure is improved

	int	nlevels = levels->nlevels;
	Level	*level = levels->level;

	int pgd = -1;
	for (int l=0; l<nlevels; l++) {
		int *gridId = level[l].gridId;
		for (int lg=0; lg<level[l].ngrids; lg++) {
			if (gridId[lg] > pgd) {
				pgd = gridId[lg];
			} else {
				return 1;
			}
		}
	}
	return 0;
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

long int GetGridBlockStart(Grid *grid, int *blockID) {
	// Computes the starting global grid index for given block
	int	dimension = grid->topo->dimension;
	int	**range = grid->range;
	int	*un = grid->un;
	int	id = blockID[0];
	int	jd = blockID[1];
	int	kd = blockID[2];
	
//	int	crange[MAX_DIMENSION][2]; // Coordinate ranges
//	for (int dim=0; dim<dimension; dim++)
//		for (int i=0; i<2; i++)
//			crange[dim][i] = range[dim][blockID[dim]+i];
	
	long int start = (range[1][jd]-1)*(un[0]) 
		+ (range[0][id]-1)*(range[1][jd+1]-range[1][jd]);
	if (dimension == 3) {
		start = start*(range[2][kd+1]-range[2][kd]) + (range[2][kd]-1)*un[1]*un[0];
	}
	
//	long int start = (crange[1][0]-1)*(un[0]) 
//		+ (crange[0][0]-1)*(crange[1][1]-crange[1][0]);
//	if (dimension == 3) {
//		start = start*(crange[2][1]-crange[2][0]) + (crange[2][0]-1)*un[1]*un[0];
//	}

	return start;
}

long int GetBlockStart(Grids *grids, Level *level, int *blockID) {
	// Compute the starting global level index for given block
	
	int	dimension = grids->topo->dimension;
	int	ngrids = level->ngrids;
	int	*gridID = level->gridId;

	long int start = 0;
	for (int i=0; i<dimension; i++)	{
		if (blockID[i] < 0) {
			start = -1;
			break;
		}
	}
	
	if (start<0) return start;

	int	crange[MAX_DIMENSION][2]; // Coordinate ranges
	for (int i=0; i<MAX_DIMENSION; i++)
		for (int j=0; j<2; j++)
			crange[i][j] = 0;
	
	for (int lg=0; lg<ngrids; lg++) {
		Grid	*grid = grids->grid+gridID[lg];
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
	// Get global index range in this block for each grid in the list of gridID
	
	Grid	*grid = grids->grid;
	int	*blockid = grids->topo->blockID;
	int	ngrids = level->ngrids;
	int	*gridID = level->gridId;
	
	long int *ranges = level->ranges;

	ranges[0] = GetBlockStart(grids, level, blockid);
	for (int lg=1; lg<ngrids+1; lg++) {
		ranges[lg] = ranges[lg-1] + grid[gridID[lg-1]].tln;
	}	
}

void GetGridRanges(Grids *grids, Level *level) {
	// Get global grid index range in this block for each grid in this level
	
	Grid	*grid = grids->grid;
	int	*blockid = grids->topo->blockID;
	int	ngrids = level->ngrids;
	int	*gridID = level->gridId;
	
	long int (*granges)[2] = level->granges;

	for (int lg=0; lg<ngrids; lg++) {
		granges[lg][0] = GetGridBlockStart(grid+gridID[lg], blockid);
		granges[lg][1] = granges[lg][0] + grid[gridID[lg]].tln;
	}	
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
	level->granges	= NULL;
	level->inc	= NULL;
	level->bcindices= NULL;
	level->ebcindices= NULL;
	level->cbcindices= NULL;
	level->is	= NULL;
}

void InitializeLevels(Levels *levels) {
	levels->A = NULL;
	levels->b = NULL;
	levels->u = NULL;
	levels->level = NULL;
}

void CreateBCindicesFromNblock(int targetlg, int dim, int dir, Grids *grids, Level *level, Nblock *nblock, BCindices *bcindices) {
	// Creates BCindices corresponding to a Nblock for targeted grid
	// i.e., create neighbor grid indices info from neighbor grid block info
	
	Grid	*grid = grids->grid;
	int	g = level->gridId[targetlg];
	
	int	*iln = bcindices->ln;
	int	*gln = nblock->ln;
	int	*iblockID = bcindices->blockID;
	int	*gblockID = nblock->blockID;
	long int *startindex = &(bcindices->bcStartIndex);
	long int *gstartindex = &(bcindices->bcGStartIndex);
	
	bcindices->sbcindices = NULL;	
	bcindices->rank = nblock->rank;
	for (int i=0; i<MAX_DIMENSION; i++) {
		iblockID[i] = gblockID[i]; 
		iln[i] = gln[i];
	}
	*startindex = GetBlockGridStart(grids, level, iblockID, targetlg);
	if (*startindex >=0) {
		*gstartindex = GetGridBlockStart(grid+g, iblockID);
	} else {
		*gstartindex = -1;
	}
	
	int dimension = grids->topo->dimension;
	long int *inc = bcindices->bcInc;
	int *ln = nblock->ln;
	GetGridIncrements(dimension, ln, inc);
	
	if (*startindex >= 0) {
		*startindex += ((1-dir)*inc[dim]*(ln[dim]-1));
		*gstartindex += ((1-dir)*inc[dim]*(ln[dim]-1));
	} else {
		for(int i=0; i<dimension; i++) inc[i] = 0;
	}
	
}

void CreateBCindices(Grids *grids, Level *level, int dim, int dir, int targetlg) {
	
	Grid	*grid = grids->grid;
	int	g = level->gridId[targetlg];
	
	CreateBCindicesFromNblock(targetlg, dim, dir, grids, level, 
		&(grid[g].nblock[dim][dir]), &(level->bcindices[targetlg][dim][dir]));
//	level->bcindices[targetlg][dim][dir].rank = grid[g].nblock[dim][dir].rank;
//	int *iblockID = level->bcindices[targetlg][dim][dir].blockID;
//	int *gblockID = grid[g].nblock[dim][dir].blockID;
//	long int *startindex = &(level->bcindices[targetlg][dim][dir].bcStartIndex);
//	long int *gstartindex = &(level->bcindices[targetlg][dim][dir].bcGStartIndex);
//
//	for (int i=0; i<MAX_DIMENSION; i++) iblockID[i] = gblockID[i]; 
//	*startindex = GetBlockGridStart(grids, level, iblockID, targetlg);
//	if (*startindex >=0) {
//		*gstartindex = GetGridBlockStart(grid+g, iblockID);
//	} else {
//		*gstartindex = -1;
//	}
//
//	int dimension = grids->topo->dimension;
//	long int *inc = level->bcindices[targetlg][dim][dir].bcInc;
//	int *ln = grid[g].nblock[dim][dir].ln;
//	GetGridIncrements(dimension, ln, inc);
//	
//	if (*startindex >= 0) {
//		*startindex += ((1-dir)*inc[dim]*(ln[dim]-1));
//		*gstartindex += ((1-dir)*inc[dim]*(ln[dim]-1));
//	} else {
//		for(int i=0; i<dimension; i++) inc[i] = 0;
//	}
	
	// Create second neighbor indices
	level->bcindices[targetlg][dim][dir].sbcindices = malloc(sizeof(BCindices));
	BCindices	*secbcindices = level->bcindices[targetlg][dim][dir].sbcindices;
	Nblock		*snblock = grid[g].nblock[dim][dir].snblock;
	
	CreateBCindicesFromNblock(targetlg, dim, dir, grids, level, snblock, secbcindices);
//	secbcindices->rank = snblock->rank;
//	int *siblockID = secbcindices->blockID;
//	int *sgblockID = snblock->blockID;
//	long int *sstartindex = &(secbcindices->bcStartIndex);
//	long int *sgstartindex = &(secbcindices->bcGStartIndex);
//	
//	for (int i=0; i<MAX_DIMENSION; i++) siblockID[i] = sgblockID[i]; 
//	*sstartindex = GetBlockGridStart(grids, level, siblockID, targetlg);
//	if (*sstartindex >=0) {
//		*sgstartindex = GetGridBlockStart(grid+g, siblockID);
//	} else {
//		*sgstartindex = -1;
//	}
//	long int *sinc = secbcindices->bcInc;
//	int *sln = snblock->ln;
//	GetGridIncrements(dimension, sln, sinc);
//	
//	if (*sstartindex >= 0) {
//		*startindex += ((1-dir)*sinc[dim]*(sln[dim]-1));
//		*sgstartindex += ((1-dir)*sinc[dim]*(sln[dim]-1));
//	} else {
//		for(int i=0; i<dimension; i++) sinc[i] = 0;
//	}
	
}

void CreateEBCindices(Grids *grids, Level *level, int dim, int jdir, int kdir, int targetlg) {
	
	Grid	*grid = grids->grid;
	int	g = level->gridId[targetlg];

	Nblock		*eblock = &(grid[g].eblock[dim][jdir][kdir]);
	BCindices	*ebcindices = &(level->ebcindices[targetlg][dim][jdir][kdir]);

	ebcindices->rank = eblock->rank;
	int	*iln = ebcindices->ln;
	int	*gln = eblock->ln;
	int *iblockID = ebcindices->blockID;
	int *gblockID = eblock->blockID;
	long int *startindex = &(ebcindices->bcStartIndex);
	long int *gstartindex = &(ebcindices->bcGStartIndex);

	for (int i=0; i<MAX_DIMENSION; i++) {
		iblockID[i] = gblockID[i]; 
		iln[i] = gln[i];
	}
	*startindex = GetBlockGridStart(grids, level, iblockID, targetlg);
	if (*startindex >=0) {
		*gstartindex = GetGridBlockStart(grid+g, iblockID);
	} else {
		*gstartindex = -1;
	}

	int dimension = grids->topo->dimension;
	long int *inc = ebcindices->bcInc;
	int *ln = eblock->ln;
	GetGridIncrements(dimension, ln, inc);
	
	if (*startindex >= 0) {
		int jdim = (dim+1)%3;
		int kdim = (dim+2)%3;
		long int temp = ((1-jdir)*inc[jdim]*(ln[jdim]-1)+(1-kdir)*inc[kdim]*(ln[kdim]-1));
		*startindex += temp;
		*gstartindex += temp;
	} else {
		for(int i=0; i<dimension; i++) inc[i] = 0;
	}
}

void CreateCBCindices(Grids *grids, Level *level, int idir, int jdir, int kdir, int targetlg) {
	
	Grid	*grid = grids->grid;
	int	g = level->gridId[targetlg];
	
	Nblock		*cblock = &(grid[g].cblock[idir][jdir][kdir]);
	BCindices	*cbcindices = &(level->cbcindices[targetlg][idir][jdir][kdir]);

	cbcindices->rank = cblock->rank;
	int	*iln = cbcindices->ln;
	int	*gln = cblock->ln;
	int *iblockID = cbcindices->blockID;
	int *gblockID = cblock->blockID;
	long int *startindex = &(cbcindices->bcStartIndex);
	long int *gstartindex = &(cbcindices->bcGStartIndex);

	for (int i=0; i<MAX_DIMENSION; i++) {
		iblockID[i] = gblockID[i]; 
		iln[i] = gln[i];
	}
	*startindex = GetBlockGridStart(grids, level, iblockID, targetlg);
	if (*startindex >=0) {
		*gstartindex = GetGridBlockStart(grid+g, iblockID);
	} else {
		*gstartindex = -1;
	}

	int dimension = grids->topo->dimension;
	long int *inc = cbcindices->bcInc;
	int *ln = cblock->ln;
	GetGridIncrements(dimension, ln, inc);
	
	if (*startindex >= 0) {
		long int temp = ((1-idir)*inc[0]*(ln[0]-1)+(1-jdir)*inc[1]*(ln[1]-1));
		if (dimension == 3) temp += (1-kdir)*inc[2]*(ln[2]-1);
		*startindex += temp;
		*gstartindex += temp;
	} else {
		for(int i=0; i<dimension; i++) inc[i] = 0;
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
	level->granges = malloc((ngrids)*sizeof(long int[2]));
	GetGridRanges(grids, level);

	level->inc = malloc(ngrids*sizeof(long int[MAX_DIMENSION]));

	for (int lg=0; lg<ngrids; lg++) {
		GetGridIncrements(dimension, grid[gridId[lg]].ln, level->inc[lg]);
	}
	
	if (ngrids > 1) {
		level->is  = malloc(ngrids*sizeof(IS));
		for (int lg=0; lg<ngrids; lg++) {
			GetSubIS(lg, grid, level, level->is+lg);
		}
	}
	
	level->bcindices = malloc(ngrids*sizeof(BCindices[MAX_DIMENSION][2]));
	for (int lg=0; lg<ngrids; lg++) {
		for (int i=0; i<dimension; i++) {
			for (int j=0; j<2; j++) {
				CreateBCindices(grids, level, i, j, lg);
			}
		}
	}

	level->cbcindices = malloc(ngrids*sizeof(BCindices[2][2][2]));
	for (int lg=0; lg<ngrids; lg++) {
		for (int i=0; i<2; i++) {
			for (int j=0; j<2; j++) {
				for (int k=0; k<dimension-1; k++) {
					CreateCBCindices(grids, level, i, j, k, lg);
				}
			}
		}
	}
	
	if (dimension == 2) return;
	level->ebcindices = malloc(ngrids*sizeof(BCindices[MAX_DIMENSION][2][2]));
	for (int lg=0; lg<ngrids; lg++) {
		for (int dim=0; dim<3; dim++) {
			for (int j=0; j<2; j++) {
				for (int k=0; k<2; k++) {
					CreateEBCindices(grids, level, dim, j, k, lg);
				}
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
	levels->dimension = grids->topo->dimension;	
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
	ierr = LevelsGridIDsCheck(levels);
	if (ierr) {
		PetscBarrier(PETSC_NULL);
		pERROR_MSG("Grids are not assigned in a strictly ascending manner to Levels");
		return 1;
	}
	levels->res= malloc((grids->ngrids-1)*sizeof(Mat));
	levels->A = malloc(levels->nlevels*sizeof(Mat));
	levels->b = malloc(levels->nlevels*sizeof(Vec));
	levels->u = malloc(levels->nlevels*sizeof(Vec));
	int	prob = levels->prob;
	double	eps = levels->eps;
	for (int i=0;i<levels->nlevels;i++) {
		fillLevel(grids, levels->level+i);
		levels->level[i].prob = prob;
		levels->level[i].eps = eps;
	}
	return 0;
}

void DestroyLevels(Levels *levels) {
	// Free the memory allocated to indices
	
	if (!levels) return;
	if (levels->res) {
		for (int l=0; l<levels->nlevels; l++) {
			int ngrids = levels->level[l].ngrids;
			for (int lg=0; lg<ngrids-1; lg++) {
				int g = levels->level[l].gridId[lg];
				MatDestroy(levels->res+g);
			}
		}
		for (int l=0; l<levels->nlevels-1; l++) {
			int ngrids = levels->level[l].ngrids;
			int g = levels->level[l].gridId[ngrids-1];
			MatDestroy(levels->res+g);
		}
		free(levels->res);
	}
	if (levels->A) {
		for (int l=0; l<levels->nlevels; l++) {
			MatDestroy(levels->A+l);
		}
		free(levels->A);
	}
	if (levels->b) {
		for (int l=0; l<levels->nlevels; l++) {
			VecDestroy(levels->b+l);
		}
		free(levels->b);
	}
	if (levels->u) {
		for (int l=0; l<levels->nlevels; l++) {
			VecDestroy(levels->u+l);
		}
		free(levels->u);
	}
	if (levels->level) {
		for (int l=0;l<levels->nlevels;l++) {
			if (levels->level[l].gridId) free(levels->level[l].gridId);
			if (levels->level[l].ranges) free(levels->level[l].ranges);
			if (levels->level[l].granges) free(levels->level[l].granges);
			if (levels->level[l].inc) free(levels->level[l].inc);
			if (levels->level[l].bcindices) {
				for (int lg=0; lg<levels->level[l].ngrids; lg++) {
				for (int dim=0; dim<levels->dimension; dim++) {
				for (int dir=0; dir<2; dir++) {
					BCindices *tbcindices = &(levels->level[l].bcindices[lg][dim][dir]);
					if(tbcindices->sbcindices) free(tbcindices->sbcindices);
				}
				}
				}
				free(levels->level[l].bcindices);
			}
			if (levels->level[l].ebcindices) free(levels->level[l].ebcindices);
			if (levels->level[l].cbcindices) free(levels->level[l].cbcindices);
			if (levels->level[l].is) {
				int ngrids = levels->level[l].ngrids;
				for (int lg=0; lg<ngrids; lg++) {
					ISDestroy(levels->level[l].is+lg);
				}
				free(levels->level[l].is);
			}
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

