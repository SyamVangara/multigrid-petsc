#include "mesh.h"

#define ERROR_MSG(message) (fprintf(stderr,"ERROR: %s:%d: %s\n",__FILE__,__LINE__,(message)))
#define ERROR_RETURN(message) {ERROR_MSG(message);return ierr;}
#define CHKERR_PRNT(message) {if(ierr==1) {ERROR_MSG(message);}}
#define CHKERR_RETURN(message) {if(ierr==1) {ERROR_RETURN(message);}}
#define PI 3.14159265358979323846

#define pERROR_MSG(message) (PetscPrintf(PETSC_COMM_WORLD,"ERROR: %s:%d: %s\n",__FILE__,__LINE__,(message)))
#define pERROR_RETURN(message) {pERROR_MSG(message);return ierr;}
#define pCHKERR_PRNT(message) {if(ierr==1) {pERROR_MSG(message);}}
#define pCHKERR_RETURN(message) {if(ierr==1) {pERROR_RETURN(message);}}

#define METRICS(i,j,k) (metrics->data[metrics->nk*((i)*metrics->nj+(j))+(k)])
#define isGRIDtoGLOBAL(l,i,j) (IsGridToGlobal[l].data[((i)*IsGridToGlobal[l].nj+(j))])
#define isGLOBALtoGRID(l,i,j) (IsGlobalToGrid[l].data[((i)*IsGlobalToGrid[l].nj+(j))])

static int CommCost3D(int *l, int *n) {
	// Computes communicstion cost based on number of points in each direction
	// and number of directionwise ranks/blocks
	
	return (l[0]-1)*n[1]*n[2] + n[0]*(l[1]-1)*n[2] + n[0]*n[1]*(l[2]-1); 
}

static int ComputeNInterfaces3D(int *l) {
	// Compute number of interfaces to communicate data
		
	return (l[0]-1)*l[1]*l[2] + l[0]*(l[1]-1)*l[2] + l[0]*l[1]*(l[2]-1); 
}

static void ComputeLoad(int *n, int *l, double *load) {
	// Compute the maximum and minimum loads, and corresponding load factors
	
	int ones[3], q[3];

	for (int i=0; i<3; i++) {
		q[i] = n[i]/l[i];
		if (q[i]*l[i] == n[i])
			ones[i] = 0;
		else
			ones[i] = 1;
	}

	load[0] = (q[0]+ones[0])*(q[1]+ones[1])*(q[2]+ones[2]);
	load[1] = q[0]*q[1]*q[2];
	load[2] = (double)load[0]/(double)load[1];
}

static void sort3Index(int *n, int *nsorted, int *index) {
	// Sorts three integers in ascending order
	// Alos, outputs indices of new positions
	
	if (n[0]>n[1]) {
		if (n[0]>n[2]) {
			index[0] = 2;
			nsorted[2] = n[0];
			if (n[1] > n[2]) {
				index[1] = 1;
				index[2] = 0;
			        nsorted[1] = n[1];
			        nsorted[0] = n[2];
			} else {
				index[2] = 1;
				index[1] = 0;
			        nsorted[1] = n[2];
			        nsorted[0] = n[1];
			}
		} else {
			index[0] = 1;
			index[2] = 2;
			index[1] = 0;
			nsorted[1] = n[0];
			nsorted[2] = n[2];
			nsorted[0] = n[1];
		}
	} else {
		if (n[1]>n[2]) {
			index[1] = 2;
			nsorted[2] = n[1];
			if (n[0]>n[2]) {
				index[0] = 1;
				index[2] = 0;
				nsorted[1] = n[0];
				nsorted[0] = n[2];
			} else {
				index[2] = 1;
				index[0] = 0;
				nsorted[1] = n[2];
				nsorted[0] = n[0];
			}
		} else {
			index[0] = 0;
			index[2] = 2;
			index[1] = 1;
			nsorted[1] = n[1];
			nsorted[2] = n[2];
			nsorted[0] = n[0];
		}
	}
}

static void sort3(int *n, int *nsorted) {
	// Sorts three integers in ascending order
	
	if (n[0]>n[1]) {
		if (n[0]>n[2]) {
			nsorted[2] = n[0];
			if (n[1] > n[2]) {
			        nsorted[1] = n[1];
			        nsorted[0] = n[2];
			} else {
			        nsorted[1] = n[2];
			        nsorted[0] = n[1];
			}
		} else {
			nsorted[1] = n[0];
			nsorted[2] = n[2];
			nsorted[0] = n[1];
		}
	} else {
		if (n[1]>n[2]) {
			nsorted[2] = n[1];
			if (n[0]>n[2]) {
				nsorted[1] = n[0];
				nsorted[0] = n[2];
			} else {
				nsorted[1] = n[2];
				nsorted[0] = n[0];
			}
		} else {
			nsorted[1] = n[1];
			nsorted[2] = n[2];
			nsorted[0] = n[0];
		}
	}
}

//void MetricsUniform(void *mesh, double x, double y, double *metrics) {
//	//Computes following metrics at (x,y)
//	//
//	//metrics[0] = (xi_x)^2 + (xi_y)^2
//	//metrics[1] = (eta_x)^2 + (eta_y)^2
//	//metrics[2] = (xi_xx) + (xi_yy)
//	//metrics[3] = (eta_xx) + (eta_yy)
//	//metrics[4] = (xi_x)(eta_x) + (xi_y)(eta_y)
//	
//	metrics[0] = 1.0;
//	metrics[1] = 1.0;
//	metrics[2] = 0.0;
//	metrics[3] = 0.0; 
//	metrics[4] = 0.0;
//}
//
//void MetricsNonUniform1(void *mesh1, double x, double y, double *metrics) {
//	//Computes following metrics at (x,y)
//	//
//	//metrics[0] = (xi_x)^2 + (xi_y)^2
//	//metrics[1] = (eta_x)^2 + (eta_y)^2
//	//metrics[2] = (xi_xx) + (xi_yy)
//	//metrics[3] = (eta_xx) + (eta_yy)
//	//metrics[4] = (xi_x)(eta_x) + (xi_y)(eta_y)
//	//
//	//bounds[0] - Lower bound of "x"
//	//bounds[1] - Upper bound of "x"
//	//bounds[2] - Lower bound of "y"
//	//bounds[3] - Upper bound of "y"
//	
//	double	temp;
//	double	*bounds;
//	Mesh	*mesh;
//	int	procs, rank;
//	
//	MPI_Comm_size(PETSC_COMM_WORLD, &procs);
//	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
//	
//	mesh = (Mesh*)mesh1;
//	bounds = mesh->bounds;
//	temp = ((bounds[3]-bounds[2])*(bounds[3]-bounds[2])-(bounds[3]-y)*(bounds[3]-y));
//	metrics[0] = 1.0;
//	metrics[1] = 4.0/(PI*PI*temp);
//	metrics[2] = 0.0;
//	metrics[3] = (-2.0*(bounds[3]-y))/(PI*sqrt(temp*temp*temp)); 
//	metrics[4] = 0.0;
//}
//
//void MetricsNonUniform2(void *mesh1, double x, double y, double *metrics) {
//	//Computes following metrics at (x,y)
//	//
//	//metrics[0] = (xi_x)^2 + (xi_y)^2
//	//metrics[1] = (eta_x)^2 + (eta_y)^2
//	//metrics[2] = (xi_xx) + (xi_yy)
//	//metrics[3] = (eta_xx) + (eta_yy)
//	//metrics[4] = (xi_x)(eta_x) + (xi_y)(eta_y)
//	//
//	//bounds[0] - Lower bound of "x"
//	//bounds[1] - Upper bound of "x"
//	//bounds[2] - Lower bound of "y"
//	//bounds[3] - Upper bound of "y"
//	
//	double	temp;
//	double	*bounds;
//	Mesh	*mesh;
//	int	procs, rank;
//	
//	MPI_Comm_size(PETSC_COMM_WORLD, &procs);
//	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
//	
//	mesh = (Mesh*)mesh1;
//	bounds = mesh->bounds;
//	temp = ((exp(2)-1)*(exp(2)-1))/(((y-bounds[2])*(exp(2)-1)+(bounds[3]-bounds[2]))*((y-bounds[2])*(exp(2)-1)+(bounds[3]-bounds[2])));
//	metrics[0] = 1.0/((bounds[1]-bounds[0])*(bounds[1]-bounds[0]));
//	metrics[1] = 0.25*temp;
//	metrics[2] = 0.0;
//	metrics[3] = (-0.5)*temp; 
//	metrics[4] = 0.0;
//}

int Compute_coord(Grid *grid) {
	// computes coords in each direction of the structured grid
	
	int	ierr = 0;
	double	length, d[MAX_DIMENSION];

	int	*n;
	double	**coord;
	int	*type;
	int	dimension = grid->topo->dimension;
	Topo	*topo = grid->topo;
	n 	= grid->n;
	coord	= grid->coord;
	type	= grid->topo->gridtype;
	
	for(int i=0;i<dimension;i++){
		if (n[i]<MIN_POINTS) {
			PetscPrintf(PETSC_COMM_WORLD,
				"ERROR: %s:%d: Need at least %d points in each direction\n",
				__FILE__,__LINE__,MIN_POINTS);
			return 1;
		}
		coord[i][0] = topo->bounds[i][0]; //Lower bound
		coord[i][n[i]-1] = topo->bounds[i][1]; //Upper bound
		
		length = (coord[i][n[i]-1]-coord[i][0]);
		double tmp_d = (length/(double)(n[i]-1));
		double tmp_eta;
		d[i] = 0.0;
		for(int j=1;j<n[i]-1;j++){
			if (type[i] == 0) coord[i][j] = coord[i][j-1] + tmp_d;
			if (type[i] == 1) coord[i][j] = grid->topo->bounds[i][1]-length*(cos(PI*0.5*(j/(double)(n[i]-1))));
			if (type[i] == 2) {
				tmp_eta = (j/(double)(n[i]-1));
				coord[i][j] = grid->topo->bounds[i][0]+length*((exp(2*tmp_eta)-1)/(exp(2)-1));
			}
			d[i] = fmax(d[i],fabs(coord[i][j]-coord[i][j-1])); 
		}
		d[i] = fmax(d[i],fabs(coord[i][n[i]-2]-coord[i][n[i]-1])); 
	}

	grid->h = 0.0;
	for (int i=0;i<dimension;i++){
		grid->h += d[i]*d[i];
	}
	grid->h = sqrt(grid->h);

	return ierr;
}

int MetricCoefficients2D(Array2d *metrics, double **coord, ArrayInt2d *IsGlobalToGrid, IsRange *range, double *bounds, int dimension, void (*MetricCoefficientsFunc)(double *metricsAtPoint, double *bounds, double *lengths, double x, double y)) {
	//This is a metric coefficients computing shell
	//Note: Metric coefficients at each point excluding BC points are computed
	//
	//metrics.indices[(n[0]-2)*(n[1]-2)][5]
	//metrics.data[(n[0]-2)*(n[1]-2)][5]
	//range[0:1] = start to end global index of unknowns in this process
	//
	//bounds[0] - Lower bound of "x"
	//bounds[1] - Upper bound of "x"
	//bounds[2] - Lower bound of "y"
	//bounds[3] - Upper bound of "y"
	//
	//(*MetricCoefficientsFunc) compute and return metric coefficients at (x,y)
	
	double	lengths[dimension];
	int	igrid, jgrid;
	int	ierr = 0;
	int	range0[2];

	range0[0] = range[0].start;	
	range0[1] = range[0].end;	
	metrics->ni = range0[1]-range0[0];
	metrics->nj = 5;
	metrics->data = malloc(metrics->ni*metrics->nj*sizeof(double));if (metrics->data==NULL) ERROR_MSG("malloc failed");

	for(int i=0;i<dimension;i++) {
		lengths[i] = bounds[i*2+1] - bounds[i*2];
	}
	
	for(int i=range0[0];i<range0[1];i++) {
		igrid = isGLOBALtoGRID(0,i,0);
		jgrid = isGLOBALtoGRID(0,i,1);
		(*MetricCoefficientsFunc)(((metrics->data)+((i-range0[0])*metrics->nj)),bounds,lengths,coord[0][jgrid+1],coord[1][igrid+1]);
	}
	
	return ierr;
}

static void CheckAndAssign(int *ijk, int *n, int *commcost, double *maxload, double *loadfac, int *nInterfaces, int *l) {
	// - Check if given factorization (ijk) is better based on provided communication
	// cost, maximum load and load factor
	// - If yes, then that factorization is assigned in ascending order and 
	// returned in "l" along with corresponding communication cost, maximum load and
	// load factor	
	
	int ijks[3], tcost, tInterfaces;
	double load[3];
		
	sort3(ijk, ijks);
	tcost =  CommCost3D(ijks, n);
	ComputeLoad(n, ijks, load);
	tInterfaces = ComputeNInterfaces3D(ijks);
	
	if (tcost < *commcost ||
	   (tcost == *commcost &&
	   (load[0] < *maxload ||
	   (load[0] == *maxload &&
	   (load[2] < *loadfac ||
	   (load[2] == *loadfac &&
	    tInterfaces < *nInterfaces)))))) {
		l[0] = ijks[0];
		l[1] = ijks[1];
		l[2] = ijks[2];
		*commcost = tcost;
		*maxload = load[0];
		*loadfac = load[2];
		*nInterfaces = tInterfaces;
	}
}

static void factorize2(int p, int i, int *nsorted, int *l, double *para) {
	// Factorize "p*i" into three factors (with "i" being one of them)
	// such that communication cost, maximum load and load factors
	// are reduced in that order of priority.
		
	double	temp = ((double)p*nsorted[1])/((double)nsorted[2]); 
	int	mopt = (int) floor(sqrt(temp));
	int	sqrtm = (int) floor(sqrt((double)p));
	int	k;
	int	cost = (int) para[0];
	double	maxload = para[1];
	double	loadfac = para[2];
	int	nInterfaces = (int) para[3];

	for (int j=mopt; j>0; j--) {
		k = p/j;
		if (k*j == p) {
			int ijk[3] = {i, j, k};
			CheckAndAssign(ijk, nsorted, &cost, &maxload, &loadfac, &nInterfaces, l);
			break;
		}
	}
	for (int j=mopt+1; j<sqrtm+1; j++) {
		k = p/j;
		if (k*j == p) {
			int ijk[3] = {i, j, k};
			CheckAndAssign(ijk, nsorted, &cost, &maxload, &loadfac, &nInterfaces, l);
			break;
		}
	}
	para[0] = (double) cost;
	para[1] = maxload;
	para[2] = loadfac;
	para[3] = (double) nInterfaces;
}

int factorize3(int procs, int *nsorted, int *l, double *para) {
	// Factorize "procs" into three factors such that 
	// communication cost, maximum load and load factors
	// are reduced in that order of priority.
	
	for (int i=0; i<3; i++)	l[i] = procs;
	para[0] = (double) CommCost3D(l, nsorted); // intialize with maximum cost
	para[1] = (double) nsorted[0]*nsorted[1]*nsorted[2];
	para[2] = (double) para[1];
	para[3] = (double) ComputeNInterfaces3D(l);

	double temp = (double)(procs*nsorted[0]*nsorted[0])/(double)(nsorted[1]*nsorted[2]);
	int lopt = (int) floor(pow(temp,(1.0/3.0)));
	int cubethp = (int) floor(pow(procs,(1.0/3.0)));
	int fac;

	for (int i=lopt; i>0; i--) {
		fac = procs/i;
		if (fac*i == procs) {
			factorize2(fac, i, nsorted, l, para);
			break;
		}
	}
	for (int i=lopt+1; i<cubethp+1; i++) {
		fac = procs/i;
		if (fac*i == procs) {
			factorize2(fac, i, nsorted, l, para);
			break;
		}
	}

	return 0;
}

int factorize(Grid *grid) {
	// Factorize "procs" depending on the dimension such that 
	// communication cost, maximum load and load factors
	// are reduced in that order of priority.
	
	int	procs, rank;
	int	ierr = 0;
	
	MPI_Comm_size(MPI_COMM_WORLD, &procs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	int dimension = grid->topo->dimension;
	int *dimProcs = grid->topo->dimProcs;
	int *n = grid->n;
	double *para = grid->para;; 
	
	int nbc = 2; // Assuming 2 boundary points; 1 on each side
	int nsorted[MAX_DIMENSION], index[3];

	if (dimension == 2) n[2] = nbc+1;
	sort3Index(n, nsorted, index); 
	for (int i=0; i<3; i++) nsorted[i] = nsorted[i]-nbc;

	int l[3];
	if (dimension == 2) {
		l[0] = 1;
		for (int i=1; i<3; i++)	l[i] = procs;
		para[0] = (double) CommCost3D(l, nsorted); // intialize with maximum cost
		para[1] = (double) nsorted[0]*nsorted[1]*nsorted[2];
		para[2] = (double) para[1];
		para[3] = (double) ComputeNInterfaces3D(l);
		factorize2(procs, 1, nsorted, l, para);
	} else if (dimension == 3) {
		ierr = factorize3(procs, nsorted, l, para);
	}

	if (l[0] > nsorted[0] || l[1] > nsorted[1] || l[2] > nsorted[2]) {
		pERROR_MSG("Rerun with different no. of procs or mesh sizes");
		return 1;
	};
	for (int i=0; i<3; i++)	dimProcs[i] = l[index[i]];

//	PetscPrintf(PETSC_COMM_WORLD,"Chosen: \n");
//	for (int i=0; i<dimension; i++)
//		PetscPrintf(PETSC_COMM_WORLD,"p[%d] = %d, ", i, dimProcs[i]);
//	PetscPrintf(PETSC_COMM_WORLD,"CommCost = %d, MaxLoad = %d, LoadFactor = %lf, nInterfaces = %d\n", (int)para[0], (int)para[1], para[2], (int)para[3]);
	
	return ierr;
}

int split_domain(Grid *grid) {
	// Split the domain for given factorization of total
	// no. of procs.
	
	int	procs, rank;
	
	MPI_Comm_size(MPI_COMM_WORLD, &procs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	int dimension = grid->topo->dimension;
	int *l = grid->topo->dimProcs;
	int *n = grid->n;
	int (*range)[2] = grid->range;
	int *blockID = grid->topo->blockID;
	int nbc = 2; // Assuming 2 boundary points; 1 on each side

	blockID[2] = rank/(l[0]*l[1]);
	blockID[1] = (rank-blockID[2]*l[0]*l[1])/l[0];
	blockID[0] = (rank-blockID[2]*l[0]*l[1]-blockID[1]*l[0]);
	
	int rem[3], quo[3];
	for (int i=0; i<dimension; i++) {
		rem[i] = (n[i]-nbc)%l[i];
		quo[i] = (n[i]-nbc)/l[i];
	}

	for (int i=0; i<dimension; i++) {
		for (int j=0; j<2; j++)
			range[i][j] = 1 + quo[i]*(blockID[i]+j) + ((blockID[i]+j < rem[i]) ? blockID[i]+j : rem[i]) ;
	}
	for (int i=0; i<dimension; i++) {
		if (blockID[i] == 0) range[i][0] = range[i][0] - 1;
		if (blockID[i] == l[i]-1) range[i][1] = range[i][1] + 1;
	}
	
	return 0;
}

int ReadTopo(Topo *topo) {
	// Reads topological information from input
	
	int		ierr = 0;
	PetscBool	set;	

	PetscOptionsGetInt(NULL, NULL, "-dim", &(topo->dimension), &set);
	if (!set) {
		pERROR_MSG("Dimension of the problem not set");
		pERROR_MSG("Set '-dim n' for n-dimension");
		return 1;
	} else if (topo->dimension > 3 || topo->dimension < 2) {
		pERROR_MSG("Only 2D or 3D is valid");
		return 1;
	}
	
	int		nmax;
	nmax = topo->dimension;
	
	PetscOptionsGetIntArray(NULL, NULL, "-grid_type", topo->gridtype, &nmax, &set);
	if (!set || nmax != topo->dimension) {
		pERROR_MSG("Grid type is not set properly");
		pERROR_MSG("Set '-grid_type n1,n2,n3' for n1, n2 and n3 type spacing in each direction");
		return 1;
	}
	nmax = topo->dimension*2;
	double	temp[MAX_DIMENSION*2];
	PetscOptionsGetRealArray(NULL, NULL, "-bounds", temp, &nmax, &set);
	if (!set || nmax != topo->dimension*2) {
		pERROR_MSG("Grid bounds are not set properly");
		pERROR_MSG("Set '-bounds l0,u0,l1,u1,l2,u2' for lower and upper bounds in each of the three dimensions");
		return 1;
	}
	for (int i=0; i<topo->dimension; i++) {
		for (int j=0; j<2; j++) {
			topo->bounds[i][j] = temp[i*2+j];
		}
	}

	return ierr;
}

int create_coarse_grid(Grid *topgrid, Grid *botgrid, int *cfactor) {
	// Creates botgrid by coarsening topgrid using coarsening factors (cfactor)

	int	ierr=0;

	int	dimension = topgrid->topo->dimension;
	char	dir[4] = "ijk";

	botgrid->topo = topgrid->topo;

	for (int i=0;i<dimension;i++) {
		int temp;
		temp = (topgrid->n[i]-1)/cfactor[i];
		if (temp*cfactor[i] != topgrid->n[i]-1 || temp==0) {
			PetscPrintf(PETSC_COMM_WORLD,"ERROR: %s:%d: Coarsening along %.1s-th direction failed\n",__FILE__,__LINE__,dir+i);
			PetscPrintf(PETSC_COMM_WORLD,"ERROR: %s:%d: Change no. of points along %.1s-th direction or no. of grids\n",__FILE__,__LINE__,dir+i);
			return 1;
		}
		botgrid->n[i] = temp+1;

		for (int j=0;j<2;j++)
			botgrid->range[i][j] = (topgrid->range[i][j]+cfactor[i]-1)/cfactor[i];
	}

	ierr = malloc2dY(&(botgrid->coord), dimension, botgrid->n);
	if (ierr) {
		pERROR_MSG("Mesh memory allocation failed");
		return 1;
	}
	
	for (int i=0;i<dimension;i++) {
		for (int j=0;j<botgrid->n[i];j++) {
			botgrid->coord[i][j] = topgrid->coord[i][cfactor[i]*j];
		}
	}
	
	double	d[MAX_DIMENSION];
	for (int i=0;i<dimension;i++) {
		d[i] = 0.0;
		for (int j=1;j<botgrid->n[i];j++) {
			d[i] = fmax(d[i],fabs(botgrid->coord[i][j]-botgrid->coord[i][j-1])); 
		}
	}

	botgrid->h = 0.0;
	for (int i=0;i<dimension;i++){
		botgrid->h += d[i]*d[i];
	}
	botgrid->h = sqrt(botgrid->h);
	
	double	load[3];
	int	nreduced[3], nbc = 2; // Default no. of BC points in a direction

	if (dimension == 2) botgrid->n[2] = 1 + nbc;
	for (int i=0;i<3;i++) nreduced[i] = botgrid->n[i]-nbc;
	ComputeLoad(nreduced, botgrid->topo->dimProcs, load);

	botgrid->para[0] = (double) CommCost3D(botgrid->topo->dimProcs, nreduced);
	botgrid->para[1] = load[0];
	botgrid->para[2] = load[2];
	botgrid->para[3] = topgrid->para[3];
	return 0;
}

int create_coarse_grids(Grids *grids) {
	
	int	ierr=0;
	int	ngrids = grids->ngrids; 
	
	for (int i=0; i<ngrids-1; i++) {
		Grid	*topgrid = grids->grid+i;
		Grid	*botgrid = grids->grid+i+1;
		int	*cfactor = grids->cfactor[i];
		botgrid->id = i+1;
		ierr = create_coarse_grid(topgrid, botgrid, cfactor);
		if (ierr) {
			PetscPrintf(PETSC_COMM_WORLD,"ERROR: %s:%d: Coarsening of grid-%d failed\n",__FILE__,__LINE__,i+1);
			return 1;
		}
	}

	return 0;	
}

int CreateGrids(Grids *grids) {
	// Creates all grids
	// =================
	//
	// Reads information
	// Allocates memory
	
	int ierr = 0;
	
	grids->topo = malloc(sizeof(Topo));
	Topo	*topo = grids->topo;
	
	ierr = ReadTopo(topo); if(ierr) return 1;

	PetscBool	set;	
	PetscOptionsGetInt(NULL, NULL, "-ngrids", &(grids->ngrids), &set);
	if (!set) {
		pERROR_MSG("Number of grids for MG not set");
		pERROR_MSG("Set '-ngrids n' for n grids");
		return 1;
	}
	grids->grid = malloc(grids->ngrids*sizeof(Grid));

	int	nmax;
	Grid	*grid = grids->grid;

	grid->topo = grids->topo;
	grid->id = 0;
	nmax = topo->dimension;

	PetscOptionsGetIntArray(NULL, NULL, "-npts", grid->n, &nmax, &set);
	if (!set || nmax != topo->dimension) {
		pERROR_MSG("No. of grid points not set properly");
		pERROR_MSG("Set '-npts n0,n1,n2' for no. of grid points in each of the three dimensions");
		return 1;
	}
	ierr = malloc2dY(&(grid->coord), topo->dimension, grid->n);
	if (ierr) {
		pERROR_MSG("Mesh memory allocation failed");
		return 1;
	}

	for (int i=0;i<MAX_GRIDS;i++)
		for (int j=0;j<MAX_DIMENSION;j++)
			grids->cfactor[i][j] = 2; // Default coarsening factor "2"

	ierr = Compute_coord(grid); pCHKERR_RETURN("Coordinate computation failed");
	ierr = factorize(grid); pCHKERR_RETURN("Factorization of total no. of procs failed");
	ierr = split_domain(grid); pCHKERR_RETURN("Domain splitting failed");
	ierr = create_coarse_grids(grids); pCHKERR_RETURN("Coarse grids creation failed");
//	if (mesh->type == UNIFORM) mesh->MetricCoefficients = &MetricsUniform;
//	if (mesh->type == NONUNIFORM1) mesh->MetricCoefficients = &MetricsNonUniform1;
//	if (mesh->type == NONUNIFORM2) mesh->MetricCoefficients = &MetricsNonUniform2;

	return ierr;
}

void DestroyGrids(Grids *grids) {
	// Deallocates memory
	
	for (int i=0;i<grids->ngrids;i++)
		if (grids->grid[i].coord) free2dArray(&(grids->grid[i].coord));
	if (grids->grid) free(grids->grid);
	if (grids->topo) free(grids->topo);

}
