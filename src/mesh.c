#include "mesh.h"

#define ERROR_MSG(message) (fprintf(stderr,"Error:%s:%d: %s\n",__FILE__,__LINE__,(message)))
#define ERROR_RETURN(message) {ERROR_MSG(message);return ierr;}
#define CHKERR_PRNT(message) {if(ierr==1) {ERROR_MSG(message);}}
#define CHKERR_RETURN(message) {if(ierr==1) {ERROR_RETURN(message);}}
#define PI 3.14159265358979323846

#define METRICS(i,j,k) (metrics->data[metrics->nk*((i)*metrics->nj+(j))+(k)])
#define isGRIDtoGLOBAL(l,i,j) (IsGridToGlobal[l].data[((i)*IsGridToGlobal[l].nj+(j))])
#define isGLOBALtoGRID(l,i,j) (IsGlobalToGrid[l].data[((i)*IsGlobalToGrid[l].nj+(j))])
//#define isGRIDtoGLOBAL(i,j) (IsGridToGlobal.data[((i)*IsGridToGlobal.nj+(j))])
//#define isGLOBALtoGRID(i,j) (IsGlobalToGrid.data[((i)*IsGlobalToGrid.nj+(j))])

double TransformFunc(double *bounds, double length, double xi) {
	//Transformation function from computational to physical space
	//
	//bounds - lower and upper bounds of physical coordinate
	//length = (bounds[1]-bounds[0])
	//
	//x or y = T(xi)
	
	double val;
//	val = bounds[1]-length*(cos(PI*0.5*xi));
	val = xi;
	return val;
}

void MetricsUniform(void *mesh, double x, double y, double *metrics) {
	//Computes following metrics at (x,y)
	//
	//metrics[0] = (xi_x)^2 + (xi_y)^2
	//metrics[1] = (eta_x)^2 + (eta_y)^2
	//metrics[2] = (xi_xx) + (xi_yy)
	//metrics[3] = (eta_xx) + (eta_yy)
	//metrics[4] = (xi_x)(eta_x) + (xi_y)(eta_y)
	
	metrics[0] = 1.0;
	metrics[1] = 1.0;
	metrics[2] = 0.0;
	metrics[3] = 0.0; 
	metrics[4] = 0.0;
}

void MetricsNonUniform1(void *mesh1, double x, double y, double *metrics) {
	//Computes following metrics at (x,y)
	//
	//metrics[0] = (xi_x)^2 + (xi_y)^2
	//metrics[1] = (eta_x)^2 + (eta_y)^2
	//metrics[2] = (xi_xx) + (xi_yy)
	//metrics[3] = (eta_xx) + (eta_yy)
	//metrics[4] = (xi_x)(eta_x) + (xi_y)(eta_y)
	//
	//bounds[0] - Lower bound of "x"
	//bounds[1] - Upper bound of "x"
	//bounds[2] - Lower bound of "y"
	//bounds[3] - Upper bound of "y"
	
	double	temp;
	double	*bounds;
	Mesh	*mesh;
	int	procs, rank;
	
	MPI_Comm_size(PETSC_COMM_WORLD, &procs);
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	
	mesh = (Mesh*)mesh1;
	bounds = mesh->bounds;
	temp = ((bounds[3]-bounds[2])*(bounds[3]-bounds[2])-(bounds[3]-y)*(bounds[3]-y));
	metrics[0] = 1.0;
	metrics[1] = 4.0/(PI*PI*temp);
	metrics[2] = 0.0;
	metrics[3] = (-2.0*(bounds[3]-y))/(PI*sqrt(temp*temp*temp)); 
	metrics[4] = 0.0;
}

void MetricsNonUniform2(void *mesh1, double x, double y, double *metrics) {
	//Computes following metrics at (x,y)
	//
	//metrics[0] = (xi_x)^2 + (xi_y)^2
	//metrics[1] = (eta_x)^2 + (eta_y)^2
	//metrics[2] = (xi_xx) + (xi_yy)
	//metrics[3] = (eta_xx) + (eta_yy)
	//metrics[4] = (xi_x)(eta_x) + (xi_y)(eta_y)
	//
	//bounds[0] - Lower bound of "x"
	//bounds[1] - Upper bound of "x"
	//bounds[2] - Lower bound of "y"
	//bounds[3] - Upper bound of "y"
	
	double	temp;
	double	*bounds;
	Mesh	*mesh;
	int	procs, rank;
	
	MPI_Comm_size(PETSC_COMM_WORLD, &procs);
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	
	mesh = (Mesh*)mesh1;
	bounds = mesh->bounds;
	temp = ((exp(2)-1)*(exp(2)-1))/(((y-bounds[2])*(exp(2)-1)+(bounds[3]-bounds[2]))*((y-bounds[2])*(exp(2)-1)+(bounds[3]-bounds[2])));
	metrics[0] = 1.0/((bounds[1]-bounds[0])*(bounds[1]-bounds[0]));
	metrics[1] = 0.25*temp;
	metrics[2] = 0.0;
	metrics[3] = (-0.5)*temp; 
	metrics[4] = 0.0;
}

int UniformMesh(double ***pcoord, int *n, double *bounds, double *h, int dimension) {
	
	int ierr = 0;
	
	//Assign memory
	ierr = malloc2dY(pcoord,dimension,n); CHKERR_RETURN("malloc failed");
	//Compute uniform grid in each dimension
	for(int i=0;i<dimension;i++){
		if (n[i]<2) {ierr=1; ERROR_RETURN("Need at least 2 points in each direction");}
		(*pcoord)[i][0] = bounds[i*2]; //Lower bound
		(*pcoord)[i][n[i]-1] = bounds[i*2+1]; //Upper bound
		
		h[i] = ((*pcoord)[i][n[i]-1]-(*pcoord)[i][0])/(n[i]-1); //Spacing
		for(int j=1;j<n[i]-1;j++){
			(*pcoord)[i][j] = (*pcoord)[i][j-1] + h[i];
		}
	}

	return ierr;
}

int Coords(Mesh *mesh, MeshType type) {
	// computes coords in each direction of the structured grid
	//
	// MeshType type: {UNIFORM, NONUNIFORM}
	
	int	ierr = 0;
	double	length, d[MAX_DIMENSION];

	int	*n;
	double	**coord, *bounds;
	
	n 	= mesh->n;
	bounds	= mesh->bounds;
	coord	= mesh->coord;
	for(int i=0;i<1;i++){
		if (n[i]<2) {ierr=1; ERROR_RETURN("Need at least 2 points in each direction");}
		coord[i][0] = bounds[i*2]; //Lower bound
		coord[i][n[i]-1] = bounds[i*2+1]; //Upper bound
		
		d[i] = (coord[i][n[i]-1]-coord[i][0])/(n[i]-1); //Spacing
		for(int j=1;j<n[i]-1;j++){
			coord[i][j] = coord[i][j-1] + d[i];
		}
	}

	for(int i=1;i<2;i++){
		if (n[i]<2) {ierr=1; ERROR_RETURN("Need at least 2 points in each direction");}
		coord[i][0] = bounds[i*2]; //Lower bound
		coord[i][n[i]-1] = bounds[i*2+1]; //Upper bound
		
		length = (coord[i][n[i]-1]-coord[i][0]);
		double tmp_d = (length/(double)(n[i]-1));
		double tmp_eta;
		d[i] = 0.0;
		for(int j=1;j<n[i]-1;j++){
			if (type == NONUNIFORM1) coord[i][j] = bounds[i*2+1]-length*(cos(PI*0.5*(j/(double)(n[i]-1))));
			if (type == NONUNIFORM2) {
				tmp_eta = (j/(double)(n[i]-1));
				coord[i][j] = bounds[i*2]+length*((exp(2*tmp_eta)-1)/(exp(2)-1));
			}
			if (type == UNIFORM) coord[i][j] = coord[i][j-1] + tmp_d;
//			coord[i][j] = (*Transform)(&(bounds[i*2]), length, j/(double)(n[i]-1));
			d[i] = fmax(d[i],fabs(coord[i][j]-coord[i][j-1])); 
		}
		d[i] = fmax(d[i],fabs(coord[i][n[i]-2]-coord[i][n[i]-1])); 
	}

	for(int i=2;i<DIMENSION;i++){
		if (n[i]<2) {ierr=1; ERROR_RETURN("Need at least 2 points in each direction");}
		coord[i][0] = bounds[i*2]; //Lower bound
		coord[i][n[i]-1] = bounds[i*2+1]; //Upper bound
		
		d[i] = (coord[i][n[i]-1]-coord[i][0])/(n[i]-1); //Spacing
		for(int j=1;j<n[i]-1;j++){
			coord[i][j] = coord[i][j-1] + d[i];
		}
	}
	
	mesh->h = 0.0;
	for (int i=0;i<DIMENSION;i++){
		mesh->h += d[i]*d[i];
	}
	mesh->h = sqrt(mesh->h);

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

void CreateMesh(Mesh *mesh, int meshflag) {
	// Creates mesh 
	// Allocates memory
	// Computes coords of a structured mesh
	// Assigns metric coefficients computing function
	int ierr = 0;

	// Reads mesh inputs
	PetscBool	set;	
	PetscOptionsGetInt(NULL, NULL, "-dim", &(mesh->dimension), &set);
	if (!set) {
		PetscPrintf(PETSC_COMM_WORLD, "ERROR: Dimension of the problem not set!\n"); 
		PetscFinalize();
		MPI_Finalize();
		return 0;
	}
	int		nmax;
	nmax = mesh->dimension;
	PetscOptionsGetIntArray(NULL, NULL, "-npts", mesh->n, &nmax, &set);
	if (!set || nmax != mesh->dimension) {
		PetscPrintf(PETSC_COMM_WORLD, "ERROR: No. of mesh points not set properly!\n"); 
		PetscFinalize();
		MPI_Finalize();
		return 0;
	}
	int	meshflag;
	PetscOptionsGetInt(NULL, NULL, "-mesh", &meshflag, &set);
	if (!set) {
		PetscPrintf(PETSC_COMM_WORLD, "ERROR: Mesh type is not set properly!\n"); 
		PetscFinalize();
		MPI_Finalize();
		return 0;
	}
	nmax = mesh->dimension*2;
	PetscOptionsGetIntArray(NULL, NULL, "-bounds", mesh->bounds, &nmax, &set);
	if (!set || nmax != mesh->dimension*2) {
		PetscPrintf(PETSC_COMM_WORLD, "ERROR: Mesh bounds are not set properly!\n"); 
		PetscFinalize();
		MPI_Finalize();
		return 0;
	}
	ierr = malloc2dY(&(mesh->coord), mesh.dimension, mesh->n); CHKERR_PRNT("malloc failed");
	if (meshflag == 0) mesh->type = UNIFORM;
	if (meshflag == 1) mesh->type = NONUNIFORM1;
	if (meshflag == 2) mesh->type = NONUNIFORM2;
	ierr = Coords(mesh, type); CHKERR_PRNT("Meshing failed");
	if (mesh->type == UNIFORM) mesh->MetricCoefficients = &MetricsUniform;
	if (mesh->type == NONUNIFORM1) mesh->MetricCoefficients = &MetricsNonUniform1;
	if (mesh->type == NONUNIFORM2) mesh->MetricCoefficients = &MetricsNonUniform2;

//	return ierr;
}

void DestroyMesh(Mesh *mesh) {
	// Deallocates memory
	
	if (mesh->coord != NULL) free2dArray(&(mesh->coord));

}
