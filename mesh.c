#include "mesh.h"

#define ERROR_MSG(message) (fprintf(stderr,"Error:%s:%d: %s\n",__FILE__,__LINE__,(message)))
#define ERROR_RETURN(message) {ERROR_MSG(message);return ierr;}
#define CHKERR_PRNT(message) {if(ierr==1) {ERROR_MSG(message);}}
#define CHKERR_RETURN(message) {if(ierr==1) {ERROR_RETURN(message);}}
#define PI 3.14159265358979323846

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

int NonUniformMeshY(double ***pcoord, int *n, double *bounds, double *h, int dimension, double (*Transform)(double *bounds, double range, double s) ) {
	
	int	ierr = 0;
	double	range, d[dimension];
	
	//Assign memory
	ierr = malloc2dY(pcoord,dimension,n); CHKERR_RETURN("malloc failed");
	//Compute uniform grid in each dimension
	for(int i=0;i<1;i++){
		if (n[i]<2) {ierr=1; ERROR_RETURN("Need at least 2 points in each direction");}
		(*pcoord)[i][0] = bounds[i*2]; //Lower bound
		(*pcoord)[i][n[i]-1] = bounds[i*2+1]; //Upper bound
		
		d[i] = ((*pcoord)[i][n[i]-1]-(*pcoord)[i][0])/(n[i]-1); //Spacing
		for(int j=1;j<n[i]-1;j++){
			(*pcoord)[i][j] = (*pcoord)[i][j-1] + d[i];
		}
	}

	for(int i=1;i<2;i++){
		if (n[i]<2) {ierr=1; ERROR_RETURN("Need at least 2 points in each direction");}
		(*pcoord)[i][0] = bounds[i*2]; //Lower bound
		(*pcoord)[i][n[i]-1] = bounds[i*2+1]; //Upper bound
		
		range = ((*pcoord)[i][n[i]-1]-(*pcoord)[i][0]);
		d[i] = 0.0;
		for(int j=1;j<n[i]-1;j++){
			(*pcoord)[i][j] = (*Transform)(&(bounds[i*2]), range, j/(double)(n[i]-1));
			d[i] = fmax(d[i],fabs((*pcoord)[i][j]-(*pcoord)[i][j-1])); 
		}
		d[i] = fmax(d[i],fabs((*pcoord)[i][n[i]-2]-(*pcoord)[i][n[i]-1])); 
	}

	for(int i=2;i<dimension;i++){
		if (n[i]<2) {ierr=1; ERROR_RETURN("Need at least 2 points in each direction");}
		(*pcoord)[i][0] = bounds[i*2]; //Lower bound
		(*pcoord)[i][n[i]-1] = bounds[i*2+1]; //Upper bound
		
		d[i] = ((*pcoord)[i][n[i]-1]-(*pcoord)[i][0])/(n[i]-1); //Spacing
		for(int j=1;j<n[i]-1;j++){
			(*pcoord)[i][j] = (*pcoord)[i][j-1] + d[i];
		}
	}
	
	*h = 0.0;
	for (int i=0;i<dimension;i++){
		*h = *h + d[i]*d[i];
	}
	*h = sqrt(*h);

	return ierr;
}
