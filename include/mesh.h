



#ifndef _MESH_H_
#define _MESH_H_

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <petscksp.h>

#include "array.h"
#include "problem.h"

#define PI 3.14159265358979323846
#define MAX_DIMENSION 3
#define MIN_POINTS 3
#define	MAX_GRIDS 40

//typedef enum {UNIFORM, NONUNIFORM1, NONUNIFORM2} MeshType;

typedef struct {
	int	dimension;
	int	gridtype[MAX_DIMENSION];
	int	blockID[MAX_DIMENSION];
	int	dimProcs[MAX_DIMENSION];
	double	bounds[MAX_DIMENSION][2];
} Topo;

typedef struct {
	int	dimension;
	int	gridtype[MAX_DIMENSION];
	int	blockID[MAX_DIMENSION];
	int	dimProcs[MAX_DIMENSION];
	double	bounds[MAX_DIMENSION][2];
	int	n[MAX_DIMENSION]; // Grids in each direction
	int	cfactor[MAX_DIMENSION]; // Coarsening factor for next grid
	double	range[MAX_DIMENSION][2]; // Range of grid points
	double	**coord; // Coordinates in each direction
	double	h; // Grid characteristic length
	double	para[4]; // Domain splitting quality measures
	void	(*MetricCoefficients)(void *mesh, double x, double y, double *metrics);
} Mesh;

typedef struct {
	Topo	*topo; //Topology info
	int	id; // grid id; 0:finest, ngrids-1:coarsest
	int	n[MAX_DIMENSION]; // No. of grid points in each direction
	int	ln[MAX_DIMENSION]; // No. of local points in each direction
	int	**range; // Range of grid points in rank and direction wise
	double	**coord; // Coordinates in each direction
	double	h; // Grid characteristic length
	double	para[4]; // Domain splitting quality measures
	void	(*MetricCoefficients)(void *mesh, double x, double y, double *metrics);
} Grid;

typedef struct {
	int	ngrids; // Total number of grids
	int	cfactor[MAX_GRIDS][MAX_DIMENSION]; // Coarsening factors for all grids
	Topo	*topo; // Topology info
	Grid	*grid; // Store each grid
} Grids;

int CreateGrids(Grids *grids);
void DestroyGrids(Grids *grids);

//int CreateMesh(Mesh *mesh);
//void DestroyMesh(Mesh *mesh);

#endif
