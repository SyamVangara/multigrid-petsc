//#include "matbuild.h"
#include "solver.h"

#define METRICS(i,j) (metrics.data[(i)*metrics.nj+(j)])
//#define METRICS(i,j,k) (metrics.data[metrics.nk*((i)*metrics.nj+(j))+(k)])
#define PMETRICS(i) ((metrics.data)+((i)*metrics.nj))
#define F(i,j) (f.data[((i)*f.nj+(j))])
#define U(i,j) (u.data[((i)*u.nj+(j))])

static int ipow(int base, int exp) {
	// Computes integer power of a real or integer
	//
	// base - real or int
	// exp - int
	// output: base^(exp)

	int result = 1;
	while (exp) {
		if (exp & 1)
			result *= base;
		exp >>= 1;
		base *= base;
	}
	return result;
}

void GridId(Indices *indices) {
	// Computes the gridId range and number grids in each level
	
	int q, count;
	
	q = 0;
	for (int l=0;l<indices->levels;l++) {
		indices->level[l].grids = 1;
		q += 1;
	}
	indices->level[indices->levels-1].grids += (indices->totalGrids-q);
	
	count = 0;
	for (int l=0;l<indices->levels;l++) {
		indices->level[l].gridId = malloc(indices->level[l].grids*sizeof(int));
		for (int g=0;g<indices->level[l].grids;g++) {
			indices->level[l].gridId[g] = count;
			count += 1;
		}
	}
}

void CreateIndexMaps(Mesh *mesh, Indices *indices) {
	// Allocate memory to global-to-grid and grid-to-global maps
	
	int	temp, g;
	int	totaln, n0, n1, fn0, fn1;
	int	factor;
	const	int	M = 3; // i,j,g
	
	factor = indices->coarseningFactor;
	fn0 = mesh->n[0];
	fn1 = mesh->n[1];
	for (int l=0;l<indices->levels;l++) {
		totaln = 0;
		for (int lg=0;lg<indices->level[l].grids;lg++) {
			g = indices->level[l].gridId[lg];
			temp = ipow(factor,g);
			n0 = (fn0-1)/temp - 1;
			n1 = (fn1-1)/temp - 1;
			totaln = totaln + n0*n1; 
			CreateArrayInt2d(n1, n0,&(indices->level[l].grid[lg]));
		}
		CreateArrayInt2d(totaln, M, &(indices->level[l].global));
	}
}

void DestroyIndexMaps(Indices *indices) {
	// Free the memory of global-to-grid and grid-to-global maps	
	
	for (int l=0;l<indices->levels;l++) {
		for (int g=0;g<indices->level[l].grids;g++) {
			DeleteArrayInt2d(&(indices->level[l].grid[g]));
		}
		DeleteArrayInt2d(&(indices->level[l].global));
	}
}

void SetUpIndices(Mesh *mesh, Indices *indices) {
	// Get the GridId range, then allocate memory
	
	int	procs;
	MPI_Comm_size(PETSC_COMM_WORLD, &procs);
	
	indices->level = malloc(indices->levels*sizeof(Level));
	GridId(indices);	
	for (int i=0;i<indices->levels;i++) {
		indices->level[i].ranges = malloc((procs+1)*sizeof(int));
		indices->level[i].grid = malloc(indices->level[i].grids*sizeof(ArrayInt2d));
		indices->level[i].h = malloc(indices->level[i].grids*sizeof(double[2]));
	}
	CreateIndexMaps(mesh, indices);
	for (int l=0;l<indices->levels;l++) {
		for (int lg=0;lg<indices->level[l].grids;lg++) {
			indices->level[l].h[lg][0] = 1.0/(indices->level[l].grid[lg].ni+1);
			indices->level[l].h[lg][1] = 1.0/(indices->level[l].grid[lg].nj+1);
		}
	}
}

void DestroyIndices(Indices *indices) {
	// Free the memory allocated to indices
	
	DestroyIndexMaps(indices);
	for (int l=0;l<indices->levels;l++) {
		free(indices->level[l].grid);
		free(indices->level[l].h);
		free(indices->level[l].gridId);
		free(indices->level[l].ranges);
	}
	free(indices->level);
}

void GetRanges(int *ranges, int totaln) {
	// Computes the ranges of indices for all processes for a given total number of indices	
	//
	// length of ranges = procs + 1
	// range[p] - index start in process "p"
	// range[p+1] - (index end + 1) in process "p"
	
	int	remainder, quotient;
	int	procs, rank;
	
	MPI_Comm_size(PETSC_COMM_WORLD, &procs);
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	
	remainder = (totaln)%procs;
	quotient  = (totaln)/procs;
	ranges[0] = 0;
	for (int p=0;p<procs;p++) {
		if (p<remainder) {
			ranges[p+1] = ranges[p] + (quotient + 1);
		}
		else {
			ranges[p+1] = ranges[p] + quotient;
		}
	}
}

void mappingThroughGrids(Indices *indices) {
	// Maps global indices to grid unknowns and vice-versa
	// Mapping style: places points of any grid corresponding to same (x, y) physical coordinates next to each other
	
	int	count, fcount;
	int	factor, gfactor;
	int	ig, jg;
	int	gfine, g;
	int	*ranges;
	ArrayInt2d	a, b, fine;
	
	int	procs, rank;
	
	MPI_Comm_size(PETSC_COMM_WORLD, &procs);
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	
	factor = indices->coarseningFactor;
	for (int l=0;l<indices->levels;l++) {
		count = 0; // counts number of points mapped
		a 	= indices->level[l].global;
		ranges 	= indices->level[l].ranges;
		gfine 	= indices->level[l].gridId[0];
		fine 	= indices->level[l].grid[0];

		// Compute the ranges of fine grid indices for processes
		GetRanges(ranges, fine.ni*fine.nj);
		fcount = 0; // counts number of fine grid points mapped
		rank = 0; // Keeps track of the rank
		for (int i=0;i<fine.ni;i++) {
			for (int j=0;j<fine.nj;j++) {
				for (int lg=0;lg<indices->level[l].grids;lg++) {
					g = indices->level[l].gridId[lg];
					gfactor = ipow(factor, g-gfine);
					if ((i+1)%gfactor!=0 || (j+1)%gfactor!=0) continue;
					ig = (i+1)/gfactor-1;
					jg = (j+1)/gfactor-1;
					b = indices->level[l].grid[lg];
					b.data[ig*b.nj+jg] = count;
					a.data[count*a.nj+0] = ig;
					a.data[count*a.nj+1] = jg;
					a.data[count*a.nj+2] = g;
					count = count + 1;
				}
				fcount += 1;
				if (fcount == ranges[rank+1]) {
					ranges[rank+1] = count; // update the number of points in "rank"
					rank += 1;
				}
			}
		}
	}
}

void mappingGridAfterGrid(Indices *indices) {
	// Maps global indices to grid unknowns and vice-versa
	// Mapping style: places all the points of a grid together
	
	int	count;
	int	g;
	int	*ranges;
	ArrayInt2d	a, b;
	
	for (int l=0;l<indices->levels;l++) {
		count = 0;
		a = indices->level[l].global;
		ranges 	= indices->level[l].ranges;
		for (int lg=0;lg<indices->level[l].grids;lg++) {
			g = indices->level[l].gridId[lg];
			b = indices->level[l].grid[lg];
			for (int i=0;i<b.ni;i++) {
				for (int j=0;j<b.nj;j++) {
					b.data[i*b.nj+j] = count;
					a.data[count*a.nj+0] = i;
					a.data[count*a.nj+1] = j;
					a.data[count*a.nj+2] = g;
					count = count + 1;
				}
			}
		}
		GetRanges(ranges, count);
	}

}

void mapping(Indices *indices, int mappingStyleflag) {
	// Maps global indices to grid unknowns and vice-versa
	
       	if (mappingStyleflag == 0) {
		mappingGridAfterGrid(indices);
	} else if (mappingStyleflag == 1) {
		mappingThroughGrids(indices);
	} else {
		printf("Unknown indices mapping style\n");
	}
}


void SetUpOperator(Indices *indices, Operator *op) {
	// Allocates memory to Operator struct
	int	order; // order of grid transfer operators
	int	stencilSize;
	
	order = indices->coarseningFactor; // order of grid transfer operators is same as the coarsening factor
	op->totalGrids = indices->totalGrids;
	op->res = malloc((op->totalGrids-1)*sizeof(ArrayInt2d));
	op->pro = malloc((op->totalGrids-1)*sizeof(ArrayInt2d));
	stencilSize = 1;
	for (int i=0;i<op->totalGrids-1;i++) {
		stencilSize = (stencilSize+1)*order-1;
		CreateArray2d(stencilSize, stencilSize, op->res+i);
		CreateArray2d(stencilSize, stencilSize, op->pro+i);
	}

}

void DestroyOperator(Operator *op) {
	// Free the memory in Operator struct
	
	for (int i=0;i<op->totalGrids-1;i++) {
		DeleteArray2d(op->res+i);
		DeleteArray2d(op->pro+i);
	}
	free(op->res);
	free(op->pro);
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

void ProlongationOperator(Array2d pro) {
	// Builds prolongation 2D stencilwise operator (pro)
	// Stencil size: pro.ni x pro.nj
	if(pro.ni != 3 || pro.nj != 3) {printf("Error in ProlongationOperator\n"); return;}
	for (int i=0;i<pro.ni;i++) {
 		pro.data[i*pro.nj]= 0.5 - 0.25*fabs(1-i);
 		pro.data[i*pro.nj+1]= 1.0 - 0.5*fabs(1-i);
 		pro.data[i*pro.nj+2]= 0.5 - 0.25*fabs(1-i);
	}
}

void RestrictionOperator(Array2d res) {
	// Builds Restriction 2D stencilwise operator (res)
	// Stencil size: res.ni x res.nj
	
	if(res.ni != 3 || res.nj != 3) {printf("Error in RestrictionOperator\n"); return;}
	for (int i=0;i<res.nj;i++) {
 		res.data[i*res.nj]= 0.0;
 		res.data[i*res.nj+1]= 0.0;
 		res.data[i*res.nj+2]= 0.0;
	}
	res.data[res.nj+1] = 1.0;
}

void GridTransferOperators(Operator op, Indices indices) {
	// Builds stencilwise grid transfer operators between any two grids that have "x" grids inbetween
	// where x = {0,...,levels-2}
	if (op.totalGrids < 2) return;
	RestrictionOperator(op.res[0]);
	ProlongationOperator(op.pro[0]);
	
	GridTransferOperator(op.res, indices.coarseningFactor, op.totalGrids);
	GridTransferOperator(op.pro, indices.coarseningFactor, op.totalGrids);
}


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

