#include "myhead.h"
void snglscan(MPI_Comm comm, int iam, float a, int root, float *b) {
	//Ïàµ±ÓÚMPI_scan
	MPI_Comm newcomm;
	int color, key;
	if (iam <= root) {
		color = 0;
	}
	else {
		color = MPI_UNDEFINED;
	}
	key = iam;
	MPI_Comm_split(comm, color, key, &newcomm);
	if (iam <= root) {
		MPI_Reduce(&a, b, 1, MPI_FLOAT, MPI_SUM, root, newcomm);
		MPI_Comm_free(&newcomm);
	}
	return;
}
