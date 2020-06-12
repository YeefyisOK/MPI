#include "myhead.h"
void diagonal(int m, int n, int lda, MPI_Datatype *newtp, MPI_Datatype *rect) {
	MPI_Datatype oldtype[2];
	int blklen[2];
	MPI_Aint displs[2];
	MPI_Type_vector(m, n, lda, MPI_FLOAT, rect);
	oldtype[0] = *rect;
	oldtype[1] = MPI_UB;

	blklen[0] = 1;
	blklen[1] = 1;

	displs[0] = 0;
	displs[1] = sizeof(float)*(m*lda + n);
	MPI_Type_create_struct(2, blklen, displs, oldtype, newtp);
}