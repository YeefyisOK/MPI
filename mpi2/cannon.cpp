#include "myhead.h"
void cannon(MPI_Comm rowcom, MPI_Comm colcom, int p, int myrow, int mycol, int m, int k, int n,
	float *a, int lda, float *b, int ldb, float *c, int ldc, float *at, int ldaw, float *bt, int ldbw) {
	//abt是临时数组
	int i, j;
	MPI_Datatype atp, btp, attp, bttp;
	MPI_Status st;
	typemat(k, n, lda, &atp);
	MPI_Type_commit(&atp);
	typemat(m, k, ldaw, &attp);
	MPI_Type_commit(&attp);
	typemat(k, n, ldb, &btp);
	MPI_Type_commit(&btp);
	typemat(k, n, ldbw, &bttp);
	MPI_Type_commit(&bttp);

	int l = myrow;
	int front = (myrow - 1 + p) % p;//上一行
	int next = (myrow + 1 + p) % p;//下一行
	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++) {
			c[i*ldc + j] = 0.0;
		}
	}
	for (i = 0; i < p; i++) {
		if (mycol == l) scopy(m, k, a, lda, at, ldaw);
		MPI_Bcast(at, 1, attp, l, rowcom);

		gemmm(m, k, n, at, ldaw, b, ldb, c, ldc);
		if (i == p - 1) continue;
		MPI_Sendrecv(b, 1, btp, front, 1, bt, 1, bttp, next, 1, colcom, &st);

		scopy(k, n, bt, ldbw, b, ldb);
		l = (l + 1) % p;

	}

}