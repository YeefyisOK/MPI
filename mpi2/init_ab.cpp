#include "myhead.h"
#include <iostream>
using namespace std;
//矩阵A按行分块
void init_a(int m, int k, int lda, float *a, int iam) {
	int i, j, offside;
	offside = iam * m;
	for (i = 0; i < m; i++) {
		for (j = 0; j < k; j++) {
			a[i*lda+j] = offside + i + j;
		}
	}

	return;
}
//矩阵B按列分块
void init_b(int k, int n, int ldb, float *b, int iam) {
	int i, j, offside;
	offside = iam * n;
	for (i = 0; i < k; i++) {
		for (j = 0; j < n; j++) {
			b[i*ldb+j] = 1.0 - 2.0*((i + j + offside) % 2);
		}
	}
	return;
}

void init_c(int m, int n, int ldc, float *c, int iam) {
	int i, j, offside;
	offside = iam * n;
	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++) {
			c[i*ldc + j] = 1111;
		}
	}
	return;
}
void matmul(int m, int k, int n,int lda, float *a, int ldb, float * b, int ldc, float *c) {
	int i, j, l;
	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++) {
			c[i*ldc + j] = 0.0;
			for (l = 0; l < k; l++) {
				c[i*ldc + j] += a[i*lda + l]*b[l*ldb + j];
			}
		}
	}
	return;
}
void rcmatmul(MPI_Comm comm, int np, int iam, int m, int k, int n,
	int lda, float * a, int ldb, float * b, int ldc, float * c, int ldw, float * w) {
	int i, front, next, l,j;
	MPI_Datatype rectb, rectw;
	MPI_Status st;
	MPI_Type_vector(k, n, ldb, MPI_FLOAT, &rectb);
	MPI_Type_vector(k, n, ldw, MPI_FLOAT, &rectw);

	MPI_Type_commit(&rectb);
	MPI_Type_commit(&rectw);
	l = iam * n;
	front = (np + iam - 1) % np;
	next = (iam + 1) % np;

	for (i = 0; i < np - 1; i++) {
		if (i % 2 == 0) {
			matmul(m, k, n, lda, a, ldb, b, ldc, &c[l]);		
			MPI_Sendrecv(b, 1, rectb, front, 1, 
						w, 1, rectw, next, 1, comm, &st);

		}
		else {
			matmul(m, k, n, lda, a, ldw, w, ldc, &c[l]);
			MPI_Sendrecv(w, 1, rectw, front, 1, 
						b, 1, rectb, next, 1, comm, &st);

		}
		l += n;
		if (np * n == l)	
			l = 0;
	}
	
	if ((np - 1) % 2 == 0)
		matmul(m, k, n, lda, a, ldb, b, ldc, &c[l]);
	else
		matmul(m, k, n, lda, a, ldw, w, ldc, &c[l]);
	
	MPI_Type_free(&rectb);
	MPI_Type_free(&rectw);
	return;
}