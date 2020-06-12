#include "myhead.h"

//矩阵A按行分块
void init_a2(int m_, int k, int lda, float *a, int iam) {
	int i, j, offside;
	offside = iam * m_;
	for (i = 0; i < m_; i++) {
		for (j = 0; j < k; j++) {
			a[i*lda + j] = offside + i + j;
		}
	}
	return;
}
//矩阵B按行分块
void init_b2(int k_, int n, int ldb, float *b, int iam) {
	int i, j, offside;
	offside = iam * k_;
	for (i = 0; i < k_; i++) {
		for (j = 0; j < n; j++) {
			b[i*ldb + j] = 1.0 - 2.0*((i + j + offside) % 2);
		}
	}
	return;
}
void matmul2(int m_, int k_, int n, int lda, float *a, int ldb, float * b, int ldc, float *c) {
	/*
	c=c+a*b
	m k n矩阵维数
	m是a的行数
	k是b的行数
	n是b的列数
	*/
	int i, j, l;//三重循环
	for (int i = 0; i < m_; i++) {
		for (int j = 0; j < n; j++) {
			for (int l = 0; l < k_; l++) {
				c[i*ldc + j] += a[i*lda + l] * b[l*ldb + j];
			}
		}
	}
	return;
}

void hanghang(MPI_Comm comm, int np, int iam, int m, int k, int n,
	int lda, float * a, int ldb, float * b, int ldc, float * c, int ldw, float * w) {
	int i, j, front, next, l;
	MPI_Datatype rectb, rectw;
	MPI_Status st;
	int m_ = m / np;
	int k_ = k / np;
	int n_ = n / np;
	MPI_Type_vector(k_, n , ldb, MPI_FLOAT, &rectb);
	MPI_Type_vector(k_, n ,ldw, MPI_FLOAT, &rectw);

	MPI_Type_commit(&rectb);
	MPI_Type_commit(&rectw);
	for (i = 0; i < m_; i++) {
		for (j = 0; j < n; j++) {
			c[i*ldc + j] = 0.0;
		}
	}

	front = (np + iam - 1) % np;
	next = (iam + 1) % np;
	int r = 0;
	for (i = 0; i < np - 1; i++) {

		if (i % 2 == 0) {
			matmul2(m_, k_, n, lda, &a[i*k_], ldb, b, ldc, &c[0]);
			MPI_Sendrecv(b, 1, rectb, front, 1,
						w, 1, rectw, next, 1, comm, &st);

			//cout << "proc=" << iam << " c[1][1] =" << c[1] << " c[1][2] =" << c[2] <<
				//" c[1][3]= " << c[3] << " c[1][4]= " << c[4] << endl;
		}
		else {
			matmul2(m_, k_, n, lda, &a[i*k_], ldw, w, ldc, &c[0]);
			MPI_Sendrecv(w, 1, rectw, front, 1, 
						b, 1, rectb, next, 1, comm, &st);

		}
	}

	if ((np - 1) % 2 == 0)
		matmul2(m_, k_, n, lda, &a[(np - 1)*k_], ldb, b, ldc, &c[0]);
	else
		matmul2(m_, k_, n, lda, &a[i*k_], ldw, w, ldc, &c[0]);
	MPI_Type_free(&rectb);
	MPI_Type_free(&rectw);
	return;
}