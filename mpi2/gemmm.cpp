#include "myhead.h"

void gemmm(int m, int k, int n, float *a, int lda, float *b, int ldb, float *c, int ldc) {
	/*
	c=c+a*b
	m k n����ά��
	m��a������
	k��b������
	n��b������
	*/
	int i, j, l;//����ѭ��
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			for (int l = 0; l < k; l++) {
				c[i*ldc + j] += a[i*lda + l] * b[l*ldb + j];
			}
		}
	}	
	return;
}

//�����������
void typemat(int m, int n, int lda, MPI_Datatype *newtp) {
	MPI_Type_vector(m, n, lda, MPI_FLOAT, newtp);
	return;
}
void scopy(int m, int n, float *a, int lda, float *b, int ldb) {
	int i, j;
	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++) {
			b[i*ldb + j] = a[i*lda + j];
		}
	}
	return;
}
//i+jΪż��ʱΪ1
void setinitab(int p, int myrow, int mycol, int m, int k, int n,
	float *a, int lda, float *b, int ldb) {
	int i, j, offsizea, offsizeb;
	offsizea = m * myrow + k * mycol;
	offsizeb = k * myrow + n * mycol;

	for (i = 0; i < m; i++) {
		for (j = 0; j < k; j++) {
			a[i*lda + j] = i + j + offsizea;
		}
	}
	for (i = 0; i < k; i++) {
		for (j = 0; j < n; j++) {
			b[i*ldb + j] = 1.0 - 2.0*((i + j + offsizeb) % 2);
		}
	}
}