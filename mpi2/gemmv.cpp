#include "myhead.h"

#include <iostream>
using namespace std;
void gemmv(int n, int m, float *a, int lda, float *x, float *y) {
	//���� ����������   y=ax+y
	//�ĳ� x=Ax A m*n x n*1 
	int i, j;
	
	for (int i = 0; i < n; i++) {
		y[i] = 0;
		for (int j = 0; j < m; j++) {
			y[i] += a[i*lda + j] * x[j];
		}
	}
	return;
}