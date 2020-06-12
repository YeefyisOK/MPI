#include "myhead.h"
#include <stdlib.h>
#include <iostream>
using namespace std;

//jacobi����  x=(I-A)x+b 
//�����зֿ�
void iteration(MPI_Comm comm, int np, int iam, int n,
	int en, float *a, int lda, float *b, float *x, int num) {
	//ÿ�������������=n/np  enÿ������������� num��������
	int i, j, *rc;
	float *y;
	rc = (int*)malloc(np * sizeof(int));
	for (int i = 0; i < np; i++) {
		rc[i] = en;
	}
	y = (float*)malloc(n * sizeof(float));
	for (int i = 0; i < num; i++) {
		if (iam == 0) {
			for (int j = 0; j < n; j++) {
				y[j] = b[j];
			}
		}
		else {
			for (int j = 0; j < n; j++) {
				y[j] = 0;
			}
		}
		gemmv(n, en, a, lda, x, y);
		/*
		for (int k = 0; k < n; k++) {
			cout <<"y"<< k << "x��ֵΪ" << y[k] << " ";
		}
		cout << endl;*/
		MPI_Reduce_scatter(y, x, rc, MPI_FLOAT, MPI_SUM, comm);
	}
	free(y);
	free(rc);
	//�������Ax=b,A�ǶԽǾ���aii=1/2,bi=i
}


void init_ax(int m, int n, int lda, float *a, int iam,float *x) {
	//k*n��Ϊ n*m
	int i, j, offside;
	offside = iam * m;
	for (i = 0; i < n; i++) {
		for (j = 0; j < m; j++) {
			a[i*lda + j] = i + j + offside;
		}
	}
	int offisidex = iam * m;
	for (int i = 0; i < m; i++) {
		x[i] = 1.0 - 2.0*((i+offisidex) % 2);
	}
	return;
}
/*
if (iam == 1) {
	//���A
	std::cout << iam << std::endl;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			std::cout << a[i*lda+j] << " ";
		}
		std::cout << std::endl;
	}
	//���xx
	for (int i = 0; i < n; i++) {
		std::cout << x[i] << " ";

		std::cout << std::endl;
	}
}*/

void kaoshi(MPI_Comm comm, int p, int iam, int n,
	int m, float *a, int lda, float *b, float *x) {
	//A��n*m��x��m*1 y,b n*1, 
	float *y;
	y = (float*)malloc(n * sizeof(float));

	float *tl;
	tl = (float*)malloc(m * sizeof(float));

	MPI_Status st;
	//�㲥b
	if (iam == 0) {
		for (int i = 1; i < p ; i++) {
			
			MPI_Send(&b[i*m], m, MPI_FLOAT,i, 5, comm);
		}
	}
	else {
		//�������̽���b����b0λ��
		MPI_Recv(&b[0], m, MPI_FLOAT, 0, 5, comm, &st);
	}

	//A*x
	gemmv(n, m, &a[0], lda, &x[0], &y[0]);
	//����ֵ �õ�y��һ����
	int offiside = iam * m;
	for (int i = 0; i < m; i++) {
		x[i] = y[offiside+i];
	}
	/*
	if (iam == 0) {
		//���l
		for (int i = 0; i < n; i++) {
			std::cout << y[i] << " ";

			std::cout << std::endl;
		}
	}*/
	//��y��Ӧλ�ø����н���
	for (int i = 0; i < p; i++) {
		if (iam != i) {

			//MPI_Send(&y[i*m], m, MPI_FLOAT, i, 4, comm);
			MPI_Sendrecv(&y[i*m], m, MPI_FLOAT, i, 4,
				&tl[0], m, MPI_FLOAT, i, 4, comm, &st);
			for (int j = 0; j < m; j++) {
				x[j] += tl[j];
			}

		}
	}
	/*
	//���ն�Ӧλ�õ�y,��ӵõ�l
	for (int i = 0; i < p; i++) {
		if (iam != i) {
			MPI_Recv(&tl[0], m, MPI_FLOAT, i, 4, comm, &st);
			for (int j = 0; j < m; j++) {
				x[j] += tl[j];
			}
		}
	}*/
	
	

	for (int j = 0; j < m; j++) {
		x[j] += b[j];
	}
	/*
	for (int j = 0; j < m; j++) {
		x[j] = l[j];
	}*/
	free(y);
	free(tl);
}