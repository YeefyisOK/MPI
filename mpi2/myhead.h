#pragma once
#include <mpi.h>

#define maxrow 100
#define maxcol 100
#define maxa 101

typedef struct {
	int a;
	float b[2];
	char c[3];
}abc;
/*
void mybegin(int *argc, char ***argv, MPI_Comm comm, int np, int iam);
void myend();
void gemmm(int m, int k, int n, float *a, int lda, float *b, int ldb, float *c, int ldc);
void typemat(int m, int n, int lda, MPI_Datatype *newtp);*/
void ring(int m, int *n, MPI_Comm comm, int np, int iam);
void datatype(char which, MPI_Datatype *newtype); 
void mpistruct(MPI_Datatype *newtp);
void diagonal(int m, int n, int lda, MPI_Datatype *newtp, MPI_Datatype *rect);
//矩阵乘法 行列分块
void init_a(int m, int k, int lda, float *a, int iam);
void init_b(int k, int n, int ldb, float *b, int iam);
void init_c(int m, int n, int ldc, float *c, int iam);
void matmul(int m, int k, int n, int lda, float *a, int ldb, float * b, int ldc, float *c);
void rcmatmul(MPI_Comm comm, int np, int iam, int m, int k, int n,
	int lda, float * a, int ldb, float * b, int ldc, float *c, int ldw, float * w);

//矩阵乘法 行行分块
void init_a2(int m_, int k, int lda, float *a, int iam); 
void init_b2(int k_, int n, int ldb, float *b, int iam);
void matmul2(int m_, int k_, int n, int lda, float *a, int ldb, float * b, int ldc, float *c);
void hanghang(MPI_Comm comm, int np, int iam, int m, int k, int n,
	int lda, float * a, int ldb, float * b, int ldc, float * c, int ldw, float * w);



void proc2d(MPI_Comm comm, int np, int iam, int p,
	int q, MPI_Comm *rowcom, MPI_Comm *colcom, int *rowid, int *colid);
void snglscan(MPI_Comm comm, int iam, float a, int root, float *b);
void gemmv(int m, int n, float *a, int lda, float *x, float *y);
void iteration(MPI_Comm comm, int np, int iam, int n,
	int en, float *a, int lda, float *b, float *x, int num);
void gemmm(int m, int k, int n, float *a, int lda, float *b, int ldb, float *c, int ldc);
void typemat(int m, int n, int lda, MPI_Datatype *newtp);
void scopy(int m, int n, float *a, int lda, float *b, int ldb); 
void setinitab(int p, int myrow, int mycol, int m, int k, int n,
	float *a, int lda, float *b, int ldb);
void cannon(MPI_Comm rowcom, MPI_Comm colcom, int p, int myrow, int mycol, int m, int k, int n,
	float *a, int lda, float *b, int ldb, float *c, int ldc, float *at, int ldaw, float *bt, int ldbw);




void init_ax(int m, int n, int lda, float *a, int iam, float *x);
void kaoshi(MPI_Comm comm, int p, int iam, int n,
	int m, float *a, int lda, float *b, float *x);