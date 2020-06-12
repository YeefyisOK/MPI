#include"myhead.h"
#include <iostream>
void proc2d(MPI_Comm comm, int np, int iam, int p,
	int q, MPI_Comm *rowcom, MPI_Comm *colcom, int *rowid, int *colid) {
	int color, key, pxq;
	MPI_Comm valcom;

	pxq = p * q;
	if (np < pxq) {
		return;
	}
	if (iam < pxq) {
		color = 0;
	}
	else {
		color = MPI_UNDEFINED;
	}
	//С��p*q�Ĳ�Ҫ�����ڵĲ�����
	key = iam;
	MPI_Comm_split(comm, color, key, &valcom);
	if (valcom == MPI_COMM_NULL) return;

	/*�γ���ͨѶ�� 0 1 2��һ��*/
	color = iam / q;
	MPI_Comm_split(valcom, color, key, rowcom);
	/*�γ���ͨѶ�� 0 3 6��һ��*/
	color = iam % q;
	MPI_Comm_split(valcom, color, key, colcom);

	MPI_Comm_rank(*colcom, rowid);//��ͨѶ�Ӳ�����id
	MPI_Comm_rank(*rowcom, colid);
	
}