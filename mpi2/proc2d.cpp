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
	//小于p*q的才要，大于的不用了
	key = iam;
	MPI_Comm_split(comm, color, key, &valcom);
	if (valcom == MPI_COMM_NULL) return;

	/*形成行通讯子 0 1 2在一行*/
	color = iam / q;
	MPI_Comm_split(valcom, color, key, rowcom);
	/*形成列通讯子 0 3 6在一列*/
	color = iam % q;
	MPI_Comm_split(valcom, color, key, colcom);

	MPI_Comm_rank(*colcom, rowid);//列通讯子产生行id
	MPI_Comm_rank(*rowcom, colid);
	
}