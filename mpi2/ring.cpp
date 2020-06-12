#include "myhead.h"
void ring(int m, int *n, MPI_Comm comm, int np, int iam) {
	MPI_Status st;
	int front = (np + iam - 1) % np;
	int next = (iam + 1) % np;
	/*
	if (iam % 2 == 0) {
		MPI_Send(&m, 1, MPI_INT, next, 1, comm);
		MPI_Recv(n, 1, MPI_INT, front, 1, comm, &st);
	}
	else {
		MPI_Recv(n, 1, MPI_INT, front, 1, comm, &st);
		MPI_Send(&m, 1, MPI_INT, next, 1, comm);
	}
	*/
	/*
	if (iam == 0)
		front = MPI_PROC_NULL;
	else if(iam==np-1)
		next = MPI_PROC_NULL;

	if (iam == 0)
		MPI_Send(&m, 1, MPI_INT, next, 1, comm);
	else if (iam == np - 1)
		MPI_Recv(n, 1, MPI_INT, front, 1, comm, &st);
	else
		MPI_Sendrecv(&m, 1, MPI_INT, next, 1, n, 1, MPI_INT, front,1, comm, &st);//挺好的
	*/
	/*
	//一起等待
	MPI_Status sts[2];
	MPI_Request reqs[2];
	MPI_Isend(&m, 1, MPI_INT, next, 1, comm, &reqs[0]);
	MPI_Irecv(n, 1, MPI_INT, front, 1, comm, &reqs[1]);
	MPI_Waitall(2, reqs, sts);
	*/
	//分别等待
	MPI_Request sreq,rreq;
	MPI_Isend(&m, 1, MPI_INT, next, 1, comm, &sreq);
	MPI_Irecv(n, 1, MPI_INT, front, 1, comm, &rreq);
	//MPI_Wait(&sreq, &st);
	MPI_Wait(&rreq, &st);

	//MPI_Request_free(reqs);这个为什么会出现问题
	return;
}