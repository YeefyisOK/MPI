#include "myhead.h"
void datatype(char which,MPI_Datatype *newtype) {
	//用哪个字符表示数据类型
	int count = 2;
	int stride = 25,length=3;
	int lengths[2] = { 3,2 }, displs[2] = { 0,5 };
	if (which == 'c')
		MPI_Type_contiguous(count, MPI_INT, newtype);
	else if (which == 'v')
		MPI_Type_vector(count, length, stride, MPI_INT, newtype);
	else if (which == 'i')
		MPI_Type_indexed(count, lengths, displs, MPI_INT, newtype);

	return;

}