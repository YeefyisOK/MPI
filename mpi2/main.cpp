#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <mpi.h>
#include "myhead.h"

using namespace std;
#define maxnp 16
int main(int argc, char **argv)
{
	MPI_Comm comm;
	int np, iam;
	int m, k, n, narray[100], marray[101], matr[11][25];
	MPI_Datatype newtp, rect;
	MPI_Status st;
	abc x[10];
	MPI_Aint sizeabc, extnewtp, extrect, lb, lb1;//什么数据类型
	int i, j;
	//float a[31][57], b[53][59], c[31][61], w[51][53],u[37][41];
	float a[44][44], b[53][59], c[44][61], w[51][53], u[37][41];

	//cannon w u两个临时的工作空间
	int rcounts[maxnp];//gatherv 接受数目数组
	int displs[maxnp];//gatherv 接受位移数组

	FILE* fp;
	int bnp[3];//读取矩阵
	typedef struct floatint
	{
		float a;
		int m;
	};
	floatint mxl, resmxl;

	//MPI_Split的参数
	MPI_Comm rowcom, colcom;
	int rowid, colid;
	//雅可比迭代的参数
	float xx[44];
	float rhs[44];
	int en;

	//组操作的参数
	MPI_Group grp1, newgrp;
	int ranks[10];
	int gnp = 111, giam = 111;

	int p;

	//修改这里就是修改整个的
#define chkkaoshi
	/*进入MPI环境*/
	//mybegin(&argc, &argv, comm, np, iam);
	/*整个这4句，第一句没有& & 就是mybegin*/
	MPI_Init(&argc, &argv);
	MPI_Comm_dup(MPI_COMM_WORLD, &comm);
	MPI_Comm_size(comm, &np);//num of process
	MPI_Comm_rank(comm, &iam);//当前进程号
	/*函数主体*/
#ifdef chkhanghang
	if (iam == 0) {

		bnp[0] = 8;
		bnp[1] = 16;
		bnp[2] = 16;
		/*
		bnp[0] = 2 * 4;
		bnp[1] = 3 * 4;
		bnp[2] = 4 * 4;*/
	}
	MPI_Bcast(bnp, 3, MPI_INT, 0, comm);
	m = bnp[0];
	k = bnp[1];
	n = bnp[2];

	init_a2(m / np, k, 57, &a[0][0], iam);
	init_b2(k / np, n, 59, &b[0][0], iam);
	if (iam == 0) {
		cout << iam << endl;
		for (i = 0; i < m / np; i++) {
			for (j = 0; j < k; j++) {
				cout << a[i][j] << " ";
			}
			cout << endl;
		}
	}
	/*
	if (iam == 1) {
		cout << iam << endl;
		for (i = 0; i < k / np; i++) {
			for (j = 0; j < n; j++) {
				cout << b[i][j] << " ";
			}
			cout << endl;
		}
	}*/
	hanghang(comm, np, iam, m, k, n, 57, &a[0][0], 59, &b[0][0],
		61, &c[0][0], 53, &w[0][0]);
	/*
	if (iam == 3) {
		cout << iam << endl;
		for (i = 0; i < m/np; i++) {
			for (j = 0; j < n; j++) {
				cout << c[i][j] << " ";
			}
			cout << endl;
		}
	}
	*/
	cout << "proc=" << iam << " c[1][1] =" << c[1][1] << " c[1][2] =" << c[1][2] <<
		" c[1][3]= " << c[1][3] << " c[1][4]= " << c[1][4] << endl;

#endif //chkhanghang
#ifdef chkcannon
	p = 3;
	if (np < 9) return 0;
	proc2d(comm, np, iam, p, p, &rowcom, &colcom, &rowid, &colid);
	if (iam == 0) {
		/*
		fp = fopen("inputmkn.txt", "r");
		i = fscanf_s(fp, "%*[^\n]%*c");
		i = fscanf_s(fp, "%*[^\n]%*c %d %d %d", &bnp[0],&bnp[1], &bnp[2]);
		*/
		bnp[0] = 11;
		bnp[1] = 11;
		bnp[2] = 12;
	}
	MPI_Bcast(bnp, 3, MPI_INT, 0, comm);
	m = bnp[0];
	k = bnp[1];
	n = bnp[2];
	if (iam < 9) {
		setinitab(p, rowid, colid, m, k, n, &a[0][0], 57, &b[0][0], 59);
		cannon(rowcom, colcom, p, rowid, colid, m, k, n, &a[0][0], 57,
			&b[0][0], 59, &c[0][0], 61, &w[0][0], 53, &u[0][0], 41);
		cout << "proc=" << iam << " rowid=" << rowid << " colid=" << colid << " c[1][1] =" << c[1][1] << " c[1][2] =" << c[1][2] <<
			" c[1][3]= " << c[1][3] << " c[1][4]= " << c[1][4] << endl;
	}
#endif //chkcannon
#ifdef chkgroup
	MPI_Comm_group(comm, &grp1);
	ranks[0] = 1;
	ranks[1] = 3;

	MPI_Group_incl(grp1, 2, ranks, &newgrp);
	if (iam == 1 || iam == 3) {
		MPI_Group_size(newgrp, &gnp);
		MPI_Group_rank(newgrp, &giam);
	}
	std::cout << "the proc " << iam << ",group" << giam << endl;

#endif //chkgroup

#ifdef chkkaoshi
	/*
	n = 8;
	m = 2;*/
	n = 11*np;
	m = 11;
	p = np;
	if (iam == 0) {
		for (int i = 0; i < n; i++) {
			rhs[i] = 1.0 - 2.0*(i % 2);
		}
	}
	init_ax(m, n, 44,&a[0][0], iam, &xx[0]);

	kaoshi(comm, p, iam, n, m, &a[0][0], 44, &rhs[0], &xx[0]);
	
	if (iam ==1) {
		//输出x
		for (i = 0; i < m; i++) {
			std::cout << xx[i] << " ";
			std::cout << std::endl;
		}
	}
	
	/*
	
	if (iam == 1) {
		//输出A
		std::cout << iam << std::endl;
		for (i = 0; i < n; i++) {
			for (j = 0; j < m; j++) {
				std::cout << a[i][j] << " ";
			}
			std::cout << std::endl;
		}
		//输出B
		for (i = 0; i < n; i++) {
			std::cout << xx[i] << " ";

			std::cout << std::endl;
		}
	}
	*/


#endif //chkkaoshi

#ifdef chkiteration
	//迭代求解Ax=b,A是对角矩阵，aii=1/2,bi=i
	en = 5;
	n = en * np;
	for (int i = 0; i < n; i++) rhs[i] = i;
	int offside = iam * en;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < en; j++) {
			a[i][j] = 0;
			if (i == (j + offside)) a[i][j] = 0.5;
		}
	}

	/*jacobi迭代  x=(I-A)x+b    */
	iteration(comm, np, iam, n, en, &a[0][0], 57, rhs, xx, 10);
	std::cout << "iam=" << iam << " x0=" << xx[0] << " x1=" << xx[1]
		<< " x2=" << xx[2] << endl;

#endif //chkiteration

#ifdef chksngl

	a[0][0] = iam + 1.0;
	snglscan(comm, iam, a[0][0], 2, &b[0][0]);
	if (iam == 2) {
		std::cout << "sum:" << b[0][0] << endl;
	}
	
	MPI_Scan(&a[0][0], &b[0][0], 1, MPI_FLOAT, MPI_SUM, comm);
	std::cout << "Each proc value =" << b[0][0]
		<< " in " << iam << endl;
#endif // chksngl

#ifdef chkproc2d
	if (np < 12) return 0;	
	proc2d(comm, np, iam, 3, 4, &rowcom, &colcom, &rowid, &colid);
	std::cout << "进程号:" << iam<<" rowid:"<<rowid<<" colid: "<<colid << endl;
	
#endif // chkproc2d

#ifdef chkreduce
	mxl.a = (iam + 1) * 20;
	mxl.m = iam;
	MPI_Allreduce(&mxl, &resmxl, 1, MPI_FLOAT_INT, MPI_MAXLOC, comm);
	cout << "proc=" << iam << "最大值 =" << resmxl.a << " loc=" << resmxl.m<< endl;

#endif // chkreduce

#ifdef chkmatmul
	if (iam == 0) {
		/*
		float a[31][57], b[53][59], c[31][61], w[51][53],u[37][41];

		*/
		bnp[0] = 2;
		bnp[1] = 3;
		bnp[2] = 4;
	}
	MPI_Bcast(bnp, 3, MPI_INT, 0, comm);
	m = bnp[0];
	k = bnp[1];
	n = bnp[2];
	init_a(m, k, 57, &a[0][0], iam);
	init_b(k, n, 59, &b[0][0], iam);
	init_c(m, n*np, 61, &c[0][0], iam);

	rcmatmul(comm, np, iam, m, k, n, 57, &a[0][0], 59, &b[0][0], 
		61, &c[0][0], 53, &w[0][0]);


	cout << "proc=" << iam <<" c[1][1] =" << c[1][1] << " c[1][2] =" << c[1][2] <<
		" c[1][3]= " << c[1][3] << " c[1][4]= " << c[1][4] << endl;
#endif // chkmatmul


#ifdef  chkgather
	if (iam == 0) {
		for (i = 0; i < 31; i++) {
			for (j = 0; j < 31; j++) {
				a[i][j] = i + j;
			}
		}
	}
	MPI_Bcast(&a[0][0], np*5, MPI_FLOAT, 0, comm);
	/*
	//把第0行20个数广播
	cout << "proc=" << iam << "a[0][0] =" << a[0][0] << " a[0][1] =" << a[0][1] <<
		" a[0][2]= " << a[0][2] << endl;
	MPI_Scatter(&a[0][0], 5, MPI_FLOAT, &a[1][0], 5, MPI_FLOAT, 0, comm);
	//MPI_Scatter0进程的20个数5个5个的分给了其他所有进程的第一行

	MPI_Allgather(&a[1][0], 5, MPI_FLOAT, &a[2][0], 5, MPI_FLOAT, comm);
	//MPI_Allgather每个进程的5个5个的数都被收集给每一个进程的第二行
	cout << "proc=" << iam << "a[2][0] =" << a[2][0] << " a[2][1] =" << a[2][1] <<
		" a[2][2]= " << a[2][2] << endl;*/
	MPI_Alltoall(&a[0][0], 5, MPI_FLOAT, &a[1][0], 5, MPI_FLOAT, comm);
	cout << "proc=" << iam << "a[1][0] =" << a[1][0] << " a[0][1] =" << a[1][1] <<
		" a[0][2]= " << a[1][2] << endl;
	a[1][1] = iam + 0.1;
	
	MPI_Reduce(&a[1][0], &a[0][0], 1, MPI_FLOAT, MPI_SUM, 0, comm); 
	cout << "proc=" << iam << "a[0][0] =" << a[0][0] << " a[0][1] =" << a[0][1] <<
		" a[0][2]= " << a[0][2] << endl;
	//MPI_Reduce求和,结果为30
	/*
	narray[0] = (iam + 1) * 20;
	narray[1] = iam;
	MPI_Reduce(narray, marray, 1, MPI_2INT, MPI_MAXLOC, 0, comm);
	i = marray[1];

	cout << "proc=" << iam << "最大值 =" << marray[0] << " loc=" << i<<endl;
	*/

	/*
		<< " a[0][3] =" << a[0][3] <<
		" a[0][4]= " << a[0][4] << " a[0][5] =" << a[0][5] <<endl;*/
	/*
	for (i = 0; i < 31; i++) {
		for (j = 0; j < 31; j++) {
			a[i][j] = i + j;
		}
	}
	j = iam * 5;
	MPI_Gather(&a[0][j], 3, MPI_FLOAT, &a[1][0], 3, MPI_FLOAT, 0, comm);
	if(iam==0)
		cout << "proc=" << iam << "a[1][0] =" << a[1][0] << " a[1][1] =" << a[1][1] <<
		" a[1][2]= " << a[1][2] << " a[1][3] =" << a[1][3] <<
		" a[1][4]= " << a[1][4] << " a[1][5] =" << a[1][5] << endl;
	MPI_Scatter(&a[1][0], 3, MPI_FLOAT, &a[0][0], 3, MPI_FLOAT, 0, comm);

	cout << "proc=" << iam << "a[0][0] =" << a[0][0] << " a[0][1] =" << a[0][1] <<
		" a[0][2]= " << a[0][2] <<endl;

	for (i = 0; i < np; i++) {
		rcounts[i] = 3;
		displs[i] = 5 * i;
	}
	MPI_Gatherv(&a[0][j], 3, MPI_FLOAT, &a[1][0], rcounts, displs, MPI_FLOAT, 0, comm);
	cout << "proc=" << iam << "a[1][0] =" << a[1][0] << " a[1][1] =" << a[1][1] <<
		" a[1][2]= " << a[1][2] << " a[1][5] =" << a[1][5] <<
		" a[1][6]= " << a[1][6] << " a[1][7] =" << a[1][7] << endl;
		*/
#endif //  chkgather

#ifdef  chkbcast
	if (iam == 0) {
		for (i = 0; i < 31; i++) {
			for (j = 0; j < 31; j++) {
				a[i][j] = i + j;
			}
		}
	}
	MPI_Bcast(&a[0][0], 5, MPI_FLOAT, 0, comm);
	cout << "proc=" << iam << "a[0][0] =" << a[0][0] << " a[0][1] =" << a[0][1] <<
	" a[0][2]= " << a[0][2] << " a[0][3] =" << a[0][3] << endl;

#endif //  chkbcast


#ifdef  chkdiagnoal
	diagonal(2, 3, 57, &newtp, &rect);
	MPI_Type_commit(&newtp);
	MPI_Type_commit(&rect);
	MPI_Type_get_extent(newtp, &lb, &extnewtp);
	MPI_Type_get_extent(rect, &lb1, &extrect);
	cout << "extrect=" << extrect << "extnewtp=" << extnewtp << endl;
	if (iam == 0) {
		for (i = 0; i < 31; i++) {
			for (j = 0; j < 31; j++) {
				a[i][j] = i + j;
			}
		}
		MPI_Send(a, 3, newtp, 1, 5, comm);

	}
	if (iam == 1) {

		MPI_Recv(a, 3, newtp, 0, 5, comm, &st);

		cout << "a[0][0] =" << a[0][0] << " a[0][1] =" << a[0][1] << 
			" a[1][0]= " << a[1][0] << " a[1][1] =" << a[1][1] <<
			" a[2][0] =" << a[2][0] << " a[2][1] =" << a[2][1] <<
			" a[2][3]= " << a[2][3] << " a[2][4] =" << a[2][4] <<
			endl;
	}

	MPI_Type_free(&newtp);
	MPI_Type_free(&rect);
#endif //  chkdiagnoal


#ifdef chkstruct
	mpistruct(&newtp);
	MPI_Type_commit(&newtp);
	if (iam == 0) {
		for (m = 0; m < 10; m++) {
			x[m].a = m;
			x[m].b[0] = 20.0*(m + 1);
			x[m].b[1] = 30.0*(m + 1);
			x[m].c[0] = 'a' + 3 * m;
			x[m].c[1] = 'b' + 3 * m;
			x[m].c[2] = 'c' + 3 * m;

		}
		MPI_Send(x, 3, newtp, 1, 5, comm);
		sizeabc = sizeof(abc);//什么数据类型
		//MPI_Type_extent(newtp, &extnewtp);

		MPI_Type_get_extent(newtp, &lb,&extnewtp);
		cout << "sizeof=" << sizeabc << "extnewtp=" << extnewtp << endl;
	}	
	if (iam == 1) {
		MPI_Recv(x, 3, newtp, 0, 5, comm, &st);
		cout << "x0的结果" << x[0].a << " " << x[0].b[0] << " " << x[0].b[1] <<
			" " << x[0].c[0] << " " << x[0].c[1] << " " << x[0].c[2] << endl;
		cout << "x1的结果" << x[1].a << " " << x[1].b[0] << " " << x[1].b[1] <<
			" " << x[1].c[0] << " " << x[1].c[1] << " " << x[1].c[2] << endl;
	}

	MPI_Type_free(&newtp);

#endif //chkstruct

#ifdef chkdatatype
	for (m = 0; m < 100; m++) narray[m] = m;
	char a = 'i';
	datatype(a, &newtp);

	MPI_Type_commit(&newtp);
	if (iam == 0) {
		/*
		for (m = 0; m < 10; m++) {
			for (n = 0; n < 25; n++) {
				matr[m][n] = m + n;
			}
			MPI_Send(matr, 1, newtp, 1, 5, comm);
		}
		*/
		MPI_Send(narray, 1, newtp, 1, 5, comm);
	}
	if (iam == 1) {
		//MPI_Recv(matr, 1, newtp, 0, 5, comm, &st);
		MPI_Recv(marray, 1, newtp, 0, 5, comm, &st);
		cout << "data from proc" << iam << "is" << marray[0] << "," << marray[1] << ","
			<< marray[2] << "," << marray[3] << "," << marray[4] << "," << marray[5]<<"," << marray[6]
			<< endl;
		/*
		cout << "data from proc" << iam << "is" << matr[0][0] << "," << matr[0][1] << ","
			<< matr[0][2] << "," << matr[1][0] << "," << matr[1][1] << "," << matr[1][2]
			<< endl;
		for (m = 0; m < 10; m++) {
			for (n = 0; n < 25; n++) {
				cout<<matr[m][n]<<" ";
			}
			cout << endl;
		}*/
	}
	MPI_Type_free(&newtp);
#endif // chkdatatype


	//m = iam + 3;
	//std::cout << "helloworld from proc" << iam << ",m=" << m << endl;
#ifdef chkring
	m = iam;
	n = 100;
	ring(m, &n, comm, np, iam);
	cout << "helloworld from proc" << iam << ",n=" << n << endl;
#endif
	/*结束*/
	MPI_Finalize();
}