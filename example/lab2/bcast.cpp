#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include <iostream>

int main(int argc, char* argv[])
{
	double x=120;
	int ProcRank, ProcNum;
	MPI_Status Status;
	// инициализация
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD,&ProcRank);
	// подготовка данных
	// рассылка данных на все процессы
	if(ProcRank==0)
		x=100;
	std::cout<<x<<"Current rank is"<<ProcRank<<std::endl;
	//MPI_Bcast(&x, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Finalize();
}
