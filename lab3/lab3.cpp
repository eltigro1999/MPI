#include <mpi.h>
#include <iostream>

void Ring(const int& CommSize, const int& ProcRank);

int main(int argc, char* argv[]) {
	int ProcRank;
	int CommSize;
	int RoundsAmount=10;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &CommSize);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	for(int Round=0; Round<RoundsAmount; ++Round) {
		Ring(CommSize, ProcRank);
		if(ProcRank==0)
		std::cout<<std::endl;	//This string is needed just for going to a new line
	}
	MPI_Finalize();
	return 0;
}

void Ring(const int& CommSize, const int& ProcRank) {
	int messageSize=CommSize;
	int* messageToSend=new int[messageSize];
	for(int i=0; i<messageSize; ++i) {
		messageToSend[i]=ProcRank*10;	//Every Process has the message in it which is just sequence of ProcRank*10
	}
	int messageToReceive=2;	// 2 is Random number here
	for(int i=0; i<CommSize; ++i) {
		//every process has to wait each other in here
		MPI_Barrier(MPI_COMM_WORLD);
		// process 0 won't be shown first as it has to receive a message from the process n-1
		if(i!=0 && i==ProcRank) {
			std::cout<<"Current rank is "<<ProcRank<<". Message is "<<messageToReceive<<". Received from "<<messageToReceive/10<<std::endl;
		}
		//scattering data
		MPI_Scatter(messageToSend, 1, MPI_INT, &messageToReceive, 1, MPI_INT, i, MPI_COMM_WORLD);
	}
	//Process 0 gets the data from the process n-1 so now it can be shown
	if(ProcRank==0) 
		std::cout<<"Current rank is "<<ProcRank<<". Message is "<<messageToReceive<<". Received from "<<messageToReceive/10<<std::endl;
	delete[] messageToSend;
}
