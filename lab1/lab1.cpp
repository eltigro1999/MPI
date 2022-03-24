#include <mpi.h>
#include <iostream>

void Ring(const int& ProcAmount, const int& CurrentRank, MPI_Status* status);

int main(int argc, char *argv[]) {
	int ProcAmount=0, ProcRank=0;
	int RoundsAmount=5;
	MPI_Init(&argc, &argv);
	MPI_Status status;
	MPI_Comm_size(MPI_COMM_WORLD, &ProcAmount);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	for(int Round=0; Round<RoundsAmount; Round++) {
		Ring(ProcAmount, ProcRank, &status);
	}
	MPI_Finalize();
	return 0;
}

void Ring(const int& ProcAmount, const int& CurrentRank, MPI_Status* status) {
	int SendTo=0;
	int ReceivedFrom=0;
	if (CurrentRank==0) {
		SendTo=1;
		MPI_Send(&CurrentRank, 1, MPI_INT, SendTo, CurrentRank, MPI_COMM_WORLD);
		std::cout<<"Message from "<<CurrentRank<<" To "<<SendTo<<" was sent. Current rank is "<<CurrentRank<<std::endl;
		ReceivedFrom=ProcAmount-1;
		MPI_Recv(&ReceivedFrom, 1, MPI_INT, ReceivedFrom, MPI_ANY_TAG, MPI_COMM_WORLD, status);
		std::cout<<"Message from "<<ReceivedFrom<<" To "<<CurrentRank<<" was received. Current rank is "<<CurrentRank<<std::endl;
	} else {
		if (CurrentRank==ProcAmount-1) {
			SendTo=0;
			ReceivedFrom=ProcAmount-2;
		} else {
			SendTo=CurrentRank+1;
			ReceivedFrom=CurrentRank-1;
		}

		MPI_Recv(&ReceivedFrom, 1, MPI_INT, ReceivedFrom, MPI_ANY_TAG, MPI_COMM_WORLD, status);
		std::cout<<"Message from "<<ReceivedFrom<<" To "<<CurrentRank<<" was received. Current rank is "<< CurrentRank<<std::endl;
		MPI_Send(&CurrentRank, 1, MPI_INT, SendTo, 0, MPI_COMM_WORLD);
		std::cout<<"Message from "<<CurrentRank<<" To "<<SendTo<<" was sent. Current rank is "<<CurrentRank<<std::endl;
	}
}
