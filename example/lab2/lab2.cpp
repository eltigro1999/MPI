#include <mpi.h>
#include <iostream>

void Ring(const int& ProcAmount, const int& ProcRank);

int main(int argc, char* argv[]) {
	int ProcAmount=0;
	int ProcRank=0;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcAmount);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	int RoundsAmount=10;
	for(int Round=0; Round<RoundsAmount; ++Round) {
		if(ProcRank==0) {
			std::cout<<"Current round is "<<Round<<std::endl;
		}
		Ring(ProcAmount, ProcRank);
	}
	MPI_Finalize();
	return 0;
}

void Ring(const int& ProcAmount, const int& ProcRank) {
	int Message=0;
	for (int CurrentRank=0; CurrentRank<ProcAmount; ++CurrentRank) {
		if(CurrentRank==0) {
			Message = CurrentRank;
		} else if (ProcRank==CurrentRank){
			std::cout<<"Current rank is "<<CurrentRank<<". Message was received from "<<Message<<std::endl;
			Message+=1;
		}
		MPI_Bcast(&Message, 1, MPI_INT, CurrentRank, MPI_COMM_WORLD);
	}
	if (ProcRank==0) {
		std::cout<<"Current rank is "<<ProcRank<<". Message was received from "<<Message<<std::endl;
	}
}
