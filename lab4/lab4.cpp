#include <mpi.h>
#include <iostream>
 
using namespace std;
 
struct complex {
	int Re;
	int Im;
};
 
void MtrxMultiplic(complex*& A, complex*& B, complex*& C, int& Size, int ProcNum, int ProcRank) {
	MPI_Datatype COMPLEX_NUMBER;
	MPI_Type_contiguous(2, MPI_INT, &COMPLEX_NUMBER);
	MPI_Type_commit(&COMPLEX_NUMBER);

	int dim = Size;
	int i, j, k, p, ind;
	complex temp;
	MPI_Status Status;
	int ProcPartSize = dim / ProcNum;
	int ProcPartElem = ProcPartSize * dim;
	complex* bufA = new complex[ProcPartElem];
	complex* bufB = new complex[ProcPartElem];
	complex* bufC = new complex[ProcPartElem];
	
	int ProcPart = dim / ProcNum, part = ProcPart * dim;
	if(ProcRank==0) {
		std::cout<<"Proc part size = "<<ProcPartSize<<std::endl;
		std::cout<<"Proc part elem = "<<ProcPartElem<<std::endl;
		std::cout<<"Proc part = "<<ProcPart<<std::endl;
		std::cout<<"Part = "<<part<<std::endl;
	}
	MPI_Scatter(A, part, COMPLEX_NUMBER, bufA, part, COMPLEX_NUMBER, 0, MPI_COMM_WORLD); //разбивают на равные части и посылает данные каждому процессу, что бы воссоздать полную матрицу
	MPI_Scatter(B, part, COMPLEX_NUMBER, bufB, part, COMPLEX_NUMBER, 0, MPI_COMM_WORLD);
 	if(ProcRank==1) {
		std::cout<<"Second part"<<std::endl;
		for(int i=0; i<ProcPartElem; ++i) {
			std::cout<<bufA[i].Re<<" "<<bufA[i].Im<<std::endl;
		}
		for(int i=0; i<ProcPartElem; ++i) {
			std::cout<<bufB[i].Re<<" "<<bufB[i].Im<<std::endl;
		}
	}
	temp.Re = 0;
	temp.Im = 0;
	for (i = 0; i < ProcPartSize; i++) {
		for (j = 0; j < ProcPartSize; j++) {
			for (k = 0; k < dim; k++) {
				temp.Re += bufA[i * dim + k].Re * bufB[j * dim + k].Re - bufA[i * dim + k].Im * bufB[j * dim + k].Im;
				temp.Im += bufA[i * dim + k].Re * bufB[j * dim + k].Im + bufB[j * dim + k].Re * bufA[i * dim + k].Im;
			}
 
			bufC[i * dim + j + ProcPartSize * ProcRank].Re = temp.Re;
			bufC[i * dim + j + ProcPartSize * ProcRank].Im = temp.Im;
 
			temp.Re = 0;
			temp.Im = 0;
		}
	}
 
	int NextProc; int PrevProc;
	for (p = 1; p < ProcNum; p++) {
		NextProc = ProcRank + 1;
		if (ProcRank == ProcNum - 1)
			NextProc = 0;
		PrevProc = ProcRank - 1;
		if (ProcRank == 0)
			PrevProc = ProcNum - 1;
		MPI_Sendrecv_replace(bufB, part, MPI_DOUBLE, NextProc, 0, PrevProc, 0, MPI_COMM_WORLD, &Status); //посланное сообщение замещается полученным
 
		temp.Re = 0;
		temp.Im = 0;
  	if(ProcRank==1) {
		std::cout<<"First part"<<std::endl;
		for(int i=0; i<ProcPartElem; ++i) {
			std::cout<<bufA[i].Re<<" "<<bufA[i].Im<<std::endl;
		}
		for(int i=0; i<ProcPartElem; ++i) {
			std::cout<<bufB[i].Re<<" "<<bufB[i].Im<<std::endl;
		}
	}
		for (i = 0; i < ProcPartSize; i++) {
			for (j = 0; j < ProcPartSize; j++) {
				for (k = 0; k < dim; k++) {
					temp.Re += bufA[i * dim + k].Re * bufB[j * dim + k].Re - bufA[i * dim + k].Im * bufB[j * dim + k].Im;
					temp.Im += bufA[i * dim + k].Re * bufB[j * dim + k].Im + bufB[j * dim + k].Re * bufA[i * dim + k].Im;
				}
				if (ProcRank - p >= 0)
					ind = ProcRank - p;
				else ind = (ProcNum - p + ProcRank);
 
				bufC[i * dim + j + ind * ProcPartSize].Re = temp.Re;
				bufC[i * dim + j + ind * ProcPartSize].Im = temp.Im;
 
				temp.Re = 0;
				temp.Im = 0;
			}
		}
	}
 
	MPI_Gather(bufC, ProcPartElem, COMPLEX_NUMBER, C, ProcPartElem, COMPLEX_NUMBER, 0, MPI_COMM_WORLD); //производит сборку от всех процессов в один массив C
 
	delete[]bufA;
	delete[]bufB;
	delete[]bufC;
}
 
int main(int argc, char* argv[]) {
	int ProcNum, ProcRank;
	int N = 4, A = 2;
 
	complex* matrix1 = new complex[N * N];
	complex* matrix2 = new complex[N * N];
	complex* result = new complex[N * N];
	complex* resultNotParallel = new complex[N * N];
 
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
 
 
	if (ProcRank == 0) {
 
		for (int i = 0; i < N * N; i++) {
			matrix1[i].Re = 1 + rand() % 10;
			matrix1[i].Im = 1 + rand() % 10;
 
			matrix2[i].Re = 1 + rand() % 10;
			matrix2[i].Im = 1 + rand() % 10;
 
			result[i].Re = 0;
			result[i].Im = 0;
		}
 
		cout << "Matrix1: " << endl;
		for (int i = 0; i < N * N; i++) {
			cout << "(" << matrix1[i].Re << " ," << matrix1[i].Im << ") ";
			if ((i + 1) % N == 0)
				cout << endl;
		}
 
		cout << "Matrix2: " << endl;
		for (int i = 0; i < N * N; i++) {
			cout << "(" << matrix2[i].Re << " ," << matrix2[i].Im << ") ";
			if ((i + 1) % N == 0)
				cout << endl;
		}
 
	}
 
	MtrxMultiplic(matrix1, matrix2, result, N, ProcNum, ProcRank);
	MPI_Barrier(MPI_COMM_WORLD);
 
	if (ProcRank == 0) {
		cout << "result: ";
		for (int i = 0; i < N * N; i++) {
			cout << "[" << result[i].Re << " ," << result[i].Im << "] ";
			if ((i + 1) % N == 0)
				cout << endl;
		}
 
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				resultNotParallel[i * N + j].Re = 0;
				resultNotParallel[i * N + j].Im = 0;
				for (int k = 0; k < N; k++) {
					resultNotParallel[i * N + j].Re += matrix1[i * N + k].Re * matrix2[j * N + k].Re - matrix1[i * N + k].Im * matrix2[j * N + k].Im;
					resultNotParallel[i * N + j].Im += matrix1[i * N + k].Re * matrix2[j * N + k].Im + matrix2[j * N + k].Re * matrix1[i * N + k].Im;
				}
			}
		}
 
		cout << "NotParallel: ";
		for (int i = 0; i < N * N; i++) {
			cout << "[" << resultNotParallel[i].Re << " ," << resultNotParallel[i].Im << "] ";
			if ((i + 1) % N == 0)
				cout << endl;
		}
	}
 
	MPI_Finalize();
 
	return 0;
}

