#include <mpi.h>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <vector>

const int N=2;

class Complex {
public:
	int Re;
	int Im;
	Complex operator*(const Complex& complex2) const {
		Complex result;
		result.Re=this->Re*complex2.Re-this->Im*complex2.Im;
		result.Im=this->Re*complex2.Im+this->Im*complex2.Re;
		return result;
	}
	Complex operator+(const Complex& complex2) const {
		Complex result;
		result.Re=this->Re*complex2.Re;
		result.Im=this->Im*complex2.Im;
		return result;
	}
	Complex& operator+=(const Complex& complex2) {
		this->Re+=complex2.Re;
		this->Im+=complex2.Im;
		return *this;
	}
	friend std::ostream& operator<<(std::ostream& os, const Complex& complex);
};

std::ostream& operator<<(std::ostream& os, const Complex& complex) {
	os<< "("<<complex.Re<<", "<<complex.Im<<")";
		return os;
	}

void InitializeMatrix(Complex*& matrix);
void MultiplyMatrices(const Complex* matrix1, const Complex* matrix2, Complex*& result, const int& N);
void ShowMatrix(const Complex* matrix);

int main(int argc, char* argv[]) {
//	srand(time(NULL));
	int ProcRank;
	int CommSize;
	int A=7;
	if(A<=0) {
		std::cout<<"No matrices to multiply"<<std::endl;
		return 0;
	}
	Complex** matrices=new Complex*[A];
	Complex* matrix;
	for(int i=0; i<A; ++i) {
		matrix=new Complex[N*N];
		InitializeMatrix(matrix);
		matrices[i]=matrix;
		ShowMatrix(matrices[i]);
	}
	Complex* matrix1=new Complex[N*N];
	Complex* result=new Complex[N*N];
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &CommSize);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	MPI_Status Status;
	MPI_Datatype MPI_Complex;
	MPI_Type_contiguous(2, MPI_INT, &MPI_Complex);
	MPI_Type_commit(&MPI_Complex);
	int send_to;
	int recv_from;
	if(A==1) {
		ShowMatrix(matrices[0]);
	}
	int round=1;
	do {
		MPI_Barrier(MPI_COMM_WORLD);

		if(ProcRank==0) {
			send_to=1;
			MPI_Send(matrices[0], N*N, MPI_Complex, send_to, 0, MPI_COMM_WORLD);	
		} else {
			recv_from=ProcRank-1;
			MPI_Recv(matrix1, N*N, MPI_Complex, ProcRank-1, 0, MPI_COMM_WORLD, &Status);
			MultiplyMatrices(matrices[ProcRank],matrix1, result, N);
			if(ProcRank==CommSize-1) send_to=0;
			else send_to=ProcRank+1;
			MPI_Send(result, N*N, MPI_Complex, send_to, 0, MPI_COMM_WORLD);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		if(ProcRank==0) {
			recv_from=CommSize-1;
			MPI_Recv(result, N*N, MPI_Complex, recv_from, 0, MPI_COMM_WORLD, &Status);
		}
		++round;
	} while (CommSize*round<=A);
	if(ProcRank==0) {
		std::cout<<"Result:"<<std::endl;
		ShowMatrix(result);
	}
	for(int i=0; i<A; ++i) {
		delete[] matrices[i];
	}
	delete[] matrices;
	MPI_Type_free(&MPI_Complex);
	MPI_Finalize();
}

void InitializeMatrix(Complex*& matrix) {
	for(int i=0; i<N*N; ++i) {
		matrix[i].Re=rand()%10+1;
		matrix[i].Im=rand()%10+1;
	}
}

void MultiplyMatrices(const Complex* matrix1, const Complex* matrix2, Complex*& result, const int& N) {
	for(int i=0; i<N*N; ++i) {
		result[i].Re=0;
		result[i].Im=0;
	}
	for(int i=0; i<N; ++i) {
		for(int j=0; j<N; ++j) {
			for(int k=0; k<N; ++k) {
				result[i*N+j]+=matrix1[i*N+k]*matrix2[k*N+j];
			}
		}
	}
}

void ShowMatrix(const Complex* matrix) {
	std::cout<<std::endl;
	for(int i=0; i<N; ++i) {
		for(int j=0; j<N; ++j) {
			std::cout<< matrix[i*N+j] <<", ";
		}
		std::cout<<std::endl;
	}
	std::cout<<std::endl;
}
