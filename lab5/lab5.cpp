#include "mpi.h"
#include <iostream>
#include <vector>
 
#define BASE 10 // Система счисления
#define MIN_LENGTH_FOR_KARATSUBA 4 // Числа короче умножаются квадратичным алгоритмом
 
using namespace std;
 
void finalize(vector<int>& res) {
    for (auto i = 0; i < res.size(); ++i) {
        res[i + 1] += res[i] / BASE;
        res[i] %= BASE;
    }
}
 
vector<int> naive_mul(const vector<int>& x, const vector<int>& y) {
    auto len = x.size();
    vector<int> res(2 * len);
 
    for (auto i = 0; i < len; ++i) {
        for (auto j = 0; j < len; ++j) {
            res[i + j] += x[i] * y[j];
        }
    }
 
    return res;
}
 
vector<int> karatsuba_mul(const vector<int>& x, const vector<int>& y) {
    auto len = x.size();
    vector<int> res(2 * len);
 
    if (len <= MIN_LENGTH_FOR_KARATSUBA) {
        return naive_mul(x, y);
    }
 
    auto k = len / 2;
 
    vector<int> Xr{ x.begin(), x.begin() + k };
    vector<int> Xl{ x.begin() + k, x.end() };
    vector<int> Yr{ y.begin(), y.begin() + k };
    vector<int> Yl{ y.begin() + k, y.end() };
 
    vector<int> P1 = karatsuba_mul(Xl, Yl);
    vector<int> P2 = karatsuba_mul(Xr, Yr);
 
    vector<int> Xlr(k);
    vector<int> Ylr(k);
 
    for (int i = 0; i < k; ++i) {
        Xlr[i] = Xl[i] + Xr[i];
        Ylr[i] = Yl[i] + Yr[i];
    }
 
    vector<int> P3 = karatsuba_mul(Xlr, Ylr);
 
    for (auto i = 0; i < len; ++i) {
        P3[i] -= P2[i] + P1[i];
    }
 
    for (auto i = 0; i < len; ++i) {
        res[i] = P2[i];
    }
 
    for (auto i = len; i < 2 * len; ++i) {
        res[i] = P1[i - len];
    }
 
    for (auto i = k; i < len + k; ++i) {
        res[i] += P3[i - k];
    }
 
 
    return res;
}
 
void doMul(int*& buffer, int bufferSize) {
    vector<int> first(bufferSize / 2);
    vector<int> second(bufferSize / 2);
    vector<int> result;
 
    for (auto i = 0; i < bufferSize / 2; i++) {
        first[i] = buffer[bufferSize / 2 - 1 - i];
        second[i] = buffer[bufferSize - 1 - i];
    }
 
    result = karatsuba_mul(first, second);
    finalize(result);
 
    for (auto i = 0; i < bufferSize; i++)
        buffer[i] = 0;
 
    for (auto i = 0; i < bufferSize; i++)
        buffer[i] = result[bufferSize - 1 - i];
 
    return;
}
 
int main(int argc, char* argv[]) {
    int ProcNum, ProcRank;
    int A = 8, N = 6;
 
    int currentBufSize = N * 2;
    int blockMultiplier = 1;
 
    int loop = 0;
 
    int* numbers = new int[A * N]; //|N|N|N|N|... A раз
    for (auto i = 0; i < A * N; i++)
        numbers[i] = rand() % 10;
 
 
    int temp = A;
    while (temp != 1) {
        temp = temp / 2;
        loop += 1;
    }
 
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
 
    MPI_Status Status;
    MPI_Datatype PAIR_LONG_NUMBER;
    MPI_Type_contiguous(N * 2, MPI_INT, &PAIR_LONG_NUMBER);
    MPI_Type_commit(&PAIR_LONG_NUMBER);
    int* recvBuffer = new int[currentBufSize];
 
    MPI_Group group;
    MPI_Comm_group(MPI_COMM_WORLD, &group);
 
    int* ranks = new int[A / 2];
    for (int i = 0; i < A / 2; i++)
        ranks[i] = i;
 
    MPI_Group my_group;
    MPI_Group_incl(group, A / 2, ranks, &my_group);
 
    MPI_Comm COMM;
    MPI_Comm_create(MPI_COMM_WORLD, my_group, &COMM);
 
    if (ProcRank == 0) {
        for (auto i = 0; i < A * N; i++) {
            if (i != 0 && i % N == 0)
                cout << endl;
 
            cout << numbers[i];
        }
        cout << endl;
    }
 
    for (auto i = 0; i < loop; i++) {
 
        MPI_Scatter(numbers, blockMultiplier, PAIR_LONG_NUMBER, recvBuffer, blockMultiplier, PAIR_LONG_NUMBER, 0, COMM);
        doMul(recvBuffer, currentBufSize);
 
        MPI_Gather(recvBuffer, blockMultiplier, PAIR_LONG_NUMBER, numbers, blockMultiplier, PAIR_LONG_NUMBER, 0, COMM);
 
        MPI_Barrier(COMM);
 
        if (ProcRank == 0) {
            for (auto i = 0; i < A * N; i++) {
                if (i != 0 && i % currentBufSize == 0)
                    cout << endl;
 
                cout << numbers[i];
            }
            cout << endl;
        }
 
        blockMultiplier = blockMultiplier * 2;
        currentBufSize = currentBufSize * 2;
        recvBuffer = new int[currentBufSize];
 
    }
 
    MPI_Type_free(&PAIR_LONG_NUMBER);
    delete[] recvBuffer;
    MPI_Finalize();
 
    delete[] numbers;
 
    return 0;
}

