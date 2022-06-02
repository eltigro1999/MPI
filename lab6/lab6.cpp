#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "mpi.h"
#include <time.h>
#include <string>
 
using namespace std;
int RecvRank;
 
int ProcNum;
int ProcRank;
int GroupProcRank;
const int N = 5, A = 10;
 
 
// Убирает все 0 в старших разрядах
int change_res_size(int* num, int size) { // num - число которое надо отформатировать, size - разрядность числа
 
    while (num[size - 1] == 0 && size > 1) {
        size--;
    }
    return size; // возвращает новое значение разрядности
}
 
 
// Функция показа числа
 
void show_num(int* num, int size) {
    cout << endl;
    size = change_res_size(num, size);
    int count = (size % 3);
 
    for (int i = size - 1; i >= 0; i--) {
 
        printf("%d", num[i]);
        if (count == 3)
        {
            printf(" ");
            count = 0;
        }
        count++;
 
    }
 
};
 
// Функция создания рандомно числа разрядности N
int* createNum() {
 
    int* res = new int[N]; // N - глобальное число
    for (int i = 0; i < N; i++) {
        do {
            res[N - 1 - i] = rand() % 10;
        } while (i == 0 && res[N - 1] == 0); // Проверка что самый старший разряд не 0, то есть что все числа одного разряда
    }
 
    return res; //возвращает ссылку на число
}
 
 
/// Функция для создания оболочки числа, для удобного сложения
int* create_Null(int size) {
    int* res = new int[size];
    for (int i = 0; i < size; i++)
    {
        res[i] = 0;
    }
    return res; //возвращает ссылку на число
}
 
 
// Функция которая форматирует массив(которое позиционируется как число), таким образом чтобы в ячейках не было чисел больших 10, 
//а соответствеено лишние разряды переностcя на старше разряды
void normilize_number(int* number, int size) { // number - число, которое надо нормализовать, size - разрядность
    for (int i = 0; i < size; ++i) {
        if (number[i] > 9) {
            number[i + 1] += number[i] / 10;
            number[i] %= 10;
        }
    }
}
 
// Функция умножения чисел представленных ввиде массива
int* multyply(int* a, int* b, int size_a, int size_b) { //a,b - числа массивы, size_a,size_b - раразрядности чисел
    int res_size = size_a + size_b;
    int* res = create_Null(res_size);
    for (int Bi = 0; Bi < size_b; Bi++) {
        for (int Ai = 0; Ai < size_a; Ai++) {
            res[Ai + Bi] += a[Ai] * b[Bi];
        }
 
    }
    normilize_number(res, res_size);
    return res; //возвращает результат умножения чисел
}
 
 
int main(int argc, char* argv[])
{
    srand(time(0));
 
    MPI_Datatype LongInt; // Тип длинного числа передаваемого 0 - ым процессом
    MPI_Datatype LongIntForResFromProc; // Тип длинного числа передаваемого i - ым процессом в качестве результата выполнения умножения
    MPI_Comm Topology;
 
    MPI_Init(&argc, &argv);
 
    int ProcNum;
    int ProcRank;
    MPI_Status Status;
 
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
 
 
    MPI_Type_contiguous(N, MPI_INT, &LongInt);
    MPI_Type_contiguous(2 * N, MPI_INT, &LongIntForResFromProc);
 
    MPI_Type_commit(&LongInt);
    MPI_Type_commit(&LongIntForResFromProc);
 
 
    int ProcessInActive = ProcNum - 1; //Так как 0 в промежуточном умножении не действует
 
    while (A / ProcessInActive < 2) { //Вычисляет количество задействованных процессов
        ProcessInActive--;
 
    }
    ProcessInActive++;
 
    int* indexes = new int[ProcessInActive];
    int* edges = new int[2 * ProcessInActive];
 
    for (int i = 0; i < 2 * ProcessInActive; i++) {
 
        if (i < ProcessInActive - 1) {
            edges[i] = i + 1;
        }
        else
        {
            edges[i] = 0;
        }
 
        if (i == 0) {
            indexes[i] = ProcessInActive;
        }
        else if (i < ProcessInActive)
        {
            indexes[i] = ProcessInActive + i;
        }
    }
    MPI_Graph_create(MPI_COMM_WORLD, ProcessInActive, indexes, edges,
        1, &Topology);
 
    int cycleforProc = A;
 
    int cycle = ProcRank * 2;
 
    if (A % 2 != 0) {
 
        cycleforProc = A - 1;
    }
 
 
 
    if (ProcRank == 0)
    {
 
        int evenA = A;
 
        int* res;
 
 
 
        int res_size = 1;
        res_size = 0;
        int count = 0;
 
        if (A % 2 != 0) {
 
            evenA = A - 1;
            res = createNum();
            show_num(res, N);
            res_size = N;
            count = 1;
        }
        else
        {
            res = new int[2 * N];
        }
 
 
 
        int procNum = 1;
 
        while (count < A)
        {
            int* num1 = createNum();
            show_num(num1, N);
            int* num2 = createNum();
            show_num(num2, N);
 
 
            if (procNum == ProcNum) {
                procNum = 1;
            }
 
 
            MPI_Send(num1, 1, LongInt, procNum, 0, Topology);
            MPI_Send(num2, 1, LongInt, procNum, 0, Topology);
 
            procNum++;
 
            if (procNum == ProcessInActive) {
 
                for (int i = 1; i < ProcessInActive; i++) {
 
                    int* resFromProc = new int[N + N];
                    MPI_Recv(resFromProc, 1, LongIntForResFromProc, i, 1, Topology, &Status);
 
                    if (res_size == 0) {
 
                        res = resFromProc;
 
                    }
                    else {
 
                        res = multyply(resFromProc, res, 2 * N, res_size);
                    }
 
                    res_size += 2 * N;
 
                }
            }
 
            count += 2;
 
            delete[] num1;
            delete[] num2;
        }
 
 
        int smth = (A / 2) % (ProcessInActive - 1);
 
        for (int i = 1; i <= smth; i++) {
 
            int* resFromProc = new int[N + N];
            MPI_Recv(resFromProc, 1, LongIntForResFromProc, i, 1, Topology, &Status);
 
            if (res_size == 0) {
 
                res = resFromProc;
 
            }
            else {
 
                res = multyply(resFromProc, res, 2 * N, res_size);
            }
 
            res_size += 2 * N;
 
        }
 
        cout << endl;
        show_num(res, res_size);
        cout << endl;
        delete[] res;
    }
    else {
        while (cycle <= cycleforProc) {
 
            int* a = new int[N];
 
            MPI_Recv(a, 1, LongInt, 0, 0, Topology, &Status);
 
            int* b = new int[N];
            MPI_Recv(b, 1, LongInt, 0, 0, Topology, &Status);
 
            int* res = multyply(a, b, N, N); //Умножение 
            MPI_Send(res, 1, LongIntForResFromProc, 0, 1, Topology);
 
 
            cycle += 2 * (ProcNum - 1);
 
            delete[] a;
            delete[] b;
            delete[] res;
        }
    }
    MPI_Type_free(&LongInt);
    MPI_Type_free(&LongIntForResFromProc);
    MPI_Finalize();
 
    return 0;
}

