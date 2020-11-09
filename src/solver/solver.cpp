#include <chrono>
#include <cstring>
#include <thread>
#include <iostream>

#include <mpi.h>

#include "solver.hpp"
#include "../matrix/matrix.hpp"

using std::this_thread::sleep_for;
using std::chrono::microseconds;
using std::cout;
using std::endl;
using std::flush;
using std::stod;
using std::stoi;

void solveSeq(int rows, int cols, int iterations, double td, double h, int sleep, double * matrix) {
    double c, l, r, t, b;
    
    double h_square = h * h;

    double * linePrevBuffer = new double[cols];
    double * lineCurrBuffer = new double[cols];

    for(int k = 0; k < iterations; k++) {

        memcpy(linePrevBuffer, matrix, cols * sizeof(double));
        for(int i = 1; i < rows - 1; i++) {

            memcpy(lineCurrBuffer, matrix + i * cols, cols * sizeof(double));
            for(int j = 1; j < cols - 1; j++) {
                c = lineCurrBuffer[j];
                t = linePrevBuffer[j];
                b = getMatrixCell(matrix, i + 1, j, cols);
                l = lineCurrBuffer[j - 1];
                r = lineCurrBuffer[j + 1];

                sleep_for(microseconds(sleep));
                matrix[i * cols + j] = c * (1.0 - 4.0 * td / h_square) + (t + b + l + r) * (td / h_square);
            }

            memcpy(linePrevBuffer, lineCurrBuffer, cols * sizeof(double));
        }
    }
}

void solvePar(int rows, int cols, int iterations, double td, double h, int sleep, double * matrix) {
    int rank;
    int processCount;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &processCount);

    double c, l, r, t, b;
    double h_square = h * h;

    int *sendcounts;            // array describing how many elements to send to each process
    int *displs;                // array describing the displacements where each segment begins
    int rem = (rows*cols)%processCount; // elements remaining after division among processes
    int sum = 0;                // Sum of counts. Used to calculate displacements
    int rec_buf_size = (rows*cols/processCount) + 1;
    double rec_buf[rec_buf_size];          // buffer where the received data should be stored

    sendcounts = new int[processCount];
    displs = new int[processCount];

    // calculate send counts and displacements
    for (int i = 0; i < processCount; i++) {
        sendcounts[i] = (rows*cols)/processCount;
        if (rem > 0) {
            sendcounts[i]++;
            rem--;
        }

        displs[i] = sum;
        sum += sendcounts[i];
    }

    MPI_Scatterv(matrix, sendcounts, displs, MPI_DOUBLE, &rec_buf, 100, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // print what each process received
    // printf("%d: ", rank);
    // for (int i = 0; i < sendcounts[rank]; i++) {
    //     printf("%f\t", rec_buf[i]);
    // }
    // printf("\n");


    double * processPrevBuffer = new double[sendcounts[rank]];
    double * processCurrBuffer = new double[sendcounts[rank]];
    double * nextRankValues;

    int processIndexStart = 0;
    for(int i = 0; i < rank; i++)
    {
        processIndexStart += sendcounts[i];
    }

    // std::cout << "process " << rank << " starts at " << processIndexStart << "\n";


    for(int k = 0; k < iterations; k++) 
    {
        memcpy(processPrevBuffer, rec_buf, sendcounts[rank] * sizeof(double));
        memcpy(processCurrBuffer, rec_buf, sendcounts[rank] * sizeof(double));

        // if(rank > 0)
        // {
        //     MPI_Send(rec_buf, sendcounts[rank], MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD); 
        // }

        // if(rank < processCount - 1)
        // {
        //     nextRankValues = new double[sendcounts[rank + 1]];
        //     MPI_Recv(&nextRankValues, sendcounts[rank + 1], MPI_DOUBLE, rank + 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // }

        
        for(int i = 0; i < sendcounts[rank]; i++)
        {
            if((processIndexStart + i) % cols == 0 || (processIndexStart + i + 1) % cols == 0) continue;

            c = processCurrBuffer[i];
            t = processPrevBuffer[i];
            b = matrix[processIndexStart + i];
            l = matrix[processIndexStart - 1];
            r = matrix[processIndexStart + 1];

            sleep_for(microseconds(sleep));
            rec_buf[i] = c * (1.0 - 4.0 * td / h_square) + (t + b + l + r) * (td / h_square);
        }

        memcpy(processPrevBuffer, processCurrBuffer, sendcounts[rank] * sizeof(double));
    }

    MPI_Gatherv(rec_buf, sendcounts[rank], MPI_DOUBLE, matrix, sendcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if(rank == 0)
    {
        cout << "-----  PARALLEL  -----" << endl << flush;
        printRowOrderMatrix(matrix, rows, cols);
    }
}
