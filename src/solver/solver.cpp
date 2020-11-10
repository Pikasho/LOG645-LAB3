#include <chrono>
#include <cstring>
#include <thread>
#include <iostream>
#include <math.h>

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

    int rowsPerProcess = ceil((double)rows / processCount);
    int *sendcounts;            
    int *displs;                
    int rem = (rows*cols)%processCount; 
    int sum = 0;                


    int rec_buf_size = (rowsPerProcess * cols);
    double rec_buf[rec_buf_size];          

    int remainingRows = rows;

    sendcounts = new int[processCount];
    displs = new int[processCount];

    // calculate send counts and displacements
    for (int i = 0; i < processCount; i++) {
        sendcounts[i] = 0;

        if(remainingRows > rowsPerProcess)
        {
            sendcounts[i] = rowsPerProcess * cols;
            remainingRows = remainingRows - rowsPerProcess;
            displs[i] = sum;
            sum += sendcounts[i];
        }
        else
        {
            sendcounts[i] = remainingRows * cols;
            displs[i] = sum;
            sum += sendcounts[i];
            break;
        }
    }



    double * processPrevBuffer = new double[sendcounts[rank]];
    double * processCurrBuffer = new double[sendcounts[rank]];

    int processIndexStart = 0;
    for(int i = 0; i < rank; i++)
    {
        processIndexStart += sendcounts[i];
    }

    // std::cout << "process " << rank << " starts at " << processIndexStart << "\n";


    for(int k = 0; k < iterations; k++) 
    {
        MPI_Scatterv(matrix, sendcounts, displs, MPI_DOUBLE, rec_buf, 100, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        for(int i = 0; i < sendcounts[rank]; i++)
        {
            if((processIndexStart + i) % cols == 0 || (processIndexStart + i + 1) % cols == 0) continue;
            if(processIndexStart + i < cols || processIndexStart + i > rows * cols - cols) continue;

            c = matrix[processIndexStart + i];
            t = matrix[processIndexStart + i - cols];
            b = matrix[processIndexStart + i + cols];
            l = matrix[processIndexStart + i - 1];
            r = matrix[processIndexStart + i + 1];

            sleep_for(microseconds(sleep));
            rec_buf[i] = c * (1.0 - 4.0 * td / h_square) + (t + b + l + r) * (td / h_square);
        }
        MPI_Gatherv(rec_buf, sendcounts[rank], MPI_DOUBLE, matrix, sendcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }


    if(rank == 0)
    {
        cout << "-----  PARALLEL  -----" << endl << flush;
        printRowOrderMatrix(matrix, rows, cols);
    }
}
