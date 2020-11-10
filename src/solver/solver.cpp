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

                // if(k == 0)
                // {
                //     std::cout << t << "    \n";
                // }

                // if(k == 1)
                // {
                //     std::cout << b << "    \n";
                // }

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

    int rowsPerProcess = floor((double)rows / processCount);         
    int sum = 0;    
    int *sendcounts;            
    int *displs;                         

    int remainingRows = rows;

    sendcounts = new int[processCount];
    displs = new int[processCount];

    // calculate send counts and displacements
    for (int i = 0; i < processCount; i++) {
        if(i != processCount -1)
        {
            sendcounts[i] = rowsPerProcess * cols;
            remainingRows -= rowsPerProcess;
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

    // for(int i = 0; i < processCount; i++)
    // {
    //     std::cout << "\n rank: " << rank << " " << sendcounts[i];
    // }


    int processIndexStart = 0;
    for(int i = 0; i < rank; i++)
    {
        processIndexStart += sendcounts[i];
    }


    int rec_buf_size = sendcounts[rank];
    double rec_buf[rec_buf_size];    
    double calcBuf[sendcounts[rank]];

    for(int k = 0; k < iterations; k++) 
    {
        MPI_Scatterv(matrix, sendcounts, displs, MPI_DOUBLE, rec_buf, sendcounts[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

        int prevProcessCapacity = 0;
        if(rank > 0) prevProcessCapacity = sendcounts[rank - 1];

        int nextProcessCapacity = 0;
        if(rank < processCount - 1) nextProcessCapacity = sendcounts[rank + 1];

        double prevRowBuf[prevProcessCapacity];

        double nextRowBuf[nextProcessCapacity];

        memcpy(calcBuf, rec_buf, sendcounts[rank] * sizeof(double));

        // send the rows in this rank to previous rank 
        if(rank < processCount - 1)
        {
            MPI_Send(rec_buf, sendcounts[rank], MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);

            MPI_Recv(nextRowBuf, sendcounts[rank + 1], MPI_DOUBLE, rank + 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        if(rank > 0)
        {
            MPI_Recv(prevRowBuf, sendcounts[rank - 1], MPI_DOUBLE, rank - 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            MPI_Send(rec_buf, sendcounts[rank], MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
        }

        for(int i = 0; i < sendcounts[rank]; i++)
        {
            if((processIndexStart + i) % cols == 0 || (processIndexStart + i + 1) % cols == 0) continue;
            if(processIndexStart + i < cols || processIndexStart + i > rows * cols - cols) continue;

            c = rec_buf[i];
            

            if(i - cols > 0)
            {
                t = rec_buf[i - cols];
            }
            else
            {
                t = prevRowBuf[i];
            }

            if(i + cols < sendcounts[rank])
            {
                b = rec_buf[i + cols];
            }
            else
            {
                b = nextRowBuf[i];
            }

            // std::cout << "rank: " << rank << "   " << sendcounts[rank] << "\n";
            

            l = rec_buf[i - 1];
            r = rec_buf[i + 1];
            


            sleep_for(microseconds(sleep));
            calcBuf[i] = c * (1.0 - 4.0 * td / h_square) + (t + b + l + r) * (td / h_square);
        }

        MPI_Gatherv(calcBuf, sendcounts[rank], MPI_DOUBLE, matrix, sendcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }



    if(rank == 0)
    {
        cout << "-----  PARALLEL  -----" << endl << flush;
        printRowOrderMatrix(matrix, rows, cols);
    }
}
