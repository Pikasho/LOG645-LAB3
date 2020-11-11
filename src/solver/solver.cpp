/*  Solver
*
* Fonctions permettant de solve le problème de transfert de chaleur sur une plaque 2D
*
* Philippe Bolduc
* Rémi Thibault
*/

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

// Function to solve the problem sequentially
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

// Function to solve the problem in parallel using MPI
void solvePar(int rows, int cols, int iterations, double td, double h, int sleep, double * matrix) {
    int rank;           // rank of the current process
    int processCount;   // number of processes running
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &processCount);

    double c, l, r, t, b;  
    double h_square = h * h;

    int rowsPerProcessFloored = floor((double)rows / processCount); // number of rows we will give to each process 
    if(rowsPerProcessFloored == 0)  rowsPerProcessFloored = 1;      // if the number of rows per process was floored to 0, we set it to 1 manually

    // variables we will use for Scatterv 
    int sum = 0;    
    int *sendCounts;            
    int *displs;                         
    sendCounts = new int[processCount];
    displs = new int[processCount];

    // set remainingRows to the number of rows
    int remainingRows = rows;

    // Logic for setting up how we will split the rows in the processes
    // sendCounts will contain how many elements each rank will be attributed
    // displs will contain the displacements
    for (int i = 0; i < processCount; i++)
    {
        if(remainingRows <= 0)
        {
            sendCounts[i] = 0;
            continue;
        }

        // every process except the last will contained the floored rowsPerProcess
        if(i != processCount - 1)
        {
            sendCounts[i] = rowsPerProcessFloored * cols;
            remainingRows -= rowsPerProcessFloored;
            displs[i] = sum;
            sum += sendCounts[i];
        }

        // the last proceess will contain the remaining rows
        else
        {
            sendCounts[i] = remainingRows * cols;
            displs[i] = sum;
            sum += sendCounts[i];
        }
    }

    // Keep track of the the matrix index of the first element in the process
    int processIndexStart = 0;
    for(int i = 0; i < rank; i++)
    {
        processIndexStart += sendCounts[i];
    }


    int recBufSize = sendCounts[rank];  // size of the receiving buffer
    double recBuf[recBufSize];          // buffer to contain the scattered elements
    double calcBuf[sendCounts[rank]];   // buffer to contain the calculated values

    // Scatter the elements in all the processes
    MPI_Scatterv(matrix, sendCounts, displs, MPI_DOUBLE, recBuf, sendCounts[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Start of the iterations
    for(int k = 0; k < iterations; k++) 
    {
        // setup the capacity for communications with previous and next processes
        int prevProcessCapacity = sendCounts[rank];
        if(rank > 0) prevProcessCapacity = sendCounts[rank - 1];

        // setup the capacity for communications with previous and next processes
        int nextProcessCapacity = sendCounts[rank];
        if(rank < processCount - 1) nextProcessCapacity = sendCounts[rank + 1];

        double prevRowBuf[prevProcessCapacity]; // initialize the receiving buffers for the communications
        double nextRowBuf[nextProcessCapacity]; // initialize the receiving buffers for the communications
        memcpy(calcBuf, recBuf, sendCounts[rank] * sizeof(double));

        if(rank < processCount - 1)
        {
            // send to the next process
            MPI_Send(recBuf, sendCounts[rank], MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);

            // receive from the next process
            MPI_Recv(nextRowBuf, sendCounts[rank + 1], MPI_DOUBLE, rank + 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        if(rank > 0)
        {
            // send to the previous process
            MPI_Recv(prevRowBuf, sendCounts[rank - 1], MPI_DOUBLE, rank - 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            // receive from the previous process
            MPI_Send(recBuf, sendCounts[rank], MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
        }

        // for each value in the current process
        for(int i = 0; i < sendCounts[rank]; i++)
        {
            if((processIndexStart + i) % cols == 0 || (processIndexStart + i + 1) % cols == 0) continue;    // make sure we ignore the first and last cols
            if(processIndexStart + i < cols || processIndexStart + i > rows * cols - cols) continue;        // make sure we ignore the first and last rows

            // [i, j, k]
            c = recBuf[i];

            // [i - 1, j, k]
            if(i - cols >= 0)
            {
                t = recBuf[i - cols];
            }
            else
            {
                int index = i - cols + prevProcessCapacity;
                t = prevRowBuf[index];
            }

            // [i + 1, j, k]
            if(i + cols < sendCounts[rank])
            {
                b = recBuf[i + cols];
            }
            else
            {
                int index = i + cols - sendCounts[rank];
                b = nextRowBuf[index];
            }
        
            l = recBuf[i - 1];  //[i, j - 1, k]
            r = recBuf[i + 1];  // [i, j + 1, k]
            
            sleep_for(microseconds(sleep)); //sleep
            calcBuf[i] = c * (1.0 - 4.0 * td / h_square) + (t + b + l + r) * (td / h_square);   // formula
        }

        // cpy the calcBuf into recBuf
        memcpy(recBuf, calcBuf, sendCounts[rank] * sizeof(double));
    }

    // gather from all the process into the main matrix
    MPI_Gatherv(calcBuf, sendCounts[rank], MPI_DOUBLE, matrix, sendCounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // free some pointers
    free(sendCounts);
    free(displs);
}
