#include <chrono>
#include <cstring>
#include <thread>

#include <mpi.h>

#include "solver.hpp"
#include "../matrix/matrix.hpp"

using std::memcpy;
using std::chrono::microseconds;
using std::this_thread::sleep_for;

void solveSeq(int rows, int cols, int iterations, double td, double h, int sleep, double **matrix)
{
    double c, l, r, t, b;

    double h_square = h * h;

    double *linePrevBuffer = new double[cols];
    double *lineCurrBuffer = new double[cols];

    for (int k = 0; k < iterations; k++)
    {

        memcpy(linePrevBuffer, matrix[0], cols * sizeof(double));
        for (int i = 1; i < rows - 1; i++)
        {

            memcpy(lineCurrBuffer, matrix[i], cols * sizeof(double));
            for (int j = 1; j < cols - 1; j++)
            {
                c = lineCurrBuffer[j];
                t = linePrevBuffer[j];
                b = matrix[i + 1][j];
                l = lineCurrBuffer[j - 1];
                r = lineCurrBuffer[j + 1];

                sleep_for(microseconds(sleep));
                matrix[i][j] = c * (1.0 - 4.0 * td / h_square) + (t + b + l + r) * (td / h_square);
            }
            memcpy(linePrevBuffer, lineCurrBuffer, cols * sizeof(double));
        }
    }
}

void solvePar(int rows, int cols, int iterations, double td, double h, int sleep, double **matrix)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // int *rowOrderMatrix = NULL; // same as "matrix" but rearranged into a row major order matrix

    // int previousLocalEntry; // matrix entry from [i][j-1]
    // int localEntry;         // matrix entry from [i][j]

    // if (rank == 0)
    // {
    //     rowOrderMatrix = allocateRowOrderMatrix(8, 8);                // allocate the memory for the matrix
    //     fillRowOrderMatrix(rowOrderMatrix, initialValue, rows, cols); // fills the matrix with the initial value
    //     gettimeofday(&timestamp_s, NULL);                             // sets the execution's start time
    // }

    // // Scatter the matrix values into 64 processes with 1 matrix-entry each
    // // The data will be scattered in a way that each process can know the row and col of it's value based on his rank
    // MPI_Scatter(rowOrderMatrix, 1, MPI_INT, &localEntry, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // int row = floor(rank / rows); // calculates the row of the entry in the process based on the rank
    // int col = rank - row * rows;  // calculates the col of the entry in the process based on teh rank

    // double c, l, r, t, b;

    // double h_square = h * h;

    // for (int k = 0; k < iterations; ++k)
    // {
    //     for (int i = 0; i < rows; ++i)
    //     {
    //         for (int j = 0; j < cols; ++j)
    //         {
    //             c = matrix[i][j];
    //             l = matrix[i - 1][j];
    //             r = matrix[i + 1][j];
    //             t = matrix[i][j - 1];
    //             b = matrix[i][j + 1];

    //             sleep_for(microseconds(500000));

    //             matrix[i][j] = (1 - 4 * td / h_square) * c + (td / h_square) * (l + r + t + b);
    //         }
    //     }
    // }

    // // retrieve the final data and gathers it in the process at rank = 0
    // MPI_Gather(&localEntry, 1, MPI_INT, matrix, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (0 != rank)
    {
        deallocateMatrix(rows, matrix);
    }

    sleep_for(microseconds(500000));
}
