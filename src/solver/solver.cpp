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

    double c, l, r, t, b;

    double h_square = h * h;

    for (int k = 0; k < iterations; ++k)
    {
        for (int i = 0; i < rows; ++i)
        {
            for (int j = 0; j < cols; ++j)
            {
                c = matrix[i][j];
                l = matrix[i - 1][j];
                r = matrix[i + 1][j];
                t = matrix[i][j - 1];
                b = matrix[i][j + 1];
                matrix[i][j] = (1 - 4 * td / h_square) * c + (td / h_square) * (l + r + t + b);
                sleep_for(microseconds(500000));
            }
        }
    }

    if (0 != rank)
    {
        deallocateMatrix(rows, matrix);
    }

    //sleep_for(microseconds(500000));
}
