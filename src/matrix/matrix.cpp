/*  Matrix
*
* Fonctions utilitaires pour des matrices stockées dans un tableau 2d ou 1d
*
* Philippe Bolduc
* Rémi Thibault
*/

#include "matrix.hpp"
#include <stdlib.h>
#include <stdio.h>

double ** allocateMatrix(int rows, int cols) {
    double ** matrix = new double*[rows];

    for(int i = 0; i < rows; i++) {
        matrix[i] = new double[cols];
    }

    return matrix;
}

void deallocateMatrix(int rows, double ** matrix) {
    for(int i = 0; i < rows; i++) {
        delete(matrix[i]);
        matrix[i] = nullptr;
    }

    delete(matrix);
    *matrix = nullptr;
}

void fillMatrix(int rows, int cols, double ** matrix) {
     for(int row = 0; row < rows; row++) {
        for(int col = 0; col < cols; col++) {
            matrix[row][col] = row * (rows - row - 1) * col * (cols - col - 1);
        }
    }
}

// **********************************
// NEW FUNCTIONS FOR ROW ORDER MATRIX
// **********************************

// Allocate the memory for a row order matrix
double* allocateRowOrderMatrix(int rows, int cols)
{
    double* matrix = (double*)malloc((rows * cols) * sizeof(double));
    return matrix;
}

// Fill a row order matrix with the correct values for the problem
void fillRowOrderMatrix(double *rowOrderMatrix, int rows, int cols)
{
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            rowOrderMatrix[i * cols + j] = i * (rows - i - 1) * j * (cols - j - 1);
        }
    }
}

// returns the matric cell at [i, j] in a row order matrix
double getMatrixCell(double* matrix, int i, int j, int cols)
{
    return matrix[i * cols + j];
}