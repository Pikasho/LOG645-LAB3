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

double* allocateRowOrderMatrix(int rows, int cols)
{
    double* matrix = (double*)malloc((rows * cols) * sizeof(double));
    return matrix;
}

void fillRowOrderMatrix(int *rowOrderMatrix, int value, int rows, int cols)
{
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            rowOrderMatrix[rows * i + j] = value;
        }
    }
}

int getRankfromIndex(int rows, int row, int col)
{
    return (rows * row) + col;
}

double getMatrixCell(double* matrix, int i, int j, int cols)
{
    return matrix[i * cols + j];
}

void printRowOrderMatrix(int rows, int cols, int *matrix)
{
    // prints the result matrix
    int count = 0;
    for (int i = 0; i < rows * cols; i++)
    {
        printf("%12d   ", matrix[i]);
        count++;

        if (count == rows)
        {
            printf("\n");
            count = 0;
        }
    }

    printf("\n");
}
