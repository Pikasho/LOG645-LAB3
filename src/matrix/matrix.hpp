#ifndef MATRIX_HPP
#define MATRIX_HPP

double ** allocateMatrix(int rows, int cols);
void deallocateMatrix(int rows, double ** matrix);

void fillMatrix(int rows, int cols, double ** matrix);

double* allocateRowOrderMatrix(int rows, int cols);
void fillRowOrderMatrix(double *rowOrderMatrix, int rows, int cols);
int getRankfromIndex(int rows, int row, int col);
double getMatrixCell(double* matrix, int i, int j, int cols);
void printRowOrderMatrix(double* matrix, int rows, int cols);

#endif
