#pragma once

void cholesky_decomposition_basic(double* A, double* L, int n);
void printMatrix(double* a, int n);
void Cholesky_Decomposition(double* a, double* l, int n, int subMatrixSize, int rowBias, int colBias);
void Cholesky_Decomposition_Block(double* a, double* l, int n, int r);
void Cholesky_Decomposition(double* A, double* L, int n);
void multiply_matrix(double* a, double* b, double* c, int n);
void transpose_matrix(double* a, double* at, int n);
