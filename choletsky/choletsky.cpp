// choletsky.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include "choletsky.h"
using namespace std;






void fillZeroes(double* a, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            a[i * n + j] = 0;
        }
    }
}


void cholesky_decomposition_basic(double* A, double* L, int n)
{
    fillZeroes(L, n);

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j <= i; j++)
        {
            int sum = 0;

            if (j == i)
            {
                for (int k = 0; k < j; k++)
                    sum += pow(L[j * n + k], 2);
                L[j * n + j] = sqrt(A[j * n + j] -
                    sum);
            }
            else
            {
                for (int k = 0; k < j; k++)
                    sum += (L[i * n + k] * L[j * n + k]);
                L[i * n + j] = (A[i * n + j] - sum) /
                    L[j * n + j];
            }
        }
    }
}

void printMatrix(double* a, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            cout << a[i * n + j] << " ";
        }
        cout << endl;
    }
}


using namespace std;
void Cholesky_Decomposition(double* a, double* l, int n, int subMatrixSize, int rowBias, int colBias)
{
    for (int i = 0; i < subMatrixSize; i++)
    {
        for (int j = 0; j <= i; j++)
        {
            int sum = 0;

            if (j == i)
            {
                for (int k = 0; k < j; k++)
                    sum += pow(l[(j + rowBias) * n + (k + colBias)], 2);
                l[(j + rowBias) * n + (j + colBias)] = sqrt(a[(j + rowBias) * n + (j + colBias)] -
                    sum);
            }
            else
            {
                for (int k = 0; k < j; k++)
                    sum += (l[(i + rowBias) * n + (k + colBias)] * l[(j + rowBias) * n + (k + colBias)]);
                l[(i + rowBias) * n + (j + colBias)] = (a[(i + rowBias) * n + (j + colBias)] - sum) /
                    l[(j + rowBias) * n + (j + colBias)];
            }
        }
    }
}

void solve_system(double* a, double* L, int N, int subMatrixHeight, int subMatrixWidth, int rowBiasL,
    int colBiasL, int rowBiasL2, int colBiasL2)
{
    double sum;
#pragma omp parallel private(sum)
    {
        int m, n, r;
        int coo1 = rowBiasL2 * N + colBiasL2;
        int coo2 = rowBiasL * N + colBiasL;

#pragma omp for
        for (int j = 0; j < subMatrixHeight; j++)
        {
            m = j * N + coo2;


            for (int i = 0; i < subMatrixWidth; i++)
            {
                n = coo1 + i;
                sum = 0;
#pragma simd
                for (int k = 0; k < i; k++)
                {
                    sum += L[i * N + coo1 + k] * L[k + m];
                }

                L[m + i] = (a[m + i] - sum) / L[i * N + coo1 + i];
            }
        }
    }
}

void mult_l(double* l, double* Res, int size, int subMatrixHeight, int subMatrixWidth, int rowBiasA,
    int colBiasA, int rowBiasB, int colBiasB)
{
    int u, v;
#pragma omp parallel private(u,v)
    {
#pragma omp for
        for (int i = 0; i < subMatrixHeight; i++)
        {
            u = (i + rowBiasA) * size + colBiasA;
            for (int j = 0; j < subMatrixHeight; j++)
            {

                v = (j + rowBiasB) * size + (colBiasB);
#pragma simd
                for (int t = 0; t < subMatrixWidth; t++)
                {
                    Res[i * subMatrixHeight + j] += l[u + t] * l[v + t];
                }
            }
        }
    }
}

void sub_from_a(double* a, double* b, int size, int sub_matrix_size, int row_bias, int col_bias)
{
#pragma omp parallel
    {
        int m, n;
        int coo = row_bias * size + col_bias;
#pragma omp for
        for (int i = 0; i < sub_matrix_size; i++)
        {
            m = i * sub_matrix_size;
            n = i * size + coo;
            for (int j = 0; j < sub_matrix_size; j++)
            {
                a[j + n] -= b[m + j];
            }
        }
    }
}

void Cholesky_Decomposition_Block(double* a, double* l, int n, int r)
{
    int blockCount = n / r;
    int currMatrixSize = n;

    for (int i = 0; i <= blockCount; i++)
    {
        if (currMatrixSize - r <= 0)
        {
            Cholesky_Decomposition(a, l, n, currMatrixSize, i * r, i * r);
            // cout << "Block #" << i<< ", last decomposition"<<endl;
            //  printMatrix(l, n);
            break;
        }

        Cholesky_Decomposition(a, l, n, r, i * r, i * r);
        // cout << "Block #" << i << ", decomposition, l" << endl;
        // printMatrix(l, n);

        solve_system(a, l, n, currMatrixSize - r, r, (i + 1) * r, i * r, i * r, i * r);
        //  cout << "Block #" << i << ", solve system, l" << endl;
        // printMatrix(l, n);

        double* L21L21 = new double[(currMatrixSize - r) * (currMatrixSize - r)]();
        mult_l(l, L21L21, n, currMatrixSize - r, r, (i + 1) * r, i * r, (i + 1) * r, i * r);
        // cout << "Block #" << i << ", multiplication res, l21l21" << endl;
        // printMatrix(L21L21, (currMatrixSize - r));

        sub_from_a(a, L21L21, n, currMatrixSize - r, (i + 1) * r, (i + 1) * r);
        // cout << "Block #" << i << ",diff, a" << endl;
        //  printMatrix(a, n);
        currMatrixSize = currMatrixSize - r;

        delete[] L21L21;
    }
}

void Cholesky_Decomposition(double* A, double* L, int n)
{
    Cholesky_Decomposition_Block(A, L, n, 25);
}

void GenerateTestMatrix(double* a, int n = 3)
{
    a[0] = 4;
    a[1] = 12;
    a[2] = -16;
    a[3] = 12;
    a[4] = 37;
    a[5] = -43;
    a[6] = -16;
    a[7] = -43;
    a[8] = 98;
}

void GenerateTestMatrixAnswer(double* a, int n = 3)
{
    a[0] = 2;
    a[1] = 0;
    a[2] = 0;
    a[3] = 6;
    a[4] = 1;
    a[5] = 0;
    a[6] = -8;
    a[7] = 5;
    a[8] = 3;
}


void generate_imatrix(double* a, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            a[i * n + j] = 1;
        }
    }
}

void transpose_matrix(double* a, double* at, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            at[i * n + j] = a[j * n + i];
        }
    }
}

void generate_random_matrix(double* a, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            a[i * n + j] = rand() % n;
        }
    }
}


void multiply_matrix(double* a, double* b, double* c, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            c[i * n + j] = 0;
            for (int k = 0; k < n; k++)
            {
                c[i * n + j] += (a[i * n + k] * b[k * n + j]);
            }
        }
    }
}

void decomposition_test()
{
    int n = 3;
    double* a = new double[n * n];
    double* expected = new double[n * n];
    double* actual = new double[n * n];
    GenerateTestMatrix(a, n);
    printMatrix(a, n);
    GenerateTestMatrixAnswer(expected, n);
    cholesky_decomposition_basic(a, actual, n);
    cout << "Actual" << endl;
    printMatrix(actual, n);
    cout << "Expected" << endl;
    printMatrix(expected, n);
}

void multiplication_test()
{
    int n = 3;
    double* a = new double[n * n];
    double* expected = new double[n * n]();
    double* decomposed = new double[n * n]();
    double* actual = new double[n * n]();
    GenerateTestMatrix(a, n);
    printMatrix(a, n);
    GenerateTestMatrix(expected, n);
    cholesky_decomposition_basic(a, decomposed, n);
    GenerateTestMatrixAnswer(expected, n);
    mult_l(decomposed, actual, n, n, n, 0, 0, 0, 0);
    cout << "Actual" << endl;
    printMatrix(actual, n);
    cout << "Expected" << endl;
    printMatrix(expected, n);
}

void multiplication_test_step()
{
    int n = 3;
    double* a = new double[n * n];
    double* expected = new double[3 * 3]();
    double* decomposed = new double[n * n]();
    double* decomposed2 = new double[n * n]();
    double* actual = new double[3 * 3]();
    GenerateTestMatrix(a, n);
    printMatrix(a, n);
    GenerateTestMatrix(expected, n);
    cholesky_decomposition_basic(a, decomposed, n);
    GenerateTestMatrixAnswer(decomposed2, n);
    mult_l(decomposed, actual, n, 3, 2, 0, 1, 0, 1);
    cout << "Decomposed" << endl;
    printMatrix(decomposed2, 3);
    cout << "Actual" << endl;
    printMatrix(actual, 3);
    cout << "Expected" << endl;
    printMatrix(expected, 3);
}

void generate_test_matrix5(double* a, int n)
{
    a[0] = 1;
    a[1] = 1;
    a[2] = 2;
    a[3] = 4;
    a[4] = 7;
    a[5] = 1;
    a[6] = 5;
    a[7] = 8;
    a[8] = 14;
    a[9] = 23;
    a[10] = 2;
    a[11] = 8;
    a[12] = 22;
    a[13] = 41;
    a[14] = 65;
    a[15] = 4;
    a[16] = 14;
    a[17] = 41;
    a[18] = 93;
    a[19] = 162;
    a[20] = 7;
    a[21] = 23;
    a[22] = 65;
    a[23] = 162;
    a[24] = 319;
}


void generate_test_matrix_answer5(double* a, int n)
{
    a[0] = 1;
    a[1] = 0;
    a[2] = 0;
    a[3] = 0;
    a[4] = 0;
    a[5] = 1;
    a[6] = 2;
    a[7] = 0;
    a[8] = 0;
    a[9] = 0;
    a[10] = 2;
    a[11] = 3;
    a[12] = 3;
    a[13] = 0;
    a[14] = 0;
    a[15] = 4;
    a[16] = 5;
    a[17] = 6;
    a[18] = 4;
    a[19] = 0;
    a[20] = 7;
    a[21] = 8;
    a[22] = 9;
    a[23] = 10;
    a[24] = 5;
}

void decomposition_block_test()
{
    int n = 5;
    double* a = new double[n * n];
    double* expected = new double[n * n]();
    double* actual = new double[n * n]();
    generate_test_matrix5(a, n);
    printMatrix(a, n);
    generate_test_matrix_answer5(expected, n);
    Cholesky_Decomposition_Block(a, actual, n, 2);
    cout << "Actual" << endl;
    printMatrix(actual, n);
    cout << "Expected" << endl;
    printMatrix(expected, n);
}

void decomposition_block_test2()
{
    int n = 3;
    double* a = new double[n * n];
    double* expected = new double[n * n]();
    double* actual = new double[n * n]();
    GenerateTestMatrix(a, n);
    printMatrix(a, n);
    GenerateTestMatrixAnswer(expected, n);
    Cholesky_Decomposition_Block(a, actual, n, 2);
    cout << "Actual" << endl;
    printMatrix(actual, n);
    cout << "Expected" << endl;
    printMatrix(expected, n);
}


