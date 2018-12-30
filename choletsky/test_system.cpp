#include "stdafx.h"
#include <cstdlib>
#include "choletsky.h"
#include <omp.h>
#include <iostream>
#include <ctime>
#include "grad.h"
#include <set>
#include <cassert>
#include <fstream>
#include <string>
const double MAX_VAL = 100.0;
void generateRandomSymmetricMatrix(double * a, int n)
{
    double * l = new double[n*n]();
    double *lt = new double[n*n]();
    srand(time(NULL));
    for (int i = 0; i <n; i++)
    {
       for(int j = 0; j <= i; j++)
       {
           l[i*n + j] = rand() % 10000000000 + 1;

       }
    }

    transpose_matrix(l, lt, n);
    multiply_matrix(l, lt, a, n);
    delete[] l;
    delete[] lt;
}

void choletsky_seq_test()
{
    int n = 10;
    int maxn = 1000;
    omp_set_num_threads(1);
    for(int i = n ; i <= maxn; i+=10)
    {
        double full_time = 0;
        for(int j = 0; j < 5; j++)
        {
            double * a = new double[i*i]();
            double * l = new double[i*i]();
            generateRandomSymmetricMatrix(a, i);
            double start = omp_get_wtime();
            Cholesky_Decomposition(a, l, i);
            double time = omp_get_wtime() - start;
            full_time += time;
            delete[] a;
            delete[] l;
        }
        
        std::cout << i<< ";"  <<  full_time / 5.0 * 1000 << std::endl;
   
    }
}

double next()
{
    return ((double)rand() / (double)RAND_MAX);
}
void InitializeMatrix(int N, int NZ, CRSMatrix &mtx)
{
    mtx.n = N;
    mtx.nz = NZ;
    mtx.val.resize(mtx.nz);
    mtx.colIndex.resize(mtx.nz);
    mtx.rowPtr.resize(mtx.n + 1);
}

void CreateCRSMatrix(CRSMatrix& crsMatrix, int n, int nz)
{
    crsMatrix.n = n;
    crsMatrix.m = n;
    crsMatrix.nz = nz;
    crsMatrix.val.resize(crsMatrix.nz);
    crsMatrix.colIndex.resize(crsMatrix.nz);
    crsMatrix.rowPtr.resize(crsMatrix.n + 1);
}

void InitUpperTriangleFullCRSPositiveDefiniteMatrix(CRSMatrix& crsMatrix,
    int n, int countElementsInRow)
{
    /* Creating CRS matrix */
    CRSMatrix tmpCRS;
    int tmpnNz = n * countElementsInRow;
    CreateCRSMatrix(tmpCRS, n, tmpnNz);

    /* Initializing CRS matrix */
    /* Filling colIndex vector (diagonal element must be specified) */
    for (int i = 0; i < n; ++i)
    {
        tmpCRS.colIndex[i * countElementsInRow + countElementsInRow - 1] = i; // diagonal element 
        for (int j = 0; j < countElementsInRow - 1; ++j)
        {
            int col = n / countElementsInRow * j;
            if (col != i)
            {
                tmpCRS.colIndex[i * countElementsInRow + j] = col;
            }
            else
            {
                tmpCRS.colIndex[i * countElementsInRow + j] = col + 1;
            }
        }
    }

    /* Filling non-zero values */
    int intervalSize = 10;
    for (int i = 0; i < tmpCRS.nz; ++i)
    {
        tmpCRS.val[i] = rand() % intervalSize - intervalSize / 2 +
            (rand() % 9 * 0.1);
        if (tmpCRS.val[i] == 0.0)
        {
            tmpCRS.val[i] += 1;
        }
    }

    /* Filling rows vector */
    tmpCRS.rowPtr[0] = 0;
    for (int i = 1; i < tmpCRS.n + 1; ++i)
    {
        tmpCRS.rowPtr[i] = tmpCRS.rowPtr[i - 1] + countElementsInRow;
    }

    /* Add diagonal predominance */
    for (int i = 0; i < n; ++i)
    {
        for (int j = tmpCRS.rowPtr[i]; j < tmpCRS.rowPtr[i + 1]; ++j)
        {
            if (tmpCRS.colIndex[j] == i)
            {
                tmpCRS.val[j] = (intervalSize / 2) * n + rand() % intervalSize;
            }
        }
    }

    /* Extract upper diagonal from tmpCRS */
    int nz = 0;
    crsMatrix.n = n;
    crsMatrix.m = n;
    crsMatrix.rowPtr.resize(n + 1);
    for (int i = 0; i < n; ++i)
    {
        crsMatrix.rowPtr[i] = crsMatrix.rowPtr[0] + nz;
        for (int j = tmpCRS.rowPtr[i]; j < tmpCRS.rowPtr[i + 1]; ++j)
        {
            if (tmpCRS.colIndex[j] >= i)
            {
                ++nz;
                crsMatrix.val.push_back(tmpCRS.val[j]);
                crsMatrix.colIndex.push_back(tmpCRS.colIndex[j]);
            }
        }
    }
    crsMatrix.nz = nz;
    crsMatrix.rowPtr[n] = nz;
   
}
void loadFromMtxFile(const std::string& fileName, CRSMatrix& a)
{
    std::ifstream infile(fileName);

    int n, m, nz;
    infile >> n >> m >> nz;
    assert(n == m);

    a.n = n;
    a.m = m;
    a.nz = nz;

    //int colind = 0;
    int row, col, prevRow = -1;
    double val;
    while (infile >> col >> row >> val) {
        if (prevRow < row - 1) {
            int cur = static_cast<int>(a.val.size());
            a.rowPtr.push_back(cur);
            prevRow = row - 1;
        }
        a.val.push_back(val);
        a.colIndex.push_back(col - 1);
    }
    a.rowPtr.push_back(nz);
}
CRSMatrix Generate(size_t n, size_t nz)
{
    CRSMatrix result;
    result.n = n;
    result.m = n;
    result.nz = nz;
    srand(time(NULL));
    std::set<std::pair<size_t, size_t>> pairs;
    for (int i = 0; i < n; ++i)
    {
        pairs.emplace(i, i);
    }
    while (pairs.size() <= (nz + n) / 2)
    {
      
        int i = rand() % (n/2);
        int j = n/2 + rand() % (n/2);
        if (i < j)
        {
            pairs.emplace(i, j);
        }
       
    }

    int idx = 0;
    int cur_row = -1;

    for (auto& pair : pairs)
    {
        while (cur_row != pair.first)
        {
            result.rowPtr.push_back(idx);
            ++cur_row;
        }

        result.colIndex.push_back(pair.second);
        result.val.push_back((10.0 * rand() / RAND_MAX) - 5.0);
        if (pair.first == pair.second)
        {
            result.val.back() += 5.00 * ((nz / n) + 3);
        }
        ++idx;
    }
    while (result.rowPtr.size() != (n + 1))
    {
        result.rowPtr.push_back(idx);
    }

    return result;
}
void GenerateRegularCRS(int seed, int N, int cntInRow, CRSMatrix& mtx)
{
    int i, j, k, f, c;

    srand(seed);

    int notNull = cntInRow * N;
    InitializeMatrix(N, notNull, mtx);

    for (i = 0; i < N; i++)
    {
        // Формируем номера столбцов в строке i
        for (j = 0; j < cntInRow; j++)
        {
            do
            {
                mtx.colIndex[i * cntInRow + j] = rand() % N;
                f = 0;
                for (k = 0; k < j; k++)
                    if (mtx.colIndex[i * cntInRow + j] == mtx.colIndex[i * cntInRow + k])
                        f = 1;
            } while (f == 1);
        }
        // Сортируем номера столцов в строке i
        for (j = 0; j < cntInRow - 1; j++)
            for (k = 0; k < cntInRow - 1; k++)
                if (mtx.colIndex[i * cntInRow + k] > mtx.colIndex[i * cntInRow + k + 1])
                {
                    int tmp = mtx.colIndex[i * cntInRow + k];
                    mtx.colIndex[i * cntInRow + k] = mtx.colIndex[i * cntInRow + k + 1];
                    mtx.colIndex[i * cntInRow + k + 1] = tmp;
                }
    }

    // Заполняем массив значений
    for (i = 0; i < cntInRow * N; i++)
        mtx.val[i] = next() * MAX_VAL;

    // Заполняем массив индексов строк
    c = 0;
    for (i = 0; i <= N; i++)
    {
        mtx.rowPtr[i] = c;
        c += cntInRow;
    }
}

void GenerateVector(int seed, int N, double * vec)
{
    srand(seed);
    
    for (int i = 0; i < N; i++)
    {
        vec[i] = next() * MAX_VAL;
    }
}

void grad_test()
{
    int n = 100;
    int maxn = 10000000;
    int cnt = 5;
    double eps = 0.001;
    int iter = 10000;
    omp_set_num_threads(2);        
        
        int runCount = 3;
        for (int i = n; i <= maxn; i *= 10, cnt *= 10) {
            double full_time = 0;
            int avg_count = 0;
            for (int j = 0; j < runCount; j++)
            {
               // std::cout << "generate matrix" << std::endl;
                CRSMatrix a = Generate(i, cnt);
               // std::cout << " matrix generated" << std::endl;
                double* v = new double[i];
                double* x = new double[i];
                GenerateVector(time(NULL), i, v);
              //  std::cout << " vector generated" << std::endl;
                int realCount = 0;
                double start = omp_get_wtime();
                SLE_Solver_CRS(a, v, eps, iter, x, realCount);
                double time = omp_get_wtime() - start;
                full_time += time;
                avg_count += realCount;

                delete[] v;
                delete[] x;
               // std::cout << j << std::endl;
            }
            std::cout << i << ";" << full_time / runCount << ";" << avg_count / runCount << std::endl;
        }
        

   
}
void test(std::string fileName, double eps)
{
    CRSMatrix a;
    loadFromMtxFile(fileName, a);

    auto* b = new double[a.n];
    auto* x = new double[a.n];

    for (int i = 0; i < a.n; ++i) {
        x[i] = rand() % 100 - 50;
    }
    Mult(a, x, b);

    for (int num_threads = 1; num_threads < 5; num_threads *= 2) {
        omp_set_num_threads(num_threads);

        memset(x, 0, a.n * sizeof(double));

        int iter;
        double start = omp_get_wtime();
        SLE_Solver_CRS(a, b, eps, 100000, x, iter);
        double time = omp_get_wtime() - start;
        std::cout << fileName << ", " << eps << ", " << iter << ", " << num_threads << ", " << time << std::endl;
    }
}

int main()
{
    omp_set_num_threads(1);

    std::cout << "fileName, eps, iter, threads, time" << std::endl;
     /*test("bcsstk01.mtx", 0.0001);
      test("bcsstk01.mtx", 0.0000001);
      test("bcsstk05.mtx", 0.0001);
      test("bcsstk05.mtx", 0.0000001);
      test("bcsstk10.mtx", 0.0001);
      test("bcsstk10.mtx", 0.0000001);
     test("bcsstk13.mtx", 0.0001);
      test("bcsstk13.mtx", 0.0000001);
      test("bcsstk16.mtx", 0.0001);
      test("bcsstk16.mtx", 0.0000001);*/
    test("tmt_sym.mtx", 0.0001);

    return 0;
   
}