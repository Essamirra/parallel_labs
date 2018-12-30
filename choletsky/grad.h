#pragma once
#include <vector>

struct CRSMatrix
{
    int n; // Число строк в матрице 
    int m; // Число столбцов в матрице 
    int nz; // Число ненулевых элементов в разреженной симметричной матрице, лежащих не ниже главной диагонали 
    std::vector<double> val; // Массив значений матрицы по строкам 
    std::vector<int> colIndex; // Массив номеров столбцов 
    std::vector<int> rowPtr; // Массив индексов начала строк 
};
void SLE_Solver_CRS(CRSMatrix & A, double * b, double eps, int max_iter, double * x, int & count);
void Mult(const CRSMatrix& a, double* x, double* ax);