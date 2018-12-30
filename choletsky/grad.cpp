#include "stdafx.h"
#include "grad.h"
#include <iostream>


bool isEnough(double * x, double *x1, int n, double eps)
{
    double sum = 0;
    #pragma omp parallel for reduction(+ :sum)
    for (int i = 0; i < n; i++)
    {
        sum += (x[i] - x1[i])*(x[i] - x1[i]);
    }
    return  sqrt(sum) < eps;

}

double ScalarMult(double * a, double * b, int n)
{
    
    double sum = 0;
    #pragma omp parallel for reduction(+ :sum)
    for (int i = 0; i < n;i++)
    {
        sum += (a[i] * b[i]);
    }
    return sum;

}

double CalculateAlpha(double * ap, double* r, double* p, int n)
{
    double r1 = ScalarMult(r, r, n);
    double app = ScalarMult(ap, p, n);
   
   
    
    return  r1 / app;
}


void Mult(const CRSMatrix& a, double* x, double* ax)
{
    for (int i = 0; i < a.n; i++)
        ax[i] = 0;
    
    for (int i = 0; i < a.n; i++)
    {
        #
        for (int j = a.rowPtr[i] + 1; j < a.rowPtr[i + 1];
            j++)
        {

            ax[i] += a.val[j] * x[a.colIndex[j]];

            ax[a.colIndex[j]] += a.val[j] * x[i];
        }

        ax[i] += a.val[a.rowPtr[i]] *
            x[a.colIndex[a.rowPtr[i]]];
    }

    

}

void Init(double* p, double* r, double* b, const CRSMatrix& a, double* x, int n)
{
    double * ax = new double[n]();
    
    Mult(a, x, ax);
    //    for (int i = 0; i < 5; i++)
    //       std::cout << "ax i = " << ax[i]<< std::endl;
#pragma omp parallel for
    for(int i = 0; i < n; i++)
    {
         r[i] = b[i] - ax[i];
         p[i] = b[i] - ax[i];
    }   
   // for (int i = 0; i < 5; i++)
  //      std::cout << "r i = " << r[i] << std::endl;
    delete[] ax;
}
void swap(double * &a, double * &b)
{
    double * t;
    t = a;
    a = b;
    b = t;

}


double CalculateBeta(double* r1, double* r, int n)
{
    return ScalarMult(r1, r1, n) / ScalarMult(r, r, n);
}

void SLE_Solver_CRS(CRSMatrix & A, double * b, double eps, int max_iter, double * x, int & count)
{
    int n = A.n;
    double * p = new double[n]();
    double * r = new double[n]();
    double * x1 = new double[n]();
    double * r1 = new double[n]();
    double * p1 = new double[n]();
    double * Ap = new double[n]();    
    count = 0;
   
    Init(r, p, b, A, x, n);
    while (count <= max_iter)
    {
       
        Mult(A, p, Ap);
     //   for (int i = 0; i < 5; i++)
     //       std::cout << "Ap i = " << Ap[i] << std::endl;
        
        double alpha = CalculateAlpha(Ap, r, p,n);
      //  std::cout << "Alpha = " << alpha << std::endl;

        for(int i = 0; i < n; i++)
        {
            x1[i] = x[i] + alpha * p[i];
        }
       
        
        if (isEnough(x, x1, n, eps))
        {
            for (int i = 0; i < n; i++)
                x[i] = x1[i];
           
           return;
        }            
        
        #pragma omp parallel for
        for(int i = 0; i < n; i++)
        {
            r1[i] = r[i] - alpha * Ap[i];
        }
        
        double beta = CalculateBeta(r1, r, n);
        #pragma omp parallel for
        for(int i = 0; i < n; i++)
        {
            p1[i] = r1[i] + beta * p[i];
        }
       
        for (int i = 0; i < n; i++)
            x[i] = x1[i];
        swap(r1, r);
        swap(p1, p);
           
        count++;
    }    
    delete[] p;
    delete[] r;
    delete[] x1;
    delete[] Ap;
    delete[] r1;
    delete[] p1;


}

void fillTestA(CRSMatrix & a)
{

    a.n = 5;
    a.m = a.n;
    a.nz = 9;
    a.val.resize(a.nz);
    a.colIndex.resize(a.nz);
    a.rowPtr.resize(a.n + 1);

    a.val[0] = 4;
    a.val[1] = 1;
    a.val[2] = 2;
    a.val[3] = 0.5;
    a.val[4] = 2;
    a.val[5] = 0.5;
    a.val[6] = 3;
    a.val[7] = 0.625;
    a.val[8] = 16;

    a.colIndex[0] = 0;
    a.colIndex[1] = 1;
    a.colIndex[2] = 2;
    a.colIndex[3] = 3;
    a.colIndex[4] = 4;
    a.colIndex[5] = 1;
    a.colIndex[6] = 2;
    a.colIndex[7] = 3;
    a.colIndex[8] = 4;

    a.rowPtr[0] = 0;
    a.rowPtr[1] = 5;
    a.rowPtr[2] = 6;
    a.rowPtr[3] = 7;
    a.rowPtr[4] = 8;
    a.rowPtr[5] = 9;

    
}

void fillTestB(double * b)
{
    
}

//int main()
//{
//    CRSMatrix a;
//    fillTestA(a);
//    double* b = new double[a.n];
//    b[0] = 17;
//    b[1] = 3;
//    b[2] = 7;
//    b[3] = 6;
//    b[4] = 12;
//   
//    
//    double * x = new double[a.n]();
//    x[0] = 1;
//    
//    int count;
//    SLE_Solver_CRS(a, b, 0.00000001, 1000, x, count);
//    for (int i = 0; i < 5; i++)
//        std::cout << "x i = " << x[i];
//    return 0;
//}