#include <iostream>
#include <iomanip>
#include <cstring>
//#include <lapack>
using namespace std;

#include "TriMatrix.h"


#define F77NAME(x) x##_
extern "C" {
    void F77NAME(dgemv)(char* TRANS, const int* M, const int* N,
                double* alpha, double* A, const int* LDA, double* X,
                const int& INCX, double* beta, double* C, const int& INCY);
    void F77NAME(dgetrf)(const int* M, const int* N, double* A,
                const int* lda, int* Piv, int info);
    void F77NAME(dgetrs)(char* TRANS, const int* M, double* A,
                const int* lda, int* Piv, double* b, const int* ldb, int info);
    }

TriMatrix::TriMatrix(double nu, double theta, int mSize) //constructor
{
    //setting temporary storage for TriMatrix
    double* upper;
    double* lower;
    double* diag;
    double* LHS_t;
    double* RHS_t;

    upper = new double[mSize-1];
    lower = new double[mSize-1];
    diag = new double[mSize];
    LHS_t = new double[mSize*mSize];
    RHS_t = new double[mSize*mSize];

    double d_nu= 2.0 * nu;
    double identity = 1.0;

    //constructing LHS trimatix
    //inputting values for matricies
    for (int i = 1; i < mSize-1; i++)
    {
        upper[i]= -theta * nu;
        lower[i-1] = -theta * nu;
    }

    upper[0]=0;
    lower[mSize-1] =0;

    for (int i = 1; i < (mSize-1); i++)
    {
        diag[i] = identity + theta * d_nu;
    }

    diag[0] = identity;
    diag[mSize-1] = identity;

    for (int i=0; i<mSize*mSize; i++)
    {
        if (i%mSize == 0)
        {
            LHS_t[i]=diag[i/mSize];
        }
        else if (i%mSize == 1)
        {
            LHS_t[i]=lower[(i-1)/mSize];
        }
        else if (i%mSize == mSize-1)
        {
            LHS_t[i]=upper[(i-mSize+1)/mSize];
        }
        else
        {
            LHS_t[i]=0;
        }
    };

    //pre LU decomposition
    int* Piv = new int[mSize];
    int info;
/*
    void F77NAME(dgetrf)(const int* M, const int* N, double* A,
                const int* lda, int* Piv, int info);
*/
    const int* mSize_t = &mSize;
    cout << mSize << endl;
    cout << LHS_t << endl;
    cout << Piv << endl;
    cout << info << endl;

    F77NAME(dgetrf)(mSize_t, mSize_t, LHS_t, mSize_t, Piv, info);
    cout << 1<< endl;
    //setting variables to private class
    this -> LHS = LHS_t;

    //constructing RHS trimatrix
    upper = new double[mSize-1];
    lower = new double[mSize-1];
    diag = new double[mSize];

    //inputting values for matricies
    for (int i = 1; i < mSize-1; i++)
    {
        upper[i]= (1 - theta) * nu;
        lower[i-1] = (1 - theta) * nu;
    }

    upper[0]=0;
    lower[mSize-1] =0;
    cout << 2<< endl;
    for (int i = 1; i < (mSize-1); i++)
    {
        diag[i] = identity - (1 - theta) * d_nu;
    }

    diag[0] = identity;
    diag[mSize-1] = identity;

    for (int i=0; i<mSize*mSize; i++)
    {
        if (i%mSize == 0)
        {
            RHS_t[i]=diag[i/mSize];
        }
        else if (i%mSize == 1)
        {
            RHS_t[i]=lower[(i-1)/mSize];
        }
        else if (i%mSize == mSize-1)
        {
            RHS_t[i]=upper[(i-mSize+1)/mSize];
        }
        else
        {
            RHS_t[i]=0;
        }
    }
    cout << 3<< endl;
    //setting variables to private class
    this -> RHS = RHS_t;

}

TriMatrix::~TriMatrix() //destructor of TriMatrix
{
    delete[] RHS;
    delete[] LHS;
}

//matrix multiplication function
void TriMatrix::matrixMultiplication (vector<double> &U, double ini_con_1, double ini_con_2)
{
    //temporary storage for the solution
    double* U_temp;

    U_temp = new double[U.size()];
    cout << 4<< endl;
    //calculating RHS of the equation (I-(1-theta)*nu*L)U^k
    //assigning values considering boundary condition
    for (int i=0; i < U.size(); i++)
    {
        if (i==0)
        {
            U_temp[i] = ini_con_1;
        }
        else if (i == U.size()-1)
        {
            U_temp[i] = ini_con_2;
        }
        else
        {
            U_temp[i] = U[i];
        };
    };

    int matSize = sizeof(U_temp);
    char TRANS = 'N';
    double one = 1.0;
    double zero = 0.0;
    int* Piv = new int[matSize];
    int info;
    cout << 5<< endl;
    //performing RHS first to obtain Ax=b format
    F77NAME(dgemv)(&TRANS, &matSize, &matSize, &one, RHS, &matSize, U_temp, 1, &zero, U_temp, 1);

    //then solving equation in Ax=b format
    F77NAME(dgetrs)(&TRANS, &matSize, LHS, &matSize, Piv, U_temp, &matSize, info);

    this -> U_h = U_temp;

    //replacing current values with the processed values
    for (unsigned int i=0; i < U.size(); i++)
    {
        U[i] = U_h[i];
    };
    cout << 6<< endl;


};
