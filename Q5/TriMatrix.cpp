#include <iostream>
#include <iomanip>
#include <cstring>

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

//splitted TriMatrix function in Q4 for two processors
TriMatrix::TriMatrix(double nu, double theta, int mSize, int my_rank) //constructor
{
    double d_nu= 2.0 * nu;
    double identity = 1.0;
    int* Piv = new int[mSize];
    int info;

    // constructing the half matricies for each processor
    if (my_rank == 0) //for processor 0
    {
        double* upper0;
        double* lower0;
        double* diag0;
        double* LHS_t0;
        double* RHS_t0;

        upper0 = new double[mSize-1];
        lower0 = new double[mSize-1];
        diag0 = new double[mSize];
        LHS_t0 = new double[mSize*mSize];
        RHS_t0 = new double[mSize*mSize];

        //constructing LHS trimatix
        //inputting values for matricies
        for (int i = 0; i < mSize-1; i++)
        {
            upper0[i]= -theta * nu;
            lower0[i] = -theta * nu;
        }

        upper0[0]=0;

        for (int i = 1; i < mSize; i++)
        {
            diag0[i] = identity + theta * d_nu;
        }

        diag0[0] = identity;
        //constructing full matrix in order to use LAPACK and BLAS
        for (int i=0; i<mSize*mSize; i++)
        {
            if (i%mSize == 0)
            {
                LHS_t0[i]=diag0[i/mSize];
            }
            else if (i%mSize == 1)
            {
                LHS_t0[i]=lower0[(i-1)/mSize];
            }
            else if (i%mSize == mSize-1)
            {
                LHS_t0[i]=upper0[(i-mSize+1)/mSize];
            }
            else
            {
                LHS_t0[i]=0;
            }
        };

        //pre LU decomposition using LAPACK
        F77NAME(dgetrf)(&mSize, &mSize, LHS_t0, &mSize, Piv, info);

        //setting variables to private class
        this -> LHS0 = LHS_t0;

        //constructing RHS trimatrix
        upper0 = new double[mSize-1];
        lower0 = new double[mSize-1];
        diag0 = new double[mSize];

        //inputting values for matricies
        for (int i = 0; i < mSize-1; i++)
        {
            upper0[i]= (1 - theta) * nu;
            lower0[i] = (1 - theta) * nu;
        }

        upper0[0]=0;

        for (int i = 1; i < mSize; i++)
        {
            diag0[i] = identity - (1 - theta) * d_nu;
        }

        diag0[0] = identity;
        //constructing full matrix in order to use LAPACK and BLAS
        for (int i=0; i<mSize*mSize; i++)
        {
            if (i%mSize == 0)
            {
                RHS_t0[i]=diag0[i/mSize];
            }
            else if (i%mSize == 1)
            {
                RHS_t0[i]=lower0[(i-1)/mSize];
            }
            else if (i%mSize == mSize-1)
            {
                RHS_t0[i]=upper0[(i-mSize+1)/mSize];
            }
            else
            {
                RHS_t0[i]=0;
            }
        }

        //setting variables to private class
        this -> RHS0 = RHS_t0;
    }
    else //for processor 1, same steps as above
    {
        double* upper1;
        double* lower1;
        double* diag1;
        double* LHS_t1;
        double* RHS_t1;

        upper1 = new double[mSize-1];
        lower1 = new double[mSize-1];
        diag1 = new double[mSize];
        LHS_t1 = new double[mSize*mSize];
        RHS_t1 = new double[mSize*mSize];

        //constructing LHS trimatix
        //inputting values for matricies
        for (int i = 0; i < mSize-1; i++)
        {
            upper1[i]= -theta * nu;
            lower1[i] = -theta * nu;
        }

        lower1[mSize-2]=0;

        for (int i = 0; i < mSize-1; i++)
        {
            diag1[i] = identity + theta * d_nu;
        }

        diag1[mSize-1] = identity;

        for (int i=0; i<mSize*mSize; i++)
        {
            if (i%mSize == 0)
            {
                LHS_t1[i]=diag1[i/mSize];
            }
            else if (i%mSize == 1)
            {
                LHS_t1[i]=lower1[(i-1)/mSize];
            }
            else if (i%mSize == mSize-1)
            {
                LHS_t1[i]=upper1[(i-mSize+1)/mSize];
            }
            else
            {
                LHS_t1[i]=0;
            }
        };

        //pre LU decomposition using LAPACK
        F77NAME(dgetrf)(&mSize, &mSize, LHS_t1, &mSize, Piv, info);

        //setting variables to private class
        this -> LHS1 = LHS_t1;

        //constructing RHS trimatrix
        upper1 = new double[mSize-1];
        lower1 = new double[mSize-1];
        diag1 = new double[mSize];

        //inputting values for matricies
        for (int i = 0; i < mSize-1; i++)
        {
            upper1[i] = (1 - theta) * nu;
            lower1[i] = (1 - theta) * nu;
        }

        lower1[mSize-2]=0;

        for (int i = 1; i < mSize; i++)
        {
            diag1[i] = identity - (1 - theta) * d_nu;
        }

        diag1[0] = identity;
        //constructing full matrix in order to use LAPACK and BLAS
        for (int i=0; i<mSize*mSize; i++)
        {
            if (i%mSize == 0)
            {
                RHS_t1[i]=diag1[i/mSize];
            }
            else if (i%mSize == 1)
            {
                RHS_t1[i]=lower1[(i-1)/mSize];
            }
            else if (i%mSize == mSize-1)
            {
                RHS_t1[i]=upper1[(i-mSize+1)/mSize];
            }
            else
            {
                RHS_t1[i]=0;
            }
        }

        //setting variables to private class
        this -> RHS1 = RHS_t1;
    }
}

TriMatrix::~TriMatrix() //destructor of TriMatrix
{
    delete[] U_h0;
    delete[] RHS0;
    delete[] LHS0;
    delete[] U_h1;
    delete[] RHS1;
    delete[] LHS1;
}

//matrix multiplication function for two processors
void TriMatrix::matrixMultiplication (vector<double> &U, double ini_con, int my_rank)
{
    int matSize = U.size();
    char TRANS = 'N';
    double one = 1.0;
    double zero = 0.0;
    int* Piv = new int[matSize];
    int info;

    if (my_rank==0) //processor 0
    {
        //temporary storage for the solution
        double* U_temp0;

        U_temp0 = new double[U.size()];

        //calculating RHS of the equation (I-(1-theta)*nu*L)U^k
        //assigning values considering boundary condition
        for (unsigned int i=0; i < U.size(); i++)
        {
            if (i==0)
            {
                U_temp0[i] = ini_con;
            }
            else
            {
                U_temp0[i] = U[i];
            };
        };

        //performing RHS(Bb) first to obtain Ax=b format from Ax=Bb format
        F77NAME(dgemv)(&TRANS, &matSize, &matSize, &one, RHS0, &matSize, U_temp0, 1, &zero, U_temp0, 1);

        //then solving x, equation in Ax=b format
        F77NAME(dgetrs)(&TRANS, &matSize, LHS0, &matSize, Piv, U_temp0, &matSize, info);

        this -> U_h0 = U_temp0;

        //replacing current values with the processed values
        for (unsigned int i=0; i < U.size(); i++)
        {
            U[i] = U_h0[i];
        };
    }
    else //processor 1
    {
        //temporary storage for the solution
        double* U_temp1;

        U_temp1 = new double[U.size()];

        //calculating RHS of the equation (I-(1-theta)*nu*L)U^k
        //assigning values considering boundary condition
        for (unsigned int i=0; i < U.size(); i++)
        {
            if (i==U.size()-1)
            {
                U_temp1[i] = ini_con;
            }
            else
            {
                U_temp1[i] = U[i];
            };
        };

        //performing RHS(Bb) first to obtain Ax=b format from Ax=Bb format
        F77NAME(dgemv)(&TRANS, &matSize, &matSize, &one, RHS0, &matSize, U_temp1, 1, &zero, U_temp1, 1);

        //then solving x, equation in Ax=b format
        F77NAME(dgetrs)(&TRANS, &matSize, LHS0, &matSize, Piv, U_temp1, &matSize, info);

        this -> U_h1 = U_temp1;

        //replacing current values with the processed values
        for (unsigned int i=0; i < U.size(); i++)
        {
            U[i] = U_h1[i];
        };
    }
};
