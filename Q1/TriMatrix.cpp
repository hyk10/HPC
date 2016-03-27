#include <iostream>
#include <iomanip>
#include <cstring>

using namespace std;

#include "TriMatrix.h"

TriMatrix::TriMatrix(double nu, int mSize) //constructor
{
    //setting temporary storage for TriMatrix
    double* upper;
    double* lower;
    double* diag;
    upper = new double[mSize-1];
    lower = new double[mSize-1];
    diag = new double[mSize];

    double d_nu= 2.0 * nu;
    double identity = 1.0;

    //inputting values for matricies
    for (int i = 1; i < mSize-1; i++)
    {
        upper[i]= nu;
        lower[i-1] = nu;
    }

    upper[0]=0;
    lower[mSize-1] =0;

    for (int i = 1; i < (mSize-1); i++)
    {
        diag[i] = identity - d_nu;
    }

    diag[0] = identity;
    diag[mSize-1] = identity;

    //setting variables to private class
    this -> mUpper = upper;
    this -> mLower = lower;
    this -> mDiag = diag;

}

TriMatrix::~TriMatrix() //destructor of TriMatrix
{
    delete[] mDiag;
    delete[] mLower;
    delete[] mUpper;
}

//overloaded function
double& TriMatrix::operator()(unsigned int i, unsigned int j)
{
    double ijthentry;

    //returning relevant value for specific matrix position
    if (i==j)
    {
        ijthentry = mDiag[i-1];
    }
    else if (i==(j-1))
    {
        ijthentry = mUpper[i-1];
    }
    else if (i==(j+1))
    {
        ijthentry = mLower[i-2];
    }
    else
    {
        ijthentry = 0;
    }
    return ijthentry;
};

//matrix multiplication function
void TriMatrix::matrixMultiplication (vector<double> &U, double ini_con_1, double ini_con_2)
{
    //temporary storage for the solution
    double* U_temp;
    U_temp = new double[U.size()];

    //assigning values considering boundary condition
    for (unsigned int i=0; i < U.size(); i++)
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
            U_temp[i] = this->operator()(i+1,i)*U[i-1] + this->operator()(i+1,i+1)*U[i] + this->operator()(i+1,i+2)*U[i+1];
        };
    };

    //inputting final value to the private variable
    this -> U_h = U_temp;

    //replacing current values with the processed values
    for (unsigned int i=0; i < U.size(); i++)
    {
        U[i] = U_h[i];
    };
};
