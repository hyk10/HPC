#include <iostream>
#include <iomanip>
#include <cstring>

using namespace std;

#include "TriMatrix.h"

TriMatrix::TriMatrix(double nu, double theta, int mSize) //constructor
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

    //setting variables to private class
    this -> mUpper_lhs = upper;
    this -> mLower_lhs = lower;
    this -> mDiag_lhs = diag;

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

    for (int i = 1; i < (mSize-1); i++)
    {
        diag[i] = identity - (1 - theta) * d_nu;
    }

    diag[0] = identity;
    diag[mSize-1] = identity;

    //setting variables to private classvU_heat
    this -> mUpper = upper;
    this -> mLower = lower;
    this -> mDiag = diag;

}

TriMatrix::~TriMatrix() //destructor of TriMatrix
{
    delete[] mDiag_lhs;
    delete[] mLower_lhs;
    delete[] mUpper_lhs;
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

    //calculating RHS of the equation (I-(1-theta)*nu*L)U^k
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


    double* upper;
    double* lower;
    double* diag;

    upper = mUpper_lhs;
    lower = mLower_lhs;
    diag = mDiag_lhs;

    // performing forward elimination
    for (unsigned int i = 1; i < U.size(); i++)
    {
        diag[i] = diag[i] - lower[i]*upper[i-1]/diag[i-1];
        U_temp[i] = U[i] - lower[i]*U[i-1]/diag[i-1];
    };


    //performing back substitution
    //the bottom value
    U_temp[U.size()-1] = U_temp[U.size()-1]/diag[U.size()-1];
    //the rest
    for (int i=U.size()-2; i >= 0; i--)
    {
        U_temp[i]=(U_temp[i]-upper[i]*U_temp[i+1])/diag[i];
    };

    this -> U_h = U_temp;

    //replacing current values with the processed values
    for (unsigned int i=0; i < U.size(); i++)
    {
        U[i] = U_h[i];
    };

};
