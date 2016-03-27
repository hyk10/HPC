#ifndef TRIMATRIX_H_INCLUDED
#define TRIMATRIX_H_INCLUDED

#include <vector>

class TriMatrix
{
    public:
        TriMatrix(double nu, double theta, int mSize); //constructor

        ~TriMatrix(); //destructor

        //overloaded function
        double& operator()(unsigned int i, unsigned int j);

        //function which operates matrix multiplication
        void matrixMultiplication(std::vector<double> &U, double ini_con_1, double ini_con_2);



    private:
        double* mDiag_lhs;//mSize]; // diagonal value
        double* mUpper_lhs;//mSize-1]; // upper value
        double* mLower_lhs;//mSize-1]; // lower value
        double* mDiag;//mSize]; // diagonal value
        double* mUpper;//mSize-1]; // upper value
        double* mLower;//mSize-1]; // lower value
        double* U_h;

};


#endif // TRIMATRIX_H_INCLUDED
