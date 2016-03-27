#ifndef TRIMATRIX_H_INCLUDED
#define TRIMATRIX_H_INCLUDED

#include <vector>

class TriMatrix
{
    public:
        TriMatrix(double nu, int mSize); //constructor

        ~TriMatrix(); //destructor

        //overloaded function to return value of row i, column j
        double& operator()(unsigned int i, unsigned int j);

        //function which operates matrix multiplication
        void matrixMultiplication(std::vector<double> &U, double ini_con_1, double ini_con_2);



    private:
        double* mDiag;// diagonal value
        double* mUpper;// upper value
        double* mLower;// lower value
        double* U_h;// calculated U

};


#endif // TRIMATRIX_H_INCLUDED
