#ifndef TRIMATRIX_H_INCLUDED
#define TRIMATRIX_H_INCLUDED

#include <vector>

class TriMatrix
{
    public:
        TriMatrix(double nu, double theta, int mSize, int my_rank); //constructor

        ~TriMatrix(); //destructor

        //function which operates matrix multiplication
        void matrixMultiplication(std::vector<double> &U, double ini_con, int my_rank);

    private:
        double* U_h0;
        double* RHS0;
        double* LHS0;
        double* U_h1;
        double* RHS1;
        double* LHS1;

};


#endif // TRIMATRIX_H_INCLUDED
