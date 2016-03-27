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
        //first half, to be processed in processor 0
        double* U_h0; //solution of the next time step
        double* RHS0; // right hand side matrix
        double* LHS0; // left hand side matrix
        //second half, to be processed in processor 1
        double* U_h1; //solution of the next time step
        double* RHS1; // right hand side matrix
        double* LHS1; // left hand side matrix

};


#endif // TRIMATRIX_H_INCLUDED
