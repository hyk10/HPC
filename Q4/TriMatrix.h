#ifndef TRIMATRIX_H_INCLUDED
#define TRIMATRIX_H_INCLUDED

#include <vector>

class TriMatrix
{
    public:
        TriMatrix(double nu, double theta, int mSize); //constructor

        ~TriMatrix(); //destructor

        //function which operates matrix multiplication
        void matrixMultiplication(std::vector<double> &U, double ini_con_1, double ini_con_2);

    private:
        double* U_h; //solution of the next time step
        double* RHS; // right hand side matrix
        double* LHS; // left hand side matrix

};


#endif // TRIMATRIX_H_INCLUDED
