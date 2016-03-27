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
        double* mDiag_lhs; // diagonal value of LHS matrix
        double* mUpper_lhs; // upper value of LHS matrix
        double* mLower_lhs; //lower value of LHS matrix
        double* mDiag; //diagonal value of RHS matrix
        double* mUpper; //upper value of RHS matrix
        double* mLower; //lower value of RHS matrix
        double* U_h; //solution of the next time step

};


#endif // TRIMATRIX_H_INCLUDED
