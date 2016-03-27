#include "TriMatrix.h"
#include <iostream>
#include <iomanip>
#include <cstring>
#include <vector>

using namespace std;

double L=1.0;
double gamma_0=0.0;
double gamma_1=0.0;
double d_t=0.001;
double alpha=1.0;
int N_x=20;
double d_x=L/N_x;
double nu=alpha*d_t/(d_x*d_x);
vector<double> vU_heat(N_x+1);

int main()
{
    //constructing matrix
    TriMatrix IpnuL = TriMatrix(nu,N_x+1);

    //setting initial U
    for (unsigned int i = 0; i < vU_heat.size(); i++)
    {
        vU_heat[i]=i*d_x*(1-i*d_x);
    };

    //processing multiplication
    for(unsigned int i =0;; i++)
    {
        IpnuL.matrixMultiplication(vU_heat, gamma_0, gamma_1);

        //checking convergence
        if (vU_heat[N_x/2] < 0.000001)
        {
            return i;
            continue;
        };
    };
};
