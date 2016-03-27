#include "TriMatrix.h"
#include <iostream>
#include <iomanip>
#include <cstring>
#include <vector>
#include <cmath>
#include <fstream>

using namespace std;



int main()
{
    double L;
    int N_x;
    double T;
    double N_t;
    double alpha;

    //inputting variables
    cout << "Enter L value in double format(eg 1.0) :" << endl;
    cin >> L;
    cout << "Enter number of discretised domain in integer format (eg 20) :" << endl;
    cin >> N_x;
    cout << "Enter target time T in double format (eg 5.0) :" << endl;
    cin >> T;
    cout << "Enter number of time step in double format (eg 5000.0) :" << endl;
    cin >> N_t;
    cout << "Enter alpha value in double format (eg 1.0):" << endl;
    cin >> alpha;

    //boundary conditions
    double gamma_0 = 0.0;
    double gamma_1 = 0.0;

    //setting initial time to 0
    double refT = 0.0;

    //making vector to store solutions
    vector<double> vU_heat(N_x+1);

    double d_x = L/N_x;
    double d_t = T/N_t;
    double nu = alpha*d_t/(d_x*d_x);
    double deviation;
    double rmse;
    const double pi = 3.1415926535897;

    //text files to store time and corresponding mid point value
    ofstream vOut1("refT.txt");
    ofstream vOut2("Nx_2.txt");

    //constructing matrix
    TriMatrix IpnuL = TriMatrix(nu,N_x+1);

    //setting initial U
    for (unsigned int i = 0; i < vU_heat.size(); i++)
    {
        vU_heat[i] = sin(pi*i*d_x/L);
    };

    //processing multiplication
    while(refT < T)
    {
        vOut1 << refT << "," ;
        vOut2 << vU_heat[floor(N_x/2)+1] << "," ;
        IpnuL.matrixMultiplication(vU_heat, gamma_0, gamma_1);
        refT+= d_t;


    };

    //calculating deviation*deviation
    for (unsigned int i = 0; i < vU_heat.size(); i++)
    {
        deviation += (vU_heat[i]-sin(pi*i*d_x/L)*exp(-alpha*pi*pi*refT/(L*L)))*(vU_heat[i]-sin(pi*i*d_x/L)*exp(-alpha*pi*pi*refT/(L*L)));
    };
    //calculating root mean square error
    rmse = sqrt(deviation/(N_x+1));

    vOut1.close();
    vOut2.close();

    //storing time step, descritised d_x and rmse
    ofstream vOut3("dtdx.txt");
    vOut3 << d_t << "," << d_x << "," << rmse;
    vOut3.close();
};
