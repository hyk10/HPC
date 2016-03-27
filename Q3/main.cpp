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
    double L=1.0;
    int N_x=20;
    double T=5.0;
    double N_t=5000.0;
    double alpha=1.0;
    double theta=0.5;

    //cout << "Enter L value in double format(eg 1.0) :" << endl;
    //cin >> L;
    //cout << "Enter number of discretised domain in integer format (eg 20) :" << endl;
    //cin >> N_x;
    //cout << "Enter target time T in double format (eg 5.0) :" << endl;
    //cin >> T;
    //cout << "Enter number of time step in double format (eg 5000.0) :" << endl;
    //cin >> N_t;
    //cout << "Enter alpha value in double format (eg 1.0):" << endl;
    //cin >> alpha;
    //cout << "Enter theta value(rad) in double format (eg 0.5):" << endl;
    //cin >> theta;

    double gamma_0 = 0.0;
    double gamma_1 = 0.0;
    double refT = 0.0;
    vector<double> vU_heat(N_x+1);
    double d_x=L/N_x;
    double d_t=T/N_t;
    double nu=alpha*d_t/(d_x*d_x);
    const double pi = 3.1415926535897;
    double deviation;
    double rmse;

    //constructing matrix
    TriMatrix IpnuL = TriMatrix(nu, theta, N_x+1);

    //setting initial U
    for (unsigned int i = 0; i < vU_heat.size(); i++)
    {
        vU_heat[i] = sin(pi*i*d_x/L);
    };

    //processing multiplication
    while(refT < T)
    {

        IpnuL.matrixMultiplication(vU_heat, gamma_0, gamma_1);
        refT+= d_t;
        cout << refT << endl;
    };



    for (unsigned int i = 0; i < vU_heat.size(); i++)
    {
        deviation += (vU_heat[i]-sin(pi*i*d_x/L)*exp(-alpha*pi*pi*refT/(L*L)))*(vU_heat[i]-sin(pi*i*d_x/L)*exp(-alpha*pi*pi*refT/(L*L)));

        //cout << sin(pi*i*d_x/L)*exp(-alpha*pi*pi*refT/(L*L)) << endl;
    };

    rmse = sqrt(deviation/(N_x+1));
    cout << rmse << endl;

    // for plotting
    ofstream vOut1("dtdx.txt");
    vOut1 << d_t << "," << d_x << "," << rmse << endl;
    vOut1.close();

};
