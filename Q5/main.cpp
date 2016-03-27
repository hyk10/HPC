#include "TriMatrix.h"
#include <iostream>
#include <iomanip>
#include <cstring>
#include <vector>
#include <cmath>

using namespace std;

#include "mpi.h"

int main(int argc,char *argv[])
{
    int p;
    int my_rank;

    MPI_Init(&argc, &argv); // Initialise MPI program
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); // Get rank
    MPI_Comm_size(MPI_COMM_WORLD, &p); // Get number of procs

    double L;
    int N_x;
    double T;
    double N_t;
    double alpha;
    double theta;
    double d_x;
    double d_t;
    double nu;

   //boundary condition
    double gamma_0 = 0.0;
    double gamma_1 = 0.0;

    //setting initial time to 0
    double refT = 0.0;

    //empty vector for
    vector<double> vU_heat(N_x+1);
    const double pi = 3.1415926535897;

    //variables to copy values through communicator
    double var0;
    double var1;
    int sizeOfV;

    //empty vectors to operate in different processors
    vector<double> vU_heat0;
    vector<double> vU_heat1;


    if (my_rank == 0)
    {
        //taking arguments through terminal
        cout << "Enter L value in double format(eg 1.0) :" << endl;
        cin >> L;
        cout << "Enter even number of discretised domain in integer format (eg 20) :" << endl;
        cin >> N_x;
        cout << "Enter target time T in double format (eg 5.0) :" << endl;
        cin >> T;
        cout << "Enter number of time step in double format (eg 5000.0) :" << endl;
        cin >> N_t;
        cout << "Enter alpha value in double format (eg 1.0):" << endl;
        cin >> alpha;
        cout << "Enter theta value(rad) in double format (eg 0.5):" << endl;
        cin >> theta;

        d_x=L/N_x;
        d_t=T/N_t;
        nu=alpha*d_t/(d_x*d_x);

        //setting size of matrix to divide
        sizeOfV=floor(N_x/2)+2;

        //sending variables to processor 1
        MPI_Send(&L, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
        MPI_Send(&N_x, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
        MPI_Send(&T, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
        MPI_Send(&N_t, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
        MPI_Send(&alpha, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
        MPI_Send(&theta, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
        MPI_Send(&d_x, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
        MPI_Send(&d_t, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
        MPI_Send(&nu, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
        MPI_Send(&sizeOfV, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);



        for (unsigned int i = 0; i < sizeOfV; i++)
        {
            vU_heat0[i] = sin(pi*i*d_x/L);
        };

        ////constructing matrix
        //TriMatrix IpnuL0 = TriMatrix(nu, theta, sizeOfV, 0);
    }
    else
    {
        //receiving variables from processor 0
        MPI_Recv(&L, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&N_x, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&T, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&N_t, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&alpha, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&theta, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&d_x, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&d_t, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&nu, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&sizeOfV, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        for (unsigned int i = sizeOfV-1; i < 2*sizeOfV; i++)
        {
            //initial conditions
            vU_heat1[i] = sin(pi*i*d_x/L);
        };

        ////constructing matrix
        //TriMatrix IpnuL1 = TriMatrix(nu, theta, sizeOfV, 1);
    }

    for (unsigned int i = 0; i < N_t; i++)
    {
        //constructing matrix
        TriMatrix IpnuL0 = TriMatrix(nu, theta, sizeOfV, 0);
        TriMatrix IpnuL1 = TriMatrix(nu, theta, sizeOfV, 1);
        if (my_rank==0)
        {
            //processing half of the matrix in processor 0
            IpnuL0.matrixMultiplication(vU_heat0, gamma_0, 0);
            var0 = vU_heat0[sizeOfV-3];

            //sending and receiving the bottom value of first half of U
            MPI_Send(&var0, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
            MPI_Recv(&var1, 1, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            vU_heat0[sizeOfV-1] = var1;
        }
        else
        {
            //processing the other half of matrix in processor 1
            IpnuL1.matrixMultiplication(vU_heat1, gamma_1, 1);
            var1 = vU_heat1[sizeOfV-3];
            //sending and receiving the bottom value of the other half of U
            MPI_Recv(&var0, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(&var1, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
            vU_heat1[0] = var0;
        }
    }

    //combining all values in processor 0
    if (my_rank == 0)
    {
        for (unsigned int i = 0; i < sizeOfV; i++)
        {
            MPI_Recv(&vU_heat1[i], 1, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        };


        for (unsigned int i = 0; i < sizeOfV; i++)
        {
            vU_heat0.push_back(vU_heat1[i]);
        };

        for (unsigned int i = 0; i < vU_heat0.size(); i++)
        {
            //showing final values through terminal
            cout << vU_heat0[i] << endl;
        };
    }
    else
    {
        for (unsigned int i = 0; i < sizeOfV; i++)
        {
            MPI_Send(&vU_heat1[i], 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
        };
    }
    MPI_Finalize(); // Finish MPI program
};
