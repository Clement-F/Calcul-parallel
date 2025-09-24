// Compilation:
//   mpicxx MPI_FD1D.cpp
// Execution (replace 'N' and 'L' with numbers of spatial/time steps):
//   mpirun -np 2 ./a.out 'N' 'L'

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
// #include <mpi.h>

//  var global

int n =100;                     // taille matrices
int Nt =100, Nx= 100;           // nb time step,  nb space step
double Length_x =1;             // dimension carre
double Time =1;                 // Temps d arret
double dx = Length_x/Nx;        // space step
double dt = Time/Nt;            // time step

double Length_y =Length_x;

using namespace std;

double f(double x, double y)
{
  return 0;
}


double U0(double x, double y)   // gaussian centrer en a/2, b/2
{
  double r2 = (x-0.5*Length_x)*(x-0.5*Length_x) + (y-0.5*Length_y)*(y-0.5*Length_y); // carré du rayon au centre du carre
  double a =100;    // parametre de centrage de la gaussienne
  return exp(-a*r2);
}

int main(int argc, char* argv[]){

    vector<double> U (Nx*Nx);     // vecteur de la solution
    
    // init U
    double x,y;
    for(int i; i<Nx;i++)
    {
        
        for(int j;j<Nx;j++)
        {
            U[i*Nx+j] = U0()
        } 
    }


    return 0;
}
