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

int n;                     // taille matrices
int Nt,Nx;           // nb time step,  nb space step
double a,b;             // dimension carre
double Time;                 // Temps d arret
double dx;        // space step
double dt;            // time step


using namespace std;

double f(double x, double y)
{
  return 0;
}


double U0(double x, double y)   // gaussian centrer en a/2, b/2
{
  double r2 = (x-0.5*a)*(x-0.5*a) + (y-0.5*b)*(y-0.5*b); // carré du rayon au centre du carre
  double a =100;    // parametre de centrage de la gaussienne
  return exp(-a*r2);
}

int main(int argc, char* argv[]){

    // Paramètres du problème
    if (argc != 4) {
        cout << "Il faut entrer 3 arguments, ici il y en a " << argc - 1 << endl;
        return 1;
    }

    int Nx = atoi(argv[1]);
    double a = atoi(argv[2]);
    double b = atoi(argv[3]);
    double dx = (b-a)/(Nx+1);

    if (a <= 0) {
        cout << "a doit être strictement positif!" << endl;
        return 1;
    }
    else if (b <= 0) {
        cout << "b doit être strictement positif!" << endl;
    }

    vector<double> U (Nx*Nx);     // vecteur de la solution
    
    // init U
    double x,y;
    for(int i; i<Nx;i++)
    {
        
        for(int j;j<Nx;j++)
        {
            U[i*Nx+j];
        } 
    }


    return 0;
}
