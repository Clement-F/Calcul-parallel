// Compilation:
//   mpicxx MPI_FD1D.cpp
// Execution (replace 'N' and 'L' with numbers of spatial/time steps):
//   mpirun -np 2 ./a.out 'N' 'L'

#include <iostream>
#include <fstream>
#include <filesystem>
#include <cmath>
#include <math.h>
#include <cstring>
#include <vector>
// #include <mpi.h>

//  var global

int n;                     // taille matrices
int Nt,Nx,Ny;           // nb time step,  nb space step
double a,b;             // dimension carre
double Time;                 // Temps d arret
double dx,dy;        // space step
double dt;            // time step
double U0, alpha;

using namespace std;

double f(double x, double y)
{
  return 0;
}

void move(const string source, const string destination)
{
  char cmd[100];

  strcpy(cmd,"mv ");
  strcat(cmd,source.c_str());
  strcat(cmd," ");
  strcat(cmd,destination.c_str());
  system(cmd);
}

// double U0(double x, double y)   // gaussian centrer en a/2, b/2
// {
//   double r2 = (x-0.5*a)*(x-0.5*a) + (y-0.5*b)*(y-0.5*b); // carré du rayon au centre du carre
//   double a =100;    // parametre de centrage de la gaussienne
//   return exp(-a*r2);
// }

double V(double y)
{
    return 1-cos(2*M_PI *y/b);
}

int main(int argc, char* argv[]){

    // Paramètres du problème

    Nx = 100;   Ny = 100;   Nt = 100;
    a = 1;      b = 1;      
    dx = a/Nx;  dy= b/Ny;
    Time = 1;   dt =Time/Nt;
    U0 = 1;

    cout<<" parametres du probleme : \n";
    cout<<" parametre d'espace : Nx ="<<Nx<<", Ny ="<<Ny<<", a="<<a<<", b="<<b<<'\n';
    cout<<" parametre de temps : Nt ="<<Nt<<", Time ="<<Time<<'\n';
    cout<<" parametre de bords : U0 ="<<U0<<", alpha="<<alpha<<'\n';


    if (a <= 0) {
        cout << "a doit être strictement positif!" << endl;
        return 1;
    }
    else if (b <= 0) {
        cout << "b doit être strictement positif!" << endl;
    }

    vector<double> U ((Nx+2)*(Ny+2));           // vecteur de la solution au temps t
    vector<double> U_Next ((Nx+2)*(Ny+2));      // vecteur de la solution au temps t+dt
    
    // init U
    double x,y;
    for(int i=1; i<Nx;i++)
    {
        x = i*dx;
        U[i*Ny]     = U0;
        U[i*Ny+Ny]  = U0;

        for(int j=1;j<Ny;j++)
        {
            y = j*dy;
            U[i*Ny+j] = U0;
        } 
    }

    for(int j=0;j<Ny+1;j++)
    {
        y = j*dy;
        U[j] = U0 *(1+ alpha*V(y));
        U[(Nx+1)*(Ny) +j] = U0;
    }

    // scheme
    double t=0;
    for(int l=1;l<=Nt;l++)
    {
        t=t+dt;
        double progress = round(double(l)/Nt*10000)/100;
        cout<<"progress : "<<progress<<"%  t="<<t<<'\n';
        for(int i=1; i<Nx;i++)
        {
            x = i*dx;
            for(int j=1;j<Ny;j++)
            {
                y = j*dy;
                U_Next[i*Nx+j] = 0.25*(U[(i+1)*Nx+j] + U_Next[(i-1)*Nx+j]+ U[i*Nx+(j+1)]+ U[i*Nx+(j-1)]) - 0.25* dx*dx *f(x,y);   
            } 
        }
        U.swap(U_Next);
    }

    // save to file


    std::ofstream myfile;
    myfile.open("U_sol.txt");
    for (int i=0; i<= (Nx+2)*(Ny+2) ;i++ ){myfile<<U[i]<<"\n";}
    myfile.close();         

    return 0;
}
