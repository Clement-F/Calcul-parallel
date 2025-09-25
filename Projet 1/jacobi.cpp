// Compilation
//    g++ jacobi.cpp
// Exécution (remplacer Nx par le nombre de points intérieurs du maillage)
//    ./a.out Nx Ny a b

#include <iostream>
#include <fstream>
#include <cmath>
#include <math.h>
#include <vector>

using namespace std;

// Paramètres du domaine et de la discrétisation
double a,b;          // domaine d'étude [0,a]\times[0,b]
double U0 = 1;     // Conditions aux bords 1
double alpha = 0.5;    // Conditions aux bords 2

int Nt,Nx,Ny;        // # pas de temps, # pas d'espace horizontal, # pas d'espace vertical
double dx;           // pas d'espace horizontal
double dy;           // pas d'espace vertical
double dt;           // pas de temps
int Nmax = 100000;     // nombre d'itérations max de l'algorithme


// Second membre de l'équation de Poisson
double f(double x, double y)
{
    return 0;
}

double V(double y, const double b)
{
    // Conditions aux limites (coin en haut à gauche)
    return 1 - cos(2*M_PI*y/b);
}


int main(int argc, char* argv[])
{
    // Paramètres du problème
    if (argc != 5) {
        cout << "Il faut entrer 4 arguments: 'Nx' 'Ny' 'a' 'b', ici il y en a " << argc - 1 << endl;
        return 1;
    }

    Nx = atoi(argv[1]);
    Ny = atoi(argv[2]);
    a = atof(argv[3]);
    b = atof(argv[4]);

    dx = a/(Nx+1);
    dy = b/(Ny+1);

    if (a <= 0) {
        cout << "a doit être strictement positif!" << endl;
        return 1;
    }
    else if (b <= 0) {
        cout << "b doit être strictement positif!" << endl;
        return 1;
    }

    vector<double> x(Nx+2, 0.0);
    vector<double> y(Ny+2, 0.0);
    for (int i=0; i<Nx+2; i++) {x[i] = i*dx;}
    for (int j=0; j<Ny+2; j++) {y[j] = j*dy;}

    // Initilisation
    vector<double> u( (Nx+2)*(Ny+2), U0);
    for (int i=0; i<Nx+2; i++) {
        u[i*(Ny+2)] = U0;                       // bas du domaine
        u[i*(Ny+2)+(Ny+1)] = U0;                // haut du domaine
    }

    for (int j=0; j<Ny+2; j++) {
        u[(Nx+1)*(Ny+2) + j] = U0;              // droite du domaine
        u[j] = U0 * (1.0 + alpha*V(y[j], b));   // gauche du domaine
    }


    vector<double> uNew((Nx+2)*(Ny+2), 0.0);
    // Schéma
    for (int l=0; l<Nmax; l++) {
        for (int i=1; i<Nx+1; i++) {
            for (int j=1; j<Ny+1; j++) {
                uNew[i*(Ny+2) + j] = 0.25 * ( u[(i+1)*(Ny+2) + j] + u[(i-1)*(Ny+2) + j] + u[i*(Ny+2) + j+1] + u[i*(Ny+2) + j-1] ) - 0.25 * dx*dx * f(x[i], y[j]); 
            }
        }
        u.swap(uNew);
        if (l == 0 || l == 1 || l == 10 || l == Nmax-1) {
            cout << "Itération " << l << " échantillon points intérieurs/extérieurs:" << endl;
            for (int i = 1; i <= 3; i++) {
                for (int j = 1; j <= 3; j++) {
                    cout << "u[" << i << "," << j << "] = " << u[i*(Ny+2) + j] << "  ";
                }
            cout << endl;
            }
        }   
    }


    // Sauvegarder dans un .txt pour la visualisation
    std::ofstream myfile;
    myfile.open("u_sol.txt");
    for (int i=0; i<Nx+2; i++) {
        for (int j=0; j<Ny+2; j++) {
            myfile << u[i*(Ny+2) + j] << "\n";
        }
    }
    // for (int k=0; k<(Nx+2)*(Ny+2); k++) {
    //     myfile << u[k] << "\n";
    // }
    myfile.close();
    return 0;
}
