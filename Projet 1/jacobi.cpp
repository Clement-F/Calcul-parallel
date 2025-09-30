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
double a, b;         // Domaine d'étude (0,a)x(0,b)
double U0 = 1;       // Conditions aux bords 1
double alpha = 0.5;  // Conditions aux bords 2

int Nx, Ny;          // # pas d'espace horizontal, # pas d'espace vertical
double dx;           // Pas d'espace horizontal
double dy;           // Pas d'espace vertical
int Nmax = 100000;   // Nombre d'itérations max de l'algorithme
double tol = 1e-6;   // Tolérance de l'algorithme (critère d'arrêt)


// Solution particulière régulière (pour la validation du code)
double u_ex(double x, double y)
{
    return sin(M_PI*x/a)*sin(M_PI*y/b);
}


// Second membre de l'équation de Poisson
double f(double x, double y)
{
    return -M_PI*M_PI * (1/(a*a) + 1/(b*b)) * u_ex(x, y);
    //return 0;
}

// Conditions aux limites
double V(double y, const double b)
{
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

    double cst = 1/(2/(dx*dx) + 2/(dy*dy));

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
        u[i*(Ny+2)] = u_ex(x[i], y[0]);             // bas
        u[i*(Ny+2) + (Ny+1)] = u_ex(x[i], y[Ny+1]);
        //u[i*(Ny+2)] = U0;                         // bas du domaine
        //u[i*(Ny+2)+(Ny+1)] = U0;                  // haut du domaine
    }

    for (int j=0; j<Ny+2; j++) {
        u[0*(Ny+2) + j] = u_ex(x[0], y[j]);         // 0 ???
        u[(Nx+1)*(Ny+2) + j] = u_ex(x[Nx+1], y[j]);
        //u[(Nx+1)*(Ny+2) + j] = U0;                // droite du domaine
        //u[j] = U0 * (1.0 + alpha*V(y[j], b));     // gauche du domaine
    }

    vector<double> uNew = u;
    // Schéma
    int iteration = 0;              // Nombre d'itérations
    double laplacien;               // Laplacien approché
    double maxResidu = tol + 1.0;   // Norme infinie du résidu
    double residu;                  // Résidu || Ax^l
    while (iteration < Nmax && maxResidu > tol) {
        maxResidu = 0.0;
        for (int i=1; i<Nx+1; i++) {
            for (int j=1; j<Ny+1; j++) {
                // Mise à jour
                uNew[i*(Ny+2) + j] = cst * ( 1/(dx*dx) * (u[(i+1)*(Ny+2) + j] + u[(i-1)*(Ny+2) + j]) + 1/(dy*dy) * (u[i*(Ny+2) + j+1] + u[i*(Ny+2) + j-1]) - f(x[i], y[j]));
                
                // Calcul de la norme du résidu
                laplacien = 1/(dx*dx) * (u[(i+1)*(Ny+2) + j] - 2*u[i*(Ny+2) + j] + u[(i-1)*(Ny+2) + j]) + 1/(dy*dy) * (u[(i*(Ny+2) + j+1)] - 2*u[i*(Ny+2) + j] + u[i*(Ny+2) + j-1]);
                residu = fabs( -laplacien + f(x[i], y[j]));
                maxResidu = max(maxResidu, residu);
            }
        }
        u.swap(uNew);
        iteration++;
        // Suivi de la convergence
        if (iteration == 1 || iteration % 10000 == 0 || iteration == Nmax) {
            cout << "Itération " << iteration << ", résidu max: " << maxResidu << endl;
        }
    }

    if (maxResidu <= tol) {
        cout << "Convergence après " << iteration << " itérations." << endl;
        cout << "Résidu max final: " << maxResidu << "." << endl;
    }



    // Calcul erreur L^inf sur le maillage
    double err_inf = 0;
    for (int i=0; i<Nx+2; i++) {
        for (int j=0; j<Ny+2; j++) {
            err_inf = max(err_inf, fabs(u[i*(Ny+2) + j] - u_ex(x[i], y[j])));
        }
    }

    cout << "h = max(dx, dy) = " << max(dx, dy) << endl;
    cout << "Erreur L^inf sur le maillage: " << err_inf << endl;


    // Sauvegarder dans un .txt pour la visualisation
    std::ofstream myfile;
    myfile.open("u_sol.txt");
    for (int i=0; i<Nx+2; i++) {
        for (int j=0; j<Ny+2; j++) {
            myfile << u[i*(Ny+2) + j] << "\n";
        }
    }
    myfile.close();
    cout << "Sauvegarde dans u_sol.txt effectuée." << endl;
    return 0;
}
