// Compilation
//    mpicxx jacobi.cpp
// Exécution (remplacer Nx par le nombre de points intérieurs du maillage)
//    mpirun -np k ./a.out Nx Ny a b

#include <iostream>
#include <fstream>
#include <mpi.h>
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
int Nmax = 10000;   // Nombre d'itérations max de l'algorithme
double tol = 1e-6;   // Tolérance de l'algorithme (critère d'arrêt)
bool EtudeErr = true; // Résolution du problème original (false) OU Etude de l'erreur sur solution particulère (true)


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


    if (a <= 0) {
        cout << "a doit être strictement positif!" << endl;
        return 1;
    }
    else if (b <= 0) {
        cout << "b doit être strictement positif!" << endl;
        return 1;
    }

    dx = a/(Nx+1);
    dy = b/(Ny+1);
    double cst = 1.0/(2.0/(dx*dx) + 2.0/(dy*dy));


    MPI_Init(&argc, &argv);

    int nbTask;
    int myRank;

    MPI_Comm_size(MPI_COMM_WORLD, &nbTask);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    vector<double> x(Nx+2, 0.0);
    vector<double> y(Ny+2, 0.0);
    for (int i=0; i<Nx+2; i++) {x[i] = i*dx;}
    for (int j=0; j<Ny+2; j++) {y[j] = j*dy;}



    // Subdivision pour le calcul parallèle
    int taille_bloc = Nx / nbTask;  // Nx points pour nbTask process
    int reste = Nx % nbTask;        // Reste des points (si Nx non divisible par nbTask)

    int Nx_local = taille_bloc + (myRank < reste ? 1 : 0);          // Tous les process ont taille_bloc points et les premiers (myRank < reste) récupèrent tous un point
    int decalage = myRank * taille_bloc + min(myRank, reste);       // Nombre de points déjà attribués avant ce process
    
    int i_start = 1 + decalage;
    int i_end = i_start + Nx_local - 1;

    // Initilisation & Conditions de bords
    vector<double> u( (Nx_local+2)*(Ny+2), U0);
    // Bas/Haut du domaine
    for (int i_local=1; i_local<=Nx_local; i_local++) {
        int i_global = decalage + i_local;
        if (EtudeErr) {
            u[i_local*(Ny+2)] = u_ex(x[i_global], y[0]);              // bas du domaine
            u[i_local*(Ny+2) + Ny+1] = u_ex(x[i_global], y[Ny+1]);    // haut du domaine
        } else {
            u[i_local*(Ny+2)] = U0;                                   // bas du domaine
            u[i_local*(Ny+2) + Ny+1] = U0;                            // haut du domaine
        }
    }
    // Bord gauche (seulement rang 0)
    if (myRank == 0) {
        for (int j=0; j<Ny+2; j++) {
            if (EtudeErr) {
                u[0*(Ny+2) + j] = u_ex(x[0], y[j]);
            } else {
                u[0*(Ny+2) + j] = U0*(1.0 + alpha*V(y[j], b));
            }
        }
    }
    // Bord droit (seulement rang nbTask-1)
    if (myRank == nbTask-1) {
        for (int j=0; j<Ny+2; j++) {
            if (EtudeErr) {
                u[(Nx_local+1)*(Ny+2) + j] = u_ex(x[Nx+1], y[j]);
            } else {
                u[(Nx_local+1)*(Ny+2) + j] = U0;
            }
        }
    }
    vector<double> uNew = u;


    // Schéma
    int iteration = 0;              // Nombre d'itérations
    double laplacien;               // Laplacien approché
    double maxResidu = tol + 1.0;   // Norme infinie du résidu
    double residu;                  // Résidu || Ax^l
    double maxResiduLoc = tol + 1.0;
    while (iteration < Nmax && maxResidu > tol) {
        maxResiduLoc = 0.0;
        if (myRank > 0) {
            MPI_Sendrecv(&u[(Ny+2)], Ny+2, MPI_DOUBLE, myRank-1, 0, &u[0], Ny+2, MPI_DOUBLE, myRank-1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if (myRank < nbTask-1) {
            MPI_Sendrecv(&u[Nx_local*(Ny+2)], Ny+2, MPI_DOUBLE, myRank+1, 1, &u[(Nx_local+1)*(Ny+2)], Ny+2, MPI_DOUBLE, myRank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        for (int i=1; i<=Nx_local; i++) {
            int i_global = decalage + i;
            for (int j=1; j<Ny+1; j++) {
                // Mise à jour
                uNew[i*(Ny+2) + j] = cst * ( 1/(dx*dx) * (u[(i+1)*(Ny+2) + j] + u[(i-1)*(Ny+2) + j]) + 1/(dy*dy) * (u[i*(Ny+2) + j+1] + u[i*(Ny+2) + j-1]) - f(x[i_global], y[j]));
                
                // Calcul de la norme du résidu
                laplacien = 1/(dx*dx) * (u[(i+1)*(Ny+2) + j] - 2*u[i*(Ny+2) + j] + u[(i-1)*(Ny+2) + j]) + 1/(dy*dy) * (u[(i*(Ny+2) + j+1)] - 2*u[i*(Ny+2) + j] + u[i*(Ny+2) + j-1]);
                residu = fabs( -laplacien + f(x[i_global], y[j]));
                maxResiduLoc = max(maxResiduLoc, residu);
            }
        }
        u.swap(uNew);

        MPI_Allreduce(&maxResiduLoc, &maxResidu, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        iteration++;
        // Suivi de la convergence
        if (myRank == 0) {
            if (iteration == 1 || iteration % 10000 == 0 || iteration == Nmax) {
                cout << "Itération " << iteration << ", résidu max: " << maxResidu << endl;
            }
        }
    }
    if (maxResidu <= tol) {
        if (myRank == 0) {
            cout << "Convergence après " << iteration << " itérations." << endl;
            cout << "Résidu max final: " << maxResidu << "." << endl;
        }
    }


    
    // Calcul erreur L^inf sur le maillage
    if (EtudeErr) {
        double err_inf = 0;
        double err_inf_loc = 0.0;
        for (int i_local=1; i_local<=Nx_local; i_local++) {
            int i_global = decalage + i_local;
            for (int j=0; j<Ny+2; j++) {
                err_inf_loc = max(err_inf_loc, fabs(u[i_local*(Ny+2) + j] - u_ex(x[i_global], y[j])));
            }
        }
        MPI_Reduce(&err_inf_loc, &err_inf, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        if (myRank == 0) {
            cout << "h = max(dx, dy) = " << max(dx, dy) << endl;
            cout << "Erreur L^inf sur le maillage: " << err_inf << endl;
        }
    }


    // Sauvegarder dans un .txt pour visualiser
    if (myRank == 0) {
        vector<double> u_global((Nx+2)*(Ny+2), 0.0);

        // Copie de la partie locale à myRank=0
        for (int i_local=1; i_local<=Nx_local; i_local++) {
            int i_global = decalage + i_local;
            for (int j=0; j<Ny+2; j++) {
                u_global[i_global*(Ny+2)+j] = u[i_local*(Ny+2)+j]; 
            }
        }

        // Réception des autres process (myRank != 0)
        for (int p=1; p<nbTask; p++) {
            int taille_bloc = Nx / nbTask;
            int reste = Nx % nbTask;
            int Nx_p = taille_bloc + (p < reste ? 1 : 0);
            int decalage_p = p * taille_bloc + min(p, reste);

            vector<double> temp((Nx_p+2)*(Ny+2)); // Vecteur temporaire pour recevoir
            MPI_Recv(temp.data(), (Nx_p+2)*(Ny+2), MPI_DOUBLE, p, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            // Ne prendre que les cellules intérieurs (éviter les cellules temporaires [0, Nx_p+1] et garder uniquement [1, Nx_p])
            for (int i_local=1; i_local<=Nx_p; i_local++) {
                int i_global = decalage_p + i_local;
                for (int j=0; j<Ny+2; j++) {
                    u_global[(i_global)*(Ny+2) + j] = temp[i_local*(Ny+2) + j];
                }
            }
        }

        // Écriture du fichier
        std::ofstream myfile;
        myfile.open("u_sol.txt");
        for (int i=0; i<Nx+2; i++) {
            for (int j=0; j<Ny+2; j++) {
                myfile << u_global[i*(Ny+2) + j] << "\n";
            }
        }
        myfile.close();
        cout << "Sauvegarde dans u_sol.txt effectuée." << endl;
    } else {
        MPI_Send(u.data(), (Nx_local+2)*(Ny+2), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }

    MPI_Finalize();

    return 0;
}