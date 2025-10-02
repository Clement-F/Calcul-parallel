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

// Paramètres du domaine
double a, b;          // Domaine d'étude (0,a)x(0,b)
double U0 = 1;        // Conditions aux bords 1
double alpha = 0.5;   // Conditions aux bords 2

// Paramètres de discrétisation
int Nx, Ny;           // # pas d'espace horizontal, # pas d'espace vertical
double dx;            // Pas d'espace horizontal
double dy;            // Pas d'espace vertical
int Nmax = 10000;     // Nombre d'itérations max de l'algorithme
double tol = 1e-6;    // Tolérance de l'algorithme (critère d'arrêt)
bool EtudeErr = false; // Résolution du problème original (false) OU Etude de l'erreur sur solution particulière (true)

// Solution particulière régulière (pour la validation du code)
double u_ex(double x, double y)
{
    return sin(M_PI*x/a)*sin(M_PI*y/b);
}


// Second membre de l'équation de Poisson
double f(double x, double y)
{
    if (EtudeErr) {
        return -M_PI*M_PI * (1/(a*a) + 1/(b*b)) * u_ex(x, y);
    } else {
        return 0;
    }
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

    // Initilisation
    vector<double> u( (Nx+2)*(Ny+2), 0.0);
    for (int i=0; i<Nx+2; i++) {
        if (EtudeErr) {
            u[i*(Ny+2)] = u_ex(x[i], y[0]);             // bas du domaine
            u[i*(Ny+2) + (Ny+1)] = u_ex(x[i], y[Ny+1]); // haut du domaine
        } else {
            u[i*(Ny+2)] = U0;                           // bas du domaine
            u[i*(Ny+2)+(Ny+1)] = U0;                    // haut du domaine
        }
    }

    for (int j=0; j<Ny+2; j++) {
        if (EtudeErr) {
            u[0*(Ny+2) + j] = u_ex(x[0], y[j]);         // droite du domaine
            u[(Nx+1)*(Ny+2) + j] = u_ex(x[Nx+1], y[j]); // gauche du domaine
        } else {
            u[(Nx+1)*(Ny+2) + j] = U0;                  // droite du domaine
            u[j] = U0 * (1.0 + alpha*V(y[j], b));       // gauche du domaine
        }
    }
    vector<double> uNew = u;    // Pour itérer
    vector<double> u0 = u;      // Pour écrire dans le .txt, ne pas modifier


    // Schéma
    int iteration = 0;              // Nombre d'itérations
    double laplacien;               // Laplacien approché
    double maxResidu = tol + 1.0;   // Norme infinie du résidu
    double residu;                  // Résidu || Ax^l
    double maxResiduLoc = tol + 1.0;

    int taille_bloc = Nx / nbTask;  // Subdivision en x
    int reste = Nx % nbTask;        // Si la division n'est pas ronde, on a un reste != 0 et on distribue

    int Nx_local = taille_bloc + (myRank < reste ? 1 : 0);          // Tous les process ont taille_bloc points et les premiers (myRank < reste) récupèrent tous un point
    int decalage = myRank * taille_bloc + min(myRank, reste);       // nombre de points déjà attribués avant ce process
    
    int i_start = 1 + decalage;
    int i_end = i_start + Nx_local - 1;

    while (iteration < Nmax && maxResidu > tol) {
        maxResiduLoc = 0.0;
        if (myRank > 0) {
            MPI_Sendrecv(&u[i_start*(Ny+2)], Ny+2, MPI_DOUBLE, myRank-1, 0, &u[(i_start-1)*(Ny+2)], Ny+2, MPI_DOUBLE, myRank-1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if (myRank < nbTask-1) {
            MPI_Sendrecv(&u[i_end*(Ny+2)], Ny+2, MPI_DOUBLE, myRank+1, 1, &u[(i_end+1)*(Ny+2)], Ny+2, MPI_DOUBLE, myRank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        for (int i=i_start; i<=i_end; i++) {
            for (int j=1; j<Ny+1; j++) {
                // Mise à jour
                uNew[i*(Ny+2) + j] = cst * ( 1/(dx*dx) * (u[(i+1)*(Ny+2) + j] + u[(i-1)*(Ny+2) + j]) 
                + 1/(dy*dy) * (u[i*(Ny+2) + j+1] + u[i*(Ny+2) + j-1]) - f(x[i], y[j]));
                
                // Calcul de la norme du résidu
                laplacien = 1/(dx*dx) * (u[(i+1)*(Ny+2) + j] - 2*u[i*(Ny+2) + j] + u[(i-1)*(Ny+2) + j]) + 1/(dy*dy) * (u[(i*(Ny+2) + j+1)] - 2*u[i*(Ny+2) + j] + u[i*(Ny+2) + j-1]);
                residu = fabs( -laplacien + f(x[i], y[j]));
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
    if (myRank == 0) {
        if (maxResidu <= tol) {
            cout << "Convergence après " << iteration << " itérations." << endl;
            cout << "Résidu max final: " << maxResidu << "." << endl;
        }
    }

    // Calcul erreur L^inf sur le maillage (si EtudeErr = true)
    if (EtudeErr) {
        double err_inf = 0;
        for (int i=0; i<Nx+2; i++) {
            for (int j=0; j<Ny+2; j++) {
                err_inf = max(err_inf, fabs(u[i*(Ny+2) + j] - u_ex(x[i], y[j])));
            }
        }
        if (myRank == 0) {
            cout << "h = max(dx, dy) = " << max(dx, dy) << endl;
            cout << "Erreur L^inf sur le maillage: " << err_inf << endl;
        }
    }

    // Sauvegarder dans un .txt pour visualiser
    if (myRank == 0) {
        // Copie de la partie locale à myRank=0
        for (int i=i_start; i<=i_end; i++) {
            for (int j=0; j<Ny+2; j++) {
                u0[i*(Ny+2)+j] = u[i*(Ny+2)+j]; 
            }
        }

        // Réception des autres process (myRank != 0)
        for (int p=1; p<nbTask; p++) {
            int taille_bloc = Nx / nbTask;
            int reste = Nx % nbTask;
            int Nx_p = taille_bloc + (p < reste ? 1 : 0);
            int decalage = p * taille_bloc + min(p, reste);
            int i_start_p = 1 + decalage;
            int i_end_p = i_start_p + Nx_p - 1;

            vector<double> temp(Nx_p*(Ny+2)); // Vecteur temporaire pour recevoir
            MPI_Recv(temp.data(), Nx_p*(Ny+2), MPI_DOUBLE, p, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (int i=0; i<Nx_p; i++) {
                for (int j=0; j<Ny+2; j++) {
                    u0[(i_start_p+i)*(Ny+2) + j] = temp[i*(Ny+2) + j];
                }
            }
        }


        // Écriture du fichier
        std::ofstream myfile;
        myfile.open("u_sol.txt");
        for (int i=0; i<Nx+2; i++) {
            for (int j=0; j<Ny+2; j++) {
                myfile << u0[i*(Ny+2) + j] << "\n";
            }
        }
        myfile.close();
        cout << "Sauvegarde dans u_sol.txt effectuée." << endl;
    } else {
        int Nx_local = i_end - i_start + 1;
        MPI_Send(&u[i_start*(Ny+2)], Nx_local*(Ny+2), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }

    MPI_Finalize();

    return 0;
}
