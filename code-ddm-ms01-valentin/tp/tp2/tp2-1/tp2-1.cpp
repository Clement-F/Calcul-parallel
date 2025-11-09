#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h>
#include <femtool.hpp>
#include <string>
#include <fstream>
#include <sstream>

bool carre_0025 = true;
bool carre_001 = true;
bool carre_0005 = true;
bool carre_00025 = true;

double err_relative = 1e-6;

int main() {
    std::ofstream iter("iterations_PCG.txt");
    std::cout << "Erreur relative recherchée: ||u_h(k) - u_h||/||u_h|| < 1e-6" << std::endl;
    if (carre_0025) {
        // Lire le maillage
        Mesh2D Omega;
        Read(Omega, "carre-0025.mesh");
        double h = 0.025;
        std::cout << "h = " << h << std::endl;
        
        // Espace Elements finis
        auto Vh = FeSpace(Omega);
        
        // Matrice du système
        auto M = Mass(Vh);
        auto K = Stiffness(Vh);
        auto A = M + K;
        
        auto uex = [](const R3& x){return std::cos(10*M_PI*x[0]);};
        auto u = Vh(uex);
        
        // Résolution par femtool
        auto b = A(u);
        auto InvA = Inv(A);
        auto uh = InvA(b);

        // Résolution par PCG
        std::cout << "Resolution par PCG" << std::endl;
        int maxIter = 250;
        int k_stop = maxIter;
        std::vector<double> Conv(maxIter);
        for (int k=0; k<maxIter; k++) {
            auto upcg = PCGSolver(A, b, k);
            auto err = u - upcg;
            Conv[k] = ((K(err)|err) + (M(err)|err))/((K(uh)|uh)+ (M(uh)|uh));
            if (Conv[k] < err_relative) {
                k_stop = k;
                std::cout << "Tolerance atteinte avec k = " << k_stop << " iterations" << std::endl;
                break;
            }
        }

        iter << h << "\t" << k_stop << "\n";
        
        std::string filename = "convergence_0025.txt";
        std::ofstream myfile;
        myfile.open(filename);
        for (int i=1; i<maxIter; i++){
            myfile << Conv[i] << "\n";
        }
        std::cout << "Sauvegarde dans " << filename << ": ok." << std::endl;
    }
    if (carre_001) {
        // Lire le maillage
        Mesh2D Omega;
        Read(Omega, "carre-001.mesh");
        double h = 0.01;
        std::cout << "h = " << h << std::endl;
        
        // Espace Elements finis
        auto Vh = FeSpace(Omega);
        
        // Matrice du système
        auto M = Mass(Vh);
        auto K = Stiffness(Vh);
        auto A = M + K;
        
        auto uex = [](const R3& x){return std::cos(10*M_PI*x[0]);};
        auto u = Vh(uex);
        
        // Résolution par femtool
        auto b = A(u);
        auto InvA = Inv(A);
        auto uh = InvA(b);
        
        // Résolution par PCG
        std::cout << "Resolution par PCG" << std::endl;
        int maxIter = 250;
        int k_stop = maxIter;
        std::vector<double> Conv(maxIter);
        for (int k=0; k<maxIter; k++) {
            auto upcg = PCGSolver(A, b, k);
            auto err = u - upcg;
            Conv[k] = ((K(err)|err) + (M(err)|err))/((K(uh)|uh)+ (M(uh)|uh));
            if (Conv[k] < err_relative) {
                k_stop = k;
                std::cout << "Tolerance atteinte avec k = " << k_stop << " iterations" << std::endl;
                break;
            }
        }
        
        iter << h << "\t" << k_stop << "\n";

        std::string filename = "convergence_001.txt";
        std::ofstream myfile;
        myfile.open(filename);
        for (int i=1; i<maxIter; i++){
            myfile << Conv[i] << "\n";
        }
        std::cout << "Sauvegarde dans " << filename << ": ok." << std::endl;
    }
    if (carre_0005) {
        // Lire le maillage
        Mesh2D Omega;
        Read(Omega, "carre-0005.mesh");
        double h = 0.005;
        std::cout << "h = " << h << std::endl;
        
        // Espace Elements finis
        auto Vh = FeSpace(Omega);
        
        // Matrice du système
        auto M = Mass(Vh);
        auto K = Stiffness(Vh);
        auto A = M + K;
        
        auto uex = [](const R3& x){return std::cos(10*M_PI*x[0]);};
        auto u = Vh(uex);
        
        // Résolution par femtool
        auto b = A(u);
        auto InvA = Inv(A);
        auto uh = InvA(b);
        
        // Résolution par PCG
        std::cout << "Resolution par PCG" << std::endl;
        int maxIter = 250;
        int k_stop = maxIter;
        std::vector<double> Conv(maxIter);
        for (int k=0; k<maxIter; k++) {
            auto upcg = PCGSolver(A, b, k);
            auto err = u - upcg;
            Conv[k] = ((K(err)|err) + (M(err)|err))/((K(uh)|uh)+ (M(uh)|uh));
            if (Conv[k] < err_relative) {
                k_stop = k;
                std::cout << "Tolerance atteinte avec k = " << k_stop << " iterations" << std::endl;
                break;
            }
        }

        iter << h << "\t" << k_stop << "\n";
        
        std::string filename = "convergence_0005.txt";
        std::ofstream myfile;
        myfile.open(filename);
        for (int i=1; i<maxIter; i++){
            myfile << Conv[i] << "\n";
        }
        std::cout << "Sauvegarde dans " << filename << ": ok." << std::endl;
    }
    if (carre_00025) {
        // Lire le maillage
        Mesh2D Omega;
        Read(Omega, "carre-00025.mesh");
        double h = 0.0025;
        std::cout << "h = " << h << std::endl;
        
        // Espace Elements finis
        auto Vh = FeSpace(Omega);
        
        // Matrice du système
        auto M = Mass(Vh);
        auto K = Stiffness(Vh);
        auto A = M + K;
        
        auto uex = [](const R3& x){return std::cos(10*M_PI*x[0]);};
        auto u = Vh(uex);
        
        // Résolution par femtool
        auto b = A(u);
        auto InvA = Inv(A);
        auto uh = InvA(b);
        
        // Résolution par PCG
        std::cout << "Resolution par PCG" << std::endl;
        int maxIter = 250;
        int k_stop = maxIter;
        std::vector<double> Conv(maxIter);
        for (int k=0; k<maxIter; k++) {
            auto upcg = PCGSolver(A, b, k);
            auto err = u - upcg;
            Conv[k] = ((K(err)|err) + (M(err)|err))/((K(uh)|uh)+ (M(uh)|uh));
            if (Conv[k] < err_relative) {
                k_stop = k;
                std::cout << "Tolerance atteinte avec k = " << k_stop << " iterations" << std::endl;
                break;
            }
        }

        iter << h << "\t" << k_stop << "\n";
        
        std::string filename = "convergence_00025.txt";
        std::ofstream myfile;
        myfile.open(filename);
        for (int i=1; i<maxIter; i++){
            myfile << Conv[i] << "\n";
        }
        std::cout << "Sauvegarde dans " << filename << ": ok." << std::endl;
    }
    iter.close();
    return 0;
}