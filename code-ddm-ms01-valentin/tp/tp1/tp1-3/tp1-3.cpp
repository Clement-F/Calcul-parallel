#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h>
#include <femtool.hpp>
#include <string>
#include <fstream>
#include <sstream>

bool PlotSol = false;

bool carre_0025 = true;
bool carre_001 = true;
bool carre_0005 = true;
bool carre_00025 = true;

void cg_iter(const CooMatrix<double>& A,
             const std::vector<double>& b,
             std::vector<double>& x,
             std::vector<double>& r,
             std::vector<double>& p,
             double& r2)
{
    auto Ap = A*p;
    double pAp = (Ap|p);
    if (std::abs(pAp) < 1e-20) return;
    double alpha = r2 / pAp;
    x += alpha*p;
    r -= alpha*Ap;
    double r2_new = (r|r);
    double beta = r2_new / r2;
    p = r + beta*p;
    r2 = r2_new;
}

int main() {
    // Fichier unique pour le nombre d'itérations nécessaires pour atteindre || A*(x_k)-b || < 1e-6
    std::ofstream iter("iterations.txt");

    // Résolution sur chaque carré
    if (carre_0025) { 
        // Carré h=0.025
        // Lecture du maillage
        double h=0.025;
        Mesh2D Omega;
        Read(Omega, "carre-0025.mesh");


        // Elements Finis
        auto Vh = FeSpace(Omega);

        auto M = Mass(Vh);
        auto K = Stiffness(Vh);
        auto A = M + K;

        // Solution exacte sur le maillage
        auto uex = [](const R3& x){return std::cos(10*M_PI*(x[0]+x[1]));};
        auto u = Vh(uex);

        // Solution approchée par inverse directe
        auto b = A(u);
        auto InvA = Inv(A);
        auto uh = InvA(b);

        // Affichage de la solution
        if (PlotSol) {
            auto err1 = u - uh;
            Plot(Vh,err1,"err1");
        }


        // Etude de la convergence du Gradient Conjugué
        //int maxIter = 250;
        //int k_stop = maxIter;
        //auto ucg = u;
        //std::vector<long double> Conv(maxIter, 0.0);
        //for (int k=1; k<maxIter; k++) {     
        //    ucg = cgsolve(A, b, k);
        //    auto err = uh - ucg;
        //    Conv[k] = ((K(err)|err) + (M(err)|err))/((K(u)|u) + (M(u)|u));//Norm(u-uh)/Norm(u);
        //    if (Conv[k] < 1e-6) {
        //        std::cout << "h = " << h << "; k = " << k << std::endl;
        //        k_stop = k;
        //        break;
        //    }
        //} 
        int maxIter = 250;
        int k_stop = maxIter;
        auto ucg = std::vector<double>(b.size(), 0.0);
        auto r = b - A*ucg;
        auto p = r;
        double r2 = (r|r);

        std::vector<double> Conv(maxIter, 0.0);
        for (int k=1; k<maxIter; k++) {     
            cg_iter(A, b, ucg, r, p, r2);
            //ucg = cgsolve(A, b, k);
            auto err = uh - ucg;
            Conv[k] = ((K(err)|err) + (M(err)|err))/((K(u)|u) + (M(u)|u));//Norm(u-uh)/Norm(u);
            if (Conv[k] < 1e-6) {
                std::cout << "h = " << h << "; k = " << k << std::endl;
                k_stop = k;
                break;
            }
        }
        iter << h << "\t" << k_stop << "\n";

        // Sauvegarde
        std::string filename = "convergence_0025.txt";
        std::ofstream myfile;
        myfile.open(filename);
        for (int i=1; i<maxIter; i++){
            myfile << Conv[i] << "\n";
        }
        std::cout << "Sauvegarde dans " << filename << ": ok." << std::endl;
        myfile.close();
    }
    if (carre_001) { 
        // Carré h=0.01
        // Lecture du maillage
        double h = 0.01;
        Mesh2D Omega;
        Read(Omega, "carre-001.mesh");

        // Elements Finis
        auto Vh = FeSpace(Omega);
        auto M = Mass(Vh);
        auto K = Stiffness(Vh);
        auto A = M + K;

        // Solution exacte sur le maillage
        auto uex = [](const R3& x){return std::cos(10*M_PI*(x[0]+x[1]));};
        auto u = Vh(uex);

        // Solution approchée par inverse directe
        auto b = A(u);
        auto InvA = Inv(A);
        auto uh = InvA(b);

        // Etude de la convergence du Gradient Conjugué
        int maxIter = 250;
        int k_stop = maxIter;
        auto ucg = std::vector<double>(b.size(), 0.0);
        auto r = b - A*ucg;
        auto p = r;
        double r2 = (r|r);

        std::vector<double> Conv(maxIter, 0.0);
        for (int k=1; k<maxIter; k++) {     
            cg_iter(A, b, ucg, r, p, r2);
            //ucg = cgsolve(A, b, k);
            auto err = uh - ucg;
            Conv[k] = ((K(err)|err) + (M(err)|err))/((K(u)|u) + (M(u)|u));//Norm(u-uh)/Norm(u);
            if (Conv[k] < 1e-6) {
                std::cout << "h = " << h << "; k = " << k << std::endl;
                k_stop = k;
                break;
            }
        }

        iter << h << "\t" << k_stop << "\n";

        // Sauvegarde
        std::string filename = "convergence_001.txt";
        std::ofstream myfile;
        myfile.open(filename);
        for (int i=1; i<maxIter; i++){
            myfile << Conv[i] << "\n";
        }
        std::cout << "Sauvegarde dans " << filename << ": ok." << std::endl;
        myfile.close();
        }
    if (carre_0005) {
        // Carré h=0.005
        // Lecture du maillage
        double h = 0.005;
        Mesh2D Omega;
        Read(Omega, "carre-0005.mesh");

        // Elements Finis
        auto Vh = FeSpace(Omega);
        auto M = Mass(Vh);
        auto K = Stiffness(Vh);
        auto A = M + K;

        // Solution exacte sur le maillage
        auto uex = [](const R3& x){return std::cos(10*M_PI*(x[0]+x[1]));};
        auto u = Vh(uex);

        // Solution approchée par inverse directe
        auto b = A(u);
        auto InvA = Inv(A);
        auto uh = InvA(b);

        // Etude de la convergence du Gradient Conjugué
        //int maxIter = 250;
        //int k_stop = maxIter;
        //auto ucg = u;
        //std::vector<long double> Conv(maxIter, 0.0);
        //for (int k=1; k<maxIter; k++) {     
        //    ucg = cgsolve(A, b, k);
        //    auto err = uh - ucg;
        //    Conv[k] = ((K(err)|err) + (M(err)|err))/((K(u)|u) + (M(u)|u));//Norm(u-uh)/Norm(u);
        //    if (Conv[k] < 1e-6) {
        //        std::cout << "h = " << h << "; k = " << k << std::endl;
        //        k_stop = k;
        //        break;
        //    }
        //}
        int maxIter = 250;
        int k_stop = maxIter;
        auto ucg = std::vector<double>(b.size(), 0.0);
        auto r = b - A*ucg;
        auto p = r;
        double r2 = (r|r);

        std::vector<double> Conv(maxIter, 0.0);
        for (int k=1; k<maxIter; k++) {     
            cg_iter(A, b, ucg, r, p, r2);
            //ucg = cgsolve(A, b, k);
            auto err = uh - ucg;
            Conv[k] = ((K(err)|err) + (M(err)|err))/((K(u)|u) + (M(u)|u));//Norm(u-uh)/Norm(u);
            if (Conv[k] < 1e-6) {
                std::cout << "h = " << h << "; k = " << k << std::endl;
                k_stop = k;
                break;
            }
        }

        iter << h << "\t" << k_stop << "\n";
        
        // Sauvegarde
        std::string filename = "convergence_0005.txt";
        std::ofstream myfile;
        myfile.open(filename);
        for (int i=1; i<maxIter; i++){
            myfile << Conv[i] << "\n";
        }
        std::cout << "Sauvegarde dans " << filename << ": ok." << std::endl;
        myfile.close();
        }
    if (carre_00025) { 
        // Carré h=0.0025

        // Lecture du maillage
        double h = 0.0025;
        Mesh2D Omega;
        Read(Omega, "carre-00025.mesh");

        // Elements Finis
        auto Vh = FeSpace(Omega);
        auto M = Mass(Vh);
        auto K = Stiffness(Vh);
        auto A = M + K;

        // Solution exacte sur le maillage
        auto uex = [](const R3& x){return std::cos(10*M_PI*(x[0]+x[1]));};
        auto u = Vh(uex);

        // Solution approchée par inverse directe
        auto b = A(u);
        auto InvA = Inv(A);
        auto uh = InvA(b);

        // Etude de la convergence du Gradient Conjugué
        //int maxIter = 250;
        //int k_stop = maxIter;
        //auto ucg = u;
        //std::vector<long double> Conv(maxIter, 0.0);
        //for (int k=1; k<maxIter; k++) {     
        //    ucg = cgsolve(A, b, k);
        //    auto err = uh - ucg;
        //    Conv[k] = ((K(err)|err) + (M(err)|err))/((K(u)|u) + (M(u)|u));//Norm(u-uh)/Norm(u);
        //    if (Conv[k] < 1e-6) {
        //        std::cout << "h = " << h << "; k = " << k << std::endl;
        //        k_stop = k;
        //        break;
        //    }
        //}
        int maxIter = 250;
        int k_stop = maxIter;
        auto ucg = std::vector<double>(b.size(), 0.0);
        auto r = b - A*ucg;
        auto p = r;
        double r2 = (r|r);

        std::vector<double> Conv(maxIter, 0.0);
        for (int k=1; k<maxIter; k++) {     
            cg_iter(A, b, ucg, r, p, r2);
            //ucg = cgsolve(A, b, k);
            auto err = uh - ucg;
            Conv[k] = ((K(err)|err) + (M(err)|err))/((K(u)|u) + (M(u)|u));//Norm(u-uh)/Norm(u);
            if (Conv[k] < 1e-6) {
                std::cout << "h = " << h << "; k = " << k << std::endl;
                k_stop = k;
                break;
            }
        }


        iter << h << "\t" << k_stop << "\n";
        
        // Sauvegarde
        std::string filename = "convergence_00025.txt";
        std::ofstream myfile;
        myfile.open(filename);
        for (int i=1; i<maxIter; i++){
            myfile << Conv[i] << "\n";
        }
        std::cout << "Sauvegarde dans " << filename << ": ok." << std::endl;
        myfile.close();
        }
    iter.close();
    return 0;
}