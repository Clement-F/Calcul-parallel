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
#include <cstdlib>
#include <algorithm>

// #include <mpi.h>

//  var global

int n;                      // taille matrices
int Nt,Nx,Ny;               // nb time step,  nb space step
double a,b;                 // dimension carre
double Time;                // Temps d arret
double dx,dy;               // space step
double dt;                  // time step
double U0, alpha;
int nb_savefile =10;

using namespace std;

// solution exacte
double u_ex(double x, double y)
{
    return sin(M_PI*x/a)*sin(M_PI*y/b);
}

// second membre
double f(double x, double y)
{
    return -M_PI*M_PI * (1/(a*a) + 1/(b*b)) * u_ex(x, y);
    // return 0;
}

// terme de bord
double V(double y)
{
    return 1-cos(2*M_PI *y/b);
}

// store une solution dans un fichier
void sol_to_file(vector<double> U, string name){
    std::ofstream myfile;
    myfile.open(name+".txt");
    for (int i=0; i<Nx+2; i++) {
        for (int j=0; j<Ny+2; j++) {
            myfile << U[i*(Ny+2) + j] << "\n";
        }
    }
    myfile.close();
}


// store un vecteur quelconque dans un fichier 
void save_to_file(vector<double> U, string name){
    std::ofstream myfile;
    myfile.open(name+".txt");
    for (int j=0; j<U.size(); j++) {
        myfile << U[j] << "\n";
    }
    myfile.close();
}

// donne l'erreur d approximation par Gauss Seidel 
// necessite d'input les parametre globaux avant d'appeler la fonction
double GS(){

    // Paramètres du problème

    dx = a/(Nx+1);  dy = b/(Ny+1);
    double dx2= dx*dx;

    // donne les parametre qui lui ont ete donne
    cout<<" parametres du probleme : \n";
    cout<<" parametre d'espace : Nx ="<<Nx<<", Ny ="<<Ny<<", a="<<a<<", b="<<b<<", dx="<<dx<<", dy="<<dy<<'\n';
    cout<<" parametre de temps : Nt ="<<Nt<<", Time ="<<Time<<'\n';
    cout<<" parametre de bords : U0 ="<<U0<<",  alpha="<<alpha<<'\n';



    if (a <= 0) {
        cout << "a doit être strictement positif!" << endl;
    }
    else if (b <= 0) {
        cout << "b doit être strictement positif!" << endl;
    }

    vector<double> U ((Nx+2)*(Ny+2));                   // vecteur de la solution au temps t
    vector<double> U_Next ((Nx+2)*(Ny+2));              // vecteur de la solution au temps t+dt
    
    // init U
    double x,y;
    for(int i=0; i<Nx+2;i++){
        x = i*dx;
        U[i*(Ny+2)]=        u_ex(x,0);      U_Next[i*(Ny+2)]=        u_ex(x,0);         // bas
        U[i*(Ny+2)+(Ny+1)]= u_ex(x,b);      U_Next[i*(Ny+2)+(Ny+1)]= u_ex(x,b);         // haut
        // cout<<u_ex(x,0)<<u_ex(x,b)<<'\n';
    }

    for(int j=0; j<Ny+2;j++){
        y = j*dy;
        U[j] =                u_ex(0,y);    U_Next[j] =                u_ex(0,y);       // gauche
        U[(Nx+1)*(Ny+2) +j] = u_ex(a,y);    U_Next[(Nx+1)*(Ny+2) +j] = u_ex(a,y);       // droite
        // cout<<u_ex(0,y)<<u_ex(a,y)<<'\n';
    }

    for(int i=1;i<Nx+1;i++){ for(int j=1;j<Ny+1;j++){
        U[i*(Ny+2)+j] = 0;
    }}

    // scheme
    double t=0;
    for(int l=1;l<=Nt;l++)
    {   
        t=t+dt;
        double progress = round(double(l)/Nt*10000)/100;
        if(progress - int(progress)<10e-7) cout<<"progress : "<<progress<<"%  t="<<t<<'\n';
        for(int i=1; i<Nx+1;i++)
        {
            x = i*dx;
            for(int j=1;j<Ny+1;j++)
            {
                y = j*dy;
                U_Next[i*(Ny+2)+j] = 0.25*(U[(i+1)*(Ny+2)+j] + U[i*(Ny+2)+(j+1)]);              // termes non update de GS
                U_Next[i*(Ny+2)+j] += 0.25*( U_Next[(i-1)*(Ny+2)+j]+ U_Next[i*(Ny+2)+(j-1)]);   // termes update de GS
                U_Next[i*(Ny+2)+j] += -0.25* dx2 *f(x,y);                                       // terme de bord
            } 
        }
        U.swap(U_Next);

    
    }

    
    vector<double> U_diff((Ny+2)*(Nx+2)); // U_diff = U_GS - U_ex

    // calcul l'erreur en norme inf de GS
    // definition de U_diff
    for(int i=0;i<Nx+2;i++)
    {
        x = i*dx;
        for(int j=0;j<Ny+2;j++)
        {
            y = j*dy;
            
            U_diff[i*(Ny+2)+j] = abs( U[i*(Ny+2)+j] - u_ex(x,y));

            if(U_diff[i*(Ny+2)+j] >0.5)
            {
                cout<<"( x:"<<x<<", y:"<<y<<") "<<i<<", "<<j<<"\n";
                cout<<U_diff[i*(Ny+2)+j]<<'\n'<<U[i*(Ny+2)+j]<<'\n'<<u_ex(x,y)<<'\n';
            }

        }
    }

    cout<<"norme infini \n";
    // calcul de la norme
    cout<<*max_element(U_diff.begin() , U_diff.end())<<" "<<endl;

    // save to file
    // sol_to_file(U_diff,"U_sol");   

    return (*max_element(U_diff.begin() , U_diff.end()));
}

int main (int argc, char* argv[])
{

    // on assigne les variables globales
    Nt =10000;
    a = 1;          b  = 1;      
    Time = 1;       dt = Time/Nt; 
    U0 = 1;         alpha = 0.5;

    vector<double> norme_inf;

    for(int k=2;k<=6;k++)
    {
        cout<<"============= "<<k<<" ===================== \n";
        Nx = int(pow(2,k));  Ny = Nx;
        dx = a/(Nx+1);  dy = b/(Ny+1);
        Nt = 2*Ny*Nx +100;      dt = Time/Nt; 
        norme_inf.push_back( GS());

    }

    save_to_file(norme_inf,"norme_inf");   


}