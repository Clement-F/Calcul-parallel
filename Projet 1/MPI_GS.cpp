// Compilation:
//   mpicxx MPI_GS.cpp
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
int Nx_n, start, end, diff;
double a,b;                 // dimension carre
double Time;                // Temps d arret
double dx,dy;               // space step
double dt;                  // time step
double U0, alpha;
int nb_savefile =10;

using namespace std;

// double U0(double x, double y)   // gaussian centrer en a/2, b/2
// {
//   double r2 = (x-0.5*a)*(x-0.5*a) + (y-0.5*b)*(y-0.5*b); // carré du rayon au centre du carre
//   double a =100;    // parametre de centrage de la gaussienne
//   return exp(-a*r2);
// }
double u_ex(double x, double y)
{
    return sin(M_PI*x/a)*sin(M_PI*y/b);
}


double f(double x, double y)
{
    return -M_PI*M_PI * (1/(a*a) + 1/(b*b)) * u_ex(x, y);
    // return 0;
}

bool is_red(int i,int j)
{
    if(i+j %2 ==0) return true;
}


double V(double y)
{
    return 1-cos(2*M_PI *y/b);
}


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



void save_to_file(vector<double> U, string name){
    std::ofstream myfile;
    myfile.open(name+".txt");
    for (int j=0; j<U.size(); j++) {
        myfile << U[j] << "\n";
    }
    myfile.close();
}

int main (int argc, char* argv[]){

    MPI_Init(&argc, &argv);

    int nbTask;
    int myRank;
    MPI_Comm_size(MPI_COMM_WORLD, &nbTask);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    // Paramètres du problème

    a = 1;          b  = 1;      
    Time = 1;       dt = Time/Nt; 
    U0 = 1;         alpha = 0.5;
    Nx = 100;       Ny = Nx;
    dx = a/(Nx+1);  dy = b/(Ny+1);
    Nt = 2*Ny*Nx +100;      dt = Time/Nt; 

    Nx_n = floor(double(Nx)/nbTask);
    start = myRank*(Nx_n+1); end = min((Nx_n+1)*(myRank+1)-1, N-1);
    diff = end-start +3;


    dx = a/(Nx+1);  dy = b/(Ny+1);
    double dx2= dx*dx;

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

    // solution du process
    vector<double> U_loc  ((Nx_n+2)*(Ny+2));            // vecteur de la solution sur le process au temps t
    vector<double> U_loc_Next  ((Nx_n+2)*(Ny+2));       // vecteur de la solution sur le process aut temps t+dt
    
    

    // init U
    double x,y;
    for(int i=0; i<Nx_n+2;i++){
        x = (i+myRank*Nx_n)*dx;
        U_loc[i*(Ny+2)]=        u_ex(x,0);      U_loc_Next[i*(Ny+2)]=        u_ex(x,0);         // bas
        U_loc[i*(Ny+2)+(Ny+1)]= u_ex(x,b);      U_loc_Next[i*(Ny+2)+(Ny+1)]= u_ex(x,b);         // haut
        // cout<<u_ex(x,0)<<u_ex(x,b)<<'\n';
    }

    for(int j=0; j<Ny+2;j++){
        y = j*dy;
        if(myRank==1)       {U[j] =                u_ex(0,y);    U_Next[j] =                u_ex(0,y);}       // gauche
        if(myRank==nbTask-1){U[(Nx+1)*(Ny+2) +j] = u_ex(a,y);    U_Next[(Nx+1)*(Ny+2) +j] = u_ex(a,y);}       // droite
        // cout<<u_ex(0,y)<<u_ex(a,y)<<'\n';
    }

    for(int i=start;i<end;i++){ for(int j=1;j<Ny+1;j++){
        U[i*(Ny+2)+j] = 0;
    }}





    // scheme
    double t=0;
    for(int l=1;l<=Nt;l++)
    {   
        t=t+dt;

        for(int i=1; i<Nx_n+1;i++)
        {
            x = i*dx;
            for(int j=1;j<Ny+1;j++)
            {
                y = j*dy;
                if(is_red(i,j))
                { 
                // cout<<"( x:"<<x<<", y:"<<y<<") "<<i<<", "<<j<<"\n";
                U_loc_Next[i*(Ny+2)+j] = 0.25*(U_loc[(i+1)*(Ny+2)+j] + U_loc[(i-1)*(Ny+2)+j]+ U_loc[i*(Ny+2)+(j+1)]+ U_loc[i*(Ny+2)+(j-1)]) - 0.25* dx2 *f(x,y);
                }
            } 
        }
        
        for(int i=1; i<Nx_n+1;i++)
        {
            x = i*dx;
            for(int j=1;j<Ny+1;j++)
            {
                y = j*dy;
                if(is_red(i,j))
                { 
                // cout<<"( x:"<<x<<", y:"<<y<<") "<<i<<", "<<j<<"\n";
                U_loc_Next[i*(Ny+2)+j] = 0.25*(U_loc[(i+1)*(Ny+2)+j] + U_loc[(i-1)*(Ny+2)+j]+ U_loc[i*(Ny+2)+(j+1)]+ U_loc[i*(Ny+2)+(j-1)]) - 0.25* dx2 *f(x,y);
                }
            } 
        }
        U.swap(U_Next);

    
    }

    
    vector<double> U_diff((Ny+2)*(Nx+2));

    for(int i=0;i<Nx+2;i++)
    {
        x = i*dx;
        for(int j=0;j<Ny+2;j++)
        {
            y = j*dy;
            
            U_diff[i*(Ny+2)+j] = abs( U[i*(Ny+2)+j] - u_ex(x,y));
            if(U_diff[i*(Ny+2)+j] >0.5)
            {
                // cout<<"( x:"<<x<<", y:"<<y<<") "<<i<<", "<<j<<"\n";
                // cout<<U_diff[i*(Ny+2)+j]<<'\n'<<U[i*(Ny+2)+j]<<'\n'<<u_ex(x,y)<<'\n';
            }

        }
    }

    cout<<"norme infini \n";
    cout<<*max_element(U_diff.begin() , U_diff.end())<<" "<<endl;

    // save to file
    // sol_to_file(U_diff,"U_sol");   

    return (*max_element(U_diff.begin() , U_diff.end()));
}
