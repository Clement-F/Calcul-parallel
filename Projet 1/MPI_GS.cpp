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

    // Subdivision pour le calcul parallèle
    int taille_bloc = Nx / nbTask;  // Nx points pour nbTask process
    int reste = Nx % nbTask;        // Reste des points (si Nx non divisible par nbTask)

    int Nx_local = taille_bloc + (myRank < reste ? 1 : 0);          // Tous les process ont taille_bloc points et les premiers (myRank < reste) récupèrent tous un point
    int decalage = myRank * taille_bloc + min(myRank, reste);       // Nombre de points déjà attribués avant ce process
    
    int i_start = 1 + decalage;
    int i_end = i_start + Nx_local - 1;
    

    // solution du probleme
    vector<double> U_loc  ((Nx_local+2)*(Ny+2));            // vecteur de la solution du probleme au temps t
    
    // init U
    double x,y;
    for(int i=0; i<Nx_local+2;i++){
        int i_global = decalage + i;
        x = (i_global)*dx;
        U_loc[i*(Ny+2)]=        u_ex(x,0);         // bas
        U_loc[i*(Ny+2)+(Ny+1)]= u_ex(x,b);         // haut
        // cout<<u_ex(x,0)<<u_ex(x,b)<<'\n';
    }

    if(myRank==0){          for(int j=0; j<Ny+2;j++){y = j*dy;  U_loc[j] = u_ex(0,y);}}                     // gauche 
    if(myRank==nbTask-1){   for(int j=0; j<Ny+2;j++){y = j*dy;  U_loc[(Nx_local+1)*(Ny+2) +j] = u_ex(a,y);}}  // droite
    

    
    vector<double> U_loc_Next =U_loc;       // vecteur de la solution du probleme au temps t+dt
    vector<double> U_send_left(Ny+2);
    vector<double> U_send_right(Ny+2);
    // scheme
    double t=0;
    for(int l=1;l<=Nt;l++)
    {   
        t=t+dt;

        // update de rouge
        // phase de com 
        // envoie de données superflue
        if(myRank!=nbTask-1){
            for(int j=0; j<Ny+1; j++){ U_send_left[j] = U_loc[Nx_local*(Ny+2)+j];}
            MPI_Isend (U_send_left,Ny+2,MPI_DOUBLE,myRank+1,0,MPI_COMM_WORLD,&reqSendLeft);
            MPI_Irecv(&U_loc[Ny+2],Ny+2,MPI_DOUBLE,myRank-1,0,MPI_COMM_WORLD,&reqRecvLeft);
        } 

        if(myRank!=0){
            for(int j=0; j<Ny+1; j++){ U_send_right[j] = U_loc[j];}
            MPI_Isend (U_send_right,Ny+2,MPI_DOUBLE,myRank+1,0,MPI_COMM_WORLD,&reqSendLeft);
            MPI_Irecv (&U_loc[Nx_local*(Ny+2)],Ny+2,MPI_DOUBLE,myRank-1,0,MPI_COMM_WORLD,&reqRecvLeft);
        } 


        for(int i=1; i<Nx_local+1;i++)
        {
            x = i*dx;
            for(int j=1;j<Ny+1;j++)
            {
                y = j*dy;
                if(is_red(i,j))
                { 
                // cout<<"( x:"<<x<<", y:"<<y<<") "<<i<<", "<<j<<"\n";
                U_loc_Next[i*(Ny+2)+j] = 0.25*(U_loc[(i+1)*(Ny+2)+j] + U_loc[i*(Ny+2)+(j+1)]);      
                U_loc_Next[i*(Ny+2)+j] += 0.25*(U_loc[(i-1)*(Ny+2)+j]+ U_loc[i*(Ny+2)+(j-1)]);      
                U_loc_Next[i*(Ny+2)+j] += -0.25*dx2 *f(x,y);                                       
                }
            } 
        }

        // update de noir
        // phase de com noir
        // envoie de données superflue        

        if(myRank!=nbTask-1){
            for(int j=0; j<Ny+1; j++){ U_send_left[j] = U_loc[Nx_local*(Ny+2)+j];}
            MPI_Isend (U_send_left,Ny+2,MPI_DOUBLE,myRank+1,0,MPI_COMM_WORLD,&reqSendLeft);
            MPI_Irecv(&U_loc[Ny+2],Ny+2,MPI_DOUBLE,myRank-1,0,MPI_COMM_WORLD,&reqRecvLeft);
        } 
        
        if(myRank!=0){
            for(int j=0; j<Ny+1; j++){ U_send_right[j] = U_loc[j];}
            MPI_Isend (U_send_right,Ny+2,MPI_DOUBLE,myRank+1,0,MPI_COMM_WORLD,&reqSendRight);
            MPI_Irecv (&U_loc[Nx_local*(Ny+2)],Ny+2,MPI_DOUBLE,myRank-1,0,MPI_COMM_WORLD,&reqRecvRight);
        } 

        for(int i=1; i<Nx_local+1;i++)
        {
            x = i*dx;
            for(int j=1;j<Ny+1;j++)
            {
                y = j*dy;
                if(not is_red(i,j))
                { 
                U_loc_Next[i*(Ny+2)+j] = 0.25*(U_loc_Next[(i+1)*(Ny+2)+j] + U_loc_Next[i*(Ny+2)+(j+1)]);    
                U_loc_Next[i*(Ny+2)+j] += 0.25*(U_loc_Next[(i-1)*(Ny+2)+j]+ U_loc_Next[i*(Ny+2)+(j-1)]);    
                U_loc_Next[i*(Ny+2)+j] += -0.25* dx2 *f(x,y);                                               
                }
            } 
        }
        U_loc.swap(U_loc_Next);

    
    }

    // Sauvegarder dans un .txt pour visualiser
    if (myRank == 0) {
        vector<double> U_global((Nx+2)*(Ny+2), 0.0);

        // Copie de la partie locale à myRank=0
        for (int i_local=1; i_local<=Nx_local; i_local++) {
            int i_global = decalage + i_local;
            for (int j=0; j<Ny+2; j++) {
                U_global[i_global*(Ny+2)+j] = U_loc[i_local*(Ny+2)+j]; 
            }
        }

        for (int j=0; j<Ny+2; j++) {
            U_global[0*(Ny+2) + j] = U_loc[0*(Ny+2) + j];
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
                    U_global[(i_global)*(Ny+2) + j] = temp[i_local*(Ny+2) + j];
                }
            }
            // Prendre en compte la cellule Nx_p pour p = nbTask-1 pour avoir les conditions de bord droit du domaine
            if (p == nbTask - 1) {
                for (int j=0; j<Ny+2; j++) {
                    U_global[(Nx+1)*(Ny+2) + j] = temp[(Nx_p+1)*(Ny+2) + j];
                }
            }
        }
        if (nbTask == 1) {
            for (int j=0; j<Ny+2; j++) {
                U_global[(Nx+1)*(Ny+2) + j] = U_loc[(Nx_local+1)*(Ny+2)+j];
            }
        }

        // Écriture du fichier
        std::ofstream myfile;
        myfile.open("u_sol.txt");
        for (int i=0; i<Nx+2; i++) {
            for (int j=0; j<Ny+2; j++) {
                myfile << U_global[i*(Ny+2) + j] << "\n";
            }
        }
        myfile.close();
        cout << "Sauvegarde dans u_sol.txt effectuée." << endl;
    } else {
        MPI_Send(U_loc.data(), (Nx_local+2)*(Ny+2), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
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
    MPI_Finalize();

    return 0;
}
