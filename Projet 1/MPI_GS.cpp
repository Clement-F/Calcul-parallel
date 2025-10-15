// Compilation:
//   mpicxx MPI_GS.cpp
// Execution 
//   mpirun -np 2 ./a.out 

#include <iostream>
#include <fstream>
#include <cmath>
#include <math.h>
#include <cstring>
#include <vector>
#include <cstdlib>
#include <algorithm>

#include <mpi.h>

//  var global

int n;                      // taille matrices
int Nt,Nx,Ny;               // nb time step,  nb space step
double a,b;                 // dimension carre
double Time;                // Temps d arret
double dx,dy;               // space step
double dt;                  // time step
double U0, alpha;
int nb_savefile =10;

int size_left_r,size_left_n,size_right_r,size_right_n;

using namespace std;

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
    if((i+j)%2 ==0){return true;}
    else{ return false;}
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


int main(int argc, char* argv[])
{

    MPI_Init(&argc, &argv);
    int nbTask;
    int myRank;
    MPI_Comm_size(MPI_COMM_WORLD, &nbTask);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    // Paramètres du problème


    a = 1;          b  = 1;      
    Time = 1;       dt = Time/Nt; 
    U0 = 1;         alpha = 0.5;
    Nx = 40;        Ny = Nx;
    dx = a/(Nx+1);  dy = b/(Ny+1);
    Nt = 10000;      dt = Time/Nt; 

    dx = a/(Nx+1);  dy = b/(Ny+1);
    double dx2= dx*dx;

    cout<<"des "<<nbTask<<" je suis le process "<<myRank<<'\n';
    if(myRank ==0)
    {
    cout<<"\n================================================== \n";
    cout<<"parametres du probleme : \n";
    cout<<"parametre d'espace : Nx ="<<Nx<<", Ny ="<<Ny<<", a="<<a<<", b="<<b<<", dx="<<dx<<", dy="<<dy<<'\n';
    cout<<"parametre de temps : Nt ="<<Nt<<", Time ="<<Time<<'\n';
    cout<<"parametre de bords : U0 ="<<U0<<",  alpha="<<alpha<<'\n';
    cout<<"================================================== \n \n";
    }


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
    
    int i_start = 1+ decalage;
    int i_end = i_start + Nx_local -1;

    cout<<"process "<<myRank<<" operates between "<<i_start<<" and "<<i_end<<" it has "<<Nx_local<<" nodes and "<<decalage<<" nodes have already been affected \n";
    
    bool isLeft_red =((i_start-1) %2==0);
    bool isRight_red=((i_end+1)  %2==0);

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
    if(myRank==(nbTask-1)){ for(int j=0; j<Ny+2;j++){y = j*dy;  U_loc[(Nx_local+1)*(Ny+2) +j] = u_ex(a,y);}}  // droite
    
    cout<<"valeurs de U_loc init attribuee pour process "<<myRank<<'\n';
    
    vector<double> U_loc_Next =U_loc;       // vecteur de la solution du probleme au temps t+
    
    int size_flank = ((Ny+2)/2)+1;
    // if(isLeft_red){  size_left_r = ceil(Ny/2);  size_left_n = floor(Ny/2);} else {size_left_r = floor(Ny/2);  size_left_n = ceil(Ny/2);}
    // if(isRight_red){ size_right_r = ceil(Ny/2); size_right_n = floor(Ny/2);}else {size_right_r = floor(Ny/2); size_right_n = ceil(Ny/2);}


    vector<double> U_left_R(size_flank), U_right_R(size_flank), U_left_N(size_flank), U_right_N(size_flank); 
    vector<double> U_left(Ny+2), U_right(Ny+2); 

    cout<<'\n';

    // scheme
    double t=0;
    for(int l=1;l<=Nt;l++)
    {   
        t=t+dt;
        if(myRank ==0){
        double progress = round(double(l)/Nt*100000)/1000;
        if(((progress - floor(progress))<10e-9) && int(progress)%10==0 ) cout<<"progress : "<<progress<<"%  t="<<t<<'\n';}
        
        // =============================================================================
        // =============================================================================

        // update de rouge
        // phase de com 
        // envoie de données superflue


        MPI_Request reqSendLeft, reqRecvLeft;
        MPI_Request reqSendRight, reqRecvRight;
        
        // =============================================================================
    
        if(myRank>0){
            for(int i=0;     i<Ny+2;     i++){
                if(is_red(i_start,i))
                { 
                    U_left[i]=U_loc[Ny+2+i];

                    if(i%2==0){ U_left_R[i/2]    = U_left[i]; if(l==1 && i<=1) cout<<myRank<<" S Left R "<<i<<" \n";}
                    if(i%2==1){ U_left_R[(i+1)/2]= U_left[i]; if(l==1 && i<=1) cout<<myRank<<" S Left R "<<i<<" \n";}
                } 
                // cout<<"R "<<i<<" "<<Ny+2+j<<'\n';
            }

            // if(l==1) cout<<"left "<<myRank<<" to "<<myRank-1<<'\n';
            MPI_Isend(&U_left_R[0], size_flank,MPI_DOUBLE,myRank-1,0,MPI_COMM_WORLD,&reqSendLeft);
            MPI_Irecv(&U_left_R[0], size_flank,MPI_DOUBLE,myRank-1,1,MPI_COMM_WORLD,&reqRecvLeft);

            MPI_Wait(&reqRecvLeft, MPI_STATUS_IGNORE);

            for(int i=0;     i<Ny+2;     i++){
                if(is_red(i_start-1,i))
                { 
                    if(i%2==0) U_left[i] = U_left_R[i/2];       if(l==1 && i<=1) cout<<"R Left R \n";
                    if(i%2==1) U_left[i] = U_left_R[(i+1)/2];   if(l==1 && i<=1) cout<<"R Left R \n";

                    U_loc[i]=U_left[i];
                }
            }
        } 

        // =============================================================================

        if(myRank<nbTask-1){
            for(int i=0;     i<Ny+2;     i++){
                if(is_red(i_end,i))
                {
                     U_right[i]=U_loc[(Nx_local)*(Ny+2)+i];

                    if(i%2==0){ U_right_R[i/2]    = U_right[i];  if(l==1 && i<=1) cout<<myRank<<" S Right R "<<i<<"\n";}
                    if(i%2==1){ U_right_R[(i+1)/2]= U_right[i];  if(l==1 && i<=1) cout<<myRank<<" S Right R "<<i<<"\n";}
                }
            }
            // if(l==1) cout<<"right "<<myRank<<" to "<<myRank+1<<'\n';
            MPI_Isend (&U_right_R[0], size_flank,MPI_DOUBLE,myRank+1,1,MPI_COMM_WORLD,&reqSendRight);
            MPI_Irecv (&U_right_R[0], size_flank,MPI_DOUBLE,myRank+1,0,MPI_COMM_WORLD,&reqRecvRight);

            MPI_Wait(&reqRecvRight, MPI_STATUS_IGNORE);
            for(int i=0;     i<Ny+2;     i++){
                if(is_red(i_end+1,i))
                {
                    if(i%2==0) U_right[i] = U_right_R[i/2];     if(l==1 && i<=1) cout<<"R Right R \n";
                    if(i%2==1) U_right[i] = U_right_R[(i+1)/2]; if(l==1 && i<=1) cout<<"R Right R \n";

                    U_loc[(Nx_local+1)*(Ny+2)+i]=U_right[i];
                } 
            }
        }
        
        // cout<<myRank<<" a envoyé et reçu ses données \n";

        // =============================================================================
        // =============================================================================

        for(int i=1; i<Nx_local+1;i++)
        {
            x = (i+decalage)*dx;
            for(int j=1;j<Ny+1;j++)
            {
                y = j*dy;
                // if(l==1)cout<<"Red ,"<<myRank<<" ,"<<is_red(i,j)<<"  ( x:"<<x<<", y:"<<y<<") "<<i+decalage<<", "<<j<<"\n";
                if(is_red(i,j)){ 
                U_loc_Next[i*(Ny+2)+j] = 0.25*(U_loc[(i+1)*(Ny+2)+j] + U_loc[i*(Ny+2)+(j+1)]);      
                U_loc_Next[i*(Ny+2)+j] += 0.25*(U_loc[(i-1)*(Ny+2)+j]+ U_loc[i*(Ny+2)+(j-1)]);      
                U_loc_Next[i*(Ny+2)+j] += -0.25*dx2 *f(x,y);                                       
                }
            } 
        }
    
        // =============================================================================
        // =============================================================================

        // cout<<myRank<<" a calculé pour les noeuds rouges \n";

        // update de noir
        // phase de com noir
        // envoie de données superflue        

        if(myRank>0){
            for(int i=0;     i<Ny+2;     i++){
                if(not is_red(i_start,i))
                {
                    U_left[i]=U_loc[Ny+2+i]; 

                    if(i%2==0){ U_left_N[i/2]    = U_left[i]; if(l==1 && i<=1) cout<<myRank<<" S Left N "<<i<<"\n";}
                    if(i%2==1){ U_left_N[(i+1)/2]= U_left[i]; if(l==1 && i<=1) cout<<myRank<<" S Left N "<<i<<"\n";}
                }
                // cout<<"N "<<i<<" "<<Ny+2+j<<'\n';
            }

            // if(l==1) cout<<"left "<<myRank<<" to "<<myRank-1<<'\n';
            MPI_Isend(&U_left_N[0], size_flank,MPI_DOUBLE,myRank-1,0,MPI_COMM_WORLD,&reqSendLeft);
            MPI_Irecv(&U_left_N[0], size_flank,MPI_DOUBLE,myRank-1,1,MPI_COMM_WORLD,&reqRecvLeft);

            MPI_Wait(&reqRecvLeft, MPI_STATUS_IGNORE);

            for(int i=0;     i<Ny+2;     i++){
                if(not is_red(i_start-1,i))
                {

                    if(i%2==0) U_left[i] = U_left_N[i/2];       if(l==1 && i<=1) cout<<"R Left N \n";
                    if(i%2==1) U_left[i] = U_left_N[(i+1)/2];   if(l==1 && i<=1) cout<<"R Left N \n";

                    U_loc[i]=U_left[i];

                }
            }
        } 
        
        // =============================================================================

        if(myRank<nbTask-1){
            for(int i=0;     i<Ny+2;     i++)
            {
                if(not is_red(i_end,i))
                {
                    U_right[i]=U_loc[(Nx_local)*(Ny+2)+i];
                     
                    if(i%2==0){ U_right_N[i/2]    = U_right[i];  if(l==1 && i<=1) cout<<myRank<<" S Right N "<<i<<"\n";}
                    if(i%2==1){ U_right_N[(i+1)/2]= U_right[i];  if(l==1 && i<=1) cout<<myRank<<" S Right N "<<i<<"\n";}
                }
            }

            // if(l==1) cout<<"right "<<myRank<<" to "<<myRank+1<<'\n';
            MPI_Isend (&U_right_N[0], size_flank,MPI_DOUBLE,myRank+1,1,MPI_COMM_WORLD,&reqSendRight);
            MPI_Irecv (&U_right_N[0], size_flank,MPI_DOUBLE,myRank+1,0,MPI_COMM_WORLD,&reqRecvRight);

            MPI_Wait(&reqRecvRight, MPI_STATUS_IGNORE);

            for(int i=0;     i<Ny+2;     i++)
            {
                if(not is_red(i_end+1,i))
                { 
                    if(i%2==0) U_right[i] = U_right_N[i/2];     if(l==1 && i<=1) cout<<"R Right N \n";
                    if(i%2==1) U_right[i] = U_right_N[(i+1)/2]; if(l==1 && i<=1) cout<<"R Right N \n";

                    U_loc[(Nx_local+1)*(Ny+2)+i]=U_right[i];

                }
            }
        }
    
        // =============================================================================
        // =============================================================================
    
        // cout<<myRank<<" a envoyé et reçu ses données updates \n";

        for(int i=1; i<Nx_local+1;i++)
        {
            x = (i+decalage)*dx;
            for(int j=1;j<Ny+1;j++)
            {
                y = j*dy;
                if(not is_red(i,j)){ 
                // cout<<"Black ,"<<myRank<<" ,"<<is_red(i,j)<<"  ( x:"<<x<<", y:"<<y<<") "<<i<<", "<<j<<"\n";
                U_loc_Next[i*(Ny+2)+j] = 0.25*(U_loc_Next[(i+1)*(Ny+2)+j] + U_loc_Next[i*(Ny+2)+(j+1)]);    
                U_loc_Next[i*(Ny+2)+j] += 0.25*(U_loc_Next[(i-1)*(Ny+2)+j]+ U_loc_Next[i*(Ny+2)+(j-1)]);    
                U_loc_Next[i*(Ny+2)+j] += -0.25* dx2 *f(x,y);                                               
                }
            } 
        }
        U_loc.swap(U_loc_Next);
    
        // =============================================================================
        // =============================================================================

        // cout<<myRank<<" a fini l iteration "<<l<<"\n";
    }

// ===========================================================================================
// ===========================================================================================

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
        myfile.open("U_sol.txt");
        for (int i=0; i<Nx+2; i++) {
            for (int j=0; j<Ny+2; j++) {
                myfile << U_global[i*(Ny+2) + j] << "\n";
            }
        }
        myfile.close();
        cout << "Sauvegarde dans u_sol.txt effectuée." << endl;

    
    vector<double> U_diff((Ny+2)*(Nx+2));

    for(int i=0;i<Nx+2;i++)
    {
        x = i*dx;
        for(int j=0;j<Ny+2;j++)
        {
            y = j*dy;
            
            U_diff[i*(Ny+2)+j] = abs( U_global[i*(Ny+2)+j] - u_ex(x,y));
            if(U_diff[i*(Ny+2)+j] >0.5)
            {
                // cout<<"( x:"<<x<<", y:"<<y<<") "<<i<<", "<<j<<"\n";
                // cout<<U_diff[i*(Ny+2)+j]<<'\n'<<U[i*(Ny+2)+j]<<'\n'<<u_ex(x,y)<<'\n';
            }

        }
    }

    cout<<"=========================== \nnorme infini \n"<<*max_element(U_diff.begin() , U_diff.end())<<" \n=========================== \n";

    // save to file
    // sol_to_file(U_diff,"U_sol");   
    } else {
        MPI_Send(U_loc.data(), (Nx_local+2)*(Ny+2), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }


    
    cout<<myRank<<" has ended its watch \n";
    MPI_Finalize();

    return 0;
}
