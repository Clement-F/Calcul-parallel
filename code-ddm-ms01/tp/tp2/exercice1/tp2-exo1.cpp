#include <cmath>
#include <math.h>
#include <vector>
#include <femtool.hpp>
#include <iostream>
#include <fstream>

const double PI = 3.141592653589793238463; 

int main(){

    std::cout<<" \n begin \n";

  // Instantiation of a 2D domain
  Mesh2D Omega;

    std::cout<<" \n instantiation \n";

  // Loading a 2D mesh
  Read(Omega,"mesh.mesh");
  
    std::cout<<" \n read \n";

  // Assembly of a finite element space over Omega
  auto Vh   = FeSpace(Omega);
  
    std::cout<<" \n assembly of Vh \n";
  
  // Function x = (x1,x2) -> cos(omega*x1)
  auto F    = [](const R3& x){return std::cos(10*PI*(x[0]+x[1]));};

    std::cout<<" \n lambda func \n";

  // Manufactured "exact solution" obtained
  // by nodal evaluation of f at the degrees of freedom of Vh
  auto uex   = Vh(F);

    std::cout<<" \n assignment of ue \n";

  // Assembly of the mass matrix for Vh
  auto M    = Mass(Vh);

  // Assembly of finite element matrix of the operator
  // -Delta + 1 with Neumann BC
  auto K    = Stiffness(Vh);
  auto A    = K+M;

  // Manufactured "right hand side"
  auto b    = A(ue);

//   int size =30;
//   std::vector<long double> Conv(size);
  long double Conv;

  
//   auto I = Identity(Vh); auto Q = I;

  auto uh   = PCGSolver(A,b,10);  

  // Error between discrete solution and
  // (manufactured) exact solution
  auto err  = ue-uh;

  // Displaying L2 norm of the error
  // Conv[i] = ((K(err)|err) + (M(err)|err)  )/((K(ue)|ue) + (M(ue)|ue));
  
  Conv = (K(err)|err)/(K(ue)|ue);
  
  std::cout << "pour k ="<<10<<", |ue-uh|^2/|ue|^2 = ";
  std::cout <<  Conv << "\n";  

  std::pair< Mesh2DPart,CooMatrix<double> > Partition;
  Partition = Partition4(Omega);
  Plot(Partition.first, "partition_mesh");

//   Plot(Vh,ue,"ue");

//   // Plotting error with vizir4
//   Plot(Vh,err,"err");
}