#include <cmath>
#include <femtool.hpp>
#include<iostream>

const double PI = 3.141592653589793238463; 

int main(){

    std::cout<<" \n begin \n";

  // Instantiation of a 2D domain
  Mesh2D Omega;

    std::cout<<" \n instantiation \n";

  // Loading a 2D mesh
  Read(Omega,"tp1-1.mesh");
  
    std::cout<<" \n read \n";

  // Assembly of a finite element space over Omega
  auto Vh   = FeSpace(Omega);
  
    std::cout<<" \n assembly of Vh \n";
  
  // Function x = (x1,x2) -> cos(omega*x1)
  auto F    = [](const R3& x){return std::cos(10*PI*(x[0]+x[1]));};

    std::cout<<" \n lambda func \n";

  // Manufactured "exact solution" obtained
  // by nodal evaluation of f at the degrees of freedom of Vh
  auto ue   = Vh(F);

    std::cout<<" \n assignment of ue \n";

  // Plotting exact solution with vizir4 
  Plot(Vh,ue,"tp1-1-output");



}