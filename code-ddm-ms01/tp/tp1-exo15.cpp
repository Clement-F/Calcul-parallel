#include <cmath>
#include <femtool.hpp>
#include<iostream>

const double PI = 3.141592653589793238463; 

int main(){

///// =================== TP1 Question 3 =======================


  // Instantiation of a 2D domain
    std::cout<<" \n instantiation \n";
  Mesh3D Omega;


  // Loading a 2D mesh
    std::cout<<" \n read \n";
  Read(Omega,"tp1-2.msh");
  

  // Assembly of a finite element space over Omega
    std::cout<<" \n assembly of Vh \n";
  auto Vh   = FeSpace(Omega);
  

    std::cout<<" \n assembly of boundary \n";
  auto temp   = Boundary(Vh); FeSpace Wh = temp.first; CooMatrix B = temp.second;
  
  
  // Function x = (x1,x2) -> cos(omega*x1)
    std::cout<<" \n lambda func \n";
  auto F    = [](const R3& x){return std::cos(5*PI*(x[0]+x[1]+x[2]));};


  // Manufactured "exact solution" obtained
  // by nodal evaluation of f at the degrees of freedom of Vh
    std::cout<<" \n assignment of ue \n";
  auto ue   = Vh(F);


  // Plotting exact solution with vizir4 
    std::cout<<" \n plot \n";
  Plot(Vh,ue,"tp1-2-output");



}