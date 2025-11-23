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
  Read(Omega,"mesh_fin.mesh");
  
  std::cout<<" \n read \n";
  auto Vh   = FeSpace(Omega);
  
  auto temp = Partition4(Omega,0);
  auto Gamma = temp.first[0]; auto R = temp.second;

  // std::cout<<R;

  auto tbl = Restrict(R,0);

  auto [Uh,P] = Restrict(Vh,Gamma,tbl);
  // std::cout<<P;

  
  auto F    = [](const R3& x){return std::cos(7*PI*(x[0]+x[1]));};
  
  // solution sur le maillage restreint /// fonctionne
  std::cout<<" \n assignment of ue sur Vh \n";
  auto ue   = Uh(F);
  std::cout<<" \n plot mesh\n";
  Plot(Uh,ue,"output-retreint-mesh");

  
  // solution restreinte      /// fonctionne pas 
  std::cout<<" \n assignment of ue sur Vh \n";
  auto fh   = Vh(F);  std::vector<double> f(dim(Vh),0);
  for(const auto& [j,k,m_jk]: P){if(int(m_jk) !=0)f[j]=fh[j];}
  std::cout<<" \n plot mesh\n";
  Plot(Vh,f,"output-function-mesh");


}