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

  Mesh2D Gamma;

  Mesh2D G = Vh.mesh();
  std::cout<<int(dim(Vh));
}