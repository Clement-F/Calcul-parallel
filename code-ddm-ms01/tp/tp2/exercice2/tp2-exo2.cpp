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

  std::pair< Mesh2DPart,CooMatrix<double> > Partition;
  Partition = Partition4(Omega,2);
  Plot(Partition.first, "partition_mesh");
}