#include <cmath>
#include <math.h>
#include <vector>
#include <femtool.hpp>
#include <iostream>
#include <fstream>

const double PI = 3.141592653589793238463; 

int main(){


  // Instantiation of a 2D domain
  Mesh2D Omega;

  Read(Omega,"mesh00025.mesh");
  auto Vh   = FeSpace(Omega);
  
  auto F    = [](const R3& x){return std::cos(10*PI*(x[1]));};

  auto ue   = Vh(F);

  auto M    = Mass(Vh);
  auto K    = Stiffness(Vh);
  auto A    = K+M;
  auto b    = A(ue);

  int size =200;
  std::vector<long double> Conv(size);
  std::vector<double> err;

  for(int i=1; i<size;i++){
  // int k = int(pow(1.5,i));
  int k=i;

  auto uh   = cgsolve(A,b,k,1e-06);  

  // Error between discrete solution and
  // (manufactured) exact solution
  err  = ue-uh;

  // Displaying L2 norm of the error  
  Conv[i] = (K(err)|err)/(K(ue)|ue);
  std::cout << "pour k ="<<k<<", |ue-uh|^2/|ue|^2 = ";
  std::cout <<  Conv[i] << "\n";  
  }

  
  std::ofstream myfile;
  myfile.open("convergence.txt");
  for (int i=1; i<size; i++) {
    myfile << Conv[i] <<"\n";
  }
  
  myfile.close();
  // Plotting exact solution with vizir4 
  Plot(Vh,ue,"ue");

  // Plotting error with vizir4
  Plot(Vh,err,"err");
}