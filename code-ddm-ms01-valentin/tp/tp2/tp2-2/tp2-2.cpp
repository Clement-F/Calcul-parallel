#include <cmath>
#include <math.h>
#include <vector>
#include <femtool.hpp>
#include "partition.hpp"
#include <iostream>
#include <fstream>

using Mesh2DPart = std::vector<Mesh2D>;

const double PI = 3.141592653589793238463;
bool question1_2 = false;
bool question4 = true;

int main(){

  Mesh2D Omega;
  Read(Omega,"carre-001.mesh");

  if (question1_2) {
    std::pair<Mesh2DPart,CooMatrix<double>> Partition;
    Partition = Partition16(Omega);
    Plot(Partition.first, "partition");
  }
  if (question4) {
    std::pair<Mesh2DPart,CooMatrix<double>> Partition;
    std::size_t nl = 2;
    Partition = Partition16(Omega, nl);
    Plot(Partition.first, "partition_recouvrement");
    
    // Vérification qu'il y ai bien recouvrement : on regarde les sommets partagés par les maillages
    int commun = 0;
    const double eps = 1e-12;

    auto& part = Partition.first;

    for (int p = 0; p < 16; ++p) {
      for (const auto& v1 : part[p].nodes()) {
        bool partage = false;
        for (int q = 0; q < 16; ++q) {
          if (q == p) continue;
          for (const auto& v2 : part[q].nodes()) {
            if (fabs(v1[0]-v2[0])<eps && fabs(v1[1]-v2[1])<eps && fabs(v1[2]-v2[2])<eps) {
              partage = true;
              break;
            }
          }
          if (partage) break;
        }
        if (partage) commun++;
      }
    }

    std::cout << "Sommets partages (par coordonnees) : " << commun << std::endl;
  }
}