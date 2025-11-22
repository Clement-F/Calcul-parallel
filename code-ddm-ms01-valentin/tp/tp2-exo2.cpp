#include <cmath>
#include <math.h>
#include <vector>
#include <femtool.hpp>
#include "partition.hpp"
#include <iostream>
#include <fstream>

using Mesh2DPart = std::vector<Mesh2D>;

const double PI = 3.141592653589793238463;
bool question1_2 = true;
bool question4 = true;

int main(){

  Mesh2D Omega;
  Read(Omega,"maillages/carre-001.mesh");

  if (question1_2) {
    std::pair<Mesh2DPart,CooMatrix<double>> Part4;
    Part4 = Partition4(Omega);
    Plot(Part4.first, "partition4");

    std::pair<Mesh2DPart,CooMatrix<double>> Part16;
    Part16 = Partition16(Omega);
    Plot(Part16.first, "partition16");
  }
  if (question4) {
    std::pair<Mesh2DPart,CooMatrix<double>> Partition4_recouvrement;
    std::size_t nl4 = 2;
    Partition4_recouvrement = Partition4(Omega, nl4);
    Plot(Partition4_recouvrement.first, "partition4_recouvrement");


    std::pair<Mesh2DPart,CooMatrix<double>> Partition16_recouvrement;
    std::size_t nl16 = 2;
    Partition16_recouvrement = Partition16(Omega, nl16);
    Plot(Partition16_recouvrement.first, "partition16_recouvrement");
    
    // Vérification qu'il y ai bien recouvrement : on regarde les sommets partagés par les maillages

    int commun4 = 0;
    const double eps4 = 1e-12;
    auto& part4 = Partition4_recouvrement.first;
    for (int p = 0; p < 4; ++p) {
      for (const auto& v1 : part4[p].nodes()) {
        bool partage = false;
        for (int q = 0; q < 4; ++q) {
          if (q == p) continue;
          for (const auto& v2 : part4[q].nodes()) {
            if (fabs(v1[0]-v2[0])<eps4 && fabs(v1[1]-v2[1])<eps4 && fabs(v1[2]-v2[2])<eps4) {
              partage = true;
              break;
            }
          }
          if (partage) break;
        }
        if (partage) commun4++;
      }
    }
    std::cout << "Sommets partages (par coordonnees, partition 4) : " << commun4 << std::endl;
    

    int commun16 = 0;
    const double eps16 = 1e-12;
    auto& part16 = Partition16_recouvrement.first;
    for (int p = 0; p < 16; ++p) {
      for (const auto& v1 : part16[p].nodes()) {
        bool partage = false;
        for (int q = 0; q < 16; ++q) {
          if (q == p) continue;
          for (const auto& v2 : part16[q].nodes()) {
            if (fabs(v1[0]-v2[0])<eps16 && fabs(v1[1]-v2[1])<eps16 && fabs(v1[2]-v2[2])<eps16) {
              partage = true;
              break;
            }
          }
          if (partage) break;
        }
        if (partage) commun16++;
      }
    }

    std::cout << "Sommets partages (par coordonnees, partition 16) : " << commun16 << std::endl;
  }
}