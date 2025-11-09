#include <iostream>
#include <fstream>
#include <sstream>
#include <memory>
#include <vector>
#include <filesystem>

using Mesh2DPart = std::vector<Mesh2D>;

std::pair<Mesh2DPart, CooMatrix<double>> 
 Partition4(const Mesh2D& Omega) {

    Mesh2DPart Meshes(4);
    int N = int(Omega.size());

    CooMatrix<double> Q;
}