#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h>
#include <vector>
#include <femtool.hpp>
#include "partition.hpp"
#include <iostream>
#include <fstream>

using Mesh2DPart = std::vector<Mesh2D>;

int main() {
    Mesh2D Omega;
    Read(Omega,"maillages/carre-001.mesh");

    FeSpace Vh = FeSpace(Omega);

    std::pair<Mesh2DPart,CooMatrix<double>> Part4;
    auto [Sigma, Q] = Partition4(Omega);

    auto f = [](const R3& x){return std::cos(7*M_PI*(x[0]+x[1]));};

    std::vector<std::size_t> tbl(Sigma[0].size());
    for (std::size_t j=0; j<Sigma[0].size(); ++j) {
        auto& elemGamma = Sigma[0][j];
        for (std::size_t k=0; k<Omega.size(); ++k) {
            if (elemGamma == Omega[k]) {
                tbl[j] = k;
                break;
            }
        }
    }
    auto Uh = Restrict(Vh, Sigma[0], tbl).first;
    auto fh_r = Uh(f);

    Plot(Uh, fh_r, "restriction-output");

    return 0;
}