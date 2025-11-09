#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h>
#include <femtool.hpp>


int main(){

    Mesh3D Omega;
    Read(Omega, "tp1-2.mesh");

    FeSpace Vh = FeSpace(Omega);

    auto bord_pair = Boundary(Vh);
    FeSpace Wh = bord_pair.first;
    CooMatrix B = bord_pair.second;
    // FeSpace Wh(bord);

    auto F = [](const R3& x){return std::cos(5*M_PI*(x[0] + x[1] + x[2]));};
    auto uex = Wh(F);

    Plot(Wh, uex,"output");

    return 0;
}