#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h>
#include <femtool.hpp>

int main(){

    Mesh2D Omega;

    Read(Omega, "tp1-1.mesh");

    FeSpace Vh = FeSpace(Omega);

    // Function x = (x1,x2) -> cos(omega*x1)
    auto F = [](const R3& x){return std::cos(10*M_PI*(x[0]+x[1]));};

    // Evaluating nodal values of f at the degrees of freedom of Vh
    auto uex = Vh(F);

    Plot(Vh,uex,"output");

    return 0;
}