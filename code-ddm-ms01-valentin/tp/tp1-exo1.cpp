#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h>
#include <femtool.hpp>

bool question2 = true;
bool question3 = true;

int main(){
    if (question2) {
        Mesh2D Omega;
        Read(Omega, "maillages/tp1-1.mesh");

        FeSpace Vh = FeSpace(Omega);

        auto F = [](const R3& x){return std::cos(10*M_PI*(x[0]+x[1]));};
        auto uex = Vh(F);

        Plot(Vh,uex,"tp1-1-output");
        std::cout << "Question 2: ok" << std::endl;
    }
    if (question3) {
        Mesh3D Omega;
        Read(Omega, "maillages/tp1-2.mesh");

        FeSpace Vh = FeSpace(Omega);

        auto bord_pair = Boundary(Vh);
        FeSpace Wh = bord_pair.first;
        CooMatrix B = bord_pair.second;

        auto F = [](const R3& x){return std::cos(5*M_PI*(x[0] + x[1] + x[2]));};
        auto uex = Wh(F);

        Plot(Wh, uex,"tp1-2-output");
        std::cout << "Question 3: ok" << std::endl;
    }
    return 0;
}