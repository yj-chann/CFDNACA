#pragma once
#include "../Solver/EulerSolver2D.h"
#include "../Core/Type.h"

class RungeKutta3 {
public:
    RungeKutta3(int nx, int ny);
    void step(EulerSolver2D& solver);

private:
    Field2D<StateVec> U_n;
    Field2D<StateVec> dU;
};
