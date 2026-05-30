#include "RungeKutta3.h"

RungeKutta3::RungeKutta3(int nx, int ny)
    : U_n(nx, ny), dU(nx, ny) {
}

void RungeKutta3::step(EulerSolver2D& solver) {
    double dt = solver.computeTimeStep();

    // Save initial state U^n
    for (int j = 0; j < solver.ny; ++j) {
        for (int i = 0; i < solver.nx; ++i) {
            U_n(i, j) = solver.U(i, j);
        }
    }

    // ==========================================
    // RK3 Stage 1: U^(1) = U^n + dt * Q(U^n)
    // ==========================================
    solver.computeFluxResidual(solver.U, dU);
    for (int j = 0; j < solver.ny; ++j) {
        for (int i = 0; i < solver.nx; ++i) {
            StateVec dtQ = dU(i, j) * (dt / solver.Volumes(i, j));
            solver.U(i, j) = U_n(i, j) + dtQ;
        }
    }

    // ==========================================
    // RK3 Stage 2: U^(2) = 3/4 U^n + 1/4 (U^(1) + dt * Q(U^(1)))
    // ==========================================
    solver.computeFluxResidual(solver.U, dU);
    for (int j = 0; j < solver.ny; ++j) {
        for (int i = 0; i < solver.nx; ++i) {
            StateVec dtQ = dU(i, j) * (dt / solver.Volumes(i, j));
            solver.U(i, j) = U_n(i, j) * 0.75 + (solver.U(i, j) + dtQ) * 0.25;
        }
    }

    // ==========================================
    // RK3 Stage 3: U^{n+1} = 1/3 U^n + 2/3 (U^(2) + dt * Q(U^(2)))
    // ==========================================
    solver.computeFluxResidual(solver.U, dU);
    for (int j = 0; j < solver.ny; ++j) {
        for (int i = 0; i < solver.nx; ++i) {
            StateVec dtQ = dU(i, j) * (dt / solver.Volumes(i, j));
            solver.U(i, j) = U_n(i, j) * (1.0 / 3.0) + (solver.U(i, j) + dtQ) * (2.0 / 3.0);
        }
    }
}