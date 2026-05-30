#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include <vector>

#include "Core/config.h"
#include "Core/Type.h"
#include "IO/MeshLoader.h"
#include "IO/TecplotWriter.h"
#include "Solver/EulerSolver2D.h"
#include "TimeIntegration/RungeKutta3.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main() {
    // ---------------------------------------------------------
    // 1. Setup Grid Dimensions
    // ---------------------------------------------------------
    const int Nx = 151;
    const int Ny = 101;
    const int nx = Nx - 1;
    const int ny = Ny - 1;

    // Relative paths assume the executable is run from the project root folder.
    const std::string meshFile = "data/mesh.plt";
    const std::string initialOut = "output/initialize.plt";
    const std::string solutionOut = "output/solution.plt";

    // ---------------------------------------------------------
    // 2. Initialize Solver & Time Integrator
    // ---------------------------------------------------------
    int spatialOrder = 2;
    double maxCFL = 0.5;
    EulerSolver2D solver(nx, ny, spatialOrder, maxCFL);
    RungeKutta3 rk3(nx, ny);

    Field2D<Point2D> nodes(Nx, Ny);
    try {
        MeshLoader::loadMesh(meshFile, Nx, Ny, solver.Volumes, solver.NormalsXi,
            solver.NormalsEta, solver.WallNormals, solver.FarfieldNormals);
        nodes = TecplotWriter::readNodes(meshFile, Nx, Ny);
    }
    catch (const std::exception& e) {
        std::cerr << "Initialization failed: " << e.what() << "\n";
        return 1;
    }

    // ---------------------------------------------------------
    // 3. Define Free-stream Conditions (NACA0012 Transonic)
    // ---------------------------------------------------------
    //double Mach_inf = 0.8;
    //double alpha_deg = 1.25;
    //double alpha_rad = alpha_deg * M_PI / 180.0;

    //double rho_inf = 1.0;
    //double p_inf = 1.0;
    double a_inf = std::sqrt(Config::GAMMA * Config::p_inf / Config::rho_inf);
    double V_inf = Config::Mach_inf * a_inf;

    double u_inf = V_inf * std::cos(Config::alpha_rad);
    double v_inf = V_inf * std::sin(Config::alpha_rad);

    solver.initialize(Config::rho_inf, u_inf, v_inf, Config::p_inf);
    TecplotWriter::exportSolution(initialOut, Nx, Ny, nodes, solver.U);

    // ---------------------------------------------------------
    // 4. Time-Stepping & Monitoring Loop
    // ---------------------------------------------------------
    const int max_iterations = 100000;
    const int output_interval = 100;
    const double convergence_tolerance = 1e-8;

    std::cout << "\nStarting Euler Solver Iterations...\n";
    std::cout << std::setw(10) << "Iter"
        << std::setw(20) << "L2_Residual_Rho" << "\n";
    std::cout << std::string(30, '-') << "\n";

    for (int iter = 1; iter <= max_iterations; ++iter) {

        // Copy old density state to compute residual
        std::vector<double> rho_old(nx * ny);
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                rho_old[j * nx + i] = solver.U(i, j).rho;
            }
        }

        // Advance one time step using the external Time Integrator
        rk3.step(solver);

        // Compute L2 Norm of density residual
        double l2_res = 0.0;
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                double drho = solver.U(i, j).rho - rho_old[j * nx + i];
                l2_res += drho * drho;
            }
        }
        l2_res = std::sqrt(l2_res / (nx * ny));

        if (iter % output_interval == 0 || iter == 1) {
            std::cout << std::setw(10) << iter
                << std::setw(20) << std::scientific << std::setprecision(6) << l2_res << "\n";
        }

        if (l2_res < convergence_tolerance) {
            std::cout << "\nConverged at iteration " << iter << "!\n";
            break;
        }
    }

    // ---------------------------------------------------------
    // 5. Final Output Generation
    // ---------------------------------------------------------
    TecplotWriter::exportSolution(solutionOut, Nx, Ny, nodes, solver.U);

    return 0;
}