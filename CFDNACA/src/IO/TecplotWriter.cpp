#include "TecplotWriter.h"
#include <iostream>
#include <fstream>
#include <stdexcept>

Field2D<Point2D> TecplotWriter::readNodes(const std::string& filename, int Nx, int Ny) {
    Field2D<Point2D> nodes(Nx, Ny);
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open " + filename + " for node extraction.");
    }
    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            file >> nodes(i, j).x >> nodes(i, j).y;
        }
    }
    return nodes;
}

void TecplotWriter::exportSolution(const std::string& filename, int Nx, int Ny,
    const Field2D<Point2D>& nodes, const Field2D<StateVec>& U) {
    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Warning: Could not open " << filename << " for writing.\n";
        return;
    }

    out << "VARIABLES = \"X\", \"Y\", \"Density\", \"U\", \"V\", \"Pressure\", \"Mach\"\n";
    out << "ZONE I=" << Nx << ", J=" << Ny << ", DATAPACKING=BLOCK, "
        << "VARLOCATION=([3,4,5,6,7]=CELLCENTERED)\n";

    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) out << nodes(i, j).x << "\n";
    }
    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) out << nodes(i, j).y << "\n";
    }

    int nx = Nx - 1;
    int ny = Ny - 1;

    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) out << U(i, j).rho << "\n";
    }
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) out << U(i, j).u() << "\n";
    }
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) out << U(i, j).v() << "\n";
    }
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) out << U(i, j).p() << "\n";
    }
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            double vel = std::sqrt(U(i, j).u() * U(i, j).u() + U(i, j).v() * U(i, j).v());
            out << vel / U(i, j).a() << "\n";
        }
    }
    out.close();
    std::cout << "Successfully exported solution to " << filename << "\n";
}