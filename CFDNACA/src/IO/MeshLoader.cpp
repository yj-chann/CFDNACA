#include "MeshLoader.h"
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <cmath>

void MeshLoader::loadMesh(const std::string& filename, int Nx, int Ny,
    Field2D<double>& Volumes,
    Field2D<FaceNormal>& NormalsXi,
    Field2D<FaceNormal>& NormalsEta,
    std::vector<FaceNormal>& WallNormals,
    std::vector<FaceNormal>& FarfieldNormals)
{
    Field2D<Point2D> nodes(Nx, Ny);
    std::ifstream file(filename);

    if (!file.is_open()) {
        throw std::runtime_error("Error: Could not open " + filename);
    }

    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            file >> nodes(i, j).x >> nodes(i, j).y;
        }
    }
    file.close();
    std::cout << "Successfully loaded " << Nx * Ny << " nodes from " << filename << "\n";

    int nx_cells = Nx - 1;
    int ny_cells = Ny - 1;

    for (int j = 0; j < ny_cells; ++j) {
        for (int i = 0; i < nx_cells; ++i) {
            Point2D sw = nodes(i, j);
            Point2D se = nodes(i + 1, j);
            Point2D nw = nodes(i, j + 1);
            Point2D ne = nodes(i + 1, j + 1);

            double diag1_x = ne.x - sw.x;
            double diag1_y = ne.y - sw.y;
            double diag2_x = nw.x - se.x;
            double diag2_y = nw.y - se.y;

            Volumes(i, j) = 0.5 * std::abs(diag1_x * diag2_y - diag1_y * diag2_x);
        }
    }

    for (int j = 0; j < ny_cells; ++j) {
        for (int i = 0; i < nx_cells; ++i) {
            Point2D bottom = nodes(i + 1, j);
            Point2D top = nodes(i + 1, j + 1);

            double dx = top.x - bottom.x;
            double dy = top.y - bottom.y;
            double length = std::sqrt(dx * dx + dy * dy);

            NormalsXi(i, j) = { dy / length, -dx / length, length };
        }
    }

    for (int j = 0; j < ny_cells; ++j) {
        for (int i = 0; i < nx_cells; ++i) {
            Point2D right = nodes(i + 1, j + 1);
            Point2D left = nodes(i, j + 1);

            double dx = left.x - right.x;
            double dy = left.y - right.y;
            double length = std::sqrt(dx * dx + dy * dy);

            NormalsEta(i, j) = { dy / length, -dx / length, length };
        }
    }

    WallNormals.resize(nx_cells);
    for (int i = 0; i < nx_cells; ++i) {
        Point2D right = nodes(i + 1, 0);
        Point2D left = nodes(i, 0);
        double dx = left.x - right.x;
        double dy = left.y - right.y;
        double length = std::sqrt(dx * dx + dy * dy);
        WallNormals[i] = { dy / length, -dx / length, length };
    }

    FarfieldNormals.resize(nx_cells);
    for (int i = 0; i < nx_cells; ++i) {
        Point2D right = nodes(i + 1, ny_cells);
        Point2D left = nodes(i, ny_cells);
        double dx = left.x - right.x;
        double dy = left.y - right.y;
        double length = std::sqrt(dx * dx + dy * dy);
        FarfieldNormals[i] = { dy / length, -dx / length, length };
    }
}