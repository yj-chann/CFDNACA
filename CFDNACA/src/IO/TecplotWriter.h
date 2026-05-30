#pragma once
#include <string>
#include "../Core/Type.h"

class TecplotWriter {
public:
    static Field2D<Point2D> readNodes(const std::string& filename, int Nx, int Ny);
    static void exportSolution(const std::string& filename, int Nx, int Ny,
        const Field2D<Point2D>& nodes, const Field2D<StateVec>& U);
};
