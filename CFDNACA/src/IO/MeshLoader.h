#pragma once
#include <string>
#include <vector>
#include "../Core/Type.h"

class MeshLoader {
public:
    static void loadMesh(const std::string& filename, int Nx, int Ny,
        Field2D<double>& Volumes,
        Field2D<FaceNormal>& NormalsXi,
        Field2D<FaceNormal>& NormalsEta,
        std::vector<FaceNormal>& WallNormals,
        std::vector<FaceNormal>& FarfieldNormals);
};