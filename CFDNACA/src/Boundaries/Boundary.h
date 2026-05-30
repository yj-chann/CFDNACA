#pragma once
#include "../Core/Type.h"

class BoundaryConditions {
public:
    // Calculates boundary state and returns the physical flux directly
    static StateVec computeSolidWallFlux(const StateVec& U_star, const FaceNormal& inward_normal);
    static StateVec computeFarFieldFlux(const StateVec& U_star, const StateVec& U_inf, const FaceNormal& outward_normal);
}; 
