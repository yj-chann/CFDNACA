#pragma once
#include "../Core/Type.h"

class BoundaryConditions {
public:
    // Calculates boundary state and returns the physical flux directly
    static StateVec computeSolidWallFlux(const StateVec& U_star, const FaceNormal& inward_normal);
    static StateVec computeFarFieldFlux(const StateVec& U_star, const StateVec& U_inf, const FaceNormal& outward_normal);

    // Viscous N-S Wall Boundary
    static StateVec computeViscousWallFlux(
        const StateVec& U_wall,
        const Gradient2D& grad_u,
        const Gradient2D& grad_v,
        const Gradient2D& grad_T,
        const FaceNormal& inward_normal,
        double mu_wall
    );
}; 
