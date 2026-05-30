#pragma once
#include "../Core/Type.h"

class ViscousFlux {
public:
    // Computes the viscous flux at a face interface
    static StateVec computeFlux(
        const StateVec& state_face,
        const Gradient2D& grad_u,
        const Gradient2D& grad_v,
        const Gradient2D& grad_T,
        const FaceNormal& normal,
        double mu_face
    );
}; 
