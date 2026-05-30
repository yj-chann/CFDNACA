#pragma once
#include "../Core/Type.h"

class RoeFlux {
public:
    static StateVec computeFlux(const StateVec& UL, const StateVec& UR, const FaceNormal& normal);
private:
    static StateVec physicalFlux(const StateVec& U, double nx, double ny);
};
