#pragma once
#include "EulerSolver2D.h"
#include "../Core/config.h"
#include "../Numerics/ViscousFlux.h"

class NavierStokesSolver2D : public EulerSolver2D {
public:
    // Raw nodal coordinates required for accurate Green-Gauss contour integrals
    Field2D<Point2D> Nodes;

    NavierStokesSolver2D(int num_cells_x, int num_cells_y, int order, double cfl, const Field2D<Point2D>& meshNodes);

    // Overridden methods to apply viscous N-S physics
    double computeTimeStep() override;
    void computeFluxResidual(const Field2D<StateVec>& state_in, Field2D<StateVec>& residualOut) const override;

private:

    // Helper to compute local dynamic viscosity based on scaling parameters
    double computeViscosity(double T) const;

    // Geometric and Interpolation Helpers
    Point2D getCellPos(int i, int j) const;
    StateVec interpolateToNode(const Field2D<StateVec>& U_in, int ni, int nj) const;
};