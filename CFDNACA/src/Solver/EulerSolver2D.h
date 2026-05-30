#pragma once
#include "../Core/Type.h"
#include <vector>

class EulerSolver2D {
public:
    int nx, ny;
    int spatialOrder; // 0 = Constant, 1 = 1st-Order Unlimited, 2 = TVD Minmod
    double maxCFL;
    StateVec U_inf;

    Field2D<StateVec> U;
    Field2D<double> Volumes;

    Field2D<FaceNormal> NormalsXi;
    Field2D<FaceNormal> NormalsEta;
    std::vector<FaceNormal> WallNormals;
    std::vector<FaceNormal> FarfieldNormals;

    EulerSolver2D(int num_cells_x, int num_cells_y, int order, double cfl);
    virtual ~EulerSolver2D() = default;

    void initialize(double rho_inf, double u_inf, double v_inf, double p_inf);
    virtual double computeTimeStep();

    // Evaluates spatial fluxes based on a given state and outputs to a residual array
    virtual void computeFluxResidual(const Field2D<StateVec>& state_in, Field2D<StateVec>& residualOut) const;

private:
    double minmod(double a, double b) const;
    StateVec minmod(const StateVec& a, const StateVec& b) const;
    StateVec reconstructXi(const Field2D<StateVec>& state_in, int i, int j, bool isRight) const;
    StateVec reconstructEta(const Field2D<StateVec>& state_in, int i, int j, bool isRight) const;
    StateVec extrapolateToWall(const Field2D<StateVec>& state_in, int i) const;
    StateVec extrapolateToFarfield(const Field2D<StateVec>& state_in, int i) const;
    int wrapXi(int i) const;
};
