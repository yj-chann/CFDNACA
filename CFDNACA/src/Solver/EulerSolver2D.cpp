#include "EulerSolver2D.h"
#include "../Numerics/RoeFlux.h"
#include "../Boundaries/Boundary.h"
#include <algorithm>

EulerSolver2D::EulerSolver2D(int num_cells_x, int num_cells_y, int order, double cfl)
    : nx(num_cells_x), ny(num_cells_y), spatialOrder(order), maxCFL(cfl),
    U(nx, ny), Volumes(nx, ny), NormalsXi(nx, ny), NormalsEta(nx, ny),
    WallNormals(nx), FarfieldNormals(nx) {
}

void EulerSolver2D::initialize(double rho_inf, double u_inf, double v_inf, double p_inf) {
    double rhoE_inf = p_inf / Config::GAMMA_MINUS_ONE + 0.5 * rho_inf * (u_inf * u_inf + v_inf * v_inf);
    U_inf = StateVec(rho_inf, rho_inf * u_inf, rho_inf * v_inf, rhoE_inf);

    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            U(i, j) = U_inf;
        }
    }
}

int EulerSolver2D::wrapXi(int i) const {
    if (i < 0) return i + nx;
    if (i >= nx) return i - nx;
    return i;
}

// Scalar minmod function
double EulerSolver2D::minmod(double a, double b) const {
    if (a * b <= 0.0) {
        return 0.0; // Opposite signs or zero, return 0
    }
    if (a > 0.0) {
        return std::min(a, b);
    }
    return std::max(a, b);
}

// Vectorized minmod function applied to the StateVec
StateVec EulerSolver2D::minmod(const StateVec& a, const StateVec& b) const {
    StateVec result;
    result.rho = minmod(a.rho, b.rho);
    result.rhou = minmod(a.rhou, b.rhou);
    result.rhov = minmod(a.rhov, b.rhov);
    result.rhoE = minmod(a.rhoE, b.rhoE);
    return result;
}

// cell center i 0 -> nx-1 j 0 -> ny-1 
StateVec EulerSolver2D::reconstructXi(const Field2D<StateVec>& state_in, int i, int j, bool isRight) const {
    if (spatialOrder == 0) {
        // 0th-Order: Piecewise Constant
        return isRight ? state_in(wrapXi(i + 1), j) : state_in(wrapXi(i), j);
    }
    else if (spatialOrder == 1) {
        // 1st-Order: Unlimited (Your original logic)
        if (isRight) {
            int ip1 = wrapXi(i + 1);
            int ip2 = wrapXi(i + 2);
            return state_in(ip1, j) * 1.5 - state_in(ip2, j) * 0.5;
        }
        else {
            int i_curr = wrapXi(i);
            int im1 = wrapXi(i - 1);
            return state_in(i_curr, j) * 1.5 - state_in(im1, j) * 0.5;
        }
    }
    else {
        // 2nd-Order: TVD with Minmod Limiter
        if (isRight) {
            int i_curr = wrapXi(i);
            int ip1 = wrapXi(i + 1);
            int ip2 = wrapXi(i + 2);

            StateVec U_ip1 = state_in(ip1, j);
            StateVec dU_fw = state_in(ip2, j) - U_ip1;
            StateVec dU_bw = U_ip1 - state_in(i_curr, j);

            return U_ip1 - minmod(dU_fw, dU_bw) * 0.5;
        }
        else {
            int im1 = wrapXi(i - 1);
            int i_curr = wrapXi(i);
            int ip1 = wrapXi(i + 1);

            StateVec U_i = state_in(i_curr, j);
            StateVec dU_fw = state_in(ip1, j) - U_i;
            StateVec dU_bw = U_i - state_in(im1, j);

            return U_i + minmod(dU_fw, dU_bw) * 0.5;
        }
    }
}

StateVec EulerSolver2D::reconstructEta(const Field2D<StateVec>& state_in, int i, int j, bool isRight) const {
    if (spatialOrder == 0) {
        // 0th-Order: Piecewise Constant
        return isRight ? state_in(i, j + 1) : state_in(i, j);
    }
    else if (spatialOrder == 1) {
        // 1st-Order: Unlimited (Your original logic)
        if (isRight) {
            StateVec U_jp1 = state_in(i, j + 1);
            StateVec U_jp2 = (j + 2 < ny) ? state_in(i, j + 2) : U_jp1;
            return U_jp1 - (U_jp2 - U_jp1) * 0.5;
        }
        else {
            StateVec U_j = state_in(i, j);
            StateVec U_jm1 = (j - 1 >= 0) ? state_in(i, j - 1) : U_j;
            return U_j + (U_j - U_jm1) * 0.5;
        }
    }
    else {
        // 2nd-Order: TVD with Minmod Limiter
        if (isRight) {
            StateVec U_j = state_in(i, j);
            StateVec U_jp1 = state_in(i, j + 1);
            StateVec U_jp2 = (j + 2 < ny) ? state_in(i, j + 2) : U_jp1;

            StateVec dU_fw = U_jp2 - U_jp1;
            StateVec dU_bw = U_jp1 - U_j;

            return U_jp1 - minmod(dU_fw, dU_bw) * 0.5;
        }
        else {
            StateVec U_jp1 = state_in(i, j + 1);
            StateVec U_j = state_in(i, j);
            StateVec U_jm1 = (j - 1 >= 0) ? state_in(i, j - 1) : U_j;

            StateVec dU_fw = U_jp1 - U_j;
            StateVec dU_bw = U_j - U_jm1;

            return U_j + minmod(dU_fw, dU_bw) * 0.5;
        }
    }
}

StateVec EulerSolver2D::extrapolateToWall(const Field2D<StateVec>& state_in, int i) const {
    if (!spatialOrder) return state_in(i, 0);
    else return state_in(i, 0) * 1.5 - state_in(i, 1) * 0.5;
}

StateVec EulerSolver2D::extrapolateToFarfield(const Field2D<StateVec>& state_in, int i) const {
    if (!spatialOrder) return state_in(i, ny - 1);
    else return state_in(i, ny - 1) * 1.5 - state_in(i, ny - 2) * 0.5;
}

double EulerSolver2D::computeTimeStep() {
    double dt_min = 1e9;
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            double u = U(i, j).u();
            double v = U(i, j).v();
            double a = U(i, j).a();
            double vol = Volumes(i, j);

            double ds_xi = NormalsXi(i, j).length;
            double ds_eta = NormalsEta(i, j).length;

            double xi_x = NormalsXi(i, j).nx / vol * ds_xi;
            double xi_y = NormalsXi(i, j).ny / vol * ds_xi;
            double grad_xi = std::sqrt(xi_x * xi_x + xi_y * xi_y);
            double lambda_xi = std::abs(xi_x * u + xi_y * v) + a * grad_xi;

            double eta_x = NormalsEta(i, j).nx / vol * ds_eta;
            double eta_y = NormalsEta(i, j).ny / vol * ds_eta;
            double grad_eta = std::sqrt(eta_x * eta_x + eta_y * eta_y);
            double lambda_eta = std::abs(eta_x * u + eta_y * v) + a * grad_eta;

            double dt_local = maxCFL / (lambda_xi + lambda_eta);
            dt_min = std::min(dt_min, dt_local);
        }
    }
    return dt_min;
}

void EulerSolver2D::computeFluxResidual(const Field2D<StateVec>& state_in, Field2D<StateVec>& residualOut) const {
    // 0. Zero out the residual array
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            residualOut(i, j) = StateVec();
        }
    }

    // Face LOOP!
    // 1. Interior Fluxes in Xi direction (i+1/2)
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            StateVec UL = reconstructXi(state_in, i, j, false);
            StateVec UR = reconstructXi(state_in, i, j, true);
            StateVec flux = RoeFlux::computeFlux(UL, UR, NormalsXi(i, j));

            residualOut(i, j) = residualOut(i, j) - flux;
            residualOut(wrapXi(i + 1), j) = residualOut(wrapXi(i + 1), j) + flux;
        }
    }

	// 2. Interior Fluxes in Eta direction (j+1/2)
    for (int j = 0; j < ny - 1; ++j) {
        for (int i = 0; i < nx; ++i) {
            StateVec UL = reconstructEta(state_in, i, j, false);
            StateVec UR = reconstructEta(state_in, i, j, true);
            StateVec flux = RoeFlux::computeFlux(UL, UR, NormalsEta(i, j));

            residualOut(i, j) = residualOut(i, j) - flux;
            residualOut(i, j + 1) = residualOut(i, j + 1) + flux;
        }
    }

    // 3. Apply Solid Boundary Fluxes (Direct flux addition)
    for (int i = 0; i < nx; ++i) {
        StateVec U_star = extrapolateToWall(state_in, i);
        StateVec wall_flux = BoundaryConditions::computeSolidWallFlux(U_star, WallNormals[i]);
        residualOut(i, 0) = residualOut(i, 0) + wall_flux;
    }

    // 4. Apply Far-field Boundary Fluxes (Direct flux subtraction)
    for (int i = 0; i < nx; ++i) {
        StateVec U_star = extrapolateToFarfield(state_in, i);
        StateVec far_flux = BoundaryConditions::computeFarFieldFlux(U_star, U_inf, FarfieldNormals[i]);
        residualOut(i, ny - 1) = residualOut(i, ny - 1) - far_flux;
    }
}