#include "NavierStokesSolver2D.h"
#include "../Boundaries/Boundary.h"
#include <algorithm>
#include <cmath>

NavierStokesSolver2D::NavierStokesSolver2D(int num_cells_x, int num_cells_y, int order, double cfl, const Field2D<Point2D>& meshNodes)
    : EulerSolver2D(num_cells_x, num_cells_y, order, cfl),
    Nodes(meshNodes) {
}


double NavierStokesSolver2D::computeViscosity(double T) const {
    // Define Sutherland Constants
    const double S = 110.4;          // Dimensional Sutherland constant for air [K]
    const double T_ref = 288.15;     // Dimensional free-stream temperature [K] (e.g., Sea Level)
	const double mu_ref = 1.7894e-5;  // Dimensional reference viscosity at T_ref [kg/(m*s)]

    // Apply Sutherland's Law
    // Formula: mu/mu_inf = (T/T_ref)^1.5 * (T_ref + S) / (T + S)

    double mu = mu_ref * std::pow(T/T_ref, 1.5) * (T_ref + S) / (T + S);

    return mu;
}




Point2D NavierStokesSolver2D::getCellPos(int i, int j) const {
    // Averages the 4 surrounding nodes to find the cell center coordinate
    int ni = i;
    int nip1 = (i + 1 < Nodes.nx) ? i + 1 : 0; // Node index wrapping 
    int nj = j;
    int njp1 = j + 1;

    double cx = 0.25 * (Nodes(ni, nj).x + Nodes(nip1, nj).x + Nodes(nip1, njp1).x + Nodes(ni, njp1).x);
    double cy = 0.25 * (Nodes(ni, nj).y + Nodes(nip1, nj).y + Nodes(nip1, njp1).y + Nodes(ni, njp1).y);
    return Point2D(cx, cy);
}

// ni , nj are node indices (0 to nx for ni, 0 to ny for nj) - note that nodes are (nx+1) x (ny+1)
StateVec NavierStokesSolver2D::interpolateToNode(const Field2D<StateVec>& U_in, int ni, int nj) const {
    // 1. Calculate horizontal weighting factor (alpha) - needed for both boundary and interior
    int ci_L = wrapXi(ni - 1);
    int ci_R = wrapXi(ni);

    Point2D p_L = Nodes(ni - 1, nj);
    Point2D p_R = Nodes(ni + 1, nj);
    Point2D p_node = Nodes(ni, nj);

    double d_L = std::hypot(p_L.x - p_node.x, p_L.y - p_node.y);
    double d_R = std::hypot(p_R.x - p_node.x, p_R.y - p_node.y);
    double alpha = d_R / (d_L + d_R + 1e-14); // Avoid div by zero

    // 2. Boundary mapped node interpolation
    if (nj == 0) { // Solid Wall No-slip Isothermal Boundary mapping
        // Horizontally interpolate pressure and density for the first two cell rows above the wall
        double p_row0 = alpha * U_in(ci_L, 0).p() + (1.0 - alpha) * U_in(ci_R, 0).p();
        double p_row1 = alpha * U_in(ci_L, 1).p() + (1.0 - alpha) * U_in(ci_R, 1).p();

        double rho_row0 = alpha * U_in(ci_L, 0).rho + (1.0 - alpha) * U_in(ci_R, 0).rho;
        double rho_row1 = alpha * U_in(ci_L, 1).rho + (1.0 - alpha) * U_in(ci_R, 1).rho;

        // Apply Boundary Conditions
        double p_node = 1.5 * p_row0 - 0.5 * p_row1; // Extrapolated pressure

        // Extrapolate density (If strictly isothermal, replace this with rho_node = p_node / (R * T_wall))
        double rho_node = 1.5 * rho_row0 - 0.5 * rho_row1;

        // Reconstruct and return conservative StateVec   
        double rhoE_node = p_node / Config::GAMMA_MINUS_ONE;
        return StateVec(rho_node, 0.0, 0.0, rhoE_node);
    }

    if (nj >= ny) return U_inf; // Simplistic Far-field mapping

    // 3. Interior Node Interpolation (nj > 0 && nj < ny)
    int cj_B = nj - 1;
    int cj_T = nj;

    // Vertical weighting factor (beta)
    Point2D p_B = Nodes(ni, nj - 1);
    Point2D p_T = Nodes(ni, nj + 1);
    double d_B = std::hypot(p_B.x - p_node.x, p_B.y - p_node.y);
    double d_T = std::hypot(p_T.x - p_node.x, p_T.y - p_node.y);
    double beta = d_T / (d_B + d_T + 1e-14);

    StateVec U_LB = U_in(ci_L, cj_B);
    StateVec U_RB = U_in(ci_R, cj_B);
    StateVec U_LT = U_in(ci_L, cj_T);
    StateVec U_RT = U_in(ci_R, cj_T);

    // Combined Bilinear Weights
    double w_LB = alpha * beta;
    double w_RB = (1.0 - alpha) * beta;
    double w_LT = alpha * (1.0 - beta);
    double w_RT = (1.0 - alpha) * (1.0 - beta);

    // Primitive variable distance-weighted interpolation
    double rho_node = w_LB * U_LB.rho + w_RB * U_RB.rho + w_LT * U_LT.rho + w_RT * U_RT.rho;
    double u_node = w_LB * U_LB.u() + w_RB * U_RB.u() + w_LT * U_LT.u() + w_RT * U_RT.u();
    double v_node = w_LB * U_LB.v() + w_RB * U_RB.v() + w_LT * U_LT.v() + w_RT * U_RT.v();
    double pressure_node = w_LB * U_LB.p() + w_RB * U_RB.p() + w_LT * U_LT.p() + w_RT * U_RT.p();

    // Reconstruct and return conservative StateVec 
    double rhoE_node = pressure_node / Config::GAMMA_MINUS_ONE + 0.5 * rho_node * (u_node * u_node + v_node * v_node);

    return StateVec(rho_node, rho_node * u_node, rho_node * v_node, rhoE_node);
}

double NavierStokesSolver2D::computeTimeStep() {
    double dt_min = 1e9;

    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            double rho = U(i, j).rho;
            double u = U(i, j).u();
            double v = U(i, j).v();
            double a = U(i, j).a();
            double mu = computeViscosity(U(i, j).T());
            double vol = Volumes(i, j);

            // Fetch nodes for geometric metrics differences
            int ni = i;
            int nip1 = (i + 1 < Nodes.nx) ? i + 1 : 0;

            double dy_xi = Nodes(nip1, j).y - Nodes(ni, j).y;
            double dy_eta = Nodes(ni, j + 1).y - Nodes(ni, j).y;
            double dx_xi = Nodes(nip1, j).x - Nodes(ni, j).x;
            double dx_eta = Nodes(ni, j + 1).x - Nodes(ni, j).x;

            // Geometric Metrics 
            double xi_x = dy_eta / vol;
            double xi_y = -dx_eta / vol;
            double eta_x = -dy_xi / vol;
            double eta_y = dx_xi / vol;

            double grad_xi_sq = xi_x * xi_x + xi_y * xi_y;
            double grad_eta_sq = eta_x * eta_x + eta_y * eta_y;

            // Viscous constraint factor
            double visc_factor = 2.0 * std::max(4.0 / 3.0, Config::GAMMA / Config::PRANDTL) * (mu / rho);

            // Constrained Spectral Radii
            double lambda_xi = std::abs(xi_x * u + xi_y * v) + a * std::sqrt(grad_xi_sq) + visc_factor * grad_xi_sq;
            double lambda_eta = std::abs(eta_x * u + eta_y * v) + a * std::sqrt(grad_eta_sq) + visc_factor * grad_eta_sq;

            double dt_local = maxCFL / (lambda_xi + lambda_eta);
            dt_min = std::min(dt_min, dt_local);
        }
    }
    return dt_min;
}

void NavierStokesSolver2D::computeFluxResidual(const Field2D<StateVec>& state_in, Field2D<StateVec>& residualOut) const {
    // 0. Zero out the residual array
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            residualOut(i, j) = StateVec();
        }
    }
    


    // 1. Get pure inviscid residual along Xi direction (i+1/2 faces)



    // 1. Get pure inviscid residual along Eta direction (j+1/2 faces)


    // 3. Viscous Fluxes along Xi direction (i+1/2 faces)
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {

        }
    }

    // 4. Viscous Fluxes along Eta direction (j+1/2 faces)
    for (int j = 0; j < ny - 1; ++j) {
        for (int i = 0; i < nx; ++i) {

        }
    }

    // 5. Viscous Wall boundary condition (at j=0)
    for (int i = 0; i < nx; ++i) {

    }

    // 6. Apply Far-field Boundary Fluxes (Direct flux subtraction)
    for (int i = 0; i < nx; ++i) {

    }

}
