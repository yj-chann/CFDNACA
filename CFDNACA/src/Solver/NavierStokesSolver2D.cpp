#include "NavierStokesSolver2D.h"
#include "../Boundaries/Boundary.h"
#include <algorithm>
#include <cmath>

NavierStokesSolver2D::NavierStokesSolver2D(int num_cells_x, int num_cells_y, int order, double cfl, const Field2D<Point2D>& meshNodes)
    : EulerSolver2D(num_cells_x, num_cells_y, order, cfl),
    Nodes(meshNodes) {
}


double NavierStokesSolver2D::computeViscosity(double T) const {
    // Non-dimensional free-stream viscosity derived from the Reynolds number definition
    double a_inf = std::sqrt(Config::GAMMA * Config::p_inf / Config::rho_inf);
    double mu_inf = (Config::rho_inf * (Config::Mach_inf * a_inf) * Config::L_REF) / Config::REYNOLDS;
    // Note: Assuming constant viscosity scaled to free-stream. 
    // Sutherland's law can be implemented here using local T if a dimensional T_ref is provided.
    return mu_inf;
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

StateVec NavierStokesSolver2D::interpolateToNode(const Field2D<StateVec>& U_in, int ni, int nj) const {
    // Boundary mapped node interpolation
    if (nj == 0) { // Solid Wall No-slip Isothermal Boundary mapping
        int ci = wrapXi(ni);
        double p_wall = 1.5 * U_in(ci, 0).p() - 0.5 * U_in(ci, 1).p(); // Pressure Extrapolation
        double T_w = Config::p_inf / (Config::rho_inf * Config::R);    // Isothermal (Free-stream Temp)
        double rho_wall = p_wall / (Config::R * T_w);
        return StateVec(rho_wall, 0.0, 0.0, p_wall / Config::GAMMA_MINUS_ONE); // u=0, v=0
    }
    if (nj >= ny) return U_inf; // Simplistic Far-field mapping

    // Fetch the 4 surrounding cell centers
    int ci_L = wrapXi(ni - 1);
    int ci_R = wrapXi(ni);
    int cj_B = nj - 1;
    int cj_T = nj;

    // Weighting factors: Inverse distance weights (alpha, beta)
    Point2D p_L = getCellPos(ci_L, cj_B);
    Point2D p_R = getCellPos(ci_R, cj_B);
    Point2D p_node = Nodes(ni, nj);

    double d_L = std::hypot(p_L.x - p_node.x, p_L.y - p_node.y);
    double d_R = std::hypot(p_R.x - p_node.x, p_R.y - p_node.y);
    double alpha = d_R / (d_L + d_R + 1e-14); // Avoid div by zero

    Point2D p_B = getCellPos(ci_R, cj_B);
    Point2D p_T = getCellPos(ci_R, cj_T);
    double d_B = std::hypot(p_B.x - p_node.x, p_B.y - p_node.y);
    double d_T = std::hypot(p_T.x - p_node.x, p_T.y - p_node.y);
    double beta = d_T / (d_B + d_T + 1e-14);

    StateVec U_LB = U_in(ci_L, cj_B);
    StateVec U_RB = U_in(ci_R, cj_B);
    StateVec U_LT = U_in(ci_L, cj_T);
    StateVec U_RT = U_in(ci_R, cj_T);

    // Bilinear distance-weighted interpolation
    return U_LB * (alpha * beta) +
        U_RB * ((1.0 - alpha) * beta) +
        U_LT * (alpha * (1.0 - beta)) +
        U_RT * ((1.0 - alpha) * (1.0 - beta));
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
    // 1. Get pure inviscid residual
    EulerSolver2D::computeFluxResidual(state_in, residualOut);

    // 2. Viscous Fluxes along Xi direction (i+1/2 faces)
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            int ip1 = wrapXi(i + 1);
            StateVec U_L = state_in(i, j);
            StateVec U_R = state_in(ip1, j);

            // Interpolate to Top and Bottom nodes defining the face
            int n_i = (i + 1 < Nodes.nx) ? i + 1 : 0;
            StateVec N_top = interpolateToNode(state_in, n_i, j + 1);
            StateVec N_bot = interpolateToNode(state_in, n_i, j);

            // Difference Operators: Delta_i (Cell R - Cell L), Delta_j (Node Top - Node Bot)
            double du_i = U_R.u() - U_L.u();   double du_j = N_top.u() - N_bot.u();
            double dv_i = U_R.v() - U_L.v();   double dv_j = N_top.v() - N_bot.v();
            double dT_i = U_R.T() - U_L.T();   double dT_j = N_top.T() - N_bot.T();

            Point2D c_L = getCellPos(i, j);
            Point2D c_R = getCellPos(ip1, j);

            double dx_i = c_R.x - c_L.x;       double dy_i = c_R.y - c_L.y;
            double dx_j = Nodes(n_i, j + 1).x - Nodes(n_i, j).x;
            double dy_j = Nodes(n_i, j + 1).y - Nodes(n_i, j).y;

            // Green-Gauss Contour Area Omega
            double Omega = dx_i * dy_j - dx_j * dy_i;
            Omega = (std::abs(Omega) < 1e-14) ? 1e-14 : Omega;

            // Gradient Reconstruction Formulas
            Gradient2D g_u((du_i * dy_j - du_j * dy_i) / Omega, (du_j * dx_i - du_i * dx_j) / Omega);
            Gradient2D g_v((dv_i * dy_j - dv_j * dy_i) / Omega, (dv_j * dx_i - dv_i * dx_j) / Omega);
            Gradient2D g_T((dT_i * dy_j - dT_j * dy_i) / Omega, (dT_j * dx_i - dT_i * dx_j) / Omega);

            StateVec U_face = (U_L + U_R) * 0.5;
            double mu_face = computeViscosity(U_face.T());

            StateVec v_flux = ViscousFlux::computeFlux(U_face, g_u, g_v, g_T, NormalsXi(i, j), mu_face);

            // Subtract Viscous flux (from Euler formulation) -> add to residualOut L, subtract from R
            residualOut(i, j) = residualOut(i, j) + v_flux;
            residualOut(ip1, j) = residualOut(ip1, j) - v_flux;
        }
    }

    // 3. Viscous Fluxes along Eta direction (j+1/2 faces)
    for (int j = 0; j < ny - 1; ++j) {
        for (int i = 0; i < nx; ++i) {
            StateVec U_B = state_in(i, j);
            StateVec U_T = state_in(i, j + 1);

            int nip1 = (i + 1 < Nodes.nx) ? i + 1 : 0;
            StateVec N_L = interpolateToNode(state_in, i, j + 1);
            StateVec N_R = interpolateToNode(state_in, nip1, j + 1);

            // Delta_i (Node Right - Node Left), Delta_j (Cell Top - Cell Bot)
            double du_i = N_R.u() - N_L.u();   double du_j = U_T.u() - U_B.u();
            double dv_i = N_R.v() - N_L.v();   double dv_j = U_T.v() - U_B.v();
            double dT_i = N_R.T() - N_L.T();   double dT_j = U_T.T() - U_B.T();

            Point2D c_B = getCellPos(i, j);
            Point2D c_T = getCellPos(i, j + 1);

            double dx_i = Nodes(nip1, j + 1).x - Nodes(i, j + 1).x;
            double dy_i = Nodes(nip1, j + 1).y - Nodes(i, j + 1).y;
            double dx_j = c_T.x - c_B.x;
            double dy_j = c_T.y - c_B.y;

            double Omega = dx_i * dy_j - dx_j * dy_i;
            Omega = (std::abs(Omega) < 1e-14) ? 1e-14 : Omega;

            Gradient2D g_u((du_i * dy_j - du_j * dy_i) / Omega, (du_j * dx_i - du_i * dx_j) / Omega);
            Gradient2D g_v((dv_i * dy_j - dv_j * dy_i) / Omega, (dv_j * dx_i - dv_i * dx_j) / Omega);
            Gradient2D g_T((dT_i * dy_j - dT_j * dy_i) / Omega, (dT_j * dx_i - dT_i * dx_j) / Omega);

            StateVec U_face = (U_B + U_T) * 0.5;
            double mu_face = computeViscosity(U_face.T());

            StateVec v_flux = ViscousFlux::computeFlux(U_face, g_u, g_v, g_T, NormalsEta(i, j), mu_face);

            residualOut(i, j) = residualOut(i, j) + v_flux;
            residualOut(i, j + 1) = residualOut(i, j + 1) - v_flux;
        }
    }

    // 4. Viscous Wall boundary condition override (at j=0)
    for (int i = 0; i < nx; ++i) {
        // Strip out the purely inviscid Euler slip-wall flux added by base class
        StateVec dummy_U_star = state_in(i, 0);
        residualOut(i, 0) = residualOut(i, 0) - BoundaryConditions::computeSolidWallFlux(dummy_U_star, WallNormals[i]);

        // Isothermal, No-Slip condition rules
        double p_wall = 1.5 * state_in(i, 0).p() - 0.5 * state_in(i, 1).p(); // Pressure Extrapolation
        double T_wall = Config::p_inf / (Config::rho_inf * Config::R);
        double rho_wall = p_wall / (Config::R * T_wall);

        Point2D c_B = getCellPos(i, 0);
        Point2D n_B = Point2D((Nodes(i, 0).x + Nodes((i + 1 < Nodes.nx) ? i + 1 : 0, 0).x) * 0.5,
            (Nodes(i, 0).y + Nodes((i + 1 < Nodes.nx) ? i + 1 : 0, 0).y) * 0.5); // Face midpoint

        // Half-cell gradient approximation for the wall
        double dy_w = c_B.y - n_B.y + 1e-14;
        double du_dy = (state_in(i, 0).u() - 0.0) / dy_w;
        double dv_dy = (state_in(i, 0).v() - 0.0) / dy_w;
        double dT_dy = (state_in(i, 0).T() - T_wall) / dy_w;

        Gradient2D g_u_w(0.0, du_dy);
        Gradient2D g_v_w(0.0, dv_dy);
        Gradient2D g_T_w(0.0, dT_dy);

        double mu_wall = computeViscosity(T_wall);
        StateVec U_wall(rho_wall, 0.0, 0.0, p_wall / Config::GAMMA_MINUS_ONE);

        StateVec ns_wall_flux = BoundaryConditions::computeViscousWallFlux(
            U_wall, g_u_w, g_v_w, g_T_w, WallNormals[i], mu_wall);

        // Add exact pressure component directly onto the boundary flux
        ns_wall_flux.rhou += p_wall * WallNormals[i].nx * WallNormals[i].length;
        ns_wall_flux.rhov += p_wall * WallNormals[i].ny * WallNormals[i].length;

        // Apply completed viscous + pressure wall flux
        residualOut(i, 0) = residualOut(i, 0) + ns_wall_flux;
    }
}