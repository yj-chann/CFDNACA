#include "NavierStokesSolver2D.h"
#include "../Boundaries/Boundary.h"
#include "../Numerics/RoeFlux.h"
#include "../Numerics/ViscousFlux.h"
#include <algorithm>
#include <cmath>

NavierStokesSolver2D::NavierStokesSolver2D(int num_cells_x, int num_cells_y, int order, double cfl, const Field2D<Point2D>& meshNodes)
    : EulerSolver2D(num_cells_x, num_cells_y, order, cfl , meshNodes){
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

// i , 1/2
StateVec NavierStokesSolver2D::extrapolateToWall(const Field2D<StateVec>& state_in, int i) const {
    if (!spatialOrder) {
        double p = state_in(i, 0).p();
        // isothermal rho = p / (R * T_wall)
        double rho = p / (Config::R * Config::Tw);
        double rhoE = p / Config::GAMMA_MINUS_ONE;
        return StateVec(rho, 0.0, 0.0, rhoE);
    }
    else{
        double p = (state_in(i, 0).p() * 1.5 - state_in(i, 1).p() * 0.5);
        double rho = p / (Config::R * Config::Tw);
        double rhoE = p / Config::GAMMA_MINUS_ONE;
        return StateVec(rho, 0.0, 0.0, rhoE);
    }
}




Point2D NavierStokesSolver2D::getCellPos(int i, int j) const {
    // Averages the 4 surrounding nodes to find the cell center coordinate
    int ni = i;
    int nip1 = i + 1;
    int nj = j;
    int njp1 = j + 1;

    double cx = 0.25 * (Nodes(ni, nj).x + Nodes(nip1, nj).x + Nodes(nip1, njp1).x + Nodes(ni, njp1).x);
    double cy = 0.25 * (Nodes(ni, nj).y + Nodes(nip1, nj).y + Nodes(nip1, njp1).y + Nodes(ni, njp1).y);
    return Point2D(cx, cy);
}


// ni , nj are node indices (0 to nx for ni, 0 to ny for nj) - note that nodes are (nx+1) x (ny+1)
// i+1/2,j+1/2
StateVec NavierStokesSolver2D::interpolateToNode(const Field2D<StateVec>& U_in, int ni, int nj) const {
    // 1. Calculate horizontal weighting factor (alpha) - needed for both boundary and interior
    int ci_L = wrapXi(ni - 1);
    int ci_R = wrapXi(ni);

 
    int ni_L = (ni - 1) == -1 ? nx - 1 : ni-1;
    int ni_R = (ni + 1) == nx + 1 ? 1 : ni+1;

    Point2D p_L = Nodes(ni_L, nj);
    Point2D p_R = Nodes(ni_R, nj);
    Point2D p_node = Nodes(ni, nj);

    double d_L = std::hypot(p_L.x - p_node.x, p_L.y - p_node.y);
    double d_R = std::hypot(p_R.x - p_node.x, p_R.y - p_node.y);
    double alpha = d_R / (d_L + d_R);

    // 2. Boundary mapped node interpolation
    if (nj == 0) { // Solid Wall No-slip Isothermal Boundary mapping
        double p_row0 = alpha * U_in(ci_L, 0).p() + (1.0 - alpha) * U_in(ci_R, 0).p();
        double p_node = p_row0;

        if (spatialOrder) {
            double p_row1 = alpha * U_in(ci_L, 1).p() + (1.0 - alpha) * U_in(ci_R, 1).p();
            p_node = p_row0 * 1.5 - p_row1 * 0.5;
        }

        // isothermal rho_node = p_node / (R * T_wall)
        double rho_node = p_node / (Config::R * Config::Tw);

        // Reconstruct and return conservative StateVec   
        double rhoE_node = p_node / Config::GAMMA_MINUS_ONE;
        return StateVec(rho_node, 0.0, 0.0, rhoE_node);


    }

    if (nj >= ny) { // Far-field Characteristic Boundary Mapping
        // Horizontally interpolate U for the last two cell rows below the far-field
        StateVec U_row_ny1 = U_in(ci_L, ny - 1) * alpha + U_in(ci_R, ny - 1) * (1.0 - alpha);
        StateVec U_row_ny2 = U_in(ci_L, ny - 2) * alpha + U_in(ci_R, ny - 2) * (1.0 - alpha);

        StateVec U_star;
        if (!spatialOrder) {
            U_star = U_row_ny1;
        }
        else {
            U_star = U_row_ny1 * 1.5 - U_row_ny2 * 0.5;
        }

        FaceNormal n_L = FarfieldNormals[ci_L];
        FaceNormal n_R = FarfieldNormals[ci_R];
        double n_x = alpha * n_L.nx + (1.0 - alpha) * n_R.nx;
        double n_y = alpha * n_L.ny + (1.0 - alpha) * n_R.ny;
        double n_len = std::hypot(n_x, n_y);
        n_x /= (n_len);
        n_y /= (n_len);

        double Vn_inf = U_inf.u() * n_x + U_inf.v() * n_y;
        double a_inf = U_inf.a();
        double Vn_star = U_star.u() * n_x + U_star.v() * n_y;
        double a_star = U_star.a();

        double R1, R2, R3, R4;

        // lambda 1
        if (Vn_inf - a_inf < 0) { R1 = Vn_inf - 2.0 * a_inf / Config::GAMMA_MINUS_ONE; }
        else { R1 = Vn_star - 2.0 * a_star / Config::GAMMA_MINUS_ONE; }

        // lambda 2 (Entropy)
        if (Vn_inf < 0) { R2 = U_inf.p() / std::pow(U_inf.rho, Config::GAMMA); }
        else { R2 = U_star.p() / std::pow(U_star.rho, Config::GAMMA); }

        // lambda 3 (Tangential Velocity)
        if (Vn_inf < 0) { R3 = U_inf.v() * n_x - U_inf.u() * n_y; } // Vt_inf
        else { R3 = U_star.v() * n_x - U_star.u() * n_y; } // Vt_star

        // lambda 4
        if (Vn_inf + a_inf < 0) { R4 = Vn_inf + 2.0 * a_inf / Config::GAMMA_MINUS_ONE; }
        else { R4 = Vn_star + 2.0 * a_star / Config::GAMMA_MINUS_ONE; }

        // Solve for boundary primitives from R1, R2, R3, R4
        double Vn_b = 0.5 * (R1 + R4);
        double a_b = 0.25 * Config::GAMMA_MINUS_ONE * (R4 - R1);
        double Vt_b = R3;

        double rho_b = std::pow((a_b * a_b) / (Config::GAMMA * R2), 1.0 / Config::GAMMA_MINUS_ONE);
        double p_b = rho_b * (a_b * a_b) / Config::GAMMA;

        double u_b = Vn_b * n_x - Vt_b * n_y;
        double v_b = Vn_b * n_y + Vt_b * n_x;

        double rhoE_b = p_b / Config::GAMMA_MINUS_ONE + 0.5 * rho_b * (u_b * u_b + v_b * v_b);

        return StateVec(rho_b, rho_b * u_b, rho_b * v_b, rhoE_b);
    }

    // 3. Interior Node Interpolation (nj > 0 && nj < ny)
    int cj_B = nj - 1;
    int cj_T = nj;

    // Vertical weighting factor (beta)
    Point2D p_B = Nodes(ni, nj - 1);
    Point2D p_T = Nodes(ni, nj + 1);
    double d_B = std::hypot(p_B.x - p_node.x, p_B.y - p_node.y);
    double d_T = std::hypot(p_T.x - p_node.x, p_T.y - p_node.y);
    double beta = d_T / (d_B + d_T);

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
            int nip1 = i + 1;

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
    // 0. Zero out the residual array (cellcenter based)
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            residualOut(i, j) = StateVec();
        }
    }

	// 1. Get pure inviscid residual along Xi direction (i+1/2 faces)
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            StateVec UL = reconstructXi(state_in, i, j, false);
            StateVec UR = reconstructXi(state_in, i, j, true);
            StateVec flux = RoeFlux::computeFlux(UL, UR, NormalsXi(i, j));

            // Inviscid fluxes are subtracted from left cell, added to right cell
            residualOut(i, j) = residualOut(i, j) - flux;
            residualOut(wrapXi(i + 1), j) = residualOut(wrapXi(i + 1), j) + flux;
        }
    }

    // 2. Get pure inviscid residual along Eta direction (j+1/2 faces)
    for (int j = 0; j < ny - 1; ++j) {
        for (int i = 0; i < nx; ++i) {
            StateVec UL = reconstructEta(state_in, i, j, false);
            StateVec UR = reconstructEta(state_in, i, j, true);
            StateVec flux = RoeFlux::computeFlux(UL, UR, NormalsEta(i, j));

            residualOut(i, j) = residualOut(i, j) - flux;
            residualOut(i, j + 1) = residualOut(i, j + 1) + flux;
        }
    }

    // 3. Viscous Fluxes along Xi direction (i+1/2 faces)
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            int i_next = wrapXi(i + 1);
            int ni = i + 1;

            Point2D c_curr = getCellPos(i, j);
            Point2D c_next = getCellPos(i_next, j);
            const FaceNormal& xiNormal = NormalsXi(i, j);

            // Contour integral differences
            double d_ix = c_next.x - c_curr.x;
            double d_iy = c_next.y - c_curr.y;
            double d_jx = -xiNormal.ny * xiNormal.length;
            double d_jy = xiNormal.nx * xiNormal.length;
            double area = d_ix * d_jy - d_jx * d_iy;

            StateVec U_curr = state_in(i, j);
            StateVec U_next = state_in(i_next, j);
            StateVec U_bottom = interpolateToNode(state_in, ni, j);
            StateVec U_top = interpolateToNode(state_in, ni, j + 1);

            // Green-Gauss Gradient Lambda
            auto calcGrad = [&](double phi_curr, double phi_next, double phi_bottom, double phi_top) {
                double d_i_phi = phi_next - phi_curr;
                double d_j_phi = phi_top - phi_bottom;
                Gradient2D g;
                g.dx = (d_i_phi * d_jy - d_j_phi * d_iy) / area;
                g.dy = (d_j_phi * d_ix - d_i_phi * d_jx) / area;
                return g;
                };

            Gradient2D g_u = calcGrad(U_curr.u(), U_next.u(), U_bottom.u(), U_top.u());
            Gradient2D g_v = calcGrad(U_curr.v(), U_next.v(), U_bottom.v(), U_top.v());
            Gradient2D g_T = calcGrad(U_curr.T(), U_next.T(), U_bottom.T(), U_top.T());

            StateVec U_face = (U_curr + U_next) * 0.5; // Average state at face
            double mu_face = computeViscosity(U_face.T());
            StateVec visc_flux = ViscousFlux::computeFlux(U_face, g_u, g_v, g_T, xiNormal, mu_face);

            // Viscous fluxes have the opposite sign mapping compared to inviscid fluxes
            residualOut(i, j) = residualOut(i, j) + visc_flux;
            residualOut(i_next, j) = residualOut(i_next, j) - visc_flux;
        }
    }

    // 4. Viscous Fluxes along Eta direction (j+1/2 faces)
    for (int j = 0; j < ny - 1; ++j) {
        for (int i = 0; i < nx; ++i) {
            int i_right = (i + 1 < Nodes.nx) ? i + 1 : 0;

            Point2D c_curr = getCellPos(i, j);
            Point2D c_next = getCellPos(i, j + 1);
            const FaceNormal& etaNormal = NormalsEta(i, j);

            double d_ix = etaNormal.ny * etaNormal.length;
            double d_iy = -etaNormal.nx * etaNormal.length;
            double d_jx = c_next.x - c_curr.x;
            double d_jy = c_next.y - c_curr.y;
            double area = d_ix * d_jy - d_jx * d_iy;

            StateVec U_curr = state_in(i, j);
            StateVec U_next = state_in(i, j + 1);
            StateVec U_left = interpolateToNode(state_in, i, j + 1);
            StateVec U_right = interpolateToNode(state_in, i_right, j + 1);

            auto calcGrad = [&](double phi_curr, double phi_next, double phi_left, double phi_right) {
                double d_i_phi = phi_right - phi_left;
                double d_j_phi = phi_next - phi_curr;
                Gradient2D g;
                g.dx = (d_i_phi * d_jy - d_j_phi * d_iy) / area;
                g.dy = (d_j_phi * d_ix - d_i_phi * d_jx) / area;
                return g;
                };

            Gradient2D g_u = calcGrad(U_curr.u(), U_next.u(), U_left.u(), U_right.u());
            Gradient2D g_v = calcGrad(U_curr.v(), U_next.v(), U_left.v(), U_right.v());
            Gradient2D g_T = calcGrad(U_curr.T(), U_next.T(), U_left.T(), U_right.T());

            StateVec U_face = (U_curr + U_next) * 0.5;
            double mu_face = computeViscosity(U_face.T());
            StateVec visc_flux = ViscousFlux::computeFlux(U_face, g_u, g_v, g_T, etaNormal, mu_face);

            residualOut(i, j) = residualOut(i, j) + visc_flux;
            residualOut(i, j + 1) = residualOut(i, j + 1) - visc_flux;
        }
    }

    // 5. Viscous Wall boundary condition (at j=1/2)
    for (int i = 0; i < nx; ++i) {
        int i_right = i + 1;

        Point2D c_curr = getCellPos(i, 0);
        Point2D n_left = Nodes(i, 0);
        Point2D n_right = Nodes(i_right, 0);
        Point2D face_center((n_left.x + n_right.x) * 0.5, (n_left.y + n_right.y) * 0.5);

        double d_ix = n_right.x - n_left.x;
        double d_iy = n_right.y - n_left.y;
        double d_jx = c_curr.x - face_center.x; // Delta j mapped from wall center to cell center
        double d_jy = c_curr.y - face_center.y;
        double area = d_ix * d_jy - d_jx * d_iy;

        StateVec U_curr = state_in(i, 0);
        StateVec U_wall = extrapolateToWall(state_in, i);
        StateVec U_left = interpolateToNode(state_in, i, 0);
        StateVec U_right = interpolateToNode(state_in, i_right, 0);

   

        auto calcGrad = [&](double phi_curr, double phi_wall, double phi_left, double phi_right) {
            double d_i_phi = phi_right - phi_left;
            double d_j_phi = phi_curr - phi_wall;
            Gradient2D g;
            g.dx = (d_i_phi * d_jy - d_j_phi * d_iy) / area;
            g.dy = (d_j_phi * d_ix - d_i_phi * d_jx) / area;
            return g;
            };

        Gradient2D g_u = calcGrad(U_curr.u(), U_wall.u(), U_left.u(), U_right.u());
        Gradient2D g_v = calcGrad(U_curr.v(), U_wall.v(), U_left.v(), U_right.v());
        Gradient2D g_T = calcGrad(U_curr.T(), U_wall.T(), U_left.T(), U_right.T());

        double mu_wall = computeViscosity(U_wall.T());

        // ComputeCompleteNSWallFlux returns: (wall_inviscid_flux - wall_viscous_flux)
        StateVec wall_flux = BoundaryConditions::computeCompleteNSWallFlux(U_wall, g_u, g_v, g_T, WallNormals[i], mu_wall);
        residualOut(i, 0) = residualOut(i, 0) + wall_flux; // Add combined physical flux inward to cell
    }

    // 6. Apply Far-field Boundary Fluxes (Direct flux subtraction)
    for (int i = 0; i < nx; ++i) {
        StateVec U_star = extrapolateToFarfield(state_in, i);
        StateVec far_flux = BoundaryConditions::computeFarFieldFlux(U_star, U_inf, FarfieldNormals[i]);
        residualOut(i, ny - 1) = residualOut(i, ny - 1) - far_flux;
    }
}