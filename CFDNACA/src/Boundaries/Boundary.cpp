#include "Boundary.h"
#include <cmath>

StateVec BoundaryConditions::computeSolidWallFlux(const StateVec& U_star, const FaceNormal& inward_n) {
    double nx = inward_n.nx;
    double ny = inward_n.ny;

    double rho_star = U_star.rho;
    double u_star = U_star.u();
    double v_star = U_star.v();
    double p_star = U_star.p();
    double a_star = U_star.a();

    double Vn_star = u_star * nx + v_star * ny;

    // Riemann Invariants & Isentropic relations
    double term = (GAMMA_MINUS_ONE * GAMMA_MINUS_ONE * std::pow(Vn_star - (2.0 * a_star) / GAMMA_MINUS_ONE, 2.0)
        * std::pow(rho_star, GAMMA)) / (4.0 * GAMMA * p_star);

    double rho_b = std::pow(term, 1.0 / GAMMA_MINUS_ONE);
    double p_b = p_star * std::pow(rho_b / rho_star, GAMMA);

    // Direct Flux Assembly for Solid Wall
    StateVec wall_flux;
    wall_flux.rho = 0.0;
    wall_flux.rhou = p_b * inward_n.nx * inward_n.length;
    wall_flux.rhov = p_b * inward_n.ny * inward_n.length;
    wall_flux.rhoE = 0.0;

    return wall_flux;
}

StateVec BoundaryConditions::computeFarFieldFlux(const StateVec& U_star, const StateVec& U_inf, const FaceNormal& outward_n) {
    double nx = outward_n.nx;
    double ny = outward_n.ny;

    double Vn_inf = U_inf.u() * nx + U_inf.v() * ny;
    double a_inf = U_inf.a();
    double Vn_star = U_star.u() * nx + U_star.v() * ny;
    double a_star = U_star.a();

    double R1, R2, R3, R4;

    // lambda 1
    if (Vn_inf - a_inf < 0) { R1 = Vn_inf - 2.0 * a_inf / GAMMA_MINUS_ONE; }
    else { R1 = Vn_star - 2.0 * a_star / GAMMA_MINUS_ONE; }

    // lambda 2 (Entropy)
    if (Vn_inf < 0) { R2 = U_inf.p() / std::pow(U_inf.rho, GAMMA); }
    else { R2 = U_star.p() / std::pow(U_star.rho, GAMMA); }

    // lambda 3 (Tangential Velocity)
    if (Vn_inf < 0) { R3 = U_inf.v() * nx - U_inf.u() * ny; } // Vt_inf
    else { R3 = U_star.v() * nx - U_star.u() * ny; } // Vt_star

    // lambda 4
    if (Vn_inf + a_inf < 0) { R4 = Vn_inf + 2.0 * a_inf / GAMMA_MINUS_ONE; }
    else { R4 = Vn_star + 2.0 * a_star / GAMMA_MINUS_ONE; }

    // Solve for boundary primitives from R1, R2, R3, R4
    double Vn_b = 0.5 * (R1 + R4);
    double a_b = 0.25 * GAMMA_MINUS_ONE * (R4 - R1);
    double Vt_b = R3;

    double rho_b = std::pow((a_b * a_b) / (GAMMA * R2), 1.0 / GAMMA_MINUS_ONE);
    double p_b = rho_b * (a_b * a_b) / GAMMA;

    double u_b = Vn_b * nx - Vt_b * ny;
    double v_b = Vn_b * ny + Vt_b * nx;

    double rhoE_b = p_b / GAMMA_MINUS_ONE + 0.5 * rho_b * (u_b * u_b + v_b * v_b);

    // Direct Flux Assembly for Far-field
    StateVec far_flux;
    far_flux.rho = rho_b * Vn_b * outward_n.length;
    far_flux.rhou = (rho_b * u_b * Vn_b + p_b * outward_n.nx) * outward_n.length;
    far_flux.rhov = (rho_b * v_b * Vn_b + p_b * outward_n.ny) * outward_n.length;
    far_flux.rhoE = (rhoE_b + p_b) * Vn_b * outward_n.length;

    return far_flux;
}