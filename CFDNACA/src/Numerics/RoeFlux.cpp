#include "RoeFlux.h"
#include <algorithm>

StateVec RoeFlux::physicalFlux(const StateVec& U, double nx, double ny) {
    double u = U.u();
    double v = U.v();
    double p = U.p();
    double Vn = u * nx + v * ny;

    return StateVec(
        U.rho * Vn,
        U.rhou * Vn + p * nx,
        U.rhov * Vn + p * ny,
        (U.rhoE + p) * Vn
    );
}

StateVec RoeFlux::computeFlux(const StateVec& UL, const StateVec& UR, const FaceNormal& normal) {
    double nx = normal.nx;
    double ny = normal.ny;

    double uL = UL.u(), vL = UL.v(), pL = UL.p(), HL = UL.H();
    double uR = UR.u(), vR = UR.v(), pR = UR.p(), HR = UR.H();

    double sq_rhoL = std::sqrt(UL.rho);
    double sq_rhoR = std::sqrt(UR.rho);
    double Rv = sq_rhoR / sq_rhoL;

    double rho_tilde = sq_rhoL * sq_rhoR;
    double u_tilde = (uL + Rv * uR) / (1.0 + Rv);
    double v_tilde = (vL + Rv * vR) / (1.0 + Rv);
    double H_tilde = (HL + Rv * HR) / (1.0 + Rv);
    double a_tilde = std::sqrt(Config::GAMMA_MINUS_ONE * (H_tilde - 0.5 * (u_tilde * u_tilde + v_tilde * v_tilde)));

    double u_breve = u_tilde * nx + v_tilde * ny;
    double lambda1 = u_breve;
    double lambda2 = u_breve + a_tilde;
    double lambda3 = u_breve - a_tilde;

    double d_rho = UR.rho - UL.rho;
    double d_u = uR - uL;
    double d_v = vR - vL;
    double d_p = pR - pL;

    double abs_l1 = std::abs(lambda1);
    double abs_l2 = std::abs(lambda2);
    double abs_l3 = std::abs(lambda3);
    double a_tilde2 = a_tilde * a_tilde;

    double beta1 = abs_l1 * (d_rho - d_p / a_tilde2);
    double beta2 = (abs_l2 / (2.0 * a_tilde2)) * (d_p + rho_tilde * a_tilde * (nx * d_u + ny * d_v));
    double beta3 = (abs_l3 / (2.0 * a_tilde2)) * (d_p - rho_tilde * a_tilde * (nx * d_u + ny * d_v));

    double beta4 = beta1 + beta2 + beta3;
    double beta5 = a_tilde * (beta2 - beta3);
    double beta7 = abs_l1 * rho_tilde * (d_v * nx - d_u * ny);

    StateVec dissipation;
    dissipation.rho = beta4;
    dissipation.rhou = u_tilde * beta4 + nx * beta5 - ny * beta7;
    dissipation.rhov = v_tilde * beta4 + ny * beta5 + nx * beta7;
    dissipation.rhoE = H_tilde * beta4 + (u_tilde * nx + v_tilde * ny) * beta5
        + (v_tilde * nx - u_tilde * ny) * beta7 - (a_tilde2 * beta1) / Config::GAMMA_MINUS_ONE;

    StateVec FL = physicalFlux(UL, nx, ny);
    StateVec FR = physicalFlux(UR, nx, ny);

    StateVec flux;
    flux.rho = 0.5 * (FL.rho + FR.rho - dissipation.rho) * normal.length;
    flux.rhou = 0.5 * (FL.rhou + FR.rhou - dissipation.rhou) * normal.length;
    flux.rhov = 0.5 * (FL.rhov + FR.rhov - dissipation.rhov) * normal.length;
    flux.rhoE = 0.5 * (FL.rhoE + FR.rhoE - dissipation.rhoE) * normal.length;

    return flux;
}