#include "ViscousFlux.h"
#include "../Core/config.h"

StateVec ViscousFlux::computeFlux(
    const StateVec& U_face,
    const Gradient2D& g_u,
    const Gradient2D& g_v,
    const Gradient2D& g_T,
    const FaceNormal& normal,
    double mu
) {
    double nx = normal.nx;
    double ny = normal.ny;
    double ds = normal.length;

    // Velocity components at face
    double u = U_face.u();
    double v = U_face.v();

    // Thermal conductivity k = (mu * Cp) / Pr
    double k = (mu * Config::CP) / Config::PRANDTL;

    // Stokes' Hypothesis
    double lambda = -(2.0 / 3.0) * mu;
    double div_V = g_u.dx + g_v.dy; // u_x + v_y

    // Stress Tensor Components
    double tau_xx = 2.0 * mu * g_u.dx + lambda * div_V;
    double tau_yy = 2.0 * mu * g_v.dy + lambda * div_V;
    double tau_xy = mu * (g_u.dy + g_v.dx);

    // Heat Flux (Fourier's Law)
    double q_x = -k * g_T.dx;
    double q_y = -k * g_T.dy;

    // Assemble Viscous Flux Vector H_v
    StateVec Hv;
    Hv.rho = 0.0;
    Hv.rhou = (tau_xx * nx + tau_xy * ny) * ds;
    Hv.rhov = (tau_xy * nx + tau_yy * ny) * ds;

    // Energy equation viscous terms
    double work_heat_x = u * tau_xx + v * tau_xy - q_x;
    double work_heat_y = u * tau_xy + v * tau_yy - q_y;
    Hv.rhoE = (work_heat_x * nx + work_heat_y * ny) * ds;

    return Hv;
}