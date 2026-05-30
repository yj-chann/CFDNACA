#pragma once

namespace Config {
    // ---------------------------------------------------------
    // Math Constants
    // ---------------------------------------------------------
    constexpr double M_PI = 3.14159265358979323846;

    // ---------------------------------------------------------
    // Fluid & Thermodynamic Constants
    // ---------------------------------------------------------
    constexpr double GAMMA = 1.4;
    constexpr double GAMMA_MINUS_ONE = GAMMA - 1.0;

    constexpr double R = 287.05;
    constexpr double CP = (GAMMA * R) / GAMMA_MINUS_ONE;
    constexpr double PRANDTL = 0.72;

    // ---------------------------------------------------------
    // Free-Stream Aerodynamic Conditions
    // ---------------------------------------------------------
    constexpr double Mach_inf = 0.8;
    constexpr double alpha_deg = 1.25;
    constexpr double alpha_rad = alpha_deg * M_PI / 180.0;

    constexpr double rho_inf = 1.0;
    constexpr double p_inf = 1.0;

    // ---------------------------------------------------------
    // Viscous & Reference Parameters
    // ---------------------------------------------------------
    // Re is required to compute non-dimensional viscosity!
    constexpr double REYNOLDS = 6.5e6; // Typical for a NACA validation case
    constexpr double L_REF = 1.0;      // Chord length

    // ---------------------------------------------------------
    // Solver Control Parameters
    // ---------------------------------------------------------
    constexpr double MAX_CFL = 0.8;
    constexpr int SPATIAL_ORDER = 2;
}
