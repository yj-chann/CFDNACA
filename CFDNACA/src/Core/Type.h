#pragma once
#include <cmath>
#include <vector>
#include"config.h"



// Physical constants
//constexpr double GAMMA = 1.4;
//constexpr double GAMMA_MINUS_ONE = GAMMA - 1.0;

// State vector: [rho, rho*u, rho*v, rho*E]
struct StateVec {
    double rho, rhou, rhov, rhoE;

    StateVec() : rho(0), rhou(0), rhov(0), rhoE(0) {}
    StateVec(double r, double ru, double rv, double rE) : rho(r), rhou(ru), rhov(rv), rhoE(rE) {}

    // Primitive variables extraction
    double u() const { return rhou / rho; }
    double v() const { return rhov / rho; }
    double p() const { return Config::GAMMA_MINUS_ONE * (rhoE - 0.5 * rho * (u() * u() + v() * v())); }
    double H() const { return (rhoE + p()) / rho; }
    double a() const { return std::sqrt(Config::GAMMA * p() / rho); }
    double T() const { return p() / (rho * Config::R); }

    // Operator overloads for vector math
    StateVec operator-(const StateVec& rhs) const {
        return StateVec(rho - rhs.rho, rhou - rhs.rhou, rhov - rhs.rhov, rhoE - rhs.rhoE);
    }
    StateVec operator+(const StateVec& rhs) const {
        return StateVec(rho + rhs.rho, rhou + rhs.rhou, rhov + rhs.rhov, rhoE + rhs.rhoE);
    }
    StateVec operator*(double scalar) const {
        return StateVec(rho * scalar, rhou * scalar, rhov * scalar, rhoE * scalar);
    }
};

struct FaceNormal {
    double nx, ny, length;
};


// Struct to hold spatial gradients (dx, dy) for Green-Gauss
struct Gradient2D {
    double dx, dy;
    Gradient2D() : dx(0.0), dy(0.0) {}
    Gradient2D(double dx, double dy) : dx(dx), dy(dy) {}
};

// 2D Array wrapper for computational grid variables
template <typename T>
class Field2D {
public:
    int nx, ny;
    std::vector<T> data;
    Field2D(int nx, int ny) : nx(nx), ny(ny), data(nx* ny) {}
    T& operator()(int i, int j) { return data[j * nx + i]; }
    const T& operator()(int i, int j) const { return data[j * nx + i]; }
};

struct Point2D {
    double x, y;
    Point2D() : x(0), y(0) {}
    Point2D(double x, double y) : x(x), y(y) {}
}; 
